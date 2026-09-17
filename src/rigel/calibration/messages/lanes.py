"""The LEVEL LANES — the gDNA lane and one RNA lane per strand — and the class they are built on.

       Gate: ``tests/calibration/test_transfer_rna_lanes.py``, ``test_transfer_policy.py``

A level is a population's density as an ABSOLUTE profile over ``u = log(rho / rho_ref)`` on the solve
grid, so it needs no map and no recipient and crosses the faces composition cannot (`LevelLane`).
`gdna_lane` serves every directed face left without a composition rule (`faces.Faces`); `rna_lanes`
serve each strand's faces read off the flag bits. Both are built by `TransferPolicy.prepare` from the
chain as the builders read it and the library's coordinates, which they take as arguments.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from . import Levels, Received
from .faces import NONE, RowTable, side_of
from .transfer_rows import (
    EPS,
    blur_row,
    count_logvar,
    flux_level,
    hop_price,
    intersect,
    level_of_profile,
    lower_side,
    poisson_level,
    profile_of_level,
    read_column,
    rna_level_of_profile,
    rna_row_of_level,
    strand_bits,
)

if TYPE_CHECKING:
    from .faces import Faces
    from .transfer import _Chain, _Library

__all__ = ["LevelLane", "gdna_lane", "rna_lanes"]


_FIELD = {"gdna": "level_gdna", "pos": "level_rna_pos", "neg": "level_rna_neg"}


def gdna_lane(c: _Chain, own: RowTable, faces: Faces, rho_ref: float) -> "LevelLane | None":
    """THE LEVEL LANE: the default of every directed face that has no composition rule. A gDNA level
    is absolute (a profile over ``u = log(rho / rho_ref)``), so it needs no map and no recipient: it
    crosses the faces composition cannot (strand changes, termini both ways, the AMBIG complex) and
    the EMPTY node — a piece with no total, which many exon pieces are — which forwards it
    unchanged. Every full node's own level: a gene edge's crossing as a Poisson level (structurally
    pure gDNA), any other node's own profile read through its total. ``rho_ref`` is the lane's
    coordinate, the library's structurally pure gDNA density (`TransferPolicy.library`); ``None``
    when the library has no positive density to serve as one."""
    empty = ~(c.n_u > 0.0) | ~(c.a_g > 0.0)
    if not rho_ref > 0.0:
        return None
    u = c.lam
    gene_edge = np.zeros(c.n, bool)
    bnd = np.flatnonzero(c.is_bnd)
    lo, hi_ = c.left[bnd], c.right[bnd]
    gene_edge[bnd] = ((lo >= 0) & c.is_intergenic[np.maximum(lo, 0)]) | (
        (hi_ >= 0) & c.is_intergenic[np.maximum(hi_, 0)]
    )
    own_level = RowTable(c.n, c.K)
    for x in np.flatnonzero(~empty & ~c.is_intergenic):
        if gene_edge[x]:
            own_level[x] = poisson_level(u, c.n_u[x], c.a_g[x], rho_ref)
        elif own[x] is not None and np.ptp(own[x]) > EPS:
            own_level[x] = level_of_profile(own[x], c.lam, u, c.n_u[x], c.a_g[x], rho_ref)
    # the lane's faces, per (destination, side): every face with no composition rule whose two ends
    # are not terminals — a terminal neither receives nor sends on the lane: nothing is imputed at
    # structurally pure gDNA, and it has no own level to send; its far side is another locus
    face = np.zeros((c.n, 2), bool)
    for side, nbr in enumerate((c.left, c.right)):
        s = np.maximum(nbr, 0)
        face[:, side] = (
            (nbr >= 0) & ~c.is_intergenic & ~c.is_intergenic[s] & (faces.kind[:, side] == NONE)
        )
    return LevelLane("gdna", u, c.lam, rho_ref, c.n_u, c.a_g, empty, own_level, face)


def rna_lanes(c: _Chain, own: RowTable, library: _Library) -> dict:
    """THE RNA LEVEL LANES, one per strand (the both-stranded locus).
    FACES from the flag bits: strand ``s``'s level crosses a face iff the boundary carries none of
    ``s``'s four bits and both nodes admit ``s``; across ``s``'s OWN junction it enters ``s``'s
    intron — the crossing IS the intron's unspliced population — and not ``s``'s exon; a terminus of
    ``s`` stops ``s`` both ways. TWO-SIDED only between an intron of ``s`` and its own boundary (one
    shared unspliced population; the intron test is PER STRAND). SOURCES: a single-strand node's
    own claim read as its live strand's RNA level, and the certified flux at each of an exon's
    junctions as that strand's level at the exon — the junction's estimate of the exon's RNA
    abundance, priced by THE NODE PAIR (the junction's spliced count at its route rate, whole-strand
    units, against the exon's count on the column that strand reads on, over the PROTOCOL'S SHARE of
    the exon's RNA opportunity, ``kappa_read · a_r`` — a column is that much opportunity for the strand's
    RNA to be counted on it, so the two densities are in one unit and the pair's agreement is priced as
    counting alone at every kappa; `EQUATIONS.md` §12), kept LOWER-SIDED because a two-sided
    estimate over-claims at the probe cliff; kept per FACE so the solve can tell which face's
    composition already carries it. An EMPTY exon piece beside a lit junction is a source too: its
    level is priced on its zero count — counting alone — and the piece emits it with the flux's own
    witness, the pooled spliced count on the pooled route opportunity, so the next full node prices
    the hop as a full exon prices its flux. The coordinate is the library's one RNA coordinate (`_Library.rho_rna`) — a level
    is absolute and its coordinate only an origin, so both strands share it and a strand with no
    single-strand exon of its own still builds its flux levels — and whether the split is a witness at
    all is the library's verdict too. Nothing pooled, no constant."""
    n_u, a_r, cnt = c.n_u, c.a_r, c.cnt
    empty = ~(n_u > 0.0) | ~(a_r > 0.0)
    single = ~(c.fp & c.fn)
    kappa = None if c.strand is None else float(c.strand[0])
    # the protocol's read rate: the share of a strand's RNA that reads on the column it reads on
    kappa_read = 0.5 if kappa is None else max(kappa, 1.0 - kappa)
    split_live = library.split_live
    lanes = {}
    for name, free, col in (("pos", c.fp, 0), ("neg", c.fn, 1)):
        all_bits, _sj_bits, term_bits = strand_bits[name]
        # the lane's WITNESS is the count of the reads strand-``s`` RNA produces: its own genome-strand
        # column when the library reads sense, the other under an antisense protocol (`read_column`).
        # ⛔ Under a strongly antisense protocol the own column holds almost nothing, so reading it
        # prices every hop as counting on an empty witness and blurs every floor away.
        col_read = read_column(col, kappa)
        intron_s = ~c.is_bnd & free & ~c.exon_of[name]
        # a face crosses when its boundary carries none of the strand's bits; across the strand's own
        # (non-terminus) bits it still enters the strand's intron, and every face into that intron is
        # two-sided — per (destination, side), the source being the neighbour on that side
        face = np.zeros((c.n, 2), bool)
        two_sided = np.zeros((c.n, 2), bool)
        dest = np.arange(c.n)
        for side, nbr in enumerate((c.left, c.right)):
            s = np.maximum(nbr, 0)
            ok = (nbr >= 0) & ~c.is_intergenic & free[s] & free
            b = np.where(c.is_bnd[s], s, dest)
            reg = np.where(c.is_bnd[s], dest, s)
            f = c.flags[b].astype(np.int64)
            crossing = (f & all_bits) == 0
            into_own_intron = ~crossing & ((f & term_bits) == 0) & intron_s[reg]
            face[:, side] = ok & (crossing | into_own_intron)
            two_sided[:, side] = ok & ((crossing & intron_s[reg]) | into_own_intron)
        rho_ref = float(library.rho_rna)
        own_level = RowTable(c.n, c.K)
        flux = RowTable((c.n, 2), c.K)  # the junction flux levels, per (exon, side of its junction)
        flux_witness = RowTable(c.n, 2)  # per empty flux source: (spliced count, route opportunity)
        if rho_ref > 0.0:
            for x in np.flatnonzero(free):
                parts = []
                if not empty[x] and single[x] and own[x] is not None and np.ptp(own[x]) > EPS:
                    parts.append(
                        rna_level_of_profile(own[x], c.lam, c.lam, n_u[x], a_r[x], rho_ref)
                    )
                c_sum = a_sum = 0.0
                if c.is_exon[x]:
                    # the junction's flux is a measurement of THIS exon's RNA whether or not the piece
                    # holds a fragment of its own: at an EMPTY piece — one shorter than a
                    # fragment, or dark — the level is built too, priced on the piece's zero count (the
                    # counting rule every hop pays), and the piece emits it with the flux's witness
                    for b in (c.left[x], c.right[x]):
                        if b < 0 or not c.is_bnd[b]:
                            continue
                        hi = 1 if c.left[x] == b else 0
                        c_j, r_j = float(c.sj_count[hi][b, col]), float(c.route_rate[hi][b, col])
                        if not (c_j > 0.0 and r_j > 0.0):
                            continue
                        # the junction's rate is whole-strand; the column count is priced on the
                        # protocol's share of the exon's opportunity, so the units agree (the exon's
                        # total density at kappa = ½, the column's at kappa → 1, the strand's own share
                        # at a both-stranded exon) — `EQUATIONS.md` §12
                        v = hop_price(c_j, c_j / r_j, cnt[x, col_read], kappa_read * a_r[x])
                        fl = flux_level(c.lam, c_j, r_j, rho_ref, v)
                        parts.append(fl)
                        flux[x, side_of(int(b), int(x))] = fl
                        c_sum += c_j
                        a_sum += c_j / r_j
                if parts:
                    own_level[x] = intersect(parts)
                    if empty[x]:
                        flux_witness[x] = (c_sum, a_sum)
        # the two columns as contiguous arrays: the native pass reads them as they are
        lanes[name] = LevelLane(
            name,
            c.lam,
            c.lam,
            rho_ref,
            np.ascontiguousarray(cnt[:, col_read]),
            a_r,
            empty,
            own_level,
            face,
            two_sided=two_sided,
            total=n_u,
            flux=flux,
            other=np.ascontiguousarray(cnt[:, 1 - col_read]) if split_live else None,
            flux_witness=flux_witness,
        )
    return lanes


class LevelLane:
    """ONE population's LEVEL LANE — the gDNA lane, or one strand's RNA lane — for one sweep. A level
    travels as an ABSOLUTE profile over ``u = log(rho / rho_ref)`` on the solve grid, so it needs no
    map and no recipient: an EMPTY node forwards it unchanged, a full node emits the INTERSECTION of
    its own level and what it holds — ⛔ bounds INTERSECT, they do not multiply: a product of
    one-sided claims ratchets along a chain into a hard bound at the noisiest node's mode — and a
    full recipient prices the hop and takes the level as a LOWER bound.

    The two kinds are one class because they differ by PARAMETERS and nothing else: the WITNESS the
    hop's price reads — every node's ``count`` on its opportunity ``a`` (the totals on the gDNA lane;
    on an RNA lane the strand's read column, and where the library's strand channel is live the
    column split's asymmetry against ``other``, the other column) — and ``two_sided``, the faces the
    WHOLE profile crosses (none on the gDNA lane: every gDNA hop is a lower bound; an intron and its
    own boundary on an RNA lane: one shared unspliced population). ``total`` is every node's total,
    through which a held level is read back as the node's composition (`row`); ``own_level`` each
    node's own level (a `RowTable`); ``faces`` the directed faces the lane serves; ``flux``, per
    ``(exon, side)`` — the junction's priced estimate of the exon's RNA, kept per FACE (a `RowTable`
    over the node's two sides) so the solve can tell which face's composition already carries it;
    ``flux_witness``, per EMPTY node
    whose own level is a flux level, the witness that level travels with — the pooled spliced count
    on the pooled route opportunity of the node's lit junctions (an empty node has no count of its
    own to stamp a level with), a `RowTable` of ``(count, opportunity)`` pairs. Every array is in the
    layout the native pass reads (`tables`)."""

    __slots__ = (
        "population",
        "field",
        "u",
        "lam",
        "rho_ref",
        "count",
        "total",
        "a",
        "empty",
        "own_level",
        "face",
        "two_sided",
        "flux",
        "other",
        "flux_witness",
    )

    def __init__(
        self,
        population: str,
        u,
        lam,
        rho_ref,
        count,
        a,
        empty,
        own_level,
        face,
        *,
        two_sided=None,
        total=None,
        flux=None,
        other=None,
        flux_witness=None,
    ):
        n = len(own_level)
        self.population, self.field = population, _FIELD[population]
        self.u, self.lam, self.rho_ref = u, lam, float(rho_ref)
        self.count, self.a, self.empty = count, a, empty
        self.total = count if total is None else total
        self.own_level = own_level
        # the faces the lane serves and the two-sided ones, per (destination, side)
        self.face = np.asarray(face, bool).reshape(n, 2)
        self.two_sided = (
            np.zeros((n, 2), bool) if two_sided is None else np.asarray(two_sided, bool)
        )
        # the junction flux levels, per (exon, side of its junction)
        self.flux = RowTable((n, 2), self.own_level.rows.shape[1]) if flux is None else flux
        self.other = other
        self.flux_witness = RowTable(n, 2) if flux_witness is None else flux_witness

    def tables(self) -> tuple:
        """The lane as the native pass reads it (`native.transfer_pass`'s ``lanes`` entry): which
        `Received` table it writes, its faces, its two-sided faces, its emptiness, the own levels with
        their mask, the witness counts and opportunities, the other column or ``None``, and the flux
        witnesses with their mask — the arrays themselves, nothing copied."""
        return (
            Received.LANES.index(self.field),
            self.face,
            self.two_sided,
            self.empty,
            self.own_level.rows,
            self.own_level.mask,
            self.count,
            self.a,
            self.other,
            self.flux_witness.rows,
            self.flux_witness.mask,
        )

    def serves(self, s: int, x: int) -> bool:
        """Does the lane carry its level across the face into ``x`` from ``s``?"""
        return bool(self.face[int(x), side_of(int(s), int(x))])

    def flux_at(self, x: int, side: int):
        """The junction flux level ``x`` holds at its ``side`` (0: from its left junction, 1: its right),
        or ``None``."""
        return self.flux[int(x), int(side)]

    def witness(self, y: int):
        """The strand's RNA count at ``y`` and its Poisson variance, read from the column split: the
        asymmetry ``count − other`` (gDNA splits evenly, so it cancels; the other strand's RNA reads on
        the other column) over the protocol's strand contrast ``|1 − 2κ|`` — a factor common to every
        node, so it cancels from every ratio the price takes and is left out. A non-positive asymmetry
        is a DARK node: no RNA of this strand is measurable there."""
        c_r, c_o = float(self.count[y]), float(self.other[y])
        return c_r - c_o, c_r + c_o

    def emit(self, s: int, x: int, levels: Levels) -> bool:
        """What ``s`` sends toward ``x`` on the lane, written into row ``x`` of ``levels`` — the pass's
        table, whose row ``s`` is what ``s`` holds from its far side. An empty node forwards what it
        holds unchanged — unless it is itself a FLUX SOURCE, whose level (intersected with what it
        holds) travels with the flux's witness; a full node the INTERSECTION of its own level and what
        it holds — the own level WHOLE across a two-sided face, its lower side everywhere else —
        stamped with its own witness. Returns whether anything was sent."""
        own = self.own_level[s]
        far = bool(levels.present[s])
        if self.empty[s] and own is None:
            if far:
                levels.forward(s, x)
            return far
        if own is not None and not self.two_sided[int(x), side_of(int(s), int(x))]:
            own = lower_side(own)
        parts = [p for p in (own, levels.profile[s] if far else None) if p is not None]
        if not parts:
            return False
        if self.empty[s]:
            n, a = self.flux_witness[s]
            levels.write(x, intersect(parts), n, a)
            return True
        rna_count = rna_var = None
        if self.other is not None:
            rna_count, rna_var = self.witness(s)
        levels.write(x, intersect(parts), self.count[s], self.a[s], rna_count, rna_var)
        return True

    def receive(self, levels: Levels, s: int, x: int) -> None:
        """What a FULL recipient ``x`` holds after the hop from ``s`` — row ``x`` re-priced in place from
        what `emit` wrote there: across a TWO-SIDED face the whole profile, everywhere else its lower
        side — and on EVERY face the hop's price: both witness counts' counting plus the disagreement,
        beyond its own counting, between the two nodes' estimates of the population's abundance
        (`hop_price` on the counts; the owner's rule per hop, nothing pooled). On an RNA lane the
        witness is the column split's asymmetry (`witness`), not the column count: the column holds
        gDNA's half, which jumps with every probe edge whether or not this strand's RNA is there, while
        the asymmetry is this strand's RNA alone. So a dark recipient (no measurable RNA of the strand)
        agrees with a dark claim and the claim arrives whole — the perfectly dark host intron's "no RNA
        of mine here" that resolves the tilt at an antisense exon's boundaries under capture — while a
        lit recipient disagrees with a claim from a dimmer node by the cliff between them and blurs it
        away. ⛔ There is no counting-only exemption: a faint intron whose neighbour's probe enriches it
        many-fold would otherwise carry its claim whole across the cliff. Where the library's strand
        channel is dead the column count is the witness."""
        x = int(x)
        n_sent, a_sent = float(levels.count[x]), float(levels.opportunity[x])
        if self.other is None or not levels.has_witness[x]:
            v = hop_price(n_sent, a_sent, self.count[x], self.a[x])
        else:
            v = float(count_logvar(n_sent) + count_logvar(self.count[x]))
            n_s, v_s = float(levels.rna_count[x]), float(levels.rna_count_var[x])
            n_x, v_x = self.witness(x)
            if n_s > 0.0 and n_x > 0.0:
                r = (n_x / float(self.a[x])) / (n_s / a_sent)
                v += max(0.0, float(np.log(r)) ** 2 - (v_s / (n_s * n_s) + v_x / (n_x * n_x)))
        p = (
            levels.profile[x]
            if self.two_sided[x, side_of(int(s), x)]
            else lower_side(levels.profile[x])
        )
        p = blur_row(p, self.u, v) if v > 0.0 else p
        rna_count = rna_var = None
        if self.other is not None:
            rna_count, rna_var = self.witness(x)
        levels.write(x, p, self.count[x], self.a[x], rna_count, rna_var)

    def row(self, profile, x: int):
        """A held level's profile read as ``x``'s composition row through its own total — a pure
        coordinate change (the level was priced on arrival): "at least this much gDNA" is a floor on
        the gDNA share; "at least this much RNA" of the node's live strand is a ceiling on it."""
        read = profile_of_level if self.population == "gdna" else rna_row_of_level
        return read(profile, self.u, self.lam, self.total[x], self.a[x], self.rho_ref)
