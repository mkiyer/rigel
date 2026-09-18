"""The LEVEL LANES — the gDNA lane and one RNA lane per strand — the class they are built on.

       Gate: ``tests/calibration/test_transfer_rna_lanes.py``, ``test_transfer_policy.py``

A level is a population's density as an ABSOLUTE profile over ``u = log(rho / rho_ref)`` on the solve
grid, so it needs no map and no recipient and crosses the faces composition cannot (`LevelLane`). The
gDNA lane serves every directed face left without a composition rule (`faces.Faces`); each RNA lane
serves its strand's faces read off the flag bits. Both are built by `native.transfer_prepare`
(`native/transfer_kernel.cpp`'s `gdna_lane` and `rna_lane`) into the tables `TransferPolicy.prepare`
allocates on a `LevelLane` — its faces, own levels, junction flux levels and flux witnesses; its
coordinates, witnesses and emptiness are the chain's arrays.
"""

from __future__ import annotations

import numpy as np

from . import Levels, Received
from .faces import RowTable, side_of
from .transfer_rows import blur_row, count_logvar, hop_price, intersect, lower_side

__all__ = ["LevelLane"]


_FIELD = {"gdna": "level_gdna", "pos": "level_rna_pos", "neg": "level_rna_neg"}


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
    through which the solve reads a held level back as the node's composition; ``own_level`` each
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
