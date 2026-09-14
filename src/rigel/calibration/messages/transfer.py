"""TransferPolicy — the COMPOSITION TRANSFER policy on the two-phase backbone: every node states its
own claim, every directed face carries a recipient's rule, and the two passes do the rest.

       Gate: ``tests/calibration/test_transfer_policy.py``

Every message is a COMPOSITION PROFILE — a max-normalised log-likelihood over the solve grid of the
destination's gDNA share — carried across ONE face by a derived map from `transfer_rows`, or a
population's LEVEL — an absolute profile over the log density — carried where composition cannot
cross, so no level ever crosses a capture cliff and no constant is anywhere. The policy has three
parts, and ``prepare`` is a table of contents: one named BUILDER per shipped message.

* Every node's OWN CLAIM (`_claims`): an intron's factory profile (its density against the
  intergenic background — the context's ``factory_rows``, the very array ψ adds as its λ-factor); an
  exon's or a boundary's own strand profile
  where the solver's derived strand deadband declares the channel live (``strand``); an
  intergenic|exon edge's gDNA COUNT (the level lane: the edge's crossing is structurally pure gDNA).
  ⛔ A claim is data only — never a belief, which already holds the prior and the neighbours.
* The RECIPIENT's rule per directed face: absent = STOP (composition cannot cross: the recipient
  holds silence); the identity = FORWARD (the two objects share one unspliced population exactly); a
  map = MODIFY (the face's arithmetic with its counting width and, where two witnesses exist, the
  pair's own discrepancy). The rules are TYPED TABLES (`faces.Faces`): every directed face is one of a
  node's two sides — it hears from its left neighbour or its right — so a rule is a KIND and its
  parameters at ``(destination, side)``, five kinds in all, and the passes read the table. The rules
  ARE the shipped messages, each built by the function named, which writes its rows into the table:

  - `_splice_faces` — the intron|exon face. intron ⇄ boundary: FORWARD both ways (one shared
    unspliced population; into the intron only under a shared single strand). boundary → exon at a
    LICENSED face (no terminus, the same strand set): the splice-in face map, the certified flux
    capping the claimable gDNA share, widened by the face's counting. exon → boundary at that face:
    the splice-in map read backwards, marginalised over the face's spliced-to-unspliced ratio.
  - `_edge_level` — intergenic|exon edge → exon: THE EDGE'S LEVEL, one-sided: the exon has at least
    the edge's gDNA density, at the count's Poisson width; nothing above (no local witness prices
    capture's enrichment of the interior), a zero count vacuous (darkness is not absence).
  - `_terminus_rules` — the exon|exon TERMINUS boundary ⇄ its OUTSIDE exon: the licence counts the
    spliced crossing, ``f_b = f_O (U_b + S_b) / U_b`` — the outside exon's message travels the
    splice-out map, the boundary's the face map with the spliced density. And a TERMINUS boundary →
    the region INSIDE it (exon|exon and exon|intron alike): THE LEVEL RULE. Composition cannot cross
    a terminus (new transcription starts or ends there), the gDNA LEVEL can: the boundary's OWN
    strand profile — a measurement, never what it holds — carried through the level-kept map (the
    inside's gDNA share is the boundary's share times its crossing density, over the inside's own
    total), its shape preserved, blurred by the two totals' counting and by the pair's own
    discrepancies: the excess of the totals' disagreement over counting, and, where the inside exon's
    strand channel is live, the excess of the two strand modes' disagreement over counting. Value
    kept, nothing pooled, no hypothesis chosen for the discrepancy. With no own claim the boundary
    still sends what every library measures: its crossing's total bounds the inside's gDNA density
    above. At an sj+terminus boundary a junction sharing the terminus does not change which flank is
    inside; its flux belongs to the flank on the junction's exon side.
  - `_alternative_splice_site` — the ALTERNATIVE SPLICE SITE ⇄ both flanks: the intron-side flank
    shares the full unspliced crossing, the exon-of-both flank the crossing plus the leaving isoform
    measured as the route flux; each pair widened by its own disagreement beyond counting, nothing
    pooled.
  - `lanes.gdna_lane` — EVERY OTHER DIRECTED FACE (a strand change, termini pointing both ways, the AMBIG
    complex, and every face into or out of an EMPTY node — a piece with no total) carries THE LEVEL
    LANE: the gDNA level as an ABSOLUTE profile over the log density, which needs no map and no
    recipient. An empty node forwards it unchanged; a full node emits the intersection of its own
    level's lower side and the level it holds; the recipient prices the hop (both totals' counting
    plus the abundance discrepancy beyond it) and takes it as a LOWER bound — a level that crosses a
    face says "at least this much gDNA" and nothing more. The only faces with no rule at all lead
    into intergenic regions (structurally pure gDNA) or off the chain.
  - `lanes.rna_lanes` — THE RNA LEVEL LANES: one lane per strand, its FACES from the flag bits (strand
    ``s`` crosses a face iff the boundary carries none of ``s``'s bits and both nodes admit ``s``;
    across ``s``'s own junction it enters ``s``'s INTRON — the intron test is per strand, from
    ``BlockContext.exon_pos`` / ``exon_neg`` — and not ``s``'s exon; a terminus of ``s`` stops ``s``),
    TWO-SIDED only between an intron of ``s`` and its own boundary (one shared unspliced population:
    the whole profile), lower-only everywhere else, and EVERY hop priced by the pair — both counts'
    counting plus the disagreement between the two nodes' estimates of the strand's abundance
    (`lanes.LevelLane.witness`). SOURCES: a single-strand node's own claim read as its live strand's RNA
    level; the certified flux at each of an exon's junctions as that strand's level at the exon (one
    hop, boundary → exon; two junctions pay their pair's disagreement beyond counting) — at an EMPTY
    exon piece too, which emits it with the flux's own witness. DELIVERED at AMBIG nodes as one row
    over ψ's (λ, θ) cube (`PsiMessage.cube_rows`): a lower bound on RNA+ is an upper bound on the
    gDNA share through the node's own strand counts — the side the gDNA lane cannot give (the bracket
    theorem, gated); and at SINGLE-STRAND nodes as a CEILING on the gDNA share, read only from a face
    that sent no composition (`_PreparedTransfer._ceilings`). The tilt needs no lane of its own.

* The two passes and the solve (`_PreparedTransfer.propagate` / `solve`): a node SENDS two
  things apart — its own claim (a measurement) and what it holds from its far side (an imputation) —
  and the recipient's rule decides what to do with each: a composition rule composes them (a witness
  product: profiles add) and maps the product; a level rule reads the measurement only, because an
  imputation is never re-issued as a level; a lane intersects them as bounds; either may stop — so a
  claim travels as far as the faces admit it, each hop charging its own counting width, and no node
  ever hears its own claim back (the forward pass composes only what came from the left, the backward
  pass only what came from the right). At the solve the two held profiles add, the two held levels
  intersect and join them as a constraint, and ψ fuses the row with the slot's own evidence and the
  prior.

The laws the policy keeps: the sender publishes its claim unchanged; the recipient decides; a no-claim
stays a no-claim — a flat profile, an absent factory, a context with no rows all deliver silence,
never a zero-filled channel; a message is built from the source's claim and the recipient's constants
and observations, never the recipient's belief.

THE LIBRARY (`TransferPolicy.library`, once per sweep over the whole chain, `_Library`): the three
level lanes' coordinates — the library's structurally pure gDNA density, and each strand's unspliced
density over its single-strand exons — and whether the strand split is a live RNA witness anywhere.
They are ratios of sums over structurally selected slots and one boolean, read from observations and
geometry only, and they are the ONLY things a message knows about slots outside its own block.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..simplex_logodds import CubeRow, strand_row_logodds
from . import BlockContext, ChainView, PsiMessage, Received
from .faces import EDGE, FORWARD, LEVEL, SPLICE_OUT, TRANSPORT, Faces, fuse, norm
from .lanes import gdna_lane, rna_lanes
from .transfer_rows import (
    EPS,
    SJ_FLAGS,
    boundary_shares_strand,
    count_logvar,
    edge_level_row,
    face_is_licensed,
    face_map_lambda,
    intersect,
    junction_exon_side,
    junction_flanks,
    level_bound_row,
    level_map_lambda,
    lower_side,
    outside_flank,
    read_column,
)

__all__ = ["TransferPolicy"]

#: the three populations' lanes: the `Message` field each population's level travels on

_RNA = ("pos", "neg")


@dataclass(frozen=True, slots=True)
class _Library:
    """What the whole library says, reduced once per sweep and read by every block — the only
    cross-block information a message may use. ``rho_gdna``: the structurally pure gDNA density, the
    gDNA lane's coordinate (``0.0``: no positive density anywhere, so no level lane can be built).
    ``rho_rna``: the RNA lanes' one coordinate — the unspliced RNA density over the library's
    single-strand exons, both strands pooled, or over every exon when no single-strand exon has counts
    (a level is absolute and the coordinate only its origin, so one serves both strands and a strand with
    no single-strand exon of its own still builds its flux levels). ``split_live``: the strand split is a witness of a strand's RNA somewhere — the
    derived deadband is open and some single-strand exon has counts."""

    rho_gdna: float
    rho_rna: float
    split_live: bool


class TransferPolicy:
    """``strand = (kappa, od_g, od_r)`` is the library's fitted strand model, which the exon's and the
    boundary's own claims need; ``None`` leaves every strand claim off. The intron factory's per-slot
    profiles arrive on the context (``factory_rows``, the sweep's own λ-factor on its own grid);
    ``None`` there means the factory has nothing to claim; the faces and the lanes are built as with
    one."""

    name = "transfer"

    def __init__(self, strand: tuple[float, float, float] | None = None):
        self._strand = None if strand is None else tuple(float(x) for x in strand)

    # ── the library: the three lanes' coordinates and the strand witness's liveness, once ───────────
    def library(self, view: ChainView) -> _Library:
        is_bnd = np.asarray(view.is_boundary, bool)
        is_exon = np.asarray(view.is_exon_region, bool)
        fp, fn = np.asarray(view.free_pos, bool), np.asarray(view.free_neg, bool)
        is_intergenic = ~is_bnd & ~is_exon & ~fp & ~fn
        n_u = np.asarray(view.n_slot, np.float64)
        a_g = np.asarray(view.eff_gdna, np.float64)
        a_r = np.asarray(view.eff_rna, np.float64)
        cnt = np.asarray(view.unspliced_count, np.float64)
        # the gDNA lane's coordinate: the library's structurally pure gDNA density. It is a REFERENCE
        # for a log axis, so any positive density serves; when the intergenic count is zero (a
        # zero-gDNA library with no intergenic reads) the library-wide density stands in. ⛔ Without
        # that fallback the coordinate is zero, and the gDNA lane plus every RNA lane hanging off it
        # is never built at all — a gDNA-free library then gets no level messages anywhere.
        pure = is_intergenic & (a_g > 0.0)
        rho = float(n_u[pure].sum() / a_g[pure].sum()) if a_g[pure].sum() > 0.0 else 0.0
        if not rho > 0.0:
            rho = float(n_u[a_g > 0.0].sum() / max(a_g[a_g > 0.0].sum(), EPS))
        if not rho > 0.0:
            rho = 0.0
        kappa = None if self._strand is None else float(self._strand[0])
        # the split is a witness of a strand's RNA only where the library's strand channel is live —
        # the derived deadband's verdict (`ChainView.strand_live`): with it open, every counted
        # single-strand exon's strand precision is positive; with it shut, none is (an intron's
        # factory precision joins its strand term, so introns cannot stand for the channel)
        single = fp != fn
        split_live = (
            kappa is not None
            and bool(view.strand_live)
            and bool(np.any(is_exon & single & (n_u > 0.0)))
        )
        # the RNA lanes' coordinate: the unspliced RNA density over the single-strand exons, each strand
        # read on its own column (`read_column`), both strands pooled; every exon when none has counts
        num = den = 0.0
        for free, col in ((fp, 0), (fn, 1)):
            sel = is_exon & free & single & (a_r > 0.0)
            num += float(cnt[sel, read_column(col, kappa)].sum())
            den += float(a_r[sel].sum())
        if den <= 0.0:
            sel = is_exon & (a_r > 0.0)
            num, den = float(cnt[sel].sum()), float(a_r[sel].sum())
        rho_rna = num / den if den > 0.0 else 0.0
        return _Library(rho, rho_rna, split_live)

    # ── phase 0: every node's own claim, and the recipient's rule per directed face ──────────────
    def prepare(self, ctx: BlockContext, library: _Library) -> "_PreparedTransfer":
        # a chain with no coarse intron has no factory rows: a factory with nothing to say, and the
        # rest of the layer — the faces, the lanes — exactly as with one
        src = ctx.factory_rows
        if src is None:
            src = np.zeros((int(ctx.n_slots), int(ctx.n_grid)))
        chain = _Chain(ctx, np.asarray(src, np.float64), self._strand)
        own = _claims(chain)
        # the recipient's rule per DIRECTED face, written into the face table by the builder that owns
        # that kind of face (the faces are disjoint: a second rule at a face is refused); the level lane
        # serves every face left without one
        faces = Faces(chain.lam, chain.left, chain.right)
        _splice_faces(chain, faces)
        _edge_level(chain, own, faces)
        _terminus_rules(chain, faces)
        _alternative_splice_site(chain, faces)
        # each lane exists iff its OWN coordinate does: a gDNA-free library has no gDNA lane and still
        # has its RNA lanes
        lanes: dict = {}
        gdna = gdna_lane(chain, own, faces, library.rho_gdna)
        if gdna is not None:
            lanes["gdna"] = gdna
        lanes.update(rna_lanes(chain, own, library))
        site = _SolveSite(chain.fp & chain.fn, {"pos": chain.fp, "neg": chain.fn})
        return _PreparedTransfer(own, faces, lanes, site)


# ══ THE CHAIN AS THE BUILDERS READ IT ═══════════════════════════════════════════════════════════════


class _Chain:
    """One sweep's chain as the builders read it: every array of the context the messages need,
    unpacked once — the node classes, the geometry, the counts and opportunities per node and per
    face, whether each node has own composition evidence — plus the factory's per-slot profiles ``src`` and the
    library's fitted strand model ``strand`` (``None``: every strand claim off)."""

    __slots__ = (
        "src",
        "strand",
        "n",
        "K",
        "lam",
        "is_bnd",
        "is_exon",
        "is_intron",
        "is_intergenic",
        "exon_of",
        "fp",
        "fn",
        "left",
        "right",
        "flags",
        "n_u",
        "n_s",
        "flux",
        "a_g",
        "a_r",
        "has_own_composition",
        "cnt",
        "belief",
        "route_rate",
        "sj_count",
    )

    def __init__(self, ctx: BlockContext, src: np.ndarray, strand):
        self.src, self.strand = src, strand
        self.n, self.K = ctx.n_slots, src.shape[1]
        self.lam = np.linspace(-float(ctx.logodds_window), float(ctx.logodds_window), self.K)
        self.is_bnd = np.asarray(ctx.is_boundary, bool)
        self.is_exon = np.asarray(ctx.is_exon_region, bool)
        self.fp, self.fn = np.asarray(ctx.free_pos, bool), np.asarray(ctx.free_neg, bool)
        self.is_intron = ~self.is_bnd & ~self.is_exon & (self.fp | self.fn)
        self.is_intergenic = ~self.is_bnd & ~self.is_exon & ~self.fp & ~self.fn
        #: the exon bit PER STRAND: a region that admits ``s`` and carries no exon of ``s`` is ``s``'s
        #: intron whatever the other strand does there
        self.exon_of = {
            "pos": np.asarray(ctx.exon_pos, bool),
            "neg": np.asarray(ctx.exon_neg, bool),
        }
        self.left, self.right = np.asarray(ctx.left, np.int64), np.asarray(ctx.right, np.int64)
        self.flags = np.asarray(ctx.boundary_flags, np.uint16)
        self.n_u = np.asarray(ctx.n_slot, np.float64)
        self.n_s = np.asarray(ctx.spliced_slot, np.float64)
        self.flux = np.asarray(ctx.sj_count, np.float64).sum(axis=1)
        self.a_g = np.asarray(ctx.eff_gdna, np.float64)
        self.a_r = np.asarray(ctx.eff_rna, np.float64)
        self.has_own_composition = np.asarray(ctx.has_own_composition, bool)
        self.cnt = np.asarray(ctx.unspliced_count, np.float64)
        self.belief = np.asarray(ctx.belief_fg, np.float64)
        #: per face by transcript strand, ``(low end, high end)``: the flux that crosses an exon's LEFT
        #: face belongs to the sj whose HIGH end is that boundary and vice versa, so one column serves
        #: both directions of a face
        self.route_rate = (
            np.asarray(ctx.route_rate_lo, np.float64),
            np.asarray(ctx.route_rate_hi, np.float64),
        )
        self.sj_count = (
            np.asarray(ctx.sj_count_lo, np.float64),
            np.asarray(ctx.sj_count_hi, np.float64),
        )

    def other_flank(self, b, e):
        return self.left[b] if self.right[b] == e else self.right[b]

    def intron_exon_pair(self, b):
        """``(exon, intron)`` across boundary ``b``, or ``(None, None)`` when it is not that face."""
        lo, hi_ = self.left[b], self.right[b]
        if lo < 0 or hi_ < 0:
            return None, None
        if self.is_exon[lo] and self.is_intron[hi_]:
            return lo, hi_
        if self.is_exon[hi_] and self.is_intron[lo]:
            return hi_, lo
        return None, None

    def strand_profile(self, x):
        """A slot's OWN strand profile: its unspliced count split by genome strand, the variance frozen
        at its incoming belief (a source-side read); the mode is data only."""
        kappa, od_g, od_r = self.strand
        f_ref = float(self.belief[x]) if np.isfinite(self.belief[x]) else 0.5
        return norm(
            strand_row_logodds(
                self.lam, self.cnt[x, 0], self.cnt[x, 1], bool(self.fp[x]), kappa, od_g, od_r, f_ref
            )
        )

    def strand_mode(self, y):
        """A node's own strand MODE where its channel is live: the gDNA share its unspliced split implies
        under the library's κ and that share's log-odds variance by the delta method — data only, the
        pair-discrepancy witness of the terminus rule and the alternative splice site. ``None`` at a
        vertex mode, which has no log-odds."""
        kappa = self.strand[0]
        nn = self.cnt[y].sum()
        p = self.cnt[y, 0] / nn
        ks = kappa if self.fp[y] else 1.0 - kappa
        f = (p - ks) / (0.5 - ks)
        if not 0.0 < f < 1.0:
            return None
        v_log = p * (1.0 - p) / nn / (p - ks) ** 2
        return f, v_log / (1.0 - f) ** 2


# ══ THE BUILDERS — one per shipped message ══════════════════════════════════════════════════════════


def _claims(c: _Chain) -> list:
    """Every node's OWN CLAIM — data only, never a belief: an intron's factory profile where the
    factory has one; an exon's own strand profile where the derived strand deadband declares the
    channel live; a single-strand boundary's own strand profile likewise (an AMBIG boundary's split
    constrains only the tilt, never the gDNA level — the Schur complement the local solve applies)."""
    own: list = [None] * c.n
    for i in np.flatnonzero(c.is_intron):
        if np.ptp(c.src[i]) > EPS:
            own[i] = norm(c.src[i])
    if c.strand is not None:
        live, single = c.has_own_composition, c.fp != c.fn
        for e in np.flatnonzero(c.is_exon & single & live):
            own[e] = c.strand_profile(e)
        for b in np.flatnonzero(c.is_bnd & live & single):
            own[b] = c.strand_profile(b)
    return own


def _splice_faces(c: _Chain, faces: Faces) -> None:
    """THE INTRON|EXON FACE. intron → boundary: FORWARD (one shared unspliced population); boundary →
    intron: FORWARD under a shared single strand. At a LICENSED face (no terminus, the same strand
    set, nothing depleted): boundary → exon by the splice-in face map — the certified flux capping the
    claimable gDNA share, widened by the face's counting — and exon → boundary by the splice-in map
    read backwards, marginalised over the face's spliced-to-unspliced ratio."""
    rr = tuple(r.sum(axis=1) for r in c.route_rate)
    sc = tuple(s.sum(axis=1) for s in c.sj_count)
    for b in np.flatnonzero(c.is_bnd):
        e, i = c.intron_exon_pair(b)
        if e is None:
            continue
        faces.set(i, b, FORWARD)
        if boundary_shares_strand(c.fp[b], c.fn[b], c.fp[i], c.fn[i]):
            faces.set(b, i, FORWARD)
        hi = 1 if c.left[e] == b else 0
        if not face_is_licensed(c.flags[b], c.fp[e], c.fn[e], c.fp[i], c.fn[i]):
            continue
        if not (c.n_u[b] > 0 and c.a_g[b] > 0 and c.a_r[b] > 0 and c.a_g[e] > 0 and c.a_r[e] > 0):
            continue  # a depleted face: silence
        le = face_map_lambda(
            c.lam, c.n_u[b], c.a_g[b], c.a_r[b], c.a_g[e], c.a_r[e], float(rr[hi][b])
        )
        s = float(sc[hi][b])
        faces.set(b, e, TRANSPORT, row=le, n_u=c.n_u[b], n_s=s)
        faces.set(e, b, SPLICE_OUT, n_u=c.n_u[b], n_s=s, a_b=c.a_g[b], a_x=c.a_g[e])


def _edge_level(c: _Chain, own: list, faces: Faces) -> None:
    """THE INTERGENIC|EXON EDGE → its exon: A LEVEL, one-sided. The edge's own claim is its gDNA
    COUNT — a marker profile; the rule reads the
    count, and nothing is held at an edge — and the exon converts it through its own total: the exon
    has at least the edge's gDNA density, at the count's own Poisson width; nothing above (no local
    witness prices a probed interior's enrichment over its edge), and a zero count is vacuous
    (darkness under capture is not absence)."""
    for e in np.flatnonzero(c.is_exon):
        for b in (c.left[e], c.right[e]):
            if b < 0 or not c.is_bnd[b]:
                continue
            o = c.other_flank(b, e)
            if not (
                o >= 0 and c.is_intergenic[o] and c.a_g[b] > 0 and c.a_g[e] > 0 and c.n_u[e] > 0
            ):
                continue
            own[b] = np.zeros(c.K)  # the level claim (the rule reads the count, not the row)
            faces.set(b, e, EDGE, row=edge_level_row(c.lam, c.n_u[b], c.n_u[e], c.a_g[b], c.a_g[e]))


def _terminus_rules(c: _Chain, faces: Faces) -> None:
    """THE TERMINUS BOUNDARY. With its OUTSIDE exon (exon|exon, a shared single strand): composition by
    the spliced-crossing licence — the outside exon's message travels the splice-out map, the
    boundary's the face map with the spliced density. Into the region INSIDE it (exon|exon and
    exon|intron alike): THE LEVEL RULE — the boundary's OWN profile through the level-kept map, blurred
    by the two totals' counting and the pair's own discrepancies (the totals' disagreement beyond
    counting; where both strand channels are live, the two strand modes' disagreement beyond counting);
    with no own claim, the crossing total's upper bound. At an sj+terminus boundary the junction's flux
    belongs to the flank on the junction's exon side: the outside map when that flank is outside, the
    inside's total when it is inside (the RNA joining there is measured; only the terminus's own
    transcription is not)."""
    for b in np.flatnonzero(c.is_bnd & (c.left >= 0) & (c.right >= 0)):
        lo, hi_ = c.left[b], c.right[b]
        o, i = outside_flank(c.flags[b], lo, hi_)
        if o is None:
            continue
        ex_side = junction_exon_side(c.flags[b], lo, hi_) if int(c.flags[b]) & SJ_FLAGS else None
        flux_b = float(c.flux[b]) if ex_side is not None else 0.0
        if (
            c.is_exon[lo]
            and c.is_exon[hi_]
            and boundary_shares_strand(c.fp[b], c.fn[b], c.fp[o], c.fn[o])
            and c.n_u[b] > 0
            and c.a_g[b] > 0
            and c.a_g[o] > 0
        ):
            s_out = c.n_s[b] + (flux_b if ex_side == o else 0.0)
            le5 = face_map_lambda(
                c.lam, c.n_u[b], c.a_g[b], c.a_g[b], c.a_g[o], c.a_g[o], s_out / c.a_g[b]
            )
            faces.set(o, b, SPLICE_OUT, n_u=c.n_u[b], n_s=s_out, a_b=c.a_g[b], a_x=c.a_g[o])
            faces.set(b, o, TRANSPORT, row=le5, n_u=c.n_u[b], n_s=s_out)
        # ── THE LEVEL RULE: the boundary's gDNA level into the inside region ─────────────────────
        if not (c.is_exon[i] and (c.is_exon[o] or c.is_intron[o])):
            continue
        if not (c.n_u[b] > 0 and c.n_u[i] > 0 and c.a_g[b] > 0 and c.a_g[i] > 0):
            continue
        if not boundary_shares_strand(c.fp[b], c.fn[b], c.fp[i], c.fn[i]):
            continue  # a strand change is the level lane's
        density_b = c.n_u[b] / c.a_g[b]
        total_b = c.n_u[b] + c.n_s[b]
        # (i) the totals per opportunity, beyond counting
        total_i = total_b + (flux_b if ex_side == i else 0.0)
        r = (c.n_u[i] / c.a_g[i]) / (total_i / c.a_g[b])
        v_pair = max(0.0, float(np.log(r)) ** 2 - (1.0 / c.n_u[i] + 1.0 / total_i))
        # (ii) the two strand modes where both channels are live: the inside's share predicted from
        # the crossing alone (the flux is not a crossing), against its own mode, beyond counting
        r_mode = (c.n_u[i] / c.a_g[i]) / (total_b / c.a_g[b])
        if c.strand is not None and c.has_own_composition[b] and c.has_own_composition[i]:
            mode_b = c.strand_mode(b)
            mode_i = None if mode_b is None else c.strand_mode(i)
            if mode_i is not None:
                (f_b, v_b), (f_i, v_i) = mode_b, mode_i
                f_pred = min(
                    f_b / r_mode, 1.0 - 1e-9
                )  # the level kept: the boundary's share over r
                d = np.log(f_i / (1.0 - f_i)) - np.log(f_pred / (1.0 - f_pred))
                v_pair += max(0.0, float(d * d) - (v_b + v_i + 1.0 / c.n_u[i] + 1.0 / total_b))
        v_level = float(count_logvar(c.n_u[b]) + count_logvar(c.n_u[i])) + v_pair
        # a level is made from the sender's MEASUREMENT only: its own profile through the level-kept
        # map; with none, the crossing total's upper bound. What it holds — an imputation — never
        # crosses a level face (`Faces.apply`, LEVEL)
        faces.set(
            b,
            i,
            LEVEL,
            row=level_map_lambda(c.lam, density_b, c.a_g[i], c.n_u[i]),
            row2=level_bound_row(c.lam, density_b, c.a_g[i], c.n_u[i], v_level),
            var=v_level,
        )


def _alternative_splice_site(c: _Chain, faces: Faces) -> None:
    """THE ALTERNATIVE SPLICE SITE — an exon|exon boundary carrying a junction and no terminus — and
    both its flanks, both ways: the intron-side flank shares the full unspliced crossing, the
    exon-of-both flank the crossing plus the isoform that splices out here, measured as the route
    flux; each pair widened by its own strand modes' disagreement beyond counting (the discrepancy
    rule PER PAIR, nothing pooled). Needs the strand model: the width is a strand-mode witness."""
    if c.strand is None:
        return
    for b in np.flatnonzero(c.is_bnd & (c.left >= 0) & (c.right >= 0)):
        lo, hi_ = c.left[b], c.right[b]
        if not (c.is_exon[lo] and c.is_exon[hi_]):
            continue
        c_side, e_side = junction_flanks(c.flags[b], lo, hi_)
        if c_side is None or not (c.n_u[b] > 0 and c.a_g[b] > 0):
            continue
        for x, s_out in ((e_side, c.n_s[b] + c.flux[b]), (c_side, c.n_s[b])):
            if not (boundary_shares_strand(c.fp[b], c.fn[b], c.fp[x], c.fn[x]) and c.a_g[x] > 0):
                continue
            width = 0.0
            if c.has_own_composition[b] and c.has_own_composition[x]:
                mode_b = c.strand_mode(b)
                mode_x = None if mode_b is None else c.strand_mode(x)
                if mode_x is not None:
                    (f_b, v_b), (f_x, v_x) = mode_b, mode_x
                    lo_b, lo_x = np.log(f_b / (1.0 - f_b)), np.log(f_x / (1.0 - f_x))
                    v_ratio = s_out / (c.n_u[b] * (c.n_u[b] + s_out)) if s_out > 0 else 0.0
                    d = lo_b - lo_x - np.log((c.n_u[b] + s_out) / c.n_u[b])
                    width = max(0.0, d * d - (v_b + v_x + v_ratio))
            le7 = face_map_lambda(
                c.lam, c.n_u[b], c.a_g[b], c.a_g[b], c.a_g[x], c.a_g[x], s_out / c.a_g[b]
            )
            faces.set(
                x, b, SPLICE_OUT, n_u=c.n_u[b], n_s=s_out, a_b=c.a_g[b], a_x=c.a_g[x], width=width
            )
            faces.set(b, x, TRANSPORT, row=le7, n_u=c.n_u[b], n_s=s_out, width=width)


class _SolveSite:
    """What the solve's two RNA deliveries read at a node: the AMBIG mask and each strand's admitting
    nodes."""

    __slots__ = ("ambig", "free")

    def __init__(self, ambig, free):
        self.ambig, self.free = ambig, free


class _PreparedTransfer:
    """One sweep's working object: every node's own claim, the rules per directed face, the lanes by
    population (``"gdna"``, ``"pos"``, ``"neg"`` — any may be absent) and the site. The passes' state
    is the backbone's table, never a copy here."""

    def __init__(self, own, faces: "Faces | None", lanes: dict | None = None, site=None):
        self.own = own
        self.faces = faces
        self.lanes: dict = {} if lanes is None else dict(lanes)
        self.site = site

    # ── phase 1: propagate — the recipient's kernel, run by the backbone in chain order ──────────
    def propagate(self, received: Received, *, backward: bool):
        faces = self.faces
        if self.own is None or not (
            faces.any() or any(ln.face.any() for ln in self.lanes.values())
        ):
            return None  # nothing to say anywhere: every node holds silence
        # each lane's faces, its table and its emptiness, unpacked once rather than per node
        lanes = tuple(
            (ln, ln.face, getattr(received, ln.field), ln.empty.tolist())
            for ln in self.lanes.values()
        )
        kind, apply, own = faces.kind, faces.apply, self.own
        composition, has_composition = received.composition, received.has_composition

        def receive(s: int, i: int):
            """``s`` sends two things apart — its own claim and what it holds from its far side (its own
            row of this pass's table, written one step earlier) — and the face decides: a composition
            rule composes and maps (MODIFY), passes (FORWARD) or reads the measurement only (a level
            face); each lane whose faces include this one carries its population's level; absent all (a
            structural pure-gDNA neighbour) STOP, and row ``i`` stays silent. A face is ``(i, side)``:
            the side of ``i`` that ``s`` is on."""
            side = 0 if s < i else 1
            if kind[i, side]:
                far = composition[s] if has_composition[s] else None
                out = apply(s, i, own[s], far)
                if out is not None and out.max() - out.min() > EPS:
                    composition[i] = norm(out)
                    has_composition[i] = True
            for lane, face, levels, empty in lanes:
                if face[i, side] and lane.emit(s, i, levels) and not empty[i]:
                    lane.receive(levels, s, i)

        return receive

    # ── phase 2: solve — the policy's half ───────────────────────────────────────────────────────
    def solve(self, from_left: Received, from_right: Received) -> PsiMessage:
        if self.own is None:
            return PsiMessage.silent()
        # the two held compositions add (independent witnesses about one slot), in table order: left,
        # right — an absent one adds nothing
        rows = np.where(from_left.has_composition[:, None], from_left.composition, 0.0) + np.where(
            from_right.has_composition[:, None], from_right.composition, 0.0
        )
        fused = from_left.has_composition | from_right.has_composition
        gdna = self.lanes.get("gdna")
        if gdna is not None:
            # a held gDNA level read as the node's composition row through its own total; two bounds on
            # one density intersect (the tighter wins) and join the compositions as one more witness
            held = (from_left.level_gdna, from_right.level_gdna)
            for i in np.flatnonzero((held[0].present | held[1].present) & ~gdna.empty):
                bound = intersect([gdna.row(lv.profile[i], i) for lv in held if lv.present[i]])
                rows[i] = rows[i] + bound if fused[i] else bound
                fused[i] = True
        rows[fused] -= rows[fused].max(axis=1, keepdims=True)  # `fuse`: add, then re-normalise
        live = bool(fused.any())
        live = self._ceilings(from_left, from_right, rows) or live
        cube = self._cube_rows(from_left, from_right)
        if not live and not cube:
            return PsiMessage.silent()
        return PsiMessage(lam_rows=rows if live else None, cube_rows=cube or None)

    def _ceilings(self, from_left, from_right, rows) -> bool:
        """THE UPPER SIDE AT SINGLE-STRAND NODES: an RNA level of the node's live strand says "at most
        this much gDNA". ⛔ Read ONLY from a face that sent no composition — a held level on that
        side, and the node's own junction flux at that face — because a licensed face's splice-in map
        already carries the flux as its cap and a composition already carries its sender's witnesses,
        so reading them again counts the same evidence twice. Bounds intersect; the row joins the
        node's other witnesses. Returns whether anything was added."""
        site = self.site
        if site is None:
            return False
        added = False
        for name in _RNA:
            rl = self.lanes.get(name)
            if rl is None:
                continue
            for i in np.flatnonzero(site.free[name] & ~site.ambig & ~rl.empty):
                bounds = []
                for side, t in ((0, from_left), (1, from_right)):
                    if not t.has_neighbour[i] or t.has_composition[i]:
                        continue
                    lv = getattr(t, rl.field)
                    if lv.present[i]:
                        bounds.append(lv.profile[i])
                    fx = rl.flux_at(i, side)
                    if fx is not None:
                        bounds.append(fx)
                if not bounds:
                    continue
                row = rl.row(intersect(bounds), i)
                if np.ptp(row) <= EPS:
                    continue
                rows[i] = fuse([rows[i], row]) if np.ptp(rows[i]) > EPS else row
                added = True
        return added

    def _cube_rows(self, from_left, from_right) -> dict:
        """THE DELIVERY AT AMBIG NODES: the held RNA levels — both sides intersected, plus the node's
        OWN flux level (the spliced claim's one hop, boundary → exon, read at the exon; an AMBIG node
        has no own strand claim, so its own level is the flux alone) — as a :class:`CubeRow`, the row's
        ingredients, which ψ evaluates at its own θ nodes. The tilt needs no lane of its own: both
        strands' bounds constrain it through the shares."""
        site = self.site
        pos, neg = self.lanes.get("pos"), self.lanes.get("neg")
        if site is None or pos is None or neg is None:
            return {}
        out = {}
        for i in np.flatnonzero(site.ambig & ~pos.empty):
            profiles = {}
            for rl in (pos, neg):
                bounds = [
                    getattr(t, rl.field).profile[i]
                    for t in (from_left, from_right)
                    if getattr(t, rl.field).present[i]
                ]
                if rl.own_level[i] is not None:
                    bounds.append(lower_side(rl.own_level[i]))
                if bounds:
                    profiles[rl.population] = intersect(bounds)
            if profiles:
                out[int(i)] = CubeRow(
                    profile_pos=profiles.get("pos"),
                    profile_neg=profiles.get("neg"),
                    u=pos.u,
                    total=float(pos.total[i]),
                    opportunity=float(pos.a[i]),
                    rho_ref=float(pos.rho_ref),
                )
        return out
