"""TransferPolicy — the COMPOSITION TRANSFER policy on the two-phase backbone: every node states its
own claim, every directed face carries a recipient's rule, and the two passes do the rest.

       Gate: ``tests/calibration/test_transfer_policy.py``

Every message is a COMPOSITION PROFILE — a max-normalised log-likelihood over the solve grid of the
destination's gDNA share — carried across ONE face by a derived map from `transfer_rows`, or a
population's LEVEL — an absolute profile over the log density — carried where composition cannot
cross, so no level ever crosses a capture cliff and no constant is anywhere. The policy has three
parts, and ``prepare`` is a table of contents: one named BUILDER per shipped message.

* Every node's OWN CLAIM (`claims`): an intron's factory profile (its density against the
  intergenic background — the context's ``factory_rows``, the very array ψ adds as its λ-factor); an
  exon's or a boundary's own strand profile
  where the library's strand protocol decision declares the channel live (``strand``); an
  intergenic|exon edge's gDNA COUNT (the level lane: the edge's crossing is structurally pure gDNA).
  ⛔ A claim is data only — never a belief, which already holds the prior and the neighbours.
* The RECIPIENT's rule per directed face: absent = STOP (composition cannot cross: the recipient
  holds silence); the identity = FORWARD (the two objects share one unspliced population exactly); a
  map = MODIFY (the face's arithmetic with its counting width and, where two witnesses exist, the
  pair's own discrepancy). The rules are TYPED TABLES (`faces.Faces`): every directed face is one of a
  node's two sides — it hears from its left neighbour or its right — so a rule is a KIND and its
  parameters at ``(destination, side)``, five kinds in all, and the passes read the table. The rules
  ARE the shipped messages, each built by the builder named — a function of `native/prepare_kernel.cpp`,
  which builds a whole block's claims, rules and lanes in one call (`native.transfer_prepare`) into the
  tables ``prepare`` allocates:

  - `splice_faces` — the intron|exon face. intron ⇄ boundary: FORWARD both ways (one shared
    unspliced population; into the intron only under a shared single strand). boundary → exon at a
    LICENSED face (no terminus, the same strand set): the splice-in face map, the certified flux
    capping the claimable gDNA share, widened by the face's counting. exon → boundary at that face:
    the splice-in map read backwards, marginalised over the face's spliced-to-unspliced ratio.
  - `edge_level` — intergenic|exon edge → exon: THE EDGE'S LEVEL, one-sided: the exon has at least
    the edge's gDNA density, at the count's Poisson width; nothing above (no local witness prices
    capture's enrichment of the interior), a zero count vacuous (darkness is not absence).
  - `terminus_rules` — the exon|exon TERMINUS boundary ⇄ its OUTSIDE exon: the licence counts the
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
  - `alternative_splice_site` — the ALTERNATIVE SPLICE SITE ⇄ both flanks: the intron-side flank
    shares the full unspliced crossing, the exon-of-both flank the crossing plus the leaving isoform
    measured as the route flux; each pair widened by its own disagreement beyond counting, nothing
    pooled.
  - `gdna_lane` — EVERY OTHER DIRECTED FACE (a strand change, termini pointing both ways, the AMBIG
    complex, and every face into or out of an EMPTY node — a piece with no total) carries THE LEVEL
    LANE: the gDNA level as an ABSOLUTE profile over the log density, which needs no map and no
    recipient. An empty node forwards it unchanged; a full node emits the intersection of its own
    level's lower side and the level it holds; the recipient prices the hop (both totals' counting
    plus the abundance discrepancy beyond it) and takes it as a LOWER bound — a level that crosses a
    face says "at least this much gDNA" and nothing more. The only faces with no rule at all lead
    into intergenic regions (structurally pure gDNA) or off the chain.
  - `rna_lane` — THE RNA LEVEL LANES: one lane per strand, its FACES from the flag bits (strand
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
stays a no-claim — a flat profile or an absent factory row is no claim on that channel, never a
zero-filled one; a message is built from the source's claim and the recipient's constants
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

from ...native import transfer_pass, transfer_prepare
from ..simplex_logodds import CubeRow
from . import BlockContext, ChainView, PsiMessage, Received
from .faces import Faces, RowTable, fuse, norm
from .lanes import LevelLane
from .transfer_rows import _MARGINAL_NODES, EPS, intersect, lower_side, read_column

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
    library's protocol preserves strand and some single-strand exon has counts."""

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
        # the protocol decision (`ChainView.strand_live`): where the protocol preserves strand, every
        # counted single-strand exon's strand precision is positive; where it does not, none is (an intron's
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
        return _prepare(_Chain(ctx, np.asarray(src, np.float64), self._strand), library)


def _prepare(c: _Chain, library: _Library) -> "_PreparedTransfer":
    """THE BUILDERS IN ONE NATIVE CALL (`native.transfer_prepare`): the tables allocated here — the
    claims, the face tables, each lane's arrays in the layout the pass reads — and written there. The lanes'
    coordinates, witnesses and emptiness are the chain's arrays; the kernel writes the faces, the own levels,
    the flux levels and the flux witnesses; each lane exists iff its OWN coordinate does (a gDNA-free library
    has no gDNA lane and still has its RNA lanes). Gates: the transfer gates hold the tables to independent
    recomputes (``test_transfer_faces.py``, ``test_transfer_policy.py``, ``test_transfer_rna_lanes.py``)."""
    n, K, lam = c.n, c.K, c.lam
    own = RowTable(n, K)
    faces = Faces(lam, c.left, c.right)
    kappa = None if c.strand is None else float(c.strand[0])
    split_live = library.split_live

    def lane(name, rho_ref, count, a, other=None):
        empty = ~(c.n_u > 0.0) | ~(a > 0.0)
        return LevelLane(
            name,
            lam,
            lam,
            rho_ref,
            count,
            a,
            empty,
            RowTable(n, K),
            np.zeros((n, 2), bool),
            two_sided=np.zeros((n, 2), bool),
            total=c.n_u,
            flux=RowTable((n, 2), K),
            other=other,
            flux_witness=RowTable(n, 2),
        )

    def tables(ln):
        return (
            ln.face,
            ln.two_sided,
            ln.own_level.rows,
            ln.own_level.mask,
            ln.flux.rows,
            ln.flux.mask,
            ln.flux_witness.rows,
            ln.flux_witness.mask,
        )

    lanes: dict = {}
    if library.rho_gdna > 0.0:
        lanes["gdna"] = lane("gdna", library.rho_gdna, c.n_u, c.a_g)
    for name, col in (("pos", 0), ("neg", 1)):
        col_read = read_column(col, kappa)
        lanes[name] = lane(
            name,
            library.rho_rna,
            np.ascontiguousarray(c.cnt[:, col_read]),
            c.a_r,
            other=np.ascontiguousarray(c.cnt[:, 1 - col_read]) if split_live else None,
        )
    strand = (False, 0.0, 0.0, 0.0) if c.strand is None else (True, *c.strand)
    faces.n_rows = transfer_prepare(
        lam=lam,
        is_boundary=c.is_bnd,
        is_exon=c.is_exon,
        free_pos=c.fp,
        free_neg=c.fn,
        exon_pos=c.exon_of["pos"],
        exon_neg=c.exon_of["neg"],
        left=c.left,
        right=c.right,
        flags=c.flags,
        n_u=c.n_u,
        n_s=c.n_s,
        a_g=c.a_g,
        a_r=c.a_r,
        cnt=c.cnt,
        belief=c.belief,
        has_own_composition=c.has_own_composition,
        flux=c.flux,
        route_rate_lo=c.route_rate[0],
        route_rate_hi=c.route_rate[1],
        sj_count_lo=c.sj_count[0],
        sj_count_hi=c.sj_count[1],
        src=c.src,
        has_strand=strand[0],
        kappa=strand[1],
        od_g=strand[2],
        od_r=strand[3],
        rho_gdna=float(library.rho_gdna),
        rho_rna=float(library.rho_rna),
        own=own.rows,
        own_mask=own.mask,
        faces=(
            faces.kind,
            faces.row,
            faces.row2,
            faces.n_u,
            faces.n_s,
            faces.a_b,
            faces.a_x,
            faces.width,
            faces.var,
            faces.rows,
        ),
        gdna=None if "gdna" not in lanes else tables(lanes["gdna"]),
        pos=tables(lanes["pos"]),
        neg=tables(lanes["neg"]),
    )
    site = _SolveSite(c.fp & c.fn, {"pos": c.fp, "neg": c.fn})
    return _PreparedTransfer(own, faces, lanes, site)


# ══ THE CHAIN AS THE BUILDERS READ IT ═══════════════════════════════════════════════════════════════


class _Chain:
    """One sweep's chain as the builders read it: every array of the context the messages need,
    unpacked once and handed to `native.transfer_prepare` by name — the node classes, the geometry, the
    counts and opportunities per node and per face, whether each node has own composition evidence — plus
    the factory's per-slot profiles ``src`` and the library's fitted strand model ``strand`` (``None``:
    every strand claim off)."""

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


class _SolveSite:
    """What the solve's two RNA deliveries read at a node: the AMBIG mask and each strand's admitting
    nodes."""

    __slots__ = ("ambig", "free")

    def __init__(self, ambig, free):
        self.ambig, self.free = ambig, free


class _PreparedTransfer:
    """One sweep's working object: every node's own claim (a `RowTable`), the rules per directed face
    (`Faces`), the lanes by population (``"gdna"``, ``"pos"``, ``"neg"`` — any may be absent) and the
    site. The passes' state is the backbone's table, never a copy here."""

    def __init__(self, own, faces: "Faces | None", lanes: dict | None = None, site=None):
        self.own = own
        self.faces = faces
        self.lanes: dict = {} if lanes is None else dict(lanes)
        self.site = site

    # ── phase 1, natively: the whole pass on the table in one call ─────────────────────────────────
    def tables(self) -> dict:
        """The tables the native pass reads, by the kernel's argument names — the builders' own, shared
        by both passes: every node's own claim with its mask, the faces' arrays and their row store,
        each lane's arrays (`LevelLane.tables`). Nothing is copied."""
        faces = self.faces
        return dict(
            lam=faces.lam,
            own=self.own.rows,
            own_mask=self.own.mask,
            f_kind=faces.kind,
            f_row=faces.row,
            f_row2=faces.row2,
            f_n_u=faces.n_u,
            f_n_s=faces.n_s,
            f_a_b=faces.a_b,
            f_a_x=faces.a_x,
            f_width=faces.width,
            f_var=faces.var,
            f_rows=faces.rows[: faces.n_rows],
            marginal_nodes=_MARGINAL_NODES,
            lanes=[ln.tables() for ln in self.lanes.values()],
        )

    def run_pass(self, received: Received, seq, nbr, terminal, *, backward: bool) -> None:
        """PHASE 1 on the table in one native call (`native.transfer_pass`): the same hops, in the same
        order, as ``propagate``'s kernel run by the backbone; nothing to say anywhere is a no-op. Gate:
        ``tests/calibration/test_pass_kernel.py``."""
        faces = self.faces
        if self.own is None or not (
            faces.any() or any(ln.face.any() for ln in self.lanes.values())
        ):
            return
        levels = [
            (
                lv.present,
                lv.profile,
                lv.count,
                lv.opportunity,
                lv.has_witness,
                lv.rna_count,
                lv.rna_count_var,
            )
            for lv in (received.level_gdna, received.level_rna_pos, received.level_rna_neg)
        ]
        transfer_pass(
            seq=np.ascontiguousarray(seq, np.int64),
            nbr=np.ascontiguousarray(nbr, np.int64),
            terminal=np.ascontiguousarray(terminal, bool),
            has_composition=received.has_composition,
            composition=received.composition,
            levels=levels,
            **self.tables(),
        )

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
