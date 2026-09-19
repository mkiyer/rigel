"""TransferPolicy — the COMPOSITION TRANSFER policy on the two-phase backbone: every node states its
own claim, every directed face carries a recipient's rule, and the two passes do the rest.

       Gate: ``tests/calibration/test_transfer_policy.py``

Every message is a COMPOSITION PROFILE — a max-normalised log-likelihood over the solve grid of the
destination's gDNA share — carried across ONE face by a derived map (the row constructors of
`native/transfer_rows.h`), or a population's LEVEL — an absolute profile over the log density — carried
where composition cannot cross, so no level ever crosses a capture cliff and no constant is anywhere. The
policy has three parts, all of them the kernel's (`native/transfer_kernel.h`, run per block by
`native/solve_kernel.cpp`); ``prepare_block`` is a table of contents: one named BUILDER per shipped
message.

* Every node's OWN CLAIM (``claims``): an intron's factory profile (its density against the
  intergenic background — the factory's row for it, the very row ψ adds as its λ-factor); an exon's or
  a boundary's own strand profile where the library's strand protocol decision declares the channel
  live (``strand``); an intergenic|exon edge's gDNA COUNT (the level lane: the edge's crossing is
  structurally pure gDNA). ⛔ A claim is data only — never a belief, which already holds the prior and
  the neighbours.
* The RECIPIENT's rule per directed face: absent = STOP (composition cannot cross: the recipient
  holds silence); the identity = FORWARD (the two objects share one unspliced population exactly); a
  map = MODIFY (the face's arithmetic with its counting width and, where two witnesses exist, the
  pair's own discrepancy). The rules are TYPED TABLES: every directed face is one of a node's two
  sides — it hears from its left neighbour or its right — so a rule is a KIND and its parameters at
  ``(destination, side)``, five kinds in all, and the passes read the table. The rules ARE the shipped
  messages, each built by the builder named:

  - ``splice_faces`` — the intron|exon face. intron ⇄ boundary: FORWARD both ways (one shared
    unspliced population; into the intron only under a shared single strand). boundary → exon at a
    LICENSED face (no terminus, the same strand set): the splice-in face map, the certified flux
    capping the claimable gDNA share, widened by the face's counting. exon → boundary at that face:
    the splice-in map read backwards, marginalised over the face's spliced-to-unspliced ratio.
  - ``edge_level`` — intergenic|exon edge → exon: THE EDGE'S LEVEL, one-sided: the exon has at least
    the edge's gDNA density, at the count's Poisson width; nothing above (no local witness prices
    capture's enrichment of the interior), a zero count vacuous (darkness is not absence).
  - ``terminus_rules`` — the exon|exon TERMINUS boundary ⇄ its OUTSIDE exon: the licence counts the
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
  - ``alternative_splice_site`` — the ALTERNATIVE SPLICE SITE ⇄ both flanks: the intron-side flank
    shares the full unspliced crossing, the exon-of-both flank the crossing plus the leaving isoform
    measured as the route flux; each pair widened by its own disagreement beyond counting, nothing
    pooled.
  - ``gdna_lane`` — EVERY OTHER DIRECTED FACE (a strand change, termini pointing both ways, the AMBIG
    complex, and every face into or out of an EMPTY node — a piece with no total) carries THE LEVEL
    LANE: the gDNA level as an ABSOLUTE profile over the log density, which needs no map and no
    recipient. An empty node forwards it unchanged; a full node emits the intersection of its own
    level's lower side and the level it holds; the recipient prices the hop (both totals' counting
    plus the abundance discrepancy beyond it) and takes it as a LOWER bound — a level that crosses a
    face says "at least this much gDNA" and nothing more. The only faces with no rule at all lead
    into intergenic regions (structurally pure gDNA) or off the chain.
  - ``rna_lane`` — THE RNA LEVEL LANES: one lane per strand, its FACES from the flag bits (strand
    ``s`` crosses a face iff the boundary carries none of ``s``'s bits and both nodes admit ``s``;
    across ``s``'s own junction it enters ``s``'s INTRON — the intron test is per strand, from the
    chain's ``exon_pos`` / ``exon_neg`` — and not ``s``'s exon; a terminus of ``s`` stops ``s``),
    TWO-SIDED only between an intron of ``s`` and its own boundary (one shared unspliced population:
    the whole profile), lower-only everywhere else, and EVERY hop priced by the pair — both counts'
    counting plus the disagreement between the two nodes' estimates of the strand's abundance.
    SOURCES: a single-strand node's own claim read as its live strand's RNA level; the certified flux
    at each of an exon's junctions as that strand's level at the exon (one hop, boundary → exon; two
    junctions pay their pair's disagreement beyond counting) — at an EMPTY exon piece too, which emits
    it with the flux's own witness. DELIVERED at AMBIG nodes as one row over ψ's (λ, θ) cube
    (`simplex_logodds.CubeRows`): a lower bound on RNA+ is an upper bound on the gDNA share through
    the node's own strand counts — the side the gDNA lane cannot give (the bracket theorem, gated);
    and at SINGLE-STRAND nodes as a CEILING on the gDNA share, read only from a face that sent no
    composition (the solve). The tilt needs no lane of its own.

* The two passes and the solve (``pass_block``, ``solve_block``): a node SENDS two things apart — its
  own claim (a measurement) and what it holds from its far side (an imputation) — and the recipient's
  rule decides what to do with each: a composition rule composes them (a witness product: profiles
  add) and maps the product; a level rule reads the measurement only, because an imputation is never
  re-issued as a level; a lane intersects them as bounds; either may stop — so a claim travels as far
  as the faces admit it, each hop charging its own counting width, and no node ever hears its own
  claim back (the forward pass composes only what came from the left, the backward pass only what came
  from the right). At the solve the two held profiles add, the two held levels intersect and join them
  as a constraint, and ψ fuses the row with the slot's own evidence and the prior.

The laws the policy keeps: the sender publishes its claim unchanged; the recipient decides; a no-claim
stays a no-claim — a flat profile or an absent factory row is no claim on that channel, never a
zero-filled one; a message is built from the source's claim and the recipient's constants
and observations, never the recipient's belief.

THE LIBRARY (`TransferPolicy.library`, once per sweep over the whole chain, `_Library`): the three
level lanes' coordinates — the library's structurally pure gDNA density, and each strand's unspliced
density over its single-strand exons — and whether the strand split is a live RNA witness anywhere.
They are ratios of sums over structurally selected slots and one boolean, read from observations and
geometry only, and they are the ONLY things a message knows about slots outside its own block. The
gates hold the kernel's tables to independent recomputes through `native.transfer_prepare`,
`transfer_pass` and `transfer_solve` (`tests/calibration/_transfer_harness.py`).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ...native import transfer_rows as _rows
from . import ChainView

__all__ = ["TransferPolicy"]


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
    profiles are the kernel's own rows; a factory with nothing to claim leaves the faces and the lanes
    exactly as with one."""

    name = "transfer"

    def __init__(self, strand: tuple[float, float, float] | None = None):
        self.strand = None if strand is None else tuple(float(x) for x in strand)

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
            den = float(a_g[a_g > 0.0].sum())
            rho = float(n_u[a_g > 0.0].sum()) / den if den > 0.0 else 0.0
        kappa = None if self.strand is None else float(self.strand[0])
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
        # read on its own column (the kernel's `read_column`: its own under a sense protocol, the other
        # under an antisense one), both strands pooled; every exon when none has counts
        num = den = 0.0
        for free, col in ((fp, 0), (fn, 1)):
            sel = is_exon & free & single & (a_r > 0.0)
            col_read = int(
                _rows.read_column(col, kappa is not None, 0.0 if kappa is None else kappa)
            )
            num += float(cnt[sel, col_read].sum())
            den += float(a_r[sel].sum())
        if den <= 0.0:
            sel = is_exon & (a_r > 0.0)
            num, den = float(cnt[sel].sum()), float(a_r[sel].sum())
        rho_rna = num / den if den > 0.0 else 0.0
        return _Library(rho, rho_rna, split_live)
