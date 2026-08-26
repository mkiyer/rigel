"""rigel.calibration.structural_claims — the first pass's substrate, derived from structure alone.

       Gate: ``tests/calibration/test_structural_claims.py``

Classifies every chain slot by what the ANNOTATION alone constrains its UNSPLICED population to be —
no counts, no solver, a pure function of the statics recomputed at calibration init (an owner ruling:
per-BOUNDARY structural flags were deliberately dropped from the index schema, and a persisted
calibration-facing artifact would need a reach digest no existing hash provides).

Four disjoint classes, each carrying a CLAIM that certified per-slot truth can test with no solver —
which is what makes the substrate falsifiable rather than a guess:

===================  ===============================================================================
``intergenic``       no RNA strand admissible (:func:`~rigel.calibration.region_geometry.g1_locked`,
                     both kinds): the composition is certain, truth must have no RNA fragments at all
``ss_intron_region``   exactly one strand admissible and no exon membership: contiguous mature RNA
                     does not fit inside, so truth must have no mature fragments — what remains is
                     the deconvolution problem (gDNA plus not-yet-spliced RNA)
``ss_intron_boundary`` exactly one strand continuous and no contiguous exon across: an unspliced
                     crossing has no mature term, same claim one dimension down
``solvable_exon``    single-stranded exonic REGION with at least one flanking BOUNDARY that is a
                     splice site with no contiguous exon across it. The claim is the FLANK's: its
                     unspliced crossings have no mature term, and its sj flux cannot be gDNA — so the
                     exon's own gDNA level is reachable without imputing anything from far away
===================  ===============================================================================

⭐ The set is a function of the two admissibility bits plus the mature-crossing bit, so ``|T| ∈
{1,2,3}`` always and no fourth population can be expressed here. An AMBIG slot (both strands
admissible) appears in NO class: on unstranded data it has neither a strand channel nor a
structurally clean neighbour, and excluding it from the TRAINING substrate is what keeps the prior
honest — it does not mean those slots go unanswered later.

⚠ **"solvable" here is STRUCTURAL and is not the sweep's runtime ``solvable`` mask** (``(free_pos |
free_neg) & (count > 0)`` — admits RNA and has data, nearly the opposite population). The two share a
word because both answer "can this slot be solved?", but this one answers it from the annotation and
that one from the data.

⚠ The mature-crossing bit at a BOUNDARY is the per-strand AND over its two flanks' exon membership —
a MAY-cross over-approximation, so a boundary where two different transcripts' exons abut reads
active even though no single transcript is contiguous there. That direction only UNDER-admits (a
claimed slot's claim still holds); the converse cannot happen, because a transcript contiguous across
a boundary puts its exon on both flanks.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..types import Strand
from .region_chain import BOUNDARY, REGION, RegionChain
from .region_geometry import RegionStatics, g1_locked
from .splice_graph import (
    FLAG_TES_NEG,
    FLAG_TES_POS,
    FLAG_TSS_NEG,
    FLAG_TSS_POS,
    is_splice_site,
)

__all__ = ["StructuralClaims", "build_structural_claims"]


@dataclass(frozen=True, slots=True)
class StructuralClaims:
    """Per-slot membership of the four structural classes (all arrays length ``n_slots``, bool).

    ``exon_flank_left`` / ``exon_flank_right`` name WHICH flank licenses a ``solvable_exon`` slot
    (both may be true); the flank slot itself is ``chain.left`` / ``chain.right``. They are exactly
    the flanks whose claim `structural_claims_audit.py` scores against certified truth.
    """

    n_slots: int
    intergenic: np.ndarray
    ss_intron_region: np.ndarray
    ss_intron_boundary: np.ndarray
    solvable_exon: np.ndarray
    exon_flank_left: np.ndarray
    exon_flank_right: np.ndarray
    #: ⭐ the COMPLETENESS bit, per licensed flank (a strict subset of ``exon_flank_*``): True where the
    #: flank's account of the exon's RNA is COMPLETE — no transcript-end bit faces the exon, so every
    #: molecule overlapping the exon passes through this flank (the theorem is in the builder). A
    #: complete flank's transfer is a TWO-SIDED estimate (estimate-grade for landscape training); an
    #: incomplete one licenses only the ONE-SIDED at-least-this-much-RNA bound (bound-grade).
    exon_flank_left_complete: np.ndarray
    exon_flank_right_complete: np.ndarray

    @property
    def claimed(self) -> np.ndarray:
        """The substrate: the union of the four classes. Everything else is not this pass's business."""
        return (
            self.intergenic | self.ss_intron_region | self.ss_intron_boundary | self.solvable_exon
        )


def interface_masks(claims: "StructuralClaims") -> tuple:
    """The three stage-0 masks the message layer consumes, under POLICY-NEUTRAL names — the
    backbone carries them without knowing the flank concept (its vocabulary firewall,
    `tests/calibration/test_sweep_backbone.py`): per-exon LEFT/RIGHT interface certification
    (every route arriving at that interface is certified, and no terminus admits molecules the
    flux cannot see), and the claimed ss-intron boundary class."""
    import numpy as np

    return (
        np.asarray(claims.exon_flank_left_complete, bool),
        np.asarray(claims.exon_flank_right_complete, bool),
        np.asarray(claims.ss_intron_boundary, bool),
    )


def build_structural_claims(chain: RegionChain, statics: RegionStatics) -> StructuralClaims:
    """Derive the four classes from the statics — O(n_slots) array math, nothing read but structure."""
    kind = np.asarray(chain.kind)
    is_region = kind == REGION
    is_boundary = kind == BOUNDARY
    fp = np.asarray(statics.free_pos, bool)
    fn = np.asarray(statics.free_neg, bool)
    mrna = np.asarray(statics.mrna_active_pos, bool) | np.asarray(statics.mrna_active_neg, bool)

    intergenic = g1_locked(fp, fn)
    single_stranded = fp ^ fn
    ss_intron_region = single_stranded & is_region & ~mrna
    ss_intron_boundary = single_stranded & is_boundary & ~mrna

    sj = is_splice_site(statics.boundary_flags, Strand.POS) | is_splice_site(
        statics.boundary_flags, Strand.NEG
    )
    flank_qualifies = is_boundary & sj & ~mrna
    # gather through the chain's adjacency: a -1 terminal indexes the appended sentinel, which is False
    at = np.concatenate([flank_qualifies, [False]])
    exonic = single_stranded & is_region & mrna
    left = np.asarray(chain.left)
    right = np.asarray(chain.right)
    exon_flank_left = exonic & at[left]
    exon_flank_right = exonic & at[right]

    # ── THE COMPLETENESS THEOREM ─────────────────────────────────────────────────────────────────────
    # A molecule overlapping the exon either CONTAINS the licensed flank (crosses it contiguously or
    # splices through it — both measured AT the flank) or ORIGINATES exactly at it: the exon's interior
    # holds no boundary, so a molecule that misses the flank must have its end there, and a transcript
    # end is a graph bit. So a flank's account of the exon's RNA is COMPLETE ⇔ the flank carries no
    # transcript-end bit FACING the exon: a genomic-LOW end at a LEFT flank, a genomic-HIGH end at a
    # RIGHT one. ⭐ The both-strand masks below are EXACTLY the strand-specific check at a licensed
    # flank: an other-strand end in the facing set would extend INTO the exon, making it AMBIG and
    # therefore unlicensed — so no strand plumbing is needed, and no test can distinguish the forms.
    flags = np.asarray(statics.boundary_flags, np.uint16)
    low_end = np.concatenate([(flags & (FLAG_TSS_POS | FLAG_TES_NEG)) != 0, [False]])
    high_end = np.concatenate([(flags & (FLAG_TES_POS | FLAG_TSS_NEG)) != 0, [False]])

    return StructuralClaims(
        n_slots=int(chain.n_slots),
        intergenic=intergenic,
        ss_intron_region=ss_intron_region,
        ss_intron_boundary=ss_intron_boundary,
        solvable_exon=exon_flank_left | exon_flank_right,
        exon_flank_left=exon_flank_left,
        exon_flank_right=exon_flank_right,
        exon_flank_left_complete=exon_flank_left & ~low_end[left],
        exon_flank_right_complete=exon_flank_right & ~high_end[right],
    )
