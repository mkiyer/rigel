"""SRD v1 :class:`CalibrationResult` — the public calibration output.

Carries the three Empirical-Bayes-shrunk fragment-length models
(``rna_fl_model``, ``gdna_fl_model``, ``global_fl_model``), the
library-wide gDNA pool fraction ``pi_pool``, and per-category
diagnostic counts.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Literal

import numpy as np

from ..frag_length_model import FragmentLengthModel


GdnaFlQuality = Literal["good", "weak", "fallback"]


@dataclass(frozen=True)
class CalibrationResult:
    """SRD v1 calibration output."""

    # ---- Models (the only product the rest of the pipeline consumes) ----
    gdna_fl_model: FragmentLengthModel
    rna_fl_model: FragmentLengthModel
    global_fl_model: FragmentLengthModel
    gdna_fl_quality: GdnaFlQuality

    # ---- Library-level signal passed through ----
    strand_specificity: float

    # ---- Per-(category, strand) counts (shape (N_CATEGORIES, 3)) ----
    # See `rigel.calibration._categorize` for the row/column semantics:
    #   rows    = FragmentCategory  (INTERGENIC, INTRONIC,
    #                                EXON_CONTAINED, EXON_INCOMPATIBLE)
    #   columns = FragmentStrand    (SENSE, ANTISENSE, AMBIG)
    # SRD v3 Phase 1: SENSE/ANTISENSE are in the **transcript frame**
    # for transcript-overlapping rows.  For the INTERGENIC row the
    # labels are reused as a pure naming convention based on the read's
    # genomic strand (SENSE = read on "+", ANTISENSE = read on "-");
    # downstream consumers should treat the INTERGENIC row only as a
    # 50/50 sanity check.  AMBIG marks transcript-overlapping fragments
    # where transcripts on BOTH strands overlap.
    # Counts only unique-mapper UNSPLICED fragments with a known read
    # strand; non-UNSPLICED and AMBIG-read-strand fragments are
    # excluded from calibration entirely.
    category_counts: np.ndarray  # int64[N_CATEGORIES, 3]
    n_multimap_excluded: int
    n_spliced: int
    """Unique-mapper fragments with ``splice_type != UNSPLICED`` (any of
    SPLICED_ANNOT / SPLICED_UNANNOT / SPLICED_IMPLICIT / SPLICE_ARTIFACT).
    Excluded from calibration; surfaced for the QC warning that fires
    when fewer than 100 spliced fragments are seen (``RNA_FL`` then
    collapses to ``global_FL``).
    """

    # ---- Pool diagnostics ----
    n_pool: int
    """Number of unique-mapper UNSPLICED fragments in the calibration
    pool (INTERGENIC ∪ INTRONIC ∪ EXON_INCOMPATIBLE) **after** dropping
    out-of-range lengths. This is the denominator that actually fed the
    1-D mixture EM. Pre-drop count is ``n_pool + n_pool_dropped_out_of_range``.
    """
    pi_pool: float
    mixture_converged: bool
    mixture_iterations: int

    n_pool_intronic_strand_sense: int
    n_pool_intronic_strand_antisense: int
    n_pool_intronic_strand_ambig: int
    """Transcript-frame strand-asymmetry diagnostic for the INTRONIC
    bucket (SRD v3 Phase 1).  In a stranded library, gDNA produces a
    50/50 sense/antisense split (sonicated dsDNA); nRNA produces a
    SS/(1-SS) split.  An INTRONIC sense fraction much closer to the
    library SS than to 0.5 is the direct signal of nascent-RNA
    contamination of the calibration pool.  Phase 2 of v3 will consume
    the strand axis directly via a joint (FL × strand) mixture; until
    then, see ``docs/calibration/srd_v2_phase2plus_handoff.md`` §7a
    for the underlying limitation and ``docs/calibration/srd_v3_early_plan.md``
    for the planned resolution.
    """

    n_pool_dropped_out_of_range: int
    """Pool fragments dropped because ``frag_length`` was outside
    ``[0, max_size]``. With SRD v2 Phase 1 sourcing length from
    ``genomic_footprint`` (always >= 0), this counts only fragments
    longer than ``max_size``. These are biologically anomalous (true
    chimeras, unannotated splicing, read-through into intergenic,
    annotation errors) and intentionally excluded from FL training; the
    count is reported here purely as telemetry. Typically <0.5% of pool;
    investigate if it spikes.
    """

    # ---- Config echo (for reproducibility) ----
    exon_fit_tolerance_bp: int
    fl_prior_ess: float

    # ---- Free-form extras (warnings, debug strings) ----
    extra: dict[str, Any] = field(default_factory=dict)

    def to_summary_dict(self) -> dict[str, Any]:
        """JSON-serializable summary for ``summary.json``.

        ``category_counts`` is flattened row-major:
        ``flat[i * N_FRAGMENT_STRANDS + j]`` = count of fragments in
        ``(FragmentCategory(i), FragmentStrand(j))``.
        """
        return {
            "gdna_fl_quality": self.gdna_fl_quality,
            "strand_specificity": float(self.strand_specificity),
            "category_counts": [int(c) for c in np.asarray(self.category_counts).ravel()],
            "category_counts_shape": list(np.asarray(self.category_counts).shape),
            "n_multimap_excluded": int(self.n_multimap_excluded),
            "n_spliced": int(self.n_spliced),
            "n_pool": int(self.n_pool),
            "n_pool_intronic_strand_sense": int(self.n_pool_intronic_strand_sense),
            "n_pool_intronic_strand_antisense": int(self.n_pool_intronic_strand_antisense),
            "n_pool_intronic_strand_ambig": int(self.n_pool_intronic_strand_ambig),
            "n_pool_dropped_out_of_range": int(self.n_pool_dropped_out_of_range),
            "pi_pool": float(self.pi_pool),
            "mixture_converged": bool(self.mixture_converged),
            "mixture_iterations": int(self.mixture_iterations),
            "gdna_fl_mean": (
                round(self.gdna_fl_model.mean, 2)
                if self.gdna_fl_model is not None
                else None
            ),
            "exon_fit_tolerance_bp": int(self.exon_fit_tolerance_bp),
            "fl_prior_ess": float(self.fl_prior_ess),
            "extra": dict(self.extra),
        }
