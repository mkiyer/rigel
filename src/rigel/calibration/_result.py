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

    # ---- Per-(category, strand) counts (shape (N_CATEGORIES, 4)) ----
    # See `rigel.calibration._categorize` for the row/column semantics:
    #   rows    = FragmentCategory  (INTERGENIC, INTRONIC,
    #                                EXON_CONTAINED, EXON_INCOMPATIBLE)
    #   columns = StrandLabel       (NONE, POS, NEG, AMBIG)
    # Counts only unique-mapper UNSPLICED fragments; non-UNSPLICED
    # fragments are excluded from calibration entirely and are tallied
    # separately under `n_spliced`.
    category_counts: np.ndarray  # int64[N_CATEGORIES, 4]
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
    pi_pool: float
    mixture_converged: bool
    mixture_iterations: int

    n_pool_intronic_strand_pos: int
    n_pool_intronic_strand_neg: int
    """Strand-asymmetry diagnostic for the INTRONIC bucket.  Extreme
    asymmetry given high library strand specificity is the signal of
    nascent-RNA pollution of the calibration pool (a structural
    identifiability limit of the 1-D mixture; see
    ``docs/calibration/srd_v2_phase2plus_handoff.md`` §7a).
    """

    n_pool_dropped_out_of_range: int
    """Pool fragments dropped because ``frag_length`` was outside
    ``[0, max_size]``. With SRD v2 Phase 1 sourcing length from
    ``genomic_footprint`` (always >= 0), this counts only fragments
    longer than ``max_size``. Expected to be <0.1% of the pool.
    """

    # ---- Config echo (for reproducibility) ----
    exon_fit_tolerance_bp: int
    fl_prior_ess: float

    # ---- Free-form extras (warnings, debug strings) ----
    extra: dict[str, Any] = field(default_factory=dict)

    def to_summary_dict(self) -> dict[str, Any]:
        """JSON-serializable summary for ``summary.json``.

        ``category_counts`` is flattened row-major:
        ``flat[i * N_STRAND_LABELS + j]`` = count of fragments in
        ``(FragmentCategory(i), StrandLabel(j))``.
        """
        return {
            "gdna_fl_quality": self.gdna_fl_quality,
            "strand_specificity": float(self.strand_specificity),
            "category_counts": [int(c) for c in np.asarray(self.category_counts).ravel()],
            "category_counts_shape": list(np.asarray(self.category_counts).shape),
            "n_multimap_excluded": int(self.n_multimap_excluded),
            "n_spliced": int(self.n_spliced),
            "n_pool": int(self.n_pool),
            "n_pool_intronic_strand_pos": int(self.n_pool_intronic_strand_pos),
            "n_pool_intronic_strand_neg": int(self.n_pool_intronic_strand_neg),
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
