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

    # ---- Per-category counts (length 7; matches FragmentCategory enum) ----
    category_counts: np.ndarray  # int64[7]
    n_multimap_excluded: int

    # ---- Pool diagnostics ----
    n_pool: int
    pi_pool: float
    mixture_converged: bool
    mixture_iterations: int
    n_intergenic_unique: int
    """Unique-mapper unspliced fragments with no transcript candidates.

    These are dropped before reaching the buffer (they never appear in
    ``category_counts[INTERGENIC]``) but their FL distribution is folded
    into Pool B from ``frag_length_models.intergenic``. See
    ``calibration._simple.calibrate_gdna``.
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
        """JSON-serializable summary for ``summary.json``."""
        return {
            "gdna_fl_quality": self.gdna_fl_quality,
            "strand_specificity": float(self.strand_specificity),
            "category_counts": [int(c) for c in self.category_counts],
            "n_multimap_excluded": int(self.n_multimap_excluded),
            "n_pool": int(self.n_pool),
            "n_intergenic_unique": int(self.n_intergenic_unique),
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
