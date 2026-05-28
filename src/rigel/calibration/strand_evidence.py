"""Strand-channel gDNA evidence (PR 07 Phase 2).

This module defines :class:`StrandGdnaEvidence`, the per-region orientation
evidence input to the Fractional Mixture Allocation fusion engine
(``fusion.py``). It replaces the gDNA estimates produced by the legacy
``RegionGdnaEstimate`` tier ladder with a continuous, schema-clean view that
decouples three concerns:

* **Structural applicability** — whether the region overlaps an annotated
  unique-strand feature (`TS_POS` xor `TS_NEG`). This is the only gate that
  is binary. It is derived strictly from annotation and never from observed
  counts.
* **Library-level confidence** (`global_info_scale`) — a continuous scalar in
  ``[0, 1]`` derived from the CI on the library-wide signed strand contrast.
  When the library is near unstranded or under-sampled, this shrinks to 0.
* **Region-level effective support** (`region_info_gain`) — a continuous
  factor in ``[0, 1]`` derived from the per-region effective sample size
  ``eff_support_r = T_sense_r + T_anti_r`` (fractional mass, not raw fragment
  count). Multi-mapper-diluted regions therefore lose strand information
  smoothly.

The product ``g_lib · g_region · applicable`` is the continuous gain factor
used by the fusion engine to interpolate the per-region effective crosstalk
``kappa_eff`` toward the neutral value ``0.5`` when evidence is weak.

See ``docs/newcalib/pr07_strand_density_integration_v2.md`` for the full
mathematical derivation and the rationale for the eff-support vs raw-support
split.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .fractional_evidence import TS_AMBIG, TS_NEG, TS_NONE, TS_POS
from .region_count_ledger import RegionCountLedger
from .strand_deconv import StrandRegionCounts
from .strand_summary import StrandSummary

__all__ = [
    "STRAND_INFO_PSEUDOCOUNT",
    "StrandGdnaEvidence",
    "build_strand_gdna_evidence",
    "compute_global_info_scale",
]


# Smoothing pseudocount for the per-region depth-driven information gain
# ``g_region(eff_support) = eff_support / (eff_support + STRAND_INFO_PSEUDOCOUNT)``.
# A value of 4.0 implies half-weight at ~4 fragment-equivalents of effective
# support. Continuous and monotone; no threshold.
STRAND_INFO_PSEUDOCOUNT: float = 4.0

# Numerical floor for the standard error in the library-level gain factor
# ``g_library = c^2 / (c^2 + (z*se)^2 + eps)``.
_GLOBAL_INFO_SE_EPS: float = 1.0e-12

# Normal-approximation z-score for the library-level CI (95%).
_GLOBAL_INFO_Z: float = 1.959964


@dataclass(frozen=True, slots=True)
class StrandGdnaEvidence:
    """Per-region orientation evidence for the FMA fusion engine.

    See module docstring for the role of each field. All array fields are
    aligned to the region table row order (length ``R``).
    """

    n_sense: np.ndarray              # float64[R] — T_sense_r (fractional mass)
    n_anti: np.ndarray               # float64[R] — T_anti_r  (fractional mass)
    n_total: np.ndarray              # float64[R] — T_r
    support_total: np.ndarray        # uint64[R]  — N_contained + N_left + N_right
    eff_support: np.ndarray          # float64[R] — T_sense + T_anti
    kappa: np.ndarray                # float64[R] — per-region crosstalk (global broadcast in v2)
    applicable: np.ndarray           # bool[R]    — STRUCTURAL only
    structural_absent: np.ndarray    # bool[R]    — == ~applicable
    p_r1_sense: float
    global_info_scale: float         # in [0, 1], library-level
    region_info_gain: np.ndarray     # float64[R] in [0, 1], depth-driven


def compute_global_info_scale(strand_summary: StrandSummary) -> float:
    """Library-level information gain factor in ``[0, 1]``.

    ``g_library = c^2 / (c^2 + (z*se)^2 + eps)`` where ``c`` is the absolute
    signed strand contrast and ``z*se`` is the normal-approximation CI margin
    at 95 %. Returns 0 when the standard error is infinite (no observations).
    """
    c = abs(float(strand_summary.signed_strand_contrast))
    se = float(strand_summary.signed_strand_contrast_se)
    if not np.isfinite(se):
        return 0.0
    margin = _GLOBAL_INFO_Z * se
    c_sq = c * c
    denom = c_sq + margin * margin + _GLOBAL_INFO_SE_EPS
    if denom <= 0.0:
        return 0.0
    return float(c_sq / denom)


def build_strand_gdna_evidence(
    *,
    strand_counts: StrandRegionCounts,
    strand_summary: StrandSummary,
    ledger: RegionCountLedger,
    ts_class: np.ndarray,
    pseudocount: float = STRAND_INFO_PSEUDOCOUNT,
) -> StrandGdnaEvidence:
    """Assemble per-region strand evidence from raw inputs.

    ``applicable`` is computed strictly from ``ts_class`` (annotation): only
    ``TS_POS`` and ``TS_NEG`` regions have a sense/antisense channel; ``TS_NONE``
    and ``TS_AMBIG`` regions are structurally absent and the fusion engine
    will fall back to density-only allocation for them.
    """
    n_sense = np.asarray(strand_counts.k_sense, dtype=np.float64)
    n_anti = np.asarray(strand_counts.k_antisense, dtype=np.float64)
    n_total = np.asarray(strand_counts.n_total, dtype=np.float64)

    ts = np.asarray(ts_class)
    # STRUCTURAL ONLY — derived from annotation, never from counts.
    applicable = (ts == TS_POS) | (ts == TS_NEG)
    structural_absent = (ts == TS_NONE) | (ts == TS_AMBIG)
    # Defensive: any unknown ts code is treated as structurally absent.
    assert np.all(applicable | structural_absent), (
        "ts_class contains values outside {TS_POS, TS_NEG, TS_NONE, TS_AMBIG}"
    )

    support_total = (
        ledger.contained_unspliced_support.astype(np.uint64)
        + ledger.boundary_left_unspliced_support.astype(np.uint64)
        + ledger.boundary_right_unspliced_support.astype(np.uint64)
    )

    eff_support = n_sense + n_anti
    pc = float(pseudocount)
    region_info_gain = eff_support / (eff_support + pc)

    p_r1_sense = float(strand_summary.p_r1_sense)
    kappa_lib = float(min(p_r1_sense, 1.0 - p_r1_sense))
    kappa = np.full(n_total.shape, kappa_lib, dtype=np.float64)

    g_lib = compute_global_info_scale(strand_summary)

    return StrandGdnaEvidence(
        n_sense=n_sense,
        n_anti=n_anti,
        n_total=n_total,
        support_total=support_total,
        eff_support=eff_support,
        kappa=kappa,
        applicable=applicable,
        structural_absent=structural_absent,
        p_r1_sense=p_r1_sense,
        global_info_scale=g_lib,
        region_info_gain=region_info_gain,
    )
