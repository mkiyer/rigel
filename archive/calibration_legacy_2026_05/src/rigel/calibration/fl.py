"""rigel.calibration.fl — the v6 fragment-length distribution surface.

**One owner, one policy, one product.**

* **Owner.**  This module is the sole owner of finalized FL
  distributions on the v6 calibration path.  No other module finalises
  histograms, applies priors, or interprets raw count vectors.
* **Policy.**  Empirical-Bayes Dirichlet shrinkage of ``rna`` and
    ``gdna`` toward ``global_`` with a shared maximum ``prior_ess`` and
    an evidence-adaptive cap. ``global_`` is the unconditional anchor
    (no prior).
* **Product.**  :class:`FLModels` carries raw diagnostic/calibration
    distributions plus contrast-regularized scoring distributions for EM.
    Downstream code must choose the surface that matches its statistical role.

Producers (BAM scanner + calibration C++ accumulator) emit raw
``int64`` count vectors and stop there.  ``build_fl_models`` consumes
three vectors and returns one ``FLModels``.

See ``docs/calibration/m7_implementation_plan.md``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from ..frag_length_model import FragmentLengthModel


__all__ = [
    "FLModels",
    "Quality",
    "POOL_EB_PRIOR_ESS",
    "POOL_SCORING_PRIOR_ESS",
    "POOL_QUALITY_GOOD_THRESHOLD",
    "POOL_QUALITY_WEAK_THRESHOLD",
    "build_fl_models",
]


Quality = Literal["good", "weak", "fallback"]

#: Pool with at least this many fragments uses pure-empirical FL.
POOL_QUALITY_GOOD_THRESHOLD: int = 5000

#: Pool with at least this much positive evidence but below ``GOOD`` is
#: EB-shrunk to ``global_``. The default keeps non-empty low-contamination
#: pools informative while still allowing callers to request stricter fallback.
POOL_QUALITY_WEAK_THRESHOLD: int = 1

#: Effective sample size of the global Dirichlet prior in the ``weak`` branch.
POOL_EB_PRIOR_ESS: float = 1000.0

#: Effective sample size used to decide how much RNA-vs-gDNA FL contrast is
#: reliable for scoring. This regularizes the contrast, not the raw diagnostic
#: pools: if either side lacks support, FL scoring falls back to the shared
#: global surface and contributes no class-specific likelihood ratio.
POOL_SCORING_PRIOR_ESS: float = 200.0

_SCORING_PMF_FLOOR: float = 1e-300


# ---------------------------------------------------------------------------
# Public type
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class FLModels:
    """Finalized FL distributions and score-side contrast surfaces."""

    global_: FragmentLengthModel  # unconditional anchor (no prior)
    rna: FragmentLengthModel  # SPLICED, EB-shrunk to global_
    gdna: FragmentLengthModel  # gDNA pool, EB-shrunk to global_
    rna_scoring: FragmentLengthModel  # RNA FL surface used by EM scoring
    gdna_scoring: FragmentLengthModel  # gDNA FL surface used by EM scoring

    rna_quality: Quality
    gdna_quality: Quality
    rna_fl_reliability: float
    gdna_fl_reliability: float
    fl_contrast_weight: float
    fl_scoring_prior_ess: float

    # Pool totals are float64 — the fractional accumulator emits
    # weighted mass, not integer counts. Existing call sites consume
    # these as numeric quantities only.
    n_global: float
    n_rna: float
    n_gdna: float

    def to_summary_dict(self) -> dict[str, object]:
        return {
            "n_global": float(self.n_global),
            "n_rna": float(self.n_rna),
            "n_gdna": float(self.n_gdna),
            "rna_quality": self.rna_quality,
            "gdna_quality": self.gdna_quality,
            "global_fl_mean": float(self.global_.mean),
            "rna_fl_mean": float(self.rna.mean),
            "gdna_fl_mean": float(self.gdna.mean),
            "rna_scoring_fl_mean": float(self.rna_scoring.mean),
            "gdna_scoring_fl_mean": float(self.gdna_scoring.mean),
            "rna_fl_reliability": float(self.rna_fl_reliability),
            "gdna_fl_reliability": float(self.gdna_fl_reliability),
            "fl_contrast_weight": float(self.fl_contrast_weight),
            "fl_scoring_prior_ess": float(self.fl_scoring_prior_ess),
        }


# ---------------------------------------------------------------------------
# The single EB primitive
# ---------------------------------------------------------------------------


def _finalize(
    counts: np.ndarray,
    max_size: int,
    *,
    prior_counts: np.ndarray | None,
    prior_ess: float | None,
) -> FragmentLengthModel:
    """Build and finalise one ``FragmentLengthModel`` from a raw count vector.

    Used three times in :func:`build_fl_models` (once per channel).
    """
    fl = FragmentLengthModel(max_size=max_size)
    n = min(fl.counts.size, counts.size)
    fl.counts[:n] = counts[:n].astype(np.float64, copy=False)
    fl._total_weight = float(fl.counts.sum())
    if prior_counts is None:
        fl.finalize()
    else:
        fl.finalize(
            prior_counts=prior_counts.astype(np.float64, copy=False),
            prior_ess=prior_ess,
        )
    return fl


def _adaptive_prior_ess(prior_ess: float, pool_weight: float) -> float:
    """Cap weak-pool shrinkage so global FL cannot swamp scarce evidence."""
    if prior_ess <= 0.0 or pool_weight <= 0.0:
        return 0.0
    return float(min(prior_ess, np.log1p(pool_weight)))


def _aligned_counts(counts: np.ndarray, max_size: int) -> np.ndarray:
    """Align a raw FL count vector to ``max_size + 1`` bins."""
    out = np.zeros(max_size + 1, dtype=np.float64)
    n = min(out.size, counts.size)
    out[:n] = counts[:n].astype(np.float64, copy=False)
    if counts.size > out.size:
        out[max_size] += float(np.sum(counts[out.size :], dtype=np.float64))
    return out


def _empirical_pmf_or_global(counts: np.ndarray, global_pmf: np.ndarray) -> np.ndarray:
    total = float(np.sum(counts, dtype=np.float64))
    if total <= 0.0:
        return global_pmf.copy()
    return counts / total


def _pool_reliability(quality: Quality, pool_weight: float, scoring_prior_ess: float) -> float:
    if quality == "fallback" or pool_weight <= 0.0:
        return 0.0
    if scoring_prior_ess <= 0.0:
        return 1.0
    return float(pool_weight / (pool_weight + scoring_prior_ess))


def _normalize_scoring_pmf(pmf: np.ndarray) -> np.ndarray:
    out = np.maximum(np.asarray(pmf, dtype=np.float64), _SCORING_PMF_FLOOR)
    total = float(np.sum(out, dtype=np.float64))
    if total <= 0.0 or not np.isfinite(total):
        raise ValueError("FL scoring PMF has invalid total mass.")
    return out / total


def _build_scoring_models(
    *,
    rna_counts: np.ndarray,
    gdna_counts: np.ndarray,
    global_model: FragmentLengthModel,
    rna_quality: Quality,
    gdna_quality: Quality,
    n_rna: float,
    n_gdna: float,
    max_size: int,
    scoring_prior_ess: float,
) -> tuple[FragmentLengthModel, FragmentLengthModel, float, float, float]:
    """Build contrast-regularized RNA/gDNA FL surfaces for EM scoring."""
    global_pmf = _normalize_scoring_pmf(global_model.pmf)
    aligned_rna = _aligned_counts(rna_counts, max_size)
    aligned_gdna = _aligned_counts(gdna_counts, max_size)

    rna_reliability = _pool_reliability(rna_quality, n_rna, scoring_prior_ess)
    gdna_reliability = _pool_reliability(gdna_quality, n_gdna, scoring_prior_ess)
    contrast_weight = min(rna_reliability, gdna_reliability)

    rna_empirical = _empirical_pmf_or_global(aligned_rna, global_pmf)
    gdna_empirical = _empirical_pmf_or_global(aligned_gdna, global_pmf)
    rna_score_pmf = _normalize_scoring_pmf(
        contrast_weight * rna_empirical + (1.0 - contrast_weight) * global_pmf
    )
    gdna_score_pmf = _normalize_scoring_pmf(
        contrast_weight * gdna_empirical + (1.0 - contrast_weight) * global_pmf
    )

    return (
        FragmentLengthModel.from_pmf(rna_score_pmf, max_size=max_size),
        FragmentLengthModel.from_pmf(gdna_score_pmf, max_size=max_size),
        rna_reliability,
        gdna_reliability,
        contrast_weight,
    )


def _classify_and_build(
    counts: np.ndarray,
    global_counts: np.ndarray,
    global_model: FragmentLengthModel,
    *,
    max_size: int,
    prior_ess: float,
    good_threshold: int,
    weak_threshold: int,
) -> tuple[FragmentLengthModel, Quality]:
    """One classifier; called twice (once for RNA, once for gDNA)."""
    # Pool size is fractional under the cutover; compare on weighted mass.
    n = float(counts.sum())
    if n >= good_threshold:
        return _finalize(counts, max_size, prior_counts=None, prior_ess=None), "good"
    if n > 0.0 and n >= weak_threshold:
        return (
            _finalize(
                counts,
                max_size,
                prior_counts=global_counts,
                prior_ess=_adaptive_prior_ess(prior_ess, n),
            ),
            "weak",
        )
    # Identity-share with the global anchor — no recomputation.
    return global_model, "fallback"


# ---------------------------------------------------------------------------
# The single owner
# ---------------------------------------------------------------------------


def build_fl_models(
    *,
    global_counts: np.ndarray,
    rna_counts: np.ndarray,
    gdna_counts: np.ndarray,
    max_size: int,
    prior_ess: float = POOL_EB_PRIOR_ESS,
    scoring_prior_ess: float = POOL_SCORING_PRIOR_ESS,
    good_threshold: int = POOL_QUALITY_GOOD_THRESHOLD,
    weak_threshold: int = POOL_QUALITY_WEAK_THRESHOLD,
) -> FLModels:
    """Build the three finalized FL distributions under one EB policy.

    Parameters
    ----------
    global_counts, rna_counts, gdna_counts
        Raw ``float64`` FL weight vectors emitted by the scanner /
        calibration accumulator. See ``_fl_sources.py`` for the
        canonical extractors.
    max_size
        Maximum FL bin (size of the FL histograms minus one).
    prior_ess
        Effective sample size of the Dirichlet prior used in the
        ``weak`` quality branch.
    scoring_prior_ess
        Effective sample size for the joint RNA-vs-gDNA FL contrast used by
        EM scoring. This does not alter raw diagnostic pool models.
    good_threshold, weak_threshold
        Pool-size cutoffs for the quality classifier.
    """
    n_global = float(global_counts.sum())
    n_rna = float(rna_counts.sum())
    n_gdna = float(gdna_counts.sum())

    global_ = _finalize(global_counts, max_size, prior_counts=None, prior_ess=None)

    rna, rna_q = _classify_and_build(
        rna_counts,
        global_counts,
        global_,
        max_size=max_size,
        prior_ess=prior_ess,
        good_threshold=good_threshold,
        weak_threshold=weak_threshold,
    )
    gdna, gdna_q = _classify_and_build(
        gdna_counts,
        global_counts,
        global_,
        max_size=max_size,
        prior_ess=prior_ess,
        good_threshold=good_threshold,
        weak_threshold=weak_threshold,
    )

    rna_scoring, gdna_scoring, rna_reliability, gdna_reliability, contrast_weight = (
        _build_scoring_models(
            rna_counts=rna_counts,
            gdna_counts=gdna_counts,
            global_model=global_,
            rna_quality=rna_q,
            gdna_quality=gdna_q,
            n_rna=n_rna,
            n_gdna=n_gdna,
            max_size=max_size,
            scoring_prior_ess=scoring_prior_ess,
        )
    )

    return FLModels(
        global_=global_,
        rna=rna,
        gdna=gdna,
        rna_scoring=rna_scoring,
        gdna_scoring=gdna_scoring,
        rna_quality=rna_q,
        gdna_quality=gdna_q,
        rna_fl_reliability=rna_reliability,
        gdna_fl_reliability=gdna_reliability,
        fl_contrast_weight=contrast_weight,
        fl_scoring_prior_ess=float(scoring_prior_ess),
        n_global=n_global,
        n_rna=n_rna,
        n_gdna=n_gdna,
    )
