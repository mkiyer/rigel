"""rigel.calibration.fl — the v6 fragment-length distribution surface.

**One owner, one policy, one product.**

* **Owner.**  This module is the sole owner of finalized FL
  distributions on the v6 calibration path.  No other module finalises
  histograms, applies priors, or interprets raw count vectors.
* **Policy.**  Empirical-Bayes Dirichlet shrinkage of ``rna`` and
  ``gdna`` toward ``global_`` with a single shared ``prior_ess``.
  ``global_`` is the unconditional anchor (no prior).
* **Product.**  :class:`FLModels` carries three finalized
  :class:`FragmentLengthModel` instances + per-pool quality flags +
  raw-count totals.  This is the only finalized FL surface downstream
  code (EM scorer, ``summary.json``) sees.

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
    "POOL_QUALITY_GOOD_THRESHOLD",
    "POOL_QUALITY_WEAK_THRESHOLD",
    "build_fl_models",
]


Quality = Literal["good", "weak", "fallback"]

#: Pool with at least this many fragments uses pure-empirical FL.
POOL_QUALITY_GOOD_THRESHOLD: int = 5000

#: Pool with at least this many but below ``GOOD`` is EB-shrunk to ``global_``.
POOL_QUALITY_WEAK_THRESHOLD: int = 200

#: Effective sample size of the global Dirichlet prior in the ``weak`` branch.
POOL_EB_PRIOR_ESS: float = 1000.0


# ---------------------------------------------------------------------------
# Public type
# ---------------------------------------------------------------------------

@dataclass(frozen=True, slots=True)
class FLModels:
    """The three finalized FL distributions of v6 calibration.

    Sole product of the calibration FL pipeline; sole input the EM
    scorer needs.
    """

    global_: FragmentLengthModel       # unconditional anchor (no prior)
    rna:     FragmentLengthModel       # SPLICED, EB-shrunk to global_
    gdna:    FragmentLengthModel       # gDNA pool, EB-shrunk to global_

    rna_quality:  Quality
    gdna_quality: Quality

    n_global: int
    n_rna:    int
    n_gdna:   int

    def to_summary_dict(self) -> dict[str, object]:
        return {
            "n_global":       int(self.n_global),
            "n_rna":          int(self.n_rna),
            "n_gdna":         int(self.n_gdna),
            "rna_quality":    self.rna_quality,
            "gdna_quality":   self.gdna_quality,
            "global_fl_mean": float(self.global_.mean),
            "rna_fl_mean":    float(self.rna.mean),
            "gdna_fl_mean":   float(self.gdna.mean),
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
    n = int(counts.sum())
    if n >= good_threshold:
        return _finalize(counts, max_size, prior_counts=None, prior_ess=None), "good"
    if n >= weak_threshold:
        return (
            _finalize(counts, max_size, prior_counts=global_counts, prior_ess=prior_ess),
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
    rna_counts:    np.ndarray,
    gdna_counts:   np.ndarray,
    max_size:      int,
    prior_ess:      float = POOL_EB_PRIOR_ESS,
    good_threshold: int   = POOL_QUALITY_GOOD_THRESHOLD,
    weak_threshold: int   = POOL_QUALITY_WEAK_THRESHOLD,
) -> FLModels:
    """Build the three finalized FL distributions under one EB policy.

    Parameters
    ----------
    global_counts, rna_counts, gdna_counts
        Raw ``int64`` FL count vectors emitted by the scanner /
        calibration accumulator.  See ``_fl_sources.py`` for the
        canonical extractors.
    max_size
        Maximum FL bin (size of the FL histograms minus one).
    prior_ess
        Effective sample size of the Dirichlet prior used in the
        ``weak`` quality branch.
    good_threshold, weak_threshold
        Pool-size cutoffs for the quality classifier.
    """
    n_global = int(global_counts.sum())
    n_rna    = int(rna_counts.sum())
    n_gdna   = int(gdna_counts.sum())

    global_ = _finalize(global_counts, max_size, prior_counts=None, prior_ess=None)

    rna, rna_q = _classify_and_build(
        rna_counts, global_counts, global_,
        max_size=max_size, prior_ess=prior_ess,
        good_threshold=good_threshold, weak_threshold=weak_threshold,
    )
    gdna, gdna_q = _classify_and_build(
        gdna_counts, global_counts, global_,
        max_size=max_size, prior_ess=prior_ess,
        good_threshold=good_threshold, weak_threshold=weak_threshold,
    )

    return FLModels(
        global_=global_, rna=rna, gdna=gdna,
        rna_quality=rna_q, gdna_quality=gdna_q,
        n_global=n_global, n_rna=n_rna, n_gdna=n_gdna,
    )
