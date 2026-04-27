"""SRD v1 Pass 3 — Empirical-Bayes FL model construction.

Builds three :class:`FragmentLengthModel` instances:

* ``global_fl`` — finalized from the unconditional global FL histogram.
  Acts as the Dirichlet prior anchor for the other two models.
* ``rna_fl`` — finalized from the spliced (SPLICED-category) histogram,
  Dirichlet-smoothed toward ``global_fl`` with effective sample size
  ``prior_ess``. Collapses to ``global_fl`` when SPLICED counts are
  negligible (QC-fail libraries).
* ``gdna_fl`` — finalized from the Pass 2 mixture-recovered histogram,
  Dirichlet-smoothed toward ``global_fl`` with the same ESS. Falls back
  to ``global_fl`` when the mixture skipped or failed.

The smoothing primitive is ``FragmentLengthModel.finalize(prior_counts=...,
prior_ess=...)`` which already implements the Dirichlet shrinkage; no
new EB code is needed here.
"""

from __future__ import annotations

import numpy as np

from ..frag_length_model import FragmentLengthModel


def _materialize(counts: np.ndarray, max_size: int) -> FragmentLengthModel:
    """Build an unfinalized :class:`FragmentLengthModel` from a count vector."""
    fl = FragmentLengthModel(max_size=max_size)
    n = min(fl.counts.size, counts.size)
    fl.counts[:n] = counts[:n]
    fl._total_weight = float(fl.counts.sum())
    return fl


def build_global_fl(global_counts: np.ndarray, max_size: int) -> FragmentLengthModel:
    """Finalize the global FL model (no prior — it *is* the prior anchor)."""
    fl = _materialize(global_counts, max_size)
    fl.finalize()
    return fl


def build_rna_fl(
    spliced_counts: np.ndarray,
    global_counts: np.ndarray,
    max_size: int,
    *,
    prior_ess: float,
) -> FragmentLengthModel:
    """Finalize ``RNA_FL`` from spliced fragments, EB-shrunk to global FL."""
    fl = _materialize(spliced_counts, max_size)
    fl.finalize(prior_counts=global_counts, prior_ess=prior_ess)
    return fl


def build_gdna_fl(
    gdna_counts: np.ndarray | None,
    global_counts: np.ndarray,
    max_size: int,
    *,
    prior_ess: float,
) -> FragmentLengthModel:
    """Finalize ``gDNA_FL`` from mixture output, EB-shrunk to global FL.

    When ``gdna_counts`` is ``None`` or sums to zero, returns a model
    finalized purely from the prior — i.e. ``gDNA_FL = global_FL``.
    """
    counts = (
        np.asarray(gdna_counts, dtype=np.float64)
        if gdna_counts is not None
        else np.zeros(max_size + 1, dtype=np.float64)
    )
    fl = _materialize(counts, max_size)
    fl.finalize(prior_counts=global_counts, prior_ess=prior_ess)
    return fl
