"""rigel.calibration._orchestrator — top-level ``calibrate(...)``.

Composes the M3–M7 building blocks into a single entry point that
replaces the SRD-v1 :func:`rigel.calibration.calibrate_gdna`.

The function is deliberately thin: it builds the FL models once,
uses the resulting gDNA-FL mean to compute the three global gDNA
densities, and assembles a :class:`CalibrationResult` with an empty
:class:`PriorTable`.  Per-``MultiLocus`` priors are filled in later
by the pipeline via :func:`assemble_priors` →
:meth:`CalibrationResult.with_priors`.

There is no recomputation: the FL models built here are passed
through to :func:`build_calibration_result` so it does not rebuild
them.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from ._fl_sources import (
    extract_gdna_counts,
    extract_global_counts,
    extract_rna_counts,
)
from ._result_v6 import CalibrationResult, build_calibration_result
from .density_global import compute_global_densities
from .fl import POOL_EB_PRIOR_ESS, build_fl_models
from .scan_payload import CalibrationScanPayload

if TYPE_CHECKING:
    from ..frag_length_model import FragmentLengthModels
    from ..index import TranscriptIndex


__all__ = ["calibrate"]


def calibrate(
    *,
    index:        "TranscriptIndex",
    payload:      CalibrationScanPayload,
    scan_trained: "FragmentLengthModels",
    fl_prior_ess: float = POOL_EB_PRIOR_ESS,
) -> CalibrationResult:
    """Run the v6 calibration pipeline end-to-end (sans per-locus priors).

    Parameters
    ----------
    index
        The loaded :class:`TranscriptIndex`.  Must carry a non-``None``
        ``region_df`` (built by M1 / loaded by M2).  Stale indexes raise
        :class:`RuntimeError` with a rebuild instruction.
    payload
        The :class:`CalibrationScanPayload` produced by the C++ scanner
        (M3).  Carries the 8-state mask histogram, per-region counts,
        boundary-flux counters, and exclusion totals.
    scan_trained
        The scanner-trained :class:`FragmentLengthModels`.  Used as the
        source of raw global and SPLICED-annotated count vectors that
        feed :func:`build_fl_models`.
    fl_prior_ess
        Empirical-Bayes evidence strength for the FL Dirichlet
        shrinkage.  Defaults to :data:`POOL_EB_PRIOR_ESS`.

    Returns
    -------
    CalibrationResult
        With ``prior_table = PriorTable.empty()``; the caller is
        expected to build the per-``MultiLocus`` prior table via
        :func:`rigel.calibration.locus_prior.assemble_priors` and swap
        it in via :meth:`CalibrationResult.with_priors`.
    """
    if index.region_df is None:
        raise RuntimeError(
            "Index has no region table. Rebuild the index "
            "(rigel index --fasta ... --gtf ...). Older indexes "
            "are not supported."
        )

    fl_models = build_fl_models(
        global_counts = extract_global_counts(scan_trained),
        rna_counts    = extract_rna_counts(scan_trained),
        gdna_counts   = extract_gdna_counts(payload),
        max_size      = scan_trained.max_size,
        prior_ess     = fl_prior_ess,
    )

    global_densities = compute_global_densities(
        index.region_df,
        payload,
        gdna_fl_mean=float(fl_models.gdna.mean),
    )

    return build_calibration_result(
        payload          = payload,
        scan_trained     = scan_trained,
        global_densities = global_densities,
        fl_models        = fl_models,
        fl_prior_ess     = fl_prior_ess,
    )
