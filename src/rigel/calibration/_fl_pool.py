"""rigel.calibration._fl_pool — Pool fragment-length models for v6.

Per ``docs/calibration/calibration_v6_plan.md`` §2.9 and
``docs/calibration/m7_implementation_plan.md`` §3.

Three FL pools, one new + two passthrough:

* ``gdna``        — derived from ``payload.fl_hist[2,3,4]`` (mask 2 INTRON_ONLY
                    ∪ 3 EXON_INTRON ∪ 4 INTERGENIC_ONLY — all unspliced ⇒
                    genomic span == FL).  Quality-classified into
                    ``good`` / ``weak`` / ``fallback``.
* ``rna_spliced`` — passthrough of ``scan_trained.rna_model`` (annotated-spliced).
* ``global_``     — passthrough of ``scan_trained.global_model``.

**Spliced-genomic-span trap:** ``fl_hist[1]`` (mask 1 = EXON_ONLY) is
predominantly *spliced* mRNA whose genomic span is the distance between
the outer read ends — NOT the fragment length.  Histogramming this column
as a FL distribution overestimates FL by 100s of bp.  M7 explicitly
rejects mask 1 as a gDNA FL source; the scanner-trained ``rna_model``
(SPLICED_ANNOT-bin) is the only valid mRNA FL.

There is no scan-trained nRNA FL model: the scanner cannot disambiguate
unspliced RNA from gDNA at scan time.  Downstream code that needs an
"unspliced RNA" FL distribution should use ``global_`` as a stand-in.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np

from ..frag_length_model import FragmentLengthModel, FragmentLengthModels
from ._fl_empirical_bayes import build_gdna_fl, build_global_fl
from .scan_payload import CalibrationScanPayload


__all__ = [
    "PoolFLModels",
    "PoolQuality",
    "POOL_QUALITY_GOOD_THRESHOLD",
    "POOL_QUALITY_WEAK_THRESHOLD",
    "POOL_EB_PRIOR_ESS",
    "compute_pool_fl_models",
]


PoolQuality = Literal["good", "weak", "fallback"]


#: Pool with at least this many fragments uses pure-empirical FL.
POOL_QUALITY_GOOD_THRESHOLD = 5000

#: Pool with at least this many but below ``GOOD`` is EB-shrunk to global.
POOL_QUALITY_WEAK_THRESHOLD = 200

#: Effective sample size of the global Dirichlet prior in the ``weak`` branch.
POOL_EB_PRIOR_ESS = 1000.0


@dataclass(frozen=True, slots=True)
class PoolFLModels:
    """The three pool FL models + provenance / quality diagnostics.

    Field discipline:

    * ``gdna``        — newly derived from the scan payload (one of three
                         quality branches, see :func:`compute_pool_fl_models`).
    * ``rna_spliced`` — object identity with ``scan_trained.rna_model``.
    * ``global_``     — object identity with ``scan_trained.global_model``.
    """

    gdna: FragmentLengthModel
    rna_spliced: FragmentLengthModel
    global_: FragmentLengthModel

    gdna_quality: PoolQuality
    gdna_n_fragments: int
    gdna_fl_mean: float
    gdna_eb_ess: float
    gdna_used_global_fallback: bool

    #: Mask-N counts (N ∈ {0, 1, 5, 6, 7}) that did NOT contribute to gDNA pool.
    n_pool_annotation_gap: dict[str, int]

    def to_summary_dict(self) -> dict[str, object]:
        return {
            "gdna_quality":              self.gdna_quality,
            "gdna_n_fragments":          int(self.gdna_n_fragments),
            "gdna_fl_mean":              float(self.gdna_fl_mean),
            "gdna_eb_ess":               float(self.gdna_eb_ess),
            "gdna_used_global_fallback": bool(self.gdna_used_global_fallback),
            "rna_spliced_fl_mean":       float(self.rna_spliced.mean),
            "global_fl_mean":            float(self.global_.mean),
            "n_pool_annotation_gap":     {k: int(v) for k, v in self.n_pool_annotation_gap.items()},
        }


# ---------------------------------------------------------------------------
# Builder
# ---------------------------------------------------------------------------

#: gDNA-bearing masks (all unspliced ⇒ genomic span == FL).
_GDNA_MASKS: tuple[int, ...] = (0b010, 0b011, 0b100)  # INTRON_ONLY, EXON_INTRON, INTERGENIC_ONLY

#: Masks NOT contributing to gDNA pool (annotation gap diagnostic).
_NON_GDNA_MASKS: tuple[int, ...] = (0b000, 0b001, 0b101, 0b110, 0b111)


def compute_pool_fl_models(
    payload: CalibrationScanPayload,
    scan_trained: FragmentLengthModels,
) -> PoolFLModels:
    """Build the three v6 pool FL models from the scan payload + trained models.

    Pure function: no I/O, no mutation of inputs.  The ``rna_spliced``
    and ``global_`` pools are object-identity passthroughs of the
    scanner-trained models.
    """
    max_size = int(scan_trained.max_size)

    # 1. Sum gDNA-bearing fl_hist columns over the unspliced masks.
    gdna_hist = np.zeros(payload.fl_hist.shape[1], dtype=np.int64)
    for m in _GDNA_MASKS:
        gdna_hist += payload.fl_hist[m]
    n_gdna = int(gdna_hist.sum())

    global_counts = np.asarray(scan_trained.global_model.counts, dtype=np.float64)

    # 2. Quality classifier → branch.
    if n_gdna >= POOL_QUALITY_GOOD_THRESHOLD:
        # Pure empirical (Laplace-smoothed; no global prior shrinkage).
        gdna_model = FragmentLengthModel.from_counts(
            gdna_hist.astype(np.float64), max_size=max_size
        )
        quality: PoolQuality = "good"
        eb_ess = 0.0
        fallback = False

    elif n_gdna >= POOL_QUALITY_WEAK_THRESHOLD:
        gdna_model = build_gdna_fl(
            gdna_hist.astype(np.float64),
            global_counts=global_counts,
            max_size=max_size,
            prior_ess=POOL_EB_PRIOR_ESS,
        )
        quality = "weak"
        eb_ess = POOL_EB_PRIOR_ESS
        fallback = False

    else:
        # < POOL_QUALITY_WEAK_THRESHOLD ⇒ identity copy of global.
        gdna_model = build_global_fl(global_counts, max_size=max_size)
        quality = "fallback"
        eb_ess = 0.0
        fallback = True

    # 3. Annotation-gap diagnostic.
    n_pool_annotation_gap = {
        f"mask_{m}": int(payload.fl_hist[m].sum()) for m in _NON_GDNA_MASKS
    }

    return PoolFLModels(
        gdna=gdna_model,
        rna_spliced=scan_trained.rna_model,
        global_=scan_trained.global_model,
        gdna_quality=quality,
        gdna_n_fragments=n_gdna,
        gdna_fl_mean=float(gdna_model.mean),
        gdna_eb_ess=eb_ess,
        gdna_used_global_fallback=fallback,
        n_pool_annotation_gap=n_pool_annotation_gap,
    )
