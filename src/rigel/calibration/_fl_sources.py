"""rigel.calibration._fl_sources — raw FL count vector extractors (private).

The **only** Python code that knows where v6 raw FL count histograms
live.  Three one-line functions, one per channel:

* ``extract_global_counts``  — every observed unique-mapper fragment
                                (from the BAM scanner).
* ``extract_rna_counts``     — SPLICED-ANNOT fragments only
                                (from the BAM scanner; sole authoritative
                                RNA source — covers mRNA AND nRNA).
* ``extract_gdna_counts``    — UNSPLICED ∩ {INTRON_ONLY, EXON_INTRON,
                                INTERGENIC_ONLY}
                                (from the calibration accumulator).

The 3-bit region mask is referenced exactly once in this module
(``extract_gdna_counts``) and nowhere else on the v6 path.  After M7+a
(C++ ``fl_hist`` collapse to 3×N), only these three function bodies
change.

See ``docs/calibration/m7_implementation_plan.md`` §4.2.
"""

from __future__ import annotations

import numpy as np

from ..frag_length_model import FragmentLengthModels
from .scan_payload import (
    MASK_EXON,
    MASK_INTERGENIC,
    MASK_INTRON,
    CalibrationScanPayload,
)


__all__ = ["extract_global_counts", "extract_rna_counts", "extract_gdna_counts"]


def extract_global_counts(scan_trained: FragmentLengthModels) -> np.ndarray:
    """Raw global FL histogram from the scanner (every observed fragment)."""
    return np.asarray(scan_trained.global_model.counts, dtype=np.int64)


def extract_rna_counts(scan_trained: FragmentLengthModels) -> np.ndarray:
    """Raw SPLICED-ANNOT FL histogram from the scanner.

    SPLICED fragments are the sole authoritative RNA source: they are
    pure RNA by construction (no genomic-DNA fragment can carry an
    annotated splice junction) and the same FL distribution applies to
    both mature mRNA and nascent RNA.
    """
    from ..splice import SpliceType
    return np.asarray(
        scan_trained.category_models[SpliceType.SPLICED_ANNOT].counts,
        dtype=np.int64,
    )


def extract_gdna_counts(payload: CalibrationScanPayload) -> np.ndarray:
    """Raw gDNA-pool FL histogram from the calibration accumulator.

    Sums mask rows {INTRON_ONLY, EXON_INTRON, INTERGENIC_ONLY}.  All
    three are unspliced ⇒ ``frag_end - frag_start`` IS the fragment
    length.  Mask EXON_ONLY (0b001) is excluded — it predominantly
    contains spliced fragments whose genomic span is **not** the FL.
    """
    h = payload.fl_hist
    return (
        h[MASK_INTRON]
        + h[MASK_EXON | MASK_INTRON]
        + h[MASK_INTERGENIC]
    ).astype(np.int64, copy=False)
