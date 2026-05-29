"""rigel.calibration._fl_sources \u2014 raw FL count vector extractors (private).

The **only** Python code that knows where v6 raw FL count vectors live.
Three one-line functions, one per channel:

* ``extract_global_counts`` \u2014 every observed unique-mapper fragment
                                (from the BAM scanner).
* ``extract_rna_counts``    \u2014 SPLICED-ANNOT fragments only
                                (from the BAM scanner; sole authoritative
                                RNA source \u2014 covers mRNA AND nRNA).
* ``extract_gdna_counts``   \u2014 INTERGENIC + INTRONIC FL pools
                                (from the calibration accumulator).

All three return ``float64`` vectors. After the fractional cutover,
``extract_gdna_counts`` aggregates the four non-EXON FL pools.
"""

from __future__ import annotations

import numpy as np

from ..frag_length_model import FragmentLengthModels
from .fractional_evidence import gdna_fl_mass
from .scan_payload import CalibrationScanPayload


__all__ = ["extract_global_counts", "extract_rna_counts", "extract_gdna_counts"]


def extract_global_counts(scan_trained: FragmentLengthModels) -> np.ndarray:
    """Raw global FL histogram from the scanner (every observed fragment)."""
    return np.asarray(scan_trained.global_model.counts, dtype=np.float64)


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
        dtype=np.float64,
    )


def extract_gdna_counts(payload: CalibrationScanPayload) -> np.ndarray:
    """Aggregate gDNA FL mass: INTERGENIC + INTRONIC pools (both compartments)."""
    return gdna_fl_mass(payload)
