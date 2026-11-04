"""Region signature encoding: the 4-bit exon/intron × strand annotation (layer 0).

Every genome position carries a 4-bit signature ``{intron_pos, intron_neg, exon_pos, exon_neg}`` recording which
transcript features cover it; a calibration region (:mod:`rigel.calibration.region_arrays`) is a maximal interval of
constant signature. This module holds the bit constants, the ``RegionType`` / ``RegionStrand`` enums and the coarse
classes derived from a signature: the type (exon wins over intron), the transcript-strand class and the per-strand RNA
masks. The bit values are a wire encoding, and the module imports nothing else in the package.
"""

from __future__ import annotations

from enum import IntEnum, IntFlag

import numpy as np

# ---------------------------------------------------------------------------
# Signature bits (canonical 4-bit encoding)
# ---------------------------------------------------------------------------

BIT_INTRON_POS = 0x8
BIT_INTRON_NEG = 0x4
BIT_EXON_POS = 0x2
BIT_EXON_NEG = 0x1
#: Number of distinct 4-bit signatures (0x0 .. 0xF).
N_SIGNATURES = 16


class RegionType(IntEnum):
    """Coarse region type derived from a signature (exon wins over intron)."""

    INTERGENIC = 0
    INTRON = 1
    EXON = 2


class RegionStrand(IntFlag):
    """Coarse strand flags for a region.

    ``AMBIG`` is the bitwise OR ``POS | NEG``: both transcript strands are
    represented over the region.
    """

    NONE = 0
    POS = 1
    NEG = 2
    AMBIG = 3  # POS | NEG


# The int8 transcript-strand class :class:`~rigel.calibration.region_arrays.RegionArrays` carries: aliases of
# :class:`RegionStrand`, which mirrors :class:`rigel.types.Strand`, so a value in the Strand convention (such as a
# splice junction's strand) compares against TS_POS / TS_NEG with no conversion.
TS_NONE: int = int(RegionStrand.NONE)  # 0
TS_POS: int = int(RegionStrand.POS)  # 1
TS_NEG: int = int(RegionStrand.NEG)  # 2
TS_AMBIG: int = int(RegionStrand.AMBIG)  # 3


# ---------------------------------------------------------------------------
# Coarse derivations
# ---------------------------------------------------------------------------


def coarse_type_array(signature: np.ndarray) -> np.ndarray:
    """Map a signature array to its uint8 :class:`RegionType` (exon > intron).

    Returns ``0`` (INTERGENIC) / ``1`` (INTRON) / ``2`` (EXON) per region: the region-type axis of the gDNA
    fragment-length pools, the same convention as the C++ accumulator's ``fl_pool_idx``.
    """
    sig = np.asarray(signature)
    has_exon = (sig & (BIT_EXON_POS | BIT_EXON_NEG)) != 0
    has_intron = (sig & (BIT_INTRON_POS | BIT_INTRON_NEG)) != 0
    out = np.full(sig.shape, int(RegionType.INTERGENIC), dtype=np.uint8)
    out[has_intron] = int(RegionType.INTRON)
    out[has_exon] = int(RegionType.EXON)  # exon wins over intron
    return out


def transcript_strand_class(signature: np.ndarray) -> np.ndarray:
    """Map a uint8 signature array to its int8 transcript-strand class.

    Returns one of ``{TS_NONE=0, TS_POS=1, TS_NEG=2, TS_AMBIG=3}`` per region.
    ``TS_NONE`` covers intergenic regions; ``TS_AMBIG`` covers regions with
    transcripts on both strands.
    """
    sig = np.asarray(signature)
    has_pos = (sig & (BIT_INTRON_POS | BIT_EXON_POS)) != 0
    has_neg = (sig & (BIT_INTRON_NEG | BIT_EXON_NEG)) != 0
    out = np.zeros(sig.shape, dtype=np.int8)
    out[has_pos & ~has_neg] = TS_POS
    out[~has_pos & has_neg] = TS_NEG
    out[has_pos & has_neg] = TS_AMBIG
    return out


def rna_active_strands(signature: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Per-strand RNA-active masks from a signature array — ``(pos, neg)``.

    RNA can be present wherever a transcript is, exon or intron (RNA inside an intron is RNA that has not spliced
    there), so the mask on strand ``s`` is the union of that strand's exon and intron bits. It is the solver's
    transcript-continuity gate ``free_s``: a region's mask is its own bits, and a boundary's is the AND of its two
    flanks, because both must carry the strand for an unspliced fragment to cross.
    """
    sig = np.asarray(signature)
    return (sig & (BIT_EXON_POS | BIT_INTRON_POS)) != 0, (
        sig & (BIT_EXON_NEG | BIT_INTRON_NEG)
    ) != 0


__all__ = [
    "BIT_INTRON_POS",
    "BIT_INTRON_NEG",
    "BIT_EXON_POS",
    "BIT_EXON_NEG",
    "N_SIGNATURES",
    "RegionType",
    "RegionStrand",
    "TS_NONE",
    "TS_POS",
    "TS_NEG",
    "TS_AMBIG",
    "coarse_type_array",
    "transcript_strand_class",
    "rna_active_strands",
]
