"""Region signature encoding — the 4-bit exon/intron x strand annotation.

Every genome position carries a 4-bit *signature* recording which transcript features cover it:
``{intron_pos, intron_neg, exon_pos, exon_neg}``. A calibration *region* (see
:mod:`rigel.calibration.region_arrays`) is a maximal interval over which this signature is constant.

This module is the pure encoding layer: the bit constants, the ``RegionType`` / ``RegionStrand``
enums, and the small functions that derive a coarse class from a signature. It has no tunable
parameters — the bit values and enum integers are a wire encoding, not heuristics — and it depends
on nothing else in the package, so every other calibration module may import it freely.

The derivations here are the vocabulary the rest of calibration is written in: the coarse type
(exon wins over intron), the transcript-strand class, and the two per-strand activity masks that
say which RNA populations the annotation admits at a position.
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


# Transcript-strand class for a region's signature — the int8 encoding used by
# :class:`rigel.calibration.region_arrays.RegionArrays`. These are ALIASES of
# :class:`RegionStrand` (which mirrors :class:`rigel.types.Strand`): there is ONE
# strand-value convention across the codebase — NONE=0, POS=1, NEG=2, AMBIG=3.
# All consumers compare against these named constants (equality, no sign
# arithmetic), so a value in the Strand convention — e.g. the accumulator's
# per-boundary splice-junction strand (Strand POS=1 / NEG=2) — routes directly
# against TS_POS / TS_NEG with no conversion.
#
# NONE and AMBIG both lack a single transcript strand but are NOT
# interchangeable for the strand channel:
#   * TS_NONE  — no transcript (intergenic). gDNA is unstranded, so an arbitrary
#                sense assignment is safe (neutral). Stays in the strand model.
#   * TS_AMBIG — transcripts on BOTH strands (overlapping opposite-strand
#                annotations). Every read is sense for one and antisense for the
#                other, so there is no valid sense split. AMBIG regions are
#                excluded from strand deconvolution and recovered by density,
#                boundary-sweep imputation and the global fallback.
TS_NONE: int = int(RegionStrand.NONE)  # 0
TS_POS: int = int(RegionStrand.POS)  # 1
TS_NEG: int = int(RegionStrand.NEG)  # 2
TS_AMBIG: int = int(RegionStrand.AMBIG)  # 3


# ---------------------------------------------------------------------------
# Coarse derivations
# ---------------------------------------------------------------------------


def coarse_type_array(signature: np.ndarray) -> np.ndarray:
    """Map a signature array to its uint8 :class:`RegionType` (exon > intron).

    Returns ``0`` (INTERGENIC) / ``1`` (INTRON) / ``2`` (EXON) per region — the
    region-type axis of the gDNA fragment-length pools, and the same convention as
    the C++ accumulator's ``fl_pool_idx``.
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


def nrna_active_strands(signature: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Per-strand nascent-RNA-active masks from a signature array — ``(pos, neg)``.

    Nascent RNA is present wherever a transcript is (exon OR intron), so ``nrna_active`` on strand
    ``s`` is the union of that strand's exon and intron bits. This is the boundary solver's
    transcript-continuity gate ``free_s``: a region's ``nrna_active`` is its own bits; a boundary's
    is the AND of its two flanks, because both flanks must carry the strand for an unspliced
    fragment to cross.
    """
    sig = np.asarray(signature)
    return (sig & (BIT_EXON_POS | BIT_INTRON_POS)) != 0, (
        sig & (BIT_EXON_NEG | BIT_INTRON_NEG)
    ) != 0


def mrna_active_strands(signature: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Per-strand mature-RNA-active masks from a signature array — ``(pos, neg)``.

    Mature RNA lives only in exons, so ``mrna_active`` on strand ``s`` is that strand's exon bit,
    hence ``mrna_active_s`` implies ``nrna_active_s``. A region's ``mrna_active`` is its own exon
    bits; a boundary's is the AND of its two flanks, because an unspliced fragment crosses as mature
    only over contiguous exon — otherwise it would be spliced. This selects the region's solver
    prior: ``nrna_active`` without ``mrna_active`` means the annotation admits RNA here only inside
    an intron, which takes the near-zero nascent prior.
    """
    sig = np.asarray(signature)
    return (sig & BIT_EXON_POS) != 0, (sig & BIT_EXON_NEG) != 0


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
    "nrna_active_strands",
    "mrna_active_strands",
]
