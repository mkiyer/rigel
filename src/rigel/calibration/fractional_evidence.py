"""Fractional evidence — the single Python interface to the new payload.

Owns:

* Signature-derived per-region masks (``is_intergenic``, ``is_intron_only``,
  ``is_exon_any``, ``has_both_strands``).
* ``transcript_strand_class(signature)`` int8 column (TS_NONE / TS_POS /
  TS_NEG / TS_AMBIG).
* ``sense_antisense_split(...)`` and the ``SenseSplit`` dataclass.
* ``FractionalEvidenceView`` — cached per-payload helpers.
* FL pool helpers: ``gdna_fl_mass``, ``rna_candidate_fl_mass``, ``pool``,
  ``pool_total``.

No legacy ``ORIENT_*`` constants. Strand-relative analysis goes through
``sense_antisense_split`` exclusively.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from .signature import (
    BIT_EXON_NEG,
    BIT_EXON_POS,
    BIT_INTRON_NEG,
    BIT_INTRON_POS,
    CHANNEL_STRAND_NEG,
    CHANNEL_STRAND_POS,
    COMPARTMENT_BOUNDARY_LEFT,
    COMPARTMENT_BOUNDARY_RIGHT,
    COMPARTMENT_CONTAINED,
    FL_POOL_EXONIC_BOUNDARY,
    FL_POOL_EXONIC_CONTAINED,
    FL_POOL_INTERGENIC_BOUNDARY,
    FL_POOL_INTERGENIC_CONTAINED,
    FL_POOL_INTRONIC_BOUNDARY,
    FL_POOL_INTRONIC_CONTAINED,
    FL_POOL_NAMES,
    N_COMPARTMENTS,
    N_FL_POOLS,
    N_SIGNATURES,
    N_SPLICE_STATES,
    SPLICE_SPLICED,
    SPLICE_UNSPLICED,
    channel_index,
)

if TYPE_CHECKING:
    from .scan_payload import CalibrationScanPayload


# Transcript-strand class for a region's signature.
TS_NONE: int = 0
TS_POS: int = 1
TS_NEG: int = -1
TS_AMBIG: int = 2


# ---------------------------------------------------------------------------
# Per-region signature derivations
# ---------------------------------------------------------------------------


def transcript_strand_class(signature: np.ndarray) -> np.ndarray:
    """Map a uint8 signature array to int8 transcript-strand class.

    Returns one of {TS_NONE=0, TS_POS=+1, TS_NEG=-1, TS_AMBIG=2} per region.
    """
    sig = np.asarray(signature)
    has_pos = (sig & (BIT_INTRON_POS | BIT_EXON_POS)) != 0
    has_neg = (sig & (BIT_INTRON_NEG | BIT_EXON_NEG)) != 0
    out = np.zeros(sig.shape, dtype=np.int8)
    out[has_pos & ~has_neg] = TS_POS
    out[~has_pos & has_neg] = TS_NEG
    out[has_pos & has_neg] = TS_AMBIG
    return out


def is_intergenic(signature: np.ndarray) -> np.ndarray:
    """Regions with no exon or intron bit set."""
    return np.asarray(signature) == 0


def is_intron_only(signature: np.ndarray) -> np.ndarray:
    """Regions with intron bits set and no exon bits."""
    sig = np.asarray(signature)
    has_intron = (sig & (BIT_INTRON_POS | BIT_INTRON_NEG)) != 0
    has_exon = (sig & (BIT_EXON_POS | BIT_EXON_NEG)) != 0
    return has_intron & ~has_exon


def is_exon_any(signature: np.ndarray) -> np.ndarray:
    """Regions with any exon bit set (exon-only or mixed exon+intron)."""
    return (np.asarray(signature) & (BIT_EXON_POS | BIT_EXON_NEG)) != 0


def has_both_strands(signature: np.ndarray) -> np.ndarray:
    """Regions with both transcript strands represented (ambiguous strand)."""
    sig = np.asarray(signature)
    has_pos = (sig & (BIT_INTRON_POS | BIT_EXON_POS)) != 0
    has_neg = (sig & (BIT_INTRON_NEG | BIT_EXON_NEG)) != 0
    return has_pos & has_neg


# ---------------------------------------------------------------------------
# Sense / antisense split (the only orient API)
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class SenseSplit:
    """Per-region (sense, antisense, ambiguous) float32 mass triple."""

    sense: np.ndarray  # float32[R]
    antisense: np.ndarray  # float32[R]
    ambiguous: np.ndarray  # float32[R]


def sense_antisense_split(
    region_counts: np.ndarray,
    ts_class: np.ndarray,
    *,
    compartment: int,
    splice: int,
) -> SenseSplit:
    """Return per-region (sense, antisense, ambiguous) float32 mass.

    For each region:

    * ``TS_POS``: sense = mass on POS channel; antisense = mass on NEG channel.
    * ``TS_NEG``: sense = mass on NEG channel; antisense = mass on POS channel.
    * ``TS_AMBIG`` or ``TS_NONE``: both channels go to ambiguous; sense and
      antisense are zero.
    """
    if compartment < 0 or compartment >= N_COMPARTMENTS:
        raise ValueError(f"compartment must be in [0, {N_COMPARTMENTS}); got {compartment}")
    if splice < 0 or splice >= N_SPLICE_STATES:
        raise ValueError(f"splice must be in [0, {N_SPLICE_STATES}); got {splice}")
    pos = region_counts[:, channel_index(compartment, splice, CHANNEL_STRAND_POS)]
    neg = region_counts[:, channel_index(compartment, splice, CHANNEL_STRAND_NEG)]
    ts = np.asarray(ts_class)

    sense = np.zeros(pos.shape, dtype=np.float32)
    antisense = np.zeros(pos.shape, dtype=np.float32)
    ambiguous = np.zeros(pos.shape, dtype=np.float32)

    is_pos = ts == TS_POS
    is_neg = ts == TS_NEG
    is_amb = (ts == TS_AMBIG) | (ts == TS_NONE)

    sense[is_pos] = pos[is_pos]
    antisense[is_pos] = neg[is_pos]
    sense[is_neg] = neg[is_neg]
    antisense[is_neg] = pos[is_neg]
    ambiguous[is_amb] = pos[is_amb] + neg[is_amb]

    return SenseSplit(sense=sense, antisense=antisense, ambiguous=ambiguous)


# ---------------------------------------------------------------------------
# FL pool helpers
# ---------------------------------------------------------------------------


_POOL_BY_NAME: dict[str, int] = {name: idx for idx, name in enumerate(FL_POOL_NAMES)}


def pool(payload: "CalibrationScanPayload", name: str) -> np.ndarray:
    """Return the float64[1024] FL pool by name."""
    idx = _POOL_BY_NAME.get(name)
    if idx is None:
        raise ValueError(f"unknown FL pool name {name!r}; expected one of {FL_POOL_NAMES}")
    return payload.fl_pool_mass[idx]


def pool_total(payload: "CalibrationScanPayload", name: str) -> float:
    """Return the float64 total mass of the named FL pool."""
    idx = _POOL_BY_NAME.get(name)
    if idx is None:
        raise ValueError(f"unknown FL pool name {name!r}; expected one of {FL_POOL_NAMES}")
    return float(payload.fl_pool_total[idx])


def gdna_fl_mass(payload: "CalibrationScanPayload") -> np.ndarray:
    """Aggregate gDNA FL mass: INTERGENIC + INTRONIC pools (both compartments).

    Returns float64[1024].
    """
    fp = payload.fl_pool_mass
    return (
        fp[FL_POOL_INTERGENIC_CONTAINED]
        + fp[FL_POOL_INTERGENIC_BOUNDARY]
        + fp[FL_POOL_INTRONIC_CONTAINED]
        + fp[FL_POOL_INTRONIC_BOUNDARY]
    ).astype(np.float64, copy=False)


def rna_candidate_fl_mass(payload: "CalibrationScanPayload") -> np.ndarray:
    """Aggregate likely-RNA FL mass: EXONIC pools (both compartments).

    Returns float64[1024]. Diagnostic only; the trained RNA FL model lives
    on ``scan_trained.category_models[SPLICED_ANNOT]``.
    """
    fp = payload.fl_pool_mass
    return (fp[FL_POOL_EXONIC_CONTAINED] + fp[FL_POOL_EXONIC_BOUNDARY]).astype(
        np.float64, copy=False
    )


# ---------------------------------------------------------------------------
# FractionalEvidenceView — cached per-payload helpers
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class FractionalEvidenceView:
    """Cached per-payload helpers keyed on the payload + region signatures.

    The view holds references (no copies) to the payload's
    ``region_counts`` and to the int8 ``ts_class`` derived from the
    region signatures. Channel-slice helpers return float32 views into
    ``region_counts``.
    """

    payload: "CalibrationScanPayload"
    signature: np.ndarray  # uint8[R]
    ts_class: np.ndarray  # int8[R]
    mask_intergenic: np.ndarray  # bool[R]
    mask_intron_only: np.ndarray  # bool[R]
    mask_exon_any: np.ndarray  # bool[R]
    mask_both_strands: np.ndarray  # bool[R]

    @classmethod
    def from_payload(
        cls,
        payload: "CalibrationScanPayload",
        signature: np.ndarray,
    ) -> "FractionalEvidenceView":
        sig = np.ascontiguousarray(signature, dtype=np.uint8)
        if sig.shape != (payload.n_regions,):
            raise ValueError(
                f"signature shape {sig.shape} does not match payload n_regions "
                f"({payload.n_regions})."
            )
        if np.any(sig >= N_SIGNATURES):
            raise ValueError("signature contains values outside [0, 15].")
        return cls(
            payload=payload,
            signature=sig,
            ts_class=transcript_strand_class(sig),
            mask_intergenic=is_intergenic(sig),
            mask_intron_only=is_intron_only(sig),
            mask_exon_any=is_exon_any(sig),
            mask_both_strands=has_both_strands(sig),
        )

    # ----- channel slices ---------------------------------------------------

    def channel(self, compartment: int, splice: int, strand: int) -> np.ndarray:
        """Return float32[R] view of the named channel."""
        return self.payload.region_counts[:, channel_index(compartment, splice, strand)]

    def contained_unspliced(self, strand: int | None = None) -> np.ndarray:
        return self._compartment_view(COMPARTMENT_CONTAINED, SPLICE_UNSPLICED, strand)

    def boundary_left_unspliced(self, strand: int | None = None) -> np.ndarray:
        return self._compartment_view(COMPARTMENT_BOUNDARY_LEFT, SPLICE_UNSPLICED, strand)

    def boundary_right_unspliced(self, strand: int | None = None) -> np.ndarray:
        return self._compartment_view(COMPARTMENT_BOUNDARY_RIGHT, SPLICE_UNSPLICED, strand)

    def contained_spliced(self, strand: int | None = None) -> np.ndarray:
        return self._compartment_view(COMPARTMENT_CONTAINED, SPLICE_SPLICED, strand)

    def region_total_mass(self) -> np.ndarray:
        """Per-region sum across all 12 channels (float32[R])."""
        return self.payload.region_counts.sum(axis=1)

    def _compartment_view(self, compartment: int, splice: int, strand: int | None) -> np.ndarray:
        if strand is None:
            pos = self.payload.region_counts[
                :, channel_index(compartment, splice, CHANNEL_STRAND_POS)
            ]
            neg = self.payload.region_counts[
                :, channel_index(compartment, splice, CHANNEL_STRAND_NEG)
            ]
            return pos + neg
        return self.channel(compartment, splice, strand)

    # ----- sense / antisense ------------------------------------------------

    def split(self, *, compartment: int, splice: int = SPLICE_UNSPLICED) -> SenseSplit:
        """Per-region (sense, antisense, ambiguous) mass for the given channel slice."""
        return sense_antisense_split(
            self.payload.region_counts,
            self.ts_class,
            compartment=compartment,
            splice=splice,
        )


__all__ = [
    "TS_NONE",
    "TS_POS",
    "TS_NEG",
    "TS_AMBIG",
    "transcript_strand_class",
    "is_intergenic",
    "is_intron_only",
    "is_exon_any",
    "has_both_strands",
    "SenseSplit",
    "sense_antisense_split",
    "pool",
    "pool_total",
    "gdna_fl_mass",
    "rna_candidate_fl_mass",
    "FractionalEvidenceView",
    "N_FL_POOLS",
    "FL_POOL_NAMES",
]
