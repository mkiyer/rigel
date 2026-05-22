"""Fine-region signature and fractional-channel layout helpers."""

from __future__ import annotations

from .regions import RegionStrand, RegionType

BIT_INTRON_POS = 0x8
BIT_INTRON_NEG = 0x4
BIT_EXON_POS = 0x2
BIT_EXON_NEG = 0x1
N_SIGNATURES = 16

COMPARTMENT_CONTAINED = 0
COMPARTMENT_BOUNDARY_LEFT = 1
COMPARTMENT_BOUNDARY_RIGHT = 2
N_COMPARTMENTS = 3

SPLICE_UNSPLICED = 0
SPLICE_SPLICED = 1
N_SPLICE_STATES = 2

CHANNEL_STRAND_POS = 0
CHANNEL_STRAND_NEG = 1
N_CHANNEL_STRANDS = 2

CHAN_CONTAINED = COMPARTMENT_CONTAINED
CHAN_BOUNDARY_LEFT = COMPARTMENT_BOUNDARY_LEFT
CHAN_BOUNDARY_RIGHT = COMPARTMENT_BOUNDARY_RIGHT
N_CHANNELS = N_COMPARTMENTS * N_SPLICE_STATES * N_CHANNEL_STRANDS

# FL pool enum (mirrored in src/rigel/native/calibration/region_signature.h).
# Indexed via fl_pool_index(signature, compartment). Boundary_left and
# boundary_right collapse to a single BOUNDARY pool: the receiving
# region's coarse class determines the pool, and per-side discrimination
# already happened via mass routing.
FL_POOL_INTERGENIC_CONTAINED = 0
FL_POOL_INTERGENIC_BOUNDARY = 1
FL_POOL_INTRONIC_CONTAINED = 2
FL_POOL_INTRONIC_BOUNDARY = 3
FL_POOL_EXONIC_CONTAINED = 4
FL_POOL_EXONIC_BOUNDARY = 5
N_FL_POOLS = 6

FL_POOL_NAMES: tuple[str, ...] = (
    "INTERGENIC_CONTAINED",
    "INTERGENIC_BOUNDARY",
    "INTRONIC_CONTAINED",
    "INTRONIC_BOUNDARY",
    "EXONIC_CONTAINED",
    "EXONIC_BOUNDARY",
)


__all__ = [
    "BIT_INTRON_POS",
    "BIT_INTRON_NEG",
    "BIT_EXON_POS",
    "BIT_EXON_NEG",
    "N_SIGNATURES",
    "COMPARTMENT_CONTAINED",
    "COMPARTMENT_BOUNDARY_LEFT",
    "COMPARTMENT_BOUNDARY_RIGHT",
    "N_COMPARTMENTS",
    "SPLICE_UNSPLICED",
    "SPLICE_SPLICED",
    "N_SPLICE_STATES",
    "CHANNEL_STRAND_POS",
    "CHANNEL_STRAND_NEG",
    "N_CHANNEL_STRANDS",
    "CHAN_CONTAINED",
    "CHAN_BOUNDARY_LEFT",
    "CHAN_BOUNDARY_RIGHT",
    "N_CHANNELS",
    "FL_POOL_INTERGENIC_CONTAINED",
    "FL_POOL_INTERGENIC_BOUNDARY",
    "FL_POOL_INTRONIC_CONTAINED",
    "FL_POOL_INTRONIC_BOUNDARY",
    "FL_POOL_EXONIC_CONTAINED",
    "FL_POOL_EXONIC_BOUNDARY",
    "N_FL_POOLS",
    "FL_POOL_NAMES",
    "pack_signature",
    "validate_signature",
    "has_intron_flag",
    "has_exon_flag",
    "coarse_type_from_signature",
    "coarse_strand_from_signature",
    "is_ambiguous_signature",
    "channel_index",
    "channel_tuple",
    "fl_pool_index",
]


def pack_signature(
    *,
    intron_pos: bool = False,
    intron_neg: bool = False,
    exon_pos: bool = False,
    exon_neg: bool = False,
) -> int:
    """Pack four fine-region flags into the canonical uint8 signature."""
    signature = 0
    if intron_pos:
        signature |= BIT_INTRON_POS
    if intron_neg:
        signature |= BIT_INTRON_NEG
    if exon_pos:
        signature |= BIT_EXON_POS
    if exon_neg:
        signature |= BIT_EXON_NEG
    return signature


def validate_signature(signature: int) -> int:
    """Return ``signature`` as int after validating the 4-bit range."""
    value = int(signature)
    if value < 0 or value >= N_SIGNATURES:
        raise ValueError(f"signature must be in [0, 15]; got {signature!r}")
    return value


def has_intron_flag(signature: int) -> bool:
    """Return whether any intron flag is present."""
    value = validate_signature(signature)
    return bool(value & (BIT_INTRON_POS | BIT_INTRON_NEG))


def has_exon_flag(signature: int) -> bool:
    """Return whether any exon flag is present."""
    value = validate_signature(signature)
    return bool(value & (BIT_EXON_POS | BIT_EXON_NEG))


def coarse_type_from_signature(signature: int) -> int:
    """Derive the coarse region type from a fine signature.

    Exon wins over intron, matching the current coarse partition.
    """
    value = validate_signature(signature)
    if value & (BIT_EXON_POS | BIT_EXON_NEG):
        return int(RegionType.EXON)
    if value & (BIT_INTRON_POS | BIT_INTRON_NEG):
        return int(RegionType.INTRON)
    return int(RegionType.INTERGENIC)


def coarse_strand_from_signature(signature: int) -> int:
    """Derive coarse strand flags from all active fine signature bits."""
    value = validate_signature(signature)
    has_pos = bool(value & (BIT_INTRON_POS | BIT_EXON_POS))
    has_neg = bool(value & (BIT_INTRON_NEG | BIT_EXON_NEG))
    if has_pos and has_neg:
        return int(RegionStrand.AMBIG)
    if has_pos:
        return int(RegionStrand.POS)
    if has_neg:
        return int(RegionStrand.NEG)
    return int(RegionStrand.NONE)


def is_ambiguous_signature(signature: int) -> bool:
    """Return whether the derived coarse strand is ambiguous."""
    return coarse_strand_from_signature(signature) == int(RegionStrand.AMBIG)


def channel_index(compartment: int, splice_idx: int, strand_idx: int) -> int:
    """Pack ``(compartment, splice, strand)`` into the 12-channel index."""
    compartment = int(compartment)
    splice_idx = int(splice_idx)
    strand_idx = int(strand_idx)
    if compartment < 0 or compartment >= N_COMPARTMENTS:
        raise ValueError(f"compartment must be in [0, {N_COMPARTMENTS}); got {compartment}")
    if splice_idx < 0 or splice_idx >= N_SPLICE_STATES:
        raise ValueError(f"splice_idx must be in [0, {N_SPLICE_STATES}); got {splice_idx}")
    if strand_idx < 0 or strand_idx >= N_CHANNEL_STRANDS:
        raise ValueError(f"strand_idx must be in [0, {N_CHANNEL_STRANDS}); got {strand_idx}")
    return compartment * (N_SPLICE_STATES * N_CHANNEL_STRANDS) + splice_idx * 2 + strand_idx


def channel_tuple(channel: int) -> tuple[int, int, int]:
    """Unpack a 12-channel index into ``(compartment, splice, strand)``."""
    channel = int(channel)
    if channel < 0 or channel >= N_CHANNELS:
        raise ValueError(f"channel must be in [0, {N_CHANNELS}); got {channel}")
    compartment = channel // (N_SPLICE_STATES * N_CHANNEL_STRANDS)
    rem = channel % (N_SPLICE_STATES * N_CHANNEL_STRANDS)
    splice_idx = rem // N_CHANNEL_STRANDS
    strand_idx = rem % N_CHANNEL_STRANDS
    return compartment, splice_idx, strand_idx


def fl_pool_index(signature: int, compartment: int) -> int:
    """Map ``(signature, compartment)`` to one of the six FL pools.

    Boundary_left and boundary_right collapse to a single BOUNDARY pool;
    contained stays distinct. The receiving region's coarse class
    (derived from ``signature``) selects the row:

    - INTERGENIC (``signature == 0x0``) → pools 0 (contained), 1 (boundary)
    - INTRON (``{0x4, 0x8, 0xC}``)      → pools 2 (contained), 3 (boundary)
    - EXON (any signature with an exon bit) → pools 4 (contained), 5 (boundary)
    """
    coarse = coarse_type_from_signature(signature)
    compartment = int(compartment)
    if compartment < 0 or compartment >= N_COMPARTMENTS:
        raise ValueError(f"compartment must be in [0, {N_COMPARTMENTS}); got {compartment}")
    contained_offset = 0 if compartment == COMPARTMENT_CONTAINED else 1
    return int(coarse) * 2 + contained_offset
