"""The deposit-behaviour digest's fixture, run against the executable specification rather than the C++.

``deposit_digest`` in ``scripts/design/_scan_cache.py`` computes the cache key's behavioural half from the native
accumulator, because the key must certify what the production scanner deposited. This module runs the same fixture over
the specification, for two tests:

* the two must produce the same digest, or the key certifies an artifact nothing writes;
* the specification's rule can be perturbed in-process, which the compiled one cannot, so this is how a changed deposit
  rule is shown to move the digest.

The fixture is duplicated on purpose, and the agreement test catches drift; importing the cache module's fixture would let
both sides share any mistake in it.
"""

from __future__ import annotations

import numpy as np

#: Same geometry as ``_scan_cache.deposit_digest``: two annotated sj, a short region whose far
#: boundary a fragment may or may not reach, and regions wide enough that containment is reachable.
REGION_BOUNDS = [0, 60, 200, 260, 1000, 1060, 1120, 2000, 2400]
REGION_TYPES = [0, 2, 2, 1, 2, 2, 1, 2]
SJ = [(0, 260, 1000, 1), (0, 1120, 2000, 1)]

#: Contained / one boundary / three boundaries / spliced-both-cross / spliced-neither-crosses / two sj
#: and boundaries / one sj and one boundary. Every branch of the deposit rule appears at least once, so
#: a change confined to any single branch still moves the digest.
FRAGMENTS = (
    (10, 50, ()),
    (30, 150, ()),
    (30, 280, ()),
    (150, 1060, ((260, 1000),)),
    (210, 1040, ((260, 1000),)),
    (150, 2100, ((260, 1000), (1120, 2000))),
    (1030, 2050, ((1120, 2000),)),
)


def reference_deposit_digest(accumulator_cls, partition_cls) -> str:
    """The digest over the specification, byte-compatible with ``_scan_cache.deposit_digest``."""
    from _scan_cache import _digest, _schema_names

    partition = partition_cls.from_region_bounds(
        [REGION_BOUNDS], region_types=[REGION_TYPES], sj=SJ
    )
    accumulator = accumulator_cls(partition, max_fragment_length=1000)
    for start, end, introns in FRAGMENTS:
        accumulator.deposit(0, start, end, observed_introns=introns, sj_strand=1)

    tally = accumulator.tally
    parts: list[bytes] = []
    for name in sorted(_schema_names()):
        value = getattr(tally, name.split("__")[0], None)
        if isinstance(value, np.ndarray):
            parts.append(name.encode())
            parts.append(np.ascontiguousarray(value).tobytes())
    return _digest(*parts)
