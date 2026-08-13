"""The deposit-behaviour digest's fixture, run against the SPECIFICATION rather than the C++.

``rigel.scan_cache.deposit_digest`` computes the cache key's behavioural half from the NATIVE
accumulator, because the cache holds what the production scanner deposited and the key must certify
that. This module is the same fixture over the executable specification, and it exists for two gates:

* the two must produce the SAME digest — otherwise the key certifies an artifact nothing writes;
* the specification's rule can be PERTURBED in-process, which the compiled one cannot, so this is how
  "a changed deposit rule moves the digest" is actually falsified.

⛔ The fixture is duplicated deliberately and the duplication is GATED: if the two drift, the agreement
test fails loudly. Importing the src fixture and re-running it here would test nothing, because both
sides would then share whatever mistake it contained.
"""

from __future__ import annotations

import numpy as np

#: Same geometry as ``rigel.scan_cache.deposit_digest``: two annotated junctions, a short region whose far
#: line a fragment may or may not reach, and regions wide enough that containment is reachable.
CUTS = [0, 60, 200, 260, 1000, 1060, 1120, 2000, 2400]
REGION_TYPES = [0, 2, 2, 1, 2, 2, 1, 2]
JUNCTIONS = [(0, 260, 1000, 1), (0, 1120, 2000, 1)]

#: Contained / one line / three lines / spliced-both-cross / spliced-neither-crosses / two junctions and
#: lines / one junction and one line. ⭐ Every branch of the deposit rule appears at least once, so a
#: change confined to any single branch still moves the digest.
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
    """The digest over the specification, byte-compatible with ``scan_cache.deposit_digest``."""
    from rigel.scan_cache import _digest, _schema_names

    partition = partition_cls.from_cuts([CUTS], region_types=[REGION_TYPES], junctions=JUNCTIONS)
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
