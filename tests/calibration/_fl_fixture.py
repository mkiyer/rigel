"""A hand-checkable substrate for the gDNA length law's two-pool contrast: five regions on one reference.

Intergenic 10 kb | exon 200 | intron 8 kb | exon 200 | intergenic 10 kb, every off-target region holding exactly its
uniform expectation at ``RHO`` per contained start, and the five pools at the laws that implies: the contained pools
and the spliced RNA pool at one law, the two crossing pools tilted by their crossing opportunity.

Deliberately NOT a real ``AccumulatorPayload``: the length laws consume four payload fields, and building the full
scanner object here would couple this fixture to the scan schema for no gain. The dataclass below carries exactly the
consumed surface, and a schema change that widens that surface will fail these tests loudly at the attribute, which is
the right failure.
"""

from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace

import numpy as np
import pandas as pd

from rigel.calibration.effective_length import contained_eff_length
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.signature import BIT_EXON_POS, BIT_INTRON_POS

MAX_SIZE = 100
RHO = 0.05


@dataclass(frozen=True)
class _PayloadLike:
    pool_lengths: np.ndarray
    region_contained_count: np.ndarray
    max_length: int
    deposited_lengths: np.ndarray


#: A coarse type's signature (0 intergenic / 1 intron / 2 exon): one strand's bit types a region.
_SIGNATURE = {0: 0, 1: BIT_INTRON_POS, 2: BIT_EXON_POS}


def regions_of(lengths, types, mask=()) -> RegionArrays:
    """One reference tiled by regions of these *lengths* and coarse *types*, with *mask*'s ``[start, end)`` runs
    unobservable: the frame `fl.py` reads, as `RegionArrays.from_index` builds it from an index."""
    lengths = np.asarray(lengths, dtype=np.int64)
    bounds = np.concatenate([[0], np.cumsum(lengths)])
    frame = pd.DataFrame(
        {
            "ref_name": "chr1",
            "start": bounds[:-1],
            "end": bounds[1:],
            "signature": [_SIGNATURE[int(t)] for t in types],
        }
    )
    runs = pd.DataFrame({"ref": "chr1", "start": [s for s, _ in mask], "end": [e for _, e in mask]})
    return RegionArrays.from_frame(frame, {"chr1": 0}, runs if len(runs) else None)


def build_fixture():
    """``(payload_like, opportunity_like, regions)``."""
    L = np.arange(MAX_SIZE + 1, dtype=np.float64)
    pmf = np.exp(-0.5 * ((L - 50.0) / 12.0) ** 2)
    pmf[0] = 0.0
    pmf /= pmf.sum()
    mu = float((pmf * L).sum())

    region_lengths = np.array([10_000.0, 200.0, 8_000.0, 200.0, 10_000.0])
    region_types = np.array([0, 2, 1, 2, 0], dtype=np.uint8)  # 0 ig / 1 intron / 2 exon
    regions = regions_of(region_lengths, region_types)
    e_contained = contained_eff_length(region_lengths, pmf)
    counts = np.zeros((5, 2), dtype=np.float64)
    for i in (0, 2, 4):  # the off-target regions carry exactly the uniform expectation
        counts[i] = RHO * e_contained[i] / 2.0
    counts[1] = counts[3] = 500.0  # exon counts are RNA-dominated and unread by the contrast

    crossing_tilt = pmf * np.clip(L - 1.0, 0.0, None)
    n_cross = (
        2.0 * RHO * (mu - 1.0)
    )  # two boundaries per crossing pool, each at the uniform expectation
    pool_lengths = np.zeros((5, MAX_SIZE + 1), dtype=np.float64)
    pool_lengths[0] = pmf * float(counts[0].sum() + counts[4].sum())  # intergenic contained
    pool_lengths[1] = pmf * float(counts[2].sum())  # intronic contained
    pool_lengths[2] = pool_lengths[3] = (
        crossing_tilt / crossing_tilt.sum() * n_cross
    )  # the two crossing pools
    pool_lengths[4] = pmf * 1000.0  # spliced RNA

    payload = _PayloadLike(
        pool_lengths=pool_lengths,
        region_contained_count=counts,
        max_length=MAX_SIZE,
        deposited_lengths=pool_lengths.sum(axis=0),
    )
    # the opportunity object: per-pool A(L) and the total T(L), toy-scaled
    total = np.full(MAX_SIZE + 1, 28_400.0)
    a_contained = np.clip(region_lengths[[0, 2, 4]].sum() - L + 1.0, 0.0, None)
    a_cross = np.full(MAX_SIZE + 1, 2.0) * np.clip(L - 1.0, 0.0, None)
    pools = (a_contained, a_contained * 0.4, a_cross, a_cross)
    opp = SimpleNamespace(pools=pools, total=total, combined_probability=lambda: sum(pools) / total)
    return payload, opp, regions
