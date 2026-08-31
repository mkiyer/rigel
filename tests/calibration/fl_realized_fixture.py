"""A hand-checkable substrate for the realized gDNA law: five regions, four boundaries, one knob.

The knob is ``boundary_excess`` — how many times the uniform-field expectation each exon-flanking
boundary's unspliced count carries. 1.0 says "no capture" (the on-target correction must then vanish
identically); large values say "strong capture".

⚠ Deliberately NOT a real ``AccumulatorPayload``: the realized-law machinery consumes six payload
fields, and building the full scanner object here would couple this fixture to the scan schema for no
gain. The dataclass below carries exactly the consumed surface, and a schema change that widens that
surface will fail these tests loudly at the attribute, which is the right failure.
"""

from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace

import numpy as np

from rigel.calibration.gdna_density import contained_opportunity

MAX_SIZE = 100
RHO = 0.05


@dataclass(frozen=True)
class _PayloadLike:
    pool_lengths: np.ndarray
    region_contained_count: np.ndarray
    boundary_unspliced_count: np.ndarray
    ref_region_offsets: np.ndarray
    ref_boundary_offsets: np.ndarray
    max_length: int


def _pmf() -> np.ndarray:
    L = np.arange(MAX_SIZE + 1, dtype=np.float64)
    p = np.exp(-0.5 * ((L - 50.0) / 12.0) ** 2)
    p[0] = 0.0
    return p / p.sum()


def build_fixture(boundary_excess: float):
    """(payload_like, opp_like, region_lengths, region_types, rna_pmf, uniform_counts)."""
    pmf = _pmf()
    mu = float((pmf * np.arange(pmf.size)).sum())

    # regions: intergenic 10k | exon 200 | intron 8k | exon 200 | intergenic 10k
    region_lengths = np.array([10_000.0, 200.0, 8_000.0, 200.0, 10_000.0])
    region_types = np.array([0, 2, 1, 2, 0], dtype=np.uint8)  # 0 ig / 1 intron / 2 exon

    e_contained = contained_opportunity(pmf, region_lengths)
    counts = np.zeros((5, 2), dtype=np.float64)
    for i in (0, 2, 4):  # the off-target regions carry exactly the uniform expectation
        counts[i] = RHO * e_contained[i] / 2.0
    counts[1] = counts[3] = 500.0  # exon counts are RNA-dominated and unread by the machinery

    # boundaries: ig|exon, exon|intron, intron|exon, exon|ig — each at `boundary_excess` x uniform
    per_boundary = boundary_excess * RHO * (mu - 1.0)
    boundary_counts = np.full((4, 2), per_boundary / 2.0, dtype=np.float64)

    L = np.arange(MAX_SIZE + 1, dtype=np.float64)
    crossing_tilt = pmf * np.clip(L - 1.0, 0.0, None)
    pool_lengths = np.zeros((5, MAX_SIZE + 1), dtype=np.float64)
    pool_lengths[0] = pmf * float(counts[0].sum() + counts[4].sum())  # intergenic contained
    pool_lengths[1] = pmf * float(counts[2].sum())  # intronic contained
    n_cross = 2.0 * per_boundary
    pool_lengths[2] = crossing_tilt / crossing_tilt.sum() * n_cross  # intron|exon
    pool_lengths[3] = crossing_tilt / crossing_tilt.sum() * n_cross  # intergenic|exon
    pool_lengths[4] = pmf * 1000.0  # spliced RNA, unread here

    payload = _PayloadLike(
        pool_lengths=pool_lengths,
        region_contained_count=counts,
        boundary_unspliced_count=boundary_counts,
        ref_region_offsets=np.array([0, 5], dtype=np.int64),
        ref_boundary_offsets=np.array([0, 4], dtype=np.int64),
        max_length=MAX_SIZE,
    )

    # the opportunity object: per-pool A(L) and the total T(L), toy-scaled
    total = np.full(MAX_SIZE + 1, 28_400.0)
    a_contained = np.clip(region_lengths[[0, 2, 4]].sum() - L + 1.0, 0.0, None)
    a_cross = np.full(MAX_SIZE + 1, 2.0) * np.clip(L - 1.0, 0.0, None)
    opp = SimpleNamespace(pools=(a_contained, a_contained * 0.4, a_cross, a_cross), total=total)

    uniform_counts = pmf * float(pool_lengths[:4].sum())
    return payload, opp, region_lengths, region_types, pmf.copy(), uniform_counts
