"""A hand-checkable substrate for the realized gDNA law: five regions, four boundaries, one knob.

The knob is ``boundary_excess`` — how many times the uniform-field expectation each exon-flanking
boundary's unspliced count carries. 1.0 says "no capture" (the on-target correction must then vanish
identically); large values say "strong capture". ``rna_rate`` adds RNA to the intron and to the two
boundaries it flanks, at the expected counts: ``rna_rate`` per unit of the intron's contained opportunity,
and ``rna_rate·(mu - 1)`` crossing each boundary, the same crossing count the estimator's odds assume.

Deliberately NOT a real ``AccumulatorPayload``: the realized-law machinery consumes six payload
fields, and building the full scanner object here would couple this fixture to the scan schema for no
gain. The dataclass below carries exactly the consumed surface, and a schema change that widens that
surface will fail these tests loudly at the attribute, which is the right failure.
"""

from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace

import numpy as np
import pandas as pd

from rigel.calibration.gdna_density import contained_opportunity
from rigel.calibration.region_arrays import RegionArrays, boundary_reaches, observable_pieces
from rigel.calibration.signature import BIT_EXON_POS, BIT_INTRON_POS

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


def _pmf() -> np.ndarray:
    L = np.arange(MAX_SIZE + 1, dtype=np.float64)
    p = np.exp(-0.5 * ((L - 50.0) / 12.0) ** 2)
    p[0] = 0.0
    return p / p.sum()


def build_fixture(boundary_excess: float, rna_rate: float = 0.0, mask=()):
    """(payload_like, opp_like, regions, rna_pmf, uniform_counts). *mask* holds the genome's unobservable
    ``[start, end)`` runs: every count is over the observable starts alone, and a boundary with a masked base
    beside it is crossed by nothing observable. A run within a fragment of a boundary but not beside it would bind
    the boundary's crossing support and is refused."""
    pmf = _pmf()
    mu = float((pmf * np.arange(pmf.size)).sum())

    # regions: intergenic 10k | exon 200 | intron 8k | exon 200 | intergenic 10k
    region_lengths = np.array([10_000.0, 200.0, 8_000.0, 200.0, 10_000.0])
    region_types = np.array([0, 2, 1, 2, 0], dtype=np.uint8)  # 0 ig / 1 intron / 2 exon
    regions = regions_of(region_lengths, region_types, mask)
    reach_lo, reach_hi = boundary_reaches(regions)
    beside = (reach_lo == 0.0) | (reach_hi == 0.0)
    if np.any(~beside & ((reach_lo < MAX_SIZE) | (reach_hi < MAX_SIZE))):
        raise ValueError("a masked run within a fragment of a boundary binds its crossing support")

    piece_region, piece_len = observable_pieces(regions)
    e_contained = np.bincount(
        piece_region, weights=contained_opportunity(pmf, piece_len), minlength=5
    )
    counts = np.zeros((5, 2), dtype=np.float64)
    for i in (0, 2, 4):  # the off-target regions carry exactly the uniform expectation
        counts[i] = RHO * e_contained[i] / 2.0
    counts[1] = counts[3] = 500.0  # exon counts are RNA-dominated and unread by the machinery
    # the RNA law is the gDNA law here, so the intron has one contained opportunity for both
    rna_in_intron = rna_rate * e_contained[2]
    counts[2] += rna_in_intron / 2.0

    # boundaries: ig|exon, exon|intron, intron|exon, exon|ig — each at `boundary_excess` x uniform
    per_boundary = boundary_excess * RHO * (mu - 1.0)
    boundary_counts = np.full((4, 2), per_boundary / 2.0, dtype=np.float64)
    rna_crossing = rna_rate * (mu - 1.0)
    # exon|intron and intron|exon: the intron's RNA crosses the two boundaries it flanks
    boundary_counts[[1, 2]] += rna_crossing / 2.0
    boundary_counts[beside] = 0.0

    L = np.arange(MAX_SIZE + 1, dtype=np.float64)
    crossing_tilt = pmf * np.clip(L - 1.0, 0.0, None)
    pool_lengths = np.zeros((5, MAX_SIZE + 1), dtype=np.float64)
    pool_lengths[0] = pmf * float(counts[0].sum() + counts[4].sum())  # intergenic contained
    pool_lengths[1] = pmf * float(counts[2].sum())  # intronic contained
    n_cross = 2.0 * per_boundary
    # intron|exon, its gDNA and the intron's RNA crossings
    pool_lengths[2] = crossing_tilt / crossing_tilt.sum() * (n_cross + 2.0 * rna_crossing)
    pool_lengths[3] = crossing_tilt / crossing_tilt.sum() * n_cross  # intergenic|exon
    pool_lengths[4] = pmf * 1000.0  # spliced RNA, unread here

    payload = _PayloadLike(
        pool_lengths=pool_lengths,
        region_contained_count=counts,
        boundary_unspliced_count=boundary_counts,
        ref_region_offsets=np.array([0, 5], dtype=np.int64),
        ref_boundary_offsets=np.array([0, 4], dtype=np.int64),
        max_length=MAX_SIZE,
        deposited_lengths=pool_lengths.sum(axis=0),
    )

    # the opportunity object: per-pool A(L) and the total T(L), toy-scaled
    total = np.full(MAX_SIZE + 1, 28_400.0)
    a_contained = np.clip(region_lengths[[0, 2, 4]].sum() - L + 1.0, 0.0, None)
    a_cross = np.full(MAX_SIZE + 1, 2.0) * np.clip(L - 1.0, 0.0, None)
    pools = (a_contained, a_contained * 0.4, a_cross, a_cross)
    opp = SimpleNamespace(pools=pools, total=total, combined_probability=lambda: sum(pools) / total)

    uniform_counts = pmf * float(pool_lengths[:4].sum())
    return payload, opp, regions, pmf.copy(), uniform_counts
