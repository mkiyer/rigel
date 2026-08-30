"""End-to-end: calibrate() recovers a planted gDNA strand overdispersion (docs/em_strand/03 §4).

Builds a payload of many INTRON regions of one + gene (pure gDNA in them) whose contained
sense/antisense counts are drawn from BetaBinom(½, od_true), runs the real calibrator, and checks the
fitted gdna_strand_overdispersion. ⚠ Intron, not intergenic: since 2026-08-29 the fit is the away-half
moment over GENIC seeds and an intergenic region — with no gene strand to orient by — is not a seed. This exercises the full Phase-2 path: substrate → count clue →
seed extraction → fit → result.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from _synthetic import make_gdna_fl_pmf, make_strand_models

from rigel.calibration import calibrate
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.signature import BIT_INTRON_POS
from rigel.config import CalibrationConfig
from rigel.scan_payload import (
    N_FRAGMENT_POOLS,
    AccumulatorPayload,
    DeferredFragments,
    GapCensus,
    ScanQC,
)

_KAPPA = 0.95


def _consistent_strand_model(overdispersion, n_sj=800, depth=60, seed=7):
    """A strand model whose SPLICED seeds carry the SAME overdispersion as the planted gDNA one.

    ⛔ **A library has ONE technical process, and since 2026-08-30 calibrate reconciles the two components
    against each other** (`gdna_strand.reconcile_overdispersions`), so a fixture that plants `od_g = 0.2`
    beside a single junction sitting exactly at κ — which reads `od_r = 0` with real information — is
    asserting a library that cannot exist, and the reconciliation correctly splits the difference. Planting
    the same dispersion in both components is what makes this an end-to-end RECOVERY test again.
    """
    from rigel.strand_model import SJStrandTable, StrandModel, StrandModels
    from rigel.types import Strand

    rng = np.random.default_rng(seed)
    if overdispersion > 0:
        conc = (1.0 - overdispersion) / overdispersion
        p = rng.beta(_KAPPA * conc, (1.0 - _KAPPA) * conc, size=n_sj)
    else:
        p = np.full(n_sj, _KAPPA)
    n_sense = rng.binomial(depth, p).astype(np.int64)
    table = SJStrandTable(
        ref_id=np.zeros(n_sj, dtype=np.int32),
        start=np.arange(n_sj, dtype=np.int64) * 1000,
        end=np.arange(n_sj, dtype=np.int64) * 1000 + 100,
        motif_strand=np.full(n_sj, int(Strand.POS), dtype=np.int8),
        n_sense=n_sense,
        n_antisense=np.full(n_sj, depth, dtype=np.int64) - n_sense,
    )
    return StrandModels(exonic_spliced=StrandModel.from_sj_table(table))


_STRAND_MODEL = make_strand_models(_KAPPA, 200)
_FRAG_LEN = 50  # the delta the gDNA pmf sits at; every fixture fragment is this long


def _intron_betabinom_payload(n_regions, depth, overdispersion, seed):
    """A 1-reference payload of ``n_regions`` intron(+) regions; contained gDNA ~ BetaBinom(½, od).

    ⚠ **The boundary axis is empty of COUNTS but not of ROWS.** One reference with ``k`` regions owns
    ``k − 1`` boundaries, and the payload must carry them or the chain builder refuses it. Leaving the
    counts at zero is what makes this test isolate the CONTAINED-region seed arm: with no crossing
    fragments there are no boundary seeds, so a recovered overdispersion can only have come from the regions.
    """
    rng = np.random.default_rng(seed)
    a = 0.5 * (1.0 - overdispersion) / overdispersion if overdispersion > 0 else 1e9
    p = rng.beta(a, a, size=n_regions)
    pos = rng.binomial(depth, p)
    neg = depth - pos

    contained = np.stack([pos, neg], axis=1).astype(np.uint32)
    n_boundaries = n_regions - 1
    quantum = 1.0 / _FRAG_LEN

    def region_zeros(dtype):
        return np.zeros((n_regions, 2), dtype=dtype)

    def boundary_zeros(dtype):
        return np.zeros((n_boundaries, 2), dtype=dtype)

    def flat(rows, dtype):
        """A single-column bank — the length moments and the conserved mass carry no strand axis."""
        return np.zeros(rows, dtype=dtype)

    payload = AccumulatorPayload(
        region_bounds=np.arange(n_regions + 1, dtype=np.int64) * 100,
        ref_region_bound_offsets=np.array([0, n_regions + 1], dtype=np.int64),
        ref_region_offsets=np.array([0, n_regions], dtype=np.int64),
        ref_boundary_offsets=np.array([0, n_boundaries], dtype=np.int64),
        ref_sj_offsets=np.array([0, 0], dtype=np.int64),
        region_contained_count=contained,
        # ⚠ ONE column: the length moments carry no strand axis, so the two are summed.
        region_contained_inv_opportunity_sum=(
            contained.sum(axis=1).astype(np.uint64) * np.uint64(quantum)
        ),
        region_start_count=contained.astype(np.uint32),
        region_end_count=contained.astype(np.uint32),
        region_span_count=region_zeros(np.uint32),
        boundary_unspliced_count=boundary_zeros(np.uint32),
        boundary_unspliced_inv_length_sum=flat(n_boundaries, np.float64),
        # ⚠ ONE value per boundary — the conserved mass has no strand axis. Zero here is a real state and
        # not a stub: this fixture deposits no crossings at all, so there is no mass to conserve.
        boundary_unspliced_mass=np.zeros(n_boundaries, dtype=np.float64),
        boundary_spliced_count=boundary_zeros(np.uint32),
        boundary_spliced_mass=np.zeros(n_boundaries, dtype=np.float64),
        sj_count=np.zeros((0, 2), dtype=np.uint32),
        sj_inv_length_sum=np.zeros(0, dtype=np.float64),
        sj_mass=np.zeros(0, dtype=np.float64),
        pool_lengths=np.zeros((N_FRAGMENT_POOLS, 201), dtype=np.int64),
        deposited_lengths=np.zeros(201, dtype=np.uint32),
        # ⚠ Nothing was deferred here, and that is a real state, not a stub: this fixture has no
        # annotated intron in any mate gap. The two empty spellings live on the classes so a
        # hand-built payload cannot get the `[0]`-not-`[]` offset boundary wrong.
        deferred=DeferredFragments.empty(),
        gap_resolution=GapCensus.zeros(),
        qc=ScanQC(
            deposited=int(contained.sum()),
            dropped_too_long=0,
            dropped_empty=0,
            dropped_strand_undefined=0,
            deferred_undetermined_gap=0,
            unannotated_introns=0,
            contradictory_sj_strand=0,
            introns_absorbed=0,
        ),
        max_length=200,
        n_refs=1,
    )
    starts = np.arange(n_regions, dtype=np.int64) * 100
    region_df = pd.DataFrame(
        {
            "region_id": np.arange(n_regions, dtype=np.int64),
            "ref_name": pd.array(["chr1"] * n_regions, dtype="string"),
            "start": starts,
            "end": starts + 100,
            "length": np.full(n_regions, 100, dtype=np.int64),
            "signature": np.full(
                n_regions, BIT_INTRON_POS, dtype=np.uint8
            ),  # intron ⇒ count-observable, genic
        }
    )
    return payload, RegionArrays.from_frame(region_df, {"chr1": 0})


def _calibrate(payload, ra, strand_model=None):
    return calibrate(
        payload=payload,
        region_arrays=ra,
        strand_model=strand_model if strand_model is not None else _STRAND_MODEL,
        gdna_fl_pmf=make_gdna_fl_pmf(mean=_FRAG_LEN),
        rna_fl_pmf=make_gdna_fl_pmf(mean=_FRAG_LEN),
        config=CalibrationConfig(),
    )


@pytest.mark.parametrize("od_true", [0.05, 0.10, 0.20])
def test_calibrate_recovers_overdispersion(od_true):
    payload, ra = _intron_betabinom_payload(
        n_regions=400, depth=150, overdispersion=od_true, seed=int(od_true * 1000)
    )
    # ⭐ ONE library, ONE technical process: the spliced seeds carry the same dispersion as the genomic
    # ones, so the two components agree and the reconciliation is a no-op — which is what makes this an
    # end-to-end RECOVERY test rather than a test of how far the components get pulled toward each other.
    result = _calibrate(payload, ra, _consistent_strand_model(od_true))
    assert result.gdna_strand_overdispersion == pytest.approx(od_true, rel=0.25, abs=0.02)
    assert result.rna_strand_overdispersion == pytest.approx(od_true, rel=0.30, abs=0.02)


def test_calibrate_binomial_gdna_floors_to_zero():
    """Non-overdispersed (50/50) gDNA → the identifiability gate floors od to 0 (Binomial)."""
    payload, ra = _intron_betabinom_payload(n_regions=400, depth=150, overdispersion=0.0, seed=7)
    assert _calibrate(payload, ra).gdna_strand_overdispersion < 0.02


def test_a_sj_free_library_calibrates(od_true=0.10):
    """⚠ ``n_sj == 0`` is legal and must not be confused with "no sj flux": this payload's
    references are all single-region-signature intergenic, so the graph has no sj boundary at all.
    ``calibrate`` defaults to an empty sj axis and the result carries a length-0 array."""
    payload, ra = _intron_betabinom_payload(
        n_regions=400, depth=150, overdispersion=od_true, seed=11
    )
    result = _calibrate(payload, ra)
    assert result.n_sj == 0
    assert result.count_rna_sj.shape == (0,)
    assert result.n_boundaries == result.n_regions - 1
