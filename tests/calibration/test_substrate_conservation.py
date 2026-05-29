"""Substrate mass conservation on a real scanned payload."""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.region_arrays import RegionArrays, boundary_region_indices
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.config import BamScanConfig, EMConfig, PipelineConfig
from rigel.pipeline import scan_and_buffer
from rigel.sim import ReadSimConfig, Scenario

SEED = 4321


@pytest.fixture(scope="module")
def scanned():
    import tempfile
    from pathlib import Path

    work = Path(tempfile.mkdtemp())
    sc = Scenario("subcons", genome_length=6000, seed=SEED, work_dir=work / "subcons")
    sc.add_gene("g1", "+", [{"t_id": "t1", "exons": [(300, 600), (900, 1200)], "abundance": 60}])
    sc.add_gene("g2", "-", [{"t_id": "t2", "exons": [(3000, 3300), (3700, 4000)], "abundance": 40}])
    result = sc.build_oracle(
        n_fragments=300, sim_config=ReadSimConfig(frag_mean=200, frag_std=30, seed=SEED)
    )
    config = PipelineConfig(em=EMConfig(seed=SEED), scan=BamScanConfig(sj_strand_tag="auto"))
    _, _, _, _, payload = scan_and_buffer(str(result.bam_path), result.index, config.scan)
    ra = RegionArrays.from_region_df(result.index.region_df, result.index.ref_name_to_id)
    yield payload, ra
    sc.cleanup()


def test_contained_mass_equals_count(scanned):
    payload, ra = scanned
    sub = CalibrationSubstrate.from_payload(payload, ra)
    # Per region, contained RNA+spliced mass equals the raw 4-channel total.
    raw_total = payload.region_contained.sum(axis=1).astype(np.float64)
    got = sub.contained.mass_unspliced + sub.contained.mass_spliced
    np.testing.assert_allclose(got, raw_total)


def test_boundary_attribution_no_double_count_no_loss(scanned):
    payload, ra = scanned
    sub = CalibrationSubstrate.from_payload(payload, ra)

    # Independent attribution via the inverse (boundary→region) map: mass_left[b]
    # belongs to b's left region, mass_right[b] to its right region; terminals (-1)
    # are off-edge and unattributed.
    lr, rr = boundary_region_indices(payload.ref_region_offsets, payload.ref_boundary_offsets)
    ml = payload.boundary_mass_left.astype(np.float64)
    mr = payload.boundary_mass_right.astype(np.float64)
    attributed = ml[lr >= 0].sum() + mr[rr >= 0].sum()

    sub_attributed = float(
        sub.left.mass_unspliced.sum()
        + sub.left.mass_spliced.sum()
        + sub.right.mass_unspliced.sum()
        + sub.right.mass_spliced.sum()
    )
    # The substrate (region→boundary map) attributes exactly the non-terminal mass.
    np.testing.assert_allclose(sub_attributed, attributed, rtol=1e-6, atol=1e-6)

    # No overcounting: attributed mass cannot exceed the total boundary mass.
    total = float(ml.sum() + mr.sum())
    assert sub_attributed <= total + 1e-6
