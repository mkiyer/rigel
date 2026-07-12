"""CI guard for the calibration oracle (scripts/debug/oracle.py).

The oracle's correctness rests on ONE identity: partitioning the sim BAM by true fragment origin and running
the production accumulator on each partition reproduces the full payload (the accumulator deposits each
fragment independently, so the per-origin parts must sum to the whole). This test builds a tiny gDNA + mRNA
+ nascent scenario and asserts that identity holds — if the accumulator ever changes in a way that breaks
per-fragment linearity, this fails loudly rather than letting a silently-wrong oracle percolate.
"""

import sys
from pathlib import Path

import numpy as np
import pytest

from rigel.config import PipelineConfig
from rigel.sim import Scenario, ReadSimConfig, GDNAConfig

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts" / "debug"))
from oracle import ORIGINS, OracleTruth  # noqa: E402


@pytest.fixture(scope="module")
def oracle_scenario(tmp_path_factory):
    wd = tmp_path_factory.mktemp("oracle_scn")
    sc = Scenario("orc", genome_length=6000, seed=11, work_dir=wd / "sim")
    sc.add_gene("g1", "+", [{"t_id": "t1", "exons": [(400, 700), (1200, 1500)], "abundance": 60}])
    sc.add_gene("g2", "-", [{"t_id": "t2", "exons": [(3000, 3300), (3800, 4100)], "abundance": 40}])
    result = sc.build_oracle(
        n_rna_fragments=1500,
        gdna_fraction=1.0,
        nrna_abundance=20.0,
        sim_config=ReadSimConfig(
            frag_mean=180,
            frag_std=30,
            frag_min=80,
            frag_max=400,
            read_length=90,
            strand_specificity=0.99,
            seed=11,
        ),
        gdna_config=GDNAConfig(abundance=0.0, frag_mean=200, frag_std=40),
    )
    return result


def test_oracle_validates_and_partitions_sum_to_full(oracle_scenario, tmp_path):
    # from_bam runs the sum-to-full + channel-sanity + fragment-accounting gates internally (raises on any).
    orc = OracleTruth.from_bam(
        str(oracle_scenario.bam_path), oracle_scenario.index, PipelineConfig(), tmp_path, "orc"
    )

    # region_contained partitions sum to full EXACTLY (integer channels).
    rc_full = np.asarray(orc.full.region_contained, np.int64)
    rc_sum = sum(np.asarray(orc.parts[k].region_contained, np.int64) for k in ORIGINS)
    assert np.array_equal(rc_sum, rc_full)

    # gDNA is never spliced: zero spliced (ch2,3) contained mass in the gdna partition.
    assert np.asarray(orc.parts["gdna"].region_contained, np.int64)[:, 2:].sum() == 0

    # every read accounted for; gDNA and mRNA partitions both non-empty (scenario has both).
    assert orc.read_counts["gdna"] > 0 and orc.read_counts["mrna"] > 0

    # the true gDNA fraction is a valid fraction wherever there is unspliced mass.
    fg, tot = orc.region_true_fg()
    assert np.all((fg[tot > 0] >= 0) & (fg[tot > 0] <= 1))


def test_oracle_override_conserves_node_mass(oracle_scenario, tmp_path):
    """The override masses (gdna + rna, spliced-inclusive) must equal the full node mass per region —
    the conservation the EM prior relies on."""
    from rigel.calibration.region_arrays import RegionArrays
    from rigel.calibration.substrate import CalibrationSubstrate

    orc = OracleTruth.from_bam(
        str(oracle_scenario.bam_path), oracle_scenario.index, PipelineConfig(), tmp_path, "orc2"
    )
    ra = RegionArrays.from_region_df(
        oracle_scenario.index.region_df, oracle_scenario.index.ref_name_to_id
    )
    ov = orc.override_masses(ra)
    full_sub = CalibrationSubstrate.from_payload(orc.full, ra)
    total = (
        np.asarray(full_sub.contained.mass_unspliced)
        + np.asarray(full_sub.contained.mass_spliced)
        + np.asarray(full_sub.left.mass_unspliced)
        + np.asarray(full_sub.left.mass_spliced)
        + np.asarray(full_sub.right.mass_unspliced)
        + np.asarray(full_sub.right.mass_spliced)
    )
    got = (
        ov["mass_gdna_contained"]
        + ov["mass_rna_contained"]
        + ov["mass_gdna_left"]
        + ov["mass_rna_left"]
        + ov["mass_gdna_right"]
        + ov["mass_rna_right"]
    )
    assert np.allclose(got, total, atol=1e-3)
