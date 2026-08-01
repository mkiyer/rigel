"""CI guard for the calibration oracle (``tests/calibration/_oracle.py``).

The oracle's correctness rests on ONE identity: partitioning the sim BAM by true fragment origin and running
the production accumulator on each partition reproduces the full payload (the accumulator deposits each
fragment independently, so the per-origin parts must sum to the whole). This test builds a tiny gDNA + mRNA
+ nascent scenario and asserts that identity holds — if the accumulator ever changes in a way that breaks
per-fragment linearity, this fails loudly rather than letting a silently-wrong oracle percolate.
"""

import numpy as np
import pytest

from rigel.config import PipelineConfig
from rigel.sim import Scenario, ReadSimConfig, GDNAConfig

from _oracle import _BANKS, ORIGINS, OracleTruth  # noqa: E402


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

    # ⭐ EVERY bank on all three axes sums to full EXACTLY — no tolerance anywhere, because every bank
    # is an integer count. The predecessor could only be exact on two of its four arrays.
    for bank in _BANKS:
        full = np.asarray(getattr(orc.full, bank), np.int64)
        parts = sum(np.asarray(getattr(orc.parts[k], bank), np.int64) for k in ORIGINS)
        np.testing.assert_array_equal(parts, full, err_msg=f"{bank} does not sum to full")

    # gDNA is never spliced — on the contiguous-edge spliced bank AND on the junction axis.
    assert np.asarray(orc.parts["gdna"].edge_spliced_count, np.int64).sum() == 0
    assert np.asarray(orc.parts["gdna"].sj_count, np.int64).sum() == 0

    # ⚠ The scenario must actually EXERCISE the RNA-only banks, or "gDNA is zero there" is vacuous.
    assert np.asarray(orc.full.sj_count, np.int64).sum() > 0

    # every read accounted for; gDNA and mRNA partitions both non-empty (scenario has both).
    assert orc.read_counts["gdna"] > 0 and orc.read_counts["mrna"] > 0

    # the true gDNA fraction is a valid fraction wherever there is contained mass.
    fg, tot = orc.node_true_fg()
    assert np.all((fg[tot > 0] >= 0) & (fg[tot > 0] <= 1))


def test_oracle_override_conserves_mass_on_EACH_AXIS_SEPARATELY(oracle_scenario, tmp_path):
    """The override masses must equal the full object count — checked **per axis**, not pooled.

    ⚠ Pooling the two axes into one total would let an error on the node axis cancel an equal and
    opposite one on the edge axis, which is exactly the class of mistake a three-axis schema makes
    possible. ``E`` and ``N`` differ by only ``n_refs``, so such a cancellation is not far-fetched.
    """
    from rigel.calibration.region_arrays import RegionArrays
    from rigel.calibration.substrate import CalibrationSubstrate

    orc = OracleTruth.from_bam(
        str(oracle_scenario.bam_path), oracle_scenario.index, PipelineConfig(), tmp_path, "orc2"
    )
    ra = RegionArrays.from_frame(
        oracle_scenario.index.nodes_df, oracle_scenario.index.ref_name_to_id
    )
    ov = orc.override_masses(ra)
    full = CalibrationSubstrate.from_payload(orc.full, ra)

    # NODE axis: a node's contained population holds no spliced molecule, so its total is one bank.
    np.testing.assert_allclose(
        ov["mass_gdna_node"] + ov["mass_rna_node"],
        np.asarray(full.node_contained.count, np.float64).sum(1),
    )
    # EDGE axis: unspliced + spliced, because mass_rna_edge is spliced-inclusive.
    np.testing.assert_allclose(
        ov["mass_gdna_edge"] + ov["mass_rna_edge"],
        np.asarray(full.edge_unspliced.count, np.float64).sum(1)
        + np.asarray(full.edge_spliced.count, np.float64).sum(1),
    )
    # JUNCTION axis: never deconvolved — the flux verbatim.
    np.testing.assert_allclose(
        ov["mass_rna_junction"], np.asarray(full.junction.count, np.float64).sum(1)
    )


def test_the_oracle_result_is_a_VALID_CalibrationResult(oracle_scenario, tmp_path):
    """⭐ The perfect-calibration lever must actually construct. ``override_masses`` returns exactly the
    fields ``dataclasses.replace`` needs, so a rename in the schema that it missed would surface here
    rather than in whatever A/B first tried to use it."""
    import dataclasses

    from rigel.calibration.region_arrays import RegionArrays
    from rigel.calibration.result import CalibrationResult
    from rigel.config import CalibrationConfig

    orc = OracleTruth.from_bam(
        str(oracle_scenario.bam_path), oracle_scenario.index, PipelineConfig(), tmp_path, "orc3"
    )
    ra = RegionArrays.from_frame(
        oracle_scenario.index.nodes_df, oracle_scenario.index.ref_name_to_id
    )
    ov = orc.override_masses(ra)
    n, e, j = orc.full.n_nodes, orc.full.n_edges, orc.full.n_sj
    blank = CalibrationResult(
        mass_gdna_node=np.zeros(n),
        mass_rna_node=np.zeros(n),
        mass_gdna_edge=np.zeros(e),
        mass_rna_edge=np.zeros(e),
        mass_rna_spliced_edge=np.zeros(e),
        mass_rna_junction=np.zeros(j),
        gdna_node_eff_len=np.ones(n),
        gdna_edge_eff_len=np.ones(e),
        rna_node_eff_len=np.ones(n),
        rna_edge_eff_len=np.ones(e),
        gdna_density_global=0.0,
        rna_sense_frac=0.5,
        gdna_strand_overdispersion=0.0,
        rna_strand_overdispersion=0.0,
        n_nodes=n,
        n_edges=e,
        n_junctions=j,
        config=CalibrationConfig(),
    )
    truth = dataclasses.replace(blank, **ov)  # __post_init__ re-validates every axis
    assert truth.mass_rna_junction.sum() > 0
