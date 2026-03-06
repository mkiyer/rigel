"""
test_batch_parity.py — Stress tests: C++ batch_locus_em vs Python fallback.

These tests run BOTH the C++ batch path and the Python per-locus loop on
identical inputs and assert exact numerical parity.  This catches:

  - Stride / column-count mismatches (the 6-vs-8 bug)
  - Off-by-one errors in array flattening
  - gDNA boundary handling (last component in the buffer)
  - Multi-locus accumulation ordering
  - Edge cases: empty loci, single-transcript loci, large loci

Every test builds a scenario, aligns, then runs the pipeline twice:
  (1) C++ batch path (default)
  (2) Python fallback (HULK_FORCE_PYTHON=1)
and asserts exact float equality on all output arrays.
"""

import os
import shutil

import numpy as np
import pytest

from hulkrna.config import EMConfig, PipelineConfig, BamScanConfig
from hulkrna.pipeline import run_pipeline
from hulkrna.sim import GDNAConfig, Scenario, SimConfig

pytestmark = pytest.mark.skipif(
    shutil.which("minimap2") is None or shutil.which("samtools") is None,
    reason="minimap2 and/or samtools not found in PATH",
)

# Fixed seed and simulation parameters for reproducibility
_SEED = 42
_N_FRAGS = 500
_SIM = SimConfig(
    frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
    read_length=100, strand_specificity=1.0, seed=_SEED,
)
_SIM_SS90 = SimConfig(
    frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
    read_length=100, strand_specificity=0.9, seed=_SEED,
)


# =====================================================================
# Helper: run pipeline via C++ batch AND Python, assert parity
# =====================================================================


def _run_and_compare(result, *, em_iterations=1000, include_multimap=False,
                     atol=0.0):
    """Run pipeline twice (C++ batch vs Python fallback), assert parity.

    Parameters
    ----------
    result : ScenarioResult
        Output of ``scenario.build()``.
    em_iterations : int
        EM iteration budget.
    include_multimap : bool
        Whether to include multimapped reads.
    atol : float
        Absolute tolerance for count comparison.  Use 0.0 for exact parity.
    """
    em_cfg = EMConfig(iterations=em_iterations, seed=_SEED)
    scan_cfg = BamScanConfig(
        sj_strand_tag="ts",
        include_multimap=include_multimap,
    )
    config = PipelineConfig(em=em_cfg, scan=scan_cfg)
    index = result.index

    # --- C++ batch path (default) ---
    saved_env = os.environ.get("HULK_FORCE_PYTHON")
    os.environ.pop("HULK_FORCE_PYTHON", None)
    try:
        pr_cpp = run_pipeline(result.bam_path, index, config=config)
    finally:
        if saved_env is not None:
            os.environ["HULK_FORCE_PYTHON"] = saved_env

    # --- Python fallback path ---
    os.environ["HULK_FORCE_PYTHON"] = "1"
    try:
        pr_py = run_pipeline(result.bam_path, index, config=config)
    finally:
        os.environ.pop("HULK_FORCE_PYTHON", None)
        if saved_env is not None:
            os.environ["HULK_FORCE_PYTHON"] = saved_env

    # --- Assert parity ---
    _assert_parity(pr_cpp, pr_py, index=index, atol=atol)
    return pr_cpp, pr_py


def _assert_parity(pr_cpp, pr_py, *, index, atol=0.0):
    """Assert that C++ and Python paths produce identical results."""
    df_cpp = pr_cpp.estimator.get_counts_df(index)
    df_py = pr_py.estimator.get_counts_df(index)

    # Sort by transcript_id for stable comparison
    df_cpp = df_cpp.sort_values("transcript_id").reset_index(drop=True)
    df_py = df_py.sort_values("transcript_id").reset_index(drop=True)

    # Core count columns that must match exactly
    count_cols = [
        "mrna", "mrna_unambig", "mrna_em",
        "mrna_high_conf", "mrna_spliced",
        "nrna", "rna_total",
    ]
    for col in count_cols:
        np.testing.assert_allclose(
            df_cpp[col].values, df_py[col].values,
            atol=atol,
            err_msg=f"Parity fail: transcript-level '{col}'",
        )

    # TPM and posterior_mean should also match
    np.testing.assert_allclose(
        df_cpp["tpm"].values, df_py["tpm"].values,
        atol=atol,
        err_msg="Parity fail: transcript-level 'tpm'",
    )
    np.testing.assert_allclose(
        df_cpp["posterior_mean"].values, df_py["posterior_mean"].values,
        atol=atol,
        err_msg="Parity fail: transcript-level 'posterior_mean'",
    )

    # gDNA EM totals
    np.testing.assert_allclose(
        pr_cpp.estimator.gdna_em_count,
        pr_py.estimator.gdna_em_count,
        atol=atol,
        err_msg="Parity fail: gDNA EM count",
    )

    # nRNA EM totals
    np.testing.assert_allclose(
        pr_cpp.estimator.nrna_em_count,
        pr_py.estimator.nrna_em_count,
        atol=atol,
        err_msg="Parity fail: nRNA EM count",
    )

    # Gene-level counts
    gdf_cpp = pr_cpp.estimator.get_gene_counts_df(index)
    gdf_py = pr_py.estimator.get_gene_counts_df(index)
    gdf_cpp = gdf_cpp.sort_values("gene_id").reset_index(drop=True)
    gdf_py = gdf_py.sort_values("gene_id").reset_index(drop=True)
    for col in ["mrna", "mrna_unambig", "mrna_em", "mrna_high_conf",
                 "mrna_spliced", "nrna", "rna_total"]:
        np.testing.assert_allclose(
            gdf_cpp[col].values, gdf_py[col].values,
            atol=atol,
            err_msg=f"Parity fail: gene-level '{col}'",
        )


# =====================================================================
# Test classes
# =====================================================================


class TestSingleGene:
    """Simplest case: one gene, no ambiguity."""

    def test_single_exon(self, tmp_path):
        sc = Scenario("s1e", genome_length=5000, seed=_SEED,
                      work_dir=tmp_path / "s1e")
        sc.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(500, 1000)], "abundance": 50},
        ])
        sc.add_gene("g_h", "-", [
            {"t_id": "t_h", "exons": [(3000, 3300)], "abundance": 20},
        ])
        result = sc.build(n_fragments=_N_FRAGS, sim_config=_SIM)
        _run_and_compare(result)

    def test_multi_exon(self, tmp_path):
        sc = Scenario("s3e", genome_length=5000, seed=_SEED,
                      work_dir=tmp_path / "s3e")
        sc.add_gene("g1", "+", [
            {"t_id": "t1",
             "exons": [(200, 400), (700, 900), (1200, 1400)],
             "abundance": 50},
        ])
        sc.add_gene("g_h", "-", [
            {"t_id": "t_h", "exons": [(3000, 3300)], "abundance": 20},
        ])
        result = sc.build(n_fragments=_N_FRAGS, sim_config=_SIM)
        _run_and_compare(result)


class TestTwoGenes:
    """Two separate genes — multiple loci but no cross-locus ambiguity."""

    def test_two_genes_separate(self, tmp_path):
        sc = Scenario("2g", genome_length=10000, seed=_SEED,
                      work_dir=tmp_path / "2g")
        sc.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(500, 800), (1200, 1500)],
             "abundance": 50},
        ])
        sc.add_gene("g2", "-", [
            {"t_id": "t2", "exons": [(4000, 4300), (4700, 5000)],
             "abundance": 30},
        ])
        result = sc.build(n_fragments=_N_FRAGS, sim_config=_SIM)
        _run_and_compare(result)

    def test_asymmetric_abundance(self, tmp_path):
        sc = Scenario("asym", genome_length=10000, seed=_SEED,
                      work_dir=tmp_path / "asym")
        sc.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(500, 1000)], "abundance": 100},
        ])
        sc.add_gene("g2", "-", [
            {"t_id": "t2", "exons": [(4000, 4500)], "abundance": 5},
        ])
        result = sc.build(n_fragments=_N_FRAGS, sim_config=_SIM)
        _run_and_compare(result)


class TestMultiIsoform:
    """Single gene with multiple isoforms — stresses component indexing."""

    def test_two_isoforms(self, tmp_path):
        sc = Scenario("2iso", genome_length=8000, seed=_SEED,
                      work_dir=tmp_path / "2iso")
        sc.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(500, 800), (1200, 1500)],
             "abundance": 40},
            {"t_id": "t2", "exons": [(500, 800), (1800, 2100)],
             "abundance": 40},
        ])
        sc.add_gene("g_h", "-", [
            {"t_id": "t_h", "exons": [(5000, 5300)], "abundance": 20},
        ])
        result = sc.build(n_fragments=_N_FRAGS, sim_config=_SIM)
        _run_and_compare(result)

    def test_three_isoforms_varying_abundance(self, tmp_path):
        sc = Scenario("3iso", genome_length=10000, seed=_SEED,
                      work_dir=tmp_path / "3iso")
        sc.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(500, 800), (1200, 1400), (2000, 2300)],
             "abundance": 60},
            {"t_id": "t2", "exons": [(500, 800), (1600, 1800), (2000, 2300)],
             "abundance": 30},
            {"t_id": "t3", "exons": [(500, 800), (2000, 2300)],
             "abundance": 10},
        ])
        sc.add_gene("g_h", "-", [
            {"t_id": "t_h", "exons": [(6000, 6300)], "abundance": 20},
        ])
        result = sc.build(n_fragments=_N_FRAGS, sim_config=_SIM)
        _run_and_compare(result)

    def test_five_isoforms(self, tmp_path):
        """Five isoforms — pushes component count."""
        sc = Scenario("5iso", genome_length=12000, seed=_SEED,
                      work_dir=tmp_path / "5iso")
        isoforms = []
        for i in range(5):
            mid = 1200 + i * 400
            isoforms.append({
                "t_id": f"t{i+1}",
                "exons": [(500, 800), (mid, mid + 200), (3500, 3800)],
                "abundance": 20 + i * 10,
            })
        sc.add_gene("g1", "+", isoforms)
        sc.add_gene("g_h", "-", [
            {"t_id": "t_h", "exons": [(8000, 8300)], "abundance": 20},
        ])
        result = sc.build(n_fragments=_N_FRAGS, sim_config=_SIM)
        _run_and_compare(result)


class TestParalogs:
    """Paralogs with multimapped reads — all reads are NH>1."""

    def test_identical_paralogs(self, tmp_path):
        sc = Scenario("par_id", genome_length=12000, seed=_SEED,
                      work_dir=tmp_path / "par_id")
        sc.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(500, 1000)], "abundance": 40},
        ])
        sc.add_gene("g2", "+", [
            {"t_id": "t2", "exons": [(5000, 5500)], "abundance": 40},
        ])
        sc.genome.edit(5000, sc.genome[500:1000])
        sc.add_gene("g_h", "+", [
            {"t_id": "t_h", "exons": [(8000, 8300), (8700, 9000)],
             "abundance": 50},
        ])
        result = sc.build(n_fragments=_N_FRAGS, sim_config=_SIM)
        _run_and_compare(result, include_multimap=True)

    def test_distinguishable_paralogs(self, tmp_path):
        sc = Scenario("par_dist", genome_length=12000, seed=_SEED,
                      work_dir=tmp_path / "par_dist")
        sc.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(500, 800), (1200, 1500)],
             "abundance": 40},
        ])
        sc.add_gene("g2", "+", [
            {"t_id": "t2", "exons": [(5000, 5300), (5700, 6000)],
             "abundance": 40},
        ])
        sc.genome.edit(5000, sc.genome[500:1500])
        sc.add_gene("g_h", "+", [
            {"t_id": "t_h", "exons": [(8000, 8300), (8700, 9000)],
             "abundance": 50},
        ])
        result = sc.build(n_fragments=_N_FRAGS, sim_config=_SIM)
        _run_and_compare(result, include_multimap=True)


class TestGDNA:
    """gDNA contamination — attacks the gDNA component (last in buffer)."""

    def test_light_gdna(self, tmp_path):
        gdna = GDNAConfig(abundance=5.0)
        sc = Scenario("gdna_lo", genome_length=10000, seed=_SEED,
                      work_dir=tmp_path / "gdna_lo")
        sc.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(500, 800), (1200, 1500)],
             "abundance": 50},
        ])
        sc.add_gene("g2", "-", [
            {"t_id": "t2", "exons": [(4000, 4300), (4700, 5000)],
             "abundance": 30},
        ])
        result = sc.build(n_fragments=_N_FRAGS, sim_config=_SIM,
                          gdna_config=gdna)
        _run_and_compare(result)

    def test_heavy_gdna(self, tmp_path):
        gdna = GDNAConfig(abundance=50.0)
        sc = Scenario("gdna_hi", genome_length=10000, seed=_SEED,
                      work_dir=tmp_path / "gdna_hi")
        sc.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(500, 800), (1200, 1500)],
             "abundance": 30},
        ])
        sc.add_gene("g2", "-", [
            {"t_id": "t2", "exons": [(4000, 4300), (4700, 5000)],
             "abundance": 30},
        ])
        result = sc.build(n_fragments=_N_FRAGS, sim_config=_SIM,
                          gdna_config=gdna)
        _run_and_compare(result)


class TestStrandSpecificity:
    """Imperfect strand specificity — tests strand column routing."""

    def test_ss_090(self, tmp_path):
        sc = Scenario("ss90", genome_length=8000, seed=_SEED,
                      work_dir=tmp_path / "ss90")
        sc.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(500, 800), (1200, 1500)],
             "abundance": 40},
        ])
        sc.add_gene("g2", "-", [
            {"t_id": "t2", "exons": [(3000, 3300), (3700, 4000)],
             "abundance": 30},
        ])
        result = sc.build(n_fragments=_N_FRAGS, sim_config=_SIM_SS90)
        _run_and_compare(result)


class TestEdgeCases:
    """Edge cases and boundary conditions."""

    def test_zero_abundance_gene(self, tmp_path):
        """Gene with 0 abundance — no reads, zero counts everywhere."""
        sc = Scenario("zero", genome_length=8000, seed=_SEED,
                      work_dir=tmp_path / "zero")
        sc.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(500, 800), (1200, 1500)],
             "abundance": 50},
        ])
        sc.add_gene("g_zero", "-", [
            {"t_id": "t_zero", "exons": [(3000, 3300), (3700, 4000)],
             "abundance": 0},
        ])
        result = sc.build(n_fragments=_N_FRAGS, sim_config=_SIM)
        _run_and_compare(result)

    def test_low_em_iterations(self, tmp_path):
        """Very few EM iterations — early-stop parity."""
        sc = Scenario("lowiter", genome_length=8000, seed=_SEED,
                      work_dir=tmp_path / "lowiter")
        sc.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(500, 800), (1200, 1500)],
             "abundance": 40},
            {"t_id": "t2", "exons": [(500, 800), (1800, 2100)],
             "abundance": 40},
        ])
        sc.add_gene("g_h", "-", [
            {"t_id": "t_h", "exons": [(5000, 5300)], "abundance": 20},
        ])
        result = sc.build(n_fragments=_N_FRAGS, sim_config=_SIM)
        _run_and_compare(result, em_iterations=5)


class TestManyLoci:
    """Many independent loci — stresses the batch loop."""

    def test_eight_genes(self, tmp_path):
        sc = Scenario("8g", genome_length=50000, seed=_SEED,
                      work_dir=tmp_path / "8g")
        for i in range(8):
            base = 1000 + i * 5000
            strand = "+" if i % 2 == 0 else "-"
            if i % 3 == 0:
                exons = [(base, base + 500)]
            elif i % 3 == 1:
                exons = [(base, base + 300), (base + 700, base + 1000)]
            else:
                exons = [(base, base + 200), (base + 500, base + 700),
                         (base + 1000, base + 1200)]
            sc.add_gene(f"g{i+1}", strand, [
                {"t_id": f"t{i+1}", "exons": exons,
                 "abundance": 20 + i * 5},
            ])
        result = sc.build(n_fragments=_N_FRAGS, sim_config=_SIM)
        _run_and_compare(result)


class TestCombined:
    """Scenarios combining multiple challenging features at once."""

    def test_paralogs_with_gdna(self, tmp_path):
        """Identical paralogs + gDNA contamination."""
        gdna = GDNAConfig(abundance=20.0)
        sc = Scenario("par_gdna", genome_length=12000, seed=_SEED,
                      work_dir=tmp_path / "par_gdna")
        sc.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(500, 1000)], "abundance": 40},
        ])
        sc.add_gene("g2", "+", [
            {"t_id": "t2", "exons": [(5000, 5500)], "abundance": 40},
        ])
        sc.genome.edit(5000, sc.genome[500:1000])
        sc.add_gene("g_h", "+", [
            {"t_id": "t_h", "exons": [(8000, 8300), (8700, 9000)],
             "abundance": 50},
        ])
        result = sc.build(n_fragments=_N_FRAGS, sim_config=_SIM,
                          gdna_config=gdna)
        _run_and_compare(result, include_multimap=True)

    def test_isoforms_gdna_imperfect_strand(self, tmp_path):
        """Multi-isoform + gDNA + imperfect strand specificity."""
        gdna = GDNAConfig(abundance=10.0)
        sc = Scenario("combo", genome_length=10000, seed=_SEED,
                      work_dir=tmp_path / "combo")
        sc.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(500, 800), (1200, 1500)],
             "abundance": 40},
            {"t_id": "t2", "exons": [(500, 800), (1800, 2100)],
             "abundance": 20},
        ])
        sc.add_gene("g2", "-", [
            {"t_id": "t3", "exons": [(4000, 4300), (4700, 5000)],
             "abundance": 30},
        ])
        result = sc.build(n_fragments=_N_FRAGS, sim_config=_SIM_SS90,
                          gdna_config=gdna)
        _run_and_compare(result)

    def test_nrna_scenario(self, tmp_path):
        """Nascent RNA scenario — tests nRNA component handling."""
        sc = Scenario("nrna", genome_length=10000, seed=_SEED,
                      work_dir=tmp_path / "nrna")
        sc.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(500, 800), (1200, 1500)],
             "abundance": 40},
        ])
        sc.add_gene("g2", "-", [
            {"t_id": "t2", "exons": [(4000, 4300), (4700, 5000)],
             "abundance": 30},
        ])
        result = sc.build(n_fragments=_N_FRAGS, sim_config=_SIM_SS90,
                          nrna_abundance=20)
        _run_and_compare(result)
