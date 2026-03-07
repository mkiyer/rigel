"""
test_golden_output.py — Bit-exact regression tests for pipeline output stability.

PURPOSE
-------
After C++ or type-narrowing changes, this test suite verifies that the
pipeline produces **identical** numerical output on a battery of diverse
oracle scenarios.  Each scenario captures the full transcript-level,
gene-level, and locus-level output DataFrames plus scalar aggregates,
serialized to a golden CSV file.  Any drift is flagged immediately.

USAGE
-----
Normal run (compare against golden files):
    pytest tests/test_golden_output.py -v

Regenerate golden files after an intentional change:
    pytest tests/test_golden_output.py -v --update-golden

Only run a specific scenario:
    pytest tests/test_golden_output.py -v -k "multi_isoform_3tx"

DESIGN
------
- Uses oracle BAM (no minimap2/samtools dependency)
- Fixed seeds for complete reproducibility
- Captures every numerical column from get_counts_df, get_gene_counts_df,
  get_loci_df, plus gdna_em_count, nrna_em_count scalars
- Comparison is bit-exact (atol=0) on float64 representation
- Golden files stored in tests/golden/ as feather files for lossless
  float round-trip; human-readable TSV mirrors alongside for diffing
"""

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from rigel.config import BamScanConfig, EMConfig, PipelineConfig
from rigel.pipeline import run_pipeline
from rigel.sim import GDNAConfig, Scenario, SimConfig

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

GOLDEN_DIR = Path(__file__).parent / "golden"
SEED = 42
PIPELINE_SEED = 42
N_FRAGS = 1000  # enough to exercise EM meaningfully

# Tolerance for golden comparison.
# The C++ EM solver has inherent ULP-level (~1e-16) non-determinism in
# floating-point accumulation order, and near-zero quantities (e.g. gdna
# when there is no contamination, ~1e-296) can wander by several percent
# in relative terms while remaining scientifically meaningless.
# rtol=1e-12 catches real regressions on meaningful values; atol=1e-10
# absorbs noise on effectively-zero quantities.
RTOL = 1e-12
ATOL = 1e-10

# Standard simulation configs
_SIM_SS100 = SimConfig(
    frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
    read_length=100, strand_specificity=1.0, seed=SEED,
)
_SIM_SS90 = SimConfig(
    frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
    read_length=100, strand_specificity=0.9, seed=SEED,
)
_SIM_SS65 = SimConfig(
    frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
    read_length=100, strand_specificity=0.65, seed=SEED,
)


def _pipeline_config(seed=PIPELINE_SEED):
    return PipelineConfig(
        em=EMConfig(seed=seed),
        scan=BamScanConfig(sj_strand_tag="auto"),
    )


# ---------------------------------------------------------------------------
# Fixture: --update-golden flag (registered in conftest.py)
# ---------------------------------------------------------------------------


@pytest.fixture
def update_golden(request):
    return request.config.getoption("--update-golden")


# ---------------------------------------------------------------------------
# Core: capture pipeline output as DataFrames
# ---------------------------------------------------------------------------


def _capture_output(result, index):
    """Run the pipeline and capture all output DataFrames + scalars.

    Returns a dict with keys: 'transcript_df', 'gene_df', 'loci_df', 'scalars'.
    """
    pr = run_pipeline(result.bam_path, index, config=_pipeline_config())
    est = pr.estimator

    transcript_df = est.get_counts_df(index)
    gene_df = est.get_gene_counts_df(index)
    loci_df = est.get_loci_df()

    # Sort for deterministic ordering
    transcript_df = transcript_df.sort_values("transcript_id").reset_index(drop=True)
    gene_df = gene_df.sort_values("gene_id").reset_index(drop=True)
    if len(loci_df) > 0:
        loci_df = loci_df.sort_values("locus_id").reset_index(drop=True)

    scalars = {
        "gdna_em_count": float(est.gdna_em_count),
        "nrna_em_count": float(est.nrna_em_count),
        "n_fragments": int(pr.stats.n_fragments),
    }

    return {
        "transcript_df": transcript_df,
        "gene_df": gene_df,
        "loci_df": loci_df,
        "scalars": scalars,
    }


# ---------------------------------------------------------------------------
# Golden file I/O
# ---------------------------------------------------------------------------


def _golden_path(scenario_name, suffix):
    return GOLDEN_DIR / f"{scenario_name}{suffix}"


def _save_golden(scenario_name, output):
    """Save golden output to feather + TSV files."""
    GOLDEN_DIR.mkdir(parents=True, exist_ok=True)

    for key in ("transcript_df", "gene_df", "loci_df"):
        df = output[key]
        feather_path = _golden_path(scenario_name, f"_{key}.feather")
        tsv_path = _golden_path(scenario_name, f"_{key}.tsv")
        df.to_feather(str(feather_path))
        df.to_csv(str(tsv_path), sep="\t", index=False)

    scalars_path = _golden_path(scenario_name, "_scalars.json")
    with open(scalars_path, "w") as f:
        json.dump(output["scalars"], f, indent=2, sort_keys=True)


def _load_golden(scenario_name):
    """Load golden output from feather files."""
    output = {}
    for key in ("transcript_df", "gene_df", "loci_df"):
        feather_path = _golden_path(scenario_name, f"_{key}.feather")
        output[key] = pd.read_feather(str(feather_path))

    scalars_path = _golden_path(scenario_name, "_scalars.json")
    with open(scalars_path) as f:
        output["scalars"] = json.load(f)

    return output


def _golden_exists(scenario_name):
    return _golden_path(scenario_name, "_scalars.json").exists()


# ---------------------------------------------------------------------------
# Comparison
# ---------------------------------------------------------------------------

# Columns to compare with bit-exact precision
_TRANSCRIPT_NUMERIC_COLS = [
    "effective_length",
    "mrna", "mrna_unambig", "mrna_em", "mrna_high_conf", "mrna_spliced",
    "nrna", "rna_total", "tpm", "posterior_mean",
]

_GENE_NUMERIC_COLS = [
    "effective_length",
    "mrna", "mrna_unambig", "mrna_em", "mrna_high_conf", "mrna_spliced",
    "nrna", "rna_total", "tpm",
]

_LOCI_NUMERIC_COLS = [
    "n_em_fragments", "mrna", "nrna", "gdna", "gdna_rate", "gdna_init",
]


def _compare_df(actual, expected, numeric_cols, label):
    """Compare two DataFrames column-by-column with informative errors."""
    assert len(actual) == len(expected), (
        f"{label}: row count mismatch: {len(actual)} vs {len(expected)}"
    )
    for col in numeric_cols:
        if col not in actual.columns:
            continue
        a = actual[col].values.astype(np.float64)
        e = expected[col].values.astype(np.float64)
        # Handle NaN equality
        nan_mask = np.isnan(e)
        if nan_mask.any():
            assert np.isnan(a[nan_mask]).all(), (
                f"{label}.{col}: NaN mismatch"
            )
            a = a[~nan_mask]
            e = e[~nan_mask]
        np.testing.assert_allclose(
            a, e, atol=ATOL, rtol=RTOL,
            err_msg=f"{label}: column '{col}' differs",
        )


def _compare_output(actual, expected, scenario_name):
    """Full comparison of pipeline output against golden data."""
    # Transcript-level
    _compare_df(
        actual["transcript_df"], expected["transcript_df"],
        _TRANSCRIPT_NUMERIC_COLS, f"{scenario_name}/transcript",
    )

    # Gene-level
    _compare_df(
        actual["gene_df"], expected["gene_df"],
        _GENE_NUMERIC_COLS, f"{scenario_name}/gene",
    )

    # Loci-level
    _compare_df(
        actual["loci_df"], expected["loci_df"],
        _LOCI_NUMERIC_COLS, f"{scenario_name}/loci",
    )

    # Scalars
    for key in ("gdna_em_count", "nrna_em_count", "n_fragments"):
        np.testing.assert_allclose(
            actual["scalars"][key], expected["scalars"][key],
            atol=ATOL, rtol=RTOL,
            err_msg=f"{scenario_name}: scalar '{key}' differs",
        )


# ---------------------------------------------------------------------------
# Core test runner
# ---------------------------------------------------------------------------


def _run_golden_test(scenario_name, result, index, update_golden):
    """Run pipeline, compare or save golden output."""
    output = _capture_output(result, index)

    if update_golden:
        _save_golden(scenario_name, output)
        pytest.skip(f"Golden files updated for '{scenario_name}'")

    if not _golden_exists(scenario_name):
        _save_golden(scenario_name, output)
        pytest.skip(
            f"Golden files created for '{scenario_name}' (first run). "
            f"Re-run to verify stability."
        )

    expected = _load_golden(scenario_name)
    _compare_output(output, expected, scenario_name)


# =====================================================================
# SCENARIO DEFINITIONS
# =====================================================================
#
# Each scenario exercises a different axis of complexity.  Together they
# cover: single/multi-exon genes, multi-isoform EM, antisense overlap,
# gDNA contamination, nRNA, imperfect strand specificity, many-loci
# batching, extreme abundance ratios, and combinations thereof.
# =====================================================================


class TestSingleGene:
    """Baseline: one gene, no ambiguity — pure deterministic routing."""

    def test_single_exon_clean(self, tmp_path, update_golden):
        sc = Scenario("golden_1e", genome_length=5000, seed=SEED,
                      work_dir=tmp_path / "golden_1e")
        sc.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(500, 1200)], "abundance": 100},
        ])
        sc.add_gene("g_h", "-", [
            {"t_id": "t_h", "exons": [(3000, 3500)], "abundance": 30},
        ])
        result = sc.build_oracle(n_fragments=N_FRAGS, sim_config=_SIM_SS100)
        _run_golden_test("single_exon_clean", result, result.index, update_golden)
        sc.cleanup()

    def test_multi_exon_spliced(self, tmp_path, update_golden):
        sc = Scenario("golden_3e", genome_length=8000, seed=SEED,
                      work_dir=tmp_path / "golden_3e")
        sc.add_gene("g1", "+", [
            {"t_id": "t1",
             "exons": [(200, 500), (1000, 1300), (2000, 2300)],
             "abundance": 80},
        ])
        sc.add_gene("g_h", "-", [
            {"t_id": "t_h", "exons": [(5000, 5400)], "abundance": 20},
        ])
        result = sc.build_oracle(n_fragments=N_FRAGS, sim_config=_SIM_SS100)
        _run_golden_test("multi_exon_spliced", result, result.index, update_golden)
        sc.cleanup()


class TestMultiIsoform:
    """Multi-isoform genes — stresses EM transcript disambiguation."""

    def test_two_isoforms_equal(self, tmp_path, update_golden):
        sc = Scenario("golden_2iso_eq", genome_length=8000, seed=SEED,
                      work_dir=tmp_path / "golden_2iso_eq")
        sc.add_gene("g1", "+", [
            {"t_id": "t1",
             "exons": [(500, 800), (1200, 1500)],
             "abundance": 50},
            {"t_id": "t2",
             "exons": [(500, 800), (1800, 2100)],
             "abundance": 50},
        ])
        sc.add_gene("g_h", "-", [
            {"t_id": "t_h", "exons": [(5000, 5300)], "abundance": 20},
        ])
        result = sc.build_oracle(n_fragments=N_FRAGS, sim_config=_SIM_SS100)
        _run_golden_test("multi_isoform_2tx_equal", result, result.index, update_golden)
        sc.cleanup()

    def test_three_isoforms_skewed(self, tmp_path, update_golden):
        """3 isoforms with 10:3:1 abundance ratio — tests EM convergence."""
        sc = Scenario("golden_3iso", genome_length=10000, seed=SEED,
                      work_dir=tmp_path / "golden_3iso")
        sc.add_gene("g1", "+", [
            {"t_id": "t1",
             "exons": [(500, 800), (1200, 1400), (2000, 2300)],
             "abundance": 100},
            {"t_id": "t2",
             "exons": [(500, 800), (1600, 1800), (2000, 2300)],
             "abundance": 30},
            {"t_id": "t3",
             "exons": [(500, 800), (2000, 2300)],
             "abundance": 10},
        ])
        sc.add_gene("g_h", "-", [
            {"t_id": "t_h", "exons": [(6000, 6300)], "abundance": 20},
        ])
        result = sc.build_oracle(n_fragments=N_FRAGS, sim_config=_SIM_SS100)
        _run_golden_test("multi_isoform_3tx_skewed", result, result.index, update_golden)
        sc.cleanup()

    def test_five_isoforms(self, tmp_path, update_golden):
        """5 isoforms — high component count per locus."""
        sc = Scenario("golden_5iso", genome_length=12000, seed=SEED,
                      work_dir=tmp_path / "golden_5iso")
        isoforms = []
        for i in range(5):
            mid = 1200 + i * 400
            isoforms.append({
                "t_id": f"t{i+1}",
                "exons": [(500, 800), (mid, mid + 200), (3500, 3800)],
                "abundance": 20 + i * 15,
            })
        sc.add_gene("g1", "+", isoforms)
        sc.add_gene("g_h", "-", [
            {"t_id": "t_h", "exons": [(8000, 8300)], "abundance": 20},
        ])
        result = sc.build_oracle(n_fragments=N_FRAGS, sim_config=_SIM_SS100)
        _run_golden_test("multi_isoform_5tx", result, result.index, update_golden)
        sc.cleanup()


class TestAntisense:
    """Overlapping antisense genes — tests strand separation."""

    def test_antisense_overlap(self, tmp_path, update_golden):
        sc = Scenario("golden_as", genome_length=8000, seed=SEED,
                      work_dir=tmp_path / "golden_as")
        sc.add_gene("g_sense", "+", [
            {"t_id": "t_sense",
             "exons": [(500, 800), (1200, 1500)],
             "abundance": 60},
        ])
        sc.add_gene("g_anti", "-", [
            {"t_id": "t_anti",
             "exons": [(400, 700), (1100, 1400)],
             "abundance": 40},
        ])
        result = sc.build_oracle(n_fragments=N_FRAGS, sim_config=_SIM_SS90)
        _run_golden_test("antisense_overlap_ss90", result, result.index, update_golden)
        sc.cleanup()

    def test_antisense_contained(self, tmp_path, update_golden):
        """Antisense gene fully contained within sense gene's exon."""
        sc = Scenario("golden_as_cont", genome_length=8000, seed=SEED,
                      work_dir=tmp_path / "golden_as_cont")
        sc.add_gene("g_sense", "+", [
            {"t_id": "t_sense", "exons": [(500, 2000)], "abundance": 80},
        ])
        sc.add_gene("g_anti", "-", [
            {"t_id": "t_anti", "exons": [(700, 1500)], "abundance": 30},
        ])
        result = sc.build_oracle(n_fragments=N_FRAGS, sim_config=_SIM_SS90)
        _run_golden_test("antisense_contained_ss90", result, result.index, update_golden)
        sc.cleanup()


class TestGDNA:
    """gDNA contamination at various levels."""

    def test_gdna_light(self, tmp_path, update_golden):
        gdna = GDNAConfig(abundance=10.0, frag_mean=350, frag_std=100,
                          frag_min=100, frag_max=1000)
        sc = Scenario("golden_gdna_lo", genome_length=10000, seed=SEED,
                      work_dir=tmp_path / "golden_gdna_lo")
        sc.add_gene("g1", "+", [
            {"t_id": "t1",
             "exons": [(500, 800), (1200, 1500)],
             "abundance": 60},
        ])
        sc.add_gene("g2", "-", [
            {"t_id": "t2",
             "exons": [(4000, 4300), (4700, 5000)],
             "abundance": 40},
        ])
        result = sc.build_oracle(n_fragments=N_FRAGS, sim_config=_SIM_SS100,
                                  gdna_config=gdna)
        _run_golden_test("gdna_light", result, result.index, update_golden)
        sc.cleanup()

    def test_gdna_heavy(self, tmp_path, update_golden):
        gdna = GDNAConfig(abundance=80.0, frag_mean=350, frag_std=100,
                          frag_min=100, frag_max=1000)
        sc = Scenario("golden_gdna_hi", genome_length=10000, seed=SEED,
                      work_dir=tmp_path / "golden_gdna_hi")
        sc.add_gene("g1", "+", [
            {"t_id": "t1",
             "exons": [(500, 800), (1200, 1500)],
             "abundance": 40},
        ])
        sc.add_gene("g2", "-", [
            {"t_id": "t2",
             "exons": [(4000, 4300), (4700, 5000)],
             "abundance": 30},
        ])
        result = sc.build_oracle(n_fragments=N_FRAGS, sim_config=_SIM_SS100,
                                  gdna_config=gdna)
        _run_golden_test("gdna_heavy", result, result.index, update_golden)
        sc.cleanup()


class TestNRNA:
    """Nascent RNA — tests nRNA component handling in EM."""

    def test_nrna_moderate(self, tmp_path, update_golden):
        sc = Scenario("golden_nrna", genome_length=10000, seed=SEED,
                      work_dir=tmp_path / "golden_nrna")
        sc.add_gene("g1", "+", [
            {"t_id": "t1",
             "exons": [(500, 800), (1200, 1500)],
             "abundance": 50},
        ])
        sc.add_gene("g2", "-", [
            {"t_id": "t2",
             "exons": [(4000, 4300), (4700, 5000)],
             "abundance": 40},
        ])
        result = sc.build_oracle(n_fragments=N_FRAGS, sim_config=_SIM_SS90,
                                  nrna_abundance=30)
        _run_golden_test("nrna_moderate_ss90", result, result.index, update_golden)
        sc.cleanup()

    def test_nrna_heavy(self, tmp_path, update_golden):
        sc = Scenario("golden_nrna_hi", genome_length=10000, seed=SEED,
                      work_dir=tmp_path / "golden_nrna_hi")
        sc.add_gene("g1", "+", [
            {"t_id": "t1",
             "exons": [(500, 800), (1200, 1500)],
             "abundance": 40},
        ])
        sc.add_gene("g2", "-", [
            {"t_id": "t2",
             "exons": [(4000, 4300), (4700, 5000)],
             "abundance": 30},
        ])
        result = sc.build_oracle(n_fragments=N_FRAGS, sim_config=_SIM_SS90,
                                  nrna_abundance=70)
        _run_golden_test("nrna_heavy_ss90", result, result.index, update_golden)
        sc.cleanup()


class TestStrandSpecificity:
    """Imperfect strand specificity — tests strand column routing."""

    def test_low_strand_specificity(self, tmp_path, update_golden):
        sc = Scenario("golden_ss65", genome_length=10000, seed=SEED,
                      work_dir=tmp_path / "golden_ss65")
        sc.add_gene("g1", "+", [
            {"t_id": "t1",
             "exons": [(500, 800), (1200, 1500)],
             "abundance": 50},
            {"t_id": "t2",
             "exons": [(500, 800), (1800, 2100)],
             "abundance": 30},
        ])
        sc.add_gene("g2", "-", [
            {"t_id": "t3",
             "exons": [(4000, 4300), (4700, 5000)],
             "abundance": 40},
        ])
        result = sc.build_oracle(n_fragments=N_FRAGS, sim_config=_SIM_SS65)
        _run_golden_test("strand_ss65_multi_iso", result, result.index, update_golden)
        sc.cleanup()


class TestManyLoci:
    """Many independent loci — stresses the C++ batch loop."""

    def test_ten_genes(self, tmp_path, update_golden):
        sc = Scenario("golden_10g", genome_length=60000, seed=SEED,
                      work_dir=tmp_path / "golden_10g")
        for i in range(10):
            base = 1000 + i * 5000
            strand = "+" if i % 2 == 0 else "-"
            if i % 3 == 0:
                exons = [(base, base + 600)]
            elif i % 3 == 1:
                exons = [(base, base + 300), (base + 700, base + 1000)]
            else:
                exons = [(base, base + 200), (base + 500, base + 700),
                         (base + 1100, base + 1300)]
            sc.add_gene(f"g{i+1}", strand, [
                {"t_id": f"t{i+1}", "exons": exons,
                 "abundance": 15 + i * 8},
            ])
        result = sc.build_oracle(n_fragments=N_FRAGS, sim_config=_SIM_SS100)
        _run_golden_test("many_loci_10genes", result, result.index, update_golden)
        sc.cleanup()

    def test_mixed_loci_with_isoforms(self, tmp_path, update_golden):
        """6 genes, 3 with multiple isoforms — mixed complexity."""
        sc = Scenario("golden_mixed", genome_length=50000, seed=SEED,
                      work_dir=tmp_path / "golden_mixed")
        # Gene 1: single exon
        sc.add_gene("g1", "+", [
            {"t_id": "t1", "exons": [(500, 1200)], "abundance": 40},
        ])
        # Gene 2: two isoforms
        sc.add_gene("g2", "-", [
            {"t_id": "t2a",
             "exons": [(6000, 6300), (6700, 7000)],
             "abundance": 30},
            {"t_id": "t2b",
             "exons": [(6000, 6300), (7400, 7700)],
             "abundance": 20},
        ])
        # Gene 3: three isoforms
        sc.add_gene("g3", "+", [
            {"t_id": "t3a",
             "exons": [(12000, 12300), (12700, 12900), (13400, 13700)],
             "abundance": 50},
            {"t_id": "t3b",
             "exons": [(12000, 12300), (13100, 13300), (13400, 13700)],
             "abundance": 25},
            {"t_id": "t3c",
             "exons": [(12000, 12300), (13400, 13700)],
             "abundance": 10},
        ])
        # Gene 4: single exon, opposite strand
        sc.add_gene("g4", "-", [
            {"t_id": "t4", "exons": [(20000, 20800)], "abundance": 60},
        ])
        # Gene 5: wide intron
        sc.add_gene("g5", "+", [
            {"t_id": "t5",
             "exons": [(26000, 26300), (30000, 30300)],
             "abundance": 35},
        ])
        # Gene 6: low abundance
        sc.add_gene("g6", "-", [
            {"t_id": "t6",
             "exons": [(38000, 38300), (38700, 39000)],
             "abundance": 5},
        ])
        result = sc.build_oracle(n_fragments=N_FRAGS, sim_config=_SIM_SS100)
        _run_golden_test("mixed_loci_isoforms", result, result.index, update_golden)
        sc.cleanup()


class TestCombinedStress:
    """Combined stress: multi-isoform + gDNA + nRNA + imperfect strand."""

    def test_full_combo_moderate(self, tmp_path, update_golden):
        gdna = GDNAConfig(abundance=20.0, frag_mean=350, frag_std=100,
                          frag_min=100, frag_max=1000)
        sc = Scenario("golden_combo_mod", genome_length=12000, seed=SEED,
                      work_dir=tmp_path / "golden_combo_mod")
        sc.add_gene("g1", "+", [
            {"t_id": "t1",
             "exons": [(500, 800), (1200, 1500)],
             "abundance": 50},
            {"t_id": "t2",
             "exons": [(500, 800), (1800, 2100)],
             "abundance": 25},
        ])
        sc.add_gene("g2", "-", [
            {"t_id": "t3",
             "exons": [(5000, 5300), (5700, 6000)],
             "abundance": 40},
        ])
        result = sc.build_oracle(n_fragments=N_FRAGS, sim_config=_SIM_SS90,
                                  gdna_config=gdna, nrna_abundance=20)
        _run_golden_test("combo_moderate", result, result.index, update_golden)
        sc.cleanup()

    def test_full_combo_extreme(self, tmp_path, update_golden):
        """Extreme: heavy gDNA + heavy nRNA + low strand spec + 3 isoforms."""
        gdna = GDNAConfig(abundance=80.0, frag_mean=350, frag_std=100,
                          frag_min=100, frag_max=1000)
        sc = Scenario("golden_combo_ext", genome_length=12000, seed=SEED,
                      work_dir=tmp_path / "golden_combo_ext")
        sc.add_gene("g1", "+", [
            {"t_id": "t1",
             "exons": [(500, 800), (1200, 1400), (2000, 2300)],
             "abundance": 60},
            {"t_id": "t2",
             "exons": [(500, 800), (1600, 1800), (2000, 2300)],
             "abundance": 20},
            {"t_id": "t3",
             "exons": [(500, 800), (2000, 2300)],
             "abundance": 10},
        ])
        sc.add_gene("g2", "-", [
            {"t_id": "t4",
             "exons": [(5000, 5300), (5700, 6000)],
             "abundance": 30},
        ])
        result = sc.build_oracle(n_fragments=N_FRAGS, sim_config=_SIM_SS65,
                                  gdna_config=gdna, nrna_abundance=60)
        _run_golden_test("combo_extreme", result, result.index, update_golden)
        sc.cleanup()


class TestEdgeCases:
    """Boundary conditions and edge cases."""

    def test_extreme_abundance_ratio(self, tmp_path, update_golden):
        """100:1 abundance ratio between isoforms."""
        sc = Scenario("golden_ratio", genome_length=8000, seed=SEED,
                      work_dir=tmp_path / "golden_ratio")
        sc.add_gene("g1", "+", [
            {"t_id": "t_major",
             "exons": [(500, 800), (1200, 1500)],
             "abundance": 200},
            {"t_id": "t_minor",
             "exons": [(500, 800), (1800, 2100)],
             "abundance": 2},
        ])
        sc.add_gene("g_h", "-", [
            {"t_id": "t_h", "exons": [(5000, 5300)], "abundance": 30},
        ])
        result = sc.build_oracle(n_fragments=N_FRAGS, sim_config=_SIM_SS100)
        _run_golden_test("extreme_abundance_ratio", result, result.index, update_golden)
        sc.cleanup()

    def test_zero_abundance_gene(self, tmp_path, update_golden):
        """Gene with 0 abundance — must produce zero counts everywhere."""
        sc = Scenario("golden_zero", genome_length=8000, seed=SEED,
                      work_dir=tmp_path / "golden_zero")
        sc.add_gene("g1", "+", [
            {"t_id": "t1",
             "exons": [(500, 800), (1200, 1500)],
             "abundance": 50},
        ])
        sc.add_gene("g_zero", "-", [
            {"t_id": "t_zero",
             "exons": [(3000, 3300), (3700, 4000)],
             "abundance": 0},
        ])
        result = sc.build_oracle(n_fragments=N_FRAGS, sim_config=_SIM_SS100)
        _run_golden_test("zero_abundance", result, result.index, update_golden)
        sc.cleanup()

    def test_many_short_exons(self, tmp_path, update_golden):
        """Gene with 6 short exons — tests many-splice-junction scoring."""
        sc = Scenario("golden_short_exons", genome_length=10000, seed=SEED,
                      work_dir=tmp_path / "golden_short_exons")
        sc.add_gene("g1", "+", [
            {"t_id": "t1",
             "exons": [(500, 650), (1000, 1150), (1500, 1650),
                       (2000, 2150), (2500, 2650), (3000, 3150)],
             "abundance": 80},
        ])
        sc.add_gene("g_h", "-", [
            {"t_id": "t_h", "exons": [(6000, 6500)], "abundance": 30},
        ])
        result = sc.build_oracle(n_fragments=N_FRAGS, sim_config=_SIM_SS100)
        _run_golden_test("many_short_exons", result, result.index, update_golden)
        sc.cleanup()

    def test_wide_intron(self, tmp_path, update_golden):
        """Gene with very wide intron — tests gDNA footprint scoring."""
        sc = Scenario("golden_wide_intron", genome_length=30000, seed=SEED,
                      work_dir=tmp_path / "golden_wide_intron")
        sc.add_gene("g1", "+", [
            {"t_id": "t1",
             "exons": [(500, 800), (20000, 20300)],
             "abundance": 50},
        ])
        sc.add_gene("g_h", "-", [
            {"t_id": "t_h", "exons": [(25000, 25500)], "abundance": 30},
        ])
        result = sc.build_oracle(n_fragments=N_FRAGS, sim_config=_SIM_SS100)
        _run_golden_test("wide_intron", result, result.index, update_golden)
        sc.cleanup()

    def test_higher_fragment_count(self, tmp_path, update_golden):
        """3000 fragments — tests scaling beyond the default N_FRAGS."""
        sc = Scenario("golden_3k", genome_length=10000, seed=SEED,
                      work_dir=tmp_path / "golden_3k")
        sc.add_gene("g1", "+", [
            {"t_id": "t1",
             "exons": [(500, 800), (1200, 1500)],
             "abundance": 50},
            {"t_id": "t2",
             "exons": [(500, 800), (1800, 2100)],
             "abundance": 30},
        ])
        sc.add_gene("g2", "-", [
            {"t_id": "t3",
             "exons": [(4000, 4300), (4700, 5000)],
             "abundance": 40},
        ])
        result = sc.build_oracle(n_fragments=3000, sim_config=_SIM_SS100)
        _run_golden_test("higher_frag_count_3k", result, result.index, update_golden)
        sc.cleanup()
