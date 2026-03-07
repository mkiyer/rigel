"""
Shared fixtures and helpers for aligner-dependent scenario tests.

These tests require minimap2 and samtools in PATH.  They use
``Scenario.build()`` (FASTQ → alignment) rather than ``build_oracle()``.
"""

import logging
import shutil

import pytest

from rigel.config import EMConfig, PipelineConfig, BamScanConfig
from rigel.pipeline import run_pipeline
from rigel.sim import GDNAConfig, Scenario, SimConfig, run_benchmark

logger = logging.getLogger(__name__)

pytestmark = pytest.mark.skipif(
    shutil.which("minimap2") is None or shutil.which("samtools") is None,
    reason="minimap2 and/or samtools not found in PATH",
)


# =====================================================================
# Parameter grids (match oracle scenarios)
# =====================================================================

GDNA_LEVELS = [0, 20, 100]
STRAND_LEVELS = [0.65, 0.9, 1.0]
NRNA_LEVELS = [0, 30, 70]

STRESS_COMBOS = [
    (20, 30, 0.9),
    (100, 30, 0.65),
]
STRESS_IDS = ["g20n30s90", "g100n30s65"]

N_FRAGMENTS = 500
SIM_SEED = 42
PIPELINE_SEED = 42


# =====================================================================
# Helpers
# =====================================================================


def sim_config(*, strand_specificity: float = 1.0, seed: int = SIM_SEED):
    return SimConfig(
        frag_mean=200, frag_std=30, frag_min=80, frag_max=450,
        read_length=100, strand_specificity=strand_specificity, seed=seed,
    )


def gdna_config(abundance: float) -> GDNAConfig | None:
    if abundance == 0:
        return None
    return GDNAConfig(
        abundance=abundance, frag_mean=350, frag_std=100,
        frag_min=100, frag_max=1000,
    )


def build_and_run(scenario, *, n_fragments=N_FRAGMENTS,
                  gdna_abundance=0, strand_specificity=1.0,
                  nrna_abundance=0, include_multimap=False,
                  scenario_name=""):
    """Build via FASTQ+alignment and run the pipeline."""
    sc = sim_config(strand_specificity=strand_specificity)
    gdna = gdna_config(gdna_abundance)
    result = scenario.build(
        n_fragments=n_fragments, sim_config=sc,
        gdna_config=gdna, nrna_abundance=nrna_abundance,
    )
    config = PipelineConfig(
        em=EMConfig(seed=PIPELINE_SEED),
        scan=BamScanConfig(sj_strand_tag="ts", include_multimap=include_multimap),
    )
    pr = run_pipeline(result.bam_path, result.index, config=config)
    bench = run_benchmark(result, pr, scenario_name=scenario_name)
    logger.info("\n%s", bench.summary())
    return bench


# =====================================================================
# Assertion helpers
# =====================================================================


def assert_alignment(bench, min_rate=0.70):
    """Aligner may lose some reads — lower threshold than oracle."""
    assert bench.n_fragments > 0, "No fragments entered the pipeline"
    assert bench.alignment_rate > min_rate, (
        f"Low alignment rate: {bench.alignment_rate:.1%}"
    )


def assert_accountability(bench, tolerance=5):
    sd = bench.stats_dict
    mm_extra = sd.get("n_multimapper_alignments", 0) - sd.get(
        "n_multimapper_groups", 0
    )
    n_gated_out = sd.get("n_gated_out", 0)
    effective_fragments = bench.n_fragments - mm_extra
    total = (
        bench.total_observed + bench.n_nrna_pipeline
        + bench.n_gdna_pipeline + bench.n_chimeric
        + n_gated_out
    )
    assert abs(total - effective_fragments) <= tolerance, (
        f"Accountability gap: {abs(total - effective_fragments):.0f} "
        f"(total={total:.0f}, eff_frags={effective_fragments}, "
        f"mm_extra={mm_extra})"
    )


def assert_negative_control(bench, ctrl_id="t_ctrl", *,
                             gdna_abundance=0, strand_specificity=1.0):
    ctrl = next(t for t in bench.transcripts if t.t_id == ctrl_id)
    max_fp = 5
    if gdna_abundance > 0:
        max_fp += min(gdna_abundance, 60)
    if strand_specificity < 0.9:
        ss_gap = 0.9 - strand_specificity
        max_fp += round(ss_gap * 200 + gdna_abundance * ss_gap * 8)
    assert ctrl.observed <= max_fp, (
        f"Negative control {ctrl_id}: {ctrl.observed:.0f} counts "
        f"(limit={max_fp}, gdna={gdna_abundance}, "
        f"ss={strand_specificity})"
    )


def assert_gdna_accuracy(bench, gdna_abundance, max_rel_err=0.55):
    if gdna_abundance == 0:
        return
    for ta in bench.transcripts:
        if ta.t_id in ("t_ctrl", "t_helper"):
            continue
        assert ta.observed <= ta.expected + bench.n_gdna_expected + 5, (
            f"{ta.t_id}: observed={ta.observed:.0f} exceeds "
            f"RNA({ta.expected}) + gDNA({bench.n_gdna_expected})"
        )
    if bench.n_gdna_expected > 10:
        rel_err = bench.gdna_abs_diff / bench.n_gdna_expected
        assert rel_err < max_rel_err, (
            f"gDNA rel error: pipeline={bench.n_gdna_pipeline:.0f}, "
            f"expected={bench.n_gdna_expected}, rel={rel_err:.2f}"
        )
