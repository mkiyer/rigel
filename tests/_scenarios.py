"""The scenario suites' shared grid, builder and assertion helpers.

``tests/scenarios/`` builds every scenario with ``Scenario.build_oracle()``, whose alignments are perfect by
construction, so a failure is the pipeline's and never an aligner's, and the suite needs no external tools.
``tests/scenarios_aligned/`` goes through ``Scenario.build()`` — FASTQ, then a real alignment — so it sees the losses
and the multimapping an oracle BAM cannot produce, and it passes the looser tolerances that calls for. Each assertion
helper encodes a tolerance that has to mean the same thing in every scenario that uses it.
"""

import logging

from rigel.config import BamScanConfig, EMConfig, PipelineConfig
from rigel.pipeline import run_pipeline
from rigel.sim import GDNAConfig, ReadSimConfig

from _benchmark import run_benchmark

logger = logging.getLogger(__name__)

GDNA_LEVELS = [0, 20, 100]
STRAND_LEVELS = [0.65, 0.9, 1.0]
NRNA_LEVELS = [0, 30, 70]

# Stress combos: one moderate and one extreme
STRESS_COMBOS = [
    (20, 30, 0.9),
    (100, 30, 0.65),
]
STRESS_IDS = ["g20n30s90", "g100n30s65"]

N_FRAGMENTS = 500
SIM_SEED = 42
PIPELINE_SEED = 42

#: How many standard deviations of a fractional split a bound allows (:func:`split_noise`): a one-sided 3-sigma excess,
#: about one draw in 740, with the seeds fixed besides.
SPLIT_NOISE_Z = 3.0


def split_noise(n_fragments: int) -> float:
    """The most an unbiased fractional split can misplace, in fragments, at ``SPLIT_NOISE_Z``: two pools sharing ``m``
    fragments at a share ``f`` exchange a binomial ``sqrt(m·f(1 − f)) ≤ sqrt(m)/2`` of them, and no pool shares more
    than the library's ``n_fragments``."""
    return SPLIT_NOISE_Z * n_fragments**0.5 / 2.0


def sim_config(*, strand_specificity: float = 1.0, seed: int = SIM_SEED):
    """The read simulation every scenario uses."""
    return ReadSimConfig(
        frag_mean=200,
        frag_std=30,
        frag_min=80,
        frag_max=450,
        read_length=100,
        strand_specificity=strand_specificity,
        seed=seed,
    )


def gdna_config(abundance: float) -> GDNAConfig | None:
    """The gDNA simulation at *abundance*, or None when it is 0."""
    if abundance == 0:
        return None
    return GDNAConfig(
        abundance=abundance,
        frag_mean=350,
        frag_std=100,
        frag_min=100,
        frag_max=1000,
    )


def build_and_run(
    scenario,
    *,
    aligned=False,
    n_fragments=N_FRAGMENTS,
    gdna_abundance=0,
    strand_specificity=1.0,
    nrna_abundance=0,
    include_multimap=True,
    scenario_name="",
    em_config=None,
):
    """Build *scenario* — an oracle BAM, or FASTQ and an alignment when *aligned* — and run the pipeline (the shipped
    EM unless *em_config* names another); returns the benchmark against ground truth."""
    build = scenario.build if aligned else scenario.build_oracle
    result = build(
        n_fragments=n_fragments,
        sim_config=sim_config(strand_specificity=strand_specificity),
        gdna_config=gdna_config(gdna_abundance),
        nrna_abundance=nrna_abundance,
    )
    config = PipelineConfig(
        em=em_config if em_config is not None else EMConfig(seed=PIPELINE_SEED),
        scan=BamScanConfig(
            sj_strand_tag="ts" if aligned else "auto", include_multimap=include_multimap
        ),
    )
    pr = run_pipeline(result.bam_path, result.index, config=config)
    bench = run_benchmark(result, pr, scenario_name=scenario_name)
    logger.info("\n%s", bench.summary())
    return bench


def assert_alignment(bench, min_rate=0.95):
    """Fragments entered the pipeline and the alignment rate clears *min_rate* (an oracle BAM is near perfect)."""
    assert bench.n_fragments > 0, "No fragments entered the pipeline"
    assert bench.alignment_rate > min_rate, f"Low alignment rate: {bench.alignment_rate:.1%}"


def assert_accountability(bench, tolerance=5):
    """Every fragment is accounted for: annotated, synthetic (nascent), gDNA or gated out."""
    n_gated_out = bench.stats_dict.get("n_gated_out", 0)
    # n_synthetic_observed already includes the synthetic spans' counts,
    # so n_nrna_pipeline is not added separately.
    total = bench.total_observed + bench.n_synthetic_observed + bench.n_gdna_pipeline + n_gated_out
    gap = abs(total - bench.n_fragments)
    assert gap <= tolerance, (
        f"Accountability gap: {gap:.0f} (total={total:.0f}, n_fragments={bench.n_fragments})"
    )


def assert_transcript_accuracy(bench, max_abs_diff=3, exclude_ids=("t_ctrl", "t_helper")):
    """Per-transcript mRNA accuracy (skips controls/helpers)."""
    for ta in bench.transcripts:
        if ta.t_id in exclude_ids:
            continue
        assert ta.abs_diff <= max_abs_diff, (
            f"{ta.t_id}: expected={ta.expected}, observed={ta.observed:.0f}, diff={ta.abs_diff:.0f}"
        )


def assert_negative_control(bench, ctrl_id="t_ctrl", *, gdna_abundance=0, strand_specificity=1.0):
    """A negative-control transcript must carry near-zero mRNA counts.

    The control expresses nothing, so every fragment it holds is one the EM split onto it from a pool it shares
    fragments with: gDNA landing on its exons, and at imperfect strand specificity a neighbour's mis-stranded RNA. An
    unbiased split gives it none on average and at most :func:`split_noise` of the library at ``SPLIT_NOISE_Z``,
    whatever the contamination or the strand specificity, which only change how many fragments are shared. Measured
    over every scenario sweep, the control has held at most one fragment.
    """
    ctrl = next(t for t in bench.transcripts if t.t_id == ctrl_id)
    max_fp = split_noise(bench.n_fragments)
    assert ctrl.observed <= max_fp, (
        f"Negative control {ctrl_id}: {ctrl.observed:.0f} counts "
        f"(limit={max_fp:.0f}, gdna={gdna_abundance}, "
        f"ss={strand_specificity})"
    )


def assert_gdna_accuracy(bench, gdna_abundance, max_rel_err=0.60):
    """No transcript holds more than its RNA and all the gDNA, and the gDNA count is within *max_rel_err*."""
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


def assert_nrna_detected(bench, nrna_abundance):
    """When nascent RNA is present and not swamped by gDNA, the pipeline detects some."""
    if nrna_abundance == 0 or bench.n_nrna_expected <= 5:
        return
    if bench.n_gdna_expected > 3 * bench.n_nrna_expected:
        return
    assert bench.n_nrna_pipeline > 0, f"No nRNA detected but expected {bench.n_nrna_expected:.0f}"
