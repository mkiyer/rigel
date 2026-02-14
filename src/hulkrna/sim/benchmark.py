"""
hulkrna.sim.benchmark — Accuracy benchmarking for simulation scenarios.

Compares observed counts from ``run_pipeline`` against ground-truth
fragment assignments from the read simulator.  Produces per-transcript
accuracy metrics suitable for both deterministic correctness testing
and statistical accuracy benchmarking.

Usage
-----
>>> result = scenario.build(n_fragments=100)
>>> pipeline_result = run_pipeline(result.bam_path, result.index, sj_strand_tag="ts")
>>> bench = run_benchmark(result, pipeline_result)
>>> assert bench.all_exact, bench.summary()
"""

from dataclasses import dataclass, field

import numpy as np

from ..pipeline import PipelineResult

__all__ = [
    "TranscriptAccuracy",
    "BenchmarkResult",
    "run_benchmark",
]


@dataclass(slots=True)
class TranscriptAccuracy:
    """Per-transcript accuracy metrics.

    Attributes
    ----------
    t_id : str
        Transcript identifier.
    t_index : int
        Transcript index in the count array.
    expected : int
        Ground-truth fragment count (from simulation read names).
    observed : float
        Pipeline-assigned count (summed across all 12 count types).
    exact_match : bool
        True if ``observed == expected``.
    abs_diff : float
        ``|observed - expected|``.
    rel_error : float
        ``abs_diff / expected``, or 0.0 if ``expected == 0``.
    """

    t_id: str
    t_index: int
    expected: int
    observed: float
    exact_match: bool
    abs_diff: float
    rel_error: float

    def __str__(self) -> str:
        return (
            f"{self.t_id}: expected={self.expected}, observed={self.observed:.0f}, "
            f"diff={self.abs_diff:+.0f}, rel_err={self.rel_error:.1%}"
            f"{' ✓' if self.exact_match else ''}"
        )


@dataclass
class BenchmarkResult:
    """Aggregate accuracy metrics for a simulation scenario.

    Attributes
    ----------
    scenario_name : str
        Name of the scenario.
    n_simulated : int
        Total fragments passed to the read simulator.
    n_simulated_per_transcript : dict[str, int]
        Ground truth from FASTQ — fragments simulated per transcript.
    n_aligned_per_transcript : dict[str, int]
        Ground truth from BAM — fragments that appear in the BAM
        (aligned or present as read names, regardless of mapping).
    n_fragments : int
        Fragments seen by the pipeline (from stats.n_fragments).
    n_intergenic : int
        Fragments classified as intergenic.
    n_chimeric : int
        Fragments classified as chimeric.
    n_gdna_expected : int
        Ground-truth gDNA fragment count (from simulation read names).
    n_gdna_pipeline : float
        gDNA count reported by the pipeline (``counter.gdna_total``).
    transcripts : list[TranscriptAccuracy]
        Per-transcript accuracy breakdown.
    stats_dict : dict
        Full PipelineStats as a dictionary for inspection.
    """

    scenario_name: str
    n_simulated: int
    n_simulated_per_transcript: dict[str, int]
    n_aligned_per_transcript: dict[str, int]
    n_fragments: int
    n_intergenic: int
    n_chimeric: int
    n_gdna_expected: int = 0
    n_gdna_pipeline: float = 0.0
    transcripts: list[TranscriptAccuracy] = field(default_factory=list)
    stats_dict: dict = field(default_factory=dict)

    # -- Derived properties ---------------------------------------------------

    @property
    def all_exact(self) -> bool:
        """True if every transcript has an exact count match."""
        return all(t.exact_match for t in self.transcripts)

    @property
    def total_expected(self) -> int:
        """Sum of expected counts across all transcripts."""
        return sum(t.expected for t in self.transcripts)

    @property
    def total_observed(self) -> float:
        """Sum of observed counts across all transcripts."""
        return sum(t.observed for t in self.transcripts)

    @property
    def total_abs_error(self) -> float:
        """Sum of absolute differences across all transcripts."""
        return sum(t.abs_diff for t in self.transcripts)

    @property
    def max_rel_error(self) -> float:
        """Maximum relative error across all transcripts."""
        if not self.transcripts:
            return 0.0
        return max(t.rel_error for t in self.transcripts)

    @property
    def n_counted(self) -> float:
        """Total counts assigned by the pipeline."""
        return self.total_observed

    @property
    def alignment_rate(self) -> float:
        """Fraction of simulated fragments that entered the pipeline."""
        if self.n_simulated == 0:
            return 0.0
        return self.n_fragments / self.n_simulated

    @property
    def counting_rate(self) -> float:
        """Fraction of pipeline fragments that received a count."""
        if self.n_fragments == 0:
            return 0.0
        return self.total_observed / self.n_fragments

    @property
    def gdna_abs_diff(self) -> float:
        """Absolute difference between expected and pipeline gDNA counts."""
        return abs(self.n_gdna_pipeline - self.n_gdna_expected)

    # -- Display --------------------------------------------------------------

    def summary(self) -> str:
        """Human-readable summary string."""
        lines = [
            f"Benchmark: {self.scenario_name}",
            f"  Simulated: {self.n_simulated} fragments",
            f"  Pipeline fragments: {self.n_fragments} "
            f"(alignment rate: {self.alignment_rate:.1%})",
            f"  Intergenic: {self.n_intergenic}, Chimeric: {self.n_chimeric}",
            f"  Total expected: {self.total_expected}, "
            f"Total observed: {self.total_observed:.0f}",
        ]
        if self.n_gdna_expected > 0 or self.n_gdna_pipeline > 0:
            lines.append(
                f"  gDNA: expected={self.n_gdna_expected}, "
                f"pipeline={self.n_gdna_pipeline:.0f}, "
                f"diff={self.gdna_abs_diff:.0f}"
            )
        lines.extend([
            f"  All exact: {self.all_exact}, "
            f"Total abs error: {self.total_abs_error:.0f}, "
            f"Max rel error: {self.max_rel_error:.1%}",
            "  Per-transcript:",
        ])
        for t in self.transcripts:
            lines.append(f"    {t}")
        return "\n".join(lines)


def run_benchmark(
    scenario_result,
    pipeline_result: PipelineResult,
    *,
    scenario_name: str = "",
) -> BenchmarkResult:
    """Compare pipeline counts against simulation ground truth.

    Parameters
    ----------
    scenario_result : ScenarioResult
        The simulation artifacts (with BAM, FASTQ, transcripts).
    pipeline_result : PipelineResult
        Output of ``run_pipeline()``.
    scenario_name : str
        Optional label for the benchmark result.

    Returns
    -------
    BenchmarkResult
    """
    # Ground truth: actual simulated fragments per transcript (from FASTQ)
    simulated_counts = scenario_result.ground_truth_from_fastq()

    # Ground truth: fragments present in BAM (may be fewer if unmapped)
    aligned_counts = scenario_result.ground_truth_counts()

    # Ground truth: gDNA fragment count
    n_gdna_expected = scenario_result.ground_truth_gdna_count()

    # Observed: per-transcript total counts from the pipeline
    t_counts = pipeline_result.counter.t_counts  # (N_t, 12) array
    observed_per_t = t_counts.sum(axis=1)  # total per transcript

    # Pipeline gDNA count
    n_gdna_pipeline = float(pipeline_result.counter.gdna_total)

    # Build t_id → t_index mapping
    t_id_to_index = {t.t_id: t.t_index for t in scenario_result.transcripts}

    # Stats
    stats = pipeline_result.stats

    # Per-transcript accuracy
    transcript_results: list[TranscriptAccuracy] = []
    for t in scenario_result.transcripts:
        expected = simulated_counts.get(t.t_id, 0)
        observed = float(observed_per_t[t.t_index])
        abs_diff = abs(observed - expected)
        rel_error = abs_diff / expected if expected > 0 else 0.0
        exact = observed == expected

        transcript_results.append(TranscriptAccuracy(
            t_id=t.t_id,
            t_index=t.t_index,
            expected=expected,
            observed=observed,
            exact_match=exact,
            abs_diff=abs_diff,
            rel_error=rel_error,
        ))

    return BenchmarkResult(
        scenario_name=scenario_name or "unnamed",
        n_simulated=scenario_result.n_simulated,
        n_simulated_per_transcript=simulated_counts,
        n_aligned_per_transcript=aligned_counts,
        n_fragments=stats.n_fragments,
        n_intergenic=stats.n_intergenic,
        n_chimeric=stats.n_chimeric,
        n_gdna_expected=n_gdna_expected,
        n_gdna_pipeline=n_gdna_pipeline,
        transcripts=transcript_results,
        stats_dict=stats.to_dict(),
    )
