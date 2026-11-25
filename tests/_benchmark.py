"""Score ``run_pipeline`` counts against a scenario's read-name truth, per transcript and per pool.

The truth is read off the simulated read names (``rigel.sim.read_name``): from the R1 FASTQ when the scenario wrote
one, which holds every fragment aligned or not, else from the oracle BAM.

Usage
-----
>>> result = scenario.build(n_fragments=100)
>>> pipeline_result = run_pipeline(result.bam_path, result.index)
>>> bench = run_benchmark(result, pipeline_result)
>>> assert bench.all_exact, bench.summary()
"""

import gzip
from collections import Counter
from collections.abc import Iterable
from dataclasses import dataclass, field

import pysam

from rigel.pipeline import PipelineResult
from rigel.sim.read_name import parse_origin


def _read_names(scenario_result) -> Iterable[str]:
    """One read name per simulated fragment: the R1 FASTQ's when it exists, else the BAM's."""
    fastq = scenario_result.fastq_r1
    if fastq is not None and fastq.exists():
        with gzip.open(fastq, "rt") if str(fastq).endswith(".gz") else open(fastq) as handle:
            yield from (line for number, line in enumerate(handle) if number % 4 == 0)
        return
    with pysam.AlignmentFile(str(scenario_result.bam_path), "rb", check_sq=False) as bam:
        yield from {read.query_name for read in bam}


def simulated_truth(scenario_result) -> tuple[dict[str, int], Counter]:
    """The scenario's truth: mature-RNA fragments per transcript id, and fragments per origin kind (``annotated``,
    ``nascent``, ``gdna``)."""
    by_transcript: Counter[str] = Counter()
    by_kind: Counter[str] = Counter()
    for name in _read_names(scenario_result):
        origin = parse_origin(name)
        by_kind[origin.kind] += 1
        if origin.kind == "annotated" and origin.transcript_id is not None:
            by_transcript[origin.transcript_id] += 1
    return dict(by_transcript), by_kind


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
        Pipeline-assigned count (summed across all count columns).
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
    n_fragments : int
        Fragments seen by the pipeline (from stats.n_fragments).
    n_intergenic : int
        Fragments classified as intergenic.
    n_gdna_expected : int
        Ground-truth gDNA fragment count (from simulation read names).
    n_gdna_pipeline : float
        gDNA count reported by the pipeline: the intergenic fragments plus ``estimator.gdna_em_count``.
    transcripts : list[TranscriptAccuracy]
        Per-transcript accuracy breakdown.
    stats_dict : dict
        Full PipelineStats as a dictionary for inspection.
    """

    scenario_name: str
    n_simulated: int
    n_fragments: int
    n_intergenic: int
    n_gdna_expected: int = 0
    n_gdna_pipeline: float = 0.0
    n_nascent_expected: int = 0
    n_nascent_pipeline: float = 0.0
    n_synthetic_observed: float = 0.0
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
        """Sum of observed counts across all transcripts (annotated RNA only)."""
        return sum(t.observed for t in self.transcripts)

    @property
    def total_rna_observed(self) -> float:
        """Annotated plus synthetic-span counts: every RNA fragment the pipeline assigned."""
        return self.total_observed + self.n_synthetic_observed

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
    def alignment_rate(self) -> float:
        """Fraction of simulated fragments that entered the pipeline."""
        if self.n_simulated == 0:
            return 0.0
        return self.n_fragments / self.n_simulated

    @property
    def gdna_abs_diff(self) -> float:
        """Absolute difference between expected and pipeline gDNA counts."""
        return abs(self.n_gdna_pipeline - self.n_gdna_expected)

    @property
    def nascent_abs_diff(self) -> float:
        """Absolute difference between expected and pipeline nascent RNA counts."""
        return abs(self.n_nascent_pipeline - self.n_nascent_expected)

    # -- Display --------------------------------------------------------------

    def summary(self) -> str:
        """Human-readable summary string."""
        boundaries = [
            f"Benchmark: {self.scenario_name}",
            f"  Simulated: {self.n_simulated} fragments",
            f"  Pipeline fragments: {self.n_fragments} (alignment rate: {self.alignment_rate:.1%})",
            f"  Intergenic: {self.n_intergenic}",
            f"  Total expected: {self.total_expected}, Total observed: {self.total_observed:.0f}",
        ]
        if self.n_gdna_expected > 0 or self.n_gdna_pipeline > 0:
            boundaries.append(
                f"  gDNA: expected={self.n_gdna_expected}, "
                f"pipeline={self.n_gdna_pipeline:.0f}, "
                f"diff={self.gdna_abs_diff:.0f}"
            )
        if self.n_nascent_expected > 0 or self.n_nascent_pipeline > 0:
            boundaries.append(
                f"  nascent RNA: expected={self.n_nascent_expected}, "
                f"pipeline={self.n_nascent_pipeline:.0f}, "
                f"diff={self.nascent_abs_diff:.0f}"
            )
        boundaries.extend(
            [
                f"  All exact: {self.all_exact}, "
                f"Total abs error: {self.total_abs_error:.0f}, "
                f"Max rel error: {self.max_rel_error:.1%}",
                "  Per-transcript:",
            ]
        )
        for t in self.transcripts:
            boundaries.append(f"    {t}")
        return "\n".join(boundaries)


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
    simulated_counts, by_kind = simulated_truth(scenario_result)
    n_gdna_expected = by_kind["gdna"]
    n_nascent_expected = by_kind["nascent"]

    t_counts = pipeline_result.estimator.t_counts
    observed_per_t = t_counts.sum(axis=1)
    stats = pipeline_result.stats

    # gDNA: the intergenic fragments plus the EM's genic and intergenic gDNA.
    n_gdna_pipeline = float(
        stats.n_intergenic
        + pipeline_result.estimator.gdna_em_count
        + pipeline_result.estimator.gdna_intergenic_em_count
    )

    # Nascent RNA: the synthetic spans' counts.
    n_nascent_pipeline = float(pipeline_result.estimator.synthetic_count)

    transcript_results: list[TranscriptAccuracy] = []
    annotated_t_indices: set[int] = set()
    for t in scenario_result.transcripts:
        expected = simulated_counts.get(t.t_id, 0)
        observed = float(observed_per_t[t.t_index])
        abs_diff = abs(observed - expected)
        rel_error = abs_diff / expected if expected > 0 else 0.0
        exact = observed == expected
        annotated_t_indices.add(t.t_index)

        transcript_results.append(
            TranscriptAccuracy(
                t_id=t.t_id,
                t_index=t.t_index,
                expected=expected,
                observed=observed,
                exact_match=exact,
                abs_diff=abs_diff,
                rel_error=rel_error,
            )
        )

    # The counts on every index row outside the scenario's transcripts, the synthetic spans, so every fragment is
    # accounted for.
    n_synthetic_observed = sum(
        float(observed_per_t[i]) for i in range(len(observed_per_t)) if i not in annotated_t_indices
    )

    return BenchmarkResult(
        scenario_name=scenario_name or "unnamed",
        n_simulated=scenario_result.n_simulated,
        n_fragments=stats.n_fragments,
        n_intergenic=stats.n_intergenic,
        n_gdna_expected=n_gdna_expected,
        n_gdna_pipeline=n_gdna_pipeline,
        n_nascent_expected=n_nascent_expected,
        n_nascent_pipeline=n_nascent_pipeline,
        n_synthetic_observed=n_synthetic_observed,
        transcripts=transcript_results,
        stats_dict=stats.to_dict(),
    )
