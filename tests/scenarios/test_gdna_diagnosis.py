"""Pure-mRNA baseline regression for gDNA over-absorption.

Scenario: 10 kb random genome with one 2-exon transcript T1(+):
  - Exon 1: [1000, 2000) — 1000 bp
  - Exon 2: [4000, 5000) — 1000 bp
  - Intron: [2000, 4000) — 2000 bp
  - Total exonic: 2000 bp, total span: 4000 bp (40% of genome)

Locks the control case: with a pure-mRNA library the pipeline must attribute
~all fragments to T1 and call ~0 gDNA. (The wider gDNA/nRNA/strand sweep that
once lived here was an investigation harness with no assertions; it has been
removed — use scripts/debug/ tooling for that.)
"""

from rigel.config import (
    BamScanConfig,
    EMConfig,
    PipelineConfig,
)
from rigel.pipeline import run_pipeline
from rigel.sim import GDNAConfig, Scenario, ReadSimConfig, run_benchmark

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
N_FRAGMENTS = 1000
SIM_SEED = 42
PIPELINE_SEED = 42
GENOME_LENGTH = 10_000
TRANSCRIPT_ABUNDANCE = 100.0  # Dominant transcript


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _sim_config(*, strand_specificity: float = 1.0, seed: int = SIM_SEED):
    return ReadSimConfig(
        frag_mean=200,
        frag_std=30,
        frag_min=80,
        frag_max=450,
        read_length=100,
        strand_specificity=strand_specificity,
        seed=seed,
    )


def _gdna_config(abundance: float) -> GDNAConfig | None:
    if abundance == 0:
        return None
    return GDNAConfig(
        abundance=abundance,
        frag_mean=350,
        frag_std=100,
        frag_min=100,
        frag_max=1000,
    )


def _build_scenario(tmp_path):
    """Build the 10kb scenario with one 2-exon transcript."""
    scenario = Scenario(
        "gdna_diag",
        genome_length=GENOME_LENGTH,
        seed=SIM_SEED,
        work_dir=tmp_path,
    )
    scenario.add_gene(
        gene_id="TRAPS: no-magic-numbers",
        strand="+",
        transcripts=[
            {
                "t_id": "T1",
                "exons": [(1000, 2000), (4000, 5000)],
                "abundance": TRANSCRIPT_ABUNDANCE,
            },
        ],
    )
    return scenario


def _run_diagnostic(
    tmp_path,
    *,
    gdna_abundance: float = 0,
    nrna_abundance: float = 0,
    strand_specificity: float = 1.0,
    n_fragments: int = N_FRAGMENTS,
    label: str = "",
):
    """Run the full pipeline and return detailed diagnostic info.

    Returns (bench, pipeline_result, scenario_result, diagnostics_dict).
    """
    scenario = _build_scenario(tmp_path)
    sc = _sim_config(strand_specificity=strand_specificity)
    gdna = _gdna_config(gdna_abundance)

    result = scenario.build_oracle(
        n_fragments=n_fragments,
        sim_config=sc,
        gdna_config=gdna,
        nrna_abundance=nrna_abundance,
    )

    config = PipelineConfig(
        em=EMConfig(seed=PIPELINE_SEED),
        scan=BamScanConfig(sj_strand_tag="auto"),
    )
    pr = run_pipeline(result.bam_path, result.index, config=config)
    bench = run_benchmark(result, pr, scenario_name=label or "gdna_diag")
    return bench, pr, result


class TestGDNADiagnosis:
    """Regression: a pure-mRNA library must not siphon mRNA into gDNA."""

    def test_pure_mrna_baseline(self, tmp_path):
        """Pure mRNA (no gDNA, no nRNA, perfect strand).

        The control case: every fragment is mRNA for T1, so the pipeline
        should attribute ~all fragments to T1 and call ~0 gDNA / nRNA.
        """
        bench, pr, result = _run_diagnostic(
            tmp_path,
            gdna_abundance=0,
            nrna_abundance=0,
            strand_specificity=1.0,
            label="pure_mRNA_ss1.0",
        )
        t1 = next(t for t in bench.transcripts if t.t_id == "T1")
        assert t1.abs_diff <= 5, f"T1 mRNA: expected={t1.expected}, observed={t1.observed:.0f}"
        assert bench.n_gdna_pipeline <= 3, f"Spurious gDNA: {bench.n_gdna_pipeline:.0f}"
