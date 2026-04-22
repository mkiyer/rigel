"""
Diagnostic test: reproduce and root-cause the gDNA over-absorption bug.

Scenario: 10kb random genome with one 2-exon transcript T1(+):
  - Exon 1: [1000, 2000) — 1000 bp
  - Exon 2: [4000, 5000) — 1000 bp
  - Intron:  [2000, 4000) — 2000 bp
  - Total exonic: 2000 bp, total span: 4000 bp (40% of genome)

Tests sweep mRNA-only, mRNA+gDNA, mRNA+nRNA, mRNA+gDNA+nRNA conditions
with detailed diagnostic output at every stage of the EB initialization
and EM pipeline.

Goal: understand why the gDNA model siphons reads from mRNA and nRNA.
"""

import logging

from rigel.config import (
    BamScanConfig,
    EMConfig,
    PipelineConfig,
)
from rigel.pipeline import run_pipeline
from rigel.sim import GDNAConfig, Scenario, SimConfig, run_benchmark

logger = logging.getLogger(__name__)

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
    return SimConfig(
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
        gene_id="G1",
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

    # --- Collect diagnostic info ---
    diag = _collect_diagnostics(pr, result, bench, label)
    return bench, pr, result, diag


def _collect_diagnostics(pr, result, bench, label):
    """Extract and format all diagnostic information."""
    est = pr.estimator
    stats = pr.stats

    diag = {}

    # --- Ground truth ---
    gt_mrna = result.ground_truth_auto()
    gt_gdna = result.ground_truth_gdna_auto()
    gt_nrna = result.ground_truth_nrna_auto()
    diag["gt_mrna"] = gt_mrna
    diag["gt_gdna"] = gt_gdna
    diag["gt_nrna"] = gt_nrna

    # --- Pipeline results ---
    diag["pipeline_mrna"] = {t.t_id: t.observed for t in bench.transcripts}
    diag["pipeline_gdna"] = bench.n_gdna_pipeline
    diag["pipeline_nrna"] = bench.n_nrna_pipeline

    # --- Per-transcript unambig counts ---
    uc = est.unambig_counts  # (N_t, 6) array
    for t in result.transcripts:
        ti = t.t_index
        diag[f"unambig_counts_{t.t_id}"] = uc[ti].tolist()

    # --- Strand model ---
    sm = pr.strand_models
    diag["strand_specificity"] = sm.strand_specificity
    diag["exonic_spliced_ss"] = sm.exonic_spliced.strand_specificity
    diag["exonic_ss"] = sm.exonic.strand_specificity
    diag["intergenic_ss"] = sm.intergenic.strand_specificity

    # --- Transcript geometry ---
    for t in result.transcripts:
        ti = t.t_index
        if est._exonic_lengths is not None:
            diag[f"exonic_length_{t.t_id}"] = float(est._exonic_lengths[ti])
        if est._transcript_spans is not None:
            diag[f"transcript_span_{t.t_id}"] = float(est._transcript_spans[ti])
        diag[f"t_eff_len_{t.t_id}"] = float(est._t_eff_len[ti])

    # --- Locus results ---
    diag["locus_results"] = est.locus_results

    # --- Pipeline stats ---
    diag["n_fragments"] = stats.n_fragments
    diag["n_intergenic"] = stats.n_intergenic
    diag["n_chimeric"] = stats.n_chimeric
    diag["n_gdna_em"] = stats.n_gdna_em

    # --- gDNA contamination rate ---
    diag["gdna_contamination_rate"] = est.gdna_contamination_rate
    diag["gdna_total"] = est.gdna_em_count
    diag["gdna_em_count"] = est.gdna_em_count

    return diag


def _print_diagnostics(diag, label=""):
    """Pretty-print diagnostics to the logger."""
    lines = [
        "",
        "=" * 70,
        f"  DIAGNOSTIC REPORT: {label}",
        "=" * 70,
        "",
        "--- Ground Truth ---",
        f"  mRNA fragments per transcript: {diag['gt_mrna']}",
        f"  gDNA fragments: {diag['gt_gdna']}",
        f"  nRNA fragments: {diag['gt_nrna']}",
        "",
        "--- Pipeline Output ---",
        f"  mRNA pipeline: {diag['pipeline_mrna']}",
        f"  gDNA pipeline: {diag['pipeline_gdna']:.1f}",
        f"  nRNA pipeline: {diag['pipeline_nrna']:.1f}",
        "",
        "--- Stats ---",
        f"  n_fragments: {diag['n_fragments']}",
        f"  n_intergenic: {diag['n_intergenic']}",
        f"  n_chimeric: {diag['n_chimeric']}",
        f"  gdna_total: {diag['gdna_total']:.1f}",
        f"  gdna_em_count: {diag['gdna_em_count']:.1f}",
        f"  gdna_contamination_rate: {diag['gdna_contamination_rate']:.4f}",
        "",
        "--- Strand Model ---",
        f"  strand_specificity: {diag['strand_specificity']:.4f}",
        f"  exonic_spliced_ss: {diag['exonic_spliced_ss']:.4f}",
        f"  exonic_ss: {diag['exonic_ss']:.4f}",
        f"  intergenic_ss: {diag['intergenic_ss']:.4f}",
    ]

    # Per-transcript details
    for key in sorted(diag.keys()):
        if key.startswith("unambig_counts_"):
            t_id = key.replace("unambig_counts_", "")
            lines.append("")
            lines.append(f"--- Transcript {t_id} ---")
            lines.append(f"  unambig_counts: {diag[key]}")

            lines.append(
                f"  geometry: exonic_len={diag.get(f'exonic_length_{t_id}', '?')}, "
                f"span={diag.get(f'transcript_span_{t_id}', '?')}, "
                f"eff_len={diag.get(f't_eff_len_{t_id}', '?')}"
            )

    # Locus results
    lines.append("")
    lines.append("--- Locus EM Results ---")
    for lr in diag.get("locus_results", []):
        lines.append(
            f"  Locus {lr['locus_id']}: "
            f"n_t={lr['n_transcripts']}, n_genes={lr['n_genes']}, "
            f"n_em_frags={lr['n_em_fragments']}, "
            f"rna_total(mRNA+nRNA)={lr['rna_total']:.1f}, "
            f"gDNA={lr['gdna']:.1f}, alpha_gdna={lr['alpha_gdna']:.2f}, alpha_rna={lr['alpha_rna']:.2f}"
        )

    # Error analysis
    lines.append("")
    lines.append("--- Error Summary ---")
    for t_id, expected in diag["gt_mrna"].items():
        observed = diag["pipeline_mrna"].get(t_id, 0)
        diff = observed - expected
        lines.append(f"  {t_id}: expected={expected}, observed={observed:.1f}, diff={diff:+.1f}")
    gdna_diff = diag["pipeline_gdna"] - diag["gt_gdna"]
    lines.append(
        f"  gDNA: expected={diag['gt_gdna']}, "
        f"pipeline={diag['pipeline_gdna']:.1f}, diff={gdna_diff:+.1f}"
    )
    nrna_diff = diag["pipeline_nrna"] - diag["gt_nrna"]
    lines.append(
        f"  nRNA: expected={diag['gt_nrna']}, "
        f"pipeline={diag['pipeline_nrna']:.1f}, diff={nrna_diff:+.1f}"
    )

    lines.append("=" * 70)
    logger.warning("\n".join(lines))


# ---------------------------------------------------------------------------
# Test cases
# ---------------------------------------------------------------------------


class TestGDNADiagnosis:
    """Diagnostic suite for gDNA over-absorption investigation."""

    def test_pure_mrna_baseline(self, tmp_path):
        """Pure mRNA (no gDNA, no nRNA, perfect strand).

        This is the control: all fragments are mRNA for T1.
        The pipeline should attribute ~all fragments to T1 mRNA.
        gDNA and nRNA should be ~0.
        """
        bench, pr, result, diag = _run_diagnostic(
            tmp_path,
            gdna_abundance=0,
            nrna_abundance=0,
            strand_specificity=1.0,
            label="pure_mRNA_ss1.0",
        )
        _print_diagnostics(diag, "pure_mRNA_ss1.0")

        # All fragments should be mRNA
        t1 = next(t for t in bench.transcripts if t.t_id == "T1")
        assert t1.abs_diff <= 5, f"T1 mRNA: expected={t1.expected}, observed={t1.observed:.0f}"
        assert bench.n_gdna_pipeline <= 3, f"Spurious gDNA: {bench.n_gdna_pipeline:.0f}"

    def test_mrna_plus_gdna_ss1(self, tmp_path):
        """mRNA + gDNA contamination with perfect strandedness.

        This is the primary bug reproduction: gDNA should be separated
        cleanly by the strand signal (SS=1.0), but the EB prior may
        be over-estimating the gDNA rate.
        """
        bench, pr, result, diag = _run_diagnostic(
            tmp_path,
            gdna_abundance=50,
            nrna_abundance=0,
            strand_specificity=1.0,
            label="mRNA+gDNA_ss1.0",
        )
        _print_diagnostics(diag, "mRNA+gDNA_ss1.0")

        t1 = next(t for t in bench.transcripts if t.t_id == "T1")

        # Key diagnostic: how much mRNA was siphoned to gDNA?
        mrna_loss = t1.expected - t1.observed
        gdna_excess = bench.n_gdna_pipeline - bench.n_gdna_expected

        logger.warning(
            f"\n  >>> mRNA loss = {mrna_loss:.0f} frags "
            f"({mrna_loss / max(t1.expected, 1) * 100:.1f}%)\n"
            f"  >>> gDNA excess = {gdna_excess:.0f} frags "
            f"({gdna_excess / max(bench.n_gdna_expected, 1) * 100:.1f}%)"
        )

    def test_mrna_plus_gdna_ss09(self, tmp_path):
        """mRNA + gDNA with SS=0.9 (realistic strandedness)."""
        bench, pr, result, diag = _run_diagnostic(
            tmp_path,
            gdna_abundance=50,
            nrna_abundance=0,
            strand_specificity=0.9,
            label="mRNA+gDNA_ss0.9",
        )
        _print_diagnostics(diag, "mRNA+gDNA_ss0.9")

        t1 = next(t for t in bench.transcripts if t.t_id == "T1")
        mrna_loss = t1.expected - t1.observed
        gdna_excess = bench.n_gdna_pipeline - bench.n_gdna_expected

        logger.warning(
            f"\n  >>> mRNA loss = {mrna_loss:.0f} frags "
            f"({mrna_loss / max(t1.expected, 1) * 100:.1f}%)\n"
            f"  >>> gDNA excess = {gdna_excess:.0f} frags "
            f"({gdna_excess / max(bench.n_gdna_expected, 1) * 100:.1f}%)"
        )

    def test_mrna_plus_nrna_ss1(self, tmp_path):
        """mRNA + nRNA (no gDNA) with perfect strandedness.

        nRNA produces unspliced reads across the intron. The pipeline
        should correctly separate mRNA from nRNA.
        """
        bench, pr, result, diag = _run_diagnostic(
            tmp_path,
            gdna_abundance=0,
            nrna_abundance=30,
            strand_specificity=1.0,
            label="mRNA+nRNA_ss1.0",
        )
        _print_diagnostics(diag, "mRNA+nRNA_ss1.0")

        t1 = next(t for t in bench.transcripts if t.t_id == "T1")
        mrna_loss = t1.expected - t1.observed
        nrna_diff = bench.n_nrna_pipeline - bench.n_nrna_expected

        logger.warning(
            f"\n  >>> mRNA diff = {-mrna_loss:.0f} frags\n  >>> nRNA diff = {nrna_diff:.0f} frags"
        )

    def test_mrna_gdna_nrna_ss1(self, tmp_path):
        """Full 3-way: mRNA + gDNA + nRNA with perfect strandedness.

        This is the worst-case scenario from the benchmark: all three
        sources compete and the gDNA model may siphon from both.
        """
        bench, pr, result, diag = _run_diagnostic(
            tmp_path,
            gdna_abundance=50,
            nrna_abundance=30,
            strand_specificity=1.0,
            label="mRNA+gDNA+nRNA_ss1.0",
        )
        _print_diagnostics(diag, "mRNA+gDNA+nRNA_ss1.0")

        t1 = next(t for t in bench.transcripts if t.t_id == "T1")
        mrna_loss = t1.expected - t1.observed
        gdna_excess = bench.n_gdna_pipeline - bench.n_gdna_expected
        nrna_diff = bench.n_nrna_pipeline - bench.n_nrna_expected

        logger.warning(
            f"\n  >>> mRNA loss = {mrna_loss:.0f} frags\n"
            f"  >>> gDNA excess = {gdna_excess:.0f} frags\n"
            f"  >>> nRNA diff = {nrna_diff:.0f} frags"
        )

    def test_mrna_gdna_nrna_ss09(self, tmp_path):
        """Full 3-way with SS=0.9 (realistic library)."""
        bench, pr, result, diag = _run_diagnostic(
            tmp_path,
            gdna_abundance=50,
            nrna_abundance=30,
            strand_specificity=0.9,
            label="mRNA+gDNA+nRNA_ss0.9",
        )
        _print_diagnostics(diag, "mRNA+gDNA+nRNA_ss0.9")

        t1 = next(t for t in bench.transcripts if t.t_id == "T1")
        mrna_loss = t1.expected - t1.observed
        gdna_excess = bench.n_gdna_pipeline - bench.n_gdna_expected
        nrna_diff = bench.n_nrna_pipeline - bench.n_nrna_expected

        logger.warning(
            f"\n  >>> mRNA loss = {mrna_loss:.0f} frags\n"
            f"  >>> gDNA excess = {gdna_excess:.0f} frags\n"
            f"  >>> nRNA diff = {nrna_diff:.0f} frags"
        )

    def test_high_gdna_ss1(self, tmp_path):
        """Heavy gDNA contamination (abundance=200) with perfect strand.

        Extreme case: gDNA dominates. Does the pipeline over-absorb mRNA?
        """
        bench, pr, result, diag = _run_diagnostic(
            tmp_path,
            gdna_abundance=200,
            nrna_abundance=0,
            strand_specificity=1.0,
            label="high_gDNA_ss1.0",
        )
        _print_diagnostics(diag, "high_gDNA_ss1.0")

        t1 = next(t for t in bench.transcripts if t.t_id == "T1")
        mrna_loss = t1.expected - t1.observed
        gdna_excess = bench.n_gdna_pipeline - bench.n_gdna_expected

        logger.warning(
            f"\n  >>> mRNA loss = {mrna_loss:.0f} frags "
            f"({mrna_loss / max(t1.expected, 1) * 100:.1f}%)\n"
            f"  >>> gDNA excess = {gdna_excess:.0f} frags"
        )

    def test_mrna_plus_gdna_ss065(self, tmp_path):
        """mRNA + gDNA with low strandedness (SS=0.65).

        The strand signal is weak here. The pipeline relies more on
        density estimation. How well does it perform?
        """
        bench, pr, result, diag = _run_diagnostic(
            tmp_path,
            gdna_abundance=50,
            nrna_abundance=0,
            strand_specificity=0.65,
            label="mRNA+gDNA_ss0.65",
        )
        _print_diagnostics(diag, "mRNA+gDNA_ss0.65")

        t1 = next(t for t in bench.transcripts if t.t_id == "T1")
        mrna_loss = t1.expected - t1.observed
        gdna_excess = bench.n_gdna_pipeline - bench.n_gdna_expected

        logger.warning(
            f"\n  >>> mRNA loss = {mrna_loss:.0f} frags "
            f"({mrna_loss / max(t1.expected, 1) * 100:.1f}%)\n"
            f"  >>> gDNA excess = {gdna_excess:.0f} frags"
        )
