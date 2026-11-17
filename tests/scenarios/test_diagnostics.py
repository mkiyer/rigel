"""Focused end-to-end scenarios for particular algorithm behaviours,
rather than the full parameter sweeps the other scenario files run: minor-isoform
count collapse when one isoform is a strict exonic subset of another, read
retention for unspliced and spliced genes as strand specificity falls, and a
pure-mRNA control library in which the pipeline must attribute nearly every
fragment to RNA and call almost no gDNA.
"""

from rigel.config import EMConfig, PipelineConfig, BamScanConfig
from rigel.pipeline import run_pipeline
from rigel.sim import Scenario

from _benchmark import run_benchmark, simulated_truth

from _scenarios import sim_config, gdna_config, SIM_SEED, PIPELINE_SEED

N_FRAGMENTS = 1000

#: How many binomial standard deviations of the antisense share a strand-noise bound allows: a one-sided 3-sigma excess,
#: about one draw in 740, so a sweep's handful of arms fails by chance about once in a few hundred runs, and the seeds
#: are fixed besides.
STRAND_NOISE_Z = 3.0
GENOME_LENGTH = 10_000
TRANSCRIPT_ABUNDANCE = 100.0  # Dominant transcript


# ── isoform collapse when one isoform is an exonic subset of the other ─


class TestIsoformCollapse:
    """Minor isoform count collapse when t2 is a strict exonic subset of t1.

    Validates that effective-length normalization keeps t2 error < 50%.
    Also verifies recovery when t2 has its own unique exon.
    """

    def _make_two_isoform_scenario(self, tmp_path, t1_abundance, t2_abundance):
        sc = Scenario("iso_diag", genome_length=5000, seed=SIM_SEED, work_dir=tmp_path / "iso_diag")
        sc.add_gene(
            "g1",
            "+",
            [
                {
                    "t_id": "t1",
                    "exons": [(200, 500), (1000, 1300), (2000, 2300)],
                    "abundance": t1_abundance,
                },
                {"t_id": "t2", "exons": [(200, 500), (2000, 2300)], "abundance": t2_abundance},
            ],
        )
        return sc

    def test_equal_abundance(self, tmp_path):
        """1:1 isoform ratio: t2 should be within 50% of truth."""
        sc = self._make_two_isoform_scenario(tmp_path, 100, 100)
        try:
            result = sc.build_oracle(n_fragments=2000, sim_config=sim_config())
            gt, _ = simulated_truth(result)
            pr = run_pipeline(
                result.bam_path,
                result.index,
                config=PipelineConfig(
                    em=EMConfig(seed=PIPELINE_SEED),
                    scan=BamScanConfig(sj_strand_tag="auto"),
                ),
            )
            bench = run_benchmark(result, pr, scenario_name="iso_1_1")

            t2 = next(t for t in bench.transcripts if t.t_id == "t2")
            t2_err = abs(t2.observed - gt.get("t2", 0))
            truth_t2 = gt.get("t2", 0)
            assert t2_err < truth_t2 * 0.5, (
                f"t2 error {t2_err:.0f} exceeds 50% of truth ({truth_t2})"
            )
        finally:
            sc.cleanup()

    def test_unique_t2_exon_recovers(self, tmp_path):
        """When t2 has its own unique exon, the EM should recover."""
        sc = Scenario(
            "iso_unique_t2", genome_length=8000, seed=SIM_SEED, work_dir=tmp_path / "iso_unique_t2"
        )
        sc.add_gene(
            "g1",
            "+",
            [
                {"t_id": "t1", "exons": [(200, 500), (1000, 1300), (2000, 2300)], "abundance": 100},
                {"t_id": "t2", "exons": [(200, 500), (2000, 2300), (3500, 3800)], "abundance": 100},
            ],
        )
        try:
            result = sc.build_oracle(n_fragments=2000, sim_config=sim_config())
            gt, _ = simulated_truth(result)
            pr = run_pipeline(
                result.bam_path,
                result.index,
                config=PipelineConfig(
                    em=EMConfig(seed=PIPELINE_SEED),
                    scan=BamScanConfig(sj_strand_tag="auto"),
                ),
            )
            bench = run_benchmark(result, pr, scenario_name="iso_unique_t2")

            t1 = next(t for t in bench.transcripts if t.t_id == "t1")
            t2 = next(t for t in bench.transcripts if t.t_id == "t2")
            total_err = abs(t1.observed - gt["t1"]) + abs(t2.observed - gt["t2"])
            assert total_err < 0.5 * (gt["t1"] + gt["t2"]), (
                f"Total error {total_err:.0f} too high with unique exons"
            )
        finally:
            sc.cleanup()


# ── unspliced and spliced genes as strand specificity falls ────────────


class TestUnsplicedLowStrand:
    """At low SS, unspliced genes may lose reads to the gDNA component.

    Spliced genes survive because splice-junction strand evidence
    is independent of library strandedness.
    """

    def test_ss_gradient_unspliced(self, tmp_path):
        """Sweep SS for an unspliced gene in a library with no gDNA: at every SS it keeps its reads, to within what the
        strand channel's sampling noise could misread as gDNA.

        RNA's antisense rate is ``p = 1 − SS`` and gDNA's is ½, so ``n`` fragments whose antisense share sits ``z``
        binomial standard deviations above ``p`` read as a gDNA share of ``z·sqrt(p(1 − p)/n) / (½ − p)``. That is
        the most a zero-gDNA library can lose to a phantom gDNA call at ``z = STRAND_NOISE_Z``: about a quarter of the
        gene at SS 0.65, a twentieth at 0.9, and nothing at 1.0, where the strand channel has no noise.
        """
        sc = Scenario("ss_grad", genome_length=5000, seed=SIM_SEED, work_dir=tmp_path / "ss_grad")
        sc.add_gene(
            "g1",
            "+",
            [
                {"t_id": "t1", "exons": [(500, 1500)], "abundance": 100},
            ],
        )
        # Multi-exon gene so the strand model trains: the calibration requires spliced
        # reads (a pure single-exon library calibrates as unstranded), and the
        # trained κ_rna is what lets the unspliced gene's reads deconvolve as RNA.
        sc.add_gene(
            "g_train",
            "+",
            [
                {"t_id": "t_train", "exons": [(2500, 2800), (3200, 3500)], "abundance": 60},
            ],
        )

        for ss in [0.65, 0.9, 1.0]:
            result = sc.build_oracle(
                n_fragments=2000,
                sim_config=sim_config(strand_specificity=ss),
            )
            gt, _ = simulated_truth(result)
            pr = run_pipeline(
                result.bam_path,
                result.index,
                config=PipelineConfig(
                    em=EMConfig(seed=PIPELINE_SEED),
                    scan=BamScanConfig(sj_strand_tag="auto"),
                ),
            )
            bench = run_benchmark(result, pr, scenario_name=f"ss_{ss}")
            t1 = next(t for t in bench.transcripts if t.t_id == "t1")
            n = gt.get("t1", 0)
            p = 1.0 - ss
            phantom = STRAND_NOISE_Z * (p * (1.0 - p) / n) ** 0.5 / (0.5 - p)
            assert t1.observed >= n * (1.0 - phantom) - 1.0, (
                f"t1 too low at ss={ss}: {t1.observed:.0f} vs truth {n}; strand noise allows a "
                f"{phantom:.1%} phantom gDNA share"
            )

        sc.cleanup()

    def test_spliced_survives_low_ss(self, tmp_path):
        """Spliced gene retains reads even at SS=0.65."""
        sc = Scenario(
            "ss_spliced", genome_length=5000, seed=SIM_SEED, work_dir=tmp_path / "ss_spliced"
        )
        sc.add_gene(
            "g1",
            "+",
            [
                {"t_id": "t1", "exons": [(200, 500), (1000, 1300)], "abundance": 100},
            ],
        )
        try:
            result = sc.build_oracle(
                n_fragments=2000,
                sim_config=sim_config(strand_specificity=0.65),
            )
            gt, _ = simulated_truth(result)
            pr = run_pipeline(
                result.bam_path,
                result.index,
                config=PipelineConfig(
                    em=EMConfig(seed=PIPELINE_SEED),
                    scan=BamScanConfig(sj_strand_tag="auto"),
                ),
            )
            bench = run_benchmark(result, pr, scenario_name="ss_0.65_spliced")
            t1 = next(t for t in bench.transcripts if t.t_id == "t1")
            assert t1.observed >= gt.get("t1", 0) * 0.3, (
                f"Spliced gene lost too many reads at ss=0.65: "
                f"{t1.observed:.0f} vs truth {gt.get('t1', 0)}"
            )
        finally:
            sc.cleanup()


# ── a pure-mRNA library must not be over-absorbed into gDNA ────────────


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
    sc = sim_config(strand_specificity=strand_specificity)
    gdna = gdna_config(gdna_abundance)

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
    """A pure-mRNA library must not siphon mRNA into gDNA.

    One 2-exon transcript T1(+) on a 10 kb random genome: exon 1 at [1000, 2000),
    exon 2 at [4000, 5000), intron [2000, 4000). The 4000 bp span is 40% of the
    genome, so an over-absorbing gDNA model has plenty of room to show itself.
    """

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
