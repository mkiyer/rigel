"""The report substrate ``rigel quant`` writes, and the contract the report builder reads it under.

``summary.json`` keeps the fragment-length section stats-only; the per-bin histograms live in a tidy
``(category, length, count)`` companion feather whose per-category sums reconcile with the summary's observation
counts. The reported categories are the anchor, the two models and the five pure pools, and the modelled aggregates
are exact sums of the primitives beneath them, so two files describing one library agree.
"""

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pytest
import pandas as pd

from rigel import cli
from rigel.config import BamScanConfig, EMConfig, PipelineConfig
from rigel.pipeline import run_pipeline
from rigel.sim import GDNAConfig, ReadSimConfig, Scenario

SEED = 42


# ---------------------------------------------------------------------------
# Unit: the fragment-length split helper
# ---------------------------------------------------------------------------


def _fl_models(max_size: int = 1000):
    """An ``FLModels`` built straight from pool histograms, the shape production produces.

    Built through the EB kernel rather than through ``build_fl_models``, because that entry point takes
    a whole payload, which a unit test of the report should not fabricate.
    """
    from rigel.calibration.fl import fl_models_from_histograms
    from rigel.scan_payload import (
        N_FRAGMENT_POOLS,
        POOL_DNA_INTERGENIC,
        POOL_DNA_INTERGENIC_EXON,
        POOL_DNA_INTRONIC,
        POOL_DNA_INTRON_EXON,
        POOL_RNA_SPLICED,
    )

    pools = np.zeros((N_FRAGMENT_POOLS, max_size + 1), dtype=np.float64)
    pools[POOL_DNA_INTERGENIC, 100] = 2.0
    pools[POOL_DNA_INTRONIC, 150] = 1.0
    pools[POOL_DNA_INTRON_EXON, 120] = 3.0
    pools[POOL_DNA_INTERGENIC_EXON, 130] = 4.0
    pools[POOL_RNA_SPLICED, 200] = 5.0

    unconditional = np.zeros(max_size + 1, dtype=np.float64)
    unconditional[[100, 150, 200, 300]] = [2.0, 1.0, 3.0, 1.0]

    return fl_models_from_histograms(
        global_counts=unconditional,
        rna_counts=pools[POOL_RNA_SPLICED],
        gdna_counts=pools[[POOL_DNA_INTERGENIC, POOL_DNA_INTRONIC]].sum(axis=0),
        max_size=max_size,
        pool_counts=pools,
    )


def test_fragment_length_report_splits_lean_summary_from_histograms():
    summary, hist_df = cli.fragment_length_report(_fl_models())

    # summary is lean: stats only, never raw per-bin arrays
    assert "global" in summary
    for cat, s in summary.items():
        assert "histogram" not in s, f"{cat} still carries raw histogram"
        assert set(s) == {
            "n_observations",
            "mean",
            "std",
            "median",
            "mode",
            "max_size",
            "overflow_count",
            "overflow_fraction",
        }
    assert summary["global"]["n_observations"] == 7

    # tidy histogram table reconciles with the summary observation counts
    assert list(hist_df.columns) == ["category", "length", "count"]
    assert str(hist_df["length"].dtype) == "int32"
    for cat in hist_df["category"].unique():
        fsum = float(hist_df.loc[hist_df["category"] == cat, "count"].sum())
        assert abs(fsum - summary[cat]["n_observations"]) < 0.5


def test_the_reported_categories_are_the_anchor_the_two_models_and_five_pure_pools():
    """The report's output contract: the anchor, the two models, and the five pure pools.

    There is no per-``SpliceType`` histogram (``unspliced``, ``spliced_annot``, …) in the set: the five pure
    pools are a structural gDNA/RNA classification in the accumulator's one frame, and the two ``*_exon``
    rows are on-target gDNA, which nothing else surfaces.
    """
    summary, _ = cli.fragment_length_report(_fl_models())
    assert set(summary) == {
        "global",
        "gdna",
        "rna",
        "pool_gdna_intergenic",
        "pool_gdna_intronic",
        "pool_gdna_intron_exon",
        "pool_gdna_intergenic_exon",
        "pool_rna_spliced",
    }


def test_the_modelled_pools_are_sums_of_the_reported_pure_pools():
    """The aggregates and the primitives must reconcile, or the report shows two truths.

    ``rna`` is ``pool_rna_spliced``; ``gdna`` is the two contained rows and excludes the two crossing
    ones, because summing four differently-tilted pools reports a gDNA mean far above what the pure
    intergenic pool measures.
    """
    summary, _ = cli.fragment_length_report(_fl_models())
    n = {k: v["n_observations"] for k, v in summary.items()}
    assert n["rna"] == n["pool_rna_spliced"]
    assert n["gdna"] == n["pool_gdna_intergenic"] + n["pool_gdna_intronic"]
    assert n["gdna"] < n["gdna"] + n["pool_gdna_intron_exon"] + n["pool_gdna_intergenic_exon"], (
        "the on-target crossing pools are empty in this fixture, so their exclusion is untested"
    )


def _raw_model(counts: dict, max_size: int = 100):
    from rigel.frag_length_model import FragmentLengthModel

    hist = np.zeros(max_size + 1, dtype=np.float64)
    for length, n in counts.items():
        hist[length] = n
    return FragmentLengthModel(max_size=max_size, counts=hist)


def test_a_length_summary_reports_the_censored_overflow_apart():
    """The censored overflow bin is surfaced as a QC diagnostic, and kept out of the mean and the mode."""
    summary = cli.length_summary(_raw_model({50: 90.0, 100: 10.0}))
    assert summary["mode"] == 50
    assert summary["mean"] == pytest.approx(50.0, abs=0.01)
    assert summary["overflow_count"] == pytest.approx(10.0)
    assert summary["overflow_fraction"] == pytest.approx(0.1)
    assert (summary["n_observations"], summary["max_size"]) == (100, 100)


def test_a_length_histogram_is_trimmed_to_its_occupied_bins():
    lengths, counts = cli.length_histogram(_raw_model({50: 1.0, 60: 1.0}))
    assert lengths.tolist() == list(range(50, 61))  # 50 through 60 inclusive
    assert counts.sum() == 2.0
    empty_lengths, empty_counts = cli.length_histogram(_raw_model({}))
    assert empty_lengths.size == empty_counts.size == 0


# ---------------------------------------------------------------------------
# Integration: the full quant summary.json + fragment_lengths.feather
# ---------------------------------------------------------------------------


def _build_scenario(work_dir: Path):
    sc = Scenario("summary_report", genome_length=12000, seed=SEED, work_dir=work_dir)
    sc.add_gene(
        "g1",
        "+",
        [
            {"t_id": "t1", "exons": [(500, 900), (1500, 1900), (2500, 2900)], "abundance": 80},
            {"t_id": "t2", "exons": [(500, 900), (2500, 2900)], "abundance": 40},
        ],
    )
    sc.add_gene(
        "g2",
        "-",
        [
            {"t_id": "t3", "exons": [(6000, 6400), (7000, 7400)], "abundance": 50},
        ],
    )
    sim = ReadSimConfig(
        frag_mean=200,
        frag_std=40,
        frag_min=80,
        frag_max=450,
        read_length=100,
        strand_specificity=0.9,
        seed=SEED,
    )
    return sc, sc.build_oracle(
        n_rna_fragments=2000,
        sim_config=sim,
        gdna_fraction=0.2,
        gdna_config=GDNAConfig(),
        nrna_abundance=15.0,
    )


def _quant_args(result, out: Path, *flags: str) -> argparse.Namespace:
    """The arguments ``rigel quant`` hands its writer, resolved as ``quant_command`` resolves them."""
    args = cli.build_parser().parse_args(
        ["quant", "--bam", str(result.bam_path), "--index", str(result.index_dir), "-o", str(out)]
        + list(flags)
    )
    cli.resolve_quant_args(args, cli.build_quant_defaults())
    return args


@pytest.fixture(scope="module")
def quant_run(tmp_path_factory):
    """One ``rigel quant`` run's outputs, as its writer writes them: ``(pr, result, args, out, summary)``."""
    work = tmp_path_factory.mktemp("summary_report")
    out = work / "out"
    out.mkdir()
    sc, result = _build_scenario(work / "work")
    cfg = PipelineConfig(
        em=EMConfig(seed=SEED, assignment_mode="fractional", n_threads=1),
        scan=BamScanConfig(sj_strand_tag="auto"),
    )
    pr = run_pipeline(result.bam_path, result.index, config=cfg)
    args = _quant_args(result, out)
    cli.write_quant_outputs(pr, result.index, out, args)
    yield pr, result, args, out, json.loads((out / "summary.json").read_text())
    sc.cleanup()


def _manual_block(name: str, indent: str = r"\s+") -> set[str]:
    """The keys `docs/MANUAL.md` lists for the summary.json block *name*."""
    manual = (Path(__file__).resolve().parents[1] / "docs" / "MANUAL.md").read_text()
    listed = re.search(rf'"{name}": \{{\n(.*?)\n  \}}', manual, re.S)
    assert listed, f"docs/MANUAL.md lists no summary.json {name} block"
    return set(re.findall(rf'^{indent}"(\w+)":', listed.group(1), re.M))


def test_the_schema_version_is_one_constant_the_writer_and_the_report_share(quant_run, tmp_path):
    """The report targets the version the writer writes, so a bump of the writer's moves the report's with it: the
    bumped summary reads clean and the older one warns. A copy restated in the report would stay behind and warn on
    the bumped one."""
    from rigel.report.substrate import load_substrate

    pr, result, args, out, summary = quant_run

    def schema_warnings(run_dir):
        return [w for w in load_substrate(run_dir).warnings if "schema_version" in w]

    assert summary["schema_version"] == cli.SUMMARY_SCHEMA_VERSION == 9
    assert schema_warnings(out) == []
    bumped = tmp_path / "bumped"
    bumped.mkdir()
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(cli, "SUMMARY_SCHEMA_VERSION", cli.SUMMARY_SCHEMA_VERSION + 1)
        cli.write_quant_outputs(pr, result.index, bumped, args)
        assert json.loads((bumped / "summary.json").read_text())["schema_version"] == 10
        assert schema_warnings(bumped) == [] and schema_warnings(out) != []


def test_the_mismatch_rate_is_the_scans_counts_and_their_ratio(quant_run):
    """The library's mismatch rate, read on its unique unspliced hits."""
    pr, _result, _args, _out, summary = quant_run
    rate = summary["fragment_stats"]["mismatch_rate"]
    assert rate["mismatches"] == pr.stats.n_unique_unspliced_mismatches
    assert rate["bases"] == pr.stats.n_unique_unspliced_bases > 0
    assert rate["rate"] == round(rate["mismatches"] / rate["bases"], cli.SUMMARY_DECIMALS)


def test_every_read_pair_has_one_outcome_and_every_record_one_status(quant_run):
    """Each partition sums to its total as the scan counted it; the transcripts are the annotated ones, the synthetic
    spans counted apart as the nascent RNA, and the default alpha gates nothing."""
    pr, result, _args, _out, summary = quant_run
    pairs, records = summary["read_pairs"], summary["alignment_stats"]["records"]
    placed = pairs["placed"]["unique"] + pairs["placed"]["multimapping"]
    assert pairs["total"] == pr.stats.n_read_pairs > 0
    assert pairs["total"] == (
        pairs["excluded"] + sum(pairs["withheld"].values()) + pairs["excluded_by_option"] + placed
    )
    assert sum(pairs["unpaired_by_reason"].values()) == pairs["withheld"]["unpaired"]
    assert placed == summary["fragment_stats"]["total"]
    assert sum(records.values()) == summary["alignment_stats"]["total_reads"] == pr.stats.total
    is_synthetic = result.index.t_df["is_synthetic"].to_numpy(dtype=bool)
    assert summary["quantification"]["n_nascent"] == int(is_synthetic.sum()) > 0
    assert summary["quantification"]["n_transcripts"] == int((~is_synthetic).sum())
    assert summary["quantification"]["n_gated_out"] == 0


def test_the_length_block_is_stats_only_and_the_companion_table_reconciles(quant_run):
    """``fragment_length`` carries no per-bin histogram; ``fragment_lengths.feather`` does, and each category's
    counts sum to the summary's observations."""
    _pr, _result, _args, out, summary = quant_run
    fl = summary["fragment_length"]
    assert "global" in fl
    for d in fl.values():
        assert "histogram" not in d
    flf = pd.read_feather(out / "fragment_lengths.feather")
    assert list(flf.columns) == ["category", "length", "count"]
    for cat in flf["category"].unique():
        fsum = float(flf.loc[flf["category"] == cat, "count"].sum())
        assert abs(fsum - fl[cat]["n_observations"]) < 0.5


def test_the_splice_breakdown_carries_the_blacklist_provenance(quant_run):
    """Implicit, artifact and blacklisted counts beside the blacklist's size and whether one was loaded; this
    fixture's index has none, so detection is off."""
    splice = quant_run[4]["fragment_stats"]["splice"]
    assert set(splice) == {
        "unspliced",
        "spliced_annotated",
        "spliced_implicit",
        "splice_artifact",
        "sj_blacklisted",
        "sj_blacklist_size",
        "sj_blacklist_loaded",
    }
    assert splice["sj_blacklist_loaded"] is False
    assert splice["sj_blacklist_size"] == 0


def test_the_sj_strand_tag_counts_are_the_scans(quant_run):
    """The junction-strand tag, read on every spliced record: what each record's tags gave."""
    pr, _result, _args, _out, summary = quant_run
    tags = summary["alignment_stats"]["sj_strand_tag"]
    assert tags == {
        "tagged": pr.stats.n_sj_tag_tagged,
        "untagged": pr.stats.n_sj_tag_untagged,
        "conflicting": pr.stats.n_sj_tag_conflicting,
        "invalid": pr.stats.n_sj_tag_invalid,
    }
    assert tags["tagged"] > 0


def test_the_strand_block_carries_the_sj_table_and_the_contamination_diagnostic(quant_run):
    """The per-sj table's total is the 2x2's total: one observation per qualified fragment."""
    strand = quant_run[4]["strand_model"]
    assert {"n_sj", "n_observations", "depth_median", "depth_max"} <= set(strand["sj"])
    assert strand["sj"]["n_observations"] == strand["n_training_fragments"]
    assert {
        "exonic_all_specificity",
        "exonic_all_p_r1_sense",
        "exonic_all_n_fragments",
        "contamination_gap",
    } <= set(strand["diagnostics"])


def test_the_calibration_track_is_written_as_feather_and_bedgraph(quant_run):
    pr, _result, _args, out, _summary = quant_run
    assert pr.calibration_track is not None
    track = pd.read_feather(out / "calibration_track.feather")
    assert list(track.columns) == [
        "ref",
        "start",
        "end",
        "count_gdna",
        "count_rna",
        "gdna_density",
        "fraction_gdna",
    ]
    assert len(track) == pr.calibration.n_regions
    assert (track["fraction_gdna"] >= 0).all() and (track["fraction_gdna"] <= 1).all()
    bg = (out / "calibration_track.bedgraph").read_text().splitlines()
    assert bg[0].startswith("track type=bedGraph")
    assert len(bg) == len(track) + 1  # header + one line per region


def test_the_calibration_block_is_the_library_scalars_the_manual_lists(quant_run):
    """Calibration's capture answer is the bulk density it reads capture against (null when no landscape was fitted),
    the weight its reading is raised to, the saturated level over the bulk at that weight and the junction cap fitted
    above that level; the user manual's listing of the block is the same set."""
    written = set(quant_run[4]["calibration"])
    assert written == {
        "gdna_density_global",
        "gdna_bulk_density",
        "gdna_capture_weight",
        "gdna_capture_saturation",
        "junction_capture_cap",
        "fraction_rna_sense",
        "n_regions",
        "n_boundaries",
        "n_sj",
    }
    assert _manual_block("calibration") == written


def test_the_observable_genome_block_is_the_indexs_record_the_manual_lists(quant_run):
    """The index's record of its mask over every reference's bases."""
    _pr, result, _args, _out, summary = quant_run
    observable = summary["observable_genome"]
    assert observable == cli.observable_genome_summary(result.index)
    assert observable["observable_bases"] + observable["unobservable_bases"] == sum(
        result.index.ref_lengths.values()
    )
    assert _manual_block("observable_genome", indent="    ") == set(observable)


def test_the_length_tables_block_is_the_pipelines_decision_the_manual_lists(quant_run):
    """The EM's two length tables: the decision the pipeline made, as written."""
    pr, _result, _args, _out, summary = quant_run
    tables = summary["length_tables"]
    d = pr.length_tables.decision
    assert (
        tables["weight"] == round(d.weight, cli.SUMMARY_DECIMALS) and 0.0 <= tables["weight"] <= 1.0
    )
    assert tables["kappa"] == round(d.kappa, cli.SUMMARY_DECIMALS)
    assert tables["gdna_table_mean"] == round(
        float(pr.length_tables.gdna @ np.arange(pr.length_tables.gdna.size)), 2
    )
    assert _manual_block("length_tables") == set(tables)


def test_every_intergenic_fragment_is_gdna():
    """An intergenic fragment is unspliced, since a spliced hit the scan keeps has a transcript holding its junctions,
    so the gDNA count takes every one, and the three pools' fractions sum to one."""
    from types import SimpleNamespace

    q = cli.quantification_pools(
        SimpleNamespace(n_intergenic=100),
        count_annotated=600.0,
        count_nascent=100.0,
        count_gdna_em=200.0,
        count_gdna_intergenic_em=0.0,
    )
    assert set(q) == {
        "count_annotated",
        "count_nascent",
        "count_gdna",
        "count_gdna_em",
        "fraction_annotated",
        "fraction_nascent",
        "fraction_gdna",
        "fraction_gdna_em",
        "count_intergenic",
        "fraction_intergenic",
    }
    assert (q["count_gdna"], q["count_gdna_em"]) == (300.0, 200.0)
    assert q["fraction_gdna"] == round(300 / 1000, cli.SUMMARY_DECIMALS)
    assert q["fraction_annotated"] + q["fraction_nascent"] + q["fraction_gdna"] == pytest.approx(
        1.0
    )


def test_the_strand_block_carries_the_junction_diagnostics_and_one_overdispersion(tmp_path):
    """The strand model's junction fit is surfaced: κ pooled over every junction and over the genuine ones, each
    class's share of the junctions and of their fragments, and whether the strand channel is live. The two
    overdispersions, both 0 by policy, are one field that says the channel is binomial."""
    out = tmp_path / "out"
    out.mkdir()
    sc, result = _build_scenario(tmp_path / "work")
    cfg = PipelineConfig(
        em=EMConfig(seed=SEED, assignment_mode="fractional", n_threads=1),
        scan=BamScanConfig(sj_strand_tag="auto"),
    )
    pr = run_pipeline(result.bam_path, result.index, config=cfg)
    args = _quant_args(result, out)
    cli.write_quant_outputs(pr, result.index, out, args)
    summary = json.loads((out / "summary.json").read_text())
    strand = summary["strand_model"]
    assert set(strand["kappa"]) == {"pooled", "genuine"}
    assert set(strand["junction_classes"]) == {"by_junction", "by_fragment"}
    for shares in strand["junction_classes"].values():
        assert set(shares) == {"genuine", "artifact", "reversed"}
        assert sum(shares.values()) == pytest.approx(1.0, abs=1e-6)
    assert strand["strand_live"] is True  # the scenario is stranded at 0.9
    assert strand["strand_overdispersion"] == {"value": 0.0, "model": "binomial"}
    # calibration's density at full precision: rounded to 8 places, a human-genome density below ~15 fragments reads 0
    assert (
        summary["calibration"]["gdna_density_global"]
        == float(pr.calibration.gdna_density_global)
        > 0
    )
    assert not {"gdna_strand_overdispersion", "rna_strand_overdispersion"} & set(
        summary["calibration"]
    )


def test_a_dead_strand_channel_publishes_the_kappa_calibration_reads_and_no_classes():
    """Where the strand channel is dead the three classes are not identified — at κ = ½ an artifact junction looks like
    any other — so the fit's genuine κ and class shares are noise, and calibration reads the pooled counts instead
    (``StrandModel.genuine_n_same``). The block must say what calibration used: genuine equals pooled, no classes."""
    from rigel.strand_model import SJStrandTable, StrandModel, StrandModels

    rng = np.random.default_rng(
        5
    )  # a dead channel whose fit reads 98 % artifact and genuine κ 0.82 from noise
    depth = rng.integers(2, 12, size=60)
    sense = rng.binomial(depth, 0.5)
    table = SJStrandTable(
        ref_id=np.zeros(depth.size, dtype=np.int32),
        start=np.arange(depth.size, dtype=np.int64) * 1000,
        end=np.arange(depth.size, dtype=np.int64) * 1000 + 100,
        motif_strand=np.ones(depth.size, dtype=np.int8),
        n_sense=sense.astype(np.int64),
        n_antisense=(depth - sense).astype(np.int64),
    )
    block = cli.strand_fit_summary(StrandModels(exonic_spliced=StrandModel.from_sj_table(table)))
    assert block["strand_live"] is False
    assert block["kappa"]["genuine"] == block["kappa"]["pooled"]
    assert block["junction_classes"] == {"by_junction": None, "by_fragment": None}


def test_the_locus_statistics_are_one_row_per_solved_locus_in_locus_order(tmp_path):
    """Ordered by locus, never by wall time, so two runs write the same rows in the same order. The second locus is
    the heavier, so an order by time would put it first."""
    sc = Scenario("locus_order", genome_length=12000, seed=SEED, work_dir=tmp_path / "work")
    sc.add_gene("light", "+", [{"t_id": "t1", "exons": [(500, 900), (1500, 1900)], "abundance": 2}])
    sc.add_gene(
        "heavy",
        "+",
        [
            {"t_id": "t2", "exons": [(6000, 6400), (7000, 7400), (8000, 8400)], "abundance": 400},
            {"t_id": "t3", "exons": [(6000, 6400), (8000, 8400)], "abundance": 300},
            {"t_id": "t4", "exons": [(6200, 6400), (7000, 7400), (8000, 8300)], "abundance": 200},
        ],
    )
    sim = ReadSimConfig(frag_mean=200, frag_std=40, read_length=100, seed=SEED)
    result = sc.build_oracle(
        n_rna_fragments=20000, sim_config=sim, gdna_fraction=0.2, gdna_config=GDNAConfig()
    )
    cfg = PipelineConfig(
        em=EMConfig(seed=SEED, assignment_mode="fractional", n_threads=1),
        scan=BamScanConfig(sj_strand_tag="auto"),
        emit_locus_stats=True,
    )
    stats = run_pipeline(result.bam_path, result.index, config=cfg).estimator.locus_stats
    sc.cleanup()
    assert [s["locus_idx"] for s in stats] == [0, 1]


def test_a_reused_output_directory_keeps_no_optional_file_this_run_did_not_write(tmp_path):
    """The TSV mirrors and the locus statistics are written only on request, so a rerun into the same directory
    without them removes the earlier run's copies: a stale quant.tsv would disagree with the quant.feather beside it."""
    sc, result = _build_scenario(tmp_path / "work")
    cfg = PipelineConfig(
        em=EMConfig(seed=SEED, assignment_mode="fractional", n_threads=1),
        scan=BamScanConfig(sj_strand_tag="auto"),
        emit_locus_stats=True,
    )
    pr = run_pipeline(result.bam_path, result.index, config=cfg)
    out = tmp_path / "out"
    out.mkdir()

    cli.write_quant_outputs(pr, result.index, out, _quant_args(result, out, "--tsv"))
    assert (out / "quant.tsv").exists() and (out / "locus_stats.feather").exists()
    pr.estimator.locus_stats = None  # a run that collected no locus stats
    cli.write_quant_outputs(pr, result.index, out, _quant_args(result, out))
    assert not sorted(out.glob("*.tsv")) and not (out / "locus_stats.feather").exists()


def test_the_nascent_table_writes_exactly_the_columns_the_manual_lists(tmp_path):
    """``nascent_quant.feather`` is ``get_nascent_quant_df``: one row per synthetic span, and the user manual's table
    names exactly its columns, in the file's order, so the two cannot drift apart."""
    sc, result = _build_scenario(tmp_path / "work")
    cfg = PipelineConfig(
        em=EMConfig(seed=SEED, assignment_mode="fractional", n_threads=1),
        scan=BamScanConfig(sj_strand_tag="auto"),
    )
    pr = run_pipeline(result.bam_path, result.index, config=cfg)
    nascent = pr.estimator.get_nascent_quant_df(result.index)
    t_df = result.index.t_df
    assert sorted(nascent["nascent_id"]) == sorted(t_df.loc[t_df["is_synthetic"], "t_id"])
    assert len(nascent) > 0 and (nascent["n_annotated"] >= 1).all()
    manual = (Path(__file__).resolve().parents[1] / "docs" / "MANUAL.md").read_text()
    section = re.search(
        r"### nascent_quant\.feather / nascent_quant\.tsv\n(.*?)\n### ", manual, re.S
    )
    assert section, "docs/MANUAL.md has no nascent_quant section"
    rows = re.findall(r"^\| (`[^|]*`) \|", section.group(1), re.M)
    listed = [name for cell in rows for name in re.findall(r"`(\w+)`", cell)]
    assert listed == list(nascent.columns), "the MANUAL lists the columns in the file's order"
    sc.cleanup()
