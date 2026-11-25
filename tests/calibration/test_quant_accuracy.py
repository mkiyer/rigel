"""Tests for ``scripts/design/quant_accuracy.py``: the end-to-end accuracy table and the prior-injection ceiling above
it.

The instrument answers two questions with one scorer: how accurate the tool is, and how much of that error a perfect
calibration prior would remove. Both answers are differences between an estimate and a truth column, so every way of
getting them wrong produces a plausible number: the wrong truth column, an injection that never fired, a false-positive
channel filtered out before scoring. Each test covers one of those.
"""

from __future__ import annotations

import dataclasses
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from rigel.config import EMConfig, PipelineConfig
import _shared
import quant_accuracy as qa
from _prior_toy import build_toy


# ── the toy ──────────────────────────────────────────────────────────────────────────────────────


@pytest.fixture(scope="module")
def toy(tmp_path_factory):
    """A contaminated toy: three genes, staggered isoforms, gDNA at 1x the RNA fragment count."""
    return build_toy(tmp_path_factory)


@pytest.fixture(scope="module")
def toy_oracle(toy, tmp_path_factory):
    from _oracle import OracleTruth

    return OracleTruth.from_bam(
        str(toy.bam_path),
        toy.index,
        PipelineConfig(),
        tmp_path_factory.mktemp("qa_split"),
        "qa",
    )


def _run(toy, arm, oracle, seed=EMConfig().seed):
    """One full ``run_pipeline`` under one arm; returns ``(PipelineResult, fired)``. The seed is the instrument's, the
    shipped one: under the default ``assignment_mode="sample"`` the hard assignment is a draw, so byte-identity needs
    one seed.
    """
    from rigel.pipeline import run_pipeline

    cfg = PipelineConfig()
    cfg = dataclasses.replace(cfg, em=dataclasses.replace(cfg.em, seed=seed))
    restore, fired = qa.install_arm(arm, oracle)
    try:
        result = run_pipeline(str(toy.bam_path), toy.index, cfg)
    finally:
        restore()
    return result, fired


def _quant(toy, arm, oracle, seed=EMConfig().seed):
    """:func:`_run`'s transcript table; returns ``(counts_df, fired)``."""
    result, fired = _run(toy, arm, oracle, seed)
    return result.estimator.get_counts_df(toy.index), fired


# ── Gate 1: noop is byte-identical end to end, and the wrapper really ran ────────────────────────


def test_the_noop_arm_reproduces_base_byte_identically_through_the_whole_pipeline(toy, toy_oracle):
    """The injection reaches the EM, whose per-locus draw is order-sensitive, so a wrapper that perturbed anything (even
    a summation order) would add that perturbation to every ``oracle`` number. ``noop`` is not an early return: it
    builds O in full and discards it, so this covers the whole wrapper path. The real ``oracle`` arm must then move the
    EM's gDNA count and break the match, showing the comparison can fail.
    """
    base_result, base_fired = _run(toy, "base", None)
    noop_result, noop_fired = _run(toy, "noop", toy_oracle)
    assert base_fired["n"] >= 1 and noop_fired["n"] >= 1
    base = base_result.estimator.get_counts_df(toy.index)
    noop = noop_result.estimator.get_counts_df(toy.index)

    # Two standards. Every count column is integer-derived and must match exactly, since integer addition is
    # associative. `posterior_mean` is float-derived and varies by a few ulp between identical runs, because the
    # conserved-mass banks are float64 and the per-worker merge order varies. The tolerance sits far above that spread
    # (with headroom for deeper libraries, since re-association error grows with the number of additions) and far below
    # the injection it must catch, and the oracle arm below shows it is tight enough.
    exact = [c for c in base.columns if c != "posterior_mean"]

    def compare(other):
        pd.testing.assert_frame_equal(base[exact], other[exact], check_exact=True)
        pd.testing.assert_series_equal(
            base["posterior_mean"], other["posterior_mean"], check_exact=False, rtol=1e-9, atol=0.0
        )

    compare(noop)
    # The EM's gDNA count per locus, the component the injected field names: integer-derived, so exact.
    base_gdna = base_result.estimator.get_loci_df(toy.index)["count_gdna"]
    pd.testing.assert_series_equal(
        base_gdna, noop_result.estimator.get_loci_df(toy.index)["count_gdna"], check_exact=True
    )

    perturbed_result, _ = _run(toy, "oracle", toy_oracle)
    # The witness is the EM's gDNA count, not the transcript `count` column: on this toy the mass the oracle's gDNA
    # count moves lands between gDNA and nascent RNA, while the annotated transcripts move by under one fragment, which
    # one seeded draw need not flip.
    assert not perturbed_result.estimator.get_loci_df(toy.index)["count_gdna"].equals(base_gdna), (
        "taking the oracle's gDNA count changed no locus's gDNA count on a contaminated toy — the "
        "injection does not reach the EM and every ceiling this script reports would be zero"
    )
    with pytest.raises(AssertionError):
        compare(perturbed_result.estimator.get_counts_df(toy.index))


def test_base_is_the_shipped_config_and_a_pinned_arm_reproduces_itself(toy):
    """``base`` must be the shipped configuration, so the instrument's seed is the shipped one, and two runs of one arm
    must return one answer, or no arm comparison is attributable.
    """
    a, _ = _quant(toy, "base", None)
    b, _ = _quant(toy, "base", None)
    assert np.array_equal(a["count"].to_numpy(), b["count"].to_numpy()), (
        "two runs of one arm disagree — no arm comparison here is attributable"
    )


def test_an_injection_that_never_fires_raises_rather_than_reporting_no_effect(
    tmp_path, monkeypatch
):
    """ "The perfect prior changed nothing" is a publishable conclusion, so a wrapper that was never called must not be
    able to produce it; that happens when ``quant_from_buffer`` returns early on an empty unit set. ``run_pipeline`` is
    stubbed to a no-op, the shape of that early return: the pipeline completes, the wrapper is installed, and
    ``assemble_priors`` is never reached.
    """
    cond = tmp_path / "c1"
    cond.mkdir()
    pd.DataFrame({"transcript_id": ["a"], "observed_annotated_fragments": [1]}).to_csv(
        cond / "truth_abundances.tsv", sep="\t", index=False
    )
    (cond / "truth_summary.json").write_text(
        '{"origin_counts": {"gdna": 0, "annotated": 1, "nascent": 0}}'
    )
    monkeypatch.setattr(qa, "run_pipeline", lambda *a, **k: None)
    monkeypatch.setattr(qa, "load_oracle", lambda *a, **k: object())
    with pytest.raises(RuntimeError, match="never wrapped-and-called"):
        qa.run_condition("oracle", tmp_path, None, "c1", PipelineConfig(), tmp_path)


# ── Gate 3: the truth column is the realised fragment count, not the molar abundance ─────────────


def test_the_scorer_reads_the_realised_fragment_column_and_refuses_the_molar_table(tmp_path):
    """The panel ships two truth tables holding different quantities. A condition's ``truth_abundances.tsv`` gives the
    realised observed fragment count (``observed_annotated_fragments``) beside the pre-capture molar abundance
    (``pre_capture_annotated_abundance``); the suite-level ``truth_abundances_nascent_<label>.tsv`` gives the molar abundance
    alone, as ``annotated_abundance``. Scoring a fragment estimate against a molar truth would charge the tool for hybrid
    capture, which it never claims to invert.

    Both tables are written by the simulator's own writers and read as the instrument reads them. The scorer must score
    the realised count exactly, refuse the molar table, and see a molar truth placed in the realised column as error.
    """
    from rigel.transcript import Transcript

    # the writers are the whole-genome grid's, a lab script under `scripts/sim/`
    from _sim_truth import write_molar_truth, write_post_capture_truth
    from rigel.types import Interval, Strand

    transcripts = []
    for t_id, molecules in (("a", 40.0), ("b", 10.0), ("c", 5.0)):
        transcript = Transcript(
            ref="chr1",
            strand=Strand.POS,
            exons=[Interval(0, 500)],
            t_id=t_id,
            g_id=f"g_{t_id}",
            abundance=molecules,
        )
        transcript.compute_length()
        transcripts.append(transcript)
    # `c` has molecules and was never captured, so it realised no fragment
    reads = ["a"] * 1000 + ["b"] * 500
    fastq = tmp_path / "reads.fq"
    fastq.write_text("".join(f"@{t}:0-100:f:{i}/1\nACGT\n+\nIIII\n" for i, t in enumerate(reads)))
    cond = tmp_path / "cond"
    write_post_capture_truth(
        transcripts,
        cond / "truth_abundances.tsv",
        cond / "truth_fragment_lengths.tsv",
        cond / "truth_summary.json",
        fastq_path=fastq,
    )
    write_molar_truth(transcripts, tmp_path / "truth_abundances_nascent_none.tsv")
    realised = pd.read_csv(cond / "truth_abundances.tsv", sep="\t")
    molar = pd.read_csv(tmp_path / "truth_abundances_nascent_none.tsv", sep="\t")

    quant = pd.DataFrame(
        {
            "transcript_id": ["a", "b", "c"],
            "count": [1000.0, 500.0, 0.0],
            "count_em": [0.0, 0.0, 0.0],
            "tpm": [0.0, 0.0, 0.0],
            "effective_length": [100.0, 100.0, 100.0],
        }
    )
    exact = qa.score_transcripts(quant, realised)
    assert exact["count_abs_err"] == pytest.approx(0.0)
    assert exact["fp_mass"] == pytest.approx(0.0)

    with pytest.raises(KeyError, match="observed_annotated_fragments"):
        qa.score_transcripts(quant, molar)

    wrong = qa.score_transcripts(
        quant,
        realised.assign(observed_annotated_fragments=realised["pre_capture_annotated_abundance"]),
    )
    assert wrong["count_abs_err"] > 1000.0, (
        "a molar truth scored as a fragment count produced almost no error — the two units are not "
        "distinguishable to this scorer and the wrong file could be read silently"
    )


# ── Gate 4: the two false-positive channels are separated, and neither hides ─────────────────────


def test_false_positive_and_false_negative_mass_are_separate_and_reconcile():
    """``fp_mass`` (estimate on a transcript the simulator gave zero fragments) and ``fn_mass`` (real RNA assigned
    nowhere) are different failures with different fixes, which a single ``Σ|Δ|`` cannot separate. Flipping one
    transcript from silent to expressed must move its mass out of ``fp_mass``; a scorer keyed on the estimate rather
    than the truth would keep it there.
    """
    truth = pd.DataFrame(
        {"transcript_id": ["a", "b", "c"], "observed_annotated_fragments": [100.0, 0.0, 50.0]}
    )
    quant = pd.DataFrame(
        {
            "transcript_id": ["a", "b", "c"],
            "count": [120.0, 30.0, 0.0],
            "count_em": [0.0, 0.0, 0.0],
            "tpm": [0.0, 0.0, 0.0],
            "effective_length": [100.0, 100.0, 100.0],
        }
    )
    s = qa.score_transcripts(quant, truth)
    assert s["fp_mass"] == pytest.approx(30.0) and s["fp_n"] == 1
    assert s["fn_mass"] == pytest.approx(50.0) and s["fn_n"] == 1
    assert s["count_over"] - s["count_under"] == pytest.approx(s["count_net_err"])
    assert s["count_over"] + s["count_under"] == pytest.approx(s["count_abs_err"])

    lifted = truth.assign(observed_annotated_fragments=[100.0, 25.0, 50.0])
    s2 = qa.score_transcripts(quant, lifted)
    assert s2["fp_mass"] == pytest.approx(0.0), (
        "a transcript the truth now expresses is still counted as a false positive — fp_mass is "
        "keyed on the estimate, not on the truth"
    )


def test_the_gene_row_collapses_isoform_ambiguity_and_sees_real_gene_error():
    """The gene row answers "is the remaining transcript error isoform ambiguity, or something else?", which it can only
    do if moving mass between a gene's own isoforms is invisible to it while moving mass out of the gene is not. Both
    halves are asserted, since a scorer that saw neither would also report zero on the first.
    """
    truth = pd.DataFrame(
        {
            "transcript_id": ["a1", "a2", "b1"],
            "gene_id": ["A", "A", "B"],
            "observed_annotated_fragments": [100.0, 0.0, 50.0],
        }
    )
    swapped = pd.DataFrame(
        {"transcript_id": ["a1", "a2", "b1"], "count": [0.0, 100.0, 50.0]}
    )  # the isoform call is completely wrong; the gene call is perfect
    tx = qa.score_transcripts(swapped.assign(count_em=0.0, tpm=0.0, effective_length=100.0), truth)
    gene = qa.score_genes(swapped, truth)
    assert tx["count_abs_err"] == pytest.approx(200.0)
    assert gene["count_abs_err"] == pytest.approx(0.0), (
        "a pure isoform swap inside one gene shows up at gene level — the grouping is not collapsing "
        "isoform ambiguity and the gene row answers nothing the transcript row did not"
    )

    leaked = pd.DataFrame({"transcript_id": ["a1", "a2", "b1"], "count": [60.0, 0.0, 90.0]})
    gene2 = qa.score_genes(leaked, truth)
    assert gene2["count_abs_err"] == pytest.approx(80.0), (
        "40 fragments moved from gene A to gene B and the gene row did not notice — it is summing "
        "away real error, not ambiguity"
    )


def test_the_library_gdna_fraction_counts_intergenic_and_the_em_only_view_does_not():
    """``n_intergenic`` counts fragments that reached no locus; they never enter the EM, so ``gdna_em_count`` excludes
    them. Off capture gDNA is genome-uniform and more than half of it is intergenic, so a library gDNA fraction that
    drops them from the numerator while RNA stays in the denominator reads far low against truth. The deliverable
    includes the intergenic fragments, the EM-only diagnostic does not, and the two must differ.
    """
    annotated, gdna_em, intergenic = 5_000_000.0, 2_330_992.0, 2_601_271.0
    quant = pd.DataFrame({"transcript_id": ["a"], "count": [annotated]})
    result = _FakeResult(synthetic_count=0.0, gdna_em_count=gdna_em, n_intergenic=intergenic)
    summary = {"origin_counts": {"gdna": 5_000_000, "annotated": 5_000_000, "nascent": 0}}
    lib = qa.score_library(result, quant, summary)

    # The expectations are derived from the fixture's own inputs, never pasted from a run.
    want_deliverable = (gdna_em + intergenic) / (annotated + gdna_em + intergenic)
    want_em_only = gdna_em / (annotated + gdna_em)
    assert lib["gdna_frac_true"] == pytest.approx(0.5)
    assert lib["gdna_frac_est"] == pytest.approx(want_deliverable), (
        "the deliverable fraction is not counting intergenic fragments as gDNA"
    )
    assert lib["gdna_em_frac_est"] == pytest.approx(want_em_only)
    # the EM-only view reads far low against a truth of 0.5
    assert want_em_only < 0.33 and want_deliverable > 0.49
    assert lib["gdna_frac_est"] != pytest.approx(lib["gdna_em_frac_est"]), (
        "the two views are equal, so one of them is not being computed — the whole point is that "
        "they differ by exactly the population the EM never sees"
    )

    # with no intergenic fragments the two views must coincide exactly
    same = qa.score_library(
        _FakeResult(synthetic_count=0.0, gdna_em_count=2_330_992.0, n_intergenic=0.0),
        quant,
        summary,
    )
    assert same["gdna_frac_est"] == pytest.approx(same["gdna_em_frac_est"]), (
        "with zero intergenic fragments the two denominators are identical and the two numbers must "
        "be too — they are not, so one of them has a second difference in it"
    )


@dataclasses.dataclass
class _FakeEstimator:
    synthetic_count: float
    gdna_em_count: float
    gdna_intergenic_em_count: float = 0.0


@dataclasses.dataclass
class _FakeStats:
    n_intergenic: float


class _FakeResult:
    """Enough of a ``PipelineResult`` for ``score_library``: the four numbers it reads."""

    def __init__(self, synthetic_count, gdna_em_count, n_intergenic):
        self.estimator = _FakeEstimator(synthetic_count, gdna_em_count)
        self.stats = _FakeStats(n_intergenic)


def test_every_intergenic_fragment_is_gdna_in_the_library_split():
    """An intergenic fragment is unspliced (the scan withholds a spliced hit no transcript holds), so every one is
    gDNA, as in ``cli.py``'s ``fraction_gdna``."""
    annotated, gdna_em, intergenic = 5_000.0, 2_000.0, 3_000.0
    quant = pd.DataFrame({"transcript_id": ["a"], "count": [annotated]})
    summary = {"origin_counts": {"gdna": 4_900, "annotated": 5_100, "nascent": 0}}
    lib = qa.score_library(_FakeResult(0.0, gdna_em, n_intergenic=intergenic), quant, summary)
    assert lib["gdna_frac_est"] == pytest.approx(
        (gdna_em + intergenic) / (annotated + gdna_em + intergenic)
    )
    assert "n_intergenic_rna" not in lib and lib["n_intergenic_gdna"] == intergenic
    assert lib["rna_est"] == annotated


def test_the_count_metrics_do_not_move_when_the_length_model_does():
    """The count rows must be free of the fragment-length model, or the instrument would price the length models too and
    neither number would be attributable. Changing one transcript's effective length by 100x must not move them.
    """
    truth = pd.DataFrame(
        {"transcript_id": ["a", "b"], "observed_annotated_fragments": [100.0, 50.0]}
    )
    quant = pd.DataFrame(
        {
            "transcript_id": ["a", "b"],
            "count": [120.0, 30.0],
            "count_em": [0.0, 0.0],
            "tpm": [500000.0, 500000.0],
            "effective_length": [100.0, 100.0],
        }
    )
    a = qa.score_transcripts(quant, truth)
    b = qa.score_transcripts(quant.assign(effective_length=[10000.0, 100.0]), truth)
    for field in ("count_abs_err", "count_net_err", "fp_mass", "fn_mass", "spearman"):
        assert a[field] == pytest.approx(b[field], nan_ok=True), field
    assert a["tpm_abs_err"] != pytest.approx(b["tpm_abs_err"]), (
        "the TPM row did not move when an effective length moved 100x — it is not reading the "
        "length model at all and its caption is wrong"
    )


# ── Gate 5: a missing condition is a refusal, never a quiet average ──────────────────────────────


def test_the_report_refuses_arms_with_different_row_sets(tmp_path):
    """A shard that died would leave surviving rows that aggregate into a stratum total reading like a complete panel,
    so arms with different row sets are refused."""
    import json

    full = tmp_path / "full.jsonl"
    short = tmp_path / "short.jsonl"
    rows = [
        {
            "arm": "a",
            "condition": f"gdna_g{g}_ss_0.99_nascent_none_capture_on",
            "axis": ax,
            "count_abs_err": 1.0,
            "fp_mass": 1.0,
            "fn_mass": 1.0,
            "tpm_abs_err": 1.0,
            "spearman": 1.0,
            "mard": 0.0,
            "gdna_frac_true": 0.5,
            "gdna_frac_est": 0.5,
            # the pool-level fields table ⑥ reads, part of the `library` row's real schema
            "gdna_est": 1.0,
            "gdna_true": 1.0,
            "nascent_est": 0.0,
            "nascent_true": 0.0,
            "annotated_est": 1.0,
            "annotated_true": 1.0,
            "n_intergenic_gdna": 0.0,
        }
        for g in ("01", "50")
        for ax in ("transcript", "library")
    ]
    full.write_text("\n".join(json.dumps(r) for r in rows) + "\n")
    short.write_text("\n".join(json.dumps(r) for r in rows[:2]) + "\n")
    qa.report([full])  # one arm alone is fine
    with pytest.raises(SystemExit, match="different row set"):
        qa.report([full, short])


# ── Gate 6: a missing oracle cache is a refusal, never an invented truth ─────────────────────────


def test_a_missing_or_stale_oracle_cache_aborts(toy, tmp_path):
    """The oracle arms exist to inject truth, so falling back to anything else (a rescan under a different scan config,
    an empty payload) would inject something else under the name ``oracle``. ``read_scan_cache`` refuses a payload whose
    ``reach_digest`` does not describe this index; this checks the refusal propagates. Sum-to-full is checked against
    the condition's live scan cache, so partitions of another scan are refused.
    """
    import _oracle_arms

    from rigel.pipeline import scan_and_buffer
    from _scan_cache import write_scan_cache

    bam, cfg = str(toy.bam_path), PipelineConfig()
    with pytest.raises(SystemExit, match="no valid scan or oracle cache"):
        qa.load_oracle(toy.index, cfg, tmp_path, tmp_path, "absent")

    scan = cfg.scan
    _stats, strand_model, _buf, payload = scan_and_buffer(bam, toy.index, scan)
    orc = tmp_path / "oracle_cache"
    _oracle_arms.load_or_build_oracle(bam, toy.index, cfg, tmp_path / "w", "t", payload, orc)

    def cache_scan():
        write_scan_cache(
            tmp_path / "scan_cache" / "t",
            payload=payload,
            strand_model=strand_model,
            index=toy.index,
            bam=bam,
            scan_config=scan,
        )

    cache_scan()
    qa.load_oracle(toy.index, cfg, tmp_path, orc, "t")  # the live scan's partitions load
    payload.region_contained_count[0, 0] += 1
    cache_scan()
    with pytest.raises(AssertionError, match="oracle INVALID"):
        qa.load_oracle(toy.index, cfg, tmp_path, orc, "t")


def test_an_oracle_part_records_the_bam_it_was_split_from(toy, tmp_path):
    """An oracle part is scanned from a split file in a scratch directory named for its condition, which the next build
    of a same-named condition, in any panel, overwrites. The part is a function of the BAM it was split from, so its
    cache records that BAM, and a reused split file is refused.
    """
    import _oracle_arms

    from rigel.pipeline import scan_and_buffer
    from _scan_cache import read_scan_cache

    bam, cfg = str(toy.bam_path), PipelineConfig()
    scan = cfg.scan
    _stats, _sm, _buf, payload = scan_and_buffer(bam, toy.index, scan)
    orc, work = tmp_path / "oracle_cache", tmp_path / "w"
    _oracle_arms.load_or_build_oracle(bam, toy.index, cfg, work, "t", payload, orc)
    parts = sorted(m.parent for m in orc.glob("t/*/manifest.json"))
    splits = sorted(work.glob("*.bam"))
    assert parts and splits
    for split in splits:  # what another condition's build with this tag leaves at the same paths
        split.write_bytes(b"another condition's reads")
    for part in parts:
        recorded = json.loads((part / "manifest.json").read_text())["bam"]
        assert recorded["path"] == str(Path(bam).resolve()), part.name
        read_scan_cache(part, toy.index, scan)


# ── Gate 7: the ruler arm hands the EM the lengths it names, and its noop is inert ──────────────────


def _quant_ruler(toy, arm, factor, seed=EMConfig().seed):
    """One full ``run_pipeline`` with the EM's effective lengths under ``install_ruler_arm``, the toy's one reference
    carrying gDNA, so nothing is held."""
    from rigel.pipeline import run_pipeline

    cfg = PipelineConfig()
    cfg = dataclasses.replace(cfg, em=dataclasses.replace(cfg.em, seed=seed))
    restore, fired = qa.install_ruler_arm(arm, factor, np.zeros(len(factor), dtype=bool))
    try:
        result = run_pipeline(str(toy.bam_path), toy.index, cfg)
    finally:
        restore()
    return result.estimator.get_counts_df(toy.index), fired


def test_the_ruler_arm_hands_the_em_the_lengths_it_names_and_its_noop_is_inert(toy):
    """The arm claims the EM divided by ``fl × factor``, so the test reads the EM's own published length
    (``em_effective_length``) rather than the wrapper's return value, which would stay healthy if the length were
    dropped before the solver. The factor differs between isoforms of one gene, so the arm must move the counts; the
    noop builds the same factor and hands back the shipped lengths, so it must reproduce ``base`` by the noop test's two
    standards. On this capture-OFF toy the shipped lengths are the plain ones."""
    base, _ = _quant(toy, "base", None)
    ids = toy.index.t_df["t_id"].to_numpy()
    factor = np.where(np.isin(ids, ["t1b", "t3b"]), 0.25, 1.0)

    sub, fired = _quant_ruler(toy, "oracle_ruler", factor)
    plain = base.set_index("transcript_id")["em_effective_length"]
    want = plain * pd.Series(factor, index=ids).reindex(plain.index)
    got = sub.set_index("transcript_id")["em_effective_length"].reindex(plain.index)
    np.testing.assert_allclose(got.to_numpy(), want.to_numpy(), rtol=1e-12, atol=0.0)
    assert fired["ruler"] >= 1 and fired["max_abs_delta"] > 0.0
    assert not sub["count"].equals(base["count"]), (
        "a quarter of the length on two isoforms moved no count — the lengths never reached the split"
    )

    noop, fired_noop = _quant_ruler(toy, "oracle_ruler_noop", factor)
    assert fired_noop["ruler"] >= 1 and fired_noop["truth_delta"] > 0.0
    assert fired_noop["max_abs_delta"] == 0.0
    exact = [c for c in base.columns if c != "posterior_mean"]
    pd.testing.assert_frame_equal(base[exact], noop[exact], check_exact=True)
    pd.testing.assert_series_equal(
        base["posterior_mean"], noop["posterior_mean"], check_exact=False, rtol=1e-9, atol=0.0
    )


def test_the_ruler_guard_refuses_an_arm_that_could_not_have_measured():
    """The arm must move the lengths wherever the truth differs from the shipped ruler, its noop must not move them, and
    a shrinkage that never ran is no measurement. Where the truth is the shipped ruler (capture-OFF: the plain length on
    both sides) the arm is the identity and says so rather than refusing."""
    fired = dict(n=1, ruler=1, max_abs_delta=3.0, truth_delta=3.0)
    qa.check_ruler_arm("c", "oracle_ruler", fired)
    qa.check_ruler_arm("c", "oracle_ruler", dict(fired, max_abs_delta=0.0, truth_delta=0.0))
    with pytest.raises(RuntimeError, match="did not move"):
        qa.check_ruler_arm("c", "oracle_ruler", dict(fired, max_abs_delta=0.0))
    with pytest.raises(RuntimeError, match="NOOP"):
        qa.check_ruler_arm("c", "oracle_ruler_noop", dict(fired, max_abs_delta=1e-9))
    with pytest.raises(RuntimeError, match="never ran"):
        qa.check_ruler_arm("c", "oracle_ruler", dict(fired, ruler=0))


@pytest.fixture(scope="module")
def toy_zero(tmp_path_factory):
    """The zero-gDNA toy: no reference carries gDNA."""
    from _prior_toy import build_toy_zero_gdna

    return build_toy_zero_gdna(tmp_path_factory)


def _ruler_condition(toy, suite, n_gdna, monkeypatch):
    """``oracle_ruler`` installed as ``run_condition`` installs it, on a panel whose manifest names the toy's one
    reference as genomic, under a condition that drew ``n_gdna`` gDNA fragments, with the simulator's factor stood in
    by a probed spike-in's (about 1000) on every annotated transcript. Returns the run's transcript table, its fired
    counter and the factor."""
    from rigel.pipeline import run_pipeline

    (ref,) = toy.index.ref_name_to_id
    suite.mkdir(exist_ok=True)
    (suite / "manifest.json").write_text(
        json.dumps(
            {"gdna": {"genomic_refs": [ref]}, "conditions": [{"name": "c", "n_gdna": n_gdna}]}
        )
    )
    factor = np.where(toy.index.t_df["is_synthetic"].to_numpy(bool), 1.0, 1000.0)
    monkeypatch.setattr(qa, "capture_truth_factor", lambda *a: factor)
    cfg = PipelineConfig()
    restore, fired = qa.install_condition_arm(
        "oracle_ruler", suite, toy.index, "c", cfg, None, None
    )
    try:
        result = run_pipeline(str(toy.bam_path), toy.index, cfg)
    finally:
        restore()
    return result.estimator.get_counts_df(toy.index), fired, factor


def test_the_ruler_arm_keeps_the_shipped_length_where_the_condition_drew_no_gdna(
    toy, toy_zero, tmp_path, monkeypatch
):
    """The arm leaves the locus gDNA length as shipped, and its factor is anchored on that length's scale: where gDNA
    was drawn, the shipped gDNA length reads gDNA's own capture. Where none was drawn (an RNA-only spike-in reference,
    or every reference of a zero-gDNA condition) nothing measured capture, the shipped gDNA length is the uncaptured
    one, and a probed transcript at ``fl × factor`` would sit three orders of magnitude above it: the E-step then
    hands the uncaptured components every fragment they explain alike (on the ladder, about every probed spike-in
    fragment to gDNA). There the arm must hand the EM the shipped lengths, read off the EM's own published length and
    counts; where gDNA was drawn it substitutes as before. Which references drew gDNA is read off the manifest."""
    base, _ = _quant(toy_zero, "base", None)
    held, fired, _ = _ruler_condition(toy_zero, tmp_path / "zero", 0, monkeypatch)
    assert fired["ruler"] >= 1
    np.testing.assert_array_equal(held["em_effective_length"], base["em_effective_length"])
    exact = [c for c in base.columns if c != "posterior_mean"]
    pd.testing.assert_frame_equal(base[exact], held[exact], check_exact=True)

    shipped, _ = _quant(toy, "base", None)
    sub, fired, factor = _ruler_condition(toy, tmp_path / "contaminated", 1, monkeypatch)
    ids = toy.index.t_df["t_id"].to_numpy()
    plain = shipped.set_index("transcript_id")["em_effective_length"]
    want = plain * pd.Series(factor, index=ids).reindex(plain.index)
    got = sub.set_index("transcript_id")["em_effective_length"].reindex(plain.index)
    np.testing.assert_allclose(got.to_numpy(), want.to_numpy(), rtol=1e-12, atol=0.0)
    assert fired["max_abs_delta"] > 0.0


def test_the_references_without_gdna_are_the_manifests_never_the_counts():
    """A spike-in reference outside ``gdna.genomic_refs`` holds no gDNA in any condition, and a zero-gDNA condition
    holds none on any reference; a condition the manifest does not list is refused rather than read as one."""
    t_df = pd.DataFrame({"ref": ["chr1", "chr1", "ERCC-1", "chr2"]})
    index = type("I", (), {"t_df": t_df, "ref_name_to_id": {"chr1": 0, "ERCC-1": 1, "chr2": 2}})()
    manifest = {
        "gdna": {"genomic_refs": ["chr1", "chr2"]},
        "conditions": [{"name": "g05", "n_gdna": 5}, {"name": "g00", "n_gdna": 0}],
    }
    assert qa.gdna_free_transcripts(manifest, index, "g05").tolist() == [False, False, True, False]
    assert qa.gdna_free_transcripts(manifest, index, "g00").tolist() == [True] * 4
    with pytest.raises(KeyError):
        qa.gdna_free_transcripts(manifest, index, "g50")


def test_the_no_capture_arm_replaces_the_readers_answer_and_nothing_else(monkeypatch):
    """The arm's wrapper sits on the one place the result takes its capture reading. A reading with capture in it (a
    saturated level of 50, efficiencies of 7 and 3) comes out as no capture at all: every efficiency exactly 1 and the
    saturated level 1, with calibration's bulk density kept, since the bulk is not the reader's answer about capture.
    After ``restore`` the reading is the reader's again."""
    import importlib

    CALIBRATE = importlib.import_module("rigel.calibration.calibrate")
    read = dict(
        gdna_bulk_density=2e-4,
        gdna_capture_saturation=50.0,
        gdna_capture_efficiency_region=np.array([7.0, 1.0, 0.2]),
        gdna_capture_efficiency_boundary=np.array([3.0, 0.5]),
    )
    monkeypatch.setattr(CALIBRATE, "_capture_fields", lambda *a, **k: dict(read))
    restore, fired = qa.install_no_capture()
    try:
        got = CALIBRATE._capture_fields()
    finally:
        restore()
    assert fired["n"] == 1
    assert got["gdna_bulk_density"] == 2e-4 and got["gdna_capture_saturation"] == 1.0
    np.testing.assert_array_equal(got["gdna_capture_efficiency_region"], np.ones(3))
    np.testing.assert_array_equal(got["gdna_capture_efficiency_boundary"], np.ones(2))
    assert CALIBRATE._capture_fields()["gdna_capture_saturation"] == 50.0, (
        "restore left the wrapper installed"
    )


def test_the_no_capture_arm_reaches_the_result_and_its_guard_refuses_one_it_did_not(toy):
    """End to end the arm's claim is read off the run's own calibration result, not the wrapper's return value, which
    would stay healthy if the result were assembled from another reading. A result holding any efficiency other than
    1, or a saturated level above 1, is refused as no measurement of reading no capture."""
    from rigel.pipeline import run_pipeline

    restore, fired = qa.install_no_capture()
    try:
        result = run_pipeline(str(toy.bam_path), toy.index, PipelineConfig())
    finally:
        restore()
    assert fired["n"] >= 1
    cal = result.calibration
    qa.check_no_capture("toy", cal)
    c = np.array(cal.gdna_capture_efficiency_region)
    c[0] = 2.0
    with pytest.raises(RuntimeError, match="no_capture"):
        qa.check_no_capture("toy", dataclasses.replace(cal, gdna_capture_efficiency_region=c))
    if cal.gdna_bulk_density is not None:
        with pytest.raises(RuntimeError, match="no_capture"):
            qa.check_no_capture("toy", dataclasses.replace(cal, gdna_capture_saturation=3.0))


def test_the_capture_truth_reads_1_at_the_uncaptured_bulk():
    """The shipped efficiencies are yields against the uncaptured bulk, so the EM's lengths put an unprobed transcript
    at its plain length and a probed one above it, on the scale of the locus gDNA length the arm leaves as shipped. The
    simulator's factor is anchored on the same class: the unprobed transcripts' median reads 1, and a probed one keeps
    its whole yield above them. Anchored on the fully probed class instead, every unprobed transcript's length shrank by
    the saturation against the shipped gDNA length, and the arm read several times the RNA truth under capture. A
    transcript with no length is the identity; capture-OFF every factor is 1."""
    factor = np.array([1000.0, 1000.0, 1200.0, 2.0, 2.0, 2.0, 50.0, np.nan])
    probed = np.array([1.0, 0.95, 0.92, 0.0, 0.0, 0.0, 0.1, 0.0])
    plain = np.array([100.0] * 7 + [0.0])
    got = qa.anchor_capture_factor(factor, probed, plain)
    np.testing.assert_allclose(got, [500.0, 500.0, 600.0, 1.0, 1.0, 1.0, 25.0, 1.0], rtol=1e-12)
    off = qa.anchor_capture_factor(np.full(4, 3.0), np.zeros(4), np.full(4, 100.0))
    assert np.array_equal(off, np.ones(4))


def test_a_panel_that_probes_every_transcript_keeps_the_simulators_own_scale():
    """With no uncaptured transcript there is no class to read the bulk off; the simulator's factor carries the
    off-target weight in both its counts, so it is already a yield against the bulk and is kept as it is."""
    factor = np.array([1000.0, 800.0, 20.0])
    got = qa.anchor_capture_factor(factor, np.array([1.0, 0.9, 0.05]), np.full(3, 100.0))
    np.testing.assert_allclose(got, factor, rtol=1e-12)


# ── Gate 8: arms scored under different assignment modes are refused, never compared ──────────────


def test_the_report_refuses_arms_scored_under_different_assignment_modes(tmp_path):
    """``--set em.assignment_mode=fractional`` changes every count, so a fractional arm beside a sampled one would
    report the mode's effect as the arm's. Every row carries its mode and the report refuses a mix; a row without the
    stamp counts as a different mode."""
    row = dict(
        condition="gdna_g05_ss_0.99_nascent_mid_capture_on", axis="transcript", count_abs_err=1.0
    )
    a, b = tmp_path / "a.jsonl", tmp_path / "b.jsonl"
    a.write_text(json.dumps(dict(row, assignment_mode="fractional")) + "\n")
    b.write_text(json.dumps(dict(row, assignment_mode="sample")) + "\n")
    with pytest.raises(SystemExit, match="assignment mode"):
        qa.report([a, b])
    b.write_text(json.dumps(row) + "\n")
    with pytest.raises(SystemExit, match="assignment mode"):
        qa.report([a, b])


# ── the markdown report: what a release is judged on ────────────────────────────────────────────
#
# The report renders numbers the arm files already hold, so every way of getting it wrong produces a plausible-looking
# table: the wrong gDNA column, a ratio invented where the truth is zero, a pool folded into another. These tests cover
# those.


def _arm_rows(tmp_path, **over):
    """One condition, three axes, written as an arm jsonl, with numbers chosen so every quantity the report derives is
    distinguishable from every other."""
    cond = "gdna_g50_ss_0.99_nascent_mid_capture_off"
    lib = {
        "condition": cond,
        "axis": "library",
        "arm": "a",
        "assignment_mode": "fractional",
        "annotated_est": 4_000.0,
        "annotated_true": 4_100.0,
        "nascent_est": 900.0,
        "nascent_true": 1_000.0,
        "rna_est": 4_900.0,
        "rna_true": 5_100.0,
        #: EM-only. The comparable total is this plus the intergenic fragments.
        "gdna_est": 2_000.0,
        "n_intergenic_gdna": 3_100.0,
        "gdna_true": 4_900.0,
        "gdna_frac_est": 0.51,
        "gdna_frac_true": 0.49,
    }
    lib.update(over)
    axis = {
        "condition": cond,
        "arm": "a",
        "assignment_mode": "fractional",
        "n_tx": 10,
        "n_expressed": 6,
        "n_detected": 7,
        "count_true": 4_100.0,
        "count_est": 4_000.0,
        "count_abs_err": 410.0,
        "count_net_err": -100.0,
        "count_over": 155.0,
        "count_under": 255.0,
        "fp_mass": 33.0,
        "fp_n": 2,
        "fn_mass": 11.0,
        "fn_n": 1,
        "spearman": 0.9,
        "pearson_log2": 0.8,
        "mard": 0.25,
    }
    tx = dict(axis, axis="transcript", median_rel_err=0.12, tpm_abs_err=1.0, tpm_fp=0.5)
    gene = dict(axis, axis="gene", n_tx=4, count_abs_err=41.0)
    p = tmp_path / "qa_arm.jsonl"
    p.write_text("\n".join(json.dumps(r) for r in (lib, tx, gene)) + "\n")
    return p, cond


def _render(tmp_path, **over):
    src, _cond = _arm_rows(tmp_path, **over)
    out = tmp_path / "report.md"
    qa.markdown_report(src, out)
    return out.read_text()


def test_the_report_scores_gdna_including_the_intergenic_fragments(tmp_path):
    """``gdna_est`` is the EM's number and excludes fragments that reached no locus, while ``gdna_true`` counts every
    gDNA fragment. Off capture the intergenic share is the larger one, so printing ``gdna_est`` against that truth would
    read as a large under-call that is not there. The comparable estimate is ``gdna_est + n_intergenic_gdna``.
    """
    md = _render(tmp_path)
    assert "5,100" in md, "the gDNA estimate did not include the intergenic fragments"
    assert "| 2,000 |" not in md, "the report printed the EM's gDNA alone against the full truth"
    assert "+4.1 %" in md, f"expected (5,100-4,900)/4,900; got:\n{md}"


def test_a_percent_error_is_refused_where_the_truth_is_zero(tmp_path):
    """The `g00` rung's true gDNA count is exactly 0, and a percent error there would invent a scale. The raw count is
    the whole answer; the ratio is `n/a`."""
    md = _render(tmp_path, gdna_true=0.0)
    assert "n/a" in md
    assert "inf" not in md.lower() and "nan %" not in md.lower()


def test_the_three_pools_are_reported_apart_and_the_synthetic_one_by_name(tmp_path):
    """The split is on ``is_synthetic``, and the report must say so: a single-exon annotated transcript
    (``is_single_exon``) belongs to the annotated pool. Folding the two RNA pools together would hide the false-positive
    channel the nascent row exposes."""
    md = _render(tmp_path)
    assert "SYNTHETIC entities only" in md
    assert "is_synthetic`, never `is_single_exon" in md
    for n in ("900", "1,000", "4,000", "4,100"):  # nascent est/true, annotated est/true
        assert n in md, f"pool figure {n} missing"


def test_the_transcript_and_gene_sections_score_the_annotated_pool_alone(tmp_path):
    """The transcript table a user reads holds no synthetic entity (``get_counts_df`` drops them), so the transcript and
    gene truth must equal the library's annotated truth, not the RNA total; scoring against ``rna_true`` would credit
    the tool for nascent mass not in the table."""
    md = _render(tmp_path)
    src, cond = _arm_rows(tmp_path)
    rows = qa._load(src)
    assert rows[(cond, "transcript")]["count_true"] == rows[(cond, "library")]["annotated_true"]
    assert rows[(cond, "transcript")]["count_true"] != rows[(cond, "library")]["rna_true"]
    assert "10.00 %" in md, "Σ|Δ| as a share of the annotated truth (410/4,100) is missing"


def test_every_table_in_the_report_is_rectangular(tmp_path):
    """A ragged markdown table renders as a wrong table rather than an error, and the escaped pipes in a `Σ|Δ|` header
    are cell content, not delimiters."""
    md = _render(tmp_path)
    for block in md.split("\n\n"):
        lines = [x for x in block.splitlines() if x.startswith("|")]
        if len(lines) < 2:
            continue
        widths = {len(x.replace("\\|", "").split("|")) for x in lines}
        assert len(widths) == 1, f"ragged table {widths}: {lines[0][:80]}"


def test_the_deferred_stratum_is_marked_wherever_it_appears(tmp_path):
    """Unstranded × capture-ON is reported on every benchmark, is never a development target, and carries most of the
    error, so a reader who misses the mark would read a pooled total as the tool's accuracy."""
    lib, tx, gene = (json.loads(x) for x in _arm_rows(tmp_path)[0].read_text().splitlines())
    for r in (lib, tx, gene):
        r["condition"] = "gdna_g50_ss_0.50_nascent_mid_capture_on"
    src = tmp_path / "def.jsonl"
    src.write_text("\n".join(json.dumps(r) for r in (lib, tx, gene)) + "\n")
    out = tmp_path / "d.md"
    qa.markdown_report(src, out)
    md = out.read_text()
    assert md.count("DEFERRED") >= 4, "the deferred stratum is unmarked in at least one section"


def test_the_report_presents_the_false_negative_mass_beside_the_under_assignment(tmp_path):
    """``fn_mass`` is the truth on the expressed transcripts the estimate gives exactly zero. A component the EM left at
    its floor is assigned nothing, so the column fires under fractional assignment too, and leaving it out would hide
    real RNA the tool reports as absent. It is a part of ``count_under``, which stays beside it.
    """
    md = _render(tmp_path)
    assert "cannot fire" not in md, "the report still says a false negative cannot be measured"
    for section in ("## 2.", "## 3."):
        table = [x for x in md.split(section)[1].split("## ")[0].splitlines() if x.startswith("|")]
        head, row = (
            [c.strip() for c in x.replace("\\|", "").split("|")] for x in (table[0], table[-1])
        )
        at = head.index("false-negative mass")
        assert (row[at], row[at + 1]) == ("11", "1"), "fn_mass and fn_n are not under their header"
        assert row[head.index("under-assigned")] == "255", "count_under is missing"
        assert row[head.index("false-positive mass")] == "33", "fp_mass moved from under its header"


def test_the_report_does_not_call_the_truth_tables_row_count_a_transcript_count(tmp_path):
    """The truth table carries one row per synthetic nascent entity as well as per annotated transcript, and those rows
    are zero on both sides, so they enter no figure. Printing the row count as "transcripts" would claim thousands of
    scored transcripts that were never scored. The scored sets are ``expressed`` and ``detected``.
    """
    md = _render(tmp_path)
    body = md.split("## 2.")[1]
    assert "| 10 |" not in body, "the truth table's row count was reported as a transcript count"
    assert "| 4 |" not in body, "the gene axis reported its grouping-key count"
    assert "expressed" in body and "detected" in body


def _arm_file(tmp_path, gdna_frac_true: dict) -> Path:
    """``_arm_rows``'s three rows once per condition, each library row carrying its own true share."""
    base = [json.loads(x) for x in _arm_rows(tmp_path)[0].read_text().splitlines()]
    rows = [
        {**r, "condition": cond, **({"gdna_frac_true": f} if r["axis"] == "library" else {})}
        for cond, f in gdna_frac_true.items()
        for r in base
    ]
    path = tmp_path / "multi.jsonl"
    path.write_text("\n".join(json.dumps(r) for r in rows) + "\n")
    return path


def _gdna_pool_blocks(md: str) -> dict:
    """The gDNA pool table's stratum headings, each with the conditions filed under it."""
    table = md.split("### gDNA (EM + intergenic)\n\n")[1].split("\n\n")[0]
    blocks, head = {}, None
    for line in table.splitlines():
        if line.startswith("| **"):
            head = line.split("**")[1]
            blocks[head] = []
        elif line.startswith("| `"):
            blocks[head].append(line.split("`")[1])
    return blocks


def test_every_gdna_rung_is_labelled_from_its_truth_not_from_a_map_of_known_rungs(tmp_path):
    """The level column is read off the row's own truth, ``gdna_frac_true``, so a panel with any rung (the test
    chromosome's ``g25``, the depth panels' ``g001`` and ``g01``) gets a report. A hand map of the ladder's rung tokens
    would raise on ``g25``, and digits read from the name would print ``g001`` as 1 %.
    """
    rungs = {"g00": (0.0, "0 %"), "g001": (0.001, "0.1 %"), "g05": (0.0500004, "5 %"),
             "g25": (0.25, "25 %"), "g98": (0.98, "98 %")}  # fmt: skip
    conds = {f"gdna_{k}_ss_0.99_nascent_mid_capture_off": f for k, (f, _label) in rungs.items()}
    out = tmp_path / "rungs.md"
    qa.markdown_report(_arm_file(tmp_path, conds), out)
    md = out.read_text()
    for (k, (_f, label)), cond in zip(rungs.items(), conds):
        assert f"| `{cond}` | {label} |" in md, f"rung {k} is not labelled {label}"


def test_a_strand_specificity_in_neither_half_is_reported_apart_never_pooled(tmp_path):
    """ss 0.50 is the unstranded half and ss 0.99 the stranded half, judged against different bars. The test panel's ss
    0.70 is neither, so it is reported apart rather than moving the half a policy must win on. One definition
    (`_shared.strandedness`) serves every instrument, and a name that spells no specificity raises rather than landing
    in a half.
    """
    assert _shared.strandedness("gdna_g50_ss_0.50_nascent_mid_capture_on") == "unstranded"
    assert _shared.strandedness("gdna_g50_ss_0.99_nascent_mid_capture_on") == "stranded"
    assert _shared.strandedness("gdna_g50_ss_0.70_nascent_file_capture_on") == "ss 0.70"
    with pytest.raises(ValueError):
        _shared.strandedness("gdna_g50_nascent_mid_capture_on")

    c50, c70 = (
        "gdna_g50_ss_0.50_nascent_file_capture_off",
        "gdna_g50_ss_0.70_nascent_file_capture_off",
    )
    assert _shared.strata([c50, c70]) == [*_shared.STRATA, ("ss 0.70", "capture OFF")]
    out = tmp_path / "apart.md"
    qa.markdown_report(_arm_file(tmp_path, {c50: 0.5, c70: 0.5}), out)
    md = out.read_text()
    assert _gdna_pool_blocks(md) == {
        "unstranded × capture OFF": [c50],
        "ss 0.70 × capture OFF": [c70],
    }
    assert "| ss 0.70 × capture OFF |" in md.split("## 4.")[1], (
        "the roll-up dropped the apart stratum"
    )


# ── Gate 16: the allocation arm's weights must reach every transcript the simulator gave RNA ──────


def _alloc_truth_table():
    """Two annotated transcripts and the synthetic span of each of their genes, as the simulator writes them: a
    synthetic span's fragments are ``observed_nascent_fragments`` on its own row, and its ``observed_annotated_fragments`` is
    zero because a synthetic span emits no mature RNA. A single-exon annotated transcript carries both columns."""
    return pd.DataFrame(
        {
            "transcript_id": ["t_spliced", "t_single_exon", "nascent_gene_a", "nascent_gene_b"],
            "observed_annotated_fragments": [900.0, 100.0, 0.0, 0.0],
            "observed_nascent_fragments": [0.0, 7.0, 300.0, 50.0],
        }
    )


class _AllocIndex:
    """The three columns ``truth_weights`` reads off an index."""

    def __init__(self, t_ids):
        self.t_df = pd.DataFrame(
            {"t_id": list(t_ids), "t_index": np.arange(len(t_ids), dtype=np.int64)}
        )
        self.num_transcripts = len(t_ids)


def _weighted_counts(toy, weights):
    """The toy's transcript table under ``oracle_alloc_seed``'s injection of ``weights`` (``None``: no injection), with
    the guard's counter."""
    from rigel.pipeline import run_pipeline

    restore, fired = (lambda: None), {"n": 0}
    if weights is not None:
        restore, fired = qa.install_truth_weights(weights)
    try:
        result = run_pipeline(str(toy.bam_path), toy.index, PipelineConfig())
    finally:
        restore()
    return result.estimator.get_counts_df(toy.index)["count"].to_numpy(), fired


def test_the_allocation_arm_moves_the_em_and_its_guard_counts_only_weights_the_em_reads(toy):
    """``oracle_alloc_seed`` sets ``rna_prior_weight`` on the estimator's batch call, which the pipeline never passes,
    so weights must move the EM's answer, or every allocation ceiling the arm reports is zero. The native prior shares a
    locus whose weights sum to zero by the evidence rule, the shipped one, so weights on no transcript the EM solves are
    a silent ``base``: the guard counts only a call whose weights reach a solved transcript, and ``run_condition``
    refuses a count of zero.
    """
    n = int(toy.index.num_transcripts)
    base, _ = _weighted_counts(toy, None)
    uniform, fired = _weighted_counts(toy, np.ones(n))
    assert fired["n"] == 1
    assert not np.array_equal(uniform, base), (
        "uniform weights did not move the EM: the injection is lost"
    )
    missed, fired = _weighted_counts(toy, np.zeros(n))
    assert fired["n"] == 0, "weights the EM never reads were counted as an injection"
    np.testing.assert_array_equal(missed, base)  # the silent base the guard exists to refuse


def test_the_allocation_arm_never_weighs_an_rna_emitting_transcript_at_zero():
    """The arm must weigh every RNA component, and a synthetic span is one.

    ``--arm oracle_alloc_seed`` claims to hand the EM the true per-transcript allocation. Reading a mature column alone
    would give every synthetic nascent entity a weight of zero, which is not a true allocation but a suppression of the
    component, and the arm would then appear to remove the nascent over-call. A transcript the simulator gave RNA
    fragments must receive a positive weight, and one it gave none must receive zero.
    """
    truth = _alloc_truth_table()
    index = _AllocIndex(truth["transcript_id"])
    w = qa.truth_weights(truth, index)

    emits = (truth["observed_annotated_fragments"] + truth["observed_nascent_fragments"]).to_numpy()
    assert np.all(w[emits > 0] > 0.0), (
        "a transcript the simulator gave RNA fragments received a weight of ZERO — a zero weight is "
        "a suppression, not an allocation"
    )
    assert np.all(w[emits == 0] == 0.0), "a silent transcript received prior mass"
    # the weights are the RNA the simulator realised, mature and nascent alike
    np.testing.assert_allclose(w, [900.0, 107.0, 300.0, 50.0], rtol=1e-12)


def test_the_allocation_weights_refuse_a_truth_table_that_carries_no_nascent_column():
    """A truth table with a mature column and no nascent one cannot express a synthetic span's allocation, so the arm
    refuses rather than weighting every synthetic span at zero. The suite-level molar table is refused too: its
    ``annotated_abundance`` / ``nascent_abundance`` pair counts molecules, not the fragments a condition realised."""
    truth = _alloc_truth_table().drop(columns=["observed_nascent_fragments"])
    with pytest.raises(SystemExit, match="nascent"):
        qa.truth_weights(truth, _AllocIndex(truth["transcript_id"]))
    molar = _alloc_truth_table().set_axis(
        ["transcript_id", "annotated_abundance", "nascent_abundance"], axis=1
    )
    with pytest.raises(SystemExit, match="nascent"):
        qa.truth_weights(molar, _AllocIndex(molar["transcript_id"]))


# ── an aligned input, scored against the same truth, its aligner's losses counted ────────────


def _templates(bam_path):
    """The BAM's records grouped by template, in file order: ``[(name, [AlignedSegment, ...]), ...]``."""
    import pysam

    groups = []
    with pysam.AlignmentFile(str(bam_path), "rb", check_sq=False) as fh:
        header = fh.header.to_dict()
        for r in fh.fetch(until_eof=True):
            if not groups or groups[-1][0] != r.query_name:
                groups.append((r.query_name, []))
            groups[-1][1].append(r)
    return header, groups


def _write(path, header, groups):
    import pysam

    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for _name, recs in groups:
            for r in recs:
                out.write(r)
    return path


def _unmap(r):
    """The record as an aligner writes a read it could not place, its pair kept (`outSAMunmapped Within`)."""
    r.flag = (r.flag | 0x4 | 0x8) & ~0x2
    return r


def _origin_counts(groups):
    from rigel.sim.read_name import parse_origin

    counts = {"gdna": 0, "annotated": 0, "nascent": 0}
    for name, _recs in groups:
        counts[parse_origin(name).kind] += 1
    return counts


def test_the_aligner_loss_counts_unmapped_templates_by_their_true_origin(toy, tmp_path):
    """An aligned input is scored against the simulator's truth, which counts every fragment it drew, so a row read from
    it must carry what the aligner left out, by origin: on the ladder's suite genome part of the gDNA is drawn over N
    runs no aligner can place, and read against the full truth that would look like Rigel's error. A template is lost
    only when none of its records is mapped; one placed mate is the aligner's delivery and Rigel's to judge. Counting
    records instead of templates, or reading the origin from anything but the read name, would move the counts.
    """
    from rigel.sim.read_name import parse_origin

    header, groups = _templates(toy.bam_path)
    expected = {"gdna": 0, "annotated": 0, "nascent": 0}
    half = None
    for i, (name, recs) in enumerate(groups):
        if i % 5 == 0:
            for r in recs:
                _unmap(r)
            expected[parse_origin(name).kind] += 1
        elif half is None and len(recs) >= 2:
            _unmap(recs[0])  # one mate unplaced: the template still reached Rigel
            half = name
    assert half is not None and min(expected.values()) > 0, (
        "the toy must lose templates of every origin"
    )
    bam = _write(tmp_path / "aligned.bam", header, groups)
    got = qa.aligner_loss(bam, _origin_counts(groups))
    assert got == {f"unaligned_{k}": v for k, v in expected.items()}


def test_an_aligned_bam_that_is_not_the_simulators_is_refused(toy, tmp_path):
    """The aligned input must hold the simulator's templates, each once and in one group: a missing template (a
    truncated alignment) or one split into two runs of records (an input that is not collated) would shift the aligner's
    losses silently."""
    header, groups = _templates(toy.bam_path)
    counts = _origin_counts(groups)
    short = _write(tmp_path / "short.bam", header, groups[:-1])
    with pytest.raises(SystemExit, match="not the simulator's"):
        qa.aligner_loss(short, counts)
    name, recs = next((n, r) for n, r in groups if len(r) >= 2)
    split = [g for g in groups if g[0] != name] + [(name, recs[:1])]
    split.insert(0, (name, recs[1:]))
    with pytest.raises(SystemExit, match="not the simulator's"):
        qa.aligner_loss(_write(tmp_path / "split.bam", header, split), counts)


def test_another_input_bam_is_scored_by_the_truth_free_arms_alone_and_must_exist(tmp_path):
    """The oracle arms inject truth built from the simulator's own alignments (the oracle cache), so on another input
    their substitution would be another scan's and is refused. ``base`` and ``no_capture`` substitute no truth and
    score any BAM of the reads. A condition without the named BAM is refused by name before anything runs, never
    scored as a shorter panel."""
    for c in ("c1", "c2"):
        (tmp_path / c).mkdir()
        (tmp_path / c / qa.ORACLE_BAM).write_bytes(b"")
    (tmp_path / "c1" / "aligned.bam").write_bytes(b"")
    with pytest.raises(SystemExit, match="c2"):
        qa.check_input_bam(tmp_path, ["c1", "c2"], "base", "aligned.bam")
    qa.check_input_bam(tmp_path, ["c1"], "base", "aligned.bam")  # every condition has it
    qa.check_input_bam(tmp_path, ["c1"], "no_capture", "aligned.bam")  # no truth substituted
    with pytest.raises(SystemExit, match="base"):
        qa.check_input_bam(tmp_path, ["c1"], "oracle", "aligned.bam")
    qa.check_input_bam(tmp_path, ["c1", "c2"], "oracle", qa.ORACLE_BAM)  # the simulator's: any arm
