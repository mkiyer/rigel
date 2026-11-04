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
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from rigel.config import EMConfig, PipelineConfig
from _prior_toy import build_toy

_MODULES: dict = {}


def _load_sibling(name: str):
    """Import a ``scripts/design/`` instrument by path. Registered in ``sys.modules`` before
    execution, or ``@dataclass`` fails resolving its own ``__module__``."""
    import sys

    key = name[:-3]
    if key not in _MODULES:
        path = Path(__file__).resolve().parents[2] / "scripts" / "design" / name
        spec = importlib.util.spec_from_file_location(key, path)
        module = importlib.util.module_from_spec(spec)
        sys.modules[key] = module
        _MODULES[key] = module
        spec.loader.exec_module(module)
    return _MODULES[key]


QA = _load_sibling("quant_accuracy.py")


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
    restore, fired = QA.install_arm(arm, oracle)
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


def test_the_noop_arm_reproduces_BASE_byte_identically_through_the_whole_pipeline(toy, toy_oracle):
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


def test_base_IS_the_shipped_config_and_a_pinned_arm_reproduces_itself(toy):
    """``base`` must be the shipped configuration, so the instrument's seed is the shipped one, and two runs of one arm
    must return one answer, or no arm comparison is attributable.
    """
    a, _ = _quant(toy, "base", None)
    b, _ = _quant(toy, "base", None)
    assert np.array_equal(a["count"].to_numpy(), b["count"].to_numpy()), (
        "two runs of one arm disagree — no arm comparison here is attributable"
    )


def test_an_injection_that_never_fires_RAISES_rather_than_reporting_no_effect(
    tmp_path, monkeypatch
):
    """ "The perfect prior changed nothing" is a publishable conclusion, so a wrapper that was never called must not be
    able to produce it; that happens when ``quant_from_buffer`` returns early on an empty unit set. ``run_pipeline`` is
    stubbed to a no-op, the shape of that early return: the pipeline completes, the wrapper is installed, and
    ``assemble_priors`` is never reached.
    """
    cond = tmp_path / "c1"
    cond.mkdir()
    pd.DataFrame({"transcript_id": ["a"], "observed_mrna_fragments": [1]}).to_csv(
        cond / "truth_abundances.tsv", sep="\t", index=False
    )
    (cond / "truth_summary.json").write_text('{"origin_counts": {"gdna": 0, "mrna": 1, "nrna": 0}}')
    monkeypatch.setattr(QA, "run_pipeline", lambda *a, **k: None)
    monkeypatch.setattr(QA, "load_oracle", lambda *a, **k: object())
    with pytest.raises(RuntimeError, match="never wrapped-and-called"):
        QA.run_condition("oracle", tmp_path, None, "c1", PipelineConfig(), tmp_path)


# ── Gate 3: the truth column is the realised fragment count, not the molar abundance ─────────────


def test_the_scorer_reads_the_realised_fragment_column_and_refuses_the_molar_table(tmp_path):
    """The panel ships two truth tables holding different quantities. A condition's ``truth_abundances.tsv`` gives the
    realised observed fragment count (``observed_mrna_fragments``) beside the pre-capture molar abundance
    (``pre_capture_mrna_abundance``); the suite-level ``truth_abundances_nrna_<label>.tsv`` gives the molar abundance
    alone, as ``mrna_abundance``. Scoring a fragment estimate against a molar truth would charge the tool for hybrid
    capture, which it never claims to invert.

    Both tables are written by the simulator's own writers and read as the instrument reads them. The scorer must score
    the realised count exactly, refuse the molar table, and see a molar truth placed in the realised column as error.
    """
    from rigel.sim.truth import write_molar_truth, write_post_capture_truth
    from rigel.transcript import Transcript
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
    write_molar_truth(transcripts, tmp_path / "truth_abundances_nrna_none.tsv")
    realised = pd.read_csv(cond / "truth_abundances.tsv", sep="\t")
    molar = pd.read_csv(tmp_path / "truth_abundances_nrna_none.tsv", sep="\t")

    quant = pd.DataFrame(
        {
            "transcript_id": ["a", "b", "c"],
            "count": [1000.0, 500.0, 0.0],
            "count_em": [0.0, 0.0, 0.0],
            "tpm": [0.0, 0.0, 0.0],
            "effective_length": [100.0, 100.0, 100.0],
        }
    )
    exact = QA.score_transcripts(quant, realised)
    assert exact["count_abs_err"] == pytest.approx(0.0)
    assert exact["fp_mass"] == pytest.approx(0.0)

    with pytest.raises(KeyError, match="observed_mrna_fragments"):
        QA.score_transcripts(quant, molar)

    wrong = QA.score_transcripts(
        quant, realised.assign(observed_mrna_fragments=realised["pre_capture_mrna_abundance"])
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
        {"transcript_id": ["a", "b", "c"], "observed_mrna_fragments": [100.0, 0.0, 50.0]}
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
    s = QA.score_transcripts(quant, truth)
    assert s["fp_mass"] == pytest.approx(30.0) and s["fp_n"] == 1
    assert s["fn_mass"] == pytest.approx(50.0) and s["fn_n"] == 1
    assert s["count_over"] - s["count_under"] == pytest.approx(s["count_net_err"])
    assert s["count_over"] + s["count_under"] == pytest.approx(s["count_abs_err"])

    lifted = truth.assign(observed_mrna_fragments=[100.0, 25.0, 50.0])
    s2 = QA.score_transcripts(quant, lifted)
    assert s2["fp_mass"] == pytest.approx(0.0), (
        "a transcript the truth now expresses is still counted as a false positive — fp_mass is "
        "keyed on the estimate, not on the truth"
    )


def test_the_GENE_row_collapses_isoform_ambiguity_and_still_sees_real_gene_error():
    """The gene row answers "is the remaining transcript error isoform ambiguity, or something else?", which it can only
    do if moving mass between a gene's own isoforms is invisible to it while moving mass out of the gene is not. Both
    halves are asserted, since a scorer that saw neither would also report zero on the first.
    """
    truth = pd.DataFrame(
        {
            "transcript_id": ["a1", "a2", "b1"],
            "gene_id": ["A", "A", "B"],
            "observed_mrna_fragments": [100.0, 0.0, 50.0],
        }
    )
    swapped = pd.DataFrame(
        {"transcript_id": ["a1", "a2", "b1"], "count": [0.0, 100.0, 50.0]}
    )  # the isoform call is completely wrong; the gene call is perfect
    tx = QA.score_transcripts(swapped.assign(count_em=0.0, tpm=0.0, effective_length=100.0), truth)
    gene = QA.score_genes(swapped, truth)
    assert tx["count_abs_err"] == pytest.approx(200.0)
    assert gene["count_abs_err"] == pytest.approx(0.0), (
        "a pure isoform swap inside one gene shows up at gene level — the grouping is not collapsing "
        "isoform ambiguity and the gene row answers nothing the transcript row did not"
    )

    leaked = pd.DataFrame({"transcript_id": ["a1", "a2", "b1"], "count": [60.0, 0.0, 90.0]})
    gene2 = QA.score_genes(leaked, truth)
    assert gene2["count_abs_err"] == pytest.approx(80.0), (
        "40 fragments moved from gene A to gene B and the gene row did not notice — it is summing "
        "away real error, not ambiguity"
    )


def test_the_library_gDNA_fraction_COUNTS_INTERGENIC_and_the_EM_only_view_does_not():
    """``n_intergenic`` counts fragments that reached no locus; they never enter the EM, so ``gdna_em_count`` excludes
    them. Off capture gDNA is genome-uniform and more than half of it is intergenic, so a library gDNA fraction that
    drops them from the numerator while RNA stays in the denominator reads far low against truth. The deliverable
    includes the intergenic fragments, the EM-only diagnostic does not, and the two must differ.
    """
    mrna, gdna_em, intergenic = 5_000_000.0, 2_330_992.0, 2_601_271.0
    quant = pd.DataFrame({"transcript_id": ["a"], "count": [mrna]})
    result = _FakeResult(synthetic_count=0.0, gdna_em_count=gdna_em, n_intergenic=intergenic)
    summary = {"origin_counts": {"gdna": 5_000_000, "mrna": 5_000_000, "nrna": 0}}
    lib = QA.score_library(result, quant, summary)

    # The expectations are derived from the fixture's own inputs, never pasted from a run.
    want_deliverable = (gdna_em + intergenic) / (mrna + gdna_em + intergenic)
    want_em_only = gdna_em / (mrna + gdna_em)
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
    same = QA.score_library(
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
    mrna, gdna_em, intergenic = 5_000.0, 2_000.0, 3_000.0
    quant = pd.DataFrame({"transcript_id": ["a"], "count": [mrna]})
    summary = {"origin_counts": {"gdna": 4_900, "mrna": 5_100, "nrna": 0}}
    lib = QA.score_library(_FakeResult(0.0, gdna_em, n_intergenic=intergenic), quant, summary)
    assert lib["gdna_frac_est"] == pytest.approx(
        (gdna_em + intergenic) / (mrna + gdna_em + intergenic)
    )
    assert "n_intergenic_rna" not in lib and lib["n_intergenic_gdna"] == intergenic
    assert lib["rna_est"] == mrna


def test_the_count_metrics_do_not_move_when_the_length_model_does():
    """The count rows must be free of the fragment-length model, or the instrument would price the length models too and
    neither number would be attributable. Changing one transcript's effective length by 100x must not move them.
    """
    truth = pd.DataFrame({"transcript_id": ["a", "b"], "observed_mrna_fragments": [100.0, 50.0]})
    quant = pd.DataFrame(
        {
            "transcript_id": ["a", "b"],
            "count": [120.0, 30.0],
            "count_em": [0.0, 0.0],
            "tpm": [500000.0, 500000.0],
            "effective_length": [100.0, 100.0],
        }
    )
    a = QA.score_transcripts(quant, truth)
    b = QA.score_transcripts(quant.assign(effective_length=[10000.0, 100.0]), truth)
    for field in ("count_abs_err", "count_net_err", "fp_mass", "fn_mass", "spearman"):
        assert a[field] == pytest.approx(b[field], nan_ok=True), field
    assert a["tpm_abs_err"] != pytest.approx(b["tpm_abs_err"]), (
        "the TPM row did not move when an effective length moved 100x — it is not reading the "
        "length model at all and its caption is wrong"
    )


# ── Gate 5: a missing condition is a refusal, never a quiet average ──────────────────────────────


def test_the_report_REFUSES_arms_with_different_row_sets(tmp_path):
    """A shard that died would leave surviving rows that aggregate into a stratum total reading like a complete panel,
    so arms with different row sets are refused."""
    import json

    full = tmp_path / "full.jsonl"
    short = tmp_path / "short.jsonl"
    rows = [
        {
            "arm": "a",
            "condition": f"gdna_g{g}_ss_0.99_nrna_none_capture_on",
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
            "nrna_est": 0.0,
            "nrna_true": 0.0,
            "mrna_est": 1.0,
            "mrna_true": 1.0,
            "n_intergenic_gdna": 0.0,
        }
        for g in ("01", "50")
        for ax in ("transcript", "library")
    ]
    full.write_text("\n".join(json.dumps(r) for r in rows) + "\n")
    short.write_text("\n".join(json.dumps(r) for r in rows[:2]) + "\n")
    QA.report([full])  # one arm alone is fine
    with pytest.raises(SystemExit, match="different row set"):
        QA.report([full, short])


# ── Gate 6: a missing oracle cache is a refusal, never an invented truth ─────────────────────────


def test_a_missing_or_stale_oracle_cache_ABORTS(toy, tmp_path):
    """The oracle arms exist to inject truth, so falling back to anything else (a rescan under a different scan config,
    an empty payload) would inject something else under the name ``oracle``. ``read_scan_cache`` refuses a payload whose
    ``reach_digest`` does not describe this index; this checks the refusal propagates. Sum-to-full is checked against
    the condition's live scan cache, so partitions of another scan are refused.
    """
    import _oracle_arms

    from rigel.pipeline import scan_and_buffer
    from rigel.scan_cache import write_scan_cache

    bam, cfg = str(toy.bam_path), PipelineConfig()
    with pytest.raises(SystemExit, match="no valid scan or oracle cache"):
        QA.load_oracle(bam, toy.index, cfg, tmp_path, tmp_path, "absent")

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
    QA.load_oracle(bam, toy.index, cfg, tmp_path, orc, "t")  # the live scan's partitions load
    payload.region_contained_count[0, 0] += 1
    cache_scan()
    with pytest.raises(AssertionError, match="oracle INVALID"):
        QA.load_oracle(bam, toy.index, cfg, tmp_path, orc, "t")


def test_an_oracle_part_records_the_BAM_it_was_SPLIT_from(toy, tmp_path):
    """An oracle part is scanned from a split file in a scratch directory named for its condition, which the next build
    of a same-named condition, in any panel, overwrites. The part is a function of the BAM it was split from, so its
    cache records that BAM, and a reused split file is refused.
    """
    import _oracle_arms

    from rigel.pipeline import scan_and_buffer
    from rigel.scan_cache import read_scan_cache

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
    """One full ``run_pipeline`` with the EM's effective lengths under ``install_ruler_arm``."""
    from rigel.pipeline import run_pipeline

    cfg = PipelineConfig()
    cfg = dataclasses.replace(cfg, em=dataclasses.replace(cfg.em, seed=seed))
    restore, fired = QA.install_ruler_arm(arm, factor)
    try:
        result = run_pipeline(str(toy.bam_path), toy.index, cfg)
    finally:
        restore()
    return result.estimator.get_counts_df(toy.index), fired


def test_the_ruler_arm_hands_the_EM_the_lengths_it_names_and_its_noop_is_inert(toy):
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
    QA.check_ruler_arm("c", "oracle_ruler", fired)
    QA.check_ruler_arm("c", "oracle_ruler", dict(fired, max_abs_delta=0.0, truth_delta=0.0))
    with pytest.raises(RuntimeError, match="did not move"):
        QA.check_ruler_arm("c", "oracle_ruler", dict(fired, max_abs_delta=0.0))
    with pytest.raises(RuntimeError, match="NOOP"):
        QA.check_ruler_arm("c", "oracle_ruler_noop", dict(fired, max_abs_delta=1e-9))
    with pytest.raises(RuntimeError, match="never ran"):
        QA.check_ruler_arm("c", "oracle_ruler", dict(fired, ruler=0))


def test_the_capture_truth_is_anchored_on_the_fully_probed_class():
    """The shipped ruler's convention puts a fully probed transcript at 1 and every other below it, so the simulator's
    factor is divided by the fully probed class's median; a transcript with no length is the identity; and with no
    probed class (capture-OFF) every factor is exactly 1. The fixture's median over every transcript differs from the
    probed class's, so an anchor read off the wrong pool shows."""
    factor = np.array([10.0, 10.0, 12.0, 1.0, 2.0, 3.0, np.nan])
    probed = np.array([1.0, 0.95, 0.92, 0.0, 0.1, 0.2, 0.0])
    plain = np.array([100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 0.0])
    got = QA.anchor_capture_factor(factor, probed, plain)
    np.testing.assert_allclose(got, [1.0, 1.0, 1.2, 0.1, 0.2, 0.3, 1.0], rtol=1e-12)
    off = QA.anchor_capture_factor(np.full(4, 3.0), np.zeros(4), np.full(4, 100.0))
    assert np.array_equal(off, np.ones(4))


# ── Gate 8: arms scored under different assignment modes are refused, never compared ──────────────


def test_the_report_REFUSES_arms_scored_under_different_assignment_modes(tmp_path):
    """``--set em.assignment_mode=fractional`` changes every count, so a fractional arm beside a sampled one would
    report the mode's effect as the arm's. Every row carries its mode and the report refuses a mix; a row without the
    stamp counts as a different mode."""
    row = dict(
        condition="gdna_g05_ss_0.99_nrna_mid_capture_on", axis="transcript", count_abs_err=1.0
    )
    a, b = tmp_path / "a.jsonl", tmp_path / "b.jsonl"
    a.write_text(json.dumps(dict(row, assignment_mode="fractional")) + "\n")
    b.write_text(json.dumps(dict(row, assignment_mode="sample")) + "\n")
    with pytest.raises(SystemExit, match="assignment mode"):
        QA.report([a, b])
    b.write_text(json.dumps(row) + "\n")
    with pytest.raises(SystemExit, match="assignment mode"):
        QA.report([a, b])


# ── the markdown report: what a release is judged on ────────────────────────────────────────────
#
# The report renders numbers the arm files already hold, so every way of getting it wrong produces a plausible-looking
# table: the wrong gDNA column, a ratio invented where the truth is zero, a pool folded into another. These tests cover
# those.


def _arm_rows(tmp_path, **over):
    """One condition, three axes, written as an arm jsonl, with numbers chosen so every quantity the report derives is
    distinguishable from every other."""
    cond = "gdna_g50_ss_0.99_nrna_mid_capture_off"
    lib = {
        "condition": cond,
        "axis": "library",
        "arm": "a",
        "assignment_mode": "fractional",
        "mrna_est": 4_000.0,
        "mrna_true": 4_100.0,
        "nrna_est": 900.0,
        "nrna_true": 1_000.0,
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
    qa = _load_sibling("quant_accuracy.py")
    src, _cond = _arm_rows(tmp_path, **over)
    out = tmp_path / "report.md"
    qa.markdown_report(src, out)
    return out.read_text()


def test_the_report_scores_gDNA_INCLUDING_the_intergenic_fragments(tmp_path):
    """``gdna_est`` is the EM's number and excludes fragments that reached no locus, while ``gdna_true`` counts every
    gDNA fragment. Off capture the intergenic share is the larger one, so printing ``gdna_est`` against that truth would
    read as a large under-call that is not there. The comparable estimate is ``gdna_est + n_intergenic_gdna``.
    """
    md = _render(tmp_path)
    assert "5,100" in md, "the gDNA estimate did not include the intergenic fragments"
    assert "| 2,000 |" not in md, "the report printed the EM's gDNA alone against the full truth"
    assert "+4.1 %" in md, f"expected (5,100-4,900)/4,900; got:\n{md}"


def test_a_percent_error_is_REFUSED_where_the_truth_is_zero(tmp_path):
    """The `g00` rung's true gDNA count is exactly 0, and a percent error there would invent a scale. The raw count is
    the whole answer; the ratio is `n/a`."""
    md = _render(tmp_path, gdna_true=0.0)
    assert "n/a" in md
    assert "inf" not in md.lower() and "nan %" not in md.lower()


def test_the_three_pools_are_reported_APART_and_the_synthetic_one_by_name(tmp_path):
    """The split is on ``is_synthetic``, and the report must say so: a single-exon annotated transcript
    (``is_single_exon``) belongs to the annotated pool. Folding the two RNA pools together would hide the false-positive
    channel the nascent row exposes."""
    md = _render(tmp_path)
    assert "SYNTHETIC entities only" in md
    assert "is_synthetic`, never `is_single_exon" in md
    for n in ("900", "1,000", "4,000", "4,100"):  # nascent est/true, annotated est/true
        assert n in md, f"pool figure {n} missing"


def test_the_transcript_and_gene_sections_score_the_ANNOTATED_POOL_ALONE(tmp_path):
    """The transcript table a user reads holds no synthetic entity (``get_counts_df`` drops them), so the transcript and
    gene truth must equal the library's annotated truth, not the RNA total; scoring against ``rna_true`` would credit
    the tool for nascent mass not in the table."""
    md = _render(tmp_path)
    qa = _load_sibling("quant_accuracy.py")
    src, cond = _arm_rows(tmp_path)
    rows = qa._load(src)
    assert rows[(cond, "transcript")]["count_true"] == rows[(cond, "library")]["mrna_true"]
    assert rows[(cond, "transcript")]["count_true"] != rows[(cond, "library")]["rna_true"]
    assert "10.00 %" in md, "Σ|Δ| as a share of the annotated truth (410/4,100) is missing"


def test_every_table_in_the_report_is_RECTANGULAR(tmp_path):
    """A ragged markdown table renders as a wrong table rather than an error, and the escaped pipes in a `Σ|Δ|` header
    are cell content, not delimiters."""
    md = _render(tmp_path)
    for block in md.split("\n\n"):
        lines = [x for x in block.splitlines() if x.startswith("|")]
        if len(lines) < 2:
            continue
        widths = {len(x.replace("\\|", "").split("|")) for x in lines}
        assert len(widths) == 1, f"ragged table {widths}: {lines[0][:80]}"


def test_the_DEFERRED_stratum_is_MARKED_wherever_it_appears(tmp_path):
    """Unstranded × capture-ON is reported on every benchmark, is never a development target, and carries most of the
    error, so a reader who misses the mark would read a pooled total as the tool's accuracy."""
    qa = _load_sibling("quant_accuracy.py")
    lib, tx, gene = (json.loads(x) for x in _arm_rows(tmp_path)[0].read_text().splitlines())
    for r in (lib, tx, gene):
        r["condition"] = "gdna_g50_ss_0.50_nrna_mid_capture_on"
    src = tmp_path / "def.jsonl"
    src.write_text("\n".join(json.dumps(r) for r in (lib, tx, gene)) + "\n")
    out = tmp_path / "d.md"
    qa.markdown_report(src, out)
    md = out.read_text()
    assert md.count("DEFERRED") >= 4, "the deferred stratum is unmarked in at least one section"


def test_the_report_does_NOT_present_a_column_that_cannot_fire(tmp_path):
    """``fn_mass`` needs an estimate of exactly zero, which a posterior under fractional assignment essentially never
    is, so it reads 0.0 everywhere, and a column of zeros labelled "false-negative mass" would read as a perfect score
    for something never measured. ``count_under`` is the live quantity and is reported instead.
    """
    md = _render(tmp_path)
    assert "false-negative" not in md.split("## 2.")[1], "a column that cannot fire was reported"
    assert "under-assigned" in md and "over-assigned" in md
    assert "255" in md, "count_under (the live under-assignment) is missing"


def test_the_report_does_NOT_call_the_truth_table_s_ROW_COUNT_a_transcript_count(tmp_path):
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


def test_every_gDNA_RUNG_is_labelled_from_its_truth_not_from_a_map_of_known_rungs(tmp_path):
    """The level column is read off the row's own truth, ``gdna_frac_true``, so a panel with any rung (the test
    chromosome's ``g25``, the depth panels' ``g001`` and ``g01``) gets a report. A hand map of the ladder's rung tokens
    would raise on ``g25``, and digits read from the name would print ``g001`` as 1 %.
    """
    rungs = {"g00": (0.0, "0 %"), "g001": (0.001, "0.1 %"), "g05": (0.0500004, "5 %"),
             "g25": (0.25, "25 %"), "g98": (0.98, "98 %")}  # fmt: skip
    conds = {f"gdna_{k}_ss_0.99_nrna_mid_capture_off": f for k, (f, _label) in rungs.items()}
    out = tmp_path / "rungs.md"
    QA.markdown_report(_arm_file(tmp_path, conds), out)
    md = out.read_text()
    for (k, (_f, label)), cond in zip(rungs.items(), conds):
        assert f"| `{cond}` | {label} |" in md, f"rung {k} is not labelled {label}"


def test_a_strand_specificity_in_NEITHER_half_is_reported_APART_never_pooled(tmp_path):
    """ss 0.50 is the unstranded half and ss 0.99 the stranded half, judged against different bars. The test panel's ss
    0.70 is neither, so it is reported apart rather than moving the half a policy must win on. One definition
    (`_shared.strandedness`) serves every instrument, and a name that spells no specificity raises rather than landing
    in a half.
    """
    shared = importlib.import_module("_shared")
    assert shared.strandedness("gdna_g50_ss_0.50_nrna_mid_capture_on") == "unstranded"
    assert shared.strandedness("gdna_g50_ss_0.99_nrna_mid_capture_on") == "stranded"
    assert shared.strandedness("gdna_g50_ss_0.70_nrna_file_capture_on") == "ss 0.70"
    with pytest.raises(ValueError):
        shared.strandedness("gdna_g50_nrna_mid_capture_on")

    c50, c70 = "gdna_g50_ss_0.50_nrna_file_capture_off", "gdna_g50_ss_0.70_nrna_file_capture_off"
    assert shared.strata([c50, c70]) == [*shared.STRATA, ("ss 0.70", "capture OFF")]
    out = tmp_path / "apart.md"
    QA.markdown_report(_arm_file(tmp_path, {c50: 0.5, c70: 0.5}), out)
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
    synthetic span's fragments are ``observed_nrna_fragments`` on its own row, and its ``observed_mrna_fragments`` is
    zero because a synthetic span emits no mature RNA. A single-exon annotated transcript carries both columns."""
    return pd.DataFrame(
        {
            "transcript_id": ["t_spliced", "t_single_exon", "nrna_gene_a", "nrna_gene_b"],
            "observed_mrna_fragments": [900.0, 100.0, 0.0, 0.0],
            "observed_nrna_fragments": [0.0, 7.0, 300.0, 50.0],
        }
    )


class _AllocIndex:
    """The three columns ``truth_weights`` reads off an index."""

    def __init__(self, t_ids):
        self.t_df = pd.DataFrame(
            {"t_id": list(t_ids), "t_index": np.arange(len(t_ids), dtype=np.int64)}
        )
        self.num_transcripts = len(t_ids)


def test_the_allocation_arm_does_NOT_hand_a_weight_of_ZERO_to_an_RNA_EMITTING_transcript():
    """The arm must weigh every RNA component, and a synthetic span is one.

    ``--arm oracle_alloc_seed`` claims to hand the EM the true per-transcript allocation. Reading a mature column alone
    would give every synthetic nascent entity a weight of zero, which is not a true allocation but a suppression of the
    component, and the arm would then appear to remove the nascent over-call. A transcript the simulator gave RNA
    fragments must receive a positive weight, and one it gave none must receive zero.
    """
    truth = _alloc_truth_table()
    index = _AllocIndex(truth["transcript_id"])
    w = QA.truth_weights(truth, index)

    emits = (truth["observed_mrna_fragments"] + truth["observed_nrna_fragments"]).to_numpy()
    assert np.all(w[emits > 0] > 0.0), (
        "a transcript the simulator gave RNA fragments received a weight of ZERO — a zero weight is "
        "a suppression, not an allocation"
    )
    assert np.all(w[emits == 0] == 0.0), "a silent transcript received prior mass"
    # the weights are the RNA the simulator realised, mature and nascent alike
    np.testing.assert_allclose(w, [900.0, 107.0, 300.0, 50.0], rtol=1e-12)


def test_the_allocation_weights_REFUSE_a_truth_table_that_carries_no_nascent_column():
    """A truth table with a mature column and no nascent one cannot express a synthetic span's allocation, so the arm
    refuses rather than weighting every synthetic span at zero. The suite-level molar table is refused too: its
    ``mrna_abundance`` / ``nrna_abundance`` pair counts molecules, not the fragments a condition realised."""
    truth = _alloc_truth_table().drop(columns=["observed_nrna_fragments"])
    with pytest.raises(SystemExit, match="nascent"):
        QA.truth_weights(truth, _AllocIndex(truth["transcript_id"]))
    molar = _alloc_truth_table().set_axis(
        ["transcript_id", "mrna_abundance", "nrna_abundance"], axis=1
    )
    with pytest.raises(SystemExit, match="nascent"):
        QA.truth_weights(molar, _AllocIndex(molar["transcript_id"]))


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

    counts = {"gdna": 0, "mrna": 0, "nrna": 0}
    for name, _recs in groups:
        counts[parse_origin(name).kind] += 1
    return counts


def test_the_aligner_loss_counts_UNMAPPED_TEMPLATES_by_their_true_origin(toy, tmp_path):
    """An aligned input is scored against the simulator's truth, which counts every fragment it drew, so a row read from
    it must carry what the aligner left out, by origin: on the ladder's suite genome part of the gDNA is drawn over N
    runs no aligner can place, and read against the full truth that would look like Rigel's error. A template is lost
    only when none of its records is mapped; one placed mate is the aligner's delivery and Rigel's to judge. Counting
    records instead of templates, or reading the origin from anything but the read name, would move the counts.
    """
    from rigel.sim.read_name import parse_origin

    header, groups = _templates(toy.bam_path)
    expected = {"gdna": 0, "mrna": 0, "nrna": 0}
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
    got = QA.aligner_loss(bam, _origin_counts(groups))
    assert got == {f"unaligned_{k}": v for k, v in expected.items()}


def test_an_aligned_BAM_that_is_not_the_simulators_is_REFUSED(toy, tmp_path):
    """The aligned input must hold the simulator's templates, each once and in one group: a missing template (a
    truncated alignment) or one split into two runs of records (an input that is not collated) would shift the aligner's
    losses silently."""
    header, groups = _templates(toy.bam_path)
    counts = _origin_counts(groups)
    short = _write(tmp_path / "short.bam", header, groups[:-1])
    with pytest.raises(SystemExit, match="not the simulator's"):
        QA.aligner_loss(short, counts)
    name, recs = next((n, r) for n, r in groups if len(r) >= 2)
    split = [g for g in groups if g[0] != name] + [(name, recs[:1])]
    split.insert(0, (name, recs[1:]))
    with pytest.raises(SystemExit, match="not the simulator's"):
        QA.aligner_loss(_write(tmp_path / "split.bam", header, split), counts)


def test_another_input_BAM_is_scored_by_the_BASE_arm_alone_and_must_exist(tmp_path):
    """The oracle arms inject truth built from the simulator's own alignments (the oracle cache), so on another input
    their substitution would be another scan's and is refused. A condition without the named BAM is refused by name
    before anything runs, never scored as a shorter panel."""
    for c in ("c1", "c2"):
        (tmp_path / c).mkdir()
        (tmp_path / c / QA.ORACLE_BAM).write_bytes(b"")
    (tmp_path / "c1" / "aligned.bam").write_bytes(b"")
    with pytest.raises(SystemExit, match="c2"):
        QA.check_input_bam(tmp_path, ["c1", "c2"], "base", "aligned.bam")
    QA.check_input_bam(tmp_path, ["c1"], "base", "aligned.bam")  # every condition has it
    with pytest.raises(SystemExit, match="base"):
        QA.check_input_bam(tmp_path, ["c1"], "oracle", "aligned.bam")
    QA.check_input_bam(tmp_path, ["c1", "c2"], "oracle", QA.ORACLE_BAM)  # the simulator's: any arm
