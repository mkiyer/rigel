"""Falsification gates for ``scripts/design/quant_accuracy.py`` — the end-to-end accuracy table and
the prior-injection ceiling above it.

The instrument answers two questions with one scorer: how accurate is the tool, and how much of that
error would a PERFECT calibration prior remove. Both answers are differences between an estimate and a
truth column, so every way of getting them wrong produces a *plausible* number — the wrong truth
column, an injection that never fired, a false-positive channel that is invisible because the rows are
filtered out before scoring. These gates are those ways.

⛔ **EACH GATE CARRIES ITS OWN PERTURBATION.** A gate that has never been watched to fire has not been
written yet.
"""

from __future__ import annotations

import dataclasses
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from rigel.calibration.priors import LocusPriors
from rigel.config import PipelineConfig
from rigel.sim import GDNAConfig, ReadSimConfig, Scenario

_MODULES: dict = {}


def _load_sibling(name: str):
    """Import a ``scripts/design/`` instrument by path. ⚠ Registered in ``sys.modules`` BEFORE
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
    wd = tmp_path_factory.mktemp("qa_sim")
    sc = Scenario("qa", genome_length=9000, seed=17, work_dir=wd)
    sc.add_gene(
        "g1",
        "+",
        [
            {"t_id": "t1", "exons": [(600, 1100), (1800, 2300)], "abundance": 60},
            {"t_id": "t1b", "exons": [(1000, 1100), (1800, 1900)], "abundance": 30},
        ],
    )
    sc.add_gene("g2", "-", [{"t_id": "t2", "exons": [(4000, 4500), (5200, 5700)], "abundance": 40}])
    sc.add_gene(
        "g3",
        "+",
        [
            {"t_id": "t3", "exons": [(6600, 7100), (7800, 8300)], "abundance": 25},
            {"t_id": "t3b", "exons": [(6600, 7100), (7400, 7440), (7800, 8300)], "abundance": 12},
        ],
    )
    return sc.build_oracle(
        n_rna_fragments=4000,
        gdna_fraction=1.0,
        nrna_abundance=15.0,
        sim_config=ReadSimConfig(
            frag_mean=180,
            frag_std=30,
            frag_min=80,
            frag_max=400,
            read_length=90,
            strand_specificity=0.99,
            seed=17,
        ),
        gdna_config=GDNAConfig(abundance=0.0, frag_mean=240, frag_std=45),
    )


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


def _quant(toy, arm, oracle, seed=QA.DEFAULT_EM_SEED):
    """One full ``run_pipeline`` under one arm; returns ``(counts_df, fired)``.

    ⛔ **The seed is PINNED and that is not tidiness.** ``EMConfig.seed`` ships as ``None`` with
    ``assignment_mode="sample"``, so the EM's hard assignment is an unseeded categorical draw: two
    back-to-back runs of the IDENTICAL pipeline differ here on 4 transcripts by up to 43 fragments
    (measured). Without the pin the byte-identity gate below could never pass and its failure would
    say nothing about the injection.
    """
    from rigel.pipeline import run_pipeline

    cfg = PipelineConfig()
    cfg = dataclasses.replace(cfg, em=dataclasses.replace(cfg.em, seed=seed))
    restore, fired = QA.install_arm(arm, oracle)
    try:
        result = run_pipeline(str(toy.bam_path), toy.index, cfg)
    finally:
        restore()
    return result.estimator.get_counts_df(toy.index), fired


# ── GATE 1: noop is byte-identical END TO END, and the wrapper really ran ────────────────────────


def test_the_noop_arm_reproduces_BASE_byte_identically_through_the_whole_pipeline(toy, toy_oracle):
    """⛔ TRAPS: byte-identity-gate, and this is the harness's own falsification. The injection reaches the EM, whose
    per-locus draw is order-sensitive; if the wrapper perturbed anything at all — even the order a
    float was summed in — every ``oracle`` number would be that perturbation plus the effect.

    ⭐ **``noop`` is not an early return.** It builds O in full and discards it, so this gate covers
    the whole wrapper path and not an ``if`` (TRAPS: could-the-arm-have-fired). The perturbation makes it take one field
    and asserts the result stops matching — proving the comparison could have failed.
    """
    base, base_fired = _quant(toy, "base", None)
    noop, noop_fired = _quant(toy, "noop", toy_oracle)
    assert base_fired["n"] >= 1 and noop_fired["n"] >= 1

    # ⭐⭐ TWO STANDARDS, AND THEY ARE THE NUMERIC CONVENTION. Every COUNT column is integer-derived and
    # must match EXACTLY — measured over 8 shipped multi-threaded runs, their run-to-run spread is
    # 0.000e+00, because integer addition is associative. `posterior_mean` is float-derived and wanders
    # by **1.503e-15** (~7 ulp) between identical runs, because the conserved-mass banks are float64 and
    # the per-worker merge order is a data-dependent race.
    #
    # ⛔ The tolerance is DERIVED and GATED, not chosen: it sits ~6 orders above that measured spread
    # (headroom for a library far deeper than this toy, since re-association scales with the number of
    # additions) and ~9 orders below the injection it must catch — and the perturbation arm below proves
    # it is still tight enough, because taking ONE oracle field must break the comparison.
    # ⚠ Owner ruling 2026-08-11: the tool is not bit-reproducible and tests validate within a tolerance.
    exact = [c for c in base.columns if c != "posterior_mean"]
    pd.testing.assert_frame_equal(base[exact], noop[exact], check_exact=True)
    pd.testing.assert_series_equal(
        base["posterior_mean"], noop["posterior_mean"], check_exact=False, rtol=1e-9, atol=0.0
    )

    saved = QA._ARM_FIELDS["noop"]
    try:
        QA._ARM_FIELDS["noop"] = ("gdna_prior_count",)
        perturbed, _ = _quant(toy, "noop", toy_oracle)
    finally:
        QA._ARM_FIELDS["noop"] = saved
    assert not perturbed["count"].equals(base["count"]), (
        "taking the oracle's gDNA prior changed no transcript count on a contaminated toy — the "
        "injection does not reach the EM and every ceiling this script reports would be zero"
    )


def test_the_UNSEEDED_shipped_config_is_not_reproducible_and_that_is_why_the_seed_is_pinned(toy):
    """⛔⛔ **A PROPERTY OF THE SHIPPED TOOL, PINNED HERE BECAUSE EVERY A/B ABOVE DEPENDS ON IT.**
    ``EMConfig.seed`` defaults to ``None`` and ``assignment_mode`` to ``"sample"``, so two runs of the
    identical pipeline on the identical BAM return different transcript counts. Any end-to-end arm
    comparison run on the default config therefore reports sampling noise on top of its effect, and a
    byte-identity ``noop`` is impossible.

    ⭐ This asserts BOTH halves: unseeded runs differ, and the same runs with a pinned seed do not. If
    the default ever becomes deterministic the first half fails — which is a result to record, not a
    test to widen (TRAPS: waive-with-a-measurement).
    """
    from rigel.pipeline import run_pipeline

    unseeded = [
        run_pipeline(str(toy.bam_path), toy.index, PipelineConfig())
        .estimator.get_counts_df(toy.index)["count"]
        .to_numpy()
        for _ in range(2)
    ]
    assert not np.array_equal(unseeded[0], unseeded[1]), (
        "the unseeded pipeline reproduced itself exactly — if EMConfig.seed no longer defaults to "
        "None this gate has become a measurement of the new default and must be re-recorded"
    )
    a, _ = _quant(toy, "base", None)
    b, _ = _quant(toy, "base", None)
    assert np.array_equal(a["count"].to_numpy(), b["count"].to_numpy()), (
        "pinning the seed did not make the pipeline reproducible — no arm comparison here is "
        "attributable"
    )


def test_an_injection_that_never_fires_RAISES_rather_than_reporting_no_effect(
    tmp_path, monkeypatch
):
    """⛔ TRAPS: an-ablation-that-never-ran. "The perfect prior changed nothing" is a publishable conclusion, so the
    one thing that must not be able to produce it is a wrapper that was never called — which happens
    for real when ``quant_from_buffer`` returns early on an empty unit set.

    ⚠ ``run_pipeline`` is stubbed to a no-op, which is exactly the shape of that early return: the
    pipeline completes, the wrapper is installed, and ``assemble_priors`` is never reached.
    """
    cond = tmp_path / "c1"
    cond.mkdir()
    pd.DataFrame({"transcript_id": ["a"], "mrna_abundance": [1.0]}).to_csv(
        cond / "truth_abundances.tsv", sep="\t", index=False
    )
    (cond / "truth_summary.json").write_text('{"origin_counts": {"gdna": 0, "mrna": 1, "nrna": 0}}')
    monkeypatch.setattr(QA, "run_pipeline", lambda *a, **k: None)
    monkeypatch.setattr(QA, "load_oracle", lambda *a, **k: object())
    with pytest.raises(RuntimeError, match="never wrapped-and-called"):
        QA.run_condition("oracle", tmp_path, None, "c1", PipelineConfig(), tmp_path)


# ── GATE 2: each single-array arm replaces exactly its own field ─────────────────────────────────


def test_each_oracle_arm_replaces_ITS_field_and_leaves_the_other_two_shipped(toy, toy_oracle):
    """⭐ The three single-array arms exist to say WHICH of the prior's three numbers carries the
    value, and that only works if each one is surgical. This drives the wrapper directly with two
    recognisably different priors and checks, field by field, that exactly the named ones moved.

    ⛔ The perturbation is the whole point of the loop: every arm is checked against every field, so
    an arm that quietly replaced all three (or none) fails on the fields it should not have touched.
    """
    import rigel.calibration.priors as PRIORS
    from rigel.calibration.region_arrays import RegionArrays

    ra = RegionArrays.from_index(toy.index)
    shipped = LocusPriors(np.array([1.0, 2.0]), np.array([3.0, 4.0]), np.array([5.0, 6.0]))
    oracle_p = LocusPriors(np.array([10.0, 20.0]), np.array([30.0, 40.0]), np.array([50.0, 60.0]))
    calls = {"n": 0}

    def fake_assemble(cal, ra, ml):
        calls["n"] += 1
        return oracle_p if calls["n"] % 2 == 0 else shipped

    original = PRIORS.assemble_priors
    PRIORS.assemble_priors = fake_assemble
    try:
        for arm, fields in QA._ARM_FIELDS.items():
            if arm == "base":
                continue
            calls["n"] = 0
            restore, fired = QA.install_arm(arm, toy_oracle)
            try:
                out = PRIORS.assemble_priors(_FakeCal(), ra, [])
            finally:
                restore()
            assert fired["n"] == 1
            for f in ("gdna_prior_count", "rna_prior_count", "gdna_eff_len"):
                want = oracle_p if f in fields else shipped
                assert np.array_equal(getattr(out, f), getattr(want, f)), (
                    f"arm {arm!r}: field {f} came from the wrong prior"
                )
    finally:
        PRIORS.assemble_priors = original


@dataclasses.dataclass
class _FakeCal:
    """Enough of a ``CalibrationResult`` for ``dataclasses.replace`` inside the wrapper."""

    mass_gdna_region: object = None
    mass_rna_region: object = None
    mass_gdna_edge: object = None
    mass_rna_edge: object = None
    mass_rna_spliced_edge: object = None
    mass_rna_junction: object = None


# ── GATE 3: the truth column is the REALISED fragment count, not the molar abundance ─────────────


def test_the_scorer_reads_the_realised_fragment_column_and_a_molar_one_would_be_visible():
    """⛔ The panel ships TWO truth tables and they are different quantities. The per-condition
    ``truth_abundances.tsv`` gives the realised OBSERVED fragment count; the suite-level
    ``truth_abundances_nrna_none.tsv`` gives the pre-capture MOLAR abundance. Scoring a fragment
    estimate against a molar truth charges the tool for hybrid capture, which it never claims to
    invert.

    ⭐ The perturbation feeds the scorer a molar-scaled truth (the same shape, a different unit) and
    asserts the error explodes — so a future edit that reaches for the wrong file cannot pass quietly.
    """
    truth = pd.DataFrame(
        {
            "transcript_id": ["a", "b", "c"],
            "mrna_abundance": [1000.0, 500.0, 0.0],
        }
    )
    quant = pd.DataFrame(
        {
            "transcript_id": ["a", "b", "c"],
            "count": [1000.0, 500.0, 0.0],
            "count_em": [0.0, 0.0, 0.0],
            "tpm": [0.0, 0.0, 0.0],
            "effective_length": [100.0, 100.0, 100.0],
        }
    )
    exact = QA.score_transcripts(quant, truth)
    assert exact["count_abs_err"] == pytest.approx(0.0)
    assert exact["fp_mass"] == pytest.approx(0.0)

    molar = truth.assign(mrna_abundance=truth["mrna_abundance"] / 1500.0 * 4.6)
    wrong = QA.score_transcripts(quant, molar)
    assert wrong["count_abs_err"] > 1000.0, (
        "a molar truth scored as a fragment count produced almost no error — the two units are not "
        "distinguishable to this scorer and the wrong file could be read silently"
    )


# ── GATE 4: the two false-positive channels are separated, and neither hides ─────────────────────


def test_false_positive_and_false_negative_mass_are_separate_and_reconcile():
    """⭐ ``fp_mass`` (estimate on a transcript the simulator gave zero fragments) and ``fn_mass``
    (real RNA assigned nowhere) are different failures with different fixes, and a single ``Σ|Δ|``
    cannot separate them. ⛔ The perturbation flips one transcript from silent to expressed and
    asserts its mass MOVES OUT of ``fp_mass`` — a scorer that keyed on the estimate rather than on
    the truth would keep it there.
    """
    truth = pd.DataFrame({"transcript_id": ["a", "b", "c"], "mrna_abundance": [100.0, 0.0, 50.0]})
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

    lifted = truth.assign(mrna_abundance=[100.0, 25.0, 50.0])
    s2 = QA.score_transcripts(quant, lifted)
    assert s2["fp_mass"] == pytest.approx(0.0), (
        "a transcript the truth now expresses is still counted as a false positive — fp_mass is "
        "keyed on the estimate, not on the truth"
    )


def test_the_GENE_row_collapses_isoform_ambiguity_and_still_sees_real_gene_error():
    """⭐⭐ The gene row exists to answer "is the surviving transcript error isoform ambiguity, or
    something else?" — and it can only answer that if swapping mass BETWEEN a gene's own isoforms is
    invisible to it while moving mass OUT of the gene is not. Both halves are asserted here, because a
    gene scorer that saw neither would also report zero on the first case.
    """
    truth = pd.DataFrame(
        {
            "transcript_id": ["a1", "a2", "b1"],
            "gene_id": ["A", "A", "B"],
            "mrna_abundance": [100.0, 0.0, 50.0],
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
    """⛔⛔ **THIS IS A DEFECT THAT SHIPPED IN THIS FILE AND WAS CAUGHT BY A SANITY CHECK, NOT BY A
    GATE.** ``n_intergenic`` counts fragments that reached no locus; they never enter the EM, so
    ``gdna_em_count`` excludes them. Off capture, gDNA is genome-uniform and MORE THAN HALF of it is
    intergenic — measured at `g50 ss0.50 capture_off`, 2,601,271 intergenic against 2,330,992
    EM-assigned. The first version of this scorer left them out of the numerator while RNA stayed in
    the denominator: it read **0.3151 against a truth of 0.5000** and fabricated a "systematic
    off-capture EM under-call" that does not exist. With them it reads **0.4933**.

    ⭐ The fixture reproduces that exact shape, and asserts BOTH views: the deliverable includes the
    intergenic fragments, the EM-only diagnostic does not, and they must not be equal — otherwise one
    of the two is not being computed.
    """
    mrna, gdna_em, intergenic = 5_000_000.0, 2_330_992.0, 2_601_271.0
    quant = pd.DataFrame({"transcript_id": ["a"], "count": [mrna]})
    result = _FakeResult(nrna_em_count=0.0, gdna_em_count=gdna_em, n_intergenic=intergenic)
    summary = {"origin_counts": {"gdna": 5_000_000, "mrna": 5_000_000, "nrna": 0}}
    lib = QA.score_library(result, quant, summary)

    # ⚠ the expectations are DERIVED from the fixture's own inputs, not pasted from a run. The first
    # version pasted 0.4932 off the real `g50` condition, whose mRNA total is not this one's, and
    # failed for a reason that had nothing to do with the behaviour under test.
    want_deliverable = (gdna_em + intergenic) / (mrna + gdna_em + intergenic)
    want_em_only = gdna_em / (mrna + gdna_em)
    assert lib["gdna_frac_true"] == pytest.approx(0.5)
    assert lib["gdna_frac_est"] == pytest.approx(want_deliverable), (
        "the deliverable fraction is not counting intergenic fragments as gDNA"
    )
    assert lib["gdna_em_frac_est"] == pytest.approx(want_em_only)
    # and the defect's own signature: the EM-only view reads far LOW against a truth of 0.5
    assert want_em_only < 0.33 and want_deliverable > 0.49
    assert lib["gdna_frac_est"] != pytest.approx(lib["gdna_em_frac_est"]), (
        "the two views are equal, so one of them is not being computed — the whole point is that "
        "they differ by exactly the population the EM never sees"
    )

    # the perturbation: with NO intergenic fragments the two views must coincide exactly
    same = QA.score_library(
        _FakeResult(nrna_em_count=0.0, gdna_em_count=2_330_992.0, n_intergenic=0.0), quant, summary
    )
    assert same["gdna_frac_est"] == pytest.approx(same["gdna_em_frac_est"]), (
        "with zero intergenic fragments the two denominators are identical and the two numbers must "
        "be too — they are not, so one of them has a second difference in it"
    )


@dataclasses.dataclass
class _FakeEstimator:
    nrna_em_count: float
    gdna_em_count: float


@dataclasses.dataclass
class _FakeStats:
    n_intergenic: float


class _FakeResult:
    """Enough of a ``PipelineResult`` for ``score_library`` — the three numbers it reads."""

    def __init__(self, nrna_em_count, gdna_em_count, n_intergenic):
        self.estimator = _FakeEstimator(nrna_em_count, gdna_em_count)
        self.stats = _FakeStats(n_intergenic)


def test_the_count_metrics_do_not_move_when_the_length_model_does():
    """⛔ TRAPS: price-the-halves-separately. The count rows must be free of the fragment-length model, or this
    instrument would silently price the length models too and neither number would be attributable.
    The perturbation changes ONE transcript's effective length by 100x."""
    truth = pd.DataFrame({"transcript_id": ["a", "b"], "mrna_abundance": [100.0, 50.0]})
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


# ── GATE 5: a missing condition is a refusal, never a quiet average ──────────────────────────────


def test_the_report_REFUSES_arms_with_different_row_sets(tmp_path):
    """⛔ TRAPS: byte-identity-gate's first recorded lie was an arm with ZERO rows scoring "32/32 IDENTICAL"
    because the comparison looped over one arm's keys. Here the equivalent is a shard that died: the
    surviving rows would aggregate into a stratum total that reads like a complete panel."""
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
            # ⭐ the POOL-LEVEL fields table ⑥ reads. They are part of the `library` row's real schema,
            # so a fixture without them was testing the report against a row shape that never ships.
            "gdna_est": 1.0,
            "gdna_true": 1.0,
            "nrna_est": 0.0,
            "nrna_true": 0.0,
            "mrna_est": 1.0,
            "mrna_true": 1.0,
            "n_intergenic": 0.0,
        }
        for g in ("01", "50")
        for ax in ("transcript", "library")
    ]
    full.write_text("\n".join(json.dumps(r) for r in rows) + "\n")
    short.write_text("\n".join(json.dumps(r) for r in rows[:2]) + "\n")
    QA.report([full])  # one arm alone is fine
    with pytest.raises(SystemExit, match="different row set"):
        QA.report([full, short])


# ── GATE 6: a missing oracle cache is a refusal, never an invented truth ─────────────────────────


def test_a_missing_or_stale_oracle_cache_ABORTS(toy, tmp_path):
    """⛔ The oracle arms exist to inject TRUTH. Falling back to anything else — a rescan under a
    different scan config, an empty payload — would inject something else under the name ``oracle``.
    ``read_scan_cache`` already refuses a payload whose ``reach_digest`` does not describe this index;
    this asserts the refusal is propagated rather than swallowed."""
    with pytest.raises(SystemExit, match="no valid oracle cache"):
        QA.load_oracle(str(toy.bam_path), toy.index, PipelineConfig(), tmp_path, "absent")
