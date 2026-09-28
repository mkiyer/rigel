"""Falsification gates for ``scripts/design/prior_vs_oracle.py`` — the instrument that scores
calibration's endpoint (``LocusPriors``) against the origin-split oracle.

Everything the instrument prints is a difference between two per-locus arrays of fragment counts,
so every way of getting it wrong is a way of getting a plausible number: a lever that never fired,
a reference that quietly became the arm, a weight that rewards the solver for declining to answer,
a projection that loses mass off the end of a locus. These gates are those ways, and each carries
its own perturbation in the same test (TRAPS: perturb-every-gate). The scenario is a
single-reference toy, which is enough for the scoring and plumbing they cover — arithmetic over two
per-locus arrays plus the override lever — and deliberately not enough to judge the deposit path,
which has its own gates in ``tests/native/`` and its truth-scored instruments on the panel.
"""

from __future__ import annotations

import dataclasses
import importlib.util
from pathlib import Path

import numpy as np
import pytest

from rigel.config import PipelineConfig
from _prior_toy import build_toy, build_toy_zero_gdna

_MODULES: dict = {}


def _load_sibling(name: str):
    """Import a ``scripts/design/`` instrument by path.

    ``scripts/`` is not a package and must not become one. The module is registered in
    ``sys.modules`` BEFORE execution: ``@dataclass`` resolves its own class's ``__module__`` through
    that table and fails at class-definition time otherwise.
    """
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


PV = _load_sibling("prior_vs_oracle.py")


# ── the toys ─────────────────────────────────────────────────────────────────────────────────────


@pytest.fixture(scope="module")
def toy(tmp_path_factory):
    return build_toy(tmp_path_factory)


@pytest.fixture(scope="module")
def toy_zero_gdna(tmp_path_factory):
    return build_toy_zero_gdna(tmp_path_factory)


def _measure(result, tmp_path, tag):
    return PV.measure_condition(
        str(result.bam_path),
        result.index,
        PipelineConfig(),
        tmp_path,
        tag,
        oracle_cache=None,
    )


@pytest.fixture(scope="module")
def measured(toy, tmp_path_factory):
    return _measure(toy, tmp_path_factory.mktemp("pv_run"), "toy")


@pytest.fixture(scope="module")
def measured_zero(toy_zero_gdna, tmp_path_factory):
    return _measure(toy_zero_gdna, tmp_path_factory.mktemp("pv0_run"), "toy_zero")


# ── GATE 1: the noop arm is byte-identical, AND it notices when it is not ────────────────────────


def _nudged_prior(measured, field, index, delta):
    """``assemble_priors`` re-run with ``+delta`` added to one element of one mass array."""
    cal = _rebuild_calibration(measured)
    arr = np.array(getattr(cal, field), dtype=np.float64, copy=True)
    arr[index] += delta
    return PV.PRIORS.assemble_priors(
        dataclasses.replace(cal, **{field: arr}), measured.region_arrays, measured.multi_loci
    )


def _moved(a, b) -> bool:
    return any(not np.array_equal(getattr(a, f), getattr(b, f)) for f in PV.PRIOR_FIELDS)


def test_the_noop_arm_is_byte_identical_and_the_lever_resolves_a_PICOFRAGMENT(measured):
    """TRAPS: byte-identity-gate. The whole instrument rests on ``dataclasses.replace`` being an inert way to swap
    the six mass arrays. If the replace dropped a field, O would differ from P for a reason that is
    not deconvolution error, and that bug would BE the headline number.

    The perturbation site has to be chosen by perturbation, and that is the gate's real content.
    Nudging ``argmax(count_gdna_region)`` reads "no effect", because on any genome the largest gDNA
    region is intergenic and ``_project_regions_to_loci`` drops every region overlapping no locus —
    correct behaviour, and it would retire the gate as broken (TRAPS: could-the-arm-have-fired). So
    this asserts both directions: in-locus moves, intergenic does not.

    The resolution: 1e-12 fragments at an in-locus region moves the prior by 1.0e-12, and one ULP
    (~3e-14 on a mass of a few hundred) does not, because the projection is a plain summation and a
    perturbation below the summand's own rounding is absorbed. 1e-12 fragments is twelve orders of
    magnitude below the unit the prior is denominated in, so the lever is not the limiting factor.
    """
    assert all(measured.noop_identical.values()), measured.noop_identical

    from rigel.calibration.priors import _RNA_SIGNATURE_BITS

    sig = np.asarray(measured.region_arrays.signature).astype(np.int64)
    in_locus = (sig & _RNA_SIGNATURE_BITS) != 0
    mass = np.asarray(measured.calibration.count_gdna_region, np.float64)
    inside = int(np.argmax(np.where(in_locus, mass, -1.0)))
    outside = int(np.argmax(np.where(~in_locus, mass, -1.0)))
    assert mass[inside] > 0.0 and mass[outside] > 0.0, (
        "the toy has no gDNA at an in-locus region or no gDNA at an intergenic one — one half of this "
        "gate would be vacuous"
    )

    base = measured.priors["P"]
    assert _moved(_nudged_prior(measured, "count_gdna_region", inside, 1e-12), base), (
        "1e-12 fragments at an in-locus region changed no prior — the lever cannot resolve an override"
    )
    # The intergenic direction is asserted on the COUNT field only: the locus projection dropping
    #   intergenic regions from the count is what this gate protects. ``gdna_eff_len`` is the
    #   contraction's business — it reads the result's reference density and the locus's own objects —
    #   and has its own gates in ``test_priors.py``; asserting all of ``PRIOR_FIELDS`` here conflated
    #   the two.
    nudged_out = _nudged_prior(measured, "count_gdna_region", outside, 1.0)
    count_moved = not np.array_equal(nudged_out.gdna_count, base.gdna_count)
    assert not count_moved, (
        "a whole fragment at an INTERGENIC region changed a locus COUNT prior — the locus "
        "projection is no longer dropping regions that overlap no locus"
    )


def test_the_prior_reads_exactly_the_two_gDNA_mass_fields_and_provably_NONE_of_the_RNA_ones(
    measured,
):
    """TRAPS: an-ablation-that-never-ran, applied to the override itself. The two gDNA arrays
    ``override_masses`` writes must reach the prior — an override landing on a field nothing reads is an
    override that never ran, and it would read as "calibration is already correct on that channel".

    And the four RNA arrays must provably NOT reach it. Calibration's RNA count does not enter the EM:
    ``pipeline.em_pseudocounts`` reads the gDNA count against the EM's own count of the locus's
    fragments, so the two stages never have to agree on which fragments are spliced. An RNA field that
    moved the prior would be that count leaking back in. The rule otherwise lives only in a docstring;
    this is what keeps it true. Each field is nudged at its largest IN-LOCUS site, so a field that does
    not move the prior is not read, never merely nudged where the projection drops it.
    """
    base = measured.priors["P"]
    reads = {}
    for field in PV.OVERRIDE_FIELDS:
        site = _biggest_in_locus_site(measured, field)
        reads[field] = _moved(_nudged_prior(measured, field, site, 1.0), base)
    assert {f for f, m in reads.items() if m} == {"count_gdna_region", "count_gdna_boundary"}, (
        f"the override fields the prior reads are {sorted(f for f, m in reads.items() if m)} — an RNA "
        "field reaching the prior is calibration's RNA count back in the EM, and a gDNA field missing "
        "is an override that never ran"
    )


def test_an_override_that_stops_writing_a_field_ABORTS_rather_than_scoring(measured, monkeypatch):
    """TRAPS: an-ablation-that-never-ran. If ``override_masses`` were changed to stop writing one of the six
    arrays, O would silently keep the SHIPPED value there and would be a hybrid of truth and estimate
    — which reads as "calibration is better than we thought" and is the most flattering possible bug.
    """
    cal = _rebuild_calibration(measured)
    full = measured.oracle.override_masses(measured.region_arrays)
    for dropped in PV.OVERRIDE_FIELDS:
        partial = {k: v for k, v in full.items() if k != dropped}
        monkeypatch.setattr(
            type(measured.oracle), "override_masses", lambda self, ra, _p=partial: dict(_p)
        )
        with pytest.raises(RuntimeError, match="OVERRIDE_FIELDS"):
            PV.oracle_priors(measured.oracle, cal, measured.region_arrays, measured.multi_loci)


# ── GATE 2: the lever COULD have fired ───────────────────────────────────────────────────────────


def test_the_oracle_lever_actually_MOVES_the_prior(measured):
    """TRAPS: could-the-arm-have-fired. "P equals O" would be the headline result of the whole campaign, so the
    one thing that must not produce it is a lever that did nothing. On a toy with real gDNA the two
    priors must differ at a substantial number of loci and by a substantial total.

    Stated as a floor on the COUNT of differing loci as well as on the total, because a single
    enormous locus difference and a broad small one are different findings and only one of them
    proves the lever reaches the whole array.
    """
    p, o = measured.priors["P"], measured.priors["O"]
    differ = ~np.isclose(p.gdna_count, o.gdna_count, rtol=1e-9, atol=1e-9)
    assert differ.sum() >= 2, f"the oracle lever moved {differ.sum()} loci — it barely fired"
    assert o.gdna_count.sum() > 0.0, "the oracle prior claims no gDNA on a contaminated toy"


# ── GATE 3: the two arms cannot be scored against a different locus partition ────────────────────


def test_scoring_against_a_DIFFERENT_locus_partition_raises(measured):
    """The locus partition is a function of the SCORING stage, not of the index — ``build_multi_loci``
    unions transcripts linked by scored fragments — so two runs of the pipeline can produce different
    numbers of loci. Index-aligning two such arrays would silently compare locus 7 of one run with
    locus 7 of another, which is not a small error.
    """
    p = measured.priors["P"].gdna_count
    o = measured.priors["O"].gdna_count
    PV.score_arm(p, o)  # the honest call
    with pytest.raises(ValueError, match="different locus partitions"):
        PV.score_arm(p, o[:-1])


# ── GATE 4: every sharded run writes its shards to a directory of its own ───────────────────────


def test_two_sharded_runs_never_share_a_shard_directory_and_every_shard_gets_the_set(
    tmp_path, monkeypatch
):
    """Two A/B runs sharing a work dir once wrote their shards to one fixed ``_shards``, so each could
    merge the other's rows. The shards are faked (each writes an empty row list), so this checks the
    plumbing only: two runs, two directories, and every shard command carries the run's ``--set``.

    Perturbation: a fixed shard directory makes the two runs' directories equal.
    """
    import subprocess
    import sys

    seen: list[list[str]] = []

    class _Shard:
        returncode = 0

        def __init__(self, cmd, **_kw):
            seen.append(cmd)
            Path(cmd[cmd.index("--json") + 1]).write_text("[]")

        def communicate(self):
            return "", None

    monkeypatch.setattr(subprocess, "Popen", _Shard)
    argv = [
        "prior_vs_oracle.py",
        "--suite",
        str(tmp_path / "suite"),
        "--work-dir",
        str(tmp_path),
        "--jobs",
        "2",
        "--set",
        "scan.total_threads=1",
        "--conditions",
        "a",
        "b",
    ]
    dirs = []
    for _run in range(2):
        seen.clear()
        monkeypatch.setattr(sys, "argv", argv)
        assert PV.main() == 0
        assert len(seen) == 2 and all("scan.total_threads=1" in cmd for cmd in seen)
        dirs.append({Path(cmd[cmd.index("--json") + 1]).parent for cmd in seen})
    assert len(dirs[0]) == len(dirs[1]) == 1 and dirs[0] != dirs[1]


# ── GATE 5: the ZERO-gDNA control ────────────────────────────────────────────────────────────────


def test_at_zero_gDNA_the_ORACLE_prior_is_identically_zero_and_the_shipped_one_is_scored_against_it(
    measured_zero,
):
    """THE OWNER-REQUIRED ZERO CONTROL. With no gDNA in the library the oracle's gDNA mass is
    exactly 0 at every object, so ``O.gdna_count`` must be exactly 0 at every locus — not small,
    not floored, zero. Anything the SHIPPED prior puts there is a false positive with nothing to
    cancel it, which is the only reading of that arm that is unambiguous.

    The perturbation is the other direction and it is what makes the assertion non-vacuous: hand
    the same assembler a single fabricated gDNA fragment and the prior must come off zero. A gate that
    only ever sees zeros cannot tell "correct" from "the array is not wired".
    """
    o = measured_zero.priors["O"]
    assert measured_zero.n_loci > 0, (
        "the zero-gDNA toy has no loci — the zero below is a sum over nothing"
    )
    assert float(np.asarray(o.gdna_count).sum()) == 0.0, (
        "the oracle prior claims gDNA in a library that has none"
    )

    cal = _rebuild_calibration(measured_zero)
    truth = measured_zero.oracle.override_masses(measured_zero.region_arrays)
    seeded = np.array(truth["count_gdna_region"], copy=True)
    # one gDNA fragment, seeded at the largest RNA region
    i = int(np.argmax(np.asarray(truth["count_rna_region"])))
    seeded[i] = 1.0
    with_one = PV.PRIORS.assemble_priors(
        dataclasses.replace(cal, **{**truth, "count_gdna_region": seeded}),
        measured_zero.region_arrays,
        measured_zero.multi_loci,
    )
    assert float(np.asarray(with_one.gdna_count).sum()) > 0.0, (
        "one fabricated gDNA fragment produced a prior of exactly zero — the zero above is the "
        "array being dead, not the library being clean"
    )


# ── GATE 6: a capture that never happened is an error, not a zero ────────────────────────────────


def test_a_run_that_never_reaches_assemble_priors_RAISES(measured, monkeypatch):
    """TRAPS: an-ablation-that-never-ran. ``quant_from_buffer`` returns early when there are no EM units, and a
    silently-absent capture would read as "a condition with no loci and therefore no error" — the most
    flattering possible failure of the harness.
    """
    monkeypatch.setattr(PV, "quant_from_buffer", lambda *a, **k: None)
    with pytest.raises(RuntimeError, match="never called"):
        PV.capture_priors(None, None, None, None, None, None, None, PipelineConfig())


# ── GATE 7: the aggregate re-derives its rates, never averages them ──────────────────────────────


def test_the_stratum_aggregate_is_a_RATIO_OF_SUMS_not_a_mean_of_ratios():
    """TRAPS: never-pool-the-strata's third way. A panel condition at 10 M fragments and one at 10 k are not two
    equally-informative opinions about a rate, and averaging their ``rel`` values gives the shallow
    one equal say. The aggregate must recompute from the summed totals.

    The two constructed rows differ by 1,000x in depth and have opposite-signed errors, so the mean
    of ratios and the ratio of sums are far apart and a gate that confused them could not pass by
    coincidence.
    """
    deep = PV.ArmScore(
        n_loci=1,
        n_claiming=1,
        total_arm=1_100_000.0,
        total_ref=1_000_000.0,
        net_err=100_000.0,
        abs_err=100_000.0,
        over_call=100_000.0,
        under_call=0.0,
    )
    shallow = PV.ArmScore(
        n_loci=1,
        n_claiming=1,
        total_arm=100.0,
        total_ref=1_000.0,
        net_err=-900.0,
        abs_err=900.0,
        over_call=0.0,
        under_call=900.0,
    )
    agg = PV._agg([deep, shallow])
    assert agg.rel == pytest.approx(100_900.0 / 1_001_000.0)
    assert agg.rel != pytest.approx((deep.rel + shallow.rel) / 2.0, rel=1e-3)
    assert agg.net_err == pytest.approx(99_100.0)
    assert agg.over_call == pytest.approx(100_000.0) and agg.under_call == pytest.approx(900.0)


# ── GATE 8: the directional split is reported and reconciles ─────────────────────────────────────


def test_over_and_under_call_are_reported_separately_and_reconcile(measured):
    """The library-level figure is ``|Σ(a − a*)|`` and the per-locus answer is ``Σ|a − a*|``; when a
    large under-call sits next to a large over-call the first flatters the second by whatever
    ``cancellation`` reports. Both halves must therefore exist and must reconcile exactly::

        over − under == net        over + under == abs
    """
    s = PV.score_arm(measured.priors["P"].gdna_count, measured.priors["O"].gdna_count)
    assert s.over_call - s.under_call == pytest.approx(s.net_err, rel=1e-9, abs=1e-6)
    assert s.over_call + s.under_call == pytest.approx(s.abs_err, rel=1e-9, abs=1e-6)
    assert s.over_call >= 0.0 and s.under_call >= 0.0


# ── helpers ──────────────────────────────────────────────────────────────────────────────────────


def _in_locus_regions(measured) -> np.ndarray:
    """``bool[N]`` — regions the locus projection keeps. An intergenic region carries no exon/intron bit
    and is DROPPED, so a perturbation there is inert by design (see the noop gate)."""
    from rigel.calibration.priors import _RNA_SIGNATURE_BITS

    sig = np.asarray(measured.region_arrays.signature).astype(np.int64)
    return (sig & _RNA_SIGNATURE_BITS) != 0


def _biggest_in_locus_site(measured, field):
    """The largest element of ``field`` that the locus projection actually reaches.

    Boundary-indexed arrays are selected through the SHIPPED ``_boundary_locus_shares`` rather than a local
    rule: a locus's boundaries are the boundaries that TOUCH its regions, which is exactly the decision that
    function exists to make. Restating it here would drift from the code under test.
    """
    from rigel.calibration.priors import _boundary_locus_shares

    arr = np.asarray(getattr(measured.calibration, field), np.float64)
    keep = _in_locus_regions(measured)
    if "boundary" in field:
        e_idx, _lid, _w = _boundary_locus_shares(
            measured.region_arrays, measured.multi_loci, len(measured.multi_loci)
        )
        keep = np.zeros(int(measured.calibration.n_boundaries), dtype=bool)
        keep[e_idx] = True
    if arr.shape[0] == keep.shape[0]:
        ranked = np.where(keep.reshape((-1,) + (1,) * (arr.ndim - 1)), arr, -np.inf)
    else:  # the sj axis — never projected, and never read by the prior
        ranked = arr
    site = np.unravel_index(int(np.argmax(ranked)), arr.shape)
    return site


def _rebuild_calibration(m):
    """The condition's own ``CalibrationResult``, reconstructed from the shipped masses the
    measurement kept rather than re-running ``calibrate``: the perturbations below only need an object
    whose six mass arrays are P's, and re-solving would take seconds per gate and could drift."""
    return m.calibration
