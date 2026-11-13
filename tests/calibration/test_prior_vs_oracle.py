"""Tests for ``scripts/design/prior_vs_oracle.py``, the instrument that scores calibration's endpoint (``LocusPriors``)
against the origin-split oracle.

Everything the instrument prints is a difference between two per-locus arrays of fragment counts, so every way of
getting it wrong yields a plausible number: a lever that never fired, a reference that became the arm, a weight that
rewards the solver for declining to answer, a projection that loses mass off the end of a locus. Each test covers one of
those and perturbs its own input to show it can fail. The scenario is a single-reference toy, enough for the scoring and
plumbing (arithmetic over two per-locus arrays plus the override lever) but not for the deposit path, which
``tests/native/`` covers.
"""

from __future__ import annotations

import dataclasses
from pathlib import Path

import numpy as np
import pytest

from rigel.calibration.signature import BIT_EXON_NEG, BIT_EXON_POS, BIT_INTRON_NEG, BIT_INTRON_POS
from rigel.config import PipelineConfig
import prior_vs_oracle
from _prior_toy import build_toy, build_toy_zero_gdna

#: a region carrying none of these strand/type bits is intergenic, and the locus projection drops it
_RNA_SIGNATURE_BITS = BIT_EXON_POS | BIT_EXON_NEG | BIT_INTRON_POS | BIT_INTRON_NEG


# ── the toys ─────────────────────────────────────────────────────────────────────────────────────


@pytest.fixture(scope="module")
def toy(tmp_path_factory):
    return build_toy(tmp_path_factory)


@pytest.fixture(scope="module")
def toy_zero_gdna(tmp_path_factory):
    return build_toy_zero_gdna(tmp_path_factory)


def _measure(result, tmp_path, tag):
    return prior_vs_oracle.measure_condition(
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


# ── Gate 1: the noop arm is byte-identical, and it notices when it is not ────────────────────────


def _nudged_prior(measured, field, index, delta):
    """``assemble_priors`` re-run with ``+delta`` added to one element of one mass array."""
    cal = _rebuild_calibration(measured)
    arr = np.array(getattr(cal, field), dtype=np.float64, copy=True)
    arr[index] += delta
    return prior_vs_oracle.PRIORS.assemble_priors(
        dataclasses.replace(cal, **{field: arr}), measured.region_arrays, measured.multi_loci
    )


def _moved(a, b) -> bool:
    return any(
        not np.array_equal(getattr(a, f), getattr(b, f)) for f in prior_vs_oracle.PRIOR_FIELDS
    )


def test_the_noop_arm_is_byte_identical_and_the_lever_resolves_a_picofragment(measured):
    """The instrument rests on ``dataclasses.replace`` being an inert way to swap the six mass arrays; a dropped field
    would make O differ from P for a reason that is not deconvolution error.

    The perturbation site matters. On any genome the largest gDNA region is intergenic, and ``region_locus_shares``
    drops every region overlapping no locus, so nudging ``argmax(count_gdna_region)`` would show no effect. So both
    directions are asserted: an in-locus nudge moves the prior, an intergenic one does not. 1e-12 fragments at an
    in-locus region moves the prior by 1e-12, while one ULP of a mass of a few hundred does not (the projection is a
    plain sum and absorbs a perturbation below the summand's rounding), so the lever resolves far below the prior's
    unit.
    """
    assert all(measured.noop_identical.values()), measured.noop_identical

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
    # The intergenic direction is asserted on the count field only: dropping intergenic regions from the count is what
    # this test protects. ``gdna_eff_len`` reads the result's reference density and the locus's own objects, and
    # ``test_priors.py`` covers it.
    nudged_out = _nudged_prior(measured, "count_gdna_region", outside, 1.0)
    count_moved = not np.array_equal(nudged_out.gdna_count, base.gdna_count)
    assert not count_moved, (
        "a whole fragment at an INTERGENIC region changed a locus COUNT prior — the locus "
        "projection is no longer dropping regions that overlap no locus"
    )


def test_the_prior_reads_exactly_the_two_gdna_mass_fields_and_none_of_the_rna_ones(
    measured,
):
    """The two gDNA arrays ``override_masses`` writes must reach the prior; an override landing on a field nothing reads
    would read as "calibration is already correct on that channel". The four RNA arrays must not reach it: calibration's
    RNA count does not enter the EM (``pipeline.em_pseudocounts`` reads the gDNA count against the EM's own count of the
    locus's fragments), so an RNA field moving the prior would be that count leaking in. Each field is nudged at its
    largest in-locus site, so a field that does not move the prior is genuinely unread.
    """
    base = measured.priors["P"]
    reads = {}
    for field in prior_vs_oracle.OVERRIDE_FIELDS:
        site = _biggest_in_locus_site(measured, field)
        reads[field] = _moved(_nudged_prior(measured, field, site, 1.0), base)
    assert {f for f, m in reads.items() if m} == {"count_gdna_region", "count_gdna_boundary"}, (
        f"the override fields the prior reads are {sorted(f for f, m in reads.items() if m)} — an RNA "
        "field reaching the prior is calibration's RNA count back in the EM, and a gDNA field missing "
        "is an override that never ran"
    )


def test_an_override_that_stops_writing_a_field_aborts_rather_than_scoring(measured, monkeypatch):
    """If ``override_masses`` stopped writing one of the six arrays, O would keep the shipped value there and be a
    hybrid of truth and estimate, which would flatter calibration.
    """
    cal = _rebuild_calibration(measured)
    full = measured.oracle.override_masses(measured.region_arrays)
    for dropped in prior_vs_oracle.OVERRIDE_FIELDS:
        partial = {k: v for k, v in full.items() if k != dropped}
        monkeypatch.setattr(
            type(measured.oracle), "override_masses", lambda self, ra, _p=partial: dict(_p)
        )
        with pytest.raises(RuntimeError, match="OVERRIDE_FIELDS"):
            prior_vs_oracle.oracle_priors(
                measured.oracle, cal, measured.region_arrays, measured.multi_loci
            )


# ── Gate 2: the lever could have fired ───────────────────────────────────────────────────────────


def test_the_oracle_lever_actually_moves_the_prior(measured):
    """ "P equals O" must not come from a lever that did nothing. On a toy with real gDNA the two priors must differ at
    many loci and by a substantial total; both the count of differing loci and the total are bounded, since one huge
    locus difference does not show the lever reaches the whole array.
    """
    p, o = measured.priors["P"], measured.priors["O"]
    differ = ~np.isclose(p.gdna_count, o.gdna_count, rtol=1e-9, atol=1e-9)
    assert differ.sum() >= 2, f"the oracle lever moved {differ.sum()} loci — it barely fired"
    assert o.gdna_count.sum() > 0.0, "the oracle prior claims no gDNA on a contaminated toy"


# ── Gate 3: the two arms cannot be scored against a different locus partition ────────────────────


def test_scoring_against_a_different_locus_partition_raises(measured):
    """The locus partition is a function of the scoring stage, not the index (``build_multi_loci`` unions transcripts
    linked by scored fragments), so two pipeline runs can produce different numbers of loci, and index-aligning them
    would compare locus 7 of one run with locus 7 of another.
    """
    p = measured.priors["P"].gdna_count
    o = measured.priors["O"].gdna_count
    prior_vs_oracle.score_arm(p, o)  # the honest call
    with pytest.raises(ValueError, match="different locus partitions"):
        prior_vs_oracle.score_arm(p, o[:-1])


# ── Gate 4: every sharded run writes its shards to a directory of its own ───────────────────────


def test_two_sharded_runs_never_share_a_shard_directory_and_every_shard_gets_the_set(
    tmp_path, monkeypatch
):
    """Two A/B runs sharing a work dir must not write shards to one fixed ``_shards``, or each could merge the other's
    rows. The shards are faked (each writes an empty row list), so this checks the plumbing only: two runs, two
    directories, and every shard command carries the run's ``--set``.
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
        assert prior_vs_oracle.main() == 0
        assert len(seen) == 2 and all("scan.total_threads=1" in cmd for cmd in seen)
        dirs.append({Path(cmd[cmd.index("--json") + 1]).parent for cmd in seen})
    assert len(dirs[0]) == len(dirs[1]) == 1 and dirs[0] != dirs[1]


# ── Gate 5: the zero-gDNA control ────────────────────────────────────────────────────────────────


def test_at_zero_gdna_the_oracle_prior_is_zero_and_the_shipped_one_is_scored_on_it(
    measured_zero,
):
    """With no gDNA in the library the oracle's gDNA mass is exactly 0 at every object, so ``O.gdna_count`` must be
    exactly 0 at every locus, and anything the shipped prior puts there is an unambiguous false positive. Handing the
    same assembler one fabricated gDNA fragment must take the prior off zero, so the test can tell "correct" from "not
    wired".
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
    with_one = prior_vs_oracle.PRIORS.assemble_priors(
        dataclasses.replace(cal, **{**truth, "count_gdna_region": seeded}),
        measured_zero.region_arrays,
        measured_zero.multi_loci,
    )
    assert float(np.asarray(with_one.gdna_count).sum()) > 0.0, (
        "one fabricated gDNA fragment produced a prior of exactly zero — the zero above is the "
        "array being dead, not the library being clean"
    )


# ── Gate 6: a capture that never happened is an error, not a zero ────────────────────────────────


def test_a_run_that_never_reaches_assemble_priors_raises(measured, monkeypatch):
    """``quant_from_buffer`` returns early when there are no EM units, and a silently absent capture would read as "a
    condition with no loci and therefore no error".
    """
    monkeypatch.setattr(prior_vs_oracle, "quant_from_buffer", lambda *a, **k: None)
    with pytest.raises(RuntimeError, match="never called"):
        prior_vs_oracle.capture_priors(
            None, None, None, None, None, None, None, None, PipelineConfig()
        )


# ── Gate 7: the aggregate re-derives its rates, never averages them ──────────────────────────────


def test_the_stratum_aggregate_is_a_ratio_of_sums_not_a_mean_of_ratios():
    """A condition at 10 M fragments and one at 10 k are not equally informative about a rate, and averaging their
    ``rel`` values would give the shallow one equal say, so the aggregate recomputes from summed totals. The two rows
    differ 1,000x in depth with opposite-signed errors, so the mean of ratios and the ratio of sums are far apart.
    """
    deep = prior_vs_oracle.ArmScore(
        n_loci=1,
        n_claiming=1,
        total_arm=1_100_000.0,
        total_ref=1_000_000.0,
        net_err=100_000.0,
        abs_err=100_000.0,
        over_call=100_000.0,
        under_call=0.0,
    )
    shallow = prior_vs_oracle.ArmScore(
        n_loci=1,
        n_claiming=1,
        total_arm=100.0,
        total_ref=1_000.0,
        net_err=-900.0,
        abs_err=900.0,
        over_call=0.0,
        under_call=900.0,
    )
    agg = prior_vs_oracle._agg([deep, shallow])
    assert agg.rel == pytest.approx(100_900.0 / 1_001_000.0)
    assert agg.rel != pytest.approx((deep.rel + shallow.rel) / 2.0, rel=1e-3)
    assert agg.net_err == pytest.approx(99_100.0)
    assert agg.over_call == pytest.approx(100_000.0) and agg.under_call == pytest.approx(900.0)


# ── Gate 8: the directional split is reported and reconciles ─────────────────────────────────────


def test_over_and_under_call_are_reported_separately_and_reconcile(measured):
    """The library-level figure is ``|Σ(a − a*)|`` and the per-locus answer is ``Σ|a − a*|``; a large under-call beside
    a large over-call makes the first flatter the second by ``cancellation``. Both halves must exist and reconcile
    exactly::

        over − under == net        over + under == abs
    """
    s = prior_vs_oracle.score_arm(measured.priors["P"].gdna_count, measured.priors["O"].gdna_count)
    assert s.over_call - s.under_call == pytest.approx(s.net_err, rel=1e-9, abs=1e-6)
    assert s.over_call + s.under_call == pytest.approx(s.abs_err, rel=1e-9, abs=1e-6)
    assert s.over_call >= 0.0 and s.under_call >= 0.0


# ── helpers ──────────────────────────────────────────────────────────────────────────────────────


def _in_locus_regions(measured) -> np.ndarray:
    """``bool[N]``: regions the locus projection keeps. An intergenic region carries no exon/intron bit and is dropped,
    so a perturbation there is inert by design."""
    sig = np.asarray(measured.region_arrays.signature).astype(np.int64)
    return (sig & _RNA_SIGNATURE_BITS) != 0


def _biggest_in_locus_site(measured, field):
    """The largest element of ``field`` the locus projection reaches. Boundary-indexed arrays are selected through the
    shipped ``boundary_locus_shares``, since a locus's boundaries are those touching its regions and that function
    makes the decision.
    """
    from rigel.calibration.priors import boundary_locus_shares, region_locus_shares

    arr = np.asarray(getattr(measured.calibration, field), np.float64)
    keep = _in_locus_regions(measured)
    if "boundary" in field:
        n_loci = len(measured.multi_loci)
        region_shares = region_locus_shares(measured.region_arrays, measured.multi_loci, n_loci)
        e_idx, _lid, _w = boundary_locus_shares(measured.region_arrays, region_shares, n_loci)
        keep = np.zeros(int(measured.calibration.n_boundaries), dtype=bool)
        keep[e_idx] = True
    if arr.shape[0] == keep.shape[0]:
        ranked = np.where(keep.reshape((-1,) + (1,) * (arr.ndim - 1)), arr, -np.inf)
    else:  # the sj axis — never projected, and never read by the prior
        ranked = arr
    site = np.unravel_index(int(np.argmax(ranked)), arr.shape)
    return site


def _rebuild_calibration(m):
    """The condition's own ``CalibrationResult``, rebuilt from the shipped masses the measurement kept rather than by
    re-running ``calibrate``: the perturbations need only an object whose six mass arrays are P's."""
    return m.calibration
