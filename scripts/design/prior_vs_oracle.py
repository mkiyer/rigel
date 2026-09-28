#!/usr/bin/env python
"""Is ``LocusPriors`` -- what the EM takes from calibration -- right?

Calibration does not ship a number; it ships two float64 arrays indexed by ``multi_locus_id``:
``gdna_count``, its conserved count of the locus's gDNA fragments, and ``gdna_eff_len``, the gDNA
component's effective length. The EM reads ``gdna_count`` against its OWN count of the locus's fragments
(``pipeline.em_pseudocounts``: the odds are ``G : N - G``, with ``N`` the locus's units plus its
deterministic fragments), so the composition the EM is handed is ``gdna_count / N`` over a denominator
calibration does not supply: its error IS the count's error, divided by a number every arm shares, and
there is no separate composition or scale to score. Calibration's RNA count does not reach the EM at all,
so there is no RNA arm. This instrument scores the gDNA count, per condition and per stratum, against the
origin-split oracle. Three arms: ``P`` is the shipped prior, ``O`` is the same assembler fed the true
per-object masses (``OracleTruth.override_masses``, the one lever that exists -- O is never an
estimator), and ``S`` is O with the crossing term converted by gDNA's own true share instead of the pooled
one. ``P - O`` is calibration's error and ``O - S`` the pooled share's part of the assembler's
(``ISSUES: the-pooled-q-in-the-gdna-count``). The count's error is reported in fragments (``sum |dA|``,
additive). ``gdna_eff_len`` reads none of the six override fields, so O carries P's length
and it is not scored here; its truth is ``ruler_vs_truth.py``'s. Every arm runs in the drained frame; the
drain's spliced-gDNA leak is reported beside the numbers as ``gdna_spliced_leak`` and the lift's attribution
error as ``n_ambiguous``. No per-locus EM runs -- the
pipeline is stopped after its scoring stage.

The panel's default paths, the six override fields and the stratum readers this and the other oracle
instruments share live in ``_shared`` (``DEFAULT_SUITE``, ``DEFAULT_INDEX``, ``OVERRIDE_FIELDS``, ``stratum``,
``is_zero_gdna``), so no instrument loads another to read them. Gates:
``tests/calibration/test_prior_vs_oracle.py``.

Usage::

    python scripts/design/prior_vs_oracle.py --suite DIR --oracle-cache DIR/oracle_cache --jobs 6
    python scripts/design/prior_vs_oracle.py --conditions gdna_g50_ss_0.50_nrna_mid_capture_on
    python scripts/design/prior_vs_oracle.py --index INDEX --work-dir SCRATCH
    python scripts/design/prior_vs_oracle.py --set scan.total_threads=1   # a pinned, reproducible run
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import os
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402

_REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO / "tests" / "calibration"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from _shared import DEFAULT_INDEX, DEFAULT_SUITE, OVERRIDE_FIELDS, is_zero_gdna, set_field, strata, stratum  # noqa: E402,F401

from _oracle import ORIGINS, OracleTruth, _split_bam  # noqa: E402

import rigel.calibration.priors as PRIORS  # noqa: E402
from rigel.calibration.calibrate import calibrate  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402
from rigel.pipeline import (  # noqa: E402
    _drain_side_buffer,
    _native_detect_sj_tag,
    quant_from_buffer,
    scan_and_buffer,
)
from rigel.scan_cache import ScanCacheKeyError, read_scan_cache, write_scan_cache  # noqa: E402

#: The two ``LocusPriors`` fields, in the order every table prints them.
PRIOR_FIELDS = ("gdna_count", "gdna_eff_len")



# ── the scoring ──────────────────────────────────────────────────────────────────────────────────


@dataclass(frozen=True, slots=True)
class ArmScore:
    """One arm's per-locus gDNA count scored against one reference's, over one locus selection.

    Every field is a fragment total or a locus count, and the two rates ``rel`` and ``cancellation`` are
    derived from them — so the totals ADD across a partition of the loci and the two rates do not.
    """

    n_loci: int  #: loci in the selection. never the panel's locus count
    n_claiming: int  #: loci where EITHER arm or reference puts a nonzero gDNA count
    total_arm: float  #: Σ a, the arm's own fragment total
    total_ref: float  #: Σ a*, the reference's
    net_err: float  #: Σ (a − a*). What a library-level figure would see
    abs_err: float  #: Σ |a − a*|. What the per-locus answer is
    over_call: float  #: Σ (a − a*)+
    under_call: float  #: Σ (a* − a)+

    @property
    def rel(self) -> float:
        """``Σ|a − a*| / Σ a*`` — the absolute error as a fraction of the reference total."""
        return self.abs_err / self.total_ref if self.total_ref > 0 else float("nan")

    @property
    def cancellation(self) -> float:
        """``Σ|err| / |net|``. Large means a big under-call is sitting next to a big over-call and
        any library-level summary of this arm is flattering it."""
        return self.abs_err / abs(self.net_err) if self.net_err != 0.0 else float("inf")


def score_arm(arm: np.ndarray, ref: np.ndarray, select: np.ndarray | None = None) -> ArmScore:
    """Score one arm's per-locus gDNA count (a fragment count) against a reference's."""
    a = np.asarray(arm, np.float64)
    r = np.asarray(ref, np.float64)
    if a.shape != r.shape:
        raise ValueError(
            f"arm has {a.shape[0]} loci, reference has {r.shape[0]} — these are different locus "
            "partitions. A prior scored against another run's loci is not a small error."
        )
    if select is not None:
        a, r = a[select], r[select]
    d = a - r
    return ArmScore(
        n_loci=int(a.shape[0]),
        n_claiming=int(np.count_nonzero((a > 0.0) | (r > 0.0))),
        total_arm=float(a.sum()),
        total_ref=float(r.sum()),
        net_err=float(d.sum()),
        abs_err=float(np.abs(d).sum()),
        over_call=float(np.maximum(d, 0.0).sum()),
        under_call=float(np.maximum(-d, 0.0).sum()),
    )


# ── the arms ─────────────────────────────────────────────────────────────────────────────────────


def capture_priors(buffer, index, strand_models, fl, region_arrays, stats, calibration,
                   pipeline_config):
    """Run the production quant path far enough to get ``(multi_loci, LocusPriors)``, then STOP.

    The loci are the production loci, not a re-derivation. ``build_multi_loci`` builds
    connected components of transcripts linked by SCORED fragments, so the locus partition is a
    function of the scoring stage and cannot be reconstructed from the index alone. Wrapping
    :func:`~rigel.calibration.priors.assemble_priors` — which ``quant_from_buffer`` imports
    function-locally, so patching the module attribute is picked up at call time — takes both objects
    from the call production itself makes (TRAPS: a-test-that-redefines).

    The sentinel exception is what makes this affordable: the per-locus EM is the single most
    expensive stage and this instrument does not read its output. An experiment that injects the
    oracle prior and re-quantifies needs the EM and must not reuse this path.
    """

    class _StopAfterPriors(Exception):
        pass

    captured: dict = {}
    original = PRIORS.assemble_priors

    def _wrapper(cal, ra, multi_loci):
        captured["multi_loci"] = multi_loci
        captured["priors"] = original(cal, ra, multi_loci)
        raise _StopAfterPriors

    PRIORS.assemble_priors = _wrapper
    try:
        quant_from_buffer(
            buffer, index, strand_models, fl, region_arrays, stats, calibration,
            em_config=pipeline_config.em, scoring=pipeline_config.scoring,
        )
    except _StopAfterPriors:
        pass
    finally:
        PRIORS.assemble_priors = original
    if "priors" not in captured:
        # TRAPS: an-ablation-that-never-ran. ``quant_from_buffer`` returns early on an empty unit set, and a
        # silently-absent capture would read here as a condition with no loci rather than as a
        # harness that never fired.
        raise RuntimeError(
            "assemble_priors was never called — quant_from_buffer returned before the prior stage "
            "(no EM units, or no multi-loci). This is not a condition with zero error."
        )
    return captured["multi_loci"], captured["priors"]


def oracle_priors(oracle: OracleTruth, calibration, region_arrays, multi_loci):
    """O — the prior a perfect deconvolution would produce, and the ``noop`` gate beside it.

    The ONE lever: ``override_masses`` swaps the six mass arrays for the origin-split truth and
    changes nothing else. The ``noop`` arm re-injects the SHIPPED masses through the identical
    ``dataclasses.replace`` and must come back byte-identical (TRAPS: byte-identity-gate) — an override that silently
    landed on the wrong field, or a ``replace`` that dropped a field, would otherwise look like a
    result. Returns ``(O, noop)``.
    """
    override = oracle.override_masses(region_arrays)
    missing = set(OVERRIDE_FIELDS) - set(override)
    extra = set(override) - set(OVERRIDE_FIELDS)
    if missing or extra:
        raise RuntimeError(
            f"override_masses no longer writes OVERRIDE_FIELDS: missing={sorted(missing)} "
            f"extra={sorted(extra)}. The noop gate would be testing a different set than the arm."
        )
    o = PRIORS.assemble_priors(
        dataclasses.replace(calibration, **override), region_arrays, multi_loci
    )
    noop = PRIORS.assemble_priors(
        dataclasses.replace(calibration, **{f: getattr(calibration, f) for f in OVERRIDE_FIELDS}),
        region_arrays,
        multi_loci,
    )
    return o, noop


def share_priors(oracle: OracleTruth, calibration, region_arrays, multi_loci):
    """S — the O arm, with the gDNA count's crossing term converted by gDNA's OWN true per-boundary share.

    WHY THIS ARM EXISTS. ``assemble_priors`` converts the gDNA count's crossing term by ONE pooled
    share, ``q = mass / count`` off the mixture — RNA's where RNA dominates a boundary (the length moved to
    gDNA's own conserved share and the count did not, ``ISSUES: the-pooled-q-in-the-gdna-length``). Where
    the two components' shares differ, the gDNA count is off by ``q / share_g`` at that boundary, and
    whatever the pooled share moves off gDNA it moves onto RNA, so no gate on the locus total can see it.
    Only a per-component comparison can.

    ``O − S`` is therefore the pooled share's own contribution, isolated.

    It is the shipped function, run once with one input varied — ``boundary_mass_per_crossing`` set to
    gDNA's true share — never a re-implementation. gDNA's share is the only one the arm needs:
    calibration's RNA count does not reach the EM, so the pooled share's RNA half has no consumer.

    The shares are MEASURED off the origin split (``OracleTruth.component_shares``), never derived
    from a pmf — see that method for why an analytic share would make this a model arm.
    """
    shares = oracle.component_shares()
    truth_cal = dataclasses.replace(calibration, **oracle.override_masses(region_arrays))
    return PRIORS.assemble_priors(
        dataclasses.replace(truth_cal, boundary_mass_per_crossing=shares["gdna"]),
        region_arrays, multi_loci,
    ), shares


def eff_len_inflation(calibration, region_arrays, multi_loci) -> dict:
    """Is ``gdna_eff_len``'s span the locus's genomic extent, each start counted once?

    ``assemble_priors`` clamps ``gdna_eff_len`` to ``span = Σ share·(S_region + M_boundary)``, ``M`` gDNA's
    conserved share at each boundary (``gdna_boundary_conserved_len``): a crossing start's unit split over the
    boundaries its fragment crosses, so the span over a locus of pieces shorter than a fragment is still its
    extent plus the fragments straddling its two ends, never a fragment length per interior boundary. The EM
    divides the gDNA component's abundance by this array, so an inflation here is a direct scale error on one
    of the two numbers calibration ships.

    Reports the ratio to the locus's GENOMIC span, mass-weighted by the gDNA count, so the number is
    what the consumer feels rather than what an unweighted locus average would say
    (``TRAPS: weight-it-like-the-consumer``).
    """
    # Regions and boundaries are projected on their OWN axes, exactly as `assemble_priors` does — no boundary is
    # folded onto a flank region. Re-deriving the fold here would measure a span the assembler no longer
    # builds (`TRAPS: a-test-that-redefines`).
    n_loci = len(multi_loci)
    region_support = np.maximum(np.asarray(calibration.gdna_region_eff_len, np.float64), 0.0)
    boundary_support = np.maximum(np.asarray(calibration.gdna_boundary_conserved_len, np.float64), 0.0)
    proj = PRIORS._project_regions_to_loci(
        region_arrays, multi_loci, n_loci,
        {
            "region_only": region_support,
            "genomic": np.asarray(region_arrays.region_size_bp, np.float64),
        },
    )
    e_idx, e_lid, e_w = PRIORS._boundary_locus_shares(region_arrays, multi_loci, n_loci)
    proj["support"] = proj["region_only"] + PRIORS._sum_by_locus(
        e_idx, e_lid, e_w, boundary_support, n_loci
    )
    live = proj["genomic"] > 0
    ratio = np.divide(proj["support"], proj["genomic"], out=np.ones_like(proj["support"]), where=live)
    region_ratio = np.divide(proj["region_only"], proj["genomic"],
                           out=np.ones_like(proj["support"]), where=live)
    return {
        "median_support_over_genomic": float(np.median(ratio[live])) if live.any() else float("nan"),
        "median_region_over_genomic": float(np.median(region_ratio[live])) if live.any() else float("nan"),
        "total_support": float(proj["support"].sum()),
        "total_genomic": float(proj["genomic"].sum()),
    }


# ── one condition ────────────────────────────────────────────────────────────────────────────────


@dataclass
class ConditionResult:
    """Everything one condition produced. The arrays are kept so a gate can re-score them."""

    condition: str
    n_loci: int
    priors: dict  #: arm name -> LocusPriors
    noop_identical: dict  #: field -> bool
    #: is gdna_eff_len clamped by an incidence sum? See :func:`eff_len_inflation`
    eff_len: dict
    drain: dict  #: the measured drain caveat
    library: dict  #: condition-level totals, for the header row
    seconds: float
    #: The three inputs a gate needs to RE-DERIVE any arm above rather than trust the arrays it was
    #: handed. Kept deliberately: a falsification test that can only read the outputs can check they
    #: are self-consistent and never that they are the outputs of the shipped code.
    oracle: OracleTruth = None
    region_arrays: object = None
    multi_loci: list = None
    calibration: object = None


def _oracle_parts(bam, index, scan, work_dir, tag, cache_root):
    """The three origin partitions, from the cache when it is valid — otherwise split and scan.

    Keyed by the SHIPPED ``read_scan_cache``, never a home-made key: ``reach`` is covered by no other
    hash, so a rebuilt index would verify clean against one. Same argument as
    ``_oracle_arms.load_or_build_oracle``, which this deliberately mirrors.
    """
    if cache_root is not None:
        dirs = {k: Path(cache_root) / tag / k for k in ORIGINS}
        try:
            return {k: read_scan_cache(dirs[k], index, scan).payload for k in ORIGINS}
        except (FileNotFoundError, KeyError, ScanCacheKeyError):
            pass
    paths, _counts = _split_bam(bam, Path(work_dir), tag)
    parts = {}
    for origin in ORIGINS:
        _s, strand_model, _b, payload = scan_and_buffer(paths[origin], index, scan)
        parts[origin] = payload
        if cache_root is not None:
            d = Path(cache_root) / tag / origin
            write_scan_cache(d, payload=payload, strand_model=strand_model, index=index,
                             bam=paths[origin], scan_config=scan)
    return parts


def _calibrate_and_prior(payload, strand_model, buffer, stats, index, ra, pipeline_config):
    """calibrate → score → ``(calibration, multi_loci, LocusPriors)`` on ONE payload."""
    from rigel.calibration.splice_graph import build_boundary_flags_array, build_sj_geometry_arrays
    from rigel.pipeline import library_fl_models

    fl = library_fl_models(payload, index)
    cal = calibrate(
        payload=payload,
        region_arrays=ra,
        strand_model=strand_model,
        gdna_fl_pmf=fl.gdna_pmf,
        rna_fl_pmf=fl.rna_pmf,
        config=pipeline_config.calibration,
        sj=build_sj_geometry_arrays(index),
        boundary_flags=build_boundary_flags_array(index),
    )
    multi_loci, priors = capture_priors(
        buffer, index, strand_model, fl, ra, stats, cal, pipeline_config
    )
    return cal, multi_loci, priors


def measure_condition(bam, index, pipeline_config, work_dir, tag, *, oracle_cache=None) -> ConditionResult:
    """Scan once, drain to the production frame, build P / O / S."""
    start = time.perf_counter()
    scan = dataclasses.replace(pipeline_config.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    ra = RegionArrays.from_index(index)

    stats, strand_model, buffer, payload = scan_and_buffer(bam, index, scan)
    n_held = int(payload.deferred.n_fragments)
    # the drained frame: drained at the production seed; the oracle partitions are drained by
    # replaying the whole's choices, and sum-to-full then validates the lift.
    lift: dict = {}
    payload = _drain_side_buffer(
        payload, index, strand_model, seed=pipeline_config.second_pass_seed, _lift=lift
    )
    parts = _oracle_parts(bam, index, scan, work_dir, tag, oracle_cache)
    oracle = OracleTruth.from_cached_parts(payload, parts, lift)  # raises if the lift breaks a bank

    cal, multi_loci, p_arm = _calibrate_and_prior(
        payload, strand_model, buffer, stats, index, ra, pipeline_config
    )
    o_arm, noop = oracle_priors(oracle, cal, ra, multi_loci)
    s_arm, _shares = share_priors(oracle, cal, ra, multi_loci)
    eff_len = eff_len_inflation(cal, ra, multi_loci)

    noop_identical = {
        f: bool(np.array_equal(getattr(noop, f), getattr(p_arm, f))) for f in PRIOR_FIELDS
    }

    # the drained-frame report: the leak is production's own behaviour, recorded beside the
    # numbers it rides with (`ISSUES: drain-contaminates-certified-rna`); ``n_ambiguous`` bounds the
    # lift's origin attribution. ``gdna_spliced_leak`` is ``None`` only when nothing was held.
    drain: dict = {
        "n_held": n_held,
        "n_ambiguous": int(oracle.n_ambiguous),
        "gdna_spliced_leak": (
            None if oracle.gdna_spliced_leak is None
            else int(sum(oracle.gdna_spliced_leak.values()))
        ),
    }

    g = float(np.asarray(oracle.parts["gdna"].region_start_count, np.int64).sum())
    r = float(
        np.asarray(oracle.parts["mrna"].region_start_count, np.int64).sum()
        + np.asarray(oracle.parts["nrna"].region_start_count, np.int64).sum()
    )
    return ConditionResult(
        condition=tag,
        n_loci=len(multi_loci),
        priors={"P": p_arm, "O": o_arm, "S": s_arm},
        noop_identical=noop_identical,
        eff_len=eff_len,
        drain=drain,
        library={
            "true_gdna_fragments": g,
            "true_rna_fragments": r,
            "true_f_gdna": g / (g + r) if (g + r) > 0 else 0.0,
        },
        seconds=time.perf_counter() - start,
        oracle=oracle,
        region_arrays=ra,
        multi_loci=multi_loci,
        calibration=cal,
    )


# ── reporting ────────────────────────────────────────────────────────────────────────────────────


def _agg(scores):
    """Sum a list of ``ArmScore`` into one. The rates are RE-DERIVED from the summed totals, never
    averaged: a mean of ratios over conditions of wildly different depth is a number with no
    consumer (TRAPS: never-pool-the-strata's third way)."""
    scores = [s for s in scores if s is not None]
    if not scores:
        return None
    return ArmScore(
        n_loci=sum(s.n_loci for s in scores),
        n_claiming=sum(s.n_claiming for s in scores),
        total_arm=sum(s.total_arm for s in scores),
        total_ref=sum(s.total_ref for s in scores),
        net_err=sum(s.net_err for s in scores),
        abs_err=sum(s.abs_err for s in scores),
        over_call=sum(s.over_call for s in scores),
        under_call=sum(s.under_call for s in scores),
    )


def _selections(conds) -> tuple:
    """Every selection every table prints: one per stratum, never pooled, g00 excluded (it is read per
    condition). One list, so a stratum appears on every table at once; a stratum apart from both strand
    halves gets its own row (`_shared.strata`)."""
    return tuple(
        (" x ".join(st), (lambda c, st=st: stratum(c) == st and not is_zero_gdna(c)))
        for st in strata(conds)
    )


def _rel(x: float) -> str:
    """``rel`` in 8 columns. Switches to scientific below 1e-3 rather than rounding to ``0.000``: a
    small residual is not zero, and a fixed 3-decimal format would print it as one."""
    if not np.isfinite(x):
        return f"{'nan':>8}"
    return f"{x:>8.1e}" if 0.0 < abs(x) < 1e-3 else f"{x:>8.4f}"


def report(rows: list[dict], settings=()) -> None:
    """The whole report, from the per-condition JSON — the only report path there is.

    It reads JSON in the serial case too, and that is the point: a `--jobs 1` run and a `--jobs 6`
    run print through the same table builder from the same source, so the shard merge is not a
    special case and nothing is produced by two code paths (`TRAPS: a-test-that-redefines`).
    """
    rows = sorted(rows, key=lambda r: r["condition"])
    selections = _selections(r["condition"] for r in rows)
    print()
    print("=" * 104)
    print("  ⭐⭐⭐ CALIBRATION'S ENDPOINT vs THE ORACLE — LocusPriors, in FRAGMENTS")
    config = " ".join(f"--set {spec}" for spec in settings) or "the shipped defaults"
    print(f"  {len(rows)} conditions   config: {config}   the drained frame (the drain's leak is "
          "reported below)")
    print("=" * 104)

    # ── the gate first. A table read before its gate is a table nobody checked. ──
    bad = [(r["condition"], f) for r in rows for f, ok in r["noop_identical"].items() if not ok]
    if bad:
        print(f"\n  ⛔⛔ NOOP GATE FAILED on {len(bad)} (condition, field) pairs — the override "
              "plumbing is NOT inert and every number below is unattributable:")
        for c, f in bad[:10]:
            print(f"       {c}  {f}")
        raise SystemExit(2)
    print(f"\n  ✅ noop gate: re-injecting the shipped masses reproduces P byte-identically on all "
          f"{len(rows)} x {len(PRIOR_FIELDS)} arrays (TRAPS: byte-identity-gate)")

    def arm_table(title: str, key: str, note: str = "") -> None:
        print()
        print(f"  {title}")
        if note:
            print(f"  {note}")
        print(f"    {'stratum':<26} {'ref total':>14} {'arm total':>14} {'Σ|Δ|':>14} "
              f"{'rel':>8} {'net':>14} {'canc':>7}")
        print("    " + "-" * 102)
        for label, sel in selections:
            s = _agg([ArmScore(**r[key]) for r in rows if sel(r["condition"])])
            if s is None:
                print(f"    {label:<26} {'(empty)':>14}")
                continue
            print(f"    {label:<26} {s.total_ref:>14,.0f} {s.total_arm:>14,.0f} "
                  f"{s.abs_err:>14,.0f} {_rel(s.rel)} {s.net_err:>+14,.0f} "
                  f"{s.cancellation:>7.1f}")

    arm_table("① P vs O — CALIBRATION'S OWN ERROR (a perfect deconvolution, same assembler) · gdna_count",
              "P_vs_O")
    arm_table("② O vs S — THE POOLED SHARE'S OWN CONTRIBUTION, ISOLATED · gdna_count",
              "O_vs_S",
              "⛔ What the pooled share moves off gDNA it moves onto RNA, so no gate on the locus total "
              "can see it (ISSUES: the-pooled-q-in-the-gdna-count).")

    # ── ③ is gdna_eff_len clamped by an incidence sum? ──
    print()
    print("  ③ ⭐ gdna_eff_len's CLAMP — is `span` a genomic extent or an INCIDENCE sum?")
    print(f"    {'stratum':<26} {'support/genomic':>16} {'regions only':>12} {'Σ support':>16} "
          f"{'Σ genomic':>16}")
    print("    " + "-" * 92)
    for label, sel in selections:
        sub_rows = [r["eff_len"] for r in rows if sel(r["condition"]) and r.get("eff_len")]
        if not sub_rows:
            continue
        print(f"    {label:<26} "
              f"{float(np.median([x['median_support_over_genomic'] for x in sub_rows])):>16.2f} "
              f"{float(np.median([x['median_region_over_genomic'] for x in sub_rows])):>12.2f} "
              f"{sum(x['total_support'] for x in sub_rows):>16,.0f} "
              f"{sum(x['total_genomic'] for x in sub_rows):>16,.0f}")
    print("    ⚠ `support/genomic` well above 1 means every interior boundary is adding ~mu_g − 1 to the "
          "locus's clamp.")
    print("    The EM divides the gDNA component's abundance by gdna_eff_len, so this is a direct "
          "scale error on a shipped number.")

    # ── ④ the drained-frame report ──
    print()
    print("  ④ THE DRAINED FRAME — the lift's attribution bound, and production's certified-RNA leak")
    d = [r["drain"] for r in rows]
    amb, held = sum(x["n_ambiguous"] for x in d), sum(x["n_held"] for x in d)
    leak = sum(x["gdna_spliced_leak"] or 0 for x in d)
    n_leak = sum(1 for x in d if (x["gdna_spliced_leak"] or 0) > 0)
    print(f"     {held:,} held fragments drained at the production seed; lift ambiguity "
          f"{amb:,} ({amb / max(held, 1):.2%}) bounds the per-origin truth attribution.")
    print(f"     ⚠ the drain deposited {leak:,} spliced records into the gdna partition on "
          f"{n_leak}/{len(d)} conditions — production behaviour, RECORDED "
          "(ISSUES: drain-contaminates-certified-rna), never an oracle refusal.")

    # ── ⑤ the per-condition ladder, because a stratum total hides the shape ──
    print()
    print("  ⑤ PER CONDITION — gdna_count, P and O as fragment totals")
    print(f"    {'condition':<44} {'true f_g':>9} {'O_g':>13} {'P_g':>13} {'P/O':>8} "
          f"{'Σ|Δ|':>13} {'rel':>8}")
    print("    " + "-" * 116)
    for r in rows:
        g = ArmScore(**r["P_vs_O"])
        ratio = g.total_arm / g.total_ref if g.total_ref > 0 else float("nan")
        print(f"    {r['condition']:<44} {r['library']['true_f_gdna']:>9.4f} {g.total_ref:>13,.0f} "
              f"{g.total_arm:>13,.0f} {ratio:>8.3f} {g.abs_err:>13,.0f} {_rel(g.rel)}")

    print()
    print(f"  total wall clock {sum(r['seconds'] for r in rows):,.0f} s over {len(rows)} conditions")


def to_json(results: list[ConditionResult]) -> list[dict]:
    out = []
    for r in results:
        row: dict = {
            "condition": r.condition,
            "stratum": list(stratum(r.condition)),
            "n_loci": r.n_loci,
            "noop_identical": r.noop_identical,
            "drain": r.drain,
            "eff_len": r.eff_len,
            "library": r.library,
            "seconds": r.seconds,
        }
        p, o, sa = r.priors["P"], r.priors["O"], r.priors["S"]
        for ref_name, ref, arm in (
            ("P_vs_O", o.gdna_count, p.gdna_count),
            ("O_vs_S", sa.gdna_count, o.gdna_count),
        ):
            row[ref_name] = dataclasses.asdict(score_arm(arm, ref))
        out.append(row)
    return out


# ── main ─────────────────────────────────────────────────────────────────────────────────────────


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--suite", type=Path, default=DEFAULT_SUITE)
    ap.add_argument("--index", type=Path, default=DEFAULT_INDEX)
    ap.add_argument("--conditions", nargs="*", default=None)
    ap.add_argument("--oracle-cache", type=Path, default=None,
                    help="defaults to <suite>/oracle_cache when that directory exists")
    ap.add_argument("--work-dir", type=Path,
                    default=Path(os.environ.get("RIGEL_SCRATCH", "/tmp")) / "rigel_prior_oracle")
    ap.add_argument("--json", type=Path, default=None)
    ap.add_argument("--jobs", type=int, default=1,
                    help="run this many conditions CONCURRENTLY by re-invoking on shards. The "
                         "conditions are independent, so this changes no number.")
    ap.add_argument("--set", dest="settings", action="append", default=[], metavar="SECTION.FIELD=VALUE",
                    help="override one PipelineConfig field, typed from the field; repeatable "
                         "(e.g. --set scan.total_threads=1 for a reproducible run)")
    args = ap.parse_args()

    names = args.conditions or sorted(
        p.name for p in args.suite.iterdir() if (p / "sim_oracle.bam").is_file()
    )
    if not names:
        raise SystemExit(f"no conditions with a sim_oracle.bam under {args.suite}")
    cache = args.oracle_cache
    if cache is None and (args.suite / "oracle_cache").is_dir():
        cache = args.suite / "oracle_cache"

    if args.jobs > 1 and len(names) > 1:
        # Shards, not threads: conditions share nothing but a read-only index and cache, and
        # re-invoking the single-process path keeps the measured code byte-for-byte the serial one.
        # OMP_NUM_THREADS=1 is forced at import so the workers do not fight. Each run shards into its
        # own directory, so two concurrent runs sharing a work dir never read each other's shards.
        shards = [s for s in (names[i:: args.jobs] for i in range(args.jobs)) if s]
        args.work_dir.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=args.work_dir, prefix="shards_") as td:
            procs, outs = [], []
            for i, sh in enumerate(shards):
                o = Path(td) / f"{i}.json"
                outs.append(o)
                cmd = [sys.executable, str(Path(__file__).resolve()),
                       "--suite", str(args.suite), "--index", str(args.index),
                       "--work-dir", str(Path(td) / f"shard{i}"), "--json", str(o),
                       "--conditions", *sh]
                if cache is not None:
                    cmd += ["--oracle-cache", str(cache)]
                for spec in args.settings:
                    cmd += ["--set", spec]
                procs.append(subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                              stderr=subprocess.STDOUT, text=True))
            rc = 0
            for i, pr in enumerate(procs):
                out, _ = pr.communicate()
                if pr.returncode != 0:
                    rc = pr.returncode
                    print(f"  ⛔ shard {i} FAILED (rc={pr.returncode}):\n{out}", flush=True)
                else:
                    print(f"  shard {i}: {len(shards[i])} conditions ok", flush=True)
            if rc:
                # TRAPS: an-ablation-that-never-ran's shape — a short output file reads as a complete panel.
                raise SystemExit("a shard failed; refusing to report a partial panel")
            merged = [row for o in outs for row in json.loads(o.read_text())]
        if args.json is not None:
            args.json.write_text(json.dumps(merged, indent=1))
        report(merged, args.settings)
        return 0

    index = TranscriptIndex.load(str(args.index))
    pipeline_config = PipelineConfig()
    for spec in args.settings:
        pipeline_config = set_field(pipeline_config, spec)
    args.work_dir.mkdir(parents=True, exist_ok=True)
    results = []
    for name in names:
        bam = str(args.suite / name / "sim_oracle.bam")
        print(f"  … {name}", flush=True)
        results.append(measure_condition(
            bam, index, pipeline_config, args.work_dir, name,
            oracle_cache=cache,
        ))
    payload = to_json(results)
    if args.json is not None:
        args.json.write_text(json.dumps(payload, indent=1))
    report(payload, args.settings)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
