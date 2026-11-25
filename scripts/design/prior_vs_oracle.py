#!/usr/bin/env python
"""Is ``LocusPriors`` -- what the EM takes from calibration -- right?

Calibration hands the EM, per ``multi_locus_id``, ``gdna_count`` (the locus's conserved gDNA fragment count),
``fragment_count`` (every fragment it deposited on the same objects) and ``gdna_eff_len``. The EM reads the gDNA count
as a share of the fragment count (odds ``G : T - G``, ``pipeline.em_pseudocounts``), and the fragment count is each
object's whole tally, which no split of it moves, so the gDNA count's error is the whole error of the split. This
scores ``gdna_count`` in fragments (``Σ|Δ|``), per condition and per stratum, over three arms:

- ``P``: the shipped prior;
- ``O``: the same assembler fed the true per-object masses (``OracleTruth.override_masses``);
- ``S``: ``O`` with the crossing term converted by gDNA's own true share instead of the pooled one.

``P − O`` is calibration's error and ``O − S`` the pooled share's. ``O`` keeps ``P``'s ``gdna_eff_len``, so the
length is not scored here (``ruler_vs_truth.py`` scores it). Every arm runs in the drained frame, with the drain's
spliced-gDNA leak and the lift's attribution bound reported beside the numbers. The pipeline stops at the prior
stage; no per-locus EM runs.

Gate: ``tests/calibration/test_prior_vs_oracle.py``.

Usage::

    python scripts/design/prior_vs_oracle.py --suite DIR --oracle-cache DIR/oracle_cache --jobs 6
    python scripts/design/prior_vs_oracle.py --conditions gdna_g50_ss_0.50_nascent_mid_capture_on
    python scripts/design/prior_vs_oracle.py --index INDEX --work-dir SCRATCH
    python scripts/design/prior_vs_oracle.py --set scan.total_threads=1   # a pinned, reproducible run
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import os
import time
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402

from _shared import (  # noqa: E402
    DEFAULT_INDEX,
    DEFAULT_SUITE,
    ORACLE_BAM,
    OVERRIDE_FIELDS,
    SCRATCH,
    add_set_option,
    check_override_fields,
    configured,
    oracle_cache_dir,
    panel_conditions,
    run_json_shards,
    set_args,
    stratum,
    stratum_selections,
)

from _oracle import OracleTruth  # noqa: E402
from _oracle_arms import load_or_build_oracle  # noqa: E402

import rigel.calibration.priors as PRIORS  # noqa: E402
from rigel.calibration.calibrate import calibrate  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402
from rigel.pipeline import (  # noqa: E402
    drain_side_buffer,
    length_models,
    length_tables,
    quant_from_buffer,
    scan_and_buffer,
    transcript_geometry,
)

#: The ``LocusPriors`` fields: the noop gate holds every one byte-identical.
PRIOR_FIELDS = ("gdna_count", "fragment_count", "gdna_eff_len")

# ── the scoring ──────────────────────────────────────────────────────────────────────────────────


@dataclass(frozen=True, slots=True)
class ArmScore:
    """One arm's per-locus gDNA count scored against a reference's, over one locus selection. The fields add
    across a partition of the loci; the rates ``rel`` and ``cancellation`` are derived from them."""

    n_loci: int  #: loci in the selection, not the panel's locus count
    n_claiming: int  #: loci where the arm or the reference puts a nonzero gDNA count
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
        """``Σ|err| / |net|``. Large means under- and over-calls cancel, so a library-level summary flatters
        the arm."""
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


def capture_priors(buffer, index, strand_models, tables, region_arrays, stats, calibration,
                   junction_cap, geometry, pipeline_config):
    """Run the production quant path far enough to get ``(multi_loci, LocusPriors)``, then stop.

    The loci depend on the scoring stage, so they are taken from production's own call by wrapping
    :func:`~rigel.calibration.priors.assemble_priors` (imported function-locally by ``quant_from_buffer``, so the
    patch is seen at call time). A sentinel exception then skips the per-locus EM.
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
            buffer, index, strand_models, tables, region_arrays, stats, calibration, junction_cap, geometry,
            em_config=pipeline_config.em, scoring=pipeline_config.scoring,
        )
    except _StopAfterPriors:
        pass
    finally:
        PRIORS.assemble_priors = original
    if "priors" not in captured:
        # ``quant_from_buffer`` returns early on an empty unit set
        raise RuntimeError(
            "assemble_priors was never called — quant_from_buffer returned before the prior stage "
            "(no EM units, or no multi-loci). This is not a condition with zero error."
        )
    return captured["multi_loci"], captured["priors"]


def oracle_priors(oracle: OracleTruth, calibration, region_arrays, multi_loci):
    """O, the prior a perfect deconvolution would produce, and the ``noop`` gate beside it. Returns ``(O, noop)``.

    ``override_masses`` swaps the six mass arrays for the origin-split truth. ``noop`` re-injects the shipped
    masses through the identical ``dataclasses.replace`` and must come back byte-identical.
    """
    override = oracle.override_masses(region_arrays)
    check_override_fields(override)
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
    """S, the O arm with the gDNA count's crossing term converted by gDNA's own true per-boundary share.

    ``assemble_priors`` converts the crossing term by one pooled share, ``q = mass / count`` off the mixture, so
    where the components' shares differ the gDNA count is off by ``q / share_g`` and the locus total cannot show
    it. This runs the shipped function with ``boundary_mass_per_crossing`` set to gDNA's true share, measured off
    the origin split (``OracleTruth.component_shares``).
    """
    shares = oracle.component_shares()
    truth_cal = dataclasses.replace(calibration, **oracle.override_masses(region_arrays))
    return PRIORS.assemble_priors(
        dataclasses.replace(truth_cal, boundary_mass_per_crossing=shares["gdna"]),
        region_arrays, multi_loci,
    ), shares


# ── one condition ────────────────────────────────────────────────────────────────────────────────


@dataclass
class ConditionResult:
    """Everything one condition produced. The arrays are kept so a gate can re-score them."""

    condition: str
    n_loci: int
    priors: dict  #: arm name -> LocusPriors
    noop_identical: dict  #: field -> bool
    drain: dict  #: the measured drain caveat
    library: dict  #: condition-level totals, for the header row
    seconds: float
    #: the inputs a gate needs to re-derive any arm from the shipped code
    oracle: OracleTruth = None
    region_arrays: object = None
    multi_loci: list = None
    calibration: object = None


def _calibrate_and_prior(payload, strand_model, buffer, stats, index, ra, pipeline_config):
    """calibrate → score → ``(calibration, multi_loci, LocusPriors)`` on one payload."""
    from rigel.calibration.splice_graph import build_boundary_flags_array, build_sj_geometry_arrays
    from rigel.pipeline import library_fl_models

    fl = library_fl_models(payload, index, ra)
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
    # the EM's length tables and the junction cap, as production decides them
    tables = length_tables(fl, cal, ra)
    junction_cap, geometry = transcript_geometry(index, length_models(tables)[0], cal, ra, tables)
    multi_loci, priors = capture_priors(
        buffer, index, strand_model, tables, ra, stats, cal, junction_cap, geometry, pipeline_config,
    )
    return cal, multi_loci, priors


def measure_condition(bam, index, pipeline_config, work_dir, tag, *, oracle_cache=None) -> ConditionResult:
    """Scan once, drain to the production frame, build P / O / S."""
    start = time.perf_counter()
    scan = pipeline_config.scan
    ra = RegionArrays.from_index(index)

    stats, strand_model, buffer, payload = scan_and_buffer(bam, index, scan)
    n_held = int(payload.deferred.n_fragments)
    # drained at the production seed; the oracle partitions replay the whole's choices, and sum-to-full checks it
    lift: dict = {}
    payload = drain_side_buffer(
        payload, index, strand_model, regions=ra, seed=pipeline_config.second_pass_seed, lift=lift
    )
    # sum-to-full is re-run on the loaded parts, cached or rebuilt; it raises if the lift breaks a bank
    oracle = load_or_build_oracle(
        bam, index, pipeline_config, work_dir, tag, payload,
        Path(work_dir) / "oracle_cache" if oracle_cache is None else oracle_cache, lift=lift,
    )

    cal, multi_loci, p_arm = _calibrate_and_prior(
        payload, strand_model, buffer, stats, index, ra, pipeline_config
    )
    o_arm, noop = oracle_priors(oracle, cal, ra, multi_loci)
    s_arm, _shares = share_priors(oracle, cal, ra, multi_loci)

    noop_identical = {
        f: bool(np.array_equal(getattr(noop, f), getattr(p_arm, f))) for f in PRIOR_FIELDS
    }

    # ``gdna_spliced_leak`` is ``None`` only when nothing was held
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
        np.asarray(oracle.parts["annotated"].region_start_count, np.int64).sum()
        + np.asarray(oracle.parts["nascent"].region_start_count, np.int64).sum()
    )
    return ConditionResult(
        condition=tag,
        n_loci=len(multi_loci),
        priors={"P": p_arm, "O": o_arm, "S": s_arm},
        noop_identical=noop_identical,
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
    """Sum a list of ``ArmScore`` into one; the rates are re-derived from the summed totals, never averaged
    over conditions of different depth."""
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


def _rel(x: float) -> str:
    """``rel`` in 8 columns, scientific below 1e-3 so a small residual never prints as zero."""
    if not np.isfinite(x):
        return f"{'nan':>8}"
    return f"{x:>8.1e}" if 0.0 < abs(x) < 1e-3 else f"{x:>8.4f}"


def _check_noop_gate(rows: list[dict]) -> None:
    """Refuse the report when re-injecting the shipped masses did not reproduce ``P`` byte-identically."""
    bad = [(r["condition"], f) for r in rows for f, ok in r["noop_identical"].items() if not ok]
    if bad:
        print(f"\n  ⛔⛔ NOOP GATE FAILED on {len(bad)} (condition, field) pairs — the override "
              "plumbing is NOT inert and every number below is unattributable:")
        for c, f in bad[:10]:
            print(f"       {c}  {f}")
        raise SystemExit(2)
    print(f"\n  ✅ noop gate: re-injecting the shipped masses reproduces P byte-identically on all "
          f"{len(rows)} x {len(PRIOR_FIELDS)} arrays (TRAPS: byte-identity-gate)")


def _print_arm_table(rows: list[dict], title: str, key: str, note: str = "") -> None:
    """One arm pair's ``gdna_count`` scores, summed per stratum."""
    print()
    print(f"  {title}")
    if note:
        print(f"  {note}")
    print(f"    {'stratum':<26} {'ref total':>14} {'arm total':>14} {'Σ|Δ|':>14} "
          f"{'rel':>8} {'net':>14} {'canc':>7}")
    print("    " + "-" * 102)
    for label, sel in stratum_selections(r["condition"] for r in rows):
        s = _agg([ArmScore(**r[key]) for r in rows if sel(r["condition"])])
        if s is None:
            print(f"    {label:<26} {'(empty)':>14}")
            continue
        print(f"    {label:<26} {s.total_ref:>14,.0f} {s.total_arm:>14,.0f} "
              f"{s.abs_err:>14,.0f} {_rel(s.rel)} {s.net_err:>+14,.0f} "
              f"{s.cancellation:>7.1f}")


def _print_drained_frame(rows: list[dict]) -> None:
    """③: the lift's attribution bound and the drain's spliced-gDNA leak."""
    print()
    print("  ③ THE DRAINED FRAME — the lift's attribution bound, and production's certified-RNA leak")
    d = [r["drain"] for r in rows]
    amb, held = sum(x["n_ambiguous"] for x in d), sum(x["n_held"] for x in d)
    leak = sum(x["gdna_spliced_leak"] or 0 for x in d)
    n_leak = sum(1 for x in d if (x["gdna_spliced_leak"] or 0) > 0)
    print(f"     {held:,} held fragments drained at the production seed; lift ambiguity "
          f"{amb:,} ({amb / max(held, 1):.2%}) bounds the per-origin truth attribution.")
    print(f"     ⚠ the drain deposited {leak:,} spliced records into the gdna partition on "
          f"{n_leak}/{len(d)} conditions — production behaviour, RECORDED "
          "(issue #17), never an oracle refusal.")


def _print_per_condition(rows: list[dict]) -> None:
    """④: ``P`` and ``O`` per condition, since a stratum total hides the shape."""
    print()
    print("  ④ PER CONDITION — gdna_count, P and O as fragment totals")
    print(f"    {'condition':<44} {'true f_g':>9} {'O_g':>13} {'P_g':>13} {'P/O':>8} "
          f"{'Σ|Δ|':>13} {'rel':>8}")
    print("    " + "-" * 116)
    for r in rows:
        g = ArmScore(**r["P_vs_O"])
        ratio = g.total_arm / g.total_ref if g.total_ref > 0 else float("nan")
        print(f"    {r['condition']:<44} {r['library']['true_f_gdna']:>9.4f} {g.total_ref:>13,.0f} "
              f"{g.total_arm:>13,.0f} {ratio:>8.3f} {g.abs_err:>13,.0f} {_rel(g.rel)}")


def report(rows: list[dict], settings=()) -> None:
    """The whole report, from the per-condition JSON rows; serial and sharded runs share this one path."""
    rows = sorted(rows, key=lambda r: r["condition"])
    print()
    print("=" * 104)
    print("  ⭐⭐⭐ CALIBRATION'S ENDPOINT vs THE ORACLE — LocusPriors, in FRAGMENTS")
    config = " ".join(f"--set {spec}" for spec in settings) or "the shipped defaults"
    print(f"  {len(rows)} conditions   config: {config}   the drained frame (the drain's leak is "
          "reported below)")
    print("=" * 104)
    _check_noop_gate(rows)
    _print_arm_table(rows, "① P vs O — CALIBRATION'S OWN ERROR (a perfect deconvolution, same assembler) · gdna_count",
                     "P_vs_O")
    _print_arm_table(rows, "② O vs S — THE POOLED SHARE'S OWN CONTRIBUTION, ISOLATED · gdna_count",
                     "O_vs_S",
                     "⛔ What the pooled share moves off gDNA it moves onto RNA, so no gate on the locus total "
                     "can see it (issue #28).")
    _print_drained_frame(rows)
    _print_per_condition(rows)
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
    ap.add_argument("--work-dir", type=Path, default=SCRATCH / "rigel_prior_oracle")
    ap.add_argument("--json", type=Path, default=None)
    ap.add_argument("--jobs", type=int, default=1,
                    help="run this many conditions CONCURRENTLY by re-invoking on shards. The "
                         "conditions are independent, so this changes no number.")
    add_set_option(ap, "override one PipelineConfig field, typed from the field; repeatable "
                       "(e.g. --set scan.total_threads=1 for a reproducible run)")
    args = ap.parse_args()

    names = panel_conditions(args.suite, args.conditions)
    if not names:
        raise SystemExit(f"no conditions with a sim_oracle.bam under {args.suite}")
    cache = oracle_cache_dir(args.suite, args.oracle_cache)

    if args.jobs > 1 and len(names) > 1:
        # each shard works in a directory beside its output
        payload = run_json_shards(
            Path(__file__).resolve(), names, args.jobs, args.work_dir,
            lambda out, shard: ["--suite", str(args.suite), "--index", str(args.index),
                                "--work-dir", str(out.with_suffix("")), "--json", str(out),
                                "--conditions", *shard,
                                *([] if cache is None else ["--oracle-cache", str(cache)]),
                                *set_args(args.settings)],
        )
    else:
        index = TranscriptIndex.load(str(args.index))
        pipeline_config = configured(args.settings)
        args.work_dir.mkdir(parents=True, exist_ok=True)
        results = []
        for name in names:
            print(f"  … {name}", flush=True)
            results.append(measure_condition(
                str(args.suite / name / ORACLE_BAM), index, pipeline_config, args.work_dir, name,
                oracle_cache=cache,
            ))
        payload = to_json(results)
    if args.json is not None:
        args.json.write_text(json.dumps(payload, indent=1))
    report(payload, args.settings)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
