"""Re-scan a panel's oracle caches after a deposit-rule change, and GATE the rebuild on byte-identity.

⭐⭐⭐ **WHY THIS EXISTS.** Re-scanning is the first IRREVERSIBLE step of a deposit-rule change: the old
caches are overwritten and cannot be recovered. Everything before it is a code revert. So the rebuild
carries its own falsification — for every condition, every bank the change was NOT supposed to touch must
come back **byte-identical**, and the symmetric difference must be the SAME one rename on every condition.

    old cache          fresh scan          verdict
    ---------          ----------          -------
    31 shared banks    31 shared banks     ⛔ every one byte-identical, or STOP
    bank A only        bank B only         ⛔ the SAME (A, B) on every condition, or STOP

⛔ **THE ORDERING IS THE WHOLE POINT.** ``write_scan_cache`` overwrites ``payload.npz`` in place, so the
comparison baseline is destroyed by the very act of rebuilding. The old arrays are therefore read into
memory BEFORE the scan runs, and a condition that fails the gate is reported without its cache being
written — the stale one survives for inspection.

⛔ **NO EXPECTED DELTA IS HARD-CODED.** The tool does not know which bank was renamed; it DERIVES the
delta from the first condition it rebuilds and then requires every later condition to reproduce it
exactly. A side effect that appears only under capture, or only at ``g98``, is a mismatch against the
first condition rather than a mismatch against a constant nobody re-checks.

⚠ The manifest's deposit-sensitive scalars are compared too — ``qc`` (``deposited``, every ``dropped_*``)
and ``gap_resolution`` — because an arbitration change moves those while leaving every array plausible.
``payload_schema_digest`` is excluded: moving it is the point.

    python scripts/design/rescan_panels.py --index IDX --suite SUITE [--conditions A B] [--jobs N]
    python scripts/design/rescan_panels.py --self-test        # perturbs the comparator, no I/O

⚠ Split BAMs are deleted after their scan unless ``--keep-splits``; the shipped oracle path leaves them,
which is ~800 MB per condition and ~35 GB across a 44-condition rebuild.
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tests" / "calibration"))

from _oracle import ORIGINS, _split_bam  # noqa: E402

from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402
from rigel.pipeline import _native_detect_sj_tag, scan_and_buffer  # noqa: E402
from rigel.scan_cache import write_scan_cache  # noqa: E402

#: The four payloads cached per condition: the production scan plus the origin split it is validated
#: against. ``_main`` first so a condition's cheapest check runs before its 105 s split.
PAYLOADS = ("_main", *ORIGINS)

#: Manifest sub-dictionaries that a DEPOSIT or ARBITRATION change moves. Compared alongside the arrays.
DEPOSIT_SENSITIVE_SCALARS = ("qc", "gap_resolution")


# ────────────────────────────────────────────────────────────────────────────────────────────────────
# The comparator — pure, and the only thing --self-test exercises.
# ────────────────────────────────────────────────────────────────────────────────────────────────────


def compare(old_arrays, old_scalars, new_arrays, new_scalars, expect_changed=()):
    """Compare two payload snapshots. Returns ``(only_old, only_new, differing, changed_as_expected)``.

    ``differing`` lists every shared key whose value is not byte-identical — arrays by
    :func:`numpy.array_equal` after a shape check, scalars by equality. ⛔ Byte-identity, never a
    tolerance: these are integer tallies and the whole claim is that they did not move.

    ⭐ ``expect_changed`` names the banks a DEPOSIT-RULE change is supposed to move. They are reported
    in ``changed_as_expected`` instead of ``differing``, so the gate keeps its real teeth — *nothing
    moved EXCEPT what was named* — rather than being switched off for the whole rebuild.
    ⛔ A named bank that did NOT move is reported too, and it is the more interesting failure: it means
    the rule change did not reach the data (`TRAPS: could-the-arm-have-fired`).
    """
    ko, kn = set(old_arrays), set(new_arrays)
    expect_changed = set(expect_changed)
    differing, changed_as_expected = [], []
    for k in sorted((ko & kn) & expect_changed):
        a, b = old_arrays[k], new_arrays[k]
        if a.dtype != b.dtype:
            # ⭐ A dtype change IS a change, and a numeric-convention change is exactly that. Reporting
            # it as "unchanged" would fail the gate for the one reason it was told to expect.
            changed_as_expected.append(f"{k}: dtype {a.dtype} -> {b.dtype} (expected)")
        elif a.shape != b.shape:
            changed_as_expected.append(f"{k}: shape {a.shape} -> {b.shape} (expected)")
        elif not np.array_equal(a, b):
            n = int((a != b).sum())
            changed_as_expected.append(f"{k}: {n} of {a.size} elements moved (expected)")
        else:
            changed_as_expected.append(f"{k}: ⛔ UNCHANGED — the rule change did not reach this bank")
    for k in sorted((ko & kn) - expect_changed):
        a, b = old_arrays[k], new_arrays[k]
        if a.shape != b.shape:
            differing.append(f"{k}: shape {a.shape} -> {b.shape}")
        elif a.dtype != b.dtype:
            differing.append(f"{k}: dtype {a.dtype} -> {b.dtype}")
        elif not np.array_equal(a, b):
            n = int((a != b).sum())
            differing.append(f"{k}: {n} of {a.size} elements differ")
    for group in DEPOSIT_SENSITIVE_SCALARS:
        a, b = old_scalars.get(group), new_scalars.get(group)
        if a != b:
            differing.append(f"{group}: {a!r} -> {b!r}")
    return sorted(ko - kn), sorted(kn - ko), differing, changed_as_expected


def _snapshot(cache_dir: Path):
    """Read a cache's arrays and deposit-sensitive scalars into memory. ``None`` if it is not there.

    ⛔ Arrays are materialised (``np.load`` is lazy over the zip), because the file is about to be
    overwritten by the rebuild and a lazy handle would then read the NEW bytes.
    """
    npz, manifest = cache_dir / "payload.npz", cache_dir / "manifest.json"
    if not npz.is_file() or not manifest.is_file():
        return None
    with np.load(npz) as z:
        arrays = {k: np.array(z[k]) for k in z.files}
    scalars = json.loads(manifest.read_text()).get("payload_scalars", {})
    return arrays, {g: scalars.get(g) for g in DEPOSIT_SENSITIVE_SCALARS}


def _fresh(payload):
    """The same two views, taken from a payload object rather than from disk."""
    import dataclasses

    from rigel.scan_payload import AccumulatorPayload

    arrays, scalars = {}, {}
    for field in dataclasses.fields(AccumulatorPayload):
        value = getattr(payload, field.name)
        if isinstance(value, np.ndarray):
            arrays[field.name] = np.ascontiguousarray(value)
        elif dataclasses.is_dataclass(value):
            nested = {}
            for sub in dataclasses.fields(value):
                sub_value = getattr(value, sub.name)
                if isinstance(sub_value, np.ndarray):
                    arrays[f"{field.name}__{sub.name}"] = np.ascontiguousarray(sub_value)
                else:
                    nested[sub.name] = sub_value
            scalars[field.name] = nested
        else:
            scalars[field.name] = value
    return arrays, {g: scalars.get(g) for g in DEPOSIT_SENSITIVE_SCALARS}


# ────────────────────────────────────────────────────────────────────────────────────────────────────
# The rebuild
# ────────────────────────────────────────────────────────────────────────────────────────────────────


def rebuild_condition(index, suite: Path, condition: str, work_dir: Path, keep_splits: bool,
                      flat: bool = False, expect_changed=()):
    """Rebuild one condition's payloads, gating each against the cache it replaces.

    Returns a dict per payload: ``{only_old, only_new, differing, had_baseline, seconds}``. ⛔ A payload
    whose ``differing`` is non-empty is NOT written — the stale cache survives so the difference can be
    read off disk afterwards.

    ⭐ TWO CACHE LAYOUTS, one rebuild. The oracle layout is ``oracle_cache/<cond>/{_main,gdna,mrna,nrna}``
    and needs the origin split; ``flat`` is ``scan_cache/<cond>/`` holding the production payload alone,
    which is what the pilot panel uses. Only the directory shape and the payload list differ, so they
    share the gate rather than growing a second tool with a second comparator.
    """
    bam = str(suite / condition / "sim_oracle.bam")
    payloads = ("_main",) if flat else PAYLOADS
    cache_root = (suite / "scan_cache" / condition) if flat else (suite / "oracle_cache" / condition)
    scan = dataclasses.replace(
        PipelineConfig().scan, sj_strand_tag=_native_detect_sj_tag(bam)
    )

    out: dict[str, dict] = {}
    split_paths = None
    try:
        for name in payloads:
            t0 = time.perf_counter()
            if name == "_main":
                source = bam
            else:
                if split_paths is None:
                    split_paths, _counts = _split_bam(bam, work_dir, condition)
                source = split_paths[name]

            cache_dir = cache_root if flat else cache_root / name
            baseline = _snapshot(cache_dir)  # ⛔ BEFORE the write, or the baseline is the new bytes

            _stats, strand_model, _buffer, payload = scan_and_buffer(source, index, scan)
            new_arrays, new_scalars = _fresh(payload)

            if baseline is None:
                record = {"only_old": [], "only_new": [], "differing": [], "expected": [],
                          "had_baseline": False}
            else:
                only_old, only_new, differing, expected = compare(
                    *baseline, new_arrays, new_scalars, expect_changed
                )
                record = {
                    "only_old": only_old,
                    "only_new": only_new,
                    "differing": differing,
                    "expected": expected,
                    "had_baseline": True,
                    "n_shared": len(set(baseline[0]) & set(new_arrays)),
                    # ⛔ A COMPARISON AGAINST AN ALL-ZERO BASELINE PROVES NOTHING, and this panel
                    # manufactures them: ``g00``'s gdna partition and every ``nrna_none`` condition's
                    # nrna partition are empty BY CONSTRUCTION. Zeros match zeros whatever the deposit
                    # rule is, so those payloads are counted as VACUOUS, never as passes
                    # (`TRAPS: could-the-arm-have-fired` — a comparison that could not have differed is not a control).
                    "live_banks": sum(
                        1 for k, v in baseline[0].items() if v.size and np.any(v)
                    ),
                }
            record["seconds"] = round(time.perf_counter() - t0, 1)

            if not record["differing"]:
                write_scan_cache(
                    cache_dir,
                    payload=payload,
                    strand_model=strand_model,
                    index=index,
                    bam=source,
                    scan_config=scan,
                )
                record["written"] = True
            else:
                record["written"] = False
            out[name] = record
    finally:
        if split_paths is not None and not keep_splits:
            for p in split_paths.values():
                Path(p).unlink(missing_ok=True)
    return out


def _delta_key(record) -> tuple:
    """The (removed, added) bank sets a condition reported — the thing every condition must agree on."""
    return (tuple(record["only_old"]), tuple(record["only_new"]))


# ────────────────────────────────────────────────────────────────────────────────────────────────────
# --self-test: PERTURB THE COMPARATOR. A gate nobody broke is a gate nobody knows fires.
# ────────────────────────────────────────────────────────────────────────────────────────────────────


def self_test() -> int:
    base_a = {"x": np.arange(10, dtype=np.uint64), "y": np.ones(4, dtype=np.float64)}
    base_s = {"qc": {"deposited": 7}, "gap_resolution": {"a": 1}}

    def fresh():
        return ({k: v.copy() for k, v in base_a.items()}, json.loads(json.dumps(base_s)))

    checks: list[tuple[str, bool]] = []

    # ✅ identical must be silent — the gate has to be able to PASS, or it proves nothing when it does
    o, n, d, _e = compare(*fresh(), *fresh())
    checks.append(("identical => no difference", (o, n, d) == ([], [], [])))

    # ⛔ a 1-ULP nudge in a float bank
    a2, s2 = fresh()
    a2["y"][2] = np.nextafter(a2["y"][2], 2.0)
    o, n, d, _e = compare(*fresh(), a2, s2)
    checks.append(("1-ULP float nudge => caught", len(d) == 1 and d[0].startswith("y:")))

    # ⛔ a single count off by one
    a2, s2 = fresh()
    a2["x"][0] += 1
    o, n, d, _e = compare(*fresh(), a2, s2)
    checks.append(("integer +1 => caught", len(d) == 1 and d[0].startswith("x:")))

    # ⛔ a shape change and a dtype change are not value comparisons and must be caught separately
    a2, s2 = fresh()
    a2["x"] = np.arange(11, dtype=np.uint64)
    o, n, d, _e = compare(*fresh(), a2, s2)
    checks.append(("shape change => caught", len(d) == 1 and "shape" in d[0]))

    a2, s2 = fresh()
    a2["x"] = base_a["x"].astype(np.int64)
    o, n, d, _e = compare(*fresh(), a2, s2)
    checks.append(("dtype change => caught", len(d) == 1 and "dtype" in d[0]))

    # ⛔ an added and a removed bank land in the symmetric difference, not in `differing`
    a2, s2 = fresh()
    del a2["y"]
    a2["z"] = np.zeros(3)
    o, n, d, _e = compare(*fresh(), a2, s2)
    checks.append(("bank removed+added => symmetric difference", (o, n, d) == (["y"], ["z"], [])))

    # ⛔ a deposit-sensitive SCALAR moving while every array holds — the arbitration failure mode
    a2, s2 = fresh()
    s2["qc"]["deposited"] = 8
    o, n, d, _e = compare(*fresh(), a2, s2)
    checks.append(("qc scalar moved => caught", len(d) == 1 and d[0].startswith("qc:")))

    a2, s2 = fresh()
    s2["gap_resolution"]["a"] = 2
    o, n, d, _e = compare(*fresh(), a2, s2)
    checks.append(("gap_resolution moved => caught", len(d) == 1 and d[0].startswith("gap_res")))

    # ⛔ an EMPTY baseline must not read as agreement (arm_identity.py's recorded lie: zero rows
    # once scored "32/32 IDENTICAL")
    o, n, d, _e = compare({}, base_s, *fresh())
    checks.append(("empty baseline => every bank is `only_new`", n == ["x", "y"] and d == []))

    # ⭐ --expect-changed: a NAMED bank that moved is expected, not a failure — but the gate must keep
    # its teeth for every OTHER bank in the same comparison.
    a2, s2 = fresh()
    a2["x"][0] += 1
    a2["y"][1] += 1
    o, n, d, e = compare(*fresh(), a2, s2, expect_changed={"x"})
    checks.append((
        "expect-changed: named bank moves => expected, unnamed still fails",
        len(e) == 1 and e[0].startswith("x:") and "expected" in e[0]
        and len(d) == 1 and d[0].startswith("y:"),
    ))

    # ⛔⛔ AND THE ONE THAT MATTERS MOST: a bank you SAID would move and did not. That is the rule
    # change failing to reach the data, and it must not read as a clean pass.
    o, n, d, e = compare(*fresh(), *fresh(), expect_changed={"x"})
    checks.append((
        "expect-changed: named bank did NOT move => reported UNCHANGED",
        len(e) == 1 and "UNCHANGED" in e[0] and d == [],
    ))

    width = max(len(name) for name, _ in checks)
    for name, ok in checks:
        print(f"  {'PASS' if ok else 'FAIL'}  {name:<{width}}", flush=True)
    failed = [name for name, ok in checks if not ok]
    print(f"\n{len(checks) - len(failed)}/{len(checks)} comparator perturbations fire", flush=True)
    return 1 if failed else 0


# ────────────────────────────────────────────────────────────────────────────────────────────────────


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--index", type=Path, default=None)
    ap.add_argument("--suite", type=Path, default=None)
    ap.add_argument("--conditions", nargs="*", default=None)
    ap.add_argument("--work-dir", type=Path, default=Path(os.environ.get("RIGEL_SCRATCH", "/tmp")))
    ap.add_argument("--json", type=Path, default=None)
    ap.add_argument("--keep-splits", action="store_true")
    ap.add_argument(
        "--expect-changed", nargs="*", default=[],
        help="banks a DEPOSIT-RULE change is supposed to move. Named explicitly so the gate keeps its "
             "teeth: nothing may move EXCEPT these, and one of these that does NOT move is reported "
             "too, because that means the rule change never reached the data.",
    )
    ap.add_argument(
        "--flat", action="store_true",
        help="the pilot layout: scan_cache/<cond>/ holding the production payload alone, with no "
             "origin split. Same gate, different directory shape.",
    )
    ap.add_argument(
        "--jobs", type=int, default=1,
        help="rebuild this many conditions CONCURRENTLY by re-invoking this script on shards. The "
             "conditions are independent BAMs, so this changes no number.",
    )
    ap.add_argument("--self-test", action="store_true", help="perturb the comparator; no I/O")
    args = ap.parse_args()

    if args.self_test:
        return self_test()
    if args.index is None or args.suite is None:
        ap.error("--index and --suite are required unless --self-test")

    names = args.conditions or sorted(
        p.name for p in args.suite.iterdir() if (p / "sim_oracle.bam").is_file()
    )

    if args.jobs > 1 and len(names) > 1:
        shards = [names[i :: args.jobs] for i in range(args.jobs)]
        shards = [s for s in shards if s]
        # ⛔ PER-INVOCATION, not a shared name. Two concurrent runs sharing a `--work-dir` used the same
        # `rescan_shards/`, and the first to finish `rmtree`d it out from under the other's shards —
        # which then died writing their JSON. The rebuild itself was unharmed (a payload is written only
        # when it compares identical), but every shard's VERDICT was lost, so 8 conditions were
        # re-scanned with no gate report. The isolation is the fix; the pid makes it unique.
        tmp = args.work_dir / f"rescan_shards_{os.getpid()}"
        tmp.mkdir(parents=True, exist_ok=True)
        procs = []
        for i, shard in enumerate(shards):
            out = tmp / f"shard{i}.json"
            cmd = [sys.executable, __file__, "--index", str(args.index), "--suite", str(args.suite),
                   "--work-dir", str(args.work_dir), "--json", str(out), "--conditions", *shard]
            if args.keep_splits:
                cmd.append("--keep-splits")
            if args.flat:
                cmd.append("--flat")
            if args.expect_changed:
                cmd += ["--expect-changed", *args.expect_changed]
            procs.append((subprocess.Popen(cmd), out))
        results: dict[str, dict] = {}
        lost = []
        for proc, out in procs:
            proc.wait()
            if out.is_file():
                results.update(json.loads(out.read_text()))
            else:
                lost.append(out.name)
        shutil.rmtree(tmp, ignore_errors=True)
        # ⛔ A SHARD THAT LEFT NO RESULT IS NOT A SHARD THAT PASSED. Its conditions may well have been
        # re-scanned — and if they were, they compared identical, because nothing else is written — but
        # the verdict is gone, and a silently short aggregate would read as a clean run over fewer
        # conditions. Name them so the operator can re-gate or accept the loss knowingly.
        if lost:
            print(f"⛔ {len(lost)} shard(s) produced no result file ({', '.join(lost)}). Their "
                  f"conditions are NOT covered by the verdict below.", flush=True)
    else:
        index = TranscriptIndex.load(str(args.index))
        args.work_dir.mkdir(parents=True, exist_ok=True)
        results = {}
        for name in names:
            results[name] = rebuild_condition(
                index, args.suite, name, args.work_dir, args.keep_splits, args.flat,
                args.expect_changed,
            )
            if args.json is None:
                _report_one(name, results[name])

    if args.json is not None:
        args.json.write_text(json.dumps(results, indent=1))

    return _verdict(results)


def _report_one(condition: str, payloads: dict) -> None:
    secs = sum(p["seconds"] for p in payloads.values())
    print(f"\n{condition}   ({secs:.0f}s)", flush=True)
    for name, r in payloads.items():
        if not r["had_baseline"]:
            print(f"  {name:<6} no prior cache — built, nothing to compare", flush=True)
            continue
        state = "⛔ DIFFERS" if r["differing"] else "identical"
        print(f"  {name:<6} {r['n_shared']:>2} shared banks {state}"
              f"   -{r['only_old']} +{r['only_new']}", flush=True)
        for d in r["differing"]:
            print(f"           {d}", flush=True)
        for e in r.get("expected", ()):
            print(f"           {e}", flush=True)


def _verdict(results: dict) -> int:
    """⛔ THE GATE. Every shared bank identical, and the SAME rename delta on every condition."""
    print("\n" + "=" * 96, flush=True)
    bad_arrays, deltas = [], {}
    compared = vacuous = 0
    for condition, payloads in sorted(results.items()):
        for name, r in payloads.items():
            if r["differing"]:
                bad_arrays.append(f"{condition}/{name}: {'; '.join(r['differing'])}")
            if r["had_baseline"]:
                compared += 1
                if not r.get("live_banks"):
                    vacuous += 1
                deltas.setdefault(_delta_key(r), []).append(f"{condition}/{name}")

    inert = [f"{c}/{n}: {e}" for c, ps in sorted(results.items()) for n, r in ps.items()
             for e in r.get("expected", ()) if "UNCHANGED" in e]

    n_written = sum(p.get("written", False) for ps in results.values() for p in ps.values())
    print(f"{len(results)} conditions, {n_written} payloads written, {compared} compared to a baseline",
          flush=True)
    if vacuous:
        print(f"⚠ {vacuous} of those baselines were ALL-ZERO (g00's gdna, every nrna_none's nrna) — "
              f"zeros match zeros whatever the deposit rule is. {compared - vacuous} comparisons "
              f"carry the verdict.", flush=True)

    if compared - vacuous == 0:
        print("⛔ NO PAYLOAD WITH A NON-EMPTY BASELINE — this run proves nothing. An empty comparison "
              "is not a pass.", flush=True)
        return 1

    print(f"\nrename deltas observed ({len(deltas)} distinct):", flush=True)
    for (removed, added), where in sorted(deltas.items(), key=lambda kv: -len(kv[1])):
        print(f"  -{list(removed)} +{list(added)}   on {len(where)} payloads", flush=True)
        if len(where) <= 3:
            for w in where:
                print(f"      {w}", flush=True)

    ok = True
    if bad_arrays:
        ok = False
        print(f"\n⛔ {len(bad_arrays)} payloads have a bank that MOVED — the deposit change had a side "
              f"effect. These were NOT written:", flush=True)
        for b in bad_arrays[:20]:
            print(f"   {b}", flush=True)
    if inert:
        ok = False
        print(f"\n⛔ {len(inert)} payload(s) have a bank named in --expect-changed that did NOT move. "
              f"The rule change did not reach the data (`TRAPS: could-the-arm-have-fired`):", flush=True)
        for i in inert[:10]:
            print(f"   {i}", flush=True)
    if len(deltas) > 1:
        ok = False
        print("\n⛔ the rename delta is NOT the same on every condition — a bank appeared or vanished "
              "on some conditions only.", flush=True)

    print("\n" + ("✅ GATE PASSED — every shared bank byte-identical, one delta everywhere"
                  if ok else "⛔ GATE FAILED — STOP, do not treat these caches as valid"), flush=True)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
