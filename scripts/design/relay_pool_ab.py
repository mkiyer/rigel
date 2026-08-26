"""THE RELAY A/B, PER POOL, PER OBJECT, AGAINST ORIGIN-SPLIT TRUTH — message propagation OFF vs ON.

⭐⭐⭐ **WHAT REGRESSIONS DOES MESSAGE PROPAGATION CAUSE, AND WHERE?** One row per condition per arm,
**NEVER COLLAPSED** — the panel total hides a sign flip between strata
(TRAPS: never-pool-the-strata), and an A/B whose only output is a total cannot say which stratum moved.

**THE POOLS.** Truth comes from the oracle cache's origin partitions, which are three separate scans::

    gdna    genomic DNA
    nrna    UNANNOTATED SYNTHETIC NASCENT RNA — its own pool, and the one the panel held at 0
    mrna    annotated RNA

⛔⛔ **CALIBRATION ESTIMATES TWO POOLS, NOT THREE, AND THAT IS AXIOM 0 RATHER THAN A GAP IN THIS FILE.**
The solver's populations are ``{gDNA, RNA+, RNA-}``; "nascent" and "annotated" are not species and are
not a degree of freedom there. So the ESTIMATE columns are gDNA and RNA, while the TRUTH columns carry
all three — the nascent/annotated split is an EM ENTITY-level quantity and is table B's
(``--table pipeline``), not this one's. Printing the three truth pools here anyway is deliberate: it is
what tells you whether a gDNA over-call is eating nascent RNA or annotated RNA.

**THE TWO ERRORS, and they answer different questions.**

===========  ====================================================================================
``net``      ``sum_i (est_i - true_i)`` over every object — the SIGNED total. A net near zero with
             a large absolute error means the tool is misplacing mass, not mis-sizing the pool
``abs``      ``sum_i |est_i - true_i|`` — the total misplaced mass, in FRAGMENTS. This is the one
             that cannot cancel, and it is the one to rank on
===========  ====================================================================================

⚠ **`nrna = 0` ON A CONDITION IS NOT A MEASUREMENT OF THE NASCENT POOL.** The 16-condition ladder held
``nrna: ratios: [0.0]`` on every row until 2026-08-18, so every nascent number it produced was that
zero restated. Rows whose truth nascent pool is exactly 0 are STAMPED, so a reader cannot mistake an
empty pool for a pool that was measured and found small.

⛔ Both arms are run IN ONE PROCESS off the SAME cached payload, so the only thing that differs is
``CalibrationConfig.message_propagation`` — no re-scan, no second truth source, no reseeding.
⛔ The relay arm ASSERTS IT RAN: ``_uni`` is written only under ``RelayPolicy``, and the muted arm must
reproduce ``f_g == fg_loc`` exactly. An arm that silently did not switch is
TRAPS: an-ablation-that-never-ran, which has already cost this project a 314-second run reported as
"all arms byte-identical".

Usage::

    python scripts/design/relay_pool_ab.py                       # the whole ladder, both arms
    python scripts/design/relay_pool_ab.py --conditions NAME ...
    python scripts/design/relay_pool_ab.py --self-test           # no I/O

⚠ SINCE 2026-08-25 the certified-flux stream is a MESSAGE (relay-only): the "off" arm
(SilentPolicy) carries NO anchor, so the off/on delta includes the stream and is not comparable
with numbers recorded before commit cb2268f1.
"""

from __future__ import annotations

import argparse
import dataclasses
import importlib.util
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402


def _sibling(name: str):
    key = name[:-3]
    if key not in sys.modules:
        spec = importlib.util.spec_from_file_location(key, Path(__file__).resolve().parent / name)
        module = importlib.util.module_from_spec(spec)
        sys.modules[key] = module
        spec.loader.exec_module(module)
    return sys.modules[key]


OC = _sibling("object_composition.py")

from rigel.calibration import calibrate  # noqa: E402
from rigel.calibration.region_chain import BOUNDARY, REGION, build_region_chain  # noqa: E402
from rigel.config import CalibrationConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402
from rigel.scan_cache import calibration_inputs, read_scan_cache  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tests"))
from calibration._oracle import ORIGINS, OracleTruth  # noqa: E402

DEFAULT_SUITE = OC.DEFAULT_SUITE
DEFAULT_INDEX = OC.DEFAULT_INDEX
_EPS = 1.0e-12

#: The three ORIGIN pools, in report order. ⭐ ``nrna`` is the unannotated synthetic nascent pool and
#: is deliberately named apart from ``mrna``: they are different entities to the EM, and lumping them
#: is what makes a nascent regression invisible.
POOLS = ("gdna", "nrna", "mrna")


def pool_truth(parts, region_arrays, chain) -> dict[str, np.ndarray]:
    """Per-slot TRUE count for each origin pool, on the chain the solver works in.

    ⭐ ``slot_counts`` projects ``region_contained`` at a REGION and ``boundary_unspliced`` at a
    BOUNDARY — exactly the mixture ψ deconvolves — so a truth built from the origin partitions and an
    estimate built from the full payload are on ONE basis. Any other projection would compare two
    different populations and call the difference an error.
    """
    return {k: OC.slot_counts(parts[k], region_arrays, chain) for k in POOLS}


def arm(payload, kw, *, messages: bool, policy: str = "relay", injected_priors=None) -> np.ndarray:
    """One arm's per-slot ``f_g``. ⛔ Returns the FINAL belief, which is what `assemble_priors` reads.

    ⭐ ``injected_priors`` is for a TOY reference — a chromosome small enough to have no library-level
    statistics of its own. `calibrate` REFUSES a library with zero spliced unique mappers (*"a real
    RNA-seq library always carries spliced reads"*), which is right about real data and is exactly the
    state of the method-development test chromosome while every transcript on it is single-exon
    (`docs/TESTING.md` §0a). Harvest them from a cached donor condition — `toy_harness.harvest` — and
    the toy is measured against ITS OWN truth with the library's parameters supplied, which is the
    contract the toy harness already runs under.
    """
    cfg = dataclasses.replace(
        CalibrationConfig(), message_propagation=messages, message_policy=policy
    )
    debug: dict = {}
    call = {k: v for k, v in kw.items() if k != "payload"}
    if injected_priors is not None:
        call["injected_priors"] = injected_priors
    calibrate(payload=payload, config=cfg, _debug=debug, **call)
    cap = debug["capture"]
    # ⛔ THE ARM RAN, AND BOTH DIRECTIONS ARE CHECKED. `_uni` is written only at messages/relay.py,
    #    i.e. only under RelayPolicy; and muted, ψ carries each slot's own evidence alone, so the
    #    final belief must BE the message-free local solve, bit for bit.
    relay_ran = "_uni" in cap
    currency_ran = "_currency" in cap
    if relay_ran != (messages and policy == "relay"):
        raise AssertionError(
            f"messages={messages} policy={policy} but the relay {'ran' if relay_ran else 'did not run'} "
            "— this arm is not the arm it claims to be (`_uni` is written only under RelayPolicy)."
        )
    if currency_ran != (messages and policy == "currency"):
        raise AssertionError(
            f"messages={messages} policy={policy} but the currency policy "
            f"{'ran' if currency_ran else 'did not run'} — this arm is not the arm it claims to be "
            "(`_currency` is written only under CurrencyPolicy)."
        )
    f_g = np.asarray(cap["f_g"], np.float64)
    if not messages:
        loc = np.asarray(cap["fg_loc"], np.float64)
        if not np.array_equal(f_g, loc):
            raise AssertionError(
                "messages=False but the final belief differs from the message-free local solve — "
                "SilentPolicy delivered something, so the muted arm is not muted."
            )
    return f_g


def score(f_g: np.ndarray, truth: dict[str, np.ndarray], sel: np.ndarray) -> dict:
    """The two errors, per estimated pool, over one axis of the chain.

    ⛔ ``net`` is SIGNED and ``abs`` is not, and reporting only one of them is how a misplacement is
    read as a sizing error (or the reverse). A tool that puts every fragment in the wrong object but
    the right pool has ``net = 0``.
    """
    true_g = truth["gdna"][sel]
    true_r = truth["mrna"][sel] + truth["nrna"][sel]
    mass = true_g + true_r
    est_g = f_g[sel] * mass
    est_r = (1.0 - f_g[sel]) * mass
    out = {"n_obj": int(sel.sum()), "mass": float(mass.sum())}
    for name, est, true in (("gdna", est_g, true_g), ("rna", est_r, true_r)):
        d = est - true
        out[f"{name}_est"] = float(est.sum())
        out[f"{name}_true"] = float(true.sum())
        out[f"{name}_net"] = float(d.sum())
        out[f"{name}_abs"] = float(np.abs(d).sum())
    return out


#: label -> (message_propagation, message_policy). ⭐ "off"/"on" are the standing benchmark's two
#: arms; "currency" is the Stage-3 policy under development, opt-in via --arms so every existing
#: caller and artifact keeps its exact shape.
ARMS: dict[str, tuple[bool, str]] = {
    "off": (False, "relay"),
    "on": (True, "relay"),
    "currency": (True, "currency"),
}


def measure(index, region_arrays, suite: Path, oracle_cache: Path, condition: str,
            injected_priors=None, arms: tuple[str, ...] = ("off", "on")) -> list[dict]:
    """One condition, the requested arms, three axes. ⛔ ONE cached payload, so the arms differ by
    exactly the config values `ARMS` names."""
    cache = read_scan_cache(Path(suite) / "scan_cache" / condition, index)
    kw = calibration_inputs(cache, index)
    chain = build_region_chain(
        cache.payload.ref_region_offsets, cache.payload.ref_boundary_offsets
    )
    root = Path(oracle_cache) / condition
    parts = {k: read_scan_cache(root / k, index).payload for k in ORIGINS}
    # ⭐ sum-to-full runs as a HARD gate on every condition — a cached partition that does not
    #   reconstruct the scan's own read is a silently wrong truth, invisible in every number below.
    OracleTruth.from_parts(cache.payload, parts)
    truth = pool_truth(parts, region_arrays, chain)

    kind = np.asarray(chain.kind)
    axes = {
        "REGION": kind == REGION,
        "BOUNDARY": kind == BOUNDARY,
        "ALL": np.ones(int(chain.n_slots), bool),
    }
    rows = []
    for label in arms:
        messages, policy = ARMS[label]
        f_g = arm(cache.payload, kw, messages=messages, policy=policy,
                  injected_priors=injected_priors)
        for axis, sel in axes.items():
            row = {"condition": condition, "arm": label, "axis": axis}
            row.update(score(f_g, truth, sel))
            for k in POOLS:
                row[f"true_{k}"] = float(truth[k][sel].sum())
            rows.append(row)
    return rows


def _f(x: float, w: int = 13) -> str:
    return f"{x:>{w},.0f}"


def report(rows: list[dict], axis: str) -> None:
    """The full table for one axis — one line per condition per arm, never collapsed."""
    sel = [r for r in rows if r["axis"] == axis]
    if not sel:
        return
    print()
    print("=" * 168)
    print(f"⭐⭐ RELAY A/B — axis = {axis}.  Counts in FRAGMENTS, against the origin-split oracle.")
    print("=" * 168)
    # ⭐ THE TRUTH IS THREE POOLS AND THE ESTIMATE IS TWO, and the table says so rather than hiding it:
    #   `nrna_true` / `mrna_true` are printed with no estimate beside them because calibration has none
    #   to give (AXIOM 0). They are here to say WHICH RNA a gDNA over-call is eating.
    print(
        f"   {'condition':<40} {'arm':<4} {'gdna_est':>13} {'gdna_true':>13} {'gdna_net':>13} "
        f"{'gdna_abs':>13} {'rna_est':>13} {'rna_true':>13} {'rna_net':>13} {'rna_abs':>13} "
        f"{'|nrna_true':>13} {'mrna_true':>13}"
    )
    by_cond: dict[str, list[dict]] = {}
    for r in sel:
        by_cond.setdefault(r["condition"], []).append(r)
    for cond, pair in by_cond.items():
        for r in sorted(pair, key=lambda r: r["arm"]):
            print(
                f"   {cond:<40} {r['arm']:<4} {_f(r['gdna_est'])} {_f(r['gdna_true'])} "
                f"{_f(r['gdna_net'])} {_f(r['gdna_abs'])} {_f(r['rna_est'])} {_f(r['rna_true'])} "
                f"{_f(r['rna_net'])} {_f(r['rna_abs'])} {_f(r['true_nrna'])} {_f(r['true_mrna'])}"
            )
        off = next(r for r in pair if r["arm"] == "off")
        # ⚠ the nascent stamp: an empty pool is not a measured pool
        nas = off["true_nrna"]
        stamp = (
            "   ⚠ nascent pool is EXACTLY 0 on this condition — NOT a measurement of it"
            if nas <= 0.0
            else f"   nascent truth {nas:,.0f} fragments"
        )
        for other in sorted((r for r in pair if r["arm"] != "off"), key=lambda r: r["arm"]):
            marks = []
            for p in ("gdna", "rna"):
                base = off[f"{p}_abs"]
                ratio = other[f"{p}_abs"] / base if base > 0 else float("nan")
                verdict = "HELPS" if ratio < 0.98 else ("HURTS" if ratio > 1.02 else "flat")
                marks.append(f"{p} {ratio:6.3f}x {verdict}")
            print(
                f"   {'':<40} {'->':<4} abs-error ratio {other['arm']}/off:  "
                + "   ".join(marks) + stamp
            )
    print()


def self_test() -> int:
    """⛔ Perturb every comparator, with NO I/O. A gate that has never been shown to fail is not a gate."""
    ok = fail = 0

    def check(name: str, cond: bool) -> None:
        nonlocal ok, fail
        if cond:
            ok += 1
        else:
            fail += 1
            print(f"   ⛔ {name}")

    n = 6
    truth = {
        "gdna": np.array([10.0, 0.0, 5.0, 0.0, 20.0, 0.0]),
        "nrna": np.array([0.0, 4.0, 0.0, 0.0, 0.0, 6.0]),
        "mrna": np.array([0.0, 6.0, 5.0, 8.0, 0.0, 4.0]),
    }
    allsel = np.ones(n, bool)
    mass = truth["gdna"] + truth["nrna"] + truth["mrna"]

    # ① a PERFECT f_g must give zero error on both pools and both error kinds
    perfect = np.where(mass > 0, truth["gdna"] / np.maximum(mass, _EPS), 0.0)
    s = score(perfect, truth, allsel)
    check("perfect f_g leaves gdna_abs > 0", abs(s["gdna_abs"]) < 1e-9)
    check("perfect f_g leaves rna_abs > 0", abs(s["rna_abs"]) < 1e-9)
    check("perfect f_g leaves gdna_net > 0", abs(s["gdna_net"]) < 1e-9)

    # ② the two errors are DIFFERENT statistics — a compensating pair cancels in net, not in abs
    comp = perfect.copy()
    comp[0] -= 0.5  # object 0 under-calls gDNA by 5
    comp[4] += 0.25  # object 4 over-calls gDNA by 5
    s2 = score(comp, truth, allsel)
    check("a compensating pair did not cancel in net", abs(s2["gdna_net"]) < 1e-9)
    check("a compensating pair also cancelled in abs — abs is not absolute", s2["gdna_abs"] > 9.99)

    # ③ the two pools are complementary: every fragment is gDNA or RNA, so the nets are opposite
    check("gdna_net and rna_net are not opposite", abs(s2["gdna_net"] + s2["rna_net"]) < 1e-9)

    # ④ truth totals must reproduce the partition exactly
    check("gdna_true wrong", abs(s["gdna_true"] - truth["gdna"].sum()) < 1e-9)
    check("rna_true wrong", abs(s["rna_true"] - (truth["nrna"] + truth["mrna"]).sum()) < 1e-9)

    # ⑤ a SELECTION must actually restrict — scoring a subset is not scoring everything
    half = np.array([True, True, True, False, False, False])
    check("selection ignored", score(perfect, truth, half)["mass"] < s["mass"])

    # ⑥ an all-half f_g on a pure-gDNA object is an error of exactly half its mass
    pure = {"gdna": np.array([100.0]), "nrna": np.array([0.0]), "mrna": np.array([0.0])}
    sh = score(np.array([0.5]), pure, np.ones(1, bool))
    check("half-call on pure gDNA is not 50", abs(sh["gdna_abs"] - 50.0) < 1e-9)
    check("half-call on pure gDNA has wrong sign", sh["gdna_net"] < 0)

    print(f"\n   self-test: {ok} passed, {fail} failed")
    return 1 if fail else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--suite", type=Path, default=DEFAULT_SUITE)
    ap.add_argument("--index", type=Path, default=DEFAULT_INDEX)
    ap.add_argument("--oracle-cache", type=Path, default=None)
    ap.add_argument("--conditions", nargs="*", default=None)
    ap.add_argument("--donor", type=Path, default=None,
                    help="a cached condition dir whose LIBRARY-LEVEL priors are injected — required "
                         "for a toy reference with no spliced reads of its own (docs/TESTING.md §0a)")
    ap.add_argument("--donor-index", type=Path, default=None,
                    help="the donor's index, when it differs from --index (a toy has its own)")
    ap.add_argument("--out", type=Path, default=None,
                    help="write every (condition x arm x axis) row as TSV — the machine-readable form "
                         "of the standing benchmark, and what `benchmark_report.py` renders")
    ap.add_argument("--arms", nargs="+", default=["off", "on"], choices=sorted(ARMS),
                    help="which arms to run — 'off'/'on' are the standing benchmark; 'currency' is "
                         "the Stage-3 policy under development (an 'off' baseline is required)")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()

    if args.self_test:
        return self_test()

    from rigel.calibration.region_arrays import RegionArrays

    index = TranscriptIndex.load(args.index)
    region_arrays = RegionArrays.from_index(index)
    oracle = args.oracle_cache or (args.suite / "oracle_cache")
    conds = args.conditions or sorted(p.name for p in (args.suite / "scan_cache").iterdir())
    injected = None
    if args.donor is not None:
        TH = _sibling("toy_harness.py")
        donor_index = TranscriptIndex.load(str(args.donor_index or args.index))
        injected = TH.harvest(args.donor, donor_index).priors
        print(f"⭐ library-level priors injected from the donor {args.donor.name} "
              f"(kappa={injected.rna_sense_frac:.6f}) — the toy cannot fit them for itself")
        # ⛔⛔ **THE DONOR MUST MATCH THE CONDITION'S AXES, AND A MISMATCH IS SILENT** — it produced a
        #    spurious catastrophe the first time this ran (2026-08-19): a `ss_0.99` donor injected into
        #    the `ss_0.50` conditions told the solver the library was stranded, and the zero control
        #    reported **82,581** false-positive gDNA fragments that a matching donor puts at **0**.
        #    Nothing in the output says "wrong donor", so it is refused here instead.
        def _axes(name: str) -> tuple[str, str]:
            return ("ss_0.99" if "ss_0.99" in name else "ss_0.50",
                    "capture_on" if "capture_on" in name else "capture_off")
        d_axes = _axes(args.donor.name)
        bad = [c for c in conds if _axes(c) != d_axes]
        if bad:
            raise SystemExit(
                f"⛔ the donor {args.donor.name} is {d_axes[0]} / {d_axes[1]}, but {len(bad)} of the "
                f"requested conditions are not — e.g. {bad[0]}. The injected priors are LIBRARY-level "
                "(kappa above all), so a mismatched donor silently answers a different library: run "
                "one donor per stratum, or pass --conditions for this donor's stratum only."
            )
        print()

    rows: list[dict] = []
    for c in conds:
        try:
            rows.extend(measure(index, region_arrays, args.suite, oracle, c, injected,
                                 arms=tuple(args.arms)))
            print(f"   ✔ {c}", flush=True)
        except Exception as exc:  # noqa: BLE001
            print(f"   ⛔ {c}: {type(exc).__name__}: {exc}", flush=True)
    for axis in ("REGION", "BOUNDARY", "ALL"):
        report(rows, axis)
    if args.out is not None and rows:
        cols = list(dict.fromkeys(k for r in rows for k in r))
        with open(args.out, "w") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in rows:
                fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")
        print(f"\n→ {args.out}  ({len(rows)} rows: condition x arm x axis)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
