#!/usr/bin/env python
"""⭐⭐⭐ **WHAT DOES ONE MESSAGE OPERATOR DO, PER SLOT, ON A REAL CHAIN?** — the fast half of the byte-identity gate.

Two policies, one real 70,176-slot chain, one process, **every output array compared element by element**.

⭐ **It is strictly stronger per condition than the panel, and strictly weaker across conditions.** The
panel reports 1,872 SCALARS, so an error that cancels between two slots is invisible there and visible
here; the panel in turn covers 36 conditions and this covers one. Run this first, on the condition with the
most structure, then run the panel. ⛔ It is also the honest answer to TRAPS: all-small-singly-large-jointly — *when every single
ablation is small and the joint one is large, go one stage upstream* — because an aggregate cannot tell a
switch that moves nothing from one that moves two slots in opposite directions.

⭐⭐ **THE VERDICT IT WAS BUILT FOR, and the reason it exists at all.** It gated the backbone restructure
against the solver it replaced: ``HeadPolicy`` reproduced the shipped answer on **421,056 output elements
and 18,245,830 diagnostic elements, zero differences**, and ``SilentPolicy`` reproduced muting psi's four
imputed channels on all 421,056 — which also PROVES the relay reaches the answer ONLY through those four
channels, since one arm runs the whole relay and discards it while the other never runs it at all. The
predecessor is deleted; the machinery is kept because pricing the switches one at a time needs exactly it.

What it compares
----------------
Two ``sweep.solve_chain`` arms on the SAME inputs:

* the six :class:`NodeBelief` arrays — ``f_pos``, ``f_neg``, ``f_g``, ``var_pos``, ``var_neg``, ``var_gdna``
* the two chain PROJECTIONS both feed ``CalibrationResult``
* every shared key of the diagnostics ``_capture``, because the dissect loop's instruments read it and a
  change that silently drops one of its keys breaks them with no error
* ⭐ the five BACKBONE ASSERTIONS, as counts beside their ELIGIBLE sets — because an assertion reporting
  zero violations where its predicate can never fire is not evidence of anything (TRAPS: could-the-arm-have-fired)

⛔ **TRAPS: byte-identity-gate's first lie is gated: this refuses to pass if it compared nothing.** It asserts both arms produced
the same key set and that the number of elements compared is nonzero, so "0 differences" cannot mean "0
comparisons".

Usage
-----
    # the live invariant: SilentPolicy is the message-free floor
    python scripts/design/backbone_parity.py --suite .../ladder --index .../rigel_index

    # what ONE operator does, per slot
    python scripts/design/backbone_parity.py --suite ... --index ... --arm-b no_peel

⚠ Both arms run in ONE process against ONE set of inputs, which is what makes this an identity test rather
than a reproducibility test.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from rigel.calibration import sweep as SW  # noqa: E402
from rigel.calibration.messages.head import HeadPolicy, HeadSwitches  # noqa: E402
from rigel.calibration.messages.silent import SilentPolicy  # noqa: E402

#: capture keys that CANNOT match by construction, with the reason. ⛔ Keep this list empty unless the
#: reason is structural — every entry is a hole in the gate.
_EXPECTED_ABSENT = {
    # the backbone publishes its assertion counts; the predecessor had no assertions to count.
    "backbone_assertions",
}


def _cmp(a, b):
    """Element-wise identity, not closeness. Returns ``(n_elements, n_differing, max_abs_delta)``."""
    if a is None and b is None:
        return 0, 0, 0.0
    if (a is None) != (b is None):
        return 1, 1, float("inf")
    if isinstance(a, dict) or isinstance(b, dict):
        return 0, 0, 0.0  # compared key-by-key by the caller where it matters
    try:
        x = np.asarray(a, np.float64)
        y = np.asarray(b, np.float64)
    except (TypeError, ValueError):
        return 1, int(a != b), 0.0
    if x.shape != y.shape:
        return max(x.size, y.size), max(x.size, y.size), float("inf")
    same = (x == y) | (np.isnan(x) & np.isnan(y))
    n_diff = int((~same).sum())
    d = 0.0
    if n_diff:
        with np.errstate(invalid="ignore"):
            dd = np.abs(y - x)[~same]
            dd = dd[np.isfinite(dd)]
            d = float(dd.max()) if dd.size else float("inf")
    return int(x.size), n_diff, d


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--suite", type=Path, required=True)
    ap.add_argument("--index", type=Path, required=True)
    ap.add_argument("--condition", default="gdna_g50_ss_0.50_nrna_none_capture_on")
    ap.add_argument("--oracle-cache", type=Path, default=None)
    ap.add_argument("--work-dir", type=Path, default=Path("/tmp/rigel_backbone_parity"))
    ap.add_argument("--arm-a", default="head", help="'head', 'silent', or 'no_<switch>'")
    ap.add_argument("--arm-b", default="silent", help="'head', 'silent', or 'no_<switch>'")
    args = ap.parse_args()

    from rigel.config import CalibrationConfig, PipelineConfig  # noqa: PLC0415
    from rigel.index import TranscriptIndex  # noqa: PLC0415

    # ── capture ONE real set of node_sweep inputs by intercepting the production call ──────────────────
    grabbed: list[dict] = []
    orig = SW.solve_chain

    def spy(chain, statics, geometry, belief, region_arrays, **kw):
        if len(grabbed) < 2:
            grabbed.append(
                {
                    "chain": chain,
                    "statics": statics,
                    "geometry": geometry,
                    "belief": belief,
                    "region_arrays": region_arrays,
                    "kw": dict(kw),
                }
            )
        return orig(chain, statics, geometry, belief, region_arrays, **kw)

    # ⛔⛔ TRAPS: an-ablation-that-never-ran, and it fired on the first run: ``import rigel.calibration.calibrate as CAL`` binds the
    # re-exported FUNCTION, not the module, so patching an attribute on it patches nothing and the spy
    # reads as "never called". Go through ``sys.modules``, and RAISE if the patch did not fire.
    CAL = sys.modules["rigel.calibration.calibrate"]
    CAL.solve_chain = spy
    cond = args.suite / args.condition
    print(f"\n   running one calibration to capture real inputs: {args.condition}", flush=True)
    index = TranscriptIndex.load(str(args.index))
    from scripts.design import pass0_vs_oracle as P0  # noqa: PLC0415

    P0.measure_condition(
        bam=str(cond / "sim_oracle.bam"),
        index=index,
        pipeline_config=PipelineConfig(),
        calibration_config=CalibrationConfig(),
        work_dir=args.work_dir,
        tag=args.condition,
        truth_pmfs=None,
        oracle_cache=args.oracle_cache,
    )
    CAL.solve_chain = orig
    if not grabbed:
        raise SystemExit("⛔ never captured a solve_chain call — the spy did not fire (TRAPS.md an-ablation-that-never-ran)")
    g = grabbed[0]
    kw = {k: v for k, v in g["kw"].items() if k not in ("policy", "_capture")}
    n = int(g["chain"].n_slots)
    print(f"   captured a chain of {n:,} slots (prior-free pass)", flush=True)

    def run(fn, **extra):
        cap: dict = {}
        out = fn(
            g["chain"], g["statics"], g["geometry"], g["belief"], g["region_arrays"],
            _capture=cap, **kw, **extra,
        )
        return out, cap

    def policy_for(spec: str):
        """``head`` | ``silent`` | ``no_<switch>``. ⛔ An unknown switch name RAISES rather than being
        silently ignored — an arm that turns nothing off scores identical and reads as "inert" (TRAPS: an-ablation-that-never-ran)."""
        if spec == "silent":
            return SilentPolicy()
        if spec == "head":
            return HeadPolicy()
        if spec.startswith("no_"):
            name = spec[3:]
            if name not in HeadSwitches().names():
                raise SystemExit(
                    f"⛔ no such switch {name!r}. Available: {', '.join(HeadSwitches().names())}"
                )
            return HeadPolicy(HeadSwitches(**{name: False}))
        raise SystemExit(f"⛔ unknown arm {spec!r} — use 'head', 'silent' or 'no_<switch>'")

    pa, pb = policy_for(args.arm_a), policy_for(args.arm_b)
    if pa.name == pb.name:
        raise SystemExit(f"⛔ both arms are {pa.name!r} — this would compare a run against itself (TRAPS: could-the-arm-have-fired)")
    print(f"\n   ARM A: {pa.name}", flush=True)
    old, cap_old = run(SW.solve_chain, policy=pa)
    print(f"   ARM B: {pb.name}", flush=True)
    new, cap_new = run(SW.solve_chain, policy=pb)

    rows, tot_el, tot_diff = [], 0, 0
    for f in ("f_pos", "f_neg", "f_g", "var_pos", "var_neg", "var_gdna"):
        el, nd, d = _cmp(getattr(old, f), getattr(new, f))
        rows.append((f"belief.{f}", el, nd, d))
        tot_el += el
        tot_diff += nd

    sub = g["region_arrays"]
    for nm, fn in (("chain_node_deconv", SW.chain_node_deconv), ("chain_edge_deconv", SW.chain_edge_deconv)):
        try:
            a = fn(g["chain"], old, sub)
            b = fn(g["chain"], new, sub)
        except AttributeError:
            continue  # region_arrays is not the substrate object; the belief comparison covers it
        for f in ("gdna_mass", "rna_mass", "gdna_frac", "rna_pos_frac", "rna_neg_frac"):
            el, nd, d = _cmp(getattr(a, f), getattr(b, f))
            rows.append((f"{nm}.{f}", el, nd, d))
            tot_el += el
            tot_diff += nd

    # ── the diagnostics capture: the dissect loop reads it, so a dropped key is a real regression ──────
    ka, kb = set(cap_old), set(cap_new)
    missing = sorted(ka - kb)
    added = sorted((kb - ka) - _EXPECTED_ABSENT)
    cap_el = cap_diff = 0
    cap_bad = []
    for k in sorted(ka & kb):
        a, b = cap_old[k], cap_new[k]
        if isinstance(a, list) and isinstance(b, list):
            if len(a) != len(b):
                cap_bad.append((f"_capture[{k}] (list len {len(a)} vs {len(b)})", 1, 1, float("inf")))
                continue
            for j, (xa, xb) in enumerate(zip(a, b)):
                for kk in sorted(set(xa) & set(xb)) if isinstance(xa, dict) else []:
                    el, nd, d = _cmp(xa[kk], xb[kk])
                    cap_el += el
                    cap_diff += nd
                    if nd:
                        cap_bad.append((f"_capture[{k}][{j}][{kk}]", el, nd, d))
            continue
        if isinstance(a, dict) and isinstance(b, dict):
            for kk in sorted(set(a) & set(b)):
                el, nd, d = _cmp(a[kk], b[kk])
                cap_el += el
                cap_diff += nd
                if nd:
                    cap_bad.append((f"_capture[{k}][{kk}]", el, nd, d))
            for kk in sorted(set(a) - set(b)):
                cap_bad.append((f"_capture[{k}][{kk}] MISSING", 1, 1, float("inf")))
            continue
        el, nd, d = _cmp(a, b)
        cap_el += el
        cap_diff += nd
        if nd:
            cap_bad.append((f"_capture[{k}]", el, nd, d))

    print()
    print(f"   {'array':<34} {'elements':>12} {'differing':>10} {'max |delta|':>14}")
    print("   " + "-" * 74)
    for nm, el, nd, d in rows:
        flag = "  ⛔" if nd else ""
        print(f"   {nm:<34} {el:>12,} {nd:>10,} {d:>14.6g}{flag}")
    print("   " + "-" * 74)
    print(f"   {'OUTPUT TOTAL':<34} {tot_el:>12,} {tot_diff:>10,}")
    print(f"   {'_capture TOTAL':<34} {cap_el:>12,} {cap_diff:>10,}")

    if cap_bad:
        print("\n   ⛔ capture differences (the dissect loop reads these):")
        for nm, el, nd, d in cap_bad[:40]:
            print(f"      {nm:<52} {nd:>9,}/{el:<10,} max {d:.6g}")
        if len(cap_bad) > 40:
            print(f"      … and {len(cap_bad) - 40} more")
    if missing:
        print(f"\n   ⚠ diagnostic keys only ARM A publishes: {missing}")
        print("      ⭐ EXPECTED for an ablation (a silent or switched-off operator publishes nothing);")
        print("         a BUG for an identity gate, where both arms must publish the same keys.")
    if added:
        print(f"\n   ⚠ diagnostic keys only ARM B publishes: {added}")

    # ── THE FIVE BACKBONE ASSERTIONS, as counts — and TRAPS: could-the-arm-have-fired is why the ELIGIBLE column is here ───────────
    # An assertion reporting 0 violations on a substrate where its predicate can never fire is not evidence
    # of anything. So print what each one could have caught beside what it did.
    # ⛔ BOTH arms, because a policy that sends nothing skips every check on a message and its table would
    # otherwise read as "the assertion holds" when the truth is "the assertion never ran" — TRAPS: could-the-arm-have-fired one level up.
    aa = cap_old.get("backbone_assertions") or {}
    ab = cap_new.get("backbone_assertions") or {}
    if aa or ab:
        print()
        print("   ⭐ THE FIVE BACKBONE ASSERTIONS — violations / eligible, BOTH arms")
        print(f"   {'assertion':<32} {pa.name[:16]:>17} {pb.name[:16]:>17}   verdict")
        print("   " + "-" * 92)
        for k in sorted(set(aa) | set(ab)):
            cells, worst = [], None
            for d in (aa, ab):
                if k not in d:
                    cells.append("        not run")
                    continue
                v, e = d[k]["violations"], d[k]["eligible"]
                cells.append(f"{v:>8,}/{e:<8,}")
                if e and v:
                    worst = max(worst or 0.0, 100.0 * v / e)
                elif e == 0 and worst is None:
                    worst = -1.0
            if worst is None:
                verdict = "✅ holds"
            elif worst < 0:
                verdict = "⚠ NO ELIGIBLE SLOTS — said nothing here (TRAPS: could-the-arm-have-fired)"
            else:
                waived = k in SW._KNOWN_VIOLATIONS
                verdict = f"{'⛔ WAIVED' if waived else '⛔ UNWAIVED'}, up to {worst:.2f}% of eligible"
            print(f"   {k:<32} {cells[0]:>17} {cells[1]:>17}   {verdict}")

    ok = tot_diff == 0 and cap_diff == 0 and not missing and not added and not cap_bad
    if tot_el == 0:
        print("\n   ⛔ COMPARED NOTHING — this gate would have passed vacuously (TRAPS.md byte-identity-gate)")
        return 1
    print()
    if ok:
        print(f"   ✅ {pa.name} and {pb.name} are BYTE-IDENTICAL: {tot_el:,} output elements and "
              f"{cap_el:,} capture elements, zero differences")
        print("      ⚠ For an ABLATION arm that is not good news — it means the operator is INERT on this")
        print("        condition, and TRAPS: hard-labels-miss-soft-change says treat an identical result as NO EVIDENCE, not as no change.")
    else:
        print(f"   ⭐ {pa.name} vs {pb.name}: {tot_diff:,} of {tot_el:,} output elements differ, "
              f"{cap_diff:,} of {cap_el:,} capture elements")
        print("      ⛔ For an IDENTITY gate any difference is a bug. For an ablation arm this is the")
        print("         measurement — now score it PER STRATUM on the panel (TRAPS: panel-before-src), never pooled.")
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
