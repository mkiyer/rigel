#!/usr/bin/env python
"""What does one message operator do, per slot, on a real chain? Two policies run through
``sweep.solve_chain`` on the same captured inputs, in one process, and every output array is
compared element by element: the six ``RegionBelief`` arrays, the two chain projections that feed
``CalibrationResult``, every shared key of the diagnostics ``_capture`` (the dissect loop reads it,
so a silently dropped key is a regression), and the backbone assertions as violation counts beside
their eligible sets, because zero violations where the predicate can never fire is not evidence.
It is strictly stronger than the panel per condition — an aggregate is a handful of scalars, so an
error that cancels between two slots is invisible there and visible here — and covers one condition
where the panel covers the ladder, so run it first and the panel after. ``transfer`` against
``silent`` is the per-slot view of the shipped policy; ``transfer`` against a prototype arm
(``module:<file.py>:<arm>``, the ``ARMS`` table `policy_prototype.py` loads) is the per-slot view
of one mechanism, and a prototype identical to the shipped policy must score byte-identical. The
policy is the arm, named on both sides and printed before each run, so there is no separate
``--messages`` switch; the one calibration at the top exists only to capture ``solve_chain``'s
inputs and the shipped policy instance. Identity is bit-equality, never a tolerance, and it refuses
to pass if it compared nothing or if both arms are the same policy. For an ablation arm a
byte-identical result means the operator is inert here and is no evidence, not "no change".

Usage::

    python scripts/design/backbone_parity.py --suite .../ladder --index .../rigel_index   # shipped transfer vs silent, per slot
    python scripts/design/backbone_parity.py --suite ... --index ... --arm-b module:proto.py:my_arm   # one prototype mechanism
    python scripts/design/backbone_parity.py --suite ... --index ... --condition <name> --oracle-cache <dir>
"""

from __future__ import annotations

import argparse
import dataclasses
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from rigel.calibration import sweep as SW  # noqa: E402
from rigel.calibration.messages.silent import SilentPolicy  # noqa: E402

#: capture keys that cannot match by construction. Keep this set empty unless the reason is
#: structural — every entry is a hole in the gate.
_EXPECTED_ABSENT: set[str] = set()


def _cmp(a, b, where: str = ""):
    """Element-wise identity, not closeness. Returns ``(n_elements, n_differing, max_abs_delta)``.

    A type this function does not understand raises, naming the key, rather than falling through to
    ``!=``: a dataclass holding arrays would answer that through a generated ``__eq__`` that cannot
    reduce an array to one bool, and nothing in the suite reads ``_capture``, so a silent guess here
    leaves the instrument dead with a green suite.
    """
    if a is None and b is None:
        return 0, 0, 0.0
    if (a is None) != (b is None):
        return 1, 1, float("inf")
    # a dict is walked key by key here as well as by the caller: ``main`` walks the top level of
    #   ``_capture`` itself, but a dict reached through a list element or a dataclass field has no
    #   such caller, and returning "0 elements, 0 differing" for it would hide data that differed.
    if isinstance(a, dict) or isinstance(b, dict):
        if not (isinstance(a, dict) and isinstance(b, dict)):
            return 1, 1, float("inf")
        el = nd = 0
        d = 0.0
        for kk in sorted(set(a) & set(b)):
            e, n, dd = _cmp(a[kk], b[kk], f"{where}[{kk}]")
            el += e
            nd += n
            d = max(d, dd)
        one_sided = set(a) ^ set(b)
        if one_sided:
            el += len(one_sided)
            nd += len(one_sided)
            d = float("inf")
        return el, nd, d
    # a dataclass is walked field by field, because its members carry the numbers; comparing the
    #   objects would ask a generated ``__eq__`` to reduce an array to one bool, and a field-wise walk
    #   keeps the element count honest (all-``None`` members contribute 0 comparisons, not 1).
    if dataclasses.is_dataclass(a) or dataclasses.is_dataclass(b):
        if type(a) is not type(b):
            return 1, 1, float("inf")
        el = nd = 0
        d = 0.0
        for f in dataclasses.fields(a):
            e, n, dd = _cmp(getattr(a, f.name), getattr(b, f.name), f"{where}.{f.name}")
            el += e
            nd += n
            d = max(d, dd)
        return el, nd, d
    if isinstance(a, (str, bytes)) or isinstance(b, (str, bytes)):
        return 1, int(a != b), 0.0
    try:
        x = np.asarray(a, np.float64)
        y = np.asarray(b, np.float64)
    except (TypeError, ValueError) as exc:
        raise TypeError(
            f"⛔ _cmp cannot compare {where or '<unnamed>'}: {type(a).__name__} vs {type(b).__name__}. "
            "Add a branch for it — falling through to `!=` is what killed this instrument once."
        ) from exc
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
    ap.add_argument("--condition", default="gdna_g50_ss_0.50_nrna_mid_capture_on")
    ap.add_argument("--oracle-cache", type=Path, default=None)
    ap.add_argument("--work-dir", type=Path, default=Path("/tmp/rigel_backbone_parity"))
    ap.add_argument("--arm-a", default="transfer", help="'transfer', 'silent', or 'module:<file.py>:<arm>'")
    ap.add_argument("--arm-b", default="silent", help="'transfer', 'silent', or 'module:<file.py>:<arm>'")
    args = ap.parse_args()

    from rigel.config import CalibrationConfig, PipelineConfig  # noqa: PLC0415
    from rigel.index import TranscriptIndex  # noqa: PLC0415

    # ── capture ONE real set of region_sweep inputs by intercepting the production call ──────────────────
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

    # ``import rigel.calibration.calibrate as CAL`` binds the re-exported FUNCTION, not the module, so
    # patching an attribute on it patches nothing and the spy reads as "never called". Go through
    # ``sys.modules``, and RAISE if the patch did not fire (`TRAPS: an-ablation-that-never-ran`).
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

    # the captured policy is the shipped transfer policy with its intron-row memo and strand model;
    # a prototype arm is built on the same two, so the two arms differ by the mechanism alone.
    _shipped = g["kw"].get("policy")

    def policy_for(spec: str):
        """``transfer`` | ``silent`` | ``module:<file.py>:<arm>``. An unknown arm raises rather than
        being silently ignored — an arm that changes nothing scores identical and reads as inert."""
        if spec == "silent":
            return SilentPolicy()
        if spec == "transfer":
            if _shipped is None or getattr(_shipped, "name", None) == "silent":
                raise SystemExit(
                    "⛔ the captured shipped policy is not the transfer policy — run under the shipped config"
                )
            return _shipped
        if spec.startswith("module:"):
            _, path, arm = spec.split(":", 2)
            from scripts.design import policy_prototype as PP  # noqa: PLC0415

            arms = PP.load_arms(Path(path))
            if arm not in arms:
                raise SystemExit(f"⛔ {path} defines no arm {arm!r}; it has {sorted(arms)}")
            return arms[arm](strand=_shipped._strand)
        raise SystemExit(f"⛔ unknown arm {spec!r} — use 'transfer', 'silent' or 'module:<file.py>:<arm>'")

    pa, pb = policy_for(args.arm_a), policy_for(args.arm_b)
    if pa.name == pb.name:
        raise SystemExit(f"⛔ both arms are {pa.name!r} — this would compare a run against itself (TRAPS: could-the-arm-have-fired)")
    print(f"\n   ARM A: {pa.name}", flush=True)
    old, cap_old = run(SW.solve_chain, policy=pa)
    print(f"   ARM B: {pb.name}", flush=True)
    new, cap_new = run(SW.solve_chain, policy=pb)

    rows, tot_el, tot_diff = [], 0, 0
    for f in ("f_pos", "f_neg", "f_g", "var_gdna"):
        el, nd, d = _cmp(getattr(old, f), getattr(new, f), f"belief.{f}")
        rows.append((f"belief.{f}", el, nd, d))
        tot_el += el
        tot_diff += nd

    sub = g["region_arrays"]
    for nm, fn in (("chain_region_deconv", SW.chain_region_deconv), ("chain_boundary_deconv", SW.chain_boundary_deconv)):
        try:
            a = fn(g["chain"], old, sub)
            b = fn(g["chain"], new, sub)
        except AttributeError:
            continue  # region_arrays is not the substrate object; the belief comparison covers it
        for f in ("gdna_mass", "rna_mass", "gdna_frac", "rna_pos_frac", "rna_neg_frac"):
            el, nd, d = _cmp(getattr(a, f), getattr(b, f), f"{nm}.{f}")
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
                    el, nd, d = _cmp(xa[kk], xb[kk], f"_capture[{k}][{j}][{kk}]")
                    cap_el += el
                    cap_diff += nd
                    if nd:
                        cap_bad.append((f"_capture[{k}][{j}][{kk}]", el, nd, d))
            continue
        if isinstance(a, dict) and isinstance(b, dict):
            for kk in sorted(set(a) & set(b)):
                el, nd, d = _cmp(a[kk], b[kk], f"_capture[{k}][{kk}]")
                cap_el += el
                cap_diff += nd
                if nd:
                    cap_bad.append((f"_capture[{k}][{kk}]", el, nd, d))
            for kk in sorted(set(a) - set(b)):
                cap_bad.append((f"_capture[{k}][{kk}] MISSING", 1, 1, float("inf")))
            continue
        el, nd, d = _cmp(a, b, f"_capture[{k}]")
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
        # VALUE differences are printed BEFORE structural ones: a policy ablation produces dozens of
        #   "key MISSING" rows by construction and the list is truncated at 40, so a real numeric
        #   difference sorted below them would appear only as a count in the TOTAL line.
        n_val = sum(1 for _n, _e, _d, dd in cap_bad if np.isfinite(dd))
        cap_bad.sort(key=lambda r: (0 if np.isfinite(r[3]) else 1, -r[2]))
        # "structural" is every row whose delta is not finite: a missing key, a shape mismatch, a None
        #   on one side only, or a type mismatch.
        print(f"\n   ⛔ capture differences (the dissect loop reads these): {n_val:,} carry a VALUE "
              f"delta, {len(cap_bad) - n_val:,} are structural (a missing key, shape, None or type)")
        for nm, el, nd, d in cap_bad[:40]:
            print(f"      {nm:<52} {nd:>9,}/{el:<10,} max {d:.6g}")
        if len(cap_bad) > 40:
            print(f"      … and {len(cap_bad) - 40} more (VALUE deltas are listed first)")
    if missing:
        print(f"\n   ⚠ diagnostic keys only ARM A publishes: {missing}")
        print("      ⭐ EXPECTED for an ablation (a silent or switched-off operator publishes nothing);")
        print("         a BUG for an identity gate, where both arms must publish the same keys.")
    if added:
        print(f"\n   ⚠ diagnostic keys only ARM B publishes: {added}")

    # ── the backbone assertions, as violation counts beside their ELIGIBLE sets ──────────────────────
    # An assertion reporting 0 violations where its predicate can never fire is not evidence of
    # anything, so print what each one could have caught beside what it did — for BOTH arms, because a
    # policy that sends nothing skips every check on a message and would otherwise read as "holds"
    # when the truth is "never ran" (`TRAPS: could-the-arm-have-fired`).
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
