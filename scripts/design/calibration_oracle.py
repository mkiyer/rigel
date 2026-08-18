"""THE CALIBRATION ORACLE — per-REGION and per-BOUNDARY ground-truth composition, CERTIFIED.

⭐⭐⭐ **THE ONE TABLE CALIBRATION IS DEBUGGED AGAINST, AND EVERY ROW OF IT IS GATED BEFORE IT MAY BE
USED.** Owner, 2026-08-18: calibration is showing catastrophic error, so nothing downstream of the
accumulator can be trusted until the ground truth itself is proven. This instrument produces, for every
slot of the chain:

    count        the ACCUMULATOR's own tally (the mixture calibration deconvolves)
    n_gdna       TRUE genomic-DNA fragments        } realized, from the origin-split oracle BAMs,
    n_nrna       TRUE unannotated nascent fragments } each scanned by the SAME accumulator and
    n_mrna       TRUE annotated-RNA fragments       } re-validated against the full scan
    true_f_g     n_gdna / (n_gdna + n_nrna + n_mrna)

and REFUSES to emit it unless every gate below passes. ⛔ An oracle that is merely *plausible* is how a
calibration bug gets debugged against a truth bug and both survive.

**THE GATES — what "100 % sure" means here, stated as the checks that enforce it. Named, never
numbered** (`tests/test_no_jargon_labels.py`):

==============================  ======================================================================
**sum-to-full**                 the three origin partitions re-scanned through the accumulator
                                reconstruct the full scan's every bank (`OracleTruth.from_parts`); a
                                cached partition that does not describe this index/scan is refused by
                                the cache layer itself
**partition-projects-exactly**  at every slot ``n_gdna + n_nrna + n_mrna == count`` — the projection
                                used for truth is the projection used for the estimate, or the two
                                tables describe different populations
**gdna-field-uniformity**       ⭐⭐ the gDNA field is UNIFORM where physics says it must be
                                (capture-OFF, per genomic reference): every slot's expected gDNA count
                                is ``rho_ref x E_g(slot)``. Scored as Poisson z per slot and as a
                                ratio-z per object CLASS — a per-class bias is the signature of
                                exactly the counting/arithmetic bug this campaign hunts
**exact-zeros**                 a reference with no simulated gDNA, or a ``g00`` condition, must carry
                                a gdna partition that is identically 0 — counted VACUOUS, never
                                "pass", when the whole bank is empty
**nascent-in-annotation**       ``n_nrna > 0`` only where a transcript span admits RNA
                                (``free_pos | free_neg``); nascent in intergenic space is a split or
                                projection bug
==============================  ======================================================================

⚠ **WHAT IS *NOT* YET CERTIFIED, SAID PLAINLY RATHER THAN IMPLIED:** the per-object ANALYTIC RNA
expectation (abundance x per-transcript opportunity). The RNA truth here is *realized* — certified as
accumulator-consistent (sum-to-full, partition-projects-exactly) and annotation-consistent
(nascent-in-annotation), and it inherits fragment-level truth
from the simulator's read names via the origin split. The analytic cross-check is the natural next
gate; the gDNA arm (**gdna-field-uniformity**) is fully analytic already because a uniform field needs no per-transcript
model.

Output: a per-slot ``.npz`` beside the report, for `calibration_walk.py` and any stepwise dissection
to read — one truth source, derived once, gated once.

Usage::

    python scripts/design/calibration_oracle.py --condition NAME     # certify one condition
    python scripts/design/calibration_oracle.py                      # certify the whole ladder
    python scripts/design/calibration_oracle.py --self-test          # perturb every gate, no I/O
"""

from __future__ import annotations

import argparse
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

from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.region_chain import REGION, build_region_chain  # noqa: E402
from rigel.calibration.region_geometry import build_region_geometry, build_region_statics  # noqa: E402
from rigel.calibration.splice_graph import build_boundary_flags_array, build_sj_geometry_arrays  # noqa: E402
from rigel.calibration.substrate import CalibrationSubstrate  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402
from rigel.scan_cache import calibration_inputs, read_scan_cache  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tests"))
from calibration._oracle import ORIGINS, OracleTruth  # noqa: E402

DEFAULT_SUITE = OC.DEFAULT_SUITE
DEFAULT_INDEX = OC.DEFAULT_INDEX
_EPS = 1.0e-12
#: |z| above which a slot is flagged under **gdna-field-uniformity**. 4 sigma two-sided is ~6e-5 expected false flags per
#: 70k slots ~ 4 slots; the gate is on the FLAG RATE and the CLASS ratio, not on any single slot.
_Z_FLAG = 4.0
#: Poisson slots below this expectation are pooled into their class rather than z-scored singly —
#: the normal approximation is not honest there, and a per-slot z on lambda = 0.3 flags noise.
_MIN_LAMBDA = 5.0


# ── the gates, arrays in / verdict out, so the self-test can feed synthetic inputs ────────────────


def gate_partition(count: np.ndarray, n_g, n_n, n_m, atol: float = 1e-6) -> dict:
    """**partition-projects-exactly** — the origin partition must reproduce the accumulator count at EVERY slot."""
    gap = np.abs((n_g + n_n + n_m) - count)
    worst = float(gap.max()) if gap.size else 0.0
    return {"gate": "gate partition-projects-exactly", "ok": bool(worst <= atol), "worst_gap": worst,
            "n_bad": int((gap > atol).sum())}


def gate_uniformity(n_g: np.ndarray, eff_g: np.ndarray, classes: np.ndarray,
                    select: np.ndarray) -> dict:
    """**gdna-field-uniformity** — Poisson uniformity of the realized gDNA field against ``rho x E_g``, per class.

    ⭐ The CLASS ratio-z is the counting-bug detector: a deposit-rule or divisor bug is per-KIND
    (regions vs boundaries vs a signature class), so it shows as one class sitting off the shared
    rate while per-slot flags stay unremarkable. ⛔ VACUOUS when the selected field is empty.
    """
    sel = select & (eff_g > 0.0)
    tot_n, tot_e = float(n_g[sel].sum()), float(eff_g[sel].sum())
    if tot_n <= 0.0:
        return {"gate": "gate gdna-field-uniformity", "ok": True, "vacuous": True}
    rho = tot_n / tot_e
    lam = rho * eff_g
    zable = sel & (lam >= _MIN_LAMBDA)
    z = (n_g[zable] - lam[zable]) / np.sqrt(lam[zable])
    flag_rate = float(np.mean(np.abs(z) > _Z_FLAG)) if z.size else 0.0
    rows = []
    ok = True
    for cls in sorted(set(classes[sel].tolist())):
        m = sel & (classes == cls)
        cn, ce = float(n_g[m].sum()), float(rho * eff_g[m].sum())
        cz = (cn - ce) / np.sqrt(max(ce, _EPS))
        rows.append({"class": cls, "n": cn, "expected": ce, "ratio": cn / max(ce, _EPS), "z": cz})
        # ⚠ class totals are 1e5-1e6 fragments, so a real bias is hundreds of sigma; 6 allows the
        #   slight non-independence of the shared rho-hat without admitting any real bias.
        if abs(cz) > 6.0:
            ok = False
    return {"gate": "gate gdna-field-uniformity", "ok": ok and flag_rate < 1e-3, "vacuous": False,
            "rho": rho, "slot_flag_rate": flag_rate, "n_z_scored": int(z.size), "classes": rows}


def gate_zero(n: np.ndarray, select: np.ndarray, what: str) -> dict:
    """**exact-zeros** — a bank that must be empty, checked as EXACTLY zero. Vacuous if nothing is selected."""
    if not select.any():
        return {"gate": f"gate exact-zeros:{what}", "ok": True, "vacuous": True}
    bad = float(n[select].sum())
    return {"gate": f"gate exact-zeros:{what}", "ok": bool(bad == 0.0), "vacuous": False, "mass": bad}


def gate_nascent_scope(n_n: np.ndarray, rna_admissible: np.ndarray) -> dict:
    """**nascent-in-annotation** — nascent RNA may exist only where the annotation admits ANY RNA."""
    bad = float(n_n[~rna_admissible].sum())
    return {"gate": "gate nascent-in-annotation", "ok": bool(bad == 0.0), "out_of_scope_mass": bad}


# ── the derivation ────────────────────────────────────────────────────────────────────────────────


def derive(index, region_arrays, suite: Path, condition: str) -> tuple[dict, list[dict]]:
    """One condition's certified per-slot truth table plus its gate verdicts."""
    sj = build_sj_geometry_arrays(index)
    bflags = build_boundary_flags_array(index)
    cache = read_scan_cache(Path(suite) / "scan_cache" / condition, index)
    kw = calibration_inputs(cache, index)
    chain = build_region_chain(cache.payload.ref_region_offsets, cache.payload.ref_boundary_offsets)
    statics = build_region_statics(chain, region_arrays, bflags)
    geom = build_region_geometry(
        chain, CalibrationSubstrate.from_payload(cache.payload, region_arrays),
        region_arrays, sj, kw["gdna_fl_pmf"], kw["rna_fl_pmf"], None,
    )
    label = OC.strata(chain, statics, geom, region_arrays)["label"]

    # **sum-to-full** — sum-to-full runs as a hard gate inside from_parts; an exception IS the verdict.
    root = Path(suite) / "oracle_cache" / condition
    parts = {k: read_scan_cache(root / k, index).payload for k in ORIGINS}
    OracleTruth.from_parts(cache.payload, parts)

    count = OC.slot_counts(cache.payload, region_arrays, chain)
    n_g = OC.slot_counts(parts["gdna"], region_arrays, chain)
    n_n = OC.slot_counts(parts["nrna"], region_arrays, chain)
    n_m = OC.slot_counts(parts["mrna"], region_arrays, chain)

    kind = np.asarray(chain.kind)
    obj = np.asarray(chain.obj_idx, np.int64)
    is_region = kind == REGION
    left = np.clip(np.asarray(chain.left, np.int64), 0, int(chain.n_slots) - 1)
    # a slot's reference: a REGION's own; a BOUNDARY sits between two regions of ONE reference,
    # so its left flank's region names it.
    ref = np.where(is_region, region_arrays.ref_id[np.clip(obj, 0, len(region_arrays.ref_id) - 1)], -1)
    ref = np.where(is_region, ref, ref[left])
    eff_g = np.asarray(geom.eff_gdna, np.float64)
    rna_ok = np.asarray(statics.free_pos, bool) | np.asarray(statics.free_neg, bool)

    capture_on = condition.endswith("_capture_on")
    gdna_refs = np.array(sorted({int(r) for r in ref[n_g > 0]}), np.int64)

    verdicts = [gate_partition(count, n_g, n_n, n_m)]
    if capture_on:
        verdicts.append({"gate": "gate gdna-field-uniformity", "ok": True, "vacuous": True,
                         "note": "capture-ON: the field is deliberately non-uniform; gate not applicable"})
    else:
        for r in gdna_refs:
            v = gate_uniformity(n_g, eff_g, label, ref == r)
            v["ref"] = int(r)
            verdicts.append(v)
        if gdna_refs.size == 0:
            verdicts.append({"gate": "gate gdna-field-uniformity", "ok": True, "vacuous": True,
                             "note": "no gDNA anywhere (a zero condition)"})
    # refs that carry no gDNA at all must be EXACTLY zero (ERCC backbone, and every ref at g00)
    verdicts.append(gate_zero(n_g, ~np.isin(ref, gdna_refs), "gdna outside genomic refs"))
    verdicts.append(gate_nascent_scope(n_n, rna_ok))

    mass = n_g + n_n + n_m
    table = {
        "condition": condition, "kind": kind, "obj": obj, "stratum": label.astype(str),
        "ref": ref, "eff_g": eff_g, "eff_r": np.asarray(geom.eff_rna, np.float64),
        "count": count, "n_gdna": n_g, "n_nrna": n_n, "n_mrna": n_m,
        "true_f_g": np.where(mass > 0, n_g / np.maximum(mass, _EPS), 0.0),
        "live": mass > 0,
    }
    return table, verdicts


def report(condition: str, verdicts: list[dict]) -> tuple[bool, bool]:
    """Two certification LEVELS, because the gates certify two different things.

    ⭐ **COMPOSITION-certified** (**sum-to-full** **partition-projects-exactly** **exact-zeros** **nascent-in-annotation**): ``true_f_g`` is sound — it is realized counts against
    realized counts through one accumulator, and no opportunity model enters it anywhere.
    ⭐ **FIELD-certified** (**gdna-field-uniformity** as well): the deposit geometry also matches the opportunity model, so a
    DENSITY (``n/E``) may be trusted too. ⛔ First run of this gate on real data (2026-08-18) FAILED
    field certification: every BOUNDARY class reads **2-3 % below** ``rho x E_cross`` while every
    REGION class is exact — reproduced on two refs and three independent simulations. ``true_f_g`` is
    untouched by that; every density calibration divides out is not, so the table is stamped with
    both verdicts rather than one.
    """
    comp_ok, field_ok = True, True
    print(f"\n== {condition}")
    for v in verdicts:
        ok = v["ok"]
        if v["gate"].startswith("gate gdna-field"):
            field_ok &= ok
        else:
            comp_ok &= ok
        mark = "✔" if ok else "⛔"
        extra = ""
        if v.get("vacuous"):
            extra = "  (VACUOUS — the selected field is empty; this is not a pass)"
        elif v["gate"].startswith("gate gdna-field") and not v.get("vacuous"):
            extra = (f"  rho={v['rho']:.6g}  slot flags {v['slot_flag_rate']:.2e} "
                     f"over {v['n_z_scored']:,} z-scored slots")
        ref_tag = f"  ref {v['ref']}" if "ref" in v else ""
        print(f"   {mark} {v['gate']:<34}{ref_tag}{extra}")
        if v["gate"].startswith("gate gdna-field") and not v.get("vacuous"):
            for c in v["classes"]:
                flag = "" if abs(c["z"]) <= 6.0 else "   ⛔ CLASS BIAS"
                print(f"        {c['class']:<18} n {c['n']:>13,.0f}  expected {c['expected']:>13,.0f}"
                      f"  ratio {c['ratio']:.4f}  z {c['z']:+8.2f}{flag}")
    return comp_ok, field_ok


# ── self-test: every gate shown to FIRE on the defect it exists for ───────────────────────────────


def self_test() -> int:
    ok = fail = 0

    def check(name, cond):
        nonlocal ok, fail
        if cond:
            ok += 1
        else:
            fail += 1
            print(f"   ⛔ {name}")

    rng = np.random.default_rng(7)
    n = 4000
    eff = rng.uniform(50, 5000, n)
    cls = np.array(["R intron", "R exon", "B exon|exon", "B exon|intron"])[np.arange(n) % 4]
    rho = 0.05
    n_g = rng.poisson(rho * eff).astype(float)
    sel = np.ones(n, bool)

    v = gate_uniformity(n_g, eff, cls, sel)
    check("uniformity passes on a genuinely uniform field", v["ok"] and not v["vacuous"])

    # ① the counting-bug shape: ONE CLASS deposited at half weight
    bad = n_g.copy()
    bad[cls == "B exon|exon"] *= 0.5
    check("a half-weight class fires the class ratio-z",
          not gate_uniformity(bad, eff, cls, sel)["ok"])
    # ② a wrong DIVISOR on one class (eff doubled) fires the same way
    bad_e = eff.copy()
    bad_e[cls == "R intron"] *= 2.0
    check("a doubled divisor on one class fires", not gate_uniformity(n_g, bad_e, cls, sel)["ok"])
    # ③ vacuous is not a pass
    check("an empty field reports VACUOUS", gate_uniformity(np.zeros(n), eff, cls, sel)["vacuous"])

    n_n = rng.poisson(2.0, n).astype(float)
    n_m = rng.poisson(10.0, n).astype(float)
    count = n_g + n_n + n_m
    check("partition gate passes when exact", gate_partition(count, n_g, n_n, n_m)["ok"])
    check("one lost fragment fires the partition gate",
          not gate_partition(count, n_g - (np.arange(n) == 17), n_n, n_m)["ok"])

    admissible = np.arange(n) % 3 != 0
    n_n2 = np.where(admissible, n_n, 0.0)
    check("nascent-scope passes when confined", gate_nascent_scope(n_n2, admissible)["ok"])
    leak = n_n2.copy()
    leak[0] = 1.0  # slot 0 is inadmissible
    check("one intergenic nascent fragment fires", not gate_nascent_scope(leak, admissible)["ok"])

    z = np.zeros(n)
    check("zero gate passes on exact zeros", gate_zero(z, ~admissible, "t")["ok"])
    z[3] = 1e-9
    check("1e-9 of forbidden mass fires the zero gate", not gate_zero(z, ~admissible, "t")["ok"])
    check("zero gate on empty selection is VACUOUS", gate_zero(z, np.zeros(n, bool), "t")["vacuous"])

    print(f"\n   self-test: {ok} passed, {fail} failed")
    return 1 if fail else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--suite", type=Path, default=DEFAULT_SUITE)
    ap.add_argument("--index", type=Path, default=DEFAULT_INDEX)
    ap.add_argument("--condition", default=None)
    ap.add_argument("--out", type=Path, default=None,
                    help="write the certified per-slot table as .npz (default: <suite>/oracle_cache/<condition>/slot_truth.npz)")
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()
    if args.self_test:
        return self_test()

    index = TranscriptIndex.load(args.index)
    region_arrays = RegionArrays.from_index(index)
    conds = [args.condition] if args.condition else sorted(
        p.name for p in (args.suite / "scan_cache").iterdir()
    )
    bad = 0
    for c in conds:
        try:
            table, verdicts = derive(index, region_arrays, args.suite, c)
        except Exception as exc:  # noqa: BLE001 — **sum-to-full** failures surface here and must be a verdict
            print(f"\n== {c}\n   ⛔ gate sum-to-full / cache validation FAILED: {type(exc).__name__}: {exc}")
            bad += 1
            continue
        comp_ok, field_ok = report(c, verdicts)
        if comp_ok:
            out = args.out or (args.suite / "oracle_cache" / c / "slot_truth.npz")
            table["field_certified"] = bool(field_ok)
            np.savez_compressed(out, **{k: v for k, v in table.items() if k != "condition"})
            level = "COMPOSITION + FIELD" if field_ok else "COMPOSITION only (⛔ field gate failed — densities not certified)"
            print(f"   → table written [{level}]: {out}")
        if not comp_ok:
            bad += 1
            print("   ⛔ NOT CERTIFIED — no table written. Fix the gate, not the gate's bound.")
        elif not field_ok:
            bad += 1
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
