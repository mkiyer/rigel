"""Is the calibration result itself right, scored against an oracle calibration?

``P = calibrate(...)`` off the cached scan is compared with ``O``, the same ``CalibrationResult`` with only its six
deconvolved arrays replaced by the origin-split truth (`OracleTruth`, drained frame, sum-to-full gated); no EM and
no re-scan. ``O`` keeps ``P``'s efficiencies and lengths, so the capture-contracted length is invisible here: read
`ruler_vs_truth.py` beside a capture-ON arm. ``--set SECTION.FIELD=VALUE`` (repeatable) applies to both arms, so a
run prices one config value. The report, per stratum and tagged in scope or deferred: the pool ledger in
fragments, the composition per axis, the library gDNA fraction, the signed per-object error histogram, the
conditions ranked worst first (take the worst in-scope row), and a per-condition table.

Usage::

    python scripts/design/calibration_vs_oracle.py                           # the whole ladder
    python scripts/design/calibration_vs_oracle.py --conditions <name>       # one condition
    python scripts/design/calibration_vs_oracle.py --jobs 4                  # sharded, one report path
    python scripts/design/calibration_vs_oracle.py --set calibration.message_policy=silent   # price a policy on both arms
    python scripts/design/calibration_vs_oracle.py --set calibration.sweep_logodds_step=0.1   # any config value, both arms
    python scripts/design/calibration_vs_oracle.py --json rows.json          # write the rows and exit
    python scripts/design/calibration_vs_oracle.py --self-test               # no I/O
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import os
import tempfile
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402

from _shared import (  # noqa: E402
    DEFAULT_INDEX,
    DEFAULT_SUITE,
    OVERRIDE_FIELDS,
    SCRATCH,
    STRATA,
    add_set_option,
    check_override_fields,
    configured,
    oracle_cache_dir,
    panel_conditions,
    pool_ledger,
    refuses,
    run_json_shards,
    scope,
    set_args,
    set_field,
    sibling,
    stratum,
    stratum_selections,
)

OA = sibling("_oracle_arms.py")

from rigel.calibration import calibrate  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.substrate import CalibrationSubstrate  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402
from _scan_cache import calibration_inputs, read_scan_cache  # noqa: E402
from _oracle import OracleTruth  # noqa: E402

#: The two axes ``CalibrationResult`` deconvolves; the sj axis is certified RNA and is never split.
AXES = OA.AXES


#: Signed per-object error bucket edges in fragments, log-spaced, mirrored about an exact-zero bucket of its own.
_HIST_EDGES = (1e-9, 0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0)


def signed_histogram(delta: np.ndarray) -> dict:
    """The count and summed |error| per bucket of ``delta = gDNA_arm − gDNA_true`` (fragments, per object);
    positive is an over-call."""
    d = np.asarray(delta, np.float64)
    edges = np.asarray(_HIST_EDGES, np.float64)
    out = {"n_objects": int(d.size), "buckets": []}
    labels: list[tuple[str, np.ndarray]] = []
    for i in range(edges.size - 1, 0, -1):
        labels.append((f"< -{edges[i - 1]:g}" if i == edges.size - 1 else
                       f"[-{edges[i]:g}, -{edges[i - 1]:g})",
                       (d <= -edges[i - 1]) & ((d > -edges[i]) if i < edges.size - 1 else True)))
    labels.append(("EXACT 0", np.abs(d) < edges[0]))
    for i in range(1, edges.size):
        labels.append((f"[{edges[i - 1]:g}, {edges[i]:g})" if i < edges.size - 1 else
                       f"> +{edges[i - 1]:g}",
                       (d >= edges[i - 1]) & ((d < edges[i]) if i < edges.size - 1 else True)))
    for name, sel in labels:
        out["buckets"].append(
            {"label": name, "n": int(sel.sum()), "abs_err": float(np.abs(d[sel]).sum())}
        )
    out["n_under"] = int((d < -edges[0]).sum())
    out["n_over"] = int((d > edges[0]).sum())
    out["n_exact"] = int((np.abs(d) <= edges[0]).sum())
    return out


# ── one condition ────────────────────────────────────────────────────────────────────────────────


def score_axis_row(p_arm, o_arm, axis: str) -> dict:
    """``P`` against ``O`` on one axis: the ``AxisScore``, the signed per-object histogram over the objects with
    mass, and the gate that RNA mirrors gDNA per object."""
    def arr(arm, population):
        return np.asarray(getattr(arm, f"count_{population}_{axis}"), np.float64)

    gp, go, rp, ro = arr(p_arm, "gdna"), arr(o_arm, "gdna"), arr(p_arm, "rna"), arr(o_arm, "rna")
    row = dataclasses.asdict(OA.score_axis(gp, rp, go, ro))
    # the full distribution; objects with no mass have no answer and are excluded
    live = go + ro > 0.0
    row["hist"] = signed_histogram(gp[live] - go[live])
    # Δ_RNA ≡ −Δ_gDNA per object, since both arms carry the same per-object total; recorded as a gate
    row["rna_mirrors_gdna_max_dev"] = float(
        np.abs((rp - ro)[live] + (gp - go)[live]).max() if live.any() else 0.0
    )
    return row


def measure_condition(index, region_arrays, pipeline_config, suite: Path, oracle_cache: Path,
                      condition: str) -> dict:
    """calibrate once, build P and O, gate them, and score. One JSON-able row."""
    start = time.perf_counter()
    cache = read_scan_cache(Path(suite) / "scan_cache" / condition, index)
    lift: dict = {}
    kw = calibration_inputs(cache, index, lift_out=lift)
    # the drained frame: P, the substrate and the oracle's full side all read this payload, never `cache.payload`
    payload = kw["payload"]
    p_arm = calibrate(config=pipeline_config.calibration, **kw)

    # the cached parts are drained by replaying the whole's choices, and sum-to-full validates that lift
    oracle = OracleTruth.from_cached_parts(payload, OA.oracle_parts(oracle_cache, condition, index), lift)
    override = oracle.override_masses(region_arrays)
    check_override_fields(override)
    o_arm = dataclasses.replace(p_arm, **override)

    # both arms must carry the payload's own per-object totals, per axis
    substrate = CalibrationSubstrate.from_payload(payload, region_arrays)
    OA.check_same_basis("P", p_arm, substrate)
    OA.check_same_basis("O", o_arm, substrate)

    row = {
        "condition": condition,
        "stratum": list(stratum(condition)),
        "gdna_spliced_leak": oracle.gdna_spliced_leak,
        "lift_n_ambiguous": oracle.n_ambiguous,
        "seconds": time.perf_counter() - start,
        "library_f_gdna_P": OA.library_f_gdna(p_arm),
        "library_f_gdna_O": OA.library_f_gdna(o_arm),
        "axes": {axis: score_axis_row(p_arm, o_arm, axis) for axis in AXES},
    }
    # ── the pool ledger, in fragments ──
    truth = pool_ledger(Path(suite) / condition)
    row["pools"] = {
        "true_gdna": truth["gdna"],
        "true_mrna": truth["mrna"],
        "true_nrna": truth["nrna"],
        "true_rna": truth["mrna"] + truth["nrna"],
        # conserved fragment counts (each axis converted by its population's mass per crossing), not incidences
        "P_gdna": p_arm.library_gdna_fragments,
        "P_rna": p_arm.library_rna_fragments,
        "O_gdna": o_arm.library_gdna_fragments,
        "O_rna": o_arm.library_rna_fragments,
    }
    return row


# ── reporting ────────────────────────────────────────────────────────────────────────────────────


def _selections(rows: list[dict]) -> list:
    """The selections every table prints: one per stratum, tagged with its scope, zero-gDNA conditions excluded."""
    return stratum_selections((r["condition"] for r in rows), label=lambda st: f"{' x '.join(st)}  [{scope(st)}]")


def _agg_axis(scores: list[dict]) -> dict | None:
    """Sum ``AxisScore`` dicts into one (``None`` if empty); ``mwae`` is re-derived from the sums."""
    scores = [s for s in scores if s is not None]
    if not scores:
        return None
    out = {
        k: sum(s[k] for s in scores)
        for k in ("n_scored", "mass", "net_err", "abs_err", "over_call", "under_call")
    }
    out["mwae"] = out["abs_err"] / out["mass"] if out["mass"] > 0.0 else 0.0
    return out


def _print_gates(rows: list[dict]) -> None:
    print()
    print("=" * 118)
    print("  ⭐⭐⭐ CALIBRATION vs ORACLE CALIBRATION — the 0.8.0 metric")
    print(f"  {len(rows)} conditions   no solver, no EM, no BAM re-scan")
    print("=" * 118)
    print()
    print("  ✅ GATE  override_masses writes exactly the override field set")
    print("  ✅ GATE  P and O are on the payload's own per-object totals, both axes (check_same_basis)")
    worst_mirror = max((r["axes"][a]["rna_mirrors_gdna_max_dev"] for r in rows for a in AXES),
                       default=0.0)
    print(f"  ✅ GATE  Δ_RNA ≡ −Δ_gDNA per object on both axes (max deviation {worst_mirror:.3e}), so")
    print("           ONE signed number per object describes the whole error — reported once, not twice")


def _print_pool_ledger(rows: list[dict], selections) -> None:
    """⓪: what went in and what calibration says, in fragments."""
    print()
    print("  ⓪ ⭐⭐⭐ THE POOL LEDGER — FRAGMENT COUNTS, not ratios. What went in, and what calibration says.")
    print("     ⛔ THREE pools on the truth side, TWO on the answer side: `calibrate` splits gDNA from RNA")
    print("        and CANNOT split mature from nascent — that is the EM's job (`quant_accuracy.py`).")
    print(f"    {'condition / stratum':<44} {'true gDNA':>12} {'true mRNA':>12} {'true nRNA':>10} "
          f"{'P gDNA':>13} {'P RNA':>13} {'ΔgDNA vs O':>13} {'P/O gDNA':>9}")
    print("    " + "-" * 142)
    for title, sub in selections:
        if not sub:
            continue
        t = {k: sum(r["pools"][k] for r in sub) for k in
             ("true_gdna", "true_mrna", "true_nrna", "P_gdna", "P_rna", "O_gdna")}
        d = t["P_gdna"] - t["O_gdna"]
        ratio = t["P_gdna"] / t["O_gdna"] if t["O_gdna"] > 0 else float("nan")
        print(f"    {title:<44} {t['true_gdna']:>12,.0f} {t['true_mrna']:>12,.0f} "
              f"{t['true_nrna']:>10,.0f} {t['P_gdna']:>13,.0f} {t['P_rna']:>13,.0f} "
              f"{d:>+13,.0f} {ratio:>9.4f}")
    if all(r["pools"]["true_nrna"] == 0.0 for r in rows):
        print("    ⛔ EVERY condition has true nRNA = 0 — this panel cannot resolve the nascent axis at")
        print("       all. A gDNA-vs-nascent-vs-annotated verdict needs nascent-bearing conditions.")
    print("    ⚠ `P/O gDNA` is the like-for-like comparison. The simulator's own count is the OUTER")
    print("      reference and differs from O by the accumulator's retention, which is a different axis.")


def _print_composition(selections) -> None:
    """①: P against O, mass-weighted, per axis."""
    for axis in AXES:
        print()
        print(f"  ① COMPOSITION on the {axis.upper()} axis — P against O, mass-weighted")
        print(f"    {'stratum':<38} {'objects':>9} {'mass':>14} {'Σ|Δ gDNA|':>14} "
              f"{'mwae':>8} {'net':>14} {'over':>13} {'under':>13}")
        print("    " + "-" * 130)
        for title, sub in selections:
            a = _agg_axis([r["axes"][axis] for r in sub])
            if a is None:
                continue
            print(f"    {title:<38} {a['n_scored']:>9,} {a['mass']:>14,.0f} {a['abs_err']:>14,.0f} "
                  f"{a['mwae']:>8.4f} {a['net_err']:>+14,.0f} {a['over_call']:>13,.0f} "
                  f"{a['under_call']:>13,.0f}")


def _print_library_fraction(selections) -> None:
    """②: the library gDNA fraction per stratum, a ratio of sums over conditions of different depth."""
    print()
    print("  ② THE LIBRARY gDNA FRACTION — the one-number thermometer, per stratum")
    print(f"    {'stratum':<38} {'f_gdna P':>10} {'f_gdna O':>10} {'|Δ|':>10}")
    print("    " + "-" * 72)
    for title, sub in selections:
        if not sub:
            continue
        p = sum(r["library_f_gdna_P"] * r["axes"]["region"]["mass"] for r in sub)
        o = sum(r["library_f_gdna_O"] * r["axes"]["region"]["mass"] for r in sub)
        w = sum(r["axes"]["region"]["mass"] for r in sub)
        p, o = (p / w, o / w) if w > 0 else (float("nan"), float("nan"))
        print(f"    {title:<38} {p:>10.4f} {o:>10.4f} {abs(p - o):>10.4f}")


def _print_histograms(rows: list[dict], selections) -> None:
    """③: the full signed per-object error distribution, per axis."""
    for axis in AXES:
        print()
        print(f"  ③ ⭐⭐⭐ THE SIGNED ERROR DISTRIBUTION on the {axis.upper()} axis — zero IS the truth.")
        print("     Δ = gDNA_shipped − gDNA_true, in FRAGMENTS, per object. + is an OVER-call (gDNA")
        print("     claimed where the truth was RNA). Objects carrying no mass are excluded.")
        labels = [b["label"] for b in rows[0]["axes"][axis]["hist"]["buckets"]]
        print(f"    {'stratum':<38} " + " ".join(f"{lab:>13}" for lab in labels))
        print("    " + "-" * (38 + 14 * len(labels)))
        for title, sub in selections:
            if not sub:
                continue
            tot = [sum(r["axes"][axis]["hist"]["buckets"][i]["n"] for r in sub)
                   for i in range(len(labels))]
            print(f"    {title:<38} " + " ".join(f"{n:>13,}" for n in tot))
        print()
        print(f"    {'stratum':<38} {'under':>13} {'EXACT':>13} {'over':>13} {'objects':>13} "
              f"{'Σ|Δ| frags':>15}")
        print("    " + "-" * 110)
        for title, sub in selections:
            if not sub:
                continue
            h = [r["axes"][axis]["hist"] for r in sub]
            print(f"    {title:<38} {sum(x['n_under'] for x in h):>13,} "
                  f"{sum(x['n_exact'] for x in h):>13,} {sum(x['n_over'] for x in h):>13,} "
                  f"{sum(x['n_objects'] for x in h):>13,} "
                  f"{sum(r['axes'][axis]['abs_err'] for r in sub):>15,.0f}")


def _print_ranking(rows: list[dict]) -> None:
    """④: the conditions ranked by total Σ|Δ|, worst first, with the worst in-scope one named."""
    print()
    print("  ④ ⭐⭐⭐ SCENARIOS RANKED BY TOTAL Σ|Δ| IN FRAGMENTS (region + boundary) — worst first.")
    print("     ⛔ The scope tag is the ranking. The worst row overall is the DEFERRED stratum every")
    print("        time; take the worst IN-SCOPE row instead.")
    print(f"    {'#':>3} {'condition':<44} {'scope':<9} {'Σ|Δ| region':>13} {'Σ|Δ| bnd':>13} "
          f"{'Σ|Δ| TOTAL':>13} {'over':>12} {'under':>12}")
    print("    " + "-" * 126)
    ranked = sorted(
        rows, key=lambda r: -(r["axes"]["region"]["abs_err"] + r["axes"]["boundary"]["abs_err"])
    )
    for i, r in enumerate(ranked, 1):
        reg, bnd = r["axes"]["region"], r["axes"]["boundary"]
        print(f"    {i:>3} {r['condition']:<44} {scope(r['stratum']):<9} {reg['abs_err']:>13,.0f} "
              f"{bnd['abs_err']:>13,.0f} {reg['abs_err'] + bnd['abs_err']:>13,.0f} "
              f"{reg['over_call'] + bnd['over_call']:>12,.0f} "
              f"{reg['under_call'] + bnd['under_call']:>12,.0f}")
    in_scope = [r for r in ranked if scope(r["stratum"]) == "IN SCOPE"]
    if in_scope:
        w = in_scope[0]
        print(f"    ⭐ WORST IN-SCOPE: {w['condition']}  "
              f"Σ|Δ| = {w['axes']['region']['abs_err'] + w['axes']['boundary']['abs_err']:,.0f} fragments")


def _print_per_condition(rows: list[dict]) -> None:
    """⑤: one line per condition, ranked within a stratum only."""
    print()
    print("  ⑤ PER CONDITION — rank within a stratum, never across one")
    print(f"    {'condition':<44} {'mwae reg':>9} {'mwae bnd':>9} {'lib P':>8} {'lib O':>8} {'s':>6}")
    print("    " + "-" * 89)
    for r in sorted(rows, key=lambda x: (tuple(x["stratum"]), x["condition"])):
        print(f"    {r['condition']:<44} {r['axes']['region']['mwae']:>9.4f} "
              f"{r['axes']['boundary']['mwae']:>9.4f} {r['library_f_gdna_P']:>8.4f} "
              f"{r['library_f_gdna_O']:>8.4f} {r['seconds']:>6.1f}")
    print()
    print(f"  total wall clock {sum(r['seconds'] for r in rows):.1f} s over {len(rows)} conditions")


def report(rows: list[dict]) -> None:
    """The whole report, from the per-condition JSON rows; serial and sharded runs share this one path."""
    # each stratum's title and its rows
    selections = [(title, [r for r in rows if pred(r["condition"])]) for title, pred in _selections(rows)]
    _print_gates(rows)
    _print_pool_ledger(rows, selections)
    _print_composition(selections)
    _print_library_fraction(selections)
    _print_histograms(rows, selections)
    _print_ranking(rows)
    _print_per_condition(rows)


# ── the self-test: every comparator perturbed, with no I/O ───────────────────────────────────────


def _toy_calibration(n_regions: int = 24, n_boundaries: int = 20, n_sj: int = 4, seed: int = 7):
    """A minimal well-formed ``CalibrationResult`` for the self-test. No index, no payload, no I/O."""
    from rigel.calibration.result import CalibrationResult

    rng = np.random.default_rng(seed)
    ones_b = np.ones(n_boundaries)
    gr = rng.uniform(1.0, 9.0, n_regions)
    gb = rng.uniform(1.0, 9.0, n_boundaries)
    return CalibrationResult(
        count_gdna_region=gr,
        count_rna_region=rng.uniform(1.0, 9.0, n_regions),
        count_gdna_boundary=gb,
        count_rna_boundary=rng.uniform(1.0, 9.0, n_boundaries),
        count_rna_spliced_boundary=np.zeros(n_boundaries),
        boundary_mass_per_crossing=ones_b.copy(),
        count_rna_sj=np.zeros(n_sj),
        boundary_spliced_mass_per_crossing=ones_b.copy(),
        sj_mass_per_crossing=np.ones(n_sj),
        sj_src_region=np.zeros(n_sj, dtype=np.int64),
        sj_dst_region=np.ones(n_sj, dtype=np.int64),
        gdna_region_eff_len=rng.uniform(50.0, 500.0, n_regions),
        gdna_boundary_eff_len=np.full(n_boundaries, 120.0),
        gdna_boundary_conserved_len=np.full(n_boundaries, 120.0),
        gdna_frac_region=np.full(n_regions, 0.5),
        rna_pos_frac_region=np.full(n_regions, 0.25),
        rna_neg_frac_region=np.full(n_regions, 0.25),
        gdna_frac_boundary=np.full(n_boundaries, 0.5),
        rna_pos_frac_boundary=np.full(n_boundaries, 0.25),
        rna_neg_frac_boundary=np.full(n_boundaries, 0.25),
        gdna_density_global=0.1,
        gdna_bulk_density=None,
        gdna_capture_saturation=1.0,
        gdna_capture_efficiency_region=np.ones(n_regions),
        gdna_capture_efficiency_boundary=np.ones(n_boundaries),
        rna_sense_frac=0.5,
        n_regions=n_regions,
        n_boundaries=n_boundaries,
        n_sj=n_sj,
    )


def _check_score_axis(cal, check) -> None:
    """``score_axis`` of an arm against itself is exactly zero, and a one-ULP nudge makes it nonzero."""
    s = OA.score_axis(cal.count_gdna_region, cal.count_rna_region,
                      cal.count_gdna_region, cal.count_rna_region)
    check("score_axis(P, P) is exactly 0", s.abs_err == 0.0 and s.mwae == 0.0)
    nudged = np.array(cal.count_gdna_region, copy=True)
    nudged[3] = np.nextafter(nudged[3], np.inf)
    # the RNA side moves the opposite way, so the per-object total (the basis) is preserved
    rna_n = np.array(cal.count_rna_region, copy=True)
    rna_n[3] = np.nextafter(rna_n[3], -np.inf)
    s1 = OA.score_axis(nudged, rna_n, cal.count_gdna_region, cal.count_rna_region)
    check("score_axis resolves a ONE-ULP nudge", s1.abs_err > 0.0)


def _check_override_gate(cal, check) -> None:
    """The override field-set gate passes the real set and refuses a dropped field and an extra one."""
    full = {f: getattr(cal, f) for f in OVERRIDE_FIELDS}
    check("field-set gate passes the real set", not refuses(check_override_fields, full, exc=RuntimeError))
    for label, mutated in (("a DROPPED field", {k: v for k, v in list(full.items())[1:]}),
                           ("an EXTRA field", {**full, "gdna_density_global": 0.0})):
        check(f"field-set gate refuses {label}", refuses(check_override_fields, mutated, exc=RuntimeError))


def _check_aggregate(check) -> None:
    """The aggregate is a ratio of sums, not a mean of ratios; two scores of very different mass tell them apart."""
    a = {"n_scored": 1, "mass": 1.0, "net_err": 0.0, "abs_err": 1.0,
         "over_call": 0.0, "under_call": 0.0, "mwae": 1.0}
    b = {"n_scored": 1, "mass": 999.0, "net_err": 0.0, "abs_err": 0.0,
         "over_call": 0.0, "under_call": 0.0, "mwae": 0.0}
    agg = _agg_axis([a, b])
    check("aggregate mwae is Σ|err| / Σmass", np.isclose(agg["mwae"], 1.0 / 1000.0))
    check("aggregate mwae is NOT the mean of the rates", not np.isclose(agg["mwae"], 0.5))
    check("aggregate returns None on an empty selection", _agg_axis([]) is None)


def _check_histogram(check) -> None:
    """The signed histogram places a known vector exactly, keeps the sign, and drops no object."""
    known = np.array([0.0, 0.0, -5.0, -500.0, +2.0, +2000.0, -0.05, +0.05, +20000.0, -20000.0])
    h = signed_histogram(known)
    check("histogram conserves every object", sum(b["n"] for b in h["buckets"]) == known.size)
    check("histogram counts the exact zeros", h["n_exact"] == 2)
    check("histogram splits under/over correctly", (h["n_under"], h["n_over"]) == (4, 4))
    check("histogram Σ|err| equals Σ|input|",
          np.isclose(sum(b["abs_err"] for b in h["buckets"]), np.abs(known).sum()))
    # and moves an object between buckets when its value crosses an edge
    lo = signed_histogram(np.array([9.99]))
    hi = signed_histogram(np.array([10.01]))
    check("histogram resolves a bucket EDGE crossing",
          [b["n"] for b in lo["buckets"]] != [b["n"] for b in hi["buckets"]])
    check("histogram is SIGN-sensitive: +x and -x land in mirrored buckets",
          [b["n"] for b in signed_histogram(np.array([+7.0]))["buckets"]]
          == [b["n"] for b in signed_histogram(np.array([-7.0]))["buckets"]][::-1])


def _check_pool_ledger(check) -> None:
    """The pool ledger reads the file, and a pool it cannot find reads 0."""
    with tempfile.TemporaryDirectory() as td:
        (Path(td) / "truth_summary.json").write_text(
            json.dumps({"origin_counts": {"gdna": 7.0, "mrna": 11.0, "nrna": 13.0}})
        )
        led = pool_ledger(Path(td))
        check("pool ledger reads all three origin pools",
              (led["gdna"], led["mrna"], led["nrna"]) == (7.0, 11.0, 13.0))
        (Path(td) / "truth_summary.json").write_text(json.dumps({"origin_counts": {"gdna": 7.0}}))
        check("pool ledger reports a MISSING pool as 0, not as absent",
              pool_ledger(Path(td)) == {"gdna": 7.0, "mrna": 0.0, "nrna": 0.0})


def _check_scope(check) -> None:
    """Every stratum of the panel is covered exactly once by the release scope."""
    scopes = [scope(st) for st in STRATA]
    check("the scope table covers all four strata", len(STRATA) == 4 and "APART" not in scopes)
    check("three strata are IN SCOPE and one is DEFERRED",
          sorted(scopes) == ["DEFERRED", "IN SCOPE", "IN SCOPE", "IN SCOPE"])


def _check_set(check) -> None:
    """``--set``: the parser types from the field, touches nothing else, and refuses what it cannot spell."""
    base = PipelineConfig()
    knob = set_field(base, "calibration.sweep_block_slots=90")
    check("--set types an int field from the field",
          knob.calibration.sweep_block_slots == 90 and type(knob.calibration.sweep_block_slots) is int)
    # two settings in sequence: the first must survive the second
    knob2 = set_field(knob, "calibration.sweep_logodds_step=0.1")
    check("--set leaves every other field identical, and an earlier --set survives a later one",
          knob2.calibration.sweep_block_slots == 90
          and knob2.calibration.sweep_logodds_step == 0.1
          and dataclasses.replace(knob2.calibration, sweep_block_slots=base.calibration.sweep_block_slots,
                                  sweep_logodds_step=base.calibration.sweep_logodds_step)
          == base.calibration
          and dataclasses.replace(knob2, calibration=base.calibration) == base)
    check("--set fills an `int | None` field as an int",
          set_field(base, "calibration.sweep_block_slots=1000").calibration.sweep_block_slots == 1000)
    check("--set refuses an unknown field", refuses(set_field, base, "calibration.no_such_field=1"))
    check("--set refuses an unknown section", refuses(set_field, base, "nowhere.sweep_block_slots=1"))
    check("--set refuses a value the field's type cannot take",
          refuses(set_field, base, "calibration.sweep_block_slots=sixty"))
    # a ``Literal`` field (the EM's mode and warm start) is filled by membership
    check("--set fills a Literal field by membership",
          set_field(base, "em.mode=map").em.mode == "map"
          and set_field(base, "em.warm_start=uniform").em.warm_start == "uniform")
    check("--set refuses a value outside a Literal field's set", refuses(set_field, base, "em.mode=maximum"))


def self_test() -> int:
    """Perturb every comparator and require it to fire."""
    checks: list[tuple[str, bool]] = []

    def check(name: str, ok: bool):
        checks.append((name, bool(ok)))

    cal = _toy_calibration()
    _check_score_axis(cal, check)
    _check_override_gate(cal, check)
    _check_aggregate(check)
    _check_histogram(check)
    _check_pool_ledger(check)
    _check_scope(check)
    _check_set(check)

    width = max(len(n) for n, _ in checks)
    for name, ok in checks:
        print(f"  {'✅' if ok else '⛔'} {name:<{width}}")
    n_ok = sum(1 for _n, ok in checks if ok)
    print(f"\n  {n_ok}/{len(checks)} self-test gates fired as specified")
    return 0 if n_ok == len(checks) else 1


# ── cli ──────────────────────────────────────────────────────────────────────────────────────────


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--suite", type=Path, default=DEFAULT_SUITE)
    ap.add_argument("--index", type=Path, default=DEFAULT_INDEX)
    ap.add_argument("--conditions", nargs="*", default=None)
    ap.add_argument("--oracle-cache", type=Path, default=None,
                    help="defaults to <suite>/oracle_cache when that directory exists")
    ap.add_argument("--work-dir", type=Path, default=SCRATCH / "rigel_cal_vs_oracle")
    ap.add_argument("--json", type=Path, default=None, help="write the per-condition rows and exit")
    ap.add_argument("--jobs", type=int, default=1)
    add_set_option(ap, "override one PipelineConfig field on BOTH arms, typed from the field; repeatable "
                       "(e.g. --set calibration.sweep_logodds_step=0.1)")
    ap.add_argument("--self-test", action="store_true", help="perturb every comparator; no I/O")
    args = ap.parse_args()

    if args.self_test:
        return self_test()

    names = panel_conditions(args.suite, args.conditions)
    if not names:
        raise SystemExit(f"⛔ no condition with a sim_oracle.bam under {args.suite}")
    cache = oracle_cache_dir(args.suite, args.oracle_cache)
    if cache is None:
        raise SystemExit("⛔ --oracle-cache is required; this script refuses to invent a truth")

    if args.jobs > 1:
        # merged in ``names`` order, so a sharded run reports as a serial one
        merged = run_json_shards(
            Path(__file__).resolve(), names, args.jobs, args.work_dir,
            lambda out, shard: ["--suite", str(args.suite), "--index", str(args.index),
                                "--oracle-cache", str(cache), "--json", str(out), "--conditions", *shard,
                                *set_args(args.settings)],
        )
        order = {name: i for i, name in enumerate(names)}
        merged.sort(key=lambda r: order[r["condition"]])
        return _emit(merged, args.json)

    t0 = time.perf_counter()
    index = TranscriptIndex.load(str(args.index))
    region_arrays = RegionArrays.from_index(index)
    # applied to both arms: P and O share one payload and differ only in the six deconvolved arrays
    pipeline_config = configured(args.settings)
    if args.json is None:
        for spec in args.settings:
            print(f"⭐ --set {spec} on BOTH arms")
        print(f"index + region arrays loaded in {time.perf_counter() - t0:.2f} s", flush=True)

    rows = []
    for name in names:
        rows.append(
            measure_condition(index, region_arrays, pipeline_config, args.suite, cache, name)
        )
        if args.json is None:
            print(f"  scored {name}  ({rows[-1]['seconds']:.1f} s)", flush=True)

    return _emit(rows, args.json)


def _emit(rows: list[dict], json_path: Path | None) -> int:
    """The per-condition rows, written to ``json_path`` when one is given, else reported."""
    if json_path is not None:
        json_path.parent.mkdir(parents=True, exist_ok=True)
        json_path.write_text(json.dumps(rows))
        return 0
    report(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
