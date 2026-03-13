#!/usr/bin/env python3
"""Three-way comparison: Original Baseline vs Phase 1 vs T+N+2.

Produces comprehensive analysis with gDNA accuracy tables, nRNA accuracy,
mRNA regression check, and deep root cause analysis.
"""
import csv
import sys
from collections import defaultdict

import numpy as np

TRANSCRIPTS = ["TA1", "TA2", "TA3", "TA4", "TB1", "TB2", "TCF", "TCR", "TD", "TE"]
NRNA_LABELS = ["NTA1", "NTA2", "NTA34", "NTB12", "NTCR", "NTE"]

PATTERN_LABELS = {
    0: "TA1_only_no_nrna",
    1: "TA4_only_no_nrna",
    2: "mixed_A_no_nrna",
    3: "convergent_CF_CR_no_nrna",
    4: "TA1_mod_nrna",
    5: "TA4_mod_nrna",
    6: "mixed_A_mod_nrna",
    7: "convergent_mod_nrna",
    8: "TA1_heavy_nrna",
    9: "TA4_heavy_nrna",
    10: "mixed_A_heavy_nrna",
    11: "convergent_heavy_nrna",
    12: "everything_moderate",
    13: "TA1+TCF+heavy_nrna",
    14: "TA3+TCR+heavy_nrna",
    15: "nrna_only_no_mrna",
}


def classify_pattern(row):
    ta1 = float(row["TA1"])
    ta2 = float(row["TA2"])
    ta3 = float(row["TA3"])
    ta4 = float(row["TA4"])
    tcf = float(row["TCF"])
    tcr = float(row["TCR"])
    nta1 = float(row["NTA1"])
    nta34 = float(row["NTA34"])
    ntcr = float(row["NTCR"])
    total_mrna = ta1 + ta2 + ta3 + ta4 + tcf + tcr
    total_nrna = nta1 + float(row["NTA2"]) + nta34 + float(row["NTB12"]) + ntcr
    if total_mrna == 0 and total_nrna > 0:
        return 15
    if ta1 == 500 and ta2 == 300 and ta3 == 200 and ta4 == 100:
        return 12
    if ta1 == 1000 and tcf == 1000:
        return 13
    if ta3 == 1000 and tcr == 1000:
        return 14
    if total_nrna == 0:
        if ta1 == 1000 and ta2 == 0:
            return 0
        if ta4 == 1000 and ta1 == 0:
            return 1
        if ta1 == 500 and ta2 == 300:
            return 2
        if tcf == 500 and tcr == 500:
            return 3
    else:
        if ta1 == 1000 and nta1 == 200:
            return 4
        if ta4 == 1000 and nta34 == 200:
            return 5
        if ta1 == 500 and nta1 == 200:
            return 6
        if tcf == 500 and ntcr == 200:
            return 7
        if ta1 == 1000 and nta1 == 1000:
            return 8
        if ta4 == 1000 and nta34 == 1000:
            return 9
        if ta1 == 500 and nta1 == 1000:
            return 10
        if tcf == 500 and ntcr == 1000:
            return 11
    return -1


def classify_nrna_level(row):
    vals = [float(row[k]) for k in NRNA_LABELS]
    total = sum(vals)
    if total == 0:
        return "none"
    pat = classify_pattern(row)
    if pat in (4, 5, 6, 7):
        return "moderate"
    if pat in (8, 9, 10, 11):
        return "heavy"
    if pat in (12, 13, 14, 15):
        return "complex"
    return "unknown"


def load_sweep(path):
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)
    for r in rows:
        r["_pattern"] = classify_pattern(r)
        r["_nrna_level"] = classify_nrna_level(r)
        r["_ss"] = float(r["strand_specificity"])
        r["_gdna_frac"] = float(r["gdna_fraction"])
    return rows


def safe_rel(obs, exp):
    if exp == 0:
        return None
    return (obs - exp) / exp


def pct(v):
    if v is None:
        return "n/a"
    return f"{v:+.1f}%"


def main():
    if len(sys.argv) < 4:
        print("Usage: deep_compare_three_way.py <baseline.tsv> <phase1.tsv> <tn2.tsv>")
        sys.exit(1)

    labels = ["Baseline", "Phase 1", "T+N+2"]
    datasets = [load_sweep(p) for p in sys.argv[1:4]]

    for i, (label, ds) in enumerate(zip(labels, datasets)):
        print(f"{label}: {len(ds)} runs ({sys.argv[i+1]})")

    ss_vals = [0.5, 0.75, 0.9, 1.0]
    gdna_fracs = [0.1, 0.3, 0.5]
    nrna_levels = ["none", "moderate", "heavy", "complex"]

    # ==================================================================================
    # SECTION 1: gDNA ACCURACY — THE PRIMARY TARGET
    # ==================================================================================
    print("\n" + "=" * 120)
    print("  1. gDNA ACCURACY (mean signed relative error) — THREE-WAY COMPARISON")
    print("=" * 120)

    for nlevel in nrna_levels:
        print(f"\n  nRNA level: {nlevel}")
        header = "  gdna_f"
        for ss in ss_vals:
            header += f"     ss={ss}                          "
        print(header)
        sub = "        "
        for _ in ss_vals:
            sub += f"{'Baseline':>10s}{'Phase1':>10s}{'T+N+2':>10s}   "
        print(sub)
        print("  " + "-" * 110)

        for gf in gdna_fracs:
            line = f"     {gf}"
            for ss in ss_vals:
                vals = []
                for ds in datasets:
                    matching = [r for r in ds if r["_nrna_level"] == nlevel
                                and abs(r["_ss"] - ss) < 0.01
                                and abs(r["_gdna_frac"] - gf) < 0.01]
                    if matching:
                        rels = []
                        for r in matching:
                            exp = float(r["gdna_expected"])
                            obs = float(r["gdna_observed"])
                            rv = safe_rel(obs, exp)
                            if rv is not None:
                                rels.append(rv * 100)
                        vals.append(np.mean(rels) if rels else None)
                    else:
                        vals.append(None)
                for v in vals:
                    line += f"  {pct(v):>8s}"
                line += "   "
            print(line)

    # ==================================================================================
    # SECTION 2: gDNA ACCURACY AT SS≥0.75 — SUCCESS CRITERIA CHECK
    # ==================================================================================
    print("\n" + "=" * 120)
    print("  2. gDNA ACCURACY AT SS≥0.75 — SUCCESS CRITERIA: |mean error| < 15%")
    print("=" * 120)

    for label, ds in zip(labels, datasets):
        matching = [r for r in ds if r["_ss"] >= 0.75 and r["_gdna_frac"] > 0]
        rels = []
        for r in matching:
            exp = float(r["gdna_expected"])
            obs = float(r["gdna_observed"])
            rv = safe_rel(obs, exp)
            if rv is not None:
                rels.append(rv * 100)
        mean_abs = np.mean(np.abs(rels))
        mean_signed = np.mean(rels)
        median = np.median(rels)
        p25 = np.percentile(rels, 25)
        p75 = np.percentile(rels, 75)
        n_under_15 = sum(1 for r in rels if abs(r) < 15)
        n_total = len(rels)
        print(f"\n  {label}:")
        print(f"    Mean |error|:   {mean_abs:.1f}%")
        print(f"    Mean signed:    {mean_signed:+.1f}%")
        print(f"    Median:         {median:+.1f}%")
        print(f"    IQR:            [{p25:+.1f}%, {p75:+.1f}%]")
        print(f"    Within ±15%:    {n_under_15}/{n_total} ({100*n_under_15/n_total:.0f}%)")
        pass_fail = "PASS" if mean_abs < 15 else "FAIL"
        print(f"    SUCCESS CRITERIA: {pass_fail}")

    # ==================================================================================
    # SECTION 3: nRNA ACCURACY — REGRESSION CHECK
    # ==================================================================================
    print("\n" + "=" * 120)
    print("  3. nRNA ACCURACY (mean signed relative error) — REGRESSION CHECK")
    print("=" * 120)

    for nlevel in ["moderate", "heavy", "complex"]:
        print(f"\n  nRNA level: {nlevel}")
        header = "  gdna_f"
        for ss in ss_vals:
            header += f"     ss={ss}                          "
        print(header)
        sub = "        "
        for _ in ss_vals:
            sub += f"{'Baseline':>10s}{'Phase1':>10s}{'T+N+2':>10s}   "
        print(sub)
        print("  " + "-" * 110)

        for gf in [0.0, 0.1, 0.3, 0.5]:
            line = f"     {gf}"
            for ss in ss_vals:
                vals = []
                for ds in datasets:
                    matching = [r for r in ds if r["_nrna_level"] == nlevel
                                and abs(r["_ss"] - ss) < 0.01
                                and abs(r["_gdna_frac"] - gf) < 0.01]
                    if matching:
                        rels = []
                        for r in matching:
                            exp = float(r["nrna_expected"])
                            obs = float(r["nrna_observed"])
                            rv = safe_rel(obs, exp)
                            if rv is not None:
                                rels.append(rv * 100)
                        vals.append(np.mean(rels) if rels else None)
                    else:
                        vals.append(None)
                for v in vals:
                    line += f"  {pct(v):>8s}"
                line += "   "
            print(line)

    # ==================================================================================
    # SECTION 4: mRNA ACCURACY — REGRESSION CHECK
    # ==================================================================================
    print("\n" + "=" * 120)
    print("  4. TOTAL mRNA ACCURACY (mean |relative error|) — REGRESSION CHECK")
    print("=" * 120)

    for nlevel in nrna_levels:
        print(f"\n  nRNA level: {nlevel}")
        header = "  gdna_f"
        for ss in ss_vals:
            header += f"     ss={ss}                          "
        print(header)
        sub = "        "
        for _ in ss_vals:
            sub += f"{'Baseline':>10s}{'Phase1':>10s}{'T+N+2':>10s}   "
        print(sub)
        print("  " + "-" * 110)

        for gf in [0.0, 0.1, 0.3, 0.5]:
            line = f"     {gf}"
            for ss in ss_vals:
                vals = []
                for ds in datasets:
                    matching = [r for r in ds if r["_nrna_level"] == nlevel
                                and abs(r["_ss"] - ss) < 0.01
                                and abs(r["_gdna_frac"] - gf) < 0.01]
                    if matching:
                        rels = []
                        for r in matching:
                            exp = float(r["total_mrna_expected"])
                            obs = float(r["total_mrna_observed"])
                            rv = safe_rel(obs, exp)
                            if rv is not None:
                                rels.append(abs(rv) * 100)
                        vals.append(np.mean(rels) if rels else None)
                    else:
                        vals.append(None)
                for v in vals:
                    if v is not None:
                        line += f"  {v:>7.1f}%"
                    else:
                        line += f"  {'n/a':>7s} "
                line += "   "
            print(line)

    # ==================================================================================
    # SECTION 5: DEEP gDNA ANALYSIS — BY PATTERN AND CONDITION
    # ==================================================================================
    print("\n" + "=" * 120)
    print("  5. gDNA ERROR BY EXPRESSION PATTERN (SS≥0.75, all gdna fracs)")
    print("=" * 120)

    for pat_id in sorted(PATTERN_LABELS.keys()):
        pat_name = PATTERN_LABELS[pat_id]
        print(f"\n  Pattern {pat_id}: {pat_name}")
        header = "  gdna_f"
        for ss in [0.75, 0.9, 1.0]:
            header += f"     ss={ss}                          "
        print(header)
        sub = "        "
        for _ in [0.75, 0.9, 1.0]:
            sub += f"{'Baseline':>10s}{'Phase1':>10s}{'T+N+2':>10s}   "
        print(sub)
        print("  " + "-" * 90)

        for gf in gdna_fracs:
            line = f"     {gf}"
            for ss in [0.75, 0.9, 1.0]:
                vals = []
                for ds in datasets:
                    matching = [r for r in ds if r["_pattern"] == pat_id
                                and abs(r["_ss"] - ss) < 0.01
                                and abs(r["_gdna_frac"] - gf) < 0.01]
                    if matching:
                        rels = []
                        for r in matching:
                            exp = float(r["gdna_expected"])
                            obs = float(r["gdna_observed"])
                            rv = safe_rel(obs, exp)
                            if rv is not None:
                                rels.append(rv * 100)
                        vals.append(np.mean(rels) if rels else None)
                    else:
                        vals.append(None)
                for v in vals:
                    line += f"  {pct(v):>8s}"
                line += "   "
            print(line)

    # ==================================================================================
    # SECTION 6: nRNA-gDNA CROSSTALK ANALYSIS
    # ==================================================================================
    print("\n" + "=" * 120)
    print("  6. nRNA-gDNA CROSSTALK: nRNA error DIRECTION vs gDNA error DIRECTION (SS≥0.75)")
    print("=" * 120)
    print("  (Positive crosstalk = nRNA over-counted while gDNA under-counted, or vice versa)")

    for label, ds in zip(labels, datasets):
        print(f"\n  {label}:")
        # Categorize runs
        cats = {"nrna_over_gdna_under": 0, "nrna_under_gdna_over": 0,
                "both_over": 0, "both_under": 0, "other": 0}
        total = 0
        for r in ds:
            if r["_ss"] < 0.75 or r["_gdna_frac"] == 0 or r["_nrna_level"] == "none":
                continue
            nrna_exp = float(r["nrna_expected"])
            nrna_obs = float(r["nrna_observed"])
            gdna_exp = float(r["gdna_expected"])
            gdna_obs = float(r["gdna_observed"])
            if nrna_exp == 0 or gdna_exp == 0:
                continue
            nrna_err = (nrna_obs - nrna_exp) / nrna_exp
            gdna_err = (gdna_obs - gdna_exp) / gdna_exp
            total += 1
            if nrna_err > 0.05 and gdna_err < -0.05:
                cats["nrna_over_gdna_under"] += 1
            elif nrna_err < -0.05 and gdna_err > 0.05:
                cats["nrna_under_gdna_over"] += 1
            elif nrna_err > 0.05 and gdna_err > 0.05:
                cats["both_over"] += 1
            elif nrna_err < -0.05 and gdna_err < -0.05:
                cats["both_under"] += 1
            else:
                cats["other"] += 1
        if total > 0:
            print(f"    nRNA↑ gDNA↓ (crosstalk): {cats['nrna_over_gdna_under']:>4d}/{total} ({100*cats['nrna_over_gdna_under']/total:.0f}%)")
            print(f"    nRNA↓ gDNA↑ (inv-cross): {cats['nrna_under_gdna_over']:>4d}/{total} ({100*cats['nrna_under_gdna_over']/total:.0f}%)")
            print(f"    Both↑ (absorption fail): {cats['both_over']:>4d}/{total} ({100*cats['both_over']/total:.0f}%)")
            print(f"    Both↓ (mass-loss):       {cats['both_under']:>4d}/{total} ({100*cats['both_under']/total:.0f}%)")
            print(f"    Within ±5% (accurate):   {cats['other']:>4d}/{total} ({100*cats['other']/total:.0f}%)")

    # ==================================================================================
    # SECTION 7: OVERESTIMATION ANALYSIS — T+N+2 gDNA FLIP
    # ==================================================================================
    print("\n" + "=" * 120)
    print("  7. gDNA OVERESTIMATION ANALYSIS — WHERE T+N+2 FLIPS FROM UNDER TO OVER")
    print("=" * 120)

    # For each condition, show where Phase 1 was negative (underestimate) and T+N+2 is positive (overestimate)
    flip_count = 0
    flip_worse = 0
    flip_better = 0
    flip_details = []

    base = datasets[0]
    ph1 = datasets[1]
    tn2 = datasets[2]

    for i in range(len(base)):
        if base[i]["_gdna_frac"] == 0:
            continue
        gdna_exp_b = float(base[i]["gdna_expected"])
        if gdna_exp_b == 0:
            continue
        gdna_err_b = (float(base[i]["gdna_observed"]) - gdna_exp_b) / gdna_exp_b * 100
        gdna_err_p = (float(ph1[i]["gdna_observed"]) - float(ph1[i]["gdna_expected"])) / float(ph1[i]["gdna_expected"]) * 100
        gdna_err_t = (float(tn2[i]["gdna_observed"]) - float(tn2[i]["gdna_expected"])) / float(tn2[i]["gdna_expected"]) * 100

        # Where baseline was underestimating and T+N+2 flipped to overestimate
        if gdna_err_b < -5 and gdna_err_t > 5:
            flip_count += 1
            # Is T+N+2 error larger in magnitude?
            if abs(gdna_err_t) > abs(gdna_err_b):
                flip_worse += 1
            else:
                flip_better += 1
            flip_details.append({
                "ss": base[i]["_ss"],
                "gf": base[i]["_gdna_frac"],
                "nrna": base[i]["_nrna_level"],
                "pattern": PATTERN_LABELS.get(base[i]["_pattern"], "?"),
                "base": gdna_err_b,
                "ph1": gdna_err_p,
                "tn2": gdna_err_t,
            })

    n_total_gdna = sum(1 for r in base if r["_gdna_frac"] > 0 and float(r["gdna_expected"]) > 0)
    print(f"\n  Runs where baseline underestimated (>5%) but T+N+2 overestimates (>5%): {flip_count}/{n_total_gdna}")
    if flip_count:
        print(f"    Magnitude worse (|T+N+2| > |base|): {flip_worse}")
        print(f"    Magnitude better (|T+N+2| < |base|): {flip_better}")

        # Group by condition
        by_cond = defaultdict(list)
        for d in flip_details:
            key = (d["ss"], d["gf"], d["nrna"])
            by_cond[key].append(d)

        print(f"\n  {'SS':>5s} {'gdna':>5s} {'nRNA':>10s} {'count':>6s} {'base_mean':>10s} {'ph1_mean':>10s} {'tn2_mean':>10s}")
        print("  " + "-" * 60)
        for (ss, gf, nrna), items in sorted(by_cond.items()):
            n = len(items)
            bm = np.mean([d["base"] for d in items])
            pm = np.mean([d["ph1"] for d in items])
            tm = np.mean([d["tn2"] for d in items])
            print(f"  {ss:>5.2f} {gf:>5.1f} {nrna:>10s} {n:>6d} {bm:>+9.1f}% {pm:>+9.1f}% {tm:>+9.1f}%")

    # ==================================================================================
    # SECTION 8: SUMMARY SCORECARD
    # ==================================================================================
    print("\n" + "=" * 120)
    print("  8. SUMMARY SCORECARD")
    print("=" * 120)

    def calc_metric(ds, filter_fn, metric_fn):
        vals = []
        for r in ds:
            if filter_fn(r):
                v = metric_fn(r)
                if v is not None:
                    vals.append(v)
        return np.mean(vals) if vals else None

    metrics = [
        ("mRNA accuracy (no nRNA, gdna=0)",
         lambda r: r["_nrna_level"] == "none" and r["_gdna_frac"] == 0,
         lambda r: abs(safe_rel(float(r["total_mrna_observed"]), float(r["total_mrna_expected"])) or 0) * 100),
        ("mRNA accuracy (SS≥0.75, gdna=0)",
         lambda r: r["_ss"] >= 0.75 and r["_gdna_frac"] == 0,
         lambda r: abs(safe_rel(float(r["total_mrna_observed"]), float(r["total_mrna_expected"])) or 0) * 100),
        ("mRNA accuracy (SS≥0.75, all gdna)",
         lambda r: r["_ss"] >= 0.75,
         lambda r: abs(safe_rel(float(r["total_mrna_observed"]), float(r["total_mrna_expected"])) or 0) * 100),
        ("gDNA |error| (SS≥0.75)",
         lambda r: r["_ss"] >= 0.75 and r["_gdna_frac"] > 0,
         lambda r: abs((safe_rel(float(r["gdna_observed"]), float(r["gdna_expected"])) or 0) * 100)),
        ("gDNA signed error (SS≥0.75)",
         lambda r: r["_ss"] >= 0.75 and r["_gdna_frac"] > 0,
         lambda r: (safe_rel(float(r["gdna_observed"]), float(r["gdna_expected"])) or 0) * 100),
        ("gDNA |error| (SS=1.0)",
         lambda r: abs(r["_ss"] - 1.0) < 0.01 and r["_gdna_frac"] > 0,
         lambda r: abs((safe_rel(float(r["gdna_observed"]), float(r["gdna_expected"])) or 0) * 100)),
        ("gDNA |error| (SS=0.9)",
         lambda r: abs(r["_ss"] - 0.9) < 0.01 and r["_gdna_frac"] > 0,
         lambda r: abs((safe_rel(float(r["gdna_observed"]), float(r["gdna_expected"])) or 0) * 100)),
        ("gDNA |error| (SS=0.75)",
         lambda r: abs(r["_ss"] - 0.75) < 0.01 and r["_gdna_frac"] > 0,
         lambda r: abs((safe_rel(float(r["gdna_observed"]), float(r["gdna_expected"])) or 0) * 100)),
        ("nRNA signed error (SS≥0.75, mod, gdna=0)",
         lambda r: r["_ss"] >= 0.75 and r["_nrna_level"] == "moderate" and r["_gdna_frac"] == 0,
         lambda r: (safe_rel(float(r["nrna_observed"]), float(r["nrna_expected"])) or 0) * 100),
        ("nRNA signed error (SS≥0.75, mod, all gdna)",
         lambda r: r["_ss"] >= 0.75 and r["_nrna_level"] == "moderate",
         lambda r: (safe_rel(float(r["nrna_observed"]), float(r["nrna_expected"])) or 0) * 100),
        ("nRNA signed error (SS≥0.75, heavy, all gdna)",
         lambda r: r["_ss"] >= 0.75 and r["_nrna_level"] == "heavy",
         lambda r: (safe_rel(float(r["nrna_observed"]), float(r["nrna_expected"])) or 0) * 100),
        ("nRNA |error| (SS=0.5)",
         lambda r: abs(r["_ss"] - 0.5) < 0.01 and r["_nrna_level"] != "none",
         lambda r: abs((safe_rel(float(r["nrna_observed"]), float(r["nrna_expected"])) or 0) * 100)),
    ]

    print(f"\n  {'Metric':<50s} {'Baseline':>10s} {'Phase 1':>10s} {'T+N+2':>10s} {'Δ(T+N+2-Base)':>14s}")
    print("  " + "-" * 100)
    for name, filt, metric in metrics:
        vals = []
        for ds in datasets:
            v = calc_metric(ds, filt, metric)
            vals.append(v)
        delta = None
        if vals[0] is not None and vals[2] is not None:
            delta = vals[2] - vals[0]
        line = f"  {name:<50s}"
        for v in vals:
            line += f"  {v:>8.1f}%" if v is not None else f"  {'n/a':>8s} "
        if delta is not None:
            line += f"  {delta:>+12.1f}%"
        print(line)

    # ==================================================================================
    # SECTION 9: gDNA DIRECTION SHIFT ANALYSIS
    # ==================================================================================
    print("\n" + "=" * 120)
    print("  9. gDNA ERROR DIRECTION DISTRIBUTION (SS≥0.75)")
    print("=" * 120)

    for label, ds in zip(labels, datasets):
        matching = [r for r in ds if r["_ss"] >= 0.75 and r["_gdna_frac"] > 0]
        under_20 = 0
        under_10 = 0
        within_5 = 0
        over_10 = 0
        over_20 = 0
        rels = []
        for r in matching:
            exp = float(r["gdna_expected"])
            obs = float(r["gdna_observed"])
            if exp == 0:
                continue
            err = (obs - exp) / exp * 100
            rels.append(err)
            if err < -20:
                under_20 += 1
            elif err < -10:
                under_10 += 1
            elif err <= 5:
                within_5 += 1
            elif err <= 20:
                over_10 += 1
            else:
                over_20 += 1
        n = len(rels)
        print(f"\n  {label} (n={n}):")
        print(f"    < -20%:   {under_20:>4d} ({100*under_20/n:.0f}%)")
        print(f"    -20..-10: {under_10:>4d} ({100*under_10/n:.0f}%)")
        print(f"    -10..+5:  {within_5:>4d} ({100*within_5/n:.0f}%)")
        print(f"    +5..+20:  {over_10:>4d} ({100*over_10/n:.0f}%)")
        print(f"    > +20%:   {over_20:>4d} ({100*over_20/n:.0f}%)")

    # ==================================================================================
    # SECTION 10: PATTERN 15 (nRNA-only) STRESS TEST
    # ==================================================================================
    print("\n" + "=" * 120)
    print("  10. PATTERN 15 STRESS TEST: nRNA only (no mRNA expressed)")
    print("=" * 120)
    print(f"\n  {'SS':>5s} {'gdna':>5s}  ", end="")
    for label in labels:
        print(f"  {'mrna':>6s} {'nrna%':>8s} {'gdna%':>8s}", end="")
    print()
    print("  " + "-" * 100)

    for ss in ss_vals:
        for gf in [0.0, 0.1, 0.3, 0.5]:
            line = f"  {ss:>5.2f} {gf:>5.1f}  "
            for ds in datasets:
                matching = [r for r in ds if r["_pattern"] == 15
                            and abs(r["_ss"] - ss) < 0.01
                            and abs(r["_gdna_frac"] - gf) < 0.01]
                if matching:
                    r = matching[0]
                    mrna_obs = float(r["total_mrna_observed"])
                    nrna_exp = float(r["nrna_expected"])
                    nrna_obs = float(r["nrna_observed"])
                    nrna_rel = safe_rel(nrna_obs, nrna_exp)
                    gdna_exp = float(r["gdna_expected"])
                    gdna_obs = float(r["gdna_observed"])
                    gdna_rel = safe_rel(gdna_obs, gdna_exp) if gdna_exp > 0 else None
                    line += f"  {mrna_obs:>6.0f} {pct(nrna_rel*100 if nrna_rel is not None else None):>8s} {pct(gdna_rel*100 if gdna_rel is not None else None):>8s}"
                else:
                    line += f"  {'n/a':>6s} {'n/a':>8s} {'n/a':>8s}"
            print(line)

    print("\n" + "=" * 120)
    print("  END OF ANALYSIS")
    print("=" * 120)


if __name__ == "__main__":
    main()
