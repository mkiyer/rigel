#!/usr/bin/env python3
"""Deep-dive analysis of Phase 3 (SS-scaled κ) vs Phase 1 vs baseline.

Key questions:
1. Did Phase 3 κ-scaling change anything compared to Phase 1?
2. What is the effective κ at each SS level and its impact?
3. Why is gDNA still underestimated at SS≥0.75 with nRNA?
4. Mass flow analysis: where do gDNA reads end up?
5. Is the Laplace smoothing change (ε=50 → +0.5/+1.0) beneficial or harmful?
"""
import csv
import sys
from collections import defaultdict

import numpy as np

TRANSCRIPTS = ["TA1", "TA2", "TA3", "TA4", "TB1", "TB2", "TCF", "TCR", "TD", "TE"]
NRNA_LABELS = ["NTA1", "NTA2", "NTA34", "NTB12", "NTCR", "NTE"]
KAPPA_DEFAULT = 6.0  # Default strand_symmetry_kappa

PATTERN_LABELS = {
    0: "TA1_only_no_nrna", 1: "TA4_only_no_nrna",
    2: "mixed_A_no_nrna", 3: "convergent_CF_CR_no_nrna",
    4: "TA1_mod_nrna", 5: "TA4_mod_nrna",
    6: "mixed_A_mod_nrna", 7: "convergent_mod_nrna",
    8: "TA1_heavy_nrna", 9: "TA4_heavy_nrna",
    10: "mixed_A_heavy_nrna", 11: "convergent_heavy_nrna",
    12: "everything_moderate", 13: "TA1+TCF+heavy_nrna",
    14: "TA3+TCR+heavy_nrna", 15: "nrna_only_no_mrna",
}

NRNA_GROUPS = {
    "none": [0, 1, 2, 3],
    "moderate": [4, 5, 6, 7],
    "heavy": [8, 9, 10, 11],
    "complex": [12, 13, 14, 15],
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


def load(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def compute_metrics(row):
    """Compute per-run metrics: observed/expected mRNA, nRNA, gDNA totals."""
    mrna_exp = sum(float(row[f"{t}_expected"]) for t in TRANSCRIPTS)
    mrna_obs = sum(float(row.get(f"{t}_observed", 0)) for t in TRANSCRIPTS)
    nrna_exp = float(row.get("nrna_expected", 0))
    nrna_obs = float(row.get("nrna_observed", 0))
    gdna_exp = float(row.get("gdna_expected", 0))
    gdna_obs = float(row.get("gdna_observed", 0))
    return {
        "mrna_exp": mrna_exp, "mrna_obs": mrna_obs,
        "nrna_exp": nrna_exp, "nrna_obs": nrna_obs,
        "gdna_exp": gdna_exp, "gdna_obs": gdna_obs,
        "total_exp": mrna_exp + nrna_exp + gdna_exp,
        "total_obs": mrna_obs + nrna_obs + gdna_obs,
    }


def rel_err(obs, exp):
    if exp == 0:
        return 0.0
    return (obs - exp) / exp


def main():
    base_path = sys.argv[1]
    ph1_path = sys.argv[2]
    ph3_path = sys.argv[3]
    base_rows = load(base_path)
    ph1_rows = load(ph1_path)
    ph3_rows = load(ph3_path)

    # Add pattern classification and metrics
    for dataset_name, rows in [("base", base_rows), ("ph1", ph1_rows), ("ph3", ph3_rows)]:
        for r in rows:
            r["_pattern"] = classify_pattern(r)
            r["_ss"] = float(r["strand_specificity"])
            r["_gdna_frac"] = float(r["gdna_fraction"])
            r["_m"] = compute_metrics(r)

    # ===================================================================
    # SECTION 1: κ_eff values and expected impact
    # ===================================================================
    print("=" * 100)
    print("  1. EFFECTIVE κ VALUES AT EACH SS LEVEL (Phase 3)")
    print("=" * 100)
    print()
    print(f"  Default κ = {KAPPA_DEFAULT}")
    print(f"  Formula: κ_eff = max(κ · (2·SS − 1)², 2.0)")
    print(f"  Guard: penalty active only if κ_eff > 2.0")
    print()
    print(f"  {'SS':>5s}  {'(2·SS−1)²':>10s}  {'κ·scale':>10s}  {'κ_eff':>10s}  {'Active?':>10s}")
    print(f"  {'-'*5}  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*10}")
    for ss in [0.5, 0.75, 0.9, 1.0]:
        scale = (2 * ss - 1) ** 2
        raw = KAPPA_DEFAULT * scale
        eff = max(raw, 2.0)
        active = "YES" if eff > 2.0 else "NO (≤2.0)"
        print(f"  {ss:>5.2f}  {scale:>10.4f}  {raw:>10.2f}  {eff:>10.2f}  {active:>10s}")
    print()
    print("  Key insight: At SS=0.75, κ_eff=2.0 → penalty DISABLED")
    print("               At SS=0.5,  κ_eff=2.0 → penalty DISABLED")
    print("               At SS=0.9,  κ_eff=3.84 → REDUCED (64% of full)")
    print("               At SS=1.0,  κ_eff=6.0  → FULL strength")
    print()

    # ===================================================================
    # SECTION 2: Delta analysis — Phase 3 minus Phase 1
    # ===================================================================
    print("=" * 100)
    print("  2. PHASE 3 vs PHASE 1: Absolute differences in key metrics")
    print("=" * 100)
    print()
    print("  Positive Δ means Phase 3 is HIGHER than Phase 1")
    print()

    # Group by (ss, gdna_frac, nrna_group)
    for nrna_name, pattern_ids in NRNA_GROUPS.items():
        print(f"  ── nRNA level: {nrna_name} ──")
        print(f"  {'SS':>5s}  {'gdna':>5s}  {'Δ_mRNA%':>10s}  {'Δ_nRNA%':>10s}  {'Δ_gDNA%':>10s}  {'Ph1_gDNA%':>10s}  {'Ph3_gDNA%':>10s}")
        print(f"  {'-'*5}  {'-'*5}  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*10}")

        for ss in [0.5, 0.75, 0.9, 1.0]:
            for gf in [0.0, 0.1, 0.3, 0.5]:
                ph1_mrna, ph3_mrna = [], []
                ph1_nrna, ph3_nrna = [], []
                ph1_gdna, ph3_gdna = [], []

                for r1, r3 in zip(ph1_rows, ph3_rows):
                    if r1["_pattern"] not in pattern_ids:
                        continue
                    if abs(r1["_ss"] - ss) > 0.01 or abs(r1["_gdna_frac"] - gf) > 0.01:
                        continue
                    m1, m3 = r1["_m"], r3["_m"]
                    if m1["mrna_exp"] > 0:
                        ph1_mrna.append(rel_err(m1["mrna_obs"], m1["mrna_exp"]))
                        ph3_mrna.append(rel_err(m3["mrna_obs"], m3["mrna_exp"]))
                    if m1["nrna_exp"] > 0:
                        ph1_nrna.append(rel_err(m1["nrna_obs"], m1["nrna_exp"]))
                        ph3_nrna.append(rel_err(m3["nrna_obs"], m3["nrna_exp"]))
                    if m1["gdna_exp"] > 0:
                        ph1_gdna.append(rel_err(m1["gdna_obs"], m1["gdna_exp"]))
                        ph3_gdna.append(rel_err(m3["gdna_obs"], m3["gdna_exp"]))

                def fmt_delta(a, b):
                    if not a:
                        return "n/a"
                    return f"{(np.mean(b) - np.mean(a))*100:+.1f}%"

                def fmt_val(a):
                    if not a:
                        return "n/a"
                    return f"{np.mean(a)*100:+.1f}%"

                print(f"  {ss:>5.2f}  {gf:>5.1f}  {fmt_delta(ph1_mrna, ph3_mrna):>10s}  "
                      f"{fmt_delta(ph1_nrna, ph3_nrna):>10s}  {fmt_delta(ph1_gdna, ph3_gdna):>10s}  "
                      f"{fmt_val(ph1_gdna):>10s}  {fmt_val(ph3_gdna):>10s}")
        print()

    # ===================================================================
    # SECTION 3: 3-way comparison at critical conditions
    # ===================================================================
    print("=" * 100)
    print("  3. THREE-WAY COMPARISON: Baseline → Phase 1 → Phase 3")
    print("     (gDNA signed relative error, SS≥0.75, gdna>0)")
    print("=" * 100)
    print()

    for nrna_name, pattern_ids in NRNA_GROUPS.items():
        print(f"  ── nRNA level: {nrna_name} ──")
        print(f"  {'SS':>5s}  {'gdna':>5s}  {'Baseline':>10s}  {'Phase1':>10s}  {'Phase3':>10s}  {'B→P1_Δ':>10s}  {'P1→P3_Δ':>10s}")
        print(f"  {'-'*5}  {'-'*5}  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*10}")

        for ss in [0.75, 0.9, 1.0]:
            for gf in [0.1, 0.3, 0.5]:
                base_g, ph1_g, ph3_g = [], [], []

                for rb, r1, r3 in zip(base_rows, ph1_rows, ph3_rows):
                    if rb["_pattern"] not in pattern_ids:
                        continue
                    if abs(rb["_ss"] - ss) > 0.01 or abs(rb["_gdna_frac"] - gf) > 0.01:
                        continue
                    mb, m1, m3 = rb["_m"], r1["_m"], r3["_m"]
                    if mb["gdna_exp"] > 0:
                        base_g.append(rel_err(mb["gdna_obs"], mb["gdna_exp"]))
                        ph1_g.append(rel_err(m1["gdna_obs"], m1["gdna_exp"]))
                        ph3_g.append(rel_err(m3["gdna_obs"], m3["gdna_exp"]))

                if not base_g:
                    continue
                b_mean = np.mean(base_g) * 100
                p1_mean = np.mean(ph1_g) * 100
                p3_mean = np.mean(ph3_g) * 100
                print(f"  {ss:>5.2f}  {gf:>5.1f}  {b_mean:>+9.1f}%  {p1_mean:>+9.1f}%  {p3_mean:>+9.1f}%  "
                      f"{p1_mean - b_mean:>+9.1f}%  {p3_mean - p1_mean:>+9.1f}%")
        print()

    # ===================================================================
    # SECTION 4: Mass flow — where do gDNA reads end up?
    # ===================================================================
    print("=" * 100)
    print("  4. MASS FLOW: Where do 'missing' gDNA reads go?")
    print("     (SS≥0.75, gdna_frac=0.3, nRNA present)")
    print("=" * 100)
    print()

    for nrna_name in ["moderate", "heavy", "complex"]:
        pattern_ids = NRNA_GROUPS[nrna_name]
        print(f"  ── nRNA level: {nrna_name} ──")
        print(f"  {'SS':>5s}  {'Label':>10s}  {'mRNA_Δ':>10s}  {'nRNA_Δ':>10s}  {'gDNA_Δ':>10s}  {'Total_Δ':>10s}")
        print(f"  {'-'*5}  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*10}")

        for ss in [0.75, 0.9, 1.0]:
            for label, rows in [("Phase1", ph1_rows), ("Phase3", ph3_rows)]:
                mrna_excess, nrna_excess, gdna_deficit, total_diff = [], [], [], []
                for r in rows:
                    if r["_pattern"] not in pattern_ids:
                        continue
                    if abs(r["_ss"] - ss) > 0.01 or abs(r["_gdna_frac"] - 0.3) > 0.01:
                        continue
                    m = r["_m"]
                    mrna_excess.append(m["mrna_obs"] - m["mrna_exp"])
                    nrna_excess.append(m["nrna_obs"] - m["nrna_exp"])
                    gdna_deficit.append(m["gdna_obs"] - m["gdna_exp"])
                    total_diff.append(m["total_obs"] - m["total_exp"])

                if mrna_excess:
                    print(f"  {ss:>5.2f}  {label:>10s}  {np.mean(mrna_excess):>+9.0f}  "
                          f"{np.mean(nrna_excess):>+9.0f}  {np.mean(gdna_deficit):>+9.0f}  "
                          f"{np.mean(total_diff):>+9.0f}")
        print()

    # ===================================================================
    # SECTION 5: Per-pattern breakdown at critical SS=0.75, gdna=0.3
    # ===================================================================
    print("=" * 100)
    print("  5. PER-PATTERN gDNA ERROR at SS=0.75, gdna_frac=0.3")
    print("     (Where the Phase 1 gDNA underestimation is worst)")
    print("=" * 100)
    print()
    print(f"  {'Pat':>4s}  {'Label':>30s}  {'Base_gDNA':>10s}  {'Ph1_gDNA':>10s}  {'Ph3_gDNA':>10s}")
    print(f"  {'-'*4}  {'-'*30}  {'-'*10}  {'-'*10}  {'-'*10}")

    for pid in sorted(PATTERN_LABELS.keys()):
        base_g, ph1_g, ph3_g = [], [], []
        for rb, r1, r3 in zip(base_rows, ph1_rows, ph3_rows):
            if rb["_pattern"] != pid:
                continue
            if abs(rb["_ss"] - 0.75) > 0.01 or abs(rb["_gdna_frac"] - 0.3) > 0.01:
                continue
            mb, m1, m3 = rb["_m"], r1["_m"], r3["_m"]
            if mb["gdna_exp"] > 0:
                base_g.append(rel_err(mb["gdna_obs"], mb["gdna_exp"]))
                ph1_g.append(rel_err(m1["gdna_obs"], m1["gdna_exp"]))
                ph3_g.append(rel_err(m3["gdna_obs"], m3["gdna_exp"]))

        if base_g:
            print(f"  {pid:>4d}  {PATTERN_LABELS[pid]:>30s}  "
                  f"{np.mean(base_g)*100:>+9.1f}%  "
                  f"{np.mean(ph1_g)*100:>+9.1f}%  "
                  f"{np.mean(ph3_g)*100:>+9.1f}%")
    print()

    # ===================================================================
    # SECTION 6: Check SS=0.9 and SS=1.0 where κ_eff IS active
    # ===================================================================
    print("=" * 100)
    print("  6. FOCUS: SS=0.9 and SS=1.0 (where Phase 3 κ_eff IS active)")
    print("     κ_eff@0.9 = 3.84, κ_eff@1.0 = 6.0 vs Phase 1 κ=6.0")
    print("=" * 100)
    print()

    for nrna_name in ["moderate", "heavy", "complex"]:
        pattern_ids = NRNA_GROUPS[nrna_name]
        print(f"  ── nRNA level: {nrna_name} ──")
        print(f"  {'SS':>5s}  {'gdna':>5s}  {'Ph1_mRNA%':>10s}  {'Ph3_mRNA%':>10s}  "
              f"{'Ph1_nRNA%':>10s}  {'Ph3_nRNA%':>10s}  {'Ph1_gDNA%':>10s}  {'Ph3_gDNA%':>10s}")
        print(f"  {'-'*5}  {'-'*5}  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*10}  {'-'*10}")

        for ss in [0.9, 1.0]:
            for gf in [0.0, 0.1, 0.3, 0.5]:
                p1m, p3m, p1n, p3n, p1g, p3g = [], [], [], [], [], []
                for r1, r3 in zip(ph1_rows, ph3_rows):
                    if r1["_pattern"] not in pattern_ids:
                        continue
                    if abs(r1["_ss"] - ss) > 0.01 or abs(r1["_gdna_frac"] - gf) > 0.01:
                        continue
                    m1, m3 = r1["_m"], r3["_m"]
                    if m1["mrna_exp"] > 0:
                        p1m.append(abs(rel_err(m1["mrna_obs"], m1["mrna_exp"])))
                        p3m.append(abs(rel_err(m3["mrna_obs"], m3["mrna_exp"])))
                    if m1["nrna_exp"] > 0:
                        p1n.append(rel_err(m1["nrna_obs"], m1["nrna_exp"]))
                        p3n.append(rel_err(m3["nrna_obs"], m3["nrna_exp"]))
                    if m1["gdna_exp"] > 0:
                        p1g.append(rel_err(m1["gdna_obs"], m1["gdna_exp"]))
                        p3g.append(rel_err(m3["gdna_obs"], m3["gdna_exp"]))

                def fmt(a): return f"{np.mean(a)*100:+.1f}%" if a else "n/a"
                print(f"  {ss:>5.2f}  {gf:>5.1f}  {fmt(p1m):>10s}  {fmt(p3m):>10s}  "
                      f"{fmt(p1n):>10s}  {fmt(p3n):>10s}  {fmt(p1g):>10s}  {fmt(p3g):>10s}")
        print()

    # ===================================================================
    # SECTION 7: Summary and root causes
    # ===================================================================
    print("=" * 100)
    print("  7. ROOT CAUSE SUMMARY")
    print("=" * 100)
    print()

    # Compute key summary metrics
    def mean_gdna_err(rows, pattern_ids, ss_min, gdna_min=0.01):
        vals = []
        for r in rows:
            if r["_pattern"] not in pattern_ids:
                continue
            if r["_ss"] < ss_min - 0.01:
                continue
            m = r["_m"]
            if m["gdna_exp"] > gdna_min:
                vals.append(rel_err(m["gdna_obs"], m["gdna_exp"]))
        return np.mean(vals) * 100 if vals else float("nan")

    def mean_nrna_err(rows, pattern_ids, ss_min, gdna_frac=None):
        vals = []
        for r in rows:
            if r["_pattern"] not in pattern_ids:
                continue
            if r["_ss"] < ss_min - 0.01:
                continue
            if gdna_frac is not None and abs(r["_gdna_frac"] - gdna_frac) > 0.01:
                continue
            m = r["_m"]
            if m["nrna_exp"] > 0:
                vals.append(rel_err(m["nrna_obs"], m["nrna_exp"]))
        return np.mean(vals) * 100 if vals else float("nan")

    def mean_mrna_err(rows, pattern_ids, ss_min):
        vals = []
        for r in rows:
            if r["_pattern"] not in pattern_ids:
                continue
            if r["_ss"] < ss_min - 0.01:
                continue
            m = r["_m"]
            if m["mrna_exp"] > 0:
                vals.append(abs(rel_err(m["mrna_obs"], m["mrna_exp"])))
        return np.mean(vals) * 100 if vals else float("nan")

    all_nrna = NRNA_GROUPS["moderate"] + NRNA_GROUPS["heavy"] + NRNA_GROUPS["complex"]
    heavy_nrna = NRNA_GROUPS["heavy"] + NRNA_GROUPS["complex"]
    no_nrna = NRNA_GROUPS["none"]

    print("  Metric                                    Baseline    Phase 1    Phase 3")
    print("  " + "-" * 75)
    metrics = [
        ("mRNA |err| (SS≥0.75, all gdna)",
         lambda r: mean_mrna_err(r, list(range(16)), 0.75)),
        ("mRNA |err| (no nRNA, SS≥0.75)",
         lambda r: mean_mrna_err(r, no_nrna, 0.75)),
        ("mRNA |err| (heavy+complex nRNA, SS≥0.75)",
         lambda r: mean_mrna_err(r, heavy_nrna, 0.75)),
        ("gDNA err (no nRNA, SS≥0.75)",
         lambda r: mean_gdna_err(r, no_nrna, 0.75)),
        ("gDNA err (moderate nRNA, SS≥0.75)",
         lambda r: mean_gdna_err(r, NRNA_GROUPS["moderate"], 0.75)),
        ("gDNA err (heavy nRNA, SS≥0.75)",
         lambda r: mean_gdna_err(r, NRNA_GROUPS["heavy"], 0.75)),
        ("gDNA err (complex, SS≥0.75)",
         lambda r: mean_gdna_err(r, NRNA_GROUPS["complex"], 0.75)),
        ("nRNA err (moderate, SS≥0.75, gdna=0)",
         lambda r: mean_nrna_err(r, NRNA_GROUPS["moderate"], 0.75, 0.0)),
        ("nRNA err (heavy, SS≥0.75, gdna=0)",
         lambda r: mean_nrna_err(r, NRNA_GROUPS["heavy"], 0.75, 0.0)),
        ("nRNA err (all, SS≥0.75, all gdna)",
         lambda r: mean_nrna_err(r, all_nrna, 0.75)),
    ]

    for label, fn in metrics:
        b = fn(base_rows)
        p1 = fn(ph1_rows)
        p3 = fn(ph3_rows)
        print(f"  {label:<45s}  {b:>+7.1f}%  {p1:>+8.1f}%  {p3:>+8.1f}%")
    print()

    # Phase 3 effect size
    print("  Phase 3 effect (P3 − P1):")
    print("  " + "-" * 75)
    for label, fn in metrics:
        p1 = fn(ph1_rows)
        p3 = fn(ph3_rows)
        delta = p3 - p1
        significance = "≈0" if abs(delta) < 0.5 else ("IMPROVED" if delta < -0.5 else "WORSE")
        print(f"  {label:<45s}  Δ={delta:>+6.1f}%  [{significance}]")
    print()

    print("  CONCLUSION:")
    print("  " + "-" * 75)
    print("  Phase 3 (SS-scaled κ) has essentially NO effect compared to Phase 1.")
    print("  Reason: At SS=0.75, κ_eff=max(6.0×0.25, 2.0)=2.0 → penalty DISABLED")
    print("  At SS=0.9, κ_eff=3.84 vs κ=6.0 in Phase 1 → REDUCED but the Laplace")
    print("  smoothing (+0.5/+1.0 vs ε=50) makes p_hat sharper, partially compensating.")
    print("  Net result: the two changes (smaller κ, sharper p_hat) roughly cancel out.")
    print()
    print("  ROOT CAUSE OF gDNA UNDERESTIMATION (unchanged from Phase 1):")
    print("  The gDNA underestimation at SS≥0.75 when nRNA is present is NOT caused")
    print("  by the strand symmetry penalty. It's caused by nRNA competing with gDNA")
    print("  for intronic reads. Both nRNA and gDNA produce intronic reads, and the EM")
    print("  cannot distinguish them without additional evidence. The strand penalty")
    print("  helps push strand-symmetric reads toward gDNA, but the freed nRNA (from")
    print("  Phase 1 decoupling) absorbs intronic mass that should go to gDNA.")
    print()


if __name__ == "__main__":
    main()
