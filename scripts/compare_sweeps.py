#!/usr/bin/env python3
"""Compare two sweep_results.tsv files side-by-side.

Usage:
    python scripts/compare_sweeps.py <baseline.tsv> <phase1.tsv>
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


def nrna_level(pattern):
    if pattern in (0, 1, 2, 3):
        return "none"
    elif pattern in (4, 5, 6, 7):
        return "moderate"
    elif pattern in (8, 9, 10, 11):
        return "heavy"
    elif pattern in (12, 13, 14, 15):
        return "complex"
    return "unknown"


def load_data(tsv_path):
    with open(tsv_path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def compute_metrics(rows):
    """Compute per-cell metrics for the sweep grid."""
    # Group by (nrna_level, ss, gdna_frac)
    mrna_err = defaultdict(list)
    nrna_err = defaultdict(list)
    gdna_err = defaultdict(list)
    pattern15 = []

    for row in rows:
        ss = float(row["strand_specificity"])
        gf = float(row["gdna_fraction"])
        pat = classify_pattern(row)
        nl = nrna_level(pat)

        # Total mRNA error
        total_exp = float(row.get("total_mrna_expected", 0))
        total_obs = float(row.get("total_mrna_observed", 0))
        if total_exp > 0:
            rel = abs(total_obs - total_exp) / total_exp
            mrna_err[(nl, ss, gf)].append(rel)

        # nRNA error
        nrna_exp = float(row.get("nrna_expected", 0))
        nrna_obs = float(row.get("nrna_observed", 0))
        if nrna_exp > 0:
            nrna_rel = (nrna_obs - nrna_exp) / nrna_exp
            nrna_err[(nl, ss, gf)].append(nrna_rel)

        # gDNA error
        gdna_exp = float(row.get("gdna_expected", 0))
        gdna_obs = float(row.get("gdna_observed", 0))
        if gdna_exp > 0:
            gdna_rel = (gdna_obs - gdna_exp) / gdna_exp
            gdna_err[(nl, ss, gf)].append(gdna_rel)

        # Pattern 15
        if pat == 15:
            pattern15.append({
                "ss": ss, "gf": gf,
                "mrna_exp": total_exp, "mrna_obs": total_obs,
                "nrna_exp": nrna_exp, "nrna_obs": nrna_obs,
                "gdna_exp": gdna_exp, "gdna_obs": gdna_obs,
            })

    return mrna_err, nrna_err, gdna_err, pattern15


def fmt_pct(val):
    return f"{val:+.1f}%" if val is not None else "   n/a"


def mean_or_none(lst):
    return np.mean(lst) if lst else None


SS_VALS = [0.5, 0.75, 0.9, 1.0]
GF_VALS = [0.0, 0.1, 0.3, 0.5]
NRNA_LEVELS = ["none", "moderate", "heavy", "complex"]


def print_comparison_table(title, base_data, ph1_data, signed=False):
    print(f"\n{'='*90}")
    print(f"  {title}")
    print(f"{'='*90}")
    for nl in NRNA_LEVELS:
        print(f"\n  nRNA level: {nl}")
        print(f"  {'gdna_f':>8}  ", end="")
        for ss in SS_VALS:
            print(f"  {'ss='+str(ss):^24}", end="")
        print()
        print(f"  {'':>8}  ", end="")
        for ss in SS_VALS:
            print(f"  {'base':>10}  {'ph1':>10}", end="")
        print()
        print("  " + "-" * 106)
        for gf in GF_VALS:
            print(f"  {gf:>8.1f}  ", end="")
            for ss in SS_VALS:
                key = (nl, ss, gf)
                bv = base_data.get(key, [])
                pv = ph1_data.get(key, [])
                if signed:
                    b = mean_or_none(bv)
                    p = mean_or_none(pv)
                else:
                    b = np.mean([abs(x) for x in bv]) if bv else None
                    p = np.mean([abs(x) for x in pv]) if pv else None
                bs = f"{b*100:+.1f}%" if b is not None else "   n/a"
                ps = f"{p*100:+.1f}%" if p is not None else "   n/a"
                print(f"  {bs:>10}  {ps:>10}", end="")
            print()


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <baseline.tsv> <phase1.tsv>")
        sys.exit(1)

    base_rows = load_data(sys.argv[1])
    ph1_rows = load_data(sys.argv[2])

    print(f"Baseline: {len(base_rows)} runs")
    print(f"Phase 1:  {len(ph1_rows)} runs")

    b_mrna, b_nrna, b_gdna, b_p15 = compute_metrics(base_rows)
    p_mrna, p_nrna, p_gdna, p_p15 = compute_metrics(ph1_rows)

    print_comparison_table(
        "TOTAL mRNA ACCURACY (mean |relative error|)",
        b_mrna, p_mrna, signed=False)

    print_comparison_table(
        "nRNA ACCURACY (mean signed relative error)",
        b_nrna, p_nrna, signed=True)

    print_comparison_table(
        "gDNA ACCURACY (mean signed relative error)",
        b_gdna, p_gdna, signed=True)

    # Pattern 15 summary
    print(f"\n{'='*90}")
    print("  PATTERN 15 STRESS TEST: nRNA only (no mRNA expressed)")
    print(f"{'='*90}")
    print(f"\n  {'SS':>4} {'gdna':>5}  {'---BASELINE---':>32}  {'---PHASE 1---':>32}")
    print(f"  {'':>4} {'frac':>5}  {'mrna_obs':>8} {'nrna_rel%':>10} {'gdna_rel%':>10}  {'mrna_obs':>8} {'nrna_rel%':>10} {'gdna_rel%':>10}")
    print("  " + "-" * 85)

    # Group pattern 15 by (ss, gf)
    b15 = defaultdict(list)
    p15 = defaultdict(list)
    for d in b_p15:
        b15[(d["ss"], d["gf"])].append(d)
    for d in p_p15:
        p15[(d["ss"], d["gf"])].append(d)

    for ss in SS_VALS:
        for gf in GF_VALS:
            bd = b15.get((ss, gf), [])
            pd = p15.get((ss, gf), [])
            if not bd and not pd:
                continue
            # Average
            b_mrna_obs = np.mean([d["mrna_obs"] for d in bd]) if bd else 0
            b_nrna_rel = np.mean([(d["nrna_obs"] - d["nrna_exp"]) / d["nrna_exp"] for d in bd]) if bd else 0
            b_gdna_rel = np.mean([(d["gdna_obs"] - d["gdna_exp"]) / d["gdna_exp"] for d in bd if d["gdna_exp"] > 0]) if any(d["gdna_exp"] > 0 for d in bd) else None

            p_mrna_obs = np.mean([d["mrna_obs"] for d in pd]) if pd else 0
            p_nrna_rel = np.mean([(d["nrna_obs"] - d["nrna_exp"]) / d["nrna_exp"] for d in pd]) if pd else 0
            p_gdna_rel = np.mean([(d["gdna_obs"] - d["gdna_exp"]) / d["gdna_exp"] for d in pd if d["gdna_exp"] > 0]) if any(d["gdna_exp"] > 0 for d in pd) else None

            bg = f"{b_gdna_rel*100:+.1f}%" if b_gdna_rel is not None else "     n/a"
            pg = f"{p_gdna_rel*100:+.1f}%" if p_gdna_rel is not None else "     n/a"

            print(f"  {ss:>4} {gf:>5.1f}  {b_mrna_obs:>8.1f} {b_nrna_rel*100:>+10.1f}% {bg:>10}  {p_mrna_obs:>8.1f} {p_nrna_rel*100:>+10.1f}% {pg:>10}")

    # Summary deltas
    print(f"\n{'='*90}")
    print("  SUMMARY: Key metric changes (Phase 1 vs Baseline)")
    print(f"{'='*90}")

    comparisons = [
        ("mRNA accuracy (no nRNA, gdna=0)", "none", 0.0, False),
        ("mRNA accuracy (moderate nRNA, SS≥0.75, gdna=0)", "moderate", 0.0, False),
        ("mRNA accuracy (heavy nRNA, SS≥0.75, gdna=0)", "heavy", 0.0, False),
        ("mRNA accuracy (complex, SS≥0.75, gdna=0)", "complex", 0.0, False),
    ]

    for label, nl, gf, signed in comparisons:
        b_vals = []
        p_vals = []
        for ss in [0.75, 0.9, 1.0]:
            b_vals.extend(b_mrna.get((nl, ss, gf), []))
            p_vals.extend(p_mrna.get((nl, ss, gf), []))
        if nl == "none":
            for ss in SS_VALS:
                b_vals.extend(b_mrna.get((nl, ss, gf), []))
                p_vals.extend(p_mrna.get((nl, ss, gf), []))
        bm = np.mean([abs(x) for x in b_vals]) * 100 if b_vals else 0
        pm = np.mean([abs(x) for x in p_vals]) * 100 if p_vals else 0
        delta = pm - bm
        print(f"  {label:55s}  base={bm:5.1f}%  ph1={pm:5.1f}%  Δ={delta:+.1f}%")

    # nRNA comparisons
    nrna_comparisons = [
        ("nRNA accuracy (moderate, SS≥0.75, gdna=0)", "moderate", [0.75, 0.9, 1.0], [0.0]),
        ("nRNA accuracy (moderate, SS≥0.75, all gdna)", "moderate", [0.75, 0.9, 1.0], GF_VALS),
        ("nRNA accuracy (heavy, SS≥0.75, gdna=0)", "heavy", [0.75, 0.9, 1.0], [0.0]),
        ("nRNA accuracy (all, SS=0.5)", "ALL", [0.5], GF_VALS),
    ]
    print()
    for label, nl, ss_list, gf_list in nrna_comparisons:
        b_vals = []
        p_vals = []
        nls = NRNA_LEVELS[1:] if nl == "ALL" else [nl]
        for n in nls:
            for ss in ss_list:
                for gf in gf_list:
                    b_vals.extend(b_nrna.get((n, ss, gf), []))
                    p_vals.extend(p_nrna.get((n, ss, gf), []))
        bm = np.mean(b_vals) * 100 if b_vals else 0
        pm = np.mean(p_vals) * 100 if p_vals else 0
        print(f"  {label:55s}  base={bm:+6.1f}%  ph1={pm:+6.1f}%")

    # gDNA worst cases
    print()
    gdna_comparisons = [
        ("gDNA accuracy (no nRNA, SS≥0.75)", "none", [0.75, 0.9, 1.0], [0.1, 0.3, 0.5]),
        ("gDNA accuracy (moderate nRNA, SS≥0.75)", "moderate", [0.75, 0.9, 1.0], [0.1, 0.3, 0.5]),
        ("gDNA accuracy (all nRNA, SS=0.5)", "ALL", [0.5], [0.1, 0.3, 0.5]),
    ]
    for label, nl, ss_list, gf_list in gdna_comparisons:
        b_vals = []
        p_vals = []
        nls = NRNA_LEVELS if nl == "ALL" else [nl]
        for n in nls:
            for ss in ss_list:
                for gf in gf_list:
                    b_vals.extend(b_gdna.get((n, ss, gf), []))
                    p_vals.extend(p_gdna.get((n, ss, gf), []))
        bm = np.mean(b_vals) * 100 if b_vals else 0
        pm = np.mean(p_vals) * 100 if p_vals else 0
        print(f"  {label:55s}  base={bm:+6.1f}%  ph1={pm:+6.1f}%")


if __name__ == "__main__":
    main()
