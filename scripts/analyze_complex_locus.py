#!/usr/bin/env python3
"""Analyze complex locus ablation sweep results.

Produces diagnostic tables for:
  1. mRNA accuracy by scenario × SS × gdna_fraction
  2. nRNA accuracy (expected vs observed)
  3. gDNA accuracy (expected vs observed)
  4. Per-transcript error patterns
  5. Worst-case failure modes
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

TSV = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(
    "/Users/mkiyer/Downloads/rigel_runs/complex_locus/sweep_results.tsv")

TRANSCRIPTS = ["TA1", "TA2", "TA3", "TA4", "TB1", "TB2", "TCF", "TCR", "TD", "TE"]
NRNA_LABELS = ["NTA1", "NTA2", "NTA34", "NTB12", "NTCR", "NTE"]

# Pattern labels for readability
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
    """Classify a row into a pattern index based on mRNA/nRNA values."""
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
    if ta1 == 1000 and tcf == 1000 and nta1 == 1000:
        return 13
    if ta3 == 1000 and tcr == 1000 and nta34 == 1000:
        return 14

    # nRNA level
    nrna_level = "none"
    if total_nrna > 0 and total_nrna <= 1000:
        nrna_level = "mod"
    elif total_nrna > 1000:
        nrna_level = "heavy"

    # Expression pattern
    if ta1 > 0 and ta2 == 0 and ta3 == 0 and ta4 == 0 and tcf == 0:
        base = {("none", 0): 0, ("mod", 4): 4, ("heavy", 8): 8}
        return base.get((nrna_level, {0: 0, 4: 4, 8: 8}.get(
            0 if nrna_level == "none" else 4 if nrna_level == "mod" else 8, 0)), 0)
    if ta4 > 0 and ta1 == 0 and ta2 == 0 and ta3 == 0 and tcf == 0:
        return {"none": 1, "mod": 5, "heavy": 9}[nrna_level]
    if ta1 > 0 and ta2 > 0 and tcf == 0:
        return {"none": 2, "mod": 6, "heavy": 10}[nrna_level]
    if tcf > 0 and tcr > 0 and ta1 == 0:
        return {"none": 3, "mod": 7, "heavy": 11}[nrna_level]

    return -1  # unknown


def load_data(tsv_path):
    with open(tsv_path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def safe_float(v, default=0.0):
    try:
        return float(v) if v else default
    except (ValueError, TypeError):
        return default


def section(title):
    w = 100
    print(f"\n{'=' * w}")
    print(f"  {title}")
    print(f"{'=' * w}")


def subsection(title):
    print(f"\n{'─' * 90}")
    print(f"  {title}")
    print(f"{'─' * 90}")


def analyze_mrna_accuracy(rows):
    """mRNA total accuracy by nRNA_level × SS × gdna_fraction."""
    section("1. TOTAL mRNA ACCURACY (relative error)")

    # Group by (nrna_level, ss, gdna_frac)
    # nrna_level: none (patterns 0-3), moderate (4-7), heavy (8-11), mixed (12-15)
    def nrna_level(r):
        pidx = classify_pattern(r)
        if pidx <= 3:
            return "none"
        elif pidx <= 7:
            return "moderate"
        elif pidx <= 11:
            return "heavy"
        else:
            return "complex"

    groups = defaultdict(list)
    for r in rows:
        nl = nrna_level(r)
        ss = r["strand_specificity"]
        gf = r["gdna_fraction"]
        rel = safe_float(r["total_mrna_rel_err"])
        groups[(nl, ss, gf)].append(rel)

    for nl in ["none", "moderate", "heavy", "complex"]:
        subsection(f"nRNA level: {nl}")
        ss_vals = sorted(set(r["strand_specificity"] for r in rows), key=float)
        gf_vals = sorted(set(r["gdna_fraction"] for r in rows), key=float)

        header = f"{'gdna_frac':>10}"
        for ss in ss_vals:
            header += f"  {'ss='+ss:>12}"
        print(header)
        print("-" * len(header))

        for gf in gf_vals:
            line = f"{gf:>10}"
            for ss in ss_vals:
                vals = groups.get((nl, ss, gf), [])
                if vals:
                    mean_err = sum(vals) / len(vals)
                    max_err = max(vals)
                    line += f"  {mean_err*100:>5.1f}%/{max_err*100:>4.0f}%"
                else:
                    line += f"  {'n/a':>12}"
            print(line)
        print(f"  (format: mean%/max%)")


def analyze_nrna_accuracy(rows):
    """nRNA accuracy by scenario."""
    section("2. NASCENT RNA (nRNA) ACCURACY")

    subsection("nRNA abs error by nRNA_level × SS × gdna_fraction")
    # Only rows where nRNA was expected
    nrna_rows = [r for r in rows if safe_float(r["nrna_expected"]) > 0]

    def nrna_level(r):
        pidx = classify_pattern(r)
        if pidx <= 3:
            return "none"
        elif pidx <= 7:
            return "moderate"
        elif pidx <= 11:
            return "heavy"
        else:
            return "complex"

    groups = defaultdict(list)
    for r in nrna_rows:
        nl = nrna_level(r)
        ss = r["strand_specificity"]
        gf = r["gdna_fraction"]
        exp = safe_float(r["nrna_expected"])
        obs = safe_float(r["nrna_observed"])
        # Signed error (positive = overestimate)
        err = obs - exp
        rel = err / exp if exp > 0 else 0
        groups[(nl, ss, gf)].append((err, rel, exp, obs))

    for nl in ["moderate", "heavy", "complex"]:
        subsection(f"nRNA level: {nl}")
        ss_vals = sorted(set(r["strand_specificity"] for r in nrna_rows), key=float)
        gf_vals = sorted(set(r["gdna_fraction"] for r in nrna_rows), key=float)

        header = f"{'gdna_frac':>10}"
        for ss in ss_vals:
            header += f"  {'ss='+ss:>14}"
        print(header)
        print("-" * len(header))

        for gf in gf_vals:
            line = f"{gf:>10}"
            for ss in ss_vals:
                data = groups.get((nl, ss, gf), [])
                if data:
                    mean_rel = sum(d[1] for d in data) / len(data)
                    mean_abs = sum(d[0] for d in data) / len(data)
                    line += f"  {mean_rel*100:>+6.1f}% ({mean_abs:>+.0f})"
                else:
                    line += f"  {'n/a':>14}"
            print(line)
        print(f"  (format: mean_rel% (mean_abs_err))")


def analyze_gdna_accuracy(rows):
    """gDNA accuracy by scenario."""
    section("3. GENOMIC DNA (gDNA) ACCURACY")

    # Only rows where gDNA was expected
    gdna_rows = [r for r in rows if safe_float(r["gdna_expected"]) > 0]

    if not gdna_rows:
        print("  No runs with gDNA expected > 0")
        return

    def nrna_level(r):
        pidx = classify_pattern(r)
        if pidx <= 3:
            return "none"
        elif pidx <= 7:
            return "moderate"
        elif pidx <= 11:
            return "heavy"
        else:
            return "complex"

    subsection("gDNA signed error by nRNA_level × SS × gdna_fraction")
    groups = defaultdict(list)
    for r in gdna_rows:
        nl = nrna_level(r)
        ss = r["strand_specificity"]
        gf = r["gdna_fraction"]
        exp = safe_float(r["gdna_expected"])
        obs = safe_float(r["gdna_observed"])
        err = obs - exp
        rel = err / exp if exp > 0 else 0
        groups[(nl, ss, gf)].append((err, rel, exp, obs))

    for nl in ["none", "moderate", "heavy", "complex"]:
        subsection(f"gDNA | nRNA level: {nl}")
        ss_vals = sorted(set(r["strand_specificity"] for r in gdna_rows), key=float)
        gf_vals = sorted(set(r["gdna_fraction"] for r in gdna_rows), key=float)

        header = f"{'gdna_frac':>10}  {'n':>4}"
        for ss in ss_vals:
            header += f"  {'ss='+ss:>16}"
        print(header)
        print("-" * len(header))

        for gf in gf_vals:
            for ss in ss_vals:
                data = groups.get((nl, ss, gf), [])
                if data:
                    break
            else:
                continue

            line = f"{gf:>10}"
            n_total = 0
            for ss in ss_vals:
                data = groups.get((nl, ss, gf), [])
                n_total = max(n_total, len(data))

            line = f"{gf:>10}  {n_total:>4}"
            for ss in ss_vals:
                data = groups.get((nl, ss, gf), [])
                if data:
                    mean_rel = sum(d[1] for d in data) / len(data)
                    mean_exp = sum(d[2] for d in data) / len(data)
                    mean_obs = sum(d[3] for d in data) / len(data)
                    line += f"  {mean_rel*100:>+6.1f}% {mean_obs:>5.0f}/{mean_exp:>5.0f}"
                else:
                    line += f"  {'n/a':>16}"
            print(line)
        print(f"  (format: mean_rel% obs/exp)")


def analyze_per_transcript(rows):
    """Per-transcript error analysis to find systematic biases."""
    section("4. PER-TRANSCRIPT ERROR ANALYSIS")

    # Group by transcript × gdna_fraction × SS
    for t_id in TRANSCRIPTS:
        # Only rows where this transcript was expected
        t_rows = [r for r in rows if safe_float(r[f"{t_id}_expected"]) > 0]
        if not t_rows:
            continue

        subsection(f"Transcript {t_id} (n={len(t_rows)} runs with expected > 0)")

        # By SS × gdna_frac
        groups = defaultdict(list)
        for r in t_rows:
            ss = r["strand_specificity"]
            gf = r["gdna_fraction"]
            exp = safe_float(r[f"{t_id}_expected"])
            obs = safe_float(r[f"{t_id}_observed"])
            signed_err = obs - exp
            rel = signed_err / exp if exp > 0 else 0
            groups[(ss, gf)].append((signed_err, rel, exp, obs))

        ss_vals = sorted(set(r["strand_specificity"] for r in t_rows), key=float)
        gf_vals = sorted(set(r["gdna_fraction"] for r in t_rows), key=float)

        header = f"{'gdna_frac':>10}"
        for ss in ss_vals:
            header += f"  {'ss='+ss:>14}"
        print(header)
        print("-" * len(header))

        for gf in gf_vals:
            line = f"{gf:>10}"
            for ss in ss_vals:
                data = groups.get((ss, gf), [])
                if data:
                    mean_rel = sum(d[1] for d in data) / len(data)
                    mean_signed = sum(d[0] for d in data) / len(data)
                    line += f"  {mean_rel*100:>+6.1f}%({mean_signed:>+5.0f})"
                else:
                    line += f"  {'n/a':>14}"
            print(line)
        print(f"  (format: mean_rel%(mean_signed_err))")


def analyze_false_positives(rows):
    """Find cases where unexpressed transcripts get non-zero estimates."""
    section("5. FALSE POSITIVES (unexpressed transcripts with non-zero estimates)")

    fp_cases = defaultdict(list)
    for r in rows:
        for t_id in TRANSCRIPTS:
            exp = safe_float(r[f"{t_id}_expected"])
            obs = safe_float(r[f"{t_id}_observed"])
            if exp == 0 and obs > 0.5:  # threshold for FP
                fp_cases[t_id].append({
                    "obs": obs,
                    "ss": r["strand_specificity"],
                    "gf": r["gdna_fraction"],
                    "pattern": classify_pattern(r),
                    "tb1": r["TB1"],
                    "tb2": r["TB2"],
                })

    if not fp_cases:
        print("  No false positives detected (all unexpressed transcripts estimated at 0)")
        return

    for t_id in TRANSCRIPTS:
        cases = fp_cases.get(t_id, [])
        if not cases:
            continue
        subsection(f"{t_id}: {len(cases)} false positive cases")

        # Summarize by SS × gdna_frac
        groups = defaultdict(list)
        for c in cases:
            groups[(c["ss"], c["gf"])].append(c["obs"])

        ss_vals = sorted(set(c["ss"] for c in cases), key=float)
        gf_vals = sorted(set(c["gf"] for c in cases), key=float)

        header = f"{'gdna_frac':>10}"
        for ss in ss_vals:
            header += f"  {'ss='+ss:>14}"
        print(header)
        print("-" * len(header))

        for gf in gf_vals:
            line = f"{gf:>10}"
            for ss in ss_vals:
                data = groups.get((ss, gf), [])
                if data:
                    mean_obs = sum(data) / len(data)
                    line += f"  {len(data):>3}× avg={mean_obs:>5.1f}"
                else:
                    line += f"  {'—':>14}"
            print(line)

        # What patterns generate FPs?
        pat_counts = defaultdict(int)
        for c in cases:
            pat_counts[c["pattern"]] += 1
        print(f"\n  Top patterns generating FPs for {t_id}:")
        for pidx, cnt in sorted(pat_counts.items(), key=lambda x: -x[1])[:5]:
            print(f"    pattern {pidx} ({PATTERN_LABELS.get(pidx, '?')}): {cnt} cases")


def analyze_worst_cases(rows):
    """Find the worst failures overall."""
    section("6. WORST CASES (top 20 by total_mrna_rel_err)")

    sorted_rows = sorted(rows, key=lambda r: -safe_float(r["total_mrna_rel_err"]))

    header = (f"{'pattern':>8} {'SS':>5} {'gdna_f':>7} {'TB1':>4} {'TB2':>4} "
              f"{'mrna_err%':>10} {'nrna_err':>9} {'gdna_err':>9} "
              f"{'mrna_exp':>9} {'mrna_obs':>9}")
    print(header)
    print("-" * len(header))

    for r in sorted_rows[:20]:
        pidx = classify_pattern(r)
        mrna_rel = safe_float(r["total_mrna_rel_err"]) * 100
        nrna_err = safe_float(r["nrna_observed"]) - safe_float(r["nrna_expected"])
        gdna_err = safe_float(r["gdna_observed"]) - safe_float(r["gdna_expected"])
        mrna_exp = safe_float(r["total_mrna_expected"])
        mrna_obs = safe_float(r["total_mrna_observed"])
        print(f"{pidx:>8} {r['strand_specificity']:>5} {r['gdna_fraction']:>7} "
              f"{r['TB1']:>4} {r['TB2']:>4} "
              f"{mrna_rel:>+9.1f}% {nrna_err:>+9.0f} {gdna_err:>+9.0f} "
              f"{mrna_exp:>9.0f} {mrna_obs:>9.1f}")


def analyze_nrna_only(rows):
    """Analyze the stress-test pattern: nRNA only, no mRNA."""
    section("7. STRESS TEST: nRNA ONLY (pattern 15 — no mRNA expressed)")

    p15_rows = [r for r in rows if classify_pattern(r) == 15]
    if not p15_rows:
        print("  No pattern-15 runs found")
        return

    header = (f"{'SS':>5} {'gdna_f':>7} {'TB1':>4} {'TB2':>4} "
              f"{'mrna_exp':>9} {'mrna_obs':>9} {'nrna_exp':>9} {'nrna_obs':>9} "
              f"{'gdna_exp':>9} {'gdna_obs':>9}")
    print(header)
    print("-" * len(header))

    for r in sorted(p15_rows, key=lambda r: (float(r["strand_specificity"]),
                                               float(r["gdna_fraction"]),
                                               float(r["TB1"]))):
        print(f"{r['strand_specificity']:>5} {r['gdna_fraction']:>7} "
              f"{r['TB1']:>4} {r['TB2']:>4} "
              f"{safe_float(r['total_mrna_expected']):>9.0f} "
              f"{safe_float(r['total_mrna_observed']):>9.1f} "
              f"{safe_float(r['nrna_expected']):>9.0f} "
              f"{safe_float(r['nrna_observed']):>9.1f} "
              f"{safe_float(r['gdna_expected']):>9.0f} "
              f"{safe_float(r['gdna_observed']):>9.1f}")


def analyze_convergent(rows):
    """Analyze convergent gene pairs (TCF/TCR) where exons overlap on opposite strands."""
    section("8. CONVERGENT GENE ACCURACY (TCF + strand, TCR - strand)")

    # Patterns 3, 7, 11 have only TCF+TCR active
    conv_rows = [r for r in rows if classify_pattern(r) in (3, 7, 11)]
    if not conv_rows:
        print("  No convergent-only runs found")
        return

    header = (f"{'pattern':>8} {'SS':>5} {'gdna_f':>7} {'TB1':>4} {'TB2':>4} "
              f"{'TCF_exp':>8} {'TCF_obs':>8} {'TCF_err%':>9} "
              f"{'TCR_exp':>8} {'TCR_obs':>8} {'TCR_err%':>9}")
    print(header)
    print("-" * len(header))

    for r in sorted(conv_rows, key=lambda r: (classify_pattern(r),
                                                float(r["strand_specificity"]),
                                                float(r["gdna_fraction"]))):
        pidx = classify_pattern(r)
        tcf_exp = safe_float(r["TCF_expected"])
        tcf_obs = safe_float(r["TCF_observed"])
        tcf_err = ((tcf_obs - tcf_exp) / tcf_exp * 100) if tcf_exp > 0 else 0
        tcr_exp = safe_float(r["TCR_expected"])
        tcr_obs = safe_float(r["TCR_observed"])
        tcr_err = ((tcr_obs - tcr_exp) / tcr_exp * 100) if tcr_exp > 0 else 0
        # Only show a sample (first 40)
        print(f"{pidx:>8} {r['strand_specificity']:>5} {r['gdna_fraction']:>7} "
              f"{r['TB1']:>4} {r['TB2']:>4} "
              f"{tcf_exp:>8.0f} {tcf_obs:>8.1f} {tcf_err:>+8.1f}% "
              f"{tcr_exp:>8.0f} {tcr_obs:>8.1f} {tcr_err:>+8.1f}%")


def analyze_gdna_fraction_effect(rows):
    """Quantify how much gDNA contamination biases mRNA estimates."""
    section("9. gDNA CONTAMINATION EFFECT ON mRNA")

    subsection("Mean signed mRNA error (obs - exp) by gdna_fraction × SS")

    groups = defaultdict(list)
    for r in rows:
        ss = r["strand_specificity"]
        gf = r["gdna_fraction"]
        exp = safe_float(r["total_mrna_expected"])
        obs = safe_float(r["total_mrna_observed"])
        if exp > 0:
            signed = obs - exp
            groups[(ss, gf)].append(signed)

    ss_vals = sorted(set(r["strand_specificity"] for r in rows), key=float)
    gf_vals = sorted(set(r["gdna_fraction"] for r in rows), key=float)

    header = f"{'gdna_frac':>10}"
    for ss in ss_vals:
        header += f"  {'ss='+ss:>12}"
    print(header)
    print("-" * len(header))

    for gf in gf_vals:
        line = f"{gf:>10}"
        for ss in ss_vals:
            data = groups.get((ss, gf), [])
            if data:
                mean_signed = sum(data) / len(data)
                line += f"  {mean_signed:>+12.1f}"
            else:
                line += f"  {'n/a':>12}"
        print(line)
    print(f"  (positive = mRNA overestimated due to gDNA leakage)")

    subsection("Where does gDNA mass go? (obs-exp for each component)")
    for gf in gf_vals:
        if gf == "0.0":
            continue
        print(f"\n  gdna_fraction = {gf}:")
        header2 = f"  {'SS':>5}  {'mrna_Δ':>10}  {'nrna_Δ':>10}  {'gdna_Δ':>10}  {'sum_Δ':>10}  {'gdna_exp':>10}"
        print(header2)
        for ss in ss_vals:
            mrna_delta = []
            nrna_delta = []
            gdna_delta = []
            gdna_exps = []
            for r in rows:
                if r["strand_specificity"] != ss or r["gdna_fraction"] != gf:
                    continue
                mrna_d = safe_float(r["total_mrna_observed"]) - safe_float(r["total_mrna_expected"])
                nrna_d = safe_float(r["nrna_observed"]) - safe_float(r["nrna_expected"])
                gdna_d = safe_float(r["gdna_observed"]) - safe_float(r["gdna_expected"])
                mrna_delta.append(mrna_d)
                nrna_delta.append(nrna_d)
                gdna_delta.append(gdna_d)
                gdna_exps.append(safe_float(r["gdna_expected"]))

            if mrna_delta:
                m = sum(mrna_delta) / len(mrna_delta)
                n = sum(nrna_delta) / len(nrna_delta)
                g = sum(gdna_delta) / len(gdna_delta)
                ge = sum(gdna_exps) / len(gdna_exps)
                print(f"  {ss:>5}  {m:>+10.1f}  {n:>+10.1f}  {g:>+10.1f}  {m+n+g:>+10.1f}  {ge:>10.0f}")


def analyze_tb_overlap(rows):
    """Analyze how TB1/TB2 (antisense overlap with gene A) affects accuracy."""
    section("10. ANTISENSE OVERLAP EFFECT (TB1/TB2 on Gene A accuracy)")

    # Compare runs with TB=0 vs TB>0 for patterns that have gene A active
    a_patterns = [0, 2, 4, 6, 8, 10, 12, 13]  # patterns with gene A active
    a_rows = [r for r in rows if classify_pattern(r) in a_patterns]

    groups = defaultdict(list)
    for r in a_rows:
        ss = r["strand_specificity"]
        gf = r["gdna_fraction"]
        tb_status = "TB_off" if (r["TB1"] == "0" and r["TB2"] == "0") else "TB_on"
        # Total A-locus error
        a_err = 0
        for t in ["TA1", "TA2", "TA3", "TA4"]:
            a_err += abs(safe_float(r[f"{t}_observed"]) - safe_float(r[f"{t}_expected"]))
        groups[(ss, gf, tb_status)].append(a_err)

    ss_vals = sorted(set(r["strand_specificity"] for r in a_rows), key=float)
    gf_vals = sorted(set(r["gdna_fraction"] for r in a_rows), key=float)

    header = f"{'gdna_frac':>10} {'TB':>6}"
    for ss in ss_vals:
        header += f"  {'ss='+ss:>10}"
    print(header)
    print("-" * len(header))

    for gf in gf_vals:
        for tb in ["TB_off", "TB_on"]:
            line = f"{gf:>10} {tb:>6}"
            for ss in ss_vals:
                data = groups.get((ss, gf, tb), [])
                if data:
                    mean = sum(data) / len(data)
                    line += f"  {mean:>10.1f}"
                else:
                    line += f"  {'n/a':>10}"
            print(line)
    print(f"  (values = mean total A-locus abs error)")


def main():
    rows = load_data(TSV)
    print(f"Loaded {len(rows)} runs from {TSV}")

    analyze_mrna_accuracy(rows)
    analyze_nrna_accuracy(rows)
    analyze_gdna_accuracy(rows)
    analyze_per_transcript(rows)
    analyze_false_positives(rows)
    analyze_worst_cases(rows)
    analyze_nrna_only(rows)
    analyze_convergent(rows)
    analyze_gdna_fraction_effect(rows)
    analyze_tb_overlap(rows)


if __name__ == "__main__":
    main()
