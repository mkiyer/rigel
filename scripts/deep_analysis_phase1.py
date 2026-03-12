#!/usr/bin/env python3
"""Deep-dive analysis of Phase 1 vs baseline sweep results.

Focuses on:
1. SS=0.5 mass distribution (the fundamental identifiability issue)
2. gDNA underestimation at SS≥0.75 (new issue in Phase 1)
3. nRNA overestimation at SS≥0.75 with gDNA (nRNA absorbs gDNA mass)
4. Pattern 15 mRNA false positives (stress test improvement)
5. Per-transcript false positive analysis
"""
import csv
import sys
from collections import defaultdict

import numpy as np

TRANSCRIPTS = ["TA1", "TA2", "TA3", "TA4", "TB1", "TB2", "TCF", "TCR", "TD", "TE"]
NRNA_LABELS = ["NTA1", "NTA2", "NTA34", "NTB12", "NTCR", "NTE"]

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


def main():
    base_path = sys.argv[1]
    ph1_path = sys.argv[2]
    base = load(base_path)
    ph1 = load(ph1_path)

    # ===================================================================
    # SECTION 1: Mass conservation analysis at SS=0.5
    # ===================================================================
    print("=" * 100)
    print("  1. MASS DISTRIBUTION AT SS=0.5 (fundamental identifiability issue)")
    print("=" * 100)
    print()
    print("At SS=0.5, strand info = 0 → gDNA and nRNA are indistinguishable")
    print("by strand evidence. The EM must rely solely on intronic vs exonic")
    print("evidence to separate nRNA from mRNA+gDNA.")
    print()
    hdr = f"{'pat':>4} {'gf':>4} | {'mrna_exp':>8} {'B_obs':>7} {'P_obs':>7} | {'nrna_exp':>8} {'B_obs':>7} {'P_obs':>7} | {'gdna_exp':>8} {'B_obs':>7} {'P_obs':>7} | {'B_tot':>7} {'P_tot':>7}"
    print(hdr)
    print("-" * len(hdr))

    for i, (b, p) in enumerate(zip(base, ph1)):
        ss = float(b["strand_specificity"])
        pat = classify_pattern(b)
        if ss != 0.5 or pat not in (4, 8, 15):
            continue
        gf = float(b["gdna_fraction"])
        if gf not in (0.0, 0.3):
            continue
        mrna_e = float(b["total_mrna_expected"])
        mrna_b = float(b["total_mrna_observed"])
        mrna_p = float(p["total_mrna_observed"])
        nrna_e = float(b["nrna_expected"])
        nrna_b = float(b["nrna_observed"])
        nrna_p = float(p["nrna_observed"])
        gdna_e = float(b["gdna_expected"])
        gdna_b = float(b["gdna_observed"])
        gdna_p = float(p["gdna_observed"])
        tot_b = mrna_b + nrna_b + gdna_b
        tot_p = mrna_p + nrna_p + gdna_p
        pl = PATTERN_LABELS.get(pat, "?")[:4]
        print(
            f"{pat:>4} {gf:>4.1f} | {mrna_e:>8.0f} {mrna_b:>7.0f} {mrna_p:>7.0f} | "
            f"{nrna_e:>8.0f} {nrna_b:>7.0f} {nrna_p:>7.0f} | "
            f"{gdna_e:>8.0f} {gdna_b:>7.0f} {gdna_p:>7.0f} | "
            f"{tot_b:>7.0f} {tot_p:>7.0f}"
        )

    # ===================================================================
    # SECTION 2: gDNA underestimation at SS≥0.75 — where does mass go?
    # ===================================================================
    print()
    print("=" * 100)
    print("  2. gDNA UNDER-ESTIMATION AT SS≥0.75 (NEW Phase 1 issue)")
    print("=" * 100)
    print()
    print("Phase 1 removed the β prior on η which constrained nRNA<->mRNA ratio.")
    print("Now nRNA freely absorbs gDNA's intronic reads at high SS, since both")
    print("nRNA and gDNA produce intronic reads but only nRNA has strand signal.")
    print()

    # For each nRNA level and SS, show gDNA error breakdown
    gdna_analysis = defaultdict(list)
    for i, (b, p) in enumerate(zip(base, ph1)):
        ss = float(b["strand_specificity"])
        gf = float(b["gdna_fraction"])
        pat = classify_pattern(b)
        if gf == 0.0:
            continue
        nrna_e = float(b["nrna_expected"])
        if nrna_e == 0:
            nl = "none"
        elif nrna_e < 3000:
            nl = "moderate"
        elif nrna_e < 8000:
            nl = "heavy"
        else:
            nl = "complex"

        gdna_e = float(b["gdna_expected"])
        gdna_b = float(b["gdna_observed"])
        gdna_p = float(p["gdna_observed"])
        nrna_b = float(b["nrna_observed"])
        nrna_p = float(p["nrna_observed"])

        gdna_analysis[(nl, ss, gf)].append({
            "gdna_exp": gdna_e,
            "gdna_base": gdna_b,
            "gdna_ph1": gdna_p,
            "nrna_exp": nrna_e,
            "nrna_base": nrna_b,
            "nrna_ph1": nrna_p,
        })

    print("  Key insight: gDNA is now systematically UNDER-estimated at SS≥0.75")
    print("  while nRNA is OVER-estimated. The freed nRNA component absorbs gDNA mass.")
    print()
    hdr = f"{'nRNA_lvl':>10} {'SS':>5} {'gf':>5} | {'gDNA base%':>11} {'gDNA ph1%':>11} {'Δ':>6} | {'nRNA base%':>11} {'nRNA ph1%':>11} {'Δ':>6}"
    print(hdr)
    print("-" * len(hdr))

    for nl in ["none", "moderate", "heavy"]:
        for ss in [0.75, 0.9, 1.0]:
            for gf in [0.1, 0.3, 0.5]:
                data = gdna_analysis.get((nl, ss, gf), [])
                if not data:
                    continue
                gdna_b_rel = np.mean([(d["gdna_base"] - d["gdna_exp"]) / d["gdna_exp"] for d in data]) * 100
                gdna_p_rel = np.mean([(d["gdna_ph1"] - d["gdna_exp"]) / d["gdna_exp"] for d in data]) * 100
                delta_g = gdna_p_rel - gdna_b_rel
                if data[0]["nrna_exp"] > 0:
                    nrna_b_rel = np.mean([(d["nrna_base"] - d["nrna_exp"]) / d["nrna_exp"] for d in data]) * 100
                    nrna_p_rel = np.mean([(d["nrna_ph1"] - d["nrna_exp"]) / d["nrna_exp"] for d in data]) * 100
                    delta_n = nrna_p_rel - nrna_b_rel
                    print(
                        f"{nl:>10} {ss:>5} {gf:>5.1f} | {gdna_b_rel:>+10.1f}% {gdna_p_rel:>+10.1f}% {delta_g:>+5.0f} | "
                        f"{nrna_b_rel:>+10.1f}% {nrna_p_rel:>+10.1f}% {delta_n:>+5.0f}"
                    )
                else:
                    print(
                        f"{nl:>10} {ss:>5} {gf:>5.1f} | {gdna_b_rel:>+10.1f}% {gdna_p_rel:>+10.1f}% {delta_g:>+5.0f} | "
                        f"{'n/a':>11} {'n/a':>11} {'':>6}"
                    )

    # ===================================================================
    # SECTION 3: Pattern 15 improvement (nRNA-only stress test)
    # ===================================================================
    print()
    print("=" * 100)
    print("  3. PATTERN 15 STRESS TEST: nRNA only (no mRNA) — mRNA false positives")
    print("=" * 100)
    print()
    print("This is the critical test from RC3: when only nRNA is expressed,")
    print("the old Beta MAP prevented η from reaching 1.0, forcing mRNA FPs.")
    print("Phase 1 removes this constraint.")
    print()
    hdr = f"{'SS':>5} {'gf':>5} {'TB':>3} | {'B mrna_FP':>10} {'P mrna_FP':>10} {'improv%':>8} | {'B nrna%':>9} {'P nrna%':>9} | {'B gdna%':>9} {'P gdna%':>9}"
    print(hdr)
    print("-" * len(hdr))

    for i, (b, p) in enumerate(zip(base, ph1)):
        pat = classify_pattern(b)
        if pat != 15:
            continue
        ss = float(b["strand_specificity"])
        gf = float(b["gdna_fraction"])
        tb1 = float(b["TB1"])
        tb2 = float(b["TB2"])
        tb = f"{int(tb1)}/{int(tb2)}"
        mrna_b = float(b["total_mrna_observed"])
        mrna_p = float(p["total_mrna_observed"])
        nrna_e = float(b["nrna_expected"])
        nrna_b = float(b["nrna_observed"])
        nrna_p = float(p["nrna_observed"])
        gdna_e = float(b["gdna_expected"]) if float(b["gdna_expected"]) > 0 else None
        gdna_b = float(b["gdna_observed"])
        gdna_p = float(p["gdna_observed"])

        improv = ((mrna_b - mrna_p) / mrna_b * 100) if mrna_b > 1 else 0
        nrna_b_r = (nrna_b - nrna_e) / nrna_e * 100
        nrna_p_r = (nrna_p - nrna_e) / nrna_e * 100
        if gdna_e:
            gdna_b_r = f"{(gdna_b - gdna_e) / gdna_e * 100:+.0f}%"
            gdna_p_r = f"{(gdna_p - gdna_e) / gdna_e * 100:+.0f}%"
        else:
            gdna_b_r = f"{gdna_b:.0f}"
            gdna_p_r = f"{gdna_p:.0f}"

        if gf in (0.0, 0.1, 0.5) and tb1 == 0 and tb2 == 500:
            print(
                f"{ss:>5} {gf:>5.1f} {tb:>3} | {mrna_b:>10.1f} {mrna_p:>10.1f} {improv:>+7.0f}% | "
                f"{nrna_b_r:>+8.0f}% {nrna_p_r:>+8.0f}% | {gdna_b_r:>9} {gdna_p_r:>9}"
            )

    # ===================================================================
    # SECTION 4: mRNA accuracy improvement summary
    # ===================================================================
    print()
    print("=" * 100)
    print("  4. mRNA ACCURACY: Baseline vs Phase 1 (delta table)")
    print("=" * 100)
    print()

    mrna_base = defaultdict(list)
    mrna_ph1 = defaultdict(list)
    for i, (b, p) in enumerate(zip(base, ph1)):
        ss = float(b["strand_specificity"])
        gf = float(b["gdna_fraction"])
        pat = classify_pattern(b)
        nrna_e = float(b["nrna_expected"])
        if nrna_e == 0:
            nl = "none"
        elif nrna_e < 3000:
            nl = "moderate"
        elif nrna_e < 8000:
            nl = "heavy"
        else:
            nl = "complex"
        mrna_e = float(b["total_mrna_expected"])
        if mrna_e > 0:
            mrna_base[(nl, ss, gf)].append(abs(float(b["total_mrna_observed"]) - mrna_e) / mrna_e)
            mrna_ph1[(nl, ss, gf)].append(abs(float(p["total_mrna_observed"]) - mrna_e) / mrna_e)

    hdr = f"{'nRNA_lvl':>10} {'gf':>5} | {'SS=0.5':>20} | {'SS=0.75':>20} | {'SS=0.9':>20} | {'SS=1.0':>20}"
    print(hdr)
    print(f"{'':>10} {'':>5} | {'base':>8} {'ph1':>8} {'Δ':>4} | {'base':>8} {'ph1':>8} {'Δ':>4} | {'base':>8} {'ph1':>8} {'Δ':>4} | {'base':>8} {'ph1':>8} {'Δ':>4}")
    print("-" * 110)
    for nl in ["none", "moderate", "heavy", "complex"]:
        for gf in [0.0, 0.1, 0.3, 0.5]:
            parts = []
            parts.append(f"{nl:>10} {gf:>5.1f} |")
            for ss in [0.5, 0.75, 0.9, 1.0]:
                bv = mrna_base.get((nl, ss, gf), [])
                pv = mrna_ph1.get((nl, ss, gf), [])
                if bv:
                    bm = np.mean(bv) * 100
                    pm = np.mean(pv) * 100
                    d = pm - bm
                    parts.append(f" {bm:>7.1f}% {pm:>7.1f}% {d:>+3.0f} |")
                else:
                    parts.append(f" {'n/a':>8} {'n/a':>8} {'':>4} |")
            print("".join(parts))

    # ===================================================================
    # SECTION 5: Root cause breakdown
    # ===================================================================
    print()
    print("=" * 100)
    print("  5. ROOT CAUSE ANALYSIS")
    print("=" * 100)
    print()

    # Count improvements vs regressions
    improved = 0
    regressed = 0
    unchanged = 0
    total = 0

    for i, (b, p) in enumerate(zip(base, ph1)):
        mrna_e = float(b["total_mrna_expected"])
        if mrna_e == 0:
            continue
        total += 1
        base_err = abs(float(b["total_mrna_observed"]) - mrna_e) / mrna_e
        ph1_err = abs(float(p["total_mrna_observed"]) - mrna_e) / mrna_e
        if ph1_err < base_err - 0.005:
            improved += 1
        elif ph1_err > base_err + 0.005:
            regressed += 1
        else:
            unchanged += 1

    print(f"  Total runs with mRNA expected > 0: {total}")
    print(f"  Improved (>0.5% lower error):     {improved} ({improved/total*100:.1f}%)")
    print(f"  Regressed (>0.5% higher error):    {regressed} ({regressed/total*100:.1f}%)")
    print(f"  Unchanged (within ±0.5%):          {unchanged} ({unchanged/total*100:.1f}%)")

    # Break regressions by SS
    print()
    print("  Regressions by SS:")
    for ss in [0.5, 0.75, 0.9, 1.0]:
        reg_count = 0
        reg_examples = []
        for i, (b, p) in enumerate(zip(base, ph1)):
            if float(b["strand_specificity"]) != ss:
                continue
            mrna_e = float(b["total_mrna_expected"])
            if mrna_e == 0:
                continue
            base_err = abs(float(b["total_mrna_observed"]) - mrna_e) / mrna_e
            ph1_err = abs(float(p["total_mrna_observed"]) - mrna_e) / mrna_e
            if ph1_err > base_err + 0.005:
                reg_count += 1
                pat = classify_pattern(b)
                gf = float(b["gdna_fraction"])
                reg_examples.append((pat, gf, base_err * 100, ph1_err * 100))
        if reg_count > 0:
            print(f"    SS={ss}: {reg_count} regressions")
            # Show worst 3
            reg_examples.sort(key=lambda x: x[3] - x[2], reverse=True)
            for pat, gf, be, pe in reg_examples[:3]:
                print(f"      pat={pat} ({PATTERN_LABELS.get(pat, '?')[:25]:25s}) gf={gf} base={be:.1f}% → ph1={pe:.1f}%")

    # Categorize the gDNA underestimation issue
    print()
    print("  gDNA absorption by nRNA (nRNA over-estimate → gDNA under-estimate):")
    print("  This is Phase 1's primary new failure mode.")
    print()
    print("  At SS≥0.75 with gDNA>0 and nRNA>0:")
    for ss in [0.75, 0.9, 1.0]:
        n_cases = 0
        gdna_base_errs = []
        gdna_ph1_errs = []
        nrna_base_errs = []
        nrna_ph1_errs = []
        for i, (b, p) in enumerate(zip(base, ph1)):
            if float(b["strand_specificity"]) != ss:
                continue
            gf = float(b["gdna_fraction"])
            nrna_e = float(b["nrna_expected"])
            gdna_e = float(b["gdna_expected"])
            if gf == 0 or nrna_e == 0:
                continue
            n_cases += 1
            gdna_base_errs.append((float(b["gdna_observed"]) - gdna_e) / gdna_e)
            gdna_ph1_errs.append((float(p["gdna_observed"]) - gdna_e) / gdna_e)
            nrna_base_errs.append((float(b["nrna_observed"]) - nrna_e) / nrna_e)
            nrna_ph1_errs.append((float(p["nrna_observed"]) - nrna_e) / nrna_e)

        if n_cases > 0:
            print(f"    SS={ss} (n={n_cases}): gDNA base={np.mean(gdna_base_errs)*100:+.1f}% → ph1={np.mean(gdna_ph1_errs)*100:+.1f}%  |  "
                  f"nRNA base={np.mean(nrna_base_errs)*100:+.1f}% → ph1={np.mean(nrna_ph1_errs)*100:+.1f}%")

    # ===================================================================
    # SECTION 6: Diagnosis — why does nRNA absorb gDNA?
    # ===================================================================
    print()
    print("=" * 100)
    print("  6. DIAGNOSIS: Why does nRNA absorb gDNA mass?")
    print("=" * 100)
    print()
    print("  Mechanism: gDNA fragments span the entire genome uniformly, including")
    print("  intronic regions. nRNA components also cover introns. Without the Beta")
    print("  prior constraint, the freed nRNA component greedily absorbs intronic")
    print("  gDNA reads because:")
    print("    (a) nRNA has higher effective length in introns than gDNA (nRNA spans")
    print("        are targeted, gDNA is genome-wide)")
    print("    (b) The strand symmetry penalty on gDNA discourages gDNA from absorbing")
    print("        reads that have strand signal, but at SS=0.5 this penalty is zero")
    print("    (c) At SS≥0.75, intronic gDNA reads DO have strand signal, making them")
    print("        look more like nRNA than gDNA in the likelihood")
    print()
    print("  This is exactly the Phase 3 issue identified in the plan:")
    print("  'The strand symmetry penalty assumes balanced strand coverage implies gDNA.'")
    print("  Phase 3 (scale strand penalty by SS) should help by making the penalty")
    print("  stronger at high SS, more correctly penalizing the asymmetric gDNA component.")
    print()
    print("  Summary of inter-component mass flow:")
    print("    - Baseline: β prior restricts nRNA → nRNA under-estimated → mass goes to gDNA (over-est)")
    print("    - Phase 1:  no β prior → nRNA over-estimated → mass comes FROM gDNA (under-est)")
    print("    - The truth lies in between: need a mechanism that correctly separates")
    print("      intronic gDNA from intronic nRNA without an artificial β constraint")


if __name__ == "__main__":
    main()
