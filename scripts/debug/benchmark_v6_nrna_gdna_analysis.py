#!/usr/bin/env python3
"""Analyze nRNA and gDNA estimation quality from benchmark v6.

Reads summary.json (pool-level mRNA/nRNA/gDNA counts) and per-transcript CSVs
(nrna_truth column) to assess how well Rigel estimates nascent RNA and genomic DNA.
"""

import json
import pandas as pd
import numpy as np
from pathlib import Path

OUTDIR = Path("/Users/mkiyer/Downloads/rigel_runs/benchmark_output_v6_prune")


def load_summary_json():
    with open(OUTDIR / "summary.json") as f:
        return json.load(f)


def load_per_transcript(condition, aligner):
    p = OUTDIR / condition / f"per_transcript_counts_{aligner}.csv"
    if not p.exists():
        return None
    return pd.read_csv(p)


def section(title):
    bar = "=" * 90
    print(f"\n{bar}\n{title}\n{bar}\n")


def subsection(title):
    print(f"\n--- {title} ---\n")


# ═══════════════════════════════════════════════════════════
# 1. POOL-LEVEL TRIPARTITE ACCURACY (all conditions, all tools)
# ═══════════════════════════════════════════════════════════
def analyze_pool_level(data):
    section("1. POOL-LEVEL TRIPARTITE ACCURACY: mRNA / nRNA / gDNA")

    print("How well does Rigel estimate the total fragment counts for each pool?")
    print("(Salmon/kallisto do not model nRNA or gDNA, so only Rigel tools are shown.)\n")

    header = (
        f"{'Condition':<40} {'Aligner':<10} {'Tool':<25} "
        f"{'Pool':<6} {'Truth':>12} {'Predicted':>12} {'Error%':>9}"
    )
    print(header)
    print("-" * len(header))

    # Collect aggregates for summary
    rows = []

    for r in data:
        ds = r["dataset_name"]
        aligner = r["aligner"]
        n_mrna = r["n_mrna_truth"]
        n_nrna = r["n_nrna_truth"]
        n_gdna = r["n_gdna_truth"]

        pc = r.get("pool_counts", {})
        for tool, c in sorted(pc.items()):
            if "rigel" not in tool:
                continue
            mrna_pred = c["mature_rna"]
            nrna_pred = c["nascent_rna"]
            gdna_pred = c["genomic_dna"]

            mrna_err = (mrna_pred - n_mrna) / n_mrna * 100 if n_mrna > 0 else 0
            nrna_err = (nrna_pred - n_nrna) / n_nrna * 100 if n_nrna > 0 else 0
            if n_gdna > 0:
                gdna_err = (gdna_pred - n_gdna) / n_gdna * 100
            else:
                gdna_err = None  # can't compute % error when truth=0

            print(f"{ds:<40} {aligner:<10} {tool:<25} {'mRNA':<6} {n_mrna:>12} {mrna_pred:>12.0f} {mrna_err:>+8.2f}%")
            print(f"{'':>40} {'':>10} {'':>25} {'nRNA':<6} {n_nrna:>12} {nrna_pred:>12.0f} {nrna_err:>+8.2f}%")
            if gdna_err is not None:
                print(f"{'':>40} {'':>10} {'':>25} {'gDNA':<6} {n_gdna:>12} {gdna_pred:>12.0f} {gdna_err:>+8.2f}%")
            else:
                print(f"{'':>40} {'':>10} {'':>25} {'gDNA':<6} {n_gdna:>12} {gdna_pred:>12.0f}  (truth=0)")

            rows.append({
                "dataset": ds, "aligner": aligner, "tool": tool,
                "mrna_truth": n_mrna, "mrna_pred": mrna_pred, "mrna_err_pct": mrna_err,
                "nrna_truth": n_nrna, "nrna_pred": nrna_pred, "nrna_err_pct": nrna_err,
                "gdna_truth": n_gdna, "gdna_pred": gdna_pred,
                "gdna_err_pct": gdna_err if gdna_err is not None else float("nan"),
            })

    return pd.DataFrame(rows)


# ═══════════════════════════════════════════════════════════
# 2. POOL-LEVEL SUMMARY BY CONDITION TYPE
# ═══════════════════════════════════════════════════════════
def analyze_pool_summary(pool_df):
    section("2. POOL-LEVEL SUMMARY BY gDNA LEVEL AND ALIGNER")

    # Focus on rigel_default (pruning doesn't matter here)
    for aligner in ["oracle", "minimap2"]:
        for tool_pat in ["rigel_default"]:
            subset = pool_df[
                (pool_df["aligner"] == aligner) & pool_df["tool"].str.contains(tool_pat)
            ]
            if subset.empty:
                continue

            subsection(f"{tool_pat}_{aligner}")

            # By gDNA level
            for gdna_label in ["none", "low", "high"]:
                gdna_rows = subset[subset["dataset"].str.contains(f"gdna_{gdna_label}")]
                if gdna_rows.empty:
                    continue

                avg_mrna = gdna_rows["mrna_err_pct"].mean()
                avg_nrna = gdna_rows["nrna_err_pct"].mean()
                gdna_errs = gdna_rows["gdna_err_pct"].dropna()
                avg_gdna = gdna_errs.mean() if len(gdna_errs) > 0 else float("nan")

                print(f"  gDNA={gdna_label:<6}: "
                      f"mRNA err={avg_mrna:>+7.2f}%  "
                      f"nRNA err={avg_nrna:>+7.2f}%  "
                      f"gDNA err={'N/A' if np.isnan(avg_gdna) else f'{avg_gdna:>+7.2f}%'}")

            # Overall averages
            avg_mrna = subset["mrna_err_pct"].mean()
            avg_nrna = subset["nrna_err_pct"].mean()
            gdna_errs = subset["gdna_err_pct"].dropna()
            avg_gdna = gdna_errs.mean() if len(gdna_errs) > 0 else float("nan")
            print(f"  {'OVERALL':<12}: "
                  f"mRNA err={avg_mrna:>+7.2f}%  "
                  f"nRNA err={avg_nrna:>+7.2f}%  "
                  f"gDNA err={'N/A' if np.isnan(avg_gdna) else f'{avg_gdna:>+7.2f}%'}")


# ═══════════════════════════════════════════════════════════
# 3. nRNA LEAKAGE ANALYSIS: Where do nRNA fragments go?
# ═══════════════════════════════════════════════════════════
def analyze_nrna_leakage(data):
    section("3. nRNA FRAGMENT ACCOUNTING: Where do the fragments go?")

    print("Total fragments = mRNA + nRNA + gDNA. How does prediction compare to truth?")
    print("(Positive error = overestimation, negative = underestimation)\n")

    for r in data:
        ds = r["dataset_name"]
        aligner = r["aligner"]
        n_mrna = r["n_mrna_truth"]
        n_nrna = r["n_nrna_truth"]
        n_gdna = r["n_gdna_truth"]
        total_truth = n_mrna + n_nrna + n_gdna

        pc = r.get("pool_counts", {})
        tool = None
        for t in sorted(pc.keys()):
            if "rigel_default" in t:
                tool = t
                break
        if not tool:
            continue

        c = pc[tool]
        total_pred = c["mature_rna"] + c["nascent_rna"] + c["genomic_dna"]

        delta_mrna = c["mature_rna"] - n_mrna
        delta_nrna = c["nascent_rna"] - n_nrna
        delta_gdna = c["genomic_dna"] - n_gdna
        delta_total = total_pred - total_truth

        # Compute fraction of nRNA truth that "leaked" to mRNA (mRNA overestimate ÷ nRNA truth)
        nrna_to_mrna_leak = delta_mrna / n_nrna * 100 if n_nrna > 0 and delta_mrna > 0 else 0
        # Fraction of nRNA truth that "leaked" to gDNA
        nrna_to_gdna_leak = delta_gdna / n_nrna * 100 if n_nrna > 0 and delta_gdna > 0 else 0

        print(f"{ds} / {aligner} / {tool}")
        print(f"  Total frags: truth={total_truth:>12}  pred={total_pred:>12.0f}  Δ={delta_total:>+10.0f}")
        print(f"  mRNA:  Δ={delta_mrna:>+10.0f}  ({delta_mrna/n_mrna*100:>+6.2f}% of mRNA truth)")
        print(f"  nRNA:  Δ={delta_nrna:>+10.0f}  ({delta_nrna/n_nrna*100:>+6.2f}% of nRNA truth)")
        if n_gdna > 0:
            print(f"  gDNA:  Δ={delta_gdna:>+10.0f}  ({delta_gdna/n_gdna*100:>+6.2f}% of gDNA truth)")
        else:
            print(f"  gDNA:  Δ={delta_gdna:>+10.0f}  (truth=0, pred={c['genomic_dna']:.0f})")
        if nrna_to_mrna_leak > 0:
            print(f"  nRNA→mRNA leak: {nrna_to_mrna_leak:.2f}% of nRNA truth absorbed by mRNA")
        if nrna_to_gdna_leak > 0:
            print(f"  nRNA→gDNA leak: {nrna_to_gdna_leak:.2f}% of nRNA truth absorbed by gDNA")
        print()


# ═══════════════════════════════════════════════════════════
# 4. gDNA PHANTOM ESTIMATION: What happens when gDNA truth = 0?
# ═══════════════════════════════════════════════════════════
def analyze_gdna_phantom(data):
    section("4. gDNA PHANTOM PREDICTION: When truth=0, how much gDNA does Rigel predict?")

    for r in data:
        ds = r["dataset_name"]
        if "gdna_none" not in ds:
            continue
        aligner = r["aligner"]
        n_gdna = r["n_gdna_truth"]
        n_total = r["n_rna_fragments"] + n_gdna

        pc = r.get("pool_counts", {})
        for tool in sorted(pc.keys()):
            if "rigel" not in tool:
                continue
            gdna_pred = pc[tool]["genomic_dna"]
            pct_of_total = gdna_pred / n_total * 100
            print(f"  {ds:<40} {aligner:<10} {tool:<25}: "
                  f"gDNA pred={gdna_pred:>10.0f} ({pct_of_total:.2f}% of total)")


# ═══════════════════════════════════════════════════════════
# 5. gDNA ACCURACY BY CONTAMINATION LEVEL
# ═══════════════════════════════════════════════════════════
def analyze_gdna_by_level(data):
    section("5. gDNA ACCURACY BY CONTAMINATION LEVEL")

    for gdna_label in ["low", "high"]:
        subsection(f"gDNA = {gdna_label}")
        for r in data:
            ds = r["dataset_name"]
            if f"gdna_{gdna_label}" not in ds:
                continue
            aligner = r["aligner"]
            n_gdna = r["n_gdna_truth"]
            pc = r.get("pool_counts", {})

            for tool in sorted(pc.keys()):
                if "rigel" not in tool:
                    continue
                gdna_pred = pc[tool]["genomic_dna"]
                err = (gdna_pred - n_gdna) / n_gdna * 100 if n_gdna > 0 else 0
                print(f"  {ds:<40} {aligner:<10} {tool:<25}: "
                      f"truth={n_gdna:>10}  pred={gdna_pred:>10.0f}  err={err:>+7.2f}%")


# ═══════════════════════════════════════════════════════════
# 6. PER-TRANSCRIPT nRNA ACCURACY
# ═══════════════════════════════════════════════════════════
def analyze_per_transcript_nrna():
    section("6. PER-TRANSCRIPT nRNA ACCURACY")

    print("The per-transcript CSV contains 'nrna_truth' but Rigel tool columns are mRNA only.")
    print("However, we can compute per-transcript nRNA accuracy from the model's nrna_abundance")
    print("column (which is the GTF annotation, NOT the EM prediction).\n")
    print("For per-transcript nRNA prediction quality, we need to examine the nRNA pool-level")
    print("statistics, which are covered in sections above.\n")

    # What we CAN do: check if the nRNA truth distribution correlates with mRNA error
    # i.e., does high nRNA cause mRNA overestimation?
    for aligner in ["oracle", "minimap2"]:
        cond = "gdna_low_ss_0.90_nrna_default"
        ptx = load_per_transcript(cond, aligner)
        if ptx is None:
            continue

        subsection(f"nRNA Impact on mRNA Error ({cond} / {aligner})")

        truth_mrna = ptx["mrna_truth"].values.astype(float)
        truth_nrna = ptx["nrna_truth"].values.astype(float)

        meta = {
            "transcript_id", "gene_id", "gene_name",
            "mrna_abundance", "nrna_abundance", "mrna_truth", "nrna_truth",
        }
        tools = [c for c in ptx.columns if c not in meta]

        expr = truth_mrna > 0
        nrna_ratio = np.zeros_like(truth_mrna)
        nrna_ratio[expr] = truth_nrna[expr] / truth_mrna[expr]

        bins = [
            ("no nRNA (ratio=0)", nrna_ratio == 0),
            ("low nRNA (0<r<=1)", (nrna_ratio > 0) & (nrna_ratio <= 1)),
            ("mid nRNA (1<r<=5)", (nrna_ratio > 1) & (nrna_ratio <= 5)),
            ("high nRNA (r>5)", nrna_ratio > 5),
        ]

        for tool in sorted(tools):
            if "rigel" not in tool:
                continue
            obs = ptx[tool].values.astype(float)
            print(f"\n  Tool: {tool}")
            print(f"  {'nRNA Stratum':<25} {'N':>6} {'MAE':>10} {'MedAE':>10} {'Bias':>10}")

            for name, mask in bins:
                combined = expr & mask
                n = combined.sum()
                if n == 0:
                    continue
                t = truth_mrna[combined]
                o = obs[combined]
                ae = np.abs(t - o)
                bias = (o - t).mean()
                print(f"  {name:<25} {n:>6} {ae.mean():>10.2f} {np.median(ae):>10.2f} {bias:>+10.2f}")


# ═══════════════════════════════════════════════════════════
# 7. TOTAL FRAGMENT BUDGET ACCURACY
# ═══════════════════════════════════════════════════════════
def analyze_fragment_budget(data):
    section("7. TOTAL FRAGMENT BUDGET: Does Rigel account for all fragments?")

    print("If mRNA+nRNA+gDNA_pred ≈ total_fragments_truth, the model accounts for all reads.\n")

    for aligner in ["oracle", "minimap2"]:
        subsection(f"Aligner: {aligner}")
        print(f"  {'Condition':<40} {'Tool':<25} {'Total_truth':>12} {'Total_pred':>12} {'Δ':>10} {'Δ%':>8}")
        print("  " + "-" * 110)

        for r in data:
            if r["aligner"] != aligner:
                continue
            ds = r["dataset_name"]
            total_truth = r["n_mrna_truth"] + r["n_nrna_truth"] + r["n_gdna_truth"]
            pc = r.get("pool_counts", {})

            for tool in sorted(pc.keys()):
                if "rigel_default" not in tool:
                    continue
                c = pc[tool]
                total_pred = c["mature_rna"] + c["nascent_rna"] + c["genomic_dna"]
                delta = total_pred - total_truth
                pct = delta / total_truth * 100
                print(f"  {ds:<40} {tool:<25} {total_truth:>12} {total_pred:>12.0f} {delta:>+10.0f} {pct:>+7.2f}%")


# ═══════════════════════════════════════════════════════════
# 8. SUMMARY TABLE: Compact accuracy across all 3 pools
# ═══════════════════════════════════════════════════════════
def summary_table(data):
    section("8. COMPACT SUMMARY: Mean absolute pool-level error by aligner and gDNA level")

    # Group by aligner, gdna_level, compute mean abs error for each pool
    rows = []
    for r in data:
        ds = r["dataset_name"]
        aligner = r["aligner"]
        n_mrna = r["n_mrna_truth"]
        n_nrna = r["n_nrna_truth"]
        n_gdna = r["n_gdna_truth"]

        if "gdna_none" in ds:
            gdna_level = "none"
        elif "gdna_low" in ds:
            gdna_level = "low"
        elif "gdna_high" in ds:
            gdna_level = "high"
        else:
            gdna_level = "unknown"

        pc = r.get("pool_counts", {})
        for tool in sorted(pc.keys()):
            if "rigel_default" not in tool:
                continue
            c = pc[tool]
            rows.append({
                "aligner": aligner,
                "gdna_level": gdna_level,
                "tool": tool,
                "mrna_abs_err_pct": abs(c["mature_rna"] - n_mrna) / n_mrna * 100 if n_mrna > 0 else 0,
                "nrna_abs_err_pct": abs(c["nascent_rna"] - n_nrna) / n_nrna * 100 if n_nrna > 0 else 0,
                "gdna_abs_err_pct": abs(c["genomic_dna"] - n_gdna) / n_gdna * 100 if n_gdna > 0 else float("nan"),
                "mrna_err_pct": (c["mature_rna"] - n_mrna) / n_mrna * 100 if n_mrna > 0 else 0,
                "nrna_err_pct": (c["nascent_rna"] - n_nrna) / n_nrna * 100 if n_nrna > 0 else 0,
                "gdna_err_pct": (c["genomic_dna"] - n_gdna) / n_gdna * 100 if n_gdna > 0 else float("nan"),
            })

    df = pd.DataFrame(rows)

    for aligner in ["oracle", "minimap2"]:
        subsection(f"Aligner: {aligner}")
        sub = df[df["aligner"] == aligner]

        print(f"  {'gDNA Level':<10} {'mRNA Err%':>12} {'nRNA Err%':>12} {'gDNA Err%':>12}")
        print("  " + "-" * 50)

        for gdna_level in ["none", "low", "high"]:
            g = sub[sub["gdna_level"] == gdna_level]
            mrna = g["mrna_err_pct"].mean()
            nrna = g["nrna_err_pct"].mean()
            gdna_vals = g["gdna_err_pct"].dropna()
            gdna = gdna_vals.mean() if len(gdna_vals) > 0 else float("nan")
            gdna_str = f"{gdna:>+11.2f}%" if not np.isnan(gdna) else "      N/A"
            print(f"  {gdna_level:<10} {mrna:>+11.2f}% {nrna:>+11.2f}% {gdna_str}")

        overall_mrna = sub["mrna_err_pct"].mean()
        overall_nrna = sub["nrna_err_pct"].mean()
        gdna_vals = sub["gdna_err_pct"].dropna()
        overall_gdna = gdna_vals.mean() if len(gdna_vals) > 0 else float("nan")
        gdna_str = f"{overall_gdna:>+11.2f}%" if not np.isnan(overall_gdna) else "      N/A"
        print(f"  {'OVERALL':<10} {overall_mrna:>+11.2f}% {overall_nrna:>+11.2f}% {gdna_str}")


# ═══════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════
def main():
    print("=" * 90)
    print("BENCHMARK V6: nRNA AND gDNA ESTIMATION QUALITY ANALYSIS")
    print(f"Data: {OUTDIR}")
    print("=" * 90)

    data = load_summary_json()
    print(f"\nConditions: {len(data)} entries in summary.json")

    pool_df = analyze_pool_level(data)
    analyze_pool_summary(pool_df)
    analyze_nrna_leakage(data)
    analyze_gdna_phantom(data)
    analyze_gdna_by_level(data)
    analyze_per_transcript_nrna()
    analyze_fragment_budget(data)
    summary_table(data)

    print("\n" + "=" * 90)
    print("ANALYSIS COMPLETE")
    print("=" * 90)


if __name__ == "__main__":
    main()
