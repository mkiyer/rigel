#!/usr/bin/env python3
"""Analyze seed 101 benchmark results from summary.json."""
import json
import sys
from collections import defaultdict

SUMMARY = "/Users/mkiyer/Downloads/hulkrna_runs/bench_2026_02_19/seed_101/summary.json"

with open(SUMMARY) as f:
    data = json.load(f)

TOOLS = ["hulkrna", "hulkrna_mm", "salmon", "kallisto"]
GENE_TOOLS = ["hulkrna", "hulkrna_mm", "salmon", "kallisto", "htseq"]

# ── Transcript-level ──────────────────────────────────────────────
print("=" * 80)
print("TRANSCRIPT-LEVEL ANALYSIS")
print("=" * 80)

# Collect per-condition rows
rows = []
for r in data:
    region = r["region"]
    gdna = r["gdna_label"]
    ss = r["strand_specificity"]
    tm = r.get("transcript_metrics", {})
    for tool in TOOLS:
        if tool in tm:
            m = tm[tool]
            rows.append(dict(
                region=region, gdna=gdna, ss=ss, tool=tool,
                mae=m["mean_abs_error"], rmse=m["rmse"],
                pearson=m["pearson"], spearman=m["spearman"],
                time=m["elapsed_sec"],
            ))

# ── Overall means ──
print("\n--- Overall Means (transcript-level) ---")
agg = defaultdict(lambda: defaultdict(list))
for row in rows:
    for k in ("mae", "rmse", "pearson", "spearman", "time"):
        agg[row["tool"]][k].append(row[k])

header = f"{'Tool':<14} {'MAE':>10} {'RMSE':>10} {'Pearson':>10} {'Spearman':>10} {'Time(s)':>10}"
print(header)
print("-" * len(header))
for tool in TOOLS:
    vals = agg[tool]
    n = len(vals["mae"])
    print(f"{tool:<14} {sum(vals['mae'])/n:>10.1f} {sum(vals['rmse'])/n:>10.1f} "
          f"{sum(vals['pearson'])/n:>10.4f} {sum(vals['spearman'])/n:>10.4f} "
          f"{sum(vals['time'])/n:>10.2f}")

# ── MAE by gDNA level ──
print("\n--- MAE by gDNA Level (transcript) ---")
gdna_order = ["none", "low", "moderate", "high"]
gdna_agg = defaultdict(lambda: defaultdict(list))
for row in rows:
    gdna_agg[(row["gdna"], row["tool"])]["mae"].append(row["mae"])

header = f"{'gDNA':<12}" + "".join(f" {t:>14}" for t in TOOLS)
print(header)
print("-" * len(header))
for g in gdna_order:
    line = f"{g:<12}"
    for t in TOOLS:
        vals = gdna_agg[(g, t)]["mae"]
        if vals:
            line += f" {sum(vals)/len(vals):>14.1f}"
        else:
            line += f" {'N/A':>14}"
    print(line)

# ── MAE by Strand Specificity ──
print("\n--- MAE by Strand Specificity (transcript) ---")
ss_order = [0.95, 0.99, 1.0]
ss_agg = defaultdict(lambda: defaultdict(list))
for row in rows:
    ss_agg[(row["ss"], row["tool"])]["mae"].append(row["mae"])

header = f"{'SS':<12}" + "".join(f" {t:>14}" for t in TOOLS)
print(header)
print("-" * len(header))
for s in ss_order:
    line = f"{s:<12.2f}"
    for t in TOOLS:
        vals = ss_agg[(s, t)]["mae"]
        if vals:
            line += f" {sum(vals)/len(vals):>14.1f}"
        else:
            line += f" {'N/A':>14}"
    print(line)

# ── MAE by Region ──
print("\n--- MAE by Region (transcript) ---")
region_agg = defaultdict(lambda: defaultdict(list))
for row in rows:
    region_agg[(row["region"], row["tool"])]["mae"].append(row["mae"])

regions = sorted(set(r["region"] for r in data))
header = f"{'Region':<22}" + "".join(f" {t:>14}" for t in TOOLS)
print(header)
print("-" * len(header))
for reg in regions:
    line = f"{reg:<22}"
    for t in TOOLS:
        vals = region_agg[(reg, t)]["mae"]
        if vals:
            line += f" {sum(vals)/len(vals):>14.1f}"
        else:
            line += f" {'N/A':>14}"
    print(line)

# ── Pearson by Region ──
print("\n--- Pearson r by Region (transcript) ---")
pearson_agg = defaultdict(lambda: defaultdict(list))
for row in rows:
    pearson_agg[(row["region"], row["tool"])]["pearson"].append(row["pearson"])

header = f"{'Region':<22}" + "".join(f" {t:>14}" for t in TOOLS)
print(header)
print("-" * len(header))
for reg in regions:
    line = f"{reg:<22}"
    for t in TOOLS:
        vals = pearson_agg[(reg, t)]["pearson"]
        if vals:
            line += f" {sum(vals)/len(vals):>14.4f}"
        else:
            line += f" {'N/A':>14}"
    print(line)

# ── Win rate ──
print("\n--- Win Rate: hulkrna vs others (transcript MAE, lower=better) ---")
wins = defaultdict(int)
total = 0
for r in data:
    tm = r.get("transcript_metrics", {})
    if "hulkrna" not in tm:
        continue
    h_mae = tm["hulkrna"]["mean_abs_error"]
    total += 1
    for other in ["salmon", "kallisto"]:
        if other in tm:
            if h_mae < tm[other]["mean_abs_error"]:
                wins[f"hulkrna < {other}"] += 1
            elif h_mae == tm[other]["mean_abs_error"]:
                wins[f"hulkrna == {other}"] += 1
            else:
                wins[f"hulkrna > {other}"] += 1

print(f"Total conditions: {total}")
for label in sorted(wins):
    print(f"  {label}: {wins[label]}/{total} ({100*wins[label]/total:.1f}%)")

# ── hulkrna vs hulkrna_mm ──
print("\n--- hulkrna vs hulkrna_mm (multimap) ---")
mm_wins = defaultdict(int)
mm_total = 0
for r in data:
    tm = r.get("transcript_metrics", {})
    if "hulkrna" in tm and "hulkrna_mm" in tm:
        mm_total += 1
        h = tm["hulkrna"]["mean_abs_error"]
        hmm = tm["hulkrna_mm"]["mean_abs_error"]
        if h < hmm:
            mm_wins["hulkrna < hulkrna_mm"] += 1
        elif h > hmm:
            mm_wins["hulkrna > hulkrna_mm"] += 1
        else:
            mm_wins["hulkrna == hulkrna_mm"] += 1

print(f"Total: {mm_total}")
for label in sorted(mm_wins):
    print(f"  {label}: {mm_wins[label]}/{mm_total} ({100*mm_wins[label]/mm_total:.1f}%)")

# ── Gene-level ────────────────────────────────────────────────────
print("\n" + "=" * 80)
print("GENE-LEVEL ANALYSIS")
print("=" * 80)

gene_rows = []
for r in data:
    region = r["region"]
    gdna = r["gdna_label"]
    ss = r["strand_specificity"]
    gm = r.get("gene_metrics", {})
    for tool in GENE_TOOLS:
        if tool in gm:
            m = gm[tool]
            gene_rows.append(dict(
                region=region, gdna=gdna, ss=ss, tool=tool,
                mae=m["mean_abs_error"], rmse=m["rmse"],
                pearson=m["pearson"], spearman=m["spearman"],
            ))

# Gene overall means
print("\n--- Overall Means (gene-level) ---")
gagg = defaultdict(lambda: defaultdict(list))
for row in gene_rows:
    for k in ("mae", "rmse", "pearson", "spearman"):
        gagg[row["tool"]][k].append(row[k])

header = f"{'Tool':<14} {'MAE':>10} {'RMSE':>10} {'Pearson':>10} {'Spearman':>10}"
print(header)
print("-" * len(header))
for tool in GENE_TOOLS:
    vals = gagg[tool]
    n = len(vals["mae"])
    if n == 0:
        continue
    print(f"{tool:<14} {sum(vals['mae'])/n:>10.1f} {sum(vals['rmse'])/n:>10.1f} "
          f"{sum(vals['pearson'])/n:>10.4f} {sum(vals['spearman'])/n:>10.4f}")

# Gene MAE by region
print("\n--- Gene MAE by Region ---")
greg_agg = defaultdict(lambda: defaultdict(list))
for row in gene_rows:
    greg_agg[(row["region"], row["tool"])]["mae"].append(row["mae"])

header = f"{'Region':<22}" + "".join(f" {t:>14}" for t in GENE_TOOLS)
print(header)
print("-" * len(header))
for reg in regions:
    line = f"{reg:<22}"
    for t in GENE_TOOLS:
        vals = greg_agg[(reg, t)]["mae"]
        if vals:
            line += f" {sum(vals)/len(vals):>14.1f}"
        else:
            line += f" {'N/A':>14}"
    print(line)

# Gene win rate
print("\n--- Gene Win Rate: hulkrna vs others ---")
gwins = defaultdict(int)
gtotal = 0
for r in data:
    gm = r.get("gene_metrics", {})
    if "hulkrna" not in gm:
        continue
    h_mae = gm["hulkrna"]["mean_abs_error"]
    gtotal += 1
    for other in ["salmon", "kallisto", "htseq"]:
        if other in gm:
            if h_mae < gm[other]["mean_abs_error"]:
                gwins[f"hulkrna < {other}"] += 1
            elif h_mae == gm[other]["mean_abs_error"]:
                gwins[f"hulkrna == {other}"] += 1
            else:
                gwins[f"hulkrna > {other}"] += 1

print(f"Total conditions: {gtotal}")
for label in sorted(gwins):
    print(f"  {label}: {gwins[label]}/{gtotal} ({100*gwins[label]/gtotal:.1f}%)")

# ── Worst conditions for hulkrna ──
print("\n" + "=" * 80)
print("WORST CONDITIONS FOR HULKRNA (transcript MAE)")
print("=" * 80)

# Compute ratio hulkrna_mae / best_other_mae
worst = []
for r in data:
    tm = r.get("transcript_metrics", {})
    if "hulkrna" not in tm:
        continue
    h_mae = tm["hulkrna"]["mean_abs_error"]
    best_other = float("inf")
    best_tool = ""
    for other in ["salmon", "kallisto"]:
        if other in tm and tm[other]["mean_abs_error"] < best_other:
            best_other = tm[other]["mean_abs_error"]
            best_tool = other
    if best_other > 0:
        ratio = h_mae / best_other
    else:
        ratio = 0
    worst.append(dict(
        region=r["region"], gdna=r["gdna_label"], ss=r["strand_specificity"],
        hulkrna_mae=h_mae, best_other_mae=best_other,
        best_tool=best_tool, ratio=ratio,
    ))

worst.sort(key=lambda x: -x["ratio"])
print(f"\n{'Region':<22} {'gDNA':<10} {'SS':>6} {'hulkrna':>10} {'Best':>10} {'Who':<10} {'Ratio':>8}")
print("-" * 78)
for w in worst[:20]:
    print(f"{w['region']:<22} {w['gdna']:<10} {w['ss']:>6.2f} "
          f"{w['hulkrna_mae']:>10.1f} {w['best_other_mae']:>10.1f} "
          f"{w['best_tool']:<10} {w['ratio']:>8.2f}x")

# ── Best conditions for hulkrna ──
print("\n" + "=" * 80)
print("BEST CONDITIONS FOR HULKRNA (transcript MAE)")
print("=" * 80)
worst.sort(key=lambda x: x["ratio"])
print(f"\n{'Region':<22} {'gDNA':<10} {'SS':>6} {'hulkrna':>10} {'Best':>10} {'Who':<10} {'Ratio':>8}")
print("-" * 78)
for w in worst[:20]:
    print(f"{w['region']:<22} {w['gdna']:<10} {w['ss']:>6.2f} "
          f"{w['hulkrna_mae']:>10.1f} {w['best_other_mae']:>10.1f} "
          f"{w['best_tool']:<10} {w['ratio']:>8.2f}x")

# ── Speed comparison ──
print("\n" + "=" * 80)
print("SPEED COMPARISON (seconds)")
print("=" * 80)
time_agg = defaultdict(list)
for row in rows:
    time_agg[row["tool"]].append(row["time"])

for tool in TOOLS:
    vals = time_agg[tool]
    if vals:
        print(f"{tool:<14}  mean={sum(vals)/len(vals):.2f}  "
              f"min={min(vals):.2f}  max={max(vals):.2f}  "
              f"total={sum(vals):.1f}")

print("\nDone.")
