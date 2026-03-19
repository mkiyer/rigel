#!/usr/bin/env python3
"""Root cause analysis of rigel false negatives vs salmon/kallisto."""
import pandas as pd
import numpy as np

outdir = "/Users/mkiyer/Downloads/rigel_runs/benchmark_output_v5_all"
cond = "gdna_low_ss_0.90_nrna_default"
ptx = pd.read_csv(f"{outdir}/{cond}/per_transcript_counts_oracle.csv")

truth = ptx["mrna_truth"].values.astype(float)
rigel = ptx["rigel_oracle"].values.astype(float)
salmon = ptx["salmon"].values.astype(float)
kallisto = ptx["kallisto"].values.astype(float)
nrna_truth = ptx["nrna_truth"].values.astype(float)

# Count isoforms per gene
iso_counts = ptx.groupby("gene_id").size()
ptx["n_iso"] = ptx["gene_id"].map(iso_counts)

# ── 1. Rigel false negatives ──
print("=" * 100)
print("1. RIGEL FALSE NEGATIVES (truth > 10, rigel < 0.5)")
print("=" * 100)
mask = (truth > 10) & (rigel < 0.5)
print(f"Count: {mask.sum()} / {(truth > 10).sum()} expressed transcripts with truth > 10")
subset = ptx[mask].copy()
subset["truth"] = truth[mask]
subset["rigel"] = rigel[mask]
subset["salmon_pred"] = salmon[mask]
subset["kallisto_pred"] = kallisto[mask]
subset["nrna"] = nrna_truth[mask]
subset = subset.sort_values("truth", ascending=False)

print(f"\n{'Transcript':<30} {'Gene':<15} {'Truth':>8} {'nRNA':>8} {'Rigel':>8} "
      f"{'Salmon':>8} {'Kallisto':>8} {'Isoforms':>8}")
print("-" * 110)
for _, r in subset.head(30).iterrows():
    print(f"{r['transcript_id']:<30} {r['gene_name']:<15} {r['truth']:>8.0f} "
          f"{r['nrna']:>8.0f} {r['rigel']:>8.1f} {r['salmon_pred']:>8.1f} "
          f"{r['kallisto_pred']:>8.1f} {r['n_iso']:>8d}")

# ── 2. Characterize the false negatives ──
print("\n" + "=" * 100)
print("2. CHARACTERIZATION OF FALSE NEGATIVES")
print("=" * 100)
fn = ptx[mask].copy()
fn["truth"] = truth[mask]
fn["nrna"] = nrna_truth[mask]
fn["nrna_ratio"] = fn["nrna"] / fn["truth"]
print(f"\nMedian isoforms per gene: {fn['n_iso'].median():.0f}")
print(f"Mean isoforms per gene: {fn['n_iso'].mean():.1f}")
print(f"Median nRNA/mRNA ratio: {fn['nrna_ratio'].median():.1f}")
print(f"Mean nRNA/mRNA ratio: {fn['nrna_ratio'].mean():.1f}")
print(f"Median truth count: {fn['truth'].median():.0f}")
print(f"Pct in multi-isoform genes (n_iso >= 2): {100 * (fn['n_iso'] >= 2).mean():.1f}%")

# By isoform count
for iso_range, lo, hi in [("1", 1, 1), ("2", 2, 2), ("3-5", 3, 5), ("6+", 6, 999)]:
    sub = fn[(fn["n_iso"] >= lo) & (fn["n_iso"] <= hi)]
    print(f"  n_iso={iso_range}: {len(sub)} false negatives")

# ── 3. Isoform competition: for multi-isoform genes with a FN, where did counts go? ──
print("\n" + "=" * 100)
print("3. ISOFORM COMPETITION: Where did rigel allocate the missing counts?")
print("=" * 100)
fn_genes = fn[fn["n_iso"] >= 2]["gene_id"].unique()
print(f"Multi-isoform genes with at least one FN: {len(fn_genes)}")

examples = []
for gene_id in fn_genes[:30]:
    gene_txs = ptx[ptx["gene_id"] == gene_id].copy()
    gene_txs["truth"] = truth[ptx["gene_id"] == gene_id]
    gene_txs["rigel"] = rigel[ptx["gene_id"] == gene_id]
    gene_txs["salmon_pred"] = salmon[ptx["gene_id"] == gene_id]
    gene_txs["nrna"] = nrna_truth[ptx["gene_id"] == gene_id]

    gene_truth_total = gene_txs["truth"].sum()
    gene_rigel_total = gene_txs["rigel"].sum()
    gene_salmon_total = gene_txs["salmon_pred"].sum()

    if gene_truth_total < 50:
        continue

    examples.append({
        "gene_id": gene_id,
        "gene_name": gene_txs["gene_name"].iloc[0],
        "n_iso": len(gene_txs),
        "truth_total": gene_truth_total,
        "rigel_total": gene_rigel_total,
        "salmon_total": gene_salmon_total,
        "rigel_err": gene_rigel_total - gene_truth_total,
        "salmon_err": gene_salmon_total - gene_truth_total,
    })

examples_df = pd.DataFrame(examples).sort_values("truth_total", ascending=False)
print(f"\n{'Gene':<15} {'Isos':>5} {'Truth':>10} {'Rigel':>10} {'Salmon':>10} "
      f"{'Rigel_err':>10} {'Salmon_err':>10}")
print("-" * 80)
for _, r in examples_df.head(20).iterrows():
    print(f"{r['gene_name']:<15} {r['n_iso']:>5} {r['truth_total']:>10.0f} "
          f"{r['rigel_total']:>10.0f} {r['salmon_total']:>10.0f} "
          f"{r['rigel_err']:>+10.0f} {r['salmon_err']:>+10.0f}")

# ── 4. Gene-level accuracy comparison ──
print("\n" + "=" * 100)
print("4. GENE-LEVEL VS TRANSCRIPT-LEVEL: Does rigel recover at gene level?")
print("=" * 100)

ptx_agg = ptx.copy()
ptx_agg["_truth"] = truth
ptx_agg["_rigel"] = rigel
ptx_agg["_salmon"] = salmon
ptx_agg["_kallisto"] = kallisto
gene_agg = ptx_agg.groupby("gene_id")[["_truth", "_rigel", "_salmon", "_kallisto"]].sum()
gene_truth = gene_agg["_truth"]
gene_rigel = gene_agg["_rigel"]
gene_salmon = gene_agg["_salmon"]
gene_kallisto = gene_agg["_kallisto"]

expr_genes = gene_truth > 0
print(f"Expressed genes: {expr_genes.sum()}")
for tool_name, gene_pred in [("rigel", gene_rigel), ("salmon", gene_salmon), ("kallisto", gene_kallisto)]:
    ae = np.abs(gene_pred[expr_genes] - gene_truth[expr_genes])
    re = ae / gene_truth[expr_genes]
    print(f"  {tool_name:>10}: MAE={ae.mean():.2f}, MedAE={ae.median():.2f}, "
          f"MedRE={re.median()*100:.1f}%, Pearson={np.corrcoef(gene_truth[expr_genes], gene_pred[expr_genes])[0,1]:.6f}")

# ── 5. Spearman paradox investigation ──
print("\n" + "=" * 100)
print("5. SPEARMAN PARADOX: Why does salmon have higher Spearman than rigel?")
print("=" * 100)

def _spearman(x, y):
    """Simple Spearman correlation via ranks."""
    from pandas import Series
    rx = Series(x).rank()
    ry = Series(y).rank()
    return np.corrcoef(rx, ry)[0, 1]

expr_mask = truth > 0
r_sp = _spearman(truth[expr_mask], rigel[expr_mask])
s_sp = _spearman(truth[expr_mask], salmon[expr_mask])
k_sp = _spearman(truth[expr_mask], kallisto[expr_mask])
print(f"Spearman (expressed only):")
print(f"  rigel:    {r_sp:.4f}")
print(f"  salmon:   {s_sp:.4f}")
print(f"  kallisto: {k_sp:.4f}")

# Check if it's driven by zero-truth transcripts
r_sp_all = _spearman(truth, rigel)
s_sp_all = _spearman(truth, salmon)
k_sp_all = _spearman(truth, kallisto)
print(f"\nSpearman (all transcripts including zeros):")
print(f"  rigel:    {r_sp_all:.4f}")
print(f"  salmon:   {s_sp_all:.4f}")
print(f"  kallisto: {k_sp_all:.4f}")

# Check false positive rates
print(f"\nFalse positive counts (truth=0, pred>0.5):")
zero_mask = truth == 0
for name, pred in [("rigel", rigel), ("salmon", salmon), ("kallisto", kallisto)]:
    fp = (pred[zero_mask] > 0.5).sum()
    total_zero = zero_mask.sum()
    print(f"  {name:>10}: {fp:>6} / {total_zero} ({100*fp/total_zero:.2f}%)")

# Spearman among highly expressed only
high_mask = truth > 50
r_sp_high = _spearman(truth[high_mask], rigel[high_mask])
s_sp_high = _spearman(truth[high_mask], salmon[high_mask])
print(f"\nSpearman (truth > 50 only, n={high_mask.sum()}):")
print(f"  rigel:    {r_sp_high:.4f}")
print(f"  salmon:   {s_sp_high:.4f}")

# ── 6. Quantify the isoform redistribution problem ──
print("\n" + "=" * 100)
print("6. ISOFORM REDISTRIBUTION QUANTIFICATION")
print("=" * 100)

# For each expressed gene with 2+ isoforms, compare isoform proportion accuracy
multi_genes = ptx[ptx["n_iso"] >= 2]["gene_id"].unique()
rigel_prop_errors = []
salmon_prop_errors = []

for gene_id in multi_genes:
    gene_mask = ptx["gene_id"] == gene_id
    t = truth[gene_mask]
    if t.sum() < 10:
        continue
    t_prop = t / t.sum()
    r_total = rigel[gene_mask].sum()
    s_total = salmon[gene_mask].sum()
    if r_total > 0:
        r_prop = rigel[gene_mask] / r_total
        rigel_prop_errors.append(np.abs(r_prop - t_prop).mean())
    if s_total > 0:
        s_prop = salmon[gene_mask] / s_total
        salmon_prop_errors.append(np.abs(s_prop - t_prop).mean())

print(f"Multi-isoform genes with truth >= 10 frags: {len(rigel_prop_errors)}")
print(f"Mean isoform proportion error:")
print(f"  rigel:  {np.mean(rigel_prop_errors):.4f}")
print(f"  salmon: {np.mean(salmon_prop_errors):.4f}")
print(f"Median isoform proportion error:")
print(f"  rigel:  {np.median(rigel_prop_errors):.4f}")
print(f"  salmon: {np.median(salmon_prop_errors):.4f}")
