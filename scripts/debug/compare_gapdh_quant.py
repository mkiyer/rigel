"""Compare oracle vs minimap2 rigel quant for GAPDH debug scenario."""
import pandas as pd
import json

base = "/Users/mkiyer/Downloads/rigel_runs/gapdh_debug/gdna_r0_ss_0.95_nrna_file"

# Load truth abundances
truth = pd.read_csv(
    "/Users/mkiyer/Downloads/rigel_runs/gapdh_debug/truth_abundances_nrna_file.tsv",
    sep="\t",
)
truth_gapdh = truth[truth["gene_name"] == "GAPDH"]

# Load quant results
oracle_q = pd.read_feather(f"{base}/rigel_oracle/quant.feather")
mm2_q = pd.read_feather(f"{base}/rigel_minimap2/quant.feather")

# Merge on transcript_id
gapdh_tx = truth_gapdh[["transcript_id", "mrna_abundance", "nrna_abundance",
                         "spliced_length"]].copy()
gapdh_tx = gapdh_tx.rename(columns={"mrna_abundance": "truth_mrna",
                                     "nrna_abundance": "truth_nrna"})

# Get oracle estimates for GAPDH transcripts
oracle_gapdh = oracle_q[oracle_q["transcript_id"].isin(gapdh_tx["transcript_id"])]
mm2_gapdh = mm2_q[mm2_q["transcript_id"].isin(gapdh_tx["transcript_id"])]

# Also check for non-GAPDH transcripts that got abundance
oracle_nongapdh = oracle_q[
    (oracle_q["mrna"] > 0.1) &
    (~oracle_q["transcript_id"].isin(gapdh_tx["transcript_id"]))
]
mm2_nongapdh = mm2_q[
    (mm2_q["mrna"] > 0.1) &
    (~mm2_q["transcript_id"].isin(gapdh_tx["transcript_id"]))
]

# Compute truth fragment counts per transcript
# The sim distributes n_rna=50000 fragments by abundance * eff_length
# We need to compute expected fragment counts
total_abund = gapdh_tx["truth_mrna"].sum()
gapdh_tx["truth_frac"] = gapdh_tx["truth_mrna"] / total_abund

print("=" * 80)
print("GAPDH Debug: Oracle vs Minimap2 Quantification")
print("=" * 80)

print("\n--- Truth Abundances ---")
for _, r in gapdh_tx.iterrows():
    print(f"  {r['transcript_id']:30s}  mrna={r['truth_mrna']:8.1f}  "
          f"length={r['spliced_length']:5d}  frac={r['truth_frac']:.4f}")

print("\n--- Oracle Quant (GAPDH transcripts) ---")
oracle_merge = gapdh_tx.merge(
    oracle_gapdh[["transcript_id", "mrna", "nrna", "mrna_unambig", "mrna_em"]],
    on="transcript_id", how="left"
).fillna(0)
for _, r in oracle_merge.iterrows():
    print(f"  {r['transcript_id']:30s}  truth={r['truth_mrna']:8.1f}  "
          f"est={r['mrna']:8.1f}  "
          f"unambig={r['mrna_unambig']:8.1f}  "
          f"em={r['mrna_em']:8.1f}  "
          f"nrna={r['nrna']:8.1f}")

print("\n--- Minimap2 Quant (GAPDH transcripts) ---")
mm2_merge = gapdh_tx.merge(
    mm2_gapdh[["transcript_id", "mrna", "nrna", "mrna_unambig", "mrna_em"]],
    on="transcript_id", how="left"
).fillna(0)
for _, r in mm2_merge.iterrows():
    print(f"  {r['transcript_id']:30s}  truth={r['truth_mrna']:8.1f}  "
          f"est={r['mrna']:8.1f}  "
          f"unambig={r['mrna_unambig']:8.1f}  "
          f"em={r['mrna_em']:8.1f}  "
          f"nrna={r['nrna']:8.1f}")

# Compute total GAPDH
print("\n--- GAPDH Gene-Level Summary ---")
truth_total = gapdh_tx["truth_mrna"].sum()
oracle_total = oracle_merge["mrna"].sum()
mm2_total = mm2_merge["mrna"].sum()
print(f"  Truth total mRNA:    {truth_total:10.1f}")
print(f"  Oracle total mRNA:   {oracle_total:10.1f}  (error: {oracle_total - truth_total:+.1f})")
print(f"  Minimap2 total mRNA: {mm2_total:10.1f}  (error: {mm2_total - truth_total:+.1f})")

# Check leakage to non-GAPDH
print(f"\n--- Non-GAPDH Leakage ---")
print(f"  Oracle:   {len(oracle_nongapdh)} transcripts with mrna > 0.1")
if len(oracle_nongapdh) > 0:
    oracle_leak = oracle_nongapdh.nlargest(10, "mrna")
    for _, r in oracle_leak.iterrows():
        print(f"    {r['transcript_id']:30s}  g={r['gene_name']:15s}  mrna={r['mrna']:8.1f}")

print(f"  Minimap2: {len(mm2_nongapdh)} transcripts with mrna > 0.1")
if len(mm2_nongapdh) > 0:
    mm2_leak = mm2_nongapdh.nlargest(20, "mrna")
    for _, r in mm2_leak.iterrows():
        print(f"    {r['transcript_id']:30s}  g={r['gene_name']:15s}  mrna={r['mrna']:8.1f}")

# Total fragment accounting
print(f"\n--- Fragment Accounting ---")
oracle_total_mrna = oracle_q["mrna"].sum()
mm2_total_mrna = mm2_q["mrna"].sum()
oracle_total_nrna = oracle_q["nrna"].sum()
mm2_total_nrna = mm2_q["nrna"].sum()
print(f"  Oracle:   mRNA={oracle_total_mrna:.0f}  nRNA={oracle_total_nrna:.0f}  total={oracle_total_mrna + oracle_total_nrna:.0f}")
print(f"  Minimap2: mRNA={mm2_total_mrna:.0f}  nRNA={mm2_total_nrna:.0f}  total={mm2_total_mrna + mm2_total_nrna:.0f}")

# Summary JSON stats
for label, path in [("Oracle", f"{base}/rigel_oracle/summary.json"),
                     ("Minimap2", f"{base}/rigel_minimap2/summary.json")]:
    with open(path) as f:
        s = json.load(f)
    stats = s.get("stats", s.get("scan_stats", {}))
    chim = stats.get("n_chimeric", 0)
    multi = stats.get("multimapping", 0)
    intergenic = stats.get("n_intergenic", 0)
    print(f"\n  {label} scan stats: chimeric={chim}, multimapping={multi}, intergenic={intergenic}")
