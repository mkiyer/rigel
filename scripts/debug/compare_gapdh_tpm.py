"""Compare oracle vs minimap2 GAPDH quant using TPM."""
import pandas as pd

base = "/Users/mkiyer/Downloads/rigel_runs/gapdh_debug/gdna_r0_ss_0.95_nrna_file"
truth = pd.read_csv(
    "/Users/mkiyer/Downloads/rigel_runs/gapdh_debug/truth_abundances_nrna_file.tsv",
    sep="\t",
)
oracle_q = pd.read_feather(f"{base}/rigel_oracle/quant.feather")
mm2_q = pd.read_feather(f"{base}/rigel_minimap2/quant.feather")

gapdh_tids = truth[truth["gene_name"] == "GAPDH"]["transcript_id"].tolist()

# TPM comparison
print(f"{'Transcript':30s} {'Truth':>8s} {'Oracle':>8s} {'MM2':>8s} {'Ora.Err':>8s} {'MM2.Err':>8s}")
print("-" * 80)
for tid in gapdh_tids:
    t_row = truth[truth["transcript_id"] == tid].iloc[0]
    t_mrna = t_row["mrna_abundance"]
    o = oracle_q[oracle_q["transcript_id"] == tid]
    m = mm2_q[mm2_q["transcript_id"] == tid]
    o_tpm = o["tpm"].iloc[0] if len(o) > 0 else 0.0
    m_tpm = m["tpm"].iloc[0] if len(m) > 0 else 0.0
    o_err = o_tpm - t_mrna
    m_err = m_tpm - t_mrna
    print(f"  {tid:30s} {t_mrna:8.1f} {o_tpm:8.1f} {m_tpm:8.1f} {o_err:+8.1f} {m_err:+8.1f}")

# Fragment count comparison
print()
print(f"{'Transcript':30s} {'TruthTPM':>8s} {'Ora.Cnt':>8s} {'MM2.Cnt':>8s}")
print("-" * 60)
for tid in gapdh_tids:
    t_row = truth[truth["transcript_id"] == tid].iloc[0]
    t_mrna = t_row["mrna_abundance"]
    o = oracle_q[oracle_q["transcript_id"] == tid]
    m = mm2_q[mm2_q["transcript_id"] == tid]
    o_cnt = o["mrna"].iloc[0] if len(o) > 0 else 0.0
    m_cnt = m["mrna"].iloc[0] if len(m) > 0 else 0.0
    print(f"  {tid:30s} {t_mrna:8.1f} {o_cnt:8.1f} {m_cnt:8.1f}")

# Non-GAPDH leakage (minimap2)
print()
print("Non-GAPDH leakage (mm2, tpm > 1):")
mm2_leak = mm2_q[(mm2_q["tpm"] > 1) & (~mm2_q["transcript_id"].isin(gapdh_tids))]
for _, r in mm2_leak.nlargest(10, "tpm").iterrows():
    print(f"  {r['transcript_id']:30s}  gene={r['gene_name']:15s}  tpm={r['tpm']:8.1f}")

# Totals
oracle_gapdh_tpm = oracle_q[oracle_q["transcript_id"].isin(gapdh_tids)]["tpm"].sum()
mm2_gapdh_tpm = mm2_q[mm2_q["transcript_id"].isin(gapdh_tids)]["tpm"].sum()
truth_gapdh_tpm = truth[truth["gene_name"] == "GAPDH"]["mrna_abundance"].sum()

print()
print(f"Total GAPDH TPM:   truth={truth_gapdh_tpm:.0f}  oracle={oracle_gapdh_tpm:.0f}  mm2={mm2_gapdh_tpm:.0f}")
print(f"Total all TPM:     oracle={oracle_q['tpm'].sum():.0f}  mm2={mm2_q['tpm'].sum():.0f}")

# Chimeric fragments
print()
print(f"Oracle fragments: {oracle_q['mrna'].sum():.0f} total mRNA")
print(f"MM2 fragments: {mm2_q['mrna'].sum():.0f} total mRNA")

# Non-GAPDH fragment leakage (count space)
oracle_gapdh_cnt = oracle_q[oracle_q["transcript_id"].isin(gapdh_tids)]["mrna"].sum()
mm2_gapdh_cnt = mm2_q[mm2_q["transcript_id"].isin(gapdh_tids)]["mrna"].sum()
mm2_nongapdh_cnt = mm2_q[~mm2_q["transcript_id"].isin(gapdh_tids)]["mrna"].sum()
print(f"Oracle GAPDH counts: {oracle_gapdh_cnt:.0f}")
print(f"MM2 GAPDH counts: {mm2_gapdh_cnt:.0f}")
print(f"MM2 non-GAPDH counts (leakage): {mm2_nongapdh_cnt:.0f}")
