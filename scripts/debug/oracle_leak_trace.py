"""Trace the DOWNSTREAM gDNA siphon genomically, with the calibration held at ORACLE (validated
scripts/debug/oracle.py). Per locus: leak = gdna_prior_count (TRUE gDNA under oracle calib) - assigned gDNA
(EM). Ranks loci by leak, shows the transcript structure (single-exon / multi-iso / nascent) and where the
leaked gDNA lands (nascent siphon / mature). Answers: WHERE (which regions/transcripts) and is it a
systematic anti-DNA bias?
"""
import os
import sys
os.environ.setdefault("OMP_NUM_THREADS", "1")
from dataclasses import replace as dc
import dataclasses
from pathlib import Path
import numpy as np
import pandas as pd
sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from oracle import OracleTruth
from rigel.index import TranscriptIndex
from rigel.config import PipelineConfig
from rigel.pipeline import scan_and_buffer, quant_from_buffer, _native_detect_sj_tag
from rigel.calibration import calibrate
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.splice import SpliceType

S = os.environ.get("RIGEL_SUITE", "/Users/mkiyer/Downloads/rigel_runs/quick_3to1_5mb")
COND = sys.argv[1] if len(sys.argv) > 1 else "gdna_gdna300_ss_0.99_nrna_none_capture_on"
TOPN = int(sys.argv[2]) if len(sys.argv) > 2 else 15
bam = f"{S}/{COND}/sim_oracle.bam"
index = TranscriptIndex.load(f"{S}/rigel_index")
cfg = PipelineConfig()
ra = RegionArrays.from_index(index)
wd = Path(os.environ.get("RIGEL_SCRATCH", "/tmp")) / "rigel_oracle_split"

orc = OracleTruth.from_bam(bam, index, cfg, wd, f"{os.path.basename(S)}_{COND}")
override = orc.override_masses(ra)

sc = dc(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))
stats, sm, flm, buffer, payload = scan_and_buffer(bam, index, sc)
fl = build_fl_models(global_counts=flm.global_model.counts,
                     rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                     gdna_counts=gdna_fl_mass(payload), max_size=flm.max_size)
cal = calibrate(payload=payload, region_arrays=ra, strand_model=sm,
                gdna_fl_pmf=fl.gdna_pmf, rna_fl_pmf=fl.rna_pmf, config=cfg.calibration)
cal = dataclasses.replace(cal, **override)  # ORACLE calibration
est, _ = quant_from_buffer(buffer, index, sm, fl, ra, stats, cal, payload,
                           em_config=cfg.em, scoring=cfg.scoring)

# per-transcript EM mass + structure
em = np.asarray(est.em_counts, np.float64).sum(axis=1)
tdf = index.t_df
lid = np.asarray(est.locus_id_per_transcript, np.int64)
is_syn = tdf["is_synthetic"].to_numpy(bool)
is_nr = tdf["is_nrna"].to_numpy(bool)
n_ex = tdf["n_exons"].to_numpy(int)
se = (n_ex == 1) & ~is_syn      # single-exon (mature+nascent, real)
me = (n_ex > 1) & ~is_syn & ~is_nr  # multi-exon mature
abund_arr = pd.to_numeric(tdf["abundance"], errors="coerce").fillna(0).to_numpy()
ref_arr = tdf["ref"].to_numpy()
start_arr = tdf["start"].to_numpy()
end_arr = tdf["end"].to_numpy()
strand_arr = np.array([str(s) for s in tdf["strand"].to_numpy()])
gid_arr = tdf["g_id"].to_numpy()

lr = pd.DataFrame(est.locus_results)
print(f"=== {COND} — DOWNSTREAM leak trace (ORACLE calibration) ===")
print(f"locus_results cols: {list(lr.columns)}")
tot_prior = lr["gdna_prior_count"].sum() if "gdna_prior_count" in lr else float("nan")
tot_assigned = lr["gdna"].sum() if "gdna" in lr else float("nan")
print(f"Σ gDNA prior(true)={tot_prior:,.0f}  Σ gDNA assigned={tot_assigned:,.0f}  "
      f"total leak(prior-assigned)={tot_prior-tot_assigned:+,.0f}")

# per-locus aggregation
rows = []
for _, r in lr.iterrows():
    L = int(r["locus_id"])
    m = lid == L
    tx = np.where(m)[0]
    leak = float(r.get("gdna_prior_count", 0)) - float(r.get("gdna", 0))
    mtx = tx[~is_syn[tx]]  # non-synthetic (real) transcripts in the locus
    strands = set(strand_arr[mtx].tolist())
    npos = int(np.sum(strand_arr[mtx] == "+")) + int(np.sum(strand_arr[mtx] == "1"))
    nneg = int(np.sum(strand_arr[mtx] == "-")) + int(np.sum(strand_arr[mtx] == "-1"))
    ambig = (npos > 0 and nneg > 0)
    rows.append(dict(
        locus=L, ref=str(ref_arr[tx][0]) if tx.size else "?",
        start=int(start_arr[tx].min()) if tx.size else -1,
        end=int(end_arr[tx].max()) if tx.size else -1,
        prior_g=float(r.get("gdna_prior_count", 0)), assigned_g=float(r.get("gdna", 0)), leak=leak,
        siphon=float(em[tx][is_syn[tx]].sum()), mature=float(em[tx][~is_syn[tx]].sum()),
        n_tx=int(r.get("n_transcripts", tx.size)), n_se=int(se[tx].sum()), n_me=int(me[tx].sum()),
        n_syn=int(is_syn[tx].sum()),
        eff_g=float(r.get("gdna_eff_len_em", 0)), enable_g=bool(r.get("enable_gdna", True)),
        n_genes=int(r.get("n_genes", 0)), npos=npos, nneg=nneg, ambig=ambig,
        n_gid=int(np.unique(gid_arr[mtx]).size) if mtx.size else 0,
        abund=float(abund_arr[tx][me[tx]].sum()) if tx.size else 0.0,
    ))
D = pd.DataFrame(rows).sort_values("leak", ascending=False)
# CENSUS: do AMBIG (opposite-strand overlap) loci EXIST here, and do they hold gDNA / leak?
nA = int(D.ambig.sum())
print(f"\nAMBIG census (ALL loci): {nA}/{len(D)} loci are opposite-strand-overlap (AMBIG); "
      f"they hold prior_g={D.loc[D.ambig,'prior_g'].sum():,.0f} gDNA and leak {D.loc[D.ambig,'leak'].sum():+,.0f} "
      f"(single-strand loci hold {D.loc[~D.ambig,'prior_g'].sum():,.0f}, leak {D.loc[~D.ambig,'leak'].sum():+,.0f})")
pos = D[D.leak > 0]
print(f"\ndirectional leak: Σ(+)={pos.leak.sum():,.0f} over {len(pos)} loci; "
      f"top {TOPN} = {pos.head(TOPN).leak.sum()/max(pos.leak.sum(),1)*100:.0f}% of it")
# where does the leak go? on nrna_none, all siphon is leaked; mature excess = leaked too
print(f"leak destination (all loci): Σ siphon(nascent)={D.siphon.sum():,.0f}  "
      f"(true nascent=0 ⇒ pure siphon)")
# structural pattern of the leaking loci
def _cls(row):
    if row.n_me >= 2: return "multi-iso"
    if row.n_me == 1 and row.n_se == 0: return "single-mRNA"
    if row.n_se >= 1 and row.n_me == 0: return "single-exon-only"
    return "mixed"
pos = pos.copy(); pos["cls"] = pos.apply(_cls, axis=1)
print("\nleak by locus CLASS:")
print(pos.groupby("cls").agg(n=("leak", "size"), leak=("leak", "sum"),
                             siphon=("siphon", "sum")).sort_values("leak", ascending=False).to_string())
# *** THE KEY QUESTION: is the leak in AMBIG (opposite-strand overlap) loci? ***
pos_all = D[D.leak > 0].copy()
pos_all["strand_cls"] = np.where(pos_all.ambig, "AMBIG(both-strand)",
                                 np.where(pos_all.npos > 0, "POS-only", "NEG-only"))
print("\n*** leak by STRAND composition (AMBIG = transcripts on BOTH strands overlap) ***")
print(pos_all.groupby("strand_cls").agg(
    n=("leak", "size"), leak=("leak", "sum"), siphon=("siphon", "sum"),
    prior_g=("prior_g", "sum")).sort_values("leak", ascending=False).to_string())
tot = pos_all.leak.sum()
print(f"AMBIG share of directional leak: {100*pos_all[pos_all.ambig].leak.sum()/max(tot,1):.0f}%  "
      f"of siphon: {100*pos_all[pos_all.ambig].siphon.sum()/max(pos_all.siphon.sum(),1):.0f}%")
# cross-tab AMBIG x #genes (multi-gene overlap)
print("\nleak by (ambig, n_genes>1):")
pos_all["multigene"] = pos_all.n_gid > 1
print(pos_all.groupby(["ambig", "multigene"]).agg(n=("leak", "size"), leak=("leak", "sum"),
      siphon=("siphon", "sum")).to_string())
# is gDNA ever DISABLED on a leaking locus (forced leak)?
n_dis = int((~D["enable_g"]).sum())
dis_leak = float(D.loc[~D["enable_g"], "leak"].clip(lower=0).sum())
print(f"\ngDNA-DISABLED loci: {n_dis}  (their +leak = {dis_leak:,.0f}) — a forced leak if >0")
# SCALING: does the leak FRACTION grow with isoform/nascent-shadow count?
D2 = D[D.prior_g > 5000].copy()  # loci with meaningful gDNA
D2["leak_frac"] = D2.leak / D2.prior_g
print("\nleak FRACTION (leak/prior_g) by #mature-isoforms (n_me), loci with prior_g>5k:")
print(D2.groupby("n_me").agg(n=("leak_frac", "size"), mean_leak_frac=("leak_frac", "mean"),
                             tot_leak=("leak", "sum")).to_string())
if len(D2) > 3:
    print(f"corr(leak_frac, n_me)={np.corrcoef(D2.n_me, D2.leak_frac)[0,1]:+.2f}  "
          f"corr(leak_frac, n_syn)={np.corrcoef(D2.n_syn, D2.leak_frac)[0,1]:+.2f}")
print(f"\nTOP {TOPN} leaking loci:")
cols = ["locus", "ref", "start", "end", "prior_g", "leak", "siphon", "mature",
        "n_tx", "n_me", "npos", "nneg", "ambig", "n_gid", "cls"]
with pd.option_context("display.width", 200, "display.max_columns", 20):
    print(pos.head(TOPN)[cols].to_string(index=False,
          formatters={c: "{:,.0f}".format for c in ["prior_g", "assigned_g", "leak", "siphon", "mature", "abund"]}))
