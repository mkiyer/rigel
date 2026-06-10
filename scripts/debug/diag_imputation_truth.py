"""AUDIT: is the POST-TRANSPORT per-region gDNA mass accurate vs oracle truth?

Compares the calibration's per-region gDNA fraction of unspliced mass — measured AFTER the
boundary-flux transport (priors._transport_boundary_flux re-attributes boundary-crossing gDNA into
regions) — to the ORACLE truth (read origin: 'gdna*'=gDNA, 'GENE*'=RNA). The earlier version of
this audit read count_gdna_frac (the CONTAINED prior mean, pre-transport) and so wrongly omitted the
crossing gDNA the transport moves into the exon. This version audits the mass the EM actually gets.
"""
import os, numpy as np, pysam
from dataclasses import replace as _dc_replace
from rigel.index import TranscriptIndex
from rigel.config import PipelineConfig
from rigel.pipeline import scan_and_buffer, _native_detect_sj_tag
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration import calibrate
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.calibration.priors import _transport_boundary_flux
from rigel.splice import SpliceType

SUITE="/Users/mkiyer/Downloads/rigel_runs/gdna_benchmark_5mb"
COND=os.environ.get("DIAG_COND","gdna_gdna1000_ss_0.99_nrna_none_capture_on")
idx=TranscriptIndex.load(f"{SUITE}/rigel_index")
df=idx.region_df
EXON_BITS=0b0011
sig=df["signature"].to_numpy(); observable=(sig & EXON_BITS)==0
starts={}; ids={}
for ref,g in df.groupby("ref_name"):
    starts[ref]=g["start"].to_numpy(); ids[ref]=g["region_id"].to_numpy()
R=len(df)
true_gdna=np.zeros(R); true_tot=np.zeros(R)
bam=f"{SUITE}/{COND}/sim_oracle.bam"
sf=pysam.AlignmentFile(bam,"rb")
for r in sf.fetch(until_eof=True):
    if r.is_secondary or r.is_supplementary or r.is_unmapped or not r.is_read1: continue
    if r.cigartuples and any(op==3 for op,_ in r.cigartuples): continue  # spliced → skip (unspliced audit)
    ref=r.reference_name
    if ref not in starts: continue
    i=np.searchsorted(starts[ref],r.reference_start,side="right")-1
    if i<0: continue
    rid=int(ids[ref][i]); true_tot[rid]+=1.0
    if r.query_name.startswith("gdna"): true_gdna[rid]+=1.0
sf.close()

cfg=PipelineConfig(); scan=_dc_replace(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))
stats,sm,flm,buffer,payload=scan_and_buffer(bam,idx,scan); sm.finalize()
ra=RegionArrays.from_region_df(idx.region_df,idx.ref_name_to_id)
fl=build_fl_models(global_counts=flm.global_model.counts,rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,gdna_counts=gdna_fl_mass(payload),max_size=flm.max_size)
sub=CalibrationSubstrate.from_payload(payload,ra)
cal=calibrate(payload=payload,region_arrays=ra,strand_model=sm,gdna_fl_pmf=fl.gdna_pmf,config=cfg.calibration)

length=np.asarray(ra.region_size_bp,dtype=np.float64)
gdna_region=_transport_boundary_flux(cal.mass_gdna_contained,cal.mass_gdna_left,cal.mass_gdna_right,length,cal.gdna_boundary_len,np.asarray(ra.ref_id))
rna_region=cal.mass_rna_contained+cal.mass_rna_left+cal.mass_rna_right
spliced_region=(np.asarray(sub.contained.mass_spliced)+np.asarray(sub.left.mass_spliced)+np.asarray(sub.right.mass_spliced))
rna_unspl=np.maximum(rna_region-spliced_region,0.0)
with np.errstate(divide="ignore",invalid="ignore"):
    est_frac=np.where(gdna_region+rna_unspl>0, gdna_region/(gdna_region+rna_unspl), np.nan)
    true_frac=np.where(true_tot>0, true_gdna/true_tot, np.nan)

# also the pre-transport contained prior mean for comparison
from rigel.calibration.density_model import node_gdna_density
from rigel.calibration.effective_length import region_eff_length, boundary_eff_length
from rigel.calibration.strand_balance import fit_strand_balance
rel=region_eff_length(ra.region_size_bp,fl.gdna_pmf); flm_=boundary_eff_length(fl.gdna_pmf)
kappa=float(fit_strand_balance(sm).rna_sense_frac)
nd=node_gdna_density(sub,ra,rel,flm_,rna_sense_frac=kappa,strand_clean_prior_weight=cfg.calibration.strand_clean_prior_weight)
pre=np.asarray(nd.count_gdna_frac)

exon=(~observable)&(true_tot>30)
w=true_tot[exon]
print(f"\n=== {COND}  kappa={kappa:.4f}  (exon regions n={exon.sum()}) ===")
print(f"  TRUE gDNA frac of unspliced (oracle):    mass-wtd={np.average(true_frac[exon],weights=w):.3f}")
print(f"  PRE-transport (count_gdna_frac):         mass-wtd={np.average(pre[exon],weights=w):.3f}")
print(f"  POST-transport (gdna_region fraction):   mass-wtd={np.average(est_frac[exon],weights=w):.3f}")
err=est_frac[exon]-true_frac[exon]
print(f"  POST-transport bias vs truth: mean={np.average(err,weights=w):+.3f}  MAE={np.average(np.abs(err),weights=w):.3f}")
