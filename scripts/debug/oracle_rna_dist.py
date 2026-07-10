"""Quantify the RNA over-attribution (cal - oracle) on the FULL per-node mass basis (contained + boundary
sides), matching what oracle_rna overrides. Break down by strand class + region class + mass concentration.
This is the correct basis (pass_trace's contained-region reg_true_g misses ~half the gDNA)."""
import os, sys
os.environ.setdefault("OMP_NUM_THREADS", "1")
from dataclasses import replace as dc
import numpy as np
sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from rigel.index import TranscriptIndex
from rigel.config import PipelineConfig
from rigel.pipeline import scan_and_buffer, _native_detect_sj_tag
from rigel.calibration import calibrate
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.signature import BIT_EXON_POS, BIT_EXON_NEG, BIT_INTRON_POS, BIT_INTRON_NEG
from rigel.splice import SpliceType
from _metrics import oracle_node_masses

S = "/Users/mkiyer/Downloads/rigel_runs/quick_3to1_5mb"
COND = sys.argv[1] if len(sys.argv) > 1 else "gdna_gdna300_ss_0.99_nrna_none_capture_on"
bam = f"{S}/{COND}/sim_oracle.bam"
index = TranscriptIndex.load(f"{S}/rigel_index")
cfg = PipelineConfig()
scan = dc(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))
stats, sm, flm, buffer, payload = scan_and_buffer(bam, index, scan)
sm.finalize()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
fl = build_fl_models(global_counts=flm.global_model.counts,
                     rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                     gdna_counts=gdna_fl_mass(payload), max_size=flm.max_size)
cal = calibrate(payload=payload, region_arrays=ra, strand_model=sm,
                gdna_fl_pmf=fl.gdna_pmf, rna_fl_pmf=fl.rna_pmf, config=cfg.calibration)
om = oracle_node_masses(bam, ra, index)

# per-region RNA mass on the full basis (contained + both boundary sides) — cal vs oracle
def rtot(m):
    return np.asarray(m["mass_rna_contained"]) + np.asarray(m["mass_rna_left"]) + np.asarray(m["mass_rna_right"])
def gtot(m):
    return np.asarray(m["mass_gdna_contained"]) + np.asarray(m["mass_gdna_left"]) + np.asarray(m["mass_gdna_right"])
cal_r = rtot({k: getattr(cal, k) for k in om}); orc_r = rtot(om)
cal_g = gtot({k: getattr(cal, k) for k in om}); orc_g = gtot(om)
over_rna = cal_r - orc_r                       # >0 = RNA OVER-attributed (= the leak driver)

sig = np.asarray(ra.signature).astype(np.int64)
is_exon = (sig & (BIT_EXON_POS | BIT_EXON_NEG)) != 0
has_p = (sig & (BIT_EXON_POS | BIT_INTRON_POS)) != 0
has_n = (sig & (BIT_EXON_NEG | BIT_INTRON_NEG)) != 0
scl = np.where(has_p & has_n, "AMBIG", np.where(has_p, "POS", np.where(has_n, "NEG", "NONE")))

print(f"=== {COND}: RNA over-attribution (cal-oracle, full basis) ===")
print(f"total RNA: cal={cal_r.sum():,.0f} oracle={orc_r.sum():,.0f}  OVER(+)={over_rna.sum():+,.0f}")
print(f"total gDNA: cal={cal_g.sum():,.0f} oracle={orc_g.sum():,.0f}")
pos = np.maximum(over_rna, 0); neg = np.maximum(-over_rna, 0)
print(f"directional RNA over={pos.sum():,.0f}  under={neg.sum():,.0f}  net={over_rna.sum():+,.0f}")
print("\nRNA over-attribution by strand class (regions, directional +):")
for k in ["AMBIG", "POS", "NEG", "NONE"]:
    m = scl == k
    print(f"  {k:6} n={m.sum():>4}  over={pos[m].sum():>10,.0f} ({100*pos[m].sum()/max(pos.sum(),1):4.0f}%)  "
          f"exon_over={pos[m & is_exon].sum():>10,.0f}")
print(f"\nexon vs non-exon over: exon={pos[is_exon].sum():,.0f}  non-exon={pos[~is_exon].sum():,.0f}")
order = np.argsort(pos)[::-1]
cum = np.cumsum(pos[order]) / max(pos.sum(), 1)
print("mass concentration of RNA over-attribution:")
for tn in [5, 10, 20, 50, 100]:
    print(f"  top {tn:>3}: {100*cum[tn-1]:4.0f}%")
print(f"  ({int((pos>0).sum())} regions over-attribute RNA)")
print("\nTOP 12 RNA-over-attributing regions:")
print(f"  {'reg':>5} {'scl':6} {'exon':4} {'cal_r':>10} {'orc_r':>10} {'over':>10} {'cal_g':>10} {'orc_g':>10}")
for i in order[:12]:
    print(f"  {i:>5} {scl[i]:6} {str(bool(is_exon[i])):5} {cal_r[i]:>10,.0f} {orc_r[i]:>10,.0f} "
          f"{over_rna[i]:>+10,.0f} {cal_g[i]:>10,.0f} {orc_g[i]:>10,.0f}")
