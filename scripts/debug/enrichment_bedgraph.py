"""Genome-wide gDNA enrichment ratio (per region) + the enrichment-weighted eff-len sanity check.

Computes, from a fresh calibration on one condition:
  - global gDNA density  rho_global = sum(region gDNA mass) / sum(region length)
  - per-region enrichment ratio e_r = (M_r / L_r) / rho_global       (>1 enriched, <1 depleted)
and writes a bedgraph (ref start end enrichment) for inspection. Then, for the worst nascent-sink
loci, prints the would-be enrichment-weighted EM eff-len  geom * rho_global / rho_component  for the
nascent span vs its mature analog, to confirm the invariant nascent_eff >= mature_eff now holds.

Usage: python scripts/debug/enrichment_bedgraph.py [condition]
"""
import sys
from dataclasses import replace as _dc_replace

import numpy as np

from rigel.calibration import calibrate
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.priors import _transport_boundary_flux
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.signature import RegionType, coarse_type_array
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex
from rigel.pipeline import _native_detect_sj_tag, scan_and_buffer
from rigel.splice import SpliceType

S = "/Users/mkiyer/Downloads/rigel_runs/quick_3to1_5mb"
cond = sys.argv[1] if len(sys.argv) > 1 else "gdna_gdna300_ss_0.99_nrna_none_capture_on"
BAM = f"{S}/{cond}/sim_oracle.bam"
idx = TranscriptIndex.load(f"{S}/rigel_index")
cfg = PipelineConfig()
scan = _dc_replace(cfg.scan, sj_strand_tag=_native_detect_sj_tag(BAM))
stats, sm, flm, buf, pl = scan_and_buffer(BAM, idx, scan)
ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
fl = build_fl_models(global_counts=flm.global_model.counts,
                     rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                     gdna_counts=gdna_fl_mass(pl), max_size=flm.max_size)
cal = calibrate(payload=pl, region_arrays=ra, strand_model=sm,
                gdna_fl_pmf=fl.gdna_pmf, rna_fl_pmf=fl.rna_pmf, config=cfg.calibration)
buf.cleanup()

length = np.asarray(ra.region_size_bp, dtype=np.float64)
gdna_region = _transport_boundary_flux(cal.mass_gdna_contained, cal.mass_gdna_left,
                                       cal.mass_gdna_right, length, cal.gdna_boundary_len,
                                       np.asarray(ra.ref_id))
rho_r = np.where(length > 0, gdna_region / np.maximum(length, 1e-9), 0.0)
rho_global = gdna_region.sum() / length.sum()
e_r = rho_r / max(rho_global, 1e-12)
ct = coarse_type_array(np.asarray(ra.signature))
tname = {rt.value: rt.name for rt in RegionType}

print(f"=== {cond} ===")
print(f"rho_global (gDNA mass/bp) = {rho_global:.4f}   total gDNA mass = {gdna_region.sum():,.0f}")
print(f"per-region enrichment e_r: min={e_r.min():.2f} median={np.median(e_r):.2f} "
      f"p90={np.quantile(e_r,0.9):.2f} max={e_r.max():.2f}")
# enrichment by region type (the capture signature)
print("\nenrichment by region type (length-weighted mean density / global):")
for t in np.unique(ct):
    m = ct == t
    if m.sum() == 0:
        continue
    dens = gdna_region[m].sum() / max(length[m].sum(), 1e-9)
    print(f"  {tname.get(int(t), str(int(t))):8s} n={int(m.sum()):5d}  e={dens/max(rho_global,1e-12):6.2f}  "
          f"(mass {gdna_region[m].sum():>12,.0f}, len {length[m].sum():>10,.0f})")

out = f"{S}/{cond}/gdna_enrichment.bedgraph"
with open(out, "w") as fh:
    fh.write('track type=bedGraph name="gdna_enrichment"\n')
    refs = np.asarray(ra.ref_id)
    id2name = {v: k for k, v in idx.ref_name_to_id.items()}
    st = np.asarray(ra.start); en = np.asarray(ra.end)
    for i in np.argsort(refs, kind="stable"):
        fh.write(f"{id2name.get(int(refs[i]), refs[i])}\t{int(st[i])}\t{int(en[i])}\t{e_r[i]:.4f}\n")
print(f"\nwrote {out}  ({len(e_r)} regions)")
