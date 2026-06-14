"""Dissect the gDNA effective-length (IPR) error: calibration per-region gDNA mass vs oracle.

gdna_eff_len is the IPR of the calibration's per-region deconvolved gDNA mass. If it is too long
(under-contracted), the calibration is over-SPREADING gDNA across the locus's regions (over-calling
depleted introns / under-concentrating on capture-enriched exons). This compares, per locus, the
eff-len the calibration produces against the eff-len it WOULD produce if its per-region gDNA equalled
the oracle truth — using the identical priors.py formula — and breaks down the worst loci by region.

Usage: python scripts/debug/dissect_efflen.py [condition]
"""
import dataclasses as dc
import sys

import numpy as np
import pandas as pd
import pysam

from rigel.calibration import calibrate
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.priors import _project_regions_to_loci, _transport_boundary_flux
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.signature import RegionType, coarse_type_array
from rigel.config import PipelineConfig
from rigel.frag_length_model import FragmentLengthModel
from rigel.index import TranscriptIndex
from rigel.locus import build_multi_loci
from rigel.pipeline import (
    _native_detect_sj_tag, _score_fragments, _setup_geometry_and_estimator, scan_and_buffer,
)
from rigel.scan import PipelineStats
from rigel.sim.read_name import parse_origin
from rigel.splice import SpliceType

S = "/Users/mkiyer/Downloads/rigel_runs/quick_3to1_5mb"
cond = sys.argv[1] if len(sys.argv) > 1 else "gdna_gdna300_ss_0.99_nrna_none_capture_on"
bam = f"{S}/{cond}/sim_oracle.bam"
idx = TranscriptIndex.load(f"{S}/rigel_index")
cfg = PipelineConfig()
scan = dc.replace(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))
stats, sm, flm, buf, pl = scan_and_buffer(bam, idx, scan)
ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
fl = build_fl_models(global_counts=flm.global_model.counts,
                     rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                     gdna_counts=gdna_fl_mass(pl), max_size=flm.max_size)
cal = calibrate(payload=pl, region_arrays=ra, strand_model=sm,
                gdna_fl_pmf=fl.gdna_pmf, rna_fl_pmf=fl.rna_pmf, config=cfg.calibration)
rna_fl = FragmentLengthModel.from_pmf(fl.rna_pmf, fl.max_size)
gdna_fl = FragmentLengthModel.from_pmf(fl.gdna_pmf, fl.max_size)
_g, est = _setup_geometry_and_estimator(idx, rna_fl, cfg.em, calibration=cal, region_arrays=ra)
em_data = _score_fragments(buf, idx, sm, rna_fl, gdna_fl, PipelineStats(), est, cfg.scoring, 1 << 30, None)
ml = build_multi_loci(em_data, idx)

# Calibration per-region transported gDNA (exactly as assemble_priors builds it) + FL-aware geom.
length = np.asarray(ra.region_size_bp, dtype=np.float64)
gdna_region = _transport_boundary_flux(cal.mass_gdna_contained, cal.mass_gdna_left,
                                       cal.mass_gdna_right, length, cal.gdna_boundary_len,
                                       np.asarray(ra.ref_id))
geom = np.asarray(cal.gdna_geom_len, dtype=np.float64)

# Oracle per-region gDNA mass: assign each gDNA read1 (unspliced; midpoint) to its region.
df = idx.region_df
starts = {r: g["start"].to_numpy() for r, g in df.groupby("ref_name")}
ids = {r: g["region_id"].to_numpy() for r, g in df.groupby("ref_name")}
R = len(df)
orac_g = np.zeros(R)
with pysam.AlignmentFile(bam, "rb") as b:
    for r in b.fetch(until_eof=True):
        if r.is_secondary or r.is_supplementary or r.is_unmapped or not r.is_read1:
            continue
        ref = r.reference_name
        if ref not in starts:
            continue
        mid = (r.reference_start + (r.reference_end or r.reference_start)) // 2
        i = int(np.searchsorted(starts[ref], mid, side="right") - 1)
        if i < 0:
            continue
        if parse_origin(r.query_name).kind == "gdna":
            orac_g[int(ids[ref][i])] += 1.0

ct = coarse_type_array(np.asarray(ra.signature))

def eff_len_from(g_region):
    """Reproduce priors.py eff_len, projected to loci, for an arbitrary per-region gDNA mass."""
    sq = np.where(geom > 0.0, g_region**2 / np.maximum(geom, 1e-9), 0.0)
    proj = _project_regions_to_loci(ra, ml, len(ml), {"g": g_region, "sq": sq, "span": geom})
    G, support, span = proj["g"], proj["sq"], proj["span"]
    safe = np.maximum(span, 1e-9)
    return np.minimum((G + 1.0) ** 2 / (support + (2.0 * G + 1.0) / safe), span), G

el_cal, G_cal = eff_len_from(gdna_region)
el_orac, G_orac = eff_len_from(orac_g)

# Join with the measured leak per locus.
plt = pd.read_csv(f"{S}/net_flow_per_locus.tsv", sep="\t")
flk = plt[(plt.gdna_label == "gdna300") & (plt.ss == 0.99) & (plt.capture == "on")] \
    .groupby("locus_id", as_index=False).net_gdna_to_rna.mean()
leak = flk.set_index("locus_id").net_gdna_to_rna

ids_loc = [m.multi_locus_id for m in ml]
d = pd.DataFrame({"locus_id": ids_loc, "eff_cal": el_cal[ids_loc], "eff_orac": el_orac[ids_loc],
                  "G_cal": G_cal[ids_loc], "G_orac": G_orac[ids_loc]})
d["eff_ratio"] = d.eff_cal / np.maximum(d.eff_orac, 1e-9)
d = d.merge(leak.rename("leak"), on="locus_id", how="left")
print(f"=== gDNA eff-len: calibration vs oracle-truth (same priors.py formula), {cond} ===")
print(f"GLOBAL mass-weighted eff_cal/eff_orac (by G_cal): "
      f"{np.sum(d.eff_cal * d.G_cal) / np.sum(d.eff_orac * d.G_cal):.2f}  "
      f"(>1 ⇒ calibration eff-len too long ⇒ gDNA over-spread)")
top = d.sort_values("leak", ascending=False).head(12)
print("\ntop-12 leak loci:  eff_cal / eff_orac / ratio | G_cal / G_orac | leak")
for _, r in top.iterrows():
    print(f"  L{int(r.locus_id):3d}: {r.eff_cal:8.0f} / {r.eff_orac:8.0f} / {r.eff_ratio:4.2f} "
          f"| {r.G_cal:8.0f} / {r.G_orac:8.0f} | {r.leak:7.0f}")

# Region breakdown for the worst-by-leak locus: where does calibration mis-place gDNA?
worst = int(top.iloc[0].locus_id)
blocks = [(blk.ref_id, blk.start, blk.end) for blk in next(m for m in ml if m.multi_locus_id == worst).loci]
sel = np.zeros(R, bool)
for rid in range(R):
    rr = int(ra.ref_id[rid]); rs = int(ra.start[rid]); re = int(ra.end[rid])
    for bref, bs, be in blocks:
        if rr == bref and min(be, re) > max(bs, rs):
            sel[rid] = True
print(f"\n=== worst locus L{worst}: per-region gDNA calib vs oracle (density = g/geom) ===")
print("region  type     geom    g_cal   g_orac  dens_cal  dens_orac")
tname = {rt.value: rt.name for rt in RegionType}
for rid in np.flatnonzero(sel):
    t = tname.get(int(ct[rid]), str(int(ct[rid])))
    dc_ = gdna_region[rid] / max(geom[rid], 1e-9); do_ = orac_g[rid] / max(geom[rid], 1e-9)
    print(f"  {rid:5d} {t:7s} {geom[rid]:7.0f} {gdna_region[rid]:8.0f} {orac_g[rid]:8.0f} "
          f"{dc_:9.3f} {do_:9.3f}")
