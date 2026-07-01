"""Does the ê(z) enrichment transfer ENGAGE on the benchmark capture conditions? (enrich_w>0?)

Loads a pre-built benchmark index, scans a condition's sim_oracle.bam, runs calibrate() with the
node_sweep capture hook, and reports enrich_w + ρ_global/ρ_floor + the per-rtype solved-vs-? gДНА.
If enrich_w≈0 under capture, the enriched exons fall back to the floor/global (→ the gdna300 exon
gDNA→RNA leak); if >0, ê(z) is present but compressed (the predictor-fix problem).
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np

from rigel.index import TranscriptIndex
from rigel.config import BamScanConfig, CalibrationConfig
from rigel.pipeline import scan_and_buffer
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.calibrate import calibrate
from rigel.calibration.node_chain import build_node_chain, REGION
from rigel.splice import SpliceType

SUITE = Path.home() / "Downloads/rigel_runs/quick_3to1_5mb"
IDX = SUITE / "rigel_index"


def probe(cond):
    bam = SUITE / cond / "sim_oracle.bam"
    idx = TranscriptIndex.load(IDX)
    _st, sm, flm, _buf, pl = scan_and_buffer(str(bam), idx, BamScanConfig())
    fl = build_fl_models(global_counts=flm.global_model.counts,
                         rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                         gdna_counts=gdna_fl_mass(pl), max_size=flm.max_size)
    ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
    from rigel.calibration.strand_balance import fit_strand_balance
    kfit = float(fit_strand_balance(sm).rna_sense_frac)
    cap = {}
    calmod = sys.modules["rigel.calibration.calibrate"]
    orig = calmod.node_sweep
    calmod.node_sweep = lambda *a, **k: orig(*a, _capture=cap, **k)
    try:
        calibrate(pl, ra, sm, fl.gdna_pmf, fl.rna_pmf, CalibrationConfig())
    finally:
        calmod.node_sweep = orig
    print(f"\n=== {cond}  (fit κ={kfit:.3f}) ===")
    print(f"  enrich_w = {cap.get('enrich_w'):.4f}   rho_global={cap.get('rho_global'):.4f}")
    # per-rtype solved f_g distribution (exon nodes are where the capture leak lives)
    from rigel.calibration.bp_solver import _node_region_type
    chain = build_node_chain(pl.ref_region_offsets, pl.ref_boundary_offsets)
    nrt, _ = _node_region_type(chain, ra)
    kind = np.asarray(chain.kind); fg = np.asarray(cap["f_g"])
    enrich_apply = np.asarray(cap["enrich_apply"])
    for rt, name in [(0, "intergenic"), (1, "intron"), (2, "exon")]:
        m = (kind == REGION) & (nrt == rt) & (np.asarray(cap["mass_l"]) > 0)
        if m.sum():
            print(f"  {name:>10}: n={int(m.sum()):4d}  f_g median={np.median(fg[m]):.3f}  "
                  f"mean={fg[m].mean():.3f}  (enrich_apply on {int((enrich_apply & m).sum())})")


if __name__ == "__main__":
    conds = sys.argv[1:] or [
        "gdna_gdna300_ss_0.50_nrna_none_capture_on",
        "gdna_gdna300_ss_0.99_nrna_none_capture_on",
        "gdna_gdna300_ss_0.50_nrna_none_capture_off",
    ]
    for c in conds:
        probe(c)
