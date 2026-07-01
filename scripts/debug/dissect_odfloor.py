"""Deep dissection of the gDNA=0 + nascent AMBIG false-gDNA (test_ambig_no_false_gdna_from_nascent).

Reproduces the exact test scenario and traces the failure end-to-end:
  * PASS 1 (no KDE): the AMBIG node f_g, and WHY the pure-RNA nodes solve to f_g≈0.02 (strand-only vs
    local vs the messages) — the suspected od_r floor that seeds the KDE's RNA-residual mode.
  * PASS 2 (KDE): the AMBIG node's decomposition — strand-only / local(+KDE) / final, the messages, and the
    KDE prior term evaluated at the node's grid — to see exactly what pulls it to 0.376.
"""
from __future__ import annotations

import dataclasses
import sys

import numpy as np

from rigel.calibration.calibrate import calibrate
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.signature import TS_AMBIG
from rigel.config import PipelineConfig
from rigel.pipeline import _native_detect_sj_tag, scan_and_buffer
from rigel.sim import ReadSimConfig, Scenario
from rigel.splice import SpliceType
from rigel.calibration.node_chain import build_node_chain, REGION
import tempfile


def _scenario():
    sc = Scenario("ambig_reg", genome_length=30000, seed=7, work_dir=tempfile.mkdtemp())
    sc.add_gene("gA", "+", [{"t_id": "TA", "exons": [(1000, 1500), (4000, 6000)], "abundance": 100}])
    sc.add_gene("gB", "-", [{"t_id": "TB", "exons": [(5000, 7000), (10000, 10500)], "abundance": 100}])
    sc.add_gene("s1", "+", [{"t_id": "S1", "exons": [(12000, 12500), (13500, 14000), (15000, 15500)],
                             "abundance": 120}])
    sc.add_gene("s2", "-", [{"t_id": "S2", "exons": [(17000, 17500), (18500, 19000), (20000, 20500)],
                             "abundance": 120}])
    res = sc.build_oracle(
        n_fragments=8000,
        sim_config=ReadSimConfig(frag_mean=250, frag_std=50, frag_min=80, frag_max=600,
                                 read_length=100, strand_specificity=0.99, seed=7),
        gdna_config=None, nrna_abundance=30.0)
    idx, bam = res.index, str(res.bam_path)
    cfg = PipelineConfig()
    scan = dataclasses.replace(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    _st, sm, flm, _buf, pl = scan_and_buffer(bam, idx, scan)
    ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
    fl = build_fl_models(global_counts=flm.global_model.counts,
                         rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                         gdna_counts=gdna_fl_mass(pl), max_size=flm.max_size)
    return pl, ra, sm, fl, cfg


def _run(pl, ra, sm, fl, cfg, *, kde: bool):
    caps = []
    cm = sys.modules["rigel.calibration.calibrate"]
    orig = cm.node_sweep
    def wrap(*a, **k):
        cap = {}
        r = orig(*a, _capture=cap, **k)
        caps.append(cap)
        return r
    cm.node_sweep = wrap
    cal_cfg = dataclasses.replace(cfg.calibration, gdna_prior_enable=kde)
    try:
        res = calibrate(pl, ra, sm, fl.gdna_pmf, fl.rna_pmf, cal_cfg)
    finally:
        cm.node_sweep = orig
    return res, caps


def main():
    pl, ra, sm, fl, cfg = _scenario()
    chain = build_node_chain(pl.ref_region_offsets, pl.ref_boundary_offsets)
    kind = np.asarray(chain.kind); ridx = np.asarray(chain.ref_idx)
    sc_arr = np.asarray(ra.strand_class)
    ambig_regions = np.flatnonzero(sc_arr == TS_AMBIG)
    # chain node(s) for the AMBIG region(s)
    ambig_nodes = [int(n) for n in np.where(kind == REGION)[0] if ridx[n] in ambig_regions]
    print(f"AMBIG region idx={list(ambig_regions)}  chain node(s)={ambig_nodes}")

    res1, caps1 = _run(pl, ra, sm, fl, cfg, kde=False)
    res2, caps2 = _run(pl, ra, sm, fl, cfg, kde=True)

    def ambig_frac(res):
        g = np.asarray(res.mass_gdna_contained)[ambig_regions]
        r = np.asarray(res.mass_rna_contained)[ambig_regions]
        return float(g.sum() / max(g.sum() + r.sum(), 1e-9))
    print(f"\nAMBIG contained gDNA fraction: PASS1(no KDE)={ambig_frac(res1):.3f}  PASS2(KDE)={ambig_frac(res2):.3f}")

    cap1 = caps1[0]  # pass1 (only one sweep)
    cap2 = caps2[-1]  # pass2 (the LAST sweep = the KDE re-solve)
    nd = ambig_nodes[0]
    fgs = lambda c, k: float(np.asarray(c[k])[nd])
    print("\n--- AMBIG node decomposition ---")
    for tag, c in (("PASS1", cap1), ("PASS2", cap2)):
        print(f"  {tag}: fg_strand={fgs(c,'fg_strand'):.3f}  fg_loc={fgs(c,'fg_loc'):.3f}  "
              f"fg_final={fgs(c,'f_g'):.3f}  | gMSG(mode={fgs(c,'mode_g'):+.2f},prec={fgs(c,'prec_g'):.2f})  "
              f"pMSG(mode={fgs(c,'mode_p')+fgs(c,'mode_n'):+.2f},prec={fgs(c,'prec_p')+fgs(c,'prec_n'):.2f})")

    # WHY do the pure-RNA nodes solve to f_g≈0.02? strand-only fg for every solved REGION node.
    print("\n--- pure-RNA region nodes: strand-only vs local f_g (PASS1) ---")
    fg_strand = np.asarray(cap1["fg_strand"]); fg_loc = np.asarray(cap1["fg_loc"]); fg_fin = np.asarray(cap1["f_g"])
    solv = np.asarray(cap1["solvable"])
    for n in np.where((kind == REGION) & solv)[0][:14]:
        print(f"  node {int(n):3d} ref{int(ridx[n]):3d}: fg_strand={fg_strand[n]:.3f} "
              f"fg_loc={fg_loc[n]:.3f} fg_final={fg_fin[n]:.3f}")


if __name__ == "__main__":
    main()
