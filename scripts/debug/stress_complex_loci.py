"""Stress the calibration solver with progressively complex overlapping-opposite-strand loci.

These loci create AMBIG (both-strand) regions that the per-node 1-D solver cannot resolve locally
(the strand is uninformative at an AMBIG node) but the 2-D simplex sweep CAN, by propagating the
per-strand RNA:gDNA odds in from strand-informative single-strand neighbours. Compares, per AMBIG
region, the estimated gDNA fraction of contained-unspliced mass: production (1-D fusion,
use_propagation=False) vs the sweep (2-D propagation, use_propagation=True) vs the oracle.

Loci (well separated; + and − seed genes give strand training + intergenic/intron regions):
  L1 cross      +[5000,7000)  −[6000,8000)                         → +only|AMBIG|−only
  L2 nested     +[10000,16000)  −exons[(11000,12000),(14000,15000)] → AMBIG patches flanked by +only
  L3 interleave +exons[(20000,22000),(24000,26000)] −[(21000,23000),(25000,27000)]
  L4 chain      +[30000,32000) −[31000,34000) +[33000,36000)       → AMBIG|−only|AMBIG chain
  L5 dense      +exons[(40000,41000),(42000,43000),(44000,45000)] −[(40500,41500),(42500,43500),(44500,45500)]
"""
import dataclasses, functools, sys, tempfile
import numpy as np, pysam
import rigel.calibration.simplex_sweep as ssm
from rigel.sim import Scenario, ReadSimConfig, GDNAConfig
from rigel.pipeline import scan_and_buffer, _native_detect_sj_tag
from rigel.config import PipelineConfig
from rigel.calibration.calibrate import calibrate
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.signature import TS_AMBIG, TS_POS, TS_NEG, TS_NONE
from rigel.splice import SpliceType
_orig = ssm.deconv_regions_sweep
SC = {TS_POS: "POS", TS_NEG: "NEG", TS_NONE: "NONE", TS_AMBIG: "AMBIG"}

K_SWEEP = int(sys.argv[1]) if len(sys.argv) > 1 else 60


def make_scenario(work):
    sc = Scenario("stress", genome_length=60000, seed=11, work_dir=work)
    sc.add_gene("A", "+", [{"t_id": "A1", "exons": [(5000, 7000)], "abundance": 100}])
    sc.add_gene("B", "-", [{"t_id": "B1", "exons": [(6000, 8000)], "abundance": 100}])
    sc.add_gene("C", "+", [{"t_id": "C1", "exons": [(10000, 16000)], "abundance": 100}])
    sc.add_gene("D", "-", [{"t_id": "D1", "exons": [(11000, 12000), (14000, 15000)], "abundance": 100}])
    sc.add_gene("E", "+", [{"t_id": "E1", "exons": [(20000, 22000), (24000, 26000)], "abundance": 100}])
    sc.add_gene("F", "-", [{"t_id": "F1", "exons": [(21000, 23000), (25000, 27000)], "abundance": 100}])
    sc.add_gene("G", "+", [{"t_id": "G1", "exons": [(30000, 32000)], "abundance": 100}])
    sc.add_gene("H", "-", [{"t_id": "H1", "exons": [(31000, 34000)], "abundance": 100}])
    sc.add_gene("I", "+", [{"t_id": "I1", "exons": [(33000, 36000)], "abundance": 100}])
    sc.add_gene("J", "+", [{"t_id": "J1", "exons": [(40000, 41000), (42000, 43000), (44000, 45000)],
                            "abundance": 100}])
    sc.add_gene("K", "-", [{"t_id": "K1", "exons": [(40500, 41500), (42500, 43500), (44500, 45500)],
                            "abundance": 100}])
    # strand-training seeds (well separated, both strands, multi-exon → spliced reads)
    sc.add_gene("s1", "+", [{"t_id": "S1", "exons": [(50000, 50500), (51500, 52000), (53000, 53500)],
                             "abundance": 120}])
    sc.add_gene("s2", "-", [{"t_id": "S2", "exons": [(55000, 55500), (56500, 57000), (58000, 58500)],
                             "abundance": 120}])
    return sc


def oracle(bam, starts, ids, R):
    g, n, m = np.zeros(R), np.zeros(R), np.zeros(R)
    with pysam.AlignmentFile(bam, "rb") as b:
        for r in b.fetch(until_eof=True):
            if r.is_secondary or r.is_supplementary or r.is_unmapped or not r.is_read1:
                continue
            if r.cigartuples and any(op == 3 for op, _ in r.cigartuples):
                continue
            ref = r.reference_name
            if ref not in starts:
                continue
            i = int(np.searchsorted(starts[ref], r.reference_start, side="right") - 1)
            if i < 0:
                continue
            rid = int(ids[ref][i]); nm = r.query_name
            (g if nm.startswith("gdna") else n if nm.startswith("nrna") else m)[rid] += 1
    return g, n, m


with tempfile.TemporaryDirectory() as work:
    sc = make_scenario(work)
    res = sc.build_oracle(
        n_fragments=60000,
        sim_config=ReadSimConfig(frag_mean=250, frag_std=50, frag_min=80, frag_max=600,
                                 read_length=100, strand_specificity=0.99, seed=11),
        gdna_config=GDNAConfig(abundance=120, frag_mean=350, frag_std=100, frag_min=100, frag_max=1000),
        nrna_abundance=25.0,
    )
    idx = res.index; bam = str(res.bam_path); df = idx.region_df
    starts = {r: gg["start"].to_numpy() for r, gg in df.groupby("ref_name")}
    ids = {r: gg["region_id"].to_numpy() for r, gg in df.groupby("ref_name")}
    R = len(df)
    scan = dataclasses.replace(PipelineConfig().scan, sj_strand_tag=_native_detect_sj_tag(bam))
    st, sm, flm, buf, pl = scan_and_buffer(bam, idx, scan)
    ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
    fl = build_fl_models(global_counts=flm.global_model.counts,
                         rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                         gdna_counts=gdna_fl_mass(pl), max_size=flm.max_size)
    og, on, om = oracle(bam, starts, ids, R)
    sc_arr = np.asarray(ra.strand_class); start = df["start"].to_numpy(); end = df["end"].to_numpy()

    def est(use_prop, K=None):
        ssm.deconv_regions_sweep = functools.partial(_orig, n_grid=K) if K else _orig
        cfg = dataclasses.replace(PipelineConfig().calibration, use_propagation=use_prop)
        try:
            r = calibrate(pl, ra, sm, fl.gdna_pmf, fl.rna_pmf, cfg)
        finally:
            ssm.deconv_regions_sweep = _orig
        gm = np.asarray(r.mass_gdna_contained); rm = np.asarray(r.mass_rna_contained)
        return np.where(gm + rm > 0, gm / np.maximum(gm + rm, 1e-9), np.nan)

    e_prod = est(False)
    e_sweep = est(True, K_SWEEP)
    amb = np.flatnonzero(sc_arr == TS_AMBIG)
    print(f"complex stress loci — AMBIG regions (gDNA=120, nRNA=25, ss=0.99); sweep K={K_SWEEP}")
    print(f"  {'region':>14} {'cls':>5} {'g':>4} {'n':>4} {'m':>4} {'orcl_g':>7} {'2term':>6} "
          f"| {'PROD(1D)':>8} {'SWEEP(2D)':>9} {'1D_err':>7} {'2D_err':>7}")
    tot1 = tot2 = 0.0
    for i in amb:
        tot = og[i] + on[i] + om[i]
        if tot < 10:
            continue
        orc = og[i] / tot
        tt = (og[i] + on[i]) / tot
        e1 = abs(e_prod[i] - orc); e2 = abs(e_sweep[i] - orc)
        tot1 += e1 * tot; tot2 += e2 * tot
        print(f"  {f'{start[i]}-{end[i]}':>14} {SC[int(sc_arr[i])]:>5} {og[i]:>4.0f} {on[i]:>4.0f} "
              f"{om[i]:>4.0f} {orc:>7.3f} {tt:>6.3f} | {e_prod[i]:>8.3f} {e_sweep[i]:>9.3f} "
              f"{e1:>7.3f} {e2:>7.3f}")
    print(f"\n  mass-weighted AMBIG |err|:  PROD(1D)={tot1:,.0f}   SWEEP(2D)={tot2:,.0f}   "
          f"({'2D better' if tot2 < tot1 else '1D better'})")
    sc.cleanup()
