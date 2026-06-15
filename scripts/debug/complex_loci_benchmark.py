"""Complex-locus benchmark SET for the calibration simplex solver.

A library of progressively complex overlapping-opposite-strand loci that stress the deconvolution —
AMBIG (both-strand) regions the per-node 1-D solver cannot resolve locally, which the 2-D simplex
sweep resolves by propagating per-strand RNA:gDNA odds in from strand-informative neighbours. Each
locus is defined in LOCAL coordinates and placed in its own genome window; the harness runs the 1-D
(production) and 2-D (sweep) deconvolutions and scores per-locus mass-weighted AMBIG-node error vs the
oracle. Extend by appending to LOCI; re-run to track the solver improving.

  python scripts/debug/complex_loci_benchmark.py [K_sweep=60] [gdna=120] [nrna=25]
"""
import dataclasses, functools, json, os, sys, tempfile
import numpy as np, pysam
import rigel.calibration.simplex_sweep as ssm
from rigel.sim import Scenario, ReadSimConfig, GDNAConfig
from rigel.pipeline import scan_and_buffer, _native_detect_sj_tag
from rigel.config import PipelineConfig
from rigel.calibration.calibrate import calibrate
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.signature import TS_AMBIG
from rigel.splice import SpliceType
_orig = ssm.deconv_regions_sweep

# The battery is data-driven: a JSON list of loci (name, description, genes[gene,strand,exons,abundance])
# in LOCAL coords [0, WIN). Extend the set by editing complex_loci_battery.json. Falls back to a tiny
# inline set if the JSON is absent.
WIN = 10000
_JSON = os.path.join(os.path.dirname(__file__), "complex_loci_battery.json")
if os.path.exists(_JSON):
    _spec = json.load(open(_JSON))
    LOCI = [(d["name"], d.get("description", ""),
             [(g["gene"], g["strand"], [tuple(e) for e in g["exons"]], g["abundance"]) for g in d["genes"]])
            for d in _spec]
else:
    LOCI = [("cross", "simple +/- overlap",
             [("a+", "+", [(1000, 3000)], 100), ("a-", "-", [(2000, 4000)], 100)])]


def build(work, gdna, nrna):
    sc = Scenario("cxbench", genome_length=(len(LOCI) + 4) * WIN, seed=13, work_dir=work)
    spans = {}
    for li, (name, _desc, genes) in enumerate(LOCI):
        base = (li + 1) * WIN
        spans[name] = (base, base + WIN)
        for gid, strand, exons, ab in genes:
            sc.add_gene(f"{name}_{gid}", strand, [{"t_id": f"{name}_{gid}_t",
                        "exons": [(base + s, base + e) for s, e in exons], "abundance": ab}])
    # strand-training seeds in the trailing window
    sbase = (len(LOCI) + 2) * WIN
    sc.add_gene("seed1", "+", [{"t_id": "S1", "exons": [(sbase, sbase + 500), (sbase + 1500, sbase + 2000),
                (sbase + 3000, sbase + 3500)], "abundance": 150}])
    sc.add_gene("seed2", "-", [{"t_id": "S2", "exons": [(sbase + 5000, sbase + 5500),
                (sbase + 6500, sbase + 7000), (sbase + 8000, sbase + 8500)], "abundance": 150}])
    return sc, spans


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


def main():
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    gdna = float(sys.argv[2]) if len(sys.argv) > 2 else 120.0
    nrna = float(sys.argv[3]) if len(sys.argv) > 3 else 25.0
    with tempfile.TemporaryDirectory() as work:
        sc, spans = build(work, gdna, nrna)
        res = sc.build_oracle(n_fragments=max(80000, len(LOCI) * 5000),
            sim_config=ReadSimConfig(frag_mean=250, frag_std=50, frag_min=80, frag_max=600,
                                     read_length=100, strand_specificity=0.99, seed=13),
            gdna_config=GDNAConfig(abundance=gdna, frag_mean=350, frag_std=100, frag_min=100, frag_max=1000),
            nrna_abundance=nrna)
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

        def est(use_prop, Kg=None):
            # K is now driven by CalibrationConfig.sweep_n_grid (calibrate's hybrid passes it through).
            extra = {"sweep_n_grid": Kg} if Kg else {}
            cfg = dataclasses.replace(
                PipelineConfig().calibration, use_propagation=use_prop, **extra
            )
            r = calibrate(pl, ra, sm, fl.gdna_pmf, fl.rna_pmf, cfg)
            gm = np.asarray(r.mass_gdna_contained); rm = np.asarray(r.mass_rna_contained)
            return np.where(gm + rm > 0, gm / np.maximum(gm + rm, 1e-9), np.nan)

        e1 = est(False); e2 = est(True, K)
        print(f"COMPLEX-LOCI BENCHMARK (gDNA={gdna:.0f} nRNA={nrna:.0f} ss=0.99, sweep K={K})")
        print(f"  mass-weighted AMBIG |f_g − oracle| per locus:\n  {'locus':>12} {'#amb':>4} "
              f"{'1D_err':>9} {'2D_err':>9}   verdict")
        T1 = T2 = 0.0
        amb = sc_arr == TS_AMBIG
        for name, _desc, _g in LOCI:
            lo, hi = spans[name]
            sel = np.flatnonzero(amb & (start >= lo) & (start < hi))
            l1 = l2 = w = 0.0
            for i in sel:
                tot = og[i] + on[i] + om[i]
                if tot < 10 or not np.isfinite(e1[i]) or not np.isfinite(e2[i]):
                    continue
                orc = og[i] / tot
                l1 += abs(e1[i] - orc) * tot; l2 += abs(e2[i] - orc) * tot; w += tot
            if w == 0:
                print(f"  {name:>12} {len(sel):>4}   (no AMBIG mass)"); continue
            T1 += l1; T2 += l2
            v = "2D wins" if l2 < l1 * 0.95 else ("1D wins" if l1 < l2 * 0.95 else "~tie")
            print(f"  {name:>12} {len(sel):>4} {l1:>9,.0f} {l2:>9,.0f}   {v}")
        print(f"  {'TOTAL':>12}      {T1:>9,.0f} {T2:>9,.0f}   "
              f"2D/1D = {T2/max(T1,1):.2f} ({'2D better' if T2 < T1 else '1D better'})")
        sc.cleanup()


if __name__ == "__main__":
    main()
