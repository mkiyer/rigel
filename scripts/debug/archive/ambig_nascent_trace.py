"""Interrogate the failing test_ambig_no_false_gdna_from_nascent (gDNA=0 + nascent, ss0.99).

The phantom gDNA appears at BOTH the AMBIG exon AND the single-strand POS/NEG exons in a ZERO-gDNA
library. Ablation of the sweep terms (this script) localises the cause:

  - ê(z) enrichment is NOT firing (enrich_sig=False, ρ_global≈0.01) — refuted as the driver.
  - The strand-only init is CORRECT: single-strand exon f_g=0.000, AMBIG f_g=1.0 (the {0,0,1} default).
  - Killing the per-strand RNA IMPUTATION message drops the single-strand-exon phantom hard
    (ALL.exon 0.38→0.17): the RNA message pulls a confident single-strand exon's f_pos DOWN toward its
    lower-density intron/boundary neighbour (mature RNA is concentrated in the exon and is NOT smooth),
    manufacturing gDNA. Killing the gDNA message makes it WORSE (the gDNA msg was helping).
  - The AMBIG node sits at the uniform-simplex marginal (~0.38): balanced counts are strand-indistinct
    from gDNA, the gentle global (N≈1 at zero-gDNA) can't pin it, and its single-strand flanks — which
    SHOULD impute "low gDNA" — are themselves RNA-message-inflated → a dirty cascade.

    OMP_NUM_THREADS=1 python scripts/debug/ambig_nascent_trace.py
"""
from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
import dataclasses
import tempfile

import numpy as np

import rigel.calibration.bp_solver as bp
from rigel.calibration.calibrate import calibrate
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.signature import TS_AMBIG, coarse_type_from_signature
from rigel.config import PipelineConfig
from rigel.pipeline import _native_detect_sj_tag, scan_and_buffer
from rigel.sim import ReadSimConfig, Scenario
from rigel.splice import SpliceType


def build(work):
    sc = Scenario("ambig_reg", genome_length=30000, seed=7, work_dir=work)
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
        gdna_config=None, nrna_abundance=30.0,
    )
    cfg = PipelineConfig()
    scan = dataclasses.replace(cfg.scan, sj_strand_tag=_native_detect_sj_tag(str(res.bam_path)))
    _st, sm, flm, _buf, pl = scan_and_buffer(str(res.bam_path), res.index, scan)
    ra = RegionArrays.from_region_df(res.index.region_df, res.index.ref_name_to_id)
    fl = build_fl_models(global_counts=flm.global_model.counts,
                         rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                         gdna_counts=gdna_fl_mass(pl), max_size=flm.max_size)
    return sc, pl, ra, sm, fl, cfg


def ambig_frac(result, ra):
    sig = np.asarray(ra.signature)
    sc_arr = np.asarray(ra.strand_class)
    g = np.asarray(result.mass_gdna_contained)
    r = np.asarray(result.mass_rna_contained)
    exon = np.array([coarse_type_from_signature(int(s)) == 2 for s in sig])
    out = {}
    for cls, lab in [(TS_AMBIG, "AMBIG.exon")]:
        m = (sc_arr == cls) & exon
        gg, rr = g[m].sum(), r[m].sum()
        out[lab] = gg / max(gg + rr, 1e-9)
    # all exon phantom (gDNA=0 library ⇒ every gDNA call is phantom)
    out["ALL.exon"] = g[exon].sum() / max(g[exon].sum() + r[exon].sum(), 1e-9)
    return out


def _run(pl, ra, sm, fl, cfg, *, kill_messages=False, kill_global=False,
         kill_rna_msg=False, kill_gdna_msg=False, force_odg_eq_odr=False, force_odg_zero=False,
         scale_global=1.0):
    """Calibrate with optional ablation of imputation messages, the global prior, or the gDNA strand od."""
    import rigel.calibration.simplex_sweep as ss
    orig_msg = bp._message
    orig_glp = bp._global_logprior
    orig_ll = bp._local_loglik
    orig_mix = ss._mixture_strand_loglik
    if kill_messages:
        bp._message = lambda *a, **k: (0.0, 0.0)
    if kill_global:
        bp._global_logprior = lambda fgg, *a, **k: np.zeros((np.asarray(a[0]).shape[0], fgg.shape[0]))
    elif scale_global != 1.0:
        bp._global_logprior = lambda *a, **k: scale_global * orig_glp(*a, **k)  # scale the global pseudocount N
    if kill_rna_msg or kill_gdna_msg:
        def ll_wrap(*a, **k):
            if kill_rna_msg:
                k["rna_imp_frac"] = None
                k["rna_imp_count"] = None
            if kill_gdna_msg:
                k["gdna_imp_frac"] = None
                k["gdna_imp_count"] = None
            return orig_ll(*a, **k)
        bp._local_loglik = ll_wrap
    if force_odg_eq_odr or force_odg_zero:
        def mix_wrap(u_pos, n, f_g, f_pos, f_neg, kappa, od_g, od_r):
            new_odg = 0.0 if force_odg_zero else od_r
            return orig_mix(u_pos, n, f_g, f_pos, f_neg, kappa, new_odg, od_r)
        ss._mixture_strand_loglik = mix_wrap
    try:
        return calibrate(pl, ra, sm, fl.gdna_pmf, fl.rna_pmf, cfg.calibration)
    finally:
        bp._message = orig_msg
        bp._global_logprior = orig_glp
        bp._local_loglik = orig_ll
        ss._mixture_strand_loglik = orig_mix


def precision_report(pl, ra, sm, fl, cfg):
    """Capture the EXACT precision values at the single-strand POS exon (u_pos≈5, u_neg≈622) on the
    converged pass, and watch f_g move as each ψ term is added: strand → +global → +gDNA msg → +RNA msg."""
    from rigel.calibration.simplex_sweep import _fg_median, _fg_var, _local_loglik, _simplex_lattice

    cap = {}
    orig_ll = bp._local_loglik

    def ll_spy(u_pos, u_neg, *a, **k):
        n = float(u_pos[0]) + float(u_neg[0])
        # single-strand POS exon: ~627 unspliced, +strand allowed, −strand forbidden.
        if 600 <= n <= 700 and bool(a[2][0]) and not bool(a[3][0]):  # a=(spl_pos,spl_neg,allow_pos,allow_neg,...)
            cap["u_pos"] = float(u_pos[0])
            cap["u_neg"] = float(u_neg[0])
            cap["args"] = (u_pos.copy(), u_neg.copy(), a)
            cap["kw"] = dict(k)
        return orig_ll(u_pos, u_neg, *a, **k)

    bp._local_loglik = ll_spy
    calibrate(pl, ra, sm, fl.gdna_pmf, fl.rna_pmf, cfg.calibration)
    bp._local_loglik = orig_ll
    if "args" not in cap:
        print("  (target single-strand exon node not captured)")
        return

    u_pos, u_neg, a = cap["args"]
    kw = cap["kw"]
    n_grid = cfg.calibration.sweep_n_grid
    lat = _simplex_lattice(int(n_grid))
    fgg = lat[2]
    kappa = a[4]  # (spl_pos, spl_neg, allow_pos, allow_neg, kappa, od_g, od_r, lattice)
    odg, odr = a[5], a[6]
    mu_g = float(kw["gdna_imp_frac"][0])
    n_g = float(kw["gdna_imp_count"][0])
    mu_p = float(kw["rna_imp_frac"][0][0])
    n_p = float(kw["rna_imp_count"][0][0])

    def solve(**terms):
        psi = _local_loglik(u_pos, u_neg, a[0], a[1], a[2], a[3], kappa, odg, odr, lat,
                            strand_obs=kw.get("strand_obs"), **terms)
        return float(_fg_median(psi, fgg)[0]), float(_fg_var(psi, fgg)[0])

    fg_s, var_s = solve()
    fg_sg, _ = solve(global_logprior=kw["global_logprior"])
    fg_sgg, _ = solve(global_logprior=kw["global_logprior"],
                      gdna_imp_frac=kw["gdna_imp_frac"], gdna_imp_count=kw["gdna_imp_count"])
    fg_full, _ = solve(global_logprior=kw["global_logprior"],
                       gdna_imp_frac=kw["gdna_imp_frac"], gdna_imp_count=kw["gdna_imp_count"],
                       rna_imp_frac=kw["rna_imp_frac"], rna_imp_count=kw["rna_imp_count"])

    N = cap["u_pos"] + cap["u_neg"]
    contrast = (2 * kappa - 1) ** 2
    print("\n=== PRECISION at the single-strand POS exon (the confident-RNA node) ===")
    print(f"  counts: u_pos={cap['u_pos']:.0f} u_neg={cap['u_neg']:.0f}  N={N:.0f}")
    print(f"  strand: κ={kappa:.4f}  (2κ−1)²={contrast:.3f}  ⇒ strand N_eff = N·(2κ−1)² = {N*contrast:.0f}")
    print(f"  overdispersion: od_g={odg:.4g}  od_r={odr:.4g}")
    # the strand likelihood var (Gaussian) at f_g=0 vs f_g=0.167 — shows whether overdispersion flattens it
    rscale = kappa * (1.0 - kappa)
    for fgt in (0.0, 0.167):
        fp_t = 1.0 - fgt
        p = 0.5 * fgt + kappa * fp_t
        var = N * p * (1 - p) + (N * fgt) ** 2 * 0.25 * odg + (N * fp_t) ** 2 * rscale * odr
        pen = (cap["u_pos"] - N * p) ** 2 / max(var, 1e-9)
        print(f"    f_g={fgt:.3f}: pred mean(u_pos)={N*p:.1f} (obs {cap['u_pos']:.0f})  var={var:.1f}  "
              f"Gaussian penalty={pen:.2f}")
    print(f"          strand-only posterior var(f_g)={var_s:.4g}  ⇒ precision 1/var = {1/max(var_s,1e-12):.0f}")
    print(f"  gDNA msg:  μ_g={mu_g:.3f}  N_g={n_g:.2f}")
    print(f"  RNA+ msg:  μ_p={mu_p:.3f}  N_p={n_p:.2f}   wall-force N_p·(1−μ_p)={n_p*(1-mu_p):.2f}")
    print("             (μ_p<1 ⇒ binom pseudo-count adds (N_p·(1−μ_p))·log(f_g) → −∞ wall pushing f_g UP)")
    print("  f_g median as terms accumulate:")
    print(f"      strand only           : {fg_s:.3f}")
    print(f"      + global              : {fg_sg:.3f}")
    print(f"      + gDNA msg            : {fg_sgg:.3f}")
    print(f"      + RNA msg (= full)    : {fg_full:.3f}")


def main():
    with tempfile.TemporaryDirectory() as work:
        sc, pl, ra, sm, fl, cfg = build(work)

        configs = [
            ("full sweep", dict()),
            ("no MESSAGES", dict(kill_messages=True)),
            ("no GLOBAL", dict(kill_global=True)),
            ("no RNA-msg", dict(kill_rna_msg=True)),
            ("no gDNA-msg", dict(kill_gdna_msg=True)),
            ("no msg+glob", dict(kill_messages=True, kill_global=True)),
            ("od_g:=od_r", dict(force_odg_eq_odr=True)),
            ("od_g:=0", dict(force_odg_zero=True)),
            ("od_g:=od_r+nomsg", dict(force_odg_eq_odr=True, kill_messages=True)),
            ("global×10", dict(scale_global=10.0)),
            ("global×30", dict(scale_global=30.0)),
            ("global×30+od_g:=od_r", dict(scale_global=30.0, force_odg_eq_odr=True)),
        ]
        rows = {}
        kappa = None
        for label, kw in configs:
            res = _run(pl, ra, sm, fl, cfg, **kw)
            kappa = res.rna_sense_frac
            rows[label] = ambig_frac(res, ra)

        print(f"\n=== gDNA=0 + nascent=30, ss0.99 (κ={kappa:.4f}) — ablation of the sweep terms ===")
        print("  (strand-only init: single-strand exon f_g=0.000 CORRECT, AMBIG f_g=1.0 default)")
        print(f"  {'config':>14} {'AMBIG.exon':>11} {'ALL.exon':>9}   (test wants AMBIG < 0.08)")
        for label, _ in configs:
            print(f"  {label:>14} {rows[label]['AMBIG.exon']:>11.3f} {rows[label]['ALL.exon']:>9.3f}")

        precision_report(pl, ra, sm, fl, cfg)
        sc.cleanup()


if __name__ == "__main__":
    main()
