"""Monte-Carlo + fold validation of the HYBRID enrichment-corrected density message mode
(docs/CARRY_FORWARD.md §9) against the two predecessors it unifies.

THE CLAIM. A region↔boundary message provides two soft targets to the dst's λ-fold: a gDNA factor on
``log f_g`` and an RNA-total factor on ``log f_r``. Write ``M_c^pred = ρ_c^src·E_c^dst`` (source density in the
dst frame), ``_den = ΣM_c^pred`` (imputed total), ``md`` = observed dst total, ``logR = log(ρ_tot^dst/ρ_tot^src)``
(the enrichment jump the NPMLE projection measures as ``mu_proj[dst]−mu_proj[src]``). The three modes:

    density      t_c = log(M_c^pred / md)                 anchored to observed md; Σexp(t)=_den/md ≠ 1 at a cliff
    composition  t_c = log(M_c^pred / _den)               Σexp(t)=1 (consistent) but total is the imputed _den
    HYBRID       t_c = log(M_c^pred / md) + logR           = density + logR; Σexp(t)=R·_den/md → 1 iff src accurate

The log-ODDS ``t_g − t_r = log(M_g/M_r)`` is the cliff-invariant eff-length shift for BOTH composition and hybrid
(the ``+logR`` cancels in the difference). The difference is the CONSISTENCY of the pair, which governs how the
fold weights the observed-md anchor:
  * accurate source  → hybrid Σexp=1 → behaves like composition (confident, cliff-correct);
  * inaccurate source → hybrid Σexp≠1 → the fold compromises toward the observed md (regularized), while
    composition sits confidently at the WRONG transported λ.

We PROVE both by simulating a genuine cliff (e_src ≠ e_dst), feeding each mode's factors through the REAL fold
(`bp_solver._fold_lambda`), and comparing the folded f_g to ground truth — across FL pairs, cliff magnitudes,
cliff directions, and injected source-belief errors. Run in the rigel env (OMP_NUM_THREADS=1).
"""

from __future__ import annotations

import math

import numpy as np

from rigel.calibration.bp_solver import _fold_lambda
from rigel.calibration.effective_length import (
    boundary_side_eff_length,
    region_eff_length,
)


def make_pmf(kind: str, m: float, s: float, L: int = 900) -> np.ndarray:
    x = np.arange(1, L)
    if kind == "gauss":
        p = np.exp(-0.5 * ((x - m) / s) ** 2)
    elif kind == "gamma":
        k = (m / s) ** 2
        p = x ** (k - 1) * np.exp(-x / (s * s / m))
    elif kind == "unif":
        p = ((x >= m - s) & (x <= m + s)).astype(float)
    elif kind == "bimodal":
        p = np.exp(-0.5 * ((x - m) / s) ** 2) + 0.7 * np.exp(-0.5 * ((x - 2.2 * m) / s) ** 2)
    else:
        raise ValueError(kind)
    p = np.maximum(p, 0.0)
    pmf = np.zeros(L)
    pmf[1:] = p / p.sum()
    return pmf


def _fold_fg(t_g, pr_g, t_r, pr_r, *, mu0=0.0, var0=100.0):
    """Run the REAL two-factor λ-fold from a prior N(mu0,var0) and return (f_g, folded_var). The message alone
    determines the outcome when var0 large, isolating the MODE arithmetic exactly as it enters the solver."""
    factors = []
    if pr_g > 0:
        factors.append((True, float(t_g), float(pr_g)))
    if pr_r > 0:
        factors.append((False, float(t_r), float(pr_r)))
    mu, var = _fold_lambda(
        mu0, var0, factors, L=10.0, coarse_k=33, fine_k=33, sigma_cov=6.0, refine=3
    )
    return 1.0 / (1.0 + math.exp(-mu)), var


def modes(rho_g, rho_r, EgD, ErD, md, logR):
    """The three (t_g, t_r) target pairs. rho_* are the SOURCE densities (belief-derived); md observed dst."""
    Mg = rho_g * EgD
    Mr = rho_r * ErD
    den = max(Mg + Mr, 1e-12)
    lg, lr = math.log(max(Mg, 1e-12)), math.log(max(Mr, 1e-12))
    lmd = math.log(md)
    return {
        "density": (lg - lmd, lr - lmd),
        "composition": (lg - math.log(den), lr - math.log(den)),
        "hybrid": (lg - lmd + logR, lr - lmd + logR),
    }


def run(gpmf, rpmf, *, Ls=3000.0, Rfl=3000.0, d_g=0.4, d_n=0.6, e_src=1.0, e_dst=1.0, belief_err=0.0,
        r_use="est", pr_g=8.0, pr_r=8.0):
    """One edge: intron REGION (contained, e_src) -> intron-exon BOUNDARY (crossing, e_dst), clean (gDNA+nascent).
    r_use: 'true' = exact e_dst/e_src; 'est' = the observed total-density ratio (what the NPMLE denoises)."""
    EgS = float(region_eff_length(np.array([Ls]), gpmf)[0])
    ErS = float(region_eff_length(np.array([Ls]), rpmf)[0])
    EgD = float(boundary_side_eff_length(gpmf, np.array([Rfl]))[0])
    ErD = float(boundary_side_eff_length(rpmf, np.array([Rfl]))[0])

    # observed masses (enrichment scales the observed density)
    MgS, MrS = e_src * d_g * EgS, e_src * d_n * ErS
    sm = MgS + MrS
    MgD, MrD = e_dst * d_g * EgD, e_dst * d_n * ErD
    md = MgD + MrD
    fg_dst_true = MgD / md

    # source belief (optionally corrupted), and its imputed densities
    fg_src_true = MgS / sm
    lam_src = math.log(fg_src_true / (1 - fg_src_true)) + belief_err
    fg_src_b = 1.0 / (1.0 + math.exp(-lam_src))
    rho_g = fg_src_b * sm / EgS
    rho_r = (1.0 - fg_src_b) * sm / ErS

    # enrichment ratio logR: the NPMLE projects log(mass/eff_g). For a clean edge the composition is uniform so
    # the observed total-density ratio (on a COMMON gDNA-eff basis) recovers log(e_dst/e_src) exactly.
    if r_use == "true":
        logR = math.log(e_dst / e_src)
    else:
        logR = math.log((md / EgD) / (sm / EgS))

    out, var = {}, {}
    for name, (t_g, t_r) in modes(rho_g, rho_r, EgD, ErD, md, logR).items():
        fg, fv = _fold_fg(t_g, pr_g, t_r, pr_r)
        out[name] = abs(fg - fg_dst_true)
        var[name] = fv
    return out, fg_dst_true, var


def run_corrupt_local(gpmf, rpmf, *, e_dst=100.0, belief_err=1.0, r_use="est", local_pr=6.0, msg_pr=8.0):
    """Test E: the dst already has a GOOD local belief at the TRUE f_g (precision local_pr). A WRONG message
    (source error) folds in at precision msg_pr. Which mode corrupts the good local belief LEAST? This mirrors
    the sim: exons have their own strand/global belief, and a confident-wrong exon-source message overrides it."""
    EgS = float(region_eff_length(np.array([3000.0]), gpmf)[0])
    ErS = float(region_eff_length(np.array([3000.0]), rpmf)[0])
    EgD = float(boundary_side_eff_length(gpmf, np.array([3000.0]))[0])
    ErD = float(boundary_side_eff_length(rpmf, np.array([3000.0]))[0])
    d_g, d_n = 0.4, 0.6
    MgS, MrS = d_g * EgS, d_n * ErS
    sm = MgS + MrS
    MgD, MrD = e_dst * d_g * EgD, e_dst * d_n * ErD
    md = MgD + MrD
    fg_dst_true = MgD / md
    lam_dst_true = math.log(fg_dst_true / (1 - fg_dst_true))
    fg_src_true = MgS / sm
    lam_src = math.log(fg_src_true / (1 - fg_src_true)) + belief_err
    fg_src_b = 1.0 / (1.0 + math.exp(-lam_src))
    rho_g = fg_src_b * sm / EgS
    rho_r = (1.0 - fg_src_b) * sm / ErS
    logR = math.log(e_dst) if r_use == "true" else math.log((md / EgD) / (sm / EgS))
    res = {}
    for name, (t_g, t_r) in modes(rho_g, rho_r, EgD, ErD, md, logR).items():
        # start from the good local belief N(lam_dst_true, 1/local_pr), fold the wrong message
        fg, _ = _fold_fg(t_g, msg_pr, t_r, msg_pr, mu0=lam_dst_true, var0=1.0 / local_pr)
        res[name] = abs(fg - fg_dst_true)
    return res


def main():
    fls = [
        ("gDNA=RNA g200", make_pmf("gauss", 200, 40), make_pmf("gauss", 200, 40)),
        ("gDNA300/RNA150", make_pmf("gauss", 300, 50), make_pmf("gauss", 150, 30)),
        ("gDNA120/RNA350", make_pmf("gauss", 120, 25), make_pmf("gauss", 350, 60)),
        ("gamma300/180", make_pmf("gamma", 300, 90), make_pmf("gamma", 180, 60)),
        ("unif250/150", make_pmf("unif", 250, 100), make_pmf("unif", 150, 60)),
        ("bimodal/g220", make_pmf("bimodal", 150, 30), make_pmf("gauss", 220, 45)),
    ]
    for r_use in ("true", "est"):
        print("=" * 100)
        print(f"A[R={r_use}]. ACCURATE SOURCE across a 100x cliff. |Δf_g| vs oracle (diffuse prior = mode only):")
        print(f"{'FL pair':16s} {'fg_true':>8s} {'density':>9s} {'composit':>9s} {'HYBRID':>9s}")
        agg = {"density": [], "composition": [], "hybrid": []}
        for name, g, r in fls:
            o, ft, _ = run(g, r, e_src=1.0, e_dst=100.0, belief_err=0.0, r_use=r_use)
            for k in agg:
                agg[k].append(o[k])
            print(f"{name:16s} {ft:8.4f} {o['density']:9.4f} {o['composition']:9.4f} {o['hybrid']:9.4f}")
        print(f"{'MEAN':16s} {'':8s} {np.mean(agg['density']):9.4f} {np.mean(agg['composition']):9.4f} "
              f"{np.mean(agg['hybrid']):9.4f}")

    print("\n" + "=" * 100)
    print("B. CLIFF-INVARIANCE (accurate source, gDNA300/RNA150), R=true vs est: |Δf_g| at each cliff:")
    g, r = make_pmf("gauss", 300, 50), make_pmf("gauss", 150, 30)
    print(f"{'e_dst/e_src':>12s} | {'density':>9s} {'comp':>7s} {'HYB(true)':>9s} {'HYB(est)':>9s}")
    for e_dst in (1.0, 10.0, 100.0, 1000.0, 0.1, 0.01):
        ot, _, _ = run(g, r, e_src=1.0, e_dst=e_dst, belief_err=0.0, r_use="true")
        oe, _, _ = run(g, r, e_src=1.0, e_dst=e_dst, belief_err=0.0, r_use="est")
        print(f"{e_dst:12.2f} | {ot['density']:9.4f} {ot['composition']:7.4f} {ot['hybrid']:9.4f} "
              f"{oe['hybrid']:9.4f}")

    print("\n" + "=" * 100)
    print("C. SOURCE-BELIEF ERROR, gDNA300/RNA150, 100x cliff, R=true. |Δf_g| and folded VARIANCE (confidence):")
    print(f"{'belief_err':>10s} | {'dens Δ':>7s} {'comp Δ':>7s} {'HYB Δ':>7s} | {'dens σ²':>8s} {'comp σ²':>8s} {'HYB σ²':>8s}")
    for be in (-2.0, -1.0, 0.0, 1.0, 2.0):
        o, _, v = run(g, r, e_src=1.0, e_dst=100.0, belief_err=be, r_use="true")
        print(f"{be:10.2f} | {o['density']:7.4f} {o['composition']:7.4f} {o['hybrid']:7.4f} | "
              f"{v['density']:8.4f} {v['composition']:8.4f} {v['hybrid']:8.4f}")

    print("\n" + "=" * 100)
    print("D. SOURCE ERROR averaged over FL pairs + cliffs + error levels (R=true), diffuse prior:")
    for r_use in ("true", "est"):
        agg = {"density": [], "composition": [], "hybrid": []}
        for name, g, r in fls:
            for e_dst in (1.0, 30.0, 100.0, 300.0):
                for be in (-1.5, -0.75, 0.75, 1.5):
                    o, _, _ = run(g, r, e_src=1.0, e_dst=e_dst, belief_err=be, r_use=r_use)
                    for k in agg:
                        agg[k].append(o[k])
        print(f"  R={r_use:4s}: density={np.mean(agg['density']):.4f}  composition={np.mean(agg['composition']):.4f}"
              f"  HYBRID={np.mean(agg['hybrid']):.4f}   (mean |Δf_g| over {len(agg['hybrid'])} cases)")

    print("\n" + "=" * 100)
    print("E. CORRUPT-A-GOOD-LOCAL-BELIEF (the sim mechanism). dst has a GOOD local belief at the TRUE f_g;")
    print("   a WRONG source message folds in. Lower = protects the good belief better. R=true, msg_pr=8, local_pr=6:")
    print(f"{'FL pair':16s} {'be':>5s} | {'density':>9s} {'composit':>9s} {'HYBRID':>9s}")
    agg = {"density": [], "composition": [], "hybrid": []}
    for name, g, r in fls:
        for be in (1.0, 2.0):
            res = run_corrupt_local(g, r, e_dst=100.0, belief_err=be, r_use="true")
            for k in agg:
                agg[k].append(res[k])
            print(f"{name:16s} {be:5.1f} | {res['density']:9.4f} {res['composition']:9.4f} {res['hybrid']:9.4f}")
    print(f"{'MEAN':16s} {'':5s} | {np.mean(agg['density']):9.4f} {np.mean(agg['composition']):9.4f} "
          f"{np.mean(agg['hybrid']):9.4f}")

    print("\n" + "=" * 100)
    print("F. THE UNIFICATION TEST — cliff sweep with source error (gDNA300/RNA150, be=1.0, R=true, diffuse).")
    print("   Claim: hybrid → DENSITY at logR≈0 (small exon cliff = the safe clean win) and → COMPOSITION at")
    print("   big logR (intron cliff = the shift). So it is safe where pure composition over-commits:")
    print(f"{'e_dst/e_src':>12s} {'logR':>7s} | {'density':>9s} {'composit':>9s} {'HYBRID':>9s}   note")
    for e_dst in (1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 30.0, 100.0):
        o, _, _ = run(g, r, e_src=1.0, e_dst=e_dst, belief_err=1.0, r_use="true")
        note = "hyb≈dens" if abs(o["hybrid"] - o["density"]) < abs(o["hybrid"] - o["composition"]) else "hyb≈comp"
        print(f"{e_dst:12.2f} {math.log(e_dst):7.2f} | {o['density']:9.4f} {o['composition']:9.4f} "
              f"{o['hybrid']:9.4f}   {note}")


if __name__ == "__main__":
    main()
