"""P1d — is the GRAFT the missing mirror of M10? A code-level + measured audit.

Answers, with numbers:
  A. STRUCTURAL — does the graft assert phi == 1 identically (no share, no share variance)?
     Exact replay checks against the shipped arrays, both the sum and the peel's twin.
  B. THE PREMISE — Var(log phi) with the ORACLE capture step, decomposed by
       (i) the junction's spliced COUNT, with the EXACT trigamma psi'(n) Poisson baseline (not 1/n)
       (ii) the boundary's SPLICED SHARE w_mu -- the test that says WHERE the term belongs:
            spliced-arm-only (=> Var(log phi) ~ w_mu^2 * const, and RIGEL_XVAR's placement is right)
            vs whole-sum (=> Var(log phi) flat in w_mu, and a SHARE on tp/tn is right)
       (iii) the OTHER junction's flux share -- the ledger's proposed prior-free observable
  C. WHAT THE GRAFT ACTUALLY CHARGES -- the shipped v_mu (1/SP_mass, NOT 1/n_spl) vs the measured premise.
  D. XVAR REACH -- what fraction of the delivered RNA precision the probe can even touch.
  E. THE MODEL-SIDE candidate: phi_hat = (tp+tn)/nu_dst, where nu_dst is the destination exon's fused RNA
     level ALREADY computed by _peel_share at every node on every face. Includes its precision provenance
     (own / far / M11) -- the circularity audit.

    OMP_NUM_THREADS=1 python scratchpad/x2_graft_mirror.py [cond ...]
"""

from __future__ import annotations

import dataclasses
import importlib
import sys
from pathlib import Path

import numpy as np
from scipy.special import polygamma

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from selfsolve_diag import _scan_and_truth  # noqa: E402

from rigel.calibration.bp_solver import REGION  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
_EPS = 1e-9
CONDS = sys.argv[1:] or [
    "gdna_gdna300_ss_0.99_nrna_present_capture_off",
    "gdna_gdna300_ss_0.50_nrna_present_capture_off",
    "gdna_gdna300_ss_0.99_nrna_present_capture_on",
    "gdna_gdna300_ss_0.50_nrna_none_capture_on",
]

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)


def qbin(x, k=5):
    """Quantile bins -- data-defined edges, no magic thresholds."""
    q = np.quantile(x, np.linspace(0.0, 1.0, k + 1))
    q[0] -= 1e-12
    q[-1] += 1e-12
    return np.clip(np.digitize(x, q[1:-1]), 0, k - 1), q


def summ(v):
    return f"n={v.size:>5d} var={np.var(v):>7.3f} sd={np.std(v):>6.3f} med={np.median(v):>7.3f}"


for cond in CONDS:
    print("=" * 118)
    print(f"CONDITION  {cond}")
    print("=" * 118)
    inp = _scan_and_truth(
        SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache"
    )
    dbg: dict = {}
    calmod.calibrate(
        inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
        np.asarray(inp["rna_fl_pmf"]),
        dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg,
    )
    chain, cap, geo = dbg["chain"], dbg["capture"], dbg["geometry"]
    us = cap["_uni_static"]
    kind = np.asarray(chain.kind)
    idx = np.asarray(chain.ref_idx, np.int64)
    isr = kind == REGION
    n = kind.shape[0]

    def pool(k):
        a = np.asarray(inp["region_pools"][k], float)
        b = np.asarray(inp["boundary_pools"][k], float)
        return np.where(isr, a[np.clip(idx, 0, a.shape[0] - 1)], b[np.clip(idx, 0, b.shape[0] - 1)])

    G = pool("gdna_pos") + pool("gdna_neg")
    Ru = pool("mat_uns_pos") + pool("nas_uns_pos") + pool("mat_uns_neg") + pool("nas_uns_neg")
    Rs = pool("mat_spl") + pool("nas_spl")
    M, E_g, E_r = us["M"], us["E_g"], us["E_r"]
    li, ri = us["left"], us["right"]
    is_bnd, is_exon, fp, fn = us["is_bnd"], us["is_exon"], us["fp"], us["fn"]
    SP = (us["SP_l"], us["SP_r"])
    SN = (us["SN_l"], us["SN_r"])
    NSP = (us["spl_n_pos_l"], us["spl_n_pos_r"])
    NSN = (us["spl_n_neg_l"], us["spl_n_neg_r"])
    ESP = (np.asarray(geo.eff_spl_left, float), np.asarray(geo.eff_spl_right, float))
    spl_p_f = tuple(np.where(SP[k] > _EPS, SP[k] / np.maximum(ESP[k], _EPS), 0.0) for k in (0, 1))
    spl_n_f = tuple(np.where(SN[k] > _EPS, SN[k] / np.maximum(ESP[k], _EPS), 0.0) for k in (0, 1))
    fwd = (us["fwd_g"], us["fwd_p"], us["fwd_n"])
    bwd = (us["bwd_g"], us["bwd_p"], us["bwd_n"])
    pin = cap["_pin"][-2:]      # the LAST rho-iteration: [df=0 (left msg), df=1 (right msg)]
    lvl = cap["_lvl"][-4:]      # same iteration: [df0-pos, df0-neg, df1-pos, df1-neg]
    uni = cap["_uni"][-1]

    # ══ A. STRUCTURAL — does the graft assert phi == 1? ════════════════════════════════════════════
    print("\n-- A. STRUCTURAL replay of both edges (exact, against the shipped arrays) --")
    for pe in pin:
        df = pe["df"]
        sf = 1 - df
        src = pe["src"]
        rel = fwd if df == 0 else bwd
        rg_s, rp_s, rn_s = rel[0][src], rel[1][src], rel[2][src]
        r = np.where(rg_s > _EPS, pe["tg"] / np.maximum(rg_s, _EPS), np.nan)
        gr = pe["graft"]
        gp = np.where(gr, spl_p_f[sf][src], 0.0)
        gn = np.where(gr, spl_n_f[sf][src], 0.0)
        # A1 the GRAFT: tp == (rho_nu_src + rho_mu_src)*r  -- coefficient 1 on BOTH arms, no share
        mgr = gr & fp & np.isfinite(r) & (pe["tp"] > _EPS)
        exp_p = (rp_s + gp) * r
        e1 = np.abs(pe["tp"][mgr] - exp_p[mgr]) / np.maximum(np.abs(exp_p[mgr]), _EPS)
        mgn = gr & fn & np.isfinite(r) & (pe["tn"] > _EPS)
        exp_n = (rn_s + gn) * r
        e1n = np.abs(pe["tn"][mgn] - exp_n[mgn]) / np.maximum(np.abs(exp_n[mgn]), _EPS)
        # A2 the PEEL: tp == rho_R_src*r*w   -- the share IS applied
        pl = is_bnd & is_exon[src] & pe["valid"]
        w_p = lvl[2 * df + 0]["w"]
        mp = pl & fp & np.isfinite(r) & (rp_s > _EPS)
        exp_pl = rp_s * r * w_p
        e2 = np.abs(pe["tp"][mp] - exp_pl[mp]) / np.maximum(np.abs(exp_pl[mp]), _EPS)
        print(f"   df={df}  GRAFT edges={int(mgr.sum()):>5d}  max|tp/((rho_nu+rho_mu)r) - 1| = "
              f"{(e1.max() if e1.size else np.nan):.3e}   (neg arm {int(mgn.sum()):>5d}: "
              f"{(e1n.max() if e1n.size else np.nan):.3e})")
        print(f"   df={df}  PEEL  edges={int(mp.sum()):>5d}  max|tp/(rho_R*r*w)      - 1| = "
              f"{(e2.max() if e2.size else np.nan):.3e}   -> the share IS applied on the peel")
        print(f"   df={df}  peel share w: med={np.median(w_p[mp]):.3f}  "
              f"Var(log w) med={np.median(lvl[2 * df]['w'][mp] * 0 + 1):.0f}"
              f" (placeholder)  |  graft share applied: IDENTICALLY 1.0")

    # ══ B. THE PREMISE, with the ORACLE capture step ══════════════════════════════════════════════
    rho_g = np.where(E_g > _EPS, G / np.maximum(E_g, _EPS), np.nan)
    rho_R_ex = np.where(E_r > _EPS, (Ru + Rs) / np.maximum(E_r, _EPS), np.nan)
    rho_nu_b = np.where(E_r > _EPS, Ru / np.maximum(E_r, _EPS), np.nan)

    rows = []  # (phi, w_mu, n_spl, s_other, spl_mass, phi_hat, po,pf,pm, nu)
    for df, nbr in ((0, li), (1, ri)):
        face = 1 - df  # the SOURCE boundary's face that looks at the exon
        oth = ri if df == 0 else li
        pe = pin[df]
        nu_p, nu_n = lvl[2 * df + 0]["nu"], lvl[2 * df + 1]["nu"]
        pv = lvl[2 * df + 0]
        for i in np.flatnonzero(is_exon):
            b = nbr[i]
            if b < 0 or not is_bnd[b]:
                continue
            mu = (SP[face][b] + SN[face][b]) / max(ESP[face][b], _EPS)
            nspl = NSP[face][b] + NSN[face][b]
            smass = SP[face][b] + SN[face][b]
            if not (mu > _EPS) or not np.isfinite(rho_R_ex[i]) or rho_R_ex[i] <= _EPS:
                continue
            if not (rho_g[b] > _EPS and rho_g[i] > _EPS) or not (rho_nu_b[b] >= 0):
                continue
            step = rho_g[b] / rho_g[i]
            phi = (rho_nu_b[b] + mu) / (rho_R_ex[i] * step)
            w_mu = mu / max(rho_nu_b[b] + mu, _EPS)
            ob = oth[i]
            mu_o, n_o = np.nan, 0.0
            if ob >= 0 and is_bnd[ob]:
                fo = 1 - (1 - df)  # the OTHER boundary's face looking back at the exon
                mu_o = (SP[fo][ob] + SN[fo][ob]) / max(ESP[fo][ob], _EPS)
                n_o = NSP[fo][ob] + NSN[fo][ob]
            nud = nu_p[i] + nu_n[i]
            phih = (pe["tp"][i] + pe["tn"][i]) / max(nud, _EPS) if nud > _EPS else np.nan
            # the DENSITY share the graft's measured spliced arm supplies, in the DELIVERED claim,
            # and the PRECISION share the same arm supplies -- the mode/precision mismatch test.
            d_spl = (spl_p_f[face][b] + spl_n_f[face][b])
            d_tot = pe["tp"][i] + pe["tn"][i]
            rel = fwd if df == 0 else bwd
            rr = pe["tg"][i] / max(rel[0][b], _EPS) if rel[0][b] > _EPS else np.nan
            w_dens = (d_spl * rr) / max(d_tot, _EPS) if np.isfinite(rr) and d_tot > _EPS else np.nan
            pprec = pe["tpp"][i] + pe["tpn"][i]
            w_prec = pe["spl_prec"][i] / max(pprec, _EPS) if pprec > _EPS else np.nan
            rows.append((phi, w_mu, nspl, mu_o, smass, phih,
                         pv["po"][i], pv["pf"][i], pv["pm"][i], nud, mu, w_dens, w_prec,
                         rel[1][b] + rel[2][b], 1.0 / max(pe["tpg"][i], _EPS), n_o))
    A = np.asarray([r_ for r_ in rows if np.isfinite(r_[0]) and r_[0] > 0], float)
    if A.shape[0] < 20:
        print("   (too few edges)")
        continue
    phi, w_mu, nspl, mu_o, smass, phih = A[:, 0], A[:, 1], A[:, 2], A[:, 3], A[:, 4], A[:, 5]
    po, pf, pm = A[:, 6], A[:, 7], A[:, 8]
    mu_this, w_dens, w_prec = A[:, 10], A[:, 11], A[:, 12]
    lp = np.log(phi)
    print(f"\n-- B. Var(log phi), ORACLE capture step --  {summ(lp)}")

    print("\n   B(i) by junction SPLICED COUNT  [EXACT trigamma baseline]")
    print(f"   {'n_spl bin':<22}{'edges':>7}{'med n':>8}{'Var(log phi)':>14}"
          f"{'psi1(n)':>10}{'1/n':>9}{'excess':>9}")
    bi, q = qbin(nspl)
    for k in range(5):
        m = bi == k
        if m.sum() < 5:
            continue
        v = np.var(lp[m])
        p1 = float(np.mean(polygamma(1, np.maximum(nspl[m], 1.0))))
        print(f"   [{q[k]:>8.1f},{q[k+1]:>8.1f}]{'':<3}{int(m.sum()):>7d}"
              f"{np.median(nspl[m]):>8.1f}{v:>14.3f}{p1:>10.4f}"
              f"{float(np.mean(1.0/np.maximum(nspl[m],1.0))):>9.4f}{v - p1:>9.3f}")

    print("\n   B(ii) by SPLICED SHARE w_mu  -- WHERE does the term belong?")
    print(f"   {'w_mu bin':<22}{'edges':>7}{'med w_mu':>10}{'Var(log phi)':>14}"
          f"{'/w_mu^2':>12}{'/(1-w_mu)^2':>13}")
    bi, q = qbin(w_mu)
    for k in range(5):
        m = bi == k
        if m.sum() < 5:
            continue
        v = np.var(lp[m])
        wm = float(np.median(w_mu[m]))
        print(f"   [{q[k]:>8.3f},{q[k+1]:>8.3f}]{'':<3}{int(m.sum()):>7d}{wm:>10.3f}"
              f"{v:>14.3f}{v/max(wm*wm,1e-6):>12.2f}{v/max((1-wm)**2,1e-6):>13.2f}")

    ok_o = np.isfinite(mu_o) & (mu_o > _EPS) & (mu_this > _EPS)
    if ok_o.sum() > 20:
        s_share = mu_this[ok_o] / (mu_this[ok_o] + mu_o[ok_o])
        print("\n   B(iii) by the OTHER junction's flux share  s = mu_this/(mu_this+mu_other)"
              "   [the ledger 5.2 candidate observable]")
        print(f"   {'s bin':<22}{'edges':>7}{'med s':>8}{'Var(log phi)':>14}{'med log phi':>13}")
        bi, q = qbin(s_share)
        lpo = lp[ok_o]
        for k in range(5):
            m = bi == k
            if m.sum() < 5:
                continue
            print(f"   [{q[k]:>8.3f},{q[k+1]:>8.3f}]{'':<3}{int(m.sum()):>7d}"
                  f"{np.median(s_share[m]):>8.3f}{np.var(lpo[m]):>14.3f}{np.median(lpo[m]):>13.3f}")

    print("\n   B(iv) ** THE MODE/PRECISION MISMATCH ** -- the grafted spliced arm's share of the")
    print("        delivered RNA DENSITY vs its share of the delivered RNA PRECISION")
    okw = np.isfinite(w_dens) & np.isfinite(w_prec)
    print(f"   {'':<22}{'p25':>10}{'median':>10}{'p75':>10}")
    for nm, v in (("density share w_dens", w_dens[okw]), ("precision share w_prec", w_prec[okw])):
        print(f"   {nm:<22}{np.percentile(v,25):>10.3f}{np.percentile(v,50):>10.3f}"
              f"{np.percentile(v,75):>10.3f}")
    print(f"   ratio  w_dens/w_prec  median = "
          f"{np.median(w_dens[okw]/np.maximum(w_prec[okw],1e-9)):.2f}x"
          f"   -> a variance placed on the SPLICED ARM ONLY reaches "
          f"{np.median(w_prec[okw])*100:.0f}% of the delivered confidence,")
    print(f"      while the premise error contaminates {np.median(w_dens[okw])*100:.0f}% "
          f"of the delivered density.")

    # ══ C. WHAT THE GRAFT CHARGES vs the premise ══════════════════════════════════════════════════
    print("\n-- C. what the graft is actually CHARGED (its own count term) --")
    v_shipped = 1.0 / np.maximum(smass, _EPS)     # the code: precision = SP mass  => v = 1/mass
    v_ledger = 1.0 / np.maximum(nspl, _EPS)       # the ledger's stated 1/n_spl
    v_trig = polygamma(1, np.maximum(nspl, 1.0))  # the EXACT trigamma of the count
    print(f"   median spliced MASS  {np.median(smass):>9.1f}   median spliced COUNT {np.median(nspl):>9.1f}"
          f"   mass/count = {np.median(smass)/max(np.median(nspl),1e-9):.3f}")
    print(f"   v charged (1/mass, SHIPPED)  med={np.median(v_shipped):.5f}  mean={v_shipped.mean():.5f}")
    print(f"   v if 1/n_spl (ledger text)   med={np.median(v_ledger):.5f}  mean={v_ledger.mean():.5f}")
    print(f"   v if psi'(n_spl) (exact)     med={np.median(v_trig):.5f}  mean={v_trig.mean():.5f}")
    print(f"   MEASURED premise Var(log phi) = {np.var(lp):.4f}"
          f"   =>  over-confidence  {np.var(lp)/max(np.median(v_shipped),1e-12):>9.1f}x (vs shipped)")
    print(f"   {'n_spl bin':<22}{'edges':>7}{'Var(log phi)':>14}{'med v charged':>15}{'over-conf':>11}")
    bi, q = qbin(nspl)
    for k in range(5):
        m = bi == k
        if m.sum() < 5:
            continue
        vv, vc = np.var(lp[m]), float(np.median(v_shipped[m]))
        print(f"   [{q[k]:>8.1f},{q[k+1]:>8.1f}]{'':<3}{int(m.sum()):>7d}{vv:>14.3f}"
              f"{vc:>15.5f}{vv/max(vc,1e-12):>10.0f}x")

    # ══ G. the THIRD stream: the graft moves the lambda MODE at a precision that never saw it ══════
    rnu_src = A[:, 13]
    okg = rnu_src > _EPS
    dlam = np.log((rnu_src[okg] + mu_this[okg]) / np.maximum(rnu_src[okg], _EPS))
    ct = np.asarray(uni["c_tau"], float)
    exg = np.zeros(n, bool)
    for pe in pin:
        exg |= pe["graft"]
    live = exg & (ct > _EPS)
    print("\n-- G. the COMPOSITION (lambda) stream: the graft sets its MODE but never its PRECISION --")
    print(f"   the graft shifts the message's lambda mode by  -log(1 + rho_mu/rho_nu):"
          f"  med={np.median(dlam):.3f} nats  p75={np.percentile(dlam,75):.3f}  n={dlam.size}")
    print(f"   the premise error enters lambda at FULL weight (d lambda = -log phi)"
          f"  =>  Var = {np.var(lp):.3f}")
    print(f"   the STATED Var(lambda) = 1/c_tau at graft-fed exons:"
          f"  med={np.median(1.0/np.maximum(ct[live],_EPS)):.4f}  n={int(live.sum())}"
          f"   -> over-confident {np.var(lp)*np.median(ct[live]):.0f}x")

    # ══ D. XVAR REACH ═════════════════════════════════════════════════════════════════════════════
    print("\n-- D. RIGEL_XVAR reach: what fraction of the delivered RNA precision it can touch --")
    tot_gr = 0
    for pe in pin:
        gr = pe["graft"]
        sp = pe["spl_prec"]
        m = gr & (pe["tpp"] + pe["tpn"] > 0)
        frac_mode = sp[m] / np.maximum(pe["tpp"][m] + pe["tpn"][m], _EPS)
        mm = gr & (pe["tmp"] + pe["tmn"] > 0) if "tmp" in pe else None
        print(f"   df={pe['df']}  graft edges={int(gr.sum()):>5d}   "
              f"_spc/(tpp+tpn): med={np.median(frac_mode):.3f} "
              f"p25={np.percentile(frac_mode,25):.3f} p75={np.percentile(frac_mode,75):.3f}")
        del mm
        tot_gr += int(gr.sum())
    # the relay fires the SAME graft edge set, and XVAR never reaches it
    relay_graft = int((is_exon & (li >= 0) & is_bnd[np.clip(li, 0, n - 1)]).sum()
                      + (is_exon & (ri >= 0) & is_bnd[np.clip(ri, 0, n - 1)]).sum())
    print(f"   graft edge FIRINGS per sweep:  combine (XVAR-reachable) = {tot_gr//len(pin)*len(pin)}"
          f" x {len([1])}   |   relay (NOT reachable) = {relay_graft}")
    print(f"   -> XVAR reaches {tot_gr} of {tot_gr + relay_graft} graft firings "
          f"({100.0*tot_gr/max(tot_gr+relay_graft,1):.0f}%), and the relay's firing COMPOUNDS along the chain")

    # ══ E. the MODEL-SIDE candidate + its circularity audit ═══════════════════════════════════════
    okh = np.isfinite(phih) & (phih > 0)
    if okh.sum() > 20:
        lh = np.log(phih[okh])
        c = np.corrcoef(lh, lp[okh])[0, 1]
        print("\n-- E. model-side phi_hat = (tp+tn)/nu_dst  (nu_dst ALREADY computed by _peel_share) --")
        print(f"   {summ(lh)}   corr(log phi_hat, log phi_oracle) = {c:>6.3f}   n={int(okh.sum())}")
        pt = po + pf + pm
        okp = pt > 0
        print(f"   nu_dst precision provenance at graft destinations (share of total precision):")
        print(f"      own  po/pt  med={np.median((po/np.maximum(pt,_EPS))[okp]):.3f}"
              f"   far  pf/pt  med={np.median((pf/np.maximum(pt,_EPS))[okp]):.3f}"
              f"   M11  pm/pt  med={np.median((pm/np.maximum(pt,_EPS))[okp]):.3f}")
        print(f"      fraction of edges where M11 is >90% of nu's precision: "
              f"{float(np.mean((pm/np.maximum(pt,_EPS))[okp] > 0.9)):.3f}  (circularity risk)")

    # ══ F. the two-junction MoM, with the EXACT trigamma ══════════════════════════════════════════
    gaps, npair = [], []
    for i in np.flatnonzero(is_exon):
        lb, rb = li[i], ri[i]
        if lb < 0 or rb < 0 or not is_bnd[lb] or not is_bnd[rb]:
            continue
        ml = (SP[1][lb] + SN[1][lb]) / max(ESP[1][lb], _EPS)
        mr = (SP[0][rb] + SN[0][rb]) / max(ESP[0][rb], _EPS)
        nl = NSP[1][lb] + NSN[1][lb]
        nr = NSP[0][rb] + NSN[0][rb]
        if ml > _EPS and mr > _EPS and nl > 0 and nr > 0:
            gaps.append(np.log(ml / mr))
            npair.append((nl, nr))
    gaps = np.asarray(gaps)
    if gaps.size > 20:
        npair = np.asarray(npair, float)
        pois_trig = float(np.mean(polygamma(1, np.maximum(npair[:, 0], 1.0))
                                  + polygamma(1, np.maximum(npair[:, 1], 1.0))))
        pois_inv = float(np.mean(1.0 / npair[:, 0] + 1.0 / npair[:, 1]))
        vg = float(np.var(gaps))
        print("\n-- F. two-junction MoM estimator (ledger 5.1), Poisson part by EXACT trigamma --")
        print(f"   pairs={gaps.size:>5d}  Var(gap)={vg:>7.3f}  "
              f"E[psi1(n_l)+psi1(n_r)]={pois_trig:>7.4f}  E[1/n_l+1/n_r]={pois_inv:>7.4f}")
        print(f"   omega_trigamma = {max(0.0,(vg-pois_trig))/2:>7.4f}   "
              f"omega_1overn = {max(0.0,(vg-pois_inv))/2:>7.4f}   "
              f"TRUE Var(log phi) = {np.var(lp):>7.4f}")

    # ══ H. the PER-EDGE estimator implied by B(iii): a single-observation second moment ════════════
    # d = log(rho_mu_this / rho_mu_other) is ONE draw from a distribution of variance 2*omega, so
    #     v_P1d = max(0, d^2 - psi'(n_this) - psi'(n_other)) / 2
    # is the method-of-moments second moment of a single observation -- the same construction as M8's
    # (log r)^2 and M7's b-hat^2. No tuned constant, no pooling, prior-free, and BOTH inputs are in scope.
    n_o = A[:, 15]
    okd = ok_o & (n_o > 0) & (nspl > 0)
    if okd.sum() > 20:
        d = np.log(mu_this[okd] / np.maximum(mu_o[okd], _EPS))
        pois = polygamma(1, np.maximum(nspl[okd], 1.0)) + polygamma(1, np.maximum(n_o[okd], 1.0))
        v_p1d = np.maximum(d * d - pois, 0.0) / 2.0
        e2 = lp[okd] ** 2
        print("\n-- H. the PER-EDGE estimator  v_P1d = max(0, d^2 - psi'(n_this) - psi'(n_other))/2 --")
        print(f"   edges={int(okd.sum()):>5d}   E[(log phi)^2]={e2.mean():>7.3f}   "
              f"E[v_P1d]={v_p1d.mean():>7.3f}   z2 = E[dsq]/E[v] = {e2.mean()/max(v_p1d.mean(),1e-9):>6.2f}")
        print(f"   {'v_P1d quintile':<22}{'edges':>7}{'med v_P1d':>11}{'E[(log phi)^2]':>16}{'z2':>8}")
        bi, q = qbin(v_p1d)
        for k in range(5):
            m = bi == k
            if m.sum() < 5:
                continue
            print(f"   [{q[k]:>8.3f},{q[k+1]:>8.3f}]{'':<3}{int(m.sum()):>7d}"
                  f"{np.median(v_p1d[m]):>11.3f}{e2[m].mean():>16.3f}"
                  f"{e2[m].mean()/max(v_p1d[m].mean(),1e-9):>8.2f}")
        print(f"   pooled-omega comparison: a single omega={np.var(lp):.3f} for every edge would give"
              f"  z2 by quintile = "
              + " ".join(f"{e2[qbin(v_p1d)[0]==k].mean()/np.var(lp):.2f}" for k in range(5)))
