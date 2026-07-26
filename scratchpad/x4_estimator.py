"""P1d re-audit — the two-junction estimator under the FULL Poisson (trigamma) subtraction.

Starting point: a COPY of scratchpad/x1_graft.py (same oracle frame, same edge enumeration), extended with

  (A) the FULL Poisson accounting on BOTH sides of the comparison.
      log phi = log[rho_nu(b) + rho_mu(b)] - log rho_R(exon) - log rho_g(b) + log rho_g(exon)
      so Var(log phi) carries FIVE independent count terms, of which the ledger's estimator subtracts ONE
      (the junction's spliced count) and only in its 1/n approximation:
          w_nu^2 * psi'(n_nu^b)   the boundary's unspliced-RNA count
        + w_mu^2 * psi'(n_spl)    the junction's spliced count          <- the only one subtracted today
        + psi'(n_R^exon)          the exon's own RNA count
        + psi'(n_g^b) + psi'(n_g^exon)   the ORACLE capture step's two gDNA counts
      Everything uses the EXACT trigamma psi'(n) (project precedent: M11 `residual_level`), with the 1/n
      form printed alongside for contrast.
      CAPTURE-OFF: the true frame step is 1 EXACTLY, so phi is formed with step==1 and the two gDNA counts
      VANISH rather than being subtracted. That is the clean target, it needs no gDNA at all (so the
      gdna_none and gdna1 conditions become measurable), and TABLE A-OFF cross-checks it against the
      oracle-gDNA route minus its two psi'(n_g) terms.

  (B) the divide-by-2's independence assumption, measured DIRECTLY: phi_left and phi_right against the SAME
      exon's oracle RNA density, their Pearson correlation, and the implied divisor 2*(1-rho). Both raw and
      after removing the SHARED Poisson terms (the exon's own RNA count -- and, on capture-ON, its gDNA
      count -- are common to both sides and manufacture correlation that is not premise correlation).

  (C) SELECTION. The estimator is computable only on exons with TWO live flanking junctions. Var(log phi) is
      reported on: all graft edges / exon has 2 live junctions / 1 / the edge carries no spliced flux.

  (G) candidate per-junction observables, all prior-free, and whether they are computable at a SINGLE edge.

    OMP_NUM_THREADS=1 python scratchpad/x4_estimator.py [cond ...]
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
from rigel.calibration.node_geometry import _node_region_type  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
_EPS = 1e-9

CAP_OFF = [
    "gdna_gdna300_ss_0.99_nrna_present_capture_off",
    "gdna_gdna300_ss_0.50_nrna_present_capture_off",
    "gdna_gdna300_ss_0.99_nrna_none_capture_off",
    "gdna_gdna300_ss_0.50_nrna_none_capture_off",
    "gdna_gdna100_ss_0.99_nrna_present_capture_off",
    "gdna_gdna100_ss_0.50_nrna_present_capture_off",
    "gdna_gdna100_ss_0.99_nrna_none_capture_off",
    "gdna_gdna100_ss_0.50_nrna_none_capture_off",
    "gdna_gdna5_ss_0.50_nrna_present_capture_off",
    "gdna_gdna1_ss_0.50_nrna_present_capture_off",
    "gdna_none_ss_0.99_nrna_present_capture_off",
    "gdna_none_ss_0.50_nrna_present_capture_off",
    "gdna_none_ss_0.99_nrna_none_capture_off",
    "gdna_none_ss_0.50_nrna_none_capture_off",
]
CONTRAST = [
    "gdna_gdna300_ss_0.99_nrna_present_capture_on",
    "gdna_gdna300_ss_0.50_nrna_none_capture_on",
    "gdna_gdna100_ss_0.50_nrna_present_capture_verystrong",
]
CONDS = sys.argv[1:] or (CAP_OFF + CONTRAST)

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)


def tri(n):
    """EXACT trigamma log-count variance, floored at one fragment (M11's `k >= 1` limit: psi'(1)=1.6449)."""
    return polygamma(1, np.maximum(np.asarray(n, float), 1.0))


def inv(n):
    """the 1/n approximation the shipped model and the ledger's estimator use."""
    return 1.0 / np.maximum(np.asarray(n, float), 1.0)


def cc(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 5 or np.std(a[ok]) < 1e-12 or np.std(b[ok]) < 1e-12:
        return np.nan
    return float(np.corrcoef(a[ok], b[ok])[0, 1])


def reg_resid(y, x):
    """premise variance left AFTER a 1-parameter linear regression of y on x (and the R^2)."""
    y, x = np.asarray(y, float), np.asarray(x, float)
    ok = np.isfinite(y) & np.isfinite(x)
    if ok.sum() < 10 or np.std(x[ok]) < 1e-12:
        return np.nan, np.nan
    b = np.polyfit(x[ok], y[ok], 1)
    r = y[ok] - np.polyval(b, x[ok])
    v0 = float(np.var(y[ok]))
    return float(np.var(r)), (1.0 - float(np.var(r)) / v0 if v0 > 0 else np.nan)


ROWS: dict[str, list] = {k: [] for k in "SAOBCDEFGHLI"}

for cond in CONDS:
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

    def pool(k, isr=isr, idx=idx, inp=inp):
        a = np.asarray(inp["region_pools"][k], float)
        b = np.asarray(inp["boundary_pools"][k], float)
        return np.where(isr, a[np.clip(idx, 0, a.shape[0] - 1)], b[np.clip(idx, 0, b.shape[0] - 1)])

    G = pool("gdna_pos") + pool("gdna_neg")
    Ru = pool("mat_uns_pos") + pool("nas_uns_pos") + pool("mat_uns_neg") + pool("nas_uns_neg")
    Rs = pool("mat_spl") + pool("nas_spl")
    M, E_g, E_r = us["M"], us["E_g"], us["E_r"]
    li, ri = us["left"], us["right"]
    is_bnd, is_exon = us["is_bnd"], us["is_exon"]
    SPf, SNf = (us["SP_l"], us["SP_r"]), (us["SN_l"], us["SN_r"])
    SNPf = (us["spl_n_pos_l"], us["spl_n_pos_r"])   # spliced fragment COUNTS per face
    SNNf = (us["spl_n_neg_l"], us["spl_n_neg_r"])
    ESP = (np.asarray(geo.eff_spl_left, float), np.asarray(geo.eff_spl_right, float))
    rt, _ = _node_region_type(chain, ra)
    capoff = cond.endswith("capture_off")

    rho_g = np.where(E_g > _EPS, G / np.maximum(E_g, _EPS), np.nan)
    rho_R_ex = np.where(E_r > _EPS, (Ru + Rs) / np.maximum(E_r, _EPS), np.nan)
    rho_nu_b = np.where(E_r > _EPS, Ru / np.maximum(E_r, _EPS), np.nan)
    rho_M = np.where(E_r > _EPS, M / np.maximum(E_r, _EPS), np.nan)  # OBSERVABLE total-mass density

    rec: list[dict] = []
    for face, nbr in ((1, li), (0, ri)):
        for i in np.flatnonzero(is_exon):
            b = nbr[i]
            if b < 0 or not is_bnd[b]:
                continue
            n_spl = float(SNPf[face][b] + SNNf[face][b])
            mu = (SPf[face][b] + SNf[face][b]) / max(ESP[face][b], _EPS)
            if not np.isfinite(rho_R_ex[i]) or rho_R_ex[i] <= _EPS:
                continue
            if not np.isfinite(rho_nu_b[b]):
                continue
            claim = rho_nu_b[b] + mu
            if not (claim > _EPS):
                continue
            step_ok = bool(rho_g[b] > _EPS and rho_g[i] > _EPS)
            step = (rho_g[b] / rho_g[i]) if step_ok else np.nan
            rec.append(dict(
                i=int(i), b=int(b), face=face, mu=float(mu), n_spl=n_spl, claim=float(claim),
                step=float(step) if step_ok else np.nan,
                lphi_g=float(np.log(claim / (rho_R_ex[i] * step))) if step_ok else np.nan,
                lphi_1=float(np.log(claim / rho_R_ex[i])),
                w_nu=float(rho_nu_b[b] / claim), w_mu=float(mu / claim),
                n_nu=float(Ru[b]), n_gb=float(G[b]), n_gi=float(G[i]), n_Rex=float(Ru[i] + Rs[i]),
                # OBSERVABLE single-edge candidates (no oracle, no second junction)
                o_exon=float(mu / rho_M[i]) if rho_M[i] > _EPS else np.nan,
                o_bnd=float(mu / rho_M[b]) if rho_M[b] > _EPS else np.nan,
                # the graft's OWN share weight w_mu (M2) — here with the ORACLE boundary RNA in place of
                # the solver's, to test the CONCEPT. Shares nothing with the exon-side denominator of phi.
                o_wmu=float(mu / (rho_nu_b[b] + mu)),
                live=bool(n_spl > 0.0 and mu > _EPS),
            ))
    if not rec:
        continue
    R = {k: np.asarray([r[k] for r in rec]) for k in rec[0]}
    KEY = "lphi_1" if capoff else "lphi_g"   # capture-OFF: step==1 is EXACT, no gDNA counts needed

    def claim_pois(f, k=slice(None)):
        """the CLAIM's own Poisson variance: boundary unspliced-RNA count + junction spliced count."""
        t_mu = np.where(R["n_spl"][k] > 0, R["w_mu"][k] ** 2 * f(np.maximum(R["n_spl"][k], 1.0)), 0.0)
        return R["w_nu"][k] ** 2 * f(R["n_nu"][k]) + t_mu

    def shared_pois(f, k=slice(None)):
        """the terms COMMON to an exon's two edges: its own RNA count (+ its gDNA count on capture-ON)."""
        s = f(R["n_Rex"][k])
        return s if capoff else s + f(R["n_gi"][k])

    def side_pois(f, k=slice(None)):
        """everything NOT shared: the claim, plus (capture-ON only) this boundary's gDNA count."""
        return claim_pois(f, k) + (0.0 if capoff else f(R["n_gb"][k]))

    def full_pois(f, k=slice(None)):
        return side_pois(f, k) + shared_pois(f, k)

    def prem(mask, f=tri, key=None):
        key = key or KEY
        x = R[key][mask]
        if x.size < 5:
            return x.size, np.nan, np.nan, np.nan
        v = float(np.var(x))
        p = float(np.mean(full_pois(f, mask))) if key == KEY else float(
            np.mean(claim_pois(f, mask) + shared_pois(f, mask))
        )
        return x.size, v, p, v - p

    m_live = R["live"] & np.isfinite(R[KEY])
    m_all = np.isfinite(R[KEY])

    sok = np.isfinite(R["step"])
    ROWS["S"].append((cond, int(m_all.sum()),
                      float(np.median(R["step"][sok])) if sok.sum() else np.nan,
                      float(np.std(np.log(R["step"][sok]))) if sok.sum() else np.nan,
                      float(np.mean(tri(R["n_gb"]) + tri(R["n_gi"]))),
                      float(np.median(R["n_gb"])), float(np.median(R["n_Rex"])),
                      float(np.median(R["n_spl"][m_live])) if m_live.sum() else np.nan))

    nA, vA, pA, qA = prem(m_live, tri)
    _, _, pA1, qA1 = prem(m_live, inv)
    ROWS["A"].append((cond, nA, vA, pA1, qA1, pA, qA, float(np.median(R[KEY][m_live]))))

    if capoff:  # the cross-check: oracle-gDNA route MINUS its two psi'(n_g), vs the exact step==1 route
        mg = m_live & np.isfinite(R["lphi_g"])
        if mg.sum() >= 5:
            vg = float(np.var(R["lphi_g"][mg]))
            pg = float(np.mean(claim_pois(tri, mg) + tri(R["n_Rex"][mg])
                               + tri(R["n_gb"][mg]) + tri(R["n_gi"][mg])))
            ROWS["O"].append((cond, int(mg.sum()), vg, pg, vg - pg, qA, (vg - pg) - qA,
                              float(np.mean(tri(R["n_gb"][mg]) + tri(R["n_gi"][mg])))))

    # ── (A/B) the estimator, on exons with TWO live flanking junctions ─────────────────────────────────
    byex: dict[int, dict] = {}
    for k in np.flatnonzero(m_live):
        byex.setdefault(int(R["i"][k]), {})[int(R["face"][k])] = k
    pairs = [(v[1], v[0]) for v in byex.values() if 0 in v and 1 in v]
    kl = np.asarray([p[0] for p in pairs], int)
    kr = np.asarray([p[1] for p in pairs], int)
    if len(pairs) >= 5:
        # E1 — the LEDGER's statistic: the two junctions' SPLICED densities only
        g1 = np.log(R["mu"][kl] / R["mu"][kr])
        v1 = float(np.var(g1))
        p1i = float(np.mean(inv(R["n_spl"][kl]) + inv(R["n_spl"][kr])))
        p1t = float(np.mean(tri(R["n_spl"][kl]) + tri(R["n_spl"][kr])))
        # E2 — the FULL graft claim differenced: log(phi_l/phi_r). The exon's own counts CANCEL.
        g2 = R[KEY][kl] - R[KEY][kr]
        v2 = float(np.var(g2))
        p2t = float(np.mean(side_pois(tri, kl) + side_pois(tri, kr)))
        p2i = float(np.mean(side_pois(inv, kl) + side_pois(inv, kr)))
        # E1-LOOSE — x1_graft.py's EXACT gap loop (no is_bnd test, no rho_R(exon) test, both faces of the
        # flanking boundaries), to reconcile against the ledger's published Var(gap).
        gl = []
        for i2 in np.flatnonzero(is_exon):
            lb, rb = li[i2], ri[i2]
            if lb < 0 or rb < 0:
                continue
            ml_ = (SPf[1][lb] + SNf[1][lb]) / max(ESP[1][lb], _EPS)
            mr_ = (SPf[0][rb] + SNf[0][rb]) / max(ESP[0][rb], _EPS)
            if ml_ > _EPS and mr_ > _EPS:
                gl.append((np.log(ml_ / mr_),
                           float(SNPf[1][lb] + SNNf[1][lb]), float(SNPf[0][rb] + SNNf[0][rb])))
        if len(gl) >= 5:
            ga = np.asarray([g[0] for g in gl])
            na, nb = np.asarray([g[1] for g in gl]), np.asarray([g[2] for g in gl])
            ROWS["L"].append((cond, len(gl), float(np.var(ga)), float(np.median(np.abs(ga))),
                              float(np.mean(inv(na) + inv(nb))), float(np.mean(tri(na) + tri(nb))),
                              max(0.0, float(np.var(ga)) - float(np.mean(inv(na) + inv(nb)))) / 2,
                              max(0.0, float(np.var(ga)) - float(np.mean(tri(na) + tri(nb)))) / 2))
        # the TARGET on the very subset the estimator is fitted on
        _, _, _, q2j = prem(m_live & np.isin(np.arange(R["i"].size), np.concatenate([kl, kr])))
        ROWS["B"].append((cond, len(pairs), v1, p1i, max(0.0, v1 - p1i) / 2, p1t,
                          max(0.0, v1 - p1t) / 2, v2, p2t, max(0.0, v2 - p2t) / 2,
                          max(0.0, v2 - p2i) / 2, q2j, qA))

        xl, xr = R[KEY][kl], R[KEY][kr]
        Vl, Vr, Cv = float(np.var(xl)), float(np.var(xr)), float(np.cov(xl, xr, ddof=0)[0, 1])
        rr = Cv / np.sqrt(max(Vl * Vr, _EPS))
        sh_p = float(np.mean(shared_pois(tri, kl)))
        Pl = Vl - float(np.mean(side_pois(tri, kl))) - sh_p
        Pr = Vr - float(np.mean(side_pois(tri, kr))) - sh_p
        Pc = Cv - sh_p
        rp = Pc / np.sqrt(max(Pl * Pr, _EPS)) if (Pl > 0 and Pr > 0) else np.nan
        dp = 2.0 * (1.0 - rp) if np.isfinite(rp) else np.nan
        ROWS["C"].append((cond, len(pairs), Vl, Vr, Cv, rr, Pl, Pr, Pc, rp,
                          2.0 * (1.0 - rr), dp,
                          (max(0.0, v2 - p2t) / dp) if (np.isfinite(dp) and dp > 0) else np.nan,
                          0.5 * (Pl + Pr)))

    # ── (C) SELECTION ─────────────────────────────────────────────────────────────────────────────────
    nlive = np.zeros(n, int)
    for k in np.flatnonzero(R["live"]):
        nlive[int(R["i"][k])] += 1
    exn = nlive[R["i"].astype(int)]
    out = [cond, int(m_all.sum()), int(m_live.sum())]
    for sel in (m_all, m_live, m_live & (exn == 2), m_live & (exn == 1), m_all & ~R["live"]):
        nn, vv, pp, qq = prem(sel)
        out += [nn, qq, (float(np.median(R[KEY][sel])) if nn else np.nan)]
    ROWS["D"].append(tuple(out))

    row = [cond]
    for lo, hi in ((0, 30), (30, 100), (100, 300), (300, 1000), (1000, 1e18)):
        nn, vv, pp, qq = prem(m_live & (R["n_spl"] >= lo) & (R["n_spl"] < hi))
        row += [nn, qq]
    ROWS["E"].append(tuple(row))

    # ── (G) candidate per-junction observables ────────────────────────────────────────────────────────
    x = R[KEY]
    if len(pairs) >= 5:
        sl = R["mu"][kl] / (R["mu"][kl] + R["mu"][kr])
        two = np.concatenate([kl, kr])
        share = np.concatenate([sl, 1.0 - sl])
        ROWS["F"].append((cond, int(m_live.sum()), len(pairs),
                          cc(x[two], np.log(np.maximum(share, 1e-12))),
                          cc(np.abs(x[two]), np.log(np.maximum(share, 1e-12))),
                          cc(x[m_live], np.log(np.maximum(R["n_spl"][m_live], 1.0))),
                          cc(np.abs(x[m_live]), np.log(np.maximum(R["n_spl"][m_live], 1.0))),
                          cc(x[m_live], np.log(np.maximum(R["o_exon"][m_live], 1e-12))),
                          cc(np.abs(x[m_live]), np.log(np.maximum(R["o_exon"][m_live], 1e-12))),
                          cc(x[m_live], np.log(np.maximum(R["o_wmu"][m_live], 1e-12))),
                          cc(np.abs(x[m_live]), np.log(np.maximum(R["o_wmu"][m_live], 1e-12)))))
        # ── (H) how much of the premise variance can ONE fitted slope on an observable remove? ────────
        pool_shared = float(np.mean(shared_pois(tri, two)))
        r_sh, R2_sh = reg_resid(x[two], np.log(np.maximum(share, 1e-12)))
        mv = np.flatnonzero(m_live)
        r_ct, R2_ct = reg_resid(x[mv], np.log(np.maximum(R["n_spl"][mv], 1.0)))
        r_ex, R2_ex = reg_resid(x[mv], np.log(np.maximum(R["o_exon"][mv], 1e-12)))
        r_bd, R2_bd = reg_resid(x[mv], np.log(np.maximum(R["o_wmu"][mv], 1e-12)))
        p_live = float(np.mean(full_pois(tri, mv)))
        # the SHARE observable, restricted to the same 2j pairs, but scored on the live-edge scale for a
        # like-for-like read against the single-edge candidates
        ROWS["H"].append((cond,
                          float(np.var(x[two])) - float(np.mean(full_pois(tri, two))),
                          r_sh - float(np.mean(full_pois(tri, two))), R2_sh,
                          float(np.var(x[mv])) - p_live,
                          r_ct - p_live, R2_ct, r_ex - p_live, R2_ex, r_bd - p_live, R2_bd))
        # ── (I) is this a BIAS with a structural fix, or genuine noise? ───────────────────────────────
        # If log phi tracks log(share) with slope ~1, the graft is not "uncertain", it is systematically
        # dropping the flux that entered through the OTHER junction. Test alternative claim forms.
        bfit = np.polyfit(np.log(np.maximum(share, 1e-12)), x[two], 1)
        rn = np.concatenate([R["claim"][kl] - R["mu"][kl], R["claim"][kr] - R["mu"][kr]])  # rho_nu
        mu_t = np.concatenate([R["mu"][kl], R["mu"][kr]])          # this junction
        mu_o = np.concatenate([R["mu"][kr], R["mu"][kl]])          # the OTHER junction
        den = np.concatenate([R["claim"][kl] / np.exp(x[kl]), R["claim"][kr] / np.exp(x[kr])])
        pshared = float(np.mean(full_pois(tri, two)))
        alts = []
        for cl in (rn + mu_t, rn + np.maximum(mu_t, mu_o), rn + mu_t + mu_o, rn + 2.0 * mu_t):
            lz = np.log(np.maximum(cl, 1e-30) / den)
            alts += [float(np.var(lz)) - pshared, float(np.median(lz))]
        ROWS["I"].append((cond, len(two), float(bfit[0]), float(bfit[1]), *alts))
        # stratify by the two-junction FLUX SHARE — mean (bias) and premise variance per stratum
        rowg = [cond]
        for lo, hi in ((0.0, 0.2), (0.2, 0.4), (0.4, 0.6), (0.6, 0.8), (0.8, 1.01)):
            s = two[(share >= lo) & (share < hi)]
            if s.size >= 5:
                v = float(np.var(x[s]))
                rowg += [s.size, float(np.mean(x[s])), v - float(np.mean(full_pois(tri, s)))]
            else:
                rowg += [s.size, np.nan, np.nan]
        ROWS["G"].append(tuple(rowg))


def show(title, hdr, rows, fmt):
    print("\n" + title)
    print(hdr)
    print("-" * max(len(hdr), 40))
    for r in rows:
        print(fmt(r))
    off = [r for r in rows if r[0].endswith("capture_off")]
    if len(off) > 1:  # MEAN over the capture-OFF conditions (where M8 is identically 0)
        m = [off[0][0].replace(off[0][0], "MEAN over %d capture-OFF" % len(off))]
        for j in range(1, len(off[0])):
            v = np.asarray([r[j] for r in off], float)
            v = v[np.isfinite(v)]
            m.append(float(v.mean()) if v.size else np.nan)
        print(fmt(tuple(m)).replace(m[0][5:34].ljust(34), m[0][:34].ljust(34)))
    print()


sh = lambda c: c[5:].replace("_nrna", "").replace("_capture", "")  # noqa: E731
f3 = lambda v: ("     nan" if not np.isfinite(v) else f"{v:>8.3f}")  # noqa: E731

show("TABLE S — the substrate: counts and the oracle capture step (median over graft edges)",
     f"{'condition':<34}{'edges':>7}{'med step':>10}{'sd(log step)':>14}{'E[2*psi_g]':>12}"
     f"{'med n_g(b)':>12}{'med n_R(ex)':>13}{'med n_spl':>11}",
     ROWS["S"],
     lambda r: f"{sh(r[0]):<34}{r[1]:>7}{r[2]:>10.3f}{r[3]:>14.3f}{r[4]:>12.3f}"
               f"{r[5]:>12.0f}{r[6]:>13.0f}{r[7]:>11.0f}")

show("TABLE A — (A) THE TARGET Var(log phi) UNDER THE FULL POISSON SUBTRACTION   [live graft edges]\n"
     "          capture-OFF uses the EXACT step==1 frame; capture-ON uses the oracle gDNA step.",
     f"{'condition':<34}{'edges':>7}{'Var raw':>10}{'E[P] 1/n':>10}{'prem 1/n':>10}"
     f"{'E[P] tri':>10}{'prem tri':>10}{'tri-1/n':>10}{'med log phi':>13}",
     ROWS["A"],
     lambda r: f"{sh(r[0]):<34}{r[1]:>7}{r[2]:>10.3f}{r[3]:>10.3f}{r[4]:>10.3f}"
               f"{r[5]:>10.3f}{r[6]:>10.3f}{r[6] - r[4]:>10.3f}{r[7]:>13.3f}")

show("TABLE O — (A) CROSS-CHECK, capture-OFF: the ORACLE-gDNA route with its two psi'(n_g) SUBTRACTED "
     "should equal the EXACT step==1 route",
     f"{'condition':<34}{'edges':>7}{'Var(oracle)':>12}{'E[P] tri':>10}{'prem(oracle)':>14}"
     f"{'prem(step1)':>13}{'diff':>8}{'of which 2psi_g':>17}",
     ROWS["O"],
     lambda r: f"{sh(r[0]):<34}{r[1]:>7}{r[2]:>12.3f}{r[3]:>10.3f}{r[4]:>14.3f}"
               f"{r[5]:>13.3f}{r[6]:>8.3f}{r[7]:>17.3f}")

show("TABLE B — (A) THE ESTIMATOR.  E1 = ledger's spliced-only gap; E2 = the full graft claim differenced "
     "(the exon's own counts CANCEL there)",
     f"{'condition':<34}{'pairs':>6}{'Var(g1)':>9}{'P 1/n':>7}{'om1 1/n':>9}{'P tri':>7}{'om1 tri':>9}"
     f"{'Var(g2)':>9}{'P2 tri':>8}{'om2 tri':>9}{'T(2j)':>8}{'T(live)':>9}"
     f"{'om1/T2j':>9}{'om2/T2j':>9}{'om1/Tlive':>11}",
     ROWS["B"],
     lambda r: f"{sh(r[0]):<34}{r[1]:>6}{r[2]:>9.3f}{r[3]:>7.3f}{r[4]:>9.3f}{r[5]:>7.3f}{r[6]:>9.3f}"
               f"{r[7]:>9.3f}{r[8]:>8.3f}{r[9]:>9.3f}{r[11]:>8.3f}{r[12]:>9.3f}"
               f"{r[6] / r[11] if r[11] > 0 else np.nan:>9.2f}"
               f"{r[9] / r[11] if r[11] > 0 else np.nan:>9.2f}"
               f"{r[6] / r[12] if r[12] > 0 else np.nan:>11.2f}")

show("TABLE C — (B) THE DIVIDE-BY-2: corr(log phi_left, log phi_right), MEASURED DIRECTLY",
     f"{'condition':<34}{'pairs':>6}{'Var_L':>8}{'Var_R':>8}{'Cov':>8}{'rho_raw':>9}"
     f"{'prem_L':>8}{'prem_R':>8}{'premCov':>9}{'rho_prem':>10}{'2(1-r)raw':>11}"
     f"{'2(1-r)prem':>12}{'om2_corr':>10}{'T(2j)':>8}",
     ROWS["C"],
     lambda r: f"{sh(r[0]):<34}{r[1]:>6}{r[2]:>8.3f}{r[3]:>8.3f}{r[4]:>8.3f}{r[5]:>9.3f}"
               f"{r[6]:>8.3f}{r[7]:>8.3f}{r[8]:>9.3f}{f3(r[9])[1:]:>10}{r[10]:>11.3f}"
               f"{f3(r[11])[:8]:>12}{f3(r[12]):>10}{r[13]:>8.3f}")

show("TABLE L — RECONCILIATION: x1_graft.py's EXACT gap loop (the ledger's published statistic)",
     f"{'condition':<34}{'pairs':>7}{'Var(gap)':>10}{'med|gap|':>10}{'E[1/n]':>9}{'E[psi]':>9}"
     f"{'om 1/n':>9}{'om tri':>9}",
     ROWS["L"],
     lambda r: f"{sh(r[0]):<34}{r[1]:>7}{r[2]:>10.3f}{r[3]:>10.3f}{r[4]:>9.3f}{r[5]:>9.3f}"
               f"{r[6]:>9.3f}{r[7]:>9.3f}")

show("TABLE D — (C) SELECTION: premise Var(log phi) by subset. Is the 2-junction subset representative?\n"
     "          triplets are (n, premise Var, median log phi).",
     f"{'condition':<34}{'Eall':>6}{'Elive':>6}"
     + "".join(f"{a:>7}{b:>8}{c:>8}" for a, b, c in
               (("n_all", "P_all", "med"), ("n_liv", "P_live", "med"), ("n_2j", "P_2j", "med"),
                ("n_1j", "P_1j", "med"), ("n_0spl", "P_0spl", "med"))),
     ROWS["D"],
     lambda r: f"{sh(r[0]):<34}{r[1]:>6}{r[2]:>6}"
               + "".join(f"{r[3 + 3 * j]:>7}{f3(r[4 + 3 * j])[:7]:>8}{f3(r[5 + 3 * j])[:7]:>8}"
                         for j in range(5)))

show("TABLE E — (§5.2 SHAPE) premise variance by the junction's spliced COUNT (trigamma-subtracted)",
     f"{'condition':<34}"
     + "".join(f"{lab:>17}" for lab in ("n<30", "30-100", "100-300", "300-1000", ">1000")),
     ROWS["E"],
     lambda r: f"{sh(r[0]):<34}"
               + "".join(f"{r[1 + 2 * j]:>7}{f3(r[2 + 2 * j])[:8]:>10}" for j in range(5)))

show("TABLE F — (G) candidate observables. SHARE needs the second junction (2j only); COUNT / vs-EXON / "
     "w_mu are single-edge computable.",
     f"{'condition':<34}{'live':>6}{'2j':>5}"
     f"{'share:x':>9}{'share:|x|':>11}{'count:x':>9}{'count:|x|':>11}"
     f"{'vsEXON:x':>10}{'vsEXON:|x|':>12}{'w_mu:x':>9}{'w_mu:|x|':>11}",
     ROWS["F"],
     lambda r: f"{sh(r[0]):<34}{r[1]:>6}{r[2]:>5}"
               + "".join(f3(v)[:9].rjust(w) for v, w in zip(r[3:], (9, 11, 9, 11, 10, 12, 9, 11))))

show("TABLE G — (G) stratified by the junction's two-junction FLUX SHARE: the MEAN of log phi is the BIAS, "
     "the premise variance is the spread",
     f"{'condition':<34}"
     + "".join(f"{lab:>21}" for lab in ("share<0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", ">0.8")),
     ROWS["G"],
     lambda r: f"{sh(r[0]):<34}"
               + "".join(f"{r[1 + 3 * j]:>5}{f3(r[2 + 3 * j])[:7]:>8}{f3(r[3 + 3 * j])[:7]:>8}"
                         for j in range(5)))
print("\nTABLE G column triplets are (n, mean log phi, premise Var).")

show("TABLE I — (I) BIAS-or-NOISE. The share regression's slope/intercept, then the premise Var and median "
     "log phi under FOUR claim forms, on the 2-junction edges.",
     f"{'condition':<34}{'n':>5}{'slope':>7}{'icpt':>7}"
     + "".join(f"{a:>9}{b:>8}" for a, b in (("V:this", "med"), ("V:max", "med"),
                                            ("V:sum", "med"), ("V:2x", "med"))),
     ROWS["I"],
     lambda r: f"{sh(r[0]):<34}{r[1]:>5}{r[2]:>7.2f}{r[3]:>7.2f}"
               + "".join(f"{f3(r[4 + 2 * j])[:9]:>9}{f3(r[5 + 2 * j])[:8]:>8}" for j in range(4)))

show("TABLE H — (G) how much premise variance ONE fitted slope removes.  SHARE is 2j-only; COUNT / vsEXON /"
     " vsBND are single-edge computable on every live graft edge.",
     f"{'condition':<34}{'P(2j)':>8}{'resid|share':>12}{'R2':>7}"
     f"{'P(live)':>9}{'resid|count':>12}{'R2':>7}{'resid|vsEXON':>14}{'R2':>7}"
     f"{'resid|w_mu':>13}{'R2':>7}",
     ROWS["H"],
     lambda r: f"{sh(r[0]):<34}{r[1]:>8.3f}{f3(r[2])[:8]:>12}{f3(r[3])[:7]:>7}"
               f"{r[4]:>9.3f}{f3(r[5])[:8]:>12}{f3(r[6])[:7]:>7}{f3(r[7])[:8]:>14}{f3(r[8])[:7]:>7}"
               f"{f3(r[9])[:8]:>13}{f3(r[10])[:7]:>7}")
