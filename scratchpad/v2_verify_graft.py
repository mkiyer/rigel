"""ADVERSARIAL VERIFICATION, round 2.

  W1. T4's DENOMINATOR. The report compares Var(log phi) (a premise about the WHOLE RNA claim, per its
      own finding 4) against the SPLICED ARM's own count variance 1/SP_mass. If the term belongs on
      the whole claim, the honest comparison is against the DELIVERED RNA variance 1/(tpp+tpn).
      Measure both, and the resulting z2.
  W2. T6's DENOMINATOR. The report uses the FUSED c_tau (both messages) against a PER-MESSAGE premise
      variance. Recompute with the per-message tau (from the _dl capture's tau_post) and with the
      pre-gate ttau.
  W3. T3b ROBUSTNESS. Var(log phi) by w_mu bin is a non-robust statistic on a heavy tail (their own
      caveat 2). Repeat with a robust scale (IQR/1.349)^2 and a bootstrap CI on the variance.
  W4. FIRING COUNT. Combine grafts over ALL rho-iterations vs the relay's.
  W5. Does the graft's own count precision reach the destination at all, or is it swamped?

    OMP_NUM_THREADS=1 python scratchpad/v2_verify_graft.py
"""

from __future__ import annotations

import dataclasses
import importlib
import sys
from pathlib import Path

import numpy as np

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
rng = np.random.default_rng(0)


def qbin(x, k=5):
    q = np.quantile(x, np.linspace(0.0, 1.0, k + 1))
    q[0] -= 1e-12
    q[-1] += 1e-12
    return np.clip(np.digitize(x, q[1:-1]), 0, k - 1), q


def boot_var(x, B=2000):
    idx = rng.integers(0, x.size, size=(B, x.size))
    v = np.var(x[idx], axis=1)
    return np.percentile(v, 2.5), np.percentile(v, 97.5)


for cond in CONDS:
    print("=" * 108)
    print(f"CONDITION {cond}")
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
    E_g, E_r = us["E_g"], us["E_r"]
    li, ri = us["left"], us["right"]
    is_bnd, is_exon = us["is_bnd"], us["is_exon"]
    SP, SN = (us["SP_l"], us["SP_r"]), (us["SN_l"], us["SN_r"])
    ESP = (np.asarray(geo.eff_spl_left, float), np.asarray(geo.eff_spl_right, float))
    pin = cap["_pin"][-2:]
    uni = cap["_uni"][-1]
    dl = cap["_dl"][-2:] if "_dl" in cap else None

    # ── W4 firing count ───────────────────────────────────────────────────────────────────────
    tot_all = sum(int(p["graft"].sum()) for p in cap["_pin"])
    tot_last = sum(int(p["graft"].sum()) for p in pin)
    relay = int((is_exon & (li >= 0) & is_bnd[np.clip(li, 0, n - 1)]).sum()
                + (is_exon & (ri >= 0) & is_bnd[np.clip(ri, 0, n - 1)]).sum())
    print(f"  W4. graft firings/sweep: combine ALL rho-iters = {tot_all}  (last iter only = {tot_last})"
          f"   relay = {relay}   -> XVAR reach = {100*tot_all/(tot_all+relay):.0f}%"
          f"   (report said {100*tot_last/(tot_last+relay):.0f}%)")

    rho_g = np.where(E_g > _EPS, G / np.maximum(E_g, _EPS), np.nan)
    rho_R = np.where(E_r > _EPS, Ru / np.maximum(E_r, _EPS), np.nan)

    rows = []
    for df, nbr in ((0, li), (1, ri)):
        face = 1 - df
        pe = pin[df]
        for i in np.flatnonzero(is_exon):
            b = nbr[i]
            if b < 0 or not is_bnd[b]:
                continue
            mu = (SP[face][b] + SN[face][b]) / max(ESP[face][b], _EPS)
            smass = SP[face][b] + SN[face][b]
            if not (mu > _EPS) or not (rho_R[i] > _EPS):
                continue
            if not (rho_g[b] > _EPS and rho_g[i] > _EPS) or not (rho_R[b] >= 0):
                continue
            ph = (rho_R[b] + mu) / (rho_R[i] * rho_g[b] / rho_g[i])
            if not (np.isfinite(ph) and ph > 0):
                continue
            tp_tot = pe["tpp"][i] + pe["tpn"][i]       # DELIVERED RNA mode precision (pre-DL)
            rows.append((np.log(ph), smass, tp_tot, pe["spl_prec"][i],
                         mu / max(rho_R[b] + mu, _EPS), df, i))
    A = np.asarray(rows, float)
    lp, smass, tptot, spc, w_mu = A[:, 0], A[:, 1], A[:, 2], A[:, 3], A[:, 4]

    # ── W1 the honest over-confidence ─────────────────────────────────────────────────────────
    ok = tptot > 0
    v_spl = 1.0 / np.maximum(smass, _EPS)
    v_del = 1.0 / np.maximum(tptot[ok], _EPS)
    print(f"  W1. n={A.shape[0]}  Var(log phi)={np.var(lp):.4f}")
    print(f"      vs the SPLICED ARM's own v (report's T4): med v={np.median(v_spl):.5f}"
          f"  -> {np.var(lp)/np.median(v_spl):>7.1f}x   [REPORT]")
    print(f"      vs the DELIVERED RNA claim's v = 1/(tpp+tpn): med v={np.median(v_del):.5f}"
          f"  -> {np.var(lp)/np.median(v_del):>7.1f}x   [HONEST, per their own finding 4]")
    print(f"      per-edge z2 = E[(log phi)^2]*E[tpp+tpn] "
          f"= {np.mean(lp[ok]**2)*np.mean(tptot[ok]):.1f} ; "
          f"median-based = {np.var(lp)*np.median(tptot[ok]):.1f}")
    bi, q = qbin(smass[ok])
    lpo, vo = lp[ok], v_del
    print(f"      {'spl MASS bin':<22}{'edges':>7}{'Var(logphi)':>13}{'med v_spl':>11}"
          f"{'med v_deliv':>13}{'oc SPL':>9}{'oc DELIV':>10}")
    for k in range(5):
        m = bi == k
        if m.sum() < 5:
            continue
        vv = np.var(lpo[m])
        print(f"      [{q[k]:>8.1f},{q[k+1]:>8.1f}]{'':<3}{int(m.sum()):>7d}{vv:>13.3f}"
              f"{np.median(v_spl[ok][m]):>11.5f}{np.median(vo[m]):>13.5f}"
              f"{vv/np.median(v_spl[ok][m]):>8.0f}x{vv/np.median(vo[m]):>9.0f}x")

    # ── W2 T6 per-message tau ─────────────────────────────────────────────────────────────────
    ct = np.asarray(uni["c_tau"], float)
    exg = pin[0]["graft"] | pin[1]["graft"]
    live = exg & (ct > _EPS)
    print(f"  W2. FUSED 1/c_tau at graft-fed exons  med={np.median(1.0/np.maximum(ct[live],_EPS)):.4f}"
          f"  n={int(live.sum())}   -> report's over-conf {np.var(lp)*np.median(ct[live]):.0f}x")
    if dl is not None:
        for d in dl:
            tp_ = np.asarray(d["tau_post"], float)
            lv = exg & (tp_ > _EPS)
            print(f"      df={d['df']} PER-MESSAGE 1/tau_post med="
                  f"{np.median(1.0/np.maximum(tp_[lv],_EPS)):.4f} n={int(lv.sum())}"
                  f"   -> over-conf {np.var(lp)*np.median(tp_[lv]):.0f}x")
    print(f"      frac of graft-fed exons with c_tau == 0 (NO lambda claim at all): "
          f"{float(np.mean(ct[exg] <= _EPS)):.3f}")

    # ── W3 T3b robustness ─────────────────────────────────────────────────────────────────────
    if np.percentile(w_mu, 90) - np.percentile(w_mu, 10) > 0.05:
        print("  W3. Var(log phi) by w_mu quintile -- with bootstrap CI and a ROBUST scale")
        bi, q = qbin(w_mu)
        print(f"      {'w_mu bin':<20}{'n':>6}{'med w_mu':>10}{'Var':>9}{'boot 95% CI':>22}"
              f"{'robust (IQR/1.349)^2':>22}{'rob/w_mu^2':>12}")
        for k in range(5):
            m = bi == k
            if m.sum() < 5:
                continue
            x = lp[m]
            lo, hi = boot_var(x)
            iqr = np.percentile(x, 75) - np.percentile(x, 25)
            rob = (iqr / 1.349) ** 2
            wm = float(np.median(w_mu[m]))
            print(f"      [{q[k]:>6.3f},{q[k+1]:>6.3f}]{'':<3}{int(m.sum()):>6d}{wm:>10.3f}"
                  f"{np.var(x):>9.3f}   [{lo:>7.3f},{hi:>7.3f}]{rob:>22.3f}"
                  f"{rob/max(wm*wm,1e-9):>12.2f}")
    else:
        print(f"  W3. w_mu has no spread here (p10={np.percentile(w_mu,10):.3f} "
              f"p90={np.percentile(w_mu,90):.3f}) -- test unavailable")
