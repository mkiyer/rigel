"""TASK 1 — THE CONTAMINATION PATH.

For every node in the hyperprior's TRAINING SUBSTRATE (region x (single|gonly)), how much of the composition
precision it ends the sweep with came from messages whose relay chain passes through an EXCLUDED node
(AMBIG region, or boundary)?

The accounting is EXACT, and needs no code change, because every relay accumulator is ADDITIVE with a
purely multiplicative per-hop damping:

    tau[i]  = tau_own[i] + ttau,     ttau = _damp(tau[s], s2t)          (bp_solver `_relay`)
    pg[i]   = pg_own[i]  + tpg,      tpg  = _damp(pg[s],  s2t)

so  delivered_i = accum_i - own_i  and the tainted part of the delivered precision is
delivered_i * share[s].  Hence the recursion along the chain order

    taint_i = T[i]*own_i + (accum_i - own_i) * share[s] ,   share_i = taint_i / accum_i

is the share-preserving apportionment, and it is exact wherever the stream has no non-relay addition
(true for tau and for pg; the graft adds a spliced COUNT to the RNA streams only, and that count's source
is the boundary itself, i.e. also excluded).

The combine's per-message composition precisions are read directly off the solver's own `_dl` capture
(`tau_post` == the atau / btau the combine fuses), so nothing here is re-derived offline.

    OMP_NUM_THREADS=1 python scratchpad/ho_taint.py [cond ...]
"""

from __future__ import annotations

import dataclasses
import importlib
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from flagship_interrogate import _oracle_per_node  # noqa: E402
from selfsolve_diag import _scan_and_truth  # noqa: E402

from rigel.calibration.bp_solver import REGION  # noqa: E402
from rigel.calibration.node_geometry import _node_region_type  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
_EPS = 1e-9
CONDS = sys.argv[1:] or [
    "gdna_gdna300_ss_0.99_nrna_present_capture_off",
    "gdna_gdna300_ss_0.99_nrna_present_capture_on",
    "gdna_gdna300_ss_0.50_nrna_present_capture_on",
    "gdna_gdna100_ss_0.50_nrna_present_capture_verystrong",
]

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)


def _share(accum, own, order, nbr, T):
    """share[i] = tainted fraction of node i's relay accumulator, along one direction."""
    n = accum.shape[0]
    sh = np.zeros(n)
    for i in order:
        s = nbr[i]
        up = 0.0
        if s >= 0 and accum[s] > 0.0:
            up = max(accum[i] - own[i], 0.0) * sh[s]
        t = (own[i] if T[i] else 0.0) + up
        sh[i] = t / accum[i] if accum[i] > _EPS else 0.0
    return np.clip(sh, 0.0, 1.0)


ALLROWS = []
for cond in CONDS:
    inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
    dbg: dict = {}
    calmod.calibrate(
        inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
        np.asarray(inp["rna_fl_pmf"]),
        dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg,
    )
    chain, cap, st = dbg["chain"], dbg["capture"], dbg["statics"]
    us, uni = cap["_uni_static"], cap["_uni"][-1]
    dls = cap["_dl"][-2:]
    dl_l = next(d for d in dls if d["df"] == 0)
    dl_r = next(d for d in dls if d["df"] == 1)

    isr = np.asarray(chain.kind) == REGION
    fp, fn = np.asarray(st.free_pos, bool), np.asarray(st.free_neg, bool)
    amb, single, gonly = fp & fn, fp ^ fn, ~fp & ~fn
    mass = np.asarray(cap["mass_global"], float)
    eff = np.asarray(cap["eff_global"], float)
    live = (eff > 1e-9) & (mass > 1e-12)
    SUB = live & isr & (single | gonly)          # the hyperprior's training substrate
    EXCL = ~(isr & (single | gonly))             # AMBIG regions + every boundary
    rt, _ = _node_region_type(chain, ra)

    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
    T_ = Gp + Gn + Rp + Rn
    fo = np.where(T_ > _EPS, (Gp + Gn) / np.maximum(T_, _EPS), np.nan)
    fg = np.asarray(cap["f_g"], float)
    err = np.abs(fg - fo)

    order = us["order"].tolist()
    left, right = us["left"], us["right"]
    tau_own, pg_own = us["tau_own"], us["pg_own"]

    out = {}
    for tag, T in (("EXCL(ambig+bnd)", EXCL), ("AMBIG-region", isr & amb), ("boundary", ~isr)):
        shf_t = _share(us["fwd_tau"], tau_own, order, left, T)
        shb_t = _share(us["bwd_tau"], tau_own, order[::-1], right, T)
        shf_p = _share(us["fwd_pg"], pg_own, order, left, T)
        shb_p = _share(us["bwd_pg"], pg_own, order[::-1], right, T)
        sl, sr = np.clip(left, 0, len(order) - 1), np.clip(right, 0, len(order) - 1)
        at, bt = dl_l["tau_post"], dl_r["tau_post"]
        ct = at + bt
        tau_msg_taint = np.where(ct > _EPS, (at * shf_t[sl] + bt * shb_t[sr]) / np.maximum(ct, _EPS), 0.0)
        ap, bp = uni["apg"], uni["bpg"]
        cp = ap + bp
        pg_msg_taint = np.where(cp > _EPS, (ap * shf_p[sl] + bp * shb_p[sr]) / np.maximum(cp, _EPS), 0.0)
        # message-borne share of the node's TOTAL composition precision
        msg_frac_tau = np.where(ct + tau_own > _EPS, ct / np.maximum(ct + tau_own, _EPS), 0.0)
        msg_frac_pg = np.where(cp + pg_own > _EPS, cp / np.maximum(cp + pg_own, _EPS), 0.0)
        out[tag] = (tau_msg_taint, pg_msg_taint, msg_frac_tau, msg_frac_pg, ct, cp)

    print(f"\n{'=' * 118}\n{cond}\n{'=' * 118}")
    print(f"  substrate nodes={int(SUB.sum()):,}  mass={mass[SUB].sum():,.0f}   "
          f"excluded nodes={int((EXCL & live).sum()):,}  mass={mass[EXCL & live].sum():,.0f}")
    hdr = (f"  {'substrate stratum':<22}{'n':>7}{'msg%tau':>9}{'msg%pg':>8}"
           f"{'taint tau':>11}{'taint pg':>10}{'  |  tau taint by src: excl/ambR/bnd':<38}{'mwae':>8}")
    print(hdr)
    rows = [("ALL substrate", SUB)]
    for k, nm in ((2, "exon"), (1, "intron"), (0, "intergenic")):
        rows.append(("  " + nm, SUB & (rt == k)))
    rows.append(("  single-strand", SUB & single))
    rows.append(("  gonly(structural)", SUB & gonly))
    for lab, m in rows:
        m = m & np.isfinite(fo)
        if not m.any():
            continue
        w = mass[m]
        tt, tp, mft, mfp, ct, cp = out["EXCL(ambig+bnd)"]
        a_t = out["AMBIG-region"][0][m]
        b_t = out["boundary"][0][m]
        print(f"  {lab:<22}{int(m.sum()):>7}{np.average(mft[m], weights=w):>9.1%}"
              f"{np.average(mfp[m], weights=w):>8.1%}"
              f"{np.average(tt[m], weights=w):>11.1%}{np.average(tp[m], weights=w):>10.1%}"
              f"  |  {np.average(tt[m], weights=w):>6.1%} /{np.average(a_t, weights=w):>7.1%} /"
              f"{np.average(b_t, weights=w):>7.1%}      "
              f"{np.average(err[m], weights=w):>8.4f}")

    # ── does the contamination CORRELATE with the substrate node's own error? ─────────────────────────
    tt = out["EXCL(ambig+bnd)"][0]
    at_ = out["AMBIG-region"][0]
    m = SUB & np.isfinite(fo) & (mass > _EPS)
    print(f"\n  contamination vs |f_g - oracle| on the substrate  (n={int(m.sum()):,})")
    for nm, x in (("taint share (excl)", tt), ("taint share (AMBIG-region only)", at_)):
        q = np.quantile(x[m], [0.25, 0.5, 0.75])
        print(f"    {nm:<34} quartiles {q[0]:.3f} / {q[1]:.3f} / {q[2]:.3f}")
        b = np.digitize(x[m], [0.001, 0.05, 0.2, 0.5])
        line = "      by taint bin  "
        for k in range(5):
            s = b == k
            if s.sum() < 5:
                line += f"{'-':>13}"
                continue
            line += f"{np.average(err[m][s], weights=mass[m][s]):>8.4f}({s.sum():>4})"
        print(line)
        xc, yc = x[m], err[m]
        if xc.std() > 1e-9:
            print(f"      corr(taint, |err|) = {np.corrcoef(xc, yc)[0, 1]:+.3f}")
    ALLROWS.append((cond, float(np.average(tt[m], weights=mass[m])),
                    float(np.average(at_[m], weights=mass[m]))))

print(f"\n{'=' * 70}\nSUMMARY — substrate composition-message taint share (mass-weighted)")
print(f"  {'condition':<50}{'excl':>9}{'ambigR':>9}")
for c, a, b in ALLROWS:
    print(f"  {c[5:]:<50}{a:>9.1%}{b:>9.1%}")
