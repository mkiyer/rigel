"""P1d / P1e FIRING DIAGNOSTIC — characterise both terms on REAL cfRNA (and one toy for contrast).

No oracle. Everything here is either read straight out of ``_debug["capture"]`` or recomputed from it
with the EXACT expressions the solver uses (`bp_solver._transport`'s P1e block, `graft_premise_logvar`).

    OMP_NUM_THREADS=1 python scratchpad/p1de_firing.py --input cfrna:LBX0190
    OMP_NUM_THREADS=1 python scratchpad/p1de_firing.py --input toy
    OMP_NUM_THREADS=1 python scratchpad/p1de_firing.py --input cfrna:LBX0190 --no-p1d

``--no-p1d`` monkeypatches ``bp_solver.graft_premise_logvar`` to return a POOLED value of 0 (P1d
neutralised, per-edge diagnostic array untouched) and restores it afterwards. Nothing under src/ is
modified.
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import os
import pickle
import sys
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")

from rigel.config import PipelineConfig  # noqa: E402

import rigel.calibration.bp_solver as bps  # noqa: E402
from rigel.calibration.calibrate import calibrate  # noqa: E402
from rigel.calibration.node_chain import build_node_chain  # noqa: E402
from rigel.calibration.node_geometry import _node_region_type  # noqa: E402

CF = Path("/Users/mkiyer/Downloads/rigel_runs/cfrna/_calib_cache")
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
TOY_COND = "gdna_gdna300_ss_0.50_nrna_present_capture_off"
EPS = 1.0e-9  # bp_solver._EPS


# ───────────────────────────────────────────────────────────────────────────── loading

def load(spec: str):
    if spec.startswith("cfrna:"):
        d = pickle.load(open(CF / f"{spec.split(':', 1)[1]}.pkl", "rb"))
        return (
            d["payload"],
            d["region_arrays"],
            d["strand_model"],
            np.asarray(d["gdna_fl_pmf"]),
            np.asarray(d["rna_fl_pmf"]),
        )
    from selfsolve_diag import _scan_and_truth

    from rigel.calibration.region_arrays import RegionArrays
    from rigel.index import TranscriptIndex

    index = TranscriptIndex.load(str(SUITE / "rigel_index"))
    cfg = PipelineConfig()
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    cond = spec.split(":", 1)[1] if ":" in spec else TOY_COND
    inp = _scan_and_truth(
        SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache"
    )
    return (
        inp["payload"],
        ra,
        inp["strand_model"],
        np.asarray(inp["gdna_fl_pmf"]),
        np.asarray(inp["rna_fl_pmf"]),
    )


# ───────────────────────────────────────────────────────────────────────────── helpers

def q(x, ps=(10, 25, 50, 75, 90, 99)):
    x = np.asarray(x, np.float64)
    if x.size == 0:
        return {f"p{p}": float("nan") for p in ps}
    return {f"p{p}": float(np.percentile(x, p)) for p in ps}


def fmtq(d):
    return "  ".join(f"{k}={v:.4g}" for k, v in d.items())


def top_decile_share(x):
    x = np.asarray(x, np.float64)
    if x.size == 0 or x.sum() <= 0:
        return float("nan")
    k = max(1, int(np.ceil(0.1 * x.size)))
    return float(np.sort(x)[::-1][:k].sum() / x.sum())


# ───────────────────────────────────────────────────────────────────────────── P1e replay

def p1e_replay(pin, us, n_node):
    """Reproduce `_transport`'s P1e block EXACTLY for one captured message direction.

    Returns a dict of per-destination-node arrays. Every line below is a transcription of
    bp_solver.py lines 1051-1102.
    """
    tg, tp, tn = pin["tg"], pin["tp"], pin["tn"]
    tpg, tpp, tpn = pin["tpg"], pin["tpp"], pin["tpn"]
    valid, src, s2t = pin["valid"], pin["src"], pin["s2t"]
    M, E_g, E_r = us["M"], us["E_g"], us["E_r"]
    og, op, on = us["og"], us["op"], us["on"]

    _p3 = np.stack([tpg, tpp, tpn], axis=-1)
    _sup = _p3 > 0.0
    _mc = np.where(_sup, np.stack([tg, tp, tn], axis=-1), np.stack([og, op, on], axis=-1))
    _mc = _mc * np.stack([E_g, E_r, E_r], axis=-1)
    _S = _mc.sum(axis=-1)
    _okc = valid & (_S > EPS) & (M > EPS)
    _al = _mc / np.maximum(_S, EPS)[..., None]
    _vc = np.where(_sup, 1.0 / np.maximum(_p3, EPS), 0.0)
    _s2c = (np.where(np.isfinite(s2t), s2t, 0.0) + 1.0 / np.maximum(pin["n_src"], EPS))[..., None]
    _sv = np.where(_sup, _s2c + _al * np.maximum(_vc - _s2c, 0.0), 0.0)
    _aSa = np.sum(_al * _sv, axis=-1)
    _dlt = np.where(_okc, np.log(np.maximum(M, EPS) / np.maximum(_S, EPS)), 0.0)
    _den = _aSa + 1.0 / np.maximum(n_node, EPS)
    _b2 = np.maximum(_dlt * _dlt - _den, 0.0)
    _b2 = np.where(_dlt < 0.0, _b2, 0.0)
    _b2 = np.where(_sup.any(axis=-1), _b2, 0.0)
    return {
        "sup_any": _sup.any(axis=-1),
        "okc": _okc,
        "S": _S,
        "delta": _dlt,
        "aSa": _aSa,
        "den": _den,
        "b2": _b2,
        "prec": np.stack([tpg, tpp, tpn], axis=-1),
        "src": src,
        "valid": valid,
        "graft": pin["graft"],
    }


# ───────────────────────────────────────────────────────────────────────────── main

ap = argparse.ArgumentParser()
ap.add_argument("--input", default="toy")
ap.add_argument("--refit", type=int, default=1)
ap.add_argument("--no-p1d", action="store_true")
ap.add_argument("--json-out", default="")
a = ap.parse_args()

payload, ra, sm, gpmf, rpmf = load(a.input)
cc = dataclasses.replace(PipelineConfig().calibration, calib_refit_iters=a.refit)

_orig_glv = bps.graft_premise_logvar
if a.no_p1d:
    bps.graft_premise_logvar = lambda *x, **k: (_orig_glv(*x, **k)[0], 0.0)

# VALIDATION HOOK. `mismatch_deflate` is called ONLY inside `_transport`, 7× per invocation, and the
# first three calls receive tpg / tpp / tpn AFTER the P1e damping. Recording them lets us check the
# replay below against the live solver EXACTLY, rather than trusting a transcription.
_orig_md = bps.mismatch_deflate
_md_seen: list = []


def _md_spy(pr, *x, **k):
    if len(_md_seen) % 7 < 3:
        _md_seen.append(np.asarray(pr, np.float64).copy())
    else:
        _md_seen.append(None)
    return _orig_md(pr, *x, **k)


bps.mismatch_deflate = _md_spy

dbg: dict = {}
t0 = time.perf_counter()
try:
    calibrate(payload, ra, sm, gpmf, rpmf, cc, _debug=dbg)
finally:
    bps.graft_premise_logvar = _orig_glv
    bps.mismatch_deflate = _orig_md
wall = time.perf_counter() - t0

cap = dbg["capture"]
us = cap["_uni_static"]
glv = cap["_glv"]
pins = cap["_pin"]

is_bnd = np.asarray(us["is_bnd"], bool)
is_reg = ~is_bnd
is_exon = np.asarray(us["is_exon"], bool)
fp = np.asarray(cap["free_pos"], bool)
fn = np.asarray(cap["free_neg"], bool)
n_node = np.where(is_reg, us["n_unspl_l"], us["n_unspl_l"] + us["n_unspl_r"])
n = is_bnd.size

chain = build_node_chain(payload.ref_region_offsets, payload.ref_boundary_offsets)
ntype, _ = _node_region_type(chain, ra)  # 0 intergenic, 1 intron, 2 exon, -1 boundary
kindname = np.array(["intergenic", "intron", "exon", "bnd"])[np.where(is_bnd, 3, ntype)]
strandname = np.where(fp & fn, "ambig", np.where(fp | fn, "single", "none"))
cls = np.char.add(np.char.add(kindname, "/"), strandname)

hdr = f"{a.input}  refit={a.refit}  P1d={'OFF' if a.no_p1d else 'ON'}  wall={wall:.1f}s"
print("=" * 100)
print(hdr)
print(f"nodes={n}  region={int(is_reg.sum())}  bnd={int(is_bnd.sum())}  exon={int(is_exon.sum())}")
print(f"strand: single={int((strandname == 'single').sum())} ambig={int((strandname == 'ambig').sum())}"
      f" none={int((strandname == 'none').sum())}")
print("=" * 100)

out: dict = {"input": a.input, "no_p1d": a.no_p1d, "n_nodes": int(n)}

# ═════════════════════════════════════════════════════════════════════ P1d
print("\n### P1d  (graft_premise_logvar)")
graft_any = np.zeros(n, bool)
n_graft_edges = 0
graft_ok_frac = {}
for pin in pins:
    graft_any |= pin["graft"]
    n_graft_edges += int(pin["graft"].sum())
print(f"graft EDGES (over both transports) = {n_graft_edges};  distinct graft-receiving nodes = "
      f"{int(graft_any.sum())}")
out["n_graft_edges"] = n_graft_edges
out["n_graft_nodes"] = int(graft_any.sum())
_ok_any = np.zeros(n, bool)
for g in glv:
    _ok_any |= np.asarray(g["ok"], bool)
_ge_any = sum(int((pin["graft"] & _ok_any).sum()) for pin in pins)
print(f"graft EDGES with a live pair on EITHER strand = {_ge_any} / {n_graft_edges} = "
      f"{_ge_any / max(1, n_graft_edges):.4f}")
out["graft_edge_pair_frac_either"] = _ge_any / max(1, n_graft_edges)

for g in glv:
    s = "pos" if g["strand"] == 0 else "neg"
    ok = np.asarray(g["ok"], bool)
    d = np.asarray(g["d"], np.float64)
    nz = np.asarray(g["noise"], np.float64)
    d2 = d[ok] ** 2
    nz_ok = nz[ok]
    # Q1: n_pairs as a fraction of GRAFT EDGES (ω is applied per graft edge, indexed by destination)
    ge_ok = sum(int((pin["graft"] & ok).sum()) for pin in pins)
    print(f"\n-- strand {s}: omega = {g['omega']:.6g}   n_pairs = {g['n_pairs']}")
    print(f"   n_pairs / n_region_nodes      = {g['n_pairs'] / max(1, int(is_reg.sum())):.4f}")
    print(f"   n_pairs / n_exon_nodes        = {g['n_pairs'] / max(1, int(is_exon.sum())):.4f}")
    print(f"   graft EDGES with a live pair  = {ge_ok} / {n_graft_edges} = "
          f"{ge_ok / max(1, n_graft_edges):.4f}   <<< Q1")
    print("   pair-node kind mix: " + "  ".join(
        f"{k}={int((kindname[ok] == k).sum())}" for k in ("exon", "intron", "intergenic", "bnd")))
    # Q2: noise subtraction
    Ed2, En = float(g["Ed2"]), float(g["Enoise"])
    surv = (Ed2 - En) / Ed2 if Ed2 > 0 else float("nan")
    print(f"   E[d^2] = {Ed2:.6g}   E[noise] = {En:.6g}   surviving share = {surv:.4f}   <<< Q2")
    print(f"   (omega = max(0, E[d2]-E[noise])/2 = {max(0.0, Ed2 - En) * 0.5:.6g})")
    print(f"   d^2   {fmtq(q(d2))}")
    print(f"   noise {fmtq(q(nz_ok))}")
    trunc = float((d2 <= nz_ok).mean()) if d2.size else float("nan")
    zero = float((d2 <= 0.0).mean()) if d2.size else float("nan")
    print(f"   per-pair truncation rate (d^2 <= noise) = {trunc:.4f}   (d == 0 exactly: {zero:.4f})")
    # Q3: concentration
    tds = top_decile_share(d2)
    # inverse participation ratio of d^2 — how many pairs the mean E[d^2] is EFFECTIVELY built from
    ipr = float(d2.sum() ** 2 / max((d2 ** 2).sum(), 1e-300)) if d2.size else float("nan")
    print(f"   top-decile share of sum(d^2) = {tds:.4f}   <<< Q3")
    print(f"   top-1% share of sum(d^2)     = "
          f"{float(np.sort(d2)[::-1][:max(1, d2.size // 100)].sum() / max(d2.sum(), 1e-300)):.4f}")
    print(f"   effective n behind E[d^2] (IPR) = {ipr:.1f}  of {d2.size}  "
          f"({ipr / max(1, d2.size):.4f})")
    out[f"p1d_{s}"] = {
        "omega": float(g["omega"]), "n_pairs": int(g["n_pairs"]),
        "graft_edge_pair_frac": ge_ok / max(1, n_graft_edges),
        "Ed2": Ed2, "Enoise": En, "surviving_share": surv,
        "trunc_rate": trunc, "zero_rate": zero, "top_decile_share": tds, "ipr": ipr,
        "d2_q": q(d2), "noise_q": q(nz_ok),
    }

# ═════════════════════════════════════════════════════════════════════ P1e
print("\n\n### P1e  (conservation surprise)")
rep = [p1e_replay(pin, us, n_node) for pin in pins]

VAL = np.concatenate([r["valid"] for r in rep])
SUP = np.concatenate([r["sup_any"] for r in rep])
OKC = np.concatenate([r["okc"] for r in rep])
B2 = np.concatenate([r["b2"] for r in rep])
DLT = np.concatenate([r["delta"] for r in rep])
ASA = np.concatenate([r["aSa"] for r in rep])
DEN = np.concatenate([r["den"] for r in rep])
PREC = np.concatenate([r["prec"] for r in rep], axis=0)
GRAFT = np.concatenate([r["graft"] for r in rep])
SRC = np.concatenate([r["src"] for r in rep])
CLS = np.concatenate([cls, cls])
PEEL = np.concatenate([
    (is_bnd & is_exon[pin["src"]] & pin["valid"]) for pin in pins
])
EDGE = np.where(GRAFT, "graft", np.where(PEEL, "peel", "plain"))

live = VAL & SUP
fire = B2 > 0.0
# precision REMOVED by the damping, summed over the three mode-fusion components (tpg/tpp/tpn).
# (the measurement-stream precisions tmg/tmp/tmn are damped identically but are not in the capture)
rem = np.where(fire[:, None], PREC - PREC / (1.0 + PREC * B2[:, None]), 0.0)
REM = rem.sum(axis=1)
PSUM = PREC.sum(axis=1)

# ── VALIDATION: replay vs the LIVE solver, exactly. ``_md_seen``'s last len(pins)*7 entries belong
# to the captured (last) sweep; entries 7k+0/1/2 are tpg/tpp/tpn AFTER the P1e damping.
_tail = _md_seen[-len(pins) * 7:]
_werr = 0.0
for _k, _pin in enumerate(pins):
    _pre = np.stack([_pin["tpg"], _pin["tpp"], _pin["tpn"]], axis=-1)
    _b2k = rep[_k]["b2"]
    _exp = _pre / (1.0 + _pre * _b2k[:, None])
    _got = np.stack(_tail[7 * _k: 7 * _k + 3], axis=-1)
    _werr = max(_werr, float(np.max(np.abs(_exp - _got) / np.maximum(np.abs(_got), 1e-12))))
print(f"[replay validation] max rel |replayed damped precision - live| = {_werr:.3e}  "
      f"({'EXACT' if _werr < 1e-12 else 'MISMATCH — do not trust the P1e numbers'})")
out["replay_max_rel_err"] = _werr

print(f"messages: total={VAL.size}  valid={int(VAL.sum())}  valid&supplied={int(live.sum())}  "
      f"valid&okc={int((VAL & OKC).sum())}")
print("\nQ4  firing rate")
print(f"   b2>0 / valid            = {int(fire.sum())} / {int(VAL.sum())} = "
      f"{fire.sum() / max(1, VAL.sum()):.4f}")
print(f"   b2>0 / valid&supplied   = {int(fire.sum())} / {int(live.sum())} = "
      f"{fire.sum() / max(1, live.sum()):.4f}   <<< Q4")
print(f"   delta<0 / valid&okc     = "
      f"{float((DLT[VAL & OKC] < 0).mean()) if (VAL & OKC).any() else float('nan'):.4f}  "
      f"(the scope gate)")
_d0 = (DLT < 0) & live
print(f"   of delta<0 messages, share that clear the noise floor (d^2>den) = "
      f"{float(fire[_d0].mean()) if _d0.any() else float('nan'):.4f}")
print(f"   precision on firing msgs / precision on all live msgs = "
      f"{PSUM[fire].sum() / max(PSUM[live].sum(), 1e-300):.4f}")
print(f"   TOTAL precision removed / total live precision        = "
      f"{REM.sum() / max(PSUM[live].sum(), 1e-300):.4f}   <<< damping mass share")
print(f"   sum(removed) = {REM.sum():.6g}   sum(live precision) = {PSUM[live].sum():.6g}")
print(f"   top-decile of firing msgs carries {top_decile_share(REM[fire]):.4f} of the removal")
out["p1e"] = {
    "n_valid": int(VAL.sum()), "n_live": int(live.sum()), "n_fire": int(fire.sum()),
    "fire_rate_live": float(fire.sum() / max(1, live.sum())),
    "fire_rate_valid": float(fire.sum() / max(1, VAL.sum())),
    "removed": float(REM.sum()), "prec_live": float(PSUM[live].sum()),
    "removed_frac": float(REM.sum() / max(PSUM[live].sum(), 1e-300)),
    "prec_on_firing_frac": float(PSUM[fire].sum() / max(PSUM[live].sum(), 1e-300)),
}


def bias_row(mask, label):
    m = mask & fire
    nm = int(m.sum())
    if nm == 0:
        print(f"   {label:<26} n=0")
        return None
    d = DLT[m]
    mu, sd = float(d.mean()), float(d.std())
    bs = mu * mu / max(mu * mu + sd * sd, 1e-300)
    w = REM[m].sum()
    print(f"   {label:<26} n={nm:<9} E[d]={mu:+.4f}  sd={sd:.4f}  bias={bs:.3f}  "
          f"E[b2]={float(B2[m].mean()):.4g}  rem={w:.4g} ({w / max(REM.sum(), 1e-300) * 100:5.1f}%)")
    return {"n": nm, "Ed": mu, "sd": sd, "bias_share": bs, "removed": float(w)}


print("\nQ5  bias share of delta  (bias = E[d]^2 / E[d^2]);  'rem' = precision removed")
print("  -- over FIRING messages (delta<0 by construction — read with that selection in mind)")
out["p1e_bias_fire"] = {"ALL": bias_row(np.ones_like(fire), "ALL firing")}
for e in ("graft", "peel", "plain"):
    out["p1e_bias_fire"][e] = bias_row(EDGE == e, f"edge={e}")
for c in sorted(set(CLS.tolist())):
    r = bias_row(CLS == c, f"dst={c}")
    if r:
        out["p1e_bias_fire"][f"dst={c}"] = r

print("  -- over ALL live+okc messages (UNCONDITIONAL: is delta systematic before the scope gate?)")
base = live & OKC
for label, m in (("ALL", base), ("graft", base & (EDGE == "graft")),
                 ("peel", base & (EDGE == "peel")), ("plain", base & (EDGE == "plain"))):
    if m.sum() == 0:
        continue
    d = DLT[m]
    mu, sd = float(d.mean()), float(d.std())
    print(f"   {label:<26} n={int(m.sum()):<9} E[d]={mu:+.4f}  sd={sd:.4f}  "
          f"bias={mu * mu / max(mu * mu + sd * sd, 1e-300):.3f}  "
          f"P(d<0)={float((d < 0).mean()):.3f}")
    out.setdefault("p1e_bias_all", {})[label] = {
        "n": int(m.sum()), "Ed": mu, "sd": sd,
        "bias_share": mu * mu / max(mu * mu + sd * sd, 1e-300),
        "p_neg": float((d < 0).mean()),
    }

print("\nQ6  magnitude of b2 and the attenuation  p -> p/(1+p*b2)")
print(f"   b2 over firing    {fmtq(q(B2[fire]))}")
print(f"   |delta| firing    {fmtq(q(np.abs(DLT[fire])))}")
print(f"   den(=aSa+1/n_dst) {fmtq(q(DEN[fire]))}")
print(f"   aSa  firing       {fmtq(q(ASA[fire]))}")
out["p1e_b2_q"] = q(B2[fire])
out["p1e_delta_q"] = q(np.abs(DLT[fire]))
for j, cname in enumerate(("gDNA", "RNA+", "RNA-")):
    pj = PREC[:, j]
    m = fire & (pj > 0)
    if m.sum() == 0:
        print(f"   attenuation {cname}: n=0")
        continue
    att = 1.0 / (1.0 + pj[m] * B2[m])
    print(f"   attenuation {cname:<5} n={int(m.sum()):<9} {fmtq(q(att))}   <<< Q6")
    out.setdefault("p1e_atten", {})[cname] = {"n": int(m.sum()), **q(att)}

# ═════════════════════════════════════════════════════════════════════ Q7, mechanical half
# The CLEAN isolation of the hypothesis. P1d damps ONLY the graft edges' RNA precisions and touches no
# DENSITY, so on the captured state ``S``, ``alpha`` and ``delta`` are IDENTICAL with and without it —
# only ``_vc = 1/p`` moves, hence ``aSa``, hence ``den``, hence ``b2``. Undo the damping algebraically
# (p = p'/(1 - p'*omega), the exact inverse of p' = p/(1+p*omega)) and recompute P1e on the SAME
# messages. No re-solve, no belief-change confound.
if not a.no_p1d:
    om = {0: float(glv[0]["omega"]), 1: float(glv[1]["omega"])}
    om_p, om_n = om[0], om[1]
    b2u_l, precu_l = [], []
    for _k, pin in enumerate(pins):
        gm = pin["graft"]
        tpp_u = np.where(gm & (pin["tpp"] > 0) & (pin["tpp"] * om_p < 1.0),
                         pin["tpp"] / np.maximum(1.0 - pin["tpp"] * om_p, 1e-300), pin["tpp"])
        tpn_u = np.where(gm & (pin["tpn"] > 0) & (pin["tpn"] * om_n < 1.0),
                         pin["tpn"] / np.maximum(1.0 - pin["tpn"] * om_n, 1e-300), pin["tpn"])
        pin_u = dict(pin, tpp=tpp_u, tpn=tpn_u)
        r_u = p1e_replay(pin_u, us, n_node)
        b2u_l.append(r_u["b2"])
        precu_l.append(r_u["prec"])
    B2U = np.concatenate(b2u_l)
    PRECU = np.concatenate(precu_l, axis=0)
    fireu = B2U > 0.0
    remu = np.where(fireu[:, None], PRECU - PRECU / (1.0 + PRECU * B2U[:, None]), 0.0).sum(axis=1)
    print("\nQ7-MECH  counterfactual: SAME captured messages, P1d's precision damping undone")
    print("   delta unchanged by construction (P1d touches no density) — verified: "
          "S identical, so only aSa/den/b2 move")
    for lab, m in (("ALL", np.ones_like(fire)), ("graft edges", GRAFT),
                   ("non-graft", ~GRAFT)):
        print(f"   {lab:<12} fire {int((fire & m).sum()):>7} -> {int((fireu & m).sum()):>7}   "
              f"removed {REM[m].sum():.4g} -> {remu[m].sum():.4g}   "
              f"mean b2(fire) {float(B2[fire & m].mean()) if (fire & m).any() else 0:.4g} -> "
              f"{float(B2U[fireu & m].mean()) if (fireu & m).any() else 0:.4g}")
    out["q7_mech"] = {
        "n_fire": int(fire.sum()), "n_fire_noP1d": int(fireu.sum()),
        "removed": float(REM.sum()), "removed_noP1d": float(remu.sum()),
        "n_fire_graft": int((fire & GRAFT).sum()),
        "n_fire_graft_noP1d": int((fireu & GRAFT).sum()),
        "removed_graft": float(REM[GRAFT].sum()),
        "removed_graft_noP1d": float(remu[GRAFT].sum()),
    }

if a.json_out:
    Path(a.json_out).write_text(json.dumps(out, indent=1, default=float))
    print(f"\n[json -> {a.json_out}]")
