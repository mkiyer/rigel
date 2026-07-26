"""RC3 — the 'STRUCTURAL ABSENCE DELIVERED AS A CONFIDENT CLAIM' bug-class scan.

Enumerates every place in `_relay` / `_transport` / the combine where a quantity can be structurally ZERO
(or saturated at a wall) and is nonetheless handed to ψ as a MODE with LIVE precision, then measures each
one's fire-rate + the destination error mass, over the whole suite.

    OMP_NUM_THREADS=1 python scratchpad/rc3_hunt.py [--conds a,b] [--top 0]
"""

from __future__ import annotations

import argparse
import sys

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
from rc3_replay import CLS, SUITE, build, load, oracle, validate  # noqa: E402

_EPS = 1e-9
LOGEPS = np.log(_EPS)

ap = argparse.ArgumentParser()
ap.add_argument("--conds", default=None)
ap.add_argument("--out", default="/tmp/rc3_sites.npz")
a = ap.parse_args()

conds = (
    a.conds.split(",")
    if a.conds
    else sorted(d.name for d in SUITE.iterdir() if (d / "sim_oracle.bam").exists())
)

# per-site accumulators: name -> [n_fired, err_mass_fired, n_fired_fullrank, err_mass_fullrank, node_mass]
acc: dict[str, np.ndarray] = {}
tot = np.zeros(4)  # nodes, err_mass, fullrank nodes, fullrank err_mass
hurt_tot = np.zeros(2)
extra: dict[str, list] = {}


def hit(name, m, err_mass, full, hurt, mass):
    v = acc.setdefault(name, np.zeros(7))
    m = np.asarray(m, bool)
    v[0] += m.sum()
    v[1] += err_mass[m].sum()
    v[2] += (m & full).sum()
    v[3] += err_mass[m & full].sum()
    v[4] += mass[m].sum()
    v[5] += (m & hurt).sum()
    v[6] += err_mass[m & hurt].sum()


for ci, cond in enumerate(conds):
    ctx = load(cond)
    S = build(ctx)
    fid = max(validate(ctx, S).values())
    assert fid == 0.0, f"replay drift {fid}"
    fo = oracle(ctx)
    mass, solved, self_fg = S["mass"], S["solved"], S["fg_loc"]
    ok = np.isfinite(fo) & (mass > _EPS)
    err = np.where(ok, np.abs(solved - fo), 0.0)
    serr = np.where(ok, np.abs(self_fg - fo), 0.0)
    em = err * mass
    full = (S["tau_own"] > _EPS) & ok
    hurt = full & (err > serr + 0.02)
    tot += [ok.sum(), em[ok].sum(), full.sum(), em[full].sum()]
    hurt_tot += [hurt.sum(), em[hurt].sum()]

    A, B = S["A"], S["B"]
    cg, cp, cn, cR = S["cg"], S["cp"], S["cn"], S["cR"]
    cm_g, cm_p, cm_n, c_tau = S["cm_g"], S["cm_p"], S["cm_n"], S["c_tau"]
    v_own_r, v_own_g, v_own_lam = S["v_own_r"], S["v_own_g"], S["v_own_lam"]

    # ── S1  THE PEEL: a fully-consumed peel emits "no RNA continues past here" at live precision ────────
    for tag, D in (("L", A), ("R", B)):
        for s, pre, post, mp, mm in (
            ("p", "tp_pre_peel", "tp_post_peel", "tpp", "tmp"),
            ("n", "tn_pre_peel", "tn_post_peel", "tpn", "tmn"),
        ):
            fired = D["peel"] & (D[pre] > _EPS) & (D[post] <= _EPS)
            live = fired & ((D[mp] > 0.0) | (D[mm] > 0.0))
            hit(f"S1 peel-zero emitted @live prec ({s})", live, em, full, hurt, mass)
            # partial peel: the M3 u-weight (derived, unwired)
            part = D["peel"] & (D[pre] > _EPS) & (D[post] > _EPS)
            u = np.where(part, D[pre] / np.maximum(D[post], _EPS), np.nan)
            extra.setdefault("peel_u", []).append(u[part])
            hit(f"S1b peel u>3 (M3 invalid) ({s})", part & (u > 3.0), em, full, hurt, mass)

    # ── S2  MODE FLOORS: a component whose fused density is 0 but whose precision is live ───────────────
    hit("S2 mo_g at log(_EPS) & cm_g>0", (cg <= _EPS) & (cm_g > 0), em, full, hurt, mass)
    hit("S2 mo_p at log(_EPS) & cm_p>0", (cp <= _EPS) & (cm_p > 0), em, full, hurt, mass)
    hit("S2 mo_n at log(_EPS) & cm_n>0", (cn <= _EPS) & (cm_n > 0), em, full, hurt, mass)
    hit("S2 mo_R at log(_EPS) & c_tau>0", (cR <= _EPS) & (c_tau > 0), em, full, hurt, mass)
    # the same on a LIVE strand only (a dead strand's ψ message is a no-op through f_act)
    hit("S2* mo_p floor, +strand LIVE", (cp <= _EPS) & (cm_p > 0) & S["fp_a"], em, full, hurt, mass)
    hit("S2* mo_n floor, −strand LIVE", (cn <= _EPS) & (cm_n > 0) & S["fn_a"], em, full, hurt, mass)

    # ── S3  DL's contradiction mask is gated on `known` — it does nothing where v_own = ∞ ───────────────
    for tag, D in (("L", A), ("R", B)):
        for s, c, vo, mm in (("g", "contra_g", v_own_g, "tmg_pre_dl"),
                             ("p", "contra_p", v_own_r, "tmp_pre_dl"),
                             ("n", "contra_n", v_own_r, "tmn_pre_dl")):
            m = D[c] & ~np.isfinite(vo) & (D[mm] > 0.0)
            hit(f"S3 contradicted msg SURVIVES (v_own=inf) ({s})", m, em, full, hurt, mass)

    # ── S4  the θ tilt message: saturated at the wall / built from nothing ──────────────────────────────
    amb, thp = S["is_amb"], S["th_prec"]
    hit("S4 θ msg pinned at |τ|=1 (one arm 0)", amb & (thp > 0) & (np.abs(S["tau_tilt"]) >= 1.0 - 1e-12)
        & ((cp <= _EPS) | (cn <= _EPS)), em, full, hurt, mass)
    hit("S4b θ msg = 0 from cR=0", amb & (thp > 0) & (cR <= _EPS), em, full, hurt, mass)

    # ── S5  _pin_v substitutes the DESTINATION's own density for an uninformed component ────────────────
    for tag, D in (("L", A), ("R", B)):
        # substituting the UNSOLVED-DEFAULT own gDNA (f_g=1 at ZERO precision) as if it were a real level
        unsolved_default = (S["tau_own"] <= _EPS) & ~S["struct"] & (S["og"] > _EPS)
        m = D["pin_sub_g"] & unsolved_default & ((D["tpp"] > 0) | (D["tpn"] > 0))
        hit("S5 pin uses unsolved-default own gDNA", m, em, full, hurt, mass)

    # ── S6  the graft's MEASUREMENT precision is the spliced MASS, not the spliced COUNT ────────────────
    for tag, D in (("L", A), ("R", B)):
        g = D["sp_mass"] > _EPS
        extra.setdefault("spl_mass_over_count", []).append(
            (D["sp_mass"][g] / np.maximum(D["sp_count"][g], _EPS))
        )
        extra.setdefault("spl_prec", []).append(D["spc"][g])

    # ── S7  the λ-emission gate: how often does it fire (the LANDED template fix), split by CAUSE ───────
    for tag, D in (("L", A), ("R", B)):
        k = D["lam_gate_killed"]
        hit("S7 λ-gate fired: source RNA-free", k & ((D["tp"] + D["tn"]) <= _EPS), em, full, hurt, mass)
        hit("S7 λ-gate fired: source gDNA-free", k & (D["tg"] <= _EPS), em, full, hurt, mass)
    hit("S7b combine λ-gate fired", (S["c_tau_pre"] > 0) & (c_tau <= 0), em, full, hurt, mass)

    # ── S8  _pin_v pushes a PARTIAL message onto the missing component's vertex ─────────────────────────
    # A component the message does not carry AND the destination has no own density for is supplied by
    # NOBODY, so `s` omits it and k renormalizes the surviving component to the f_c = 1 vertex.
    og_, op_, on_ = S["og"], S["op"], S["on"]
    for tag, D in (("L", A), ("R", B)):
        vtx = (D["pin_sub_p"] & (op_ <= _EPS)) & (D["pin_sub_n"] & (on_ <= _EPS)) & (D["tpg"] > 0)
        hit("S8 pin → f_g=1 vertex (RNA absent both)", vtx, em, full, hurt, mass)
        vtx2 = D["pin_sub_g"] & (og_ <= _EPS) & ((D["tpp"] > 0) | (D["tpn"] > 0))
        hit("S8b pin → f_R=1 vertex (gDNA absent)", vtx2, em, full, hurt, mass)
        # the SATURATION form: the raw (pre-pin) message claims more than the node's whole mass
        raw = (D["tg_pre_pin"] * S["E_g"] + (D["tp_pre_pin"] + D["tn_pre_pin"]) * S["E_r"]) \
            / np.maximum(S["M"], _EPS)
        extra.setdefault("pin_raw_sum_fc", []).append(raw[(D["tpg"] > 0) | (D["tpp"] > 0) | (D["tpn"] > 0)])
        hit("S8c pin scale |log k| > 1 (saturated msg)", np.abs(np.log(np.maximum(D["k_pin"], 1e-30))) > 1.0,
            em, full, hurt, mass)

    # ── S9  a MATURE-derived measurement precision delivered where mature structurally cannot exist ─────
    st = ctx["dbg"]["statics"]
    mrp = np.asarray(st.mrna_active_pos, bool)
    mrn = np.asarray(st.mrna_active_neg, bool)
    hit("S9 mature meas prec at non-mature node (+)", (cm_p > 0) & ~mrp & S["fp_a"], em, full, hurt, mass)
    hit("S9 mature meas prec at non-mature node (−)", (cm_n > 0) & ~mrn & S["fn_a"], em, full, hurt, mass)

    # ── S10  MISSING FRAME: one side has no ρ_tot ⇒ r is UNDEFINED and the code passes r = 1.0 (identity),
    # i.e. the source's density is delivered in the SOURCE's own enrichment frame, at full precision.
    for tag, D, val in (("L", A, S["vl"]), ("R", B, S["vr"])):
        nof = val & (D["r"] == 1.0)
        anyp = (D["tpg"] > 0) | (D["tpp"] > 0) | (D["tpn"] > 0) | (D["tmg"] > 0) | (D["ttau"] > 0)
        hit("S10 r=1 identity frame (undefined) @live", nof & anyp, em, full, hurt, mass)

    print(f"  [{ci + 1:>2}/{len(conds)}] {cond}", flush=True)

print(f"\nnodes={tot[0]:,.0f}  err-mass={tot[1]:,.0f}   full-rank nodes={tot[2]:,.0f} "
      f"err-mass={tot[3]:,.0f}   HURT nodes={hurt_tot[0]:,.0f} err-mass={hurt_tot[1]:,.0f}")
print(f"\n{'site':<44}{'fired':>9}{'err-mass':>12}{'%all':>7}"
      f"{'fullrk':>8}{'fr-em':>11}{'hurt':>7}{'hurt-em':>10}{'%hurt':>7}")
for k in sorted(acc):
    v = acc[k]
    print(f"{k:<44}{v[0]:>9,.0f}{v[1]:>12,.0f}{v[1] / tot[1]:>6.1%}"
          f"{v[2]:>8,.0f}{v[3]:>11,.0f}{v[5]:>7,.0f}{v[6]:>10,.0f}{v[6] / hurt_tot[1]:>6.1%}")

for k, lst in extra.items():
    x = np.concatenate([np.asarray(z) for z in lst])
    x = x[np.isfinite(x)]
    if x.size:
        q = np.percentile(x, [50, 75, 90, 99])
        print(f"\n{k}: n={x.size:,}  median={q[0]:.4g} p75={q[1]:.4g} p90={q[2]:.4g} p99={q[3]:.4g} "
              f"max={x.max():.4g}")
