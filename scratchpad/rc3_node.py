"""RC3 — full per-message provenance for one node, straight off the bit-faithful combine replay."""

from __future__ import annotations

import argparse
import sys

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
from rc3_replay import CLS, build, load, oracle, validate  # noqa: E402

_EPS = 1e-9
ap = argparse.ArgumentParser()
ap.add_argument("cond")
ap.add_argument("nodes", nargs="+", type=int)
a = ap.parse_args()

ctx = load(a.cond)
S = build(ctx)
assert max(validate(ctx, S).values()) == 0.0
fo = oracle(ctx)
us = S["us"]

for i in a.nodes:
    cl = CLS[int(S["cls"][i])]
    dof = "AMBIG" if S["is_amb"][i] else ("single" if S["fp_a"][i] ^ S["fn_a"][i] else "locked")
    print(f"\n=== node {i} [{cl}/{dof}] cond={a.cond}")
    print(f"  mass={S['mass'][i]:,.1f}  M={S['M'][i]:,.1f} E_g={S['E_g'][i]:.1f} E_r={S['E_r'][i]:.1f}")
    print(f"  oracle={fo[i]:.4f}  self-solve={S['fg_loc'][i]:.4f}  SOLVED={S['solved'][i]:.4f}")
    print(f"  tau_own={S['tau_own'][i]:.4g}  struct={bool(S['struct'][i])}  "
          f"v_own_g={S['v_own_g'][i]:.4g} v_own_r={S['v_own_r'][i]:.4g} v_own_lam={S['v_own_lam'][i]:.4g}")
    print(f"  own densities  og={S['og'][i]:.5g} op={S['op'][i]:.5g} on={S['on'][i]:.5g}   "
          f"own f_g·M/E_g check={S['fg_loc'][i] * S['M'][i] / S['E_g'][i]:.5g}")
    print(f"  ψ RECEIVES: mo_g={np.exp(S['mo_g'][i]):.4f}/{S['cm_g'][i]:.4g}  "
          f"mo_p={np.exp(S['mo_p'][i]):.4f}/{S['cm_p'][i]:.4g}  mo_n={np.exp(S['mo_n'][i]):.4f}/{S['cm_n'][i]:.4g}")
    print(f"              λ f_g_eq={1 / (1 + np.exp(-S['lam_msg'][i])):.4f}/{S['c_tau'][i]:.4g} "
          f"(pre-gate {S['c_tau_pre'][i]:.4g})   θ={S['th_msg'][i]:+.3f}/{S['th_prec'][i]:.4g} "
          f"τ_tilt={S['tau_tilt'][i]:+.3f}")
    for tag, D, srcarr in (("LEFT ", S["A"], S["sl"]), ("RIGHT", S["B"], S["sr"])):
        s = int(srcarr[i])
        valid = (S["li"] if tag == "LEFT " else S["ri"])[i] >= 0
        if not valid:
            print(f"  {tag} msg: (no neighbour)")
            continue
        scl = CLS[int(S["cls"][s])]
        print(f"  {tag} msg  src={s} [{scl}] src_oracle={fo[s]:.4f} src_solved={S['solved'][s]:.4f} "
              f"src_self={S['fg_loc'][s]:.4f} src_tau={S['tau_own'][s]:.4g}")
        print(f"        r={D['r'][i]:.4g} s2t={D['s2t'][i]:.4g} graft={bool(D['graft'][i])} "
              f"peel={bool(D['peel'][i])}")
        print(f"        relayed src rho: g={us['fwd_g' if tag == 'LEFT ' else 'bwd_g'][s]:.5g} "
              f"p={us['fwd_p' if tag == 'LEFT ' else 'bwd_p'][s]:.5g} "
              f"n={us['fwd_n' if tag == 'LEFT ' else 'bwd_n'][s]:.5g}")
        print(f"        pre-peel tp={D['tp_pre_peel'][i]:.5g} tn={D['tn_pre_peel'][i]:.5g}  "
              f"mature_p={D['mature_p'][i]:.5g} mature_n={D['mature_n'][i]:.5g}  "
              f"post-peel tp={D['tp_post_peel'][i]:.5g} tn={D['tn_post_peel'][i]:.5g}")
        print(f"        pre-pin  tg={D['tg_pre_pin'][i]:.5g} tp={D['tp_pre_pin'][i]:.5g} "
              f"tn={D['tn_pre_pin'][i]:.5g}   k_pin={D['k_pin'][i]:.5g} "
              f"(sub g={bool(D['pin_sub_g'][i])} p={bool(D['pin_sub_p'][i])} n={bool(D['pin_sub_n'][i])})")
        fgm = D["tg"][i] * S["E_g"][i] / max(S["M"][i], _EPS)
        fpm = D["tp"][i] * S["E_r"][i] / max(S["M"][i], _EPS)
        fnm = D["tn"][i] * S["E_r"][i] / max(S["M"][i], _EPS)
        print(f"        POST-PIN implied fractions: f_g={fgm:.4f} f_p={fpm:.4f} f_n={fnm:.4f} "
              f"(Σ={fgm + fpm + fnm:.4f})")
        print(f"        graft spliced mass={D['sp_mass'][i]:.4g}/{D['sn_mass'][i]:.4g} "
              f"count={D['sp_count'][i]:.4g}/{D['sn_count'][i]:.4g} -> prec _spc={D['spc'][i]:.4g} "
              f"_snc={D['snc'][i]:.4g}")
        print(f"        prec: mode-fuse g={D['tpg'][i]:.4g} p={D['tpp'][i]:.4g} n={D['tpn'][i]:.4g}   "
              f"meas g={D['tmg'][i]:.4g} p={D['tmp'][i]:.4g} n={D['tmn'][i]:.4g}   τ={D['ttau'][i]:.4g}")
        print(f"        pre-DL meas: g={D['tmg_pre_dl'][i]:.4g} p={D['tmp_pre_dl'][i]:.4g} "
              f"n={D['tmn_pre_dl'][i]:.4g}   contra g={bool(D['contra_g'][i])} p={bool(D['contra_p'][i])} "
              f"n={bool(D['contra_n'][i])}   λ-gate killed={bool(D['lam_gate_killed'][i])}")
