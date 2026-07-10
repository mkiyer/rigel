"""Phase 0 — isolate & lock the precision-state mechanism (NO production change).

Validates the Q2 fix in isolation, on the real AMBIG-locus geometry (gDNA=0, nascent, ss0.99):

  1. Var(f_g) from the lattice: balanced AMBIG  ≫  confident single-strand (the strand-resolving-power claim).
  2. Var_own = (M/E)²·Var(f_g): AMBIG (source) ≫ boundary (dest).
  3. message strength  N_OLD = ρ²/(σ²_bio+pois)   vs   N_NEW = ρ²/(Var_own + σ²_bio + pois)  → N_NEW ≪ N_OLD.
  4. apply the AMBIG→boundary gDNA message to the boundary, at N_OLD and N_NEW, via BOTH forms
     (count-binom, density-Gaussian) → does the boundary stay f_g≈0? (the form decision, §7).

    OMP_NUM_THREADS=1 python scripts/debug/phase0_precision_state_probe.py
"""
from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ambig_nascent_trace import build  # noqa: E402

import rigel.calibration.bp_solver as bp  # noqa: E402
from rigel.calibration.calibrate import calibrate  # noqa: E402
from rigel.calibration.node_chain import REGION  # noqa: E402
from rigel.calibration.simplex_sweep import (  # noqa: E402
    _binom_pseudo,
    _fg_median,
    _fg_var,
    _local_loglik,
    _simplex_lattice,
)

_EPS = 1e-9
_Z = np.array([0.0])
_T = np.array([True])
_F = np.array([False])
_MSG_PSEUDO = 1.0


def lattice_solve(u_pos, u_neg, kappa, odg, odr, *, single_strand, msg=None):
    """Solve one node; optional gDNA message msg=(form, mu, N, M, E). Returns (f_g median, Var(f_g))."""
    lat = _simplex_lattice(60)
    fgg = lat[2]
    allow_p, allow_n = _T, (_F if single_strand else _T)
    so = _T if single_strand else None
    psi = _local_loglik(np.array([float(u_pos)]), np.array([float(u_neg)]), _Z, _Z, allow_p, allow_n,
                        kappa, odg, odr, lat, strand_obs=so)
    if msg is not None:
        form, mu, N, M, E = msg
        if form == "binom":
            mu_c = float(np.clip(mu, 1e-6, 1 - 1e-6))
            psi = psi + _binom_pseudo(fgg[None, :], np.array([[mu_c]]), np.array([[N]]))
        elif form == "gauss":
            rho_msg = mu * M / max(E, _EPS)            # the source density the fraction mu corresponds to
            tau_density = N / max(rho_msg * rho_msg, _EPS)  # N effective counts ⇒ density precision
            dens = fgg * M / max(E, _EPS)
            psi = psi + (-0.5 * tau_density * (dens - rho_msg) ** 2)[None, :]
    return float(_fg_median(psi, fgg)[0]), float(_fg_var(psi, fgg)[0])


def main():
    with tempfile.TemporaryDirectory() as work:
        sc, pl, ra, sm, fl, cfg = build(work)
        cap = {}
        orig_sweep = bp.node_sweep
        orig_seed = bp._gdna_seed_estimate

        def sweep_spy(chain, statics, geometry, belief, *a, **k):
            cap.update(chain=chain, geom=geometry, statics=statics)
            return orig_sweep(chain, statics, geometry, belief, *a, **k)

        def seed_spy(*a, **k):
            out = orig_seed(*a, **k)
            cap["gdna_vm"] = out[1]
            return out

        bp.node_sweep = sweep_spy
        bp._gdna_seed_estimate = seed_spy
        calibrate.__globals__["node_sweep"] = sweep_spy
        res = calibrate(pl, ra, sm, fl.gdna_pmf, fl.rna_pmf, cfg.calibration)
        bp.node_sweep = orig_sweep
        bp._gdna_seed_estimate = orig_seed
        calibrate.__globals__["node_sweep"] = orig_sweep

        kappa = res.rna_sense_frac
        odg = res.gdna_strand_overdispersion
        odr = res.rna_strand_overdispersion
        chain, geom, st = cap["chain"], cap["geom"], cap["statics"]
        gdna_vm = cap["gdna_vm"]
        kind = np.asarray(chain.kind)
        ref_idx = np.asarray(chain.ref_idx)
        is_reg = kind == REGION
        start = np.asarray(ra.start)
        amb_r = int(np.where(start == 5000)[0][0])
        amb = int(np.where(is_reg & (ref_idx == amb_r))[0][0])
        bnd = int(np.asarray(chain.left)[amb])  # boundary 5000 = AMBIG's left neighbour

        def geomvals(node, face=0):
            M = (geom.mass_left if face == 0 else geom.mass_right)[node]
            Eg = (geom.eff_gdna_left if face == 0 else geom.eff_gdna_right)[node]
            return float(M), float(Eg)

        M_amb, E_amb = geomvals(amb)
        M_bnd, E_bnd = geomvals(bnd)
        up_amb, un_amb = float(st.u_pos[amb]), float(st.u_neg[amb])
        up_bnd, un_bnd = float(st.u_pos[bnd]), float(st.u_neg[bnd])

        print(f"\n=== Phase 0: precision-state mechanism (κ={kappa:.4f}, od_g={odg:.3f}, od_r={odr:.3f}) ===")
        print(f"  AMBIG node  : u_pos={up_amb:.0f} u_neg={un_amb:.0f}  M={M_amb:.0f} E_gdna={E_amb:.0f}")
        print(f"  boundary    : u_pos={up_bnd:.0f} u_neg={un_bnd:.0f}  M={M_bnd:.0f} E_gdna={E_bnd:.0f}")

        # (1) Var(f_g) from the lattice
        fg_amb, var_amb = lattice_solve(up_amb, un_amb, kappa, odg, odr, single_strand=False)
        fg_bnd, var_bnd = lattice_solve(up_bnd, un_bnd, kappa, odg, odr, single_strand=True)
        print("\n  [1] Var(f_g) from the lattice (strand-resolving power):")
        print(f"      AMBIG (balanced)        f_g={fg_amb:.3f}  Var(f_g)={var_amb:.5f}")
        print(f"      boundary (191:2)        f_g={fg_bnd:.3f}  Var(f_g)={var_bnd:.5f}")
        print(f"      ratio AMBIG/boundary = {var_amb/max(var_bnd,1e-9):.0f}×   "
              f"{'PASS' if var_amb > 10*var_bnd else 'CHECK'} (expect ≫1)")

        # (2) Var_own = (M/E)²·Var(f_g)
        var_own_amb = (M_amb / max(E_amb, _EPS)) ** 2 * var_amb
        var_own_bnd = (M_bnd / max(E_bnd, _EPS)) ** 2 * var_bnd
        print("\n  [2] Var_own(ρ_g) = (M/E)²·Var(f_g):")
        print(f"      AMBIG    = ({M_amb:.0f}/{E_amb:.0f})²·{var_amb:.4f} = {var_own_amb:.4f}")
        print(f"      boundary = ({M_bnd:.0f}/{E_bnd:.0f})²·{var_bnd:.5f} = {var_own_bnd:.5f}")

        # (3) message strength AMBIG→boundary (use the AMBIG's strand-solved f_g for ρ_src)
        rho_src = fg_amb * M_amb / max(E_amb, _EPS)
        s2_bio = float(gdna_vm.predict(np.array([max(rho_src, _EPS)]))[0])
        pois = (rho_src + _MSG_PSEUDO / max(E_amb, _EPS)) / max(E_amb, _EPS)
        N_old = rho_src ** 2 / max(s2_bio + pois, _EPS)
        N_new = rho_src ** 2 / max(var_own_amb + s2_bio + pois, _EPS)
        mu_c = float(np.clip(rho_src * E_bnd / max(M_bnd, _EPS), 0.0, 1.0))
        print("\n  [3] message strength AMBIG→boundary:")
        print(f"      ρ_src={rho_src:.3f}  σ²_bio={s2_bio:.5f}  pois={pois:.6f}  →  μ_c(at boundary)={mu_c:.3f}")
        print(f"      N_OLD = ρ²/(σ²_bio+pois)            = {N_old:.2f}   (today — the corruptor)")
        print(f"      N_NEW = ρ²/(Var_own+σ²_bio+pois)    = {N_new:.2f}   "
              f"{'PASS' if N_new < 0.2 * N_old else 'CHECK'} (expect ≪ N_OLD)")

        # (4) apply to the boundary, both N and both forms → does it stay f_g≈0?
        print("\n  [4] boundary f_g after the gDNA message (boundary alone solves to f_g≈0):")
        print(f"      {'config':>22} {'count-binom':>12} {'density-Gauss':>13}")
        for tag, N in (("N_OLD (today)", N_old), ("N_NEW (Var_own fix)", N_new)):
            fb = lattice_solve(up_bnd, un_bnd, kappa, odg, odr, single_strand=True,
                               msg=("binom", mu_c, N, M_bnd, E_bnd))[0]
            fgs = lattice_solve(up_bnd, un_bnd, kappa, odg, odr, single_strand=True,
                                msg=("gauss", mu_c, N, M_bnd, E_bnd))[0]
            print(f"      {tag:>22} {fb:>12.3f} {fgs:>13.3f}")
        print("\n  DECISION: with N_NEW, does count-binom keep f_g≈0? If yes → keep count-binom (discreteness).")
        print("            If the residual wall keeps f_g high at N_NEW → density-Gaussian needed.")
        sc.cleanup()


if __name__ == "__main__":
    main()
