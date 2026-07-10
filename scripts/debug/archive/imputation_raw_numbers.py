"""The imputation, laid bare — the FAITHFUL pass-by-pass trajectory of the boundary node at 5000.

The boundary holds 191 clean +RNA crossing fragments (u_pos=2, u_neg=191) → strand alone says f_g≈0. Yet it
ends phantom'd. Single-message reconstruction is unreliable (Gauss-Seidel updates neighbours mid-pass), so this
records, at EVERY _solve of that node, the actual messages it receives (gDNA μ_g/N_g, RNA+ μ_p/N_p) and the
f_g it is assigned — so we watch exactly when and by which message the phantom enters.

    OMP_NUM_THREADS=1 python scripts/debug/imputation_raw_numbers.py
"""
from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
import tempfile
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parent))
from ambig_nascent_trace import build  # noqa: E402

import rigel.calibration.bp_solver as bp  # noqa: E402
from rigel.calibration.simplex_sweep import _fg_median, _simplex_lattice  # noqa: E402


def main():
    with tempfile.TemporaryDirectory() as work:
        sc, pl, ra, sm, fl, cfg = build(work)
        lat = _simplex_lattice(int(cfg.calibration.sweep_n_grid))
        fgg = lat[2]

        traj = []
        state = {"pass": 0}
        orig_ll = bp._local_loglik
        orig_nd = bp.node_densities

        def nd_spy(*a, **k):
            state["pass"] += 1
            return orig_nd(*a, **k)

        def ll_spy(u_pos, u_neg, *a, **k):
            psi = orig_ll(u_pos, u_neg, *a, **k)
            # boundary 5000: ~2 pos / ~191 neg crossing, single-strand + (allow_pos, not allow_neg)
            if u_pos[0] < 5 and 185 <= u_neg[0] <= 196 and bool(a[2][0]) and not bool(a[3][0]):
                gf = k.get("gdna_imp_frac")
                gc = k.get("gdna_imp_count")
                rf = k.get("rna_imp_frac")
                rc = k.get("rna_imp_count")
                traj.append((
                    state["pass"],
                    float(gf[0]) if gf is not None else 0.0,
                    float(gc[0]) if gc is not None else 0.0,
                    float(rf[0][0]) if rf is not None else 0.0,
                    float(rc[0][0]) if rc is not None else 0.0,
                    float(_fg_median(psi, fgg)[0]),
                ))
            return psi

        bp._local_loglik = ll_spy
        bp.node_densities = nd_spy
        from rigel.calibration.calibrate import calibrate as _cal
        _cal(pl, ra, sm, fl.gdna_pmf, fl.rna_pmf, cfg.calibration)
        bp._local_loglik = orig_ll
        bp.node_densities = orig_nd

        print("\n=== boundary 5000 (u_pos=2, u_neg=191 → strand alone: f_g≈0) — every solve, in order ===")
        print("  Each row = one _solve of this node. μ_g/N_g = gDNA message; μ_p/N_p = +RNA message (the")
        print("  imputed FRACTION and its pseudo-count). f_g = what the node is assigned after that message.")
        print(f"\n  {'pass':>4} {'μ_g':>7} {'N_g':>7} {'μ_p':>7} {'N_p':>7} {'→ f_g':>8}   note")
        prev = 0.0
        for p, mg, ng, mp, npc, fg in traj:
            note = ""
            if fg - prev > 0.05:
                note = f"PHANTOM injected (+{fg - prev:.2f})"
            elif fg < 0.05:
                note = "clean (~0)"
            # which message would push f_g up? a gDNA μ_g>0, or an RNA μ_p<1 (low + frac ⇒ high f_g)
            drivers = []
            if mg > 0.05 and ng > 0.5:
                drivers.append(f"gDNA μ_g={mg:.2f}@N={ng:.1f}")
            if mp < 0.9 and npc > 0.5:
                drivers.append(f"RNA μ_p={mp:.2f}@N={npc:.1f} (low + ⇒ ↑f_g)")
            print(f"  {p:>4} {mg:>7.3f} {ng:>7.2f} {mp:>7.3f} {npc:>7.2f} {fg:>8.3f}   "
                  f"{note}{'  ['+'; '.join(drivers)+']' if drivers else ''}")
            prev = fg
        sc.cleanup()


if __name__ == "__main__":
    main()
