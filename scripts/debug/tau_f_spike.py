"""Phase 0 — tau_f spike characterization (does the (M/E)^2 jacobian blow up at short-flank boundaries? P4).

Calls the CURRENT bp_solver._message directly, sweeping the DESTINATION effective length E_dst down toward a
tiny-exon boundary (E[min(l,R)] -> R ~ 10-50bp) while a full fragment's mass crosses the seam. Reports the
emitted message precision tau_f and whether it spikes past the destination's own one-fragment binomial
resolution (the count-space cap the plan's GATE 1c will enforce). Establishes the P4 baseline + whether the
binomial var_floor saves it at the simplex walls (it should NOT — f(1-f)->0 there).

    OMP_NUM_THREADS=1 python scripts/debug/tau_f_spike.py
"""
from __future__ import annotations
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np

from rigel.calibration.bp_solver import _message


class _StubVarMean:
    """sigma^2_bio(mu) = const (isolate the jacobian, not the var~mean)."""
    def __init__(self, s2=0.01):
        self.s2 = s2
    def predict(self, mu):
        return np.full(np.asarray(mu, float).shape, self.s2)


def run():
    vm = _StubVarMean(s2=0.01)
    # a confident source: density rho_src with plenty of source mass
    rho_src = 2.0          # source gDNA density (frag/bp-ish)
    eff_src = 250.0        # a normal-length source
    mass_src = 500.0       # well-sampled source
    rho_dst_cur = 2.0
    M_dst = 30.0           # a full fragment's worth of crossing mass at the dst seam

    print("=" * 92)
    print("TAU_F SPIKE — _message into a shrinking destination flank (M_dst=30, source confident)")
    print("=" * 92)
    print(f"  {'E_dst':>7} {'mu_f':>7} {'(M/E)^2':>12} {'tau_f':>14} {'1fragRes=1/(f(1-f)/M)':>22} {'> cap?':>7}")
    for eff_dst in (250.0, 100.0, 50.0, 25.0, 10.0, 5.0):
        mu_f, tau_f = _message(rho_src, eff_src, eff_dst, M_dst, mass_src, rho_dst_cur, vm)
        jac = (M_dst / max(eff_dst, 1e-9)) ** 2
        # the destination's honest one-fragment binomial precision at this mu_f
        one_frag = 1.0 / max(mu_f * (1 - mu_f) / max(M_dst, 1e-9), 1.0 / M_dst**2)
        flag = "SPIKE" if tau_f > 2 * one_frag else ""
        print(f"  {eff_dst:>7.1f} {mu_f:>7.3f} {jac:>12.1f} {tau_f:>14.3g} {one_frag:>22.3g} {flag:>7}")

    print("\n  --- at the SIMPLEX WALL (mu_f -> 1: a near-pure-gDNA tiny-exon boundary) ---")
    print("  (binomial var_floor = f(1-f)(1/Ms+1/Md) -> 0 as mu_f->1, so it canNOT bound tau_f here)")
    print(f"  {'E_dst':>7} {'rho_src':>8} {'mu_f':>7} {'tau_f':>14}")
    for rho_hi, eff_dst in ((8.0, 25.0), (8.0, 10.0), (20.0, 10.0), (20.0, 5.0)):
        # high source density vs short flank pushes mu_f toward 1 (the wall)
        mu_f, tau_f = _message(rho_hi, eff_src, eff_dst, M_dst, mass_src, rho_dst_cur, vm)
        print(f"  {eff_dst:>7.1f} {rho_hi:>8.1f} {mu_f:>7.3f} {tau_f:>14.3g}")


if __name__ == "__main__":
    run()
