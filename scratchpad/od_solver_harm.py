"""Direct measurement of the harm at the CEILING: what does the psi solve return for a
pure-gDNA node (exact 50/50) and a pure-RNA control, as od goes to the ceiling?

Mirrors tests/calibration/test_bp_solver.py::test_pure_gdna_node_confident_at_near_binomial_od
but sweeps od up to the proposal's value (0.2 = Beta(2,2) ceiling).
"""
from __future__ import annotations

import numpy as np

from rigel.calibration.simplex_logodds import _solve_nodes_logodds_all


def solve(u_pos, u_neg, od, kappa):
    z = np.zeros(1)
    n = float(u_pos + u_neg)
    return float(
        _solve_nodes_logodds_all(
            np.array([float(u_pos)]),
            np.array([float(u_neg)]),
            np.array([True]),
            np.array([False]),
            np.array([n]),
            z,
            kappa=kappa,
            od_g=od,
            od_r=od,
            n_grid=80,
        ).gdna_frac[0]
    )


ODS = [0.0345, 0.10, 0.1429, 0.20]
LBL = ["a=14 (today)", "a=4.5 (C/2)", "a=3 (rejected)", "a=2 (CEILING = the proposal)"]

for kappa in (0.7, 0.9, 0.99):
    print(f"\n=== kappa = {kappa} ===")
    for tag, (up, un) in [
        ("pure gDNA N=1000 (500/500), truth f_g=1", (500, 500)),
        ("pure gDNA N=100  (50/50),   truth f_g=1", (50, 50)),
        ("pure gDNA N=20   (10/10),   truth f_g=1", (10, 10)),
        (f"pure RNA  N=1000 (kappa),   truth f_g=0", (int(1000 * kappa), 1000 - int(1000 * kappa))),
        (f"pure RNA  N=100  (kappa),   truth f_g=0", (int(100 * kappa), 100 - int(100 * kappa))),
    ]:
        vals = [solve(up, un, od, kappa) for od in ODS]
        print(f"  {tag:42s} " + "  ".join(f"{l.split()[0]}:{v:.4f}" for l, v in zip(LBL, vals)))

print("\nlegend: a=14 -> od=0.0345 (shipped) | a=4.5 -> 0.10 (C/2) | a=3 -> 0.1429 (measured-bad,"
      " reverted in d6f79c53) | a=2 -> 0.20 (the ceiling = the proposal)")
