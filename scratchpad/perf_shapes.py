"""Report the real per-solve shapes `_solve_nodes_logodds_all` sees, at both scales.

The 1-D and 2-D solvers each materialise several (m, K) / (m, K, K_t) temporaries; whether that is a
cache-resident kernel or a gigabyte of streamed memory is the whole question for cluster C, and it cannot
be read off the toy.
"""
from __future__ import annotations

import argparse
import os
import pickle
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")

import rigel.calibration.simplex_logodds as slo  # noqa: E402

CF = Path("/Users/mkiyer/Downloads/rigel_runs/cfrna/_calib_cache")
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
TOY_COND = "gdna_gdna300_ss_0.50_nrna_present_capture_off"

ap = argparse.ArgumentParser()
ap.add_argument("--input", default="toy")
ap.add_argument("--refit", type=int, default=1)
a = ap.parse_args()

_orig_ss, _orig_amb, _orig_rg = (
    slo._solve_nodes_logodds,
    slo._solve_ambig_logodds,
    slo._regrid_global,
)
log: list[str] = []


def _ss(u_pos, *rest, **kw):
    log.append(f"  1-D single-strand: m={np.asarray(u_pos).shape[0]:>8}  K={kw['n_grid']}"
               f"  -> (m,K) f64 = {np.asarray(u_pos).shape[0] * kw['n_grid'] * 8 / 2**20:8.1f} MB/temp")
    return _orig_ss(u_pos, *rest, **kw)


def _amb(u_pos, *rest, **kw):
    m, K = np.asarray(u_pos).shape[0], kw["n_grid"]
    Kt = kw.get("n_tilt") or K
    log.append(f"  2-D AMBIG batch  : m={m:>8}  K={K} Kt={Kt}"
               f"  -> (m,K,Kt) f32 = {m * K * Kt * 4 / 2**20:8.1f} MB/cube")
    return _orig_amb(u_pos, *rest, **kw)


def _rg(glp, n_from, n_to, L):
    if glp is not None and int(n_from) != int(n_to):
        g = np.asarray(glp)
        log.append(f"  _regrid_global   : {g.shape} -> ({g.shape[0]},{n_to})"
                   f"  = {g.shape[0] * n_to * 8 / 2**20:8.1f} MB/temp")
    return _orig_rg(glp, n_from, n_to, L)


slo._solve_nodes_logodds, slo._solve_ambig_logodds, slo._regrid_global = _ss, _amb, _rg

sys.argv = ["perf_scale", "--input", a.input, "--refit", str(a.refit), "--repeat", "0"]
exec(compile(Path("/Users/mkiyer/proj/rigel/scratchpad/perf_scale.py").read_text(), "perf_scale", "exec"))

seen: dict[str, int] = {}
for line in log:
    seen[line] = seen.get(line, 0) + 1
print(f"\n=== solve shapes ({a.input}, refit={a.refit}) — {len(log)} calls ===")
for line, c in seen.items():
    print(f"x{c:<3}{line}")
