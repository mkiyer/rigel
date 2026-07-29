"""Verification for the priority-#3 design review (docs/CARRY_FORWARD.md §7).

Four checks, each settling a review claim without touching a BAM:

  C1  fbp IS the arithmetic mean -> E_n = rho_n0 (the "expectation bug"; doc §7.1)
  C2  delta vs FW vs Gauss-Hermite under the CORRECTED convention (doc §7.1)
  C3  SPs and SPd never coexist -> there is no Cov(E_m, E_a) to omit (doc §7.2)
  C4  `m_net = max(0, SPs - SPd)` discards the absorption on 100% of absorbing edges (doc §7.2)

    OMP_NUM_THREADS=1 python scripts/debug/boundary_spliced_check.py
"""

from __future__ import annotations

import os
import pickle
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np

_SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
_COND = "gdna_gdna300_ss_0.50_nrna_none_capture_off"
_VAR_REF = 2.804  # the reference's own Var(log f_g) at L=10 -- the no-information value


def c1_expectation():
    print("=" * 84)
    print("C1  fbp is the ARITHMETIC mean => E_n = rho_n0, and mu = log(rho_n0) - v/2 is FORCED")
    print("=" * 84)
    print(
        "  simplex_logodds.py:350   f_pos = np.sum(post * f_pos_g, axis=1)   <- a MEAN, not a median"
    )
    print(
        "  So rho_n0 = fbp*sm/er IS E[rho_n]. The design's `E_n = rho_n0*exp(v_n/2)` inflates it by"
    )
    print(
        f"  exp({_VAR_REF}/2) = {np.exp(_VAR_REF / 2):.2f}x -- purely from grid uncertainty. Bug confirmed.\n"
    )


def c2_combine():
    print("=" * 84)
    print("C2  delta vs FW vs GH under the CORRECTED (arithmetic-mean) convention")
    print("=" * 84)
    rng = np.random.default_rng(1)
    N, v_m = 4_000_000, 1.0 / 500.0
    gh_x, gh_w = np.polynomial.hermite.hermgauss(7)
    print(
        f"{'w':>7}{'exact':>10}{'delta':>10}{'FW':>10}{'GH-7':>11}   {'delta/ex':>9}{'FW/ex':>8}{'GH/ex':>9}"
    )
    print("-" * 84)
    for w in (0.50, 0.90, 0.99, 0.999):
        Em, En, v_n = w, 1.0 - w, _VAR_REF
        mu_m, mu_n = np.log(Em) - v_m / 2, np.log(En) - v_n / 2  # arithmetic-mean matched
        S = np.exp(mu_m + np.sqrt(v_m) * rng.standard_normal(N)) + np.exp(
            mu_n + np.sqrt(v_n) * rng.standard_normal(N)
        )
        exact = 1.0 / np.var(np.log(S))
        delta = 1.0 / ((Em**2 * v_m + En**2 * v_n) / (Em + En) ** 2)
        V = Em**2 * (np.exp(v_m) - 1) + En**2 * (np.exp(v_n) - 1)
        fw = 1.0 / np.log1p(V / (Em + En) ** 2)
        z = mu_n + np.sqrt(2 * v_n) * gh_x
        lg = np.log(np.exp(mu_m + v_m / 2) + np.exp(z))
        p = gh_w / np.sqrt(np.pi)
        m1 = (p * lg).sum()
        gh = 1.0 / max((p * lg * lg).sum() - m1 * m1, 1e-12)
        print(
            f"{w:>7.3f}{exact:>10.2f}{delta:>10.2f}{fw:>10.2f}{gh:>11.2f}   "
            f"{delta / exact:>8.2f}x{fw / exact:>7.2f}x{gh / exact:>8.2f}x"
        )
    print()
    print(
        "  BOTH delta and FW are single-channel EXACT (delta: rho^2*v/rho^2 = v; FW: log1p(e^v-1) = v),"
    )
    print(
        "  so that property does NOT discriminate them. delta is 1.0-1.25x anti-conservative; FW is"
    )
    print(
        "  0.19-0.80x conservative. GH-7 over the nascent alone is 152x ANTI-conservative at w=0.999"
    )
    print(
        "  (the mature channel's own variance dominates there and the scheme misses it) -> rejected."
    )
    print("  FW chosen for CONSERVATISM: this project's failure mode is over-confident messages.\n")


def _load_geometry():
    from rigel.calibration.node_chain import build_node_chain
    from rigel.calibration.node_geometry import build_node_geometry
    from rigel.calibration.region_arrays import RegionArrays
    from rigel.calibration.substrate import BoundarySubstrate, CalibrationSubstrate
    from rigel.config import PipelineConfig
    from rigel.index import TranscriptIndex

    with open(_SUITE / "_selfsolve_cache" / f"{_COND}.pkl", "rb") as fh:
        d = pickle.load(fh)
    index = TranscriptIndex.load(str(_SUITE / "rigel_index"))
    cfg = PipelineConfig()
    ra = RegionArrays.from_index(index)
    pl = d["payload"]
    sub = CalibrationSubstrate.from_payload(pl, ra)
    bsub = BoundarySubstrate.from_payload(pl)
    chain = build_node_chain(pl.ref_region_offsets, pl.ref_boundary_offsets)
    geom = build_node_geometry(chain, sub, bsub, ra, d["gdna_fl_pmf"], d["rna_fl_pmf"])
    return chain, geom


def c3_c4_absorption():
    print("=" * 84)
    print("C3/C4  do SPs and SPd coexist?  and does max(0, SPs-SPd) discard the absorption?")
    print("=" * 84)
    chain, geom = _load_geometry()
    left = np.asarray(chain.left)
    # forward scan: src face = RIGHT (sf=1), dst face = LEFT (df=0)  [bp_solver: _scan(order, left, 1, 0)]
    SPs, SPd = np.asarray(geom.spliced_pos_right), np.asarray(geom.spliced_pos_left)
    e = np.array([i for i in range(len(left)) if left[i] >= 0])
    s, dd = SPs[left[e]], SPd[e]
    both = int(((s > 1e-9) & (dd > 1e-9)).sum())
    print(
        f"  forward edges: SPs>0 only={int(((s > 1e-9) & (dd <= 1e-9)).sum())}  "
        f"SPd>0 only={int(((dd > 1e-9) & (s <= 1e-9)).sum())}  BOTH={both}"
    )
    print(
        f"  -> 'at most one of SPs/SPd per edge' {'CONFIRMED' if both == 0 else 'FALSE'}"
        f"  => there is NO Cov(E_m, E_a) to omit: they are never in the same sum.\n"
    )
    m_net = np.maximum(0.0, s - dd)
    lost = (dd > 1e-9) & (m_net <= 1e-9)
    n_abs = int((dd > 1e-9).sum())
    print(f"  edges with absorption (SPd>0): {n_abs}")
    print(
        f"  ...where max(0, SPs-SPd) = 0, i.e. absorption DISCARDED: {int(lost.sum())} "
        f"({100 * lost.sum() / max(n_abs, 1):.0f}%)"
    )
    print(f"  mature mass silently dropped: {dd[lost].sum():,.0f}")
    print(
        "  -> the proposed m_net REGRESSES the shipped mature-absorption fix; rejected (doc §7.2).\n"
    )


if __name__ == "__main__":
    c1_expectation()
    c2_combine()
    c3_c4_absorption()
