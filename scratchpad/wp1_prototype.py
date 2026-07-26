"""WEIGHTED RESCALE — derivation prototype, offline, against the solver's real PRE-PIN state.

THE PROBLEM. A message arrives at node i claiming per-component densities rho_c. The node observed M
fragments, and the identity  sum_c rho_c * E_c = M  must hold. The claim violates it (measured 1.26-1.33x on
the failing arms). Today `_pin_v` restores the identity by scaling EVERY component by the same k = M/S -- so
when one component is a good measurement and another is a bad guess, the good one is punished equally.

THE DERIVATION. Write the claim's error in log space, rho_c = rho_c^true * exp(eps_c). Choose the correction
a_c (rho_c -> rho_c * exp(a_c)) that is MOST LIKELY under the error model, subject to the identity:

    minimize   a^T Sigma^{-1} a        subject to   sum_c m_c exp(a_c) = M ,   m_c = rho_c * E_c

Linearizing the constraint about a = 0 gives  alpha^T a = delta,  alpha_c = m_c / S,  delta = log(M/S), and
the Lagrangian solution is

    a  =  delta * Sigma alpha / (alpha^T Sigma alpha)                                            (*)

so the DIRECTION of the correction is Sigma*alpha and only its MAGNITUDE is set by the constraint. The
magnitude is then re-solved EXACTLY (not linearized) along that direction, because the identity is an
identity: find mu with  sum_c m_c exp(mu * s_c) = M ,  s_c = (Sigma alpha)_c.

THE ERROR MODEL. Every component of a message shares the reframe r and the source's count, so those are
COMMON-mode; whatever is left in a component's variance is its OWN:

    Sigma = sigma_cm^2 * 11^T + diag(w),   sigma_cm^2 = sigma^2_transfer + 1/n_src,   w_c = v_c - sigma_cm^2

    ==> s_c = sigma_cm^2 + alpha_c * w_c

TWO EXACT LIMITS, and they are the whole point:
  * w -> 0 (the ONLY error is the shared frame): s_c is constant, so every component moves by the same
    factor. That is EXACTLY today's `_pin_v` -- so the shipped operator is the zero-independent-variance
    limit of this one, not a different thing.
  * sigma_cm^2 -> 0 and w_R -> 0 (no frame error, and the RNA arm is a MEASUREMENT): s_R = 0, so the RNA
    claim does not move at all and the gDNA absorbs the entire residual.

Capture-OFF drives the second limit (r ~ 1 => sigma_cm^2 ~ 0) exactly when the graft makes the RNA arm a
measurement (w_R small) -- which is the failing regime. Capture-ON drives the first (Var(log r) large).

⚠ WHERE THE VIOLATION IS CREATED (measured, and it corrects the first reading of this file). Each MESSAGE is
already pinned inside `_transport`, so each individually satisfies the identity. The combine then fuses the
two messages PER COMPONENT with per-component precisions -- gDNA weighted by the gDNA precisions, RNA by the
RNA precisions -- and a per-component weighted average of two claims that each sum to M does NOT sum to M.
**The FUSED claim is never pinned at all**, which is precisely what `_pin_v`'s own docstring says should
happen ("applied at EVERY node rather than only at the final combine"). So the variants below are applied to
the FUSED claim, and there is no common-scale term left to model: both messages were pinned to the same M, so
the residual violation is pure composition DISAGREEMENT between them.

    OMP_NUM_THREADS=1 python scratchpad/wp1_prototype.py
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
    "gdna_gdna300_ss_0.50_nrna_present_capture_off",
    "gdna_gdna300_ss_0.50_nrna_none_capture_off",
    "gdna_gdna100_ss_0.50_nrna_present_capture_off",
    "gdna_gdna300_ss_0.99_nrna_present_capture_off",
    "gdna_gdna1_ss_0.50_nrna_present_capture_on",
    "gdna_gdna300_ss_0.50_nrna_present_capture_on",
    "gdna_gdna300_ss_0.99_nrna_present_capture_on",
    "gdna_gdna300_ss_0.50_nrna_none_capture_on",
    "gdna_gdna100_ss_0.50_nrna_present_capture_verystrong",
]


def solve_mu(m, s, M, iters=60):
    """Exact magnitude along the fixed direction s: the unique root of sum_c m_c exp(mu s_c) = M.

    g(mu) is non-decreasing whenever every s_c >= 0, so a bracketed bisection is exact and unconditionally
    stable -- no step control, no cap, no tuned constant."""
    def g(mu):
        return np.sum(m * np.exp(np.clip(mu[:, None] * s, -50.0, 50.0)), axis=1)

    lo = np.full(M.shape, -1.0)
    hi = np.full(M.shape, 1.0)
    for _ in range(60):
        need_lo, need_hi = g(lo) > M, g(hi) < M
        if not (need_lo.any() or need_hi.any()):
            break
        lo = np.where(need_lo, lo * 2.0, lo)
        hi = np.where(need_hi, hi * 2.0, hi)
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        up = g(mid) < M
        lo, hi = np.where(up, mid, lo), np.where(up, hi, mid)
    return 0.5 * (lo + hi)


index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
NAMES = ["shipped", "fused_common", "fused_weighted", "fused_rna_exact"]
print(f"{'condition':<46}{'n':>5}{'reads':>11}"
      + "".join(f"{x:>16}" for x in NAMES) + f"{'claim/obs':>11}")
print("-" * (46 + 5 + 11 + 16 * len(NAMES) + 11))
for cond in CONDS:
    inp = _scan_and_truth(
        SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache"
    )
    dbg: dict = {}
    calmod.calibrate(
        inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
        np.asarray(inp["rna_fl_pmf"]),
        dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg,
    )
    chain, cap, st = dbg["chain"], dbg["capture"], dbg["statics"]
    us, uni = cap["_uni_static"], cap["_uni"][-1]
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
    T = Gp + Gn + Rp + Rn
    fo = np.where(T > _EPS, (Gp + Gn) / np.maximum(T, _EPS), np.nan)
    M, E_g, E_r = us["M"], us["E_g"], us["E_r"]
    rt, _ = _node_region_type(chain, ra)
    kind = np.asarray(chain.kind)
    CLS = {0: "intergenic", 1: "intron", 2: "exon"}
    cls = np.array([CLS.get(int(rt[j]), "?") if kind[j] == REGION else "boundary" for j in range(len(M))])
    fp, fn = np.asarray(st.free_pos, bool), np.asarray(st.free_neg, bool)

    # the FUSED claim, exactly as the solver forms it, and its per-component fused precisions
    dens = np.stack([uni["cg"], uni["cp"], uni["cn"]], axis=1)
    prec = np.stack([uni["pg"], uni["pp"], uni["pn"]], axis=1)  # the fused per-component precisions
    E = np.stack([E_g, E_r, E_r], axis=1)
    m = dens * E
    S = m.sum(axis=1)
    ok = (S > _EPS) & (M > _EPS)
    alpha = m / np.maximum(S, _EPS)[:, None]
    v = np.where(prec > 0.0, 1.0 / np.maximum(prec, _EPS), 0.0)
    sup = prec > 0.0

    out = {"shipped": dens}
    k = np.where(ok, M / np.maximum(S, _EPS), 1.0)
    out["fused_common"] = dens * k[:, None]
    for name, sdir in (
        ("fused_weighted", alpha * v),                        # Sigma = diag(v); no common term survives
        ("fused_rna_exact", np.stack([np.ones_like(M), np.zeros_like(M), np.zeros_like(M)], axis=1)),
    ):
        sd = np.where(sup, sdir, 0.0)
        dead = ~(sd > _EPS).any(axis=1)
        sd = np.where(dead[:, None], np.ones_like(sd), sd)
        mu = solve_mu(m, sd, np.where(ok, M, S))
        out[name] = dens * np.exp(np.clip(mu[:, None] * sd, -50.0, 50.0))

    sel = ((cls == "exon") & ~(fp & fn) & np.isfinite(fo) & (M > _EPS)
           & np.asarray(cap["solvable"], bool) & ok)
    row = f"{cond[5:]:<46}{int(sel.sum()):>5}{M[sel].sum():>11,.0f}"
    for kname in NAMES:
        d3 = out[kname]
        num = d3[:, 0] * E_g
        den = num + (d3[:, 1] + d3[:, 2]) * E_r
        f = np.where(den > _EPS, num / np.maximum(den, _EPS), np.nan)
        m2 = sel & np.isfinite(f)
        row += f"{np.average(np.abs(f[m2] - fo[m2]), weights=M[m2]):>16.4f}"
    row += f"{np.average((S / np.maximum(M, _EPS))[sel], weights=M[sel]):>10.2f}x"
    print(row)
