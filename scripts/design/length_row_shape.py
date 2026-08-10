#!/usr/bin/env python
"""⭐⭐⭐ **WHAT SHAPE IS THE LENGTH ROW, AND WHICH REBUILD OF IT GOES INERT WHEN THE GAP GOES TO ZERO?**
— no solver runs, everything off the origin-split oracle caches.

⛔ **THE QUESTION, EXACTLY.** ``length_channel_census.py`` established that the channel FIRES and that its
claim is false at ``g00`` (bias +0.66, median |Δ| = 1.0000 against a truth of exactly 0, precision pinned
at ``1/Var(λ grid)``). It never said what the row LOOKS like, and two hypotheses were live with different
repairs: **a mode near the MIDDLE of the grid regardless of the data** — a broad hump, repaired by the
precision (``EQUATIONS.md`` §3d) — or **a MONOTONE RAMP whose argmax is a grid END**, a hard 0-or-1 vote
repaired only by the row's FUNCTIONAL FORM, because a 2:1 mass split between the two ends averages to the
~0.5 the library reports without any slot ever believing 0.5.

⭐ **PART A settles that.** Ⓐ2 is the decisive table; Ⓐ4 splits the row into ``−½·quad`` and ``−½·log det``
so the argmax has an owner, and Ⓐ7 puts the identifying gap Δ next to every verdict. ⭐ **PART B prices the
rebuilds** — each arm a different ROW on the same substrate, scored against the origin-split oracle's TRUE
per-slot ``f_g`` plus the implied library ``f_gdna``, with the live-slot fraction beside every score
because an arm that wins by being inert is not a fix (`TRAPS: could-the-arm-have-fired`).

⛔⛔ **THE DECISIVE NULL IS NOT THE ONE-pmf NULL.** One pmf on both arms zeroes Δ *and* ΔV, and the shipped
discrimination gate then makes EVERY arm inert — free, re-verified, and vacuous. The null that separates
the arms is **Δ = 0 with ΔV ≠ 0**: gDNA's moments with ``m1``/``m2`` overwritten by RNA's, leaving
``q1``/``q2``/``q12`` gDNA's. The gate still passes (the q's differ), so a well-formed row MUST be flat and
the shipped one is not. Table Ⓑ2 is that arm.

Usage::

    python scripts/design/length_row_shape.py --panel ladder --conditions gdna_g00_ss_0.50_..._capture_on
    python scripts/design/length_row_shape.py --panel flgap_short --conditions ... --drain
    python scripts/design/length_row_shape.py --panel ladder --conditions ... --perturb logdet_sign
"""

from __future__ import annotations

import argparse
import dataclasses
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402

_REPO = Path(__file__).resolve().parents[2]
for _p in (_REPO / "scripts" / "design", _REPO / "tests" / "calibration"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from length_channel_census import (  # noqa: E402
    _BINS,
    build_channel,
    fisher_information,
    slot_channels,  # noqa: F401  — imported so a reader sees the census is this file's source of truth
    slot_truth,
)

RUNS = Path.home() / "Downloads" / "rigel_runs"
INDEX = RUNS / "suite" / "rigel_index"

#: the percentiles Ⓐ1/Ⓐ5 report — the two ends plus the quartiles, so a distribution is on the record
#: rather than a mean. Not tunable.
_PCT = (0, 1, 25, 50, 75, 99, 100)

#: ⛔ THE PERTURBATIONS. Each deliberately breaks the fixed code so a named gate must fire.
_PERTURB = ("none", "logdet_sign", "as_is_scale", "null_tilt", "empty_live", "var_identity",
            "delta0_mute")

_GATES: list[tuple[str, bool, str]] = []


def _gate(name: str, ok: bool, detail: str) -> bool:
    _GATES.append((name, bool(ok), detail))
    print(f"    {'✅ PASS' if ok else '⛔ FAIL'}  {name:<32} {detail}")
    return bool(ok)


# ────────────────────────────────────────────────────────────────────────────────────────────────
#  THE ROW, RE-DERIVED SO ITS TWO TERMS CAN BE SEPARATED
# ────────────────────────────────────────────────────────────────────────────────────────────────
def row_terms(m_g, m_r, count, d_obs, s_obs, fg_grid, perturb: str = "none") -> dict:
    """The SHIPPED ``length_loglik`` algebra, term by term — gated byte-identical to it by Ⓐ4.

    ⛔ Every line mirrors ``length_likelihood.length_loglik``; the only difference is that the pieces are
    RETURNED instead of summed, so ``ll = −½(quad + log det)`` can be attributed. A decomposition that
    re-derived the row its own way would be splitting a different function
    (`TRAPS: a-test-that-redefines`), which is why Ⓐ4 prints the max |Δ| against the shipped row first.
    """
    f8 = np.float64
    fg = np.asarray(fg_grid, f8)[None, :]
    n = np.asarray(count, f8)[:, None]
    d = np.asarray(d_obs, f8)[:, None]
    s = np.asarray(s_obs, f8)[:, None]

    def mix(g, r):
        return fg * np.asarray(g, f8)[:, None] + (1.0 - fg) * np.asarray(r, f8)[:, None]

    mu_d, mu_s = mix(m_g.m1, m_r.m1), mix(m_g.m2, m_r.m2)
    v_dd = n * (mix(m_g.q1, m_r.q1) - mu_d * mu_d)
    v_ss = n * (mix(m_g.q2, m_r.q2) - mu_s * mu_s)
    v_ds = n * (mix(m_g.q12, m_r.q12) - mu_d * mu_s)
    r_d, r_s = d - n * mu_d, s - n * mu_s
    del mu_d, mu_s
    det = v_dd * v_ss - v_ds * v_ds

    disc = np.zeros(n.shape, dtype=bool)
    for a, b in ((m_g.m1, m_r.m1), (m_g.m2, m_r.m2), (m_g.q1, m_r.q1),
                 (m_g.q2, m_r.q2), (m_g.q12, m_r.q12)):
        disc |= (np.broadcast_to(np.asarray(a, f8), n.shape[:1])
                 != np.broadcast_to(np.asarray(b, f8), n.shape[:1]))[:, None]

    live = ((n > 0.0) & disc
            & (np.asarray(m_g.eff, f8)[:, None] > 0.0) & (np.asarray(m_r.eff, f8)[:, None] > 0.0)
            & (det > 0.0) & (v_dd > 0.0) & (v_ss > 0.0))
    if perturb == "empty_live":  # ⛔ break: no slot may speak
        live = np.zeros_like(live)
    live_all = live.all(axis=1)
    safe = np.where(live, det, 1.0)
    quad = (v_ss * r_d * r_d - 2.0 * v_ds * r_d * r_s + v_dd * r_s * r_s) / safe
    logdet = np.log(safe)
    if perturb == "logdet_sign":  # ⛔ break: the heteroscedastic term flipped
        logdet = -logdet
    return {"quad": quad, "logdet": logdet, "live": live, "live_all": live_all, "det": det,
            "v_dd": v_dd, "v_ss": v_ss, "v_ds": v_ds, "r_d": r_d, "r_s": r_s, "n": n}


def _perdraw(m_g, m_r, pi):
    """The PER-DRAW mixture covariance at a scalar / per-slot / grid ``pi`` → ``(v_dd, v_ss, v_ds, det)``."""
    f8 = np.float64
    mu_d = pi * np.asarray(m_g.m1, f8) + (1 - pi) * np.asarray(m_r.m1, f8)
    mu_s = pi * np.asarray(m_g.m2, f8) + (1 - pi) * np.asarray(m_r.m2, f8)
    v_dd = pi * np.asarray(m_g.q1, f8) + (1 - pi) * np.asarray(m_r.q1, f8) - mu_d * mu_d
    v_ss = pi * np.asarray(m_g.q2, f8) + (1 - pi) * np.asarray(m_r.q2, f8) - mu_s * mu_s
    v_ds = pi * np.asarray(m_g.q12, f8) + (1 - pi) * np.asarray(m_r.q12, f8) - mu_d * mu_s
    return v_dd, v_ss, v_ds, v_dd * v_ss - v_ds * v_ds


def _gap(m_g, m_r):
    f8 = np.float64
    return (np.asarray(m_g.m1, f8) - np.asarray(m_r.m1, f8),
            np.asarray(m_g.m2, f8) - np.asarray(m_r.m2, f8))


def _s_grid(T, m_g, m_r, pi0):
    """``S(pi)`` on the whole grid with the covariance FROZEN at ``pi0`` (scalar or per-slot).

    ⭐ ``S`` is EXACTLY quadratic in ``pi`` by construction — ``r(pi) = x − N·mu(pi)`` is affine — so the
    completed square gives the mode and the curvature from ONE derivation. Returns ``(S, ok)``.
    """
    v_dd, v_ss, v_ds, det = _perdraw(m_g, m_r, pi0)
    ok = (det > 0.0) & (T["n"][:, 0] > 0.0)
    scale = np.where(ok, T["n"][:, 0] * det, 1.0)[:, None]
    r_d, r_s = T["r_d"], T["r_s"]
    S = (v_ss[:, None] * r_d * r_d - 2.0 * v_ds[:, None] * r_d * r_s
         + v_dd[:, None] * r_s * r_s) / scale
    return np.where(ok[:, None], S, 0.0), ok


def _gls(m_g, m_r, count, d_obs, s_obs, pi_ref):
    """One Fisher-scoring step from ``pi_ref``: ``(pihat, I_pi, ok)``, UNCLIPPED.

    ``I_pi = N·Δ' V(pi_ref)^-1 Δ`` is ``EQUATIONS.md`` §3d's information, arriving here from the mode's
    own derivation rather than beside it. ``pihat = pi_ref + U/I`` with ``U = Δ' V^-1 r(pi_ref)``.
    """
    f8 = np.float64
    n = np.asarray(count, f8)
    v_dd, v_ss, v_ds, det = _perdraw(m_g, m_r, pi_ref)
    mu_d = pi_ref * np.asarray(m_g.m1, f8) + (1 - pi_ref) * np.asarray(m_r.m1, f8)
    mu_s = pi_ref * np.asarray(m_g.m2, f8) + (1 - pi_ref) * np.asarray(m_r.m2, f8)
    r_d, r_s = np.asarray(d_obs, f8) - n * mu_d, np.asarray(s_obs, f8) - n * mu_s
    dd, ds = _gap(m_g, m_r)
    ok = (det > 0.0) & (n > 0.0)
    sd = np.where(ok, det, 1.0)
    u = (v_ss * dd * r_d - v_ds * (dd * r_s + ds * r_d) + v_dd * ds * r_s) / sd
    i_pi = n * (v_ss * dd * dd - 2.0 * v_ds * dd * ds + v_dd * ds * ds) / sd
    with np.errstate(divide="ignore", invalid="ignore"):
        pihat = np.asarray(pi_ref, f8) + np.where(i_pi > 0.0, u / np.where(i_pi > 0.0, i_pi, 1.0), np.nan)
    return pihat, np.where(ok, i_pi, 0.0), ok


# ────────────────────────────────────────────────────────────────────────────────────────────────
#  THE ARMS
# ────────────────────────────────────────────────────────────────────────────────────────────────
#: ⛔ TWO REFERENCE pi's APPEAR BELOW AND BOTH ARE **DECISIONS OWED**, not free choices:
#:   * ``frozen_half`` freezes the covariance at pi = ½ — the grid's own centre and the unique fixed point
#:     of the gDNA↔RNA relabelling involution. Not tuned, but still a choice: V_lin(½) = (V_g+V_r)/2 is a
#:     defensible competitor and they differ at O(|Δ|²).
#:   * ``frozen_step`` freezes it at the DATA's own one-step estimate, CLIPPED to [0,1] — the parameter's
#:     own domain, not a threshold. The clip is what makes V0 a covariance at all.
#: Neither introduces a numeric constant; both are owner calls this instrument prices rather than settles.
#: ``disp_profile``'s exponent is ``dof = k − 1 = 1`` (two moment coordinates, one parameter) — a
#: DIMENSION COUNT, pinned by the requirement that the curvature reduce to §3d's I_pi.
_DOF = 1.0

ARMS = ("as_is", "no_logdet", "frozen_half", "frozen_step", "gauss_I", "m1_only",
        "disp_profile", "null_onepmf")


def build_arm(name, m_g, m_r, count, d_obs, s_obs, fg_grid, lam_grid, kind, perturb="none"):
    """One arm's ROW, ``(n_slots, K)``, gated to the SHIPPED live set so live fractions compare."""
    from rigel.calibration.node_chain import EDGE

    if name == "null_onepmf":
        m_r = m_g  # one pmf on both arms — the established null, re-verified per condition
    T = row_terms(m_g, m_r, count, d_obs, s_obs, fg_grid, perturb=perturb)
    keep = T["live_all"][:, None] & T["live"]

    def out(r):
        return np.where(keep, r, 0.0)

    shipped = -0.5 * (T["quad"] + T["logdet"])
    if name in ("as_is", "null_onepmf"):
        if perturb == "as_is_scale" and name == "as_is":  # ⛔ break: the arm is no longer the arm
            shipped = shipped * (1.0 + 1e-12)
        return out(shipped), T
    if name == "no_logdet":
        return out(-0.5 * T["quad"]), T
    if name == "frozen_half":
        S, ok = _s_grid(T, m_g, m_r, 0.5)
        return out(np.where(ok[:, None], -0.5 * S, 0.0)), T
    if name == "disp_profile":
        # R5: profile an unknown mis-specification variance s² ≥ 0 on ITS OWN boundary.
        # 1 + s²hat(pi) = max(1, S(pi)/dof) ⇒ the row below. C¹ in pi, no tuned quantity.
        S, ok = _s_grid(T, m_g, m_r, 0.5)
        r = np.where(S <= _DOF, -0.5 * S,
                     -0.5 * (_DOF * np.log(np.maximum(S, 1e-300) / _DOF) + _DOF))
        return out(np.where(ok[:, None], r, 0.0)), T
    if name == "frozen_step":
        ph, _i, _ok = _gls(m_g, m_r, count, d_obs, s_obs, 0.5)
        ph = np.where(np.isfinite(ph), np.clip(ph, 0.0, 1.0), 0.5)
        S, ok = _s_grid(T, m_g, m_r, ph)
        return out(np.where(ok[:, None], -0.5 * S, 0.0)), T
    if name == "gauss_I":
        # §3d verbatim: an explicit Gaussian in λ, mode = the channel's OWN argmax (what the census
        # scores), precision = the analytic I_lam. A PRECISION-only repair, by construction.
        i_lam = fisher_information(m_g, m_r, count, fg_grid)
        mode = np.asarray(lam_grid, np.float64)[np.argmax(out(shipped), axis=1)]
        lam = np.asarray(lam_grid, np.float64)[None, :]
        return out(-0.5 * i_lam[:, None] * (lam - mode[:, None]) ** 2), T
    if name == "m1_only":
        vd = np.where(T["v_dd"] > 0.0, T["v_dd"], 1.0)
        one = -0.5 * (T["r_d"] * T["r_d"] / vd + np.log(vd))
        return out(np.where((kind == EDGE)[:, None], one, shipped)), T
    raise ValueError(name)


def _delta0(m_g, m_r):
    """gDNA's CENTRAL covariance carried onto RNA's means — Δ = 0 with ΔV ≠ 0, the decisive null.

    ⛔ **The obvious construction is inadmissible and the positive-control gate caught it.** Overwriting
    ``m1``/``m2`` with RNA's while keeping gDNA's RAW ``q``'s implies ``V_dd = q1_g − m1_r²`` and
    ``V_ss = q2_g − m2_r²`` — and on a real-gap panel (flgap_short, mu_g 127.65 vs mu_r 203.66) the second
    of those is NEGATIVE, so the shipped ``v_ss > 0`` gate silences every slot and the table read
    "✅ INERT" for a reason that has nothing to do with the row's shape.

    ⭐ So shift the mean and keep the CENTRAL covariance: ``q' = V_g + mu_r mu_r'``. Then ``V_g' = V_g``
    exactly (``ΔV = V_g − V_r``, unchanged and nonzero) and ``Δ = 0`` exactly, on every panel.
    """
    f8 = np.float64
    m1g, m2g = np.asarray(m_g.m1, f8), np.asarray(m_g.m2, f8)
    m1r, m2r = np.asarray(m_r.m1, f8), np.asarray(m_r.m2, f8)
    return type(m_g)(
        m1=m1r, m2=m2r,
        q1=np.asarray(m_g.q1, f8) - m1g * m1g + m1r * m1r,
        q2=np.asarray(m_g.q2, f8) - m2g * m2g + m2r * m2r,
        q12=np.asarray(m_g.q12, f8) - m1g * m2g + m1r * m2r,
        eff=m_g.eff)


# ────────────────────────────────────────────────────────────────────────────────────────────────
def _pfmt(x, w=12, p=4):
    return " ".join(f"{float(np.percentile(x, q)):>{w}.{p}g}" for q in _PCT)


def _split(row, sel, weight, k):
    """Mass fractions of the argmax at the LOWEST grid point, the HIGHEST, and the INTERIOR."""
    i = np.argmax(row, axis=1)[sel]
    w = weight[sel]
    tot = max(float(w.sum()), 1.0)
    return (float(w[i == 0].sum()) / tot, float(w[i == k - 1].sum()) / tot,
            float(w[(i > 0) & (i < k - 1)].sum()) / tot)


def _blab(lo, hi):
    return "ALL" if (lo, hi) == (1, 10**9) else (
        f"{lo}" if lo == hi else (f"{lo}+" if hi > 10**8 else f"{lo}-{hi}"))


def analyse(panel, condition, drain, perturb):  # noqa: C901
    from rigel.calibration.density_deconv import density_factor_precision
    from rigel.calibration.length_likelihood import build_slot_moments, length_loglik
    from rigel.calibration.node_chain import EDGE, NODE
    from rigel.calibration.region_arrays import RegionArrays
    from rigel.calibration.simplex_logodds import _logodds_grid
    from rigel.config import PipelineConfig
    from rigel.index import TranscriptIndex
    from rigel.pipeline import _drain_side_buffer, _native_detect_sj_tag
    from rigel.scan_cache import read_scan_cache

    from _oracle import ORIGINS, OracleTruth

    suite = RUNS / "suite" / panel
    bam = str(suite / condition / "sim_oracle.bam")
    index = TranscriptIndex.load(str(INDEX))
    cfg = PipelineConfig()
    ra = RegionArrays.from_index(index)
    scan = dataclasses.replace(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))

    cache = suite / "oracle_cache" / condition
    payload = read_scan_cache(cache / "_main", index, scan).payload
    parts = {k: read_scan_cache(cache / k, index, scan).payload for k in ORIGINS}
    if drain:
        from rigel.pipeline import scan_and_buffer
        from rigel.second_pass import drain as sp_drain, lift_choices

        _st, sm, _buf, _p = scan_and_buffer(bam, index, scan)
        lift: dict = {}
        payload = _drain_side_buffer(payload, index, sm, seed=cfg.second_pass_seed, _lift=lift)
        lifted, n_amb = lift_choices(lift["undrained"], [parts[k] for k in ORIGINS], lift["choices"])
        parts = {k: sp_drain(parts[k], ch, node_types=lift["node_types"], junctions=lift["junctions"])
                 for k, ch in zip(ORIGINS, lifted)}
        arm_note = f"DRAINED (production) — lift_ambiguous {int(n_amb):,}"
        oracle = OracleTruth(full=payload, parts=parts,
                             read_counts={k: -1 for k in ORIGINS}, n_ambiguous=int(n_amb))
    else:
        arm_note = "UNDRAINED cached payload (diagnostic arm — ⛔ PRODUCTION DRAINS)"
        oracle = OracleTruth.from_parts(payload, parts)

    chain, geometry, fl, ll_ship = build_channel(payload, index, ra, cfg)
    kind = np.asarray(chain.kind)
    count = np.asarray(geometry.unspliced_count, np.float64).sum(axis=1)
    d_obs = np.asarray(geometry.unspliced_inv_length_sum, np.float64)
    s_obs = np.asarray(geometry.unspliced_length_sum, np.float64)
    fg_true, tot_true = slot_truth(oracle, chain, ra)
    lam_grid, fg_grid = _logodds_grid(int(cfg.calibration.sweep_n_grid),
                                      float(cfg.calibration.sweep_logodds_window))
    k = int(fg_grid.shape[0])
    m_g = build_slot_moments(chain, ra, fl.gdna_pmf)
    m_r = build_slot_moments(chain, ra, fl.rna_pmf)
    mu_g = float((np.arange(fl.gdna_pmf.shape[0]) * fl.gdna_pmf).sum())
    mu_r = float((np.arange(fl.rna_pmf.shape[0]) * fl.rna_pmf).sum())
    all_mass = float(tot_true.sum())
    true_lib = float(np.nansum(fg_true * tot_true)) / max(all_mass, 1.0)

    print()
    print("=" * 116)
    print(f"  ⭐⭐⭐ THE LENGTH ROW'S SHAPE, AND THE RESCUE ARMS PRICED — {panel}/{condition}")
    print(f"  ARM: {arm_note}")
    print(f"  PERTURBATION: {perturb}")
    print(f"  mu_g {mu_g:.2f} bp   mu_r {mu_r:.2f} bp   gap {mu_g - mu_r:+.2f} bp   K={k}   "
          f"lam in [{lam_grid[0]:.1f},{lam_grid[-1]:.1f}]   fg_grid[0]={fg_grid[0]:.6g} "
          f"fg_grid[-1]={fg_grid[-1]:.6g}")
    print(f"  TRUE library f_gdna on the channel's own population = {true_lib:.6f}")
    print("=" * 116)

    T = row_terms(m_g, m_r, count, d_obs, s_obs, fg_grid, perturb=perturb)
    keep = T["live_all"][:, None] & T["live"]
    row = np.where(keep, -0.5 * (T["quad"] + T["logdet"]), 0.0)
    q_row = np.where(keep, -0.5 * T["quad"], 0.0)
    l_row = np.where(keep, -0.5 * T["logdet"], 0.0)
    live = T["live_all"]
    ptp = np.ptp(row, axis=1)
    ar = np.argmax(row, axis=1)
    am = fg_grid[ar]

    # ── THE GATES, FIRST AND ALWAYS ──
    print("\n  ⛔ GATES (each verified failing before it was trusted — see --perturb)")
    g0 = len(_GATES)
    dec = float(np.abs(row - np.asarray(ll_ship, np.float64)).max())
    _gate("A4-decomposition-identity", dec <= 1e-9,
          f"max|−½(quad+logdet) − shipped length_loglik| = {dec:.3e}   (need ≤ 1e-9)")
    as_is, _t = build_arm("as_is", m_g, m_r, count, d_obs, s_obs, fg_grid, lam_grid, kind, perturb)
    aid = float(np.abs(as_is - np.asarray(ll_ship, np.float64)).max())
    _gate("as_is-IS-the-shipped-arm", aid == 0.0, f"max|as_is − shipped| = {aid:.3e}   (need EXACTLY 0)")
    m_same = build_slot_moments(chain, ra, fl.gdna_pmf)
    nul = length_loglik(m_same, m_same, count, d_obs, s_obs, fg_grid)
    if perturb == "null_tilt":  # ⛔ break: the null speaks
        nul = nul + 1e-9 * fg_grid[None, :]
    nmax = float(np.abs(np.ptp(nul, axis=1)).max())
    _gate("one-pmf-null-exactly-inert", nmax == 0.0, f"max|ptp| = {nmax:.3e}   (need EXACTLY 0)")
    live_mass = float(tot_true[live].sum())
    _gate("live-set-non-empty", live_mass > 0.0,
          f"live slots {int(live.sum()):,}/{live.size:,}, live mass {live_mass:,.0f}/{all_mass:,.0f} "
          f"({100 * live_mass / max(all_mass, 1):.1f} %)")
    vg = _perdraw(m_g, m_g, 1.0)
    vr = _perdraw(m_r, m_r, 1.0)
    dd, ds = _gap(m_g, m_r)
    pig, n1 = fg_grid[None, :], T["n"]
    pq = pig * (1 - pig)
    pred = {"v_dd": n1 * (pig * vg[0][:, None] + (1 - pig) * vr[0][:, None] + pq * (dd * dd)[:, None]),
            "v_ss": n1 * (pig * vg[1][:, None] + (1 - pig) * vr[1][:, None] + pq * (ds * ds)[:, None]),
            "v_ds": n1 * (pig * vg[2][:, None] + (1 - pig) * vr[2][:, None] + pq * (dd * ds)[:, None])}
    if perturb == "var_identity":  # ⛔ break: drop the between-component term
        pred["v_ss"] = n1 * (pig * vg[1][:, None] + (1 - pig) * vr[1][:, None])
    # ⛔ AN EMPTY LIVE SET MUST **FAIL** THIS GATE, NOT PASS IT VACUOUSLY. `--perturb empty_live` found
    # that hole: with no live slot the reduction has nothing to reduce, and the obvious repair
    # (`.max(initial=0)`) would have printed max|Δ| = 0 and a green PASS on a channel that had been
    # silenced entirely (`TRAPS: could-the-arm-have-fired`).
    if live.any():
        vid = max(float(np.abs(pred[j][live] - T[j][live]).max()) for j in pred)
        vrel = vid / max(float(np.abs(T["v_ss"][live]).max()), 1e-300)
        _gate("A6-mixture-variance-identity", vrel <= 1e-9,
              f"max|Δ| = {vid:.3e} (rel {vrel:.3e}) over live slots   (need rel ≤ 1e-9)")
    else:
        _gate("A6-mixture-variance-identity", False,
              "NOT EVALUABLE — the live set is empty, so the identity has nothing to check")
    del pred
    # ⭐ THE POSITIVE CONTROL FOR Ⓑ2, and it was a GATE HOLE until `--perturb delta0_mute` found it:
    # with every arm muted at Δ = 0 the decisive-null table reads "✅ INERT" on every row and looks like
    # a clean bill of health while testing nothing (`TRAPS: could-the-arm-have-fired`). The SHIPPED arm
    # must SPEAK there — that is the defect Ⓑ2 exists to price — so a silent as_is voids the table.
    m_g0 = _delta0(m_g, m_r)
    r0_as_is, _t0 = build_arm("as_is", m_g0, m_r, count, d_obs, s_obs, fg_grid, lam_grid, kind)
    if perturb == "delta0_mute":  # ⛔ break: pretend every arm is inert at a zero gap
        r0_as_is = np.zeros_like(r0_as_is)
    d0max = float(np.ptp(r0_as_is, axis=1).max())
    _gate("delta0-positive-control", d0max > 0.0,
          f"shipped row at Δ=0, ΔV≠0: max|ptp| = {d0max:.4e}   (need > 0, else Ⓑ2 tests nothing)")
    if not all(g[1] for g in _GATES[g0:]):
        print("\n  ⛔⛔ A GATE FAILED — refusing to print the measurement. A split that is not the row, an "
              "arm that is not the arm, a null that speaks, or an empty live set makes every number "
              "below unreadable.")
        return False

    # ── Ⓐ1 ROW AMPLITUDE ──
    hdr = " ".join(f"{'p' + str(q):>12}" for q in _PCT)
    print("\n  Ⓐ1 ROW AMPLITUDE ptp(row) in NATS over live slots")
    print(f"    {'slot':<6} {'N':<8} {'slots':>9} {'mass':>13} {hdr}")
    print("    " + "-" * 126)
    for nm, tsel in (("NODE", kind == NODE), ("EDGE", kind == EDGE)):
        for lo, hi in (*_BINS, (1, 10**9)):
            s = live & tsel & (count >= lo) & (count <= hi)
            if not s.any():
                continue
            print(f"    {nm:<6} {_blab(lo, hi):<8} {int(s.sum()):>9,} "
                  f"{float(tot_true[s].sum()):>13,.0f} {_pfmt(ptp[s])}")

    # ── Ⓐ2 THE DECISIVE TABLE ──
    print("\n  Ⓐ2 ⛔⛔ THE DECISIVE TABLE — ARGMAX DISTRIBUTION over fg_grid, MASS-WEIGHTED by the slot's TRUE")
    print("      fragment count. 'a middle mode' predicts INTERIOR; 'a monotone ramp' predicts fg_lo + fg_hi.")
    print(f"    {'slot':<6} {'N':<8} {'slots':>9} {'mass':>13} {'at fg_lo':>10} {'at fg_hi':>10} "
          f"{'INTERIOR':>10} {'mass-w argmax':>15} {'verdict':>11}")
    print("    " + "-" * 106)
    for nm, tsel in (("NODE", kind == NODE), ("EDGE", kind == EDGE), ("all", np.ones_like(kind, bool))):
        for lo, hi in (*_BINS, (1, 10**9)):
            s = live & tsel & (count >= lo) & (count <= hi)
            if not s.any():
                continue
            f_lo, f_hi, f_in = _split(row, s, tot_true, k)
            mw = float((am[s] * tot_true[s]).sum()) / max(float(tot_true[s].sum()), 1.0)
            print(f"    {nm:<6} {_blab(lo, hi):<8} {int(s.sum()):>9,} {float(tot_true[s].sum()):>13,.0f} "
                  f"{f_lo:>10.4f} {f_hi:>10.4f} {f_in:>10.4f} {mw:>15.6f} "
                  f"{('GRID EDGE' if (f_lo + f_hi) > f_in else 'MIDDLE'):>11}")

    # ── Ⓐ3 MONOTONICITY ──
    print("\n  Ⓐ3 MONOTONICITY — the mass fraction of live rows with NO interior stationary point")
    print(f"    {'slot':<6} {'slots':>9} {'mass':>13} {'mono (mass)':>13} {'mono (slots)':>14}")
    print("    " + "-" * 62)
    for nm, tsel in (("NODE", kind == NODE), ("EDGE", kind == EDGE), ("all", np.ones_like(kind, bool))):
        s = live & tsel
        if not s.any():
            continue
        df = np.diff(row[s], axis=1)
        mono = (df >= 0).all(axis=1) | (df <= 0).all(axis=1)
        w = tot_true[s]
        print(f"    {nm:<6} {int(s.sum()):>9,} {float(w.sum()):>13,.0f} "
              f"{float(w[mono].sum()) / max(float(w.sum()), 1.0):>13.4f} {float(mono.mean()):>14.4f}")

    # ── Ⓐ4 TERM DECOMPOSITION ──
    print("\n  Ⓐ4 TERM DECOMPOSITION — ll = −½·quad + −½·log det. WHICH TERM OWNS THE ARGMAX?")
    aq, al = np.argmax(q_row, axis=1), np.argmax(l_row, axis=1)
    print(f"    {'slot':<6} {'N':<8} {'slots':>9} {'ptp −½quad p50':>15} {'ptp −½logdet p50':>17} "
          f"{'quad owns (mass)':>17} {'logdet owns (mass)':>19}")
    print("    " + "-" * 110)
    for nm, tsel in (("NODE", kind == NODE), ("EDGE", kind == EDGE)):
        for lo, hi in (*_BINS, (1, 10**9)):
            s = live & tsel & (count >= lo) & (count <= hi)
            if not s.any():
                continue
            w = tot_true[s]
            tw = max(float(w.sum()), 1.0)
            print(f"    {nm:<6} {_blab(lo, hi):<8} {int(s.sum()):>9,} "
                  f"{float(np.median(np.ptp(q_row[s], axis=1))):>15.4g} "
                  f"{float(np.median(np.ptp(l_row[s], axis=1))):>17.4g} "
                  f"{float(w[aq[s] == ar[s]].sum()) / tw:>17.4f} "
                  f"{float(w[al[s] == ar[s]].sum()) / tw:>19.4f}")

    # ── Ⓐ5 STANDARDISED RESIDUALS + THE CLIPPED-GLS DIAGNOSTIC ──
    print("\n  Ⓐ5 STANDARDISED MEAN RESIDUALS AT pi = ½, live slots with N ≥ 50. Many sd ⇒ the model does not")
    print("      fit and quad is a RAMP, not a peak. sqrt(S) is the Mahalanobis residual at 1 dof (a correct")
    print("      model gives ~1); sqrt(I_pi) is how many sd the WHOLE pi axis can move the predicted mean.")
    jm = int(np.argmin(np.abs(fg_grid - 0.5)))
    z_d = T["r_d"][:, jm] / np.sqrt(np.maximum(T["v_dd"][:, jm], 1e-300))
    z_s = T["r_s"][:, jm] / np.sqrt(np.maximum(T["v_ss"][:, jm], 1e-300))
    S_half, _ok = _s_grid(T, m_g, m_r, 0.5)
    ph, i_pi_h, _ok2 = _gls(m_g, m_r, count, d_obs, s_obs, 0.5)
    print(f"    {'slot':<6} {'stat':<15} {'slots':>9} {hdr}")
    print("    " + "-" * 122)
    for nm, tsel in (("NODE", kind == NODE), ("EDGE", kind == EDGE)):
        s = live & tsel & (count >= 50)
        if not s.any():
            continue
        for lbl, v in (("z_d", z_d[s]), ("z_s", z_s[s]),
                       ("sqrt(S) at ½", np.sqrt(np.maximum(S_half[s, jm], 0.0))),
                       ("sqrt(I_pi)", np.sqrt(np.maximum(i_pi_h[s], 0.0)))):
            print(f"    {nm:<6} {lbl:<15} {int(s.sum()):>9,} {_pfmt(v)}")
    print("\n    ⭐ CLIPPED GLS — mass fraction whose UNCONSTRAINED vertex pihat lies OUTSIDE [0,1]; the Δ→0")
    print("      derivation says that is the SAME EVENT as an endpoint argmax, and 'agreement' tests it.")
    print(f"    {'slot':<6} {'slots':>9} {'mass':>13} {'pihat outside [0,1]':>21} {'endpoint argmax':>17} "
          f"{'agreement':>11}")
    print("    " + "-" * 84)
    endp = (ar == 0) | (ar == k - 1)
    outside = ~((ph >= 0.0) & (ph <= 1.0))
    for nm, tsel in (("NODE", kind == NODE), ("EDGE", kind == EDGE), ("all", np.ones_like(kind, bool))):
        s = live & tsel
        if not s.any():
            continue
        w = tot_true[s]
        tw = max(float(w.sum()), 1.0)
        print(f"    {nm:<6} {int(s.sum()):>9,} {float(w.sum()):>13,.0f} "
              f"{float(w[outside[s]].sum()) / tw:>21.4f} {float(w[endp[s]].sum()) / tw:>17.4f} "
              f"{float(w[outside[s] == endp[s]].sum()) / tw:>11.4f}")

    # ── Ⓐ6 THE VARIANCE GEOMETRY ──
    print("\n  Ⓐ6 VARIANCE GEOMETRY — det V_g vs det V_r and where log det V(pi) peaks. −½logdet is CONVEX in")
    print("      pi so its max is at an ENDPOINT: pi=1 iff det V_g < det V_r, else pi=0. Its ptp is the closed")
    print("      form ½|log(det V_g/det V_r)|, and at an EDGE that is ONE number for the whole genome.")
    ell = np.zeros(int(chain.n_slots))
    obj = np.asarray(chain.obj_idx, np.int64)
    ell[kind == NODE] = np.asarray(ra.region_size_bp, np.float64)[obj[kind == NODE]]
    nl = np.flatnonzero(live & (kind == NODE))
    el = np.flatnonzero(live & (kind == EDGE))
    picks = []
    if el.size:
        picks.append(("EDGE (scalar)", int(el[0])))
    if nl.size:
        picks.append((f"NODE median ell={ell[nl[np.argsort(ell[nl])[nl.size // 2]]]:.0f}",
                      int(nl[np.argsort(ell[nl])[nl.size // 2]])))
    print(f"    {'slot set':<26} {'det V_g':>14} {'det V_r':>14} {'ratio g/r':>11} "
          f"{'argmax_pi logdetV':>19} {'argmin_pi':>12} {'−½logdet ptp':>14}")
    print("    " + "-" * 116)
    for lbl, i in picks:
        one = {f: np.full(1, float(np.asarray(getattr(m_g, f))[i])) for f in
               ("m1", "m2", "q1", "q2", "q12", "eff")}
        two = {f: np.full(1, float(np.asarray(getattr(m_r, f))[i])) for f in
               ("m1", "m2", "q1", "q2", "q12", "eff")}
        _a, _b, _c, dt = _perdraw(type(m_g)(**one), type(m_r)(**two), fg_grid[None, :])
        lg = np.log(np.maximum(dt[0], 1e-300))
        dg, dr = float(vg[3][i]), float(vr[3][i])
        print(f"    {lbl:<26} {dg:>14.8g} {dr:>14.8g} {dg / max(dr, 1e-300):>11.6f} "
              f"{fg_grid[int(np.argmax(lg))]:>19.6g} {fg_grid[int(np.argmin(lg))]:>12.6g} "
              f"{0.5 * float(np.abs(np.log(dg / max(dr, 1e-300)))):>14.4e}")
    print("\n    −½logdet's OWN argmax distribution (mass-weighted), for comparison with Ⓐ2:")
    print(f"    {'slot':<6} {'at fg_lo':>10} {'at fg_hi':>10} {'INTERIOR':>10} "
          f"{'slots det V_g<det V_r':>23} {'ptp spread over slots':>23}")
    print("    " + "-" * 88)
    for nm, tsel in (("NODE", kind == NODE), ("EDGE", kind == EDGE)):
        s = live & tsel
        if not s.any():
            continue
        f_lo, f_hi, f_in = _split(l_row, s, tot_true, k)
        p = np.ptp(l_row[s], axis=1)
        print(f"    {nm:<6} {f_lo:>10.4f} {f_hi:>10.4f} {f_in:>10.4f} "
              f"{float(np.mean(vg[3][s] < vr[3][s])):>23.4f} {float(np.ptp(p)):>23.4e}")

    # ── Ⓐ7 DELTA ITSELF ──
    print("\n  Ⓐ7 Δ = (m1_g−m1_r, m2_g−m2_r) — THE IDENTIFYING QUANTITY, next to every verdict above.")
    v05 = _perdraw(m_g, m_r, 0.5)
    mah_hdr = "D_V^-1_D"
    print(f"    {'slot set':<26} {'m1_g':>11} {'m1_r':>11} {'D_d':>12} {'m2_g':>10} {'m2_r':>10} "
          f"{'D_s bp':>10} {mah_hdr:>11} {'N*':>10}")
    print("    " + "-" * 118)
    for lbl, i in picks:
        a, b = float(dd[i]), float(ds[i])
        vd, vs_, vx, dt0 = (float(v05[0][i]), float(v05[1][i]), float(v05[2][i]), float(v05[3][i]))
        mah = (vs_ * a * a - 2 * vx * a * b + vd * b * b) / dt0 if dt0 > 0 else float("nan")
        nstar = 4.0 / mah if mah > 0 else float("nan")
        print(f"    {lbl:<26} {float(np.asarray(m_g.m1)[i]):>11.6g} {float(np.asarray(m_r.m1)[i]):>11.6g} "
              f"{a:>12.4g} {float(np.asarray(m_g.m2)[i]):>10.4f} {float(np.asarray(m_r.m2)[i]):>10.4f} "
              f"{b:>+10.4f} {mah:>11.4g} {nstar:>10.4g}")
    print("      N* = 4/(Δ'V^-1Δ) is the per-slot COUNT needed for sd(pihat) ≤ ½ — the depth the gap buys.")

    # ── PART B ──
    print("\n  Ⓑ1 ⭐⭐⭐ THE ARM SCOREBOARD — every row scored against the origin-split oracle's TRUE per-slot")
    print("      f_g, mass-weighted, PLUS the implied library f_gdna against the condition's own truth.")
    print("      ⛔ live% is beside every score: an arm that wins by being inert is not a fix.")
    print("      ⛔ THE |Δ| COLUMNS CANNOT SEE AN AMPLITUDE REPAIR: `disp_profile` is a strictly monotone")
    print("      transform of `frozen_half`'s S(pi), so they share an argmax BY CONSTRUCTION. What separates")
    print("      them is 'ptp cw-mean' — the amplitude psi actually adds, 43.64 nats for the shipped row at")
    print("      g00 against a 4.30-nat Beta(½,½) reference.")
    print(f"    {'arm':<13} {'live slots':>11} {'live mass%':>11} {'ptp cw-mean':>12} {'ptp p50':>10} "
          f"{'ptp p99':>11} {'bias':>9} {'mean|Δ|':>9} {'med|Δ|':>9} {'lib f_gdna':>11} {'lib err':>9} "
          f"{'END mass':>9}")
    print("    " + "-" * 146)
    scored = np.isfinite(fg_true)
    arm_rows = {}
    for nm in ARMS:
        r, _t = build_arm(nm, m_g, m_r, count, d_obs, s_obs, fg_grid, lam_grid, kind, perturb)
        lv = np.ptp(r, axis=1) > 0.0
        arm_rows[nm] = (r, lv)
        s = lv & scored
        if not s.any():
            print(f"    {nm:<13} {0:>11,} {0.0:>10.1f}%  ⛔ INERT on every slot — scores nothing")
            continue
        a2 = fg_grid[np.argmax(r, axis=1)]
        w = tot_true[s]
        tw = max(float(w.sum()), 1.0)
        err = a2[s] - fg_true[s]
        f_lo, f_hi, _fi = _split(r, s, tot_true, k)
        p = np.ptp(r[s], axis=1)
        libf = float((a2[s] * w).sum()) / tw
        cw = float((np.ptp(r, axis=1) * count)[lv].sum()) / max(float(count[lv].sum()), 1.0)
        print(f"    {nm:<13} {int(s.sum()):>11,} {100 * tw / max(all_mass, 1):>10.1f}% {cw:>12.4g} "
              f"{float(np.median(p)):>10.4g} {float(np.percentile(p, 99)):>11.4g} "
              f"{float((err * w).sum()) / tw:>+9.4f} {float((np.abs(err) * w).sum()) / tw:>9.4f} "
              f"{float(np.median(np.abs(err))):>9.4f} {libf:>11.6f} {libf - true_lib:>+9.4f} "
              f"{f_lo + f_hi:>9.4f}")
    print(f"    (TRUE library f_gdna = {true_lib:.6f}; 'END mass' = mass whose argmax is a grid END)")

    # ── Ⓑ2 THE DECISIVE NULL ──
    print("\n  Ⓑ2 ⛔⛔ THE DECISIVE NULL — Δ = 0 with ΔV ≠ 0 (gDNA's CENTRAL covariance on RNA's means). The")
    print("      shipped gate still passes because q1/q2/q12 differ, so a WELL-FORMED row MUST be exactly flat.")
    print(f"    {'arm':<13} {'live slots':>11} {'max|ptp|':>13} {'count-w mean ptp':>18} "
          f"{'END-argmax mass':>17} {'verdict':>12}")
    print("    " + "-" * 92)
    for nm in ARMS:
        if nm == "as_is":
            r0 = r0_as_is  # the gated positive control, reused so the table and the gate agree
        else:
            r0, _t0 = build_arm(nm, m_g0, m_r, count, d_obs, s_obs, fg_grid, lam_grid, kind)
        if perturb == "delta0_mute":  # ⛔ break: pretend every arm is inert at a zero gap
            r0 = np.zeros_like(r0)
        p0 = np.ptp(r0, axis=1)
        lv0 = p0 > 0.0
        mx = float(p0.max())
        cw = float((p0 * count).sum()) / max(float(count[lv0].sum()), 1.0) if lv0.any() else 0.0
        eg = _split(r0, lv0 & scored, tot_true, k) if (lv0 & scored).any() else (0.0, 0.0, 0.0)
        print(f"    {nm:<13} {int(lv0.sum()):>11,} {mx:>13.4e} {cw:>18.4g} {eg[0] + eg[1]:>17.4f} "
              f"{('✅ INERT' if mx == 0.0 else '⛔ SPEAKS'):>12}")
    print("      ⚠ A residue at 1e-13…1e-7 is NOT model structure — it is float cancellation, in the mixture")
    print("        (`fg*a + (1−fg)*a != a` bitwise, scaling with N and m2) and in this table's own re-centring")
    print("        of q1/q2/q12. ⛔ It still MATTERS: `density_factor_precision`'s live test is ptp > 1e-12, so")
    print("        anything above that reads LIVE and its near-uniform posterior returns tau = 1/Var(lam grid),")
    print("        the grid's own width sold as evidence. The cure is Δ made EXACTLY 0, not a wider tolerance.")

    # ── Ⓑ3 the one-pmf null per arm ──
    print("\n  Ⓑ3 the ONE-pmf null per arm (Δ = 0 AND ΔV = 0) — enforced by the shipped gate for every arm;")
    print("      re-verified per condition because it is free and its failure would void everything above.")
    print("    max|ptp|:  " + "  ".join(
        f"{nm}={float(np.abs(np.ptp(build_arm(nm, m_same, m_same, count, d_obs, s_obs, fg_grid, lam_grid, kind)[0], axis=1)).max()):.1e}"
        for nm in ARMS))

    # ── Ⓑ4 the precision each arm declares ──
    print("\n  Ⓑ4 THE PRECISION EACH ARM DECLARES — the tabulated tau, against §3d's analytic I_lam.")
    floor = 1.0 / float(np.var(lam_grid))
    info = fisher_information(m_g, m_r, count, fg_grid)
    print(f"      grid floor 1/Var(lam) = {floor:.6f};  median analytic I_lam over live = "
          f"{float(np.median(info[live])):.6g}")
    print(f"    {'arm':<13} {'med tau':>12} {'p99 tau':>12} {'% at floor':>12}")
    print("    " + "-" * 54)
    for nm in ARMS:
        r, lv = arm_rows[nm]
        if not lv.any():
            print(f"    {nm:<13} {'—':>12} {'—':>12} {'—':>12}")
            continue
        tt = density_factor_precision(r, lam_grid)
        print(f"    {nm:<13} {float(np.median(tt[lv])):>12.6g} {float(np.percentile(tt[lv], 99)):>12.6g} "
              f"{100 * float(np.mean(np.abs(tt[lv] - floor) < 1e-4)):>11.1f}%")
    return True


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--panel", default="ladder")
    ap.add_argument("--conditions", nargs="+", default=["gdna_g00_ss_0.50_nrna_none_capture_on"])
    ap.add_argument("--drain", action="store_true",
                    help="DRAIN as production does. ⛔ re-scans the BAM (~1-2 min/condition).")
    ap.add_argument("--perturb", choices=_PERTURB, default="none",
                    help="deliberately break the fixed code so a named gate must fire")
    args = ap.parse_args()
    ok = True
    for cond in args.conditions:
        ok = analyse(args.panel, cond, args.drain, args.perturb) and ok
    print("\n  ⛔ GATE LEDGER")
    for nm, g, detail in _GATES:
        print(f"    {'PASS' if g else 'FAIL':<5} {nm:<32} {detail}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
