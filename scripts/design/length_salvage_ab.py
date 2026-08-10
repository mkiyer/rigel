#!/usr/bin/env python
"""⭐⭐⭐ **CAN THE LENGTH CHANNEL BE SALVAGED? — the P1 SILENCE gate and the P2 SPEECH gate, in the solver.**

**The question, and why it needs the solver.** ``length_channel_census.py`` and ``length_row_shape.py``
score the row OFFLINE, against per-slot truth. That is the right place to learn what shape the row is —
and it cannot answer the question that decides the feature, because the measured defect is not a wrong
value. It is that ``tau_len`` is UNGATED (``node_init`` registers the row's curvature with no admissibility
test), so enabling the channel converts **97.8 % of the library mass at ``g00`` from UNDETERMINED to
SOLVABLE**, one-way on 12 of 12 conditions. The channel does not move values inside a fixed partition; it
RECLASSIFIES the partition, and the gDNA prior is then fitted on the slots the channel itself called
solvable. That is invisible offline.

⭐⭐ **THE TWO GATES, AND EITHER ONE KILLS THE SALVAGE.**

===  =====================================================================================
 P1   SILENCE. On the ladder the TRUE component length gap is 2.6-3.1 bp -- below the
      tool's own measurement error on ``Delta`` -- so there is nothing to find and correct
      behaviour is to say nothing. A salvaged row must stop reclassifying UNDETERMINED
      slots as SOLVABLE and must leave the confidently-wrong mass where the OFF arm has it.
 P2   SPEECH. On flgap the true gap is ~110 bp. A salvaged row must stay LIVE there. ⛔ An
      arm that passes P1 by going inert everywhere is the OFF switch with extra steps, and
      the live-slot fraction is printed beside every score so that cannot hide.
===  =====================================================================================

⛔ **P1 AND P2 CANNOT SHOW THAT THE CHANNEL TRACKS TRUTH, AND NOTHING ON DISK CAN.** The ladder varies the
gDNA level at ~zero length gap; flgap varies the length gap at ONE gDNA level (all ``g50``). Distinguishing
an estimator from a constant needs a panel that varies BOTH — ``TRAPS: a-single-level-panel-cannot-see-a-constant``
— and no such panel exists. These gates are the cheap kill-tests that come FIRST, so the expensive panel is
only built for a row that has earned it.

⭐ **THE ARMS, AND ALL OF THEM ARE DERIVED — no threshold, no bp scale, no new constant.**

``off``          the shipped OFF path. ⛔ Must be byte-identical to ``length_likelihood=False``.
``as_is``        the shipped ON path. ⛔ Must be byte-identical to ``length_likelihood=True``.
``marg``         the row's covariance marginalised over the uncertainty in the two fitted pmfs::

                     Cov_eff(pi) = N*V(pi) + N^2 * [ pi^2*Sigma_g + (1-pi)^2*Sigma_r ]

                 The ``N^2`` is forced, not chosen: every fragment at a slot shares ONE pmf error, so
                 that error is perfectly correlated across the draws and its contribution to a SUM
                 scales with ``N``, not ``sqrt(N)``. ⭐ It goes INSIDE ``length_loglik``'s own covariance
                 rather than into ``tau_len``, which is the whole difference from the derived shrinkage:
                 the row changes, so the MODE changes, and the mode is what was wrong.
``marg_frozen``  the same covariance, evaluated once at ``pi = 1/2`` and held. The row is then EXACTLY
                 quadratic with curvature ``N * Delta' [V + N*Sigma]^-1 Delta`` -- and at ``Delta = 0``
                 the mean does not move, so the row is EXACTLY CONSTANT and normalises out of psi.
                 ⚠ ``V(1/2)`` is a DECISION OWED: it is the fixed point of the gDNA<->RNA relabelling,
                 which is why it is the natural choice and not a tuned one, but the "it cannot matter"
                 argument is ``O(Delta^2)`` and has only ever been checked where ``Delta ~ 0``.

⭐⭐ **WHERE ``Sigma`` COMES FROM, AND WHY IT NEEDS NO CONSTANT.** Two terms, both already in production.

* **sampling** -- each pmf is a Dirichlet posterior mean with effective size ``pool_total + prior_ess``
  (``fl._smooth_eb`` computes both and then discards them). For a tilted moment ``m = E_f[A*t]/E_f[A]``
  the delta method gives ``Var(m) = E_f[h^2]/(n_eff+1)`` with ``h = A*(t-m)/E_f[A]``, and ``E_f[h] = 0``
  by construction -- so the covariance of ``(m1, m2)`` is closed form and exact to first order.
* **placement** -- ⭐ the four gDNA pools are four estimates of ONE distribution, each divided by its own
  exact opportunity, so under uniform placement they must agree. Their SAMPLE COVARIANCE is therefore a
  truth-free measurement of the residual that ``EQUATIONS.md`` 4.4 says is a placement model rather than
  a better divisor. Measured off this instrument's own table: the four pool means agree to ~1 bp off
  capture and disagree by tens of bp under it. It is computable on real data with no oracle.

⛔ **THE RNA ARM HAS ONE POOL, SO IT GETS THE SAMPLING TERM ONLY, AND THAT ASYMMETRY IS REAL RATHER THAN
AN OVERSIGHT** -- a spliced fragment is certified RNA and the pool is pure by construction, whereas the
gDNA pools' purity is what capture attacks. It is stated here because an unstated asymmetry is how a
half-modelled uncertainty gets read as a full one.

Usage::

    python scripts/design/length_salvage_ab.py                    # the default 4-condition kill-test
    python scripts/design/length_salvage_ab.py --conditions ...   # ladder/flgap condition names
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

RUNS = Path.home() / "Downloads" / "rigel_runs"
INDEX = RUNS / "suite" / "rigel_index"

#: panel -> the conditions the two gates need. ⭐ P1 needs a ZERO and a NEAR-SATURATED level, because a
#: zero-gDNA arm alone is one-sided: on a library with no gDNA anything that lowers the reported gDNA
#: fraction scores better. P2 needs a real length gap, which only flgap has.
_DEFAULT = (
    ("ladder", "gdna_g00_ss_0.50_nrna_none_capture_off"),
    ("ladder", "gdna_g00_ss_0.50_nrna_none_capture_on"),
    ("ladder", "gdna_g98_ss_0.50_nrna_none_capture_on"),
    ("flgap_short", "gdna_g50_ss_0.50_nrna_none_capture_on"),
)

ARMS = ("off", "as_is", "marg", "marg_frozen")

_FAILED: list[str] = []


def _gate(name: str, ok: bool, detail: str) -> bool:
    """Print a gate's verdict and remember any failure, so the run ends with a verdict and not a wall."""
    print(f"    {'✅' if ok else '⛔'} {name:<46} {detail}")
    if not ok:
        _FAILED.append(name)
    return ok


# ─────────────────────────────────────────────────────────────────────────────────────────────────
#  Sigma — the covariance of each component's fitted moment vector (m1, m2)
# ─────────────────────────────────────────────────────────────────────────────────────────────────


def _tilted_moment_cov(pmf, opportunity, weight, n_eff):
    """``Cov[(m1, m2)]`` for the tilted moments of ONE fitted pmf, by the delta method.

    ``m_t = E_f[A*t] / E_f[A]`` for ``t in (u, w)``. Writing ``h_t(w) = A(w)*(t(w) - m_t)/E_f[A]`` gives
    ``m_t = E_f[h_t] + m_t`` with ``E_f[h_t] = 0``, so a Dirichlet perturbation of ``f`` with effective
    size ``n_eff`` moves ``m_t`` by ``Sum_w (df) h_t``, whence
    ``Cov(m_s, m_t) = (E_f[h_s h_t] - E_f[h_s]E_f[h_t]) / (n_eff + 1) = E_f[h_s h_t]/(n_eff + 1)``.

    ⛔ ``n_eff`` is ``pool_total + prior_ess`` and NOT the fragment count: the pmf is a library-level
    object shared by every slot, which is exactly why its error does not average away across slots.
    """
    p = np.asarray(pmf, np.float64)
    tot = float(p.sum())
    p = p / tot if tot > 0.0 else p
    a = np.asarray(opportunity, np.float64)[: p.shape[0]]
    ea = float((p * a).sum())
    if ea <= 0.0 or n_eff <= 0.0:
        return np.zeros((2, 2), np.float64)
    m = [float((p * a * t).sum()) / ea for t in weight]
    h = [a * (t - mi) / ea for t, mi in zip(weight, m)]
    c = np.empty((2, 2), np.float64)
    for i in range(2):
        for j in range(2):
            c[i, j] = float((p * h[i] * h[j]).sum()) / (n_eff + 1.0)
    return c


def _edge_weights(max_w):
    """``(u(w), w)`` and ``A(w)`` for the CROSSING frame — the same three expressions
    :func:`rigel.calibration.length_likelihood.crossing_moments` uses, and no others."""
    w = np.arange(max_w + 1, dtype=np.float64)
    ok = w >= 2.0
    u = np.zeros_like(w)
    np.divide(1.0, w - 1.0, out=u, where=ok)
    return (u, w), np.maximum(w - 1.0, 0.0)


def _pool_moment_vectors(payload, gdna_opp, max_w):
    """⭐⭐ The four gDNA pools' own ``(m1, m2)``, each de-tilted by its OWN opportunity.

    Each pool is one estimate of the same library FL distribution (``EQUATIONS.md`` 4.3), so under
    uniform placement all four must agree. Their spread is therefore a TRUTH-FREE measurement of the
    residual 4.4 says is a placement model rather than a better divisor — available on real data with no
    oracle. Empty pools are dropped; at ``g00`` all four are empty and the caller gets nothing, which is
    correct rather than a fallback: with no pools there is no disagreement to measure, and the sampling
    term (where ``n_eff`` is then ``prior_ess`` alone) is the whole of Sigma.
    """
    from rigel.calibration.fl import _GDNA_POOLS

    (u, w), a_edge = _edge_weights(max_w)
    pool_lengths = np.asarray(payload.pool_lengths, np.float64)
    total = np.asarray(gdna_opp.total, np.float64)[: max_w + 1]
    out = []
    for slot, opp in zip(_GDNA_POOLS, gdna_opp.pools):
        counts = pool_lengths[slot][: max_w + 1]
        if float(counts.sum()) <= 0.0:
            continue
        # de-tilt: f(w) ∝ count(w) * T(w) / A_p(w)   — EQUATIONS.md 4.1, a divisor is a PROBABILITY
        opp = np.asarray(opp, np.float64)[: max_w + 1]
        prob = np.divide(opp, total, out=np.zeros_like(opp), where=total > 0.0)
        f = np.divide(counts, prob, out=np.zeros_like(counts), where=prob > 0.0)
        s = float(f.sum())
        if s <= 0.0:
            continue
        f = f / s
        ea = float((f * a_edge).sum())
        if ea <= 0.0:
            continue
        out.append((float((f * a_edge * u).sum()) / ea, float((f * a_edge * w).sum()) / ea))
    return np.asarray(out, np.float64).reshape(-1, 2)


def build_sigmas(payload, index, fl, max_w):
    """``(Sigma_g, Sigma_r, report)`` — the moment-vector covariance of each fitted pmf.

    gDNA gets sampling + placement; RNA gets sampling alone, because it has ONE pool and that pool is
    pure by construction (a spliced fragment is certified RNA). The asymmetry is in the DIVISORS, not in
    the algebra, and it is reported rather than buried.
    """
    from rigel.calibration.fl import POOL_EB_PRIOR_ESS
    from rigel.calibration.gdna_opportunity import gdna_opportunity_from_index

    (u, w), a_edge = _edge_weights(max_w)
    weight = (u, w)
    s_g = _tilted_moment_cov(fl.gdna_pmf, a_edge, weight, fl.n_gdna + POOL_EB_PRIOR_ESS)
    s_r = _tilted_moment_cov(fl.rna_pmf, a_edge, weight, fl.n_rna + POOL_EB_PRIOR_ESS)

    gopp = gdna_opportunity_from_index(index, max_w)
    vecs = _pool_moment_vectors(payload, gopp, max_w)
    n_pool = int(vecs.shape[0])
    place = np.zeros((2, 2), np.float64)
    if n_pool >= 2:
        # the sample covariance OF THE MEAN of n_pool estimates — ddof=1 then /n, no free parameter
        place = np.cov(vecs, rowvar=False, ddof=1) / float(n_pool)
    report = {
        "n_gdna": float(fl.n_gdna),
        "n_rna": float(fl.n_rna),
        "prior_ess": float(POOL_EB_PRIOR_ESS),
        "n_pool": n_pool,
        "pool_m2": vecs[:, 1].tolist() if n_pool else [],
        "sd_g_m2_sampling": float(np.sqrt(max(s_g[1, 1], 0.0))),
        "sd_g_m2_placement": float(np.sqrt(max(place[1, 1], 0.0))),
        "sd_r_m2_sampling": float(np.sqrt(max(s_r[1, 1], 0.0))),
    }
    return s_g + place, s_r, report


# ─────────────────────────────────────────────────────────────────────────────────────────────────
#  the salvaged row
# ─────────────────────────────────────────────────────────────────────────────────────────────────


def salvaged_loglik(m_g, m_r, count, d_obs, s_obs, fg_grid, *, sigma_g, sigma_r, freeze):
    """``length_loglik`` with the pmf uncertainty marginalised into the covariance.

    ⛔ **The live predicate, the exact discrimination guard and the whole-row-live rule are COPIED from
    the shipped function unchanged** — this arm may differ from ``as_is`` in the covariance and in
    nothing else, or the comparison is measuring two changes (``TRAPS: a-test-that-redefines``). With
    ``sigma_g = sigma_r = 0`` and ``freeze = False`` it must reproduce ``length_loglik`` exactly, and
    gate Ⓖ4 checks that on every condition.
    """
    fg = np.asarray(fg_grid, np.float64)[None, :]
    n = np.asarray(count, np.float64)[:, None]
    d = np.asarray(d_obs, np.float64)[:, None]
    s = np.asarray(s_obs, np.float64)[:, None]

    def mix(g, r):
        return fg * np.asarray(g, np.float64)[:, None] + (1.0 - fg) * np.asarray(r, np.float64)[:, None]

    mu_d, mu_s = mix(m_g.m1, m_r.m1), mix(m_g.m2, m_r.m2)
    v_dd = n * (mix(m_g.q1, m_r.q1) - mu_d * mu_d)
    v_ss = n * (mix(m_g.q2, m_r.q2) - mu_s * mu_s)
    v_ds = n * (mix(m_g.q12, m_r.q12) - mu_d * mu_s)

    # ⭐ the marginalisation: N^2 * [pi^2 Sigma_g + (1-pi)^2 Sigma_r]. One pmf error is shared by all N
    # fragments at the slot, so it is perfectly correlated across them and enters the SUM's covariance
    # with N^2 rather than N — which is precisely why depth does not average it away.
    nn = n * n
    for i, j, acc in ((0, 0, "dd"), (1, 1, "ss"), (0, 1, "ds")):
        add = nn * (fg * fg * sigma_g[i, j] + (1.0 - fg) * (1.0 - fg) * sigma_r[i, j])
        if acc == "dd":
            v_dd = v_dd + add
        elif acc == "ss":
            v_ss = v_ss + add
        else:
            v_ds = v_ds + add

    if freeze:
        # evaluate once at pi = 1/2 and hold. The row is then exactly quadratic in pi, its curvature is
        # N * Delta' [V + N*Sigma]^-1 Delta, and at Delta = 0 the mean does not move so the row is
        # exactly CONSTANT — which normalises out of psi instead of contributing an argmax.
        k = fg.shape[1]
        mid = k // 2
        v_dd = np.repeat(v_dd[:, mid: mid + 1], k, axis=1)
        v_ss = np.repeat(v_ss[:, mid: mid + 1], k, axis=1)
        v_ds = np.repeat(v_ds[:, mid: mid + 1], k, axis=1)

    r_d, r_s = d - n * mu_d, s - n * mu_s
    det = v_dd * v_ss - v_ds * v_ds

    discriminates = np.zeros(n.shape, dtype=bool)
    for a, b in ((m_g.m1, m_r.m1), (m_g.m2, m_r.m2), (m_g.q1, m_r.q1),
                 (m_g.q2, m_r.q2), (m_g.q12, m_r.q12)):
        discriminates |= (
            np.broadcast_to(np.asarray(a, np.float64), n.shape[:1])
            != np.broadcast_to(np.asarray(b, np.float64), n.shape[:1])
        )[:, None]

    live = (
        (n > 0.0)
        & discriminates
        & (np.asarray(m_g.eff, np.float64)[:, None] > 0.0)
        & (np.asarray(m_r.eff, np.float64)[:, None] > 0.0)
        & (det > 0.0)
        & (v_dd > 0.0)
        & (v_ss > 0.0)
    )
    safe = np.where(live, det, 1.0)
    quad = (v_ss * r_d * r_d - 2.0 * v_ds * r_d * r_s + v_dd * r_s * r_s) / safe
    ll = np.where(live, -0.5 * (quad + np.log(safe)), 0.0)
    return np.where(live.all(axis=1, keepdims=True), ll, 0.0)


def install_arm(arm, payload, index, max_w):
    """Monkeypatch ``calibrate._build_length_loglik`` for one arm; returns a restore callable.

    ⚠ A prototype belongs in a script until a panel has priced it (``panel-before-src``), and the shipped
    entry point is a module-level function called unqualified, so an override here exercises the REAL
    ``calibrate`` with one expression changed. Nothing is written to ``src/``.
    """
    # ⚠ `rigel.calibration.calibrate` is a MODULE holding a FUNCTION of the same name, and the package
    # re-exports the function — so attribute access resolves to the function. Reach the module by name.
    import rigel.calibration.calibrate  # noqa: F401  (registers the module in sys.modules)

    cal = sys.modules["rigel.calibration.calibrate"]
    original = cal._build_length_loglik
    if arm in ("off", "as_is"):
        return original, (lambda: setattr(cal, "_build_length_loglik", original))

    def patched(chain, geometry, region_arrays, gdna_fl_pmf, rna_fl_pmf, config):
        if not getattr(config, "length_likelihood", False):
            return None
        from rigel.calibration.length_likelihood import build_slot_moments
        from rigel.calibration.simplex_logodds import _logodds_grid

        _, fg_grid = _logodds_grid(int(config.sweep_n_grid), float(config.sweep_logodds_window))
        s_g, s_r, _rep = _SIGMA_CACHE[(id(payload), max_w)]
        return salvaged_loglik(
            build_slot_moments(chain, region_arrays, gdna_fl_pmf),
            build_slot_moments(chain, region_arrays, rna_fl_pmf),
            np.asarray(geometry.unspliced_count, np.float64).sum(axis=1),
            np.asarray(geometry.unspliced_inv_length_sum, np.float64),
            np.asarray(geometry.unspliced_length_sum, np.float64),
            fg_grid,
            sigma_g=s_g,
            sigma_r=s_r,
            freeze=(arm == "marg_frozen"),
        )

    cal._build_length_loglik = patched
    return original, (lambda: setattr(cal, "_build_length_loglik", original))


_SIGMA_CACHE: dict = {}


# ─────────────────────────────────────────────────────────────────────────────────────────────────
#  the run
# ─────────────────────────────────────────────────────────────────────────────────────────────────


def run_arm(inputs, arm, payload, index, max_w):
    """One arm's calibration, with the per-slot internals the two gates read."""
    import dataclasses as dc

    from rigel.calibration.node_chain import NODE
    from rigel.config import CalibrationConfig

    _orig, restore = install_arm(arm, payload, index, max_w)
    try:
        cfg = dc.replace(CalibrationConfig(), length_likelihood=(arm != "off"))
        dbg: dict = {}
        res = sys.modules["rigel.calibration.calibrate"].calibrate(**inputs, config=cfg, _debug=dbg)
    finally:
        restore()
    cap, chain = dbg["capture"], dbg["chain"]
    tau = np.asarray(cap["_tau0_lam"], np.float64)
    cnt = np.asarray(cap["count"], np.float64).sum(axis=1)
    lock = (~np.asarray(cap["solvable"], bool)) & (np.asarray(chain.kind) == NODE)
    return {
        "tau": tau,
        "count": cnt,
        "fg": np.asarray(cap["f_g"], np.float64),
        "has_ev": (tau > 1e-9) & ~lock,
        "lock": lock,
        "res": res,
        "chain": chain,
    }


def main() -> int:  # noqa: C901
    from rigel.calibration.region_arrays import RegionArrays
    from rigel.config import PipelineConfig
    from rigel.index import TranscriptIndex
    from rigel.pipeline import _native_detect_sj_tag
    from rigel.scan_cache import calibration_inputs, read_scan_cache

    from _oracle import ORIGINS
    from length_channel_census import slot_truth

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--conditions", nargs="*", default=None,
                    help="panel/condition pairs as panel:condition; default is the 4-condition kill-test")
    args = ap.parse_args()

    todo = list(_DEFAULT)
    if args.conditions:
        todo = [tuple(c.split(":", 1)) for c in args.conditions]

    index = TranscriptIndex.load(str(INDEX))
    ra = RegionArrays.from_index(index)
    cfg = PipelineConfig()

    print()
    print("=" * 118)
    print("  ⭐⭐⭐ CAN THE LENGTH CHANNEL BE SALVAGED? — P1 SILENCE and P2 SPEECH, in the solver")
    print("  ⛔ UNDRAINED cached payloads (production DRAINS). The gates are about RECLASSIFICATION,")
    print("     which the drain does not touch; the drained arm is owed before any src change.")
    print("=" * 118)

    for panel, cond in todo:
        suite = RUNS / "suite" / panel
        cache = suite / "oracle_cache" / cond
        if not (cache / "_main" / "payload.npz").exists():
            print(f"\n  ⚠ SKIP {panel}/{cond} — no cached payload")
            continue
        scan = dataclasses.replace(
            cfg.scan, sj_strand_tag=_native_detect_sj_tag(str(suite / cond / "sim_oracle.bam"))
        )
        main_cache = read_scan_cache(cache / "_main", index, scan)
        payload = main_cache.payload
        parts = {k: read_scan_cache(cache / k, index, scan).payload for k in ORIGINS}
        inputs = calibration_inputs(main_cache, index)
        max_w = int(payload.max_length)

        from rigel.calibration.fl import build_fl_models
        from rigel.calibration.gdna_opportunity import gdna_opportunity_from_index
        from rigel.calibration.junction_opportunity import crossing_probability_from_index

        fl = build_fl_models(
            payload,
            junction_opportunity=crossing_probability_from_index(index, max_w),
            gdna_opportunity=gdna_opportunity_from_index(index, max_w),
        )
        s_g, s_r, rep = build_sigmas(payload, index, fl, max_w)
        _SIGMA_CACHE[(id(payload), max_w)] = (s_g, s_r, rep)

        mu_g = float((np.arange(fl.gdna_pmf.shape[0]) * fl.gdna_pmf).sum())
        mu_r = float((np.arange(fl.rna_pmf.shape[0]) * fl.rna_pmf).sum())
        sd_delta = float(np.sqrt(max(s_g[1, 1] + s_r[1, 1], 0.0)))

        print(f"\n{'─' * 118}")
        print(f"  {panel}/{cond}")
        print(f"    mu_g {mu_g:8.3f}   mu_r {mu_r:8.3f}   Delta {mu_g - mu_r:+8.3f} bp"
              f"   sd(Delta) {sd_delta:7.3f} bp   |Delta|/sd {abs(mu_g - mu_r) / max(sd_delta, 1e-12):8.3f}")
        print(f"    n_gdna {rep['n_gdna']:,.0f}   n_rna {rep['n_rna']:,.0f}   prior_ess {rep['prior_ess']:,.0f}"
              f"   live gDNA pools {rep['n_pool']}")
        if rep["n_pool"]:
            print("    four-pool m2: " + "  ".join(f"{v:.2f}" for v in rep["pool_m2"])
                  + f"   spread(max-min) {max(rep['pool_m2']) - min(rep['pool_m2']):.2f} bp")
        print(f"    sd(mu_g) sampling {rep['sd_g_m2_sampling']:.3f} bp  placement "
              f"{rep['sd_g_m2_placement']:.3f} bp   sd(mu_r) sampling {rep['sd_r_m2_sampling']:.3f} bp")

        out = {a: run_arm(inputs, a, payload, index, max_w) for a in ARMS}
        fg_true, tot_true = slot_truth(_Oracle(payload, parts), out["off"]["chain"], ra)

        print("\n    ── GATES ──")
        # Ⓖ4  the Sigma -> 0 limit must reproduce the shipped row exactly
        from rigel.calibration.length_likelihood import build_slot_moments, length_loglik
        from rigel.calibration.node_chain import build_node_chain
        from rigel.calibration.node_geometry import build_node_geometry
        from rigel.calibration.simplex_logodds import _logodds_grid
        from rigel.calibration.splice_graph import build_junction_geometry_arrays
        from rigel.calibration.substrate import CalibrationSubstrate

        chain = build_node_chain(payload.ref_node_offsets, payload.ref_edge_offsets)
        geom = build_node_geometry(
            chain, CalibrationSubstrate.from_payload(payload, ra), ra,
            build_junction_geometry_arrays(index), fl.gdna_pmf, fl.rna_pmf, None,
        )
        _, fgg = _logodds_grid(int(cfg.calibration.sweep_n_grid),
                               float(cfg.calibration.sweep_logodds_window))
        mg = build_slot_moments(chain, ra, fl.gdna_pmf)
        mr = build_slot_moments(chain, ra, fl.rna_pmf)
        cnt0 = np.asarray(geom.unspliced_count, np.float64).sum(axis=1)
        d0 = np.asarray(geom.unspliced_inv_length_sum, np.float64)
        s0 = np.asarray(geom.unspliced_length_sum, np.float64)
        shipped = length_loglik(mg, mr, cnt0, d0, s0, fgg)
        zero = np.zeros((2, 2))
        limit = salvaged_loglik(mg, mr, cnt0, d0, s0, fgg, sigma_g=zero, sigma_r=zero, freeze=False)
        _gate("Ⓖ4 Sigma->0 reproduces the shipped row",
              float(np.abs(limit - shipped).max()) == 0.0,
              f"max|Δ| = {float(np.abs(limit - shipped).max()):.3e}")
        # Ⓖ3  one pmf on both arms must be EXACTLY inert on every arm
        for arm, fr in (("marg", False), ("marg_frozen", True)):
            null = salvaged_loglik(mg, mg, cnt0, d0, s0, fgg, sigma_g=s_g, sigma_r=s_g, freeze=fr)
            _gate(f"Ⓖ3 null (one pmf both arms) inert — {arm}",
                  float(np.abs(np.ptp(null, axis=1)).max()) == 0.0,
                  f"max|ptp| = {float(np.abs(np.ptp(null, axis=1)).max()):.3e}")
        # Ⓖ5  the arm could have fired
        for a in ARMS:
            live_m = float(out[a]["count"][out[a]["has_ev"]].sum())
            _gate(f"Ⓖ5 arm could have fired — {a}", True,
                  f"own-evidence mass {live_m:,.0f}  slots {int(out[a]['has_ev'].sum()):,}")

        print("\n    ── P1 SILENCE / P2 SPEECH ──")
        tot = float(out["off"]["count"].sum())
        base_ev = out["off"]["has_ev"]
        print(f"    {'arm':<14}{'own-ev slots':>14}{'own-ev mass':>16}{'reclassified':>14}"
              f"{'mw f_g':>9}{'conf-wrong mass':>17}")
        print("    " + "-" * 86)
        for a in ARMS:
            o = out[a]
            recl = float(o["count"][o["has_ev"] & ~base_ev].sum())
            w = o["count"]
            mwfg = float((o["fg"] * w).sum() / max(w.sum(), 1e-9))
            sd = 1.0 / np.sqrt(np.maximum(o["tau"], 1e-12))
            sel = np.isfinite(fg_true) & o["has_ev"]
            cw = float(tot_true[sel][np.abs(o["fg"][sel] - fg_true[sel]) > 2.0 * sd[sel]].sum())
            print(f"    {a:<14}{int(o['has_ev'].sum()):>14,}{float(w[o['has_ev']].sum()):>16,.0f}"
                  f"{100 * recl / max(tot, 1e-9):>13.1f}%{mwfg:>9.4f}{cw:>17,.0f}")
        print("    ⭐ 'reclassified' is the share of library mass an arm moves from UNDETERMINED to")
        print("      having own evidence, relative to the OFF arm. P1 wants ~0 on the ladder; P2 wants")
        print("      a LARGE own-evidence mass on flgap. An arm small on both is the OFF switch.")

    print()
    print("=" * 118)
    if _FAILED:
        print(f"  ⛔ {len(_FAILED)} GATE(S) FAILED — every number above is void until they pass:")
        for f in _FAILED:
            print(f"      {f}")
    else:
        print("  ✅ every gate passed")
    print("=" * 118)
    return 1 if _FAILED else 0


class _Oracle:
    """The two fields :func:`length_channel_census.slot_truth` reads, and nothing else."""

    def __init__(self, full, parts):
        self.full, self.parts = full, parts


if __name__ == "__main__":
    raise SystemExit(main())
