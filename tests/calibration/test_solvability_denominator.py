"""⛔⛔ THE RANKING COLUMN MUST NOT BE GAMEABLE BY THE SOLVER KNOWING LESS — `TRAPS.md` A12, A12b.

`solvability_audit`'s headline scores the SOLVABLE population, and "solvable" is a BOOLEAN on a
CONTINUOUS quantity: `own_composition_logvar` treats any `tau_lam > 1e-9` as own evidence. On an
unstranded library the strand arm carries a genuinely nonzero but physically nil precision —
measured `I ≈ Var(κ̂)·N_eff/(p(1−p))`, i.e. roughly the node's depth over the library's spliced
depth — so that boolean flips on fitting noise at some conditions and not others.

⭐ **That is not fixable in the solver and three designs prove it** (`SESSION_HANDOFF`, `TRAPS` A12b):
a resolving-power floor at `1/(2L)²` was derived, implemented and REFUTED (τ is continuous across
the region — no empty interval, so any floor is a tuned constant); subtracting exactly `Var(κ̂)`
instead of `σ²_d` OPENS the zero-gDNA control; and propagating `Var(κ̂)` into the Schur denominator
is inert because the node's own binomial noise dominates it by 87×–5,179×.

⛔ **So the headline needs a FIXED-DENOMINATOR companion, which is exactly what A12 already
prescribes and what the table did not carry:** ``mwae`` over ALL live objects, and the raw
``Σ|err|`` in fragments. Neither has a denominator the solver can move by declining to answer.

===  ===========================================================================================
D1   ``summarise`` emits ``all_mwae`` and ``abs_err``
D2   ⭐ they are FIXED-DENOMINATOR — shrinking the solvable set leaves both BIT-IDENTICAL, while
     every existing headline field moves. This is the property the whole file exists for
D3   ``all_mwae`` is the honest mass-weighted mean over the live population (brute-forced)
D4   ONE HOME for "has own composition evidence" — the instruments import the solver's predicate
     instead of restating ``1e-9`` (`TRAPS.md` A11)
===  ===========================================================================================
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.node_init import has_own_composition_evidence, own_composition_logvar


def _fixture(n=400, seed=3):
    """A synthetic scored population: truth, prediction, mass, and a τ spanning the noise region."""
    rng = np.random.default_rng(seed)
    total = rng.uniform(10.0, 5000.0, n)
    f_true = rng.uniform(0.0, 1.0, n)
    f_pred = np.clip(f_true + rng.normal(0.0, 0.15, n), 0.0, 1.0)
    err = np.abs(f_pred - f_true) * total
    # τ spanning six decades either side of the 1e-9 boolean, which is the whole point
    tau = 10.0 ** rng.uniform(-14.0, 2.0, n)
    return dict(total=total, f_true=f_true, f_pred=f_pred, err=err, tau=tau,
                live=np.ones(n, bool), z=rng.normal(0.0, 3.0, n), gap=rng.normal(0.0, 0.4, n))


# ── D4 — one home ───────────────────────────────────────────────────────────────────────────────


def test_D4_the_evidence_predicate_has_ONE_home_and_the_instruments_import_it():
    """⭐ `TRAPS.md` A11: a gate that re-derives a definition cannot detect drift in it, so the
    definition must live in ONE place with every consumer importing it. The home is production —
    the predicate is a production concept and ``scripts/`` is deliberately not importable.

    ⛔ Three instruments used to each restate ``_EPS = 1.0e-9`` beside a comment saying it must
    match the solver. Changing the solver would have moved none of them."""
    tau = np.array([0.0, 1e-12, 1e-9, 2e-9, 1e-4, 1.0])
    got = has_own_composition_evidence(tau)
    # it is EXACTLY the predicate `own_composition_logvar` applies — read off that function, not
    # restated here: where it says "no evidence" the variance is infinite.
    v_fg, v_fr = own_composition_logvar(np.full(tau.shape, 0.4), tau, np.zeros(tau.shape, bool))
    assert np.array_equal(got, np.isfinite(v_fg))
    assert np.array_equal(got, np.isfinite(v_fr))


def test_D4_perturbation_a_DIFFERENT_predicate_stops_matching_the_solver():
    """⚠ The falsification for D4 — if a consumer picked its own cut, the identity above breaks."""
    tau = np.array([0.0, 1e-12, 1e-9, 2e-9, 1e-4, 1.0])
    theirs = tau > 1e-6  # a plausible, wrong, home-made floor
    v_fg, _ = own_composition_logvar(np.full(tau.shape, 0.4), tau, np.zeros(tau.shape, bool))
    assert not np.array_equal(theirs, np.isfinite(v_fg))


# ── D1/D3 — the two fixed-denominator fields ────────────────────────────────────────────────────


def _summarise(fx, det):
    import importlib.util
    import sys
    from pathlib import Path

    key = "solvability_audit"
    if key not in sys.modules:
        p = Path(__file__).resolve().parents[2] / "scripts" / "design" / "solvability_audit.py"
        spec = importlib.util.spec_from_file_location(key, p)
        m = importlib.util.module_from_spec(spec)
        sys.modules[key] = m
        spec.loader.exec_module(m)
    SA = sys.modules[key]
    n = fx["total"].shape[0]
    a = dict(
        determined=det, total=fx["total"], err=fx["f_pred"] - fx["f_true"], live=fx["live"],
        # ⚠ z must be WIDE enough that |z| ≥ 2 selects a real subset, or D2's "the gameable columns
        # moved" control is vacuous — `conf_wrong` would be 0 in both arms and prove nothing (A14).
        z=fx["z"], sd=np.full(n, 0.5), gap=fx["gap"], sd_lam=np.full(n, 1.0),
        f_true=fx["f_true"], ladder={"fg_loc": fx["f_pred"], "f_g": fx["f_pred"]},
    )
    a["err"] = (fx["f_pred"] - fx["f_true"]) * fx["total"]
    return SA.summarise(a)


def test_D1_D3_summarise_emits_all_mwae_and_abs_err_and_they_are_correct():
    """⭐ ``all_mwae`` is Σ(mass·|Δf_g|) / Σ(mass) over every LIVE object; ``abs_err`` is Σ|Δ gDNA
    mass| in fragments. Brute-forced against the fixture, not against the implementation."""
    fx = _fixture()
    det = has_own_composition_evidence(fx["tau"])
    s = _summarise(fx, det)
    assert "all_mwae" in s and "abs_err" in s
    err = np.abs(fx["f_pred"] - fx["f_true"]) * fx["total"]
    assert s["all_mwae"] == pytest.approx(float(err.sum() / fx["total"].sum()), rel=1e-12)
    assert s["abs_err"] == pytest.approx(float(err.sum()), rel=1e-12)


# ── D2 — the property the whole file exists for ─────────────────────────────────────────────────


def test_D2_the_new_columns_are_BIT_IDENTICAL_when_the_solvable_set_SHRINKS():
    """⭐⭐⭐ THE GATE. A12's destruction control, in miniature: make the solver "know less" by
    shrinking the determined set, and check which columns move.

    Every existing headline field is defined over the determined population, so all of them move.
    ``all_mwae`` and ``abs_err`` are defined over the LIVE population, so they must be **bit**
    identical — that is precisely what makes them safe to rank a panel on, and it is why the g75
    row could mis-rank the ladder without them."""
    fx = _fixture()
    wide = has_own_composition_evidence(fx["tau"])
    narrow = fx["tau"] > 1e-3  # the same solver, declining to answer more often
    assert narrow.sum() < wide.sum() and narrow.sum() > 0, (wide.sum(), narrow.sum())

    a, b = _summarise(fx, wide), _summarise(fx, narrow)
    # the fixed-denominator pair: BIT identical
    assert a["all_mwae"] == b["all_mwae"]
    assert a["abs_err"] == b["abs_err"]
    # and the gameable ones genuinely moved, so this is not a vacuous comparison (`TRAPS.md` A14)
    moved = [k for k in ("solvable_mass_share", "solvable_mwae", "conf_wrong_objects")
             if a[k] != b[k]]
    assert len(moved) == 3, (moved, a, b)
