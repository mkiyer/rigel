"""The fragment-length composition likelihood — gates U1–U5.

     (P2) · Module: `rigel.calibration.length_likelihood`

⛔ **U1 IS VERIFIED BY EXACT ENUMERATION OVER INTEGER START POSITIONS — never against the module's own
closed form.**: a validator that calls the builder's own helper validates
nothing, and re-deriving by a *different algorithm* is what caught both real bugs in the index work. The
closed form here is a cumulative-sum identity; the oracle below is a literal loop over every start
position and every fragment length, with no tolerance beyond float round-off.
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.effective_length import (
    UNBOUNDED_REACH,
    contained_eff_length,
    crossing_eff_length,
)
from rigel.calibration.length_likelihood import (
    contained_moments,
    crossing_moments,
    length_loglik,
)

_UNB = np.full(1, UNBOUNDED_REACH)


def _pmf(pairs, max_len=64):
    """A pmf from explicit ``(length, mass)`` pairs — exact, no distributional slack."""
    p = np.zeros(max_len, dtype=np.float64)
    for w, m in pairs:
        p[w] += m
    return p / p.sum()


# --- the oracle: brute force over start positions ---------------------------------------------------


def _brute_contained(node_len: int, pmf: np.ndarray):
    """Enumerate EVERY (start, length) placement whose fragment lies wholly inside ``[0, node_len)``.

    Returns the five tilted moments plus ``E[A]``, accumulated one placement at a time. This is the
    definition; the module computes the same thing by a cumulative-sum identity.
    """
    tot = 0.0
    acc = dict(u=0.0, w=0.0, uu=0.0, ww=0.0, uw=0.0)
    for w in range(1, pmf.shape[0]):
        f = pmf[w]
        if f == 0.0:
            continue
        starts = 0
        for s in range(0, node_len + 1):  # a fragment [s, s+w) fits iff s + w <= node_len
            if s + w <= node_len:
                starts += 1
        if starts == 0:
            continue
        mass = f * starts
        tot += mass
        u = 1.0 / w  # the accumulator's node weight: inv_length_quantum(L)
        acc["u"] += mass * u
        acc["w"] += mass * w
        acc["uu"] += mass * u * u
        acc["ww"] += mass * w * w
        acc["uw"] += mass * u * w
    if tot == 0.0:
        return dict(m1=0.0, m2=0.0, q1=0.0, q2=0.0, q12=0.0, eff=0.0)
    return dict(
        m1=acc["u"] / tot,
        m2=acc["w"] / tot,
        q1=acc["uu"] / tot,
        q2=acc["ww"] / tot,
        q12=acc["uw"] / tot,
        eff=tot,
    )


def _brute_crossing(pmf: np.ndarray):
    """Enumerate every offset ``a ∈ [1, w−1]`` placing a length-``w`` fragment across a 0-bp line."""
    tot = 0.0
    acc = dict(u=0.0, w=0.0, uu=0.0, ww=0.0, uw=0.0)
    for w in range(1, pmf.shape[0]):
        f = pmf[w]
        if f == 0.0:
            continue
        starts = sum(1 for a in range(1, w))  # bases on BOTH sides ⇒ w − 1 offsets
        if starts == 0:
            continue
        mass = f * starts
        tot += mass
        u = 1.0 / (w - 1.0)  # the accumulator's edge weight: inv_length_quantum(L − 1)
        acc["u"] += mass * u
        acc["w"] += mass * w
        acc["uu"] += mass * u * u
        acc["ww"] += mass * w * w
        acc["uw"] += mass * u * w
    return dict(
        m1=acc["u"] / tot,
        m2=acc["w"] / tot,
        q1=acc["uu"] / tot,
        q2=acc["ww"] / tot,
        q12=acc["uw"] / tot,
        eff=tot,
    )


_PMFS = {
    "point mass w=8": _pmf([(8, 1.0)]),
    "two point 5/20": _pmf([(5, 0.4), (20, 0.6)]),
    "three point": _pmf([(3, 0.2), (11, 0.5), (29, 0.3)]),
    "wide": _pmf([(w, 1.0) for w in range(2, 40)]),
}


# --- U1 ---------------------------------------------------------------------------------------------


@pytest.mark.parametrize("pmf_name", list(_PMFS))
@pytest.mark.parametrize("ell", [1, 5, 10, 25, 40, 151])
def test_contained_moments_match_brute_force(pmf_name, ell):
    """U1a — the node frame, against a literal loop over start positions."""
    pmf = _PMFS[pmf_name]
    got = contained_moments(np.array([float(ell)]), pmf)
    want = _brute_contained(ell, pmf)
    for name in ("m1", "m2", "q1", "q2", "q12", "eff"):
        np.testing.assert_allclose(
            getattr(got, name)[0],
            want[name],
            rtol=1e-12,
            atol=1e-12,
            err_msg=f"{name} at ell={ell}",
        )


@pytest.mark.parametrize("pmf_name", list(_PMFS))
def test_crossing_moments_match_brute_force(pmf_name):
    """U1b — the edge frame, same treatment."""
    pmf = _PMFS[pmf_name]
    got = crossing_moments(pmf)
    want = _brute_crossing(pmf)
    for name in ("m1", "m2", "q1", "q2", "q12", "eff"):
        np.testing.assert_allclose(
            float(getattr(got, name)), want[name], rtol=1e-12, atol=1e-12, err_msg=name
        )


def test_node_cross_moment_is_exactly_one():
    """⭐ ``u(w)·w = 1`` at a node, so ``q12 ≡ 1`` for any pmf and any node length — a structural
    identity of the deposit rule, not a coincidence. If it ever fails, the node deposit weight has
    stopped being ``1/L``."""
    for pmf in _PMFS.values():
        m = contained_moments(np.array([5.0, 25.0, 151.0, 4000.0]), pmf)
        live = m.eff > 0
        np.testing.assert_allclose(m.q12[live], 1.0, rtol=1e-12)


# --- U2: one implementation of the effective length, not two ----------------------------------------


@pytest.mark.parametrize("pmf_name", list(_PMFS))
def test_eff_matches_the_solver_divisor(pmf_name):
    """U2 — ``moments.eff`` must be byte-identical to what the solver divides by.

    ⛔: the prose said "the AVERAGE", the code followed the prose, and a
    sibling docstring had the right formula the whole time. If these two ever disagree there are two
    implementations of one quantity and one of them is wrong.
    """
    pmf = _PMFS[pmf_name]
    ell = np.array([1.0, 5.0, 25.0, 151.0, 1000.0])
    np.testing.assert_array_equal(contained_moments(ell, pmf).eff, contained_eff_length(ell, pmf))
    np.testing.assert_allclose(
        float(crossing_moments(pmf).eff), float(crossing_eff_length(pmf, _UNB, _UNB)[0]), rtol=1e-12
    )


# --- U3: the likelihood is maximised at the truth ---------------------------------------------------


def _grid(k=201):
    return np.linspace(1e-4, 1 - 1e-4, k)


def _simulate(pi_true, n, moments_g, moments_r, seed):
    """Draw ``n`` landed fragments from the mixture and return the two observed sums.

    ⚠ Sampling from the TILTED pmf, which is what "landed here" means — sampling the raw pmf would test
    a model the accumulator does not implement.
    """
    rng = np.random.default_rng(seed)
    k_g = int(rng.binomial(n, pi_true))
    d = s = 0.0
    for tilted, k in ((moments_g, k_g), (moments_r, n - k_g)):
        if k:
            w = rng.choice(tilted["w"], size=k, p=tilted["p"])
            d += float((1.0 / w).sum())
            s += float(w.sum())
    return d, s


def _tilted_contained(ell, pmf):
    w = np.arange(pmf.shape[0], dtype=np.float64)
    a = np.maximum(ell - w + 1.0, 0.0)
    t = pmf * a
    return {"w": w[t > 0], "p": (t / t.sum())[t > 0]}


@pytest.mark.parametrize("pi_true", [0.2, 0.5, 0.8])
def test_loglik_is_maximised_near_the_truth(pi_true):
    """U3 — with genuinely different lengths, the argmax over the grid recovers ``pi``."""
    ell = 600.0
    pmf_g, pmf_r = _pmf([(60, 1.0)], 700), _pmf([(300, 1.0)], 700)
    # ⚠ point masses give a SINGULAR covariance (no spread to measure), so use spread pmfs here —
    # the degenerate case is U5's job, not U3's.
    pmf_g = _pmf([(50, 0.3), (60, 0.4), (70, 0.3)], 700)
    pmf_r = _pmf([(280, 0.3), (300, 0.4), (320, 0.3)], 700)
    mg = contained_moments(np.array([ell]), pmf_g)
    mr = contained_moments(np.array([ell]), pmf_r)
    grid = _grid()
    hits = []
    for seed in range(24):
        d, s = _simulate(
            pi_true, 400, _tilted_contained(ell, pmf_g), _tilted_contained(ell, pmf_r), seed
        )
        ll = length_loglik(mg, mr, np.array([400.0]), np.array([d]), np.array([s]), grid)
        hits.append(grid[int(np.argmax(ll[0]))])
    assert abs(float(np.mean(hits)) - pi_true) < 0.05, (
        f"argmax averaged {np.mean(hits):.3f} against a true {pi_true}"
    )


def test_the_log_det_term_is_present_and_dominates_at_low_count():
    """⭐ **PERTURBATION TEST (P2d).** U3 passes with the ``−½ log det`` term deleted, because at
    ``N = 400`` the quadratic swamps it. That is not evidence the term is optional.

    **The scaling argument.** Across the grid the quadratic's variation goes as ``N·Δμ²/v`` while
    ``log det``'s is ``O(1)`` in ``N``. So the term is negligible at high depth and decisive at low —
    and low is where the genome lives: the median node's entire evidence is ~0.6 contained + ~4.1
    spanning fragments, and 80.5 % of partition nodes carry none at all.

    **The test.** Feed the observation that makes the residual identically zero at ``pi0``. A
    quadratic-only form then peaks at exactly ``pi0`` for every ``N``. The full likelihood does not,
    because ``Sigma`` depends on ``pi`` — and the gap must shrink monotonically as ``N`` grows.

    ⚠ **What this also documents: the Gaussian is asymptotic in ``N``, and at ``N = 1`` it is not a
    trustworthy LOCATION.** The measured displacement is 0.32 at ``N = 1``, 0.05 at 5, 0.004 at 50. The
    true single-draw likelihood of a two-component mixture is bimodal and no Gaussian can represent it;
    what the term is groping toward is that one fragment really is either gDNA or RNA. Recorded here
    rather than papered over, and scored as its own low-count stratum in the P2 A/B.
    """
    ell = 600.0
    mg = contained_moments(np.array([ell]), _pmf([(50, 0.3), (60, 0.4), (70, 0.3)], 700))
    mr = contained_moments(np.array([ell]), _pmf([(280, 0.3), (300, 0.4), (320, 0.3)], 700))
    grid = _grid(4001)
    pi0 = 0.5
    mu_d = pi0 * mg.m1[0] + (1 - pi0) * mr.m1[0]
    mu_s = pi0 * mg.m2[0] + (1 - pi0) * mr.m2[0]

    shifts = []
    for n in (1.0, 5.0, 50.0, 400.0):
        ll = length_loglik(mg, mr, np.array([n]), np.array([n * mu_d]), np.array([n * mu_s]), grid)
        shifts.append(abs(float(grid[int(np.argmax(ll[0]))]) - pi0))

    assert shifts[0] > 0.15, (
        f"at N=1 the heteroscedastic term moved the peak by only {shifts[0]:.4f} — "
        "the -0.5*log(det) term is missing or inert"
    )
    assert shifts == sorted(shifts, reverse=True), (
        f"the displacement must decay monotonically with N; got {shifts}"
    )
    assert shifts[-1] < 0.01, f"at N=400 the term should be negligible; got {shifts[-1]:.4f}"


# --- U4: identical components carry no information --------------------------------------------------


def test_identical_pmfs_give_a_flat_loglik():
    """U4 — if the two components have the SAME length distribution the channel cannot tell them apart,
    and must say so with a flat row rather than a confident one.

    ⭐ This is the module's own null. A term that moved here would be reading noise as composition.
    """
    pmf = _PMFS["three point"]
    m = contained_moments(np.array([200.0]), pmf)
    ll = length_loglik(m, m, np.array([50.0]), np.array([1.7]), np.array([600.0]), _grid())
    assert float(np.ptp(ll[0])) < 1e-9, "identical components produced a non-flat likelihood"


def test_a_flat_row_carries_no_precision():
    """And a flat row must read as ZERO evidence through the precision reader the solver uses."""
    from rigel.calibration.density_deconv import density_factor_precision

    pmf = _PMFS["three point"]
    m = contained_moments(np.array([200.0]), pmf)
    lam = np.linspace(-10, 10, 201)
    ll = length_loglik(m, m, np.array([50.0]), np.array([1.7]), np.array([600.0]), _grid())
    tau = density_factor_precision(ll, lam)
    np.testing.assert_allclose(tau, [0.0])


# --- U5: nothing to say, said as nothing ------------------------------------------------------------


@pytest.mark.parametrize(
    "label,count,eff_kill",
    [
        ("no count", 0.0, None),
        ("no gDNA opportunity", 10.0, "g"),
        ("no RNA opportunity", 10.0, "r"),
    ],
)
def test_inert_where_it_cannot_speak(label, count, eff_kill):
    """U5 — zero count, or zero opportunity for either component, must give a flat row and no nan.

    ⛔. A node shorter than one RNA fragment has RNA opportunity exactly 0
    on **21.7 % of chr22 nodes** against the measured pure pools — not a corner case.
    """
    pmf_g = _pmf([(20, 0.5), (30, 0.5)], 700)
    pmf_r = _pmf([(400, 0.5), (500, 0.5)], 700)
    ell = np.array([50.0])  # every RNA fragment is longer than the node ⇒ eff_r == 0
    mg = contained_moments(ell, pmf_g if eff_kill != "g" else pmf_r)
    mr = contained_moments(ell, pmf_r)
    ll = length_loglik(mg, mr, np.array([count]), np.array([0.4]), np.array([250.0]), _grid())
    assert np.all(np.isfinite(ll)), f"{label}: non-finite log-likelihood"
    np.testing.assert_allclose(ll, 0.0, err_msg=f"{label}: expected an inert (flat, zero) row")


def test_identical_components_are_EXACTLY_inert_not_nearly_inert():
    """⭐ **The null must be exact, not approximate — a third defect found by the P2 A/B.**

    With identical components the row is constant in ``pi``, but computed as a difference of large floats
    it lands flat only to ~1e-11. `density_deconv.density_factor_precision` tests ``ptp > 1e-12``, so such
    a row reads as LIVE and its near-uniform posterior returns ``tau = 1/Var(uniform over the lambda
    grid)`` — **the grid's own width sold as composition evidence**, which is precisely what that
    function's docstring says must never happen. Measured on the chr22 pilot before the fix: 689 slots
    with a spurious tau, max 0.02902 against ``1/Var(lambda)`` = 0.029016.

    ⚠ ``assert_allclose`` would pass at 1e-11 and miss this entirely. The assertion is **exact
    equality**, and it has to be.
    """
    from rigel.calibration.density_deconv import density_factor_precision

    pmf = _PMFS["wide"]
    m = contained_moments(np.array([200.0, 40.0, 5.0]), pmf)
    ll = length_loglik(
        m,
        m,
        np.array([50.0, 7.0, 1.0]),
        np.array([1.7, 0.3, 0.1]),
        np.array([600.0, 90.0, 12.0]),
        _grid(),
    )
    assert np.array_equal(ll, np.zeros_like(ll)), "the null must be EXACTLY zero, not merely small"
    tau = density_factor_precision(ll, np.linspace(-10, 10, 201))
    np.testing.assert_array_equal(tau, np.zeros_like(tau))
