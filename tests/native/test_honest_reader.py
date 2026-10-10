"""``native/honest_reader.h`` against its executable specification (``_honest_reader_reference.py``).

The standalone binding ``_honest_reader_impl`` takes the per-slot factors packed as arrays — the same inputs
the block solve assembles in place from its arena after the last refit — so the kernel is held to the
reference curve by curve, and its mode search, threading and analytic limits are held on the same fixture.
"""

import math

import numpy as np
import pytest

from rigel import _honest_reader_impl as HR
from rigel.calibration.simplex_logodds import _TILT_NODES
from rigel.config import CONSTANTS

from . import _honest_reader_reference as REF

COARSE_STEP = 0.2
WINDOW_SD = math.sqrt(2.0 * -math.log(np.finfo(np.float64).eps))
RULES = dict(coarse_step=COARSE_STEP, window_sd=WINDOW_SD, tilt_nodes=int(_TILT_NODES))
LAM = np.arange(-10.0, 10.0 + 1e-9, 0.2)
X = np.linspace(-14.0, 3.0, int(CONSTANTS.landscape.grid_points))


def _landscape(rng):
    """A bimodal landscape on the production grid, floored as the fit floors its density."""
    p = 0.8 * np.exp(-0.5 * ((X + 6.0) / 0.6) ** 2) + 0.2 * np.exp(-0.5 * ((X + 1.0) / 0.4) ** 2)
    p = p / p.sum() + 1e-12
    return np.log(p / p.sum())


def _level(rng, lo_centre, hi_centre):
    """A held level: a concave curve on its own grid with an origin."""
    m = int(rng.integers(5, 30))
    a = float(rng.uniform(lo_centre, hi_centre))
    g = np.linspace(a - rng.uniform(1.0, 4.0), a + rng.uniform(1.0, 4.0), m)
    w = float(rng.uniform(0.3, 2.0))
    vals = -0.5 * ((g - a) / w) ** 2 * float(rng.uniform(0.5, 3.0))
    return g, vals, float(np.exp(rng.uniform(-3.0, 3.0)))


def _comp(rng, live=True):
    if not live:
        return np.zeros(LAM.size)
    c, w = float(rng.uniform(-4.0, 4.0)), float(rng.uniform(0.5, 3.0))
    if rng.uniform() < 0.3:  # a cliff, the shape that moves the peak off the own term's mode
        return -np.abs(LAM - c) * float(rng.uniform(1.0, 6.0))
    return -0.5 * ((LAM - c) / w) ** 2 * float(rng.uniform(0.5, 4.0))


def _slots(rng, n=48):
    """A fixture covering every branch: pure-DNA, no admission, single-strand with and without reads, both-
    strand with and without reads, with and without each held level, live and dead compositions."""
    slots = []
    for i in range(n):
        kind = i % 6
        reads = rng.uniform() < 0.6
        u = float(rng.integers(0, 40)) if reads else 0.0
        v = float(rng.integers(0, 40)) if reads else 0.0
        if kind == 5 and reads:
            u, v = float(rng.integers(100, 400)), float(rng.integers(100, 400))
        eg = float(np.exp(rng.uniform(-3.0, 8.0)))
        er = 0.0 if kind == 0 else float(np.exp(rng.uniform(-3.0, 8.0)))
        q = (
            float(rng.uniform(0.55, 0.99))
            if rng.uniform() < 0.5
            else float(rng.uniform(0.01, 0.45))
        )
        both = kind in (3, 4, 5)
        admits = kind != 1
        slots.append(
            dict(
                u=u,
                v=v,
                eg=eg,
                er=er,
                q=q,
                both=both,
                admits=admits,
                lam=LAM,
                comp=_comp(rng, live=rng.uniform() < 0.8),
                dna=_level(rng, -9.0, 0.0) if rng.uniform() < 0.5 else None,
                pos=_level(rng, -2.0, 6.0) if (both or kind == 2) and rng.uniform() < 0.5 else None,
                neg=_level(rng, -2.0, 6.0) if both and rng.uniform() < 0.5 else None,
            )
        )
    return slots


def _pack(slots):
    n = len(slots)
    cols = {k: np.array([s[k] for s in slots], float) for k in ("u", "v", "eg", "er", "q")}
    flags = {k: np.array([s[k] for s in slots], np.uint8) for k in ("both", "admits")}
    comp = np.ascontiguousarray(np.stack([s["comp"] for s in slots]))
    levels = {}
    for k in ("dna", "pos", "neg"):
        width = max([s[k][0].size for s in slots if s[k] is not None] + [1])
        grid = np.full((n, width), np.nan)
        vals = np.full((n, width), np.nan)
        origin = np.ones(n)
        has = np.zeros(n, np.uint8)
        for i, s in enumerate(slots):
            if s[k] is not None:
                g, v_, o = s[k]
                grid[i, : g.size], vals[i, : g.size], origin[i], has[i] = g, v_, o, 1
        levels[k] = (grid, vals, origin, has)
    args = [
        cols["u"],
        cols["v"],
        cols["eg"],
        cols["er"],
        cols["q"],
        flags["both"],
        flags["admits"],
        LAM,
        comp,
    ]
    for k in ("dna", "pos", "neg"):
        args += list(levels[k])
    return args


def _native_modes(slots, logP, *, stride, n_threads=1, **rules):
    r = dict(RULES, **rules)
    mode = np.zeros(len(slots))
    mx = np.zeros(len(slots))
    HR.honest_modes(
        *_pack(slots),
        X,
        logP,
        mode,
        mx,
        r["coarse_step"],
        r["window_sd"],
        r["tilt_nodes"],
        True,
        True,
        1.0,
        False,
        False,
        int(stride),
        int(n_threads),
    )
    return mode, mx


def _native_curve(slots, i):
    out = np.zeros(X.size)
    HR.honest_logL(
        *_pack(slots),
        X,
        int(i),
        out,
        RULES["coarse_step"],
        RULES["window_sd"],
        RULES["tilt_nodes"],
        True,
        True,
        1.0,
        False,
        False,
    )
    return out


@pytest.fixture(scope="module")
def fixture():
    rng = np.random.default_rng(20261010)
    slots = _slots(rng)
    return slots, _landscape(rng)


def test_the_kernel_reproduces_the_reference_curve_at_every_slot(fixture):
    """Curve by curve, on every branch of the model: the evaluator is the specification to rounding."""
    slots, _ = fixture
    rho = np.exp(X)
    covered = set()
    for i, s in enumerate(slots):
        want = REF.log_L(s, rho, **RULES)
        got = _native_curve(slots, i)
        np.testing.assert_allclose(
            got, want, rtol=1e-9, atol=1e-7, err_msg=f"slot {i}: {s['u']}/{s['v']} both={s['both']}"
        )
        covered.add(
            (
                s["er"] <= 0.0 or not s["admits"],
                s["both"],
                s["u"] + s["v"] > 0,
                s["pos"] is not None or s["neg"] is not None,
            )
        )
    assert len(covered) >= 8, f"the fixture covers too few branches: {sorted(covered)}"


def test_the_kernel_reproduces_the_reference_mode_at_every_slot(fixture):
    slots, logP = fixture
    got, _ = _native_modes(slots, logP, stride=0)
    want = np.array([REF.mode(REF.log_L(s, np.exp(X), **RULES), X, logP) for s in slots])
    np.testing.assert_allclose(got, want, rtol=0.0, atol=1e-7)


def test_the_zero_read_search_finds_the_full_grid_mode(fixture):
    """The shipped stride (``CONSTANTS.calibration.reader_search_stride``) searches a slot with no read
    through the seeds and the windows; it must land on exactly the mode the full grid gives, at every slot."""
    slots, logP = fixture
    full, _ = _native_modes(slots, logP, stride=0)
    searched, _ = _native_modes(slots, logP, stride=CONSTANTS.calibration.reader_search_stride)
    assert any(s["u"] + s["v"] == 0 for s in slots)
    np.testing.assert_array_equal(searched, full)


def test_thread_count_changes_no_value(fixture):
    slots, logP = fixture
    one, mx1 = _native_modes(
        slots, logP, stride=CONSTANTS.calibration.reader_search_stride, n_threads=1
    )
    four, mx4 = _native_modes(
        slots, logP, stride=CONSTANTS.calibration.reader_search_stride, n_threads=4
    )
    np.testing.assert_array_equal(one, four)
    np.testing.assert_array_equal(mx1, mx4)


def _bare(**kw):
    s = dict(
        u=0.0,
        v=0.0,
        eg=1.0,
        er=0.0,
        q=0.9,
        both=False,
        admits=True,
        lam=LAM,
        comp=np.zeros(LAM.size),
        dna=None,
        pos=None,
        neg=None,
    )
    s.update(kw)
    return s


def test_flat_evidence_reads_the_priors_own_mode(fixture):
    """A slot with no read, a negligible opportunity and nothing delivered has a flat evidence, and its
    mode is the landscape's own refined peak."""
    _, logP = fixture
    slots = [_bare(eg=1e-12)]
    curve = _native_curve(slots, 0)
    assert np.all(np.abs(curve) < 1e-6)
    got, _ = _native_modes(slots, logP, stride=0)
    assert got[0] == pytest.approx(REF.mode(np.zeros(X.size), X, logP), abs=1e-12)


def test_pure_dna_is_the_poisson_closed_form_with_the_composition_at_its_gdna_end():
    """``Er = 0`` (and a slot admitting no RNA strand) is the pure-DNA limit: the two columns' Poisson
    log-likelihoods at ``ρEg/2`` each, plus the composition row's all-gDNA end when it is live, plus the held
    DNA level at ``ρ`` — in closed form, no integral."""
    u, v, eg = 37.0, 63.0, 250.0
    rho = np.exp(X)
    comp = -0.5 * ((LAM - 1.5) / 0.7) ** 2
    g = np.linspace(-9.0, -1.0, 11)
    vals = -0.5 * ((g + 4.0) / 0.8) ** 2
    for slot in (
        _bare(u=u, v=v, eg=eg, er=0.0, comp=comp, dna=(g, vals, 2.5)),
        _bare(u=u, v=v, eg=eg, er=40.0, admits=False, comp=comp, dna=(g, vals, 2.5)),
    ):
        want = (u + v) * np.log(rho * eg / 2.0) - rho * eg - REF.gammaln(u + 1) - REF.gammaln(v + 1)
        want = want + comp[-1] + np.interp(X - math.log(2.5), g, vals)
        np.testing.assert_allclose(_native_curve([slot], 0), want, rtol=1e-12, atol=1e-9)


def test_pure_dna_under_a_flat_prior_reads_the_count_over_the_opportunity():
    """With a flat prior the pure-DNA mode is ``log(n/Eg)`` exactly in the continuum; on the lattice the
    parabola through the three nearest nodes misses it by at most the square of the grid step."""
    n, eg = 100.0, 100.0
    got, _ = _native_modes([_bare(u=40.0, v=60.0, eg=eg)], np.zeros(X.size), stride=0)
    step = float(X[1] - X[0])
    assert abs(got[0] - math.log(n / eg)) < step**2
