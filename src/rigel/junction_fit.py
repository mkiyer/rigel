"""
rigel.junction_fit — the genuine junctions' strand rate κ, fitted over the per-junction strand table.

In a stranded library an annotated splice junction is one of three classes: genuine RNA, whose reads fall on the
minority strand at rate κ; a splice artifact — gDNA the aligner wrote as spliced — at ½; or reversed, at 1 − κ. The
fit maximises that three-class binomial mixture over κ and the class shares together, then reads κ as the shipped
``Beta(1, 1)`` posterior mean over the reads the maximum credits to RNA (:func:`genuine_sense_fraction`).

Every numerical setting of the search is ``rigel.config.CONSTANTS.junction_fit``; this module holds no constant of
its own. The one consumer is :attr:`rigel.strand_model.StrandModel.genuine_n_same`. Gates:
``tests/test_strand_model.py``'s ``TestGenuineKappa``.
"""

import numpy as np

from .config import CONSTANTS, JunctionFitConstants


def genuine_sense_fraction(
    n_sense, n_total, constants: JunctionFitConstants = CONSTANTS.junction_fit
) -> float:
    """The genuine junctions' sense fraction κ: the shipped ``Beta(1, 1)`` posterior mean, ``(minority + 1)/(reads +
    2)``, over the reads the three-class mixture credits to RNA at its maximum (:func:`_genuine_mixture`). Where every
    junction is genuine that is exactly ``(n_same + 1)/(n_obs + 2)``; it is never 0. The table is oriented so the
    genuine class's minority rate lies below ½."""
    s = np.asarray(n_sense, dtype=np.float64)
    n = np.asarray(n_total, dtype=np.float64)
    keep = n > 0
    s, n = s[keep], n[keep]
    flip = s.sum() > 0.5 * n.sum()
    k = n - s if flip else s
    kn, cnt = np.unique(np.stack([k, n], axis=1), axis=0, return_counts=True)
    _, _, minority, reads, _ = _genuine_mixture(
        kn[:, 0], kn[:, 1], cnt.astype(np.float64), constants
    )
    kappa = (minority + 1.0) / (reads + 2.0)
    return 1.0 - kappa if flip else kappa


def _genuine_mixture(
    k: np.ndarray,
    n: np.ndarray,
    c: np.ndarray,
    constants: JunctionFitConstants = CONSTANTS.junction_fit,
) -> tuple[float, np.ndarray, float, float, float]:
    """The three-class binomial mixture's maximum likelihood. Junction group g (``c_g`` junctions with ``k_g`` of
    ``n_g`` reads on the minority side) is genuine RNA (rate κ), a splice artifact — misaligned gDNA (½) — or reversed
    (1 − κ). The class shares are solved exactly at each κ (:func:`_class_shares`); κ maximises that profile, by a scan
    over ``u = log(κ/(½ − κ))`` and Brent's method in its best cell — bounded work, where EM crawls without end as κ
    nears ½ (the classes merge) and can stop at a saddle. Returns ``(κ, shares, minority, reads, log-likelihood)``,
    ``minority`` and ``reads`` being the RNA classes' expected counts at the maximum."""
    from scipy.optimize import minimize_scalar

    half = n * np.log(0.5)

    def profile(u):
        kappa = 0.5 / (1.0 + float(np.exp(-u)))
        lp = np.stack(
            [
                k * np.log(kappa) + (n - k) * np.log1p(-kappa),
                half,
                k * np.log1p(-kappa) + (n - k) * np.log(kappa),
            ]
        )
        top = lp.max(axis=0)
        f = np.exp(lp - top)
        value, shares = _class_shares(f, c, constants)
        return value + float(c @ top), shares, f

    floor = constants.kappa_floor
    end = float(np.log((0.5 - floor) / floor))
    grid = np.linspace(-end, end, constants.kappa_scan_points)
    scanned = [profile(u)[0] for u in grid]
    i = int(np.argmax(scanned))
    found = minimize_scalar(
        lambda u: -profile(u)[0],
        bounds=(grid[max(i - 1, 0)], grid[min(i + 1, grid.size - 1)]),
        method="bounded",
        options={"xatol": constants.kappa_tolerance},
    )
    u = float(found.x) if -found.fun >= scanned[i] else float(grid[i])
    loglik, shares, f = profile(u)
    with np.errstate(divide="ignore", invalid="ignore"):
        resp = shares[:, None] * f / (shares @ f)
    minority = float(c @ (resp[0] * k + resp[2] * (n - k)))
    reads = float(c @ ((resp[0] + resp[2]) * n))
    return 0.5 / (1.0 + float(np.exp(-u))), shares, minority, reads, loglik


def _class_shares(
    f: np.ndarray, c: np.ndarray, constants: JunctionFitConstants = CONSTANTS.junction_fit
) -> tuple[float, np.ndarray]:
    """``max_w Σ_g c_g·log(Σ_k w_k f_kg)`` over the 3-simplex, and the maximising shares. The objective is concave in
    ``w``, so the KKT conditions are necessary and sufficient: at the maximum ``Σ_g c_g f_kg/(w·f)_g`` equals ``Σc`` for
    every share in use and is at most ``Σc`` for every share at 0. A vertex, then an edge, that meets them is the
    maximum; only an interior maximum needs the log-barrier ascent. ``f`` is (3, G), each column scaled by its own
    constant (which only shifts the objective). A candidate on which some junction has no probability scores −∞, so the
    arithmetic warnings it raises are silenced, never acted on."""
    with np.errstate(divide="ignore", invalid="ignore", over="ignore", under="ignore"):
        return _class_shares_unguarded(f, c, constants)


def _class_shares_unguarded(
    f: np.ndarray, c: np.ndarray, constants: JunctionFitConstants
) -> tuple[float, np.ndarray]:
    from scipy.optimize import brentq

    def value(w):
        return float(c @ np.log(w @ f))

    total = float(c.sum())
    slack = constants.kkt_slack * total

    def satisfies(w, unused):
        mult = f @ (c / (w @ f))
        return bool(np.all(np.isfinite(mult[unused])) and np.all(mult[unused] <= total + slack))

    candidates = []
    for k in range(3):
        vertex = np.eye(3)[k]
        if satisfies(vertex, [j for j in range(3) if j != k]):
            return value(vertex), vertex
        candidates.append(vertex)
    margin = constants.edge_margin
    for a_, b_, third in ((0, 1, 2), (0, 2, 1), (1, 2, 0)):  # an edge t·e_a + (1 − t)·e_b
        dd = f[a_] - f[b_]

        def edge_slope(t, dd=dd, a_=a_, b_=b_):
            return float(c @ (dd / (t * f[a_] + (1.0 - t) * f[b_])))

        if edge_slope(margin) > 0.0 > edge_slope(1.0 - margin):
            t = brentq(edge_slope, margin, 1.0 - margin, xtol=margin)
            e = np.zeros(3)
            e[a_], e[b_] = t, 1.0 - t
            if satisfies(e, [third]):
                return value(e), e
            candidates.append(e)
    d = f[:2] - f[2]  # the interior: free coordinates (w0, w1), w2 = 1 − w0 − w1
    w = np.full(3, 1.0 / 3.0)
    mu = constants.barrier_start * total
    while True:

        def phi(v, mu=mu):
            return value(v) + mu * float(np.log(v).sum())

        for _ in range(
            constants.newton_steps
        ):  # Newton on the barrier objective, kept strictly inside
            m = w @ f
            q = d / m
            g = q @ c + mu * (1.0 / w[:2] - 1.0 / w[2])
            h = -(q * c) @ q.T - mu * (np.diag(1.0 / w[:2] ** 2) + 1.0 / w[2] ** 2)
            step = -np.linalg.solve(h, g)
            dw = np.array([step[0], step[1], -step[0] - step[1]])
            shrink = dw < 0.0
            t = (
                min(1.0, constants.boundary_fraction * float(np.min(-w[shrink] / dw[shrink])))
                if shrink.any()
                else 1.0
            )
            now = phi(w)
            while t > constants.line_search_floor and phi(w + t * dw) < now:
                t *= constants.line_search_shrink
            if t <= constants.line_search_floor:
                break
            w = w + t * dw
            if abs(float(g @ step)) * t < constants.newton_tolerance * max(1.0, abs(now)):
                break
        if mu <= constants.barrier_floor:
            break
        mu = max(mu * constants.barrier_shrink, constants.barrier_floor)
    candidates.append(w)
    best = max(candidates, key=value)
    return value(best), best
