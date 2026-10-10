"""The honest capture reader — the executable SPECIFICATION of ``native/honest_reader.h``.

This module is the reference implementation of the per-slot capture reader the block solve runs after its
last refit's final ψ (``native.solve_blocks``, ``reader=``), and it is the authority: the native evaluator is
required to reproduce every value of its log-evidence curve and every mode (``test_honest_reader.py``). It is
the archived NumPy prototype the reader was A/B'd with, lattices and all, written out once more so the
suite carries it.

THE MODEL. One slot of the chain holds ``u`` and ``v`` unspliced fragments on its two strand columns, a gDNA
opportunity ``Eg`` and an RNA opportunity ``Er``, and the protocol's column probability ``q`` for its RNA
strand. gDNA is unstranded and lands on either column with probability ½; RNA of amount ``r`` (fragments per
unit RNA opportunity, times ``Er``) lands on the RNA strand's column with probability ``q``. The log evidence
for a gDNA density ``ρ`` integrates the RNA amount out under the factors the sweep delivered to the slot::

    log L(x) = log ∫ Pois(u; ρEg/2 + q·r) · Pois(v; ρEg/2 + (1−q)·r) · C(log(ρEg/r))
                   · H_pos(r·f_pos / Er) · H_neg(r·f_neg / Er) · r^(−1/2) dr  +  D(log ρ),     ρ = e^x,

``C`` the delivered composition row on the λ lattice (the log-odds of gDNA to RNA mass), ``D`` the held DNA
level at ``ρ`` (a neighbour's statement of this slot's gDNA density, at its own origin), ``H_pos`` and
``H_neg`` the held RNA levels at their own amounts, every delivered curve held constant beyond its ends, and
``r^(−1/2)`` the amount's reference measure. A both-strand slot (``both``) integrates the share ``s`` of
its RNA on the positive strand over the arcsine continuum plus the pure atoms its witness bits admit — a
strand with no held RNA level admits the pure-other-strand atom — normalised by ``−log 3``; a slot with no
RNA opportunity (``Er = 0`` or no strand admitted) is the pure-DNA limit with the composition read at its
all-gDNA end. The weight the solve publishes is ``exp(x*)`` relative to the largest, ``x*`` the mode of
``log L(x) + logP(x)`` under the fitted landscape, the lattice argmax refined by the parabola through its two
neighbours (``mode``). No reference density, no background, no detector, no clip.

THE LATTICES, each a derivation already in the tree. The RNA amount ``y = log r``: the union of the count
solver's 0.2-nat lattice on ``[y_hi − 20, y_hi]``, ``y_hi = log(n + 12√n + 12)``, with a window of
``⌈window_sd⌉`` nodes either side at step ``h = min(0.2, 1/√(n+1))`` around the own term's mode ``y*`` (the
right-hand zero of ``r·d/dr log f``, by bisection, centred once per slot at the slot's own ``q``) and the same
window around the full integrand's maximum on the first lattice; trapezoid weights from the spacings, every
node clipped to ``[y_lo − 5, y_hi + 1]``. The strand share: ψ's 24 midpoint nodes in ``φ ∈ (0, π/2)`` with
``s = cos²φ`` and the arcsine measure ``2dφ/π``, plus a window at step ``sd_s`` around the observed share
when the columns resolve it, plus a window around the density nodes' maxima clustered at ``sd_s``, the
endpoints included; trapezoid in ``φ``. ``window_sd = √(2T)``, ``T = −log ε₆₄``, is ψ's window rule.
"""

import math

import numpy as np
from scipy.special import gammaln, logsumexp


def rna_mode_y(a, u, v, q, y_lo, y_hi):
    """Per density node, the RNA amount (log) maximising the own term with its ``r^(−1/2)`` measure: the
    unique right-hand zero of ``r·d/dr log f``, by bisection in ``y = log r``."""
    lo = np.full_like(a, y_lo)
    hi = np.full_like(a, y_hi)
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        r = np.exp(mid)
        g = -r + 0.5
        if u > 0:
            g = g + u * q * r / (a / 2 + q * r)
        if v > 0:
            g = g + v * (1 - q) * r / (a / 2 + (1 - q) * r)
        up = g > 0
        lo = np.where(up, mid, lo)
        hi = np.where(up, hi, mid)
    return 0.5 * (lo + hi)


def log_L(slot, rho, *, coarse_step, window_sd, tilt_nodes):
    """The log evidence at every density in ``rho`` (1-D), under the lattices above.

    ``slot`` is a dict: ``u, v, eg, er, q`` (floats), ``both, admits`` (bools), ``lam`` and ``comp`` (the λ
    lattice and the composition row on it), and ``dna, pos, neg`` — each ``None`` or ``(grid, vals, origin)``,
    a held level on its own grid (nats of log density for ``dna``, nats of log amount for the RNA levels)
    with its origin."""
    u, v, eg, er, q = slot["u"], slot["v"], slot["eg"], slot["er"], slot["q"]
    lam, comp = np.asarray(slot["lam"], float), np.asarray(slot["comp"], float)
    rho = np.asarray(rho, float)
    dna = rho * eg
    n = u + v
    sd = 1.0 / math.sqrt(n + 1.0)
    h = min(coarse_step, sd)
    lg = gammaln(u + 1) + gammaln(v + 1)
    comp_live = np.ptp(comp) > 0
    dna_lv, pos_lv, neg_lv = slot["dna"], slot["pos"], slot["neg"]

    def dna_term(rho_):
        if dna_lv is None:
            return 0.0
        g, vals, origin = dna_lv
        return np.interp(np.log(rho_) - math.log(origin), g, vals)

    if er <= 0.0 or not slot["admits"]:
        with np.errstate(divide="ignore"):
            own = (
                (u * np.log(dna / 2) if u > 0 else 0.0)
                + (v * np.log(dna / 2) if v > 0 else 0.0)
                - dna
                - lg
            )
        return own + (comp[-1] if comp_live else 0.0) + dna_term(rho)

    y_hi = math.log(n + 12.0 * math.sqrt(n) + 12.0)
    y_lo = y_hi - 20.0
    coarse = np.arange(y_lo, y_hi + coarse_step, coarse_step)
    log_dna = np.log(dna)
    offs = np.arange(-math.ceil(window_sd), math.ceil(window_sd) + 1) * h
    y_star = rna_mode_y(dna, u, v, q, y_lo - 5.0, y_hi + 1.0) if h < coarse_step else None

    def lattice(centres):
        parts = [np.broadcast_to(coarse, (dna.size, coarse.size))]
        for c in centres:
            parts.append(c[:, None] + offs[None, :])
        y = np.sort(np.concatenate(parts, axis=1), axis=1)
        y = np.clip(y, y_lo - 5.0, y_hi + 1.0)
        dy = np.diff(y, axis=1)
        logw = np.empty_like(y)
        with np.errstate(divide="ignore"):
            logw[:, 1:-1] = np.log(0.5 * (dy[:, :-1] + dy[:, 1:]))
            logw[:, 0] = np.log(0.5 * dy[:, 0])
            logw[:, -1] = np.log(0.5 * dy[:, -1])
        return y, logw

    a2 = dna[:, None] / 2.0

    def log_f(y, qq, fpos, fneg):
        r = np.exp(y)
        out = -dna[:, None] - r - lg + 0.5 * y
        if u > 0:
            out = out + u * np.log(a2 + qq * r)
        if v > 0:
            out = out + v * np.log(a2 + (1.0 - qq) * r)
        if comp_live:
            out = out + np.interp(log_dna[:, None] - y, lam, comp)
        for lv, frac in ((pos_lv, fpos), (neg_lv, fneg)):
            if lv is None:
                continue
            if frac == 0.0:
                out = out + lv[1][0]
            else:
                g, vals, origin = lv
                out = out + np.interp(y, g + math.log(origin) + math.log(er) - math.log(frac), vals)
        return out

    def inner(qq, fpos, fneg):
        centres = [] if y_star is None else [y_star]
        y, logw = lattice(centres)
        f = log_f(y, qq, fpos, fneg)
        if y_star is not None:
            peak = y[np.arange(y.shape[0]), np.argmax(f, axis=1)]
            y, logw = lattice(centres + [peak])
            f = log_f(y, qq, fpos, fneg)
        return logsumexp(f + logw, axis=1)

    if not slot["both"]:
        out = inner(q, 1.0, 0.0)
    else:
        m = int(tilt_nodes)
        kwin = np.arange(-math.ceil(window_sd), math.ceil(window_sd) + 1)
        sd_s = 1.0 / (math.sqrt(n + 1.0) * max(abs(2 * q - 1), 1.0 / math.sqrt(n + 1.0)))

        def s_window(centre):
            s_win = centre + kwin * sd_s
            return s_win[(s_win > 0.0) & (s_win < 1.0)]

        def tilt(phi):
            phi = np.unique(np.concatenate([[0.0], phi, [math.pi / 2]]))
            parts = []
            for p in phi:
                s_ = math.cos(p) ** 2
                parts.append(inner(q * s_ + (1.0 - q) * (1.0 - s_), s_, 1.0 - s_))
            K = np.stack(parts, axis=1)
            dphi = np.diff(phi)
            wphi = np.empty_like(phi)
            wphi[1:-1] = 0.5 * (dphi[:-1] + dphi[1:])
            wphi[0], wphi[-1] = 0.5 * dphi[0], 0.5 * dphi[-1]
            with np.errstate(divide="ignore"):
                return (
                    phi,
                    K,
                    logsumexp(K + np.log(wphi)[None, :], axis=1) + math.log(2.0 / math.pi),
                )

        phi = (np.arange(m) + 0.5) * (math.pi / 2) / m
        if n > 0 and abs(2 * q - 1) * math.sqrt(n + 1) > 1.0:
            s_star = (u / n - (1 - q)) / (2 * q - 1)
            phi = np.concatenate([phi, np.arccos(np.sqrt(s_window(s_star)))])
        phi, K, mixed = tilt(phi)
        peaks = np.cos(phi[np.argmax(K, axis=1)]) ** 2
        centres = np.unique(np.round(peaks / sd_s)) * sd_s
        extra = (
            np.unique(np.concatenate([s_window(float(c)) for c in centres]))
            if centres.size
            else np.zeros(0)
        )
        if extra.size:
            phi, K, mixed = tilt(np.concatenate([phi, np.arccos(np.sqrt(extra))]))
        hyps = [mixed]
        if neg_lv is None:
            hyps.append(inner(q, 1.0, 0.0))
        if pos_lv is None:
            hyps.append(inner(1.0 - q, 0.0, 1.0))
        out = logsumexp(np.stack(hyps, axis=1), axis=1) - math.log(3.0)
    return out + dna_term(rho)


def mode(logL, x, logP):
    """The posterior mode in log density: the lattice argmax of ``logL + logP`` refined by the parabola
    through its two neighbours (at an end node, the node itself)."""
    lp = np.asarray(logL, float) + np.asarray(logP, float)
    k = int(np.argmax(lp))
    if 0 < k < lp.size - 1:
        a, b, c = lp[k - 1], lp[k], lp[k + 1]
        den = a - 2.0 * b + c
        off = 0.5 * (a - c) / den if den < 0 else 0.0
        return float(x[k] + off * (x[k + 1] - x[k]))
    return float(x[k])
