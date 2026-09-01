"""PROTOTYPE (outside src/) — the arcsine MAGNITUDE coordinate for the simplex solver.

    ⚠ **A SANDBOX FILE, and it is PRESERVED EVIDENCE rather than a proposal.** The A/B it exists for
    RAN (2026-09-01) and the coordinate did NOT ship — `ISSUES: arcsine-magnitude-coordinate` carries
    the verdict, and the derivation doc beside this file carries the argument. It is kept because
    ① its self-test gates are what REFUTED the issue's headline premise and re-deriving them is a
    day's work, and ② if the owner rules to try rung B (phi-native message delivery), this is where
    that starts. ⛔ Nothing may cite it: not the source, not a test, not a permanent doc.

`ISSUES: arcsine-magnitude-coordinate`, step 3. The derivation is the doc beside this file.

The whole change, per the derivation:
  §1  the base measure conditioned on the slot total is `df/(f(1-f)) = dlam` — flat in lambda. Under
      `f = sin^2(phi)` the conversion `log(dlam/dphi) = -1/2 log f - 1/2 log(1-f) + log 4` is EXACTLY
      minus the two written Jeffreys reference halves, so under phi BOTH vanish: no Jacobian AND no
      reference is written.
  §2  flat-in-phi pushes forward to Beta(1/2,1/2) — the Berger-Bernardo f_g marginal. The reference
      becomes the grid measure itself. The theta axis is untouched (its own cancellation is
      independent: the lambda<->phi Jacobian is theta-free).
  §4  MIDPOINT nodes, uniform in phi on [0, pi/2]: `phi_k = (k+1/2)*(pi/2)/K`. Uniform weights are
      then the reference exactly, every node is strictly interior (so every log-rate term is finite —
      no endpoint quadrature, no floor, no epsilon), and the OUTER BINS reach the vertices: a
      posterior crossing 1/2 at the last bin's far edge reads out f = 1 exactly.

Three module-level patches implement it (`install`, an idempotent context manager):
  A `_logodds_grid`  -> the phi grid, returned as `(logit(f), f)` so every row-builder (the fitted
                        arms, the intron-factory NB rows, the certified-flux `lam_rows`) evaluates
                        at the new f values, and every lambda-coordinate consumer stays consistent.
  B `_gdna_arm` / `_rna_arm` -> drop the reference half, keep the fitted logP (§1).
  C `_posterior_median_fg`   -> interpolate the histogram CDF in PHI (the coordinate whose measure is
                        uniform), not in lambda. Bin edges are phi midpoints; the outer bins are
                        [0, dphi] and [pi/2 - dphi, pi/2], inside the domain by construction.

⛔ The NOOP arm (`install(arm="noop")`) patches with the ORIGINAL implementations and MUST be
byte-identical to unpatched production. That is this file's falsification: a harness that changes an
answer by existing cannot price anything.
"""

from __future__ import annotations

import contextlib
import sys as _sys

import numpy as np
from scipy.special import expit, logit

from rigel.calibration import region_init as _ri_mod
from rigel.calibration import simplex_logodds as S
from rigel.calibration import sweep as _sweep_mod

_HALF_PI = 0.5 * np.pi


def phi_grid(n_grid: int):
    """The MIDPOINT uniform-phi lattice (§4): `phi_k = (k+1/2)*dphi`, `dphi = (pi/2)/K`, ascending.

    Every node is strictly interior, so `f in (0,1)` and every log-rate term is finite — the
    endpoint singularity the lambda grid handled with a bracket, and that a node-on-endpoint grid
    would handle with quadrature, does not arise. Returns `(phi, f)`.
    """
    K = int(n_grid)
    dphi = _HALF_PI / K
    phi = (np.arange(K, dtype=np.float64) + 0.5) * dphi
    return phi, np.sin(phi) ** 2


def _phi_logodds_grid(n_grid: int, L: float = S._DEFAULT_L):
    """Patch A. The signature `_logodds_grid` has, returning the SAME pair meaning: `(lambda, f)`.

    `L` is accepted and IGNORED — a bounded coordinate has no bracket, which is the point
    (`ISSUES: psi-lambda-bracket-unshipped` dissolves). Callers that grid rows off `f` get the phi
    nodes; callers that need the lambda coordinate of those nodes get it exactly.
    """
    phi, fg = phi_grid(n_grid)
    return logit(fg), fg


def _phi_gdna_arm(lam, global_logprior):
    """Patch B (gDNA). No reference term: under phi the measure IS the reference (§1-§2). A fitted
    `logP_g`, pre-evaluated on this f grid, is the arm's whole content when present."""
    if global_logprior is None:
        return np.zeros((1, np.asarray(lam).shape[0]), np.float64)
    return np.asarray(global_logprior, np.float64)


def _phi_rna_arm(lam, rna_logprior=None):
    """Patch B (RNA-total). The exact mirror of :func:`_phi_gdna_arm`."""
    if rna_logprior is None:
        return np.zeros((1, np.asarray(lam).shape[0]), np.float64)
    return np.asarray(rna_logprior, np.float64)


def _phi_posterior_median_fg(post, lam, fg):
    """Patch C. The posterior 1/2-quantile read off the CDF, interpolated in PHI.

    Same construction as production's, and for the same reason (`_posterior_median_fg`'s own
    argument): the interpolation must happen in the coordinate whose measure is uniform, or a bin's
    midpoint is not the image of its midpoint and a one-hot posterior comes back biased toward 1/2.
    Under phi that coordinate is phi. The median is transform-equivariant, so `sin^2` of the phi
    quantile IS the f quantile, exactly.

    ⭐ The outer bins differ from production's mirrored half-bins: the phi domain is BOUNDED, so the
    outer edges are the domain ends 0 and pi/2 — no mirroring outside the domain, and the vertices
    are reachable read-outs (a CDF crossing at the last bin's far edge returns f = 1 exactly).
    """
    p = np.asarray(post, np.float64)
    f = np.asarray(fg, np.float64)
    x = np.arcsin(np.sqrt(np.clip(f, 0.0, 1.0)))  # phi of each node
    edges = np.empty(x.shape[0] + 1, np.float64)
    edges[1:-1] = 0.5 * (x[:-1] + x[1:])
    edges[0] = 0.0
    edges[-1] = _HALF_PI
    tot = p.sum(axis=1, keepdims=True)
    cdf = np.concatenate(
        [np.zeros((p.shape[0], 1)), np.cumsum(p / np.where(tot > 0.0, tot, 1.0), axis=1)], axis=1
    )
    k = np.clip((cdf < 0.5).sum(axis=1), 1, x.shape[0])
    lo = np.take_along_axis(cdf, (k - 1)[:, None], 1)[:, 0]
    hi = np.take_along_axis(cdf, k[:, None], 1)[:, 0]
    span = hi - lo
    t = np.where(span > 0.0, (0.5 - lo) / np.where(span > 0.0, span, 1.0), 0.5)
    return np.sin(edges[k - 1] + t * (edges[k] - edges[k - 1])) ** 2


# `_logodds_grid` is imported BY NAME into these modules, so each one's global must be rebound too.
# ⚠ `rigel.calibration.calibrate` resolves to the FUNCTION through the package's re-export, so the
#   module object is taken from sys.modules by its dotted name.
_calibrate_module = _sys.modules["rigel.calibration.calibrate"]
_GRID_IMPORTERS = (S, _sweep_mod, _ri_mod, _calibrate_module)
# (`rna_anchor` imports it inside the function body, so it picks the patched module attribute up.)
_ORIG = {
    "grid": S._logodds_grid,
    "gdna_arm": S._gdna_arm,
    "rna_arm": S._rna_arm,
    "median": S._posterior_median_fg,
}


@contextlib.contextmanager
def install(arm: str = "arcsine"):
    """Patch the three surfaces for the duration of the block.

    ``arm="noop"`` reinstalls the ORIGINALS through the identical mechanism — the harness's own
    falsification, which must be byte-identical to unpatched production.
    """
    if arm not in ("arcsine", "noop"):
        raise ValueError(f"unknown arm {arm!r} — expected 'arcsine' or 'noop'")
    new = (
        {
            "grid": _phi_logodds_grid,
            "gdna_arm": _phi_gdna_arm,
            "rna_arm": _phi_rna_arm,
            "median": _phi_posterior_median_fg,
        }
        if arm == "arcsine"
        else dict(_ORIG)
    )
    saved = [(m, m._logodds_grid) for m in _GRID_IMPORTERS]
    saved_s = (S._gdna_arm, S._rna_arm, S._posterior_median_fg)
    try:
        for m in _GRID_IMPORTERS:
            m._logodds_grid = new["grid"]
        S._gdna_arm = new["gdna_arm"]
        S._rna_arm = new["rna_arm"]
        S._posterior_median_fg = new["median"]
        yield
    finally:
        for m, fn in saved:
            m._logodds_grid = fn
        S._gdna_arm, S._rna_arm, S._posterior_median_fg = saved_s


def self_test() -> int:
    """Gates on the derivation's own claims — each perturbed, per `falsification_needs_perturbation`."""
    ok = True

    def check(name, cond):
        nonlocal ok
        print(f"   {'✔' if cond else '✘'} {name}")
        ok = ok and bool(cond)

    # ① the grid: midpoint nodes, strictly interior, uniform in phi, reaching far past the lambda bracket.
    phi, f = phi_grid(256)
    d = np.diff(phi)
    check("phi nodes uniformly spaced", np.allclose(d, d[0], rtol=0, atol=1e-15))
    check("phi nodes strictly interior to [0, pi/2]", bool(phi[0] > 0 and phi[-1] < _HALF_PI))
    check("f strictly interior to (0,1)", bool(f[0] > 0.0 and f[-1] < 1.0))
    check(f"extreme f ({f[0]:.3e}) closer to the vertex than sigma(-10) ({expit(-10.0):.3e})",
          bool(f[0] < expit(-10.0)))
    check("f symmetric about 1/2 (the composition has no preferred component)",
          np.allclose(f + f[::-1], 1.0, atol=1e-15))

    # ② the measure: uniform weights in phi ARE Beta(1/2,1/2) in f. Perturbation: uniform-in-LAMBDA
    #    weights on the same nodes must FAIL the same test.
    from scipy.stats import beta as _beta
    edges = np.concatenate([[0.0], 0.5 * (phi[:-1] + phi[1:]), [_HALF_PI]])
    f_edges = np.sin(edges) ** 2
    exact = np.diff(_beta.cdf(f_edges, 0.5, 0.5))
    unif = np.full(phi.shape[0], 1.0 / phi.shape[0])
    check(f"uniform-phi weights == Beta(1/2,1/2) bin mass (max dev {np.abs(unif-exact).max():.2e})",
          np.allclose(unif, exact, atol=1e-12))
    lam_nodes = logit(f)
    w_lam = np.diff(np.concatenate([[lam_nodes[0]], 0.5 * (lam_nodes[:-1] + lam_nodes[1:]),
                                    [lam_nodes[-1]]]))
    w_lam = w_lam / w_lam.sum()
    check("PERTURBED (uniform-in-lambda weights on the same nodes) does NOT match Beta(1/2,1/2)",
          not np.allclose(w_lam, exact, atol=1e-12))

    # ③ the reference cancellation, numerically: psi_lambda(with reference) on lambda-uniform nodes and
    #    psi_phi(no reference) on phi-uniform nodes must integrate the SAME model. Test on the pure
    #    reference: both must return the Beta(1/2,1/2) median = 1/2, and both must return the same
    #    median for a tilted likelihood.
    K = 20001
    lam_u = np.linspace(-40.0, 40.0, K)
    f_lam = expit(lam_u)
    phi_u, f_phi = phi_grid(K)
    for name, loglik in (("flat", lambda z: np.zeros_like(z)),
                         ("Binom(7,10)", lambda z: 7 * np.log(z) + 3 * np.log1p(-z)),
                         ("Binom(60,60)", lambda z: 60 * np.log(z))):
        p_l = np.exp(loglik(f_lam) + 0.5 * np.log(f_lam) + 0.5 * np.log1p(-f_lam))
        p_l /= p_l.sum()
        m_l = _ORIG["median"](p_l[None, :], lam_u, f_lam)[0]
        p_p = np.exp(loglik(f_phi) - np.max(loglik(f_phi)))
        p_p /= p_p.sum()
        m_p = _phi_posterior_median_fg(p_p[None, :], logit(f_phi), f_phi)[0]
        check(f"reference cancellation, {name}: lambda-with-ref {m_l:.9f} == phi-no-ref {m_p:.9f}",
              abs(m_l - m_p) < 2e-6)
    # PERTURBATION: keeping the reference under phi must BREAK the agreement. ⚠ It has to be tested on
    # an ASYMMETRIC likelihood — under phi the retained reference gives Beta(1,1), whose median is 1/2
    # like Beta(1/2,1/2)'s, so the symmetric case cannot detect the very error it is meant to catch.
    _ll = 7 * np.log(f_phi) + 3 * np.log1p(-f_phi)
    p_bad = np.exp(_ll + 0.5 * np.log(f_phi) + 0.5 * np.log1p(-f_phi) - np.max(_ll))
    p_bad /= p_bad.sum()
    m_bad = _phi_posterior_median_fg(p_bad[None, :], logit(f_phi), f_phi)[0]
    check(f"PERTURBED (reference kept under phi) shifts Binom(7,10): {m_bad:.6f} != 0.693176",
          abs(m_bad - 0.693176) > 1e-3)

    # ④ ⛔ THE VERTEX, MEASURED RATHER THAN ASSERTED — and it REFUTES the reachability claim the issue
    #    and this file's first draft both made. Two facts, each gated:
    #    (a) the posterior MEDIAN is bounded away from the vertex by half an outer bin in ANY
    #        coordinate. For the last bin, the CDF crossing fraction is t = 1 - 0.5/m <= 1/2 for any
    #        bin mass m <= 1, so the read-out cannot pass the bin's MIDPOINT however concentrated the
    #        posterior is. The vertex is unreachable for a median read-out, not for lambda.
    #    (b) in f-space the lambda grid is FINER near the vertices than uniform-phi, not coarser:
    #        sigma(+-L) already sits 4.5e-5 from the vertex.
    phi_s, f_s = phi_grid(60)
    lam_s, fl_s = S._logodds_grid(60, 10.0)
    one_hot = np.zeros((1, 60))
    one_hot[0, -1] = 1.0
    hi_phi = _phi_posterior_median_fg(one_hot, logit(f_s), f_s)[0]
    hi_lam = _ORIG["median"](one_hot, lam_s, fl_s)[0]
    dphi = _HALF_PI / 60
    check(f"(a) phi's max possible read-out is 1 - sin^2(dphi/2) = {1-np.sin(dphi/2)**2:.9f}, "
          f"attained by a one-hot outer bin ({hi_phi:.9f})",
          abs(hi_phi - (1 - np.sin(dphi / 2) ** 2)) < 1e-12)
    tilted = np.zeros((1, 60))
    tilted[0, -1], tilted[0, 0] = 0.999, 0.001
    check("(a) and NO posterior mass distribution passes it (t = 1 - 0.5/m <= 1/2)",
          _phi_posterior_median_fg(tilted, logit(f_s), f_s)[0] <= hi_phi + 1e-15)
    check(f"(b) ⛔ REFUTES 'phi reaches the vertex, lambda cannot': lambda reads {hi_lam:.9f}, "
          f"CLOSER to 1 than phi's {hi_phi:.9f} at K=60",
          bool(hi_lam > hi_phi))
    check(f"(b) the lambda grid is finer in f near the vertex: top f-gap lambda "
          f"{fl_s[-1]-fl_s[-2]:.3e} vs phi {f_s[-1]-f_s[-2]:.3e}",
          bool((fl_s[-1] - fl_s[-2]) < (f_s[-1] - f_s[-2])))
    mid_l = np.abs(fl_s - 0.5).argmin()
    mid_p = np.abs(f_s - 0.5).argmin()
    check(f"⭐ but phi is finer MID-simplex — the reallocation: f-gap at 1/2 lambda "
          f"{fl_s[mid_l+1]-fl_s[mid_l]:.4f} vs phi {f_s[mid_p+1]-f_s[mid_p]:.4f}",
          bool((f_s[mid_p + 1] - f_s[mid_p]) < (fl_s[mid_l + 1] - fl_s[mid_l])))
    check(f"⭐ and the MEASURED shortfall (0.15-0.29) is ~4 decades past either grid's vertex "
          f"resolution ({1-fl_s[-1]:.1e} / {1-hi_phi:.1e}) — neither grid is what pins it",
          bool((1 - hi_phi) < 1e-3 and (1 - fl_s[-1]) < 1e-3))

    # ⑤ the patch mechanism: `install` restores exactly, and "noop" is the identity.
    before = (S._logodds_grid, S._gdna_arm, S._rna_arm, S._posterior_median_fg,
              _sweep_mod._logodds_grid, _ri_mod._logodds_grid, _calibrate_module._logodds_grid)
    with install("arcsine"):
        inside = S._logodds_grid is _phi_logodds_grid and _sweep_mod._logodds_grid is _phi_logodds_grid
    after = (S._logodds_grid, S._gdna_arm, S._rna_arm, S._posterior_median_fg,
             _sweep_mod._logodds_grid, _ri_mod._logodds_grid, _calibrate_module._logodds_grid)
    check("install('arcsine') patches every importer", inside)
    check("install restores every surface exactly", before == after)
    with install("noop"):
        noop_ok = (S._logodds_grid is _ORIG["grid"] and S._gdna_arm is _ORIG["gdna_arm"]
                   and S._posterior_median_fg is _ORIG["median"])
    check("install('noop') installs the ORIGINALS (the harness's own falsification)", noop_ok)
    try:
        with install("bogus"):
            pass
        check("an unknown arm RAISES", False)
    except ValueError:
        check("an unknown arm RAISES", True)

    print(f"\n{'✔ all gates pass' if ok else '✘ FAILURES'}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(self_test())
