"""How much of the gDNA/RNA information does the accumulator's stored pair actually keep?

       Log:

THE QUESTION
    Each object stores an integer ``count`` and a fixed-point density ``Sum 1/placements``. Whether
    that pair is the right one has never been derived -- /TRAPS: pure-and-length-censored are open, and
    §1 fact 5's ranking was measured at ONE fragment-length setting. The accumulator's deposit weight
    must be fixed BEFORE the scan, because the pure length pools are built in the same pass, so the
    choice has to be robust across every library this tool will ever see -- not tuned to cfRNA.

THE CEILING, AND WHY THIS IS ANSWERABLE EXACTLY
    Fragments arrive as a Poisson process. At one object the number of length-``w`` fragments in
    population ``p`` is Poisson with rate ``opp_p(w) * (rho_g f_g(w) + rho_r f_r(w))``, independently
    across ``w``. So the FULL LENGTH HISTOGRAM is the sufficient statistic and its Fisher information
    has a closed form -- no simulation, no approximation:

        I[a][b]  =  sum_{w,p}  opp_p(w) * f_a(w) * f_b(w) / (rho_g f_g(w) + rho_r f_r(w))

    Any stored summary is a projection of that histogram, so

        efficiency(set)  =  Var_full(phi_g) / Var_set(phi_g)   in [0, 1]

    is exactly "what fraction of the available information does this storage choice keep?".  That is
    the decision metric: if the shipped pair already keeps ~all of it, storing more is pointless; if
    it keeps a third, the rest is on the table.

WHY A THIRD MOMENT COULD HELP, PRECISELY
    Given the count, the fully efficient 1-D statistic is ``sum h*(w)`` with the score weight

        h*(w)  =  (f_g(w) - f_r(w)) / (phi f_g(w) + (1 - phi) f_r(w))

    which depends on the fitted length models and so is NOT available during the pass. But a stored
    triple ``(N, sum 1/L, sum L)`` spans ``{1, 1/w, w}``, and a linear combination of those CAN be
    formed AFTER the models are fit. So the triple's efficiency measures how well the optimal score
    is approximated in that fixed basis. That is what this script computes.

WHAT IS APPROXIMATED, AND HOW IT IS CHECKED
    ``I_full`` is exact (Poisson). A summary set's information uses the Gaussian/known-covariance
    approximation, which FLATTERS a heavy-tailed statistic like ``sum L`` -- so every headline cell
    is re-measured by Monte Carlo in ``--validate``, and the opportunity formulas themselves are
    verified by exact enumeration of integer start positions.

Usage::

    python scripts/design/observable_efficiency.py                 # the grid
    python scripts/design/observable_efficiency.py --validate      # + Monte Carlo check
    python scripts/design/observable_efficiency.py --shapes        # + lognormal / normal families
"""

from __future__ import annotations

import argparse
import itertools

import numpy as np

MAX_W = 3000
W = np.arange(MAX_W + 1, dtype=np.float64)

#: The owner's span: 50 bp is highly degraded material, 300 bp the long end. Applied to gDNA and RNA
#: INDEPENDENTLY -- there is no rule that RNA is longer, and assuming one is how a tool overfits to
#: cfRNA.
MEANS = (50.0, 75.0, 100.0, 150.0, 200.0, 250.0, 300.0)

#: Coefficient of variation: a narrow peak through to a very broad one. 1.0 is exponential-tailed.
CVS = (0.15, 0.35, 0.60, 1.00)

#: Node lengths spanning the real partition: 1 bp nodes exist, the median is 151 bp, the mean 2,970.
NODE_LENGTHS = (25, 50, 100, 151, 250, 400, 1000, 3000)


# ---------------------------------------------------------------------------
# fragment-length families
# ---------------------------------------------------------------------------


def fl_pmf(mean: float, cv: float, family: str = "gamma") -> np.ndarray:
    """A discrete fragment-length pmf on 1..MAX_W with this mean and coefficient of variation.

    Three families, because an earlier pass of this derivation reached the WRONG ranking from a
    single assumed shape: the gamma result reversed once the same moments were carried by a measured
    histogram. Shape is a variable here, not a background assumption.
    """
    from scipy.stats import gamma, lognorm, norm

    sd = mean * cv
    edges = np.arange(0, MAX_W + 2, dtype=np.float64)
    if family == "gamma":
        cdf = gamma.cdf(edges, a=(mean / sd) ** 2, scale=sd**2 / mean)
    elif family == "lognormal":
        sigma2 = np.log1p(cv**2)
        cdf = lognorm.cdf(edges, s=np.sqrt(sigma2), scale=mean / np.sqrt(1.0 + cv**2))
    elif family == "normal":
        cdf = norm.cdf(edges, loc=mean, scale=sd)
    else:
        raise ValueError(f"unknown family {family!r}")
    p = np.diff(cdf)
    p[0] = 0.0  # a zero-length fragment does not exist
    total = p.sum()
    if total <= 0:
        raise ValueError(f"{family}(mean={mean}, cv={cv}) put no mass on 1..{MAX_W}")
    return p / total


# ---------------------------------------------------------------------------
# the deposit rule's own geometry
# ---------------------------------------------------------------------------


def opportunity(population: str, node_len: float) -> np.ndarray:
    """Admissible start-position count per fragment length -- the deposit rule, per population."""
    if population == "contained":
        return np.maximum(node_len - W + 1.0, 0.0)
    if population == "spanning":
        return np.maximum(W - node_len - 1.0, 0.0)
    if population == "crossing":
        return np.maximum(W - 1.0, 0.0)
    raise ValueError(population)


#: Fixed, GLOBAL histogram bin edges -- chosen once and used for every object and every grid cell.
#: Tuning these per library (or per cell) would be cheating: the accumulator must pick its bins before
#: it has seen a fragment, exactly as it must pick a deposit weight.
def hist_edges(n_bins: int, scheme: str = "geometric", top: float = 1000.0) -> np.ndarray:
    """`n_bins` edges over plausible fragment lengths, with the last bin absorbing the overflow."""
    if scheme == "geometric":
        inner = np.geomspace(20.0, top, n_bins)
    elif scheme == "uniform":
        inner = np.linspace(20.0, top, n_bins)
    else:
        raise ValueError(scheme)
    return np.concatenate(([0.0], inner, [float(MAX_W + 1)]))[: n_bins + 1]


def weight(kind: str, population: str, node_len: float = 0.0) -> np.ndarray:
    """The per-fragment quantity summed into a channel.

    ``invL`` is the accumulator's own rule today, and its two cases are NOT one quantity: ``1/L`` at a
    node, ``1/(L-1)`` at a 0-bp line.

    ``invA`` is the RECIPROCAL OPPORTUNITY -- the derivation in
    ``scripts/design/node_density_derivation.py``. It is the unique weight for which E[sum h] is
    proportional to the start density with a length-distribution-independent constant, and the shipped
    edge rule is its ``A(w) = w - 1`` special case. At a node it is ``1/(l - L + 1)`` when contained and
    ``1/(L - l - 1)`` when spanning, neither of which has ever been tried.
    """
    if kind == "count":
        return np.ones_like(W)
    if kind == "sumL":
        return W.copy()
    if kind.startswith("hist"):
        # "hist<scheme><B>_<k>" -- the indicator of bin k of a B-bin fixed histogram. A histogram
        # channel stores an integer COUNT: no fixed point, no scale, no rounding, no overflow scheme.
        spec, k = kind[4:].split("_")
        scheme = "uniform" if spec[0] == "u" else "geometric"
        n_bins = int(spec[1:])
        edges = hist_edges(n_bins, scheme)
        lo, hi = edges[int(k)], edges[int(k) + 1] if int(k) + 1 < edges.size else float(MAX_W + 1)
        return ((W >= lo) & (W < hi)).astype(np.float64)
    out = np.zeros_like(W)
    if kind == "invA":
        A = opportunity(population, node_len)
        np.divide(1.0, A, out=out, where=A > 0)
        return out
    if kind == "invL":
        denom = (W - 1.0) if population == "crossing" else W
        np.divide(1.0, denom, out=out, where=denom > 0)
        return out
    raise ValueError(kind)


def verify_opportunity_counts() -> int:
    """EXACT brute force over integer start positions. Deterministic: no sampling, no tolerance.

    Checks the three opportunity formulas every divisor below is built from. Returns the mismatch
    count so a caller can fail on it.
    """
    bad = 0
    for node_len in (1, 2, 7, 25, 151, 400):
        for w in (1, 2, 3, 26, 150, 151, 152, 401, 900):
            starts = np.arange(-w - 2, node_len + w + 2)
            ends = starts + w
            checks = (
                ("contained", int(np.sum((starts >= 0) & (ends <= node_len))), max(node_len - w + 1, 0)),
                ("spanning", int(np.sum((starts < 0) & (ends > node_len))), max(w - node_len - 1, 0)),
                ("crossing", int(np.sum((starts < 0) & (ends > 0))), max(w - 1, 0)),
            )
            for name, got, want in checks:
                if got != want:
                    bad += 1
                    print(f"    MISMATCH {name} node={node_len} w={w}: enum {got} != formula {want}")
    return bad


# ---------------------------------------------------------------------------
# information
# ---------------------------------------------------------------------------

#: A "frame" is the set of populations that co-occur at ONE object, with its geometry.
FRAMES = {**{f"node {ell} bp": ("contained", "spanning", float(ell)) for ell in NODE_LENGTHS},
          "contiguous edge": ("crossing", None, 0.0)}


def _populations(frame):
    a, b, geom = FRAMES[frame]
    return [p for p in (a, b) if p is not None], geom


#: Candidate storage choices, as (population, channel) lists. `_expand` fills in the frame's
#: populations so one definition covers both node populations and the single edge population.
CANDIDATE_SETS = {
    "SHIPS: count + Sum1/L": [("*", "count"), ("*", "invL")],
    "count + Sum1/L + SumL": [("*", "count"), ("*", "invL"), ("*", "sumL")],
    "DERIVED: count + Sum1/A": [("*", "count"), ("*", "invA")],
    "count + Sum1/A + SumL": [("*", "count"), ("*", "invA"), ("*", "sumL")],
    "Sum1/A alone (both pops)": [("*", "invA")],
    "count + Sum1/A + Sum1/L": [("*", "count"), ("*", "invA"), ("*", "invL")],
    "ALL FOUR": [("*", "count"), ("*", "invA"), ("*", "invL"), ("*", "sumL")],
    # -- same bytes, spent on a length histogram instead of moments (uint32 bin vs uint64 moment) --
    "hist 3 bins  (= 12 B)": [("*", f"histg3_{k}") for k in range(3)],
    "hist 5 bins  (= 20 B)": [("*", f"histg5_{k}") for k in range(5)],
    "hist 7 bins  (= 28 B)": [("*", f"histg7_{k}") for k in range(7)],
    "hist 7 uniform (= 28 B)": [("*", f"histu7_{k}") for k in range(7)],
}


def _expand(channels, populations):
    out = []
    for pop, kind in channels:
        for p in (populations if pop == "*" else [pop]):
            if p in populations:
                out.append((p, kind))
    return out


def moments(pmf: np.ndarray, populations, geom: float) -> dict:
    """Per population, ``E[opp*h_j]`` and ``E[opp*h_j*h_k]`` for every channel pair.

    Precomputed once per (pmf, frame) so the grid sweep is matrix algebra rather than integration.
    """
    kinds = ("count", "invL", "invA", "sumL") + tuple(
        f"hist{p}{b}_{k}" for p in "gu" for b in (3, 5, 7) for k in range(b)
    )
    out = {}
    for pop in populations:
        opp = opportunity(pop, geom)
        h = {k: weight(k, pop, geom) for k in kinds}
        out[pop] = {
            "mean": {k: float(np.sum(pmf * opp * h[k])) for k in kinds},
            "cross": {
                (j, k): float(np.sum(pmf * opp * h[j] * h[k]))
                for j in kinds
                for k in kinds
            },
        }
    return out


def _var_phi(fisher: np.ndarray, rho_g: float, rho_r: float) -> float:
    """Delta-method Var(phi_g) from the 2x2 information about (rho_g, rho_r), phi = rho_g/(rho_g+rho_r)."""
    if not np.all(np.isfinite(fisher)) or np.linalg.cond(fisher) > 1e12:
        return float("inf")
    total = rho_g + rho_r
    grad = np.array([rho_r, -rho_g]) / total**2
    return float(grad @ np.linalg.inv(fisher) @ grad)


def var_full(pmf_g, pmf_r, populations, geom, rho_g, rho_r) -> float:
    """EXACT Var(phi_g) from the full length histogram -- the ceiling. Poisson, closed form."""
    fisher = np.zeros((2, 2))
    for pop in populations:
        opp = opportunity(pop, geom)
        lam = rho_g * pmf_g + rho_r * pmf_r
        live = (opp > 0) & (lam > 0)
        if not live.any():
            continue
        base = opp[live] / lam[live]
        fisher[0, 0] += float(np.sum(base * pmf_g[live] * pmf_g[live]))
        fisher[0, 1] += float(np.sum(base * pmf_g[live] * pmf_r[live]))
        fisher[1, 1] += float(np.sum(base * pmf_r[live] * pmf_r[live]))
    fisher[1, 0] = fisher[0, 1]
    return _var_phi(fisher, rho_g, rho_r)


def var_set(mom_g, mom_r, channels, rho_g, rho_r) -> float:
    """Gaussian/known-covariance Var(phi_g) from a stored summary set.

    ⚠ This APPROXIMATION flatters a heavy-tailed channel; ``--validate`` re-measures by Monte Carlo.
    Channels with no opportunity at this geometry are dropped -- an empty channel carries no
    information, which is not the same as the set being uninformative.

    ⚠ **The channels are STANDARDISED before the linear algebra**, and that is not cosmetic. The raw
    scales differ by ~1e10 (``E[opp*w^2]`` against ``E[opp/w^2]``), so a plain condition-number guard
    rejects perfectly informative sets purely for being badly scaled -- which showed up as ADDING a
    channel losing information, an impossibility that is how the bug was caught. Fisher information
    is invariant under an invertible linear change of the observables, so dividing each channel by
    its own sd changes no number that survives.
    """
    live = [
        (pop, kind)
        for pop, kind in channels
        if rho_g * mom_g[pop]["cross"][(kind, kind)] + rho_r * mom_r[pop]["cross"][(kind, kind)] > 0
    ]
    if len(live) < 2:
        return float("inf")
    M = np.array([[mom_g[p]["mean"][k], mom_r[p]["mean"][k]] for p, k in live])
    n = len(live)
    Sigma = np.zeros((n, n))
    for j, (pj, kj) in enumerate(live):
        for k, (pk, kk) in enumerate(live):
            if pj != pk:
                continue  # the populations are disjoint sets of fragments
            Sigma[j, k] = rho_g * mom_g[pj]["cross"][(kj, kk)] + rho_r * mom_r[pj]["cross"][(kj, kk)]

    scale = np.sqrt(np.diag(Sigma))
    if not np.all(np.isfinite(scale)) or np.any(scale <= 0):
        return float("inf")
    corr = Sigma / np.outer(scale, scale)  # a correlation matrix: diagonal 1, entries in [-1, 1]
    Mz = M / scale[:, None]
    # pinv, not solve: two channels can be exactly collinear (a degenerate but legitimate cell), and
    # the pseudo-inverse gives that set the information it genuinely has instead of rejecting it.
    fisher = Mz.T @ np.linalg.pinv(corr, rcond=1e-10) @ Mz
    return _var_phi(fisher, rho_g, rho_r)


# ---------------------------------------------------------------------------
# the sweep
# ---------------------------------------------------------------------------


def sweep(family: str, phi_g: float, rho_tot: float):
    """Efficiency of every candidate set over the whole (mean, cv) x (mean, cv) grid, per frame."""
    specs = list(itertools.product(MEANS, CVS))
    pmfs = {spec: fl_pmf(spec[0], spec[1], family) for spec in specs}
    rho_g, rho_r = phi_g * rho_tot, (1.0 - phi_g) * rho_tot

    results = {frame: {name: [] for name in CANDIDATE_SETS} for frame in FRAMES}
    for frame in FRAMES:
        populations, geom = _populations(frame)
        mom = {spec: moments(pmfs[spec], populations, geom) for spec in specs}
        for sg, sr in itertools.product(specs, specs):
            if sg == sr:
                continue  # identical pools carry no information about the split, by construction
            vf = var_full(pmfs[sg], pmfs[sr], populations, geom, rho_g, rho_r)
            if not np.isfinite(vf) or vf <= 0:
                continue
            for name, channels in CANDIDATE_SETS.items():
                chans = _expand(channels, populations)
                if not chans:
                    continue
                vs = var_set(mom[sg], mom[sr], chans, rho_g, rho_r)
                results[frame][name].append(vf / vs if np.isfinite(vs) and vs > 0 else 0.0)
    return results


def print_sweep(results, title: str) -> None:
    print(f"\n{'=' * 118}\n{title}\n{'=' * 118}")
    names = list(CANDIDATE_SETS)
    print(f"{'frame':>17} " + "".join(f"{n[:24]:>25}" for n in names))
    print(f"{'':>17} " + "".join(f"{'median   min   p05':>25}" for _ in names))
    for frame in FRAMES:
        row = f"{frame:>17} "
        for name in names:
            vals = np.asarray(results[frame][name], dtype=float)
            if vals.size == 0:
                row += f"{'--':>25}"
                continue
            row += f"{np.median(vals):>11.3f}{np.min(vals):>7.3f}{np.percentile(vals, 5):>7.3f}"
        print(row)
    print("\n  efficiency = Var(full length histogram) / Var(this stored set); 1.000 = loses nothing")
    print("  'min' and 'p05' are over the whole gDNA x RNA fragment-length grid -- the robustness")
    print("  numbers, since the deposit weight is fixed before any library is seen.")


# ---------------------------------------------------------------------------
# Monte Carlo validation
# ---------------------------------------------------------------------------


def validate(family: str, phi_g: float, rho_tot: float, n_trials: int, seed: int) -> None:
    """Re-measure a few cells by simulation, because the Gaussian approximation flatters `sum L`.

    Simulates the actual Poisson fragment process at one object, forms each stored set, and solves
    the 2x2 moment system per trial. Reports the realised sd against the predicted one.
    """
    print(f"\n{'=' * 118}\nMONTE CARLO CHECK -- does the Gaussian approximation hold? "
          f"({n_trials} trials/cell)\n{'=' * 118}")
    rng = np.random.default_rng(seed)
    rho_g, rho_r = phi_g * rho_tot, (1.0 - phi_g) * rho_tot
    cells = [
        ((100.0, 0.35), (200.0, 0.35), "node 400 bp"),
        ((200.0, 0.35), (100.0, 0.35), "node 400 bp"),   # the inverted case: gDNA longer
        ((100.0, 1.00), (200.0, 0.15), "node 400 bp"),   # broad vs narrow
        ((100.0, 0.35), (200.0, 0.35), "contiguous edge"),
        ((200.0, 1.00), (100.0, 0.15), "contiguous edge"),
    ]
    print(f"{'gDNA':>13}{'RNA':>13}{'frame':>17}{'set':>28}{'predicted sd':>15}{'realised sd':>14}"
          f"{'robust':>12}{'ratio':>8}")
    for sg, sr, frame in cells:
        populations, geom = _populations(frame)
        pmf_g, pmf_r = fl_pmf(*sg, family), fl_pmf(*sr, family)
        mom_g, mom_r = moments(pmf_g, populations, geom), moments(pmf_r, populations, geom)
        for name in ("SHIPS: count + Sum1/L", "DERIVED: count + Sum1/A", "count + Sum1/A + SumL"):
            chans = _expand(CANDIDATE_SETS[name], populations)
            pred = np.sqrt(var_set(mom_g, mom_r, chans, rho_g, rho_r))
            est = _mc_estimates(pmf_g, pmf_r, chans, populations, geom, rho_g, rho_r, n_trials, rng)
            if est.size <= 8:
                continue
            # sd AND a robust scale: they diverge exactly when a heavy tail is doing the damage,
            # which is the thing this check exists to detect.
            got = float(np.std(est))
            robust = float(np.subtract(*np.percentile(est, [75, 25]))) / 1.349
            ratio = got / pred if np.isfinite(pred) and pred > 0 else np.nan
            print(f"{sg[0]:>7.0f}/{sg[1]:<5.2f}{sr[0]:>7.0f}/{sr[1]:<5.2f}{frame:>17}"
                  f"{name[:26]:>28}{pred:>15.4f}{got:>14.4f}{robust:>12.4f}{ratio:>8.2f}")
    print("\n  ratio > 1 means the Gaussian prediction was OPTIMISTIC for that set -- i.e. the stored")
    print("  choice is WORSE in reality than the grid above claims. A realised sd far above the")
    print("  robust (IQR) scale means a heavy tail, which is the failure mode `Sum L` is suspected of.")


def _mc_estimates(pmf_g, pmf_r, channels, populations, geom, rho_g, rho_r, n_trials, rng):
    """Simulate the object `n_trials` times and solve for phi_g each time.

    ⚠ **No clipping and GLS, both deliberately.** A first version clipped phi_hat to [0,1] and solved
    by unweighted least squares; the clip truncates the spread and the wrong solver is less efficient
    than the Fisher calculation assumes, so the simulation came out 2x BETTER than the theory it
    exists to test -- flattering exactly the heavy-tailed channel under suspicion. The estimator here
    is the one the Fisher number describes: generalised least squares at the true covariance, and the
    raw unbounded estimate.
    """
    mom_g = moments(pmf_g, populations, geom)
    mom_r = moments(pmf_r, populations, geom)
    live = list(channels)
    M = np.array([[mom_g[p]["mean"][k], mom_r[p]["mean"][k]] for p, k in live])
    if np.linalg.matrix_rank(M) < 2:
        return np.empty(0)
    n = len(live)
    Sigma = np.zeros((n, n))
    for j, (pj, kj) in enumerate(live):
        for k, (pk, kk) in enumerate(live):
            if pj == pk:
                Sigma[j, k] = rho_g * mom_g[pj]["cross"][(kj, kk)] + rho_r * mom_r[pj]["cross"][(kj, kk)]
    scale = np.sqrt(np.diag(Sigma))
    Winv = np.linalg.pinv(Sigma / np.outer(scale, scale), rcond=1e-10)
    Mz = M / scale[:, None]
    solver = np.linalg.pinv(Mz.T @ Winv @ Mz) @ Mz.T @ Winv  # GLS: rho_hat = solver @ y_standardised

    lam = {p: (rho_g * pmf_g + rho_r * pmf_r) * opportunity(p, geom) for p in populations}
    out = []
    for _ in range(n_trials):
        draws = {p: rng.poisson(lam[p]) for p in populations}
        y = np.array([float(np.sum(draws[p] * weight(k, p))) for p, k in live]) / scale
        rho_hat = solver @ y
        tot = rho_hat.sum()
        if tot != 0.0:
            out.append(rho_hat[0] / tot)  # unbounded: truncating here is what hid the tail
    return np.asarray(out)


def check_monotonicity(family: str, phi_g: float, rho_tot: float) -> int:
    """A superset of channels can NEVER carry less information than a subset. Assert it.

    This is the harness's own falsification test, and it earned its place: the first run reported
    ``count + Sum1/L + SumL`` at 0.000 where ``count + Sum1/L`` scored 0.832, which is impossible and
    was a covariance-conditioning failure rather than a result. Without this check the write-up would
    have argued from it. Returns the number of violations.
    """
    rho_g, rho_r = phi_g * rho_tot, (1.0 - phi_g) * rho_tot
    nested = [
        ([("*", "count"), ("*", "invL")], [("*", "count"), ("*", "invL"), ("*", "sumL")]),
        ([("*", "count"), ("*", "invA")], [("*", "count"), ("*", "invA"), ("*", "sumL")]),
        ([("*", "count"), ("*", "invA")], [("*", "count"), ("*", "invA"), ("*", "invL")]),
        ([("*", "count")], [("*", "count"), ("*", "invA")]),
        ([("contained", "count"), ("contained", "invA")], [("*", "count"), ("*", "invA")]),
    ]
    specs = [(m, c) for m in (50.0, 150.0, 300.0) for c in (0.15, 0.6, 1.0)]
    pmfs = {s: fl_pmf(s[0], s[1], family) for s in specs}
    bad = 0
    for frame in FRAMES:
        populations, geom = _populations(frame)
        mom = {s: moments(pmfs[s], populations, geom) for s in specs}
        for sg, sr in itertools.product(specs, specs):
            if sg == sr:
                continue
            for small, big in nested:
                vs = var_set(mom[sg], mom[sr], _expand(small, populations), rho_g, rho_r)
                vb = var_set(mom[sg], mom[sr], _expand(big, populations), rho_g, rho_r)
                if np.isfinite(vs) and vb > vs * 1.001:
                    bad += 1
                    if bad <= 5:
                        print(f"    VIOLATION {frame} gDNA{sg} RNA{sr}: subset var {vs:.4g} < "
                              f"superset var {vb:.4g}")
    return bad


def hardest_cells(family: str, phi_g: float, rho_tot: float, frame: str, name: str, n: int = 6):
    """The grid cells where one storage choice loses the most -- what a robustness number is made of."""
    populations, geom = _populations(frame)
    specs = list(itertools.product(MEANS, CVS))
    pmfs = {s: fl_pmf(s[0], s[1], family) for s in specs}
    mom = {s: moments(pmfs[s], populations, geom) for s in specs}
    rho_g, rho_r = phi_g * rho_tot, (1.0 - phi_g) * rho_tot
    chans = _expand(CANDIDATE_SETS[name], populations)
    rows = []
    for sg, sr in itertools.product(specs, specs):
        if sg == sr:
            continue
        vf = var_full(pmfs[sg], pmfs[sr], populations, geom, rho_g, rho_r)
        if not np.isfinite(vf) or vf <= 0:
            continue
        vs = var_set(mom[sg], mom[sr], chans, rho_g, rho_r)
        rows.append((vf / vs if np.isfinite(vs) and vs > 0 else 0.0, sg, sr))
    rows.sort(key=lambda r: r[0])
    print(f"\n  hardest cells for '{name}' at {frame}:")
    for eff, sg, sr in rows[:n]:
        print(f"    efficiency {eff:6.3f}   gDNA mean {sg[0]:5.0f} cv {sg[1]:4.2f}   "
              f"RNA mean {sr[0]:5.0f} cv {sr[1]:4.2f}"
              f"{'    <- SAME MEAN' if sg[0] == sr[0] else ''}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--phi", type=float, default=0.5, help="true gDNA share of the start density")
    ap.add_argument("--rho", type=float, default=0.05, help="total start density, fragments per bp")
    ap.add_argument("--validate", action="store_true", help="Monte Carlo check of the approximation")
    ap.add_argument("--shapes", action="store_true", help="repeat over lognormal and normal families")
    ap.add_argument("--trials", type=int, default=4000)
    ap.add_argument("--seed", type=int, default=20260730)
    args = ap.parse_args()

    print("exact enumeration of the opportunity counts (deterministic, no sampling):")
    bad = verify_opportunity_counts()
    print(f"  {'ALL MATCH' if bad == 0 else f'{bad} MISMATCHES'}")
    if bad:
        raise SystemExit("opportunity formulas are wrong; every number below would be meaningless")

    print("\nharness self-check -- a superset of channels must never carry less information:")
    violations = check_monotonicity("gamma", args.phi, args.rho)
    print(f"  {'NONE' if violations == 0 else f'{violations} VIOLATIONS'}")
    if violations:
        raise SystemExit("the information calculation is inconsistent; the grid below is meaningless")

    print(f"\ngrid: gDNA and RNA fragment lengths INDEPENDENTLY over means {MEANS} bp")
    print(f"      x cv {CVS}, both directions (no assumption that RNA is longer)")
    print(f"      {len(MEANS) * len(CVS)} pools -> {len(MEANS) * len(CVS) * (len(MEANS) * len(CVS) - 1)} ordered pairs per frame")

    families = ("gamma", "lognormal", "normal") if args.shapes else ("gamma",)
    for family in families:
        print_sweep(
            sweep(family, args.phi, args.rho),
            f"{family.upper()} fragment lengths   (true phi_g={args.phi}, rho_tot={args.rho}/bp)",
        )

    for frame in ("node 151 bp", "contiguous edge"):
        hardest_cells("gamma", args.phi, args.rho, frame, "SHIPS: count + Sum1/L")
        hardest_cells("gamma", args.phi, args.rho, frame, "DERIVED: count + Sum1/A")
        hardest_cells("gamma", args.phi, args.rho, frame, "count + Sum1/A + SumL")

    if args.validate:
        validate("gamma", args.phi, args.rho, args.trials, args.seed)


if __name__ == "__main__":
    main()
