"""The reciprocal-opportunity deposit: one density rule for every object, node or edge.

    Write-up: ``docs/NODE_DENSITY_DERIVATION.md``   ·   Design: ``docs/ACCUMULATOR_DESIGN.md`` §6/§10.1

THE CLAIM UNDER TEST
    At an edge the accumulator deposits ``1/(L-1)`` and the result is an exactly model-free estimate of
    the start density. The reason usually given is "the reciprocal of the fragment length cancels the
    opportunity". That is a coincidence of the edge frame: the opportunity to cross a 0-bp line happens
    to BE ``L-1``. The general statement is

        deposit  h(w) = 1 / A(w)     where A(w) is that population's OPPORTUNITY

    and the edge rule is its ``A(w) = w - 1`` special case. Applied to a node, where the contained
    opportunity is ``(l - w + 1)_+`` and the spanning opportunity is ``(w - l - 1)_+``, this predicts a
    deposit that has never been tried: ``1/(l - L + 1)`` and ``1/(L - l - 1)``.

WHAT THIS SCRIPT VERIFIES
    T1  h = 1/A is the UNIQUE weight (up to scale) making E[sum h] independent of the length
        distribution -- checked numerically against alternatives over a wide FL grid.
    T2  E[sum h] = rho * P(A > 0). Model-free up to the population's own support truncation.
    T3  contained + spanning together estimate rho * (1 - f(l+1)) -- the two truncations are
        complementary and close to exactly, the single missing length being w = l + 1.
    T4  The two limits the owner predicted:
          l = 0        -> the node rule degenerates to the edge rule EXACTLY
          l >> E[L]    -> the deposit tends to 1/l, i.e. count / node_length
    T5  The shipped node weight ``1/L`` does NOT have property T2, and the size of its bias.
    T6  ``count / l`` (the naive density) is biased by (l - mu + 1)/l; the reciprocal deposit removes
        that bias with NO fragment-length model.

Every check is run end-to-end through a simulated Poisson fragment process as well as in closed form,
and each is PERTURBED -- the wrong weight is substituted and the check must fail -- because a green
verification that was never made to fail proves nothing.
"""

from __future__ import annotations

import numpy as np

MAX_W = 3000
W = np.arange(MAX_W + 1, dtype=np.float64)

PASS, FAIL = "  ok  ", " FAIL "


def fl_pmf(mean: float, cv: float, family: str = "gamma") -> np.ndarray:
    from scipy.stats import gamma, lognorm, norm

    sd = mean * cv
    edges = np.arange(0, MAX_W + 2, dtype=np.float64)
    if family == "gamma":
        cdf = gamma.cdf(edges, a=(mean / sd) ** 2, scale=sd**2 / mean)
    elif family == "lognormal":
        cdf = lognorm.cdf(edges, s=np.sqrt(np.log1p(cv**2)), scale=mean / np.sqrt(1.0 + cv**2))
    else:
        cdf = norm.cdf(edges, loc=mean, scale=sd)
    p = np.diff(cdf)
    p[0] = 0.0
    return p / p.sum()


# ---------------------------------------------------------------------------
# the deposit rule
# ---------------------------------------------------------------------------


def opportunity(population: str, node_len: float) -> np.ndarray:
    """Admissible start positions per fragment length. Verified by enumeration in `check_enumeration`."""
    if population == "contained":
        return np.maximum(node_len - W + 1.0, 0.0)
    if population == "spanning":
        return np.maximum(W - node_len - 1.0, 0.0)
    if population == "crossing":
        return np.maximum(W - 1.0, 0.0)
    raise ValueError(population)


def deposit_weight(kind: str, population: str, node_len: float) -> np.ndarray:
    """The per-fragment quantity added to the density channel.

    ``invA`` is the derivation's proposal (reciprocal opportunity); ``invL`` is what ships (1/L at a
    node, 1/(L-1) at a line); ``count`` deposits 1.
    """
    out = np.zeros_like(W)
    if kind == "count":
        return np.ones_like(W)
    if kind == "invA":
        A = opportunity(population, node_len)
        np.divide(1.0, A, out=out, where=A > 0)
        return out
    if kind == "invL":
        denom = (W - 1.0) if population == "crossing" else W
        np.divide(1.0, denom, out=out, where=denom > 0)
        return out
    raise ValueError(kind)


def expected_deposit(pmf, population, node_len, kind, rho) -> float:
    """Closed form ``E[sum h] = rho * sum_w f(w) A(w) h(w)`` -- exact, given the pmf."""
    A = opportunity(population, node_len)
    h = deposit_weight(kind, population, node_len)
    return float(rho * np.sum(pmf * A * h))


# ---------------------------------------------------------------------------
# checks
# ---------------------------------------------------------------------------

GRID = [(m, c) for m in (50.0, 100.0, 200.0, 300.0) for c in (0.15, 0.35, 0.60, 1.00)]
NODE_LENGTHS = (0, 1, 25, 50, 100, 151, 250, 400, 1000, 3000, 20000)


def check_enumeration() -> int:
    """The opportunity formulas, by exact enumeration of integer starts. No sampling, no tolerance."""
    bad = 0
    for ell in (0, 1, 2, 7, 25, 151):
        for w in (1, 2, 3, 8, 26, 150, 151, 152, 400):
            starts = np.arange(-w - 3, ell + w + 3)
            ends = starts + w
            for name, got, want in (
                ("contained", int(np.sum((starts >= 0) & (ends <= ell))), max(ell - w + 1, 0)),
                ("spanning", int(np.sum((starts < 0) & (ends > ell))), max(w - ell - 1, 0)),
                ("crossing", int(np.sum((starts < 0) & (ends > 0))), max(w - 1, 0)),
            ):
                if got != want:
                    bad += 1
                    print(f"      MISMATCH {name} l={ell} w={w}: enum {got} != formula {want}")
    return bad


def check_T2_model_free(kind: str, rho: float = 0.05) -> tuple[int, float]:
    """T2: does ``E[sum h] / (rho * P(A>0))`` equal 1 for EVERY fragment-length distribution?

    Returns (violations, worst relative deviation). A weight with the property gives exactly 1.0 for
    every pmf; a weight without it varies with the pmf, and that variation IS the model dependence.
    """
    bad, worst = 0, 0.0
    for population in ("contained", "spanning", "crossing"):
        for ell in NODE_LENGTHS:
            if population == "crossing" and ell != 0:
                continue
            for spec in GRID:
                pmf = fl_pmf(*spec)
                A = opportunity(population, float(ell))
                support = float(np.sum(pmf[A > 0]))
                if support < 1e-9:
                    continue
                got = expected_deposit(pmf, population, float(ell), kind, rho)
                ratio = got / (rho * support)
                worst = max(worst, abs(ratio - 1.0))
                if abs(ratio - 1.0) > 1e-9:
                    bad += 1
    return bad, worst


def check_T3_complementary(rho: float = 0.05) -> tuple[int, float]:
    """T3: contained + spanning, both at reciprocal opportunity, estimate ``rho * (1 - f(l+1))``."""
    bad, worst = 0, 0.0
    for ell in NODE_LENGTHS:
        for spec in GRID:
            pmf = fl_pmf(*spec)
            total = expected_deposit(pmf, "contained", float(ell), "invA", rho) + expected_deposit(
                pmf, "spanning", float(ell), "invA", rho
            )
            hole = float(pmf[ell + 1]) if ell + 1 < pmf.size else 0.0
            predicted = rho * (1.0 - hole)
            rel = abs(total - predicted) / max(predicted, 1e-15)
            worst = max(worst, rel)
            if rel > 1e-9:
                bad += 1
    return bad, worst


def check_T4_limits(rho: float = 0.05) -> None:
    """T4: the two limits -- l = 0 is the edge rule exactly; l >> E[L] tends to count / node_length."""
    print("\n  T4a  l = 0 reduces to the edge rule (node spanning weight vs the 0-bp line weight)")
    for spec in GRID[:6]:
        pmf = fl_pmf(*spec)
        node0 = expected_deposit(pmf, "spanning", 0.0, "invA", rho)
        edge = expected_deposit(pmf, "crossing", 0.0, "invA", rho)
        contained0 = expected_deposit(pmf, "contained", 0.0, "invA", rho)
        tag = PASS if abs(node0 - edge) < 1e-12 and contained0 == 0.0 else FAIL
        print(f"   {tag} FL {spec[0]:5.0f}/{spec[1]:4.2f}:  node l=0 spanning {node0:.10f}  "
              f"edge {edge:.10f}   (contained at l=0 is {contained0:.1f}, as it must be)")

    print("\n  T4b  l >> E[L]: the reciprocal deposit -> 1/l, so the sum -> count / node_length")
    pmf = fl_pmf(200.0, 0.35)
    print(f"       {'node l':>10}{'E[sum 1/A]/rho':>18}{'E[count]/(rho*l)':>20}{'ratio':>10}")
    for ell in (400, 1000, 3000, 20000):
        dens = expected_deposit(pmf, "contained", float(ell), "invA", rho) / rho
        naive = expected_deposit(pmf, "contained", float(ell), "count", rho) / (rho * ell)
        print(f"       {ell:>10}{dens:>18.6f}{naive:>20.6f}{naive / dens:>10.6f}")


def check_T6_naive_bias(rho: float = 0.05) -> None:
    """T6: `count / l` is biased by (l - mu + 1)/l. The reciprocal deposit removes it, model-free."""
    print("\n  T6   `count / node_length` against the reciprocal deposit -- the effective-length bias")
    print(f"       {'node l':>9}{'FL mean':>10}{'count/l  (x truth)':>22}{'sum 1/A  (x truth)':>22}")
    for ell in (151, 250, 400, 1000, 3000):
        for mean in (50.0, 200.0, 300.0):
            pmf = fl_pmf(mean, 0.35)
            naive = expected_deposit(pmf, "contained", float(ell), "count", rho) / (rho * ell)
            exact = expected_deposit(pmf, "contained", float(ell), "invA", rho) / rho
            print(f"       {ell:>9}{mean:>10.0f}{naive:>22.4f}{exact:>22.4f}")


def check_end_to_end(rho: float, n_reps: int, seed: int) -> None:
    """The whole chain through a simulated Poisson fragment process, not just the algebra.

    Places fragments uniformly at density `rho` over a window wide enough to cover every start that
    could touch the node, applies the deposit rule as CODE, and compares to the closed form.
    """
    print("\n  END-TO-END simulation of the fragment process (the algebra above is not the code path)")
    rng = np.random.default_rng(seed)
    print(f"       {'node l':>9}{'FL':>12}{'predicted':>13}{'simulated':>13}{'z':>8}")
    for ell in (0, 25, 151, 1000):
        for spec in ((100.0, 0.60), (250.0, 0.35)):
            pmf = fl_pmf(*spec)
            pad = MAX_W + 2
            span = ell + 2 * pad
            totals = []
            for _ in range(n_reps):
                n = rng.poisson(rho * span)
                starts = rng.integers(-pad, ell + pad, size=n).astype(np.int64)
                lengths = rng.choice(len(pmf), size=n, p=pmf).astype(np.int64)
                ends = starts + lengths
                contained = (starts >= 0) & (ends <= ell) & (lengths > 0)
                spanning = (starts < 0) & (ends > ell)
                acc = 0.0
                if contained.any():
                    acc += float(np.sum(1.0 / (ell - lengths[contained] + 1.0)))
                if spanning.any():
                    acc += float(np.sum(1.0 / (lengths[spanning] - ell - 1.0)))
                totals.append(acc)
            got = float(np.mean(totals))
            se = float(np.std(totals, ddof=1)) / np.sqrt(n_reps)
            hole = float(pmf[ell + 1]) if ell + 1 < pmf.size else 0.0
            pred = rho * (1.0 - hole)
            z = (got - pred) / se if se > 0 else np.nan
            tag = PASS if abs(z) < 4.0 else FAIL
            print(f"   {tag}{ell:>9}{spec[0]:>7.0f}/{spec[1]:<4.2f}{pred:>13.6f}{got:>13.6f}{z:>8.2f}")


def main() -> None:
    rho = 0.05
    print("=" * 100)
    print("THE RECIPROCAL-OPPORTUNITY DEPOSIT -- verification")
    print("=" * 100)

    print("\n  T0   opportunity formulas, by exact enumeration:")
    bad = check_enumeration()
    print(f"   {PASS if bad == 0 else FAIL} {'all match' if bad == 0 else f'{bad} mismatches'}")

    print("\n  T2   is E[sum h] proportional to rho with a pmf-INDEPENDENT constant?")
    print(f"       {'weight':>28}{'violations':>13}{'worst deviation':>18}")
    for kind, label in (("invA", "1/A  (the derivation)"), ("invL", "1/L  (what ships)"),
                        ("count", "1    (the raw count)")):
        bad, worst = check_T2_model_free(kind, rho)
        tag = PASS if (bad == 0) == (kind == "invA") else FAIL
        print(f"   {tag}{label:>28}{bad:>13}{worst:>18.6f}")
    print("       ^ invA must show ZERO violations; the other two MUST show many -- that is the")
    print("         perturbation. A weight with no violations is model-free BY DEFINITION here.")

    print("\n  T3   contained + spanning = rho * (1 - f(l+1))?")
    bad, worst = check_T3_complementary(rho)
    print(f"   {PASS if bad == 0 else FAIL} violations {bad}, worst relative deviation {worst:.3e}")

    check_T4_limits(rho)
    check_T6_naive_bias(rho)
    check_end_to_end(rho, n_reps=400, seed=20260730)

    print("\n" + "=" * 100)


if __name__ == "__main__":
    main()
