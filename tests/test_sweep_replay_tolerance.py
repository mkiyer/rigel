"""The replay's tolerance report (`scripts/profiling/sweep_replay.py --tolerance`): the comparator's own
falsification runs here as a suite gate, and the derived budget is checked against brute force — ψ read
out from a log-density rounded to float32, a perturbation of every term at once bounded by ``ε₃₂·A``,
must move no slot past the float32 budget. Run without strand overdispersion, where the strand term's
magnitude scales with the count as the bound says (overdispersion caps it at ``c_κ/od``, which only
loosens the bound). Gate for the derivation in `sweep_replay`'s docstring.
"""

from __future__ import annotations

import importlib.util
import pathlib

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]


def _replay():
    spec = importlib.util.spec_from_file_location(
        "sweep_replay", ROOT / "scripts" / "profiling" / "sweep_replay.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_the_comparator_fires_its_own_gates():
    assert _replay().self_test() == 0


def _read_out(SL, psi, lam, fg) -> float:
    post = np.exp(psi - SL._lse(psi, axis=1, keepdims=True))
    return float(np.clip(SL._posterior_median_fg(post, lam, fg), 0.0, 1.0)[0])


def test_the_budget_covers_term_and_intermediate_rounding_and_its_constants_are_load_bearing():
    """Brute force against the derivation, at the two places an implementation rounds. A slot's ψ is
    its strand term plus the arms plus a λ-factor row; here the row is the strand profile of a slot of
    the OPPOSITE composition, so the two large terms disagree and cancel at the mode — the case the
    bound prices. (i) Every term rounded to float32 before the sum, one rounding unit of each term's
    magnitude; (ii) the strand term evaluated in float32 arithmetic, so the mean ``n·p`` is rounded
    before ``u − n·p`` cancels. Both must stay inside the float32 budget at every strength; the
    achieved ratio is recorded as orders inside it, a bound being a bound."""
    from rigel.calibration import simplex_logodds as SL

    sr = _replay()
    kappa, window, K, od = 0.99, 10.0, 138, 0.0  # no overdispersion: the count scaling is real
    lam, fg = SL._logodds_grid(K, window)
    F = np.float32

    def strand(u_pos, n, share, dtype):
        f_act = (1.0 - fg)[None, :]
        return SL._mixture_strand_loglik(
            np.asarray([u_pos], dtype)[:, None],
            np.asarray([n], dtype)[:, None],
            fg.astype(dtype)[None, :],
            f_act.astype(dtype),
            np.zeros((1, K), dtype),
            kappa,
            od,
            od,
            np.asarray([share], dtype)[:, None],
            np.asarray([1.0 - share], dtype)[:, None],
            np.asarray([0.0], dtype)[:, None],
        ).astype(np.float64)

    worst = 0.0
    for n in (10.0, 1.0e3, 1.0e5, 1.0e6):
        for share, other in ((0.001, 0.999), (0.999, 0.001), (0.02, 0.98), (0.3, 0.7)):
            u_a = n * (0.5 * share + kappa * (1.0 - share))
            u_b = n * (0.5 * other + kappa * (1.0 - other))
            arms = SL._gdna_arm(lam, None) + SL._rna_arm(lam)
            terms = [strand(u_a, n, share, np.float64), arms, strand(u_b, n, other, np.float64)]
            exact = _read_out(SL, sum(terms), lam, fg)
            rounded_terms = sum(np.asarray(np.asarray(t, F), np.float64) for t in terms)
            f32_terms = _read_out(SL, rounded_terms, lam, fg)
            f32_inner = _read_out(SL, strand(u_a, n, share, F) + arms + terms[2], lam, fg)
            b_frac, _b_var = sr.budget(n, kappa, window, sr.EPS32, K)
            for got in (f32_terms, f32_inner):
                move = abs(got - exact)
                worst = max(worst, move / b_frac)
                assert move <= b_frac, (n, share, move, b_frac)
    assert worst > 0.0, "nothing moved under float32 rounding: the brute force tested nothing"
    # the achieved worst ratio sits orders inside the bound (rounding errors on neighbouring cells are
    # incoherent and the read-out averages them); a bound is what it is — this records the scale
    assert worst < 1e-2, f"the brute force came within 1 % of the bound: re-derive it ({worst:.3e})"
