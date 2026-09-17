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


def _read_out(SL, psi, lam) -> float:
    from scipy.special import logsumexp

    post = np.exp(psi - logsumexp(psi, axis=1, keepdims=True))
    return float(np.clip(SL.posterior_median_fg(post, lam), 0.0, 1.0)[0])


def test_the_budget_covers_term_and_intermediate_rounding_and_its_constants_are_load_bearing():
    """Brute force against the derivation, at the two places an implementation rounds. A slot's ψ is
    its strand term plus the arms plus a λ-factor row; here the row is the strand profile of a slot of
    the OPPOSITE composition, so the two large terms disagree and cancel at the mode — the case the
    bound prices. (i) Every term rounded to float32 before the sum, one rounding unit of each term's
    magnitude; (ii) the strand term evaluated in float32 arithmetic, so the mean ``n·p`` is rounded
    before ``u − n·p`` cancels. Both must stay inside the float32 budget at every strength; the
    achieved ratio is recorded as orders inside it, a bound being a bound."""
    import sys
    from pathlib import Path

    from rigel.calibration import simplex_logodds as SL

    sys.path.insert(0, str(Path(__file__).resolve().parent / "calibration"))
    from _psi_reference import jeffreys_arms, strand_loglik_mixture

    sr = _replay()
    kappa, window, K, od = 0.99, 10.0, 138, 0.0  # no overdispersion: the count scaling is real
    lam, fg = SL._logodds_grid(K, window)
    F = np.float32

    def strand(u_pos, n, share, dtype):
        f_act = (1.0 - fg)[None, :]
        return strand_loglik_mixture(
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
            arms = jeffreys_arms(lam)[None, :]
            terms = [strand(u_a, n, share, np.float64), arms, strand(u_b, n, other, np.float64)]
            exact = _read_out(SL, sum(terms), lam)
            rounded_terms = sum(np.asarray(np.asarray(t, F), np.float64) for t in terms)
            f32_terms = _read_out(SL, rounded_terms, lam)
            f32_inner = _read_out(SL, strand(u_a, n, share, F) + arms + terms[2], lam)
            b_frac, _b_var = sr.budget(n, kappa, window, sr.EPS32, K)
            for got in (f32_terms, f32_inner):
                move = abs(got - exact)
                worst = max(worst, move / b_frac)
                assert move <= b_frac, (n, share, move, b_frac)
    assert worst > 0.0, "nothing moved under float32 rounding: the brute force tested nothing"
    # the achieved worst ratio sits orders inside the bound (rounding errors on neighbouring cells are
    # incoherent and the read-out averages them); a bound is what it is — this records the scale
    assert worst < 1e-2, f"the brute force came within 1 % of the bound: re-derive it ({worst:.3e})"


def test_an_ambig_slot_is_held_to_the_float64_budget():
    """ψ is float64 at every slot, AMBIG ones included, so the report must hold a move at an AMBIG slot to
    the float64 budget: a 1e-9 move there is orders past it and must read BEYOND, where a float32 unit
    would call it within."""
    import dataclasses
    from types import SimpleNamespace

    sr = _replay()

    @dataclasses.dataclass
    class Belief:
        f_g: np.ndarray

    n = 4
    # every slot AMBIG: both strands admissible
    statics = SimpleNamespace(free_pos=np.ones(n, bool), free_neg=np.ones(n, bool))
    geometry = SimpleNamespace(
        unspliced_count=np.full((n, 2), 50.0), spliced_count=np.zeros((n, 2))
    )
    args = (None, statics, geometry, None, None)
    kwargs = {"rna_sense_frac": 0.99, "logodds_window": 10.0, "n_grid": 138}
    expected = Belief(f_g=np.full(n, 0.25))
    moved = expected.f_g.copy()
    moved[2] += 1e-9
    b_frac, _ = sr.budget(100.0, 0.99, 10.0, sr.EPS64, 138)
    assert 1e-9 > b_frac, "the move must exceed the float64 budget for this gate to mean anything"
    report = sr.tolerance_report(expected, Belief(f_g=moved), args, kwargs)
    line = next(row for row in report if "f_g" in row)
    assert "BEYOND the budget" in line, line
