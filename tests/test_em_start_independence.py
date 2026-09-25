"""The locus EM's answer does not depend on where it starts.

SQUAREM jumps ahead along the EM's own path. Where a jump carried a shrinking component below the floor, the
solver once CLAMPED it there, and a component at the floor takes no responsibility in the E-step and no share of
the evidence-proportional prior, so it never came back: which of the components sharing fragments survived was
decided by the warm start. The solver now shortens the jump toward the plain double step instead
(``em_solver.cpp``'s ``backtracked_squarem_step``), so only the EM's own step can take a component to the floor.

The fixtures are three real loci from the ladder's g00 ss.99 capture-ON condition — every array the locus EM
reads, the locus's transcripts renumbered 0..n-1, and the two pseudocounts ``pipeline.em_pseudocounts`` gave
them — each of which the clamp forked: the coverage and the uniform start ended 7 to 72 fragments apart, in the
gDNA component or a transcript. Solved from both starts at 10,000 iterations and δ 1e-9 they now agree to within
1.3e-5 fragments; the bound below sits between the two by orders of magnitude on each side. Only the (locus,
mode) pairs the clamp forked are listed, so every case can fire.

The ``prior`` start is left out on purpose: with no per-component weight it starts every RNA component at exactly
zero, the absorbing state by design (the warm start's own comment in ``em_solver.cpp``).
"""

from pathlib import Path

import numpy as np
import pytest

from rigel.config import EMConfig, TranscriptGeometry
from rigel.estimator import AbundanceEstimator

FIXTURES = Path(__file__).parent / "data" / "em_start"
PARTITION_KEYS = (
    "offsets",
    "t_indices",
    "log_liks",
    "coverage_weights",
    "count_cols",
    "is_spliced",
    "gdna_log_liks",
    "locus_t_indices",
    "locus_count_cols",
)
# The clamp forked these pairs between the coverage and the uniform start (7 to 72 fragments apart).
FORKED = [
    ("locus303", "map"),
    ("locus751", "map"),
    ("locus751", "vbem"),
    ("locus803", "map"),
    ("locus803", "vbem"),
]
# Fragments. The repaired solver's two answers agree to 1.3e-5 at most; the clamp's were 7 or more apart.
AGREE = 1e-3


def _solve(locus: str, mode: str, warm_start: str):
    """One captured locus through the shipped batch EM: each component's fractional count (transcripts, then the
    gDNA component) and the solver's per-locus profile."""
    fx = dict(np.load(FIXTURES / f"g00_ss099_on_{locus}.npz"))
    n_t = int(fx["t_eff_len_em"].size)
    geometry = TranscriptGeometry(
        effective_lengths=np.maximum(fx["t_eff_len_em"], 1.0),
        effective_lengths_em=fx["t_eff_len_em"],
    )
    config = EMConfig(
        mode=mode,
        warm_start=warm_start,
        assignment_mode="fractional",
        iterations=10000,
        convergence_delta=1e-9,
        n_threads=1,
    )
    rc = AbundanceEstimator(n_t, geometry=geometry, em_config=config)
    rc.unambig_counts[:] = fx["unambig"]
    gdna, _, _ = rc.run_batch_locus_em_partitioned(
        [tuple(fx[k] for k in PARTITION_KEYS)],
        [np.arange(n_t, dtype=np.int32)],
        np.array([float(fx["gdna_prior"])]),
        rna_prior_count=np.array([float(fx["rna_prior"])]),
        gdna_eff_len=np.array([float(fx["gdna_eff_len"])]),
        em_iterations=config.iterations,
        em_convergence_delta=config.convergence_delta,
        emit_locus_stats=True,
    )
    return np.append(rc.em_counts.sum(axis=1), gdna), rc.locus_stats[0]


@pytest.mark.parametrize("locus,mode", FORKED)
def test_the_coverage_and_the_uniform_start_reach_one_answer(locus, mode):
    coverage, _ = _solve(locus, mode, "coverage")
    uniform, _ = _solve(locus, mode, "uniform")
    moved = 0.5 * np.abs(coverage - uniform).sum()
    assert moved < AGREE, f"{locus} {mode}: the two starts end {moved:.4g} fragments apart"


@pytest.mark.parametrize("locus", sorted({locus for locus, _ in FORKED}))
def test_each_captured_locus_exercises_the_backtrack(locus):
    """The repair's own path runs on every fixture: some extrapolation step was shortened rather than taken."""
    halvings = sum(
        _solve(locus, mode, warm)[1]["squarem_backtrack_count"]
        for mode in ("map", "vbem")
        for warm in ("coverage", "uniform")
    )
    assert halvings > 0
