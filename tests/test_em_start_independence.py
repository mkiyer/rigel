"""SQUAREM does not let the start choose the locus EM's answer.

SQUAREM jumps ahead along the EM's own path. A jump that clamped a shrinking component at the floor would let the
warm start decide which of the components sharing fragments survives, since a component at the floor takes no
responsibility in the E-step and no share of the evidence-proportional prior. The solver shortens the jump toward
the plain double step instead (``em_squarem.h``'s ``backtracked_squarem_step``), so only the EM's own step can take
a component to the floor.

The EM has one warm start, the coverage seed, so the other starts are made through it. The seed hands each ambiguous
fragment to its candidates in proportion to their coverage weights, the gDNA candidate's being 1, and the weights
enter nothing else. So a start that raises one transcript's weights a hundredfold leans the fragments it shares toward
that transcript, and one that lowers every transcript's leans every unspliced fragment toward gDNA.

The fixtures are three captured loci (every array the locus EM reads, the transcripts renumbered 0..n-1, and the two
pseudocounts ``pipeline.em_pseudocounts`` gave them). Listed are the (locus, mode) pairs whose plain, unaccelerated EM
reaches one answer from every leaning start, and on each a clamping jump forks the coverage seed from some leaning
start by whole fragments, so every case can fire. VBEM at locus751 is not listed: two of its transcripts compete for
the same fragments and whichever starts ahead keeps them, under the plain EM as well, so VBEM itself has two answers
there and the start rightly picks one.
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
)
LOCI = ("locus303", "locus751", "locus803")
MODES = ("map", "vbem")
PAIRS = [
    ("locus303", "map"),
    ("locus303", "vbem"),
    ("locus751", "map"),
    ("locus803", "map"),
    ("locus803", "vbem"),
]
# A hundredfold: under a clamping jump every listed pair forks at a lean of 30 (locus303 does not at 10), and VBEM's
# own answer at locus303 and locus803 moves, under the plain EM too, only between a lean of 1e4 and 1e5.
LEAN = 100.0
# Fragments: orders of magnitude below the whole fragments a clamping jump forks the starts by, and above the
# solver's own convergence residue at 10,000 iterations and δ 1e-9.
AGREE = 1e-3


def _fixture(locus: str) -> dict:
    return dict(np.load(FIXTURES / f"g00_ss099_on_{locus}.npz"))


def _leans(locus: str) -> list:
    """What a start can lean on: gDNA, then each transcript."""
    return ["gdna", *range(int(_fixture(locus)["t_eff_len_em"].size))]


def _solve(locus: str, mode: str, lean=None):
    """One captured locus through the shipped batch EM, from the coverage seed or from the seed leaning on ``lean``:
    each component's fractional count (transcripts, then the gDNA component) and the solver's per-locus profile."""
    fx = _fixture(locus)
    weights = fx["coverage_weights"].astype(np.float64)
    if lean == "gdna":
        weights /= LEAN
    elif lean is not None:
        weights[fx["t_indices"] == lean] *= LEAN
    fx["coverage_weights"] = weights.astype(np.float32)
    n_t = int(fx["t_eff_len_em"].size)
    geometry = TranscriptGeometry(
        effective_lengths=np.maximum(fx["t_eff_len_em"], 1.0),
        effective_lengths_em=fx["t_eff_len_em"],
    )
    config = EMConfig(
        mode=mode,
        assignment_mode="fractional",
        iterations=10000,
        convergence_delta=1e-9,
        n_threads=1,
    )
    rc = AbundanceEstimator(n_t, geometry=geometry, em_config=config)
    rc.unambig_counts[:] = fx["unambig"]
    gdna, _, _ = rc.run_batch_locus_em_partitioned(
        [
            (
                *(fx[k] for k in PARTITION_KEYS),
                # no unit of a fixture locus has an intergenic place
                np.full(len(fx["offsets"]) - 1, -np.inf, dtype=np.float32),
            )
        ],
        [np.arange(n_t, dtype=np.int32)],
        np.array([float(fx["gdna_prior"])]),
        rna_prior_count=np.array([float(fx["rna_prior"])]),
        gdna_eff_len=np.array([float(fx["gdna_eff_len"])]),
        emit_locus_stats=True,
    )
    return np.append(rc.em_counts.sum(axis=1), gdna), rc.locus_stats[0]


@pytest.mark.parametrize("locus,mode", PAIRS)
def test_every_start_reaches_one_answer(locus, mode):
    coverage, _ = _solve(locus, mode)
    for lean in _leans(locus):
        leaning, _ = _solve(locus, mode, lean)
        moved = 0.5 * np.abs(coverage - leaning).sum()
        assert moved < AGREE, (
            f"{locus} {mode}: the start leaning on {lean} ends {moved:.4g} fragments from the coverage seed's"
        )


@pytest.mark.parametrize("locus", LOCI)
def test_each_captured_locus_exercises_the_backtrack(locus):
    """The backtrack runs on every fixture: some extrapolation step is shortened rather than taken."""
    halvings = sum(
        _solve(locus, mode, lean)[1]["squarem_backtrack_count"]
        for mode in MODES
        for lean in (None, *_leans(locus))
    )
    assert halvings > 0
