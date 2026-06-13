"""Isolated unit tests for the simplex per-node solver (``calibration.simplex.solve_node``).

Validates the rev-2 frame (``docs/calibration/propagation_simplex_plan.md``) at the per-node level,
before the graph/worklist are wired: the pie partition is correct, over-subtraction is structurally
impossible (the 4030271-type catastrophe cannot recur), single-strand collapses to a strand result
(no-regression), the sided spliced flux lifts RNA, and mass is conserved.
"""

import numpy as np

from rigel.calibration.signature import TS_AMBIG, TS_NEG, TS_NONE, TS_POS
from rigel.calibration.simplex import init_from_signature, region_adjacency, solve_node


def _f(sol, i=0):
    return float(sol.f_rna_pos[i]), float(sol.f_rna_neg[i]), float(sol.f_g[i])


def test_mass_conservation():
    """The pie always sums to 1 on any active node, whatever the evidence."""
    sol = solve_node(
        np.array([50.0, 0.0, 30.0]), np.array([10.0, 0.0, 30.0]), kappa=0.9,
        spliced_pos=np.array([5.0, 0.0, 0.0]), count_gdna_frac=np.array([0.3, 0.5, 0.8]),
        count_precision=20.0, allow_pos=True, allow_neg=np.array([True, True, True]),
    )
    s = sol.f_rna_pos + sol.f_rna_neg + sol.f_g
    assert np.allclose(s[[0, 2]], 1.0, atol=1e-6)  # active nodes
    assert s[1] == 0.0  # empty node → all zero


def test_single_strand_pure_sense_is_rna():
    """A stranded POS node read almost entirely in sense (κ high) → nearly all RNA (f_g≈0)."""
    sol = solve_node(
        100.0, 1.0, kappa=0.99, allow_pos=True, allow_neg=False, gdna_prior_count=0.5,
    )
    f_pos, f_neg, f_g = _f(sol)
    assert f_neg == 0.0  # disallowed strand pinned
    assert f_g < 0.1
    assert f_pos > 0.9


def test_single_strand_balanced_is_gdna():
    """A stranded POS node read 50/50 (κ high) → balanced means gDNA (f_g≈1), not RNA."""
    sol = solve_node(
        50.0, 50.0, kappa=0.99, allow_pos=True, allow_neg=False, gdna_prior_count=0.5,
    )
    _f_pos, f_neg, f_g = _f(sol)
    assert f_neg == 0.0
    assert f_g > 0.85


def test_ambig_balanced_no_oversubtraction():
    """The 4030271-type catastrophe cannot recur. Balanced AMBIG reads (no spliced, no count) are
    genuinely ambiguous — t=0 fits gDNA (p=½) AND a balanced bidirectional RNA mix (f₊=f₋ ⇒ t=0)
    equally, so strand alone leaves f_g WIDE, leaning gDNA via the weak prior. The point is the
    simplex CANNOT drive f_g→0 by spurious RNA (the old count-cascade hit f_g=0.00 on a 94%-gDNA
    region); the count evidence is what resolves the ambiguity (next test)."""
    sol = solve_node(
        200.0, 200.0, kappa=0.99, allow_pos=True, allow_neg=True, gdna_prior_count=0.5,
    )
    _f_pos, _f_neg, f_g = _f(sol)
    assert f_g > 0.4  # never zeroed; leans gDNA but honestly wide


def test_ambig_balanced_resolved_by_count():
    """What strand can't resolve (balanced AMBIG = gDNA vs bidirectional RNA), the count does: a high
    count_gdna_frac drives the balanced AMBIG node to gDNA — the real 4030271 resolution path."""
    sol = solve_node(
        200.0, 200.0, kappa=0.99, count_gdna_frac=0.94, count_precision=150.0,
        allow_pos=True, allow_neg=True, gdna_prior_count=0.5,
    )
    _f_pos, _f_neg, f_g = _f(sol)
    assert f_g > 0.8


def test_ambig_strand_resolved_rna():
    """An AMBIG node read strongly +-skewed with κ high → +-strand RNA (f_pos high, f_g low)."""
    sol = solve_node(
        180.0, 20.0, kappa=0.99, allow_pos=True, allow_neg=True, gdna_prior_count=0.5,
    )
    f_pos, f_neg, f_g = _f(sol)
    assert f_pos > f_neg
    assert f_g < 0.3
    assert f_pos > 0.6


def test_unstranded_no_evidence_leans_gdna():
    """Unstranded (κ=½), no spliced, no count: RNA-vs-gDNA is unidentifiable, so the node rests
    above the uninformative centroid (1/3) toward gDNA — the honest 'unknown' the weak gDNA prior
    biases. (A derived, stronger gDNA prior sharpens this later.)"""
    sol = solve_node(
        100.0, 100.0, kappa=0.5, allow_pos=True, allow_neg=True, gdna_prior_count=0.5,
    )
    _f_pos, _f_neg, f_g = _f(sol)
    assert f_g > 1.0 / 3.0


def test_spliced_lower_bound_lifts_rna_when_strand_blind():
    """With no strand info (κ=½), a strong +-sided spliced flux is a lower bound on +-strand RNA →
    lifts f_pos and pushes f_g down (the sided mature evidence)."""
    base = solve_node(100.0, 100.0, kappa=0.5, allow_pos=True, allow_neg=True, gdna_prior_count=0.5)
    lifted = solve_node(
        100.0, 100.0, kappa=0.5, spliced_pos=80.0, allow_pos=True, allow_neg=True, gdna_prior_count=0.5,
    )
    assert float(lifted.f_rna_pos[0]) > float(base.f_rna_pos[0])
    assert float(lifted.f_rna_pos[0]) > float(lifted.f_rna_neg[0])
    assert float(lifted.f_g[0]) < float(base.f_g[0])


def test_count_evidence_pulls_fg():
    """The count term pulls f_g toward count_gdna_frac (here under κ=½ where strand is silent)."""
    lo = solve_node(100.0, 100.0, kappa=0.5, count_gdna_frac=0.1, count_precision=200.0,
                    allow_pos=True, allow_neg=True)
    hi = solve_node(100.0, 100.0, kappa=0.5, count_gdna_frac=0.9, count_precision=200.0,
                    allow_pos=True, allow_neg=True)
    assert float(lo.f_g[0]) < 0.4
    assert float(hi.f_g[0]) > 0.6


# --- increment 2: init-from-signature + chain adjacency -----------------------------------------


def test_init_from_signature_vectors():
    """The plan §2 table: NONE→fully-specified seed (f_g=1); POS/NEG→one strand active, strand-seedable;
    AMBIG→both active, not strand-seedable."""
    ts = np.array([TS_NONE, TS_POS, TS_NEG, TS_AMBIG])
    init = init_from_signature(ts)
    assert list(init.allow_pos) == [False, True, False, True]
    assert list(init.allow_neg) == [False, False, True, True]
    assert list(init.intergenic_seed) == [True, False, False, False]
    assert init.f_g_fixed[0] == 1.0 and np.all(np.isnan(init.f_g_fixed[1:]))
    assert list(init.strand_seedable) == [False, True, True, False]


def test_region_adjacency_chain_breaks_at_reference():
    """prev/next follow genomic order within a reference and are −1 at reference edges."""
    ref = np.array([0, 0, 0, 1, 1])  # two references: a 3-region chain and a 2-region chain
    prev, nxt = region_adjacency(ref)
    assert list(prev) == [-1, 0, 1, -1, 3]
    assert list(nxt) == [1, 2, -1, 4, -1]
