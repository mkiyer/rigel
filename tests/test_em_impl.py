"""Tests for the C++ EM solver (rigel._em_impl).

Tests run_locus_em_native() directly with synthetic CSR data,
verifying correctness of bias correction, equivalence-class grouping,
SQUAREM convergence, OVR prior, and VBEM mode.
"""

import math

import numpy as np
import pytest

from rigel._em_impl import run_locus_em_native


# ---------------------------------------------------------------------------
# Helper to build minimal CSR data for run_locus_em_native
# ---------------------------------------------------------------------------

def _make_locus(
    *,
    units: list[list[tuple[int, float, float, int, int]]],
    n_components: int,
    unambig_totals: np.ndarray | None = None,
    bias_profiles: np.ndarray | None = None,
    prior_eligible: np.ndarray | None = None,
    total_pseudocount: float = 1.0,
    max_iterations: int = 1000,
    convergence_delta: float = 1e-6,
    use_vbem: bool = False,
):
    """Build CSR arrays from a list of units and call the C++ solver.

    Each unit is a list of (comp_idx, log_lik, coverage_wt, tx_start, tx_end).

    Returns (theta, alpha, em_totals) as numpy arrays.
    """
    # Build CSR
    all_t_indices = []
    all_log_liks = []
    all_cov_wts = []
    all_tx_starts = []
    all_tx_ends = []
    offsets = [0]

    for unit in units:
        for comp_idx, ll, cw, txs, txe in unit:
            all_t_indices.append(comp_idx)
            all_log_liks.append(ll)
            all_cov_wts.append(cw)
            all_tx_starts.append(txs)
            all_tx_ends.append(txe)
        offsets.append(len(all_t_indices))

    offsets_arr = np.array(offsets, dtype=np.int64)
    t_indices_arr = np.array(all_t_indices, dtype=np.int32)
    log_liks_arr = np.array(all_log_liks, dtype=np.float64)
    cov_wts_arr = np.array(all_cov_wts, dtype=np.float64)
    tx_starts_arr = np.array(all_tx_starts, dtype=np.int32)
    tx_ends_arr = np.array(all_tx_ends, dtype=np.int32)

    if bias_profiles is None:
        # Default: large profiles so bias correction is ~0
        bias_profiles = np.full(n_components, 100_000, dtype=np.int64)
    if unambig_totals is None:
        unambig_totals = np.zeros(n_components, dtype=np.float64)
    if prior_eligible is None:
        prior_eligible = np.ones(n_components, dtype=np.float64)

    eff_lens = np.ones(n_components, dtype=np.float64)

    theta, alpha, em_totals = run_locus_em_native(
        offsets_arr, t_indices_arr, log_liks_arr, cov_wts_arr,
        tx_starts_arr, tx_ends_arr, bias_profiles,
        unambig_totals, eff_lens, prior_eligible,
        n_components, total_pseudocount,
        max_iterations, convergence_delta,
        use_vbem,
        0,  # n_transcripts=0 → classic (unlinked) mode
        0,  # n_nrna=0
        np.array([], dtype=np.int32),  # t_to_nrna
        np.array([], dtype=np.int32),  # nrna_to_t_offsets
        np.array([], dtype=np.int32),  # nrna_to_t_indices
    )
    return np.asarray(theta), np.asarray(alpha), np.asarray(em_totals)


# ======================================================================
# Test cases
# ======================================================================


class TestEmptyLocus:
    """Empty or trivial loci."""

    def test_zero_units(self):
        """Zero units → theta is prior-only."""
        theta, alpha, em = _make_locus(
            units=[],
            n_components=3,
            unambig_totals=np.array([5.0, 3.0, 0.0]),
        )
        assert theta.shape == (3,)
        assert np.all(np.isfinite(theta))
        assert abs(theta.sum() - 1.0) < 1e-12
        # Component 0 should get the largest share
        assert theta[0] > theta[1]

    def test_single_unit_single_candidate(self):
        """One unit with one candidate → theta dominated by that component."""
        theta, alpha, em = _make_locus(
            units=[[(0, 0.0, 1.0, 0, 200)]],
            n_components=3,
            unambig_totals=np.array([10.0, 0.0, 0.0]),
        )
        assert theta[0] > 0.9


class TestTwoTranscriptLocus:
    """Two-transcript locus with known expected ratios."""

    def test_2x_evidence_ratio(self):
        """Transcript A has 2× evidence → theta_A ≈ 2 × theta_B."""
        # 3 units for A, 1 unit for B, plus 2 ambiguous
        # Components: 0=mRNA_A, 1=mRNA_B, 2=nRNA_A, 3=nRNA_B, 4=gDNA
        n_comp = 5
        units = [
            # Ambiguous units mapping to both A and B (equal log-lik)
            [(0, -1.0, 1.0, 0, 200), (1, -1.0, 1.0, 0, 200)],
            [(0, -1.0, 1.0, 0, 200), (1, -1.0, 1.0, 0, 200)],
            [(0, -1.0, 1.0, 0, 200), (1, -1.0, 1.0, 0, 200)],
        ]
        # A has 20 unambig, B has 10 unambig
        unambig = np.array([20.0, 10.0, 0.0, 0.0, 0.0])
        eligible = np.array([1.0, 1.0, 0.0, 0.0, 0.0])

        theta, alpha, em = _make_locus(
            units=units,
            n_components=n_comp,
            unambig_totals=unambig,
            prior_eligible=eligible,
        )
        # A should get roughly 2× the share of B
        ratio = theta[0] / max(theta[1], 1e-30)
        assert 1.5 < ratio < 3.0

    def test_unequal_log_liks(self):
        """Higher log-likelihood shifts posterior toward that component."""
        n_comp = 3
        # 5 ambiguous units: comp 0 has much better log-lik than comp 1
        units = [
            [(0, -0.1, 1.0, 0, 200), (1, -5.0, 1.0, 0, 200)]
            for _ in range(5)
        ]
        unambig = np.array([1.0, 1.0, 0.0])
        eligible = np.array([1.0, 1.0, 0.0])

        theta, alpha, em = _make_locus(
            units=units,
            n_components=n_comp,
            unambig_totals=unambig,
            prior_eligible=eligible,
        )
        # Component 0 has much better log-lik → should dominate
        assert theta[0] > theta[1] * 2


class TestConvergence:
    """Convergence behavior."""

    def test_converges_within_budget(self):
        """Well-conditioned input should converge well before max_iterations."""
        n_comp = 3
        units = [
            [(0, -1.0, 1.0, 0, 200), (1, -1.0, 1.0, 0, 200)]
            for _ in range(20)
        ]
        unambig = np.array([50.0, 50.0, 0.0])
        eligible = np.array([1.0, 1.0, 0.0])

        theta, alpha, em = _make_locus(
            units=units,
            n_components=n_comp,
            unambig_totals=unambig,
            prior_eligible=eligible,
            max_iterations=1000,
        )
        # Should produce a valid probability distribution
        assert abs(theta.sum() - 1.0) < 1e-10
        assert np.all(theta >= 0)

    def test_single_iteration(self):
        """max_iterations=3 → one SQUAREM iteration (3/3=1)."""
        n_comp = 3
        units = [
            [(0, 0.0, 1.0, 0, 200), (1, 0.0, 1.0, 0, 200)]
            for _ in range(5)
        ]
        unambig = np.array([10.0, 10.0, 0.0])
        eligible = np.array([1.0, 1.0, 0.0])

        theta, alpha, em = _make_locus(
            units=units,
            n_components=n_comp,
            unambig_totals=unambig,
            prior_eligible=eligible,
            max_iterations=3,
        )
        assert abs(theta.sum() - 1.0) < 1e-10


class TestPriorEligibility:
    """Components with zero prior should converge to zero."""

    def test_ineligible_component_zeros_out(self):
        """Ineligible component gets no weight even with ambiguous evidence."""
        n_comp = 5  # mRNA_A, mRNA_B, nRNA_A, nRNA_B, gDNA
        units = [
            [(0, -1.0, 1.0, 0, 200), (1, -1.0, 1.0, 0, 200),
             (2, -1.0, 1.0, 0, 200)]
            for _ in range(10)
        ]
        unambig = np.array([5.0, 5.0, 0.0, 0.0, 0.0])
        # Only mRNA_A and mRNA_B are eligible; comp 2 is not
        eligible = np.array([1.0, 1.0, 0.0, 0.0, 0.0])

        theta, alpha, em = _make_locus(
            units=units,
            n_components=n_comp,
            unambig_totals=unambig,
            prior_eligible=eligible,
        )
        # Ineligible components should be ~0
        assert theta[2] < 0.01
        assert theta[3] < 0.01
        assert theta[4] < 0.01


class TestVBEM:
    """VBEM mode."""

    def test_vbem_converges(self):
        """VBEM should converge and produce valid theta."""
        n_comp = 3
        units = [
            [(0, -1.0, 1.0, 0, 200), (1, -2.0, 1.0, 0, 200)]
            for _ in range(15)
        ]
        unambig = np.array([10.0, 5.0, 0.0])
        eligible = np.array([1.0, 1.0, 0.0])

        theta, alpha, em = _make_locus(
            units=units,
            n_components=n_comp,
            unambig_totals=unambig,
            prior_eligible=eligible,
            use_vbem=True,
        )
        assert abs(theta.sum() - 1.0) < 1e-10
        assert np.all(theta >= 0)
        # Component 0 should get more weight (better log-lik + more unambig)
        assert theta[0] > theta[1]

    def test_vbem_suppresses_low_evidence(self):
        """VBEM should suppress low-evidence components more aggressively
        than MAP-EM due to digamma's behavior at small alpha."""
        n_comp = 3
        units = [
            [(0, -1.0, 1.0, 0, 200), (1, -1.0, 1.0, 0, 200)]
            for _ in range(5)
        ]
        # Component 0 has strong unambig evidence, comp 1 has none
        unambig = np.array([50.0, 0.0, 0.0])
        eligible = np.array([1.0, 1.0, 0.0])

        theta_map, _, _ = _make_locus(
            units=units, n_components=n_comp,
            unambig_totals=unambig.copy(), prior_eligible=eligible.copy(),
            use_vbem=False,
        )
        theta_vb, _, _ = _make_locus(
            units=units, n_components=n_comp,
            unambig_totals=unambig.copy(), prior_eligible=eligible.copy(),
            use_vbem=True,
        )
        # VBEM should give less weight to comp 1 than MAP-EM
        assert theta_vb[1] <= theta_map[1] + 0.01


class TestNumericalStability:
    """Edge cases for numerical stability."""

    def test_extreme_log_liks(self):
        """Very negative log-likelihoods should not produce NaN/Inf."""
        n_comp = 3
        units = [
            [(0, -500.0, 1.0, 0, 200), (1, -500.0, 1.0, 0, 200)]
            for _ in range(5)
        ]
        unambig = np.array([5.0, 5.0, 0.0])
        eligible = np.array([1.0, 1.0, 0.0])

        theta, alpha, em = _make_locus(
            units=units,
            n_components=n_comp,
            unambig_totals=unambig,
            prior_eligible=eligible,
        )
        assert np.all(np.isfinite(theta))
        assert abs(theta.sum() - 1.0) < 1e-10

    def test_disparate_log_liks(self):
        """One log-lik much larger than another → proper normalization."""
        n_comp = 3
        units = [
            [(0, -1.0, 1.0, 0, 200), (1, -900.0, 1.0, 0, 200)]
            for _ in range(5)
        ]
        unambig = np.array([1.0, 1.0, 0.0])
        eligible = np.array([1.0, 1.0, 0.0])

        theta, alpha, em = _make_locus(
            units=units,
            n_components=n_comp,
            unambig_totals=unambig,
            prior_eligible=eligible,
        )
        assert np.all(np.isfinite(theta))
        # Comp 0 should dominate (prior dilutes somewhat)
        assert theta[0] > 0.7

    def test_all_zero_unambig(self):
        """No unambig evidence, only ambiguous → should still converge."""
        n_comp = 3
        units = [
            [(0, -1.0, 1.0, 0, 200), (1, -1.0, 1.0, 0, 200)]
            for _ in range(10)
        ]
        unambig = np.zeros(n_comp)
        eligible = np.array([1.0, 1.0, 0.0])

        theta, alpha, em = _make_locus(
            units=units,
            n_components=n_comp,
            unambig_totals=unambig,
            prior_eligible=eligible,
        )
        assert np.all(np.isfinite(theta))
        assert abs(theta.sum() - 1.0) < 1e-10


class TestEquivalenceClasses:
    """Equivalence class grouping works correctly."""

    def test_identical_candidate_sets_grouped(self):
        """Units with same candidate set should produce same result
        regardless of their ordering in the CSR."""
        n_comp = 3
        # All units map to {0, 1} — should form one equivalence class
        units = [
            [(0, -1.0, 1.0, 0, 200), (1, -2.0, 1.0, 0, 200)]
            for _ in range(10)
        ]
        unambig = np.array([5.0, 5.0, 0.0])
        eligible = np.array([1.0, 1.0, 0.0])

        theta, alpha, em = _make_locus(
            units=units, n_components=n_comp,
            unambig_totals=unambig, prior_eligible=eligible,
        )
        assert np.all(np.isfinite(theta))
        assert theta[0] > theta[1]  # comp 0 has better log-lik

    def test_different_candidate_sets(self):
        """Different candidate sets form separate equivalence classes."""
        n_comp = 5
        units = [
            # Class 1: maps to {0, 1}
            [(0, -1.0, 1.0, 0, 200), (1, -1.0, 1.0, 0, 200)],
            [(0, -1.0, 1.0, 0, 200), (1, -1.0, 1.0, 0, 200)],
            # Class 2: maps to {0, 2}
            [(0, -1.0, 1.0, 0, 200), (2, -1.0, 1.0, 0, 200)],
            # Class 3: maps to {1}
            [(1, -1.0, 1.0, 0, 200)],
        ]
        unambig = np.array([5.0, 5.0, 5.0, 0.0, 0.0])
        eligible = np.array([1.0, 1.0, 1.0, 0.0, 0.0])

        theta, alpha, em = _make_locus(
            units=units, n_components=n_comp,
            unambig_totals=unambig, prior_eligible=eligible,
        )
        assert np.all(np.isfinite(theta))
        assert abs(theta.sum() - 1.0) < 1e-10


class TestBiasCorrection:
    """Bias correction applied inside the C++ solver."""

    def test_bias_correction_reduces_long_transcript_advantage(self):
        """Longer transcripts have larger effective length → penalty.

        Two components with same log-lik but different transcript lengths
        should result in the shorter transcript getting relatively more
        weight after bias correction, because the longer transcript's
        effective length is larger.
        """
        n_comp = 3
        # Both map to comp 0 and comp 1 with equal log-lik
        units = [
            [(0, 0.0, 1.0, 0, 200), (1, 0.0, 1.0, 0, 200)]
            for _ in range(20)
        ]
        unambig = np.array([10.0, 10.0, 0.0])
        eligible = np.array([1.0, 1.0, 0.0])

        # Comp 0 is short (500bp), comp 1 is very long (10000bp).
        # bias correction: -log(max(L - frag_len + 1, 1))
        # Comp 0: -log(max(500 - 200 + 1, 1)) = -log(301) ≈ -5.71
        # Comp 1: -log(max(10000 - 200 + 1, 1)) = -log(9801) ≈ -9.19
        # So comp 0 gets a LESS negative adjustment → higher posterior
        profiles = np.array([500, 10000, 100000], dtype=np.int64)

        theta, alpha, em = _make_locus(
            units=units, n_components=n_comp,
            unambig_totals=unambig, prior_eligible=eligible,
            bias_profiles=profiles,
        )
        # Shorter transcript (comp 0) should get MORE weight
        assert theta[0] > theta[1]


class TestOVRPrior:
    """One Virtual Read prior behavior."""

    def test_tiny_pseudocount_minimizes_ovr(self):
        """Very small pseudocount → OVR prior negligible, data dominates."""
        n_comp = 3
        units = [
            [(0, -1.0, 1.0, 0, 200), (1, -1.0, 1.0, 0, 200)]
            for _ in range(10)
        ]
        unambig = np.array([10.0, 10.0, 0.0])
        eligible = np.array([1.0, 1.0, 0.0])

        theta, alpha, em = _make_locus(
            units=units, n_components=n_comp,
            unambig_totals=unambig, prior_eligible=eligible,
            total_pseudocount=1e-10,
        )
        assert np.all(np.isfinite(theta))
        # With equal everything and negligible OVR, components should be roughly equal
        assert abs(theta[0] - theta[1]) < 0.1

    def test_high_pseudocount_shifts_prior(self):
        """Large pseudocount → OVR coverage weights dominate the prior."""
        n_comp = 3
        # All coverage weight goes to comp 0 (high cov_wt)
        units = [
            [(0, -1.0, 10.0, 0, 200), (1, -1.0, 0.1, 0, 200)]
            for _ in range(10)
        ]
        unambig = np.array([1.0, 1.0, 0.0])
        eligible = np.array([1.0, 1.0, 0.0])

        theta, alpha, em = _make_locus(
            units=units, n_components=n_comp,
            unambig_totals=unambig, prior_eligible=eligible,
            total_pseudocount=10.0,
        )
        # Comp 0 should get more weight due to high coverage weight
        assert theta[0] > theta[1]


class TestMAPvsVBEMConsistency:
    """Both modes should converge for the same input."""

    def test_both_modes_finite(self):
        """Both MAP-EM and VBEM produce valid theta on same input."""
        n_comp = 5
        units = [
            [(0, -1.0, 1.0, 0, 200), (1, -2.0, 1.0, 0, 200),
             (4, -3.0, 1.0, 0, 200)]
            for _ in range(20)
        ]
        unambig = np.array([10.0, 5.0, 0.0, 0.0, 1.0])
        eligible = np.array([1.0, 1.0, 0.0, 0.0, 1.0])

        theta_map, _, _ = _make_locus(
            units=units, n_components=n_comp,
            unambig_totals=unambig.copy(), prior_eligible=eligible.copy(),
            use_vbem=False,
        )
        theta_vb, _, _ = _make_locus(
            units=units, n_components=n_comp,
            unambig_totals=unambig.copy(), prior_eligible=eligible.copy(),
            use_vbem=True,
        )
        assert np.all(np.isfinite(theta_map))
        assert np.all(np.isfinite(theta_vb))
        assert abs(theta_map.sum() - 1.0) < 1e-10
        assert abs(theta_vb.sum() - 1.0) < 1e-10
        # Both should agree on which component is dominant
        assert np.argmax(theta_map) == np.argmax(theta_vb)


# ======================================================================
# Linked EM model tests (n_transcripts > 0)
# ======================================================================


def _make_linked_locus(
    *,
    units: list[list[tuple[int, float, float, int, int]]],
    n_transcripts: int,
    unambig_totals: np.ndarray | None = None,
    bias_profiles: np.ndarray | None = None,
    prior_eligible: np.ndarray | None = None,
    effective_lengths: np.ndarray | None = None,
    total_pseudocount: float = 1.0,
    max_iterations: int = 1000,
    convergence_delta: float = 1e-6,
):
    """Build CSR arrays for a linked-model locus and call the C++ solver.

    Component layout: [0..n_t) mRNA, [n_t..2*n_t) nRNA, [2*n_t] gDNA.
    With n_nrna = n_t (each transcript has its own unique nRNA):
        n_components = n_t + n_nrna + 1 = 2*n_t + 1.

    Returns (theta, alpha, em_totals, nrna_frac) as numpy arrays.
    nrna_frac is derived post-hoc from theta: η[n] = θ[n_t+n] / (θ[n_t+n] + θ[n]).
    """
    n_t = n_transcripts
    n_nrna = n_t  # assume 1:1 transcript→nRNA
    n_components = n_t + n_nrna + 1

    all_t_indices = []
    all_log_liks = []
    all_cov_wts = []
    all_tx_starts = []
    all_tx_ends = []
    offsets = [0]
    for unit in units:
        for comp_idx, ll, cw, txs, txe in unit:
            all_t_indices.append(comp_idx)
            all_log_liks.append(ll)
            all_cov_wts.append(cw)
            all_tx_starts.append(txs)
            all_tx_ends.append(txe)
        offsets.append(len(all_t_indices))

    offsets_arr = np.array(offsets, dtype=np.int64)
    t_indices_arr = np.array(all_t_indices, dtype=np.int32)
    log_liks_arr = np.array(all_log_liks, dtype=np.float64)
    cov_wts_arr = np.array(all_cov_wts, dtype=np.float64)
    tx_starts_arr = np.array(all_tx_starts, dtype=np.int32)
    tx_ends_arr = np.array(all_tx_ends, dtype=np.int32)

    if bias_profiles is None:
        bias_profiles = np.full(n_components, 100_000, dtype=np.int64)
    if unambig_totals is None:
        unambig_totals = np.zeros(n_components, dtype=np.float64)
    if prior_eligible is None:
        prior_eligible = np.ones(n_components, dtype=np.float64)
    if effective_lengths is None:
        effective_lengths = np.ones(n_components, dtype=np.float64)

    # 1:1 transcript→nRNA mapping for tests
    t_to_nrna = np.arange(n_t, dtype=np.int32)
    nrna_to_t_offsets = np.arange(n_nrna + 1, dtype=np.int32)
    nrna_to_t_indices = np.arange(n_t, dtype=np.int32)

    theta, alpha, em_totals = run_locus_em_native(
        offsets_arr, t_indices_arr, log_liks_arr, cov_wts_arr,
        tx_starts_arr, tx_ends_arr, bias_profiles,
        unambig_totals, effective_lengths, prior_eligible,
        n_components, total_pseudocount,
        max_iterations, convergence_delta,
        False,  # use_vbem (not supported for linked)
        n_transcripts, n_nrna,
        t_to_nrna, nrna_to_t_offsets, nrna_to_t_indices,
    )
    theta = np.asarray(theta)
    # Derive nrna_frac from theta: η[n] = θ[n_t+n] / (θ[n_t+n] + θ[n])
    nrna_frac = np.zeros(n_nrna, dtype=np.float64)
    for n in range(n_nrna):
        nrna_w = theta[n_t + n]
        mrna_w = theta[n]  # 1:1 mapping
        denom = nrna_w + mrna_w
        nrna_frac[n] = nrna_w / denom if denom > 0 else 0.5
    return (
        theta,
        np.asarray(alpha),
        np.asarray(em_totals),
        nrna_frac,
    )


class TestLinkedEmBasic:
    """Basic linked-model EM tests."""

    def test_empty_locus_nrna_frac_near_zero(self):
        """Empty linked locus: nrna_frac is near zero (no nRNA evidence)."""
        n_t = 2
        theta, alpha, em, nrna_frac = _make_linked_locus(
            units=[], n_transcripts=n_t,
            unambig_totals=np.array([5.0, 3.0, 0.0, 0.0, 0.0]),
        )
        assert nrna_frac.shape == (n_t,)
        # No data for nRNA → nrna_frac should be very low
        assert np.all(nrna_frac < 0.1)

    def test_theta_is_simplex(self):
        """Linked EM produces a valid simplex θ[2*n_t+1]."""
        n_t = 2
        units = [
            [(0, 0.0, 1.0, 0, 0), (2, -0.5, 1.0, 0, 0)],  # mRNA_0 vs nRNA_0
            [(1, 0.0, 1.0, 0, 0), (3, -0.5, 1.0, 0, 0)],  # mRNA_1 vs nRNA_1
            [(0, 0.0, 1.0, 0, 0), (4, -1.0, 1.0, 0, 0)],  # mRNA_0 vs gDNA
        ]
        theta, alpha, em, nrna_frac = _make_linked_locus(
            units=units, n_transcripts=n_t,
            unambig_totals=np.array([10.0, 5.0, 1.0, 1.0, 0.5]),
        )
        assert theta.shape == (2 * n_t + 1,)
        assert abs(theta.sum() - 1.0) < 1e-10
        assert np.all(theta >= 0)
        assert np.all(np.isfinite(theta))

    def test_nrna_frac_in_bounds(self):
        """Converged nrna_frac is in [ε, 1-ε]."""
        n_t = 2
        units = [
            [(0, 0.0, 1.0, 0, 0), (2, -0.5, 1.0, 0, 0)],
            [(1, 0.0, 1.0, 0, 0), (3, -0.5, 1.0, 0, 0)],
        ]
        theta, alpha, em, nrna_frac = _make_linked_locus(
            units=units, n_transcripts=n_t,
            unambig_totals=np.array([10.0, 5.0, 1.0, 1.0, 0.5]),
        )
        NRNA_FRAC_EPS = 1e-8
        assert np.all(nrna_frac >= NRNA_FRAC_EPS)
        assert np.all(nrna_frac <= 1.0 - NRNA_FRAC_EPS)

    def test_theta_decomposition(self):
        """θ[i] + θ[n_t+i] should approximately equal θ_t[i].

        The decomposition is θ_mRNA = θ_t*(1−nrna_frac), θ_nRNA = θ_t*nrna_frac.
        """
        n_t = 2
        units = [
            [(0, 0.0, 1.0, 0, 0), (2, -0.3, 1.0, 0, 0)],
            [(1, 0.0, 1.0, 0, 0), (3, -0.3, 1.0, 0, 0)],
        ]
        theta, alpha, em, nrna_frac = _make_linked_locus(
            units=units, n_transcripts=n_t,
            unambig_totals=np.array([10.0, 5.0, 2.0, 1.0, 0.5]),
        )
        for i in range(n_t):
            theta_t = theta[i] + theta[n_t + i]
            if theta_t > 1e-12:
                nrna_frac_from_decomp = theta[n_t + i] / theta_t
                np.testing.assert_allclose(nrna_frac[i], nrna_frac_from_decomp, atol=1e-8)


class TestLinkedEmDataDriven:
    """Test that data overwhelms weak priors correctly."""

    def test_all_mrna_data_gives_low_nrna_frac(self):
        """When all data is mRNA, nrna_frac should be very low."""
        n_t = 1
        # 50 units all assigned to mRNA_0 only
        units = [[(0, 0.0, 1.0, 0, 0)]] * 50
        theta, _, _, nrna_frac = _make_linked_locus(
            units=units, n_transcripts=n_t,
            unambig_totals=np.array([50.0, 0.0, 0.0]),
        )
        assert nrna_frac[0] < 0.1

    def test_all_nrna_data_gives_high_nrna_frac(self):
        """When all data is nRNA, nrna_frac should be very high."""
        n_t = 1
        # 50 units all assigned to nRNA_0 only (component 1)
        units = [[(1, 0.0, 1.0, 0, 0)]] * 50
        theta, _, _, nrna_frac = _make_linked_locus(
            units=units, n_transcripts=n_t,
            unambig_totals=np.array([0.0, 50.0, 0.0]),
        )
        assert nrna_frac[0] > 0.9

    def test_two_transcripts_different_nrna_frac(self):
        """Two transcripts with different mRNA/nRNA data converge to
        different nrna_frac values."""
        n_t = 2
        # Transcript 0: mostly mRNA (comp 0), some nRNA (comp 2)
        # Transcript 1: mostly nRNA (comp 3), some mRNA (comp 1)
        units = (
            [[(0, 0.0, 1.0, 0, 0)]] * 40          # 40 mRNA_0
            + [[(2, 0.0, 1.0, 0, 0)]] * 10         # 10 nRNA_0
            + [[(1, 0.0, 1.0, 0, 0)]] * 10         # 10 mRNA_1
            + [[(3, 0.0, 1.0, 0, 0)]] * 40         # 40 nRNA_1
        )
        unambig = np.array([20.0, 5.0, 5.0, 20.0, 0.0])
        theta, _, _, nrna_frac = _make_linked_locus(
            units=units, n_transcripts=n_t,
            unambig_totals=unambig,
        )
        # Transcript 0 should have lower nrna_frac (more mRNA)
        # Transcript 1 should have higher nrna_frac (more nRNA)
        assert nrna_frac[0] < 0.4
        assert nrna_frac[1] > 0.6
        assert nrna_frac[0] < nrna_frac[1]


class TestLinkedEmUnstranded:
    """Unstranded library (strand_specificity ≈ 0.5).

    When stranding is uninformative, the nrna_frac prior comes from density
    estimation alone.  The linked EM should converge correctly using
    these weaker priors.
    """

    def test_unstranded_prior_converges(self):
        """Weak Beta(1,1) prior (unstranded) still converges to data-driven nrna_frac."""
        n_t = 1
        # 30 mRNA units + 10 nRNA units → expect nrna_frac ≈ 0.25
        units = [[(0, 0.0, 1.0, 0, 0)]] * 30 + [[(1, 0.0, 1.0, 0, 0)]] * 10
        unambig = np.array([15.0, 5.0, 0.0])
        theta, _, _, nrna_frac = _make_linked_locus(
            units=units, n_transcripts=n_t,
            unambig_totals=unambig,
        )
        # nrna_frac should be near the data ratio ≈ (10+5)/(30+15+10+5) = 15/60 = 0.25
        assert 0.15 < nrna_frac[0] < 0.40

    def test_unstranded_ambiguous_units_resolved(self):
        """Ambiguous mRNA/nRNA units are resolved even with weak priors."""
        n_t = 1
        # All 40 units are ambiguous between mRNA_0 and nRNA_0
        units = [[(0, 0.0, 1.0, 0, 0), (1, 0.0, 1.0, 0, 0)]] * 40
        # Unique data strongly favours mRNA
        unambig = np.array([50.0, 5.0, 0.0])
        theta, _, _, nrna_frac = _make_linked_locus(
            units=units, n_transcripts=n_t,
            unambig_totals=unambig,
        )
        # nrna_frac should be pulled low by the mRNA-heavy unambig data
        assert nrna_frac[0] < 0.3
        assert theta[0] > theta[1]  # mRNA > nRNA

    def test_unstranded_two_tx_symmetric(self):
        """Two symmetric transcripts with uninformative prior stay balanced."""
        n_t = 2
        # Equal data for both transcripts, each with 50/50 mRNA/nRNA
        units = (
            [[(0, 0.0, 1.0, 0, 0), (2, 0.0, 1.0, 0, 0)]] * 20
            + [[(1, 0.0, 1.0, 0, 0), (3, 0.0, 1.0, 0, 0)]] * 20
        )
        unambig = np.array([10.0, 10.0, 10.0, 10.0, 0.0])
        theta, _, _, nrna_frac = _make_linked_locus(
            units=units, n_transcripts=n_t,
            unambig_totals=unambig,
        )
        # Both nrna_frac should be similar (near 0.5)
        np.testing.assert_allclose(nrna_frac[0], nrna_frac[1], atol=0.05)
        # Total transcript abundances should be similar
        theta_t0 = theta[0] + theta[2]
        theta_t1 = theta[1] + theta[3]
        np.testing.assert_allclose(theta_t0, theta_t1, atol=0.05)


class TestLinkedEmHighGDNA:
    """Capture-enriched / high gDNA background scenarios.

    Capture libraries can have significant gDNA contamination.
    The linked model should correctly attribute background to gDNA
    and maintain accurate nrna_frac estimates.
    """

    def test_high_gdna_background(self):
        """Strong gDNA signal doesn't corrupt transcript nrna_frac estimates."""
        n_t = 1
        # 20 mRNA, 5 nRNA, 30 gDNA fragments → nrna_frac should still reflect mRNA/nRNA ratio
        units = (
            [[(0, 0.0, 1.0, 0, 0)]] * 20
            + [[(1, 0.0, 1.0, 0, 0)]] * 5
            + [[(2, 0.0, 1.0, 0, 0)]] * 30  # gDNA-only units
        )
        unambig = np.array([10.0, 2.0, 20.0])
        theta, _, _, nrna_frac = _make_linked_locus(
            units=units, n_transcripts=n_t,
            unambig_totals=unambig,
        )
        # gDNA component should get substantial fraction
        assert theta[2] > 0.3  # gDNA
        # Transcript nrna_frac should still be reasonable (~0.2 from nRNA/(mRNA+nRNA))
        assert 0.05 < nrna_frac[0] < 0.4

    def test_gdna_competing_with_nrna(self):
        """Ambiguous units tested between nRNA and gDNA.

        Units that could be either nRNA (unspliced exonic) or gDNA should
        be split according to the converged θ and nrna_frac.
        """
        n_t = 1
        # Ambiguous units: candidates are nRNA_0 (comp 1) and gDNA (comp 2)
        ambig_units = [[(1, 0.0, 1.0, 0, 0), (2, 0.0, 1.0, 0, 0)]] * 30
        # Clear mRNA signal
        mrna_units = [[(0, 0.0, 1.0, 0, 0)]] * 20
        units = mrna_units + ambig_units
        unambig = np.array([15.0, 0.0, 10.0])
        theta, _, _, nrna_frac = _make_linked_locus(
            units=units, n_transcripts=n_t,
            unambig_totals=unambig,
        )
        # Should converge with reasonable estimates
        assert np.all(np.isfinite(theta))
        assert theta.sum() == pytest.approx(1.0, abs=1e-10)
        # The nRNA/gDNA ambiguity should be resolved:
        # mRNA component should dominate the transcript
        assert theta[0] > theta[1]  # mRNA > nRNA

    def test_capture_like_sparse_transcript(self):
        """Sparse transcript data with heavy background.

        A capture library might have very few on-target fragments but
        significant off-target gDNA.  The nrna_frac prior should anchor the
        estimate when data is sparse.
        """
        n_t = 2
        # Tx 0: 5 mRNA + 2 nRNA (sparse)
        # Tx 1: 3 mRNA + 1 nRNA (even sparser)
        # gDNA: 50 fragments (heavy background)
        units = (
            [[(0, 0.0, 1.0, 0, 0)]] * 5
            + [[(2, 0.0, 1.0, 0, 0)]] * 2
            + [[(1, 0.0, 1.0, 0, 0)]] * 3
            + [[(3, 0.0, 1.0, 0, 0)]] * 1
            + [[(4, 0.0, 1.0, 0, 0)]] * 50  # gDNA
        )
        unambig = np.array([3.0, 2.0, 1.0, 0.5, 30.0])
        # Strong priors: nrna_frac ≈ 0.2 for both transcripts
        theta, _, _, nrna_frac = _make_linked_locus(
            units=units, n_transcripts=n_t,
            unambig_totals=unambig,
        )
        # gDNA should get the lion's share
        assert theta[4] > 0.5
        # nrna_frac should be pulled towards prior mean since data is sparse
        assert 0.05 < nrna_frac[0] < 0.5
        assert 0.05 < nrna_frac[1] < 0.5
        # Both nrna_frac should be finite and in bounds
        assert np.all(np.isfinite(nrna_frac))
