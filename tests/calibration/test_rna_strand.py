"""The strand likelihood's symmetry between the two components' overdispersions.

Both overdispersions are 0 in production (binomial by policy), but the strand term keeps its od arguments so an od can
return, and these gate what it must then do: with both components Beta-Binomial at one od, an unstranded region
(κ = ½) is uninformative and a balanced region grows more gDNA-like as the library becomes more stranded. An od on
the gDNA side alone breaks that — it pulls balanced unstranded regions toward RNA, a composition claim manufactured
out of a nuisance parameter — which is why the two components share one value.
"""

from __future__ import annotations

import numpy as np
import pytest

from _psi_reference import strand_loglik


def _od_for_beta(a: float) -> float:
    """The intra-class correlation of a symmetric ``Beta(a, a)``: ``1/(2a + 1)``."""
    return 1.0 / (2.0 * a + 1.0)


def _decoded_gdna_frac(sense, antisense, kappa, *, gdna_od, rna_od, n_grid=4000):
    """Posterior median gDNA fraction of one region under a FLAT count prior (strand-only deconv).

    Mirrors the strand module (the per-region strand branch): a weak prior ×
    the strand likelihood, isolating the strand likelihood's effect on the deconvolution.
    """
    grid = np.linspace(1e-6, 1.0 - 1e-6, n_grid)
    ll = strand_loglik(
        grid,
        sense,
        antisense,
        kappa,
        gdna_strand_overdispersion=gdna_od,
        rna_strand_overdispersion=rna_od,
    )
    w = np.exp(ll - ll.max())
    w /= w.sum()
    return float(np.interp(0.5, np.cumsum(w), grid))


# --------------------------------------------------------------------------- deconv symmetry


def test_unstranded_is_uninformative_with_symmetric_overdispersion():
    """κ = ½, equal gDNA/RNA overdispersion ⇒ a balanced region deconvolves to gdna_frac ≈ ½ (flat)."""
    od = _od_for_beta(3.0)
    frac = _decoded_gdna_frac(50, 50, 0.5, gdna_od=od, rna_od=od)
    assert frac == pytest.approx(0.5, abs=0.02)


def test_asymmetric_overdispersion_biases_unstranded_toward_rna():
    """An asymmetric pair (gDNA Beta-Binomial, RNA Binomial) spuriously pulls a balanced unstranded
    region toward RNA. Symmetric overdispersion removes the pull, which is why the RNA side is
    modelled at all."""
    od = _od_for_beta(3.0)
    asym = _decoded_gdna_frac(50, 50, 0.5, gdna_od=od, rna_od=0.0)
    symm = _decoded_gdna_frac(50, 50, 0.5, gdna_od=od, rna_od=od)
    assert asym < 0.4  # materially pulled toward RNA
    assert symm == pytest.approx(0.5, abs=0.02)  # pull removed


def test_graded_information_balanced_region_more_gdna_as_library_stranded():
    """A balanced (50/50) region reads more gDNA-like as κ rises from ½ (unstranded) toward 1.

    At κ = ½ it is uninformative (½); as the library becomes more stranded, a *symmetric* split
    looks increasingly like the symmetric gDNA component, so gdna_frac rises monotonically — the
    'unstranded → weakly → strongly stranded' information gradient.
    """
    od = _od_for_beta(3.0)
    kappas = [0.5, 0.6, 0.7, 0.8, 0.9, 0.99]
    fracs = [_decoded_gdna_frac(50, 50, k, gdna_od=od, rna_od=od) for k in kappas]
    assert fracs[0] == pytest.approx(0.5, abs=0.02)
    assert all(b >= a - 1e-6 for a, b in zip(fracs, fracs[1:]))  # monotone non-decreasing
    assert fracs[-1] > fracs[0] + 0.2  # materially more gDNA at high κ


def test_rna_overdispersion_zero_recovers_binomial_decode():
    """rna_strand_overdispersion = 0 ⇒ strand_loglik is byte-identical to the gDNA-only formula."""
    grid = np.linspace(1e-6, 1 - 1e-6, 200)
    base = strand_loglik(grid, 30, 10, 0.9, gdna_strand_overdispersion=0.1)
    with_rna0 = strand_loglik(
        grid, 30, 10, 0.9, gdna_strand_overdispersion=0.1, rna_strand_overdispersion=0.0
    )
    np.testing.assert_array_equal(base, with_rna0)
