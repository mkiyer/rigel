"""The strand protocol decision: whether the library's strand split is a witness at all.

The pre-sweep self-solve records each slot's own composition evidence ``tau_lam`` from the binomial strand
likelihood, computed per block by the kernel's ``strand_evidence`` (``native/solve_kernel.cpp``) and scaled by the
library-level factor defined here. The strand count informs only the mixed sense rate ``p``, so it pins ``f_g`` at a
single-strand slot and cancels out of ``f_g`` at an ambiguous one; the term is therefore gated to single-strand slots
and is zero on an unstranded library. Structural certainty of a pure-gDNA object is a separate predicate,
``region_geometry.g1_locked``. Layer 6; it imports nothing from the backbone that consumes it.
"""

from __future__ import annotations

import numpy as np
from scipy.special import betaln

__all__ = [
    "strand_discriminability",
]


def strand_discriminability(kappa, n_rna_obs) -> float:
    """The library's strand discriminability: ``4(κ−½)²`` if the protocol preserves strand, else exactly 0.

    It is the library-level factor of the strand Fisher information (the kernel's ``strand_evidence``) and the
    message layer's test of whether the strand split is a witness (``messages.ChainView.strand_live``).

    The decision is a Bayes factor on the spliced sense counts (``n_same`` of ``N``, from
    :func:`strand_balance.fit_strand_balance`). H0, unstranded: read 1's strand is independent of the transcript's,
    so κ = ½ exactly and the split is ``Binomial(N, ½)``. H1: κ free under the fit's ``Beta(1, 1)`` prior. Then

        ln BF₁₀ = N·ln 2 + ln B(a, b),   a = n_same + 1 = κ̂·(N + 2),   b = n_opp + 1 = (1 − κ̂)·(N + 2),

    and the channel is live iff ``BF₁₀ > 1`` (equal prior odds, so no threshold constant). For large ``N``,
    ``ln BF₁₀ ≈ ½·[z² − ln(2N/π)]`` with ``z = (κ̂ − ½)/√(¼/N)``: the free parameter's Occam penalty grows with the
    sample, so a sampling fluctuation of an unstranded library is not read as a protocol. gDNA must stay out of this
    test: its strand mean is ½ by symmetry, and a gDNA term would switch off a gDNA-free library's channel. With no
    spliced observation the two hypotheses tie and the channel is dead.

    Gate: ``tests/calibration/test_region_init.py``.
    """
    n = float(n_rna_obs)
    if not n > 0.0:
        return 0.0
    a = kappa * (n + 2.0)
    b = (1.0 - kappa) * (n + 2.0)
    ln_bf = n * np.log(2.0) + float(betaln(a, b))
    return 4.0 * (kappa - 0.5) ** 2 if ln_bf > 0.0 else 0.0
