"""rigel.calibration.region_init — the strand PROTOCOL DECISION and the OWN-EVIDENCE predicate.

The message-free self-solve ("pass-0") deconvolves each slot's unspliced fragment mass into
``(f_pos, f_neg, f_g)`` before any message passing and records what the slot's OWN data says about that
split as ``tau_lam``, the slot's λ-axis composition evidence — the DATA's Fisher information and never a
prior's, from two sources: the STRAND deconvolution (a Beta-Binomial that is RANK-1, so it informs only
``p``; at a single-strand slot the tilt is structurally locked and the strand PINS ``f_g``, at an AMBIG
slot the tilt is free and the strand cancels out of ``f_g`` — the Schur marginal — so the strand term is
gated to single-strand slots; identically zero on an unstranded library, where the protocol decision
:func:`strand_discriminability` reads the spliced split as κ = ½ exactly) and the INTRON FACTORY (the
curvature of the density deconvolution's per-slot λ-factor). Both are the kernel's, per block inside the
solve (`native/solve_kernel.cpp`: ``strand_evidence``, ``factor_precision_row``; bound for the gates as
`native.transfer_rows.strand_evidence` / `factor_precision`), and the sweep's diagnostics capture
publishes ``tau_lam`` and the factory's part ``tau_fac`` per slot (`blocks.SweepCapture`).

What lives here is what the kernel is told and what the instruments read: the protocol decision — the
library's one verdict on whether its strand split is a witness at all — and the ONE predicate on
``tau_lam``, `has_own_composition_evidence`. Structural certainty is not recorded here: the one predicate
for a pure-gDNA object is `region_geometry.g1_locked`, which the instruments import.

Layer: LAYER 6. It imports nothing from the backbone that consumes it.
"""

from __future__ import annotations

import numpy as np
from scipy.special import betaln

from ..native import transfer_rows as _rows

__all__ = [
    "has_own_composition_evidence",
    "strand_discriminability",
]

#: the own-evidence predicate's guard — the kernel's, so the two homes are one
_EPS = float(_rows.OWN_EVIDENCE_EPS)


def has_own_composition_evidence(tau_lam) -> np.ndarray:
    """THE ONE DEFINITION of "this slot has own composition evidence", and it lives here so every
    consumer imports it instead of restating the number.

    ``tau_lam`` is the λ-axis Fisher precision summed over the sources (the kernel's, published by the
    sweep's capture); anything above the divide-by-zero guard is a live channel. That is the whole
    content — the predicate is read off the solver's behaviour, not chosen. The kernel's own liveness
    test for the transfer policy's claims is ``tau_lam > 0``, and its ``has_composition`` reads this
    same guard; the two agree wherever a positive ``tau_lam`` exceeds it.

    ⛔ It is NOT a resolving-power test and must not become one. ``τ`` is continuous across the
    interesting region, so a floor on it is a tuned constant. On an unstranded library the strand arm
    is exactly zero — not by a floor on ``τ`` but because :func:`strand_discriminability` decides the
    PROTOCOL on the spliced 2×2 and reads κ = ½ exactly there (a sampling excursion of κ̂ is not a
    protocol). The consumer's defence against a weak live channel is a FIXED-DENOMINATOR score, not a
    tighter bound here (``solvability_audit.summarise``'s ``all_mwae`` / ``abs_err``, gated in
    ``test_solvability_audit.py``).

    The home is production, rather than each instrument restating the constant beside a comment saying
    it must match the solver, because the predicate is a production concept and ``scripts/`` is
    deliberately not importable.
    """
    return np.asarray(tau_lam, np.float64) > _EPS


def strand_discriminability(kappa, n_rna_obs) -> float:
    """The library's strand DISCRIMINABILITY: ``4(κ−½)²``, the strand Fisher information's library-level
    factor (the kernel's ``strand_evidence``), where the PROTOCOL preserves strand, and exactly 0 where it
    does not. The same for every slot: a positive value is the one condition under which any counted
    single-strand slot has a live strand channel, so it is also what the message layer reads as "is the
    strand split a witness at all" (`messages.ChainView.strand_live`). One definition, used by both;
    gated in ``tests/calibration/test_region_init.py``.

    Whether a protocol preserves strand is a decision between two hypotheses on the spliced 2×2 the
    strand fit read (`strand_balance.fit_strand_balance`: ``n_same`` sense reads of ``N``). H0, the
    unstranded protocol: read 1's strand is independent of the transcript's, so every sj reads ½ and the
    pooled split is ``Binomial(N, ½)`` with no free parameter — κ = ½ EXACTLY, not approximately. H1: κ
    free under the fit's own ``Beta(1, 1)`` prior. The Bayes factor is closed-form,

        ln BF₁₀ = N·ln 2 + ln B(a, b),   a = n_same + 1 = κ̂·(N + 2),   b = n_opp + 1 = (1 − κ̂)·(N + 2),

    and the channel is live iff ``BF₁₀ > 1`` (equal prior odds; no constant). For large ``N``,
    ``ln BF₁₀ ≈ ½·[z² − ln(2N/π)]`` with ``z = (κ̂ − ½)/√(¼/N)``: the free parameter's Occam penalty
    grows with the sample, so the excursion a sampling fluctuation must clear to be read as a protocol
    grows as ``√ln N`` — 3.8σ at three million spliced fragments — while a stranded library clears it by
    10⁴–10⁶ nats. The form this replaced, an unbiased estimate of ``(κ−½)²`` floored at zero, was
    positive on 32 % of genuinely unstranded libraries (a χ²₁ above its mean): a coin toss, not a
    deadband. gDNA enters nowhere — its strand mean is ½ by symmetry and needs no observation, and a
    term in its count had switched every gDNA-free library's channel off. With no spliced observation
    the two hypotheses have equal marginal likelihood and the channel is dead."""
    n = float(n_rna_obs)
    if not n > 0.0:
        return 0.0
    a = kappa * (n + 2.0)
    b = (1.0 - kappa) * (n + 2.0)
    ln_bf = n * np.log(2.0) + float(betaln(a, b))
    return 4.0 * (kappa - 0.5) ** 2 if ln_bf > 0.0 else 0.0
