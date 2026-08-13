"""Strand likelihood — the two-component gDNA/RNA strand log-likelihood of one region.

⭐⭐ **THIS IS AN EXECUTABLE REFERENCE, NOT A PRODUCTION PATH — and saying so is the point.** Nothing in
``src/`` calls it. What ships is the THREE-component generalization in
``simplex_logodds._mixture_strand_loglik``, and this two-component special case is what that form is gated
against: ``tests/calibration/test_strand_likelihood_reference.py`` asserts the collapse when one RNA strand
is dead, so a drift between them fires a test instead of a comment being quietly wrong.

⛔ This docstring used to say "Used by the per-region strand module (`strand_deconv`)", and `strand_deconv`
does not import it. It read as dead code for that reason, and the guard the source claimed elsewhere
("collapses to ``strand_loglik`` exactly — the no-regression guard") **did not exist until 2026-08-07**.

A region's unspliced fragments are a mix of **gDNA** and **RNA**, both **Beta-Binomial overdispersed**:

* gDNA is unstranded (oriented-sense rate ½) with intra-class correlation
  ``gdna_strand_overdispersion`` (fitted in :mod:`gdna_strand`).
* RNA is stranded at ``rna_sense_frac`` (κ) with intra-class correlation
  ``rna_strand_overdispersion`` (fitted from the per-sj SJ strand table in :mod:`gdna_strand`).

Both components carry a fitted overdispersion with the **same default prior**, so they collapse to
the same distribution under sparse data and an unstranded region (κ = ½, equal overdispersions) is
**uninformative** — the symmetry that an earlier gDNA-only Beta-Binomial broke (the ``−½·log var``
term preferred the lower-variance component, pulling balanced regions toward RNA).

Given the observed (sense, antisense) split, this returns the strand log-likelihood over a grid
of candidate gDNA fractions ``gdna_frac``, which the per-region strand module multiplies into the
count prior. Only POS/NEG regions carry a valid sense split; AMBIG (both strands) and intergenic
(no transcript) regions are strand-uninformative and left to the count clue.
"""

from __future__ import annotations

import numpy as np


def strand_loglik(
    gdna_frac: np.ndarray,
    sense: float,
    antisense: float,
    rna_sense_frac: float,
    *,
    gdna_strand_overdispersion: float = 0.0,
    rna_strand_overdispersion: float = 0.0,
) -> np.ndarray:
    """Two-component gDNA/RNA strand log-likelihood of one region over a ``gdna_frac`` grid (§6).

    Of ``N = sense + antisense`` discrete unspliced fragments, a fraction ``gdna_frac`` are gDNA
    (oriented-sense rate ½, intra-class correlation ``gdna_strand_overdispersion``) and
    ``1 − gdna_frac`` are RNA (oriented-sense rate ``rna_sense_frac``, intra-class correlation
    ``rna_strand_overdispersion``). The mixture sense count has mean ``N·p`` and a variance in
    three parts — the Binomial mixture variance plus the **excess variance each component's shared
    per-region sense rate contributes**, scaled by that component's own ``μ_c(1−μ_c)``: the
    ``N·gdna_frac`` gDNA fragments (mean ½) add ``(N·gdna_frac)²·¼·gdna_strand_overdispersion`` and
    the ``N·(1−gdna_frac)`` RNA fragments (mean κ) add
    ``(N·(1−gdna_frac))²·κ(1−κ)·rna_strand_overdispersion``. Normal moment approximation::

        p   = ½·gdna_frac + rna_sense_frac·(1 − gdna_frac);  mean = N·p
        var = N·p·(1 − p)
            + (N·gdna_frac)²·¼·gdna_strand_overdispersion
            + (N·(1 − gdna_frac))²·κ(1 − κ)·rna_strand_overdispersion
        loglik(gdna_frac) = −½·(sense − mean)² / var  −  ½·log(var)

    Limits: ``gdna_frac → 1`` ⇒ Beta-Binomial(½, od_gdna); ``gdna_frac → 0`` ⇒ Beta-Binomial(κ,
    od_rna); both ``od → 0`` ⇒ the Binomial mixture (recovers the pre-overdispersion behaviour
    exactly). **Symmetry:** at ``κ = ½`` the RNA scale ``κ(1−κ)`` equals the gDNA ¼, so with
    ``od_gdna = od_rna`` the variance is flat in ``gdna_frac`` (the means coincide and only the
    ``g² + (1−g)²`` scaling depends on ``gdna_frac``) — an unstranded region is uninformative. Each
    component's excess variance uses its own mean ``μ_c(1−μ_c)``, consistent with the moment fit in
    :mod:`gdna_strand`; the normal-moment vs exact-mixture discrepancy is small and ~constant in
    ``N``.
    """
    # Fully elementwise in (sense, antisense, gdna_frac), so it broadcasts: scalar (sense, antisense)
    # with a 1-D grid returns one region's curve; column (sense, antisense) of shape (K, 1) with a row
    # grid of shape (1, n_grid) returns the whole (K, n_grid) batch at once (used by strand_deconv).
    n = sense + antisense
    p = 0.5 * gdna_frac + rna_sense_frac * (1.0 - gdna_frac)
    mean = n * p
    rna_var_scale = rna_sense_frac * (1.0 - rna_sense_frac)  # κ(1−κ); the RNA component's μ(1−μ)
    var = (
        n * p * (1.0 - p)
        + (n * gdna_frac) ** 2 * 0.25 * gdna_strand_overdispersion
        + (n * (1.0 - gdna_frac)) ** 2 * rna_var_scale * rna_strand_overdispersion
    )
    var = np.maximum(var, 1.0e-9)
    return -0.5 * (sense - mean) ** 2 / var - 0.5 * np.log(var)


__all__ = ["strand_loglik"]
