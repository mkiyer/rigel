"""The log-density 1-D/2-D per-region solver — the single production per-region solve driving
``sweep.solve_chain`` (the memory-prohibitive 2-simplex lattice it replaced is retired).

The latent magnitude dof is the
gDNA-vs-RNA **log-odds** ``λ = logit(f_g) = log ρ_g − log ρ_rna`` (log-odds bounds the 5–6-decade ρ_g
range and resolves both ``f_g→0`` and ``f_g→1`` vertices, which the uniform linear lattice cannot). We
grid ``λ`` on a FIXED ``[−L, L]`` window (no region-adaptivity) and read out the linear fraction
``f_g = σ(λ)``. ``O(m·K)`` per region (vs the lattice's ``O(m·K²)`` 2-simplex), so it is genome-scale
tractable.

The ``ψ`` integrand is ``strand + (gDNA arm) + (RNA arm) + the imputation messages``, where each **arm** is
that component group's fitted log-rate prior when we have one, else the **Jeffreys reference**
``+½·log f`` (``_JEFFREYS_REF``). Derivation, review, and the resolved design:
 (§10 is authoritative).

Three facts that determine this file's shape:

1. **Omitting a component's term is not "no prior"** — the grid's own measure supplies one. A bare uniform-λ
   grid IS Haldane per component ⇒ Beta(0,0) on the composition: improper at BOTH vertices, a vertex
   amplifier. There is no third option. (This is what shipped before; it hid while ψ stayed symmetric.)
2. **The composition is a TWO-GROUP split on the λ axis** — gDNA vs RNA-**total** — which is what calibration
   models. The per-strand tilt is a NUISANCE parameter. On the two-group axis the measure residual is
   **exactly zero**: each ``logP`` is a density in LOG-rate, so its linear-rate conversion ``−log ρ`` cancels
   ``log σ'(λ)`` exactly, once per group. **No Jacobian is written.**
3. **The tilt is gridded as ``θ = arcsin(τ)``, NOT ``τ``.** The Berger–Bernardo reference prior for this model
   (``f_g`` of interest, ``τ`` nuisance — the two are information-ORTHOGONAL, ``I_{f_g,τ} = 0`` exactly) has a
   ``(1−τ²)^{−½}`` tilt conditional. Under ``θ = arcsin(τ)`` the Jacobian ``|dτ/dθ| = cos θ = (1−τ²)^{½}``
   cancels it **identically**, so the tilt term is **exactly 0** and the reference collapses to ONE expression
   for both region classes: ``ψ_ref = ½·log f_g + ½·log(1−f_g)``. No class branch, no endpoint singularity, no
   quadrature weights. θ is to the tilt what λ is to ``f_g``: the coordinate the geometry asks for.
   *(This vanishing is a property of the BB reference specifically — a Dirichlet(½,¼,¼) reference would leave
   a residual ``−¼·log(1−τ²)``.)*

There is NO spliced term: ``mass_spliced`` is consumed only by the returned ``rna_mass``, never by ψ. That is
correct — at a sj mature RNA *splices*, so the unspliced crossing mass is gDNA + nascent, a channel
genuinely disjoint from the (directly observed, already-pure-RNA) spliced mass.

Single-strand regions (exactly one of ``allow_pos`` / ``allow_neg``) are an exact 1-D solve over ``λ``; AMBIG
regions (both set) marginalize the tilt on a 2-D ``(λ, θ)`` grid (``_solve_ambig_logodds``).
``_solve_regions_logodds_all`` dispatches between the two. Structurally RNA-free regions (neither strand live —
intergenic / TSS / TES) have no composition DOF and never reach either solver: ``sweep.solve_chain`` gates
them out via ``solvable``, so no reference is applied to a region whose composition is known structurally.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import expit, log_expit

from .region_chain import RegionDeconv

# Public surface consumed by sweep / messages / region_geometry. The remaining private helpers stay importable
# for tests but are not part of the module's external API.
__all__ = [
    "CompositionPriors",
    "_logodds_grid",
    "_solve_regions_logodds_all",
]

_EPS = 1.0e-9

# The reference exponent for an UNFITTED component group, as a density in LOG-rate.
#
# DERIVED, not tuned, by two agreeing routes: (a) Jeffreys for a Poisson rate — `g ~ Poisson(ρE)` ⇒
# `I(ρ) ∝ 1/ρ` ⇒ `p(ρ) ∝ ρ^(−½)`, which as a LOG-rate density is `ρ^(+½)` ⇒ `+½·log f`; (b) the
# Berger–Bernardo reference prior for the composition with `f_g` of interest and the tilt as nuisance, whose
# `f_g` marginal is Beta(½,½) — the SAME `+½·log f_g + ½·log(1−f_g)`.
#
# Its ONLY job is to make ψ PROPER (Beta(½,½) integrates; Beta(0,0) does not). A fitted `logP_g`/`logP_r`
# is ADDED to it, never substituted for it — the reference is the MEASURE ψ is written against, not an
# information claim to be superseded. ⛔ This comment said "must be REPLACED … never added on top" until
# 2026-08-15, which contradicted `_gdna_arm`'s own code and the argument in its docstring.
#
# ⚠ A DECLARED CHOICE, not forced by the likelihood: the observed-data Fisher information for f_g is
# `∝ n(½−κ)²` = EXACTLY 0 on unstranded libraries, where the strand term is bit-flat and the posterior simply
# IS this reference. Licensed as the "structural Jeffreys" prior; §10.5
# records the known cost — it forbids the simplex vertices, where some truth genuinely lives.
_JEFFREYS_REF = 0.5

#: ⛔⛔ **THE CERTIFIED-RNA CLAIM IS A LOWER BOUND, AND ψ HAS ALWAYS APPLIED IT AS A TWO-SIDED GAUSSIAN.**
#: the message policy states the premise in its own words — ``rho_R(exon) >= rho_nu(B) + rho_mu(B)``, because
#: the exon may also hold molecules that never touch that boundary — "and it uses it as an equality". Three
#: operators price that inequality as a VARIANCE and none prices it as a DIRECTION, which is `TRAPS.md`
#: TRAPS: a-variance-cannot-fix-a-bias: a variance cannot move a mode toward truth.
#: ⭐ ``True`` selects the one-sided form: no penalty when the destination holds MORE RNA than the bound,
#: full penalty when it holds less. ⛔ ``False`` is today's behaviour and is BYTE-IDENTICAL by
#: construction — :func:`_rna_residual` then returns its input difference unmodified.
#: Set by ``ladder_arm_ab.py --arm onesided_rna``.
ONE_SIDED_RNA = [False]


def _rna_residual(log_f, mode, one_sided=None):
    """The residual the RNA imputed-message penalty is built on — ONE HOME for both ψ paths.

    Returns ``log_f - mode``, or its negative part when :data:`ONE_SIDED_RNA` is set. The message asserts
    ``log_f >= mode`` (the destination holds AT LEAST the RNA the bound accounts for), so only
    ``log_f < mode`` is a contradiction and only that side may be penalised.

    ⚠ Dtype is preserved deliberately: the AMBIG cube is float32 and the 1-D path float64, and the clamp
    must not promote either. ⛔ With the flag off this is exactly ``log_f - mode``, which is what makes
    the default byte-identical rather than approximately so.
    """
    d = log_f - mode
    if one_sided is None:
        if not ONE_SIDED_RNA[0]:
            return d
        return np.minimum(d, d.dtype.type(0.0))
    # per-slot (stage 4d): the mask arrives already broadcast to ``d``'s leading axis; ``where``
    # preserves the dtype, so the AMBIG cube's float32 path stays float32.
    return np.where(one_sided, np.minimum(d, d.dtype.type(0.0)), d)


# f_g ∈ [σ(−10), σ(10)] = [4.5e-5, 1−4.5e-5]. A pure STATE-SPACE bracket: the widest f_g the grid can
# represent, NOT an accuracy knob — but that is a PROPERTY OF A PROPER ψ, not of this constant. It holds
# because both arms are now always written (`_JEFFREYS_REF`): under Beta(½,½) ~0.9% of the reference's mass
# lies outside L=10, and the answer is L-invariant. An improper ψ (either arm omitted) has plateau mass
# growing linearly in L, and then L silently sets the prior strength — which is what the `+0.5·λ` ramp was.
# **L-invariance is the acceptance test for PRIOR-FREE ψ**, where it holds to seven digits.
#
# ⛔⛔ **AND IT IS SCOPED TO PRIOR-FREE ψ BECAUSE THE SHIPPED PIPELINE FAILS IT — MEASURED 2026-08-17. This
# comment claimed the property unconditionally, for the whole file, until then.** Holding the lattice
# spacing `dlam` fixed at 0.3390 and widening ONLY the bracket, Σ|Δ| in object-incidence fragments
# (`pass0_vs_oracle.score_axis`, region+boundary):
#
#     condition                  L=10 (shipped)   L=20 K=119/511    L=40 K=237/1021   resolution-only L=10
#     g00 ss_0.50 capture_off         1,898,257   614,587  0.3238x  592,135  0.3119x  1,926,127  1.0147x
#     g50 ss_0.50 capture_off           116,807   114,974  0.9843x  114,960  0.9842x    117,343  1.0046x
#     g98 ss_0.99 capture_off           119,332   116,162  0.9734x  116,152  0.9733x    119,767  1.0036x
#
# It SATURATES by L=40 and the resolution-only control (K raised at L=10) moves the OTHER way, so the
# effect is the BRACKET and not the lattice. ⭐ The mechanism is the FITTED landscape, which is why
# prior-free ψ is untouched: `landscape.logprior` evaluates at `log rho = log f_c + log M − log E`, and ψ
# can only offer `f_c ∈ [σ(−L), σ(L)]` — at g00 `σ(−10) = 4.540e-05` sits 370x ABOVE the median density
# the fitted prior points at (1.225e-07), so the low end of the bracket is a wall the prior is pushing
# against rather than empty state space. ⛔ Do not read the table as licence to raise L: the derivation
# and the repair are the owner's, and no arm here has priced what a wider bracket costs elsewhere.
#
# NB: production does not read this default — `sweep.solve_chain` threads `logodds_window` (=10.0)
# explicitly, from `CalibrationConfig.sweep_logodds_window`.
_DEFAULT_L = 10.0

# Cache-tiling target for BOTH per-region solves, as a working-set size rather than a row count — the per-row
# footprint differs ~7× between the 1-D grid (K f64) and the 2-D cube (K·K_t f32), so no single row count
# serves both. `_block_rows` turns it into rows.
#
# NOT a model parameter. Every region solves INDEPENDENTLY and every reduction in both solvers is WITHIN a row
# (the ψ logsumexp, the moment sums, the CDF cumsum, and the `post @ log f` gemv), so the block size cannot
# reach the arithmetic — verified bitwise for all five reduction kinds at block sizes from 64 to 65,536.
# It is purely a memory knob, and at genome scale it is the dominant one: the 1-D path runs 357,739
# single-strand regions at K=256, i.e. a **699 MB temporary** per intermediate, ~10 of them live, streamed
# from DRAM. Blocking makes the same arithmetic run out of cache.
_SOLVE_BLOCK_BYTES = 1 << 20


def _block_rows(cells_per_row: int, itemsize: int) -> int:
    return max(1, _SOLVE_BLOCK_BYTES // max(1, int(cells_per_row) * int(itemsize)))


def _lse(a, axis, keepdims=False):
    """Lean numpy log-sum-exp — a drop-in for ``scipy.special._lse(a, axis, keepdims)`` without the
    scipy wrapper overhead (arg validation, ``b``/``return_sign`` handling), which the profiler flagged as
    ~9 s of pure per-call cost across the AMBIG solve. Same max-shift stabilisation; the ``m→0`` guard makes
    an all-``-inf`` slice give ``log(0) = -inf`` and avoids ``(-inf)-(-inf) = nan``."""
    m = np.max(a, axis=axis, keepdims=True)
    m = np.where(np.isfinite(m), m, 0.0)
    with np.errstate(
        divide="ignore"
    ):  # all-(-inf) slice ⇒ log(0) = -inf (correct); suppress the warning
        # Accumulate the sum in float64 even for a float32 cube (a 60-/3600-wide reduction loses too much in
        # f32), then return the input dtype — keeps the f64 path byte-identical, gives the f32 cube an
        # accurate normalizer.
        s = np.sum(np.exp(a - m), axis=axis, keepdims=True, dtype=np.float64)
        r = (m + np.log(s)).astype(m.dtype, copy=False)
    return r if keepdims else np.squeeze(r, axis=axis)


# ──────────────────────────────────────────────────────────────────────────────────────────────────────
# ⭐ THE THREE-COMPONENT STRAND LIKELIHOOD — folded in from `simplex.py`, which was 55 boundaries with ZERO
# public names and exactly ONE importer: this file. A module whose whole surface is one private function
# used in one place is a file a reader has to open to learn nothing, and the flat pile was made of those.
# ⚠ Its TWO-component special case lives in `strand_likelihood.strand_loglik` and is an executable
# REFERENCE, not dead code — `test_strand_likelihood_reference.py` gates that this generalization collapses
# onto it when one RNA strand is dead. That gate did not exist until 2026-08-07 even though `simplex.py`'s
# docstring asserted it did ("the no-regression guard"), which is why the reference read as unused.
# ──────────────────────────────────────────────────────────────────────────────────────────────────────


def _mixture_strand_loglik(
    u_pos, n, f_g, f_pos, f_neg, kappa, od_g, od_r, f_g_ref, f_pos_ref, f_neg_ref
):
    """Three-component gDNA/RNA₊/RNA₋ strand loglik — :func:`strand_loglik` generalized to two RNA strands.

    Broadcasts ``(u_pos, n)`` of shape ``(regions, 1)`` against the lattice ``(f_*)`` of shape
    ``(1, P)`` → ``(regions, P)``. Mean ``N·p`` with ``p = ½·f_g + κ·f₊ + (1−κ)·f₋``.

    **Count-zero-information freeze**: the mean stays LIVE in
    the solved composition ``(f_g, f_pos, f_neg)`` — the legitimate strand channel — but the variance is
    evaluated at the fixed REFERENCE composition ``(f_g_ref, f_pos_ref, f_neg_ref)`` (per-region scalars,
    broadcast). This keeps the heteroscedastic precision (the count still sets a composition-aware variance
    via the reference) while removing the ``f_g``-tilt of the normalizer, so the raw count can no longer
    manufacture a composition preference toward the variance-minimum when the mean degenerates (κ→½). The
    reference is a NEUTRAL structural default at init and the incoming belief in the sweep.
    """
    p = (
        0.5 * f_g + kappa * f_pos + (1.0 - kappa) * f_neg
    )  # LIVE mean channel (the composition solved)
    mean = n * p
    rscale = kappa * (1.0 - kappa)  # κ(1−κ): each RNA strand's μ(1−μ)
    # Variance at the REFERENCE composition (NOT the solved f_g) — the freeze.
    p_ref = 0.5 * f_g_ref + kappa * f_pos_ref + (1.0 - kappa) * f_neg_ref
    var = (
        n * p_ref * (1.0 - p_ref)
        + (n * f_g_ref) ** 2 * 0.25 * od_g
        + (n * f_pos_ref) ** 2 * rscale * od_r
        + (n * f_neg_ref) ** 2 * rscale * od_r
    )
    var = np.maximum(var, _EPS)
    return -0.5 * (u_pos - mean) ** 2 / var - 0.5 * np.log(var)


def _log_fg(lam):
    """``log f_g = log σ(λ)`` — computed via ``scipy.special.log_expit`` (stable: it never forms
    ``1−σ(λ)``, so it stays exact in the depleted tail where the naive ``log(clip(σ(λ), ε))`` underflows)."""
    return log_expit(np.asarray(lam, dtype=np.float64))


def _log1m_fg(lam):
    """``log(1 − f_g) = log σ(−λ) = log f_rna`` — exact via ``log_expit(−λ)``."""
    return log_expit(-np.asarray(lam, dtype=np.float64))


def _logodds_grid(n_grid: int, L: float = _DEFAULT_L):
    """The fixed log-odds lattice: ``λ`` uniform on ``[−L, L]`` (``K = n_grid`` points, ascending) and
    the matching ``f_g = σ(λ)`` (also ascending). Returns ``(lam, fg)``, each length ``K``."""
    lam = np.linspace(-float(L), float(L), int(n_grid))
    return lam, expit(lam)


def _posterior_median_fg(post, lam, fg):
    """Per-region point estimate of ``f_g``: the posterior's ½-QUANTILE, read off the CDF.

    Transform-invariant and robust to the SKEW of the ``f_g`` posterior, which is why ``f_g`` is a median
    and not a mean — measured, the mean is worse at BOTH simplex vertices at every depth, and the vertices
    carry 49–83 % of in-scope error (`vertex_ceiling.py --arm psi_mean`: 1.352 / 1.573 / 3.756 on the three
    in-scope strata, 1.801 on the zero control).

    ⭐⭐⭐ **IT IS A CONTINUOUS QUANTILE, NOT THE GRID POINT WHERE THE CDF FIRST REACHES ½.** The grid mass
    is treated as a histogram — bin edges at the midpoints — and the crossing bin is interpolated. ⛔ The
    predecessor returned ``fg[(cw < 0.5).sum()]``, which SNAPS to a lattice point and therefore carries up
    to half a grid step of quantisation (Δf_g ≈ 0.085 at ``n_grid`` 60, 0.0196 at ``n_grid_ss`` 256). That
    snap was not noise to be averaged away:

    * on an EVIDENCE-FREE object the posterior IS the reference — symmetric, so its ½-quantile is exactly
      ½ — and the snapped version returned the adjacent grid point instead. ⭐ Exact here by the symmetry
      of ``σ``: ``σ(δ) + σ(−δ) = 1``, so the edge between the two central bins is exactly 0.5;
    * at ``κ = ½`` the strand term is bit-flat, so EVERY slot is in that state and the snap was a constant
      ±0.0423 (``n_grid`` 60) that did not shrink with depth — it was the whole closure defect there;
    * and it propagates: `_compose` builds the RNA fractions from ``1 − f_g``, so a snapped ``f_g`` snaps
      the whole composition. `test_relay_mass_rescale`'s ``R_own = 0.5`` is what caught it.

    ⚠ **NOT the "sub-grid MODE" this docstring used to refuse**, and the distinction is the point: a mode
    is an argmax and can chase a single spike, which is what would under-call a skewed or vertex-near
    posterior. A quantile cannot — it is monotone in the CDF and reads the same half-mass crossing the
    snapped version was approximating.

    ⭐⭐⭐ **THE INTERPOLATION IS DONE ON λ, NOT ON f_g, AND THAT IS WHAT KEEPS IT TRANSFORM-INVARIANT.**
    ``σ`` is monotone, so the ½-quantile in ``λ`` maps through it to the ½-quantile in ``f_g`` exactly —
    median equivariance, which is the whole reason ``f_g`` is a median. ⛔ Interpolating in ``f``-space
    instead looks equivalent and is not: the σ grid is highly NON-uniform (spacing ~1e-5 at the ends,
    ~0.085 in the middle), so a bin's midpoint in ``f`` is not the image of its midpoint in ``λ``, and a
    posterior concentrated on one grid point comes back **biased toward ½** — measured **2.71e-03** at
    ``n_grid`` 60 and 1.48e-04 at 256, on every interior point. On ``λ`` the same posterior returns its
    own grid point to **2.2e-16**. ⚠ That case is not synthetic: an unsolved slot's fed-back belief
    produces a one-hot posterior.

    ``post``: (m,K) normalized posterior; ``lam`` the UNIFORM log-odds grid; ``fg`` = σ(λ). Returns (m,)."""
    p = np.asarray(post, np.float64)
    x = np.asarray(lam, np.float64)
    # histogram edges on the UNIFORM λ lattice — the two outer half-bins mirrored.
    edges = np.empty(x.shape[0] + 1, np.float64)
    edges[1:-1] = 0.5 * (x[:-1] + x[1:])
    edges[0] = x[0] - 0.5 * (x[1] - x[0])
    edges[-1] = x[-1] + 0.5 * (x[-1] - x[-2])
    tot = p.sum(axis=1, keepdims=True)
    cdf = np.concatenate(
        [np.zeros((p.shape[0], 1)), np.cumsum(p / np.where(tot > 0.0, tot, 1.0), axis=1)], axis=1
    )
    k = np.clip((cdf < 0.5).sum(axis=1), 1, x.shape[0])
    lo = np.take_along_axis(cdf, (k - 1)[:, None], 1)[:, 0]
    hi = np.take_along_axis(cdf, k[:, None], 1)[:, 0]
    span = hi - lo
    t = np.where(span > 0.0, (0.5 - lo) / np.where(span > 0.0, span, 1.0), 0.5)
    return expit(edges[k - 1] + t * (edges[k] - edges[k - 1]))


def _compose(f_g, w_pos, allow_pos, allow_neg):
    """ψ's composition, as the MAP from its parameters → ``(f_pos, f_neg)``.

    ⭐⭐⭐ **THE COMPOSITION HAS TWO DEGREES OF FREEDOM, NOT THREE, AND THIS IS WHERE THAT IS WRITTEN
    DOWN.** ψ solves a point on the 2-simplex, parametrised by

        λ  — the gDNA-vs-RNA LEVEL, read out as ``f_g`` (the posterior median over the λ grid)
        θ  — the RNA-internal TILT, a SHARE with no absolute scale, read out as ``w_pos``

    and the composition is their image::

        f_pos = (1 − f_g)·w_pos        f_neg = (1 − f_g)·(1 − w_pos)

    ⛔⛔ **CLOSURE IS THEREFORE STRUCTURAL — IT CANNOT FAIL**, because the map lands on the simplex.
    ``f_g + f_pos + f_neg = 1`` identically (measured: ``|SUM − 1| ≤ 1.11e-16`` on both ψ paths at every
    κ and depth), so no consumer has to check it and no arithmetic has to be renormalised.

    ⚠ **WHAT THIS REPLACED, AND WHY IT WAS WRONG.** ψ used to read out all THREE coordinates
    independently — ``f_g`` as the posterior MEDIAN over λ, and ``f_pos``/``f_neg`` as posterior MEANS of
    the grid quantity ``1 − f_g``. Mixing a quantile with expectations does not land on the simplex:

        SUM = median(f_g) + (1 − mean(f_g)) = 1 + median(f_g) − mean(f_g)

    i.e. **the closure error was exactly the posterior's SKEW** (verified to 5.8e-15), and only 74–77 % of
    real objects closed. ⛔ It was never a decision: the median for ``f_g`` was argued for, and the RNA
    fractions were never chosen at all — they fell out as expectations of a grid array.

    ⛔ **The repair is NOT "take means everywhere".** That closes too, by linearity of expectation, and was
    measured on all 16 ladder conditions at **1.352 / 1.573 / 3.756** on the three in-scope strata and
    **1.801** on the zero control: the median is closer to the truth at both simplex vertices, where
    49–83 % of in-scope error lives. ⛔ Nor is it renormalising three numbers at publication, which would
    make a 15 %-short object indistinguishable from a solved one. Nothing here is rescaled: ``f_g`` and the
    RNA total are exact complements *by parametrisation*, and the tilt is estimated as a SHARE because a
    share is what it is.

    ``w_pos`` is the + strand's share of RNA — the RNA-mass-weighted posterior share on an AMBIG slot.
    ⛔ **ADMISSIBILITY IS ENFORCED HERE, NOT BY THE CALLER**, because a share that ignores it loses RNA
    silently: with only the + strand free, ``w_pos = ½`` would place half the RNA on a forbidden strand,
    where it is zeroed and simply vanishes. The tilt of a single-strand slot is structurally locked, so
    the admissible strand takes the whole RNA total whatever ``w_pos`` says. ⚠ A slot with NEITHER strand
    admissible has no RNA to place and returns ``(0, 0)`` — its composition is ``f_g`` alone, which is the
    honest statement and not a closure failure (nothing dispatches such a slot to a ψ solve anyway).
    """
    fr = 1.0 - np.asarray(f_g, np.float64)
    ap = np.asarray(allow_pos, bool)
    an = np.asarray(allow_neg, bool)
    # ⛔ ``w_pos`` is a SHARE, so it is clamped to [0,1] here rather than trusted. Unclamped, a share
    #   outside the range produces a NEGATIVE fraction that still sums to 1 — a composition that passes a
    #   closure check and is nonsense. ⚠ Not reachable from the shipped AMBIG caller (it is a ratio of two
    #   non-negative expectations), which is exactly why the constraint belongs where it can be seen.
    w = np.clip(np.asarray(w_pos, np.float64), 0.0, 1.0)
    w = np.where(ap & an, w, np.where(ap, 1.0, 0.0))
    return np.where(ap, fr * w, 0.0), np.where(an, fr * (1.0 - w), 0.0)


def _slice_rows(a, msk):
    """One prior array's rows for a block, or ``None``. ⚠ Kept module-level so
    :class:`CompositionPriors` can use it; the solver's own ``_s`` closure is the same operation."""
    return None if a is None else np.asarray(a)[msk]


@dataclass(frozen=True, slots=True)
class CompositionPriors:
    """ψ's two fitted composition arms, carried as ONE object.

    ⭐⭐ **Why a pair rather than two parameters.** ψ's solvers already take eighteen arguments, and the
    two arms are one concept: the fitted population density for each component of the gDNA-vs-RNA split.
    Replacing ``global_logprior`` with this adds no parameter and gives the concept a name.

    ⛔⛔ **AND IT MAKES A WHOLE CLASS OF BUG UNREPRESENTABLE.** Each arm has to be row-SLICED per solve
    block and, on the fine single-strand grid, RE-GRIDDED from the coarse λ lattice. Those are two
    separate operations applied at two separate call sites, so a second arm threaded by hand could
    easily be sliced and not regridded — which would not raise, would not change shape, and would
    silently evaluate one component's prior on the wrong lattice. :meth:`select` and :meth:`regrid`
    apply to BOTH members at once, so the two arms cannot drift apart.

    ``None`` on either member means that arm is **not fitted** and takes its derived reference — which
    is not the same as "no term" (see :func:`_gdna_arm`).

    ⛔ The reference LOCATION (a per-slot mean ``m_i`` moving the reference off ½ — the structural
    0.75 assertion and the measured intron override) was a third member here until 2026-08-24, when
    the owner refuted the concept: a location is a prior assertion, pass-0's job is to LEARN the
    prior, and where the strand channel is dead the term was the entire answer at any depth. The
    reference is now the symmetric Jeffreys measure and nothing else; background information enters
    as LIKELIHOOD terms (the density λ-factor), whose precision scales with counts.
    """

    gdna: np.ndarray | None = None
    rna: np.ndarray | None = None

    def select(self, msk) -> "CompositionPriors":
        """Both arms restricted to one solve block's rows."""
        return CompositionPriors(_slice_rows(self.gdna, msk), _slice_rows(self.rna, msk))

    def regrid(self, n_from: int, n_to: int, L: float) -> "CompositionPriors":
        """Both arms carried from the coarse λ lattice to the fine single-strand one."""
        return CompositionPriors(
            _regrid_global(self.gdna, n_from, n_to, L),
            _regrid_global(self.rna, n_from, n_to, L),
        )


#: The prior-free solve — both arms take their derived reference. ⭐ A first-class configuration, not a
#: degenerate one: pass-0 runs here by design.
_NO_PRIORS = CompositionPriors()


def _regrid_global(glp, n_from, n_to, L):
    """Interpolate a ``(m, n_from)`` global-logprior (evaluated on the σ(λ) grid at ``n_from``) onto the
    ``n_to`` σ(λ) grid (Fix 3). The single-strand solve runs on a finer grid than the AMBIG cube, so the
    coarse-grid global prior is regridded to feed it; the global is smooth in ``f_g``, so linear
    interpolation is exact to interpolation accuracy. ``None`` / equal grids ⇒ passthrough (bit-identical)."""
    if glp is None or int(n_from) == int(n_to):
        return glp
    _, fc = _logodds_grid(int(n_from), L)
    _, ff = _logodds_grid(int(n_to), L)
    j = np.clip(np.searchsorted(fc, ff), 1, int(n_from) - 1)
    x0, x1 = fc[j - 1], fc[j]
    t = np.clip((ff - x0) / np.maximum(x1 - x0, _EPS), 0.0, 1.0)  # (n_to,)
    g = np.asarray(glp, np.float64)
    return g[:, j - 1] + t[None, :] * (g[:, j] - g[:, j - 1])


def _gdna_arm(lam, global_logprior):
    """The gDNA arm of ψ over the λ grid → broadcastable to ``(m, K)``.

    ``_JEFFREYS_REF`` reference ``+½·log f_g`` → ``(1, K)``. When a prior IS fitted, the fitted
    ``logP_g(log ρ_g)`` (pre-evaluated on THIS ``f_g`` grid → ``(m, K)``) is **ADDED** to it.

    ⚠ It used to REPLACE the reference, on the argument that the reference "carries no information and must
    not be double-counted". That argument is wrong twice over. The reference is not an information claim to
    be superseded — it is the **measure** ψ is written against, and it is the ONLY term bounding this arm at
    ``f_g → 0`` (`_rna_arm` bounds the other vertex and is never replaced, so the two arms were not even
    treated alike). Deleting it left ψ improper at the vertex the fitted prior most often points at, which is
    the documented region-1055 crush. Bayes composes a prior with a measure by ADDITION in log space; there is
    no double-count to avoid.

    ``None`` means "not fitted", **not** "no term"."""
    ref = _JEFFREYS_REF * _log_fg(lam)[None, :]
    if global_logprior is None:
        return ref
    return ref + np.asarray(global_logprior, np.float64)


def _rna_arm(lam, rna_logprior=None):
    """The RNA-**total** arm of ψ over the λ grid → broadcastable to ``(m, K)``.

    The ``_JEFFREYS_REF`` reference ``+½·log(1 − f_g)`` → ``(1, K)``, plus a fitted ``logP_r`` when one is
    supplied — the exact mirror of :func:`_gdna_arm`. ``None`` means "not fitted", **not** "no term".

    This is the **two-group** arm (gDNA vs RNA-total): the per-strand split is the nuisance tilt, integrated
    out on the θ axis, and needs no prior of its own.

    ⛔⛔ **NOTHING FITS ``logP_r`` YET, AND THE COST OF THAT IS MEASURED.** Until one does, the reference
    alone bounds the ``f_g → 1`` vertex — and unlike its gDNA twin it is never swamped by evidence, so it
    is a FIXED repulsion of **3.107 nats** at ``f_g = 0.999`` relative to ``f_g = ½`` (a 22:1 handicap).
    Measured 2026-08-15 on the 16-condition ladder: objects whose TRUE ``f_g ≥ 0.999`` carry **49–83 %** of
    all calibration error on the three in-scope strata, read **0.13–0.23 below** the vertex. ⭐ The
    parameter exists now so that an estimator can close that asymmetry; the socket is not speculative
    surface, it is the repair's landing point."""
    ref = _JEFFREYS_REF * _log1m_fg(lam)[None, :]
    if rna_logprior is None:
        return ref
    return ref + np.asarray(rna_logprior, np.float64)


def _tilt_grid(n_tilt: int) -> np.ndarray:
    """The RNA-internal tilt grid as the ANGLE ``θ ∈ [−π/2, π/2]`` (``K_t`` points), with ``τ = sin θ``.

    Gridding θ (not τ) is what makes the Berger–Bernardo tilt conditional ``(1−τ²)^{−½}`` vanish identically:
    ``|dτ/dθ| = cos θ = (1−τ²)^{½}`` cancels it exactly, so **no tilt term is written at all** and the ψ
    reference is the same expression for AMBIG as for single-strand regions (module docstring §3). It also
    removes the endpoint singularity outright — no clipping, no Gauss–Jacobi weights, no constant.

    Resolution follows the reference measure rather than being uniform in τ: at ``K_t=60`` the τ-spacing is
    ~0.053 near balanced tilt and ~0.0014 near strand-purity (vs a flat 0.034). That is the intended trade —
    it spends grid on the strand-pure boundaries, where distinguishing a pure strand from a small antisense leak is
    the high-stakes call, and economizes on the balanced middle, where the distinction rarely matters.

    ``θ = ±π/2`` ⇒ ``τ = ±1`` ⇒ all RNA on one strand; ``θ = 0`` ⇒ balanced. Only AMBIG regions integrate it."""
    return np.linspace(-0.5 * np.pi, 0.5 * np.pi, int(n_tilt))


def _single_strand_mask(allow_pos, allow_neg) -> np.ndarray:
    """The regions the 1-D (Phase-1) solver is valid for: exactly one strand live (tilt determined)."""
    ap = np.asarray(allow_pos, bool)
    an = np.asarray(allow_neg, bool)
    return ap ^ an


def _ambig_mask(allow_pos, allow_neg) -> np.ndarray:
    """AMBIG regions (both strands live) — the Phase-2 2-D ``(λ, τ)`` path."""
    return np.asarray(allow_pos, bool) & np.asarray(allow_neg, bool)


def _local_loglik_logodds(
    u_pos,
    u_neg,
    allow_pos,
    allow_neg,
    kappa,
    od_g,
    od_r,
    lam,
    fg,
    f_g_ref,
    f_pos_ref,
    f_neg_ref,
    priors: "CompositionPriors | None" = None,
    gdna_imp_mode=None,
    gdna_imp_prec=None,
    rna_imp_mode=None,
    rna_imp_prec=None,
    rna_one_sided=None,
    lam_logprior=None,
    lam_imp_mode=None,
    lam_imp_prec=None,
):
    """ψ over the log-odds grid for single-strand regions (strand mixture, the two arms, imputation), evaluated
    at ``f_g = σ(λ)`` with the live strand carrying ``f_active = 1 − f_g``. Returns ``(m, K)``.

    ψ = strand + ``_gdna_arm`` + ``_rna_arm`` + messages. **No Jacobian** — on the two-group axis the log-rate
    conversions cancel ``log σ'(λ)`` exactly (module docstring §2). Both arms are ALWAYS written: a fitted
    ``logP`` if we have one, else the ``_JEFFREYS_REF`` reference. Omitting one is not neutral.

    ``priors`` carries each arm already evaluated on THIS ``fg`` grid → ``(m, K)``; a ``None`` member ⇒ that arm
    takes its reference (a PRIOR-FREE solve is not a REFERENCE-FREE solve).

    ``f_g_ref`` / ``f_pos_ref`` / ``f_neg_ref`` (per-region ``(m,)``) are the count-zero-information freeze
    reference: the strand mixture's variance is evaluated at THIS fixed
    composition — not the grid ``f_g`` being integrated — so the count sets precision, not composition."""
    ap = np.asarray(allow_pos, bool)
    an = np.asarray(allow_neg, bool)
    pos_live = (ap & ~an)[:, None]  # (m,1)
    neg_live = (an & ~ap)[:, None]
    fg2 = fg[None, :]  # (1,K)
    f_act = 1.0 - fg2
    f_pos = np.where(pos_live, f_act, 0.0)  # (m,K)
    f_neg = np.where(neg_live, f_act, 0.0)
    n = (u_pos + u_neg)[:, None]
    # ── strand mixture — variance FROZEN at the reference composition (count-zero-info §2). Broadcasts
    #    the (m,1) reference against the (m,K) grid. ──
    psi = _mixture_strand_loglik(
        u_pos[:, None],
        n,
        fg2,
        f_pos,
        f_neg,
        kappa,
        od_g,
        od_r,
        np.asarray(f_g_ref, np.float64)[:, None],
        np.asarray(f_pos_ref, np.float64)[:, None],
        np.asarray(f_neg_ref, np.float64)[:, None],
    )
    # ── the two composition arms: gDNA and RNA-total. Each is its fitted logP if we have one, else the
    #    derived Jeffreys reference. ALWAYS both — see `_gdna_arm` / `_rna_arm`. Together they make ψ proper
    #    (Beta(½,½) when neither is fitted); the gDNA arm alone would leave f_g→1 unbounded, and the RNA arm
    #    alone would leave f_g→0 unbounded. ──
    _p = priors or _NO_PRIORS
    # ⛔ No location term: the reference is the SYMMETRIC Jeffreys measure and asserts nothing
    #    (owner refutation, 2026-08-24 — see CompositionPriors). Background information enters below
    #    as the intron-factory λ-factor, a LIKELIHOOD whose precision scales with counts.
    psi = psi + _gdna_arm(lam, _p.gdna) + _rna_arm(lam, _p.rna)
    # ── the gDNA INTRON FACTORY λ-factor: a per-region (m,K) log-likelihood on
    #    the λ axis, ``log NegBinom(f_g·C; ρ_bg·E_g, α_eff)``, ADDED as its own term rather than folded into
    #    the gDNA arm — it is a per-region LIKELIHOOD, not a population density, so it does not belong inside
    #    a term whose units are `log P(log ρ)`. It deconvolves confident gDNA from
    #    introns against the intergenic background; zero on non-intron regions ⇒ a no-op there. ──
    if lam_logprior is not None:
        psi = psi + np.asarray(lam_logprior, np.float64)
    # ── imputation messages: LOG-FRACTION Gaussians (the overhaul). The mode is a log-FRACTION target
    #    (``log`` of the imputed fraction, built in ``_scan``); evaluated against ``log f_c(λ)``. No clip —
    #    an off-grid target (source denser than the dst can hold) is a bounded monotone pull toward the
    #    boundary, governed by precision (D-plan P6, verify-don't-clip). ──
    log_fg = _log_fg(lam)[None, :]  # log f_g = log σ(λ)
    log_fact = _log1m_fg(lam)[None, :]  # log(1−f_g) = log f_active (the single live strand)
    if gdna_imp_mode is not None and gdna_imp_prec is not None:
        m_ = np.asarray(gdna_imp_mode, np.float64)[:, None]
        p_ = np.asarray(gdna_imp_prec, np.float64)[:, None]
        psi = psi - 0.5 * p_ * (log_fg - m_) ** 2
    if rna_imp_mode is not None and rna_imp_prec is not None:
        # single-strand: the live strand carries f_active = 1−f_g; the per-strand precision gates which
        # message applies (the dead strand's prec is 0 → no-op). Both evaluate against log f_active.
        _os = None if rna_one_sided is None else np.asarray(rna_one_sided, bool)[:, None]
        for ms, ps in ((rna_imp_mode[0], rna_imp_prec[0]), (rna_imp_mode[1], rna_imp_prec[1])):
            psi = (
                psi
                - 0.5
                * np.asarray(ps, np.float64)[:, None]
                * _rna_residual(log_fact, np.asarray(ms, np.float64)[:, None], _os) ** 2
            )
    # ── the SINGLE-λ composition message (the the-single-lambda-combine rank-1 fix): ONE Gaussian on the log-odds grid variable λ
    #    DIRECTLY (not on log f_c) — the one gDNA-vs-RNA-total DOF, so ψ counts it ONCE, not twice
    # Enrichment-invariant: λ carries no reframe. ──
    if lam_imp_mode is not None and lam_imp_prec is not None:
        lm_ = np.asarray(lam_imp_mode, np.float64)[:, None]
        lp_ = np.asarray(lam_imp_prec, np.float64)[:, None]
        psi = psi - 0.5 * lp_ * (lam[None, :] - lm_) ** 2
    # ── NO change-of-variable Jacobian, and that is deliberate: a fitted `logP` is a density in LOG-rate,
    #    so its conversion to a linear-rate density (−log f_c, up to a constant) cancels log σ'(λ) exactly,
    #    ONCE PER COMPONENT — which is why the cancellation keeps holding as each arm acquires a fitted
    #    prior. ⛔ The reference IS written, always, in both arms; this comment read "and NO reference
    #    prior … ψ_λ = strand + logP_g + logP_r, bare" until 2026-08-15, describing a retired design in
    #    which a fitted prior REPLACED the reference. ⇒ ψ_λ = strand + (ref + logP_g) + (ref + logP_r). ──
    # ⭐ ``f_pos``/``f_neg`` are LOCAL to the strand mixture and are no longer returned: the caller builds
    #   the composition from the PARAMETERS (:func:`_compose`), so it never needed the grid arrays.
    return psi


def _solve_regions_logodds(
    u_pos,
    u_neg,
    allow_pos,
    allow_neg,
    f_g_ref,
    f_pos_ref,
    f_neg_ref,
    *,
    kappa,
    od_g,
    od_r,
    n_grid,
    L: float = _DEFAULT_L,
    priors: "CompositionPriors | None" = None,
    gdna_imp_mode=None,
    gdna_imp_prec=None,
    rna_imp_mode=None,
    rna_imp_prec=None,
    rna_one_sided=None,
    lam_logprior=None,
    lam_imp_mode=None,
    lam_imp_prec=None,
) -> RegionDeconv:
    """The log-odds 1-D per-region solve for SINGLE-STRAND regions.

    Read-out: ``f_g`` = posterior median over the ``λ`` grid, and ``f_pos``/``f_neg`` are its IMAGE under
    :func:`_compose` — a single-strand slot's tilt is structurally locked, so the RNA total ``1 − f_g``
    goes entirely to the admissible strand and the composition closes by construction.
    ``*_frac_var`` = the grid-moment ``Var(log f_c)``. The dead strand is locked-certain (var 0); zero-mass
    regions report 0. ``f_g_ref``/``f_pos_ref``/``f_neg_ref`` (per-region) are the count-zero-info variance
    freeze reference (§2). AMBIG regions are out of contract — masked out."""
    lam, fg = _logodds_grid(int(n_grid), L)
    psi = _local_loglik_logodds(
        u_pos,
        u_neg,
        allow_pos,
        allow_neg,
        kappa,
        od_g,
        od_r,
        lam,
        fg,
        f_g_ref,
        f_pos_ref,
        f_neg_ref,
        priors=priors,
        gdna_imp_mode=gdna_imp_mode,
        gdna_imp_prec=gdna_imp_prec,
        rna_imp_mode=rna_imp_mode,
        rna_imp_prec=rna_imp_prec,
        rna_one_sided=rna_one_sided,
        lam_logprior=lam_logprior,
        lam_imp_mode=lam_imp_mode,
        lam_imp_prec=lam_imp_prec,
    )
    ap = np.asarray(allow_pos, bool)
    an = np.asarray(allow_neg, bool)
    post = np.exp(psi - _lse(psi, axis=1, keepdims=True))  # (m,K)
    # f_g posterior median (fg ascending ⇒ cumulative CDF directly)
    f_g = _posterior_median_fg(post, lam, fg)
    # ⭐⭐ THE COMPOSITION IS THE MAP FROM THE PARAMETERS — see :func:`_compose`. A single-strand slot has
    #    ONE degree of freedom (λ); its tilt is structurally locked, so the RNA total goes entirely to the
    #    admissible strand and the simplex closes by construction.
    # precision state = Var(log f_c), moment-matched on the grid. The dead strand is locked-certain
    # (f=0) → var 0. Capping is AUTOMATIC: the send prec_log = 1/(var+σ²+pois) ≤ 1/(σ²+pois).
    Lg = _log_fg(lam)
    mLg = post @ Lg
    var_g = np.maximum(post @ (Lg * Lg) - mLg * mLg, 0.0)
    La = _log1m_fg(lam)
    mLa = post @ La
    var_act = np.maximum(post @ (La * La) - mLa * mLa, 0.0)
    var_pos = np.where(ap & ~an, var_act, 0.0)
    var_neg = np.where(an & ~ap, var_act, 0.0)
    # ⛔ ``_compose`` runs on the CLIPPED ``f_g``, so the RNA total is its exact complement and the clip
    #   on ``f_pos``/``f_neg`` — which is what used to break the simplex — has nothing left to do.
    active = (u_pos + u_neg) > 0.0
    f_g = np.where(active, np.clip(f_g, 0.0, 1.0), 0.0)
    f_pos, f_neg = _compose(f_g, 0.0, ap, an)
    f_pos = np.where(active, f_pos, 0.0)
    f_neg = np.where(active, f_neg, 0.0)
    var_g = np.where(active, var_g, 0.0)
    var_pos = np.where(active, var_pos, 0.0)
    var_neg = np.where(active, var_neg, 0.0)
    return RegionDeconv(
        gdna_frac=f_g,
        rna_pos_frac=f_pos,
        rna_neg_frac=f_neg,
        gdna_frac_var=var_g,
        rna_pos_frac_var=var_pos,
        rna_neg_frac_var=var_neg,
    )


def _solve_ambig_logodds(
    u_pos,
    u_neg,
    f_g_ref,
    f_pos_ref,
    f_neg_ref,
    *,
    kappa,
    od_g,
    od_r,
    n_grid,
    L: float = _DEFAULT_L,
    n_tilt: int | None = None,
    priors: "CompositionPriors | None" = None,
    gdna_imp_mode=None,
    gdna_imp_prec=None,
    rna_imp_mode=None,
    rna_imp_prec=None,
    rna_one_sided=None,
    lam_logprior=None,
    lam_imp_mode=None,
    lam_imp_prec=None,
    theta_imp_mode=None,
    theta_imp_prec=None,
) -> RegionDeconv:
    """The 2-D ``(λ, θ)`` solve for AMBIG regions (both strands live). Grids the gDNA-vs-RNA-total log-odds
    ``λ`` (outer, ``K = n_grid``) and the tilt ANGLE ``θ = arcsin(τ)`` (inner, ``K_t = n_tilt`` or
    ``n_grid``), evaluates ψ on the ``(m, K, K_t)`` cube, and **marginalizes θ** (``logsumexp``) for the
    ``f_g`` read-out. ``f_g`` = posterior median over the θ-marginal λ-posterior; ``f_pos``/``f_neg`` = means
    over the full 2-D posterior.

    **No Jacobian and no tilt term are written** — and that is the point of the θ coordinate. The
    Berger–Bernardo tilt conditional ``(1−τ²)^{−½}`` is cancelled identically by ``|dτ/dθ| = cos θ``
    (``_tilt_grid``), and on the two-group λ axis the log-rate conversions cancel ``log σ'(λ)``. So ψ here is
    the SAME expression as the 1-D path: strand + ``_gdna_arm`` + ``_rna_arm`` + messages. The 1-DOF/AMBIG
    reference asymmetry is closed **identically**, not approximately.

    ``priors``' members are ``(m, K)`` on the σ(λ) grid (broadcast over θ). The cube is only
    materialized for the AMBIG subset (the caller masks); ``K·K_t`` is the per-region cost."""
    lam, fg = _logodds_grid(int(n_grid), L)  # (K,)
    Kt = int(n_tilt) if n_tilt else int(n_grid)
    theta = _tilt_grid(Kt)  # (Kt,) the ANGLE; τ = sin θ
    tau = np.sin(theta)  # (Kt,) the tilt itself — smooth and exact across the whole θ domain
    f_act = (1.0 - fg)[:, None]  # (K,1)
    f_pos_kt = f_act * (1.0 + tau[None, :]) / 2.0  # (K,Kt) f64 (reused f64 in the moment sums)
    f_neg_kt = f_act * (1.0 - tau[None, :]) / 2.0  # (K,Kt) f64
    n = u_pos + u_neg
    # ── float32 cube (authorized small-tolerance perf): the (m,K,Kt) log-posterior and its per-cell inputs
    #    are f32 — halves the cube's memory, ~2x the elementwise exp/log. Every REDUCTION below accumulates
    #    in f64 and the (m,K) marginals are lifted to f64, so medians / means / variances keep precision. ──
    F = np.float32
    fg32 = fg.astype(F)
    fpk = f_pos_kt.astype(F)  # (K,Kt) f32 cube grid
    fnk = f_neg_kt.astype(F)
    # ── strand mixture over the cube: (m,1,1)×(1,K,1)×(1,K,Kt) → (m,K,Kt) f32 ──
    psi = _mixture_strand_loglik(
        np.asarray(u_pos, F)[:, None, None],
        np.asarray(n, F)[:, None, None],
        fg32[None, :, None],
        fpk[None, :, :],
        fnk[None, :, :],
        kappa,
        od_g,
        od_r,
        np.asarray(f_g_ref, F)[:, None, None],
        np.asarray(f_pos_ref, F)[:, None, None],
        np.asarray(f_neg_ref, F)[:, None, None],
    )
    # ── LOG-fraction grids (the overhaul): log f_g (τ-independent) + log f_pos/f_neg over the cube,
    #    floored at one pseudo-fragment 1/(n+1) (TRAPS: no-prior-means-haldane: the τ=±1 boundaries have f_s=0 → log(0); the count floor
    #    keeps it finite + consistent with pois_log). ──
    log_fg_grid = _log_fg(lam)  # (K,) f64 = log f_g (moments use f64)
    log_fg32 = log_fg_grid.astype(F)  # (K,) f32 for the cube message
    frac_floor = (1.0 / (n + 1.0)).astype(F)[:, None, None]  # (m,1,1) f32
    # ``log ∘ max ≡ max ∘ log``, BITWISE, because numpy's float32 ``log`` is monotone — verified
    # exhaustively over ALL 1,065,353,217 float32 values in [0,1], which is the entire domain both
    # arguments live in (a fraction and ``1/(n+1)``). So the log is taken on the (K,K_t) GRID and the
    # (m,K,K_t) cube only ever sees a ``maximum``: ~140× fewer transcendental evaluations, same bits.
    # ``τ = ±1`` puts one strand's fraction at exactly 0 ⇒ ``log 0 = −inf``, which the floor then discards.
    with np.errstate(divide="ignore"):
        log_fpk, log_fnk, log_floor = np.log(fpk), np.log(fnk), np.log(frac_floor)
    log_fpos = np.maximum(log_fpk[None, :, :], log_floor)  # (m,K,Kt) f32
    log_fneg = np.maximum(log_fnk[None, :, :], log_floor)
    # ── the two composition arms (θ-independent — they live on the λ axis, which is exactly what makes the
    #    tilt a nuisance). Identical call to the 1-D path; see `_gdna_arm` / `_rna_arm`. ──
    _p = priors or _NO_PRIORS
    # ⛔ No location term (owner refutation, 2026-08-24 — see CompositionPriors): the reference is
    #   the symmetric Jeffreys measure alone, identical to the 1-D path.
    psi += np.asarray(_gdna_arm(lam, _p.gdna) + _rna_arm(lam, _p.rna), F)[:, :, None]
    # ── the gDNA INTRON FACTORY λ-factor (θ-independent — it lives on the λ axis), ADDED like the arms; the
    #    [:, :, None] broadcast makes it constant across the tilt, so θ is integrated out cleanly. ──
    if lam_logprior is not None:
        psi += np.asarray(lam_logprior, F)[:, :, None]
    # ── ⭐ the FRAGMENT-LENGTH λ-factor. θ-independent — the length channels do not depend on the strand
    #    tilt at all — so it broadcasts across the cube and θ integrates out cleanly. ⭐ **That
    #    independence is precisely why this source speaks on an AMBIG region where the strand term cannot**:
    #    the Schur complement that zeroes a rank-1-in-θ term does not apply to a term with no θ. ──
    # ── gDNA LOG-fraction message on log f_g (τ-independent) ──
    if gdna_imp_mode is not None and gdna_imp_prec is not None:
        mo = np.asarray(gdna_imp_mode, F)[:, None, None]
        pr = np.asarray(gdna_imp_prec, F)[:, None, None]
        psi -= F(0.5) * pr * (log_fg32[None, :, None] - mo) ** 2
    # ── per-strand RNA LOG-fraction messages on log f_pos/log f_neg (τ-dependent — inside the cube) ──
    if rna_imp_mode is not None and rna_imp_prec is not None:
        for log_f, ms, ps in (
            (log_fpos, rna_imp_mode[0], rna_imp_prec[0]),
            (log_fneg, rna_imp_mode[1], rna_imp_prec[1]),
        ):
            _os = None if rna_one_sided is None else np.asarray(rna_one_sided, bool)[:, None, None]
            psi -= (
                F(0.5)
                * np.asarray(ps, F)[:, None, None]
                * _rna_residual(log_f, np.asarray(ms, F)[:, None, None], _os) ** 2
            )
    # ── the SINGLE-λ composition message on λ DIRECTLY (θ-INDEPENDENT — it lives on the λ axis, which is
    #    exactly what makes the tilt a nuisance): one Gaussian, ψ counts the g-vs-R DOF ONCE. ──
    if lam_imp_mode is not None and lam_imp_prec is not None:
        lm_ = np.asarray(lam_imp_mode, F)[:, None, None]
        lp_ = np.asarray(lam_imp_prec, F)[:, None, None]
        psi -= F(0.5) * lp_ * (lam.astype(F)[None, :, None] - lm_) ** 2
    # ── the TILT message on θ (λ-INDEPENDENT — the separate strand-tilt DOF an AMBIG region needs; not part of
    #    the g-vs-R double-count): a Gaussian on the θ = arcsin(τ) grid. ──
    if theta_imp_mode is not None and theta_imp_prec is not None:
        tm_ = np.asarray(theta_imp_mode, F)[:, None, None]
        tp_ = np.asarray(theta_imp_prec, F)[:, None, None]
        psi -= F(0.5) * tp_ * (theta.astype(F)[None, None, :] - tm_) ** 2
    # θ-marginal λ-posterior (m,K) — lift to f64 so the posterior median + moments are full-precision.
    psi_lam = _lse(psi, axis=2).astype(np.float64)
    post_lam = np.exp(psi_lam - _lse(psi_lam, axis=1, keepdims=True))
    f_g = _posterior_median_fg(post_lam, lam, fg)
    # precision state = Var(log f_g) over the θ-marginal λ-posterior (TRAPS: two-gaussians-one-latent).
    mLg = post_lam @ log_fg_grid
    var_g = np.maximum(post_lam @ (log_fg_grid * log_fg_grid) - mLg * mLg, 0.0)
    # Var(log f_pos/neg) over the FULL 2-D posterior (f32 cube; sums accumulate in f64), and the TILT
    # SHARE that :func:`_compose` maps into the composition.
    flat = psi.reshape(psi.shape[0], -1)
    post2d = np.exp(flat - _lse(flat, axis=1, keepdims=True)).reshape(psi.shape)  # (m,K,Kt) f32
    fp_grid = fpk[None, :, :]
    fn_grid = fnk[None, :, :]
    # ⭐ ``f_pos_kt + f_neg_kt = 1 − f_g`` POINTWISE, so these two sum to the posterior MEAN of the RNA
    #   total and their RATIO is the RNA-mass-weighted posterior mean of the + share — which is the
    #   natural estimator of θ, since the tilt matters in proportion to how much RNA is there.
    m_pos = np.sum(post2d * fp_grid, axis=(1, 2), dtype=np.float64)
    m_neg = np.sum(post2d * fn_grid, axis=(1, 2), dtype=np.float64)
    rna = m_pos + m_neg
    w_pos = np.where(rna > 0.0, m_pos / np.where(rna > 0.0, rna, 1.0), 0.5)
    mLp = np.sum(post2d * log_fpos, axis=(1, 2), dtype=np.float64)
    mLn = np.sum(post2d * log_fneg, axis=(1, 2), dtype=np.float64)
    var_pos = np.maximum(
        np.sum(post2d * log_fpos * log_fpos, axis=(1, 2), dtype=np.float64) - mLp * mLp, 0.0
    )
    var_neg = np.maximum(
        np.sum(post2d * log_fneg * log_fneg, axis=(1, 2), dtype=np.float64) - mLn * mLn, 0.0
    )
    # ⭐ AMBIG: BOTH strands are admissible, so the tilt share splits the RNA total between them.
    #   The clip is on ``f_g`` ALONE — ``f_pos``/``f_neg`` are its image under :func:`_compose` and are
    #   in [0,1] by construction, so clipping them independently (which is what used to break the
    #   simplex) has nothing left to do.
    active = n > 0.0
    f_g = np.where(active, np.clip(f_g, 0.0, 1.0), 0.0)
    f_pos, f_neg = _compose(f_g, w_pos, np.ones_like(active), np.ones_like(active))
    f_pos = np.where(active, f_pos, 0.0)
    f_neg = np.where(active, f_neg, 0.0)
    var_g = np.where(active, var_g, 0.0)
    var_pos = np.where(active, var_pos, 0.0)
    var_neg = np.where(active, var_neg, 0.0)
    return RegionDeconv(
        gdna_frac=f_g,
        rna_pos_frac=f_pos,
        rna_neg_frac=f_neg,
        gdna_frac_var=var_g,
        rna_pos_frac_var=var_pos,
        rna_neg_frac_var=var_neg,
    )


def _solve_regions_logodds_all(
    u_pos,
    u_neg,
    allow_pos,
    allow_neg,
    mass_unspl,
    mass_spliced,
    *,
    kappa,
    od_g,
    od_r,
    n_grid,
    L: float = _DEFAULT_L,
    n_tilt: int | None = None,
    n_grid_ss: int | None = None,
    priors: "CompositionPriors | None" = None,
    gdna_imp_mode=None,
    gdna_imp_prec=None,
    rna_imp_mode=None,
    rna_imp_prec=None,
    rna_one_sided=None,
    lam_logprior=None,
    lam_imp_mode=None,
    lam_imp_prec=None,
    theta_imp_mode=None,
    theta_imp_prec=None,
    fg_ref=None,
    fpos_ref=None,
    fneg_ref=None,
) -> RegionDeconv:
    """The full per-region log-odds dispatcher (Phase 3 #1): routes single-strand regions to the 1-D
    ``λ`` solve (:func:`_solve_regions_logodds`) and AMBIG regions to the 2-D ``(λ, τ)`` solve
    (:func:`_solve_ambig_logodds`), scattering both into full-length arrays. G1 / zero-mass regions
    report 0 (``sweep.solve_chain`` keeps their signature-binary init via the ``solvable`` write-back). A
    drop-in for the lattice ``_local_loglik``+``_region_marginals`` pair: same ψ terms — the log-density
    log-fraction Gaussian messages + the global prior — evaluated on the ``σ(λ)`` log-odds grid.

    All array inputs are full length ``m``; ``priors``' members are ``(m, K)`` on the σ(λ) grid;
    ``gdna_imp_*`` are ``(m,)``; ``rna_imp_*`` are 2-tuples of ``(m,)``. Each is sub-indexed per class."""
    m = int(np.asarray(u_pos).shape[0])
    ap_all = np.asarray(allow_pos, bool)
    an_all = np.asarray(allow_neg, bool)
    # Count-zero-info variance-freeze reference (§2, TRAPS: measure-the-ceiling-first). Supplied by the sweep as the incoming belief; at
    # init (None) use the structural-neutral default: f_g=½ with the remaining ½ split among the LIVE
    # strands (single-strand → ½ on its strand; AMBIG → ¼ each). Location is prior/likelihood-set; the
    # reference only fixes the variance/precision.
    if fg_ref is None or fpos_ref is None or fneg_ref is None:
        nlive = ap_all.astype(np.float64) + an_all.astype(np.float64)
        half = np.where(nlive > 0.0, 0.5 / np.maximum(nlive, 1.0), 0.0)
        fg_ref = np.full(m, 0.5)
        fpos_ref = np.where(ap_all, half, 0.0)
        fneg_ref = np.where(an_all, half, 0.0)
    else:
        fg_ref = np.asarray(fg_ref, np.float64)
        fpos_ref = np.asarray(fpos_ref, np.float64)
        fneg_ref = np.asarray(fneg_ref, np.float64)
    out = {k: np.zeros(m, dtype=np.float64) for k in ("fg", "fp", "fn", "vg", "vp", "vn")}
    # Skip EMPTY regions — no per-strand counts AND no unspliced/spliced mass. Both per-class solvers zero
    # every output for an inactive region (gdna/rna_mass = f_g·M = (1−f_g)·M + S = 0 when all are 0), so an
    # empty region's solve is identical to the zero-initialized `out` — skipping is BIT-IDENTICAL. At genome
    # scale most region/boundary regions carry no fragments (unexpressed genes, intergenic deserts), so this
    # is the dominant cost saver, not a slice artifact. (A spliced-only region has signal ⇒ still solved.)
    signal = (
        np.asarray(u_pos, np.float64)
        + np.asarray(u_neg, np.float64)
        + np.asarray(mass_unspl, np.float64)
        + np.asarray(mass_spliced, np.float64)
    ) > 0.0
    ss = _single_strand_mask(allow_pos, allow_neg) & signal
    amb = _ambig_mask(allow_pos, allow_neg) & signal

    def _s(a, msk):
        return None if a is None else np.asarray(a)[msk]

    def _sp(pair, msk):
        return None if pair is None else (np.asarray(pair[0])[msk], np.asarray(pair[1])[msk])

    def _scatter(msk, dc):
        out["fg"][msk] = dc.gdna_frac
        out["fp"][msk] = dc.rna_pos_frac
        out["fn"][msk] = dc.rna_neg_frac
        out["vg"][msk] = dc.gdna_frac_var
        out["vp"][msk] = dc.rna_pos_frac_var
        out["vn"][msk] = dc.rna_neg_frac_var

    if bool(ss.any()):
        # Single-strand regions solve on the FINE 1-D grid (Fix 3, n_grid_ss); the coarse-grid global prior is
        # regridded onto it. AMBIG keeps the coarse n_grid (the expensive 2-D cube). n_grid_ss=None ⇒ n_grid.
        # Tiled into row blocks for the same reason as the AMBIG cube below — see `_SOLVE_BLOCK_BYTES`. The
        # regrid rides inside the loop, so its (block, K_ss) temporaries are tiled too.
        k_ss = int(n_grid_ss) if n_grid_ss else int(n_grid)
        ss_idx = np.where(ss)[0]
        rows = _block_rows(k_ss, 8)
        for s0 in range(0, ss_idx.size, rows):
            bidx = ss_idx[s0 : s0 + rows]
            _scatter(
                bidx,
                _solve_regions_logodds(
                    u_pos[bidx],
                    u_neg[bidx],
                    allow_pos[bidx],
                    allow_neg[bidx],
                    fg_ref[bidx],
                    fpos_ref[bidx],
                    fneg_ref[bidx],
                    kappa=kappa,
                    od_g=od_g,
                    od_r=od_r,
                    n_grid=k_ss,
                    L=L,
                    priors=(priors or _NO_PRIORS).select(bidx).regrid(n_grid, k_ss, L),
                    gdna_imp_mode=_s(gdna_imp_mode, bidx),
                    gdna_imp_prec=_s(gdna_imp_prec, bidx),
                    rna_imp_mode=_sp(rna_imp_mode, bidx),
                    rna_imp_prec=_sp(rna_imp_prec, bidx),
                    rna_one_sided=_s(rna_one_sided, bidx),
                    lam_logprior=_regrid_global(_s(lam_logprior, bidx), n_grid, k_ss, L),
                    lam_imp_mode=_s(lam_imp_mode, bidx),
                    lam_imp_prec=_s(lam_imp_prec, bidx),
                ),
            )
    if bool(amb.any()):
        # The 2-D (λ,τ) cube is (B,K,K_t); materialized for ALL ambig regions at once it is ~O(m·K²) (the
        # memory the lattice OOM'd on). AMBIG regions solve independently, so tile the subset into row blocks
        # — bit-identical results, peak memory bounded to one (rows, K, K_t) cube.
        amb_idx = np.where(amb)[0]
        rows = _block_rows(int(n_grid) * (int(n_tilt) if n_tilt else int(n_grid)), 4)
        for s0 in range(0, amb_idx.size, rows):
            bidx = amb_idx[s0 : s0 + rows]
            _scatter(
                bidx,
                _solve_ambig_logodds(
                    u_pos[bidx],
                    u_neg[bidx],
                    fg_ref[bidx],
                    fpos_ref[bidx],
                    fneg_ref[bidx],
                    kappa=kappa,
                    od_g=od_g,
                    od_r=od_r,
                    n_grid=n_grid,
                    L=L,
                    n_tilt=n_tilt,
                    priors=(priors or _NO_PRIORS).select(bidx),
                    gdna_imp_mode=_s(gdna_imp_mode, bidx),
                    gdna_imp_prec=_s(gdna_imp_prec, bidx),
                    rna_imp_mode=_sp(rna_imp_mode, bidx),
                    rna_imp_prec=_sp(rna_imp_prec, bidx),
                    rna_one_sided=_s(rna_one_sided, bidx),
                    lam_logprior=_s(lam_logprior, bidx),
                    lam_imp_mode=_s(lam_imp_mode, bidx),
                    lam_imp_prec=_s(lam_imp_prec, bidx),
                    theta_imp_mode=_s(theta_imp_mode, bidx),
                    theta_imp_prec=_s(theta_imp_prec, bidx),
                ),
            )
    return RegionDeconv(
        gdna_frac=out["fg"],
        rna_pos_frac=out["fp"],
        rna_neg_frac=out["fn"],
        gdna_frac_var=out["vg"],
        rna_pos_frac_var=out["vp"],
        rna_neg_frac_var=out["vn"],
    )
