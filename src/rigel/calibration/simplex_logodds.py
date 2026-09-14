"""The log-density per-slot solver on the ``(λ, θ)`` grid — ONE solve for every slot class, driving
``sweep.solve_chain``.

The latent magnitude dof is the gDNA-vs-RNA log-odds ``λ = logit(f_g) = log ρ_g − log ρ_rna``: log-odds
bounds the 5–6-decade ρ_g range and resolves both the ``f_g→0`` and the ``f_g→1`` vertex, which a uniform
linear lattice cannot. ``λ`` is gridded on a FIXED ``[−L, L]`` window (no region-adaptivity) and the
linear fraction is read out as ``f_g = σ(λ)``. ``O(m·K)`` per region, so it is genome-scale tractable.

The ``ψ`` integrand is ``strand + (gDNA arm) + (RNA arm) + the λ-factor rows``, where each arm is that
component group's fitted log-rate prior when there is one, else the Jeffreys reference ``+½·log f``
(``_JEFFREYS_REF``).

Three facts that determine this file's shape:

1. Omitting a component's term is not "no prior" — the grid's own measure supplies one. A bare
   uniform-λ grid IS Haldane per component ⇒ Beta(0,0) on the composition: improper at BOTH vertices, a
   vertex amplifier. There is no third option.
2. The composition is a TWO-GROUP split on the λ axis — gDNA against RNA-total — which is what
   calibration models. The per-strand tilt is a nuisance parameter. On the two-group axis the measure
   residual is exactly zero: each ``logP`` is a density in LOG-rate, so its linear-rate conversion
   ``−log ρ`` cancels ``log σ'(λ)`` exactly, once per group. No Jacobian is written.
3. The tilt is gridded as ``θ = arcsin(τ)``, not as ``τ``. The Berger–Bernardo reference prior for this
   model (``f_g`` of interest, ``τ`` nuisance — the two are information-ORTHOGONAL, ``I_{f_g,τ} = 0``
   exactly) has a ``(1−τ²)^{−½}`` tilt conditional. Under ``θ = arcsin(τ)`` the Jacobian
   ``|dτ/dθ| = cos θ = (1−τ²)^{½}`` cancels it identically, so the tilt term is exactly 0 and the
   reference collapses to ONE expression for both region classes:
   ``ψ_ref = ½·log f_g + ½·log(1−f_g)``. No class branch, no endpoint singularity, no measure weights.
   θ is to the tilt what λ is to ``f_g``: the coordinate the geometry asks for. *(The vanishing is a
   property of this reference specifically — a Dirichlet(½,¼,¼) reference would leave a residual
   ``−¼·log(1−τ²)``.)*
4. The θ nodes are not a fixed lattice: at fixed λ the strand term is an exact Gaussian in τ whose θ
   peak narrows as ``n^{−½}`` (0.005 rad at 50k fragments), so ψ places its nodes across each
   ``(slot, λ)``'s own peak and weights them as the trapezoid rule (:func:`_tilt_window`) — the
   marginal is then exact at every depth with a DERIVED node count (``_TILT_NODES``), where a fixed
   lattice's sum was a comb. The weights written are the quadrature's (``log h``) and, with 5, the
   continuum's share of the reference mass (``−log π``) — never a tilt density.
5. The AMBIG tilt's hypothesis space is {pure +, pure −, mixed} at equal reference weight — THE TILT
   ATOM. Presence per strand is discrete, so beside the continuum's nodes ψ carries two columns at
   ``τ = ±1`` exactly (a single-strand solve inside the cube, no tilt parameter), and the continuum's
   weights carry ``−log π`` so the three hypotheses' masses are equal wherever the strand term is flat.
   At a strand-pure slot the atom explains the split with no parameter and the read-out lands at the
   cap, where the continuum's median sat below it; a held RNA level on a strand (a `CubeRow` profile) is
   a certified witness of that strand's RNA and rules the OTHER strand's atom out. Gated in
   ``tests/calibration/test_vertex_reference.py``; the cost at an unwitnessed both-strand slot is the
   atom's and is recorded where it was measured.

There is NO spliced term: ``mass_spliced`` is consumed only by the returned ``rna_mass``, never by ψ. That
is correct — at a sj mature RNA *splices*, so the unspliced crossing mass is gDNA plus RNA that has not
spliced there, a channel genuinely disjoint from the (directly observed, already-pure-RNA) spliced mass.

One solver, :func:`_solve_logodds`, over the ``(λ, θ)`` cube in float64. A single-strand region (exactly
one of ``allow_pos`` / ``allow_neg``) has its tilt fixed by its live strand, so it is the ``K_t = 1`` case —
a 1-D solve over ``λ`` at the 1-D cost — and AMBIG regions (both set) marginalise the tilt on the θ grid.
``_solve_regions_logodds_all`` runs the two classes on ONE λ lattice (ruled 2026-09-13: a finer
single-strand grid with a regrid between the two measured worse than one lattice) and tiles the rows so the
working set stays in cache. Structurally RNA-free regions (neither
strand live — intergenic / TSS / TES) have no composition dof and never reach the solver:
``sweep.solve_chain`` gates them out via ``solvable``, so no reference is applied to a region whose
composition is known structurally.
"""

from __future__ import annotations

import numpy as np
from scipy.special import expit, log_expit

from dataclasses import dataclass

from .region_chain import RegionDeconv

# Public surface consumed by sweep / messages / region_geometry. The remaining private helpers stay importable
# for tests but are not part of the module's external API.
__all__ = [
    "CubeRow",
    "_logodds_grid",
    "_solve_regions_logodds_all",
    "strand_row_logodds",
]

_EPS = 1.0e-9

# The reference exponent for an UNFITTED component group, as a density in LOG-rate.
#
# DERIVED, not tuned, by two agreeing routes: (a) Jeffreys for a Poisson rate — `g ~ Poisson(ρE)` ⇒
# `I(ρ) ∝ 1/ρ` ⇒ `p(ρ) ∝ ρ^(−½)`, which as a LOG-rate density is `ρ^(+½)` ⇒ `+½·log f`; (b) the
# Berger–Bernardo reference prior for the composition with `f_g` of interest and the tilt as nuisance, whose
# `f_g` marginal is Beta(½,½) — the SAME `+½·log f_g + ½·log(1−f_g)`.
#
# Its ONLY job is to make ψ proper (Beta(½,½) integrates; Beta(0,0) does not). A fitted `logP_g`/`logP_r`
# is ADDED to it, never substituted for it — the reference is the MEASURE ψ is written against, not an
# information claim to be superseded.
#
# A declared choice, not forced by the likelihood: the observed-data Fisher information for f_g is
# `∝ n(½−κ)²`, exactly 0 on an unstranded library, where the strand term is bit-flat and the posterior
# simply IS this reference. Its known cost is that it forbids the simplex vertices, where some truth
# genuinely lives.
_JEFFREYS_REF = 0.5

# f_g ∈ [σ(−10), σ(10)] = [4.5e-5, 1−4.5e-5]. A pure STATE-SPACE bracket: the widest f_g the grid can
# represent, NOT an accuracy knob — but that is a property of a PROPER ψ, not of this constant. It holds
# because both arms are always written (`_JEFFREYS_REF`): under Beta(½,½) a fraction of a percent of the
# reference's mass lies outside L=10, and the answer is L-invariant. An improper ψ (either arm omitted)
# has plateau mass growing linearly in L, and then L silently sets the prior strength.
# L-invariance is the acceptance test for a prior-free ψ, where it holds to seven digits.
#
# ⛔ It is scoped to prior-free ψ, because with a FITTED landscape installed the pipeline fails it:
# widening only the bracket (at fixed lattice spacing) moves the answer, and a resolution-only control
# moves it the other way, so the effect is the bracket and not the lattice. The mechanism is the fitted
# prior — `landscape.logprior` evaluates at `log rho = log f_c + log M − log E` and ψ can only offer
# `f_c ∈ [σ(−L), σ(L)]`, so on a gDNA-poor library σ(−10) sits well ABOVE the density the prior points
# at and the low end of the bracket is a wall the prior pushes against rather than empty state space.
# `landscape.required_logodds_window` is the derived demand. Do not read any of this as licence to widen
# L here: nothing in this file has priced what a wider bracket costs elsewhere.
#
# NB: production does not read this default — `sweep.solve_chain` threads `logodds_window` (=10.0)
# explicitly, from `CalibrationConfig.sweep_logodds_window`.
_DEFAULT_L = 10.0

# Cache-tiling target for BOTH per-region solves, as a working-set size rather than a row count — the per-row
# footprint differs ``K_t + 2`` = 26× between a single-strand row (``K``) and an AMBIG row (``K·(K_t + 2)``),
# both float64, so no single row count serves both. `_block_rows` turns it into rows.
#
# NOT a model parameter. Every region solves independently and every reduction in both solvers is within
# a row (the ψ logsumexp, the moment sums, the CDF cumsum, and the `post @ log f` gemv), so the block size
# cannot reach the arithmetic — verified bitwise for all five reduction kinds across a wide range of block
# sizes. It is purely a memory knob, and at genome scale it is the dominant one: unblocked, the 1-D path's
# intermediates are hundreds of MB each, about ten of them live, streamed from DRAM. Blocking makes the
# same arithmetic run out of cache.
_SOLVE_BLOCK_BYTES = 1 << 20


def _block_rows(cells_per_row: int, itemsize: int) -> int:
    return max(1, _SOLVE_BLOCK_BYTES // max(1, int(cells_per_row) * int(itemsize)))


def _row_moment(post, g):
    """``Σ_k post[i, k]·g[k]`` per row — a grid moment, as a per-row pairwise sum and never as a BLAS
    matrix-vector product, whose kernel (and rounding) changes with the row count. That is what makes
    ψ's read-out CHUNK-EXACT: the same row gives the same bits whatever rows share the call."""
    return np.sum(post * g[None, :], axis=1)


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
        r = m + np.log(np.sum(np.exp(a - m), axis=axis, keepdims=True))
    return r if keepdims else np.squeeze(r, axis=axis)


# ──────────────────────────────────────────────────────────────────────────────────────────────────────
# THE THREE-COMPONENT STRAND LIKELIHOOD. Its two-component special case lives in
# `strand_likelihood.strand_loglik` and is an executable REFERENCE, not dead code:
# `test_strand_likelihood_reference.py` gates that this generalization collapses onto it when one RNA
# strand is dead, so the two cannot drift apart unnoticed.
# ──────────────────────────────────────────────────────────────────────────────────────────────────────


def _mixture_strand_loglik(
    u_pos, n, f_g, f_pos, f_neg, kappa, od_g, od_r, f_g_ref, f_pos_ref, f_neg_ref
):
    """Three-component gDNA/RNA₊/RNA₋ strand loglik — :func:`strand_loglik` generalized to two RNA strands.

    Broadcasts ``(u_pos, n)`` of shape ``(regions, 1)`` against the lattice ``(f_*)`` of shape
    ``(1, P)`` → ``(regions, P)``. Mean ``N·p`` with ``p = ½·f_g + κ·f₊ + (1−κ)·f₋``.

    The count-zero-information freeze: the mean stays live in the solved composition
    ``(f_g, f_pos, f_neg)``, which is the legitimate strand channel, but the variance is evaluated at the
    fixed reference composition ``(f_g_ref, f_pos_ref, f_neg_ref)`` (per-region scalars, broadcast). That
    keeps the heteroscedastic precision — the count still sets a composition-aware variance through the
    reference — while removing the ``f_g``-tilt of the normalizer, so a raw count can no longer
    manufacture a composition preference toward the variance-minimum when the mean degenerates (κ→½). The
    reference is a neutral structural default at init and the incoming belief in the sweep.
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


def strand_row_logodds(lam, u_pos, u_neg, live_pos, kappa, od_g, od_r, f_ref):
    """One single-strand region's own strand log-likelihood over the log-odds grid ``lam``
    (``f_g = sigma(lam)``, the live RNA strand carrying ``1 - f_g``), max-normalised.

    The same term the local solve uses (:func:`_mixture_strand_loglik`), with the same
    count-zero-information freeze: the variance is evaluated at the reference composition ``f_ref``, the
    slot's incoming belief, so at ``kappa = 1/2`` the mean is constant in ``f_g`` and the row is exactly
    flat — an unstranded library says nothing, structurally. ``live_pos`` names the live strand: a ``-``
    region's RNA reads sense at rate ``1 - kappa``. This is the row a region publishes about itself, and
    the exon-to-boundary message consumes it; it carries no prior and no message."""
    lam = np.asarray(lam, np.float64)
    fg = expit(lam)
    f_ref = float(np.clip(f_ref, _EPS, 1.0 - _EPS))
    zero = np.zeros_like(fg)
    f_pos, f_neg = ((1.0 - fg), zero) if live_pos else (zero, (1.0 - fg))
    ref_pos, ref_neg = ((1.0 - f_ref), 0.0) if live_pos else (0.0, (1.0 - f_ref))
    row = _mixture_strand_loglik(
        float(u_pos),
        float(u_pos) + float(u_neg),
        fg,
        f_pos,
        f_neg,
        float(kappa),
        float(od_g),
        float(od_r),
        f_ref,
        ref_pos,
        ref_neg,
    )
    return row - row.max()


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


def _posterior_median_fg(post, lam):
    """Per-region point estimate of ``f_g``: the posterior's ½-QUANTILE, read off the CDF.

    Transform-invariant and robust to the skew of the ``f_g`` posterior, which is why ``f_g`` is a median
    and not a mean: the mean is measurably worse at both simplex vertices at every depth, and the vertices
    are where most of the in-scope calibration error lives.

    It is a continuous quantile, not the grid point where the CDF first reaches ½. The grid mass is
    treated as a histogram, with bin edges at the midpoints, and the crossing bin is interpolated.
    Returning ``fg[(cw < 0.5).sum()]`` instead snaps to a lattice point and so carries up to half a grid
    step of quantisation, which is not noise to be averaged away:

    * on an evidence-free object the posterior IS the reference — symmetric, so its ½-quantile is exactly
      ½ — and the snapped form returns the adjacent grid point instead. The exactness here is the symmetry
      of ``σ``: ``σ(δ) + σ(−δ) = 1``, so the edge between the two central bins is exactly 0.5;
    * at ``κ = ½`` the strand term is bit-flat, so every slot is in that state and the snap is a constant
      offset that does not shrink with depth — the whole closure defect there;
    * and it propagates: `_compose` builds the RNA fractions from ``1 − f_g``, so a snapped ``f_g`` snaps
      the whole composition.

    This is not a sub-grid MODE, and the distinction is the point: a mode is an argmax and can chase a
    single spike, which is what would under-call a skewed or vertex-near posterior. A quantile cannot —
    it is monotone in the CDF and reads the same half-mass crossing the snapped form approximates.

    The interpolation is done on λ, not on ``f_g``, and that is what keeps it transform-invariant.
    ``σ`` is monotone, so the ½-quantile in ``λ`` maps through it to the ½-quantile in ``f_g`` exactly:
    median equivariance, the whole reason ``f_g`` is a median. Interpolating in ``f``-space looks
    equivalent and is not — the σ grid is highly non-uniform, spacing ~1e-5 at the ends against ~0.085 in
    the middle, so a bin's midpoint in ``f`` is not the image of its midpoint in ``λ`` and a posterior
    concentrated on one grid point comes back biased toward ½, where on ``λ`` it returns its own grid
    point to machine precision. That case
    is not synthetic: an unsolved slot's fed-back belief produces a one-hot posterior.

    ``post``: (m,K) normalized posterior; ``lam`` the uniform log-odds grid. Returns (m,)."""
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

    The composition has TWO degrees of freedom, not three, and this is where that is written down. ψ
    solves a point on the 2-simplex, parametrised by

        λ  — the gDNA-vs-RNA LEVEL, read out as ``f_g`` (the posterior median over the λ grid)
        θ  — the RNA-internal TILT, a SHARE with no absolute scale, read out as ``w_pos``

    and the composition is their image::

        f_pos = (1 − f_g)·w_pos        f_neg = (1 − f_g)·(1 − w_pos)

    Closure is therefore structural and cannot fail, because the map lands on the simplex:
    ``f_g + f_pos + f_neg = 1`` identically, to float64 rounding, on both ψ paths at every κ and depth, so
    no consumer has to check it and no arithmetic has to be renormalised.

    Reading the three coordinates out independently does not close. Taking ``f_g`` as the posterior
    median over λ while ``f_pos``/``f_neg`` are posterior MEANS of the grid quantity ``1 − f_g`` mixes a
    quantile with expectations:

        SUM = median(f_g) + (1 − mean(f_g)) = 1 + median(f_g) − mean(f_g)

    so the closure error is exactly the posterior's SKEW, and most real objects miss.

    The repair is not "take means everywhere". That closes too, by linearity of expectation, and is
    measurably worse: the median is closer to the truth at both simplex vertices, where most of the
    in-scope error lives. Nor is it renormalising three numbers at
    publication, which would make a badly short object indistinguishable from a solved one. Nothing here
    is rescaled: ``f_g`` and the RNA total are exact complements *by parametrisation*, and the tilt is
    estimated as a share because a share is what it is.

    ``w_pos`` is the + strand's share of RNA — the RNA-mass-weighted posterior share on an AMBIG slot.

    Admissibility is enforced here, not by the caller, because a share that ignores it loses RNA
    silently: with only the + strand free, ``w_pos = ½`` would place half the RNA on a forbidden strand,
    where it is zeroed and simply vanishes. The tilt of a single-strand slot is structurally locked, so
    the admissible strand takes the whole RNA total whatever ``w_pos`` says. A slot with neither strand
    admissible has no RNA to place and returns ``(0, 0)``: its composition is ``f_g`` alone, which is the
    honest statement and not a closure failure, and nothing dispatches such a slot to a ψ solve anyway.
    """
    fr = 1.0 - np.asarray(f_g, np.float64)
    ap = np.asarray(allow_pos, bool)
    an = np.asarray(allow_neg, bool)
    # ``w_pos`` is a share, so it is clamped to [0,1] here rather than trusted. Unclamped, a share
    #   outside the range produces a negative fraction that still sums to 1 — a composition that passes a
    #   closure check and is nonsense. It is not reachable from the AMBIG caller, being a ratio of two
    #   non-negative expectations, which is exactly why the constraint belongs where it can be seen.
    w = np.clip(np.asarray(w_pos, np.float64), 0.0, 1.0)
    w = np.where(ap & an, w, np.where(ap, 1.0, 0.0))
    return np.where(ap, fr * w, 0.0), np.where(an, fr * (1.0 - w), 0.0)


def _gdna_arm(lam, global_logprior):
    """The gDNA arm of ψ over the λ grid → broadcastable to ``(m, K)``.

    ``_JEFFREYS_REF`` reference ``+½·log f_g`` → ``(1, K)``. When a prior is fitted, the fitted
    ``logP_g(log ρ_g)`` (pre-evaluated on THIS ``f_g`` grid → ``(m, K)``) is ADDED to it.

    ⛔ A fitted prior must never REPLACE the reference. The reference is not an information claim to be
    superseded, it is the measure ψ is written against, and it is the only term bounding this arm at
    ``f_g → 0`` (`_rna_arm` bounds the other vertex and is never replaced, so replacing here would not
    even treat the two arms alike). Dropping it leaves ψ improper at exactly the vertex a fitted gDNA
    prior most often points at. Bayes composes a prior with a measure by addition in log space; there is
    no double-count to avoid.

    ``None`` means "not fitted", not "no term"."""
    ref = _JEFFREYS_REF * _log_fg(lam)[None, :]
    if global_logprior is None:
        return ref
    return ref + np.asarray(global_logprior, np.float64)


def _rna_arm(lam):
    """The RNA-total arm of ψ over the λ grid → ``(1, K)``: the ``_JEFFREYS_REF`` reference
    ``+½·log(1 − f_g)``, the exact mirror of :func:`_gdna_arm`'s reference.

    This is the two-group arm (gDNA vs RNA-total): the per-strand split is the nuisance tilt, integrated
    out on the θ axis, and needs no prior of its own.

    Nothing fits ``logP_r``, and the cost of that is known: the reference alone bounds the ``f_g → 1``
    vertex, and unlike its gDNA twin it is never swamped by evidence — a fixed repulsion of about 3.1
    nats at ``f_g = 0.999`` relative to ``f_g = ½``, roughly a 22:1 handicap. Objects whose true ``f_g``
    sits at that vertex carry most of the calibration error on the in-scope strata and read below the
    vertex. A fitted RNA arm would land here as a second argument, the mirror of ``global_logprior``;
    the unfed socket for it was removed (2026-09-13) rather than carried."""
    return _JEFFREYS_REF * _log1m_fg(lam)[None, :]


@dataclass(frozen=True)
class CubeRow:
    """The RNA level lanes' delivery at an AMBIG slot, as what it is made of: per strand the held level
    profile over ``u = log(ρ/ρ_ref)`` on the solve grid (either may be absent), the slot's total and its
    RNA opportunity, and the lanes' one reference density (a level is absolute; the coordinate is only
    its origin, and both strands share it). ψ evaluates it at its own θ nodes (:meth:`at`) —
    there is no θ lattice for a row to be built on, so nothing is interpolated. A one-sided profile stays
    one-sided through the map (it is monotone in each share), so "at least this much RNA+" arrives as a
    wall in the cube and no parametric summary is made."""

    profile_pos: np.ndarray | None
    profile_neg: np.ndarray | None
    u: np.ndarray
    total: float
    opportunity: float
    rho_ref: float

    def at(self, fg, tau) -> np.ndarray:
        """The row over ψ's cells: at each ``(λ, θ)`` the strand's share ``f_s = (1 − f_g)(1 ± τ)/2``
        implies the density ``f_s·n/a_r``, and the held profile is read at ``log(ρ_s/ρ_ref_s)`` — the
        `profile_of_level` map on the λ axis with the tilt inside. ``fg`` is ``(K,)``; ``tau`` is
        ``(K_t,)`` or ``(K, K_t)``; the result is ``(K, K_t)``, max-normalised."""
        fg = np.asarray(fg, np.float64)
        tau = np.asarray(tau, np.float64)
        if tau.ndim == 1:
            tau = tau[None, :]
        f_act = (1.0 - fg)[:, None]
        out = np.zeros(np.broadcast_shapes(f_act.shape, tau.shape))
        scale = float(self.total) / float(self.opportunity)
        for prof, sign in ((self.profile_pos, 1.0), (self.profile_neg, -1.0)):
            if prof is None:
                continue
            prof = np.asarray(prof, np.float64)
            with np.errstate(divide="ignore"):
                u_s = np.log(f_act * (1.0 + sign * tau) / 2.0 * scale) - np.log(float(self.rho_ref))
            out += np.interp(
                u_s, np.asarray(self.u, np.float64), prof, left=prof[0], right=prof[-1]
            )
        return out - out.max()


# The θ quadrature's truncation: the strand term's mass outside a window is below double precision.
# DERIVED, not tuned — ``erfc(√T) ≈ e^{−T}/√(πT) < ε₆₄`` at ``T = −log ε₆₄`` — so nothing about the window
# is a choice; a wider one adds nodes where the integrand is zero to the last bit.
_T_NATS = -np.log(np.finfo(np.float64).eps)

# The node count that resolves the peak inside its window. The interior window is ``2√(2T)·σ_θ`` wide and
# the trapezoid rule's error on a Gaussian of width ``σ_θ`` at spacing ``h`` is ``2·e^{−2π²(σ_θ/h)²}``,
# below ``e^{−T}`` once ``h ≤ σ_θ·π√2/√T`` — so ``K_t − 1 ≥ 2T/π``. At a domain end the window is only
# ``T^{¼}`` widths per side and this count over-resolves it. DERIVED from ``_T_NATS``; the gate
# `test_vertex_reference.test_the_derived_node_count_is_converged` reads it out against 60 nodes.
_TILT_NODES = int(np.ceil(2.0 * _T_NATS / np.pi)) + 1


def _tilt_window(u_pos, n, fg, kappa, od_g, od_r, f_g_ref, f_pos_ref, f_neg_ref):
    """The θ nodes and log-weights that integrate the strand term EXACTLY at every depth: per
    ``(slot, λ)`` a window across the term's peak, ``_TILT_NODES`` uniform nodes in θ inside it, the
    trapezoid weights. Returns ``(theta (m, K, K_t), log_weight (m, K, K_t))``.

    At fixed λ the strand term is an exact Gaussian in the tilt ``τ`` (the variance is frozen at the
    reference): ``p = ½ + a(λ)·τ`` with ``a = (1 − f_g)(κ − ½)``, so its centre is ``τ̂ = d/a`` for the
    observed contrast ``d = u₊/n − ½`` and its width ``σ_τ = σ_p/|a|`` with ``σ_p = √V/n``. Its θ peak is
    ``σ_τ/cos θ̂`` wide — 0.005 rad at 50k fragments against a 60-node lattice's 0.053 step — so a fixed
    lattice's sum is a comb across λ, a factor between 1 and ``e^{−(h/2)²/2σ_τ²}`` chosen by where the peak
    happens to fall, and a deep AMBIG slot's λ posterior is a set of spikes at arbitrary λ (the ladder's
    recorded K_t 30 failure was one 25k-fragment slot whose peak fell between nodes; 60 held it by a coin
    toss). The window is where the term lies within ``_T_NATS`` of its maximum ON the domain,

        τ_m = clip(τ̂, −1, 1),   ρ = √((τ_m − τ̂)² + 2σ_τ²T),   [τ_lo, τ_hi] = [τ̂ − ρ, τ̂ + ρ] ∩ [−1, 1],

    one closed form for an interior peak, a peak on the boundary and a peak beyond it. In ``θ = arcsin τ``
    the integrand is smooth and EVEN about ``±π/2``, so the trapezoid rule with weight ½ at both window
    ends is spectrally accurate in every regime: at an interior end the integrand is ``e^{−T}`` of its
    peak and the weight is immaterial; at a domain end it is exactly the periodic trapezoid rule's weight
    on the reflected window. The node spacing ``h(slot, λ)`` enters ψ as ``log h`` — the window scales with
    ``σ_τ(λ) ∝ 1/(1 − f_g)``, and that scaling is the λ-dependence the fixed lattice got wrong. A slot
    with no strand information (``κ = ½``, ``n = 0``, or a window wider than the domain) gets the whole
    domain: the uniform lattice with trapezoid weights, exactly where the lattice was right.
    """
    m, K, Kt = u_pos.shape[0], fg.shape[0], _TILT_NODES
    # the frozen variance, the same expression `_mixture_strand_loglik` evaluates (per slot)
    p_ref = 0.5 * f_g_ref + kappa * f_pos_ref + (1.0 - kappa) * f_neg_ref
    rscale = kappa * (1.0 - kappa)
    var = (
        n * p_ref * (1.0 - p_ref)
        + (n * f_g_ref) ** 2 * 0.25 * od_g
        + (n * f_pos_ref) ** 2 * rscale * od_r
        + (n * f_neg_ref) ** 2 * rscale * od_r
    )
    var = np.maximum(var, _EPS)
    counted = n > 0.0
    n_safe = np.where(counted, n, 1.0)
    d = np.where(counted, u_pos / n_safe - 0.5, 0.0)  # (m,)
    sig_p = np.where(counted, np.sqrt(var) / n_safe, np.inf)
    a = (1.0 - fg)[None, :] * (kappa - 0.5)  # (1, K)
    informative = counted[:, None] & (a != 0.0)  # (m, K)
    a_safe = np.where(a != 0.0, a, 1.0)
    tau_hat = np.where(informative, d[:, None] / a_safe, 0.0)
    with np.errstate(over="ignore", invalid="ignore"):
        sig_tau = np.where(informative, sig_p[:, None] / np.abs(a_safe), np.inf)
        rho = np.sqrt((np.clip(tau_hat, -1.0, 1.0) - tau_hat) ** 2 + 2.0 * sig_tau**2 * _T_NATS)
    tau_lo = np.where(informative, np.maximum(-1.0, tau_hat - rho), -1.0)
    tau_hi = np.where(informative, np.minimum(1.0, tau_hat + rho), 1.0)
    th_lo = np.arcsin(tau_lo)
    h = (np.arcsin(tau_hi) - th_lo) / (Kt - 1)  # (m, K)
    theta = th_lo[:, :, None] + h[:, :, None] * np.arange(Kt, dtype=np.float64)[None, None, :]
    log_weight = np.broadcast_to(np.log(h)[:, :, None], (m, K, Kt)).copy()
    log_weight[:, :, 0] -= np.log(2.0)
    log_weight[:, :, -1] -= np.log(2.0)
    return theta, log_weight


def _single_strand_mask(allow_pos, allow_neg) -> np.ndarray:
    """The regions the 1-D solver is valid for: exactly one strand live, so the tilt is determined."""
    ap = np.asarray(allow_pos, bool)
    an = np.asarray(allow_neg, bool)
    return ap ^ an


def _ambig_mask(allow_pos, allow_neg) -> np.ndarray:
    """AMBIG regions (both strands live) — the 2-D ``(λ, θ)`` path."""
    return np.asarray(allow_pos, bool) & np.asarray(allow_neg, bool)


def _psi(
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
    lam,
    fg,
    ambig: bool = False,
    gdna_logprior=None,
    lam_logprior=None,
    cube_rows=None,
):
    """ψ over the ``(λ, θ)`` cube for ``m`` slots — strand + ``_gdna_arm`` + ``_rna_arm`` + the λ-factor
    rows (+ the delivered cube rows) + the θ quadrature's log-weights — as ``(m, K, K_t)`` in float64,
    with the two strand-fraction grids it was evaluated on, ``(f_pos, f_neg)``, and the tilt ``tau`` they
    were built from. A single-strand call (``ambig=False``) has the tilt of each slot's live strand
    (``τ = ±1``), a ``(m, K, 1)`` cube and no weight. An AMBIG call places the θ nodes across each slot's
    strand term (:func:`_tilt_window`, ``K_t = _TILT_NODES``) for the MIXED hypothesis, appends the two
    PURE hypotheses as the columns ``τ = +1, −1`` (the tilt atom: the three at equal reference weight,
    a delivered level on a strand ruling the other strand's atom out) and evaluates each delivered
    :class:`CubeRow` at every column (``cube_rows``: one record or ``None`` per slot), so the cube is
    ``(m, K, K_t + 2)``. :func:`_solve_logodds` reads ψ out; the vertex-reference gates read ψ itself."""
    ap = np.asarray(allow_pos, bool)
    an = np.asarray(allow_neg, bool)
    u_pos = np.asarray(u_pos, np.float64)
    n = u_pos + np.asarray(u_neg, np.float64)
    fg = np.asarray(fg, np.float64)
    f_g_ref = np.asarray(f_g_ref, np.float64)
    f_pos_ref = np.asarray(f_pos_ref, np.float64)
    f_neg_ref = np.asarray(f_neg_ref, np.float64)
    if not ambig:
        # the tilt of a single-strand slot is its live strand: the cube is (m, K, 1)
        tau = np.where(ap & ~an, 1.0, -1.0)[:, None, None]
        log_weight = None
    else:
        theta, log_weight = _tilt_window(
            u_pos, n, fg, float(kappa), float(od_g), float(od_r), f_g_ref, f_pos_ref, f_neg_ref
        )
        # THE TILT ATOM: the hypothesis space is {pure +, pure −, mixed} at equal reference weight —
        # two more columns at τ = ±1 exactly (a single-strand solve inside the cube, no tilt parameter,
        # log-weight 0) beside the mixed continuum, whose trapezoid weights carry −log π so that its
        # measure is the uniform density over the domain (its window's share of it). Where the data
        # are pure the atom explains them with no parameter and wins the Occam contest against the
        # continuum's ∝ σ_θ; where they are not it sits e^{−(1−τ̂)²/2σ_τ²} below the peak and vanishes.
        m = u_pos.shape[0]
        tau = np.concatenate(
            [np.sin(theta), np.ones((m, fg.shape[0], 1)), -np.ones((m, fg.shape[0], 1))], axis=2
        )  # τ = sin θ, exact across the domain; then the two atoms
        log_weight = np.concatenate(
            [log_weight - np.log(np.pi), np.zeros((m, fg.shape[0], 2))], axis=2
        )
    f_act = (1.0 - fg)[None, :, None]  # (1, K, 1)
    f_pos = f_act * (1.0 + tau) / 2.0  # (m|1, K, K_t)
    f_neg = f_act * (1.0 - tau) / 2.0
    psi = _mixture_strand_loglik(
        u_pos[:, None, None],
        n[:, None, None],
        fg[None, :, None],
        f_pos,
        f_neg,
        kappa,
        od_g,
        od_r,
        f_g_ref[:, None, None],
        f_pos_ref[:, None, None],
        f_neg_ref[:, None, None],
    )
    psi = psi + (_gdna_arm(lam, gdna_logprior) + _rna_arm(lam))[:, :, None]
    if lam_logprior is not None:
        psi = psi + np.asarray(lam_logprior, np.float64)[:, :, None]
    if cube_rows is not None:
        for j, row in enumerate(cube_rows):
            if row is not None:
                psi[j] += row.at(fg, tau[j])
                if ambig:
                    # THE WITNESS: a held RNA level on a strand certifies that the strand carries RNA,
                    # so the hypothesis that all of the slot's RNA is on the OTHER strand is out. A level
                    # on a strand says nothing against that strand's own atom; nothing delivered leaves
                    # both atoms standing.
                    if row.profile_neg is not None:
                        psi[j, :, -2] = -np.inf
                    if row.profile_pos is not None:
                        psi[j, :, -1] = -np.inf
    if log_weight is not None:
        psi = psi + log_weight
    return psi, f_pos, f_neg, tau


def _solve_logodds(
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
    ambig: bool = False,
    gdna_logprior=None,
    lam_logprior=None,
    cube_rows=None,
) -> RegionDeconv:
    """THE per-slot solve, one for every slot class, over the ``(λ, θ)`` grid: the gDNA-vs-RNA log-odds
    ``λ`` (outer, ``K = n_grid``) and the tilt ANGLE ``θ = arcsin(τ)`` (inner, ``K_t = _TILT_NODES`` nodes
    placed per slot and λ by :func:`_tilt_window`). ψ is the
    same expression for every slot — strand + ``_gdna_arm`` + ``_rna_arm`` + the λ-factor rows (+ the
    delivered :class:`CubeRow` where there is one) — evaluated on the ``(m, K, K_t)`` cube in float64, and read
    out once: ``f_g`` is the posterior median over the θ-marginal λ-posterior, ``f_pos`` / ``f_neg`` its
    image under :func:`_compose` with the tilt share ``w_pos`` the RNA-mass-weighted posterior share, and
    ``Var(log f_g)`` is a grid moment over the λ-marginal — the one precision the tool reads (the
    landscape prior's training weight). One exp over the cube serves both read-outs, and no
    transcendental is taken over the cube after it.

    A SINGLE-STRAND slot is the ``K_t = 1`` case and not a second solver: its tilt is fixed by its live
    strand (``τ = +1`` where only ``+`` is admissible, ``−1`` where only ``−`` is), so the caller passes
    ``ambig=False`` and the cube is ``(m, K, 1)`` — the 1-D solve over λ at the 1-D cost. Every slot of a
    call must then be single-strand; an AMBIG call (both strands admissible) takes the windowed θ nodes.
    The caller (:func:`_solve_regions_logodds_all`) separates the two classes.

    No Jacobian and no tilt term are written, and that is the point of the θ coordinate: the
    Berger–Bernardo tilt conditional ``(1−τ²)^{−½}`` is cancelled identically by ``|dτ/dθ| = cos θ``,
    and on the two-group λ axis the log-rate conversions cancel ``log σ'(λ)``; the only weights in ψ are
    the θ quadrature's own (``_tilt_window``). Both arms
    are always written — a fitted ``logP`` where there is one, else the ``_JEFFREYS_REF`` reference —
    because omitting one is not neutral (the module docstring). The strand mixture's variance is frozen
    at the reference composition (``f_g_ref`` / ``f_pos_ref`` / ``f_neg_ref``, per slot), so the count
    sets precision and not composition. Zero-mass slots report 0.

    THE READ-OUT IS CHUNK-EXACT: every reduction runs per row in a fixed order (the θ and λ logsumexps,
    the moment sums as per-row sums and never BLAS), so a slot's answer is a function of its own inputs
    and the grid and never of which rows share the call; ψ is made contiguous first, since a reduction
    over a non-contiguous array follows the strides. One ulp is a different number, and the locus solve
    tiles the rows freely (gate: ``test_sweep.test_the_psi_solve_is_chunk_exact_so_a_block_split_moves_no_number``).
    """
    lam, fg = _logodds_grid(int(n_grid), L)
    ap = np.asarray(allow_pos, bool)
    an = np.asarray(allow_neg, bool)
    n = np.asarray(u_pos, np.float64) + np.asarray(u_neg, np.float64)
    psi, f_pos, f_neg, tau = _psi(
        u_pos,
        u_neg,
        ap,
        an,
        f_g_ref,
        f_pos_ref,
        f_neg_ref,
        kappa=kappa,
        od_g=od_g,
        od_r=od_r,
        lam=lam,
        fg=fg,
        ambig=ambig,
        gdna_logprior=gdna_logprior,
        lam_logprior=lam_logprior,
        cube_rows=cube_rows,
    )
    psi = np.ascontiguousarray(psi)
    m = psi.shape[0]
    # ── ONE posterior over the cube; the λ read-out on its θ-marginal ─────────────────────────────
    flat = psi.reshape(m, -1)
    post = np.exp(flat - _lse(flat, axis=1, keepdims=True)).reshape(psi.shape)  # (m, K, K_t)
    post_lam = post.sum(axis=2)  # the θ-marginal, per row
    f_g = _posterior_median_fg(post_lam, lam)
    log_fg = _log_fg(lam)
    m_lg = _row_moment(post_lam, log_fg)
    var_g = np.maximum(_row_moment(post_lam, log_fg * log_fg) - m_lg * m_lg, 0.0)
    # ── the tilt share on the whole cube: the RNA-mass-weighted posterior share of the + strand ────
    m_pos = np.sum(post * f_pos, axis=(1, 2))
    m_neg = np.sum(post * f_neg, axis=(1, 2))
    rna = m_pos + m_neg
    w_pos = np.where(rna > 0.0, m_pos / np.where(rna > 0.0, rna, 1.0), 0.5)
    # ── the composition as the image of its two parameters ────────────────────────────────────────
    active = n > 0.0
    f_g = np.where(active, np.clip(f_g, 0.0, 1.0), 0.0)
    f_pos, f_neg = _compose(f_g, w_pos, ap, an)
    return RegionDeconv(
        gdna_frac=f_g,
        rna_pos_frac=np.where(active, f_pos, 0.0),
        rna_neg_frac=np.where(active, f_neg, 0.0),
        gdna_frac_var=np.where(active, var_g, 0.0),
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
    gdna_logprior=None,
    lam_logprior=None,
    fg_ref=None,
    fpos_ref=None,
    fneg_ref=None,
    cube_rows=None,
) -> RegionDeconv:
    """The per-region dispatcher: runs :func:`_solve_logodds` on the single-strand regions (a one-cell
    tilt) and on the AMBIG regions (the windowed θ nodes, ``K_t = _TILT_NODES``), both on the one ``λ`` lattice,
    scattering both into full-length arrays. Structurally pure-gDNA and zero-mass regions report 0, and
    ``sweep.solve_chain`` keeps their signature-binary init through the ``solvable`` write-back.

    All array inputs are full length ``m``; ``gdna_logprior`` is ``(m, K)`` on the σ(λ) grid;
    ``lam_logprior`` is ``(m, K)``. Each is sub-indexed per class.
    ``cube_rows`` is ``{slot: CubeRow}`` for AMBIG slots (the RNA level lanes' delivery), evaluated at
    each slot's own θ nodes inside its ψ; ``None`` or an absent slot changes nothing."""
    m = int(np.asarray(u_pos).shape[0])
    ap_all = np.asarray(allow_pos, bool)
    an_all = np.asarray(allow_neg, bool)
    # The count-zero-information variance-freeze reference, supplied by the sweep as the incoming belief.
    # At init (None) it is the structural-neutral default: f_g=½ with the remaining ½ split among the live
    # strands (single-strand → ½ on its strand, AMBIG → ¼ each). The location is prior- and
    # likelihood-set; the reference only fixes the variance, hence the precision.
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
    out = {k: np.zeros(m, dtype=np.float64) for k in ("fg", "fp", "fn", "vg")}
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

    def _scatter(msk, dc):
        out["fg"][msk] = dc.gdna_frac
        out["fp"][msk] = dc.rna_pos_frac
        out["fn"][msk] = dc.rna_neg_frac
        out["vg"][msk] = dc.gdna_frac_var

    if bool(ss.any()):
        # Single-strand regions: the 1-D solve (a one-cell tilt), tiled into row blocks for the same
        # reason as the AMBIG cube below — see `_SOLVE_BLOCK_BYTES`.
        ss_idx = np.where(ss)[0]
        rows = _block_rows(int(n_grid), 8)
        for s0 in range(0, ss_idx.size, rows):
            bidx = ss_idx[s0 : s0 + rows]
            _scatter(
                bidx,
                _solve_logodds(
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
                    n_grid=n_grid,
                    L=L,
                    gdna_logprior=_s(gdna_logprior, bidx),
                    lam_logprior=_s(lam_logprior, bidx),
                ),
            )
    if bool(amb.any()):
        # The 2-D (λ,θ) cube is (B,K,K_t); materialised for every AMBIG region at once it would be
        # O(m·K·K_t). AMBIG regions solve independently, so the subset is tiled into row blocks —
        # bit-identical results, peak memory bounded to one (rows, K, K_t) cube.
        amb_idx = np.where(amb)[0]
        rows = _block_rows(int(n_grid) * _TILT_NODES, 8)
        for s0 in range(0, amb_idx.size, rows):
            bidx = amb_idx[s0 : s0 + rows]
            delivered = None
            if cube_rows:
                delivered = [cube_rows.get(int(slot)) for slot in bidx]
                if all(r is None for r in delivered):
                    delivered = None
            _scatter(
                bidx,
                _solve_logodds(
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
                    n_grid=n_grid,
                    L=L,
                    ambig=True,
                    gdna_logprior=_s(gdna_logprior, bidx),
                    lam_logprior=_s(lam_logprior, bidx),
                    cube_rows=delivered,
                ),
            )
    return RegionDeconv(
        gdna_frac=out["fg"],
        rna_pos_frac=out["fp"],
        rna_neg_frac=out["fn"],
        gdna_frac_var=out["vg"],
    )
