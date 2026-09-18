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
   ``(slot, λ)``'s own peak and weights them as the trapezoid rule (`psi_kernel.cpp`'s ``tilt_window``) — the
   marginal is then exact at every depth with a DERIVED node count (``_TILT_NODES``), where a fixed
   lattice's sum was a comb. The weights written are the quadrature's (``log h``) and, with 5, the
   continuum's share of the reference mass (``−log π``) — never a tilt density.
5. The AMBIG tilt's hypothesis space is {pure +, pure −, mixed} at equal reference weight — THE TILT
   ATOM. Presence per strand is discrete, so beside the continuum's nodes ψ carries two columns at
   ``τ = ±1`` exactly (a single-strand solve inside the cube, no tilt parameter), and the continuum's
   weights carry ``−log π`` so the three hypotheses' masses are equal wherever the strand term is flat.
   At a strand-pure slot the atom explains the split with no parameter and the read-out lands at the
   cap, where the continuum's median sat below it; a held RNA level on a strand (a `CubeRows` profile) is
   a certified witness of that strand's RNA and rules the OTHER strand's atom out. Gated in
   ``tests/calibration/test_vertex_reference.py``; the cost at an unwitnessed both-strand slot is the
   atom's and is recorded where it was measured.

There is NO spliced term: ``mass_spliced`` is consumed only by the returned ``rna_mass``, never by ψ. That
is correct — at a sj mature RNA *splices*, so the unspliced crossing mass is gDNA plus RNA that has not
spliced there, a channel genuinely disjoint from the (directly observed, already-pure-RNA) spliced mass.

One solver, native (`native.psi_solve`, ``src/rigel/native/psi_kernel.cpp``), over the ``(λ, θ)`` cube in
float64, one slot at a time on the slot's own cube. A single-strand region (exactly one of ``allow_pos`` /
``allow_neg``) has its tilt fixed by its live strand, so it is the one-column case — a 1-D solve over ``λ``
at the 1-D cost — and AMBIG regions (both set) marginalise the tilt over the windowed θ nodes and the two
atoms. :func:`_solve_regions_logodds_all` is the dispatcher: the reference defaults, the signal mask, the
delivered cube rows packed, one native call, both classes on ONE λ lattice (ruled 2026-09-13: a finer
single-strand grid with a regrid between the two measured worse than one lattice). Structurally RNA-free
regions (neither strand live — intergenic / TSS / TES) have no composition dof and never reach the solver:
``sweep.solve_chain`` gates them out via ``solvable``, so no reference is applied to a region whose
composition is known structurally. The gates read ψ through the same code: :func:`psi_cube` (the cube),
:func:`posterior_median_fg` (the read-out's quantile) and :func:`compose` (the composition).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import expit

from ..native import psi_compose, psi_cube_native, psi_posterior_median, psi_solve
from .region_chain import RegionDeconv

# Public surface consumed by sweep / messages / region_geometry, and the pieces the gates read.
__all__ = [
    "CubeRows",
    "_logodds_grid",
    "_solve_regions_logodds_all",
    "compose",
    "posterior_median_fg",
    "psi_cube",
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

# Cache-tiling target for the row-tiled fits (the landscape's kernels, the capture efficiency), as a
# working-set size rather than a row count; `_block_rows` turns it into rows. NOT a model parameter: every
# reduction those fits make is within a row, so the block size cannot reach the arithmetic. It is purely a
# memory knob.
_SOLVE_BLOCK_BYTES = 1 << 20


def _block_rows(cells_per_row: int, itemsize: int) -> int:
    return max(1, _SOLVE_BLOCK_BYTES // max(1, int(cells_per_row) * int(itemsize)))


def _logodds_grid(n_grid: int, L: float = _DEFAULT_L):
    """The fixed log-odds lattice: ``λ`` uniform on ``[−L, L]`` (``K = n_grid`` points, ascending) and
    the matching ``f_g = σ(λ)`` (also ascending). Returns ``(lam, fg)``, each length ``K``."""
    lam = np.linspace(-float(L), float(L), int(n_grid))
    return lam, expit(lam)


@dataclass(frozen=True)
class CubeRows:
    """The RNA level lanes' delivery at the AMBIG slots of a block, as a TABLE — one row per delivered slot:
    per strand the held level profile over ``u = log(ρ/ρ_ref)`` on the solve grid with its presence bit
    (either may be absent), the slot's total and its RNA opportunity, and the lanes' one reference density
    (a level is absolute; the coordinate is only its origin, and both strands share it); ``u`` is the one
    grid every row is on. ψ evaluates a row at its own θ nodes — at each ``(λ, θ)`` the strand's share
    ``f_s = (1 − f_g)(1 ± τ)/2`` implies the density ``f_s·n/a_r``, and the held profile is read at
    ``log(ρ_s/ρ_ref)``, the `profile_of_level` map on the λ axis with the tilt inside — so there is no θ
    lattice for a row to be built on and nothing is interpolated in θ. A one-sided profile stays one-sided
    through the map (it is monotone in each share), so "at least this much RNA+" arrives as a wall in the
    cube and no parametric summary is made. Built by the policy's solve (`native.transfer_solve`), read by
    ψ (`native.psi_solve`); the arrays ARE the kernels' arguments."""

    slot: np.ndarray  # (d,) int64 — the block's slot each row belongs to
    profile_pos: np.ndarray  # (d, K)
    has_pos: np.ndarray  # (d,) bool
    profile_neg: np.ndarray  # (d, K)
    has_neg: np.ndarray  # (d,) bool
    u: np.ndarray  # (K,)
    total: np.ndarray  # (d,)
    opportunity: np.ndarray  # (d,)
    rho_ref: np.ndarray  # (d,)

    #: the per-row arrays, in field order
    PER_ROW = (
        "slot",
        "profile_pos",
        "has_pos",
        "profile_neg",
        "has_neg",
        "total",
        "opportunity",
        "rho_ref",
    )

    @classmethod
    def blank(cls, d: int, u) -> CubeRows:
        """``d`` rows with nothing delivered yet — the table a solve writes into."""
        u = np.ascontiguousarray(u, np.float64)
        K = u.shape[0]
        return cls(
            np.zeros(d, np.int64),
            np.zeros((d, K)),
            np.zeros(d, bool),
            np.zeros((d, K)),
            np.zeros(d, bool),
            u,
            np.zeros(d),
            np.zeros(d),
            np.zeros(d),
        )

    def __len__(self) -> int:
        return int(self.slot.shape[0])

    def select(self, keep) -> CubeRows:
        """The rows ``keep`` (a mask or index array) as a table on the same grid."""
        return CubeRows(**{f: getattr(self, f)[keep] for f in self.PER_ROW}, u=self.u)

    def shifted(self, offset: int) -> CubeRows:
        """The same rows with their slots re-keyed by ``offset`` — a block's table as the chain's."""
        return CubeRows(
            **{f: getattr(self, f) for f in self.PER_ROW if f != "slot"},
            slot=self.slot + int(offset),
            u=self.u,
        )

    @classmethod
    def concat(cls, parts: list) -> CubeRows:
        return cls(
            **{f: np.concatenate([getattr(q, f) for q in parts]) for f in cls.PER_ROW}, u=parts[0].u
        )

    @property
    def nbytes(self) -> int:
        return sum(getattr(self, f).nbytes for f in self.PER_ROW) + self.u.nbytes

    def kernel_args(self) -> dict:
        """The table as ψ's kernel takes it, by argument name."""
        return dict(
            cube_slot=self.slot,
            cube_pos=self.profile_pos,
            cube_has_pos=self.has_pos,
            cube_neg=self.profile_neg,
            cube_has_neg=self.has_neg,
            cube_u=self.u,
            cube_total=self.total,
            cube_opportunity=self.opportunity,
            cube_rho=self.rho_ref,
        )


# The θ quadrature's truncation: the strand term's mass outside a window is below double precision.
# DERIVED, not tuned — ``erfc(√T) ≈ e^{−T}/√(πT) < ε₆₄`` at ``T = −log ε₆₄`` — so nothing about the window
# is a choice; a wider one adds nodes where the integrand is zero to the last bit. The kernel carries the
# same constant (`psi_kernel.cpp`); this one derives the node count below.
_T_NATS = -np.log(np.finfo(np.float64).eps)

# The node count that resolves the peak inside its window. The interior window is ``2√(2T)·σ_θ`` wide and
# the trapezoid rule's error on a Gaussian of width ``σ_θ`` at spacing ``h`` is ``2·e^{−2π²(σ_θ/h)²}``,
# below ``e^{−T}`` once ``h ≤ σ_θ·π√2/√T`` — so ``K_t − 1 ≥ 2T/π``. At a domain end the window is only
# ``T^{¼}`` widths per side and this count over-resolves it. DERIVED from ``_T_NATS``; the gate
# `test_vertex_reference.test_the_derived_node_count_is_converged` reads it out against 60 nodes.
_TILT_NODES = int(np.ceil(2.0 * _T_NATS / np.pi)) + 1


def _cube_args(cube_rows, u) -> dict:
    """The delivered rows as ψ's kernel takes them — the table's arrays, or an empty table's on this grid."""
    if cube_rows is None:
        return CubeRows.blank(0, u).kernel_args()
    if cube_rows.u.shape != u.shape or not np.array_equal(cube_rows.u, u):
        raise ValueError("the delivered rows are not on the solve grid")
    return cube_rows.kernel_args()


def _reference_composition(allow_pos, allow_neg, fg_ref, fpos_ref, fneg_ref):
    """The count-zero-information variance-freeze reference, per slot: the incoming belief where the sweep
    supplies one; at init (``None``) the structural-neutral default — ``f_g = ½`` with the remaining ½ split
    among the live strands (single-strand → ½ on its strand, AMBIG → ¼ each). The location is prior- and
    likelihood-set; the reference only fixes the variance, hence the precision."""
    ap, an = np.asarray(allow_pos, bool), np.asarray(allow_neg, bool)
    if fg_ref is None or fpos_ref is None or fneg_ref is None:
        nlive = ap.astype(np.float64) + an.astype(np.float64)
        half = np.where(nlive > 0.0, 0.5 / np.maximum(nlive, 1.0), 0.0)
        return np.full(ap.shape[0], 0.5), np.where(ap, half, 0.0), np.where(an, half, 0.0)
    return (
        np.ascontiguousarray(fg_ref, np.float64),
        np.ascontiguousarray(fpos_ref, np.float64),
        np.ascontiguousarray(fneg_ref, np.float64),
    )


def _prior(rows, K: int):
    """A per-slot ``(m, K)`` log-prior on the solve grid as the kernel reads it, or ``None`` for none."""
    if rows is None:
        return None
    rows = np.ascontiguousarray(rows, np.float64)
    if rows.ndim != 2 or rows.shape[1] != K:
        raise ValueError(
            f"a ψ prior must be (m, K) on the solve grid; got {rows.shape} for K = {K}"
        )
    return rows


def psi_cube(
    u_pos,
    u_neg,
    allow_pos,
    allow_neg,
    fg_ref,
    fpos_ref,
    fneg_ref,
    *,
    kappa,
    od_g,
    od_r,
    lam,
    ambig: bool,
    gdna_logprior=None,
    lam_logprior=None,
    cube_rows=None,
    n_tilt: int | None = None,
):
    """ψ ITSELF over the ``(λ, θ)`` cube for ``m`` slots of one class — the strand term + the two Jeffreys
    arms (``_JEFFREYS_REF``) + the fitted gDNA prior + the λ-factor rows (+ the delivered cube rows) + the
    θ quadrature's log-weights — as ``(m, K, C)`` in float64, with the two strand-fraction grids it was
    evaluated on, ``(f_pos, f_neg)``, and the tilt ``tau`` they were built from. A single-strand call
    (``ambig=False``) has one column, the tilt of each slot's live strand (``τ = ±1``) and no weight. An
    AMBIG call places the θ nodes across each slot's strand term (``n_tilt`` of them, ``_TILT_NODES`` by
    derivation) for the MIXED hypothesis and appends the two PURE hypotheses as the columns ``τ = +1, −1``
    (the tilt atom: the three at equal reference weight, a delivered level on a strand ruling the other
    strand's atom out). The same code the solve runs (`native.psi_cube_native`); the gates read it here."""
    lam = np.ascontiguousarray(lam, np.float64)
    K = lam.shape[0]
    u_pos = np.ascontiguousarray(u_pos, np.float64)
    u_neg = np.ascontiguousarray(u_neg, np.float64)
    ap, an = np.ascontiguousarray(allow_pos, bool), np.ascontiguousarray(allow_neg, bool)
    fg_ref, fpos_ref, fneg_ref = _reference_composition(ap, an, fg_ref, fpos_ref, fneg_ref)
    n_tilt = int(_TILT_NODES if n_tilt is None else n_tilt)
    m, C = u_pos.shape[0], (n_tilt + 2 if ambig else 1)
    psi, f_pos, f_neg, tau = (np.zeros((m, K, C)) for _ in range(4))
    psi_cube_native(
        u_pos=u_pos,
        u_neg=u_neg,
        allow_pos=ap,
        allow_neg=an,
        fg_ref=fg_ref,
        fpos_ref=fpos_ref,
        fneg_ref=fneg_ref,
        kappa=float(kappa),
        od_g=float(od_g),
        od_r=float(od_r),
        lam=lam,
        gdna_logprior=_prior(gdna_logprior, K),
        lam_logprior=_prior(lam_logprior, K),
        **_cube_args(cube_rows, lam),
        n_tilt=n_tilt,
        ambig=bool(ambig),
        out_psi=psi,
        out_fpos=f_pos,
        out_fneg=f_neg,
        out_tau=tau,
    )
    return psi, f_pos, f_neg, tau


def posterior_median_fg(post, lam):
    """Per-slot point estimate of ``f_g``: the posterior's ½-QUANTILE, read off the CDF on the uniform λ
    lattice — a continuous quantile (the grid mass as a histogram with edges at the midpoints, the crossing
    bin interpolated ON λ, then mapped through σ: median equivariance, which is why ``f_g`` is a median
    and not a mean; `DESIGN.md` §6c). ``post``: ``(m, K)``; ``lam``: ``(K,)``. Returns ``(m,)``."""
    post = np.ascontiguousarray(post, np.float64)
    out = np.zeros(post.shape[0])
    psi_posterior_median(post, np.ascontiguousarray(lam, np.float64), out)
    return out


def compose(f_g, w_pos, allow_pos, allow_neg):
    """ψ's composition as the MAP from its two parameters — ``f_g`` (the level, a median) and ``w_pos``
    (the + strand's share of the RNA) — onto the simplex: ``f_pos = (1 − f_g)·w``, ``f_neg = (1 − f_g)·(1 −
    w)``, so closure is structural; the share is clamped to ``[0, 1]`` and restricted to the admissible
    strands (a single-strand slot's whole RNA sits on its live strand whatever the share says; a slot with
    neither strand has no RNA to place). Returns ``(f_pos, f_neg)``."""
    f_g = np.ascontiguousarray(f_g, np.float64)
    ap, an = np.ascontiguousarray(allow_pos, bool), np.ascontiguousarray(allow_neg, bool)
    f_pos, f_neg = np.zeros(f_g.shape[0]), np.zeros(f_g.shape[0])
    psi_compose(f_g, np.ascontiguousarray(w_pos, np.float64), ap, an, f_pos, f_neg)
    return f_pos, f_neg


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
    n_tilt: int | None = None,
) -> RegionDeconv:
    """THE per-slot solve for every slot of a block, dispatched: one native call (`native.psi_solve`) over
    the slots that admit a strand and carry a fragment, each on its own ``(λ, θ)`` cube — ψ built as
    :func:`psi_cube` describes and read out once: ``f_g`` the posterior median over the θ-marginal
    λ-posterior, ``f_pos`` / ``f_neg`` its image under :func:`compose` with the tilt share ``w_pos`` the
    RNA-mass-weighted posterior share, and ``Var(log f_g)`` a grid moment over the λ-marginal — the one
    precision the tool reads (the landscape prior's training weight). Every slot is solved on its own, so
    the read-out is chunk-exact (gate: ``test_sweep.test_the_psi_solve_is_chunk_exact_so_a_block_split_moves_no_number``).

    No Jacobian and no tilt term are written, and that is the point of the θ coordinate: the Berger–Bernardo
    tilt conditional ``(1−τ²)^{−½}`` is cancelled identically by ``|dτ/dθ| = cos θ``, and on the two-group λ
    axis the log-rate conversions cancel ``log σ'(λ)``; the only weights in ψ are the θ quadrature's own.
    Both arms are always written — a fitted ``logP`` where there is one, else the ``_JEFFREYS_REF``
    reference — because omitting one is not neutral (the module docstring). The strand mixture's variance
    is frozen at the reference composition (``fg_ref`` / ``fpos_ref`` / ``fneg_ref``, per slot, the
    incoming belief; the structural-neutral default at init), so the count sets precision and not
    composition. Zero-count slots report 0.

    All array inputs are full length ``m``; ``gdna_logprior`` and ``lam_logprior`` are ``(m, K)`` on the
    σ(λ) grid; ``cube_rows`` is the RNA level lanes' delivery at the AMBIG slots (:class:`CubeRows`, on this
    grid), evaluated at each slot's own θ nodes inside its ψ; ``None`` or an absent slot changes nothing. EMPTY
    slots — no per-strand count and no unspliced or spliced mass — are not solved: at genome scale most
    slots carry no fragments, and their zeros are the solve's own answer. ``n_tilt`` is the derived
    ``_TILT_NODES`` unless a gate asks for another count."""
    u_pos = np.ascontiguousarray(u_pos, np.float64)
    u_neg = np.ascontiguousarray(u_neg, np.float64)
    m = u_pos.shape[0]
    ap, an = np.ascontiguousarray(allow_pos, bool), np.ascontiguousarray(allow_neg, bool)
    fg_ref, fpos_ref, fneg_ref = _reference_composition(ap, an, fg_ref, fpos_ref, fneg_ref)
    signal = (
        u_pos + u_neg + np.asarray(mass_unspl, np.float64) + np.asarray(mass_spliced, np.float64)
    ) > 0.0
    slots = np.flatnonzero((ap | an) & signal).astype(np.int64)
    lam, _fg = _logodds_grid(int(n_grid), L)
    out = {k: np.zeros(m) for k in ("fg", "fp", "fn", "vg")}
    if slots.size:
        psi_solve(
            slots=slots,
            u_pos=u_pos,
            u_neg=u_neg,
            allow_pos=ap,
            allow_neg=an,
            fg_ref=fg_ref,
            fpos_ref=fpos_ref,
            fneg_ref=fneg_ref,
            kappa=float(kappa),
            od_g=float(od_g),
            od_r=float(od_r),
            lam=lam,
            gdna_logprior=_prior(gdna_logprior, lam.shape[0]),
            lam_logprior=_prior(lam_logprior, lam.shape[0]),
            **_cube_args(cube_rows, lam),
            n_tilt=int(_TILT_NODES if n_tilt is None else n_tilt),
            out_fg=out["fg"],
            out_fpos=out["fp"],
            out_fneg=out["fn"],
            out_var=out["vg"],
        )
    return RegionDeconv(
        gdna_frac=out["fg"],
        rna_pos_frac=out["fp"],
        rna_neg_frac=out["fn"],
        gdna_frac_var=out["vg"],
    )
