"""The per-slot composition solve ψ on the ``(λ, θ)`` grid, one solve for every slot class; layer 3.

The magnitude coordinate is the gDNA-vs-RNA log-odds ``λ = logit(f_g) = log ρ_g − log ρ_rna``, gridded on a fixed
``[−L, L]`` window and read out as ``f_g = σ(λ)``: log-odds spans the decades of ρ_g and resolves both vertices,
which a linear lattice cannot. The integrand is ``strand + gDNA arm + RNA arm + the λ-factor rows``. Each arm carries
the Jeffreys reference ``+½·log f`` on its own fraction (the kernel's ``JEFFREYS_REF``), and the gDNA arm adds the
fitted landscape's log-rate prior when there is one. The reference is derived two ways: Jeffreys for a Poisson rate
(``I(ρ) ∝ 1/ρ`` ⇒ ``p(ρ) ∝ ρ^(−½)``, which is ``+½·log f`` as a log-rate density), and the Berger–Bernardo reference
with ``f_g`` of interest, whose marginal is Beta(½,½). It must always be written: a bare uniform-λ grid is Beta(0,0),
improper at both vertices. On an unstranded library the strand term is flat and the posterior is this reference.

The facts that shape the solve:

1. The composition is a two-group split, gDNA against total RNA, with the per-strand tilt ``τ`` a nuisance. Each
   ``logP`` is a density in log-rate, so its conversion ``−log ρ`` cancels ``log σ'(λ)`` exactly; no Jacobian is
   written.
2. The tilt is gridded as ``θ = arcsin(τ)``. The reference's tilt conditional is ``(1−τ²)^{−½}`` (``f_g`` and ``τ``
   are information-orthogonal), and ``|dτ/dθ| = cos θ`` cancels it, so the reference is
   ``½·log f_g + ½·log(1−f_g)`` for every slot class, with no endpoint singularity.
3. At fixed λ the strand term is Gaussian in τ with a θ peak narrowing as ``n^{−½}``, so ψ places its θ nodes across
   each ``(slot, λ)``'s own peak and weights them by the trapezoid rule (``tilt_window`` in ``psi_kernel.h``). The
   marginal is exact at every depth with the derived node count ``_TILT_NODES``.
4. A both-strand slot's tilt hypotheses are {pure +, pure −, mixed} at equal reference weight (the tilt atom): ψ
   carries two columns at ``τ = ±1`` beside the continuum, whose weights carry ``−log π``. A held RNA level on one
   strand (a :class:`CubeRows` profile) witnesses that strand's RNA and rules out the other strand's atom.

ψ never reads the spliced mass: RNA that splices at a sj jumps it rather than crossing it, so the unspliced mass is
gDNA plus RNA unspliced there, disjoint from the spliced channel. The solver is native (``native.psi_solve``,
``native/psi_kernel.h``), one slot at a time on its own cube: a single-strand slot is the one-column case, and a
both-strand slot marginalises the tilt over its θ nodes and the two atoms. Structurally RNA-free slots never reach
it; ``sweep.solve_chain`` gates them out.

Gate: ``tests/calibration/test_vertex_reference.py``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import expit

from ..config import CONSTANTS
from ..native import psi_solve
from .region_chain import RegionDeconv

# Public surface consumed by sweep / messages / region_geometry.
__all__ = [
    "CubeRows",
    "_logodds_grid",
    "_solve_regions_logodds_all",
]


def _block_rows(cells_per_row: int, itemsize: int) -> int:
    return max(
        1, CONSTANTS.resources.solve_block_bytes // max(1, int(cells_per_row) * int(itemsize))
    )


# The window L bounds f_g to [σ(−L), σ(L)], [4.5e-5, 1−4.5e-5] at L = 10: the widest f_g the grid can represent.
# Without a fitted prior ψ is proper (both arms carry JEFFREYS_REF) and the answer does not depend on L; an improper ψ
# would let L set the prior strength. With a fitted landscape the low end can bind: ψ reads the curve at
# log rho = log f_c + log M − log E, and on a gDNA-poor library σ(−L) can sit above the density the prior points at.
# `landscape.required_logodds_window` derives the window that demand needs.
def _logodds_grid(n_grid: int, L: float):
    """The fixed log-odds lattice: ``λ`` uniform on ``[−L, L]`` (``K = n_grid`` points, ascending) and
    the matching ``f_g = σ(λ)`` (also ascending). Returns ``(lam, fg)``, each length ``K``."""
    lam = np.linspace(-float(L), float(L), int(n_grid))
    return lam, expit(lam)


@dataclass(frozen=True)
class CubeRows:
    """The RNA level lanes' delivery at a block's both-strand slots, one row per delivered slot.

    Each row holds, per strand, the held level profile over ``u = log(ρ/ρ_ref)`` on the solve grid with its presence
    bit, plus the slot's total, its RNA opportunity and the lanes' reference density (shared by both strands); ``u``
    is the one grid every row is on. ψ reads a row at its own θ nodes: at each ``(λ, θ)`` the strand share
    ``f_s = (1 − f_g)(1 ± τ)/2`` implies the density ``f_s·n/a_r``, at which the profile is read, so nothing is
    interpolated in θ. The map is monotone in each share, so a one-sided profile arrives as a wall in the cube. Built
    by the policy's solve (``native/transfer_kernel.h``) and read by ``native.psi_solve``; the arrays are the kernels'
    arguments.
    """

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

    @classmethod
    def concat(cls, parts: list) -> CubeRows:
        return cls(
            **{f: np.concatenate([getattr(q, f) for q in parts]) for f in cls.PER_ROW}, u=parts[0].u
        )

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


# The θ window's truncation, in nats: the strand term's mass outside the window is below double precision, since
# ``erfc(√T) ≈ e^{−T}/√(πT) < ε₆₄`` at ``T = −log ε₆₄``. The kernel carries the same constant (``psi_kernel.h``).
_T_NATS = -np.log(np.finfo(np.float64).eps)

# The node count that resolves the peak inside its window: the interior window is ``2√(2T)·σ_θ`` wide and the
# trapezoid rule's error on a Gaussian at spacing ``h`` is ``2·e^{−2π²(σ_θ/h)²}``, below ``e^{−T}`` once
# ``h ≤ σ_θ·π√2/√T``, so ``K_t − 1 ≥ 2T/π``.
_TILT_NODES = int(np.ceil(2.0 * _T_NATS / np.pi)) + 1


def _cube_args(cube_rows, u) -> dict:
    """The delivered rows as ψ's kernel takes them — the table's arrays, or an empty table's on this grid."""
    if cube_rows is None:
        return CubeRows.blank(0, u).kernel_args()
    if cube_rows.u.shape != u.shape or not np.array_equal(cube_rows.u, u):
        raise ValueError("the delivered rows are not on the solve grid")
    return cube_rows.kernel_args()


def _reference_composition(allow_pos, allow_neg, fg_ref, fpos_ref, fneg_ref):
    """The per-slot composition the strand variance is frozen at: the incoming belief where the sweep supplies one,
    else ``f_g = ½`` with the other ½ split among the live strands. It sets the precision, not the location."""
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
    """A per-slot ``(m, K)`` log-prior row array on the solve grid as the kernel reads it, or ``None``."""
    if rows is None:
        return None
    rows = np.ascontiguousarray(rows, np.float64)
    if rows.ndim != 2 or rows.shape[1] != K:
        raise ValueError(
            f"a ψ prior must be (m, K) on the solve grid; got {rows.shape} for K = {K}"
        )
    return rows


def _arm(gdna_prior, gdna_support, m: int):
    """The fitted gDNA arm as the kernel takes it, ``(log_rho, logP, mass, eff)``, or ``None`` for the reference alone.

    ``(log_rho, logP)`` is the landscape's curve, a density in log-rate, and ``(mass, eff)`` the per-slot support it
    is read on; ψ evaluates it at each cell's ``log ρ_c = log f_g + log M − log E``, ends held constant. Both arrive as
    arrays, never as the landscape, which sits in a higher layer.
    """
    if gdna_prior is None:
        if gdna_support is not None:
            raise ValueError(
                "a gDNA support without a curve: the arm is the curve read on the support"
            )
        return None
    if gdna_support is None:
        raise ValueError("the fitted gDNA arm needs its per-slot support (mass, eff) to be read on")
    log_rho, logP = (np.ascontiguousarray(a, np.float64) for a in gdna_prior)
    mass, eff = (np.ascontiguousarray(a, np.float64) for a in gdna_support)
    if log_rho.ndim != 1 or log_rho.shape != logP.shape or log_rho.shape[0] < 1:
        raise ValueError("the gDNA arm's curve must be two (G,) arrays")
    if mass.shape != (m,) or eff.shape != (m,):
        raise ValueError(
            f"the gDNA arm's support must be (mass, eff) of shape ({m},); got {mass.shape} and {eff.shape}"
        )
    return log_rho, logP, mass, eff


def _solve_regions_logodds_all(
    u_pos,
    u_neg,
    allow_pos,
    allow_neg,
    *,
    kappa,
    n_grid,
    L: float,
    gdna_prior=None,
    gdna_support=None,
    row_logprior=None,
    fg_ref=None,
    fpos_ref=None,
    fneg_ref=None,
    cube_rows=None,
    n_tilt: int | None = None,
    n_threads: int = 1,
) -> RegionDeconv:
    """Solve ψ for every slot of a block that admits a strand and carries a fragment, in one native call.

    Each slot is solved on its own ``(λ, θ)`` cube (``native/psi_kernel.h``) and read out as ``f_g``, the posterior
    median of the θ-marginal λ posterior; ``f_pos``/``f_neg``, its image on the admissible strands; and
    ``Var(log f_g)``, a grid moment over the λ marginal. Slots with no count report 0. Because every slot is solved on
    its own, the result is the same for any block split and any thread count.

    All arrays are length ``m``. ``gdna_prior`` is the landscape's curve ``(log_rho, logP)`` and ``gdna_support`` the
    per-slot ``(mass, eff)`` it is read on (:func:`_arm`). ``row_logprior`` is the policy's delivered composition
    rows, ``(m, K)`` on the λ grid; ``cube_rows`` is the RNA level lanes' delivery (:class:`CubeRows`). The strand
    variance is frozen at ``fg_ref``/``fpos_ref``/``fneg_ref`` (:func:`_reference_composition`), so the count sets
    precision, not composition. ``n_tilt`` defaults to ``_TILT_NODES``; ``n_threads`` is the kernel's thread budget
    (0 is every core).
    """
    u_pos = np.ascontiguousarray(u_pos, np.float64)
    u_neg = np.ascontiguousarray(u_neg, np.float64)
    m = u_pos.shape[0]
    ap, an = np.ascontiguousarray(allow_pos, bool), np.ascontiguousarray(allow_neg, bool)
    fg_ref, fpos_ref, fneg_ref = _reference_composition(ap, an, fg_ref, fpos_ref, fneg_ref)
    signal = (u_pos + u_neg) > 0.0
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
            lam=lam,
            gdna=_arm(gdna_prior, gdna_support, m),
            row_logprior=_prior(row_logprior, lam.shape[0]),
            **_cube_args(cube_rows, lam),
            n_tilt=int(_TILT_NODES if n_tilt is None else n_tilt),
            out_fg=out["fg"],
            out_fpos=out["fp"],
            out_fneg=out["fn"],
            out_var=out["vg"],
            n_threads=int(n_threads),
        )
    return RegionDeconv(
        gdna_frac=out["fg"],
        rna_pos_frac=out["fp"],
        rna_neg_frac=out["fn"],
        gdna_frac_var=out["vg"],
    )
