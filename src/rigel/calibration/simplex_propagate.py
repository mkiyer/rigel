"""Forward–backward gDNA-density propagation — the rev-2 simplex solver (increment 3).

Spec: ``docs/calibration/propagation_message_passing.md`` (§6 dry-run). The deconvolution is inference on
a **chain** of nodes (a locus is linear → a tree → belief propagation is exact in two sweeps,
order-independent). The design review (issue A) established the key factoring:

* **Only gDNA propagates.** gDNA contamination density ``ρ_g`` is piecewise-uniform by *enrichment class*
  (on-target exon vs off-target intron/intergenic — capture enriches exons), so it is the **transferable**
  quantity. RNA does not propagate (its abundance is local, varying >10000× between neighbours). So the
  cross-node coupling is a 1-D smoothness prior on ``ρ_g`` **within an enrichment class**, and the message
  is scalar — a Gaussian on ``ρ_g``. The non-Gaussian, simplex-constrained work (the strand mixture, the
  one-sided spliced bound, the gDNA prior) is all **local** and handled by the per-node grid solve
  (:func:`simplex.solve_node`).
* **State is fractional** (the pie ``(f_rna₊, f_rna₋, f_g)``): a propagated density can only *redistribute*
  a node's observed mass, never inject — the capping / no-injection safety.

Pipeline: (1) each node's **local** ``ρ_g`` self-estimate + precision (intergenic seed ``f_g=1``;
single-strand → strand posterior; AMBIG → none); (2) a **Kalman/RTS smoother** (forward then backward —
the exact two-sweep BP for a scalar chain) over per-``(reference, enrichment-class)`` chains fills every
node's ``ρ_g`` from the seeds, precision-weighted; (3) each node solves its pie via
:func:`simplex.solve_node` with the smoothed ``ρ_g`` as the **count evidence** (``f_g ≈ ρ_g·L/U``) plus its
own local strand/prior. Boundary sides inherit their region's smoothed ``ρ_g`` (gDNA density is shared
across a locus). Returns ``(regions, left, right)`` :class:`NodeDeconv`, the same interface as the
opt-in count-cascade ``propagation.py`` (this is built standalone; wiring into ``calibrate`` is a later
increment).

**Placeholders pending derivation (no magic widths):** the RTS process variance ``Q`` (the per-hop
coupling variance) defaults to the median local observation variance — the count ``var~mean`` model is the
proper derivation (increment 4); ``gdna_prior_count`` (the weak gDNA prior) likewise. The forward–backward
*mechanics* validated here do not depend on their exact values.
"""

from __future__ import annotations

import numpy as np

from .run_fill import same_ref_left_right
from .signature import RegionType, TS_AMBIG, TS_NEG, TS_NONE, TS_POS, coarse_type_array
from .simplex import init_from_signature, solve_node
from .strand_deconv import NodeDeconv, strand_posterior_gdna_frac

__all__ = ["propagate_simplex"]

_EPS = 1.0e-9
_DIFFUSE_VAR = 1.0e12  # diffuse prior variance for an unobserved chain start


def _local_gdna_density(substrate, region_arrays, gdna_eff_len, kappa, od_g, od_r, n_grid):
    """Per-region **local** gDNA-density estimate ``ρ_g`` + precision (the chain's observations).

    Seeds emit a density; the rest emit nothing (precision 0) and are filled by the smoother.
    * intergenic (``TS_NONE``): ``f_g=1`` ⇒ ``ρ_g = U/L``, precision ``L²/U`` (Poisson on ``U``).
    * single-strand (``TS_POS/NEG``): the strand posterior gives ``f_g`` (+ variance) ⇒ ``ρ_g = f_g·U/L``,
      precision ``1/Var(ρ_g)`` — automatically small at κ≈½ / thin (a flat posterior), so an
      uninformative single-strand node is not a confident seed.
    * AMBIG: no sense split ⇒ no local density (precision 0).
    """
    ts = np.asarray(region_arrays.strand_class)
    c = substrate.contained
    u_pos = c.n_unspliced_pos.astype(np.float64)
    u_neg = c.n_unspliced_neg.astype(np.float64)
    U = u_pos + u_neg
    L = np.maximum(np.asarray(gdna_eff_len, dtype=np.float64), _EPS)
    rho = np.full(U.shape, np.nan, dtype=np.float64)
    prec = np.zeros(U.shape, dtype=np.float64)

    inter = (ts == TS_NONE) & (U > 0.0)
    rho[inter] = U[inter] / L[inter]
    prec[inter] = L[inter] ** 2 / np.maximum(U[inter], _EPS)

    for strand in (TS_POS, TS_NEG):
        idx = np.flatnonzero((ts == strand) & (U > 0.0))
        if idx.size:
            sense = (u_neg if strand == TS_NEG else u_pos)[idx]
            anti = U[idx] - sense
            g_q, var_g = strand_posterior_gdna_frac(
                sense, anti, kappa,
                gdna_strand_overdispersion=od_g, rna_strand_overdispersion=od_r, n_grid=n_grid,
            )
            rho[idx] = g_q * U[idx] / L[idx]
            var_rho = (U[idx] / L[idx]) ** 2 * np.maximum(var_g, _EPS)
            prec[idx] = 1.0 / np.maximum(var_rho, _EPS)
    return rho, prec, U, L


def _rts_smooth(y: np.ndarray, r: np.ndarray, process_var: float) -> tuple[np.ndarray, np.ndarray]:
    """Scalar Kalman filter + RTS smoother over one ordered chain — the exact two-sweep BP.

    State = ``ρ_g`` under a random-walk smoothness prior (process variance ``process_var`` per hop);
    observations ``y`` with precision ``r`` (``r=0`` ⇒ no observation at that node). Returns the smoothed
    mean + variance per node. The forward (filter) and backward (smoother) passes are the two sweeps; the
    result is independent of which end starts (a tree has a unique fixed point).
    """
    k = y.shape[0]
    m_pred = np.empty(k)
    p_pred = np.empty(k)
    m_filt = np.empty(k)
    p_filt = np.empty(k)
    m, p = 0.0, _DIFFUSE_VAR
    for i in range(k):
        mp, pp = m, p + process_var  # predict (random walk)
        m_pred[i], p_pred[i] = mp, pp
        if r[i] > 0.0:
            p = 1.0 / (1.0 / pp + r[i])
            m = p * (mp / pp + r[i] * y[i])
        else:
            m, p = mp, pp
        m_filt[i], p_filt[i] = m, p
    ms = m_filt.copy()
    ps = p_filt.copy()
    for i in range(k - 2, -1, -1):
        c = p_filt[i] / max(p_pred[i + 1], _EPS)
        ms[i] = m_filt[i] + c * (ms[i + 1] - m_pred[i + 1])
        ps[i] = p_filt[i] + c * c * (ps[i + 1] - p_pred[i + 1])
    return ms, ps


def _smooth_density(rho_obs, rho_prec, ref_id, on_target, process_var):
    """Run the RTS smoother over each ``(reference, enrichment-class)`` chain → ``(ρ_g, precision)``.

    Enrichment-class chains (on-target exon vs off-target intron/intergenic) keep ``ρ_g`` comparable
    within a chain without a capture factor (issue A): exon densities inform exon nodes, off-target
    densities inform off-target nodes. Region indices within a ``(ref, class)`` group are already in
    genomic order, so the chain hop is the genomic neighbour.
    """
    r = rho_obs.shape[0]
    rho_post = np.zeros(r, dtype=np.float64)
    prec_post = np.zeros(r, dtype=np.float64)
    ref = np.asarray(ref_id)
    cls = on_target.astype(np.int64)
    keys = ref.astype(np.int64) * 2 + cls
    for key in np.unique(keys):
        idx = np.flatnonzero(keys == key)  # ascending ⇒ genomic order
        y = np.where(np.isnan(rho_obs[idx]), 0.0, rho_obs[idx])
        ms, ps = _rts_smooth(y, rho_prec[idx], process_var)
        rho_post[idx] = np.clip(ms, 0.0, None)
        prec_post[idx] = 1.0 / np.maximum(ps, _EPS)
    return rho_post, prec_post


def _solve_view(view, allow_pos, allow_neg, rho_post, prec_post, eff_len, kappa, od_g, od_r,
                gdna_prior_count, n_grid):
    """Solve one view's pie (region-contained or a boundary side) given the smoothed ρ_g as count evidence."""
    u_pos = view.n_unspliced_pos.astype(np.float64)
    u_neg = view.n_unspliced_neg.astype(np.float64)
    U = u_pos + u_neg
    L = np.maximum(np.asarray(eff_len, dtype=np.float64), _EPS)
    mass_unspl = np.asarray(view.mass_unspliced, dtype=np.float64)
    mass_spliced = np.asarray(view.mass_spliced, dtype=np.float64)
    # smoothed density → count evidence on f_g: f_g ≈ ρ_g·L/U, var(f_g) = (L/U)²/prec_post.
    with np.errstate(divide="ignore", invalid="ignore"):
        count_frac = np.where(U > 0.0, np.clip(rho_post * L / np.maximum(U, _EPS), 0.0, 1.0), 0.0)
        count_prec = np.where(U > 0.0, prec_post * (U / L) ** 2, 0.0)
    sol = solve_node(
        u_pos, u_neg, kappa=kappa, count_gdna_frac=count_frac, count_precision=count_prec,
        allow_pos=allow_pos, allow_neg=allow_neg, strand_od_gdna=od_g, strand_od_rna=od_r,
        gdna_prior_count=gdna_prior_count, n_grid=n_grid,
    )
    f_g = np.where(U > 0.0, sol.f_g, 0.0)
    return NodeDeconv(
        gdna_mass=f_g * mass_unspl, rna_mass=(1.0 - f_g) * mass_unspl + mass_spliced,
        gdna_frac=f_g, gdna_frac_var=sol.f_g_var,
    )


def _side_allow(ts_self, ts_other):
    """A boundary side's active RNA strands = the union of its two flanks' signatures.

    NONE↔NONE → neither (gDNA only); POS↔NONE → +; POS↔NEG → both (AMBIG side). Mirrors the
    ``_side_strand_orientation`` observability without re-deriving the sense split (the smoothed density
    is the gDNA evidence; the side's own strand resolves the rest in ``solve_node``).
    """
    def has_pos(ts):
        return (ts == TS_POS) | (ts == TS_AMBIG)

    def has_neg(ts):
        return (ts == TS_NEG) | (ts == TS_AMBIG)

    return has_pos(ts_self) | has_pos(ts_other), has_neg(ts_self) | has_neg(ts_other)


def propagate_simplex(
    substrate, region_arrays, *, gdna_region_eff_len, gdna_boundary_side_eff_len, rna_sense_frac,
    gdna_strand_overdispersion=0.0, rna_strand_overdispersion=0.0, gdna_prior_count=0.5,
    process_var=None, n_grid=60,
):
    """Forward–backward gDNA-density propagation → ``(regions, left, right)`` :class:`NodeDeconv`.

    See the module docstring. ``process_var`` (the RTS per-hop coupling variance) defaults to the median
    local observation variance — a flagged placeholder for the count ``var~mean`` model (increment 4).
    """
    ts = np.asarray(region_arrays.strand_class)
    ref_id = np.asarray(region_arrays.ref_id)
    sig = np.asarray(region_arrays.signature)
    kappa = float(rna_sense_frac)
    od_g, od_r = gdna_strand_overdispersion, rna_strand_overdispersion

    rho_obs, rho_prec, _U, _L = _local_gdna_density(
        substrate, region_arrays, gdna_region_eff_len, kappa, od_g, od_r, n_grid
    )
    if process_var is None:
        seen = rho_prec > 0.0
        process_var = float(np.median(1.0 / rho_prec[seen])) if seen.any() else 1.0
    on_target = coarse_type_array(sig) == int(RegionType.EXON)
    rho_post, prec_post = _smooth_density(rho_obs, rho_prec, ref_id, on_target, process_var)

    init = init_from_signature(ts)
    regions = _solve_view(
        substrate.contained, init.allow_pos, init.allow_neg, rho_post, prec_post,
        gdna_region_eff_len, kappa, od_g, od_r, gdna_prior_count, n_grid,
    )

    # Boundary sides inherit their region's smoothed ρ_g (shared gDNA density across the locus); the side's
    # active RNA strands = the union of its two flanks (the gene-edge / opposite-strand rule).
    _l_same, _r_same = same_ref_left_right(ref_id)
    r = ts.shape[0]
    ts_prev = np.empty(r, dtype=ts.dtype)
    ts_next = np.empty(r, dtype=ts.dtype)
    ts_prev[0] = TS_NONE
    ts_next[-1] = TS_NONE
    if r > 1:
        ts_prev[1:] = np.where(_l_same[1:], ts[:-1], TS_NONE)
        ts_next[:-1] = np.where(_r_same[:-1], ts[1:], TS_NONE)
    side_eff = np.asarray(gdna_boundary_side_eff_len, dtype=np.float64)

    lp, ln = _side_allow(ts, ts_prev)
    left = _solve_view(substrate.left, lp, ln, rho_post, prec_post, side_eff, kappa, od_g, od_r,
                       gdna_prior_count, n_grid)
    rp, rn = _side_allow(ts, ts_next)
    right = _solve_view(substrate.right, rp, rn, rho_post, prec_post, side_eff, kappa, od_g, od_r,
                        gdna_prior_count, n_grid)
    return regions, left, right
