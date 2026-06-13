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

from .density_model import count_observable_masks, density_variance_curve, node_gdna_density
from .run_fill import same_ref_left_right
from .signature import RegionType, TS_AMBIG, TS_NEG, TS_NONE, TS_POS, coarse_type_array
from .simplex import init_from_signature, solve_node
from .strand_deconv import NodeDeconv, strand_posterior_gdna_frac

__all__ = ["deconv_regions_simplex", "propagate_simplex"]

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


def _rts_smooth(y: np.ndarray, r: np.ndarray, process_var) -> tuple[np.ndarray, np.ndarray]:
    """Scalar Kalman filter + RTS smoother over one ordered chain — the exact two-sweep BP.

    State = ``ρ_g`` under a random-walk smoothness prior; observations ``y`` with precision ``r``
    (``r=0`` ⇒ no observation at that node). ``process_var`` is the per-hop process variance ``Q`` —
    a **scalar** (constant coupling, increment 3) or a **per-node array** (``Q[i]`` = noise on the hop
    into node ``i``; the increment-4 ``var~mean`` coupling). Returns the smoothed mean + variance per
    node. The forward (filter) and backward (smoother) passes are the two sweeps; the result is
    independent of which end starts (a tree has a unique fixed point).
    """
    k = y.shape[0]
    q = np.broadcast_to(np.asarray(process_var, dtype=np.float64), (k,))
    m_pred = np.empty(k)
    p_pred = np.empty(k)
    m_filt = np.empty(k)
    p_filt = np.empty(k)
    m, p = 0.0, _DIFFUSE_VAR
    for i in range(k):
        mp, pp = m, p + q[i]  # predict (random walk; per-hop process variance)
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
    genomic order, so the chain hop is the genomic neighbour. ``process_var`` may be a scalar or a
    per-node array (sliced per chain).
    """
    r = rho_obs.shape[0]
    rho_post = np.zeros(r, dtype=np.float64)
    prec_post = np.zeros(r, dtype=np.float64)
    ref = np.asarray(ref_id)
    cls = on_target.astype(np.int64)
    keys = ref.astype(np.int64) * 2 + cls
    q_arr = np.broadcast_to(np.asarray(process_var, dtype=np.float64), (r,))
    for key in np.unique(keys):
        idx = np.flatnonzero(keys == key)  # ascending ⇒ genomic order
        y = np.where(np.isnan(rho_obs[idx]), 0.0, rho_obs[idx])
        ms, ps = _rts_smooth(y, rho_prec[idx], q_arr[idx])
        rho_post[idx] = np.clip(ms, 0.0, None)
        prec_post[idx] = 1.0 / np.maximum(ps, _EPS)
    return rho_post, prec_post


def _solve_view(view, allow_pos, allow_neg, rho_post, beta, eff_len, kappa, od_g, od_r,
                gdna_prior_count, n_grid):
    """Solve one view's pie: the node's own strand + the smoothed count density at trust ``β``.

    PHASE-1 count trust: ``solve_node``'s count term gets the **fixed effective precision ``β``** (the
    successor to the old hard-capped ``I₀``), not the count's own variance — so the strand
    (curvature ``N·(2κ−1)²``, vanishing at κ=½) governs single-strand nodes where it is informative and
    the count cleanly takes over where the strand is silent (κ=½ / AMBIG). See count_trust_design.md.
    """
    u_pos = view.n_unspliced_pos.astype(np.float64)
    u_neg = view.n_unspliced_neg.astype(np.float64)
    U = u_pos + u_neg
    L = np.maximum(np.asarray(eff_len, dtype=np.float64), _EPS)
    mass_unspl = np.asarray(view.mass_unspliced, dtype=np.float64)
    mass_spliced = np.asarray(view.mass_spliced, dtype=np.float64)
    with np.errstate(divide="ignore", invalid="ignore"):
        count_frac = np.where(U > 0.0, np.clip(rho_post * L / np.maximum(U, _EPS), 0.0, 1.0), 0.0)
    count_prec = np.where(U > 0.0, float(beta), 0.0)  # fixed count trust β (phase-1)
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


def deconv_regions_simplex(
    substrate, region_arrays, count_gdna_frac, *, rna_sense_frac, gdna_strand_overdispersion=0.0,
    rna_strand_overdispersion=0.0, count_trust_beta=10.0, gdna_prior_count=0.0, n_grid=60,
):
    """PHASE-1 per-region simplex deconvolution — the per-node β combine, **no spatial propagation**.

    Each region's pie is solved by :func:`simplex.solve_node` combining its **own strand** (``u_pos/u_neg``;
    intrinsic precision ``N·(2κ−1)²``, vanishing at κ=½) with the **count signal** ``count_gdna_frac`` (the
    splice-upgraded, count-mean-bias-corrected gDNA fraction) at the **fixed trust ``β=count_trust_beta``**.
    The strand governs single-strand nodes where informative; the count governs at κ=½ / AMBIG. This is the
    simplex/pie successor to the old ``w=I/(I+I₀)`` combine (count precision hard-capped at ``I₀``) — it does
    **not** spatially smooth (the count-density RTS smeared unrelated exons; strand→AMBIG propagation is a
    later phase). See docs/calibration/count_trust_design.md.
    """
    ts = np.asarray(region_arrays.strand_class)
    c = substrate.contained
    u_pos = c.n_unspliced_pos.astype(np.float64)
    u_neg = c.n_unspliced_neg.astype(np.float64)
    U = u_pos + u_neg
    mass_unspl = np.asarray(c.mass_unspliced, dtype=np.float64)
    mass_spliced = np.asarray(c.mass_spliced, dtype=np.float64)
    cf = np.clip(np.asarray(count_gdna_frac, dtype=np.float64), 0.0, 1.0)

    # Default = the count signal (the count governs where the strand is undefined: AMBIG, intergenic, κ=½).
    # A region is STRAND-OBSERVABLE iff its transcript strand is defined (POS/NEG); for those the simplex
    # solve combines the OWN strand (precision N·(2κ−1)², vanishing at κ=½) with the count at trust β —
    # the strand governs where informative. AMBIG (overlapping opposite strands → no valid sense) and
    # intergenic stay count-only, exactly as the production combine (w=0 there).
    f_g = np.where(U > 0.0, cf, 0.0)
    f_g_var = np.zeros_like(f_g)
    obs = np.flatnonzero(((ts == TS_POS) | (ts == TS_NEG)) & (U > 0.0))
    if obs.size:
        allow_pos = ts[obs] == TS_POS
        allow_neg = ts[obs] == TS_NEG
        sol = solve_node(
            u_pos[obs], u_neg[obs], kappa=float(rna_sense_frac), count_gdna_frac=cf[obs],
            count_precision=float(count_trust_beta), allow_pos=allow_pos, allow_neg=allow_neg,
            strand_od_gdna=gdna_strand_overdispersion, strand_od_rna=rna_strand_overdispersion,
            gdna_prior_count=gdna_prior_count, n_grid=n_grid,
        )
        f_g[obs] = sol.f_g
        f_g_var[obs] = sol.f_g_var
    return NodeDeconv(
        gdna_mass=f_g * mass_unspl, rna_mass=(1.0 - f_g) * mass_unspl + mass_spliced,
        gdna_frac=f_g, gdna_frac_var=f_g_var,
    )


def _coupling_process_var(substrate, region_arrays, gdna_region_eff_len, gdna_fl_mean):
    """Per-node RTS process variance ``Q`` from the triplet ``var~mean`` curve (increment 4).

    The chain coupling is gDNA-density continuity; its per-hop variance is how much ``ρ_g`` varies
    node-to-node, which the **triplet disagreement** (``{left boundary, count-observable contained, right
    boundary}``) estimates directly (the density gradient across a region's span). We fit
    :func:`density_variance_curve` on the triplet observations and read ``σ²_density`` at each node's local
    density ``μ`` (the count module's imputed density — available *before* propagation, so non-circular).
    Returns ``(Q, fallback_scalar)``: ``Q`` is ``NaN`` where the curve is undefined (the caller fills those
    with ``fallback_scalar``, the median local observation variance — the increment-3 default).
    """
    sig = np.asarray(region_arrays.signature)
    ref_id = np.asarray(region_arrays.ref_id)
    r = sig.shape[0]
    L = np.maximum(np.asarray(gdna_region_eff_len, dtype=np.float64), _EPS)
    inv_fl = 1.0 / gdna_fl_mean if gdna_fl_mean > 0.0 else 0.0

    def total(view):
        return view.n_unspliced_pos.astype(np.float64) + view.n_unspliced_neg.astype(np.float64)

    region_obs, boundary_obs = count_observable_masks(sig, ref_id)
    left_same, right_same = same_ref_left_right(ref_id)
    left_anchor = np.zeros(r, dtype=bool)
    right_anchor = np.zeros(r, dtype=bool)
    if r > 1:
        left_anchor[1:] = boundary_obs[:-1] & left_same[1:]
        right_anchor[:-1] = boundary_obs[:-1] & right_same[:-1]
    d_left = np.where(left_anchor, total(substrate.left) * inv_fl, np.nan)
    d_right = np.where(right_anchor, total(substrate.right) * inv_fl, np.nan)
    contained_density = np.where(region_obs & (L > _EPS), total(substrate.contained) / L, np.nan)

    # μ for the curve lookup = the count module's local density estimate (pre-propagation, non-circular).
    nd = node_gdna_density(substrate, region_arrays, L, gdna_fl_mean, need_count_variance=False)
    sigma_d2 = density_variance_curve(
        np.asarray(nd.density, dtype=np.float64), d_left=d_left, d_right=d_right,
        left_ok=left_anchor, right_ok=right_anchor,
        contained=contained_density, contained_ok=(region_obs & (L > _EPS)),
    )
    return sigma_d2


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
    substrate, region_arrays, *, count_gdna_frac, gdna_region_eff_len, gdna_boundary_side_eff_len,
    gdna_fl_mean, rna_sense_frac, gdna_strand_overdispersion=0.0, rna_strand_overdispersion=0.0,
    gdna_prior_count=0.5, count_trust_beta=10.0, process_var=None, n_grid=60,
):
    """Per-node simplex deconvolution with the count signal at trust ``β`` → ``(regions, left, right)``.

    PHASE-1 (count_trust_design.md). The **count signal** is the splice-upgraded ``count_gdna_frac`` (the
    count-mean-bias-corrected gDNA fraction), converted to a density and **smoothed along
    (reference, enrichment-class) chains** by the RTS (the propagation of the count magnitude). Each node's
    pie is then solved by :func:`simplex.solve_node`, combining the node's **own strand** (intrinsic
    precision ``N·(2κ−1)²``) with that smoothed count at the **fixed trust ``β = count_trust_beta``** — so
    the strand governs where informative and the count takes over where the strand is silent (κ=½ / AMBIG).

    (Phase 2 will propagate single-strand neighbours' *strand-derived* densities into AMBIG nodes; phase
    3–4 make ``β`` per-node / derived. The RTS ``Q`` is the triplet ``var~mean`` coupling variance;
    ``process_var`` overrides it for tests/ablation.)
    """
    ref_id = np.asarray(region_arrays.ref_id)
    sig = np.asarray(region_arrays.signature)
    ts = np.asarray(region_arrays.strand_class)
    kappa = float(rna_sense_frac)
    od_g, od_r = gdna_strand_overdispersion, rna_strand_overdispersion
    count_gdna_frac = np.asarray(count_gdna_frac, dtype=np.float64)

    # The count signal as a density (the propagated magnitude), + the var~mean coupling/observation variance.
    c = substrate.contained
    mass = np.asarray(c.mass_unspliced, dtype=np.float64)
    u = c.n_unspliced_pos.astype(np.float64) + c.n_unspliced_neg.astype(np.float64)
    L = np.maximum(np.asarray(gdna_region_eff_len, dtype=np.float64), _EPS)
    sigma_d2 = _coupling_process_var(substrate, region_arrays, gdna_region_eff_len, gdna_fl_mean)
    fin = np.isfinite(sigma_d2) & (sigma_d2 > 0.0)
    fallback_var = float(np.median(sigma_d2[fin])) if fin.any() else 1.0
    rho_count = np.where(u > 0.0, np.clip(count_gdna_frac, 0.0, 1.0) * mass / L, np.nan)
    prec_obs = np.where(u > 0.0, 1.0 / np.where(fin, sigma_d2, fallback_var), 0.0)
    q = process_var if process_var is not None else np.where(fin, sigma_d2, fallback_var)
    on_target = coarse_type_array(sig) == int(RegionType.EXON)
    rho_post, _prec_post = _smooth_density(rho_count, prec_obs, ref_id, on_target, q)

    init = init_from_signature(ts)
    regions = _solve_view(
        c, init.allow_pos, init.allow_neg, rho_post, count_trust_beta,
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
    left = _solve_view(substrate.left, lp, ln, rho_post, count_trust_beta, side_eff, kappa, od_g, od_r,
                       gdna_prior_count, n_grid)
    rp, rn = _side_allow(ts, ts_next)
    right = _solve_view(substrate.right, rp, rn, rho_post, count_trust_beta, side_eff, kappa, od_g, od_r,
                        gdna_prior_count, n_grid)
    return regions, left, right
