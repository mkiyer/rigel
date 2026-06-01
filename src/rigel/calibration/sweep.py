"""AMBIG boundary↔region sweep imputation (PR 4b §I.2–I.5, D7 leg 2).

A strand-ambiguous (AMBIG) region cannot use the strand channel, so its gDNA
exposure ``ω`` is imputed from strand-**decodable** neighbours by propagating
gDNA-density evidence along the **alternating region↔boundary chain**
(``… B → R → B → R …``), decaying smoothly across each boundary by a purity
weight. Every node carries a gDNA Gamma pair ``(α, β) = (gDNA count, effective
length)`` whose ratio estimates the local density ``ρ_0·ω``:

* **region node** ``(M_g_contained, region_eff_len)`` — contained gDNA, FL-corrected
  contained exposure (:mod:`rigel.calibration.effective_length`).
* **boundary node** ``(gDNA crossing flux, μ_FL)`` — crossing gDNA, the gДNA FL
  mean (region-independent).

The two passes pool decayed neighbour evidence; the immediately adjacent boundary
nodes (one decay each) dominate (nearest-neighbour priority). The imputed AMBIG
``ω`` is the **same** closed-form Gamma posterior used elsewhere, fed the swept
``(α, β)`` (doc 03 §4). It only overwrites AMBIG-region exposure; decodable
regions keep their PR 4a exposure. The gDNA **mass** attribution (D1, per side)
is untouched — the sweep borrows *density*, not mass.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .signature import TS_AMBIG

# Unit pseudocounts for the boundary purity weight (Q6: weakly-informative,
# one pseudo-fragment each; not cliffs). See PR04b §I.3.
_TRAFFIC_PSEUDOCOUNT = 1.0
_SPLICED_PSEUDOCOUNT = 1.0


@dataclass(frozen=True, slots=True)
class SweepResult:
    """Exposure after the AMBIG sweep: AMBIG rows imputed, others copied."""

    omega: np.ndarray  # float64[R]
    log_omega_var: np.ndarray  # float64[R]
    n_ambig: int  # AMBIG regions whose exposure was imputed


def _internal_boundary_nodes(
    ref_id: np.ndarray,
    alloc_left,
    alloc_right,
    sub_left,
    sub_right,
) -> tuple[np.ndarray, np.ndarray]:
    """Per-region gDNA crossing count + purity weight of the boundary to its RIGHT.

    Indexed by the **left** region ``r``: entry ``r`` describes the internal
    boundary between region ``r`` and ``r + 1`` (defined iff ``r + 1`` is in the
    same reference). The crossing flux/π_g come from the two adjacent per-side
    views (PR04b §I.2/§I.3); the two sides of a contiguous crossing carry the
    same fragments, so each quantity is the average of the two sides.
    """
    r_total = ref_id.shape[0]
    has_right = np.zeros(r_total, dtype=bool)
    has_right[:-1] = ref_id[:-1] == ref_id[1:]

    # left side of the boundary = region r's RIGHT view; right side = (r+1)'s LEFT view.
    u_left = sub_right.n_unspliced.astype(np.float64)  # flux, region r side
    s_left = sub_right.n_spliced.astype(np.float64)
    pi_left = alloc_right.pi_g
    u_right = np.empty(r_total, dtype=np.float64)
    s_right = np.empty(r_total, dtype=np.float64)
    pi_right = np.empty(r_total, dtype=np.float64)
    u_right[:-1] = sub_left.n_unspliced[1:].astype(np.float64)
    s_right[:-1] = sub_left.n_spliced[1:].astype(np.float64)
    pi_right[:-1] = alloc_left.pi_g[1:]
    u_right[-1] = 0.0
    s_right[-1] = 0.0
    pi_right[-1] = 0.0

    u_b = 0.5 * (u_left + u_right)
    s_b = 0.5 * (s_left + s_right)
    purity = 0.5 * (pi_left + pi_right)
    gdna_flux = 0.5 * (pi_left * u_left + pi_right * u_right)

    traffic = u_b / (u_b + _TRAFFIC_PSEUDOCOUNT)
    rna_clean = _SPLICED_PSEUDOCOUNT / (_SPLICED_PSEUDOCOUNT + s_b)
    weight = traffic * rna_clean * purity
    weight[~has_right] = 0.0
    gdna_flux[~has_right] = 0.0
    return gdna_flux, weight


def sweep_ambig_exposure(
    substrate,
    region_arrays,
    *,
    alloc_contained,
    alloc_left,
    alloc_right,
    region_eff_len: np.ndarray,
    mu_fl: float,
    rho_0: float,
    exposure_dispersion: float,
    base_omega: np.ndarray,
    base_log_omega_var: np.ndarray,
) -> SweepResult:
    """Impute AMBIG-region exposure via the alternating region↔boundary sweep."""
    ts_class = np.asarray(substrate.ts_class)
    ref_id = np.asarray(region_arrays.ref_id)
    ref_offsets = np.asarray(region_arrays.ref_offsets, dtype=np.int64)
    inv_dispersion = 1.0 / exposure_dispersion

    omega = np.array(base_omega, dtype=np.float64, copy=True)
    log_omega_var = np.array(base_log_omega_var, dtype=np.float64, copy=True)
    if not np.any(ts_class == TS_AMBIG) or mu_fl <= 0.0:
        return SweepResult(omega=omega, log_omega_var=log_omega_var, n_ambig=0)

    # Region nodes: gDNA contained count + FL-corrected contained exposure.
    a_reg = np.asarray(alloc_contained.m_g_unspl, dtype=np.float64)
    b_reg = np.asarray(region_eff_len, dtype=np.float64)
    # Boundary-to-the-right nodes (indexed by left region): crossing gDNA + weight.
    a_bnd, weight = _internal_boundary_nodes(
        ref_id, alloc_left, alloc_right, substrate.left, substrate.right
    )

    r_total = ts_class.shape[0]
    from_left_a = np.zeros(r_total, dtype=np.float64)
    from_left_b = np.zeros(r_total, dtype=np.float64)
    from_right_a = np.zeros(r_total, dtype=np.float64)
    from_right_b = np.zeros(r_total, dtype=np.float64)

    for f in range(ref_offsets.shape[0] - 1):
        start, end = int(ref_offsets[f]), int(ref_offsets[f + 1])
        if end <= start:
            continue
        # Forward: crossing a boundary scales the inflow AND the boundary's own
        # node by the purity weight w (run ← w·(run + boundary_node)). Weighting
        # both α and β preserves each contribution's density (α/β) while scaling
        # its evidence weight by reliability — so a choked boundary fades out
        # rather than injecting μ_FL exposure with no gDNA. Then record the
        # region's inflow and inject its own (full-weight) contained node.
        run_a = run_b = 0.0
        for r in range(start, end):
            if r > start:
                w = weight[r - 1]  # boundary between R[r-1] and R[r]
                run_a = w * (run_a + a_bnd[r - 1])
                run_b = w * (run_b + mu_fl)
            from_left_a[r] = run_a
            from_left_b[r] = run_b
            run_a += a_reg[r]
            run_b += b_reg[r]
        # Reverse.
        run_a = run_b = 0.0
        for r in range(end - 1, start - 1, -1):
            if r < end - 1:
                w = weight[r]  # boundary between R[r] and R[r+1]
                run_a = w * (run_a + a_bnd[r])
                run_b = w * (run_b + mu_fl)
            from_right_a[r] = run_a
            from_right_b[r] = run_b
            run_a += a_reg[r]
            run_b += b_reg[r]

    alpha_swept = a_reg + from_left_a + from_right_a
    beta_swept = b_reg + from_left_b + from_right_b
    ambig = ts_class == TS_AMBIG
    alpha_post = inv_dispersion + alpha_swept[ambig]
    omega[ambig] = alpha_post / (inv_dispersion + rho_0 * beta_swept[ambig])
    log_omega_var[ambig] = 1.0 / alpha_post
    return SweepResult(omega=omega, log_omega_var=log_omega_var, n_ambig=int(ambig.sum()))


__all__ = ["SweepResult", "sweep_ambig_exposure"]
