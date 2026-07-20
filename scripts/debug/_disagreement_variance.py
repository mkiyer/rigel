"""RETIRED total-density σ²_imp helpers — relocated out of production (`bp_solver`) into this debug util
(Phase-0 cleanup, `docs/calibration/background_reference_implementation_plan.md`). These estimated the old
adjacent-node disagreement variance that the NPMLE projection σ²_transfer replaced; kept ONLY for the historical
σ²_transfer investigation harnesses that still import them. Not part of any production path.
"""

from __future__ import annotations

import numpy as np

from rigel.calibration.node_chain import BOUNDARY, NodeChain
from rigel.calibration.node_geometry import NodeGeometry

_EPS = 1.0e-12

__all__ = [
    "_poisson_moment_var",
    "_adjacent_edges",
    "_adjacent_log_density_residuals",
    "adjacent_disagreement_variance",
]


def _poisson_moment_var(resid, ns, nd) -> float:
    """The Poisson disagreement-variance moment estimator (the total-density σ²_imp).
    ``Var(d) = σ²_imp + 1/n_i + 1/n_j`` ⇒ subtract each pair's known Poisson sampling and
    inverse-variance-weight-average (``w = nᵢnⱼ/(nᵢ+nⱼ)``, harmonic — 0 at zero count, no threshold).
    ``resid`` must already be oriented to a single sign convention; it is median-centered here."""
    resid = np.asarray(resid, float)
    ns = np.asarray(ns, float)
    nd = np.asarray(nd, float)
    ok = np.isfinite(resid) & (ns > _EPS) & (nd > _EPS)
    resid, ns, nd = resid[ok], ns[ok], nd[ok]
    if resid.size < 2:
        return 1.0
    dc = resid - np.median(resid)  # remove the systematic frame offset
    w = (ns * nd) / (ns + nd)
    den = float(np.sum(w))
    if den <= _EPS:
        return 1.0
    return float(max(np.sum(w * (dc * dc - (1.0 / ns + 1.0 / nd))) / den, 0.0))


def _adjacent_edges(chain: NodeChain):
    """Adjacent boundary↔region edges as (src, dst, src_is_boundary), each undirected edge once (i→right[i]).
    src presents its RIGHT face to dst; dst its LEFT face."""
    kind = np.asarray(chain.kind)
    right = np.asarray(chain.right)
    src = np.arange(kind.shape[0])
    dst = right
    m = dst >= 0
    src, dst = src[m], dst[m]
    return src, dst, kind[src] == BOUNDARY


def _adjacent_log_density_residuals(chain: NodeChain, geometry: NodeGeometry):
    """Oriented adjacent boundary↔region log gDNA-density residuals + the two source/dst gDNA counts. ``ρ =
    mass / eff_gdna`` oriented boundary→region so the systematic frame offset is one median-removable mode."""
    ML = np.asarray(geometry.mass_left)
    MR = np.asarray(geometry.mass_right)
    EGL = np.asarray(geometry.eff_gdna_left)
    EGR = np.asarray(geometry.eff_gdna_right)
    src, dst, s_bnd = _adjacent_edges(chain)
    n_i, e_i, n_j, e_j = MR[src], EGR[src], ML[dst], EGL[dst]
    ok = (n_i > _EPS) & (n_j > _EPS) & (e_i > _EPS) & (e_j > _EPS)
    n_i, e_i, n_j, e_j, s_bnd = n_i[ok], e_i[ok], n_j[ok], e_j[ok], s_bnd[ok]
    lr_i = np.log(n_i / e_i)
    lr_j = np.log(n_j / e_j)
    resid = np.where(s_bnd, lr_i - lr_j, lr_j - lr_i)  # orient boundary→region (one mode)
    return resid, n_i, n_j


def adjacent_disagreement_variance(chain: NodeChain, geometry: NodeGeometry) -> float:
    """The v1 total-density σ²_imp: the mean adjacent-node naive-gDNA log-density disagreement variance
    (Poisson sampling removed) — the population-average message process variance."""
    resid, n_i, n_j = _adjacent_log_density_residuals(chain, geometry)
    return _poisson_moment_var(resid, n_i, n_j)
