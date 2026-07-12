"""Phase-2 gDNA-density mixture prior — the nonparametric ``P(log ρ_g)`` trained on the Phase-1 single-strand
solutions (design: ``docs/calibration/PHASE2_gdna_mixture_fit_design.md``).

Two objects:

* :func:`build_training_substrate` (P2.0) — extract the training-node set from a solved chain: per-node
  ``log ρ_g`` (gDNA density) + a reliability ``weight``. Training nodes are the single-strand + structural
  region nodes (incl. zero-count intergenic/intronic — the depleted-floor anchor) and the clean-exon
  boundary crossings; **AMBIG nodes are excluded** (they are the target — their ``f_g`` is the unknown, solved
  in Phase 3 from the trained prior). This keeps the fit non-circular.
* :class:`GdnaDensityPrior` — a weighted Gaussian KDE in log space with a swappable bandwidth estimator
  (Silverman / likelihood-CV / fixed), pre-evaluated on a fine grid; the per-node solve uses its true-kernel
  :meth:`GdnaDensityPrior.logpdf_kernel` (real quadratic tails), NOT the clamped :meth:`logpdf` interpolation.

This is the PRODUCTION Phase-2 gDNA-density prior: it is trained on the pass-1 solved belief in
:func:`calibrate.calibrate` and added to the per-node ψ in pass 2 via ``bp_solver._kde_logprior``.
``scripts/debug/plot_gdna_prior.py`` visualises the fit against the oracle.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .node_chain import BOUNDARY, REGION
from .node_geometry import _node_region_type, node_densities

__all__ = ["TrainingSubstrate", "GdnaDensityPrior", "build_training_substrate"]

_EPS = 1.0e-9

# node-kind codes for the substrate (diagnostics / plot colouring)
KIND_INTERGENIC = 0
KIND_INTRON = 1
KIND_EXON = 2
KIND_BOUNDARY = 3
KIND_NAMES = {0: "intergenic", 1: "intron", 2: "exon", 3: "boundary"}


@dataclass(frozen=True, slots=True)
class TrainingSubstrate:
    """The training-node set: one row per training node.

    ``log_rho`` — the node's gDNA density ``log ρ_g = log(f_g·M/E_gdna)``, floored at the min-observable
    ``1/E_gdna`` (a finite count can never give an exactly-zero density). ``weight`` — **UNIT** (all 1.0):
    each solved node is one observation of a genomic location's density; precision is NOT used as a weight
    (it correlates with the density and would bias the distribution's shape — design §8e). ``node_kind`` —
    the ``KIND_*`` code. ``node_index`` — the chain node id (cross-reference to the belief/geometry).
    ``log_rho_std`` — the per-node log-density noise scale (= ``sqrt(Var(log f_g) + 1/(gcount+1))``); the KDE
    bandwidth is floored at its weighted median so a tight cluster is not resolved into spurious modes."""

    log_rho: np.ndarray
    weight: np.ndarray
    node_kind: np.ndarray
    node_index: np.ndarray
    log_rho_std: np.ndarray

    @property
    def n(self) -> int:
        return int(self.log_rho.shape[0])

    @property
    def n_eff(self) -> float:
        """Kish effective sample size ``(Σw)²/Σw²`` — the weighting-aware count used by the bandwidth rules."""
        w = self.weight
        s2 = float(np.sum(w) ** 2)
        return s2 / max(float(np.sum(w * w)), _EPS)


def _clean_exon_boundary(chain, region_arrays, boundary_substrate):
    """(clean_exon_bnd, exon_on_right) per chain node — the intron/intergenic↔exon crossings, exon-facing side.
    Mirrors the cleanliness test in ``bp_solver._gdna_seed_estimate``; the two live in different modules, so a
    small shared-helper dedup remains as an optional follow-up (both compute the same clean-boundary mask)."""
    kind = np.asarray(chain.kind)
    is_bnd = kind == BOUNDARY
    idx = np.asarray(chain.ref_idx, dtype=np.int64)
    _, rtype = _node_region_type(chain, region_arrays)
    R = rtype.shape[0]
    blr = np.asarray(boundary_substrate.left_region, dtype=np.int64)
    brr = np.asarray(boundary_substrate.right_region, dtype=np.int64)
    Bn = blr.shape[0]
    bi_ = np.clip(idx, 0, Bn - 1)
    lt = np.where((blr[bi_] >= 0) & is_bnd, rtype[np.clip(blr[bi_], 0, R - 1)], -1)
    rt = np.where((brr[bi_] >= 0) & is_bnd, rtype[np.clip(brr[bi_], 0, R - 1)], -1)
    left_clean = (lt == 0) | (lt == 1)
    right_clean = (rt == 0) | (rt == 1)
    clean_exon_bnd = is_bnd & ((left_clean & (rt == 2)) | (right_clean & (lt == 2)))
    exon_on_right = clean_exon_bnd & (rt == 2)
    return clean_exon_bnd, exon_on_right


def build_training_substrate(
    chain,
    belief,
    geometry,
    statics,
    region_arrays,
    boundary_substrate,
    *,
    min_eff_length: float = 0.0,
    include_boundaries: bool = False,
) -> TrainingSubstrate:
    """Extract the Phase-2 training-node set from a solved chain (P2.0).

    The unified solver (structure + strand deconvolution + message/belief propagation + the extremely weak
    gDNA floor) has ALREADY solved every node; this just READS its output. Training nodes = every SOLVED
    non-AMBIG node — intergenic sinks + single-strand introns + single-strand exons (and, if
    ``include_boundaries``, the clean-exon boundary crossings), INCLUDING zero-count intergenic/intronic (the
    depleted-floor anchor: a zero count over a large eff-length is a real, if imprecise, observation that the
    density is low). AMBIG nodes (both strands free) are excluded — they are the target. Unsolved nodes
    (``Var(log f_g)=∞`` — no strand and no mass) are excluded (no observation).

    The training-node density is the solve's own ``ρ_g = f_g·M/E_gdna`` (:func:`node_densities`), floored at
    the min-observable ``1/E_gdna``. **Every training node carries UNIT weight** — NOT precision weight:
    precision correlates with the density (the depleted floor is low-count/low-precision, the enriched mode
    high-count/high-precision), so precision-weighting would systematically down-weight the whole floor MODE
    and bias the distribution's SHAPE (measured — design §8e). The solve's log-density noise scale
    (``log_rho_std``) is kept only to floor the KDE bandwidth, not to weight.

    ``min_eff_length`` excludes region nodes whose gDNA eff-length ``E_gdna`` is below it — a region shorter
    than a fragment cannot CONTAIN one, so ``E_gdna→0`` and its contained-density ``1/E`` blows up spuriously
    (design §8e). The caller passes the gDNA mean fragment length; such tiny regions' gDNA is still observed
    via boundary crossings, so nothing is lost from the (undefined) contained density."""
    kind = np.asarray(chain.kind)
    is_reg = kind == REGION
    is_bnd = kind == BOUNDARY
    node_rtype, _ = _node_region_type(chain, region_arrays)

    fp = np.asarray(statics.free_pos, bool)
    fn = np.asarray(statics.free_neg, bool)
    is_ambig = (
        fp & fn
    )  # both strands free (G3) — EXCLUDED (target, not training; non-circular split)

    dens = node_densities(belief, geometry)
    fg = np.asarray(belief.f_g, dtype=np.float64)
    var_g = np.asarray(
        belief.var_gdna, dtype=np.float64
    )  # Var(log f_g); ∞ ⇒ unsolved (no observation)
    solved = np.isfinite(var_g)
    Ml = np.asarray(geometry.mass_left, dtype=np.float64)
    Mr = np.asarray(geometry.mass_right, dtype=np.float64)
    EGl = np.maximum(np.asarray(geometry.eff_gdna_left, dtype=np.float64), _EPS)
    EGr = np.maximum(np.asarray(geometry.eff_gdna_right, dtype=np.float64), _EPS)

    def _std(
        gcount,
    ):  # per-node log-density noise scale √(Var(log f_g)+1/(gcount+1)) — for the bandwidth floor
        return np.sqrt(np.maximum(var_g, 0.0) + 1.0 / (np.maximum(gcount, 0.0) + 1.0))

    # ---- region training nodes: every SOLVED non-AMBIG region large enough to contain a fragment
    #      (E_gdna ≥ min_eff_length), INCLUDING zero-count intergenic/intronic (the depleted-floor anchor). ----
    reg_train = is_reg & ~is_ambig & solved & (EGl >= float(min_eff_length))
    rho_reg = np.maximum(np.asarray(dens.rho_g_left, dtype=np.float64), 1.0 / EGl)
    kind_reg = np.where(
        node_rtype == 2, KIND_EXON, np.where(node_rtype == 1, KIND_INTRON, KIND_INTERGENIC)
    )

    # ---- boundary training nodes (clean intron/intergenic↔exon crossing; exon-facing side) — off by default;
    #      the crossing-density normalization biases these low (dissection §8c). The compute is skipped
    #      entirely when off (bnd_train is then an all-False mask and the collect section is a no-op on it);
    #      toggle ``include_boundaries`` to experiment. ----
    if include_boundaries:
        clean_exon_bnd, exon_on_right = _clean_exon_boundary(
            chain, region_arrays, boundary_substrate
        )
        M_bnd = np.where(exon_on_right, Mr, Ml)
        E_bnd = np.where(exon_on_right, EGr, EGl)
        rho_bnd = np.where(exon_on_right, np.asarray(dens.rho_g_right), np.asarray(dens.rho_g_left))
        rho_bnd = np.maximum(rho_bnd, 1.0 / np.maximum(E_bnd, _EPS))
        bnd_train = (
            is_bnd & clean_exon_bnd & solved & (M_bnd > 0.0) & (E_bnd >= float(min_eff_length))
        )
    else:
        M_bnd = np.zeros(int(chain.n_nodes), dtype=np.float64)
        rho_bnd = np.ones(int(chain.n_nodes), dtype=np.float64)
        bnd_train = np.zeros(int(chain.n_nodes), dtype=bool)

    # ---- collect ----
    node_idx = np.arange(int(chain.n_nodes))
    keep = np.zeros(int(chain.n_nodes), dtype=bool)
    log_rho = np.zeros(int(chain.n_nodes), dtype=np.float64)
    weight = np.zeros(int(chain.n_nodes), dtype=np.float64)
    std = np.zeros(int(chain.n_nodes), dtype=np.float64)
    nkind = np.full(int(chain.n_nodes), -1, dtype=np.int64)

    log_rho[reg_train] = np.log(rho_reg[reg_train])
    std[reg_train] = _std(fg * Ml)[reg_train]
    nkind[reg_train] = kind_reg[reg_train]
    keep |= reg_train

    log_rho[bnd_train] = np.log(rho_bnd[bnd_train])
    std[bnd_train] = _std(fg * M_bnd)[bnd_train]
    nkind[bnd_train] = KIND_BOUNDARY
    keep |= bnd_train

    weight[keep] = (
        1.0  # UNIT weight (design §8e) — NOT precision; noise is handled by the bandwidth
    )
    keep &= np.isfinite(log_rho)
    return TrainingSubstrate(
        log_rho=log_rho[keep],
        weight=weight[keep],
        node_kind=nkind[keep],
        node_index=node_idx[keep],
        log_rho_std=std[keep],
    )


# ---------------------------------------------------------------------------
# P2.1 — the nonparametric fit.
# ---------------------------------------------------------------------------


def _weighted_kde_logpdf(x_eval, x, w, h, *, chunk: int = 4096, eval_chunk: int = 16384):
    """log of the weighted Gaussian KDE ``Σ w_j φ((x−x_j)/h)/(h·Σw)`` at ``x_eval``. Returns ``(n_eval,)``.

    Tiles BOTH axes: the query axis (n_eval) in ``eval_chunk`` blocks and the sample axis (n_samp) in
    ``chunk`` blocks, so the ``(n_eval, n_samp)`` product never materialises. n_eval reaches ~10^8 at genome
    scale (nodes × solve grid); without the query tiling this allocated TiB and OOM-crashed. Both are pure
    memory knobs — the log-sum-exp is exact regardless of tiling, so the result is unchanged."""
    x = np.asarray(x, dtype=np.float64)
    w = np.asarray(w, dtype=np.float64)
    xe = np.asarray(x_eval, dtype=np.float64)
    W = max(float(np.sum(w)), _EPS)
    lw = np.log(np.maximum(w, _EPS)) - np.log(W) - np.log(max(h, _EPS)) - 0.5 * np.log(2.0 * np.pi)
    hh = max(float(h), _EPS)
    out = np.empty(xe.shape[0], dtype=np.float64)
    for e0 in range(0, xe.shape[0], eval_chunk):
        xeb = xe[e0 : e0 + eval_chunk][:, None]  # (b, 1)
        # accumulate logsumexp over sample chunks (running max for stability)
        m = np.full(xeb.shape[0], -np.inf)
        s = np.zeros(xeb.shape[0])
        for s0 in range(0, x.shape[0], chunk):
            xj = x[s0 : s0 + chunk][None, :]
            lwj = lw[s0 : s0 + chunk][None, :]
            z = -0.5 * ((xeb - xj) / hh) ** 2 + lwj  # (b, chunk)
            cmax = z.max(axis=1)
            both = np.maximum(m, cmax)
            s = s * np.exp(m - both) + np.sum(np.exp(z - both[:, None]), axis=1)
            m = both
        out[e0 : e0 + eval_chunk] = m + np.log(np.maximum(s, _EPS))
    return out


def _weighted_median(v, w):
    """Weighted median of ``v`` with weights ``w`` — CONTINUOUS (interpolated) so a machine-ε change in the
    weights or values moves the result by machine-ε, not discretely.

    The old ``vs[searchsorted(cw, 0.5·tot)]`` picked the sample AT the 0.5-mass crossing, a step function of
    the inputs: a 1e-15 nudge could flip the crossing to the adjacent sample and jump the result — this fed the
    KDE bandwidth floor and was a cross-process nondeterminism amplifier. Interpolate the value at cumulative
    mass 0.5 on the mid-rank CDF (``(cumsum−½w)/tot``) instead; ``np.interp`` is continuous in the data."""
    v = np.asarray(v, dtype=np.float64)
    w = np.asarray(w, dtype=np.float64)
    if v.size == 0:
        return 0.0
    order = np.argsort(v, kind="stable")
    vs, ws = v[order], w[order]
    tot = float(np.sum(ws))
    if tot <= _EPS:
        return float(np.median(v))
    cw_mid = (np.cumsum(ws) - 0.5 * ws) / tot  # cumulative weight at each sample's mass midpoint
    return float(np.interp(0.5, cw_mid, vs))


def _silverman_bandwidth(x, w):
    """Silverman's rule with a weighting-aware (Kish) effective count and a robust scale
    (``min(std, IQR/1.349)``)."""
    x = np.asarray(x, dtype=np.float64)
    w = np.asarray(w, dtype=np.float64)
    n_eff = float(np.sum(w) ** 2) / max(float(np.sum(w * w)), _EPS)
    mu = float(np.sum(w * x) / max(float(np.sum(w)), _EPS))
    std = float(np.sqrt(max(np.sum(w * (x - mu) ** 2) / max(float(np.sum(w)), _EPS), 0.0)))
    # weighted IQR
    order = np.argsort(x)
    xs, ws = x[order], w[order]
    cw = np.cumsum(ws) / max(float(np.sum(ws)), _EPS)
    q1 = float(np.interp(0.25, cw, xs))
    q3 = float(np.interp(0.75, cw, xs))
    scale = min(std, (q3 - q1) / 1.349) if q3 > q1 else std
    scale = scale if scale > _EPS else max(std, _EPS)
    return 0.9 * scale * max(n_eff, 1.0) ** (-0.2)


def _lscv_bandwidth(x, w, candidates, *, max_n: int = 4000, seed_stride: int | None = None):
    """Leave-one-out likelihood cross-validation: pick ``h`` maximising ``Σ_i w_i·log P̂_{−i}(x_i)``. O(n²)
    per candidate, so sub-sample to ``max_n`` points (deterministic stride — no RNG in the solve path)."""
    x = np.asarray(x, dtype=np.float64)
    w = np.asarray(w, dtype=np.float64)
    n = x.shape[0]
    if n > max_n:
        stride = seed_stride or max(1, n // max_n)
        sel = np.arange(0, n, stride)
        x, w = x[sel], w[sel]
        n = x.shape[0]
    W = max(float(np.sum(w)), _EPS)
    d = x[:, None] - x[None, :]  # (n,n)
    np.fill_diagonal(d, np.inf)  # leave-one-out: drop self
    best_h, best_ll = float(candidates[0]), -np.inf
    for h in candidates:
        # P̂_{−i}(x_i) = Σ_{j≠i} w_j φ(d_ij/h) / (h·(W−w_i))
        k = np.exp(-0.5 * (d / max(h, _EPS)) ** 2) / (max(h, _EPS) * np.sqrt(2.0 * np.pi))
        dens = (k * w[None, :]).sum(axis=1) / np.maximum(W - w, _EPS)
        ll = float(np.sum(w * np.log(np.maximum(dens, _EPS))))
        if ll > best_ll:
            best_ll, best_h = ll, float(h)
    return best_h


def _find_modes(x_grid, logP_grid):
    """Local maxima of ``logP`` (diagnostic only — the solve never hard-assigns a mode). Returns
    ``[(x, logP)]`` sorted by height descending."""
    p = np.asarray(logP_grid)
    interior = np.where((p[1:-1] > p[:-2]) & (p[1:-1] >= p[2:]))[0] + 1
    modes = [(float(x_grid[i]), float(p[i])) for i in interior]
    modes.sort(key=lambda t: t[1], reverse=True)
    return modes


@dataclass(frozen=True, slots=True)
class GdnaDensityPrior:
    """The fitted ``P(log ρ_g)`` — a weighted log-space Gaussian KDE, pre-evaluated on a fine grid.

    :meth:`logpdf` (linear interpolation of ``logP_grid``) is the Phase-3 hook: added per grid point to the
    per-node ψ, replacing the Gaussian ``bp_solver._global_logprior`` currently builds. The remaining fields
    are diagnostics for the P2.2 plotting framework."""

    x_grid: np.ndarray
    logP_grid: np.ndarray
    bandwidth: float
    n_eff: float
    train_x: np.ndarray
    train_w: np.ndarray
    train_kind: np.ndarray
    modes: tuple
    #: MIXTURE-BRIDGE weight ε∈[0,1) (Fix 1; `CalibrationConfig.gdna_prior_mixture_bridge`). The KDE is
    #: estimated from clean (unimodal) region nodes, so it has a deep VALLEY between the depleted and enriched
    #: modes; a node whose current-belief gDNA density is a spatial MIXTURE (a capture boundary, a
    #: sparse-probe region) lands in that valley by construction and collapses to f_g≈0. `_kde_logprior`
    #: mixes the KDE with a uniform "any-mixture" bridge over the observed density support at weight ε — this
    #: floors the valley (no collapse) while leaving the KDE's real tails OUTSIDE the support intact (so the
    #: high-density false-positive suppression is unchanged). ε=0 ⇒ no bridge (bit-exact legacy KDE). The
    #: level is robust (the peak/valley gap is ~10² nats), so any small ε defeats the collapse cliff. Design:
    #: `docs/calibration/boundary_kde_valley_collapse_and_simplex_precision.md`.
    mixture_bridge: float = 0.0
    #: Bridge SUPPORT TRIM (percent; `CalibrationConfig.calib_kde_bridge_trim_pct`) — the mixture bridge is
    #: bounded to the ``[trim, 100−trim]`` percentiles of the training gDNA-density support (a robustness trim so
    #: outliers don't set the bridge range). Consumed by `bp_solver._kde_logprior`.
    bridge_trim_pct: float = 0.5

    def logpdf(self, log_rho) -> np.ndarray:
        """``log P̂(log ρ_g)`` by linear interpolation; the flat boundary value extrapolates outside the grid
        (= "no prior information beyond the observed range" — the tilt / messages / residual then decide).

        NOTE: the clamped tails make this UNSUITABLE for the per-node SOLVE — a high-density node whose
        implied ``log ρ_g`` falls in the (constant) upper tail sees no penalty and drifts to ``f_g≈0.5`` →
        catastrophic false-positive gDNA. Use :meth:`logpdf_kernel` (real quadratic tails) in the solve;
        keep this interpolation only for plotting / diagnostics."""
        x = np.asarray(log_rho, dtype=np.float64)
        return np.interp(
            x,
            self.x_grid,
            self.logP_grid,
            left=float(self.logP_grid[0]),
            right=float(self.logP_grid[-1]),
        )

    def logpdf_kernel(self, log_rho) -> np.ndarray:
        """``log P̂(log ρ_g)`` by the TRUE weighted-Gaussian kernel sum — real quadratic tails (unlike the
        clamped :meth:`logpdf` interpolation). The solve MUST use this: a genuine ``−½((x−xᵢ)/h)²`` tail
        penalises implausibly-high gDNA density instead of drifting to a constant, which is what removes the
        density-prior false-positive gDNA (verified: clamped tails give ~585k FP on gDNA-free unstranded data
        vs ~2k with real tails). Same kernel/bandwidth/weights the fit uses, so it agrees with the plotted
        curve inside the observed range."""
        x = np.asarray(log_rho, dtype=np.float64)
        return _weighted_kde_logpdf(x.ravel(), self.train_x, self.train_w, self.bandwidth).reshape(
            x.shape
        )

    @classmethod
    def fit(
        cls,
        substrate: TrainingSubstrate,
        *,
        bandwidth="silverman",
        n_grid: int = 1024,
        pad: float = 1.0,
        floor_log_rho: float | None = None,
        floor_weight: float = 0.0,
        n_lscv: int = 40,
        mixture_bridge: float = 0.0,
        bridge_trim_pct: float = 0.5,
    ) -> "GdnaDensityPrior":
        """Fit the KDE.

        ``bandwidth`` — ``'silverman'`` | ``'lscv'`` | a float (the fixed ``h``). ``floor_log_rho`` /
        ``floor_weight`` — optionally seed the depleted mode with a virtual sample at ``log ρ_floor``
        (design §2.4). ``pad`` — grid margin (in units of the final bandwidth) beyond the sample range.
        ``mixture_bridge`` — the Fix-1 valley-fill weight ε (stored on the prior; consumed by
        ``bp_solver._kde_logprior``).
        """
        x = np.asarray(substrate.log_rho, dtype=np.float64).copy()
        w = np.asarray(substrate.weight, dtype=np.float64).copy()
        kind = np.asarray(substrate.node_kind).copy()
        std = np.asarray(substrate.log_rho_std, dtype=np.float64).copy()
        if floor_log_rho is not None and floor_weight > 0.0:
            x = np.append(x, float(floor_log_rho))
            w = np.append(w, float(floor_weight))
            kind = np.append(kind, KIND_INTERGENIC)
            std = np.append(std, float(np.median(std)) if std.size else 0.0)
        if x.shape[0] == 0:
            raise ValueError("empty training substrate — no training nodes")

        if isinstance(bandwidth, str):
            if bandwidth == "silverman":
                h = _silverman_bandwidth(x, w)
            elif bandwidth == "lscv":
                h0 = _silverman_bandwidth(x, w)
                cand = h0 * np.geomspace(0.25, 4.0, int(n_lscv))
                h = _lscv_bandwidth(x, w, cand)
            else:
                raise ValueError(f"unknown bandwidth estimator {bandwidth!r}")
        else:
            h = float(bandwidth)
        # NOISE FLOOR (no magic number): the KDE must not resolve below each node's own Poisson
        # log-density uncertainty, else a tight (uniform-gDNA) cluster fractures into spurious modes. Floor h
        # at the weighted median per-node sampling std. This is what makes the estimators robust off-capture.
        h_floor = _weighted_median(std, w) if std.size else 0.0
        h = max(float(h), float(h_floor), _EPS)

        lo = float(np.min(x)) - pad * h
        hi = float(np.max(x)) + pad * h
        if hi <= lo:
            hi = lo + 1.0
        x_grid = np.linspace(lo, hi, int(n_grid))
        logP = _weighted_kde_logpdf(x_grid, x, w, h)
        n_eff = float(np.sum(w) ** 2) / max(float(np.sum(w * w)), _EPS)
        return cls(
            x_grid=x_grid,
            logP_grid=logP,
            bandwidth=h,
            n_eff=n_eff,
            train_x=x,
            train_w=w,
            train_kind=kind,
            modes=tuple(_find_modes(x_grid, logP)),
            mixture_bridge=float(mixture_bridge),
            bridge_trim_pct=float(bridge_trim_pct),
        )
