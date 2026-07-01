"""Phase-2 gDNA-density mixture prior — the nonparametric ``P(log ρ_g)`` trained on the Phase-1 single-strand
solutions (design: ``docs/calibration/PHASE2_gdna_mixture_fit_design.md``).

Two objects:

* :func:`build_training_substrate` (P2.0) — extract the teacher set from a solved chain: per-node
  ``log ρ_g`` (gDNA density) + a reliability ``weight``. Teachers are the single-strand + structural region
  nodes and the clean-exon boundary crossings; **AMBIG nodes are excluded** (they are the students — their
  ``f_g`` is the unknown, so they cannot teach). This keeps the fit non-circular.
* :class:`GdnaDensityPrior` (P2.1) — a weighted Gaussian KDE in log space with a swappable bandwidth
  estimator (Silverman / likelihood-CV / fixed), pre-evaluated on a fine grid and consumed by the per-node
  solve via :meth:`GdnaDensityPrior.logpdf` (the drop-in for ``bp_solver._global_logprior``'s Gaussian).

Nothing here is wired into the solve yet (that is P2.4); the immediate consumer is the P2.2 plotting
framework, which visualises the fit against the oracle before any bandwidth is chosen for production.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .bp_solver import BOUNDARY, REGION, _node_region_type, node_densities

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
    """The teacher set: one row per teaching node.

    ``log_rho`` — the node's gDNA density ``log ρ_g = log(f_g·M/E_gdna)``, floored at the min-observable
    ``1/E_gdna`` (a finite count can never give an exactly-zero density). ``weight`` — the reliability
    ``strand_discriminability·precision`` (§1 of the design; no new constant). ``node_kind`` — the
    ``KIND_*`` code. ``node_index`` — the chain node id (cross-reference to the belief/geometry).
    ``log_rho_std`` — the per-node sampling std of ``log ρ_g`` (= ``sqrt(Var(log f_g) + 1/(gcount+1))``, the
    density-noise scale); the KDE bandwidth is floored at its weighted median so a tight (uniform-gDNA)
    cluster is not resolved into spurious modes."""

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
    Mirrors the cleanliness test in ``bp_solver._gdna_seed_estimate`` (kept local to avoid a shipped-file
    refactor at this stage; consolidate when Phase 2 lands)."""
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


def _exon_spliced_residual(chain, geometry, statics, clean_exon_bnd, node_rtype):
    """``ρ_resid = clip(M_unspliced − ρ_mature·E_rna, 0) / E_gdna`` per EXON region node — the STRAND-FREE,
    gDNA-specific density (the mature RNA subtracted via the flanking clean-exon boundary's motif-stranded
    spliced mass). NaN on non-exon nodes and exons with no flanking spliced. Mirrors the ``rho_spliced``
    computation in :func:`bp_solver.fit_enrichment_transfer` (the strand-immune estimator promoted here to
    the unstranded exon TEACHER density, per design D2)."""
    left = np.asarray(chain.left)
    right = np.asarray(chain.right)
    fp = np.asarray(statics.free_pos, bool)  # single-strand: motif strand is + iff free_pos (else −)
    SPl, SPr = geometry.spliced_pos_left, geometry.spliced_pos_right
    SNl, SNr = geometry.spliced_neg_left, geometry.spliced_neg_right
    ESPl, ESPr = geometry.eff_spl_left, geometry.eff_spl_right
    ERl = np.asarray(geometry.eff_rna_left, dtype=np.float64)
    EGl = np.maximum(np.asarray(geometry.eff_gdna_left, dtype=np.float64), _EPS)
    Ml = np.asarray(geometry.mass_left, dtype=np.float64)
    out = np.full(int(chain.n_nodes), np.nan, dtype=np.float64)
    for i in np.where(node_rtype == 2)[0]:
        sl, sr = (SPl, SPr) if fp[i] else (SNl, SNr)  # spliced on this exon's motif strand
        m_spl, e_spl = 0.0, _EPS
        lb = int(left[i])
        if lb >= 0 and clean_exon_bnd[lb] and sr[lb] > m_spl:  # exon is lb's RIGHT region → lb's right side
            m_spl, e_spl = float(sr[lb]), float(ESPr[lb])
        rb = int(right[i])
        if rb >= 0 and clean_exon_bnd[rb] and sl[rb] > m_spl:  # exon is rb's LEFT region → rb's left side
            m_spl, e_spl = float(sl[rb]), float(ESPl[rb])
        if m_spl > 0.0:
            rho_mature = m_spl / max(e_spl, _EPS)
            out[i] = max(float(Ml[i]) - rho_mature * float(ERl[i]), 0.0) / float(EGl[i])
    return out


def build_training_substrate(
    chain,
    belief,
    geometry,
    statics,
    region_arrays,
    boundary_substrate,
    kappa,
    *,
    include_boundaries: bool = True,
) -> TrainingSubstrate:
    """Extract the Phase-2 teacher set from a solved chain (P2.0).

    Teachers = non-AMBIG region nodes with mass (intergenic / single-strand introns / single-strand exons)
    plus, if ``include_boundaries``, the clean-exon boundary crossings (exon-facing side). AMBIG nodes (both
    strands free) are excluded — they are the students. The density is ``ρ_g = f_g·M/E_gdna`` (from
    :func:`node_densities`), floored at the min-observable ``1/E_gdna``. The weight is
    ``strand_factor·precision`` where ``strand_factor = (2κ−1)²`` on single-strand nodes (fades to 0 at
    κ→½, so unstranded data leans on the structural intergenic teachers) and ``1`` on structural nodes
    (intergenic sinks, clean-exon boundaries), and ``precision = 1/(Var(log f_g) + 1/(gcount+1))``."""
    kind = np.asarray(chain.kind)
    is_reg = kind == REGION
    is_bnd = kind == BOUNDARY
    node_rtype, _ = _node_region_type(chain, region_arrays)

    fp = np.asarray(statics.free_pos, bool)
    fn = np.asarray(statics.free_neg, bool)
    is_ss = fp ^ fn  # exactly one strand free (G2 single-strand — strand-derived density)
    is_ambig = fp & fn  # both free (G3) — EXCLUDED (students; non-circular teacher/target split)

    dens = node_densities(belief, geometry)
    fg = np.asarray(belief.f_g, dtype=np.float64)
    var_g = np.maximum(np.asarray(belief.var_gdna, dtype=np.float64), 0.0)  # Var(log f_g)
    Ml = np.asarray(geometry.mass_left, dtype=np.float64)
    Mr = np.asarray(geometry.mass_right, dtype=np.float64)
    EGl = np.maximum(np.asarray(geometry.eff_gdna_left, dtype=np.float64), _EPS)
    EGr = np.maximum(np.asarray(geometry.eff_gdna_right, dtype=np.float64), _EPS)
    w_str = float((2.0 * kappa - 1.0) ** 2)

    def _precision(gcount):
        return 1.0 / (var_g + 1.0 / (np.maximum(gcount, 0.0) + 1.0))

    def _std(gcount):  # per-node sampling std of log ρ_g (= 1/sqrt(precision))
        return np.sqrt(var_g + 1.0 / (np.maximum(gcount, 0.0) + 1.0))

    clean_exon_bnd, exon_on_right = _clean_exon_boundary(chain, region_arrays, boundary_substrate)

    # ---- region teachers (contained face; left == right for regions) ----
    # Density = a STRAND ⊕ ρ_resid blend (design D2): the strand-solved density (weight (2κ−1)²) blended with
    # the strand-free spliced-residual (weight 1−(2κ−1)², EXONS with flanking spliced only), so the enriched
    # exon teachers survive unstranded (κ→½) where the strand solve is worthless. Structural intergenic sinks
    # (neither strand free) keep strand-factor 1 (their M/E IS the gDNA density). denom→0 (unstranded exon with
    # NO spliced, or an unstranded intron) ⇒ dropped (no reliable teacher). weight = denom·precision.
    is_exon = node_rtype == 2
    rho_strand = np.maximum(np.asarray(dens.rho_g_left, dtype=np.float64), 1.0 / EGl)
    rho_resid_raw = _exon_spliced_residual(chain, geometry, statics, clean_exon_bnd, node_rtype)
    avail = is_exon & np.isfinite(rho_resid_raw) & (rho_resid_raw > 0.0)
    rho_resid = np.where(avail, np.maximum(rho_resid_raw, 1.0 / EGl), 0.0)
    w_strand_factor = np.where(is_ss, w_str, 1.0)          # 1 on structural intergenic; (2κ−1)² on SS
    w_spl_factor = (1.0 - w_strand_factor) * avail.astype(np.float64)
    denom = w_strand_factor + w_spl_factor
    rho_reg = np.where(denom > _EPS,
                       (w_strand_factor * rho_strand + w_spl_factor * rho_resid) / np.maximum(denom, _EPS),
                       rho_strand)
    reg_teacher = is_reg & (Ml > 0.0) & ~is_ambig & (denom > _EPS)
    w_reg = denom * _precision(fg * Ml)
    kind_reg = np.where(is_exon, KIND_EXON, np.where(node_rtype == 1, KIND_INTRON, KIND_INTERGENIC))

    # ---- boundary teachers (clean intron/intergenic↔exon crossing; exon-facing side) ----
    M_bnd = np.where(exon_on_right, Mr, Ml)
    E_bnd = np.where(exon_on_right, EGr, EGl)
    rho_bnd = np.where(exon_on_right, np.asarray(dens.rho_g_right), np.asarray(dens.rho_g_left))
    rho_bnd = np.maximum(rho_bnd, 1.0 / np.maximum(E_bnd, _EPS))
    w_bnd = _precision(fg * M_bnd)  # structural crossing → full trust
    bnd_teacher = is_bnd & clean_exon_bnd & (M_bnd > 0.0) & bool(include_boundaries)

    # ---- collect ----
    node_idx = np.arange(int(chain.n_nodes))
    keep = np.zeros(int(chain.n_nodes), dtype=bool)
    log_rho = np.zeros(int(chain.n_nodes), dtype=np.float64)
    weight = np.zeros(int(chain.n_nodes), dtype=np.float64)
    std = np.zeros(int(chain.n_nodes), dtype=np.float64)
    nkind = np.full(int(chain.n_nodes), -1, dtype=np.int64)

    log_rho[reg_teacher] = np.log(rho_reg[reg_teacher])
    weight[reg_teacher] = w_reg[reg_teacher]
    std[reg_teacher] = _std(fg * Ml)[reg_teacher]
    nkind[reg_teacher] = kind_reg[reg_teacher]
    keep |= reg_teacher

    log_rho[bnd_teacher] = np.log(rho_bnd[bnd_teacher])
    weight[bnd_teacher] = w_bnd[bnd_teacher]
    std[bnd_teacher] = _std(fg * M_bnd)[bnd_teacher]
    nkind[bnd_teacher] = KIND_BOUNDARY
    keep |= bnd_teacher

    keep &= np.isfinite(log_rho) & (weight > 0.0)
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


def _weighted_kde_logpdf(x_eval, x, w, h, *, chunk: int = 4096):
    """log of the weighted Gaussian KDE ``Σ w_j φ((x−x_j)/h)/(h·Σw)`` at ``x_eval`` (chunked over samples so
    the ``(n_eval, n_samp)`` product never materialises whole). Returns ``(n_eval,)``."""
    x = np.asarray(x, dtype=np.float64)
    w = np.asarray(w, dtype=np.float64)
    xe = np.asarray(x_eval, dtype=np.float64)
    W = max(float(np.sum(w)), _EPS)
    lw = np.log(np.maximum(w, _EPS)) - np.log(W) - np.log(max(h, _EPS)) - 0.5 * np.log(2.0 * np.pi)
    # accumulate logsumexp over sample chunks (running max for stability)
    m = np.full(xe.shape[0], -np.inf)
    s = np.zeros(xe.shape[0])
    for s0 in range(0, x.shape[0], chunk):
        xj = x[s0 : s0 + chunk][None, :]
        lwj = lw[s0 : s0 + chunk][None, :]
        z = -0.5 * ((xe[:, None] - xj) / max(h, _EPS)) ** 2 + lwj  # (n_eval, chunk)
        cmax = z.max(axis=1)
        both = np.maximum(m, cmax)
        s = s * np.exp(m - both) + np.sum(np.exp(z - both[:, None]), axis=1)
        m = both
    return m + np.log(np.maximum(s, _EPS))


def _weighted_median(v, w):
    """Weighted median of ``v`` with weights ``w``."""
    v = np.asarray(v, dtype=np.float64)
    w = np.asarray(w, dtype=np.float64)
    order = np.argsort(v)
    vs, ws = v[order], w[order]
    cw = np.cumsum(ws)
    tot = cw[-1] if cw.size else 0.0
    if tot <= _EPS:
        return float(np.median(v)) if v.size else 0.0
    return float(vs[np.searchsorted(cw, 0.5 * tot)])


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
    teacher_x: np.ndarray
    teacher_w: np.ndarray
    teacher_kind: np.ndarray
    modes: tuple

    def logpdf(self, log_rho) -> np.ndarray:
        """``log P̂(log ρ_g)`` by linear interpolation; the flat boundary value extrapolates outside the grid
        (= "no prior information beyond the observed range" — the tilt / messages / residual then decide)."""
        x = np.asarray(log_rho, dtype=np.float64)
        return np.interp(
            x,
            self.x_grid,
            self.logP_grid,
            left=float(self.logP_grid[0]),
            right=float(self.logP_grid[-1]),
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
    ) -> "GdnaDensityPrior":
        """Fit the KDE.

        ``bandwidth`` — ``'silverman'`` | ``'lscv'`` | a float (the fixed ``h``). ``floor_log_rho`` /
        ``floor_weight`` — optionally seed the depleted mode with a virtual sample at ``log ρ_floor``
        (design §2.4). ``pad`` — grid margin (in units of the final bandwidth) beyond the sample range.
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
            raise ValueError("empty training substrate — no teacher nodes")

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
            teacher_x=x,
            teacher_w=w,
            teacher_kind=kind,
            modes=tuple(_find_modes(x_grid, logP)),
        )
