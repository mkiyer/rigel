"""Pass-0 gDNA KDE prototype — the UNIFICATION, validated on real cached data (no production change).

The KDE (`GdnaDensityPrior.fit` + `logpdf_kernel`) and the projection (`bp_solver._kde_logprior`:
log ρ_g = log(f_g·M/E) → KDE kernel → prior) are ALREADY pass-agnostic. The only pass-specific piece is
the training substrate. This builds the **pass-0 confident substrate** (structural gDNA, no solve) and
reuses the production KDE + projection:

    confident set (f_g=1 by structure) = intergenic regions · single-strand intron regions ·
      intergenic↔exon boundaries · SJ (single-strand intron↔exon) boundaries with spliced=0

It fits `GdnaDensityPrior` on that substrate and plots (a) the fitted KDE over the confident-density
histogram — does it model the depleted+enriched landscape — and (b) the PROJECTION function
`observed density → prior-preferred f_g` (argmax over the f_g grid of the KDE-kernel + Jeffreys prior,
the exact `_kde_logprior` integrand). That projection curve is the deliverable: it shows at what observed
density the prior calls gDNA vs RNA.

    python scripts/debug/pass0_kde_prototype.py --cache CACHE.pkl[,...] --out DIR
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from rigel.calibration.gdna_density_prior import GdnaDensityPrior, TrainingSubstrate
from rigel.calibration.node_chain import REGION, build_node_chain
from rigel.calibration.node_geometry import build_node_geometry, build_node_statics
from rigel.calibration.signature import (
    BIT_EXON_NEG,
    BIT_EXON_POS,
    BIT_INTRON_NEG,
    BIT_INTRON_POS,
    TS_NEG,
    TS_POS,
    transcript_strand_class,
)
from rigel.calibration.substrate import BoundarySubstrate, CalibrationSubstrate

sys.path.insert(0, str(Path(__file__).parent))
from calib_cache import load  # noqa: E402
from pass0_kde_landscape import _sj_boundary_masks  # noqa: E402

_EPS = 1.0e-9
_EXON = BIT_EXON_POS | BIT_EXON_NEG
_INTRON = BIT_INTRON_POS | BIT_INTRON_NEG


def _build_structures(inp):
    """Chain / geometry / statics / substrates directly from the cached payload — NO sweep (fast)."""
    payload, ra = inp["payload"], inp["region_arrays"]
    sub = CalibrationSubstrate.from_payload(payload, ra)
    bsub = BoundarySubstrate.from_payload(payload)
    chain = build_node_chain(payload.ref_region_offsets, payload.ref_boundary_offsets)
    geom = build_node_geometry(chain, sub, bsub, ra, inp["gdna_fl_pmf"], inp["rna_fl_pmf"])
    st = build_node_statics(chain, sub, bsub, ra)
    return chain, geom, st, sub, bsub, ra


def build_pass0_substrate(chain, geom, ra, bsub, *, min_eff_length: float) -> TrainingSubstrate:
    """THE NEW PIECE: the pass-0 confident-node substrate (observed density, f_g=1 by structure).
    Everything downstream (GdnaDensityPrior.fit, _kde_logprior projection) is the unchanged production KDE."""
    kind, ref = np.asarray(chain.kind), np.asarray(chain.ref_idx, np.int64)
    is_reg = kind == REGION
    is_bnd = ~is_reg
    egl, egr = np.asarray(geom.eff_gdna_left), np.asarray(geom.eff_gdna_right)
    msl, msr = np.asarray(geom.mass_left), np.asarray(geom.mass_right)
    eff = np.where(is_reg, egl, 0.5 * (egl + egr))
    mass = np.where(is_reg, msl, msl + msr)
    dens = np.maximum(
        mass / np.maximum(eff, _EPS), 1.0 / np.maximum(eff, _EPS)
    )  # f_g=1, floored 1/E

    sig = np.asarray(ra.signature).astype(np.int64)
    sc = transcript_strand_class(sig)
    R = sig.shape[0]
    ri = np.clip(ref, 0, R - 1)
    sig_r = np.where(is_reg, sig[ri], 0)
    single = (sc[ri] == TS_POS) | (sc[ri] == TS_NEG)
    intergenic_reg = is_reg & (sig_r == 0)
    ss_intron_reg = is_reg & single & ((sig_r & _INTRON) != 0) & ((sig_r & _EXON) == 0)

    B = np.asarray(bsub.left_region).shape[0]
    bi = np.clip(ref, 0, B - 1)
    lr = np.asarray(bsub.left_region, np.int64)
    rr = np.asarray(bsub.right_region, np.int64)
    sig_l = np.where(lr >= 0, sig[np.clip(lr, 0, None)], 0)
    sig_rr = np.where(rr >= 0, sig[np.clip(rr, 0, None)], 0)
    ex_l, ex_r = (sig_l & _EXON) != 0, (sig_rr & _EXON) != 0
    interg_exon_b = ((sig_l == 0) & ex_r) | ((sig_rr == 0) & ex_l)
    spl_b = np.asarray(bsub.left.mass_spliced, float) + np.asarray(bsub.right.mass_spliced, float)
    sj0_b = _sj_boundary_masks(bsub, ra) & (spl_b <= _EPS)  # SJ spliced=0
    interg_exon_bnd = is_bnd & interg_exon_b[bi]
    sj0_bnd = is_bnd & sj0_b[bi]

    keep = (
        (intergenic_reg | ss_intron_reg | interg_exon_bnd | sj0_bnd)
        & (mass > 0.0)
        & (eff >= min_eff_length)
    )
    log_rho = np.log(dens[keep])
    gcount = mass[keep]  # f_g=1 ⇒ gDNA count ≈ observed mass
    std = np.sqrt(1.0 / (gcount + 1.0))  # no var_g pre-solve — count noise only (bandwidth floor)
    nkind = np.where(is_reg[keep], 0, 1).astype(np.int64)
    return TrainingSubstrate(
        log_rho=log_rho,
        weight=np.ones(int(keep.sum())),
        node_kind=nkind,
        node_index=np.where(keep)[0],
        log_rho_std=std,
    ), dict(
        intergenic_reg=intergenic_reg,
        ss_intron_reg=ss_intron_reg,
        interg_exon_bnd=interg_exon_bnd,
        sj0_bnd=sj0_bnd,
        dens=dens,
        is_reg=is_reg,
        sig_r=sig_r,
    )


def _projection_curve(prior: GdnaDensityPrior, log_d_grid, fg_grid):
    """The `_kde_logprior` projection, prior-only (no strand/messages). For each observed log-density d the
    prior over f_g is ``p(f_g|d) ∝ KDE.logpdf_kernel(log(fg·d)) + Jeffreys(1/(1−fg))`` — the KDE mixture
    evaluated at the candidate gDNA density (spatial-likelihood weighting: nearby kernels dominate). We read
    the POSTERIOR MEAN (what the solve integrates), NOT the argmax — the mean is smooth; the argmax cliffs.
    Returns ``(mean, p10, p90)`` per d (the band shows the spread — a node is a MIXTURE, not a point)."""
    fg = np.clip(fg_grid, _EPS, 1 - _EPS)
    jeff = -np.log1p(-fg)  # (K,)
    mean = np.empty_like(log_d_grid)
    p10 = np.empty_like(log_d_grid)
    p90 = np.empty_like(log_d_grid)
    for i, ld in enumerate(log_d_grid):
        psi = prior.logpdf_kernel(np.log(fg) + ld) + jeff
        w = np.exp(psi - psi.max())
        w /= w.sum()
        mean[i] = float((w * fg).sum())
        cdf = np.cumsum(w)
        p10[i] = fg[int(np.searchsorted(cdf, 0.10))]
        p90[i] = fg[int(np.searchsorted(cdf, 0.90))]
    return mean, p10, p90


def run(cache_path: str, outdir: Path, fl_mean: float = 200.0):
    inp = load(cache_path)
    chain, geom, st, sub, bsub, ra = _build_structures(inp)
    substrate, masks = build_pass0_substrate(chain, geom, ra, bsub, min_eff_length=fl_mean)
    prior = GdnaDensityPrior.fit(substrate, bandwidth="silverman", mixture_bridge=0.01)
    cond = Path(cache_path).stem

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    # A: fitted KDE over the confident-density histogram
    x = substrate.log_rho / np.log(10)  # log10 for display
    axes[0].hist(
        x, bins=60, density=True, alpha=0.45, color="tab:blue", label=f"confident (n={substrate.n})"
    )
    xs = np.linspace(x.min(), x.max(), 400)
    logp = prior.logpdf_kernel(xs * np.log(10))  # kernel is in natural-log units
    axes[0].plot(
        xs, np.exp(logp) * np.log(10), "k-", lw=1.8, label=f"pass-0 KDE (h={prior.bandwidth:.2f})"
    )
    axes[0].set_title(f"{cond}: pass-0 KDE over confident set")
    axes[0].set_xlabel("log10 gDNA density")
    axes[0].legend(fontsize=8)
    # B: the projection function observed-density -> prior-preferred f_g. Span the confident range and
    #    1.5 decades ABOVE it (to show what the prior does with high-density exon nodes beyond the support).
    ld_grid = np.linspace(x.min() - 0.5, x.max() + 1.5, 120) * np.log(10)
    fg_grid = np.linspace(1e-3, 1 - 1e-3, 400)
    fg_mean, fg_p10, fg_p90 = _projection_curve(prior, ld_grid, fg_grid)
    xd = ld_grid / np.log(10)
    axes[1].fill_between(xd, fg_p10, fg_p90, color="red", alpha=0.15, label="10–90% posterior")
    axes[1].plot(xd, fg_mean, "r-", lw=2, label="posterior mean f_g")
    axes[1].legend(fontsize=8)
    axes[1].set_title("projection: observed density → posterior f_g (smooth)")
    axes[1].set_xlabel("log10 observed density")
    axes[1].set_ylabel("prior-preferred f_g")
    axes[1].set_ylim(-0.05, 1.05)
    axes[1].axhline(0.5, color="gray", lw=0.5, ls=":")
    fig.tight_layout()
    p = outdir / f"pass0_kde_{cond}.png"
    fig.savefig(p, dpi=110)
    plt.close(fig)
    print(
        f"{cond:10s} n_confident={substrate.n} bandwidth={prior.bandwidth:.3f} "
        f"median_log10dens={np.median(x):.2f} -> {p}",
        flush=True,
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cache", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)
    for c in args.cache.split(","):
        run(c, outdir)


if __name__ == "__main__":
    main()
