"""Three figures: (1) σ²_bio fit with LOESS added (vs monotone + free); (2) the ê(z) enrichment-transfer
fit + the z-compression root cause; (3) dispersion-vs-disagreement (is μ or |Δlogρ| the real driver)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

from rigel.calibration.variance_model import MonotoneVarMean, MonotoneMean, _setup_spline, _select_lambda
from rigel.index import TranscriptIndex
from rigel.config import BamScanConfig, CalibrationConfig
from rigel.pipeline import scan_and_buffer
from rigel.calibration.substrate import BoundarySubstrate, CalibrationSubstrate
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.node_chain import build_node_chain, REGION
from rigel.calibration import bp_solver as bp
from rigel.calibration.strand_balance import fit_strand_balance
from rigel.calibration.signature import coarse_type_from_signature
from rigel.splice import SpliceType

_EPS = 1e-12
ROOT = Path.home() / "Downloads/rigel_runs/quick_3to1_5mb"
SCR = Path("/private/tmp/claude-503/-Users-mkiyer-proj-rigel/9e5820a2-4635-49bb-9a62-f9debab6095c/scratchpad")
_BISQUARE_C = 6.0


def _loess(x, y, xq, frac=0.4, robust_iters=2):
    """Faithful historical robust LOESS (density_model.py @200143d4): local-linear, tricube kernel,
    adaptive k-NN bandwidth, bisquare(6·MAD) robustness."""
    x = np.asarray(x, float); y = np.asarray(y, float); n = x.shape[0]
    k = int(np.clip(int(frac * n), 2, n)); rob = np.ones(n)

    def fit_at(q):
        out = np.empty(q.shape[0])
        for i, xi in enumerate(q):
            d = np.abs(x - xi); h = np.partition(d, k - 1)[k - 1]
            w = np.clip(1.0 - (d / max(h, 1e-12)) ** 3, 0.0, 1.0) ** 3 * rob
            sw = w.sum()
            if sw <= 0.0:
                out[i] = np.average(y, weights=rob) if rob.sum() > 0 else y.mean(); continue
            mx = (w * x).sum() / sw; my = (w * y).sum() / sw
            sxx = (w * (x - mx) ** 2).sum()
            b = (w * (x - mx) * (y - my)).sum() / sxx if sxx > 1e-12 else 0.0
            out[i] = my + b * (xi - mx)
        return out
    for _ in range(robust_iters):
        resid = y - fit_at(x); mad = np.median(np.abs(resid))
        if mad <= 0.0:
            break
        rob = np.clip(1.0 - (resid / (_BISQUARE_C * mad)) ** 2, 0.0, 1.0) ** 2
    return fit_at(np.asarray(xq, float))


def _free_fit_offset(cls, mean, raw, offset, weight=None, *, k=18, degree=3, robust_iters=2, n_lambda=40, population_spread=False):
    mean = np.asarray(mean, float); raw = np.asarray(raw, float); off = np.asarray(offset, float)
    wt = np.ones_like(off) if weight is None else np.asarray(weight, float)
    ok = (np.isfinite(mean) & (mean > _EPS) & np.isfinite(raw) & (raw >= 0) & np.isfinite(off) & (off >= 0) & (wt > 0))
    mean, raw, off, wt = mean[ok], raw[ok], off[ok], wt[ok]
    if mean.size < max(k, 8):
        return cls._constant_offset(mean, raw, off, degree, wt)
    o = np.argsort(mean); mean, raw, off, wt = mean[o], raw[o], off[o], wt[o]
    x = np.log(mean); z = raw - off
    kn, B, P, xlo, xhi = _setup_spline(x, k, degree); lam, edf = _select_lambda(B, z, P, mean.size, n_lambda)
    total = np.maximum(off, _EPS); coeffs = None
    for it in range(robust_iters + 1):
        w = wt / np.maximum(total, _EPS) ** 2
        if it > 0:
            r = (z - B @ coeffs) * np.sqrt(w); s = 1.4826 * float(np.median(np.abs(r - np.median(r)))) + _EPS
            u = r / (4.685 * s); w = w * np.where(np.abs(u) < 1.0, (1.0 - u ** 2) ** 2, 0.0)
        coeffs = np.linalg.solve(B.T @ (w[:, None] * B) + lam * P, B.T @ (w * z))
        if it == robust_iters:
            break
        total = np.maximum(B @ coeffs, 0.0) + off
    return cls(knots=kn, degree=degree, coeffs=coeffs, x_lo=xlo, x_hi=xhi, fit_mean=mean, fit_var=np.maximum(z, 0), edf=float(edf), lam=float(lam))


_orig = MonotoneVarMean.fit_offset


def setup(cond):
    idx = TranscriptIndex.load(str(ROOT / "rigel_index"))
    st, sm, flm, buf, pl = scan_and_buffer(str(ROOT / cond / "sim_oracle.bam"), idx, BamScanConfig())
    fl = build_fl_models(global_counts=flm.global_model.counts, rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts, gdna_counts=gdna_fl_mass(pl), max_size=flm.max_size)
    ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
    sub = CalibrationSubstrate.from_payload(pl, ra); bsub = BoundarySubstrate.from_payload(pl)
    chain = build_node_chain(pl.ref_region_offsets, pl.ref_boundary_offsets)
    geom = bp.build_node_geometry(chain, sub, bsub, ra, fl.gdna_pmf, fl.rna_pmf)
    statics = bp.build_node_statics(chain, sub, bsub, ra)
    kappa = float(fit_strand_balance(sm).rna_sense_frac)
    b = bp.init_beliefs(chain, sub, bsub, ra, rna_sense_frac=kappa, n_grid=60, statics=statics)
    b, _ = bp.node_sweep(chain, statics, geom, b, ra, bsub, rna_sense_frac=kappa, n_grid=60, max_outer=CalibrationConfig().sweep_max_outer)
    return dict(idx=idx, sm=sm, ra=ra, sub=sub, bsub=bsub, chain=chain, geom=geom, statics=statics, kappa=kappa, belief=b)


def gdna_edges(S):
    dens = bp.node_densities(S["belief"], S["geom"])
    m, r, o = bp._edge_varmean(S["chain"], np.asarray(dens.rho_g_left), np.asarray(dens.rho_g_right),
                               np.asarray(S["geom"].eff_gdna_left), np.asarray(S["geom"].eff_gdna_right), np.asarray(S["statics"].strand_obs, bool))
    mean = np.concatenate(m); raw = np.concatenate(r); off = np.concatenate(o)
    ok = (mean > _EPS) & (raw >= 0) & (off >= 0)
    return mean[ok], raw[ok], off[ok], dens


# ============ FIGURE 1: σ²_bio with LOESS ============
CONDS = [("Capture-ON, stranded (ss0.99)", "gdna_gdna300_ss_0.99_nrna_none_capture_on"),
         ("Capture-ON, UNSTRANDED (ss0.50)", "gdna_gdna300_ss_0.50_nrna_none_capture_on"),
         ("Capture-OFF control (ss0.99)", "gdna_gdna300_ss_0.99_nrna_none_capture_off")]
flag_S = None
fig, axes = plt.subplots(1, 3, figsize=(21, 7.0))
for ax, (title, cond) in zip(axes, CONDS):
    MonotoneVarMean.fit_offset = _orig
    S = setup(cond)
    if cond.endswith("ss_0.99_nrna_none_capture_on"):
        flag_S = S
    mean, raw, off, dens = gdna_edges(S)
    z = raw - off
    mono = bp.fit_gdna_varmean(S["chain"], dens, S["geom"], S["statics"])
    MonotoneVarMean.fit_offset = classmethod(_free_fit_offset)
    free = MonotoneVarMean.fit_offset(mean, raw, off)
    MonotoneVarMean.fit_offset = _orig
    o = np.argsort(mean)
    loess_curve = lambda g: np.maximum(_loess(np.log(mean[o]), z[o], np.log(g)), 0.0)
    qs = np.quantile(mean, np.linspace(0, 1, 15)); bx, bmed, blo, bhi = [], [], [], []
    for a, bb in zip(qs[:-1], qs[1:]):
        mm = (mean >= a) & (mean < bb) if bb < qs[-1] else (mean >= a) & (mean <= bb)
        if mm.sum() >= 3:
            bx.append(np.median(mean[mm])); bmed.append(np.median(z[mm])); blo.append(np.percentile(z[mm], 25)); bhi.append(np.percentile(z[mm], 75))
    bx, bmed, blo, bhi = map(np.array, (bx, bmed, blo, bhi))
    grid = np.logspace(np.log10(max(mean.min(), 1e-5)), np.log10(mean.max()), 300)
    pm, pf, pl_ = mono.predict(grid), free.predict(grid), loess_curve(grid)
    ytop = max(np.max(bhi) * 1.15, pm.max(), pf.max(), pl_.max()) * 1.05
    ax.scatter(mean, np.clip(z, -2, ytop), s=5, alpha=0.07, color="0.55", rasterized=True)
    ax.fill_between(bx, blo, bhi, color="0.7", alpha=0.4, label="data IQR (per μ-bin)")
    ax.plot(bx, bmed, "o-", color="black", lw=2.2, ms=5, label="empirical median (true shape)", zorder=5)
    ax.plot(grid, pm, color="#d62728", lw=2.6, ls="--", label=f"MONOTONE (production) edf={mono.edf:.1f}", zorder=6)
    ax.plot(grid, pf, color="#1f77b4", lw=2.6, label=f"FREE spline edf={free.edf:.1f}", zorder=6)
    ax.plot(grid, pl_, color="#2ca02c", lw=2.6, ls=":", label="LOESS (span 0.4, robust)", zorder=6)
    ax.set_xscale("log"); ax.set_xlabel("μ = ½(ρ_dst+ρ_src)  [log]"); ax.set_ylabel("σ²_bio  /  z=(Δlogρ)²−offset")
    ax.set_ylim(-2, ytop); ax.set_title(f"{title}\nκ={S['kappa']:.3f}, n={mean.size:,}", fontsize=10.5)
    ax.grid(True, which="both", alpha=0.2); ax.legend(fontsize=8.5, loc="upper left", framealpha=0.9)
fig.suptitle("σ²_bio(μ): production MONOTONE vs FREE spline vs LOESS — all three are 1-D μ-smoothers", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.94]); fig.savefig(SCR / "fig1_sigma2_loess.png", dpi=130); plt.close(fig)
print("saved fig1_sigma2_loess.png")

# ============ FIGURE 3: dispersion vs disagreement (flagship) ============
S = flag_S
mean, raw, off, dens = gdna_edges(S)
dlog = np.sqrt(np.maximum(raw, 0.0))  # |Δlogρ| (raw = dlog²)
fig, ax = plt.subplots(1, 2, figsize=(15, 6.2))
sc0 = ax[0].scatter(mean, np.clip(raw, 0, 60), c=np.clip(dlog, 0, 8), s=10, alpha=0.5, cmap="viridis")
ax[0].set_xscale("log"); ax[0].set_xlabel("μ = ½(ρ_dst+ρ_src)  [log]"); ax[0].set_ylabel("(Δlogρ)²  (clipped 60)")
ax[0].set_title("dispersion vs μ (the fit's predictor)\ncolored by |Δlogρ| — at a given μ, raw scatters with the disagreement"); ax[0].grid(alpha=0.2)
fig.colorbar(sc0, ax=ax[0], label="|Δlogρ| (disagreement)")
sc1 = ax[1].scatter(dlog, np.clip(raw, 0, 60), c=np.log10(np.maximum(mean, 1e-5)), s=10, alpha=0.5, cmap="plasma")
ax[1].set_xlabel("|Δlogρ|  (cross-node disagreement)"); ax[1].set_ylabel("(Δlogρ)²  (clipped 60)")
ax[1].set_title("dispersion vs |Δlogρ| (the disagreement)\ncolored by log μ — tight law; μ varies freely ALONG it"); ax[1].grid(alpha=0.2)
fig.colorbar(sc1, ax=ax[1], label="log₁₀ μ")
fig.suptitle("Is the dispersion a function of the LEVEL (μ) or the DISAGREEMENT (|Δlogρ|)?  [flagship cap-ON ss0.99]", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.93]); fig.savefig(SCR / "fig3_disagreement.png", dpi=130); plt.close(fig)
print("saved fig3_disagreement.png")

# ============ FIGURE 2: ê(z) enrichment transfer (flagship) ============
cap = {}
_omm = MonotoneMean.fit.__func__
def _wrap(cls, x, y, weight=None, *, recal_weight=None, **k):
    cap["z"] = np.asarray(x, float).copy(); cap["rho"] = np.asarray(y, float).copy()
    cap["rw"] = None if recal_weight is None else np.asarray(recal_weight, float).copy()
    return _omm(cls, x, y, weight, recal_weight=recal_weight, **k)
MonotoneMean.fit = classmethod(_wrap)
rho_global, gdna_vm, var_mean = bp._gdna_seed_estimate(S["chain"], S["statics"], S["geom"], S["ra"], S["bsub"], S["belief"].f_g, S["kappa"])
ehat, z_all, sig2lvl, w_enr = bp.fit_enrichment_transfer(S["chain"], S["statics"], S["geom"], S["ra"], S["bsub"], S["belief"].f_g, S["kappa"], float(rho_global))
MonotoneMean.fit = classmethod(_omm)
# contained interior density T = Ml/EGl over the single-strand EXON fit set (the ~true level, workflow corr 0.99)
kind = np.asarray(S["chain"].kind); is_reg = kind == REGION
nrt, rtype = bp._node_region_type(S["chain"], S["ra"])
EGl = np.asarray(S["geom"].eff_gdna_left); Ml = np.asarray(S["geom"].mass_left)
T = Ml / np.maximum(EGl, _EPS)
ss_exon = is_reg & (nrt == 2) & np.asarray(S["statics"].strand_obs, bool) & (Ml > 0) & np.isfinite(z_all)
zE = z_all[ss_exon]; TE = T[ss_exon]
zf = cap["z"]; rho = cap["rho"]
fig, ax = plt.subplots(1, 2, figsize=(15, 6.2))
gz = np.logspace(np.log10(max(zf.min(), 1e-5)), np.log10(zf.max()), 200)
ax[0].scatter(zf, rho, s=14, alpha=0.5, color="0.4", label="fit points (z, ρ_g)")
ax[0].plot(gz, ehat.predict(gz), color="#d62728", lw=2.8, label=f"ê(z) fitted (w_enrich={ehat.w_enrich:.2f}, scale={ehat.scale:.2f})")
ax[0].axhline(rho_global, color="0.5", ls=":", label=f"ρ_global={rho_global:.3f}")
ax[0].set_xscale("log"); ax[0].set_yscale("log"); ax[0].set_xlabel("z = boundary-crossing gDNA density (the predictor)")
ax[0].set_ylabel("ρ_g  (response: exon-body gDNA density)"); ax[0].set_title("ê(z) enrichment-transfer fit"); ax[0].grid(True, which="both", alpha=0.2); ax[0].legend(fontsize=8.5)
lim = [min(zE.min(), TE.min()) * 0.7, max(zE.max(), TE.max()) * 1.4]
ax[1].scatter(TE, zE, s=14, alpha=0.5, color="#1f77b4")
ax[1].plot(lim, lim, "k--", lw=1.5, label="z = contained (no compression)")
ax[1].set_xscale("log"); ax[1].set_yscale("log"); ax[1].set_xlim(lim); ax[1].set_ylim(lim)
ax[1].set_xlabel("contained interior density M/E  (≈ true level, corr 0.99 w/ true E)")
ax[1].set_ylabel("z = boundary-crossing density (the predictor)")
ax[1].set_title("ROOT CAUSE: z COMPRESSES — at high enrichment z flattens\nwhile the true (contained) level keeps rising → ê(z) plateaus"); ax[1].grid(True, which="both", alpha=0.2); ax[1].legend(fontsize=8.5)
fig.suptitle("Enrichment transfer ê(z): the fit (left) and why the PREDICTOR z is impoverished (right)  [flagship cap-ON ss0.99]", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.93]); fig.savefig(SCR / "fig2_ehat.png", dpi=130); plt.close(fig)
print(f"saved fig2_ehat.png  (n_fit={zf.size}, n_ss_exon={int(ss_exon.sum())})")
