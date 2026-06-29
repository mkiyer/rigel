"""The same 3 fit diagnostics on REAL VCaP (cached scan, genome-scale). fig1: σ²_bio monotone/free/LOESS;
fig2: ê(z) + z-compression; fig3: dispersion vs μ vs |Δlogρ|. Honest about the real-data complications
(eff_gdna spans 8 orders of mag → most edges at-floor; κ≈0 → strand signal dead)."""
import matplotlib, time, pickle
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

from rigel.calibration.variance_model import MonotoneVarMean, MonotoneMean, _setup_spline, _select_lambda
from rigel.index import TranscriptIndex
from rigel.config import CalibrationConfig
from rigel.calibration.substrate import BoundarySubstrate, CalibrationSubstrate
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.node_chain import build_node_chain, REGION
from rigel.calibration import bp_solver as bp
from rigel.calibration.strand_balance import fit_strand_balance
from rigel.splice import SpliceType

_EPS = 1e-12
SCR = Path("/private/tmp/claude-503/-Users-mkiyer-proj-rigel/9e5820a2-4635-49bb-9a62-f9debab6095c/scratchpad")
VCAP_PKL = SCR / "vcap_scan_fixed.pkl"
VCAP_IDX = Path.home() / "Downloads/rigel_runs/refs/rigel_index_v5"
rng = np.random.default_rng(0)
_BISQUARE_C = 6.0


def _loess(x, y, xq, frac=0.4, robust_iters=2):
    x = np.asarray(x, float); y = np.asarray(y, float); n = x.shape[0]
    k = int(np.clip(int(frac * n), 2, n)); rob = np.ones(n)
    def fit_at(q):
        out = np.empty(q.shape[0])
        for i, xi in enumerate(q):
            d = np.abs(x - xi); h = np.partition(d, k - 1)[k - 1]
            w = np.clip(1.0 - (d / max(h, 1e-12)) ** 3, 0.0, 1.0) ** 3 * rob; sw = w.sum()
            if sw <= 0.0:
                out[i] = y.mean(); continue
            mx = (w * x).sum() / sw; my = (w * y).sum() / sw; sxx = (w * (x - mx) ** 2).sum()
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
print("[vcap] loading cached scan ...", flush=True)
sm, flm, pl = pickle.load(open(VCAP_PKL, "rb"))
idx = TranscriptIndex.load(str(VCAP_IDX))
fl = build_fl_models(global_counts=flm.global_model.counts, rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts, gdna_counts=gdna_fl_mass(pl), max_size=flm.max_size)
ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
sub = CalibrationSubstrate.from_payload(pl, ra); bsub = BoundarySubstrate.from_payload(pl)
chain = build_node_chain(pl.ref_region_offsets, pl.ref_boundary_offsets)
geom = bp.build_node_geometry(chain, sub, bsub, ra, fl.gdna_pmf, fl.rna_pmf)
statics = bp.build_node_statics(chain, sub, bsub, ra)
kappa = float(fit_strand_balance(sm).rna_sense_frac)
print(f"[vcap] setup done; kappa={kappa:.5f}; nodes={len(np.asarray(chain.kind)):,}. sweeping (max_outer=1) ...", flush=True)
t0 = time.time()
b0 = bp.init_beliefs(chain, sub, bsub, ra, rna_sense_frac=kappa, n_grid=60, statics=statics)
belief, _ = bp.node_sweep(chain, statics, geom, b0, ra, bsub, rna_sense_frac=kappa, n_grid=60, max_outer=1)
print(f"[vcap] sweep done in {time.time()-t0:.0f}s", flush=True)
dens = bp.node_densities(belief, geom)

# gDNA edges (all single-strand) + off-floor flag
m, r, o = bp._edge_varmean(chain, np.asarray(dens.rho_g_left), np.asarray(dens.rho_g_right), np.asarray(geom.eff_gdna_left), np.asarray(geom.eff_gdna_right), np.asarray(statics.strand_obs, bool))
mean = np.concatenate(m); raw = np.concatenate(r); off = np.concatenate(o)
ok = (mean > _EPS) & (raw >= 0) & (off >= 0)
mean, raw, off = mean[ok], raw[ok], off[ok]
z = raw - off
print(f"[vcap] gDNA edges: {mean.size:,}  μ∈[{mean.min():.2e},{mean.max():.2e}]", flush=True)
mono = bp.fit_gdna_varmean(chain, dens, geom, statics)
MonotoneVarMean.fit_offset = classmethod(_free_fit_offset)
free = MonotoneVarMean.fit_offset(mean, raw, off)
MonotoneVarMean.fit_offset = _orig
# LOESS on a subsample (O(n²) robust pass infeasible at full n)
nss = min(4000, mean.size); ss = rng.choice(mean.size, nss, replace=False)
oss = np.argsort(mean[ss]); xl = np.log(mean[ss][oss]); zl = z[ss][oss]

# ---- FIG 1: σ²_bio ----
qs = np.quantile(mean, np.linspace(0, 1, 15)); bx, bmed, blo, bhi = [], [], [], []
for a, bb in zip(qs[:-1], qs[1:]):
    mm = (mean >= a) & (mean < bb) if bb < qs[-1] else (mean >= a) & (mean <= bb)
    if mm.sum() >= 3:
        bx.append(np.median(mean[mm])); bmed.append(np.median(z[mm])); blo.append(np.percentile(z[mm], 25)); bhi.append(np.percentile(z[mm], 75))
bx, bmed, blo, bhi = map(np.array, (bx, bmed, blo, bhi))
grid = np.logspace(np.log10(max(mean.min(), 1e-7)), np.log10(mean.max()), 250)
pm, pf = mono.predict(grid), free.predict(grid)
pl_ = np.maximum(_loess(xl, zl, np.log(grid)), 0.0)
ytop = max(np.max(bhi) * 1.15, pm.max(), pf.max(), float(np.nanmax(pl_))) * 1.05
sc = rng.choice(mean.size, min(8000, mean.size), replace=False)
fig, ax = plt.subplots(figsize=(11, 7))
ax.scatter(mean[sc], np.clip(z[sc], -2, ytop), s=5, alpha=0.10, color="0.55", rasterized=True)
ax.fill_between(bx, blo, bhi, color="0.7", alpha=0.4, label="data IQR (per μ-bin)")
ax.plot(bx, bmed, "o-", color="black", lw=2.2, ms=4, label="empirical median (true shape)", zorder=5)
ax.plot(grid, pm, color="#d62728", lw=2.6, ls="--", label=f"MONOTONE (production) edf={mono.edf:.1f}", zorder=6)
ax.plot(grid, pf, color="#1f77b4", lw=2.6, label=f"FREE spline edf={free.edf:.1f}", zorder=6)
ax.plot(grid, pl_, color="#2ca02c", lw=2.6, ls=":", label="LOESS (span 0.4, n=4k subsample)", zorder=6)
ax.set_xscale("log"); ax.set_xlabel("μ = ½(ρ_dst+ρ_src)  [gDNA density, log]"); ax.set_ylabel("σ²_bio  /  z=(Δlogρ)²−offset")
ax.set_ylim(-2, ytop); ax.grid(True, which="both", alpha=0.2); ax.legend(fontsize=9, loc="upper left")
ax.set_title(f"REAL VCaP — σ²_bio(μ): monotone vs free vs LOESS\nκ={kappa:.4f} (≈unstranded), n_edges={mean.size:,}  [eff_gdna spans ~8 orders of mag]", fontsize=11)
fig.tight_layout(); fig.savefig(SCR / "vcap_fig1_sigma2.png", dpi=130); plt.close(fig); print("saved vcap_fig1_sigma2.png", flush=True)

# ---- FIG 3: dispersion vs μ vs |Δlogρ| ----
dlog = np.sqrt(np.maximum(raw, 0.0))
fig, axx = plt.subplots(1, 2, figsize=(15, 6.2))
s2 = rng.choice(mean.size, min(12000, mean.size), replace=False)
c0 = axx[0].scatter(mean[s2], np.clip(raw[s2], 0, 40), c=np.clip(dlog[s2], 0, 6), s=8, alpha=0.4, cmap="viridis")
axx[0].set_xscale("log"); axx[0].set_xlabel("μ  [log]"); axx[0].set_ylabel("(Δlogρ)² (clip 40)"); axx[0].grid(alpha=0.2)
axx[0].set_title("dispersion vs μ (predictor), colored by |Δlogρ|"); fig.colorbar(c0, ax=axx[0], label="|Δlogρ|")
c1 = axx[1].scatter(dlog[s2], np.clip(raw[s2], 0, 40), c=np.log10(np.maximum(mean[s2], 1e-7)), s=8, alpha=0.4, cmap="plasma")
axx[1].set_xlabel("|Δlogρ| (disagreement)"); axx[1].set_ylabel("(Δlogρ)² (clip 40)"); axx[1].grid(alpha=0.2)
axx[1].set_title("dispersion vs |Δlogρ|, colored by log μ"); fig.colorbar(c1, ax=axx[1], label="log₁₀ μ")
fig.suptitle("REAL VCaP — is dispersion driven by the LEVEL (μ) or the DISAGREEMENT (|Δlogρ|)?", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.93]); fig.savefig(SCR / "vcap_fig3_disagreement.png", dpi=130); plt.close(fig); print("saved vcap_fig3_disagreement.png", flush=True)

# ---- FIG 2: ê(z) ----
print("[vcap] fitting enrichment transfer (Python region loops — may take ~30s) ...", flush=True)
cap = {}
_omm = MonotoneMean.fit.__func__
def _wrap(cls, x, y, weight=None, *, recal_weight=None, **k):
    cap["z"] = np.asarray(x, float).copy(); cap["rho"] = np.asarray(y, float).copy()
    return _omm(cls, x, y, weight, recal_weight=recal_weight, **k)
MonotoneMean.fit = classmethod(_wrap)
t0 = time.time()
rho_global, gdna_vm, var_mean = bp._gdna_seed_estimate(chain, statics, geom, ra, bsub, belief.f_g, kappa)
ehat, z_all, sig2lvl, w_enr = bp.fit_enrichment_transfer(chain, statics, geom, ra, bsub, belief.f_g, kappa, float(rho_global))
MonotoneMean.fit = classmethod(_omm)
print(f"[vcap] ê fit done in {time.time()-t0:.0f}s", flush=True)
kind = np.asarray(chain.kind); is_reg = kind == REGION
nrt, rtype = bp._node_region_type(chain, ra)
EGl = np.asarray(geom.eff_gdna_left); Ml = np.asarray(geom.mass_left); T = Ml / np.maximum(EGl, _EPS)
ss_exon = is_reg & (nrt == 2) & np.asarray(statics.strand_obs, bool) & (Ml > 0) & np.isfinite(z_all)
zE = z_all[ss_exon]; TE = T[ss_exon]
zf = cap.get("z", np.zeros(0)); rho = cap.get("rho", np.zeros(0))
fig, ax = plt.subplots(1, 2, figsize=(15, 6.2))
if zf.size:
    gz = np.logspace(np.log10(max(zf[zf > 0].min(), 1e-7)), np.log10(zf.max()), 200)
    ax[0].scatter(zf, np.maximum(rho, 1e-7), s=12, alpha=0.4, color="0.4", label=f"fit points (n={zf.size})")
    ax[0].plot(gz, ehat.predict(gz), color="#d62728", lw=2.8, label=f"ê(z) (w_enrich={ehat.w_enrich:.2f}, scale={ehat.scale:.2f})")
    ax[0].axhline(rho_global, color="0.5", ls=":", label=f"ρ_global={rho_global:.3g}")
    ax[0].set_xscale("log"); ax[0].set_yscale("log"); ax[0].set_xlabel("z = boundary-crossing gDNA density"); ax[0].set_ylabel("ρ_g (response)")
    ax[0].set_title("ê(z) enrichment-transfer fit"); ax[0].grid(True, which="both", alpha=0.2); ax[0].legend(fontsize=8.5)
if zE.size:
    pos = (zE > 0) & (TE > 0)
    lim = [min(zE[pos].min(), TE[pos].min()) * 0.7, max(zE[pos].max(), TE[pos].max()) * 1.4]
    ax[1].scatter(TE[pos], zE[pos], s=10, alpha=0.4, color="#1f77b4")
    ax[1].plot(lim, lim, "k--", lw=1.5, label="z = contained (no compression)")
    ax[1].set_xscale("log"); ax[1].set_yscale("log"); ax[1].set_xlim(lim); ax[1].set_ylim(lim)
    ax[1].set_xlabel("contained interior density M/E (≈ true level)"); ax[1].set_ylabel("z = boundary-crossing density")
    ax[1].set_title("z vs contained — does z compress at high enrichment?"); ax[1].grid(True, which="both", alpha=0.2); ax[1].legend(fontsize=8.5)
fig.suptitle(f"REAL VCaP — enrichment transfer ê(z) (n_ss_exon={int(ss_exon.sum())})", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.93]); fig.savefig(SCR / "vcap_fig2_ehat.png", dpi=130); plt.close(fig); print("saved vcap_fig2_ehat.png", flush=True)
print("DONE", flush=True)
