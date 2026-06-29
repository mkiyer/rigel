"""CAUSAL test: does dropping ONLY the monotone-increasing constraint on σ²_bio (a free GCV-penalized
P-spline that CAN fit the hump) reduce the enriched-exon under-call on the flagship? Minimal change —
same setup/GCV-λ/offset/robustness, only _fit_monotone → unconstrained penalized WLS inside the σ²_bio fit.
"""
import sys
from pathlib import Path
import numpy as np
from scipy.interpolate import BSpline

import rigel.calibration.variance_model as vm
from rigel.calibration.variance_model import MonotoneVarMean, _setup_spline, _select_lambda
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
COND = sys.argv[1] if len(sys.argv) > 1 else "gdna_gdna300_ss_0.99_nrna_none_capture_on"
BAM = ROOT / COND / "sim_oracle.bam"
print(f"### COND={COND}")

# ---- free (non-monotone) penalized-WLS variant of MonotoneVarMean.fit_offset (isolates the σ²_bio fit) ----
_orig_fit_offset = MonotoneVarMean.fit_offset.__func__


def _free_fit_offset(cls, mean, raw, offset, weight=None, *, k=18, degree=3,
                     robust_iters=2, n_lambda=40, population_spread=False):
    mean = np.asarray(mean, float); raw = np.asarray(raw, float); off = np.asarray(offset, float)
    if population_spread:
        wt = np.ones_like(off); robust_iters = 0
    else:
        wt = np.ones_like(off) if weight is None else np.asarray(weight, float)
    ok = (np.isfinite(mean) & (mean > _EPS) & np.isfinite(raw) & (raw >= 0.0)
          & np.isfinite(off) & (off >= 0.0) & np.isfinite(wt) & (wt > 0.0))
    mean, raw, off, wt = mean[ok], raw[ok], off[ok], wt[ok]
    if mean.size < max(k, 8):
        return cls._constant_offset(mean, raw, off, degree, wt)
    order = np.argsort(mean); mean, raw, off, wt = mean[order], raw[order], off[order], wt[order]
    x = np.log(mean); z = raw - off
    if float(x[-1]) - float(x[0]) < _EPS:
        return cls._constant_offset(mean, raw, off, degree, wt)
    kn, B, P, x_lo, x_hi = _setup_spline(x, k, degree)
    lam, edf = _select_lambda(B, z, P, mean.size, n_lambda)

    def free_solve(Bm, y, w, Pm, lm):
        A = Bm.T @ (w[:, None] * Bm) + lm * Pm
        return np.linalg.solve(A, Bm.T @ (w * y))   # UNCONSTRAINED penalized WLS (no monotone cumsum-exp)

    total = np.maximum(off, _EPS); coeffs = None
    for it in range(robust_iters + 1):
        w = wt if population_spread else wt / np.maximum(total, _EPS) ** 2
        if it > 0:
            r = (z - B @ coeffs) * np.sqrt(w)
            s_mad = 1.4826 * float(np.median(np.abs(r - np.median(r)))) + _EPS
            u = r / (4.685 * s_mad)
            w = w * np.where(np.abs(u) < 1.0, (1.0 - u ** 2) ** 2, 0.0)
        coeffs = free_solve(B, z, w, P, lam)
        if it == robust_iters:
            break
        total = np.maximum(B @ coeffs, 0.0) + off
    return cls(knots=kn, degree=degree, coeffs=coeffs, x_lo=x_lo, x_hi=x_hi,
               fit_mean=mean, fit_var=np.maximum(z, 0.0), edf=float(edf), lam=float(lam))


def setup():
    idx = TranscriptIndex.load(str(ROOT / "rigel_index"))
    st, sm, flm, buf, pl = scan_and_buffer(str(BAM), idx, BamScanConfig())
    fl = build_fl_models(global_counts=flm.global_model.counts,
                         rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                         gdna_counts=gdna_fl_mass(pl), max_size=flm.max_size)
    ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
    sub = CalibrationSubstrate.from_payload(pl, ra); bsub = BoundarySubstrate.from_payload(pl)
    chain = build_node_chain(pl.ref_region_offsets, pl.ref_boundary_offsets)
    geom = bp.build_node_geometry(chain, sub, bsub, ra, fl.gdna_pmf, fl.rna_pmf)
    statics = bp.build_node_statics(chain, sub, bsub, ra)
    kappa = float(fit_strand_balance(sm).rna_sense_frac)
    return idx, ra, sub, bsub, chain, geom, statics, kappa


def run(idx, ra, sub, bsub, chain, geom, statics, kappa):
    b0 = bp.init_beliefs(chain, sub, bsub, ra, rna_sense_frac=kappa, n_grid=60, statics=statics)
    b, _ = bp.node_sweep(chain, statics, geom, b0, ra, bsub, rna_sense_frac=kappa,
                         n_grid=60, max_outer=CalibrationConfig().sweep_max_outer)
    dens = bp.node_densities(b, geom)
    vmfit = bp.fit_gdna_varmean(chain, dens, geom, statics)
    return b, dens, vmfit


S = setup()
idx = S[0]; ra = S[1]; sub = S[2]; chain = S[4]
kind = np.asarray(chain.kind); refidx = np.asarray(chain.ref_idx, np.int64)
is_reg = kind == REGION; R = len(idx.region_df)
sig = np.asarray(ra.signature)
tcl = np.array([coarse_type_from_signature(int(s)) for s in sig])
ntype = np.where(is_reg, tcl[np.clip(refidx, 0, R - 1)], -1)

print("=== MONOTONE (production) ===")
MonotoneVarMean.fit_offset = classmethod(_orig_fit_offset)
b_m, dens_m, vm_m = run(*S)
print("=== FREE (non-monotone) ===")
MonotoneVarMean.fit_offset = classmethod(_free_fit_offset)
b_f, dens_f, vm_f = run(*S)
MonotoneVarMean.fit_offset = classmethod(_orig_fit_offset)

# (1) does FREE fit the hump? evaluate at transition μ and high μ.
for mu_q, nm in ((1.5, "transition μ≈1.5 (true raw~30)"), (20.0, "high μ≈20 (true raw~2.2)")):
    print(f"  σ²_bio @ {nm}: monotone={vm_m.predict(np.array([mu_q]))[0]:.2f}  free={vm_f.predict(np.array([mu_q]))[0]:.2f}")

# (2)+(3) converged f_g on ENRICHED exon region nodes (high gDNA density) — the under-call locus.
fg_m = np.asarray(b_m.f_g); fg_f = np.asarray(b_f.f_g)
rho_m = 0.5 * (np.asarray(dens_m.rho_g_left) + np.asarray(dens_m.rho_g_right))
exon = is_reg & (ntype == 2)
enr = exon & (rho_m > 5.0)   # enriched exon regions (ρ_g > 5 = probe-targeted)
print(f"\nENRICHED exon region nodes (ρ_g>5): n={int(enr.sum())}")
print(f"  mean f_g: monotone={fg_m[enr].mean():.3f}  free={fg_f[enr].mean():.3f}  (under-call ⇒ higher=better)")
print(f"  median f_g: monotone={np.median(fg_m[enr]):.3f}  free={np.median(fg_f[enr]):.3f}")
# aggregate: total calibrated gDNA mass over region nodes (production under-calls vs the 3:1 library)
mug = 0.5 * (np.asarray(S[5].mass_left) + np.asarray(S[5].mass_right))  # S[5]=geom; per-node unspliced mass
gtot_m = float(np.sum(fg_m[is_reg] * mug[is_reg])); rtot_m = float(np.sum((1 - fg_m[is_reg]) * mug[is_reg]))
gtot_f = float(np.sum(fg_f[is_reg] * mug[is_reg])); rtot_f = float(np.sum((1 - fg_f[is_reg]) * mug[is_reg]))
print(f"\n  library gDNA frac over region unspliced mass: monotone={gtot_m/max(gtot_m+rtot_m,1e-9):.3f}  "
      f"free={gtot_f/max(gtot_f+rtot_f,1e-9):.3f}  (production under-calls; higher = toward truth)")
print(f"  σ²_bio max over μ-grid: monotone={vm_m.predict(np.array([0.5,1.0,1.5,2.0,5.0,20.0])).max():.2f}  "
      f"free={vm_f.predict(np.array([0.5,1.0,1.5,2.0,5.0,20.0])).max():.2f}")
