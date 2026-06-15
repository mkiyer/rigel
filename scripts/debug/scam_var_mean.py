"""Var~mean fitter iteration surface — the production monotone P-spline (variance_model.MonotoneVarMean)
on real pipeline data, with the DIRECT vs IMPUTATION split and a comparison to LOESS / isotonic / power-law.

This is the surface for refining the fitter (K, robustness, locality/adaptive span). It calls the SAME
``variance_model`` the calibration uses, so a change there is reflected here. Validates the §17 claim that
the IMPUTATION curve sits ABOVE the DIRECT curve (imputed nodes are properly humbler).

Usage: python scripts/debug/scam_var_mean.py [condition] [out.png]
"""
import dataclasses
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from rigel.calibration.density_model import _loess
from rigel.calibration.effective_length import boundary_eff_length, region_eff_length
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.calibration.variance_model import (
    MonotoneVarMean, fit_direct_varmean, fit_imputation_varmean, varmean_points,
)
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex
from rigel.pipeline import _native_detect_sj_tag, scan_and_buffer
from rigel.splice import SpliceType

SUITE = "/Users/mkiyer/Downloads/rigel_runs/quick_3to1_5mb"
cond = sys.argv[1] if len(sys.argv) > 1 else "gdna_gdna300_ss_0.99_nrna_none_capture_on"
out = sys.argv[2] if len(sys.argv) > 2 else "/tmp/scam_var_mean.png"
bam = f"{SUITE}/{cond}/sim_oracle.bam"
idx = TranscriptIndex.load(f"{SUITE}/rigel_index")
cfg = PipelineConfig()
scan = dataclasses.replace(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))
_st, sm, flm, _buf, pl = scan_and_buffer(bam, idx, scan)
ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
fl = build_fl_models(global_counts=flm.global_model.counts,
                     rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                     gdna_counts=gdna_fl_mass(pl), max_size=flm.max_size)
sub = CalibrationSubstrate.from_payload(pl, ra)
reg_el = region_eff_length(ra.region_size_bp, fl.gdna_pmf)
fl_mean = boundary_eff_length(fl.gdna_pmf)

pts = varmean_points(sub, ra, reg_el, fl_mean)
n_dir = int(pts.region_observable.sum())
n_imp = int((~pts.region_observable).sum())
direct = fit_direct_varmean(pts)
imp = fit_imputation_varmean(pts)
allfit = MonotoneVarMean.fit(pts.mean, pts.raw_var)
print(f"=== var~mean points: {pts.mean.size} ({n_dir} DIRECT region-observable, {n_imp} IMPUTATION) ===")
print(f"DIRECT     : λ={direct.lam:.3g} edf={direct.edf:.2f}")
print(f"IMPUTATION : λ={imp.lam:.3g} edf={imp.edf:.2f}")
# the key §17 check: imputation variance >= direct variance (imputed nodes humbler)
g = np.logspace(np.log10(max(pts.mean.min(), 1e-3)), np.log10(pts.mean.max()), 50)
ratio = imp.predict(g) / np.maximum(direct.predict(g), 1e-12)
print(f"IMPUTATION/DIRECT variance ratio over the range: median {np.median(ratio):.2f} "
      f"(>1 ⇒ imputation humbler, as intended)")

fig, ax = plt.subplots(1, 2, figsize=(15, 6))
dm, im = pts.region_observable, ~pts.region_observable
for a in ax:
    a.scatter(pts.mean[dm], pts.raw_var[dm], s=10, c="tab:green", alpha=0.4, label=f"DIRECT pts ({n_dir})")
    a.scatter(pts.mean[im], pts.raw_var[im], s=10, c="tab:red", alpha=0.4, label=f"IMPUTATION pts ({n_imp})")
    a.set_xscale("log")
    a.set_yscale("log")
    a.set_xlabel("mean gDNA density μ")
    a.set_ylabel("raw var")
gg = np.logspace(np.log10(pts.mean.min()), np.log10(pts.mean.max()), 200)
ax[0].plot(gg, direct.predict(gg), "tab:green", lw=2.5, label="DIRECT mono-Pspline")
ax[0].plot(gg, imp.predict(gg), "tab:red", lw=2.5, label="IMPUTATION mono-Pspline")
ax[0].set_title("DIRECT vs IMPUTATION (production variance_model)")
ax[0].legend(fontsize=8, loc="upper left")
# comparison panel: the combined fit vs LOESS / power-law
x, y = np.log(pts.mean), np.log(pts.raw_var)
ax[1].plot(gg, allfit.predict(gg), "b-", lw=2.5, label=f"SCAM all (edf={allfit.edf:.1f})")
ax[1].plot(gg, np.exp(_loess(x, y, np.log(gg), 0.4)), color="0.3", lw=1.3, label="LOESS 0.4")
b, a0 = np.polyfit(x, y, 1)
ax[1].plot(gg, np.exp(a0) * gg**b, "k--", lw=1.6, label=f"power-law β={b:.2f}")
ax[1].set_title("SCAM vs LOESS / power-law")
ax[1].legend(fontsize=8, loc="upper left")
fig.suptitle(f"var~mean (monotone P-spline) — {cond}", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(out, dpi=120, bbox_inches="tight")
print(f"wrote {out}")
