"""Visualize the variance~mean model for the count-module imputation (Phase 4-var).

For every non-observable region with BOTH a left and a right observable boundary anchor, the two
anchor densities (d_L, d_R) are a PAIRED estimate of the SAME node's gDNA density — so their
disagreement is genuine (conflation-free) imputation variance. We collect
(mean=(d_L+d_R)/2, var=¼(d_L−d_R)²) across all such nodes and fit variance ~ mean, which denoises the
noisy per-node 2-point variance into a smooth curve. Plots capture-on vs capture-off to show that the
imputation variance grows under capture (the regime that needs it).
"""
import os, numpy as np
from dataclasses import replace as _dc_replace
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from rigel.index import TranscriptIndex
from rigel.config import PipelineConfig
from rigel.pipeline import scan_and_buffer, _native_detect_sj_tag
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.calibration.density_model import count_observable_masks
from rigel.calibration.effective_length import region_eff_length, boundary_eff_length
from rigel.splice import SpliceType

SUITE = "/Users/mkiyer/Downloads/rigel_runs/gdna_benchmark_5mb"
idx = TranscriptIndex.load(f"{SUITE}/rigel_index")


def anchor_pairs(cond):
    """Return (mean, var) for non-observable regions with both L and R observable anchors."""
    bam = f"{SUITE}/{cond}/sim_oracle.bam"
    cfg = PipelineConfig(); scan = _dc_replace(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    stats, sm, flm, buffer, payload = scan_and_buffer(bam, idx, scan); sm.finalize()
    ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
    fl = build_fl_models(global_counts=flm.global_model.counts,
                         rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                         gdna_counts=gdna_fl_mass(payload), max_size=flm.max_size)
    sub = CalibrationSubstrate.from_payload(payload, ra)
    fl_mean = boundary_eff_length(fl.gdna_pmf); inv_fl = 1.0 / fl_mean
    sig = np.asarray(ra.signature); ref_id = np.asarray(ra.ref_id)
    region_obs, bnd_obs = count_observable_masks(sig, ref_id)
    r = sig.shape[0]
    left = (np.asarray(sub.left.n_unspliced_pos) + np.asarray(sub.left.n_unspliced_neg)).astype(float)
    right = (np.asarray(sub.right.n_unspliced_pos) + np.asarray(sub.right.n_unspliced_neg)).astype(float)
    l_anchor = np.zeros(r, bool); r_anchor = np.zeros(r, bool)
    if r > 1:
        l_anchor[1:] = bnd_obs[:-1] & (ref_id[1:] == ref_id[:-1])
        r_anchor[:-1] = bnd_obs[:-1] & (ref_id[:-1] == ref_id[1:])
    sel = (~region_obs) & l_anchor & r_anchor & (left > 0) & (right > 0)
    d_L = left[sel] * inv_fl; d_R = right[sel] * inv_fl
    mean = 0.5 * (d_L + d_R); var = 0.25 * (d_L - d_R) ** 2
    return mean, var


fig, axes = plt.subplots(1, 2, figsize=(13, 5.2), sharex=False)
for ax, cond, title in [
    (axes[0], "gdna_gdna1000_ss_0.99_nrna_none_capture_off", "capture OFF (gdna1000)"),
    (axes[1], "gdna_gdna1000_ss_0.99_nrna_none_capture_on", "capture ON (gdna1000)"),
]:
    m, v = anchor_pairs(cond)
    ax.scatter(m, v, s=6, alpha=0.25, color="#4477aa", label=f"{len(m)} 2-anchor nodes")
    # binned-mean trend (the variance~mean fit), in log-mean bins
    order = np.argsort(m); ms, vs = m[order], v[order]
    pos = ms > 0
    bins = np.geomspace(max(ms[pos].min(), 1e-3), ms.max() + 1e-9, 12)
    bidx = np.digitize(ms, bins)
    bx, by = [], []
    for b in range(1, len(bins)):
        msk = bidx == b
        if msk.sum() >= 5:
            bx.append(ms[msk].mean()); by.append(vs[msk].mean())
    ax.plot(bx, by, "-o", color="#cc3311", lw=2, label="binned variance~mean")
    # reference: Poisson-equivalent var of the density (mean·inv_fl) — shows NB >> Poisson under capture
    ax.set_title(title); ax.set_xlabel("mean anchor density  (½(d_L+d_R))")
    ax.set_ylabel("imputation variance  ¼(d_L−d_R)²"); ax.legend(loc="upper left", fontsize=8)
    ax.set_xscale("log"); ax.set_yscale("log")
    # quick fit var = a·mean^b (log-log slope) over binned points
    if len(bx) >= 3:
        bxa, bya = np.log(np.array(bx)), np.log(np.maximum(np.array(by), 1e-12))
        b_slope, b_int = np.polyfit(bxa, bya, 1)
        ax.text(0.05, 0.05, f"binned fit: var ~ mean^{b_slope:.2f}", transform=ax.transAxes,
                fontsize=9, color="#cc3311")
        print(f"{title}: n={len(m)}  var~mean^{b_slope:.2f}  "
              f"median mean={np.median(m):.3g} median var={np.median(v):.3g}")
fig.suptitle("Phase 4-var: variance~mean from 2-anchor boundary disagreements (paired, conflation-free)")
fig.tight_layout()
out = "/tmp/variance_mean.png"; fig.savefig(out, dpi=110)
print(f"saved {out}")
