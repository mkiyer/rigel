"""var~mean span intuition: gDNA + RNA disagreement points with LOESS at several spans + a power-law.

For both the gDNA (triplet boundary/contained disagreement) and the RNA (splice-junction-pair) var~mean,
overlay the raw (μ, raw_var) points with the robust LOESS at spans {0.3, 0.5, 0.75, 0.9} and a
parameter-free global log-log power-law (var = α·μ^β). Use this to choose a span that is smooth/monotonic,
or to decide the power-law (no span) is enough.

Usage: python scripts/debug/plot_var_mean_span.py [condition] [out.png]
"""
import dataclasses
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from rigel.calibration.density_model import _loess, count_observable_masks
from rigel.calibration.effective_length import boundary_eff_length, region_eff_length
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.rna_variance import rna_spliced_variance
from rigel.calibration.run_fill import same_ref_left_right
from rigel.calibration.substrate import CalibrationSubstrate
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex
from rigel.pipeline import _native_detect_sj_tag, scan_and_buffer
from rigel.splice import SpliceType

SUITE = "/Users/mkiyer/Downloads/rigel_runs/quick_3to1_5mb"
SPANS = [0.3, 0.5, 0.75, 0.9]
SPAN_C = ["tab:blue", "tab:green", "tab:orange", "tab:red"]


def gdna_points(sub, ra, reg_el, fl_mean):
    sig = np.asarray(ra.signature)
    ref_id = np.asarray(ra.ref_id)
    R = sig.shape[0]
    L = np.maximum(reg_el, 1e-9)
    inv_fl = 1.0 / fl_mean

    def total(v):
        return v.n_unspliced_pos.astype(float) + v.n_unspliced_neg.astype(float)

    rco, bco = count_observable_masks(sig, ref_id)
    ls, rs = same_ref_left_right(ref_id)
    la = np.zeros(R, bool)
    rb = np.zeros(R, bool)
    if R > 1:
        la[1:] = bco[:-1] & ls[1:]
        rb[:-1] = bco[:-1] & rs[:-1]
    d_left = np.where(la, total(sub.left) * inv_fl, np.nan)
    d_right = np.where(rb, total(sub.right) * inv_fl, np.nan)
    c_ok = rco & (L > 1e-9)
    contained = np.where(c_ok, total(sub.contained) / L, np.nan)
    obs = np.stack([d_left, d_right, contained])
    msk = np.stack([la, rb, c_ok])
    k = msk.sum(0).astype(float)
    mu = np.where(msk, np.nan_to_num(obs), 0.0).sum(0) / np.maximum(k, 1.0)
    dev2 = np.where(msk, (np.nan_to_num(obs) - mu) ** 2, 0.0).sum(0)
    raw = np.where(k > 1.0, dev2 / np.maximum(k - 1.0, 1e-9), np.nan) / np.maximum(k, 1.0)
    sel = (k >= 2.0) & np.isfinite(mu) & (mu > 1e-9) & np.isfinite(raw) & (raw > 1e-9)
    return mu[sel], raw[sel]


def panel(ax, mu, raw, title):
    if mu.size < 5:
        ax.set_title(f"{title} (only {mu.size} points)")
        return
    x, y = np.log(mu), np.log(raw)
    ax.scatter(mu, raw, s=10, c="0.6", alpha=0.5, label=f"{mu.size} pts")
    xq = np.logspace(np.log10(mu.min()), np.log10(mu.max()), 100)
    for span, col in zip(SPANS, SPAN_C):
        yq = np.exp(_loess(x, y, np.log(xq), span))
        ax.plot(xq, yq, color=col, lw=1.6, label=f"LOESS span={span}")
    b, a = np.polyfit(x, y, 1)  # power-law var = e^a · μ^b
    ax.plot(xq, np.exp(a) * xq**b, "k--", lw=2, label=f"power-law β={b:.2f}")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("mean density μ")
    ax.set_ylabel("raw_var")
    ax.set_title(title)
    ax.legend(fontsize=7)


def main():
    cond = sys.argv[1] if len(sys.argv) > 1 else "gdna_gdna300_ss_0.99_nrna_none_capture_on"
    out = sys.argv[2] if len(sys.argv) > 2 else f"{SUITE}/{cond}/var_mean_span.png"
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
    reg_el_rna = region_eff_length(ra.region_size_bp, fl.rna_pmf)
    fl_mean_rna = boundary_eff_length(fl.rna_pmf)

    g_mu, g_raw = gdna_points(sub, ra, reg_el, fl_mean)
    rv = rna_spliced_variance(sub, ra, reg_el_rna, fl_mean_rna)
    r_ok = np.isfinite(rv.fit_mu) & (rv.fit_mu > 1e-9) & np.isfinite(rv.fit_var) & (rv.fit_var > 1e-9)

    fig, ax = plt.subplots(1, 2, figsize=(15, 6))
    fig.suptitle(f"var~mean span sensitivity — {cond}", fontsize=11)
    panel(ax[0], g_mu, g_raw, f"gDNA (triplet disagreement, {g_mu.size} pts)")
    panel(ax[1], rv.fit_mu[r_ok], rv.fit_var[r_ok], f"RNA (splice-junction pairs, {int(r_ok.sum())} pts)")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out, dpi=120)
    print(f"wrote {out}")
    print(f"gDNA: {g_mu.size} pts, power-law β={np.polyfit(np.log(g_mu),np.log(g_raw),1)[0]:.2f}")
    print(f"RNA:  {int(r_ok.sum())} pts, power-law β="
          f"{np.polyfit(np.log(rv.fit_mu[r_ok]),np.log(rv.fit_var[r_ok]),1)[0]:.2f}")


if __name__ == "__main__":
    main()
