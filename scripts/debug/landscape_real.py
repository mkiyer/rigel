"""Real-data density landscape (no oracle) — visualize what the pass-0 gDNA KDE must model.

Drives off a `calib_cache` pickle (scan-free), runs production `calibrate` with the `_debug` hook to get
the solved chain, and plots the log-density distribution of the **pass-0 confident set** (the user's
inclusion list — nodes whose unspliced mass is gDNA with no solve needed):

    intergenic regions · single-strand intron regions · intergenic↔exon boundaries · SJ (intron↔exon
    single-strand) boundaries

against the exon regions (the RNA-confounded contrast). Colored by the CURRENT production calibrated f_g
(not truth — real STAR BAMs have no oracle). The question this answers: is the gDNA-density landscape a
clean bimodal, a hump-in-a-mountain, or a smear — i.e. what must the projection KDE (KDE × locality
kernel) represent.

    python scripts/debug/landscape_real.py --cache CACHE.pkl --out DIR
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from rigel.calibration.node_chain import REGION
from rigel.calibration.signature import (
    BIT_EXON_NEG,
    BIT_EXON_POS,
    BIT_INTRON_NEG,
    BIT_INTRON_POS,
    TS_NEG,
    TS_POS,
    transcript_strand_class,
)

sys.path.insert(0, str(Path(__file__).parent))
from calib_cache import calibrate_from_cache, load  # noqa: E402
from pass0_kde_landscape import _sj_boundary_masks  # noqa: E402

_EPS = 1.0e-9
_EXON = BIT_EXON_POS | BIT_EXON_NEG
_INTRON = BIT_INTRON_POS | BIT_INTRON_NEG


def analyze(cache_path: str) -> dict:
    inp = load(cache_path)
    dbg: dict = {}
    res = calibrate_from_cache(inp, _debug=dbg)
    chain, geom, bsub, ra = (
        dbg["chain"],
        dbg["geometry"],
        dbg["boundary_substrate"],
        dbg["region_arrays"],
    )
    belief = dbg["belief"]
    kind, ref = np.asarray(chain.kind), np.asarray(chain.ref_idx, np.int64)
    is_reg = kind == REGION
    is_bnd = ~is_reg

    egl, egr = np.asarray(geom.eff_gdna_left), np.asarray(geom.eff_gdna_right)
    msl, msr = np.asarray(geom.mass_left), np.asarray(geom.mass_right)
    eff_global = np.where(is_reg, egl, 0.5 * (egl + egr))
    mass_global = np.where(is_reg, msl, msl + msr)
    dens = mass_global / np.maximum(eff_global, _EPS)
    cal_fg = np.asarray(belief.f_g)

    sig = np.asarray(ra.signature).astype(np.int64)
    sc = transcript_strand_class(sig)
    R = sig.shape[0]
    ri = np.clip(ref, 0, R - 1)
    sig_r = np.where(is_reg, sig[ri], 0)
    sc_r = np.where(is_reg, sc[ri], 0)
    single_reg = is_reg & ((sc_r == TS_POS) | (sc_r == TS_NEG))

    # region-node confident classes
    intergenic_reg = is_reg & (sig_r == 0)
    ss_intron_reg = single_reg & ((sig_r & _INTRON) != 0) & ((sig_r & _EXON) == 0)
    exon_reg = is_reg & ((sig_r & _EXON) != 0)

    # boundary-node confident classes (per boundary index → chain)
    B = np.asarray(bsub.left_region).shape[0]
    bi = np.clip(ref, 0, B - 1)
    lr = np.asarray(bsub.left_region, np.int64)
    rr = np.asarray(bsub.right_region, np.int64)
    sig_l = np.where(lr >= 0, sig[np.clip(lr, 0, None)], 0)
    sig_rr = np.where(rr >= 0, sig[np.clip(rr, 0, None)], 0)
    exon_l, exon_r = (sig_l & _EXON) != 0, (sig_rr & _EXON) != 0
    interg_l, interg_r = sig_l == 0, sig_rr == 0
    intergenic_exon_b = (interg_l & exon_r) | (interg_r & exon_l)  # per boundary
    sj_b = _sj_boundary_masks(bsub, ra)  # per boundary (single-strand intron↔exon)
    spl_b = np.asarray(bsub.left.mass_spliced, float) + np.asarray(bsub.right.mass_spliced, float)

    intergenic_exon_bnd = is_bnd & intergenic_exon_b[bi]
    sj_bnd = is_bnd & sj_b[bi]
    spliced = np.where(is_bnd, spl_b[bi], 0.0)

    active = mass_global > 0.0
    return dict(
        cond=Path(cache_path).stem,
        dens=dens,
        cal_fg=cal_fg,
        spliced=spliced,
        rna_sense_frac=float(res.rna_sense_frac),
        gdna_density_global=float(res.gdna_density_global),
        intergenic_reg=intergenic_reg & active,
        ss_intron_reg=ss_intron_reg & active,
        intergenic_exon_bnd=intergenic_exon_bnd & active,
        sj_bnd=sj_bnd & active,
        exon_reg=exon_reg & active,
    )


def plot(res: dict, outdir: Path) -> Path:
    ld = np.log10(np.maximum(res["dens"], 1e-6))
    sj0 = res["sj_bnd"] & (res["spliced"] <= _EPS)
    confident = [
        ("intergenic reg", res["intergenic_reg"], "tab:gray"),
        ("ss intron reg", res["ss_intron_reg"], "tab:green"),
        ("intergenic↔exon bnd", res["intergenic_exon_bnd"], "tab:purple"),
        ("SJ bnd (spliced=0)", sj0, "tab:blue"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(17, 4.6))
    # A: the pass-0 confident training distribution (what the KDE fits) — modality?
    lo = min(ld[m].min() for _, m, _ in confident if m.any())
    hi = max(ld[m].max() for _, m, _ in confident if m.any())
    bins = np.linspace(lo, hi, 50)
    for label, m, c in confident:
        if m.any():
            axes[0].hist(ld[m], bins=bins, alpha=0.5, label=f"{label} (n={int(m.sum())})", color=c)
    axes[0].set_title("pass-0 confident set: gDNA-density landscape")
    axes[0].set_xlabel("log10 density")
    axes[0].legend(fontsize=7)
    # B: SJ (spliced=0) density vs production calibrated f_g
    axes[1].scatter(
        ld[sj0],
        res["cal_fg"][sj0],
        s=7,
        alpha=0.4,
        c=res["cal_fg"][sj0],
        cmap="coolwarm_r",
        vmin=0,
        vmax=1,
    )
    axes[1].set_title("SJ spliced=0: density vs PRODUCTION f_g")
    axes[1].set_xlabel("log10 density")
    axes[1].set_ylabel("calibrated f_g")
    axes[1].set_ylim(-0.05, 1.05)
    # C: exon regions (RNA-confounded contrast)
    ex = res["exon_reg"]
    if ex.any():
        axes[2].hist(
            ld[ex], bins=50, alpha=0.6, color="tab:orange", label=f"exon reg (n={int(ex.sum())})"
        )
    axes[2].set_title("exon regions (contrast)")
    axes[2].set_xlabel("log10 density")
    axes[2].legend(fontsize=7)
    fig.suptitle(
        f"{res['cond']}   κ(rna_sense_frac)={res['rna_sense_frac']:.3f}   "
        f"gdna_density_global={res['gdna_density_global']:.2e}",
        fontsize=11,
    )
    fig.tight_layout()
    p = outdir / f"real_landscape_{res['cond']}.png"
    fig.savefig(p, dpi=110)
    plt.close(fig)
    return p


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cache", required=True, help="comma-separated calib_cache .pkl paths")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)
    for c in args.cache.split(","):
        res = analyze(c)
        p = plot(res, outdir)
        print(
            f"{res['cond']:10s} κ={res['rna_sense_frac']:.3f} gdna_glob={res['gdna_density_global']:.2e} "
            f"| n: intergenic={int(res['intergenic_reg'].sum())} ss_intron={int(res['ss_intron_reg'].sum())} "
            f"interg↔exon={int(res['intergenic_exon_bnd'].sum())} "
            f"SJ_spl0={int((res['sj_bnd'] & (res['spliced'] <= _EPS)).sum())} -> {p}",
            flush=True,
        )


if __name__ == "__main__":
    main()
