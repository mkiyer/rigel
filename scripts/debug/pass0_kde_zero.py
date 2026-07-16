"""Count–length view of the confident gDNA set — the discreteness + the zero problem, one root cause.

Density ρ = count/E. Both are discrete/bounded: at low integer count over a SHORT boundary eff-length
(~fragment length), ρ jumps in visible steps, and ρ can NEVER be 0 (a 0 count floors to the KDE's 1/E).
This tool plots the confident set in (log E, log count) space (INCLUDING count=0 nodes, which the density
substrate drops) to show: (a) the discrete count bands, (b) boundaries (short E) vs regions (long E),
(c) how many confident nodes are truly count=0 — the zero-inflation π and the reason a density-only KDE
cannot represent absence.

    python scripts/debug/pass0_kde_zero.py --cache C.pkl[,...] --out DIR
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
from calib_cache import load  # noqa: E402
from pass0_kde_landscape import _sj_boundary_masks  # noqa: E402
from pass0_kde_prototype import _build_structures  # noqa: E402

_EPS = 1.0e-9
_EXON = BIT_EXON_POS | BIT_EXON_NEG
_INTRON = BIT_INTRON_POS | BIT_INTRON_NEG


def confident(inp):
    """Per confident node: integer unspliced count, gDNA eff-length E, is_boundary — INCLUDING count=0."""
    chain, geom, st, _sub, bsub, ra = _build_structures(inp)
    kind, ref = np.asarray(chain.kind), np.asarray(chain.ref_idx, np.int64)
    is_reg = kind == REGION
    is_bnd = ~is_reg
    egl, egr = np.asarray(geom.eff_gdna_left), np.asarray(geom.eff_gdna_right)
    eff = np.where(is_reg, egl, 0.5 * (egl + egr))
    count = np.asarray(st.u_pos, float) + np.asarray(st.u_neg, float)  # integer unspliced count

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
    sj0_b = _sj_boundary_masks(bsub, ra) & (spl_b <= _EPS)
    interg_exon_bnd = is_bnd & interg_exon_b[bi]
    sj0_bnd = is_bnd & sj0_b[bi]

    m = intergenic_reg | ss_intron_reg | interg_exon_bnd | sj0_bnd
    return count[m], np.maximum(eff[m], _EPS), is_bnd[m]


def run(cache_path, outdir: Path):
    cond = Path(cache_path).stem
    count, eff, is_bnd = confident(load(cache_path))
    n = count.shape[0]
    zero = count <= 0.0
    pi_zero = float(zero.mean())
    nz = ~zero
    dens = count[nz] / eff[nz]
    ld = np.log10(dens)
    lE = np.log10(eff[nz])
    lc = np.log10(count[nz])
    bnd = is_bnd[nz]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    # A: (log E, log count) — the discreteness + boundary/region split. Horizontal bands = integer counts.
    axes[0].scatter(lE[~bnd], lc[~bnd], s=5, alpha=0.25, color="tab:green", label="region")
    axes[0].scatter(lE[bnd], lc[bnd], s=5, alpha=0.25, color="tab:blue", label="boundary")
    for k in (1, 2, 3, 5, 10):
        axes[0].axhline(np.log10(k), color="gray", lw=0.4, ls=":")
    axes[0].set_title(f"{cond}: confident nodes in (E, count) space")
    axes[0].set_xlabel("log10 gDNA eff-length E")
    axes[0].set_ylabel("log10 count")
    axes[0].legend(fontsize=8)
    # B: fine-binned density — discrete atoms (k/E stripes) at low count / short E (boundaries).
    axes[1].hist(ld[~bnd], bins=200, alpha=0.5, color="tab:green", label="region")
    axes[1].hist(ld[bnd], bins=200, alpha=0.5, color="tab:blue", label="boundary")
    axes[1].set_title("density = count/E (fine bins — discrete atoms?)")
    axes[1].set_xlabel("log10 density")
    axes[1].legend(fontsize=8)
    fig.suptitle(
        f"{cond}   n_confident={n}   count=0 (π_zero)={pi_zero:.3f}   "
        f"count=1: {float((count == 1).mean()):.3f}  count≤3: {float((count <= 3).mean()):.3f}",
        fontsize=11,
    )
    fig.tight_layout()
    p = outdir / f"count_length_{cond}.png"
    fig.savefig(p, dpi=110)
    plt.close(fig)
    med_bE = np.median(eff[is_bnd]) if is_bnd.any() else np.nan
    med_rE = np.median(eff[~is_bnd]) if (~is_bnd).any() else np.nan
    print(
        f"{cond:10s} n={n} π_zero={pi_zero:.3f} count=1:{float((count == 1).mean()):.3f} "
        f"count≤3:{float((count <= 3).mean()):.3f} medE(bnd)={med_bE:.0f} medE(reg)={med_rE:.0f} -> {p}",
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
