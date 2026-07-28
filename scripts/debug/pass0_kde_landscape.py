"""Pass-0 KDE landscape diagnostic — can we detect the gDNA enriched/depleted bimodality at INIT
(no solve), and is the boundary a cleaner window than the region?

Design dialogue: `docs/calibration/archive/node_prior_design.md` §4.4 / Phase C. The thesis (user's):

  * A single-stranded **intron↔exon splice-junction** boundary is a structural gDNA window: its
    UNSPLICED crossing mass EXCLUDES mature RNA (mature is spliced away at the junction), so it is
    nascent+gDNA — and nascent is a small fraction in real biology — i.e. ≈ gDNA even for EXPRESSED loci.
  * A REGION's unspliced mass is mature+nascent+gDNA, so at realistic (low) gDNA it is RNA-DOMINATED and
    its density landscape is confounded (the "enriched mode" is highly-expressed RNA, not captured gDNA).
  * So for UNSTRANDED data the boundary is our best enrichment detector; the boundary SPLICED mass is the
    RNA-level proxy (spliced≈0 ⇒ the adjacent exon is an RNA-free measurement).

This tool MEASURES those claims against oracle truth, across gDNA levels (the sims are high-gDNA = easy;
the hard case is low gDNA where RNA dominates). For each condition it emits, per node:
  * density = mass_global / eff_global (the exact basis the gDNA-density KDE sees),
  * oracle TRUE gDNA fraction of the unspliced mass,
  * node class: single-stranded intron↔exon SJ boundary vs exon region; boundary spliced mass.

Outputs a per-condition PNG (density landscape, coloured by oracle gDNA-fraction) + a summary TSV.
Built on the validated `oracle.OracleTruth` + production `calibrate` (same wiring as calib_pool_benchmark).

    OMP_NUM_THREADS=1 python scripts/debug/pass0_kde_landscape.py --conditions a,b --outdir DIR
"""

from __future__ import annotations

import argparse
import os
import sys
from dataclasses import replace as dc
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from rigel.calibration import calibrate
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.node_chain import REGION
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.signature import (
    BIT_EXON_NEG,
    BIT_EXON_POS,
    nrna_active_strands,
)
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex
from rigel.pipeline import _native_detect_sj_tag, scan_and_buffer
from rigel.splice import SpliceType

sys.path.insert(0, str(Path(__file__).parent))
from oracle import OracleTruth  # noqa: E402

_UNSPL = (0, 1)
_EPS = 1.0e-9


def _sj_boundary_masks(bsub, ra):
    """Per-boundary (length B): single-stranded intron↔exon splice-junction mask + which live strand.
    SJ on strand s ⇔ nrna-active on s on BOTH flanks (continuity) AND exactly ONE flank is exon on s
    (exon XOR intron across the seam). Single-stranded ⇔ exactly one of ±SJ holds."""
    sig = np.asarray(ra.signature).astype(np.int64)
    lr = np.asarray(bsub.left_region, np.int64)
    rr = np.asarray(bsub.right_region, np.int64)
    sig_l = np.where(lr >= 0, sig[np.clip(lr, 0, None)], 0)
    sig_r = np.where(rr >= 0, sig[np.clip(rr, 0, None)], 0)
    nrp_l, nrn_l = nrna_active_strands(sig_l)
    nrp_r, nrn_r = nrna_active_strands(sig_r)
    ex_p_l, ex_p_r = (sig_l & BIT_EXON_POS) != 0, (sig_r & BIT_EXON_POS) != 0
    ex_n_l, ex_n_r = (sig_l & BIT_EXON_NEG) != 0, (sig_r & BIT_EXON_NEG) != 0
    sj_pos = (nrp_l & nrp_r) & (ex_p_l ^ ex_p_r)  # +strand intron↔exon junction
    sj_neg = (nrn_l & nrn_r) & (ex_n_l ^ ex_n_r)
    return sj_pos ^ sj_neg  # single-stranded SJ (exactly one strand); XOR excludes ambiguous double


def analyze(suite: Path, cond: str, index, cfg, work_dir: Path) -> dict:
    """Scan + calibrate + oracle for one condition; return per-node arrays (density, oracle f_g, class)."""
    bam = str(suite / cond / "sim_oracle.bam")
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    sc = dc(cfg.scan, sj_strand_tag=_native_detect_sj_tag(bam))
    _stats, sm, flm, _buf, payload = scan_and_buffer(bam, index, sc)
    fl = build_fl_models(
        global_counts=flm.global_model.counts,
        rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
        gdna_counts=gdna_fl_mass(payload),
        max_size=flm.max_size,
    )
    dbg: dict = {}
    calibrate(
        payload=payload,
        region_arrays=ra,
        strand_model=sm,
        gdna_fl_pmf=fl.gdna_pmf,
        rna_fl_pmf=fl.rna_pmf,
        config=cfg.calibration,
        _debug=dbg,
    )
    chain, geom, bsub = dbg["chain"], dbg["geometry"], dbg["boundary_substrate"]
    kind, ref = np.asarray(chain.kind), np.asarray(chain.ref_idx, np.int64)
    is_reg = kind == REGION
    is_bnd = ~is_reg

    # per-node density on the KDE basis (node_sweep: region = contained; boundary = averaged two-side).
    egl, egr = np.asarray(geom.eff_gdna_left), np.asarray(geom.eff_gdna_right)
    msl, msr = np.asarray(geom.mass_left), np.asarray(geom.mass_right)
    eff_global = np.where(is_reg, egl, 0.5 * (egl + egr))
    mass_global = np.where(is_reg, msl, msl + msr)
    dens = mass_global / np.maximum(eff_global, _EPS)

    # oracle TRUE gDNA fraction of the UNSPLICED mass, per boundary and per region, mapped onto the chain.
    def _bnd_uns(part):
        pl = np.asarray(orc.parts[part].boundary_mass_left, np.float64)[:, _UNSPL].sum(1)
        pr = np.asarray(orc.parts[part].boundary_mass_right, np.float64)[:, _UNSPL].sum(1)
        return pl + pr

    orc = OracleTruth.from_bam(
        bam, index, cfg, work_dir, cond, full_payload=payload, boundary_mass_tol=0.5
    )
    gb, mb, nb = _bnd_uns("gdna"), _bnd_uns("mrna"), _bnd_uns("nrna")
    tot_b = gb + mb + nb
    ofg_b = np.where(tot_b > 0, gb / np.maximum(tot_b, _EPS), np.nan)  # per boundary index
    ofg_r, _tot_r = orc.region_true_fg()  # per region index
    # boundary spliced mass (RNA-level proxy at the junction), per boundary index.
    spl_b = np.asarray(bsub.left.mass_spliced, np.float64) + np.asarray(
        bsub.right.mass_spliced, np.float64
    )

    R = ofg_r.shape[0]
    B = ofg_b.shape[0]
    ri = np.clip(ref, 0, R - 1)
    bi = np.clip(ref, 0, B - 1)
    oracle_fg = np.where(is_reg, ofg_r[ri], ofg_b[bi])
    spliced = np.where(is_bnd, spl_b[bi], 0.0)

    # node classes
    sj_b = _sj_boundary_masks(bsub, ra)  # per boundary index
    is_sj = is_bnd & sj_b[bi]
    sig_chain = np.where(is_reg, np.asarray(ra.signature).astype(np.int64)[ri], 0)
    is_exon_reg = is_reg & ((sig_chain & (BIT_EXON_POS | BIT_EXON_NEG)) != 0)

    active = mass_global > 0.0
    return dict(
        cond=cond,
        dens=dens,
        oracle_fg=oracle_fg,
        spliced=spliced,
        is_sj=is_sj & active,
        is_exon_reg=is_exon_reg & active,
    )


def _panel(ax, x, y, title, xlabel):
    ax.scatter(x, y, s=6, alpha=0.35, c=y, cmap="coolwarm_r", vmin=0, vmax=1)
    ax.set_title(title, fontsize=9)
    ax.set_xlabel(xlabel, fontsize=8)
    ax.set_ylabel("oracle gDNA fraction", fontsize=8)
    ax.set_ylim(-0.05, 1.05)
    ax.axhline(0.5, color="gray", lw=0.5, ls=":")


def _hist(ax, groups, title, xlabel):
    for label, vals, color in groups:
        if len(vals):
            ax.hist(vals, bins=40, alpha=0.5, label=f"{label} (n={len(vals)})", color=color)
    ax.set_title(title, fontsize=9)
    ax.set_xlabel(xlabel, fontsize=8)
    ax.legend(fontsize=7)


def plot_condition(res: dict, outdir: Path):
    ld = np.log10(np.maximum(res["dens"], 1e-6))
    sj, ex = res["is_sj"], res["is_exon_reg"]
    sj_spl0 = sj & (res["spliced"] <= _EPS)
    sj_spl = sj & (res["spliced"] > _EPS)
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    _panel(
        axes[0, 0],
        ld[sj],
        res["oracle_fg"][sj],
        "SJ boundary: log-density vs oracle f_g",
        "log10 density",
    )
    _hist(
        axes[0, 1],
        [("SJ spliced=0", ld[sj_spl0], "tab:blue"), ("SJ spliced>0", ld[sj_spl], "tab:orange")],
        "SJ boundary density (bimodality?)",
        "log10 density",
    )
    _panel(
        axes[1, 0],
        ld[ex],
        res["oracle_fg"][ex],
        "exon region: log-density vs oracle f_g",
        "log10 density",
    )
    _hist(
        axes[1, 1],
        [("exon regions", ld[ex], "tab:green")],
        "exon region density (RNA-confounded?)",
        "log10 density",
    )
    fig.suptitle(res["cond"], fontsize=11)
    fig.tight_layout()
    p = outdir / f"landscape_{res['cond']}.png"
    fig.savefig(p, dpi=110)
    plt.close(fig)
    return p


def summarize(res: dict) -> dict:
    sj, ex = res["is_sj"], res["is_exon_reg"]
    ofg, spl = res["oracle_fg"], res["spliced"]
    sj0 = sj & (spl <= _EPS)

    def med(mask):
        v = ofg[mask]
        v = v[np.isfinite(v)]
        return float(np.median(v)) if len(v) else np.nan

    def puref(mask):  # fraction with oracle gDNA-fraction > 0.8 (gDNA-pure)
        v = ofg[mask]
        v = v[np.isfinite(v)]
        return float((v > 0.8).mean()) if len(v) else np.nan

    return dict(
        condition=res["cond"],
        n_sj=int(sj.sum()),
        n_sj_spliced0=int(sj0.sum()),
        n_exon_reg=int(ex.sum()),
        sj_median_oracle_fg=med(sj),
        sj_spliced0_median_oracle_fg=med(sj0),
        exon_reg_median_oracle_fg=med(ex),
        sj_frac_gdna_pure=puref(sj),
        exon_reg_frac_gdna_pure=puref(ex),
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--suite", default="/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
    ap.add_argument("--conditions", required=True, help="comma-separated condition names")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    suite = Path(args.suite)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    index = TranscriptIndex.load(str(suite / "rigel_index"))
    cfg = PipelineConfig()
    work_dir = Path(os.environ.get("RIGEL_SCRATCH", "/tmp")) / "rigel_pass0_landscape_split"

    rows = []
    for i, cond in enumerate(args.conditions.split(","), 1):
        print(f"[{i}] {cond} ...", flush=True)
        res = analyze(suite, cond, index, cfg, work_dir)
        png = plot_condition(res, outdir)
        rows.append(summarize(res))
        print(f"    -> {png}", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(outdir / "pass0_landscape_summary.tsv", sep="\t", index=False)
    pd.set_option("display.width", 220)
    print("\n=== pass-0 landscape: oracle gDNA-purity of the candidate training nodes ===")
    print(df.to_string(index=False, float_format=lambda x: f"{x:.3f}"))


if __name__ == "__main__":
    main()
