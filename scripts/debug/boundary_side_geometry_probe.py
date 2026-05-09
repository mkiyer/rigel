"""Boundary-flux diagnostic for the rho_ex / rho_ig deficit.

Re-runs the BAM scan against the v3 region partition for ONE simulation
condition, then reports three orthogonal slices of the EXON-INTRON
density estimator:

  (A) rho_ex split by side eligibility:
        - internal:        bf_left & bf_right  (both sides eligible)
        - terminal:        bf_left ^ bf_right  (only one side eligible)
      A persistent gap here would rescue the "terminal exons" hypothesis.

  (B) per-EXON observed/expected ratio binned by exon width.
      Concentration on small exons => double-cross interaction.
      Flat across widths => aligner / FL bias.

  (C) FL distribution of mask=EXON|INTRON (the boundary-crossers) vs
      mask=INTRON-only (a near-pure gDNA population at high contam).
      Length-distribution mismatch => systematic soft-clip bias near
      annotated boundaries.

Run:
    conda activate rigel && python scripts/debug/boundary_side_geometry_probe.py \
        [--bam PATH] [--index PATH] [--K 0|3] [--out scratch/probe.txt]

Defaults to the highest-contamination synthetic condition under
/Users/mkiyer/Downloads/rigel_runs/sim_synthetic.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

from rigel.config import BamScanConfig
from rigel.index import TranscriptIndex
from rigel.calibration.scan_payload import (
    MASK_EXON, MASK_INTRON, MASK_INTERGENIC,
)
from rigel.pipeline import scan_and_buffer


SIM_ROOT = Path("/Users/mkiyer/Downloads/rigel_runs/sim_synthetic")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--bam",
                   default=str(SIM_ROOT / "gdna_high_ss_0.99_nrna_none" / "annotated.bam"))
    p.add_argument("--index", default=str(SIM_ROOT / "rigel_index"))
    p.add_argument("--K", type=int, default=3,
                   help="boundary tolerance K to test (default 3)")
    p.add_argument("--out", default=None,
                   help="optional path to dump full text report")
    return p.parse_args()


def fl_summary(hist: np.ndarray) -> tuple[np.ndarray, float, int]:
    """Return (pmf, mean, n) from a 1-D count histogram."""
    n = int(hist.sum())
    if n == 0:
        return np.zeros_like(hist, dtype=np.float64), 0.0, 0
    pmf = hist.astype(np.float64) / float(n)
    mean = float((np.arange(len(hist), dtype=np.float64) * pmf).sum())
    return pmf, mean, n


def b_cross_from_pmf(pmf: np.ndarray, K: int) -> float:
    """Σ_ℓ pmf[ℓ] · max(ℓ - 2*max(K,1) + 1, 0)."""
    q = max(K, 1)
    ells = np.arange(len(pmf), dtype=np.float64)
    width = np.maximum(ells - 2.0 * q + 1.0, 0.0)
    return float((pmf * width).sum())


def l_eff_contained(spans: np.ndarray, pmf: np.ndarray) -> np.ndarray:
    """FL-PMF-weighted contained effective length per region.

    L_eff(s) = Σ_{ell <= s} pmf[ell] * (s - ell + 1)
    """
    max_ell = len(pmf)
    ells = np.arange(max_ell, dtype=np.float64)
    out = np.zeros(len(spans), dtype=np.float64)
    for i, s in enumerate(spans):
        sf = float(s)
        valid = ells <= sf
        out[i] = float((pmf[valid] * (sf - ells[valid] + 1.0)).sum())
    return out


def section(title: str) -> str:
    return f"\n{'=' * 78}\n{title}\n{'=' * 78}"


def main() -> int:
    args = parse_args()
    out_lines: list[str] = []

    def log(msg: str = "") -> None:
        out_lines.append(msg)
        print(msg)

    log(f"[probe] BAM   : {args.bam}")
    log(f"[probe] INDEX : {args.index}")
    log(f"[probe] K     : {args.K}")

    log("\n[probe] loading index ...")
    index = TranscriptIndex.load(args.index)
    region_df = index.region_df.reset_index(drop=True)
    log(f"[probe] {len(region_df):,} regions in partition")

    log("[probe] running scan_and_buffer ...")
    scan_cfg = BamScanConfig(splicing_anchor_tolerance=args.K)
    _stats, _strand, _fl_models, _buf, payload = scan_and_buffer(
        args.bam, index, scan_cfg,
    )
    if payload is None:
        log("ERROR: payload is None; index lacks regions.feather")
        return 1
    log(f"[probe] payload n_observed = {payload.n_observed:,}")

    # --------------------------------------------------------------
    # Global mask population
    # --------------------------------------------------------------
    log(section("GLOBAL MASK POPULATION"))
    mask_names = {
        0: "(no annotation)",
        MASK_EXON: "EXON",
        MASK_INTRON: "INTRON",
        MASK_EXON | MASK_INTRON: "EXON|INTRON  (boundary-crossing)",
        MASK_INTERGENIC: "INTERGENIC",
        MASK_EXON | MASK_INTERGENIC: "EXON|INTERGENIC",
        MASK_INTRON | MASK_INTERGENIC: "INTRON|INTERGENIC",
        7: "EXON|INTRON|INTERGENIC",
    }
    total = int(payload.global_counts.sum())
    for m in range(8):
        c = int(payload.global_counts[m])
        if c == 0:
            continue
        log(f"  mask={m}  {mask_names[m]:<40}  n={c:>10,}  ({100*c/total:5.2f}%)")

    # --------------------------------------------------------------
    # FL distributions (Diagnostic C)
    # --------------------------------------------------------------
    log(section("(C) FRAGMENT-LENGTH DISTRIBUTION COMPARISON"))
    pmf_xi, mean_xi, n_xi = fl_summary(payload.fl_hist[MASK_EXON | MASK_INTRON])
    pmf_in, mean_in, n_in = fl_summary(payload.fl_hist[MASK_INTRON])
    pmf_ig, mean_ig, n_ig = fl_summary(payload.fl_hist[MASK_INTERGENIC])
    log(f"  mask=EXON|INTRON  (boundary)  : n={n_xi:>10,}  mean_FL={mean_xi:7.2f}")
    log(f"  mask=INTRON       (intron)    : n={n_in:>10,}  mean_FL={mean_in:7.2f}")
    log(f"  mask=INTERGENIC   (intergenic): n={n_ig:>10,}  mean_FL={mean_ig:7.2f}")

    if n_xi > 0 and (n_in > 0 or n_ig > 0):
        def qs(pmf):
            cdf = np.cumsum(pmf)
            return [int(np.searchsorted(cdf, q)) for q in (0.10, 0.25, 0.50, 0.75, 0.90)]
        log(f"  quantiles (10/25/50/75/90 bp):")
        log(f"     EXON|INTRON  : {qs(pmf_xi)}")
        log(f"     INTRON       : {qs(pmf_in)}")
        log(f"     INTERGENIC   : {qs(pmf_ig)}")
        if n_in > 0:
            log(f"  mean ratio  EXON|INTRON / INTRON     = "
                f"{mean_xi / max(mean_in, 1e-9):.4f}")
        if n_ig > 0:
            log(f"  mean ratio  EXON|INTRON / INTERGENIC = "
                f"{mean_xi / max(mean_ig, 1e-9):.4f}")
        log(f"  >>> if << 1: boundary-crossers are systematically SHORTER")
        log(f"      than interior gDNA fragments — consistent with aligner")
        log(f"      soft-clipping near exon-intron edges.")

    # gDNA FL proxy for B_cross
    if n_ig >= 5000:
        gdna_pmf = pmf_ig
        gdna_label = "INTERGENIC"
    elif n_in >= 5000:
        gdna_pmf = pmf_in
        gdna_label = "INTRON"
    else:
        log("WARN: no large-enough pure-gDNA pool; using INTERGENIC PMF "
            "(may be noisy)")
        gdna_pmf = pmf_ig
        gdna_label = "INTERGENIC"
    b_cross = b_cross_from_pmf(gdna_pmf, args.K)
    log(f"\n  gDNA FL proxy = {gdna_label}; B_cross(K={args.K}) = "
        f"{b_cross:.4f} bp/frag")

    # --------------------------------------------------------------
    # Reference rho_ig and rho_in
    # --------------------------------------------------------------
    log(section("REFERENCE rho_ig and rho_in (containment density)"))

    types = region_df["type"].to_numpy()
    spans = (region_df["end"].to_numpy().astype(np.int64)
             - region_df["start"].to_numpy().astype(np.int64))
    leff = l_eff_contained(spans, gdna_pmf)
    prc = payload.per_region_counts

    is_ig = (types == 0)
    is_in = (types == 1)
    is_ex = (types == 2)

    n_ig_frag = int(prc[is_ig, MASK_INTERGENIC].sum())
    L_ig = float(leff[is_ig].sum())
    rho_ig = n_ig_frag / max(L_ig, 1e-9)
    log(f"  rho_ig    = {rho_ig:.6e}   (n={n_ig_frag:,}, L_eff={L_ig:,.0f})")

    n_in_frag = int(prc[is_in, MASK_INTRON].sum())
    L_in = float(leff[is_in].sum())
    rho_in = n_in_frag / max(L_in, 1e-9)
    log(f"  rho_in    = {rho_in:.6e}   (n={n_in_frag:,}, L_eff={L_in:,.0f})")

    # --------------------------------------------------------------
    # (A) rho_ex split by side eligibility
    # --------------------------------------------------------------
    log(section("(A) rho_ex SPLIT BY SIDE-ELIGIBILITY OF EXON ROW"))

    bfl = region_df["boundary_flux_left"].to_numpy().astype(np.int64)
    bfr = region_df["boundary_flux_right"].to_numpy().astype(np.int64)
    sides = bfl + bfr  # 0/1/2

    u_left = payload.u_left
    u_right = payload.u_right

    log(f"  {'category':<40} {'n_regions':>9}  {'count':>10}  "
        f"{'L_eff':>15}  {'rho':>13}  {'rho/rho_ig':>12}")

    def _summary(label: str, mask: np.ndarray) -> None:
        sel = is_ex & mask
        cnt = int((u_left[sel] * bfl[sel] + u_right[sel] * bfr[sel]).sum())
        leff_sum = float((sides[sel].astype(np.float64) * b_cross).sum())
        nreg = int(sel.sum())
        rho = cnt / max(leff_sum, 1e-9)
        log(f"  {label:<40} {nreg:>9,}  {cnt:>10,}  "
            f"{leff_sum:>15,.1f}  {rho:>13.6e}  "
            f"{rho/max(rho_ig,1e-12):>12.4f}")

    _summary("ALL eligible (sides>0)",                sides > 0)
    _summary("INTERNAL (sides==2; both eligible)",    sides == 2)
    _summary("TERMINAL (sides==1; one side only)",    sides == 1)
    _summary("INELIGIBLE (sides==0)",                 sides == 0)

    log("\n  >>> If TERMINAL rho_ex << INTERNAL rho_ex, the geometry")
    log("      hypothesis is supported. If they are comparable,")
    log("      look at (B) and (C).")

    # Also break down by which single side is eligible (for terminals).
    log("\n  Terminal-side breakdown (sides==1):")
    sel_tL = is_ex & (sides == 1) & (bfl == 1)
    sel_tR = is_ex & (sides == 1) & (bfr == 1)
    for label, sel in (("only LEFT  eligible", sel_tL),
                       ("only RIGHT eligible", sel_tR)):
        nreg = int(sel.sum())
        cnt = int((u_left[sel] * bfl[sel] + u_right[sel] * bfr[sel]).sum())
        leff_sum = float((sides[sel].astype(np.float64) * b_cross).sum())
        rho = cnt / max(leff_sum, 1e-9)
        log(f"    {label}: n_regions={nreg:>6,}  count={cnt:>8,}  "
            f"rho={rho:.6e}  rho/rho_ig={rho/max(rho_ig,1e-12):.4f}")

    # --------------------------------------------------------------
    # (B) per-EXON observed/expected ratio by width bin
    # --------------------------------------------------------------
    log(section("(B) PER-EXON OBSERVED/EXPECTED RATIO BY WIDTH BIN"))

    eligible = is_ex & (sides > 0)
    ex_widths = spans[eligible]
    ex_obs = (u_left[eligible] * bfl[eligible]
              + u_right[eligible] * bfr[eligible]).astype(np.float64)
    ex_exp = sides[eligible].astype(np.float64) * b_cross * rho_ig

    bins = [0, 100, 200, 400, 800, 1600, 3200, 10**9]
    bin_labels = [
        f"{bins[i]}-{bins[i+1]-1 if bins[i+1]<10**9 else 'inf'}"
        for i in range(len(bins) - 1)
    ]
    bidx = np.digitize(ex_widths, bins) - 1
    log(f"  {'width_bp_bin':<14} {'n_regions':>10} {'sum_obs':>12} "
        f"{'sum_exp':>14} {'obs/exp':>9} {'mean_obs/region':>16}")
    for b in range(len(bin_labels)):
        sel = bidx == b
        n = int(sel.sum())
        if n == 0:
            continue
        s_obs = float(ex_obs[sel].sum())
        s_exp = float(ex_exp[sel].sum())
        ratio = s_obs / max(s_exp, 1e-9)
        per_reg = s_obs / max(n, 1)
        log(f"  {bin_labels[b]:<14} {n:>10,} {s_obs:>12,.0f} "
            f"{s_exp:>14,.1f} {ratio:>9.4f} {per_reg:>16.3f}")
    log("  >>> Flat ratio across width bins => uniform deficit (FL/aligner).")
    log("      Ratio that grows with width  => small-exon double-cross interaction.")

    # --------------------------------------------------------------
    # Headline
    # --------------------------------------------------------------
    log(section("HEADLINE"))
    n_ex_frag = int(ex_obs.sum())
    L_ex = float((sides[eligible].astype(np.float64) * b_cross).sum())
    rho_ex = n_ex_frag / max(L_ex, 1e-9)
    log(f"  rho_ex           = {rho_ex:.6e}")
    log(f"  rho_ig           = {rho_ig:.6e}")
    log(f"  rho_in           = {rho_in:.6e}")
    log(f"  rho_ex / rho_ig  = {rho_ex / max(rho_ig, 1e-12):.4f}")
    log(f"  rho_ex / rho_in  = {rho_ex / max(rho_in, 1e-12):.4f}")

    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True)
        Path(args.out).write_text("\n".join(out_lines) + "\n")
        log(f"\n[probe] wrote report -> {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
