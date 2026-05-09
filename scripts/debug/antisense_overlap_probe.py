"""Diagnose antisense-overlap impact on the EXON-INTRON density estimator.

The genic clustering in `_emit_layout` (rigel.index) is STRAND-AGNOSTIC:
any (+)/(-) strand overlap merges into one `_GenicSpan`. Within the
span, an EXON region is flagged wherever ANY-strand has exon coverage.
This means a (+)-strand intron position that overlaps a (-)-strand
exon becomes an EXON region in the partition.

Consequences hypothesized:

  H1. Antisense-only-overlapping (AMBIG) EXON regions inflate the
      EXON denominator (sides * B_cross) without proportionate
      boundary-crossing numerator, because boundary flux at the
      "internal" boundaries of an AMBIG EXON region is driven by the
      LESS-EXPRESSED strand's gene structure, not the dominant one.

  H2. AMBIG EXON regions partition INTRON runs of one strand into
      shorter pieces, possibly removing genuine boundary-crossing
      events that would otherwise be counted.

This script:
  - loads the index and reports strand composition of EXON regions;
  - splits ρ_ex by region strand:  POS-only, NEG-only, AMBIG;
  - reports # of distinct genic clusters affected and how many
    distinct genes contribute to each AMBIG cluster.
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

import numpy as np

from rigel.config import BamScanConfig
from rigel.index import TranscriptIndex
from rigel.calibration.regions import RegionStrand, RegionType
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
    p.add_argument("--K", type=int, default=3)
    return p.parse_args()


def fl_summary(hist: np.ndarray) -> tuple[np.ndarray, float, int]:
    n = int(hist.sum())
    if n == 0:
        return np.zeros_like(hist, dtype=np.float64), 0.0, 0
    pmf = hist.astype(np.float64) / float(n)
    mean = float((np.arange(len(hist), dtype=np.float64) * pmf).sum())
    return pmf, mean, n


def b_cross_from_pmf(pmf: np.ndarray, K: int) -> float:
    q = max(K, 1)
    ells = np.arange(len(pmf), dtype=np.float64)
    return float((pmf * np.maximum(ells - 2.0 * q + 1.0, 0.0)).sum())


def l_eff_contained(spans: np.ndarray, pmf: np.ndarray) -> np.ndarray:
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

    print(f"[probe] BAM   : {args.bam}")
    print(f"[probe] INDEX : {args.index}")
    print(f"[probe] K     : {args.K}")

    print("\n[probe] loading index ...")
    index = TranscriptIndex.load(args.index)
    region_df = index.region_df.reset_index(drop=True)
    t_df = index.t_df

    print(section("STRAND COMPOSITION OF THE PARTITION"))
    types = region_df["type"].to_numpy()
    strands = region_df["strand"].to_numpy()
    spans = (region_df["end"].to_numpy(np.int64)
             - region_df["start"].to_numpy(np.int64))

    type_names = {0: "INTERGENIC", 1: "INTRON", 2: "EXON"}
    strand_names = {
        int(RegionStrand.NONE):  "NONE",
        int(RegionStrand.POS):   "POS",
        int(RegionStrand.NEG):   "NEG",
        int(RegionStrand.AMBIG): "AMBIG",
    }
    print(f"  {'type':<11} {'strand':<6} {'n_regions':>10} {'sum_bp':>14} {'mean_bp':>10}")
    for tval in (0, 1, 2):
        for sval in (0, 1, 2, 3):
            sel = (types == tval) & (strands == sval)
            n = int(sel.sum())
            if n == 0:
                continue
            ssum = int(spans[sel].sum())
            print(f"  {type_names[tval]:<11} {strand_names[sval]:<6} "
                  f"{n:>10,} {ssum:>14,} {ssum/n:>10,.0f}")

    # AMBIG EXON regions are the suspected pain point. Quantify by
    # comparing partition row counts (already printed above).
    print(section("PARTITION-LEVEL ANTISENSE FOOTPRINT"))
    n_ambig_ex = int(((types == 2) & (strands == int(RegionStrand.AMBIG))).sum())
    n_ex = int((types == 2).sum())
    n_ambig_in = int(((types == 1) & (strands == int(RegionStrand.AMBIG))).sum())
    n_in = int((types == 1).sum())
    bp_ambig_ex = int(spans[(types == 2) & (strands == int(RegionStrand.AMBIG))].sum())
    bp_ex = int(spans[types == 2].sum())
    print(f"  AMBIG EXON   rows: {n_ambig_ex:>4} / {n_ex:>4} "
          f"({100*n_ambig_ex/max(n_ex,1):.1f}%)  bp share={100*bp_ambig_ex/max(bp_ex,1):.1f}%")
    print(f"  AMBIG INTRON rows: {n_ambig_in:>4} / {n_in:>4} "
          f"({100*n_ambig_in/max(n_in,1):.1f}%)")

    # Run the BAM scan
    print("\n[probe] running scan_and_buffer ...")
    payload = scan_and_buffer(args.bam, index,
                              BamScanConfig(splicing_anchor_tolerance=args.K))[-1]
    if payload is None:
        print("ERROR: payload is None")
        return 1

    pmf_ig, mean_ig, n_ig = fl_summary(payload.fl_hist[MASK_INTERGENIC])
    b_cross = b_cross_from_pmf(pmf_ig, args.K)
    leff = l_eff_contained(spans, pmf_ig)

    is_ig = (types == 0); is_in = (types == 1); is_ex = (types == 2)
    prc = payload.per_region_counts

    n_ig_frag = int(prc[is_ig, MASK_INTERGENIC].sum())
    L_ig = float(leff[is_ig].sum())
    rho_ig = n_ig_frag / max(L_ig, 1e-9)

    n_in_frag = int(prc[is_in, MASK_INTRON].sum())
    L_in = float(leff[is_in].sum())
    rho_in = n_in_frag / max(L_in, 1e-9)

    print(f"\n  rho_ig = {rho_ig:.6e},  rho_in = {rho_in:.6e}")
    print(f"  B_cross(K={args.K}) = {b_cross:.4f} bp/frag")

    # --------------------------------------------------------------
    # rho_ex split by EXON region strand
    # --------------------------------------------------------------
    print(section("rho_ex SPLIT BY EXON REGION STRAND"))

    bfl = region_df["boundary_flux_left"].to_numpy(np.int64)
    bfr = region_df["boundary_flux_right"].to_numpy(np.int64)
    sides = bfl + bfr
    u_left = payload.u_left
    u_right = payload.u_right

    print(f"  {'strand':<8} {'n_eligible':>11} {'count':>10} "
          f"{'L_eff':>14} {'rho_ex':>13} {'rho_ex/rho_ig':>14}")
    rows = []
    for sval in (int(RegionStrand.POS), int(RegionStrand.NEG),
                 int(RegionStrand.AMBIG)):
        sel = is_ex & (strands == sval) & (sides > 0)
        n = int(sel.sum())
        if n == 0:
            continue
        cnt = int((u_left[sel] * bfl[sel] + u_right[sel] * bfr[sel]).sum())
        L_eff = float((sides[sel].astype(np.float64) * b_cross).sum())
        rho = cnt / max(L_eff, 1e-9)
        rows.append((strand_names[sval], n, cnt, L_eff, rho))
        print(f"  {strand_names[sval]:<8} {n:>11,} {cnt:>10,} "
              f"{L_eff:>14,.1f} {rho:>13.6e} {rho/max(rho_ig,1e-12):>14.4f}")

    # AMBIG denominator share
    sel_amb = is_ex & (strands == int(RegionStrand.AMBIG)) & (sides > 0)
    sel_all_eligible = is_ex & (sides > 0)
    L_amb = float((sides[sel_amb].astype(np.float64) * b_cross).sum())
    L_all = float((sides[sel_all_eligible].astype(np.float64) * b_cross).sum())
    if L_all > 0:
        print(f"\n  AMBIG share of EXON denominator: "
              f"L_amb/L_all = {L_amb/L_all:.4f} ({100*L_amb/L_all:.1f}%)")

    # --------------------------------------------------------------
    # Per-region containment density on EXON rows split by strand:
    # this checks whether AMBIG regions are 'genuinely transcribed'
    # vs 'protected by capture'. If AMBIG EXON regions have
    # n[INTRON-mask] > 0 fragments contained, that's a tell that they
    # behave more like INTRON than EXON.
    # --------------------------------------------------------------
    print(section("EXON-row CONTAINMENT counts (mask-conditioned)"))
    print(f"  {'strand':<8} {'n_regions':>10} {'n_EXON_frag':>14} "
          f"{'n_INTRON_frag':>14} {'n_BOUND_frag':>14}")
    for sval in (int(RegionStrand.POS), int(RegionStrand.NEG),
                 int(RegionStrand.AMBIG)):
        sel = is_ex & (strands == sval)
        n = int(sel.sum())
        if n == 0:
            continue
        n_ex = int(prc[sel, MASK_EXON].sum())
        n_in_ = int(prc[sel, MASK_INTRON].sum())
        n_bd = int(prc[sel, MASK_EXON | MASK_INTRON].sum())
        print(f"  {strand_names[sval]:<8} {n:>10,} {n_ex:>14,} "
              f"{n_in_:>14,} {n_bd:>14,}")

    # --------------------------------------------------------------
    # Without AMBIG: re-compute headline.
    # --------------------------------------------------------------
    print(section("HEADLINE WITH AND WITHOUT AMBIG EXON ROWS"))
    sel_ALL = is_ex & (sides > 0)
    sel_NA  = is_ex & (sides > 0) & (strands != int(RegionStrand.AMBIG))
    for label, sel in (("ALL eligible", sel_ALL),
                       ("non-AMBIG only", sel_NA),
                       ("AMBIG only", sel_amb)):
        cnt = int((u_left[sel] * bfl[sel] + u_right[sel] * bfr[sel]).sum())
        L_eff = float((sides[sel].astype(np.float64) * b_cross).sum())
        rho = cnt / max(L_eff, 1e-9)
        print(f"  {label:<18}  rho_ex={rho:.6e}  rho_ex/rho_ig="
              f"{rho/max(rho_ig,1e-12):.4f}  (n={int(sel.sum())} regions, "
              f"count={cnt}, L_eff={L_eff:,.0f})")

    return 0


if __name__ == "__main__":
    sys.exit(main())
