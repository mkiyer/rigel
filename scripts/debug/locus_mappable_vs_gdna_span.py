"""Correlate per-locus (mappable_bp / gdna_span) with gdna_rate.

Hypothesis: the EM's gDNA bias-correction denominator uses
``gdna_span + gdna_flank`` (intron-inclusive genomic span), while the
calibration's λ_G is estimated per mappable bp.  In loci with a large
unmappable fraction, gdna_span >> mappable_bp causes the gDNA component
to be over-penalized and its rate under-estimated.

This script:
  1. Loads loci.feather + quant.feather for a chosen run,
  2. Reconstructs per-locus merged genomic intervals from member
     transcripts (same algorithm as locus.build_loci),
  3. Overlaps those intervals with regions.feather to sum
     mappable_effective_length per locus,
  4. Correlates mappable_bp/gdna_span ratio with gdna_rate and
     stratifies by locus size.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd


def reconstruct_locus_intervals(tx: pd.DataFrame, locus_ids: pd.Series) -> dict:
    """Return {locus_id: [(ref, start, end), ...]} merged genomic intervals."""
    merged = {}
    df = tx.assign(locus_id=locus_ids.values)
    df = df[df["locus_id"] >= 0]
    for lid, sub in df.groupby("locus_id"):
        sub = sub.sort_values(["ref", "start"])
        intervals = []
        prev_ref, prev_s, prev_e = None, None, None
        for ref, s, e in zip(sub["ref"].values, sub["start"].values, sub["end"].values):
            if prev_ref is None:
                prev_ref, prev_s, prev_e = ref, int(s), int(e)
                continue
            if ref != prev_ref or s > prev_e:
                intervals.append((prev_ref, prev_s, prev_e))
                prev_ref, prev_s, prev_e = ref, int(s), int(e)
            else:
                prev_e = max(prev_e, int(e))
        intervals.append((prev_ref, prev_s, prev_e))
        merged[int(lid)] = intervals
    return merged


def sum_mappable_bp_per_locus(locus_intervals: dict, regions: pd.DataFrame) -> pd.DataFrame:
    """For each locus, sum overlapping regions' mappable_effective_length and length."""
    # Pre-index regions by ref for fast interval overlap.
    rows = []
    by_ref = {ref: g.sort_values("start").reset_index(drop=True) for ref, g in regions.groupby("ref")}
    for lid, ivals in locus_intervals.items():
        mbp_sum = 0.0
        len_sum = 0  # sum of region.length overlapped (with clipping)
        span_sum = 0  # sum of locus interval spans
        for ref, s, e in ivals:
            span_sum += e - s
            rg = by_ref.get(ref)
            if rg is None:
                continue
            # Binary search for regions with start < e and end > s
            starts = rg["start"].values
            ends = rg["end"].values
            # Regions that might overlap: end > s and start < e
            i0 = np.searchsorted(ends, s, side="right")
            i1 = np.searchsorted(starts, e, side="left")
            if i0 >= i1:
                continue
            # For each, compute overlap fraction and apportion mappable_bp
            rstart = starts[i0:i1]
            rend = ends[i0:i1]
            rmbp = rg["mappable_effective_length"].values[i0:i1]
            rlen = rg["length"].values[i0:i1]
            overlap = np.minimum(rend, e) - np.maximum(rstart, s)
            overlap = np.clip(overlap, 0, None)
            frac = np.where(rlen > 0, overlap / rlen, 0.0)
            mbp_sum += float((rmbp * frac).sum())
            len_sum += int(overlap.sum())
        rows.append({"locus_id": lid, "mappable_bp": mbp_sum, "overlap_len": len_sum, "span": span_sum})
    return pd.DataFrame(rows)


def main():
    run_dir = Path(
        sys.argv[1] if len(sys.argv) > 1
        else "/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna/runs/human/mctp_vcap_rna20m_dna80m"
    )
    index_dir = Path(
        sys.argv[2] if len(sys.argv) > 2
        else "/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna/refs/human/rigel_index"
    )

    print(f"[load] run_dir={run_dir}")
    loci = pd.read_feather(run_dir / "rigel" / "loci.feather")
    quant = pd.read_feather(run_dir / "rigel" / "quant.feather")
    tx = pd.read_feather(index_dir / "transcripts.feather")
    regions = pd.read_feather(index_dir / "regions.feather")
    print(f"[load] {len(loci)} loci, {len(quant)} transcripts in quant, {len(tx)} in index, {len(regions)} regions")

    # Join quant → tx to attach locus_id.
    tx = tx.merge(
        quant[["transcript_id", "locus_id"]],
        left_on="t_id", right_on="transcript_id", how="left"
    )
    tx["locus_id"] = tx["locus_id"].fillna(-1).astype(int)
    print(f"[join] tx with locus_id: {(tx['locus_id'] >= 0).sum()}")

    print("[reconstruct] merged locus intervals ...")
    locus_intervals = reconstruct_locus_intervals(tx, tx["locus_id"])
    print(f"[reconstruct] {len(locus_intervals)} loci")

    print("[overlap] sum mappable_bp per locus ...")
    mbp = sum_mappable_bp_per_locus(locus_intervals, regions)

    # Join with loci.feather
    df = loci.merge(mbp, on="locus_id", how="left")
    df["mbp_ratio"] = df["mappable_bp"] / df["locus_span_bp"].clip(lower=1)
    df["size_bin"] = pd.cut(
        df["locus_span_bp"],
        bins=[0, 1e4, 1e5, 1e6, 1e7, 1e10],
        labels=["<10kb", "10-100kb", "100kb-1Mb", "1-10Mb", ">10Mb"],
    )

    # Correlation (only loci with non-trivial EM signal)
    has_data = (df["n_em_fragments"] >= 10) & df["mappable_bp"].notna()
    sub = df[has_data].copy()
    print(f"\n[overall] n loci with ≥10 EM fragments: {len(sub)}")
    print(f"  mbp_ratio: mean={sub['mbp_ratio'].mean():.3f}  median={sub['mbp_ratio'].median():.3f}")
    print(f"  gdna_rate: mean={sub['gdna_rate'].mean():.3f}  median={sub['gdna_rate'].median():.3f}")
    r_p = sub[["mbp_ratio", "gdna_rate"]].corr().iloc[0, 1]
    r_s = sub[["mbp_ratio", "gdna_rate"]].corr(method="spearman").iloc[0, 1]
    print(f"  Pearson(mbp_ratio, gdna_rate)  = {r_p:+.4f}")
    print(f"  Spearman(mbp_ratio, gdna_rate) = {r_s:+.4f}")

    print("\n[by size_bin] weighted by n_em_fragments:")
    grp = sub.groupby("size_bin", observed=True).agg(
        n_loci=("locus_id", "size"),
        frags=("n_em_fragments", "sum"),
        gdna_rate_mean=("gdna_rate", "mean"),
        mbp_ratio_mean=("mbp_ratio", "mean"),
        mbp_ratio_median=("mbp_ratio", "median"),
        gdna_rate_wmean=("gdna_rate", lambda s: np.average(s, weights=sub.loc[s.index, "n_em_fragments"])),
    )
    print(grp.to_string())

    # Focus on mega-loci (span > 10 Mb) — print top 5 by fragments
    print("\n[top mega-loci by n_em_fragments]")
    mega = sub[sub["locus_span_bp"] > 1e7].sort_values("n_em_fragments", ascending=False).head(5)
    cols = ["locus_id", "locus_span_bp", "mappable_bp", "mbp_ratio",
            "n_em_fragments", "mrna", "nrna", "gdna", "gdna_rate", "gdna_prior"]
    print(mega[cols].to_string(index=False))

    # Save full table for user
    out = run_dir / "rigel" / "locus_mbp_diagnostic.tsv"
    df.to_csv(out, sep="\t", index=False)
    print(f"\n[save] {out}")


if __name__ == "__main__":
    main()
