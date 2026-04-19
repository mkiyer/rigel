"""Phase 0 probe — alignable Zarr landscape sanity check.

Goals (per docs/mappability/unified_proposal_v4.md, Phase 0):

1. Confirm the alignable Zarr store opens and exposes the documented API.
2. For the existing rigel index regions, compute the per-region effective
   length E_i = sum of fractional_mappability across the region at
   read_length=100 — the proposed default for v4.
3. Compare against the existing binary uniquely-mappable fraction
   (the v3-era signal we are replacing).
4. Spot-check a handful of regions (centromeric, exonic, repeat-rich)
   for biological plausibility.
5. Profile throughput: is sequential per-chromosome scanning fast
   enough for index-build time on ~700k regions?

This script is read-only and produces stdout diagnostics + an optional
TSV dump of per-region E_i for downstream review. Run it from the
``alignable`` conda environment::

    conda activate alignable
    python scripts/debug/alignable_landscape_probe.py \
        --zarr /scratch/.../alignable/grch38_star \
        --rigel-index /scratch/.../hulkrna/refs/human/rigel_index \
        --read-length 100 \
        --tsv-out /tmp/region_eff_length.tsv
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

import alignable
from alignable.api import AlignableStore


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--zarr", required=True, help="Path to alignable Zarr store or output dir")
    p.add_argument(
        "--rigel-index",
        required=True,
        help="Path to a rigel index directory containing regions.feather",
    )
    p.add_argument(
        "--read-length", type=int, default=100, help="Read-length bin (default 100)"
    )
    p.add_argument(
        "--max-regions",
        type=int,
        default=0,
        help="If >0, sample this many regions (for quick runs). 0 = all regions",
    )
    p.add_argument(
        "--cache-size",
        type=int,
        default=16,
        help="alignable.open(cache_size=...) — number of decompressed chunks to cache",
    )
    p.add_argument(
        "--tsv-out",
        type=str,
        default="",
        help="Optional path to write per-region (ref, start, end, length, eff_length, "
        "binary_mappable_length, fi) TSV",
    )
    p.add_argument(
        "--spot-check",
        type=int,
        default=10,
        help="Number of regions to print for visual inspection at each percentile",
    )
    return p.parse_args()


def load_regions(index_dir: Path) -> pd.DataFrame:
    rf = index_dir / "regions.feather"
    if not rf.exists():
        raise SystemExit(f"regions.feather not found in {index_dir}")
    df = pd.read_feather(rf)
    keep = [c for c in ("ref", "start", "end", "length") if c in df.columns]
    df = df[keep].copy()
    if "length" not in df.columns:
        df["length"] = df["end"] - df["start"]
    return df


def compute_effective_lengths(
    store: AlignableStore,
    regions: pd.DataFrame,
    read_length: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Return (eff_length, binary_mappable_length) per region.

    eff_length = sum(fractional_mappability) over region (1/NH semantics).
    binary_mappable_length = sum(binary_mappability) (uniquely-mappable bases).

    Regions are processed in (ref, start) order so that adjacent regions
    hit the same Zarr chunks and benefit from the LRU cache.
    """
    chrom_set = set(store.chromosomes.keys())
    n = len(regions)
    eff = np.zeros(n, dtype=np.float64)
    binm = np.zeros(n, dtype=np.float64)
    skipped = 0

    # Sort by (ref, start) for cache locality, but remember original order.
    order = np.lexsort((regions["start"].values, regions["ref"].values))

    last_chrom = None
    for ii, idx in enumerate(order):
        ref = regions["ref"].iat[idx]
        start = int(regions["start"].iat[idx])
        end = int(regions["end"].iat[idx])
        if ref not in chrom_set:
            skipped += 1
            continue
        if start >= end:
            continue
        # Clamp to chromosome length defensively.
        clen = store.chromosomes[ref]
        s = max(0, start)
        e = min(end, clen)
        if s >= e:
            continue
        try:
            fm = store.fractional_mappability(ref, s, e, read_length)
            bm = store.binary_mappability(ref, s, e, read_length)
        except Exception as exc:
            if last_chrom != ref:
                print(f"  query failure on {ref}:{s}-{e}: {exc}")
                last_chrom = ref
            continue
        eff[idx] = float(fm.sum())
        binm[idx] = float(bm.sum())

    if skipped:
        print(f"  skipped {skipped:,} regions on chromosomes not present in the Zarr")
    return eff, binm


def main() -> int:
    args = parse_args()

    print(f"Opening alignable store: {args.zarr}")
    # NOTE: alignable.open() does not currently forward kwargs to
    # AlignableStore (README documents cache_size= but open() ignores it).
    # Construct directly to honour --cache-size. Reported upstream.
    store = AlignableStore(args.zarr, cache_size=args.cache_size)
    print(
        f"  aligner={store.aligner} v{store.aligner_version}  "
        f"status={store.status}  read_length_bins={store.read_length_bins}  "
        f"frag_len_mean={store.frag_len_mean}  error_rate={store.error_rate}"
    )
    if args.read_length not in store.read_length_bins:
        raise SystemExit(
            f"read_length={args.read_length} not in store bins {store.read_length_bins}"
        )

    regions = load_regions(Path(args.rigel_index))
    print(f"Loaded {len(regions):,} regions from {args.rigel_index}/regions.feather")
    print(f"  ref values: {regions['ref'].nunique():,} unique")

    if args.max_regions > 0 and len(regions) > args.max_regions:
        regions = regions.sample(n=args.max_regions, random_state=0).reset_index(drop=True)
        print(f"  sampled {len(regions):,} regions for this probe run")

    print(f"\nComputing per-region E_i at read_length={args.read_length} ...")
    t0 = time.perf_counter()
    eff, binm = compute_effective_lengths(store, regions, args.read_length)
    dt = time.perf_counter() - t0
    print(
        f"  done in {dt:.1f}s ({len(regions) / dt:,.0f} regions/sec)  "
        f"cache={store.cache_info()}"
    )

    L = regions["length"].astype(np.float64).values
    fi = np.divide(eff, L, out=np.zeros_like(eff), where=L > 0)
    fi_bin = np.divide(binm, L, out=np.zeros_like(binm), where=L > 0)

    print("\n=== Per-region length distribution ===")
    print(pd.Series(L).describe(percentiles=[0.05, 0.25, 0.5, 0.75, 0.95]).round(1))

    print("\n=== Effective length E_i = sum(fractional_mappability) ===")
    print(pd.Series(eff).describe(percentiles=[0.05, 0.25, 0.5, 0.75, 0.95]).round(1))

    print("\n=== f_i (mean fractional mappability) ===")
    print(pd.Series(fi).describe(percentiles=[0.05, 0.25, 0.5, 0.75, 0.95]).round(3))

    print("\n=== f_i_binary (mean binary mappability — old uniquely-mappable signal) ===")
    print(pd.Series(fi_bin).describe(percentiles=[0.05, 0.25, 0.5, 0.75, 0.95]).round(3))

    print("\n=== Coverage summary ===")
    total_L = float(L.sum())
    total_eff = float(eff.sum())
    total_binm = float(binm.sum())
    print(f"  total region span:                {total_L:,.0f} bp")
    print(
        f"  total effective length (E):       {total_eff:,.0f} bp "
        f"({100 * total_eff / max(total_L, 1):.1f}% of span)"
    )
    print(
        f"  total uniquely-mappable bases:    {total_binm:,.0f} bp "
        f"({100 * total_binm / max(total_L, 1):.1f}% of span)"
    )
    delta = total_eff - total_binm
    print(
        f"  multi-mapper credit (E - binary): {delta:,.0f} bp "
        f"({100 * delta / max(total_L, 1):.2f}% of span) "
        "<- bases recovered by going continuous"
    )

    # Sanity: f_i should be >= f_i_binary at every region
    bad = (fi + 1e-6 < fi_bin).sum()
    print(f"  regions where f_i < f_i_binary (should be 0): {bad}")

    # Effective-length floor diagnostic from v4 (E >= L_bar fragment length).
    # alignable's frag_len_mean is the synthetic fragment length; rigel's
    # actual L_bar comes from real data. We use the alignable value as a
    # ballpark for this probe.
    floor = float(store.frag_len_mean)
    admitted = int((eff >= floor).sum())
    print(
        f"\n=== Effective-length floor check (E >= {floor:.0f} bp ~= L_bar) ==="
    )
    print(
        f"  admitted regions: {admitted:,} / {len(regions):,} "
        f"({100 * admitted / len(regions):.1f}%)"
    )
    print(
        f"  admitted effective span: {float(eff[eff >= floor].sum()):,.0f} bp "
        f"({100 * float(eff[eff >= floor].sum()) / max(total_eff, 1):.1f}% of total E)"
    )

    # Spot check
    n_spot = args.spot_check
    if n_spot > 0:
        print(f"\n=== Spot-check: {n_spot} highest-E and {n_spot} lowest-E regions ===")
        df = regions.copy()
        df["E"] = eff
        df["f_i"] = fi
        df["f_i_binary"] = fi_bin
        cols = ["ref", "start", "end", "length", "E", "f_i", "f_i_binary"]
        print("  Top E:")
        print(df.nlargest(n_spot, "E")[cols].to_string(index=False))
        print("  Bottom-E (with length > 1 kb to avoid trivial cases):")
        big = df[df["length"] > 1000]
        if len(big) > 0:
            print(big.nsmallest(n_spot, "E")[cols].to_string(index=False))

    if args.tsv_out:
        outp = Path(args.tsv_out)
        outp.parent.mkdir(parents=True, exist_ok=True)
        df = regions.copy()
        df["E"] = eff
        df["binary_mappable_length"] = binm
        df["f_i"] = fi
        df["f_i_binary"] = fi_bin
        df.to_csv(outp, sep="\t", index=False)
        print(f"\nWrote per-region table to {outp}")

    summary = {
        "zarr": args.zarr,
        "read_length": args.read_length,
        "n_regions": int(len(regions)),
        "elapsed_sec": dt,
        "regions_per_sec": len(regions) / max(dt, 1e-9),
        "total_span_bp": total_L,
        "total_eff_bp": total_eff,
        "total_binary_bp": total_binm,
        "fi_median": float(np.median(fi)),
        "fi_binary_median": float(np.median(fi_bin)),
        "admitted_at_floor": admitted,
        "admit_floor_bp": floor,
        "cache": str(store.cache_info()),
    }
    print("\n=== JSON summary ===")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
