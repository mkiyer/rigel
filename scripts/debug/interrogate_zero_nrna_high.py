#!/usr/bin/env python3
"""Interrogate the gdna_zero_ss_0.50_nrna_high synthetic-24 condition.

The goal is to understand why true nRNA fragments receive zero nRNA assignments.
This script uses the annotated BAM as the fragment-level source of truth and joins
read-name origins back to Rigel's transcript/nRNA/locus outputs.
"""

from __future__ import annotations

import argparse
import math
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pysam

from rigel.sim.truth import parse_origin


AF_MRNA_BIT = 0x02
AF_GDNA_BIT = 0x04
AF_NRNA_BIT = 0x08

DEFAULT_BASE = Path("/Users/mkiyer/Downloads/rigel_runs/sim_synthetic_24")
DEFAULT_CONDITION = "gdna_zero_ss_0.50_nrna_high"
DEFAULT_OUT = Path("results/synthetic_24_deep_analysis/gdna_zero_ss_0.50_nrna_high")


def assigned_pool(zf_value: int) -> str:
    if zf_value & AF_GDNA_BIT:
        return "gdna"
    if zf_value & AF_NRNA_BIT:
        return "nrna"
    if zf_value & AF_MRNA_BIT:
        return "mrna"
    return "unresolved"


def get_tag(read: pysam.AlignedSegment, tag: str, default: Any) -> Any:
    try:
        return read.get_tag(tag)
    except KeyError:
        return default


def parse_gtf_attributes(text: str) -> dict[str, str]:
    return {match.group(1): match.group(2) for match in re.finditer(r'(\S+)\s+"([^"]*)"', text)}


def parse_gtf_exons(gtf_path: Path) -> dict[str, tuple[tuple[int, int], ...]]:
    exons: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with open(gtf_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            attrs = parse_gtf_attributes(fields[8])
            transcript_id = attrs.get("transcript_id")
            if transcript_id:
                exons[transcript_id].append((int(fields[3]) - 1, int(fields[4])))
    return {transcript_id: tuple(sorted(items)) for transcript_id, items in exons.items()}


def merge_blocks(blocks: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not blocks:
        return []
    blocks = sorted(blocks)
    merged = [blocks[0]]
    for start, end in blocks[1:]:
        prev_start, prev_end = merged[-1]
        if start <= prev_end:
            merged[-1] = (prev_start, max(prev_end, end))
        else:
            merged.append((start, end))
    return merged


def overlap_len(left: tuple[int, int], right: tuple[int, int]) -> int:
    start = max(left[0], right[0])
    end = min(left[1], right[1])
    return max(0, end - start)


def classify_blocks(
    blocks: list[tuple[int, int]],
    exons: tuple[tuple[int, int], ...],
) -> tuple[str, int, int]:
    merged = merge_blocks(blocks)
    aligned = sum(end - start for start, end in merged)
    if aligned == 0:
        return "unaligned", 0, 0
    exon_overlap = 0
    for block in merged:
        exon_overlap += sum(overlap_len(block, exon) for exon in exons)
    if exon_overlap == 0:
        return "intronic", aligned, exon_overlap
    if exon_overlap == aligned:
        return "exon_contained", aligned, exon_overlap
    return "exon_intron", aligned, exon_overlap


def first_primary(records: list[pysam.AlignedSegment]) -> pysam.AlignedSegment | None:
    for read in records:
        if not read.is_secondary and not read.is_supplementary:
            return read
    return None


def primary_blocks(records: list[pysam.AlignedSegment]) -> list[tuple[int, int]]:
    blocks: list[tuple[int, int]] = []
    for read in records:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        blocks.extend(read.get_blocks())
    return merge_blocks(blocks)


def summarize_numeric(values: list[float]) -> dict[str, float]:
    if not values:
        return {"n": 0, "min": math.nan, "q01": math.nan, "median": math.nan, "q99": math.nan, "max": math.nan, "mean": math.nan}
    arr = np.asarray(values, dtype=float)
    return {
        "n": int(arr.size),
        "min": float(np.min(arr)),
        "q01": float(np.quantile(arr, 0.01)),
        "median": float(np.median(arr)),
        "q99": float(np.quantile(arr, 0.99)),
        "max": float(np.max(arr)),
        "mean": float(np.mean(arr)),
    }


def scan_annotated_bam(
    bam_path: Path,
    exons_by_tx: dict[str, tuple[tuple[int, int], ...]],
    origin_to_nrna_id: dict[str, str],
    nrna_locus_by_id: dict[str, int],
    *,
    max_examples_per_pool: int = 8,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    examples = []
    example_counts: Counter[str] = Counter()

    frag_id = -1

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        current_name: str | None = None
        group: list[pysam.AlignedSegment] = []

        def flush_group() -> None:
            nonlocal frag_id
            if not group:
                return
            frag_id += 1
            primary = first_primary(group)
            if primary is None:
                return
            origin = parse_origin(primary.query_name)
            zf_value = int(get_tag(primary, "ZF", 0))
            pool = assigned_pool(zf_value)
            origin_tx = origin.transcript_id or ""
            blocks = primary_blocks(group)
            region_class, aligned_bases, exon_bases = classify_blocks(
                blocks,
                exons_by_tx.get(origin_tx, ()),
            )
            nrna_id = origin_to_nrna_id.get(origin_tx, "")
            inferred_locus = nrna_locus_by_id.get(nrna_id, -1)
            row = {
                "frag_id": frag_id,
                "qname": primary.query_name,
                "origin": origin.kind,
                "origin_tx": origin_tx,
                "origin_nrna_id": nrna_id,
                "inferred_origin_locus": inferred_locus,
                "assigned_pool": pool,
                "assigned_tx": str(get_tag(primary, "ZT", ".")),
                "assigned_gene": str(get_tag(primary, "ZG", ".")),
                "zf": zf_value,
                "posterior": float(get_tag(primary, "ZW", 0.0)),
                "zc": str(get_tag(primary, "ZC", ".")),
                "n_candidates": int(get_tag(primary, "ZN", 0)),
                "splice_type": str(get_tag(primary, "ZS", ".")),
                "tag_locus": int(get_tag(primary, "ZL", -1)),
                "region_class": region_class,
                "aligned_bases": aligned_bases,
                "exon_bases": exon_bases,
                "blocks": ";".join(f"{start}-{end}" for start, end in blocks),
            }
            rows.append(row)
            if origin.kind == "nrna" and example_counts[pool] < max_examples_per_pool:
                examples.append(row.copy())
                example_counts[pool] += 1

        for read in bam:
            if current_name is None:
                current_name = read.query_name
            if read.query_name != current_name:
                flush_group()
                group = []
                current_name = read.query_name
            group.append(read)
        flush_group()

    return pd.DataFrame(rows), pd.DataFrame(examples)


def add_fraction(df: pd.DataFrame, group_cols: list[str], count_col: str = "count") -> pd.DataFrame:
    df = df.copy()
    totals = df.groupby(group_cols, dropna=False)[count_col].transform("sum")
    df["fraction"] = np.where(totals > 0, df[count_col] / totals, np.nan)
    return df


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, default=DEFAULT_BASE)
    parser.add_argument("--condition", default=DEFAULT_CONDITION)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()

    base = args.base.resolve()
    condition_dir = base / args.condition
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    transcripts = pd.read_csv(base / "rigel_index" / "transcripts.tsv", sep="\t")
    nrna_quant = pd.read_csv(condition_dir / "rigel_out" / "nrna_quant.tsv", sep="\t")
    loci = pd.read_csv(condition_dir / "rigel_out" / "loci.tsv", sep="\t")
    quant = pd.read_csv(condition_dir / "rigel_out" / "quant.tsv", sep="\t")
    truth = pd.read_csv(base / "truth_abundances_nrna_high.tsv", sep="\t")
    exons_by_tx = parse_gtf_exons(base / "reference" / "genes.gtf")

    t_id_by_index = dict(zip(transcripts["t_index"].astype(int), transcripts["t_id"].astype(str)))
    origin_to_nrna_id = {}
    for row in transcripts.itertuples(index=False):
        if bool(row.is_synthetic):
            continue
        nrna_idx = int(row.nrna_t_index)
        if nrna_idx >= 0:
            origin_to_nrna_id[str(row.t_id)] = t_id_by_index.get(nrna_idx, "")
    nrna_locus_by_id = dict(zip(nrna_quant["nrna_id"].astype(str), nrna_quant["locus_id"].astype(int)))

    fragment_rows, examples = scan_annotated_bam(
        condition_dir / "annotated.bam",
        exons_by_tx,
        origin_to_nrna_id,
        nrna_locus_by_id,
    )
    fragment_rows.to_csv(out_dir / "fragment_rows.tsv", sep="\t", index=False)
    examples.to_csv(out_dir / "example_true_nrna_fragments.tsv", sep="\t", index=False)

    origin_pool = (
        fragment_rows.groupby(["origin", "assigned_pool"], dropna=False)
        .size()
        .reset_index(name="count")
    )
    origin_pool = add_fraction(origin_pool, ["origin"])
    origin_pool.to_csv(out_dir / "origin_pool_counts.tsv", sep="\t", index=False)

    true_nrna = fragment_rows[fragment_rows["origin"] == "nrna"].copy()
    nrna_by_region = (
        true_nrna.groupby(["region_class", "assigned_pool", "zc"], dropna=False)
        .size()
        .reset_index(name="count")
        .sort_values(["region_class", "assigned_pool", "zc"])
    )
    nrna_by_region = add_fraction(nrna_by_region, ["region_class"])
    nrna_by_region.to_csv(out_dir / "true_nrna_by_region_pool_zc.tsv", sep="\t", index=False)

    nrna_by_locus = (
        true_nrna.groupby(["inferred_origin_locus", "assigned_pool"], dropna=False)
        .size()
        .reset_index(name="count")
    )
    nrna_by_locus = add_fraction(nrna_by_locus, ["inferred_origin_locus"])
    nrna_by_locus = nrna_by_locus.sort_values(["count"], ascending=False)
    nrna_by_locus.to_csv(out_dir / "true_nrna_by_locus_pool.tsv", sep="\t", index=False)

    nrna_by_origin_tx = (
        true_nrna.groupby(["origin_tx", "origin_nrna_id", "inferred_origin_locus", "assigned_pool"])
        .size()
        .reset_index(name="count")
    )
    nrna_by_origin_tx = add_fraction(nrna_by_origin_tx, ["origin_tx"])
    nrna_by_origin_tx = nrna_by_origin_tx.sort_values("count", ascending=False)
    nrna_by_origin_tx.to_csv(out_dir / "true_nrna_by_origin_tx_pool.tsv", sep="\t", index=False)

    posterior_rows = []
    for keys, group in true_nrna.groupby(["region_class", "assigned_pool", "zc"], dropna=False):
        stats = summarize_numeric(group["posterior"].tolist())
        posterior_rows.append(
            {
                "region_class": keys[0],
                "assigned_pool": keys[1],
                "zc": keys[2],
                **stats,
                "median_n_candidates": float(group["n_candidates"].median()),
                "min_n_candidates": int(group["n_candidates"].min()),
                "max_n_candidates": int(group["n_candidates"].max()),
            }
        )
    posterior_summary = pd.DataFrame(posterior_rows).sort_values(
        ["region_class", "assigned_pool", "zc"]
    )
    posterior_summary.to_csv(out_dir / "true_nrna_winner_posterior_summary.tsv", sep="\t", index=False)

    top_loci = nrna_by_locus.groupby("inferred_origin_locus")["count"].sum().nlargest(10)
    locus_detail = loci[loci["locus_id"].isin(top_loci.index)].copy()
    locus_detail["true_nrna_fragments"] = locus_detail["locus_id"].map(top_loci).fillna(0).astype(int)
    locus_detail.to_csv(out_dir / "top_true_nrna_locus_outputs.tsv", sep="\t", index=False)

    top_locus_ids = set(top_loci.index.astype(int).tolist())
    nrna_quant[nrna_quant["locus_id"].isin(top_locus_ids)].to_csv(
        out_dir / "top_locus_nrna_quant.tsv",
        sep="\t",
        index=False,
    )
    quant[quant["locus_id"].isin(top_locus_ids)].merge(
        truth[["transcript_id", "mrna_abundance", "nrna_abundance"]],
        on="transcript_id",
        how="left",
    ).sort_values(["locus_id", "count"], ascending=[True, False]).to_csv(
        out_dir / "top_locus_quant.tsv",
        sep="\t",
        index=False,
    )

    print(f"wrote diagnostics to {out_dir}")
    print("origin -> assigned pool")
    print(origin_pool.to_string(index=False))
    print("\ntrue nRNA by region / assigned pool / ZC")
    print(nrna_by_region.to_string(index=False))
    print("\ntop true-nRNA loci")
    print(locus_detail.to_string(index=False))
    print("\nexample true nRNA fragments")
    print(examples.to_string(index=False))


if __name__ == "__main__":
    main()