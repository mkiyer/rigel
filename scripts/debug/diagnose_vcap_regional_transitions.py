#!/usr/bin/env python3
"""Diagnose per-fragment VCaP pool transitions after regional exposure."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import pysam


BASE = Path("/Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m")
GDNA_FLOWCELL = "H7MFFDSXY"
RNA_FLOWCELL = "C6EL5ANXX"


def truth_source(qname: str) -> str:
    fields = qname.split(":")
    flowcell = fields[2] if len(fields) > 2 else ""
    if flowcell == GDNA_FLOWCELL:
        return "gdna"
    if flowcell == RNA_FLOWCELL:
        return "rna"
    return "unknown"


def tag(read: pysam.AlignedSegment, name: str, default=""):
    try:
        return read.get_tag(name)
    except KeyError:
        return default


def pred_pool(read: pysam.AlignedSegment) -> str:
    zf = int(tag(read, "ZF", 0))
    if zf & 0x04:
        return "gdna"
    if zf & 0x08:
        return "nrna"
    if zf & 0x02:
        return "mrna"
    return "unresolved"


def iter_primary(path: Path):
    with pysam.AlignmentFile(str(path), "rb") as bam:
        for read in bam:
            if read.is_read2 or read.is_secondary or read.is_supplementary:
                continue
            yield read


def locus_lookup(loci_path: Path) -> pd.DataFrame:
    loci = pd.read_feather(loci_path)
    keep = [
        "locus_id",
        "locus_span_bp",
        "n_transcripts",
        "n_em_fragments",
        "mrna",
        "nrna",
        "gdna",
        "gdna_eff_len",
        "gdna_eff_len_unweighted",
        "gdna_eff_len_weight_ratio",
        "gdna_prior_count",
        "gdna_prior_count_regional",
    ]
    keep = [col for col in keep if col in loci.columns]
    return loci[keep].set_index("locus_id")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--before", type=Path, default=BASE / "regional_off_v3_with_mm" / "annotated.bam")
    parser.add_argument("--after", type=Path, default=BASE / "regional_auto_v3_with_mm" / "annotated.bam")
    parser.add_argument("--after-loci", type=Path, default=BASE / "regional_auto_v3_with_mm" / "loci.feather")
    parser.add_argument("--out-dir", type=Path, default=BASE / "regional_v3_confusion")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    loci = locus_lookup(args.after_loci)
    transition_counts: Counter[tuple[str, str, str]] = Counter()
    zc_counts: Counter[tuple[str, str, str]] = Counter()
    zs_counts: Counter[tuple[str, str, str]] = Counter()
    locus_counts: Counter[int] = Counter()
    locus_transition_counts: Counter[tuple[int, str, str]] = Counter()
    zw_by_transition: defaultdict[tuple[str, str], list[float]] = defaultdict(list)

    before_iter = iter_primary(args.before)
    after_iter = iter_primary(args.after)
    n = 0
    mismatched = 0
    for before_read, after_read in zip(before_iter, after_iter):
        n += 1
        if before_read.query_name != after_read.query_name:
            mismatched += 1
            continue
        truth = truth_source(before_read.query_name)
        before_pool = pred_pool(before_read)
        after_pool = pred_pool(after_read)
        transition_counts[(truth, before_pool, after_pool)] += 1

        if truth == "gdna" and before_pool == "gdna" and after_pool in {"mrna", "nrna"}:
            zc_counts[(before_pool, after_pool, str(tag(after_read, "ZC", "missing")))] += 1
            zs_counts[(before_pool, after_pool, str(tag(after_read, "ZS", "missing")))] += 1
            lid = int(tag(after_read, "ZL", -1))
            if lid >= 0:
                locus_counts[lid] += 1
                locus_transition_counts[(lid, before_pool, after_pool)] += 1
            try:
                zw_by_transition[(before_pool, after_pool)].append(float(tag(after_read, "ZW", np.nan)))
            except (TypeError, ValueError):
                pass

    transition_rows = [
        {"truth": t, "before": b, "after": a, "count": c}
        for (t, b, a), c in transition_counts.items()
    ]
    transitions = pd.DataFrame(transition_rows).sort_values(["truth", "count"], ascending=[True, False])
    transitions.to_csv(args.out_dir / "pool_transition_counts.tsv", sep="\t", index=False)

    zc = pd.DataFrame(
        [
            {"before": b, "after": a, "ZC": zc_value, "count": c}
            for (b, a, zc_value), c in zc_counts.items()
        ]
    ).sort_values("count", ascending=False)
    zc.to_csv(args.out_dir / "gdna_before_correct_after_rna_by_zc.tsv", sep="\t", index=False)
    zs = pd.DataFrame(
        [
            {"before": b, "after": a, "ZS": zs_value, "count": c}
            for (b, a, zs_value), c in zs_counts.items()
        ]
    ).sort_values("count", ascending=False)
    zs.to_csv(args.out_dir / "gdna_before_correct_after_rna_by_zs.tsv", sep="\t", index=False)

    locus_rows = []
    for lid, count in locus_counts.most_common():
        row = {"locus_id": lid, "new_false_gdna_to_rna": count}
        if lid in loci.index:
            row.update(loci.loc[lid].to_dict())
        row["new_false_to_mrna"] = locus_transition_counts[(lid, "gdna", "mrna")]
        row["new_false_to_nrna"] = locus_transition_counts[(lid, "gdna", "nrna")]
        locus_rows.append(row)
    locus_df = pd.DataFrame(locus_rows)
    locus_df.to_csv(args.out_dir / "gdna_before_correct_after_rna_by_locus.tsv", sep="\t", index=False)

    zw_rows = []
    for transition, values in zw_by_transition.items():
        arr = np.asarray(values, dtype=np.float64)
        zw_rows.append({
            "before": transition[0],
            "after": transition[1],
            "n": int(arr.size),
            "mean": float(np.mean(arr)),
            "median": float(np.median(arr)),
            "p05": float(np.quantile(arr, 0.05)),
            "p95": float(np.quantile(arr, 0.95)),
        })
    pd.DataFrame(zw_rows).to_csv(args.out_dir / "gdna_before_correct_after_rna_zw.tsv", sep="\t", index=False)

    print(f"primary read1 compared: {n:,}; mismatched order: {mismatched:,}")
    print("top transitions")
    print(transitions.head(20).to_string(index=False))
    print("top new false loci")
    print(locus_df.head(20).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())