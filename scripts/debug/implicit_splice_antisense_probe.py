"""Interrogate implicit-splice calls and antisense overlap in synthetic oracle data.

The native BAM annotator currently stamps splice type code 3 (SPLICED_IMPLICIT)
as ``ZS=unknown``. This diagnostic therefore infers likely implicit-splice
fragments from three signals:

* ``ZS=unknown``
* at least one EM candidate (``ZN > 0`` / ``ZC != '.'``)
* no CIGAR ``N`` and no splice-blacklist marker (``ZB == 0``)

It then independently checks three geometry levels:

* strict full-intron containment in the unsequenced paired-end gap,
* any intron overlap in that gap,
* any intron overlap across the full outer genomic fragment span.

The last one matches the current native projection rule
``genomic_span - projected_tx_length > 0``. The qname provides the simulation
truth origin, so we can tabulate true-gDNA implicit-like fragments that Rigel
assigned to RNA.
"""

from __future__ import annotations

import argparse
import json
import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import pysam

BASE_DEFAULT = Path("/Users/mkiyer/Downloads/rigel_runs/sim_synthetic")
RESULTS_DEFAULT = Path("results")

AF_MRNA_BIT = 0x02
AF_GDNA_BIT = 0x04
AF_NRNA_BIT = 0x08
AF_INTERGENIC_BIT = 0x20

INTERVAL_EXON = 0
INTERVAL_TRANSCRIPT = 1
REGION_INTERGENIC = 0
REGION_INTRON = 1
REGION_EXON = 2

_QNAME_COORD_RE = re.compile(r"^[^:]+:(\d+)-(\d+):")


@dataclass(frozen=True)
class IntronIndex:
    starts_by_ref: dict[str, np.ndarray]
    ends_by_ref: dict[str, np.ndarray]
    strands_by_ref: dict[str, np.ndarray]


@dataclass(frozen=True)
class FragmentRecord:
    qname: str
    origin: str
    truth_fl: int | None
    genomic_span: int
    gap_start: int | None
    gap_end: int | None
    gap_len: int
    has_cigar_n: bool
    zf: int
    zn: int
    zb: int
    zs: str
    zc: str
    assignment: str
    contained_introns: int
    contained_intron_bp_max: int
    overlapped_introns: int
    intron_overlap_bp: int
    max_single_intron_overlap_bp: int
    span_overlapped_introns: int
    span_intron_overlap_bp: int
    span_max_single_intron_overlap_bp: int
    contained_pos: int
    contained_neg: int

    @property
    def has_contained_intron(self) -> bool:
        return self.contained_introns > 0

    @property
    def has_intron_overlap(self) -> bool:
        return self.overlapped_introns > 0

    @property
    def has_span_intron_overlap(self) -> bool:
        return self.span_overlapped_introns > 0

    @property
    def is_inferred_implicit(self) -> bool:
        return (
            self.zs == "unknown"
            and self.zn > 0
            and self.zc != "."
            and self.zb == 0
            and not self.has_cigar_n
        )

    @property
    def is_geom_candidate(self) -> bool:
        return self.zn > 0 and self.zc != "." and self.zb == 0 and not self.has_cigar_n


@dataclass(frozen=True)
class AntisenseSummary:
    metric: str
    value: float | int | str


def _get_tag(read: pysam.AlignedSegment, tag: str, default: int | str) -> int | str:
    try:
        return read.get_tag(tag)
    except KeyError:
        return default


def _read_groups(bam: pysam.AlignmentFile) -> Iterable[tuple[str, list[pysam.AlignedSegment]]]:
    current_name: str | None = None
    current_reads: list[pysam.AlignedSegment] = []
    for read in bam.fetch(until_eof=True):
        if read.is_secondary or read.is_supplementary:
            continue
        qname = read.query_name
        if current_name is None:
            current_name = qname
        if qname != current_name:
            yield current_name, current_reads
            current_name = qname
            current_reads = []
        current_reads.append(read)
    if current_name is not None:
        yield current_name, current_reads


def _truth_origin(qname: str) -> str:
    if qname.startswith("gdna:"):
        return "gdna"
    if qname.startswith("nrna:"):
        return "nrna"
    return "rna"


def _truth_fragment_length(qname: str) -> int | None:
    if qname.startswith("gdna:"):
        parts = qname.split(":")
        if len(parts) >= 3 and "-" in parts[2]:
            start_s, end_s = parts[2].split("-", 1)
            return int(end_s) - int(start_s)
        return None
    match = _QNAME_COORD_RE.match(qname)
    if match:
        return int(match.group(2)) - int(match.group(1))
    return None


def _assignment_label(zf: int) -> str:
    if zf & AF_GDNA_BIT:
        return "gdna_intergenic" if zf & AF_INTERGENIC_BIT else "gdna_em"
    if zf & AF_NRNA_BIT:
        return "nrna"
    if zf & AF_MRNA_BIT:
        return "mrna"
    return "unresolved"


def _has_cigar_n(reads: list[pysam.AlignedSegment]) -> bool:
    for read in reads:
        if read.cigartuples is None:
            continue
        if any(op == 3 for op, _length in read.cigartuples):
            return True
    return False


def _read_span(read: pysam.AlignedSegment) -> tuple[int, int] | None:
    if read.is_unmapped or read.reference_start < 0 or read.reference_end is None:
        return None
    return int(read.reference_start), int(read.reference_end)


def _paired_gap(reads: list[pysam.AlignedSegment]) -> tuple[int | None, int | None, int]:
    spans = sorted(span for read in reads if (span := _read_span(read)) is not None)
    if len(spans) < 2:
        return None, None, 0
    merged: list[list[int]] = []
    for start, end in spans:
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    if len(merged) < 2:
        return None, None, 0
    gap_start = int(merged[0][1])
    gap_end = int(merged[1][0])
    return gap_start, gap_end, max(gap_end - gap_start, 0)


def _genomic_span(reads: list[pysam.AlignedSegment], truth_fl: int | None) -> int:
    spans = [span for read in reads if (span := _read_span(read)) is not None]
    if not spans:
        return int(truth_fl or 0)
    return max(end for _start, end in spans) - min(start for start, _end in spans)


def _fragment_bounds(reads: list[pysam.AlignedSegment]) -> tuple[int | None, int | None]:
    spans = [span for read in reads if (span := _read_span(read)) is not None]
    if not spans:
        return None, None
    return min(start for start, _end in spans), max(end for _start, end in spans)


def build_intron_index(sj_df: pd.DataFrame) -> IntronIndex:
    sj_unique = sj_df[["ref", "start", "end", "strand"]].drop_duplicates()
    starts_by_ref: dict[str, np.ndarray] = {}
    ends_by_ref: dict[str, np.ndarray] = {}
    strands_by_ref: dict[str, np.ndarray] = {}
    for ref_name, sub in sj_unique.groupby("ref", sort=False):
        sub = sub.sort_values(["start", "end", "strand"])
        starts_by_ref[str(ref_name)] = sub["start"].to_numpy(np.int64)
        ends_by_ref[str(ref_name)] = sub["end"].to_numpy(np.int64)
        strands_by_ref[str(ref_name)] = sub["strand"].to_numpy(np.int8)
    return IntronIndex(starts_by_ref, ends_by_ref, strands_by_ref)


def contained_introns(
    intron_index: IntronIndex,
    ref_name: str,
    gap_start: int | None,
    gap_end: int | None,
) -> tuple[int, int, int, int]:
    if gap_start is None or gap_end is None or gap_end <= gap_start:
        return 0, 0, 0, 0
    starts = intron_index.starts_by_ref.get(ref_name)
    if starts is None:
        return 0, 0, 0, 0
    ends = intron_index.ends_by_ref[ref_name]
    strands = intron_index.strands_by_ref[ref_name]
    left = int(np.searchsorted(starts, gap_start, side="left"))
    right = int(np.searchsorted(starts, gap_end, side="left"))
    if right <= left:
        return 0, 0, 0, 0
    sub_ends = ends[left:right]
    ok = sub_ends <= gap_end
    if not np.any(ok):
        return 0, 0, 0, 0
    ok_ends = sub_ends[ok]
    ok_starts = starts[left:right][ok]
    ok_strands = strands[left:right][ok]
    lengths = ok_ends - ok_starts
    return (
        int(ok.sum()),
        int(lengths.max(initial=0)),
        int(np.count_nonzero(ok_strands == 1)),
        int(np.count_nonzero(ok_strands == 2)),
    )


def intron_overlap(
    intron_index: IntronIndex,
    ref_name: str,
    gap_start: int | None,
    gap_end: int | None,
) -> tuple[int, int, int]:
    if gap_start is None or gap_end is None or gap_end <= gap_start:
        return 0, 0, 0
    starts = intron_index.starts_by_ref.get(ref_name)
    if starts is None:
        return 0, 0, 0
    ends = intron_index.ends_by_ref[ref_name]
    right = int(np.searchsorted(starts, gap_end, side="left"))
    if right == 0:
        return 0, 0, 0
    sub_starts = starts[:right]
    sub_ends = ends[:right]
    ok = sub_ends > gap_start
    if not np.any(ok):
        return 0, 0, 0
    overlap_bp = np.minimum(sub_ends[ok], gap_end) - np.maximum(sub_starts[ok], gap_start)
    overlap_bp = overlap_bp[overlap_bp > 0]
    if len(overlap_bp) == 0:
        return 0, 0, 0
    return int(len(overlap_bp)), int(overlap_bp.sum()), int(overlap_bp.max(initial=0))


def fragment_record_from_reads(
    qname: str,
    reads: list[pysam.AlignedSegment],
    intron_index: IntronIndex,
) -> FragmentRecord | None:
    primary = [read for read in reads if not read.is_unmapped]
    if not primary:
        return None
    first = primary[0]
    ref_name = first.reference_name
    if any(read.reference_name != ref_name for read in primary):
        return None
    truth_fl = _truth_fragment_length(qname)
    gap_start, gap_end, gap_len = _paired_gap(primary)
    n_introns, max_intron_bp, n_pos, n_neg = contained_introns(
        intron_index, str(ref_name), gap_start, gap_end
    )
    n_overlapped, overlap_bp, max_overlap_bp = intron_overlap(
        intron_index, str(ref_name), gap_start, gap_end
    )
    frag_start, frag_end = _fragment_bounds(primary)
    span_n_overlapped, span_overlap_bp, span_max_overlap_bp = intron_overlap(
        intron_index, str(ref_name), frag_start, frag_end
    )
    zf = int(_get_tag(first, "ZF", 0))
    zn = int(_get_tag(first, "ZN", 0))
    zb = int(_get_tag(first, "ZB", 0))
    zs = str(_get_tag(first, "ZS", "missing"))
    zc = str(_get_tag(first, "ZC", "."))
    return FragmentRecord(
        qname=qname,
        origin=_truth_origin(qname),
        truth_fl=truth_fl,
        genomic_span=_genomic_span(primary, truth_fl),
        gap_start=gap_start,
        gap_end=gap_end,
        gap_len=gap_len,
        has_cigar_n=_has_cigar_n(primary),
        zf=zf,
        zn=zn,
        zb=zb,
        zs=zs,
        zc=zc,
        assignment=_assignment_label(zf),
        contained_introns=n_introns,
        contained_intron_bp_max=max_intron_bp,
        overlapped_introns=n_overlapped,
        intron_overlap_bp=overlap_bp,
        max_single_intron_overlap_bp=max_overlap_bp,
        span_overlapped_introns=span_n_overlapped,
        span_intron_overlap_bp=span_overlap_bp,
        span_max_single_intron_overlap_bp=span_max_overlap_bp,
        contained_pos=n_pos,
        contained_neg=n_neg,
    )


def _summary_stats(values: list[int | float]) -> dict[str, float]:
    if not values:
        return {"q50": float("nan"), "q90": float("nan"), "q99": float("nan"), "max": float("nan")}
    arr = np.asarray(values, dtype=np.float64)
    return {
        "q50": float(np.quantile(arr, 0.50)),
        "q90": float(np.quantile(arr, 0.90)),
        "q99": float(np.quantile(arr, 0.99)),
        "max": float(arr.max()),
    }


def scan_condition(
    *,
    base: Path,
    condition: dict[str, object],
    intron_index: IntronIndex,
    collect_examples: int,
) -> tuple[dict[str, object], list[dict[str, object]], list[dict[str, object]]]:
    name = str(condition["name"])
    bam_path = base / name / "annotated.bam"
    if not bam_path.exists():
        raise FileNotFoundError(bam_path)

    counters: Counter[str] = Counter()
    by_origin_assignment: Counter[tuple[str, str, str]] = Counter()
    by_origin_zs: Counter[tuple[str, str]] = Counter()
    by_origin_zc: Counter[tuple[str, str]] = Counter()
    gdna_implicit_truth_fl: list[int] = []
    gdna_implicit_genomic_span: list[int] = []
    rna_implicit_truth_fl: list[int] = []
    rna_implicit_genomic_span: list[int] = []
    gdna_contained_gap_lengths: list[int] = []
    rna_contained_gap_lengths: list[int] = []
    examples: list[dict[str, object]] = []

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for qname, reads in _read_groups(bam):
            rec = fragment_record_from_reads(qname, reads, intron_index)
            if rec is None:
                counters["skipped_unmapped_or_trans"] += 1
                continue
            counters["n_fragments"] += 1
            counters[f"origin_{rec.origin}"] += 1
            counters[f"assignment_{rec.assignment}"] += 1
            by_origin_zs[(rec.origin, rec.zs)] += 1
            by_origin_zc[(rec.origin, rec.zc)] += 1

            if rec.has_contained_intron:
                counters["n_contained_intron_any"] += 1
                counters[f"origin_{rec.origin}_contained_intron_any"] += 1
                if rec.origin == "gdna":
                    gdna_contained_gap_lengths.append(rec.gap_len)
                elif rec.origin in {"rna", "nrna"}:
                    rna_contained_gap_lengths.append(rec.gap_len)

            if rec.has_intron_overlap:
                counters["n_intron_overlap_any"] += 1
                counters[f"origin_{rec.origin}_intron_overlap_any"] += 1

            if rec.has_span_intron_overlap:
                counters["n_span_intron_overlap_any"] += 1
                counters[f"origin_{rec.origin}_span_intron_overlap_any"] += 1

            if rec.is_geom_candidate and rec.has_contained_intron:
                counters["n_geom_candidate_contained"] += 1
                counters[f"origin_{rec.origin}_geom_candidate_contained"] += 1
                if not rec.is_inferred_implicit:
                    counters["n_geom_candidate_missed"] += 1
                    counters[f"origin_{rec.origin}_geom_candidate_missed"] += 1

            if rec.is_inferred_implicit:
                counters["n_inferred_implicit"] += 1
                counters[f"origin_{rec.origin}_inferred_implicit"] += 1
                counters[f"origin_{rec.origin}_inferred_implicit_assign_{rec.assignment}"] += 1
                by_origin_assignment[(rec.origin, "implicit", rec.assignment)] += 1
                if rec.has_contained_intron:
                    counters["n_inferred_implicit_geom_ok"] += 1
                    counters[f"origin_{rec.origin}_inferred_implicit_geom_ok"] += 1
                else:
                    counters["n_inferred_implicit_no_contained_intron"] += 1
                    counters[f"origin_{rec.origin}_inferred_implicit_no_contained_intron"] += 1
                if rec.has_intron_overlap:
                    counters["n_inferred_implicit_intron_overlap_ok"] += 1
                    counters[f"origin_{rec.origin}_inferred_implicit_intron_overlap_ok"] += 1
                else:
                    counters["n_inferred_implicit_no_intron_overlap"] += 1
                    counters[f"origin_{rec.origin}_inferred_implicit_no_intron_overlap"] += 1
                if rec.has_span_intron_overlap:
                    counters["n_inferred_implicit_span_intron_overlap_ok"] += 1
                    counters[f"origin_{rec.origin}_inferred_implicit_span_intron_overlap_ok"] += 1
                else:
                    counters["n_inferred_implicit_no_span_intron_overlap"] += 1
                    counters[f"origin_{rec.origin}_inferred_implicit_no_span_intron_overlap"] += 1
                if rec.origin == "gdna":
                    if rec.assignment in {"mrna", "nrna"}:
                        counters["n_true_gdna_implicit_to_rna"] += 1
                    if rec.truth_fl is not None:
                        gdna_implicit_truth_fl.append(rec.truth_fl)
                    gdna_implicit_genomic_span.append(rec.genomic_span)
                    if len(examples) < collect_examples:
                        examples.append(_example_row(name, rec))
                elif rec.origin in {"rna", "nrna"}:
                    if rec.truth_fl is not None:
                        rna_implicit_truth_fl.append(rec.truth_fl)
                    rna_implicit_genomic_span.append(rec.genomic_span)

            if rec.zc == "ambig_opp_strand":
                counters["n_ambig_opp_strand"] += 1
                counters[f"origin_{rec.origin}_ambig_opp_strand"] += 1
                counters[f"origin_{rec.origin}_ambig_opp_assign_{rec.assignment}"] += 1

    gdna_fl = _summary_stats(gdna_implicit_truth_fl)
    gdna_span = _summary_stats(gdna_implicit_genomic_span)
    rna_fl = _summary_stats(rna_implicit_truth_fl)
    rna_span = _summary_stats(rna_implicit_genomic_span)
    gdna_gap = _summary_stats(gdna_contained_gap_lengths)
    rna_gap = _summary_stats(rna_contained_gap_lengths)

    n_gdna = max(int(counters["origin_gdna"]), 1)
    n_implicit = max(int(counters["n_inferred_implicit"]), 1)
    n_gdna_implicit = int(counters["origin_gdna_inferred_implicit"])
    n_gdna_implicit_den = max(n_gdna_implicit, 1)
    row: dict[str, object] = {
        "condition": name,
        "n_fragments": int(counters["n_fragments"]),
        "n_true_rna": int(counters["origin_rna"] + counters["origin_nrna"]),
        "n_true_gdna": int(counters["origin_gdna"]),
        "n_inferred_implicit": int(counters["n_inferred_implicit"]),
        "n_inferred_implicit_geom_ok": int(counters["n_inferred_implicit_geom_ok"]),
        "n_inferred_implicit_no_contained_intron": int(
            counters["n_inferred_implicit_no_contained_intron"]
        ),
        "inferred_implicit_geom_ok_frac": float(
            counters["n_inferred_implicit_geom_ok"] / n_implicit
        ),
        "n_inferred_implicit_intron_overlap_ok": int(
            counters["n_inferred_implicit_intron_overlap_ok"]
        ),
        "n_inferred_implicit_no_intron_overlap": int(
            counters["n_inferred_implicit_no_intron_overlap"]
        ),
        "inferred_implicit_intron_overlap_ok_frac": float(
            counters["n_inferred_implicit_intron_overlap_ok"] / n_implicit
        ),
        "n_inferred_implicit_span_intron_overlap_ok": int(
            counters["n_inferred_implicit_span_intron_overlap_ok"]
        ),
        "n_inferred_implicit_no_span_intron_overlap": int(
            counters["n_inferred_implicit_no_span_intron_overlap"]
        ),
        "inferred_implicit_span_intron_overlap_ok_frac": float(
            counters["n_inferred_implicit_span_intron_overlap_ok"] / n_implicit
        ),
        "n_geom_candidate_contained": int(counters["n_geom_candidate_contained"]),
        "n_geom_candidate_missed": int(counters["n_geom_candidate_missed"]),
        "n_true_gdna_contained_intron_any": int(counters["origin_gdna_contained_intron_any"]),
        "n_true_gdna_intron_overlap_any": int(counters["origin_gdna_intron_overlap_any"]),
        "n_true_gdna_span_intron_overlap_any": int(
            counters["origin_gdna_span_intron_overlap_any"]
        ),
        "n_true_gdna_inferred_implicit": n_gdna_implicit,
        "true_gdna_inferred_implicit_frac": float(n_gdna_implicit / n_gdna),
        "n_true_gdna_implicit_to_rna": int(counters["n_true_gdna_implicit_to_rna"]),
        "true_gdna_implicit_to_rna_frac_of_gdna": float(
            counters["n_true_gdna_implicit_to_rna"] / n_gdna
        ),
        "true_gdna_implicit_to_rna_frac_of_gdna_implicit": float(
            counters["n_true_gdna_implicit_to_rna"] / n_gdna_implicit_den
        ),
        "n_true_gdna_implicit_to_mrna": int(
            counters["origin_gdna_inferred_implicit_assign_mrna"]
        ),
        "n_true_gdna_implicit_to_nrna": int(
            counters["origin_gdna_inferred_implicit_assign_nrna"]
        ),
        "n_true_gdna_implicit_to_gdna_em": int(
            counters["origin_gdna_inferred_implicit_assign_gdna_em"]
        ),
        "n_true_gdna_implicit_to_gdna_intergenic": int(
            counters["origin_gdna_inferred_implicit_assign_gdna_intergenic"]
        ),
        "n_true_gdna_implicit_unresolved": int(
            counters["origin_gdna_inferred_implicit_assign_unresolved"]
        ),
        "n_true_rna_inferred_implicit": int(
            counters["origin_rna_inferred_implicit"] + counters["origin_nrna_inferred_implicit"]
        ),
        "n_ambig_opp_strand": int(counters["n_ambig_opp_strand"]),
        "n_true_gdna_ambig_opp_strand": int(counters["origin_gdna_ambig_opp_strand"]),
        "n_true_rna_ambig_opp_strand": int(
            counters["origin_rna_ambig_opp_strand"] + counters["origin_nrna_ambig_opp_strand"]
        ),
        "n_true_gdna_ambig_opp_to_rna": int(
            counters["origin_gdna_ambig_opp_assign_mrna"]
            + counters["origin_gdna_ambig_opp_assign_nrna"]
        ),
        "n_true_gdna_ambig_opp_to_gdna": int(
            counters["origin_gdna_ambig_opp_assign_gdna_em"]
            + counters["origin_gdna_ambig_opp_assign_gdna_intergenic"]
        ),
    }
    for prefix, stats in [
        ("gdna_implicit_truth_fl", gdna_fl),
        ("gdna_implicit_genomic_span", gdna_span),
        ("rna_implicit_truth_fl", rna_fl),
        ("rna_implicit_genomic_span", rna_span),
        ("gdna_contained_gap_len", gdna_gap),
        ("rna_contained_gap_len", rna_gap),
    ]:
        for key, value in stats.items():
            row[f"{prefix}_{key}"] = value

    long_rows: list[dict[str, object]] = []
    for (origin, zs), count in sorted(by_origin_zs.items()):
        long_rows.append({"condition": name, "axis": "ZS", "origin": origin, "label": zs, "n": count})
    for (origin, zc), count in sorted(by_origin_zc.items()):
        long_rows.append({"condition": name, "axis": "ZC", "origin": origin, "label": zc, "n": count})
    for (origin, subset, assignment), count in sorted(by_origin_assignment.items()):
        long_rows.append(
            {
                "condition": name,
                "axis": f"assignment_{subset}",
                "origin": origin,
                "label": assignment,
                "n": count,
            }
        )
    return row, long_rows, examples


def _example_row(condition: str, rec: FragmentRecord) -> dict[str, object]:
    return {
        "condition": condition,
        "qname": rec.qname,
        "truth_fl": rec.truth_fl,
        "genomic_span": rec.genomic_span,
        "gap_len": rec.gap_len,
        "contained_introns": rec.contained_introns,
        "max_contained_intron_bp": rec.contained_intron_bp_max,
        "overlapped_introns": rec.overlapped_introns,
        "intron_overlap_bp": rec.intron_overlap_bp,
        "max_single_intron_overlap_bp": rec.max_single_intron_overlap_bp,
        "span_overlapped_introns": rec.span_overlapped_introns,
        "span_intron_overlap_bp": rec.span_intron_overlap_bp,
        "span_max_single_intron_overlap_bp": rec.span_max_single_intron_overlap_bp,
        "assignment": rec.assignment,
        "zf": rec.zf,
        "zn": rec.zn,
        "zs": rec.zs,
        "zc": rec.zc,
    }


def iter_conditions(
    base: Path,
    *,
    include_all_ss: bool,
    condition_names: set[str] | None,
) -> Iterable[dict[str, object]]:
    manifest = json.loads((base / "manifest.json").read_text())
    for condition in manifest["conditions"]:
        name = str(condition["name"])
        if condition_names is not None and name not in condition_names:
            continue
        if int(condition["n_gdna"]) <= 0:
            continue
        if not include_all_ss and float(condition["strand_specificity"]) != 0.99:
            continue
        yield condition


def _merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [list(intervals[0])]
    for start, end in intervals[1:]:
        if start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(int(start), int(end)) for start, end in merged]


def _intersection_length(a: list[tuple[int, int]], b: list[tuple[int, int]]) -> int:
    i = j = 0
    total = 0
    while i < len(a) and j < len(b):
        lo = max(a[i][0], b[j][0])
        hi = min(a[i][1], b[j][1])
        if hi > lo:
            total += hi - lo
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return int(total)


def _length(intervals: list[tuple[int, int]]) -> int:
    return int(sum(end - start for start, end in intervals))


def build_strand_unions(
    intervals_df: pd.DataFrame,
    transcript_keep: set[int],
    interval_type: int,
) -> dict[str, dict[int, list[tuple[int, int]]]]:
    out: dict[str, dict[int, list[tuple[int, int]]]] = defaultdict(lambda: {1: [], 2: []})
    sub = intervals_df[
        (intervals_df["interval_type"] == interval_type)
        & (intervals_df["t_index"].isin(transcript_keep))
        & (intervals_df["strand"].isin([1, 2]))
    ]
    for row in sub.itertuples(index=False):
        out[str(row.ref)][int(row.strand)].append((int(row.start), int(row.end)))
    return {
        ref: {strand: _merge_intervals(vals) for strand, vals in by_strand.items()}
        for ref, by_strand in out.items()
    }


def antisense_summary(base: Path) -> list[AntisenseSummary]:
    t_df = pd.read_feather(base / "rigel_index" / "transcripts.feather")
    intervals = pd.read_feather(base / "rigel_index" / "intervals.feather")
    regions = pd.read_feather(base / "rigel_index" / "regions.feather")
    ref_lengths = pd.read_csv(base / "rigel_index" / "ref_lengths.tsv", sep="\t")
    genome_bp = int(ref_lengths["length"].sum())

    keep_df = t_df[(~t_df["is_nrna"].astype(bool)) & (~t_df["is_synthetic"].astype(bool))]
    transcript_keep = set(int(x) for x in keep_df["t_index"].to_numpy())

    transcript_unions = build_strand_unions(intervals, transcript_keep, INTERVAL_TRANSCRIPT)
    exon_unions = build_strand_unions(intervals, transcript_keep, INTERVAL_EXON)

    transcript_pos_bp = transcript_neg_bp = transcript_antisense_bp = 0
    exon_pos_bp = exon_neg_bp = exon_antisense_bp = 0
    for ref, by_strand in transcript_unions.items():
        pos = by_strand.get(1, [])
        neg = by_strand.get(2, [])
        transcript_pos_bp += _length(pos)
        transcript_neg_bp += _length(neg)
        transcript_antisense_bp += _intersection_length(pos, neg)
    for ref, by_strand in exon_unions.items():
        pos = by_strand.get(1, [])
        neg = by_strand.get(2, [])
        exon_pos_bp += _length(pos)
        exon_neg_bp += _length(neg)
        exon_antisense_bp += _intersection_length(pos, neg)

    opposite_pairs, transcripts_with_opposite = transcript_pair_overlap_counts(keep_df)

    region_width = regions["end"].to_numpy(np.int64) - regions["start"].to_numpy(np.int64)
    region_ambig_tx = (regions["tx_pos_bp"].to_numpy(np.int64) > 0) & (
        regions["tx_neg_bp"].to_numpy(np.int64) > 0
    )
    region_ambig_exon = (regions["exon_pos_bp"].to_numpy(np.int64) > 0) & (
        regions["exon_neg_bp"].to_numpy(np.int64) > 0
    )

    rows: list[AntisenseSummary] = [
        AntisenseSummary("n_original_transcripts", int(len(keep_df))),
        AntisenseSummary("n_original_genes", int(keep_df["g_id"].nunique())),
        AntisenseSummary("n_opposite_strand_transcript_pairs_overlapping", int(opposite_pairs)),
        AntisenseSummary("n_transcripts_with_opposite_strand_overlap", int(transcripts_with_opposite)),
        AntisenseSummary("genome_bp", genome_bp),
        AntisenseSummary("transcript_pos_union_bp", transcript_pos_bp),
        AntisenseSummary("transcript_neg_union_bp", transcript_neg_bp),
        AntisenseSummary("transcript_antisense_overlap_bp", transcript_antisense_bp),
        AntisenseSummary("transcript_antisense_overlap_frac_genome", transcript_antisense_bp / genome_bp),
        AntisenseSummary("exon_pos_union_bp", exon_pos_bp),
        AntisenseSummary("exon_neg_union_bp", exon_neg_bp),
        AntisenseSummary("exon_antisense_overlap_bp", exon_antisense_bp),
        AntisenseSummary("exon_antisense_overlap_frac_genome", exon_antisense_bp / genome_bp),
        AntisenseSummary("region_rows_with_both_tx_strands", int(region_ambig_tx.sum())),
        AntisenseSummary("region_bp_rows_with_both_tx_strands", int(region_width[region_ambig_tx].sum())),
        AntisenseSummary("region_rows_with_both_exon_strands", int(region_ambig_exon.sum())),
        AntisenseSummary("region_bp_rows_with_both_exon_strands", int(region_width[region_ambig_exon].sum())),
    ]

    for rtype, label in [
        (REGION_INTERGENIC, "intergenic"),
        (REGION_INTRON, "intron"),
        (REGION_EXON, "exon"),
    ]:
        mask = regions["type"].to_numpy(np.int64) == rtype
        rows.append(
            AntisenseSummary(
                f"region_bp_{label}_rows_with_both_tx_strands",
                int(region_width[mask & region_ambig_tx].sum()),
            )
        )
        rows.append(
            AntisenseSummary(
                f"region_bp_{label}_rows_with_both_exon_strands",
                int(region_width[mask & region_ambig_exon].sum()),
            )
        )
    return rows


def transcript_pair_overlap_counts(t_df: pd.DataFrame) -> tuple[int, int]:
    opposite_pairs = 0
    with_overlap: set[int] = set()
    for ref_name, sub in t_df.groupby("ref", sort=False):
        pos = sub[sub["strand"] == 1][["t_index", "start", "end"]].sort_values("start")
        neg = sub[sub["strand"] == 2][["t_index", "start", "end"]].sort_values("start")
        neg_starts = neg["start"].to_numpy(np.int64)
        neg_ends = neg["end"].to_numpy(np.int64)
        neg_t = neg["t_index"].to_numpy(np.int64)
        for row in pos.itertuples(index=False):
            right = int(np.searchsorted(neg_starts, int(row.end), side="left"))
            if right == 0:
                continue
            overlaps = neg_ends[:right] > int(row.start)
            n = int(np.count_nonzero(overlaps))
            if n:
                opposite_pairs += n
                with_overlap.add(int(row.t_index))
                with_overlap.update(int(x) for x in neg_t[:right][overlaps])
    return opposite_pairs, len(with_overlap)


def write_outputs(
    *,
    out_dir: Path,
    condition_rows: list[dict[str, object]],
    long_rows: list[dict[str, object]],
    examples: list[dict[str, object]],
    antisense_rows: list[AntisenseSummary],
) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(condition_rows).to_csv(
        out_dir / "implicit_splice_condition_summary.tsv", sep="\t", index=False
    )
    pd.DataFrame(long_rows).to_csv(
        out_dir / "implicit_splice_tag_breakdown.tsv", sep="\t", index=False
    )
    pd.DataFrame(examples).to_csv(
        out_dir / "implicit_splice_true_gdna_examples.tsv", sep="\t", index=False
    )
    pd.DataFrame(
        [{"metric": row.metric, "value": row.value} for row in antisense_rows]
    ).to_csv(out_dir / "antisense_overlap_summary.tsv", sep="\t", index=False)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, default=BASE_DEFAULT)
    parser.add_argument("--out-dir", type=Path, default=RESULTS_DEFAULT)
    parser.add_argument("--include-all-ss", action="store_true")
    parser.add_argument("--conditions", nargs="*", default=None)
    parser.add_argument("--examples", type=int, default=25)
    args = parser.parse_args()

    base = args.base
    condition_names = set(args.conditions) if args.conditions else None
    sj_df = pd.read_feather(base / "rigel_index" / "sj.feather")
    intron_index = build_intron_index(sj_df)

    condition_rows: list[dict[str, object]] = []
    long_rows: list[dict[str, object]] = []
    examples: list[dict[str, object]] = []
    for condition in iter_conditions(
        base,
        include_all_ss=bool(args.include_all_ss),
        condition_names=condition_names,
    ):
        row, long, ex = scan_condition(
            base=base,
            condition=condition,
            intron_index=intron_index,
            collect_examples=max(int(args.examples) - len(examples), 0),
        )
        condition_rows.append(row)
        long_rows.extend(long)
        examples.extend(ex)
        print(
            f"{row['condition']}: implicit={row['n_inferred_implicit']:,}, "
            f"true-gDNA implicit={row['n_true_gdna_inferred_implicit']:,}, "
            f"true-gDNA implicit->RNA={row['n_true_gdna_implicit_to_rna']:,}, "
            f"geom_ok={row['inferred_implicit_geom_ok_frac']:.3%}"
        )

    antisense_rows = antisense_summary(base)
    write_outputs(
        out_dir=args.out_dir,
        condition_rows=condition_rows,
        long_rows=long_rows,
        examples=examples,
        antisense_rows=antisense_rows,
    )
    print(f"Wrote outputs under {args.out_dir}")


if __name__ == "__main__":
    main()
