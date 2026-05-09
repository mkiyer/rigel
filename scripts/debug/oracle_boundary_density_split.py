"""Split oracle gDNA boundary densities by boundary-side class.

This diagnostic reads the perfect-alignment ``sim_oracle.bam`` files from the
synthetic benchmark and the index ``regions.feather``. It counts true gDNA
boundary-crossing *events* directly from qname-encoded oracle fragment spans,
then estimates separate per-side densities for:

* internal exon-intron sides: EXON region has intron neighbors on both sides
* terminal exon-intron sides: EXON region has one intron neighbor and one outer side
* exon-intergenic sides: EXON side adjacent to INTERGENIC

The scanner's boundary-flux counters are event counters, not unique-fragment
counters, so this script intentionally lets one fragment contribute multiple
boundary events when it crosses multiple boundary sides.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import pysam

BASE_DEFAULT = Path("/Users/mkiyer/Downloads/rigel_runs/sim_synthetic")
REGION_TYPE_INTERGENIC = 0
REGION_TYPE_INTRON = 1
REGION_TYPE_EXON = 2

SIDE_CLASSES = (
    "exon_intron_internal",
    "exon_intron_terminal",
    "exon_intergenic",
    "exon_genome_edge",
    "exon_other",
)
SIDE_CLASS_TO_CODE = {name: i for i, name in enumerate(SIDE_CLASSES)}


@dataclass(frozen=True)
class SideIndex:
    """Sorted boundary positions and side-class codes per reference."""

    positions_by_ref: dict[str, np.ndarray]
    class_codes_by_ref: dict[str, np.ndarray]
    exon_starts_by_ref: dict[str, np.ndarray]
    exon_ends_by_ref: dict[str, np.ndarray]
    side_counts: dict[str, int]


@dataclass(frozen=True)
class OracleConditionResult:
    condition: str
    n_gdna: int
    ref_length_bp: int
    expected_rho: float
    mean_fl: float
    b_cross: float
    physical_event_counts: dict[str, int]
    observable_event_counts: dict[str, int]
    physical_densities: dict[str, float]
    observable_densities: dict[str, float]
    physical_density_se: dict[str, float]
    observable_density_se: dict[str, float]
    n_crossing_fragments: int
    n_multi_boundary_fragments: int
    n_observable_crossing_fragments: int
    n_observable_multi_boundary_fragments: int
    event_multiplicity: dict[int, int]
    rigel_rho_exon_intron: float | None
    rigel_rho_intergenic: float | None
    rigel_rho_intron: float | None


def _neighbor_type(region_df: pd.DataFrame, idx: int, side: str) -> int | None:
    row = region_df.iloc[idx]
    if side == "left":
        if idx == 0:
            return None
        neighbor = region_df.iloc[idx - 1]
        if neighbor.ref_name != row.ref_name or int(neighbor.end) != int(row.start):
            return None
        return int(neighbor.type)
    if idx == len(region_df) - 1:
        return None
    neighbor = region_df.iloc[idx + 1]
    if neighbor.ref_name != row.ref_name or int(neighbor.start) != int(row.end):
        return None
    return int(neighbor.type)


def _classify_exon_side(region_df: pd.DataFrame, idx: int, side: str) -> str:
    neighbor = _neighbor_type(region_df, idx, side)
    opposite = _neighbor_type(region_df, idx, "right" if side == "left" else "left")
    if neighbor == REGION_TYPE_INTRON:
        if opposite == REGION_TYPE_INTRON:
            return "exon_intron_internal"
        return "exon_intron_terminal"
    if neighbor == REGION_TYPE_INTERGENIC:
        return "exon_intergenic"
    if neighbor is None:
        return "exon_genome_edge"
    return "exon_other"


def build_side_index(region_df: pd.DataFrame) -> SideIndex:
    region_df = region_df.sort_values(["ref_name", "start", "end"]).reset_index(drop=True)
    side_rows: list[tuple[str, int, int, int, str]] = []
    for idx, row in enumerate(region_df.itertuples(index=False)):
        if int(row.type) != REGION_TYPE_EXON:
            continue
        side_rows.append(
            (
                str(row.ref_name),
                int(row.start),
                int(row.start),
                int(row.end),
                _classify_exon_side(region_df, idx, "left"),
            )
        )
        side_rows.append(
            (
                str(row.ref_name),
                int(row.end),
                int(row.start),
                int(row.end),
                _classify_exon_side(region_df, idx, "right"),
            )
        )

    side_df = pd.DataFrame(
        side_rows, columns=["ref", "position", "exon_start", "exon_end", "side_class"]
    )
    side_counts = {name: int((side_df["side_class"] == name).sum()) for name in SIDE_CLASSES}
    positions_by_ref: dict[str, np.ndarray] = {}
    class_codes_by_ref: dict[str, np.ndarray] = {}
    exon_starts_by_ref: dict[str, np.ndarray] = {}
    exon_ends_by_ref: dict[str, np.ndarray] = {}
    for ref_name, sub in side_df.groupby("ref", sort=False):
        sub = sub.sort_values("position")
        positions_by_ref[str(ref_name)] = sub["position"].to_numpy(np.int64)
        exon_starts_by_ref[str(ref_name)] = sub["exon_start"].to_numpy(np.int64)
        exon_ends_by_ref[str(ref_name)] = sub["exon_end"].to_numpy(np.int64)
        class_codes_by_ref[str(ref_name)] = np.array(
            [SIDE_CLASS_TO_CODE[name] for name in sub["side_class"]], dtype=np.int16
        )
    return SideIndex(
        positions_by_ref=positions_by_ref,
        class_codes_by_ref=class_codes_by_ref,
        exon_starts_by_ref=exon_starts_by_ref,
        exon_ends_by_ref=exon_ends_by_ref,
        side_counts=side_counts,
    )


def parse_gdna_qname(qname: str) -> tuple[str, int, int] | None:
    if not qname.startswith("gdna:"):
        return None
    parts = qname.split(":")
    if len(parts) < 3:
        return None
    ref_name = parts[1]
    start_s, end_s = parts[2].split("-", 1)
    return ref_name, int(start_s), int(end_s)


def iter_manifest_conditions(
    base: Path,
    *,
    include_all_ss: bool = False,
    condition_names: set[str] | None = None,
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


def load_ref_length(base: Path) -> int:
    ref_lengths = pd.read_csv(base / "rigel_index" / "ref_lengths.tsv", sep="\t")
    return int(ref_lengths["length"].sum())


def load_rigel_densities(base: Path, condition: str) -> tuple[float | None, float | None, float | None]:
    summary_path = base / condition / "rigel_out" / "summary.json"
    if not summary_path.exists():
        return None, None, None
    summary = json.loads(summary_path.read_text())
    densities = summary["calibration"]["global_densities"]
    return (
        float(densities["EXON-INTRON"]["rho"]),
        float(densities["INTERGENIC"]["rho"]),
        float(densities["INTRON"]["rho"]),
    )


def _aligned_blocks(reads: list[pysam.AlignedSegment]) -> list[tuple[int, int]]:
    blocks: list[tuple[int, int]] = []
    for read in reads:
        if read.is_unmapped:
            continue
        blocks.extend((int(start), int(end)) for start, end in read.get_blocks())
    return blocks


def _side_observed_by_read_blocks(
    *,
    exon_start: int,
    exon_end: int,
    blocks: list[tuple[int, int]],
    q: int,
) -> bool:
    for block_start, block_end in blocks:
        overlap = min(block_end, exon_end) - max(block_start, exon_start)
        if overlap >= q:
            return True
    return False


def _densities_from_counts(
    *,
    event_counts: dict[str, int],
    side_index: SideIndex,
    b_cross: float,
) -> tuple[dict[str, float], dict[str, float]]:
    densities: dict[str, float] = {}
    density_se: dict[str, float] = {}
    for name in SIDE_CLASSES:
        denom = float(side_index.side_counts[name]) * b_cross
        count = event_counts[name]
        densities[name] = float(count / denom) if denom > 0.0 else float("nan")
        density_se[name] = float(np.sqrt(count) / denom) if denom > 0.0 and count > 0 else 0.0
    return densities, density_se


def scan_oracle_condition(
    *,
    base: Path,
    condition: dict[str, object],
    side_index: SideIndex,
    ref_length_bp: int,
    splicing_anchor_tolerance: int,
) -> OracleConditionResult:
    condition_name = str(condition["name"])
    bam_path = base / str(condition["oracle_bam"])
    q = max(int(splicing_anchor_tolerance), 1)
    physical_event_counts_arr = np.zeros(len(SIDE_CLASSES), dtype=np.int64)
    observable_event_counts_arr = np.zeros(len(SIDE_CLASSES), dtype=np.int64)
    event_multiplicity: Counter[int] = Counter()
    lengths: list[int] = []
    n_gdna = 0
    n_crossing_fragments = 0
    n_multi_boundary_fragments = 0
    n_observable_crossing_fragments = 0
    n_observable_multi_boundary_fragments = 0

    current_qname: str | None = None
    current_reads: list[pysam.AlignedSegment] = []

    def process_group(qname: str, reads: list[pysam.AlignedSegment]) -> None:
        nonlocal n_gdna, n_crossing_fragments, n_multi_boundary_fragments
        nonlocal n_observable_crossing_fragments, n_observable_multi_boundary_fragments
        parsed = parse_gdna_qname(qname)
        if parsed is None:
            return
        ref_name, start, end = parsed
        frag_len = end - start
        if frag_len <= 0:
            return
        n_gdna += 1
        lengths.append(frag_len)

        positions = side_index.positions_by_ref.get(ref_name)
        if positions is None:
            return
        class_codes = side_index.class_codes_by_ref[ref_name]
        exon_starts = side_index.exon_starts_by_ref[ref_name]
        exon_ends = side_index.exon_ends_by_ref[ref_name]
        lo = start + q
        hi = end - q
        if hi < lo:
            return
        left = int(np.searchsorted(positions, lo, side="left"))
        right = int(np.searchsorted(positions, hi, side="right"))
        if right <= left:
            return

        hit_codes = class_codes[left:right]
        physical_bincount = np.bincount(hit_codes, minlength=len(SIDE_CLASSES))
        physical_event_counts_arr[:] += physical_bincount.astype(np.int64, copy=False)
        n_events = int(right - left)
        n_crossing_fragments += 1
        if n_events >= 2:
            n_multi_boundary_fragments += 1
        event_multiplicity[n_events] += 1

        blocks = _aligned_blocks(reads)
        observable_codes: list[int] = []
        for side_idx in range(left, right):
            if _side_observed_by_read_blocks(
                exon_start=int(exon_starts[side_idx]),
                exon_end=int(exon_ends[side_idx]),
                blocks=blocks,
                q=q,
            ):
                observable_codes.append(int(class_codes[side_idx]))
        if observable_codes:
            obs_bincount = np.bincount(observable_codes, minlength=len(SIDE_CLASSES))
            observable_event_counts_arr[:] += obs_bincount.astype(np.int64, copy=False)
            n_obs_events = len(observable_codes)
            n_observable_crossing_fragments += 1
            if n_obs_events >= 2:
                n_observable_multi_boundary_fragments += 1

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam:
            if read.is_secondary or read.is_supplementary:
                continue
            if current_qname is None:
                current_qname = read.query_name
            if read.query_name != current_qname:
                process_group(current_qname, current_reads)
                current_qname = read.query_name
                current_reads = []
            current_reads.append(read)
        if current_qname is not None:
            process_group(current_qname, current_reads)

    lengths_arr = np.asarray(lengths, dtype=np.float64)
    b_cross_values = np.maximum(lengths_arr - 2.0 * float(q) + 1.0, 0.0)
    b_cross = float(b_cross_values.mean()) if lengths_arr.size else 0.0
    mean_fl = float(lengths_arr.mean()) if lengths_arr.size else 0.0
    expected_rho = float(n_gdna) / float(ref_length_bp)

    physical_event_counts = {
        name: int(physical_event_counts_arr[code]) for name, code in SIDE_CLASS_TO_CODE.items()
    }
    observable_event_counts = {
        name: int(observable_event_counts_arr[code]) for name, code in SIDE_CLASS_TO_CODE.items()
    }
    physical_densities, physical_density_se = _densities_from_counts(
        event_counts=physical_event_counts, side_index=side_index, b_cross=b_cross
    )
    observable_densities, observable_density_se = _densities_from_counts(
        event_counts=observable_event_counts, side_index=side_index, b_cross=b_cross
    )

    rho_ex, rho_ig, rho_in = load_rigel_densities(base, condition_name)
    return OracleConditionResult(
        condition=condition_name,
        n_gdna=n_gdna,
        ref_length_bp=ref_length_bp,
        expected_rho=expected_rho,
        mean_fl=mean_fl,
        b_cross=b_cross,
        physical_event_counts=physical_event_counts,
        observable_event_counts=observable_event_counts,
        physical_densities=physical_densities,
        observable_densities=observable_densities,
        physical_density_se=physical_density_se,
        observable_density_se=observable_density_se,
        n_crossing_fragments=n_crossing_fragments,
        n_multi_boundary_fragments=n_multi_boundary_fragments,
        n_observable_crossing_fragments=n_observable_crossing_fragments,
        n_observable_multi_boundary_fragments=n_observable_multi_boundary_fragments,
        event_multiplicity=dict(event_multiplicity),
        rigel_rho_exon_intron=rho_ex,
        rigel_rho_intergenic=rho_ig,
        rigel_rho_intron=rho_in,
    )


def summarize_side_geometry(side_index: SideIndex) -> None:
    print("Boundary side classes from regions.feather:")
    for name in SIDE_CLASSES:
        print(f"  {name:<24s} {side_index.side_counts[name]:8d}")
    ei_sides = side_index.side_counts["exon_intron_internal"] + side_index.side_counts[
        "exon_intron_terminal"
    ]
    terminal_fraction = side_index.side_counts["exon_intron_terminal"] / ei_sides
    print(f"  exon-intron sides total   {ei_sides:8d}")
    print(f"  terminal EI side fraction {terminal_fraction:8.3%}")
    print()


def print_condition_result(result: OracleConditionResult) -> None:
    print(f"Condition: {result.condition}")
    print(
        f"  n_gdna={result.n_gdna:,}  expected_rho=n/ref={result.expected_rho:.6f}  "
        f"mean_FL={result.mean_fl:.2f}  B_cross_emp={result.b_cross:.2f}"
    )
    if result.rigel_rho_exon_intron is not None:
        print(
            f"  Rigel summary: rho_ig={result.rigel_rho_intergenic:.6f}  "
            f"rho_in={result.rigel_rho_intron:.6f}  "
            f"rho_exon_intron={result.rigel_rho_exon_intron:.6f}"
        )
    print("  Physical-fragment oracle densities by side class:")
    for name in ["exon_intron_internal", "exon_intron_terminal", "exon_intergenic"]:
        rho = result.physical_densities[name]
        se = result.physical_density_se[name]
        ratio = rho / result.expected_rho if result.expected_rho > 0.0 else float("nan")
        print(
            f"    {name:<24s} events={result.physical_event_counts[name]:8d}  "
            f"rho={rho:.6f}  ratio_to_expected={ratio:.4f}  se={se:.6f}"
        )
    print("  Scanner-observable densities by side class:")
    for name in ["exon_intron_internal", "exon_intron_terminal", "exon_intergenic"]:
        rho = result.observable_densities[name]
        se = result.observable_density_se[name]
        ratio = rho / result.expected_rho if result.expected_rho > 0.0 else float("nan")
        print(
            f"    {name:<24s} events={result.observable_event_counts[name]:8d}  "
            f"rho={rho:.6f}  ratio_to_expected={ratio:.4f}  se={se:.6f}"
        )
    all_ei_events = (
        result.physical_event_counts["exon_intron_internal"]
        + result.physical_event_counts["exon_intron_terminal"]
    )
    all_observable_ei_events = (
        result.observable_event_counts["exon_intron_internal"]
        + result.observable_event_counts["exon_intron_terminal"]
    )
    print(
        f"  Physical fragments crossing >=1 boundary side: {result.n_crossing_fragments:,}; "
        f">=2 sides: {result.n_multi_boundary_fragments:,} "
        f"({result.n_multi_boundary_fragments / max(result.n_crossing_fragments, 1):.3%})"
    )
    print(
        f"  Scanner-observable fragments crossing >=1 side: "
        f"{result.n_observable_crossing_fragments:,}; >=2 sides: "
        f"{result.n_observable_multi_boundary_fragments:,} "
        f"({result.n_observable_multi_boundary_fragments / max(result.n_observable_crossing_fragments, 1):.3%})"
    )
    print(f"  Physical exon-intron events: {all_ei_events:,}")
    print(f"  Observable exon-intron events: {all_observable_ei_events:,}")
    print()


def rows_for_output(results: list[OracleConditionResult], side_index: SideIndex) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for result in results:
        for mode, counts, densities, density_se in [
            (
                "physical_fragment",
                result.physical_event_counts,
                result.physical_densities,
                result.physical_density_se,
            ),
            (
                "scanner_observable",
                result.observable_event_counts,
                result.observable_densities,
                result.observable_density_se,
            ),
        ]:
            for name in SIDE_CLASSES:
                side_count = side_index.side_counts[name]
                if side_count == 0:
                    continue
                rho = densities[name]
                rows.append(
                    {
                        "condition": result.condition,
                        "mode": mode,
                        "side_class": name,
                        "side_count": side_count,
                        "event_count": counts[name],
                        "density": rho,
                        "density_se": density_se[name],
                        "expected_rho": result.expected_rho,
                        "ratio_to_expected": rho / result.expected_rho,
                        "mean_fl": result.mean_fl,
                        "b_cross_empirical": result.b_cross,
                        "n_gdna": result.n_gdna,
                        "rigel_rho_exon_intron": result.rigel_rho_exon_intron,
                        "rigel_rho_intergenic": result.rigel_rho_intergenic,
                        "rigel_rho_intron": result.rigel_rho_intron,
                    }
                )
    return pd.DataFrame(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, default=BASE_DEFAULT)
    parser.add_argument("--splicing-anchor-tolerance", type=int, default=0)
    parser.add_argument("--all-ss", action="store_true", help="Include SS=0.50 duplicate conditions")
    parser.add_argument("--conditions", nargs="*", default=None)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/oracle_boundary_density_split.tsv"),
        help="TSV output path, relative to the repo unless absolute.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    base = args.base.expanduser().resolve()
    condition_names = set(args.conditions) if args.conditions else None
    region_df = pd.read_feather(base / "rigel_index" / "regions.feather")
    side_index = build_side_index(region_df)
    ref_length_bp = load_ref_length(base)
    summarize_side_geometry(side_index)

    conditions = list(
        iter_manifest_conditions(base, include_all_ss=args.all_ss, condition_names=condition_names)
    )
    if not conditions:
        raise SystemExit("No positive-gDNA conditions matched the requested filters.")

    results = [
        scan_oracle_condition(
            base=base,
            condition=condition,
            side_index=side_index,
            ref_length_bp=ref_length_bp,
            splicing_anchor_tolerance=args.splicing_anchor_tolerance,
        )
        for condition in conditions
    ]
    for result in results:
        print_condition_result(result)

    output = args.output
    if not output.is_absolute():
        output = Path.cwd() / output
    output.parent.mkdir(parents=True, exist_ok=True)
    rows_for_output(results, side_index).to_csv(output, sep="\t", index=False)
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
