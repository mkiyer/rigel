"""Run the reference accumulator on a REAL BAM against a REAL index, and cross-check it.

    Reference: ``tests/native/_accumulator_reference.py``   ·   Plan: ``docs/IMPLEMENTATION_PLAN.md`` §5 (S2)

⚠ WHY. The spec matrix exercises the reference on hand-built fixtures of two to six nodes. That cannot
catch a defect that only appears at 1.04 M nodes and 404,168 junctions — a per-reference offset that
drifts, an index space that wraps, a junction lookup that silently misses. The native accumulator is
about to be gated on byte-identity to this file, so the file has to be correct on real data first.

WHAT IS CHECKED. Every number the reference reports is re-derived here **by a different method** and
compared. The reference locates crossed lines with two ``searchsorted`` calls per segment and takes the
resulting index range; this script walks each segment's cuts with ``bisect`` and counts them one at a
time. Agreement between the two is not automatic — that is the point.

    OMP_NUM_THREADS=1 python scripts/design/reference_on_real_data.py INDEX BAM [--limit N]
"""

from __future__ import annotations

import argparse
import bisect
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pysam

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tests"))

from native._accumulator_reference import (  # noqa: E402
    DENSITY_SCALE,
    Accumulator,
    DepositOutcome,
    FragmentPool,
    Partition,
    density_quantum,
)

from rigel.calibration.splice_graph import (  # noqa: E402
    EDGE_KIND_JUNCTION,
    build_junction_edge_arrays,
    build_node_partition_arrays,
)
from rigel.index import TranscriptIndex  # noqa: E402
from rigel.types import Strand  # noqa: E402


def build_partition(index) -> Partition:
    """The real index, in the reference's own terms."""
    cuts, cut_offsets, node_types = build_node_partition_arrays(index)
    arrays = build_junction_edge_arrays(index)
    edges = index.edges_df
    is_junction = edges["kind"].to_numpy(np.uint8) == EDGE_KIND_JUNCTION
    n_junction = int(is_junction.sum())
    if arrays.acceptor_cut.shape[0] != n_junction:
        raise SystemExit("junction CSR disagrees with edges.feather on the junction count")
    return Partition(
        cut_positions=cuts,
        ref_cut_offsets=cut_offsets,
        node_types=node_types,
        ref_node_offsets=_offsets(cut_offsets, per_ref=1),
        ref_edge_offsets=_offsets(cut_offsets, per_ref=2),
        junction_offsets=arrays.offsets,
        junction_acceptor_cut=arrays.acceptor_cut,
        junction_strand=arrays.strand,
    )


def _offsets(cut_offsets: np.ndarray, per_ref: int) -> np.ndarray:
    """Node (``per_ref=1``) or contiguous-edge (``per_ref=2``) CSR offsets from the cut offsets.

    A reference contributing ``c`` cuts owns ``c − 1`` nodes and ``c − 2`` lines; one contributing none
    owns neither, which is why the subtraction is clamped at zero rather than applied blindly.
    """
    counts = np.diff(cut_offsets.astype(np.int64))
    sizes = np.maximum(counts - per_ref, 0) * (counts > 0)
    out = np.zeros(cut_offsets.shape[0], np.int64)
    np.cumsum(sizes, out=out[1:])
    return out


def fragment_paths(bam: str, name_to_ref_id: dict[str, int], limit: int | None):
    """Stream ``(ref_id, lo, hi, introns, align_strand, motif_strand)`` from a name-collated BAM.

    The path is the design's: blocks joined across the mate gap, broken at CIGAR ``N``. Introns are
    de-duplicated on ``(start, end)`` — ⚠ the scanner reads the ``XS`` tag once per RECORD, so a pair
    where read 1 carries it and read 2 does not yields the same intron twice, and crediting both would
    double-count the junction.
    """
    af = pysam.AlignmentFile(bam, "rb")
    group, current = [], None

    def emit(records):
        if not records:
            return None
        by_ref = defaultdict(list)
        introns_by_ref = defaultdict(set)
        motif = Strand.NONE
        reverse_r1 = None
        for r in records:
            if r.is_unmapped or r.is_secondary or r.is_supplementary:
                continue
            by_ref[r.reference_id].append((r.reference_start, r.reference_end))
            position = r.reference_start
            for op, length in r.cigartuples or []:
                if op in (0, 2, 7, 8):  # M D = X consume the reference
                    position += length
                elif op == 3:  # N — an intron
                    introns_by_ref[r.reference_id].add((position, position + length))
                    position += length
            if r.has_tag("XS"):
                tag = r.get_tag("XS")
                motif = Strand.POS if tag == "+" else Strand.NEG if tag == "-" else motif
            if r.is_read1 or reverse_r1 is None:
                reverse_r1 = r.is_reverse
        if not by_ref:
            return None
        ref_id = max(by_ref, key=lambda k: len(by_ref[k]))
        name = af.get_reference_name(ref_id)
        if name not in name_to_ref_id:
            return None
        lo = min(a for a, _ in by_ref[ref_id])
        hi = max(b for _, b in by_ref[ref_id])
        align = Strand.NEG if reverse_r1 else Strand.POS
        return (name_to_ref_id[name], lo, hi, sorted(introns_by_ref[ref_id]), align, motif)

    emitted = 0
    for record in af.fetch(until_eof=True):
        if record.query_name != current:
            path = emit(group)
            if path is not None:
                yield path
                emitted += 1
                if limit and emitted >= limit:
                    return
            group, current = [], record.query_name
        group.append(record)
    path = emit(group)
    if path is not None:
        yield path


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("index")
    ap.add_argument("bam")
    ap.add_argument("--limit", type=int, default=200_000, help="fragments to process (0 = all)")
    ap.add_argument("--max-fragment-length", type=int, default=1000)
    args = ap.parse_args()

    index = TranscriptIndex.load(args.index)
    partition = build_partition(index)
    print(f"index      {args.index}")
    print(
        f"partition  {partition.n_nodes:,} nodes  {partition.n_edges:,} contiguous edges  "
        f"{partition.n_junctions:,} junction edges  {partition.cut_positions.size:,} cuts"
    )
    print(f"bam        {args.bam}\n")

    name_to_ref_id = {name: i for i, name in enumerate(index.ref_names)}
    acc = Accumulator(partition, max_fragment_length=args.max_fragment_length)

    # independent expectations, re-derived per fragment by a DIFFERENT method (bisect, one cut at a time)
    expect_crossings = expect_density = 0
    expect_junctions = 0
    accepted = 0
    lengths: list[int] = []

    t0 = time.perf_counter()
    for ref_id, lo, hi, introns, align, motif in fragment_paths(
        args.bam, name_to_ref_id, args.limit or None
    ):
        outcome = acc.deposit(
            ref_id, lo, hi, introns=introns, align_strand=align, motif_strand=motif
        )
        if outcome is not DepositOutcome.DEPOSITED:
            continue
        accepted += 1
        crossings, length, junctions = _expected(partition, ref_id, lo, hi, introns, motif)
        expect_crossings += crossings
        expect_density += crossings * density_quantum(length - 1) if length >= 2 else 0
        expect_junctions += junctions
        lengths.append(length)
    elapsed = time.perf_counter() - t0

    t = acc.tally
    print(
        f"processed {accepted:,} accepted fragments in {elapsed:.1f}s "
        f"({1e6 * elapsed / max(accepted, 1):.1f} us/fragment, pure Python)\n"
    )

    print("QC")
    for key, value in sorted(t.qc.items()):
        print(f"  {key:<28} {value:>12,}")

    print("\nCROSS-CHECKS — the reference vs an independent re-derivation")
    checks = [
        ("start_count total == accepted", int(t.node_start_count.sum()), accepted),
        (
            "edge crossings",
            int(t.edge_unspliced_count.sum()) + int(t.edge_spliced_count.sum()),
            expect_crossings,
        ),
        (
            "edge density",
            int(t.edge_unspliced_density.sum()) + int(t.edge_spliced_density.sum()),
            expect_density,
        ),
        ("junction crossings", int(t.junction_count.sum()), expect_junctions),
    ]
    ok = True
    for label, got, want in checks:
        good = got == want
        ok &= good
        print(f"  {'OK ' if good else 'FAIL'}  {label:<34} {got:>16,}  vs {want:>16,}")

    print("\nWHAT LANDED")
    contained = int(t.node_contained_count.sum())
    spanning = int(t.node_spanning_count.sum())
    print(
        f"  contained            {contained:>12,}  ({100 * contained / max(accepted, 1):5.1f} % of accepted)"
    )
    print(f"  spanning             {spanning:>12,}")
    print(f"  nodes with any count {int((t.node_contained_count.sum(1) > 0).sum()):>12,}")
    print(f"  edges with any count {int((t.edge_unspliced_count.sum(1) > 0).sum()):>12,}")
    print(f"  junctions used       {int((t.junction_count.sum(1) > 0).sum()):>12,}")

    print("\nFRAGMENT-LENGTH POOLS (mean L, from the histogram)")
    bins = np.arange(t.pool_lengths.shape[1])
    for pool in FragmentPool:
        row = t.pool_lengths[pool]
        total = int(row.sum())
        mean = float((row * bins).sum() / total) if total else float("nan")
        print(f"  {pool.name:<22} {total:>12,}  mean L {mean:8.1f}")

    if lengths:
        arr = np.asarray(lengths)
        print(f"\nall accepted fragments: mean L {arr.mean():.1f}  median {np.median(arr):.0f}")
    density_total = int(t.edge_unspliced_density.sum()) + int(t.edge_spliced_density.sum())
    count_total = int(t.edge_unspliced_count.sum()) + int(t.edge_spliced_count.sum())
    if density_total:
        print(
            f"pooled count/density + 1 = {DENSITY_SCALE * count_total / density_total + 1:.1f}  "
            f"(the crossing population's own mean L, weighted by lines crossed)"
        )

    raise SystemExit(0 if ok else 1)


def _expected(partition, ref_id, lo, hi, introns, motif):
    """Crossings, L and annotated-junction count for one fragment, by bisect rather than searchsorted."""
    first, last = int(partition.ref_cut_offsets[ref_id]), int(partition.ref_cut_offsets[ref_id + 1])
    cuts = partition.cut_positions[first:last].tolist()
    lo, hi = max(lo, cuts[0]), min(hi, cuts[-1])
    introns = [(s, e) for s, e in introns if lo <= s < e <= hi]
    length = (hi - lo) - sum(e - s for s, e in introns)

    crossings, cursor = 0, lo
    segments = []
    for s, e in introns:
        if s > cursor:
            segments.append((cursor, s))
        cursor = max(cursor, e)
    if hi > cursor:
        segments.append((cursor, hi))
    for a, b in segments:
        i = bisect.bisect_right(cuts, a)
        while i < len(cuts) and cuts[i] < b:  # one cut at a time, no range arithmetic
            crossings += 1
            i += 1

    junctions = 0
    for s, e in introns:
        d = bisect.bisect_left(cuts, s)
        a = bisect.bisect_left(cuts, e)
        if d >= len(cuts) or cuts[d] != s or a >= len(cuts) or cuts[a] != e:
            continue
        for k in range(
            int(partition.junction_offsets[first + d]),
            int(partition.junction_offsets[first + d + 1]),
        ):
            if int(partition.junction_acceptor_cut[k]) != first + a:
                continue
            if motif != Strand.NONE and int(partition.junction_strand[k]) != int(motif):
                continue
            junctions += 1
            break
    return crossings, length, junctions


if __name__ == "__main__":
    main()
