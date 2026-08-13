"""THE S3 GATE ON REAL DATA — the native accumulator against the specification, at full human scale.

    Gate:  ``tests/native/test_accumulator_native_parity.py`` (the same comparison, on fixtures)
    Spec:  ``tests/native/_accumulator_reference.py``
     step 6

⚠ WHY THIS EXISTS SEPARATELY FROM THE UNIT GATE. The unit gate drives a seven-cut partition. It cannot see
a defect that only appears at **1,043,881 nodes and 404,168 junctions**: a per-reference offset that drifts,
a cut index that wraps at int32, a junction CSR slice that is off by one reference. Every such bug this
project has had was invisible on a fixture — a ref-id mismatch once dropped **476,719 of 476,732 fragments
inside deposit()** while every golden test passed.

⭐ AND THE SLICING IS RE-DERIVED HERE BY A DIFFERENT ROUTE, DELIBERATELY. ``AccumulatorSet::set_junctions``
slices the flat junction CSR per reference in C++; this script does the same arithmetic in numpy and feeds
the result to a per-reference ``Accumulator``. If the two disagree the junction banks disagree, which is
exactly the class of error a validator that called the builder's own helper would report as agreement.

    OMP_NUM_THREADS=1 python scripts/design/native_parity_on_real_data.py INDEX BAM [--limit N]

⚠ The reference costs ~500 us/fragment (pure Python, and that is expected of a specification), so
``--limit`` is how this stays runnable. Raise it until the numbers stop moving.
"""

from __future__ import annotations

import argparse
import dataclasses
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tests"))

from native._accumulator_reference import (  # noqa: E402
    Accumulator as ReferenceAccumulator,
    Tally,
)
from reference_on_real_data import build_partition, fragment_paths  # noqa: E402

from rigel._bam_impl import Accumulator as NativeAccumulator  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402


#: Which flat axis each Tally array is indexed by. This is the whole reason the comparison needs care: the
#: reference holds one flat array per quantity across every reference, the native class holds one
#: accumulator per reference, and the three axes have different per-reference lengths.
_AXIS = {
    "node_contained_count": "node",
    "node_contained_inv_opportunity_sum": "node",
    "node_spanning_count": "node",
    "node_spanning_inv_length_sum": "node",
    "node_spanning_length_sum": "node",
    "node_start_count": "node",
    "edge_unspliced_count": "edge",
    "edge_unspliced_inv_length_sum": "edge",
    "edge_spliced_count": "edge",
    "edge_spliced_inv_length_sum": "edge",
    "edge_spliced_length_sum": "edge",
    "sj_count": "sj",
    "sj_inv_length_sum": "sj",
    "sj_mass": "sj",
    "sj_length_sum": "sj",
}


def ref_sj_offsets(partition) -> np.ndarray:
    """Per-reference offsets into the junction axis, from the CSR alone.

    The CSR is keyed by the flat donor cut index and references are cut-major, so a reference's junctions
    are the contiguous slot range ``[sj_offsets[c0], sj_offsets[c1])``. That is also why the junction-edge
    id can BE the slot: the flat slot order is already the per-reference banks concatenated in order.
    """
    return partition.sj_offsets[partition.ref_cut_offsets]


def native_for_ref(partition, ref: int, max_length: int) -> NativeAccumulator:
    """One native accumulator for reference ``ref``, with the junction CSR sliced and rebased."""
    c0, c1 = int(partition.ref_cut_offsets[ref]), int(partition.ref_cut_offsets[ref + 1])
    n0, n1 = int(partition.ref_node_offsets[ref]), int(partition.ref_node_offsets[ref + 1])
    accumulator = NativeAccumulator(
        cuts=np.ascontiguousarray(partition.cut_positions[c0:c1], dtype=np.int64),
        node_types=np.ascontiguousarray(partition.node_types[n0:n1], dtype=np.uint8),
        max_length=max_length,
    )
    j0, j1 = int(partition.sj_offsets[c0]), int(partition.sj_offsets[c1])
    accumulator.set_junctions(
        np.ascontiguousarray(partition.sj_offsets[c0 : c1 + 1] - j0, dtype=np.int32),
        np.ascontiguousarray(partition.sj_boundary_right[j0:j1] - c0, dtype=np.int32),
        np.ascontiguousarray(partition.sj_strand[j0:j1], dtype=np.int8),
    )
    return accumulator


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("index")
    ap.add_argument("bam")
    ap.add_argument("--limit", type=int, default=50_000, help="fragments to process (0 = all)")
    ap.add_argument("--max-fragment-length", type=int, default=1000)
    args = ap.parse_args()

    index = TranscriptIndex.load(args.index)
    partition = build_partition(index)
    sj_offsets = ref_sj_offsets(partition)
    print(f"index      {args.index}")
    print(
        f"partition  {partition.n_nodes:,} nodes  {partition.n_edges:,} contiguous edges  "
        f"{partition.n_sj:,} junction edges  {partition.cut_positions.size:,} cuts"
    )
    print(f"bam        {args.bam}\n")

    name_to_ref_id = {name: i for i, name in enumerate(index.ref_names)}
    reference = ReferenceAccumulator(partition, max_fragment_length=args.max_fragment_length)
    natives: dict[int, NativeAccumulator] = {}
    outcomes: dict[str, int] = defaultdict(int)
    disagreements = 0
    n = 0

    t0 = time.perf_counter()
    for ref_id, start, end, introns, align_strand, sj_strand in fragment_paths(
        args.bam, name_to_ref_id, args.limit or None
    ):
        want = reference.deposit(
            ref_id, start, end, introns=introns, align_strand=align_strand, sj_strand=sj_strand
        )
        if ref_id not in natives:
            natives[ref_id] = native_for_ref(partition, ref_id, args.max_fragment_length)
        got = natives[ref_id].deposit(
            start=start,
            end=end,
            introns=introns,
            align_strand=align_strand,
            sj_strand=sj_strand,
        )
        outcomes[want.value] += 1
        n += 1
        if got != want.value:
            disagreements += 1
            if disagreements <= 5:
                print(
                    f"  OUTCOME DISAGREES  ref={ref_id} [{start},{end}) introns={introns} "
                    f"align={align_strand} sj={sj_strand}: native {got!r} vs reference {want.value!r}"
                )
    elapsed = time.perf_counter() - t0

    print(
        f"processed {n:,} fragments over {len(natives)} references in {elapsed:.1f}s "
        f"({1e6 * elapsed / max(n, 1):.0f} us/fragment, dominated by the pure-Python reference)\n"
    )
    print("OUTCOMES")
    for key, value in sorted(outcomes.items()):
        print(f"  {key:<28} {value:>12,}")

    # ── the comparison ────────────────────────────────────────────────────────────────────────────────
    print("\nBYTE-IDENTITY, per Tally field, over every reference that received a fragment")
    failures = []
    for field in dataclasses.fields(Tally):
        # `qc` and `pool_lengths` are library-wide, not per-object, so they are summed below rather than
        # sliced by reference. Anything else missing from _AXIS is a new Tally field with no comparison —
        # which must fail loudly rather than drop out of the gate.
        if field.name in ("qc", "pool_lengths"):
            continue
        if field.name not in _AXIS:
            raise SystemExit(
                f"Tally field {field.name!r} has no axis in _AXIS, so this gate would silently skip it. "
                f"Add it (or add it to the library-wide list) before trusting this run."
            )
        expected_flat = getattr(reference.tally, field.name)
        axis = _AXIS[field.name]
        bad = 0
        total = 0
        for ref_id, native in sorted(natives.items()):
            actual = getattr(native, field.name)
            if axis == "node":
                lo, hi = partition.ref_node_offsets[ref_id], partition.ref_node_offsets[ref_id + 1]
            elif axis == "edge":
                lo, hi = partition.ref_edge_offsets[ref_id], partition.ref_edge_offsets[ref_id + 1]
            else:
                lo, hi = sj_offsets[ref_id], sj_offsets[ref_id + 1]
            expected = expected_flat[int(lo) : int(hi)]
            if actual.shape != expected.shape:
                failures.append(f"{field.name} ref {ref_id}: shape {actual.shape} != {expected.shape}")
                continue
            if actual.dtype != expected.dtype:
                failures.append(f"{field.name} ref {ref_id}: dtype {actual.dtype} != {expected.dtype}")
                continue
            bad += int(np.count_nonzero(np.asarray(actual) != np.asarray(expected)))
            total += int(np.asarray(expected).size)
        status = "OK  " if bad == 0 else "FAIL"
        if bad:
            failures.append(f"{field.name}: {bad:,} cells differ")
        print(f"  {status}  {field.name:<24} {total:>14,} cells   nonzero {int((expected_flat != 0).sum()):>10,}")

    # The pools are library-wide, so they are summed over references rather than sliced.
    pools = np.zeros_like(reference.tally.pool_lengths)
    for native in natives.values():
        pools += np.asarray(native.pool_lengths)
    pools_ok = np.array_equal(pools, reference.tally.pool_lengths)
    print(f"  {'OK  ' if pools_ok else 'FAIL'}  pool_lengths (summed over references)")
    if not pools_ok:
        failures.append("pool_lengths differ")

    # QC counters likewise sum over references — except the ones the reference charges before it knows
    # which reference the fragment is on, which is why they are compared as a total and not per reference.
    qc = {key: 0 for key in reference.tally.qc}
    for native in natives.values():
        for key, value in native.qc.items():
            qc[key] += value
    print("\nQC DENOMINATORS")
    for key in sorted(qc):
        mark = "OK  " if qc[key] == reference.tally.qc[key] else "FAIL"
        print(f"  {mark}  {key:<28} native {qc[key]:>12,}   reference {reference.tally.qc[key]:>12,}")
        if qc[key] != reference.tally.qc[key]:
            failures.append(f"qc[{key}]: {qc[key]} != {reference.tally.qc[key]}")

    print()
    if disagreements:
        failures.append(f"{disagreements:,} per-fragment outcome disagreements")
    if failures:
        print("⛔ NOT BYTE-IDENTICAL")
        for line in failures[:20]:
            print(f"  {line}")
        raise SystemExit(1)
    print(f"✅ BYTE-IDENTICAL on {n:,} real fragments across {len(natives)} references")


if __name__ == "__main__":
    main()
