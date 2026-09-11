"""Does native parity hold on real data at full human scale? The native accumulator against the
executable specification, fragment by fragment, on a real BAM.

The unit gate (`tests/native/test_accumulator_native_parity.py`) drives the same comparison on a
seven-bound fixture, and a fixture cannot see a defect that only appears at a million regions and
hundreds of thousands of sj: a per-reference offset that drifts, an index that wraps at int32, a sj CSR
slice off by one reference. This streams fragment paths from a name-collated BAM, deposits each into the
pure-Python reference (`tests/native/_accumulator_reference.py`, which wins over any document about the
accumulator) and into one native accumulator per reference, and demands byte-identity per fragment
outcome, per `Tally` field over every reference that received a fragment (each flat reference array
sliced by its axis), on the library-wide pools summed over references, and on the QC denominators. The
per-reference sj CSR slicing is re-derived in numpy rather than taken from the builder's own helper, so a
disagreement there shows in the sj banks instead of being reported as agreement. A `Tally` field with no
axis in `_AXIS` refuses to run rather than dropping out of the gate. Real data is a test input here and
never a design input. The reference is pure Python and slow, so `--limit` is how this stays runnable;
raise it until the numbers stop moving.

Usage::

    OMP_NUM_THREADS=1 python scripts/design/native_parity_on_real_data.py INDEX BAM
    OMP_NUM_THREADS=1 python scripts/design/native_parity_on_real_data.py INDEX BAM --limit 0   # every fragment
    python scripts/design/native_parity_on_real_data.py INDEX BAM --limit 200000 --max-fragment-length 1000
"""

from __future__ import annotations

import argparse
import dataclasses
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pysam

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "tests"))

from native._accumulator_reference import (  # noqa: E402
    Accumulator as ReferenceAccumulator,
    Partition,
    Tally,
)
from rigel._bam_impl import Accumulator as NativeAccumulator  # noqa: E402
from rigel.calibration.splice_graph import (  # noqa: E402
    EDGE_KIND_SJ,
    build_region_partition_arrays,
    build_sj_arrays,
)
from rigel.index import TranscriptIndex  # noqa: E402
from rigel.types import Strand  # noqa: E402


#: Which flat axis each Tally array is indexed by: the reference holds one flat array per quantity
#: across every reference, the native class holds one accumulator per reference, and the three axes
#: have different per-reference lengths.
_AXIS = {
    "region_contained_count": "region",
    "region_contained_inv_opportunity_sum": "region",
    "region_spanning_count": "region",
    "region_spanning_inv_length_sum": "region",
    "region_spanning_length_sum": "region",
    "region_start_count": "region",
    "boundary_unspliced_count": "boundary",
    "boundary_unspliced_inv_length_sum": "boundary",
    "boundary_spliced_count": "boundary",
    "boundary_spliced_inv_length_sum": "boundary",
    "boundary_spliced_length_sum": "boundary",
    "sj_count": "sj",
    "sj_inv_length_sum": "sj",
    "sj_mass": "sj",
    "sj_length_sum": "sj",
}



def build_partition(index) -> Partition:
    """The real index as the reference's ``Partition``; refuses if the sj CSR and `edges_df` disagree on the sj count."""
    region_bounds, region_bound_offsets, region_types = build_region_partition_arrays(index)
    arrays = build_sj_arrays(index)
    boundaries = index.edges_df
    is_sj = boundaries["kind"].to_numpy(np.uint8) == EDGE_KIND_SJ
    n_sj = int(is_sj.sum())
    if arrays.boundary_right.shape[0] != n_sj:
        raise SystemExit("sj CSR disagrees with edges.feather on the sj count")
    return Partition(
        region_bounds=region_bounds,
        ref_region_bound_offsets=region_bound_offsets,
        region_types=region_types,
        ref_region_offsets=_offsets(region_bound_offsets, per_ref=1),
        ref_boundary_offsets=_offsets(region_bound_offsets, per_ref=2),
        sj_offsets=arrays.offsets,
        sj_boundary_right=arrays.boundary_right,
        sj_strand=arrays.strand,
    )


def _offsets(region_bound_offsets: np.ndarray, per_ref: int) -> np.ndarray:
    """Region (``per_ref=1``) or contiguous-boundary (``per_ref=2``) CSR offsets from the region-bound offsets.

    A reference contributing ``c`` region bounds owns ``c - 1`` regions and ``c - 2`` boundaries; one
    contributing none owns neither, which is why the subtraction is clamped at zero.
    """
    counts = np.diff(region_bound_offsets.astype(np.int64))
    sizes = np.maximum(counts - per_ref, 0) * (counts > 0)
    out = np.zeros(region_bound_offsets.shape[0], np.int64)
    np.cumsum(sizes, out=out[1:])
    return out


def fragment_paths(bam: str, name_to_ref_id: dict[str, int], limit: int | None):
    """Stream ``(ref_id, lo, hi, introns, align_strand, sj_strand)`` fragment paths from a name-collated BAM.

    Blocks are joined across the mate gap and broken at CIGAR ``N``; introns are de-duplicated on
    ``(start, end)``, because a pair whose two records both carry the same ``N`` would otherwise credit
    the sj twice. Stops after ``limit`` paths when ``limit`` is given.
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


def ref_sj_offsets(partition) -> np.ndarray:
    """Per-reference offsets into the sj axis, from the CSR alone.

    The CSR is keyed by the flat left region-bound index and references are region-bound-major, so a
    reference's sj are the contiguous slot range ``[sj_offsets[c0], sj_offsets[c1])``; the flat slot order
    is already the per-reference banks concatenated in order.
    """
    return partition.sj_offsets[partition.ref_region_bound_offsets]


def native_for_ref(partition, ref: int, max_length: int) -> NativeAccumulator:
    """One native accumulator for reference ``ref``, with the sj CSR sliced and rebased."""
    c0, c1 = int(partition.ref_region_bound_offsets[ref]), int(partition.ref_region_bound_offsets[ref + 1])
    n0, n1 = int(partition.ref_region_offsets[ref]), int(partition.ref_region_offsets[ref + 1])
    accumulator = NativeAccumulator(
        region_bounds=np.ascontiguousarray(partition.region_bounds[c0:c1], dtype=np.int64),
        region_types=np.ascontiguousarray(partition.region_types[n0:n1], dtype=np.uint8),
        max_length=max_length,
    )
    j0, j1 = int(partition.sj_offsets[c0]), int(partition.sj_offsets[c1])
    accumulator.set_sj(
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
        f"partition  {partition.n_regions:,} regions  {partition.n_boundaries:,} contiguous boundaries  "
        f"{partition.n_sj:,} sj boundaries  {partition.region_bounds.size:,} region_bounds"
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
            if axis == "region":
                lo, hi = partition.ref_region_offsets[ref_id], partition.ref_region_offsets[ref_id + 1]
            elif axis == "boundary":
                lo, hi = partition.ref_boundary_offsets[ref_id], partition.ref_boundary_offsets[ref_id + 1]
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
        for boundary in failures[:20]:
            print(f"  {boundary}")
        raise SystemExit(1)
    print(f"✅ BYTE-IDENTICAL on {n:,} real fragments across {len(natives)} references")


if __name__ == "__main__":
    main()
