"""Run the reference accumulator on a REAL BAM against a REAL index, and cross-check it.

    Reference: ``tests/native/_accumulator_reference.py`` · (S2)

⚠ WHY. The spec matrix exercises the reference on hand-built fixtures of two to six regions. That cannot
catch a defect that only appears at 1.04 M regions and 404,168 sj — a per-reference offset that
drifts, an index space that wraps, a sj lookup that silently misses. The native accumulator is
about to be gated on byte-identity to this file, so the file has to be correct on real data first.

WHAT IS CHECKED. Every number the reference reports is re-derived here **by a different method** and
compared. The reference locates crossed boundaries with two ``searchsorted`` calls per segment and takes the
resulting index range; this script walks each segment's region_bounds with ``bisect`` and counts them one at a
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
    Accumulator,
    DepositOutcome,
    FragmentPool,
    Partition,
    _normalise_introns,
    _segments,
)

from rigel.calibration.splice_graph import (  # noqa: E402
    EDGE_KIND_SJ,
    build_sj_arrays,
    build_region_partition_arrays,
)
from rigel.index import TranscriptIndex  # noqa: E402
from rigel.types import Strand  # noqa: E402


def build_partition(index) -> Partition:
    """The real index, in the reference's own terms."""
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
    """Region (``per_ref=1``) or contiguous-boundary (``per_ref=2``) CSR offsets from the region_bound offsets.

    A reference contributing ``c`` region_bounds owns ``c − 1`` regions and ``c − 2`` boundaries; one contributing none
    owns neither, which is why the subtraction is clamped at zero rather than applied blindly.
    """
    counts = np.diff(region_bound_offsets.astype(np.int64))
    sizes = np.maximum(counts - per_ref, 0) * (counts > 0)
    out = np.zeros(region_bound_offsets.shape[0], np.int64)
    np.cumsum(sizes, out=out[1:])
    return out


def fragment_paths(bam: str, name_to_ref_id: dict[str, int], limit: int | None):
    """Stream ``(ref_id, lo, hi, introns, strand, sj_strand)`` from a name-collated BAM.

    The path is the design's: blocks joined across the mate gap, broken at CIGAR ``N``. Introns are
    de-duplicated on ``(start, end)`` — ⚠ the scanner reads the ``XS`` tag once per RECORD, so a pair
    where read 1 carries it and read 2 does not yields the same intron twice, and crediting both would
    double-count the sj.
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


#: ⛔ The battery ``--self-check`` runs. Every case is an intron pathology on which the naive
#: ``(hi − lo) − Σ intron`` disagrees with the true segment total, or on which filtering an intron
#: disagrees with clipping it. Four of the six caught the harness itself on 2026-07-29.
_SELF_CHECK_REGION_BOUNDS = [0, 100, 200, 300, 400, 500, 600, 1000]
_SELF_CHECK_CASES = [
    ("overlapping introns (the MO_3021 chr8 case)", 150, 500, [(210, 260), (240, 300)]),
    ("nested introns", 150, 500, [(200, 400), (250, 300)]),
    ("wide nested — the naive L goes NEGATIVE", 150, 500, [(150, 480), (160, 470)]),
    ("abutting introns — malformed, merge", 150, 500, [(200, 300), (300, 400)]),
    ("intron straddling the reference end", 900, 1000, [(950, 1200)]),
    ("clean disjoint introns (the control)", 150, 500, [(200, 250), (300, 350)]),
]


def self_check() -> int:
    """Cross-check this harness's own re-derivation against the reference, on intron pathologies.

    ⚠ This exists because the harness is the ONLY gate the native accumulator has, and it was wrong. A
    re-derivation that shares a bug with the thing it checks is worse than no re-derivation, because it
    reports agreement. Real data cannot catch it: the discriminating fragment is ~1 in 875,670, so the
    prefix a ``--limit`` run sees is very likely clean.
    """
    types = [[0] * (len(_SELF_CHECK_REGION_BOUNDS) - 1)]
    partition = Partition.from_region_bounds([_SELF_CHECK_REGION_BOUNDS], region_types=types)
    print("SELF-CHECK — this harness's re-derivation vs the reference, on intron pathologies\n")
    print(f"  {'case':44s} {'reference':>10s} {'harness':>8s} {'absorbed':>9s}")
    failures = 0
    for label, lo, hi, introns in _SELF_CHECK_CASES:
        # the reference's L, through the reference's own definition: clip, normalise, sum the segments
        clipped_lo, clipped_hi = max(lo, _SELF_CHECK_REGION_BOUNDS[0]), min(hi, _SELF_CHECK_REGION_BOUNDS[-1])
        merged, absorbed = _normalise_introns(introns, clipped_lo, clipped_hi)
        reference_length = sum(b - a for a, b in _segments(clipped_lo, clipped_hi, merged))
        _, harness_length, _ = _expected(partition, 0, lo, hi, sorted(introns), Strand.NONE)
        ok = reference_length == harness_length
        failures += not ok
        print(
            f"  {label:44s} {reference_length:>10d} {harness_length:>8d} {absorbed:>9d}"
            f"{'' if ok else '   *** DISAGREE ***'}"
        )
    print(
        f"\n  {len(_SELF_CHECK_CASES) - failures}/{len(_SELF_CHECK_CASES)} agree"
        + ("" if failures else " — the harness may be trusted on real data")
    )
    return failures


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("index", nargs="?")
    ap.add_argument("bam", nargs="?")
    ap.add_argument("--limit", type=int, default=200_000, help="fragments to process (0 = all)")
    ap.add_argument("--max-fragment-length", type=int, default=1000)
    ap.add_argument(
        "--self-check",
        action="store_true",
        help="run only the intron-pathology battery against the reference; no index or BAM needed",
    )
    args = ap.parse_args()

    if args.self_check:
        raise SystemExit(1 if self_check() else 0)
    if not args.index or not args.bam:
        ap.error("index and bam are required unless --self-check is given")

    # ⛔ The harness checks itself BEFORE it is used to check anything else.
    if self_check():
        raise SystemExit("self-check failed: this harness cannot be trusted to gate the reference")
    print()

    index = TranscriptIndex.load(args.index)
    partition = build_partition(index)
    print(f"index      {args.index}")
    print(
        f"partition  {partition.n_regions:,} regions  {partition.n_boundaries:,} contiguous boundaries  "
        f"{partition.n_sj:,} sj boundaries  {partition.region_bounds.size:,} region_bounds"
    )
    print(f"bam        {args.bam}\n")

    name_to_ref_id = {name: i for i, name in enumerate(index.ref_names)}
    acc = Accumulator(partition, max_fragment_length=args.max_fragment_length)

    # independent expectations, re-derived per fragment by a DIFFERENT method (bisect, one region_bound at a time)
    expect_crossings = expect_density = 0
    expect_sj = 0
    accepted = 0
    lengths: list[int] = []

    t0 = time.perf_counter()
    for ref_id, lo, hi, introns, align, motif in fragment_paths(
        args.bam, name_to_ref_id, args.limit or None
    ):
        outcome = acc.deposit(
            ref_id, lo, hi, introns=introns, align_strand=align, sj_strand=motif
        )
        if outcome is not DepositOutcome.DEPOSITED:
            continue
        accepted += 1
        crossings, length, sj = _expected(partition, ref_id, lo, hi, introns, motif)
        expect_crossings += crossings
        # ⛔ `inv_length_quantum(p)` was `round(2^32 / p)` — the fixed-point spelling of `1/p`. The
        # layer went at `94d283c0`; the deposit is the reciprocal itself, in float64.
        expect_density += crossings / (length - 1) if length >= 2 else 0
        expect_sj += sj
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
        ("start_count total == accepted", int(t.region_start_count.sum()), accepted),
        (
            "boundary crossings",
            int(t.boundary_unspliced_count.sum()) + int(t.boundary_spliced_count.sum()),
            expect_crossings,
        ),
        (
            "boundary density",
            int(t.boundary_unspliced_inv_length_sum.sum()) + int(t.boundary_spliced_inv_length_sum.sum()),
            expect_density,
        ),
        ("sj crossings", int(t.sj_count.sum()), expect_sj),
    ]
    ok = True
    for label, got, want in checks:
        good = got == want
        ok &= good
        print(f"  {'OK ' if good else 'FAIL'}  {label:<34} {got:>16,}  vs {want:>16,}")

    print("\nWHAT LANDED")
    contained = int(t.region_contained_count.sum())
    spanning = int(t.region_spanning_count.sum())
    print(
        f"  contained            {contained:>12,}  ({100 * contained / max(accepted, 1):5.1f} % of accepted)"
    )
    print(f"  spanning             {spanning:>12,}")
    print(f"  regions with any count {int((t.region_contained_count.sum(1) > 0).sum()):>12,}")
    print(f"  boundaries with any count {int((t.boundary_unspliced_count.sum(1) > 0).sum()):>12,}")
    print(f"  sj used       {int((t.sj_count.sum(1) > 0).sum()):>12,}")

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
    # ⛔ FLOAT, NOT INT, AND NO SCALE FACTOR. `inv_length_sum` was a FIXED-POINT integer bank divided by
    # `INV_LENGTH_SCALE` (2^32); the whole fixed-point layer went at `94d283c0` under ONE NUMERIC
    # CONVENTION, so the bank is float64 in real units. `int()` on it now TRUNCATES a density toward
    # zero — a sum of 1/L terms is < 1 on any short region — so the cast is as dead as the scale.
    density_total = float(t.boundary_unspliced_inv_length_sum.sum()) + float(
        t.boundary_spliced_inv_length_sum.sum()
    )
    count_total = int(t.boundary_unspliced_count.sum()) + int(t.boundary_spliced_count.sum())
    if density_total:
        print(
            f"pooled count/density + 1 = {count_total / density_total + 1:.1f}  "
            f"(the crossing population's own mean L, weighted by boundaries crossed)"
        )

    raise SystemExit(0 if ok else 1)


def _expected(partition, ref_id, lo, hi, introns, motif):
    """Crossings, L and annotated-sj count for one fragment, by bisect rather than searchsorted.

    ⛔ ``length`` is the TOTAL OF THE SEGMENTS, never ``(hi − lo) − Σ intron``. The two differ the moment
    two introns overlap, and this function used to use the naive form — reproducing, inside the harness,
    the exact bug the reference was fixed for (TRAPS: measure-the-ceiling-first). It agreed on real data only because a discriminating
    fragment is roughly 1 in 875,670. ``--self-check`` now pins the six cases where it matters.

    ⚠ Introns are CLIPPED to ``[lo, hi)``, not filtered by it. Filtering drops an intron that straddles
    the reference end instead of truncating it, which is a second way to disagree with the reference — the
    two behaviours differ by 50 bp on the ``[900,1000)`` + intron ``[950,1200)`` case.
    """
    first, last = int(partition.ref_region_bound_offsets[ref_id]), int(partition.ref_region_bound_offsets[ref_id + 1])
    region_bounds = partition.region_bounds[first:last].tolist()
    lo, hi = max(lo, region_bounds[0]), min(hi, region_bounds[-1])
    introns = [(max(s, lo), min(e, hi)) for s, e in introns]
    introns = sorted((s, e) for s, e in introns if s < e)

    crossings, cursor = 0, lo
    segments = []
    for s, e in introns:
        if s > cursor:
            segments.append((cursor, s))
        cursor = max(cursor, e)
    if hi > cursor:
        segments.append((cursor, hi))
    length = sum(b - a for a, b in segments)
    for a, b in segments:
        i = bisect.bisect_right(region_bounds, a)
        while i < len(region_bounds) and region_bounds[i] < b:  # one region_bound at a time, no range arithmetic
            crossings += 1
            i += 1

    sj = 0
    for s, e in introns:
        d = bisect.bisect_left(region_bounds, s)
        a = bisect.bisect_left(region_bounds, e)
        if d >= len(region_bounds) or region_bounds[d] != s or a >= len(region_bounds) or region_bounds[a] != e:
            continue
        for k in range(
            int(partition.sj_offsets[first + d]),
            int(partition.sj_offsets[first + d + 1]),
        ):
            if int(partition.sj_boundary_right[k]) != first + a:
                continue
            if motif != Strand.NONE and int(partition.sj_strand[k]) != int(motif):
                continue
            sj += 1
            break
    return crossings, length, sj


if __name__ == "__main__":
    main()
