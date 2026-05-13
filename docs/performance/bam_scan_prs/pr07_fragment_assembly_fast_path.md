# PR 07: Add A Common-Case Fragment Assembly Fast Path

## Summary

Avoid `unordered_map` and `std::set` in `build_fragment` for ordinary paired-end
fragments. Keep the current general implementation as a fallback for complex
records.

## Motivation

`build_fragment` currently handles all cases with generic containers:

- `std::unordered_map<int64_t, vector<pair<int32_t, int32_t>>>` for exons by
  `(ref_id, strand)`
- `std::set<tuple<int32_t, int32_t, int32_t, int32_t>>` for introns
- sorting and moving into `AssembledFragment`

Most fragments are ordinary paired-end alignments with one R1 location and one
R2 location, often on the same reference and with small exon/intron counts. The
generic path pays allocation and hashing/tree costs even when a stack-sized or
small-vector path would be enough.

## Current Code

- Fragment assembly: [src/rigel/native/bam_scanner.cpp](../../../src/rigel/native/bam_scanner.cpp)
- Call site: `BamScanner::process_qname_group_threaded` in the same file

## Proposed Scope

Add a fast path for the common case and fallback to the current implementation
when conditions are not met.

Candidate fast-path guard:

- `r1_reads.size() <= 1`
- `r2_reads.size() <= 1`
- at least one side exists
- no supplementary multi-record bundle in the hit
- total exon block count is small enough for stack or reserved vector storage

The fast path should:

1. Compute the effective strand for R1 and flipped R2.
2. Append exon blocks into a small local vector of `ExonBlock`.
3. Append intron blocks into a small local vector of `IntronBlock`.
4. Sort exons by `(ref_id, start, end, strand)`.
5. Merge adjacent/overlapping exons only when they share `(ref_id, strand)`.
6. Sort and unique introns.
7. Return `AssembledFragment` with the same `nm` as the generic path.

## Implementation Steps

1. Rename the current implementation to `build_fragment_generic`.
2. Add `try_build_fragment_fast` returning `std::optional<AssembledFragment>` or
   a boolean plus output parameter.
3. Implement `build_fragment` as fast path first, generic fallback second.
4. Add debug-only assertions in tests or local builds that fast and generic
   paths match for representative cases.
5. Keep behavior identical for multimappers and supplementary records by falling
   back.

## Tests

Add tests around fragment assembly if there is a suitable existing native test
hook. If not, test through scanner-level fixtures:

- ordinary paired-end unspliced fragment
- paired-end spliced fragment
- R2 strand flip behavior
- overlapping mate exons that need merging
- different refs or chimeric cases should fallback or match generic behavior
- multimapper groups should remain unchanged

Run:

```bash
conda activate rigel
pip install --no-build-isolation -e .
pytest tests/test_resolution.py tests/test_cross_chunk.py tests/test_pipeline_routing.py -v
pytest tests/test_golden_output.py -v
```

## Benchmark Plan

Run the current best no-spill scan profile:

```bash
python scripts/profiling/scan_profile.py \
  --bam /Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/annotated.bam \
  --index /Users/mkiyer/Downloads/rigel_runs/refs/rigel_index \
  --outdir /Users/mkiyer/Downloads/rigel_runs/scan_profile_pr07 \
  --name-prefix pr07 \
  --n-scan-threads 4 \
  --n-decomp-threads 2 \
  --chunk-size 1000000 \
  --max-memory-gib 12
```

If possible, add counters for fast-path hit rate and fallback reason counts.

## Acceptance Criteria

- Golden outputs are unchanged.
- VCAP scan stats are unchanged: fragment counts, intergenic counts, chimera
  counts, strand observations, and FL observations.
- Fast-path hit rate is high enough to justify the added code path.
- Native samples show lower allocator activity around `build_fragment` and
  `parse_cigar` vector movement.

## Risks

- R2 strand flip convention is easy to get wrong. Keep explicit tests.
- The generic path handles complex supplementary/multimapper bundles; do not
  force those through the fast path.
- Different ordering could affect deterministic output if not matched exactly.

## Non-Goals

- Do not change multimapper pairing.
- Do not change CIGAR parsing or splice blacklist behavior in this PR.
