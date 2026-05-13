# PR 10: Buffer Dtype Diet

**Position in roadmap:** Sixth, but can run in parallel with PR07/PR08 if
carefully isolated. This PR touches scan buffer representation, not EM
payload representation.

## Summary

Reduce the logical size of `FragmentBuffer` chunks by narrowing safe
columns and removing or late-materializing columns that are no longer
consumed by the quantification path.

## Motivation

The scan buffer is still a large object. With a 4 GiB cap, the VCAP run
keeps 3.9 GB resident and spills 14 chunks. The full logical buffer is
about 6.8 GB. Lowering the cap helps RSS, but the underlying chunk format
is wider than necessary.

Narrowing columns reduces resident chunks, spill volume, Arrow read/write
time, and cache pressure during scoring.

## Current code

- Chunk schema and dtypes: [src/rigel/buffer.py](../../../src/rigel/buffer.py)
- Native accumulator fields: [src/rigel/native/resolve_context.h](../../../src/rigel/native/resolve_context.h)
- Scoring input ABI: [src/rigel/native/scoring.cpp](../../../src/rigel/native/scoring.cpp)
- Spill/load schema: [src/rigel/buffer.py](../../../src/rigel/buffer.py)

## Candidate dtype changes

Implement only after proving range bounds with assertions/tests.

| Field | Current | Candidate | Condition |
|---|---:|---:|---|
| `frag_lengths` | int32 | uint16 or int16 | `max_frag_length <= 65535`, sentinel not needed |
| `exon_bp` | int32 | uint16 | read-length bounded |
| `intron_bp` | int32 | uint16 | read/genomic footprint bounded; verify |
| `read_length` | uint32 | uint16 | aligned read length <= 65535 |
| `num_hits` | uint16 | keep | already narrow |
| `nm` | uint16 | keep | already narrow |
| `frag_id` | int64 | maybe int32 | only if annotation/writeback supports large-run fallback |

Potentially unused per-fragment columns:

- `exon_bp_pos`
- `exon_bp_neg`
- `tx_bp_pos`
- `tx_bp_neg`

These appear stored and spilled but not consumed by the current scoring
path. They may be calibration leftovers or diagnostic hooks; audit before
removal.

## Implementation steps

1. Add a schema audit script or test that summarizes min/max values for
   candidate columns on small fixtures and optionally on the VCAP run.
2. Grep and confirm consumers for `exon_bp_pos/neg` and `tx_bp_pos/neg`.
   If only storage/spill/tests reference them, remove them from the
   production chunk schema or gate them behind a diagnostics flag.
3. Narrow one group of fields at a time. Start with `read_length`,
   `exon_bp`, and `intron_bp` because they are naturally read-length
   bounded.
4. Update native accumulator vector types in
   [src/rigel/native/resolve_context.h](../../../src/rigel/native/resolve_context.h)
   and the zero-copy finalizer output.
5. Update `_FinalizedChunk.from_raw(...)`, `memory_bytes`,
   `to_scoring_arrays()`, `_spill_chunk(...)`, and `_load_chunk(...)` in
   [src/rigel/buffer.py](../../../src/rigel/buffer.py).
6. Update `ChunkPtrs` and nanobind casts in
   [src/rigel/native/scoring.cpp](../../../src/rigel/native/scoring.cpp).
   Promote to int32 locally only where arithmetic needs it.
7. Add overflow checks at append/finalize boundaries. Fail loudly if a
   value cannot fit the narrowed dtype.
8. Preserve backward-compatible spill loading only if old temporary files
   can be observed across process versions. Otherwise, document that spill
   files are ephemeral and version-local.
9. Recompute `memory_bytes` and profile JSON byte metrics.

## Tests

```bash
conda activate rigel
pip install --no-build-isolation -e .
pytest tests/test_buffer.py tests/test_resolution.py -v
pytest tests/test_pipeline_smoke.py tests/test_pipeline_routing.py -v
pytest tests/test_golden_output.py -v
pytest tests/ -q
```

Add focused tests:

- Round-trip a chunk through spill/load and assert dtypes.
- Force values at dtype boundaries.
- Force overflow and assert a clear error.
- Verify scoring arrays still produce identical routing on fixtures.

## Benchmark plan

Run scan-only and full profiles:

```bash
conda activate rigel
python scripts/profiling/scan_profile.py \
  --bam /Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/annotated.bam \
  --index /Users/mkiyer/Downloads/rigel_runs/refs/rigel_index \
  --outdir /Users/mkiyer/Downloads/rigel_runs/scan_profile_pr10_dtype_diet \
  --name-prefix pr10 --n-scan-threads 8 --n-decomp-threads 2 \
  --max-memory-gib 2 --qname-batch-size 512
```

Then run the full staged profiler. Compare buffer memory MB, spill count,
spill read/write time, peak RSS, and `fragment_router_scan`.

## Acceptance criteria

- Golden outputs unchanged.
- Buffer `memory_bytes` drops by at least 15% on VCAP; stretch target 25%.
- Spill file size drops proportionally.
- No scan or scoring wall-time regression above 5%.
- Overflow paths are tested and produce clear messages.

## Risks

- Some real datasets may exceed uint16 assumptions. Keep the checks
  explicit and consider fallback dtypes if needed.
- Removing strand-aware columns may break diagnostics not covered by the
  main pipeline. Audit thoroughly and preserve behind a flag if unsure.
- Native/Python dtype drift can create silent corruption. Add dtype tests.

## Non-goals

- Do not change `ScoredFragments` payload dtypes here; that is PR07.
- Do not alter resolver semantics or fragment classification.