# PR 10: Buffer Column Diet

**Position in roadmap:** Sixth. Touches the scan-buffer ABI; safe to
schedule in parallel with PR07/PR08 because it does not touch the
scoring CSR or partition.

## Summary

Two cleanups, in this order:

1. **Drop unused per-fragment columns.** `score_chunk()` consumes 13
   buffer columns. Four columns — `exon_bp_pos`, `exon_bp_neg`,
   `tx_bp_pos`, `tx_bp_neg` — are not in that set. If the audit
   confirms they have no production consumer, remove them from the
   chunk schema (and the spill format).
2. **Narrow safe column dtypes.** `frag_lengths`, `exon_bp`,
   `intron_bp`, `read_length` are all bounded by realistic alignment
   sizes. Move them from `int32` / `uint32` to `uint16` with explicit
   overflow guards.

## Motivation

The scan buffer logical size is ~6.8 GB on the VCAP run; 3.9 GB stays
resident at the 4 GiB cap. Lower caps are now cheap (PR06), but the
buffer itself can shrink. Removing four `int32` columns at ~32 M rows
saves ~512 MB; narrowing four more by half saves another ~512 MB. Spill
writes shrink proportionally.

## Verified ground truth

* `_FinalizedChunk` declares 21 columns in
  [src/rigel/buffer.py](../../../src/rigel/buffer.py) (~lines 138–152).
* `score_chunk()` in
  [src/rigel/native/scoring.cpp](../../../src/rigel/native/scoring.cpp)
  (~line 1128) takes a 13-tuple: `t_off, t_ind, f_len, e_bp, i_bp,
  s_type, e_str, fc, f_id, r_len, g_fp, g_sta, nm`.
* The four `_pos` / `_neg` columns are **not** in that 13-tuple.
* They are produced by the resolver and stored on the chunk; whether
  any other consumer exists must be confirmed by `grep` before
  deletion.

## Current code

* Buffer schema and spill / load: [src/rigel/buffer.py](../../../src/rigel/buffer.py).
* Native accumulator: [src/rigel/native/resolve_context.h](../../../src/rigel/native/resolve_context.h).
* Scoring ABI: [src/rigel/native/scoring.cpp](../../../src/rigel/native/scoring.cpp).

## Step 1 — Audit and drop unused columns

1. `grep -rn "exon_bp_pos\|exon_bp_neg\|tx_bp_pos\|tx_bp_neg"` across
   `src/`, `tests/`, and `scripts/`. Catalogue every consumer.
2. Classify each consumer: production (kept), test-only (kept or
   updated), diagnostic (move behind a flag), dead (removed).
3. If the only consumers are storage / spill / dead-code branches,
   remove the columns from `_FinalizedChunk`, the native accumulator,
   the zero-copy finalizer, the spill writer, and the spill loader.
4. If a real diagnostic consumer exists, gate the columns on
   `BamScanConfig.emit_strand_diagnostics: bool = False`. Default off.

This is the highest-leverage piece of the PR.

## Step 2 — Narrow safe dtypes

| Field | Current | New | Bound check |
|---|---:|---:|---|
| `frag_lengths` | int32 | uint16 | ≤ 65535 (mate-pair span) |
| `exon_bp` | int32 | uint16 | ≤ 65535 (read length) |
| `intron_bp` | int32 | uint16 | ≤ 65535 (read footprint) |
| `read_length` | uint32 | uint16 | ≤ 65535 |

Add explicit overflow guards at the *append* boundary in the native
accumulator. On overflow, raise a clear error naming the column and
the offending value. Do not silently saturate.

`frag_id` stays `int64`. `num_hits` and `nm` stay `uint16`. Do not touch
fields that already match their realistic range.

Document the dtype contract in a single comment block at the top of the
schema in `buffer.py` so future additions follow it.

## Implementation steps

1. Audit `_pos` / `_neg` consumers (Step 1).
2. Land Step 1 as commit A. Validate goldens.
3. Land Step 2 as commit B. Validate goldens.
4. Update `_FinalizedChunk.memory_bytes` and `to_scoring_arrays()`.
5. Update spill writer / loader to match. Document that spill files
   are version-local (already true in practice).
6. Update PR05's `scoring_csr.candidate_bytes` keys if applicable.
7. Update `ChunkPtrs` and nanobind casts in `scoring.cpp`. Promote
   to `int32` *locally* in scoring arithmetic where signed math
   matters; leave the storage narrow.

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

* Spill→load round-trip preserves dtypes.
* Append at the dtype boundary works (e.g. `read_length = 65535`).
* Append over the boundary raises a clear error naming the column.
* If Step 1 is gated, both `emit_strand_diagnostics=True` and
  `=False` paths produce equivalent quant outputs.

## Benchmark plan

After PR05 (for `--scan-buffer-size`):

```bash
python scripts/profiling/scan_profile.py \
  --bam /Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/annotated.bam \
  --index /Users/mkiyer/Downloads/rigel_runs/refs/rigel_index \
  --outdir /Users/mkiyer/Downloads/rigel_runs/scan_profile_pr10_dtype_diet \
  --name-prefix pr10 --threads 8 --scan-bgzf-threads 2 \
  --scan-buffer-size 2 --scan-read-name-batch-size 512
```

Then full staged. Compare buffer `memory_bytes`, spill count, spill
read/write time, peak RSS.

## Acceptance criteria

* Golden outputs unchanged.
* Buffer `memory_bytes` drops by ≥ 25% on VCAP (Step 1 + Step 2
  combined).
* Spill file size drops proportionally.
* No scan or scoring wall-time regression > 5%.
* Overflow guards have explicit tests.
* No `int32` storage remains for fields whose realistic max < 65536.

## Risks

* Real datasets may legitimately exceed `uint16` for some unusual
  column. Guards must produce clear errors, not silent corruption.
* Removing strand-aware columns may break niche diagnostics. Step 1
  exists to catch this with the audit.
* Native / Python dtype drift is the classic silent-corruption source.
  Type asserts in tests are mandatory.

## Non-goals

* No `ScoredFragments` payload changes (PR07).
* No resolver semantics changes.
* No spill format / compression changes beyond what dtype narrowing
  forces.
