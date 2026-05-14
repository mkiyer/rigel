# PR 10: Buffer Column Diet

**Status:** implemented and profiled on the VCAP RNA20M + gDNA20M workload.

## Summary

PR10 removes dead scan-buffer/scorer fields and narrows only the buffer
columns that survived real-data bounds checks.

Implemented changes:

- Dropped `intron_bp` from the scan buffer, spill format, and native
  scorer chunk ABI. The scorer accepted it positionally but never read it.
- Dropped the stale SRD-v2 strand-aware diagnostic columns from the scan
  buffer and spill format: `exon_bp_pos`, `exon_bp_neg`, `tx_bp_pos`,
  `tx_bp_neg`.
- Narrowed `exon_bp` and `read_length` to uint16 with native append guards.
- Kept `frag_lengths` as int32 after real data produced a transcript-space
  fragment length of 68,466 bp, above the uint16 limit.
- Left direct `ResolvedFragment` introspection fields intact, including
  `intron_bp` and the four strand-aware diagnostics.

The scorer chunk ABI is now a 12-tuple:

```text
t_offsets, t_indices, frag_lengths, exon_bp, splice_type, exon_strand,
fragment_classes, frag_id, read_length, genomic_footprint, genomic_start, nm
```

## Audit Results

| Field | Scope before PR10 | Audit result | PR10 action |
|---|---|---|---|
| `tx_start`, `tx_end` | Native scorer outputs removed in PR07 | Full symbol/search audit found no live downstream consumers; stale docs only | Remain removed |
| `intron_bp` | Per-candidate scan-buffer/spill/scorer ABI column | Positional scorer input, but no reads in scoring arithmetic or router construction | Removed from scan buffer, spill, and scorer ABI |
| `exon_bp_pos`, `exon_bp_neg` | Per-fragment scan-buffer/spill columns plus direct resolver fields | No production scan-buffer consumer; direct resolver/debug tests still consume public fields | Removed from scan buffer/spill only |
| `tx_bp_pos`, `tx_bp_neg` | Per-fragment scan-buffer/spill columns plus direct resolver fields | No production scan-buffer consumer; direct resolver/debug tests still consume public fields | Removed from scan buffer/spill only |
| `frag_lengths` | Per-candidate int32 buffer/scorer field | Attempted uint16 narrowing failed on VCAP at value 68,466 | Kept int32 |
| `exon_bp` | Per-candidate int32 buffer/scorer field | Bounded by aligned read overlap; native guard covers overflow | Narrowed to uint16 |
| `read_length` | Per-fragment uint32 buffer/scorer field | Aligned block span; native guard covers overflow | Narrowed to uint16 |
| `sj_strand`, `merge_criteria` | Per-fragment storage/debug fields | Not hot memory drivers; part of `BufferedFragment` debug/duck-typed surface | Kept |
| `frag_id` | Per-fragment int64 | Must support large BAMs | Kept int64 |
| `num_hits`, `nm` | Per-fragment uint16 | Already appropriately narrow | Kept uint16 |

The removed logical payload on the VCAP run is roughly 1.56 GiB before
compression and allocator effects: 4 bytes for each of 260.1 M candidate
`intron_bp` entries, plus 16 bytes for each of 32.0 M buffered fragments
across the four stale diagnostic columns. The realized RSS improvement is
smaller because the 2 GiB async spill cap already limits resident chunks.

## Real-Data Sentinel

The first VCAP PR10 profile intentionally tested the proposed uint16
`frag_lengths` narrowing. It aborted during native append with:

```text
Fragment buffer column 'frag_lengths' cannot store value 68466 as uint16
```

That was the correct failure mode. PR10 now documents and enforces the
safer contract: `frag_lengths` remains int32, while `exon_bp` and
`read_length` are the only newly narrowed columns.

## Performance Result

Comparison against the PR07 VCAP profile at the same command line:

```bash
python scripts/profiling/profiler.py \
  --bam /Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/annotated.bam \
  --index /Users/mkiyer/Downloads/rigel_runs/refs/rigel_index \
  --outdir /Users/mkiyer/Downloads/rigel_runs/profile_pr10_buffer_dtype_diet \
  --stages --threads 8 --scan-buffer-size 2 --memory-interval 250
```

| Metric | PR07 | PR10 | Delta |
|---|---:|---:|---:|
| Wall time | 89.15 s | 87.58 s | -1.57 s |
| Peak RSS | 12,086.9 MB | 11,626.2 MB | -460.7 MB |
| RSS after scan | 7,044.2 MB | 6,260.6 MB | -783.7 MB |
| RSS after router scan | 11,906.6 MB | 10,686.2 MB | -1,220.4 MB |
| Scan buffer current | 1,999.5 MB | 1,920.1 MB | -79.4 MB |
| Scan buffer peak | 2,113.1 MiB | 2,124.4 MiB | +11.3 MiB |
| Chunks spilled | 23 | 18 | -5 |
| Scoring CSR bytes | 3.94 GiB | 3.94 GiB | unchanged |
| Partition bytes | 3.94 GiB | 3.94 GiB | unchanged |

The small scan-buffer peak increase is expected noise around the 2 GiB
spill cap and different in-memory/spilled chunk mix. The meaningful PR10
signal is lower RSS after scan/router, lower final peak RSS, and fewer
spilled chunks with unchanged quantification payload sizes.

## Validation

Rebuilt native extensions after C++ edits:

```bash
conda activate rigel
pip install --no-build-isolation -e .
```

Focused checks:

```bash
pytest tests/test_buffer.py tests/test_pipeline_routing.py \
  tests/test_resolution.py tests/test_ndarray_util.py -q
# 111 passed
```

High-level checks:

```bash
pytest tests/test_pipeline_smoke.py tests/test_pipeline_routing.py \
  tests/test_golden_output.py tests/test_profiler.py -q
# 37 passed
```

Full suite:

```bash
pytest tests/ -q
# 1068 passed
```

During full-suite validation, an order-sensitive routing test exposed a
separate native scorer bug: the multimapper gDNA log-sum-exp accumulator
used an infinity sentinel inside `_scoring_impl`, which is compiled with
fast-math. PR10 replaces that sentinel with an explicit initialized flag.

## Final Buffer Dtype Contract

- `frag_id`: int64.
- `t_offsets`, `t_indices`, `frag_lengths`, `genomic_footprint`,
  `genomic_start`: int32.
- `exon_bp`, `read_length`, `num_hits`, `nm`: uint16.
- Classification and strand flags: uint8.

The append boundary in `FragmentAccumulator` is the guard point for
narrowed integer fields. Overflow raises a clear runtime error naming the
column and offending value; values are never silently saturated.

## Non-Goals

- No scoring CSR or partition payload changes; PR07 already handled those.
- No direct resolver API cleanup. Public/debug `ResolvedFragment` fields are
  intentionally retained even when the scan buffer no longer stores them.
- No spill compatibility promise for older intermediate chunk files.
