# PR 04: Use T-Aligned Fragment-Length Vectors

## Summary

Change fragment-length computation from a transcript-id keyed hash map to a
vector aligned with `cr.t_inds`.

## Motivation

`compute_frag_lengths` is hot in native samples. Its current API returns
`std::unordered_map<int32_t, int32_t>`, and `ResolvedFragment::from_core` then
looks up every `t_ind` in that map to build the final parallel array. This pays
hash allocation and hash lookup costs in the scanner hot path.

The result is already consumed as a vector parallel to `t_inds`, so the map is an
unnecessary intermediate representation.

## Current Code

- `compute_frag_lengths`: [src/rigel/native/resolve_context.h](../../../src/rigel/native/resolve_context.h)
- `ResolvedFragment::from_core`: [src/rigel/native/resolve_context.h](../../../src/rigel/native/resolve_context.h)
- `FragmentAccumulator::append`: [src/rigel/native/resolve_context.h](../../../src/rigel/native/resolve_context.h)

## Proposed Scope

Change `RawResolveResult` so it stores:

```cpp
std::vector<int32_t> frag_lengths;
```

where `frag_lengths.size() == t_inds.size()` and each entry corresponds to the
same-position transcript candidate.

## Implementation Steps

1. Add `frag_lengths` to `RawResolveResult` if it does not already exist.
2. Replace `frag_length_map` writes in `_resolve_core` with a call to
   `compute_frag_lengths_aligned(exons, introns, cr.t_inds, cr.frag_lengths,
   scratch)`.
3. Implement the aligned helper:
   - clear and reserve output to `t_inds.size()`
   - for empty input, leave output empty
   - for one exon block, append the same positive length for every `t_ind`
   - for multi-block fragments, append computed transcript-space length or `-1`
4. Change `ResolvedFragment::from_core` to move `cr.frag_lengths` directly.
5. Remove `frag_length_map` from the hot path. Keep a compatibility wrapper only
   if Python-exposed APIs still require dict-like behavior.
6. Rebuild the native extension.

## Tests

Add or update tests for:

- unique fragment length for a single candidate
- ambiguous fragment length across multiple candidates
- empty `t_inds` intergenic fragments
- synthetic nRNA candidates excluded from FL training as before

Run:

```bash
conda activate rigel
pip install --no-build-isolation -e .
pytest tests/test_resolution.py tests/test_frag_length_model.py tests/test_transcript_space_fl.py -v
pytest tests/test_pipeline_smoke.py tests/test_golden_output.py -v
```

## Benchmark Plan

Run the no-spill scan profiler at the current best setting:

```bash
python scripts/profiling/scan_profile.py \
  --bam /Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/annotated.bam \
  --index /Users/mkiyer/Downloads/rigel_runs/refs/rigel_index \
  --outdir /Users/mkiyer/Downloads/rigel_runs/scan_profile_pr04 \
  --name-prefix pr04 \
  --n-scan-threads 4 \
  --n-decomp-threads 2 \
  --chunk-size 1000000 \
  --max-memory-gib 12
```

## Acceptance Criteria

- Golden outputs are unchanged.
- Fragment-length observation counts are identical to baseline.
- Native samples show reduced time in `compute_frag_lengths` and less allocator
  activity around that function.

## Risks

- Some Python compatibility methods expose fragment lengths as a dict. Preserve
  those methods by constructing the dict lazily from `t_inds` and the aligned
  vector only when called from Python.
- Missing lengths must remain represented as `-1` where the old map lookup would
  fail.

## Non-Goals

- Do not narrow fragment-length dtype in this PR. That is a separate memory PR.
- Do not change fragment-length model training policy.
