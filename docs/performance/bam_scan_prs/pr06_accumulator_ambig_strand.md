# PR 06: Store `ambig_strand` In The Fragment Accumulator

## Summary

Stop recomputing `ambig_strand` during chunk finalization. `_resolve_core`
already computes the value for each resolved fragment; the accumulator should
store it directly.

## Motivation

`FragmentResolver::_resolve_core` computes `cr.ambig_strand` after candidate
resolution. `ResolvedFragment::from_core` copies that value. But
`FragmentAccumulator::finalize_zero_copy` scans the CSR candidate arrays again to
recompute `ambig_strand` for every fragment in the chunk.

This is redundant work on every emitted chunk and requires passing
`t_strand_arr` into finalization.

## Current Code

- `_resolve_core` ambig calculation:
  [src/rigel/native/resolve_context.h](../../../src/rigel/native/resolve_context.h)
- `ResolvedFragment::from_core`:
  [src/rigel/native/resolve_context.h](../../../src/rigel/native/resolve_context.h)
- `FragmentAccumulator::finalize_zero_copy`:
  [src/rigel/native/resolve_context.h](../../../src/rigel/native/resolve_context.h)

## Proposed Scope

Add `std::vector<uint8_t> ambig_strand_` to `FragmentAccumulator` and push
`r.ambig_strand` in `append`. Finalization should move that vector directly into
the output dictionary.

## Implementation Steps

1. Add `ambig_strand_` to `FragmentAccumulator`.
2. Reserve it in `reserve(n_fragments, n_candidates)`.
3. Push `static_cast<uint8_t>(r.ambig_strand)` in `append`.
4. Update `finalize_zero_copy` to return `ambig_strand_` directly.
5. Update `finalize` similarly if it remains in use for Python-side appends.
6. Keep the `t_strand_arr` parameter temporarily if needed for ABI stability,
   but stop using it in the hot path.
7. Rebuild native extensions.

## Tests

Add a focused test for a fragment with mixed-strand candidates to ensure the
stored value matches the previous recomputed value.

Run:

```bash
conda activate rigel
pip install --no-build-isolation -e .
pytest tests/test_buffer.py tests/test_resolution.py tests/test_opp_strand_order.py -v
pytest tests/test_golden_output.py -v
```

## Benchmark Plan

This is a small optimization, so the main benchmark is no regression:

```bash
python scripts/profiling/scan_profile.py \
  --bam /Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/annotated.bam \
  --index /Users/mkiyer/Downloads/rigel_runs/refs/rigel_index \
  --outdir /Users/mkiyer/Downloads/rigel_runs/scan_profile_pr06 \
  --name-prefix pr06 \
  --n-scan-threads 4 \
  --n-decomp-threads 2 \
  --chunk-size 1000000 \
  --max-memory-gib 12
```

## Acceptance Criteria

- Golden outputs are unchanged.
- `ambig_strand` arrays are byte-for-byte identical to baseline on test chunks.
- `finalize_zero_copy` no longer scans `t_indices` to compute ambig strand.

## Risks

- Intergenic fragments with zero candidates must still get `ambig_strand = 0`.
- Any Python-only accumulator path must append the correct stored value.

## Non-Goals

- Do not change how `ambig_strand` is defined.
- Do not narrow or remove other buffer columns in this PR.
