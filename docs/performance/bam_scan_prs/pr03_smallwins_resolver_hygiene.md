# PR 03: Small-Wins Resolver Hygiene

## Summary

Three independent micro-cleanups in the BAM scanner / resolver inner
loop. None of them is large enough to justify its own PR; together they
trim a few percent of scan wall time and remove three sources of
duplicated work.

Each sub-change has its own tests and acceptance criteria. They can land
together or be split into separate commits within the same PR; do not
split into separate PRs.

* **(a)** Store `ambig_strand` on the accumulator instead of recomputing
  it from per-record fields after `finalize_zero_copy`.
* **(b)** Cheap CIGAR pre-scan before paying for `read_sj_strand`
  (lazy SJ strand parsing).
* **(c)** Replace `std::unordered_map<int32_t, int32_t> frag_lengths`
  with a t-aligned `std::vector<int32_t>`.

## Sub-change (a): Store ambig_strand on the accumulator

### Motivation

`ambig_strand` is computed per-record during `append`, then recomputed
after `finalize_zero_copy` from the assembled vectors. The second
computation is redundant and shows up in samples.

### Current Code

* Accumulator: [src/rigel/native/bam_scanner.cpp](../../../src/rigel/native/bam_scanner.cpp) (`FragmentAccumulator`)

### Change

Add `std::vector<uint8_t> ambig_strand_;` to `FragmentAccumulator`. Push
`r.ambig_strand` in `append`. Return it from `finalize_zero_copy`.
Remove the recomputation block downstream.

### Tests

```bash
pytest tests/test_bam_tag_parsing.py tests/test_strand_model.py \
       tests/test_orient_routing.py -v
pytest tests/test_golden_output.py -v
```

### Acceptance

* All tests pass; golden outputs unchanged.
* Recomputation block deleted (verify by `git grep` for the removed
  function name returning zero matches).

## Sub-change (b): Lazy SJ strand parsing

### Motivation

`read_sj_strand` is called for every record but is only meaningful when
the read actually contains a reference skip ('N') in its CIGAR. The
function pays for tag lookup and string parsing on every record; for the
~85% of records without an 'N' op, this work is wasted.

### Current Code

* SJ strand parsing: in [src/rigel/native/bam_scanner.cpp](../../../src/rigel/native/bam_scanner.cpp)

### Change

Add a tiny helper:

```cpp
inline bool cigar_has_ref_skip(const bam1_t* b) noexcept {
    const uint32_t* cigar = bam_get_cigar(b);
    const uint32_t  n = b->core.n_cigar;
    for (uint32_t i = 0; i < n; ++i) {
        if (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP) return true;
    }
    return false;
}
```

Call it before `read_sj_strand`; skip the call when false. The function
is called inline by the scanner; do not introduce a callback indirection.

### Tests

```bash
pytest tests/test_splice.py tests/test_implicit_splice.py \
       tests/test_strand_model.py -v
pytest tests/test_golden_output.py -v
```

### Acceptance

* All tests pass; golden outputs unchanged.
* On the VCAP profile, `read_sj_strand` self-time drops by ≥ 50%.

## Sub-change (c): T-aligned fragment-length vector

### Motivation

`RawResolveResult` carries a `std::unordered_map<int32_t, int32_t>
frag_lengths` keyed by transcript ID. For every fragment, we allocate a
small hash map, populate it with K entries, then iterate. K is small
(typically 1–4); a sorted parallel vector aligned to the existing
`transcripts` field is faster and zero-allocation when reused.

### Current Code

* Result type: [src/rigel/native/resolve_context.h](../../../src/rigel/native/resolve_context.h)

### Change

```cpp
struct RawResolveResult {
    std::vector<int32_t> transcripts;   // already exists, sorted
    std::vector<int32_t> frag_lengths;  // NEW: same length, same order
    // ...
};
```

Look-ups become `frag_lengths[i]` for the i-th transcript instead of
`frag_lengths_map[transcripts[i]]`. Update the producer in
`_resolve_core` to push in lock-step with `transcripts`. Update every
consumer to index by position.

This sub-change is a representation-only refactor; no semantic change.

### Tests

```bash
pytest tests/test_resolution.py tests/test_frag_length_model.py \
       tests/test_transcript_space_fl.py tests/test_pipeline_routing.py -v
pytest tests/test_golden_output.py -v
```

### Acceptance

* All tests pass; golden outputs unchanged.
* No `unordered_map<int32_t, int32_t>` constructions remain in the
  resolver hot path (`git grep` audit).

## Combined Benchmark Plan

```bash
python scripts/profiling/scan_profile.py \
  --bam /Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/annotated.bam \
  --index /Users/mkiyer/Downloads/rigel_runs/refs/rigel_index \
  --outdir /Users/mkiyer/Downloads/rigel_runs/scan_profile_pr03 \
  --name-prefix pr03 \
  --n-scan-threads 8 --n-decomp-threads 2 \
  --chunk-size 1000000 --max-memory-gib 12
```

Combined acceptance: ≥ 3% wall-time reduction on `s8 d2 nospill` after
PR 01 has landed. The win is small, but it is an honest small win across
three locations the resolver visits on every record.

## Risks

Each sub-change is local and individually reviewable. The largest is
sub-change (c) because it touches every consumer of `frag_lengths`;
`grep` audit is mandatory before merge.
