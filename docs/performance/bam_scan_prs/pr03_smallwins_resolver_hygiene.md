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
* **(b)** One-pass CIGAR parse before paying for `read_sj_strand`
  (lazy SJ strand parsing without a duplicate CIGAR scan).
* **(c)** Replace `std::unordered_map<int32_t, int32_t> frag_lengths`
  with a t-aligned `std::vector<int32_t>`.

## Sub-change (a): Store ambig_strand on the accumulator

### Motivation

`ambig_strand` is computed by `_resolve_core` and copied into
`ResolvedFragment`, then recomputed after `finalize_zero_copy` from the
assembled vectors. The second computation is redundant and shows up in
samples.

### Current Code

* Accumulator: [src/rigel/native/resolve_context.h](../../../src/rigel/native/resolve_context.h) (`FragmentAccumulator`)

### Change

`_resolve_core` already computes `ResolvedFragment::ambig_strand`. Add
`std::vector<uint8_t> ambig_strand_;` to `FragmentAccumulator`, reserve
it, push `r.ambig_strand` in `append`, and return it from both
`finalize` and `finalize_zero_copy`. Keep the `t_strand_arr` parameter
temporarily if removing it would churn bindings, but stop using it for
this column.

### Tests

```bash
conda activate rigel
pip install --no-build-isolation -e .
pytest tests/test_bam_tag_parsing.py tests/test_strand_model.py \
       tests/test_orient_routing.py -v
pytest tests/test_golden_output.py -v
```

### Acceptance

* All tests pass; golden outputs unchanged.
* Recomputation block deleted from `finalize` and `finalize_zero_copy`.
  The only remaining `ambig_strand` computation should be in the
  resolver itself.

## Sub-change (b): Lazy SJ strand parsing

### Motivation

`read_sj_strand` is called for every record but is only meaningful when
the read actually contains a reference skip ('N') in its CIGAR. The
function pays for tag lookup and string parsing on every record; for the
~85% of records without an 'N' op, this work is wasted.

### Current Code

* SJ strand parsing: in [src/rigel/native/bam_scanner.cpp](../../../src/rigel/native/bam_scanner.cpp)

### Change

Do not add a second CIGAR pass. Parse the CIGAR once with
`STRAND_NONE`, then read `XS`/`ts` only if `rec.sjs` is non-empty and
patch the already-built junction entries:

```cpp
parse_cigar(b, mapped_ref_id, STRAND_NONE, rec.exons, rec.sjs);
rec.sj_strand = STRAND_NONE;
if (!rec.sjs.empty()) {
  rec.sj_strand = read_sj_strand(b, sj_tag_mode);
  for (auto& sj : rec.sjs) sj.strand = rec.sj_strand;
}
```

Apply the same pattern in both `parse_bam_record` and the annotated-BAM
writer path. This avoids both the aux-tag scan on unspliced records and
the duplicate CIGAR pre-scan.

### Tests

```bash
conda activate rigel
pip install --no-build-isolation -e .
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
frag_length_map` keyed by transcript ID. For every fragment, we allocate
a small hash map, populate it with K entries, then iterate. K is small
(typically 1–4); a parallel vector aligned to the existing `t_inds`
field is faster and zero-allocation when reused.

### Current Code

* Result type: [src/rigel/native/resolve_context.h](../../../src/rigel/native/resolve_context.h)

### Change

```cpp
struct RawResolveResult {
  std::vector<int32_t> t_inds;        // already exists, sorted
  std::vector<int32_t> frag_lengths;  // NEW: same length, same order
  // ...
};
```

Look-ups become `frag_lengths[i]` for the i-th `t_inds` entry instead
of `frag_length_map[t_inds[i]]`. Update `compute_frag_lengths` into a
`compute_frag_lengths_aligned(..., cr.t_inds, cr.frag_lengths, scratch)`
producer, update `ResolvedFragment::from_core` to move the vector
directly, and update the legacy `resolve()` tuple builder to reconstruct
the Python dict from `t_inds[i]` / `frag_lengths[i]`.

This sub-change is a representation-only refactor; no semantic change.

### Tests

```bash
conda activate rigel
pip install --no-build-isolation -e .
pytest tests/test_resolution.py tests/test_frag_length_model.py \
       tests/test_transcript_space_fl.py tests/test_pipeline_routing.py -v
pytest tests/test_golden_output.py -v
```

### Acceptance

* All tests pass; golden outputs unchanged.
* No `frag_length_map` field or `unordered_map<int32_t, int32_t>`
  construction remains in the fragment-length hot path (`git grep`
  audit). Other resolver maps used for index lookups may remain.

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
