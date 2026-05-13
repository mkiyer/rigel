# PR 05: Parse Splice-Junction Strand Tags Lazily

## Summary

Avoid scanning BAM auxiliary tags for splice-junction strand on unspliced
records. Only read `XS` or `ts` when the CIGAR contains a reference skip (`N`).

## Motivation

Live samples showed reader-thread cost in `read_sj_strand`, `bam_aux_get`,
`skip_aux`, and `memchr`. The current `parse_bam_record` reads `NM`, `NH`, `HI`,
and then calls `read_sj_strand` before CIGAR parsing for every mapped record.

Most records are unspliced. For those records, `sj_strand` is not needed because
there are no CIGAR splice junctions to annotate.

## Current Code

- `read_sj_strand`: [src/rigel/native/bam_scanner.cpp](../../../src/rigel/native/bam_scanner.cpp)
- `parse_bam_record`: [src/rigel/native/bam_scanner.cpp](../../../src/rigel/native/bam_scanner.cpp)
- Annotated BAM path duplicates similar logic in the same file.

## Proposed Scope

Add a lazy path that keeps the output identical for spliced records and avoids
aux-tag lookup for unspliced records.

Two viable designs:

1. Cheap pre-scan: add `cigar_has_ref_skip(const bam1_t*)`, read `XS/ts` only if
   true, then call existing `parse_cigar`.
2. Lazy callback: modify `parse_cigar` to call a strand-provider callback only
   when it encounters the first `CIG_REF_SKIP`.

The callback design avoids a second CIGAR pass, but the pre-scan design is
simpler and still likely much cheaper than `bam_aux_get` on every record.

## Implementation Steps

1. Add a small helper to detect whether a CIGAR contains `CIG_REF_SKIP`.
2. In `parse_bam_record`, set `rec.sj_strand = STRAND_NONE` unless the helper
   returns true.
3. If the CIGAR is spliced, call `read_sj_strand` and pass that into
   `parse_cigar`.
4. Apply the same change to the annotated-BAM code path so `ZB`/junction tags
   remain consistent.
5. Keep `NM`, `NH`, and `HI` parsing unchanged.

## Tests

Add tests or extend existing BAM tag tests for:

- unspliced record with `XS` tag: no behavior change in resolved fragment, and
  no spurious `sj_strand` requirement
- spliced record with `XS` tag: same `sj_strand` as baseline
- spliced record with `ts` tag on reverse reads: same flip behavior as baseline
- auto-detected `XS`/`ts` behavior unchanged

Run:

```bash
conda activate rigel
pip install --no-build-isolation -e .
pytest tests/test_bam_tag_parsing.py tests/test_splice.py tests/test_strand_model.py -v
pytest tests/test_pipeline_smoke.py tests/test_golden_output.py -v
```

## Benchmark Plan

Run scan profiler and inspect reader-thread samples:

```bash
python scripts/profiling/scan_profile.py \
  --bam /Users/mkiyer/Downloads/rigel_runs/vcap_rna20m_gdna20m/annotated.bam \
  --index /Users/mkiyer/Downloads/rigel_runs/refs/rigel_index \
  --outdir /Users/mkiyer/Downloads/rigel_runs/scan_profile_pr05 \
  --name-prefix pr05 \
  --n-scan-threads 4 \
  --n-decomp-threads 2 \
  --chunk-size 1000000 \
  --max-memory-gib 12
```

Then capture a short `sample` and confirm `read_sj_strand`/`bam_aux_get` frames
are reduced in the reader thread.

## Acceptance Criteria

- Golden outputs are unchanged.
- Strand model observations are identical on the VCAP scan.
- Reader-thread sample shows lower `read_sj_strand` cost.
- Wall time does not regress.

## Risks

- The `ts` tag reverse-read flip must remain exactly as implemented today.
- The annotated BAM path must not diverge from the scan path.

## Non-Goals

- Do not change how `NM`, `NH`, or `HI` tags are parsed in this PR.
- Do not change splice blacklist behavior.
