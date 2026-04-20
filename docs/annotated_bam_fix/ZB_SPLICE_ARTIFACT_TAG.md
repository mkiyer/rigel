# Per-record splice-artifact tag (`ZB:i`) in the annotated BAM

**Status:** Proposed — ready to implement
**Scope:** Pass 2 only (`BamAnnotationWriter`). Pass 1 (`BamScanner`) is **not** modified.
**Contract:** record-level semantic. `ZB:i` = number of splice junctions in this record's CIGAR that were dropped by the splice-artifact blacklist.

---

## Background

Rigel loads an alignable-derived splice-artifact blacklist
([splice_blacklist.py](../../src/rigel/splice_blacklist.py)) into the C++
`FragmentResolver` during `TranscriptIndex.load()`. During BAM scanning, every
CIGAR-derived junction is tested against the blacklist via
`filter_blacklisted_sjs`
([bam_scanner.cpp L460](../../src/rigel/native/bam_scanner.cpp)); junctions
whose anchors fall inside the observed artifact envelope are removed before
fragment assembly. Pass-1 stats (`n_sj_observed`, `n_sj_blacklisted`) are
aggregated but never surfaced per read. Today, users have no way to see which
reads triggered an artifact drop.

## Goal

Stamp each record written by the annotated-BAM writer with
`ZB:i <int>` = the number of its CIGAR junctions that matched the blacklist.
- `ZB:i:0` means the record's splice pattern is clean (or the record is
  unspliced, or no blacklist is loaded).
- `ZB:i:N` with `N > 0` identifies a read whose CIGAR contained `N` blacklisted
  junctions. Those junctions were treated as unspliced during scoring.

The tag is **record-level** (property of a single CIGAR), not fragment-level.
A mate whose CIGAR has 0 blacklisted junctions gets `ZB:i:0` even if its mate
recorded `ZB:i:2`. This matches the semantics of the blacklist filter, which
runs per-read in `parse_record` ([bam_scanner.cpp L585](../../src/rigel/native/bam_scanner.cpp)).

## Contract details

| Aspect | Behaviour |
| --- | --- |
| Tag name | `ZB` |
| Tag type | `i` (signed int; htslib stores smallest fitting width) |
| Presence | Always stamped on records that receive any other `Z*` annotation tag; never stamped on filtered pass-through records (QCFAIL / unmapped / unpaired / duplicates-when-skipped) |
| Value | `int` count of blacklisted CIGAR junctions in this record |
| When no blacklist is loaded | `ZB:i:0` uniformly (we still stamp it, to keep the schema uniform) |
| Relationship to `ZS` | `ZS` reports the fragment's splice class after filtering. `ZB>0` with `ZS=unspliced` is the common case where the only junction(s) were blacklisted |

---

## Where the work lands

Pass 2 only — `BamAnnotationWriter::write` in
[bam_scanner.cpp](../../src/rigel/native/bam_scanner.cpp) around L1937. The
existing writer already:

1. Re-reads the input BAM and parses each CIGAR into a `ParsedAlignment`.
2. Calls `filter_blacklisted_sjs(rec.sjs, *ctx_, mapped_ref_id)` and
   **discards** the return value.
3. Group records by qname → dispatch to `stamp_and_write_hit`, which writes
   uniform per-hit `Z*` tags onto every record in the hit.

The only meaningful changes:

- Capture the return value of `filter_blacklisted_sjs` (it already computes
  the count for free).
- Store the per-record count on `ParsedAlignment`.
- At write time, stamp `ZB:i` per record using that count, outside the
  uniform-tag helper so we do not perturb its signature.

---

## Performance analysis

**Pass 1 — completely unaffected.** Pass 1 does not write the annotated BAM and
does not need to remember per-record blacklist counts. `BamScanner` keeps its
current aggregate-only counters. This is critical — Pass 1 is our slowest
bottleneck.

**Pass 2 — marginal cost analysis:**

| Cost | Per-record increment | Notes |
| --- | --- | --- |
| Capture `int32_t` return of `filter_blacklisted_sjs` | 1 branch + 1 store | The function already computes and returns this value; today we throw it away. **Actual marginal cost: zero.** |
| One extra field on `ParsedAlignment` | ~2 bytes (pack into existing padding) | `uint16_t n_sj_blacklisted` fits into the padding after `flag` (uint16). **Actual marginal cost: zero bytes.** |
| One `bam_aux_update_int(r, "ZB", count)` per record | One aux append (~10 ns) | Negligible next to the 11 `Z*` tags already stamped and the `sam_write1` compression call. |
| Pipeline-level cost when `--annotated-bam` is not set | **None** | Pass 2 runs only when the user opts in. |
| Pipeline-level cost when blacklist is not loaded | The branch `ctx_->has_sj_blacklist()` already short-circuits the filter. `n_sj_blacklisted` is always 0. One extra `bam_aux_update_int` per record writing the constant 0. |

**Bottom line:** no change to Pass-1 hot path. In Pass 2 the blacklist work is
already paid for; we just stop discarding the answer. Expected end-to-end
Pass-2 slowdown: ≤ 1% (driven entirely by the extra tag append, not by any new
parsing).

### Optimization choices we considered and rejected

1. **Stamp ZB only when `> 0`.** Saves ~2-3 bytes per clean record in the
   output BAM. Rejected: makes downstream parsing conditional, adds a branch
   per record. The compressed BAM size delta is effectively lost after BGZF.
2. **Suppress ZB entirely when no blacklist is loaded.** Rejected for the same
   uniformity reason. Downstream tools get a stable schema.
3. **Move the stamping into `stamp_and_write_hit`.** Would force the helper to
   accept a per-record `int` vector and couple hit semantics with a record
   property. Rejected as a clarity regression.
4. **Persist ZB through the annotation table (Pass 1 → Pass 2).** Rejected:
   would add real cost to the hot path (parallel workers serializing per-record
   counts) for data we can recompute trivially in Pass 2 without even a second
   CIGAR parse (the Pass-2 writer is already parsing CIGARs).
5. **Use `ZB:c` (int8)** instead of `ZB:i`. htslib's `bam_aux_update_int`
   already picks the smallest signed type that fits, so on disk a count of 0
   or 1 serializes as a single byte regardless. `ZB:i` is the API type; the
   wire type auto-shrinks. Nothing to gain from a manual `:c`.

---

## Implementation steps

### Step 1 — `ParsedAlignment` gains a per-record counter

File: [src/rigel/native/bam_scanner.cpp](../../src/rigel/native/bam_scanner.cpp)

Add a single `uint16_t` field. Place it adjacent to `flag` so the struct
padding absorbs it — zero effective memory growth.

```cpp
struct ParsedAlignment {
    int32_t ref_id;
    int32_t ref_start;
    int32_t mate_ref_id;
    int32_t mate_ref_start;
    uint16_t flag;
    uint16_t n_sj_blacklisted;  // number of CIGAR junctions dropped by blacklist
    int32_t nm;
    // ...existing fields...
};
```

Initialize to 0 in the default constructor (already value-initialized).

### Step 2 — capture the count in Pass 2's per-record parse

File: [src/rigel/native/bam_scanner.cpp](../../src/rigel/native/bam_scanner.cpp)
(inside `BamAnnotationWriter::write`, around L1941).

```cpp
// Apply the same splice-junction blacklist filter that the
// scanner uses so the annotated BAM sees a consistent view
// of the evidence.  Record the drop count for the ZB tag.
if (ctx_ && ctx_->has_sj_blacklist()) {
    int32_t n_dropped =
        filter_blacklisted_sjs(rec.sjs, *ctx_, mapped_ref_id);
    rec.n_sj_blacklisted = static_cast<uint16_t>(n_dropped);
}
```

Pass 1's `parse_record` already captures the same return value for stats and
is unchanged.

### Step 3 — stamp `ZB` alongside the existing per-record tags

Rather than modify `stamp_and_write_hit` (which is intentionally uniform
across all records in a hit), add a small per-record post-pass inside the hit
loop where we already iterate `hit_raws`. Pseudocode:

```cpp
// Inside the hit loop, after choosing stamp_and_write_hit args,
// stamp ZB per record from its ParsedAlignment.
for (size_t k = 0; k < hit_raws.size(); ++k) {
    // Resolve the ParsedAlignment that produced hit_raws[k].
    bam_aux_update_int(hit_raws[k], "ZB",
        static_cast<int>(parsed[k].n_sj_blacklisted));
}
stamp_and_write_hit(hit_raws, out, hdr, ... /* unchanged */);
```

A cleaner variant (and the one we will use) is to thread an optional pointer
into the helper:

```cpp
static void stamp_and_write_hit(
    const std::vector<bam1_t*>& records,
    const std::vector<uint16_t>* zb_per_record,   // NEW; nullable
    /* ...existing params... */);
```

Inside the helper, after the existing tag writes and before `sam_write1`:

```cpp
for (size_t k = 0; k < records.size(); ++k) {
    int v = zb_per_record ? static_cast<int>((*zb_per_record)[k]) : 0;
    bam_aux_update_int(records[k], "ZB", v);
}
```

The writer builds the parallel `std::vector<uint16_t> zb_vals` at the same
time it builds `hit_raws`. Both vectors are `reserve()`-sized to `r1_reads.size() + r2_reads.size()` (typically 2), so no allocation churn.

Apply identically on the three call sites:
1. Multimap-skip (`include_multimap=False`) full-group stamp.
2. Annotated hit.
3. Intergenic fallback.
4. End-of-group orphan sweep (added during the collation fix).

### Step 4 — docstring / schema table

File: [src/rigel/annotate.py](../../src/rigel/annotate.py)

Add one row to the tag table:

```
   * - ZB
     - i  — number of splice junctions in this record's CIGAR that matched
            the splice-artifact blacklist (record-local, 0 means clean).
```

No Python code changes. The value is computed and stamped entirely in C++.

### Step 5 — tests

File: [tests/test_annotate.py](../../tests/test_annotate.py) — add to
`TestAnnotatedBamIntegration`.

**Test A — schema:** every resolved / intergenic record carries `ZB:i`, and
its value is ≥ 0.

**Test B — correctness on a known blacklisted junction.**

- Use an existing scenario fixture that produces spliced reads over a known
  junction `(ref, start, end)`.
- After `scenario.build(...)`, synthesize a single-row
  `splice_blacklist.feather` in the index directory whose anchors are large
  enough to blacklist the simulated reads.
- Reload the index and run the pipeline with `annotated_bam_path`.
- Assert at least one output record has `ZB:i >= 1`.
- Assert sum of `ZB:i` over all output records equals the Pass-1 stat
  `n_sj_blacklisted` reported in `PipelineResult.stats`. This pins Pass-1 /
  Pass-2 parity.

**Test C — no blacklist loaded:** every `ZB:i` value is `0`.

**Test D — filtered pass-through records carry no ZB.** QCFAIL/unmapped
records must not gain the tag (they carry no other `Z*` tags either).

---

## Compatibility and risk

- **On-disk index format:** unchanged. ZB is derived at Pass-2 time.
- **BAM consumers:** `ZB` is a custom `Z*` tag; unknown tools ignore it. Same
  contract as our existing `ZT`/`ZG`/... tags.
- **Golden regression files:** none include annotated BAM records, so no
  golden regeneration.
- **Interaction with the recent collation fix:** ZB is stamped per record
  after deduplication in the hit loop, so the "one record is written exactly
  once" invariant is preserved. The test
  `test_annotated_bam_preserves_records_with_multimappers` continues to
  guarantee no duplication regression.
- **Deprecation risk:** none. Introducing a new tag is additive.

## Out of scope (deliberately)

- Listing which specific junctions were blacklisted. A follow-up `ZJ:Z` string
  tag (`start1-end1,start2-end2,...`) is design-compatible but not part of
  this change.
- Exposing a Pass-1 per-fragment aggregated blacklist count through the
  quantification outputs (`quant.feather`). Not required by the request and
  would impose real cost on the hot path.

## Rollout

1. C++ change + recompile (`pip install --no-build-isolation -e .`).
2. Add docstring row in `annotate.py`.
3. Add 4 new tests in `TestAnnotatedBamIntegration`.
4. Run full suite. Expect 1056 passed (1052 + 4 new).

## Estimated change footprint

| Location | Approx. lines changed |
| --- | --- |
| `src/rigel/native/bam_scanner.cpp` | ~30 |
| `src/rigel/annotate.py` (docstring) | ~4 |
| `tests/test_annotate.py` (4 new tests) | ~120 |
| **Total** | **~155** |

No new public Python APIs. No breaking changes.
