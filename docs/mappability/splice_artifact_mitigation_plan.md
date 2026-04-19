# Splice Artifact Mitigation — Detailed Plan

**Date:** 2026-04-17
**Status:** Planning (no implementation)
**Scope:** Incorporate a genome-wide splice junction artifact blacklist produced
by the sister tool `alignable` into the Rigel index and BAM scanner so that
spurious spliced alignments (typically driven by gDNA contamination) are
demoted to unspliced fragments before they enter scoring and EM.

---

## 1. Problem Statement

Genomic DNA fragments frequently align as spliced reads when a short
alignment anchor on one side of a candidate intron is sufficient for the
aligner (minimap2 / STAR) to invoke an `N` CIGAR op rather than the
correct unspliced alignment.  These "splice artifacts":

1. are indistinguishable by alignment score from real splice events at the
   read level,
2. contaminate the `spliced_annot` / `spliced_unannot` fragment populations,
3. inflate mRNA abundance at the expense of nRNA/gDNA in the mixture model,
4. are especially harmful in high-gDNA libraries, hybrid-capture libraries,
   and low-input protocols.

Rigel already provides a partial mitigation through the
`SPLICED_ANNOT` vs `SPLICED_UNANNOT` classification.  This is a weak signal:
many artifacts match annotated junctions by coincidence (especially for
short introns or repeat-adjacent junctions), and many real novel junctions
are penalised unnecessarily.

`alignable` solves the detection problem directly.  It tiles the genome with
synthetic gDNA-like fragments at several read lengths, aligns them back to
the same genome, and records every splice junction that appears — these are
by construction **false positives** (there is no RNA in the input).  Each
artifact is characterised by the **maximum anchor length** on either side
that was observed across all false-positive alignments.  Any real alignment
with equal or smaller anchors is therefore indistinguishable from the
artifact population and should be treated as suspect.

We want to consume this blacklist once at `rigel index` time, and apply it
per fragment during BAM scan.

---

## 2. Design Philosophy

We want this feature to **simplify** the scanner, not grow it.  The two
design principles:

1. **Local, per-read filtering.**  A splice-junction anchor is a property
   of a single alignment record.  Filtering belongs in CIGAR parsing, not
   in post-assembly fragment resolution.  After filtering, the assembled
   fragment "just looks unspliced", and every downstream stage (merge,
   resolve, score, EM, calibration) needs **zero** changes.

2. **Data-driven, not heuristic.**  We already have an SJ annotation map
   keyed by `(ref, start, end, strand)` in `FragmentResolver`.  The
   blacklist is the structural twin of that map, keyed the same way.  We
   mirror the existing infrastructure rather than inventing a parallel
   path.

A nice side-effect: Task 1 makes the `SPLICED_UNANNOT` category less
important (most unannotated spurious junctions will be blacklisted), so we
can later consider simplifying the three-way `SpliceType` enum to a binary
spliced/unspliced flag.  We flag this as a candidate cleanup (§10) but do
not take it on here.

---

## 3. Blacklist Format & Aggregation

### 3.1 Input format (from `alignable`)

```
#chrom  intron_start  intron_end  name  count  strand  read_length  max_anchor_5p  max_anchor_3p
chr1    14061         92854       ...   2      .       100          67             35
chr1    14061         92854       ...   2      .       125          53             73
```

Observations from the human reference file on disk
(`/scratch/.../alignable/grch38_star/splice_blacklist.bed.gz`):

- Coordinates are 0-based half-open (intron span), matching our internal
  convention.
- `strand` is always `.`.  Artifacts are strand-agnostic: the aligner
  placed a spurious `N` op regardless of canonical donor/acceptor.
- `max_anchor_5p` and `max_anchor_3p` are **misnamed** — they denote
  *left* and *right* anchors along the read's 5′→3′ direction, not
  genomic 5′/3′.  We will rename them in Rigel to `max_anchor_left` /
  `max_anchor_right` at ingest time.
- The same junction appears on **multiple rows**, one per `read_length`.
  This is the key ambiguity the plan has to resolve.

### 3.2 Aggregation across read lengths

Given rows `(K_i, L_i, R_i)` (read length, left anchor, right anchor) for
one junction, we must collapse to a single pair `(L*, R*)` used at filter
time.  Three sensible choices:

| Rule | Semantics | Trade-off |
|------|-----------|-----------|
| `max` over all rows | Most aggressive filter | Blacklists longer real anchors |
| `max` over rows with `K_i ≤ sample_read_length` | Read-length aware | Requires knowing sample RL |
| Per-RL lookup (no collapse) | Most precise | Storage + complexity |

**Recommendation (phase 1):** aggregate with `max` across all rows per
junction.  Store a single `(start, end, max_anchor_left, max_anchor_right)`
per junction in the index.

Rationale:
- Artifacts grow more extreme at longer read lengths (the anchor can grow
  proportionally).  Taking `max` over all observed `K` yields the
  **conservative frontier** — every anchor at or below this frontier
  was demonstrably producible from gDNA at *some* tested read length.
- Real RNA-seq libraries typically have fragments whose effective aligned
  length is highly variable (trimming, overlap, soft-clipping) and not
  the nominal read length.  Per-RL refinement is fragile.
- We revisit this if benchmarks show unacceptable false-positive filtering.

**Concern logged (§11):** This rule can be overly strict for long-read
junctions seen only once at `K=125` with large anchors.  If a substantial
fraction of true spliced fragments are filtered, we reconsider.

### 3.3 Collapsed internal representation

One row per unique `(ref, start, end)`:

| Column | Dtype | Notes |
|--------|-------|-------|
| `ref` | str (categorical) | Reference name |
| `start` | int32 | 0-based inclusive intron start |
| `end` | int32 | 0-based exclusive intron end |
| `max_anchor_left` | int32 | Aggregated anchor threshold |
| `max_anchor_right` | int32 | Aggregated anchor threshold |

Strand is **intentionally absent**.  Artifacts are strand-agnostic; storing
a sentinel strand column would be dead weight and encourages incorrect
strand-sensitive comparisons.

---

## 4. Anchor Computation from CIGAR

### 4.1 Definition

For a spliced alignment with `k` `N` operations at CIGAR indices
`i_1 < i_2 < … < i_k`, the *j*-th junction's anchors are:

- `anchor_left_j`  = Σ ref-advancing bases (M/=/X/D) over CIGAR ops `[0, i_j)`
- `anchor_right_j` = Σ ref-advancing bases (M/=/X/D) over CIGAR ops `(i_j, end]`

Soft-clip (`S`) and hard-clip (`H`) are excluded.  Insertions (`I`) do not
advance the reference and do not count in anchors.  Deletions (`D`) are
treated as aligned reference positions (the base is placed against the
reference even though absent in the read) — this mirrors how the aligner
established the anchor region.

**Matching the user's examples:**

- `5M 12512N 145M` → `(anchor_left=5, anchor_right=145)`
- `10M 500N 15M 1000N 125M`:
  - junction 1 (500N): left=10, right=15+125=140
  - junction 2 (1000N): left=10+15=25, right=125

This matches the user-specified semantics exactly.

### 4.2 Implementation in `parse_cigar`

`src/rigel/native/bam_scanner.cpp` `parse_cigar()` (around L367) already
walks the CIGAR and emits `exons` and `sjs`.  We extend it to also compute
anchors:

```cpp
// Pre-pass: total ref-advancing bases in the read (the denominator
// for right-anchor arithmetic). This is just the sum of (end-start)
// over the exon blocks we are about to produce, so we compute it
// in a single extra pass — or better, emit a running tally during
// the existing loop and fix up right anchors at the end.
```

The natural pattern:

1. Walk CIGAR left → right, maintaining `left_len` = ref-advancing bases
   seen so far.
2. At each `N` op: record `anchor_left = left_len` for this junction and
   remember the junction's position in the output.
3. After the loop, `total = left_len`.
4. Second pass (or stored during first pass): `anchor_right = total -
   (left_len_at_junction + intron_len)` for each junction.

Cost: O(n_cigar_ops) and adds two `int32_t` fields per SJ entry.  Hot-path
impact is negligible.

### 4.3 Per-read, not per-fragment

Anchors are a **read-level** property.  A paired-end fragment has up to two
reads per alignment, each with its own CIGAR, each contributing its own
anchor observations.  The same genomic junction may be observed in R1 with
large anchors and in R2 with small anchors, or vice versa.

Filtering decision rule (see §5): a junction is blacklist-filtered iff
**every** read instance of that junction fails the anchor test.  If any
mate mate supplies an anchor strictly greater than the blacklist threshold,
the junction is trusted.  This is the permissive rule and matches biological
intuition: real splice events are supported when either mate anchors
deeply.

---

## 5. Filtering Algorithm

### 5.1 Per-read filter (CIGAR parsing)

Extend `ParsedAlignment.sjs` to carry `(start, end, strand, a_left,
a_right)`.  Immediately after CIGAR parsing, for each SJ `i`:

```
hit = blacklist_lookup(ref_id, start, end)
if hit is not None:
    if a_left[i]  <= hit.max_anchor_left \
    && a_right[i] <= hit.max_anchor_right:
        drop sj[i]                # artifactual
```

The dropped SJ simply does not enter `ParsedAlignment.sjs`.  The exon
blocks flanking it are **unchanged** — they remain two separate exon
intervals with a gap.

This is the only mutation.  There is no shadow "removed SJs" list, no flag,
no extra state.

### 5.2 Fragment assembly (`build_fragment`)

Already uses a `std::set<(ref, s, e, strand)>` union across R1/R2.  With
filtering done per-read, the union behaviour gives us the **permissive**
rule automatically: if R1 keeps the SJ and R2 drops it, the union includes
it.  No code changes needed.

### 5.3 Downstream (resolver, EM, calibration)

**Zero changes.**  The fragment presents fewer SJs, its `splice_type` is
re-computed by existing logic:

- If all SJs were blacklisted → `UNSPLICED` (0 introns remain) → naturally
  competes as gDNA / mature / nascent RNA.
- Partial blacklist (some SJs kept) → `SPLICED_ANNOT` or `SPLICED_UNANNOT`
  based on the kept SJs, exactly as today.

The exon blocks still reflect the alignment, but the intron assertion is
gone.  Transcripts that don't contain the blacklisted junction can compete
for the fragment, and gDNA can too (via the unspliced column in the per-
region accumulator).

---

## 6. Index Build Changes

### 6.1 CLI

```
rigel index \
    --fasta genome.fa \
    --gtf genes.gtf \
    --splice-blacklist artifacts.bed.gz \    # NEW, optional
    -o index/
```

When omitted: back-compatible no-op.

### 6.2 New module: `src/rigel/splice_blacklist.py`

Small, focused module (not a wart in `index.py`):

- `load_splice_blacklist(path, ref_to_id) -> pd.DataFrame`
  - Parse BED-ish file (supports gzip via pandas).
  - Validate columns; warn on unknown `ref` values.
  - Aggregate rows per `(ref, start, end)` via `max` over anchors (§3.2).
  - Return sorted DataFrame with schema from §3.3.
- `write_splice_blacklist(df, out_path)` — write `splice_blacklist.feather`.

### 6.3 `index.py` integration

- Add `splice_blacklist_file: Path | None = None` to `TranscriptIndex.build()`.
- If provided, load → aggregate → write `splice_blacklist.feather` alongside
  the existing `sj.feather`.
- Add `SJ_BLACKLIST_FEATHER = "splice_blacklist.feather"` constant.
- `TranscriptIndex.load()` reads the file if present, stashes as a numpy
  arrays-of-columns on the index object (or `None`), and forwards to the
  C++ resolver via a new `build_sj_blacklist_map()` method (symmetric with
  `build_sj_map()`).

### 6.4 C++ `FragmentResolver` additions

Add a second hash map mirroring `sj_map_`:

```cpp
using SJBlacklistKey = std::tuple<int32_t, int32_t, int32_t>;  // ref, s, e
struct SJBlacklistEntry { int32_t max_anchor_left, max_anchor_right; };
std::unordered_map<SJBlacklistKey, SJBlacklistEntry, ...> sj_blacklist_;
```

And a new method:

```cpp
void build_sj_blacklist_map(
    const std::vector<std::string>& refs,
    const std::vector<int32_t>& starts,
    const std::vector<int32_t>& ends,
    const std::vector<int32_t>& max_anchor_left,
    const std::vector<int32_t>& max_anchor_right);
```

Lookup helper on `FragmentResolver`:

```cpp
inline const SJBlacklistEntry* sj_blacklist_lookup(
    int32_t ref_id, int32_t start, int32_t end) const;
```

Returns `nullptr` if junction not in blacklist.  Strand-agnostic, as per
§3.3.

### 6.5 Existing `sj_map_` (annotated junctions) — untouched

The annotated SJ map serves a different purpose (classifying
annotated/unannotated and gathering transcript support).  We keep it.
The two maps answer two questions and their answers are independent:

- `sj_map_` → "do any reference transcripts contain this junction?"
- `sj_blacklist_` → "is this junction a known gDNA-derived artifact at low
   anchor?"

A single junction can match both, one, or neither.

---

## 7. BAM Scanner Changes

### 7.1 `ParsedAlignment` schema

Current:
```cpp
std::vector<std::tuple<int32_t, int32_t, int32_t>> sjs;  // (s, e, strand)
```

Extend to:
```cpp
struct SJCigarEntry {
    int32_t start, end, strand;
    int32_t anchor_left, anchor_right;
};
std::vector<SJCigarEntry> sjs;
```

(Introducing a named struct is a small readability win over the raw
`std::tuple`.)

### 7.2 `parse_cigar()` — compute anchors (§4)

Two-pass implementation is trivial and branch-free.  Keep it in
`bam_scanner.cpp`; no helper function bloat.

### 7.3 New filter hook

Immediately after `parse_cigar` in `parse_bam_record`, before the record
is returned to the caller:

```cpp
filter_blacklisted_sjs(rec.sjs, resolver);
```

Implementation:

```cpp
static inline void filter_blacklisted_sjs(
    std::vector<SJCigarEntry>& sjs,
    const FragmentResolver& r)
{
    if (sjs.empty() || !r.has_sj_blacklist()) return;
    sjs.erase(
        std::remove_if(sjs.begin(), sjs.end(),
            [&](const SJCigarEntry& sj) {
                const auto* hit = r.sj_blacklist_lookup(
                    sj.ref_id, sj.start, sj.end);
                return hit &&
                    sj.anchor_left  <= hit->max_anchor_left &&
                    sj.anchor_right <= hit->max_anchor_right;
            }),
        sjs.end());
}
```

That is the **entire change** to the scanner's hot path.  No
re-plumbing of fragment assembly, no new post-filter classification step.

### 7.4 `build_fragment()`

The existing code destructures `sjs` tuples.  With the struct change, a
two-line diff updates the accessors.  The `intron_set` insertion drops
the anchor fields (they are not needed downstream).  Nothing else changes.

### 7.5 Scan statistics

Add scanner-scope counters (cheap, per-hit):

- `n_sj_observed` — total SJ instances parsed from CIGAR
- `n_sj_blacklisted` — dropped by the filter
- `n_fragments_sj_demoted` — fragments where ≥1 SJ was dropped and
  (post-filter) no SJs remain
- `n_fragments_partially_demoted` — fragments where ≥1 SJ was dropped but
  others survived

Report in `summary.json` so users can see the blacklist's practical
impact.

---

## 8. Sim / Oracle BAM Considerations

`scripts/sim` generates oracle BAMs.  Oracle reads have *correct* spliced
alignments by construction, so the blacklist should almost never fire on
oracle data.  This is a useful invariant to test.

However, simulated gDNA reads aligned via minimap2/STAR may produce
splice artifacts.  After this feature lands, full benchmark runs with the
blacklist active should show:

- negligible change on `nrna_none, gdna_none` conditions,
- meaningful reduction in mRNA-over-estimation on `gdna_high` conditions,
- no change on `oracle`-aligned runs.

These become additional acceptance criteria (§9).

---

## 9. Testing Plan

Unit tests (`tests/test_splice_blacklist.py`, new file):

1. **CIGAR anchor computation**
   - `5M12512N145M` → `[(left=5, right=145)]`
   - `10M500N15M1000N125M` → `[(10,140), (25,125)]`
   - With soft clip: `5S10M100N50M5S` → `[(10, 50)]`
   - With deletion: `10M2D100N50M` → `[(12, 50)]`
   - No splice: `150M` → `[]`

2. **Blacklist loader**
   - Aggregates multi-RL rows via max.
   - Drops duplicate keys.
   - Warns on unknown refs.

3. **Filter semantics**
   - Non-blacklisted junction → always kept.
   - Blacklisted junction, anchors > thresholds → kept.
   - Blacklisted junction, both anchors ≤ thresholds → dropped.
   - Blacklisted junction, one anchor > threshold → kept (logical AND).
   - Strand-agnostic match (blacklist has `.`, read has `+`).

4. **Permissive paired-end rule**
   - R1 has junction with small anchors, R2 has same junction with large
     anchors → union keeps it.
   - Both reads fail → union drops it; fragment becomes UNSPLICED.

Integration tests:

5. **End-to-end scenario BAM**
   - Build a tiny blacklist containing one synthetic artifact junction.
   - Scenario BAM has matching fragments.
   - Verify in the resulting buffer: `splice_type == UNSPLICED`, exon
     blocks preserved, SJ not present.

6. **Golden regression**
   - No blacklist provided → bit-exact parity with existing golden
     outputs (back-compat invariant).

Benchmark tests (`scripts/benchmarking`):

7. Run `gdna_high_ss_1.00_nrna_none` with and without blacklist.
   - Expect mRNA RMSE reduction.
   - Expect gDNA fraction estimate to rise (artifacts re-bin as gDNA).
   - Expect no regression on `gdna_none` conditions.

---

## 10. Code Cleanup Opportunities Enabled by This Feature

Flagged for follow-up, **not** part of this change:

- **SpliceType simplification.**  Once blacklisting is effective, the
  `SPLICED_UNANNOT` category mostly becomes an artifact marker again.  The
  6-column strand × splice table could collapse to 4 columns (spliced /
  unspliced × sense / antisense).  Worth benchmarking post-landing.

- **Annotated BAM output.**  The annotate module already emits SJ metadata;
  consider exposing a new tag `BL:i:1` for fragments whose SJ was dropped.
  Useful for debugging, zero cost.

- **`parse_cigar` ergonomics.**  The function returns data through output
  parameters; refactoring to return a small POD struct would simplify
  callers and is a natural extension now that we're adding fields.

None of these is required.  They are mentioned so future contributors can
see the direction of travel.

---

## 11. Concerns & Open Questions

1. **Aggregation rule (§3.2).**  The `max` rule is conservative but may
   over-filter when only one `read_length` observation exists with large
   anchors.  Alternative: require ≥2 observations across `read_length`
   tiers before trusting a large anchor.  Decide after seeing first
   benchmark results.

2. **Strand-agnostic lookup.**  Blacklist is always strand=`.`; we key on
   `(ref, start, end)` without strand.  This is correct for the current
   `alignable` output but fragile if future versions emit strand-specific
   thresholds.  The loader should assert strand `.` and fail loudly
   otherwise, flagging a schema upgrade.

3. **Interaction with multimappers.**  A multimapping fragment's hits are
   each parsed and filtered independently.  The filtering decisions are
   per-hit, so filtering is correct.  But it is possible for only SOME hits
   to be demoted to `UNSPLICED`.  This is the desired behaviour — each
   alignment location is judged on its own anchors.

4. **Memory footprint.**  For the human reference, the raw blacklist is
   ~8 M rows; after aggregation it collapses to ~2 M unique junctions.  At
   24 bytes per entry plus hashmap overhead, ~100 MB resident in
   `FragmentResolver`.  Acceptable.  Worth monitoring during load.

5. **Anchor semantics under hard-clip.**  Hard-clipped bases are already
   absent from the record; they cannot appear in CIGAR parsing.  Soft-clip
   is excluded by definition.  This matches the `alignable` computation.
   Worth confirming directly with the `alignable` author and documenting.

6. **Should we filter when BOTH anchors fail OR EITHER anchor fails?**  The
   user specified "anchor_left ≤ L_bl *and* anchor_right ≤ R_bl" (both).
   That is the logically correct rule: an artifact is characterised by
   simultaneously small anchors on both sides.  A real junction with a
   small anchor on one side but a long anchor on the other is still well
   supported.  We codify the `AND` rule.

7. **Index format compatibility.**  Adding `splice_blacklist.feather`
   alongside the existing files is additive; old indexes continue to load
   (the loader treats its absence as "no blacklist").  Version bump of the
   index metadata is optional — recommended to stamp a
   `schema_version = N+1` so misloading is detectable.

8. **`alignable` output stability.**  The column names `max_anchor_5p` /
   `max_anchor_3p` are misleading.  We rename on ingest.  If `alignable`
   renames them upstream later, our loader should accept both.

9. **Testing without `alignable`.**  Synthetic blacklists (a single
   hand-crafted junction with controlled anchors) are sufficient for unit
   and integration tests.  Full-genome blacklists belong to the
   benchmark/real-data tier, not the unit suite.

---

## 12. Summary — What Changes, What Doesn't

**Changes (focused, minimal):**

- `src/rigel/splice_blacklist.py` (new, ~50 LOC)
- `src/rigel/index.py`: `--splice-blacklist` arg, feather write/read, C++
  forwarding (~30 LOC)
- `src/rigel/cli.py`: one new CLI flag
- `src/rigel/native/resolve_context.h`: second hash map + lookup helper
  (~40 LOC)
- `src/rigel/native/resolve.cpp`: bind new build method
- `src/rigel/native/bam_scanner.cpp`:
  - `SJCigarEntry` struct + anchor computation in `parse_cigar` (~15 LOC)
  - `filter_blacklisted_sjs` helper called from `parse_bam_record`
    (~20 LOC)
  - Stats counters
- `tests/test_splice_blacklist.py` (new)

**Unchanged:**

- Fragment assembly (`build_fragment`)
- Fragment resolution (SJ classification, chimera, frag-length)
- Scoring, EM, calibration, locus priors, region accumulator
- All existing tests remain green with no blacklist supplied.

The feature is **subtractive at the fragment level** (drop bad SJs) and
**additive at the infrastructure level** (one new map, one new feather
file, one new CLI flag).  The downstream pipeline is untouched — that is
the measure of its elegance.
