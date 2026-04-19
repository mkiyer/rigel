# Mappability-Aware Calibration Regions — V2

**Date**: 2026-04-18
**Supersedes**: `mappability_aware_regions_plan.md` (V1)
**Status**: Planning — approved to implement

---

## 1. Summary of changes from V1

V1 proposed taking an **unmappable BED** and splitting calibration regions at
its boundaries, flagging split pieces with a `mappable` boolean.  V2 keeps
V1's core philosophy (index-time integration, one-time cost, transparent to
downstream code) but refines four things in light of the actual input we now
have:

1. **Input is a *mappable* BED, not an unmappable BED.**  The `alignable`
   tool emits "keep" regions at a given floor read length.  Rigel consumes
   them directly — no inversion in user space, and the contract is aligned
   with the tool's natural output.

2. **Rigel prunes the BED itself.**  `alignable` emits *every* mappable
   interval, including slivers a few bp wide that cannot host a full
   fragment.  Leaving this to the user is fragile; Rigel knows the
   fragment-length regime it needs to support and must guarantee its own
   calibration quality.  A user-overridable minimum region length
   (default **500 bp**) is applied at index build time.

3. **Mappable regions are carried as a first-class index artifact**, not
   merely as a boolean column on `region_df`.  They live in their own
   feather file (`mappable.feather`) and get their own cgranges interval
   tree when loaded.  `region_df` is not altered structurally (no boundary
   injection, no schema change) — calibration instead joins existing
   regions against the mappable interval tree at calibration time.

4. **Density-eligibility tests fragment containment, not midpoint
   membership.**  A region qualifies for density-based λ\_G estimation only
   if it is *entirely inside* one mappable interval AND its length is
   ≥ `min_mappable_region_length`.  This is a conservative, defensible
   rule for both numerator (counts) and denominator (length) of the
   density estimate, and it doubles as the filter for the gDNA
   fragment-length model.

The net effect: one new feather file, one new cgranges index, one
conservative filter that replaces `density_eligible = length > 0`, and no
churn in the scanner, locus priors, or EM code paths.

---

## 2. Why not boundary-split `region_df` (V1's approach)?

The V1 plan proposed injecting mappable/unmappable boundaries into the
region-boundary set so that every atomic bin would be entirely mappable or
entirely unmappable.  On closer inspection this is the wrong layer:

- **Region semantics**: existing region boundaries encode *annotation
  structure* (exon/intron/transcript tiling).  Mappability is a genome
  property orthogonal to annotation.  Mixing the two produces region tables
  that are harder to reason about and bloats `region_df` with thousands of
  slivers that exist only to carry a boolean.

- **Granularity**: mappable intervals from `alignable` can be numerous
  (hundreds of thousands to millions on a vertebrate genome).  Folding them
  into the boundary set could multiply `|region_df|` by an order of
  magnitude, most of which is noise from our standpoint (tiny fragments of
  intergenic space).  The interval-tree / containment-test approach avoids
  this explosion entirely.

- **Correctness of density**: what density calibration actually wants is
  *regions it can fully trust*.  "Trust" means the entire region can host
  a gDNA fragment and be observed.  This is a **whole-region predicate**,
  not a fraction.  Boundary splitting optimises the wrong thing (coverage)
  at the cost of the right thing (simplicity and defensibility).

Instead, V2 preserves the invariants:

> `region_df` describes annotation structure.
> Mappability is a separate, sample-independent property of the genome.
> Calibration composes the two at query time.

This also leaves a clean path for future per-sample masks (capture kits,
custom exclusions) to compose the same way without touching the index.

---

## 3. Input: mappable BED from `alignable`

`alignable` emits a BED3 file of uniquely-mappable intervals at a chosen
floor read length (we use 50 bp as discussed in V1 §"Read-Length Problem").
Example (first lines, gzipped):

```
chr1    10468    10518
chr1    10601    10695
chr1    10700    10749
...
chrY    57217414 57217459
```

Characteristics observed:

- BED3 (3 columns, 0-based half-open), possibly gzipped.
- Intervals are non-overlapping and sorted per reference.
- Sizes range from tens of bp up to tens of kb.
- Large fraction of intervals are < 500 bp: these cannot reliably host an
  RNA-seq fragment whose length distribution is ~150–500 bp, and give
  unstable `count / length` density estimates.

Rigel's contract with `alignable`:

- Accept BED3 or BED3+extra-columns (extras ignored).
- Accept plain or gzipped.
- Be tolerant of reference names absent from the genome FASTA: emit a
  debug-level note, skip those lines, continue.
- Never fail the index build because the mappable BED is malformed on a
  single line — mappability is a soft dependency.

---

## 4. Minimum region length

### The problem

`alignable` emits every uniquely-mappable interval, including ones of a
handful of bases.  Two failure modes result from including them unfiltered
in density calibration:

1. **Fragment-containment failure.**  An RNA-seq fragment is ~150–500 bp;
   gDNA fragments after shearing span a similar range.  A 60-bp mappable
   island cannot host an entire fragment — any fragment overlapping it is
   also overlapping the neighbouring unmappable sequence, so its presence
   or absence tells us nothing clean about density in the mappable
   interval.

2. **Denominator noise.**  `density = count / length` is dominated by
   Poisson noise when `length ≪ E[fragment_length]`.  Including tiny
   intervals injects high-variance outliers into the low-density tail
   used for the percentile, biasing λ\_G downward.

### The rule

Discard mappable intervals shorter than `min_mappable_region_length`.

**Default: 500 bp.**  Rationale:

- Larger than the 99th-percentile fragment length for typical short-read
  RNA-seq (~400 bp) → a full fragment fits with margin.
- Coarse enough that Poisson noise in `count` is dominated by signal when
  the region carries real density.
- Empirically preserves a large fraction of the mappable genome (> 80%
  of floor-mappable bases on hg38 survive at a 500-bp minimum, because
  the mappable landscape is dominated by large contiguous intervals).
- Gives tens to hundreds of thousands of calibration intervals on a
  vertebrate genome — more than enough for a stable percentile.

**Exposed as an advanced CLI flag**, `--min-mappable-length` (dest
`min_mappable_region_length`), defaulting to 500.  We document it as an
advanced parameter: most users should never change it.  Exposing it lets us
(a) tune in evaluation, (b) cover shorter-fragment libraries (e.g. degraded
RNA, microRNA) without rebuilding `alignable` output.

Rigel performs the filter because Rigel is responsible for Rigel's
accuracy.  Pushing this to the producer would make calibration fragile
against inputs generated by other pipelines or older `alignable` versions.

---

## 5. Data model

### 5.1 New index artifact

```
index_dir/
    ...existing files...
    mappable.feather          # new, optional
    mappable.tsv              # optional mirror
```

Schema of `mappable.feather`:

| Column | Dtype | Notes |
|---|---|---|
| `ref`   | category (string categories) | Reference name, matches `ref_lengths.feather` |
| `start` | int32 | 0-based |
| `end`   | int32 | half-open |
| `length`| int32 | `end − start`, always ≥ `min_mappable_region_length` |

Sorted by `(ref, start)`.  Non-overlapping by construction (input is
already non-overlapping; length filtering preserves this).

### 5.2 No change to `region_df`

Region schema stays:
```
region_id | ref | start | end | length | exon_pos | exon_neg | tx_pos | tx_neg
```

Rationale: the mappability dimension is compositional, not intrinsic to the
annotation partition.  Forcing it into `region_df` would conflate two
concepts.

### 5.3 Runtime representation

On `TranscriptIndex.load()`:

- Read `mappable.feather` if present.
- Build a cgranges index `mappable_cr_` (one interval tree, all references
  pooled with ref_id encoding — same pattern as existing `region_cr_` and
  `sj_region_cr_`).
- Expose `idx.has_mappability()` and `idx.mappable_df` on the Python side.
- Absence of the file is fine: `has_mappability()` returns False; all
  calibration code paths degrade gracefully to current behaviour.

---

## 6. Calibration changes

The change is surgical: one predicate replaces `density_eligible = length > 0`.

### 6.1 `density_eligible` becomes "fully contained in a mappable region"

```python
# calibration.py
def _compute_mappable_mask(
    region_df: pd.DataFrame,
    mappable_cr,                      # cgranges or None
) -> np.ndarray:
    """True iff [start, end) is entirely inside a single mappable interval."""
    if mappable_cr is None:
        # No mappability info → fall back to current behaviour.
        return np.ones(len(region_df), dtype=bool)
    ...  # interval-tree lookup per region (vectorisable via batched query)
```

Then in `calibrate_gdna`:

```python
density_eligible = (region_length > 0) & mappable_mask
```

Nothing else in `calibrate_gdna` changes.  The percentile, the blend, the
downstream `e_gdna_density` all inherit the filter.

### 6.2 `_build_gdna_fl_model` gets the same treatment

The density-pathway FL selection uses `d_unspliced[eligible]` today; we
additionally AND in `mappable_mask` so the gDNA fragment-length model is
trained only from fragments in trusted regions.  This improves the FL model
on unstranded libraries as a bonus.

### 6.3 Containment test, not overlap

The eligibility predicate is **containment**, not overlap.  A region that
starts inside a mappable interval and extends past its boundary is
*rejected*.  This is the conservative choice:

- A region that straddles a mappability boundary has mixed reliability —
  part of its count comes from observable genome, part does not.  The
  `count / length` ratio cannot be cleanly interpreted.
- The number of regions we lose this way is small: most `region_df`
  regions are either entirely intergenic (large, likely to fit inside a
  large mappable block) or entirely exonic (within expressed transcripts
  which are almost always in mappable regions).

Containment is implemented via cgranges: for each region `(ref, s, e)`,
query the mappable tree; if exactly one mappable interval `[m_s, m_e)`
overlaps AND `m_s ≤ s AND e ≤ m_e`, the region is eligible.

### 6.4 Interaction with density percentile

The percentile is unchanged (still 10th by default).  We just compute it
over a smaller, higher-quality population.  This directly addresses the
motivating failure mode: previously the 10th percentile could be pulled to
zero by unmappable regions with zero counts; now it reflects the true
floor of gDNA density in trustworthy genome.

### 6.5 Logging

```
[CAL] Density pathway: 142,318 / 1,847,291 regions eligible
      (mappable & length > 0, 92% exclusion from mappability filter)
[CAL] Density pathway: λ_G=4.23e-05 (percentile=10)
```

### 6.6 `compute_region_stats` gains a `mappable` field

For reporting and for downstream diagnostics, we surface the boolean mask
in `stats["mappable"]` (all True when no blacklist was supplied).  This is
purely additive — no caller is broken.

---

## 7. Index build flow

```
rigel index \
    --fasta genome.fa \
    --gtf genes.gtf \
    -o index/ \
    --mappable mappable_50bp.bed.gz   \  # optional
    --min-mappable-length 500            # advanced, default 500
```

Build order inside `TranscriptIndex.build()`:

1. ref_lengths, transcripts, intervals, sj, regions (unchanged).
2. **Mappable regions (new, last):**
   - If `--mappable` not given → skip, no file written.
   - Else:
     - Parse BED3 (gzip-tolerant) via a new helper
       `rigel.mappable.load_mappable_bed()`.
     - Validate refs against `ref_lengths`; drop unknown refs with a count
       in the log.
     - Clip intervals to `[0, ref_length)` (defensive — spec-compliant
       producers never exceed this, but let's not crash on a bad file).
     - Sort by `(ref, start)`; merge adjacent/overlapping intervals (also
       defensive — `alignable` output is already clean, but the helper
       should be robust).
     - Drop intervals with `length < min_mappable_region_length`.
     - Persist to `mappable.feather` (+ optional TSV).
     - Log: `Mappable: N_in input → N_clean clean → N_kept kept
       (≥{L}bp), covering X.X Mb`.

`build_region_table` is not touched.  The mappable BED's only effect on
disk is the new feather file.

---

## 8. New Python module: `rigel.mappable`

Mirrors the style of `rigel.splice_blacklist` introduced for task 1 — a
small, single-purpose module that owns ingestion + validation + canonical
schema.

```
src/rigel/mappable.py
    MAPPABLE_COLUMNS = ("ref", "start", "end", "length")
    def load_mappable_bed(path, *, ref_lengths=None,
                          min_length: int = 500) -> pd.DataFrame: ...
```

The function returns the filtered DataFrame ready for feather persistence.
It accepts `ref_lengths` as an optional parameter so validation/clipping
can happen inline.  The min-length filter happens here, not in
`TranscriptIndex.build`, because the module is the single place where the
contract with `alignable` lives.

`TranscriptIndex.load()` gets a mirror method `load_mappable_index()` that
reads the feather and builds the cgranges tree; this is called lazily and
memoised on the instance.

---

## 9. Interval-tree representation

We already have two cgranges indices in the index (`region_cr_`,
`sj_region_cr_`); adding a third follows the established pattern and costs
very little.

Sketch in `TranscriptIndex.load()`:

```python
if (index_dir / MAPPABLE_FEATHER).exists():
    m_df = pd.read_feather(index_dir / MAPPABLE_FEATHER)
    ix.mappable_df = m_df
    cr = _cgranges_cls()
    for row in m_df.itertuples(index=False):
        cr.add(row.ref, row.start, row.end, 0)
    cr.index()
    ix.mappable_cr_ = cr
else:
    ix.mappable_df = None
    ix.mappable_cr_ = None
```

Calibration then uses `mappable_cr_.overlap(ref, s, e)` or a vectorised
batched query; the containment test wraps this.  **No C++ changes.**

---

## 10. Code-health improvements delivered in the same PR

This feature touches naturally adjacent code; we use the opportunity to
clean up, consistent with the user's directive.

1. **Small-BED-loader refactor.**  `splice_blacklist.py` (task 1) and
   `mappable.py` (task 2) share logic: gzip-tolerant parsing, reference
   validation, sort-merge.  Extract a private helper
   `_read_bedlike(path, expected_cols, schema)` into a new
   `rigel._bed_io` module.  Both loaders become thin wrappers.  Keeps each
   public module focused on its own schema and aggregation rules while
   removing duplicated I/O.

2. **`TranscriptIndex` feather manifest table.**  Today the class has five
   `*_FEATHER` constants, a sixth for splice blacklist, and now a seventh
   for mappability, with their load/save code sprinkled across
   `build()` / `load()`.  Introduce a small internal registry:

   ```python
   _INDEX_ARTIFACTS: tuple[_ArtifactSpec, ...] = (
       _ArtifactSpec("ref_lengths",   REF_LENGTHS_FEATHER,   required=True),
       _ArtifactSpec("transcripts",   TRANSCRIPTS_FEATHER,   required=True),
       _ArtifactSpec("intervals",     INTERVALS_FEATHER,     required=True),
       _ArtifactSpec("sj",            SJ_FEATHER,            required=True),
       _ArtifactSpec("regions",       REGIONS_FEATHER,       required=True),
       _ArtifactSpec("splice_blacklist", SJ_BLACKLIST_FEATHER, required=False),
       _ArtifactSpec("mappable",      MAPPABLE_FEATHER,      required=False),
   )
   ```
   Consolidates the "feather + optional TSV + feather_kwargs" boilerplate
   and gives us a single place to add the next optional artifact.

3. **`calibration.py` diet.**  Extract the ~30 lines of
   "mask construction" (`eligible`, `has_strand`, `density_eligible`)
   into a named helper `_build_eligibility_masks(stats, mappable_mask)`
   returning a small dataclass.  `calibrate_gdna` becomes more linear.

4. **Docs pass.**  Update `docs/MANUAL.md` and `docs/parameters.md` to
   describe `--mappable` and `--min-mappable-length`; cross-link from the
   calibration section.

Each of these is small, purely mechanical, and makes the resulting code
shorter and easier to extend.  None is load-bearing for the feature, and
each can be separated into its own commit within the PR for bisectability.

---

## 11. Test plan

New file `tests/test_mappable.py`:

### Loader tests (mirror `test_splice_blacklist.py`)
- BED3 → feather round-trip preserves contents.
- Gzipped input works.
- Unknown references → dropped with a log note, not raised.
- Overlapping/adjacent intervals → merged.
- Short intervals → dropped at `min_length` boundary.
- Empty BED → empty DataFrame with correct schema.
- Dtype narrowing (int32) preserved.

### Index integration
- `TranscriptIndex.build(..., mappable_file=bed, min_mappable_length=500)`
  writes `mappable.feather` with filtered rows.
- `TranscriptIndex.load()` builds `mappable_cr_` cgranges tree.
- `idx.has_mappability()` toggles correctly.
- Backwards compatibility: indexes built without `--mappable` load fine
  and expose `mappable_df is None`.

### Calibration semantics
- Containment test: region fully inside a mappable interval → eligible.
- Region straddling a mappability boundary → **not** eligible.
- Region in a gap → not eligible.
- Synthetic scenario: unmappable zero-count regions no longer pull the
  density percentile to zero.  Build a tiny index with a fake "centromere"
  block, scan a synthetic BAM with gDNA fragments in mappable parts only,
  and assert λ\_G\_density matches ground truth within tolerance.
- Regression: no mappable file → every test in
  `test_calibration_*` still passes unchanged.

### gDNA FL model
- With mappability filter, `_build_gdna_fl_model` selects fragments only
  from contained mappable regions on the density pathway.

### End-to-end
- `rigel index --mappable ...` followed by `rigel quant` completes on the
  scenario corpus with no regressions in golden outputs.
- For a scenario designed to have unmappable zero regions, the λ\_G
  reported in `summary.json` is measurably higher than without the filter.

### Properties / edge cases
- `--min-mappable-length 0` yields the pass-through behaviour.
- Very large `--min-mappable-length` (e.g. 10_000) → few or zero eligible
  regions → density pathway disables itself and logs appropriately; the
  strand pathway continues to work.

---

## 12. Rollout

1. Implement `rigel.mappable` + BED-loader tests.
2. Wire `--mappable` / `--min-mappable-length` through `cli.py`, `index.py`.
3. Add `mappable.feather` to the index build, lazy `mappable_cr_` on load.
4. Add the containment helper + `density_eligible` change in
   `calibration.py`.
5. Add `mappable` mask in `_build_gdna_fl_model` selection.
6. Land the three code-health refactors (BED-loader extraction, artifact
   registry, eligibility helper) as follow-up commits in the same PR.
7. Regenerate golden outputs only where the mappability filter changes
   λ\_G for the affected scenarios; document the reasoning in the commit.
8. Update `docs/MANUAL.md`, `docs/parameters.md`.

---

## 13. Scope exclusions (unchanged from V1)

- Sample-specific masks (capture kits, custom exclusions).
- De novo mappability computation.
- Structural-variant / copy-number effects.
- Read-length-specific mappability lookups per fragment (the floor
  principle from V1 §"Read-Length Problem" remains the answer).

---

## 14. Design summary

> One optional BED in → one optional feather out → one optional cgranges
> tree loaded → one-line change in the density-eligibility predicate.
> No scanner changes, no region-partition churn, no schema migrations.
> The feature enters and exits cleanly, leaving `region_df` alone and
> pushing the mappability concept into its own artifact where it belongs.
