# Mappability-Aware Calibration Regions

**Date**: 2026-04-08 (revised 2026-04-10)
**Status**: Planning

## Motivation

Rigel's density pathway estimates gDNA fragment density (λ\_G) from the
10th percentile of unspliced density across genomic regions.  A recent fix
correctly includes zero-count regions in this percentile so that
nRNA-dominated gene bodies do not masquerade as gDNA.

However, on a real human genome, 10–15% of bases are unmappable (centromeres,
telomeres, high-copy repeats, segmental duplications).  Those regions also
accumulate zero counts — because reads *cannot* align there, not because gDNA
is absent.  Including them naively would pull λ\_G → 0 even in the presence
of significant gDNA contamination.

We need the calibration density pathway to distinguish three classes of
zero-count regions:

| Region state | Example | Count | Should calibrate? |
|---|---|---|---|
| Mappable, no gDNA | Intergenic, expressed locus | 0 | **Yes** (true zero baseline) |
| Unmappable | Centromere, Alu cluster | 0 | **No** (uninformative) |
| Mappable, has gDNA | Intergenic with contamination | >0 | **Yes** (gDNA signal) |

The solution is to provide mappability information at index build time so that
calibration can exclude unmappable regions.

---

## The Read-Length Problem

Mappability is read-length-dependent.  A 50 bp read may be ambiguous where a
150 bp read maps uniquely:

| Read length | Approx. unmappable (hg38) |
|:-----------:|:--------------------------|
| 50 bp       | ~15–18%                   |
| 75 bp       | ~12–14%                   |
| 100 bp      | ~10–12%                   |
| 150 bp      | ~8–10%                    |
| 250 bp      | ~5–7%                     |

In a real RNA-seq experiment, fragments have *heterogeneous* effective read
lengths: adapter trimming, quality trimming, paired-end overlap merging, and
soft-clipping all produce a distribution of aligned lengths, not a single
value.  The question "is this region mappable?" has a different answer for
every fragment.

We cannot perform per-fragment mappability lookup for millions of fragments —
calibration needs a single, stable set of trusted regions.

### Conservative Floor Principle

**Key insight**: mappability is monotonically non-decreasing with read length.
A region that is uniquely mappable for 50 bp reads is *guaranteed* uniquely
mappable for all longer reads.  The converse is not true.

Therefore, if we build the blacklist at a short **floor read length** (e.g.,
50 bp), every region that passes the filter is reliably observable for all
fragments in the experiment, regardless of their individual read lengths.

This is the right tradeoff for density calibration:

- **Sufficient**: at 50 bp, ~82–85% of the human genome is mappable.  This
  gives tens of thousands of calibration regions — more than enough for a
  robust percentile estimate.
- **Unbiased**: gDNA fragments shear randomly and their density is spatially
  uniform across the genome.  Measuring density in any representative
  mappable subset gives an unbiased λ\_G estimate.
- **Conservative**: we exclude ~7–10% of the genome that *is* mappable at
  the actual read length.  This loses some statistical power but introduces
  no bias — we are simply discarding regions, not mis-measuring them.
- **Universal**: no per-fragment bookkeeping.  One blacklist, one boolean
  flag per region, one filter in calibration.  Correct for all read lengths
  ≥ floor.

### Does the floor approach distort λ\_G?

No, by two arguments:

1. **The numerator is correct**.  In a floor-mappable region of length L,
   *all* L bases can generate alignable fragments at any read length ≥ floor.
   So `density = count / L` is the true density.

2. **gDNA is spatially uniform**.  Excluding some mappable-at-150bp regions
   does not create a biased sample — those excluded regions would have the
   same gDNA density as the included ones.  We measure the same λ\_G either
   way, just from fewer regions.

The only scenario where the floor could be harmful is if the excluded regions
systematically differ in RNA content (e.g., if regions uniquely mappable at
150 bp contain genes that highly expressed genes while 50 bp-mappable regions
are depleted of genes).  This is unlikely: annotation-driven region boundaries
ensure both gene-dense and gene-sparse regions exist at all mappability
tiers.

### Why per-read-length effective mappable length is unnecessary

An earlier draft of this plan proposed storing per-region effective mappable
length at multiple read lengths (V2).  This is abandoned because:

1. The calibration density pathway only needs a **trustworthy set of
   regions**, not maximum genome coverage.  The floor approach provides this.
2. Fractional effective length adds complexity (which column to select? how
   to handle the distribution of read lengths in a sample?) without
   meaningful benefit — the λ\_G estimate is already unbiased with the floor.
3. Paired-end anchoring means the "effective read length" of a fragment is
   not simply the read length — one uniquely-mapped mate can rescue a mate
   in an ambiguous region.  This makes per-fragment mappability even less
   well-defined.  The floor sidesteps the entire issue.

### Recommended floor: 50 bp

50 bp is the practical minimum for unique short-read alignment.  Reads
shorter than 50 bp after trimming are typically discarded by aligners or
assigned MAPQ 0.  Using 50 bp as the floor ensures that:

- Any fragment that survives alignment QC maps uniquely in floor-mappable
  regions.
- The blacklist errs conservatively: it excludes more than strictly necessary
  for a 150 bp experiment, but never under-excludes.
- A single blacklist BED file per reference genome serves all experiments.

Users who know their minimum post-trimming read length (e.g., 100 bp for a
2×150 experiment) may generate a tighter blacklist at that length for slightly
more statistical power.  But 50 bp is the recommended default.

---

## Architecture Decision

**Option A — Runtime filtering**: Supply a blacklist at `rigel quant` time and
filter in the calibration module.

**Option B — Index-time splitting + flagging (chosen)**: Consume the blacklist
during `rigel index`, split calibration regions at blacklist boundaries, and
add a `mappable` flag to `region_df`.

**Rationale for B**:
- Calibration regions are already defined at index build time.
  Splitting preserves the invariant that each region has uniform annotation
  *and* uniform mappability.
- Index build is a one-time cost per reference genome; no per-sample overhead.
- Mappability is a property of the genome + aligner, not the sample.
- Downstream code (BAM scan, locus priors) needs no changes — the region
  partition is transparent.

---

## Implementation

### CLI

```bash
rigel index \
  --fasta genome.fa \
  --gtf genes.gtf \
  -o index/ \
  --blacklist unmappable.bed      # optional BED3 of unmappable intervals
```

When `--blacklist` is omitted, all regions are `mappable=True` (backward
compatible; identical behavior to current code).

### Index Build (`index.py`)

1. **`load_blacklist(bed_file, ref_lengths)`** — new function.
   - Parse BED3 (ref, start, end; 0-based half-open).
   - Validate coordinates against `ref_lengths`; skip unknown refs with warning.
   - Sort and merge overlapping intervals per reference.
   - Return `dict[str, list[tuple[int,int]]]` keyed by reference name.

2. **`build_region_table(transcripts, ref_lengths, blacklist=None)`** — modified.
   - Inject blacklist start/end coordinates into the boundary set alongside
     transcript/exon boundaries.
   - The existing boundary-sweep algorithm creates atomic bins between
     successive boundaries.  With blacklist boundaries injected, each atomic
     bin is *fully inside* or *fully outside* a blacklisted interval — no
     partial overlap.
   - After forming bins with annotation flags (exon\_pos, exon\_neg, tx\_pos,
     tx\_neg), mark each bin as `mappable=False` if its midpoint falls within
     a blacklisted interval.
   - Adjacent bins with identical `(flags, mappable)` signatures are merged
     as before.
   - New column in output DataFrame: `mappable` (bool, default True).

3. **`TranscriptIndex.build()`** — thread `blacklist_file` parameter through
   to `build_region_table()`.

### Region Schema

```
region_id | ref  | start  | end    | length | exon_pos | exon_neg | tx_pos | tx_neg | mappable
       42 | chr1 |  12100 |  12800 |    700 | True     | False    | True   | False  | True
       43 | chr1 |  12800 |  50200 |  37400 | False    | False    | False  | False  | False    ← centromeric
       44 | chr1 |  50200 |  51000 |    800 | False    | False    | False  | False  | True
```

### Calibration Filtering (`calibration.py`)

4. **`compute_region_stats()`** — pass through `mappable` from `region_df`
   into the returned stats dict.  Default to `True` when the column is
   absent (old index without blacklist).

5. **`calibrate_gdna()` density pathway** — change eligibility:
   ```python
   # Before
   density_eligible = region_length > 0
   # After
   density_eligible = (region_length > 0) & mappable
   ```
   This excludes unmappable zero-count regions while preserving mappable
   zero-count regions as the true zero baseline.

6. **`_build_gdna_fl_model()`** — add `mappable` to the eligible mask.

### No Changes Required

- **`region_cr` / `compute_locus_priors()`** — unmappable regions have
  `n_total=0`, so `E[gDNA]=0` and they contribute nothing to locus γ.
- **BAM scan (C++) / `RegionAccumulator`** — the region partition is
  transparent.  Unmappable regions simply accumulate zero counts.

### Logging

- Index build: `"Built N regions (M mappable, K blacklisted, X.X Mbp excluded)"`
- Calibration: `"N density-eligible regions (K unmappable excluded)"`

---

## Generating the Blacklist

### From `alignable` (recommended)

The sister tool `alignable` computes per-base mappability scores for a given
reference genome and read length.  To generate the blacklist:

```bash
# Generate mappability at floor read length (50 bp recommended)
alignable compute --fasta genome.fa --read-length 50 -o mappability_50bp.bedgraph

# Convert to blacklist BED (regions where unique mappability < threshold)
alignable blacklist --mappability mappability_50bp.bedgraph \
  --threshold 0.5 -o unmappable_50bp.bed
```

The threshold (e.g., 0.5) means: a region is blacklisted if fewer than 50%
of positions within it produce uniquely-mapping reads at 50 bp.

### Interim: ENCODE Blacklist

Until `alignable` is available, the ENCODE blacklist provides a reasonable
approximation for hg38:

```bash
# ENCODE blacklist v2 (hg38)
wget https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz
gunzip hg38-blacklist.v2.bed.gz
rigel index --fasta genome.fa --gtf genes.gtf --blacklist hg38-blacklist.v2.bed -o index/
```

The ENCODE blacklist covers ~100 Mb of problematic regions (centromeres, high
signal artefacts, satellite repeats).  It is less comprehensive than a full
mappability analysis but captures the worst offenders.

---

## Verification Plan

1. **Unit test — region splitting**: Build index with blacklist containing an
   interval that straddles an existing region → verify the region is split
   into 3 parts (pre-blacklist `mappable=True`, blacklist `mappable=False`,
   post-blacklist `mappable=True`).

2. **Unit test — calibration filtering**: Build a scenario with `mappable=False`
   regions carrying zero counts + `mappable=True` regions with real gDNA →
   verify density percentile ignores the unmappable zeros and correctly
   estimates λ\_G.

3. **Integration test — false zero baseline prevention**: Synthetic scenario
   with gDNA + fake unmappable regions → verify λ\_G is not pulled to zero.

4. **Regression**: All existing tests pass with no blacklist (all regions
   `mappable=True`, identical behavior).

5. **Manual**: Build index on hg38 with ENCODE blacklist → inspect
   `regions.tsv` for correct flagging, run `rigel quant` on a real sample,
   verify calibration output is reasonable.

---

## Relevant Source Files

| File | Change | Notes |
|------|--------|-------|
| `src/rigel/cli.py` | `--blacklist` argument | Optional Path |
| `src/rigel/index.py` | `load_blacklist()`, `build_region_table()`, `build()` | Boundary injection + `mappable` column |
| `src/rigel/calibration.py` | `density_eligible` filter | One-line change |
| `tests/test_index.py` | Region splitting tests | New tests |
| `tests/test_calibration.py` | Density filter tests | New tests |

---

## Scope Exclusions

- **Hybrid capture on/off-target**: Requires *sample-specific* region masks
  (the capture kit defines target regions), not genome-level mappability.
  This is a separate feature with a different input mechanism (BED of
  capture targets supplied at `rigel quant` time, not `rigel index` time).
- **De novo mappability computation**: Rigel does not compute mappability
  tracks itself.  That is the job of the sister tool `alignable`.
- **Structural variant / copy number effects**: Mappability in the presence
  of somatic SVs or CNVs is out of scope.
