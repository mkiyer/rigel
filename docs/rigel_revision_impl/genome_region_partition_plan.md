# Region Partition (Flattened Index) — Implementation Plan

## Goal

Add a non-overlapping, reference-complete genomic partition to the index. Each
atomic region carries four boolean flags—`exon_pos`, `exon_neg`, `tx_pos`,
`tx_neg`—from which all annotation context (intergenic, intronic, exonic,
strand-ambiguous, etc.) is derived at query time.

## Schema: `regions.feather`

| Column | Type | Description |
|---|---|---|
| `region_id` | int32 | Sequential, globally unique |
| `ref` | string | Reference name |
| `start` | int32 | 0-based start (inclusive) |
| `end` | int32 | 0-based end (exclusive) |
| `length` | int32 | `end - start` |
| `exon_pos` | bool | Region overlaps a + strand exon |
| `exon_neg` | bool | Region overlaps a – strand exon |
| `tx_pos` | bool | Region overlaps a + strand transcript span |
| `tx_neg` | bool | Region overlaps a – strand transcript span |

Derived interpretations (not stored):
- **intergenic**: all four flags False
- **intronic (+)**: `tx_pos` and not `exon_pos`
- **exonic (+)**: `exon_pos` (implies `tx_pos`)
- **strand-ambiguous**: any pos flag AND any neg flag
- **genic**: any flag True

## Algorithm: boundary sweep

The core builder is a pure function:

```python
def build_region_table(
    transcripts: list[Transcript],
    ref_lengths: dict[str, int],
) -> pd.DataFrame:
```

**Per reference:**

1. **Collect boundaries.** For each transcript on the reference, add `t.start`
   and `t.end` to a boundary set. For each exon of each transcript, add
   `e.start` and `e.end`. Also add `0` and `ref_length` as sentinel boundaries.

2. **Sort and deduplicate** into a sorted array of unique boundary coordinates.

3. **Form atomic bins.** Each consecutive pair `(boundaries[i], boundaries[i+1])`
   defines one atomic half-open region. Discard zero-length bins (duplicates
   already removed, but guard anyway).

4. **Assign flags.** For each bin, determine which of the four flags are True by
   testing whether any transcript span or exon overlaps the bin. Since the bin
   boundaries are drawn exactly from transcript/exon coordinates, a simple
   sorted-endpoint structure suffices:
   - Build sorted lists of `(start, end, strand)` for transcript spans and exon
     spans.
   - For each bin midpoint (or start), use `np.searchsorted` on the boundary
     array to find the bin range, then set flags.

5. **Assign `region_id`** sequentially across all references (first ref gets
   0..k, second gets k+1..m, etc.)

**Concrete `searchsorted` approach for flag assignment:**

Each transcript span `[start, end)` maps to a contiguous range of bins; each
exon likewise. Use `np.searchsorted` on the boundary array to find the first and
last bin index covered, then set the appropriate strand flag for all bins in that
range. This is O(n_transcripts × log(n_bins)) and straightforward.

## Code placement

All builder code goes in `src/rigel/index.py` as a new top-level function
`build_region_table()`, alongside the existing `build_genomic_intervals()`. No
new module needed for the first pass.

## Integration into `TranscriptIndex.build()`

Add a section after the genomic intervals block:

```python
# -- Calibration regions ------------------------------------------
logger.info("[START] Building calibration regions")
region_df = build_region_table(transcripts, ref_lengths)
logger.info(f"[DONE] Found {len(region_df)} calibration regions")

region_df.to_feather(output_dir / REGIONS_FEATHER, **feather_kwargs)
if write_tsv:
    region_df.to_csv(output_dir / REGIONS_TSV, sep="\t", index=False)
```

Constants at module level:
```python
REGIONS_FEATHER = "regions.feather"
REGIONS_TSV = "regions.tsv"
```

## Integration into `TranscriptIndex.load()`

Add after the existing interval/SJ loading:

```python
# -- calibration regions -------------------------------------------
regions_path = os.path.join(index_dir, REGIONS_FEATHER)
if os.path.exists(regions_path):
    self.region_df = pd.read_feather(regions_path)
    region_cr = _cgranges_cls()
    # vectorised iteration
    ...
    region_cr.index()
    self.region_cr = region_cr
else:
    self.region_df = None
    self.region_cr = None
```

Add to `__init__`:
```python
self.region_df: pd.DataFrame | None = None
self.region_cr: _cgranges_cls | None = None
```

## Test plan

Uses the existing `mini_index` fixture from `conftest.py`.

**Mini-GTF layout (0-based half-open):**
```
chr1 (2000 bp):
  g1 (+): t0 exons (99,200),(299,400),(499,600)   span [99,600)
          t1 exons (99,200),(499,600)              span [99,600)
  g2 (-): t2 exons (999,1100),(1199,1300)          span [999,1300)
chr2 (500 bp): no genes
```

**Expected region boundaries on chr1 (sorted unique from all coordinates + 0 + 2000):**
```
[0, 99, 200, 299, 400, 499, 600, 999, 1100, 1199, 1300, 2000]
```

**Expected 11 regions on chr1:**

| # | start | end | exon_pos | exon_neg | tx_pos | tx_neg | context |
|---|---|---|---|---|---|---|---|
| 0 | 0 | 99 | F | F | F | F | intergenic |
| 1 | 99 | 200 | T | F | T | F | exon+ |
| 2 | 200 | 299 | F | F | T | F | intron+ (t0 only) |
| 3 | 299 | 400 | T | F | T | F | exon+ (t0 only) |
| 4 | 400 | 499 | F | F | T | F | intron+ (t0 only) |
| 5 | 499 | 600 | T | F | T | F | exon+ |
| 6 | 600 | 999 | F | F | F | F | intergenic |
| 7 | 999 | 1100 | F | T | F | T | exon− |
| 8 | 1100 | 1199 | F | F | F | T | intron− |
| 9 | 1199 | 1300 | F | T | F | T | exon− |
| 10 | 1300 | 2000 | F | F | F | F | intergenic |

Plus 1 region on chr2: `[0, 500)` all-False (intergenic).

**Test cases:**

1. **Non-overlap**: no two regions share any base pair.
2. **Reference-complete**: total region length per ref == ref_length.
3. **Correct region count**: 12 total (11 chr1 + 1 chr2).
4. **Flag correctness**: spot-check each region.
5. **Intergenic**: regions with all flags False.
6. **Pure exonic**: regions with `exon_*` True.
7. **Intronic**: regions with `tx_pos=True`, `exon_pos=False`.
8. **Negative strand**: regions with neg flags set.
9. **Strand-ambiguous region**: separate test with overlapping pos/neg
   transcripts.
10. **Empty reference**: chr2 produces one intergenic region.
11. **cgranges lookup**: query a coordinate and verify the correct `region_id`.
12. **`region_id` sequencing**: IDs are sequential starting at 0.
13. **`exon_pos` → `tx_pos` consistency**: wherever `exon_pos` is True, `tx_pos`
    must also be True (and likewise for neg).

## Step-by-step coding sequence

1. Add `REGIONS_FEATHER` / `REGIONS_TSV` constants.
2. Implement `build_region_table(transcripts, ref_lengths) -> pd.DataFrame`.
3. Wire `build_region_table` into `TranscriptIndex.build()`.
4. Add `region_df` and `region_cr` attributes to `TranscriptIndex.__init__()`.
5. Wire region loading into `TranscriptIndex.load()`.
6. Write `tests/test_regions.py` with the test cases above.
7. Run tests to confirm correctness.

## Not in first pass (per revision plan)

- Region-to-transcript sidecar mapping
- Tiny-region merging
- Fragment-to-region overlap summarization (Workstream B/C)
- Any `_ambig` columns
