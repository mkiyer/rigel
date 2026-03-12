# Plan: Pool-Specific Fragment Length Distributions (gDNA vs RNA)

## Context

hulkrna separates gDNA from RNA using splice status and strand specificity. Fragment length is currently unused for discrimination — a single `global_model` is used for all pools. Bioanalyzer/tapestation confirms gDNA and RNA have distinct insert size distributions, providing an orthogonal discriminating axis.

**Goal**: Maintain two frozen fragment length distributions (gDNA and RNA) that produce P(gDNA|length) and P(RNA|length) likelihoods for the EM. The gDNA distribution doesn't need to be perfect — it needs to be *distinct* from the RNA distribution. The EM will combine this fragment-length signal with strand and geometry evidence.

## The Challenge

Building the RNA model is easy: SPLICED_ANNOT fragments are gold-standard RNA. Building the gDNA model is harder because most unspliced fragments are ambiguous (could be gDNA, nRNA, or even unspliced mRNA). We need a rough prior estimate of the probability that each fragment is gDNA, then use that probability as a fractional weight when training the gDNA histogram.

**Available signals for weighting (without running the EM)**:
1. **Intergenic unspliced** → weight 1.0 (100% gDNA, but can be very sparse)
2. **Antisense unspliced** → weight depends on strand specificity (SS=1.0 → all gDNA; SS=0.5 → uninformative)
3. **Sense unspliced** → weight depends on both SS and the local gDNA/RNA ratio
4. **Unexpressed loci** → effectively intergenic (no active transcription); a coverage-based signal for future enhancement

## Approach: Strand-Weighted Histogram Mixing

We know that RNA obeys strand specificity (s_RNA) while gDNA is perfectly unstranded (0.5). For any compartment with S sense and A antisense unspliced fragments:

```
R = max(0, (S - A) / (2·s_RNA - 1))     # estimated total RNA
G = (S + A) - R                           # estimated total gDNA
W_sense_gDNA  = (G · 0.5) / S            # gDNA weight for sense histogram
W_anti_gDNA   = (G · 0.5) / A            # gDNA weight for antisense histogram
```

These weights are applied to per-compartment histograms to build the final models:
```
H_gDNA = H_intergenic + W_sense·H_sense + W_anti·H_anti
H_RNA  = H_spliced + (1-W_sense)·H_sense + (1-W_anti)·H_anti
```

**Properties**:
- Zero parameters — no SS thresholds, no min_obs cutoffs
- Graceful degradation: when s_RNA → 0.5, all W → 0, models collapse to the safe anchors (H_intergenic and H_spliced)
- Uses ALL unspliced fragments, weighted by their gDNA probability
- O(1) mixing: numpy operations on histogram arrays of size 2001

## Implementation Steps

### Step 1: Collect sense/antisense unspliced histograms

We need 4 raw histograms (2 already exist, 2 new):
- `H_spliced` = `category_models[SPLICED_ANNOT]` — **exists**
- `H_intergenic` = `intergenic` — **exists**
- `H_unspliced_same_strand` — **NEW**: unspliced unique-mapper genic, `exon_strand == gene_strand`
- `H_unspliced_opp_strand` — **NEW**: unspliced unique-mapper genic, `exon_strand != gene_strand`

After strand model finalization, same/opp are mapped to sense/antisense based on `p_r1_sense`.

**Collection**: Add two new vectors to `FragLenObservations` in C++ BAM scanner. At the point where genic fragment lengths are recorded (~line 1050 of `bam_scanner.cpp`), also record `genomic_footprint` for unspliced unique-mapper fragments into same_strand or opp_strand vectors based on `exon_strand` vs `gene_strand`.

In Python `_replay_fraglen_observations()`, route into two new `FragmentLengthModel` instances.

**Files**: `src/hulkrna/native/bam_scanner.cpp`, `src/hulkrna/frag_length_model.py`, `src/hulkrna/pipeline.py`

### Step 2: Add `mix_models()` to `FragmentLengthModels`

New method that performs the histogram mixing math shown above. Also add `rna_model = FragmentLengthModel(...)` to `__init__`. The existing `gdna_model` is repurposed to hold the mixed result.

**File**: `src/hulkrna/frag_length_model.py`

### Step 3: Wire into pipeline

```python
strand_models.finalize()
frag_length_models.mix_models(
    s_rna=strand_models.strand_specificity,
    p_r1_sense=strand_models.exonic_spliced.p_r1_sense,
)
frag_length_models.finalize()  # caches LUTs for gdna_model, rna_model
```

**File**: `src/hulkrna/pipeline.py` (~line 747-748)

### Step 4: Add gDNA LUT fields to `FragmentScorer`

Add `gdna_fl_log_prob`, `gdna_fl_max_size`, `gdna_fl_tail_base` to the frozen dataclass.

**File**: `src/hulkrna/scoring.py`

### Step 5: Update `from_models()` to use mixed models

- RNA LUT: from `rna_model` (falls back to `global_model` if rna_model is empty)
- gDNA LUT: from `gdna_model` (falls back to `global_model` if gdna_model is empty)
- C++ `NativeFragmentScorer` automatically picks up the RNA LUT change (no C++ constructor changes needed)

**File**: `src/hulkrna/scoring.py`

### Step 6: Add `gdna_frag_len_log_lik()` and update gDNA scoring

- New function `gdna_frag_len_log_lik(ctx, flen)` reading the gDNA LUT fields
- `score_gdna_standalone()`: use `gdna_model` instead of `global_model`
- `scan.py:_gdna_log_lik()`: use `gdna_frag_len_log_lik()` instead of `frag_len_log_lik()`

**Files**: `src/hulkrna/scoring.py`, `src/hulkrna/scan.py`

### Step 7: Logging and tests

- Log computed weights (W_sense, W_anti), model stats (mean, mode, n_obs), and whether mixing was active
- New `tests/test_gdna_frag_length.py` covering mixing math, graceful fallback, and scoring

**Files**: `src/hulkrna/frag_length_model.py`, `src/hulkrna/pipeline.py`, `tests/test_gdna_frag_length.py`

## Files Modified

| File | Change |
|------|--------|
| `src/hulkrna/native/bam_scanner.cpp` | 2 new observation vectors (same/opp strand unspliced lengths) |
| `src/hulkrna/frag_length_model.py` | `rna_model`, `unspliced_same_strand`, `unspliced_opp_strand`; `mix_models()` |
| `src/hulkrna/scoring.py` | gDNA LUT fields, `gdna_frag_len_log_lik()`, update `from_models()` and `score_gdna_standalone()` |
| `src/hulkrna/scan.py` | `_gdna_log_lik()` → `gdna_frag_len_log_lik()` |
| `src/hulkrna/pipeline.py` | Replay new observations, call `mix_models()`, logging |
| `tests/test_gdna_frag_length.py` | **New** |

## Future Enhancements

- **Coverage-based weighting**: Loci with zero/minimal expression are functionally intergenic. Their unspliced fragments could receive weight ~1.0 for gDNA. This ties into the broader EM initialization problem and should be developed as part of a holistic initialization strategy.
- **Per-compartment separation**: Split unspliced into exonic vs intronic compartments for finer-grained gDNA fraction estimates (different gDNA/RNA ratios in exons vs introns).

## Verification

1. `pytest tests/` — existing tests pass
2. `pytest tests/test_gdna_frag_length.py` — mixing math + scoring tests
3. Run on stranded data — log shows distinct gDNA and RNA distributions
4. Compare quantification with/without — no regression, improved gDNA separation
