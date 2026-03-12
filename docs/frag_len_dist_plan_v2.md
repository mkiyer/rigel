# Plan: Pool-Specific Fragment Length Distributions (gDNA vs RNA) — v2

## Context

Rigel separates gDNA from RNA using splice status and strand specificity.
Fragment length is currently unused for discrimination between gDNA and RNA —
a single `global_model` is used for all pools. Bioanalyzer/tapestation
confirms gDNA and RNA have distinct insert size distributions, providing
an orthogonal discriminating axis.

**Goal**: Maintain two frozen fragment length distributions (gDNA and RNA)
that produce `P(gDNA|length)` and `P(RNA|length)` likelihoods for scoring.
During scoring, gDNA candidate hits use the gDNA distribution and RNA
candidate hits use the RNA distribution, so the EM can combine
fragment-length signal with strand and geometry evidence.

## Background: Fragment Types and Trustworthiness

A **fragment** is constructed from BAM records. To build reliable
training distributions we must carefully define which fragments
contribute and which are held out.

### Splice Status

| Category | Criteria | splice_type | Training use |
|----------|----------|-------------|--------------|
| Annotated spliced | CIGAR contains N-ops matching known splice junctions in the GTF | `SPLICED_ANNOT (2)` | **Gold-standard RNA anchor** — train RNA model |
| Unannotated spliced | CIGAR contains N-ops that do NOT match any known splice junction | `SPLICED_UNANNOT (1)` | **Hold out** — could be alignment artifacts or novel junctions; origin unclear |
| Unspliced | No N-ops in any BAM record | `SPLICE_UNSPLICED (0)` | Ambiguous — could be gDNA, nRNA, or unspliced mRNA. **Weighted by strand probability for gDNA/RNA training** |

### Gap Splice Junctions (Implicit Introns) and Fragment Length Ambiguity

When paired-end reads don't overlap, there is an **unsequenced interval**
(gap) between the exon blocks of R1 and R2. The existing
`compute_frag_lengths()` function in `resolve_context.h` (lines 547–633)
queries a **SJ gap index** (an interval tree of annotated introns) to
find introns that fit inside each gap.

If an unspliced fragment has gaps containing annotated introns, the
per-transcript fragment lengths will differ (each transcript may subtract
different implicit introns from the genomic footprint). When
`get_unique_frag_length()` finds that candidate transcripts disagree on
the fragment length, it returns -1 (ambiguous). These fragments are
**already excluded** from `fraglen_obs.lengths` by the existing check
`if (ufl > 0)`.

**Implication for this plan**: Unspliced fragments whose per-transcript
fragment lengths disagree (due to potential implicit introns in the
unsequenced interval) are naturally held out from fragment length
distribution training. They have ambiguous provenance — the fragment
could be gDNA (no introns, true insert = genomic footprint) or RNA
(with implicit introns, true insert = genomic footprint minus intron
size). Only fragments where all candidate transcripts agree on the
fragment length — meaning either (a) the fragment is a single exon
block with no gap, or (b) all candidates agree on the gap correction —
are used for training. These are the fragments whose insert size we
can trust.

### Chimeras

Chimeric fragments (`chimera_type != CHIMERA_NONE`) do not have their
`frag_length_map` populated (see `_resolve_core()`). They will have
`get_unique_frag_length() == -1` and are naturally excluded from
training. No special handling needed.

### Intergenic Fragments

Intergenic fragments (those that fail resolution — `!resolved`) with no
introns are ~100% gDNA. Their `genomic_footprint` serves as the fragment
length. These are collected separately in `fraglen_obs.intergenic_lengths`
and currently feed both the intergenic model and `gdna_model`.

### Summary: Fragment Inclusion for Distribution Training

| Fragment class | RNA model | gDNA model | Hold out |
|---------------|-----------|------------|----------|
| `SPLICED_ANNOT` with unambiguous length | ✅ (anchor) | | |
| `SPLICED_UNANNOT` | | | ✅ |
| `UNSPLICED` with unambiguous length, sense strand | weighted `(1 - W_sense)` | weighted `W_sense` | |
| `UNSPLICED` with unambiguous length, antisense strand | weighted `(1 - W_anti)` | weighted `W_anti` | |
| `UNSPLICED` with ambiguous length (gap introns disagree) | | | ✅ (automatic) |
| Chimeric (`chimera_type != NONE`) | | | ✅ (automatic) |
| Intergenic unspliced | | ✅ (anchor) | |
| Intergenic spliced | | | ✅ (automatic) |

## The Challenge

Building the RNA model is straightforward: `SPLICED_ANNOT` fragments are
gold-standard RNA.

Building the gDNA model is harder because most unspliced fragments are
ambiguous (could be gDNA, nRNA, or unspliced mRNA). We need a rough
prior estimate of the probability that each fragment is gDNA, then use
that probability as a fractional weight when training the gDNA histogram.

**Available signals for weighting (without running the EM)**:
1. **Intergenic unspliced** → weight 1.0 (100% gDNA)
2. **Antisense unspliced** → high gDNA probability (in stranded data,
   antisense fragments at genic loci are predominantly gDNA)
3. **Sense unspliced** → mixed gDNA/RNA; weight depends on strand
   specificity and local gDNA/RNA ratio

## Approach: Strand-Weighted Histogram Mixing

We know that RNA obeys strand specificity ($s_{RNA}$) while gDNA is
perfectly unstranded (0.5). For unspliced fragments at genic loci, we
observe the strand direction relative to the gene (sense vs antisense)
using the existing `StrandModel.strand_likelihood()` infrastructure.

For a pool of S sense and A antisense unspliced fragment counts:

$$R = \max\!\bigl(0,\; \frac{S - A}{2 \cdot s_{RNA} - 1}\bigr)$$

$$G = (S + A) - R$$

$$W_{\text{sense,gDNA}} = \frac{G \cdot 0.5}{S}$$

$$W_{\text{anti,gDNA}} = \frac{G \cdot 0.5}{A}$$

These weights are the estimated fraction of gDNA in each strand pool.
Applied to per-strand histograms, they build the final models:

$$H_{\text{gDNA}} = H_{\text{intergenic}} + W_{\text{sense}} \cdot H_{\text{sense}} + W_{\text{anti}} \cdot H_{\text{anti}}$$

$$H_{\text{RNA}} = H_{\text{spliced}} + (1 - W_{\text{sense}}) \cdot H_{\text{sense}} + (1 - W_{\text{anti}}) \cdot H_{\text{anti}}$$

**Properties**:
- Zero parameters — no thresholds or cutoffs
- Graceful degradation: when $s_{RNA} \to 0.5$, all $W \to 0$, and the
  gDNA model degrades to intergenic-only while the RNA model gains the
  spliced anchor plus all unspliced observations. This is the correct
  behavior: with no strand information, we cannot distinguish gDNA from
  RNA by strand, so the gDNA model is the intergenic distribution and
  the RNA model is the broader population.
- Uses ALL qualifying unspliced fragments, weighted by their gDNA probability
- O(1) mixing: numpy operations on histogram arrays of size `max_size + 1`

### Sense vs Antisense: Mapping from Same/Opposite Strand

In the C++ scanner, we observe `exon_strand` (R1 alignment strand) and
`gene_strand` (annotated transcript strand from `t_strand_arr`).

The existing `StrandModel.strand_likelihood(exon_strand, gene_strand)`
function (`strand_model.py` line 239) returns:
- `p_r1_sense` when `exon_strand == gene_strand` (same strand)
- `p_r1_antisense` (= `1 - p_r1_sense`) when they differ (opposite strand)

For the purpose of this plan, "same strand" and "opposite strand" are
what the C++ scanner directly observes. The mapping to biological
"sense" vs "antisense" depends on the library protocol:

- **R1-sense** (`p_r1_sense > 0.5`): same_strand ≈ sense, opp_strand ≈ antisense
- **R1-antisense** (`p_r1_sense < 0.5`): same_strand ≈ antisense, opp_strand ≈ sense

In `mix_models()`, we swap which histogram is `H_sense` and which
is `H_anti` based on `p_r1_sense`:

```python
if p_r1_sense >= 0.5:
    H_sense = H_same_strand
    H_anti  = H_opp_strand
else:
    H_sense = H_opp_strand
    H_anti  = H_same_strand
```

## Implementation Steps

### Step 1: Collect same/opp strand unspliced histograms (C++)

We need 4 raw histograms (2 already exist, 2 new):

| Histogram | Source | Status |
|-----------|--------|--------|
| `H_spliced` | `category_models[SPLICED_ANNOT]` | **Exists** |
| `H_intergenic` | `intergenic` model | **Exists** |
| `H_unspliced_same_strand` | Unspliced unique-mapper genic, `exon_strand == gene_strand` | **NEW** |
| `H_unspliced_opp_strand` | Unspliced unique-mapper genic, `exon_strand != gene_strand` | **NEW** |

#### C++ changes: `bam_scanner.cpp`

**`FragLenObservations` struct** (line 206): Add two new vectors.

```cpp
struct FragLenObservations {
    std::vector<int32_t> lengths;
    std::vector<int32_t> splice_types;
    std::vector<int32_t> intergenic_lengths;
    std::vector<int32_t> unspliced_same_strand_lengths;  // NEW
    std::vector<int32_t> unspliced_opp_strand_lengths;   // NEW
};
```

**Collection logic** (after line 1280): After the existing
`get_unique_frag_length()` block, add collection for the new vectors.

Qualifying criteria for a fragment to be collected into the new vectors:
1. `is_unique_mapper` — already inside the unique-mapper guard
2. `result.splice_type == SPLICE_UNSPLICED` — no explicit splice junctions
3. Not chimeric — already guaranteed because chimeric fragments have
   empty `frag_length_map` and `genomic_footprint` won't be used
4. `result.get_is_same_strand()` — unambiguous strand among candidate genes
5. `result.exon_strand` is POS or NEG (informative)
6. `gene_strand` (`ctx.t_strand_arr_[get_first_t_ind()]`) is POS or NEG

**What fragment length to use**: `genomic_footprint`. For unspliced
fragments intended to train the gDNA distribution, the genomic footprint
is the correct measurement — it represents the true insert size assuming
no internal introns. For the RNA component of unspliced, the footprint
is also adequate because the mixing weight dilutes these into a
distribution already anchored by the high-quality spliced fragments.

Fragments with ambiguous per-transcript lengths (gap introns that
differ across candidates) are included here using `genomic_footprint`
because:
- The gDNA interpretation of these fragments has no introns →
  `genomic_footprint` is their true length
- The RNA interpretation would need per-transcript correction, but these
  fragments contribute only a small fraction of RNA weight, and the
  spliced anchor dominates the RNA distribution

```cpp
// After existing frag_length recording (line ~1280), add:
if (result.splice_type == SPLICE_UNSPLICED &&
    result.get_is_same_strand() &&
    (result.exon_strand == STRAND_POS ||
     result.exon_strand == STRAND_NEG)) {
    int32_t t_idx = result.get_first_t_ind();
    if (t_idx >= 0 &&
        t_idx < static_cast<int32_t>(ctx.t_strand_arr_.size())) {
        int32_t gene_strand = ctx.t_strand_arr_[t_idx];
        if (gene_strand == STRAND_POS ||
            gene_strand == STRAND_NEG) {
            int32_t gfp = result.genomic_footprint;
            if (gfp > 0) {
                if (result.exon_strand == gene_strand) {
                    fraglen_obs.unspliced_same_strand_lengths
                        .push_back(gfp);
                } else {
                    fraglen_obs.unspliced_opp_strand_lengths
                        .push_back(gfp);
                }
            }
        }
    }
}
```

**`merge_fraglen_obs()`** (line 413): Add merge for the 2 new vectors:

```cpp
dst.unspliced_same_strand_lengths.insert(
    dst.unspliced_same_strand_lengths.end(),
    src.unspliced_same_strand_lengths.begin(),
    src.unspliced_same_strand_lengths.end());
dst.unspliced_opp_strand_lengths.insert(
    dst.unspliced_opp_strand_lengths.end(),
    src.unspliced_opp_strand_lengths.begin(),
    src.unspliced_opp_strand_lengths.end());
```

**Python dict export** (line 1355): Add 2 new keys to `fraglen_dict`:

```cpp
fraglen_dict["unspliced_same_strand_lengths"] =
    std::move(fraglen_obs_.unspliced_same_strand_lengths);
fraglen_dict["unspliced_opp_strand_lengths"] =
    std::move(fraglen_obs_.unspliced_opp_strand_lengths);
```

#### Files modified

- bam_scanner.cpp

---

### Step 2: Add new models and `mix_models()` to `FragmentLengthModels`

#### Python changes: `frag_length_model.py`

**`FragmentLengthModels.__init__()`**: Add 3 new `FragmentLengthModel`
instances:

```python
self.unspliced_same_strand = FragmentLengthModel(max_size=max_size)
self.unspliced_opp_strand = FragmentLengthModel(max_size=max_size)
self.rna_model = FragmentLengthModel(max_size=max_size)
```

The existing `gdna_model` is repurposed to hold the mixed gDNA result
(it currently holds intergenic-only observations; after mixing it will
hold intergenic + weighted unspliced).

**New `mix_models()` method**: Performs the histogram mixing math. This
is called AFTER `strand_models.finalize()` (so strand specificity is
known) and BEFORE `frag_length_models.finalize()` (so mixed models get
their log-probability lookup tables cached).

```python
def mix_models(self, s_rna: float, p_r1_sense: float) -> None:
    """Build gDNA and RNA distributions via strand-weighted mixing.

    Parameters
    ----------
    s_rna : float
        Strand specificity (0.5 = unstranded, 1.0 = perfectly
        stranded). From ``strand_models.strand_specificity``.
    p_r1_sense : float
        P(read 1 aligns in gene-sense direction).
        From ``strand_models.exonic_spliced.p_r1_sense``.
    """
    from .splice import SpliceType

    # 1. Map same/opp strand to sense/antisense
    if p_r1_sense >= 0.5:
        H_sense = self.unspliced_same_strand.counts.copy()
        H_anti  = self.unspliced_opp_strand.counts.copy()
    else:
        H_sense = self.unspliced_opp_strand.counts.copy()
        H_anti  = self.unspliced_same_strand.counts.copy()

    # 2. Total counts in each strand bin
    S = float(H_sense.sum())
    A = float(H_anti.sum())

    # 3. Estimate total RNA and gDNA in the unspliced pool
    denom = 2.0 * s_rna - 1.0
    if denom > 1e-6:
        R = max(0.0, (S - A) / denom)
        G = (S + A) - R
    else:
        # No strand information — cannot separate
        R = 0.0
        G = 0.0

    # 4. Per-strand gDNA weights (fraction of gDNA in each bin)
    W_sense = min(max((G * 0.5 / S) if S > 0 else 0.0, 0.0), 1.0)
    W_anti  = min(max((G * 0.5 / A) if A > 0 else 0.0, 0.0), 1.0)

    # 5. Build mixed histograms
    H_intergenic = self.intergenic.counts
    H_spliced = self.category_models[SpliceType.SPLICED_ANNOT].counts

    gdna_counts = (
        H_intergenic
        + W_sense * H_sense
        + W_anti * H_anti
    )
    rna_counts = (
        H_spliced
        + (1.0 - W_sense) * H_sense
        + (1.0 - W_anti) * H_anti
    )

    # 6. Write into models (replace counts and total_weight)
    self.gdna_model.counts = gdna_counts
    self.gdna_model._total_weight = float(gdna_counts.sum())

    self.rna_model.counts = rna_counts
    self.rna_model._total_weight = float(rna_counts.sum())

    # 7. Log diagnostics
    logger.info(
        f"Fragment length mixing: S={S:.0f} sense, A={A:.0f} anti, "
        f"s_rna={s_rna:.4f}, W_sense={W_sense:.4f}, "
        f"W_anti={W_anti:.4f}"
    )
    logger.info(
        f"  gDNA model: {self.gdna_model.total_weight:.0f} "
        f"weighted obs, mean={self.gdna_model.mean:.1f}"
    )
    logger.info(
        f"  RNA model: {self.rna_model.total_weight:.0f} "
        f"weighted obs, mean={self.rna_model.mean:.1f}"
    )
```

**Update `finalize()`**: Also finalize `rna_model`,
`unspliced_same_strand`, and `unspliced_opp_strand`.

**Update `to_dict()` and `log_summary()`**: Include new models in
serialization and logging.

#### Files modified

- frag_length_model.py

---

### Step 3: Replay new observation vectors in pipeline

#### Python changes: pipeline.py

**`_replay_fraglen_observations()`** (line 183): Add replay for the 2
new vectors. These feed into the raw `unspliced_same_strand` and
`unspliced_opp_strand` models on `FragmentLengthModels`:

```python
same_strand = fraglen_dict.get("unspliced_same_strand_lengths", [])
if len(same_strand) > 0:
    frag_length_models.unspliced_same_strand.observe_batch(
        np.asarray(same_strand, dtype=np.intp),
    )

opp_strand = fraglen_dict.get("unspliced_opp_strand_lengths", [])
if len(opp_strand) > 0:
    frag_length_models.unspliced_opp_strand.observe_batch(
        np.asarray(opp_strand, dtype=np.intp),
    )
```

#### Files modified

- pipeline.py

---

### Step 4: Wire `mix_models()` into the pipeline

#### Python changes: pipeline.py

**Between `strand_models.finalize()` and `frag_length_models.finalize()`**
(lines 735–736): Insert the mixing call:

```python
strand_models.finalize()
frag_length_models.mix_models(
    s_rna=strand_models.strand_specificity,
    p_r1_sense=strand_models.exonic_spliced.p_r1_sense,
)
frag_length_models.finalize()
```

This ordering is critical:
1. `strand_models.finalize()` — caches `p_r1_sense` and
   `strand_specificity`
2. `frag_length_models.mix_models()` — builds mixed histograms using
   strand info
3. `frag_length_models.finalize()` — caches log-probability lookup
   tables for all models including the mixed ones

#### Files modified

- pipeline.py

---

### Step 5: Add gDNA lookup table to `FragmentScorer` and C++ `NativeFragmentScorer`

Currently, a single fragment length lookup table (from `global_model`)
is used for both RNA and gDNA scoring. We need two separate lookup
tables: one from the RNA distribution and one from the gDNA distribution.

A "lookup table" here means the pre-computed `_log_prob` numpy array
(shape `(max_size + 1,)`) cached by `FragmentLengthModel.finalize()`.
Scoring a fragment length is a single array index: `log_prob[flen]`.
Lengths beyond `max_size` receive exponential tail decay.

#### Python changes: scoring.py

**`FragmentScorer` dataclass** (line 88): Rename existing fields for
clarity and add 3 new fields for gDNA:

```python
# Fragment-length lookup table: RNA distribution
fl_log_prob: np.ndarray | None
fl_max_size: int
fl_tail_base: float

# Fragment-length lookup table: gDNA distribution  (NEW)
gdna_fl_log_prob: np.ndarray | None
gdna_fl_max_size: int
gdna_fl_tail_base: float
```

**`from_models()` classmethod** (line 160): Change the lookup table
extraction:

```python
# RNA lookup table: from rna_model (mixed RNA distribution)
rna_fl = frag_length_models.rna_model
fl_log_prob = rna_fl._log_prob
fl_max_size = rna_fl.max_size
fl_tail_base = getattr(rna_fl, "_tail_base", 0.0)

# gDNA lookup table: from gdna_model (mixed gDNA distribution)
gdna_fl = frag_length_models.gdna_model
gdna_fl_log_prob = gdna_fl._log_prob
gdna_fl_max_size = gdna_fl.max_size
gdna_fl_tail_base = getattr(gdna_fl, "_tail_base", 0.0)
```

Both lookup tables are always populated — no fallback to `global_model`
needed. When strand information is poor ($s_{RNA} \to 0.5$), the mixing
math naturally degrades:
- `gdna_model` contains only intergenic observations (the best pure-gDNA
  anchor we have)
- `rna_model` contains spliced + all unspliced observations (the
  broadest RNA estimate)

This is the correct behavior for an uninformative strand scenario.

Pass the new gDNA fields to the `FragmentScorer` constructor and the
`NativeFragmentScorer` C++ constructor.

#### C++ changes: scoring.cpp

**`NativeFragmentScorer` class** (line 130): Add member variables:

```cpp
// gDNA fragment-length lookup table
std::vector<double> gdna_fl_log_prob_;
int32_t gdna_fl_max_size_;
double  gdna_fl_tail_base_;
bool    has_gdna_fl_lut_;
```

**New inline method** (after `frag_len_log_lik()` at line 179):

```cpp
inline double gdna_frag_len_log_lik(int32_t flen) const {
    if (flen <= 0 || !has_gdna_fl_lut_) return 0.0;
    if (flen <= gdna_fl_max_size_)
        return gdna_fl_log_prob_[static_cast<size_t>(flen)];
    return gdna_fl_tail_base_
           + (flen - gdna_fl_max_size_) * TAIL_DECAY_LP;
}
```

**Constructor** (line 224): Add 3 new parameters
(`gdna_fl_log_prob_obj`, `gdna_fl_max_size`, `gdna_fl_tail_base`).
Initialize and copy gDNA lookup table into C++ owned storage (same
pattern as existing RNA lookup table copy at line 296):

```cpp
// Copy gDNA fragment-length lookup table
gdna_fl_max_size_ = gdna_fl_max_size;
gdna_fl_tail_base_ = gdna_fl_tail_base;
has_gdna_fl_lut_ = false;
if (!gdna_fl_log_prob_obj.is_none()) {
    auto arr = nb::cast<f64_1d>(gdna_fl_log_prob_obj);
    const double* p = arr.data();
    int32_t n = static_cast<int32_t>(arr.shape(0));
    gdna_fl_log_prob_.assign(p, p + n);
    has_gdna_fl_lut_ = true;
}
```

**Nanobind bindings** (line 1496): Add the 3 new constructor parameters
to the `nb::init<>` template and `nb::arg()` declarations:

```cpp
nb::arg("gdna_fl_log_prob").none(),
nb::arg("gdna_fl_max_size"),
nb::arg("gdna_fl_tail_base"),
```

#### Files modified

- scoring.py
- scoring.cpp

---

### Step 6: Update gDNA scoring to use gDNA lookup table

Currently, gDNA scoring in C++ uses `frag_len_log_lik()` (the RNA
lookup table) to compute `P(fragment_length)` for gDNA candidates.
This should use the gDNA-specific lookup table instead.

#### C++ changes: scoring.cpp

**Multimapper gDNA scoring** (line 770):
```cpp
// BEFORE:
double gdna_fl = frag_len_log_lik(gfp_val);

// AFTER:
double gdna_fl = gdna_frag_len_log_lik(gfp_val);
```

**Single-hit gDNA scoring** (line 1218):
```cpp
// BEFORE:
double gdna_fl = frag_len_log_lik(genomic_footprint);

// AFTER:
double gdna_fl = gdna_frag_len_log_lik(genomic_footprint);
```

The RNA scoring path (`frag_len_log_lik()`) remains unchanged — it
continues to use the RNA lookup table for mRNA and nRNA candidate hits.

#### Test helper changes: scoring_helpers.py

**`score_gdna_standalone()`** (line 89): Change to use `gdna_model`:

```python
# BEFORE:
log_p_insert = (
    frag_length_models.global_model.log_likelihood(frag_length)
    if frag_length > 0
    else 0.0
)

# AFTER:
log_p_insert = (
    frag_length_models.gdna_model.log_likelihood(frag_length)
    if frag_length > 0
    else 0.0
)
```

#### Files modified

- scoring.cpp
- scoring_helpers.py

---

### Step 7: Logging and tests

#### Logging

`mix_models()` already includes diagnostic logging (Step 2).
Additionally:

- **`log_summary()`** in `FragmentLengthModels`: Add entries for
  `unspliced_same_strand`, `unspliced_opp_strand`, `rna_model`,
  and updated `gdna_model` stats.
- **`to_dict()`**: Include the new models for JSON serialization.

#### New test file: `tests/test_gdna_frag_length.py`

Test cases:

1. **Mixing math correctness**: Construct known same/opp strand
   histograms and verify that `mix_models()` produces the expected
   gDNA and RNA histogram counts.

2. **Strand specificity edge cases**:
   - $s_{RNA} = 1.0$ (perfect stranding): All antisense → gDNA, all
     sense → RNA. $W_{anti} = 1.0$, $W_{sense} = 0.0$.
   - $s_{RNA} = 0.5$ (no strand info): $W_{sense} = W_{anti} = 0$.
     gDNA = intergenic only. RNA = spliced + all unspliced.
   - $s_{RNA} = 0.75$ (typical): Intermediate weights.

3. **R1-sense vs R1-antisense**: Verify that `mix_models()` correctly
   swaps same/opp to sense/anti based on `p_r1_sense`.

4. **Graceful handling of empty models**: Verify that mixing works when
   `unspliced_same_strand` or `unspliced_opp_strand` has zero
   observations (e.g., single-end data or very sparse data).

5. **gDNA scoring uses gDNA lookup table**: Verify that
   `score_gdna_standalone()` uses `gdna_model` and returns different
   values than the RNA model for the same fragment length.

6. **End-to-end**: Build complete `FragmentLengthModels`, call
   `mix_models()`, `finalize()`, then verify the lookup tables on
   `rna_model` and `gdna_model` are distinct and correctly normalized.

#### Files modified / created

- frag_length_model.py (logging updates)
- `tests/test_gdna_frag_length.py` (**NEW**)

---

## Files Modified Summary

| File | Changes |
|------|---------|
| bam_scanner.cpp | 2 new vectors in `FragLenObservations`; collection logic for unspliced same/opp strand; merge; Python dict export |
| frag_length_model.py | `unspliced_same_strand`, `unspliced_opp_strand`, `rna_model` models; `mix_models()` method; updated `finalize()`, `to_dict()`, `log_summary()` |
| pipeline.py | Replay new observation vectors; call `mix_models()` between strand and frag length finalization |
| scoring.py | 3 new gDNA lookup table fields on `FragmentScorer`; `from_models()` extracts RNA and gDNA lookup tables separately; passes both to C++ |
| scoring.cpp | 4 new member variables; `gdna_frag_len_log_lik()` method; updated constructor + bindings; gDNA scoring uses gDNA lookup table |
| scoring_helpers.py | `score_gdna_standalone()` uses `gdna_model` |
| `tests/test_gdna_frag_length.py` | **NEW** — mixing math, edge cases, scoring integration |

## Verification

1. `pytest tests/` — existing tests pass
2. `pytest tests/test_gdna_frag_length.py` — mixing math + scoring tests
3. Run on stranded data — log shows distinct gDNA and RNA distributions
4. Compare quantification with/without — no regression, improved gDNA separation
```

Please save this to [docs/frag_len_dist_plan_v2.md](docs/frag_len_dist_plan_v2.md). I don't have file creation tools enabled to write it directly.

Key changes from v1:

1. **Careful fragment type definitions** — detailed treatment of which fragments are trustworthy for training, how gap splice junctions cause automatic exclusion via `get_unique_frag_length() == -1`, and why chimeras and `SPLICED_UNANNOT` are held out
2. **`genomic_footprint`** used for same/opp strand vectors (not `unique_frag_length`) — this is the true insert size for gDNA fragments and captures more data points
3. **C++ is the primary scoring path** — all scoring changes target `scoring.cpp`; the Python `score_gdna_standalone()` in `scoring_helpers.py` is test-only
4. **No fallback to `global_model`** — the mixing math naturally degrades when strand info is poor, producing correct distributions at all strand specificity levels
5. **Strand mapping** uses `p_r1_sense` to swap same/opp → sense/anti, leveraging the existing `StrandModel.strand_likelihood()` semantics
6. **"LUT" replaced** with descriptive "lookup table" terminology throughoutPlease save this to [docs/frag_len_dist_plan_v2.md](docs/frag_len_dist_plan_v2.md). I don't have file creation tools enabled to write it directly.

Key changes from v1:

1. **Careful fragment type definitions** — detailed treatment of which fragments are trustworthy for training, how gap splice junctions cause automatic exclusion via `get_unique_frag_length() == -1`, and why chimeras and `SPLICED_UNANNOT` are held out
2. **`genomic_footprint`** used for same/opp strand vectors (not `unique_frag_length`) — this is the true insert size for gDNA fragments and captures more data points
3. **C++ is the primary scoring path** — all scoring changes target `scoring.cpp`; the Python `score_gdna_standalone()` in `scoring_helpers.py` is test-only
4. **No fallback to `global_model`** — the mixing math naturally degrades when strand info is poor, producing correct distributions at all strand specificity levels
5. **Strand mapping** uses `p_r1_sense` to swap same/opp → sense/anti, leveraging the existing `StrandModel.strand_likelihood()` semantics
6. **"LUT" replaced** with descriptive "lookup table" terminology throughout