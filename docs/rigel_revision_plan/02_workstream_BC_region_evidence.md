# Workstream B/C: Regional Evidence Extraction

This is the most important early implementation target because it produces the
observable substrate on which every later purity approximation will depend.

## 1. Why This Is Feasible Now

The current scan path already exposes fragment-level fields needed for the first
pass.

From the buffer and scan pipeline we already have access to:

- fragment splice state
- fragment strand signal
- genomic start
- genomic footprint
- exonic overlap summary
- intronic overlap summary
- unambiguous intronic overlap summary

That means the first revision can accumulate region evidence without waiting for
the later EM redesign.

However, one design constraint needs to be made explicit: the buffered fragment
representation retains `genomic_start` and `genomic_footprint`, but not the
full aligned block geometry for spliced fragments. So exact atomic-region
overlap by genomic footprint is valid immediately for unspliced fragments, but
not for spliced fragments whose skipped introns would be falsely counted as
covered bases.

## 2. Implementation Goal

Build a standalone per-region evidence table from parsed fragments.

The first-pass evidence table should include at least:

- `region_id`
- total unspliced fragment count
- pos-strand unspliced count
- neg-strand unspliced count
- spliced evidence count in or near the region
- total fragment-length summary for calibration-eligible unspliced fragments
- effective-length-normalized density summaries by context

The first pass should also define an internal merged fragment-region summary for
unspliced fragments, because one fragment will often overlap multiple atomic
regions.

## 3. Recommended First-Pass Admissibility Rule

Use a conservative calibration subset first.

Recommended rule:

1. only unspliced fragments contribute directly to gDNA calibration evidence
2. spliced fragments are used only as RNA evidence, not as direct gDNA
   calibration observations
3. admissibility should be defined on the merged fragment-region summary, not on
  the raw count of atomic regions overlapped
4. context-ambiguous merged summaries are excluded from the first calibration
  pass

This is deliberately conservative. The priority is to build a trustworthy
signal, not to maximize calibration mass in the first implementation.

## 4. Suggested Data Flow

### 4.1 Region lookup structure

Build a query structure keyed by reference and interval coordinates. The current
code already uses cgranges-style interval lookup elsewhere, so the same pattern
is appropriate here.

### 4.2 Fragment mapping

For each fragment, compute the overlapped calibration regions from:

- `genomic_start`
- `genomic_footprint`

Derive the genomic interval as:

- `[genomic_start, genomic_start + genomic_footprint)`

For unspliced fragments, this is sufficient for first-pass overlap mapping.

For spliced fragments, do not use genomic footprint as if it were covered
sequence. Either:

- accumulate splice evidence during scan before block geometry is discarded, or
- treat spliced evidence as a separate local RNA signal not requiring exact
  atomic-region coverage in the first implementation

Then merge the raw region hits into one fragment-region summary according to the
admissibility rule.

### 4.3 Evidence accumulation

Maintain a per-region accumulator with fields such as:

- `n_unspliced`
- `n_unspliced_pos`
- `n_unspliced_neg`
- `n_spliced_local`
- `frag_length_sum`
- `frag_length_weight`
- optional context-specific density numerators

The exact schema can expand later. The important point is to keep the first
accumulator transparent and inspectable.

The fragment-region merge helper should also expose layer-specific overlap
summaries such as:

- `exon_pos_bp`
- `exon_neg_bp`
- `exon_ambig_bp`
- `tx_pos_bp`
- `tx_neg_bp`
- `tx_ambig_bp`

These are overlap totals on annotation layers, not disjoint bins. Consumers
must not sum exon and transcript-span overlap totals as though they were a
partition of bases.

## 5. Evidence Channels in the First Pass

### 5.1 Splice evidence

Initial recommendation:

- count spliced fragments overlapping the region footprint directly
- in the first implementation, prefer a separate scan-time or neighborhood-based
  splice evidence field rather than pretending genomic footprint yields exact
  spliced-region coverage

### 5.2 Density evidence

Initial recommendation:

- compute region-level unspliced density as `n_unspliced / effective_length`
- stratify interpretation by annotation context rather than forcing a single
  density score across all regions immediately

### 5.3 Strand evidence

Initial recommendation:

- store raw pos/neg counts first
- do not transform them into a symmetry score inside the evidence extractor
- leave Beta-Binomial modeling to the calibration layer

## 6. Step-by-Step Coding Plan

1. define a `RegionEvidence` schema or DataFrame output contract
2. define a merged `FragmentRegionSummary` contract for unspliced fragments
3. add a region-lookup helper for the calibration-region table
4. add a scan-adjacent accumulator that consumes fragment objects or buffered
  fragment arrays
5. emit one standalone region evidence table per run
6. add diagnostics to inspect evidence totals by annotation context
7. validate on synthetic data with known gDNA contamination patterns

## 7. Concrete Code Touchpoints

Most likely touchpoints:

- `src/rigel/buffer.py`
- `src/rigel/scan.py`
- `src/rigel/pipeline.py`
- a new module such as `src/rigel/calibration_evidence.py`

The strongest design option is to keep this path independent of the native EM
solver so it can be debugged with ordinary Python-side diagnostics first.

It should also remain independent of the main transcript-resolution index in the
first pass. Region lookup is for calibration evidence, not for replacing the
resolver hot path immediately.

## 8. Recommended Output Before Any Purity Modeling

Before implementing a purity approximation, the pipeline should be able to emit:

1. a region table
2. a region evidence table
3. summary diagnostics by context class

Once those exist, the right purity approximation will be much easier to judge
empirically.