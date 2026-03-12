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

## 3. Recommended First-Pass Admissibility Rule

Use a conservative calibration subset first.

Recommended rule:

1. only unspliced fragments contribute directly to gDNA calibration evidence
2. only fragments with one interpretable calibration-region assignment are used
3. context-ambiguous assignments are excluded from the first calibration pass
4. spliced fragments are used only as RNA evidence, not as direct gDNA
   calibration observations

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

Then assign the fragment to regions according to the admissibility rule.

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

## 5. Evidence Channels in the First Pass

### 5.1 Splice evidence

Initial recommendation:

- count spliced fragments overlapping the region footprint directly
- if overlap is too sparse in practice, add a separate neighborhood-based splice
  evidence field later rather than hiding that logic inside the first version

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
2. add a region-lookup helper for the calibration-region table
3. add a scan-adjacent accumulator that consumes fragment objects or buffered
   fragment arrays
4. emit one standalone region evidence table per run
5. add diagnostics to inspect evidence totals by annotation context
6. validate on synthetic data with known gDNA contamination patterns

## 7. Concrete Code Touchpoints

Most likely touchpoints:

- `src/rigel/buffer.py`
- `src/rigel/scan.py`
- `src/rigel/pipeline.py`
- a new module such as `src/rigel/calibration_evidence.py`

The strongest design option is to keep this path independent of the native EM
solver so it can be debugged with ordinary Python-side diagnostics first.

## 8. Recommended Output Before Any Purity Modeling

Before implementing a purity approximation, the pipeline should be able to emit:

1. a region table
2. a region evidence table
3. summary diagnostics by context class

Once those exist, the right purity approximation will be much easier to judge
empirically.