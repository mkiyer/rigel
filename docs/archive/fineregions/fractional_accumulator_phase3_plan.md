# Phase 3 Fractional Accumulator Cutover Plan

Date: 2026-05-22
Status: planning document, implementation not started
Parent design: fine_region_migration_design_v3.md

## 1. Purpose

Phase 3 replaces the native integer calibration accumulator with a fractional
float32 accumulator. This is a hard cutover, not a compatibility shim. The old
8-mask integer payload, orientation subchannels, and binary boundary counters go
away.

The purpose of this plan is narrower than the future gDNA density redesign:

- accumulate fractional per-region evidence correctly and efficiently;
- export raw evidence structures that are rich enough for downstream design;
- identify every current interface that will break;
- define minimal, auditable current-compatible views where they are useful;
- avoid silently reconstructing old behavior where the new evidence should force
  a design decision.

## 2. Decisions And Corrections

### 2.1 Current-compatible contained classes are unspliced-only

The current tool's coarse contained density paths do not use spliced fragments.
For any temporary current-compatible coarse view, use only the unspliced channel:

- INTERGENIC: `contained_unspliced_pos + contained_unspliced_neg` over
  `signature == 0x0`.
- INTRON: `contained_unspliced_pos + contained_unspliced_neg` over intron-only
  signatures `{0x4, 0x8, 0xC}`.
- EXON-CONTAINED: `contained_unspliced_pos + contained_unspliced_neg` over
  signatures with any exon bit.

Spliced channels remain in the payload because they are useful future evidence,
but they must not enter the current-compatible contained density numerators.

### 2.2 Boundary evidence is side mass, not an old integer event

The accumulator should emit raw fractional side mass for every region side. It
should not decide the final EXON_INTRON gDNA estimator. The old integer
`u_left/u_right` event counted a fragment if it crossed an exon-region endpoint
with enough clearance. The new accumulator records how much aligned mass landed
on the left or right side of each region under the exact routing rule.

A temporary current-compatible boundary view may use only unspliced mass on the
intronic side of EXON_INTRON boundaries. That view is intentionally not a promise
that the future estimator will use exactly that numerator.

### 2.3 Accumulator-side K is removed

`BamScanConfig.splicing_anchor_tolerance` stays public because the resolver uses
it for implicit-splice gap classification. Phase 3 removes K from:

- `CalibrationAccumulator` constructor and fields;
- `BamScanner.set_regions(...)` accumulator plumbing;
- boundary numerator routing;
- boundary denominator helpers.

The scan payload may retain `resolver_splicing_anchor_tolerance` as provenance.
It must be documented as resolver provenance only.

### 2.4 FL evidence is collected as six named pools, not a joint table

FL pools are organized by `(region_coarse_class, compartment)`. The receiving
region's own signature determines the coarse class via the §2.2 derivation
table, so no `boundary_kind` lookup is needed at accumulation time. The six
pools are:

- `INTERGENIC_CONTAINED`
- `INTERGENIC_BOUNDARY`
- `INTRONIC_CONTAINED`
- `INTRONIC_BOUNDARY`
- `EXONIC_CONTAINED`
- `EXONIC_BOUNDARY`

Each pool is a single `float64[1024]` histogram. Total accumulator FL state is
6 × 1024 = 48 KB plus per-worker overhead. Compare to the original v3 proposal
of a `[16, 12, 1024]` joint table (≈1.6 MB) or a boundary-kind joint table
(≈9.4 MB).

Aggregating into a single gDNA FL pool is one line of downstream code:
`INTERGENIC_CONTAINED + INTERGENIC_BOUNDARY + INTRONIC_CONTAINED +
INTRONIC_BOUNDARY`. Keeping them separate during accumulation lets calibration
diagnostics compare pool shapes (e.g. intergenic vs intronic contamination) and
lets future estimators reweight or drop a pool without re-running the scan.

The EXON pools are recorded but not used for gDNA FL in Phase 3. They support
future expressed-vs-unexpressed region partitioning and diagnostic comparison
against the scanner's existing SPLICED_ANNOT RNA FL distribution.

### 2.5 Only unspliced fragments contribute to FL pools

For unspliced fragments, `frag_end - frag_start` is the genomic fragment length.
For spliced fragments, that span includes introns and is not the fragment length.
Spliced fragments must not contribute to FL pools at all under this design.

FL pools are strand-collapsed because gDNA is unstranded. Strand information is
preserved in `region_counts` for other consumers but is summed before the FL
bin is incremented.

`observe(...)` receives a precomputed `fl_idx` (the unspliced fragment length
bin, or `-1`). The accumulator:

- routes the fragment's mass into `region_counts` regardless of splice/strand;
- updates FL pools only if `splice_idx == 0` and `fl_idx >= 0`;
- counts unspliced fragments without a valid FL bin in `n_fl_unavailable` for
  audit; this is expected to be rare.

Boundary_kind columns in the v4 region table remain as future-facing metadata.
The Phase 3 accumulator does not read them.

## 3. Payload Contract

### 3.1 FL pool enum

Define a 6-value FL pool index used by both native and Python code:

```text
FL_POOL_INTERGENIC_CONTAINED = 0
FL_POOL_INTERGENIC_BOUNDARY  = 1
FL_POOL_INTRONIC_CONTAINED   = 2
FL_POOL_INTRONIC_BOUNDARY    = 3
FL_POOL_EXONIC_CONTAINED     = 4
FL_POOL_EXONIC_BOUNDARY      = 5
N_FL_POOLS                   = 6
```

A helper `fl_pool_index(signature, compartment) -> int` is implemented once in
`signature.py` and mirrored in `region_signature.h`. It maps:

- `coarse_type_from_signature(sig)` → `INTERGENIC | INTRON | EXON`;
- `compartment` → `CONTAINED | BOUNDARY_LEFT | BOUNDARY_RIGHT`;
- with left/right collapsed to a single `BOUNDARY` slot.

### 3.2 Native payload

```cpp
struct CalibrationPayload {
    static constexpr int32_t kFlBins  = 1024;
    static constexpr int32_t kChannels = 12;
    static constexpr int32_t kSignatures = 16;
    static constexpr int32_t kFlPools = 6;

    // Dense output built once after worker merge.
    std::vector<float> region_counts;            // float32[R, 12]
    int64_t n_regions = 0;

    // Global mass summaries over the whole genome.
    std::array<double, 12> channel_mass{};       // sum over regions and FL
    std::array<double, 16> signature_mass{};     // sum over channels and FL

    // Six named FL pools. Unspliced, strand-collapsed mass only.
    std::array<std::array<double, 1024>, 6> fl_pool_mass{};
    std::array<double, 6> fl_pool_total{};       // matches fl_pool_mass.sum(axis=1)

    // Accountability.
    int64_t n_observed_for_calibration = 0;
    int64_t n_excluded_multimapper    = 0;
    int64_t n_excluded_chimera        = 0;
    int64_t n_excluded_splice_artifact = 0;
    int64_t n_excluded_strand_ambig   = 0;
    int64_t n_unobserved              = 0;
    int64_t n_unannotated_ref         = 0;
    int64_t n_fl_unavailable          = 0;   // unspliced fragments with fl_idx < 0

    int32_t resolver_splicing_anchor_tolerance = 0;  // resolver provenance only
};
```

Total FL state per worker is 48 KB plus 96 B for pool totals. The 16-element
`signature_mass` and 12-element `channel_mass` summaries remain because they
are cheap and useful for sanity checks against `region_counts`. There is no
joint FL table.

### 3.3 Python payload

`CalibrationScanPayload` becomes:

```python
@dataclass(frozen=True, slots=True)
class CalibrationScanPayload:
    region_counts: np.ndarray              # float32[R, 12]
    fl_pool_mass: np.ndarray               # float64[6, 1024]
    fl_pool_total: np.ndarray              # float64[6]
    channel_mass: np.ndarray               # float64[12]
    signature_mass: np.ndarray             # float64[16]

    n_observed_for_calibration: int
    n_excluded_multimapper: int
    n_excluded_chimera: int
    n_excluded_splice_artifact: int
    n_excluded_strand_ambig: int
    n_unobserved: int
    n_unannotated_ref: int
    n_fl_unavailable: int
    resolver_splicing_anchor_tolerance: int
```

Convenience views on the payload (or in `fractional_evidence.py`):

- `channel(compartment, splice_idx, strand_idx) -> np.ndarray`
- `contained_unspliced -> np.ndarray`
- `boundary_left_unspliced -> np.ndarray`
- `boundary_right_unspliced -> np.ndarray`
- `region_total_mass -> np.ndarray`
- `fl_pool(name) -> np.ndarray`  (`name` in the 6-value enum)
- `gdna_fl_mass() -> np.ndarray` (sums the four INTERGENIC + INTRONIC pools)
- `rna_candidate_fl_mass() -> np.ndarray` (sums the two EXONIC pools)

Validation rules:

- `region_counts` is C-contiguous float32 `(R, 12)`, finite, non-negative.
- `fl_pool_mass` is float64 `(6, 1024)`, finite, non-negative.
- `region_counts.sum()` is approximately `n_observed_for_calibration`.
- `channel_mass.sum()` and `signature_mass.sum()` match observed mass.
- `fl_pool_mass.sum()` equals `n_observed_unspliced_with_fl`, the total mass of
  unspliced fragments with a valid FL bin. This is strictly less than
  `n_observed_for_calibration`.
- Balance identity:
  `observed + excluded_multimapper + excluded_chimera + excluded_splice_artifact
  + excluded_strand_ambig + n_unobserved + n_unannotated_ref == n_read_names`.

## 4. Native Accumulation Algorithm

### 4.1 Scanner gates

At the single calibration observation site in `bam_scanner.cpp`:

1. Multimappers call `note_multimap()` and do not observe.
2. Pure chimeras call `note_chimera()`.
3. Splice artifacts call `note_splice_artifact()`.
4. Resolved fragments whose selected evidence strand is not POS or NEG call
   `note_strand_ambig()`.
5. Resolved fragments with no usable aligned blocks call `note_unobserved()`.
6. Fragments on references absent from the region index call
   `note_unannotated_ref()`.
7. Everything else calls `observe(...)` and must emit total mass 1.0.

This keeps the one-fragment-one-counter identity explicit.

### 4.2 Observe signature

`CalibrationAccumulator::observe(...)` should receive exactly the data it needs:

```cpp
void observe(int8_t splice_type,
             int32_t ref_id,
             int64_t frag_start,
             int64_t frag_end,
             int32_t fl_idx,
             const ExonBlock* blocks,
             int32_t n_blocks,
             int8_t fragment_strand,
             const RegionIndex& regions);
```

`frag_start` and `frag_end` remain useful for diagnostics, but mass routing uses
`blocks` and region overlaps. `fl_idx` is a trusted fragment-length bin or `-1`.

### 4.3 Mass routing

For each observed fragment:

1. Compute `total_aligned_bp = sum(block.end - block.start)`.
2. If `total_aligned_bp <= 0`, call `note_unobserved()` instead of observing.
3. Compute `splice_idx = (splice_type == SPLICE_UNSPLICED ? 0 : 1)`.
4. Compute `strand_idx = (fragment_strand == STRAND_POS ? 0 : 1)`.
5. Determine `fl_eligible = (splice_idx == 0 && fl_idx >= 0)`. If `splice_idx
   == 0 && fl_idx < 0`, increment `n_fl_unavailable` once.
6. For every block and every overlapping region `(rid, rs, re)`:
   - `overlap_bp = min(block.end, re) - max(block.start, rs)`
   - `w = overlap_bp / total_aligned_bp`
   - decide `compartment` from `cross_left`/`cross_right`;
     fully-spanned regions split `w` evenly between `BOUNDARY_LEFT` and
     `BOUNDARY_RIGHT`;
   - add the resulting mass to sparse `region_counts[rid, channel]`;
   - add the same mass to global `channel_mass` and `signature_mass`;
   - if `fl_eligible`, add the same mass to
     `fl_pool_mass[fl_pool_index(signature[rid], compartment), fl_idx]`.

For the fully-spanned case, both halves of `w` route to the same FL pool
(`coarse_class × BOUNDARY`), so the pool receives mass `w` total.

### 4.4 FL pool routing example

User-described case: a 100 bp unspliced fragment overlaps an exon region (right
side, 98 bp) and an intron region (left side, 2 bp).

- `total_aligned_bp = 100`
- For the exon region: `overlap_bp = 98`, `cross_right = true` →
  `region_counts[exon_rid, BOUNDARY_RIGHT(splice=0,strand)] += 0.98` and
  `fl_pool_mass[EXONIC_BOUNDARY, fl_idx] += 0.98`.
- For the intron region: `overlap_bp = 2`, `cross_left = true` →
  `region_counts[intron_rid, BOUNDARY_LEFT(splice=0,strand)] += 0.02` and
  `fl_pool_mass[INTRONIC_BOUNDARY, fl_idx] += 0.02`.

The INTRONIC_BOUNDARY pool received 0.02 mass at this fragment's FL bin, as
specified. Aggregating gDNA FL later sums INTERGENIC_CONTAINED +
INTERGENIC_BOUNDARY + INTRONIC_CONTAINED + INTRONIC_BOUNDARY. The EXON
overhang contributes only to EXONIC_BOUNDARY, which is excluded from the gDNA
FL aggregate in Phase 3.

A spliced fragment with the same geometry routes the same mass into
`region_counts` spliced channels, but contributes nothing to FL pools because
spliced fragments have no genomic-span FL.

## 5. Sparse Worker Architecture

Use a private sparse map per worker and one deterministic merge.

```cpp
struct RegionSlot {
    int32_t rid;
    std::array<float, 12> count;
};

class SparseRegionCounts {
public:
    std::array<float, 12>& row(int32_t rid);
    void add(int32_t rid, int channel, float mass);
    void merge_from(const SparseRegionCounts& other);
    std::vector<float> to_dense(int64_t n_regions) const;
};
```

Implementation notes:

- Use a simple power-of-two open-addressing table with `-1` empty keys and slot
  indices as values. This avoids the allocation pattern of `std::unordered_map`.
- Store touched rows as AoS `std::array<float, 12>` so one region row fits in one
  cache line.
- Reserve 65,536 buckets initially; grow at load factor around 0.7.
- Merge workers in worker index order into the merged accumulator.
- Build the dense `float32[R, 12]` only once during payload export.

This preserves the memory goal: worker memory scales with distinct touched
regions, not total genome region count.

## 6. Boundary Exposure Helper

Add a new helper in `_exposure.py`:

```python
def fractional_boundary_side_exposure(lengths_bp, gdna_fl) -> np.ndarray:
    ...
```

For one region side of a region with span `S` and fragment length `L`, the exact
per-side exposure is:

```text
sum over starts crossing that side of side_mass = min((L - 1) / 2, S / 2)
```

Therefore:

```text
E_side(S) = sum_L pmf[L] * min((L - 1) / 2, S / 2)
```

For large regions, this reduces to `(E[L] - 1) / 2`, which is half of the old
strict-crossing exposure at K=0. For small fine regions, it correctly accounts
for fully-spanned regions splitting mass 50/50 across both sides.

`boundary_crossing_exposure(K)` should not be used by fractional consumers.
It can be deleted or retained only for tests of historical behavior outside the
fractional path.

## 7. Current-Compatible Evidence Views

These helpers are not the final density estimator. They are an audit layer that
lets us inspect fractional evidence in coarse terms.

Proposed module: `rigel.calibration.fractional_evidence`.

Core masks derived from signature:

```python
is_intergenic = signature == 0x0
is_intron_only = signature in {0x4, 0x8, 0xC}
has_exon = (signature & (BIT_EXON_POS | BIT_EXON_NEG)) != 0
```

Contained views:

```python
intergenic_contained_unspliced = payload.contained_unspliced[is_intergenic]
intron_contained_unspliced = payload.contained_unspliced[is_intron_only]
exon_contained_unspliced = payload.contained_unspliced[has_exon]
```

Boundary side rows:

```python
@dataclass(frozen=True, slots=True)
class BoundarySideEvidence:
    region_id: np.ndarray
    side: np.ndarray                  # 0 left, 1 right
    signature: np.ndarray
    neighbor_signature: np.ndarray
    mass_unspliced_pos: np.ndarray
    mass_unspliced_neg: np.ndarray
    exposure: np.ndarray
```

`boundary_kind` is not carried on this view. The receiving region's `signature`
plus the `neighbor_signature` are sufficient to reconstruct any kind label
downstream. The Phase 3 accumulator and current-compatible views do not depend
on the v4 `boundary_kind_{left,right}` columns; they remain in `regions.feather`
as future-facing metadata.

A temporary EXON_INTRON current-compatible view is:

```python
is_intron_only(signature) and has_exon(neighbor_signature)
```

It uses intron-side unspliced boundary mass and
`fractional_boundary_side_exposure(region_length, gdna_fl)` as denominator.
Mixed-signature EXON_INTRON sides should be surfaced separately as diagnostics,
not silently folded into this view.

## 8. Interface Audit

| Surface | Current contract | Phase 3 action |
| --- | --- | --- |
| `CalibrationAccumulator` native | Dense integer arrays per worker plus K tolerance | Replace with sparse fractional region counts, FL marginals, and explicit counters. Remove K. |
| `BamScanner.set_regions(...)` | Takes type masks, strands, K | Keep v4 metadata arrays; remove K argument. K remains on resolver setup only. |
| `BamScanner.scan()` calibration dict | Exports `global_counts`, `per_region_counts`, `fl_hist`, `u_left/u_right`, orientation arrays | Export fractional payload schema. No legacy keys. |
| `scan_and_buffer(...)` | Validates old payload with `CalibrationScanPayload.from_scan_dict` | Validate fractional schema and pass `stats.n_read_names` for balance. |
| `CalibrationScanPayload` | 8-mask integer schema | Replace with `region_counts[R,12]`, six named FL pools (`fl_pool_mass[6,1024]`), channel/signature mass summaries, and fractional counters. |
| `_fl_sources.extract_gdna_counts` | Reads old `fl_hist` rows `{INTRON, EXON\|INTRON, INTERGENIC}` | Read `payload.gdna_fl_mass()` (sums INTERGENIC_CONTAINED + INTERGENIC_BOUNDARY + INTRONIC_CONTAINED + INTRONIC_BOUNDARY). Spliced and exonic pools are excluded by construction. |
| `_diagnostics.Diagnostics` | Named decomposition of old 8 masks | Replace with fractional diagnostics: observed/excluded counters, mass summaries by signature/channel, and per-pool FL pool totals. |
| `density_global.compute_global_densities` | Builds INTERGENIC, INTRON, EXON-INTRON, EXON-CONTAINED from integer counts and `boundary_crossing_exposure(K)` | Rewrite around fractional evidence views. Contained classes use unspliced-only. Boundary path should initially be an explicit intron-side EXON_INTRON view (`is_intron_only(signature) and has_exon(neighbor_signature)`) or fail fast until the new estimator interface is accepted. |
| `RegionalGdnaExposure.build` | Uses `PayloadArrays` old integer columns and old boundary exposure | Rewrite after the new evidence interface is agreed. Do not bridge old arrays. |
| `_arrays.PayloadArrays` | Sorted integer count/flux columns | Replace with sorted fractional channel views. Drop `type`, `strand`, `bf_left`, `bf_right` dependency from hot code; keep `signature`, neighbor signatures. Accumulator does not consume `boundary_kind_{left,right}`; keep them as table-only metadata. |
| `locus_prior.estimate_locus_gdna` | Local diagnostic count model uses old contained and boundary arrays | Rewrite or gate behind the new evidence interface. The current function should not consume compatibility shims pretending to be old counts. |
| `expected_gdna_count_global` | Projects old global densities using old boundary exposure | Rewrite once new global density objects exist. The EM prior cannot be meaningfully green until this interface is designed. |
| `CalibrationResult` | Stores old diagnostics and K/n_below_tolerance | Store fractional diagnostics and resolver K provenance. Remove `n_below_tolerance`. |
| `regions.feather` `boundary_kind_{left,right}` | Phase 1 metadata columns | Retain in the table; not consumed by Phase 3 accumulator or FL pool routing. Available for future estimators. |
| `BamScanConfig` and CLI help | Describes K as calibration boundary tolerance | Update text: K controls implicit-splice resolver slack only. |
| Tests | Many factories construct old payloads | Replace factories with fractional payload builders. Add direct accumulator routing tests before broad rewrites. |
| Goldens | Phase 1 refreshed; Phase 2 unchanged | Phase 3 goldens will move materially after the density/prior interface is complete. |

## 9. Implementation Sequence

### 9.1 Planning checkpoint

1. Land this plan.
2. **Resolved**: FL evidence ships as six named pools
   (`INTERGENIC_{CONTAINED,BOUNDARY}`, `INTRONIC_{CONTAINED,BOUNDARY}`,
   `EXONIC_{CONTAINED,BOUNDARY}`). No joint boundary-kind FL table.
3. Decide whether Phase 3 implementation should include a temporary
   current-compatible density view or intentionally stop at a fractional evidence
   interface for the next design pass.

### 9.2 Native accumulator cutover

1. Replace `CalibrationPayload` in `accumulator.h`.
2. Implement `SparseRegionCounts`.
3. Rewrite `observe(...)` to route fractional mass directly per block.
4. Add `note_strand_ambig()`, `note_unannotated_ref()`, and `note_fl_unavailable()`.
5. Remove `splicing_anchor_tolerance_`, `qualified_hits_`, and
   `n_below_tolerance`.
6. Update worker construction and merge logic in `bam_scanner.cpp`.
7. Update native payload export.

### 9.3 Python payload cutover

1. Rewrite `scan_payload.py` for fractional schema validation.
2. Add channel convenience properties.
3. Rewrite `_fl_sources.extract_gdna_counts` to use the new conservative FL pool
   or the agreed evidence view.
4. Rewrite `_diagnostics.py` to report fractional mass summaries.

### 9.4 Interface layer

1. Add `fractional_evidence.py` with signature masks, contained views, boundary
   side views, and exposure helpers.
2. Make `compute_global_densities` consume that layer, or make it fail with a
   targeted error until the new density estimator design lands.
3. Rewrite `RegionalGdnaExposure.build` and `locus_prior` only after the evidence
   interface is accepted.

### 9.5 Documentation and CLI

1. Update `BamScanConfig.splicing_anchor_tolerance` docstring.
2. Update CLI help for `--splicing-anchor-tolerance`.
3. Expand CHANGELOG: calibration payload schema breaking change.
4. Update this running migration log with each implementation checkpoint.

## 10. Test Plan

Native/focused:

- Direct fractional accumulator helper or native test binding for worked examples.
- One-fragment-one-unit invariant for single-block and multi-block fragments.
- Fully-spanned region splits boundary mass 50/50.
- Splice axis routes to spliced channels but does not enter current-compatible
  contained classes.
- Strand axis routes POS/NEG channels.
- Strand-ambiguous fragment increments `n_excluded_strand_ambig` and emits no mass.
- `splicing_anchor_tolerance` changes resolver implicit-splice classification but
  does not change accumulator mass routing.
- `fl_idx < 0` emits region mass but not FL mass.

Python payload:

- Shape/dtype/finiteness validation.
- Region mass identity.
- Balance identity.
- FL-known mass accounting with `n_fl_unavailable`.
- Convenience channel properties match `channel_index(...)`.

Evidence views:

- Contained current-compatible views use unspliced only.
- EXON_INTRON side view selects intron-only side only and reports mixed-signature
  sides separately.
- `fractional_boundary_side_exposure` equals `(L - 1)/2` for large regions and
  `S/2` for short regions under degenerate FL.

Integration:

- Native scan returns fractional payload on existing synthetic scenarios.
- Multithread worker merge equals single worker within float32 tolerance.
- Old payload keys are absent.
- Density/prior tests are rewritten or explicitly expected to fail until the new
  interface layer is implemented.

## 11. Open Questions For The Interface Design

1. Should Phase 3 ship a temporary current-compatible global density estimator,
   or should `calibrate(...)` fail after producing fractional evidence until the
   new estimator is designed?
2. **Resolved**: FL evidence is six named pools by
   `(coarse_class_of_region, contained|boundary)`. Receiving region's signature
   selects the pool. EXON pools are recorded but not used for gDNA FL in Phase 3.
3. Should mixed-signature EXON_INTRON boundaries be excluded from the first
   boundary density view, or assigned to a separate experimental class?
4. Should spliced mass with unavailable FL be excluded from all FL tables, or
   should the scanner compute a trustworthy fragment length for more spliced
   cases before Phase 4? (Phase 3 default: spliced is excluded from FL pools by
   construction.)
5. What is the first downstream object the new estimator should consume:
   raw `CalibrationScanPayload`, a `FractionalEvidenceView`, or a precomputed
   boundary-side table?

## 12. Downstream Consumers Of Calibration Evidence

The FL model is only one consumer. Phase 3 must keep the following surfaces in
mind so the new payload does not need a follow-up reshape. Each row lists the
consumer, what it currently reads, and what the Phase 3 evidence should provide.

| Consumer | Reads today | Needs from fractional evidence |
| --- | --- | --- |
| Library-wide gDNA FL model (`_fl_sources.extract_gdna_counts` → `_fl_mixture`) | `fl_hist` rows for INTRON, EXON\|INTRON, INTERGENIC | `payload.gdna_fl_mass()` = sum of the four non-exon pools. EXON pools available for diagnostics only. |
| Library-wide gDNA mixing weight `π_pool` (`_fl_empirical_bayes`) | Same FL aggregate, plus total observed counts | Same aggregate; total observed mass comes from `region_counts.sum()`. |
| Global gDNA density estimator (`density_global.compute_global_densities`) | INTERGENIC contained density, INTRON contained density, EXON-INTRON boundary numerator, EXON contained density; uses `boundary_crossing_exposure(K)` | Fractional contained densities by coarse class via `fractional_evidence`; intron-side boundary numerator via `BoundarySideEvidence`; denominators via `fractional_boundary_side_exposure`. |
| Regional gDNA exposure (`RegionalGdnaExposure.build`) | Sorted `PayloadArrays` integer count/flux columns | Sorted fractional channel views from the new `PayloadArrays`. No legacy columns. |
| Locus priors (`locus_prior.estimate_locus_gdna`, `expected_gdna_count_global`, `assemble_priors`) | Per-locus contained/boundary counts and global density projections | Per-locus fractional sums over `region_counts` rows in the locus, plus global density objects rebuilt from fractional evidence. |
| EM warm-start and locus construction | Indirect: consumes per-locus priors only | No direct calibration payload dependency; benefits only after the priors are rewritten. |
| Calibration diagnostics (`_diagnostics.Diagnostics`, `summary.json`) | 8-mask integer decomposition, `n_below_tolerance` | Per-pool FL totals, per-channel mass, per-signature mass, observed/excluded counters. `n_below_tolerance` is removed. |
| `CalibrationResult` surfaces | Old diagnostics, K provenance | Fractional diagnostics and resolver K provenance. |
| Tests and goldens | Old payload factories | Fractional payload factories; Phase 3 will move goldens materially. |

### 12.1 Interface decisions Phase 3 must lock in

1. **Canonical payload object.** `CalibrationScanPayload` stays the
   serialization-shaped object that crosses the native boundary. A separate
   `FractionalEvidenceView` (in `fractional_evidence.py`) wraps it and exposes
   masks, contained views, boundary side rows, and FL pool helpers. Downstream
   consumers depend on the view, not on raw payload arrays.
2. **gDNA FL aggregate is a method, not a stored array.** Native does not export
   a `gdna_candidate_fl_mass`. `payload.gdna_fl_mass()` sums the four non-exon
   pools on demand so the policy is visible in Python and easy to revise.
3. **Per-pool totals are first-class.** `fl_pool_total[6]` is exported so
   diagnostics can detect pool imbalance without re-summing the full 6×1024
   array.
4. **Boundary side evidence is built from `region_counts` + region metadata.**
   It is not a separate native export. `BoundarySideEvidence` is a Python view
   over the boundary_left and boundary_right channels of `region_counts`,
   joined to `signature` and `neighbor_signature` from the v4 region table.
5. **Density estimator surface is decided in Phase 4, not Phase 3.** Phase 3
   either fails fast inside `compute_global_densities` or ships a clearly
   labelled `_legacy_compatible_density_view` that uses only the unspliced
   contained and intron-side boundary mass. Either is acceptable; silently
   reconstructing the old estimator from fractional inputs is not.
6. **`boundary_kind_{left,right}` is index-only metadata.** No native code path
   reads it in Phase 3. If a future density estimator needs it, it joins
   against `regions.feather` in Python.
