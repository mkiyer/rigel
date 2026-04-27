# Plan: SRD v2 — Genomic-Footprint Fragment Lengths and Strand-Aware Per-Fragment Overlap

## Context

SRD v1 is mostly implemented (commit `7ebd391` on `main`) but has clear, internally-visible bugs in the calibration pipeline. Note on framing: real RNA-seq libraries — even the "purest" — always have some gDNA contamination; there is no real-world ground truth of 0% gDNA. The reported `gdna_fraction = 0.9%` for `mctp_vcap_rna20m_dna00m` is itself a perfectly acceptable result and we are **not** necesarily trying to drive it to zero. The bugs below are about *internal correctness of the gDNA_FL fit*, not about the headline contamination number.

- Fitted `gDNA_FL` for `mctp_vcap_rna20m_dna00m` has its **mode at bin 0 bp**, which is physically impossible regardless of the underlying gDNA fraction.
- 56% of fitted `gDNA_FL` mass is in bin 0; 4% saturated at bin 1000. Both are sentinel-value / clipping artifacts, not physical signal.
- Across the spike series, calibration **underestimates gDNA at high contamination** — separate from the bin-0 artifact, but plausibly downstream of the same broken `gDNA_FL` shape (a `gDNA_FL` whose mode is 0 has near-zero likelihood at the real fragment lengths gDNA actually occupies, so per-fragment posterior γ_f is suppressed at gDNA-contaminated loci).

The bugs are clear regardless of whether the headline `gdna_fraction` is "right" or "wrong" — a fitted FL distribution with mode = 0 cannot be correct. Fixing them improves the *shape* of `gDNA_FL`, which we expect to improve high-contamination accuracy materially while leaving the low-contamination headline number roughly where it is (good).

### Three confirmed root causes

**Cause 1 — `-1` sentinels masquerading as length 0.** `resolve_context.h::compute_frag_lengths` (lines 725–757) only inserts entries when the *transcript-projected* fragment length is `> 0`. For an `EXON_INCOMPATIBLE` fragment that overhangs an exon edge into intron/intergenic territory on candidate `T_j`, the projected `tx_e − tx_s` collapses to ≤ 0, so transcript `T_j` is omitted from the sparse map. Then `ResolvedFragment::from_core` (lines 165–169) fills the parallel `frag_lengths[i]` array with `-1` for missing entries. `_categorize.py` reads the **first candidate's** length (line 144) — which can be `-1`. `_simple.py` then does `np.clip(pool_lens, 0, max_size)` (line 134), silently converting `-1 → 0`. Bin 0 fills up; the 1-D mixture parks all of it in `gDNA_FL` because the fixed RNA_FL has near-zero likelihood there.

**Cause 2 — Wrong primitive.** Using *transcript-projected* fragment length to characterize *unspliced* fragments is semantically wrong. The transcript projection answers "if this fragment were mature mRNA from transcript T, what would its spliced length be?" — but for unspliced gDNA candidates that's a nonsense question. The right primitive is the **contiguous genomic span**: `genomic_end − genomic_start`. This already exists in the buffer as `genomic_footprint` (well-defined, never `-1` for resolved fragments) but is not used for the calibration pool length.

**Cause 3 — Categorization signals are too weak.** SRD v1 categorizes only by `max(exon_bp)` across candidates and `read_length`. It cannot distinguish "fragment falls inside an exon of one transcript while the same position is intronic on a second transcript" from "fragment overhangs a single transcript's exon by 30 bp into intergenic". It has no notion of **strand-specific overlap** (the partition's old `tx_pos / tx_neg / exon_pos / exon_neg` flags are gone). And `INTERGENIC` and `INTRONIC` Python counts are essentially always 0: the C++ resolver drops zero-candidate fragments before the buffer (resolve_context.h:1043), and INTRONIC requires "has candidates but zero exon overlap on every candidate", which is rare in practice. The pool is therefore almost entirely `EXON_INCOMPATIBLE` — exactly the bucket where Cause 1 lives.

The intergenic FL accumulator (bam_scanner.cpp:1320–1326) catches truly intergenic fragments by genomic footprint and feeds them in via `frag_length_models.intergenic`, but the user observes its FL shape is RNA-shaped (mean 236, mode 223) — suggesting these are dominated by unannotated transcription / off-target reads, not gDNA. So there's also a real biological issue: in a pure-RNA library there isn't much true gDNA signal, and any signal we *do* extract is mostly artifact.

**Beyond bug fixes**, the user's feedback also identifies missing pipeline capabilities:

- **Implicit splicing** (paired-end unsequenced gap contains an annotated intron) is not detected; such fragments are currently classified as `UNSPLICED` and can leak into the gDNA pool.
- **Splice artifacts** (CIGAR-N matches that pass the SJ blacklist) are currently downgraded to `UNSPLICED` and their disjoint aligned blocks are merged into a single "footprint" — which both inflates fragment length and treats artifact RNA as a gDNA candidate. They should be held out from calibration.

## Design goals (unchanged but reaffirmed)

- Library-agnostic, blazing fast, geometric-first-principles.
- No cross-region density; no SS-threshold magic numbers.
- Empirical-Bayes global-FL prior; graceful QC-fail degradation.
- Calibration runs on **uniquely-aligned, true-unspliced** fragments only.
- Per-locus γ prior preserved (already in main).
- One geometric source of truth: the BAM scanner; calibration is a pure consumer.

## Design: SRD v2

### Pillar 1 — Use genomic footprint as the calibration fragment length

For unspliced fragments, the contiguous genomic span (`genomic_footprint = genomic_end − genomic_start`) is the correct primitive. It is:

- **Always defined and positive** for any resolved fragment.
- **Physically meaningful** — what would have been measured on a gel.
- **Independent of transcript projection** — works identically for INTERGENIC, INTRONIC, EXON_INCOMPATIBLE, and EXON_CONTAINED fragments.
- **Already in the buffer** as `genomic_footprint` (buffer.py:142, populated in resolve_context.h:155).

The existing per-candidate `frag_lengths` array stays — it's still used by the per-locus EM scorer for transcript-relative spliced-coordinate likelihood. Calibration just stops reading it.

### Pillar 2 — Strand-aware per-fragment overlap counts in the scanner

Add four collapsed-across-transcripts, separated-by-strand fields to `RawResolveResult` and the buffer, populated inline in the existing `_resolve_core` cgranges hit walk:

```cpp
int32_t exon_bp_pos;     // bp of fragment in any (+) transcript's exons
int32_t exon_bp_neg;     // bp of fragment in any (−) transcript's exons
int32_t intron_bp_pos;   // bp of fragment in any (+) transcript's introns
int32_t intron_bp_neg;   // bp of fragment in any (−) transcript's introns
```

Computation: the existing inner loop (resolve_context.h ~lines 868–945) already iterates each cgranges hit per fragment exon-block, knows the hit's `interval_type` (EXON / TRANSCRIPT / INTERGENIC) and the `t_set`, and computes `overlap_bp`. We add: peek into `t_strand_arr_[t]` for at least one `t` in `t_set` of each strand, and route the `overlap_bp` into `exon_bp_pos` / `exon_bp_neg` / `intron_bp_pos` / `intron_bp_neg` accordingly. Cgranges intervals in this index are non-overlapping by construction (the tiling is collapsed exon/intron/intergenic), so simple addition is correct — no per-base bitmap needed.

Derived quantities (computed at categorization time, not stored):

```python
tx_bp_pos       = exon_bp_pos + intron_bp_pos
tx_bp_neg       = exon_bp_neg + intron_bp_neg
exon_bp_either  = bp covered by EXON-type cgranges intervals (any strand)
                = exon_bp_pos + exon_bp_neg − exon_bp_both     # see below
intergenic_bp   = genomic_footprint − tx_bp_either
```

Note on **strand-overlap zones**: a single cgranges EXON interval can have *both* + and − transcripts in its `t_set` when an exon position is covered by transcripts on both strands (real and common at compact mammalian loci). For simplicity we accumulate `exon_bp_pos += bp` AND `exon_bp_neg += bp` in those cases; this means `exon_bp_pos + exon_bp_neg` can exceed `genomic_footprint`, and we treat such positions as **strand-ambiguous** (see Pillar 3). We do not need to explicitly compute `exon_bp_both`.

Buffer cost: 4 × int32 = 16 bytes per fragment. For 100M fragments, ~1.6 GB. Comparable to or smaller than the previously-removed `region_evidence` payload.

### Pillar 3 — Categorization on contiguous genomic footprint, with strand-aware sub-labels

Pass 0 inputs (per fragment, all already in or about to be in the buffer):

| Field | Meaning | Source |
|---|---|---|
| `splice_type` | SPLICED_ANNOT / SPLICED_UNANNOT / UNSPLICED / SPLICE_ARTIFACT (new) | scanner |
| `num_hits` | NH tag | scanner |
| `genomic_footprint` | end − start of contiguous genomic span | scanner (existing) |
| `exon_bp_pos`, `exon_bp_neg` | strand-collapsed exon overlap bp | **scanner (new)** |
| `intron_bp_pos`, `intron_bp_neg` | strand-collapsed intron overlap bp | **scanner (new)** |

Filter: `num_hits == 1` AND `splice_type == UNSPLICED`. (Splice artifacts are held out.)

Categories (closed form, vectorized numpy):

```
Let:
  fp = genomic_footprint
  exon_either_bp_lower_bound = max(exon_bp_pos, exon_bp_neg)
  exon_either_bp_upper_bound = exon_bp_pos + exon_bp_neg
  tx_either_bp_lower_bound  = max(tx_bp_pos, tx_bp_neg)
  tol = exon_fit_tolerance_bp  (default 3)

INTERGENIC      : tx_bp_pos == 0  AND  tx_bp_neg == 0
INTRONIC        : exon_bp_pos == 0  AND  exon_bp_neg == 0  AND  not INTERGENIC
EXON_CONTAINED  : (fp − exon_either_bp_lower_bound) ≤ tol     # entirely covered by exon(s) of at least one strand
EXON_INCOMPATIBLE: otherwise (some exon overlap, but overhang > tol)
```

Strand sub-label (for `EXON_CONTAINED` and `EXON_INCOMPATIBLE`):

```
strand_pos     : exon_bp_pos > 0  AND  exon_bp_neg == 0
strand_neg     : exon_bp_neg > 0  AND  exon_bp_pos == 0
strand_ambig   : exon_bp_pos > 0  AND  exon_bp_neg > 0       # antisense-overlap zone
```

Diagnostic counts emitted in `CalibrationResult`: per-category × per-strand-sublabel counts, plus separate counts for `INTRONIC_pos / INTRONIC_neg / INTRONIC_ambig` for symmetry.

### Pillar 4 — Pool assembly and gDNA_FL fit

Pool for the 1-D mixture EM:

```
Pool = INTERGENIC ∪ INTRONIC ∪ EXON_INCOMPATIBLE
```

Pool fragment length: **`genomic_footprint`**, full stop. No `np.clip` saturation tricks; the histogram bin range is `[0, max_size]` and out-of-range values are dropped (with a diagnostic count emitted, not silently clipped). For typical libraries `genomic_footprint < 1000` for ≥ 99.9% of fragments; saturation is genuinely rare and worth surfacing rather than hiding.

`UNSPLICED_ANTISENSE_EXONIC` (now a strand sub-label of `EXON_CONTAINED`) remains **excluded from the pool** — at highly-expressed loci it is RNA-dominated unless SS ≈ 1.0 (math derivation in earlier plan iterations). It stays as a diagnostic.

Pass 1 (mixture): `pool_FL(l) = π · gDNA_FL(l) + (1−π) · RNA_FL(l)`, RNA_FL fixed (from SPLICED), `gDNA_FL` and `π` free. Empirical-Bayes Dirichlet smoothing toward global FL, ESS = 500. Quality enum (`good` / `weak` / `fallback`) propagated.

Pass 2 (per-locus γ aggregation): unchanged from current main — per-fragment posterior γ_f from the recovered FL models, summed per locus to produce `α_gDNA[ℓ] = γ_ℓ · c_base`.

### Pillar 5 — Scanner: correctly classify spliced fragments

**EXPLICIT splicing (existing).** CIGAR `N` matching an annotated SJ that passes the blacklist filter → `SPLICED_ANNOT`. CIGAR `N` matching an unannotated SJ → `SPLICED_UNANNOT`. Both excluded from the gDNA pool.

**IMPLICIT splicing (NEW).** A fragment with no CIGAR `N` but whose paired-end unsequenced gap contains an annotated intron should be reclassified. Detection: when R1's end and R2's start are non-overlapping (gap between aligned blocks > 0) and the gap genomic interval is fully spanned by an annotated intron of any transcript whose exons R1 and R2 each overlap, mark `splice_type = SPLICED_IMPLICIT`. This requires a new SJ-spans-gap query against an annotated-introns index (already exists in the resolver context). Treated identically to `SPLICED_ANNOT` for calibration: excluded from the pool and excluded from `genomic_footprint` reporting (because the genomic span includes the intron, which is an artifact of the projection, not of the fragment).

**SPLICE_ARTIFACT (NEW).** A fragment whose CIGAR-`N` match was rejected by the blacklist filter currently downgrades to `UNSPLICED` and the disjoint aligned blocks are merged into a single misleading "footprint". Per the user, for v1: classify as `SPLICE_ARTIFACT` (a new splice_type value), do **not** merge the blocks, and **hold out** these fragments from calibration entirely. The per-locus EM can still see them as multimapper-like ambiguity.

### Pillar 6 — Region partition: explicitly stays gone

The cgranges interval index already provides everything we need at scan time. We do **not** resurrect `region_df`, `region_cr`, or `RegionAccumulator`. The "regional" semantics live entirely in the per-fragment strand-aware overlap counts.

### Why this fixes the observed bugs

| Bug | Root cause | Fix |
|---|---|---|
| `gDNA_FL` mode = 0 (physically impossible) | `-1` sentinels in per-candidate `frag_lengths` clipped to 0 | Calibration uses `genomic_footprint`, never the per-candidate array |
| 4% mass at L=1000 | `np.clip(..., max_size)` saturating | No clip; out-of-range fragments dropped with diagnostic count |
| Pool dominated by `EXON_INCOMPATIBLE` artifacts | `INTERGENIC` / `INTRONIC` Python counts always 0 because the C++ resolver drops zero-candidate fragments | Scanner exposes strand-aware overlap counts on **all** unique-mapped unspliced fragments; the resolver no longer needs to drop intergenic ones (their strand-overlap counts are simply all zero) |
| Splice artifacts polluting gDNA pool | Blacklisted CIGAR-N → `UNSPLICED` with merged footprint | New `SPLICE_ARTIFACT` class held out from calibration |
| Implicit-splice RNA in gDNA pool | Paired-end gap with intron not detected | New `SPLICED_IMPLICIT` class at scan time |
| High-contamination `gdna_fraction` underestimated | Plausibly: a `gDNA_FL` with mode = 0 has near-zero likelihood at real gDNA fragment lengths, suppressing per-fragment γ_f at gDNA-rich loci | Fixed `gDNA_FL` shape → expected to improve; verify at Phase 6 |

---

# Phased execution

## Phase 0 — Audit (1 hour, no code changes)

Verify and document:

1. The exact line(s) where `compute_frag_lengths` returns the empty/sparse map; quote the conditions producing missing entries.
2. The `from_core` lines that fill `frag_lengths[i] = -1`.
3. The `np.clip` line in `_simple.py`.
4. Confirm `genomic_footprint` is populated for every resolved fragment (including those with empty `t_inds` if any survive — currently none do, but post-fix they will).
5. Locate the resolver's "drop zero-candidate fragments" point (resolve_context.h:1043 per the audit) and confirm it's the same path that prevents intergenic Python categorization.
6. Locate the SJ blacklist filter, the splice_type assignments, and the merge-block-into-footprint path.
7. Locate any existing annotated-introns index in the resolver context (needed for implicit splicing detection).
8. Quantify the fraction of pool fragments with `frag_length = -1` in a small replay of an existing buffer, if feasible.

**Deliverable:** a short (`docs/calibration/srd_v2_audit.md`) listing file:line refs and a one-line decision for each.

## Phase 1 — Fast safety patch (calibration-only)

Switch SRD v1 to use `genomic_footprint` for the pool histogram **without** scanner changes. This is a 10–20 line patch and should immediately eliminate the 56%/4% artifact bins:

**File:** `src/rigel/calibration/_categorize.py`
- Replace `frag_length[valid_first] = cand_frag_lengths[idx]` (line 144) with `frag_length = chunk.genomic_footprint.astype(np.int32)`.
- Drop the no-candidates fallback to `read_length` (line 149) — `genomic_footprint` is universally populated.
- Keep `cand_frag_lengths` available only for downstream consumers that want it (none currently in calibration).

**File:** `src/rigel/calibration/_simple.py`
- Replace `np.clip(pool_lens, 0, max_size)` (line 134) with explicit out-of-range filtering and a diagnostic counter:
  ```python
  in_range = (pool_lens >= 0) & (pool_lens <= max_size)
  n_pool_dropped_out_of_range = int(pool_lens.size - in_range.sum())
  pool_hist = np.bincount(pool_lens[in_range].astype(np.intp), minlength=n_bins).astype(np.float64)
  ```
- Add `n_pool_dropped_out_of_range` to `CalibrationResult` diagnostics.

**Tests:** rerun `mctp_vcap_rna20m_dna00m`. Acceptance: `gDNA_FL` bin-0 fraction < 1% and bin-1000 fraction < 1%; the fitted `gDNA_FL` mode lies in a physically plausible range (e.g. 100–400 bp) rather than at 0. We are explicitly **not** asserting anything about the headline `gdna_fraction` — `mctp_vcap_rna20m_dna00m` has no ground truth and 0.5–1% is a perfectly acceptable result for a real "pure RNA" library.

This phase is independent of Phases 2+ and ships value immediately.

## Phase 2 — Scanner: strand-aware per-fragment overlap counts

Add the four new fields end-to-end:

**Files:**
- `src/rigel/native/resolve_context.h`: add `exon_bp_pos`, `exon_bp_neg`, `intron_bp_pos`, `intron_bp_neg` to `RawResolveResult`. Augment the inner loop (~lines 868–945) to accumulate them from the existing `t_set` + `t_strand_arr_` data. Validate invariant `exon_bp_pos + exon_bp_neg + intron_bp_pos + intron_bp_neg + intergenic_bp = genomic_footprint` (with the strand-overlap caveat: positions exonic on both strands count in both pos and neg).
- `src/rigel/native/bam_scanner.cpp`: thread the new fields into `FragmentAccumulator` and the buffer-fill path.
- `src/rigel/buffer.py`: add the four new columns to `BufferedFragment` / `_FinalizedChunk`. Update the chunk-from-raw conversion.
- Important: **stop dropping zero-candidate fragments** in `_resolve_core` (resolve_context.h:1043). They get empty `t_inds` but valid `genomic_footprint` and zeroed strand-aware overlap counts → categorize as `INTERGENIC`. The intergenic FL accumulator can be removed once this works.

**Tests:** `tests/test_resolve_strand_overlap.py` — synthetic small annotation, hand-constructed fragments crossing exon/intron/intergenic boundaries on +/− strands, asserting exact bp counts. Include the strand-overlap zone case.

## Phase 3 — Categorization rewrite

**Files:**
- `src/rigel/calibration/_categorize.py`: rewrite per the Pillar 3 rules. Stop reading `cand_frag_lengths`, `t_offsets`, `t_indices`, `intron_bp`, `exon_bp` — read only the new collapsed columns plus `splice_type`, `num_hits`, `genomic_footprint`. Add strand sub-labels. Vectorized numpy throughout.
- `src/rigel/calibration/_result.py`: extend `category_counts` with strand sub-labels.
- `src/rigel/calibration/_simple.py`: pool = `INTERGENIC ∪ INTRONIC ∪ EXON_INCOMPATIBLE`, length = `genomic_footprint`. Remove the `frag_length_models.intergenic` merge — those fragments now flow through the buffer.

**Tests:** `tests/test_categorize_v2.py` — six synthetic categories × three strand sub-labels = 18 cases, plus the antisense-overlap zone, plus tolerance edge cases (overhang of exactly `tol` and `tol + 1`).

## Phase 4 — Implicit splicing detection

**Files:**
- `src/rigel/native/resolve_context.h`: new helper `gap_spans_annotated_intron(r1_end, r2_start, t_inds)` querying the existing intron index. Returns true if the unsequenced gap is fully spanned by an annotated intron of a transcript that R1 and R2 each overlap exonically.
- `src/rigel/native/bam_scanner.cpp`: after CIGAR-derived splice classification, if `splice_type == UNSPLICED` and the fragment has a paired-end gap > 0, run the helper. On match, set `splice_type = SPLICED_IMPLICIT` and recompute `genomic_footprint` to exclude the intron.
- `src/rigel/types.py`: add `SpliceType.SPLICED_IMPLICIT = 3`.

**Tests:** `tests/test_implicit_splicing.py` — synthetic paired-end fragments with gaps containing introns, gaps not containing introns, gaps containing introns of the wrong transcript, etc.

## Phase 5 — Splice artifact handling

**Files:**
- `src/rigel/native/bam_scanner.cpp`: when blacklist filter rejects all CIGAR `N` junctions, set `splice_type = SPLICE_ARTIFACT` (new `= 4`) instead of downgrading to `UNSPLICED`. Do not merge the disjoint aligned blocks; preserve them per-block. Decision: emit one buffer fragment with `splice_type = SPLICE_ARTIFACT` and a flag indicating block disjunction; categorization filter excludes these from calibration entirely.
- `src/rigel/calibration/_categorize.py`: filter `splice_type == UNSPLICED` only (drops SPLICE_ARTIFACT and all SPLICED_* types).

**Tests:** `tests/test_splice_artifact.py` — synthetic blacklisted SJ scenarios.

## Phase 6 — Validation

Rerun the spike series. Acceptance criteria are framed around **shape and trend**, not absolute numbers — none of these libraries have ground-truth gDNA fractions.

- `mctp_vcap_rna20m_dna00m` ("pure" RNA, real-world floor):
  - `gDNA_FL` bin-0 + bin-1000 mass < 2% (artifact removal — the actual claim being tested).
  - `gDNA_FL` mode in a physically plausible range (100–400 bp).
  - Headline `gdna_fraction` *not* a pass/fail criterion. Expect it to stay near 0.5–1%; if it changes a lot in either direction, investigate why before declaring success.
- `mctp_vcap_rna20m_dna20m`, `mctp_vcap_rna20m_dna80m` (real spikes):
  - Headline `gdna_fraction` should *increase* monotonically with the spike level.
  - The current SRD v1 underestimates at high contamination; we expect SRD v2 to do better, but we don't pretend to know the right number absent ground truth. Track relative improvement vs SRD v1 baseline.
  - `gDNA_FL` mode physically plausible across all three libraries; ideally similar across libraries (the same fragmentation protocol should produce similar shapes).
- `gdna_fl_model_quality` should be `"good"` on the spiked libraries and may legitimately be `"weak"` or `"fallback"` on the unspiked one (where there's little real gDNA signal to fit).

Profile to confirm scanner overhead from new fields < 5%; total calibration time still well under v5.

**Deliverable:** `docs/calibration/srd_v2_results.md` with before/after tables.

## Phase 7 — Cleanup

- Remove the now-unused intergenic FL accumulator path in `bam_scanner.cpp` (zero-candidate fragments now flow through the buffer).
- Remove `cand_frag_lengths` consumption in calibration (keep the field for the per-locus EM scorer).
- Update `MANUAL.md`, `METHODS.md`, `CHANGELOG.md`.

---

## Risk register

| Risk | Mitigation |
|---|---|
| Phase 1 patch surfaces a different artifact (e.g. `genomic_footprint` itself has a bug) | Phase 0 audit explicitly verifies `genomic_footprint` is universally populated and physically correct |
| Phase 2 strand-overlap counts off-by-one at exon/intron boundaries | Phase 2 invariant test: `exon_bp_pos + exon_bp_neg + intron_bp_pos + intron_bp_neg ≥ genomic_footprint` (with strand-overlap caveat allowing `>`) |
| Phase 2 scanner overhead noticeable | Profile in Phase 2 exit; the extra work is constant-factor on the existing cgranges hit walk |
| Removing zero-candidate-drop changes locus construction | Verify locus.py is unaffected — INTERGENIC fragments have empty `t_inds` and won't enter any locus's per-transcript machinery |
| Implicit-splice detection misses cases / false-positives | Synthetic test suite covers the obvious cases; real-data trace at Phase 6 catches the rest |
| Splice-artifact "don't merge blocks" requires deeper plumbing changes | Phase 5 audit; if it gets ugly, defer block-preservation and just hold-out the merged-footprint fragment with the new splice_type |
| nRNA contamination still pollutes gDNA_FL even with clean pool | Acknowledged trade-off; nRNA_FL ≈ RNA_FL assumption documented; the 1-D mixture still recovers gDNA_FL when it differs from RNA_FL |
| Pure-RNA library "no real gDNA signal to fit" causing weak gDNA_FL even after Phase 1 | Expected; `gdna_fl_model_quality = "fallback"` correctly fires; verify it does |

## Out of scope (deliberate)

- Strand-balanced overlap cluster channel (deferred).
- Cross-region density signals.
- Mappability-corrected effective length.
- Strand-model / RNA-FL training changes.
- Re-introducing `region_df` / `region_cr` / `RegionAccumulator`.