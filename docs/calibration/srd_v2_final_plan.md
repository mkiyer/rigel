# SRD v2 — Final Plan

> Source review of `srd_v2_bugs_and_new_plan.md` plus code‑level verification.
> The diagnosis is correct; this revision tightens the implementation, fixes a
> subtle accounting issue in Pillar 2, defers one over‑ambitious phase, and
> sharpens the acceptance criteria.

## 0. Review summary of the prior plan

**What the prior plan got right.**
- Root‑cause #1 (`-1` sentinels → bin 0 via `np.clip`) is verified at
  [resolve_context.h:756](src/rigel/native/resolve_context.h#L756) (`if (fl > 0) result[t] = fl;`),
  [resolve_context.h:165-168](src/rigel/native/resolve_context.h#L165-L168) (`-1` fill),
  [_categorize.py:144](src/rigel/calibration/_categorize.py#L144),
  [_simple.py:130-132](src/rigel/calibration/_simple.py#L130-L132).
- Root‑cause #2 (transcript‑projected length is the wrong primitive for
  unspliced fragments) is correct. `genomic_footprint` is universally
  populated in the buffer ([buffer.py:142](src/rigel/buffer.py#L142),
  [resolve_context.h:155](src/rigel/native/resolve_context.h#L155)).
- Root‑cause #3 (categorization signals are too weak; `INTERGENIC`/`INTRONIC`
  Python counts are essentially always 0) is verified by inspection of the
  spike series: every library reports `category_counts[5] = 0` and
  `category_counts[6] = 0`. The intergenic accumulator path lives at
  [bam_scanner.cpp:227-238](src/rigel/native/bam_scanner.cpp#L227-L238).

**What the prior plan needs to change.**

1. **Pillar 2's intron accounting is misspecified.** The cgranges index has
   only `ITYPE_EXON` and `ITYPE_TRANSCRIPT` intervals (no explicit `ITYPE_INTRON`
   — see the inner loop at
   [resolve_context.h:880-915](src/rigel/native/resolve_context.h#L880-L915)).
   `intron_bp_*` cannot be accumulated directly; it must be **derived** as
   `max(tx_bp_strand − exon_bp_strand, 0)`. The prior plan's wording
   ("route the `overlap_bp` into `intron_bp_pos`") implies a non‑existent
   interval class. Fixed below in §2.
2. **Phase 4 (implicit splicing) is out of scope for v2.** It is high
   complexity (requires per‑pair gap tracking in the scanner plus a new
   intron‑span query), it modifies the same source files as Phases 2 and 5,
   and it is **not on the critical path** for the observed bug. Defer to
   a follow‑up; document the leak as an accepted limitation for v2.
3. **Consolidate scanner‑touching work into one phase.** Phases 2 and 5
   both edit `resolve_context.h` and `bam_scanner.cpp`. Splitting them
   doubles the number of recompiles and full‑sample re‑runs. Merge.
4. **Add an explicit nRNA‑pool‑pollution caveat.** The 1‑D mixture
   `pool = π·gDNA + (1−π)·RNA_FL` cannot separate nRNA from gDNA when
   their FL distributions overlap. nRNA fragments enter the pool
   (INTRONIC, EXON_INCOMPATIBLE) and inflate `π`. This is a real
   limitation, not a bug. Document and surface a diagnostic.
5. **Strengthen acceptance criteria with quantitative pre‑computed
   targets** (see §6) so we can call success/failure unambiguously.

## 1. Goals (unchanged)

- Library‑agnostic; geometric‑first‑principles.
- One source of geometric truth: the BAM scanner.
- No magic constants tied to strand specificity.
- Empirical‑Bayes shrinkage to global FL on small/weak pools.
- Calibration runs on **uniquely‑aligned, true‑unspliced** fragments only.
- Per‑locus γ prior preserved (already in main).

## 2. Design

### Pillar 1 — `genomic_footprint` is the calibration FL primitive

For every fragment in the calibration pool, the fragment‑length value used
in the 1‑D mixture is **`chunk.genomic_footprint`** (`int32`, contiguous
genomic span in bp; defined for every resolved fragment, including
zero‑candidate ones once Pillar 2b lands).

The per‑candidate `frag_lengths` array stays in the buffer for the per‑locus
EM scorer, which needs the spliced transcript‑relative length. Calibration
just stops reading it.

### Pillar 2a — Strand‑aware per‑fragment overlap counts (scanner)

Add four `int32` fields to `RawResolveResult` and `_FinalizedChunk`:

```cpp
int32_t exon_bp_pos;   // bp of fragment overlapping ANY (+)-strand transcript's exon
int32_t exon_bp_neg;   // bp of fragment overlapping ANY (−)-strand transcript's exon
int32_t tx_bp_pos;     // bp of fragment overlapping ANY (+)-strand transcript's span
int32_t tx_bp_neg;     // bp of fragment overlapping ANY (−)-strand transcript's span
```

Note: we store `tx_bp_*` (transcript‑span overlap), **not** `intron_bp_*`,
because the cgranges index has only `EXON` and `TRANSCRIPT` interval types.
`intron_bp_strand` is derived in Python as
`max(tx_bp_strand − exon_bp_strand, 0)`.

Inner‑loop edit at [resolve_context.h:880-915](src/rigel/native/resolve_context.h#L880-L915):

```cpp
// Existing per-cgranges-hit code already computes:
//   bp = max(0, min(bend, h_end) - max(bstart, h_start))
//   t_set = t_set_data_[off : off + cnt]
// Cgranges intervals are non-overlapping by construction (collapsed
// exon/transcript tiling), so simple addition is correct.
bool any_pos = false, any_neg = false;
for (int32_t k = 0; k < cnt; k++) {
    int8_t s = t_strand_arr_[t_set_data_[off + k]];
    any_pos |= (s == STRAND_POS);
    any_neg |= (s == STRAND_NEG);
}
if (itype == ITYPE_EXON) {
    if (any_pos) cr.exon_bp_pos += bp;
    if (any_neg) cr.exon_bp_neg += bp;
}
if (itype == ITYPE_TRANSCRIPT) {  // NB: TRANSCRIPT intervals already
    if (any_pos) cr.tx_bp_pos += bp;   //     subsume the exon footprint
    if (any_neg) cr.tx_bp_neg += bp;   //     of the same strand.
}
```

**Strand‑overlap zone semantics.** When a single cgranges interval has both
+ and − transcripts in its `t_set`, `bp` is added to both pos and neg.
Therefore `exon_bp_pos + exon_bp_neg ≥ exon_bp_either` always, and Pillar 3
treats positions covered by both strands as **strand‑ambiguous** (this is
the "antisense‑overlap zone" case, common at compact mammalian loci).

**Invariant for tests.** For every fragment:
- `exon_bp_pos ≤ tx_bp_pos`
- `exon_bp_neg ≤ tx_bp_neg`
- `max(tx_bp_pos, tx_bp_neg) ≤ genomic_footprint`
- (without strand‑overlap) `tx_bp_pos + tx_bp_neg ≤ genomic_footprint + min(tx_bp_pos, tx_bp_neg)`

**Buffer cost.** 4 × int32 = 16 bytes/fragment. ~1.6 GB at 100M fragments.

### Pillar 2b — Stop dropping zero‑candidate fragments in the resolver

The current resolver early‑exits when a fragment has no transcript candidates
(this is the path that prevents `INTERGENIC` from ever being Python‑categorized
and forces the C++ `intergenic_lengths` accumulator at
[bam_scanner.cpp:237-238](src/rigel/native/bam_scanner.cpp#L237-L238)).

Stop the early‑exit. A zero‑candidate fragment now flows through the buffer
with:
- empty `t_inds` (already handled downstream by every consumer that checks
  `n_cands > 0`)
- `genomic_footprint > 0`
- `exon_bp_pos = exon_bp_neg = tx_bp_pos = tx_bp_neg = 0`
- categorized as `INTERGENIC` in Pillar 3.

The C++ `intergenic_obs / intergenic_truth / intergenic_lengths` accumulators
become redundant and are removed in Pillar 5 cleanup.

**Verify** before merging: `locus.py` and the per‑locus EM never iterate a
fragment with empty `t_inds` (they only project to transcripts via
`t_indices`). Confirmed by the existing skip‑on‑zero‑candidates pattern in
[scan.py FragmentRouter](src/rigel/scan.py); a one‑line guard there is the
worst case.

### Pillar 3 — Categorization

Inputs (per fragment, all in or about to be in the buffer):

| Field | Source |
|---|---|
| `splice_type` | scanner |
| `num_hits` | scanner |
| `genomic_footprint` | scanner (existing) |
| `exon_bp_pos`, `exon_bp_neg`, `tx_bp_pos`, `tx_bp_neg` | scanner (Pillar 2a) |

Filter: `num_hits == 1` AND `splice_type == UNSPLICED`. Splice artifacts and
all spliced classes are excluded.

Categories (vectorized; computed once per chunk):

```
fp                         = genomic_footprint
exon_bp_either_lower_bound = max(exon_bp_pos, exon_bp_neg)
tol                        = exon_fit_tolerance_bp  (default 3)
intron_bp_pos              = max(tx_bp_pos - exon_bp_pos, 0)
intron_bp_neg              = max(tx_bp_neg - exon_bp_neg, 0)

INTERGENIC       : tx_bp_pos == 0  AND  tx_bp_neg == 0
INTRONIC         : exon_bp_pos == 0  AND  exon_bp_neg == 0  AND  not INTERGENIC
EXON_CONTAINED   : (fp - exon_bp_either_lower_bound) <= tol
EXON_INCOMPATIBLE: otherwise
```

Strand sub‑label (for `EXON_CONTAINED`, `EXON_INCOMPATIBLE`, and `INTRONIC`):

```
strand_pos   : *_bp_pos > 0  AND  *_bp_neg == 0
strand_neg   : *_bp_neg > 0  AND  *_bp_pos == 0
strand_ambig : *_bp_pos > 0  AND  *_bp_neg  > 0
```

(For `INTRONIC` the `*_bp` is `intron_bp`; for the exon classes it is `exon_bp`.)

`CalibrationResult.category_counts` becomes a 2‑D array of shape
`(N_CATEGORIES, 4)` — last axis is `(none, pos, neg, ambig)` — JSON‑serialized
as a flat list with documented stride.

### Pillar 4 — Pool assembly

```
Pool = INTERGENIC ∪ INTRONIC ∪ EXON_INCOMPATIBLE
Length = chunk.genomic_footprint  (no clipping)
```

Out‑of‑range fragments (`fp < 0` or `fp > max_size`) are **dropped** with a
diagnostic counter `n_pool_dropped_out_of_range` on `CalibrationResult`. We
do **not** silently saturate to bin `max_size`.

`UNSPLICED_*_EXONIC` (now strand sub‑labels of `EXON_CONTAINED`) remain
**excluded from the pool**. The math is unchanged from v1 (antisense‑exonic
is RNA‑dominated unless SS ≈ 1.0).

> **Optional v2.1 follow‑up — SS‑gated antisense pool.** When the trained
> strand model reports `SS > 0.999` (true for the VCaP series), the
> `EXON_CONTAINED, strand_ambig` and `EXON_CONTAINED, strand_neg` (relative
> to the library protocol) buckets are nearly pure gDNA candidates and
> **dwarf** the current pool (4M vs 29k for the pure VCaP library). This
> is the single largest lever for improving pure‑library calibration
> sensitivity. Out of scope for v2 itself; design noted here so we don't
> rebuild the same plumbing later.

### Pillar 5 — Splice artifact handling (consolidated with Pillar 2 in scanner)

Add `SpliceType.SPLICE_ARTIFACT = 4`. When the SJ blacklist filter at
[bam_scanner.cpp:178-181](src/rigel/native/bam_scanner.cpp#L178-L181)
(`n_sj_blacklisted`) rejects all CIGAR `N` junctions for a fragment, set
`splice_type = SPLICE_ARTIFACT` instead of downgrading to `UNSPLICED`.

For v2, **do not** preserve per‑block disjunction in the buffer schema —
that is a deeper plumbing change. Emit one buffer fragment per the
existing path; `genomic_footprint` will reflect the merged span (which is
itself a useful diagnostic). Categorization filter
`splice_type == UNSPLICED` excludes these from calibration entirely.
The per‑locus EM continues to see them as multimapper‑like ambiguity.

### Pillar 6 — Implicit splicing (DEFERRED)

The prior plan's Phase 4 (`SPLICED_IMPLICIT` for paired‑end gaps spanning
annotated introns) is deferred to a v2.1 follow‑up. Rationale:

- Highest implementation complexity in the plan; touches the inner BAM
  per‑pair processing path.
- **Not on the critical path for the observed bin‑0 bug.** The pure‑library
  fitted gDNA mode is 0 because of the `-1` sentinel issue, not because of
  implicit splicing leak.
- A v2.1 implementation can be added without revisiting any of v2.

Document as a known v2 limitation: paired‑end fragments with an unsequenced
gap that hides an annotated intron are still classified `UNSPLICED` and may
enter the gDNA pool (modest effect; bounded by the rate of inserts > 2×read
length, typically a few %).

### Pillar 7 — Region partition stays gone

(Unchanged from prior plan.) `region_df`, `region_cr`, `RegionAccumulator`
are not coming back.

## 3. Bug → fix mapping

| Bug | Root cause | Fix |
|---|---|---|
| `gDNA_FL` mode = 0 (impossible) | `-1` sentinels → `np.clip` → bin 0 | Pillar 1: pool length = `genomic_footprint`, never per‑candidate |
| 4–6% mass at L=1000 | `np.clip(..., max_size)` saturating | Pool drops out‑of‑range with diagnostic counter |
| `INTERGENIC` / `INTRONIC` Python counts always 0 | resolver early‑exits zero‑candidate fragments | Pillar 2b: stop dropping; let them flow as INTERGENIC |
| Pool dominated by `EXON_INCOMPATIBLE` artifacts | (downstream of above) | Pillar 3: real INTERGENIC and INTRONIC populate the pool |
| Splice artifacts pollute pool | blacklisted CIGAR‑N → UNSPLICED | Pillar 5: new `SPLICE_ARTIFACT` class held out |
| High‑contamination underestimation | `gDNA_FL` shape (mode 0) suppresses γ_f | Verify at Phase 5 — expected to improve |

## 4. nRNA‑pool‑pollution: accepted v2 limitation

The 1‑D mixture `pool = π·gDNA + (1−π)·RNA_FL` cannot distinguish gDNA from
nascent RNA when their fragment‑length distributions overlap, which they
typically do (both span 100–600 bp from sonication). nRNA fragments enter the
pool (INTRONIC almost entirely; EXON_INCOMPATIBLE partly), inflating fitted
`π`.

This is a **structural identifiability limit**, not a bug — fixing it
requires a strand‑aware mixture or a third component, both of which are
out of scope for v2. We surface a new `n_pool_intronic_strand_pos` /
`n_pool_intronic_strand_neg` diagnostic so the user can spot nRNA‑heavy
libraries (extreme strand asymmetry inside the INTRONIC bucket given high SS).

## 5. Phased execution

### Phase 0 — Audit (no code changes; ~1 hour)

Verify and produce `docs/calibration/srd_v2_audit.md` with file:line refs:

1. Confirm `compute_frag_lengths` skips entries with projected `fl ≤ 0`
   ([resolve_context.h:756](src/rigel/native/resolve_context.h#L756)). ✅ already done in this plan.
2. Confirm `from_core` fills missing entries with `-1`
   ([resolve_context.h:165-168](src/rigel/native/resolve_context.h#L165-L168)). ✅
3. Confirm the `np.clip` call in `_simple.py` ([_simple.py:130-132](src/rigel/calibration/_simple.py#L130-L132)). ✅
4. Confirm `genomic_footprint` is universally populated and `≥ 0` for every
   resolved fragment (write a small one‑shot `scripts/debug/audit_genomic_footprint.py`
   that walks `mctp_vcap_rna20m_dna00m`'s buffer and reports `min`, `max`, and
   `n_zero` of `genomic_footprint`).
5. Confirm the resolver's "drop zero‑candidate fragments" path; quantify how
   many fragments would gain INTERGENIC categorization if the drop is removed
   (cross‑check against the C++ `n_intergenic_unspliced` counter).
6. Locate the SJ blacklist application path
   ([bam_scanner.cpp:289-291](src/rigel/native/bam_scanner.cpp#L289-L291)
   and the per‑read `n_sj_blacklisted` field at line 181) and the path that
   downgrades to `UNSPLICED`.
7. Quantify, on the same buffer, the fraction of pool fragments with
   `cand_frag_lengths[first] == -1`. (This is the "predicted impact" of
   Phase 1 alone.)

### Phase 1 — Calibration‑only safety patch (no scanner changes)

Two file edits, ~20 lines total, no recompile required.

**[_categorize.py](src/rigel/calibration/_categorize.py):**
- Replace `frag_length[valid_first] = cand_frag_lengths[idx]` (~line 144)
  with `frag_length = chunk.genomic_footprint.astype(np.int32, copy=False)`.
- Drop the `read_length` fallback at ~line 149.

**[_simple.py](src/rigel/calibration/_simple.py):**
- Replace the `np.clip` line with explicit out‑of‑range filtering and a
  `n_pool_dropped_out_of_range` diagnostic.
- Add the diagnostic field to `CalibrationResult`.

**Tests:** rerun the spike series. Acceptance (Phase 1 only):
- `mctp_vcap_rna20m_dna00m`: `gDNA_FL` mass in bin 0 < 1%, mass in bin
  `max_size` < 1%, mode in `[100, 400]`.
- All eight libraries: report new `pi_pool` and headline `gdna_fraction`;
  expect pure‑library `pi_pool` to drop from 0.42 toward 0.05–0.20 (real
  pool size is now what's in EXON_INCOMPATIBLE minus the `-1` artifacts).
- **No claim** about absolute `gdna_fraction` correctness.

This phase ships immediate value and is the gate for Phase 2.

### Phase 2 — Scanner overhaul (consolidated; one C++ build)

Single recompile combining Pillar 2a, Pillar 2b, and Pillar 5.

**Files:**
- [resolve_context.h](src/rigel/native/resolve_context.h): add `exon_bp_pos`,
  `exon_bp_neg`, `tx_bp_pos`, `tx_bp_neg` to `RawResolveResult`. Augment
  the inner loop (lines 880‑915) per Pillar 2a. Remove the
  zero‑candidate early exit (Pillar 2b).
- [bam_scanner.cpp](src/rigel/native/bam_scanner.cpp): thread the four new
  fields into `FragmentAccumulator` and the buffer‑fill path. Add the
  `SPLICE_ARTIFACT` classification on the blacklist‑rejection path.
- [buffer.py](src/rigel/buffer.py): add the four new columns to
  `_FinalizedChunk` and `from_raw`. Update `memory_bytes`.
- [types.py](src/rigel/types.py): add `SpliceType.SPLICE_ARTIFACT = 4`.

**Tests:**
- `tests/test_resolve_strand_overlap.py` — synthetic mini‑annotation with
  fragments crossing exon/intron/intergenic on +/− strands; assert exact
  bp counts; include the strand‑overlap zone case; verify the invariants
  in §2 Pillar 2a.
- `tests/test_zero_candidate_fragments.py` — synthetic intergenic
  fragments now reach the buffer with empty `t_inds` and zeroed strand
  counts; downstream consumers (FragmentRouter, locus builder) tolerate
  them.
- `tests/test_splice_artifact.py` — blacklisted SJ → `SPLICE_ARTIFACT`;
  not downgraded to `UNSPLICED`.

### Phase 3 — Categorization rewrite (Python only)

**Files:**
- [_categorize.py](src/rigel/calibration/_categorize.py): rewrite per
  Pillar 3. Stop reading `cand_frag_lengths`, `t_offsets`, `t_indices`,
  `intron_bp`, `exon_bp`. Read only the four new columns plus
  `splice_type`, `num_hits`, `genomic_footprint`. Add strand sub‑labels.
- [_result.py](src/rigel/calibration/_result.py): widen `category_counts`
  to `(N_CATEGORIES, 4)` last‑axis = `(none, pos, neg, ambig)`. Add
  `n_pool_dropped_out_of_range`, `n_pool_intronic_strand_pos`,
  `n_pool_intronic_strand_neg` diagnostics.
- [_simple.py](src/rigel/calibration/_simple.py): pool =
  `INTERGENIC ∪ INTRONIC ∪ EXON_INCOMPATIBLE`. Remove the
  `frag_length_models.intergenic` merge (Pillar 2b makes those fragments
  reach the buffer directly).

**Tests:**
- `tests/test_categorize_v2.py` — synthetic per‑category and per‑strand
  cases (4 categories × 4 strand sub‑labels = 16 base cases) plus
  tolerance edge (overhang of exactly `tol` and `tol + 1`) plus
  strand‑overlap zone.
- Update `tests/test_calibration_simple.py` for the new schema.

### Phase 4 — Validation (no code; data + report)

Rerun the spike series end‑to‑end:

```bash
bash /scratch/.../rigel/srd_v1_vcap_mixture/run_all.sh
```

Output to a sibling `srd_v2_vcap_mixture/default/` directory and produce
`docs/calibration/srd_v2_results.md` with side‑by‑side tables.

**Quantitative acceptance.**

| Metric | Library | Threshold |
|---|---|---|
| `gDNA_FL` bin‑0 mass | all | < 1% |
| `gDNA_FL` bin‑max mass | all | < 1% |
| `gDNA_FL` mode | all | in `[100, 400]` bp |
| `gDNA_FL` mode consistency | spiked libs | within ±50 bp across `dna01m..dna80m` |
| Headline `gdna_fraction` monotonicity | dna00m → dna80m | strictly increasing |
| Headline `gdna_fraction` at dna80m | spike truth ≈ 80% | ≥ 0.75 (improvement vs current 0.69) |
| Headline `gdna_fraction` at dna00m | "pure" RNA | report only; no pass/fail |
| Pool size growth from Pillar 2b | dna00m | report (expect 10×–100× with INTERGENIC + INTRONIC populated) |
| `n_pool_dropped_out_of_range` | all | < 0.1% of pool |
| Scanner runtime delta | full sample | < 5% |

If the dna80m headline does **not** improve, the structural nRNA limitation
is suspect. Drill in with the new `INTRONIC × strand` diagnostic.

### Phase 5 — Cleanup

- Remove the C++ `intergenic_obs / intergenic_truth / intergenic_lengths`
  accumulator path in `bam_scanner.cpp` (Pillar 2b makes it redundant).
- Remove the corresponding `frag_length_models.intergenic` plumbing in
  Python (`pipeline._replay_fraglen_observations`).
- Update `MANUAL.md`, `METHODS.md`, `CHANGELOG.md`.

## 6. Risk register (revised)

| Risk | Mitigation |
|---|---|
| Phase 1 surfaces a different artifact (e.g. `genomic_footprint` itself buggy) | Phase 0 step 4 explicitly checks |
| Pillar 2 strand‑bp accounting wrong | Invariant tests in `test_resolve_strand_overlap.py` |
| Removing zero‑candidate drop breaks downstream (locus builder, scorer) | Phase 0 step 5 enumerates all consumers; one‑line guard at most |
| Scanner overhead from 4 new fields | Profile in Phase 4; expect <5% (work is constant‑factor on existing hit walk) |
| nRNA contamination still pollutes gDNA_FL | Accepted limitation; new INTRONIC × strand diagnostic surfaces it |
| Pure‑RNA library has no real gDNA signal | Expected; `gdna_fl_quality = "weak"` or `"fallback"` correctly fires |
| `SPLICE_ARTIFACT` block‑disjunction not preserved | Deferred per §2 Pillar 5; merged span is acceptable for v2 |
| `SPLICED_IMPLICIT` still misclassifies gap‑hidden introns as UNSPLICED | Documented v2 limitation; v2.1 follow‑up |

## 7. Out of scope (deliberate)

- Implicit splicing (`SPLICED_IMPLICIT`) — v2.1.
- Splice‑artifact per‑block preservation — v2.1.
- SS‑gated antisense pool inclusion — v2.1, listed as the largest open
  lever for pure‑library sensitivity.
- Strand‑aware 2D mixture for nRNA / gDNA separation.
- Cross‑region density signals; mappability‑corrected effective length;
  region partition resurrection.
- Strand‑model / RNA‑FL training changes.

## 8. One‑page decision summary

- Keep all three diagnosed root causes; the fixes are correct.
- Ship Phase 1 first (calibration‑only patch, ~20 lines, no recompile) —
  it eliminates the bin‑0 artifact on its own.
- Consolidate scanner work (strand‑aware overlap counts + zero‑candidate
  fragments + `SPLICE_ARTIFACT`) into a single C++ phase to avoid two
  recompile/rerun cycles.
- Defer implicit‑splice detection and per‑block splice‑artifact
  preservation to v2.1.
- Tighten Pillar 2's intron accounting: store `tx_bp_*` (not `intron_bp_*`),
  derive `intron_bp_* = max(tx_bp_* − exon_bp_*, 0)` in Python.
- Acknowledge the nRNA‑pool‑pollution structural limit; expose a
  diagnostic; don't pretend the 1‑D mixture can solve it.
- Numerical acceptance criteria are in §5 Phase 4.
