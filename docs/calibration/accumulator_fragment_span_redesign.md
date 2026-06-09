# Accumulator fragment-span redesign — deposit the molecule, not the reads

**Status:** design FINALIZED + phased implementation plan (§7). 2026-06-09. All four review questions
resolved (§6). Supersedes the calibration-side `L_cross` self-calibration idea in
`density_sweep_implementation_plan.md`. **Lands after Phase 0, before Phase 1.** Ready to start
**Phase A (test scaffolding, no production change)** on go-ahead.

> **TL;DR.** The fractional accumulator deposits a fragment's **aligned read blocks**, not the
> **molecule's true genomic extent**. For a paired-end fragment the unsequenced insert gap between
> read1 and read2 is part of the molecule but is deposited as a hole between two blocks. The
> slice-crossing rule then over-counts boundary flux by ~1.2–1.56× (geometry-dependent) and
> mis-attributes mass. The fix is to deposit each fragment's **contiguous genomic span(s)** — fill
> the mate gap for unspliced molecules; cut out only true introns for spliced molecules — and
> normalize mass over the full span. This removes the bias **at the source**, so the Phase-1 density
> estimator (`flux / fl_mean`, `contained / region_eff_len`) is unbiased with no correction factor.

---

## 1. The evidence (why this is real)

Oracle-validated on a uniform-density (no-capture) genome
(`scripts/debug/phase1_fl_geometry_derivation.py`):

| gDNA FL | `fl_mean` est/true | oracle templates spanning seam / ρ·fl_mean | one-side flux / true span-count |
|--------:|:------------------:|:------------------------------------------:|:-------------------------------:|
| 351     | 1.000              | 0.999                                      | **1.558** (std 0.06)            |
| 216     | 1.000              | 0.992                                      | **1.192** (std 0.02)            |

- The FL estimate is exact and the count→density physics is exact (oracle template-span = ρ·fl_mean).
- The boundary **flux** over-counts true template-spanning by a factor that **shrinks as the mate gap
  shrinks** (1.56 → 1.19 → 1 as FL → read-sum). That is the signature of a **paired-end mate-gap**
  artifact, not a statistical or FL effect.

This also explains the field observation that the benchmark's GENE0037 (capture, short fragments) read
≈ unbiased while a long-fragment synthetic read 1.56× over: the bias magnitude *is* the read/insert
geometry.

---

## 2. Audit — what the accumulator does today

**Fragment assembly** (`bam_scanner.cpp:617 build_fragment`):
- Collects per-read aligned blocks (`ParsedAlignment.exons`, the CIGAR `M` segments) from read1 +
  read2 into `exon_dict`, and CIGAR-`N` junctions (`rec.sjs`) into `intron_set`.
- Merges blocks **only when overlapping/adjacent** (`bam_scanner.cpp:673`, `intervals[i].first <=
  cur_end`). A paired-end **insert gap** (read2.start > read1.end, no overlap) is **not** merged →
  two separate `ExonBlock`s with a hole between them. The gap is recorded as **neither** an exon nor
  an intron.
- `has_introns()` = "has an explicit CIGAR-`N` junction" — **false** for an unspliced paired fragment
  with a mate gap, and also false for an **implicit** splice (intron sitting in the mate gap, no
  CIGAR `N`).

**Deposit** (`calibration/accumulator.cpp:84 deposit`):
- Receives the block list + a `spliced` flag (`= has_introns()`, explicit-only) + `primary`.
- Expands each block into per-region **slices**; if all slices share one region → `contained += 1`;
  else **crossing**: each slice credits `flux_left` of its right boundary (if not last) and
  `flux_right` of its left boundary (if not first); mass is split per side by covered-base share over
  `L = Σ slice lengths` (the **covered** length, not the span).
- FL pooling already bins at the **span** `hi − lo` (`accumulator.cpp:147`), a prior fix — so the FL
  pmf is correct; only the flux/mass geometry is wrong.

**The bug, traced.** Unspliced paired fragment, read1 `[1100,1250]` straddling a boundary at `1200`,
read2 `[1300,1400]` in the downstream region B (boundary B|C at `1500`):
- blocks `[1100,1250]`, `[1300,1400]`; slices `A[1100,1200]`, `B[1200,1250]`, `B[1300,1400]`.
- slice0 `A` → `flux_left[1200]`. slice1 `B` → `flux_right[1200]` **and** `flux_left[1500]` (PHANTOM —
  the molecule never reaches C). slice2 `B` → `flux_right[1200]` **again** (double).
- Net: boundary `1200` right side credited twice, a phantom credit at `1500`. The design invariant
  "slices are monotonic, each side credited once" (`accumulator.cpp:189`) is **violated** because
  read1's tail and read2 both produce slices in region B with a hole between them.

Coalescing the unspliced molecule to a single span `[1100,1400]` gives slices `A[1100,1200]`,
`B[1200,1400]` → `flux_left[1200]`, `flux_right[1200]`, once each. Bias gone.

---

## 3. Correct behavior — deposit contiguous genomic spans

A fragment occupies a set of **contiguous genomic intervals** = its true molecular footprint. Deposit
*those*, and compute mass over their total length.

There is **one unifying rule**: *walk the sorted blocks and fill every inter-block gap UNLESS the gap
is (or contains) an intron.* Everything below is that rule specialized by fragment class — the only
inputs are the blocks and the set of **cut introns** (gaps to keep open).

| Fragment class | Cut introns (gaps kept open) | Resulting spans | Channel |
|---|---|---|---|
| **gDNA / unspliced** (no explicit or implicit intron) | none | **one** span `[min start, max end]` — fill every gap (mate gap = unsequenced contiguous DNA) | unspliced |
| **Explicitly spliced** (CIGAR `N`) | the CIGAR-`N` introns (`frag.introns`) | exonic segments; within-exon mate gaps filled | spliced |
| **Implicitly spliced** (annotated intron in the mate gap, no CIGAR `N`) | the PE gap(s) that matched an annotated intron | cut at the implied intron; other gaps filled | **spliced** (fix: today deposited as unspliced — a bug) |
| **Artifact splice** (`SPLICE_ARTIFACT`, blacklisted junction) | — | **HOLD OUT — not deposited at all** (decision §6.1) | — |
| **Chimeric** (multi-locus / interchromosomal) | — | already **held out** today (deposit is after the chimeric `continue`); deferred (§6.3) | — |

Mass normalization: `L = Σ contiguous-span lengths` (the molecule's footprint minus introns), **not**
the covered read length — directly the user's "accumulate relative to the full fragment length."
Mass still conserves to 1 per *deposited* fragment (sum over the regions its spans touch); with
artifacts held out the invariant is `total mass == n_deposited` ([[accumulator_mass_conservation]]).

**Why this is principled, not a band-aid.** The accumulator's job is to measure where fragment
*molecules* sit on the genome. A molecule is its full insert (for DNA/unspliced) or its exonic body
(for spliced RNA). Depositing read blocks instead conflates "what we sequenced" with "where the
molecule is." Fixing it makes both the count clue (`flux/fl_mean`, `contained/region_eff_len`) and
the mass clue exact by construction — the over-count factor → 1.0 with no calibration knob.

---

## 4. The implicit-splice machinery already exists

The resolver already builds PE gaps and tests them against annotated introns
(`resolve_context.h:788 has_implicit_splice_gap`, `:750 transcript_has_implicit_intron_in_gap`,
`GapBlock`, `SPLICE_IMPLICIT`). The gap-construction loop (`:817-828`) is exactly the
"sort blocks, find holes where `block.start > cur_end`" we need. So we do **not** invent intron
detection — we reuse it to decide, per gap, *fill vs cut*.

Key consequence: the **resolved** path has candidate transcripts (`cr.t_inds`) → can identify which
gaps are implicit introns. The **intergenic** path has no transcript → no intron is possible → all
gaps fill (the molecule is contiguous DNA). This neatly matches: the gDNA density bias is entirely in
the intergenic/unspliced population, where the fix is unconditional gap-fill.

---

## 5. Survey — where the pieces are, and where the new logic lives

The work is **code organization, not new computation**: every input the span rule needs is already
produced once per fragment. The principle is *reuse the existing iterations; add no fragment-level
loop; keep each function small and single-purpose.*

**Current per-fragment flow (the hot path), and what each stage already yields:**

| Stage | Location | Already produces |
|---|---|---|
| `build_fragment` | `bam_scanner.cpp:617` | sorted, within-read-merged `exons` (blocks) + explicit-`N` `introns`; one merge loop over blocks |
| `_resolve_core` | `resolve_context.h:~1000` | `splice_type` (incl. `SPLICE_IMPLICIT` at `:1345`), `t_inds`, chimera type; **`has_implicit_splice_gap` (`:788`) already builds the PE gaps and matches them to annotated introns** |
| `deposit_to_accumulator` | `bam_scanner.cpp:1355` | currently loops blocks into `bs/be` arrays → `deposit()` |

**The single primitive (pure, unit-tested):**

```cpp
// Coalesce sorted blocks into contiguous genomic spans, keeping a gap OPEN
// iff it is (or contains) a cut intron. O(blocks + cuts) on tiny sorted inputs.
void fragment_genomic_spans(const std::vector<ExonBlock>& sorted_blocks,
                            const std::vector<IntronBlock>& cut_introns,  // explicit ∪ implicit
                            std::vector<Span>& out_spans);   // reused scratch, no alloc in hot path
```

Rule: walk blocks; extend the current span across a gap unless a `cut_intron` lies within that gap
(then close the span and start a new one). A different `ref_id` always closes the span (interchrom\
osomal → separate spans, future-proofing chimeras for free). Unspliced ⇒ no cuts ⇒ one span.

**Assembling the `cut_introns` with no new loop:**
- *Explicit*: `frag.introns` (already built).
- *Implicit*: extend `has_implicit_splice_gap` → `implicit_splice_gaps(...)` that **records the
  matched gaps** (the loop already finds them at `:833-839`; today it discards them and returns
  `bool`). Stash the matched intervals on `RawResolveResult`/`ResolvedFragment` so the deposit reads
  them. The `bool` becomes "non-empty."

**Where it runs:** inside `deposit_to_accumulator`, **replacing** the existing `bs/be`-building loop
with `fragment_genomic_spans(...)` → `deposit()`. Same call site, same iteration count — *not a new
pass*. `splice_type` (threaded via `result`/`ig_result`) drives the channel (`spliced` for explicit +
implicit) and the `SPLICE_ARTIFACT` hold-out guard. `Accumulator::deposit` is **unchanged** (it still
takes spans), so its byte-for-byte reference is untouched.

**Files touched (small, focused):**
- `bam_scanner.cpp`: new free `fragment_genomic_spans` + rewire `deposit_to_accumulator` (channel from
  `splice_type`; artifact guard).
- `resolve_context.h` / `constants.h`: `implicit_splice_gaps` records matched intervals; carry them +
  ensure `splice_type` reaches the deposit (already on `ResolvedFragment`).
- No change to `Accumulator::deposit` or `_accumulator_reference.py`.

---

## 6. Resolved decisions

1. **Artifact splices → HOLD OUT of the accumulator.** We cannot recover the true span: a blacklisted
   `N` may be a real-but-blacklisted junction *or* a wholly incorrect alignment, and we derive the
   span from the alignment — so any reconstruction risks a false-positive assumption about the true
   molecule. The cleanest, most honest, lowest-bias rule is a single guard: `if splice_type ==
   SPLICE_ARTIFACT → do not deposit` (and contribute no FL pool mass). They are rare, so the data loss
   is negligible, and we introduce **zero** extrapolation. (Considered and rejected: deposit raw
   blocks = today's behavior, which reintroduces the multi-block over-count; deposit `[min,max]`,
   which assumes contiguity we can't justify; keep-largest-block, which is arbitrary and under-counts
   mass. Hold-out wins on simplicity and correctness.)
2. **Implicit splices → SPLICED channel + cut at the implied intron.** They are spliced reads; we
   already detect and label them `SPLICE_IMPLICIT`. Depositing them as unspliced (today) is a bug —
   it both fills across a real intron and mis-channels RNA as unspliced. Fix both together.
3. **Chimeras → confirmed already held out, deferred.** The resolved-chimeric path appends to the
   buffer and `continue`s *before* the deposit (`bam_scanner.cpp:~1505`), and intergenic fragments
   (empty `t_inds`) cannot be chimeric — so no chimera reaches the accumulator today. The span rule's
   "different ref closes the span" already makes future inclusion correct, but it stays **out of scope**
   for this work.
4. **Performance → integrate, don't bolt on.** No new fragment-level loop: the span computation
   replaces the existing `bs/be` loop in `deposit_to_accumulator`, the cut-intron set reuses
   `build_fragment`'s introns and the resolver's already-running implicit-gap loop, and
   `fragment_genomic_spans` is O(blocks+cuts) on ≤-few-element sorted lists with reused scratch (no
   allocation in the hot path). Acceptance: profile vs current on `gdna_benchmark_5mb`, require no
   measurable scan-phase regression.

---

## 7. Phased implementation plan

> **Discipline (per the user):** an initial phase puts the *entire correctness harness* in place —
> unit tests, regression guards, oracle assertions — **before any production code changes**, with the
> new tests written to the target behavior so they **fail** against today's code. Only then do
> implementation phases flip them green, one behavior at a time, each with golden regenerate + full
> suite green.

**Phase A — Correctness scaffolding (tests first; NO production change).**
- Unit tests for `fragment_genomic_spans` (table-driven): unspliced 2-block → one `[min,max]`;
  explicit-spliced → exonic spans, intron cut; implicit-intron-in-gap → cut; within-exon mate gap →
  filled; multi-ref → separate spans; artifact → caller skips.
- A C++/Python reference for the primitive (mirrors the chosen semantics; analogous to
  `_accumulator_reference.py`).
- Oracle assertion test built from `phase1_fl_geometry_derivation.py`: `flux / oracle-template-span`
  ≈ 1.0 at FL≈350 **and** FL≈180 (xfail/skip-marked until Phase B, documenting the target).
- Channel-correctness test: an implicit-splice fragment lands in the **spliced** channels (ch2/3),
  not unspliced (xfail until Phase C).
- Mass-conservation test extended for hold-out: `total mass == n_deposited` (artifacts excluded).
- Confirm the full existing suite is green first (regression baseline) and enumerate the golden files
  expected to change (gDNA/combo) vs not (RNA/nRNA-only).

**Phase B — Span primitive + unspliced coalescing** (the entire gDNA density bias).
- Implement `fragment_genomic_spans`; rewire `deposit_to_accumulator` to deposit coalesced spans;
  unspliced/intergenic ⇒ one `[min,max]` span. Recompile, regenerate golden, suite green.
- Flip the oracle over-count assertion green (1.56/1.19 → ~1.0).

**Phase C — Implicit splice → spliced channel + intron cut.**
- `implicit_splice_gaps` records matched intervals; thread to deposit; `splice_type`-driven channel;
  cut at implicit introns. Flip the channel-correctness test green. Golden regenerate, suite green.

**Phase D — Artifact hold-out.**
- Guard `SPLICE_ARTIFACT` out of the deposit + FL pool. Golden regenerate, suite green.

**Then:** Phase 1 of `density_sweep_implementation_plan.md` (density estimator) — now unbiased, no
correction factor.

## 8. Validation & acceptance (whole effort)

- Oracle: `flux / oracle-template-span → 1.0` across FL geometries; contained-ρ ≈ crossing-ρ on
  uniform density.
- Conservation: `total mass == n_deposited` (artifacts excluded) — extends the Phase-0 validator.
- Byte-for-byte: `Accumulator::deposit` + `_accumulator_reference.py` unchanged; new reference for
  `fragment_genomic_spans`.
- Golden + suite: gDNA/combo goldens regenerate, RNA/nRNA-only unchanged, full suite green;
  `evaluate_suite.py` on `gdna_benchmark_5mb` shows the leak drop carried by the accumulator fix with
  **no** calibration correction factor.
- Performance: no measurable scan-phase regression (profile).

## Reproduce
```bash
python scripts/debug/phase1_fl_geometry_derivation.py   # the over-count derivation/harness
```
