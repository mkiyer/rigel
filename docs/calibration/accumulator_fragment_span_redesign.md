# Accumulator fragment-span redesign — deposit the molecule, not the reads

**Status:** proposed plan (audit + design). 2026-06-09. Supersedes the calibration-side `L_cross`
self-calibration idea in `density_sweep_implementation_plan.md`. **Lands after Phase 0, before
Phase 1.**

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

| Fragment class | Contiguous spans to deposit | Channel |
|---|---|---|
| **gDNA / unspliced** (no explicit or implicit intron) | **one** span `[min block start, max block end]` — fill every gap (mate gap = unsequenced contiguous DNA) | unspliced |
| **Explicitly spliced** (CIGAR `N`) | the exonic segments; cut at each CIGAR-`N` intron; fill within-exon mate gaps | spliced |
| **Implicitly spliced** (annotated intron inside the mate gap, no CIGAR `N`) | cut at the implied intron; fill all other gaps | spliced |
| **Artifact splice** (blacklisted junction, currently treated unspliced) | **open question — §6** | (today: unspliced) |

Mass normalization: `L = Σ contiguous-span lengths` (the molecule's footprint minus introns), **not**
the covered read length — directly the user's "accumulate relative to the full fragment length."
Mass still conserves to 1 per fragment (sum over the regions its spans touch), so the
`total mass == n_fragments` invariant ([[accumulator_mass_conservation]]) holds.

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

## 5. Where to implement (proposed)

A new pure helper that maps `(blocks, explicit introns, candidate transcripts, splice_type)` → a
sorted list of contiguous genomic spans, then feed those spans to the **unchanged** `deposit()`:

```
spans = coalesce_fragment_spans(frag.exons, frag.introns, splice_type, candidate_introns)
acc.deposit(span_starts, span_ends, n_spans, spliced, primary)
```

- **Smallest-blast-radius option.** Put `coalesce_fragment_spans` in the scanner (next to
  `deposit_to_accumulator`, `bam_scanner.cpp`), so `Accumulator::deposit` and its **byte-for-byte
  reference** (`tests/native/_accumulator_reference.py`) are unchanged (they already take spans). The
  scanner-integration tests + golden change (intended).
- **Phasing the work within this fix:**
  1. **Unspliced coalescing** (the entire density bias): if `splice_type ∈ {UNSPLICED}` and not
     implicit → deposit `[min,max]`. Trivial, unconditional for the intergenic/gDNA path. **Highest
     value, lowest risk — could land alone first.**
  2. **Spliced spans:** cut at explicit + implicit introns, fill within-exon mate gaps. Needs the
     intron intervals (explicit from `frag.introns`; implicit recovered from the matched transcript).
  3. **Channel correctness:** an implicitly-spliced fragment is RNA — it should deposit on the
     **spliced** channel (today `has_introns()` puts it on unspliced). Decide whether to pass
     `splice_type`-derived `spliced` to the deposit so channel + span treatment agree.

Alternative (rejected for now): fold coalescing into `deposit()` itself and pass introns through the
ABI. More invasive (touches the spec + reference); revisit only if a non-scanner caller needs it.

---

## 6. Open questions for review

1. **Artifact splices** (`SPLICE_ARTIFACT`, blacklisted CIGAR `N` treated as unspliced). If the `N` is
   a true alignment artifact the molecule is contiguous → fill across it (consistent with "unspliced
   ⇒ one span"). If it is a real-but-blacklisted junction, filling deposits intronic coverage. Which
   default? (Leaning: treat as unspliced ⇒ fill, matching the rest of the pipeline's handling, and
   accept that a blacklisted junction is rare and low-impact.)
2. **Channel for implicit splices.** Move them to the spliced channel (correct biology, changes
   strand/splice statistics slightly) or keep current behavior and only fix the span? Recommend
   fixing both together for consistency.
3. **Interchromosomal / chimeric paired fragments.** A mate on a different ref must **not** be filled
   across (already separate refs; the span coalescing is per-ref). Confirm the per-ref grouping in
   `build_fragment` covers this (it keys by `(ref_id, strand)`).
4. **Performance.** `coalesce_fragment_spans` is O(blocks) on an already-sorted, tiny block list
   (≤ a few per fragment) in the hot per-fragment path; expected negligible. Profile against the
   current deposit on `gdna_benchmark_5mb`.

---

## 7. Validation & acceptance

- **Unit (new):** `coalesce_fragment_spans` — unspliced 2-block fragment → one `[min,max]` span;
  explicitly-spliced → exonic spans with intron cut; implicit-intron-in-gap → cut; within-exon mate
  gap → filled; interchromosomal → not merged.
- **Oracle regression:** rerun `phase1_fl_geometry_derivation.py` — `flux / oracle-template-span → 1.0`
  at BOTH FL≈350 and FL≈180 (the geometry dependence of the over-count disappears); contained-ρ and
  crossing-ρ agree on uniform density.
- **Conservation:** `total mass == n_fragments` still holds (the Phase-0 validation script extended).
- **Byte-for-byte:** `Accumulator::deposit` unchanged ⇒ `_accumulator_reference.py` unaffected; add a
  reference for `coalesce_fragment_spans`.
- **Golden + suite:** regenerate gDNA/combo goldens (calibration inputs shift), full suite green,
  `evaluate_suite.py` on `gdna_benchmark_5mb` shows the leak drop carry through with NO calibration
  correction factor.

## 8. Sequencing (updated master order)

Phase 0 (intergenic deposit, ✅ shipped) → **this accumulator span fix** (unspliced coalescing first,
then spliced spans) → Phase 1 (density estimator, now unbiased) → Phases 2–4 as in
`density_sweep_implementation_plan.md`.

## Reproduce
```bash
python scripts/debug/phase1_fl_geometry_derivation.py   # the over-count derivation/harness
```
