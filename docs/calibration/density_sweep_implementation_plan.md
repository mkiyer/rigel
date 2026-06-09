# Count-clue density sweep: implementation plan

**Status:** implementation plan (DESIGN → IMPLEMENTATION). 2026-06-08. Companion to
`density_sweep_audit_and_redesign.md` (the design/root-cause). This doc breaks the rework into
phases, records a **dry-run** of the core algorithm (with evidence + issues encountered), and
fixes the implementation-level decisions.

---

## 0. Decisions locked by the dry-run

`scripts/debug/impute_prototype.py` prototyped the replacement and compared two density estimators
against oracle truth, per region:

| estimator | GENE0037 exons (capture, true ~21) | antisense run (no capture, true ~15) |
|---|---|---|
| CURRENT global sweep | 0.62 (≈30× low) | 2.46 |
| A1 = side gDNA **mass** / `boundary_side_eff_len` | 6–8 (≈3× low) | 7.2 |
| **A2 = side gDNA *count* / `fl_mean`** | **17–23 (≈unbiased)** | 23 (≈1.5× high) |

**Decisions:**
1. **The density estimator is count-based: `crossing_gDNA_count / fl_mean`** (A2), not mass-based.
   The per-side *mass* under-counts (overlap-weighted partial bases); the crossing *count* over
   `fl_mean` (the boundary count effective length) is ~unbiased in the capture case that drives the
   leak. (The ~1.5× over in the no-capture case is a residual constant — see §6 open items.)
2. **Imputation is local + directional.** A non-observable region takes anchors from its
   observable-boundary sides; run interiors (no observable side) inherit the nearest anchored
   neighbour via a left→right + right→left carry. This recovered the right magnitude for single-side
   anchors and run interiors in the dry-run.
3. **Phase 0 (intergenic deposit) is a prerequisite** — until intergenic contained mass is deposited,
   intergenic regions read 0 (no baseline), and the global density / fallback are deflated.

---

## Phase 0 — Fix the intergenic accumulator deposit (C++) *(prerequisite)* — ✅ SHIPPED

> **DONE (2026-06-09).** Implemented in `src/rigel/native/bam_scanner.cpp`: factored the resolved
> deposit into a single shared `deposit_to_accumulator` lambda in `process_qname_group_threaded`
> and call it from **both** the resolved unique-mapper path and the intergenic (`resolved &&
> is_unique_mapper`) path. Confirmed `align_strand`/`sj_strand` are populated from the alignment
> blocks regardless of transcript resolution (`resolve_context.h:1034`), so intergenic fragments
> orient correctly. **Validation** (`scripts/debug/phase0_intergenic_deposit_validation.py`,
> oracle scenario): intergenic regions went 0 → contained gDNA mass matching oracle **exactly**
> (189,949 = 189,949, rel = 0.000); mass conservation exact (193,117 contained + 6,883 boundary =
> 200,000 = n_fragments). Full suite green (994 passed); golden regenerated for the 4 gDNA
> scenarios only (`gdna_light`, `gdna_heavy`, `combo_moderate`, `combo_extreme`; Δ 0.4–1.2% on
> counts, toward more gDNA recognized). The byte-for-byte `deposit()` reference was unaffected (the
> ABI is unchanged — only which fragments the scanner deposits changed). The "unit accumulator" and
> "_accumulator_reference.py" acceptance items below were therefore N/A; the oracle validation
> subsumes them.

**Root cause (confirmed).** In `src/rigel/native/bam_scanner.cpp` the accumulator `deposit()` (line
~1578-1608, `if (ws.acc_set && !frag.exons.empty())`) runs **only inside the resolved-hit path**.
Intergenic fragments (`any_hit_resolved == false`, lines ~1643-1649) only increment
`n_intergenic_{spliced,unspliced}` telemetry and are **never deposited** → all 101 intergenic regions
carry 0 contained mass (3.19 Mb / ~64% of genome), measured.

**Fix.** Deposit unresolved (intergenic) **unique-mapper, non-chimeric** fragments into the
accumulator using their alignment blocks, with the same orientation gating as the resolved path
(defined `align_strand`; for spliced also `sj_strand`). `Accumulator::deposit` already maps genomic
block positions → region indices, so an intergenic fragment lands in its intergenic region's
contained mass (or the gene-flanking boundary if it crosses). Implementation notes:
- The deposit must see the fragment's alignment blocks + `ref_id` + strand in the intergenic branch.
  Capture the same `frag.exons` / `result.align_strand` used by the resolved deposit; the intergenic
  branch is per-fragment (not per-hit), so thread the single unresolved alignment through.
- Multimappers: keep them out (or 1/NH) consistently with the resolved convention; simplest is to
  deposit only unique-mapper intergenic fragments (matches the resolved deposit's unique-mapper gate).
- This also restores intergenic **strand** information for the strand-balance / overdispersion fits.

**Tests / acceptance.**
- Unit (accumulator): a synthetic intergenic fragment increments the correct region's `contained`.
- Byte-for-byte: update `tests/native/_accumulator_reference.py` to include the intergenic deposit.
- Regression: on `gdna_benchmark_5mb`, intergenic regions now carry contained mass ≈ oracle
  (`gene37_region_boundary_autopsy.py`: R368/R383 cal ≈ true ~7.5k); global density rises off the
  deflated floor. Golden regenerated.
- Cost: a single per-fragment deposit; expect negligible runtime change (verify with the profiler).

> **Caveat:** Phase 0 changes calibration inputs → golden + benchmark numbers shift. Land it first,
> on its own, with golden regenerated and the suite green, so later phases build on conserved mass.

---

## Phase 1 — Replace the sweep with local boundary-anchored imputation

Rewrite `density_model.node_gdna_density` so the per-region gDNA **density** is:
- **observable region** (intron/intergenic): own contained gDNA `count / region_eff_len`
  (strand-cleaned). Post-Phase-0 the intergenic baseline is real.
- **non-observable region** (exon/AMBIG): the **A2** estimate from each *observable*-boundary side —
  `side_gDNA_count / fl_mean` — averaged over available sides; run interiors filled by the
  directional carry (L→R then R→L, combined). Remove the `from_left`/`from_right` undamped
  accumulation and the `w = flux/(flux+1)` conduit entirely.

Observability per side: region `r`'s **left** side (`substrate.left[r]`) is usable iff boundary
`(r-1, r)` is observable; its **right** side (`substrate.right[r]`) iff boundary `(r, r+1)` is. Reuse
`count_observable_masks`, `boundary_eff_length` (`fl_mean`), `boundary_side_eff_length`.

Strand-clean each side's unspliced counts to gDNA via the closed form (same as the contained path);
AMBIG sides fall back to gDNA-fraction 1 (or defer — see Phase 4).

### The FL-geometry constant — DERIVED (2026-06-09)

The dry-run's open "~1.5× over in no-capture" residual is now fully characterized
(`scripts/debug/phase1_fl_geometry_derivation.py`, oracle-validated on a uniform-density genome):

1. **`fl_mean` is unbiased.** Estimated gDNA `fl_mean` / true = **1.000** (Phase 0's intergenic
   contained deposits give the FL pipeline a clean, large-region sample). The residual is *not* an
   FL-estimate error.
2. **The crossing physics is exact.** Oracle gDNA *templates* spanning a seam = ρ·`fl_mean` to <1%
   (1172 vs 1174). So `#templates covering a point = ρ·fl_mean` holds — the count→density model is
   right.
3. **The residual is a deposit over-count.** The accumulator's one-side boundary **flux** =
   **c × (true template-span count)**, with **c ≈ 1.56** at gDNA FL≈350. Cause: paired-end reads are
   deposited as *two blocks* (read1, read2) with the insert gap between them; the slice-adjacency
   crossing rule (`accumulator.cpp` §4) credits a side-event whenever consecutive slices change
   region, which over-counts relative to a single template spanning the seam point (a read straddling
   a boundary plus a mate in the same region creates an extra crediting; a mate-gap that skips an
   interior region credits flanking seams asymmetrically).
4. **c is NOT universal — it tracks read/insert geometry.** Halving the gap (FL≈350→180) moves
   **c 1.56→1.19** (→1 as the gap vanishes). So c must be **self-calibrated per library**, never
   hardcoded. This also explains why the benchmark's GENE0037 (capture, short fragments) read ≈
   unbiased while the FL≈350 synthetic read 1.56× over.
5. **c is uniform within a library** (std ≈ 0.02–0.06 across all boundary types, exon-adjacent and
   intron/intergenic-adjacent alike) and **capture-independent by mechanism** (capture changes which
   fragments exist, not how a 2-block fragment is sliced).

**Resolution — self-calibrate, no magic number.** Tie the crossing estimator to the *exact* contained
estimator using count-observable nodes as the calibration anchor. Define one fitted per-library
scalar — the **empirical boundary effective length** `L_cross` (equivalently `c = L_cross / fl_mean`):

```
ρ̂_region (observable)   = strand_cleaned_contained_count / region_eff_len          # exact reference
L_cross = Σ_obs-boundary (one-side gDNA flux)  /  Σ_obs-boundary (ρ̂ · 1)            # fit, weighted by evidence
        # i.e. the flux a unit local density produces at a seam, given THIS library's read geometry
ρ̂_region (non-observable) = (boundary-side gDNA flux) / L_cross                      # apply to exons
```

`L_cross` is **fit, not assumed**: it is the slope that makes seam flux agree with the contained
density on observable nodes, so the deposit over-count cancels exactly whatever the read geometry.
Use the flux→`L_cross` path (not the per-side *mass*→`boundary_side_eff_length` path — the dry-run's
A1, biased ~3× low because a straddling fragment puts only part of its length on each side). Pool the
fit across all observable seams in the library (it is geometry-, not locus-, specific); guard against
loci with no observable seam via the Phase-4 fallback.

> **Alternative considered (not chosen): fix the over-count at the C++ source** so flux counts true
> template-spanning. Rejected for Phase 1 — it changes accumulator semantics + the byte-for-byte
> reference, the mass channel is already correct, and self-calibration is robust to *any* systematic
> seam/geometry effect, not just the mate gap. Revisit only if a second consumer needs a clean flux.

**Outputs unchanged in shape:** `NodeDensity.density` + `count_evidence` (the latter revised in
Phase 3). The joint deconvolution (`joint_deconv.py`) is untouched in Phase 1.

**Tests / acceptance.**
- Unit: a single exon between two introns recovers the (known) uniform density from either/both
  sides; a 3-region antisense run recovers the interior region from the carry.
- Derivation regression: `L_cross`-calibrated crossing-ρ / contained-ρ ≈ 1.0 on the uniform-density
  scenario at BOTH FL≈350 and FL≈180 (proves the self-calibration removes the geometry dependence).
- Worked: GENE0037 exon density 0.62 → ~20; antisense run interior → right magnitude.
- The dry-run code (`impute_prototype.py`) + `phase1_fl_geometry_derivation.py` are the references.

---

## Phase 2 — DNA-fraction estimator at validated splice boundaries

Add the second estimator **(B)** and use it **only** where the boundary is a true splice
donor/acceptor (an intron ending/starting at an exon, **same strand**), so spliced reads genuinely
flow across it and `spliced` ≠ "no RNA here":

```
boundary is splice-junction ⇔ one side has the exon bit and the other the intron bit of the SAME strand
   (e.g. ex+|in+ seam, or in+|ex+) — NOT intergenic↔exon (TSS/TES: no splicing), NOT ex+|ex- (antisense)
dna_fraction = gdna_unspliced_crossing / (gdna_unspliced_crossing + spliced_crossing)
region gDNA  = region_total_contained · dna_fraction   (applied to the WHOLE region)
```

**Why gate it (your point 2):** at a non-splice boundary we would see 0 spliced + some unspliced and
wrongly call it 100% DNA. Only annotated splice seams (validated from the signature transition) admit
a fair DNA-vs-RNA ratio. **Why it can beat density:** the fraction is a same-position DNA/RNA
comparison, robust to within-exon probe-enrichment gradients, and applies to the region's own total
counts (independent of where probes sit in the exon).

Implementation: compute the per-region DNA-fraction estimate where the bracketing seams qualify;
otherwise fall back to the A2 density estimate (always available). Keep **both** estimates and a
config/diagnostic switch so the **empirical comparison** (next) can select per-condition.

**Empirical comparison framework.** Extend `impute_prototype.py` into a sweep over conditions
(capture on/off × stranded/unstranded × alt-splice × antisense × gDNA level) that reports, per
region, `true` vs `A2-density` vs `B-fraction` vs combined, and the resulting per-locus prior error.
Decide the selection rule (e.g. fraction at qualified seams, density elsewhere) from the data.

---

## Phase 3 — Honest prior concentration + IPR effective-length uncertainty

1. **Concentration (`count_evidence`).** For a non-observable region set the count-prior concentration
   to the **observed crossing-fragment count** that supports its estimate (real evidence), not
   `density · region_eff_len`. Thin evidence ⇒ weak prior ⇒ strand governs. (`joint_deconv.py:94-97`
   reads `count_evidence`; feed it the supporting count.)
2. **IPR effective length (`priors.assemble_priors`)** — the make-or-break divisor. It is built from
   the per-region gDNA mass, so Phases 1–2 (correct mass) fix most of it automatically (gDNA mass
   concentrates on captured exons → correct contraction). Additionally:
   - propagate imputation uncertainty as **prior strength**, not a hard divisor — an unresolved region
     contributes little to both the prior count and the eff-len;
   - keep the Laplace `+1` (uninformative → geometric span); for genuinely-unidentifiable regions
     regularize toward a **local** prior length (geometric span as null; exon/spliced footprint as the
     informed-under-capture prior), never a global.

**Tests:** per-locus `gdna_eff_len_em` tracks the captured footprint on GENE0037; unstranded/no-gDNA
loci keep span (no spurious contraction).

---

## Phase 4 — Fallbacks for the unidentifiable floor

Opposite-strand overlap (AMBIG) at low strand-specificity is the genuine floor; goal is **rational**,
not accurate, behavior. Fallback order for a region with no usable anchor:
1. nearest observable boundary/region density (the run carry);
2. an **intron baseline** density (post-Phase-0, observable);
3. only as a last resort, a conservative length with **near-zero** concentration (add little to the
   prior; let the EM + strand decide) — **never** the global average (deflated; mixes depleted +
   enriched).
- **No capture:** the neighbor/baseline fallback is sound (density ~uniform).
- **Capture:** the fallback is intrinsically murky (captured vs not is unknown without probes); accept
  reduced accuracy here. This case requires *both* unstranded *and* capture to bite, so it is rare.

---

## Dry-run: issues encountered (implementation-readiness)

1. **Estimator normalization:** mass-based (A1) is biased ~3× low; **count-based (A2) is the form to
   ship.** Residual ~1.5× over in no-capture suggests cross-side contamination at the seam (the
   exon-side count includes some neighbor-origin crossings) and/or the exact `fl_mean` constant — needs
   a short FL-geometry derivation + empirical calibration; track as an acceptance item, not a blocker.
2. **Intergenic = 0 until Phase 0:** the prototype shows intergenic A1/A2 = 0 (no contained data) — the
   baseline/fallback is impossible until Phase 0 lands. Confirms phase ordering.
3. **Degenerate zero-`eff_len` regions** (e.g. GENE0037 R374, an alt-splice boundary-shift artifact)
   appear; the imputation must guard against div-by-eff_len and skip/zero them.
4. **AMBIG strand-cleaning** is undefined (no valid sense split); the prototype used gDNA-fraction 1.
   Phase 4 must decide AMBIG handling explicitly (likely defer to strand at the joint, weak prior).
5. **Per-side strand orientation:** a boundary spans two regions of possibly different strand; the
   prototype cleaned each side by the region-`r` strand. Verify this is correct for antisense seams.
6. **Run detection:** the simple "no observable side ⇒ carry from neighbors" two-pass fill works for
   short runs; verify on 4–5-region runs and conflicting left/right anchors.

---

## Sequencing, validation, and acceptance

**Order:** Phase 0 (C++, golden churn, land alone) → Phase 1 (count-based local imputation) →
Phase 2 (DNA-fraction + comparison) → Phase 3 (concentration + eff-len) → Phase 4 (fallbacks). Each
phase: unit tests + golden regenerate + `evaluate_suite.py` on `gdna_benchmark_5mb` + full suite green.

**Acceptance (end state):**
- GENE0037 exon swept density ~20 (not 0.62), `count_gf` ~0.8, joint `gdna_frac` ~0.8.
- Benchmark net gDNA→RNA leak in capture-on / high-gDNA drops sharply; **no** rise in gdna_none false
  gDNA; capture-off conditions unchanged.
- Intergenic contained mass ≈ oracle; per-locus eff-len contracts to the captured footprint.
- Per-region imputation validated against oracle across the full signature/condition matrix (Phase 2
  comparison harness), including runs and the unidentifiable floor (rational fallback).

**Risks:** (a) the residual estimator constant (manageable — bounded, validated empirically);
(b) golden churn across phases (mitigated by phase isolation); (c) AMBIG/unstranded edge cases (accept
reduced accuracy, ensure rational fallback); (d) Phase-0 runtime (expected negligible; profile).

## Open items
- ~~Exact FL-geometry constant for the count→density estimator (the ~1.5× no-capture residual).~~
  **RESOLVED (2026-06-09):** it is a paired-end deposit over-count `c = flux/template-span`, geometry-
  dependent (1.19–1.56), so **self-calibrated per library** as `L_cross` against the contained
  reference — see Phase 1 "The FL-geometry constant" + `phase1_fl_geometry_derivation.py`.
- Selection rule between A2-density and B-fraction (decide from the Phase-2 comparison).
- AMBIG handling at the joint (weak prior + strand) — Phase 4.
- Confirm the intergenic-deposit fix conserves total fragment mass (accumulator invariant).

## Reproduce
```bash
SUITE=/Users/mkiyer/Downloads/rigel_runs/gdna_benchmark_5mb
python scripts/debug/build_stress_overlap_scenario.py
python scripts/debug/impute_prototype.py --index $SUITE/rigel_index \
  --bam $SUITE/gdna_gdna400_ss_0.99_nrna_none_capture_on/sim_oracle.bam \
  --ref chr_syn --start 1659774 --end 1701056
python scripts/debug/impute_prototype.py --index /tmp/rigel_sim/antisense_run/index \
  --bam /tmp/rigel_sim/antisense_run/antisense_run_oracle.bam --ref antisense_run --start 4500 --end 8500
```
