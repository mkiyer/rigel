# Phase 3 — strand-resolved nascent removal (3-term splice fraction + per-strand carry-over): implementation plan

> **SUPERSEDED (2026-06-11) by [`sequential_calibration_redesign.md`](sequential_calibration_redesign.md).**
> The direct-3-term-on-eligible-regions approach here was measured net-negative (washed out on
> strand-observable regions; never reaches the count-routed AMBIG node — see the A/B + the toy scenario
> `scripts/debug/phase3_ambig_scenario.py`). The redesign instead restructures into a **sequential
> strand→count pipeline** (strand cleans first + communicates precision; count imputes the cleaned
> density field), which reaches the AMBIG node by construction. Kept for the design rationale (3-term
> arithmetic, per-strand identifiability, sizing).

**Status:** execution-ready plan, 2026-06-11. Resolves the open questions in
`count_mean_bias_design.md` §5.1/§8a.5 and `remaining_phases.md` §3.8 with the decisions below.
Supersedes those sketches for *implementation*; they remain the design rationale. **Do not start the
build until the two prerequisites (§9) land: the shared run-fill refactor and the FL-consistency
diagnostic.**

## 0. What it is + honest scope

The shipped 2-term boundary fraction (`splice_junction.region_splice_gdna_frac`) lumps crossing-unspliced
reads as **gDNA + nascent**, so nascent pre-mRNA is counted as gDNA. Phase 3 emits a **3-term** fraction
that moves nascent to the RNA side wherever a gDNA/RNA split of the unspliced is available:

```
f = (U_gDNA / E^gDNA) / ( U_gDNA/E^gDNA  +  (U_RNA + S)/E^RNA )      # pure gDNA / total
```

the split coming from **strand pre-cleaning** (sense/antisense vs κ deconvolves the crossing unspliced
into gDNA at sense ½ and nascent at sense κ) or from **carry-over** along a left→right sweep, with a
**per-strand `{gDNA, RNA⁺, RNA⁻}`** state for identifiability at antisense (AMBIG) overlaps.

**Honest expectations (measured, do not over-claim).** On the purpose-built `quick_3to1_5mb` suite
(antisense 0.15 + random nascent + capture) the target class is small: 146 AMBIG regions (2.13% of
genomic span), 92 AMBIG-&-non-observable (0.82%), and only **11 regions with exon overlap on both
strands** (0.05%) — the true §5.1 `E+·I−` identifiability case. And the 3-term only carries **full
weight on count-routed (AMBIG / unstranded) nodes**: in a stranded library `w=(2κ−1)²→1`, so the strand
module already classifies nascent as sense-RNA and the count-fraction 3-term gets only ~`(1−w)≈4%`
weight (blend washout); in an unstranded library the split is unavailable and the 3-term degrades to the
2-term (the irreducible floor). **So Phase 3 is correctness/robustness for complex regions, not a
headline-number mover.** We build it anyway because the human genome is dense with antisense/overlapping
loci where the method must not silently fail — "works most of the time, fails in a minority" is the
case Phase 3 closes.

## 1. Design decisions (resolved)

1. **Per-strand state from the start.** The sweep state is the vector `{gDNA, RNA⁺, RNA⁻}` (mirrors the
   accumulator's 4 channels — no new information, just don't sum away the unspliced strand). A scalar
   stopgap is *not* built; the per-strand state is where AMBIG identifiability lives and is cheap. The
   per-region prior still collapses to the scalar `gDNA/(gDNA+RNA⁺+RNA⁻)`.
2. **Carry BOTH the fraction and the density; the density is the FLOOR.** The absolute-density imputation
   (`node_gdna_density.count_gdna_frac`) is always available and is the **lower bound** on a region's
   gDNA share; the splice-derived fraction may **only raise** the estimate above it:
   ```
   g_region = max( g_density_floor,  g_splice )        # density floor, splice wins above
   ```
   Rationale: under hybrid capture the absolute density *under-counts* gDNA (boundary depletion), so the
   splice fraction normally sits above the floor and wins; the floor is an FP-safe catch that prevents a
   noisy or over-aggressive splice estimate from calling **less** gDNA than the raw density actually
   shows. This replaces today's unconditional `frac[i] = splice_frac[i]` with the `max`.
3. **Seeds = direct splice anchors only.** Because the density floor (decision 2) already covers every
   region universally, the carry-over sweep does *not* also seed from count-observable density regions
   (that would just re-derive the existing density run-fill). The sweep propagates only the
   **splice-derived, nascent-clean** signal; the density is an independent, universal floor underneath.
   (Resolves §8a.5-Q3.)
4. **Acyclic order — reorder the splice step after the overdispersion fits.** The 3-term needs κ and the
   fitted gDNA/RNA strand overdispersions to strand-split the crossing unspliced. Move the
   `region_splice_gdna_frac` call in `calibrate.py` from its current position (before the fits) to after
   `fit_gdna_strand_from_substrate` / `fit_rna_strand_from_substrate`. This is clean: those fits consume
   `node_density` (the un-upgraded density), **not** the splice-upgraded fraction, so nothing breaks.
   Order becomes: `κ → density → gDNA/RNA overdispersion → splice 3-term (consumes κ+overdispersion+carry)
   → blend`. Assert no density→clean→density loop.

## 2. The arithmetic

**The split (per boundary side).** A boundary side carries `n_unspliced_pos`, `n_unspliced_neg`
(genome-strand crossing-unspliced) and `n_spliced_*` (the mature reference). Reuse the **existing strand
Beta-Binomial posterior** (`strand_likelihood.strand_loglik`, the same machinery `strand_deconv` uses) to
get the gDNA fraction `g_bb` of the crossing unspliced from `(sense, antisense)` given κ + the fitted
overdispersions — **no new estimator**. Then `U_gDNA = g_bb·(n_pos+n_neg)`, `U_RNA = (1−g_bb)·(n_pos+n_neg)`,
and the per-strand nascent split follows from the sense orientation: `nascent⁺/nascent⁻` partition the
`(1−g_bb)` RNA mass by the side's pos/neg counts above the gDNA-symmetric ½ baseline.

**The 3-term fraction** then uses `boundary_gdna_fraction(unspliced_gdna=U_gDNA, unspliced_rna=U_RNA,
spliced=S, eff_gdna, eff_rna)` — the function already supports the 3-term signature (it takes
`unspliced_rna`); today the caller passes `unspliced_rna=0` (2-term). Phase 3 passes the deconvolved split.

**AMBIG identifiability (the core §5.1 case).** At an opposite-strand overlap the local unspliced is
```
n_pos = gDNA/2 + nascent⁺ ,   n_neg = gDNA/2 + nascent⁻        # 2 eq, 3 unknowns — under-determined
```
Carrying `nascent⁺` forward from a `+`-only neighbour (where it *was* identifiable) closes it:
`gDNA = 2·(n_pos − nascent⁺)`, then `nascent⁻ = n_neg − gDNA/2`. The per-strand state is what makes this
expressible; a scalar carry cannot represent two independent RNA flows crossing the same overlap.

## 3. The per-strand carry-over sweep

- **Run** = a maximal stretch of consecutive same-reference regions (the structure `density_model`
  already uses; **consume the shared `runfill_bidirectional` helper** from the §9 refactor — Phase 3 is
  its second caller, in fraction/state space).
- **State per region:** `{gDNA, RNA⁺, RNA⁻}` (densities, so they compose with effective lengths).
- **Seeds:** regions with a direct strand-splittable splice anchor (decision 3).
- **Carry:** bidirectional forward+reverse fill of the state along the run; `RNA⁺`/`RNA⁻` propagate along
  their transcript extents; at an AMBIG overlap the carried strand(s) combine with the local
  `unspliced_pos`/`unspliced_neg` to solve for `gDNA` (§2).
- **Collapse + floor:** per region, `g_splice = gDNA/(gDNA+RNA⁺+RNA⁻)`, then `g_region = max(g_density,
  g_splice)` (decision 2). Priority: direct splice anchor > carried splice > density floor.

## 4. Files & changes

- **`calibration/run_fill.py` (NEW — the §9 refactor, prerequisite):** `same_ref_left_right(ref_id)` and
  `runfill_bidirectional(values, ref_id)`. Phase 3 imports both.
- **`calibration/splice_junction.py`:** extend `region_splice_gdna_frac` (or a sibling
  `region_splice_gdna_frac_3term`) to (i) strand-split each anchoring side's crossing unspliced via
  `strand_loglik` (κ + overdispersions passed in), (ii) build the `{gDNA,RNA⁺,RNA⁻}` state, (iii) run the
  carry-over via `runfill_bidirectional`, (iv) emit the 3-term fraction, (v) apply the `max(density_floor,
  splice)` combination. Precompute eligibility once per boundary (the 256-entry signature LUT from the
  §9 refactor) instead of the current per-region double-classification.
- **`calibration/calibrate.py`:** reorder per decision 4; pass `rna_sense_frac`,
  `gdna_strand_overdispersion`, `rna_strand_overdispersion`, `n_grid` into the splice step.
- **No config, no CLI, no C++/accumulator change** (the substrate already exposes per-side
  `n_unspliced_pos/neg` + `n_spliced_*`). **No new constant** (§3.9).

## 5. Acyclicity

`κ → density → gDNA/RNA strand overdispersion → splice 3-term → blend → derive`. The 3-term **consumes**
already-fitted κ/overdispersion + the carried state and the un-upgraded density floor; it does **not**
refit the density or feed back into the overdispersion fits. Add an assertion/test that the splice step
reads only `node_density` (floor) + fitted scalars, never the post-splice fraction.

## 6. Unit tests

`tests/calibration/test_splice_junction.py` (+ realized-scenario extensions):
1. **3-term recovers pure gDNA** when nascent is present and a split is supplied (extends the realized
   `boundary_gdna_fraction` cases).
2. **AMBIG identifiability:** constructed `E+·I−` overlap; carrying `nascent⁺` from a `+`-only neighbour
   recovers `gDNA` and `nascent⁻` (the 2-eq/3-unknown resolution).
3. **Carry-over sweep:** a run with one seeded region propagates the split to non-seeded neighbours;
   unreachable runs keep the fallback floor.
4. **Density-floor rule:** splice fraction above the floor wins; a splice fraction driven *below* the
   floor is caught at the floor (FP-safe). q-independent.
5. **Reduction:** zero-nascent or zero-antisense input reproduces the current 2-term output (regression
   stability); unstranded (κ=½) degrades to the 2-term/floor.

## 7. Validation

- `quick_3to1_5mb` (antisense + nascent + capture), 3-pool net flow: the **stranded nascent-on**
  conditions' net gDNA→nRNA siphon on AMBIG regions should shrink toward the nascent-off baseline
  (per-node gDNA estimate matches the oracle *gDNA-only*, not gDNA+nascent). Unstranded+nascent unchanged
  (irreducible floor). No regression on non-AMBIG regions; mature-RNA FP stays controlled. Because the
  target mass is ~0.05–0.82% of span, report the metric **restricted to AMBIG / nascent regions**, not
  just the suite-wide rollup, or the effect will be invisible in the aggregate.
- Optionally materialize `complex_benchmark_25mb` (reference exists; condition dirs not generated) for a
  larger AMBIG sample — but watch disk (each 3:1 condition ≈ 1 GB; delete `sim_oracle.bam`/`*.fq.gz`
  after eval).

## 8. Risks & subtleties (flag, don't bury)

1. **The density floor can re-inject nascent.** The absolute density is computed from *raw* unspliced
   counts (gDNA + nascent), so `g_density_floor` itself contains nascent contamination at nascent-bearing
   imputed regions. When the 3-term correctly removes nascent and drops *below* that contaminated floor,
   `max(floor, splice)` re-injects it. This is the deliberate FP-safe stance (never under-call gDNA below
   observed density), but it bounds Phase 3's nascent removal. **Recommendation:** measure how often the
   floor binds on nascent regions in validation; if it materially limits the benefit, consider flooring
   at a *nascent-aware* density (the strand-split density) rather than the raw one — a follow-up, not v1.
2. **Blend washout (§0):** the 3-term moves the count fraction, which is down-weighted to ~`(1−w)` on
   strand-observable nodes. Set expectations: the visible movement is on count-routed AMBIG nodes.
3. **Multi-junction double-count** in the RNA overdispersion fit (existing, `gdna_strand.py:~364`,
   `rna_strand_bb_double_count_followup` memory) feeds the strand split; unchanged here but noted.

## 9. Prerequisites & sequencing (do these BEFORE finishing Phase 3)

Per the agreed order:
1. **Shared run-fill refactor** (`run_fill.py`): extract `same_ref_left_right` + `runfill_bidirectional`
   from the 6 current copies (canonical: `density_model.py:247-259`); vectorize the fwd/rev fill once
   (segmented `np.maximum.accumulate`-style) so all callers benefit and Phase 3 isn't the 7th copy. Plus
   the per-boundary eligibility LUT for `splice_junction`.
2. **Dead/vestigial cleanup + vocabulary** (Tier 3): unused `StrandSummary.minor_rate_*` /
   `signed_strand_contrast_se`, tests-only `signature` scalar helpers, stale "decode"→"deconv" comments,
   `__init__` docstring, `_transport_boundary_flux` placement.
3. **Efficiency** (optional, whole-genome): vectorize `_deconv_per_node` grid posterior; gate/vectorize
   the O(R²) `_loess` (inert at the default `q=0.5`).
4. **FL-consistency accuracy DIAGNOSTIC** (`count_mean_bias_design.md` §9.1): a dedicated scenario suite
   that *exposes* the `f_b·M_region` vs region `E^gDNA/E^RNA` short-exon error (vary exon length vs the
   gDNA/RNA FL gap), quantify it, **understand before fixing**. The fix (multiplying the region eff-length
   ratio into the partition) is downstream of the diagnostic.

Only after 1–2 (and the 4 diagnostic informing whether the partition needs the eff-length ratio) does the
Phase 3 build (§4) proceed.

## 10. Magic-number audit

None. The 3-term uses already-fit κ, the two strand overdispersions, and the FL effective lengths; the
carry-over reuses the shared run-fill; the floor is a parameter-free `max`. No new constant, config, or
CLI flag.
