# Effective-length audit — calibration vs EM

**Date:** 2026-04-22
**Triggered by:** user concern that the calibration estimates gDNA
density per **mappable** bp, while the EM might re-normalise per
**total** bp — which would systematically under-estimate gDNA.

This document is a read of the code as-is; **no behaviour changes
yet.**

## TL;DR

1. The locus-level **EM solver uses no length term at all.** `θ_c`
   is a pure mixing weight over components (mRNA transcripts,
   synthetic nRNA, one gDNA bin), scatter-normalised to sum to 1.
   `em_solver.cpp:1332`: "Effective lengths are all 1.0, so
   `log_eff_len = 0`." Therefore the EM itself cannot introduce a
   mappable-vs-total bias — it doesn't know about lengths.
2. The **per-fragment log-likelihoods** `p(F | c)` that feed the EM
   also carry no spatial/length term. Per-component emission is
   `p(fraglen | c) × p(splice | c) × ½ × 1/NH` with no
   `-log(transcript_eff_len)` and no `-log(genome_size)` /
   `-log(mappable_bp)`. Both RNA and gDNA are treated the same way
   in this respect. Scoring is therefore length-neutral.
3. The **calibration** *does* compute a physical rate:
   `λ_G = (Σ gDNA-attributable n_total) / (Σ mappable_bp)`, fit by
   the v4 mixture EM. Units: fragments per mappable base.
4. The **bridge** between calibration and EM is the per-locus
   mixing fraction γ_ℓ. In
   `locus.py::compute_locus_priors`:

   ```
   γ_ℓ  =  Σ_i region_e_gdna[i]     /     Σ_i region_n_total[i]
       =  Σ_i (λ_G × mappable_bp_i) /     Σ_i n_total_i
   ```

   Both numerator and denominator sum over the SAME set of
   calibration regions that overlap the locus's merged intervals.
   γ_ℓ is then mapped to `α_gDNA = γ_ℓ·c_base`, `α_RNA =
   (1-γ_ℓ)·c_base`, and used (a) as a Dirichlet pseudocount and
   (b) as the EM warm-start ratio
   `θ_gDNA = (α_gDNA/α_RNA) × θ_RNA_total`.

   **Key observation:** γ_ℓ is a *ratio* computed on the same
   region set with the same counting convention. Units (and
   length conventions) cancel: if you were to change mappable_bp
   to total_bp in *both* the numerator `λ_G·mbp_i` and in the
   calibration fit of `λ_G`, γ_ℓ would be unchanged (to first
   order — see §3 caveats). So the calibration→EM handshake is
   **not** an obvious source of underestimation.

## 1. What the EM sees (length-neutral)

### 1.1 Per-component log-lik

- RNA transcript `t` for fragment `F`:
  ```
  log L_RNA(F, t) = log p(fraglen_F | RNA_fl_model)
                  + strand_llr
                  + splice / overhang / mismatch penalties
                  + log(1/NH)
  ```
  *No* `-log(eff_len_t)`.

- gDNA bin for fragment `F` (unspliced only):
  ```
  log L_gDNA(F) = log p(fraglen_F | gDNA_fl_model)
                + gdna_splice_penalty (0 for unspliced)
                + log(½)                [unstranded]
                + log(1/NH)
  ```
  `scoring.cpp:910-913`. *No* `-log(genome_size)` or
  `-log(mappable_bp)`.

### 1.2 EM update

```
posterior(F, c)  ∝  θ_c · exp(log L(F, c))
θ_c_new          ∝  unambig_totals[c] + Σ_F posterior(F, c) + prior[c]
```

`em_solver.cpp::assign_posteriors` line 1332 explicitly zeros
`log_eff_len`. The EM is a Dirichlet-multinomial mixture with no
rate/length parameter anywhere. It cannot under- or over-count
gDNA based on a length-convention mismatch — because it doesn't
consume lengths at all.

## 2. What the calibration sees (mappable-bp rate model)

`_em.py::run_em` fits a two-component Poisson mixture on the
full-genome partition:

- **Class G** (gDNA-only): `E[count_i] = λ_G · mappable_bp_i`
- **Class R** (expressed): `E[count_i] = λ_G · mappable_bp_i + μ_R_i`

The denominator is always `mappable_bp` (≡
`region_df["mappable_effective_length"]` if present, else raw
`length`). So λ_G's units are **frags per mappable base**.

## 3. Where mappable vs total bp could still bite

Even though the γ_ℓ ratio is internally consistent, there are
three subtler ways the mappable-bp convention can introduce bias:

### 3.1 Residual bias in λ_G itself

The calibration counts `n_total_i` are derived from BAM-scan
per-region aggregates (`region_counts`). These aggregates
*already* include whatever multi-mapper / NH fractional
assignments that landed in region `i`. If region `i` has
`mappable_bp = 0` (entirely low-mappability repeats) but
`n_total_i > 0` (some multimappers happened to be scattered
there), this region contributes to the mixture EM marginal
likelihood without having a finite rate — in practice `_em.py`
will handle it via numeric safeguards, but the λ_G estimator
gets nudged upward by "un-rateable" observations. **Not obviously
large; worth quantifying.** Diagnostic: plot
`n_total_i / mappable_bp_i` in high-γ regions and measure the
residual density vs λ_G fit.

### 3.2 Non-overlap of locus intervals with calibration partition

`compute_locus_priors` sums over calibration regions that
overlap `locus.merged_intervals` (union of transcript spans).
If a locus is entirely *inside* one calibration region (common:
short loci inside an intergenic bin), only that region
contributes. If the locus straddles multiple regions, the sum
is area-weighted by each region's `mappable_bp` and `n_total`.

**Edge case:** a locus whose merged interval includes a strip of
zero-mappability sequence (a tandem repeat inside an intron)
gets no contribution from that strip (both e and n are
mappable-weighted by construction). γ_ℓ for that locus is
identical to what you'd get if the repeat weren't there.

That is *arguably correct* (we have no expected or observed
counts there) but it means: if the EM has to place fragments
that nonetheless fell into that repeat strip (via multi-mapper
allocation), there is no prior support for them being gDNA.
They'll be pushed toward RNA. **This is a plausible underestimation
pathway.** Diagnostic: for the dna80m mega-locus, measure how
many EM fragments have their genomic footprint inside
mappable_bp=0 sub-intervals.

### 3.3 Mismatch between calibration partition and locus partition

`n_total_i` counts in the calibration partition are tallied
from `region_counts` (BAM-scan output). `n_em_fragments`
counts in a locus are fragments that resolve to that locus's
transcripts or its gDNA bin. These two populations are related
but not identical — a multimapper might appear in both counts.
The γ_ℓ ratio doesn't "know" about this: it trusts that the
relative shares it measures (e_i/n_i) per region are
representative of the locus's population.

If RNA-rich loci preferentially attract NH-split fragments from
nearby intergenic regions (unlikely but possible), γ_ℓ could be
biased low. **Lowest-probability pathway; won't investigate
further unless other hypotheses fail.**

## 4. What the user proposed: "use total length everywhere"

The user floated this alternative: in the EM, extrapolate gDNA
density across the entire genomic region length rather than
mappable_bp. Two ways to implement:

**Option A.** Change calibration to fit λ_G per *total bp*.
Then:
- Per-region `region_e_gdna = λ_G_total · length_i`.
- Per-locus γ_ℓ = Σ e / Σ n_total.
- No other change.

This is mechanically equivalent to the current code up to a
constant rescaling (mappable_bp ↦ length) IF both n_total and e
use the same denominator in both calibration and prior — which
they do. **γ_ℓ is unchanged.** The only effect is λ_G's unit
interpretation.

**Option B.** Keep calibration with mappable_bp (correct
physics: fragments can only come from mappable bases), but
have the EM *also* consume mappable_bp as an effective length
for the gDNA component, so θ_gDNA and θ_RNA are rates-per-length
rather than pure mixing weights.

This is a substantial change: it would turn the locus EM from
Dirichlet-multinomial (mixing) into a Poisson-rate model (rate
× length). The RNA side would then also need `eff_len_t` per
transcript, which is currently discarded by setting
`log_eff_len = 0`. This is closer to a Salmon-style model.

**Recommendation:** Option A is cosmetic (no γ_ℓ change).
Option B is a rewrite of the EM likelihood and of the scoring
interface. Neither option is expected to fix the 10-20 pp
underestimate on its own.

## 5. What I think actually explains the underestimate

From the corrected P1 analysis (all size bins 10-20 pp under
expected), the remaining candidates in order of plausibility
are:

1. **`min(λ_G·mbp, n_total)` clip** — inflates `region_e_gdna`
   asymmetrically. Fix: drop the clip. Still the cheapest
   experiment.
2. **§3.1 λ_G residual bias** — test by pulling high-γ regions
   and recomputing λ_G via pooled MLE over just those regions.
3. **§3.2 unmappable-interior fragments** — measurable on the
   dna80m mega-locus by intersecting each EM fragment's
   genomic footprint with mappable-bp=0 sub-intervals.
4. **nRNA siphon** — growing with locus size; real but not
   dominant at small loci.

## 6. Refactor performed in this pass (unrelated to density)

Independently, the variable/column name `mrna` in the C++→Python
bridge was misleading because the C++ accumulator
`assign_posteriors::mrna_total` actually sums posteriors over
**all transcript-like components including synthetic nRNA
shadows.** This has now been renamed through the whole call
chain:

| Layer | Old | New |
|---|---|---|
| C++ arg `assign_posteriors` | `mrna_total` | `transcript_total` |
| C++ vec/ptr | `locus_mrna_vec/data` | `locus_transcript_vec/data` |
| C++ return tuple name | `locus_mrna` | `locus_transcript_total` |
| Python `pipeline._build_locus_meta` kwarg | `mrna=` | `transcript_total=` |
| `locus_results` dict key | `"mrna"` | `"transcript_total"` |

The user-facing **`loci.feather` column `mrna`** is unchanged; it
continues to report *annotated-only* mRNA mass after the
annotated/synthetic split in `AbundanceEstimator.get_loci_df`.
The P1 double-count accounting fix (total = transcript_total +
gdna) is preserved.

All 931 tests pass.
