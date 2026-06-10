# Phase 4-mean — the splice-junction gDNA-fraction (design)

**Status:** design for review. 2026-06-10. The **bias** fix for the count module's exon gDNA estimate
under hybrid capture (the dominant unstranded+capture leak driver). Composes with — and is independent
of — Phase 4-var (`count_posterior_design.md`, variance only). Supersedes the earlier "structural ratio
`r(length)`" sketch; see §2 for why that framing was wrong.

## 0. TL;DR

At a **true splice-junction boundary**, the crossing fragments split *cleanly by construction*:
- **crossing-unspliced** = gDNA (+ nascent) — mature mRNA cannot contiguously span an exon→intron edge;
- **crossing-spliced** = mature mRNA — only mature carries the junction.

So the boundary yields an observable, **ground-truth-free** gDNA fraction. FL-normalize each crossing
count to a density, take the gDNA share, and use it to partition the *neighbouring exon region's own
(capture-enriched) total* — which fixes the under-count that the absolute-density imputation suffers
under capture. This is the **gDNA-fraction** paradigm (redesign point-5), made precise and per-boundary.

Two refinements make it sharper and composable:
- **§3 eligibility** — exactly which adjacent-region signature pairs are genuine, unambiguous splice
  junctions, *proven* by exhaustive unit tests over the (pruned) 4-bit × 4-bit grid.
- **§5 the 3-term fraction** — when a gDNA/RNA split for the unspliced reads is already available (strand
  pre-cleaning, or carried forward in a sweep), the numerator becomes pure gDNA and nascent moves to the
  RNA side.

Where no boundary is eligible → **fall back** to absolute-density imputation (§6), always available.

> **Validated — 4-mean.1 (eligible 2-term direct anchor) on the 20-condition benchmark (2026-06-10).**
> Unstranded+capture mean leak **15.41% → 11.96%** (−3.45pp); the critical gdna1000 unstranded+capture
> **18.43% → 12.48%** (−5.95pp). No pure-RNA false positives introduced (gdna_none stays 0.0000), no
> regression on stranded (flagship 4.74% → 4.21%) or capture-off. Partial closure as expected — only
> ~29% of non-observable exon regions are junction-bearing; the residual needs the sweep (4-mean.2) and
> 3-term (4-mean.3).

## 1. The mechanism — boundary gDNA-fraction (2-term form)

Per boundary `b` between adjacent regions, and per exon region to impute:
- `U_b`, `S_b` — crossing **unspliced** / **spliced** fragment counts at the boundary.
- `E^gDNA = boundary_eff_length(gdna_pmf)`, `E^RNA = boundary_eff_length(rna_pmf)` — FL-marginal
  crossing effective lengths (`= fl_mean` of each pmf; `effective_length.py`).
- `M_region` — the exon region's **total** fragment count (its own, capture-enriched).

Crossing **densities** (molecules · bp⁻¹), each count ÷ its own FL crossing length:
```
ρ_gDNA = U_b / E^gDNA          ρ_RNA = S_b / E^RNA
f_b    = ρ_gDNA / (ρ_gDNA + ρ_RNA)          # not-mature (gDNA+nascent) molecular share
G_region = f_b · M_region                    # partition the region's own enriched total
```
> **Why this beats absolute-density imputation.** The fallback imputes `G_region = ρ_boundary · E_region`
> — an *absolute* density carried from the boundary, which under capture is **depleted** vs the enriched
> exon interior the probes target (→ the ~2× under-call, 0.41 vs 0.91). The fraction carries only a
> **ratio** `f_b` (robust — enrichment scales gDNA and mature crossing in lockstep) applied to the
> region's **own** enriched `M_region`. Validated: 0.41 → 0.88 (truth 0.885), `diag_spliced_ref_debias.py`.

## 2. Why the earlier "structural ratio r(length)" was overthinking

The diagnostic fit `r = U_mat/S ≈ a·L + b` over region length (Spearman 0.90). That length curve is not a
thing to learn — it is the closed-form FL-marginal effective-length ratio `E_region^gDNA / E_region^RNA`
(`region_eff_length`), which grows with length because a longer region has more interior where a fragment
fits without crossing a junction. **No fit, no learned constant** — the boundary supplies the density
ratio, the FL models supply every length correction in closed form. (Bonus: the boundary split is
ground-truth-free; the interior `U_mat` needed read-name truth to separate mature-exon-body from gDNA.)

## 3. Eligibility — which boundaries are true splice junctions (the crux)

A region's 4-bit signature: `{INTRON_POS=0x8, INTRON_NEG=0x4, EXON_POS=0x2, EXON_NEG=0x1}`. Per strand
`s` define, for a region: `E_s` = has exon on `s`, `I_s` = has intron on `s`.

**Clean splice junction on strand `s`** at boundary `(L,R)` — a pure exon_s↔intron_s transition with no
mixed state on `s`:
```
junction_s(L,R) = ( E_s(L) & ¬I_s(L) & I_s(R) & ¬E_s(R) )      # exon_s → intron_s
               ∨ ( I_s(L) & ¬E_s(L) & E_s(R) & ¬I_s(R) )      # intron_s → exon_s
```
The **exon side** (region to impute) is the side with `E_s & ¬I_s`. With
`junction_strands = {s : junction_s}` and, for the exon-side region, `exon_strands = {t : E_t}`:

> **Eligibility (unstranded, primary target):**
> **`exon_strands(exon-side region) == junction_strands ≠ ∅`**, *and* both junctions (if two) point at the
> same side.

A spliced read's `N`-gap sits at *this* boundary, so crossing-spliced reads come from exactly the strands
that splice here. For `f_b` to reference the *whole* mature exon-body in the region, the strands
contributing mature exon-body (`exon_strands`) must equal those contributing the spliced reference
(`junction_strands`). Mismatches → fallback: AMBIG exon vs single-strand junction (mature⁻ unreferenced);
antisense junction contamination; mixed `E_s & I_s` state (`junction_s` false); split exon sides.

In **stranded** data the strand module already deconvolves exons, so the count path is the small `(1−w)`
residual; with per-strand read splitting the rule relaxes to per-strand `junction_s` — a bonus, not the
target. The matched-set rule is the conservative one that holds with **no** strand info.

**The predicate** `splice_junction_eligibility(sig_L, sig_R) -> (exon_side, strands) | None` is pure, no
parameters. Over the raw 16×16 grid it returns **30 eligible** pairs (6 clean single-/both-strand
junctions + 24 antisense-overlap cases that survive the matched-set rule). Realizability (§4) prunes
further.

## 4. Realizability — the grid is smaller than 256

- **Distinct neighbours.** Adjacent regions have *different* signatures by definition → drop the diagonal
  → **240** (16 × 15) candidate pairs.
- **Splice motifs are asymmetric.** A genomic position is a splice-junction *start* (donor) **or** *end*
  (acceptor), never both, and two distinct genes never share *every* splice coordinate — so generically a
  boundary is a junction on **at most one strand**. The two both-strand matched pairs the signature
  predicate marks eligible (`E+·E- │ I+·I-` and its reflection) therefore **do not arise** for any
  biologically plausible annotation (proven empirically — Test B's generic sweep never produces them).
- **But the coordinate partition cannot see motifs.** `build_region_partition` works from exon/intron
  *coordinates*, not sequence, so if two opposite-strand genes are given *identical* coordinates it *will*
  emit `E+·E- │ I+·I-` (Test B confirms this degenerate case). Rather than special-case it, we proved the
  predicate is **correct anyway**: both strands splice and the region carries both exon bodies, so the
  combined spliced reference matches the combined exon body → matched-set eligible. So the predicate stays
  pure signature logic, needs no motif-asymmetry pruning, and **realizability is a separate, tested
  property** (Test B enumerates what generic vs coincident transcript layouts produce).

## 5. The 3-term fraction — composing with a prior deconvolution

The 2-term `f_b` (§1) lumps crossing-unspliced as gDNA **+** nascent, so it estimates a *not-mature*
fraction. Whenever a gDNA/RNA split of the unspliced reads is **already available**, the numerator becomes
pure gDNA and nascent moves to the RNA side — strictly more accurate:
```
unspliced split into  U_gDNA  and  U_RNA (nascent):
   f = (U_gDNA / E^gDNA) / ( U_gDNA/E^gDNA  +  (U_RNA + S_b)/E^RNA )      # pure gDNA / total
```
Two sources of that prior split (both "deconvolved gDNA + RNA counts before imputation"):

1. **Strand pre-cleaning of an observable boundary.** At an intron→exon boundary the *unspliced* crossing
   fragments (gDNA + nascent, both unstranded-vs-stranded distinguishable) are deconvolved by the strand
   model (total / sense / antisense) into `U_gDNA` vs `U_RNA` (nascent). Then `f = gDNA / (gDNA +
   unspliced_RNA + spliced_RNA)`.
2. **Sweep carry-over across a run of unidentifiable regions.** Proceeding left→right, once a region is
   imputed (by gDNA-fraction *or* absolute density) we have its `U_gDNA`/`U_RNA` split. At the next
   eligible boundary we **carry that fraction forward** to deconvolve the boundary's unspliced, then apply
   the same 3-term `f`. Method-agnostic in how the prior split was obtained — the same forward-state idea
   as the absolute-density bidirectional sweep, carrying a deconvolved fraction.

> **Architectural note.** This lets the count module *consume* a strand-derived split to sharpen its own
> estimate — it does **not** re-introduce the joint strand×count *product* the decoupling removed. There
> the concern was multiplying biased count into unbiased strand; here strand only *cleans the count's
> numerator* (removes nascent), which can only help. It is opt-in (used only when a split exists) and
> reduces to the 2-term form otherwise. Flagged for review as the one place strand feeds count.

### 5.1 Strand-resolved sweep state — do not lump RNA across strands

The carried RNA must be **per strand**. The accumulator is already 4-channel — `unspliced{genome+,
genome−}`, `spliced{sense, antisense}` — and the spliced strand is *reliable everywhere* (motif-oriented,
valid even in AMBIG regions); only the unspliced strand is genome orientation, ambiguous in AMBIG.

The reason it must be per-strand is **identifiability**, not bookkeeping. In an opposite-strand overlap
region (e.g. `E+·I−`), the unspliced reads decompose as
```
n_unspliced_pos = gDNA/2 + nascent⁺ ,   n_unspliced_neg = gDNA/2 + nascent⁻
```
— 2 equations, 3 unknowns: **locally under-determined** (exactly why AMBIG regions are excluded from the
strand module, "D7"). Carrying one strand's RNA forward from a neighbouring single-strand region (where
it *was* identifiable) along that transcript's extent supplies the missing constraint:
`gDNA = 2·(n_unspliced_pos − nascent⁺)`, and `nascent⁻` follows. **Lumping all unspliced RNA into one
bucket discards the per-strand structure and leaves AMBIG overlaps permanently under-determined.** A
scalar carried fraction also cannot represent two independent RNA flows (a `+` and a `−` transcript)
crossing the same overlap region — it conflates them and mis-attributes downstream.

**Design:** the sweep state is the small vector `{gDNA, RNA⁺, RNA⁻}` (mirroring the 4 channels — no new
information, just don't sum away the strand). `RNA⁺`/`RNA⁻` propagate along their respective transcript
extents; at AMBIG overlaps they combine with local `unspliced_pos`/`unspliced_neg` to solve for `gDNA`.
The per-region prior still collapses to the scalar `gDNA / (gDNA + RNA⁺ + RNA⁻)` — strand resolution is
an internal sweep concern, invisible to the prior, and `boundary_gdna_fraction` is unchanged (its
`unspliced_rna` arg receives the strand-resolved RNA total).

**Scope:** bites only for **stranded + nascent + AMBIG overlap**. Single-strand runs have one strand;
pure-unstranded can't separate strands at all (lump, accept nascent-as-gDNA, matched-set keeps the
fraction safe); the zero-nascent gDNA benchmark has empty RNA⁺/⁻ buckets. So build the state per-strand
from the start (cheap, general); the AMBIG carry-resolution activates only when nascent + strand are
present. The eligibility predicate (§3, Tests A/B) is unaffected — it already routes AMBIG *exons* to the
fallback; this is purely the sweep's unspliced deconvolution.

## 6. Fallback — absolute density (always available)

No eligible boundary ⇒ keep the current `density_model` absolute-density imputation (needs only an
observable density, not a junction — always possible). Eligibility is therefore a pure *upgrade*; no
region is left unhandled. **Anchor combination** (one vs two eligible boundaries) reuses the *same*
arithmetic as absolute density — mean / bidirectional sweep / precision weight (Phase 4-var). No new
machinery.

## 7. Unit-test plan — prove eligibility empirically

**(A) Exhaustive predicate test** — `tests/calibration/test_splice_junction_eligibility.py`. Enumerate all
256 `(sig_L, sig_R)`; assert `splice_junction_eligibility` against an independently hand-derived truth
table, plus invariants: strand-swap symmetry, L/R reflection, mixed-state veto (`E_s & I_s` ⇒ `s ∉
junction_strands`), matched-set (`exon_strands == junction_strands`).

**(B) Realized-scenario test** — `tests/calibration/test_splice_junction_realized.py`. For each
*realizable* pair (§4), **construct actual transcripts** (exons/introns on `+`/`−`, overlapping and
antisense in every combination) so the region partition yields adjacent regions with those signatures;
deposit reads of known origin (gDNA, mature⁺/⁻ via spliced, nascent) at controlled ratios; then compute
`f_b` from origin-blind crossing counts + FL models, impute `G_region = f_b · M_region`, compare to
planted truth. Assert: **eligible** pairs recover truth; **ineligible** pairs (antisense junction, AMBIG
exon, retained intron, split sides) *fail* recovery and route to fallback; **unrealizable** pairs are
recorded as un-constructable. This is where eligibility is *proven*, not asserted. Extend with a 3-term
case (§5): strand-clean / carry a split and confirm the pure-gDNA fraction is recovered with nascent
present.

## 8. Implementation sketch (after the predicate is proven)

- `calibration/splice_junction.py` (new): the predicate + `boundary_gdna_fraction(...)` (2-term and
  3-term) + per-region eligible-anchor selection. No new constants.
- `density_model.py`: eligible regions → fraction partition; ineligible → fallback. Clip only after
  locus-level aggregation (deferred-clip rule).
- Confirm the accumulator exposes crossing-unspliced vs crossing-spliced per boundary side (the substrate
  has `n_spliced_*`; add crossing-unspliced if missing — C++ + reference stay byte-identical).

## 9. Open questions / future work

1. **FL-consistency of the partition — `f_b · M_region`.** `f_b` is a *molecular* (density) fraction;
   `M_region` is a *fragment* count. A strictly consistent partition multiplies in the region's own
   `E_region^gDNA / E_region^RNA`. The diagnostic recovered truth without it (ratio ≈ 1 here), so the
   first cut uses the bare form. **Marked as a future investigation** — derive and test the region-eff-
   length correction across region lengths (short exons, where gDNA/RNA FL differences bite). Saved as a
   memory so it is not lost.
2. **Capture-robustness of `f_b`** — holds iff probes enrich exon-body and junction-spanning fragments
   proportionally; the nascent+capture benchmark is the test.

*(Resolved from the prior draft: anchor combination → reuse absolute-density machinery, §6; coverage-by-
mass → moot, these are toy correctness scenarios not real data.)*
