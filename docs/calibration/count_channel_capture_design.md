# Count-channel calibration under hybrid capture — v2 design & implementation plan

**Status:** implementation plan. 2026-06-10. v2 supersedes the v1 root-cause doc and the
strand-cleaning line of inquiry (`strand_cleaning_robustness_design.md`,
`strand_clean_global_target_design.md`, both **deferred** — they chase a minor residual). v2
incorporates the directional decisions from review: locality is in **count space**, not genomic
space; the count-prior mean is imputed only from **observable boundary crossing flux**; dispersion is
**mean/count-dependent**; imputation precision follows **anchor strength**; and the
**unspliced-fraction projection** is promoted from "refinement" to the **primary fix** for the exon
density under-estimate.

---

## 1. Root cause (confirmed)

The count overdispersion fit measures each seed's count against a **single global pooled density**
`ρ̄` (`count_dispersion._nb_mom`, `μ=ρ̄·eff_len`). Under capture the true local gDNA density varies
8–48× (on-target enriched / off-target depleted), so between-region **mean heterogeneity is booked as
overdispersion**: `α_crossing` inflates ~170,000× (0.0005 → 86 off→on; driven by the between-seed
density CV 0.04 → 1.19), and `N_eff=N/(1+αN)` collapses a typical crossing count's effective evidence
from **~515 to 0.012**. The count channel is annihilated; the prior degenerates to Jeffreys Beta(½,½);
the deconvolution leans entirely on strand. When unstranded, both clues are dead — the catastrophic
regime (gdna1000 unstranded+capture 18.7% leak; pure-RNA unstranded+capture FP 4.5%).

Evidence (scripts `scripts/debug/diag_density_locality.py`, `diag_fp_localize.py`; suite
`~/Downloads/rigel_runs/gdna_benchmark_5mb`):
- The 93k false gDNA on gdna_none unstranded+capture is **entirely in exon regions** whose prior mean
  is *exactly 0* and concentration *exactly 0* → Jeffreys → drift to ~0.5.
- The imputed exon prior mean is **0.418 at both ss0.99 and ss0.50** (truth ~0.9) — a
  strand-independent count under-estimate. The leak gap (3.7% stranded vs 18.7% unstranded) is the
  Beta-Binomial strand likelihood rescuing the crushed count prior when it can, and not when it can't.

This splits the leak into two count-side problems:
- **Problem A** — the imputed exon gDNA density is ~2× too low under capture (the dominant accuracy
  bug; §3.2 / §3.3 fix it).
- **Problem B** — the Jeffreys floor discards a correct local mean at low concentration (§3.5 fixes it).

---

## 2. Corrected conceptual model

Four principles, each correcting a v1 assumption:

1. **Locality lives in count space, not genomic space.** Genomic proximity is *not* informative under
   capture — adjacent exon/intron densities vary wildly, so a genomic "neighbourhood mean" is
   implausible. But the **count→variance relationship** is stable: nodes with similar counts share a
   dispersion regardless of where they sit. So the dispersion model is a smooth function of the count
   (DESeq2's mean-variance trend, fit by local regression over the sorted counts), and it
   automatically reduces to the single global value when counts are homogeneous (capture-off) — so
   capture-off behaviour is preserved by construction.

2. **The count-prior mean is imputed only from observable boundary crossing flux.** A non-observable
   region (or run of them) is bordered by **observable boundaries**; the fragments crossing those
   boundaries into the region are the *only* reliable local count signal. Accuracy is probe-placement
   dependent (probes far from a boundary → fewer crossing fragments) but there is no better source —
   the imputation must do its best with the crossing flux. **No genomic sweep / neighbourhood average.**

3. **The mean imputation should use the boundary unspliced *fraction*, projected onto the region's own
   count — not an absolute crossing density magnitude.** The current imputation uses
   `crossing_gdna_count / fl_mean` (a density *magnitude*), which is depleted and noisy at the boundary
   under capture — the source of the 0.418-vs-0.9 under-estimate. Instead use the **ratio**
   `f = unspliced_density / (unspliced_density + spliced_density)` at the observable boundary (a local
   ratio, robust to capture magnitude) and project it onto the region's **own** total count `C` (which
   *is* enriched on-target): `unspliced(gDNA+nascent) ≈ f · C`. Strand then splits gDNA vs nascent
   within. This leverages the region's well-measured enriched count instead of a depleted boundary
   magnitude.

4. **Concentration = baseline count dispersion (observable nodes) + imputation variance (exons).** Two
   distinct uncertainties, conflated today into one global α:
   - *Observable nodes* (count is gDNA by signature): trust ≈ Poisson, limited only by the small
     baseline overdispersion from the count-space trend (§1).
   - *Imputed exon nodes*: trust is governed by **imputation variance**, which follows the **anchor
     strength** — the number of crossing fragments at the observable boundaries used to impute it
     (more crossing fragments → tighter). This is the real reason to distrust an exon, and it is local
     and honest, unlike a global α.

---

## 3. Implementation plan (staged, each independently measurable)

Sequenced to remove catastrophic behaviour first, then fix accuracy, then refine precision. Each stage
has a validation gate on the 20-scenario benchmark before the next.

### Stage 0 — Excise the global-α overdispersion; concentration by node type
**Why first:** the global α is the catastrophe; removing it is the highest-leverage, lowest-risk move,
but it must not swing exons to over-trust their (still-biased) imputed mean.

- `calibrate.py`: stop computing `count_disp = fit_gdna_count_overdispersion(...)` and the
  `region_alpha` blend. Replace `count_evidence` with a node-type rule:
  - **observable region / observable boundary side** → `count_evidence = count_support` (the raw
    effective gDNA count; Poisson trust — these counts *are* gDNA).
  - **imputed exon region / non-observable side** → a deliberately **modest** concentration
    (placeholder until Stage 3), e.g. the anchor's effective crossing count, so the prior is weak but
    not zero. *Not* a global crush.
- Keep `effective_count()` (it stays the home for any α-limit once Stage 4 lands).
- **Gate:** unstranded+capture leaks must drop sharply (the crush is gone); observable-heavy and
  capture-off conditions unchanged; watch for any *new* over-call on stranded+capture exons (would
  indicate the modest exon concentration is too high pre-Stage-2).

### Stage 1 — Fix the degenerate floor (Problem B)
- `joint_deconv._joint_per_node`: when `count_concentration → 0`, the prior must fall back to a point
  near its **own mean** `count_gdna_frac`, not Jeffreys Beta(½,½). Concretely, replace the additive
  Jeffreys floor with a small **minimum concentration** `c₀` applied to the node's *own* mean:
  `a_c = max(conc, c₀)·gf + ε`, `b_c = max(conc, c₀)·(1−gf) + ε` (so at zero count support the prior
  is a weak Beta centred at `gf`, not at ½). `c₀` is the one new constant here — derive it as "≈1
  pseudo-count" and justify against the no-magic-numbers rule.
- **Gate:** gdna_none unstranded+capture FP (−4.5%) → ~0, with no change to informative-count nodes.

### Stage 2 — Unspliced-fraction projection imputation (Problem A, the accuracy fix)
- New imputation path in `density_model.node_gdna_density` (or a new `impute.py`): for each
  non-observable region, gather its anchoring **observable boundary sides** (already identified via
  `boundary_count_observable`). At each anchor compute the boundary unspliced fraction
  `f = mass_unspliced / (mass_unspliced + mass_spliced)` (both available on `SubstrateView`), and the
  region's imputed unspliced gDNA+nascent count `≈ f · C_region` (region's own contained count).
  Average over available anchors (anchor-strength-weighted, see Stage 3).
- The count-prior mean becomes `count_gdna_frac = (strand-cleaned gDNA part of f·C) / mass_unspliced`
  (strand removes nascent from the unspliced part, as today). For observable regions the mean path is
  unchanged (their own count is gDNA).
- Keep the old absolute-density imputation available behind a flag for A/B measurement; default to the
  fraction projection.
- **Gate:** imputed exon prior mean at gdna1000 capture-on rises from 0.418 toward ~0.9 (oracle check);
  flagship stranded+capture leak (4.6%) and unstranded+capture (post-Stage-0) both drop; capture-off
  unchanged.

### Stage 3 — Imputation variance from anchor strength (concentration for exons)
- Replace Stage 0's placeholder exon concentration with one derived from the **anchoring crossing
  count(s)**: the imputed mean `f·C` has precision set by the number of crossing fragments that
  produced `f` (binomial/√N on the ratio) and the consistency across multiple anchors. A region
  anchored by abundant, agreeing boundaries gets a tight prior; one anchored by a single sparse
  boundary stays weak (defers to strand).
- **Gate:** stranded+capture unchanged or better (strand still corrects where imputation is weak);
  unstranded+capture improves further (good imputations now trusted); no FP regression on gdna_none.

### Stage 4 — Count-space dispersion trend (baseline overdispersion, honest)
- New `count_dispersion.fit_dispersion_trend(counts, eff_lens)` → per-node baseline NB `α`, fit by
  **local regression over the sorted counts** (recommend a **Gaussian kernel in log-count space** —
  weights near neighbours more, no hard bin edges; the one tunable is the bandwidth). Apply via the
  existing `effective_count(count, α_node)` to observable nodes' concentrations (Stage 0 set these to
  raw counts; this caps them at the honest baseline). Capture-off → homogeneous counts → recovers the
  single global α (continuity with today's good non-capture behaviour).
- **Gate:** capture-off conditions still ≤1.5% (no regression from re-introducing a *small* baseline
  dispersion); capture-on unaffected (the trend is small there because it's count-local, not global).

### Stage 5 — Retire the contained/crossing type split (point 6)
- Once dispersion is count-local (Stage 4) and concentration is node-type + anchor-strength based
  (Stages 0/3), the `alpha_contained` / `alpha_crossing` distinction has no remaining role. Remove the
  split and the two-type plumbing in `count_dispersion.py` and `calibrate.py`. Confirm no behaviour
  change.

---

## 4. Data-flow summary (after all stages)

```
substrate (per-node unspliced/spliced counts + mass, from accumulator)
  → rna_sense_frac κ                                            (strand_balance)
  → per-node count-prior MEAN:
        observable nodes:  own strand-cleaned gDNA count / mass
        imputed exons:     [ f·C projection from observable anchors ] strand-cleaned / mass   (§3.2)
  → per-node count-prior CONCENTRATION:
        observable nodes:  count capped by count-space dispersion trend α(N)                  (§3.4)
        imputed exons:     anchor-strength imputation precision                               (§3.3)
  → joint deconv: Beta(mean, conc) [floor → own mean, not ½]  ×  Beta-Binomial strand         (§3.1/§3.5)
  → derive (gdna_density_global QC, geometric gDNA length)
```

No global density target anywhere; the only cross-node fit is the **count-space** dispersion trend.

---

## 5. Open decisions (no-magic-numbers discipline — settle before/within each stage)

1. **Stage 1 floor `c₀`** — minimum concentration / pseudo-count. Derive as "1 effective observation";
   confirm it only bites at near-zero support.
2. **Stage 2 anchor combination** — when a region has multiple observable anchors (left/right, or a run
   bordered on both ends), how to combine their `f` (anchor-strength-weighted mean? nearest?). And how
   `f·C` interacts with runs of consecutive non-observable regions (carry the bordering anchors inward).
3. **Stage 3 precision form** — the functional map from crossing fragment count → imputation precision
   (√N on the ratio `f`? Beta-binomial on `f`? plus a between-anchor consistency term).
4. **Stage 4 bandwidth** — the kernel bandwidth in log-count space (data-driven, e.g. a fixed fraction
   of the log-count IQR, or Silverman). The single tunable of the dispersion model; justify it.
5. Whether Stage 2's fraction projection fully replaces, or blends with, the absolute-density
   imputation (keep both behind a flag through Stage 2's gate, then decide).

## 6. Risks

- **Stage 0 over-trust:** removing α before Stage 2 lands could let a biased exon mean over-assert.
  Mitigated by the modest exon concentration placeholder and by sequencing Stage 2 immediately after.
- **Fraction projection assumes the boundary `f` transfers to the region interior.** Probe placement
  can break this (probe mid-exon → boundary `f` unrepresentative). Stage 3's anchor-strength precision
  is the safety valve (weak anchor → weak prior → strand governs). Quantify on the benchmark.
- **Nascent (untested by `nrna_none`).** The fraction projection lumps gDNA+nascent as "unspliced";
  strand splits them downstream, as today. Needs a nascent-bearing sim to validate — flag for the
  simulator backlog, do not assume.

## 7. Validation harness

Skill `calibration-benchmark` (net fragment-flow) on all 20 scenarios at each gate. Track: per-condition
net leak; imputed exon prior mean vs oracle (Problem A); gdna_none FP (Problem B); capture-off
no-regression. Unit tests per stage (floor → own mean; `N_eff→N` as α→0; dispersion trend → global
under homogeneous counts; fraction projection recovers known gDNA fraction on a synthetic node).
