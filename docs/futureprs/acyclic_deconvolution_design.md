# Acyclic Deconvolution — Calibration Redesign

**Supersedes** the Stage-1 proposal in `rho0_stabilization_design.md` and demotes the
exposure-prior line of work (`exposure_prior_findings_report.md`,
`exposure_prior_noise_floor_theory.md`) to optional. Those remain the diagnostic record
that motivated this redesign. This revision incorporates an external review (Risks A–D
below) and the design decisions: **single-pass, fixed-seed; every region/boundary gets a
decoded result; unbiased `z=0` default**.

---

## 0. Thesis

The sparse-data instability is a **feedback cycle**, not a prior-tuning problem: the E-step
decides each region's gDNA from the *modeled* global density `μ_g = ω·ρ₀·L`, then
re-estimates the density from those decisions. A closed loop cannot be stabilized by a
prior on the loop — only by opening it.

**Redesign:** deconvolve each region's gDNA from **local evidence only** (a strand-split
posterior + boundary-crossing density), then **derive** `ρ₀` and `ω` as *outputs*. The
graph becomes acyclic and the pipeline becomes **single-pass**. This:

- eliminates the runaway collapse (sparse → self-averaging sampling noise, not a runaway);
- **dissolves the exposure-dispersion prior** (`φ`, noise floor, `κ` knob) — `ω` is a
  derived ratio, not a brittle latent;
- uses the **captured exons** (highest-yield gDNA in a capture panel) instead of excluding
  them;
- exposes one **continuous confidence** knob `z` (the RNA specificity/sensitivity control);
- degrades gracefully as strand information weakens.

Two refinements keep it robust at the edges (from review): the local split is a
**strand-regularized posterior**, not a raw MLE (so thin counts don't over-commit); and the
**derived `ω` is shrunk toward 1 by local information scarcity** before export (so a
zero-evidence region cannot collapse `gdna_eff_len` and re-create Mechanism A downstream).

It removes machinery (the count-channel feedback and the dispersion prior) rather than
adding it.

---

## 1. The disease: the feedback cycle

```
ρ₀ ──► μ_g = ω·ρ₀·L ──► count-channel gDNA call ──► per-region m_g ──► ρ₀ (re-estimate)
└────────────────────────────── feedback ───────────────────────────────────┘
```

`nrna_dc`, n=2000: a small downward wobble in `ρ₀` lowers `μ_g` everywhere → the count
channel expects ≈0 gDNA → calls ≈0 → `ρ₀` collapses to 5e-5 (~800× low); mass leaks to
RNA/nRNA. With ω forced to 1 the same loop runs; with ω free the opposite edge over-calls
gDNA. A prior on `φ` cannot escape a closed loop from either side.

## 2. The cure: an acyclic, single-pass graph

```
seed ρ₀ (strand-free) ─► local deconvolution (strand-posterior + crossing) ─► Ĝ_r per node
                         (NO global-ρ₀ input)                                     │
   ┌──────────────────────────────────────────────────────────────────────────┘
   ▼
 impute non-decodable nodes (crossing → baseline fallback)  ─►  ρ₀ = ΣĜ/ΣL  ─►  ω_r (shrunk)
 (every node gets a result)                                     (derived outputs, never fed back)
```

`Ĝ_r` is estimated from node `r`'s own counts. `ρ₀`, `ω` are summaries computed **once**;
nothing downstream re-enters the deconvolution. A noisy `Ĝ_r` perturbs the aggregate `ρ₀`
by a self-averaging, *unbiased* amount — imprecise, never catastrophic. Because the seed is
fixed and the imputation is one-shot, this is a feed-forward pipeline (no EM loop).

## 3. Local deconvolution as a strand-regularized posterior *(addresses Risk A)*

Per region, unspliced fragments split sense `S` / antisense `A` (`N = S+A`; spliced
fragments are deterministic RNA, PR7). RNA carries the library sense rate `κ_rna` (from the
clean spliced channel, `fit_strand_balance`); gDNA is unstranded (κ_d = 0.5). The
**point identity** is

```
R = (S − A) / (2κ_rna − 1),     G = N − R         (κ_rna = 1, 150/50 ⇒ R=100, G=100)
```

but the shipped estimator is the **Beta-Binomial posterior** on the gDNA fraction `π_g`,
using the library strand constants (`κ_rna, ρ_r_bb` for RNA; `0.5, ρ_d_bb` for gDNA) as a
two-component likelihood with a weak prior on `π_g`. The closed form above is its
high-count mean; at low counts the posterior is **wide** and `π_g` regresses to an
uninformative split rather than over-committing ("2 antisense by chance ⇒ 100% gDNA" is
impossible — the posterior stays diffuse). Crucially, this regularizer is the **strand
constants, not `ρ₀`**, so it adds no feedback edge — acyclicity is preserved.

The posterior yields `Ĝ_r` (mean) **and its uncertainty**, both of which the rest of the
design consumes (`z` in §5, the shrinkage in §4).

It does the right thing on the cases that broke the old design:
- **Captured *unexpressed* exon** — pure gDNA ⇒ ~50/50 ⇒ `Ĝ≈N`. Highest-yield gDNA now
  flows into `ρ₀`.
- **Intronic nascent RNA** — sense-stranded ⇒ `R`, not gDNA. Strand separates nascent RNA
  from gDNA cleanly.
- **Expressed exon** — sense-heavy ⇒ small `Ĝ`, large `R`. Correct, locally, no `ρ₀` input.

## 4. Derived density and the *protected* exposure *(addresses Risk C)*

```
ρ₀  = (Σ_decodable Ĝ_r) / (Σ_decodable L_r)            # library-wide gDNA density
ω_r = shrink( (Ĝ_r/L_r) / ρ₀ ,  toward 1,  by info_r ) # exposure, protected before export
```

The downstream prior needs only the per-node gDNA mass (`alpha_gdna_add ∝ ρ₀·ω·L = Ĝ`);
`ρ₀`/`ω` are reparametrizations for QC and the enrichment report. There is **no Gamma
prior, no dispersion fit, no EM shrinkage of a latent**.

**But the derived `ω` cannot be exported raw.** A zero-/low-evidence region gives
`Ĝ_r ≈ 0 ⇒ ω_r ≈ 0 ⇒ gdna_eff_len → 0`, which makes the gDNA component infinitely
attractive in the downstream EM — Mechanism A, re-entering by a new path. So `ω_r` is
**shrunk toward 1** with strength set by the local information `info_r` (the gDNA
posterior's precision / effective count), via a **unit-information** rule
`ω_shrunk = (info_r·ω_r + 1·1)/(info_r + 1)`: zero evidence ⇒ `ω→1` (baseline density, safe,
no collapse); abundant evidence ⇒ `ω→` the estimate (capture enrichment preserved, because
captured = high count = high `info_r` = little shrinkage). This is a one-shot, acyclic
shrinkage of a *derived* output (not a fed-back latent), parameter-free (unit prior), and
continuous (no floor cliff). It is the structural guard that the old `gdna_eff_len`-floor
hack was groping toward.

## 5. The confidence knob `z` *(addresses Risk D; default per decision #2)*

`z` is a **continuous quantile** of the per-region gDNA posterior (§3), not a threshold:

- `z = 0` → posterior mean → unbiased `ρ₀` → balanced calls. **Shipped default.**
- `z > 0` → upper-quantile gDNA (fewer, purer RNA calls); intentionally over-estimated gDNA
  → larger downstream gDNA prior → more fragments removed from RNA. **RNA specificity.**
- `z < 0` → toward RNA sensitivity (accept gDNA contamination in the RNA pool).

One knob, propagating coherently from the local split through (intentionally biased) `ρ₀`
to the downstream prior. Because it shifts a posterior quantile, it is smooth and "honest"
— no masked cutoff. (The over-estimate at `z>0` is *intended*, exactly the requested
specificity behavior.)

## 6. Decodability hierarchy + the every-node guarantee *(addresses Risk B; decision #1)*

A node is decoded by the first applicable tier; **every region and boundary must receive a
result** (decision #1), so the hierarchy ends in a guaranteed fallback:

| Tier | Signal | Applies to | Needs strand? |
|---|---|---|---|
| 1 | **Annotation** | intergenic (TS_NONE → all gDNA) | no |
| 2 | **Strand-posterior** (§3) | POS/NEG exons + introns | **yes** |
| 3 | **Boundary-crossing density** (§8) | AMBIG, weak-strand, captured "splash" | no |
| 4 | **Expressed/unexpressed latent** (§7) | exons in *unstranded* libraries | no |
| F | **Baseline fallback** | isolated zero-evidence nodes | no |

Tier F sets `Ĝ_r = ρ₀·L_r` (ω=1). It is **neutral for `ρ₀`** (it contributes `ρ₀·L` over
`L`, so re-deriving `ρ₀` including imputed nodes does not move it — no circularity) and safe
downstream (ω=1, no collapse). So an isolated, dataless region degrades to "assume library
baseline," never to a failure.

- **Strand-specific library:** tiers 1+2+3 cover almost everything; `ρ₀` robust; capture
  works via tier 2 without leaning on crossings.
- **Unstranded library:** tier 2 unavailable; fall back to 1+3+4(+F), reduced accuracy (§7).

## 7. Unstranded / weak-strand mode

**The concession.** Without strand, an unspliced fragment may be gDNA *or* nascent RNA —
**not identifiable**. We accept this: mature RNA stays identifiable (spliced junctions),
nascent-RNA recovery degrades. A property of the data, stated openly.

**Seed + tier 4.** Intergenic/intronic + crossings seed `ρ₀` strand-free. Then (single
pass, **against the fixed seed** per decision #1) classify each exon by whether its total
density is statistically consistent with the seed `ρ₀`: density ≈ `ρ₀` ⇒ unexpressed
(≈pure gDNA) ⇒ absorbed into the decodable pool (bolstering `ρ₀` where intergenic/intronic
seed regions are sparse); density ≫ `ρ₀` ⇒ expressed ⇒ gDNA imputed via tier 3/F. Because
classification is against the **fixed seed**, tier 4 stays acyclic (no `ρ₀`→classify→`ρ₀`
loop). Iterating it is deferred until single-pass stability is demonstrated.

**Unstranded capture is partially recoverable** (decision #3, correcting the earlier
doc). A captured region — even unexpressed, even unstranded — receives a large
**boundary-crossing "splash"** from probes (exon-intron crossing fragments), which tier 3
imputes as gDNA. That lifts its `Ĝ_r` well above the intergenic/intronic baseline, so its
derived exposure is `ω ≫ 1`. Accuracy is lower than the stranded path (the splash is a
proxy for the in-probe peak, and nascent RNA contaminates it), but capture enrichment is
**not** lost — it is partially decoded from crossings.

## 8. The boundary-crossing density signal

Boundary-crossing fragments carry a **local, ρ₀-free** gDNA density:
- gDNA spans boundaries freely (contiguous genomic DNA);
- spliced mature RNA crosses only at annotated junctions;
- nascent RNA crosses intron/exon boundaries (unspliced), like gDNA.

So non-junction crossings are gDNA-or-nascent; with strand they deconvolve (gDNA
unstranded, nascent sense), without strand they are *assumed gDNA* (the concession). This is
the existing D1/sweep machinery, **repurposed to impute `Ĝ_r`** (not exposure). It makes
AMBIG and weak-strand nodes decodable (tier 3), reads the capture "splash" across
captured-exon boundaries (§7), and is the reason intergenic depletion in capture panels does
not blind the seed. Its imputed `Ĝ_r` **can exceed baseline** (captures enrichment).

## 9. What this replaces / simplifies

- **Delete** the ρ₀-dependent count channel as feedback (`estep._llr_count` over
  `μ_g = ω·ρ₀·L`). Deconvolution = strand-posterior + crossing.
- **Dissolve** the exposure-dispersion prior: the Gamma posterior, `update_exposure_dispersion`,
  `exposure_prior_pseudocount`, the noise-floor theory. `ω` becomes §4 (a derived, shrunk
  ratio). Net: fewer modules and constants.
- **Keep & repurpose** boundary/sweep (now imputes `Ĝ_r`), `fit_strand_balance` (`κ_rna`,
  dispersions — now the local-posterior regularizer), the density seed (tier-1/tier-4
  anchor), the decodable-node masks (extended into the §6 hierarchy).

## 10. Risks (review) and how the design answers them

- **A — low-count/unstranded variance trap.** Answered by §3: the split is a strand-
  regularized *posterior*, not an MLE; thin counts stay diffuse. (Regularizer = strand
  constants, not `ρ₀` → no feedback.)
- **B — AMBIG breakdown.** Answered by §6 tiers 3+F: crossing-density imputation, then a
  guaranteed neutral baseline fallback; every node decoded.
- **C — derived-ω collapse → Mechanism A redux.** Answered by §4: unit-information shrinkage
  of derived `ω` toward 1, protecting `gdna_eff_len`. (This is the single most important
  guard; validate it explicitly.)
- **D — `z` as a cliff.** Answered by §5: `z` is a continuous posterior quantile, default 0.

**Remaining open items (validate, don't assume):**
1. **Tier-4 single-pass stability** (the one residual feedback-shaped step, even fixed-seed):
   stress-test on an unstranded toy that tries to make it drift.
2. **Joint coverage**: enumerate region × data-regime; confirm tiers 1–4+F leave no node
   undetermined and none mis-tiered.
3. **ρ₀/ω scale split under capture**: pooling `Ĝ` from enriched exons constrains `ρ₀·ω`;
   confirm `E[ω]=1`-style normalization recovers the global scale without smuggling
   enrichment into `ρ₀`. (Open question: is the normalization explicit, or implicit in
   "ρ₀ = ΣĜ/ΣL over all decodable nodes including enriched ones"? The latter makes `ρ₀` a
   coverage-weighted average density, which is well-defined but means `ω` is relative to
   that average — acceptable, but document it.)
4. **Low-`κ_rna` continuity**: `1/(2κ_rna−1)` must degrade smoothly into the tier-3/4 path,
   no discontinuity at the stranded↔unstranded handoff.

## 11. Validation plan

1. **Circularity is gone (headline).** Strand-first deconvolution; confirm the `nrna_dc`
   collapse simply does not occur (ρ₀ stable, gDNA ≈ exact) with **no** ρ₀ prior, **no** φ
   machinery.
2. **Risk-C guard.** A region driven to zero local evidence must yield `ω→1`, not `ω→0`;
   confirm `gdna_eff_len` never collapses.
3. **Capture, stranded.** Add **silent-but-captured** genes (pure enriched gDNA — the case
   the old design under-called 689/2740); confirm recovery via tier 2.
4. **Capture, unstranded.** Same genes, `κ_rna→0.5`; confirm *partial* recovery via tier-3
   splash (`ω>1`, accuracy reduced but not lost) — decision #3.
5. **Three-regime harness** passes with the prior machinery removed.
6. **`z` sweep**: monotone, smooth RNA-purity / gDNA-recall tradeoff.
7. **Scenario + golden suites** regenerate (global behavior change).

## 12. Implementation sketch

- `estep.py` — replace `logit(prior)+llr_count+llr_strand` with the §3 strand-posterior (+
  §5 quantile); remove `_llr_count`/`μ_g`.
- `calibrate.py` — single pass: seed → local `Ĝ_r` → impute (tiers 3/4/F) → derive `ρ₀`
  (§4) → derive+shrink `ω` (§4). Remove the exposure posterior, `update_exposure_dispersion`,
  the outer EM loop.
- `exposure.py` — reduce to the §4 derived+shrunk ratio (or fold into `calibrate.py`).
- `mstep.py` — `update_rho_0` becomes the §4 aggregate over decodable `Ĝ_r`; dispersion fit
  removed.
- `density.py` — seed stays as tier-1 anchor + tier-4 reference.
- `config.py` — remove `exposure_prior_pseudocount`; add continuous `z` (default 0).
- Boundary/sweep repurposed to impute `Ĝ_r` with the tier-3 → tier-F fallback.

End state: a feed-forward calibration — strand-posterior + crossing → `Ĝ_r` → `ρ₀`,`ω` →
prior — with no EM feedback and a single optional tier-4 refinement for unstranded
libraries.

---

## 13. First exploration step (next)

Per "explore then implement," the first prototype is the **headline test**: a throwaway
script computing the §3 strand-posterior `Ĝ_r` on `nrna_dc` + a silent-captured-exon
scenario, deriving `ρ₀` strand-only, to **show the collapse does not occur** and the
captured gDNA is recovered — *before* touching production modules. Then validate the Risk-C
guard (drive a region to zero evidence, confirm `ω→1`).
