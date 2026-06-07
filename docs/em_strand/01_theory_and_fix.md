# Strand likelihood in the locus EM — theory, diagnosis, and fix

**Status:** theory + diagnostics complete (2026-06-06; §4.1–4.3 supersede an earlier reading
that was wrong on two points — see those sections). Key results: (1) the per-fragment strand
factor **is load-bearing** — it catches the antisense half of in-locus gDNA; neutralizing it
alone makes the leak *worse* (§4.2). (2) The antisense gDNA **is already inside the locus
EM** (the "anchor"), so the fix needs **no locus reconstruction** — just an aggregate
count×strand split term (B) + per-fragment neutralization for the split (C). Fix not yet
implemented.
**Scope:** the per-fragment strand factor in the post-calibration locus EM
(`scoring.cpp` + `em_solver.cpp`), and how the gDNA-vs-RNA strand discrimination
should be partitioned between the per-fragment E-step and the aggregate (region /
locus-weight) level.

This note is the durable record for the `em_strand` line of work. Diagnostic script:
`scripts/debug/em_strand_leak_stratify.py`.

---

## 1. The current formulation

Strand enters Rigel at **two altitudes**.

**Calibration (region scale, aggregate — correct).**
`calibration/strand_likelihood.py` + `joint_deconv.py` model a region's
sense/antisense *counts* as a **Beta-Binomial mixture**: gDNA contributes a
Beta-Binomial with mean `0.5` and overdispersion `rho_d`, RNA a Beta-Binomial with
mean `rna_sense_frac` and `rho_r`. The per-region gDNA fraction is the posterior of a
`count × strand` joint. This is overdispersion-aware and operates on the count
sufficient statistic. See `docs/caljointmodel/01_generative_model.md` §4.3.

**Locus EM (per fragment — the problem).**
`scoring.cpp::score_mm_alignment` emits, for every fragment:
- gDNA hit (`scoring.cpp:516`): `... + LOG_HALF + ...` — a flat `log 0.5` per
  fragment, unconditional.
- RNA candidate (`scoring.cpp:553-563`): `log(p_sense)` if the fragment orientation
  matches the transcript strand, `log(p_antisense)` if opposite, `LOG_HALF` if the
  transcript is strandless.

Those per-fragment numbers are the **only** strand signal the locus EM sees. The
calibration prior that crosses in (`priors.py` → `AggregatePrior` in `em_solver.cpp`)
is two strand-free scalars (`alpha_gdna_add`, `alpha_rna_add`). The M-step
(`em_solver.cpp::apply_grouped_prior_update`, line 710) forces the split as
`G = n_gdna + alpha_gdna_add`, `R = n_rna + alpha_rna_add` — but `n_gdna`/`n_rna` are
the **E-step raw counts**, which are themselves driven by the per-fragment strand
log-liks. So the gDNA-vs-RNA split is litigated in *both* places.

---

## 2. Why the per-fragment `log 0.5` misbehaves

`log 0.5` is **not** an illegitimate penalty in isolation — it is the honest emission
probability of a Bernoulli(½) outcome (gDNA is double-stranded; either orientation is
equally likely), and it is the *same* `0.5` for both orientations, so it cannot bias
*which* strand gDNA prefers. The unstranded limit is harmless: when `p_sense = 0.5`,
the RNA term is also `log 0.5` and **cancels** in the responsibility ratio.

The defect is what per-fragment `log 0.5` **secretly assumes**. For a hypothesis-fixed
dataset, summing per-fragment Bernoulli log-probs equals the *aggregate* log-likelihood
up to a hypothesis-independent constant:

```
sum over gDNA fragments of log 0.5  =  N·log 0.5  =  Binomial(N, k; ½) log-likelihood + const
```

That is the **rho_d = 0 (zero-overdispersion) special case** of the correct
Beta-Binomial. It asserts gDNA's sense/antisense split is a *rigid* 50/50 with no
sampling slack. Real gDNA is overdispersed (capture bias, PCR, pair-anchoring noise) —
which is exactly why `rho_d` exists. So:

> The per-fragment `log 0.5` gDNA strand factor is an **overconfident** gDNA strand
> model. Its overconfidence manifests as an effective penalty against gDNA whenever the
> local strand split deviates from 50/50 — and it grows with library strandedness
> (`p_sense → 1`), vanishing at unstranded (`p_sense → 0.5`).

Worked example (the doc's own 126/26 gDNA pile, `p_sense = 0.95`, + transcript):
- gDNA, per-fragment: `152 · log 0.5 = −105.3`
- RNA, per-fragment: `126·log 0.95 + 26·log 0.05 = −6.5 − 77.9 = −84.4`

RNA wins by 21 nats — a genuine, noise-asymmetric gDNA pile is misclassified as RNA,
because the rigid Binomial(½) gDNA model punishes *any* asymmetry while the RNA model
flexes toward the majority strand. The aggregate Beta-Binomial(½, rho_d) keeps it gDNA.

---

## 3. E-step or M-step? Strand balance is intrinsically aggregate

A per-fragment E-step factor is valid only if fragments are conditionally independent
given their component. With overdispersion they are **not**: all gDNA fragments in a
region share a latent per-region strand rate `kappa^(g) ~ Beta(alpha_d, alpha_d)`. Once
fragments share a latent rate, the joint likelihood does **not** factor into independent
per-fragment terms — you must work with the count sufficient statistic `(k_sense, N)`.
This is why overdispersion "only makes sense at the region scale."

Therefore the two jobs strand does must live at different altitudes:

- **gDNA-vs-RNA balance** ("is this pile 50/50 or skewed?") — an aggregate,
  overdispersed-count quantity. Belongs where the count sufficient statistic lives:
  the calibration region deconvolution (already done), or an aggregate
  **weight-level / M-step** term in the locus EM. Putting it per-fragment forces
  `rho_d = 0` and double-counts the calibration evidence.
- **Within-RNA transcript identity** ("did this fragment come from the + gene or the −
  gene overlapping here?") — a genuine per-fragment question. `log p_sense` vs
  `log p_antisense` legitimately belongs in the **E-step**, with no overdispersion
  issue, because it identifies a source rather than asserting a balance.

The current code conflates these. The fix is to separate them.

### Historical note — the two-component gDNA model

An earlier design split gDNA into gDNA⁺ / gDNA⁻ components with a free mixing weight.
That "worked" because it let the EM allocate gDNA asymmetrically at no penalty — i.e. it
injected gDNA strand overdispersion by brute force. It was overkill because a single
gDNA component with a **Beta-Binomial(½, rho_d)** strand model captures the same
flexibility with one scalar. Collapsing the two components was the right *computational*
call but silently reinstated the rigid-50/50 per-fragment penalty, because only
calibration picked up the Beta-Binomial replacement — the locus EM did not.

---

## 4. Empirical confirmation (2026-06-06)

Diagnostic `scripts/debug/em_strand_leak_stratify.py` on the synthetic hybrid-capture
suite (`hyb_capture_500kb`). For every gDNA fragment that **genomically overlaps a
transcript** (truth-based, not the pipeline's `ZL` — see note below), classify it as
**co-stranded** (its genomic orientation matches the overlapping transcript strand) or
**antisense**, and measure the gDNA→RNA leak rate (called RNA via annotated-BAM `ZF`
bit `0x04`) in each bucket.

| condition | co-stranded leak | antisense leak | asymmetry |
|---|---|---|---|
| ss99, capture-off | **15.8%** (n=18510) | **1.4%** (n=18432) | **+14.4 pts** |
| ss99, capture-on  | **37.4%** (n=49256) | **3.3%** (n=48924) | **+34.0 pts** |
| ss50, capture-off | 12.1% (n=18203) | 12.4% (n=18650) | −0.2 pts |
| ss50, capture-on  | 40.7% (n=49276) | 40.6% (n=48853) | +0.1 pts |

The co/anti fragment counts are ~50/50 — the gDNA itself is symmetric. Findings:

1. **The predicted signature is present and large.** At ss99 the co-stranded leak
   exceeds the antisense leak by **+14 pts (capture-off)** and **+34 pts
   (capture-on)**. The per-fragment `log p_sense` reward steals co-stranded gDNA into
   RNA; `log p_antisense` catches antisense gDNA.
2. **It vanishes when unstranded.** At ss50 the asymmetry is ≈0 (−0.2 / +0.1 pts) — the
   strand factor cancels, exactly as theory says.
3. **It scales with capture** — capture co-locates more gDNA on expressed exons, giving
   the strand factor more co-stranded gDNA to steal.

> Note on the denominator. The table above is from a pre-2026-06-06 build (the first
> measurement); current `main` leaks roughly half as much (the FL fix, main@36636f4). The
> *signature* (co ≫ anti at ss99, symmetric at ss50) is unchanged. **Also: do not use the
> pipeline's `ZL` tag to judge "in the EM."** `ZL` is the *winning-transcript* locus, so a
> fragment that wins the **gDNA** component gets `ZL = −1` even though it was fully inside a
> locus EM. An early cut using `ZL` made it look like antisense gDNA was "routed to
> intergenic"; the `ZF` pool-flag decode (§4.1–4.3) shows that was an artifact.

### 4.1 The signature is the per-fragment strand penalty (corrected)

`scripts/debug/em_strand_fate.py` decodes each gDNA fragment's `ZF` outcome by
location × strand. For **exon-overlapping** gDNA (the leak-relevant set), current-build
baseline:

| condition | EXONIC co (em-gDNA / mRNA-leak) | EXONIC anti (em-gDNA / mRNA-leak) |
|---|---|---|
| ss99 cap-on | 80% / **20%** | 99% / **1%** |
| ss50 cap-on | 80% / 20% | 80% / 20% |

Both co- and antisense exon gDNA are **`em-gDNA(0x05)` — fully inside the locus EM** (not
intergenic). The asymmetry is the **per-fragment strand factor** (confirmed by the ablation
in §4.2), acting in two opposite directions:

- **Antisense gDNA is caught** because the mRNA hypothesis for an antisense fragment gets
  `log p_antisense` (≈ −4.6 at ss99) → mRNA loses → gDNA wins (1% leak).
- **Co-stranded gDNA leaks** because the gDNA hypothesis gets `log 0.5` while the mRNA
  hypothesis gets `log p_sense` (≈ 0) → gDNA loses to mRNA (20% leak). This *is* the ρ_d=0
  penalty of §2.

At ss50 (`κ=½`) the factor cancels and both leak symmetrically (20%).

### 4.2 Ablation (corrected — both scoring paths): per-fragment strand is load-bearing

The first ablation patched only the multimapper path; **most gDNA fragments are unique
mappers** scored on a separate path (`scoring.cpp:983`) that still used `LOG_HALF`, so it
was a near-no-op. With **both** paths neutralized:

| condition (cap-on) | baseline co / anti | ablation co / anti | total leak base → abl |
|---|---|---|---|
| ss99 | 20% / 1% | **14% / 14%** | 9.9% → **13.4% (worse)** |
| ss50 | 20% / 20% | 20% / 20% | 19.1% → 19.4% (≈) |

Neutralizing **symmetrizes** the two halves: co improves (20→14%, removing the ρ_d=0
penalty) but antisense regresses badly (1→14%, losing the catch), and the regression wins
→ total leak up. So:

1. **The per-fragment strand factor is load-bearing** — it catches the antisense half. The
   original "C-alone makes it worse" prediction was correct; my interim "no-op / not the
   lever" reading was an artifact of the half-applied patch.
2. **C (neutralization) alone is not a fix** — it must be paired with an aggregate
   replacement (B) that catches *both* halves.

### 4.3 The anchor IS in the locus EM — piece A is unnecessary (corrected)

The decisive correction: §4.1 shows antisense exon gDNA are `em-gDNA(0x05)`, i.e. **inside
the locus EM**, attached to the same transcript's locus as the co-stranded fragments (the
antisense fragment's only candidate is that transcript; it is not dropped at resolution).
The earlier "anchor routed to intergenic, bring it into the locus (piece A)" conclusion was
based on the `ZL` misread (above) and is **withdrawn**.

So the locus EM already has, per single-strand locus, both `N_sense` (mRNA + co-stranded
gDNA) and `N_anti` (≈ pure gDNA) unspliced fragments. The antisense count is a strong
anchor (≈ the co-stranded gDNA count, since gDNA is 50/50). An **aggregate count×strand
term on the locus's own `(N_sense, N_anti)` will work** — no locus reconstruction needed.

**Budget vs assignment** (`em_strand_locus_diag.py`): `gdna_em ≈ gdna_prior` per locus
(locus 1: 31730 ≈ 29778; locus 5: 15918 ≈ 14689). The gDNA **budget total is ≈ spent** —
this is an *assignment* problem (the budget is spent on off-exon gDNA while co-stranded
exon gDNA is labeled mRNA), consistent with calibration being ≈ correct (the ~2–3%
undercount is a separate, minor issue).

> **Separate finding (the gating, confirmed).** `scoring.cpp:962` / `:980` emit the gDNA
> log-lik **only when a transcript candidate survives** (`new_cands > 0`, `best_t ≥ 0`). A
> gDNA fragment with no surviving candidate — **intronic** gDNA, or one whose only candidate
> is pruned — loses its gDNA hypothesis and is held out of the EM. This does **not** affect
> the exon-overlapping leak (those have a candidate), but it is a real gap worth a separate
> look: intronic gDNA within a gene span contributes nothing to the locus gDNA component.

---

## 5. The fix (architecture) — revised by the ablation

> **Superseded by [02_mstep_implementation_plan.md](02_mstep_implementation_plan.md) v3.**
> Re-derivation from the generative model replaced the "B (bolt-on aggregate term) + C
> (neutralize)" framing below with a clean latent-strand-rate model (no neutralization, no
> multi-strand special case, closed-form native update) — **and** flagged that the premise
> may be wrong: the leak is measured on hard argmax labels while abundances are soft, the
> per-locus budget is ≈ right, and a Beta(½) rate does not anchor sense-to-antisense count.
> A **gating measurement** (does the leak bias soft abundances, or is it a hard-label
> artifact?) must run before any EM change. The B/C sketch below is kept for history.

The (superseded) fix sketch was **B + C**, both at the locus level (no piece A — the anchor
is already in the locus, §4.3):

**B — aggregate count×strand split (the mechanism).** At a single-strand locus, the EM
already holds the unspliced fragments split into `N_sense` (mRNA + co-stranded gDNA) and
`N_anti` (≈ pure gDNA). gDNA is symmetric (mean ½), so the antisense count directly
estimates the co-stranded gDNA: in the `κ≈1` limit, gDNA sense ≈ `N_anti`, so mRNA ≈
`N_sense − N_anti`. Mechanically, replace the M-step gDNA-vs-RNA split
(`em_solver.cpp::apply_grouped_prior_update`, currently `G = n_gdna + alpha_gdna_add`) with
the **count × strand joint** — identical in form to `calibration/joint_deconv` — combining:
- the EM's current count evidence (`n_gdna`, `n_rna` responsibilities),
- the calibration prior (`alpha_gdna_add`/`alpha_rna_add`) as the Beta prior on `π_g`,
- `strand_loglik(π_g, N_sense, N_anti, κ_σ, ρ_d)` — the aggregate strand evidence.

Because this is the *aggregate* count model (not the per-fragment factor), it uses the
antisense fragments to infer the co-stranded gDNA **without** the ρ_d=0 penalty, catching
*both* halves. `strand_loglik`, `κ_σ`, `ρ_d` are existing calibration outputs — **no new
constants**.

**C — neutralize the per-fragment strand factor for the split** (the `strand_neutral_split`
flag, already built — but it must patch *both* scoring paths, see §4.2). With B owning the
gDNA-vs-RNA strand evidence, the per-fragment factor must stop doing so to avoid
double-counting; it still ranks +/− isoforms. C alone is net-negative (§4.2); it is only
valid *together with* B.

**Why no piece A.** §4.3: antisense gDNA is `em-gDNA(0x05)` — already in the locus EM. The
earlier plan to "route off-strand gDNA into the locus" is withdrawn; the locus already has
`N_anti`.

**Multi-strand loci** (both `+` and `−` transcripts overlap): skip the term (count-only),
exactly as calibration's AMBIG handling.

---

## 6. Next experiments

### 6.1 Build B + C and validate

Implement the count×strand joint in the M-step (reuse `strand_loglik`); enable
`strand_neutral_split` (both paths). Re-run `em_strand_leak_stratify.py` +
`em_strand_fate.py` + `bench_calibration.py`. Targets:
- **co-stranded leak ↓** (toward the antisense level) at ss99, **without** regressing the
  antisense catch — the B+C combination must beat *both* baseline (9.9%) and the C-only
  ablation (13.4%);
- **ss50 unchanged** (term vanishes when `κ=½`);
- **gdna_none FP unchanged** (no new false gDNA);
- pool gDNA fraction stays within ±2.2%.

### 6.2 Isoform-discrimination guard

A scenario with overlapping opposite-strand transcripts; confirm B+C preserves +/−
assignment (the per-fragment +/− gap is retained by construction; verify empirically).

### 6.3 (Separate) intronic-gDNA gating

Independently of the leak fix, assess the §4.3 gating finding: intronic / no-candidate
gDNA loses its gDNA hypothesis (`new_cands == 0`). Quantify how much gDNA this withholds
and whether it should enter the EM (e.g. as a gDNA-only equivalence class).

Baseline to beat (current `main`, EXONIC gDNA→RNA leak, capture-on):
**ss99 ≈ co 20% / anti 1% (total 9.9%); ss50 ≈ 20% / 20%.** C-only ablation: co 14% /
anti 14% (total 13.4%) — the bar B+C must clear.
