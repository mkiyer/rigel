# Implementation plan — gDNA strand handling in the locus EM (v3)

**Status:** **Phase 0 complete (2026-06-06) → VERDICT: do NOT build the EM strand term.** The
gating measurement (§4 RESULTS) shows the co-stranded leak does **not** bias soft abundances —
mRNA is never inflated (always slightly *under*, dominated by a gDNA-independent −3% baseline),
and the only real per-locus errors are in the *unstranded* regime where strand carries no
information. Plus `ρ_d = 0` by construction (the decode uses the Binomial limit), so the clean
latent-κ model below is provably inert. The strand-balance term is unwarranted for
quantification. §2's clean model is retained as the *correct* formulation should a future need
arise; §4 RESULTS is the operative conclusion.
**Prereq:** [01_theory_and_fix.md](01_theory_and_fix.md) §4.1–4.3.

---

## 0. Reviewer critiques → how this version resolves them

| Critique | Resolution in v3 |
|---|---|
| Multi-strand trap (global neutralization breaks both-strand loci) | **Dissolved** — no neutralization. Each component (gDNA + every transcript) carries its own strand rate; multi-strand loci are handled by the same per-component machinery, no special case (§2). |
| "Neutralization is a hack; `log 0.5` is already neutral" | **Agreed and adopted.** `log 0.5` *is* the correct gDNA orientation emission (Bernoulli mean ½); it is not removed. The clean model has no neutralization step (§2–3). |
| Python grid solver in the C++ hot loop | **Eliminated.** The clean M-step is a **closed-form Beta update** (two digamma evals), not a grid/Newton optimization. Native C++, single implementation (§5). |
| Double-counting the prior | The strand evidence enters the split in exactly one place; calibration's pseudo-counts stay the conjugate Beta prior. But see §3 — the real risk is double-counting *across stages*, which the gating measurement checks. |
| VBEM vs MAP point estimate | The latent rate is a proper variational factor → native to VBEM; MAP-EM uses its mode. Same closed form (§5). |
| Spliced fragments in N_sense/N_anti | Excluded by construction — only unspliced, gDNA-eligible fragments have a gDNA hypothesis and contribute to κ_g (§5). |

## 1. The defect, stated precisely

The current per-fragment gDNA strand factor is a **constant `log 0.5`** — the point estimate
`κ_g = ½` of a gDNA orientation rate. That is the correct *mean*, but it is fixed: it is the
`ρ_d → 0` (zero-overdispersion) limit. The theoretically clean object is a **latent gDNA
strand rate** `κ_g ~ Beta(a_d, a_d)` (mean ½, overdispersion `ρ_d`), inferred per locus,
with the RNA components carrying their own rates at `p_sense`.

## 2. The clean model (no neutralization, no multi-strand special case)

Per locus, over **unspliced** fragments. Components: gDNA + transcripts `{k}`, mixing
`π ~ Dir(α)` (α from calibration). Latent per-component strand rates:
- `κ_g ~ Beta(a_d, a_d)`, `a_d = ½(1−ρ_d)/ρ_d` (gDNA, symmetric, overdispersion `ρ_d`);
- `κ_k ~ Beta` at mean `p_sense` (rel. to transcript strand σ_k), overdispersion `ρ_r`
  — or simply **fixed** `κ_k = p_sense` (RNA strand is known from calibration).

Fragment f: `z_f ~ Cat(π)`; orientation `o_f | z_f=c ~ Bernoulli(κ_c oriented)`; length/position
as today.

**Variational EM** (`q(z)q(π)∏q(κ_c)`):
- `q(κ_g) = Beta(a_d + S_g, a_d + A_g)`, with expected counts `S_g = Σ_f r_{f,g}·[o_f sense]`,
  `A_g = Σ_f r_{f,g}·[o_f anti]`, `N_g = S_g + A_g`.
- **E-step orientation factor** for gDNA (replaces the constant `log 0.5`):
  - sense: `ψ(a_d + S_g) − ψ(2a_d + N_g)`
  - anti:  `ψ(a_d + A_g) − ψ(2a_d + N_g)`
  (ψ = digamma). As `ρ_d → 0` (`a_d → ∞`) both → `log ½` — i.e. **this reduces exactly to
  today's behavior**; the change is the finite-`ρ_d` adaptation.
- **M-step**: update `(S_g, A_g)` from responsibilities → update `q(κ_g)`. Closed form.

This is one coherent model. There is **no neutralization** (the gDNA factor is a real
emission, not zeroed) and **no multi-strand special case** (every component, on either
strand, uses the same per-component rate; a fragment simply scores against each component's
`κ_c`). The reviewer's trap and the "hack" both disappear.

## 3. The premise is in doubt — read before building

Re-deriving exposed a problem the v2 plan glossed: **the model above, with the correct
symmetric prior, will likely reproduce today's behavior and *not* reduce the leak.** Three
linked reasons:

1. **A Beta(½, ρ_d) rate does not anchor "sense count = antisense count."** Observing
   antisense gDNA shifts the predictive toward *more antisense*, not toward sense. The
   "gDNA is 50/50, so antisense implies equal sense" intuition is about the *prior mean*,
   not a per-fragment pull. With small `ρ_d`, `κ_g` is pinned at ½ and the factor ≈ `log ½`
   (today). With large `ρ_d` it *learns the observed asymmetry* — the wrong direction.
2. **The only valid aggregate lever is the gDNA *fraction* `π_g`**, estimable from the
   locus strand split via the count×strand mixture (`joint_deconv` at locus scale). But that
   is **what calibration already computes at region scale**, and the diagnostic
   (§4.3/Q2) shows `gdna_em ≈ gdna_prior` per locus — the budget is ≈ right. A locus-level
   re-estimate is therefore likely **redundant** with the calibration prior.
3. **The leak is measured on hard argmax labels; the reported abundance is the soft EM
   count** (`estimator.count_em` = `em_counts.sum`, [estimator.py:518](../../src/rigel/estimator.py#L518)).
   Pool gDNA fraction is within ±2.2% and per-locus budget is ≈ right, so the *soft
   quantification may already be correct* and the "leak" may be a hard-label QC artifact
   (leak and siphon cancel in the soft counts).

**Worked check (why soft counts can be right while hard labels leak):** at a single-strand
locus with `M` mRNA and `A` sense-gDNA (+ `A` antisense-gDNA, all confidently gDNA), the
per-fragment EM assigns each sense fragment to gDNA with responsibility
`ρ = π_g·½·FL_g / (π_g·½·FL_g + π_m·p_sense·FL_m)`. With `FL_g=FL_m`, `p_sense=1`, and
`π_g/π_m = 2A/M`, this gives `ρ = A/(M+A)` exactly — **the soft gDNA-sense count = A, correct.**
FL differences then move the *hard* label (long gDNA → mRNA, short mRNA → gDNA) while the
*soft* count stays ≈ A. So the leak is a labeling artifact unless FL/`π_g` are actually off.

## 4. Gating measurement (DO THIS FIRST — before any production code)

Decide whether there is an abundance problem to fix at all:

- **A. Per-transcript abundance error.** For the leaky conditions, compare Rigel's
  `count_em` (soft) per transcript against simulated truth. Compute per-transcript and
  per-locus relative error, and the gDNA `count_em` vs true gDNA. Question: **is mRNA
  inflated at high-gDNA loci beyond the ±2.2% pool tolerance, or are the soft counts
  accurate despite the hard-label leak?**
- **B. If abundances are accurate** → the leak is a hard-label/annotation artifact. The EM
  change is unnecessary for quantification; at most, improve the *annotation* argmax (cheap,
  no model change), and stop here.
- **C. If abundances are biased** → localize: is `π_g` (`gdna_prior_count`) wrong at those
  loci (→ calibration fix), or is it FL overlap (gDNA vs mRNA length) (→ FL model)? Only if
  the bias traces to the *strand split being mis-estimated* does an EM strand term help —
  and then §5 is the clean way to add it.

Tooling: `scripts/debug/em_strand_gating.py` (joins `quant.tsv` `count_em` vs
`truth_abundances.tsv` `observed_mrna_fragments` per transcript/locus).

### 4. RESULTS (2026-06-06) → verdict: branch B (do not build)

Total mRNA (soft `count_em`) vs truth, by condition:

| condition | mRNA est | mRNA true | inflation |
|---|---|---|---|
| **gdna_none** ss99 cap-on (control) | 96 924 | 99 996 | **−3.07%** |
| gdna_high ss99 cap-off | 94 440 | 97 536 | −3.17% |
| gdna_high ss99 cap-on | 94 004 | 99 999 | −6.00% |
| gdna_high ss50 cap-off | 94 705 | 97 526 | −2.89% |
| gdna_high ss50 cap-on | 98 114 | 99 995 | −1.88% |

**mRNA is never inflated** — it is always slightly *under* (the FP-safe direction). The
**gdna_none control already shows −3.07%** with zero gDNA, so the few-percent under-count is a
**gDNA-independent baseline** (effective-length / unresolved fragments), not leak. The gDNA
conditions add little beyond it, and in the *under* direction. **The co-stranded leak does not
bias mRNA abundance** — leak and siphon cancel in the soft counts. Branch **B**: hard-label
artifact.

Per-locus, at **ss99** Δmrna is small and **not** proportional to gDNA load (Δmrna/gdna ≈
−0.14…+0.02). The only large per-locus errors are at **ss50 cap-on** (locus 1 −12 947; loci
5/8/4 +2–5 k) — real gDNA-vs-mRNA mis-splits, but in the **unstranded** regime where strand
carries no information, so **not addressable by any strand term**. Combined with `ρ_d = 0`
(the decode is the Binomial limit), the EM strand-balance term is confirmed **unwarranted and
inert**.

**What this redirects to (separate, if they matter):** (i) the gDNA-independent ~3% mRNA
under-count (eff-len / unresolved); (ii) the unstranded-capture per-locus split (FL power /
prior, not strand); (iii) the intronic-gDNA gating gap (§4.3 box). None is a strand fix.

## 5. (Not triggered) IF §4 had said "build it" — the native design

Implement the §2 latent-`κ_g` factor entirely in C++ (`em_solver.cpp`), no Python in the
hot loop, no duplicated implementation:
- Carry a per-fragment sense/antisense bit (from `exon_str` vs locus strand) into the
  equivalence classes — one extra `int8`/bit per candidate row.
- Accumulate `S_g, A_g` in the E-step alongside `em_totals` (two extra doubles per locus).
- E-step gDNA orientation factor = the two digamma terms (a `digamma()` already exists in
  `em_solver.cpp` for VBEM). **Closed form; no grid, no Newton.**
- M-step: `q(κ_g)` is implicit in `(S_g, A_g)`; nothing to optimize.
- `a_d` from `ρ_d` (calibration). RNA rates fixed at `p_sense` (today) unless §4 shows `ρ_r`
  matters. **No new constants.**
- **VBEM:** `κ_g` is a variational factor; its `E[log κ_g]` enters the responsibility exactly
  like the existing `digamma(α)` terms — consistent with `vbem_step`. **MAP-EM:** use the
  mode of `q(κ_g)`. Same code path, one branch on the weight source.
- **Multi-strand / AMBIG:** no special case — each transcript scores against its own `κ_k`;
  gDNA against `κ_g`. (If `κ_k` is fixed at `p_sense`, a both-strand locus naturally has a
  `+`-transcript favouring sense and a `−`-transcript favouring antisense.)
- **Unspliced only:** only unspliced fragments carry a gDNA hypothesis and contribute to
  `S_g/A_g`; spliced fragments score against transcripts only.

## 6. Phased plan

- **Phase 0 (gating, no production code):** §4 abundance measurement. **Decision gate:** is
  there an abundance bias, and does it trace to the strand split? If no → stop (or annotation-only).
- **Phase 0.5 (prototype, only if Phase 0 says build):** implement §2 in a small Python EM
  on one extracted leaky locus; confirm the latent-`κ_g` factor measurably improves the
  *soft per-transcript abundance* vs truth (not just the hard label). Sweep `ρ_d` to see if
  any value helps — if none does, the model is confirmed inert and we stop.
- **Phase 1 (native):** port §5 to `em_solver.cpp`; assert byte-identical output at
  `ρ_d → 0` (the limit must equal today).
- **Phase 2 (validate):** abundance accuracy (primary), leak/siphon (secondary), ss50
  unchanged, gdna_none FP unchanged, isoform guard, pool ±2.2%.

## 7. Open questions for review

1. **Is the hard-label leak a real quantification problem?** (Phase 0.) If the soft
   abundances are already accurate, do the hard `ZF` labels matter downstream (QC, IGV,
   filtering)? If only annotation, we fix the annotation, not the EM.
2. **`ρ_d` magnitude.** Calibration fits `ρ_d`; is it large enough that the latent-`κ_g`
   factor differs materially from `log ½`? If `ρ_d ≈ 0`, the clean model ≈ today and we
   should expect no change — worth checking the fitted value directly.
3. **RNA rate.** Keep `κ_k` fixed at `p_sense`, or also make it a Beta factor (`ρ_r`)?
   Default: keep fixed unless Phase 0 shows RNA overdispersion matters.
