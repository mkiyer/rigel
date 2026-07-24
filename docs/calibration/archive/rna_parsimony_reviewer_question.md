# For a statistical reviewer: the RNA-parsimony reference term (½·log(1−f_g) vs −log(1−f_g))

**What we need:** help reconciling a discrepancy between a derived Jeffreys reference prior and a heuristic term
that was empirically load-bearing for false-positive control. We suspect the derivation assumed something our
current pipeline no longer provides, but we want the theory to close cleanly. This document is self-contained.

---

## 1. Setup

Rigel deconvolves each genomic **node** (a region of the genome, or a splice-junction boundary) into a
composition — the fraction of its **unspliced** fragment mass that is genomic-DNA contamination, `f_g ∈ [0,1]`,
vs RNA. Per node we maximize/integrate a log-posterior over `λ = logit(f_g)`:

```
ψ(λ) = strand-likelihood(λ)  +  gDNA-arm(λ)  +  RNA-arm(λ)  +  cross-node messages(λ)
```

- **strand-likelihood** — a two-component Beta-Binomial on the node's ± strand counts. Its Fisher information
  for `f_g` is `∝ N·(2κ−1)²`, where `κ` is the library's sense-fraction. **At κ=½ (an *unstranded* library) this
  is identically zero for any depth** — the strand carries no information about `f_g`.
- **gDNA-arm / RNA-arm** — each is a fitted log-prior on that component's rate when available, else a
  **reference** term. The gDNA hyperprior (a population prior on the absolute gDNA rate) is fit *after* an
  initial pass and applied only in a later refit.
- **messages** — belief-propagation between neighbouring nodes carrying each other's solved densities.

**The invariant we hold (count-zero-information).** A fragment *count* carries **zero** intrinsic information
about gDNA-vs-RNA composition; it may enter only as *precision*. The only intrinsic composition signal is the
strand tilt. Consequently, **at an unstranded node there is no likelihood evidence about `f_g` at all** — the
posterior is set entirely by the prior/reference and the messages.

---

## 2. The two competing RNA-side terms

**(A) The current, derived reference — Beta(½,½), coefficient ½.** Our reviewer-signed-off derivation
(`reference_prior_derivation.md` §10) sets the reference to the Jeffreys/Berger–Bernardo prior on the
composition simplex. After a change of variables to the log-odds `λ` (the latents are modeled as log-rates, and
a log-rate↔logit Jacobian cancels one term per component), the reference lands as

```
gDNA-arm_ref = +½·log f_g ,     RNA-arm_ref = +½·log(1 − f_g)          (the current code)
```

**(B) The shipped v0.7.1 heuristic — coefficient 1.** The previously released solver instead used

```
RNA-arm = −log(1 − f_g)  =  +log( 1 / (1 − f_g) )                       (coefficient 1)
```

bundled into its gDNA-prior evaluation. This term was **empirically load-bearing for false-positive safety**:
with it, false positives on a DNA-free library were ~2k fragments; **without it, ~585k**. The current ½ reference
does **not** reproduce that guard.

The two differ in coefficient (½ vs 1) and the released one was described in its own code as making "gDNA the
residual after a typical-magnitude RNA" — i.e. a *parsimony* pull toward `f_g → 0`, not merely a reference
measure.

---

## 3. The empirical symptom that forces the question

We measured the **initial pass** (no fitted gDNA prior yet) on a **true-zero-gDNA** library (`f_g` should be 0
everywhere), **unstranded (κ=½)**:

- Every single-strand *expressed* node (it carries a transcript, so its unspliced mass is predominantly RNA)
  solves to `strand=0.51`, `local(+reference)=0.51` — i.e. it **sits at the reference median ≈ 0.5**. With κ=½
  the strand is flat, so there is nothing to move it, and **the ½ reference's mode/median is ≈ 0.5**.
- Result: true-RNA nodes are called `f_g ≈ 0.47` (median), a large **systematic** false positive (~10⁵–10⁶
  fragments; ~29–66% of node mass confidently mis-assigned).
- On a **stranded** version of the same library (κ≈1) there is essentially no such false positive — the strand
  likelihood carries `f_g → 0`.
- With **nascent RNA present**, neighbouring boundaries also sit at ≈0.5 and *emit gDNA messages*, so the
  over-call **cascades** and becomes confident; with **nascent absent**, the nodes simply rest at the reference
  ≈0.5 with high variance. Both are rooted in the same reference behaviour.

So an *uninformative* reference gives `f_g ≈ 0.5` at a strand-silent node, whereas the truth is `f_g = 0`, and
the shipped `−log(1−f_g)` (coeff 1) was the thing that used to pull it down.

---

## 4. Our current understanding (from an internal derivation pass — please pressure-test)

Three independent internal derivations (reference-prior, rate-space, and empirical-Bayes framings) converged on
the following. **We now believe the "½ vs 1" framing was a false comparison** — they are different priors with
*opposite* effects, not weak-vs-strong versions of one reference. We want an external statistician to confirm or
refute, and especially to adjudicate the one remaining tension (§5).

- **`+½·log(1−f_g)` is the correct reference and should not change.** It is the proper Beta(½,½) Jeffreys prior:
  each component's rate-Jeffreys `p(ρ_c)∝ρ_c^{−½}`, written as a density in log-rate, is `+½·log ρ_c`; the
  log-rate→logit Jacobian cancels once per group (`reference_prior_derivation.md` §10, residual R=0), leaving
  `ψ_ref = ½·log f_g + ½·log(1−f_g)`. It is L-invariant and marginally consistent between the single-strand and
  the 3-part (f_pos,f_neg,f_g) node classes (the tilt is info-orthogonal; its θ=arcsin τ Jacobian cancels).
- **At κ=½ the posterior *is* the prior, and Beta(½,½) is symmetric ⇒ median 0.5.** No *symmetric* coefficient
  change moves it (Beta(0,0), Beta(½,½), Beta(1,1) all → 0.5); only an *asymmetric* reference moves it, which is
  a disguised "gDNA is rare" **abundance ramp** (the improper `+½·λ` this project already removed twice). So
  **f_g=0.5 at a strand-silent, zero-gDNA node is the honest output of an uninformative reference — not a bug —
  and inflating the coefficient is either inert or forbidden.**
- **`−log(1−f_g)` (the shipped "coeff 1") is NOT a stronger RNA parsimony — it is its improper opposite.** It is
  monotone *increasing* in f_g (→+∞ at the f_g→1 gDNA vertex): the change-of-variables measure of a Haldane
  (`∝1/ρ_r`) linear-rate prior; as a Beta exponent, `Beta(·,−1)`, giving `Beta(½,−½)` — non-integrable at f_g→1
  (median → 0.9998). In isolation it makes the zero-gDNA over-call *worse*. In the shipped code it did two other,
  real jobs (an intron-parsimony "bleed-stopper" that prevented the *opposite*-polarity FP — RNA hallucinated
  into introns — and the improper RNA companion of the KDE, made integrable only by the KDE's Gaussian tail).
  **The "585k→2k FP" A/B toggled the whole gDNA-arm block (KDE + ρ_bg-anchored global term), not the
  coefficient** — attributing the guard to the coefficient was a synecdoche.
- **The real regression is a *deleted term*, not the reference.** The shipped solver applied a
  `ρ_bg`-anchored gDNA-rate prior (`_global_logprior`) at **both** passes: a **two-sided**, one-pseudo-obs-capped
  Gaussian on `log f_g` with **mode = log(ρ_bg / ρ_tot)**. For dense exons in a DNA-free library (ρ_bg≈0) that
  mode sits low, breaking the flat-strand tie **downward toward RNA**. A refactor collapsed "prior-free initial
  pass" to "reference-only" and dropped it (`gdna_prior=None` ⇒ no ρ_bg at pass-0) → pure reference → 0.5.
- **Circularity note:** the gDNA hyperprior is fit on `f_g·mass` over these same κ=½ nodes, so a pass-0 stuck at
  0.5 teaches the refit a spurious half-density gDNA mode. Grounding pass-0 fixes both the pass-0 FP and this
  refit contamination.

## 5. The one question we most need adjudicated — TWO-SIDED vs ONE-SIDED ρ_bg

Restoring the ρ_bg anchor at pass-0 is our proposed fix, but there is a real statistical tension we cannot
resolve internally:

- **Breaking the flat-strand tie *downward* at true-zero requires the anchor's UPPER side** (the two-sided
  Gaussian mode at `log(ρ_bg/ρ_tot)`, which is low when ρ_bg≈0). A **one-sided lower floor** on ρ_bg (which
  `background_reference_derivation.md` §8 deliberately adopted) is **provably inert at true-zero** — it leaves the
  median at 0.5000.
- **But §8 chose one-sided *for capture-safety*:** pinning ρ_g at a capture-*depleted* ρ_bg via a two-sided
  anchor risks **false-negating genuine captured gDNA** (an `gdna300 + capture` under-call). The one-pseudo-obs
  stability cap is meant to keep the anchor weak enough to avoid this while still breaking the DNA-free tie.

**Questions:**
1. Is `+½·log(1−f_g)` indeed the correct objective/reference prior here (L-invariant, marginally consistent), or
   does the correct prior carry a different `(1−f_g)` power? (We believe it is correct — confirm/refute.)
2. Is a **two-sided, one-pseudo-obs-capped** ρ_bg-anchored Gaussian on `log f_g` (mode `log(ρ_bg/ρ_tot)`) the
   principled way to break the κ=½ tie toward RNA at true-zero — and is the stability cap *provably* sufficient
   to avoid false-negating real captured gDNA at high-gDNA+capture? Or is the one-sided floor correct and some
   *other* mechanism should break the tie?
3. Is it legitimate, under a strict "the fragment **count** never votes composition" invariant, for the prior to
   favour `f_g→0` at an expressed unstranded node via a ρ_bg-anchored *rate* prior? (We believe yes — ρ_bg is
   population/structural information, not the node's own count — but want the ruling.)
4. Does grounding pass-0 with a two-sided ρ_bg anchor introduce its *own* downward bias into the refit
   hyperprior (learning too-little gDNA under capture) — trading one contamination for another?

We can run any targeted A/B the reviewer wants — in particular the two corners that decide §5.2: the `gdna_none`
over-call and the `gdna300 + capture` under-call. We have a controlled sim suite with oracle truth.
