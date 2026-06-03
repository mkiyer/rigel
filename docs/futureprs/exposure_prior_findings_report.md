# Exposure-Prior Robustness — Findings Report

**Status:** investigation synopsis for external review.
**Scope:** the per-region gDNA *exposure* sub-model of Rigel's calibration stage, its
sparse-data brittleness, an attempted adaptive-regularization fix, and a structural
coupling (a two-sided bifurcation) the fix unmasked. Ends with the open question we
most want reviewed: **how should the regularizing prior be set, without a hand-tuned
constant?**

This document is self-contained; no Rigel-specific knowledge is assumed.

---

## 1. Background: the calibration exposure model

Rigel quantifies RNA-seq libraries that are a mixture of **mRNA**, **nascent RNA**
(nRNA, unspliced pre-mRNA), and **genomic DNA contamination** (gDNA). Before the
per-locus quantification EM, a genome-wide **calibration** stage estimates a handful of
library-level hyperparameters and, per genomic *region*, deconvolves the observed
fragment mass into a gDNA part and an RNA part.

The genome is partitioned into `R` non-overlapping **regions** (exons, introns,
intergenic spans; in a human annotation `R` ≈ 10⁵–10⁶; in our unit-test toys `R` ≈ 30).
Each region `r` has a physical length `L_r` and an observed gDNA-attributed mass `M_g,r`
(a non-negative real, since fragments are fractionally allocated).

The gDNA count model is **Gamma–Poisson** (negative binomial):

```
M_g,r  ~  Poisson( ρ₀ · ω_r · L_r )
ω_r    ~  Gamma(s, s)              E[ω_r] = 1,  Var[ω_r] = 1/s ≡ φ
```

- `ρ₀` — global gDNA density (expected gDNA fragments per bp). One scalar per library.
- `ω_r` — the **per-region exposure**: a multiplicative deviation of region `r`'s gDNA
  density from the library average. `ω_r = 1` ⇒ region as contaminated as average;
  `ω_r > 1` ⇒ locally enriched (e.g. a hybrid-capture probe target); `ω_r < 1` ⇒
  depleted (e.g. capture off-target). The prior `E[ω_r] = 1` says *"uniform unless the
  data say otherwise."*
- `φ = 1/s` — the **exposure dispersion**: the region-to-region variance of `ω`. `φ → 0`
  ⇒ perfectly uniform library (the Poisson limit); large `φ` ⇒ highly non-uniform
  (strong capture enrichment).

The per-region exposure posterior is conjugate:

```
ω̂_r = E[ω_r | M_g,r] = (s + M_g,r) / (s + ρ₀ L_r)         — a convex pull toward 1,
                                                              prior weight s, data weight ρ₀L_r.
```

`φ` is fit by an empirical-Bayes Gamma-shape MLE from the per-region posteriors
(stationary condition `log s − ψ(s) = mean_r(E[ω_r] − E[log ω_r]) − 1`, solved by
bisection; `ψ` = digamma). `ρ₀` is a closed-form weighted ratio
`ρ₀ = Σ M_g,r / Σ ω_r L_r`. These are interleaved in an outer EM with the gDNA-vs-RNA
soft allocation (E-step) and the other hyperparameters.

**Why `ω` matters downstream.** The calibration hands the per-locus EM two Dirichlet
scalars that set the gDNA-vs-RNA split, both proportional to the local exposure
`e_r = ω_r · L_r`. So a region's `ω` controls *both* how much gDNA the prior expects
there *and* (via the gDNA component's effective length `≈ Σ ω_r L_r`) how attractive the
gDNA state is in the EM. This dual role is central to §6.

---

## 2. The problem: sparse-data brittleness

On small synthetic scenarios the exposure weights are **brittle**: they collapse to
near-zero in some regions and inflate in others, even when the *true* library is uniform
(`ω ≡ 1`). Downstream this mis-routes RNA into gDNA (or vice-versa).

Concrete failing case — `nrna_dc` (`R = 33` regions, one expressed 8-exon gene + one
silent control gene, uniform gDNA contamination, ~2000 fragments):

| metric | expected | observed (pre-fix) |
|---|---|---|
| total RNA recovery | 100% | **20%** |
| gDNA fragments | 596 | ~1700 (**over-called**) |

**Mechanism.** With sparse data, two effects inflate the *apparent* dispersion of `ω`:
(a) Poisson sampling noise in `M_g,r`, and (b) deconvolution mis-allocation (the
gDNA-vs-RNA E-step is itself uncertain when coverage is thin). Both make regions *look*
more variable than they are ⇒ the MLE returns a large `φ` (small `s`) ⇒ the prior pull
toward `ω = 1` weakens ⇒ `ω` spreads further ⇒ larger apparent dispersion next
iteration. The feedback is self-reinforcing.

**It is specifically a sparse-data artifact (verified).** Replicating the same gene pair
`K` times into a larger genome (more regions *and* proportionally more fragments)
dissolves the failure:

| scale | regions | RNA recovery |
|---|---|---|
| K=1 | 33 | 20% |
| K=5 | 161 | 91% |
| K=20 | 641 | 104% |

Robustness at scale is simply the law of large numbers: with enough regions the inflated
tail of `ω` becomes harmless. Real genomes (10⁵–10⁶ regions) sit deep in the robust
regime; the unit-test toys are pathological. This is the key reassurance — **the model is
correct in the large-data limit; the question is purely how to degrade gracefully toward
the sparse limit.**

---

## 3. The attempted fix: evidence-weighted adaptive regularization

Idea: shrink the empirical dispersion signal toward 0 (i.e. `ω → 1`, uniform) by an
amount that *adapts to how much evidence supports the dispersion*, so sparse libraries
shrink hard and dense libraries keep their real structure.

Two changes to the `φ` M-step (`update_exposure_dispersion`):

1. **Evidence-weighting.** Weight each region's dispersion contribution by its gDNA mass
   `M_g,r`. Empty/sparse regions carry no exposure information, and — critically — the
   *collapsed* `ω → 0` regions (which have `M_g ≈ 0`) must not be allowed to inflate the
   signal, since that is the very feedback loop of §2.

2. **Adaptive shrinkage.** Scale the (evidence-weighted) signal `c` by
   `c_reg = c · N_eff / (N_eff + κ)`, where `N_eff` is the **Kish effective number of
   gDNA-bearing regions** `(Σ M_g,r)² / Σ M_g,r²` and `κ` is a prior-strength constant.
   Sparse (small `N_eff`) ⇒ strong shrinkage ⇒ `ω → 1`; dense (large `N_eff`) ⇒ empirical
   dispersion preserved.

`κ` is a single configuration parameter (`exposure_prior_pseudocount`); `κ = 0` recovers
the raw evidence-weighted MLE. **`κ` being a hand-set constant is the central weakness —
see §7.**

---

## 4. Validation harness and baselines

We built a three-regime harness (`scripts/debug/exposure_harness.py`) spanning the axes
that matter:

1. **SPARSE + downsampling sweep.** The `nrna_dc` toy at fragment budgets
   {500…16000}. Diagnostic: recovered `φ` and RNA recovery vs depth.
2. **DENSE-UNIFORM scale-up.** `K ∈ {1,5,20}` copies of the gene pair; uniform gDNA.
   Tests the law-of-large-numbers regime; the fix must not regress it.
3. **DENSE-CAPTURE.** Six genes, a synthetic hybrid-capture probe set over three of
   them, gDNA focally enriched at the targets. The critical guard: a fix that shrinks
   `ω → 1` everywhere would **destroy real capture enrichment**, so this regime must keep
   `ω ≫ 1` at the targets.

**Pre-fix baselines:**

- *Sparse* is **chaotic**, not a clean monotone collapse — RNA recovery bounces
  108→27→20→17→102→46% across depth, with `φ` flat (~2–6). At small scale the exact
  fragment count lands the deconvolution in different basins; sporadic.
- *Dense-uniform* robust: K=5 → 91%, K=20 → 104%.
- *Dense-capture*: the current calibration **does** recover focal enrichment —
  `ω` at captured regions mean **6.65** (max 11), at uncaptured **0.02**; gDNA recovery
  866/2740. So capture works qualitatively even pre-fix.

---

## 5. κ-sweep results

Running the fix across `κ ∈ {0,20,50,100,200}`:

| regime | metric | κ=0 | κ=20 | κ=50 | κ=100 | κ=200 |
|---|---|---|---|---|---|---|
| sparse headline (n=2000) | RNA recovery | 47% | **96%** | 127% | 128% | — |
| | gDNA (exp 596) | 1340 | **651** | 215 | 205 | — |
| dense-uniform K=5 | RNA recovery | 91% | 91% | — | 131% | — |
| dense-uniform K=20 | RNA recovery | 96% | 98% | — | 99% | — |
| dense-capture | ω captured (mean) | 6.28 | 5.39 | 4.39 | 4.26 | 2.39 |
| | gDNA recovery /2740 | 925 | 1012 | 1139 | 1256 | 1599 |

Observations:

1. **Evidence-weighting alone (κ=0) already helps** the sparse mid-range (the pre-fix
   `nrna_dc` was 20%; κ=0 gives 47%) by refusing to let the collapsed `M_g ≈ 0` regions
   drive the dispersion.
2. **Dense + capture are robust across all κ.** Capture enrichment is preserved
   (`ω` captured stays 2.4–6.3, never crushed to 1), and gDNA recovery there actually
   *improves* with stronger shrinkage. The realistic-data goal is met by a wide range of
   κ.
3. **The sparse toy is fragile.** It is *correct* near κ≈20 (RNA 96%, gDNA 651≈596) but
   flips to ~128% over-recovery by κ≈50–100, and even κ=20 vs κ=30 flips the same
   scenario (96% → 127%). The "right" κ for the toy is a knife-edge.
4. The shrinkage **over-corrects the mid-scales** at large κ (dense-uniform K=5: 91% →
   131%). A κ strong enough to tame the 33-region toy harms the 161-region case.

---

## 6. The structural finding: a two-sided gDNA bifurcation

The over-recovery at large κ is not a tuning nuisance; it exposes a real coupling.
Full-chain trace of `nrna_dc` at κ=100 (so `φ` floors and `ω → 1` everywhere):

```
ρ₀ = 5.1e-5            (true ρ₀ ≈ 0.04 — ~800× too small: COLLAPSED)
φ  = 1e-8 (floored)  ⇒ ω = 1.000 in every region
per-region called gDNA mass ≈ 0 everywhere
gDNA observed 205 (expected 596);  freed mass → nRNA 1175 (expected 828)
⇒ RNA over-recovery 128%
```

So forcing `ω → 1` makes the global gDNA density `ρ₀` **collapse to ~0**, and the gDNA
that should have been called leaks into (nascent) RNA.

**Why `ω → 1` collapses `ρ₀`.** gDNA *inside a high-RNA region* is statistically
near-invisible: both the count channel (the region's mass is dominated by RNA) and the
strand channel (the region is strongly stranded) say "RNA." The exposure *dispersion*
was the degree of freedom that let the model **place expected gDNA where it is
detectable** (introns, intergenic) instead of smearing it uniformly. With `ω ≡ 1`,
expected gDNA `= ρ₀ L_r` is spread uniformly, including across the regions where it can't
be seen; the E-step then calls ≈0 gDNA there; `ρ₀ = Σ M_g / Σ ω L` follows the called
mass down; and the loop converges to `ρ₀ → 0`.

**Two-sided trade-off.** The exposure freedom mediates *two opposing* gDNA failure modes:

- **Too much `ω` freedom** (κ small): in low-`ω` regions the gDNA component's effective
  length `Σ ω L` shrinks, making the gDNA state spuriously *attractive* in the EM ⇒ gDNA
  **over**-called (the pre-fix failure of §2).
- **Too little** (κ large, `ω → 1`): `ρ₀` collapses ⇒ gDNA **under**-called (this section).

κ≈20 simply threads the needle between two attractors — hence the knife-edge. Notably,
`ω → 1` *cures* the first failure (the trace shows `gdna_eff_len` = full physical span,
no collapse); we merely traded one failure for the other.

This is a **pre-existing latent instability** that the regularization *reveals* rather
than creates: pre-fix, `ρ₀` stayed alive only because free `ω` over-called gDNA. The
realistic regimes never approach either edge (enough evidence to stay in the basin), so
they are robust regardless.

---

## 7. The open question: how should the prior be set?

`κ` as a hand-tuned constant is **not an acceptable solution**, for a fundamental reason:
there is no single value that can be right across the range of RNA-seq protocols and
sequencing depths. The §5 sweep shows the "optimal" κ already differs between a 33-region
and a 161-region toy; across real protocols (poly-A vs total vs hybrid-capture) and
depths (1M vs 500M reads) it would have to vary by orders of magnitude, and may even
differ locus-to-locus. A global constant `κ = 20` is a fitted artifact of one toy.

What we *have* established is the right *shape* of the solution: **the prior must adapt to
data sparsity**, and sparsity has two independent axes:

- **Region sparsity** — few regions `R` (small genome / sparse annotation).
- **Coverage sparsity** — few fragments per region (low `ρ₀ L_r`, shallow gDNA).

Both must drive the prior: with little of *either*, shrink `ω → 1`; with plenty of both,
trust the empirical dispersion (and recover genuine capture structure). The current
`N_eff` (a Kish *count* of gDNA-bearing regions) captures region sparsity but **not**
coverage sparsity, and the crossover point `κ` is still a free constant.

**The desired property:** graceful behavior across the full sparsity continuum, derived
*empirically from the data's own sampling properties* — not set by a tunable knob. We are
developing a theoretical foundation for this (a noise-floor / Fisher-information account:
the dispersion estimate has a sampling standard error
`SE(φ̂) ∝ sqrt(2 / Σ_r (ρ₀ L_r)²)`, in which `Σ (ρ₀ L_r)²` is a *single* quantity that
encodes both the region count and the per-region coverage; the prior should believe `φ̂`
only to the extent it exceeds this noise floor, which makes the shrinkage parameter-free).
That derivation is the subject of a companion note and is **the part we would most value
external scrutiny on.**

There is also the separate, coupled `ρ₀`-collapse instability of §6: even a perfectly
set exposure prior that drives `ω → 1` in the sparse limit will trip the `ρ₀` collapse
unless the gDNA-density estimate is independently stabilized. A complete solution likely
needs both (a) a principled, data-derived exposure prior and (b) protection of `ρ₀`
against circular collapse.

---

## 8. Questions for the reviewer

1. Is the noise-floor framing (`SE(φ̂) ∝ sqrt(2/Σ(ρ₀L_r)²)`; shrink by signal-to-noise)
   the right principled route to a **parameter-free** adaptive prior, or is there a
   cleaner standard estimator (REML, a proper hierarchical hyperprior with the
   hyper-scale itself integrated out, a score-test shrinkage)?
2. Is there a recognized treatment for the **two-sided bifurcation** (§6) — jointly
   estimating a global rate `ρ₀` and per-unit exposures `ω` when the per-unit signal is
   confounded with a much larger competing component (RNA)? This feels like a known
   identifiability pathology.
3. Is empirically *learning* the optimal-prior law from a scenario sweep (region count ×
   coverage → optimal shrinkage) a sound fallback if the closed form is intractable, or a
   red flag that the model is mis-specified?

---

### Appendix: reproduction

- Harness: `scripts/debug/exposure_harness.py --regime {sparse,dense,capture,all} [--kappa K]`
- Full-chain trace: `scripts/debug/trace_nrna_dc.py --gdna 20 --nrna 70 --ss 0.65 --nfrag 2000`
- Implementation: `update_exposure_dispersion` in `src/rigel/calibration/mstep.py`;
  parameter `exposure_prior_pseudocount` in `src/rigel/config.py`.
- All commands run inside the `rigel` conda environment.
