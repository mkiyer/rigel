<!-- title: Global gDNA prior — count-space (Poisson-Gamma) derivation -->
# The global gDNA prior — a count-space (Poisson–Gamma) derivation

**Status:** derivation + partial validation. Replaces the ad-hoc, subset-selected, Gaussian-var~mean global
mean with a first-principles count-space rate. Core mechanism validated; application has known regressions
(§6). Prototype behind `RIGEL_GLOBAL_ALLNODES` (default off = current behavior).

## 1. The problem with the current global
The current global gDNA mean is a **mass-weighted average over a hand-picked node subset** (`self_solvable`),
and its per-node precision is a Gaussian `σ²(ρ_obs)` times the **`(eff/mass)²` Jacobian**. Three failures:
the subset selection is ad-hoc (and impossible to define for unstranded data); the Jacobian over-amplifies
the precision at high-mass nodes (over-pinning); and the mean has no graceful, length-aware path to **0** when
the truth is 0. The fix is to model the data as what it is — **discrete counts over a fixed exposure** — and
let the same machinery serve all conditions.

## 2. The node model (count + length, no Jacobian)
Each node *i* is a Poisson count of gDNA fragments over a **fixed exposure** (the gDNA effective length):
```
ĝ_i ~ Poisson(ρ · E_i),     ĝ_i = f_g_i · M_i  (the deconvolved gDNA count),   E_i = gDNA eff-len
```
With a `Gamma(a, b)` prior on the rate ρ, the node's **rate posterior** is `Gamma(a+ĝ_i, b+E_i)`:
```
mean      ρ̂_i = (a + ĝ_i)/(b + E_i)
precision τ_i = (b + E_i)² / (a + ĝ_i)
```
This is exactly the behavior asked for, derived rather than imposed:
- **zero count, long region** (`ĝ=0`, large `E`): mean → 0, precision → large ⇒ *confidently zero* (length-trust).
- **zero count, tiny region** (`ĝ=0`, small `E`): mean → `a/b` (prior), precision small ⇒ *honestly uncertain*
  (a single count is "infinitely more" than zero — the discrete nature is respected by the `a` pseudo-count).
- **high count**: precision ≈ `E²/ĝ` ⇒ the informative regions are confident.

The precision is in **rate space** — count and length only, **no `(eff/mass)²` Jacobian**.

## 3. Integrating nodes = summing Poisson sufficient statistics
Two nodes sharing a rate combine by **summing counts and summing exposures**:
`Gamma(a+ĝ_1+ĝ_2, b+E_1+E_2)`. This is the count-space analog of inverse-variance pooling (and mirrors how
count-space *messages* add). The combined precision is `(b+ΣE)²/(a+Σĝ)` — **not** the sum of the Gaussian
precisions; the Poisson sufficient statistic is the honest integrator.

## 4. The global mean = exposure-pooled rate over ALL nodes
```
ρ_global = (Σ_i ĝ_i + a) / (Σ_i E_i + b)          Var(ρ_global) = (Σĝ + a)/(ΣE + b)²
```
No node-subset selection — every node with mass contributes its `(ĝ, E)`. Three properties, all derived:

1. **Graceful zero, length-driven.** Empty intergenic contributes `(0 counts, huge E)`. The huge zero-count
   exposure sits in the denominator and pulls `ρ_global → 0`, with `Var → a/(ΣE)² →` tiny ⇒ **confidently 0**.
   This is the length-trust intuition: a genome-worth of zero-count exposure makes the pooled zero certain.
2. **It self-corrects to the truth by iteration.** Even when per-node deconvolution is phantom in a zero-gDNA
   *unstranded* library (`ĝ_i ≈ ½M_i` at flat-strand exons), the fixed point of the iteration is
   `ρ* = ρ*·(ΣE_signal / ΣE_total) ⇒ ρ* = 0` (because `ΣE_signal < ΣE_total`): each pass the intergenic
   exposure dilutes the phantom, `ρ_global` drops, it pulls the phantom exons down, repeat → **0**. This is the
   "initialize all-gDNA, let the data converge it down" design, now with a proof of where it converges.
3. **The informative regions govern.** High-count enriched regions dominate `Σĝ` (and the spread fit), so the
   estimate and its precision are set by where the signal actually is — the hybrid-capture intuition.

**Capture honesty.** Under capture `ρ_global` is the genome-*average* rate (the enriched counts spread over the
whole exposure), i.e. it under-states the enriched density. That is acceptable **iff** the global is then weak
where enrichment lives — which is the precision's job (§5): the population spread is large under capture, so the
global defers and the *enriched* nodes are governed by their own strand / imputation, not the global. A single
global cannot *be* the enriched density (that is per-region); its job is the zero-baseline + a weak tie-breaker.

## 5. The per-node prior precision = spread + mean-uncertainty
For an unknown node the predictive variance of its rate is the population spread plus the mean-estimate error:
```
σ²_prior = σ²_between(var~mean offset fit on the node rates)  +  Var(ρ_global)
```
- zero-gDNA uniform: `σ²_between ≈ 0`, `Var(ρ_global) ≈ 0` (huge ΣE) ⇒ confident ⇒ pins unknowns to 0.
- capture: `σ²_between` large (rate heterogeneity) ⇒ weak ⇒ defers to local.
- init (all nodes uncertain / few counts): `Var(ρ_global)` large ⇒ imprecise ⇒ "ridiculously imprecise at init,"
  tightening as nodes solve and `Σĝ/ΣE` sharpens.

## 6. The remaining gap — apply it in COUNT space (no Jacobian), and defer harder
The global *mean* is now count-space; the **application** still converts `(ρ_global, σ²_prior)` to a fraction
Gaussian via `τ = 1/(σ²·geom2)` (the `(eff/mass)²` Jacobian) and gates single-strand nodes by `1−(2κ−1)²`.
Measured (`RIGEL_GLOBAL_ALLNODES=1`, contained gDNA): **UNSTRAND-0 167,791 → 48,946** (−71%, the core win),
but **CAPTURE 1,694,088 → 1,501,443** (−11%) and **STRANDED-0 24 → 28,751**. Both regressions trace to the
*application*, not the mean:
- **capture −11%:** the genome-average `ρ_global` (low) leaks onto strand-RESOLVED exons because the
  `1−(2κ−1)²` defer (≈0.04 at κ=0.99) and the Jacobian-amplified precision are not enough to fully step aside.
- **stranded-0 → 28K:** ungating the global onto resolved single-strand nodes adds a small phantom where
  `ρ_global` is slightly > 0.

The principled fix is to make the application count-space too: the global enters node *j* as a **Beta/Dirichlet
pseudo-count** — `α_g = N_global · (ρ_global·E_j/M_j)` gDNA pseudo-fragments vs `α_r = N_global · (1−·)` — with
`N_global` the honest count-space confidence (from §5). This is **Jacobian-free**, and a strand-RESOLVED node's
large real strand count then dominates a small `N_global` *automatically* (no hand-tuned defer needed) — the
deference becomes emergent, the same count-currency principle as the messages. (The global, unlike a relay
message, is not forwarded, so the Dirichlet form is appropriate here — it does not hit the OQ1 density/precision
conflation.) `N_global` must scale so that under capture it is small (large spread) and in a confident zero
library it is large.

## 6b. Implementing the count-space application — what it revealed (the deep root cause)

I implemented the count-space application (behind `RIGEL_GLOBAL_ALLNODES`): mean = exposure-pooled rate over
the **self-solvable** nodes (excluding AMBIG — the global's targets — else `Σĝ` is self-sustaining: a nonzero
phantom fixed point, confirmed); application = a **Gamma prior on the gDNA count** `g_j=f_g·M`, shape `k≤1`
(allows g=0), rate `k/g_pseudo` (→∞ as `g_pseudo=ρ_global·E→0` ⇒ length-driven zero pin). Measured vs the
current Gaussian global (contained gDNA): **UNSTRAND-0 167,791 → 16,420** (−90%, the core win), **CAPTURE
1,694,088 → 1,626,280** (−4%, near-preserved via emergent deference), but **STRANDED-0 24 → 27,933**.

**The stranded-0 regression traced to a deeper, pre-existing bug — causally confirmed.** The phantom is
entirely on **strand-DECISIVE** AMBIG nodes (e.g. upos=16336, uneg=149 → obviously `+` RNA, f_g should be 0)
that nonetheless sit at **f_g=1**. Trace: the **+RNA imputation message** carries `μ=0.00, τ=1e7` and **pins
f_pos→0 ⇒ f_g→1**. That `τ=1e7` is the **(M/E)² geom2 Jacobian in the *message* precision** (mass≈16485 ⇒
(M/E)²≈3000×), which steamrolls the node's own decisive strand evidence (~1e4). The current strong Gaussian
global (also Jacobian-amplified, τ≈3.5e7) happened to **mask** this by pinning f_g=0 itself; the honest,
weaker count-space global **unmasks** it. **Causal proof:** capping the message precision at the dest count
collapses STRANDED-0 27,933→**1** and UNSTRAND-0 16,420→**9** — but worsens CAPTURE −9% (the messages carry
*real* gDNA signal there). So a naive cap is wrong.

**The unifying root cause:** the `(M/E)²` Jacobian over-amplifies precision **everywhere it converts a
density/rate belief to a fraction constraint — the global AND the imputation messages.** Fixing only the global
creates an imbalance (honest-weak global vs over-amplified messages) that surfaces the messages' version as a
phantom. **The complete fix is the count-space treatment of the WHOLE per-node solve** (the gDNA *count* as the
solve variable; messages as pseudo-counts) — the count-space solver pivot — not the global alone and not a
naive message cap. This is the precise, measured reason the global-only count-space change regresses, and it
re-converges the whole effort onto the count-space solve.

## 7. Status + next steps
- **Validated:** the exposure-pooled mean is the right object — it converges the zero-gDNA phantom down on its
  own (167K→49K and still descending across passes), count- and length-honest, no subset selection.
- **To do:** (a) the count-space (pseudo-count) application of §6 to remove the Jacobian and make deference
  emergent (fixes the capture + stranded-0 regressions); (b) confirm `σ²_between` is honestly large under
  capture; (c) the `a`,`b` Gamma-prior choice (currently `a=1`, `b=mean E`) through the no-magic-numbers review.
- **Honest scope:** capture-*unstranded* (enriched gDNA, no strand) remains the identifiability floor — the
  global is correctly weak there and the FL signal is the only true resolver (separate track).
