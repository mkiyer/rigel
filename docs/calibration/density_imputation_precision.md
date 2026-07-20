# Message precision for density imputation across the calibration chain

**Status: CLOSED — investigated, NOT pursued (2026-07-20).** A pass-0-vs-oracle diagnostic over all 32
`ambig_dense_10mb` scenarios (`scripts/debug/pass0_error_diagnostic.py`) showed the pass-0 error budget is
**80.5% under-determined** (honest high-variance — an *identifiability* floor on unstranded low-gDNA nodes, which
no message precision can fix and which the gDNA hyperprior is designed to resolve) and only **1.9%
confidently-wrong** (the sole precision-attackable slice, on 0.7% of node mass). The τ core had already moved
confidently-wrong mass from 10.4% → 1.4%. **Building the field-theoretic precision model below would chase a
≤1.9% error ceiling with large modeling + identifiability risk — so it was declined in favour of leaving pass-0
precision as-is and moving to the hyperprior.**

The **one durable takeaway retained** (a reviewer's correct catch): the honest density sampling variance is
`Var(log f_c) + 1/M` (composition + total count), **not** the per-component `1/(f_c·M)` — the latter
double-counts the split the τ core already owns. Equivalently, **effective length already enters precision through
the count** (`n = ρ·E`: a short node collects few fragments, so `1/M` already penalizes it) — there is no separate
"length" axis to add. The rest of this document is preserved as the derivation-and-review record that led to the
decision; it is **not** a design to implement.

---

*(Original framing below — retained for provenance.)*

**Status:** derivation + design proposal (2026-07-20). For prototyping, adversarial validation, and external
statistical review. Self-contained — no codebase familiarity assumed.

---

## 0. Purpose and scope

Rigel's calibration stage deconvolves each genomic node's fragment mass into a composition
**(f_g, f_r)** — the gDNA vs RNA fraction — by a belief-propagation (BP) sweep over a chain of nodes. Adjacent
nodes exchange **messages**: a source node imputes its estimate of a component's **density** onto a destination
node. Each message carries a **mode** (the imputed value) and a **precision** (how much the destination should
trust it). This document is about the **precision**.

The precision governs how confidently one node's belief flows to its neighbours. Getting it wrong is
consequential in two directions: too high → a node's error propagates as **confidently-wrong** mass that corrupts
downstream fits; too low → real information fails to propagate and nodes sit at their uninformative prior. We have
a working composition-precision model (the strand-evidence core, §2.1); the **open problem is the precision of an
imputed density magnitude** when the source and destination differ in (a) total density, (b) effective length,
and (c) composition — all of which vary enormously in real data, especially under hybrid-capture enrichment.

This doc: (§1) the geometry and notation; (§2) the **current** precision calculation, precisely; (§3) the
problem statement and the two discrepancies to model; (§4) the **derivation** of a density-imputation precision;
(§5) the **proposed solution** (concrete formulas + what changes); (§6) the prototype + validation plan; (§7)
open questions for review.

---

## 1. Geometry and notation

The genome is partitioned into a **chain of nodes** that alternate between two kinds:

- **Region** — a genomic interval of length `L` (an exon, an intron, or an intergenic stretch). It measures
  fragments **wholly contained** within it.
- **Boundary** — a genomic *point* between two regions (e.g. an intron↔exon junction). It measures **crossing
  flux**: fragments that span the point.

Each node's fragment mass is a mixture of components `c`:

- **gDNA** (`g`) — genomic DNA contamination, ~uniform genomically, strand-symmetric.
- **RNA** (`r`), split by strand into `+`/`−`; comprises **nascent** (unspliced pre-mRNA) and **mature** (spliced)
  populations. Boundaries additionally observe **spliced** fragments — a direct, motif-stranded measurement of
  mature RNA that crosses a junction as a spliced read.

**Effective length.** A component's density is `ρ_c = M_c / E_c`, where `M_c` is its fragment mass (a count) and
`E_c` its **effective length** — the expected number of start positions at which a fragment of that component's
fragment-length (FL) distribution produces an observable event at that node:

- Region, component `c`: `E_c^region = E_{FL_c}[max(0, L − ℓ + 1)]` (contained-fragment opportunity). **Note it
  depends on the component's FL distribution** and on `L`: for `L < ℓ` it is **0** (a region too short to contain
  a fragment of that length).
- Boundary, component `c`: `E_c^bnd ≈ E_{FL_c}[min(ℓ, flank)]/2` (crossing opportunity), `→ fl_mean_c/2` for long
  flanks.

Because **gDNA and RNA have different FL distributions**, `E_g ≠ E_r` at the same node — dramatically so for
short exons. This is central to §3–§4.

**Key symbols.**

| symbol | meaning |
|---|---|
| `ρ_c = M_c/E_c` | density of component `c` at a node (fragments per unit effective length) |
| `M_c`, `M = Σ_c M_c` | component mass (count) and total mass at a node |
| `f_c = M_c/M` | composition fraction; `λ = logit(f_g)` |
| `E_c^s`, `E_c^d` | source / destination effective length for component `c` |
| `n_c = ρ_c · E_c` | **component count** — the number of fragments supporting `ρ_c` (identically `M_c`) |
| `e(x)` | capture/coverage enrichment field (nucleic-acid-agnostic multiplier on density) |
| `S` | boundary spliced count (direct mature-RNA measurement, motif-stranded) |

The message is expressed in **log space** (log-density / log-fraction), where multiplicative enrichment and
frame factors become additive and precisions are approximately Gaussian.

---

## 2. The current message precision calculation

A message for component `c` (source `s` → destination `d`) is a soft Gaussian factor on the destination's
log-fraction `log f_c^d`, with **mode** `mo_c` and **precision** `pr_c`. The destination folds all incoming
factors into its belief. Only `pr_c` concerns us here.

### 2.1 The precision formula

```
    pr_c  =  1 / ( Var(log f_c^s)  +  1/M_src  +  σ²_transfer )
          =  n_eff / ( n_eff·(Var(log f_c^s) + σ²_transfer) + 1 ) ,   n_eff = M_src        (2.1)
```

Three additive variance terms:

1. **Composition variance `Var(log f_c^s)`** — how well the source knows the gDNA-vs-RNA **split**. Sourced from a
   reference-free **strand-evidence precision `τ`** (the "τ core"): `τ = I_strand + I_struct + relayed evidence`,
   where `I_strand` is the strand Beta-Binomial Fisher information (identically 0 for unstranded data — a count
   carries no intrinsic composition information) and `I_struct` marks structurally-certain nodes (intergenic
   locks). `Var(log f_g^s) = (1−f_g)²·(1/τ)` via the logit Jacobian. **This term is about the ratio, and we
   consider it solved.**

2. **Count term `1/M_src`** — the sampling (Poisson) precision. **`M_src` is the source's TOTAL facing mass**
   (`sm`, all components), not the component count.

3. **Transfer variance `σ²_transfer`** — the enrichment-crossing damping, from projecting each node's total
   log-density `log(mass/eff)` onto a fitted non-parametric density landscape (an NPMLE mixture of enrichment
   modes) and reading off `(μ_proj, var_proj)`:

   ```
       σ²_transfer  =  var_proj[d]  +  ( μ_proj[d] − μ_proj[s] )² .                          (2.2)
   ```

   The **mode gap** `(μ_proj[d] − μ_proj[s])²` is large when source and destination sit at very different total
   densities (a capture cliff); `var_proj[d]` is the within-mode ambiguity. **Belief-free, count-free, fit once.**

For the RNA message, an **honest clamp** replaces `n_eff` with the predicted destination component mass
`ρ_c·E_c^d` when the residual density is driven to the one-count floor (so a saturated absorption cannot be
laundered into a confident "no RNA"). The **spliced** measurement additionally credits precision `pr_c += S/(1 +
S·σ²_transfer)` — the direct mature count `S`, transfer-damped only (no composition variance, since spliced is
pure motif-stranded RNA).

### 2.2 What the current model does and does not capture

- ✅ **Composition (ratio) uncertainty** — via `τ` (term 1). Solid.
- ✅ **Enrichment / total-density discrepancy** — via `σ²_transfer` (term 3). This is the coverage cliff.
- ⚠️ **Count precision uses the TOTAL mass `M_src`, not the per-component count `M_c = f_c·M_src`.** This
  **over-credits** a component that is a minority of the mass (e.g. the gDNA density of an RNA-rich exon is
  credited with the full RNA+gDNA count).
- ❌ **No effective-length / resolution term.** A density measured over a **short** effective length (a boundary's
  crossing window; a tiny exon) is imputed onto a **long** frame (a whole exon body) with **no penalty for the
  frame mismatch** beyond the enrichment gap. Under real coverage non-uniformity this is the missing piece.
- ⚠️ **The spliced count `S` is credited raw** — honest as a *count* (`S = ρ_mature·E_spl`), but it inherits the
  same missing resolution term: a junction-window density imputed over a long exon body carries no frame penalty.

---

## 3. The problem statement

> **Given a source node's estimate of a component's density, with what precision should it be imputed onto a
> destination node — when the two nodes differ in total density (enrichment), effective length, and composition?**

Three concrete drivers, all large in real data:

### 3.1 Discrepancy in total density (enrichment / coverage)

Hybrid-capture enrichment (and, even without capture, ordinary coverage non-uniformity — a large, well-studied
phenomenon) makes the **absolute density** vary by 10²–10³× between adjacent nodes. A boundary might sit at
density 10 while its flanking exon sits at 1000. The **composition** is (locally) shared, but the density is not.
Imputing a density across such a cliff must be damped. *This is what `σ²_transfer` (2.2) targets today, on the
total-density landscape.*

### 3.2 Density is composition-dependent (different populations, different FL)

The density of a component is `ρ_c = M_c/E_c`, and **`E_c` depends on the component's fragment-length
distribution** (gDNA FL ≠ RNA FL). Consequences:

- Partitioning the **same** total unspliced count `M` over gDNA vs RNA yields **dramatically different densities**
  when `E_g ≠ E_r` — most severe for short exons where `region_eff_length(L, FL_c)` is small and steeply
  FL-sensitive.
- Some exons are **too short to contain a fragment** of a given population (`L < ℓ ⇒ E_c = 0`): that component's
  density is undefined and the node carries **no** information about it, even at high count.
- Therefore the density, and its precision, are **inherently per-component** — a single node has different
  reliabilities for its gDNA density and its RNA density.

### 3.3 Effective-length discrepancy (short ↔ long)

A boundary measures crossing flux over a **short** effective length (`~fl_mean/2`, ~75–100 bp); a region measures
containment over a **long** one (up to thousands of bp; or ~0 for a tiny exon). Imputation routinely crosses a
10–100× effective-length change. The intuition to formalize: **a density measured over a short window is a noisy
representative of a long frame** — and, symmetrically, a density imputed *onto* a short node is constrained by how
few fragments that node can hold. *Today's precision has no term for this.*

---

## 4. Derivation

### 4.1 The unifying identity: length enters precision **through the count**

The honest sampling precision of a log-density is the **component count**, and the component count **is** density
× effective length:

```
    Var(log ρ̂_c)  =  Var(log M_c)  =  1/M_c  =  1/(ρ_c · E_c) .                              (4.1)
```

So "the count is the honest precision" and "a short length makes the density imprecise" are the **same
statement**: a short `E_c` collects few fragments (`n_c = ρ_c·E_c` small) at a given density, hence low
precision. Worked example (owner's): a boundary at `1 count / 100 bp` and a region at `10 counts / 1000 bp` have
the **same density** (0.01) but precisions `n=1` (100%) vs `n=10` (32%) — the boundary is imprecise *because* its
short length yields `n=1`. **The remedy for §3.3's count side is simply to use the per-component count `M_c =
f_c·M` (2.1's term 2 currently uses the total `M`).**

But (4.1) is the **sampling** variance under a *uniform* density. It does **not** by itself penalize imputing a
short-window density onto a long frame — under strict uniformity, a perfectly-measured window density equals the
frame density. §3.1/§3.3 require modelling **non-uniformity** of the density field. We do so next.

### 4.2 A density-field model → the three-term imputation variance

Model the true per-component log-density field:

```
    θ_c(x)  =  μ_c  +  λ(x)  +  η_c(x) ,                                                      (4.2)
```

- `μ_c` — the global component mean.
- `λ(x)` — the **enrichment / coverage** field: nucleic-acid-agnostic (scales gDNA and RNA alike), the systematic
  large-scale non-uniformity that hybrid capture drives. This is the field the NPMLE landscape fits.
- `η_c(x)` — **residual non-uniformity**: the per-component coverage variation not captured by `λ` (mean 0,
  point-variance `σ²_η`, correlation length `ξ` on the order of a fragment length / probe footprint).

A node measures the component density by averaging `θ_c` over its effective-length window (weighted by fragment
geometry), collecting `n_c = ρ_c·E_c` fragments:

```
    θ̂_c^n  =  θ̄_c(window_n)  +  ε_c^n ,     Var(ε_c^n) = 1/n_c^n  (sampling) .              (4.3)
```

where `θ̄_c(window)` is the **window-average** of the field. For a field with point-variance `σ²_η` and
correlation length `ξ`, the window-average over effective length `E` has variance

```
    Var[ θ̄_η(window) ]  ≈  σ²_η · min(1, ξ/E)                                                (4.4)
```

— it saturates at `σ²_η` for `E ≲ ξ` (a boundary / tiny exon: no averaging) and falls as `σ²_η·ξ/E` for `E ≫ ξ`
(a long region: the non-uniformity averages down).

**Imputation error.** Predicting the destination's density from the source's, `θ_c^d − θ̂_c^s`, decomposes into
three independent contributions:

```
    θ_c^d − θ̂_c^s  =  [λ_d − λ_s]         (enrichment gap)
                    +  [η̄_c^d − η̄_c^s]    (residual non-uniformity, both windows)
                    −  ε_c^s              (source sampling)
```

giving the **three-term imputation variance**:

```
   ┌────────────────────────────────────────────────────────────────────────────────────────┐
   │  σ²_impute(c, s→d)  =   1/n_c^s               (1) source sampling  = 1/(ρ_c·E_c^s)        │
   │                      +  σ²_transfer           (2) enrichment gap   = (μ_proj gap)² + var  │
   │                      +  σ²_η·[ min(1,ξ/E_c^s) + min(1,ξ/E_c^d) ]   (3) non-uniformity     │
   └────────────────────────────────────────────────────────────────────────────────────────┘                  (4.5)
```

### 4.3 What each term captures, and the answers to §3

- **(1) Source sampling `1/(ρ_c·E_c^s)`** — the honest, **per-component**, length-dependent count precision
  (§4.1). Answers §3.3's count side and §3.2 (per-component `E_c` → per-component precision; `E_c=0` ⇒ ∞ variance
  ⇒ the message is correctly ignored). **Fixes term-2-uses-total in (2.1).**
- **(2) Enrichment gap `σ²_transfer`** — unchanged from (2.2). Answers §3.1 (the total-density / capture cliff),
  via the NPMLE landscape.
- **(3) Non-uniformity `σ²_η·[min(1,ξ/E_c^s) + min(1,ξ/E_c^d)]`** — the **new** term. Answers the residual of §3.3:
  a density measured over a short window, or imputed onto a short node, sees un-averaged coverage non-uniformity.
  **Both** effective lengths enter (inversely), and **per-component** (`E_g` vs `E_r`).

**Directionality (an explicit question raised in review).** Confidence is set primarily by the **source count**
(term 1): imputing *from* a long, high-count region is confident; *from* a short, low-count boundary is not. The
**frame** terms (2)+(3) are penalties on *either* end being short or enrichment-mismatched — so "short→long" and
"long→short" are **both** penalized whenever either endpoint is short. It is **not** one-directional: the penalty
is `1/E_c^s` (in term 1, sampling) plus `min(1,ξ/E_c^s)+min(1,ξ/E_c^d)` (in term 3, non-uniformity).

**Composition-dependence (§3.2) enters twice:** the per-component count `n_c = ρ_c·E_c` (term 1) and the
per-component effective lengths in term 3. A tiny exon with `E_g ≈ 0` (cannot hold a gDNA fragment) gives an
infinite gDNA-imputation variance — the gDNA message is suppressed, exactly as it should be, while its RNA message
(if `E_r > 0`) can still flow.

### 4.4 Note on the log-space subtlety

Under **strict local uniformity** (`σ²_η = 0`), term 3 vanishes and (4.5) reduces to `1/n_c^s + σ²_transfer` —
the count and enrichment only, with the effective length entering purely through `n_c = ρ_c·E_c`. In that regime
the raw Poisson count is exactly right and "length beyond count" would be double-counting. The non-uniformity term
is therefore **contingent on `σ²_η > 0`** — which is empirically overwhelming (coverage non-uniformity is large,
and hybrid capture makes it enormous), but the derivation is honest that term 3 is the *relaxation of uniformity*,
not a consequence of it.

---

## 5. Proposed solution

Replace the single count term `1/M_src` (and the raw spliced credit) with the three-term per-component variance
(4.5). Concretely, for each component `c` and edge `s→d`:

```
    n_c^s        =  ρ_c^s · E_c^s              # per-component source count (= M_c^s; NOT the total M)
    resol(E)     =  min(1, ξ / E)              # the resolution factor (4.4); ξ = a fixed length scale
    σ²_nonunif   =  σ²_η · ( resol(E_c^s) + resol(E_c^d) )
    σ²_impute_c  =  1/n_c^s  +  σ²_transfer  +  σ²_nonunif
    pr_c         =  1 / ( Var(log f_c^s)  +  σ²_impute_c )          # composition (τ) + density imputation
```

- **Composition (ratio) precision `Var(log f_c^s)`** — unchanged (the τ core). Kept **separate** from the density
  imputation variance: the two answer different questions (which split? vs how well is the density magnitude
  known and transported?).
- **Spliced** — the direct mature count is honest as term (1): use `n_mature = S` (it already equals
  `ρ_mature·E_spl`), but now carry terms (2)+(3) with `E_spl` as the source effective length, so a junction-window
  mature density imputed onto a long exon body is penalized by the resolution term. (Overdispersion cap, below,
  applies to `S` as well.)
- **Overdispersion (optional, axiom-free).** Molecular sampling is Beta-Binomial/NB over-dispersed; the τ core
  already deflates the strand count that way. For consistency, cap every count at its overdispersed effective
  value `n_eff = n/(1+(n−1)ω)`, so a high-depth count does not over-credit. This is depth-based (distinct from the
  length-based term 3) and parameter-free (`ω` is fitted).

### 5.1 Estimating the two scales (`σ²_η`, `ξ`)

Everything else in (4.5) is already computed (`ρ_c`, `E_c`, the NPMLE landscape). Two scales remain, and the
**estimation strategy is the main design decision**:

- **Option A — tie to the NPMLE landscape residual (self-fitting, no new constant).** `σ²_η` is, by construction,
  the density-field spread *not* explained by the enrichment landscape `λ`. The NPMLE already produces
  `var_proj[·]` (within-mode ambiguity) and the residual of each node's `log(mass/eff)` about the fitted modes.
  Estimate `σ²_η` from that residual (e.g. its variance after removing the mode structure), and take `ξ` at the
  fragment-length scale (the resolution floor). Fully self-fitting; couples term 3 to the enrichment fit.
- **Option B — one explicit fitted scalar.** Fit a single library-level `σ²_η` (and possibly `ξ`) by
  matching the observed node-to-node density disagreement to (4.5) — a method-of-moments on adjacent-node density
  residuals, stratified by effective length so the `resol(E)` shape is identified. Cleaner conceptually; adds one
  (or two) fitted scalars, not a per-case constant.

Both keep faith with the "no per-case magic numbers" constraint (a single library-level scale, fitted from the
data, is not a tuned constant). **Recommendation: prototype Option A first** (it introduces nothing new to fit),
and fall back to Option B if the landscape residual proves too coupled to the enrichment estimate.

### 5.2 What changes in the solver

1. Count term: total `M_src` → per-component `M_c^s = f_c^s·M_src` (and the honest clamp expressed per-component).
2. Add the non-uniformity term (3) with the fitted `σ²_η`, `ξ`, per component, using both `E_c^s` and `E_c^d`.
3. Spliced: carry terms (2)+(3) with `E_spl`; optionally the overdispersion cap on `S`.
4. `σ²_transfer` (term 2) unchanged.

All behind a default-OFF toggle; A/B before flipping (§6).

---

## 6. Prototype and validation plan

**Unit-level (Monte-Carlo, the ground-truth arbiter).** Simulate the density-field model (4.2): gDNA + RNA
fragments over region/boundary geometries with (i) controllable enrichment cliffs, (ii) controllable coverage
non-uniformity `σ²_η`, `ξ`, (iii) FL distributions that differ between gDNA and RNA, and (iv) short/tiny exons.
Measure the **true** imputation error variance node-to-node and check that (4.5) predicts it across the grid —
especially that the short-`E` and cliff regimes are captured, and that the uniform limit (`σ²_η=0`) reduces to
count+enrichment (no double-count). Reuse/extend `scripts/debug/cliff_hybrid_mc.py` (which already drives the real
`_fold_lambda` with configurable FL pairs and enrichment).

**Adversarial review.** (a) Verify the three terms are the *only* independent contributions (no missing
covariance between sampling and field). (b) Stress the `resol(E)` form at `E ≈ ξ` and `E ≪ ξ`. (c) Check the
estimator of `σ²_η` is not circular with the enrichment fit (Option A risk). (d) Confirm the per-component count
fix does not merely re-introduce, in a new place, the composition coupling the τ core already handles.

**System-level (the decision gate).** Pass-0 vs oracle on all 32 cached `ambig_dense_10mb` scenarios
(`scripts/debug/pass0_error_diagnostic.py`): mass-weighted `|Δf_g|` vs oracle by node type (**boundary / intron /
exon**) + **confidently-wrong mass %**. Gate to flip on: prior-free `f_g` accuracy improves (or is neutral) with
**no** node-type/condition regression, stranded controls flat, and confidently-wrong mass does not rise. Baseline
to beat (current production): exon 0.249 (unstranded 0.392), boundary 0.274, CW 0.6% (all-32 means).

**Sequencing.** This precision work is on **pass-0**, which must be maximally accurate *prior-free* (only the
inert reference + messages) before the gDNA hyperprior (NPMLE) is fit on its output and the refit is run. A
healthier pass-0 precision directly improves the hyperprior it feeds.

---

## 7. Open questions for review

1. **The three-term structure (4.5)** — is source-sampling + enrichment-gap + non-uniformity the correct and
   complete decomposition of the density-imputation variance? Any missing cross-term (e.g. sampling × field
   covariance when the same fragments inform both)?
2. **`σ²_η` estimation** — Option A (NPMLE landscape residual, self-fitting) vs Option B (one fitted scalar). Is
   the landscape residual a defensible estimator of the *residual* non-uniformity, or does it conflate with the
   enrichment field `λ` it is defined against?
3. **The resolution kernel `min(1, ξ/E)`** — is the window-average variance (4.4) the right functional form, and
   is `ξ ≈ fragment length` the right correlation scale, or should `ξ` itself be fit (Option B)?
4. **Per-component count** — replacing the total `M_src` with `M_c = f_c·M_src` is more honest, but `f_c` is a
   belief (uncertain). Should the count term use the *point-estimate* `f_c` (as proposed) or be marginalized over
   the composition posterior (coupling it to the τ variance)?
5. **Destination count** — (4.5) uses the destination effective length only in the non-uniformity term (3).
   Should the destination's own sampling (`1/n_c^d`, how few fragments it can hold) enter directly, or is that
   already the destination's own local likelihood (not the message's job)?
6. **Overdispersion** — apply the Beta-Binomial deflation to the density counts (as the τ core does for strand),
   or leave the count terms Poisson?

---

## Appendix — the one-paragraph summary

The precision of an imputed density is `1/(Var(log f_c) + σ²_impute)`, where the **composition** variance
`Var(log f_c)` (the gDNA/RNA split, from strand evidence `τ`) is kept separate from the **density-magnitude**
imputation variance `σ²_impute = 1/(ρ_c·E_c^s) + σ²_transfer + σ²_η·[min(1,ξ/E_c^s)+min(1,ξ/E_c^d)]`. The first
term is the honest per-component, length-aware sampling precision (a short effective length yields few counts);
the second is the enrichment/coverage cliff (the existing NPMLE-landscape transfer variance); the third is
residual coverage non-uniformity, which penalizes imputing across a short effective length at *either* endpoint.
Effective length enters precision through the count (`n = ρ·E`) and through the resolution kernel — never as an
ad-hoc factor — and everything is per-component because gDNA and RNA have different fragment-length distributions.
