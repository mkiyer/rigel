# Calibration precision: the over-confidence failure mode, and the fix

**Status:** diagnosis + proposed fix, for review before implementation.
**Scope:** the three places calibration estimates a *variance* (→ precision) that feeds the per-node
simplex solve, why two of them can become catastrophically over-confident, and how to make them honest.

---

## 0. Why precision is the whole game

Every node (a region or a boundary) owns a slice of **unspliced** fragment mass `M`, which the calibrator
partitions on the 2-simplex into `(f₊, f₋, f_g)` — sense-RNA / antisense-RNA / gDNA, summing to 1. The
per-node posterior is a **product of evidence terms**, each entering its log-posterior as a
precision-weighted penalty `−½·τ·(f − μ)²` (plus the strand likelihood, which is curved). Concretely the
solver maximizes (over the simplex):

```
log ψ(f₊,f₋,f_g) =  L_strand(tilt)                        # intrinsic, curved
                  −  ½·τ_g   ·(f_g − μ_g  )²              # gDNA imputation message
                  −  ½·τ₊    ·(f₊  − μ₊  )²              # +RNA imputation message
                  −  ½·τ₋    ·(f₋  − μ₋  )²              # −RNA imputation message
                  −  ½·τ_glob·(f_g − μ_glob)²            # node-class prior (global, at AMBIG/intergenic)
                  +  [spliced lower-bound, Jeffreys at single-strand]   # (omitted here)
```

This is a precision-weighted average: **the term with the largest `τ` wins.** If any `τ` is wrongly
inflated — say `1e7` when it should be `~3` — that term *pins* `f_g` to its `μ` and silently overrides
everything else, including a sharp, correct strand likelihood. That is the entire failure mode of this
document. Honest `τ` is not a nicety; it is the correctness condition for the solver.

---

## 1. The three variance models

### 1.1 The strand Beta-Binomial likelihood — `simplex._mixture_strand_loglik`

Of `N = u₊ + u₋` unspliced fragments at a node, a fraction `f_g` are gDNA (genomic strand symmetric ⇒
plus-rate ½), `f₊` are +RNA (plus-rate κ = library sense fraction) and `f₋` are −RNA (plus-rate `1−κ`).
The mixture plus-strand probability is

```
p = ½·f_g + κ·f₊ + (1−κ)·f₋  =  ½ + (κ−½)·(f₊ − f₋)  =  ½ + (κ−½)·t
```

where **`t = f₊ − f₋` is the tilt**. The second equality is an identity (substitute `f_g = 1−f₊−f₋`),
and it is the single most important fact about strand information:

> **`p` depends only on the tilt `t`, never on `f_g` directly.**

The likelihood is a Gaussian approximation to the (overdispersed) Beta-Binomial count `u₊ ~ (N, p)`:

```
mean = N·p
var  = N·p(1−p)                 # Binomial sampling
     + (N·f_g)²·¼·od_g          # gDNA overdispersion excess (μ(1−μ)=¼ at rate ½)
     + (N·f₊)²·κ(1−κ)·od_r      # +RNA overdispersion excess
     + (N·f₋)²·κ(1−κ)·od_r      # −RNA overdispersion excess
log L_strand = −½·(u₊ − mean)²/var − ½·log var
```

`od_g, od_r` are the **Beta-Binomial overdispersions**, fit once by method-of-moments from seed regions
(`gdna_strand.py`): they inflate the variance above Binomial to account for fragment clustering, so that a
finite count carries *less* than its naïve Binomial information. Mean ½ for gDNA, κ for RNA, applied
symmetrically, so an unstranded library (κ=½) is uninformative by construction.

**Derived precision (the key consequence).** Near the optimum, the Fisher information the strand carries
about the **tilt** is

```
I(t) ≈ N·(dp/dt)² / [p(1−p)]  =  N·(κ−½)² / [p(1−p)]   ∝  N·(2κ−1)²
```

— high at high depth × high strand-specificity, zero at κ=½. But because `p` depends only on `t`, this
information lives **entirely in the tilt direction**. Projected onto `f_g`:

- **Single-strand node** (`f₋ ≡ 0`): `t = f₊ = 1 − f_g`, so the tilt *is* `f_g` (1:1). The strand pins
  `f_g` directly and sharply. ✅
- **AMBIG node** (both strands free): `f_g` and `t` are independent axes. Moving along the ridge
  `f_g ↔ (balanced f₊=f₋)` keeps `t` fixed and `p = ½` unchanged ⇒ **the strand likelihood is exactly
  flat in `f_g`.** It constrains the tilt but says *nothing* about the gDNA-vs-RNA split.

So at AMBIG nodes `f_g` is determined **only** by the priors/messages of §1.2–1.3. This is not a bug — it
is the genuine identifiability structure (a balanced count is equally "gDNA" or "balanced bidirectional
RNA"). It is *why* the prior precision must be honest there: nothing else is holding `f_g`.

### 1.2 The imputation messages (gDNA + per-strand RNA) — `bp_solver._message`, `fit_*_varmean`

A neighbour exports its density to a node as a Gaussian prior on the node's fraction. For source `src` →
destination `dst`:

```
μ_f   = (ρ_src·E_dst − spliced_dst) / M_dst          # identity-density mean, known spliced removed
τ_ρ   = 1 / ( σ²_bio(ρ_dst) + ρ_src/E_src )           # density-space precision
τ_f   = τ_ρ · (M_dst / E_dst)²                        # jacobian density → fraction
```

The two pieces of `τ_ρ`:
- `ρ_src/E_src` — the **Poisson sampling floor**: the source's own density is a rate estimated from a
  finite count, so it is noisy. (Count → precision, the only role the count plays.)
- `σ²_bio(ρ_dst)` — the **biological dispersion**: how much the true density genuinely varies between
  adjacent nodes at this density level. This is the learned `var~mean` curve.

**`σ²_bio(μ)` — the Poisson-offset `var~mean` (`MonotoneVarMean.fit_offset`).** Each sweep pass, over the
frozen-snapshot densities, we form one point per directed chain edge:

```
μ      = ρ(dst)                       # the mean (density level)
raw    = (ρ(dst) − ρ(src))²           # observed squared cross-node residual
offset = ρ(dst)/E(dst) + ρ(src)/E(src)  # the COMPUTED Poisson sampling floor V_p
```

The observed residual decomposes as `E[raw] = σ²_bio(μ) + V_p`. We subtract the computed `V_p` and fit the
**excess** `z = raw − offset` as a monotone-increasing P-spline of `log μ` (SCAM reparameterization
`β = cumsum([α₀, exp(α₁..)])` guarantees monotonicity; one GCV-chosen λ; robust IRLS). `predict(μ)` returns
`max(spline(log μ), 0)`, flat-extrapolated outside the fit range. gDNA fits over all edges; RNA fits
per-strand over `free_s`-continuous edges, both strands pooled.

The intent (sound): where adjacent true densities genuinely differ a lot, `σ²_bio` is large ⇒ the message
is a weak hint; where they track, it is small ⇒ a strong hint. Fitting on the *frozen* previous-pass
snapshot is meant to stop `σ²_bio` collapsing to zero within a pass.

### 1.3 The global gDNA prior — `bp_solver.global_gdna_prior` (original)

The node-class prior for AMBIG/intergenic nodes (single-strand nodes get a Beta(½,½) Jeffreys instead). A
population baseline on `f_g`:

```
ρ_global  = Σ_{count-obs} f_g·M  /  Σ_{count-obs} E      # gDNA density over count-observable regions
σ_global  = MAD over ALL regions of (f_g·M / E)          # robust between-region density spread
μ_glob[i] = clip(ρ_global · E_i / M_i , 0, 1)            # density → this node's fraction
τ_glob[i] = (M_i / E_i)² / σ_global²                     # jacobian density → fraction precision
```

"count-observable" = signature-known non-exonic (intergenic / intron) — the regions whose gDNA can be read
without exon RNA contaminating it. Recomputed each pass.

---

## 2. How the simplex solver works — `simplex_sweep`

- **Lattice.** `_simplex_lattice(K)` enumerates all `(f₊, f₋, f_g)` with `f₊+f₋ ≤ 1` on a `K`-step
  barycentric grid (`K = sweep_n_grid = 60`). Uniform over the triangle ⇒ the softmax posterior integral
  is unbiased.
- **Per-node solve** (`_local_loglik`). Evaluate `log ψ` (§0) on every lattice point: the strand mixture
  (always, every node — including AMBIG), the sided spliced lower bound, the node-class prior (Jeffreys at
  single-strand, the global at AMBIG/intergenic), and the gDNA + per-strand RNA imputation penalties.
  Normalize over the lattice.
- **Readout.** `f_g` = posterior **median** over the lattice (median, not mean — avoids the
  strand-overdispersion skew); `f₊, f₋` = posterior means; `var_g` = posterior variance of `f_g` (the
  node's confidence).
- **Sweep.** A directional Gauss-Seidel pass L→R then R→L. Each pass: freeze the snapshot densities, fit
  the gDNA and RNA `var~mean`s on them, recompute the global prior, then visit each node pulling it toward
  one neighbour's current density per direction. Iterate `sweep_max_passes`.

The solver is therefore a **precision-weighted reconciliation** of (intrinsic strand) + (neighbour
imputation) + (population prior). Its correctness rests on each `τ` being an honest inverse-variance.

---

## 3. The failure mode — over-confident `τ` overriding the truth

Both bugs are the same disease seen at two terms. The mechanism is always: **a `τ` that should be modest
becomes `1e6`–`1e11`, so its term pins `f_g` and overrides the local strand likelihood.**

### 3.1 Global prior over-confident → capture AMBIG-exon UNDER-call

Traced node `r242` (AMBIG exon under capture; balanced count `u± = 10065/10523`; oracle `f_g = 0.976`):

```
GLOB  μ=0.001  τ=7.7e6   |  gDNA msg μ=0.304 τ=2.6  |  RNA± small   ⇒  f_g = 0.000   (WRONG; truth 0.976)
```

- The node is AMBIG and balanced ⇒ the strand likelihood is **flat in `f_g`** (§1.1) — no help.
- `μ_glob = ρ_global·E/M = 0.001`: `ρ_global` is measured on un-enriched count-observable regions, so it
  is the **un-enriched** density; divided by the capture-**inflated** `M` it implies "this node is 0.1%
  gDNA." Capture-blind mean.
- `τ_glob = (M/E)²/σ_global² = 7.7e6`: `σ_global` is a **MAD**, which is robust *against* the enriched
  high-density tail — the ~1000 low-density regions set the median, the ~150 enriched exons are discarded
  as outliers ⇒ `σ_global` tiny ⇒ `τ` enormous. False confidence.
- Result: the global (τ 7.7e6) steamrolls the correct gDNA propagation message (μ 0.304, τ 2.6) and pins
  `f_g→0`. Across 92 AMBIG exons this is the **−198K** gDNA under-call = the entire capture leak.
- **Contrast** `r798` (single-strand exon, same `μ_glob≈0.001, τ_glob≈9.6e6`): its strand is *not* flat
  (`t = f₊ = 1−f_g`), so the strand likelihood overrides the global and `f_g→1.0` (correct). Single-strand
  nodes are immune; only AMBIG, with no strand defence, is exposed.

### 3.2 RNA message over-confident → zero-gDNA AMBIG PHANTOM

Traced node `r231` (AMBIG exon, zero-gDNA library; strongly stranded count `u± = 16336/149` = 99%
positive; oracle `f_g = 0`, it is +RNA):

```
GLOB μ=0.000 τ=5.6e7  |  gDNA msg μ=0.077 τ=1.2  |  RNA+ μ=0.000 τ=1.5e11  RNA− μ=0.000 τ=1.5e11
                                                                          ⇒  f_g = 1.000   (PHANTOM)
```

- The RNA messages pin **both** `f₊→0` and `f₋→0` with `τ = 1.5e11`. Since `f₊+f₋+f_g = 1`, forcing
  `f₊=f₋=0` forces `f_g = 1` — a pure-gDNA phantom in a library with **no gDNA at all**.
- This overrides the global (which is correctly pulling `f_g→0` at τ 5.6e7) **and** the local strand
  likelihood — 16,336 positive reads carry Fisher information ≈ `N(2κ−1)² ≈ 15,800`, dwarfed by `1.5e11`.
- **Where the `1.5e11` comes from (the circularity).** `τ_f = τ_ρ·(M/E)²`, `τ_ρ = 1/(σ²_bio(ρ_dst) +
  ρ_src/E_src)`. The message queries `σ²_bio` at **`ρ_dst` — the node's OWN current RNA density.** Once
  `f₊` is even slightly crushed, `ρ_dst` is small; the monotone `σ²_bio` curve is **smallest at small μ**;
  so `σ²_bio(small) → ~0` ⇒ `τ_ρ` huge ⇒ (×`(M/E)²` at the high-mass exon) ⇒ `τ_f = 1.5e11`. Crushed
  density → tiny queried variance → enormous precision → stays crushed. The frozen snapshot only delays
  this; the crushed state persists across passes. **This is the original 2c collapse**, and it is
  *self-querying at a circular density* — distinct from §3.1, which was a wrong population statistic.

### 3.3 The unifying statement

> A density-space variance (`σ_global²`, or `σ²_bio(ρ_dst)`) is **under-estimated** — by a tail-robust
> statistic on the wrong subset (global), or by a self-referential query at the node's own collapsed
> density (RNA). The `(M/E)²` jacobian (correct in itself) then magnifies the under-estimate into a
> `τ` of `1e7`–`1e11`. Because the solver is precision-weighted, that term wins and overrides the
> intrinsic strand evidence. Single-strand nodes survive (their strand likelihood is sharp in `f_g`);
> AMBIG nodes, where the strand is flat in `f_g`, take the full hit — under-calling gDNA where the global
> is wrong, over-calling it where the RNA message is wrong.

---

## 4. The fix

### 4.1 Global prior — DONE (prototype in `node_sweep`, uncommitted)

Two changes, both removing the under-estimate (the jacobian is kept — it was never the problem):

1. **Honest mean over the self-solvable nodes.** Compute `ρ_mean` as the **mass-weighted gDNA density
   over all nodes solvable without external information** — `self_solvable = ~(f₊_free & f₋_free)`, i.e.
   everything except AMBIG: single-strand (strand likelihood determines `f_g`) + intergenic (signature
   determines it). This **includes the strand-solved enriched exons**, so capture enrichment is now *in*
   the mean instead of excluded.
2. **Honest density-dependent width.** Replace the scalar tail-robust `σ_global` with a Poisson-offset
   `σ²(ρ_obs)` from `fit_offset`, keyed on the **observed total density `ρ_obs = M/E` — a known,
   non-circular axis** (not `f_g`-derived). Small at low density (un-enriched, well-understood), large at
   high density (capture-enriched, genuinely uncertain) ⇒ the prior is humble exactly where enrichment
   lives, so the propagation governs there.

Result on the capture condition (oracle dissection, deterministic): contained-gDNA error **202K → 54K
(−73%)**, the AMBIG-exon under-call eliminated, single-strand/intergenic untouched. Residual: a mild
+24K overshoot (the RNA messages, §4.2, contribute it).

### 4.2 RNA messages — PROPOSED (the analogous fix)

Same disease ⇒ same medicine: **stop the variance under-estimate; keep the jacobian.** The RNA collapse
is the *self-querying* form (§3.2), so the analog of "use a non-circular axis" is:

> **Query `σ²_bio` at a non-circular density, not at the node's own (possibly-crushed) `ρ_dst`.** The
> natural non-circular choice is the *source/imputed* density `ρ_src` (what the message is actually
> claiming), or the observed total density — exactly as the global fix keyed its width on `ρ_obs`.

Why this is enough: the over-confidence came purely from evaluating the monotone curve at a collapsed
`ρ_dst≈0` (where the curve is smallest). Evaluated at the genuine imputed density `ρ_src`, `σ²_bio` reflects
the **real** between-region RNA dispersion — which for RNA is *large* (RNA abundance varies ~10⁴-fold, the
exact reason a fraction prior is hopeless), so `τ_f` is modest, the RNA message becomes a *gentle* hint,
and the local strand likelihood (16,336 reads) correctly wins ⇒ `f₊→1, f_g→0`, phantom gone. The mean
`μ_f` may still be wrong (0 when it should be high), but a *weak* wrong hint is harmless — the strand
dominates. (If desired we can additionally floor `σ²_bio` at the genuine population RNA-dispersion, the
exact analog of the global's self-solvable spread; querying at `ρ_src` is the minimal change.)

### 4.3 Consequences (predicted; to verify on the three cached conditions)

- **Zero-gDNA AMBIG phantom:** removed — `--kill-rna-prec` already showed the phantom drops 52K→35K when
  the RNA precision stops dominating; an honest (non-zero, non-circular) `τ` achieves the same without
  deleting the messages.
- **Capture overshoot (+24K):** reduced — the overshoot is the RNA messages over-imputing on top of the
  fixed global; honest (weaker) RNA precision should pull it back toward oracle.
- **Capture under-call recovery:** preserved — it came from the *global* fix, not the RNA messages
  (confirmed: `--kill-rna-prec` on capture *returns* the under-call, so the RNA messages help but the
  global does the heavy lifting; a weaker-but-present RNA message keeps the help, drops the overshoot).
- **Risk:** RNA imputation across exon↔intron (the nascent-crossing signal) becomes a weaker prior. That
  is *correct* — RNA abundance is genuinely unpredictable between regions, so the cross-region RNA hint
  *should* be weak and yield to direct evidence (strand, spliced). The residual 35K balanced-AMBIG phantom
  (where strand is truly flat *and* counts are balanced) is a separate, smaller readout question, not
  addressed here.

### 4.4 Do we need the same fix for the gDNA imputation message?

**Same structure, much lower urgency — recommend applying it for consistency, but it is not the active
bug.** The gDNA message has the identical `τ_ρ = 1/(σ²_bio(ρ_dst) + …)` × `(M/E)²` form, so it *can*
collapse in principle. But empirically it does not bite, for two reasons:

- **Under capture**, gDNA density genuinely varies (enrichment jumps at target boundaries), so the gDNA
  `σ²_bio` stays large and the message stays honest — in the `r242` trace the gDNA message was `τ=2.6`,
  not collapsed. The enrichment heterogeneity self-protects it.
- **In a uniform-gDNA library** the gDNA densities are genuinely near-equal, so `σ²_bio→0` *can* happen —
  but there the collapse is **benign**: every neighbour agrees on the same correct density, so an
  over-confident message imputes the *right* value. (Contrast RNA, where the sweep's own smoothing
  manufactures false agreement around a *wrong*, crushed value.)

So the gDNA message is largely self-protected by the very heterogeneity that breaks the global and the
RNA. The clean engineering choice is to route **all three** through one honest-precision rule (query the
`var~mean` on a non-circular density; never let a population/neighbour variance be under-estimated) so the
collapse is structurally impossible everywhere — but if we stage it, gDNA-message is last and lowest-risk.

---

## 5. Open implementation choices (for the working session)

1. **RNA `σ²_bio` query density:** `ρ_src` (minimal, matches the message's own claim) vs `ρ_obs = M/E`
   (matches the global fix's axis) vs both-plus-a-floor (most robust). Prototype and compare on the three
   cached conditions.
2. **Whether to unify** all three precisions behind one helper (`honest_message_tau`) now, or fix RNA
   first and fold gDNA/global in once validated.
3. **The residual 35K balanced-AMBIG phantom** (strand flat *and* counts balanced, zero-gDNA): out of
   scope here; likely a readout/prior question on the genuinely-unidentifiable subset.

Harness: `scripts/debug/dissect_regions.py` (deterministic, OMP=1, cached payloads) with
`--kill-rna-prec` / `--soften-global` probes; oracle = origin-split re-scan. The global prototype lives in
`bp_solver.node_sweep` (uncommitted).
