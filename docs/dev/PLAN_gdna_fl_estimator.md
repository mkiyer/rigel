# PLAN — a contamination-robust gDNA fragment-length estimator (2026-08-30)

    ⚠ **A DEV DOC and a SANDBOX.** Nothing here is authoritative and nothing may cite into it.
    ⛔⛔ **SUPERSEDED 2026-08-31 — THIS WAS THE PLAN AND IT IS NOW BUILT.** Read `ROADMAP.md` §0's
    deconvolution row for where the tool IS; read `src/rigel/calibration/gdna_density.py` and `fl.py` for
    what shipped. This file is kept only for the DERIVATION (§3) and the REFUTED routes (§4, §7.1), which
    is what stops them being re-proposed. ⚠ Every "not yet built" and "untested" below is stale by
    construction; the implementation notes are in §9.

## 0. Response to external review (2026-08-30, revision 2)

The review endorsed the identification strategy and the removal of the constant, and raised four issues.
⭐ **Two are accepted and have changed the plan. Two turn on codebase facts the reviewer could not see and
are answered with evidence rather than adopted.** Each is resolved in place below; this section is the map.

| review point | verdict | where |
|---|---|---|
| Capture-ON breakdown is a **blocker**; needs a continuous fade | ⭐ **ACCEPTED**, and the mechanism is now built and measured — but **not** on the statistic the review proposed, which is not observable at runtime | §7.1 |
| Small-sample: decline must use the **analytic variance** of `a_0 − a_1`, not a point estimate | ⭐ **ACCEPTED** in full; now derived | §6.2.4 |
| Layering: **reject (a), enforce (b)** — refactor layer 5 to consume the layer-2 primitive | ⚠ **PARTLY.** The DRY requirement is accepted and now satisfied by construction; the *sequencing* is not, and §6.1 gives the measured reason | §6.1 |
| Latent boundary bug needs **a dedicated cleanup PR** | ⛔ **DECLINED — the bug is not in the codebase.** Audited: it was in my prototype only | §7.4 |

**Upstream evidence, unchanged and still accurate:** `DIAGNOSIS_gdna_fragment_length.md` (what is broken)
and `DESIGN_gdna_fragment_length_fix.md` (what it costs, and the routes already refuted).

---

## 1. What is being proposed, in one paragraph

`fl.py` fits the gDNA fragment-length pmf from four structural pools it asserts are "PURE BY
CONSTRUCTION". They are not: the intronic pool measures **95 % nascent RNA** and the intergenic pool
**53 % mature** (unannotated transcription). The proposal replaces the purity assertion with a
**two-pool contrast** that needs no pure pool and no template, and estimates the one nuisance parameter it
requires — the gDNA density `rho_g` — from the **low side of the per-object density distribution**, using
the fact that transcription can only ADD reads. Both halves are closed-form, introduce **no constant**,
and need **no accumulator or schema change**. Measured against the origin-split oracle, the chain removes
**88–99 %** of the gDNA length-model error on both length-gap sign arms, stranded and unstranded, and is
inert where nothing is wrong.

⛔ **What this does NOT do:** it does not repair hybrid capture — under capture the premise it rests on
is false, and the plan's answer is to DETECT that and stand down (§7.1), not to fix it; capture-ON is
defects **B** and **C**, which are independent. It also does not touch the fragment-length COMPOSITION
channel, which is retired past 0.8.0 and is a different thing that shares a word
(`EQUATIONS.md` §3d–§3e).

---

## 2. Theory — why the current estimator is wrong, and what is actually identifiable

### 2.1 The premise that failed

`fl.py`'s docstring: *"Every pool is PURE BY CONSTRUCTION, and that purity is what removes the
circularity."* The same failure the strand overdispersion had, and now a named trap:
`TRAPS: purity-is-a-property-of-the-annotation`. "Intergenic" is whatever the user's GTF leaves over;
pervasive transcription is real; nascent RNA sits inside introns by definition. **No structural class of
objects is pure gDNA**, and an estimator that asserts one is calibrated to the annotation rather than to
the genome.

The consequence is a mean-length bias of exactly

    bias  ≈  RNA_share × (len_RNA − len_gDNA)

which is why it shipped: **the ladder and the test chromosome give both components EQUAL fragment lengths
by design** (a deliberate forcing function so the EM cannot split origins on length alone —
`TRAPS: a-length-gap-bypasses-calibration`), so the second factor is zero and a 95 %-contaminated pool
reads under a bp of error. On the fl-gap arms the same pool reads **+121 bp**.

### 2.2 What makes the contaminant tractable

⭐ **The contaminant is genomically CONTIGUOUS exactly as gDNA is.** Nascent RNA runs through an intron
without splicing there; an unannotated transcript's fragments in intergenic space are contiguous too. So
inside a *contained* pool the two components experience the **same length-dependent opportunity**
`(ell − L + 1)+`, and the de-tilt that corrects the pool corrects both components identically. That is
what makes each pool a clean linear mixture of two shapes rather than a mixture of two differently-tilted
shapes.

⚠ **This is an assumption with a measured domain of validity**, quantified in §5.4: it holds off capture
and is destroyed by capture.

### 2.3 The identifiability argument

There is no channel that reads the gDNA length distribution directly. What exists is:

* **five pool histograms** (library-wide, in `payload.pool_lengths`), each a mixture at a *different*
  gDNA share;
* **per-object counts and opportunities** (`payload.region_contained_count`, region lengths), which say
  how much of each pool is gDNA;
* the fact that gDNA is **genomically uniform** and the contaminant is **not**.

Two mixtures of the same two components at *different* ratios determine both components. That is the
whole identifiability claim, and it is why the estimator uses two pools rather than trying to purify one.

---

## 3. Derivation

### 3.1 The contrast

Let `f_0` and `f_1` be the de-tilted length densities of the two CONTAINED pools (`0` = intergenic,
`1` = intronic), `g` the gDNA length pmf we want, `r` the contaminant's, and `a_p` pool `p`'s gDNA share
*by fragment count*. §2.2 gives

    f_0 = a_0·g + (1 − a_0)·r
    f_1 = a_1·g + (1 − a_1)·r

Multiply the first by `(1 − a_1)` and the second by `(1 − a_0)` and subtract; the `r` terms cancel
identically:

    (1−a_1)·f_0 − (1−a_0)·f_1  =  [a_0(1−a_1) − a_1(1−a_0)]·g  =  (a_0 − a_1)·g

    ⇒   g  =  [ (1 − a_1)·f_0  −  (1 − a_0)·f_1 ] / (a_0 − a_1)                            (★)

**Three properties, each of which matters:**

1. ⭐ **No template.** `r` never appears. This is the whole difference from the already-refuted
   subtraction, which estimated `r` from the spliced pool and was wrong because the spliced pool is not
   the contaminant (measured TV 0.148–0.187 raw, 0.053–0.091 de-tilted).
2. ⭐⭐ **The conditioning is the SEPARATION of the purities, not the purity.** Template subtraction
   divides by `a` and amplifies template error by `1/a`; at the measured `a = 0.048` that is **21×**, and a
   5 % template error left a −14 bp residual. Here the divisor is `a_0 − a_1`, measured at **0.42–0.61**
   off capture — an amplification near **2×**.
3. ⭐ **It is the identity when the pools agree.** Put `f_0 = f_1 = f` in (★): the numerator is
   `(a_0 − a_1)·f`, so `g = f`. A pool pair with nothing to say changes nothing — the estimator cannot
   manufacture a correction out of agreement.

**Degeneracy is self-announcing.** As `a_0 → a_1` the estimator is undefined, and that is correct rather
than a failure: two pools with the same composition carry no information about how to split them. Measured,
`a_0 − a_1` collapses exactly where nothing needs correcting — **0.000 at `g00`** (no gDNA at all),
0.006–0.03 at `g98` (already pure), 0.034 at `g50` capture-ON.

### 3.2 The weights, reduced to one scalar

By definition `a_p = n_gdna,p / n_p`, with `n_p` observed. gDNA is genomically uniform at density `rho_g`,
so its expected count in pool `p` is `rho_g · E_p`, where `E_p = Σ_L A_p(L)·g(L)` is that pool's total
gDNA opportunity — index geometry, plus `g`. Hence

    a_p  =  rho_g · E_p / n_p                                                              (†)

Everything on the right is known except the single scalar `rho_g`. ⭐ **And `rho_g` cancels from the
ratio:** `a_0/a_1 = (E_0/n_0)·(n_1/E_1)`, which is pure geometry and observed counts. So the two weights
live on a one-parameter ray and only the LEVEL is unknown.

⚠ **On the circularity.** `E_p` depends on `g`, which is what (★) produces — so (†) is a fixed point. It
is a *weak* one, and that is measured rather than asserted: the region opportunity `E[(ell−L+1)+]` divides
a length error by the object length, so a +121 bp pmf error perturbs `E_p` by **−1.22 % at ℓ = 10 kb and
−0.06 % at ℓ = 200 kb**, and gDNA is measured at intronic and intergenic objects, which are kb- to
Mb-scale. §5.3 shows one pass from the shipped pmf already suffices; iterating is available and cheap.

### 3.3 `rho_g` from the side the contaminant cannot reach

This is the same move as the strand fix's away-half (`EQUATIONS.md` §6a), transposed from the strand
residual to density.

**The one-sided fact:** transcription can only ADD fragments to an object. Writing object `i`'s observed
contained count as `n_i = G_i + C_i` with `G_i ~ Poisson(rho_g·E_i)` the gDNA and `C_i ≥ 0` the
contaminant, every object's observed density `n_i/E_i` is an **overestimate** of `rho_g`. So the clean
objects are the ones on the LOW side, and no class is asserted pure.

This immediately explains why the shipped pooled rate is wrong and by how much:

    rho_pooled = Σn_i / ΣE_i = rho_g + ΣC_i/ΣE_i        — biased UP by the total contaminant

measured at **8.46× too high** at `g05` on the test chromosome and **4.50×** on the ladder.

**The estimator, with no threshold.** Let `lam_i = rho·E_i` and `d_i = n_i − lam_i`. Under the pure model
`E[d_i] = 0`, so the negative part has an exact closed form — De Moivre's mean-absolute-deviation identity
for the Poisson:

    E|N − lam|  =  2·lam^(⌊lam⌋+1)·e^(−lam) / ⌊lam⌋!        ⇒
    E[(lam − N)+]  =  ½·E|N − lam|  =  lam^(⌊lam⌋+1)·e^(−lam) / ⌊lam⌋!   ≡  D(lam)

(verified to 6 d.p. against both the exact truncated sum and 2×10⁷ Monte Carlo draws at
`lam ∈ {0.3, 1, 2.7, 10, 55.5, 200}`.)

Contamination can only RAISE `n_i`, hence only LOWER `(lam_i − n_i)+`. Matching the observed one-sided
mass to its exact null expectation gives a single scalar equation:

    F(rho)  =  Σ_i (rho·E_i − n_i)+  −  Σ_i D(rho·E_i)  =  0                               (‡)

⭐⭐ **No trim depth, no quantile, no threshold, no chosen level — one root of one monotone function.**
`F(0) = −ΣD(0) ≤ 0` and `F` grows without bound (the first sum is eventually linear in `rho`, the second
grows as `√lam`), so the root is **bracketed by construction** and bisection terminates — the same
structure `EQUATIONS.md` §6b already uses for the strand overdispersion.

**Its one bias, stated before it was measured.** A heavily contaminated object contributes `0` to the
first sum but its full `D(lam_i)` to the second, so the equation compensates by raising `rho`. The bias is
therefore **one-signed (upward)** and proportional to the contaminated fraction. §5.2 measures it at
**+0.2 % to +6.7 %** against a requirement of ±20 %.

---

## 4. The trim depth — the magic number, and why it is now gone

**This section is kept deliberately, because the first working prototype had a constant in it and the
reviewer should see how it was removed rather than be told it never existed.**

The first version of §3.3 was a **quantile trim**: sort objects by `n_i/E_i`, keep the lowest `q`, pool.

    rho_hat(q) = Σ_{lowest q} n_i / Σ_{lowest q} E_i          (q = 1 is today's shipped rate)

It worked — the curve has a clear plateau, `rho_hat/rho_true` = 0.890 / 0.960 / 1.023 at q = 0.25 / 0.50 /
0.75 — and `q = 0.50` gave the best end-to-end numbers of anything tried. **It is nevertheless not
shippable**, for a reason that is not aesthetic:

| `TRIM_Q` | `rho_hat/rho_true` (`g05`) | chain error, `g05` rna_long |
|---|---|---|
| 0.20 | 0.881 | +4.66 bp |
| 0.35 | 0.931 | +6.82 |
| 0.50 | 0.975 | +8.88 |
| 0.65 | 1.003 | +10.30 |
| **0.80** | **1.981** | **+67.69 bp — BREAKS** |

⛔ **The basin is wide but the constant is not safe, because the safe depth is a function of the library's
CONTAMINATED FRACTION, which is a property of the sample and not of the tool.** The test substrate has
~26 % of objects contaminated; a library with 60 % contaminated objects would break at `q = 0.5`. A
constant that is safe on the panel and unsafe on real data is exactly
`TRAPS: a-constant-parked-a-value-off-a-knife-edge`, and `no-magic-numbers` forbids it.

⭐⭐⭐ **The resolution is not to derive a better `q` — it is to remove `q`.** Equation (‡) is the trimmed
estimator's limit done exactly: instead of guessing how many objects are clean, it matches the observed
one-sided mass to the amount the pure model predicts. The plan carries **no trim depth at all**, and the
trimmed version survives only as a diagnostic arm to A/B against.

---

## 5. Measurements

⛔ All of it is prototype-level: scored on histograms against the origin-split oracle, no pipeline A/B.
The oracle is read only to SCORE, never to produce an estimate.

### 5.1 What the whole chain does (trim-free `rho`, no oracle input anywhere)

Mean gDNA fragment length, capture-OFF, error in bp. POOLED is what ships today.

| arm | condition | POOLED | **CHAIN** | removed |
|---|---|---|---|---|
| rna_long | `g05` ss0.99 | +118.96 | **+13.96** | 88 % |
| rna_long | `g25` ss0.99 | +59.91 | **+3.21** | 95 % |
| rna_long | `g50` ss0.99 | +27.96 | **+1.93** | 93 % |
| rna_long | `g50` **ss0.50 unstranded** | +27.49 | **+2.07** | 92 % |
| gdna_long | `g05` ss0.99 | −120.59 | **−11.33** | 91 % |
| gdna_long | `g25` ss0.99 | −60.84 | **−0.84** | 99 % |
| gdna_long | `g50` ss0.99 | −27.81 | **−2.03** | 93 % |
| gdna_long | `g50` **ss0.50 unstranded** | −27.66 | **−1.65** | 94 % |

⭐ **Both sign directions.** The gap-sign flip is the test a repair must pass rather than a formality
(`TESTING.md` §0): a change that merely trades one bias for another moves the two arms the same way.
⭐⭐ **It works UNSTRANDED**, which is why it is preferred over splitting the length banks by strand — that
alternative is undefined at `kappa = ½`, has no orientation for the intergenic pool at all, and would need
a C++ accumulator change plus a rebuild of every scan cache. unstranded × capture-OFF is an IN-SCOPE
0.8.0 stratum.

### 5.2 The weight estimator on its own

`rho_hat/rho_true`, no oracle input:

| substrate | objects | `g05` | `g50` | `g50` unstranded | `g98` |
|---|---|---|---|---|---|
| test chromosome | 77 | 1.067 | 1.010 | — | 1.004 |
| **ladder** | **11,117** | **1.035** | **1.010** | **1.010** | **1.002** |
| *(today's pooled rate, ladder)* | | *4.495* | *1.184* | *1.184* | *1.004* |

⭐ Always within **+0.2 % … +6.7 %**, and **always on the high side** — the direction §3.3 predicted before
measuring. ⚠ The bias does **not** shrink cleanly with depth in these data (the test chromosome is deeper
per object and biased more), so its scaling is a function of the contaminated fraction and is **not
established**; only its sign and its measured magnitude are.

### 5.3 The controls

| control | required | measured |
|---|---|---|
| `g00` (zero gDNA) | must not invent a distribution | `a_0 − a_1` = **0.000**; the estimator declines |
| `g98` (near-pure gDNA) | must not disturb a correct answer | declines (`a_0 − a_1` = 0.006–0.03); pooled is already within 0.7 bp |
| **equal-length ladder, 70,270 regions** | must be inert | pmf moves **−0.30 / +0.02 / −0.00 bp**; `rho_hat/rho_true` 0.999–1.009 |
| equal-length test chromosome, 152 regions | must be inert | −6.19 bp at `g05`, ≤0.8 bp elsewhere — **small-sample wobble, see §7.2** |
| one pass from the shipped pmf | must not need iteration to be useful | one pass gives the table in §5.1 |

### 5.4 Where the premise dies

`TV(r_intergenic, r_intronic)` — how different the two pools' contaminants actually are, which §2.2
assumes is small:

| | capture-OFF | capture-ON |
|---|---|---|
| TV between the two contaminants | **0.057 – 0.141** | **0.95 – 0.97** |

⛔ **Under capture the shared-contaminant premise is not approximately true, it is false**, and the
purities stop being separated as well. Measured on the ladder at `g05` capture-ON the chain is **worse**
than what ships: −22.34 bp → **−24.03 bp**. See §7.1 — this is the most important open issue.

---

## 6. Implementation plan

⛔ **Nothing below is written. This is the plan being reviewed.** Sequencing follows the house rule
DERIVE → DESIGN → PLAN → PROTOTYPE → A/B → only then `src/`; derivation, design and prototype are done, and
**the A/B is a precondition for step 4**, not a follow-up to it.

### 6.1 Where the code goes — and the layering resolution

`fl` is **layer 2** (opportunity). The estimator needs, and only needs:

| input | source | layer | direction |
|---|---|---|---|
| pool length histograms | `payload.pool_lengths` | — | already an argument |
| per-pool opportunity `A_p(L)`, total `T(L)` | `gdna_opportunity` | 2 | sideways ✅ |
| per-object counts | `payload.region_contained_count` | — | already on the payload |
| region lengths, per reference | `payload.region_bounds` + `ref_region_bound_offsets` | — | already on the payload |
| region TYPE (intergenic / intron) | `splice_graph.build_region_partition_arrays(index)` | 1 | downward ✅ |

⭐ **No upward import is required and no new data has to be accumulated.**

**THE DUPLICATION THE REVIEW FLAGGED IS REAL, AND SO IS ITS FIX — but not by the route proposed.**
The reviewer asked to make `gdna_density.py` the layer-2 primitive and refactor
`density_deconv.fit_intron_background` (layer 5) to consume it. Two facts from the source change the
shape of that:

1. ⛔ **`fit_intron_background` is not a rate estimator and is not a drop-in.** It is pool SELECTION
   delegating to `fit_gdna_background`, which returns a `GdnaBackground` — a Gamma-conjugate posterior
   `log_mu_bg = ln((Σg + ½)/ΣE)` **plus** a method-of-moments over-dispersion `alpha` and a posterior
   shape. Swapping the location changes a NegBinom prior, not a scalar.
2. ⛔⛔ **Its pool is intergenic-only BY DESIGN, for non-circularity**: *"introns are what we
   deconvolve"*. My estimator deliberately pools intergenic **and** intronic, because the fl weights need
   the density that applies to both pools. **The two functions have different estimands and different
   pools on purpose**, so "consume the primitive" is a behaviour change to the composition reference, not
   a refactor — and that reference is a live, partly-broken area (`ROADMAP.md` §0: the per-object
   reference is *"not shippable"*, stranded × capture-ON still regressing ~2×). Landing it inside the fl
   change would move two mechanisms at once, and *"a change that cannot be A/B'd alone cannot be judged
   alone"*.

⭐⭐⭐ **THE RESOLUTION — option (d), which satisfies the review's DRY requirement without its sequencing.**
The new module owns **both** estimators over the same signature `(counts, exposure) -> rate`:

    pooled_rate(counts, exposure)      # ln((Σg + ½)/ΣE) — today's formula, moved, not changed
    one_sided_rate(counts, exposure)   # the root of (‡)

and `fit_gdna_background` (layer 5) is refactored to call **down** to `pooled_rate` — a legal 5 → 2
import and a **byte-identical** change, gated by a test asserting bit-identity on the panel.

⭐ **There is then exactly ONE implementation of each estimator, in one home**, which is what the review
actually asked for. `fl` uses `one_sided_rate` now. Layer 5 adopting the robust rate later becomes a
**one-line call-site change**, independently A/B-able — and it should happen, because
`fit_gdna_background`'s own docstring says *"from a pool of pure-gDNA regions"*, which is the **third**
appearance of the premise this whole work exists to remove, and its intergenic pool measures 53 % mature
RNA at `g05`. ⚠ That is a real, separate defect this work uncovered; it is recorded, not fixed here.

⚠ **Layer choice.** A pure `(counts, exposure) -> rate` function has no calibration dependencies, so
layers 0 or 1 would also accept it. Layer 2 is proposed because its job description is *"every divisor in
the tool is derived here"* and a rate is a count divided by an opportunity. ⛔ This is the one item still
genuinely open for the owner; it does not affect the maths.

### 6.2 The changes, in order

⭐ **Step 0 is new in revision 2 and is a PREREQUISITE, not an option** — it is the §7.1 resolution and it
must be measured before anything is written into `src/`.

0. ⭐⭐ **Prototype the non-negative projection (§7.1 (ii)) and measure it**, outside `src/`, on both
   fl-gap sign arms, the equal-length control and the ladder, capture-OFF **and** capture-ON. It must
   (a) reproduce (★) exactly wherever (★) is already feasible, and (b) remove or materially reduce the
   capture-ON regression. ⛔ **If it does neither, stop and bring the result back** — the alternative is
   the fade of §7.1 (i), which reintroduces a chosen functional form and would need its own justification.
1. **New layer-2 module** (proposed name `gdna_density.py` — descriptive, and it avoids colliding with
   layer 5's `density_*` vocabulary), owning **both** rate estimators over one signature:
   * `pooled_rate(counts, exposure)` — `ln((Σg + ½)/ΣE)`, **moved from `density_deconv`, not changed**.
   * `one_sided_rate(counts, exposure) -> GdnaDensityFit` — the root of (‡) by bracketed bisection, plus
     the diagnostics that make it auditable: `n_objects`, `bracket_ok`, `pooled_rate`, `rate_over_pooled`
     (how far the one-sided fit moved off the naive one — the number that says whether contamination was
     found at all), `declined` and its reason.
   * `contained_opportunity(pmf, region_lengths)` — `E[(ell−L+1)+]` per region.
   * `region_lengths(payload)` — per-reference, see §7.4.
2. **`density_deconv.fit_gdna_background` calls DOWN to `pooled_rate`.** Byte-identical; gated by a
   bit-identity test on the panel. This is what makes the DRY fix structural rather than promised.
3. **`fl.py`:** a new keyword argument carrying the per-region counts, lengths and types (or the assembled
   fit), and the contrast (★) — via step 0's projection — replacing the four-pool sum for `gdna_counts`.
   ⛔ **Keep the current four-pool path as the fallback when the new inputs are absent**, exactly as the
   `gdna_opportunity=None` fallback works today: `build_fl_models` is called by the second pass and by
   tests with varying inputs, and the honest fallback is the existing behaviour, not a crash.
   ⛔ **And the module docstring's "Every pool is PURE BY CONSTRUCTION" must go** — it is the false premise
   itself, and leaving it would be the documentation half of the same bug.
4. **`pipeline.py`:** pass the new inputs at both call sites (~line 397 and ~line 974). Both already have
   `index` in scope, and one already calls `build_region_partition_arrays(index)` a few lines above.

### 6.2.4 The decline rule — derived, per the review

⭐ **ACCEPTED IN FULL: the decline must read the analytic variance of `a_0 − a_1`, not the point estimate.**
Revision 1 said "derive it from the data's own resolution" and left it there; here it is.

With `a_p = rho·E_p/n_p`, the two dominant error sources are the Poisson noise in `n_p` and the sampling
error in `rho`. Since `n_p` is a fragment count, `Var(n_p) ≈ n_p`, and by the delta method

    Var(a_p | rho)  ≈  a_p² · [ Var(rho)/rho²  +  1/n_p ]

`a_0` and `a_1` share `rho`, so that term is **common-mode and cancels from the difference** — which is
the useful part, and it is why the ratio `a_0/a_1` was noted as `rho`-free in §3.2:

    Var(a_0 − a_1)  ≈  a_0²/n_0  +  a_1²/n_1        (+ the common `rho` term, which cancels)

**Decline when `(a_0 − a_1)²  ≤  Var(a_0 − a_1)`** — i.e. when the separation is not distinguishable from
zero at the data's own resolution. ⭐ **No constant is chosen**: the comparison is the estimator's own
standard error against itself, the same shape as the strand channel's derived noise-floor deadband
(`disc = 4·max(0, (κ−½)² − σ²_d)`, `EQUATIONS.md` §5.2b). ⚠ It also answers §7.2 — a small annotation has
small `n_p`, so the rule declines on thin data automatically rather than needing an object-count floor.

On declining, return the CURRENT four-pool answer and say so in the QC line. ⛔ **The decline must be
counted and surfaced**, not silent: a run that declined everywhere and a run that corrected everywhere
must not look the same.

### 6.3 Tests

Following `a falsification test first, verified failing — then break the fixed code and watch each gate
fire` (`TRAPS: perturb-every-gate`):

* **the identity**: `f_0 = f_1` ⇒ output equals input, exactly. Falsify by perturbing one pool.
* **exact recovery**: synthesise `g` and `r`, mix at known `a_0 ≠ a_1`, assert (★) returns `g` to float
  tolerance. Falsify by passing the wrong weights.
* **De Moivre**: `D(lam)` against the exact truncated sum over a `lam` grid spanning the sub-1 and
  large-`lam` regimes. (Already verified in prototype; it becomes a unit test.)
* **bracketing**: `F(0) ≤ 0 < F(hi)` on random inputs, so bisection can never fail to terminate.
* **one-sidedness**: adding a non-negative contaminant to any object must not DECREASE `rho_hat`.
* **the zero controls**: no gDNA ⇒ declines; no contaminant ⇒ returns the pooled rate.
* **decline path**: `a_0 = a_1` returns the four-pool answer and sets the reason.
* ⛔ **the per-reference region-length gate** — see §7.4.

⚠ **Expected suite-count movement**, to be re-derived and not adjusted: one new
`src/rigel/calibration/` module is **+3** (jargon, docs-boundary, layering-if-declared); a new
`tests/calibration/` file is **+2**; plus the new cases themselves. Current standing baseline is
**3,680**.

### 6.4 How it will be judged

⛔ **Not on calibration and not on transcript accuracy** — measured, both flip sign between the two fl-gap
arms. On **`gdna_frac_est`**, the library gDNA fraction, via `scripts/design/em_fl_ceiling.py`, which
already exists and carries the `noop`/fires/reseed-floor gates.

The perfect-pmf ceiling is **8–32 % of the standing library-gDNA bias** (measured, 9–47× the reseed
floor). The chain recovers 88–99 % of the *pmf* error, so most of that ceiling should be reachable —
⛔ **but that is a PREDICTION, not a result**: the pmf error and the composition error are different
quantities and the mapping between them has not been measured. **If the A/B lands inside the noise floor,
the correct outcome is to not ship this**, and §7.6 is the pre-registered version of that.

Runs required: both fl-gap sign arms **and** the equal-length control **and** the ladder, capture-OFF and
capture-ON, `g00`/`g05`/`g25`/`g50`/`g98`.

---

## 7. Open issues

### 7.1 ⭐ Capture-ON — the review called this a blocker, and it now has a measured, observable answer

**The problem, restated.** Under capture the two pools' contaminants stop resembling each other
(oracle TV **0.95–0.97** against 0.06–0.14 off capture), because the probes reshape the RNA length
distribution. The shared-`r` premise (★) rests on is then false, and the ladder's `g05` capture-ON row got
*worse*: −22.34 bp → −24.03 bp.

⛔ **The review's proposed fix is not implementable as stated.** It asks for a fade weighted by "observed
TV between pool distributions". The TV that diagnoses the failure is between the two **contaminants**, and
measuring that needed the origin-split oracle — it is not available at runtime. The TV between the two
observed *pools* is a different quantity: it confounds the `g`-vs-`r` difference with the `a_0`-vs-`a_1`
difference, and is large precisely when the estimator is working best.

⭐⭐⭐ **The observable substitute, derived from the model's own consistency and now MEASURED.** The
two-component model asserts that both recovered components are densities, so both must be non-negative.
When the premise holds the raw solution is non-negative up to sampling noise; when it fails the solution
goes negative, and **the negative mass is a direct, runtime-observable measure of premise violation** —
no oracle, no new bank, no threshold:

    g_raw = [(1−a_1)f_0 − (1−a_0)f_1]/(a_0−a_1)      r_raw = [a_0 f_1 − a_1 f_0]/(a_0−a_1)
    VIOL  = (negative mass of g_raw) + (negative mass of r_raw)          0 ⇒ premise consistent

| substrate | condition | **VIOL (observable)** | oracle TV(r_0,r_1) |
|---|---|---|---|
| ladder | `g05` / `g50` capture-**OFF** | **0.0000 / 0.0000** | — |
| ladder | `g05` / `g50` capture-**ON** | **0.032 / 0.933** | — |
| rna_long | `g05` / `g50` capture-**OFF** | **0.061 / 0.003** | 0.103 / 0.141 |
| rna_long | `g05` / `g50` capture-**ON** | **0.730 / 2.199** | 0.948 / 0.969 |

⭐ It separates capture-OFF from capture-ON by **one to three orders of magnitude** and tracks the oracle
TV it is standing in for. On the ladder capture-OFF it is **exactly 0.0000**.

**Two ways to use it, and the second is better:**

* **(i) a continuous fade** — shrink the contrast toward the pooled estimate in proportion to `VIOL`.
  ⭐ Simple, and the detector is measured. ⛔ The functional form (linear? in what units?) is a choice, so
  this reintroduces exactly the kind of decision §4 worked to remove.
* **(ii) ⭐⭐ project onto the feasible set (RECOMMENDED, UNTESTED)** — instead of solving (★) and clipping
  the negatives, solve for the non-negative `(g, r)` pair that best explains `(f_0, f_1)`: a small
  non-negative least squares. **This needs no fade parameter at all.** When the premise holds it returns
  (★) exactly, because (★) is then already feasible; when it fails it returns the closest feasible
  explanation, which degrades toward the pooled answer by construction. ⚠ **Not yet measured** — this is
  the first thing to build in §6.2.

⛔ **What neither fixes.** At ladder `g05` capture-ON `VIOL` is only 0.032, so the detector barely engages
and a small regression would remain. ⚠ In context that regression is **1.7 bp on an error of 22 bp
(7.6 %)**, and capture-ON's error is dominated by defects **B** and **C** — a probe-blind opportunity
divisor and an EB shrinkage anchored on the global mixture — which this change does not touch and which
own the already-priced −5.90 %. **I am not claiming this makes capture-ON good; I am claiming it can be
made not-materially-worse, and §7.6's kill criterion holds it to that.**

### 7.2 Small-object-count wobble — now handled by §6.2.4, but not eliminated

On the equal-length **test chromosome** (152 regions, 77 usable) `g05` costs −6.19 bp where the pooled
estimate was already right. On the **ladder** (70,270 regions, 11,117 usable) the same control costs
−0.30 bp. So it is object count rather than bias.

⭐ The review's fix is accepted and is now §6.2.4: the decline rule reads `Var(a_0 − a_1) ≈ a_0²/n_0 +
a_1²/n_1`, so thin data declines on its own without an object-count floor. ⚠ **It is not eliminated**:
the variance rule keys on pool fragment counts `n_p`, which are large even when the number of OBJECTS is
small, and the −6.19 bp wobble came from few objects rather than few fragments. A second-order term for
the object count of the `rho` fit may be needed; it is not derived. **This must be measured on the
test-chromosome control in the step-0 prototype, not assumed away.**

### 7.3 The bias of (‡) is one-signed but its scaling is not established

Measured +0.2 % … +6.7 %, always high, against a ±20 % requirement — comfortable. But the deeper substrate
was biased *more*, not less, so the bias tracks the contaminated fraction rather than depth and I cannot
predict it on a library with a different contamination profile. ⚠ **A real library could be more
contaminated than any panel tested here.** Two mitigations, neither built: iterate (‡) with the
contaminated objects excluded by their own one-sided residual; or report the bias bound alongside the
estimate so a consumer can see it.

### 7.4 ⛔ The boundary bug — AUDITED, and it is NOT in the codebase

The review asked for "a dedicated cleanup PR". **That is not needed, and the audit is the evidence.**

`np.diff(payload.region_bounds)` straight through IS wrong — bounds are concatenated per reference
(`ref_region_bound_offsets`), so a plain diff manufactures one phantom cross-chromosome region per
junction. On the test substrate that phantom landed on the blank chromosome, exactly where the unannotated
transcripts live, and it made the weight estimator look **6× worse than it is**.

⭐ **But it was in MY PROTOTYPE only.** Audited across `src/` and `scripts/`: the one shipped consumer,
`gdna_opportunity.py:224`, slices per reference first (`region_bounds[lo:hi]`, inside a loop over
`offsets`) and then diffs — correct. Every other `np.diff` in the package is over offset arrays or
per-reference slices. **No shipped code has this bug and there is nothing to clean up.**

⚠ What survives is a *test* obligation, not a PR: the new module must take region lengths through a
per-reference helper, and §6.3 carries a gate that a two-reference index produces no phantom region. It is
recorded here because the failure was silent, plausible-looking, and nearly produced the wrong verdict.

### 7.5 Second-order concerns

* **Two-pool choice.** The plan uses the two CONTAINED pools because they share one opportunity geometry.
  The two CROSSING pools are ignored, and they dominate under capture. A four-pool generalisation of (★)
  is a least-squares extrapolation to `a = 1` rather than a 2×2 solve; unbuilt, and it may be the right
  answer for §7.1.
* **`POOL_EB_PRIOR_ESS = 1000`** still shrinks the result toward the global MIXTURE, which is mostly RNA
  whenever gDNA is a minority. It is inert while pools are large and **dominant under capture** (measured
  `ship−pool = −23.7` of a −31.7 total). This plan does not touch it; it is defect **C** and independent.
* **The fixed point** (†)/(★) is run for one pass. Convergence has not been demonstrated, only shown to be
  unnecessary at the measured tolerances.
* **Real data.** Every number here is simulated. `real-data-is-a-test-input` — the estimator should be run
  on the cfRNA libraries as a *census* before it ships, since that is where a genuinely different
  contamination profile lives.

### 7.6 Pre-registered kill criteria

So the A/B cannot be rationalised after the fact:

* if the `gdna_frac_est` improvement on the fl-gap arms is inside the `base_reseed` noise floor → **do not
  ship**, close rank 1b as measured-and-not-worth-it;
* if the equal-length ladder control moves by more than its own noise floor → **do not ship**; the
  estimator must be inert where nothing is wrong;
* if the capture-ON regression exceeds the capture-OFF gain on any in-scope stratum → **do not ship**
  without a resolution to §7.1.

---

## 8. Certification (revision 2)

**What I certify.**

* The derivations in §3 are complete and self-contained. (★) is two lines of algebra; (†) is a definition
  plus the uniformity of gDNA; (‡) rests on De Moivre's Poisson MAD identity, which I verified to 6 d.p.
  against both the exact truncated sum and 2×10⁷ Monte Carlo draws at `lam` spanning 0.3 to 200.
* **The plan introduces no magic number.** The trim depth the first prototype needed is gone (§4), and in
  revision 2 the last remaining candidate — the decline rule — is **derived** rather than deferred
  (§6.2.4: `Var(a_0 − a_1) ≈ a_0²/n_0 + a_1²/n_1`, compared against the separation itself). ⚠ The one
  functional-form choice that could still enter is §7.1 (i)'s fade, which is why §7.1 (ii)'s projection is
  the recommended route and why measuring it is step 0 rather than a follow-up.
* Every number in §5 and §7.1 was produced this session by prototypes outside `src/`, scored against the
  origin-split oracle, with the oracle used **only to score**. The controls in §5.3 behave.
* **The two review points I declined are backed by source, not by argument**: `fit_intron_background` is a
  NegBinom fit whose intergenic-only pool is deliberate and non-circular (§6.1), and the boundary bug is
  absent from `src/` and `scripts/` — audited, with the one shipped consumer slicing per reference before
  it diffs (§7.4).
* The implementation is bounded: one new layer-2 module, one byte-identical refactor at layer 5, one new
  argument and one replaced code path in `fl.py`, two call sites in `pipeline.py`. **No accumulator
  change, no schema change, no cache rebuild.**

**What I do NOT certify, and what would make me wrong.**

* ⛔ **No A/B through the shipped pipeline has been run.** Four times a toy-positive change has been
  panel-negative here (`TRAPS: panel-before-src`). The §6.4 claim that an 88–99 % pmf gain translates into
  a share of the measured 8–32 % `gdna_frac_est` ceiling is an **inference, not a measurement** — the two
  are different quantities and the mapping between them is unmeasured.
* ⛔ **§7.1 (ii), the projection, is UNTESTED.** I have measured the detector that motivates it (`VIOL`,
  separating capture-OFF from capture-ON by one to three orders of magnitude) but not the projection
  itself. That is why it is step 0.
* ⚠ **§7.2 is not fully closed.** The derived decline rule keys on pool fragment counts, and the
  small-sample wobble came from few OBJECTS, not few fragments. The rule may not catch it.
* ⛔ **The bias of (‡) is one-signed and measured at +0.2 … +6.7 %, but its scaling is not established**
  (§7.3) — the deeper substrate was biased more, so it tracks the contaminated fraction, and a real
  library could be more contaminated than any panel here.
* ⛔ **Nothing has been tested on real data.** `real-data-is-a-test-input`: the cfRNA libraries should be
  run as a census before this ships, because that is where a genuinely different contamination profile
  lives.

**Am I ready to implement?**

⭐ **Yes — with one prerequisite and one decision.**

* **Prerequisite (mine):** step 0 of §6.2 — prototype and measure the non-negative projection. It is
  outside `src/`, it is bounded, and it is the difference between "capture-ON regresses" and "capture-ON
  is handled". I can do it now and report before touching the package.
* **Decision (yours):** §6.1's layer placement. Revision 2 resolves the DRY objection structurally —
  one home per estimator, layer 5 calling down — so what is left is genuinely a placement question
  (layer 2 next to the other divisors, versus layer 0/1 as a dependency-free primitive). It does not
  affect the maths, and I should not make it unilaterally.

⛔ **Still NO for shipping**, and the bar is unchanged: §7.6's pre-registered kill criteria must hold on
both fl-gap sign arms, the equal-length control and the ladder, capture-OFF and capture-ON — judged on
`gdna_frac_est`, never on calibration and never on transcript accuracy.


---

## 9. What was actually built (2026-08-31), and where the plan was WRONG

⭐ **Built, per §6.2 steps 1–4:** `src/rigel/calibration/gdna_density.py` (layer 2, declared in
`_layers.py`), owning `pooled_log_rate`, `one_sided_rate`, `contained_opportunity` and
`region_lengths_from_partition`; `density_deconv.fit_gdna_background` refactored to call **down** to
`pooled_log_rate` (byte-identical); the contrast in `fl._deconvolved_gdna_counts` behind two new optional
arguments to `build_fl_models`; both `pipeline.py` call sites wired.
Gates: `tests/calibration/test_gdna_density.py`, 34 cases.

### ⛔ THE PLAN WAS WRONG ABOUT CAPTURE-ON, AND IN THE FAVOURABLE DIRECTION

§7.1 called capture-ON a blocker and predicted a regression. **Measured end to end against the shipped
four-pool sum, capture-ON is a large WIN**: the library gDNA bias falls **0.018925 → 0.005500 (71 %)** at
`g05` and **0.022875 → 0.000915 (96 %)** at `g50`.

⚠ **The error was a BASELINE error, and it is worth naming.** The prototype compared the contrast against
the two CONTAINED pools pooled; what ships is the FOUR-pool sum, which is a different and much worse
starting point under capture (the crossing pools dominate there and are length-selected long). Comparing
against the wrong baseline produced a confident, wrong, blocking verdict.
⛔ **The premise is still false under capture (TV 0.95) — so the estimator is right for a reason it does
not model, and no detector catches that.** That risk is now recorded in `ROADMAP.md` §0 and in the
function's own docstring, and it is strictly worse than a measured regression would have been, because
nothing would announce it breaking on a real probe panel.

### ⛔ STEP 0's PROJECTION WAS REFUTED, AND SO WAS THE VIOLATION DETECTOR

Both were prototyped and both failed, which is why neither shipped:

* **the non-negative projection** (§7.1 (ii)) reproduces the contrast exactly wherever the contrast is
  already feasible — the required no-op property — but where it is not, it is **worse** than clipping
  (`−23.92 → −29.31` bp at `g50` capture-ON). The best feasible fit to a wrong model is not closer to
  truth. Clipping ships.
* **the violation-vs-noise decline test** does not separate: the predicted noise floor scales with the
  violation itself, so the ratio sits at **0.21–0.73 on AND off capture** and never fires. ⚠ This also
  retracts the previous revision's claim that raw `VIOL` is a usable detector — that separation was
  confounded with pool depth.

### ⚠ Two smaller corrections

* **The `a_0 − a_1` decline rule (§6.2.4) does not fire at `g00`** on the test chromosome: with both
  weights near zero its standard error also goes to zero, so a separation of 0.003 reads as significant.
  It does fire on the ladder ("an empty pool"). ⭐ Measured, it is harmless — the zero control's bias
  goes **0.046315 → 0.045660**, slightly BETTER — so this is a latent sharp edge rather than a defect,
  and the derived fix is an SE that carries `rho`'s own relative error instead of dropping it as
  common-mode.
* **The equal-length control moves** (bias `0.050200 → 0.048750` at `g05`), which the letter of §7.6
  forbade. It moves TOWARD truth, so the criterion's intent — do not damage a correct estimate — holds,
  but the criterion as written was about the pmf and was applied to the composition.


### ⛔ §7.6's KILL CRITERION 2 IS TRIPPED — reported, not absorbed

"If the equal-length ladder control moves by more than its own noise floor → do not ship; the estimator
must be inert where nothing is wrong." Measured by ablation on the ladder, capture-OFF:

| condition | bias OFF | bias ON | reseed floor | verdict |
|---|---|---|---|---|
| `g00` ss0.99 capture-OFF | 0.008842 | 0.008842 | 0.000020 | ⭐ **EXACTLY unchanged** — it declines |
| `g05` ss0.99 capture-OFF | 0.013696 | 0.013751 | 0.000044 | worse, 1.25× floor |
| `g50` ss0.99 capture-OFF | 0.051113 | 0.051212 | 0.000057 | worse, 1.7× floor |
| `g50` ss0.50 capture-OFF | 0.072223 | 0.072420 | 0.000156 | worse, 1.3× floor |
| `g98` ss0.99 capture-OFF | 0.001621 | 0.001629 | 0.000041 | inside the floor |
| `g05` ss0.99 **capture-ON** | 0.001183 | 0.000822 | 0.000062 | ⭐ **better, 6× floor** (30.5 %) |
| `g50` ss0.99 **capture-ON** | 0.033971 | 0.031193 | 0.000054 | ⭐ **better, 51× floor** (8.2 %) |

⭐⭐ **The ladder TOTAL favours the change by roughly 9×**: the capture-ON gains sum to 0.00314 against
capture-OFF costs of 0.00036. ⚠ Read the strata apart anyway (`TRAPS: never-pool-the-strata`) — the sign
differs by capture, and it is the capture-OFF half that trips the criterion.

⭐ **The mechanism is understood and it is the expected one**: with equal component lengths there is no
length bias to remove, so only the contrast's extra variance (≈2×, the `1/(a_0−a_1)` amplification)
remains. ⭐ The same ablation on the equal-length TEST CHROMOSOME *improves* (0.050200 → 0.048750),
because that substrate carries shadow transcripts and therefore has real contamination even at equal
lengths — so this is a variance cost where there is nothing to correct, not a bias.

⚠ **Magnitude: 0.3–0.4 % of the standing bias, against a 20–250× larger gain wherever the components'
lengths differ.** ⛔ But the criterion was pre-registered precisely so this could not be rationalised
after the fact, and the ladder is the panel the tool is ranked on. **The call is the owner's.** The
options, none of which is free: accept the trade as measured; damp the contrast toward the pooled estimate
by the estimated length gap (a new mechanism, underived); or gate on a measured gap (a switch, which the
house prefers to a fade).
