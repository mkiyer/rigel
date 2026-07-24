# The GRAFT message: mode and precision — derivation, validation, and implementation plan

**Status:** derived, **Monte-Carlo validated**, ready to implement. **Date:** 2026-07-23.
**Scope:** the RNA + gDNA message a splice-junction **boundary** sends to its **exon** flank (`boundary → exon`).
**Not in scope:** the `exon → boundary` PEEL, which is a *difference* and is a genuinely different problem —
see §9 and `peel_message_precision_open.md` (to be written; the peel is **not** derived).
**Validation harness:** `scripts/debug/message_precision_mc.py` (seeded, reproducible).
**Supersedes:** `prediction_measurement_merge_derivation.md` (item E) — that derivation's *algebra* is correct
and is reproduced here as T2, but its §5 "structural guarantee" was justified with an `exon → boundary`
failure narrative that cannot occur on a graft edge (§9.1). `message_arithmetic_reconciliation.md` Part D/E
contains the same conflation and should be read only for its measured tables.

---

## 0. Grounding — RNA is RNA (owner, 2026-07-23)

There is **one RNA species**. "Mature" and "nascent" are not two entities to track; they are two **routes**
the same RNA takes at a junction:

* it **splices out** — jumps to another exon; observed in the **spliced** channel;
* it **continues contiguously** — observed **unspliced** in the next genomic region.

Pass-0 solves the **unspliced** pool. The spliced count is a direct measurement of the departing flux. All of
the arithmetic below is bookkeeping on "where did the RNA go at this junction", nothing more.

Capture multiplies the observed density of **every** component at position `x` by the same efficiency `e(x)`
(nucleic-acid-agnostic at pass-0).

---

## 1. Quantities at one junction

Boundary `b` (a point — crossing flux), exon `x` (a region — contained mass).

| symbol | definition | provenance |
|---|---|---|
| `M_b`, `n_b` | boundary unspliced crossing mass / **integer** count | observed |
| `S_b`, `n_s` | boundary spliced mass / **integer** count | observed |
| `f_g` | boundary's solved gDNA fraction of its **unspliced crossing** | solved (has `Var(f_g)`) |
| `ρ_g = f_g·M_b/E_g^b` | gDNA density crossing `b` | imputed |
| `ρ_ν = (1−f_g)·M_b/E_r^b` | RNA that **continues** through `b` | imputed |
| `ρ_μ = S_b/E_spl^b` | RNA that **splices out** at `b` | measured |
| `ρ_R = ρ_ν + ρ_μ` | **total RNA flux** through the junction | mixed |
| `w_μ = ρ_μ/ρ_R`, `w_ν = ρ_ν/ρ_R` | provenance shares, `w_μ + w_ν = 1` | — |

**The exon receives the whole flux.** A transcript passing through the exon deposits contained-unspliced
fragments there *regardless of which route it takes at the junction*. So the exon's component set is
`{gDNA, RNA}` and the boundary's — **once the spliced is included** — is the same. That component-set match
is what makes the graft clean, and it is the entire reason the graft is easy and the peel is hard.

---

## 2. The mode — transport `k`, never `f_g`

```
k          = ρ_g / ρ_R                                    the enrichment-free content ratio
f_g(x)     = k·E_g^x / (k·E_g^x + E_r^x)                  re-formed in the EXON's own frame
f_R(x)     = 1 − f_g(x)                                   by construction
```

`k` is `enrichment_frame.k_from_belief`; the re-forming is `f_g_from_k`. **Do not copy `f_g` across** — the
two nodes have different effective lengths (contained vs crossing, and gDNA-FL vs RNA-FL), so `f_g` is
frame-dependent while `k` is not.

### T1 — enrichment cancels **identically**

Capture scales every density at `b` by `e(b)`. `k` is a *ratio of two densities at the same node*, so `e(b)`
cancels exactly; and `f_g(x)` depends only on `k` and the **exon's own** effective lengths. Therefore the
graft message is **enrichment-free** and needs no reframe ratio `r` at all.

> **MC: max |Δf_g(x)| = 0.00e+00 for e ∈ [10⁻³, 10⁶].** Exact to machine precision.

This is why `Var(log r)` (design item R2) is **irrelevant on the graft** — and, by contrast, load-bearing on
the peel, where an absolute density is subtracted (§9).

### 2.1 Two factors, one statement

`f_g(x)` and `f_R(x)` are emitted as two λ-factors but derive from the **same** `k`, so they sum to 1 by
construction. This retires the incoherence measured in `message_arithmetic_reconciliation.md` §E.5 defect 2
(the two emitted fractions summed to 0.753 / 0.768 / 0.767 — the fold was reconciling a self-contradiction).

---

## 3. The precision — T2

`log k = log ρ_g − log(ρ_ν + ρ_μ)`. Delta method, with `ρ_g` and `ρ_ν` **sharing** the mass `M_b`:

```
∂ log k / ∂ log M_b = 1 − w_ν = w_μ          (M_b scales ρ_g AND ρ_ν — it partially cancels)
∂ log k / ∂ log S_b = − w_μ
∂ log k / ∂ f_g     = 1/f_g + w_ν/(1−f_g)
```

giving `Var(log k) = w_μ²(1/n_b^eff + 1/n_s^eff) + [1/f_g + w_ν/(1−f_g)]²·Var(f_g)`. **Do not implement that
form.** Reparameterize onto the solver's own coordinate `λ = logit(f_g)` (independent review, 2026-07-23).
Since `w_ν = 1 − w_μ`,

```
∂ log k / ∂ f_g = [ (1−f_g) + (1−w_μ)f_g ] / [ f_g(1−f_g) ] = (1 − w_μ·f_g) / [ f_g(1−f_g) ]
```

and with `df_g/dλ = f_g(1−f_g)` the singular denominators **cancel identically**:

```
┌────────────────────────────────────────────────────────────────────────────────────┐
│  ∂ log k / ∂λ = 1 − w_μ·f_g                                                        │
│                                                                                     │
│  Var(log k) = w_μ²·( 1/n_b^eff + 1/n_s^eff )  +  ( 1 − w_μ·f_g )²·Var(λ)           │
└────────────────────────────────────────────────────────────────────────────────────┘
```

> **MC: 3.5 % / 4.0 % / 4.2 % / 1.9 % relative error** at (f_g, w_μ) = (0.40, 0.85) / (0.40, 0.00) /
> (0.40, 1.00) / (0.85, 0.93). The two forms are algebraically identical and agree to all printed digits —
> this is a pure reparameterization, so adopting it is **behaviour-preserving** and carries no A/B risk.

The coefficient `(1 − w_μ f_g)` is **bounded in `[1−f_g, 1]`**: no division by zero, no blow-up as
`f_g → 0` or `f_g → 1`. Limits: `w_μ = 0 ⇒ coefficient 1 ⇒ Var(log k) = Var(λ)` exactly; `w_μ = 1 ⇒
coefficient (1−f_g) = f_R`, which is right because `ρ_R = ρ_μ` no longer depends on `f_g`, so
`log k = log f_g + const` and `d log f_g/dλ = 1−f_g`.

### 3.0 Which variance field to use — the §6 blocker, resolved

`Var(λ)` is **already computed**: `NodeDeconv.lam_var` = `Var[λ]` over the grid posterior, on the solver's
native log-odds lattice (`simplex_logodds`). **Use it directly.**

Do **not** use `NodeDeconv.gdna_frac_var` for this. It is documented as `Var(f_g)` (`strand_deconv.py:50`)
but computes `post_lam @ (log_fg_grid²) − mLg²` — i.e. **`Var(log f_g)`** (observed max 4.049, above the
Popoviciu bound ¼ that `Var(f_g)` could not exceed). Converting it to the logit scale requires `÷(1−f_g)²`,
which reintroduces exactly the singularity this section removes, and at `f_g → 1` — the gDNA-rich regime that
matters most. The field's docstring is wrong and should be corrected separately; `enrichment_frame.
composition_logvar` explicitly expects `Var(f_g)` and is a live trap for the same reason.

### 3.1 Why this is the right combination rule

Two *independent observations of the same quantity* → **precisions add**. Two *components that sum to the
quantity* → **variances add**, share-weighted. Mature and nascent are **addends** of `ρ_R`, not two looks at
one number. The shipped `pr += S` uses the replicate rule on an additive decomposition, which asserts that
observing the departing RNA makes us more certain about the continuing RNA. It does not.

**The share-weighting is not postulated — it falls out of the delta method.** `w_μ²` and `w_ν²` appear
because the weights are squared, so a minority component contributes *quadratically* little: a large
*relative* error on a small addend is a small *absolute* error on the sum. Precision-addition is wrong in
**both** directions — the weak part does not piggyback, and the strong part is not diluted either.

### 3.2 Two limits that pin the formula

* **No spliced flux (`w_μ = 0`)** — `Var(log k)` collapses to `Var(f_g)/[f_g(1−f_g)]² = Var(logit f_g)`
  **exactly**, because `k` is then a ratio of two densities built from the *same* mass and the counting
  cancels identically. *MC: predicted 0.069444 vs empirical 0.072067.*
* **Mature-dominated (`w_μ → 1`)** — the counting no longer cancels (the two densities no longer share a
  common factor) and the composition term degrades to the `1/f_g` part alone.

### 3.3 Counts, not mass; and overdispersion

`n_b` and `n_s` are the **integer** counts (`n_unspl_*`, `spliced_n_*`), never the fractional mass — mass sums
per-fragment shares and is split across nodes by the accumulator, so `1/mass` is not a counting variance.

Under the Gamma-Poisson (NegBinom) count model the effective count is

```
n^eff = n / (1 + (n−1)·ω)        ⇒   1/n^eff → ω   as n → ∞
```

so **`ω` is a floor on a count channel's precision that no sequencing depth removes** (`ω = 0.05` ⇒ ±22 %).

**`ω` has now been measured (2026-07-24), and it is NOT fittable per-library.** Findings:

* **The synthetic suite is Poisson by construction: `ω = 0`.** `sim/wgs_engine.py:473` allocates fragments by
  multinomial `rng.choice(...,p=probs)` at fixed abundance, no Gamma/NB mixing; measured `ω < 5e-5`, flat
  across depth. **The battery therefore CANNOT validate anything that depends on `ω`.** It validates the
  `ω`-independent parts — the mode (T1), the shares, the destination Jacobian (T5), the *relative* precisions
  — but a count-precision floor fit or accepted on `ambig_dense_10mb` is over-confident by
  `√(1+(n−1)·ω_real)` at real depths (n = 100, ω = 0.02 ⇒ 1.4×).
* **Real data: `ω ≈ 0.02` (a depth-extrapolated upper bound), not point-identified.** Within-transcript MoM
  decays monotonically with depth (0.41 → 0.08 as µ climbs) — the fingerprint of a between-junction *fixed*
  effect (splicing-efficiency / eff-length bias), not a constant counting floor. Net of it, the counting `ω`
  is ≤ 0.03. Real data has no technical replicates, so it is **not fittable** (every single-library proxy is
  either structurally binomial → 0 or biology-contaminated → upward biased).

**Decision:** ship a conservative CONSTANT `ω_spliced ≈ 0.02` in `CalibrationConfig`, alongside the strand-BB
overdispersion defaults, consumed as `n^eff = n/(1+(n−1)·ω)` (the `_compile_strand_evidence` form). Flag it for
the no-magic-numbers review — it is a *chosen upper bound*, not a fitted quantity, precisely because the data
cannot identify it; that non-identifiability is itself the justification. **The graft is well-conditioned
(a sum), so a modest `ω` mis-set changes only the absolute count-floor, not the shares or the mode** — this is
why the graft can ship on the battery despite `ω` being untestable there. The channel where `ω` is decisive is
the PEEL (§9), and that is where the "conservative constant" matters.

---

## 4. From `Var(log k)` to the message precision — T5

```
d log f_g(x) / d log k = 1 − f_g(x)          ⇒  Var(log f_g(x)) = (1 − f_g(x))² · Var(log k)
d log f_R(x) / d log k = − f_g(x)            ⇒  Var(log f_R(x)) =      f_g(x)²  · Var(log k)
precision = 1 / Var
```

> **MC: 0.03 % – 4.0 %.**

**The Jacobian is evaluated at the DESTINATION's `f_g`, not the source's.** The transport itself moves `f_g`
(measured: 0.400 → 0.195, 0.614, 0.005 in the three MC regimes) because the effective lengths differ. Using
the source's value — which is what `_scan` does today (`v_logfg = fr_s²·ev_lam`) — is wrong by **45 %, 132 %,
64 %** in those same regimes. This is a real, large, previously-unnoticed defect.

Finally add the per-edge transfer term once, if it is still in use: `+ σ²_transfer`. Note that R3 (retiring
σ²_transfer) is *justified on the graft by T1* — the mode is now enrichment-exact, so the cliff-height penalty
is a pure double-count there. It is **not** justified on the peel. R3 must therefore be scoped per direction,
not applied globally (this corrects the framing in `message_arithmetic_reconciliation.md` §B.1).

---

## 5. What this retires

| retired | replaced by |
|---|---|
| `pr += S` (precision addition of the measurement) | T2's share-weighted `Var(log k)` |
| the `n_eff` "honest clamp" | falls out — `w_μ → 1` as the imputed part vanishes |
| the source Jacobian `fr_s²·ev_lam` / `fg_s²·ev_lam` | T5's **destination** Jacobian |
| the mass-vs-count inconsistency (`SPs` mass used as a count) | `n_s` integer counts in the variance, mass only in the density |
| the two incoherent λ-factors into an exon | one `k`, two factors summing to 1 (§2.1) |
| `Var(log r)` on the graft (R2) | T1 — enrichment cancels; `r` never appears |

---

## 6. Implementation

All inside `_unified_solve` (`bp_solver.py`), the graft path only.

1. **Per boundary, precompute the parts** (composition-free where possible): `ρ_μ` per face
   (`spl_*_f[face]`), `n_s` per face (`spliced_n_*`), `n_b`, `E_g^b`, `E_r^b`.
2. **At the graft edge** (`ex_a[dst] & is_bnd_a[src]`), form
   `ρ_R(src) = ρ_ν(src) + ρ_μ(src, face)` — already implemented (the graft is applied **before** the reframe
   because it is a density measured at the source) — and the shares `w_μ, w_ν`.
3. **Compute `Var(log k)`** by T2 (§3, logit form) from the source's `f_g`, **`lam_var`** (= `Var(λ)`, §3.0 —
   NOT `gdna_frac_var`), `n_b^eff`, `n_s^eff`.
4. **Compute the destination `f_g(x)`** from `k` and the exon's `E_g, E_r` (`f_g_from_k`), then apply T5's
   Jacobian to get the two factor precisions.
5. **Emit two coherent factors** on `log f_g(x)` and `log f_R(x)`; split `f_R` across strands by the tilt as
   today (the tilt is a separate axis — `k` is defined on RNA-total).

**Variance field:** resolved in §3.0 — consume `lam_var` (`Var(λ)`) directly; the logit form of T2 needs no
conversion, and the previously-flagged `gdna_frac_var` blocker no longer applies to this design. Correcting
that field's docstring (and auditing `composition_logvar`'s callers) remains worth doing separately.

---

## 7. Gates

* flag-off **byte-identical** (`gdna_none` guard unchanged);
* the **Σ_c f_c = 1** message invariant (`unified_message_audit.py`) — the cheap 1-second regression catch;
* **stranded must not regress** (R4), per condition;
* the 32-scenario battery per condition, unstranded AND stranded;
* `message_precision_mc.py` promoted to a unit test pinning T1/T2/T5;
* goldens last.

---

## 8. Honest limits

1. **Delta method.** T2/T5 are first-order. Validated to ~4 % in the regimes above; they degrade when any
   input's relative spread approaches 1 (see §9 for where that bites — it bites the *peel*, not the graft,
   because the graft is a sum and sums are well-conditioned). Measured degradation at an extreme composition:
   at `f_g = 0.05` with `σ_λ ≈ 0.6` the error reaches **17.3 %**. Nodes with a near-vertex `f_g` AND a wide
   posterior are the graft's weakest case.

2. **Independence — stated explicitly (review, 2026-07-23).** T2 assumes `Cov(λ, log M_b) = Cov(λ, log S_b) =
   0`; the omitted term would be `2·w_μ·(1 − w_μ f_g)·Cov(λ, log M_b)`. Count-zero-information is what makes
   this defensible: the count may set *precision*, never a composition *vote*, so `λ`'s LOCATION does not read
   `M_b`. It is **not exactly zero**, though — message precisions scale with `n`, so a higher-count node is
   pulled harder toward its incoming messages, which correlates `λ` with `n` at second order. Treat as a
   documented approximation, and measure the residual correlation before relying on the graft at high
   precision.

3. **Low counts break the log-Gaussian count model — and this is NOT fixable by a Taylor term.**
   `Var(log n) = 1/n^eff` is a first-order expansion. Measured against Poisson MC it **under-states**
   (empirical 0.05430 / 0.02063 / 0.01016 at λ = 20 / 50 / 100 versus `1/λ` = 0.05 / 0.02 / 0.01), so the
   graft is **over-confident** at low counts. The excess scales as ≈`1.5/n²` — *three times* the `1/(2n²)`
   second-order term one might bolt on, and the naive analytic expansion actually gives the **wrong sign**
   (`−1/(2n²)`), because the Taylor series is asymptotic while the left tail of `log` dominates at finite `n`.
   For the NegBinom (the real model) `1/n^eff` under-states by 3–28 % with a persistent ≈2.7 % floor.
   **This is material:** measured real junction counts are `n_s` p10 = 1, p25 = 4, **median 36**, p75 = 677 —
   **34.7 % of junctions have `n_s < 10`**, where `log(count)` is not remotely Gaussian and `P(count = 0)` is
   non-negligible. No expansion is valid there, and a fitted constant would violate the no-magic-numbers
   directive. Resolve by exact 1-D quadrature on the fitted NegBinom (cheap, tabulatable) or an explicit
   low-count gate — **not** by an added Taylor term.

4. **`ω` unfitted** for the spliced channel (§3.3). Until measured, the graft is *optimistic* about `ρ_μ`.
5. **Tilt.** `k` is defined on RNA-total; the ± split rides the existing tilt axis. A per-strand `k` is not
   derived here.

---

## 9. Why the PEEL is a different problem (scope boundary)

`exon → boundary` runs the other way: the boundary receives only what **continues**, so the departing flux
must be **subtracted**:

```
ρ_ν(b) = ρ_R(x)/r − ρ_μ(b)            ← r does NOT cancel: an absolute density is subtracted
Var(log ρ_ν) = u²·σ_T² + (u−1)²·σ_μ²  ,  u = T/ρ_ν = 1/(fraction that continues)
```

**A sum carries convex weights in [0,1]; a difference carries weights ≥ 1.** Same delta method, opposite
regime: subtracting similar numbers destroys precision. Consequences:

* `Var(log r)` (R2) is **load-bearing** on the peel — it carried 92 % of the peel's variance in MC — and
  **irrelevant** on the graft (T1).
* The linearization is valid only for `ε = √(u²σ_T² + (u−1)²σ_μ²) ≲ 0.15`, and through the next band it
  **under-states** the variance — the over-confident direction, and the signature of the +49 % phantom.
* `ρ_ν < 0` is **not** physically impossible; it is a sampling outcome (owner, 2026-07-23). The `ρ_ν ≥ 0`
  constraint belongs in the **prior** — the same shape as the intron factory's `NegBinom(g)·1[g ≤ C]` — not
  as an emission gate. `P(ρ_ν<0) ≈ Φ(−log(T₀/ρ_μ)/√(σ_T²+σ_μ²))` is the honest read of "can the data exclude
  *nothing continues*?".
* Since `ε ≥ (u−1)·√ω`, `ω` sets a *floor* on the peel's error — but the measurement (§3.3, 2026-07-24)
  shows **`ω ≈ 0.02` is small; it is `u` that limits the peel, not `ω`.** At `ω = 0.02` the floor clears
  `ε ≤ 0.15` only for `u ≤ 2.06`, i.e. ~40 % of junctions (the ones where >½ the RNA continues). The median
  junction (`u = 2.3`) sits at the boundary; `u ≳ 3` junctions (p75 `u = 5.3`, only 19 % continuing) are
  hopeless at any depth — `ρ_ν = T − ρ_μ` subtracts large near-equal numbers. **So the peel is `u`-weighted,
  not globally demoted:** quantitative where `u ≲ 2`, a weak one-sided constraint where `u ≳ 3`. `ε` is
  already a per-junction precision — use it as the weight directly.

### 9.1 A correction to the prior record

`prediction_measurement_merge_derivation.md` §5 claims the additive lower bound "`ρ_r = ρ_m > 0`" prevents the
+49 % phantom. That bound is valid **only on the graft**, where the destination genuinely carries mature. The
+49 % was attributed to the `exon → boundary` mature absorption — and the two mechanisms **cannot occur on the
same edge**: measured on `gdna300_ss_0.50_nrna_present_capture_on`, `pr += S` fires on **400 edges, all with a
REGION destination**, the absorption on **399 edges, all with a BOUNDARY destination**, **overlap 0** (spliced
mass is identically zero on every region node). Item E therefore does unblock R2 — but on graft edges, for the
reason in §3.1, not the reason stated there.
