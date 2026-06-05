# ρ₀ variance & count-prior widening (documented next step)

**Status:** design note for a future PR. The *clean global density* joint model is landed
(`calibrate` reorders strand→density; `density_model.gdna_frac` and `decode_sides` strand-clean
ρ₀ on contained nodes *and* boundary crossings; the decode is joint count×strand; `confidence`
is a posterior quantile). This note records the next refinement: propagating the **uncertainty**
of the strand-cleaned global density ρ₀ into the count clue.

## 1. What's landed, and why it's principled

The global gDNA density **ρ₀ is a hyperparameter**, estimated by empirical Bayes from the
count-decodable nodes (intergenic = pure gDNA; introns = strand-cleaned via the closed-form
gDNA fraction `π_g = clip((sense_frac − κ_rna)/(½ − κ_rna), 0, 1)`). Each node's per-node decode
then **conditions** on ρ₀ and κ_rna and combines the two clues:

```
posterior(π_g | node) ∝ Beta(π_g ; density·eff_len/M, κ_c) · BB_strand(sense, antisense | π_g)
```

This is **not** vicious circularity:

- ρ₀ is global; a single node contributes ~`1/N` to it, then uses its own strand fully in its
  own decode. The `1/N` self-use is negligible (leave-one-out makes it exact if ever needed).
- The main decode targets (exonic nodes) never enter ρ₀ (only intergenic + introns do), so for
  them count and strand are **fully independent** channels.
- Conditional independence given `(ρ₀, κ_rna)` makes the product legitimate.

It also degrades correctly across the **strand-specificity spectrum**: unstranded (κ_rna = ½)
⇒ the strand BB is degenerate ⇒ count-only, and ρ₀ falls back to the raw count density
(strand can't clean it); increasingly stranded ⇒ the strand cleans ρ₀ better *and* contributes
more to each node — the synergy.

## 2. The gap: ρ₀ enters as a point estimate

Today the count clue uses ρ₀ as a **point** value. Its uncertainty does not propagate. But the
strand-cleaning is imperfect and ρ₀ is finite-sample:

- **Strand-cleaning uncertainty.** `π_g` per decodable intron is a *finite-count, finite-SS*
  estimate. As κ_rna → ½ (low strand specificity) the map `(sense_frac − κ)/(½ − κ)` has a
  large derivative, so small sense-fraction noise becomes large gDNA-fraction noise — the
  cleaned gDNA mass feeding ρ₀ is noisy.
- **Count uncertainty.** ρ₀ aggregates a finite number of gDNA fragments; sparse libraries (or
  capture, where off-target intergenic is depleted) give few count-decodable fragments.

The point-ρ₀ count clue is therefore **over-confident** exactly when ρ₀ is poorly determined
(low SS, sparse gDNA). Observed symptom: at SS 0.9 + 20 % gDNA the silent-gene negative control
leaks ~5 % (the count no longer over-masks the strand noise that the *old* contaminated ρ₀ hid).
This is the diagnostically-relevant regime — imperfect strand specificity with real gDNA — so it
is worth modeling rather than papering over with the `confidence` knob (which is the output
quantile and trades the leak for gDNA over-calls elsewhere).

## 3. The model: ρ₀ variance → count-prior widening

Carry ρ₀ as a **distribution**, not a point, and let its spread widen the count prior.

### 3.1 ρ₀ as a Gamma posterior

Treat the strand-cleaned gDNA as a Gamma–Poisson hyperposterior over the global rate:

```
G_clean = Σ_decodable  π_g,r · M_r           (strand-cleaned gDNA mass)
L_eff   = Σ_decodable  eff_len_r
ρ₀ ~ Gamma(a₀, b₀),   a₀ = G_clean_effective + ½,   b₀ = L_eff
E[ρ₀] = a₀/b₀,        Var[ρ₀] = a₀/b₀²
```

`G_clean_effective` is the **effective** cleaned count, discounted for the strand-cleaning
uncertainty: a per-node `π_g` estimated from `n_r` oriented fragments at specificity κ_rna has a
Beta-Binomial effective sample size `n_eff,r = n_r / (1 + ρ_strand·(n_r − 1))`, and its
contribution to ρ₀ is down-weighted by `(½ − κ_rna)²` (the squared cleaning gain — vanishing as
κ_rna → ½). So low SS or thin counts ⇒ small `a₀` ⇒ large `Var[ρ₀]`.

### 3.2 Widening the count prior

The count clue's gDNA-fraction mean is `π_count = ρ₀·eff_len / M`. With ρ₀ uncertain, integrate
it out (or moment-match) so the **concentration** `κ_c` of the count Beta shrinks as `Var[ρ₀]`
grows:

```
π_count   = E[ρ₀]·eff_len / M
κ_c_eff   = κ_c · E[ρ₀]² / (E[ρ₀]² + Var[ρ₀])          (≤ κ_c; → κ_c when ρ₀ is well-determined)
a_c = κ_c_eff·π_count + ½,   b_c = κ_c_eff·(1 − π_count) + ½
```

When ρ₀ is well-determined (`Var[ρ₀] ≪ E[ρ₀]²`) this reduces to today's behaviour. When ρ₀ is
uncertain (low SS, sparse gDNA), `κ_c_eff → 0` ⇒ the count prior flattens to Jeffreys ⇒ the
decode leans on whatever strand evidence the node has, instead of trusting an ill-determined
reference density. This is the correct propagation: the silent-gene call stops over-trusting a
shaky ρ₀, and the nascent/RNA calls stop being whipped by it either.

## 4. Why this is the right axis (and the `confidence` knob is not)

`confidence` is the **output** quantile (how conservatively we *report* π_g — the FP-averse
dial the user sets per sample). ρ₀ variance is an **input** uncertainty (how well we *know* the
reference density). They are orthogonal: raising `confidence` to mask an uncertain ρ₀ over-calls
gDNA everywhere (validated: default 0.8 sent antisense 37→169 and added RNA-conservation
failures). The variance model fixes the *input*, leaving `confidence` free for its real job.

## 5. Scope / acceptance for the future PR

- Add `rho0` and `rho0_var` to the calibration internals (and optionally the result), computed
  in `density_model` from the cleaned `a₀, b₀`.
- Thread `κ_c_eff` into `joint_decode._joint_per_node`.
- Acceptance: SS-0.9 + 20 % gDNA negative control returns under the 25-count limit *without*
  raising the `confidence` default; nascent stays fixed; capture recovery and the 10:1 FP
  protection are unchanged; multimap stays green.
- Non-goal: per-node leave-one-out ρ₀ (the `1/N` self-use is negligible at scale).
