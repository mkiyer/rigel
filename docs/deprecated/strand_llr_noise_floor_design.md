# Calibration v4 — Strand LLR Noise Floor Design

**Date:** 2026-04-19
**Status:** design notes — layer (1) implemented; layer (2) deferred.
**Companion doc:** `strand_llr_perfect_ss_failure.md` (root cause).

---

## 1. Context recap

The v4 calibration mixture on the mini-stress sweep is accurate to within
~9 % at `ss ≤ 0.9` but overestimates `λ_G` by 1.45 – 2.68× at `ss = 1.0`.
Root cause: the `1 − ε` clamp in `_strand_llr` with `_EPS = 1e-12` turns a
single antisense read into +26.94 nats of gDNA evidence, flipping γ→1 in
any region with fewer than ~40 sense reads and cascading into π_soft≈0.99,
σ_R inflation, and a biased μ_R. See the companion doc for numeric detail.

This doc captures the **model reasoning** behind the fix we've shipped and
the self-consistent extension we've deliberately deferred, so a future
reader (or future us) can pick up exactly where we left off.

---

## 2. What the two mixture classes really are

The calibration EM partitions **regions**, not reads, into two classes.
This matters for the strand LLR because it changes what antisense rate we
should expect under each class.

- **Class R — RNA-bearing region.** Reads are a mixture of RNA and the
  gDNA background that blankets the whole genome.  Under Class R:

  $$
  P(\text{anti} \mid R, i) \;=\; (1 - f_{G,i})\,(1 - ss_{\text{true}}) \;+\; \epsilon_{\text{bio}} \;+\; f_{G,i}\cdot 0.5
  $$

  where $f_{G,i}$ is the local gDNA contamination fraction in region $i$
  and $\epsilon_{\text{bio}}$ captures structural RNA antisense (readthrough,
  convergent transcription, antisense ncRNA, template-switching artefacts).

- **Class G — gDNA-only region.** Reads are pure gDNA background with
  $P(\text{anti} \mid G) = 0.5$ by construction.

The previous implementation used $P(\text{anti} \mid R) = 1 - ss_{\text{est}}$,
which — at `ss_est ≥ 1 − 10⁻⁶` — is **smaller than the gDNA background
contribution its own model predicts should be present in every RNA region**.
That's the mathematical self-contradiction that makes the LLR explode.

---

## 3. Three sources of irreducible antisense rate

Ordered by typical magnitude in real data:

| source                                   | typical magnitude           | scales with N_train? |
|------------------------------------------|-----------------------------|:---------------------:|
| in-region gDNA background  `0.5 · f_G,i` | 0.01 – 0.10                 | no (depends on λ_G, L_i) |
| biological/technical RNA antisense `ε_bio` | 0.001 – 0.010             | no (structural)      |
| measurement uncertainty on ss_est `ε_CI` | √(ss(1−ss)/N_train)         | **yes** — 1/√N       |

Only the third shrinks with library depth.  The other two form a floor
that **does not go away at 10⁸ reads**.  A Beta-Binomial likelihood with
posterior parameters derived from $N_{\text{train}}$ converges to the
plug-in Binomial as $N_{\text{train}} \to \infty$ and therefore cannot
protect against the structural floors — its effective sample size washes
out under real library depths.  This is the scaling trap the user flagged.

The right decomposition is: use the measurement uncertainty `ε_CI` when
the strand model is under-trained, and a small biological floor `ε_bio`
when it is well-trained, and ignore the in-region background for now
(it's a per-iteration self-consistent refinement — see §6).

---

## 4. Shipped fix — layer (1)

```python
# src/rigel/calibration/_em.py::_strand_llr  (new)
eps_floor = max(epsilon_ci, STRAND_SPECIFICITY_NOISE_FLOOR)
ss_c = float(np.clip(strand_specificity, 0.5, 1.0 - eps_floor))
log_ratio_sense = math.log(0.5 / ss_c)
log_ratio_anti  = math.log(0.5 / (1.0 - ss_c))
```

- `STRAND_SPECIFICITY_NOISE_FLOOR = 0.001` — exposed on `CalibrationConfig`,
  interpreted as "structural antisense rate under RNA".  Default 0.001
  rewards pristine stranded kits (logged `ss_est > 0.999` is common in
  user's corpus of >100 k samples).
- `epsilon_ci = 1 − UCL_{99}(ss)` — one-sided 99 % upper credible limit on
  `ss` under a Beta(1, 1) prior given the strand trainer's
  `(n_sense, n_anti)` spliced-read counts.  Computed by the strand trainer
  itself and carried on the `StrandModel` object.
- The effective floor is `max(ε_CI, ε_bio)`:
  - small N_train: `ε_CI` dominates → LLR automatically bounded by
    measurement uncertainty.
  - large N_train: `ε_bio = 0.001` dominates → LLR caps at
    `log(0.5 / 0.001) ≈ 6.2` nats.  A region needs ~9 sense reads to
    offset 1 antisense read.  Strong, informative, non-catastrophic.

Numeric reference, per antisense read, `log_ratio_anti`:

| ε        | nats   | sense reads to offset |
|---------:|-------:|----------------------:|
|  0.01    |  3.91  |  ~5.6                 |
|  0.005   |  4.60  |  ~6.6                 |
|  0.001   |  6.21  |  ~9.0                 |
|  0.0001  |  8.52  | ~12.3                 |
|  1e-12   | 26.94  | ~38.9 ← old behaviour |

---

## 5. Logging

The calibration module now emits (at INFO) per-region-EM-run:

- `ss_est` from strand model.
- `n_sense`, `n_anti`, `ε_CI` from strand model (what the trainer saw).
- `STRAND_SPECIFICITY_NOISE_FLOOR` (current config value).
- effective `ε_floor = max(ε_CI, ε_bio)` and which term dominates.
- resulting `ss_c`, `log_ratio_sense`, `log_ratio_anti`.
- final `π`, `π_soft`, `μ_R`, `σ_R`, `λ_G`.

This exposes every knob that touches the LLR and makes the floor logic
auditable from a single `grep` of the run log.

---

## 6. Deferred — layer (2): self-consistent in-region background

The proper Class-R antisense probability is

$$
P(\text{anti} \mid R, i) \;=\; (1 - f_{G,i})\,(1 - ss_c) \;+\; f_{G,i}\cdot 0.5
$$

with $f_{G,i} = \min\!\big(1,\; \lambda_G \cdot L_i / n_i\big)$ taken from the
previous EM iteration (exactly as we already do for the count LLR — the
strand LLR would join the same self-consistent outer loop).  This removes
the remaining hand-tuned component by letting the model explain its own
observed antisense reads as expected gDNA leakage.

We have deliberately **not shipped** this, because:

1. Layer (1) alone should eliminate the detonation across the sweep.
2. It adds a coupling between strand LLR and $\lambda_G$ that we'd want
   a separate convergence test for.
3. If layer (1) validation is clean, layer (2) is information-theoretic
   nice-to-have, not a correctness fix.

We should revisit layer (2) **only if** the post-(1) sweep still shows
bias at ss = 1.0 with high λ_G, or if vcap oracle σ_R doesn't drop.

---

## 7. Validation gates (same as companion doc)

Pass criteria for layer (1):

- `|λ_est/λ_true − 1| ≤ 0.15` across every ss row of the 4×5 mini-stress
  sweep, **including** ss = 1.0.
- `σ_R ≤ 0.8` in every cell.
- γ distribution unimodal in all `gdna_fraction = 0` cells.

Then vcap oracle re-run: expect σ_R 1.85 → ≤ 1.0 and −11 % gDNA bias to
close to within a few percent.

---

## 8. Open questions for the future

- Does the per-iteration `f_{G,i}` coupling (layer 2) destabilise the EM
  fixed point, or does it just tighten the estimator?
- Should `STRAND_SPECIFICITY_NOISE_FLOOR` be per-chemistry (dUTP vs
  SMARTer vs NEB Directional)?  Worth auto-detecting from the
  strand-model empirical antisense rate on high-coverage regions.
- Can we measure `ε_bio` empirically as the asymptotic antisense fraction
  in the top-coverage fully-mappable exonic bins?  That would turn the
  last remaining constant into a data-driven estimate.
