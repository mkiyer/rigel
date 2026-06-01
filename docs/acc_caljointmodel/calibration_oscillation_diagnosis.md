# Calibration φ-oscillation — root-cause diagnosis (PR 6 follow-up)

**Symptom.** On low-gDNA scenarios (`ten_genes`, the pipeline smoke fixture, the
`g*_n0_*` scenario sweeps) the calibration outer loop never converges: it settles
into a **period-2 limit cycle** where the NB dispersion `φ` ping-pongs between its
floor and ~3–4 every iteration, while `ρ₀` slowly inflates. High-gDNA scenarios
(`gdna_heavy`) converge cleanly in ~7 iterations.

Reproducer: `scripts/debug/trace_calibration_oscillation.py` (builds the smoke
oracle scenario — 3 transcripts, 5000 bp, **no gDNA, no nRNA, 0 spliced reads**,
67 intra-exon unspliced reads across 5 exon regions; 11 regions total).

## 1. The mechanism (traced fragment-level)

The φ M-step (`mstep.fit_phi`) and the exposure posterior (`exposure_posterior`)
form a positive-feedback loop **mediated by the large empty (n_u=0)
intergenic/intron regions**.

The exposure posterior gives every region
`ω = (1/φ + M_g) / (1/φ + ρ₀·L_eff)`. For an **empty** region (M_g = 0):

```
ω_empty = 1 / (1 + φ·ρ₀·L_eff)        ⇒   μ_g = ω·ρ₀·L_eff = ρ₀L_eff / (1 + φ·ρ₀·L_eff)
```

So the *modeled gDNA mean of an empty region is a function of φ* — and `fit_phi`
is then fed exactly those `μ_g` and asked to fit φ from `NB(n_u | μ_g + m_d, φ)`.
Circular, with the wrong sign:

| | empty region ω | empty region μ_g | n_u | `fit_phi` sees | → returns |
|---|---|---|---|---|---|
| E-step φ **small** (0.037) | 0.17–0.53 | **~26** | 0 | `NB(0\|26)` — catastrophic misfit | φ **large** (~3.2) |
| E-step φ **large** (3.2) | 0.002–0.009 | ~0.27 | 0 | `NB(0\|0.27)` — fine | φ **small** (floor) |

Per-region NB-NLL decomposition at the two cycle states (smoke, iters 4↔5),
`Δ = NLL@φ=4 − NLL@φ=floor` (Δ<0 ⇒ region pulls φ up):

```
iter 4 (E-step φ=0.037):  region 6 (intergenic,L=1203,n_u=0) Δ=-25.1
                          region 10(intergenic,L=1603,n_u=0) Δ=-27.3
                          region 8 (intron,    L=103, n_u=0) Δ= -5.3
                          4 exon regions (n_u>0)             Δ= +7.8 total
                          ───────────────────────────── TOTAL Δ = -50  ⇒ fit_phi → φ=3.2
iter 5 (E-step φ=3.2):    same empty regions now μ_g≈0.27    Δ≈ -0.2 total
                          4 exon regions                     Δ= +6.6 total
                          ───────────────────────────── TOTAL Δ = +6.8 ⇒ fit_phi → φ=floor
```

The two states are mutually reinforcing → a stable 2-cycle that **never**
converges. (`ρ₀` also ratchets up 0.08→0.16 across the cycle — a pure-mRNA
scenario where the true ρ₀≈0.)

## 2. Root cause — `fit_phi` is the wrong M-step

`φ` is **one** parameter playing two roles that the Poisson–Gamma identity ties
together (doc 01 §4.2): it is the NB count dispersion **and** the exposure-prior
variance (`ω ~ Gamma(1/φ, 1/φ)`). The exposure posterior uses it correctly. But
the **M-step for φ is mis-specified**:

- `fit_phi` maximizes the *marginal* NB count likelihood
  `Σ_r log NB(n_u_r | μ_g_r + μ_d_r, φ)` — **with the posterior mean `ω_r` plugged
  into `μ_g_r = ω_r·ρ₀·L_eff`.** That double-uses the exposure (ω already absorbed
  the count; φ then re-fits the residual that ω's value depends on). It is neither
  the clean marginal (which would use `ρ₀·L_eff`, no ω) nor the EM M-step, so it
  carries no monotone-ascent guarantee → it oscillates.
- The **proper EM M-step** for φ maximizes the expected complete-data
  log-likelihood, whose only φ-term is the exposure prior
  `Σ_r E_{ω|data}[ log Gamma(ω_r | s, s) ]`, `s = 1/φ`. That is the standard
  Gamma-shape MLE from the posteriors `Gamma(α_post_r, β_post_r)`:
  ```
  solve   log s − ψ(s) = mean_r( E[ω_r] − E[log ω_r] ) − 1     for s,   φ = 1/s
  E[ω_r]     = α_post/β_post                 (= ω, already computed)
  E[log ω_r] = ψ(α_post) − log(β_post)
  ```
  It uses **only the ω posteriors**, never the raw counts → no empty-region
  count-misfit feedback, and EM guarantees monotone convergence.

This is a correctable algorithm bug, **not** a case for a damping factor.

## 3. Empirical comparison (smoke scenario, 14 iterations)

| φ M-step | outcome | φ | ρ₀ |
|---|---|---|---|
| `nb_omega` (current) | **period-2 cycle, never converges**; ρ₀ inflates | floor ↔ 3.1 | 0.08→0.16 |
| `nb_active` (drop n_u=0 regions) | converges | floor | 0.079 |
| `nb_marginal` (mean = ρ₀L_eff, no ω) | slow drift, no fixed point | ~4.3↑ | ~0.08↑ |
| **`gamma` (proper EM M-step)** | **converges by iter ~5** (Δ~1e-9) | ~120 | **0.0035** |

`gamma` both converges and yields the most-correct ρ₀ (≈0, matching the
no-gDNA truth). Its large φ is *adaptive and arguably correct*: when gDNA
exposure is extremely heterogeneous (gDNA concentrated in a few regions, ~0
elsewhere), a high exposure variance is the right inference — and it makes the
count channel appropriately weak, deferring to the strand channel + prior.

## 4. Recommended fix

Replace `mstep.fit_phi` (the NB-count NLL minimiser) with the **Gamma exposure-
variance EM M-step** (§2). Inputs: the per-region `α_post = 1/φ_old + M_g` and
`β_post = 1/φ_old + ρ₀·L_eff` already formed by `exposure_posterior`. No new
parameter, no damping. `calibrate()` passes the exposure posteriors (or `M_g_tot`
+ `ρ₀` + `L_eff` + `φ_old`) to the new M-step instead of `(n_u, μ_g, m_d)`.

**Open items to validate when implementing:**
1. Confirm convergence on `ten_genes` + the failing scenario sweeps, and that it
   does **not** regress the currently-converging gDNA scenarios (`gdna_heavy`,
   `gdna_light`) — re-run the full scenario suite + check accuracy.
2. φ can land large (≈120) under extreme heterogeneity; revisit the `_PHI_MAX=100`
   ceiling (raise it, or confirm clipping to 100 is harmless for the count BF).
3. A small residual ρ₀↔φ drift remains after the mass-change converges (benign,
   Δ~1e-9); confirm it's acceptable or tighten.

## 5. Secondary degeneracies (related, not the oscillation)

The smoke scenario has **0 spliced reads**, which makes two other parameters
degenerate (independent of the φ cycle, but worth noting for "tiny scenario"
realism):
- `ε_s → 1` (clipped): `update_eps_s` is unidentifiable with `Σ n_s = 0`.
- `κ_rna, ρ_r_bb → 1e-6` (floor): `fit_strand_balance` reads the spliced channel,
  which is empty here.

`ten_genes` oscillates **with** spliced reads, confirming the φ cycle is the
count-channel/exposure feedback above, not the strand degeneracy.

## 6. Implementation + validation (the fix landed)

`fit_phi` was replaced by `update_exposure_dispersion` (the Gamma-variance EM
M-step, §2). The `phi` family was also renamed to intuitive names
(`exposure_dispersion`, `exposure_dispersion_floor`, `_EXPOSURE_DISPERSION_MAX`).

**The oscillation is fixed.** The smoke fixture and `ten_genes` now converge
(`ten_genes`: iter 13, δ=7.8e-5). Full suite improved **962→974 pass, 23→12
fail**; calibration unit tests 112/112; goldens regenerated from the converged
output.

**But validation surfaced a deeper, separate problem — NOT the oscillation:**
1. On low-gDNA / low-spliced scenarios `exposure_dispersion` pins at the
   ceiling (100). A huge dispersion makes the NB count channel (which uses
   `1/dispersion`) essentially uninformative — so the gDNA-vs-RNA call falls to
   the strand channel + prior.
2. When the strand model is *also* degenerate (no spliced reads ⇒ κ_rna→floor),
   nothing correctly identifies unspliced RNA, and it is **over-called as gDNA**.
   `single_exon`: t1 expected 339, observed 41 (87.9% error) — 298 of t1's reads
   assigned to gDNA where the truth is 0.
3. Many scenario-sweep runs still hit the 25-iter cap with δ spanning ~1e-4
   (basically converged — a few more iters / the slow exposure_dispersion drift)
   up to ~0.25 (genuinely unconverged).

These are calibration-*accuracy* issues (count-channel disabling ↔ exposure
heterogeneity, degenerate strand model on no-spliced scenarios, slow
convergence), distinct from the period-2 cycle this fix removed. They are the
next investigation (PR 7 validation / a follow-up deep-dive). Candidate threads:
whether the count channel should be down-weighted by `1/dispersion` at all when
dispersion saturates; the strand model's behaviour with no spliced reads; and
whether `max_outer_iterations` / `_EXPOSURE_DISPERSION_MAX` need raising.
