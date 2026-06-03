# ρ₀ Stabilization — Mechanism and Design

**Companion to** `exposure_prior_findings_report.md` and `exposure_prior_noise_floor_theory.md`.
**Premise (established):** the sparse-data failure is the `(ρ₀, ω)` identifiability
bifurcation, not the exposure prior. A prior on the dispersion `φ` is the wrong lever; the
gDNA density `ρ₀` collapses, and *that* is what mis-routes mass. This note pins the
collapse mechanism in the code, presents the decisive experiment, and proposes a
principled, data-derived stabilization.

---

## 1. Mechanism — why ρ₀ collapses

Two estimators of `ρ₀` coexist in the calibration:

**The seed** (`density.py::estimate_gdna_density`) — *robust*. It pools the **non-exonic**
regions (intergenic + intronic, gDNA-dominated *by annotation*), treats their unspliced
mass as gDNA, and forms a Gamma posterior with a **unit prior**:
`ρ₀ = (1 + Σ_seed M_g) / (1 + Σ_seed L_eff)`. On the `nrna_dc` toy this gives ρ₀ ≈ 0.055
(true ≈ 0.04; mildly upper-biased by intronic nascent RNA, as its own docstring notes).

**The EM update** (`mstep.py::update_rho_0`) — *fragile*. Every outer iteration it
**replaces** ρ₀ with the bare allocation ratio over **all decodable** regions, with **no
prior**:
`ρ₀ = Σ_decodable m_g / Σ_decodable(ω·L_eff)`.

The seed is used only as the initial value (`calibrate.py:105`); the EM update then throws
away all three of the seed's robustness features (non-exonic pool, mass-as-gDNA, unit
prior). The collapse follows:

1. The count channel (`estep._llr_count`) compares gDNA-only mean `μ_g = ω·ρ₀·L` against a
   mixture `μ_g + m_d_prev`. Once a region carries any RNA mass `m_d_prev`, the mixture
   fits the observed flux better ⇒ `π_g ↓` ⇒ allocated gDNA `m_g ↓`.
2. `update_rho_0` then sees `Σ m_g ↓` in the numerator, while the **exonic** regions keep
   their `ω·L_eff` in the denominator (mass there is RNA, so they contribute length but no
   gDNA). ρ₀ is diluted downward.
3. Lower ρ₀ ⇒ lower `μ_g` ⇒ the count channel prefers the mixture even more ⇒ `m_g → 0`
   everywhere ⇒ `ρ₀ → 0`. A circular collapse.

Trace (`nrna_dc`, n=2000, with the φ prior forcing ω→1): ρ₀ = **5e-5** (~800× too small),
per-region called gDNA ≈ 0 everywhere, gDNA observed 205 (expected 596), the freed mass
leaking to nRNA (1175 vs 828), RNA over-recovery 128%.

This is the second face of the bifurcation: with ω free, a few small-R regions instead
*over*-attract gDNA (gDNA over-called); with ω→1, ρ₀ collapses (gDNA under-called). The
exposure prior cannot escape it from either side.

---

## 2. The decisive experiment

Freeze ρ₀ at the seed (skip the EM update entirely); rerun `nrna_dc`, n=2000:

```
ρ₀ = 0.055 (seed)            gDNA:  expected 596, observed 596.0   ← exact
t_ctrl regions → gfrac 1.00  total RNA: expected 1404, observed 1404.0   ← 100%
alpha_gdna_add 441 / 413     nRNA: 792 (expected 828);  t_ctrl leak 27 (minor)
(was 0.4 under collapse)
```

The headline goes **20% → 100%**, gDNA **exact**, and — critically — `exp_disp` floors
(ω→1) yet gDNA is now correctly called: **with ρ₀ healthy, ω→1 is safe.** The bifurcation
is resolved by anchoring ρ₀ alone.

Across the full harness with ρ₀ anchored (κ=100):

- **Sparse sweep:** 90 / 97 / **100** / 104 / 99 / 66 % (was a chaotic 17–108%).
- **Dense-uniform:** K=1 100%, K=5 86%, K=20 87%.
- **Dense-capture:** captured ω mean **9.6** (vs ≈4.3 *without* the anchor) — enrichment
  preserved and sharpened; uncaptured ω ≈ 0.08.

And the κ brittleness is gone: the catastrophic over-recovery that made κ a knife-edge was
the collapse; once ρ₀ is anchored, larger κ is simply "more uniform," monotone and safe.
The uniform-vs-capture φ magnitudes now separate cleanly (≈0.1–0.5 vs ≈10), so the φ prior
becomes an easy secondary regularizer rather than the load-bearing fix.

---

## 3. Design — anchor ρ₀ to the gDNA-representative regions

The seed is the right estimator; the EM update must not be allowed to collapse below it.
The principle: **estimate ρ₀ only where gDNA is observable — the unexpressed / non-exonic
regions — and anchor it with the seed's Gamma prior; never let the expressed regions
(where gDNA is invisible under RNA) drag it down.** This is exactly the external review's
Step 2 ((1−γ_r)-weighting, with γ_r the expression posterior) and the "expressed/unexpressed
latent" already flagged as future work in `density.py`'s docstring.

Concretely, two staged options:

### Stage 1 — seed-anchored, signature-gated update (immediate, validated in principle)

Replace the bare ratio with a seed-prior-anchored update restricted to the non-exonic
decodable pool:

```
ρ₀ = (α_seed + Σ_{nonexonic ∩ decodable} m_g) / (β_seed + Σ_{nonexonic ∩ decodable} ω·L_eff)
```

- `(α_seed, β_seed)` is the seed's own Gamma posterior — a strong, data-derived anchor (not
  a constant), so the allocation can refine ρ₀ but cannot collapse it.
- Restricting to non-exonic removes the exonic-length dilution of the denominator (step 2
  of the collapse).
- Uses the existing `_non_exonic` mask from `density.py` — no new machinery, no new
  constant (the unit Gamma `α₀=β₀=1` already in the seed is the only prior).

This is the freeze experiment made principled: it tracks the seed when the allocation is
uninformative and refines it when the allocation carries signal.

### Stage 2 — expression latent (refinement, for high-expression / capture libraries)

Stage 1 inherits the seed's mild upper-bias: a transcribed gene's introns hold nascent RNA
that the non-exonic pool counts as gDNA. The principled remedy is a **per-region expression
posterior γ_r derived from non-circular RNA evidence** — spliced mass (gDNA cannot splice)
and sense-strand excess (gDNA is unstranded) — neither of which depends on ρ₀ or ω. Then
`(1−γ_r)`-weight the ρ₀ (and φ) pools so transcribed introns are softly excluded:

```
ρ₀ = (α_seed + Σ_r (1−γ_r) m_g) / (β_seed + Σ_r (1−γ_r) ω·L_eff)
```

This subsumes Stage 1 (the signature gate is a hard `γ_r ∈ {0,1}` prior; the latent makes
it soft and data-driven) and handles capture, where intergenic gDNA is depleted and the
representative gDNA lives near transcribed regions.

---

## 4. Composition with the φ prior, and remaining items

- **The φ prior stays, demoted.** With ρ₀ anchored, the φ regularizer's only job is to pick
  ω→1 for genuinely-uniform libraries vs ω-enriched for capture. The magnitudes are now
  well-separated, so this is robust — but it should be revisited for parameter-freeness in
  light of the `noise_floor_theory` §7 conclusions (the empirical null refuted a closed-form
  floor; with ρ₀ anchored the φ estimate is cleaner and the question may simplify).
- **Capture gDNA quantification** is still under-recovered (≈689/2740) even though the ω map
  is correct — a separate downstream item (the captured-region gDNA count, not the
  exposure shape).
- **Stage-1 risk:** restricting ρ₀ to non-exonic regions assumes a usable non-exonic
  annotation; the seed already falls back to a half-and-half density when none exists.

---

## 5. Recommendation

Implement **Stage 1** first (seed-anchored, non-exonic-gated `update_rho_0`) and validate
against the three-regime harness and the scenario/golden suites — the freeze experiment
shows it recovers the sparse failure with gDNA exact while preserving capture. Treat
**Stage 2** (the expression latent) as the principled follow-up that removes the seed's
intronic upper-bias for high-expression and capture libraries. The φ prior remains as a
secondary regularizer, to be made parameter-free separately. No floor constants are
introduced at any stage; every anchor is a data-derived Gamma prior.
