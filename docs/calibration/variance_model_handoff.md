# Variance model for the composition (unified) solver — derivation handoff

**Status:** WIP, **expected to be redone**. This is a HANDOFF that records what was derived, what was
validated, and — critically — the parameterization error found at the end that means the clean derivation is
still ahead. Read `ROADMAP.md` §3 first. Date: 2026-07-24.

**Validation harness (keep, reuse):** `scripts/debug/message_precision_mc.py` — seeded Monte-Carlo that pins
the variance laws against empirical spread. Any new derivation should be validated the same way before code.

---

## 0. Why a new variance model is needed (the conceptual root)

The old solver transported **absolute densities** between nodes and its variance model assumed **genome-wide
density uniformity** (`ρ_g ≡ ρ_bg` everywhere). Hybrid capture **breaks uniformity**: it enriches on-target
regions and depletes off-target ones by factors of 10²–10³. The new (unified) solver handles the *mode*
correctly by switching to a **composition** transport — it reframes each message by an **enrichment ratio**
`r = ρ_tot(dst)/ρ_tot(src)` that cancels capture efficiency, then normalizes by the destination’s own observed
mass (**÷M_dst**). But the **variance model was never re-derived for this composition transport.** That is the
gap.

## 1. What the message actually is (ground truth for any derivation)

**RNA is RNA** (owner). There is one RNA species; “mature/nascent” are two *routes* at a junction: it either
**splices out** (measured in the spliced channel, `ρ_μ = S_b/E_spl`) or **continues** unspliced (imputed,
`ρ_ν`). Pass-0 solves the unspliced pool. A boundary→exon message therefore carries a gDNA density and an RNA
density that is itself a **sum** `ρ_R = ρ_ν + ρ_μ` (imputed continue + measured splice).

Two message directions, NOT mirror images:
* **GRAFT** (boundary→exon): the exon receives the **whole** flux → `ρ_R = ρ_ν + ρ_μ` is a **SUM**;
  component sets match; enrichment **cancels**.
* **PEEL** (exon→boundary): the boundary receives only what **continues** → `ρ_ν = ρ_R(x)/r − ρ_μ` is a
  **DIFFERENCE**; enrichment does **not** cancel (an absolute density is subtracted).

## 2. What was derived and MC-validated (correct, but in the WRONG parameterization)

All in the **ratio (k) parameterization**, `k = ρ_g/ρ_R`, validated in `message_precision_mc.py`:

* **T1 — enrichment cancels on the matched-set graft.** The reframe + re-form-f_g is exact to machine
  precision over e ∈ [10⁻³, 10⁶]. ✅ (This part is parameterization-independent and stands.)
* **T2 — graft `Var(log k)`** (share-weighted, logit form):
  `Var(log k) = w_μ²(1/n_b + 1/n_s) + (1 − w_μ·f_g)²·Var(λ)`, `w_μ = ρ_μ/ρ_R`. MC ~4%. The **item-E
  share-weighting falls out of the delta method** (a SUM → variances add, share-weighted, squared weights),
  which is the correct rule vs the shipped `pr += S` (a replicate rule on an additive decomposition).
* **T5 — destination Jacobian:** `Var(log f_g) = (1 − f_g^DST)²·Var(log k)`, evaluated at the **destination’s**
  fraction (the old code used the source’s — wrong by 45–132%). MC ~1–4%.
* **PEEL (T3/T4):** `ρ_ν = ρ_R(x)/r − ρ_μ`; `Var(log ρ_ν) = u²·σ_T² + (u−1)²·σ_μ²`, `u = 1/(fraction that
  continues)`. A DIFFERENCE carries weights ≥ 1 (subtracting near-equal numbers destroys precision) — the
  mirror of the graft’s convex weights. Valid only for `ε = √(that) ≲ 0.15`; beyond it the linearization
  **under-states** the variance (over-confident) — the signature of the +49% phantom that killed R2.

## 3. ⛔ THE ERROR THAT SENDS IT BACK — ratio vs per-component (owner, 2026-07-24)

T2/T5 are written on `k = ρ_g/ρ_R` (a **ratio**), which is the correct object for the **shift mode**
(normalize by the imputed total). But **the unified solver uses the ÷M_dst density mode**, not the shift.
Consequences the ratio form gets wrong:

* **It is SINGULAR at the pure-gDNA anchor** (`f_g = 1 ⇒ k = ∞`). The intergenic anchor is the single most
  important source on unstranded data (the only decided node), and the ratio form breaks exactly there.
* **The owner’s decomposition is per-component, not ratio:** an intergenic anchor has **composition variance
  0** (structurally pure gDNA) but still **count variance (Poisson `1/n_b`)** and **transfer variance**. The
  ratio form makes the count cancel (`w_μ = 0 ⇒ 0`), contradicting this.

**The right framework is the per-component density variance** (this is the ÷M mode’s natural variance, and it
matches the code’s data flow — the relay transports per-component densities `cg, cp, cn`):

```
    Var(log ρ_c) = Σ_k (ρ_k/ρ_c)²·v_k          +  σ²_transfer          (once per edge)
       gDNA (imputed)         : v = 1/n_unspl + Var(log f_g^src)
       RNA-continue (imputed) : v = 1/n_unspl + Var(log f_R^src)
       RNA-splice (measured)  : v = 1/n_spliced      (no composition term — spliced IS pure RNA)
```

This **unifies everything the scattered docs were circling**:
* **Intergenic gDNA** = `1/n_b + 0` — count precision, zero composition variance. ✅ (owner point 3)
* **Graft RNA** = share-weighted `w_ν²(1/n + comp) + w_μ²(1/n_spl)` — item E, recovered, non-singular. ✅
* **Transfer** = the clean additive `+ σ²_transfer` hook. ✅ (owner point 5)

The ratio form (T2/T5) and the per-component form **agree on the matched-set graft** (why the MC passed) but
diverge on partial messages (the anchor), where ÷M is the mode the whole unified solver was built around.

**So the next derivation must:** (a) recast the graft/peel precision in the per-component density
parameterization above; (b) confirm the ÷M → `Var(log f_c)` conversion (the destination normalizer M_dst is a
common-mode across components — likely drops from the *relative* composition the ψ solve cares about, but this
must be shown, not assumed); (c) MC-validate the per-component form including the anchor limit.

## 4. σ²_transfer — the ground-up piece (owner: item 3, redo everything at once)

The current `σ²_transfer = var_proj[dst] + (μ_proj[dst] − μ_proj[src])²` is a **density-uniformity proxy**: it
penalizes a message by how far the two nodes’ *total densities* sit on the learned enrichment landscape. That
was a stand-in for the composition-transport variance we now mean to derive properly. Under the composition
solver the honest transfer term is **`Var(log r)`** — the variance of the enrichment ratio itself
(`enrichment_frame.composition_logvar` derives `Var(log ρ_tot)`, and `Var(log r) = Var(log ρ_tot^dst) +
Var(log ρ_tot^src)`). Key facts already established:
* On the **GRAFT**, enrichment cancels (T1), so `Var(log r)` does **not** enter — the transfer term should be
  ~0 there. Applying the current σ²_transfer on the graft is a **double-count**.
* On the **PEEL**, `Var(log r)` is **load-bearing** — it carried ~92% of the peel’s variance in MC.
* So the transfer term is **direction-dependent**, which the current single global σ²_transfer is not.
* **The transfer overhaul must cover the RELAY, not just the combine** (owner) — the relay accumulates context
  with the same damping; a correct transfer model replaces both in one pass.

## 5. Overdispersion — deliberately OUT (owner, 2026-07-24)

Model counts as **Poisson (ω = 0)**. Overdispersion is real on real data (measured ω ≈ 0.02, a depth-decaying
upper bound) but **we have no way to fit it** (no technical replicates) and **the synthetic suite is Poisson by
construction** so it cannot even validate an ω term (memory `synthetic_suite_is_poisson_omega_zero`). When
real data forces it, `1/n → 1/n_eff = (1+(n−1)ω)/n` is a one-line substitution. **No ω term, no config hook,
until it is proven necessary and fittable.** Precision is knowingly optimistic at high counts on real data.

## 6. The u-distribution decides the peel’s fate (measured)

`u = 1/(fraction that continues unspliced)`: real junctions have median u = 2.3, p75 = 5.3, p90 = 9.0. Even at
ω = 0.02 only ~40% of junctions (low-u, >½ continues) can be quantitative; u ≳ 3 junctions are hopeless at any
depth. So the peel is **u-weighted**, not a blanket channel and not globally demoted — use `ε` as a per-junction
precision weight, and treat `ρ_ν < 0` (a real **sampling** outcome, not an impossibility) as a truncation in
the PRIOR (like the intron factory’s `NegBinom·1[g≤C]`), not an emission gate.

## 7. What to carry forward vs discard

**Carry forward:** the per-component framework (§3); T1 (enrichment cancels on the graft); the item-E
share-weighting insight (a SUM → share-weighted variances); the peel’s difference-variance and u-weighting
(§6); the MC harness; the ω decision (§5).

**Discard / redo:** the ratio (k) parameterization of T2/T5 as the *implementation* form (keep only as the
matched-graft cross-check); the current σ²_transfer (§4); every archived precision doc
(`message_precision_derivation.md`, `transfer_variance_formal_derivation.md`,
`prediction_measurement_merge_derivation.md`, `graft_message_precision_design.md`,
`graft_precision_implementation_plan.md`, `spliced_precision_*`, `density_imputation_precision.md`, …) — all in
`archive/`, all superseded by this handoff.
