# The MESSAGE variance model — derivation (per-component density form)

**Status:** derivation for review, 2026-07-24. Built on the settled `τ_λ` FOUNDATION
(`variance_foundation_proposal.md`). Retires the density-uniformity `σ²_transfer` proxy
(`variance_model_handoff.md` §4) and the singular ratio-`k` parameterization (§3). Every law below is
Monte-Carlo-validated in `scripts/debug/message_variance_mc.py` (the arbiter; pure numpy, mirrors the template
`scripts/debug/message_precision_mc.py`). Read `ROADMAP.md` §3 and `variance_model_handoff.md` first.

---

## 0. What the message actually delivers, and in what units

The unified solver transports **per-component densities** `ρ_c` (c ∈ {g, R} at a face; the RNA splits to
{+,−} by the tilt). The ψ solve consumes, per component, a **Gaussian on `log f_c^dst`** with precision `p_c`
(`simplex_logodds._local_loglik_logodds`: `−½·p_c·(log f_c − mo_c)²`). The message mode is

```
    mo_c = log( ρ_c^src · r · E_c / M_dst ) ,      r = ρ_tot(dst)/ρ_tot(src)   (the reframe)
```

and the TRUTH is `f_c^dst = ρ_c^dst,true · E_c / M_dst`. **`M_dst` appears in both**, so it is COMMON-MODE and
cancels in the mode error:

```
    mo_c − log f_c^dst,true  =  log(ρ_c^src · r)  −  log(ρ_c^dst,true · r_true)
```

Under the imputation premise (neighbours share composition), holding the truth fixed,

```
    p_c = 1 / Var(log ρ_c^msg) ,    ρ_c^msg = ρ_c^src · r      ⟵ NO 1/n_dst, NO (1−f_g)² Jacobian.
```

This is the resolution of `variance_model_handoff.md` §3(b): **the ÷M_dst normalizer is common-mode and drops
out.** The message precision is `1/Var(log ρ_c^msg)`, and the whole model reduces to deriving
`Var(log ρ_c^msg)` per component, per direction. (Validated: M4/M5 below.)

The three components and their **source** log-variances (the "transport seed"):

| component | provenance | `v_c^src = Var(log ρ_c^src)` |
|---|---|---|
| gDNA `g` (imputed) | `f_g·M/E_g` | `Var(log f_g) + 1/n` |
| RNA-continue `ν` (imputed) | `f_R·M/E_r` | `Var(log f_R) + 1/n` |
| RNA-splice `μ` (**measured**) | `S/E_spl` | `1/n_spl` (composition-CERTAIN → no comp term) |

`Var(log f_g) = (1−f_g)²/τ_λ`, `Var(log f_R) = f_g²/τ_λ` — the FOUNDATION Jacobians of the single Schur scalar
`τ_λ`. The `1/n` is the **node's unspliced count** (`geometry.n_unspl_*`); g and ν **share it** (same `M`).
The measured spliced has its **own** count `n_spl` and zero composition variance (spliced IS pure RNA).

> **The count-zero-information split, made precise.** `M` (the node total) is common to g and ν, so it CANCELS
> in the composition ratio `ρ_g:ρ_R` and never votes the split — it only sets an overall scale. What survives
> into the composition is (i) the composition variances `Var(log f_c)`, and (ii) the **independent** spliced
> count `1/n_spl` (the only count not shared with the total). This is why the graft SUM (M2) carries the count
> as `w_μ²·1/n_spl` (spliced share only), and why the pure-gDNA anchor (M4) carries `1/n` (its gDNA IS the
> whole node, nothing to cancel against).

---

## 1. GRAFT — boundary→exon (a SUM), share-weighted [M1, M2]

The exon receives the boundary's WHOLE RNA flux, so the RNA density is a **sum** `ρ_R = ρ_ν + ρ_μ`
(imputed-continue + measured-splice), and the component sets match ({g, R}) so enrichment cancels (T1).
By the delta method on a sum (`d log ρ_R = Σ w_c d log ρ_c`, `w_c = ρ_c/ρ_R`, `Σw_c = 1`):

```
    Var(log ρ_R^msg) = w_ν²·v_ν + w_μ²·v_μ + σ²_transfer          w_ν = ρ_ν/ρ_R,  w_μ = ρ_μ/ρ_R
    Var(log ρ_g^msg) = v_g + σ²_transfer                          (gDNA transports alone)
```

The shares are **convex** (≤ 1): a minority component contributes quadratically LITTLE (the item-E
share-weighting, which falls out of the delta method — it is NOT the shipped replicate rule `pr += S`).
`σ²_transfer ≈ 0` on the graft (M5). **MC:** M2 rel 0.09–3.3% across `w_μ ∈ {0.85, 0}`; the naive unweighted
`v_ν + v_μ` is 10× too large at intermediate `w_μ`.

Degenerate check (no spliced, `w_μ = 0`): `Var(log ρ_R) → v_ν` exactly — the count-cancelling limit.

---

## 2. PEEL — exon→boundary (a DIFFERENCE), u-weighted [M3]

The boundary receives only what CONTINUES: `ρ_ν = ρ_R(x)/r − ρ_μ` — a **difference** (an absolute measured
density is subtracted, so enrichment does NOT cancel). With `T = ρ_R(x)/r`, `u = T/ρ_ν = 1/(fraction that
continues) ≥ 1`:

```
    Var(log ρ_ν^msg) = u²·Var(log T) + (u−1)²·v_μ ,    Var(log T) = Var(log ρ_R(x)) + σ²_transfer
```

A difference carries weights **≥ 1** (subtracting near-equal numbers destroys precision) — the mirror of the
graft's convex weights. `σ²_transfer` is LOAD-BEARING here (§4). **MC:** M3 rel 1.3% at `u≈2`; the
linearization under-states beyond `ε=√(fraction continuing) ≳ 0.15` (`u ≳ 3`, ~40% off at `u=10` — the
`variance_model_handoff.md` §6 "hopeless" tail, so the peel is **u-weighted**, not a blanket channel). gDNA at
a boundary just reframes: `Var(log ρ_g^msg) = v_g + σ²_transfer`.

---

## 3. σ²_transfer = Var(log r) — direction-dependent [M4, M5]

`r = ρ_tot(dst)/ρ_tot(src)` isolates the enrichment ratio `e(dst)/e(src)`; its uncertainty is additive for a
ratio of independent totals:

```
    σ²_transfer = Var(log r) = Var(log ρ_tot(dst)) + Var(log ρ_tot(src))
```

Each `Var(log ρ_tot)` is the EXISTING `enrichment_frame.composition_logvar`
(`= 1/n + [(1/E_g − 1/E_r)/B]²·Var(f_g)`) — already derived, already tested. **Direction (M5, rel 0.02%):**

* **GRAFT (matched set):** `r` multiplies BOTH components equally → common-mode → **cancels in the
  composition** → σ²_transfer ≈ 0. Applying the density-uniformity proxy here was a DOUBLE-COUNT.
* **PEEL / ANCHOR (partial message):** no matched partner to cancel `r` against → **load-bearing** (~92% of the
  peel variance in MC; on the anchor, `Δ(from r) ≈ σ²_transfer` to 0.02%).

This **retires** the density-uniformity proxy `var_proj[dst] + (μ_proj[dst]−μ_proj[src])²` and covers BOTH the
relay and the combine (one law, both places).

**Anchor limit [M4]:** a pure-gDNA source (`f_g=1`, `Var(log f_g)=0`) gives `Var(log ρ_g^msg) = 1/n +
σ²_transfer` — FINITE (rel 0.08%), where the ratio-`k` form was singular (`k = ρ_g/ρ_R → ∞`). The intergenic
anchor (the single most important source on unstranded data) is now regular.

---

## 4. ⚠ THE COMBINE — two independent messages are ≥2× (up to ~7×) OVER-confident [M6]

**A new finding, not in the handoff — SETTLED by the derivation workflow (`wf_c952640d`, 3 independent
derivations + 2 adversaries + adjudicator, all MC).** The composition is **one DOF** (`λ = logit f_g`). But the
÷M mode hands ψ **two** Gaussians from one source's belief — a gDNA message on `log f_g` AND an RNA message on
`log f_R = log(1−f_g)`. These are **RANK-1**: `d log f_g = (1−f_g)dλ`, `d log f_R = −f_g dλ`, both from the one
latent `λ`, so their correlation is **exactly −1** (proven to survive the spliced graft: all coefficient ratios
equal `−(1−a_g)/a_g`). ψ entering them as INDEPENDENT Gaussians sums their λ-Fisher:

```
    λ-Fisher from both = (1−f_g)²/v_g + f_g²/v_R = τ_T + τ_T = 2·τ_T      ⟹ EXACTLY 2× over-confident
    correct joint      = 1/Var(log k) = τ_T                              (one λ constraint = the joint truth)
```

**MC verdict (all agents + adjudicator):** the two-message stated variance is **≤ ½ the joint truth** —
**exactly 2× over-confident in the pure-composition limit, GROWING with the depth of the independent (spliced)
content, measured up to ~7× for a deep spliced measurement.** The "1.7–3×" of an earlier draft was a FLOOR, not
a range. Over-confidence is the count-zero-info FORBIDDEN direction (never confidently wrong), so this is a real
correctness defect — but a **bounded and localized** one:

* **INERT on the dominant-error unstranded/AMBIG arm.** There `τ_λ = 0` ⇒ both own-composition arms carry ZERO
  precision ⇒ no composition message is emitted at all (only count / spliced / transfer). The double-count
  cannot bite where there is no composition message. This is why the unstranded solve — the error mass — is
  unaffected, and why the mode has looked "largely correct."
* **Bites the already-healthy stranded / intron-factory arm** (`τ_λ > 0`), where a composition message IS sent.
  Installing the honest (smaller) M1–M5 variances makes those messages more confident, so the double-count can
  begin to regress the stranded arm — watch it in the A/B.

**THE SETTLED FIX — single-λ composition message (two-regime, fused precision).** Emit ONE composition
constraint per source: a Gaussian on `λ = logit f_g` with the FULL joint `Var(log k)` — the **spliced FOLDED
INTO R-total** via its share-weighted `w_μ²` count terms (§1), **NOT a separate/orthogonal spliced message** (a
separate spliced arm double-counts by 23%; dropping the ρ_g–ρ_R covariance under-states by 16–21% — both
REFUTED by the workflow). Two regimes, gated by the EXISTING `free_pos`/`free_neg`/`struct_lock`:
* **both components live** → single-λ via `Var(log k)` (the k-mode precision, count-cancelled);
* **structural anchor** (`ρ_R ≡ 0`, `k` singular) → the ÷M single-arm gDNA density mode (finite `mo_g`, M4).

The per-STRAND RNA split (pos vs neg) is a genuinely separate DOF (the tilt `θ`) and stays a separate message;
**only the g-vs-R-total pair is the rank-1 double-count.** This requires a **FUSED composition precision** (the
ρ_g–ρ_R cross-covariance = one λ constraint), NOT three independent per-component precisions — a genuine combine
**re-architecture**, not a surgical edit. It closes the foundation↔message unification
(`variance_foundation_proposal.md` §6: the same Schur formula IS the message-combine rule) and adds **no magic
number** — it REMOVES a factor.

**SEQUENCING (adjudicated — an engineering order, not a correctness compromise).** Land the validated core FIRST
— **M1–M5 + the `1/n` fusion/transport split + the `struct_lock` hard override** — since the double-count is
inert on the error-mass arm and the core is a large win over the density-uniformity proxy. Then build the
single-λ combine, validated by the SAME per-condition A/B. The single-λ combine is the **committed correct
target**, sequenced after the core — not deferred indefinitely, not implemented speculatively.

### Guards the settled laws still need (workflow, not re-derivation)
* **M1 near-wall minority component** — the log-Jacobian `Var(log f_c) = var/f_c²` UNDER-states by 36–92%
  (over-confident) for a small-AND-high-CV minority component against the `f_c→0` wall (measured `f_g=0.04`,
  CV 2.5 → 92%). Down-weight the minority-arm message there or document the knowing optimism — do not ship it
  silent. (The foundation's "10–18% near the vertex" understates this.)
* **M3 peel u/ε weighting is LOAD-BEARING** — a blanket peel is over-confident 27–40% at `u ≥ 3` (= p75 of real
  junctions). Use `u` (equivalently `ε = √(fraction continuing)`) as a per-junction precision weight, and treat
  `ρ_ν < 0` as a PRIOR TRUNCATION (like the intron factory's `NegBinom·1[g≤C]`), NOT an emission gate
  (`variance_model_handoff.md` §6).

---

## 5. What lands in code (Step 2 preview)

Pure layer (`enrichment_frame.py` → `message_arithmetic.py`), each with a closed-form unit test:
* `transport_seed_logvar(var_log_f, n)` → `Var(log f) + 1/n` (M1; the measured-spliced case is `n=n_spl`,
  `var_log_f=0`).
* `graft_sum_logvar(v_nu, v_mu, w_nu, w_mu, s2_transfer)` → M2 (share-weighted SUM).
* `peel_diff_logvar(var_log_T, v_mu, u)` → M3 (u-weighted DIFFERENCE); `var_log_T = Var(log ρ_R) + s2_transfer`.
* `transfer_logvar(logvar_tot_dst, logvar_tot_src, *, graft)` → M5 (`0` on the graft, the sum on the peel).
  Feeds off the existing `composition_logvar`.
* `divm_precision(var_log_rho_msg)` → `1/var` (M4 — the ÷M_dst conversion; no Jacobian, no `1/n_dst`).

Solver plumbing (`bp_solver._unified_solve`), from the FOUNDATION's deferred spec (plan v4 §4):
* **Fusion vs transport `1/n` split** — the fusion weight is `1/Var(log f_c)` (composition ONLY); the transport
  seed is `1/(Var(log f_c) + 1/n)`. Land byte-identically first, then flip the fusion.
* **`struct_lock` HARD OVERRIDE in `_fuse`** — a composition-certain node ADOPTS its own belief; never an `∞`
  weight (`1/0 = ∞ ⇒ nan` cascade). Interior-anchor no-nan test.

The invariants (`variance_foundation_plan.md` v4 §2 / handoff §7) are preserved: `N` enters only as power
(`τ_λ` or `1/n`); variances in log-odds / `Var(log f_c)`, never on simplex fractions; ONE Schur scalar `τ_λ`;
`struct_lock` = hard override.

---

## 6. Empirical results — Step 2b (M5 `σ²_transfer` wired) and the confirmed combine blocker

**M5 `σ²_transfer` landed** (`bp_solver._unified_solve`: the density-uniformity NPMLE proxy — which was
identically **0 in pass-0**, so pass-0 had NO transfer damping — replaced by the direction-dependent
`Var(log r)` from `composition_logvar`, 0 on the graft, load-bearing elsewhere). Full `ambig_dense_10mb` A/B
(`pass0_oracle_bench.py`, `OMP_NUM_THREADS=1`):

| condition | refit=0 base→M5 | refit=1 base→M5 |
|---|---|---|
| **aggregate** | 0.1267 → **0.1099** | 0.1234 → **0.0945** (= the legacy-with-factory target) |
| unstranded ss_0.50 | 0.1981 → **0.1328** | 0.1902 → **0.1068** |
| **unstr capON** (the landmine) | 0.2543 → **0.1761** | 0.3874 → **0.1626** (−0.225) |
| verystrong capture | 0.2985 → **0.1801** | 0.1514 → **0.1166** |
| **stranded ss_0.99** | 0.0262 → **0.0777** ⬆ REGRESS | 0.0293 → **0.0773** ⬆ REGRESS |

So `σ²_transfer` is a **large win on the error-mass arm** (unstranded × capture) but **regresses the stranded
arm** — a per-condition regression the gate forbids, so M5 is NOT gate-passing on its own.

**The regression is the TWO-MESSAGE combine, not `σ²_transfer` placement (both proven by experiment).**
* A "matched-exempt" variant (`σ²_transfer = 0` on every matched reframe, not just the graft) made stranded
  **worse** (0.0863 vs 0.0777) — so the honest σ²_transfer=0-on-matched placement (§3) is empirically
  anti-aligned with performance in the two-message mode, itself a symptom that the mode is wrong.
* Dissection (`pass0_node_dissect.py`, worst stranded scenario) shows the **strand-alone solve is good
  (boundary mwae 0.05) and the imputation MESSAGES corrupt it (→0.20)**, concentrated on **AMBIG boundaries**
  (`tau0_lam = 0`: the strand cannot pin f_g for a 2-DOF node at any κ, so messages are the only info there —
  and the two-message combine delivers wrong ones: e.g. an AMBIG boundary with oracle f_g=0.255, all RNA on
  the − strand, gets a precision-0 RNA− message and a confident wrong gDNA/RNA+ pair → solved 0.995).

**THE COMMITTED FIX (next task) — the single-λ combine, on a λ-space relay.** The composition is ONE DOF,
`λ = logit f_g`, and — the key architectural insight this session surfaced — **λ is ENRICHMENT-INVARIANT**:
under a matched reframe both ρ_g and ρ_R scale by the same `r`, so `λ = log(ρ_g/ρ_R) + log(E_g/E_r)` has `r`
cancel identically; the frame change only shifts λ by the known, enrichment-free eff-length ratio
`log(E_g^dst E_r^src / E_r^dst E_g^src)`. So the relay should carry a **λ-belief `(λ_mean, Var(λ))` per node
(+ a separate θ tilt belief for AMBIG)** rather than three per-component densities — which is both **simpler**
(no enrichment residual — the current per-component relay laments `Σ_c f_c = 75.6` at introns) and the **fix
for the M6 double-count** (one λ message, not two rank-1 copies). Each neighbour contributes one λ-constraint
with the **joint** `Var(log k)` (the T2 form `w_μ²(1/n+1/n_s) + (1−w_μ f_g)²·Var(λ_src)`, MC-validated in
`scripts/debug/message_precision_mc.py`), which needs the source's `Var(λ)` and spliced share `w_μ` carried
through the relay — so the λ-relay is what makes the joint precision computable without a cross-component
covariance. Two regimes: single-λ where both components live; the ÷M gDNA density mode at structural anchors
(`ρ_R ≡ 0`, k singular). This retires the two-message per-component combine entirely.
