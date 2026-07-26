# The MESSAGE variance model — derivation (per-component density form)

**Status: COMPLETE and LANDED (2026-07-25).** M1–M5 + the M6 single-λ fix + **M7, the cross-cliff cost** (§7,
added 2026-07-25 — the term that finished the model). Every law is MC-validated in
`scripts/debug/message_variance_mc.py` (0 failures end-to-end). The shipped message precision is

```
    p = 1 / ( Var(log f_c^src) + 1/n_src  +  σ²_transfer  +  b̂² )
```

and no prior of any kind enters it. A/B: aggregate pass-0-vs-oracle 0.1267→**0.0969** (refit=0),
0.1234→**0.0828** (refit=1). §6's "next task" is done; the live handoff is `SESSION_2026_07_25_HANDOFF_6.md`.

Original header: derivation for review, 2026-07-24. Built on the settled `τ_λ` FOUNDATION
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

---

## 7. M7 — the CROSS-CLIFF cost: SCALE (`σ²_transfer`) vs COMPOSITION (`b̂²`)

M1–M5 price a message's **sampling** error. They do not price the imputation **premise** — "my neighbour and I
share a composition" — being false, and that premise is the message's entire content. That omission was the
solver's last real defect: a message crossing a 407× enrichment cliff arrived at full confidence and
steamrolled a weak-but-correct own belief (the anchor, node 1909: an oracle-0.985 mostly-gDNA exon collapsed
to 0.51 by a flanking boundary's spliced measurement).

### 7.1 The error splits into two orthogonal defects — exactly

Holding the truth fixed, the delivered mode error is (MC-validated to machine precision, M7a/M7b):

```
    mo_c − log f_c^dst,true  =  log( s_c^src / s_c^dst,true )  +  log( r̂ / r_true )
                                └── composition-SHARE mismatch ──┘  └── the reframe's own scale noise ──┘
    ⇒   σ²_c,delivered  =  v_src,c  +  σ²_transfer  +  b_c²
```

`s_c = ρ_c/ρ_tot`. Both the destination mass `M_dst` and BOTH effective lengths cancel identically — run the
MC with deliberately mismatched masses and eff-lengths on the two sides and the identity holds to 1e-16.

**This is why the shipped `σ²_cliff = (log r)²` proxy had to go.** It charged the *whole* cliff as mismatch,
which is conservative and did recover the stranded arm — but `log r` is dominated by the mass ratio, i.e. by
ENRICHMENT, so on the extreme-capture arm it over-damped messages whose composition was in fact preserved
across a 1000× density step. Pricing the two terms separately is what lets a pure-enrichment cliff cost
nothing while a composition-mismatch cliff costs everything.

### 7.2 `b̂²` prior-free: DerSimonian–Laird against the node's own self-solve

`b_c` is a POPULATION quantity — the third information source (`CALIBRATION_ARCHITECTURE.md`) — and pass-0 has
no population prior. But the destination holds an **independent** estimate of the same composition: its
message-free self-solve (`node_init`). Two estimators of one quantity is a two-study random-effects
meta-analysis, whose between-study variance has a closed-form method-of-moments estimator with **no tuned
constant** (the coefficient is 1 because it is a second moment; the truncation at 0 is the method's own):

```
    b̂_c² = max( 0 , G_c² − v_msg,c − v_own,c )     G_c = mo_c^msg − mo_c^own
    p_eff = 1/(v_msg + b̂²) = 1 / max( v_msg , G² − v_own )
```

with `v_own` the `τ_λ` FOUNDATION variance (`node_init.own_composition_logvar`), and `v_msg` the stream's own
variance — which already contains `σ²_transfer`, so the scale noise is correctly SUBTRACTED from the gap
rather than double-counted as mismatch.

**Three regimes, no gate, no constant** — they fall out of `v_own` alone:

| `v_own` | node | behaviour |
|---|---|---|
| finite (`τ_own > 0`) | strand- or factory-solved | a CONFLICTING message is killed, an agreeing one barely touched — **the stranded arm's fix** |
| `+∞` (`τ_own = 0`) | every AMBIG node; ALL unstranded data | `b̂² ≡ 0`, the message passes **bit-identically** — where messages are the only information there is no own opinion to contradict them, so **the M5 unstranded/capture win propagates untouched** |
| `0` (`struct_lock`) | structural pure-gDNA anchor | the full `G²` is charged (inert at the combine today — such a node is never `solvable`) |

**The safety property, exact.** `p_eff = 1/max(v_msg, G²−v_own)` ⇒ a message out-weighs the own belief **iff
`|G| < √2·σ_own`** — at every node, every composition, and independently of the source's depth. This is the
owner's "pass-0 must be weak and correctable" as an inequality rather than a knob, and it is the invariant to
preserve under any future change.

**Estimator quality (MC, M7c):** recovers the true `b²` to 0.1–0.4% for real mismatches (`b ≥ 1.5`); at `b ≈ 0`
it over-estimates by exactly `E[max(0,χ²₁−1)]·(v_msg+v_own) = 0.4839·(v_msg+v_own)` — the OVER-damping
direction, and harmless, because a message that agrees with the own belief moves the fused mode nowhere.

### 7.3 Where it is applied, and three choices that are load-bearing

`bp_solver._transport`, after `_pin_v`, on all three streams (mode-fusion, measurement, composition τ).

1. **POST-`_pin_v`.** The pin rescales the message to the destination's own mass, which removes its SCALE
   claim; what remains in `G` is purely a disagreement about the SPLIT. Verified: at matched composition the
   pinned gap is 0 to 2.4e-16 at cliffs `r = 1, 10, 407, 10⁴`. A pre-pin gap would re-import `log r` and
   recreate the `(log r)²` over-damping.
2. **The composition stream gets ONE λ-axis gap** (`G_λ = G_g − G_R`, `v_own,λ = 1/τ_own`), because it delivers
   one Gaussian on λ. A per-component gap there is a category error (λ-space variances minus a log-`f_c`-space
   gap) and, since `G_λ ≡ G_g − G_R` with opposite signs, satisfies `G_λ²/4 ≤ max_c G_c² ≤ G_λ²` — i.e. it
   under-damps by 1–4×, worst exactly when the gap straddles `f_g = ½`. Always the forbidden direction.
3. **A one-sided absence is the `b̂² → ∞` limit, carried as a MASK.** A source with no RNA at all (a pure-gDNA
   seam; a fully-consumed peel) asserts `f_c = 0` to ψ — the most confidently-wrong claim in the solver.
   Reporting the excess as `+∞` (⇒ precision 0) rather than as `log(t_c/_EPS)` keeps `_EPS` a numerical
   zero-test instead of promoting it to the constant that sets the surviving precision. (Measured identical to
   4 decimals against the floored variant across the suite — but the mask form has no constant in it.)

### 7.4 What M7 does NOT do

It is **inert wherever `τ_own = 0`**, which on unstranded data is every exon and every boundary (`i_strand ≡ 0`
by the κ=½ deadband, so `τ_own = τ_factory`, nonzero only at introns). So M7 is not what recovers the
unstranded/capture arm — retiring `(log r)²` is. The two effects are disjoint and the A/B separates them
cleanly: deleting `(log r)²` recovers verystrong (0.2779→0.1864 refit=0), `b̂²` recovers stranded
(0.0792→0.0390), and neither costs the other. Protecting AMBIG nodes requires a finite `v_own` there, which
only the trained gDNA hyperprior can supply — Phase 2 (`SESSION_2026_07_25_HANDOFF_6.md` §3).

---

## 8. M8 — the GRAFT's FRAME MISLIFT: the measured spliced is already in the destination's frame

**Status: DERIVED, MC-validated, implemented, A/B-won (2026-07-25).** Found by the single-strand × capture
study (`SESSION_2026_07_25_HANDOFF_8.md`). This is the term M5's graft exemption was hiding.

### 8.1 The defect

M5 sets `σ²_transfer = 0` on a graft edge because "`r` is common-mode across the matched set `{g, R}` and
cancels in the composition". **That cancellation is real** — verified to `1.8e-15` against the shipped arrays:
the delivered exon composition is exactly

```
    f_g : f_R  =  ρ_g^src·E_g^dst : (ρ_ν^src + ρ_μ)·E_r^dst
```

with the reframe `r` and the `_pin_v` factor `k` multiplying every component alike and dropping out. But the
exemption's *premise* is false for one of the three terms. `ρ_ν` is an imputed density that travels with
`ρ_g` from the same source, so `r` genuinely is common-mode for it. `ρ_μ = S/E_spl` is not: **a spliced
fragment's two blocks lie in the flanking EXONS**, so it measures the *destination's* RNA density directly.

Measured on the suite with oracle densities and no solver in the loop:

| | `ρ_R(exon)/ρ_spl(bnd)` | gDNA step exon↔boundary |
|---|---|---|
| capture OFF | 1.02 / 1.08 | **1.03** |
| capture ON | 1.16 / 1.86 | **6.1 / 6.8** |

The RNA channel is **capture-invariant** — the graft's own value needs no correction (the required frame
factor solved per edge has median `log c` = +0.009 / −0.008 / +0.054 off-capture). What changes under capture
is the *other* side of the ratio: the boundary's gDNA density falls 6× below the exon's, and since `r`
cancels, the graft edge never reframes it. The delivered claim is therefore wrong by exactly that step —
measured delivered bias −1.6 to −2.6 nats of RNA over-claim, against +0.4 to +1.0 off-capture.

So the grafted component has **no matched partner to cancel `r` against** — precisely M5's peel /
partial-anchor case, where `σ²_transfer` is load-bearing.

### 8.2 The law

```
    Var(log ρ_μ^msg)  +=  (log r)²
```

the method-of-moments second moment of a single observation of the un-cancelled frame step, with the model's
own `r` as the estimate of that step — the same logic as M7's `b̂²`, and **no tuned constant**. It is
identically 0 at `r = 1`, which is why the shipped model is exact off-capture. Applied to the spliced
measurement's OWN precision (`1/n_spl ⊕ (log r)²`), so M2's share weighting `w_μ²` arises implicitly from the
inverse-variance fusion with the correctly-framed `ρ_ν` arm: a graft that is a minority of the RNA is damped
proportionately less. **MC-validated** (`message_variance_mc.py`, M8): rel 0.24 % / 0.02 % / 0.01 % at
`r = 1 / 2.6 / 6.1`.

Calibration of the delivered `λ` against the oracle, on graft edges into exons —
`z2 = E[Δ²]/E[v]`, 1.0 = honest:

| | z2 shipped | z2 with M8 |
|---|---|---|
| gdna300 ss0.99 nrna_none capOFF | 62.0 | **3.4** |
| gdna300 ss0.99 nrna_none capON | 57.5 | **3.8** |
| gdna300 ss0.99 present capOFF | 58.7 | **3.7** |
| gdna300 ss0.99 present capON | 31.4 | **3.0** |
| gdna100 ss0.99 present capOFF | 111.5 | **2.1** |
| gdna100 ss0.99 present capON | 83.3 | **2.8** |
| gdna100 verystrong | 310.2 | **3.8** |

### 8.3 Why it is a VARIANCE and not a mode correction

The mode correction is available and it is §5.1's "graft outside the reframe" — and it was re-measured on
today's (post-relay-pin) baseline as part of this derivation. It wins more aggregate and **reproduces §5.1's
documented failure exactly**: it regresses every zero-gDNA condition (`gdna_none` capOFF 0.3668 → 0.4369,
capON 0.3780 → 0.4917), which is `density_composition_reconciliation.md` §3.3's forbidden direction (more
confidently *wrong*). M8 alone is the only variant that improves `gdna_none` rather than damaging it, because
it can only *remove* confidence — it never moves a mode toward a wrong answer. That is the governing
principle ("pass-0 must be WEAK and CORRECTABLE") deciding between two arms that the aggregate alone would
have decided the other way. Full ablation in `SESSION_2026_07_25_HANDOFF_8.md` §4.

**A/B (32 conditions, `OMP_NUM_THREADS=1`):** aggregate **0.0926 → 0.0900** (refit=0, 15 better / 9 worse /
8 flat) and **0.0779 → 0.0700** (refit=1, 20 better / 6 worse / 6 flat).

### 8.4 Known cost, and the next refinement

`r` conflates the capture step with genuine density structure — off-capture `r ≈ 1.15–1.3` boundary→exon
because the exon really does carry more RNA per base, while the true frame step is 1.03. So M8 over-damps
off-capture (charging `(log 1.2)² = 0.033` where the truth is 0.0009) and under-damps on-capture (`(log 2.6)²
= 0.91` where the truth is `(log 6.4)² = 3.4`). The measured cost is concentrated on **capture-OFF
unstranded**, the one arm where the graft is the only RNA information and damping it is pure loss:
`gdna100_ss_0.50_capOFF` 0.0866 → 0.1129, `gdna300_ss_0.50_capOFF` 0.0504 → 0.0631. Isolating the *enrichment*
part of `r` prior-free is the open refinement.


---

## 10. M10 — the PEEL as a COMPOSITION, not a subtraction (derived + MC-validated; NOT yet landed)

**Status: DERIVED and MC-VALIDATED (2026-07-25); the wiring is measured and does NOT land yet — the blocker
is the SHARE ESTIMATOR, not the law.** Owner-directed: the subtraction is the last absolute-density operation
left from the pre-composition solver and must be retired.

### 10.1 The law

Of the RNA at a seam, a fraction continues unspliced and the rest splices away. That is a **share**, and a
share is **enrichment-free** — capture multiplies both channels alike, so `e(B)` cancels identically:

```
    w = ρ_ν / (ρ_ν + ρ_μ)         ρ_ν = e·c_ν , ρ_μ = e·c_μ   ⇒  e cancels EXACTLY
    ρ_ν^msg = ρ_R(x) · r · w                                     a SCALING, not a difference
    Var(log ρ_ν^msg) = Var(log ρ_R(x)) + σ²_transfer + w_μ²·(v_ν^B + 1/n_spl)      w_μ = 1 − w
```

Delta method on `log w = log ρ_ν − log(ρ_ν+ρ_μ)`: both partials are `±w_μ`. The weights are **CONVEX**
(`w_μ ≤ 1`) — the exact mirror of M2's graft SUM — where M3's difference carried `u = 1/w ≥ 1` and amplified.
**MC (`message_variance_mc.py` M10): rel 0.20–0.96 %.**

### 10.2 Why it matters — the conditioning result (M10b, exact to 1e-12)

A systematic log-scale error `δ` in the reframed source arrives:

| through | delivered error | at δ = 0.30 |
|---|---|---|
| the SUBTRACTION `A − ρ_μ` | `log(u·e^δ − (u−1))` → `u·δ` as δ→0 | **1.77× / 2.39× / 5.01×** for u = 2/3/10 |
| the SHARE `A·w` | `δ`, exactly, for every u and every δ | **1.00×** |

This is load-bearing because the exon-facing reframe error is **irreducible**: a boundary samples an `fl_mean`
window around a point while an exon samples its interior, so with mid-exon probes they sit at genuinely
different capture — measured **0.4–1.3 nats, and unchanged when fed the ORACLE `f_g`**. Meanwhile `u`'s p75 on
real junctions is ≈ 3. The subtraction is the ONLY place in the solver where a scale error does not cancel
(`_pin_v` cancels `r` on 87.6 % of the error).

It also retires two known defects for free: `w ∈ [0,1]` so the peel can no longer go negative (the
zero-truncation that emitted "no RNA continues past here" at a live precision), and M3's ill-conditioned
`u ≳ 3` tail ceases to exist.

### 10.3 ⚠ WHY IT DOES NOT LAND YET — the share estimator

`A·w = A − ρ_μ` **exactly** when `w = 1 − ρ_μ/A`. So the subtraction already *is* a composition peel whose
share is taken from the exon's own claim. Switching to a share therefore only helps if `w` comes from
somewhere **better** — and the obvious source, the boundary's own self-solve, is not:

| condition | `w` from the boundary self-solve | `w` TRUE | \|Δ\| |
|---|---|---|---|
| gdna300 ss0.99 present capON | 0.531 | 0.435 | 0.115 |
| **gdna300 ss0.50 none capON** | **0.664** | **0.124** | **0.548** |
| gdna100 verystrong | 0.475 | 0.463 | 0.160 |

On unstranded data the boundary says 66 % of its RNA continues when 12 % does — because there its self-solve
IS the uninformed ~0.5 default. Wired that way the A/B loses: **0.0907 / 0.0707** vs HEAD's 0.0884 / 0.0678
(10 better / 13 worse at refit=0). A purely observable upper bound (treat all unspliced mass as RNA) is worse
still (\|Δ\| 0.282 / 0.626 / 0.142).

**Conclusion, and it is the next task's specification:** M10 is the correct vehicle and its conditioning
advantage is real and proven, but it must be applied AFTER the boundary has solved its own composition from
its three sources (strand tilt, the neighbouring intron's density deconvolution, the exon's gDNA claim).
Peel-by-composition and a properly-solved boundary are one piece of work, not two.

The pure laws are landed and unit-tested (`enrichment_frame.peel_continue_share`, `peel_share_logvar`);
only the solver wiring waits.
