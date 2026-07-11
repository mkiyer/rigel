# Fix 2 — interdependent-simplex message precision (SHELVED; reimplementation reference)

**Status:** implemented + empirically evaluated 2026-07-06; **NET-NEGATIVE, not shipped.** The full
implementation is preserved on branch **`fix2-joint-simplex-precision`** (this doc is the reference to
resume it). It remains a *potential* weakness of the calibration solver (see §5), so we keep the design.

**One-line summary:** representing a node's per-component precision as three independent marginal variances
lets a message on an *uncertain* pie slice move a *confident* coupled slice (the "trojan horse"). The joint
(conditional-variance) fix down-weights such messages — but it is net-negative because for a strand-null
AMBIG node the neighbour messages ARE the load-bearing signal, and you cannot tell a *good* disagreeing
message from a *bad* one at the destination. The trojan-horse defense belongs at the SOURCE (= Fix 1, the
mixture bridge, which stops nodes from emitting the bad message), which subsumes Fix 2's goal.

Companion: `boundary_kde_valley_collapse_and_simplex_precision.md` (§3 states the problem; the top UPDATE
block records the ship decision). Diagnosis memory: `node1085_kde_valley_boundary_crush_diagnosis.md`.

---

## 1. The problem (why Fix 2 exists)

The calibration solve deconvolves each node's unspliced mass into a pie `(f_pos, f_neg, f_g)` on the
2-simplex (`f_pos + f_neg + f_g = 1`). The per-node belief stores three **marginal** posterior variances
`Var(log f_pos)`, `Var(log f_neg)`, `Var(log f_g)`. Messages between nodes are per-component log-fraction
Gaussians whose precision is the disagreement-aware `σ²_edge` (`dispersion_aware_message_precision.md`),
anchored on the destination's per-component **marginal** message-free belief.

**The trojan horse (verified on region 1085, gdna300/ss0.99/capture-on):** after the pass-2 KDE, node 1085
is *confident* about gDNA (`f_g≈0.90`, `Var(log f_g)=0.017`) but *marginally uncertain* about the RNA split
(`f_neg≈0.08`, `Var(log f_neg)=3.36`). A neighbour (a collapsed boundary) sends an RNA− message
(`f_neg≈0.29`, precision 22.6). The per-axis disagreement check passes — `residual²=1.67 < Var(log f_neg)=3.36`,
so `s²_edge≈0` and the message keeps full precision — because the *marginal* f_neg variance overstates
f_neg's freedom. But accepting `f_neg=0.29` forces `f_g` down from 0.90 to 0.62 through the simplex
constraint, and f_g is confident. **The per-component surprise metric never sees the induced move on the
confident coupled slice.** The least-precise slice is a backdoor to the most-precise one.

## 2. The implemented fix (what is on the branch)

Measure disagreement against the **conditional** variance (given the confident coupled slice) — the diagonal
of the joint precision `Σ⁻¹` — instead of the marginal. Files/edits (all env-gated `RIGEL_JOINT_DISAGREE=1`,
default off ⇒ bit-exact baseline):

* **`strand_deconv.py`** — `NodeDeconv` gains three optional log-fraction cross-covariance fields:
  `gdna_pos_cov`, `gdna_neg_cov`, `pos_neg_cov` = `Cov(log f_a, log f_b)` over the per-node posterior.
* **`simplex_logodds.py`** — both solvers compute the cross-covariances from the grid posterior:
  * single-strand (`_solve_nodes_logodds`): `Cov(log f_g, log f_active)` on the 1-D λ ridge, routed to the
    live strand (dead strand cov = 0; `cov_pn = 0`).
  * AMBIG (`_solve_ambig_logodds`): all three, from the 2-D `(λ,τ)` posterior (`log f_g` is τ-independent,
    broadcast; `log f_pos/f_neg` vary over the cube).
  * `_solve_nodes_logodds_all` scatters them into the returned `NodeDeconv`.
* **`bp_solver.py`** — `node_sweep._local_solve` returns the covs; the sweep computes, **for AMBIG
  destinations only**, the conditional variance `cv_c = max(Var(c) − max_o Cov(c,o)²/max(Var(o), δ), δ)` (the
  pairwise conditioning on the most-constraining coupled slice; `δ=1e-6` is the user-prescribed
  full-rank/rank-deficiency floor), and uses `cv_*` (not the marginal `v_*`) as the `s²_edge` anchor in the
  three message-emission blocks in `_scan`.

**Restricted to AMBIG destinations on purpose:** a single-strand node's live strand is intrinsically coupled
to f_g (`f_active = 1−f_g`), so its conditional variance is ~0 *tautologically* (not from confidence);
applying the conditioning there would wrongly crush the load-bearing nascent-in-intron imputation. The trojan
horse is exclusive to AMBIG nodes (KDE pins f_g while a strand-degenerate RNA slice stays marginally free).

## 3. Why it is net-negative (the empirical + conceptual result)

4-arm in-process A/B (target condition, per-node calibration Σ|err|; deterministic):

| arm | Σ\|err\| | region 1085 (true 0.920) |
|---|---|---|
| baseline | 34,875 | 0.624 |
| Fix 1 only | 31,020 (−11%) | 0.927 |
| **Fix 2 only** | **38,569 (+10.6%)** | 0.766 (partial defend) |
| both | 32,906 (worse than Fix 1 alone) | 0.927 |

Deliverable (gene quant): Fix 2 on top of Fix 1 changes the false-nascent by only −216 (Fix 1 alone does
−25,434), i.e. negligible.

**Conceptual reason (the important part for any reimplementation):** for a strand-null AMBIG node the
neighbour messages are the *primary* signal — there is no strand to fall back on. Fix 2 makes a confident
AMBIG node skeptical of every *disagreeing* message, but most disagreeing messages are the **correct**
cross-node gDNA-continuity correction of a wrong local belief, not trojan horses. At the destination you
cannot distinguish the two — the local belief (the thing you'd judge the message against) is precisely what
is unreliable at an AMBIG node. So destination-side skepticism blocks more good than bad. The trojan-horse
defense therefore belongs at the **source** (don't emit the bad message), which is exactly what **Fix 1**
does by un-collapsing the boundary that emitted it. Fix 1 subsumes Fix 2's goal for the current failure mode.

## 4. When to resume Fix 2

Fix 1 removes the *known* source (KDE-valley boundary collapse). Fix 2 becomes worth revisiting if a
**different** source of confident-but-wrong messages appears that Fix 1 does not cover — e.g. a real-data
regime where AMBIG nodes still receive high-precision conflicting messages after Fix 1. Symptoms to watch:
per-node calibration Σ|err| dominated by AMBIG exons whose KDE-corrected local belief is overturned by a
single neighbour message, where the source node is NOT a collapsed boundary (so Fix 1 does not help).

## 5. What to try differently (a reimplementation should not repeat the destination-skepticism trap)

* **Source-side confidence gating (preferred).** Instead of the destination distrusting disagreeing
  messages, make the *source* emit a weaker message when its belief is prior-driven rather than
  evidence-driven. The message precision already includes the source belief variance (`base_var = var[src] +
  pois`); the gap is that a KDE/prior-pinned source can be *confidently wrong*. A source-confidence term that
  distinguishes "pinned by my own strand+counts" from "pinned by the population prior" would down-weight the
  bad messages without touching the good ones. (Fix 1 already prevents the specific pathological pin.)
* **Full joint (Mahalanobis) combination**, not the pairwise-conditional approximation: represent the belief
  as a Gaussian on the simplex in ILR/ALR coordinates with a full 2×2 precision, combine messages by
  precision-matrix addition, and measure disagreement as the Mahalanobis distance under the full precision.
  This is the principled version, but note it will hit the SAME conceptual wall (§3) if applied at the
  destination — it is only worth it combined with source-side gating.
* **Numerical robustness (carry over):** strict interior simplex clamp `[δ, 1−δ]`, δ≈1e-7, before any
  log-ratio/precision math (ILR fallback if ALR blows up at vertices); a full-rank `δI` baseline seeded into
  the local precision so summed rank-deficient (single-strand / strand-blind) messages cannot leave the
  precision matrix singular. The shelved implementation already floors denominators at δ=1e-6.

## 6. Reproduce / restore

* Branch: `git switch fix2-joint-simplex-precision` (Fix 1 env-gated + Fix 2 intact).
* Toggle at runtime: `RIGEL_JOINT_DISAGREE=1` (on that branch, or if the code is restored).
* A/B tool: `scripts/debug/pass_trace.py` (per-node stage trace) — the 4-arm table in §3 is
  `RIGEL_MIX_BRIDGE ∈ {0,0.01} × RIGEL_JOINT_DISAGREE ∈ {0,1}`.
