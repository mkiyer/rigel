<!-- title: Calibration — next steps (FL deferred) -->
# Calibration — next steps (FL implementation deferred)

State as of 2026-06-20, after the honest-precision arc + the adversarial readiness review.

## 1. Where we are (re-baselined against the SHIPPED production path)
- **The honest-precision `(density, precision)` message is already shipped** (`dfc3360f`): `_message` carries
  the source density `ρ_src` and a count-and-length-aware precision `τ_ρ = 1/(σ²_bio + (n+a)/E²)`; the
  `(M/E)²` factor is the necessary density→fraction conversion, bounded by the load-bearing binomial floor.
- **The count-space-SOLVE rewrite (Dirichlet `Σα·log f` pseudo-counts) was explored and REJECTED** (6-lens
  review): a one-sided `α·log(f)` is a monotone push to the f=1 vertex (measured: μ=0.2→f_g≈0.77, opposite of
  content); it contradicts the user-locked `(density, precision)` decision; and its motivating regression
  (node463/stranded-0 +28K) was a `RIGEL_GLOBAL_ALLNODES`-ON (default-off) **artifact**. The flag-on prototype
  is **reverted**; production code is unchanged.
- **Production calibration is solid:** stranded-0 phantom = **0** (the shipped intrinsic-solve `225d65ac`),
  capture gDNA ≈ 359,766 (oracle f_g=0.84; the EM recovers it off-capture), unstranded-0 ≈ 34,886.
- **The one open calibration phantom — unstranded-0 at κ≈½ (≈34,886) — is an IDENTIFIABILITY FLOOR**, not a
  precision bug: with no strand information, gDNA is indistinguishable from balanced RNA. The only resolver is
  the **fragment-length signal** (gDNA fragments run longer than RNA) → the **FL track, DEFERRED**.

## 2. The global gDNA prior — status (valid; precision honest-enough; do NOT rewrite now)
**Valid.** The production global (the `else` path, flag off) — a Gaussian pull of the AMBIG/intergenic nodes
toward the self-solvable mean, with the var~mean spread — **works**: it gives stranded-0 = 0 and the correct
capture baseline, as a weak tie-breaker the per-fragment EM then refines. The reverted prototype (exposure-
pooled mean + Gamma-on-count, flag on) was the rejected count-space-solve exploration.

**Honest precision?** *Partially, and adequately for production.* The global's precision shares the same
`(M/E)²`-Gaussian form as the messages — it is **not** the fully count-space form. It works because the global
and message `(M/E)²` factors **mask each other** on production (both bounded; the intrinsic-solve handles the
strand-decisive nodes). The genuinely-honest global (the exposure-pooled mean — which converges to 0 in a
zero-gDNA library by intergenic-length dilution — plus a count-space precision) is **real and recorded**
(`global_prior_count_space_derivation.md`) but is part of the **deferred** count-space solve: it cannot ship
piecemeal (a count-space global with `(M/E)²` messages unmasks the message imbalance — measured).

**Revert / rewrite?** Revert: **done** (prototype removed). Rewrite: **not now.** The global's honest-precision
upgrade is bundled with the messages into the deferred count-space solve, and both inherit the **σ²_bio
honesty** fix below (the shared input). The exposure-pooled-mean insight is preserved for that work.

## 3. The next implementation step — the σ²_bio raw-density honesty fit
This is the **one remaining honest-precision piece** that is shippable on its own and improves **both** the
messages and the global (they share `σ²_bio`).

**The problem.** `σ²_bio(μ)` (the between-node biological density spread; `fit_gdna_varmean` /
`fit_rna_varmean`) is fit each pass on the **frozen previous-pass snapshot densities**, which at init are the
all-gDNA seed — so the spread is artificially tight in the early passes that set the basin → `τ_ρ` (hence the
message/global precision) is over-confident exactly when it most shapes the solution. `N_src` / the message
precision is only as honest as `σ²_bio`.

**The fix.** Fit `σ²_bio` on the **raw observable per-node densities** (`ρ_obs = M/E`, gДNK-channel
f_g-independent) — a stable, non-circular axis — rather than the relay-smoothed posterior densities. (OQ4 from
the FB review; `count_space_execution_plan.md` §1c.) The var~mean machinery (`MonotoneVarMean.fit_offset`) is
unchanged; only its **input** changes from snapshot-density to raw-density.

**Scope:** `fit_gdna_varmean` / `fit_rna_varmean` inputs in `bp_solver.node_sweep`; no new constants; keep the
1-pseudo-observation stabilizer (`_MSG_PSEUDOCOUNT`). Likely a small, in-place change.

**Honest expectation:** production is already solid, so this may move the headline metrics little — it is a
**correctness/robustness** fix (the precision becomes honest, removing the early-pass over-confidence and
making any future count-space work safe). **Measure** the per-condition signed phantom/siphon before claiming a
win; if it is neutral, it still de-risks the deferred count-space solve.

## 4. Validation harness (build BEFORE any precision change)
The flag-on/flag-off prototype debugging showed how easy it is to chase a non-production artifact. Build a
**faithful, reproducible** harness first (extend `scripts/debug/gmean_study.py` — it already has the by-origin
oracle geometry):
- **signed** per-condition node-level **phantom AND siphon** (not just |phantom|), vs the oracle;
- the **production default path** baseline as a hard assertion: stranded-0 = 0, unstranded-0 ≈ 34,886,
  capture ≈ 359,766;
- a **capture-unstranded** condition (gdna300/ss0.50) — the current 3 study scenarios miss it;
- the post-EM `net_gdna_to_rna` gate (the calibration change must propagate through quant, not just the prior);
- forbid blanket `--update-golden`: per-golden diff + one-line justification.

## 5. Deferred tracks (parked, with their design docs)
- **FL signal** (gDNA-FL vs RNA-FL per-node) — the unifying resolver for *both* the unstranded-0 κ≈½
  identifiability floor (§1) *and* the EM-blind-under-capture deficit (the stress-test dominant error). The
  single highest-leverage real-data lever; deferred by decision.
- **Full count-space solve** (exposure-pooled global + count-space messages, done **together** after σ²_bio
  honesty) — keeps the `(density, precision)` carrier; the Dirichlet form is rejected
  (`count_space_solve_impl_plan.md`, superseded). Only pursue if a measured non-uniform-gDNA need appears.

## 6. Higher-ROI direction (once calibration precision is locked)
The stress test (`net_flow` decomposition) showed the **dominant real-data error is EM-side, not calibration**:
off-capture the EM recovers a 2×-wrong prior to near-perfect, but **under capture it goes blind** (the gDNA
component eff-len vs the RNA/nascent eff-len construction). That — and the FL signal that resolves it — is the
bigger prize than further calibration-prior precision. Recommended sequence: lock calibration (§3+§4), then
take up the EM-side eff-len / FL track.
