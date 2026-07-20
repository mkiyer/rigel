# The NPMLE projection-variance model — design & production roadmap

**Status:** design, written 2026-07-16 (branch `calib-ambig-init-wip`), after the σ²_transfer derivation
([`sigma2_transfer_derivation.md`](sigma2_transfer_derivation.md)) and the projection-variance prototype
(`scripts/debug/npmle_projection_variance.py`, `npmle_deconv_plot.py`). **Owner-directed** (2026-07-16): build
the projection-variance model on our data; start from the single total-density NPMLE; make it the prior + the
message-variance facilitator for the *first* forward-backward pass; aim to **replace the painfully-derived
reference prior** ([`reference_prior_derivation.md`](reference_prior_derivation.md)) with the NPMLE. **Companions:**
[`CALIBRATION_ARCHITECTURE.md`](CALIBRATION_ARCHITECTURE.md) (count-zero-information — authoritative),
[`npmle_roadmap.md`](npmle_roadmap.md) (the shipped `GdnaRatePrior`; the KEY NEGATIVE RESULT this must not repeat).

> **The one object, three jobs.** A *single* NPMLE `P(log ρ) = Σ_j w_j·N(μ_j, h²)` fit **once** on the
> belief-free total density (`M/E`, `f_g=1`) does all three: **(1)** the gDNA prior (project `ρ_g=f_g·M/E`),
> **(2)** the message transfer variance `σ²_transfer` (project both edge endpoints; the mode-gap between them is
> the damping), **(3)** the reference — because a *proper* density with real tails is what the Jeffreys
> reference was standing in for. The count enters **nowhere** in (2): the message variance is the count-free
> mixture geometry, so a high count still buys no false confidence.

---

## 0. What the total-density NPMLE *is* — the enrichment profile (the only assumption-free invariant)

**We can assume NOTHING about the shape of the total-density landscape.** It is the convolution of the DNA and
RNA densities, and the mixture is unknown and unbounded: real RNA-seq spans no-DNA to heavily-DNA-contaminated,
so the modes, the spread, and the number of peaks are all circumstance-dependent. Do not build in any prior on
its shape.

**The one thing we CAN assume: hybrid capture is nucleic-acid-agnostic — a probe enriches whatever nucleic acid
is present, DNA and RNA alike.** So the total-density NPMLE, fit prior-free and belief-free at the very start,
is a **profile of ENRICHMENT**: its modes are probe-on / probe-off regimes, not DNA modes or RNA modes. That
profile is **robust to the mixture** — whatever species dominates, capture stamps the *same* enrichment
geometry onto it. This is the actual convolved prior, before deconvolution begins to separate DNA from RNA.

**Why this is the right foundation for `σ²_transfer` (job 2).** The dominant transfer discontinuity — the
crossing that must gag a message — is a **capture-enrichment boundary**, and because capture is agnostic, that
boundary is a discontinuity in **both** components at once. So the enrichment profile's crossing detection
applies to the gDNA message *and* the RNA message from one belief-free fit: exactly what pass-0 needs, before we
can tell DNA from RNA. (Residual RNA-specific spikiness — expression variation with no capture — is *not* in the
enrichment profile; it is refined in later per-component passes, §4. Off-capture this means the enrichment
profile can over-damp a gDNA message across an RNA-expression boundary; that is conservative, second-order, and
not the flagship regime — the flagship is unstranded **capture**, where total ≈ enrichment ≈ the shared
crossing.)

**And it justifies job 3 at pass-0:** since the enrichment profile is the landscape *both* components live on,
both the gDNA and RNA arms may legitimately take it as their pass-0 prior — retiring the Jeffreys reference from
the first pass, with deconvolution refining each arm afterward.

---

## 1. The object and its three roles

**Fit (exists — `GdnaRatePrior.fit`).** `P(log ρ)` on all nodes' total density `M/E` at `f_g=1`, τ=0:
a fixed-kernel Poisson-lognormal mixture, weights by deterministic EM. Mixture components: `μ_j` = grid points,
`w_j` = EM weights, `h` = bandwidth.

**Role 1 — gDNA prior (exists).** `logprior(fg_grid, M, E)` → `logP(log ρ_g)` at `ρ_g = f_g·M/E`, the ψ additive
term. Unchanged.

**Role 2 — message transfer variance `σ²_transfer` (NEW).** For an edge `src→dst`, project each endpoint's
observed density onto the mixture — responsibilities `r_j = softmax_j( log w_j − ½((logρ_obs−μ_j)/h)² )` — and
read the **predictive variance of the destination's rate given the source's message**:
```
    μ_proj(x) = Σ_j r_j^x·μ_j          (the denoised rate: snaps toward a mode)
    Var_proj(x) = Σ_j r_j^x·(μ_j−μ_proj)²  +  h²      (between-mode ambiguity + within-mode floor)
    σ²_transfer(src→dst)  =  Var_proj(dst)  +  ( μ_proj(dst) − μ_proj(src) )²
                             ╰─ dst ambiguity ─╯   ╰──── the CROSSING (mode gap) ────╯
```
*Measured behavior* (`npmle_projection_variance.py`, capON+nascent): same-mode edges → `σ² ≈ 2h² ≈ 0.1`;
a depleted↔enriched crossing → `gap² ≈ 38` — a ~200× separation, **invisible to the source alone** (the gap is
0 without the destination). This is the pair result of `sigma2_transfer_derivation.md`, re-derived from the
owner's projection idea and **smooth** (the responsibilities interpolate continuously through the valley, so no
step-function — the reviewer's Landmine 2 dissolves).

**Role 3 — the reference replacement (NEW; the prize).** The reference prior exists only because the RNA arm's
`logP_r` is unwritten, which makes ψ Haldane-improper at the vertices (`reference_prior_derivation.md` §4). A
**proper** NPMLE with genuine (non-clamped) tails supplies exactly the properness the Jeffreys reference
manufactured. So `logP_r` ← the NPMLE (the total-density fit at pass-0; the deconvolved RNA fit later) **retires
the `½·log(1−f_g)` reference term**. What must survive: the AMBIG **tilt** residual `R = −½·log(1−τ²)` is *pure
measure*, prior-independent (`reference_prior_derivation.md` §3.2) — it stays. The log-rate NPMLE cancels the
`f_g` Jacobian exactly (both arms), so no Jacobian is written (as `GdnaRatePrior.logprior` already notes).

---

## 2. The knob ledger — what is derived, what is a form-decision, what is the one real risk

| knob | status | how to set / test |
|---|---|---|
| **bandwidth `h`** | **the one real magic-number risk** | The crossing signal (`gap²` = mode separation) is **h-robust** — it is the distance between modes, not their width. Only the within-mode floor (`2h²`) and the fit's smoothness depend on `h`. So `h` is a **resolution knob** IF we verify accuracy is flat over a range: sweep `h∈[0.1,0.5]`, confirm (a) the oracle-strata reproduction and the flagship leak are stable, (b) the RNA/total fit is non-wiggly (needs `h≳0.2` — see §3). Accept a fixed `h` only if that holds; else adaptive `h` (§5, open). |
| **variance floor `ε_v = 2h²`** | **DERIVED** | The within-mode variance *is* the kernel resolution — a message can never be more precise than the landscape can resolve. This is the reviewer's Landmine-3 floor, and it is not a tuned ε. It also caps the max message precision at `1/(2h²)`, preventing rigid lock. |
| **projection obs-width** | tied to `h` | The observation is resolved at the landscape's own scale ⇒ obs-width `= h`. Not a separate number. |
| **count-free vs count-aware projection** | **DECISION: count-free** for `σ²_transfer` | Count-zero-information: the transfer variance is a world property, count-independent. The count enters *only* as the separate `1/M_src` term. (Count-aware projection is the source's *own-rate* precision — already covered by `1/M_src`; do not fold it into `σ²_transfer`.) |
| **the `σ²_transfer` formula (F1 vs F2)** | form-decision, prototype | F1 `Var_dst + gap²` (predictive) vs F2 `Var_src+Var_dst+gap²` (symmetric). Pick by which reproduces the oracle strata; both give the dominant `gap²`. |
| grid range / resolution, EM iters/tol | perf/numeric, derived | from the data support; unchanged from `GdnaRatePrior`. Not accuracy constants. |
| `calib_refit_iters` | existing | the deconvolution iteration (§4). |

**Net new magic numbers: at most one (`h`), and it is arguably a resolution knob.** The floor is derived; the
rest are form-decisions gated on the oracle. This is a far better ledger than the retired KDE (bandwidth +
mixture-bridge + trim + min-nodes) — and it is why this is worth doing.

---

## 3. The bandwidth ↔ variance coupling (the one hard interaction)

The RNA/total fit is wiggly at `h=0.15` (the production default) because that data is sparse and broad
(`npmle_projection_variance.py` panel A: `h=0.3/0.5` smooth it). **A wiggly NPMLE makes a wiggly `σ²_transfer`**
(panel B), because the projection reads the local mixture shape. So the bandwidth is not just a fit-cosmetics
knob here — it is load-bearing for Role 2. Resolution: `h` must be large enough that the *total-density* fit
(the belief-free substrate we actually use) is smooth (`h≳0.2` empirically), and the crossing must survive that
smoothing (it does — the modes are ~6 nats apart, `h=0.5` ≈ 1.15 nats). The sweep in §6 Phase 0 settles this.

---

## 4. The deconvolution iteration — refit as a deliberate, guarded direction (not a firewall)

The owner's "iterate and converge": start convolved (total), end deconvolved (gDNA + RNA).

```
pass 0:  fit P_total on M/E (belief-free).  Both arms' prior = P_total.  σ²_transfer from P_total's projection.
         → first forward-backward solve → belief f_g per node.
pass ≥1: deconvolve  ρ_g = f_g·M/E,  ρ_r = (1−f_g)·M/E.  Refit P_g, P_r.  gDNA arm ← P_g, RNA arm ← P_r.
         → recompute σ²_transfer PER COMPONENT from P_g / P_r.  → re-solve.  (calib_refit_iters)
```

**Why per-component refit is likely valuable (the motivation, confirmed on the oracle).** The convolved
`P_total` forces *both* components onto one crossing structure. That is right for gDNA and wrong for RNA, and
the oracle deconvolution (`npmle_deconv_plot.py`) shows exactly why: **gDNA is clean and bi/tri-modal** (capture
enrichment/depletion — a robust, predictable pattern that correlates with the probes), while **RNA is a broad
smear** over 3+ decades (expression-driven, a large dynamic range, *not* cleanly bimodal under capture). So
deconvolving gives **gDNA a crisp mode/crossing** (sharp σ²_transfer with a real gap) and **RNA an honest broad
variance** whose mode tracks expression rather than enrichment. `P_total` structurally cannot separate these;
`P_g`/`P_r` can. This is the strongest case for refit, and it is a per-component *variance* improvement, not just
a per-component *mode* one.

**The failure mode is real but not a law.** The earlier honest-σ²_imp refit *backfired* — but specifically
because it was measured on a **not-yet-peeled** (wrong) belief: adjacent wrong nodes agreed ⇒ σ² collapsed ⇒
confident wrong messages propagated (`npmle_roadmap.md` KEY NEGATIVE RESULT). That was a symptom of a broken
pipeline (bugs, bad math — much since fixed), not an inherent circularity. Refit-on-a-good-belief is honest
fixed-point iteration; refit-on-a-wrong-belief is self-reinforcing error. **Whether we are now in the first
regime is an empirical, oracle-testable question — not something to settle by assertion.** (Retracts an earlier
over-strong "firewall" claim.)

**So refit is IN, pursued deliberately, with three guards — measured, not assumed:**
1. **Belief-quality weighting (already exists).** `GdnaRatePrior` refits with each node's observation width
   `τ² = Var(log f_g)`, so weak/uncertain beliefs enter *wide* and cannot pin the fit — the self-tightening
   empirical-Bayes mechanism. The same width protects a per-component σ²_transfer refit.
2. **A belief-free floor/fallback.** The pass-0 `P_total` σ²_transfer is retained as the floor; a refit pass is
   *accepted* only if it does not regress the message-free self-solve or the oracle strata — so a divergent
   refit is caught and dropped, never silently taken.
3. **Oracle-gated convergence monitoring (in dev).** Across passes, does mwae / the σ²_transfer strata move
   *toward* the oracle (converge) or away (the old backfire)? This is the Phase-4 gate. The toys let us watch it
   directly.

If refit converges cleanly (plausible now, given the fixes), it is the better model and becomes the default. If
it diverges on a component, that component falls back to the belief-free `P_total` σ². The point is we **decide
this with the oracle**, not up front.

---

## 5. Open design decisions (to prototype in Phase 0, before wiring)

1. **`h` fixed vs adaptive.** Prefer fixed (simplest, deterministic) if the §3/§6 sweep shows flat accuracy.
   Adaptive (`h` ∝ local sparsity) only if a fixed `h` cannot be both smooth and crossing-preserving.
2. **The `σ²_transfer` formula** (F1/F2, and whether `Var_dst` belongs or only `gap²`) — pick against the
   oracle strata.
3. **Reference replacement scope at pass-0.** Does the RNA arm use `P_total` at pass-0 (retire the Jeffreys
   reference immediately), or keep the Jeffreys reference until `P_r` is deconvolved (pass-1)? The former is
   cleaner but asserts RNA looks like total (≈ gDNA under capture — rough). Prototype both against the
   L-invariance + zero-gDNA-pin tests.
4. **The percentile readout** (p10–p90) the owner mentioned — a *diagnostic* on the projection, not a core knob;
   the 2nd moment `Var_proj` is the operative summary. Keep percentiles for the QC/diagnostic surface.

---

## 6. Roadmap — prototype → production

> **What the toy suite is for (owner, 2026-07-16).** These synthetic datasets are a **development vehicle**, not
> the target — they are far-end (300 % gDNA), small, and not representative. The goal now is **the algorithm and
> the code**: get the mechanism right, the wiring right, the bugs out, everything running smoothly, with toy
> accuracy as a *confidence check that we built it properly* — NOT as the thing to optimize. **All real knob
> tuning (the bandwidth above all) happens later, on real data.** So `h` is picked "smooth enough" by eye on the
> toys and left; do not over-invest in locking it now.

**Phase 0 — belief-free harness + a smooth-enough `h` (NEXT; no production code).** One prototype that, per
condition: fits `P_total`; computes `σ²_transfer` per edge via the projection; compares to the oracle σ²_transfer
strata (`transfer_variance_diag.py` is the oracle). Pick `h` "smooth enough" (the RNA/total fit is non-wiggly,
`h≈0.25`) and confirm the crossing survives it. **Gate:** the projection `σ²_transfer` reproduces the oracle
ordering (dep-dep small, crossing large) belief-free; settle the *form* decisions §5.1–5.2 (which are about
correctness, not tuning). Deliverable: the mechanism validated end-to-end + the form-decisions locked. *(The
`h`-sweep is a light sanity check that accuracy is flat, not a tuning exercise — the real value settles on real
data.)*

**Phase 1 — the projection module (standalone, unit-tested).** A pure function/class:
`NpmleProjection(prior)` with `.project(M, E) → (mode, var)` and `.edge_sigma2(src, dst) → σ²_transfer`.
Deterministic, belief-free, no solve. **Unit tests** (portable, in `tests/`): synthetic bimodal mixture →
valley variance spike; well-separated modes → crossing `gap²`; single mode → `≈2h²`; floor honored. This is the
count-zero-info-clean primitive.

**Phase 2 — wire `σ²_transfer` into the message precision.** In `bp_solver._scan`, the four `pr=` sites
(gDNA :462, RNA :479/:492/:499): `pr = 1/(v + 1/M_src + σ²_transfer[edge])`, the edge σ² precomputed once
(belief-free) before the sweep. Keep the mature-gate-dismantled state. **Gate (physics, not suite-MAE):** the
flagship unstranded+capture leak improves; the message-free self-solve is a **floor** the messages never break;
high-nascent neutral; 5-pass refit stability (Landmine-2 check). Pin the four src hashes before/after.

**Phase 3 — replace the reference prior with the NPMLE.** Retire the `½·log(1−f_g)` RNA reference in favor of
`logP_r` (= `P_total` at pass-0 per §5.3, `P_r` after deconvolution); keep the tilt residual `R`. **Gate:** the
`reference_prior_derivation.md` §10.7 tests — L-invariance (a proper NPMLE must give a finite L→∞ limit), the
zero-gDNA pin, no vertex blow-up; goldens regen with an attributable diff.

**Phase 4 — the deconvolution iteration.** Refit `P_g`/`P_r` on the peeled belief (§4), σ²_transfer fixed
belief-free. **Gate:** convergence over 5 passes (mwae stabilizes, no oscillation); the flagship + the 24-cond
oracle benchmark.

**Phase 5 — real data.** LBX0190 (pristine cfRNA, no nascent) stops injecting gDNA; a vcap capture library's
enriched boundaries are not crushed. Only then does this become the default path.

---

## 6.5 Phase 0 RESULTS (2026-07-16, `scripts/debug/phase0_projection_sigma2.py`)

Belief-free projection σ²_transfer vs the oracle **raw** gDNA disagreement variance (the apples-to-apples
comparison — both carry sampling noise; F1≈F0≈F2, the `gap²` dominates), 7 conditions, `h`-swept.

1. **Core mechanism VALIDATED in the flagship (capON).** The projection reproduces the true (raw) gDNA transfer
   variance nearly exactly and **`h`-robustly**: dep-dep 0.33/0.33, enr-enr ~6–7 vs 6.2–6.8 (F1), MIXED ~10–13
   vs 10–13. The calibration curve tracks `y=x` above the floor, slightly **conservative** (safe). *From one
   belief-free total-density fit, the projection predicts the per-edge transfer variance where it matters.*
2. **`h` is EMPIRICALLY a resolution knob.** The enr-enr / MIXED σ² (load-bearing) is `h`-invariant and matches
   the oracle across `h∈{0.15,0.25,0.35}`; **only the floor `= h²` moves with `h`** (0.12→0.33→0.66, exactly
   `(h·ln10)²`). Confirms §2: the accuracy-relevant signal is `h`-robust; `h` only sets the floor.
3. **The floor `h²` IS the count-zero-information cap** — a message can never exceed precision `1/h²`, no matter
   the count. This is not an ad-hoc `ε_v`; it is the principled "high count ≠ high precision" ceiling, and it
   equals the landscape resolution. (Supersedes the earlier `ε_v=2h²` framing.)
4. **Form decision SETTLED: F1** (F0/F2 differ by <5% — the `gap²` dominates).
5. **Off-capture, the CONVOLVED σ² over-damps gDNA — refit is REQUIRED, not optional (empirical proof of §4).**
   True off-capture gDNA variance is 0.02–0.04 (uniform); the convolved projection gives 0.3–1.8, because the
   total-density landscape off-capture is RNA-expression, not gDNA enrichment, so it manufactures gDNA
   "crossings" from RNA structure. The **deconvolved gDNA NPMLE is unimodal off-capture ⇒ no crossings ⇒ σ² = the
   floor** — the refit removes the RNA-contamination component (0.3–1.8 → ~floor), giving the gDNA channel its own
   (clean) landscape while the RNA channel keeps its smear. This is exactly the owner's DNA-clean/RNA-smear
   argument, now measured: **the convolved pass-0 σ² is a good flagship facilitator but structurally over-damps
   gDNA off-capture, and only per-component deconvolution fixes it.** (The floor `h²` remains — it is the
   count-zero cap, not a defect.)

**Open (Phase 2):** F1 matches the *raw* oracle (includes Poisson); the message formula also carries `1/M_src`,
so σ²_transfer mildly double-counts the source Poisson. Second-order (the floor dominates same-mode; the crossing
is Poisson-insensitive) — revisit at wiring.

---

## 6.6 Phase 1+2 BUILT (2026-07-16) — the projection is in the real solve

**Production code (surgical, contained to two files):** `GdnaRatePrior` now stores the mixture `weights` and has
`.project(mass, eff) → (mu_proj, var_proj)` (chunked, belief-free); `bp_solver.node_sweep` computes the
projection once from the NPMLE and threads a `transfer_variance` toggle; `_scan` adds
`s2t = var_proj[dst] + (mu_proj[dst]−mu_proj[src])²` to all four `pr` denominators (the one-line `v → v+s2t`).
`calibrate` already passes `gdna_prior` through — no orchestration change. **Tests:** 3 new projection unit tests
(floor `= h²`, crossing dominates within-mode, determinism) + full calibration suite green (215 passed, 1 xfail;
the 1 red `test_gdna_sweep_zero_gdna_pin_and_monotone` is the pre-existing known-red artifact — **proven inert to
this change**: it runs `gdna_prior=None` ⇒ `s2t=0` ⇒ byte-identical `pr`; `git stash` reproduces the same
`0.6833`).

**A/B on the real solve (`scripts/debug/sigma2_transfer_ab.py`, per-region mwae vs oracle, projection ON vs the
shipped σ²=0):**

| regime | mwae OFF → ON | verdict |
|---|---|---|
| capOFF (all 4) | 0.099→0.033, 0.093→0.034, 0.034→0.013, 0.039→0.010 | **helps (big)** |
| capON **stranded** (ss0.99) | 0.297 → **0.038** | **helps (huge)** |
| capON **unstranded** (ss0.50, the flagship) | 0.484→0.631, 0.446→0.579 | **regresses** |
| mean | 0.213 → 0.191 (−0.022) | net better |

**Reading it.** Honest crossing-damping helps wherever the node has a fallback (capOFF: the prior/self-solve;
stranded capON: the strand — a **huge** 0.30→0.04). It **regresses the unstranded flagship**: there `I(f_g)=0`
(no strand) and the pass-0 prior is deliberately weak (n_eff≈0.15), so the node was leaning on the
*over-confident* crossing messages that σ²_transfer correctly gags — and the net gDNA under-call worsens
(−3.4M→−4.5M). This is not a bug (the mechanism is honest and stranded capON proves it works); it is the
unstranded regime being **information-starved**, and it points squarely at the fix: the **deconvolution refit**,
which gives the enriched exon a strong enriched-gDNA *prior* (a clear mode) to replace the message help the
gagging removed. The A/B thus independently motivates Phase 4 — exactly the movie the owner wants to watch.
(Mid-implementation; production accuracy is not yet the target.)

---

## 7. Why this is worth the build (the case for the elegance)

* **One fit, three jobs** — collapses the gDNA prior, the message damping, and the reference into a single
  object with one bandwidth, replacing (prior KDE + floor + global + Jeffreys reference + a separate σ²_imp).
* **Count-zero-information by construction** — the message variance is mixture geometry; the count is only
  `1/M_src`.
* **The reference-prior replacement** retires the hardest, most fragile derivation in the calibration
  (`reference_prior_derivation.md`'s 734 lines of Berger–Bernardo) with a proper, informative, data-driven prior.
* **The magic-number ledger shrinks** — from the KDE's four knobs to one (`h`, arguably a resolution knob) plus
  a *derived* floor.
