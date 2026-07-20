# NPMLE gDNA prior — implementation roadmap

**Status:** the core is SHIPPED to branch `calib-ambig-init-wip` (2026-07-14). The calibration gDNA prior —
the scalar depleted floor **and** the density KDE — is replaced by one **count-space Fixed-Kernel
Poisson-lognormal Mixture NPMLE** (`calibration/gdna_rate_prior.py::GdnaRatePrior`). It is fit ONCE before
the sweep on ALL nodes' total unspliced density (`f_g=1`, belief-free), is extremely weak (`n_eff≈0.15`
pseudo-obs), and projects onto each node's ψ via `GdnaRatePrior.logprior` (real tails + one Jeffreys).
`gdna_density_prior.py` + `bp_solver._gdna_seed_estimate`/`_floor_estimate`/`_global_logprior`/`_kde_logprior`
(~400 lines) are DELETED; 1089 tests green. Fitter/struggles: `npmle_struggles.md`; theory:
`gdna_prior_zero_handling.md`.

## Remaining work (tracked; 2026-07-14)

1. **Perf** (IN PROGRESS — dev blocker). The FIT is cheap on small suites (~0.35s, and minor genome-scale).
   The **genome-scale bottleneck is the SWEEP**: `simplex_logodds._solve_nodes_logodds_all` ≈ 42s of a 55s
   `calibrate` (`_solve_ambig_logodds` the 2-D (λ,τ) cube 18.5s, `_local_loglik_logodds` 11.8s,
   `_regrid_global` 4.2s, `_scan` 6.6s), called 3× (init + local + final). Levers: reuse the message-free ψ
   across the local+final solves (they recompute strand+global+regrid, differing only in messages — ~2× if
   split); cache the anchored regridded global once; shrink the AMBIG cube. Fit already optimized
   (grid 300→200, iters 200→150 w/ tol 1e-5; ~1.8× on small).
2. **The refit loop — SHIPPED (`calib_refit_iters=1`), a large win.** pass-0 (gentle: weak prior + the
   RNA-inflated scalar σ²_imp ⇒ messages nudge, never ruin) → solve → **refit `P(ρ)` on the peeled belief**
   (`ĝ=f_g·count`, `τ²=Var(log f_g)`) → re-solve. Oracle A/B over all 24 ambig conditions (mean `mwae_fg`):
   **0.151 → 0.118 (−22 %)**; ss0.50 0.274→0.218; ss0.99 0.028→0.018; and the worst corner
   (zero-gDNA/unstranded/nascent/capture) **0.52 → 0.28 (−47 %)**.

   **KEY NEGATIVE RESULT — an "honest" per-node message precision BACKFIRES.** `message_precision.
   adjacent_imputation_variance` (σ²_imp(ρ) = a fixed-kernel regression of the Poisson-corrected adjacent
   disagreement on the source's *solved* log-rate; deterministic, no spline) was built and ablated:
   * prior refit only → mean mwae **0.118** ✅
   * honest σ²_imp only → **0.189** ❌ (worse than baseline; the worst corner 0.52 → **0.95**)
   * both → 0.156 (the effects cancel)

   Why: measured on a **not-yet-peeled** belief, adjacent nodes agree on their *wrong* gDNA rates ⇒ σ² small
   ⇒ messages turn confident and **propagate the error**. *Honesty measured against a wrong belief is not
   honesty.* Keeping messages weak lets the prior refit do the work. Any future non-constant precision must be
   **gated on belief quality**, not merely measured. (`message_precision.py` is retained but UNUSED — delete or
   revisit with a quality gate.) The DNA+RNA belief-evolution diagnostic (both `P(ρ)` per iteration) is still
   worth building on top of the refit loop.
3. **The over-call / unstranded-RNA peeling (the main accuracy problem).** → **See
   `message_precision_and_dof.md`** (the consolidated review document: the trace showing the local PRIOR — not
   the messages — drives the over-call, the measured degrees-of-freedom double-count, and the open questions). The 24-cond oracle benchmark:
   stranded (ss0.99) near-perfect (mwae ≤ 0.08); the OLD gDNA→RNA under-call collapse is GONE. But UNSTRANDED
   (ss0.50) now OVER-calls gDNA — worst on zero-gDNA (mwae 0.42–0.58: RNA hallucinated as gDNA). Cause: with
   no strand, RNA is peeled only by messages, and (a) the message precision `σ²_imp` is too weak, (b)
   unspliced/nascent RNA has no clean message route. **Message precision `σ²_imp` is the agreed next target.**
4. **Docs (N4 tail).** `CALIBRATION_ARCHITECTURE.md`, `calibration_prior_production_reference.md`, this
   roadmap's body (§0–§6 below still describe the pre-ship plan), and repo `CLAUDE.md` still describe the
   retired 2-pass/floor/KDE — rewrite to the shipped pass-0 NPMLE regime.

The sections below (§0–§6) are the PRE-SHIP plan, kept for provenance; the shipped design differs (Fixed-Kernel
Poisson-lognormal Mixture, all-nodes total-density substrate, single pass) — see `npmle_struggles.md` §8–9.

---

## 0. The one-paragraph shape

The prior lives in three layers. NPMLE **replaces the estimation layer**, **extends the projection**, and
**reuses the solver unchanged**:

| layer | symbols | action |
|---|---|---|
| **estimation** | `_gdna_seed_estimate`, `_LogLinearVarMean`, `_fit_seed_varmean`, `_floor_estimate`, `_global_logprior`; `GdnaDensityPrior` + `build_training_substrate` + bandwidth machinery | **replace** (one NPMLE) |
| **projection** | `_kde_logprior` (bp_solver.py:400) — `log ρ_g = log(f_g·M/E)` → evaluate `logP(ρ)` on the f_g grid → `(n_nodes,K)` additive term | **extend** (swap the density evaluator for `np.interp` on the NPMLE grid) |
| **solve** | `node_sweep`, `_local_solve`, `_scan`, `_comb`, write-back, `chain_*_deconv`, `adjacent_disagreement_variance` (σ²_imp, **not** the prior) | **reuse** (untouched) |

The prior enters the solver at exactly two `_local_solve(global_lp, …)` calls (bp_solver.py:740 local,
:899 final); NPMLE only changes **how `global_lp` is built** (the block at bp_solver.py:664-724), never the
`(n_nodes,K)`-additive contract `simplex_logodds._solve_nodes_logodds_all` consumes (regridded for
single-strand via `_regrid_global`). That contract is the clean seam — the solver is not touched.

---

## 1. Performance — solved by design (the key result)

Naive per-node Poisson-EM would be a real blocker (~10⁶ confident nodes × 200-grid × 120 iters × refits ≈
30–60 min). **Two structural facts collapse it to a non-issue:**

1. **`(count, E)` collapse.** The EM responsibility `r_ij` depends only on `(k_i, E_i)`. Counts are
   tiny-integer-dominated (64–94 % at k=0, 90–99.6 % at k≤3) and boundary E is near-constant (~fragment
   length); region E bins cleanly in log space. Unique `(k, log-E-bin)` cells collapse to **~10³–10⁴** (from
   ~10⁶) — 100–1000×. Run the EM over **weighted cells** → `O(iters·cells·grid) ≈ 2.4e8` → **sub-second**,
   bit-exact up to a log-E binning error ≪ the log-ρ grid step. *Implement the collapse FIRST — it is not an
   optimization, it is the difference between prohibitive and free.*
2. **Belief-free fit set.** The confident selector reads only `statics.u_pos/u_neg` + `geometry.eff` — **no
   `f_g`** — so the *structural* fit does not change across passes. The **projection** reuses `_kde_logprior`'s
   tabulate-and-interpolate lattice verbatim (built once per pass, outside the `_scan` loop), so the NPMLE
   never multiplies the sweep's dominant cost.

**No numba/C++ needed for the NPMLE.** The one genuine sweep hotspot — the Python `_scan` forward/backward
loop (bp_solver.py:801-885, run 4×) — is **pre-existing and orthogonal**; do NOT push any NPMLE evaluation
into `_scan` (keep it a static precomputed array).

**Immediate perf hygiene** (from the prototype, mechanical): hoist the iteration-invariant `logpois` out of
the EM loop; use `simplex_logodds._lse` not `scipy.logsumexp`; **drop the RNG subsample** (the collapse
removes the need — and RNG is a cross-process nondeterminism hazard).

---

## 2. Implementation phases

**N1 — the NPMLE fitter (standalone, in `gdna_density_prior.py`).** A fixed log-ρ grid + the `(count,
log-E-bin)` weighted-cell **collapse** + deterministic Poisson-EM (warm-startable) + a `logP(log_ρ)`
evaluator. **Consolidate the confident-node selector here** — one function returning `(count, eff,
node_kind)` incl `count=0`, replacing the three overlapping masks (`_gdna_seed_estimate.struct_seed`,
`gdna_density_prior._clean_exon_boundary`, `scripts/debug/pass0_kde_zero.confident`). Test against the
prototype on a cached sample (P(ρ) agreement; collapse ratio).

**N2 — the projection (extend `_kde_logprior`).** Swap `gdna_prior.logpdf_kernel` → `np.interp` on the NPMLE
grid (NPMLE already lives on a grid, so the lattice-tabulation branch simplifies); **delete the
mixture-bridge block** (NPMLE is zero/valley-native); apply the **Jeffreys term exactly once** (see §3).

**N3 — wire into `calibrate` / `node_sweep`.** Fit the NPMLE at **pass-0** (structural, solve-free) and pass
it to both sweeps; the **evolving** step adds *solved* exon nodes to the fit set at pass ≥ 1 (sub-second via
the collapse — see §4 open question). Delete `_gdna_seed_estimate`/`_floor_estimate`/`_global_logprior` from
the `global_lp` block; `global_lp` flows unchanged into `_local_solve`. Resolve the cap-semantics question
(§3).

**N4 — cleanup + migration (see §5).** Delete the dead estimation functions; rewrite `GdnaDensityPrior` →
the NPMLE object; adapt `diagnostics.py::from_prior` + the `_capture` hook; rewrite
`test_gdna_density_prior.py`; migrate config; update the docs.

**N5 — validation (the gate; behind a flag until green).** Golden regen (`--update-golden`); **region A/B**
on the stranded/gDNA sim suites (no regression); the **AMBIG + unstranded+capture** benchmark
(`calib_pool_benchmark` — the 13 M gDNA→RNA leak and the documented **0.02→0.38 AMBIG drift**); and the real
caches (does the enriched under-call vanish, does pristine LBX0190 stop injecting gDNA).

---

## 3. Blockers & how each is handled

- **[HIGH] Jeffreys double-count.** `−log(1−f_g)` appears in `_global_logprior` (:387, floor nodes) *and*
  `_kde_logprior` (:470, non-floor), de-duped only by `floor_mask`. Unifying to one model dissolves
  `floor_mask`, so the two sources merge. **Fix:** fold Jeffreys into the single NPMLE projection, applied
  uniformly; delete both originals; add a regression assert that each node's ψ carries exactly one
  `−log(1−f_g)`.
- **[HIGH] Loss of the floor's downward anchor.** The scalar floor pins gDNA-poor AMBIG nodes → ~0; a prior
  KDE-replaces-floor attempt drifted **0.02→0.38** (false-positive) because the per-node KDE's RNA-residual
  mode sat above `ρ_global`. **The NPMLE zero-atom (count-0 low-rate mass) must reproduce this** — the
  **primary make-or-break**. **Fix:** keep the confident set rich in `count=0` structural nodes so the
  zero-atom carries real mass; N5 verifies AMBIG-poor → ~0 on the benchmark.
- **[MED] Cap semantics.** Today is hybrid: non-floor nodes get a `_GLOBAL_STAB_PREC`-capped hyperprior
  ("never overrides a node's own strand"); floor nodes get an *uncapped* data-driven density likelihood. The
  NPMLE is one term — **decide: capped hyperprior or full likelihood?** Leaning likelihood (it is the
  population model fit on confident nodes), but then the cap's stability role for low-count nodes must be
  re-derived. Prototype both on the quick benchmark.
- **[MED] Enriched-mode coverage at pass-0.** The structural confident set excludes exons (they carry RNA),
  so pass-0 P(ρ) may under-represent the exon-interior enriched mode until the pass≥1 re-fit adds *solved*
  exons. **Fix:** the evolving re-fit (§4) — and validate the enriched mode appears on a capture condition.
- **[LOW] Removed-symbol references.** `_capture` (bp_solver.py:951-973), `diagnostics.py::from_prior`, and
  `test_gdna_density_prior.py` reference `rho_global`/`rho_floor`/KDE fields. **Fix:** adapt `from_prior` to
  the NPMLE (`x_grid`/`logP_grid`/`modes` survive; `bandwidth`/`train_x` retire), fix `_capture`, rewrite the
  test module.

---

## 4. Design decisions (RESOLVED 2026-07-14)

**D1+D2 — fit ALL nodes with the current belief; precision enters as observation WIDTH, not weight; no cap.**
Decisions 1 (evolving fit set) and 2 (precision) are the same question. Resolution: **always fit `P(ρ_g)` on
ALL nodes** using the current belief (exons are the largest measurement source — essential). Safety comes
from **honest precision applied as the deconvolution WIDTH of each node's unit-mass observation** —
`σᵢ(log ρ_g) = sqrt(Var(log f_g) + 1/(f_g·k))` (the existing `build_training_substrate.log_rho_std`):
- Precision as **width** (not a per-node weight) preserves each node's unit mass but spreads it — so it does
  NOT bias the distribution shape (this is what sidesteps the §8e precision-weight objection, which was about
  *weights*), while giving the weak-start/tighten behavior.
- **Weak start**: early diffuse beliefs (large `Var(log f_g)`) ⇒ wide observations ⇒ broad `P(ρ)` ⇒ weak
  prior ⇒ cannot pin nodes to confident wrong values. **Tightens** as beliefs converge (`Var→0` ⇒ narrow ⇒
  sharp `P(ρ)`). Empirical-Bayes, self-tightening.
- ⇒ **Full likelihood, no `_GLOBAL_STAB_PREC` cap.** The cap was a band-aid for a dishonest prior strength;
  honest width replaces it. Fitting ALL nodes subsumes the confident-first route (that is this with a hard
  0/1 belief weight).
- **Zero anchor**: wide-σ (count-0) nodes are uninformative in log-rate space (`log 0=−∞`), so the "gDNA is
  scarce" anchor comes from the **count-space Poisson** treatment of low-count nodes (`Poisson(0|ρE)=e^{−ρE}`
  correctly prefers low ρ, bounded by ~1/E). → the hybrid, D3.

**D3 — the grid: log everywhere; exact-low-count / continuous-high-count (the same hybrid as precision).**
- `P(ρ)` on a **log-ρ grid**; range **derived** from the observed `(k,E)` support (`1/E_max` → max cell
  rate) — not a magic number.
- Fit collapses to **`(count, log-E-bin)` cells**: count **exact for low k** (≤ `K*`; the dominant 90–99 %
  and where discreteness bites; also the zero anchor via the Poisson path), **log-binned/Gaussian for high k**
  (rare; `Poisson→Gaussian` so lossless); **log-E bins** (`log λ = log ρ + log E` smooth in `log E` ⇒ binning
  error ≪ the log-ρ step).
- **Precision hybrid = grid hybrid**: low-count/short-E → exact count-space **Poisson** (zero anchor);
  high-count/long-E → continuous log-rate **deconvolution** (belief+count width). One threshold `K*`.
- Thresholds (`K*`, log-E bins, log-ρ grid size) trade **memory/compute, not accuracy** (finer = slower, not
  more correct, above a floor) — the robust parameterization the "no magic numbers" rule permits (they are
  perf knobs, not accuracy constants). Retire the prototype's un-derived `0.1/E_max`, `dens_max×3`, `w[:3]`.

**The remaining hard piece is a PROTOTYPE, not a decision:** cleanly fuse the count-space Poisson (zero +
low-count, exact) with the log-rate deconvolution (belief+count width, high-count) into one fit, and verify
on the caches across pass 0→1→2 that `P(ρ)` starts weak (broad) and sharpens ONLY as beliefs converge, keeps
the zero anchor, and never pins prematurely. This precedes N3 wiring.

---

## 5. Cleanup, config, tests (the two review lenses that failed structured-output — covered here)

**Dead code / shrinkage (bp_solver.py loses ~200+ lines):** delete `_gdna_seed_estimate`,
`_LogLinearVarMean`, `_fit_seed_varmean`, `_floor_estimate`, `_global_logprior`, `_KDE_LATTICE_PTS_PER_BW`
usage folded into the NPMLE eval. `gdna_density_prior.py`: `GdnaDensityPrior`/`build_training_substrate`/the
bandwidth estimators (`_silverman_bandwidth`, `_lscv_bandwidth`, `_weighted_median`, `_weighted_kde_logpdf`)
→ replaced by the NPMLE fitter + the consolidated confident-selector. Net: the calibration package **shrinks
and simplifies** — one prior model, one selector, one projection.

**Final module layout:** NPMLE fitter + confident-selector live in `gdna_density_prior.py` (renamed in spirit
to the gDNA-rate prior); the projection stays in `bp_solver._kde_logprior` (extended); the solver is
untouched. Stays within the package discipline (≤ ~25 modules, ≤ ~8 constants).

**Config (`config.py::CalibrationConfig`):** retire `gdna_prior_bandwidth`, `gdna_prior_mixture_bridge`,
`calib_kde_bridge_trim_pct`, `calib_kde_min_training_nodes`. Add `npmle_grid_size`, `npmle_em_iters`,
`npmle_em_tol`, `npmle_log_e_bins`, `npmle_min_confident_nodes` — with a `__post_init__` validator mirroring
the existing one. **Keep `sweep_n_grid` (=60)** — the f_g solve grid, orthogonal.

**Tests / docs:** rewrite `tests/calibration/test_gdna_density_prior.py` for the NPMLE; update
`test_priors.py` / `test_bp_solver.py`; regenerate goldens (behavior change). Update
`CALIBRATION_ARCHITECTURE.md` and `calibration_prior_production_reference.md` (the prior section) + the
docs/README index.

---

## 6. Validation gates (acceptance)

1. **Unit:** NPMLE `P(ρ)` matches the prototype on a cached sample; collapse is bit-exact up to log-E binning.
2. **No-regression:** golden regen reviewed for direction; region A/B on stranded/gDNA sims ~neutral.
3. **The make-or-break:** on the AMBIG + unstranded+capture benchmark, gDNA-poor AMBIG → ~0 (no 0.02→0.38
   drift) **and** the enriched under-call (the 13 M leak) shrinks.
4. **Real data:** pristine LBX0190 stops injecting gDNA; the vcap enriched boundaries are no longer crushed.

Only after 1–4 does the flag flip and the KDE/floor path get deleted.
