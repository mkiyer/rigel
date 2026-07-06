# Adaptive-refine per-node solver — design & implementation plan

**Status:** design (2026-07-06). Not implemented. Successor to the fixed-grid solve
(`log_density_1d_solver_design.md`) that Fix 3 (`fix3_single_strand_grid_quantization.md`) tuned. Focus:
the **2-D AMBIG solver**, where the accuracy fix has the largest performance/memory consequences.

> **TL;DR.** The per-node solve grids the log-odds posterior and reads `f_g` as the median snapped to the
> grid. Fix 3 de-quantized the cheap **1-D single-strand** path with a finer grid (`sweep_n_grid_single_strand
> =256`). The same quantization afflicts the **2-D AMBIG** path, but there a finer grid is a genome-scale
> memory blow-up (`O(m·K·K_t)`). The adaptive-refine solver gets the **finer-grid median at coarse-grid
> cost**: a coarse grid locates the CDF-0.5 crossing, then the posterior is re-evaluated on a small local
> fine sub-grid there and the continuous median is read off. It keeps the **median** (vertex-safe,
> transform-invariant — the parabolic mode was rejected, see Fix 3), unifies the 1-D and 2-D paths, and — via
> the validated fact that **AMBIG `f_g` quantization is λ-driven, τ-independent** — de-quantizes AMBIG
> without a fine `(λ,τ)` cube. This is the solver to harden *before* the C++ optimization pass.

---

## 1. Problem

The per-node solve (`simplex_logodds`) evaluates the posterior `ψ` over the gDNA-vs-RNA log-odds `λ =
logit(f_g)` (1-D, single-strand) or `(λ, τ)` (2-D, AMBIG; `τ` = RNA sense/antisense tilt), then reads:

* `f_g` = posterior **median** over the τ-marginal λ-posterior, **snapped** to the grid point where the CDF
  first reaches 0.5 (`fg[idx]`).
* `f_pos`, `f_neg` = posterior **means** over the full 2-D posterior.

Two coupled issues:

1. **Quantization.** A deep-count node's posterior is far narrower than a grid cell, so it concentrates at
   the nearest grid point and `f_g` snaps to the log-odds lattice (Δf_g ≈ 0.085 at `K=60`). The true off-grid
   `f_g` differs by up to ±0.04; × high node mass this was the dominant post-Fix-1 error. **Fix 3** solved it
   for single-strand with a finer 1-D grid — cheap because that path is `O(m·K)`.

2. **The 2-D cost wall.** The AMBIG solve materialises an `(m, K, K_t)` cube — `O(m·K·K_t)` compute and
   memory. At genome scale (`m` ~ 10⁶ AMBIG nodes) a fine cube is prohibitive: at `K=K_t=256` one batch is
   `8192·256·256·4 B ≈ 2 GB` (18× the `K=60` cube; the retired 2-simplex lattice OOM'd on exactly this). So
   Fix 3 deliberately left AMBIG at `K=60` — leaving its `f_g` quantized.

**Goal:** de-quantize `f_g` for BOTH paths at ~coarse-grid cost, keeping the (correct) median estimator, so
the solver is both accurate and fast enough to commit to C++.

### Validated enabling facts

* **The median is the right estimator; the mode is not** (Fix 3). A parabolic sub-grid *mode* de-quantizes
  symmetric peaks for free but under-calls skewed/edge-piled posteriors — a confident pure-gDNA node piling
  up at `f_g→1` reads ~0.75 vs the correct >0.8. The median integrates that edge mass correctly and is
  transform-invariant (`median σ(λ) = σ(median λ)`). So the target is a **de-quantized median**, not a mode.
* **AMBIG `f_g` quantization is λ-driven, τ-independent** (measured: region 236 `f_g` = 0.490 at `(λ256,τ60)`
  vs 0.471 at `(λ256,τ256)` — τ resolution is irrelevant to `f_g`; only λ matters). Because `f_g` is read
  from the **τ-marginal** λ-posterior, we can de-quantize it by refining the **λ-marginal median alone**,
  leaving τ coarse. This is what makes a cheap 2-D fix possible.

---

## 2. The adaptive-refine median

The posterior over λ is a smooth function evaluable at any λ (strand Beta-Binomial + KDE global +
Jeffreys + cross-node messages + Jacobian — all closed-form or interpolatable, incl. `GdnaDensityPrior.
logpdf_kernel` at arbitrary `log ρ_g`). So we can afford to *locate* on a coarse grid and *refine* locally.

### 2.1 1-D (single-strand)

```
1. COARSE: evaluate ψ on the coarse λ grid (K_c, e.g. 48–60); normalize → p_c; CDF C_c = cumsum(p_c).
2. LOCATE: k = first index with C_c[k] ≥ 0.5  → median lies in the cell [λ_{k-1}, λ_k].
3. REFINE: re-evaluate ψ on a local fine λ sub-grid over a window around the crossing
   (e.g. [λ_{k-2}, λ_{k+1}], r ≈ 16–32 points). Normalize the window against the coarse tail masses
   (C_c[k-2] below, 1−C_c[k+1] above — the tails are far from the median, coarse is exact enough).
   Build the refined CDF and linearly interpolate the continuous λ* where it crosses 0.5.
4. READ: f_g = σ(λ*).  Vertex-safe: if the crossing cell is the boundary cell, the fine points resolve
   where within the edge the median sits (the piled mass is integrated, not extrapolated).
```

Cost: `O(K_c + r)` per node vs `O(K_fine)` for a globally fine grid — with `K_c=48, r=16` that is ~64 ψ
evals for finer-than-`K=512` median resolution. Replaces `sweep_n_grid_single_strand=256` with a coarser
grid + refine at equal-or-better accuracy and lower cost.

### 2.2 2-D (AMBIG) — the focus

`f_g` is the median of the **τ-marginal** λ-posterior `p(λ) = Σ_τ exp ψ(λ,τ)`. Since τ resolution does not
move `f_g`, keep the coarse τ grid and refine only the λ-marginal median:

```
1. COARSE CUBE: evaluate ψ on the coarse (λ,τ) grid (K_c × K_{t,c}); marginalize τ (logsumexp over τ)
   → coarse τ-marginal p_c(λ); CDF; LOCATE the λ-median cell (as in §2.1 step 2).
2. REFINE (λ only): for each of r fine λ* in the crossing window, marginalize τ ON THE COARSE τ GRID:
   p(λ*) = logsumexp_τ ψ(λ*, τ)   — r × K_{t,c} extra ψ evals (small).
   Build the refined τ-marginal CDF; interpolate the continuous λ-median λ*; f_g = σ(λ*).
3. RNA SPLIT: f_pos / f_neg remain the full-2-D posterior means on the coarse cube (τ resolution governs
   the RNA tilt, not f_g; refine only if the AMBIG RNA split later proves grid-limited — currently the
   AMBIG error is dominated by the κ strand-flip, a SEPARATE issue, not quantization).
```

Cost: coarse cube `O(m·K_c·K_{t,c})` (unchanged, the batched cube we already have) **+** refine
`O(m·r·K_{t,c})` (r ≈ 16 extra λ-columns). Memory: the coarse cube ONLY — no fine cube is ever
materialised. This is the crux: **fine-λ `f_g` accuracy with coarse-cube memory.** Contrast the naive fix
(fine `(λ,τ)` cube): `O(m·K_f·K_{t,f})` — the ~2 GB/batch genome-scale blow-up we must avoid.

---

## 3. Performance & memory

| approach | AMBIG compute / node | AMBIG peak memory | `f_g` accuracy |
|---|---|---|---|
| coarse fixed (today) | `K_c·K_{t,c}` | coarse cube | quantized (Δ≈0.085) |
| fine fixed cube | `K_f·K_{t,f}` (≈16×) | **fine cube (~2 GB/batch)** | fine |
| **adaptive-refine** | `K_c·K_{t,c} + r·K_{t,c}` (≈1.3×) | **coarse cube** | fine-λ (median) |

The refine adds a small local loop, not a cube-size multiplier, so it is genome-scale safe on both axes and
is the natural shape for the C++ port: the coarse cube stays the batched hot kernel; the refine is a short
per-node inner loop reusing the same Beta-Binomial / KDE / logsumexp primitives at a handful of λ*. The
single-strand 1-D path collapses to the same primitive (r-point refine, no τ). Expected net effect vs Fix 3:
**equal-or-better accuracy at lower cost**, and it removes the `sweep_n_grid_single_strand=256` grid.

---

## 4. Implementation plan (phased)

Each phase env-gated + in-process A/B (`total_threads=1`) against the fixed-grid baseline; ship only on a
clean full-suite result (mature accuracy flat-to-better, no FP regression, per-node Σ|err| down).

* **Phase 0 — expose ψ as a callable.** Refactor `_local_loglik_logodds` (1-D) and the AMBIG cube builder so
  the posterior can be evaluated at an arbitrary λ vector (reusing the strand mixture, sided-spliced,
  Jeffreys, global-`_regrid`-able, message terms). Prereq for refinement; no behaviour change.
* **Phase 1 — 1-D adaptive-refine median.** Add the coarse-locate + local-refine median to
  `_solve_nodes_logodds` behind `RIGEL_SOLVER=refine`. Validate it MATCHES `sweep_n_grid_single_strand=256`
  accuracy at a coarser `K_c` + `r` refine (target Σ|err| ≈ 16.8k on the benchmark condition) and PASSES the
  vertex unit tests. On success, make it the default and retire `sweep_n_grid_single_strand`.
* **Phase 2 — 2-D AMBIG λ-marginal refine.** Refine the τ-marginal λ-median in `_solve_ambig_logodds` (§2.2);
  τ stays coarse. Validate `f_g` de-quantizes to match a fine-λ cube (region-236-class nodes) at coarse-cube
  memory. **Requires the complex AMBIG benchmark (§6).**
* **Phase 3 — unify.** One refine routine drives both paths (1-D = 2-D with no τ). Delete the fixed
  fine-grid knob; the coarse `sweep_n_grid` + refine window become the only resolution parameters.
* **Phase 4 — C++ port.** Port the coarse cube + refine loop; the refine is a bounded local loop over the
  shared kernels. Only after the algorithm is frozen and validated.

---

## 5. Risks & open questions

* **ψ callable fidelity.** The refined ψ must be byte-consistent with the coarse-grid ψ (same terms, same
  global regridding). The KDE global at arbitrary `log ρ_g` uses `logpdf_kernel` (real tails) — already
  evaluable; keep the message log-fraction targets grid-independent (they already are).
* **Tail-mass normalization.** The refine window normalizes against coarse tail masses; verify the window is
  wide enough that the excluded tails carry negligible probability near the median (widen adaptively if the
  posterior is broad / low-count).
* **Determinism.** The locate + refine must be deterministic (no data-dependent branching that reorders FP
  sums) — the calibration determinism work (`calibrate_cross_process_nondeterminism`) applies.
* **The AMBIG strand-flip is NOT a grid issue.** Region 236's residual (`f_g`≈0.49 vs true 0.53 even at fine
  λ) is the κ-driven strand-attribution flip, a separate defect (the next dissection). The refine
  de-quantizes; it will not fix the flip. Do not conflate them when reading Phase-2 results.
* **Refine window vs multimodality.** A single crossing-cell window assumes a unimodal posterior near the
  median (true for the strand posterior). If a future prior induces multimodality, locate on the coarse CDF
  still finds the median cell correctly, but validate.

---

## 6. Benchmark prerequisite (deferred)

The current suite (`quick_3to1_5mb`) is **single-strand-dominated** — only a handful of AMBIG nodes (the
top-40 remaining error is 37 single-strand / 3 AMBIG). It cannot stress-test the 2-D solver: Phase-2's
accuracy/perf/memory claims need a scenario with **many AMBIG nodes** — dense overlapping opposite-strand
transcripts, nested/antisense loci, sparse-probe capture creating mixture regions, and enough AMBIG mass
that the 2-D path dominates runtime and memory. Building that scenario (via the simulator; see
`SIMULATOR.md` and the `scripts/sim` config) is the **prerequisite** for aggressively pursuing and
validating the AMBIG adaptive-refine, and is **deferred** for now. Until it exists, Phase 2 can be prototyped
and sanity-checked on the few AMBIG nodes present (e.g. region 236's neighbourhood) but not performance- or
accuracy-validated at scale.
