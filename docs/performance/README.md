# Rigel performance & memory optimization

Status of the compute/memory optimization effort on the `rigel quant` pipeline, the profiling
infrastructure to reproduce it, and the ranked roadmap for further work.

**Scope note.** This is *engineering* performance work — making the existing algorithm run faster and
in less memory **without changing its scientific output** (verified against the golden regression at
`rtol=1e-6` and the 16-condition accuracy benchmark). It is deliberately separate from the *calibration
accuracy* effort (the redescending-prior / nascent-hallucination redesign), which changes the answer.
When those two touch the same code, accuracy wins take priority and perf work re-bases onto them.

---

## TL;DR — where we are (as of 0.6.2)

The first real-data run (human GENCODE index, a 971k-fragment targeted library) **crashed trying to
allocate 2.69 TiB**. It now runs end-to-end. Cumulative, all output-preserving:

| Metric | Before | After 0.6.2 |
|---|---|---|
| Outcome on real data | **OOM crash** | completes |
| Peak RSS (971k-frag human lib) | 2.69 TiB attempted | **8.6 GB** |
| Calibration stage wall time | ~110 s | **~66 s (−40%)** |
| End-to-end `rigel quant` | (crashed) | ~105 s |
| Deconvolution accuracy | — | unchanged (goldens `rtol=1e-6`; benchmark within run-to-run noise) |

The pipeline is now **fast and small enough for routine production** (snakemake, per-sample, ~96 GB
nodes). The next tier of wins exists but needs either a compiled inner loop or more mixed-precision —
see [Roadmap](#roadmap--biggest-opportunities).

---

## Why this became urgent

`rigel` was developed and validated on small synthetic suites (a 5 Mb genome, ~1000s of fragments).
Cost there is dominated by fixed overhead, so the *scaling* behaviour was never exercised. The first
whole-human-genome run exposed it immediately: the calibration solves over **every** region/boundary
node in the index (~1.5M nodes for GENCODE), **independent of read depth** — a tiny targeted BAM still
pays the full genome-scale calibration cost. Two things blew up:

1. **A quadratic-in-disguise allocation** — the gDNA-density KDE prior evaluated `P(log ρ_g)` at every
   `(node × solve-grid)` point (`log_rho` = **90,335,640** values) and built one dense
   `(n_eval, n_train)` = `(90M, 4096)` float64 matrix = **2.69 TiB**.
2. **Genome-scale per-node work** in pure Python/numpy that was fine at 1000s of nodes but not at 1.5M.

---

## What shipped (0.6.2)

All four are **output-preserving** (byte-identical where noted, else within the golden tolerance).
See `CHANGELOG.md` [0.6.2] and the commits below.

| # | Optimization | Where | Impact | Precision | Commit |
|---|---|---|---|---|---|
| 1 | **KDE query-axis tiling + tabulate-and-interpolate** | `gdna_density_prior._weighted_kde_logpdf`, `bp_solver._kde_logprior` | 2.69 TiB → bounded; the fix that stops the crash | tiling byte-identical; tabulation ≪ bandwidth (gated to large queries; small/golden inputs take the exact path) | `2e93021a` |
| 2 | **Lean numpy `logsumexp`** (drop `scipy.special.logsumexp`) | `simplex_logodds._lse` | ~15 s off calibration | matches scipy to ~1e-15; **byte-identical** on the f64 path | `2e93021a` |
| 3 | **In-place `psi` accumulation + CSE** | `_solve_ambig_logodds`, `_mixture_strand_loglik` | ~6% off the AMBIG solve, fewer live cube temporaries | **byte-identical** (`np.array_equal` on an 8192-node cube) | `dcc97ac6` |
| 4 | **float32 AMBIG cube + float64 reductions** | `_solve_ambig_logodds`, `_lse` | AMBIG solve 1103 → 657 ms (**~1.7×**); ≈½ the cube memory | `f_g` identical, RNA fractions ≤7.7e-7 (within `rtol=1e-6`) | `6a5792c4` |

### The KDE fix in detail (the load-bearing one)

`log_rho[node, grid] = log(M/E)[node] + log(f_g)[grid]` — every one of the ~90M query points lies on a
**bounded 1-D interval**, and `logpdf_kernel` is a **smooth 1-D function**. So:

- **Tile the query axis** in `_weighted_kde_logpdf` (in addition to the existing sample-chunk tiling) so
  the `(n_eval, n_samp)` product never materializes — bounds memory at *any* scale, exactly.
- In `_kde_logprior`, when the query set is large, **tabulate the exact kernel** on a dense lattice
  (`~16` points per KDE bandwidth) spanning the query range and **interpolate**. This preserves the
  real quadratic tails (the lattice spans the full range → no clamping, unlike the retired `logpdf`
  interpolation), turning `O(n_eval·n_train)` into `O(L·n_train + n_eval)`. Small query sets keep the
  exact per-point path, so golden scenarios stay bit-identical.

### The float32 pattern (the template for future precision work)

Store/evaluate the big `(m,K,Kt)` cube in **float32** (½ the memory, ~2× the SIMD throughput of the
elementwise `exp`/`log`), but **accumulate every reduction in float64** and keep the small `(m,K)`
marginals in float64. Concretely: `np.sum(..., dtype=np.float64)`, `_lse` accumulates its inner sum in
f64 and returns the input dtype, and the τ-marginal is lifted to f64 before the median/moment math.
Result: the posterior median (`f_g`) is *identical* and the cancellation-sensitive variances stay
accurate (≤6e-6). This mixed-precision recipe is reusable anywhere a large intermediate feeds a
reduction.

---

## The cost model — where time and memory go

Measured on the mo_3005 human library (971k fragments) at 0.6.2.

**Stage split (end-to-end ~105 s):**

| Stage | Time | Notes |
|---|---|---|
| index load | ~7 s | 457k transcripts + 5.8M splice-blacklist junctions |
| BAM scan (C++) | ~2 s | already fast |
| **calibration** | **~66 s (≈63%)** | the optimization target |
| per-locus EM (C++) | ~24 s | 9444 loci; OpenMP |
| write outputs | ~5 s | — |

**Inside calibration** (isolated cProfile, self-time; totals ~2× under cProfile):

| Function | Self-time | Nature |
|---|---|---|
| `simplex_logodds._solve_ambig_logodds` | ~34 s | the AMBIG (λ,τ) cube: float32 exp/log over `(m,K,Kt)` |
| `simplex._mixture_strand_loglik` | ~16 s | strand mixture over the cube (pass-invariant — see roadmap) |
| `bp_solver._scan` | ~15 s | **sequential per-node Python loop** (Gauss-Seidel; +22.9M `max`, 5.7M `log`/`exp` scalar calls) |
| `simplex_logodds._lse` | ~11 s | the cube's log-sum-exp reductions |

The AMBIG cube (`_solve_ambig_logodds` + `_mixture_strand_loglik` + most of `_lse`) is ~60% of
calibration; `_scan` is ~20%.

**Quant stage (~33 s), two Python hotspots worth noting** (downstream of calibration, *not* part of the
calibration-accuracy redesign, so durable targets):

- `priors._project_regions_to_loci` — ~14 s of pure-Python nested interval-overlap loops.
- `capture_eff_length._transcript_node_incidence` (`_add`) — a **1.87M-call** Python loop building the
  transcript→node incidence (annotation-only; the docstring notes it "could be precomputed at index build").

---

## Profiling infrastructure (how to reproduce)

Everything runs inside the activated `rigel` conda env.

**1. Full-pipeline per-stage + cProfile — `scripts/profiling/profiler.py`** (committed). Wraps the real
`run_pipeline` and reports per-stage wall/memory with optional `--cprofile`:

```bash
python scripts/profiling/profiler.py --bam sample.bam --index index/ --threads 8 --repeat 1 --cprofile
```

**2. Calibration-in-isolation (fast iteration loop).** The full pipeline is slow to iterate on; dump
`calibrate()`'s real inputs once, then re-run/profile calibration alone (deterministic, no scan/EM/write
overhead). This was built ad hoc this session — **promote to `scripts/profiling/` when continuing**:

```python
# dump once: monkeypatch calibrate() to pickle its (args, kwargs) then sys.exit, run run_pipeline.
# rerun: np.random.seed(0); calibrate(*args, **kwargs) under cProfile / timeit.
```

**3. Byte-identity + micro-timing for a single hot function** — the pattern that gated every win: build
fixed deterministic inputs for the target function, save a baseline `.npz`, and on each edit assert
`np.array_equal` (for byte-identical wins) or report `max|Δ|` (for tolerance wins) plus the min wall
time over N repeats. Used for `_solve_ambig_logodds` (an 8192-node cube).

**4. Accuracy guardrail — the 16-condition benchmark** (`calibration-benchmark` skill /
`scripts/sim/evaluate_suite.py`). Run **before shipping any perf change** to confirm the gDNA/nRNA/mRNA
net fragment-flow is unchanged. Force a re-quant after a calibration/quant change by deleting
`*/rigel_out/quant.feather`. This is the real-scenario complement to the goldens.

---

## Guiding principles

- **Prove output preservation, every time.** Byte-identical is the gold standard (verify with
  `np.array_equal` on the isolated function). Where a win requires it (float32), a small tolerance is
  acceptable *if* the goldens still pass at `rtol=1e-6` and the 16-condition benchmark is within its
  ~1% run-to-run noise. Never loosen the golden tolerance silently.
- **Gate approximations by size.** The tabulation and (potential) further float32 should engage only
  where they pay off (large query sets / genome scale) and leave small/golden inputs on the exact path,
  so the regression suite stays a meaningful bit-level check.
- **Measure the stage in isolation.** Genome-scale calibration cost is depth-independent; use the
  dumped-input harness, not full runs, to iterate.
- **Respect the accuracy redesign.** `simplex_logodds.py`/`simplex.py` (the solver mechanics) are safe
  to optimize and are *not* in the redescending-prior WIP; `bp_solver.py`/`calibrate.py`/
  `gdna_density_prior.py` (the prior) are — coordinate there.

---

## Roadmap — biggest opportunities

Ranked by impact × tractability. Effort/impact are rough; all should be verified against goldens +
benchmark.

### 1. Compile `_scan` (the sequential belief-propagation loop) — **biggest remaining calibration win**
- **Impact:** ~15 s (the whole `_scan` cost); calibration ~66 → ~50 s.
- **Why it needs compilation:** `_scan` is a genuine Gauss-Seidel sweep — each node reads the *running*
  belief updated by earlier nodes, so it cannot be vectorized. It's ~1.5M iterations of scalar Python
  (`math.log`, `max`, comparisons) × 4 directional passes.
- **Preferred route: C++ / nanobind.** rigel already ships 5 nanobind modules; a `_scan` kernel adds
  **no new runtime dependency** and fits the architecture. Port the per-node body (density message,
  structural emit gates, disagreement-aware precision combine) to C++; keep the Python orchestration.
  Verify against goldens (float ops may differ at ULP — within the authorized tolerance).
- **Alternative: numba** — faster to write (decorate the loop) but adds a heavy LLVM dependency to the
  conda recipe / wheels. Avoid for a distributed tool unless the C++ route proves too costly.
- **Cheaper partial win first:** hoist node-static quantities (eff-lengths, emit gates, floors) out of
  the loop and pre-filter nodes with no neighbour/signal from the scan order — trims per-iteration work
  and iteration count, byte-identically, though the Python interpreter overhead per node remains.

### 2. Quant-stage Python vectorizations — **durable, off the accuracy-redesign path**
- **Impact:** ~28 s total across the two hotspots below (quant ~33 → ~10 s).
- `priors._project_regions_to_loci` (~14 s): vectorize the region↔locus-block interval overlap +
  distribute per-array shares with `np.bincount`/`np.add.at` instead of the nested Python loops.
- `capture_eff_length._transcript_node_incidence` (~14 s, 1.87M calls): vectorize the per-ref
  `searchsorted` + build the ragged incidence with the standard `np.repeat`/offset trick; or (bigger)
  **precompute it at `rigel index` time** since it is annotation-only and sample-independent.
- **Caveat:** both feed `em_effective_length` / the per-locus prior → the goldens. Verify the
  vectorized output is set-equal (order-independent) or bit-equal to the loop before shipping.

### 3. Pass-invariant caching in the AMBIG solver — **bit-identical, compute-for-memory**
- **Impact:** ~7–10 s (post-float32) — the strand-mixture + spliced-floor `psi` is identical across the
  ~6 solve passes (init + sweep); only the message/prior terms change. Cache the invariant part and
  re-add the per-pass messages.
- **Cost:** holds the `(n_ambig, K, Kt)` cube resident across passes (~1.2 GB in float32). A speed-for-
  memory trade — reasonable on a 96 GB node, but weigh against the memory goal; consider a budget gate.

### 4. Wider mixed precision (the float32 template, item 4 above)
- **Impact:** moderate; the single-strand solver `_local_loglik_logodds` (~3 s) and some `node_sweep`
  arrays are still float64. Apply the same f32-cube/f64-reduction recipe. Diminishing returns after the
  AMBIG cube, but cheap and on-pattern.

### 5. Peak-memory attribution pass
- Current peak is 8.6 GB; the float32 cube is no longer the dominant term. **Profile the peak allocation
  site** (candidates: the loaded index incl. the 5.8M-junction blacklist; the `(n_nodes × K)` sweep
  arrays; the scoring CSR) before optimizing memory further — don't guess. `tracemalloc` or
  `memray` on the isolated calibration / quant stages will localize it.

### 6. Tunable knobs already in place
- `simplex_logodds._SOLVE_BLOCK_BYTES` (1 MB) tiles BOTH per-node solves — bit-identical at any value
  (every node solves independently and every reduction is within a row), so it is purely a cache/peak-memory
  knob. Measured at genome scale: flat from 1–16 MB, ~20 % slower unblocked, and sharply worse below
  256 KB where per-call overhead takes over. It replaced `_AMBIG_BATCH` (8192 rows), which bounded peak
  memory but sat 100× outside cache.
- The KDE lattice density (`_KDE_LATTICE_PTS_PER_BW=16`) trades interpolation accuracy for tabulation
  cost — already comfortably below any real regression.

---

## Ceiling & honest expectations

Output-preserving wins on the calibration have a ceiling: the AMBIG cube's `exp`/`log` FLOPs are
intrinsic (float32 already ~2×'d them), and `_scan` is Python-loop-bound. The realistic remaining
output-preserving headroom is roughly: `_scan` compile (~15 s) + quant vectorizations (~28 s, mostly
outside calibration) + AMBIG caching (~8 s). Beyond that, larger multipliers (coarser grids, fewer
sweep passes, a fundamentally cheaper AMBIG solve) **change the answer** and belong to the calibration
redesign, not here.

---

## Commit references (0.6.2 perf campaign)

- `2e93021a` — KDE-prior OOM fix (tiling + tabulation) + lean `logsumexp`
- `dcc97ac6` — in-place `psi` accumulation + CSE (byte-identical)
- `6a5792c4` — float32 AMBIG cube with float64 reductions
- Regression coverage: `tests/calibration/test_gdna_density_prior.py`
  (`test_weighted_kde_logpdf_query_tiling_is_exact`, `test_kde_logprior_tabulation_matches_direct_at_scale`)
