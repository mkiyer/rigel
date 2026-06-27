<!-- title: Simplex-lattice grid — cost, accuracy, and solver alternatives (deferred exploration) -->
# Simplex-lattice grid: cost, accuracy, and better-solver ideas

**Status:** ideas log (2026-06-27), DEFERRED for later exploration. Records what was *measured* about the
per-node 2-simplex lattice solve (`config.sweep_n_grid`, default 60) and the candidate alternatives, so we can
pursue them deliberately. Prompted by the question: *the grid size of 60 is arbitrary — would a larger grid
help, what does it cost, and what would a better solver look like?* The short answer (measured): a larger grid
does **not** improve flagship accuracy, costs **O(K²)** compute + memory, and is a genome-scale bottleneck. The
real wins are in the **readout** and a **semi-analytic solve**, not the grid resolution.

## 0. What the grid is

The per-node deconvolution (`simplex_sweep._local_loglik` + `_node_marginals`, driven by `bp_solver.node_sweep`)
evaluates the node's log-posterior over a discretized 2-simplex of `(f_rna+, f_rna-, f_gdna)`:
`_simplex_lattice(K)` returns `P = (K+1)(K+2)/2` lattice points. Every node's solve is a `(m_nodes × P)` array;
the readout marginalises it (`f_g` posterior **median**, `f_±` posterior means, per-component variances).

## 1. Measured cost + accuracy (flagship `gdna_gdna300_ss_0.99_nrna_none_capture_on`, 2443 nodes)

`scripts/debug/` n_grid sweep (single-thread, `tracemalloc` peak; |err| = Σ|mass_gdna_contained − oracle|):

| K | lattice pts | time | peak mem | net err | **\|err\| total** | \|err\| AMBIG | \|err\| single-strand |
|---|---|---|---|---|---|---|---|
| 40 | 861 | 7.6 s | 154 MB | −20,138 | 45,359 | 11,808 | 32,356 |
| **60** | 1,891 | 9.4 s | 335 MB | −16,639 | 45,910 | 10,901 | 33,557 |
| 90 | 4,186 | 11.8 s | 739 MB | −18,418 | 45,952 | 10,131 | 34,199 |
| 120 | 7,381 | 14.0 s | 1.3 GB | −18,977 | 48,395 | 8,924 | 37,803 |
| 180 | 16,471 | 22.1 s | 2.9 GB | −21,052 | 45,983 | 10,639 | 33,597 |
| 240 | 29,161 | 40.6 s | 5.1 GB | −20,453 | 48,290 | 9,036 | 37,449 |

**Findings:**
- **Cost is O(K²).** `P = (K+1)(K+2)/2`; the hot array and the memory both scale with `P`. At K=240 the solve
  peaks at **5.1 GB on a 2443-node toy** — at genome scale (millions of nodes) this is a hard wall.
- **Larger grid does NOT improve accuracy.** Total |err| is flat-to-worse (45K→48K, no trend). AMBIG |err|
  improves only modestly (~25% to K=120) then is noisy; the dominant single-strand |err| does not improve at
  all. The grid controls only **quantization** error; the flagship error is **continuous-model** error (the
  boundary-seam message under-stating exon-interior density, ê saturation, the strand tilt). A finer grid just
  resolves the same *biased* posterior more precisely.
- **The coarse grid's median-snap is a mild accidental regulariser** — K=60 gives the *best* net err (−16,639)
  of the whole sweep (cf. the "Jacobi under-convergence was accidental regularisation" lesson in
  `forward_backward_state.md`). So **60 is, in effect, a reasonable sweet spot**, not just an arbitrary constant.

## 2. The −1/60 quantization snap (the original trigger)

The cluster-A AMBIG nodes show a uniform `final − local = −0.0167 = exactly 1/60` — the posterior **median**
readout snaps the continuous −0.014…−0.024 propagation shift to one grid cell. Two facts bound its importance:
1. The continuous shift is real (the seam message), so the snap is a small *quantization* on top of a real
   drag — fixing the snap recovers at most ~one cell per node, sub-dominant to the continuous error.
2. The snap is the **median** readout's doing; the **mean** would not snap.

## 3. Better-solver ideas (each needs its own exploration)

### 3a. Posterior MEAN readout (with a skew correction) — zero extra grid cost
`_node_marginals` **already computes** the posterior mean `mg = post @ fgg` (it is used for `var_g`); we just
return the **median** (`fg_med`). Switching `f_g` to the mean would erase the −1/60 snap at **no extra cost**.
**Caveat (load-bearing):** the median was chosen deliberately — the docstring records that the mean "biased the
grid-MAP low (the +8.7 pt regression)" because of the strand-posterior overdispersion *skew*. So a naive mean
is a regression. The exploration: a **skew-corrected mean** (e.g. mean with a bias term estimated from the
posterior's third moment, or a trimmed/winsorised mean) that keeps the median's robustness while removing the
quantization snap. Cheap to prototype (the moments are already computed).

### 3b. Semi-analytic / reduced-dimension solve — the real O(K²)→O(K) (or O(1)) win
The per-node posterior is **mostly Gaussian**: the cross-node messages (gDNA + per-strand RNA) and the global
prior are all Gaussians in fraction space; only the **strand Beta-Binomial likelihood** is non-Gaussian, and it
constrains a *1-D* quantity (the per-strand split). Ideas:
- **1-D quadrature over `f_g` only.** Marginalise `f±` analytically (they are Gaussian given `f_g` and the
  strand split), reducing the lattice from a 2-simplex `O(K²)` to a 1-D grid `O(K)` in `f_g`.
- **Moment-matched / Laplace solve.** Approximate the strand BB by its Gaussian moments and combine all terms in
  closed form (an information-form Gaussian product), with a correction only where the BB is far from Gaussian
  (thin counts / near the simplex walls). Potentially **O(1) per node**.
- Either removes the genome-scale memory wall entirely and makes the solve trivially vectorisable/parallelisable.
- Risk: the BB skew at AMBIG/thin-count nodes is exactly where the Gaussian approx is worst — needs a guarded
  fallback (the same nodes the lattice serves today). Validate against the lattice solve for agreement.

### 3c. Adaptive lattice refinement — likely low payoff
Refine the grid only near the posterior mode (coarse global grid + local zoom). Cuts cost for a target
accuracy, but since accuracy is **continuous-model-bound** (§1), the accuracy payoff is small; mainly a
perf play, and dominated by 3b. Lower priority.

### 3d. Keep K=60 (or even drop to ~40–50) + invest in the readout/model
Given the measured flat accuracy vs K, the pragmatic near-term stance is: **do not raise K** (it only costs
memory). If perf matters, K≈40–50 is nearly as accurate and ~2× cheaper in memory; spend the effort on 3a/3b and
on the continuous-model errors (the seam-message density, ê saturation) that actually move |err|.

## 4. Performance / genome-scale notes

- The lattice solve is **already vectorised** over `(m, P)`; raising K does not parallelise error away — it just
  enlarges the array (memory bandwidth + footprint bound).
- Memory `O(m·K²)` is the binding constraint at genome scale. Chunking over nodes bounds peak memory but not
  total work; **3b (reduced-dimension solve) is the structural fix.**
- The outer var~mean fixed point and the two batched lattice solves per FB pass mean the lattice solve is the
  hot path; 3b would move the bottleneck off it.

## 5. Decision

Grid resolution is a **measured dead-end for accuracy** and a **perf liability** — not a lever for the flagship
error. Deferred ideas, in priority order when revisited: **3b (semi-analytic) > 3a (skew-corrected mean) > 3d
(keep/lower K) > 3c (adaptive)**. The flagship accuracy lever is elsewhere — the global-vs-strand authority on
single-strand exons (the 3× larger source) and the boundary-seam message density — tracked separately.
