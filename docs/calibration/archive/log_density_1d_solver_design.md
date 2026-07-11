# Log-density 1-D calibration solver — design & implementation plan

**Status:** design, awaiting GO. Supersedes the 2-simplex lattice
(`simplex.py` / `simplex_sweep.py`) as the per-node solve. Motivated by two findings that
turn out to have **one fix**: the genome-scale memory wall and the log-space-noise objection.

---

## 0. TL;DR

The per-node solve currently grids the **linear** 2-simplex `(f₊, f₋, f_g)` with
`P = (K+1)(K+2)/2` points and evaluates an `(m, P)` array. At genome scale (`m ≈ 1.5M`
nodes, `P = 1891`) that array is **23 GB each, several live at once → >100 GB** (observed on
VCaP). It is also the *wrong space*: high-throughput count data is multiplicative
(PCR doubling → log-normal), so densities must be modeled in **log space** — which the
reliability machinery (`var~mean`, `ê(z)`) already does, while the solve does not. That
inconsistency is the source of the fragility.

**The fix solves both.** Re-express the per-node solve as a **1-D quadrature over the log gDNA
density `g = log ρ_g`**, with the strand tilt marginalized separately. This:

1. collapses `O(m·K²)` → `O(m·K)` (≈ 3 GB at genome scale, ~30× less memory *and* compute);
2. puts the priors, messages, global, and grid **all in log-density space**, consistent with
   `var~mean`/`ê(z)`;
3. gives the grid resolution where the dynamic range actually is (near `f_g → 0` depleted and
   `f_g → 1` enriched), which the uniform linear lattice does not.

It does **not** change the estimand (the *linear* fraction `f_g`), the count-zero-information
principle, the strand likelihood, the FB message structure, or the chain geometry.

---

## 0b. Phase 0 findings (measured 2026-06-28) — they refine the grid

Measured on the dissect cache (6 sim conditions) + the VCaP human index (752,654 regions). Tool:
`scripts/debug/phase0_logdens_census.py`. Three results, each load-bearing:

**(1) ~80–90% of nodes are NOT AMBIG → the lattice over-paid `O(K²)` for almost everything.**
Node-class shares are *structural* (signature-derived) and stable across gDNA level / strand /
capture: single-strand ~76–80%, AMBIG ~11–16%, G1 ~4–13%. Only the AMBIG minority needs the 2-D
(tilt) treatment; the rest is exact 1-D.

| | single-strand | AMBIG | G1 |
|---|---|---|---|
| sim (all conditions) | 76.1% | 11.4% | 12.5% |
| VCaP regions | 79.8% | **15.8%** | 4.4% |

**(2) `ρ_g` spans 5–6 decades under capture → log is mandatory; and the mass is BIMODAL at the
`f_g` vertices.** Oracle `ρ_g` dynamic range: **4.6–5.2 decades under capture** (p99/p1 ≈
30,000–40,000×), 1.9–2.2 decades off-capture. And the *true* `f_g` mass distribution piles at the
**vertices**, not the middle:

| condition | mass at `f_g < 0.01` | mass in `[0.05, 0.95]` | mass at `f_g > 0.99` |
|---|---|---|---|
| capture, ss99 (flagship) | ~0% | 66.6% | **26.1%** |
| off-cap, ss99 | ~0% | 14.9% | **77.3%** |
| off-cap, ss50 | ~0% | 15.1% | **80.9%** |
| zero-gDNA | **100%** | 0% | 0% |

700–775 regions per condition sit **within one linear-lattice cell (1/60 ≈ 0.0167) of a vertex** —
i.e. 26–81% of the mass is crammed into the top 1–2 lattice points, while the uniform lattice spends
most of its 1,891 points in the empty middle. **This refutes my doc's original "grid `log ρ_g`"
choice** (which resolves `f_g→0` finely but `f_g→1` only at spacing `~Δ`): the data demands fine
resolution at *both* vertices. See §6 for the fix.

**(3) Memory @ genome scale is ~3–5 GB, runnable.** VCaP `m ≈ 1.5M`:

| | per array |
|---|---|
| 2-simplex lattice `(m, P=1891)` | 22.8 GB → ×several → **>100 GB (the wall)** |
| log-odds 1-D `(m, K=60)` | **0.72 GB** |
| AMBIG tilt `(m_ambig, K, K_t=40)` | ≤ 4.6 GB (chunkable — the one term needing node-batching) |

---

## 1. Why we solve on a grid at all

The per-node latent is a composition `(f₊, f₋, f_g)` on the 2-simplex. The posterior is a product
of non-conjugate terms — an overdispersed Beta-Binomial strand likelihood, a one-sided clipped
spliced floor, a population prior, and several neighbour-message priors. There is no closed form,
so we evaluate `ψ` on a grid and take posterior summaries (the `f_g` median, the `f±` means, the
per-component variances = the precision state). The grid is the honest way to integrate a
non-Gaussian, boundary-constrained posterior. **We keep the grid.** We change *what we grid over*
and *in what space*.

---

## 2. The memory wall (the trigger)

`simplex_sweep._local_loglik` broadcasts `(m, 1)` node stats against the `(1, P)` lattice →
`(m, P)`. `_node_marginals` then forms `post = exp(ψ − logsumexp(ψ))`, another `(m, P)`. Each outer
iteration calls both **twice** (local solve `ψ_loc` + final solve `ψ`).

| quantity | value |
|---|---|
| `K = n_grid` | 60 |
| `P = (K+1)(K+2)/2` | 1,891 |
| human nodes `m` (regions+boundaries) | ~1.5M (752,654 regions) |
| one `(m, P)` float64 array | 1.5M × 1,891 × 8 B ≈ **23 GB** |
| arrays live per outer iter | several (ψ, post, reductions) → **>100 GB** |

The toy suites (~2.4K nodes) never approach this — they hide a hard genome-scale blocker.
Chunking the `(m, P)` arrays over node-batches would bound peak memory bit-identically, but it
keeps the `O(K²)` lattice and the wrong (linear) space. We do the structural fix instead.

---

## 3. Log space: reconciling densities and fractions

> *Your concern:* 40 counts = 32 RNA + 8 gDNA. Linear `f_g = 8/40 = 0.20`. In log2,
> `32 = 2⁵`, `8 = 2³`, "log fraction" `= 3/(3+5) = 0.375`. How does log space translate to our
> fractional solve?

### 3.1 The estimand is the linear fraction; the ratio-of-logs is not a quantity

The physical gDNA fraction **is** `8/40 = 0.20` — 20% of the molecules are gDNA. The
`3/(3+5) = 0.375` "ratio of logs" is not a meaningful quantity: it does not correspond to any
mixture proportion. Molecule pools **add in linear space** (a region's mass is `(ρ_g + ρ_rna)·L`),
not in log space — `log(a+b) ≠ log a + log b`. So the estimand we report stays `f_g = ρ_g / (ρ_g
+ ρ_₊ + ρ_₋)`, linear. Log space is **not** a redefinition of the fraction.

### 3.2 What log space *is* for, and the bridge

Log space is for the **noise model and the dynamic range**, not the mixing fraction. Each
*component density* `ρ_c` (gDNA, RNA₊, RNA₋ per bp) spans orders of magnitude — a depleted
intron at `ρ_g ≈ 0.01/bp` and a capture-enriched exon at `ρ_g ≈ 16/bp` differ 1000×. PCR doubling
makes the sampling noise **multiplicative** (log-normal), so the variance is roughly constant *in
log space*. Modeling `log ρ_c` (not `ρ_c`) is what every published RNA-seq model does (DESeq2/edgeR
log-fold-changes with NB), and it is the only way to handle the dynamic range without a single
heteroskedastic mess.

The bridge between the log-density we *model* and the fraction we *report* is one line. With the
observed total density `d = M / E` (the count-zero-info scale — see §4.2):

```
log ρ_c = log d + log f_c          ⟺          f_c = exp(log ρ_c − log d)
```

So a **Gaussian prior on `log ρ_c`** (mean from `ê(z)` or a message, precision from `var~mean`) is
**exactly a Gaussian prior on `log f_c`, shifted by the known constant `log d`**. We grid `log ρ_g`,
all priors are Gaussian-in-log, and we read out the linear `f_g` at the end. Your "solve for
log(density)" instinct is precisely right, and the fraction falls out by subtracting `log d`.

### 3.3 The current solver is internally inconsistent — that is the fragility

This is the deeper point your objection exposes. Today:

- `variance_model.MonotoneVarMean` fits `σ²_bio(μ)` in **log-log** space (`lx, ly = log x, log y`).
- `fit_enrichment_transfer` / `ê(z)` fits in **log-log** space.
- But `simplex_sweep._local_loglik` applies the messages as **linear-fraction** Gaussians
  `−½·prec·(f_c − mode)²` (line 105/113) and the global as a **linear-fraction** Beta pseudo-count
  (`_binom_pseudo`), on a **uniform linear** lattice.

We fit reliability in log space and then *spend* it in linear space. A linear-fraction Gaussian
cannot represent "within a multiplicative factor of 2 of `ρ_g = 0.01`" and "within 2× of `ρ_g = 16`"
with the same machinery — but a log-density Gaussian does so natively (constant log-variance =
constant multiplicative factor). The redesign removes the mismatch: **fit in log, solve in log.**

---

## 4. The robust solver: separate the tilt from the magnitude

The 2-simplex has two degrees of freedom. They are physically different and should be modeled in
different spaces. The current code already proves they are orthogonal; we make that structural.

### 4.1 The strand likelihood constrains only the tilt (a proportion)

From `simplex._mixture_strand_loglik` (line 69), the mixture plus-strand probability is an
**identity**:

```
p = ½·f_g + κ·f₊ + (1−κ)·f₋ = ½ + (κ−½)·(f₊ − f₋) = ½ + (κ−½)·t,   t ≜ f₊ − f₋
```

The strand evidence depends **only on the tilt `t`**. Strand is a *proportion* model (which way do
the reads lean), naturally linear/logit in `t` — it does **not** go to log-density. Its precision is
the intrinsic BB Fisher information `N·(2κ−1)²` (count-zero-info: the count enters here and nowhere
else as composition information). *(Note: the overdispersion variance term, line 72–77, depends on
the individual `f_c` magnitudes, so the strand `ψ` is fully determined once both `g` and `t` are
fixed — see §5.2. The mean still depends on `t` alone.)*

### 4.2 The magnitude priors live on log-density; the grid axis is the log-odds `λ`

The gDNA-vs-RNA **magnitude** is where the priors live, and where the dynamic range is:

- the population/enrichment prior `ê(z)` / `ρ_global` → a Gaussian on `log ρ_g`;
- the gDNA neighbour message → a Gaussian on `log ρ_g`;
- the per-strand RNA messages → Gaussians on `log ρ_₊`, `log ρ_₋`;
- the spliced floor → a lower bound `log ρ_s ≥ log(spliced_s / E_spl)`.

The observed total `d = M/E` is fixed (count-zero-info: the count `M` sets the strand precision and
the scale `d`, but carries **no** intrinsic gDNA/RNA split). Given `d`, the only free magnitude dof
is the **split**, and (Phase 0 §6) the grid coordinate for it is the **log-odds**
`λ = log ρ_g − log ρ_rna = logit(f_g)`. The priors above stay on `log ρ_g` / `log ρ_±` and are
evaluated *pointwise* at each grid `λ` via `log ρ_g(λ) = log d + log σ(λ)`,
`log ρ_rna(λ) = log d + log σ(−λ)`.

### 4.3 The model

Latent: `(λ, τ)` where `λ = logit(f_g) ∈ (−∞, ∞)` is the gDNA-vs-RNA log-odds (grid: `[−L, L]`) and
`τ ∈ [−1, 1]` is the RNA-internal tilt. Derived:

```
f_g = σ(λ),   ρ_g = f_g·d,   ρ_rna = (1−f_g)·d,   ρ_₊ = ρ_rna·(1+τ)/2,   ρ_₋ = ρ_rna·(1−τ)/2
```

(`τ` is the RNA-internal tilt; the strand-likelihood tilt `t = f₊ − f₋ = (1−f_g)·τ`. The
strand-likelihood tilt and the RNA-message split are the same dof.) The posterior is

```
log P(λ, τ) =  strand_LL(λ, τ)                       # tilt (+ magnitude via overdispersion)
             + spliced_floor(λ, τ)                    # one-sided lower bounds on ρ_±
             + N(log ρ_g(λ); μ_msg_g, σ²_msg_g)       # gDNA message, Gaussian-in-log
             + N(log ρ_₊(λ,τ); μ_msg_₊, σ²_msg_₊)      # RNA₊ message
             + N(log ρ_₋(λ,τ); μ_msg_₋, σ²_msg_₋)      # RNA₋ message
```

This is the *same* hierarchy as today (`CALIBRATION_ARCHITECTURE.md` §1) — strand likelihood +
node-class prior + imputation — re-expressed with the magnitude priors in log-density and the grid
over the log-odds `λ`.

---

## 5. The 1-D reduction

### 5.1 Single-strand & intergenic nodes — exact 1-D over `λ`

These are the **majority** of nodes (Phase 0: ~80–90%). For them the tilt is **determined**, so the
solve is genuinely 1-D — no approximation:

- **single-strand** (`free_pos ^ free_neg`): `ρ_rna` is entirely on the live strand, `τ = ±1`. The
  strand LL and RNA message collapse to that one strand. Grid `λ` only.
- **intergenic / G1 sink**: `ρ_rna = 0`, `f_g = 1` (the locked `{0,0,1}` vertex). Trivial.

We evaluate `ψ(λ)` on a `K`-point `λ`-grid → `(m, K)`, take the posterior over `λ`, and read out
`f_g = σ(λ)`, `var(λ)`. `f₊`/`f₋` follow from `τ = ±1`.

### 5.2 AMBIG nodes — marginalize the tilt (cheap inner 1-D)

AMBIG nodes (`free_pos & free_neg`, both strands live — the hard capture-exon case; Phase 0:
~10–16%) carry the extra tilt dof. For each `λ` we marginalize `τ` over a small inner grid `K_t`:

```
ψ(λ) = logsumexp_τ [ strand_LL(λ, τ) + RNA_msg(λ, τ) ]   over τ ∈ [−1, 1], K_t points
```

The strand LL varies with `τ` (mean via `p`, variance via the per-strand overdispersion), so this
is a genuine integral — but it is **1-D and cheap**, and only AMBIG nodes pay it. Memory is
`(m_ambig, K, K_t)` ≤ 4.6 GB at genome scale (Phase 0 §0b.3); node-chunked or done semi-analytically
(the RNA-message part is Gaussian in `log ρ_±`; only the strand LL is non-Gaussian in `τ`). **The
entire reason the 2-simplex lattice materializes `(K+1)(K+2)/2` points for every node is to carry
this AMBIG tilt — which only the AMBIG minority needs, and only as a 1-D inner integral.**

### 5.3 Readout

The posterior over `λ` gives the precision state in the natural space:

- `f_g` = posterior **median** of `f_g = σ(λ)` (robust, as today — avoids the overdispersion skew
  that biased the grid-MAP);
- `f₊`, `f₋` = posterior means (the current-state partition for the RNA messages);
- the **precision state** is `var(λ)` — a *log-odds* uncertainty, the natural and dimensionless
  confidence that the messages should carry. This is strictly better than `var(f_g)` (which is
  scale-dependent and collapses near the vertices regardless of real confidence — exactly where
  Phase 0 shows most of the mass lives).

---

## 6. Grid construction — the robustness core (revised by Phase 0)

This is where "the scale changes a lot" is handled, and it is the single biggest robustness gain.

**Today:** a *fixed, uniform, linear* lattice on `[0, 1]`, spacing `1/60 ≈ 0.017`. Phase 0 (§0b.2)
shows this is not merely wasteful but actively wrong: 26–81% of the true mass sits within one cell of
a vertex, while the lattice spends ~1,800 of its 1,891 points in the empty middle.

**The axis: the gDNA-vs-RNA log-odds.** My original "grid `log ρ_g`" was asymmetric (resolves
`f_g→0` but not `f_g→1`); the bimodal-at-vertices data needs *both*. The fix is to grid

```
λ ≜ logit(f_g) = log( ρ_g / ρ_rna ) = log ρ_g − log ρ_rna
```

the **log-odds of gDNA vs RNA** — still "log-density" (a *difference* of log-densities, honoring the
multiplicative-noise model), and the *natural* axis: since `d = M/E` is observed (count-zero-info),
the only free magnitude dof is the **split**, which is exactly `λ`. Uniform spacing in `λ` is
geometric in `f_g` at *both* ends → fine resolution at `f_g→0` (depleted introns, zero-gDNA) **and**
`f_g→1` (captured pure-gDNA exons). `f_g = σ(λ)` at readout.

**Bonus: the log-odds bounds the range, so the grid can be FIXED — no node-adaptivity needed.** The
5–6 decade `ρ_g` range (which *would* force a node-adaptive `log ρ_g` grid) maps into a bounded `λ`
window: `f_g ∈ [5e-5, 1−5e-5]` is `λ ∈ [−10, 10]`. A fixed `λ ∈ [−L, L]` grid at `K ≈ 60`
(spacing ≈ `2L/K`) covers every regime; the **node-specific prior is evaluated pointwise on the
fixed grid**, not by moving the grid. This is strictly simpler than the node-adaptive `log ρ_g` plan.

**Priors on the `λ` grid.** Each prior stays Gaussian-in-log-density, evaluated at the grid point:

```
log ρ_g(λ)   = log d + log σ(λ),     log ρ_rna(λ) = log d + log σ(−λ)
```

so the `ê(z)`/global/gDNA-message Gaussians (on `log ρ_g`) and the RNA-message Gaussians (on
`log ρ_±`) are smooth functions of `λ`, evaluated once per grid point. The **depleted-floor +
pseudocount** (validated, 42× below `ρ_global` under capture) now anchors the *prior mean* `μ_λ` (so
zero-gDNA pins `λ → −∞` / `f_g → 0` gently) rather than serving as a grid bound — a cleaner role.

`K` and `L` are the accuracy↔cost knobs, set by verification (sweep `K`, check convergence) like the
current `n_grid = 60` — not magic constants.

---

## 7. Memory + compute

| | 2-simplex lattice (now) | log-density 1-D (proposed) |
|---|---|---|
| per-node array | `(m, P)`, `P = 1891` | `(m, K)`, `K ≈ 60` + AMBIG `(m_ambig, K, K_t)` |
| one array @ 1.5M nodes (measured) | **22.8 GB** | **0.72 GB** (+ ≤ 4.6 GB AMBIG, chunkable) |
| peak (several live) | **>100 GB** | **~3–5 GB** |
| compute / node | `O(P) = O(K²)` | `O(K)` single-strand (~84%), `O(K·K_t)` AMBIG (~16%) |
| genome-scale | **OOM** | **runnable** |

~30× memory and ~30× compute at genome scale, plus the FB scans are unchanged (already O(m)).

---

## 8. Implementation plan (phased; each phase verified before the next)

The change is behind the `node_sweep` / `_solve_nodes` interface, so the chain, geometry, statics,
init, and the FB driver are untouched in structure — we swap the per-node kernel and the message
*space*.

- **Phase 0 — measure & ground (no code change). ✅ DONE (§0b).** Findings: ~80–90% non-AMBIG
  (1-D-exact); `ρ_g` spans 5–6 decades under capture; true `f_g` mass is bimodal at the vertices
  (26–81% at `f_g>0.99`, 100% at `f_g<0.01` for zero-gDNA) → **grid the log-odds `λ`, fixed `[−L,L]`,
  no node-adaptivity** (§6); genome-scale memory ~3–5 GB. Tool: `scripts/debug/phase0_logdens_census.py`.

- **Phase 1 — log-density 1-D core (single-strand + intergenic). ✅ DONE & VERIFIED.**
  `simplex_logodds.py` (`_solve_nodes_logodds`): the `λ`-grid solve evaluating the SAME ψ terms as the
  lattice (strand / spliced / Jeffreys / global — still linear-fraction; the *space* migration is
  deferred to Phase 3) + the change-of-variable Jacobian `log(f_g(1−f_g))`. Verified via
  `scripts/debug/phase1_logodds_equiv.py` (a capture-gated probe in `node_sweep`, identical inputs):
  **per-node `|Δf_g|` → 0 as O(1/K)** (median 0.028→0.014→0.006 at K=60→120→240, flagship) — proving it
  integrates the *same* posterior (a measure/Jacobian bug would converge to a different value). At
  K=60 the aggregate differs ≲2% in most regimes; the one large case (zero-gDNA +43%) is the log-odds
  being the **better-converged quadrature** — the lattice climbs toward it as K↑ (71K→99K vs
  102K→108K), i.e. the coarse linear grid was *accidentally under-resolving* the model's phantom (same
  pattern as the FB-vs-Jacobi cliff). 188/188 calibration tests green (production path untouched).
  *Finding:* the log-odds grid is fine at the vertices (where the true mass is — Phase 0) and coarse in
  the middle (where the κ=½ phantom lands) — the inverse of the lattice; `K`/`L` tuning + the Phase-3
  log-density priors address the middle. **The faithful solver UNMASKS the zero-gDNA/κ=½ phantom the
  lattice's coarse grid hid — a model-level issue (global/Jeffreys at κ=½) for Phase 3 / the E4 work,
  not a solver bug.**

- **Phase 2 — AMBIG tilt marginalization. ✅ DONE & VERIFIED.** `_solve_ambig_logodds`: the 2-D
  `(λ, τ)` cube, SAME ψ terms (no Jeffreys — matches the lattice at AMBIG), `logsumexp` over τ, then
  the **τ-independent 2-D Jacobian** `log|J| = log f_g + 2·log(1−f_g)` (uniform `(f_pos,f_neg)` →
  `(λ,τ)`). Verified (`phase1_logodds_equiv.py`, AMBIG probe): **O(1/K) convergence** (max |Δf_g| flagship
  0.050→0.026, unstranded 0.042→0.022 at K=60→120; aggregate gDNA mass +1.8%→+0.7% flagship,
  +0.58%→+0.02% unstranded). At K=60 the AMBIG agreement is **tighter than single-strand** (max 0.050
  vs 0.133); worst cases at f_g≈0.5–0.75 (the log-odds middle). 188/188 tests green. The cube is
  materialized only for the AMBIG subset (the caller masks).

- **Phase 3 #1 — mechanical dispatch into `node_sweep`. ✅ DONE.** Config flag
  `CalibrationConfig.calibration_solver = "lattice" (default) | "logodds"` (+ `sweep_logodds_window`,
  `sweep_n_tilt`). `node_sweep` dispatches the phase-A local + phase-D final solves through a backend
  closure (`_solve_nodes_logodds_all` for logodds, the lattice otherwise); the FB α/β scan + the
  fraction-space messages + the global construction (on the solver's `f_g` grid) are SHARED — only the
  grid changes. **Lattice default byte-identical (188 tests green).** Logodds end-to-end vs lattice on
  the per-region gDNA mass (the prior input): **±2% on all gDNA-present conditions** (flagship +0.3%,
  ρ_g within ~1%), **zero-gDNA +73%** — the phantom unmasking, amplified by FB propagation (a small
  absolute on a zero-truth library; the model-level issue #2 must fix, not a dispatch bug). Ruff clean.

- **Phase 3 #2 — messages + global → log-density.** *(needs design alignment — intersects
  `message_state_separation`.)* Convert `_scan` to send a **log-density Gaussian** (`mode = log
  ρ_src_imputed`, `prec` from `var(λ)` + log-space `σ²_bio`) instead of the linear-fraction Gaussian;
  the global from `_binom_pseudo` (Beta-on-`f_g`) to a Gaussian-on-`log ρ_g`. This is the actual
  robustness payoff (fixing the §3.3 fit-in-log/solve-in-linear inconsistency) and where the zero-gDNA
  phantom gets a real shot at improvement (sharper log-density global near `f_g→0`). **Verify** the
  full 16-condition net-leak benchmark vs the shipped baseline (`@0ab636e1`) + the `(K, L)` sweep vs
  ground truth; goldens regenerated intentionally.

- **Phase 4 — genome scale.** VCaP on one chromosome → confirm ~3 GB; then the full genome end-to-end
  → confirm runnable + the **`(K, L)` sweep on real data** (now feasible) + extract the real `ê(z)`
  points (the original blocked task).

The `calibration_solver` flag keeps the lattice available for A/B through Phase 3, retired after Phase 4.

---

## 9. Verification

- **Bit-comparable, not bit-identical**: the log grid and space differ from the linear lattice, so
  the gate is a numerical tolerance on `f_g` per node + the aggregate net-leak, not byte-identity.
  We assert the *aggregate* metrics (net fragment flow per condition) move less than the
  pass-to-pass solver noise, and spot-check the known hard nodes (AMBIG capture exons,
  zero-gDNA introns).
- **Memory bound**: measured RSS on the VCaP chromosome slice < a few GB; the full genome completes.
- **No-regression on the suites**: the calibration unit tests + golden outputs (regenerated with a
  documented rationale, as for prior intentional solve changes).
- **The hard regimes** from the gateless work (`calib_honest_precision_gateless`) — off-cap
  unstranded ss50 basin chaos, the E4 zero-gDNA+nascent phantom — are *re-measured* under the new
  solver; the log grid's better extreme-resolution may help the basin stability, but this is a
  measurement, not an assumption.

---

## 10. Risks & open questions

1. **AMBIG marginalization accuracy.** The inner `t`-grid must match the 2-simplex's handling of the
   coupled (tilt × magnitude) overdispersion. Mitigation: Phase 2 verifies exactly here; `K_t` is a
   verified knob. The semi-analytic alternative (Gaussian RNA × tabulated strand-LL-in-`t`) is a
   fallback.
2. ~~**High-end (`f_g → 1`) resolution.**~~ **RESOLVED by Phase 0.** Measured: 26–81% of true mass
   sits at `f_g > 0.99` (captured/off-cap), 100% at `f_g < 0.01` (zero-gDNA) — bimodal at *both*
   vertices. The asymmetric `log ρ_g` grid is out; the **log-odds `λ`** grid (§6) resolves both ends
   geometrically and, because `λ` is bounded for `f_g ∈ [5e-5, 1−5e-5]`, a *fixed* `[−L, L]` grid
   suffices (node-adaptivity dropped). `L` is verified, not magic.
3. **Grid size `K` and window `L` are data-driven knobs, NOT the inherited magic `60`.** `sweep_n_grid
   = 60` was never calibrated; Phase 1 shows it is under-converged (the lattice under-counts its own
   K=240 answer by ~30% in zero-gDNA). The log-odds `O(m·K)` solver is what makes K-exploration on
   real data feasible at all (the lattice's `O(m·K²)` needs ~370 GB at K=240 on VCaP; the 1-D is
   ~3 GB). **Selection: on sim — sweep `(K, L)` to the accuracy-vs-oracle plateau; on real VCaP — sweep
   to the self-consistency plateau + the cost curve (a Phase-4 deliverable, now possible).** `L` sets
   how far toward the `f_g→0`/`f_g→1` vertices we resolve (Phase 0: the mass is at the vertices, so `L`
   matters as much as `K`). The inner τ-grid `K_t` can be modest and *shrinks* as `f_g→1` (the
   2-simplex resolves the tilt 61/31/7 at f_g=0/0.5/0.9 — fine when there is RNA to split, coarse when
   `ρ_rna→0`). Neither is a fixed constant; both are verified against data.
4. **Message-space change is load-bearing.** Sending log-density Gaussians changes how disagreement
   propagates (a multiplicative, not additive, notion of "close"). This is the *intended* robustness
   gain, but it interacts with the FB relay; Phase 3's benchmark is the arbiter.
5. **The strand overdispersion's magnitude dependence** (line 72–77) means `t` and `g` are not fully
   separable in the variance term — handled by evaluating on the joint `(g, t)` for AMBIG, but it
   means the AMBIG marginal is not purely "strand × independent magnitude". Documented, verified in
   Phase 2.

---

## 11. What does **not** change

- The estimand: the linear fraction `f_g` (and `f₊`, `f₋`).
- The count-zero-information principle: the count sets the strand precision + the scale `d`, nothing
  else.
- The strand Beta-Binomial likelihood and its `N·(2κ−1)²` Fisher information.
- The bipartite region↔boundary chain, the FB forward/backward structure, init, the geometry/statics,
  the spliced floor, the enrichment transfer `ê(z)`, the depleted-floor anchor, the unstranded+nascent
  caveat.
- The `var~mean` reliability machine (already log-space) — it becomes *consistent* with the solve,
  not replaced.
```
