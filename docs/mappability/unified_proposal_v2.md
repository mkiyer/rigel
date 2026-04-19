# Unified gDNA Calibration — Proposal v2

_Date: 2026-04-18 · Supersedes `unified_proposal.md`_

## Why a v2

v1 proposed a two-component Poisson mixture over regions plus an
unmapped-read pathway. After review, three concerns drove a redesign:

1. **The mixture is over-engineered.** It introduces a soft RNA-rate
   nuisance that is fragile under the 5+ orders of magnitude dynamic
   range of real RNA expression, and requires EM bookkeeping for what
   is fundamentally a one-dimensional estimation problem.
2. **The unmapped pathway is research, not infrastructure.** It is
   STAR-specific and exposes a π_RNA bias that we don't yet know how
   to bound on real data. It should be probed offline first and only
   promoted to a likelihood term once we understand the bias surface.
3. **Tile-based denominators are a nonstarter.** The region partition
   (transcript-aware, unequal lengths) is the natural quantum of the
   index and must remain so. Any new estimator must consume regions
   as-is.

This document proposes a **single, simple, region-based estimator**
that subsumes the role of P₁₀ without introducing soft mixtures, magic
percentiles, or a tile lattice. It also retires the heuristic
strand/density blend in favour of inverse-variance pooling.

---

## 1. Statement of the problem

Given the existing region partition (one row per genomic region with
columns `length`, transcript memberships, and unspliced fragment
counts), estimate a single scalar:

$$\lambda_G = \text{per-bp gDNA fragment rate}$$

under the working assumption that gDNA fragmentation is approximately
uniform across the mappable genome, while RNA is concentrated in a
small fraction of regions with arbitrary intensity.

We must:
- Treat zero-count regions as informative (large exposure × zero count
  is the strongest possible upper-bound evidence on λ_G).
- Survive the RNA dynamic range (a single highly-expressed exon must
  not dominate the estimate).
- Use a continuous mappability weight, not a hard mask.
- Combine with the strand pathway without magic blending constants.
- Stay simple enough that the entire estimator is one short function.

---

## 2. Core idea — the profile-Poisson estimator

### 2.1. Model

For region $i$ with mappable exposure $E_i = f_i \cdot S_i$ (mappable
fraction × genomic span), let $k_i$ be the unspliced fragment count.
Model:

$$k_i \sim \text{Poisson}(\lambda_G\, E_i + r_i),\qquad r_i \ge 0$$

where $r_i$ is the **non-negative RNA contribution** for region $i$.
The $r_i$ are nuisance parameters, one per region, with no shared
prior — this lets each region's RNA rate be arbitrary, sidestepping
the dynamic-range problem entirely.

### 2.2. Profile likelihood

For fixed $\lambda$, the MLE of each $r_i$ is closed form:

$$\hat r_i(\lambda) = \max(0,\; k_i - \lambda E_i)$$

Substitute back. Define the **null-dominated set**

$$\mathcal A(\lambda) = \{\,i : k_i < \lambda E_i\,\}$$

(regions where the observed count is below what the gDNA rate alone
would predict; any RNA contribution is set to zero). The profile
log-likelihood collapses to

$$\ell(\lambda) = \sum_{i \in \mathcal A(\lambda)} \big[\, k_i \log(\lambda E_i) - \lambda E_i\,\big] + \text{const}$$

Setting $d\ell/d\lambda = 0$ on a bracket where $\mathcal A(\lambda)$
is locally constant gives the closed-form M-step:

$$\boxed{\;\hat\lambda = \frac{\sum_{i \in \mathcal A(\hat\lambda)} k_i}{\sum_{i \in \mathcal A(\hat\lambda)} E_i}\;}$$

This is a **self-consistent rate over the regions whose counts fall
below what that very rate predicts**.

### 2.3. The algorithm

```
λ ← initial guess  (e.g. total counts / total exposure)
repeat:
    A ← { i : k_i < λ · E_i }
    λ_new ← Σ_{A} k_i / Σ_{A} E_i
    if |λ_new − λ| < tol: break
    λ ← λ_new
```

Convergence is monotone (the sequence is non-increasing because each
iteration narrows or holds $\mathcal A$, both of which can only lower
$\hat\lambda$) and reaches a fixed point in a small number of steps —
typically ≤ 10. The whole estimator is ~15 lines of Python.

### 2.4. Properties (and why they matter)

| Property | Mechanism |
|---|---|
| **Zeros are fully informative.** | Any region with $k_i = 0$ satisfies $0 < \lambda E_i$ for all $\lambda > 0$, so it is always in $\mathcal A$. Its exposure $E_i$ contributes to the denominator, pulling $\hat\lambda$ down. |
| **RNA dynamic range is non-issue.** | A region with arbitrarily high $k_i$ leaves $\mathcal A$ as soon as $k_i \ge \lambda E_i$. It contributes nothing to the rate estimate, regardless of how extreme its count. No outlier control, no winsorisation, no log-transform needed. |
| **Continuous mappability weight is native.** | $E_i = f_i S_i$ enters as exposure. A region with $f_i = 0$ contributes zero exposure and is silently irrelevant; a region with $f_i = 0.4$ contributes 40% of its span. No hard mask, no containment filter. |
| **No arbitrary cutoff.** | The set $\mathcal A$ is determined by the data and the estimate — there is no "10th percentile" knob. The only hyperparameters are convergence tolerance and a minimum-regions-in-$\mathcal A$ guard. |
| **Closed-form Fisher information.** | At the fixed point, $I(\hat\lambda) = \sum_{i \in \mathcal A(\hat\lambda)} E_i / \hat\lambda$. Gives a direct standard error and enables principled pooling with the strand pathway (§4). |
| **Diagnoseable.** | $|\mathcal A| / N$ and $\sum_{i \in \mathcal A} E_i / \sum_i E_i$ are interpretable summaries (fraction of regions / fraction of mappable genome that "look like pure gDNA"). |

### 2.5. Relation to the v1 proposal

This is the hard-EM (Viterbi) limit of v1's two-component Poisson
mixture. The hard limit is preferable because:
- The M-step is closed form, no soft posteriors to track.
- The "RNA component rate" $\mu_R$ is replaced by per-region $r_i$ — no
  pooling of unrelated genes' expression rates.
- Convergence proof is trivial; numerical instability impossible.
- Implementation is a one-page function rather than an EM module.

We lose the soft posterior $\pi_i$ (probability that region $i$
contains RNA), but we gain a binary indicator membership in $\mathcal
A(\hat\lambda)$ that serves the same diagnostic purpose.

### 2.6. Relation to the historical "two-state EM"

The previous Rigel architecture used a two-state EM (gDNA-only vs.
gDNA+RNA) seeded by the P₁₀ percentile. The profile-Poisson estimator
is a **strict superset of this idea** with two simplifications:

1. The latent state is replaced by an analytic indicator
   ($k_i < \lambda E_i$), eliminating the EM bookkeeping.
2. No seed is required: any positive starting $\lambda$ converges to
   the same fixed point because the iteration map is monotone and the
   problem has a unique well-defined optimum (modulo the trivial
   collapse to zero, which we guard against — see §2.7).

So this is not a regression to the old approach — it is the same
geometric idea recast in a closed form that doesn't need a seed.

### 2.7. Edge cases and guards

- **Trivial collapse.** If $\hat\lambda \to 0$ the set $\mathcal A$
  becomes "everywhere" and the iteration is stable at zero. This
  happens if every region has $k_i = 0$ (genuinely no gDNA) or if a
  small number of regions with truly extreme RNA are the only sources
  of counts and zero everywhere else is the equilibrium. Guard: if
  $|\mathcal A(\hat\lambda)| < N_\text{min}$ or $\sum_\mathcal{A} E_i$
  drops below a fraction of the total mappable exposure, fall back to
  $\hat\lambda$ from the previous iteration and emit a calibration
  warning.
- **Hybrid-capture / target-enrichment regime** (see §6) where the
  gDNA rate is non-uniform and elevated inside exons. The global
  estimator under-estimates exonic gDNA. We do **not** attempt to
  handle this in v2; it is flagged as a known limitation with a
  proposed extension (per-locus λ_G).
- **All-zero regions in the index.** Regions with $k_i = 0$ and $E_i
  = 0$ contribute neither to numerator nor denominator and are
  invisible — correct behaviour.

---

## 3. Mappability as a continuous weight

Replace the hard `compute_containment_mask` with a per-region
mappable fraction $f_i \in [0,1]$ computed at index build time:

$$f_i = \frac{|\,\text{region}_i \cap \text{mappable BED}\,|}{|\,\text{region}_i\,|}$$

Stored as a `mappable_fraction` column on `regions.feather`. Default
to $f_i = 1$ when no mappable BED is supplied (preserves behaviour
for indexes built without alignability data).

Index format bumps a major version; users rebuild. No backward
compatibility constraints per project policy.

---

## 4. Combining with the strand pathway — inverse-variance pooling

The strand pathway already produces an estimate $\hat\lambda^{(\text{strand})}$
with a closed-form Fisher information based on the opposite-strand
exonic count under SS. The density pathway (§2) produces
$\hat\lambda^{(\text{density})}$ with $I_\text{density}(\hat\lambda)
= \sum_{\mathcal A} E_i / \hat\lambda$.

Combine via standard inverse-variance pooling:

$$\hat\lambda_G = \frac{I_\text{strand}\,\hat\lambda^{(\text{strand})} + I_\text{density}\,\hat\lambda^{(\text{density})}}{I_\text{strand} + I_\text{density}}$$

Properties:
- At SS = 1.0 the strand pathway dominates ($I_\text{strand}$ is large).
- At SS = 0.5 the strand pathway's information collapses naturally
  ($N_\text{opp}$ has no discriminating power) and the density pathway
  takes over.
- No magic $w = (2\sigma - 1)^2$ blend weight. The blend emerges from
  the Fisher information of each pathway.
- Per-pathway $\hat\lambda$ and Fisher info are reported in
  `summary.json` for diagnostics and consistency checks.

---

## 5. Unmapped reads — research, not infrastructure

For v2 we **do not** add unmapped reads to the calibration likelihood.
Open questions to resolve via offline probes before any production
integration:

1. What is the empirical relationship between $N_\text{uT3}/S_\text{rep}$
   and known $\lambda_G$ across simulated and real datasets? Does it
   track truth or systematically over/under-estimate?
2. How sensitive is that ratio to library type (poly-A, ribo-depleted,
   total RNA), tissue (testis, brain, cancer), and aligner parameters
   (`--outFilterMultimapNmax`)?
3. What fraction of `uT:3` reads are reproducibly RNA in matched
   stranded data (use strand bias to estimate per-tissue $\pi_\text{RNA}$)?
4. Is the distinction between $N_\text{uT3}/S_\text{rep}$ as a **ceiling**
   on $\lambda_G$ (under $\pi_\text{RNA} \ge 0$) more useful than as a
   point estimate?

Action: a future turn will add `scripts/debug/unmapped_pathway_probe.py`
that produces per-condition tables of $(N_\text{uT3}, S_\text{rep},
\hat\lambda_G^\text{truth}, \hat\lambda_G^\text{density},
\hat\lambda_G^\text{strand})$ across the existing benchmark
conditions. We decide its role based on what that data shows.

In the meantime, the BAM scanner can begin **counting** unmapped reads
by `uT` value (cheap, harmless, exposes the data to the Python layer
for the probe). That count flows into `summary.json` as a diagnostic
field but does not enter calibration.

---

## 6. Known limitation — non-uniform gDNA (hybrid capture / WES-from-RNA)

The profile-Poisson estimator (and v1's mixture, and the historical
two-state EM, and P₁₀) all assume gDNA contributes a **globally
uniform** background density $\lambda_G$. In hybrid-capture libraries
(e.g. RNA-seq with exome-enriched probes, or RNA captured by a
gene-panel kit), gDNA contamination is **heavily concentrated inside
the capture targets**. The on-target gDNA rate may be 10–100× the
off-target rate.

Under this model, a globally-fit $\lambda_G$ is meaningless — there is
no single rate. The estimator will report the off-target rate (low),
which is uninformative for the on-target regions where the actual
gDNA-RNA mixing is happening.

Two extensions are possible, both deferred:

1. **Per-locus $\lambda_G$.** Run the profile-Poisson estimator per
   locus over its constituent regions. Requires enough regions per
   locus to be statistical (likely false for most loci).
2. **Two-rate model.** Fit a global "off-target" $\lambda_G^\text{off}$
   on intergenic regions and a global "on-target" $\lambda_G^\text{on}$
   on exonic regions. Use a prior or a known target BED to allocate
   regions. This is essentially modelling the capture protocol.

For v2 we ship the global estimator with a **diagnostic check** that
flags suspicious libraries: if $\hat\lambda^{(\text{intergenic-only})}$
differs from $\hat\lambda^{(\text{exonic-zero-counts})}$ by more than
some factor, emit a "possible target-enrichment artefact" warning.

---

## 7. Implementation plan

Each phase is small and independently mergeable.

### Phase 1 — Index: continuous mappability fraction

1. In `index.py`, compute $f_i$ for every region using the mappable
   BED (already loaded for `mappable.feather`).
2. Add `mappable_fraction` column (float32) to `regions.feather`.
3. Bump index format version. Update `tests/test_index_integrity.py`.
4. Remove `compute_containment_mask` callers (will be unused after
   Phase 2).

### Phase 2 — Calibration: profile-Poisson estimator

1. New function in `calibration.py`:
   ```python
   def density_pathway_lambda(
       counts: np.ndarray,        # k_i, shape (N,)
       exposures: np.ndarray,     # E_i = f_i * S_i, shape (N,)
       max_iter: int = 50,
       tol: float = 1e-6,
       min_regions: int = 100,
   ) -> tuple[float, float, int]:
       """Returns (lambda_hat, fisher_info, n_in_A)."""
   ```
2. Replace `_density_pathway` body (the P₁₀ block in
   `calibrate_gdna`) with a call to `density_pathway_lambda` using
   `counts = n_unspliced`, `exposures = mappable_fraction * region_length`.
3. Compute strand-pathway Fisher information at the existing
   strand-pathway estimate.
4. Replace the $w = (2\sigma - 1)^2$ blend with inverse-variance
   pooling (§4).
5. Augment `CalibrationResult` and `summary.json` with per-pathway
   $\hat\lambda$, Fisher info, and `n_in_A`. Keep
   `lambda_gdna` as the pooled estimate for backward dataflow.

### Phase 3 — Unit tests

1. Synthetic Poisson with known $\lambda_G$, varying $\pi_\text{null}$
   (fraction of regions that are null) and exposure heterogeneity.
   Verify recovery within Cramér-Rao bounds.
2. Adversarial RNA: a few regions with $r_i$ = 10⁵ × $\lambda_G E_i$;
   verify $\hat\lambda_G$ is unaffected.
3. All-zero edge case → graceful $\hat\lambda = 0$ with diagnostic.
4. Hybrid-capture-like adversarial input → estimator returns the
   intergenic rate, target-enrichment warning fires.

### Phase 4 — Validation

1. Re-run the 24-scenario `stress_calibration.py` analytical sweep.
2. Re-run the STAR pipeline-validation BAM (the unstranded + 50% gDNA
   case where the current calibration collapses to zero with the
   mappability mask). Expected outcome: $\hat\lambda_G$ within 30% of
   truth without any tuning.
3. Re-run the full hulkrna benchmark sweep. Confirm no regression on
   stranded high-quality libraries (where strand pathway dominates).

### Phase 5 — Unmapped-read research probe (separate turn)

Per §5 — out of scope for this implementation but the next planning
turn after Phase 4 lands.

---

## 8. What we are not doing (and why)

| Idea (from v1 or earlier) | Decision | Rationale |
|---|---|---|
| Soft-mixture Poisson EM with $\pi_i$ posteriors | Dropped | Hard-EM (§2) is the closed-form limit; no information loss for our use case |
| Tile-based denominators | Rejected | Region partition is the index's natural quantum; tiling breaks the architecture |
| Joint likelihood combining unmapped pathway as a likelihood term | Deferred | π_RNA misspecification is a real bias risk; needs offline probe first |
| Hard mappability containment mask | Replaced | Continuous $f_i$ exposure weight is strictly more informative |
| `density_percentile` parameter | Removed | No percentile in the new estimator; CLI flag retained for one release as a no-op with deprecation warning |
| Heuristic blend weight $w = (2\sigma-1)^2$ | Replaced | Inverse-variance pooling is principled and self-silences each pathway when uninformative |

---

## 9. Summary

The new calibration is one short function (the iteration in §2.3),
plus inverse-variance pooling with the existing strand pathway. It
uses the existing region partition (no tiles), handles zero counts as
informative evidence, is immune to the RNA dynamic range, and has no
arbitrary percentile. The unmapped-read pathway becomes a research
diagnostic to be evaluated in a separate turn before any production
integration. Hybrid-capture / target-enrichment libraries are flagged
as a known limitation with a documented extension path.
