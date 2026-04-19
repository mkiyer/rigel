# Unified gDNA Calibration — Proposal v3

_Date: 2026-04-18 · Supersedes `unified_proposal_v2.md`_

## Changes from v2

v2 settled the core estimator (profile-Poisson over the existing
region partition) and the pathway-blending strategy (inverse-variance
pooling). v3 makes two substantive additions and one clarification:

1. **Mappability becomes a continuous 1/NH-style score**, sourced from
   `alignable`'s Zarr per-base landscape rather than a binary BED.
   This aligns the exposure denominator $E_i = f_i S_i$ with rigel's
   multimapper credit-splitting convention, so $f_i$ has a clean
   probabilistic interpretation.
2. **Hybrid-capture support is promoted from "deferred limitation" to
   a first-class target-BED-aware two-rate density estimator.** This
   is needed soon and is small.
3. The unmapped-read pathway and the pooled inverse-variance blend
   from v2 are retained as written; the §2 estimator is unchanged.

---

## 1. Mappability as a continuous per-region weight

### 1.1. The right thing for $f_i$ to mean

The estimator's denominator is $E_i = f_i \cdot S_i$. For this to be
the correct **expected gDNA exposure** under a uniform fragmentation
rate $\lambda_G$, we need:

$$E_i = \int_{b \in \text{region}_i} \Pr[\text{a gDNA fragment originating at base } b \text{ contributes 1 unit of count to region } i] \, db$$

So $f_i$ is the *average per-base probability that a fragment is
"observable" inside region $i$*. The relevant per-base spectrum is:

| Base type | Per-base contribution to $f_i$ |
|---|---|
| Uniquely mappable (NH = 1) | 1.0 |
| Multimapping with $NH = K$, **rigel splits credit $1/K$** | $1/K$ (this region gets $1/K$, the other $K-1$ regions also get $1/K$ each, total = 1) |
| Above aligner tolerance (NH > NH_max → unmapped) | 0.0 |

The natural per-base score is therefore **1/NH where NH is the number
of tile placements at that base under the aligner's mapping
tolerance**, with NH = ∞ → 0. This is exactly what GENMAP/UMAP-style
tools produce, and what `alignable`'s continuous landscape is
expected to provide (per-base float in $[0, 1]$).

Defining

$$f_i = \frac{1}{S_i}\sum_{b \in \text{region}_i} \text{mappability}(b)$$

makes $f_i \in [0, 1]$ a **mean per-base mappability**, and the
exposure $E_i = f_i S_i$ reduces to the expected number of
"observable bases" in the region. This is the right denominator for
the profile-Poisson estimator, period.

### 1.2. Why the continuous score beats either BED

| Source | What it represents | Pathology |
|---|---|---|
| Uniquely-mappable BED | $1/NH = 1$ vs. 0 (binary) | Throws away multimapping bases entirely. A region that is 100% NH=2 would get $f_i = 0$ instead of the correct 0.5, badly under-counting exposure |
| Unmappable BED | NH > NH_max (binary) | Treats NH=2 and NH=1 identically as "in", over-counting exposure for repeat-rich regions |
| Continuous per-base 1/NH | The actual probability | Correct |

The continuous score is also strictly more informative: both BEDs are
projections of it (`uniquely-mappable BED = {b : score(b) = 1}`,
`unmappable BED = {b : score(b) ≈ 0}`).

### 1.3. The read-length problem

`alignable` computes mappability at a chosen read length $L$. We
don't know the data's read length at index build time. Three options:

**A. Fixed conservative read length.** Build the index at $L = 50$.
Pros: matches the existing "minimum supported read length" floor;
single value to store. Cons: $f_i$ is a lower bound on the true
exposure for typical 100–150 bp data, biasing $\hat\lambda_G$
**upward** (denominator too small).

**B. Fixed mid-range read length (recommended).** Build at $L = 100$
or a CLI-selectable value, default 100. Pros: matches the modal RNA-seq
read length, mappability differences across 75–150 bp are small for
euchromatin (the regime where the density pathway operates), tiny
implementation. Cons: significantly off for very long-read or very
short-read data.

**C. Multiple read-length tables, select at quant time.** Pros:
exact. Cons: index storage roughly $K\times$, selection plumbing,
interpolation logic, marginal gain for the typical 75–150 bp range.

**Decision: B.** Make read length a CLI flag at index time
(`--mappability-read-length`, default 100). Document the assumption
in the index manifest. Rigel emits a soft warning at quant time if
the observed median read length differs from the index-build value
by more than ~30%.

The estimator is robust to modest exposure mis-calibration: a 10%
error in $f_i$ produces a ~10% error in $\hat\lambda_G$, which is
small compared to the order-of-magnitude errors the current P₁₀
estimator can exhibit.

### 1.4. Source format priority

`alignable` can produce BED, bigwig, or Zarr. Priority for v3:

1. **Zarr (preferred when available).** Per-base landscape; computes
   $f_i$ exactly via array slicing. Ship a small loader that takes a
   region table and returns $f_i$ for each row.
2. **BigWig fallback.** Same data, file-based; load with `pyBigWig`
   if Zarr is absent.
3. **BED fallback.** If only a uniquely-mappable BED is supplied, set
   $f_i = $ fraction of region bases inside the BED — i.e. treat
   multimapping bases as fully unmappable. Document this as a coarse
   fallback that will under-estimate exposure and bias $\hat\lambda_G$
   upward in repeat-rich loci.

Index stores $f_i$ as a single `mappable_fraction` column on
`regions.feather` (float32). The source artifact (Zarr/bigwig/BED
path + read length) is recorded in the index manifest for
traceability.

### 1.5. Edge cases

- **No mappability source supplied.** Default $f_i = 1$ for every
  region — equivalent to assuming the whole genome is uniquely
  mappable. Strictly worse calibration in repeat-rich loci, but the
  estimator still runs. A prominent log message warns the user.
- **Region with $f_i = 0$.** Zero exposure → contributes nothing to
  numerator or denominator → silently irrelevant. Correct.
- **Region with $f_i$ very low (e.g. 0.05) but $k_i$ moderate.** Will
  almost certainly be in $\mathcal A$ (since $\lambda E_i$ is tiny)
  and contribute its $k_i$ to the numerator. This is the most
  pathological case: a repeat-rich region with a few mapped fragments
  that look "above background" in raw count but only because exposure
  is so low. Mitigation: optional minimum-$f_i$ threshold (default
  0.1) below which a region is excluded from $\mathcal A$. This is a
  soft guard and the only "magic constant" in the new estimator.

---

## 2. Profile-Poisson estimator (unchanged from v2)

Verbatim restatement for completeness; see v2 §2 for the full
derivation.

For region $i$ with exposure $E_i = f_i S_i$ and unspliced count $k_i$:

```
λ ← total counts / total exposure
repeat:
    A ← { i : k_i < λ · E_i  AND  f_i ≥ f_min }
    λ_new ← Σ_{A} k_i / Σ_{A} E_i
    if |λ_new − λ| < tol: break
    λ ← λ_new
```

Fisher information at the fixed point: $I(\hat\lambda) = \sum_{\mathcal A} E_i / \hat\lambda$.

---

## 3. Hybrid-capture / target-BED awareness — two-rate density estimator

### 3.1. The problem

Hybrid-capture and exome-enriched RNA-seq libraries do not satisfy
the uniform-$\lambda_G$ assumption. On-target gDNA contamination can
be 10–100× the off-target rate, because the same probes that pull
down RNA also pull down genomic fragments overlapping the target
intervals.

A single global $\lambda_G$ fit across all regions would converge to
something between the two rates and be useful for neither. In
particular, on-target regions (which is where the gDNA-vs-RNA
deconvolution actually matters) would get a $\lambda_G$ that is far
too low, so the density pathway would systematically under-correct
their gDNA contribution.

### 3.2. The fix

Accept an optional **target BED** (probe / capture interval set) at
quant time (`--target-bed`). At calibration time:

1. Mark every region with `is_target = (region overlaps target BED by
   ≥ T bp)`. Default $T = 1$ (any overlap). The target BED is loaded,
   sorted, merged once per run.
2. Run the profile-Poisson estimator (§2) **twice**, on disjoint
   subsets:
   - On-target regions → $\hat\lambda_G^{\text{on}}$, $I^{\text{on}}$
   - Off-target regions → $\hat\lambda_G^{\text{off}}$, $I^{\text{off}}$
3. The per-region gDNA rate becomes
   $\lambda_G(i) = \hat\lambda_G^{\text{on}}$ if `is_target[i]` else
   $\hat\lambda_G^{\text{off}}$.

The downstream consumers (`region_e_gdna`, locus priors) take a
**per-region** rate vector instead of a scalar. Currently they
already accept a vector multiplied by region length; the only change
is the source.

### 3.3. Behaviour without a target BED

When `--target-bed` is not supplied, `is_target = False` everywhere
and the off-target estimator runs over all regions — the v2 behaviour
is recovered exactly. This is the default for vanilla RNA-seq.

### 3.4. Diagnostic

When a target BED is supplied, report
$\hat\lambda_G^{\text{on}} / \hat\lambda_G^{\text{off}}$ in
`summary.json`. A ratio close to 1.0 indicates the hybrid-capture
artefact is absent (probably wrong BED or non-capture library); a
ratio of 10–100 confirms the artefact and validates the two-rate
treatment.

When `--target-bed` is **not** supplied but the data is hybrid
capture, we want to detect it. Heuristic: split regions into "exonic"
and "intergenic" using the existing index annotations; if
$\hat\lambda_G^{\text{exonic-zero-counts-only}} \gg \hat\lambda_G^{\text{intergenic}}$
by more than (say) 5×, emit a "possible target enrichment — consider
supplying `--target-bed`" warning. This is the same diagnostic
flagged in v2 §6, now actionable.

### 3.5. Strand pathway under hybrid capture

The strand pathway is unchanged. Hybrid-capture libraries are
typically strand-specific, so the strand pathway already dominates
the inverse-variance pool for those libraries. The two-rate density
extension exists to handle the (less common but real) **unstranded**
hybrid-capture case where the density pathway is the only signal.

### 3.6. Why this small extension is sufficient

A more general "spatially-varying $\lambda_G$" model (e.g., a
Gaussian process over genomic position) would be statistically
elegant but operationally heavy. The two-rate model captures the
essential structural feature of hybrid capture (a dichotomy of rates
between probed and unprobed regions) at minimal cost: one extra BED
file, one extra estimator call, one extra summary field. Per-locus
rates remain as a future extension if hybrid-capture libraries with
heterogeneous probe efficiency become a problem.

---

## 4. Pathway pooling (unchanged from v2)

Inverse-variance pool the strand and density pathways using each
pathway's Fisher info. Under the two-rate density estimator, the
pooling is done **per-target-class**:

- On-target regions: pool $\hat\lambda_G^{\text{on, density}}$ and
  $\hat\lambda_G^{\text{on, strand}}$ (the strand pathway is fit
  separately on on-target regions if a target BED is supplied).
- Off-target regions: pool $\hat\lambda_G^{\text{off, density}}$ and
  $\hat\lambda_G^{\text{off, strand}}$.

Without a target BED, this collapses to the single-pool v2 behaviour.

---

## 5. Implementation plan (revised)

### Phase 0 — Alignable integration sanity check

Before any code changes, write a short notebook /
`scripts/debug/alignable_landscape_probe.py` that:
- Loads the existing Zarr (or bigwig) landscape from `alignable` for
  GRCh38 at $L = 100$.
- Computes $f_i$ for the existing rigel index regions.
- Reports the distribution of $f_i$, scatter against the
  uniquely-mappable-BED-fraction, and a few sanity checks (centromere
  → $f_i \approx 0$, exonic → $f_i \approx 1$).

This validates that the per-base score is interpretable as 1/NH and
that mean-aggregation gives sensible region scores before we wire it
into the index.

### Phase 1 — Index changes

1. Add `mappable_fraction` (float32) to `regions.feather`, computed
   from the chosen mappability source (Zarr / bigwig / BED, in that
   order).
2. CLI flag at index build: `--mappability {zarr,bigwig,bed,none}`,
   `--mappability-path PATH`, `--mappability-read-length INT`
   (default 100). Record all of these in the index manifest.
3. Bump index format version. Update affected tests.
4. Remove the old `compute_containment_mask` and the binary-BED
   loading code path (or keep as the BED fallback only).

### Phase 2 — Calibration: profile-Poisson + two-rate

1. Implement `density_pathway_lambda(counts, exposures, ...)` per v2
   §2 / §2.7.
2. Add target-BED loading utility (reuse `_bed_io.read_bedlike`).
3. Compute `is_target` per region from the target BED (intersection
   with region intervals via cgranges, threshold $T$).
4. In `calibrate_gdna`, branch:
   - No target BED → single estimator over all regions (v2 path).
   - Target BED → run on-target and off-target estimators, return a
     per-region $\lambda_G$ vector.
5. Inverse-variance pool with strand pathway per §4.
6. Augment `CalibrationResult` to carry both rates and the boolean
   target mask. Update `summary.json` schema.

### Phase 3 — Tests

1. Synthetic uniform-gDNA benchmark, recovery within Cramér-Rao.
2. Synthetic adversarial RNA dynamic range.
3. Synthetic hybrid-capture: gDNA rate 10× higher in target regions
   → estimator recovers both rates within tolerance, ratio matches
   simulation.
4. No-mappability-source test → graceful $f_i = 1$ behaviour.
5. Read-length-mismatch warning fires when applicable.

### Phase 4 — Validation

1. Re-run the 24-scenario `stress_calibration.py` sweep — confirm no
   regression on uniform-gDNA conditions.
2. Re-run the STAR pipeline-validation BAM (unstranded + 50% gDNA).
3. Identify or simulate a hybrid-capture-like dataset and verify the
   two-rate estimator behaves sensibly.
4. Full hulkrna benchmark sweep.

### Phase 5 — Unmapped-read research probe (separate turn, unchanged)

Per v2 §5.

---

## 6. Outstanding questions

1. **`alignable` Zarr API.** What is the current schema (chromosome
   keying, dtype, chunk size)? The Phase 0 probe will surface this.
   Do we need to vendor a thin wrapper, or is the existing public API
   stable?
2. **Default read length.** 100 bp is a reasonable modal value, but
   recent datasets are skewing 150 bp. Worth checking the
   distribution across the `hulkrna` benchmark inputs before
   committing to a default.
3. **Minimum-$f_i$ threshold (`f_min`).** Default 0.1 is a guess. A
   quick sensitivity analysis over `[0, 0.05, 0.1, 0.2, 0.5]` on
   simulated data should pin this down.
4. **Target BED format.** Standard BED3 is sufficient; should we also
   accept a probe interval list with metadata (probe ID, capture
   efficiency)? For v3, BED3 only — extensions deferred.

---

## 7. Summary

| Topic | v3 decision |
|---|---|
| Mappability source | Continuous per-base 1/NH-style score from `alignable` Zarr (with bigwig/BED fallbacks); aggregate to per-region mean $f_i$ |
| Read length | Single CLI-configurable value at index time, default 100 bp; soft mismatch warning at quant |
| Estimator | Profile-Poisson over regions, exposure $E_i = f_i S_i$ (unchanged from v2) |
| Hybrid capture | Optional `--target-bed`; two-rate density estimator (on-target / off-target $\lambda_G$); per-region vector returned |
| Pathway pooling | Inverse-variance with per-pathway Fisher info, applied per-target-class when a target BED is supplied |
| Unmapped reads | Counted in BAM scanner as diagnostic; production integration deferred to a future turn |

The estimator core stays a one-page function. Mappability becomes
probabilistic instead of binary, matching how multimappers actually
contribute to counts. Hybrid capture is supported with one extra BED
input and one extra estimator call. Everything else from v2 stands.
