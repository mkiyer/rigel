# Count-overdispersion integration — implementation plan (clean slate)

**Status:** implementation plan + dry-run. 2026-06-09. The grounded foundation for the chosen
clue-integration solution. Supersedes the exploration in `clue_integration_reliability_design.md`
(which logs the vetoed count-reliability idea and the superseded strand-reliability weight).

---

## 1. Goal & principle

Give the count clue an **honest precision** so the joint deconvolution self-weights count vs strand
correctly across the strand-specificity × capture matrix — replacing the categorical
`defer_to_strand` zeroing with one principled mechanism.

The count-prior concentration is currently the **raw expected gDNA count** `density·eff_len`
(thousands) — the *Poisson* precision. RNA-seq / gDNA counts are **NB-overdispersed**, so the honest
precision is the **overdispersion-limited effective count** `N_eff = N / (1 + α·N)` (→ `1/α` for large
`N`), where `α` is the gDNA count overdispersion fitted from the data — the count-side twin of the
existing strand Beta-Binomial overdispersion.

**Validated** (`α` via DESeq2-style common-dispersion MoM `CV²(ρ) − mean(1/N)` on gDNA counts from
count-observable **regions** + exon-adjacent observable **boundaries**; `count_evidence = N/(1+αN)`,
no strand-reliability weight, no zeroing):

| condition | boundary α | oracle gf | deconv gf |
|---|---|---|---|
| capture / ss0.99 | 1.21 (eff-N ≈ 1) | 0.804 | **0.795** — α huge ⇒ count fades ⇒ strand governs ⇒ leak closed |
| no-capture / ss0.50 | 0.0002 (eff-N ≈ 5000) | 0.963 | **0.883** — α tiny ⇒ count carries ⇒ weak-SS floor 0.60→0.88 |

The **boundaries** carry the capture enrichment heterogeneity (the imputation source); the
intron/intergenic **regions**, uniformly depleted, do not (α≈0.04 even under capture) — so both seed
types are needed, regions for count-observable regions, boundaries for imputed ones.

## 2. Current data flow + the two couplings to untangle

```
node_gdna_density → NodeDensity{density, count_evidence(=density·eff_len, zeroed for single-strand exon),
                                region_count_observable, boundary_count_observable}
   ├─ joint_deconv._joint_per_node:   count_gdna_frac = clip(density·eff_len/mass)   (the prior MEAN)
   │                                  count_concentration = count_evidence            (the prior PRECISION)
   ├─ gdna_strand._region_seeds:      weight = clip(count_evidence/mass)              (= the MEAN, NOT precision!)
   └─ deconv_sides (region_density = node_density.density, fallback)
```

- **Coupling A (must fix):** `count_evidence` is double-duty — the deconv uses it as the *precision*,
  but `_region_seeds` uses `count_evidence/mass` as the strand-fit *weight* (which is really the count
  **mean**). Changing `count_evidence` to `N_eff` (much smaller) would collapse the strand-fit weight.
  ⇒ **decouple the mean from the precision**: add an explicit `count_gdna_frac` (the mean) to
  NodeDensity; the deconv mean and the `_region_seeds` weight both read it; `count_evidence` becomes
  the precision only.
- **Coupling B (ordering):** `α` is fit from the seed gDNA counts **and** the density, so it must run
  *after* `node_gdna_density`. ⇒ `count_evidence` cannot be finalized inside `node_gdna_density`; it is
  computed after the α fit. So `node_gdna_density` returns the raw **supporting count** `count_support`
  and the α step turns it into `count_evidence = count_support/(1+α·count_support)`.

## 3. Target data flow

```
node_gdna_density → NodeDensity{density, count_gdna_frac (MEAN), count_support (raw N),
                                region_count_observable, boundary_count_observable}
fit_gdna_count_overdispersion(substrate, region_arrays, node_density, region_eff_len, fl_mean, κ)
        → (α_region, α_boundary)                                   # MoM + shrinkage; mirrors the strand fit
count_evidence[r] = count_support[r] / (1 + α[r]·count_support[r]) # α_region if observable else α_boundary
joint_deconv: mean = node_density.count_gdna_frac ; precision = count_evidence
gdna_strand._region_seeds: weight = node_density.count_gdna_frac   # decoupled from count_evidence
```

## 4. Phased plan

- **Phase A — decouple mean/precision (bit-identical).** Add `count_gdna_frac` + `count_support` to
  NodeDensity; `deconv_regions` reads `count_gdna_frac` (instead of recomputing `density·eff_len/mass`);
  `_region_seeds` reads `count_gdna_frac` (instead of `count_evidence/mass`); keep `count_evidence`
  equal to the **current** value (`density·eff_len`, with the existing zeroing) so the run is
  bit-for-bit. Verify golden unchanged + suite green. *(This is the long-deferred explicit-`count_gdna_frac`
  cleanup; it lands risk-free here.)*
- **Phase B — count overdispersion.** Add `fit_gdna_count_overdispersion` (new module
  `count_dispersion.py`, mirroring `gdna_strand.py`); set `count_evidence = count_support/(1+α·count_support)`;
  **remove** the `defer_to_strand` zeroing (the overdispersion-honest precision replaces it; AMBIG needs
  no special-case). Add config knobs (default α + prior weight). Validate on the SS×capture matrix +
  the antisense stress; regenerate golden.
- **Phase C — validate end-to-end.** Re-quant the benchmark (+ an antisense+gDNA condition) and run the
  net-flow analysis; confirm transcript- and dataset-level gDNA↔RNA error.

## 5. Dry-run pseudocode

**`density_model.node_gdna_density`** (return mean + support; drop count_evidence/zeroing):
```python
# ... unchanged: compute `density` (local imputation), masks ...
mass = substrate.contained.mass_unspliced                      # contained unspliced mass per region
with errstate:                                                 # the count-prior MEAN (was in the deconv)
    count_gdna_frac = clip(where(mass>0, density*region_eff_len/mass, 0.0), 0, 1)
count_support = ??? # OPEN Q (§6.1): density*region_eff_len  vs  the anchoring crossing gDNA count
return NodeDensity(density, count_gdna_frac, count_support, region_count_observable, boundary_count_observable, ...)
```

**`count_dispersion.fit_gdna_count_overdispersion`** (mirror `gdna_strand`):
```python
def fit_gdna_count_overdispersion(substrate, ra, node_density, region_eff_len, fl_mean, *, κ,
                                  prior_alpha, prior_weight):
    # REGION seeds: count-observable, non-AMBIG. ρ = clean_gdna_contained / region_eff_len ; N = clean_gdna_contained
    # BOUNDARY seeds: count-observable exon-adjacent sides. ρ = clean_gdna_cross / fl_mean ; N = clean_gdna_cross
    #   (reuse the strand fit's seed collection: _region_seeds geometry + boundary_side_seeds)
    def mom(rho, N):                       # common-dispersion NB MoM (CV² minus Poisson), guarded
        cv2 = var(rho)/mean(rho)**2
        return max(cv2 - mean(1.0/N), 0.0)
    a_reg = shrink(mom(reg_rho, reg_N),  n_reg,  prior_alpha, prior_weight)   # EB shrinkage to prior
    a_bnd = shrink(mom(bnd_rho, bnd_N),  n_bnd,  prior_alpha, prior_weight)   #   when few seeds
    return a_reg, a_bnd
# shrink(â, n, a0, w) = (n·â + w·a0) / (n + w)   # OPEN Q (§6.3): exact shrinkage / DESeq2 trend
```

**`calibrate`** (insert the fit; finalize `count_evidence`):
```python
node_density = node_gdna_density(substrate, ra, region_eff_len, fl_mean, rna_sense_frac=κ)
a_reg, a_bnd = fit_gdna_count_overdispersion(substrate, ra, node_density, region_eff_len, fl_mean,
                                             κ=κ, prior_alpha=cfg.count_overdispersion_prior,
                                             prior_weight=cfg.count_overdispersion_prior_weight)
alpha = where(node_density.region_count_observable, a_reg, a_bnd)
N = node_density.count_support
count_evidence = N / (1.0 + alpha * N)          # overdispersion-limited effective count
# strand-overdispersion fit unchanged (now reads node_density.count_gdna_frac for its weight)
# deconv_regions/deconv_sides: pass count_evidence as the concentration, count_gdna_frac as the mean
```

**`joint_deconv`**: `_joint_per_node` takes `count_gdna_frac` (mean) directly instead of recomputing
`density·eff_len/mass`; `count_concentration = count_evidence` (now overdispersion-limited).
`_region_seeds`: `weight = clip(node_density.count_gdna_frac, 0, 1)`.

## 6. Open questions / decisions (resolve before/while implementing)

1. **Supporting count `count_support` for IMPUTED regions** — `density·eff_len` (the expected
   *contained* gDNA count, used in the validation, worked) vs the **anchoring crossing gDNA count**
   (the actual imputation *evidence*). They coincide under capture (α saturates `N_eff`→1/α regardless
   of `N`) but differ at no-capture/low-SS (sub-saturated). Lean: the crossing count (the honest
   evidence). **Decide empirically in Phase B** (re-check the no-cap/ss0.50 number both ways).
2. **`α_region` vs `α_boundary` assignment** — observable regions → `α_region`; imputed (exon/AMBIG)
   regions → `α_boundary`. (Confirm; it falls straight out of the masks.)
3. **α estimator + shrinkage** — start with the guarded common-dispersion MoM `CV²−Poisson` + linear
   EB shrinkage to a prior (mirroring `gdna_strand_prior_weight`). Optional later: the DESeq2 mean-trend
   `α(μ)=asymptDisp+extraPois/μ` if small-count seeds bias the common estimate (our boundaries are
   high-count, so the asymptotic α dominates — likely unnecessary).
3a. Weighting of the MoM across seeds (equal vs by-`N`); robustness to outlier seeds (trim?).
4. **Default `α` prior + prior_weight** (the few-seed fallback) — needs a principled default, not a
   magic constant. A *weak* prior (low weight) lets data dominate where seeds are plentiful; but with
   *few* seeds under capture, under-estimating α would over-trust a biased count → leak. So the default
   should be **conservative (moderate α)** so a data-poor library does not over-trust the count. Pin
   the value + weight from the SS×capture sweep. **Open.**
5. **AMBIG + capture floor** — unchanged by this work (α_boundary huge ⇒ count fades ⇒ strand absent ⇒
   ~0.5). Genuine; the DNA-fraction (path 3, defect B) remains the future de-bias. Confirm acceptable.
6. **`count_support` = 0 / run-filled regions** — a region with no anchor (density 0) → support 0 →
   count_evidence 0 (Jeffreys). A run-filled region inherits a carried density; its `count_support`
   should reflect the (indirect) carried evidence — set to the carried support or a low value. Specify.
7. **RNA twin (deferred)** — the same machinery could fit an RNA count overdispersion from
   spliced-fragment counts at boundaries; out of scope here, noted.

## 7. Tests & acceptance

- **Phase A bit-identical:** golden unchanged; full suite green; `count_gdna_frac`/`count_support`
  unit checks; `_region_seeds` weight equals the old `count_evidence/mass`.
- **Phase B unit:** `mom` recovers a known α on synthetic over-/under-dispersed counts; `N/(1+αN)`
  monotone, → 1/α for large N; `fit_gdna_count_overdispersion` returns small α on uniform seeds, large
  α on heterogeneous seeds; AMBIG excluded from the fit.
- **Phase B worked (oracle):** capture/ss0.99 exon deconv gf ≈ oracle (≈0.80) with NO zeroing/weight;
  no-cap/ss0.50 single-strand exons recover (≈0.88, not 0.60); gdna_none clean at all SS; antisense
  AMBIG correct at no-capture.
- **Golden + suite** regenerate (calibration concentration changes); full suite green.
- **Phase C:** end-to-end net-flow leak across the matrix + antisense condition.

## 8. Files touched
- `density_model.py` — NodeDensity gains `count_gdna_frac` + `count_support`; drop `count_evidence`
  field + `defer_to_strand`.
- `count_dispersion.py` (new) — `fit_gdna_count_overdispersion` (+ MoM/shrinkage helpers), mirroring
  `gdna_strand.py`; reuse the seed geometry.
- `joint_deconv.py` — `_joint_per_node` takes the mean; `_region_seeds` weight ← `count_gdna_frac`.
- `calibrate.py` — insert the α fit; compute `count_evidence = N/(1+αN)`.
- `config.py` — `count_overdispersion_prior` + `count_overdispersion_prior_weight`.
- tests + golden.
