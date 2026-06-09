# Count-overdispersion integration — implementation plan (clean slate)

**Status:** Phase A + Phase B **IMPLEMENTED** (2026-06-09; Phase A bit-identical, Phase B golden
regenerated, full suite green). Phase C (end-to-end net-flow benchmark) pending. Supersedes the
exploration in `clue_integration_reliability_design.md` (which logs the vetoed count-reliability
idea and the superseded strand-reliability weight).

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
# ... compute `density` (local imputation) + masks; CHANGE: no-anchor regions shrink to the GLOBAL
#     gDNA density (not 0) — a sensible baseline (decision §6.6) — flagged via `has_local_evidence`.
mass = substrate.contained.mass_unspliced                      # contained unspliced mass per region
with errstate:                                                 # the count-prior MEAN (was in the deconv)
    count_gdna_frac = clip(where(mass>0, density*region_eff_len/mass, 0.0), 0, 1)
# count_support = the HONEST gDNA evidence behind the estimate (decision §6.1):
#   observable region : the contained gDNA count            (clean_gf · contained_total)
#   imputed region    : the ANCHORING crossing gDNA count   (Σ exon-side clean crossing flux used)
#   no local evidence : ~0 (global-density baseline ⇒ minimal concentration ⇒ defer)
return NodeDensity(density, count_gdna_frac, count_support, region_count_observable,
                   boundary_count_observable, ...)
```

**`count_dispersion.fit_gdna_count_overdispersion`** (mirror `gdna_strand`):
```python
def fit_gdna_count_overdispersion(substrate, ra, node_density, region_eff_len, fl_mean, *, κ,
                                  gdna_density_global, prior_alpha, prior_weight):
    # Two seed sets = two SUPPORTING-COUNT TYPES (decision §6.2 — NOT arbitrary):
    #   CONTAINED: count-observable, non-AMBIG regions. N = clean_gdna_contained ; μ = ρ̄·region_eff_len
    #   CROSSING : count-observable exon-adjacent boundary sides. N = clean_gdna_cross ; μ = ρ̄·fl_mean
    #   (ρ̄ = gdna_density_global, the COMMON expectation, so residuals capture heterogeneity)
    #   Reuse the strand fit's seed geometry (_region_seeds, boundary_side_seeds).
    def nb_mom(N, mu):                     # proper NB method-of-moments (handles small counts, §6.3)
        # Var = μ + α·μ²  ⇒  α̂ = Σ[(N-μ)² - μ] / Σ μ²   (Poisson term subtracted per seed)
        return max(sum((N-mu)**2 - mu) / sum(mu**2), 0.0)
    a_contained = shrink(nb_mom(reg_N, reg_mu), W(reg), prior_alpha, prior_weight)
    a_crossing  = shrink(nb_mom(bnd_N, bnd_mu), W(bnd), prior_alpha, prior_weight)
    return a_contained, a_crossing
# shrink(â, W, a0, w) = (W·â + w·a0)/(W + w)   # EB; sparse ⇒ a0 (conservative-high) dominates (§6.4)
# W(seeds) = seed evidence weight (n_seeds or Σμ); prior = "unreliable until proven" (§6.4, OPEN)
```

**`calibrate`** (insert the fit; finalize `count_evidence`):
```python
node_density = node_gdna_density(substrate, ra, region_eff_len, fl_mean, rna_sense_frac=κ)
a_contained, a_crossing = fit_gdna_count_overdispersion(
    substrate, ra, node_density, region_eff_len, fl_mean, κ=κ,
    gdna_density_global=node_density_global_density,           # the COMMON expectation ρ̄
    prior_alpha=cfg.count_overdispersion_prior,
    prior_weight=cfg.count_overdispersion_prior_weight)
alpha = where(node_density.region_count_observable, a_contained, a_crossing)  # α of the support type
N = node_density.count_support                  # honest gDNA evidence (contained or crossing)
count_evidence = N / (1.0 + alpha * N)          # overdispersion-limited effective count
# strand-overdispersion fit unchanged (now reads node_density.count_gdna_frac for its weight)
# deconv_regions/deconv_sides: pass count_evidence as the concentration, count_gdna_frac as the mean
```

**`joint_deconv`**: `_joint_per_node` takes `count_gdna_frac` (mean) directly instead of recomputing
`density·eff_len/mass`; `count_concentration = count_evidence` (now overdispersion-limited).
`_region_seeds`: `weight = clip(node_density.count_gdna_frac, 0, 1)`.

## 6. Decisions (resolved) + the one collaborative item

1. **Supporting count = the HONEST evidence** *(decided)*. Observable region → its contained gDNA
   count; imputed region → the anchoring **crossing** gDNA count (not `density·eff_len`). For
   zero-evidence regions, see §6.6.
2. **Two α's — one model, two count TYPES** *(decided)*. `α_contained` (contained-count seeds) for
   observable regions; `α_crossing` (crossing-count seeds) for imputed regions. They differ ~30× under
   capture (crossing carries the enrichment heterogeneity — the leak-closing signal); a single pooled
   α would mis-serve both. Each node uses the α of the count *type* that supports it.
3. **α estimator** *(decided)* — the proper NB MoM `Σ[(N−μ)²−μ]/Σμ²` (μ = global-density expectation),
   which subtracts the Poisson term per seed and so is well-behaved for the **many small-count seeds
   of real data** (the toy data's high coverage hid this — the global `CV²` form would be biased
   there). DESeq2 mean-trend `α(μ)=asymptDisp+extraPois/μ` is the documented fallback if α varies with
   μ; our leak-critical regime is the high-count asymptote, so start without it and add if a
   small-count scenario shows mean-dependence.
4. **THE COLLABORATIVE ITEM — the overdispersion prior** *(RESOLVED 2026-06-09, with the user)*.
   **Decision:** fallback anchor **α₀ = 1** (the *geometric* NB — max-entropy at fixed mean;
   `N_eff → 1`, so absent dispersion evidence a count is worth ~1 pseudo-observation — the count
   analog of the Jeffreys `Beta(½,½)` floor; conservative/leak-averse). Shrink each per-type MoM
   toward the **global pooled-seed trend** (NB-MoM over all seeds combined, DESeq2-style), with α₀
   used **only** when that pool is degenerate (no usable seeds). Strength `count_overdispersion_prior_weight`
   in seed-node units, default **30** (mirrors the strand fit); negligible with the abundant seeds
   of a normal library (the distinct contained/crossing α's are preserved — that distinctness is
   what closes the leak), a guard only in the degenerate few-seed case. *Original note below kept for
   provenance.*
   **Key simplification (from the precedent survey):** edgeR (`prior.df`) and DESeq2 carry elaborate
   empirical-Bayes shrinkage *because they fit a dispersion per gene* (few replicates each, so each
   estimate is noisy and must borrow from a trend). **We fit a single *common* α per count-type,
   pooled over all seeds** — so even with the many small-count seeds of real data, the pooled NB-MoM
   `Σ[(N−μ)²−μ]/Σμ²` is well-determined by averaging, and per-feature shrinkage is *not* needed. The
   prior is therefore only a **guard for the degenerate few-seed library** (too few seeds to estimate
   α at all). Principle (rigel's strand precedent *sparse ⇒ high overdispersion*): that guard is
   **conservative-high** — "unreliable until proven reliable" — so a data-poor library defers to
   strand. EB shrinkage `α = (W·α̂ + w·a0)/(W + w)` with `W` = total seed evidence; when `W ≫ w` (the
   normal case) the data dominates and the prior is irrelevant. **To pin together:** `a0` + `w` — but
   the surface is small (it only bites at very low seed counts). Validate on a **low-coverage /
   few-seed** scenario, not just the toy high-coverage benchmark. **Do NOT pick a magic number
   unilaterally.** *Note: Phase A does not depend on this; only Phase B does.*
5. **Unstranded + capture = worst case** *(accepted)* — strand blind AND count locally biased; push
   forward, handle later if needed (the DNA-fraction / defect B is the eventual lever).
6. **Zero-evidence regions → shrink to the GLOBAL gDNA density** *(decided, per your steer)*. A region
   with no anchor (no own crossing, no run-fill reach) takes `density = gdna_density_global` (a
   sensible baseline, not 0) with a **minimal `count_support`** (it is a global fallback, not local
   evidence) ⇒ small concentration ⇒ it defers (to strand where present, else a weak global baseline).
   Run-filled regions (carry reaches) keep the carried density with the carried (attenuated) support.
7. **RNA twin — deferred** *(accepted)*. The same machinery could model RNA count overdispersion from
   spliced-fragment counts later; out of scope.

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
