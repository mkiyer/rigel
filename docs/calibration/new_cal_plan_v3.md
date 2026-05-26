# New Calibration Plan v3

Date: 2026-05-25
Status: implementation-ready design (incorporates v2 review)
Supersedes: `new_cal_plan_v2.md`, `new_cal_plan_v1.md`
Resolved against user review of v2 open questions Q1–Q10.

This document keeps the v2 latent-variable framing and module layout, but bakes in the resolved answers, promotes one v2 sub-point ("strand finds captured-but-not-expressed regions") to a first-class signal, and replaces the fixed two-pass loop with a convergence-gated iteration. A new section §7 enumerates residual issues found while sharpening the design.

---

## 0. North Star (unchanged from v2)

Per fine region `r`, calibration learns two continuous latents:

```text
e_r = P(region r contains any RNA fragments)            "is_expressed"
c_r = P(region r is enriched by the capture assay)      "is_captured"
```

Plus two derived per-region scalars consumed by downstream EM:

```text
mu_gdna(r)   = E[ D_r | data ]                          (numerator)
A_r          = 1 + c_r * (gamma_r - 1)                  (denominator scalar)
```

And two scalars for the whole library:

```text
rho_off       off-target gDNA fragment rate (frag/bp)
kappa_d       gDNA strand-balance concentration (stranded libraries only)
```

The `400` precision cap and the `rho_post / rho_ref` exposure pathway are gone. `A_r` is derived from the latent capture model, not from a ratio of density posteriors.

---

## 1. What Has To Go

Same as v2 §1 — unchanged. The cleanup is captured in PR 0. Re-stating the deletion list for completeness only:

`DENSITY_PRIOR_MIN_CV`, `DENSITY_PRIOR_MAX_PRECISION`, `compute_beta_cap`, `GammaRatePrior.beta_cap` / `cap_applied`, `DensityEvidence.relative_exposure`, `RegionExposure.from_density(max_exposure=)`, `_fit_anchor_priors` family `"ALL"`, `select_rho_ref` `"WEIGHTED_FAMILIES"`, `density_max_exposure` CLI flag, `fuse_density_and_strand` (rewritten, not kept), and the duplicate `boundary_left_leff`/`boundary_right_leff` fields (collapsed to one).

`bp_weighted_mean_exposure_over_blocks` is **kept as the v3 denominator** (resolution to Q1: fine-region granularity is high enough that bp-weighted mean is a fair first approximation; FL-convolved integral is a follow-up PR if and only if a real benchmark demands it).

---

## 2. The Latent Model

### 2.1 Variables (unchanged from v2)

```text
e_r       in [0, 1]    expression probability                per region
c_r       in [0, 1]    capture probability                   per region
gamma_r   >= 1         capture enrichment when captured      per region
rho_off                off-target gDNA rate (frag/bp)        scalar
kappa_d                gDNA strand-balance concentration     scalar (stranded only)
```

### 2.2 Generative model (unchanged from v2 §2.2)

```text
lambda_gdna(r) = rho_off * (1 + c_r * (gamma_r - 1))         = rho_off * A_r

contained_unspliced_total(r) ~ Poisson( lambda_gdna(r) * L_c(r) + nu_unspliced(r) * L_c(r) )
boundary_unspliced_total(r)  ~ Poisson( lambda_gdna(r) * L_b(r) + nu_unspliced(r) * L_b(r) )
contained_spliced_total(r)   ~ Poisson( nu_spliced(r)  * L_c(r) )
boundary_spliced_total(r)    ~ Poisson( nu_spliced(r)  * L_b(r) )

strand split of unspliced    ~ BetaBinomial( kappa_d ) for gDNA
                               Binomial( p_r1_sense )  for RNA            (stranded only)
```

`nu_unspliced(r)` and `nu_spliced(r)` are nuisance RNA rates we do not estimate. Their evidence enters through three channels:

- `boundary_spliced_total > 0` or `contained_spliced_total > 0` ⇒ `e_r := 1`.
- strand-deconvolved RNA-lower-bound > 0 ⇒ `e_r := 1` (stranded only).
- otherwise `e_r` floats and is updated by the joint posterior.

Boundary spliced mass must be subtracted from boundary unspliced mass when fitting `gamma_r`; see §3.2.

### 2.3 Three seed pools

This is where v3 diverges from v2. v2 had one off-target pool `S_off`. v3 has three pools, each with a different role and a different evidence type. The names are chosen to match the user's TODO §"4 groups based on expressed × targeted".

```text
Pool A  --  off-target background        (not expressed AND not captured)
            seeds rho_off
            intergenic + intron-only anchors
            AND no spliced mass
            AND no strand-RNA evidence (stranded only)
            AND not in top T% by contained density (top-T clip)

Pool B  --  on-target gDNA control       (not expressed AND captured)         *** key insight from Q6 ***
            seeds gamma_r and validates capture detection
            *** only available with strand-specific data ***
            exonic regions consistent with pure gDNA at confidence c
              (kappa_d strand-balance test passes, no spliced mass)
            AND boundary or contained gDNA density > rho_off
            screen_no_rna_exons() already does most of this work

Pool C  --  on-target expressed          (expressed AND captured)
            *** never used to fit any parameter ***
            it is the population we predict for, not from.
            Diagnostics only: median A_r, gDNA leak rate.
```

The fourth combination (expressed AND off-target) is not a pool. It is a residual category and is the hardest to model. Its members are simply regions that never qualified for A, B, or C.

**Why Pool B matters.** Without strand evidence, the only way to detect capture is boundary flux enrichment, which is sparse and zero-inflated (Q4). With strand evidence, we can identify exonic regions whose strand signature is consistent with pure gDNA (i.e., not expressed in this sample) but whose contained density is elevated above `rho_off`. These regions reveal `gamma_r` with no RNA confounding. They are the gold standard for capture-enrichment estimation and they exist precisely because hybrid capture panels carry many baits over genes that are not expressed in the user's sample.

This pool is what makes self-sufficient (no-BED) capture detection actually work on stranded libraries. The existing `screen_no_rna_exons` in [strand_deconv.py](src/rigel/calibration/strand_deconv.py) line 563 already implements the "exon consistent with pure gDNA" screen; v3 just consumes its output instead of using it only to refit `kappa_d`.

### 2.4 Initialization

```text
0.  build region table, ledger, FL models                             (existing)
1.  if stranded:
       fit strand_summary, strand_counts, kappa_d                     (existing)
       compute per-region strand-deconvolved gDNA/RNA estimates       (existing)
    else:
       all strand-derived quantities are None
2.  Pool A := {intergenic ∪ intron-only}
              ∩ {no spliced mass}
              ∩ {no strand-RNA evidence}                              (no-op when unstranded)
3.  trim top T% of Pool A by contained density                        (default T = 0.01)
4.  rho_off := sum(contained[A]) / sum(L_c[A])                        (fragments per bp)
5.  if stranded:
       Pool B := screen_no_rna_exons(...) refined to {density > rho_off}
    else:
       Pool B := empty
6.  for each region, evaluate Poisson tail of boundary-unspliced-minus-spliced
      flux under rho_off * L_b(r); call this the boundary enrichment test.
7.  seed c_r:
       Pool A members:          c_r := 0
       Pool B members:          c_r := 1
       boundary test rejects:   c_r := 1
       all others:              c_r := 0.5 (unknown)
8.  seed gamma_r:
       Pool B:                  gamma_r := contained_density(r) / rho_off
       boundary-test rejects:   gamma_r := boundary_density(r) / rho_off    (after spliced subtraction)
       else:                    gamma_r := 1
9.  seed e_r:
       any spliced mass:        e_r := 1
       strand-RNA-lower > 0:    e_r := 1   (stranded only)
       Pool A:                  e_r := 0
       Pool B:                  e_r := 0
       else:                    e_r := 0.5
```

There is **one** tunable threshold: `alpha_capture` (default 0.01) for the boundary Poisson test. `trim_top_fraction` (default 0.01) is the other configurable knob.

### 2.5 Iteration

v2 fixed two passes. v3 iterates until convergence, capped at `max_passes` (default 5):

```text
repeat:
  prev_rho_off  := rho_off
  prev_A_r      := A_r

  Pool A := refresh from current (e_r, c_r): require c_r < c_thresh AND e_r < e_thresh
            AND (signature is intergenic | intron-only)
            AND (no spliced mass)
            AND (strand-RNA-lower == 0 if stranded)
            AND not in current top-T% by contained density
  rho_off := sum(contained[A]) / sum(L_c[A])

  Pool B := refresh from screen_no_rna_exons() and {contained_density > rho_off}
            (stranded only)
  refit gamma_r on Pool B (preferred) and boundary-test rejects (fallback)

  re-evaluate joint posterior p(e_r, c_r, D_r | counts, rho_off, gamma_r, kappa_d)
    using the exact/approx split inherited from fuse_density_and_strand

  A_r := 1 + c_r * (gamma_r - 1)
  converged if   |rho_off - prev_rho_off| / max(prev_rho_off, eps) < tol_rho
            AND  max_r |A_r - prev_A_r|                            < tol_A
  break if converged
until max_passes
```

Defaults: `c_thresh = 0.5`, `e_thresh = 0.5`, `tol_rho = 0.01`, `tol_A = 0.05`, `max_passes = 5`.

Resolution to Q3 / Q10: iteration is the mechanism for purifying Pool A. Each pass can eject a region from Pool A (a high-nRNA intron, an unannotated transcript, a negative-control probe in intergenic space) once enough evidence has accumulated. The top-T clip is the robustness fallback for hostile data where iteration alone fails. In synthetic tests we expect 1–3 passes to suffice; the loop must report `n_passes` and `converged` in `summary.json`.

Optional refinement: refit `kappa_d` once at the top of each pass using the *new* Pool A ∪ Pool B as the seed set. The existing `estimate_kappa_d` already supports this seed-plus-self-training shape.

### 2.6 Derived outputs (unchanged from v2 §2.5)

```text
mu_gdna(r)    = E[ D_r | counts, posteriors ]
upper_gdna(r) = c-quantile of D_r posterior
A_r           = 1 + c_r * (gamma_r - 1)
rna_lower(r)  = max(observed_compatible(r) - upper_gdna(r), 0)
```

---

## 3. Implementation Details Worth Pinning

### 3.1 Boundary signal must be spliced-subtracted

Boundary flux is the workhorse for `gamma_r` estimation when Pool B is unavailable or sparse. But boundary flux contains both intron→exon (gDNA / nRNA) and exon→exon (mature spliced RNA across junctions) fragments. The 12-channel ledger already separates `boundary_unspliced` from `boundary_spliced`. v3 uses **only** the unspliced channel for the boundary enrichment test and for `gamma_r` estimation:

```text
B_gdna_candidate(r) = boundary_left_unspliced(r) + boundary_right_unspliced(r)
B_under_null(r)     ~ Poisson( rho_off * L_b(r) )
```

This is correct because `lambda_gdna * L_b` is the only term contributing to the unspliced boundary channel for an unexpressed-but-captured exon. For an expressed-AND-captured exon (Pool C), the unspliced boundary count *also* contains nRNA flux through introns, which biases `gamma_r` upward; that is acceptable because we never fit `gamma_r` from Pool C anyway. We only ever fit `gamma_r` from Pool B (no expression) and from boundary-test rejects that are not in Pool A.

### 3.2 Region geometry caveat: what is a "boundary"?

The boundary compartments in [signature.py](src/rigel/calibration/signature.py) (`COMPARTMENT_BOUNDARY_LEFT`/`RIGHT`) tag fragments that *cross* the receiving region's edges. For an exonic region this is intron→exon flux on the left side and exon→intron flux on the right side — exactly the gDNA-rich edge signal we want. For an intronic region this is exon→intron flux. v3 needs no change here, but the documentation in `density_observation.py` should be updated to make this explicit.

### 3.3 `c_r` posterior softness — start hard, soften if needed (resolution to Q4)

PR 2 ships the boundary enrichment test as a hard indicator first:

```text
c_seed(r) = 1[ p_value(r) < alpha_capture ]
```

The joint posterior in PR 3 returns continuous `c_r ∈ [0, 1]` regardless of seed hardness. The seed is only an initialization. If sparse-boundary instability appears in benchmarks, swap the hard indicator for a sigmoid on `-log p_value` in a follow-up; the seam is intentional.

### 3.4 FL is observed (resolution to Q5)

The gDNA FL model already trains on observed fragments, which under capture is the capture-conditioned distribution. v3 makes no attempt to recover a pre-capture latent FL. **Truth for benchmarking is the post-capture library**, not the simulator's pre-capture parameters. Sim infrastructure (`scripts/sim/`) must surface a `truth_post_capture.tsv` that reflects the realized abundances and realized FL distribution after the capture step. The pre-capture-truth comparison reporters should be removed from the synthetic capture benchmarks.

This is a benchmarking change, not a calibration change, but it is a prerequisite for honest acceptance gates.

### 3.5 Ordering of strand and density passes (resolution to Q6)

The v2 orchestrator runs density before strand. v3 reorders so strand runs first whenever the library is informatively stranded:

```text
fl_models
strand_summary
if stranded:
   strand_counts
   kappa_d (seed)
   region_gdna (strand-deconvolved RNA-lower)
fit_background_model       (consumes strand_rna_lower if available)
fit_capture_model          (consumes Pool B from screen_no_rna_exons if available)
latent_states              (iteration loop)
region_exposure
```

For unstranded libraries the strand block is a no-op and `Pool B = empty`. The downstream code paths must handle this gracefully and `summary.json` must report `pool_b_available: false`.

### 3.6 Per-region scalar A_r — no double counting (resolution to Q9)

`A_r` is an intensive quantity (a multiplier on a rate). Two overlapping multi-loci that share a region each multiply their own `gdna_eff_len` by that region's `A_r`. This is correct: each multi-locus has its own physical opportunity surface, and a shared region contributes its `A_r` to both surfaces. The EM never sees the same *fragment* twice; the geometry allocator already prevents that.

### 3.7 Capture state diagnostics on `quant.feather` (resolution to Q7)

No RNA exposure model. Diagnostics only:

- `quant.feather` columns: `mean_A_r`, `n_regions_captured`, `n_regions_total`.
- `gene_quant.feather`: same columns at gene level.
- `loci.feather`: `n_regions_captured`, `mean_A_r`, `n_regions_expressed_post`.
- `summary.json`: `capture_detection` block with `mode ∈ {uniform, suspected_capture, capture}`, median `gamma_r` over captured regions, `n_captured_seeded`, `p99_A_r`, `pool_b_available`, `n_passes`, `converged`.

The user is responsible for interpreting low-`mean_A_r` transcripts as off-target. Rigel does not silently translate low capture into low TPM and does not silently report off-target transcripts as zero expression — but it does not "rescue" them either. They are simply not measurable in this library by design.

### 3.8 Unstranded capture (resolution to Q8)

Unstranded capture is the hardest regime. v3 implements it on the same code paths but with `Pool B = empty` and no strand likelihood term in the joint posterior. The only signals are spliced counts (for `e_r`) and unspliced boundary flux (for `c_r`, `gamma_r`). The synthetic suite must include cells `gdna_high_ss_0.50_capture_on` and `gdna_none_ss_0.50_capture_on`. Acceptance is "best effort" — not the same bar as stranded — and the `summary.json` must include an honest `confidence: low | medium | high` label keyed off `pool_b_available` and the count of boundary-test rejects.

---

## 4. Module Layout (delta from v2 §3)

Same module layout as v2 plus one explicit responsibility shift:

```
src/rigel/calibration/
  background.py        NEW  — Pool A and rho_off (iterative)
  capture_model.py     NEW  — Pool B (from screen_no_rna_exons) + boundary enrichment test
                              + gamma_r fit (iterative)
  latent_states.py     NEW  — joint posterior p(e_r, c_r, D_r) + the iteration loop
  exposure.py          rewrite — A_r from LatentStates
  integration.py       rewrite — joint posterior cell evaluators
  density_model.py     gut — keep only the Pearson-trimmed Gamma fit utility for diagnostics
  prior.py             minor — read mu_gdna / upper_gdna / A_r from LatentStates
  _orchestrator.py     rewrite — new ordering, calls iteration driver
```

Strand modules (`strand_deconv.py`, `strand_balance.py`, `strand_summary.py`) are unchanged. The C++ scanner, signature/ledger/payload, FL models, and region geometry are all unchanged.

---

## 5. PRs (delta from v2 §4)

The PR sequencing is the same as v2 with three substantive changes:

| PR | v2 description | v3 change |
|---|---|---|
| 0 | Cleanup | unchanged |
| 1 | `background.py` | now exposes a `refresh(model, latent_states) -> BackgroundModel` method for the iteration loop |
| 2 | `capture_model.py` | now consumes a `Pool B` mask from strand deconvolution as an additional `gamma_r` training source; boundary signal must be spliced-subtracted; `refresh(model, latent_states) -> CaptureModel` for iteration |
| 3 | `latent_states.py` | replace fixed two-pass with convergence-gated loop; report `n_passes` and `converged` |
| 4 | orchestrator + exposure | re-order strand-before-density |
| 5 | prior wiring | unchanged |
| 6 | diagnostics + CLI | add `--alpha-capture`, `--background-trim-top-fraction`, `--max-calibration-passes`, `--tol-rho-off`, `--tol-a-r`; add transcript-level `mean_A_r` and `n_regions_captured` |
| 7 | tests + synthetic gates | benchmark against post-capture truth, not pre-capture (Q5) |
| 8 (deferred) | FL-convolved exposure integral | only if a real benchmark demands it (Q1) |

PR 6 also adds the `pool_b_available` and `confidence` flags to `summary.json`.

---

## 6. Acceptance Summary

Unchanged structure from v2, plus the iteration acceptance and the post-capture-truth correction:

| Condition | Acceptance |
|---|---|
| capture-off, gdna=none | byte-identical numerics to today (within float tolerance) |
| capture-off, gdna=high | byte-identical numerics to today |
| capture-on, gdna=none | `mu_gdna ≈ 0` everywhere; no invented gDNA |
| capture-on, gdna=high, ss=0.99 | gDNA→RNA leak drops ≥5× vs today; `n_passes ≤ 3` in 95% of runs |
| capture-on, gdna=high, ss=0.50 | gDNA→RNA leak drops ≥2× vs today |
| Pool B (stranded captured-not-expressed) | median `gamma_r` from Pool B is within 30% of median `gamma_r` from boundary test on overlapping exons (cross-validation) |
| Iteration convergence | `converged == True` in ≥95% of synthetic conditions within `max_passes = 5` |

Structural:

- No precision floors. No CV floors. No `400`.
- All tunables (`alpha_capture`, `trim_top_fraction`, `max_passes`, `tol_rho`, `tol_A`, `c_thresh`, `e_thresh`) statistically interpretable.
- Benchmarks compare against post-capture truth, never pre-capture sim parameters.

---

## 7. Residual Issues Found While Sharpening v3

The following are real concerns that were not in v2 and are not yet fully resolved by v3. None is a blocker. Each has a recommendation.

### R1. Boundary spliced contamination of `gamma_r` for Pool C

In an expressed-AND-captured exon, the unspliced boundary count includes intron-retained nRNA boundary flux. If a Pool B exon is mis-screened (strand evidence falsely says "no RNA" while there is in fact nRNA), `gamma_r` inflates. The strand RNA-lower test screens against mature RNA balance, not nRNA. **Recommendation**: track the ratio `contained_spliced(r) / contained_unspliced(r)` per Pool B region; reject from Pool B if it crosses a threshold (default 0.1). This is a one-line filter and a one-knob addition.

### R2. `kappa_d` refit after Pool A purification

`kappa_d` is currently fit before the iteration loop. After iteration purifies Pool A (e.g., ejects high-nRNA introns and unannotated tx loci), the seed set for `kappa_d` is *better*. **Recommendation**: refit `kappa_d` once at the top of each iteration pass using the current Pool A ∪ Pool B. This is cheap and the strand-deconv machinery already supports a `seed_mask` argument.

### R3. Identifiability collapse: unstranded + small locus + all captured + all expressed

If a locus has no anchor regions, no strand information, and all of its exons are both expressed and captured, the model has no way to separate `gamma_r` from RNA expression — the boundary test fires but cannot tell us how much of the elevated unspliced count is nRNA versus gDNA. **Recommendation**: when a locus meets this signature, set `enable_gdna = False` for that locus and emit a `gdna_unidentifiable` flag in `loci.feather`. This is conservative and avoids inventing gDNA where we cannot measure it.

### R4. nRNA in introns under strand-specific data

Under stranded data, intron-only regions with non-zero strand-RNA-lower are correctly identified as carrying RNA (nascent transcription) and are excluded from Pool A. The library does contain those fragments and they are physically gDNA-adjacent. Excluding them does not bias `rho_off` (they are not gDNA-only), but it does shrink Pool A in nRNA-rich libraries. **Recommendation**: report `pool_a_size`, `pool_a_fragments`, `pool_a_eff_length` in `summary.json` so we can see when Pool A is starved. If `pool_a_fragments < 1000` or `pool_a_eff_length < 1 Mb`, switch `fit_status` to `"sparse"` and broaden `rho_off` uncertainty downstream.

### R5. Pool B sample size under aggressive capture

In a panel with very narrow gene targeting and a sample that expresses most targeted genes, Pool B (captured-and-not-expressed exons) can be small. **Recommendation**: track `n_pool_b` and emit a `low_pool_b` confidence flag. The boundary-test fallback for `gamma_r` is the bypass when Pool B is too small.

### R6. `gamma_r` is a noisy per-region estimate

A single exon's `gamma_r` is one ratio of two Poisson counts and has large variance for short exons. v3 reports the per-region value, but per-region `A_r` propagates that noise into per-locus denominators. **Recommendation**: in PR 5, when computing the per-locus bp-weighted mean of `A_r`, weight by bp *and* by `L_c(r)` (effective length) — short noisy regions contribute less. This is a one-line change in `bp_weighted_mean_exposure_over_blocks`.

### R7. Boundary effective length `L_b(r)` is roughly the mean FL

`fractional_boundary_side_exposure` in [_exposure.py](src/rigel/calibration/_exposure.py) yields an FL-aware boundary opportunity that is bounded above by the receiving region's length. For exons shorter than the mean FL (~150–170 bp), boundary opportunity saturates and the boundary enrichment test loses power. **Recommendation**: in the boundary enrichment test, downweight regions where `L_b(r) < min_boundary_opportunity` (default 1.0 bp, already a threshold in `density_model`); these regions cannot reject the null and must abstain (`p_value = NaN`, `c_seed = 0`, gamma sourced from Pool B if available).

### R8. Iteration may oscillate at the `c_thresh = 0.5` boundary

A region with posterior `c_r` hovering near 0.5 could flip in and out of Pool A across iterations, preventing convergence. **Recommendation**: apply a hysteresis band. A region is admitted to Pool A only when `c_r < 0.4` and ejected only when `c_r > 0.6`. Same for `e_r`.

### R9. Top-T clip vs strand-RNA filter interaction

The top-T contained-density clip is currently described as operating on the initial Pool A. Once strand evidence has already filtered nRNA-bright introns, the top-T clip is mostly removing real gDNA hotspots and unannotated transcripts. **Recommendation**: keep both filters but make `T` configurable (default 0.01 with strand, 0.05 without). The CLI exposes `--background-trim-top-fraction`.

### R10. CLI knobs surface area

PR 6 introduces seven new knobs (`alpha_capture`, `trim_top_fraction`, `max_passes`, `tol_rho`, `tol_A`, `c_thresh`/`e_thresh` hysteresis). This is a lot. **Recommendation**: bundle into a `CalibrationConfig.calibration_v3` dataclass with sane defaults and one CLI flag `--calibration-profile {default, sensitive, conservative}` that picks presets. The individual knobs remain available for advanced users via `--calibration-overrides key=value,...`.

### R11. `enable_gdna` gating in PR 5 needs a sharper criterion

v2 PR 5 said: "enable gDNA when at least one region has `c_r * gamma_r > 1 + epsilon` OR contains anchor regions with `mu_gdna > 0`." This is fine, but on capture-off/no-gDNA, `mu_gdna` is dominated by Poisson noise on anchors, and we may falsely enable gDNA. **Recommendation**: gate on `sum(mu_gdna) > threshold` where `threshold = 0.5 * sqrt(sum(L_c) * rho_off)` (signal-to-noise floor at one σ). This is the same gate currently informally controlled by `_GDNA_ENABLE_UPPER_EPS = 1e-6` in `prior.py`.

### R12. Strand-deconv input depends on the calibration scan, not the EM data

`build_strand_region_counts` reads from the calibration scan payload, which is a sample of fragments. Pool B sizing therefore depends on the scan sample size. For very small samples, Pool B may be falsely small. **Recommendation**: document and surface `n_calibration_fragments_used_for_pool_b` in `summary.json`. No code change required.

### R13. Pre-capture sim truth files

The benchmarking framework currently emits `truth_abundances_*.tsv` keyed to the pre-capture simulator state. v3 acceptance gates against the post-capture realized state. **Recommendation**: extend `scripts/sim/locus_sweep.py` and the larger benchmarking framework to emit `truth_post_capture.tsv` with the realized FL distribution and realized per-transcript abundance after the capture sampler runs. This is a benchmarking PR, not a calibration PR; sequence it before PR 7 if possible.

### R14. The `rho_off == 0` deterministic-zero path

When Pool A is starved (no-gDNA libraries, fully-captured panels with no intergenic depth), `rho_off → 0`. The boundary Poisson test then degenerates: `lambda_null = 0` and every nonzero boundary count rejects the null with p = 0. We must guard against this. **Recommendation**: when `rho_off < eps_rho_off` (e.g. 1e-6 frag/bp), skip the boundary test entirely, set `c_r = 0` everywhere, `gamma_r = 1` everywhere, and emit `capture_detection.mode = "uniform"`. The library is treated as gDNA-free and the EM does not enable gDNA at any locus. This is also the safe behavior for non-capture libraries.

---

## 8. Non-Goals (unchanged from v2 §7)

- v3 does **not** recover all-transcript TPM under capture.
- v3 does **not** implement per-fragment `log W(f)` (resolution to Q2: locus-level denominator is sufficient; the older "per-fragment exposure correction" proposal is superseded).
- v3 does **not** introduce an RNA exposure model.
- v3 does **not** ship the FL-convolved exposure integral (deferred to PR 8 contingent on benchmark need).
- v3 does **not** consume a probe BED; the `--capture-bed` seam is reserved for a future opt-in path and is intentionally not shipped in v3.
