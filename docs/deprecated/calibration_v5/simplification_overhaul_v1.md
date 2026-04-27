# Plan: Replace v5 Calibration with Simple Regional Deconvolution (SRD)

## Context

Rigel's **calibration step** (`src/rigel/calibration/`) is failing on hybrid-capture data. Its two outputs are:

1. A **gDNA fragment-length model** (`gdna_fl_model`), consumed by the C++ per-locus EM scorer to tell gDNA apart from RNA by length. **This is the primary product.**
2. A **per-region gDNA estimate** (`region_e_gdna[i]`, `region_n_total[i]`), currently used by `locus.compute_locus_priors()` to seed per-locus Dirichlet warm-starts. Demoted to diagnostic / opt-in under SRD.

The current v5 implementation (~1000 LOC in `_em.py` + `_stats.py` + `_fl_model.py`) is a three-channel mixture EM over regions (density + strand + FL) with Beta-Binomial dispersion fitting and Gauss–Hermite / Gauss–Legendre quadrature. **It fails because the density channel assumes a single global gDNA rate `λ_G`, but hybrid capture creates 100–140× spatial variation.** Density LLRs (±20 nats) overwhelm the correct strand LLRs (~0.23·n nats). Patching with Negative Binomial dispersion kept adding parameters and fragility.

Design goals for the replacement:

- **Library-agnostic.** Must work well on both stranded and unstranded data — unstranded RNA-seq is too large a population to regress on.
- **Simple and blazing fast.** Closed-form per-fragment classification; at most one 1-D mixture fit. No region-level EM, no quadrature, no dispersion estimation.
- **No density comparison across regions.** Immune to hybrid-capture spatial bias.
- **No distributional assumptions on RNA expression.** Histograms only.
- **Few hyperparameters.** Nothing per-protocol.
- **Principled fragment pools.** Derived from genomic-geometry first principles (exon fit, splice compatibility, strand), not heuristics.

## Design: Simple Regional Deconvolution (SRD)

**Core shift**: stop classifying *regions* with density. Classify individual *fragments* using cheap geometric / strand signals, then aggregate per region for diagnostics.

### Signals used

| Signal | Strength | Library |
|---|---|---|
| Spliced CIGAR (N-junction) | **Hard** — 100% RNA | All |
| Exon-fit geometry (unspliced fragment must fit inside one annotated exon to be mature mRNA) | **Strong** — excludes mRNA | All |
| Strand vs tx_strand (antisense in confident tx_strand region) | **Strong** — purity ≈ (2·SS − 1) | Stranded only |
| Fragment length | **Soft** — posterior weight | All (after FL training) |

**Explicitly dropped**: regional density comparisons, Poisson/NegBin rate estimation, LogNormal expression priors, Beta-Binomial dispersion.

### Pass 0 — Fragment categorization (single walk over buffer, O(n_frags))

Walk `FragmentBuffer` once and assign each fragment exactly one category using cheap geometric checks against the annotation index.

The geometric test "does an unspliced fragment of length FL fit inside an annotated exon of length L?" has a clean definition: the fragment's genomic span `[g0, g1]` must satisfy `g0 ≥ exon_start ∧ g1 ≤ exon_end` for some exon. The per-exon **unspliced landing zone** is `max(0, L − FL + 1)` bases; when this is 0 for every overlapping exon, the fragment cannot be mature mRNA.

Categories (names per your feedback):

1. **SPLICED** — has N in CIGAR (annotated or unannotated splice). Already flagged by the BAM scanner.
2. **UNSPLICED_SENSE_EXONIC** — unspliced, fits entirely inside one annotated exon, strand matches tx_strand. (Stranded libraries only; else falls into next category.)
3. **UNSPLICED_ANTISENSE_EXONIC** — unspliced, fits entirely inside one annotated exon, strand opposite to tx_strand. **High-confidence gDNA for stranded libraries** (purity ≈ 2·SS − 1).
4. **UNSPLICED_EXONIC_AMBIG** — unspliced, fits inside an exon, but tx_strand is unknown or library is unstranded.
5. **EXON_INCOMPATIBLE** — unspliced, but does **not** fit inside any annotated exon: either overhangs an exon boundary into an intron (with a splice site the fragment "should have" used), or is longer than any containing exon. **Cannot be mature mRNA.** This is the new first-principles signal.
6. **NONEXONIC** — unspliced with no annotated exon overlap. Split for diagnostics into:
   - `INTRONIC` — wholly inside an annotated intron.
   - `INTERGENIC` — no transcript overlap.

Notes:

- `EXON_INCOMPATIBLE` and `NONEXONIC` are each a **gDNA + nRNA mix**, not pure gDNA. We handle that in Pass 1 by assuming `nRNA_FL ≈ mRNA_FL` (both reflect post-library-prep fragmentation of long precursor molecules), so the mixture decomposition still recovers a clean `gDNA_FL`.
- We reuse the exon/intron interval index already built by `index.py`. No new data structures.
- A single-pass numpy implementation using the existing `FragmentBuffer` columns plus an interval-tree lookup should easily stay under 1s per million fragments.

### Pass 1 — gDNA_FL training (one 1-D fit, all libraries)

The goal: recover `gDNA_FL` as a length histogram. Different libraries use different source pools, but the same 1-D mixture decomposition.

**Source pool selection:**

- **Stranded libraries (SS ≥ 0.6)** — use `UNSPLICED_ANTISENSE_EXONIC` directly. Purity is 2·SS − 1 (≥ 0.2 for SS ≥ 0.6, typically ≥ 0.8 for dUTP protocols). A single histogram over these fragment lengths is a near-clean gDNA_FL. Optional light cleanup: subtract the SS-proportional mRNA leakage using the known RNA_FL.
- **Unstranded libraries (SS < 0.6)** — use `EXON_INCOMPATIBLE ∪ NONEXONIC` as the "non-mRNA pool". This pool is gDNA + nRNA. Under the assumption `nRNA_FL ≈ RNA_FL`, fit the 1-D mixture:
  ```
  pool_FL(l)  =  π · gDNA_FL(l)  +  (1 − π) · RNA_FL(l)
  ```
  with `RNA_FL` fixed (from spliced reads) and `gDNA_FL` free as a Dirichlet-smoothed histogram. Recovers `gDNA_FL` and the pool-level gDNA fraction `π` in <30 iterations. Well-posed because `RNA_FL` is known; no quadrature.

- **Stranded libraries** can optionally also run the 1-D mixture on their pool for cross-validation, but the direct histogram is the default (simpler and faster).

**Failure mode handling**: if either source pool is too small (< a few hundred fragments) or the 1-D mixture fails to separate (`π → 0` or `gDNA_FL ≈ RNA_FL`), emit a diagnostic warning and fall back to a deliberately-broad default `gDNA_FL` (e.g. geometric over the observed length range). The per-locus EM still runs — with weaker FL discrimination.

### Pass 2 — Per-fragment classification (closed form)

For each fragment `f` in the buffer, compute `γ_f = P(gDNA | features(f))` with no iteration:

```
SPLICED                         → γ_f = 0                                       (hard RNA)
UNSPLICED_ANTISENSE_EXONIC      → γ_f = 2·SS − 1                                (strand purity)
UNSPLICED_SENSE_EXONIC          → γ_f from FL-LLR + strand-LLR                  (sigmoid)
UNSPLICED_EXONIC_AMBIG          → γ_f from FL-LLR                               (sigmoid)
EXON_INCOMPATIBLE / NONEXONIC   → γ_f from FL-LLR with prior shifted toward gDNA
                                  (starting log-prior = log(π_pool/(1−π_pool))
                                  from the Pass-1 mixture)
```

The FL log-likelihood-ratio is a single array lookup: `log gDNA_FL(l) − log RNA_FL(l)`. Strand LLR reuses the existing `strand_model.strand_likelihood`. No region-level term — that's the whole point.

### Pass 3 — Aggregate to regions (diagnostic / optional prior)

```
region_e_gdna[i]  = Σ_{f ∈ region i} γ_f
region_n_total[i] = number of fragments in region i
```

These are **reported** in `CalibrationResult` for diagnostics but are not load-bearing. See the "priors question" below.

### The priors question: capitalism vs warm-start

Your direct question: *"how important is it to set the prior for the EM? What if we ran the EM with no prior, letting fragments compete freely?"*

My answer:

- With SRD, the per-fragment signals available to the per-locus EM are **strong on their own**: splice hard-anchor + FL discrimination (via the newly-reliable `gDNA_FL`) + strand (stranded libs). The region prior mostly serves as a warm-start, which matters most when per-locus fragment counts are low.
- v5's region prior was load-bearing because density was doing most of the work; under SRD it is not.
- **Default recommendation: capitalism.** Set `c_base = 0` by default in `compute_locus_priors`, producing flat Dirichlet(0, 0, …) priors and letting the C++ EM run the per-locus mixture without a warm-start bias. `region_e_gdna` / `region_n_total` stay in the result for reporting, but aren't consumed.
- **Keep the machinery** behind a single config flag (`EMConfig.calibration_prior_strength: float = 0.0`). Set it to a small value (1.0) for users who want to experiment with informative warm-starts. This is one flag, easily removable later if capitalism proves universally better.

This keeps the simplification honest: the pipeline wiring stays symmetric, nothing is deleted that we might want back, but the default behavior reflects the observation that regional estimates are diagnostic, not causal.

### Density for unstranded libraries — explicitly deferred

You raised the possibility that density over `NONEXONIC` could add signal for unstranded libs. I agree in principle but am keeping it **out of scope for v1**. Rationale:

- The exon-geometry pool (`EXON_INCOMPATIBLE ∪ NONEXONIC`) already provides library-agnostic gDNA_FL training without density.
- Adding density reintroduces the exact failure mode we're trying to escape (spatial-variation misspecification under hybrid capture).
- If the v1 diagnostics show unstranded-library gDNA_FL quality is insufficient, we can add density as a **tertiary signal gated on a coverage-uniformity check** in v2. Ship simple first; measure; add only if needed.

The categorization already produces `INTRONIC` and `INTERGENIC` counts separately, so the hook is there if we want it.

## Why this works where v5 failed

| v5 failure mode | SRD behavior |
|---|---|
| Single global `λ_G` misspecified under hybrid capture (100× spatial variation) | **No density channel.** Spatial variation irrelevant. |
| Poisson / LogNormal / Beta-Binomial / κ dispersion fitting | None. Histograms + Bernoulli posteriors only. |
| 50-iteration fused-LLR EM with quadrature | Closed-form per fragment. Only the 1-D FL-mixture EM iterates (well-posed, <30 iters). |
| Brittle hyperparameters (κ bracket, quadrature nodes, eligibility floor) | One effective knob: pool-size warning threshold. |
| Regional γ derived from one global number (`λ_G`) | Regional γ is the actual sum of per-fragment posteriors; genuinely regional. |
| No signal for unstranded libraries without density | **Exon-incompatible pool**: library-agnostic first-principles gDNA source. |
| Region prior load-bearing but misfit | Region prior **optional and off by default** (capitalism). |

## Files to change

### New

- `src/rigel/calibration/_simple.py` — SRD implementation. Exports `calibrate_gdna(...)` returning `CalibrationResult`. Target ~300–400 LOC including docstrings. Sections: Pass 0 categorizer, Pass 1 FL trainer (branch on SS), Pass 2 per-fragment classifier, Pass 3 aggregator.
- `src/rigel/calibration/_fl_mixture.py` — 1-D FL mixture EM for unstranded libraries (~80 LOC). `fit_fl_mixture(pool_fl_hist, rna_fl_model, max_iter=50, tol=1e-4) -> (pi, gdna_fl_model, converged)`.
- `src/rigel/calibration/_categorize.py` — Pass-0 fragment categorizer using the annotation index. Returns a per-fragment category array (uint8 codes). This is the only file that touches the exon/intron interval indexes.
- `tests/test_calibration_simple.py` — unit tests (see Verification).
- `tests/test_fl_mixture.py` — focused tests for the 1-D mixture fit.
- `tests/test_categorize.py` — focused tests for exon-fit geometry.

### Modified

- `src/rigel/calibration/__init__.py` — export the new `calibrate_gdna` from `_simple`.
- `src/rigel/calibration/_result.py` — add fields: `n_spliced`, `n_sense_exonic`, `n_antisense_exonic`, `n_exonic_ambig`, `n_exon_incompatible`, `n_intronic`, `n_intergenic`, `n_gdna_fl_pool`, `pi_pool`, `fl_mixture_converged`. Mark `mu_R`, `sigma_R`, `kappa_G`, `kappa_R`, `mixing_pi_soft` as deprecated (`Optional[float] = None`); they stay for one release so external consumers don't break.
- `src/rigel/config.py` — add `EMConfig.calibration_prior_strength: float = 0.0` (default: capitalism).
- `src/rigel/locus.py` — in `compute_locus_priors`, multiply all output `α` values by `calibration_prior_strength` passed through. When `0.0`, returns flat zeros (C++ EM treats as uninformative).
- `src/rigel/estimator.py` / `pipeline.py` — thread `calibration_prior_strength` from config. No other logic change.

### Parked (per your earlier answer — "park in a branch, fresh implementation on main")

On a branch `legacy/calibration-v5-em`, preserve for reference:
- `src/rigel/calibration/_em.py`
- `src/rigel/calibration/_stats.py`
- `src/rigel/calibration/_fl_model.py` (old)
- `src/rigel/calibration/_calibrate.py` (old orchestrator)
- `docs/calibration/*v5*.md`, `docs/calibration/channel_fusion*.md`, `docs/calibration/density_nb_model_plan.md`, etc.

Main branch loses these files. Git history preserves the full rollback path.

### Reused primitives (no changes)

- `FragmentBuffer` — columnar per-fragment store with `region_id`, `frag_len`, `frag_strand`, splice flags. Already populated by the BAM scanner.
- `frag_length_model.FragmentLengthModel` — histogram storage and likelihood lookup; correct for both RNA_FL and gDNA_FL.
- `strand_model.StrandModel` and `strand_likelihood()` — reused for the strand LLR term in Pass 2.
- `index.py` exon/intron interval index — used in `_categorize.py`.
- C++ `_scoring_impl` hard constraint on spliced fragments (line ~693: `!st.mm_is_any_spliced && st.mm_nh_gdna > 0`) — do not touch. This is what makes `γ = 0` for SPLICED fragments in the per-locus EM as well.
- `locus.compute_locus_priors` — interface preserved; behavior gated by `calibration_prior_strength`.

## Key file references (for implementation)

- Output contract to preserve: [`_result.py`](src/rigel/calibration/_result.py) `CalibrationResult`.
- Current orchestrator to replace: [`_calibrate.py:26`](src/rigel/calibration/_calibrate.py:26) `calibrate_gdna`.
- Downstream consumer (will be gated by `c_base`): [`locus.py:125`](src/rigel/locus.py:125) `compute_locus_priors`.
- Existing FL histogram primitive to reuse: [`frag_length_model.py`](src/rigel/frag_length_model.py).
- Existing strand LLR to reuse: [`strand_model.py`](src/rigel/strand_model.py) `strand_likelihood`.
- Exon/intron interval index source: [`index.py`](src/rigel/index.py).

## Verification

### Unit tests

1. **`test_categorize`** — synthetic exon/intron layouts, exercise all categories: a fragment fitting in one exon (SENSE/ANTISENSE/AMBIG), a fragment longer than its containing exon (EXON_INCOMPATIBLE), a fragment spanning exon-intron boundary (EXON_INCOMPATIBLE), purely intronic (INTRONIC), off-annotation (INTERGENIC), spliced (SPLICED). Assert exact categorization.
2. **`test_fl_mixture`** — inject known `gDNA_FL` (e.g. shifted log-normal) into a synthetic pool with known RNA_FL + nRNA_FL. Assert recovered `gDNA_FL` within KL tolerance and recovered `π` within ±0.02.
3. **`test_calibration_simple`**:
   - **Stranded pure-RNA**: SS=0.95, zero gDNA. Expect `region_e_gdna ≈ 0`, `n_antisense_exonic` small, `gdna_fl_model` flagged as low-confidence.
   - **Stranded 20% gDNA**: SS=0.95, injected 20% gDNA with distinct FL. Expect `π_global ≈ 0.20 ± 0.02`, `gDNA_FL` matches injection.
   - **Unstranded 20% gDNA**: SS=0.5, same mixture. Expect `π_pool` within ±0.05 of true pool fraction; `gDNA_FL` matches injection within tolerance.
   - **Hybrid-capture-like**: synthetic regions with 100× density variation but uniform 10% gDNA. v5 fails this; SRD should recover uniform 10% per region (not exact — will have Poisson scatter — but no 100× bias).
   - **Degenerate empty pool**: zero EXON_INCOMPATIBLE fragments and SS=0.5. Expect graceful fallback to broad default `gDNA_FL`, diagnostic warning emitted, no crash.
   - **Contract test**: `to_summary_dict()` returns all current top-level keys (plus new diagnostic keys). v5-era fields (`mu_R` etc.) present as `None` with a deprecation log on access if feasible.

### Integration

- `pytest tests/ -v` — full suite; regenerate goldens with `--update-golden` after eyeballing diffs. Splice-anchored loci should match v5 tightly; hybrid-capture regions should differ (and SRD should be closer to truth).
- Run any stranded VCaP titration scenario in `tests/scenarios_aligned/` if still present. v5's `quant.gdna_total / strand_inferred_gdna` was 1.1–1.5×; SRD target < 1.1×.
- Run at least one unstranded scenario to confirm the exon-incompatible pool produces a sensible `gDNA_FL`.

### Smoke test

```bash
pip install --no-build-isolation -e .
pytest tests/ -v -k calibration
pytest tests/ -v
rigel quant --bam <stranded_sample.bam> --index <idx>/ -o out/
rigel quant --bam <unstranded_sample.bam> --index <idx>/ -o out/
```

## Out of scope (deliberate)

- Density signal for `NONEXONIC` in unstranded libraries (deferred to v2 if diagnostics justify).
- Changes to the C++ per-locus EM or fragment scorer.
- Changes to strand model, RNA FL model, or locus-construction code.
- Removing the deprecated `CalibrationResult` fields (`mu_R`, `sigma_R`, `kappa_G`, `kappa_R`) — one release of overlap before deletion.
- Any "capture vs. non-capture" protocol detection (we make the method immune by not using cross-region density, rather than trying to detect capture).