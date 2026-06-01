# PR 6 — Integrate: calibration → locus prior → `quant_from_buffer`

This is the capstone. PR 0–5 built the calibrator in isolation; the pipeline
still stops the moment `calibrate()` returns (two `NotImplementedError`s, at
`pipeline.py:597` and `:733`). PR 6 closes the loop:

```
CalibrationResult  ──assemble_priors──▶  per-locus gDNA/RNA split prior
       +                                          │
   FragmentBuffer ──score──▶ loci ──partition──▶ per-locus EM ──▶ AbundanceEstimator ──▶ PipelineResult
```

End-to-end green returns at PR 6; the golden / smoke tests can run again.

---

## ✅ Prerequisite resolved — PR 5.5 (accumulator mass conservation)

A prerequisite bug surfaced while confirming the accumulator's boundary
attribution: a region a fragment **encompassed** (a fully-traversed interior
slice) had its mass deposited at full value on *both* boundaries, so a spanning
fragment deposited total mass 1.5 instead of 1.0 — breaking the exposure's
`∫overlap/ℓ\,dx = L` identity and over-stating gDNA/RNA mass for short exons
spanned by contiguous fragments.

**Fixed and merged in [PR 5.5](PR055_accumulator_conservation.md).** The deposit
is now per-slice: each slice contributes `slice_mass / n_cross` to each side it
crosses (`n_cross ∈ {1,2}`), so an interior slice splits 50/50 and end slices
keep full mass. Flux stays integer; downstream is unchanged. The accumulator now
conserves per-fragment mass (locked by `test_t6b`), so PR 6 can treat the
substrate's `M_cont + M_left + M_right` as a faithful per-region fragment-mass
total — the foundation the prior formulas (§I.2) rely on.

---

## 0. Problem & contract

The calibration emits **per-region** deconvoluted mass + exposure +
hyperparameters. The locus EM wants a **per-locus** prior. The only missing
pieces are (a) the *bridge* (`assemble_priors`) and (b) the *orchestration*
(`quant_from_buffer` body + `run_pipeline` tail).

**Survived the burn — reuse as-is** (verified present): `_setup_geometry_and_estimator`
(`pipeline.py:435`), `_score_fragments` (`:469`), `_assign_locus_ids` (`:508`),
`build_multi_loci` (`locus.py:88`), `partition_and_free` (`locus_partition.py:45`),
`run_batch_locus_em_partitioned` (`estimator.py:308`), `FragmentLengthModel.from_pmf`
(`frag_length_model.py:142`).

**Recover from `fc96902` — lean skeleton only:** `_run_locus_em_partitioned`
(the 9-tuple pack + batch call + a **lean** `locus_results` dict — drop the ~20
v5 `prior_*` diagnostic columns); the region→locus overlap allocator
(`adaptive_prior.py::_project_to_loci`, extract the overlap math).

**Write fresh:** `assemble_priors` (the clean bridge below).

**Explicitly DROP** (MEMORY `calibration_no_magic_numbers`): the v5 adaptive
prior (`compute_adaptive_prior`, `MAX_ESS`, `rna_call_bias`, ESS caps), the
per-transcript gDNA *exposure factor*, the v5 `locus_results` diagnostics.

---

# Part I — Theory: the calibration → locus prior bridge

## I.1 The prior is two per-locus scalars — it does **not** touch individual transcripts

A `MultiLocus` is an independent EM subproblem: a set of transcripts + the
overlapping gDNA component. Calibration deconvolves the locus's total fragment
mass into **gDNA mass** and **RNA mass**. **Those two numbers are the prior.**

The prior's job is to **correctly split the unspliced fragments between gDNA and
RNA**. It does **not** — and cannot — say how the RNA mass divides among the
compatible transcripts: that is precisely what the EM solves. So `assemble_priors`
produces exactly **two scalars per locus**:

- `alpha_rna_add[L]`  — total deconvolved RNA mass for the locus,
- `alpha_gdna_add[L]` — deconvolved/modeled gDNA mass for the locus,

which the C++ EM consumes as the Dirichlet prior on the gDNA-vs-RNA-group split
(`AggregatePrior`, `em_solver.cpp:699`). The EM seeds the RNA group across
transcripts by current evidence and then lets the fragment likelihoods do the
attribution. **No per-transcript prior weights are computed anywhere.**

## I.2 Per-region evidence → the two per-locus scalars

Per region `r` (geometry from `RegionArrays`; mass from `CalibrationResult`):

$$
D_r \;=\; M_r^{(d,\mathrm{cont})} + M_r^{(d,L)} + M_r^{(d,R)}
\qquad\text{(deconvolved RNA mass — no ½, D1)}
$$
$$
G_r \;=\; \hat\omega_r\,\rho_0\,L_r^{\mathrm{phys}}
\qquad\text{(modeled gDNA mass; }L_r^{\mathrm{phys}}=\texttt{region\_size\_bp}[r]\text{)}
$$

Region→locus overlap fraction (locus granularity — §I.3):
$\varphi_{L,r} = \mathrm{overlap\_bp}(r,L)\,/\,L_r^{\mathrm{phys}}$.

$$
\boxed{\;\alpha^{(d)}_L = w_{\text{prior}}\sum_r \varphi_{L,r}\,D_r
\qquad
\alpha^{(g)}_L = w_{\text{prior}}\sum_r \varphi_{L,r}\,G_r\;}
$$

`w_prior` (≡ the old κ) is a single global multiplier on prior strength
(§II.7 / D-3). The masses *are* the prior; `w_prior=1` means "the prior carries
exactly the deconvolved calibration mass." No per-region inverse-variance
weighting in the default — see **D-WR** (open).

## I.3 Confirmed corrections to doc 04 §6.2 (it predates D1)

1. **No ½ on boundary mass** (D1). A crossing fragment's mass is attributed
   *fractionally* to each side; the sides are disjoint and sum to the full mass.
   So `D_r = M_cont + M_left + M_right` (full). *(Independent of the accumulator
   interior-slice fix above — that fix makes those per-side masses conserve in
   the first place.)*
2. **Physical `L` in `G_r`** (not FL-corrected `region_eff_length`), **weighted
   by exposure**. The EM's gDNA count scale is `Σ_r M_g_tot_r`, which includes
   boundary-crossers; `M_g_tot_r = ω_r ρ_0 L_r^{phys}` by the exposure mass form.
   FL-correction (which counts only fully-contained fragments) is for the
   calibration *count channel* only.
3. **Locus-granular `φ`** (not per-transcript): the C++ throws away any
   per-transcript split, and `Σ_t φ_{t,r}` would double-count shared exons.
   `φ_{L,r}` counts each region's mass once.

## I.4 Modeled gDNA — to be evaluated empirically

`G_r` uses the **modeled** `ω ρ_0 L`, not the deconvolved `M_g_cont+L+R`,
because gDNA has a forward generative model (smooth, defined on zero-fragment
regions, no double-dip), whereas RNA (no forward model) must use the deconvolved
residual `D_r`. **This is untested** — how well the model behaves is unknown
until we can run the tool. We adopt modeled for v1 and revisit after evaluation
(the deconvolved alternative gives exact conservation but `α_g=0` on empty
regions). — **D-2**.

## I.5 gDNA component effective length = exposure-weighted **physical** bp

*(Corrected — my prior FL-corrected draft was wrong; thank you for the catch.)*

After deconvolution each region's gDNA total is
`M_g_cont + M_g_left + M_g_right`, and the effective length for that total is the
**region size in nucleotides, unmodified** (not FL-corrected) — scaled by the
region's exposure `ω` (gDNA accessibility). So the EM's per-locus gDNA component
effective length is:

$$
\texttt{gdna\_eff\_len}[L] = \max\!\Big(1,\ \sum_r \varphi_{L,r}\,\hat\omega_r\,L_r^{\mathrm{phys}}\Big)
$$

This is *raw physical bp, exposure-weighted* — no FL marginalization. It is the
same quantity that appears inside `G_r` (since `G_r = ρ_0 · ω_r L_r^{phys}`), so
the gDNA prior count is just `w_prior · ρ_0 · (`exposure-weighted physical length
of the locus`)` — fully consistent with how `ρ_0` was fit
(`ρ_0 = Σ M_g_tot / Σ ω L^{phys}`). The `1.0` floor is the EM's existing
numerical default, not a tunable.

---

# Part II — Implementation plan

## II.1 New — `assemble_priors` + `LocusPriors` (`calibration/priors.py`)

```python
@dataclass(frozen=True, slots=True)
class LocusPriors:
    alpha_gdna_add: np.ndarray   # float64[n_loci]
    alpha_rna_add:  np.ndarray   # float64[n_loci]
    gdna_eff_len:   np.ndarray   # float64[n_loci]

def assemble_priors(
    calibration: CalibrationResult,
    region_arrays: RegionArrays,
    multi_loci: Sequence[MultiLocus],
    *,
    prior_weight: float,
) -> LocusPriors: ...
```

Body (≈45 lines): per-region `D_r`, `G_r` (§I.2) as vectors; overlap-project to
loci (§II.2) → the two α arrays + `gdna_eff_len` (§I.5). No per-fragment loop, no
per-transcript anything, no v5 machinery.

## II.2 Recover (lean) — overlap allocator + `_run_locus_em_partitioned`

- **Overlap allocator** (extract `_project_to_loci`): region intervals × each
  `MultiLocus.loci` interval list → `φ_{L,r}` shares. ~25 lines, in `priors.py`.
- **`_run_locus_em_partitioned`** (`pipeline.py`): pack the 9-tuples from
  `partition_and_free`; call `run_batch_locus_em_partitioned(...,
  gdna_prior_count=α_g, rna_prior_count=α_d, gdna_eff_len=..., enable_gdna=None)`
  (`enable_gdna=None` → the estimator auto-derives it, `estimator.py:398`); append
  a **lean** dict per locus to `estimator.locus_results`. `get_loci_df` reads via
  `.get(default)` (`estimator.py:901–917`), so the lean set suffices:
  `locus_id, locus_span_bp, n_transcripts, n_genes, n_em_fragments, rna_total,
  gdna, alpha_gdna_add, alpha_rna_add, enable_gdna, gdna_eff_len_em`.

## II.3 `quant_from_buffer` body (`pipeline.py`)

```
rna_fl  = FragmentLengthModel.from_pmf(fl_models.rna_pmf,  fl_models.max_size)
gdna_fl = FragmentLengthModel.from_pmf(fl_models.gdna_pmf, fl_models.max_size)
geometry, estimator = _setup_geometry_and_estimator(index, rna_fl, em_config)
em_data    = _score_fragments(buffer, index, strand_models, rna_fl, gdna_fl, ...)
multi_loci = build_multi_loci(em_data, index);  _assign_locus_ids(estimator, multi_loci)
priors     = assemble_priors(calibration, region_arrays, multi_loci,
                             prior_weight=em_config_or_calibration_weight)
partitions = partition_and_free(em_data, multi_loci)
_run_locus_em_partitioned(estimator, partitions, multi_loci, index,
                          priors.alpha_gdna_add, priors.alpha_rna_add,
                          priors.gdna_eff_len, em_config, annotations=...)
return estimator, calibration
```

No `set_em_effective_lengths` — the geometry default (FL-marginal transcript
lengths, unit exposure; `estimator.py:207`) is exactly the new design.

## II.4 `run_pipeline` tail

- **Move `buffer.cleanup()`** out of the `finally` around `calibrate()`
  (`pipeline.py:714`) to *after* `quant_from_buffer` — calibrate uses the
  *payload*, quant scans the *buffer*; the current `finally` frees it too early.
- Drop both `NotImplementedError`s; call `quant_from_buffer`, build
  `PipelineResult(estimator=..., calibration=..., ...)`. Output frames come on
  demand from `estimator.get_*_df(index)` (no new result fields).

## II.5 Signature change — `quant_from_buffer`

Add `region_arrays` (prior geometry) + `fl_models` (scorer FL + gDNA eff-len);
drop `frag_length_models` (already distilled into `fl_models` upstream). Only
`run_pipeline` + tests call it.

## II.6 Tests

- **`assemble_priors` unit** (`tests/calibration/test_priors.py`): split
  conservation for a fully-tiled single locus (`α_d+α_g = w_prior Σ_r φ(D_r+G_r)`);
  fractional `φ` for a region split across two loci; zero-fragment region
  (`D_r=0`, `ω=1` → contributes only to `α_g` and `gdna_eff_len`); transcript not
  overlapping any region → 0.
- **End-to-end smoke** (`tests/test_pipeline_smoke.py`): `run_pipeline` completes,
  `locus_results` populated, `get_{counts,gene_counts,loci}_df` well-formed,
  `gdna_em_count ≥ 0`.
- **Golden re-baseline deferred to PR 7** (goldens are stale anyway — D-8).

## II.7 Constants & heuristics — full enumeration (Q6)

| Symbol | Value | New? | Justification |
|---|---|---|---|
| `prior_weight` (≡ κ) | `1.0` | **1 new knob** | Global prior-strength multiplier; `1.0` = prior carries the deconvolved calibration mass. Expose as `EMConfig.prior_weight` (D-3). |
| `gdna_eff_len` floor | `1.0` | No | Already the EM default (`estimator.py` `gdna_eff_len=None` path / `set_em_effective_lengths`). Numerical guard. |

Everything else (`D_r`, `G_r`, `φ_{L,r}`, `gdna_eff_len`) is derived. PR 6 adds
**one** knob, defaulting to an inert `1.0`.

---

# Part III — Decisions

**Resolved (your feedback):**
- **D-1** gDNA prior uses physical locus length, **exposure-weighted** (§I.3.2). ✓
- **D-2** Try the **modeled** gDNA estimate; evaluate after we can run the tool (§I.4). ✓
- **D-3** Expose `EMConfig.prior_weight` (intuitive name for the prior multiplier),
  default `1.0`. ✓
- **D-4** Drop the ½ boundary half-split + use locus-granular `φ` (§I.3.1/3). ✓
- **D-5** `gdna_eff_len` = exposure-weighted **physical** bp, region size unmodified
  (§I.5). ✓
- **D-6** `quant_from_buffer` signature change (add `region_arrays` + `fl_models`). ✓
- **D-7** Lean `locus_results`; add diagnostics later if needed. ✓
- **D-8** Golden re-baseline anytime (goldens stale); PR 6 bar = runs + smoke. ✓

- **D-ACC** Accumulator mass conservation fixed first as a dedicated **PR 5.5**
  (merged). PR 6 builds on the conserved per-region mass. ✓
- **D-WR** **Omit** the inverse-variance weight `w_r` for v1 — the prior is the
  plain deconvolved mass (`α = prior_weight · Σ_r φ_{L,r} · mass`). Revisit only
  if evaluation shows uncertain-exposure regions distorting the prior. ✓

**All decisions are resolved. There are no open questions.**

---

## Readiness certification — cleared to implement PR 6

| Gate | Status |
|---|---|
| All design decisions resolved (D-1…D-8, D-ACC, D-WR) | ✓ |
| Prerequisite (accumulator conservation, PR 5.5) merged to `main` | ✓ |
| Downstream contract mapped (EM API, router, loci, partition, `from_pmf`) | ✓ verified against live code (`estimator.py:308`, `scan.py`, `locus.py:88`, `locus_partition.py:45`, `frag_length_model.py:142`) |
| Prior role locked: two per-locus scalars, no per-transcript weights (§I.1) | ✓ |
| Constants budget (Q6): exactly one new knob, `EMConfig.prior_weight=1.0` | ✓ |
| Recover targets identified (`fc96902`: `_run_locus_em_partitioned` lean, `_project_to_loci`) | ✓ |
| Test plan: `assemble_priors` units + end-to-end smoke; golden re-baseline → PR 7 | ✓ |

**Build/verify in the activated `rigel` env** (per `CLAUDE.md`): after wiring,
`pip install --no-build-isolation -e .` (no C++ change expected — Python-only PR,
unless `EMConfig.prior_weight` is the only edit), then `pytest tests/` should turn
the ~481 `NotImplementedError` failures green, and `ruff check`.

**PR 6 is cleared to begin.** Implementation order: (1) `assemble_priors` +
`LocusPriors` + overlap allocator; (2) recover lean `_run_locus_em_partitioned`;
(3) `quant_from_buffer` body + signature; (4) `run_pipeline` tail (move
`buffer.cleanup()`, drop both `NotImplementedError`s, build `PipelineResult`);
(5) `EMConfig.prior_weight`; (6) tests.

---

## Rollback

Pure additive integration. Revert by restoring the two `NotImplementedError`s and
the `finally: buffer.cleanup()`. The calibrator and PR 0–5 / PR 5.5 modules are
untouched.
