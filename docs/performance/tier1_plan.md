# Tier 1 Implementation Plan + `coverage_weights` Audit

**Date:** 2026-05-13
**Scope:** Bug fixes surfaced by the VCAP profile run, plus the four Tier 1 items from [profile_2026-05-13_vcap_rna20m_gdna20m.md](profile_2026-05-13_vcap_rna20m_gdna20m.md). Plus a deep audit of the `coverage_weights` candidate-array because it is the largest single contributor to peak RSS during scoring.

This plan is meant to be executed top-to-bottom in a single PR (or one PR per section). Each item has explicit acceptance criteria so we can verify the change without re-running the full benchmark.

---

## 0. Pre-flight bug fixes (already applied)

Both fixes were necessary to even produce the profile and are already in tree. Documenting here so they show up in the PR description.

### 0.1 `int8 → int32` in `build_t_to_local_locus`

**File:** [src/rigel/calibration/_locus_n_obs.py](../../src/rigel/calibration/_locus_n_obs.py)

**Symptom:** `RuntimeError: build_t_to_local_locus: MultiLocus 3 has 1497 loci; int8 binning supports ≤ 127.` on the VCAP workload.

**Root cause:** A defensive guard from a prior assumption that "the largest GENCODE cluster is 19 loci" — invalidated by real data where one MultiLocus contains 1 497 distinct contiguous Locus intervals (likely a paralog cluster or pseudo-gene-rich region).

**Fix:** Promote the bin-id dtype from `int8` to `int32` and remove the cap. Memory cost: ≤ 460 K transcripts × 4 bytes = 1.8 MB, completely negligible.

**Acceptance:** Run `pytest tests/test_locus_partition.py tests/test_partition_units.py` — must pass. Then run quant on the VCAP BAM — must not raise.

### 0.2 Stale `c_base` kwarg in profiler

**File:** [scripts/profiling/profiler.py](../../scripts/profiling/profiler.py)

**Symptom:** `AttributeError: 'CalibrationConfig' object has no attribute 'c_base'` at the end of the stage-decomposed profile.

**Root cause:** Profiler was written against an older calibration API; v6 renamed `c_base` to `nrna_weight`. Replaced and added the missing `strand_summary=` kwarg to `calibrate(...)`.

**Acceptance:** Profiler completes end-to-end on a small BAM with `--stages --cprofile`.

---

## 1. Tier 1 fixes

Listed in implementation order. Each is independently mergeable.

### 1.1 Memoise `FragmentLengthModel._build_eff_len_cache` after finalize

**File:** [src/rigel/frag_length_model.py:399](../../src/rigel/frag_length_model.py#L399)

**Current state:**
```python
def _build_eff_len_cache(self) -> tuple[np.ndarray, np.ndarray]:
    probs = self._normalized_probs()      # (max_size+1,)
    cdf = np.cumsum(probs)
    l_vals = np.arange(len(probs), dtype=np.float64)
    cmom = np.cumsum(probs * l_vals)
    return cdf, cmom
```

**Profile evidence:** 169 460 calls during the prior assembly stage; 567 K `cumsum` invocations cumulative; 2.0 s self time.

**Why it's safe to cache:** Once `finalize()` is called, the FL model is immutable through the rest of the pipeline. Calibration v6 builds `FLModels` at the end of `calibrate()` and never re-trains.

**Implementation:**

1. Add two lazy slots to `FragmentLengthModel.__init__`: `self._cdf_cache: np.ndarray | None = None`, `self._cmom_cache: np.ndarray | None = None`.
2. Modify `_build_eff_len_cache` to return cached values when populated, else compute, store, and return. Keep the cache invariant: cache is populated only when `self._finalized is True`.
3. Modify `_invalidate_cache()` (or every `observe_*` method that mutates the histogram) to null out both caches. Verify all mutators by grepping for `self.counts`, `self._prob`, and `_finalized`.
4. Verify by repeating the call: cached path must produce arrays `array_equal` with the freshly-computed ones.

**Acceptance criteria:**
- All existing unit tests pass: `pytest tests/test_fl_models.py tests/test_frag_length_model.py tests/test_transcript_space_fl.py -v`.
- Add a small regression test asserting `_build_eff_len_cache()` returns the *same* `id(...)` on repeated calls after `finalize()`.
- Re-run the VCAP stage profile; `_build_eff_len_cache` cumulative time should drop from 1.96 s to < 50 ms; `cumsum` total calls should fall from 567 K to ~1 K.

**Risk:** Trivially low. The cache is a pure function of the histogram; mutating the histogram already requires going through `observe_*` which we will instrument.

**Estimated payoff:** −2 s wall on `eb_gdna_priors`; no RSS change.

---

### 1.2 Cast EM `log_liks` and `coverage_weights` to float32

**Files:**
- [src/rigel/scored_fragments.py](../../src/rigel/scored_fragments.py) — change dtype annotations
- [src/rigel/native/scoring.cpp](../../src/rigel/native/scoring.cpp) — `v_ll`, `v_cw` push_back as float
- [src/rigel/native/em_solver.cpp](../../src/rigel/native/em_solver.cpp) — `coverage_wts` and `log_liks` ABI: accept `f32_1d`, promote to `double` only on read into the per-row buffer
- [src/rigel/locus_partition.py](../../src/rigel/locus_partition.py) — switch `scatter_candidates_f64` to `scatter_candidates_f32` for these two
- [src/rigel/scan.py](../../src/rigel/scan.py), [src/rigel/pipeline.py](../../src/rigel/pipeline.py) — wire dtype through
- [tests/conftest.py](../../tests/conftest.py), all tests using `np.ones(..., dtype=np.float64)` for these arrays — switch to float32

**Profile evidence:** `log_liks` and `coverage_weights` are 260 M × 8 B = **2.08 GB each**. Casting each to float32 saves ≈ **2 GB peak RSS**.

**Numerical analysis:**

- **`log_liks`**: Magnitudes in our pipeline are ≤ ~50 in absolute value (sum of log-strand ≈ −7, log-FL ≈ −7, overhang penalty ≈ −2.3 × N_overhang_bp, NM penalty ≈ −2.3 × N_mismatches). Float32 has 23 mantissa bits ⇒ ~7 decimal digits of precision; relative error in `log_lik` is ~1e-7. The EM kernel immediately subtracts a per-row `max_val` (line 396-399 in `em_step_kernel_range`) before exp, so the absolute error in the post-pivot exponent is < 1e-5 ⇒ relative error in posterior is < 1e-5. SQUAREM convergence delta is 1e-5 by default; this is well within tolerance.
- **`coverage_weights`**: Trapezoid model returns values in `[1.0, ∞)` but in practice ≤ ~3 (ratio of plateau width to local coverage area). Float32 is overwhelmingly precise here.

**Implementation steps:**

1. **C++ scoring write side** (scoring.cpp lines 932-933): change `std::vector<double>* v_cw` to `std::vector<float>* v_cw` and `v_ll` similarly. The push_back values become `static_cast<float>(s.cov_wt)` and `static_cast<float>(s.log_lik)`.
2. **Streaming scorer output** (around `finish()`): expose float32 ndarrays.
3. **C++ EM read side** (em_solver.cpp `build_equiv_classes` line 247-248): keep `wt_flat`/`ll_flat` as `std::vector<double>` (they live in equiv-class scratch and are read repeatedly by SQUAREM), but cast on copy from the input float32 arrays. The flat arrays are O(units × k) and not the dominant cost; conversion happens once.
4. Update the nanobind signatures (em_solver.cpp line 1086, 1971, 2598) to take `f32_1d` instead of `f64_1d` for `log_liks` and `coverage_wts`.
5. Update Python-side tests (conftest.py, test_em_impl.py, test_em_prior_weight.py, etc.) to allocate float32 — search for `np.ones(..., dtype=np.float64)` followed by `coverage_weights` or `log_liks` usage.
6. Update [tests/golden/](../../tests/golden/) regression files only if numeric drift > 1e-5 (it shouldn't, given the pivot).

**Acceptance criteria:**
- Full test suite passes (`pytest tests/ -v`). Golden output tests should pass without `--update-golden` (drift below tolerance). If they fail, check the diff is < 1e-5 relative.
- Re-run VCAP profile; `peak_rss_mb` should drop from 15 836 to ~13 800.
- EM convergence delta on a representative locus must not change by more than 10 %.

**Risk:** Medium. The ABI change ripples through many test fixtures. The EM is well-conditioned (per-row pivot, log-sum-exp, deterministic equiv-class sort) but we must validate with golden tests before merging.

**Estimated payoff:** −2 GB peak RSS; small (~5 %) speedup in scoring + EM kernel due to halved memory bandwidth.

---

### 1.3 Vectorise `partition_units_to_loci` slow path

**File:** [src/rigel/calibration/_locus_n_obs.py:111](../../src/rigel/calibration/_locus_n_obs.py#L111)

**Current state:**
```python
unit_indices = ml.unit_indices
return tuple(
    np.ascontiguousarray(unit_indices[bins == j], dtype=np.int32)
    for j in range(n_loci)
)
```

For each locus `j` in `range(n_loci)`, this allocates a `bool[n_units]` mask, fancy-indexes `unit_indices`, and casts to int32. On the largest mega-MultiLocus (`1497 loci × ~thousands of units`), this is the dominant call.

**Profile evidence:** 7 351 generator-expression invocations; 1.25 s self time.

**Replacement:**
```python
unit_indices = ml.unit_indices
order = np.argsort(bins, kind="stable")        # int64 perm
sorted_bins = bins[order]
sorted_units = unit_indices[order].astype(np.int32, copy=False)
# Per-bin slice offsets via searchsorted on a strictly-increasing range
edges = np.searchsorted(sorted_bins, np.arange(n_loci + 1, dtype=bins.dtype))
return tuple(
    np.ascontiguousarray(sorted_units[edges[j]:edges[j + 1]])
    for j in range(n_loci)
)
```

This is one O(N log N) sort + O(n_loci log N) searchsorted instead of O(n_loci × N) mask scans. Output is **identical** (same units, same per-locus order — `argsort(stable)` preserves original order within each bin).

**Acceptance criteria:**
- New unit test: random `bins` array of size 10 000 across 100 buckets — assert `set(new[j]) == set(old[j])` for every `j` and `np.array_equal(new[j], old[j])` (stable order).
- Existing tests pass: `pytest tests/test_partition_units.py tests/test_locus_partition.py tests/test_locus_rename.py`.
- Re-run VCAP profile; the partition_units_to_loci genexpr should drop from 1.25 s to < 100 ms.

**Risk:** Low. Pure numpy; preserves stable order.

**Estimated payoff:** −1.2 s wall in `eb_gdna_priors`.

---

### 1.4 Fuse `estimate_locus_gdna` + `expected_gdna_count_global` per Locus

**File:** [src/rigel/calibration/locus_prior.py](../../src/rigel/calibration/locus_prior.py)

**Profile evidence:** Each call costs ≈ 0.17 ms × 28 K calls. The two functions duplicate:
- `region_index.overlap(...)` (called twice per Locus, once each)
- `contained_exposure_clipped(starts, ends, locus_start, locus_end, gdna_fl)` (the eff_full + eff_clip_core compute is identical across both)
- `boundary_crossing_exposure(gdna_fl, splicing_anchor_tolerance=...)` (returns a scalar — easy to hoist out of the per-Locus loop entirely)

The two functions share **almost every input**; the divergence is in how they aggregate. Specifically:
- `estimate_locus_gdna` consumes `payload_arrays` (per-region observed counts) and applies locoregional shrinkage.
- `expected_gdna_count_global` ignores `payload_arrays` and projects global densities.

**Implementation steps:**

1. **Hoist `boundary_crossing_exposure` out of the per-Locus loop.** It's a pure function of `gdna_fl` and `splicing_anchor_tolerance` — both immutable here. Compute once before the multi-locus loop and pass as a scalar.
2. **Hoist FL-cache** by depending on item 1.1 above (after that lands, re-running `compute_all_transcript_eff_lens` per call is cheap; we can keep two separate functions).
3. **Introduce `_estimate_locus_combined`** that takes the union of arguments (region arrays, payload, densities, fl), runs `region_index.overlap` and `contained_exposure_clipped` exactly once, then computes both the locoregional `LocusGdnaEstimate` and the global `ExpectedGdnaPriorParts` from the shared scratch.
4. Wire `assemble_priors` to call the combined helper inside the per-Locus loop and unpack two return values.
5. Keep `estimate_locus_gdna` and `expected_gdna_count_global` exported for unit tests; the combined helper is an internal optimisation.

**Acceptance criteria:**
- Unit-level: every existing call to `estimate_locus_gdna` and `expected_gdna_count_global` produces the identical numeric output (compare new combined helper's two halves to the originals on a representative Locus). Use `pytest tests/test_density_loco.py tests/test_density_global.py tests/test_locus_partition.py`.
- Re-run VCAP profile; `eb_gdna_priors` should drop from 9.5 s to ≈ 4 s (boundary_crossing alone removes 0.33 s; the rest comes from de-duplicating the FL exposure + region overlap).

**Risk:** Medium. The two helpers have subtly different windowing (intergenic flank vs. unflanked). The combined helper must preserve both.

**Estimated payoff:** −3 s wall in `eb_gdna_priors`.

---

## 2. `coverage_weights` audit

The Tier 1 plan above casts `coverage_weights` to float32 (item 1.2). Before doing that we need to understand whether this 2 GB array carries its weight (pun intended) or whether it can be dropped entirely.

### 2.1 Where `coverage_weights` is computed

[src/rigel/native/scoring.cpp](../../src/rigel/native/scoring.cpp) — function `compute_fragment_weight` (lines 99-145).

Per fragment-candidate pair `(fragment, transcript)`, it computes the **trapezoid coverage capacity** — under uniform random sampling, the expected number of fragments of the observed length that *could* have been generated at the observed position. Edge regions of a transcript have lower sampling capacity than the plateau, so a fragment observed at the edge gets a weight `> 1` (less likely to be drawn ⇒ stronger evidence). Plateau fragments get weight `1.0`.

Storage path:
1. **Scoring** (scoring.cpp ~line 596): `cov_wt = compute_fragment_weight(tx_s, cov_end, t_len)` per (frag, candidate), then pushed to `v_cw` for the chunk.
2. **CSR build** (scan.py): `coverage_weights` becomes a per-candidate float64 array in `ScoredFragments` — 260 M × 8 B = **2.08 GB**.
3. **Locus partition** (locus_partition.py): scattered to per-locus `LocusPartition.coverage_weights` (still float64).
4. **EM** (em_solver.cpp `build_equiv_classes` line 247-248): copied from per-locus arrays into `EmEquivClass.wt_flat` (still float64).

### 2.2 Where `coverage_weights` is consumed

Exhaustive grep result. Three call sites in C++; no Python consumer:

| File / function | Line(s) | What it does | Observable effect |
|---|---|---|---|
| `em_solver.cpp::build_equiv_classes` | 237-248 | Stores `wt_flat[i*k+j] = coverage_wts[s+j]` into the per-EC dense matrix | Pre-stage for the next two |
| `em_solver.cpp::build_equiv_classes` (sort tie-break) | 281-296 | Used as a **secondary** sort key for deterministic row ordering when `ll_flat` rows are exactly equal | Determinism only — no effect on EM convergence |
| `em_solver.cpp::compute_ovr_prior_and_warm_start` | 763-806 | (a) Per-unit row-normalised coverage shares produce the **warm-start `theta_init`**; (b) component coverage totals feed the **OVR prior allocation** weights | Initial point and per-component RNA prior strength |

**Critical fact:** `coverage_weights` is **NOT used in the EM E-step inner loop** (`em_step_kernel_range`, lines 357-419). The hot loop reads `ll_flat` only. So the 2 GB pays for one warm-start pass plus an occasional sort tie-break.

### 2.3 What the warm start actually does

`compute_ovr_prior_and_warm_start` runs once at EM startup. For every equiv class, every unit row, it computes:
- A per-row sum `row_sum = Σ wt[i,j] * eligible[cidx[j]]`
- Per-cell shares `share = wt[i,j] * eligible[cidx[j]] / row_sum`
- Adds `share` into `theta_init_out[cidx[j]]` (warm start)
- Adds `share` into `coverage_totals[cidx[j]]` (used downstream for OVR prior allocation across RNA components)

If every `coverage_weight == 1.0`, the warm start collapses to **uniform fractional assignment**: each candidate gets `1/k` of the unit, exactly as if we passed `np.ones`. The OVR prior also collapses to a coverage-uniform prior allocation across RNA components.

So the question is: **does the trapezoid edge-aware weighting actually move EM convergence in a measurable way?**

### 2.4 Hypothesis

The trapezoid weighting was useful when the EM started cold (all theta uniform) and convergence was sensitive to the warm start. With Calibration v6 in place, every EM run has:
- A **strong, prior-informed gDNA pseudocount** (`alpha_gdna = η_g(ℓ)`).
- An **OVR prior** that is itself a function of `coverage_totals`.

The OVR prior is the only place where coverage-weighted information persists into the M-step. For an isolated transcript with no overlapping isoforms, every candidate has weight 1 anyway. For overlapping isoforms, the weight ratios across candidates within one fragment are *small* — the trapezoid model usually returns 1.0 for plateau fragments and ~1.05–1.5 for edge fragments. The induced asymmetry is dwarfed by the OVR prior's dependence on `total_rna_coverage`.

### 2.5 Recommended experiment (must run before deciding T1.2 vs. removal)

Three configurations on the VCAP BAM and one or two sim-suite conditions:

| Config | `coverage_weights` source | Expected outcome |
|---|---|---|
| **A — baseline** | trapezoid (current) | reference |
| **B — float32** | trapezoid, downcast to float32 | should be bit-equivalent within EM tolerance |
| **C — uniform** | hard-code `cov_wt = 1.0` in scoring.cpp | tells us whether the weighting matters at all |

For each config compute:
- Per-transcript abundance ranks (Spearman ρ vs A) — should be > 0.999 for B and we need to see how close C is.
- Per-locus EM iteration count — if C uses ≤ 10 % more iterations on average, the warm start is ornamental.
- Mass conservation and overall mRNA / nRNA / gDNA fractions — should match within 0.5 %.

**Decision rule:**
- If C is statistically equivalent to A on the metrics above, **drop `coverage_weights` entirely** (replace with hard-coded `1.0` in the EM warm start, eliminate the array end-to-end). This saves 2 GB peak RSS and removes a column from the buffer + scoring + EM pipelines.
- If C is meaningfully worse but B is equivalent to A, ship the **float32 cast** (T1.2 above) and keep the trapezoid weighting. Saves 1 GB peak RSS.
- If both B and C diverge from A, ship neither; coverage_weights is structurally important.

The experiment is cheap to set up: a single C++ flag `RIGEL_DISABLE_COV_WT` that forces `cov_wt = 1.0` in the two scoring code paths. Run on the existing benchmark fixtures (`scripts/benchmarking`).

### 2.6 Code-archaeology note

The trapezoid model appears to be a port of salmon's Equivalence Classes Coverage (ECC) weighting (salmon, McManes & Patro 2020 — though salmon stores weights as float32, not float64). The salmon team validated it as part of the *EM warm start*, with the explicit acknowledgement that their EM is otherwise bias-corrected. Rigel inherited the construct in the early "salmon-style" scoring rewrite, but the bias correction in calibration v6 (effective length, FL-aware, gDNA-aware) supersedes the warm-start argument. **My prediction is that config C will match A within tolerance**, in which case we can drop the column entirely.

---

## 3. Suggested PR layout

| PR | Items | Estimated diff size |
|---|---|---|
| **PR-1: Hotfixes** | §0.1 + §0.2 | ~50 LOC, no test changes besides one regression test |
| **PR-2: Tier 1 Python** | §1.1 + §1.3 + §1.4 | ~250 LOC, ~5 new tests |
| **PR-3: coverage_weights audit** | §2.5 experiment + decision write-up | ~30 LOC C++ flag, ~50 LOC analysis script under `scripts/benchmarking/` |
| **PR-4: Tier 1 dtype** | §1.2 *or* coverage_weights removal (whichever §2.5 selects) | ~400 LOC across C++ + Python + tests + golden file refresh |

Estimated combined wall-time saving on VCAP after all four PRs: **~7 s on quant** (46 s → ~39 s); peak RSS savings: **2–4 GB** depending on §2.5 outcome.

---

## 4. Definition of done

A Tier-1 PR is done when:

1. The full test suite passes (`pytest tests/ -v --update-golden=False`), excluding the documented pre-existing strand-LLR test failure.
2. `ruff check src/ tests/` is clean.
3. The relevant cells in [profile_2026-05-13_vcap_rna20m_gdna20m.md](profile_2026-05-13_vcap_rna20m_gdna20m.md) §1 ("Stage breakdown") and §2 ("Hotspot deep-dive") are re-validated by a follow-up profile run (no need to re-publish the doc; a JSON diff of `profile_summary.json` in the PR description is sufficient).
4. Memory and wall-time deltas match the predictions within 25 %; if not, root-cause before merging.

---

## 5. Out of scope (deferred to Tier 2+)

- Buffer dtype narrowing (Tier 2.1) — separate PR; touches the C++ buffer struct.
- Streaming scorer parallelisation (Tier 3.1) — needs an EM-input concatenation strategy.
- BAM-scan I/O re-architecture (Tier 4.2) — a multi-day project on its own.
