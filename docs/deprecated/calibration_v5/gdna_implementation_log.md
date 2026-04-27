# gDNA Option B — Implementation Log

Companion to [`gdna_length_model_options.md`](gdna_length_model_options.md) §6.

---

## Phase 1 — Scorer-side physics

**Goal**: per-hit harmonic-mean gDNA length correction inside `scoring.cpp`;
`gdna_log_liks` arrives pre-corrected at the buffer.

### Changes

- `scoring.cpp`:
  - New file-local helper `lse_update(max, sum, x)` for numerically stable
    online logsumexp. Handles the first-input (`-inf` → finite) transition
    explicitly.
  - New `gdna_flank_` int32 member on `NativeFragmentScorer`, injected via
    one new constructor argument (inserted between `gdna_fl_tail_base` and
    `t_strand_arr`; nanobind `.def(nb::init<…>)` updated accordingly).
  - `FillState`:
    - **Added**: `mm_gdna_lse_max` (double), `mm_gdna_lse_sum` (double),
      `mm_nh_gdna` (int32).
    - **Kept as bridge** for Phase 1 standalone: `mm_best_gdna_fp`,
      `mm_first_gfp` (still consumed by `v_gfp` push; deleted in Phase 3
      when `v_gfp` disappears).
    - **Deleted**: `mm_best_gdna_ll` (replaced by lse accumulator).
  - `score_mm_alignment`: per-hit gDNA block rewritten. Anchors the
    sampling window on the hit's first RNA candidate transcript
    (`t_span_[anchor_t] + gdna_flank_`), computes
    `e_h = max(L_h − gfp + 1, 1)`, and accumulates
    `gdna_fl + gdna_log_sp + LOG_HALF + log_nm − log(e_h)` via
    `lse_update`. `mm_best_gdna_fp` is tracked as the fp of the
    current-lse-max hit (Phase-3 bridge). Hits with zero RNA candidates
    (unresolved) are skipped — they cannot contribute to the gDNA
    sampling window anyway.
  - `flush_mm_group`: emits
    `final_gdna_ll = lse_max + log(lse_sum) − log(nh_gdna)` when
    `nh_gdna > 0 && !mm_is_any_spliced`, else `−∞`.
  - Non-MM path: emits
    `gdna_fl + gdna_log_sp + LOG_HALF + log_nm − log(e_h)` using
    `best_t` as anchor.

- `src/rigel/scoring.py`:
  - `FragmentScorer.from_models` reads `gdna_flank = int(gdna_fl.mean)`
    (safe on zero-weight models — returns 0) and forwards it as
    `gdna_flank=` to `NativeFragmentScorer`.

### Compile + test

- `pip install --no-build-isolation -e .` — compiled clean.
- `pytest tests/ -x --ignore=tests/test_golden_output.py` failed at
  `tests/scenarios/test_nrna_double_counting.py::TestNrnaDoubleCounting::test_full_sweep[g20_n0_s65]`
  with `rna_rel_err=0.65` (expected < 0.60).

**Diagnosis**: expected behavior for Phase 1 standalone.
`em_solver.cpp::apply_bias_correction_uniform` still subtracts
`log(gdna_span − gfp + 1)` on top of the now-already-Option-B-corrected
`gdna_log_liks[u]`. The double correction makes gDNA log-likelihood
artificially low, pushing fragments to the RNA hypothesis and
over-estimating total RNA. Phase 2 removes this redundant correction.

### Status

✅ Phase 1 complete. No blockers. Proceeding to Phase 2.

---

## Phase 2 — EM solver simplification

**Goal**: strip redundant gDNA length correction from the EM solver; the
scorer now pre-corrects `gdna_log_liks`, so `apply_bias_correction_uniform`
must skip the gDNA row and the per-locus `gdna_span` plumbing becomes
dead weight.

### Changes

- `em_solver.cpp`:
  - `apply_bias_correction_uniform`: new `int32_t n_t` parameter
    inserted before `n_candidates`. Rows with `t_indices[i] >= n_t`
    are skipped — this covers the gDNA row (local component == n_t)
    whose log-lik is already length-corrected upstream.
  - `extract_locus_sub_problem_from_partition`: dropped `int64_t gdna_span`
    parameter. `sub.bias_profiles` now sized to `n_t` (RNA only);
    deleted `sub.bias_profiles[sub.gdna_idx] = gdna_span;`.
  - gDNA sort_buf entry at the candidate-scatter step zeroes out the
    fragment-length coordinates (`tx_s=0, tx_e=0`). Value is no longer
    used for bias correction; the `pv.genomic_footprints[ui]` read was
    removed from this site (bridge field still compiled, deleted in
    Phase 3).
  - `batch_locus_em_partitioned`: dropped `i64_1d gdna_spans` parameter,
    local `gs_ptr`, and one `nb::arg("gdna_spans")` binding entry.
  - Non-partitioned `run_locus_em_native` path (test-only/legacy)
    passes `static_cast<int32_t>(nc)` to `apply_bias_correction_uniform`
    so every component retains bias correction — this entry point has
    no gDNA-row convention.

- `src/rigel/estimator.py`:
  - `run_batch_locus_em_partitioned` wrapper dropped `gdna_spans`
    parameter (+ docstring entry + one `ascontiguousarray` call).

- `src/rigel/pipeline.py`:
  - `_run_locus_em_partitioned`: dropped `gdna_flank: int = 0` kwarg.
  - `_call_batch_em`: dropped `batch_spans` positional parameter.
  - Mega-locus call site: deleted
    `np.array([locus.gdna_span + gdna_flank], dtype=np.int64)`.
  - Normal-locus batch: deleted `normal_spans` construction.
  - Top-level: deleted `gdna_flank = int(calibration.gdna_fl_model.mean)`
    block and the `gdna_flank=gdna_flank` kwarg (flank now lives inside
    `FragmentScorer.from_models` — Phase 1).

- `tests/conftest.py`:
  - `_run_and_assign` no longer constructs `gdna_spans` nor forwards it.

- `tests/test_gdna.py`:
  - `test_gdna_log_lik_determines_absorption`: test-supplied
    `gdna_log_lik` values shifted by the removed −log(eff_len) ≈ −9.19
    constant so the Option B pre-corrected inputs preserve the test's
    low/high contrast (−5.0 → −14.0, 0.0 → −9.0).
  - `test_tiny_alpha_gdna_suppresses_gdna`: explicit
    `gdna_log_lik=-9.0` added to represent a realistic pre-corrected
    "no-signal" value.

### Compile + test

- `pip install --no-build-isolation -e .` — clean compile.
- `pytest tests/ --ignore=tests/test_golden_output.py` — **1037 passed**,
  including the previously-failing
  `test_nrna_double_counting[g20_n0_s65]`.

### Status

✅ Phase 2 complete. No blockers. Proceeding to Phase 3
(delete `genomic_footprints` end-to-end).

---

## Phase 3 — Delete `genomic_footprints` plumbing

**Goal**: remove the per-unit `genomic_footprints` array entirely. With
Option B's in-scorer harmonic-mean length correction, downstream code
no longer needs per-unit genomic footprint information — it was used
solely to size the now-deleted gDNA bias profile row.

### Changes

- `scoring.cpp`:
  - Removed `FillState::v_gfp`, `mm_best_gdna_fp`, `mm_first_gfp`
    (Phase-1 bridge fields). Constructor / `reset_mm_group` updated.
  - Removed the bridge write `st.mm_best_gdna_fp = gfp_val` at the
    current-lse-max hit and the "track first member's footprint"
    fallback.
  - Removed `st.v_gfp->push_back(...)` in `flush_mm_group` and the
    non-MM emission path.
  - `StreamingScorer`: removed `v_gfp_` member, its allocation,
    wiring, deletion, `vec_to_ndarray` output slot, and the null-out
    on finish.

- `em_solver.cpp`:
  - `PartitionView`: removed `const int32_t* genomic_footprints` field.
  - Tuple parsing in `batch_locus_em_partitioned`: removed `tup[9]`
    read; shifted `locus_t_indices` / `locus_count_cols` to `tup[9]`,
    `tup[10]`. Tuple arity dropped from 12 to 11. Docstrings updated.

- `src/rigel/scored_fragments.py`: removed `genomic_footprints`
  field from `ScoredFragments` and `LocusPartition` dataclasses
  (including docstring entries).

- `src/rigel/scan.py`: dropped `genomic_footprints` from the
  `StreamingScorer.finish()` tuple unpacking and the
  `ScoredFragments(...)` constructor call.

- `src/rigel/partition.py`: removed `genomic_footprints` from the
  `UNIT_ARRAYS` scatter list and the `LocusPartition(...)` ctor.

- `src/rigel/pipeline.py`: removed `p.genomic_footprints` slot from
  the partition tuple in `_call_batch_em`.

- `src/rigel/estimator.py`: updated "12-tuples" docstring to
  "11-tuples".

- Tests (`tests/conftest.py`, `tests/test_partition.py`,
  `tests/test_estimator.py`): dropped `genomic_footprints=...`
  from every `ScoredFragments(...)` / `LocusPartition(...)` /
  partition-tuple construction site, plus the round-trip assertion
  in `test_partition`.

### Compile + test

- `pip install --no-build-isolation -e .` — clean compile.
- `pytest tests/ --ignore=tests/test_golden_output.py` — **1037 passed**.

### Status

✅ Phase 3 complete.

---

## Phase 4 — Validation

### 4a. Harmonic-mean accumulator unit tests

New file `tests/test_gdna_harmonic_length.py` — 7 focused tests that
exercise the Option B LSE accumulator end-to-end through the real
C++ `NativeFragmentScorer`:

1. **NH=1 unique unspliced** — emitted `gdna_log_lik` is finite.
2. **NH=2 identical spans** — collapses to the NH=1 value
   (`log((1/e + 1/e)/2) = −log(e)`).
3. **NH=2 dissimilar spans** — exact match to hand-computed
   `C + log((1/e₀ + 1/e₁)/2)` using a reference NH=1 unit to
   recover the per-hit constant `C`.
4. **Any-spliced in MM group** — `gdna_log_lik = −∞`.
5. **All-spliced MM group** — `gdna_log_lik = −∞`.
6. **`gdna_flank` widens window** — verifies
   `e(flank) > e(no_flank)` and the length term becomes more
   negative (wider uniform prior).
7. **NH=3 three-way harmonic mean** — three-hit case matches the
   closed-form arithmetic mean of reciprocals.

### 4b. Full test suite + goldens

- `pytest tests/test_gdna_harmonic_length.py -v` — **7 passed**.
- `pytest tests/test_golden_output.py --update-golden` — 21
  goldens regenerated (the 4 that changed were rotated-mass
  shifts in combined-stress scenarios, within ~0.1% of prior values).
- `pytest tests/` — **1065 passed** (all prior tests + 7 new + 21
  regenerated goldens that were previously skipping).

### 4c. vcap mixture titration (real data)

Re-ran `rigel quant` with the new Option B build against the
8-library vcap gDNA titration under
`/scratch/.../runs/human/mctp_vcap_rna20m_dna{00,01,02,05,10,20,40,80}m/`
(STAR-aligned, name-sorted, NH-tagged, multimappers retained).
Each library has 20M RNA reads + a variable gDNA contamination
spike-in. Outputs written to `/scratch/.../runs/human_optionb/`;
compared against the pre-Option-B baseline in `<lib>/rigel/`.

#### Headline metrics (mRNA / nRNA should be flat across titration;
gDNA should rise linearly with input)

| Library | Δ mRNA (vs 0M, old → new) | Δ nRNA (vs 0M, old → new) | gDNA slope (old → new) |
|---|---|---|---|
| dna00m | baseline | baseline | — |
| dna01m | +0.09M → +0.04M | +0.10M → +0.01M | 0.738 → 0.875 |
| dna10m | +0.73M → +0.37M | +0.80M → +0.22M | 0.752 → 0.847 |
| dna20m | +1.38M → +0.56M | +1.67M → +0.40M | 0.734 → 0.838 |
| dna40m | +2.56M → +1.04M | +3.33M → +0.78M | 0.702 → 0.803 |
| dna80m | +4.60M → +1.86M | +6.45M → +1.39M | 0.644 → 0.742 |

Interpretation:

- **mRNA drift** (should be 0 in a constant-RNA titration): cut
  from +34% (old) to +14% (new) at 80M gDNA — **2.5× smaller
  contamination-induced bias**.
- **nRNA siphon** (the dominant failure mode in synthetic
  benchmarks): cut from +1850% to +470% — **4.6× smaller**. nRNA
  no longer absorbs gDNA fragments as aggressively.
- **gDNA recovery**: per-input-read slope up from 0.64–0.76 (old)
  to 0.74–0.88 (new). Option B correctly classifies an extra
  ~10% of input gDNA that the old code misrouted into mRNA/nRNA.

The calibration's `E[gdna]` estimate remains ~10× smaller than
the observed `gdna_total` across the titration; this ratio is
stable (0.10 ± 0.01) and was present before Option B. It's a
separate calibration-layer issue, not in scope here.

Outputs:
- `scripts/benchmark/run_vcap_optionb.sh` — driver
- `scripts/benchmark/compare_vcap_optionb.py` — comparator
- `/scratch/.../runs/human_optionb/comparison.tsv` — raw table

### Status

✅ Phase 4 complete. Option B ships.

