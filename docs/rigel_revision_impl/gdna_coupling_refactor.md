# gDNA Coupling Refactor: T+N+2 → T+N+1

## Overview

Revert from the T+N+2 gDNA architecture (two strand-split gDNA components
coupled via a Beta prior in the M-step) to a T+N+1 architecture (single gDNA
component with an explicit `log(0.5)` strand probability baked into the
per-fragment likelihood).

**Zero backwards compatibility**: all deprecated code, parameters, and
documentation referring to the T+N+2 model are to be **deleted entirely**, not
deprecated or feature-gated.

---

## Rationale

### The Problem

In the T+N+2 model, gDNA fragments have **no strand term** in their
per-fragment log-likelihood (effectively `p(strand|gDNA) = 1.0`), while
mRNA/nRNA fragments incur a strand penalty (`log(SS)` for sense, `log(1-SS)`
for antisense, where SS ≈ 0.9).  This gives gDNA a **10× per-fragment
advantage** on antisense fragments compared to nRNA (`1.0` vs `0.1`).

The Beta(κ/2, κ/2) coupling in the M-step was supposed to enforce 50/50 strand
balance on gDNA, but it operates at the *aggregate* level (not per-fragment)
and fails catastrophically at achievable κ values.  Sweep data shows
+107.9% gDNA overestimation and −65.3% nRNA underestimation.

### The Fix

Add `log(0.5)` (≈ −0.693) to the gDNA per-fragment log-likelihood.  This
makes gDNA pay `p(strand|gDNA) = 0.5` per fragment, reducing the antisense
advantage from 10× to 5× (`0.5` vs `0.1`).  At realistic θ values
(nRNA >> gDNA), the abundance ratio naturally overcomes this residual
advantage.

### Mathematical Equivalence

At the true fixed point, the T+N+2 (κ = ∞), dual coupling, and T+N+1 +
log(0.5) models all produce **algebraically identical** E-step posteriors.
T+N+1 is the simplest implementation; no coupling block is needed.

---

## Component Layout Change

### Before (T+N+2)

```
Index:   [0, n_t)              → mRNA (one per transcript)
         [n_t, n_t + n_nrna)   → nRNA (one per unique nRNA span)
         [n_t + n_nrna]        → g_pos (positive-strand gDNA)
         [n_t + n_nrna + 1]    → g_neg (negative-strand gDNA)

Total:   n_t + n_nrna + 2
```

### After (T+N+1)

```
Index:   [0, n_t)              → mRNA (one per transcript)
         [n_t, n_t + n_nrna)   → nRNA (one per unique nRNA span)
         [n_t + n_nrna]        → gdna (single gDNA component)

Total:   n_t + n_nrna + 1
```

---

## File-by-File Changes

### 1. `src/rigel/native/scoring.cpp` — Add `log(0.5)` to gDNA likelihood

This is the core fix.  gDNA log-likelihood must include a strand term.

**Define constant** (near top of file, alongside other constants):

```cpp
static constexpr double LOG_HALF = -0.6931471805599453;  // log(0.5)
```

**Multi-hit path** (line ~800, inside `score_multimapper_units()`):

```cpp
// BEFORE:
double gdna_ll_val =
    gdna_fl + gdna_log_sp;

// AFTER:
double gdna_ll_val =
    gdna_fl + gdna_log_sp + LOG_HALF;
```

**Single-hit path** (line ~1245, inside `score_single_hit_units()`):

```cpp
// BEFORE:
st.gdna_ll[st.unit_cur] =
        gdna_fl + gdna_log_sp;

// AFTER:
st.gdna_ll[st.unit_cur] =
        gdna_fl + gdna_log_sp + LOG_HALF;
```

No other changes needed in scoring.cpp.  mRNA/nRNA strand scoring
(`log_p_sense_`, `log_p_antisense_`) remains unchanged.

---

### 2. `src/rigel/native/em_solver.cpp` — Remove coupling and dual gDNA

This file requires the most changes.  Every reference to `g_pos`/`g_neg`,
`gdna_pos_idx`/`gdna_neg_idx`, and `strand_symmetry_kappa` must be updated.

#### 2a. `map_em_step()` (lines 516–580)

**Remove `strand_symmetry_kappa` parameter**.  Delete the entire Beta coupling
block (lines 555–575).

```cpp
// BEFORE (line 527):
    double        strand_symmetry_kappa = 2.0)

// AFTER:
    )   // strand_symmetry_kappa parameter removed
```

**Delete the entire coupling block** (lines 555–575):

```cpp
// DELETE THIS ENTIRE BLOCK:
    // Couple g_pos and g_neg via symmetric Beta(κ/2, κ/2) prior.
    // g_pos is always nc-2, g_neg is always nc-1.
    // kappa <= 2.0 disables the coupling (fast path).
    if (strand_symmetry_kappa > 2.0 && n_components >= 2) {
        int g_pos = n_components - 2;
        int g_neg = n_components - 1;
        double R_pos = theta_new[g_pos];
        double R_neg = theta_new[g_neg];
        double R_g = R_pos + R_neg;
        if (R_g > 0.0) {
            double phi = (R_pos + strand_symmetry_kappa / 2.0 - 1.0)
                       / (R_g + strand_symmetry_kappa - 2.0);
            phi = std::clamp(phi, 0.0, 1.0);
            double new_pos = phi * R_g;
            double new_neg = (1.0 - phi) * R_g;
            total += (new_pos - theta_new[g_pos]) + (new_neg - theta_new[g_neg]);
            theta_new[g_pos] = new_pos;
            theta_new[g_neg] = new_neg;
        }
    }
```

#### 2b. `run_squarem()` (line 716)

**Remove `strand_symmetry_kappa` parameter**:

```cpp
// BEFORE (line 716):
    double        strand_symmetry_kappa = 2.0)

// AFTER:
    )   // strand_symmetry_kappa parameter removed
```

**Remove `strand_symmetry_kappa` from all `map_em_step()` calls** inside
`run_squarem()`.  There are **4 call sites** (lines 832, 838, 875, 959):

For each call, remove the trailing `,\n                        strand_symmetry_kappa);`
argument.

#### 2c. Pruning protection (lines 914–916)

```cpp
// BEFORE:
        // Never prune gDNA (last two components: g_pos and g_neg)
        if (nc >= 2) prune_mask[nc - 2] = false;
        prune_mask[nc - 1] = false;

// AFTER:
        // Never prune gDNA (last component)
        prune_mask[nc - 1] = false;
```

#### 2d. `LocusSubProblem` struct (lines 1147–1155)

```cpp
// BEFORE:
    int n_components;  // n_t + n_nrna + 2
    int gdna_pos_idx;  // = n_t + n_nrna     (positive-strand gDNA)
    int gdna_neg_idx;  // = n_t + n_nrna + 1 (negative-strand gDNA)

// AFTER:
    int n_components;  // n_t + n_nrna + 1
    int gdna_idx;      // = n_t + n_nrna     (single gDNA component)
```

#### 2e. `extract_locus_sub_problem()` (lines 1248–1250)

```cpp
// BEFORE:
    sub.n_components = n_t + n_nrna + 2;
    sub.gdna_pos_idx = n_t + n_nrna;
    sub.gdna_neg_idx = n_t + n_nrna + 1;

// AFTER:
    sub.n_components = n_t + n_nrna + 1;
    sub.gdna_idx = n_t + n_nrna;
```

#### 2f. gDNA strand routing in `extract_locus_sub_problem()` (lines 1383–1389)

Remove the strand-based routing.  All unspliced units point to the single
`gdna_idx`.

```cpp
// BEFORE:
            // Strand bit from locus_ct_arr: count_col = stype*2 + is_anti
            bool is_anti = (sub.locus_ct_arr[ui] % 2) != 0;
            int32_t gdna_comp = is_anti ? sub.gdna_neg_idx : sub.gdna_pos_idx;

// AFTER:
            int32_t gdna_comp = sub.gdna_idx;
```

#### 2g. Bias profiles in `extract_locus_sub_problem()` (lines 1438–1439)

```cpp
// BEFORE:
    sub.bias_profiles[sub.gdna_pos_idx] = locus_span;
    sub.bias_profiles[sub.gdna_neg_idx] = locus_span;

// AFTER:
    sub.bias_profiles[sub.gdna_idx] = locus_span;
```

#### 2h. Prior zeroing in `extract_locus_sub_problem()` (lines 1446–1447)

```cpp
// BEFORE:
    if (gdna_init == 0.0) {
        sub.prior[sub.gdna_pos_idx] = 0.0;
        sub.prior[sub.gdna_neg_idx] = 0.0;
    }

// AFTER:
    if (gdna_init == 0.0) {
        sub.prior[sub.gdna_idx] = 0.0;
    }
```

#### 2i. `assign_posteriors()` (lines 1512–1614)

Replace dual gDNA index tracking with single index.

```cpp
// BEFORE:
    int gdna_pos = sub.gdna_pos_idx;
    int gdna_neg = sub.gdna_neg_idx;

// AFTER:
    int gdna = sub.gdna_idx;
```

Update the nRNA check boundary (line ~1580):

```cpp
// BEFORE:
            } else if (comp < gdna_pos) {

// AFTER:
            } else if (comp < gdna) {
```

Update the gDNA locus attribution check (line ~1607):

```cpp
// BEFORE:
            if (c == gdna_pos || c == gdna_neg) {
                gdna_unit_sum += posteriors[j];
            }

// AFTER:
            if (c == gdna) {
                gdna_unit_sum += posteriors[j];
            }
```

#### 2j. `batch_locus_em()` (line 1690)

**Remove `strand_symmetry_kappa` parameter**:

```cpp
// BEFORE:
    double strand_symmetry_kappa = 2.0)

// AFTER:
    )   // strand_symmetry_kappa parameter removed
```

**Remove `strand_symmetry_kappa` from `run_squarem()` call** (line 1919):

```cpp
// BEFORE:
            EMResult result = run_squarem(
                ec_data, log_eff_len.data(),
                sub.unambig_totals.data(),
                prior.data(),
                theta_init.data(),
                nc, max_iterations, convergence_delta,
                use_vbem, prune_threshold,
                estep_thr, pool,
                strand_symmetry_kappa);

// AFTER:
            EMResult result = run_squarem(
                ec_data, log_eff_len.data(),
                sub.unambig_totals.data(),
                prior.data(),
                theta_init.data(),
                nc, max_iterations, convergence_delta,
                use_vbem, prune_threshold,
                estep_thr, pool);
```

#### 2k. nanobind module binding (line ~2363)

```cpp
// BEFORE:
          nb::arg("strand_symmetry_kappa") = 2.0,

// AFTER:
          // Remove this line entirely
```

#### 2l. Comment/docstring cleanup in `extract_locus_sub_problem()`

Update the layout comment (line ~1189):

```cpp
// BEFORE:
// Component layout: [0,n_t) mRNA + [n_t, n_t+n_nrna) nRNA + [n_t+n_nrna] g_pos + [n_t+n_nrna+1] g_neg

// AFTER:
// Component layout: [0,n_t) mRNA + [n_t, n_t+n_nrna) nRNA + [n_t+n_nrna] gdna
```

Update the block heading comment above `map_em_step()` (lines ~508–515):

```cpp
// BEFORE:
// Hierarchical MAP-EM step: theta → theta_new
// ...
// Plain MAP-EM step with symmetric Beta coupling for gDNA strand balance.
// g_pos and g_neg are the last two components (nc-2, nc-1).
// After the standard M-step accumulation, the coupling redistributes
// gDNA mass between g_pos and g_neg via a Beta(κ/2, κ/2) prior on the
// strand balance φ = θ_{g_pos} / (θ_{g_pos} + θ_{g_neg}).

// AFTER:
// MAP-EM step: theta → theta_new
//
// Plain MAP-EM with per-component Dirichlet prior.
// gDNA strand symmetry is enforced via per-fragment log(0.5) in the likelihood
// (added at scoring time), so no M-step coupling is needed.
```

---

### 3. `src/rigel/locus.py` — Simplify component layout

#### 3a. Component setup (lines 181–183)

```python
# BEFORE:
    gdna_pos_idx = n_t + n_nrna  # positive-strand gDNA component index
    gdna_neg_idx = n_t + n_nrna + 1  # negative-strand gDNA component index
    n_components = n_t + n_nrna + 2

# AFTER:
    gdna_idx = n_t + n_nrna  # single gDNA component index
    n_components = n_t + n_nrna + 1
```

#### 3b. gDNA strand routing (lines 302–307)

Remove strand-based routing; all gDNA units map to the single `gdna_idx`.

```python
# BEFORE:
        # Strand-route: sense → g_pos, antisense → g_neg
        is_anti = em_data.locus_count_cols[locus.unit_indices][valid_gdna] % 2
        gdna_lidx_arr = np.where(
            is_anti,
            gdna_neg_idx,
            gdna_pos_idx,
        ).astype(np.int32)

# AFTER:
        gdna_lidx_arr = np.full(n_gdna, gdna_idx, dtype=np.int32)
```

#### 3c. Bias profiles (lines 394–395)

```python
# BEFORE:
    bias_profiles[gdna_pos_idx] = int(locus_span)
    bias_profiles[gdna_neg_idx] = int(locus_span)

# AFTER:
    bias_profiles[gdna_idx] = int(locus_span)
```

#### 3d. Prior zeroing (lines 404–405)

```python
# BEFORE:
    if gdna_init == 0.0:
        prior[gdna_pos_idx] = 0.0
        prior[gdna_neg_idx] = 0.0

# AFTER:
    if gdna_init == 0.0:
        prior[gdna_idx] = 0.0
```

#### 3e. Component layout docstrings (lines 130–140, 130–140 in function)

Update both the class docstring and function docstring:

```python
# BEFORE:
        [n_t + n_nrna]          — g_pos (positive-strand gDNA)
        [n_t + n_nrna + 1]     — g_neg (negative-strand gDNA)

    Total components = n_transcripts + n_nrna + 2.

# AFTER:
        [n_t + n_nrna]          — gdna (single gDNA component)

    Total components = n_transcripts + n_nrna + 1.
```

And in `build_locus_em_data()`:

```python
# BEFORE:
        [n_t + n_nrna]          - g_pos (positive-strand gDNA)
        [n_t + n_nrna + 1]     - g_neg (negative-strand gDNA)

    Only UNSPLICED units get a gDNA candidate (g_pos or g_neg, depending
    on strand).

# AFTER:
        [n_t + n_nrna]          - gdna (single gDNA component)

    Only UNSPLICED units get a gDNA candidate.
```

#### 3f. Bias profiles comment block (lines ~386–392)

```python
# BEFORE:
    #   [n_t + n_nrna]         → g_pos: length = locus span
    #   [n_t + n_nrna + 1]    → g_neg: length = locus span

# AFTER:
    #   [n_t + n_nrna]         → gdna: length = locus span
```

---

### 4. `src/rigel/config.py` — Delete `strand_symmetry_kappa`

**Delete the entire field and its docstring** (lines 92–97):

```python
# DELETE:
    strand_symmetry_kappa: float = 6.0
    """Strand symmetry penalty strength for gDNA in the M-step.
    Effective κ scales by strand specificity: κ_eff = κ · (2·SS − 1)².
    Set ≤ 2.0 to disable the penalty entirely.
    """
```

Also delete the section header comment above it (line ~91):

```python
# DELETE:
    # -- gDNA strand symmetry penalty knobs --
```

---

### 5. `src/rigel/estimator.py` — Remove kappa from `batch_locus_em` call

**Remove the `strand_symmetry_kappa` argument** (line 393):

```python
# BEFORE:
            self.em_config.n_threads,
            self.em_config.strand_symmetry_kappa,
        )

# AFTER:
            self.em_config.n_threads,
        )
```

---

### 6. `src/rigel/cli.py` — Delete CLI parameter

#### 6a. Delete `_ParamSpec` entry (line 423):

```python
# DELETE:
    _ParamSpec("strand_symmetry_kappa", "em.strand_symmetry_kappa"),
```

#### 6b. Delete argparse argument (lines 824–829):

```python
# DELETE:
    adv.add_argument(
        "--strand-symmetry-kappa", dest="strand_symmetry_kappa",
        type=float, default=None,
        help="Strand symmetry penalty κ for gDNA in the M-step. "
             "Effective κ scales by strand specificity: κ_eff = κ·(2·SS−1)². "
             "Set ≤ 2.0 to disable (default: 6.0).",
    )
```

---

### 7. `src/rigel/scored_fragments.py` — Update docstrings

Update the `LocusEMInput` docstring (lines 136–139):

```python
# BEFORE:
        [n_t + n_nrna]          — g_pos (positive-strand gDNA)
        [n_t + n_nrna + 1]     — g_neg (negative-strand gDNA)

    Total components = n_transcripts + n_nrna + 2.

# AFTER:
        [n_t + n_nrna]          — gdna (single gDNA component)

    Total components = n_transcripts + n_nrna + 1.
```

---

### 8. Tests

#### 8a. `tests/test_gdna.py` (line ~925)

Delete or update the assertion checking `strand_symmetry_kappa == 6.0`:

```python
# DELETE:
    assert config.strand_symmetry_kappa == 6.0
```

#### 8b. `tests/test_em_impl.py`

Many tests already use `n_components = n_t + n_nrna + 1`.  **Search** for any
tests still using `+ 2` and update.  Also search for any tests passing
`strand_symmetry_kappa` to the C++ `batch_locus_em()` or setting up g_pos/g_neg
component arrays.

Key things to check:
- Tests that construct `prior`, `unambig_totals`, `bias_profiles` arrays
  with `n_t + n_nrna + 2` length → change to `+ 1`
- Tests that assign separate g_pos/g_neg component indices → use single
  `gdna_idx = n_t + n_nrna`
- Tests that pass `strand_symmetry_kappa` kwarg → remove it
- Golden output test baselines may need regeneration if counts shift

#### 8c. Golden outputs (`tests/golden/`)

After implementation, run:

```bash
pytest tests/test_golden_output.py --update-golden
```

Verify that the updated golden outputs show reduced gDNA and increased nRNA
as expected from the model change.

#### 8d. Any test constructing `EMConfig`

Search for `strand_symmetry_kappa` in all test files and remove.

---

### 9. Scripts

**All scripts** that reference `strand_symmetry_kappa` must have the parameter
removed.  Affected files (non-exhaustive — grep the `scripts/` directory):

- `scripts/ablation_comprehensive.yaml`
- `scripts/analyze_ablation_deep.py`
- `scripts/analyze_ablation_synthesis.py`
- `scripts/analyze_ablation_v2.py`
- `scripts/analyze_ablation.py`
- `scripts/analyze_sweep.py`
- `scripts/benchmark_example.yaml`
- `scripts/benchmark.py`
- `scripts/nrna_sweep_config.yaml`
- `scripts/sim_example.yaml`
- `scripts/sim.py`
- `scripts/synthetic_sim_sweep.py`
- Any other YAML configs or Python scripts

**Action**: Run `grep -rn strand_symmetry_kappa scripts/` and remove every
occurrence.

---

### 10. Documentation

Update or delete references to g_pos/g_neg, Beta coupling, and
strand_symmetry_kappa in:

- `docs/METHODS.md`
- `docs/parameters.md`
- `docs/MANUAL.md`
- `docs/TODO.md`
- `docs/em_dynamics_analysis.md` (if present)
- `docs/CODE_PATH.md`
- `CLAUDE.md` (component layout table if it references g_pos/g_neg)
- Any other docs referencing the old architecture

**Action**: Run `grep -rn 'g_pos\|g_neg\|strand_symmetry_kappa\|Beta.*coupling' docs/`
and update every occurrence.

---

## Implementation Order

Execute changes in this order to maintain a compilable/testable state at
each step:

1. **scoring.cpp**: Add `LOG_HALF` constant and the two `+ LOG_HALF` terms.
   This change is backward-compatible — the old T+N+2 EM still runs, just
   with slightly different likelihoods.

2. **em_solver.cpp**: All changes in one commit:
   - Update `LocusSubProblem` struct
   - Update `extract_locus_sub_problem()`
   - Remove coupling block from `map_em_step()`
   - Remove `strand_symmetry_kappa` from all function signatures
   - Update `assign_posteriors()`
   - Update pruning protection
   - Update nanobind binding

3. **locus.py**: Simplify component layout, remove strand routing.

4. **config.py**: Delete `strand_symmetry_kappa` field.

5. **estimator.py**: Remove kappa from `batch_locus_em()` call.

6. **cli.py**: Delete CLI parameter and `_ParamSpec`.

7. **scored_fragments.py**: Update docstrings.

8. **Tests**: Fix all test files, update golden outputs.

9. **Scripts and Docs**: Clean up all references.

**Recompile after step 2**: `pip install --no-build-isolation -e .`

**Run full test suite after step 8**: `pytest tests/ -v`

---

## Verification Checklist

After implementation, verify:

- [ ] `pip install --no-build-isolation -e .` compiles cleanly
- [ ] `pytest tests/ -v` — all tests pass
- [ ] `grep -rn 'g_pos\|g_neg\|gdna_pos\|gdna_neg' src/` — zero matches
- [ ] `grep -rn strand_symmetry_kappa src/ tests/ scripts/` — zero matches
- [ ] `grep -rn 'n_nrna + 2' src/ tests/` — zero matches
- [ ] `ruff check src/ tests/` — no lint errors
- [ ] `ruff format src/ tests/` — clean
- [ ] Golden output regenerated and inspected for expected gDNA reduction
- [ ] Component count in a sample run is `n_t + n_nrna + 1`

---

## Summary of Removed Artifacts

| Artifact | Location | Action |
|----------|----------|--------|
| `strand_symmetry_kappa` field | `config.py:92` | Delete |
| `--strand-symmetry-kappa` CLI arg | `cli.py:824-829` | Delete |
| `_ParamSpec("strand_symmetry_kappa", ...)` | `cli.py:423` | Delete |
| `self.em_config.strand_symmetry_kappa` | `estimator.py:393` | Remove from call |
| Beta coupling block | `em_solver.cpp:555-575` | Delete |
| `strand_symmetry_kappa` parameter (4 functions) | `em_solver.cpp` | Remove from signatures |
| `gdna_pos_idx` / `gdna_neg_idx` struct fields | `em_solver.cpp:1149-1150` | Replace with `gdna_idx` |
| `gdna_pos_idx` / `gdna_neg_idx` variables | `locus.py:181-182` | Replace with `gdna_idx` |
| Strand-routing `np.where(is_anti, ...)` | `locus.py:302-307` | Replace with `np.full` |
| `nb::arg("strand_symmetry_kappa")` binding | `em_solver.cpp:2363` | Delete |
| g_pos/g_neg docstrings | `scored_fragments.py`, `locus.py` | Update |
| All script/YAML references | `scripts/**` | Delete |
| All doc references | `docs/**` | Update |
