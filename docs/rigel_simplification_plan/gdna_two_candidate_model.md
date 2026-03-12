# gDNA Two-Candidate Model (T+N+2)

## Summary

Replace the single collapsed gDNA component (T+N+1) with two explicit
strand-specific gDNA components (T+N+2): `g_pos` (positive-strand DNA) and
`g_neg` (negative-strand DNA). This eliminates the systematic per-fragment
LOG_HALF penalty that causes ~40% gDNA underestimation at SS ≥ 0.75.

**Component layout change:**

```
Before: [0, n_t)  mRNA  |  [n_t, n_t+n_nrna)  nRNA  |  [n_t+n_nrna]  gDNA
After:  [0, n_t)  mRNA  |  [n_t, n_t+n_nrna)  nRNA  |  [n_t+n_nrna]  g_pos  |  [n_t+n_nrna+1]  g_neg
```

**Net effect on codebase**: Simplification. Removes the M-step targeted-excess
penalty (~60 lines), removes `gdna_col`/`gdna_is_anti` from equivalence
classes (~40 lines), and removes `strand_specificity` parameter threading
(~20 references). Adds ~10 lines for M-step theta coupling.

---

## Problem Statement

### The LOG_HALF Disadvantage

In `scoring.cpp:1245`, gDNA log-likelihood is computed as:

```cpp
gdna_ll = LOG_HALF + gdna_fl + gdna_log_sp   // LOG_HALF = log(0.5) = -0.693
```

Meanwhile, mRNA/nRNA use the data-driven strand model (`scoring.cpp:1025`):

```cpp
log_strand = same ? log_p_sense_ : log_p_antisense_
```

At SS=0.9: `log_p_sense = log(0.9) = -0.105`. For each sense-strand fragment,
gDNA is at a **0.588 nat disadvantage** vs. mRNA/nRNA. Over N fragments, gDNA
mass is suppressed by a factor of ~exp(0.588 × N_sense - 1.610 × N_anti),
which becomes astronomically large for even moderate fragment counts.

### Why Phase 3 Failed

Phase 3 scaled the M-step strand penalty by κ_eff = κ · (2·SS − 1)². This
had no effect because the penalty operates on E-step responsibilities that are
ALREADY corrupt: the LOG_HALF disadvantage means gDNA gets near-zero
responsibility for sense fragments before the M-step ever runs.

### T+N+2 Solution

Two explicit gDNA strand components with **deterministic** strand likelihoods:

| | Sense fragment | Antisense fragment |
|---|---|---|
| g_pos | log(1) = 0 | −∞ |
| g_neg | −∞ | log(1) = 0 |

The correct-strand gDNA component gets **no strand penalty** (log(1) = 0),
which is actually slightly better than mRNA's log(0.9) = −0.105 at SS=0.9.
Competition is then decided fairly by fragment length, abundance, and
overhang — not by the strand model asymmetry.

---

## Code Changes

### Overview

| # | File | Action | Lines affected |
|---|---|---|---|
| 1 | `src/rigel/native/scoring.cpp` | Remove LOG_HALF from gDNA LL | ~10 lines |
| 2 | `src/rigel/native/em_solver.cpp` | T+N+2 component layout, remove EC strand tracking, replace M-step penalty with theta coupling | ~250 lines |
| 3 | `src/rigel/native/constants.h` | Remove LOG_HALF constant | 1 line |
| 4 | `src/rigel/locus.py` | T+N+2 component layout, strand-routed gDNA injection | ~30 lines |
| 5 | `src/rigel/scored_fragments.py` | Update component count comment | 1 line |
| 6 | `src/rigel/estimator.py` | Collapse g_pos + g_neg in output | ~10 lines |
| 7 | `src/rigel/pipeline.py` | Remove strand_specificity parameter | ~5 lines |
| 8 | `src/rigel/config.py` | Remove strand_symmetry_kappa | ~5 lines |
| 9 | `src/rigel/cli.py` | Remove --strand-symmetry-kappa CLI flag | ~5 lines |
| 10 | `tests/` | Update all gDNA tests, regenerate golden outputs | ~200 lines |

---

### Step 1: `src/rigel/native/scoring.cpp` — Remove LOG_HALF from gDNA LL

**Goal**: gDNA log-likelihoods no longer include any strand term. The strand
routing is handled by component selection (g_pos vs g_neg) at locus
construction time.

#### 1a. Single-mapper gDNA (lines 1241–1249)

**Before:**
```cpp
st.gdna_ll[st.unit_cur] = LOG_HALF + gdna_fl + gdna_log_sp;
```

**After:**
```cpp
st.gdna_ll[st.unit_cur] = gdna_fl + gdna_log_sp;
```

#### 1b. Multi-mapper gDNA (lines 797–801)

**Before:**
```cpp
double gdna_ll_val = LOG_HALF + gdna_fl + gdna_log_sp;
```

**After:**
```cpp
double gdna_ll_val = gdna_fl + gdna_log_sp;
```

#### 1c. Remove `using rigel::LOG_HALF` import (line 30)

Only if LOG_HALF is not used elsewhere in scoring.cpp. Search for all
LOG_HALF references in the file. The mRNA/nRNA strand scoring path uses
LOG_HALF as a fallback when `has_strand` is false (lines 534, 1032):
```cpp
log_strand = LOG_HALF;  // unstranded fallback
```
**Decision**: Keep `LOG_HALF` import if the unstranded fallback still uses it.
The unstranded fallback (`has_strand == false`) sets `log_strand = LOG_HALF`
for mRNA/nRNA candidates. This is correct: for unstranded libraries, all
components including mRNA/nRNA should use LOG_HALF for their strand term.
gDNA does not need LOG_HALF because strand routing is deterministic.

**If LOG_HALF is only used in the gDNA path**: remove the import and the
constant from `constants.h`.

**If LOG_HALF is used in mRNA/nRNA unstranded fallback**: keep the constant
but remove it from gDNA scoring.

→ **Check at implementation time.** Preliminary scan shows LOG_HALF is used
in the unstranded fallback for mRNA/nRNA (lines 534, 1032), so we keep the
constant and only remove it from gDNA scoring.

---

### Step 2: `src/rigel/native/em_solver.cpp` — Core EM Changes

This is the largest change. Broken into sub-steps:

#### 2a. Remove `gdna_col` and `gdna_is_anti` from EmEquivClass (lines 100–115)

**Before:**
```cpp
struct EmEquivClass {
    std::vector<int32_t> comp_idx;
    std::vector<double>  ll_flat;
    int gdna_col;                       // REMOVE
    std::vector<uint8_t> gdna_is_anti;  // REMOVE
    int n, k;
    // ...
};
```

**After:**
```cpp
struct EmEquivClass {
    std::vector<int32_t> comp_idx;
    std::vector<double>  ll_flat;
    int n, k;
    // ...
};
```

#### 2b. Simplify `build_equiv_classes()` (lines 135–280)

Remove all gDNA-specific tracking:

- **Remove parameters**: `const uint8_t* locus_ct_arr = nullptr` (line 138),
  `int gdna_idx = -1` (line 139)
- **Remove gdna_col search** (lines 172–178): No longer needed. g_pos and
  g_neg are regular components.
- **Remove gdna_is_anti population** (lines 180–194): Strand info is encoded
  in which component the unit maps to (g_pos or g_neg), not in a side array.
- **Remove gdna_is_anti reordering** (lines 252–271): Part of the
  deterministic sort; remove the gdna_is_anti row-swap code.

The function signature simplifies to:
```cpp
static std::vector<EmEquivClass> build_equiv_classes(
    const int64_t* offsets,
    const int32_t* t_indices,
    const double*  log_liks,
    const double*  coverage_wts,
    int            n_units)
```

#### 2c. Replace M-step penalty in `map_em_step()` (lines 537–640)

**Remove entirely** (lines 555–556, 589–638):
- `strand_symmetry_kappa` parameter
- `strand_specificity` parameter
- `int gi = n_components - 1` (gDNA-is-last assumption)
- The entire kappa_eff computation
- The entire e_sense/e_anti accumulation loop
- The entire targeted-excess penalty (e_sym, e_excess, w_sym, penalized_em)

**Replace with** theta coupling (~10 lines) after the standard M-step
normalization. The coupling enforces a symmetric Beta prior on the gDNA
strand balance φ = θ_{g_pos} / (θ_{g_pos} + θ_{g_neg}):

```cpp
// Couple g_pos and g_neg via symmetric Beta(κ/2, κ/2) prior
// g_pos is always nc-2, g_neg is always nc-1
int g_pos = n_components - 2;
int g_neg = n_components - 1;
double R_pos = theta_new[g_pos];  // unnormalized (pre-norm)
double R_neg = theta_new[g_neg];
double R_g = R_pos + R_neg;
if (R_g > 0.0 && strand_symmetry_kappa > 2.0) {
    double phi = (R_pos + strand_symmetry_kappa / 2.0 - 1.0)
               / (R_g + strand_symmetry_kappa - 2.0);
    phi = std::clamp(phi, 0.0, 1.0);
    theta_new[g_pos] = phi * R_g;
    theta_new[g_neg] = (1.0 - phi) * R_g;
}
```

**Note**: This coupling must happen BEFORE normalization. The current M-step
first accumulates `theta_new[i] = unambig + em + prior`, THEN applies the
penalty, THEN normalizes. The coupling replaces the penalty step in the same
position.

**New parameter**: `strand_symmetry_kappa` is kept (controls coupling
strength). `strand_specificity` is **removed** — the scaling logic
κ_eff = κ·(2·SS−1)² is gone because the two-component model inherently
handles strand specificity correctly.

**Updated `map_em_step` signature:**
```cpp
static void map_em_step(
    const double* theta,
    const std::vector<EmEquivClass>& ec_data,
    const double* log_eff_len,
    const double* unambig_totals,
    const double* prior,
    double*       em_totals,
    double*       theta_new,
    int           n_components,
    int           estep_threads = 1,
    rigel::EStepThreadPool* pool = nullptr,
    double        strand_symmetry_kappa = 2.0)
```

#### 2d. Update `run_squarem()` (lines 760–1030)

- Remove `strand_specificity` parameter (line 775).
- Update all `map_em_step(...)` call sites to remove `strand_specificity`
  argument (lines 891, 897, 934, 1017).
- Update "never prune gDNA" logic (line 975): Currently
  `prune_mask[nc - 1] = false`. With T+N+2, two components are gDNA:
  ```cpp
  prune_mask[nc - 1] = false;  // g_neg
  prune_mask[nc - 2] = false;  // g_pos
  ```

#### 2e. Update `LocusSubProblem` struct (lines 1200–1240)

**Before:**
```cpp
int n_components;  // n_t + n_nrna + 1
int gdna_idx;      // = n_t + n_nrna
```

**After:**
```cpp
int n_components;   // n_t + n_nrna + 2
int gdna_pos_idx;   // = n_t + n_nrna
int gdna_neg_idx;   // = n_t + n_nrna + 1
```

Update the component layout comment:
```
// [0,n_t) mRNA + [n_t, n_t+n_nrna) nRNA + [n_t+n_nrna] g_pos + [n_t+n_nrna+1] g_neg
```

#### 2f. Update `extract_locus_sub_problem()` (lines 1246–1510)

- Line 1304: `sub.gdna_idx = n_t + n_nrna` → `sub.gdna_pos_idx = n_t + n_nrna; sub.gdna_neg_idx = n_t + n_nrna + 1`
- Line 1305: `sub.n_components = n_t + n_nrna + 1` → `sub.n_components = n_t + n_nrna + 2`

**gDNA candidate injection** (lines 1441–1453):

Currently each unspliced unit gets a single gDNA candidate at `gdna_comp =
sub.gdna_idx`. In T+N+2, each unspliced unit gets ONE gDNA candidate at
EITHER `gdna_pos_idx` or `gdna_neg_idx`, depending on strand:

```cpp
// count_col encoding: stype * 2 + is_anti
// Strand bit = locus_ct_arr[ui] % 2: 0 = sense, 1 = antisense
bool is_anti = (sub.locus_ct_arr[ui] % 2) != 0;
int32_t gdna_comp = is_anti ? sub.gdna_neg_idx : sub.gdna_pos_idx;
```

The existing `locus_ct_arr` already stores the strand information per unit,
so no new data flow is needed. The strand bit is extracted the same way as
the current `gdna_is_anti` logic in `build_equiv_classes`.

**Bias profiles** (line 1492):
```cpp
// Before
sub.bias_profiles[sub.gdna_idx] = locus_span;
// After
sub.bias_profiles[sub.gdna_pos_idx] = locus_span;
sub.bias_profiles[sub.gdna_neg_idx] = locus_span;
```

**Prior gating** (lines 1497–1499):
```cpp
// Before
if (gdna_init == 0.0) { sub.prior[sub.gdna_idx] = 0.0; }
// After
if (gdna_init == 0.0) {
    sub.prior[sub.gdna_pos_idx] = 0.0;
    sub.prior[sub.gdna_neg_idx] = 0.0;
}
```

**Eligible array**:
Both g_pos and g_neg must be eligible (nonzero) to participate.

#### 2g. Update `assign_posteriors()` (lines 1545–1680)

The posterior scatter function checks `comp < n_t` → mRNA, `comp < gdna_idx`
→ nRNA, else → gDNA. With T+N+2:

```cpp
int gdna_pos = sub.gdna_pos_idx;
int gdna_neg = sub.gdna_neg_idx;

if (comp < n_t) {
    // mRNA — unchanged
} else if (comp < gdna_pos) {
    // nRNA — unchanged
} else {
    // gDNA (either g_pos or g_neg)
    gdna_total += p;
}
```

The gDNA locus attribution (lines 1660–1672) currently checks
`sub.t_indices[s + j] == gdna_idx`. Update to:
```cpp
if (sub.t_indices[s + j] == gdna_pos || sub.t_indices[s + j] == gdna_neg) {
    gdna_unit_sum += posteriors[j];
}
```

#### 2h. Update `batch_locus_em()` entry point (lines 1683–2091)

- In the `build_equiv_classes` call (line ~1945): Remove `locus_ct_arr` and
  `gdna_idx` arguments.
- Remove `strand_specificity` parameter from the function signature and
  nanobind binding.
- `strand_symmetry_kappa` remains (controls coupling).

#### 2i. Update nanobind module definition (lines 2370–2430)

- Remove `nb::arg("strand_specificity") = 0.5` from `batch_locus_em` binding.
- Keep `nb::arg("strand_symmetry_kappa") = 2.0`.
- Update docstring to reflect T+N+2.

#### 2j. Update `run_locus_em_native()` entry point (lines 1035–1185)

This simpler path (used in some tests) does NOT pass gdna_idx or
locus_ct_arr to build_equiv_classes — it already uses the simplified call.
No `strand_symmetry_kappa` is passed to `run_squarem` (uses defaults).

**Changes needed**:
- No build_equiv_classes signature changes affect this path (it already
  uses the minimal signature).
- The `run_squarem` call (line 1151) uses default kappa=2.0 which means
  coupling is disabled (kappa <= 2.0 is the fast-path). This is fine for
  tests that use this entry point.

---

### Step 3: `src/rigel/native/constants.h` — Cleanup

- **Keep** `LOG_HALF`: still used in mRNA/nRNA unstranded fallback.
- No change needed unless a full audit reveals LOG_HALF is only used in gDNA.

---

### Step 4: `src/rigel/locus.py` — Python Locus Construction

The Python locus construction (`build_locus_em_data`) is the fallback/test
path. It must mirror the C++ changes.

#### 4a. Component layout (line 180–181)

**Before:**
```python
gdna_idx = n_t + n_nrna
n_components = n_t + n_nrna + 1
```

**After:**
```python
gdna_pos_idx = n_t + n_nrna
gdna_neg_idx = n_t + n_nrna + 1
n_components = n_t + n_nrna + 2
```

#### 4b. gDNA CSR injection (lines 298–317)

Currently, all unspliced units get gDNA candidates at `gdna_idx`. In T+N+2,
each unit gets a candidate at either `gdna_pos_idx` or `gdna_neg_idx`
depending on strand.

The strand bit is available from `em_data.locus_count_cols[locus.unit_indices]`:
```python
is_anti = em_data.locus_count_cols[locus.unit_indices] % 2  # 0=sense, 1=anti
```

For each valid gDNA unit:
```python
gdna_lidx_arr = np.where(
    is_anti[valid_gdna],
    gdna_neg_idx,
    gdna_pos_idx,
).astype(np.int32)
```

#### 4c. Effective lengths (line 376)

**Before:**
```python
# [n_t + n_nrna] → gDNA: length = locus span
```

**After:**
```python
# [n_t + n_nrna] → g_pos: length = locus span
# [n_t + n_nrna + 1] → g_neg: length = locus span
```

Both g_pos and g_neg get `locus_span`.

#### 4d. Prior gating (lines 392–394)

**Before:**
```python
if gdna_init == 0.0:
    prior[gdna_idx] = 0.0
```

**After:**
```python
if gdna_init == 0.0:
    prior[gdna_pos_idx] = 0.0
    prior[gdna_neg_idx] = 0.0
```

#### 4e. Update docstring and comments

Update the component layout docstring (lines 131–133) and all internal
comments that reference "+1" or "single gDNA component".

---

### Step 5: `src/rigel/scored_fragments.py` — Comment Update

Line 139: `"Total components = n_transcripts + n_nrna + 1"` → `"Total
components = n_transcripts + n_nrna + 2"`.

No structural changes — `gdna_log_liks` remains a per-unit array of
pre-computed values (without LOG_HALF now). The strand routing happens
downstream in locus construction.

---

### Step 6: `src/rigel/estimator.py` — Collapse g_pos + g_neg

The C++ `batch_locus_em` returns `(total_gdna_em, locus_mrna, locus_nrna,
locus_gdna)`. The `total_gdna_em` scalar already sums both g_pos and g_neg
(from `assign_posteriors`). The per-locus `locus_gdna` array also sums both.
No Python-side collapse is needed — the C++ output is already collapsed.

However, verify that `run_batch_locus_em()` correctly passes parameters:

- **Remove**: `strand_specificity` parameter (line 279) and its forwarding
  (line ~394).
- **Keep**: `strand_symmetry_kappa` forwarding from `em_config`.

The `gdna_em_count`, `gdna_total`, `gdna_contamination_rate` properties are
unchanged — they read from `_gdna_em_total` which is already the collapsed
sum.

---

### Step 7: `src/rigel/pipeline.py` — Remove strand_specificity

- In `_run_locus_em()` (line ~508): Remove `strand_specificity: float = 0.5`
  parameter.
- Remove the forwarding of `strand_specificity` to
  `estimator.run_batch_locus_em()`.
- Remove the computation/lookup of strand_specificity from strand_models in
  the calling code.

---

### Step 8: `src/rigel/config.py` — Remove strand_symmetry_kappa

**Remove** `strand_symmetry_kappa` from `EMConfig` (line 92). The
coupling strength κ is now a hardcoded constant in the C++ M-step. With
T+N+2, the coupling is simple enough that a fixed κ=6.0 (matching the
previous default) provides symmetric regularization without needing
per-library tuning.

**Rationale**: The κ_eff = κ · (2·SS−1)² scaling was added in Phase 3
specifically because a fixed κ was too strong at low SS. In T+N+2, the
strand handling is in the likelihood (deterministic 0/−∞ routing), not in
the prior/penalty. The coupling κ controls only how much the gDNA strand
ratio can deviate from 50:50 — a fundamentally different (and much simpler)
role than the old penalty. A fixed value is appropriate.

**Alternative**: If we prefer to keep κ tunable for future calibration
(e.g., the EB κ estimation in DESIGN_REVIEW.md Phase 2), keep the
`strand_symmetry_kappa` field in EMConfig but remove `strand_specificity`.
The CLI flag `--strand-symmetry-kappa` would remain.

→ **Decision**: Keep `strand_symmetry_kappa` in EMConfig and CLI for now.
It costs nothing and enables future calibration. Remove only
`strand_specificity`.

---

### Step 9: `src/rigel/cli.py` — CLI Changes

- **Remove**: Any `--strand-specificity` flag if it exists (it does not —
  strand_specificity is computed from the strand model, not passed via CLI).
- **Keep**: `--strand-symmetry-kappa` (per Step 8 decision).
- No other CLI changes needed.

---

### Step 10: Tests

#### 10a. Tests that must be updated

| Test file | Tests affected | Change needed |
|---|---|---|
| `test_gdna.py` | `TestLocusGDNATheta` (6 tests) | Update component count assertions, theta indexing |
| `test_gdna.py` | `TestLocusGDNAAssignment` (4 tests) | Update gDNA component indices |
| `test_gdna.py` | `TestGDNALocusAttribution` (2 tests) | Update gDNA component counting |
| `test_gdna.py` | `TestScoreGDNA` (4 tests) | Update expected gDNA LL values (no LOG_HALF) |
| `test_em_impl.py` | `TestLinkedEmHighGDNA` (3 tests) | Update component layout, expected theta |
| `test_em_impl.py` | `TestPruning.test_gdna_never_pruned` | Update to check both g_pos and g_neg |
| `test_estimator.py` | 4 gDNA tests | Update component expectations |
| `test_golden_output.py` | All ~35 scenarios | Regenerate golden outputs |
| `test_cross_chunk.py` | `test_gdna_em_count_match` | Update gDNA totals |
| `test_opp_strand_order.py` | 2 gDNA tests | Update component layout expectations |

#### 10b. Golden output regeneration

```bash
conda run -n rigel pytest tests/test_golden_output.py --update-golden -x -v
```

All ~35 golden scenarios will get new reference outputs reflecting the T+N+2
model. The numerical values will change because:
1. gDNA log-likelihoods no longer include LOG_HALF
2. Two gDNA components compete instead of one
3. M-step uses theta coupling instead of targeted-excess penalty

#### 10c. New tests to add

| Test | Purpose |
|---|---|
| `test_gdna_strand_routing` | Verify sense units → g_pos, antisense → g_neg |
| `test_gdna_theta_coupling` | Verify symmetric Beta coupling at various κ values |
| `test_gdna_no_log_half` | Verify gDNA LL = gdna_fl + gdna_log_sp (no LOG_HALF) |

---

## Code Cleanup Summary

### Removed code

| Item | Location | Why |
|---|---|---|
| `EmEquivClass.gdna_col` | em_solver.cpp:105 | Strand encoded in component identity, not side array |
| `EmEquivClass.gdna_is_anti` | em_solver.cpp:106 | Same — no longer needed |
| `gdna_idx` parameter in `build_equiv_classes` | em_solver.cpp:139 | No gDNA-specific EC handling |
| `locus_ct_arr` parameter in `build_equiv_classes` | em_solver.cpp:138 | Same |
| All gdna_col/gdna_is_anti logic in EC build | em_solver.cpp:172-194, 252-271 | ~40 lines deleted |
| Targeted-excess penalty in `map_em_step` | em_solver.cpp:600-638 | ~40 lines deleted, replaced by ~10 lines theta coupling |
| `strand_specificity` parameter | em_solver.cpp, estimator.py, pipeline.py | κ_eff scaling obsolete |
| `strand_specificity` in nanobind binding | em_solver.cpp:2416-2417 | Removed from C++ API |

### Kept code (intentionally preserved)

| Item | Location | Why |
|---|---|---|
| `LOG_HALF` constant | constants.h:76 | Still used in mRNA/nRNA unstranded fallback |
| `strand_symmetry_kappa` in EMConfig | config.py:92 | Controls coupling strength; enables future calibration |
| `--strand-symmetry-kappa` CLI flag | cli.py | Same |
| `gdna_log_liks` per-unit array | scored_fragments.py | Still needed; just no longer includes LOG_HALF |
| gDNA FL LUT in scoring.cpp | scoring.cpp:166-190 | Fragment length model unchanged |
| EB gDNA priors (`compute_eb_gdna_priors`) | locus.py:682-858 | gDNA initialization unchanged |
| `gdna_splice_penalties` | scoring.py:44-47 | Splice-type penalty unchanged |

---

## Implementation Order

The changes have dependencies. This ordering minimizes broken intermediate
states:

### Phase A: Scoring (no build needed, tests won't break yet)

1. **scoring.cpp**: Remove LOG_HALF from gDNA LL (Step 1a, 1b)

### Phase B: EM Solver (single build, core logic)

2. **em_solver.cpp — EmEquivClass**: Remove gdna_col, gdna_is_anti (Step 2a)
3. **em_solver.cpp — build_equiv_classes**: Simplify signature and body (Step 2b)
4. **em_solver.cpp — LocusSubProblem**: gdna_idx → gdna_pos_idx, gdna_neg_idx (Step 2e)
5. **em_solver.cpp — extract_locus_sub_problem**: Strand-routed gDNA injection (Step 2f)
6. **em_solver.cpp — map_em_step**: Replace penalty with theta coupling (Step 2c)
7. **em_solver.cpp — run_squarem**: Update prune protection, remove strand_specificity (Step 2d)
8. **em_solver.cpp — assign_posteriors**: Update gDNA scatter (Step 2g)
9. **em_solver.cpp — batch_locus_em**: Update build_equiv_classes call, remove strand_specificity (Step 2h)
10. **em_solver.cpp — nanobind**: Update binding (Step 2i)

### Phase C: Python (must match C++ changes)

11. **locus.py**: T+N+2 layout, strand-routed injection (Step 4)
12. **scored_fragments.py**: Comment update (Step 5)
13. **estimator.py**: Remove strand_specificity parameter (Step 6)
14. **pipeline.py**: Remove strand_specificity forwarding (Step 7)
15. **config.py**: No change (keep strand_symmetry_kappa)
16. **cli.py**: No change

### Phase D: Build + Test

17. Build: `pip install -e . --no-build-isolation`
18. Run full test suite: `pytest tests/ -x -v`
19. Fix test failures (update assertions, component counts)
20. Regenerate golden outputs: `pytest tests/test_golden_output.py --update-golden -x -v`
21. Full test suite pass: `pytest tests/ -v` — all 871+ tests green

---

## Verification and Validation

### Unit Test Verification

After implementation, every test must pass:

```bash
conda run -n rigel pytest tests/ -x -v
```

Key areas to verify:
- gDNA theta is comparable to pre-change (no catastrophic regression)
- Spliced fragments still get −∞ for both g_pos and g_neg
- Unspliced sense fragments → g_pos candidate only (not g_neg)
- Unspliced antisense fragments → g_neg candidate only (not g_pos)
- Theta coupling produces φ ≈ 0.5 when input is symmetric
- κ ≤ 2.0 disables coupling (fast path preserved)

### Ablation Validation

The primary validation uses the same 1024-run complex locus sweep
infrastructure used in Phase 1 and Phase 3 validation.

#### Step 1: Run the 1024-run sweep

```bash
mkdir -p /Users/mkiyer/Downloads/rigel_runs/complex_locus_tn2
conda run -n rigel python scripts/synthetic_sim_sweep.py \
    -c scripts/complex_locus_config.yaml \
    -o /Users/mkiyer/Downloads/rigel_runs/complex_locus_tn2 \
    -v 2>&1 | tee /Users/mkiyer/Downloads/rigel_runs/complex_locus_tn2/tn2_sweep.log
```

This uses the existing `complex_locus_config.yaml` which sweeps:
- 4 strand specificities: [0.5, 0.75, 0.9, 1.0]
- 4 gDNA fractions: [0.0, 0.1, 0.3, 0.5]
- 16 expression patterns (mRNA-only, moderate nRNA, heavy nRNA, everything)
- 4 B-locus configurations (convergent transcription on/off)
- = 1024 total runs

#### Step 2: Analyze the sweep

```bash
conda run -n rigel python scripts/analyze_complex_locus.py \
    /Users/mkiyer/Downloads/rigel_runs/complex_locus_tn2
```

#### Step 3: Compare against Phase 1 baseline

```bash
conda run -n rigel python scripts/compare_sweeps.py \
    /Users/mkiyer/Downloads/rigel_runs/complex_locus_phase1 \
    /Users/mkiyer/Downloads/rigel_runs/complex_locus_tn2 \
    --labels "Phase 1 (T+N+1)" "T+N+2"
```

#### Step 4: Three-way comparison (baseline → Phase 1 → T+N+2)

```bash
conda run -n rigel python scripts/analyze_ablation_deep.py \
    --baseline /Users/mkiyer/Downloads/rigel_runs/complex_locus \
    --phase1 /Users/mkiyer/Downloads/rigel_runs/complex_locus_phase1 \
    --phase3 /Users/mkiyer/Downloads/rigel_runs/complex_locus_tn2 \
    --labels "Baseline" "Phase 1" "T+N+2"
```

(Reuses the existing deep analysis script with relabeled directories.)

### Success Criteria

The T+N+2 model is validated if:

1. **gDNA accuracy improves at SS ≥ 0.75**: The ~40% systematic
   underestimation seen in Phase 1/Phase 3 should be substantially reduced.
   Target: |gDNA error| < 15% at SS ≥ 0.75 across all nRNA levels.

2. **mRNA accuracy does not regress**: Phase 1 achieved |mRNA error| ≈ 5.6%
   at SS ≥ 0.75. T+N+2 must stay within 1% of this.

3. **nRNA accuracy does not regress**: nRNA estimation must remain comparable
   or improve. Specifically, watch for nRNA absorption into g_pos at high
   expression loci — this would indicate the φ/nRNA identifiability problem
   (DESIGN_REVIEW.md Section 5).

4. **Unstranded (SS=0.5) cases unaffected**: At SS=0.5, both g_pos and g_neg
   receive equal fragments and the coupling drives φ → 0.5, equivalent to
   the old model. Verify no regression.

5. **All unit tests pass**: 871+ tests green, including updated golden outputs.

### Failure Modes to Watch

| Failure mode | Symptom | Mitigation |
|---|---|---|
| nRNA absorption into g_pos | nRNA underestimation at high SS + heavy nRNA | Add asymmetric constraint (DESIGN_REVIEW.md Section 5) |
| gDNA overestimation | gDNA > truth at moderate/high SS | κ too low — increase coupling strength |
| Convergence issues | SQUAREM iteration count increases | Check state vector dimensionality effect |
| EC explosion | Excessive EC count or memory | Profile EC statistics; should be similar to T+N+1 |

---

## Risk Assessment

| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| gDNA overestimation (overcorrection) | Low | Medium | κ coupling constrains g_pos/g_neg ratio |
| nRNA absorbed into g_pos | Medium | Medium | Ablation will detect; asymmetric constraint is ready |
| SQUAREM convergence slower | Low | Low | Two more components in state vector — negligible for typical loci |
| EC count increase | Low | Low | At most 2× gDNA-containing ECs; total work O(NK) unchanged |
| Test update complexity | Certain | Low | Mechanical — update indices and expected values |

---

## Appendix: Mathematical Derivation

### Per-fragment log-likelihood under T+N+2

For a fragment $f$ from unit $u$:

**mRNA candidate $t$:**
$$\ell_t(f) = \log p_{\text{strand}}(f \mid t) + \log p_{\text{FL}}(f \mid t) + \text{OH} \cdot \log p_{\text{oh}} + \log p_{\text{mm}}$$

where $\log p_{\text{strand}} = \log(\text{SS})$ if sense, $\log(1-\text{SS})$
if antisense.

**gDNA candidate $g_\text{pos}$ (positive strand):**
$$\ell_{g_\text{pos}}(f) = \begin{cases} \log p_{\text{FL}}^{\text{gdna}}(f) + \log p_{\text{splice}} & \text{if } f \text{ is sense} \\ -\infty & \text{if } f \text{ is antisense} \end{cases}$$

**gDNA candidate $g_\text{neg}$ (negative strand):**
$$\ell_{g_\text{neg}}(f) = \begin{cases} -\infty & \text{if } f \text{ is sense} \\ \log p_{\text{FL}}^{\text{gdna}}(f) + \log p_{\text{splice}} & \text{if } f \text{ is antisense} \end{cases}$$

### M-step theta coupling

After standard EM responsibility accumulation:
$$R_{g_\text{pos}} = \sum_f r_{f,g_\text{pos}}, \quad R_{g_\text{neg}} = \sum_f r_{f,g_\text{neg}}$$

Total gDNA: $R_g = R_{g_\text{pos}} + R_{g_\text{neg}}$

MAP estimate of strand balance under $\text{Beta}(\kappa/2, \kappa/2)$ prior:
$$\hat{\phi} = \frac{R_{g_\text{pos}} + \kappa/2 - 1}{R_g + \kappa - 2}$$

Apply:
$$\theta_{g_\text{pos}}^{\text{new}} = \hat{\phi} \cdot R_g, \quad \theta_{g_\text{neg}}^{\text{new}} = (1 - \hat{\phi}) \cdot R_g$$

Then normalize alongside all other components.

At $\kappa = 6$ (default), the prior is $\text{Beta}(3, 3)$, which is
moderately informative toward 0.5. This allows genuine asymmetry to emerge
from data while preventing extreme imbalance.

At $\kappa \leq 2$, the prior becomes $\text{Beta}(\leq 1, \leq 1)$ which is
improper or flat — coupling disabled.

### Equivalence to collapsed T+N+1

When $\hat{\phi} = 0.5$ exactly:
$$\theta_g = \theta_{g_\text{pos}} + \theta_{g_\text{neg}}$$

The total gDNA mass equals what T+N+1 would produce with the symmetric
strand model. T+N+2 generalizes T+N+1 by allowing $\phi \neq 0.5$.
