# Phase 3: Scale Strand Symmetry Penalty by Strand Specificity

## Summary

Scale the gDNA strand symmetry penalty in the EM M-step by the library's
strand specificity (SS). Currently the penalty operates at fixed strength
regardless of how much strand information the library actually provides.
At SS=0.5 (no strand info), the penalty wrongly penalises *all* asymmetric
gDNA — including gDNA that only appears asymmetric because RNA reads are
strand-biased. At SS≥0.75, the penalty is too weak to prevent nRNA from
absorbing gDNA's intronic reads.

The fix: replace the fixed `strand_symmetry_kappa` with an effective kappa
that scales by SS:

$$\kappa_{eff} = \max\bigl(\kappa \cdot (2 \cdot SS - 1)^2,\; 2.0\bigr)$$

Fold `strand_symmetry_pseudo` into the κ scaling, removing it as a separate
parameter.

**Root causes addressed**: RC1 (synergy with Phase 1), RC6 (gDNA mass
misallocation)

**Parameters removed**: 1 (`strand_symmetry_pseudo`) → 17 remaining

---

## Rationale

### Phase 1 exposed the gDNA–nRNA separation problem

Phase 1 (decouple nRNA) freed nRNA from the Beta prior constraint. This
dramatically improved mRNA and nRNA accuracy at SS≥0.75, but introduced a
new failure mode: **gDNA is now systematically underestimated** at SS≥0.75
when nRNA is present.

| SS   | gDNA error (baseline) | gDNA error (Phase 1) | nRNA error (baseline) | nRNA error (Phase 1) |
|------|-----------------------|----------------------|-----------------------|----------------------|
| 0.75 | +33.9%                | **-47.8%**           | -18.3%                | **+22.5%**           |
| 0.9  | +31.4%                | **-45.8%**           | -18.8%                | **+17.5%**           |
| 1.0  | +42.9%                | **-26.5%**           | -22.6%                | **+8.8%**            |

The mass flow simply inverted: baseline pushed mass from nRNA → gDNA
(overestimate gDNA); Phase 1 pushes mass from gDNA → nRNA (underestimate
gDNA).

### Why the current penalty doesn't help

The current strand symmetry penalty in `map_em_step` computes:

```
p̂ = (e_sense_gdna + ε/2) / (N_g + ε)        # gDNA strand fraction
balance = 4 · p̂ · (1 - p̂)                    # ∈ [0, 1], =1 when symmetric
w_sym = balance^(κ/2 - 1)                     # discount on asymmetric excess
```

Where:
- `κ` = `strand_symmetry_kappa` (default 6.0) — penalty strength
- `ε` = `strand_symmetry_pseudo` (default 50.0) — pseudo-count regulariser

**The problem**: κ=6.0 is applied uniformly regardless of SS. But the
optimal penalty strength depends critically on SS:

- **At SS=1.0**: Perfect strand info. Symmetric gDNA counts genuinely
  indicate gDNA. The penalty should be strong (high κ) to protect symmetric
  gDNA from nRNA absorption.
- **At SS=0.75**: Good strand info but imperfect. The penalty should be
  moderately strong.
- **At SS=0.5**: No strand info. All reads appear ~50/50 sense/anti regardless
  of source. The penalty should be disabled (κ=2) because strand balance
  carries no gDNA-specific information.

### The `strand_symmetry_pseudo` parameter is redundant

The pseudo-count `ε` serves as an evidence threshold: at low gDNA counts,
`p̂` is pulled toward 0.5 by the pseudo-count, making `balance` ≈ 1 and
`w_sym` ≈ 1. This is an indirect way to reduce penalty strength at low counts.

But this is better handled by the SS scaling itself:
- At low SS, κ_eff is naturally small → penalty is naturally weak
- At high SS, κ_eff is strong → penalty engages even at low counts
  (which is correct — at SS=1.0, even a few symmetric reads are strong
  evidence for gDNA)

The current `ε=50` was tuned to compensate for the fixed-strength penalty.
With SS-aware scaling, this compensation is no longer needed. We fold the
evidence-gating role into the κ_eff formula and remove the parameter.

### Mathematical formulation

**Current** (fixed κ, separate ε):
```
w_sym = [4 · p̂(ε) · (1 - p̂(ε))]^(κ/2 - 1)
```

**Proposed** (SS-scaled κ, no ε):
```
κ_eff = max(κ · (2·SS - 1)², 2.0)
p̂ = (e_sense + 0.5) / (N_g + 1.0)           # minimal Laplace smoothing
w_sym = [4 · p̂ · (1 - p̂)]^(κ_eff/2 - 1)
```

Behavior table:

| SS   | (2·SS-1)² | κ_eff (κ=6) | Penalty effect        |
|------|-----------|-------------|------------------------|
| 0.50 | 0.00      | 2.0         | Disabled (w_sym = 1)   |
| 0.60 | 0.04      | 2.0         | Disabled (clamped)     |
| 0.75 | 0.25      | 2.0         | Disabled (clamped)     |
| 0.80 | 0.36      | 2.16        | Very mild              |
| 0.90 | 0.64      | 3.84        | Moderate               |
| 0.95 | 0.81      | 4.86        | Strong                 |
| 1.00 | 1.00      | 6.0         | Full strength          |

Note: with the κ=2 clamp, the penalty only engages above SS≈0.79
(where `κ · (2·SS−1)² > 2`). This is desirable — below ~0.8, strand
evidence is too weak to reliably distinguish gDNA symmetry from RNA
strand noise.

---

## Changes by File

### 1. `src/rigel/config.py` — Remove `strand_symmetry_pseudo`

**Delete** the `strand_symmetry_pseudo` field from `EMConfig`:

```python
# REMOVE this field entirely:
strand_symmetry_pseudo: float = 50.0
"""Bayesian pseudo-count (α₀) for gDNA strand fraction regularisation.
Acts as a Beta(α₀, α₀) prior belief: "α₀ sense and α₀ antisense
fragments" before observing data.  Controls evidence threshold
for the symmetry penalty to engage.  Higher → more forgiving
at low fragment counts.
"""
```

**Update** the `strand_symmetry_kappa` docstring to mention SS scaling:

```python
strand_symmetry_kappa: float = 6.0
"""Strand symmetry penalty strength for gDNA in the M-step.
Effective κ is scaled by strand specificity: κ_eff = κ · (2·SS − 1)².
Set ≤ 2.0 to disable the penalty entirely.
"""
```

---

### 2. `src/rigel/native/em_solver.cpp` — Scale κ by SS, remove `strand_eps`

#### 2a. `map_em_step` — Replace `strand_eps` with `strand_specificity`

Change the signature (line ~543):

```cpp
// BEFORE:
static void map_em_step(
    ...
    double        strand_symmetry_kappa = 2.0,
    double        strand_eps = 10.0)

// AFTER:
static void map_em_step(
    ...
    double        strand_symmetry_kappa = 2.0,
    double        strand_specificity = 0.5)
```

Modify the penalty block (line ~597):

```cpp
// BEFORE:
if (strand_symmetry_kappa > 2.0 && n_components > 0) {
    ...
    if (N_g > 0.0) {
        double pseudo_half = strand_eps / 2.0;
        double p_hat = (e_sense_gdna + pseudo_half) / (N_g + strand_eps);
        double balance = 4.0 * p_hat * (1.0 - p_hat);
        w_sym = std::pow(balance, strand_symmetry_kappa / 2.0 - 1.0);
    }
    ...
}

// AFTER:
// Compute effective kappa scaled by strand specificity
double ss_scale = 2.0 * strand_specificity - 1.0;
double kappa_eff = strand_symmetry_kappa * ss_scale * ss_scale;
if (kappa_eff > 2.0 && n_components > 0) {
    ...
    if (N_g > 0.0) {
        // Laplace smoothing (replaces strand_eps pseudo-count)
        double p_hat = (e_sense_gdna + 0.5) / (N_g + 1.0);
        double balance = 4.0 * p_hat * (1.0 - p_hat);
        w_sym = std::pow(balance, kappa_eff / 2.0 - 1.0);
    }
    ...
}
```

The key changes:
1. `kappa_eff = κ · (2·SS − 1)²` replaces fixed `κ`
2. Guard condition changes from `κ > 2.0` to `κ_eff > 2.0` (penalty auto-disables at low SS)
3. `strand_eps` pseudo-count replaced by minimal Laplace smoothing `+0.5 / +1.0`

#### 2b. `run_squarem` — Replace strand_eps with strand_specificity

Same signature change. Both MAP-EM and VBEM paths call into this.

```cpp
// BEFORE:
static EMResult run_squarem(
    ...
    double        strand_symmetry_kappa = 2.0,
    double        strand_eps = 10.0)

// AFTER:
static EMResult run_squarem(
    ...
    double        strand_symmetry_kappa = 2.0,
    double        strand_specificity = 0.5)
```

Forward to `map_em_step` calls (lines ~889, ~895, ~932):

```cpp
map_em_step(state0.data(), ec_data, log_eff_len,
            unambig_totals, prior, em_totals.data(),
            state1.data(), n_components,
            estep_threads, pool,
            strand_symmetry_kappa, strand_specificity);
```

**Note on VBEM path**: Currently `vbem_step` does NOT apply the strand
penalty (pre-existing inconsistency). Phase 3 does NOT change this — VBEM
strand penalty support is out of scope. The MAP-EM path (default) is the
only one affected.

#### 2c. `batch_locus_em` — Replace `strand_eps` with `strand_specificity`

Same signature change (line ~1738):

```cpp
// BEFORE:
    double strand_symmetry_kappa = 2.0,
    double strand_eps = 10.0)

// AFTER:
    double strand_symmetry_kappa = 2.0,
    double strand_specificity = 0.5)
```

Forward to `run_squarem` (line ~1970):

```cpp
EMResult result = run_squarem(
    ec_data, log_eff_len.data(),
    sub.unambig_totals.data(),
    prior.data(),
    theta_init.data(),
    nc, max_iterations, convergence_delta,
    use_vbem, prune_threshold,
    estep_thr, pool,
    strand_symmetry_kappa, strand_specificity);
```

#### 2d. Nanobind binding — Replace `strand_eps` with `strand_specificity`

At line ~2414:

```cpp
// BEFORE:
    nb::arg("strand_symmetry_kappa") = 2.0,
    nb::arg("strand_eps") = 10.0,

// AFTER:
    nb::arg("strand_symmetry_kappa") = 2.0,
    nb::arg("strand_specificity") = 0.5,
```

---

### 3. `src/rigel/estimator.py` — Pass `strand_specificity` instead of `strand_symmetry_pseudo`

The `run_batch_locus_em` method (line ~392) needs to pass `strand_specificity`
instead of `strand_symmetry_pseudo`:

```python
# BEFORE:
    self.em_config.strand_symmetry_kappa,
    self.em_config.strand_symmetry_pseudo,

# AFTER:
    self.em_config.strand_symmetry_kappa,
    self._strand_specificity,
```

This requires adding a `strand_specificity` attribute to the estimator.
There are two options:

**Option A** (simpler): Add it as a parameter to `run_batch_locus_em`:

```python
def run_batch_locus_em(
    self,
    loci: list[Locus],
    em_data: ScoredFragments,
    index: TranscriptIndex,
    gdna_inits: np.ndarray,
    *,
    strand_specificity: float = 0.5,   # NEW
    em_iterations: int = 1000,
    ...
```

**Option B** (matching existing pattern): Set it as an attribute before
calling `run_batch_locus_em`, same as `nrna_init`:

```python
estimator.strand_specificity = strand_models.strand_specificity
```

**Recommendation**: Option A — explicit parameter, avoids hidden state.
The caller in `pipeline.py` already has `strand_specificity` available.

---

### 4. `src/rigel/pipeline.py` — Pass strand_specificity to the EM

In `_run_em_phase` (around line ~500, where `run_batch_locus_em` is called):

```python
# BEFORE:
total_gdna_em, locus_mrna, locus_nrna, locus_gdna = estimator.run_batch_locus_em(
    loci, em_data, index, gdna_inits,
    em_iterations=..., em_convergence_delta=...,
    confidence_threshold=...,
)

# AFTER:
total_gdna_em, locus_mrna, locus_nrna, locus_gdna = estimator.run_batch_locus_em(
    loci, em_data, index, gdna_inits,
    strand_specificity=strand_models.strand_specificity,
    em_iterations=..., em_convergence_delta=...,
    confidence_threshold=...,
)
```

---

### 5. `src/rigel/cli.py` — Remove `--strand-symmetry-pseudo` argument

Delete the CLI argument definition (lines ~833-837):

```python
# REMOVE:
adv.add_argument(
    "--strand-symmetry-pseudo", dest="strand_symmetry_pseudo",
    type=float, default=None,
    help="Bayesian pseudo-count (α₀) for gDNA strand fraction. "
         "Controls evidence threshold for the symmetry penalty. "
         "Higher = more forgiving at low counts (default: 50.0).",
)
```

Remove from `_PARAM_SPECS` (line ~424):

```python
# REMOVE:
_ParamSpec("strand_symmetry_pseudo", "em.strand_symmetry_pseudo"),
```

Update `--strand-symmetry-kappa` help text:

```python
adv.add_argument(
    "--strand-symmetry-kappa", dest="strand_symmetry_kappa",
    type=float, default=None,
    help="Strand symmetry penalty κ for gDNA in the M-step. "
         "Effective κ scales by strand specificity: κ_eff = κ·(2·SS−1)². "
         "Set ≤ 2.0 to disable (default: 6.0).",
)
```

---

### 6. Scripts — Remove `strand_symmetry_pseudo` references

#### `scripts/diagnostic_init_params.py` (~line 106, 134)

Remove references to `strand_symmetry_pseudo` in the diagnostic parameter
output and EM call.

#### `scripts/analyze_sweep.py` (~line 54)

Remove `strand_symmetry_pseudo` from the printed sweep parameters.

#### `scripts/benchmarking/gdna_kappa_pseudo_grid.yaml` (~line 70)

Remove `strand_symmetry_pseudo` from sweep parameters. This sweep config
was specifically designed to grid-search κ × pseudo combinations — it should
now only sweep κ.

---

## Test Changes

### `tests/test_gdna.py`

1. **Remove** assertion `assert cfg.strand_symmetry_pseudo == 50.0`
   (line ~926).

2. **Update** any test that directly passes `strand_symmetry_pseudo` or
   `strand_eps` to the C++ native function.

### `tests/test_em_impl.py`

No changes needed — this file does not reference strand symmetry parameters.

### `tests/test_bias.py`

Check for any direct `run_locus_em_native` calls that pass the old
`strand_eps` parameter. Update to `strand_specificity` if found.

### `tests/test_golden_output.py`

The golden reference files will need regeneration since the EM numerics
change when the penalty strength varies by SS. Run:

```bash
pytest tests/test_golden_output.py -v --update-golden
```

### New tests

Add a focused test in `tests/test_gdna.py`:

```python
class TestStrandPenaltySSScaling:
    """Verify strand symmetry penalty scales with strand specificity."""

    def test_penalty_disabled_at_ss_half(self):
        """At SS=0.5, κ_eff = κ·0² = 0 ≤ 2 → penalty disabled."""
        # Run EM with gDNA that has asymmetric strand counts
        # Verify gDNA mass is NOT penalised at SS=0.5

    def test_penalty_full_at_ss_one(self):
        """At SS=1.0, κ_eff = κ·1² = κ → full penalty."""
        # Run EM with gDNA that has asymmetric strand counts
        # Verify gDNA mass IS penalised at SS=1.0

    def test_penalty_scales_monotonically(self):
        """Higher SS → stronger penalty → less gDNA mass."""
        # Run same locus at SS=0.5, 0.75, 0.9, 1.0
        # Verify gDNA decreases monotonically as SS increases
        # (for a case with asymmetric gDNA strand counts)
```

---

## Validation

1. **Re-run the 1,024-run complex locus ablation**:
   ```bash
   conda run -n rigel python scripts/synthetic_sim_sweep.py \
     -c scripts/complex_locus_config.yaml \
     -o /Users/mkiyer/Downloads/rigel_runs/complex_locus_phase3 \
     -v 2>&1 | tee phase3_sweep.log
   ```

2. **Run analysis**:
   ```bash
   conda run -n rigel python scripts/analyze_complex_locus.py \
     /Users/mkiyer/Downloads/rigel_runs/complex_locus_phase3/sweep_results.tsv
   ```

3. **Compare key metrics against Phase 1**:

   | Metric | Phase 1 | Phase 3 Target |
   |--------|---------|----------------|
   | gDNA error at SS≥0.75 with nRNA (currently -48%) | -47.8% | < ±15% |
   | nRNA error at SS≥0.75 (currently +22.5%) | +22.5% | < ±10% |
   | mRNA accuracy at SS≥0.75 | 1–34% | Same or better |
   | SS=0.5 gDNA error | ~+270% | Same (penalty disabled) |
   | SS=0.5 nRNA/mRNA | -81% / 34–58% | Same (penalty disabled) |
   | Pure mRNA accuracy (no nRNA, no gDNA) | 0.0% | Unchanged |
   | Pattern 15 mRNA FP (SS≥0.9) | 220–603 | Same or better |

4. **Run full test suite**: `pytest tests/ -q`

---

## Implementation Order

1. **C++ changes**: Modify `map_em_step`, `run_squarem`, `batch_locus_em`,
   and the nanobind binding to accept `strand_specificity` instead of
   `strand_eps`. Implement the κ_eff formula.

2. **Python plumbing**: Add `strand_specificity` parameter to
   `run_batch_locus_em` in `estimator.py`. Pass from `pipeline.py`.

3. **Config cleanup**: Remove `strand_symmetry_pseudo` from `EMConfig`,
   CLI, and `_PARAM_SPECS`.

4. **Script cleanup**: Remove `strand_symmetry_pseudo` from diagnostic and
   sweep scripts.

5. **Tests**: Update `test_gdna.py`, regenerate golden outputs, add new
   SS scaling tests.

6. **Full test run**: `pytest tests/ -q` — target 871 passed, 0 failed.

---

## Risk Assessment

**Low risk**. This change is localized to a single formula in the M-step
and parameter plumbing. The code paths are well-understood from the
exploration above.

| Risk | Mitigation |
|------|-----------|
| κ_eff too weak at moderate SS, gDNA still underestimated | The formula is monotonic in SS; if needed, use a steeper exponent e.g. `(2·SS−1)^1` instead of `(2·SS−1)^2` |
| Laplace smoothing (+0.5/+1.0) too aggressive/weak | Minimal change from large `ε=50` → `ε=1`; the SS-scaling matters more than the smoothing |
| VBEM path inconsistency (no strand penalty) | Pre-existing issue, out of scope. Default mode is MAP-EM. |
| Golden output changes | Expected — regenerate via `--update-golden` |
| Sweep configs referencing `strand_symmetry_pseudo` | grep and update all YAML/Python scripts |
