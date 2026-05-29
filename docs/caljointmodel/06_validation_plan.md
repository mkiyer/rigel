# 06 — Validation Plan

Four validation layers, ordered cheapest to most expensive. A
regression at any layer blocks promotion to the next.

| Layer | Scope | Runtime | Phase |
|---|---|---|---|
| 1. Unit tests | Per-function correctness | < 30 s | 3, 4, 5 |
| 2. Synthetic end-to-end | Hyperparameter recovery | < 1 min | 4 |
| 3. Scenario tests | Paralog, hybrid-capture, … | < 5 min | 5, 6 |
| 4. Real-data benchmarks | armis2 sweeps | hours | 6 |

## 1. Unit tests (`tests/calibration/`)

### 1.1 G2 / G3 — Soft allocation

> G1 (gDNA FL pmf) was dropped from the calibration goals — FL is no
> longer a calibration channel. See
> [`00_overview.md`](00_overview.md) §1.

| Test | Asserts |
|---|---|
| `test_estep_recovers_known_pi_g` | Synthetic substrate where $\pi_r^{(g)}$ is known → recovered within 5% on regions with ≥ 20 fragments. |
| `test_estep_spliced_routes_to_rna` | Region with $n_r^{\text{s}} = 100$ → $M_r^{(d, \text{cont})}$ at least 99 (modulo $\epsilon_s$). |
| `test_estep_unstranded_ignores_strand_channel` | $\kappa_r^{\text{rna}} = 0.5$ → strand LLR within $10^{-3}$ of zero. |
| `test_estep_count_strand_dominate_in_paralog` | 126/26 strand-asymmetric gDNA straddlers (no FL signal) → $M_r^{(d, \text{cont})} < 10$ via the three-leg rescue (count NB matching $\mu_r^{(g)}$ + BB strand absorbing 126/26 + zero spliced). **Risk-flagged**: empirically weaker than the prior FL-bearing design; benchmark layer 4 must verify. |
| `test_estep_vectorized_matches_scalar` | 1000 random regions. |
| `test_mass_conservation` | $M_r^{(g)} + M_r^{(d)}$ = total fragment count to $10^{-9}$. |
| `test_boundary_half_split_conservation` | $M_r^{(g, L)} + M_r^{(g, R)}$ correctly half-split into two adjacent regions. |
| `test_betabinom_at_50_50_split` | Beta-Binomial likelihood under fixed $\kappa_d = 0.5$ for various $\rho_d^{\text{BB}}$ produces non-spurious LLR (no phantom-RNA from balanced gDNA strand splits). |
| `test_kappa_d_is_fixed_at_half` | Calibration result has no `kappa_d` field; internally always 0.5 by biological constraint. |
| `test_kappa_rna_passthrough_from_strand_model` | Per-region $\kappa_r^{\text{rna}}$ used internally equals `StrandModel.p_r1_sense` exactly; never re-fitted. |

### 1.2 G4 — Exposure

| Test | Asserts |
|---|---|
| `test_exposure_closed_form` | Gamma posterior mean matches formula. |
| `test_exposure_empty_region` | $M = 0$ → $\omega = 1$, $\log\sigma^2 = \phi$. |
| `test_exposure_variance_scaling` | $\hat{\sigma}^2_{\log\omega} = 1 / (1/\phi + M)$. |

### 1.3 M-step

| Test | Asserts |
|---|---|
| `test_rho_0_recovers_truth` | Synthetic with known $\rho_0$ → within 5% relative. |
| `test_phi_recovers_truth` | Synthetic NB with known $\phi$ → within 20% relative. |
| `test_rho_d_bb_recovers_truth` | Synthetic gDNA BB with known $\rho_d^{\text{BB}}$ (mean fixed at 0.5) → within 30% relative. |
| `test_rho_r_bb_recovers_truth` | Synthetic RNA BB with known $\rho_r^{\text{BB}}$ and given $\kappa_r^{\text{rna}}$ → within 30% relative. |
| `test_eps_s_closed_form` | Matches Beta(1,1)+counts formula exactly. |

### 1.4 Outer loop

| Test | Asserts |
|---|---|
| `test_mass_change_monotone` | Mass-change diagnostic non-increasing across iterations. |
| `test_converges_within_max_outer` | Default 25 iters is far above what's needed. |
| `test_convergence_error_on_increase` | Manually corrupting an M-step output to force ELBO regress raises `CalibrationConvergenceError`. |

### 1.5 Substrate + result invariants

| Test | Asserts |
|---|---|
| `test_substrate_invariant_violations_raise` | Each invariant in [`04_interface_contract.md`](04_interface_contract.md) §4.1 raises. |
| `test_result_invariants_enforced` | Each invariant in §5.1 raises on violation. |

### 1.6 Locus prior consumer (Phase 5)

| Test | Asserts |
|---|---|
| `test_pseudocount_formula_exact` | Matches §6.2 formulas to $10^{-10}$. |
| `test_mass_conservation_invariant` | §6.3 invariant holds. |
| `test_symmetric_paralog_locus_symmetric_pseudocount` | Two transcripts sharing region equally → identical pseudocounts. |
| `test_empty_region_only_gdna_prior_contribution` | Only $\kappa w_r \rho_0 L_{t,r}^{\text{eff}}$ contributes. |

### 1.7 Acceptance

All unit tests pass on every commit.

## 2. Synthetic end-to-end

`tests/calibration/test_e2e_synthetic.py`:

```python
def test_calibrate_recovers_known_hyperparams():
    rng = np.random.default_rng(0)
    R = 1000
    truth = {
        "rho_0": 0.005,
        "phi": 0.2,
        "rho_d_bb": 0.02,
        "rho_r_bb": 0.01,
        "eps_s": 1e-3,
    }
    substrate, omega_true, pi_g_true = simulate_calibration_substrate(R, truth, rng)

    result = calibrate(
        substrate, strand_model_fixture(),
        config=CalibrationConfig(),
    )

    assert result.converged
    assert abs(result.rho_0 - truth["rho_0"]) / truth["rho_0"] < 0.10
    assert abs(result.phi - truth["phi"]) / truth["phi"] < 0.25
    assert abs(result.rho_d_bb - truth["rho_d_bb"]) / truth["rho_d_bb"] < 0.50
    assert abs(result.rho_r_bb - truth["rho_r_bb"]) / truth["rho_r_bb"] < 0.50

    omega_err = np.abs(np.log(result.omega) - np.log(omega_true))
    well_observed = (omega_true * truth["rho_0"] * substrate.L_eff) >= 10
    assert np.median(omega_err[well_observed]) < 0.10
```

### Acceptance thresholds

| Parameter | Tolerance | Rationale |
|---|---|---|
| $\rho_0$ | 10% rel | Library-wide mean; small SE at $R = 1000$. |
| $\phi$ | 25% rel | Overdispersion harder to pin down. |
| $\rho_d^{\text{BB}}$ | 50% rel | gDNA strand dispersion (mean fixed at 0.5). |
| $\rho_r^{\text{BB}}$ | 50% rel | RNA strand dispersion (mean $\kappa_r^{\text{rna}}$ pre-computed). |
| $\rho_d^{\text{BB}}$ | 50% rel | Dispersion estimator has high SE at moderate counts. |
| $\omega_r$ median log-error | 0.10 | On regions with ≥ 10 expected gDNA fragments. |

## 3. Scenario tests

### 3.1 Paralog scenario (root motivator)

`tests/scenarios_aligned/test_paralogs.py`. Post-rewrite expected:

| Quantity | Legacy (broken) | New (expected) |
|---|---|---|
| Paralog region `mass_d_contained` | ~26 | < 5 |
| Paralog region `mass_g_contained` | ~125 | ~145 |
| Locus EM t1 / t2 split | 4 / 638 | between 400/600 and 600/400 |

**Acceptance.** Passes naturally. t1, t2 ∈ [400, 600].

### 3.2 Hybrid-capture scenario (new)

`tests/scenarios_aligned/test_hybrid_capture.py` — to be authored in Phase 6.

Setup: single annotated exon with neighbouring intronic flanks.
Library = 95% RNA from exon + 5% gDNA. gDNA captured at the same 90%
rate as the exon's RNA inside the exonic boundary, < 5% in intronic
flanks. Therefore:
- Contained gDNA mass in exon: moderate.
- Contained gDNA mass in intronic flank: ~zero.
- Boundary flux at exon edges: moderate, gDNA-dominated.

**Acceptance.** Recovered exon RNA mass within 15% of truth.
Boundary gDNA mass within 20% of truth.

### 3.3 Other scenario BAMs

All BAMs in `tests/scenarios_aligned/` revalidated; visual inspection
of outputs (goldens will be regenerated).

## 4. Synthetic benchmark sweeps

### 4.1 Sweeps

| Config | Tests |
|---|---|
| `locus_simple_baseline.yaml` | Nominal across ss × gdna × n_rna grid |
| `locus_simple_em_prior.yaml` | EM stability vs prior concentration |
| `locus_simple_strand.yaml` | Strand specificity sweep |
| `locus_simple_em_mode.yaml` | MAP vs VBEM (orthogonal) |
| `locus_simple_scoring.yaml` | Scoring penalty sweep (orthogonal) |

### 4.2 Acceptance

- **mRNA recovery**: median relative error < 0.15.
- **nRNA recovery**: nRNA siphon reduced from baseline ~37% to < 20% worst case.
- **gDNA recovery**: low-contamination relative error < 1.3× (from baseline 1.5–2.0×).

### 4.3 Comparison procedure

```bash
git checkout main
for cfg in scripts/benchmark/configs/locus_simple_*.yaml; do
  name=$(basename "$cfg" .yaml)
  python scripts/sim/locus_sweep.py -c "$cfg" -o "scratch/preburn/$name"
done

git checkout calibration-rebuild
for cfg in scripts/benchmark/configs/locus_simple_*.yaml; do
  name=$(basename "$cfg" .yaml)
  python scripts/sim/locus_sweep.py -c "$cfg" -o "scratch/postburn/$name"
done

python scripts/benchmark/analyze_golden.py --baseline scratch/preburn --candidate scratch/postburn
python scripts/benchmark/analyze_deep.py --baseline scratch/preburn --candidate scratch/postburn
```

## 5. Armis2 real-data benchmarks

### 5.1 Corner conditions

- `gdna_none_ss_1.00_nrna_none` — easiest.
- `gdna_none_ss_0.50_nrna_rand` — unstranded with nRNA.
- `gdna_high_ss_1.00_nrna_none` — high gDNA.
- `gdna_high_ss_0.50_nrna_rand` — hardest.

### 5.2 Procedure

```bash
conda activate rigel
# scripts/benchmarking/configs/calibration_rebuild.yaml defines a new named
# rigel config "calibration_rebuild_v1"

python -m scripts.benchmarking run -c scripts/benchmarking/configs/calibration_rebuild.yaml \
  --conditions gdna_none_ss_1.00_nrna_none gdna_none_ss_0.50_nrna_rand \
               gdna_high_ss_1.00_nrna_none gdna_high_ss_0.50_nrna_rand

python -m scripts.benchmarking analyze -c scripts/benchmarking/configs/calibration_rebuild.yaml \
  -o results/calibration_rebuild_v1
```

### 5.3 Acceptance

vs latest legacy run:
- mRNA-level Spearman vs truth: equal or better.
- gDNA absolute fraction: within 5% absolute for `gdna_high`, within 2% absolute for `gdna_none`.
- nRNA recovery: improved on `nrna_rand`.
- Runtime: within 2× of legacy.

## 6. Numerical-stability stress

`tests/calibration/test_numerical_stability.py`:

| Input | Expected |
|---|---|
| All $n_r^{\text{u}} = 0$ | $M = 0$, $\omega = 1$, $\log\sigma^2 = \phi$. Diagnostics finite. |
| Single region | Calibration runs. |
| All-intergenic signature | All hyperparameters finite; $\pi_r^{(g)} \to 1$. |
| Extreme imbalance: 99.99% intergenic + 0.01% exonic | No overflow. |
| All-sense or all-antisense strand data | $\rho_r^{\text{BB}}$ fitted; gDNA likelihood still centred at 0.5 (no drift). |
| Tiny $\rho_0 \to 0$ limit | $\rho_0$ floored; $\pi_r^{(g)} \to 0$ smoothly. |

## 7. Reporting

After Phase 6: `docs/caljointmodel/validation_report_2026_05.md`:

1. Unit test pass count + runtime.
2. Synthetic-recovery table.
3. Paralog scenario before/after.
4. Hybrid-capture scenario before/after.
5. Sweep-by-sweep comparison plots.
6. Armis2 condition-by-condition comparison.
7. Regressions + root causes.
8. Sign-off: "The rebuild meets all acceptance criteria for promotion to Phase 7."

## 8. Explicit non-goals for v1

- **FL as a calibration channel.** FL is consumed only by the downstream EM scorer (`scoring.py`, `frag_length_model.py`) at full per-fragment resolution. The calibrator does not accumulate per-region or per-boundary FL histograms; see [`00_overview.md`](00_overview.md) §1.
- **Mappability-aware $L_r^{\text{eff}}$.** v1 uses physical length.
- **Per-region BB dispersion.** Single global $\rho_d^{\text{BB}}$.
- **Junction-geometric RNA splice rate model.** Not needed; spliced fragments are deterministic RNA.

Tracked in `docs/caljointmodel/TODO.md` (created Phase 7).
