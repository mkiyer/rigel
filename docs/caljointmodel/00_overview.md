# 00 — Overview

## 1. The three goals of calibration

Calibration exists to serve three — and only three — downstream
goals. Everything in this design is derived from these goals.
Anything that does not contribute to one of them is out of scope.

| # | Goal | Output type | Used by |
|---|---|---|---|
| **G2** | Deconvolute each region's **contained** unspliced mass into gDNA mass and RNA mass. | Two scalars per region: $M_r^{(g, \text{cont})}, M_r^{(d, \text{cont})}$. | Locus prior pseudocounts. |
| **G3** | Deconvolute each region's **boundary** unspliced flux into gDNA mass and RNA mass, separately for the left and right boundaries. | Four scalars per region: $M_r^{(g, L)}, M_r^{(d, L)}, M_r^{(g, R)}, M_r^{(d, R)}$. | Locus prior pseudocounts (boundary-aware); hybrid-capture libraries. |
| **G4** | Estimate a per-region gDNA exposure factor with posterior uncertainty. | Two scalars per region: $\hat{\omega}_r$ (Gamma posterior mean), $\hat{\sigma}_{\log\omega_r}^2$ (delta-method variance). | Locus prior gDNA pseudocount weighting. |

> **G1 (gDNA fragment-length pmf) is intentionally not a calibration
> goal.** Earlier drafts placed a binned FL pmf $p_{\text{FL}}^{(g)}$
> here; that has been removed. The calibrator does not use FL as a
> likelihood signal — spliced fragments are FL-ambiguous, and
> per-region FL histograms add substantial substrate cost without
> improving the count/strand/spliced separation the calibrator
> actually depends on. FL continues to be applied per-fragment by
> the **downstream EM scorer** ([`scoring.py`](../../src/rigel/scoring.py),
> [`frag_length_model.py`](../../src/rigel/frag_length_model.py))
> at full resolution; that path is unaffected. Goal numbering (G2, G3,
> G4) is preserved so existing cross-references stay valid.

Spliced fragments are **not** modelled with a "splice rate" parameter
(that quantity is governed by transcript geometry and read length, not
library biochemistry). They are routed as deterministic RNA mass into
$M_r^{(d)}$; see [`01_generative_model.md`](01_generative_model.md) §3.

## 2. The substrate the calibrator consumes

The native BAM scanner already aggregates fragments into per-region and
per-boundary counters. The calibrator operates **on these aggregates
only** — there is no per-fragment iteration anywhere in the calibrator.
Per region $r$ and per boundary $b \in \{L_r, R_r\}$, the native
substrate emits:

| Sufficient statistic | Symbol | Shape |
|---|---|---|
| Unspliced fragment count | $n_r^{\text{u}}, n_b^{\text{u}}$ | scalar |
| Spliced fragment count | $n_r^{\text{s}}, n_b^{\text{s}}$ | scalar |
| Sense-strand unspliced count | $k_r^{(+)}, k_b^{(+)}$ | scalar |
| Region length | $L_r^{\text{eff}}$ | scalar |
| Region strand annotation | $\text{strand}_r$ | scalar |

No FL histogram axis. Memory budget for the calibration substrate at
$R = 200{,}000$ regions × 3 channel groups (contained, left, right) ×
12 bytes per slot ≈ 7 MB — well below previous (FL-bearing) designs.

These are the same kind of aggregates the existing
`RegionCountLedger` already produces; no FL extension is required.
See [`04_interface_contract.md`](04_interface_contract.md) §3 for the
substrate contract.

## 3. What downstream consumes

There is exactly one downstream consumer:
`compute_locus_priors_from_partitions` in `src/rigel/locus.py`. It
reads the per-region G2/G3/G4 outputs (plus the library
hyperparameter $\rho_0$) and emits the per-transcript Dirichlet
pseudocounts $(\alpha_t^{(g)}, \alpha_t^{(d)})$ for the per-locus EM.
Exact formulas are in [`04_interface_contract.md`](04_interface_contract.md) §5.

## 4. Operating principles

1. **Goal-directed.** Every line of new calibration code traces to G2–G4.
2. **Sufficient statistics, not fragments.** Calibration is $O(R)$ scalar algebra, never $O(F)$. Per-fragment scoring happens once, downstream, at the existing EM stage. Calibration does not duplicate it.
3. **Overdispersion is the default, not an option.** Count channel is Negative Binomial. Strand channel is Beta-Binomial. The Poisson and Binomial special cases would be hand-tuned approximations and we have learned to distrust those.
4. **Per-region exposure is the only continuous latent.** It has a closed-form Gamma posterior (from the NB ≡ Gamma-Poisson identity). No Newton, no Laplace, no quadrature.
5. **Library-wide hyperparameters in the M-step.** $\rho_0$, $\phi$, $\rho_d^{\text{BB}}$, $\rho_r^{\text{BB}}$, $\epsilon_s$ are jointly estimated over the library so localized pathologies cannot dominate. The gDNA strand mean $\kappa_d$ is fixed at 0.5 (biology); the per-region RNA strand mean $\kappa_r^{\text{rna}}$ is pre-computed by `StrandModel` from spliced uniques during the BAM scan.
6. **Posteriors not labels.** Outputs are masses and Gamma posteriors. Downstream weights by inverse variance. No qualitative thresholds.
7. **No constants without derivation.** Every numeric is a Bayesian prior with cited strength, a tolerance with a derived unit, or a machine-precision floor. Target: ≤ 8 total.
8. **No two-state expression latent.** Whether a region "is expressed" is a quantitative not a categorical question. The downstream consumer reads the deconvoluted RNA mass directly.

## 5. Scope: keep / delete / archive

### Keep (substrate; minor extension only)

| File | Why |
|---|---|
| `src/rigel/native/bam_scanner.cpp` | Already aggregates region/boundary counters; v1 extends with the new `CalibrationAggregates` scalars (no FL histograms). |
| `src/rigel/native/region_index.h`, `accumulator.{h,cpp}` | Native fragment routing. |
| `src/rigel/calibration/scan_payload.py` | Substrate adapter. |
| `src/rigel/calibration/_arrays.py` | Region array projections. |
| `src/rigel/calibration/region_count_ledger.py` | Region count layout. |
| `src/rigel/calibration/regions.py` | Region schema. |
| `src/rigel/calibration/fractional_evidence.py` | Signature predicates. |
| `src/rigel/strand_model.py` | Per-region RNA strand rate from spliced uniques. |
| `src/rigel/frag_length_model.py` | RNA FL pmf per strand. |

### Delete

Everything else under `src/rigel/calibration/`. Per
[`02_failure_audit.md`](02_failure_audit.md) the deletion removes ~91
hand-tuned constants and ~25 qualitative cliffs.

### Archive (for algorithmic reference)

| File | Why archive |
|---|---|
| `_fl_empirical_bayes.py` | EB Dirichlet smoothing pattern — retained for potential reuse by downstream EM. |
| `_fl_mixture.py` | Closed-form FL mixture EM — retained for potential reuse by downstream EM. |
| `fl.py` | Public FL surface — not consumed by the new calibrator. |

## 6. Boundaries are first-class

The boundary deconvolution (G3) uses **the same per-region inference
machinery** on a parallel set of sufficient statistics. Boundary flux
carries the only gDNA signal in hybrid-capture libraries (intron and
intergenic depleted; gDNA detected by fragments that splash out of
captured exons into flanking sequence). The boundary substrate
counters are populated identically by the native accumulator
(the boundary kind is already tagged); the calibrator treats each
boundary as a small "region" with its own sufficient statistics. After
inference, boundary gDNA mass is half-split into the two neighbouring
regions for downstream pseudocount construction (avoids double-count).

There is no "v2" deferral. G3 ships in v1.

## 7. The inference shape in one paragraph

For each region (and each boundary, by the same code), the unspliced
count $n_r^{\text{u}}$ is modelled as $\mathrm{NB}(\mu_r^{(g)} + \mu_r^{(d)}, \phi)$
with $\mu_r^{(g)} = \omega_r \rho_0 L_r^{\text{eff}}$ and $\mu_r^{(d)}$
the inferred RNA mean. The sense count $k_r^{(+)}$ is a
Beta-Binomial mixture under the per-region mixing proportion
$\pi_r^{(g)}$ (gDNA strand mean $\kappa_d = 0.5$ fixed by biology; the
per-region RNA strand mean $\kappa_r^{\text{rna}}$ is pre-computed by
`StrandModel`). Per-region E-step: combine the count log-evidence
ratio and the strand log-evidence ratio via sigmoid to soft-allocate
$n_r^{\text{u}}$ into $(M_r^{(g)}, M_r^{(d)})$. Spliced fragments add
directly to $M_r^{(d)}$ (deterministic RNA) and also contribute a
Poisson failsafe term against spurious gDNA mass via $\epsilon_s$.
Per-region exposure update: closed-form Gamma posterior
$\omega_r \mid M_r^{(g)} \sim \mathrm{Gamma}(1/\phi + M_r^{(g)},\, 1/\phi + \rho_0 L_r^{\text{eff}})$.
Global M-step: closed forms for $\rho_0$ and $\epsilon_s$; 1-D Newton
for the two dispersion scalars $\phi$ and $\rho_d^{\text{BB}}$, and
(optionally) $\rho_r^{\text{BB}}$. Outer EM converges in 3–10
iterations. **FL is not used as a likelihood signal in calibration**
— it is consumed per-fragment by the downstream EM scorer.

## 8. Definition of done (design phase)

Design is signed off when:
1. The four goals (§1) are agreed.
2. Sufficient-statistic substrate (§2) is agreed.
3. Generative model in [`01_generative_model.md`](01_generative_model.md), including NB count and Beta-Binomial strand, is agreed.
4. Splice-as-deterministic-RNA treatment is agreed (no `q̄` hyperparameter).
5. Inference algorithm in [`03_inference.md`](03_inference.md) is agreed.
6. Interface contract in [`04_interface_contract.md`](04_interface_contract.md), including exact downstream pseudocount formulas, is agreed.
7. Implementation phases in [`05_implementation_plan.md`](05_implementation_plan.md) are agreed.
8. Validation plan in [`06_validation_plan.md`](06_validation_plan.md) is agreed.

Only then does Phase 1 (Archive) begin.

## 9. Document index

| Doc | Purpose |
|---|---|
| [00_overview.md](00_overview.md) | Goals, substrate, scope, principles, doc index. |
| [01_generative_model.md](01_generative_model.md) | Region/boundary sufficient-statistic model; NB count, Beta-Binomial strand, deterministic spliced RNA, Poisson failsafe. |
| [02_failure_audit.md](02_failure_audit.md) | Why the legacy must be burned. |
| [03_inference.md](03_inference.md) | Per-region E-step (count + strand sigmoid mix, closed-form Gamma exposure); global M-step (closed forms + 1-D Newton scalars). |
| [04_interface_contract.md](04_interface_contract.md) | Substrate aggregates; `CalibrationConfig`, `CalibrationResult`; exact downstream pseudocount formulas. |
| [05_implementation_plan.md](05_implementation_plan.md) | Seven-phase rebuild with sub-phases per goal. |
| [06_validation_plan.md](06_validation_plan.md) | Unit, synthetic, scenario (paralog + hybrid-capture), armis2. |
