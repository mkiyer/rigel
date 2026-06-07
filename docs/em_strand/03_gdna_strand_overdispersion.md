# Design — restore the gDNA strand Beta-Binomial overdispersion (ρ_d)

**Status (2026-06-07): Phase 1 + Phase 2 SHIPPED.** Fitted (`gdna_strand.py`, pooled MoM +
identifiability gate), applied (`strand_likelihood.py` two-component variance), wired
(`calibrate.py` → `deconv_regions/sides`), surfaced (`CalibrationResult` + `summary.json`),
simulated (`GDNAConfig.gdna_strand_overdispersion`), tested (estimator + wrapper + end-to-end
calibrate recovery across {0.05,0.1,0.2} + 126/26 regression + od→0 byte-identical). Full suite
989 green, ruff clean, goldens byte-identical (the gate floors od on the tiny test scenarios).
**Two new constants** (`_MIN_SEED_NODES=30`, `_SIGNIFICANCE_Z=2`) flagged for review per
no-magic-numbers (§1.2). Phase 3 remaining: boundary-side seeds, benchmark overdispersion
conditions, stale-CLAUDE.md fix.

**Naming (no Greek):** the fitted gDNA strand overdispersion is the code/config/sim variable
**`gdna_strand_overdispersion`** — the intra-class correlation in `[0, 1)`, the exact twin of the
existing `rna_strand_overdispersion`. The math symbol `ρ_d` appears only in formulas below; no
`rho`/`kappa`/`rho_d` in any identifier.

**Resolved decisions:** (1) **single-pass fit** (§1.2), iterative refinement (§1.3) only as a
data-gated backstop. (2) **one global** `gdna_strand_overdispersion` scalar per library.
(3) **`strand_loglik` form — exact-vs-normal — settled by the numerical study in §2.1** (hybrid:
exact below a small-count threshold, normal-moment above). (4) naming as above.
**Problem:** the gDNA strand Beta-Binomial overdispersion `ρ_d` — a documented cornerstone
(`docs/caljointmodel/03_inference.md` §5.2, `01_generative_model.md` §4.3) — is **fitted
nowhere and applied nowhere** on current `main`. The acyclic redesign collapsed the
two-component Beta-Binomial strand model to a single normal-moment term and defaulted its
overdispersion to `0` (the Binomial limit). The synthetic benchmark uses `strand_kappa=None`
(50/50 Binomial gDNA) so it cannot see the gap. This document restores the fit, the
application, the simulator support, and the tests.

This is a **calibration-decode** fix (independent of the locus-EM line in 01/02 — that line
concluded *no* EM strand term; this is the separate, real gap on the gDNA side).

> **Scope note.** RNA strand stays **Binomial** by deliberate decision
> (`strand_balance.py`; `rna_strand_overdispersion` is QC-only). This doc only restores the
> **gDNA** Beta-Binomial. The asymmetry (gDNA BB, RNA Binomial) is intentional and will be
> documented in code (§5).

---

## 1. Fit `ρ_d` from seed nodes, breaking the circularity

### 1.1 Seed nodes (not intergenic-only)

The fit must use the **count-observable** nodes already identified by
`density_model.count_observable_masks` — the same seeds the density estimator trusts:
- **regions with no exon bit** → intergenic + intron-only (`region_count_observable`);
- **boundaries with no *shared* exon** → exon–intron and exon–intergenic seams
  (`boundary_count_observable`).

Intergenic-only is insufficient: under **hybrid capture**, off-target intergenic gDNA is
depleted, so intergenic regions carry few fragments; the intronic and seam nodes carry the
bulk of the observable gDNA and must contribute to the fit.

Each seed node supplies a gDNA-eligible unspliced strand split `(k⁺, n)` from the substrate
(`substrate.contained.n_unspliced_pos/neg` for regions; `substrate.left/right.n_unspliced_pos/neg`
for boundary sides). Sense is oriented by the node's transcript strand where defined
(`region_arrays.strand_class`); intergenic nodes have no strand, but the **overdispersion of
the split around ½ is orientation-independent**, so an arbitrary orientation is fine for the
fit (we measure variance of `k⁺/n` about ½, not a signed skew).

### 1.2 The circularity, and how to break it

A seed node's strand split is a **mixture**: gDNA (rate ½, overdispersed by `ρ_d`) + RNA
(nascent at introns/seams; rate `κ_rna`, stranded). Isolating the gDNA overdispersion needs
the per-node gDNA fraction `π_g` — which the strand deconvolution produces, and which (with
this fix) depends on `ρ_d`. Circular.

**Break it with the count⊥strand conditional independence the engine already relies on.**
The count clue produces a per-node gDNA fraction **independent of `ρ_d`**:
`π_g = clip(density·eff_len / mass, 0, 1)` (from `node_gdna_density`), where the density was
cleaned by the **mean-based** closed-form `gdna_frac = (sense_frac − κ_rna)/(½ − κ_rna)`
([calibrate.py:104](../../src/rigel/calibration/calibrate.py#L104)) — which uses only the gDNA
strand **mean (½)** and `κ_rna`, never `ρ_d`. So `π_g` is a clean, `ρ_d`-free mixing weight.

Given `π_g` (count), `κ_rna` (spliced-channel `StrandModel`, `ρ_d`-free), fit `ρ_d` by a
**pooled method of moments** over seed nodes `s` (implemented in
`calibration/gdna_strand.py::fit_gdna_strand_overdispersion`):

```
mean_s        = ½·π_g,s + κ_rna·(1 − π_g,s)                       # mixture sense rate
excess_var_s  = (k⁺_s − n_s·mean_s)²  −  n_s·mean_s·(1 − mean_s)  # observed − Binomial
gdna_var_s    = (π_g,s·n_s)·(π_g,s·n_s − 1)·¼                     # BetaBinom excess-var scale
gdna_strand_overdispersion = Σ_s excess_var_s / Σ_s gdna_var_s   # pooled, floored at 1e-6
```

`ρ_d` is identified from the **excess variance of `k⁺` beyond Binomial, after removing the
count-derived mean** — a single pass, fully acyclic.

**Identifiability gate (added in Phase 2).** A between-node dispersion cannot be estimated from
a handful of nodes, and a single node's squared residual is a heavy-tailed draw that will
masquerade as overdispersion. So we fit a positive `ρ_d` only when **(a)** there are at least
`_MIN_SEED_NODES = 30` seed nodes (the CLT rule-of-thumb, also what the Gaussian SE below needs)
**and (b)** the pooled excess variance exceeds Binomial sampling noise by `≥ _SIGNIFICANCE_Z = 2`
standard errors (`SE = √(Σ 2·(N·μ(1−μ))²)`, the SE of the excess-variance sum under the
`od = 0` null). Otherwise we default to the **Binomial limit `od = 0`** (`fallback_used`), which
makes the decode byte-identical to the pre-overdispersion behaviour. On real genomes (thousands
of seed nodes) a true overdispersion clears the gate by a wide margin; tiny test scenarios and
50/50 gDNA floor to 0. **These two constants are new heuristics — flagged for review per the
no-magic-numbers policy** (`30` = standard CLT threshold; `2` = standard 2σ significance). **Why MoM, not the MLE:** it is
closed-form (`O(n_seed_nodes)`, no optimizer/convolution), robust (pooling stabilizes
noisy/negative per-node terms), and — crucially — uses the **same variance decomposition the
decode applies** (`(n·gdna_frac)²·¼·od`, §2.1), so the fit and the application are mutually
consistent rather than two different approximations. Pooled denominator ≤ 0 ⇒ no identifiable
gDNA strand signal ⇒ `fallback_used`, floored to Binomial. **Verified:** recovers the planted
overdispersion across `{0.01, 0.05, 0.1, 0.2}` for pure-gDNA and RNA-contaminated seeds
(`tests/calibration/test_gdna_strand.py`).

### 1.3 Optional iterative refinement (only if needed)

The count-based `π_g` is noisier at intron/seam nodes (nascent contamination) than at
intergenic ones. If validation (§4) shows the single-pass `ρ_d` is biased, add a short
outer loop:

```
ρ_d ← 0
repeat (≤ 3×):
    deconv with current ρ_d → refined per-node π_g (joint count×strand)
    ρ_d ← refit from seeds using the refined π_g
until |Δρ_d| < tol
```

This re-introduces a *small, bounded* iteration only for `ρ_d` (everything else stays
acyclic). **Default to the single pass (§1.2);** adopt iteration only if the tests demand it.
The exit decision is data-driven, not assumed.

### 1.4 New module + result field

- **New `calibration/gdna_strand.py`** (Phase 1, done), parallel to `strand_balance.py`. Pure
  estimator + model, trivially testable:
  - `GdnaStrandModel(gdna_strand_overdispersion, n_seed_nodes, n_seed_fragments, fallback_used)`
    with `.beta_concentration()` (`a = ½(1−od)/od`).
  - `fit_gdna_strand_overdispersion(sense, total, gdna_weight, rna_sense_frac) -> GdnaStrandModel`
    — the pooled MoM over seed-node arrays.
  - `fallback_used = True` ⇒ no identifiable gDNA strand signal ⇒ `_OVERDISPERSION_FLOOR` (1e-6,
    Binomial), not silent.
- **Phase 2 wrapper** (in `calibrate.py`): extract the seed-node `(sense, total, gdna_weight)`
  arrays from the substrate / `region_arrays` / `node_density` (count-observable regions +
  boundary sides; `gdna_weight` = the count-clue gDNA fraction) and call the estimator.
- Add `gdna_strand_overdispersion: float` to `CalibrationResult` and to `summary.json`
  (observability — the value must be inspectable, unlike today).

## 2. Apply the gDNA Beta-Binomial in the decode

### 2.1 Correct the strand likelihood variance

`strand_likelihood.strand_loglik` currently applies a single shared overdispersion to the
whole mixture variance (and is fed `0`). Replace with overdispersion on the **gDNA fraction
only** (RNA is Binomial):

```
p   = ½·gdna_frac + κ_rna·(1 − gdna_frac)                 # mixture mean (unchanged)
var = N·p·(1 − p)            (Binomial mixture variance)
    + (N·gdna_frac)² · ¼ · ρ_d   (excess from the shared gDNA strand rate, Beta(½,ρ_d))
loglik(gdna_frac) = −½·(k⁺ − N·p)² / var − ½·log var
```

Limits (correctness checks): `gdna_frac → 1` ⇒ `var → N·¼·(1 + N·ρ_d)` = Beta-Binomial(½, ρ_d)
(normal approx); `gdna_frac → 0` ⇒ `var → N·κ(1−κ)` = Binomial RNA; `ρ_d → 0` ⇒ today's
Binomial mixture exactly (regression-safe). The `(N·gdna_frac)²` scaling correctly localizes
the excess variance to the gDNA fragments, so a node that is mostly RNA is not spuriously
widened.

**Resolved split (exact vs normal-moment) — settled by `scripts/debug/strand_ll_study.py`:**
- **The fit (§1) uses the exact Beta-Binomial** log-pmf. It is evaluated only at the seed
  nodes' observed counts during a 1-D optimization over the overdispersion — *not* on a
  per-grid-point grid — so its cost is bounded and exactness keeps the estimate unbiased.
- **The decode (§2) uses the corrected normal-moment everywhere — no small-`n` fallback.** The
  study compared the normal-moment posterior-median `gdna_frac` against the exact two-component
  mixture marginal across `n ∈ {8…500}` and overdispersion `∈ {0.05, 0.1, 0.2}`: the discrepancy
  is **small and ~constant in `n`** (max ≈ 0.03; ≤ 0.01 at realistic overdispersion ≤ 0.1), *not*
  a small-`n` effect — because an overdispersed Beta-Binomial's *fraction* converges to a Beta,
  so the Normal-vs-Beta shape error persists at all `n` and a threshold buys nothing. The bias
  grows with overdispersion and is in the **FP-safe direction** (slightly over-calls gDNA for
  pure-gDNA nodes). The normal-moment already *is* Beta-Binomial-aware (the variance is inflated
  by the gDNA term); the residual is the moment-matching shape error, acceptable at the
  calibration's ±2–3 % accuracy floor. (If a future need demands it, the cheaper fix than the
  O(n²) exact mixture is a Beta-fraction likelihood; not needed now.)

### 2.2 Thread `ρ_d` through

- `calibrate.py`: after `fit_gdna_strand_dispersion`, pass `ρ_d` to `deconv_regions` /
  `deconv_sides` (which forward `strand_overdispersion=ρ_d` to `strand_loglik`). This is the
  one-line wiring that is currently missing (the `=0.0` defaults at
  [joint_deconv.py:137,172](../../src/rigel/calibration/joint_deconv.py#L137) are the bug).
- The mean-based closed-form gDNA-fraction cleaning (calibrate.py:104) stays **unchanged**:
  `ρ_d` does not bias the mean (½), so the density cleaning is already correct; `ρ_d` only
  affects the per-node strand *likelihood confidence* (whether a noise-skewed gDNA node is
  held as gDNA rather than mis-called RNA — the 126/26 motivation).

## 3. Simulator — gDNA strand overdispersion

The machinery exists: `GDNAConfig.strand_kappa` + `_init_strand_regions` draw per-region
`Beta(κ/2, κ/2)` sense rates ([reads.py:251](../../src/rigel/sim/reads.py#L251)). The
intra-class correlation is **`ρ_d = 1/(κ + 1)`** (so `κ = (1 − ρ_d)/ρ_d`); `κ → ∞` ⇒ `ρ_d → 0`
(Binomial), small `κ` ⇒ strong overdispersion.

- Expose **`gdna_strand_overdispersion`** in `GDNAConfig` and the suite config
  (`scripts/sim/configs/*.yaml`, `rigel.sim.suite`/`whole_genome`) as a swept axis; convert to
  the existing `strand_kappa` Beta concentration **internally** (`strand_kappa = (1 −
  overdispersion)/overdispersion`), keeping one user-facing knob in the clear units.
- Add benchmark conditions at `gdna_strand_overdispersion ∈ {0, 0.05, 0.2}` crossed with the
  existing gDNA-level × strand-spec × capture grid (at least for `gdna_high`).
- Record the *true* `gdna_strand_overdispersion` in the condition metadata / `truth_summary.json`
  so tests can compare fitted vs true.

## 4. Tests — recover `ρ_d` and prevent the failure mode

**Unit (the estimator, `tests/calibration/`):**
1. `fit_gdna_strand_dispersion` on **synthetic pure-gDNA seeds** (draw per-node `(k⁺, n)` from
   `BetaBinom(n, ½, ρ_d)` at known `ρ_d`): recovered `ρ_d` within tolerance across a grid
   `ρ_d ∈ {0.01, 0.05, 0.1, 0.2}` and node counts/depths. **This is the core "fit the BB
   correctly" test you asked for.**
2. **Mixture identifiability:** seeds drawn from the `π_g·BetaBinom + (1−π_g)·Binom(κ_rna)`
   mixture with known `π_g`, `κ_rna`, `ρ_d` → recovered `ρ_d` within tolerance; check bias as
   nascent fraction `(1−π_g)` grows (validates the §1.2 break, and whether §1.3 iteration is
   needed).
3. **Thin-seed fallback:** few/empty seeds → `fallback_used`, `ρ_d = _BB_FLOOR`, no crash.

**Property / regression:**
4. `ρ_d → 0` ⇒ `strand_loglik` and the full decode are **byte-identical** to today (guards the
   Binomial limit; protects the deliberate RNA-Binomial path).
5. `strand_loglik` variance limits (§2.1) hold at `gdna_frac ∈ {0, 1}`.

**Integration (sim → calibrate):**
6. Simulate at known `ρ_d` (via `gdna_strand_kappa`) → calibrate → **fitted `ρ_d` ≈ true**
   across the grid (the end-to-end recovery test).
7. **The 126/26 motivation:** a pure-gDNA region with a noise-driven strand skew must be
   deconvolved as gDNA (not mis-called RNA) once `ρ_d > 0`, and *would* be mis-called under
   `ρ_d = 0`. Direct regression on the failure the BB exists to prevent.

**Benchmark:**
8. Add the overdispersed-gDNA conditions to the calibration-benchmark skill; assert (a) no
   regression on `ρ_d = 0` conditions, (b) improved gDNA/RNA separation on overdispersed
   conditions vs the `ρ_d = 0` baseline (pool fraction + the strand-skew mis-call rate).

## 5. Code clarity / reorganization (explicit ask)

- **`gdna_strand.py`** (new) mirrors `strand_balance.py`: the gDNA strand model lives in one
  obviously-named place, symmetric to the RNA strand model. `GdnaStrandModel` carries `ρ_d`.
- **`strand_likelihood.py`**: rename the parameter `strand_overdispersion →
  gdna_strand_overdispersion`; rewrite the docstring as the **two-component (gDNA Beta-Binomial +
  RNA Binomial)** model, not "single overdispersion." Make the gDNA-only excess-variance term and
  its derivation explicit.
- **One place documents the asymmetry**: gDNA = Beta-Binomial(½, `ρ_d`) *fitted*; RNA =
  Binomial (deliberate, `rna_strand_overdispersion` QC-only). Put this side-by-side in the
  `calibrate.py` module docstring and the two model modules so the design intent is
  unmissable.
- **Fix stale `CLAUDE.md`**: it claims `mstep.py` fits `ρ_d_bb`; `mstep.py` was removed in the
  acyclic redesign. Update the architecture section to the acyclic fit (`gdna_strand.py`).
- **Result/observability**: `gdna_strand_overdispersion` surfaced in `CalibrationResult` +
  `summary.json` so it is never again silently 0.

## 6. Phasing, risks, gates

- **Phase 1 — simulator + tests first (TDD).** Add `gdna_strand_kappa` to the sim/suite and
  the unit tests (§4.1–4.5). These are independent of the fit and define "correct."
- **Phase 2 — fit + apply.** `gdna_strand.py` (§1), the `strand_likelihood` variance (§2.1),
  the wiring (§2.2), result field. Land §4.6–4.7 integration tests.
- **Phase 3 — clarity + benchmark.** §5 cleanup; §4.8 benchmark conditions; CLAUDE.md.
- **Decision gate after Phase 2:** if §4.2 shows the single-pass fit is biased by nascent
  contamination, enable the §1.3 iteration; otherwise keep the single pass.

**Risks:**
- *Seed contamination by nascent RNA* at introns/seams biases `ρ_d` high (extra variance from
  the stranded nascent mean spread). Mitigated by the count-`π_g` mixture weighting (§1.2);
  measured by §4.2; iteration (§1.3) as backstop.
- *Normal-moment approximation at small `N`* — guarded by §2.1 limits + the exact-BB fallback
  for small `N`; decided by §4 tests.
- *No new magic numbers*: `ρ_d` is fitted; `κ_d = ½` is biology; `_BB_FLOOR = 1e-6` is the
  documented numerical floor. The simulator `ρ_d` grid is test config, not product code.

## 7. Decisions (resolved 2026-06-06)

1. **Single-pass fit** (§1.2); the §1.3 iteration is a data-gated backstop (gate after Phase 2,
   driven by test §4.2). ✔
2. **One global** `gdna_strand_overdispersion` scalar per library (matches the design). ✔
3. **`strand_loglik` form — settled by the §2.1 study:** fit = exact Beta-Binomial; decode =
   corrected normal-moment **everywhere, no fallback** (the normal-vs-exact discrepancy is
   small, ~constant in `n`, and FP-safe — a small-`n` threshold buys nothing). ✔
4. **Naming**: `gdna_strand_overdispersion` everywhere (twin of `rna_strand_overdispersion`); no
   Greek/`rho`/`kappa` identifiers. ✔
