# Implementation plan — RNA strand Beta-Binomial overdispersion (symmetric with gDNA)

**Status: IMPLEMENTED (2026-06-07).** Companion to the shipped gDNA strand BB work
([03](03_gdna_strand_overdispersion.md), [04](04_overdispersion_prior_and_shrinkage.md)). Both
components are now Beta-Binomial with the **same default prior** (`rna_strand_prior_alpha_beta = 3`,
weight 30, matching gDNA — symmetry chosen over FP-tightness). Shipped: shared component-agnostic
MoM core `gdna_strand._fit_overdispersion(node_mean, component_frac, component_mean)`;
`fit_rna_strand_overdispersion` / `fit_rna_strand_from_substrate` (boundary-side spliced);
`strand_loglik` RNA variance term scaled by `κ(1−κ)` (the RNA component's μ(1−μ), **not** ¼ — a real
bug caught in testing: a ¼ scale biases the fit by 4κ(1−κ)); `CalibrationResult
.rna_strand_overdispersion`; `calibrate.py` fit + decode wiring + a boundary-spliced-mean-vs-κ
diagnostic; `summary.json` field; 23 tests in `tests/calibration/test_rna_strand.py`.

**Verified:** unstranded balanced node 0.319 (asymmetric, old) → 0.500 (symmetric); graded
balanced-node gdna_frac 0.50 → 0.81 as κ 0.5→0.99; RNA overdispersion recovered at known levels.
**Measured cost (§9):** in the SS=0.65 clean-sim the junction double-counting inflates the fit to
`od≈0.092` (sim has no RNA overdispersion), softening the strand clue and costing ~6 frags of
weak-SS accuracy (tolerance widened 30→40, documented in-code). Tracked as future minor work (count
one side per boundary per fragment).

> **Note (revises §9):** the double-counting **over**-estimates the overdispersion here (a
> multi-exon fragment credits several sides with the *same* strand → positive between-side
> correlation), the opposite of §9's original "under-estimate" guess.

## 1. Motivation — the asymmetry is real and harmful (verified)

The decode currently models **gDNA strand as Beta-Binomial** (fitted overdispersion) but **RNA
strand as Binomial** (overdispersion forced to 0). That asymmetry makes the strand likelihood
**spuriously informative on unstranded data**, where it must be flat. Measured on
`strand_likelihood.strand_loglik` (κ = `rna_sense_frac` = 0.5, a balanced node, 50/100 sense):

| model | median `gdna_frac` (should be 0.5) |
|---|---|
| **current** (od_gdna = 0.143, od_rna = 0) | **0.317** ← spurious pull toward RNA |
| symmetric BB (od_gdna = od_rna = 0.143) | **0.499** ← uninformative restored |
| no overdispersion (both 0) | 0.499 |

The mechanism: the likelihood's `−½·log(var)` term prefers the **lower-variance** component;
with only gDNA overdispersed, an unstranded node looks "tighter" as RNA, so `gdna_frac` is pulled
down (loglik varies 1.36 nats across the grid when it should be flat). Making RNA Beta-Binomial
with a comparable overdispersion **symmetrizes the variance and restores uninformativeness**
(0.499). This is the user's "unstranded ⇒ uninformative; weakly→strongly stranded ⇒ graded
information" requirement — currently violated, fixed by RNA BB.

> **Verification status:** the *unstranded* case is verified above. "Weakly stranded" (κ ≈
> 0.6–0.8) and "near-perfect" (κ ≈ 0.99) graded behavior should be added as tests (§7).

## 2. Current RNA strand handling (code map)

- **`calibration/strand_balance.py`** — `fit_strand_balance(StrandModels) -> StrandBalance`:
  - `rna_sense_frac` = posterior mean `(n_same+1)/(n_obs+2)` of the **BAM-scan** `StrandModel`'s
    2×2 spliced-unique-mapper counts (global, per-library). **Used by the decode** (centers the
    strand likelihood; cleans the count density).
  - `rna_strand_overdispersion` = `1/(n_obs+3)` — a **posterior-predictive thin-count**
    overdispersion (statistical-power diagnostic), **QC-only, deliberately not in the decode**.
- **`calibration/strand_likelihood.py::strand_loglik`** — the decode's per-node strand term:
  `var = N·p(1−p) + (N·gdna_frac)²·¼·gdna_strand_overdispersion`, `p = ½·gdna_frac +
  rna_sense_frac·(1−gdna_frac)`. **No RNA overdispersion term.**
- **`calibrate.py`** — calls `fit_strand_balance`, uses `rna_sense_frac` for density cleaning and
  passes it to `deconv_regions`/`deconv_sides`; at unstranded (`|0.5−rna_sense_frac|<1e-6`) the
  *density cleaning* falls back to count-only, but `strand_loglik` is **still added** for POS/NEG
  regions (hence the §1 spurious pull).
- **`CalibrationResult`** carries `rna_sense_frac` + `gdna_strand_overdispersion`; **no
  `rna_strand_overdispersion`** field yet.

> **⚠ Prior decision to respect (PR 9, documented in `strand_balance.py`):** wiring the RNA BB
> widening into the decode was previously **measured to be negligible (rel ≤ 1e-3) and to slightly
> *worsen* the silent-gene false-positive rate** (it softens the strand clue), so it was kept out
> **by design**. **But that test used the posterior-predictive `1/(n_obs+3)` quantity applied when
> gDNA was *also* Binomial.** The situation is now different: (a) gDNA BB is live, *creating* the
> §1 asymmetry that RNA BB would cancel; (b) the proposed RNA overdispersion is a **fitted
> between-region biological** quantity, not the thin-count posterior-predictive. So PR 9's
> conclusion must be **re-evaluated**, not assumed — see §7 (the silent-gene FP benchmark is the
> gate).

## 3. The proposed fit — RNA overdispersion from boundary-side spliced counts

Spliced fragments are **pure RNA**, so their per-boundary-side sense/antisense split directly
measures RNA strand overdispersion. The data already exists, no new structures:

- **Source:** `substrate.left` / `substrate.right` (per-region boundary-side views) expose
  `n_spliced_sense` and `n_spliced_antisense` (motif-relative, already oriented — valid even in
  AMBIG regions). Populated by the accumulator: a spliced "jump" credits boundary-side flux in the
  spliced channels (`accumulator.cpp:187`, `flux_*[ch] += 1`).
- **Double-counting (accepted):** a fragment spanning K junctions credits ~K boundary sides
  (`flux += 1` per crossed slice). This inflates the effective count and mildly *under*-estimates
  the overdispersion (correlated repeats look more Binomial), but needs no new data structures and
  is a reasonable first approximation. **Document it; revisit only if the fit looks biased.**
- **Estimator:** the **same pooled MoM** as gDNA, but for the RNA component:
  - per side `s`: `sense_s = n_spliced_sense`, `n_s = n_spliced_sense + n_spliced_antisense`;
  - `mean = rna_sense_frac` (the existing global κ — keep the *mean* from the reliable BAM-scan
    `StrandModel`; fit only the *overdispersion* from the between-side spliced variance);
  - `excess_var_s = (sense_s − n_s·κ)² − n_s·κ(1−κ)`; `var_scale_s = n_s(n_s−1)/4` (pure RNA ⇒
    component fraction 1); `od_rna = Σ excess / Σ scale`, shrunk toward an RNA prior (§5).

## 4. Generalize the estimator (the one real refactor)

`gdna_strand.fit_gdna_strand_overdispersion` currently hard-codes the gDNA roles: `mean_rate =
½·weight + κ·(1−weight)` and `var_scale = (weight·N)²·¼` both use `weight = gdna_weight`. For RNA
the component mean is κ (not ½) and the component fraction is 1 (not `gdna_weight`). **Generalize
to a component-agnostic estimator:**

```
fit_strand_overdispersion(sense, total, *, component_mean, component_frac,
                          prior_overdispersion=0, prior_weight=0)
    excess_var  = (sense − total·node_mean)² − total·node_mean·(1−node_mean)   # node_mean passed in
    var_scale   = (component_frac·total)·(component_frac·total − 1)·¼
    od          = shrink(Σ excess / Σ var_scale, prior, n_seed_nodes)
```

- **gDNA** call: `node_mean = ½·gdna_weight + κ·(1−gdna_weight)`, `component_frac = gdna_weight`.
- **RNA** call: `node_mean = κ`, `component_frac = 1` (pure spliced).

Rename `gdna_strand.py` internals to be component-neutral (or keep `gdna_strand.py` + add a thin
`fit_rna_strand_from_substrate` that calls the shared core). Either way the shrinkage, clamp, and
`overdispersion_for_beta` helper are reused unchanged. **This is the bulk of the new code** and is
small (~mirrors the gDNA wrapper, reading the spliced channels of `substrate.left/right`).

## 5. Apply it in the decode — the two-component variance

`strand_loglik` gains the symmetric RNA term:

```
var = N·p·(1−p)
    + (N·gdna_frac)²·¼·gdna_strand_overdispersion
    + (N·(1−gdna_frac))²·¼·rna_strand_overdispersion        # NEW
```

`od_rna → 0` recovers today's behavior (regression-safe). §1 shows equal od restores unstranded
uninformativeness (median 0.499).

> **Subtlety to flag (not a blocker):** the per-component excess-variance decomposition
> (`g²·od_g + (1−g)²·od_r`) is exact when the two components have well-separated rates (stranded:
> gDNA ½ vs RNA κ) but is an approximation when their rate *distributions* coincide (unstranded +
> equal od). Empirically it is near-flat there (median 0.499, §1), so it is good enough; an exact
> mixture-variance form is a possible future refinement, not needed for v1.

## 6. Wiring (all mirrors the gDNA path)

- **`CalibrationConfig`** (config.py): add `rna_strand_prior_alpha_beta` (default — see below) +
  `rna_strand_prior_weight`, symmetric to the gDNA knobs, with the same `≥ 2.0` validation and the
  `od = 1/(2a+1)` conversion.
- **`calibrate.py`**: after `fit_strand_balance`, fit `rna_strand_overdispersion` from the boundary
  spliced counts (§3/§4) with the config prior; pass it into `deconv_regions`/`deconv_sides`
  (alongside `gdna_strand_overdispersion`); add to the `CalibrationResult` + debug log.
- **`joint_deconv.py`**: thread `rna_strand_overdispersion` through `deconv_regions`/`deconv_sides`/
  `_joint_per_node` → `strand_loglik` (exactly parallel to the existing `gdna_strand_overdispersion`
  param).
- **`CalibrationResult`** (result.py): add `rna_strand_overdispersion: float` field + `[0, ceil]`
  validation.
- **`cli.py`** `summary.json`: add `rna_strand_overdispersion` to the calibration block.
- **Default prior `a` for RNA:** unlike gDNA (where high overdispersion is FP-safe), an
  over-wide RNA prior *softens* the strand clue (the PR 9 risk). So the RNA prior should be
  **weaker/tighter** (larger `a`, smaller od₀) and/or a larger `prior_weight` is *not* the lever —
  pick `a` so od₀ is small (e.g. `a` giving od₀ ≈ the typical fitted RNA overdispersion). **Decide
  the RNA default `a` from the benchmark (§7), not by symmetry with gDNA's `a=3`.**

## 7. Validation (the gate is the silent-gene FP rate — PR 9's concern)

1. **Unit — unstranded uninformative:** `strand_loglik` at κ=0.5 with `od_gdna = od_rna` ⇒ median
   `gdna_frac ≈ 0.5` and near-flat loglik (the §1 check, as a test).
2. **Unit — graded information:** κ ∈ {0.5, 0.7, 0.9, 0.99} ⇒ strand information increases
   monotonically; unstranded ≈ 0.
3. **Unit — RNA overdispersion recovery:** boundary-side spliced counts drawn at known od_rna ⇒
   recovered (mirrors the gDNA recovery test); `od_rna → 0` ⇒ `strand_loglik` byte-identical.
4. **Simulator:** the sim already supports RNA `strand_specificity`; add an RNA-strand-
   overdispersion knob if we want graded-od conditions (the gDNA per-region-Beta machinery in
   `whole_genome.py` can be reused for RNA). *Optional for v1.*
5. **Benchmark — the decision gate:** re-run the calibration benchmark and confirm RNA BB **does
   not worsen the silent-gene / gdna_none false-positive rate** (the explicit PR 9 finding). If it
   does, tighten the RNA prior (§6) or keep RNA BB gated behind a config flag (default off) until
   the FP concern is resolved.
6. **Goldens:** regenerate (the new RNA term shifts the decode); confirm shifts are small +
   directional.

## 8. Effort / risk

**Medium — mostly parallel to the shipped gDNA work**, since the estimator, shrinkage, config
pattern, decode-threading, and fit-from-substrate wrapper all already exist and just need to be
made component-agnostic + an RNA path added. Concretely:

| piece | size | risk |
|---|---|---|
| Generalize estimator (component_mean/frac) + RNA substrate wrapper | small–med | low |
| `strand_loglik` RNA variance term | tiny | low (od→0 = no-op) |
| Config + result + summary + joint_deconv threading | small | low (mirror gDNA) |
| Tests (unstranded, graded, recovery) | med | low |
| **Benchmark FP re-evaluation (PR 9 gate)** | — | **the real risk** — RNA BB softens the strand clue; must confirm it doesn't regress silent-gene FP |

**Recommendation:** implement behind the existing config pattern; gate adoption on the §7.5
benchmark. The mean stays from the BAM-scan `StrandModel`; only the overdispersion is new. Land
the unstranded-uninformative win (§1) only if the silent-gene FP rate holds.

## 9. Open questions

1. **RNA prior default `a`** — tighter than gDNA's `a=3` (less FP-safe to over-widen RNA). Set
   from the benchmark.
2. **Double-counting** — accept for v1 (documented), or down-weight multi-junction fragments
   later?
3. **Mean source** — keep `rna_sense_frac` from the BAM-scan `StrandModel` (recommended), or also
   re-fit the mean from boundary spliced counts? (Recommend: keep the existing mean.)
4. **Exact vs approximate mixture variance** at unstranded (§5 subtlety) — approximate is
   empirically fine; revisit only if a test demands it.
