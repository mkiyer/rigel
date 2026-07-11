# Enriched-mode reference-density detector (eff-length contraction)

**Purpose:** set `ρ_ref` — the fully-captured gDNA density that anchors the capture eff-length contraction
(`w_n = min(ρ_n/ρ_ref, 1)`, `eff = Σ Sₙwₙ`) — from **observed node densities**, unbiased, with **no
assumption about probe locations**. Companion to `efflen_shared_reference_fix_plan.md` (which proves a
*shared* reference makes `eff(nascent) ≥ eff(mature)` — no inversion). This doc supplies the missing piece:
**how to detect the reference from the data.**
**Validated:** `scratchpad/kde_mode_scan.py` across all 16 `quick_3to1_5mb` scenarios (plot in that dir).

## Why a shared, data-detected reference

- **Stops the inversion & siphon.** A per-transcript reference (`ρ* = G_c/E_c`) contracts on *within-transcript*
  density variation — including noise/false gDNA — which breaks monotonicity and (confirmed) fires even in
  `gNONE`: 40% of transcripts contracted `<0.9×` with **no real gDNA**. A single shared reference makes the
  contraction depend on each node's density vs one library-wide level, guaranteeing `nascent ≥ mature`.
- **gDNA is stable, so a global reference is legitimate.** gDNA copy number varies a few-fold at most (even
  in cancer); RNA varies ~10⁴-fold. So "the fully-captured gDNA density" is a real, stable, library-wide
  constant — a global `ρ_ref` is principled for gDNA and would be absurd for RNA.
- **No probe assumptions.** Under capture the per-node gDNA density is **bimodal** — a depleted (off-target)
  mass and an enriched (on-target) mode. `ρ_ref` = the **enriched mode**, read from the density distribution
  itself, never from probe annotations.

## The method

1. Per-node gDNA density `ρ_n = f_g,n · mass_n / eff_gdna,n`; per-node weight = deconvolved gDNA mass
   `g_n = f_g,n · mass_n`.
2. Build a **MASS-WEIGHTED** KDE of `log ρ_n`. *Mass-weighting is the key idea:* a small captured panel
   (e.g. 50 of thousands of genes) is a tiny bump in the **count** distribution but — because enriched nodes
   carry ~100× the mass — a **dominant** peak in the **mass** distribution. So the enriched mode is prominent
   and detectable even for a small panel (the concern that motivated this).
3. Find local maxima. `ρ_ref` = the **rightmost significant peak** (the enriched mode).
4. **Unimodal ⇒ no enrichment.** If there is one mode (capture-off, or no gDNA), `ρ_ref` = that mode ⇒
   `ρ_n ≈ ρ_ref` for all n ⇒ `w ≈ 1` ⇒ **no contraction**. This is the correct, automatic fallback.

```python
x, w = log(rho[rho>0]), gmass[rho>0]           # log density, weighted by gDNA mass
km = weighted_kde(x, w)                         # MASS-weighted (enriched mode prominent)
peaks = local_maxima(km)
rho_ref = exp(x_grid[ rightmost significant peak ])   # unimodal -> the single mode
```

## Validation — 16 scenarios (`kde_mode_scan.py`)

| bucket | modes | depleted | enriched (ρ_ref) | ratio |
|---|--:|--:|--:|--:|
| g300 capOFF (ss99/ss50) | 1 | 0.60 | 0.60 | **1.0× (no contraction ✓)** |
| **g300 capON** (ss99/ss50) | 2 | 0.013 | **15–16** | **~1150–1220×** |
| gNONE capOFF | 1 | ~5e-4 | ~5e-4 | 1.0× (no contraction ✓) |
| gNONE capON | 2–4 | tiny | tiny | 7–200× (FALSE gDNA — see below) |

- **capture-OFF → unimodal → ratio 1.0× → no contraction.** Exactly right.
- **capture-ON (real gDNA) → clean bimodal**; the enriched mode is a small *count* peak but the dominant
  *mass* peak (mass-weighting recovers it); `ρ_ref` ≈ the true enrichment level, ratio ≈ the enrichment
  factor. Exactly right.
- **`gNONE capON` → capture-patterned FALSE gDNA** (bimodal, noisy). The calibration mis-deconvolves some
  capture-enriched RNA as gDNA; this is the root of the `gNONE` nascent siphon. The detector exposes it.

## Answers to two design questions

- **"Max density ⇒ all nodes except the max shrink?"** Yes — with `ρ_ref = max`, every `w = ρ_n/ρ_max ≤ 1`.
  But **don't use the raw max**: it's a one-node, high-variance estimator (a single outlier sets the whole
  library's reference). Use the KDE **enriched mode** — the robust version of "the top."
- **"Allow expansion (w>1) as well as shrinkage?"** No. The clamp `w = min(·,1)` keeps `eff = Σ Sₙwₙ ≤ fl` —
  a transcript can't be effectively *longer* than its full FL-marginal length (that's >100% capture,
  meaningless), and un-clamped it degenerates (all densities → `ρ_ref`). The clamp is also *why the mode
  beats the max*: with the enriched-mode reference, every fully-captured node reads `w=1` (no spurious shrink
  from an outlier), and only genuinely depleted nodes shrink — the meaningful (depleted-tail) distribution is
  preserved, the top is correctly flattened to "fully captured."

## Wiring plan (straightforward)

1. **Compute the global `ρ_ref` once** — mass-weighted KDE enriched mode over all node gDNA densities
   (reuse the Phase-2 `GdnaDensityPrior` KDE machinery; add the mode-detection).
2. **`transcript_capture_eff_lengths`:** replace the per-transcript `ρ* = G_c/E_c` with the global
   `ρ_ref`; `w_n = min(ρ_n/ρ_ref, 1)`; `eff_em = Σ Sₙwₙ` (clamped, ≤ fl). Keep the junction seams (shipped)
   and the multimapper-evidence shrinkage.
3. **Remove the per-transcript reference.**
4. **Global gDNA-abundance guard:** if the detector is unimodal (or the library gDNA fraction ≈ 0), skip
   contraction (`eff = fl`). This enforces "no gDNA / no capture ⇒ no contraction" globally and neutralises
   the `gNONE` false-gDNA contraction.
5. **A/B** on the soft 3-pool surplus + the nascent/mature inversion rate + the `gNONE` siphon.

## Open (basic detector; refine only if needed)

The rightmost-significant-peak rule uses a mild prominence guard (5% of the max peak height). It is clean on
the target (g300 capON: 2 well-separated modes); the `gNONE capON` false-gDNA case is noisier (2–4 modes),
but that regime is handled by the abundance guard (step 4), not by the reference. If robustness is needed
later, replace peak-finding with a 2-component mixture + separation test — but for the flagship it is not
required.
