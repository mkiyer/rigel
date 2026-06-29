<!-- title: The enrichment transfer ê(z) — why we fit it, the current monotone-spline approach, its behavior + problems, and fit alternatives -->
# The enrichment transfer ê(z): model-fit deep dive

**Status:** working design doc (2026-06-28). Records *why* calibration fits an enrichment model, *exactly how*
the current monotone P-spline works (parameters + internals), what it does well, where it misbehaves (with
plot evidence), and the candidate fixes — both within the P-spline framework and as alternative fits. Companion
to [`message_state_separation.md`](message_state_separation.md) (the precision architecture) and
[`grid_solver_ideas_and_alternatives.md`](grid_solver_ideas_and_alternatives.md). Saved fit-point datasets for
offline experimentation: `scratchpad/enrichment_fit_points/*.tsv` (sim); real VCaP extraction pending.

## 1. Why fit a model at all — the goal

The calibrator must give every node a **gDNA-density prior** `ρ_g`. By the **count-zero-info principle**
([`CALIBRATION_ARCHITECTURE.md`](CALIBRATION_ARCHITECTURE.md)) the strand likelihood alone cannot solve
gDNA-vs-RNA at an AMBIG node, and at ss≈½ it cannot solve it *anywhere* — so the node leans on a global density
prior + propagation. A single genome-wide constant is wrong under **hybrid capture**: capture concentrates gDNA
on exons (enriched, ρ_g≈16/bp on the flagship) and depletes introns/intergenic (ρ_g≈0.01/bp) — a **40×+ spatial
swing**. So we need a *spatially-varying* gDNA-density prior.

`ê(z)` is that model. It predicts a node's gDNA density from a **clean predictor** `z` that is correlated with
the local enrichment level but is itself relatively gDNA-clean. The transfer is **region ← its flanking
boundaries**: for an exon REGION (the response), `z` = the mean region-facing crossing density of its flanking
**gDNA-clean** boundaries (intron/intergenic↔exon junctions), and the response `ρ_g` = the region's own
strand/spliced-derived gDNA density. (Both clean flanking boundaries are averaged; one if only one is clean.)
The boundary crossing is the predictor because it is RNA-cleaner than the exon interior we are trying to
estimate. *(Open: a boundary-targeted transfer — predicting a boundary's gDNA from its neighbours — is not
modelled; boundaries currently get only propagation + their own solve. Deferred.)*

**The goal for the fit, precisely.** A function `ê(z)` that is:
1. **monotone increasing** — more crossing density ⇒ more gDNA;
2. **anchored at a depleted floor** as `z→0` — no enrichment ⇒ the true low-density baseline (see §4);
3. **smooth + robust** — a trustworthy curve on noisy real data, no fitting artifacts;
4. **self-collapsing** — when `z` carries no real signal, ê reduces to the flat baseline (handled downstream by
   the EB-BIC weight, §3) — capture is never *detected*, the model *adapts*.

## 2. The current approach (precise internals)

### 2.1 The data (the teacher points)
Per single-strand exon "teacher" node, one point `(z, ρ_g, w_fit)`:
- `z` = mean of `{d_right(left boundary), d_left(right boundary)}` over the gDNA-**clean** flanks, where the
  boundary face density `d = M_face / E_face` (`fit_enrichment_transfer`, bp_solver.py).
- `ρ_g` (the response) = a reliability-weighted blend of a **strand-derived** estimate (`f_g·T`, weight
  `(2κ−1)²`) and a **spliced/mature-derived** estimate (strand-immune, weight `1−(2κ−1)²`). At ss≈½ the strand
  weight → 0, so ρ_g is ~entirely the spliced-derived term — and its imperfect mature-subtraction is what leaks
  nascent RNA into ρ_g (§5).
- `w_fit` = the per-point reliability × eff-length (`denom·E_gdna`).

### 2.2 The fit — `variance_model.MonotoneMean.fit`
A **monotone-increasing P-spline in log-log space**:
- **space:** `lx = log z`, `ly = log ρ_g` (both > 0 required; the fit is on `ly`).
- **basis:** cubic B-spline (`degree = 3`), basis size `k = 18` (an *upper bound* on flexibility; λ controls the
  realized smoothness). Knots span `[log z_min, log z_max]`.
- **smoothing λ:** chosen by **GCV** over a fixed grid `λ ∈ logspace(−4, 4, 40)`. For each λ:
  `β = (BᵀB + λP)⁻¹Bᵀy`, effective dof `edf = trace((BᵀB+λP)⁻¹BᵀB)`, `GCV = n·RSS/(n−edf)²`; the λ minimizing
  GCV wins. **One global λ for the whole curve.**
- **monotonicity:** enforced by the reparameterization `β = cumsum([α₀, exp(α₁), exp(α₂), …])` — the spline
  coefficients are a cumulative sum of *non-negative* increments (`exp(·) ≥ 0`), so the fitted curve is
  non-decreasing by construction (`_fit_monotone`, penalized weighted least squares over α).
- **robustness:** 2 iterations of **bisquare-weighted IRLS** (`robust_iters = 2`) on the log-residuals
  (down-weights points with `|r| > 4.685·MAD`).
- **net-balance recalibration:** `scale = Σ(rw·y) / Σ(rw·ŷ)` so the curve's weighted mean matches the data's.
- **predict:** `ê(z) = scale · exp(spline(clip(log z, log z_min, log z_max)))`, **flat-extrapolated** outside the
  fit range.
- **fallback:** `< max(k, 8) = 18` points, or no log-range ⇒ a flat constant (the recal-weighted mean).

### 2.3 Downstream consumption
`ê(z)` is not applied raw. The per-node prior is `ρ̂ = w·ê(z) + (1−w)·ρ_global`, `N = ρ̂²/σ²_prior` (a *humble,
M-independent* hyperprior). `w = w_enrich = σ(½·ΔBIC)` is the EB-BIC slope-evidence weight (the gateless replacement
for the old permutation gate; [`message_state_separation.md`](message_state_separation.md)). The currently-shipped
anchor of that blend is `ρ_global` — which §4 argues is the wrong anchor.

## 3. Behavior — what it does well

- **Monotone + self-collapsing:** ê rises with z; off-capture the GCV fit is ~flat and `w_enrich→0`, so the
  layer reduces to the genome-wide baseline with no detector (verified: off-cap `w_enrich = 0.03`).
- **Strong enrichment is captured cleanly:** flagship log-log fit R² = 0.93 (tight rising relationship).
- **Parameter-free smoothing:** GCV picks λ; no hand-tuned smoothness constant.

## 4. Problem A — the wrong anchor (and the fix)

The blend anchors to **`ρ_global`**, which is the **mass-weighted average gDNA density over all nodes** — and
under capture it is **dragged up by the enriched exons**, so it is *not* a low-density floor. Measured (sim):

| condition | ρ_global | depleted floor (intergenic+intronic regions) | ratio |
|---|---|---|---|
| flagship (capture, ss0.99) | 0.552 | **0.013** | **42×** |
| unstranded (capture, ss0.50) | 0.520 | **0.012** | **44×** |
| off-capture | 0.552 | 0.597 | 0.9× (converge — no enrichment) |
| zero-gDNA + nascent | 0.027 | 0.013 | 2× |

**The fix (proposed).** Anchor the curve's low-z end (and the blend) to a **depleted floor**:
`floor = (Σ gDNA mass over intergenic + intronic REGIONS + pseudocount) / (Σ their gDNA eff-len)`.
- Intergenic + intronic **regions** (not boundaries) are *depleted* under capture, so their density estimates
  the true gDNA floor — the thing a no-enrichment node should fall back to.
- The **pseudocount** on the mass sum is a tiny "there is never truly zero gDNA" prior — it keeps the floor
  finite (and the log-fit stable) when the depleted regions carry ~0 gDNA (the zero-gDNA case → floor ≈
  pseudocount/Σeff, a small positive number, not 0).
- Off-capture this floor ≈ ρ_global (they converge), so the change is inert where there is no enrichment, and
  load-bearing only under capture — exactly the desired uniform behavior.

(Note: this is the floor for the *fit's* low-z anchor / the blend baseline — distinct from, and not to be
confused with, the earlier-refuted idea of replacing the per-node global *mean* with an intergenic-only value,
which collapsed the prior strength `N`. Here it stabilizes the curve and the blend baseline, not the strength.)

## 5. Problem B — the response is contaminated by nascent RNA (ACCEPTED CAVEAT)

The teacher response `ρ_g` (the strand/spliced blend) carries a real, systematic correlation with `z` even on
**zero-gDNA** data (the phantom panel's faint R²≈0.11 slope). This is **not noise** — the sim is exact — it is
**nascent pre-mRNA**: at ss≈½ the strand mixture power `(2κ−1)²→0`, so the strand cannot separate balanced
nascent *unspliced* RNA from balanced gDNA, and the spliced/mature subtraction that builds ρ_g is imperfect
(~12–23% of nascent unspliced mass leaks into ρ_g), correlated with z. **Decision (accepted limitation):** for
**unstranded + nascent** libraries we do **not** disambiguate nascent RNA from gDNA — gDNA is **overestimated**,
and nascent RNA tends to be predicted as gDNA. This is a documented caveat of ss≈½ + nascent, an identifiability
floor of the strand signal, not a fit bug. (Stranded data resolves it via the strand tilt; the spliced anchor
helps but cannot fully close the unstranded gap.)

## 6. Problem C — the spline fit has artifacts

The plots (`scratchpad/ehat_fits_loglog.png`) show the fit alternating between smooth, **piecewise-linear**, and
**parabolic/step** segments — it does not track the data cloud as a clean monotone curve. Root causes, all in
the §2.2 machinery:
- **Monotone cumsum-exp reparameterization** → where consecutive increments are ~0 then one is large, the curve
  is **flat-then-jump** (the step-like bumps). The non-decreasing constraint turns flat-then-rising data into
  visible steps.
- **Cubic B-spline (degree 3)** → the **parabolic/cubic segments** between knots.
- **One global GCV-λ** → local over-smoothing (linear stretches) coexisting with local under-smoothing (bumps);
  a single λ cannot adapt to a curve that is flat at low z and steep mid-range.
- **Bisquare IRLS (2 iters)** → down-weighted points create **local flat patches** where the curve ignores the
  cloud.
- **No low-z anchor** (§4) → the left end floats on sparse/noisy low-z teachers.

These are cosmetic on strong-enrichment data (R²=0.93 hides them) but will be load-bearing on noisy **real**
data, and they undermine trust in the curve. This is the core of the model-fit improvement effort.

## 7. Solutions WITHIN the monotone-P-spline framework

1. **Anchor the low-z end** (§4): inject the depleted floor as a fixed left-end value — e.g. a high-weight
   synthetic anchor point `(z→z_min, floor)`, or a boundary constraint on the spline's left coefficient, or a
   prior pulling `α₀` toward `log floor`. Directly removes the float.
2. **Tame the artifacts:**
   - **Stronger / adaptive smoothing:** raise the λ floor, or use **adaptive λ** (more penalty at low z where
     data is sparse); or reduce `k` (fewer knots → fewer wiggle modes).
   - **A smoother monotone parameterization:** replace cumsum-exp with **I-splines / SCAM-style monotone
     smoothing** (monotone constraint without the step-prone cumulative-exp), or **PAVA isotonic + a smoothing
     pass**.
   - **Reconsider the bisquare IRLS** — on a genuine monotone signal it mostly creates flats; a gentler robust
     loss (Huber) or none may track the cloud better.
3. **Keep the framework, swap the loss to a shape-constrained GAM** (mgcv-style monotone P-spline with a proper
   smoothness prior and REML λ) — same family, fewer artifacts.

## 8. ALTERNATIVE fits (outside the free P-spline) — to test on the saved points

The smartest way to choose is empirical: fit each to the saved `(z, ρ_g, w)` datasets (sim now; **real VCaP
next**) and compare smoothness, anchoring, robustness, and held-out fit.

1. **Parametric saturating curve (recommended first try):** a 2–3-parameter monotone form, anchored by
   construction:
   `ê(z) = floor + (ρ_max − floor) · z^h / (K^h + z^h)`  (Hill), or Michaelis–Menten (`h=1`).
   - `floor` = the depleted anchor (§4, fixed, not fitted); `ρ_max`, `K`, `h` fit by weighted NLS.
   - **Smooth by construction** (no steps/parabolas), **monotone by construction**, **anchored by
     construction**, robust (few parameters), interpretable (`K` = the z at half-saturation). Collapses to flat
     when the data has no rise (`ρ_max→floor`).
2. **Monotone GP** with mean function = `floor`: shrinks to the floor where data is sparse, smooth, gives
   uncertainty — but heavier and less interpretable.
3. **Isotonic regression (PAVA) + a smoothing kernel:** non-parametric monotone, then smooth — cheap, but still
   needs anchoring + a smoothing bandwidth.
4. **Empirical-Bayes shrinkage to the floor:** a parametric prior centred on `floor`; the posterior mean is a
   smooth shrinkage curve. Ties the anchor and the smoothness into one Bayesian object.

**Experiment harness.** `scratchpad/enrichment_fit_points/{flagship_strong,unstranded_moderate,offcap_none,
phantom_zerogdna_nascent}.tsv` hold the columns `z, rho_g, fit_weight, recal_weight` + a header with `ρ_global`
and the depleted floor. Fit each candidate, overlay on the cloud, and score: (a) does it anchor to the floor at
low z? (b) is it artifact-free? (c) does it collapse to flat on `offcap`/`phantom`? (d) held-out residual. Then
re-run on the **real VCaP** points (extraction pending — needs a human index + a scan of the 8.3 GB
`annotated.bam`) — real data is the true test, since the artifacts matter most under real noise.

## 9. Summary of the open work

| problem | status | lever |
|---|---|---|
| A. wrong anchor (ρ_global, 42× too high under capture) | proposed | depleted floor = intergenic+intronic regions + pseudocount |
| B. nascent-RNA contamination of ρ_g (unstranded) | **accepted caveat** | none — documented ss≈½ identifiability limit |
| C. spline artifacts (steps/parabolas/float) | the target | anchor + smoother monotone fit, or a parametric saturating curve |

Recommended path: **(1)** extract the real VCaP fit points; **(2)** fit the §8 candidates (esp. the anchored
Hill curve) to sim + real points offline; **(3)** adopt the cleanest anchored-smooth-robust form, validated to
hold the strong-enrichment flagship and collapse correctly off-capture.
