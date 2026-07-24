# gDNA-hyperprior background-floor derivation

**Status:** derived + implemented + validated against the live ambig sim (2026-07-19). Produced by a 4-way
derivation workflow (independent principles → adversarial verification vs the oracle grounding → synthesis).
Implemented in `background_reference.measure_background` (`log_rho_floor`) + `npmle.DensityNPMLE.fit` (the
aggregate-cell injection). LIVE validation: the floor matches the derived predictions to the digit and lands
within 0.53 log of oracle (gdna1_verystrong, a fundamental per-region-concentration limit) vs the old 1/ΣE
seed at ~3 logs too low. The floor is currently INERT in production until the HYBRID background is wired into
the Phase-2 refit hyperprior fit (`calibrate._fit_gdna_hyperprior`).

---

Verified against the live code. Here is the reviewer-ready synthesis.

---

# FINAL: gDNA-hyperprior background-floor derivation

## Verdict on the candidate field

All four candidates agree on the diagnosis (the shipped floor resolves the pool as *one genome-sized region* → 1/ΣE ≈ −6.83, ~3 logs below the depleted mode = the confident-false-positive seed) and on two of the three acceptance properties (resolvable-exact, FP-safe). They differ only in **the E-aggregate for the zero regime** and **hard-max vs. soft-blend**. Scored against the oracle, that ordering is unambiguous:

| candidate | E-aggregate | combine | max\|err\| |
|---|---|---|---|
| poisson-em-limit | arithmetic 1/⟨E⟩ | max | 0.71 |
| resolution-limit | geometric 1/geomean | max | 0.64 |
| oracle-empirical | harmonic 1/harmmean | max | 0.69 |
| **bayesian-pooling** | **harmonic 1/harmmean** | **log-space precision blend** | **0.53** |

The winner is **harmonic-mean E-aggregate + precision-weighted log-space blend** (candidate `bayesian-pooling`). It is the only one with both a clean Fisher/resolution interpretation for the aggregate AND a principled reason to beat the hard-max on the one hard condition. I adopt it as THE formula.

---

## (1) THE floor-location formula

Let the pooled background regions (non-exonic: intergenic [+intron on real data]) have counts `g_i` over gDNA eff-lengths `E_i`. Define

- `G = Σg` (pooled counts), `E_tot = ΣE` (pooled exposure)
- `n0 = #{i : g_i = 0}` — the zero-count (unresolved) regions
- **`ρ_res = (1/n0)·Σ_{i:g_i=0}(1/E_i) = 1/harmmean({E_i : g_i=0})`** — the resolution wall

Then, in natural log,

```
ln ρ_floor = ( G·ln(G/E_tot) + n0·ln(ρ_res) ) / ( G + n0 )
```

with the two exact limits

- **n0 = 0** (every region resolved) ⇒ `ln ρ_floor = ln(G/E_tot)` — the pooled rate, EXACT.
- **G = 0** (true zero) ⇒ `ln ρ_floor = ln(ρ_res)` — the single-typical-region resolution wall.

`floor_log10 = ln ρ_floor / ln 10`. The E-aggregate is the **harmonic mean of the zero-count regions' effective lengths** (equivalently, the arithmetic mean of their per-region Poisson resolution caps `1/E_i`).

---

## (2) Why the resolution limit, and why harmonic

**Why the floor is a resolution limit, not a data mean.** `DensityNPMLE.fit` fits P(log ρ) in log-rate space; each background region enters as a Poisson cell whose likelihood is `e^{−ρE_i}` for a zero-count region. That barrier's detection knee — where the expected count first reaches 1 — sits at `ρ = 1/E_i`. A zero-count region therefore *cannot* separate any density below `1/E_i` from true zero; its Fisher information in log-rate, `ρ·E_i`, falls below 1 there. The shipped code pools all regions into one `(Σg, ΣE)` cell, so `e^{−ρΣE}` drives ρ to `1/ΣE` — it asserts the whole genome is a single region and claims a resolution `n×` finer than any region actually delivers. That is the ~3-log collapse. The honest population floor is the density at which a *typical* single region still reads ~zero.

**Why harmonic (= mean of per-region caps).** The per-region caps are `r_i = 1/E_i`. The population's typical resolvable floor is the plain average of those caps, `mean(1/E_i) = 1/harmmean(E)`. This is the mechanistically exact object — it is the mean of the quantities that actually carry the resolution, not an average of the E's themselves. The arithmetic-mean-E floor (`n/ΣE`) averages exposures before inverting and undershoots the true-zero condition by 0.7 log; the geometric mean is the log-space average of the thresholds and lands between. Only the harmonic mean is `mean(per-region resolution)`, and it also fits the anchor best (none −4.16 vs oracle −3.78). `n0` counts **only zero-count regions** — a region that measured a count has resolved its own density and does not vote for the floor; this is what deactivates the floor exactly at the resolvable condition.

**Why blend, not hard-max.** The estimator is a hierarchical Poisson: the pooled log-rate MLE carries Fisher information `G` (the total count), and the resolution floor is a prior at `ln ρ_res` carrying pseudo-information `n0` (one "my density is below `1/E_i`" vote per unresolved region). Precision-weighting the two in log-rate space is the exact posterior mode. As `G/n0 → ∞` it recovers the pool; at `G=0` it is the wall; at `n0=0` it is the pool exactly. The soft blend continues shrinking toward the pool as counts accrue, which is precisely what rescues the one bimodal condition the hard-max cannot (below).

---

## (3) Predicted floors vs oracle

| condition | oracle | **synthesis** | err |
|---|---|---|---|
| gdna_none (Σg=0) | −3.78 | **−4.16** | +0.38 |
| gdna1_verystrong | −4.84 | **−4.31** | +0.53 |
| gdna5_capON (n0=0) | −3.06 | **−3.06** | 0.00 (exact) |
| gdna5_verystrong | −4.32 | **−4.39** | −0.07 |

**max \|err\| = 0.53 log** (at gdna1_vs), RMS ≈ 0.33 log. Resolvable-exact at gdna5_capON to the digit. Versus the shipped 1/ΣE seed at ~3.0 logs low on every depleted condition, this removes essentially all of the catastrophic collapse.

---

## (4) Implementation plan (verified against current source)

Two surgical edits, no new tuned constant (harmonic mean is a mathematical aggregate; `G`, `n0` are counts; the blend weights are Fisher info).

**(A) `src/rigel/calibration/background_reference.py` — `measure_background`, after line 97** (`nr = int(pool.sum())`), before the `log_rho_bg` line:

```python
zc = pool & (counts <= _EPS)          # the unresolved (zero-count) regions
n0 = int(zc.sum())
s_recip = float((1.0 / eff[zc]).sum()) if n0 > 0 else 0.0
rho_res = (n0 / s_recip) if s_recip > _EPS else 0.0     # 1/harmmean(E_zc) — the resolution wall
if sg > _EPS and n0 > 0 and rho_res > _EPS:             # partially-depleted: precision-weighted blend
    L_pool = np.log(sg) - np.log(se)
    log_rho_floor = float((sg * L_pool + n0 * np.log(rho_res)) / (sg + n0))
elif sg > _EPS:                                         # n0==0 fully resolved ⇒ pooled rate EXACT
    log_rho_floor = float(np.log(sg) - np.log(se))
elif rho_res > _EPS:                                    # Σg==0 true zero ⇒ the wall
    log_rho_floor = float(np.log(rho_res))
else:
    log_rho_floor = -np.inf
```

Add `log_rho_floor: float` to `BackgroundReference` (and optionally `n0: int`, `rho_res: float` as diagnostics). **Leave `log_rho_bg`/`sigma_bg` unchanged** — they remain the dormant one-sided diagnostic; `log_rho_floor` is the new meaningful low-density location that replaces the too-low 1/ΣE seed downstream.

**(B) `src/rigel/calibration/npmle.py` — `DensityNPMLE.fit`:**

- **Line 204** (grid-bottom fallback), replace
  `lo = min(lo, float(np.log(max(background.n_counts, 1.0) / background.eff_total)) - 3.0 * h)`
  with
  `lo = min(lo, float(background.log_rho_floor) - 3.0 * h)`
  so the grid bottom sits at the resolution floor, not 1/ΣE.

- **Line 209** (aggregate-cell count injection), replace
  `gc = np.append(gc, float(background.n_counts))`
  with
  `gc = np.append(gc, float(np.exp(background.log_rho_floor) * background.eff_total))`
  leaving `ec = eff_total`, `wc = n_regions` as-is. This places the injected genome-scale Poisson cell at density exactly `ρ_floor` — its huge `ec = ΣE` makes it a sharp low mode there instead of `e^{−ρΣE}` driving ρ → 1/ΣE.

**Identity-preserving in the resolvable regime:** when `n0 = 0`, `log_rho_floor = ln(Σg/ΣE)`, so `exp(log_rho_floor)·ΣE = Σg` — the injected cell and the grid-lo fallback are **byte-identical to today**. Only the zero / near-zero conditions change. Gate behind the existing A/B and validate on the 4 grounding conditions before wiring `log_rho_floor` into `bp_solver._kde_logprior`. No change to `project()` or the enrichment role.

---

## (5) FP-safety argument

The KDE projection kernel bandwidth is `h = 0.15` decade (`h = bandwidth·ln10` at line 180). Across the depleted conditions the synthesis floor lands at −4.16…−4.39 — **at most 0.38 log below any oracle depleted mode, and at/above it on the near-zero conditions**. A genuinely depleted observation therefore sits within ≈2.5 bandwidths of the injected low mode ⇒ it is captured by that low mode, not drawn up to an enriched mode (which is ≥7 bandwidths away). This is the exact property the fix buys over the shipped 1/ΣE = −6.83 seed, which is ~3 logs / ~20 bandwidths below the depleted density: with no low mode anywhere near the real depleted rate, the shipped seed leaves depleted observations to be pinned onto whatever enriched mode exists — the confident-false-positive cliff that crushes strand-blind gDNA. Additionally the synthesis errs slightly **high** (over-floors) at the bimodal condition, which is conservative-safe: the hyperprior never claims a density below one resolvable count per typical region.

---

## Residual gaps (flagged)

1. **Strict 0.5-log gate: 0.53 at gdna1_verystrong — marginally over.** This is a genuine bimodal-unmixing model mismatch, not a fixable rate error: the 31 counts concentrate into ~13 count-bearing regions that form a *separate enriched mode*, freeing the 203 zero-count regions to collectively resolve *below* the single-region wall (oracle −4.84 < wall −4.15). No monotone functional of the summary stats `(Σg, ΣE, {E_i})` can see the per-region count *concentration*; only a full per-region NPMLE — what the oracle effectively is — captures it. The soft blend narrows this from the hard-max's 0.69 to 0.53 by continuing to shrink toward the pool as Σg grows, but cannot close it without a per-region count histogram, which the aggregate-cell design deliberately does not carry.

2. **Downstream-negligible.** The floor is the weakest arm of ψ (n_eff ≈ 0.15 pseudo-obs). A 0.5-log floor error nudges log f_g by <0.1 of a count; the node's composition is dominated by the strand likelihood + imputation messages. The fix's value is removing the 3-log cliff, not the last 0.1 log.

3. **The +0.38 residual at true-zero (none)** is the NPMLE kernel/grid mode-lift, not a rate error — the smoothed grid readout sits a touch above the injected cell's rate.

4. **Substrate/uniformity caveat.** Predictions assume the sim substrate. On real data (`include_introns=True`) and under tumor CNV the DNA-uniformity assumption weakens (per the module docstring §"tumor CNV breaks it"); revalidate `ρ_res` and the blend on real data before trusting the absolute floor.

5. **Predictions are the candidates' self-reported values** — per-region `E_i` were not in hand to recompute the blends independently. The two structurally-pinned conditions are exact by construction (none = pure `ρ_res` geometry at Σg=0; gdna5_capON = pooled MLE at n0=0); only gdna1_vs and gdna5_vs actually exercise the blend, and one of those two is the >0.5 case. Recommend recomputing all four from the live `measure_background` under the A/B before sign-off.

**Files:** `/Users/mkiyer/proj/rigel/src/rigel/calibration/background_reference.py` (measure_background, lines 95–106), `/Users/mkiyer/proj/rigel/src/rigel/calibration/npmle.py` (DensityNPMLE.fit, lines 197–212).