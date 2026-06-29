# Dispersion-aware message precision — model-fitting state + design

**Status:** design, not implemented. Authored 2026-06-29 after a multi-scenario + real-VCaP investigation
of the calibration message-reliability model fit. Decision (with the user): **the cross-node *disagreement*,
not the mean density level, should be the message-precision estimator.** This doc is the self-contained
hand-off: the current state, the evidence, the alternatives weighed, and the recommended design with enough
detail to implement in a new session.

Constraints that bind every choice (project rules): **no capture-detection gate** (brittle on real data);
**no magic numbers / new tuned constants** (≤ ~8 total; pause before adding one); **messages communicate
densities, not fractions**; principled/derivable, not band-aids over compensating errors.

---

## 0. How to resume (orientation)

- **Solver:** the calibration is a forward-backward belief-propagation sweep over a bipartite region↔boundary
  node chain. Per-node solve = the log-density log-odds solver (`simplex_logodds`; the 2-simplex lattice is
  retired). Each node imputes its gDNA fraction `f_g` from three sources: its **strand likelihood**
  (Beta-Binomial tilt — the only intrinsic per-count signal), the **cross-node messages** (this doc), and the
  **global gDNA prior** (population baseline + enrichment transfer). Architecture: `CALIBRATION_ARCHITECTURE.md`.
- **The message machinery lives in** `src/rigel/calibration/bp_solver.py`:
  - `node_sweep` — the outer var~mean fixed point + the inner FB sweep.
  - `_scan` (inside `node_sweep`) — the sequential relay that builds + applies each message. **This is where
    the redesign lands.**
  - `_edge_varmean` / `fit_gdna_varmean` / `fit_rna_varmean` — fit the current `σ²_bio(μ)` curve. `_edge_sigma2`
    (inside `node_sweep`) queries it per edge. **These retire under the redesign (see §5.4).**
  - `_global_logprior`, `_gdna_seed_estimate` / `_fit_seed_varmean`, `fit_enrichment_transfer` — the GLOBAL
    prior + enrichment transfer ê(z). **These stay** (separate from the message precision).
- **The fit code:** `src/rigel/calibration/variance_model.py` — `MonotoneVarMean` (the monotone P-spline),
  `MonotoneMean` (the ê transfer), `BlendedVarMean` (the mean/median fork — retires).
- **Reproduce the diagnostics** (synthetic scenarios in `~/Downloads/rigel_runs/quick_3to1_5mb/`):
  - `scripts/debug/dispersion_study/edge_hump_diagnostic.py` — extract the gDNA edge population (μ, raw,
    off-floor) on a scenario; prints the per-decile hump + monotone under-fit.
  - `scripts/debug/dispersion_study/causal_free_smoother.py <cond>` — the free-vs-monotone causal test.
  - `scripts/debug/dispersion_study/plot_synthetic_fits.py` — the 3 synthetic figures (σ²_bio fit + LOESS, ê(z),
    dispersion-vs-disagreement). `plot_vcap_fits.py` — the VCaP versions (needs a fresh cached scan).
  - All write PNGs to the session scratchpad — **edit the `SCR`/`OUT` path before running.** Activate the
    `rigel` conda env first; prefix with `OMP_NUM_THREADS=4`.

---

## 1. The model-fitting state (what exists today)

The message precision a node assigns to a neighbour's gDNA message is, in `_scan` (log-density form):

```
τ  =  1 / ( vbg[src] + σ²_bio(μ) + pois )
```

- `vbg[src]` — the source node's running-belief log-variance `Var(log f_g)` (the relay carries it).
- `pois = 1/count_src` — the source's log-density sampling variance (non-detection ⇒ count→0 ⇒ pois→∞ ⇒ ~0
  precision; "zero density is not a measurement").
- **`σ²_bio(μ)`** — the biological/process dispersion, the subject of this doc. A **monotone-increasing
  P-spline in `log μ`**, `μ = ½(ρ_dst + ρ_src)` (the mean of the two endpoints' gDNA densities), fit each outer
  pass on the response `raw = (Δlog ρ_g)²` minus a Poisson offset (`MonotoneVarMean.fit_offset`,
  `fit_gdna_varmean`). The message mode `mo` is the source density re-expressed in the destination's log-f_g
  frame; the combine is a Gaussian product in log space.

Two companion fits (NOT message precision — they parameterize the GLOBAL prior; they **stay**):
- **ê(z)** (`MonotoneMean` via `fit_enrichment_transfer`) — the enrichment transfer `E[ρ_g | z]`, z = the
  boundary-crossing gDNA flux density; sets the enrichment-aware global prior on AMBIG exon nodes.
- **σ²_g / var_mean** (`_gdna_seed_estimate` → `_fit_seed_varmean`) — the population gDNA spread, sets the
  global prior's precision.

Ancillary, entangled with σ²_bio:
- **mean/median fork** (`population_spread` flag + `BlendedVarMean`) — `fit_rna_varmean` blends a
  conditional-MEAN fit and a reliability-weighted-MEDIAN fit by the enrichment weight `w`. **This is a band-aid
  for the monotone mis-fit (§2) and retires with it.**
- **Poisson offset** — `fit_offset` models `E[raw] = σ²_bio(μ) + V_p`, `V_p = 1/(ρ·E+1) + …` (the +1-pseudocount
  log-Poisson floor). Minor: under-counts sub-unit-count sampling → a small spurious low-μ σ²_bio floor.

---

## 2. The diagnosis (why the current fit is wrong)

### 2.1 The dispersion is HUMP-shaped in μ, but the model is monotone-increasing

On the controlled flagship (gdna300, ss0.99, capture-ON), the off-floor, converged-belief gDNA-edge
disagreement `(Δlog ρ_g)²` per μ-decile:

```
0.2  0.3  17.4  23.8  29.7  33.3  2.6  1.2  1.4  3.2
```

A clean **hump**: ~0 where both nodes are depleted (agree), peaking ≈ 33 in the **depleted↔enriched transition**
(intermediate μ), falling back to ~2 where both are enriched (agree — enriched nodes are similar to each other).
The monotone-increasing σ²_bio can only rise, so it pins at a low plateau and **under-fits the transition ~12×**
(message precision there ~10× too high → confident-wrong messages drag enriched exons down = the capture
under-call) while **over-fitting the high-μ tail**. Verified off-floor (not a 1/E clamp artifact) and on the
converged belief, not just init.

### 2.2 It IS the capture signature (the control)

cap-OFF (uniform gDNA, no enrichment regimes): **no hump** — ρ_g unimodal, the raw-vs-μ field flat (~20× smaller
dynamic range); the monotone fit is appropriate there (even over-fits). No regimes → no transition → no hump.
**So the hump is the capture/enrichment signature, not a generic artifact.** unstranded ss0.50 cap-ON reproduces
the hump (milder, plateau set lower by κ≈½).

### 2.3 The deeper truth: dispersion is the DISAGREEMENT, not the level — and this transfers to real data

The decisive finding (`plot_synthetic_fits.py` fig3, and **confirmed on real VCaP** `plot_vcap_fits.py` fig3):
the dispersion `raw` is a **tight function of the cross-node disagreement `|Δlog ρ_g|`** and only an *incidental*,
messy function of μ. At a fixed μ, agreeing and disagreeing edges coexist; the whole vertical scatter is
explained by `|Δlog ρ|`. **μ is the wrong predictor axis.** This is why *any* 1-D μ-smoother (monotone, free,
LOESS) is mis-specified: it can fit the aggregate conditional-mean hump but can never separate a genuinely-
agreeing edge from a seam-straddling one at the same μ. **This finding holds on BOTH synthetic and real VCaP —
the only piece of the diagnosis that transfers cleanly to real data** (see §6 for what does NOT).

### 2.4 A second, INDEPENDENT error (not addressed here): ê(z) predictor compression

Separate capture under-call, in the global prior not the message. The predictor z (boundary-crossing flux)
**compresses ~5× at high enrichment** (high-E median z≈6 vs contained interior density ≈28 vs true E≈33), so the
monotone log-log ê plateaus → over-predicts transition exons (E∈[1,10) by 2.86×) and under-predicts the most-
enriched (E≥100 by 4.2×). The fix is a richer PREDICTOR (the contained interior density, which tracks truth at
corr 0.99, fit on single-strand exons as a non-circular teacher), **not** a different curve form. **Tracked as a
follow-up; orthogonal to this doc.** (ê on VCaP is noise-dominated/flat — can't be validated there yet, §6.)

---

## 3. Alternatives weighed (and why disagreement-aware wins)

| Approach | Fits the hump? | Fixes the conflation (§2.3)? | Constraints | Verdict |
|---|---|---|---|---|
| **Monotone σ²_bio(μ)** (current) | **No** (structural) | No | ok | The error. |
| **Free / shape-free σ²_bio(μ)** spline | Yes (causally verified) | **No** (still μ-axis) | ok, small | Real partial fix — but wrong axis; tail/sparse-region instability (overshoots the thin high-μ tail). |
| **LOESS σ²_bio(μ)** (the historical impl) | Yes | **No** | reintroduces a span constant; O(n²) robust pass | Same class as free spline, no advantage; rejected. |
| **2-predictor σ²_bio(μ, \|Δlogρ\|)** | Yes | Yes | risks a detector if 2nd pred is a regime label | Collapses to disagreement-aware when the 2nd predictor is the disagreement, plus a redundant μ axis. |
| **Disagreement-aware precision** | n/a (no μ curve) | **Yes** | ok, no detector, no new constant | **Recommended.** |

Verified facts behind the table:
- The **free spline** is a genuine improvement over monotone — re-measured penalized RSS **39,288 (free) vs
  264,135 (monotone)** on the flagship edges (free ⊇ monotone feasible set ⇒ must be ≤, the earlier workflow's
  "free worse" was an implementation artifact); causal test (`causal_free_smoother.py`): drop-the-constraint
  fits the hump (σ²_bio @transition 2.9→32) and materially helps where messages dominate (unstranded library
  gDNA frac 0.388→0.459) with **no harm off-capture** (0.756→0.753, self-collapses). So a free spline is a valid
  *fallback*. But fig3 shows it fits the wrong axis.
- **LOESS** (historical `_loess`, local-linear tricube + bisquare, span 0.4, from commit 200143d4) tracks the
  hump like the free spline but is the same *class* (a μ-smoother), reintroduces the un-CV'd span constant, and
  is O(n_query·n) per eval (the sweep queries σ²_bio many times). No advantage over the free GCV spline.

**Why disagreement-aware:** it is the only candidate that addresses §2.3 (the real, transfer-able structural
truth) — and it does so with no μ curve at all, no regime detector, and no new tuned constant.

---

## 4. The recommended design — disagreement-aware (robust) message precision

### 4.1 Principle

A message src→dst proposes "dst's gDNA density ≈ ρ_src." Its reliability is the **edge process variance**
`σ²_edge` — how much true gDNA density changes across this edge. `σ²_edge` is small within a regime (smooth
imputation, trust it) and large across a regime boundary (don't impute across it). We **estimate `σ²_edge`
per-edge from the message's *surprise*** — how far it falls from the destination's own *independent* evidence,
in excess of the sampling + belief uncertainty that already explains some scatter. A surprising message (one
demanding a big jump the destination's own evidence contradicts) reveals a high-variance edge and self-silences.
This is robust/heteroskedastic Gaussian belief propagation — the edge variance is a hierarchical latent we infer
from the residual, not a curve we fit against μ.

### 4.2 The non-circular anchor (critical)

The surprise must be measured against the destination's **message-free local belief** — its strand likelihood +
global prior, with NO neighbour messages — NOT its evolving running belief. (If anchored on the running belief,
a destination already dragged toward a wrong source shrinks its own surprise → trusts more → runaway.) That
quantity already exists in `_scan`: `fg_loc` (the phase-A message-free local solve), carried as `lfg_loc` (its
log-fraction) and `pg_loc`/`vg_loc` (its precision/variance). It is fixed across the scan and recomputed each
outer pass from the node's own evidence — independent of the messages. ✔ non-circular.

### 4.3 The formula (gDNA message; per strand for RNA, symmetric)

In `_scan`, for a message from `src` into `dst` (all log-fraction quantities are already in the destination's
log-f_g frame):

```
mo        = log(ρ_src in dst frame)                  # the message mode — UNCHANGED from today
resid     = mo − lfg_loc[dst]                        # surprise vs the dst's INDEPENDENT (message-free) belief
base_var  = vbg[src] + pois                           # source posterior log-var + sampling (= 1/pr0); UNCHANGED
expected  = base_var + vg_loc[dst]                    # variance that already explains scatter (sampling+belief)
σ²_edge   = max( resid**2 − expected , 0.0 )           # per-edge process variance (method-of-moments, floored)
τ         = 1 / ( base_var + σ²_edge )                 # the message precision  ← replaces  1/(base_var + σ²_bio(μ))
# combine into the running belief (unchanged): fbg[dst] = exp((pg_loc·lfg_loc + τ·mo)/(pg_loc+τ)); vbg[dst]=1/(pg_loc+τ)
```

The single line that changes is the precision denominator: `σ²_bio(μ)` (a fitted curve) → `σ²_edge` (a per-edge
quantity from the surprise). Everything else in `_scan` is unchanged. `σ²_g`/`var_mean` (global prior) and ê(z)
are untouched.

### 4.4 Why this is correct in all three regimes (the behaviours we need)

- **Agreeing edge** (within-regime, smooth): `resid ≈ 0` ⇒ `σ²_edge = 0` ⇒ `τ = 1/base_var` (full precision) ⇒
  the smooth imputation is trusted, and the long-range relay flows. ✔
- **Conflicting message into a CONFIDENT (single-strand) destination** (`vg_loc` small): `expected` small,
  `resid²` large ⇒ `σ²_edge ≈ resid²` ⇒ `τ` small ⇒ the depleted-seam message **cannot drag the enriched exon
  down**. ✔ (the capture under-call mechanism, silenced.)
- **Conflicting message into an UNCERTAIN (AMBIG) destination** (`vg_loc` large — its local belief is just the
  global): `expected` large ⇒ `σ²_edge = max(resid² − large, 0) ≈ 0` ⇒ `τ ≈ 1/base_var` ⇒ a correct enriched
  neighbour **still lifts the AMBIG exon to its true level**. ✔ **This is the property a naive "down-weight by
  raw disagreement" would break, and why the anchor's *precision* (`vg_loc`) must enter via `expected`.**

So it self-silences exactly the harmful seam messages while preserving both the smooth relay and the AMBIG lift —
the three things that must all hold — with no regime label anywhere.

### 4.5 What this is, statistically

`σ²_edge = max(resid² − expected, 0)` is the per-edge method-of-moments estimate of the edge's excess (process)
variance, given the known sampling + belief variances. The resulting `τ` is a robust message precision: messages
are down-weighted in proportion to their surprise relative to the destination's confidence — the standard form of
robust Gaussian BP / a Student-t-like message. The data (fig3: `raw ≈ |Δlogρ|²`) is exactly the statement that
the realized residual² IS the dispersion, so estimating `σ²_edge` from it is using the right observable.

---

## 5. Scope of the change

### 5.1 Lands in `_scan` (and the symmetric RNA branches)
Three message blocks in `_scan` (`emit_g`, `emit_p`, `emit_n`): swap the `s2{g,p,n}[i]` term for the per-edge
`σ²_edge` computed from `(mo − lfg_loc[i])`, `vbg[lsrc]`, `pois`, and the destination's local variance. For RNA,
the anchor is the destination's local *per-strand* log-fraction (`lfp_loc`/`lfn_loc`) and `vp_loc`/`vn_loc`.

### 5.2 The destination local variance must be available in `_scan`
`vg_loc` / `vp_loc` / `vn_loc` are computed in phase A (`_local_solve`) and are present as `pg_loc`/`pp_loc`/
`pn_loc` (their inverses). Carry the variances (or use `1/pg_loc`) into the `_scan` closure.

### 5.3 Spliced (RNA) floor unchanged
The one-sided spliced/mature floor is a node-local term (D8), not a message — leave it.

### 5.4 RETIRES (net simplification — a whole layer goes)
- `_edge_sigma2` (the per-edge σ²_bio query), `_edge_varmean`, `fit_gdna_varmean`, `fit_rna_varmean`.
- The mean/median fork: the `population_spread` *path inside `fit_rna_varmean`* and `BlendedVarMean` (it only
  existed to lift the monotone curve toward the hump under capture — a band-aid for §2.1; measured inert
  off-capture, diverging 4.6× under capture purely for this).
- The per-pass refit of the message σ²_bio in `node_sweep`'s outer loop.

### 5.5 STAYS (do NOT remove)
- `_gdna_seed_estimate` / `_fit_seed_varmean` and `var_mean`, `σ²_g` — these parameterize the **global prior**
  precision (`_global_logprior`), a per-node term, not an edge message.
- `fit_enrichment_transfer` / ê(z) / its `σ²_level` (which *does* use `population_spread=True` — keep that path
  for the global prior even though the message fork retires).
- `MonotoneVarMean` and `MonotoneMean` classes (still used by the global-prior + ê fits).
- The source running-belief relay (`fbg[src]`) — the message mode `mo` still uses it (the relay is correct).

---

## 6. Open questions / risks (resolve during implementation)

1. **Single-residual noise.** `σ²_edge` from one residual is χ²₁-noisy. Acceptable for a one-sided protective
   down-weight, but a message could be unlucky (large residual by chance) and be over-silenced, or lucky and
   under-silenced. Options to consider: a mild shrinkage of `σ²_edge` toward 0; a light pooling that does NOT
   reintroduce a μ axis; or accept the noise (it averages out over a node's several edges in the FB combine).
   **Decide empirically with the kill-switch experiment (§7).**
2. **Over-silencing the relay (the kill-switch).** Must confirm the new τ does NOT cut genuine long-range relay
   edges (large μ-separation but small disagreement-from-anchor). §4.4 argues it doesn't (agreeing ⇒ full τ),
   but verify on the converged belief before committing.
3. **`expected`/combine double-use of `vg_loc`.** `vg_loc` enters both `expected` (the surprise threshold) and
   the final Gaussian combine (`pg_loc`). Check the derivation is self-consistent (it is the dst's belief
   variance in both roles, which is coherent, but confirm no double-counting of the dst evidence).
4. **Offset cleanup (independent, cheap).** Swap the +1-pseudocount log-Poisson offset for the delta-exact
   `1/max(count, count_floor)` (reuse the existing `_MSG_PSEUDOCOUNT/E` floor) — removes a spurious low-μ σ²_bio
   floor and reconciles the fit-time vs apply-time floor conventions. Second-order; bundle or defer.
5. **ê(z) compression (§2.4)** — separate follow-up; not in this change.
6. **Goldens WILL change** — the synthetic regression goldens must be regenerated after this lands (an intended
   solver change). Validate FIRST on the net-flow benchmark (§7).
7. **VCaP can't validate this yet (§ below).**

### Real-data (VCaP) caveats — deferred per the user (generating new real data)
On the cached VCaP slice the redesign's premise can't be *tested*: `eff_gdna` spans ~8 orders of magnitude so
~99.5% of edges floor-compress into a low-μ pile (μ tops at ~15 — the enriched downslope is never reached), and
κ fits at **0.00086**. The hump (§2.1) and ê compression (§2.4) therefore do NOT reproduce on VCaP — not a
contradiction, a coverage failure. **The κ=0.00086 question is unresolved and matters:** by the BB Fisher info
`(2κ−1)²≈1`, κ≈0 is a *fully reverse-stranded* library with near-maximal strand discriminability (fine) — but
the workflow also reported the raw sense fraction ≈0.44–0.49 (≈unstranded); those can't both be the same
quantity, so either VCaP is genuinely reverse-stranded or the κ fit collapsed on this head-of-BAM slice (then the
densities are suspect). Resolve κ before trusting any VCaP density. The ONE thing VCaP DID confirm: §2.3 (fig3) —
dispersion is the disagreement on real data too. Development stays on the controlled synthetic scenarios until
the new real data lands.

---

## 7. Validation plan (in order)

1. **Cheap pre-commit experiment — counterfactual silencing (no solver change, no golden regen).** On the
   converged flagship belief (reuse `edge_hump_diagnostic.py`), compute the proposed τ alongside the current τ
   for every gDNA edge and check: (i) transition/seam edges (high `|resid|` vs `lfg_loc`) get ≳10× *lower* trust
   (over-trust 9.5×→~1×); (ii) within-regime agreeing edges keep ~full trust (the relay survives — the
   kill-switch); (iii) AMBIG-destination edges with a correct enriched message keep high trust (the lift
   survives). If (ii)/(iii) fail, fall back to the **free shape-free spline** (§3) which carries no runaway risk.
2. **Implement** §4–§5 in `_scan`; retire §5.4.
3. **Benchmark:** the 16-condition net-flow deconvolution suite (skill `calibration-benchmark`,
   `~/Downloads/rigel_runs/gdna_benchmark_5mb`) — expect the unstranded/capture conditions (which regressed
   worst) to improve, off-capture + stranded to stay ~flat. Then regenerate the synthetic goldens
   (`pytest tests/ --update-golden`), confirm `pytest tests/` green.
4. **New real data** (user is generating) — re-run the §6 diagnostics once a sample with usable regime coverage
   + a settled κ exists; this is the real-capture validation the synthetic can't provide.

---

## 8. Evidence index (numbers cited above)

- Hump per-decile (flagship, off-floor, converged): `0.2 0.3 17.4 23.8 29.7 33.3 2.6 1.2 1.4 3.2`; monotone
  σ²_bio pinned ≈ 2–3 (production median fit) → ~12× under-fit at the peak.
- Free-vs-monotone penalized RSS on the flagship edges: **39,288 vs 264,135** (free fits the hump up-and-down:
  transition 30→32.6, high-μ 2.3→1.8; monotone flat-pinned at 15.3).
- Causal free-smoother (drop-the-constraint): flagship enriched f_g 0.706→0.719 (strand-dominated, small);
  **unstranded** median enriched f_g 0.376→0.458, library gDNA frac 0.388→0.459; cap-OFF 0.756→0.753 (no harm).
- ê(z): high-E z≈6 vs contained ≈28 vs true ≈33 (z compresses ~5×); E≥100 under-called 4.2×, E∈[1,10)
  over-called 2.86×; corr(log contained, log E)=0.99.
- VCaP: 53,494 gDNA edges, μ∈[1.2e-6, 15], κ=0.00086, ~99.5% at-floor.

Related memory: `calib_sigma2_hump_monotone_misfit.md` (the diagnosis + the free-spline correction),
`calib_logdensity_overhaul.md`, `hybrid_capture_derivation.md` (whose "bimodal seam" framing this sharpens to
hump-vs-monotone, then to disagreement-vs-level), `flagship_em_bound_not_calibration.md` (why the flagship's
own effect is small — it's EM-bound; the message fix matters most in the message-dominated regimes).
