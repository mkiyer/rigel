# gDNA-hyperprior projection derivation (the snap-to-depleted fix)

**Status:** derivation-workflow output (2026-07-19). Framing (the correct posterior / yardstick) + lever
analysis (bandwidth LOAD-BEARING, enrichment-skeleton CONDITIONAL, tempering disfavored, MIP rejected) +
FP-bounding principle + ranked prototyping plan.

> **⚠ Scope update (2026-07-19, post-A/B).** A **fitted RNA arm** (a learned `logP_r` / `π_r`) was
> prototyped and **rejected** — a full-solve A/B showed both a raw and an FP-bounded form fail (no recovery
> under message anchoring; catastrophic false positives). It is **not** a calibration option and has been
> removed from the roadmap; the yardstick below still *names* `π_r` as the (flat, reference) RNA prior in the
> generative model, which is correct and unchanged. The **live** levers are **bandwidth (S1/S2)** and the
> **enrichment-skeleton / occupancy-leveling (S3)**; see `CALIBRATION_MASTER.md` for the current plan.

## Framing — the correct generative model / posterior

I've read both files. Here is the statistical framing you asked for — the yardstick, before any lever.

---

# The correct generative model vs. the shipped ψ

## 1. Node observables and latents

A single node carries: unspliced mass `M` (≈ total count), effective length `E`, observed total density `ρ_tot^obs = M/E`, and the per-strand split `(u₊, u₋)`, `N = u₊+u₋`. The latents are two **absolute rates**: the gDNA rate `ρ_g` and the RNA (nascent-unspliced) rate `ρ_r`. Composition is the derived ratio `f_g = ρ_g/(ρ_g+ρ_r)`; total scale is `ρ = ρ_g+ρ_r`. The solver grids the log-odds `λ = logit f_g = log ρ_g − log ρ_r` (simplex_logodds §docstring).

The natural coordinates are `a = log ρ_g`, `b = log ρ_r`, so `λ = a − b` is one axis and the total scale is the other.

## 2. The correct generative model

Three independent facts generate the node:

1. **gDNA population** (the hyperprior). `ρ_g` is drawn from the fitted population density `π_g(a)` — the bimodal NPMLE (`DensityNPMLE.logP`): a tall depleted mode at `a≈−3.06·ln10` plus short enriched modes. **This is a prior over the absolute gDNA rate, not over composition.**
2. **RNA population.** `ρ_r` is drawn from `π_r(b)` — currently the flat Jeffreys reference (`_rna_arm` → `½·log(1−f_g)`), i.e. `π_r ≡ const` in `b`. gDNA-enrichment and RNA-expression are separate processes ⇒ the prior is the **product** `π_g(a)·π_r(b)`.
3. **The observation.** The total count is `M | ρ ~ Poisson(ρ·E)`, and the strand split is `(u₊,u₋) | f_g ~` the two-component Beta-Binomial (`strand_loglik`, `_mixture_strand_loglik`). The BB is informative about `f_g` **only** when `κ≠½` and `N` large; at `κ=½` (unstranded) it is exactly flat in `f_g` (Fisher info `∝ N(½−κ)² = 0`).

Joint posterior in `(a,b)`:

```
P(a,b | M,u) ∝ Poisson(M; (eᵃ+eᵇ)E) · BB(u | σ(a−b)) · π_g(a) · π_r(b)
```

The log-rate/logit Jacobian cancels exactly once per group (the module already relies on this — simplex_logodds §2), so no Jacobian is written; that part is correct.

## 3. The proper marginal over composition — where the counting variance enters

Reparametrize along fixed log-odds `λ`: put `b = a−λ`, integrate the total direction `a`. Since `f_g = σ(λ)`, the total is `ρ = eᵃ+e^{a−λ} = eᵃ/f_g`, so the Poisson total-count factor is `Poisson(M; (eᵃ/f_g)E)`. This factor, as a function of `a`, is a bump centered at

```
a*(λ) = log(f_g · M/E)          ← the "ray" the shipped code pins
```

with **width set by the node's own count**: Poisson Fisher information in log-total is `M`, so to second order

```
Poisson(M; (eᵃ/f_g)E) ≈ const · exp(−½ · M · (a − a*(λ))²)
```

— a Gaussian in `a` of **precision M** (SD `1/√M`). This is *exactly* where the node's own counting variance lives, and it is a **scale** likelihood on the total, fully consistent with count-zero-information on composition. The proper composition posterior is therefore

```
P(λ | M,u) ∝ BB(u | σ(λ)) · ∫ da  exp(−½M(a−a*(λ))²) · π_g(a) · π_r(a−λ)
```

Read structurally: **the correct arm is `π_g` convolved along the log-rate axis by the node's counting kernel (width 1/√M), then read at the ray, jointly weighted by a real `π_r` on the RNA the composition implies.**

## 4. The shipped ψ is the zero-counting-variance, flat-RNA degenerate limit

The production integrand is (simplex_logodds:267, `_gdna_arm`:155-163, npmle `logprior`:250-282):

```
ψ(λ) = BB(u|σ(λ)) + logπ_g(log(f_g·M/E)) + ½·log(1−f_g)
```

Match it to §3 by making **two** approximations:

- **(i) M → ∞.** The counting Gaussian `exp(−½M(a−a*)²)` collapses to `δ(a − a*(λ))`, so the integral becomes `π_g(a*(λ))·π_r(a*−λ)`, i.e. `logπ_g(log(f_g·M/E))`. The node's total is treated as **known exactly**, counting variance set to **zero** — the opposite extreme from noisy.
- **(ii) π_r flat.** The RNA arm degenerates to the improper reference `½·log(1−f_g)`, which bounds `f_g→1` weakly and does **nothing** at `f_g→0`. "Explain the whole node as RNA" is never penalized.

Under (i)+(ii) with the strand term silent (unstranded), `ψ(λ)` is just the height of `π_g` read along the ray `{ρ_g = f_g·M/E : f_g∈[0,1]}`. Its argmax is **the tallest mode with `μ_k ≤ log(M/E)`** — the tallest *sub-total* mode. For a node whose observed total sits at a **short enriched** mode, the **tall depleted** mode is a sub-total mode and out-votes it by the height ratio (your node 2589: peak at `f_g=0.18` hitting the depleted `−2.97`; every enriched `f_g` is ~4 logP lower ⇒ `e⁴≈50×`). The posterior median snaps `f_g` low. This is the "gDNA flows to RNA" failure, exactly.

So the shipped arm is not a wrong formula — it is the **correct posterior evaluated in the `M→∞`, `π_r≡const` corner**, which converts a prior over the *absolute rate* `ρ_g` into an *induced prior over composition* by slicing the bimodal density at a pinned total. That induced prior's mode is always the tallest reachable `π_g` peak.

## 5. Why `floor_eps=0.02` cannot fix it

The uniform interior floor (npmle:229-239) fills `−∞` valleys **between** modes so `np.interp` never returns `−inf`. But the failure is not a valley — it is the **relative height** of two present peaks (depleted vs enriched, ~50:1). Adding 2% of the mass as a flat pedestal shifts both peaks negligibly and leaves the height ratio essentially untouched. The floor fixes **support**; the bug is **relative mode mass**. Different axis.

## 6. The core reframe (and the count-zero-info tension)

The hyperprior is a prior on `ρ_g` (absolute). The shipped code turns it into a prior on `f_g` by fixing `ρ_tot = M/E` and rescaling. That conversion is legitimate **only** if you simultaneously (a) carry the `π_r` pushforward and (b) integrate the total with its Poisson counting kernel — neither of which the code does. Because it pins the total using `M/E`, the **total count leaks into composition** through the projection onto the bimodal `π_g`: a node's `M/E` decides which mode its ray lands nearest, and that manufactures a composition vote out of a count. This is precisely the count-zero-information violation flagged in your `npmle_is_enrichment_not_dna_prior` memory — density height, at a count-pinned total, becomes a composition vote.

## 7. The yardstick for the four levers

Every lever should be judged as an approximation to the **proper marginal**

```
P(λ | M,u) ∝ BB(u|σ(λ)) · ∫ da exp(−½M(a−a*(λ))²) · π_g(a) · π_r(a−λ)
```

against the **degenerate baseline** `ψ = BB + logπ_g(log(f_g M/E)) + ½log(1−f_g)`. The proper posterior exposes exactly four surfaces a lever can act on, and the framing tells you what each *is*:

- **Smoothing width of `π_g`** — the counting convolution `exp(−½M(a−a*)²)`. A fixed bandwidth `h` (npmle `fit`, `h=0.15`) is a **node-independent stand-in for the node-dependent counting kernel `1/√M`**: both smooth the same log-rate axis, and both lower a tall *sharp* depleted peak more than a broad enriched hump, shrinking the 50:1 ratio. This is why `h` is the dominant lever in this framing — it is the only place the missing counting/estimation width currently lives — and it points at an adaptive `h` that folds in `1/√M` (and the hyperprior's own estimation SD) rather than a magic constant.
- **The mode height ratio itself** — the `w_k`/peak structure of `π_g` on the composition ray.
- **The RNA arm `π_r`** — in the generative model this is a real prior on the RNA rate; in the code it is the flat Jeffreys reference `½·log(1−f_g)`. ~~A fitted `logP_r` would penalize the all-RNA escape~~ — **fitting a learned `π_r` was TESTED and REJECTED (2026-07-19, full-solve A/B):** under message anchoring it does not recover, and the FP-bounded (monotone) form manufactures gDNA everywhere. The RNA constraint is carried by the reference + the messages; there is no fitted RNA arm.
- **The total marginalization / counting variance** — restoring the `∫da` with precision `M` instead of the `δ`-pin.

Sensitivity/specificity reading of the yardstick: the `gdna_none` condition has `M/E` sitting *at* the depleted mode itself, so a*(λ) for any `f_g>0` climbs *up* the density into thinner mass — the proper posterior keeps `f_g≈0` there regardless of smoothing, because there is no enriched mode below the total to snap to. The 35 partially-enriched nodes have a genuine enriched mode below their total; smoothing + a real RNA constraint let that mode survive the depleted Everest. That asymmetry — enriched mode present-below-total vs. absent — is the structural reason the tradeoff is achievable rather than a wash.

**Files:** `/Users/mkiyer/proj/rigel/src/rigel/calibration/npmle.py` (`logprior` 250-282 = the ray interp at `log(f_g·M/E)`; `fit` 148-248 = `h`, `floor_eps`, aggregate-background cell; `project` 284-315), `/Users/mkiyer/proj/rigel/src/rigel/calibration/simplex_logodds.py` (`_gdna_arm` 155-163, `_rna_arm` 166-175, ψ assembly 267, posterior via `_posterior_median_fg` 127-136).

---

# Synthesis — the fix plan

Everything in the levers checks against the code: `logprior` (npmle.py:280-282) is the `np.interp` δ-pin at `log(f_g·M/E)`; `fit` defaults `bandwidth=0.15`, `floor_eps=0.02`, grid top `hi=max(ld)+2h` (so `ρ_g=f_g·M/E ≤ M/E` cannot extrapolate above the total — the structural FP asymmetry is real); `project` (304-314) computes exactly `r_j∝w_j·N(d;μ_j,h²)` (the total-density responsibility MIP would reuse — the §0 hook is confirmed); `_gdna_arm` (161-162) passes `global_logprior` through byte-unchanged; `_rna_arm` (173-175) is the flat `½·log(1−f_g)` reference (the unwritten `logP_r`); `_posterior_median_fg` (127-136) is the location-robust readout. Plan follows.

---

# Projection Fix — Reviewer-Ready Plan

## 0. The yardstick (restated, so every claim below is measured against it)

The proper composition posterior is

```
P(λ | M,u) ∝ BB(u|σ(λ)) · ∫da exp(−½·M·(a−a*(λ))²) · π_g(a) · π_r(a−λ)     a*(λ)=log(f_g·M/E)
```

The shipped ψ (`simplex_logodds.py:267`, `_gdna_arm` + `_rna_arm`; `npmle.logprior:280`) is the **`M→∞` δ-pin × `π_r≡const`** degenerate corner:

```
ψ(λ) = BB(u|σ(λ)) + logπ_g(log(f_g·M/E)) + ½·log(1−f_g)
```

The yardstick exposes exactly **four surfaces**; each candidate lever acts on one:

| Surface in the proper marginal | What it is | Lever(s) |
|---|---|---|
| **S1** estimation width of `π_g` | the KDE bandwidth `h` that fits the population density from Poisson-noisy per-region rates | bandwidth (change A) |
| **S2** the counting kernel `exp(−½M(a−a*)²)` | the node's own counting variance `1/M`, dropped by the δ-pin | bandwidth (change B) |
| **S3** mode-weight ratio `w_dep/w_enr` | the *fitted* height of the enriched vs depleted mode | enrichment-landscape / occupancy-leveling |
| ~~**S4** the RNA arm `π_r`~~ | ~~the discriminant~~ | **REJECTED** — fitting a learned `π_r` failed the full-solve A/B (2026-07-19); the RNA arm stays the flat reference |

The crush is a product of two independent defects on **different surfaces**: an under-smoothing/δ-pin defect (S1+S2, a *width* problem) and a self-reinforcing under-population of the enriched mode (S3, a *weight* problem). This is the crux of the load-bearing/redundant call below.

---

## 1. The correct projection, and which levers are load-bearing

**The correct projection is the un-degenerated marginal: restore the counting convolution on the absolute log-rate axis and read `π_g` through it, not at a point.** Concretely, replace `logprior`'s δ-pin `np.interp` with the analytic Gaussian-convolved mixture re-render at per-node widened variance `H² = (h·ln10)² + 1/M`, and make the fit bandwidth `h` a data-driven estimation width rather than the magic `0.15`. That is the **bandwidth lever, both changes (A fit-`h`, B counting-kernel)** — it is the *only* candidate that is a term-for-term de-approximation of the yardstick (it puts back S1 and S2 exactly, with zero new generative assumptions).

**Load-bearing vs redundant:**

- **Bandwidth (A+B) — LOAD-BEARING, primary, ship-track.** It is the first-principles fix and it alone removes the *smoothing-manufactured* fraction of the depleted advantage `(w_dep/w_enr)·√(1+σ_enr²/h²) → w_dep/w_enr` and broadens the razor-sharp depleted spike by the node's own `1/√M`. Change A does most of the work on the 35 nodes; change B is a principled second-order refinement that only bites at low count (`M≲8`), but it is the honest home of the counting variance and it is what protects the FP guard at sparse `gdna_none` nodes.

- **Enrichment-landscape (skeleton) — CONDITIONAL second lever, NOT redundant with bandwidth.** Bandwidth cannot touch S3: no amount of smoothing rebuilds an enriched mode whose EM *weight* was collapsed by pass-0 (the refit fits on `belief.f_g·mass`, so a pass-0 snap-to-RNA starves the enriched weight; `calibrate.py:104`). If, after optimal `h`, the residual honest ratio `w_dep/w_enr` still out-votes the strand+message counter-evidence (~34:1 in the correct direction, per `refit_loop_study`), you need a weight-rebuild, and the enrichment skeleton is the *principled* one: it sources the enriched-mode existence/location/occupancy from the independent TOTAL landscape and anchors the mode at the absolute `ρ_bg·e^{δ_k}` rate. Its non-redundancy is exactly that it acts on S3 (weight) while bandwidth acts on S1/S2 (width).

- **Height-normalization (tempering `P^{1/T}`) — REDUNDANT/disfavored, keep only as a bounded safety valve.** It targets the same surface (S3) as the skeleton but is *not* derived from the yardstick — `P^{1/T}` appears nowhere in the generative model; it is an ad-hoc reshaping that happens to be monotone-safe. Its one genuinely non-redundant property vs bandwidth (a monotone transform cannot merge two resolved modes, whereas over-large `h` can) is neutralized by simply respecting the `h ≤ Δ_dep-enr/2 ≈ 0.52` ceiling. Prefer the skeleton (derived, absolute-anchored) over tempering (heuristic, unproven non-regression across stranded conditions, auto-`T` perplexity heuristic brushes the no-magic-numbers directive). Reach for tempering only if bandwidth+skeleton are jointly insufficient, and then only with a swept, capped, real-data-validated `T` — never auto-`T`.

- **Mixture-integration (MIP) — REJECTED.** Its distinctive factor is `project()`'s total-density responsibility `r_k∝w_k·N(d;μ_k,h)` (npmle.py:304-314) reused as a *composition-bump weight* — this is precisely the "density votes composition" the owner reset against on 2026-07-18 (`npmle_is_enrichment_not_dna_prior`), and it manufactures a genuine coincidence FP (an RNA-only node whose *total* lands on a real enriched gDNA mode → `f_g→1`) in the unstranded+capture regime where strand cannot rescue. Its correct half (the counting-marginalized integral) is exactly bandwidth change B; the only increment MIP adds beyond that is the §0-violating directional pull. Discard.

**Bottom line:** the answer to "is bandwidth alone enough, or bandwidth + one other?" is **bandwidth is the primary and possibly sufficient fix; enrichment-landscape is a conditional second lever gated on a single cheap measurement (§4, Exp 2) — the post-optimal-`h` residual weight ratio.** Tempering and MIP are not on the ship track.

---

## 2. The exact combined mechanism

All changes isolate to `npmle.py`; `simplex_logodds` ψ-assembly is **unchanged** (`_gdna_arm` at :267/:465 consumes the pre-evaluated `(m,K)` term as-is). `project()` and its `h` are **untouched** so the signed-off σ²_transfer F1 result cannot move.

### Change A — data-driven estimation bandwidth (S1), `DensityNPMLE.fit`

Add a `bandwidth=None` sentinel. When `None`, compute the resolution-floor rule (zero new free constants):

```
g = g_hat[g_hat>0]
h_noise = median(1/sqrt(g)) / ln10                    # decade; the Poisson log-SD floor
# optional DERIVED (not tuned) upper anchor, exposed as a selector, not the default:
ld   = log((g_hat/eff)[g_hat>0])
sigma = (pctl(ld,75)-pctl(ld,25))/1.349 / ln10
h_silv = 0.9 * sigma * ld.size**(-0.2)                # Silverman AMISE
bandwidth = max(h_noise, h_silv)                       # primary default = h_noise alone
```

Everything downstream of `h = bandwidth·ln10` (npmle.py:180) is unchanged. **Applied only to the gDNA/`logprior` fit (`calibrate.gdna_prior`).** The enrichment fit (`calibrate.enrichment_prior` → `project` → σ²_transfer) stays pinned at `bandwidth=0.15`.

### Change B — restore the per-node counting kernel (S2), `DensityNPMLE.logprior`

Replace the δ-pin point read (npmle.py:280-282) with the analytic Gaussian⊛Gaussian re-render at per-node widened variance, using the stored `(weights, log_rho, bandwidth)`:

```
sigma_m = 1/sqrt(mass)                                  # (n,) nat-log counting SD
H2 = (bandwidth*ln10)**2 + sigma_m[:,None]**2          # (n,1) = h² + 1/M
a  = log(fg)[None,:] + (log(mass)-log(eff))[:,None]     # (n,K) = log(f_g·M/E)
# logP_node(a) = logsumexp_j[ log w_j - ½(a-μ_j)²/H2 - ½log(2π H2) ]
```

evaluated over the `f_g` grid → `(n,K)`, replacing the interp. This is the `∫da` of the yardstick done exactly, node-by-node. Chunk over rows like `project()` (the `(n,G,K)` cube never materializes; `G≈200`, `K≈60`). Cap `sigma_m ≤ h` (a 2-count node must not dominate the population width — this cap is the one guard the plan flags as a new heuristic and must be stress-tested, Exp 5). Keep `_posterior_median_fg` (location-robust) as the readout — load-bearing for the FP guard.

### Change C — enrichment skeleton (S3), CONDITIONAL, default-OFF

Only if Exp 2 shows the residual `w_dep/w_enr` still crushes. Three additive edits, mirroring the existing aggregate-background-cell machinery (npmle.py:207-218):

1. `DensityNPMLE.enriched_tiers(self)` — post-fit read of the already-fit TOTAL enrichment prior: render `dens_grid`, find local maxima, lowest = `μ_bg`, higher = tiers; return `[(δ_k=log_rho[peak]−μ_bg, w_k=basin-integrated weight)]`, gated `δ_k > ~h` and `w_k >` a detection threshold. `w_k` is *territory occupancy* (node-count weight from `_collapse`, not fragment mass) — the guard that keeps it shareable DNA↔RNA.
2. `calibrate.py` — `tiers = enrichment_prior.enriched_tiers()`, thread into the gDNA-hyperprior fit.
3. `DensityNPMLE.fit` — inside the `use_bg` block, after the background cell, append one skeleton cell per tier: density `= exp(log_rho_floor + δ_k)` (enriched gDNA mode = background rate × detected enrichment), `wc = w_k·n_regions` (occupancy-set height), into the same `_em_weights`. The fitted `logP` now carries enriched gDNA modes at `ρ_bg·e^{δ_k}` with genome-prevalence height.

### Config / isolation discipline

- `config.py`: `npmle_gdna_bandwidth: float|None = None` (sentinel → data-driven); `enrichment_skeleton: bool = False`; (if ever needed) `npmle_gdna_height_temper: float = 1.0`.
- `calibrate.py`: data-driven `h` and skeleton apply **only** to `gdna_prior`; `enrichment_prior` keeps `bandwidth=0.15`.
- `simplex_logodds.py`: **zero change.** Verify L-invariance and finiteness at `f_g→{0,1}` still hold (the re-rendered arm is proper).

---

## 3. The single principle that bounds FP amplification

> **Every lever may only reshape the *width, weight, or per-node smoothing* of `π_g` on the absolute log-rate axis, where `π_g` is fit on DECONVOLVED `g_hat` (never the total), and every operation must be either monotone/support-preserving OR anchored to the absolute `ρ_bg` scale. Under this constraint a lever can lower an over-dominant mode but can never *create* a mode — so it cannot vote `f_g` up where `π_g` has no mass above `ρ_bg`.**

Why this is the whole game — the **structural asymmetry** (confirmed in code): in `gdna_none`, `π_g` (fit on deconvolved `g_hat`) collapses to the single aggregate-background mode at `ρ_bg`, and the grid top is `hi=max(ld)+2h` so `ρ_g=f_g·M/E ≤ M/E` **cannot extrapolate above the total** (npmle.py:187-189). Every node's ray therefore finds only the sub-total background mode, giving `f_g ≈ ρ_bg/ρ_tot ≈ 0`. No monotone reshaping (bandwidth width, tempering) and no absolute-anchored addition (skeleton at `ρ_bg·e^{δ_k}`, which on an RNA-dominated ray sits at *small* `f_g`) can put an enriched attractor above `ρ_bg` where the population placed none.

This principle is also the clean **disqualifier for MIP**: its proximity factor reads the node's *total* to select which mode votes composition — it violates "never read composition from total," independent of the empirical FP. It is why MIP is rejected on principle, not just on a benchmark.

**Hard preconditions (not tunable knobs):** (a) `π_g` fit on deconvolved `g_hat`, never total — keeps `gdna_none`'s single-mode structure; (b) `ρ_bg` measured on pure intron/intergenic pooled at true zero (`Σg/ΣE`, `background_reference`); (c) median (location-robust) readout, not mean.

**Bounded residual risks, all covered by the diagnostic gate:** over-smoothing past `h > Δ/2` merges the depleted/enriched modes into one mis-composed blur (both A and, in dense conditions, a too-*small* `h_noise` re-sharpening the Everest) — bounded by the **mode-count / inter-mode-dip diagnostic** (refuse/clamp `h` if the dip closes below ~1 nat or the local-max count drops); counting-kernel upward smear at low-count `gdna_none` — bounded by the median readout and the `sigma_m ≤ h` cap; skeleton's `ρ_bg` over-estimate or CNV non-uniformity — real-data caveats (§5).

---

## 4. Ranked prototyping plan

All on the ambig sim (`gdna5_capON`), acceptance test = **35 partially-enriched nodes' median `f_g` → 0.63** AND **`gdna_none` median `f_g` → 0 (<~0.02)**. Cheapest-highest-signal first.

**Exp 1 — Single-node re-score curve (near-zero cost, no full run; HIGHEST signal).**
Load the fitted `gdna5_capON` `π_g` and node 2589 + the 35-node totals. Plot the gDNA arm `logP_g(f_g)` over the grid for: (a) shipped δ-pin interp; (b) `h`-sweep {0.10,0.15,0.20,0.30,0.40,0.50}; (c) counting-kernel read `σ_m=1/√M` ON vs OFF; (d) skeleton injected. Read off whether the argmax moves off `f_g≈0.18` toward the enriched mapping, and *how much* of the ~4-nat crush each change removes. This one experiment tells you the shape of the answer.

**Exp 2 — Residual-weight decomposition (cheap, analytic on the fit; the GO/NO-GO for skeleton).**
Decompose the rendered depleted:enriched ratio into `(w_dep/w_enr)·(σ_enr/h)`: measure `σ_enr, w_enr, w_dep`. Compute the **floor ratio at `h=Δ/2≈0.52`** (max safe `h`). If floor ratio ≪ the ~34:1 strand+message counter-vote → **bandwidth alone suffices**, skeleton not needed. If floor ratio still out-votes → **enrichment-landscape (Change C) is required.** This single measurement decides "bandwidth alone vs bandwidth+skeleton."

**Exp 3 — Over-smoothing / merge-cliff diagnostic (cheap).**
Across the `h`-sweep, track local-max count and inter-mode dip depth; confirm the enriched pair (−2.01/−1.64, 0.37 decade) stays resolved and a safe window `σ_enr ≲ h ≲ 0.52` exists. Establishes the shippable over-smoothing detector and the `h` ceiling. Also run on a dense/high-count condition to confirm `h_noise` does **not** drop below ~0.15 and re-sharpen the Everest (the bandwidth key_risk).

**Exp 4 — In-process A/B full solve, `gdna5_capON` (moderate; `calibration-benchmark` skill).**
Wire bandwidth A+B behind the flag; re-solve; report the 35-node refit `f_g` distribution (median → 0.63) and net fragment-flow, and confirm the well-solved single-strand nodes do not regress. Add Change C only if Exp 2 said so, then re-run.

**Exp 5 — Specificity guard, `gdna_none` full A/B (MUST PASS — acceptance gate).**
Same fixes on `gdna_none`: assert median `f_g < ~0.02` across the whole `h`-sweep (+skeleton if on). Include the FP-invariance numerical check (the ray stays monotone-decreasing in `f_g`) and the low-count stress: down-sample counts to inflate `σ_m`, verify the `σ_m≤h` cap + median readout hold `f_g≈0`.

**Exp 6 — Data-driven selector agreement, no oracle peeking (cheap).**
Compute `h_noise`, `h_silv`, and an LCV (leave-one-cell-out Poisson predictive) argmax; confirm the selected `h` lands inside the Exp-3 safe window **without** seeing the oracle. Acceptance: a purely data-driven selector picks a good `h` on `gdna_none`, `gdna5_capON`, and `gdna300` alike.

**Exp 7 — Full 24-condition net-flow regression (most expensive, LAST).**
Especially `gdna300_ss0.99` and the unstranded-resolvable conditions — the non-regression the height-normalization/bandwidth verdicts insist on. Gate any default change here, not on `gdna_none` alone.

---

## 5. Sims-now vs real-data

**Settle on sims now (the CI suite has the true-zero `gdna_none` and controlled `σ_enr`):**
- Whether the crush is smoothing-manufactured (S1/S2) vs weight-under-population (S3) — the load-bearing architecture call (Exp 1+2).
- The bandwidth mechanism itself (A data-driven `h`, B counting kernel) — pure math, validate end-to-end.
- FP-invariance / `gdna_none` specificity and the FP-bounding principle (Exp 5).
- The over-smoothing diagnostic (mode-count/dip) and the safe `h` window / merge cliff (Exp 3).
- Whether bandwidth alone suffices or the skeleton is needed (Exp 2+4) and the skeleton's detectability + `ρ_bg`-anchor FP safety on the sim's total landscape.
- The data-driven `h` selector's agreement with the oracle frontier and its downward-misfire mode in dense conditions (Exp 6).

**MUST wait for real data:**
- **The final numeric `h` selector** (`h_noise` vs `h_silv` vs LCV). Ship the zero-constant `h_noise` rule as the default; real capture chemistry sets the true enriched-mode widths/tier structure, so finalize the selector on real libraries.
- **The enrichment skeleton's strength `wc=w_k·n_regions`** (the free `n_eff` scale) and its **detectability** (can 35 gDNA-enriched nodes form a total-landscape mode distinct from the RNA-expression continuum?). Keep Change C default-OFF until validated/derived on real data; if it only ever resolves RNA-expression tiers it stays FP-safe but fails to recover — that's an acceptable, non-breaking null.
- **CNV / tumor** breaks the `ρ_g=ρ_bg·e` uniformity premise the skeleton relies on (copy-number amplification is locus-specific DNA that is not capture enrichment → shared-tier detection mis-locates enriched gDNA modes). Cannot be settled on the CNV-free sims; **hard caveat, flag skeleton for tumor libraries.**
- **Tempering `T`**: do **not** ship `T>1` as a default; its non-regression across real stranded conditions is the unproven risk, and auto-`T`'s perplexity heuristic violates the no-magic-numbers directive. Keep it off the ship track unless real-data net-flow shows the bandwidth+skeleton residual still crushes, and then only swept+capped.

**Hard constraints settled now, both sims and real:**
- `project()`'s bandwidth stays pinned at the signed-off `0.15` — do not perturb σ²_transfer F1.
- Isolation: data-driven `h` and skeleton on `gdna_prior` only, never `enrichment_prior`.

## 6. Rejected track (kept for provenance): the fitted RNA arm `logP_r`
Surface **S4** — replacing the flat `π_r` reference (`_rna_arm`, `½·log(1−f_g)`) with a **fitted** `logP_r` — was
prototyped and **rejected** on 2026-07-19. A full-solve A/B (`calib_rna_arm ∈ {off, raw, tail}`) showed: the raw
density form does **not** recover the enriched nodes (messages anchor the snap; the render-only +0.11 vanishes)
*and* adds false positives; the FP-bounded upper-tail form is monotone in `log ρ_r` ⇒ a one-directional gDNA-ward
push ⇒ catastrophic over-call. The snap is **message-and-gDNA-arm-driven**, so a local RNA arm cannot overcome it.
This track is **closed** and is not a calibration option. The RNA constraint stays the reference + the messages;
the live work is the gDNA-hyperprior *representation* (see `CALIBRATION_MASTER.md` and
`kde_vs_npmle_enriched_mode`).