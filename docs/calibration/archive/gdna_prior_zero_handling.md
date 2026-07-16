# The gDNA prior: the depleted floor, the KDE, and the zero / discreteness problem

**Status:** design + open questions, for review (branch `calib-ambig-init-wip`). Companion to
`node_prior_design.md` (the projection prior) and `calibration_prior_production_reference.md` (what ships).
**Scope:** how the gDNA-density prior behaves when a sample has **little or no gDNA**, why the current
design bolts a depleted scalar floor onto the KDE, and a parameterization question that real data forces:
should the prior be modeled over **density** (`count/E`) at all, or over the **latent rate** with a
Poisson observation of `(count, E)`?

All evidence here is from real hybrid-capture samples (three cfRNA, one ~20–25 % gDNA in-silico mix),
scan-cached (`scripts/debug/calib_cache.py`) and measured with `scripts/debug/pass0_kde_zero.py` /
`landscape_real.py`. Real STAR BAMs have **no per-node oracle** — the arguments are structural, not
truth-labeled.

---

## 1. The current architecture — a depleted scalar floor *plus* the KDE

Two priors on `f_g` are **added** in the sweep (`bp_solver.node_sweep`):

1. **`_global_logprior` — the depleted background floor, always on (both passes).** `_floor_estimate`
   pools an **exposure-weighted scalar gDNA density** `ρ_floor = (1+G)/E_tot` over ALL intergenic + intron
   **regions** (boundaries and exons excluded), *keeping zero-count regions* (0 gDNA over a large E is the
   strongest evidence gDNA is scarce → `ρ_floor → 0⁺`). `_global_logprior` applies `ρ_floor` as the per-node
   target for intergenic/intron nodes (with a data-driven density-likelihood + a Jeffreys RNA-parsimony
   term), and a weak genome-wide baseline (capped at one pseudo-observation) elsewhere.
2. **`_kde_logprior` — the density KDE, pass-2 only, ADDED on top.** `P(log ρ_g)` (real Gaussian tails) +
   the same Jeffreys term.

**Why the floor is bolted on:** the KDE + Jeffreys tilt toward gDNA (the Jeffreys `−log(1−f_g)` favors small
RNA ⇒ large `f_g`; the KDE cannot go below `1/E`), so the depleted population would be over-called. The
scalar floor drags it back. It is a *corrective patch* holding two priors in tension — the thing we want to
unify away.

---

## 2. The zero problem — and how large it is on real data

**Mechanism.** The KDE is fit on densities floored at the min-observable `1/E`, so its support is strictly
positive: it **cannot represent `ρ_g = 0`**. Project a truly gDNA-free node and it lands on the `1/E` floor
→ `f_g = 1/M > 0` → **false-positive gDNA**. A pristine sample gets gDNA injected everywhere the strand
cannot veto it. The scalar floor mitigates this (a low target) but is a fixed level, not zero-aware.

**Scale (measured).** The confident gDNA set (intergenic + single-strand introns + intergenic↔exon +
single-strand intron↔exon boundaries) is *dominated by zeros*:

| sample | fraction with **count = 0** | count ≤ 3 |
|---|---|---|
| vcap (~25 % gDNA) | **64 %** | 90 % |
| LBX0588 (cfRNA) | **82 %** | 97 % |
| LBX0190 (cfRNA, near-pristine) | **94 %** | 99.6 % |

The density-space training substrate requires `mass > 0`, so it **discards 64–94 % of the confident
nodes** — precisely the "gDNA is absent here" observations. The KDE is trained on the small remainder that
*has* gDNA and is structurally blind to absence. That is a systematic over-call by construction, and it is
worst exactly where it matters (the near-pristine cfRNA).

---

## 3. Options for zero / sub-`1/E` sparse gDNA

- **3.1 Depleted scalar floor (current).** A single pooled low level added on top. Simple, load-bearing
  ([[calibration_collapse_cleanup_locked_in]]), but not zero-aware and structurally a patch.
- **3.2 Zero-inflated KDE (spike-and-slab).** `prior(ρ_g) = π·[ρ_g→0] + (1−π)·KDE(ρ_g)`, with a δ-atom at
  zero. **π = the population fraction consistent with zero gDNA** — directly the count-0 fraction above
  (0.64 / 0.82 / 0.94), or fit by EM. Robust: pristine → `π→1` → `f_g→0`; contaminated → the KDE modes take
  over. Subsumes the floor (its job = "high π + a low mode"). *Caveat:* still a density-space model for the
  nonzero part, so it inherits the discreteness of §4.
- **3.3 Censored-below.** Treat `1/E`-floored points as "≤ 1/E" (left-censored) so the tail reaches zero.
  Works; less explicit than 3.2.
- **3.4 Poisson-rate model (§4).** The root fix — makes zero *native* rather than special-cased.

---

## 4. The parameterization question: density vs the latent rate (the count–length reframe)

Density `ρ = count / E` collapses two discrete quantities, and real data shows the collapse is not benign:

- **Zeros are unrepresentable.** `count = 0` → `ρ = 0`, which the KDE floors to `1/E`. The 64–94 % of
  confident nodes at count 0 cannot be modeled in density space at all (§2).
- **The nonzero mass is quantized.** 90–99.6 % of nonzero confident nodes have **count ≤ 3**, and
  **boundaries sit at short, near-fixed E** (median 100–260 bp ≈ a fragment length) vs regions at ~1.6–1.8
  kb. So a boundary's density is `k/E` for tiny integer `k` over an almost-constant `E` — a **comb of
  discrete atoms** (visible directly in `pass0_kde_zero.py` Panel B). The KDE's bandwidth smears that comb
  into a smooth "enriched mode" — i.e. it can mistake **quantization for enrichment**.

**The generative truth.** A node's count is `k ~ Poisson(ρ·E)` for a latent gDNA rate `ρ`, with
`ρ ~ P(ρ)` the population rate distribution we actually want as the prior. The density KDE *approximates*
`P(ρ)` by kernel-smoothing the point estimates `k/E` — which (a) cannot place mass at `ρ=0`, (b) ignores
the `E`-dependent precision (a `k/E` from `k=1,E=100` is treated like `k=100,E=10⁴`), and (c) turns
integer quantization into false structure.

**Model `P(ρ)` directly (Poisson-rate deconvolution).** Estimate the latent rate distribution such that
each confident node's `(k_i, E_i)` is `Poisson(ρ_i·E_i)` with `ρ_i ~ P(ρ)`. Then:
- **Zero is native** — `k=0` is a legitimate observation supporting `ρ→0`; the count-0 nodes *are* the
  absence signal, so zero-inflation is automatic (no atom, no π to fit).
- **Discreteness is native** — the Poisson is the correct discrete likelihood; no quantization artifact.
- **Precision is native** — high `E`/`k` → sharp; low → broad. The honest per-node precision falls out of
  the likelihood, so the separate `σ²`/bandwidth machinery is largely subsumed.

**On the literal "2D/3D KDE over (count, length)":** the `(count, E)` plane is the *observation* space
(worth plotting — it exposed all of the above), but a raw 2D KDE over it estimates `P(count, E)`, the joint
of *observations*, not the prior `P(ρ)` we need. The principled realization of the idea is the **Poisson
deconvolution to `P(ρ)`**, of which the 2D view is the input. So: the instinct is right and points at the
Poisson-rate model; the deliverable is `P(ρ)`, not a joint KDE.

**Cost / risk.** A Poisson-rate deconvolution (e.g. an EM/NPMLE mixture of rates, or a binned
Poisson-likelihood KDE) is a bigger change than the density KDE, costs more compute, and needs care at
genome scale. But it **unifies** the floor (3.1), the zero-inflation (3.2), and the precision/bandwidth
handling into one correctly-specified model.

---

## 5. Recommendation + open questions

1. **Near-term (cheap, density-space):** prototype the **zero-inflated KDE** (3.2) with `π` = the measured
   count-0 fraction, and show the projection's left side reaches `f_g→0` on near-pristine LBX0190 (does the
   false-positive floor vanish?). Low-risk, uses existing machinery.
2. **Medium-term (root fix):** prototype the **Poisson-rate `P(ρ)`** (§4) on the caches and compare to the
   density KDE — does it (a) use the count-0 majority, (b) remove the fake boundary "enriched mode," (c)
   handle zero without a floor? This is the candidate to replace the floor + KDE + zero-inflation together.
3. **Substrate weighting (deferred by §4):** the density substrate uses **unit weight** (design §8e:
   precision correlates with density, so precision-weighting biases the shape toward enriched; precision's
   correct home is the bandwidth). Worth an A/B on real data — but the Poisson-rate model makes it moot (the
   likelihood *is* the precision).

**Open questions for review.**
- Is the discreteness/comb severe enough to distort the *enriched-mode detection* we rely on, or does the
  bandwidth adequately absorb it? (The boundary comb in `pass0_kde_zero.py` Panel B suggests the former.)
- For the Poisson-rate model, what estimator at genome scale — a binned Poisson-likelihood grid over `log ρ`
  (cheap, O(nodes + grid)) or an NPMLE mixture (sharper, costlier)?
- Does the count-0 majority argue for fitting `P(ρ)` on **all** confident nodes (zeros included) from the
  start — the elegant "all nodes" the density substrate can't support because it drops the zeros?
