# The gDNA intron factory — concept → design → implementation plan

**Status:** design (2026-07-20). Promotes the concept in [`gdna_intron_factory.md`](gdna_intron_factory.md) (owner's
notes) to a worked design + implementation plan. Grounded in [`background_reference_derivation.md`](background_reference_derivation.md)
(the aggregate `ρ_bg`), [`gdna_background_floor_derivation.md`](gdna_background_floor_derivation.md) (the
resolution wall), [`CALIBRATION_ARCHITECTURE.md`](CALIBRATION_ARCHITECTURE.md) §0–1 (count-zero-information), and
**empirically pre-validated** on the cached `ambig_dense_10mb` scenarios (§7 — the numbers already track oracle).

Supersedes the "nascent factory" naming; the mechanism **peels gDNA, not RNA** (owner's key reframe), so it is
the *gDNA* factory.

---

## 0. TL;DR

Introns carry **background gDNA + (scarce) nascent RNA**. Introns are **off-target** under capture, at the *same*
depletion as intergenic regions — so the pooled intergenic gDNA rate `ρ_bg` **is the intron's own gDNA density**,
not merely a floor. For each intron we therefore deconvolve its observed unspliced count `C` (over gDNA
eff-length `E_g`) into gDNA (explained by the background) + nascent (the excess), producing a **belief on
`λ = logit(f_g)`** — the gDNA-vs-RNA-total log-odds — with **honest, count-derived precision**. We **peel gDNA**
(strand-symmetric ⇒ no strand to assign); the residual RNA's strand is left to the solver (`θ`, the tilt). The
belief enters the per-node solve as one additive `(n_nodes, K)` factor on `λ` — the *same* slot the gDNA
composition arm already occupies — and (crucially) **propagates**: a confident intron gDNA anchor flows
genomically to the adjacent exons/boundaries, where the intergenic background is *not* a valid reference.

This is the §1.3 global gDNA hyperprior **done right**: scoped to the nodes (introns) where the intergenic
reference is statistically valid, with a two-sided estimate instead of the exon-only one-sided floor.

---

## 1. The problem, and what changes

Calibration deconvolves each node's unspliced mass into `(f_pos, f_neg, f_g)`. The **hard regime is unstranded**
(`κ = ½`): the strand Beta-Binomial Fisher information is `∝ N·(2κ−1)² = 0`, so a node has **no intrinsic
composition signal** (`CALIBRATION_ARCHITECTURE.md` §1.1). Today, pass-0 is prior-free there, so an unstranded
intron sits at the reference `f_g ≈ ½` — regardless of the true DNA level. The consequences (measured, §7):

- **gDNA-free library** (`gdna_none`): intron truth `f_g = 0` (pure nascent), but pass-0 reads `≈ ½` → a
  **fabricated gDNA false-positive** that seeds the phantom NPMLE enriched mode ([`pass0_fake_enriched_mode_diagnosis`]).
- **gDNA-rich library** (`gdna300`): intron truth `f_g ≈ 0.88`, pass-0 reads `≈ ½` → a **gDNA under-call**.

The node's own composition (total density, strand) is *identical* in the two cases — the truth is opposite only
because the **DNA level** differs. No fixed prior can be right for both (`background_reference_derivation.md` §1);
**the DNA level is a measured DATA quantity.** The gDNA intron factory measures it, per intron, from the intron's
count relative to the known background — and this is exactly the missing signal that resolves both corners.

**What changes.** For **intron** REGION nodes only, add a factor on `λ` derived from the deconvolution below.
Nothing else changes (exons, boundaries, intergenic locks, the strand likelihood, the message system are all
untouched). Introns are already *solvable* (G2 single-strand / G3 AMBIG — never G1-locked, verified §7); the
factory **informs** them, it does not unlock anything.

---

## 2. Architectural placement — the §1.3 hyperprior, scoped, and the §8 reconciliation

The factory is the **global gDNA hyperprior** (`CALIBRATION_ARCHITECTURE.md` §1.3: "the population baseline gDNA
density … sways an unanchored node toward gDNA when there is global evidence of gDNA elsewhere; precision = the
robust population spread"). The *current* hyperprior (Role B, `_fit_gdna_hyperprior`) is fit on the deconvolved
gDNA of all nodes and fails under capture, because a captured exon's on-target gDNA density is `≫` the intergenic
floor — the reference is wrong there. The factory fixes the **scoping**: apply the background reference only where
it is valid.

**Why introns, and the reconciliation with `background_reference_derivation.md` §8.** §8 established that `ρ_bg`
(pooled intergenic + RNA-free intron gDNA) measures `ρ_DNA · c_off` (DNA × off-target capture efficiency), a
**one-sided lower bound** that must NEVER be used as `f_g = ρ_bg/ρ_tot` — *for an exon*. The reason is that an
exon is **on-target**: its true gDNA density is `ρ_bg · e` with enrichment `e ≫ 1`, so `ρ_bg/ρ_tot` would
false-negative a real gDNA exon.

An **intron is off-target** (`e ≈ 1`), at the *same* depletion as the intergenic pool. So for an intron the
enrichment cancels and `ρ_bg` is not a floor but the **actual gDNA density**:

```
    ρ_{g,intron}  =  ρ_bg · e_intron  ≈  ρ_bg · 1  =  ρ_bg ,     hence   f_g^{intron} ≈ ρ_bg / ρ_obs   (two-sided, valid).
```

The empirical check (§7) confirms this to the digit: intron truth `f_g` ≈ `ρ_bg/ρ_obs` across capture on/off.
**This is the whole reason the factory is intron-scoped** — introns are the nodes where §8's forbidden
`ρ_bg/ρ_tot` becomes the *correct* two-sided estimate. Exons still get gDNA only via propagation + on-target
evidence, never from `ρ_bg` directly.

**Graceful degradation is inherited from §8.** As capture strengthens, the off-target pool depletes, `ρ_bg → 0`
and `Var(log ρ_bg) = 1/Σg → ∞`; the factory's gDNA peel becomes *imprecise* (wide) rather than confidently
wrong. A genuinely DNA-free library (`Σg = 0` exactly) is the precise-near-zero case (the resolution wall,
§4) → `f_g → 0` correctly. The factory never manufactures gDNA (no false-positive) nor confidently crushes it
under depletion (the precision widens) — the §8 safety properties, restricted to the domain where the estimate
is two-sided.

---

## 3. The statistical model — deconvolving one intron

**Setup.** Intron `i`: observed integer unspliced count `C_i = n_unspliced_pos + n_unspliced_neg` (both strands
— gDNA is strand-symmetric), gDNA effective length `E_{g,i} = region_eff_length(L_i, gdna_fl_pmf)`. Observed
density `ρ_obs,i = C_i / E_{g,i}`.

**Physical decomposition.** `C_i = g_i + r_i`:
- **gDNA** `g_i`: `g_i ~ Poisson(ρ_{g,i}·E_{g,i})`, with the intron's gDNA density drawn from the background,
  `ρ_{g,i} ~ π_bg` (the intergenic gDNA density distribution — §8: intergenic ≡ intron off-target rate).
- **nascent RNA** `r_i ≥ 0`: the excess. Nascent is **scarce and region-specific** — it does *not* correlate
  with any global pattern (owner's open Q1), so it carries **no informative prior**; we use the maximally
  non-committal **flat one-sided prior** `π_r(r) ∝ 1, r ≥ 0`. (This is what dodges the "magic cutoff" the
  concept doc worried about — see below.)

**The posterior on the gDNA count**, marginalizing `r` under the flat prior (`∫₀^∞ Poisson(C; (ρ_g+ρ_r)E) dρ_r`
makes any `C ≥ g` equally likely given `g`):

```
    P(g | C_i, E_{g,i})  ∝  P_bg(g)  ·  1[0 ≤ g ≤ C_i] ,                                             (3.1)
    where   P_bg(g)  =  ∫ Poisson(g; ρ·E_{g,i}) · π_bg(ρ) dρ                                         (3.2)
```

**`P_bg` is the background *count* distribution over the intron's eff-length**, truncated at the observed total.
`f_g = g/C_i`, so (3.1) is directly a posterior on `f_g` (equivalently `λ = logit(f_g)`). **No cutoff, no tuned
constant** — the flat nascent prior + the truncated background gives a smooth posterior. This is (3.1)'s answer
to open Q1: we do **not** need an "intronic distribution", and there is no threshold; the deconvolution is the
posterior (3.1).

**The three regimes (all handled by one formula, verified §7):**

| regime | `ρ_obs` vs `ρ_bg` | posterior mode `g*` | `f_g` | meaning |
|---|---|---|---|---|
| pure gDNA (no nascent) | `ρ_obs ≈ ρ_bg` | `≈ C` (truncation binds) | `≈ 1` | all gDNA — confident |
| nascent present | `ρ_obs ≫ ρ_bg` | `≈ ρ_bg·E_g` | `ρ_bg/ρ_obs < 1` | **peel gDNA at background; excess is nascent** |
| below background (low mappability) | `ρ_obs < ρ_bg` | `= C` | `≈ 1` | all gDNA (faint locus) |
| DNA-free library | `ρ_bg ≈ 0` (Σg=0) | `≈ 0` | `≈ 0` | pure nascent RNA — **no false-positive** |

The concept doc's "how many fragments are *confidently* gDNA?" is `E[g | C]` under (3.1); the residual `C − E[g|C]`
is the (imprecisely known) nascent.

---

## 4. Precision — honest, discrete count over discrete length

The owner's hard constraint: the deconvolution must carry a precision that **honors the discrete count over
discrete length**, is **confident (high) for gDNA in most cases**, and is **imprecise about peeling nascent**
(intergenic vs intronic distributions overlap). The posterior (3.1) delivers exactly this, from **two honest
uncertainty sources, no magic number**:

1. **Per-intron Poisson (count-over-length).** In (3.2) the `Poisson(g; ρE_g)` factor makes the gDNA count
   uncertain to `±√(ρ_bg·E_g)`. A high-count intron (dense background × long intron) resolves `f_g` sharply; a
   low-count intron is diffuse — `Var(f_g) ≈ ρ_bg·E_g / C²`, the honest Fisher precision of a rate over an
   opportunity. This is `CALIBRATION_ARCHITECTURE.md` §1.3's "count enters only as statistical power".
2. **Background spread.** `π_bg` has width from (a) the pooled-mean uncertainty `σ_bg = 1/√Σg` (tiny when the
   genome-scale `Σg` is large — §7: `σ_bg ≈ 5e-4`) and (b) the per-region biological/mappability/CNV variation
   (the over-dispersion of intergenic densities). When an intron sits *just above* `ρ_bg`, this spread is what
   makes the nascent peel imprecise (the excess could be a high-density gDNA region or gDNA+nascent) — exactly
   the owner's expected "imprecise near the overlap". Far above `ρ_bg`, the spread is irrelevant and the peel is
   confident.

**Why "confident gDNA in most cases":** most introns sit at `ρ_obs ≈ ρ_bg` (nascent scarce), where (3.1)
concentrates `f_g` near 1 with the tight background precision. Nascent-bearing introns get a confident *gDNA*
floor (`g ≈ ρ_bg·E_g`, precise) and an *imprecise* nascent residual — precisely the asymmetry requested.

The factor's precision is read straight off the curvature of the tabulated `log P(f_g)` on the solve grid — the
existing solver already returns `gdna_frac_var = Var(log f_g)` and `lam_var = Var(λ)` by grid moment-matching, so
the belief is `(mode, precision)` with no extra work.

---

## 5. Strand and DOF — peel gDNA (λ), leave the tilt (θ) free

The factory constrains **only `λ`** (gDNA-vs-RNA-total). This is the owner's "peel gDNA, not RNA" made exact:
gDNA is strand-symmetric, so peeling it needs no strand; the residual RNA's `+/−` split is the tilt `θ`, resolved
elsewhere. This maps *perfectly* onto the solver's existing `(λ, θ)` coordinate model:

- **Single-strand intron (G2, 510/573 in §7):** one live strand carries `f_active = 1 − f_g`. The factory fixes
  `f_g`; the residual `(1−f_g)` lands on the inferred strand automatically. 1 DOF, fully resolved.
- **AMBIG intron (G3, 63/573):** the factory fixes `f_g` (the λ axis = "this intron is X% gDNA"); the residual
  RNA's `+/−` split (`θ`) is left to the strand likelihood (flat at `κ=½`) + neighbour messages + solver. The
  factory informs the "tilt toward gDNA" and is silent on the strand of the residual — exactly the owner's
  description. 2 DOF, λ resolved, θ free.

Mechanically this is automatic: a factor on `λ` is a `(n_nodes, K)` array over the `f_g = σ(λ)` grid; for AMBIG
the solver broadcasts it across all tilt columns (`[:, :, None]`) and integrates `θ` out by `logsumexp` — so a
λ-only factor is θ-free **by construction** for both node classes, no special-casing.

---

## 6. Count-zero-information compliance (why this is not a §0 violation)

`CALIBRATION_ARCHITECTURE.md` §0: *a fragment count carries ZERO intrinsic gDNA/RNA composition information; the
count may enter ONLY as statistical power, never a composition vote.* The factory obeys this:

- The **gDNA amount** `g ≈ ρ_bg·E_g` comes from the **external** intergenic background (`ρ_bg`), an outside
  §1.3 quantity — **not** from the intron's own count voting itself gDNA.
- The intron's count `C` enters **only** as (a) the fraction **normalizer** `f_g = g/C` (the definitional
  density↔fraction jacobian that §2 explicitly blesses — "the same `M` that *defines* the fraction") and (b) the
  **Poisson precision** (§4). It never votes the gDNA amount.

So the factory is the legitimate §1.3 channel (population baseline gDNA sways an unanchored node), scoped to the
nodes where the population baseline is the valid reference. The distinction from the forbidden count-prior is the
same as everywhere in the architecture: `ρ_bg` (external) sets the amount, `C` (own) sets normalizer + power.

---

## 7. Empirical pre-validation (cached `ambig_dense_10mb`, pass-0 vs oracle)

Before any implementation, the point estimate `f_g = ρ_bg/ρ_obs` was checked against oracle intron `f_g`
(`measure_background(include_introns=False)` for `ρ_bg`; per-intron `C/E_g` for `ρ_obs`):

| condition | `ρ_bg` (intergenic) | `Σg` (→ `σ_bg`) | intron `ρ_obs/ρ_bg` (median) | factory `f_g=ρ_bg/ρ_obs` | **oracle intron `f_g`** |
|---|---|---|---|---|---|
| gdna300 ss0.50 capOFF | 0.568 | 3.8M (5e-4) | 1.14 | **0.875** | **0.885** |
| gdna300 ss0.50 capON | 0.052 | 352k (1.7e-3) | 1.07 | **0.93** | **0.959** |
| gdna_none ss0.50 capON | 0.000 (Σg=0) | 0 (∞) | — | **≈ 0** | **0.000** (78% of introns `<0.5`) |

Intron node classes are **G1=0 / G2=510 / G3=63** in every condition — introns are solvable, never locked. The
factory's estimate tracks oracle across the capture cliff *and* both DNA extremes — the exact regime where the
current prior-free pass-0 fails (sits at `½`). The background is count-rich even under capture (`σ_bg ≤ 2e-3`), so
the peel is confident. **This is the mechanism CALIBRATION_STATUS's TL;DR is asking for**, arriving as a
first-class density factor rather than the un-landed subtraction.

---

## 8. The background distribution `π_bg` — design options + recommendation

The factory needs `P_bg(g)` = the intergenic gDNA **count** distribution over an eff-length `E_g` (3.2). Three
candidate representations, in increasing fidelity/cost. **No option introduces a tuned constant** (`ρ_bg`, the
over-dispersion, and the mixture weights are all *fit* from intergenic data).

**Option A — scalar `ρ_bg` + pure Poisson.** `π_bg = δ(ρ − ρ_bg)`, `ρ_bg` from `measure_background`, so
`P_bg(g) = Poisson(g; ρ_bg·E_g)`. Simplest; captures the per-intron Poisson precision exactly. **Misses the
per-region background spread** (CNV/mappability) → over-confident nascent peel where intergenic densities vary.
Add `σ_bg` as a lognormal jitter on `ρ_bg` to recover the pooled-mean uncertainty (still misses per-region
variation). Good first cut; the sim's intergenic is near-Poisson so this may already validate.

**Option B — Negative-Binomial (Gamma-Poisson) [RECOMMENDED].** `π_bg = Gamma(mean ρ_bg, over-dispersion φ)`,
so `P_bg(g) = NegBinom(g; mean = ρ_bg·E_g, φ)` — closed form. `ρ_bg` from `measure_background`; **`φ` fit by
method-of-moments from the intergenic per-region counts** (the excess of `Var(g_i)` over the Poisson mean — the
same style as the strand overdispersions `od_g`, `od_r`). Captures **both** honest uncertainty sources (§4):
Poisson (per-intron) + over-dispersion (per-region background variation). `Var(g) = μ + μ²/φ`. On real data with
CNV, `φ` finite widens the peel; on a Poisson-clean sim, `φ → ∞` reduces to Option A. Clean, honest,
one fitted shape parameter, no magic number.

**Option C — reuse `DensityNPMLE` fit on intergenic.** Fit `DensityNPMLE.fit(g_intergenic, E_g_intergenic)`
(the existing engine, Role-B path) → the full non-parametric `π_bg` including multimodal CNV structure, then the
factor is `DensityNPMLE.logprior(grid, C_intron, E_{g,intron})` — the *existing* `(n,K)` producer, no new math.
**Caveat:** `logprior` scores the intron's *point* density `log(f_g·C/E_g)` against the deconvolved `π_bg`
(density), so it omits the *per-intron* Poisson broadening (over-confident for low-count introns). The honest
form convolves per intron: `log Σ_j w_j · Poisson(f_g·C; ρ_j·E_g)` (a Poisson-marginalized mixture, still one
pass over the grid). Highest fidelity; most code.

**Recommendation:** **Option B (NegBinom)** as the production model — it is the smallest honest model that carries
both uncertainty sources and reuses `measure_background` for `ρ_bg`. Validate Option A first (fastest, confirms
the concept end-to-end), then B. Keep Option C in reserve for real-data CNV multimodality. This is the one
genuine **modeling choice** for the owner to ratify (per the no-new-heuristics rule) — see §9.

**Subtlety — the background pool substrate.** `measure_background(include_introns=False)` pools **intergenic
only** (avoids circularity: introns are what we deconvolve). On real data, `include_introns=True` adds
RNA-free introns for resolution (nascent is sparse, `robust_trim_mad` fences the rare transcribed intron). For
the factory, intergenic-only is the clean default; document the circularity trade-off for the real-data switch.

---

## 9. Open design decisions (owner)

1. **`π_bg` representation (§8):** ratify **Option B (NegBinom, fitted `φ`)** for production, A for first
   validation. The over-dispersion `φ` is a *fitted* quantity, not a constant — but its estimator (MoM vs ML,
   and whether to floor it) is a modeling choice worth a nod.
2. **Nascent prior:** ratify the **flat one-sided** prior (§3) — the choice that avoids a cutoff. (Alternative:
   a weak `r → 0`-favouring prior, but that needs a scale = a magic number; not recommended.)
3. **Pipeline placement & interaction with the existing hyperprior (Role B refit) — RESOLVED (owner,
   2026-07-20):** the factory resolves **intronic gDNA BEFORE the pass-0 solve** (a pre-computed λ-factor that
   enters the pass-0 ψ, so pass-0 solves introns *with* their gDNA already peeled), and the now-resolved introns
   then become part of the **initial gDNA hyperprior fit after pass-0 completes** (they join the
   `_fit_gdna_hyperprior` training substrate as clean gDNA nodes, instead of the current garbage `f_g≈½`
   unstranded introns that pollute it). This is option (c): the factory *seeds* the hyperprior with a clean
   intron substrate. No double-counting — the factory is the intron composition prior at pass-0; the hyperprior
   consumes the resulting deconvolved intron gDNA as training data, exactly as it already consumes single-strand
   exons.
4. **Does the factory apply to intergenic too?** Intergenic is already G1-locked to `{0,0,1}` (all gDNA) — the
   factory would agree (`ρ_obs ≈ ρ_bg ⇒ f_g ≈ 1`) but adds nothing. Keep the lock; scope the factory to
   `coarse_type == INTRON`.

---

## 10. Implementation plan

Phased, each phase independently validatable on the cached scenarios. Insertion points are exact (from the code
trace). **No new magic numbers**; one fitted shape parameter (`φ`, Option B) with the same status as `od_g/od_r`.

### Phase 1 — the factor primitive (pure, unit-testable)
New module `src/rigel/calibration/gdna_intron_factory.py`:
- `fit_intron_background(substrate, region_arrays, region_eff_g, *, config) -> IntronBackground` — wraps
  `measure_background(include_introns=False)` for `ρ_bg, σ_bg, log_rho_floor, Σg`, and (Option B) fits `φ` by
  MoM from the intergenic per-region counts `(n_unspliced_pos+neg)` over `region_eff_g`. Returns a small frozen
  dataclass.
- `intron_lambda_factor(background, C, E_g, fg_grid) -> ndarray[(n_intron, K)]` — the tabulated
  `log P_bg(f_g·C) + log 1[·]` (3.1) over the `f_g = σ(λ)` solve grid, using `P_bg = NegBinom(mean=ρ_bg·E_g, φ)`.
  Normalized per node (constant in `λ` is irrelevant to ψ). Vectorized; the graceful `ρ_bg → 0` / `Σg = 0`
  behaviour (fall back to the resolution wall `log_rho_floor` as the gDNA upper bound, wide precision) lives
  here.
- Unit tests: the three regimes of §3 (pure gDNA `→ f_g≈1`; nascent `→ f_g≈ρ_bg/ρ_obs`; DNA-free `→ f_g≈0`);
  precision monotone in `C` (more count ⇒ sharper); `Σg=0 ⇒ imprecise-but-low`; sums to a proper factor.

### Phase 2 — thread it into the sweep as a first-class λ arm
Following the role-separation discipline (`bp_solver.py` keeps `gdna_prior` and `enrichment_prior` as distinct
kwargs — do the same):
- `simplex_logodds`: add a `lam_logprior=None` per-node `(m, K)` argument to `_solve_nodes_logodds_all` →
  `_solve_nodes_logodds` / `_solve_ambig_logodds` → `_local_loglik_logodds`, sub-indexed per class
  (`_s(lam_logprior, ss/bidx)`) and **regridded** for the fine single-strand grid (`_regrid_global`, as
  `global_logprior` already is). **Add** it to ψ alongside the arms (do NOT fold into `global_logprior` — that
  slot *replaces* the Jeffreys reference; folding would drop the `f_g→1` bound):
  - single-strand: `psi = psi + _gdna_arm(...) + _rna_arm(lam) + lam_logprior`   (at the `+267` add site)
  - AMBIG: `psi += (_gdna_arm(...) + _rna_arm(lam) + lam_logprior)[:, :, None]`   (θ-broadcast, at `+465`)
- `bp_solver.node_sweep`: add an `intron_prior=None` kwarg; before the pass, evaluate the factor **once**
  (anchored, like every other prior) into a `(n_nodes, K)` array `lam_lp` that is zero on non-intron nodes and
  the factor on intron REGION nodes (`_node_region_type == INTRON`); pass `lam_logprior=lam_lp` through
  `_local_solve`. It participates in both the local solve *and* the message relay (so the confident intron gDNA
  belief **propagates**).
- `calibrate`: build the `IntronBackground` once (§Phase 1) and the factor; pass as `intron_prior` into
  `node_sweep` at the **pass-0** call (`_sweep`). Gate behind `config.intron_factory` (default off for a clean
  byte-identical A/B, matching how `_TAU_PRECISION` / `_COMPOSITION_MODE` landed).

### Phase 3 — reconcile with Role B (owner decision §9.3)
Once validated: make the Phase-2 `_fit_gdna_hyperprior` refit either skip introns (factory owns them) or re-fit
the enrichment on the factory-cleaned substrate. A/B the refit on/off with the factory on.

### Config
`intron_factory: bool = False`, `intron_factory_pi_bg: {"poisson","negbin","npmle"} = "negbin"` (§8). No numeric
knobs; `φ` is fitted. Reuse `background_*` config for the pool substrate.

---

## 11. Validation plan (pass-0 vs oracle, cached ambig 10mb — no full pipeline)

Reuse `scratchpad/validate_tau.py` / `validate_composition.py` harnesses (pass-0 off `_selfsolve_cache`, factory
OFF vs ON, `OMP_NUM_THREADS=1`). Metrics, per condition:
- **Intron `f_g` vs oracle** (mass-weighted `|Δf_g|`): must drop sharply on unstranded conditions (the target),
  toward the §7 point-estimate agreement, with honest precision.
- **Confidently-wrong mass** (`var<1 & |Δf_g|>0.3`): the NPMLE-corrupting quantity — must fall at introns
  (both the gdna_none false-positive and the gdna300 under-call).
- **Propagation:** exon/boundary `f_g` adjacent to informed introns should improve (gDNA flows genomically).
- **Stranded controls** (`ss0.99`): must **not** regress (the factory agrees with the strand there; it only adds
  information where strand is silent).
- **Precision honesty:** intron `var_gdna` should be small where §7 says confident (high `Σg`, `ρ_obs` far from
  `ρ_bg`) and wide where `Σg` small / `ρ_obs ≈ ρ_bg` (the overlap).
- Always include **ample single-strand nodes** (per the standing AMBIG-benchmark rule); the 573-intron mix
  (510 G2 / 63 G3) satisfies this.

Gate: intron `|Δf_g|` improves on unstranded, zero-gDNA false-positive mass → ~0, stranded controls flat, no
new confidently-wrong mass. Then hand to the owner for the Role-B reconciliation (Phase 3) and real-data
`include_introns=True` tuning.

---

## Appendix — the one-line intuition

`f_g^{intron} = ρ_bg / ρ_obs`, with a Poisson×background precision — *because an intron is off-target gDNA at the
same rate as intergenic, plus whatever nascent it happens to transcribe.* Everything above is the honest,
strand-free, propagating, count-over-length realization of that one line.

---

## 12. Extending to ALL nodes — the factor form is ENRICHMENT-dependent (owner Q, 2026-07-20)

**The question.** Can the factory be applied to exons / exon-intron boundaries, not just introns? The concern is
hybrid capture: it enriches gDNA *and* RNA above background, so at an enriched exon (ρ_obs=30, ρ_bg=0.01) we can
confidently peel only the depleted floor (~0.01) and the residual 29.99 is *unknown* (gDNA or RNA). Extending
must stay **honest** — never overconfident about an enriched node's composition.

**Why the current (two-sided) factor is UNSAFE on enriched exons.** The factory factor is
`log NegBinom(f_g·C; μ=ρ_bg·E_g, α_eff)` — a **two-sided** density with a *mode* at `f_g=ρ_bg/ρ_obs` that
penalizes `f_g` on BOTH sides. That is correct for an **off-target** node, where `ρ_g ≈ ρ_bg` (the excess is
confidently the other component). But an **on-target** exon has `ρ_g = ρ_bg·e`, `e≥1` *unknown* → `ρ_g ≥ ρ_bg`
with **no upper bound the background can see** (§8: the gDNA enrichment `e` is confounded with the RNA content in
the total density — unrecoverable). Applied to an enriched exon the two-sided factor peaks at `f_g≈0.0026` and
its NegBinom right tail is a wall (measured: **−15067 nats** at f_g=0.5 vs the peak) → it **confidently CRUSHES**
the exon to ~0. That is the §8 false-negative, and *overconfident*. **Measured through the full sweep** (exon
mass-wt |Δf_g| vs oracle, gdna300 capON): introns-only **0.272** → two-sided-on-exons **0.674** (a large
regression, the crush). Confirmed.

**The honest form for an enriched node — a ONE-SIDED lower CDF floor.** Replace the pmf with the background CDF:
`ψ_exon(f_g) = log P(g_bg ≤ f_g·C)` (the NegBinom CDF at the candidate gDNA count). This is **flat above the
floor** (measured: 0.000 nats at f_g=0.5 — *no opinion* about the 29.99 residual) and penalizes only **below**
the floor (`ρ_g` cannot be under the demonstrable background). It is the §8 one-sided log-floor exactly: it can
*lift* `f_g` off a wrong crush up to the background, but never pushes it toward a confident composition. Because
its curvature above the floor is **zero**, it *cannot* make an enriched exon overconfident. It also adds
essentially nothing to a heavily-enriched exon (the tiny floor is trivially cleared) — which is the correct,
honest behavior: **the background genuinely cannot resolve an enriched exon's composition** (§8). Measured: capON
exon 0.272 → **0.265** (safe, no crush).

**The twist — the correct form FLIPS with capture, so it cannot be chosen structurally.** Under capture-**OFF**,
exons are *also off-target* (`e≈1`, `ρ_g≈ρ_bg`), so the TWO-SIDED pin is right for them — and hugely so
(measured capOFF exon **0.249 → 0.007**). The one-sided floor, conversely, *under-constrains* those off-target
exons and **regresses** them (0.249 → 0.346). So:

| node / regime | correct form | two-sided on exon | one-sided floor on exon |
|---|---|---|---|
| intron (always off-target) | **two-sided pin** | ✓ (validated) | under-constrains |
| exon, capture-OFF (off-target) | **two-sided pin** | ✓ 0.249→**0.007** | ✗ 0.249→0.346 |
| exon, capture-ON (enriched) | **one-sided floor** | ✗ 0.272→**0.674** (crush) | ✓ 0.272→0.265 |

**Neither fixed form is universally safe on exons.** The correct form depends on the node's *enrichment* —
two-sided where `ρ_g≈ρ_bg` (off-target), one-sided where `ρ_g≥ρ_bg` (enriched) — and the per-node gDNA enrichment
is precisely the §8-confounded quantity (unreadable from total density, since a high total may be expression, not
gDNA capture). This is why the factory is **intron-scoped**: an intron is off-target *by capture design* (probes
target exons), in EVERY regime — so the two-sided pin is always valid there (no-regression across all 32
scenarios, capOFF/ON/verystrong).

**Recommendation.** Keep the *direct* factory intron- (and intergenic-) scoped. Exon gDNA should arrive by the
three §8-sanctioned routes, not a direct two-sided exon factory:
1. **Propagation** — the now-confident intron gDNA flows genomically to the flanking exon/boundary via the
   message system (gDNA is smooth). The factory supplies the *accurate source beliefs* the message model needed
   (cf. [`density_cliff_and_mature_absorption.md`] / [[composition_mode_regresses_post_tau]]); the intron→exon
   relay across the capture cliff is where that message work re-engages, now with honest sources.
2. **The one-sided background floor** in the Phase-2 hyperprior (already the §8 pinned-background component) —
   the confident lower bound on exon gDNA, safe under capture.
3. **The on-target deconvolution** (strand + spliced) + the iterative global abundance — the enriched exon gDNA
   *level*, which the background floor cannot supply.

**Optional future direct-exon term** (only if the above under-serve): an **enrichment-adaptive band**
`log P(ρ_bg ≤ ρ_g ≤ ρ_bg·e_node)` — two-sided when the node is confidently off-target (`e_node≈1`), widening to
the one-sided floor as `e_node` grows. This needs a per-node/library enrichment classifier (library-level "is
captured?" is readable from the total-density NPMLE's enriched mode; the per-position off-target vs
highly-expressed call is the hard, confounded part). Defer until a scenario demands it; the intron anchor +
propagation covers the validated need.
