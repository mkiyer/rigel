# THE θ COORDINATE — prototype, derivation, and measurements (2026-08-24)

    ⚠ **A DEV DOC.** ⛔ VOCABULARY (owner, 2026-08-24): "tilt" in this codebase means the
    RNA+ vs RNA− nuisance axis and NOTHING else. The reference prior's mean is called the
    **location term** here; an earlier draft of this doc (and the external review) called it
    "the location tilt", which collides with the solver's real tilt DOF and is retired.

     The session's answer to `NEXT_SESSION` ⑥ steps 2–3: the solve coordinate
    re-derived as θ_g = arcsin(√f_g), prototyped end-to-end against the shipped solver, and the
    density channel extended to the location term's population as the replacement voice. Nothing here is
    shipped; no default moved; every number is re-derivable from the recipe in §7.

## 1. THE DERIVATION — one identity, and it is exact

The shipped 1-D ψ is, over the λ = logit(f_g) grid:

    psi_lambda = loglik(f_g) + ½·log f_g + ½·log(1−f_g) + location

Change variables λ → θ = arcsin(√f_g). The density picks up |dλ/dθ| = 2/√(f(1−f)), and the two
Jeffreys halves cancel it **exactly**:

    exp(psi_lambda) · |dλ/dθ|  =  2 · exp( loglik + location )

⭐⭐⭐ **So a UNIFORM grid in θ carries Beta(½,½) as its own measure.** The two reference arms ARE
the Jacobian: under θ they are not written at all, and (with the refuted location term dropped)
ψ_θ is the bare likelihood — every remaining term is a likelihood or a message.

* **Fitted arms transfer verbatim**: ψ_θ = loglik + logP_g + logP_r. The log-rate-density
  conversions that cancelled `log σ'(λ)` once per group on the λ axis are absorbed by the same
  identity, so the "no Jacobian is written" property survives with a fitted landscape.
* **The tilt keeps its shipped arcsine grid**, so the 2-D AMBIG measure is flat in
  (θ_g, θ_tilt) and the Berger–Bernardo reference vanishes identically on BOTH axes — the module's
  own §3 argument, applied to the parameter of interest. ONE uniform treatment for both DOFs.
* **The vertices f_g = 0 and 1 sit AT θ = 0 and π/2** — on the grid, exactly. There is no window
  constant `L`, no state-space clip, and nothing improper anywhere: the domain is compact and the
  likelihood bounded, so ψ_θ is proper with no stabilizer term. The Haldane trap
  (`simplex_logodds` fact 1) is not dodged — it is dissolved: there is no "omitted term" state.
* **Endpoint policy** (the one new decision): the two exact-vertex nodes stand for their
  half-bins, so `log f` at θ = 0 (and `log(1−f)` at π/2) is evaluated at the half-bin's
  representative point `f(Δθ/4)` — grid-derived, count-independent, monotone; every interior node
  exact. ⛔ A count floor `log(1/(n+1))` was tried first and REJECTED: at n = 5,000 it sits at
  f = 2·10⁻⁴, COARSER than the λ window it replaces.

## 2. WHAT THE COORDINATE DOES **NOT** CHANGE — the honest half, stated first

⭐⭐⭐ **The continuum posterior is IDENTICAL — the identity above says so, and the prototype
measures it** (§4 P1: max|Δf_g| = 9.3·10⁻⁴ over random interior slots at K = 2001, and ~0.1–2 %
tracking of the `w = 0` arm across five full panel conditions). So:

* The κ = ½ pathology (a written location term never overturned at any depth) is
  **coordinate-independent** — with zero information, the answer is the prior in ANY coordinate.
  `NEXT_SESSION` ①'s "the collapse is the COORDINATE" is true of the leverage bookkeeping of a
  WRITTEN location term on the λ axis; the substantive act is REMOVING the term, and §H.10 ② already
  priced that removal on the panel (worse on every contaminated rung until a real channel
  replaces it — which is §5's job).
* The near-vertex closure rate is the posterior's honest 1/N, both coordinates (§4 P2). The
  `_rna_arm` "3.107-nat repulsion" framing measures Beta(½,½) against improper flat-λ; the θ
  solve keeps the same proper measure and therefore the same honest offset.

**What the coordinate DOES buy** — all structural, none a panel score:

1. The reference stops being a WRITTEN TERM. No arms to keep balanced, no `_location_term`
   socket for an assertion to creep back into, no `m = ½` float-noise special case, no
   "omitting a term is Haldane" trap for the next editor.
2. **Exact vertices and no `L`.** The measured L-pathology (`simplex_logodds` header: the fitted
   landscape pushing against σ(−L) costs 3.2× at `g00`, saturating only by L = 40) is a
   coordinate artefact with a coordinate-level fix: the fitted prior can point at any density and
   θ can represent it. The `required_logodds_window` machinery, the bracket-rebuild memo, and the
   `_regrid_global` window maps all become dead concepts.
3. **Honest near-vertex message precisions.** On the λ grid a vertex-concentrated posterior has
   Var(log f) computed over log σ(λ) ∈ [−10, −9.99…] — tiny, so the slot SENDS huge precision.
   The θ grid's endpoint carries the half-bin representative (−24 at K = 2001), so the same
   posterior reports honestly wide log-fraction variance. Same mechanism family as the fan-out's
   over-confidence (§H.5 ②); not separately priced yet.
4. The 1-DOF and AMBIG paths and both tilt/interest axes get ONE uniform geometry.

## 3. THE PROTOTYPE

A drop-in `solve_regions_theta_all` with `_solve_regions_logodds_all`'s exact signature and
routing (1-D single-strand / 2-D AMBIG, same `RegionDeconv` contract, same message terms, same
`_mixture_strand_loglik` / `_rna_residual` / `_compose` / median-readout machinery — the median
interpolates on the UNIFORM θ lattice and maps through sin², the same transform-invariance
argument as the shipped λ interpolation). Per-row arrays that arrive pre-evaluated on the λ grid
(fitted arms, the intron-factory factor) are re-evaluated at λ(θ) by shared-grid linear
interpolation. Installed for measurement by rebinding `_solve_regions_logodds_all` in `sweep`,
`region_init` and `region_geometry` (all three bind the name at import — the recorded harness
trap), asserted taken.

## 4. PROPERTY FALSIFICATIONS — all pass, and both perturbations fire

| | property | result |
|---|---|---|
| P1 | interior equivalence to shipped-λ (location=None), K = 2001 | max\|Δf_g\| = 9.3·10⁻⁴; injecting a stray +½ log f moves it to 6.8·10⁻² (the test can fail) |
| P2 | truth exactly at a vertex, κ = 0.99 | both close as ~1/N; θ ≤ λ at every depth (8.3·10⁻⁴ vs 1.05·10⁻³ at N = 10⁵) |
| P3 | invariance | θ: K-drift 2.3·10⁻⁵ (256 → 2001), no L exists. λ: an off-window message (target log f = −13.8) moves the answer **3.81 nats** between L = 10 and L = 20; θ lands within 0.49 nats of the L = 20 answer with no window at all |
| P4 | κ = ½, no channels | exactly 0.5000 at every depth — proper, centred, no pole, no term written |
| P5 | AMBIG 2-D equivalence (BB vanishes on both axes) | max\|Δf_g\| = 2.6·10⁻³, max\|Δw_pos\| = 7.9·10⁻⁴ vs the shipped (λ, θ_t) path |
| P6 | the review's §4.1 leverage table, no location concept | κ = 0.99/0.90/0.75 rows track the likelihood (0.090/0.133/0.202 at N = 10 → 0.003/0.003/0.006 at N = 10⁴); κ = ½ row is 0.5 at every depth — **the "never overturned at 0.747" row cannot exist because there is no term to hold it there** |

## 4b. SPARSE DATA — the owner's 0/1/2-fragment question (asked and answered 2026-08-24)

Single-strand slot, + strand live; a + read is RNA-consistent, a − read at high κ is gDNA
evidence. `f_g` readouts (θ vs shipped-λ-neutral vs shipped-λ with the 0.75 location term):

| κ | case | θ | λ neutral | λ m=0.75 |
|---|---|---|---|---|
| 0.99 | 1 read, + | 0.358 | 0.360 | 0.597 |
| 0.99 | 1 read, − | 0.823 | 0.819 | 0.900 |
| 0.99 | 2 reads, ++ | 0.265 | 0.268 | 0.445 |
| 0.99 | 2 reads, +− | 0.728 | 0.725 | 0.852 |
| 0.99 | 2 reads, −− | 0.912 | 0.910 | 0.940 |
| 0.75 | 1 read, + | 0.404 | 0.405 | 0.666 |
| 0.75 | 1 read, − | 0.655 | 0.652 | 0.831 |
| 0.50 | any counts | 0.500 | 0.500 | **0.747** |

(An EMPTY slot is zeroed by the active-mask before any solve — the ½ honest-ignorance readout is
the κ = ½ row, where the slot has mass but the strand channel carries nothing.)

Readings: the moves are exactly the likelihood's own odds — at κ = 0.99 one + read is ~2:1 for
RNA (an RNA read is +, but half of gDNA reads are + too) ⇒ 0.36, while one − read is ~50:1 for
gDNA ⇒ 0.82; the asymmetry is the evidence's, not a prior's. Two mixed reads resolve toward the
stronger read. Nothing jumps to a vertex; the posterior VARIANCE at these depths is large, so the
precision such a slot would SEND onward is honestly weak — the sparse regime feeds the
weak-message requirement rather than fighting it. The λ-neutral column agrees to ±0.003 (the
measure identity); the 0.75 location column shows what the old reference did to every one of
these reads, including holding κ = ½ at 0.747 at any depth.

## 5. ON THE PANEL — three arms at `calib_refit_iters = 0`, claimed slots, misplaced gDNA fragments

`ship` = shipped λ + location term; `w0` = shipped λ, location term neutralised; `theta` = the
prototype (no location term, no reference terms). Populations B/E as in `pass0_claimed_ab.py`; never pooled.

| condition | pol | pop | ship | w0 | theta |
|---|---|---|---|---|---|
| `g50` ss.99 ON | silent | B | **38,650** | 43,533 | 42,701 |
| `g50` ss.99 ON | silent | E | 17,359 | 17,359 | **17,145** |
| `g50` ss.99 ON | fanout | B | **41,416** | 54,971 | 54,552 |
| `g50` ss.99 ON | fanout | E | **19,106** | 20,618 | 20,535 |
| `g98` ss.99 ON | fanout | B | **54,571** | 82,888 | 81,840 |
| `g98` ss.99 ON | fanout | E | **23,910** | 28,543 | 28,318 |
| `g50` ss.50 OFF | silent | B | **68,577** | 85,742 | 85,742 |
| `g50` ss.50 OFF | silent | E | 268,498 | 268,498 | 268,508 |
| `g50` ss.50 ON | fanout | B | **187,485** | 448,749 | 447,572 |
| `g50` ss.50 ON | fanout | E | **270,741** | 319,825 | 319,172 |

Readings: θ tracks `w0` within 0.1–2 % everywhere (the measure identity, at panel scale, through
the full pipeline including fan-out) and edges it on most stranded rows (exact vertices). The gap
to `ship` is the location term's removal — §H.10 ②'s number, reproduced, NOT a property of the coordinate.
The coordinate is a foundation; the location term's only-voice role needs a real channel:

## 6. THE REPLACEMENT VOICE — the density CHANNEL at the location term's population (rec 3, prototyped)

`density_lambda_factor` extended to claimed ss-intron BOUNDARIES with the boundary's OWN frame:
`log NegBinom(f_g·C_cross; ρ_bg·E_g,cross, α_eff)` — crossing count against crossing gDNA
opportunity, ρ_bg pooled off intergenic REGION slots exactly as `fit_intron_background` selects.
The `_build_intron_prior` docstring's frame objection ("a boundary's count is a crossing and its
divisor a different formula") is answered, not overridden: both numerator and divisor are the
boundary's own. Strength scales with counts; zero counts ⇒ wide; a flat row ⇒ τ = 0.

**THE μ-CHECK FIRST — no solver, the direct test of the transported rate on the boundary crossing
frame** (`ρ_bg · E_g,cross` against certified gDNA at claimed DELIVER boundaries):

| condition | ρ_bg | n deliver | Σμ predicted | Σ certified | median ratio |
|---|---|---|---|---|---|
| `g50` ss.50 OFF | 5.13·10⁻² | 8,608 | 95,166 | 95,308 | **1.005** |
| `g50` ss.99 ON | 2.12·10⁻³ | 7,898 | 4,254 | 496,707 | **0.026** |
| `g98` ss.99 ON | 4.19·10⁻³ | 12,939 | 13,745 | 1,770,466 | **0.020** |
| `g00` (both) | ~1·10⁻⁸ | 0 | 0 | 0 | — |

⭐⭐⭐ **Capture-OFF the intergenic rate is THE right rate for a boundary crossing, to half a
percent** — the reciprocal-opportunity frame transfer works across axes. **Capture-ON it is 38–50×
too small** — §H.8 ④b's probe footprint, now measured on this exact frame.

**And the panel follows the μ-check exactly** (population B, `calib_refit_iters = 0`; `ship` and
`w0` from §5's run):

| condition | pol | ship (location) | w0 | w0+bfac | theta+bfac |
|---|---|---|---|---|---|
| `g00` ss.50 OFF | silent | 117,792 | 78,829 | **208** | **152** |
| `g00` ss.99 ON | silent | 3,601 | 2,945 | **91** | **69** |
| `g50` ss.50 OFF | silent | 68,577 | 85,742 | **35,505** | **35,294** |
| `g50` ss.99 ON | silent | **38,650** | 43,533 | 1,082,602 | 1,083,228 |
| `g50` ss.99 ON | fanout | **41,416** | 54,971 | 1,075,903 | 1,076,240 |
| `g98` ss.99 ON | fanout | **54,571** | 82,888 | 2,104,936 | 2,105,427 |
| `g50` ss.50 ON | fanout | **187,485** | 448,749 | 1,085,718 | 1,085,918 |

*(every arm at `calib_refit_iters = 0`, same process family; the `g00` ship/w0 rows re-measured
here rather than quoted from the review's shipped-refits sweep.)* Exons are untouched by the
factor under silence (E identical to `w0` per condition) and move only through messages under
fanout (small: 20,618 → 20,047 at `g50` ss.99 ON).

**READ IT AS TWO VERDICTS, one per capture state:**

* ⭐⭐⭐ **Capture-OFF, the channel is the location term's replacement and BETTER THAN it**: the zero
  controls go to ~0 (208 vs the shipped location term's 117,792 — 566×, and 91 vs 3,601 capture-ON at
  zero gDNA; requirement 3 met by construction and by measurement), and the unstranded
  contaminated stratum — where the zero-refit partition showed the location term IS the whole answer —
  improves 1.9× over the shipped location term and 2.4× over Jeffreys-only, with BOTH obligations
  improving (deliver 47,653 → 15,311, refute 38,088 → 20,194).
* ⛔⛔ **Capture-ON it must not run with the intergenic rate**: the likelihood form asserts the
  40×-wrong rate WITH count-scaled confidence, so it is far worse than the weak location form ever
  was (~26× over ship at `g50` ss.99). The gate is the enrichment detector (the ④a boolean);
  the fix is the on-rate anchor (§H.8 ④b), unchanged and now priced on this frame: 38–50×.
* θ and λ agree through the whole experiment (0.1–2 %) — the coordinate neither causes nor cures
  any of this; it just removes the term the location assertion lived in.

## 7. REBUILD RECIPE

Prototype + experiments live in this session's scratchpad (ephemeral):
`theta_solver.py` (the drop-in), `theta_properties.py` (P1–P6, self-failing), `theta_panel.py`
(the §5 table), `boundary_factor_panel.py` (§6; wraps `solve_chain` in BOTH `sweep` and
`calibrate` — the name is bound into each), `who_decides_refit0.py` (the partition re-run now
recorded in `REFERENCE_PRIOR_REVIEW.md` §4.2). Each is ~200 lines over the shipped machinery
(`pass0_claimed_ab.claimed_masks/score`, cached scans, certified `slot_truth`); nothing needs a
re-scan. The panel arms run at `calib_refit_iters = 0` off `oracle_cache/<cond>/_main`.

## 8. WHAT A SHIPPED VERSION STILL NEEDS (not done, deliberately)

* The float32 AMBIG cube path (prototype is f64) and the `_SOLVE_BLOCK_BYTES` tiling discipline.
* The fitted-arm evaluation NATIVELY on θ (the prototype interpolates from the λ grid; a real
  build evaluates `landscape.logprior` at `log ρ = log f(θ) + log M − log E` directly, which is
  where pathology-②'s L-win should materialise — untested at refits = 0).
* The reference-term test files (`test_structural_reference.py`, the location-term gates,
  `L`-invariance gates) become obsolete-by-construction and need replacement properties
  (K-invariance, endpoint policy, the P1 identity as a regression gate).
* An owner ruling: this deletes `_location_term`, `structural_reference_location`,
  `measured_reference_location`, both config flags, and the `L` window from the solve — the
  "converge and delete" is large and the panel cost of the location term's removal stands until §6's
  channel (or the landscape prior at refits > 0) covers it.
