# W12 — the θ quadrature at zero gDNA: the derivation (2026-09-13)

A working note, provisional like everything in this directory. The issue is `ISSUES: theta-quadrature-at-zero-gdna`;
the code is `simplex_logodds._tilt_grid`, `_psi` and `_solve_logodds`; the ruling home, once there is one, is
`DESIGN.md` §6b.15 and the derivation's home `EQUATIONS.md`.

## 1. What ψ integrates over θ, exactly

At an AMBIG slot ψ is evaluated on the `(λ, θ)` cube and read out on its θ-marginal, `M(λ) = Σ_k exp ψ(λ, θ_k)`
over the uniform lattice `θ_k = −π/2 + k·h`, `h = π/(K_t − 1)`, endpoints included, no weights. Only the strand
term and the delivered cube rows depend on θ; the two arms and the λ-factor rows do not.

The strand term (`_mixture_strand_loglik`) is a Gaussian quasi-likelihood in the aligned-strand rate `p` with the
variance FROZEN at the reference composition, so at fixed λ it is an exact Gaussian in the tilt `τ = sin θ`:

    p(λ, τ) = ½·f_g + κ·f₊ + (1−κ)·f₋ = ½ + a(λ)·τ,          a(λ) = (1 − f_g)(κ − ½)
    L(λ, τ) = −(u₊ − n·p)² / 2V  =  −(τ − τ̂(λ))² / 2σ_τ(λ)²  + const(λ)

    τ̂(λ)  = d / a(λ),        d = u₊/n − ½   (the observed strand contrast)
    σ_τ(λ) = σ_p / |a(λ)|,    σ_p = √V / n    (V the frozen variance; σ_p ∝ n^−½)

Both the centre and the width depend on λ: as `f_g` rises, `a` shrinks, `τ̂` moves outward and `σ_τ` widens.
For a strand-pure slot (`d = κ − ½`) the centre is `τ̂ = 1/(1 − f_g)`: ON the boundary at `f_g = 0`, BEYOND
it at every `f_g > 0`.

The reference measure is the arcsine weight (the Berger–Bernardo tilt conditional `(1 − τ²)^−½`, cancelled by
gridding θ), so the marginal the lattice approximates is

    M(λ) = ∫_{−π/2}^{π/2} e^{ψ(λ, sin θ)} dθ = ∫_{−1}^{1} e^{L(λ,τ)} · c(λ, τ) · (1 − τ²)^{−½} dτ

with `c` the delivered cube row (smooth in τ; 1 when nothing is delivered).

## 2. The lattice is the trapezoid rule in φ — with the wrong endpoint weights

With `φ = θ + π/2 ∈ [0, π]` the integrand `g(φ) = e^{ψ(λ, −cos φ)}` is the restriction of a smooth 2π-periodic
even function, so the periodic trapezoid rule is spectrally accurate on it: the correct rule on `[0, π]` is
`h·[½g(0) + g(h) + … + g(π − h) + ½g(π)]`. The plain sum weights both endpoints by 1. Its error is therefore

    S_lattice − S_trapezoid = ½·[g(0) + g(π)]     (first order in h, not spectral)

and it is not a constant factor across λ: it is large exactly where the integrand peaks ON the endpoint —
every λ with `τ̂(λ) ≥ 1`, i.e. every `f_g > 0` on a strand-pure slot — and negligible where the peak is
interior. That is a λ-dependent inflation of the marginal toward higher `f_g`: a gDNA bias, at every `K_t`, even
where the peak is resolved. The fix is the rule's own weights (½ at `θ = ±π/2`), not a tunable.

## 3. Resolution: why the bias is a λ-dependence of the error, and where it lives

Near the boundary `τ = 1 − δ²/2` (`δ = π/2 − θ`), so the Gaussian in τ becomes, in θ:

* `τ̂ = 1` (zero gDNA, strand-pure): `L = −δ⁴ / 8σ_τ²` — a QUARTIC peak of width `δ₀ ≈ (8σ_τ²)^¼ ∝ n^−¼`.
  At 500k fragments, `κ = 0.99`: `σ_p = 1.4e−4`, `σ_τ = 2.9e−4`, `δ₀ = 0.029 rad` — the recorded 0.03 rad
  against the `K_t = 60` step 0.053.
* `τ̂ = 1 + ε` (beyond): `L = L(1) − ε·δ² / 2σ_τ² − δ⁴/8σ_τ²` — a QUADRATIC peak of width
  `σ_φ = σ_τ / √ε`, `ε ≈ f_g`, NARROWING as `f_g` rises, sitting on the endpoint node.

The lattice sum over an endpoint peak of width `σ_φ` is `≈ g(0)·(1 + Σ_k e^{−(kh)²/2σ_φ²})`, while the
integral is `≈ g(0)·σ_φ·√(π/2)`. The error factor `E(λ) = h·S/M` is 1 when `σ_φ ≫ h` and grows like
`h/(1.25·σ_φ) ∝ h·√f_g / σ_τ` once `σ_φ < h`. So the θ-marginal is inflated MORE at higher `f_g`, the
λ-posterior tilts toward gDNA, and the size of the tilt is set by `h·√f_g/σ_τ` in the band `f_g ~ σ_τ` where
the true posterior still has mass (`e^{−f_g²/2σ_τ²}`): a shift of order `σ_τ` in `f_g`, i.e. `~ σ_p·n/|κ−½|`
fragments per slot — a few fragments at 50k (the recorded ≤ 4, prior-free), tens at 500k. The peak-width
analysis at `f_g = 0` alone cannot see it, because the quartic peak at zero gDNA is the WIDEST point of the
curve; the bias is the derivative of the error along λ.

Two further consequences: a finer lattice cures it only where `h < σ_φ` for every λ that carries mass, which
depends on `n` (it saturates on the ladder at 60 and is latent on a deep library's heavy exons); and any mesh
clustered at the ends in θ moves the nodes toward the quartic peak at `f_g = 0`, not toward the narrower
quadratic peaks at `f_g > 0` where the inflation lives — consistent with the refused Chebyshev mesh breaking a
`g00` row.

## 4. Candidate designs — accuracy independent of node count

The θ-integral at each `(slot, λ)` is one dimension of a Gaussian times a smooth function times the arcsine
weight. Three ways to integrate it whose error does not depend on `K_t`:

* **A. The trapezoid endpoint weights** (½ at `±π/2`). Not the whole answer — it removes the first-order term
  of §2, not the resolution term of §3 — but exact where the peak is resolved, costs nothing, and is the rule's
  own weights. Prototype first, alone, because it is the one mechanism that can be A/B'd in isolation.
* **B. Gauss–Hermite per λ.** Centre `τ̂(λ)`, scale `σ_τ(λ)`, a fixed handful of nodes: integrates the Gaussian
  factor exactly against the smooth remainder (`c(λ,τ)·(1−τ²)^−½`, interpolated from the lattice), with the
  truncation at `±1` handled by the truncated-normal weights. Accuracy is set by the smoothness of the
  remainder, not by `K_t`, and the node count falls from 60 to ~8–16 — the cube's cost lever. The delivered
  cube rows are the one complication: they arrive on the fixed θ lattice and must be read at λ-dependent nodes.
* **C. The analytic endpoint.** Beyond the boundary the marginal has the closed form
  `M ≈ e^{L(1)}·c(λ,1)·√(π σ_τ² / 2(τ̂−1))`, and at the boundary `∝ e^{L(1)}·√σ_τ·Γ(¼)/2^{5/4}`; a composite rule
  uses the lattice where resolved and the closed form where not. Exact in both limits, but a regime switch is a
  threshold, and a threshold is a tunable.

## 5. The plan, and where it stands

1. DERIVE — this note. The two recorded numbers (0.03 rad at 500k; ≤ 4 fragments at 50k) are reproduced by §3.
2. PROTOTYPE outside `src/`. **A is REFUTED (2026-09-13)**: the endpoint weights, patched into every binding of
   `_psi` (the session scratchpad's `w12/theta_endpoint.py`), moved ≤ 0.7 false-gDNA fragments on any `g00`
   row of either panel at `K_t` 60 or 30, and left the recorded failure untouched — ladder `g00 ss.99 ON`
   reads 9,821.3 → 9,821.2 at `K_t` 30 against 194.1 → 194.0 at 60; the test chromosome is identical at 30
   and 60 on every row (shallow slots, every peak resolved). So the first-order endpoint term of §2 is
   negligible at these depths and the resolution term of §3 is the whole mechanism. Next is B (or the
   analytic marginal of §4 C as one smooth formula), which changes the cube's shape and must handle the
   delivered cube rows — a design decision before the prototype.
3. A/B on both panels with both zero controls: `calibration_vs_oracle.py`, the `g00` rows under
   `--set calibration.sweep_n_tilt=30` and at the default 60, patched against unpatched, the metric per stratum;
   the panel; the profiler (the cube is the cost).
4. `src/` only after the A/B, one mechanism at a time; the ruling to `DESIGN.md` §6b.15, the derivation to
   `EQUATIONS.md`, the issue closed with its numbers.


# The design decision and its A/B (2026-09-13, the session after the derivation)

The prototype and the instruments are the session scratchpad's `w12/` (`theta_window.py` the rule patched into every binding of `_psi`; `tilt_census.py` the per-slot census of where the tilt matters; `deep_stress.py` the shared-exon toy; `quadrature_check.py` the standalone check against adaptive quadrature; `oracle_summary.py`). Nothing in `src/` has moved.

## 6. The rule (the design decision, 2026-09-13): the θ nodes follow the strand term's peak (the windowed trapezoid rule)

At fixed λ the strand term is an exact Gaussian in τ (THETA_QUADRATURE.md §1): centre τ̂(λ) = d/a(λ), width
σ_τ(λ) = σ_p/|a(λ)|, a(λ) = (1 − f_g)(κ − ½), d = u₊/n − ½, σ_p = √V/n with V the frozen variance. In
θ = arcsin τ the integrand g(θ) = exp ψ(λ, sin θ) is smooth on [−π/2, π/2] and EVEN about both ends
(sin(π − θ) = sin θ), so it is the restriction of a smooth 2π-periodic function.

The rule, per (slot, λ):

    τ_m = clip(τ̂, −1, 1)                        the maximum of the strand term ON the domain
    ρ   = √((τ_m − τ̂)² + 2 σ_τ² T)               the half-width where the term is within T nats of that maximum
    [τ_lo, τ_hi] = [τ̂ − ρ, τ̂ + ρ] ∩ [−1, 1]      one closed form for the interior, boundary and beyond regimes
    θ_k = θ_lo + k h,  h = (θ_hi − θ_lo)/(K_t − 1),  k = 0 … K_t − 1,  θ_lo/hi = arcsin τ_lo/hi
    M(λ) ≈ h · [½ g(θ_0) + g(θ_1) + … + g(θ_{K_t−2}) + ½ g(θ_{K_t−1})]

and in ψ the weights enter as log h(slot, λ) (+ log ½ at the two end nodes) added before the one exp.

Why it is exact in every regime:
* interior (|τ̂| < 1 − ρ): the window holds the peak to e^{−T}; the integrand vanishes at both window ends
  to e^{−T} of its peak, so the trapezoid rule is the periodic trapezoid rule on a smooth function and
  converges spectrally: error ≈ 2 e^{−2π² (σ_θ/h)²}, σ_θ = σ_τ / cos θ*;
* at or beyond a domain end (τ_hi = 1, say): the integrand is even about π/2, so the rule with weight ½ at
  θ = π/2 IS the periodic trapezoid rule on the reflected window [θ_lo, π − θ_lo] halved — spectral again,
  on the quartic peak (τ̂ = 1) and the one-sided quadratic peak (τ̂ > 1) alike; no case split, the
  endpoint weight is the ½ the trapezoid rule always carries (immaterial at an interior end, exact at a
  domain end);
* no strand information (κ = ½, n = 0, or ρ ≥ 1 either side): the window is the whole domain and the rule
  is the fixed lattice with trapezoid weights.

The delivered cube rows (the RNA level lanes) arrive on the fixed θ lattice and are read at the nodes by
linear interpolation — the same resolution they had, no less.

## 7. The two constants are derived, not tuned

* T = −log ε₆₄ ≈ 36.04 nats: the strand term's mass outside the window is erfc(√T) ≈ e^{−T}/√(πT) < ε₆₄ of
  the peak's integral — the truncation is below double precision.
* K_t: the window is 2√(2T) σ_θ wide in the interior, and the trapezoid error is ≤ e^{−T} when
  2π²(σ_θ/h)² ≥ T, i.e. h ≤ σ_θ π√2/√T, so K_t − 1 ≥ 2√(2T)·√T/(π√2) = 2T/π ≈ 23: **K_t = 24**. At the
  boundary the window is only T^{¼} ≈ 2.4 widths per side, so 24 nodes over-resolve it. Measured on the
  standalone check (n 500 … 500k, f_g ∈ {0, 0.3}, τ ∈ {0, 0.5, 0.9, 0.99, 1}): the λ-shape error of log M is
  ≤ 2.1e−6 nats at 24 nodes and ≤ 1e−9 at 60 (the reference's own tolerance); 12 nodes reach 0.07 nats and
  are refused. The shipped lattice at 60 reads 0.01–0.5 nats at n = 500, 8–12 at 50k, 90–130 at 500k;
  2,400 lattice nodes still read 0.1–0.4 at 500k.

## 8. What the lattice was doing (the mechanism, corrected)

The peak's width in θ is σ_θ ≈ σ_τ/cos θ* ≈ 1/(2√n |κ−½| (1−f_g) cos θ*) — at 5k fragments 0.015 rad, at
50k 0.005, at 500k 0.0015, against the 60-lattice's step 0.053. Where σ_θ < h the lattice sum is a COMB: at
each λ it is e^{−(distance of τ̂(λ) to the nearest node)²/2σ_τ²}, a factor between 1 and e^{−(h/2)²/2σ_τ²}
(e^{−100} at 50k), chosen by where τ̂(λ) = d/a(λ) happens to fall as it drifts with λ. The λ posterior of a
deep interior-tilt AMBIG slot is therefore a set of spikes at arbitrary λ. The recorded K_t 30 failure is
ONE such slot: ladder `g00 ss.99 ON`, slot 37345, n = 25,242, τ_true = −0.43, σ_θ = 0.006 rad — 9,627 false
gDNA fragments at 30 nodes, 0.0 at 60, because at 60 a node happens to sit on the peak at the right λ. So
"60 is the smallest count that holds the control" was a coin toss on one slot, not a resolution; and the
strand-purity story (the quartic peak) was not the mechanism of the recorded number.

## 9. Where the tilt matters — the census (`tilt_census.py`, both panels, per stratum)

The population: AMBIG slots (both strands admissible) whose certified RNA is on both strands with the minor
strand ≥ 5 % of it. Ladder: ~1,400–1,600 slots per condition (500–950 regions), 170–250k RNA fragments,
median depth ~35 fragments capture-OFF and 300–900 capture-ON; the predicted peak width exceeds the
half-step on 97–100 % of them, so the ladder is resolved at 60 (and at 30 but for the one slot above).
Test chromosome: 12 such slots at ~700 fragments.

Ladder, K_t 30 → 60 → 120 (the shipped lattice), the both-strand band with a 20–50 % minor strand, stranded rows:
the gDNA error is flat in K_t but for the one slot above (g00 ss.99 ON: 9,637 → 9.5 → 9.4); the TILT error is
NOT converged at 60 — g00 ss.99 OFF 1,421 → 995 → 719 RNA fragments on the wrong strand, g00 ss.99 ON 3,490 →
1,030 → 519, g05 ss.99 OFF 1,797 → 1,203 → 964 — the read-out the owner's population is judged on. The
windowed rule at 24 nodes reads what the lattice reaches only at 120: 718 / 455 / 990 on those three rows, gDNA
error identical to 0.1 fragment; on the strand-pure AMBIG band (RNA on one strand) the exact marginal reads
the gDNA error 0.3–2 % higher (g50 ss.99 ON 38,013 → 38,780) and the tilt error 7 % higher — the exact
integral realises the beyond-boundary tail factor and the posterior-MEAN tilt sits inside the boundary by
the peak's width, where the lattice's endpoint node flattered both; that is the +0.4 % on stranded ON.

The 0.8.0 metric (`calibration_vs_oracle.py`, Σ|Δ gDNA| both axes, g00 excluded), the ladder, lattice 60 (byte-
identical to `arms/2026-09-13_preport/oracle_ladder_default.json`) → window 60 → window 24: stranded OFF
249,666 → 249,703 → 249,702; stranded ON 469,164 → 471,172 → 471,263 (+0.4 %, at g50/g98 ss.99 ON); unstranded
OFF 308,530 → 308,530 → 308,529; deferred 3,622,223 → 3,619,942 → 3,622,383; the four g00 rows identical to
0.1 fragment (405.9 → 405.4 the largest move). The test chromosome: every stratum within ±3 fragments, g00
rows identical. The ladder is shallow where the tilt matters, so the rule is a numeric near-no-op there.

The deep stress (`deep_stress.py`: two spliced genes on opposite strands sharing one 10 kb exon, the shared
exon AMBIG with a junction certifying each strand at one edge, donors g00 and g50 ss.99 OFF, the donor's
priors injected, refits on), the shared exon's false gDNA in fragments, lattice 60 → window 24:

    g00  50k  minor 0.02 / 0.20 / 0.50 :     22 → 3       997 → 56      219 → 44
    g00 500k  minor 0.02 / 0.20 / 0.50 :  1,428 → 261   9,907 → 25   10,994 → 632
    g50  50k  minor 0.02 / 0.20 / 0.50 :    −18 → −19     −12 → −19     −31 → −38
    g50 500k  minor 0.02 / 0.20 / 0.50 :    465 → 29   11,448 → 11  102,076 → 19   (truth 490–498 gDNA)

and the tilt read-out's error (RNA fragments on the wrong strand) falls 6–70× at 500k (g50 balanced: 5,316 →
76). The lattice at 240 nodes roughly matches the window at 500k on gDNA but not on the tilt (1,153 vs 76).
The strand-pure rows (minor 0) are the same under every rule (the boundary peak was never the mechanism).

Node count against the rows' lattice, separated (the `window24:240` arm — 24 nodes, rows on 240): identical to
`window60:240` on every row, so the NODE COUNT is converged at 24; `window 24` differs from `window 60`
only because the config's K_t also sets the lattice the RNA level lanes build their cube rows on, and rows
at 24 differ from rows at 60 ≈ rows at 240. So the rows need their lattice (converged at 60) and the nodes
need 24: today's one knob has two jobs.

## 10. The other finding: the strand marginal's volume factor (NOT the quadrature; a new issue)

The exact θ-marginal of the strand term at an interior tilt is ∝ σ_τ(λ) = σ_p/((1 − f_g)|κ − ½|): the band
of tilts consistent with the data widens as f_g → 1, so the marginal likelihood of f_g carries a factor
1/(1 − f_g) up to where the band fills the domain (1 − f_g ≈ σ_p/|κ−½| ≈ 1/√n). Against the Beta(½,½)
reference the own solve of a balanced pure-RNA slot reads f_g ≈ 1 − 1/√n: 0.960 at n = 500, 0.9955 at 50k,
0.9985 at 500k (the standalone reference, no prior, no messages). It is the Occam factor of the nuisance
tilt — "all gDNA explains balanced strands with no parameter; all RNA needs τ tuned to ±σ_τ" — and it is
Bayes-correct under an f_g-independent arcsine measure on τ, which is what ships. On a real chain the RNA
level lanes and the landscape prior pin such a slot (ladder g00: 25k-fragment both-strand slots read 0.0);
on the mono-exon toy nothing does and a balanced 50k exon reads 0.99 at zero gDNA (prior-free or not).
The spliced toy (junctions on both strands) pins the exon: g00 500k balanced reads 632 false fragments under the exact rule, not 415k; the mono-exon toy with no junction reads 0.997 — the cap from above is the RNA level lanes' delivery, and it is what a real chain supplies. Filed as `ISSUES: strand-marginal-volume-factor` with these numbers; the
candidate is a tilt measure uniform in the observable p over the reachable band (|a(λ)| dτ), which cancels
the factor exactly and leaves the strand term flat in f_g at an interior tilt, as EQUATIONS §5.2 already
says it should be. Its own thread: derive → A/B; it changes the zero-control behaviour at strand-pure slots.

## 11. The tilt measure — the derivation (`ISSUES: strand-marginal-volume-factor`, 2026-09-13)

**The factor is intrinsic to marginalising under any proper prior on the tilt.** With `p = ½ + a(λ)τ` the
strand term at fixed λ is a Gaussian in τ of width `σ_τ = σ_p/|a|`. For ANY normalised prior density `π(τ)`
that does not depend on λ, the marginal is `∫ L π dτ ≈ L_max · π(τ̂) · √(2π) σ_p / |a(λ)|` wherever the peak is
interior — the `1/|a| ∝ 1/(1 − f_g)` volume factor. "Uniform in p over the reachable band", normalised per λ,
is `dp/(2|a|) = dτ/2`: uniform in τ, the factor intact. Normalisation cannot remove it; only a measure whose
mass scales with `|a|` can.

**Why: the strand channel identifies one number.** `u₊ ~ Binomial(n, p)` observes `p` alone; `(f_g, τ)` are not
jointly identified by it (`I_ff·I_ττ − I_fτ² = 0`: the information matrix has rank one), so the marginal
likelihood of `f_g` at an interior tilt is entirely the prior's doing. Under a proper prior it is the Occam
factor — "gDNA explains balanced strands with no parameter" — of order `√n`. The tool's design says the
opposite (`EQUATIONS.md` §5.2: strand measures the tilt and reaches gDNA only through the cap), i.e. the
marginal should be FLAT in `f_g` at an interior tilt. The measure that does that is the unnormalised Jeffreys
conditional of the observed channel, `√I_ττ = √n |a(λ)| / √(p(1−p))`: `∫ L √I_ττ dτ ≈ L_max √(2π) σ_p /
√(p̂(1−p̂))`, constant in λ. Its λ-dependence is `|a(λ)| ∝ (1 − f_g)`, so on the AMBIG cube the whole change is
ONE added term, `+ log(1 − f_g)`, and the `(p(1−p))^{−½}` shape is a λ-constant at the peak.

**Equivalently: the AMBIG reference is Beta(½, 3/2).** `∫ π(f_g, τ) dτ ∝ (1 − f_g) · Beta(½,½)(f_g)`: the
tilt's one degree of freedom, counted in the reference measure as Jeffreys counts every parameter's volume.
A single-strand slot has no tilt to count and keeps Beta(½,½); the two classes' references differ by the
parameter one of them has. This is the honest name of the change and how it should be written in ψ: the RNA
arm's exponent is `½` on a single-strand slot and `3/2` on an AMBIG one.

**The three regimes under the new measure.**
* Interior (`|τ̂| < 1`): `M' = |a|·M` is flat in λ up to the cap and the delivered rows — the design's intent.
* Beyond the boundary (strand-pure, `τ̂ = 1/(1−f_g)`): `M ∝ e^{L(1)} σ_τ /√(τ̂−1)` becomes `M' ∝ e^{L(1)} σ_p
  √(|a|/(d−a)) ∝ √((1−f_g)/f_g)`: the `1/√f_g` tail toward zero gDNA SURVIVES (it is the data's — a pure
  strand split is incompatible with `f_g > 0` beyond `σ_p`), the change is the factor `√(1−f_g)`, negligible
  at small `f_g`. The zero controls should barely move; `ISSUES: capture-on-strand-pure-ambig-undercall` is
  not this thread's.
* No strand information (`κ = ½`, or `σ_τ ≥ 1`): the strand term is flat in τ, the untruncated volume
  `σ_τ√(2π)` exceeds the domain and the marginal saturates at a constant, so `+ log(1 − f_g)` is a pure
  prior shift on every AMBIG slot — Beta(½,3/2) against Beta(½,½), median 0.17 against 0.50 at zero
  evidence. On unstranded data this is the whole effect, and the panel must price it (the unstranded ×
  capture-OFF stratum is in scope; `g00 ss.50` is a zero control). Normalising by the TRUNCATED volume
  instead (`M/V`, `V = ∫_{−1}^{1} e^{L} dτ`) would remove the saturation problem but inverts the boundary
  tail (`M/V ∝ √f_g/σ_τ`, because the arcsine weight makes `M`'s tail heavier than `V`'s) and is refused on
  paper: it would break the zero controls.

**The prototype** (`tilt_measure.py`, the session scratchpad): `_psi` patched to add `log(1 − f_g)` on AMBIG
cubes; arms `reference` / `jeffreys`; judged on the metric per stratum with both zero controls on both panels,
the census bands, and the shared-exon stress in its mono (the bare factor) and spliced (with the guards) forms.

**The Jeffreys volume is REFUTED on the ladder (2026-09-13).** `+ log(1 − f_g)` on the AMBIG cube (the arm
`jeffreys`): Σ|Δ gDNA| both axes, g00 excluded — stranded OFF 249,703 → 267,373 (+7 %), stranded ON 471,202 →
761,711 (+62 %), unstranded OFF 308,529 → 331,064 (+7 %), deferred 3,622,383 → 5,263,361 (+45 %); the four g00
rows 2–4 % better (499 → 478, 406 → 397). The test chromosome +1–3 % on every stratum. The census says where:
the gDNA-rich AMBIG slots under capture — g98 ss.99 ON strand-pure band 34,056 → 120,845, g50 ss.99 ON
38,823 → 81,151, and the both-strand bands 2–4× — exactly the regime the derivation flagged: where `a(λ)` is
small the strand term does not constrain the tilt, the untruncated volume exceeds the domain, and the term is
a prior shift toward RNA on slots that are mostly gDNA. The stress confirms the term's own claim — on the
mono toy the balanced exon's √n push falls 0.99 → 0.48 at 50k and 500k, on the spliced toy nothing moves —
so it removes the Occam factor only where nothing guards the slot, and costs a stratum where the guards
are the answer. No saturating variant survives on paper (§11: the truncated volume inverts the boundary
tail). Refused. The remaining candidate of the same derivation is the PROFILE form: the strand term enters
the λ posterior through `max_τ L(λ, τ)` (the cap, flat inside it), and the tilt is integrated only for the
delivered rows and the tilt read-out — as a cube, `ψ' = ψ + [L_max(λ) − lse_θ(L + log w)]`, a per-(slot, λ)
constant that leaves the θ-shape (hence `w_pos`) untouched. It has no volume factor, no saturation regime
and no measure; at a strand-pure slot it drops the tail's `1/√f_g` width factor, which is the family of
`ISSUES: capture-on-strand-pure-ambig-undercall`. Arm `profile`, judged the same way.

**The profile form is REFUTED too (2026-09-13)**: ladder 253,117 / 687,783 / 308,557 / 3,653,341 (+1.4 / +46 / 0 /
+0.9 %), g00 marginally better, the test chromosome +0.2 / +0.8 / 0 / 0 %, the same g98 capture-ON AMBIG bands
(19 → 131, 10 → 88); mono toy 0.997 → 0.49 at 500k balanced, spliced unchanged. The thread is closed: the width
factor is evidence the panels reward; `ISSUES: strand-marginal-volume-factor` carries the record.

## 12. The strand-pure AMBIG under-call — the dissection (`ISSUES: capture-on-strand-pure-ambig-undercall`, 2026-09-13)

**Where it is.** The per-slot census on the ladder's stranded rows (`dissect/census_slots.json`): the AMBIG slots
whose RNA is on one strand are MOSTLY gDNA — mean true `f_g` 0.53 / 0.85 / 0.89 / 0.85 / 0.36 by depth bin on
`g50 ss.99 ON`, 0.97 on `g98` — and the under-call sits on the deep ones: slots with ≥ 300 fragments carry
−23k of the −28k net on `g50 ss.99 ON` and −27k of −30k on `g98`. Boundaries carry three quarters of it.
By stage: the strand-only solve reads −97k on `g50 ss.99 ON`, the local solve (factory rows) −44k, the shipped
answer (messages, refits) −28k; off capture (`g50 ss.99 OFF`) the same structure starts at −23k strand-only
and the prior repairs it to −1.2k. So the bias originates in the own solve and the prior repairs it where
the landscape is sharp (no capture) and cannot where it is broad (capture).

**The mechanism.** At a slot whose RNA is strand-pure the truth sits AT the strand cap: with `τ = 1` the
split `p̂ = ½ + (1 − f_g)(κ − ½)` identifies `f_g`, exactly as at a single-strand slot. The AMBIG solve does
not know the tilt is pure: every `f_g` below the cap fits `p̂` with a slightly impure tilt, the marginal is
spread over `[0, cap]` (∝ the width factor and the arcsine at `τ̂(λ)`), and its median lands below the cap.
Prior-free, one slot: `f_true` 0.50 reads 0.31–0.37, 0.85 reads 0.69–0.80, 0.97 reads 0.85–0.95 (n = 30 …
30k). This is a structural bias of marginalising a tilt whose truth is at its vertex.

**The candidate the dissection names: a tilt atom.** The AMBIG tilt's hypothesis space is {pure +, pure −,
mixed} at equal reference weight; the pure hypotheses are single-strand solves inside the cube (τ = ±1, no
tilt parameter), the mixed one today's continuous θ integral with its measure normalised to the domain
(`dθ/π`). Where the data are pure the pure hypothesis explains them with no parameter and wins the Occam
contest; the smoke test reads 0.46–0.50 / 0.84–0.88 / 0.91–0.97 at the three truths above. Its cost: at an
INTERIOR tilt the pure hypothesis at `f_g = cap` ALSO explains `p̂` with no parameter, so the atom pulls a
both-strand slot toward the cap (`f_true` 0.30 at τ 0.5: 0.46 → 0.64, `w₊` 0.79 → 0.89). The ladder prices
the net (2,347 strand-pure against ~1,500 both-strand AMBIG slots on `g50 ss.99 ON`); the arm is
`tilt_atom.py`, judged on the metric, both zero controls, the census bands and the stress.

**The atom's A/B (2026-09-13).** Ladder, Σ|Δ gDNA| both axes, g00 excluded, reference → atom (equal weights) →
witnessed atom (a pure hypothesis admissible only without a delivered RNA level on the other strand):
stranded OFF 249,703 → 249,407 → 248,827; stranded ON 471,202 → 424,712 → 427,069 (−9.9 % / −9.4 %);
unstranded OFF 308,529 → 308,763 → 307,992; deferred 3,622,383 → 3,661,197 → 3,610,370. The zero controls:
`g00 ss.99 OFF` 405.5 → 445.1 → 433.3 (+28 fragments under the witnessed form, at both-strand slots pulled
toward the cap: the (0.2, 0.5] band 24 → 55, the (0.05, 0.2] band 15 → 22), `g00 ss.99 ON` 194.0 → 198.3 →
197.4, both unstranded g00 rows within a fragment. The test chromosome: −0.1 / +0.2 / −0.06 / +0.09 % under
the witnessed form. The census on `g50 ss.99 ON`: the strand-pure band 38,823 → 23,824 → 27,683 and its tilt
error 17,123 → 7,817 → 10,582; the near-pure band 14,526 → 7,884 → 9,259; the both-strand bands (0.05, 0.2]
6,005 → 5,002 → 4,835 and (0.2, 0.5] 11,301 → 11,964 → 12,246 with the tilt error 2,500 → 3,663 → 3,020 —
the witness removes about a third of the atom's cost at both-strand slots and keeps most of its win. The
spliced stress is unchanged under the witnessed form on every row (a level on each strand rules both pure
hypotheses out); the mono stress shows the plain atom's failure bare (no witness: a 20 %-minor exon 0.28 →
0.40, its tilt read as pure). The residual cost lives at both-strand slots that hold a delivered level on
only one strand.

## 13. What the atom points to — the study (2026-09-13)

**Was the tilt ever a message?** No. The ruling in `messages/__init__.py` (and `DESIGN.md` §6b.13) says the tilt
has NO lane — the two RNA levels constrain it through the shares — and the only strand-shaped row in the layer
(`strand_row_logodds`) is a composition row over λ from a single-strand exon's own strand term. The deleted
relay carried composition and flux rows. Nothing about the tilt has ever travelled, and the study says it
should not: a tilt is a RATIO of two levels at one slot, slot-specific by structure; what is conserved along
a locus is each strand's LEVEL, which the RNA lanes already carry as a lower bound.

**What a "strand-pure AMBIG slot" is, structurally** (`study/structure.json`, `g50 ss.99 ON`, regions): the
RNA's strand is an exon and the other strand is that gene's INTRON at this slot — 465 slots, 174k RNA — or
both genes have an exon here and one is silent at the locus — 214 slots, 68k RNA; and 1,484 strand-pure
BOUNDARIES carry three quarters of the under-call (−21k of −28k). For 78–81 % of pure slots the absent
strand's gene has RNA elsewhere in the same locus: "pure" is not "the other gene is silent", it is "the
other strand's RNA here could only be NASCENT" — the nascent scope ruling's own case (`DESIGN.md` §0b:
absent unless abundant evidence), which the AMBIG cube does not know: it gives the intron strand's RNA the
same reference measure as the exon strand's, and the tilt marginal then pays the Occam price of a
continuum whose truth is at its boundary.

**Presence against witness** (`study/crosstab.py`): the lanes delivered a level on neither strand at ~40 % of
pure slots, on one at ~45 %, on both at ~15 % (a delivered level is a lower bound, sometimes near zero, so
"witnessed" is coarse). The atom's win sits on pure slots of every witness state (`g50 ss.99 ON`: none
−14.5k → −5.5k, one −12.9k → −4.2k); its cost on mixed-truth slots with one witness (+2.0k → +7.2k plain,
+4.4k witnessed) — and for 98 % of those (546/558) the unwitnessed strand's gene HAS RNA elsewhere in the
locus, i.e. the presence witness exists at the locus and the lane did not deliver it to this slot. The
delivered-level witness recovers a third of the cost; the remainder is a lane reach question, not a tilt
question.

**So the atom is not a special case.** It is Axiom 0's opportunity geometry applied to the tilt: the
population set at a slot is `{gDNA} ∪ {RNA+ if free_pos} ∪ {RNA− if free_neg}`, and the AMBIG cube treats both
RNA members as equally weighted continua; but PRESENCE per strand is structural — exon of s (the gene's
RNA, if expressed) against intron of s (nascent only, sparse by ruling) — and locus-level (the gene is
expressed), which the RNA lanes carry as a lower bound. The tilt's hypothesis space {pure +, pure −, mixed}
makes presence discrete, and the witness that selects among them should be structural first (the per-strand
exon bits: a strand whose RNA here could only be nascent cannot be the pure carrier) and delivered second
(a level on a strand rules the opposite pure hypothesis out). The arm `atom_s` prices that. What it points to
beyond itself: an AMBIG slot that is one strand's exon and the other's intron is a single-strand solve with a
nascent allowance, and most AMBIG slots are that; the true both-strand solve is the exon∩exon overlap of two
expressed genes (58–108 slots per ladder row, 2–3 % of AMBIG slots). The cube is the right model for those
and an over-general one for the rest.

**The structural witness adds nothing measurable (2026-09-13).** Arm `atom_s` (the per-strand exon bits on top
of the delivered-level witness): ladder 248,789 / 426,720 / 308,050 / 3,615,015 against the witnessed atom's
248,827 / 427,069 / 307,992 / 3,610,370; `g00 ss.99 OFF` 424.5 against 433.3. By presence truth on
`g50 ss.99 ON` the two are within 2 %: pure −12,656 / −13,005 net, mixed +7,923 / +7,489. The slots that pay
are exon∩exon overlaps and boundaries, where the exon bits exclude nothing; the residual cost is the mixed-
truth slots whose second strand's gene is expressed in the locus (98 %) but reached this slot with no level —
a lane-reach question, to be filed as its own entry. The delivered-level witness alone is the landing form.

## 14. The encompassing-transcript audit (owner's case, 2026-09-13)

**The case.** A single-exon TB− (10,000–30,000) encompassing a two-exon TA+ (11,000–12,000, 19,000–20,000);
`audit_encompass.py` in the session scratchpad. The lanes' face rule says TB−'s level crosses every boundary
here (each carries only TA+'s bits, every region admits −) from its two single-strand exons, and the
patched tree confirms it on the g50 donor: the − level is held on both sides of every slot from 11,000 to
20,000 and delivered into every AMBIG cube; its wall reads 1.957/bp against true − densities of
1.92–2.22/bp; the exon∩exon slots read `f_g` 0.006 / 0.011 against 0.005 / 0.006 (0.29 / 0.44 with no level)
and the tilt 0.51 / 0.52 against 0.53 / 0.52. The witnessed atom equals the reference there.

**Four couplings that silence the RNA lanes, none of which a whole chromosome shows.** On the shipped tree
the case delivered NO level anywhere:
* A — `TransferPolicy.prepare` returned an empty layer (no faces, no lanes) when the intron factory had no
  rows, i.e. no coarse intron anywhere in the chain. Fixed: a factory with nothing to say is all-zero rows.
  Gate `test_the_lanes_are_built_when_the_intron_factory_has_no_rows`. This is also why every intron-less
  TOY (`TA_single_exon`, `one_exon`, the mono shared-exon stress) ran with no message layer at all — the mono
  toy's "unguarded" reading was the layer being off, not the lanes not reaching.
* B — the RNA lanes were built only if the gDNA lane existed, so a library with zero gDNA density had none.
  Fixed: each lane exists iff its own coordinate does (`rna_lanes` takes the grid, not the gDNA lane). Gate
  `test_the_rna_lanes_are_built_without_a_gdna_lane`.
* C — (not real) the + lane exists with a zero coordinate rather than being absent; the cube delivery was
  never blocked by it. The gate written for it (`…when_the_other_strand_has_no_coordinate`) passes and stays.
* D — OPEN: a strand's flux level — the junction's certified estimate of its exon's RNA — is built only if
  the strand's library coordinate is positive, and that coordinate is the density over the strand's
  SINGLE-STRAND exons. TA+ has none, so its certified flux is silently unused and no + level reaches its
  exons (the audit's `wit` column reads `·−`). The level is absolute, so any positive reference density is a
  coordinate; the candidate is one RNA coordinate for both strands (the pooled single-strand exon density)
  or the flux's own. The owner's call: it changes which sources exist on a real chromosome only where a
  strand has no single-strand exon, i.e. never on the ladder.
* The g00 donor delivers nothing for a fifth, ruled reason: the derived deadband (`strand_discriminability`,
  `1/N_gdna`) declares the strand channel uninformative on a library with exactly zero gDNA, so no
  single-strand exon emits a strand-derived level. The ladder's g00 rows carry a small fitted gDNA count and
  the channel is live there; a truly gDNA-free real library would not be. A ruling, not this thread's.

**The rung-0 identity gate** (`test_transfer_policy.test_an_evidence_free_transfer_is_byte_identical_to_silence`)
asserted coupling A as a floor — no factory rows ⇒ the transfer policy IS silence. With A fixed it fails by
design. The floor the tool has is `SilentPolicy`; retiring the gate is the owner's decision and the fix is
left uncommitted with it failing until ruled.
