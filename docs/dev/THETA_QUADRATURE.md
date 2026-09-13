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
