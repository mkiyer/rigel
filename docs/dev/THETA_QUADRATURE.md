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

## 5. The plan

1. DERIVE — this note. The two recorded numbers (0.03 rad at 500k; ≤ 4 fragments at 50k) are reproduced by §3.
2. PROTOTYPE outside `src/`: patch `_solve_logodds`' θ-sum (A: the endpoint weights; then B) in both bindings,
   as `w5/theta_mesh.py` did for `_tilt_grid`.
3. A/B on both panels with both zero controls: `calibration_vs_oracle.py`, the `g00` rows under
   `--set calibration.sweep_n_tilt=30` and at the default 60, patched against unpatched, the metric per stratum;
   the panel; the profiler (the cube is the cost).
4. `src/` only after the A/B, one mechanism at a time; the ruling to `DESIGN.md` §6b.15, the derivation to
   `EQUATIONS.md`, the issue closed with its numbers.
