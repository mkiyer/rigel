# DERIVATION — the arcsine magnitude coordinate

    ⚠ **A DEV DOC**: the working derivation for `ISSUES: arcsine-magnitude-coordinate`, step 2 of the
    thread (step 1's dissection record is §0). Nothing here is authoritative; when it settles, the
    derivation MOVES to `EQUATIONS.md`, the dissection verdict into the issue, and this file is deleted.
    Written 2026-09-01 on the tree at `bae76391`.

## §0 Step 1 — the dissection (2026-09-01), and the attribution VERDICT

Instruments: `worst_objects.py --arm final`, `calibration_walk.py` (certified drained-frame
`slot_truth`, COMPOSITION+FIELD), plus a scratchpad banding of the walk's own C/F beliefs by certified
`true_f_g`. Conditions: `g98 ss.99 ON` (the 0.8.0-metric worst), its OFF sibling, and
`g50 ss.50 OFF` (`solvability_audit`'s worst) as the contrast.

**The stage ladder** (Sum|f_g − true|·mass, fragments; C = message-free local solve, refit 0):

| condition | C local | D C+msgs | E C+refits | F shipped |
|---|---|---|---|---|
| g98 ss.99 ON | 299,319 | 512,131 (+212,812) | 298,597 (−722) | 456,838 (+158,240 over E) |
| g98 ss.99 OFF | 169,075 | 221,501 (+52,425) | 126,467 (−42,608) | 208,306 (+81,839 over E) |
| g50 ss.50 OFF | 1,682,749 | 221,501-scale drop | — | 189,133 |

**The banding of C's error by certified truth** (signed ≈ |err| everywhere below — one-sided):

* `g98 ss.99 ON`: slots with truth EXACTLY `f_g = 1` (`n_rna = 0`) carry **56 %** of C's error
  (−167,483 of 299,319), truth ≥ 0.99 carries **79 %** — all UNDER-call. F deepens the same band to
  −229,785 (messages add −62 k there).
* `g98 ss.99 OFF`: exact-vertex slots carry **62 %**, truth ≥ 0.99 ~69 % — same shape without capture.
* Mirror: truth-EXACTLY-0 slots show pure OVER-call (+54 / +324), matching the g00 baseline rows
  (0 under / 67,571 objects over).
* `g50 ss.50 OFF` is a DIFFERENT defect: 83 % of its C error is OVER-call at truth ∈ (0, 0.2) — the
  unstranded blind-slot default — and messages+refits then remove ~89 % of it (1.68 M → 189 k).
  That condition belongs to the message-value thread, not this one.

**The structural check**: of the 37,426 exact-vertex slots at `g98 ss.99 ON`, the only ones reading
`f_g = 1.0` are 3,219 slots with `solvable = False` — signature-locked init, not solver output.
**No solver-reached slot can read either vertex**, because no grid point represents one
(`f_g = sigma(lam)`, `lam ∈ [−10, 10]` ⇒ `f_g ∈ [4.54e-5, 1 − 4.54e-5]`).

**VERDICT — the attribution the handoff demanded is CONFIRMED, with one word changed.** The error
mass sits on NEAR-VERTEX objects, pre-message, on BOTH capture arms. The alternatives are refused by
the same tables: not the capture-ON divisor (OFF carries the same one-sided vertex band at 62 %), not
the EB anchor (same reason), not the refits (E − C = −722). Messages ADD harm on top (+158 k at ON).
⛔ **"Vertex-BLOCKED" is the word this section originally used and §3 refutes it**: the objects are
not blocked by the coordinate — no solver-reached slot reads a vertex, but the grid's representable
ceiling is `4.5e-5` from the vertex while the shortfall is `0.026`, so the block is not what pins
them. The attribution is to the near-vertex BAND; the mechanism is §7a④.

**The sharpening the dissection adds** (this constrains what the build may claim): the read-out is the
posterior MEDIAN interpolated on the grid coordinate, and a median is transform-EQUIVARIANT — so a
coordinate change alone moves nothing wherever the posterior density (against the same measure) is
unchanged. At a deep exact-vertex slot the strand channel releases `f_rna ∈ [0, w]`, `w ≈ 2/sqrt(n)`,
and the Beta(½,½)-referenced median lands at `f_rna ≈ w/4`. Measured against observed `fg_loc`:
n = 1,000 predicts 0.9842 (observed deep-slot mean 0.9845 ON); n = 55,362 predicts 0.9979 (observed
0.998). That interior shortfall is the honest posterior read and is NOT the coordinate's to fix.
What the coordinate structurally changes is exactly four things — ⛔ **(i) is REFUTED by §3 and (iv)
with it; (ii) and (iii) survive and are confirmed in §7a**: **(i)** the vertices become
representable states (today they are not, at any depth, for any evidence); **(ii)** the fixed-`L`
bracket and its measured wall (the module's stamped L=20 table: g00 error ×0.32 when the bracket
widens — the fitted landscape pushes against `sigma(−L)` as a wall) are retired by construction;
**(iii)** the grid measure becomes the reference exactly (§2), retiring the written reference terms
and the reference-location graveyard; **(iv)** claims and precisions get a coordinate in which every
share mode is in-domain (§5) — the `TRAPS: off-grid-message-mode` family's structural end. The A/B
must therefore expect its wins at the extreme bands and the zero controls, not at deep interior
shortfalls — "how much the panel moves is the A/B's to say" stands.

## §1 The measure bookkeeping (checklist point 1)

**What fact 2's "no Jacobian is written" actually is.** The model's prior content per component group
is a density in LOG-rate, `p_g(log rho_g)` and `p_r(log rho_r)`, with base measure
`d(log rho_g)·d(log rho_r)`. The solve is conditional on the slot's total density `rho`. Change
variables `(log rho_g, log rho_r) → (f, rho)` with `rho_g = f·rho`, `rho_r = (1−f)·rho`:

    d(log rho_g) d(log rho_r) = df d(rho) / (rho · f(1−f))       (the 2×2 determinant)

so at fixed `rho` the induced base measure on `f` is `df / (f(1−f)) = dlam` — **flat in λ**. That is
the whole cancellation: the λ grid's own spacing IS the conditioned base measure, so fitted arms are
pure point evaluations, and the two written reference terms `+½·log f + +½·log(1−f)` are the Jeffreys
log-rate reference densities `p_c(log rho_c) ∝ rho_c^{+½}` evaluated pointwise (the `rho^{½}` of the
total is an additive constant). Fact 2 re-derived, λ-specific as the handoff says.

**The φ form.** Put `f = sin²(phi)`, `phi ∈ [0, π/2]`. Then

    dlam/dphi = 4 / sin(2·phi)   ⇒   log(dlam/dphi) = −½·log f − ½·log(1−f) + log 4.

A grid flat in φ must carry the measure conversion `+log(dlam/dphi)` to represent the same model —
and that conversion is exactly the NEGATIVE of the two reference terms. So, writing each arm as
(reference + optional fitted logP) as today:

    psi_phi = strand + [−½·log f − ½·log(1−f)]                      (the measure conversion)
                     + [½·log f + logP_g]  + [½·log(1−f) + logP_r]   (the two arms)
            = strand + logP_g + logP_r + (likelihood rows) + (messages) + const.

**The φ-analogue of fact 2 is stronger than the λ one: neither a Jacobian NOR a reference is ever
written.** Each arm contributes its fitted `logP(log rho_c(phi))` pointwise when fitted, and NOTHING
when unfitted. The cancellation is per-component, so it keeps holding as each arm acquires a fitted
prior — same property the λ form had, minus two written terms. Every other ψ term (strand mixture,
the intron-factory NB rows, `lam_rows`, every message penalty) is a pointwise function of the
composition and never carried a Jacobian in either coordinate.

## §2 The reference in φ (checklist point 2)

The flat measure in φ pushes forward to `Beta(½,½)` in `f` (the arcsine law — `dphi ∝
df/sqrt(f(1−f))`), which is exactly the Berger–Bernardo reference's `f_g` marginal recorded in the
module (fact 3's derivation and `_JEFFREYS_REF`'s two agreeing routes). So `psi_ref` vanishes into
the grid: **the reference becomes literally the uniform measure**, as the issue conjectured — verified
against the BB structure rather than assumed:

* the BB reference factorizes as `Beta(½,½)(f_g) × arcsine(tau)`; the θ grid already supplies the
  tilt factor (fact 3: `|d tau/d theta| = cos(theta)` cancels `(1−tau²)^{−½}` identically), and the
  φ grid now supplies the magnitude factor the same way. The full 2-D reference is the uniform
  `(phi, theta)` measure; the two cancellations are independent because the λ↔φ Jacobian is
  θ-free and the τ↔θ Jacobian is λ-free — the measure factorizes, so fact 3 survives untouched.
* the tilt needs no per-strand log-rate bookkeeping in either coordinate: no fitted density is ever
  specified on a per-strand rate (the arms are the TWO-GROUP split by design), so the third dof takes
  the pure BB conditional, which the θ grid is. Unchanged.

Consequences carried by construction: the reference-location machinery stays deleted (there is no
written reference to relocate); `sigma(L)`-bracket state-space limits and the `L`-invariance
acceptance test are replaced by **no bracket at all** plus a `K`-invariance acceptance test (§4);
`ISSUES: psi-lambda-bracket-unshipped` dissolves (a bounded coordinate has no bracket).

## §3 Vertex endpoints (checklist point 3) — ⛔ THE PREMISE IS REFUTED, MEASURED

    ⛔⛔ **THIS SECTION'S FIRST DRAFT WAS WRONG AND THE PROTOTYPE'S GATES REFUTED IT.** It is kept
    with the refutation because the claim it made is the ISSUE's own headline claim ("belief cannot
    sit at a vertex … a vertex is a finite interior endpoint instead of ±∞"), and the next session
    must not re-derive it. The measurements are `arcsine_proto.py --self-test` gates ④.

**Fact (a): the posterior MEDIAN cannot reach a vertex in ANY coordinate — it is bounded away by
half of the outer bin.** For the last bin holding posterior mass `m`, the CDF crossing fraction is
`t = (0.5 − (1−m))/m = 1 − 0.5/m`, which is `≤ ½` for every `m ≤ 1`. So the read-out never passes
the outer bin's MIDPOINT however concentrated the posterior is. Under uniform-φ at `K = 60` the
ceiling is exactly `1 − sin²(dphi/2) = 0.999828662` (gated, attained). This is a property of the
median read-out, not of λ: the vertex is unreachable because the estimator is a quantile of a
continuous posterior with no atom there, and no reparameterisation changes that.

**Fact (b): in `f`-space the λ grid is FINER near the vertices than uniform-φ, not coarser.**
`sigma(±L)` already sits `4.54e-5` from the vertex. At `K = 60` the top `f`-gap is `1.83e-5` under λ
versus `1.37e-3` under φ (75× coarser), and production's one-hot read-out is `0.999954602` under λ
versus `0.999828662` under φ — **λ is CLOSER to the vertex**. The λ grid concentrates resolution at
the vertices (that is what `logit` does); uniform-φ spends it evenly.

**Fact (c): therefore the measured shortfall is not a grid artifact, by four orders of magnitude.**
The baseline's under-call is `+0.15…+0.29` at true `f_g ∈ [0.999, 1]`; both grids' vertex
resolution is `~1e-4` or finer. The C-stage belief at exact-vertex slots is `0.974` (§0), which is
`~570×` further from 1 than λ's representable ceiling. **The solver is not pinned by the coordinate.
It is answering honestly against the Beta(½,½) reference and the strand channel's release width**
(§0's `w/4` prediction, confirmed against observed deep-slot beliefs at two depths).

⭐ **What the coordinate DOES change is a resolution REALLOCATION, and that is the real trade**:
uniform-φ is ~3.1× FINER mid-simplex (`f`-gap at ½: 0.0261 vs 0.0822 at `K=60`) and ~75× coarser at
the extremes. Plus the two structural wins that survive: the measure becomes exact (§2, gated to
3.3e-15 against `Beta(½,½)` bin mass) and there is no bracket (`ISSUES: psi-lambda-bracket-unshipped`
dissolves).

⭐ **What would actually move a vertex belief** — recorded so the thread can be re-aimed rather than
repeated: an ATOM at the vertex, i.e. the zero-inflation extension of the fitted-`logP_r` socket
(`_rna_arm`'s recorded landing point: "nothing fits `logP_r` yet, and the cost of that is measured —
a FIXED repulsion of 3.107 nats at `f_g = 0.999`, a 22:1 handicap"). That repulsion is the
mechanism, it is the same size as the measured shortfall, and it is a MODEL term — not a coordinate.
A coordinate change cannot remove it because it is the reference's content, not its parametrisation.

### What the original draft argued (kept for the record — superseded by (a)–(c) above)

The endpoints `phi = 0, π/2` are ordinary grid nodes. What ψ evaluates to there decides everything,
and it comes out of the terms — derived, never clamped:

* **strand mixture**: finite at both endpoints (mean `n·p` with `p` linear in the composition;
  variance frozen at the reference composition). Pure-gDNA data (`kappa` at its null) leave the
  endpoint LIVE; strand data inconsistent with the endpoint penalize it finitely, in proportion to
  the evidence — never a veto, matching the channel's own information content (`I_f` is bounded at
  the vertices — `ISSUES: reference-prior-refuted-at-concept-level`'s recorded fact).
* **genuine likelihood rows** (the intron-factory NB factor; the certified-flux `lam_rows`): NB at
  component rate 0 is finite iff the observed count is 0 and `−inf` otherwise. **Data with counts
  forbid the vertex; zero counts reach it** — the checklist's desired behaviour, and it emerges from
  the likelihood with no new mechanism. (These rows arrive pre-evaluated on the grid: their builders
  must evaluate the φ grid's `f` values, endpoint included — re-evaluation, not interpolation, since
  the λ grid's f-range never covered the endpoints.)
* **fitted arms**: a proper fitted density on `log rho_c` tends to `−inf` as `rho_c → 0`, so a fitted
  arm's POINT value kills its component's vertex — correctly, in one precise sense: a log-rate density
  asserts the component exists at some level; it cannot say "absent". TRUE absence is a point mass —
  the zero-inflation extension of the fitted-`logP_r` socket (`_rna_arm`'s recorded repair landing
  point). That is a MODEL change, deliberately out of this build's scope; the coordinate's job is to
  make the state representable so that extension is writable at all. Meanwhile the fitted arm's mass
  near the vertex is handled by quadrature, next.
* **the endpoint-bin quadrature** — the one place point-evaluation × uniform weight is not adequate.
  Terms in `log rho_c` are unboundedly steep in φ at their component's endpoint
  (`d log rho_r/dphi = −2·tan(phi) → −∞`), so the endpoint BIN's mass must be the bin INTEGRAL, not
  the node value × Δφ (the node value can be `−inf` while the bin holds decisive mass — this exact
  case is the measured g00 landscape-vs-bracket wall, where the fitted gDNA arm's mode lies below
  the representable range and the L=20 arm recovered ×0.32). Substituting the steep term's own
  coordinate `t = log rho_c` inside the bin turns the integral into an ordinary 1-D integral of
  (fitted density) × `exp(t/2)` × (smooth terms ≈ const) — the `e^{t/2}` is the arcsine measure's
  half-power reappearing — computable by a few sub-nodes spaced by the fitted arm's own scale.
  Everything is derived from the terms; no floor, no epsilon. Interior bins keep plain node
  evaluation (every term is smooth there; midpoint quadrature to second order, same as today).
* **read-out at the endpoint**: `_posterior_median_fg`'s histogram-CDF interpolation carries over
  verbatim with `phi` as the interpolation axis (it is the same construction that made λ
  interpolation transform-invariant; the outer half-bins already mirror). A posterior concentrated
  in an endpoint bin reads out through that bin's own within-bin CDF — so `f_g` can now BE 0 or 1
  to within the endpoint bin's resolution, and exactly 0/1 when the mass sits at the node.
* **the emitted moments**: `Var(log f_c)` grid-moments hit `log 0 = −inf` at an endpoint node with
  mass. The same quadrature answers it: every emitted log-moment is the BIN-integrated expectation
  (finite — `∫ log f · f^{−½} df` converges), exact to quadrature order on interior bins and
  computed by the endpoint sub-nodes at the ends. One convention, applied uniformly, replacing
  nothing on interior-supported posteriors.

Falsifications this section fixes in advance: the `g00` panels (truth at the `f_g = 0` vertex
everywhere — today pure over-call; belief must now be able to resolve to 0, with the fitted-arm
endpoint-bin mass carrying it) and the near-pure `g98` bands (the C-stage under-call must not get
WORSE; the deep-slot interior shortfall is expected to stand until the `logP_r` socket is filled).

## §4 Grid spacing (checklist point 4)

The `dlam = 0.3390` note is measured-in-λ and does not transfer; derive the φ spacing from what the
grid must resolve:

* **the prior side**: flat-in-φ IS the reference, so uniform φ nodes give exact-to-second-order
  quadrature of the reference everywhere INCLUDING the vertices — the Beta(½,½) "divergence" in `f`
  is a coordinate artifact; in φ the density is constant. (In λ the reference decays like
  `e^{−|lam|/2}` and 0.86 % of its mass lay outside L=10 entirely.)
* **the data side**: the strand channel's information in φ is bounded by O(n) uniformly in position
  (`I_phi = 4·f(1−f)·I_f`, and `I_f` is bounded at the vertices by κ-compression), so posterior
  φ-width is ≳ `1/(2·sqrt(n))` everywhere — position-uniform, which λ does not have (λ-width diverges
  at the vertices and shrinks mid-simplex). A UNIFORM φ grid therefore has a constant
  nodes-per-posterior-width budget across slots and depths: uniform-φ is the derived choice, not an
  asserted one.
* **K**: keep `K = 60` (AMBIG axis) and `K_ss = 256` initially — one variable at a time. For
  calibration: uniform-φ at K=256 has smallest interior `f = sin²(dphi) = 3.79e-5` (vs the λ grid's
  `4.54e-5` floor) and its mid-simplex `f`-resolution is ~3× finer than λ at equal K; at K=60 the
  extreme-tail interior resolution is coarser (7.1e-4) and the endpoint bins carry the difference —
  which is what §3's bin quadrature is for. **The acceptance test is K-invariance** (double K, the
  answer must not move), replacing L-invariance; unlike L there is no second, hidden role for K to
  leak into (L silently set the prior strength when ψ was improper — a bounded coordinate with a
  proper measure has no such channel).
* **endpoint sub-nodes**: count and spacing derived from the steep terms' own scale (the fitted
  arm's sd in `log rho`; the NB rows' curvature), with their own convergence check. No constant.

⭐⭐ **WHAT THE PROTOTYPE ACTUALLY BUILT, AND WHY IT IS SIMPLER THAN THE ABOVE.** Nodes are the
MIDPOINT rule — `phi_k = (k + ½)·(π/2)/K` — not `linspace` including the endpoints. Consequences,
all derived rather than chosen: every node is strictly INTERIOR, so `f ∈ (0,1)` and **every log-rate
term is finite at every node** — the endpoint singularity §3 wanted sub-node quadrature for does not
arise, so there is no quadrature scheme, no floor and no epsilon anywhere in the change; uniform
weights are the reference EXACTLY (gated to `3.3e-15` against `Beta(½,½)` bin mass); the outer bins
are `[0, dphi]` and `[π/2 − dphi, π/2]`, wholly inside the domain, so the read-out needs no mirrored
half-bins outside it (production's λ read-out mirrors, because λ's domain is unbounded). ⚠ The
residual the midpoint grid keeps is the outermost BIN under-representing a fitted arm whose mode lies
beyond it — the same failure mode as the λ bracket, but two decades further out. ⛔ And per §3(a) it
does not buy vertex reachability either: the median's ceiling is `1 − sin²(dphi/2)`.

## §5 Messages and precisions (checklist point 5)

Inventory of every ψ input channel and what the map does to it:

| channel | today | under φ |
|---|---|---|
| strand mixture | pointwise in composition | unchanged |
| `lam_logprior` (intron-factory NB rows) | `(m, K)` rows pre-evaluated on the σ(λ) grid | re-evaluated by its builder at the φ grid's `f` (endpoint per §3) |
| `lam_rows` (certified-flux stream) | `(n_slots, K)` NB rows on the λ grid | same — re-evaluated; endpoint behaviour is the correct counts-forbid rule for free |
| fitted arms (`CompositionPriors`) | evaluated at `log rho_c` on the σ(λ) grid + reference | evaluated at `log rho_c(phi)`, NO reference term (§1); endpoint bins per §3; `_regrid_global` interpolates in φ between φ grids |
| `gdna_imp` / `rna_imp` (Gaussians in `log f_c`) | pointwise penalties | rung A: identical pointwise (endpoint value `−inf` where precision > 0 — a positive-precision log-level claim asserts presence; recorded, defensible, and superseded by rung B) |
| `lam_imp` (the single-λ composition Gaussian) | Gaussian on λ itself | rung A: pointwise `−½·p·(lam(phi) − m)²` — same values at the same `f`, `−inf` at endpoints |
| `theta_imp` (tilt) | Gaussian on θ | untouched |

**Rung B — φ-native delivery (derived, A/B'd separately from the grid swap).** A message is a Laplace
(Gaussian) summary of the source's claim, and a Laplace summary is taken IN a coordinate. The derived
coordinate is φ: it is the variance-stabilised one, so posteriors are closest to Gaussian there at
every location including near vertices — exactly where log-f and λ summaries are worst (the measured
`_KNOWN_VIOLATIONS`: ~30 % of gdna log-share modes sit 10.7 nats OFF the log-share grid on
`g50 ss.50 ON`; λ modes read ±20.7 against a ±10 grid). Mapping a claim `x ~ N(m_x, 1/tau_x)` stated
in coordinate `x ∈ {log f_c, lam}`:

    mode:      phi_m = phi(x = m_x)                      (monotone map, always lands in [0, π/2])
    precision: tau_phi = tau_x · (dx/dphi)²  at phi_m    (dlam/dphi = 4/sin(2·phi);
                                                          d log f/dphi = 2·cot(phi);
                                                          d log(1−f)/dphi = −2·tan(phi))

Every share claim in [0, 1] maps to an in-domain mode: **`TRAPS: off-grid-message-mode`'s whole
family becomes structurally unsatisfiable for level/composition claims**, and the backbone's
mode-in-grid assertion tightens from "inside an epsilon-floored open interval" to "inside a closed
interval that is the coordinate's entire range". The message laws survive the map unchanged in form:
precisions stay precisions in a declared coordinate (`TRAPS: zero-the-precision-with-the-value` is
about zeroing, not about units; single-source-may-only-reduce is precision algebra, coordinate-fixed);
`hop_logvar` is a variance in log-TOTAL-density (the reframe channel) and never touches the
composition coordinate — transport and damping stay in the claim's native coordinate, and the
re-expression happens once, at delivery into ψ.

**Emission stays contract-compatible**: the solver's outputs (`f_g` median; `Var(log f_c)` moments;
`tau_lam` consumers) are posterior FUNCTIONALS, computable on either grid (§3's bin-integrated
moments). The prototype emits byte-compatible arrays; re-basing the belief state itself on
`(phi, tau_phi)` — which de-conflates "belief at the vertex" from "infinite certainty"
(`TRAPS: a-variance-cap-asserts-certainty`'s xfail family: in λ-state those are the same state) —
is the message-policy build's business, on the settled coordinate, not this build's.

## §6 The 2-D solve (checklist point 6)

The AMBIG cube becomes `(m, K_phi, K_t)` with axis 1 gridded in φ — same shape, same f32/f64 split,
same tiling (`_block_rows` unchanged; `sweep_n_grid` semantics unchanged). Nothing on the θ axis
moves (§2). The arms, NB rows, `lam_imp`, and `gdna_imp` are θ-independent in either coordinate, so
θ still integrates out cleanly; the λ↔φ measure conversion is θ-free, so the marginalisation and the
`f_g` median on the θ-marginal carry over verbatim with φ as the axis. The cube's existing
`1/(n+1)` count floor on per-strand log-shares (the τ = ±1 edges) already meets the φ endpoint
column (`f_act = 0` there): same floor, same recorded justification, superseded at rung B the same
way as the 1-D path. The `w_pos` tilt-share read-out (a ratio of 2-D posterior means) is a
functional — unchanged.

## §7a ⭐⭐⭐ THE MEASURED VERDICT (2026-09-01) — the coordinate is NOT the lever

Prototype: `arcsine_proto.py` (scratchpad, outside `src/`) — three module patches, 23 self-test
gates all passing, each perturbed. A/B harness: `arcsine_ab.py` (walk metric, certified drained-frame
`slot_truth`), `arcsine_release_ab.py` (0.8.0's release metric,
`calibration_vs_oracle.measure_condition`), `arcsine_kinvariance.py` (the acceptance test).
⛔ **The NOOP arm is byte-identical to unpatched production on every condition and both metrics** —
the harness prices something rather than changing answers by existing.

**① The derivation's core claim is CONFIRMED** (`self_test` gate ③, three likelihoods): λ-with-a-
written-Beta(½,½)-reference and φ-with-a-flat-measure are the same model, agreeing to `2e-8` on the
posterior median. Uniform-φ weights reproduce the `Beta(½,½)` bin mass to `3.3e-15`; keeping the
reference under φ shifts an asymmetric case by `0.017` (the perturbation). §1 and §2 stand.

**② The issue's headline premise is REFUTED** — §3 above, gated. The median cannot reach a vertex in
ANY coordinate; in `f`-space λ is FINER at the vertices than uniform-φ and its read-out is CLOSER to
1 (`0.999954602` vs `0.999828662` at K=60).

**③ K-INVARIANCE — the acceptance test — separates resolution from domain** (`g98 ss.99 ON`, rung C):

| K (ambig/ss) | base err | arcsine err | arc/base | mass-wtd \|base−arcsine\| | belief@truth==1 |
|---|---|---|---|---|---|
| 60/256 | 299,319 | 292,380 | 0.9768 | 6.484e-04 | 0.974050 / 0.974729 |
| 120/512 | 300,375 | 293,294 | 0.9764 | 6.577e-04 | 0.973975 / 0.974677 |
| 240/1024 | 300,929 | 293,769 | 0.9762 | 6.686e-04 | 0.973940 / 0.974651 |

Each coordinate converges to its OWN limit (mass-wtd `|K60 − K240|`: base 1.47e-04, arcsine
1.31e-04) while the gap BETWEEN them holds at ~6.5e-04 and does not shrink with K. **That residual
is the λ BRACKET's truncated domain, not quadrature** — refining a grid cannot recover a domain it
does not cover. So the coordinate's entire leverage on the vertex band is `6.5e-04` of belief,
against a measured shortfall of `0.026` — **2.5 % of the defect.**

**④ THE DEFECT'S ACTUAL MECHANISM, MEASURED** (`g98 ss.99 ON`, 34,207 SOLVED exact-vertex slots,
5.89 M fragments): the under-call tracks the strand channel's own release width across four decades
of depth.

| depth n | slots | median n | observed 1−f_g | predicted w/4 = 1/(2·√n) | obs/pred |
|---|---|---|---|---|---|
| [10, 32) | 3,377 | 16 | 0.07932 | 0.12500 | 0.63 |
| [32, 100) | 3,898 | 63 | 0.04549 | 0.06299 | 0.72 |
| [100, 320) | 5,055 | 146 | 0.02903 | 0.04138 | 0.70 |
| [320, 1000) | 4,335 | 650 | 0.01606 | 0.01961 | 0.82 |
| [1000, 3200) | 420 | 1,179 | 0.01192 | 0.01456 | 0.82 |
| [3200, 10000) | 114 | 4,966 | 0.00858 | 0.00709 | 1.21 |
| [10000, ∞) | 35 | 14,916 | 0.00392 | 0.00409 | 0.96 |

⭐ **The shortfall scales as `1/√n`** — the signature of a POSTERIOR WIDTH, not of a grid (a grid
artifact is constant in depth, or saturates at the grid's resolution). Derivation: the strand channel
releases `f_rna ∈ [0, w]`, `w ≈ 2/√n`; under the `Beta(½,½)` reference the posterior density there is
`∝ f_rna^{−½}`, whose MEDIAN sits at `f_rna = w/4`. **The solver is answering honestly.** The error
at an exact-vertex slot is the irreducible price of a CONTINUOUS prior with no atom at the vertex,
scored against a truth that sits exactly on it.

**⑤ The panel A/B, both metrics** (`arcsine/base`; <1 is better):

| condition | walk C | walk F | release Σ\|Δ\| |
|---|---|---|---|
| `g98 ss.99 OFF` (stranded × OFF) | 0.9881 | **0.8984** | **0.8984** |
| `g98 ss.99 ON` (worst in-scope) | 0.9768 | 1.0371 | **1.0371** |
| `g50 ss.99 ON` | 1.0040 | 1.0096 | — |
| `g50 ss.50 OFF` (worst by solvability) | 0.9971 | 1.0712 | **1.0712** |
| `g05 ss.99 OFF` | 0.9878 | 1.0081 | — |
| `g05 ss.50 OFF` | 0.9999 | 1.1074 | — |
| `g00 ss.99 OFF` (zero control) | 0.9701 | 0.9781 | **0.9781** |
| `g00 ss.50 OFF` (zero control) | 1.0000 | 0.8239 | 0.8239 |
| `g00 ss.50 ON` (zero control) | — | — | 0.9247 |
| `g00 ss.99 ON` (zero control) | 1.0073 | 1.0369 | — |

⭐⭐ **THE WHOLE LADDER (all 16, release metric, noop byte-identical on every one) SPLITS CLEANLY BY
STRATUM** — this supersedes the partial reading above:

| stratum | rows (low→high gDNA) | reading |
|---|---|---|
| stranded × OFF (IN SCOPE) | 1.008 · 0.939 · 0.898 | ⭐ WINS, and more as gDNA rises |
| stranded × ON (IN SCOPE) | 0.998 · 1.010 · 1.037 | mild loss, worsening with gDNA |
| unstranded × OFF (IN SCOPE) | 1.107 · 1.071 · 1.014 | ⛔ loses EVERY row |
| unstranded × ON (DEFERRED) | 1.024 · 1.007 · 1.002 | mild loss |
| `g00` zero controls | 0.824 · 0.925 · 0.978 · 1.037 | 3 of 4 win |

⭐ **The mechanism the split suggests**: the coordinate helps where the strand LIKELIHOOD is strong
(stranded × capture-OFF — the reference matters least there and φ's ~3× finer mid-simplex resolution
is a real gain) and hurts where the likelihood is flat (every unstranded row — there the posterior
essentially IS the reference, so the grid's tail resolution moves the answer directly). A real trade,
not noise.

⭐ **The shape**: the message-free LOCAL solve is neutral-to-slightly-better nearly everywhere
(C median ≈ 0.997, range 0.970–1.007 — consistent with ③'s small, real, bracket-sourced gain), while
the SHIPPED stage carries the spread (F range 0.824–1.107). ⛔ **Losing all three IN-SCOPE unstranded
rows is what disqualifies it**, not any single row (`TRAPS: never-pool-the-strata` — the ladder total
would have hidden it).
⚠ **Rung B (φ-native message delivery, §5) is NOT measured** — rung A delivers claims pointwise in
λ's coordinate, so the F-stage spread is the expected cost of messages built for λ, and rung B is the
honest candidate to recover it. But rung B cannot exceed the coordinate's total leverage, which ③
bounds at 2.5 % of the defect.

**⑤b THE RE-AIM'S TARGET POPULATION IS HALF THE DATA** (no solver; certified `slot_truth`, whole
ladder): **50.0 % of live mass — 72.9 M of 145.8 M fragments — sits EXACTLY on a simplex vertex.**
`g98 ss.99 OFF` 73.0 %, `g98 ss.99 ON` 50.0 %, `g50` rows 19.3–38.7 %, `g05` rows 3.7–21.6 %, all
four `g00` controls 100 % (pure RNA by construction). That is the mass a continuous prior cannot
reach in any coordinate — the atom is not a corner case.

**⑥ VERDICT.** The coordinate is a sound reparameterisation — the derivation holds, the measure
becomes exact, the bracket dissolves — and it is NOT the fix for the vertex under-call, which is a
MODEL defect the coordinate cannot reach. On 0.8.0's own metric it is MIXED-SIGN across in-scope
strata (a 10 % win at `g98 ss.99 OFF`, +3.7 % and +7.1 % losses elsewhere) with the zero controls
better. ⛔ It does not clear `TRAPS: panel-before-src` and must not land as-is.

⭐ **THE RE-AIM.** The mechanism ④ names is already named in the source, with its price:
`_rna_arm`'s docstring — "**NOTHING FITS `logP_r` YET, AND THE COST OF THAT IS MEASURED** … a FIXED
repulsion of **3.107 nats** at `f_g = 0.999` relative to `f_g = ½` (a 22:1 handicap) … objects whose
TRUE `f_g ≥ 0.999` carry **49–83 %** of all calibration error … read **0.13–0.23 below** the vertex.
⭐ The parameter exists now so that an estimator can close that asymmetry; the socket is not
speculative surface, it is the repair's landing point." That is the same band, the same direction and
the same magnitude as the 2026-08-31 baseline's finding, and it is a MODEL term. The candidate is an
ATOM at the vertex — "this slot holds NO RNA" as a point mass rather than a density limit — which is
a well-posed spike-and-slab / point-null and is exactly what a continuous reference cannot express in
any coordinate. ⚠ It needs its own DERIVE → DESIGN → PLAN → PROTOTYPE → A/B, and the arcsine
coordinate is a genuine ENABLER for it (an atom needs a representable endpoint), which is the honest
reason to keep this derivation rather than delete it.

**⑦ ⛔⛔ THE RE-AIM WAS PRICED IN THE SAME SESSION AND IS ALSO NOT A WIN.** `vertex_ceiling.py`,
three conditions (`g98 ss.99 ON/OFF`, `g50 ss.50 OFF`), oracle truth pinned and the whole chain
RE-SOLVED (`noop` pins 0 objects and is byte-identical — the instrument's own falsification passing;
`vertex_free` pins 69,850–78,820 objects, `vertex_all` 365,490–391,910). Read `mwae ALL` /
`Σ|err| ALL` only; `confidently wrong` and `solv%` are contaminated by the pin, which hands objects
certainty and so moves them into the confident population:

| arm | region Σ\|err\| | boundary Σ\|err\| | net |
|---|---|---|---|
| `vertex_free` (reachable population — THE CEILING) | +18,887 (2 better/1 worse) | −20,208 (3 better/0 worse) | **≈ −1,321, a WASH** |
| `vertex_all` (looser; includes own-evidence objects) | +94,490 | −18,237 | **+76,253, WORSE** |

⭐⭐ **A perfect vertex answer at every reachable object is worth approximately nothing net.** The
boundary axis improves on every row, the region axis degrades, and they cancel. Handing a vertex
object the exact answer changes what it BROADCASTS and the relay over-propagates it — which is why
this instrument RE-SOLVES rather than substitutes (`TRAPS: substitution-understates-a-source`).
⭐ **Taken with ④ this is the session's real conclusion**: the near-vertex error is an HONEST
posterior width AND is not removable by local certainty, so the question belongs to the MESSAGE LAYER
rather than to the prior or the coordinate. ⚠ Three conditions, on `vertex_ceiling`'s
pass-0-flavoured metric — the whole-ladder run and a `calibration_vs_oracle` equivalent are OWED
before any atom is built.

## §7 The prototype and its A/B (step 3 — as executed; §7a carries the results)

Parallel solve OUTSIDE `src/`, patched in-process over `_solve_regions_logodds_all` (the ceiling
instruments' pattern), one mechanism per arm:

1. **rung 0 — equivalence**: φ grid, SAME model content, endpoints carried but with the two written
   reference terms replaced by the measure per §1, K matched, messages evaluated pointwise (rung A
   delivery). Prediction: interior beliefs track the λ solver to quantisation (median equivariance);
   differences confined to extreme bands where the bracket used to truncate. This arm is the
   falsification harness, not a candidate.
2. **rung A — endpoints + bin quadrature ON** (§3): the candidate. Scored per stratum on
   `calibration_vs_oracle.py`, BOTH zero controls (`zero_controls.py` — g00 must move TOWARD zero;
   g98's exact-vertex band must not regress), reseed floor re-derived same-session, all 16 ladder
   conditions + the 30-condition test chromosome. The fl-gap arms are not implicated (nothing
   length-adjacent moves) — recorded, not skipped silently.
3. **rung B — φ-native delivery** (§5): separate arm, after A, since it changes message evaluation.
4. `rename_identity.py --freeze/--check` wraps any landing in `src/` (not the out-of-tree A/B);
   landing is the owner's call on the A/B table.

Honest expectations, from §0's sharpening: the deep-interior shortfall at exact-vertex slots
(`w/4`-scale) does NOT move under rung 0/A — it is the honest median of the reference posterior and
its repair is the fitted-`logP_r` socket (a follow-up model issue, enabled by this coordinate). The
measurable wins rung A is accountable for: the g00 zero controls (the measured bracket wall, ×0.32
headroom at that older metric), the exact-vertex reachability under zero opposing evidence, the
extreme-band truncation, and the retirement of `L`/bracket machinery. If rung A cannot demonstrate
those on the panel, the coordinate does not ship on outcome grounds — the structural argument alone
does not clear `TRAPS: panel-before-src`.

### ⭐ WHAT WAS ACTUALLY RUN, versus the plan above

* **The patch point moved and is simpler**: not a parallel solve over `_solve_regions_logodds_all`,
  but THREE module-level patches (`_logodds_grid`, the two arms, `_posterior_median_fg`) applied
  through an idempotent context manager that also rebinds each module which imported `_logodds_grid`
  BY NAME (`sweep`, `region_init`, `calibrate`; `rna_anchor` imports inside its function body and
  picks the module attribute up). That reaches the whole solve — both solvers, every row-builder,
  the AMBIG cube — without a second copy of the solver to drift.
* **Rungs 0 and A collapsed into ONE arm**, because the midpoint grid has no endpoint nodes and so
  no bin quadrature to switch on (§4's addendum). The `noop` arm is the falsification, and it held.
* **Rung B was NOT run** — it changes message evaluation and is a separate mechanism
  (`CLAUDE.md`: one mechanism at a time). It remains derived-but-unmeasured.
* **The prediction above was RIGHT on the part that mattered and WRONG on its premise**: the
  deep-interior shortfall did not move (correct — belief at truth==1 went 0.974050 → 0.974729), but
  "the exact-vertex reachability" was never available to win, because §3(a) shows no median read-out
  reaches a vertex in any coordinate. The zero controls DID improve, as predicted from the bracket.
* **`zero_controls.py`, the 30-condition test chromosome and the reseed floor were NOT run** on the
  arcsine arm — the verdict turned on ③'s leverage bound and ⑤'s mixed sign before they could add
  anything, and running them would have priced a mechanism already disqualified. ⛔ Recorded as a
  gap, not skipped silently: if the owner revives the coordinate (e.g. with rung B), they are owed.
