<!-- title: Count-space + relay — implementation plan -->
# Count-space state + forward-backward relay — implementation plan

**Status:** detailed plan for review. Supersedes the precision dimension of `count_space_solver_design.md`
(which scoped relay *out* in its §6 — that is now retracted: relay is folded **in**). Implements the relay
analysis in `ambig_intrinsic_solve_fix.md` §7. No magic numbers introduced without an explicit flag.

**Locked decisions (user):**
- **D1 — state = per-component pseudo-counts + per-component precisions; no fractions stored or propagated**
  (`{n_pos, n_neg, n_g}`, `{τ_pos, τ_neg, τ_g}`; density = `n_c / E_c`, `E_c` fixed). Fractions survive only
  as a *readout* projection for the EM prior.
- **D2 — forward-backward** message passing (see §3.1 for the tradeoff write-up that backs this).
- **D3 — mappability deferred** (folds into a future mappability-corrected effective length, lower priority).

---

## 1. The unifying model: a latent-density state-space along the genomic chain

The cleanest frame for *both* count-space and relay is a **state-space model** (a hidden-Markov / linear-
dynamical chain), which makes relay fall out as a standard operation rather than a bolted-on rule.

- **Latent state.** Three independent latent **density** processes evolve along the genomic chain: gDNA
  density `ρ_g(x)`, sense-RNA density `ρ_+(x)`, antisense-RNA density `ρ_-(x)`. "Along the chain" = along the
  ordered region↔boundary nodes of one reference.
- **Process (transition) noise.** Between adjacent nodes the density drifts by a biological step variance
  **`σ²_bio`** — exactly the `var~mean` quantity we already fit. gDNA drifts genomically (always); a strand's
  RNA density drifts only where that strand is *continuous* (`free_s` on both endpoints) and is otherwise
  **born** (TSS) or **killed** (TES) — strand birth/death in the state-space.
- **Observation.** At each node we observe the strand-resolved unspliced counts `(u_+, u_-)` (plus the
  boundary's one-sided spliced). The observation is the **Beta-Binomial strand mixture**: `p(+ read) =
  ½f_g + κf_+ + (1−κ)f_-` with `f_c = ρ_c / Σρ`. This couples the three latent processes *at the node* and is
  **non-Gaussian** (and degenerate in `f_g` at balanced counts — the count-zero-information principle).
- **Inference = forward-backward smoothing.** Each node's posterior density = (forward-predicted state from
  the left) ⊗ (its own observation) ⊗ (backward-predicted state from the right). Exact on a path graph for
  the Gaussian parts; the non-Gaussian observation is handled by a projected (assumed-density / EP) update —
  §3.4.

**Why this is the whole answer to relay.** A zero-count node contributes *no observation*. In a smoother that
is the **prediction-only** step: the latent density is carried through it with added process noise `σ²_bio`,
i.e. the upstream signal is **relayed**, attenuated by one step of biological drift. A tiny exon (small `E` →
its own observation is near-useless, §2.2) likewise contributes almost nothing and the state passes through.
Over a run of sparse/zero nodes the latent gDNA density is *smoothed* toward the neighbourhood mean — which is
exactly the true level the user described ("the average of which is the true global gDNA level"). The current
solver instead *kills* the signal at such nodes; that is the bug, and the smoother is its principled cure.

> **OPEN QUESTION 1 (message algebra).** `count_space_solver_design.md` proposed messages as Dirichlet
> pseudo-counts entering the simplex as `Σ_c α_c·log f_c`. That form **conflates density and precision**: a
> relayed-and-discounted message necessarily *moves the mean* (Dirichlet mean ∝ α_c), so it cannot relay "the
> same density at lower confidence" — fatal for relay. The state-space frame instead keeps a **per-component
> (density, precision)** message where relay preserves the density and only decays the precision. I therefore
> **recommend the (density, precision) form** and treat the pseudo-count `n_c = ρ_c·E_c` as the *signal* and
> `τ_c` as the *separate* reliability (D1 already separates them). The valued Dirichlet properties (a zero
> channel drops out; no fraction-Jacobian) are recovered another way in §4.2. **Decision needed: confirm we
> adopt (density, precision) and retire the `Σα log f` form.**

## 2. State, precision, and the three precision sources

### 2.1 Node state (D1)
Per node `i`, per component `c ∈ {+, −, g}`:
- `n_{i,c}` — pseudo-count (the *signal*): expected fragments of component `c` = `ρ_{i,c} · E_{i,c}`.
- `τ_{i,c}` — precision on the **density** `ρ_{i,c}` (the *reliability*; `0` = no information, `∞` = locked).
- Density readout: `ρ_{i,c} = n_{i,c} / E_{i,c}`, with `E_{i,g}` the gDNA-FL eff-len and `E_{i,±}` the RNA-FL
  eff-len (the existing `NodeGeometry.eff_gdna_*` / `eff_rna_*`; fixed across passes).

`n_c` and `τ_c` are **coupled at generation, independent thereafter** — own-data sets both from the same
counts; relay forwards `ρ_c=n_c/E` and *decays* `τ_c`. This separation is the whole point of D1.

### 2.2 The three precision sources (the user's insight, written as formulas)
A node's **own-observation** density precision combines three independent terms — length, count, structure:

```
τ_own,c  =  1 / ( σ²_bio(ρ_c)            # BIOLOGICAL: between-node spread (var~mean; also the process noise)
                + (n_c + a) / E_c² )      # SAMPLING:   Poisson density variance Var(ρ)=Var(n)/E² = (n+a)/E²
structure: free_s = False  ⇒  n_c ≡ 0, τ_c ≡ ∞ (locked; relays 0 on that channel)
```

- **Length** enters through `E_c²` in the denominator: a tiny exon has tiny `E` ⇒ `(n+a)/E²` is huge ⇒
  `τ_own` → 0 ⇒ **low precision purely from short length**, even at `n=0`. ✓ (the user's central point)
- **Count** enters through `n_c`: more fragments ⇒ smaller sampling variance ⇒ higher precision.
- **Structure** enters through `free_s`: a locked channel is `τ=∞`, `n=0` — it blocks (relays 0), the
  "infinite-precision ⇒ ignore+forward-zero" half that `free_s` already implements today.

> **OPEN QUESTION 2 (the sampling pseudocount `a`).** `Var(ρ) = (n+a)/E²` needs a prior pseudo-count `a` so
> that `n=0` is *low* precision (∝`E²/a`) rather than infinite. Candidates: `a=1` (one Poisson pseudo-obs,
> matches the current `_MSG_PSEUDOCOUNT`) or `a=½` (Jeffreys). This is the one genuinely new constant; flag for
> the no-magic-numbers review. (It is *not* free to drop: `a=0` reintroduces "zero count looks infinitely
> certain.")

### 2.3 What this retires
`f_pos/f_neg/f_g` and `var_*` on `NodeBelief` → `(n_c, τ_c)` per component. The fraction lattice survives only
inside the **local observation update** (§3.4) and the **readout** (§4.3); it is never stored or sent.

## 3. The algorithm: forward-backward on the bipartite chain

### 3.1 Forward-backward vs Gauss-Seidel (the D2 tradeoff)

| | forward-backward (chosen) | Gauss-Seidel (current) |
|---|---|---|
| relay | **by construction** — the forward message *is* the cumulative left-evidence, predicted through transparent nodes | emergent, lossy, **broken at zeros** (the bug) |
| exactness | **exact on a path** (our per-reference chain is a path) given `σ²_bio` ⇒ 2 inner passes | iterative; many passes |
| to fix relay | already correct | must add message-forwarding ⇒ *reimplements FB, poorly* |
| local factor | non-Gaussian (strand BB) ⇒ projected/EP update (we already do a lattice solve) | same |
| both-flank evidence | combines forward ⊗ backward (disjoint evidence) — correct | avoided one-at-a-time (the "no-blend" heuristic) |
| cost | rewrite `node_sweep` | incremental |

The decisive point: relay is the goal, and FB delivers it natively on a path. GS would have to grow a
forwarding mechanism that *is* FB. We pay a `node_sweep` rewrite once.

> **OPEN QUESTION 3 (the "never blend two boundaries onto one region" ruling).** FB forms each node's
> posterior as `forward ⊗ observation ⊗ backward` — i.e. it **does** combine both flanks. That was a firm
> GS-era ruling (memory: "both-flanks-combine was my WRONG drift"). Under proper FB the combine is **correct**:
> the forward and backward messages carry *disjoint* evidence (left sub-chain vs right sub-chain), and the
> message a node sends *excludes* what it just received from that side, so the double-count that motivated the
> ruling cannot occur. **I believe FB supersedes that ruling — please confirm we retire it under FB.**

### 3.2 The two passes (per reference chain)
Order the chain `0..K` (alternating boundary/region). Let `predict_c(belief, edge)` be the transition:
carry density `ρ_c` unchanged, inflate variance by the edge's process noise (decay precision):
`τ_c ← 1/(1/τ_c + σ²_bio,c)`; **gate**: if the channel is not continuous across the edge (`free_s` false on
either endpoint, or a strand birth/death boundary), the message is *not* forwarded (`τ→0`, i.e. no prediction)
— gDNA always continuous; `±` only within a transcript run.

```
FORWARD:   α_0 = prior (global pseudo-count, §4.2);  for i=1..K:  α_i = predict( update(α_{i-1}, obs_{i-1}), edge_{i-1→i} )
BACKWARD:  β_K = prior;                              for i=K-1..0: β_i = predict( update(β_{i+1}, obs_{i+1}), edge_{i+1→i} )
MARGINAL:  belief_i = update( combine(α_i, β_i), obs_i )      # combine = precision-weighted product per component
```

`update(prior, obs_i)` = the local observation solve (§3.4). `combine` = per-component precision-weighted
merge of the two incoming density beliefs. A zero-count node's `update` is ~identity ⇒ `α_i ≈ predict(α_{i-1})`
⇒ **transparent relay**. A high-count node's `update` is dominated by its own observation ⇒ it **absorbs**.

### 3.3 Boundaries are first-class state-space nodes
Region↔boundary alternation is unchanged. A boundary's *observation* = its crossing counts + its one-sided,
motif-stranded spliced (a hard sense/antisense RNA pseudo-count, §4.2). A boundary is where a strand process is
**born/killed**: a TSS/TES boundary sets `free_s` for the transition on its far side, so `±` does not relay
past a transcript end while gDNA does. (All of this data is already in `NodeGeometry`/`NodeStatics` — no new
accumulator deposits, no C++; this remains the linear-genomic-chain case, distinct from the deferred
transcript-SJ graph.)

### 3.4 The non-Gaussian local update (the one approximation)
The observation factor is the strand Beta-Binomial mixture on the simplex — non-Gaussian, and degenerate in
`f_g`. We keep the existing **lattice** evaluation for `update`: evaluate, over the `(f_+, f_-, f_g)` lattice,
`strand_BB_loglik + incoming-density-priors(ρ_c=f_c·M/E_c) + global + Jeffreys + spliced`, then **project** the
posterior marginal back to a per-component `(ρ_c, τ_c)` by moment-matching (mean, variance of `ρ_c` under the
posterior). This is assumed-density filtering — the standard FB treatment of a non-Gaussian observation, and
exactly what `_solve` already computes (it reads `_fg_median`, `_axis_mean`, `_fg_var`); we repackage those as
the projected message. **The `(M/E)²` "Jacobian" is not a free amplifier here** — it is the honest statement
that a confident density on a high-count node implies a sharp fraction; over-confidence is bounded by capping
the *incoming* precision at the source's effective count (§4.2, the `N_pseudo` cap), not by a floor.

### 3.5 Outer loop (var~mean refit)
`σ²_bio` (process noise) is refit each outer pass on the current beliefs (as today). So:
`init → [fit σ²_bio per component → FORWARD → BACKWARD → readout] × until Δ<tol`. Each inner FB is exact given
`σ²_bio`; the outer loop handles the `σ²_bio` dependence. Expect few outer iterations (the inner solve is no
longer itself iterative).

> **OPEN QUESTION 4 (process-noise axis & fit set).** `σ²_bio,c(ρ)` is the per-step density drift for component
> `c`, fit by `MonotoneVarMean.fit_offset`. Today gDNA & RNA var~means are fit on adjacent-observable node
> pairs. Under FB the natural training datum is the **adjacent-node density difference** `(ρ_{i+1}−ρ_i)²` on
> nodes where both are *observed* (count>0), indexed by mean density. Confirm: per-component, difference-based,
> observed-pairs-only — and how to handle references/chains with too few observed pairs (fall back to the
> pooled `σ²_bio`).

## 4. Local factors in count space

### 4.1 Strand likelihood — unchanged
The Beta-Binomial strand mixture stays exactly as in `simplex._mixture_strand_loglik` (it is the data
likelihood; both overdispersions applied symmetrically). It is the observation factor in §3.4.

### 4.2 Priors as pseudo-counts / density-beliefs
- **Incoming relay messages**: per-component `(ρ_src, τ)` density-Gaussians (§3.2), entering the lattice
  `update` as `−½·τ_c·(f_c·M/E_c − ρ_src)²`. The incoming precision is **capped at the source's effective
  count** `N_pseudo = ρ²/(σ²_bio + (n+a)/E²)` (= `count_space_solver_design.md` §3, corrected) so a deep node
  cannot steamroll the strand.
- **Global gDNA prior**: a per-node gDNA pseudo-count at the global-mean density `ρ̄` with precision `ν_global`
  in **count units** (the between-self-solvable-node spread, §3's `fit_offset`). This is the FB chain's
  *boundary prior* `α_0`/`β_K` and the per-node fallback. (Carries forward the corrected, strong-where-it-should-
  be-strong global from `ambig_intrinsic_solve_fix.md`.)
- **Spliced RNA**: observed spliced reads are strand-known RNA → inject as a **hard** `+`/`−` density-belief
  (high precision) on the boundary's exon side. Replaces the soft clipped-Gaussian spliced floor.
- **Jeffreys**: a `½`-pseudo-count per *admissible* component (`free_s` true) — the uninformative baseline that
  keeps the lattice posterior proper when a node is otherwise unconstrained.

A zero channel (`free_s` false) contributes nothing and is locked at 0 — recovering the Dirichlet "drops out at
α=0" property that OQ1 worried about losing, without the Dirichlet's density/precision conflation.

### 4.3 Readout (EM interface unchanged)
At convergence, per region: `f_g = n_g / (n_g + n_+ + n_-)` (normalize counts — the *only* place a fraction is
formed). `chain_region_deconv` / `chain_boundary_side_deconv` consume `(n_c)` instead of `(f_c, mass)`;
`priors.assemble_priors` is **unchanged** (still the two per-locus Dirichlet scalars `α_rna_add`,
`α_gdna_add`). `gdna_density_global` likewise reads densities `n_g/E`.

## 5. Module-by-module work

| module | change |
|---|---|
| `bp_solver.py` | `NodeBelief` → `(n_c, τ_c)` per component. `node_sweep` → `forward_backward` (§3.2): `predict`/`update`/`combine` + the outer `σ²_bio` loop. Delete `has_msg_nbr`/intrinsic-solve special-case (FB makes it unconditional — *every* node gets `update`, isolated or not). `node_densities` becomes the `predict` transition. |
| `simplex_sweep.py` | `_local_loglik` → count-space priors (§4.2); `_solve` returns the projected `(ρ_c, τ_c)` per component (moment-matched marginals) instead of writing fractions. |
| `simplex.py` | `_mixture_strand_loglik` unchanged. Lattice primitives unchanged. |
| `variance_model.py` | `MonotoneVarMean` reused as **process noise** `σ²_bio,c`; fit on adjacent observed-pair density differences (OQ4). |
| `init_beliefs` | signature-binary init now emits `(n_c, τ_c)`: intergenic → gDNA pseudo-count at `ρ̄`, `τ` low; single-strand → locked off-strand (`τ=∞,n=0`); AMBIG → uninformative (`τ→0`), **no** gDNA pin (the smoother fills it, not the init). |
| `chain_region_deconv` / `chain_boundary_side_deconv` | consume `(n_c)`; emit the same `NodeDeconv` mass split. |
| `derive.py`, `priors.py`, `result.py` | density readout from `n_c/E`; EM interface unchanged. |
| C++ | **none** — same payload (counts, eff-lens, adjacency). |

## 6. Benchmark scenarios + regression tests (build FIRST — §7 Phase 0)

Sim API: `rigel.sim.Scenario(genome_length, seed).add_gene(...)` + `add_gdna`/`gdna_fraction` + `CaptureConfig`
+ `ReadSimConfig(strand_specificity=...)` (pattern: `tests/scenarios/test_overlapping_antisense.py`).

**Scenario A — tiny middle exon (relay across an RNA gap), user-specified.**
```
T1 (+): exons (1000,2000), (10000,11000)            # long intron 2000–10000 overlaps T2
T2 (−): exons (5000,6000), (7000,7050), (8000,9000) # middle exon length is the SWEEP knob
+ training transcripts elsewhere (FL + var~mean fitting)
```
Sweep the middle exon `(7000, 7000+L)` over `L ∈ {1000, 300, 100, 50, 10}`. As `L` ↓ below the fragment length
the middle region goes to zero counts and **kills** the `−`-strand relay between T2's outer exons. Run
stranded **and** unstranded. **Regression assertion:** T2's outer exons stay co-attributed to `−` RNA (and the
T1/T2 overlap regions stay correctly split) as `L→10`; under the current solver they diverge once the middle
region zeroes out.

**Scenario B — gDNA-signal killers (relay across gDNA gaps), user-specified.**
```
genome length G; many single-exon transcripts of length L_k spaced along it (the "killers");
one opposite-strand transcript with a giant intron overlapping a stretch of them; uniform gDNA background.
```
Sweep killer length `L_k` down through the fragment length. **Regression assertion:** recovered gDNA density in
the killer stretch tracks the (known, uniform) background mean within tolerance as `L_k→`small; the current
solver **underestimates** (signal dies at each killer). This is the direct test of the "sparse-gDNA
underestimate" the user flagged.

**Invariants that must still hold:** `factor-1-under-uniform-gDNA`; the existing 1050-test suite; the 3
study scenarios (zero-gDNA phantom stays 0; capture AMBIG stays ≈ oracle).

> **OPEN QUESTION 5 (scenario truth & metric).** Both scenarios need a per-region oracle (the
> `dissect_regions.py` by-origin split generalizes). Define the pass/fail metric: absolute gDNA-density error
> in the killer stretch (B) and a strand-attribution / co-membership metric for the gapped transcript (A).

## 7. Phased rollout (each phase ends green + measured)

- **Phase 0 — expose & baseline.** Build Scenarios A & B (+ oracle harness). Measure the *current* solver's
  failure (RNA relay death vs `L`; gDNA underestimate vs `L_k`). Lands the head-to-head baseline the user wants.
- **Phase 1 — state migration (no behaviour change).** `NodeBelief → (n_c, τ_c)`; readout normalizes to the
  same fractions; keep the *current* GS sweep operating on the new state via thin adapters. Suite must stay
  bit-stable (goldens unchanged). De-risks the type change in isolation.
- **Phase 2 — FB skeleton, gDNA-only, no relay yet.** Replace GS with forward-backward but *gate prediction
  off* (process noise = ∞ ⇒ no relay) ⇒ must reproduce Phase 1 within tolerance. Validates the FB plumbing.
- **Phase 3 — gDNA relay on.** Enable the gDNA `predict` transition (process noise = fitted `σ²_bio,g`).
  **Gate:** Scenario B underestimate closes; `factor-1-uniform` holds; zero-gDNA phantom stays 0; capture
  unchanged.
- **Phase 4 — per-strand RNA relay + structure birth/death.** Enable `±` prediction gated by `free_s` and
  TSS/TES birth/death. **Gate:** Scenario A relay holds as `L→10`, stranded & unstranded.
- **Phase 5 — count-space priors.** Global / spliced / Jeffreys → pseudo-counts (§4.2); `N_pseudo` cap on
  incoming precision. **Gate:** capture AMBIG ≈ oracle; unstranded divergence (Bug C) improves or is unchanged.
- **Phase 6 — retire fractions + cleanup.** Remove fraction storage everywhere; delete the GS adapters and the
  `has_msg_nbr` intrinsic-solve (subsumed by FB's unconditional `update`); doc + golden regen.

## 8. Consolidated open questions
1. **Message algebra**: (density, precision) vs Dirichlet `Σα log f` — recommend the former for relay (§1, §4.2).
2. **Sampling pseudocount `a`** in `Var(ρ)=(n+a)/E²` — `1` vs Jeffreys `½` (§2.2). The one new constant.
3. **Retire the "no-blend two boundaries" ruling** under FB? (§3.1) — I believe yes; needs your confirmation.
4. **Process-noise fit**: per-component, adjacent-observed-pair density differences, mean-indexed; sparse-chain
   fallback (§3.5).
5. **Scenario oracle + pass/fail metric** for A & B (§6).
6. **`κ` (rna_sense_frac) & the two overdispersions** are still fit upstream of the chain (calibrate.py) — does
   any of that move, or does only the per-node deconv change? (I believe upstream is unchanged; confirm.)
7. **nRNA / mature-vs-nascent**: this plan deconvolves gDNA-vs-RNA per node as today; the mature/nascent split
   remains the per-locus EM's job. Confirm the relay model does not need a 4th latent process (I believe not —
   nascent is the unspliced RNA the `±` process already carries; mature rides the spliced pseudo-count).

---

## 9. Adversarial-review verdict — REVISED sequencing (supersedes §5/§7 ordering)

A 6-lens adversarial panel + a prior external reviewer + the Phase-0 measurement converge on a **re-sequenced**
plan. The architecture (count-space `(n_c,τ_c)` state + relay-as-prediction) is endorsed; the **risk ordering in
§5/§7 was backwards**, and one load-bearing justification was overclaimed. Three findings drive the revision:

- **CRITICAL — "FB exact on a path" is FALSE here.** The observation is the non-Gaussian, f_g-degenerate strand
  BB mixture (so the update is approximate ADF/EP, not exact), and the per-strand `free_s` gate splits the chain
  into **three differently-connected sub-graphs** (gDNA / +RNA / −RNA). FB is "a principled prediction-through-
  zero message schedule," **not** an exact 2-pass solver. §3.1/§3.5 exactness claims are **demoted**. This
  removes the rewrite's headline justification → FB must earn its place on *measured* impact.
- **CRITICAL — relay may be REDUNDANT with the global prior.** At a flat/AMBIG node the forward message's
  precision decays geometrically per hop (process noise), so over a run of flat nodes it →0 and the only
  surviving term is the **global prior**. "Relay across killers" ≈ "let the global reach flat nodes" — which the
  global already does on a *uniform* background (Phase-0 measured ~95–98% recovery; `ambig_intrinsic_solve_fix.md`
  §7 measured relay's footprint ~0% under capture). So **uniform Scenario B cannot prove relay is needed.**
- **CRITICAL — the outgoing `(M/E)²` jacobian into a short-flank boundary is an unbounded amplifier** the
  `N_pseudo` cap does not reach, and the binomial floor that would bound it **vanishes at the simplex walls**
  (exactly where tiny-exon boundaries sit). (Sharpens reviewer P4: the amplifier is the tiny exon's *boundary*,
  E→R finite + mass nonzero — NOT its contained region, where E=0=M vanish together.) Needs a count-space cap on
  outgoing `τ_f`, landed with the count-space priors.

### 9.1 Received-feedback adjudication (panel consensus)
P1 splicing-DAG **rejected** (spliced is a separate channel; the genomic relay carries nascent-unspliced + gDNA,
both continuous) — but promotes a real residue: **`free_s` is a signature-bit AND, not transcript-span
continuity**, so RNA relay can be wrongly severed at antisense-overlap seams / wrongly trusted across a skip.
P2 trough-mean **upheld — strongest correctness concern**: §3.4 silently swaps the robust **median** (current
f_g readout) for a **mean** message at exactly the bimodal nodes. Fix: **matched moments for the MESSAGE, median
for the READOUT — two explicitly different operators.** P3 double-count: mechanism **rejected** (the boundary
node owns spanning mass once; proper FB excludes `obs_i`), conclusion **upheld** (keep a precision discount; do
not retire no-blend on an idealized disjoint-evidence argument). P4 **confirmed + sharpened** (see above). OQ1
**(density,precision) — agree.** OQ2 **a=1** (status quo, no new constant; the E→0 stability argument is moot —
tiny nodes have no observation). OQ3 **reframe, don't retire verbatim** (keep no-blend through GS; retire only
after FB lands + a high-spliced-junction stability test). OQ4 **agree** — fit `σ²_bio` on **RAW** per-node
densities, not relay-smoothed posteriors; add an outer-loop stability gate. OQ7 **keep 3 processes (unanimous).**

### 9.2 The decision: count-space on Gauss-Seidel FIRST; FB gated on measured relay
Count-space (state + priors) and relay (the message *schedule*) are **orthogonal axes**. Bank the count-space
value behind a **reversible GS change first**; build FB **last and only if** a measured relay impact beats both
the global prior and a cheap windowed-pooled baseline. This **merges the old Phase-1 (state) + Phase-5 (priors)**
into a single count-space-on-GS step that delivers full value and **can STOP there** (= the user's "count-space
now, defer relay" deliverable, possibly permanently if relay does not earn its place).

**Revised phase order:**
- **Pre-0 — reconcile + Step-1 removal.** Resolve the contradiction with `CALIBRATION_ARCHITECTURE.md §8`
  (AUTHORITATIVE: stateless sum-product until Step 3). Do §8's **pure-removal** first (count prior, RNA prior,
  q_rna, I₀-in-sweep, dead density-variance) — it directly kills the intron bug and is justified regardless of FB.
- **Phase 0 — scenarios that DISCRIMINATE relay from the global.** Split B into **B-uniform** (global *should*
  recover it; passing is NOT evidence relay is needed) and **B-step** (anchored high-gDNA region beside a killer
  run whose TRUE level = the anchor ≠ global mean; only genuine relay recovers it). **Ablate the global** (or
  hold at truth) in the relay gate. Add **Scenario C** (overlapping-antisense high-nascent intron, zero local
  gDNA, nonzero distal gDNA → does gDNA relay + global *launder* phantom gDNA at AMBIG introns?). Add a stand-
  alone **bimodal-node projection unit test**. Prototype a **windowed-pooled global_μ** as the relay control arm.
- **Phase 1 (merged) — count-space on the existing GS sweep.** (i) per-component `(n_c,τ_c)` state incl.
  `τ_pos/τ_neg` plumbed; (ii) factor `refit_pass(snapshot)→PassPriors` as a sweep-agnostic hook (**the door
  hinge**); (iii) count-space priors (global/spliced/Jeffreys as pseudo-counts; `N_pseudo` cap; the **outgoing-
  τ_f count-space cap**); (iv) keep `_fg_median` readout; (v) gate the self density-prior to ~0 precision when
  `E_c` < a fragment-length floor (tiny-exon transparency). **GATE:** fraction-space equivalence within tight
  `rtol` (NOT bit-stable goldens); factor-1-uniform holds; the `has_msg_nbr` phantom numbers (stranded-0 = 0,
  unstranded-0 ≈ 34,886, capture ≈ 359,766) reproduce as a **hard CI test**; capture overshoot bounded.
- **Phase 2 — FB GO/NO-GO.** Build FB **only if** B-step (global ablated) shows relay recovering the anchored
  level where global-only AND the windowed baseline cannot, by a margin that matters in downstream quant error.
  Else "defer relay" → **drop relay**, ship the count-space migration.
- **Phases 3–5 (if GO) — FB + RNA relay.** "Gated off" ≡ process-noise=∞ AND messages carry the global only.
  Matched-moment messages, median readout. Specify per-reference chain **terminal boundary conditions** + the
  **linear-vs-log (lognormal) density carrier**. **Redefine `free_s` as genomic transcript-span continuity**
  before RNA relay; gate on the antisense-overlap + alternatively-spliced Scenario-A variants. Compute the
  global mean from RAW pre-relay self-solvable observations (preserve the firewall; add Bug-C as a non-regression
  gate). Do NOT delete `has_msg_nbr`/GS adapters until FB-equivalence CI passes with them removed.
- **Phase 6 — perf/determinism/ledger.** Per-chain FB memory at genome scale (the ~31GB/chain OOM risk);
  `σ²_bio` cross-reference pooling policy; determinism under the new schedule + OpenMP; the **full new-constant
  ledger** (a, the τ_f cap multiple, the `E_c` floor, the continuity-weight threshold, the σ²_bio damping γ, the
  sparse-chain pooling threshold) through the no-magic-numbers review — not just `a`.

### 9.3 Door-open invariants (so the later FB swap is NOT a second substrate rewrite)
1. per-component `(n_c,τ_c)` for **all three** from day one (incl. `τ_pos/τ_neg`, currently unused) + a round-
   trip test so a later "simplify" cannot delete the per-strand precision channel.
2. `refit_pass(frozen_snapshot)→PassPriors` factored out — **FB must inherit this contract unchanged.**
3. message stays `(density, precision)`; discount **decays precision only, never moves the mean** (the OQ1
   invariant relay needs); never regress to Dirichlet `Σα·log f`.
4. `σ²_bio` = process noise on the **frozen previous-pass snapshot**, per-component, from **raw** observable
   densities — the exact quantity FB's `predict()` reuses.
5. keep the **no-blend** ruling and **`has_msg_nbr`** intrinsic-solve through the GS step; retire only after FB
   proves equivalence on a hard CI test.
6. message-projection (matched-moment **mean**) and readout-projection (robust **median**) are **different
   operators** — prevents the median→normalized-mean +8.7pt regression sneaking in under "bit-stable".
7. cap outgoing `τ_f` in **count space** now (the jacobian-into-short-flank amplifier) — relay leans on it harder.

### 9.4 Completeness gaps to resolve before code (no one had raised these)
linear-vs-**log** density carrier (lognormal → positive support + scale-free `σ²_bio`); FB **terminal boundary
conditions** per reference (edge nodes get one real message + one global → systematic global-pull bias);
**genome-scale per-chain memory + σ²_bio pooling**; a **direct bimodal-node** projection test; the **spliced-vs-
unspliced precision contract** at the boundary exon side; the **global-prior circularity** (computed from beliefs,
fed back into beliefs — FB-smoothed AMBIG densities can erode the self-solvable firewall → Bug C); **determinism**
under FB schedule + parallelism; the **full new-constant ledger**.
