<!-- title: Count-space on Gauss-Seidel — execution plan -->
# Count-space on Gauss-Seidel — execution plan (the "exactly what")

**Companion to** `count_space_relay_implementation_plan.md` (the *why* — design + §9 adversarial verdict).
This doc is the **execution reference**: exactly what changes, in which file/function, with which gate, for the
**count-space-on-GS** effort. Forward-backward/relay is a **gated future phase** (§7), not built here.

**Locked decisions:** count-space on the existing Gauss-Seidel sweep first (reversible, full value); FB/relay
gated on a measured non-uniform-step impact (possibly dropped); **sampling pseudocount `a = 1`**; relay deferred,
door kept open via the §6 invariants.

**Checkpoint:** clean revert point is `main@519b123a` (shipped: intrinsic-solve fix `225d65ac` + this plan +
Phase-0 harness). All work below branches forward from there.

---

## 0. Pre-0 — doc reconciliation (NO code removal; audited 2026-06-19)

`CALIBRATION_ARCHITECTURE.md §8` prescribed a "Step-1 pure removal" (count prior, RNA-odds prior, q_rna, I₀
`gdna_strand_info_scale`, `mass_u` caps, DIRECT var~mean, dead density-variance) and "stateless sum-product
until Step 3." **Audit result: every Step-1 removal target is ALREADY GONE** — the BP-sweep rebuild did them, and
boundaries are already first-class (Step 3 done). The only `deference`/`I0` matches are the *word* "deference" in
two comments. So §8's step-ordering is **superseded by the rebuild**, and the count-space plan does **not**
contradict any live code.

**Action (cheap, no code):** add a one-paragraph note to `CALIBRATION_ARCHITECTURE.md §8` stating Step-1/Step-3
are complete (per audit) and that the per-node *count-space* state object is introduced by this plan as the
next step on the already-bipartite sweep. No removal commit needed.

## 1. The count-space model (concrete definitions)

Per chain node `i`, per component `c ∈ {pos, neg, g}`, the sweep state becomes:

```
n_{i,c}   pseudo-count (signal): expected component-c fragments at node i = ρ_{i,c} · E_{i,c}
τ_{i,c}   precision on the DENSITY ρ_{i,c}  (0 = no info; locked channel ⇒ n_c≡0, τ_c≡∞)
density   ρ_{i,c} = n_{i,c} / E_{i,c}        E_{i,g}=eff_gdna (gDNA-FL), E_{i,±}=eff_rna (RNA-FL); fixed
fraction  (READOUT ONLY) f_c = n_c / Σ_c n_c
```

**Precision sources** (the three the user named — length, count, structure), with `a = 1`:

```
sampling (density) precision:   τ_samp,c = E_c² / (n_c + a)          # Var(ρ)=(n+a)/E²  → short E ⇒ low precision
message (relay-ready) precision: τ_msg   = 1 / ( σ²_bio(ρ) + (n_src + a)/E_src² )   # biological + sampling
structure:  free_s = False  ⇒  n_c ≡ 0, τ_c ≡ ∞ (locked; emits 0 on that channel)
```

**Two precision caps** (both land in this effort — the adversarial criticals #2/#3):
- **source cap (`N_pseudo`)**: a message carries at most the source's effective count
  `N_pseudo = ρ_src² / (σ²_bio + (n_src+a)/E_src²)` of evidence — a deep node can't steamroll the strand.
- **destination cap (outgoing `τ_f`)**: clamp the outgoing fraction-precision in **count space** so the
  `(M_dst/E_dst)²` jacobian into a short-flank boundary cannot emit unbounded precision:
  `τ_f ≤ 1 / max(μ_f(1−μ_f)/M_dst, 1/M_dst²)` — i.e. the fraction can't be known better than one fragment's
  binomial resolution at that node, even where the wall-vanishing floor would otherwise let it escape.
- **tiny-exon transparency**: when `E_c` is below a fragment-length-derived floor (a region too short to
  *contain* a fragment), the node's self density-prior precision → ~0 ⇒ it is a transparent relay node by
  construction (and never triggers the degenerate `M/E` lattice evaluation).

**Operator split (adversarial P2 / MAJOR-4 — load-bearing):** message-projection and readout-projection are
**different operators** and must stay so.
- **MESSAGE** (what the sweep propagates): matched moments — mean density `E[ρ_c]` and `Var[ρ_c]` → `(n_c, τ_c)`.
- **READOUT** (what the EM prior consumes): the robust **posterior median** `f_g` (`_fg_median`), NOT the
  normalized-mean `n_g/Σn` (which re-incurs the documented +8.7pt grid-MAP skew the median fixed).

## 2. Phase-0 — scenarios that DISCRIMINATE relay from the global prior

Built in `scripts/debug/relay_scenarios.py` (extend the existing harness). Each has a by-origin oracle + a
pass metric. These are the **measurement targets**; the count-space-on-GS gates (Phase 1) and the FB GO/NO-GO
(§7) reference them.

| scenario | construction | what it proves |
|---|---|---|
| **B-uniform** | killers + uniform gDNA background | global *should* recover it; passing is **expected, NOT** evidence relay is needed (the count-space gate, relay off) |
| **B-step** | anchored **high-gDNA** region beside a killer run whose TRUE level = the anchor ≠ global mean | only genuine left-evidence relay recovers it → the **FB GO/NO-GO** |
| **C** | overlapping-antisense, **high-nascent intron**, zero local gDNA, nonzero **distal** gDNA | does gDNA relay + global *launder* phantom gDNA at AMBIG introns? (the dominant real-data error mode) |
| **bimodal-unit** | one synthetic balanced-count AMBIG node, off-chain | gates the §1 projection operator directly: forwarded message must be (low-precision, prior-centered), NOT (trough-mean, inflated-var) |
| **windowed-control** | a cheap windowed/locally-pooled `global_μ` arm | if it closes B-step, FB's marginal value ≈ 0 ⇒ "defer relay" → **drop relay** |

Add a **τ_f finiteness assertion**: sweep `L ∈ {10, 50}` with the global ablated; max per-message `τ_f` must be
finite and ≤ a small multiple of the largest observed node count (guards the §1 destination cap).

## 3. Phase 1 — count-space on the existing Gauss-Seidel sweep

Sub-phased so each lands behind its own gate. The GS walk and the directional one-at-a-time solve (the no-blend
ruling) are **unchanged**; only the state, the prior terms, and the precision arithmetic change.

### 1a — state migration (rtol-equivalent; estimators UNCHANGED)
- `bp_solver.NodeBelief`: `(f_pos,f_neg,f_g, var_pos,var_neg,var_gdna)` → per-component `(n_pos,n_neg,n_g,
  τ_pos,τ_neg,τ_g)`. **Plumb `τ_pos/τ_neg` even though GS consumes only `τ_g`** (door-open invariant 1).
- Keep the CURRENT estimators (`_fg_median` for f_g, `_axis_mean` for f±) — store them as counts
  (`n_c = f_c·M`); the message/readout values are numerically the same, only the container changes.
- **Factor `refit_pass(frozen_snapshot) → PassPriors`** out of `node_sweep` (the per-pass `fit_gdna_varmean` +
  `fit_rna_varmean` + global-prior recompute) — the **door hinge** (invariant 2). `node_sweep` calls it; the
  future `forward_backward` will call it unchanged.
- **GATE 1a:** fraction-space equivalence within tight `rtol` (≤1e-9) on the 3 study scenarios + the suite;
  readout estimator unchanged; a round-trip test that `τ_pos/τ_neg` survive storage; `refit_pass` equivalence.
  *(Not "bit-stable goldens" — the normalization reorder won't be bit-identical; gate on rtol + a reviewed,
  attributable diff if goldens regen.)*

### 1b — count-space prior terms (`simplex_sweep._local_loglik`)
- Re-express the node-class priors as **pseudo-counts**: global gDNA → a gDNA pseudo-count at `ρ̄·E` with
  precision in count units; spliced → a **hard one-sided ±** RNA pseudo-count (replaces the soft clipped floor);
  Jeffreys → `a=1` per admissible component (`_PRIOR_EPS` retained as the lattice-edge floor).
- The imputation messages enter as per-component `−½·τ_c·(f_c·M/E_c − ρ_src)²` (density-Gaussian on the lattice).
- **GATE 1b:** `factor-1-under-uniform` holds; the **`has_msg_nbr` phantom numbers reproduce as a HARD CI test**
  (stranded-0 AMBIG gDNA = 0; unstranded-0 ≈ 34,886; capture ≈ 359,766); capture overshoot bounded; suite green.

### 1c — count-space precision + the two caps (`bp_solver._message`)
- Reparameterize the message precision to the §1 count form (`σ²_bio + (n+a)/E²`), the **`N_pseudo` source cap**,
  and the **outgoing `τ_f` destination cap**; the message **decays precision only, never moves the mean**
  (invariant 3 — the relay-ready form). Keep the binomial floor as a *lower* bound on variance.
- **Tiny-exon transparency:** gate the self density-prior precision → ~0 when `E_c < fragment-length floor`.
- **GATE 1c:** the τ_f finiteness assertion (§2) passes (`L∈{10,50}`, global ablated); Scenario A demonstrates
  **length-aware precision** (a short region defers, a long one asserts); capture leak ≤ the `225d65ac` baseline;
  intron de-voting demonstrated; suite green.

### 1d — readout (`chain_region_deconv` / `chain_boundary_side_deconv` / `derive` / `priors`)
- Consume `(n_c)`; emit the same `NodeDeconv` mass split. **f_g readout = the posterior median** (carried from
  the solve), NOT `n_g/Σn`. `assemble_priors` interface (two per-locus Dirichlet scalars) **unchanged**.
- **GATE 1d:** EM-prior values within rtol of 1a where estimators match; **explicitly re-run the +8.7pt
  median→mean regression case** and confirm it does not return.

**STOP POINT.** After 1a–1d the count-space migration is complete and shippable with full value (length-aware
precision, bounded capture overshoot, count-de-voting, the caps). This is the user's "count-space now, relay
deferred" deliverable. The effort can end here permanently if §7 shows relay does not earn its place.

## 4. Test gates (consolidated — all HARD CI)
1. **`has_msg_nbr` phantom CI** (from 1b on): stranded-0 = 0, unstranded-0 ≈ 34,886, capture ≈ 359,766.
2. **factor-1-under-uniform** invariant.
3. **τ_f finiteness** (`L∈{10,50}`, global ablated).
4. **+8.7pt median→mean non-regression** (readout stays median).
5. **fraction-space rtol equivalence** at 1a.
6. full suite green (1050+); goldens regen only with a reviewed, attributable diff.

## 5. Module-by-module change map
| file | change |
|---|---|
| `bp_solver.py` | `NodeBelief`→counts+precisions; `refit_pass()` hook factored out; `_solve` returns `(n_c,τ_c)` per component + the median-`f_g` for readout; `_message` count-space precision + `N_pseudo` + outgoing-`τ_f` caps + decay-precision-only; tiny-exon transparency gate; `node_densities`→`ρ_c=n_c/E_c`. GS walk + `has_msg_nbr` intrinsic-solve + no-blend **unchanged**. |
| `simplex_sweep.py` | `_local_loglik` count-space priors (global/spliced/Jeffreys as pseudo-counts); `_solve` outputs per-component matched moments (mean+var) for the MESSAGE **and** median-`f_g` for readout. |
| `simplex.py` | `_mixture_strand_loglik` unchanged (data likelihood). |
| `chain_region_deconv` / `chain_boundary_side_deconv` / `derive.py` / `priors.py` | consume `(n_c)`; median-`f_g` readout; EM interface unchanged. |
| `variance_model.py` | `MonotoneVarMean` reused; fit on **raw** frozen-snapshot densities (OQ4); GCV-λ on the **weighted** normal equations (MINOR fix). |
| C++ | **none.** |

## 6. Door-open invariants (so a later FB swap is NOT a second substrate rewrite)
1. per-component `(n_c, τ_c)` for all three from day one (incl. `τ_pos/τ_neg`) + round-trip test.
2. `refit_pass(frozen_snapshot)→PassPriors` factored out — FB inherits it unchanged.
3. message `(density, precision)`; discount decays **precision only, never the mean**; no Dirichlet `Σα·log f`.
4. `σ²_bio` = process noise on the frozen snapshot, per-component, from **raw** densities.
5. keep no-blend + `has_msg_nbr` through this effort; retire only after FB proves equivalence.
6. message-projection (matched-moment mean) ≠ readout-projection (median) — kept as two operators.
7. outgoing `τ_f` capped in count space now.

## 7. Gated future — forward-backward / relay (NOT built here)
Build FB **only if** Phase-0 **B-step (global ablated)** shows relay recovering the anchored gDNA level where
**both** global-only and the windowed-control cannot, by a margin that matters in downstream quant error. Before
building: demote the "exact on a path" claim (false — non-Gaussian degenerate observation + per-strand `free_s`
3-sub-graph); single forward + single backward ADF schedule; specify per-reference terminal boundary conditions;
decide **linear-vs-log (lognormal)** density carrier; **redefine `free_s` as genomic transcript-span continuity**
before RNA relay; compute the global mean from **raw pre-relay** self-solvable observations (preserve the
firewall; Bug-C non-regression gate). Full detail: `count_space_relay_implementation_plan.md` §9.

## 8. New-constant ledger (no-magic-numbers review)
- `a = 1` — sampling pseudocount (one Poisson pseudo-observation; **decided**).
- tiny-exon `E_c` floor — fragment-length-derived (tie to the FL pmf, e.g. the read length or `frag_min`; **flag
  for review**, not a free scalar).
- outgoing-`τ_f` cap — derived from the node's one-fragment binomial resolution (`1/M_dst²`); **not a free
  scalar**.
- (FB-only, deferred) σ²_bio damping γ, sparse-chain pooling threshold, continuity-weight threshold.
