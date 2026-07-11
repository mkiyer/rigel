# Per-node precision state: honest message send/receive (Q1+Q2 fix)

**Status:** design note, revised 2026-06-25 after the open-question review. Scoped slice of
`count_space_relay_implementation_plan.md` (the full state-space/relay plan; its adversarial review §9.2 says:
**count-space precision on the existing Gauss-Seidel sweep first**, defer the forward-backward rewrite). This
note develops *only* the precision-state + honest message send/receive — the direct cure for the phantom-gДНК
cascade (`ambig_nascent_test_root_cause`, the boundary corruption trace).

## 0. The bug, in one sentence

A node's outgoing message carries **no information about how confident the source actually is**, so an
*unsolved* AMBIG node (which knows nothing about its own gДНК) emits a gДНК message at `N≈50` pseudo-counts and
overwrites a boundary holding 191 clean +RNA fragments. **Fix:** the message strength must include the source's
own uncertainty — an unsure node necessarily speaks quietly.

## 1. The minimal change (the whole fix in one line)

The send precision today ignores the source's own composition uncertainty:

```
today:   N_src = ρ_src² / ( σ²_bio + pois )
FIX:     N_src = ρ_src² / ( Var_own,src + σ²_bio + pois )
                  with   Var_own,src = (M_src/E_src)² · Var(f_c^src)
```

`Var_own,src` is the source's own posterior variance of its component **density** `ρ_c = f_c·M/E`, obtained by
moment-matching its lattice posterior (`Var(f_c)` = the per-component posterior variance the lattice already
computes; `_fg_var` and the axis analogues). Everything else — the count-binomial message form, the Gauss-Seidel
sweep, the readout — stays. This is the entire Q2 fix; Q1/Q3/Q5/Q7 are about getting the three variance terms
in that denominator right.

**Why this is honest and safe:**
- It is a **variance blend** (Q5): the three independent variance sources *add* — source-own + communication +
  sampling. An unsure source contributes a huge `Var_own` term ⇒ `N_src → 0` *smoothly* (no cliff, no cap). A
  confident source has `Var_own ≈ 0` ⇒ `N_src` is communication-limited (today's value).
- The `(M/E)²` here sits in the **denominator** — it can only *attenuate* a message, never amplify it. This is
  the opposite of the dishonest Jacobian we deleted (which amplified in the strength). No cap is needed because
  nothing is being amplified.
- It **respects discreteness**: the message stays the count-binomial Beta pseudo-count; `Var_own` only adjusts
  its strength `N`. No continuous Gaussian approximation of the count process.

## 2. The state: value + precision, per component

Per node `i`, per component `c ∈ {+, −, g}`:
- **value** — density `ρ_{i,c}` (the currency; count `n_{i,c} = ρ_{i,c}·E_{i,c}`).
- **precision** — `Var(f_{i,c})`, the per-component posterior variance from the lattice (equivalently
  `τ_{i,c} = 1/Var(ρ_{i,c})`; `Var=∞`/`τ=0` ⇒ no info, listens; `Var=0`/`τ=∞` ⇒ locked, ignores).

This re-introduces a per-node precision — but, unlike the `var_*` we retired, it is **per-component AND
consumed** (the retired one was dead precisely because nothing read it; this is the missing wiring). Value
representation is a code-cleanliness choice (Q6), no behaviour content.

## 3. Where a node's OWN precision comes from (length, count, structure, strand)

`Var_own,c` is moment-matched from the node's local lattice posterior, which already blends:
- **Length** — short `E_c` ⇒ the density `ρ_c=n_c/E_c` has large sampling variance ⇒ low precision even at `n=0`
  (the tiny-exon/relay case).
- **Count** — more fragments ⇒ smaller sampling variance ⇒ higher precision.
- **Structure** — a forbidden strand (`free_s` false) is locked (`Var(f)=0`, `n=0`): contributes nothing, ignores
  messages on that channel.
- **Strand resolving power** *(the load-bearing one)* — at a **balanced** node the strand BB is flat in `f_g`
  (gДНК and balanced ±RNA both give `p≈½`), so `Var(f_g)` is **large** ⇒ `Var_own,g` large ⇒ the node sends a
  weak gДНК message *and* listens. At a confident single-strand node `f_g` is sharp ⇒ `Var_own,g` small ⇒ it
  anchors and is hard to move. This falls out of the existing lattice posterior — no new machinery.

## 4. SEND — sending adds communication variance (the user's distinction)

A node sends its density, with a precision degraded by the between-node biological drift `σ²_bio` (the process
noise — "how well adjacent nodes predict each other"):

```
ρ_c^msg = ρ_c^src
Var(msg) = Var_own,src + σ²_bio,c(ρ_c) + pois        # own uncertainty + communication + sampling
N_src    = ρ_src² / Var(msg)                         # the count-binom strength (effective fragments)
```

Internal precision (`Var_own`) ≠ communication: sending *adds* `σ²_bio`. The act of communicating is itself
uncertain, and that uncertainty is independent of how sure the node is about its own value.

## 5. The overdispersion is FITTED, not floored; the pseudocount is a prior (Q1, Q3, Q7)

- **`σ²_bio` is the fitted overdispersion.** We fit it as the **excess over Poisson sampling** (`raw =
  (ρ_dst−ρ_src)²` minus the Poisson offset) — so "count overdispersion is a fact of RNA-seq" *is* the thing the
  `var~mean` already fits. There is **no arbitrary floor**; "variance never zero" is earned from (a) `σ²_bio > 0`
  wherever there is real drift and (b) `Var_own` (huge at unresolved nodes). The one requirement: fit `σ²_bio`
  **honestly** — on *raw observed* densities / gДНК-clean seeds, not on the circular belief snapshot (Q4).
- **`pois` is a conjugate prior, not a magic number.** Zero counts → zero Poisson variance is the pathology; the
  principled cure is a **Gamma prior on the rate** (Poisson-Gamma = Negative Binomial): the posterior rate after
  `n` counts over exposure `E` is `Gamma(a+n, E)`, variance `(a+n)/E²`, so at `n=0` the variance is `a/E²`
  (finite, low precision = imprecise/no-info — exactly the intent). `a` is the Gamma shape (`a=1` exponential /
  one pseudo-observation, = today's `_MSG_PSEUDOCOUNT`; `a=½` Jeffreys), a *prior*, not a floor. The NB does not
  remove the need for `a` (NB variance is also 0 at `n=0`); the prior is what makes zero-count imprecise.

## 6. RECEIVE — precision-weighted combine (already the lattice solve)

The receive mechanism is *already correct*: the destination's lattice solve combines its strand likelihood (its
own observation) with the incoming messages (priors); the posterior is precision-weighted. `Var(f_c^dst)=∞`
(unsolved) ⇒ believes the messages; small ⇒ ignores them; a locked channel (RNA− at a POS node) ⇒ ignores RNA−
messages entirely. The per-strand gate becomes a special case of precision. The only reason receive misbehaves
today is the messages it receives are over-precise (§1); fix the send and the existing receive correctly
down-weights the (now-weak) garbage.

## 7. Message form: count-binom (default) vs density-Gaussian — and the discreteness argument

The message is `(ρ_c, N_src)`. How it enters the destination lattice:
- **(A) count-binomial** `_binom_pseudo(f_c, μ_c, N_src)`, `μ_c = ρ_src·E_dst/M_dst` — the status quo.
  Respects the discrete, skewed, [0,1]-bounded count process. **Risk:** the Beta density diverges at the simplex
  edge — `α·log(f) → −∞` as `f→0` — a "log-wall" of height `N·μ·|log ε|` at the edge. With `N` corrected by
  `Var_own` (§1) the wall shrinks with the message (unsure source ⇒ tiny `N` ⇒ tiny wall); Phase 0 measures
  whether the residual wall at small `N` still nudges a confident node.
- **(B) density-Gaussian** `−½·τ·(f_c·M_dst/E_dst − ρ_msg)²` — smooth, no edge wall, carries (density, precision)
  cleanly for relay; **but** treats discrete counts as continuous/symmetric (poor at low counts / near 0,1) and
  re-introduces the `(M/E)²` in the strength (honest given the §1 blend + the denominator placement, but a
  reversal of the count-space decision).

**DECISION (Phase 0, measured — `phase0_precision_state_probe.py`):** use **(B) density-Gaussian + `Var_own`**.
The discreteness-favored default (A) is **refuted by the data**: even with the `Var_own` strength fix dropping
the message from `N=42.5 → 2.31`, the count-binomial still drives the boundary to `f_g=0.90` — because the
imputed `μ_c` clips to 1.0 and the `α·log(f_g)` wall at the *correct* answer `f_g=0` is `2.31·6.9 ≈ 16`
log-units, which still beats the od-flattened strand (~3). The density-Gaussian at the same `N=2.31` → `f_g=0.017`
(clean). So **both changes are necessary**: `Var_own` (weakens the message) *and* the density-Gaussian (removes
the wall) — neither alone suffices (density-Gaussian at the old `N=42.5` also corrupts → 0.98).

**On discreteness (the concern, resolved):** the message is a *prior* — a smooth belief about a *continuous*
density — not a likelihood on counts. The discreteness lives where it belongs: in the **strength `N`**
(count-derived) and in the **observation** (the strand Beta-Binomial, untouched). The density-Gaussian prior is
evaluated only on the bounded simplex lattice (no negative-density / unbounded-tail pathology). And crucially,
**discreteness ⟹ the wall**: the only discrete-respecting message form is the Beta/binomial pseudo-count, whose
log-density *must* diverge at the simplex edge — exactly where the zero-gДНК answer sits. So a discrete message
form is structurally incompatible with not-corrupting-the-edge; the continuous Gaussian *prior* is the right
choice, with discreteness preserved in `N` and the likelihood.

**On the `(M/E)²` (now honest):** with the density-Gaussian, the curvature on `f_g` is `τ_density·(M_dst/E_dst)²`
— the honest density→fraction Jacobian — but `τ_density = N_src/ρ²` is already humbled by `Var_own,src` (§1), so
it is bounded by the source's true confidence, never a free amplifier. No cap needed (Phase 0 confirms: the
boundary's strand dominates the now-weak message).

## 8. Gauss-Seidel vs forward-backward (Q8) — and why we stay on GS

- **Gauss-Seidel (current):** iterative relaxation — sweep every node repeatedly (L→R, R→L, …), each update using
  neighbours' *current* values, until convergence. Information crawls one node per sweep; approximate; can
  oscillate (the boundary swung 0.95↔0.07 across passes).
- **Forward-backward:** a *two-pass exact* method for a chain — a forward pass accumulates *all* left-evidence
  into one running message, a backward pass all right-evidence, then each node = combine(left, own, right) (the
  Kalman smoother). Exact in two passes for a Gaussian chain; carries cumulative evidence *through* a zero-count
  node (relay) automatically.
- **Why GS for this fix:** FB is *not* exact here (our observation is the non-Gaussian strand mixture, and the
  per-strand gate splits the chain into differently-connected sub-graphs), and relay may be redundant with the
  global prior (plan §9). The precision fix (§1) works on GS as-is, and once messages are honest the
  end-anchored L→R/R→L sweep behaves (confident ends dominate, unsure middle stays quiet — no oscillation). FB
  is a later, separately-justified rewrite, gated on a *measured* relay need.

## 9. Worked scenarios

- **Boundary (2:191, single-strand +).** Self-solve: `Var_own,g` small (sharp `f_g≈0`). Receives the AMBIG's
  gДНК message: the AMBIG is unsolved ⇒ `Var_own,g(AMBIG)` huge ⇒ `N_src→~2` ⇒ the boundary **ignores it**,
  stays `f_g≈0`. (Today `N≈50` overwrites it → 0.68.) Sends a strong + message.
- **AMBIG (583:637, balanced).** `Var_own,g` huge ⇒ sends weak gДНК (no cascade) and listens. Its left boundary
  sends strong **+**, its right boundary strong **−** ⇒ `f_+,f_-` set from the flanks ⇒ `f_g = 1−f_+−f_- → ~0`.
- **Single-strand exon (5:622).** Like the boundary — high own precision, ignores garbage, anchors, relays +.
- **Tiny middle exon (relay).** Short `E` ⇒ `Var_own` huge even at a few counts ⇒ listens/relays.
- **Unstranded (κ≈½).** Strand flat in `f_g` everywhere ⇒ `Var_own,g` huge everywhere ⇒ all nodes defer to the
  global + the motif-stranded spliced (the unstranded lever).
- **Capture (gДНК-enriched AMBIG).** Legit gДНК messages come from **confident** gДНК-rich single-strand flanks
  (small `Var_own`) ⇒ they survive ⇒ AMBIG exons correctly read high gДНК. ⚠ must verify recovery preserved.

## 10. Implementation plan (meticulous, phased; GS, count-space)

Each phase ends green (full suite) + a measured gate. No new constant ships without the no-magic-numbers review.

### Phase 0 — isolate & lock the mechanism — DONE ✅ (decision: density-Gaussian)
Prototype `scripts/debug/phase0_precision_state_probe.py` (real AMBIG-locus geometry). **Measured:** Var(f_g)
AMBIG 0.060 vs boundary 0.024 (only 2.5×); but `Var_own=(M/E)²·Var(f_g)` AMBIG 0.16 vs boundary 0.0036 (44× —
the (M/E)² is the separator); message `N` 42.5→2.31 with `Var_own`; boundary `f_g` at N=2.31: count-binom 0.90
(wall persists) vs **density-Gaussian 0.017 (clean)** ⇒ density-Gaussian required. The checks below become the
Phase-1 unit-test assertions:
1. **`Var(f_c)` from the lattice** for canonical nodes — assert balanced-AMBIG `Var(f_g)` ≫ confident-single-
   strand `Var(f_g)` (the strand-resolving-power claim).
2. **`Var_own = (M/E)²·Var(f_g)`** for the AMBIG (source) vs the boundary (dest) — assert AMBIG ≫ boundary.
3. **Message strength** `N_OLD = ρ²/(σ²_bio+pois)` vs `N_NEW = ρ²/(Var_own+σ²_bio+pois)` for the AMBIG→boundary
   edge — assert `N_NEW ≪ N_OLD` (≈50 → small).
4. **Apply both** to the boundary's lattice solve, **both message forms** (count-binom, density-Gaussian) — report
   `f_g`. **Decision gate:** does count-binom + `N_NEW` keep `f_g≈0` (residual wall harmless), or is the
   density-Gaussian required? This picks the form empirically (§7).
5. Output a crisp table; these become the Phase-1 unit-test assertions.

### Phase 1 — precision state, computed & stored (no behaviour change) — DONE ✅
Implemented: `NodeBelief` + `NodeDeconv` carry per-component `(var_pos, var_neg, var_gdna)` = `Var(f_c)`;
`_solve_nodes` and `node_sweep._solve` moment-match them (`_fg_var` per axis); `_type_belief` sets the init
convention (forbidden strand `var=0`; admissible-unsolved `inf`; G2 the strand-solve variance; G1/G3 default).
Unit tests added (`test_precision_state_strand_resolution` + the init-convention asserts). **Gate met:** full
suite 1059 pass / 1 (the pre-existing `test_ambig`, unchanged), goldens bit-stable **without** `--update-golden`.
- Add per-component precision to `NodeBelief`: `Var(f_+), Var(f_-), Var(f_g)` (re-introduced, now consumed).
- After each `_solve`, moment-match the lattice posterior → per-component `(ρ_c, Var(f_c))` (the message
  projection uses **mean/variance**; the readout keeps the robust **median** — two distinct operators, plan P2).
- Init: signature-pin (forbidden strand `Var(f)=0,n=0`; admissible-unsolved `Var(f)=∞`; AMBIG gДНК `Var=∞`, no
  confident pin). Value init irrelevant.
- **Gate:** goldens bit-stable (compute+store only; messages still use the old `N`).

### Phase 2 — honest send (the core fix) — DONE ✅ (mechanism); benchmark pending (Phase 4)
- `_message` blends `Var_own,src = (M_src/E_src)²·Var(f_c^src)` into the variance, returns a density-Gaussian
  `(mode, prec)`, `prec = (M_dst/E_dst)²/(Var_own + σ²_bio + pois)`; `_local_loglik` applies `−½·prec·(f_c−mode)²`
  (the count-binom message + log-wall is gone; the global prior stays count-binom). node_sweep computes
  `Var_own` per message from the source's belief variance. Message-form tests updated.
- **Gate met (mechanism):** the toy AMBIG locus — POS exon 0.20→**0.000**, AMBIG 0.367→**0.033** (<0.08),
  NEG exon 0.30→**0.000**, boundary 5000 f_g 0.683→**0.000** (relays clean +), boundary 6000 0.217→**0.000**
  (relays clean −); `test_ambig_no_false_gdna_from_nascent` now **passes**; no oscillation. Full suite 1052
  pass / 8 golden-shift (behaviour changed — NOT yet regenerated, pending the Phase-4 net-flow gate).
- **Prior-level benchmark (`enrichment_phase3_check`):** flagship gdna300/ss0.99/cap-on AMBIG net −4K (under) →
  **+17,898 (over)**, still within the ±40K recovered band, total contained gDNA net +1,775; zero-gДНК/cap-on
  small phantom +1,886; the no-gДНК+nascent+unstranded confound (~191K) is the pre-existing accepted limitation.
  The capture flip (under→over) is the AMBIG now listening to its enriched flanks — needs the **post-EM
  net-flow benchmark** (Phase 4) to adjudicate before goldens are regenerated.

### Phase 3 — honest σ²_bio (Q3/Q4)
- Audit/ensure both channels' `var~mean` is fit on **raw observed** densities/seeds, not the circular belief
  snapshot. Gamma-prior `a` for zero-count (= status-quo `_MSG_PSEUDOCOUNT`, named as the prior).
- **Gate:** zero-gДНК phantom stays ≈0 *for the right reason*; capture AMBIG ≈ oracle (the `Var_own` term must
  not silence *legitimate* gДНК messages from confident enriched flanks); circularity check (§Q4).

### Phase 4 — full benchmark + ship
- `quick_3to1_5mb` 16-condition (esp. unstranded), the capture flagship, factor-1-under-uniform, full suite +
  goldens regen. Commit per-phase.

## 11. Open questions
- **Resolved this review:** Q1/Q5 — blend not cap; Q3 — `σ²_bio` is the fitted overdispersion, no floor; Q7 —
  pseudocount = Gamma prior (`a=1`); Q8 — stay on GS; Q6 — value representation is a code choice.
- **Resolved by Phase 0:** message form = **density-Gaussian + `Var_own`** (count-binom's log-wall corrupts even
  at the corrected small `N`; discreteness lives in `N` + the strand likelihood, not the prior). Both the
  `Var_own` strength fix and the density-Gaussian form are necessary; neither alone suffices.
- **Still decide by testing:** Q2 — message = matched mean/variance vs readout = median (validate on the
  bimodal-node test in Phase 1); Q4 — does `Var_own` alone break the `σ²_bio` circularity on the RNA channel, or
  is the raw-fit firewall needed (measure in Phase 3)?
