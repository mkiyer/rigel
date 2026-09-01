# NEXT SESSION — the coordinate thread (ruled 2026-09-01)

    ⚠ **A DEV DOC, and it is a HANDOFF.** It says where things stand, not what to build — the ranked
    list is `ROADMAP.md`. MOVE anything that settles into the permanent docs and DELETE this file.

## The ruling this session ended on (owner, 2026-09-01)

**The next build is the simplex solver's MAGNITUDE coordinate** — `ISSUES: arcsine-magnitude-coordinate`
— and it lands BEFORE the message policy, because messages deliver claims in ψ's own coordinate
(`TRAPS: off-grid-message-mode`) and a policy built on λ would be reworked. The certified-RNA leak is
PARKED (its no-leak ceiling was worth −0.6 % on the release metric at the worst condition — the
in-solve correction refused itself; two follow-ups recorded in the issue).

## Step 1 — the DISSECTION, before any derivation is trusted

Confirm the attribution at the object level: the 2026-08-31 drained-frame baseline
(`docs/dev/BASELINE_2026-08-31_drained_frame.md`) shows the solver UNDER-calling gDNA at near-pure
objects on EVERY stratum (shortfall +0.15…+0.29 at true `f_g ∈ [0.999, 1]`, 34–43 % of each stratum's
error mass; `g98 ss.99 capture-ON` is the worst in-scope condition at ~99 % under-call). The claim to
verify: that error mass sits on VERTEX-BLOCKED objects (belief pinned short of `f_g = 1` by the λ grid
and its measure), not on e.g. the capture-ON divisor (`ISSUES: capture-blind-gdna-divisor`).

    python scripts/design/worst_objects.py        # rank g98 ss.99 ON's objects by error mass
    python scripts/design/calibration_walk.py     # which STAGE introduces it (needs slot_truth — certified, drained frame)

Read the concentration curve first (concentrated ⇒ a mechanism; diffuse ⇒ systematic — here a
SYSTEMATIC shortfall in one f_g band is itself the expected signature). ⚠ `solvability_audit`'s worst
in-scope row differs (`g50 ss.50 OFF` by Σ|err|) — look at both before concluding.

## Step 2 — the DERIVATION (paper before prototype; BOTH solvers)

**The map is `src/rigel/calibration/simplex_logodds.py`'s own docstring — read it first.** Current
state: magnitude axis `λ = logit(f_g)` on a FIXED `[−L, L]` grid (`_logodds_grid`, `_DEFAULT_L`);
`f_g = σ(λ)`. Single-strand regions solve 1-D in λ; AMBIG regions marginalise on the 2-D `(λ, θ)`
grid (`_solve_ambig_logodds`). ⭐ **The TILT is already arcsine** (`θ = arcsin(τ)`, fact 3 — the
Berger–Bernardo Jacobian cancellation, landed 2026-08-24). The derivation must cover:

1. **The measure bookkeeping.** Fact 2's "no Jacobian is written" is λ-SPECIFIC (each `logP` density
   in log-rate cancels `log σ′(λ)` exactly, once per group). Under `f_g = sin²φ` that cancellation
   does not hold as written — derive the φ-analogue or write the correct Jacobian explicitly, and
   prove which on paper before any grid code.
2. **The reference in φ.** The flat measure in φ is the Jeffreys prior for the binary split — if the
   BB analysis carries over, `ψ_ref` may vanish entirely (the reference-mean graveyard moot by
   construction). Verify against the BB derivation, don't assume; the tilt's fact-3 cancellation must
   survive alongside.
3. **Vertex endpoints ON the grid.** At exact `f = 0/1` one component's rate is 0: data with counts
   there forbid the vertex through the likelihood; zero counts reach it. That is the desired
   behaviour and must come out of the derivation, never a clamp. The `g00`/near-pure behaviours on
   both panels are the falsifications.
4. **Grid spacing.** The `dlam = 0.3390` note in the module is measured-in-λ and does not transfer;
   information in φ is CONSTANT (the arcsine is the variance-stabilising transform), which suggests a
   uniform-φ grid is the derived choice — derive it, don't assert it.
5. **Messages and precisions.** Claim delivery (`TRAPS: off-grid-message-mode`), `hop_logvar`, and
   `τ_λ → τ_φ` semantics all map through the coordinate; the message laws
   (`TRAPS: zero-the-precision-with-the-value`, single-source-may-only-reduce) must survive the map.
6. **The 2-D solve.** The AMBIG `(φ, θ)` grid — memory shape (`sweep_n_grid`) and the tilt
   marginalisation.

Then: prototype OUTSIDE `src/` as a parallel solve, A/B per stratum on `calibration_vs_oracle.py`
with BOTH zero controls (and the fl-gap sign arms if anything length-adjacent moves).
⚠ `rename_identity.py --freeze` before and `--check` after — this is exactly the numerically-loaded
change it exists for.

## Cautions

* ⛔ The 2026-08-24 arcsine measurement ("equivalent-and-clean, never a panel win") was the TILT axis
  on the pre-fix tree. The magnitude-axis claim is NEW and must be measured. Expectation: the
  coordinate fixes reachability and the information pathology structurally; how much the panel moves
  is the A/B's to say.
* ⛔ The attribution can still surprise — if the dissection points at the capture-ON divisor or the
  EB anchor instead, STOP and re-rank before deriving.
* ⭐ General state: baselines re-derived 2026-08-31 in the drained frame with same-session floors
  (the snapshot beside this file — worst in-scope `g98 ss.99 ON`; relay hurts the solvable set on
  9/12 contaminated conditions; MessagePolicy ≡ silent on all 46; reseed floor ~30–2,700 transcript
  fragments per condition). The frame ruling is `DESIGN.md` §4.3; only wave-3 bank-readers still read
  pass one (`ISSUES: hygiene-ledger`). ⛔ Nothing measured before 2026-08-31 is comparable to
  current numbers.

## Parallel, cheap, no solver

`ISSUES: message-value-for-blind-slots` on the anchored twin block: is the information a blind
unstranded/AMBIG slot needs actually PRESENT in its neighbours? It decides the message policy's
ceiling before the policy is built — and the policy itself comes AFTER the coordinate.
