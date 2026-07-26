# Calibration ROADMAP — status, architecture, and the path to production

**This is the single entry point for calibration work. Read it first.** Last updated: 2026-07-25.

> **⭐ ORDER OF WORK (owner, 2026-07-25): the pass-0 SOLVER must be CORRECT before the gDNA hyperprior fit.**
> The **single-strand × capture study is DONE** — see `SESSION_2026_07_25_HANDOFF_8.md`.
> Its answer: the 10× capture degradation is 77–92 % on EXONS and is a message **MODE** defect with an exact
> mechanism — the grafted junction flux `ρ_μ` is a spliced measurement already in the destination exon's
> frame, but it is ratioed against the *boundary's* gDNA density, and since the reframe `r` cancels from the
> delivered share (verified to `1.8e-15`) the graft edge **never reframes the gDNA at all**. Under capture
> that step is 6.1–6.8× (1.03× without capture). Fixed by **M8** (`graft_frame_logvar`), which prices the
> un-cancelled step as a variance `(log r)²` on the grafted component — derived, MC-validated, A/B-won:
> **0.0926 → 0.0900 (refit=0), 0.0779 → 0.0700 (refit=1)**.
>
> **The AMBIG study is DONE — `SESSION_2026_07_25_HANDOFF_9.md` (the LIVE handoff); read its §0b FIRST, it
> CORRECTS the rest of the file.** What stands: AMBIG is **31.6 %** of suite error mass (the "≈50 %" figure is
> the stranded-only census), `τ_own = 0` on every AMBIG node, and **messages roughly double AMBIG error**
> (self 0.090 → solved 0.144–0.192) with the λ composition message the dominant defect (bit-exact ablation
> 0.1444 → 0.0986). Two candidate fixes were A/B'd and both were ~neutral; both reverted.
>
> **What was RETRACTED (owner, 2026-07-25).** The plug-in `d̂ = (p_obs−½)/(κ−½)` is **inadmissible**: κ is
> fitted (posterior `Beta(n_same+1, n_opp+1)`), and on unstranded data `|κ̂−½|` is *smaller than its own
> posterior sd* while sitting in a denominator — med `|d̂|` up to 118, 80–99 % of mass with `|d̂| > 1`. The
> derived replacement (MC-validated, M9a/M9b) never divides: the strand likelihood depends on `f_g` and `κ`
> **only** through `A = (1−f_g)·|κ−½|`, so `κ → ½` flattens it for every `f_g` and "no information" is a
> **limit**, not a guard. Its sharp-limit closed form is `p(f_g) ∝ 1/√((1−f_g)²−d²)` — an *integrable spike*
> at the bound, not a hard wall. And the shipped 2-D solve was **never degenerate**: it already computes this
> marginal. **M9d shows the residual AMBIG bias (−0.10 to −0.15 at τ=1) is the TILT PRIOR, not discarded
> information**, so the "85–91 % of the AMBIG prize" figure was an offline substitution with an inadmissible
> estimator, not an achievable prior-free target.
>
> **NEXT: the TILT MESSAGE.** What converts the constraint into an estimate is knowing τ ≈ ±1, and the
> prior-free source is the neighbours' per-strand RNA imputation — already carried as `theta_imp`, and
> measured **nearly inert** (0.1444 → 0.1417 when ablated).
>
> **⭐ CURRENT (2026-07-26): the composition peel, the MESSAGE PACKET and **P1d** are LANDED; the next task
> is P1e.** Suite **0.0865 (refit=0) / 0.0676 (refit=1)** — deliberately up from 0.0855 / 0.0671, bought
> against a 36.6 % cut in confidently-wrong mass. **Read `SESSION_2026_07_26_HANDOFF_14.md` — it is the LIVE
> handoff and has the run instructions.**
>
> **The frame for everything that remains:** the error MASS is close to its prior-free ceiling (HANDOFF_10 §3:
> 67–75 % of it is premise-limited, and ×10 on every precision moves nothing), but the **CALIBRATION of that
> error is not**. Confidently-wrong mass fell 14.4 % → **9.9 %** this session (1,777,658 → 1,186,552 reads),
> yet `z2` inside the confident quartile is still **8.4 on single-strand exons, 92.1 on AMBIG exons, 15.5
> overall** — introns (1.5) are the only honest class. A node at `z2 = 9` gets nine times too much weight in a
> hyperprior fit that trusts declared precision. **So pass-0 is nearly done on mwae and is not done on trust,
> and the remaining work is chosen for its effect on `z2`.**
>
> What landed: **M11 `residual_level`** (the seam's RNA level from the node's own mass closed against an
> imputed gDNA density — the source that resolves "there is no LEVEL at the seams that matter"), the
> **three-way level FUSE in LINEAR space** (a log fuse cannot represent "confidently near zero"), the
> `mrna_active` structural gate **removed**, **M12 `conservation_rescale`** (law kept, path superseded), and
> **the MESSAGE PACKET** — λ and θ carried explicitly and fused by their **own** precisions rather than read
> back off the density fuse. The packet alone drove 0.0889 → 0.0855 and cut AMBIG-exon confidently-wrong 45 %.
>
> **⭐ P1d IS LANDED (2026-07-26) — and its mechanism was NOT what the term was named for.**
> `graft_premise_logvar` — **ONE library-level MoM scalar** applied to every graft edge, default ON,
> `RIGEL_GLV_OFF=1` ablates it. Confidently-wrong **1,200,951 → 762,000 (−36.6 %)**, `z2` ALL
> **15.55 → 8.98**, exon-single **8.60 → 2.46**, exon-AMBIG **92.10 → 64.39**, for **+1.2 % error mass**
> (0.0855 → 0.0865 / 0.0671 → 0.0676; 9 better / 6 worse / 17 flat, 8/4/20 at refit=1).
> ⚠ A **per-edge** variant using each exon's own two-seam gap landed first and was **reverted the same day**:
> a variance from one pair is a χ²₁ (CV = √2), so it is mostly noise, and it made the LEFT seam's message
> carry a variance built from the RIGHT seam's counts — a BP violation. Pooled-only is better on every axis.
>
> **"Extrapolation" is refuted as the mechanism.** The premise residual is FLAT in the extrapolation ratio
> (1.13× over a 6.7× range) and flat in exon length. What the graft actually gets wrong is that its claim is a
> **LOWER BOUND used as an equality** — `ρ_R(exon) ≥ ρ_ν(B) + ρ_μ(B)`, asserted as `φ ≡ 1` to 2.2e-16 — and it
> fails hardest where RNA does not flow through the seam. **A single structural bit, "does this boundary carry
> a transcript TERMINUS", splits `ω̂` ≥30× (1.7–1.9 vs 0.04–0.06) and 20.8 % of edges carry 71.7 % of the
> squared deviation.** The 40× junction-count gradient the term's shape was to be fitted to is a *proxy* for
> that bit: flat inside every structural class. Full record + the refutations: `variance_ledger.md` §3.0.
>
> **That un-defers §6 below.** The missing TSS/TES in the region/boundary map was deferred as "expected to be
> low-impact on real data"; it is measured at **62–72 % of the graft's premise error**. → **P1g.**
>
> ### ⚠⚠ `ω_graft` IS A DEBT, NOT A MODEL — DO NOT LET THIS BE FORGOTTEN
>
> P1d's fitted `ω_graft` **partially compensates for a FAILURE IN OUR STRUCTURAL REPRESENTATION.** The
> region/boundary map has no TSS/TES, so the solver **cannot tell a splice junction from a transcript
> terminus** — and that distinction is the whole of the effect `ω` prices: measured `ω̂` = **1.7–1.9 at
> terminus boundaries vs 0.04–0.06 at junction-only ones, a ≥30× split**, with 20.8 % of edges carrying
> 71.7 % of the error. `ω` is a single library-wide average standing in for a bimodal quantity. It works
> because over-charging a variance is cheap and under-charging is expensive — **not because it is right.**
>
> **Three consequences that must not be overlooked:**
> 1. **It is expected to be fragile on real data.** It is fitted on ~200 exons with a 30–50 % standard
>    error, on a suite that is Poisson by construction and whose isoforms are only nested truncations plus
>    exon skips. Real annotations have alternative TSS/TES *inside* exons, mutually exclusive exons and
>    retained introns — none of which this suite contains.
> 2. **`ω` MUST be re-derived as a per-class quantity the moment TSS/TES enter the region map.** Same
>    estimating equation, two (or more) scalars keyed on the structural bit, each fitted over a
>    homogeneous population. The partial-pooling code already in `graft_premise_logvar` is the plug-in point.
> 3. **Until then, treat every `ω`-dependent result as provisional**, and never quote the fitted value as
>    if it were a measured constant of the library.
>
> **⭐ THE PASS-0 / PHASE-2 BOUNDARY IS NOW MEASURED (`SESSION_2026_07_25_HANDOFF_10.md`).** Suite state of
> play: **12.56 M of 116.7 M node-attributed fragments misassigned (10.8 %)**; unstranded carries **84 %**,
> single-strand exons **61 %**. On the largest scenario, 87.6 % of the error sits where `_pin_v` cancels the
> reframe `r`, so the delivered mode is the source's composition — and neighbours differ by
> **|Δf_g| = 0.2400**. Oracle modes remove 67–75 %; scaling every precision ×10 removes NOTHING. **That
> residual is the hyperprior's, and its precision is already honest** (`errQ1conf` 4.1 %, gDNA channel
> z2 = 0.4–1.5). The Phase-2 RISK is elsewhere: **stranded capture-OFF puts 58–76 % of its error on the
> most-confident quartile, and introns 90.2 %** — and introns are where gDNA density is measured.
>
> **Status in one line:** the message-variance model is **COMPLETE** — derived, MC-validated, independently
> verified, implemented, and A/B-won. A message's precision is
> `1/(Var(log f_c^src) + 1/n_src + σ²_transfer + b̂²)`, plus **M8's `(log r)²` on the grafted spliced
> component only**: the source's earned composition+count precision, the reframe's **scale** uncertainty
> (M5 `Var(log r)`), the **DerSimonian–Laird composition-mismatch** `b̂²` (M7), and the graft's **un-cancelled
> frame step** (M8). **Best aggregate on record: 0.0900 (refit=0) / 0.0700 (refit=1)** vs the 0.1267/0.1234
> pre-fix baseline. **NOT ready to ship** (the hyperprior refit still regresses unstranded-capON), and per the
> owner's directive the hyperprior is NOT the next task — the **AMBIG tilt message** is
> (`SESSION_2026_07_25_HANDOFF_9.md` §0b).
>
> **Update 2026-07-25 (DL cliff-term session).** `(log r)²` charged the WHOLE enrichment cliff as composition
> drift, which recovered the stranded arm but over-damped extreme capture. The delivered message error splits
> exactly into a composition-SHARE mismatch plus the reframe's scale noise, so the two are now priced
> separately; `b̂²` is estimated prior-free against the node's own self-solve (a two-study random-effects
> meta-analysis) with **no tuned constant**. Its safety property is exact: **a message can out-weigh a node's own
> belief only if it agrees within `√2·σ_own`** — the governing principle as an inequality rather than a knob.
> Attribution is clean and the two effects are disjoint: deleting `(log r)²` recovers verystrong, `b̂²` recovers
> stranded, neither costs the other. Also retired the dead NPMLE-projection σ²_transfer plumbing: **no prior of
> any kind now enters message precision.**
>
> **GOVERNING PRINCIPLE (owner):** **pass-0 must be WEAK and correctable** — an over-confident message that pins
> a node wrong is worse than a weak one slightly off. Prefer the under-confident option when unsure.
>
> **Update 2026-07-24 (post-handoff session).** Retired `_scan` + all its flags/helpers (`bp_solver.py` 1871 →
> ~730 lines); extracted the per-node self-solve into `node_init.build_node_init` (the four init sources of
> `variance_model_concepts.md`, one unit test each), behavior-preserving (byte-identical to the pre-refactor
> unified path across all 32 scenarios). Goldens regenerated to the unified default. §2 below is now historical.

The only other docs that are live (everything else is in `archive/`, kept for history, NOT to be referenced):
* `CALIBRATION_ARCHITECTURE.md` — the authoritative theory (count-zero-information; the three information
  sources). Still correct; read second.
* `unified_solver_design.md` — the target solver's architecture (the reframe + ÷M_dst mode). Its **precision /
  variance sections (§8 R1–R4) are SUPERSEDED** by `variance_model_handoff.md`; the mode design stands.
* `gdna_intron_factory_design.md` — a shipped feature (the intron gDNA factory). Live.
* **`density_composition_reconciliation.md` — ⭐ THE NEXT SUBSTANTIVE TASK** (owner-directed 2026-07-25). The
  composition regime is scale-blind: relayed claims assert more reads of a component than the node sequenced in
  total (52–71 % of nodes; p99 31–288×). Holds the measured evidence, the derivation brief, the two adjacent
  modelling gaps (§4: **no TSS/TES in the region/boundary map**; **the boundary is a slope, not a cliff — three
  enrichment ratios, not one**), and the record of what was tried and rejected.
* **`PASS0_FINISH_PLAN.md` — ⭐ THE PRIORITIZED WORK LIST for finishing pass-0.** Owner-agreed; read its §0b
  first (the state of play, and the `z2` measurement that says pass-0 is NOT finished even though the error
  mass is near its ceiling). Holds P0–P6 with status, what was refuted, and why the list is ordered by
  confidently-wrong mass rather than error mass.
* **`variance_ledger.md` — ⭐ THE STANDING AUDIT of every variance term** in a message's precision: what each
  one prices, why none of them overlap (including the proof that M7's DerSimonian–Laird cannot double-count
  with anything), and the one measured gap that is still unbuilt. **Update it whenever a variance term is
  added, removed or re-scoped.**
* `weighted_rescale_design.md` — the conservation rescale (M12) and, in §9, the MESSAGE PACKET that subsumed
  it: λ and θ carried explicitly and fused by their own precisions.
* **`SESSION_2026_07_26_HANDOFF_12.md` — the boundary/composition-peel handoff.** The
  composition peel LANDED: the three questions of HANDOFF_11 answered, M11 derived + MC-validated +
  unit-tested, the level as a three-way FUSE, the `mrna_active` gate removed, the full per-condition A/B and
  the trust view, the do-not-re-run list, and the remaining deficit localized to one regime with its
  mechanism measured.
* `SESSION_2026_07_25_HANDOFF_11.md` — the previous handoff (superseded). Its §5/§6/§7 pose the three
  questions HANDOFF_12 answers; its §8 do-not-re-run list still stands.
* **`peel_by_composition_plan.md`** — the plan for steps 1–5 with pre-registered gates; §10/§11 hold the
  execution results, including two of its own assumptions being refuted.
* `SESSION_2026_07_25_HANDOFF_10.md` — the boundary evidence base. The deep
  dive on the suite's largest error scenario: the intron factory's gDNA is excluded from the measurement
  stream (a real prior-free defect, fix validated at refit=0 but regressing refit=1 — decision needed), and
  the **count-zero-information wall measured**: on 87.6 % of that error the pin cancels `r`, and neighbouring
  nodes differ in composition by |Δf_g| = 0.24, so no prior-free message can beat it.
  `..._HANDOFF_9.md` — the AMBIG
  study: the simplex-bound result, the ceiling, the validated prior-free estimator, the two neutral A/B arms
  (both reverted), and the implementation plan. `..._HANDOFF_8.md` (the single-strand × capture result — M8,
  the 8-step measurement chain, and the four-arm ablation that chose the variance over the mode fix),
  `..._HANDOFF_7.md` (the boundary-class census — still the map; its §4–§5 study is now DONE),
  `..._HANDOFF_6.md` (whose "next task is Phase 2" is WITHDRAWN), `..._HANDOFF_5.md` and `..._HANDOFF_4.md`
  are the arc. (Handoffs 1–3 are historical.)
* `message_variance_derivation.md` — the derived + MC-validated + independently-verified message-variance laws
  (M1–M5), the M6 combine finding, and the empirical results (§6). Live.
* `variance_model_concepts.md` — the owner's spec for the **initialization** phase (the four sources).
* `variance_foundation_proposal.md` — the SETTLED foundation model (approach E, the single Schur scalar `τ_λ`),
  derived + numerically verified. `variance_foundation_plan.md` v4 — the invariants + the deferred-work spec;
  `variance_foundation_critique.md` — the adversarial critique ledger.
* `variance_model_handoff.md` — the MESSAGE variance-model derivation substrate (§3-4), the NEXT task's math,
  built on the `τ_λ` foundation.

---

## 1. What calibration is, and what “pass-0” means

Calibration deconvolves each genomic node’s **unspliced** fragment mass into a composition
**(f_rna₊, f_rna₋, f_g)** — sense-RNA / antisense-RNA / gDNA. It runs in stages:

```
   PASS-0 (prior-free)  →  fit the gDNA HYPERPRIOR  →  RE-SOLVE (with the hyperprior)
   an APPROXIMATION        (from the pass-0 result)     the actual answer
```

* **Pass-0 is not a solution — it is an approximation.** It solves each node from only the two *intrinsic*
  information sources (the strand likelihood + cross-node imputation), with **no population prior**. On
  unstranded data the strand likelihood is flat (`CALIBRATION_ARCHITECTURE.md`), so pass-0 leans almost
  entirely on cross-node message propagation.
* **The pass-0 result trains a gDNA hyperprior** — the population baseline gDNA density landscape.
* **The hyperprior is then required to re-solve**, in particular to resolve **AMBIG** (opposite-strand
  overlap) nodes, which have no intrinsic gDNA/RNA signal at all (see
  `strand_likelihood_constrains_tilt_not_fg` in memory: on AMBIG nodes the strand likelihood constrains only
  the tilt, never f_g).

Everything below is **pass-0** unless stated. The hyperprior fit is a separate, also-WIP workstream (§4).

## 2. The two-solver problem — RESOLVED (historical)

`bp_solver.py` used to contain **two** solve paths — the legacy density-transfer `_scan` (default) and the
composition `_unified_solve` (flag-gated). **As of 2026-07-24 this is resolved: `_scan` and all its flags
(`RIGEL_B1B/N4A/N4B/E2`, the `_UNIFIED` gate) and helpers were deleted; the unified solver is the sole path.**
`bp_solver.py` roughly a third of its former size (1871 → ~730 lines); the per-node INITIALIZATION self-solve
now lives in `node_init.py` (`build_node_init` — the four sources, unit-tested). The unified solver still
**loses the A/B** (measured pass-0-vs-oracle mwae: unified ~**0.15** vs the legacy-with-factory baseline) — this
is **expected and accepted on this WIP branch**; the variance model (§3) is what recovers it. Nothing ships
until it does.

## 3. ✅ RESOLVED — the variance model (historical; see `message_variance_derivation.md` for the final laws)

The conceptual shift was: the old solver compared **absolute densities** between nodes; the new solver compares
**compositions**, normalizing by an **enrichment ratio** `r` that cancels hybrid-capture enrichment/depletion.
The old variance model assumed genome-wide density **uniformity**, which capture breaks — so it had to be
rebuilt for a composition transport. It now is (`message_variance_derivation.md`, laws M1–M7, every one
MC-validated in `scripts/debug/message_variance_mc.py`, which runs 0 failures end-to-end):

```
    p_message = 1 / ( Var(log f_c^src) + 1/n_src  +  σ²_transfer  +  b̂² )
                     \__ strand ___/   \_count_/    \_ SCALE _/    \_ COMPOSITION _/
```

The last two are the cross-cliff cost, and the session's central result is that they are **different objects**
that must be priced separately. The delivered message error splits EXACTLY (to machine precision) into
`log(s_c^src/s_c^dst,true) + log(r̂/r_true)` — a composition-SHARE mismatch plus the reframe's own scale noise.
`σ²_transfer = Var(log r)` (M5) prices the scale; `b̂²` (M7) prices the imputation PREMISE ("my neighbour and I
share a composition") being false. The retired `(log r)²` proxy charged the whole cliff as mismatch, which is
why it fixed the stranded arm but over-damped extreme capture, where the composition really is preserved
across a 1000× enrichment step.

`b̂²` is a population quantity we lack prior-free — but the destination has an **independent** estimate of its
own composition: its message-free self-solve. Two estimators of one quantity is a two-study random-effects
meta-analysis, so the **DerSimonian–Laird** between-study estimator supplies it with **no tuned constant**:
`b̂² = max(0, G² − v_msg − v_own)`, giving `p_eff = 1/max(v_msg, G² − v_own)`. Its three regimes fall out of
`v_own` (from the `τ_λ` foundation) with no gate, and the middle one is the safety property, exact: **a message
can out-weigh a node's own belief only if it agrees within `√2·σ_own`.** Where a node has NO composition
evidence (`τ_own = 0` — every AMBIG node, all unstranded data) `v_own = ∞` and the term switches itself off, so
messages propagate untouched exactly where they are the only information. That inertness is deliberate — and it
is also why the remaining AMBIG error is Phase 2's problem, not this term's.

## 4. ✅ THE "BLOCKER" IS STALE — re-measured 2026-07-26, and Phase 2 is much closer than this file said

This section used to read *"⛔ THE BLOCKER — production ships `calib_refit_iters=1` and **that refit regresses
unstranded capture-ON**."* **That is no longer true.** Measured at HEAD (mass-weighted, all 32 conditions,
`refit=0 → refit=1`):

| stratum | refit 0 | refit 1 | Δ | better/worse |
|---|---|---|---|---|
| ALL 32 | 0.0865 | **0.0676** | −0.0188 | 25 / 7 |
| stranded ss_0.99 | 0.0313 | **0.0293** | −0.0020 | 10 / 2 |
| **unstranded × capON** | 0.1523 | **0.1275** | **−0.0248** | 4 / 4 |
| capture OFF | 0.0475 | **0.0319** | −0.0156 | 14 / 0 |
| verystrong | 0.1844 | **0.1264** | −0.0579 | 3 / 1 |
| gdna_none | 0.0941 | **0.0029** | −0.0912 | 9 / 0 |

**The refit improves every stratum, including the one it was recorded as regressing** — and it already did so
before P1d (pre-P1d unstranded × capON 0.1547 → 0.1001). The regression was real when it was written and has
been fixed by the intervening work (the message packet, M11, the relay pin); the note simply outlived it.

**And the second reason Phase 2 looked blocked is also wrong: AMBIG's over-confidence does not reach the
hyperprior.** `_fit_gdna_hyperprior` already selects **REGION nodes that are single-strand or structural
gDNA — no AMBIG, no boundaries** ("non-circular", `calibrate.py`). So the substrate the prior actually trains
on is exactly the classes P1d made honest:

| population | confidently-wrong | **z2 \| Q1** |
|---|---|---|
| **the hyperprior's training substrate** (region × single-strand) | 431,708 | **2.26** |
| … its exon half | 246,391 | 2.46 |
| … its intron half | 185,317 | 1.50 |
| EXCLUDED from training (AMBIG + boundary) | 330,292 | **42.04** |

**All of the remaining two-orders-of-magnitude over-confidence sits in the population the hyperprior never
sees.** So P3 (AMBIG) is a pass-0 *quality* item for the re-solve, **not** a Phase-2 blocker.

### 4b. ⭐ THE RE-SOLVE ALREADY RESETS — verified 2026-07-26, and it sharpens the objective

`calibrate.py:419-420`, inside the Phase-2 refit loop:

```python
    belief = _init_belief()            # ← FULL reset to the signature-binary default
    belief = _sweep(gdna_hyperprior)   # ← re-solve from scratch, WITH the prior
```

Nothing from pass-0 survives into the re-solve except **the fitted hyperprior itself**. The owner's option
"after fitting the hyperprior, reset all nodes and re-solve from scratch" is therefore **already what
ships** — the failure mode "an over-confident node won't budge when the prior arrives" **cannot occur through
belief inheritance.**

**What that means for pass-0's objective — and the owner's ruling (2026-07-26):**

> *"Honest precision is critical. We may adopt an iterative strategy. We may decide to use a
> precision-weighted gDNA hyperprior fit. So honest precision can be beneficial in many ways. If we push
> towards a combination of ACCURACY with honest precision, we will prevail."*

So **both** metrics are live, and they are not in competition:

* **ACCURACY on the fit substrate** (`REGION & (single | gonly)`) is what the hyperprior is built from, and
  it is the number that was NOT being tracked. Suite-wide mwae is a proxy for it, not a substitute — the
  substrate is a *subset*, and a change can move the two in opposite directions.
* **Honest precision** has three live consumers, so it is not merely instrumental: (i) it sets pass-0's own
  fused MODE, hence the substrate's accuracy; (ii) the non-additive hyperprior fit is **precision-weighted**
  (`var_gdna` weights the NPMLE — `_fit_gdna_hyperprior`); (iii) an iterative re-solve strategy remains open.

**Do NOT relitigate P1d or P4 on the grounds that the re-solve resets.** Both remain landed.

⚠ **What has NOT been re-tested at HEAD** is the specific committed Phase-2 step — feeding the hyperprior's
own λ-curvature into DL's `v_own` (~6 lines, `SESSION_2026_07_25_HANDOFF_6.md` §3). At its original HEAD it
improved stranded / verystrong / capture-off and regressed unstranded-capON 0.1702 → 0.2177. Given both
premises above have changed, **that measurement is the single highest-value next experiment** → **P2a**.

## 5. What SHIPPED (is correct and on by default) vs WIP

**Shipped / correct (default path):**
* The **gDNA intron factory** — `intron_factory=True` (default, changed this session). Deconvolves introns
  against the intergenic background NegBinom and now **carries its derived precision** as composition evidence
  (`_lambda_factor_precision`), so a factory-solved intron can actually propagate. Pass-0 vs oracle over the
  32-scenario suite: **mwae 0.1361 → 0.0949, corr 0.688 → 0.736**, 20 better / 1 worse / 11 flat; intron mwae
  0.1781 → 0.0117; every stranded scenario better-or-flat.

* The **message-variance model** (§3, M1–M7) — the sole message-precision law, no flags, no prior input.
  Pass-0-vs-oracle over the 32-scenario suite: **0.1267 → 0.0969 (refit=0), 0.1234 → 0.0828 (refit=1)**.

**Work in progress (NOT ready to ship):**
* The **gDNA hyperprior refit** (§4) — the blocker.
* **AMBIG nodes** — ⚠ **refined by `SESSION_2026_07_25_HANDOFF_9.md` §0b**: `τ_own = 0` is the *interior*
  Schur result and the constraint `f_g ≤ 1 − |d|` is real algebra, but it is an *integrable spike*, not a
  hard wall, and turning it into an estimate needs the TILT (τ ≈ ±1) — which is population knowledge unless
  it comes from the neighbours' per-strand imputation. Prior-free they are
  carried by messages alone and the DL term does not protect them. The minimal reproduction is the factor-1 bedrock toy
  (`test_gdna_sweep_factor1_ambig_recovery`, xfail): on a uniform ρ=0.5 chain the AMBIG node between two exact
  anchors reads **0.3914**. This is the designed weakness, NOT a mode defect — the shortfall shrinks
  monotonically with depth (21.7% at ρ=0.5 → 0.8% at ρ=5000), so the transported mode is right and what is
  missing is WEIGHT: messages arrive at their honest `1/n`, and ψ's uninformative reference holds the node off
  the vertex until the data earn it. Fixed by §4 (a trained prior), not by more damping.

## 6. The path to production (ordered)

0. ✅ **DONE (2026-07-24):** converge to one solver — deleted `_scan` + flags + the `_UNIFIED_*` gates;
   extracted + hardened + unit-tested the **initialization** phase (`node_init.py`, the four sources);
   regenerated goldens. `bp_solver.py` 1871 → ~730 lines. Behavior-preserving (byte-identical A/B).
1. ✅ **DONE — the variance FOUNDATION** (`variance_foundation_plan.md` v4, `variance_foundation_proposal.md`).
   The honest local composition precision is a **single Schur-marginal scalar `τ_λ`** (approach E) — a diagonal
   `(τ_λ, τ_θ)` is prohibited (the strand Fisher is rank-1). Derived by a 5-approach workflow, numerically
   validated + independently re-verified, and independently critiqued (both converge on approach E + Option B).
   Phase 1 (the strand-gate bug fix — AMBIG gets zero strand f_g credit) is landed + A/B-validated (stranded arm
   improves, no regression), committed `c6df8c50`.
2. ✅ **DONE — the MESSAGE variance model** (§3; `message_variance_derivation.md`, laws M1–M7). The M1–M5
   per-component transport laws, the single-λ combine on a three-stream relay (the M6 rank-1 double-count), and
   the cross-cliff cost split into the M5 scale term + the M7 DerSimonian–Laird composition-mismatch `b̂²`. The
   composition/sampling separation came out of the three-stream relay (τ carries composition only, the
   measurement stream the counts). Every law MC-validated; `b̂²` additionally re-derived and adversarially
   audited by a 4-agent workflow. Commits `44f1ecc6`…`1a3e0a89`.
3. ✅ **DONE — the unified solver wins the A/B**: 0.0969 (refit=0) / 0.0828 (refit=1) vs the 0.0949
   legacy-with-factory target — and note the current suite gained the hard `verystrong`/`gdna1`/`gdna5`
   scenarios since that number was set, so gate on in-run A/B deltas, not the absolute.
4. ✅ **DONE — the composition peel, the MESSAGE PACKET, and P1d.** Suite 0.0865 (r0) / 0.0676 (r1);
   confidently-wrong 1,777,658 → **762,000**; `z2` ALL 15.55 → **8.98**, exon-single 8.60 → **2.46**.

5. ⭐ **P2a — RE-TEST THE PHASE-2 λ-CURVATURE STEP AT HEAD. THIS IS THE NEXT ACTION.** ~6 lines, already
   written and measured once at a much older HEAD (`SESSION_2026_07_25_HANDOFF_6.md` §3). Both reasons it was
   shelved have since evaporated (§4): the refit no longer regresses unstranded-capON, and the hyperprior's
   training substrate is now at `z2` = 2.26. **If it holds up, Phase 2 unblocks immediately** — which is the
   stated goal. Cheapest high-value experiment on the list by a wide margin. Gate it on the trust view, not
   on mwae.

6. **P1e — the conservation SURPRISE `z² = δ²/(αᵀΣα)`.** The only damping term that is **never inert** (M7's
   DL needs `τ_own > 0` and so switches off on all unstranded data; a node's own mass is always known). It is
   also the fix for **P1d's own residual regression**: that regression is localized to unstranded × gDNA-rich
   × capture-OFF, which is P1b's population, where the RNA claim is a near-exact measurement and the **gDNA**
   claim is 47× too big — so P1d currently damps the wrong arm. `claim/obs` is the prior-free observable that
   says which arm failed, and nothing uses it as evidence.

7. ✅ **DONE — P4: the `_far` across-the-seam level estimator is DELETED.** It was a **100 % structural
   data reuse**: the node it read is *identically* the source of the destination's OTHER message (verified
   on 35,421 peel edges, every condition), so that node's evidence reached the destination twice and the
   fuse added the two precisions as if independent. Deleting it is a **trust win as well as a correctness
   fix** — on a held-fixed node set (HEAD's own confident quartile) boundary-single `z2` **5.58 → 2.98**,
   boundary-AMBIG 18.89 → 15.84, exon-AMBIG 64.39 → 52.39, ALL 8.98 → **8.28**, because the declared
   variance there was inflated **2.81×** (8.54× in the top decile). Price **+0.0007 (r0) / +0.0009 (r1)**,
   0 better / 2 worse / 30 flat. ⚠ The plan's own proposed replacement (give the left message the RIGHT
   message's delivered claim) was prototyped and is **illegal and worthless** — it does not remove the
   reuse and recovers 0.000117 of 0.000696.

   ✅ **P4b — FIXED (`_RHO_ITERS = 1`). The solver is now BP-legal.** Price: **−0.0002 at refit=0** (4
   better / 2 worse / 26 flat — removing the violation is slightly BETTER) and **+0.0011 at refit=1**,
   concentrated in unstranded × capON (+0.0104). Held-fixed z2 flat (ALL 8.53 → 8.50). With one iteration
   the frame comes from `_init_belief()`, which is belief-free, so no posterior feeds back.
   **⛔ The "elegant" pairwise-potential refactor is NOT justified by this price** — it would bid a rewrite
   of the message representation against ~0.001. A belief-free constant frame was also prototyped and is a
   worse trade (+0.0003 suite, +0.0010 fit-substrate, neutral-to-worse on held-fixed z2).

   **(historical) `_far` was NOT the last third-node dependence.** `_RHO_ITERS = 2`, and the second iteration's
   reframe faces are built from the destination's **fused posterior** — i.e. from the other message.
   Measured median `|Δlog ρ_face|` between iterations 0.0116 (stranded capOFF) → 0.1242 (unstranded capON),
   >1 % on 52.7–79.0 % of nodes. **The solver is not BP-legal yet.** The legal way to recover what `_far`
   supplied, if wanted: carry the peel share as a FUNCTION of the destination's state (a proper pairwise
   potential) rather than a plugged-in point estimate.

8. ⭐ **P1g — PUT TSS/TES IN THE REGION/BOUNDARY MAP.** The structural debt behind `ω_graft` (see the banner
   at the top of this file). Owner-flagged as a high priority and forthcoming. It reaches the index build
   (`boundary_df` today carries only `boundary_id`/`ref_name`/`position`; `t_df` already has the transcript
   spans, so the bit is derivable with no new input). On landing, **re-derive `ω_graft` per structural
   class** — same equation, one scalar per class.

9. **P3 — AMBIG exon over-confidence** (`z2 | Q1` = 64.39, 259,493 confidently-wrong reads, 34 % of what
   remains). **Demoted from "promoted": it is excluded from the hyperprior fit by construction**, so it
   blocks the *re-solve's* quality, not Phase 2's start. Established as the designed weakness, not a mode
   bug (`SESSION_2026_07_25_HANDOFF_9.md` §0b) — a trained prior is the fix, which means it is largely
   *downstream* of Phase 2 rather than upstream of it.

10. **Then the re-solve and a ship candidate.**

11. **⚠ REFUTED, do not build: per-channel enrichment ratios at the boundary face.** The physics is real and
   measured (the boundary sits on a capture SLOPE: at verystrong it is 0.125× the exon and 2113× the intron),
   but item 4's relay pin retired the motivation — substituting the ORACLE capture step for the model's `r`
   buys ~0 (−0.16/−0.17/−0.33, +0.013) and ≈4 % at verystrong, because the reframe only carries SCALE and the
   pin re-derives scale from the observed mass. What remains live is frame ASSIGNMENT of the measured spliced
   channel (`density_composition_reconciliation.md` §5.1), not ratio accuracy.

## 7. How we work (methodology — see memory `pass0_debug_iteration_loop`)

Debug by the loop: **run the full ambig_dense_10mb suite → find the worst scenario (by error MASS) → dissect
its worst nodes → trace to root cause → fix → repeat.** Compare **pass-0 vs oracle only** (`calib_refit_iters=0`)
unless testing the refit. Everything is cached (`_selfsolve_cache`), so a scenario solve is ~1 s. Standing
directives: **no magic numbers** (pause and discuss before any new constant); develop on toy/controlled
scenarios; keep the module/constant count small. **⚠ The synthetic suite is Poisson by construction** (memory
`synthetic_suite_is_poisson_omega_zero`), so it **cannot validate anything overdispersion-dependent.**

Key tools (in `scripts/debug/`): `pass0_oracle_bench.py` (THE A/B — `--arm`, `P0_REFIT`, `--report`),
`pass0_error_concentration.py` (suite → where the error is), `pass0_node_dissect.py` (exact-replay ψ channel
ablation + per-node dumps), `message_variance_mc.py` (the variance-law MC arbiter, M1–M7),
`unified_message_audit.py` (Σ_c f_c invariant). In `scratchpad/`: `dl_dissect.py` (error mass attributed by
DL-protection state / strand DOF / node class, + per-node message provenance), `dump_node.py`.
Env `RIGEL_S2T_OFF=1` disables both cliff terms for isolation; `_capture["_dl"]` publishes the per-message
gaps and the τ-stream kill.
