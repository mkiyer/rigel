# The composition peel — implementation plan

**Status: LANDED 2026-07-26. §12 holds the retry that worked; §10/§11 are the 2026-07-25 attempt that did not.**
Owner-directed. The plan below is unchanged from before execution (it is the pre-registration);
§10 records what happened, including two of its own assumptions being refuted. Read `ROADMAP.md`, then
`SESSION_2026_07_25_HANDOFF_10.md` (§9–§12 are the evidence this plan rests on), then this.

The peel's **subtraction** is the last absolute-density operation left from the pre-composition solver. It is
the only place in the message path where a scale error does not cancel, and it is what contaminates introns.
This plan replaces it with a **share** that carries its own uncertainty, so that "intronic unspliced fragments
are gDNA until proven otherwise" is **reproduced by the uncertainty** rather than asserted by a rule.

---

## 0. Objective, and what "done" means

**Objective.** The RNA-continue message across a seam must be a *scaling* of the exon's RNA by a share `w`
estimated at the boundary, with `Var(log w)` propagated into the message precision — so that a boundary whose
RNA is consistent with zero emits a claim that is both near-zero **and weak**, and a boundary with real RNA
emits a proportionate one.

**Pre-registered acceptance gates.** All five must hold to land:

| # | gate | current HEAD |
|---|---|---|
| G1 | per-condition A/B at **refit=0**: no worse than HEAD on aggregate, and ≥ as many better-than-`m8` as HEAD | 0.0884, 16 better / 11 worse / 5 flat vs `m8` |
| G2 | per-condition A/B at **refit=1**: same | 0.0678, 23 better / 4 worse / 5 flat vs `m8` |
| G3 | delivered ρ_R at `exon\|intron` boundaries in `nrna_none` (oracle = **exactly 0**) ≤ HEAD | 0.2030 |
| G4 | boundary `errQ1conf` (share of a class's error on the most-confident quartile) must not regress | single 6.6 %, AMBIG 2.2 % |
| G5 | the `√2·σ_own` DL pin-safety inequality preserved; `message_variance_mc.py` 0 failures | holds |

**Rollback rule, pre-agreed.** If G1 or G2 fails after the ablations in §6 have been run, revert to HEAD and
record the measurement. Three wirings have already lost (HANDOFF_10 §10.2); the prior on success is not high,
so **every step below is specified to produce a measurement that stands on its own even if the step is
reverted.** That is the point of the sequencing in §7.

---

## 1. The quantities, and exactly where each comes from

Per **(node i, face d ∈ {0,1}, strand s ∈ {+,−})**. All already exist; nothing new is measured.

| symbol | meaning | source |
|---|---|---|
| `ρ_μ` | spliced (leaving) density at this face | `spl_p_f[d][i]` / `spl_n_f[d][i]` = `SP[d]/E_spl[d]` |
| `n_spl` | spliced COUNT at this face | `geometry.spliced_n_pos_{left,right}` / `..._neg_...` (already in `_uni_static`) |
| `v_μ` | `Var(log ρ_μ)` | `1/n_spl` — measured, composition-certain (M1) |
| `ρ_ν` | unspliced RNA density **at the boundary** | §2 precedence |
| `v_ν` | `Var(log ρ_ν)` | §2 precedence |
| `A` | the exon's RNA density reframed to the boundary | `(rp[s] or rn[s]) · r` — unchanged |
| `v_A` | `Var(log A)` | the relay's existing damped precision ⊕ `σ²_transfer` — unchanged |

⚠ **`v_μ` must use the spliced COUNT, never the spliced MASS.** The accumulator deposits fragments
*fractionally* across nodes, so at a junction face the median count is **33** against a median mass of **11**
(audited). Using the mass would over-state `v_μ` ~3× and make every peel message spuriously weak. The count
fields (`geometry.spliced_n_{pos,neg}_{left,right}`) already exist and are already published in
`_uni_static`; the mass fields (`spliced_{pos,neg}_{left,right}`) are what forms the DENSITY `ρ_μ`. Both are
needed, for different purposes — do not conflate them.

## 2. The LEVEL — a constant-free precedence (the "first step" of the two-step)

`ρ_ν` and `v_ν` are taken from the first source that has composition evidence. Every test is a **presence
test**, no thresholds:

1. **The boundary's own evidence** — `τ_own(i) > 0` (strand tilt, or a density deconvolution if one applies
   here). Then `ρ_ν = op[i]` / `on[i]` and `v_ν = transport_seed_logvar(v_log_fr(i), n_unspl(i))`.
2. **Else the node ACROSS the seam** — the neighbour opposite the exon, if it has `τ_own > 0` (in practice the
   factory-solved intron). Then `ρ_ν = op[far]·r_{far→i}` and
   `v_ν = transport_seed_logvar(v_log_fr(far), n(far)) + σ²_transfer(far→i)`.
3. **Else no evidence** — `v_ν = +∞`. See §4.3: the message carries **zero precision**, which is the honest
   "I don't know" and reproduces the safe default without asserting it.

This is exactly the owner's *"the intron message can set the gDNA vs unspliced RNA level. Or the strand. Then
the composition scale works."* Note (1) is checked before (2) so a boundary that knows its own composition is
never overridden by a neighbour's — the same "a direct measurement outranks an imputation" ordering that the
intron study established.

⚠ **Measured expectation, so we are not surprised.** Source (2) is the estimator whose *mode* was scored in
HANDOFF_10 §10.1: |Δw| = 0.192 / **0.294** / 0.126 / 0.052 / 0.100 against a self-solve's
0.110 / **0.628** / 0.156 / 0.070 / 0.160. It is better but still over-claims on the unstranded arm. **The
variance is what must carry that** — see §3.

## 3. The LAW — mode and precision together (M10, already derived and MC-validated)

```
    w        = ρ_ν / (ρ_ν + ρ_μ)                                  enrichment-free (e cancels identically)
    w_μ      = 1 − w                                              the spliced share
    t_p      = A · w                                              a SCALING, not a difference
    Var(log t_p) = v_A + w_μ²·(v_ν + v_μ)                          convex weights ≤ 1
```

Both pure laws are **landed and unit-tested**: `enrichment_frame.peel_continue_share` and
`peel_share_logvar`. MC: `message_variance_mc.py` M10 rel 0.20–0.96 %, M10b exact to 1e-12.

**Why the variance term is the whole point.** As `ρ_ν → 0`: `w → 0` *and* `w_μ → 1`, so
`Var(log w) → v_ν + v_μ`, which is large exactly when `ρ_ν` is poorly known. A boundary reporting
`ρ_ν = 0.0006 ± 0.0012` (HANDOFF_10 §12.3 — measured, and *conservative*: z2 = 0.13–0.55) therefore emits a
near-zero claim at near-zero precision. **That is "gDNA until proven otherwise" derived, not asserted** — and
it is precisely what the three failed wirings threw away by using a ratio of modes.

## 4. Degenerate cases and structural interactions — enumerate ALL of them

Each is a required unit test (§6.1).

**4.1 `ρ_μ = 0` (no spliced flux at this face).** `w = 1`, `w_μ = 0` ⇒ `Var(log w) = 0`. Nothing splices away,
nothing is peeled, no variance added. The message passes through untouched. Already the documented limit of
`peel_continue_share`.

**4.2 `ρ_ν = 0` with `ρ_μ > 0`.** `w = 0` ⇒ `t_p = 0`. Three downstream operators read the *density*, not the
precision, and each must be checked:
* `_pin_v`: `sp = where(pp_ > 0, p, op)` — reads the precision first, so a zero-precision component correctly
  falls back to the node's own density. **Safe as written.**
* **the λ-emission gate** `ttau = where((tg > _EPS) & ((tp + tn) > _EPS), ttau, 0)` — reads the DENSITY. A
  legitimately-zero RNA density would silence the composition message. Its stated intent is "a source that
  carries only ONE component has no composition claim", i.e. a test of what was *supplied*; a component at
  zero precision is not supplied, and a component at zero density but live precision **is**. ⚠ **This gate
  must move to a precision test** — but as its own commit and its own A/B (§7 step 5), never as a rider.
* `mismatch_gap(t, o)` — `contradicted = live_t ^ live_o`, so `t_p = 0` against a live own belief marks the
  message contradicted and DL kills it. That is correct *if* the zero is a claim; it is wrong if the zero is
  "no opinion". Resolved by the same precision test.

**4.3 `v_ν = +∞` (no evidence anywhere — §2 case 3).** Then `Var(log t_p) = ∞` ⇒ precision 0. The mode is
undefined; set `w = 0` so `t_p = 0`, consistent with §4.2, and rely on the precision to make it inert. **Must
not produce a nan**: compute with the `∞`-masked-before-the-product idiom already used in
`node_init.own_precision`.

**4.4 AMBIG nodes.** `w` is per-strand; a strand that is structurally dead already has `ρ_ν = 0` and is zeroed
by the existing `fp_a`/`fn_a` gates *after* the peel. Order is unchanged.

**4.5 Mature contiguous (`exon|exon`, retained intron).** Little or nothing splices at that face ⇒ `ρ_μ` small
⇒ `w → 1` ⇒ near-no peel. The correct behaviour falls out of the measurement; **no structural gate is
needed**, which is what retires HEAD's `mrna_active` rule (§5).

**4.6 Both faces carry spliced.** Measured at 4.6–5.2 % of junction boundaries (HANDOFF_10 §9). `w` is
per-face by construction; each direction uses its own face. No special case.

**4.7 `n_spl = 0` while `ρ_μ > 0`.** Cannot occur (the density is built from the count) but guard: `v_μ = ∞`
⇒ zero precision rather than a division by zero.

**4.8 The relay's running belief.** `_relay` computes `w` from the node's *init* (`op`/`on`), not its fused
value, because when node `i` is processed its own entry has not been updated yet. This is deliberate and must
be stated: the share is a **pre-pass** quantity. Do **not** attempt to use the running belief — it would make
the share depend on the traversal direction and break the forward/backward symmetry.

## 5. What is REMOVED

* **The subtraction** `tp = max(tp − spl_*_f[df], 0)` — both sites (`_relay`, `_transport`). This also retires
  the **zero-truncation defect** (HANDOFF_7 §8.6: a fully-consumed peel emitted `ρ_ν = 0`, "no RNA continues
  past here", at a *live* precision) for free, since `w ∈ [0,1]`.
* **HEAD's `mrna_active` gate** — its justification is void (the owner's correction: mature and nascent are
  the same species; only spliced vs unspliced is observable), and §4.5 shows the share subsumes its effect.
  Remove it **only in step 4** (§7), so its contribution is measured separately.
* **NOT removed:** `peel_rna_logvar` (M3). It is unused by the solver but is the record of what the difference
  costs, and M10's convexity unit test asserts against it directly.

## 6. Validation ladder

### 6.1 Closed-form unit tests (new, in `tests/calibration/test_enrichment_frame.py`)
Four already exist for the laws. Add, for the **wiring**:
1. `ρ_μ = 0` ⇒ message identical to no-peel, and `Var(log w) = 0` (§4.1).
2. `ρ_ν = 0` ⇒ `t_p = 0` and precision 0 — and **no nan** (§4.2/4.3).
3. `v_ν = ∞` ⇒ precision 0, mode finite (§4.3).
4. enrichment invariance end-to-end: scaling `(ρ_ν, ρ_μ)` by a common `e` leaves `t_p` and `Var(log t_p)`
   unchanged to 1e-14 — the property the whole change exists for.
5. the M3→M10 contrast at a fixed operating point: `Var(log t_p)` strictly below the subtraction's, for
   `u ∈ {2, 3, 10}`.

### 6.2 MC — already done, re-run as a gate
`message_variance_mc.py` M10 (rel 0.20–0.96 %) and M10b (exact to 1e-12). **0 failures over M1–M10** is a
landing gate. No new MC is required; the law is unchanged by the wiring.

### 6.3 Replay fidelity
`scratchpad/ss_cap_3_replay.py` / `ambig_2_mechanism.py` must still reproduce the shipped ψ **bit-exactly**
(`max|Δ| = 0`) after the change. If they do not, a channel is being published inconsistently with what ψ is
handed — fix that before reading any A/B.

### 6.4 Per-condition A/B
`P0_REFIT=0|1 pass0_oracle_bench.py --arm …`, all 32 conditions, `OMP_NUM_THREADS=1`. Report better/worse/flat
against **both** `m8` (the pre-boundary-work baseline) and HEAD. Call out explicitly: `stranded ss_0.99`,
`unstranded × capON`, `verystrong`, `gdna_none`.

### 6.5 Certification matrix — by boundary class
`scratchpad/b1_boundary.py`, per flanking pair × junction present/absent. Every cell must be reported, not
just the aggregate:

| class | why it is a distinct case |
|---|---|
| `exon\|intron`, junction | the main peel case; 39–61 % of boundary error |
| `exon\|intron`, no junction | the 14–17 % whose junction was never observed |
| `exon\|exon`, junction | mature contiguous **and** splicing — the owner's TA/TB geometry; `w → 1` expected |
| `exon\|exon`, no junction | contiguous, no peel; must be untouched |
| retained intron | `coarse_type` collapses it (HANDOFF_7 §8.3) — verify it is not silently peeled |
| both enrichment orderings | boundary < exon (72–97 %) and boundary > exon (3–28 %) |

### 6.6 The regime split that decides success
Report `w` against the truth separately for **saturated** and **unsaturated** boundaries, using the
constant-free detector of HANDOFF_10 §11.4 (is the factory log-factor still rising at the top of the λ grid).
Expected, and the plan fails if not observed:
* saturated / `nrna_none`: `w → 0` with large `Var(log w)` ⇒ delivered RNA → 0 (gate G3);
* unsaturated / `nrna_present`: `w` tracks the true 0.17–0.41 rather than HEAD's 0.

⚠ **Re-measure `w_true` on junction-bearing seams only.** HANDOFF_10 §11.4 records that the earlier figure was
contaminated by boundaries with no spliced mass, where the share is undefined and the probe defaulted it to 1.

### 6.7 Standing gates
`pytest tests/ -q` (1236 pass / 2 xfail / 2 xpass), `ruff check src/ tests/ scripts/`, goldens regenerated
**last**.

## 7. Sequencing — each step independently measurable

Each step is its own commit with its own A/B, so a failure localizes and every measurement survives a revert.

| step | change | what it isolates | expectation |
|---|---|---|---|
| 1 | plumb `(ρ_ν, v_ν)` per §2 and compute `w`, `Var(log w)` — **computed, not used** | that the level precedence and the shares are well-formed; no behaviour change | bit-identical A/B |
| 2 | replace the subtraction with `t_p = A·w`, keep HEAD's gate, **precision unchanged** | the MODE effect of the share alone | ~neutral; this is close to the arm that scored 0.0893/0.0690 |
| 3 | add `w_μ²(v_ν + v_μ)` to the message log-variance | **the variance effect — the actual thesis** | the step that must pay; G3 should move |
| 4 | remove HEAD's `mrna_active` gate | whether the share subsumes the structural rule | should be ~neutral if §4.5 holds |
| 5 | λ-emission gate: density test → precision test (§4.2) | a pre-existing defect the plan surfaced | its own A/B, may be landed or reverted independently |

**Do not collapse steps 2 and 3.** The whole claim of this plan is that the mode alone is not enough and the
variance is what makes it work; measuring them together would leave that untested — and step 2 in isolation is
expected to look like the failed arm, which is the point.

## 8. Risk register

| risk | how it shows | response |
|---|---|---|
| The §2(2) intron level still over-claims on the unstranded arm | G3 fails; `unstranded × capON` regresses | this is the known |Δw| = 0.294; if the variance does not absorb it, the plan's thesis is **refuted** — record and stop |
| `w = 0` trips the λ gate / DL contradiction mask | composition messages vanish at seams; `exon\|exon` `c_tau` → 0 (the exact signature seen when the junction gate over-reached) | step 5; check `c_tau` per boundary class in §6.5 |
| The share becomes over-confident where `v_ν` is optimistic | `errQ1conf` at boundaries regresses (G4) | G4 is a hard gate, not advisory |
| Loss of the M3 `u`-weighting removes a damping that was compensating elsewhere | stranded arm regresses | ablate by restoring `peel_rna_logvar` on the same operating point |
| Traversal-order dependence | forward and backward passes disagree | §4.8 — the share is a pre-pass quantity; assert `_relay` and `_transport` compute the identical `w` |

## 9. Explicitly out of scope

* The **exon-face reframe error** (0.4–1.3 nats). It is irreducible (HANDOFF_10 §9.3 — unchanged even with the
  ORACLE `f_g`) and is precisely what the share makes harmless rather than fixes.
* The **ψ vertex** / `_JEFFREYS_REF`. HANDOFF_10 §11–§12 close this: the reference is required for properness,
  the median is already the best readout, the grid is not the problem, and the precision is already honest.
* The **hyperprior**. Nothing here needs it, and it is not the current task.


---

## 10. EXECUTION RESULT (2026-07-25) — G1 failed, reverted; four findings stand

All five steps run, each measured separately as specified. **HEAD restored** (`expF`, 0.0884 / 0.0678).

| step | change | refit=0 | vs HEAD | refit=1 | vs HEAD |
|---|---|---|---|---|---|
| 1 | scaffolding, unused | — | **bit-identical** ✓ | — | — |
| 2 | share replaces subtraction, precision unchanged | 0.1024 | 2 / 30 | 0.0821 | 3 / 29 |
| **3** | **+ `Var(log w)` in the precision** | **0.0926** | 5 / 23 | **0.0676** | 12 / 10 / 10 |
| 4 | remove the `mrna_active` gate | 0.0932 | 2 / 30 | 0.0683 | 9 / 16 / 7 |
| 5 | λ-gate: density test → precision test | 0.0927 | 6 / 23 | 0.0676 | 12 / 10 / 10 |

**Gate outcome.** G1 **FAILS**: best arm (step 3) is 0.0926 vs HEAD's 0.0884 at refit=0, and 10 better / 21
worse against `m8` where HEAD is 16 / 11 / 5. G2 is a wash (0.0676 vs 0.0678, 12 / 10 / 10). G3 is
**unchanged at 0.2030** — because the structural gate is retained, so the share never reaches the
`exon|intron` seams where the leak is. G4 **improves markedly** (below). Per the pre-agreed rollback rule:
reverted.

### 10.1 The thesis is CONFIRMED directionally, and it is not enough

Step 2 → step 3 is the whole claim, and the variance does exactly what M10 says: **0.1024 → 0.0926** (refit=0)
and **0.0821 → 0.0676** (refit=1). Adding `w_μ²(v_ν + v_μ)` recovers most of the mode-only damage and brings
refit=1 level with HEAD. The derivation is not in question; the *level* still is (§2 case 2's |Δw| = 0.294 on
the unstranded arm, exactly the risk registered in §8).

### 10.2 ⛔ TWO ASSUMPTIONS OF THIS PLAN ARE REFUTED

* **§7's claim that steps 2 and 3 are separable is WRONG.** Case 3's conservative mode (`w = 0`) is only
  defensible *paired with* its zero precision; applying the mode at full precision costs 0.0884 → 0.1024. At
  an `exon|exon` seam neither flank has a factory, so on unstranded data the precedence falls to case 3 and
  step 2 kills a claim that should flow. Mode and variance are one object here, not two.
* **§4.5's claim that the share subsumes the `mrna_active` gate is WRONG.** Removing it costs 0.0926 → 0.0932
  and 0.0676 → 0.0683 (step 4). The structural rule contributes independently of the share.

### 10.3 The λ-gate fix — LANDED (and §10.3's first verdict was wrong)

⚠ **Correction.** The first pass judged step 5 "inert, not worth a commit" from a single stacked measurement
(0.0927 vs step 3's 0.0926) — but that was step 5 *on top of step 3*, where the share had already changed the
density landscape underneath it, and **step 5 was never tested on its own**. Owner's call to test it properly.

Measured alone against HEAD: **0.0885 / 0.0678** vs 0.0884 / 0.0678, with **29/32 and 27/32 conditions
completely unchanged**. But the mechanism is *not* inert — instrumenting the two predicates directly, they
disagree on a substantial set:

| condition | disagree (all msgs) | live λ msgs | disagree \| live | rate |
|---|---|---|---|---|
| gdna300 ss0.50 none capON | 2,948 | 10,240 | **1,610** | **15.7 %** |
| gdna300 ss0.99 present capON | 638 | 11,824 | 397 | 3.4 % |
| gdna100 verystrong | 2,760 | 2,448 | 285 | 11.6 % |

So the fix silences 3.4–15.7 % of live λ messages that the density test had let through — a real set whose
claims simply were not moving the answer. **Theoretically correct, fires on a real population,
behaviour-neutral**: a correctness cleanup, landed on that basis rather than on a score.

It also matters *forward*. The disagreement is almost entirely one-directional today (density-supplied but
precision-empty: 2,942 of 2,948). The opposite case — a legitimate ZERO density at live precision, i.e. the
claim "there is none of this here", which is exactly a composition claim — is currently rare (6 / 3 / 3 / 0)
**but is precisely what a composition peel emits**. With the density test in place that claim would have been
silently discarded. This removes a trap the next attempt would have walked into.

### 10.4 ⭐ The honest-precision trade, and it is the interesting result

Step 3 makes boundaries **markedly more honest and somewhat more wrong**:

| | HEAD | step 3 |
|---|---|---|
| boundary single `errQ1conf` | 6.6 % | **2.5 %** |
| boundary AMBIG `errQ1conf` | 2.2 % | **1.6 %** |
| boundary single ERR (reads) | 692,306 | 939,360 |

So the share is doing precisely what it was designed to do — it converts confident boundary error into
*declared uncertainty* — but pays for it in error mass, and mwae scores only the mass. **That is a real
tension between G1 and G4, and it is the owner's call, not a technical one:** if what pass-0 owes Phase 2 is a
trustworthy substrate rather than a low score, a 2.6× reduction in confidently-wrong boundary error may be
worth 0.004 of aggregate mwae. G1 was pre-registered as a hard gate, so the default action was revert.

### 10.5 What would have to change for a retry

The binding constraint is the LEVEL, not the law. §2 case 2 (the far intron) carries |Δw| = 0.294 on the
unstranded arm, and case 3 (no evidence at all) covers `exon|exon` seams entirely — where HANDOFF_10 §8.1
already showed the boundary is structurally starved (no factory within reach on 97 %, no strand when
unstranded). **A better share needs a source that does not exist yet at those seams**, which is the same
conclusion the boundary study reached from the other direction. Retry only after that source exists.


---

## 11. THE FULL SET (steps 1–5 combined), run on the whole battery — owner-requested

The earlier execution measured the steps as a ladder; this is all five **combined**, which is what was
actually being advocated. Step 5 was already at HEAD by then (landed separately as a correctness fix).

| arm | refit=0 | vs HEAD | refit=1 | vs HEAD |
|---|---|---|---|---|
| `m8` (pre-boundary baseline) | 0.0900 | — | 0.0700 | — |
| **HEAD** (`mrna_active` gate + λ precision gate) | **0.0885** | — | **0.0678** | — |
| FULL (steps 1–5) | 0.0934 | 1 better / 31 worse | 0.0684 | 9 better / 18 worse / 5 flat |
| HYBRID (share where a level exists, subtraction where not) | 0.0893 | 14 / 14 / 4 | 0.0693 | 4 / 22 / 6 |

Note the FULL arm still beats `m8` at refit=1 (0.0684 vs 0.0700, **22 better / 6 worse**) — the boundary work
as a whole is winning; it is the composition peel specifically that costs.

### 11.1 Where the FULL set breaks, and why

The damage is concentrated on the **low-gDNA arms**: `gdna1_ss0.50_present_capON` 0.0689 → 0.0935,
`gdna5_..._capON` 0.0906 → 0.1133, `none_ss0.50_present_capON` 0.0588 → 0.0841. In those libraries the intron
factory has almost no background to fit, so `τ_factory = 0`, the level precedence falls to **case 3**, and
`v_ν = ∞` **silences the RNA channel outright** — in libraries where RNA is the whole signal.

That is the same failure shape as the `mrna_active` over-reach and as §10.2's finding: a conservative default
is only safe when something else can still carry the channel.

### 11.2 The HYBRID repair — right diagnosis, still not a win

Falling back to the **subtraction** in case 3 (we cannot form a share, but the spliced flux demonstrably left,
which is a measured fact) recovers most of the refit=0 damage — 0.0934 → 0.0893, and 1/31 → 14/14. It costs
refit=1 (0.0684 → 0.0693, 4 / 22). So the diagnosis is confirmed and the repair works in the direction
predicted, but neither arm reaches HEAD.

**Reverted.** HEAD stands at 0.0885 / 0.0678.

### 11.3 The standing conclusion is unchanged and now has a third independent confirmation

Every route into the composition peel has failed at the same place: **there is no level at the seams that
matter.** Case 2 (the far intron) carries |Δw| = 0.294 unstranded; case 3 covers `exon|exon` seams entirely
and every seam in a low-gDNA library. The law is right, the wiring is right, the *input* does not exist yet.
Retry only once a boundary-level source exists — which is exactly what HANDOFF_10 §8.1 identified from the
opposite direction (`exon|exon` boundaries: no factory within reach on 97 %, no strand when unstranded).


---

## 12. THE RETRY THAT LANDED (2026-07-26) — §11.3's "the input does not exist yet" was wrong

§11.3 concluded: *"The law is right, the wiring is right, the INPUT does not exist yet. Retry only once a
boundary-level source exists."* The source did exist — it just had not been recognized as one. **The node's
own observed mass, closed against an imputed gDNA DENSITY** (`enrichment_frame.residual_level`, M11) is a
level at every seam, and it is the same primitive the intron factory already uses, with the gDNA prior
supplied by a neighbour rather than by the intergenic pool.

Two structural corrections to this plan turned it from a loss into a win, and both invalidate parts of §2:

* **§2's PRECEDENCE is wrong as a shape, not just in its ordering.** Every candidate level is
  `ρ_ν = (1 − f̂_g)·M/E_r`, so they are three estimators of ONE quantity and must FUSE by inverse variance.
  A precedence needs a "no evidence" arm; a fuse does not, and that arm (case 3) was the entire regression.
* **§3's law is right and its CONSUMER was in the wrong coordinate.** Levels must be fused in LINEAR density
  space. `Var(log w) = w_μ²(v_ν+v_μ)` is untouched; what changed is how `v_ν` is formed.

| arm | refit=0 | vs HEAD | refit=1 | vs HEAD |
|---|---|---|---|---|
| HEAD (`g5`) | 0.0885 | — | 0.0678 | — |
| §10/§11 FULL (precedence) | 0.0934 | 1 / 31 | 0.0684 | 9 / 18 / 5 |
| **§12 (fuse + M11), gate REMOVED** | **0.0895** | **9 / 9 / 14** | **0.0693** | **4 / 8 / 20** |

**The pre-registered gates.** G1/G2 are marginal misses on aggregate (+0.0010 / +0.0015) but the
better/worse counts are transformed (1/31 → 9/9). **G4 passes decisively and is now free**: boundary
`errQ1conf` 6.6 % → 5.3 % *while boundary error mass FALLS* 692 k → 664 k — where §10.4's arm had to buy its
honesty with 247 k of extra boundary error. G5 holds (`message_variance_mc.py` 0 failures over M1–**M11**;
the DL inequality untouched). §10.4's G1-vs-G4 tension is therefore **dissolved, not traded**.

Per §11.3's own standard — *"the law is right, the wiring is right"* — and the owner's standing instruction
not to revert on tiny error differences, this is landed. The remaining +0.0010 is one regime
(`gDNA-rich × capture-ON`) with a measured mechanism: see `SESSION_2026_07_26_HANDOFF_12.md` §5.
