# Session handoff — the BOUNDARY work: steps 1–5 are AGREED, the revert is PROVISIONAL

**This is the LIVE handoff. START HERE.** Date: 2026-07-25. Branch `calib-ambig-init-wip`, HEAD `9cde3067`.
Suite: **0.0885 (refit=0) / 0.0678 (refit=1)**. Gates all green.

---

## 0. Read this first — the owner's position

**Steps 1–5 of `peel_by_composition_plan.md` are AGREED IN PRINCIPLE. The revert at the end of the last
session was PROVISIONAL — a clean-tree decision, not a rejection of the design.** The owner's verdict, per
step:

| step | what it is | owner's position |
|---|---|---|
| **1** | plumb the boundary's RNA *level* `(ρ_ν, v_ν)` and form the share `w`, `Var(log w)` | **KEEP.** Owner has ideas here. |
| **2** | peel by **scaling** (`×w`) instead of **subtracting** | **Confident we must scale.** The open question is *how*. Owner has experiments in mind. |
| **3** | carry `Var(log w)` into the message precision | **Absolutely required.** The open question is how to compute the uncertainty *correctly*. **Owner is working on this now.** |
| **4** | remove the `mrna_active` structural gate | **MUST remove.** It is a band-aid — a last resort that means we gave up making the algorithm correct. |
| **5** | the λ-emission gate tests **precision**, not density | **Agreed, critical.** ✅ *Already landed* (`4d3340ab`). |

> Owner: *"I agree in theory with steps 1–5. I think they give us a more solid foundation to work with and
> make our approaches elegant and consistent."*
>
> And on method: *"We shouldn't be making hard decisions about what to revert or keep based on tiny error
> differences. If it's a theoretically sound idea without a catastrophic regression, just test it across the
> entire AMBIG scenario battery and see how it performs."*

**⚠ Do not re-litigate whether to do steps 1–4. Do them. The task is to make them work.**

## 1. ▶ THE THREE QUESTIONS FOR THIS SESSION

1. **Why the regression?** Steps 1–5 combined score 0.0934 / 0.0684 against HEAD's 0.0885 / 0.0678. §5 has
   the localization already done — start there, do not re-measure it.
2. **How do you compute the precision of a density RATIO?** This is step 3's open question and the owner is
   working it in parallel. `w = ρ_ν/(ρ_ν+ρ_μ)` — the current answer is M10's delta method,
   `Var(log w) = w_μ²(v_ν + v_μ)`, which is MC-validated *given* `v_ν`. The suspicion is that `v_ν` itself is
   wrong, not the propagation.
3. **How should unspliced boundary fragments be split RNA vs DNA?** This is the level problem (§6) and it is
   what every failed attempt has bottomed out on.

## 2. Bootstrap — environment, gates, restore

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
export OMP_NUM_THREADS=1                      # REQUIRED for A/B determinism
cd /Users/mkiyer/proj/rigel
```

**Restore steps 1–4 in one command** (verified bit-identical to the `full_r0`/`full_r1` arms):

```bash
git apply scratchpad/steps_1_to_4.patch        # step 5 is already at HEAD
```

**Gates — all must stay green:**

```bash
python -m pytest tests/calibration tests/native -q   # fast loop
python -m pytest tests/ -q                           # 1236 pass, 2 xfail, 2 xpass
ruff check src/ tests/ scripts/                      # clean
python scripts/debug/message_variance_mc.py          # 0 failures over M1–M10
python -m pytest tests/ --update-golden -q           # goldens — LAST, after the solver is final
```

**The A/B** (32 conditions, cached, ~1 s/scenario, `~/Downloads/rigel_runs/ambig_dense_10mb`):

```bash
P0_REFIT=0 python scripts/debug/pass0_oracle_bench.py --arm NAME_r0
P0_REFIT=1 python scripts/debug/pass0_oracle_bench.py --arm NAME_r1
# rows accumulate in /tmp/pass0_oracle_bench.tsv; HEAD's arms are `g5_r0` / `g5_r1`
```

## 3. Where we are — the arms on record

| arm | refit=0 | refit=1 | what |
|---|---|---|---|
| `m8` | 0.0900 | 0.0700 | before any boundary work |
| **`g5` = HEAD** | **0.0885** | **0.0678** | `mrna_active` gate + λ precision gate |
| `full` | 0.0934 | 0.0684 | steps 1–5 combined |
| `hyb` | 0.0893 | 0.0693 | share where a level exists, subtraction where not |
| `s3` | 0.0926 | 0.0676 | steps 1–3 + 5, gate KEPT (best refit=1 of the peel arms) |

Note `full` still beats `m8` at refit=1 (0.0684 vs 0.0700, **22 better / 6 worse**) — the boundary work as a
whole is winning; the composition peel *specifically* is what costs. **All differences are small; treat them
as directional, not decisive** (owner's standing instruction).

Suite in reads: **12.34 M of 116.7 M node-attributed fragments misassigned (10.6 %)**.
`python scripts/debug/pass0_error_table.py` gives the full per-scenario table plus the TRUST view.

## 4. What each step IS (so the next session can reason about them, not just re-run them)

A boundary sits between two regions. RNA arriving from an exon either **splices away** or **continues
unspliced** across the seam. The solver must decide how much continues. The old code **subtracted**:
`continues = (exon RNA, rescaled) − (measured spliced density)`.

* **Step 1 — the level.** For each boundary estimate `ρ_ν` (how much unspliced RNA is there) *and* `v_ν` (its
  log-variance), then form `w = ρ_ν/(ρ_ν+ρ_μ)` and `Var(log w)`. Precedence, all presence tests: the
  boundary's own evidence (`τ_own > 0`) → the node across the seam (the factory-solved intron) → else
  `v_ν = ∞`. **Bit-identical when computed but unused** (verified).
* **Step 2 — scale, don't subtract.** `continues = (exon RNA) × w`. `w` is enrichment-free: capture
  multiplies the continuing and splicing channels alike, so `e` cancels inside the ratio. **Why it matters:**
  the boundary's enrichment is uncertain by 0.4–1.3 nats *irreducibly* (a boundary samples an `fl_mean` window
  around a point, an exon samples its whole interior — unchanged even with the ORACLE `f_g`). A subtraction
  **amplifies** that by `u = 1/w`: measured **1.77× / 2.39× / 5.01×** at δ = 0.30 for u = 2/3/10. A
  multiplication does not amplify it at all. Same information, better conditioning.
* **Step 3 — carry the uncertainty.** `Var(log w) = w_μ²(v_ν + v_μ)`, `w_μ = 1−w`. When a boundary says
  "ρ_ν = 0.0006 ± 0.0012, consistent with zero", the message becomes both near-zero **and weak**, so it cannot
  contaminate — *"gDNA until proven otherwise" arising from the uncertainty rather than a rule.*
* **Step 4 — remove the structural gate.** HEAD kills the exon's RNA claim where `mrna_active_s` is false.
  Owner: a band-aid, must go.
* **Step 5 — the λ-gate.** The λ message is a node's claim about its *composition ratio*; a source that
  supplied only one component has no such claim. "Supplied" is a statement about **precision**, not density:
  *density > 0, precision = 0* is "a number I have no confidence in" (not a claim); *density = 0, precision >
  0* is "I am confident there is **none** here" (very much a claim). Measured: the two tests disagree on
  **3.4–15.7 %** of live λ messages. ✅ Landed.

## 5. ⭐ WHY THE REGRESSION — the localization is already done

The damage from `full` is **concentrated on the low-gDNA arms**:

| condition | HEAD | full |
|---|---|---|
| `gdna1_ss0.50_present_capON` | 0.0689 | **0.0935** |
| `gdna5_ss0.50_present_capON` | 0.0906 | **0.1133** |
| `none_ss0.50_present_capON` | 0.0588 | **0.0841** |
| `gdna5_..._verystrong` | 0.1550 | 0.1666 |

**Mechanism.** In a low-gDNA library the intron factory has almost no background to fit ⇒ `τ_factory = 0` ⇒
the level precedence falls through to **case 3** ⇒ `v_ν = ∞` ⇒ the RNA channel is **silenced outright** — in
libraries where RNA is the entire signal. Same failure shape as the `mrna_active` over-reach: *a conservative
default is only safe when something else can still carry the channel.*

**Confirmed by repair.** Falling back to the subtraction in case 3 (the spliced flux demonstrably left — a
measured fact even without a level) recovers most of refit=0: **0.0934 → 0.0893**, and 1/31 → 14/14. It costs
refit=1 (4/22). Right diagnosis, not yet a win. That is the `hyb` arm.

**So question 1 has an answer already: the regression is case 3.** The next session's job is not to find it
but to fix it — i.e. to give case 3 something better than silence *and* better than the subtraction.

## 6. ⭐ THE STANDING BLOCKER — three independent confirmations

**There is no LEVEL at the seams that matter.** Every route into the peel has bottomed out here:

| route | how it failed |
|---|---|
| level from the boundary's own self-solve | |Δw| = **0.628** unstranded (ψ's Beta(½,½) leaves it at ~0.5) |
| level from the far intron | |Δw| = **0.294** unstranded — better, still over-claims |
| no level (case 3) | silences the channel; kills low-gDNA libraries |
| "gDNA until proven otherwise" as a hard default | breaks `_pin_v`'s partial-claim semantics ⇒ asserts `f_g → 1` |

And the census says why, structurally: **`exon|exon` boundaries have no factory within reach on 97 % and no
strand when unstranded** — they are the crux case with the fewest sources (HANDOFF_10 §8.1). Case 3 covers
`exon|exon` seams *entirely* and every seam in a low-gDNA library.

**This is question 3.** A boundary's unspliced crossing is `gDNA + RNA-that-continues`; splitting it needs a
source that does not exist yet at those seams. Candidate directions not yet tried:
* the exon's **gDNA** claim (which transports at 0.96–0.99 accuracy — the *good* channel) to set the gDNA
  level, leaving RNA as the residual against the observed mass;
* a genuinely two-sided boundary solve using the measured spliced/unspliced split (the owner's framing: the
  boundary *does* measure spliced vs unspliced — that composition is observed);
* the intergenic anchors, which are exact (`intergenic → boundary` reframe error **0.000**).

## 7. Question 2 — the precision of a density ratio

`w = ρ_ν/(ρ_ν+ρ_μ)`. Current answer (M10, MC-validated to 0.2–1.0 %, `enrichment_frame.peel_share_logvar`):

```
    log w = log ρ_ν − log(ρ_ν+ρ_μ) ,   ∂log w/∂log ρ_ν = 1−w = w_μ ,  ∂log w/∂log ρ_μ = −w_μ
    Var(log w) = w_μ²·( v_ν + v_μ )        convex (w_μ ≤ 1) — the mirror of M2's graft SUM
```

**The propagation is validated; `v_ν` is the suspect.** Things to check:
* the delta method is a *linearization* — near `w → 0` (exactly the RNA-free case) `log w` is unbounded and
  the linearization may understate. Compare against a direct posterior for `w` (it is a bounded [0,1]
  quantity, so its posterior mean/variance are well-defined without any log).
* `v_ν` currently = `transport_seed_logvar(v_log_fr, n)` = `Var(log f_R) + 1/n` with
  `Var(log f_R) = f_g²/τ_λ`. At `f_g → 1` that is large — correct — but it is the **λ-Jacobian** variance,
  and the measured node-level calibration is *conservative* (z2 = 0.13–0.55, HANDOFF_10 §12.3).
* `v_μ` **must use the spliced COUNT, never the mass** — the accumulator deposits fragments fractionally
  (median count 33 vs median mass 11 at a junction face); the mass would over-state `v_μ` ~3×.

## 8. ⛔ DO NOT RE-RUN — measured and settled

Everything in HANDOFF_7 §7 and HANDOFF_8 §7 still holds, plus:

| item | verdict |
|---|---|
| ψ vertex: finer grid / linear-scale grid / different point distribution | **All no.** Resolution near the vertex is Δf_g ≈ 0.0014 (15× finer than the bias); the median is transform-invariant so a coordinate change cannot move it; the measure is an explicit reference either way. HANDOFF_10 §12.1–12.2 |
| dropping `_JEFFREYS_REF` where the factory is present | Makes ψ **improper** — the median drifts 0.9835→0.9995 with L. §11.2 |
| the intron λ-gate ("direct measurement outranks imputation") | Superseded by the message fix; removed |
| the left/right message gap as a DL second study | Refuted: corr 0.028 / 0.032 / −0.116 with the mode error |
| `r_M/E_g` as a better frame ratio | Moot — `_pin_v` cancels `r` on 87.6 % of the error |
| more `_RHO_ITERS` (2→4→8) | No effect; the combine is converged |
| `gdna_fallback_admissible` wired as a hard default | 0.0906 / 0.0690, 3 better / 15 worse |
| per-channel enrichment ratios; the graft-frame fix; the RNA floor; `σ²_pin` | HANDOFF_7 §7 |

## 9. Documents, in reading order

1. **`ROADMAP.md`** — the entry point.
2. **THIS FILE** — state, the owner's position, the three questions.
3. **`peel_by_composition_plan.md`** — the plan for steps 1–5 with pre-registered gates; **§10 and §11 hold
   the execution results.** ⚠ Two of its own assumptions are refuted in §10.2 (steps 2 and 3 are *not*
   separable; the share does *not* subsume the gate).
4. **`SESSION_2026_07_25_HANDOFF_10.md`** — the boundary evidence base. §8 the `mrna_active` rule, §9 the
   two-sided reframe (validated/refuted/irreducible), §10 the share estimators, §11 the ψ vertex, §12 the
   grid + precision-honesty answer.
5. `message_variance_derivation.md` — the laws M1–M10. §10 is the composition peel.
6. `SESSION_2026_07_25_HANDOFF_9.md` §0b — the AMBIG record and its retraction.
7. `CALIBRATION_ARCHITECTURE.md` — count-zero-information, the three sources.
8. `docs/calibration/boundary_work_notes.md` — **the owner's own notes file.** Do not overwrite.

**Never reference `docs/calibration/archive/`.**

## 10. Tools

`scripts/debug/`: `pass0_oracle_bench.py` (the A/B), `pass0_error_table.py` (suite state in READS + the TRUST
`errQ1conf` view), `suite_dissect.py` (per-node table), `derive_4_boundary_classes.py`,
`message_variance_mc.py` (M1–M10 arbiter).

`scratchpad/` (untracked, **keep it**): `steps_1_to_4.patch` (**the restore**), `b1_boundary.py` (boundary
census by class + coverage), `b2_faces.py` (face-selector scoring), `p1_path.py` (exon→boundary→intron
hop-by-hop), `i1_introns.py` / `i2_intron_ablate.py` (intron factory + channel ablation), `v1_vertex.py`
(the ψ vertex), `t1..t7_*.py` (the largest-scenario dissection), `ambig_1..8_*.py`, `ss_cap_1..8_*.py`.

⚠ `ss_cap_3_replay.py` and `ambig_2_mechanism.py` are **bit-exact ψ replays** (`max|Δ| = 0`) supporting
per-channel ablation — the single most useful diagnostic in the set. Re-verify fidelity after any solver
change before trusting an A/B.

## 11. Invariants — preserve these

* The **`√2·σ_own`** DL pin-safety inequality: a message out-weighs a node's own belief only if it agrees to
  within `√2·σ_own`. Any change that breaks it is a regression even if the A/B improves.
* `Σ_c ρ_c·E_c = M` is an identity under the imputation premise, enforced at every relay hop and the combine
  with `_pin_v` semantics (a component the context does not supply is filled from the node's OWN density).
* `N` enters only as power (`τ_λ` Fisher or `1/n`), never as a composition vote.
* Variances in log-odds / `Var(log f_c)`, never on simplex fractions.
* The composition is ONE dof (single-λ message); the tilt θ is a SEPARATE dof (AMBIG only).
* **No magic numbers** — structural presence tests, derived quantities, exact limits only.
* **Pass-0 must be WEAK and CORRECTABLE.** An over-confident message that pins a node wrong is worse than a
  weak one slightly off. Prefer the under-confident option when unsure.

## 12. Vocabulary (owner's, use it)

RNA is **one species**. "Mature" and "nascent" are bookkeeping, not biology — only **SPLICED vs UNSPLICED** is
observable. A boundary can be an exon↔exon boundary that is *also* a splice junction, with RNA contiguous
across it while other RNA splices in or out. Both at once. The peel is about **spliced vs unspliced**, never
about mature vs nascent.
