# Session handoff — P1d is LANDED, and the mechanism it was named for is refuted

**This is the LIVE handoff. START HERE.** Date: 2026-07-26. Branch `calib-ambig-init-wip`.
Supersedes `SESSION_2026_07_26_HANDOFF_13.md` (still the composition-peel / message-packet record).

Suite **0.0871 (refit=0) / 0.0684 (refit=1)** — deliberately UP from 0.0855 / 0.0671, bought against a
**27.5 % cut in confidently-wrong mass**. Gates all green.

---

## 0. The one-paragraph version

P1d is built and on by default. But the term it was scoped as — *"a ~100 bp junction speaking for a ~2,100 bp
exon, a 12–21× extrapolation"* — **is not the mechanism.** The residual is flat in the extrapolation ratio and
flat in exon length. What the graft actually gets wrong is that its RNA claim is a **lower bound used as an
equality**, and it fails hardest at **transcript termini**, where RNA does not flow through the seam. One
free structural bit splits the error **≥30×**. That reopens `ROADMAP.md` §6, which had deferred TSS/TES in the
region map as "expected to be low-impact"; it is 62–72 % of the graft's premise error.

## 1. What landed

`enrichment_frame.graft_premise_logvar` — the graft's premise variance, measured per-edge from the
disagreement of an exon's **two flanking seams**:

```
    v_premise = max(0, d² − v_a − v_b) / 2 ,      d = log(ρ_μ^a / ρ_μ^b)
```

* **Method of moments**, the third use of it in the calibrator after `κ` and both strand overdispersions —
  a *fitted* quantity, not a tuned constant. No coefficient to choose: the `/2` is the **measured**
  independence of the two seams (`corr(log φ_L, log φ_R)` = 0.017 raw / −0.036 post-Poisson at capture-OFF,
  over 14 conditions — the earlier "≈0.95 correlated" claim is definitively dead), and the truncation at 0 is
  the method's own.
* **Per-edge where a second seam exists (47.7 %); the same estimator POOLED where one does not.** The pooled
  fit is on the measurably better-behaved half (premise 0.48 vs 0.60 one-seam, 2.69 no-seam), so it
  under-charges — the safe direction.
* Both fluxes are **lifted into the destination's frame first**, so a capture step common to the two seams
  cancels out of `d`; and each seam's noise carries its spliced **COUNT** (never the mass) **⊕ its lift's own
  scale sampling** (M5's source leg — the destination's leg is common to both lifts and cancels).
* Applied to the **WHOLE** grafted RNA claim, in **both** code paths. The premise is about the SUM: `Var(log φ)`
  is flat in the spliced share `w_μ` (2.02 → 1.83 across `w_μ` 0.47 → 1.00) while `Var/w_μ²` swings 5×.
* `RIGEL_GLV_OFF=1` ablates it — verified bit-identical to the pre-P1d path, 32/32.

| | mwae r0 / r1 | **confidently-wrong** | z2 exon-single | z2 exon-AMBIG | **z2 ALL** |
|---|---|---|---|---|---|
| pre-P1d (`pk2`) | 0.0855 / 0.0671 | 1,200,951 | 8.60 | 92.10 | 15.55 |
| **P1d** (`p1d`) | 0.0871 / 0.0684 | **870,245 (−27.5 %)** | **5.68** | **65.28** | **11.12** |

10 better / 7 worse / 15 flat at refit=0. Arms in `/tmp/pass0_oracle_bench.tsv`: `pk2` (pre), `p1d` (HEAD).

## 2. ⭐ What the measurement refuted — read before touching P1d again

| §5.2 / §3 claim | verdict |
|---|---|
| the mechanism is the 12–21× **extrapolation** | ⛔ **REFUTED.** Flat in the extrapolation ratio (1.13× over a 6.7× span) and flat in exon length |
| `Var(log φ)` spans **40×** across the count range, so the shape is a count law | ⛔ **a PROXY.** Flat inside every structural class (junction-only excess 0.036 → 0.073 across a 200× count range); the count tracks *terminus membership* |
| "Poisson explains almost none of it" | ⚠ **partly.** It was one term subtracted from a five-count statistic, using the spliced **MASS** (~2× the count). Full trigamma subtraction: Poisson explains **16–24 %**, residual **0.485** |
| the true `Var(log φ)` is 0.583 / 0.574 | ⚠ contains the ORACLE capture step's own two gDNA counts. Variance components across independent libraries put the counting noise at 0.1414 against the 0.0912 that was subtracted |
| the estimator over-states by 1.23× | ⛔ a **frame** conversion, not a bias: `E1` lives in μ-space, the target in claim-space, frame factor `E[w_μ²]` = 0.845–0.899. Frame-matched it is **1.04×** |
| low count ⇒ **minor isoform** ⇒ unrepresentative | ⛔ in this suite the driver is **TRANSCRIPT TERMINI**, not minor isoforms. Within the minor-isoform class the count coefficient dies (t = 0.36) |
| count merely tracks **expression** (the confound) | ⛔ also refuted — under honest centring the expression coefficient goes negative and insignificant |

**The structural partition (all 14 capture-OFF conditions, 4,763 edges):**

| class | edges | med φ | Var(log φ) | share of squared deviation |
|---|---|---|---|---|
| **A** boundary carries no junction of this exon (100 % lie inside an annotated exon, split by another transcript's TSS/TES) | 3.7 % | 0.84 | 1.66 | 9.8 % |
| **B** boundary is a transcript **TERMINUS** | 17.1 % | 0.67 | 1.99 | **61.9 %** |
| **C** minor-isoform junction | 43.5 % | 1.25 | 0.24 | 20.3 % |
| **D** sole junction into the exon | 35.7 % | 0.98 | 0.16 | 8.0 % |

## 3. ⛔ DO NOT RE-RUN — settled this session

| item | verdict |
|---|---|
| the **MODE** fix (replace the graft's claim with the dominant seam's flux, `max`) | Right physics, wrong regime. Capture-OFF **14/14 better, 0 worse**; capture-ON **8 worse** (mwae +0.0039). Under capture the graft ALREADY over-states (median φ = 2.45 — an M8 frame problem), so no bound-tightening can help. Framing the two seams first halves the damage (+0.0098 → +0.0047) but does not remove it. **Not landed** |
| `min` / `sum` / `2×` variants of that mode fix | `min` 0.514 (worse than either bound alone), `sum` +0.61 nats of over-statement, `2×` 0.503 — each fails exactly as its algebra predicts, which is what makes the `max` result structural rather than a max-of-two selection artifact |
| omitting M5's lift noise from the MoM subtraction | Scores better on trust (818 k vs 870 k confidently-wrong) and is **dishonest** — MoM books as premise error every noise it fails to subtract. Costs +0.0010 mwae to do right. Do not "recover" it |
| the two-junction gap as **≈0.95 correlated** (HANDOFF_10 §4 family) | Dead. Measured directly against the same exon: `ρ` = 0.017 raw / −0.036 post-Poisson at capture-OFF, implied divisor 1.97–2.07. (Under capture `ρ` = 0.30, divisor ≈1.39 — a correction P1d does not currently make) |
| `w_μ` as a per-edge reliability signal | R² = **0.00** on all 14 capture-OFF conditions, and not vacuously: `sd(log w_μ)` = 0.49, so it has real dynamic range. A true null |
| `μ/(M_exon/E_r)` ("vsEXON") as the covariate | Looks strong (R² 0.21–0.80) but its R² rises exactly as gDNA falls — the signature of a shared-denominator artifact, not information |
| `RIGEL_XVAR` as a size estimate | It reaches only `_transport`, not `_relay` — 3,636 of 5,454 graft firings (67 %, **not** the 50 % first computed: `_RHO_ITERS = 2` and `_transport` runs 4× per sweep while `_relay` runs 2×) |
| "M7 is inert at graft destinations" | **False.** M7 is live at 62–75 % of them and fires at 21–50 %, including on unstranded data. §2.2's `max()` no-double-count argument is untouched; only the belt-and-braces sentence was wrong |

Everything in HANDOFF_13 §5, HANDOFF_11 §8 and HANDOFF_10 §4/§9.2 still holds.

## 4. ▶ NEXT

1. **P1e — the conservation SURPRISE**, and P1d's own residual regression points straight at it. The
   regression is localized to **unstranded × gDNA-rich × capture-OFF** (`gdna100_ss0.50_none_capOFF` +0.0406,
   `gdna300_ss0.50_*_capOFF` +0.012) — which is exactly P1b's population, where the RNA claim is a near-exact
   measurement (36.6453 claimed vs 36.6734 true) and the **gDNA** claim is 47× too big. P1d damps the RNA arm
   there, which is the wrong arm. `claim/obs` is the prior-free observable that says which arm failed
   (`PASS0_FINISH_PLAN.md` §P1b: 1.00–1.01× ⇒ excellent, 1.26–1.33× ⇒ bad), and nothing uses it as evidence.
2. **P1g — put TSS/TES in the region/boundary map** (un-deferred from `ROADMAP.md` §6). P1d prices this bit
   blind, through the two seams' gap; naming it directly prices it exactly, and it is 62–72 % of the error.
   Reaches the index schema (`boundary_df` has only `boundary_id`/`ref_name`/`position` today; `t_df` has the
   transcript spans, so the bit is derivable at index-build with no new input).
3. **P3 — AMBIG exon over-confidence.** `z2|Q1` 183 → 92.1 → **65.3**, still the largest single defect.

## 5. HOW TO RUN EVERYTHING

Unchanged from `SESSION_2026_07_26_HANDOFF_13.md` §6, plus one flag. Gates:
`pytest tests/ -q` = 1248+7 pass / 2 xfail / 2 xpass · `ruff check src/ tests/ scripts/` ·
`python scripts/debug/message_variance_mc.py` = 0 failures.

| flag | what it ablates |
|---|---|
| **`RIGEL_GLV_OFF=1`** | **P1d** (bit-identical to the pre-P1d path, verified 32/32) |
| `RIGEL_S2T_OFF=1` | both cliff terms (M5 + M8) |
| `RIGEL_RNAMEAS_OFF=1` | the RNA measurement ψ factor |
| `RIGEL_M12_MSG=1` | the per-message pin → the weighted rescale |
| `RIGEL_XVAR=<v>` | the flat probe (superseded by P1d; reaches only `_transport`) |

New in `scratchpad/`: `x2_poisson.py` (the five-count Poisson decomposition), `x3_share.py` (the structural
partition + the regressions), `x4_estimator.py` (the estimator audit + the four claim forms),
`x2_graft_mirror.py` / `x3_graft_reuse.py` (the code-path replay), and the `v*_x2_audit_*.py` / `x5..x7_*.py`
verification scripts. **Keep them** — they are what refuted §5.2.

## 6. Invariants — one addition

All of HANDOFF_13 §7, plus: **a method-of-moments fit must subtract every noise source the model already
knows**, or it books them as the effect being fitted. P1d subtracts the spliced count *and* M5's lift
variance; leaving the latter out scores better and is wrong.
