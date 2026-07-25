# Session handoff — the `gdna300_ss_0.50_nrna_none_capture_on` deep dive: the count-zero-info wall, measured

**LIVE handoff. Read `docs/calibration/ROADMAP.md` first, then this.** Supersedes
`SESSION_2026_07_25_HANDOFF_9.md` (still the AMBIG record — read its §0b). Do NOT read `archive/`.
Date: 2026-07-25. Branch `calib-ambig-init-wip`, HEAD `3c9360ce` (M8: 0.0900 refit=0 / 0.0700 refit=1).

---

## 0. TL;DR

Root-caused the suite's single largest error scenario (`gdna300_ss_0.50_nrna_none_capture_on`, **1,358,610
error reads = 10.8 % of all suite error**). Two results, one landable-pending-a-decision and one that redraws
the pass-0/Phase-2 boundary:

1. **A real prior-free defect, found and fixed**: the intron factory's confident gDNA density is **excluded
   from the gDNA measurement stream** (`mg_own = where(struct_lock, …)` — measured: **0 of 572** factory
   introns seed it). Exons flanked on both sides by factory-solved introns therefore receive an unopposed
   spliced-RNA measurement. Fix wins **refit=0: 0.0900 → 0.0884, 18 better / 9 worse**, target scenario
   −5.7 %, and collapses the pathological stratum from 90 nodes/200 k reads to 12 nodes/41 k. **But refit=1
   is 6 better / 19 worse**, so it is NOT landed — see §5 for the decision needed.
2. **The wall is measured, not assumed.** On **87.6 %** of the error `_pin_v` cancels the reframe `r`
   exactly, so the delivered mode is the SOURCE's composition re-expressed in the destination's eff-lengths —
   and neighbouring nodes differ in composition by **|Δf_g| = 0.2400**. That is the count-zero-information
   wall in numbers: no prior-free message can beat it. Substituting oracle modes removes 67–75 %; scaling
   every precision ×10 removes **nothing**. So the residual here is Phase-2's, and now we know how much.

## 1. Why this scenario, and its physics

Largest single error scenario (`scripts/debug/pass0_error_table.py`). κ = 0.49921 ⇒ **the strand likelihood is
flat ⇒ τ_own = 0 at every node except introns**; no nascent RNA ⇒ every intron's true RNA is exactly 0;
capture ON ⇒ large enrichment cliffs. So the **intron factory + message propagation are the ONLY information
in the entire solve**, and every exon's answer is 100 % imported (self-solve = 0.49/0.51/0.54 everywhere —
the uninformed default).

Messages are doing real work: **self 2,247,710 → solved 1,358,610 reads**. Channel ablation (bit-exact ψ
replay, `max|Δ| = 0`): `-gdna` **+29.0 %**, `-lam` +5.6 %, `-factory` +2.7 %, `-tilt` +1.4 %, `-rna` **−3.1 %**.
The gDNA channel is the workhorse; the RNA channel is net harmful.

Error by class: exon single **903 k (66.5 %)**, exon AMBIG 155 k (11.4 %), exon+intron(x) 153 k (11.3 %),
boundary 102 k (7.5 %), intron 6 k (0.5 %). Signed error on exons **−0.1268** — the solve systematically
under-claims gDNA.

## 2. ⭐ THE DEFECT — the factory's gDNA authority never reaches the measurement stream

`bp_solver`: `mg_own = np.where(_struct, pg_own, 0.0)` seeds the gDNA MEASUREMENT stream **only at
`struct_lock` nodes**, while the RNA measurement stream is seeded at every spliced graft. Measured here:

```
struct_lock nodes (seed mg_own):     453   of which mg_own>0: 216
factory introns (tau_lam>0):         572   mg_own>0 among them:  0     <-- none
boundaries with spliced (seed RNA):  764
```

Partitioning exons by which channels arrived:

| stratum | n | ERR reads | share | mwae | signed |
|---|---|---|---|---|---|
| both gDNA + RNA | 414 | 943,417 | 69.4 % | 0.1772 | −0.1262 |
| gDNA only | 144 | 73,575 | 5.4 % | 0.1493 | **−0.0053** |
| **RNA only** | **90** | **199,973** | **14.7 %** | 0.2673 | **−0.2070** |
| NOTHING arrived | 6 | 32,732 | 2.4 % | 0.4569 | −0.4569 |

The **gDNA-only stratum is essentially unbiased** and the **RNA-only stratum is catastrophically
RNA-biased**. 65 % of RNA-only exons have a factory intron within 2 hops. The anchor case:

```
node 1075  62,195 reads  oracle 0.822  self 0.510  ship 0.305  ERR 32,164
           cm_g 0.00   cm_p+n 2.62     neighbours: intr*, boun, boun, intr*   (* = factory intron)
```

An exon flanked by two confidently-solved gDNA-rich introns gets **zero gDNA weight** and a spliced-RNA
measurement at weight 2.62. Distance is not the problem — 90.4 % of exon error is 2–3 hops from an intrinsic
source. It is a **stream-membership** problem.

**The fix** (`scratchpad/` EXP-A/B): seed `mg_own` from density-level authority, not just structural
certainty — `where(struct_lock | tau_factory > 0, pg_own, 0)`. This is principled: `pg_own =
1/(Var(log f_g) + 1/n)` is the honest precision of a gDNA DENSITY claim, and at f_g → 1 the Jacobian
`(1−f_g)²/τ_λ → 0` so it degrades smoothly to the pure count — exactly a struct_lock anchor in the limit. The
factory is a *density deconvolution*, which is what the measurement stream is for.

## 3. ⭐ THE WALL — 87.6 % of the error is imputation-premise-limited

**⚠ This corrects an intermediate finding of this session.** A per-edge decomposition suggested the gDNA
error was the transport ratio `r` (substituting the oracle capture step cut the per-edge claim error 5–8×,
and `r_M/E_g` scored 0.808 → 0.330 against the true capture step). **That measurement was of the PRE-pin
claim and is an artifact.** `_pin_v` rescales each message to the node's own mass, so when all components are
supplied the scaling is ∝ 1/r and **`r` cancels exactly**:

```
stratum                    n       ERR   share    is `r` cancelled by the pin?
both gDNA+RNA supplied   492 1,038,367  87.6%    YES — r is irrelevant there
gDNA only                163    73,782   6.2%    NO
RNA only                  12    40,823   3.4%    NO
```

This is also why §4.2's oracle-`r` substitution bought ~0, and it is consistent. What sets the mode on that
87.6 % is the imputation premise itself:

```
|f_g(dst) − f_g(src)| = 0.2400        |log f_g(dst) − log f_g(src)| = 0.3625      (exon edges, oracle)
```

and the delivered share misses by 0.2891–0.2995. **Neighbouring nodes genuinely differ in composition by
0.24, so a composition-imputation message cannot do better.** That is the count-zero-information principle
as a number.

**The ceiling** (bit-exact ψ replay, oracle modes substituted one channel at a time, precisions as shipped):

| condition | shipped | oracle mo_g | oracle mo_R | oracle both | **×10 precision** | nomsg |
|---|---|---|---|---|---|---|
| gdna300 ss0.50 none capON | 1,281,156 | **420,459** | 984,216 | 331,805 | 1,269,460 | 2,247,710 |
| gdna300 ss0.50 present capON | 1,266,170 | **311,171** | 1,165,262 | 333,214 | 1,261,139 | 2,371,171 |
| gdna100 verystrong | 1,010,313 | 692,007 | 981,299 | 876,224 | 1,068,948 | 1,402,103 |

67–75 % is available **to a better mode only**, and **no re-weighting can reach it** (×10 precision moves
nothing). Since the mode is premise-limited, that 67–75 % is the hyperprior's.

## 4. ⛔ REFUTED this session — do not re-run

| item | verdict |
|---|---|
| **Left/right message disagreement as a DL second study** (to price mismatch where τ_own = 0 and DL is inert) | **Refuted.** corr with the delivered mode error = 0.028 / 0.032 / −0.116 / 0.364, and mwae is FLAT or FALLS with the gap band. Both sides share the same `r`/graft bias, so they agree while both are wrong — the classic correlated-estimator failure. |
| **`r_M/E_g` (or any better capture-step estimate) as the frame ratio** | Scores far better against the true capture step (0.808 → 0.330 all-edges) **but is moot**: the pin cancels `r` on 87.6 % of the error. Also loses at capture-OFF (0.322 vs 0.587). |
| **More `_RHO_ITERS`** (2 → 4 → 8) | No effect: 1,281,156 → 1,282,539 → 1,281,156. The combine's lazy ρ_tot is converged. |
| **Routing the factory to exactly one stream** (EXP-C: composition stream = strand only) | **Bit-identical to EXP-B.** The double-count hypothesis for the refit=1 regression is refuted; the cause is elsewhere. |

## 5. ▶ THE DECISION NEEDED, and what is next

**EXP-B is validated at refit=0 and regresses at refit=1:**

| arm | refit=0 | | refit=1 | |
|---|---|---|---|---|
| `m8` (HEAD) | 0.0900 | — | 0.0700 | — |
| **`expB`** | **0.0884** | **18 better / 9 worse / 5 flat** | 0.0688 | **6 better / 19 worse / 7 flat** |

The refit=1 aggregate improves only because of one large win (`gdna100_verystrong` 0.3203 → 0.2548) while 19
conditions regress, concentrated on stranded × capture-ON. Every prior landing in this project won BOTH
refits, so **it is not landed** and the tree is clean at M8. The judgement call is whether a pass-0
correctness fix should be gated on a refit that is itself the acknowledged next thing to fix (ROADMAP: "the
hyperprior refit still regresses unstranded-capON, and that is now the single binding constraint"). My
recommendation: **re-test EXP-B immediately after the hyperprior work**, because the regression is an
interaction with the component being replaced.

**What this scenario says about diminishing returns.** Its precision is already honest — `errQ1conf` = 4.1 %
(share of error on the most-confident quartile) and the gDNA channel's reliability `z2` = 0.4–1.5 across
precision bands. So pass-0 here is *not* confidently wrong; it is honestly uncertain about something it
cannot know. That is exactly the substrate the hyperprior needs. **The remaining prior-free work in this
scenario is EXP-B (≈ 6 %) and the 6 "NOTHING arrived" exons (2.4 %); the other ~90 % is Phase-2's.**

⚠ **Where precision is NOT honest** (from `pass0_error_table.py`, and this is the more important Phase-2
risk): **stranded capture-OFF** scenarios put **58–76 %** of their error on the most-confident quartile, and
**introns are 90.2 %** confident-error. Introns are small in absolute error (209 k, 1.7 % of suite) but they
are precisely where gDNA density is measured for the hyperprior. **That, not the big unstranded scenarios,
is what would poison Phase 2.**

## 6. Tools

`scripts/debug/pass0_error_table.py` — the suite state of play in READS with the trust (`errQ1conf`) view.
`scratchpad/t1_char.py` (characterize + ψ ablation, bit-exact), `t2_strata.py` (channel-arrival strata +
distance-to-anchor), `t3_streams.py` (who seeds which stream), `t4_gdna_transport.py` (per-edge transport
decomposition + candidate frame ratios), `t5_lr_gap.py` (the refuted left/right statistic),
`t6_ceiling.py` (**the oracle-mode ceiling — use this before chasing any mode work**), `t7_pin_check.py`
(does the pin cancel `r`; the imputation-premise floor).
