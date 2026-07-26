# Session handoff — the largest-scenario deep dive + the INTRON fix (landed)

**LIVE handoff. Read `docs/calibration/ROADMAP.md` first, then this.** Supersedes
`SESSION_2026_07_25_HANDOFF_9.md` (still the AMBIG record — read its §0b). Do NOT read `archive/`.
Date: 2026-07-25. Branch `calib-ambig-init-wip`. **HEAD now 0.0898 refit=0 / 0.0697 refit=1** (M8 +
the intron λ gate of §6, which is the one thing LANDED this session).

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

## 6. ⭐ THE INTRON RESULT — a direct measurement was being out-voted by an imputation (LANDED)

Owner's model, and it is right: the density deconvolution delivers the **λ axis only** (gDNA vs RNA-total)
plus a precision — it does NOT assign the tilt. So a **single-strand intron** has θ structurally locked and
λ measured ⇒ a COMPLETE self-solve needing no messages; an **AMBIG intron** knows its split but not which
strand the RNA sits on ⇒ it needs messages for **θ only, never for λ**.

**What was happening.** Messages made introns worse in **16 of 16** scenario × DOF strata — the only class
where that happens (suite-wide self 0.0103 → solved 0.0133). Channel ablation isolates it to the **λ
composition message** (−12 % … −37 % when removed); `-gdna` and `-tilt` do nothing to introns.

**Why DL did not stop it.** M7's safety property lets a message out-weigh the own belief when it agrees to
within `√2·σ_own`. At a real intron:

```
factory τ_own = 0.254  ⇒  σ_own = 1.98,  √2·σ_own = 2.81
λ message |G_lam| = 2.08  <  2.81      ⇒  DL correctly PERMITS it
ψ then weights c_tau = 0.6151 vs the factory's 0.2540  ⇒  the imputation gets 70.8 %
```

DL is behaving exactly as designed. `σ_own` is large because the factory's λ variance is **conservative**:
vertex-free calibration `z2 = E[(λ_self−λ_true)²]/E[1/τ]` = **0.10–0.32** (1.0 = honest), i.e. under-confident
by 3–10×. ⚠ **Correction to an earlier reading in this session and to the owner's hypothesis: the NB factory
is NOT over-confident.** The 20–26 figures first reported were a logit-clipping artifact at `f_g = 1.0`.

**The fix (LANDED).** A node holding a *direct density measurement* of its own composition does not accept an
*imputed* composition message on that same axis — structural, no threshold, the same kind of presence test as
the strand gates. Scoped to the FACTORY half of τ (`NodeInit.tau_factory`, new), never the strand half, since
`I_strand` is itself a composition vote. **The tilt message is untouched**, which is what AMBIG introns need.

| | refit=0 | | refit=1 | |
|---|---|---|---|---|
| M8 | 0.0900 | — | 0.0700 | — |
| **+ intron gate** | **0.0898** | **10 better / 0 worse / 22 flat** | **0.0697** | **13 better / 3 worse / 16 flat** |

Intron error **244,504 → 213,636 reads (−12.6 %)**; per-scenario −37 % / −12 % on the two capture-ON
gdna300 conditions. Aggregate movement is small because introns are only 1.7 % of suite error — **the value
is substrate quality for Phase 2**, since introns are where gDNA density is measured.

**Refuted while getting there** (do not re-run): replacing the factory's grid-variance precision with the
CURVATURE at the mode — which is what `density_factor_precision`'s own docstring specifies — makes it
**worse** (intron error 5,882 → 7,099), because for `f_g → 1` introns the factor's mode sits on the λ-grid
edge and the second difference is then taken in the tail. And scaling `tau_fac` alone saturates
(6,872 → 5,753 at ×100) because ψ weights the `intron_prior` array while `v_own` reads `tau_lam` — the same
quantity down two paths that then disagree. **That split is a latent defect worth fixing properly.**

**Still open on introns:** after the λ gate, the residual message damage is the relayed **RNA measurement**
(4,341 → 3,189 when ablated, vs 3,181 for no messages at all). The same "direct measurement outranks
imputation" argument applies to it — but `th_prec = cm_p + cm_n`, so gating it would also kill the tilt
channel AMBIG introns depend on. Needs the two separated first.

## 7. ⭐⭐ THE ROOT CAUSE — an exon has no RNA claim across a splice junction (LANDED)

**Owner's call, and it was right: §6's intron λ gate was a band-aid.** The messages themselves are bad. They
travel through exon↔intron boundaries, are contaminated with the exon's RNA, over-estimate RNA, and then
poison introns — which the density deconvolution had solved accurately. Fix the message, not the victim.

**Measured, on `nrna_none`, where the oracle RNA at an exon|intron boundary AND inside the intron is EXACTLY
0** (`scratchpad/p1_path.py`):

| node kind | n | oracle ρ_R | relayed ρ_R | oracle ρ_g | relayed ρ_g | ratio |
|---|---|---|---|---|---|---|
| EXON (source) | 681 | 5.5050 | 13.4281 | 27.5491 | 25.2997 | 0.96 |
| exon\|intron BOUNDARY | 1145 | **0.0000** | **0.5280** | 14.8641 | 13.8969 | 0.98 |
| INTRON | 571 | **0.0000** | **0.0041** | 0.0524 | 0.0514 | 0.99 |

```
THE PEEL, exon -> exon|intron boundary
dir      n   ρ_R(exon)  reframed  peel amt  post-peel  oracle@bnd   peel==0
L->R   572     18.3554    7.8185    1.2189     2.0373      0.0000       17%
R->L   573     15.4216    5.9309    1.4768     1.4256      0.0000       14%
```

**gDNA transports essentially perfectly across every hop (0.96–0.99). The RNA channel alone is the
contaminant.** The peel removes only ~16 % of the reframed RNA and **does not fire at all on 14–17 % of
edges**.

**Why.** The peel subtracts the MEASURED spliced flux as a proxy for "the mature that leaves". But across an
active junction **mature RNA splices out by definition**; what continues unspliced is nascent read-through,
and an exon cannot measure how much of its own RNA that is. So the exon's RNA density is not a claim about
the other side at all — the component sets differ, which is exactly the refusal
`enrichment_frame.gdna_fallback_admissible` already describes. The honest emission is **none** (density AND
precision), not a partial subtraction.

**A/B.** `_junc` = the boundary carries spliced mass on either face; on such a peel edge the exon emits no
RNA claim:

| arm | refit=0 | | refit=1 | |
|---|---|---|---|---|
| M8 | 0.0900 | — | 0.0700 | — |
| **+ junction rule (LANDED)** | **0.0892** | 16 better / 16 worse | **0.0677** | **23 better / 6 worse / 3 flat** |
| + junction rule + §6 intron gate | 0.0888 | 18 / 14 | 0.0672 | 25 / 3 / 4 |

**The §6 intron gate was REMOVED** (with `NodeInit.tau_factory` and its tests) now that the root cause is
fixed — owner's directive. Note the trade recorded above: keeping it buys a further 0.0004 / 0.0005 and cuts
refit=1 regressions from 6 to 3. Also note intron error is 217,301 reads without it vs 208,599 at the M8
baseline and 183,148 with it — **introns alone are marginally worse than baseline**, while the suite and the
boundaries are clearly better. That is a live decision, not a settled one.

**The honesty gain is the real headline.** Boundary `errQ1conf` (share of a class's error sitting on the
most-confident quartile) collapses: **single 12.8 % → 4.3 %, AMBIG 17.7 % → 1.5 %**. The contaminated RNA
messages were producing *confidently wrong* boundaries, which is precisely what would poison the Phase-2
hyperprior fit. Suite total: **12,561,367 → 12,447,440 error reads**; stranded × capON 0.0565 → 0.0520.

**Residual, and the next step.** The leak is reduced, not closed: relayed ρ_R at the boundary 0.5280 → 0.2381
against an oracle of 0. What is left is the **14–17 % of edges where `_junc` is false** — the gate is
OBSERVATIONAL (spliced mass was measured) where it should be STRUCTURAL (the annotation says a junction is
there). `build_node_statics` already computes `mrna_active_pos/neg` and it is unused (§4.1 flags this as "a
partial mitigation available without the index work"). That is the next thing to try.

Also open, from §6: ψ weights the `intron_prior` ARRAY while `v_own` reads `tau_lam` — the same quantity down
two paths that then disagree (proved by the ×100 saturation test). Worth unifying.

## 8. ⭐⭐⭐ THE BOUNDARY RULE — an exon has no RNA claim where its MATURE cannot continue (LANDED)

**Owner's refinement of §7, and it corrects an over-reach in it.** §7 gated on *junction presence* (spliced
mass observed). That is the wrong predicate in both directions. The right one is **structural and already
computed**: `statics.mrna_active_s` = the exon bit set on BOTH flanks, i.e. *mature runs contiguously past
this seam on strand s*.

* `mrna_active_s` **TRUE** — an exon↔exon seam, including the owner's alternative-isoform case
  (`TA+ (1000,2000),(9000,10000)` splices 2000→9000 while `TB+ (1000,10000)` reads straight through). **Mature
  IS contiguous here**, so the exon's RNA legitimately transports and the measured peel removes what splices
  away. §7's junction gate wrongly silenced these: `exon|exon` junction boundaries had **`c_tau = 0.00`** —
  no composition message at all — and they are 24–36 % of boundary error.
* `mrna_active_s` **FALSE** — an exon↔intron seam. Mature splices out by definition; only nascent read-through
  continues and the exon cannot measure how much of its own RNA that is ⇒ **no RNA claim** (density AND
  precision). This also fires on the **14–17 % of exon↔intron seams whose junction was never observed**, which
  the observational gate missed entirely.

| arm | refit=0 | | refit=1 | |
|---|---|---|---|---|
| M8 | 0.0900 | — | 0.0700 | — |
| §7 junction-presence gate | 0.0892 | 16 better / 16 worse | 0.0677 | 23 / 6 / 3 |
| **§8 `mrna_active` rule (LANDED)** | **0.0884** | **16 better / 11 worse / 5 flat** | **0.0678** | **23 better / 4 worse / 5 flat** |

Better aggregate at refit=0 and **fewer regressions at both**. `exon|exon` regains its composition message
(`c_tau` 0.00 → 12.28) and the intron leak tightens further (relayed ρ_R at an exon|intron boundary
0.5280 → 0.2381 → **0.2030**, at the intron 0.0041 → 0.0028 → **0.0022**, against an oracle of exactly 0).

**Suite: 12,561,367 → 12,343,147 error reads (−1.7 %).** The honesty gain at boundaries is the headline:
`errQ1conf` **single 12.8 % → 6.6 %, AMBIG 17.7 % → 2.2 %** — the contaminated RNA messages had been
manufacturing *confidently wrong* boundaries, precisely what would poison the Phase-2 hyperprior fit.

### 8.1 The boundary census — what each class has to work with

`scratchpad/b1_boundary.py`, per flanking pair × junction, with the ORACLE composition of the crossing mass
and which of the three deconvolution sources reaches it (strand τ, factory-intron neighbour, exon gDNA):

| class (gdna300 ss0.99 present capON) | share of bnd ERR | oracle f_g | self → ship | τ>0 | facNbr | c_tau |
|---|---|---|---|---|---|---|
| exon \| intron, junction | 52.4 % | 0.6257 | 0.0403 → 0.0542 | 89 % | 100 % | 10.63 |
| **exon \| exon, junction** | **26.6 %** | 0.5621 | 0.1098 → 0.0992 | 49 % | **3 %** | 12.28 |
| exon \| intron, no junction | 11.5 % | 0.9764 | 0.0553 → 0.1010 | 90 % | 99 % | 4.26 |

**`exon|exon` boundaries are structurally starved**: both flanks are exons, so there is **no intron factory
within reach (3 %)** and on unstranded data no strand either — the owner's "the boundary is the crux" case
with the fewest sources. That is where to look next.

Note also two classes where **messages still make the boundary worse** (`exon|intron` junction 0.0403 → 0.0542;
no-junction 0.0553 → 0.1010) — the peel is fixed but something else is still contaminating those.

### 8.2 Coverage of the owner's untested regimes (measured, for the simulator question)

* **Unspliced RNA density > spliced RNA density at a junction boundary** — present but only in the tail:
  fraction > 1 is **28.7 % / 1.7 % / 9.4 % / 14.9 %** across scenarios, median ratio 0.12–0.48. The
  *very* high-nascent regime the owner describes is NOT well covered.
* **Inverted enrichment (probe over the boundary, boundary MORE enriched than the exon)** — boundary/exon
  oracle gDNA-density ratio median 0.13–0.83, but **boundary > exon on 3 % / 12 % / 22 % / 28 %** of edges. So
  it exists, thinly, and only as a by-product rather than by design.

Both argue for purpose-built scenarios rather than relying on the tail of the current suite.

## 9. ⭐ THE BOUNDARY REFRAME — what is validated, what is irreducible, and THE REMAINING TASK

**Owner's correction, and it invalidates §8's REASONING even though §8 won its A/B.** There is no difference
between mature and nascent RNA — that distinction is manufactured for our own bookkeeping. Only SPLICED vs
UNSPLICED is observable. So "mature cannot continue past an exon→intron seam" is not physics; RNA crosses
either way. `mrna_active_s` is therefore **not** a statement about what can cross. §8's gate is compensating
for a broken peel, and it must be revisited once the peel is fixed (§9.3).

**The owner's model of a junction boundary — two sides, two frames.** The side facing the EXON has a total
density INCLUDING the spliced fragments and shares the exon's content; the side facing the INTRON EXCLUDES
them (they have spliced away) and shares the intron's content. Each side's ρ_tot ratio to its neighbour should
therefore be a pure CAPTURE step. And because the boundary sits on an enrichment slope, **the peel must be
done by COMPOSITION, not by subtraction.**

### 9.1 What is VALIDATED (`scratchpad/b2_faces.py`)

The two-sided reframe IS the shipped design (`node_total_density`'s docstring states it), and **the
intron-facing side is correct**: excluding the spliced there scores 0.077–0.924 against the true capture step
while including it scores 0.499–1.875 — 3–6× worse. Confirmed, keep.

### 9.2 What is REFUTED — do not re-run

| item | verdict |
|---|---|
| Making the face selector STRUCTURAL (`faces an exon`) instead of observational (`carries spliced mass`) | **No effect.** The two agree on **94.6–96.0 %** of junction boundaries; the structural version scores marginally WORSE (0.479→0.493, 0.269→0.275, 1.139→1.145, 0.283→0.327). |
| Folding the spliced into the exon-facing density differently (`ρ_u + S/E_r`, `ρ_u·(M+S)/M`, unspliced-only) | **No variant wins** across regimes (0.494/0.413/1.155/0.576 vs shipped 0.577/0.396/1.278/0.424). And `E_spl/E_r = 1.000`, so the eff-length is not the culprit. |
| Blaming the exon-side error on a contaminated `f_g` | **Refuted.** Substituting the ORACLE `f_g` into `node_total_density` does not help: 0.577→0.614, 0.396→0.320, 1.278→1.239, 0.424→0.725. |

### 9.3 What is IRREDUCIBLE, and what actually remains

The exon-facing residual is **0.4–1.3 nats and survives perfect inputs**. That is not a formula error — it is
§4.2's confirmed physics: a boundary samples a ~`fl_mean` window around a point while an exon samples its
whole interior, so when probes sit mid-exon the ends are depleted and **the "capture step" between a point and
a region average is genuinely large and scattered** (§4.2 measured p25–p75 = 0.55–2.9 at verystrong). Without
probe placement it cannot be recovered.

**But it mostly does not matter, and that tells us exactly where to work.** §3 established that `_pin_v`
cancels `r` on **87.6 %** of the error — wherever all components are supplied. The one place scale does NOT
cancel is a **subtraction**, and the peel is the only subtraction in the solver:

```
    tp = max(tp − spl_*_f[df], 0)        an ABSOLUTE density subtracted from a differently-framed one
```

**So the owner's prescription is precisely right and is the remaining task: peel by COMPOSITION.** The
boundary's own solved state supplies both parts in ITS OWN frame, so enrichment cancels inside the ratio:

```
    w_continue = ρ_ν / (ρ_ν + ρ_μ)          ρ_ν = (1−f_g)·U/E_r  (solved),  ρ_μ = S/E_spl  (measured)
    message to the intron side = ρ_R(exon) · r · w_continue        a scale-free FRACTION, not a difference
```

This replaces the only scale-sensitive operation in the message path with a scale-free one, and it is what
makes the exon-side slope scatter harmless rather than load-bearing.

### 9.4 THE PLAN (owner asked for next steps)

1. **Implement the composition peel** — `w_continue` from the boundary's own solved composition, applied
   multiplicatively; delete the subtraction. Both sites (`_relay`, `_transport`).
2. **Derive its precision.** A ratio of two solved quantities: `Var(log w)` by the delta method from the
   boundary's own `Var(λ)` and the spliced count `1/n_spl`. MC-validate as M10 in
   `scripts/debug/message_variance_mc.py`, the way M1–M9 were.
3. **Re-test §8's `mrna_active` gate.** Its justification is void (§9 opening); if the composition peel makes
   it unnecessary — which is the expectation — remove it. Do not leave a compensation in place next to a fix.
4. **Per-condition A/B at refit=0 AND refit=1**, then certify per boundary class: `exon|exon` with and without
   a junction, `exon|intron` with and without, retained-intron, and both enrichment orderings.
5. **Then** revisit `exon|exon` boundaries, which §8.1 shows are structurally starved (no factory within
   reach on 97 %, no strand when unstranded) — they need a source, not a better transport.

**No new simulations are required for 1–4.** §8.2's coverage gaps (very-high-nascent, inverted enrichment)
matter for CERTIFYING step 4 across regimes, not for deriving or landing the fix.

## 10. ⭐ THE BOUNDARY DEFAULT AND THE TWO-STEP — measured; the blocker is ψ's vertex

**Owner:** "intronic unspliced fragments are gDNA until proven otherwise" — the default at an exon↔intron
boundary should be 100 % gDNA at low/zero precision, not 50/50. Otherwise it is a **two-step solve**: the
intron message (or the strand) sets the gDNA-vs-unspliced-RNA level, and only then does the composition scale
work.

### 10.1 Both halves are one estimator, and both are measured right in their regime

The continuing share `w = ρ_ν/(ρ_ν + ρ_μ)`, scored against truth (mass-weighted |Δ|):

| condition | `w` TRUE | from the boundary's SELF-solve | **from the far INTRON** | "gDNA default" (w=0) |
|---|---|---|---|---|
| gdna300 ss0.99 present capON | 0.391 | 0.110 | 0.192 | 0.391 |
| gdna300 ss0.50 none capON | **0.000** | **0.628** | 0.294 | **0.000** |
| gdna100 verystrong | 0.413 | 0.156 | **0.126** | 0.413 |
| gdna300 ss0.99 present capOFF | 0.168 | 0.070 | **0.052** | 0.168 |
| gdna100 ss0.50 present capON | 0.403 | 0.160 | **0.100** | 0.403 |

The owner's two prescriptions are complementary and both correct: **w = 0 is EXACT with no nascent RNA**, and
**the intron sets the level when there is**. And the intron-sourced share *is* the two-step — it returns ~0
when the intron's density deconvolution finds no RNA and rises when it does, reaching the safe default as a
LIMIT of a measurement rather than asserting it.

### 10.2 ⛔ But both wirings LOSE — do not re-run

| arm | refit=0 | refit=1 |
|---|---|---|
| HEAD (`expF`) | **0.0884** | **0.0678** |
| own belief → 100 % gDNA at zero precision (`gdna_fallback_admissible`, finally wired) | 0.0906 (3 better / 15 worse) | 0.0690 |
| two-step: composition peel with the **intron-sourced** share | 0.0893 (10 / 21) | 0.0690 (7 / 25) |

* The **fallback default** loses for a mechanical reason worth recording: zeroing `op`/`on` breaks `_pin_v`'s
  partial-claim semantics (`sp = where(pp > 0, p, op)`), so a gDNA-only message is renormalised onto the
  node's whole mass and asserts `f_g → 1` — the over-confidence that operator's own docstring warns about.
  (`gdna_fallback_admissible` had been designed, documented and unit-tested in `enrichment_frame` and never
  called; it is still never called.)
* The **two-step** loses despite a better average share, because **its residual is an OVER-claim**, and RNA
  over-claim is the damaging direction — the one thing every study this session has agreed on. HEAD's
  conservative full peel errs the safe way and wins.

### 10.3 ⭐ THE BLOCKER, and it is one number

The intron-sourced share gives **0.294 where the truth is 0.000** in `nrna_none`. That floor is not the
estimator — it is the intron's own answer: the factory self-solves an RNA-free intron to **f_g = 0.9821
against an oracle of 1.0000**. That 1.8 % phantom RNA, reframed to the seam, is the entire 0.294.

And the 1.8 % is **ψ's reference, by construction**: `_JEFFREYS_REF` makes the composition Beta(½,½), which
is proper precisely because it **forbids the simplex vertices** — `reference_prior_derivation.md` §10.5
already records this as the known cost, "it forbids the simplex vertices, where some truth genuinely lives".
A library with no nascent RNA has introns *at* the vertex.

So the chain is: ψ's proper-prior choice → intron `f_g = 0.982` not 1.000 → share 0.294 not 0 → RNA
over-claim at the seam → intron contamination. **The two-step cannot beat the conservative peel until a node
whose truth is at the vertex can reach it.** That is the next thing to solve, and it is a ψ-reference
question, not a message question.

## 11. ⭐ THE ψ VERTEX, interrogated — the bias is REAL, the obvious fix is UNAVAILABLE, and why

§10.3 traced the two-step's blocker to one number: an RNA-free intron solves to `f_g = 0.9821` against an
oracle of 1.0000. Interrogated (`scratchpad/v1_vertex.py`), four findings.

### 11.1 It IS the reference, and only at the vertex

Reconstructing the intron λ-posterior offline from the captured factory factor:

| ψ variant | median | mean f_g | mode | \|med − oracle\| |
|---|---|---|---|---|
| factory + reference (**shipped**) | 0.9744 | 0.9672 | 0.9552 | **0.0256** |
| **factory ONLY** | 0.9953 | 0.9858 | 0.9777 | **0.0047** |
| reference only | 0.5423 | 0.5000 | 0.5423 | 0.4577 |

**5.4× better without the reference at the vertex, and ~no difference away from it** (`nrna_present`, truth
0.8478: 0.0297 → 0.0293). The mechanism is exact: the reference's `e^{−λ/2}` tail shifts the median by
`2·ln 2 = 1.386` nats whenever the likelihood is flat above the mode — measured shift λ 5.36 → 3.94 = 1.42.
The READOUT is not at fault: the median (0.9744) beats both the mean (0.9672) and the mode (0.9552), so
`_posterior_median_fg` is already the best of the three. Nor is the grid: σ(L) = 0.999955, nowhere near
binding at 0.982.

### 11.2 ⛔ But dropping the reference is NOT available — the factory-only ψ is improper

| L | 6 | 10 | 14 | 20 | 30 |
|---|---|---|---|---|---|
| factory-only median | 0.9835 | 0.9954 | 0.9983 | 0.9993 | 0.9995 |

The median drifts monotonically with the window, i.e. **L silently becomes the prior strength** — exactly what
`_JEFFREYS_REF`'s docstring warns about ("an improper ψ has plateau mass growing linearly in L"). So the
reference is load-bearing and must stay. `_rna_arm`'s own docstring already said it: "the reference is what
bounds the `f_g → 1` vertex today, and it is the ONLY thing doing so."

### 11.3 The deeper reading — the likelihood SATURATES, so 0.975 is honest

That L-drift is not a numerical artifact, it is the diagnosis: **the data carries no information above
λ ≈ 4**, because `f_g = 0.98` and `f_g = 0.9999` both predict "essentially all gDNA". The factory's NB factor
is flat there. So the posterior tail is prior-determined by construction, `f_g = 0.975` is the honest median
of a proper posterior, and the "phantom" 2.5 % RNA is really the statement *"I cannot rule out 2.5 % RNA"*.

**This reframes the problem.** There is no better point estimate to find — the defect is that a SATURATED
quantity is being propagated as though it were measured. `reference_prior_derivation.md` §10.5 already
licensed the cost ("it forbids the simplex vertices, where some truth genuinely lives"); what was missing is
that pass-0 then *uses* the forbidden-vertex answer as if it were a measurement.

### 11.4 A constant-free SATURATION DETECTOR exists and works

Saturation is directly observable: **is the factory log-factor still rising at the top of the λ grid?** — a
comparison of the last two grid cells, no threshold.

| condition | oracle f_g | % saturated | f_g \| SAT | f_g \| not |
|---|---|---|---|---|
| gdna300 ss0.50 nrna_none capON | 1.0000 | 45 % | 1.0000 | 1.0000 |
| gdna300 ss0.99 nrna_none capON | 1.0000 | 48 % | 1.0000 | 1.0000 |
| gdna300 ss0.99 present capON | 0.8478 | 14 % | **0.9874** | 0.8246 |
| gdna100 ss0.50 present capON | 0.6520 | 9 % | **0.9690** | 0.6216 |
| gdna300 ss0.99 present capOFF | 0.6843 | 4 % | **0.9932** | 0.6700 |

It fires on 45–48 % of introns where the truth is at the vertex and 4–14 % where it is not — and *within* a
nascent-RNA library it picks out precisely the introns that are individually RNA-free (0.969–0.993 vs
0.62–0.82). That is a usable, constant-free flag.

⚠ One correction it forced: `w_true` at seams beside saturated introns reads 0.31–0.65, not 0 — but that
average is contaminated by boundaries with NO spliced mass, where the share is undefined and the probe
defaults it to 1. Do not read those numbers as "RNA continues past an RNA-free intron"; re-measure on
junction-bearing seams only before using them.

### 11.5 What this says to do

1. **Do not chase a better point estimate at the vertex** — §11.1–11.3 close that off; the reference is
   correct and required, and the median is already the best readout.
2. **Use the detector for PRECISION, not for the mode.** Where the factory is saturated, the node's RNA level
   is *unidentified*, so its RNA claim should carry ~zero precision rather than the precision of a median it
   cannot support. That is the governing principle applied exactly, and it is what pass-0 owes Phase 2 — an
   honest "I don't know", not a confident 2.5 %.
3. **The remaining headroom is the UNSATURATED regime.** Where the intron is genuinely RNA-bearing the
   two-step's share estimator was decent (|Δ| 0.05–0.19, §10.1) while HEAD's conservative `w = 0` is badly
   wrong (0.17–0.41). A regime-aware share — intron estimate when unsaturated, conservative when saturated —
   is the concrete next experiment, and the detector makes it constant-free.

## 12. Tools

`scripts/debug/pass0_error_table.py` — the suite state of play in READS with the trust (`errQ1conf`) view.
`scratchpad/t1_char.py` (characterize + ψ ablation, bit-exact), `t2_strata.py` (channel-arrival strata +
distance-to-anchor), `t3_streams.py` (who seeds which stream), `t4_gdna_transport.py` (per-edge transport
decomposition + candidate frame ratios), `t5_lr_gap.py` (the refuted left/right statistic),
`t6_ceiling.py` (**the oracle-mode ceiling — use this before chasing any mode work**), `t7_pin_check.py`
(does the pin cancel `r`; the imputation-premise floor).
