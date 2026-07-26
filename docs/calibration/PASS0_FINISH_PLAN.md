# Finishing pass-0 — the prioritized work list

**Owner-agreed 2026-07-26. Live status in the table below; the run instructions are in
`SESSION_2026_07_26_HANDOFF_13.md` §6.** Read `ROADMAP.md` first, then
`SESSION_2026_07_26_HANDOFF_12.md` (the live handoff), then this.

This document exists so the set does not get lost between sessions. **Update the status column in place**;
do not start a competing list.

| # | item | status |
|---|---|---|
| **P0** | change the debug loop's loss to **confidently-wrong mass** | ✅ **DONE** 2026-07-26 |
| ~~P2~~ | ~~introns at 91 % `errQ1conf`~~ | ⛔ **REFUTED** — a selection artifact, §P2 |
| **P1** | root-cause the **capture-OFF** confident error | 🔬 **ROOT-CAUSED**, §P1b |
| **P1c** | the WEIGHTED RESCALE | ⤳ **SUPERSEDED** by P1f, which makes it inert (0/0/32); law kept, path removed |
| **P1f** | **the MESSAGE PACKET** — λ and θ carried explicitly, fused by their OWN precisions | ✅ **LANDED** — `weighted_rescale_design.md` §9 |
| **P2a** | ⭐ **RE-TEST the Phase-2 λ-curvature step at HEAD** (~6 lines, `HANDOFF_6` §3) | ▶ **THE NEXT ACTION.** Both reasons it was shelved have evaporated: the refit no longer regresses unstranded-capON (0.1523 → **0.1275** at HEAD, 25/7 overall), and the hyperprior's substrate is at `z2` = **2.26**. If it holds, **Phase 2 unblocks** |
| **P1e** | the conservation **SURPRISE** — never inert | ✅ **LANDED** 2026-07-26, common-direction, **scoped to `δ < 0`**, default ON (`RIGEL_P1E_OFF=1` ablates). suite 0.0870 → **0.0841**, fit-substrate 0.0883 → **0.0845**, held-fixed z2 ALL 8.53 → **1.74**, exon-single → **0.88**. ⚠ **partly a DEBT — prices a bias as a variance**, `variance_ledger.md` §6. Historical: P1d's residual regression points straight at it: the premium must go to the arm whose premise FAILED, and `claim/obs` is the observable that says which |
| — | **`variance_ledger.md`** — every variance term, what it prices, and the no-overlap proof | ✅ **the standing audit; update it whenever a term changes** |
| **P1d** | the graft's **PREMISE** variance (was: "extrapolation" — refuted) | ✅ **LANDED** 2026-07-26, `graft_premise_logvar` — ONE library-level MoM scalar, default ON. Confidently-wrong **−36.6 %**, z2 ALL 15.55 → **8.98**, exon-single 8.60 → **2.46**, for +1.2 % mwae. (A per-edge two-seam variant landed first and was reverted same-day — χ²₁ noise + a BP violation.) The mechanism is **transcript TERMINI**, not extrapolation — `variance_ledger.md` §3.0 |
| **P3** | AMBIG exon over-confidence — `z2\|Q1` 183 → 92.1 → **64.4** | ⤵ **DEMOTED.** Still the largest single defect (259,493 CWRONG, 34 %) but **excluded from the hyperprior fit by construction**, so it blocks the RE-SOLVE's quality, not Phase 2's start. Largely downstream of the trained prior |
| **P4** | the `_far` across-the-seam level estimator | ✅ **DELETED** 2026-07-26. A 100 % structural data reuse — its node is *identically* the other message's source (35,421 peel edges). On a held-fixed node set: boundary-single `z2` **5.58 → 2.98**, ALL 8.98 → **8.28**. Price +0.0007 (r0) / +0.0009 (r1) |
| **P4b** | the last BP violation — `_RHO_ITERS = 2` | ✅ **FIXED** 2026-07-26 by `_RHO_ITERS = 1`. Price −0.0002 (r0, i.e. BETTER) / +0.0011 (r1, concentrated in unstranded × capON +0.0104). **The solver is now BP-legal.** The elegant pairwise-potential refactor is NOT justified — it would bid a message-representation rewrite against ~0.001. Historical: The 2nd iteration's reframe faces are built from the destination's **fused posterior**, i.e. from the other message. Measured median \|Δlog ρ_face\| between iterations 0.0116 (stranded capOFF) → 0.1242 (unstranded capON), >1 % on 52.7–79.0 % of nodes. Deleting `_far` removed the largest violation, not the last one |
| P5 | re-test `r` from the gDNA channel (the "moot" verdict is stale) | queued |
| P6 | fragment length as a FOURTH information source — **gated** | scoping only |
| **P1g** | ⭐ **put TSS/TES in the region/boundary map — THE STRUCTURAL DEBT BEHIND `ω_graft`** (owner: high priority, forthcoming). `ROADMAP.md` §6 deferred it as "expected to be low-impact"; it is measured at **62–72 % of the graft's premise error** and a free annotation bit splits ω̂ **≥30×**. P1d prices it blind through the two-seam gap; naming it directly would price it exactly | ▶ **UN-DEFERRED by measurement** |
| **P-hold** | hold out / reset the nodes not used to fit the hyperprior | ⛔ **MEASURED UNNECESSARY.** The refit *already resets* (`calibrate.py:420`), and holding the excluded half out moves the fit substrate by **±0.0006** — there is no contamination worth removing. A stronger version of the same idea (skip unidentified nodes, defer to the prior) was already derived, implemented and **empirically refuted** in `solve_gate_design.md`: +0.010 (r0) / +0.025 (r1) |
| **P-solv** | ⭐ **a clean `solvable` predicate**, and expose it to the prior fit | ▶ **NEXT.** Owner: *"a solvable flag is helpful for now"* — expose it; the EXCLUSION question (which nodes to drop from the fit) is deliberately deferred, because dropping a node that reports a real enriched mode could hurt badly. Evidence: `τ_own > 0` separates hard: 29.6 % of nodes at mwae **0.0227** (7.7 % of error mass) vs 70.4 % at **0.1143** (92.3 %). A composite prior-free predicate splits 45.6 % @ 0.0147 vs 54.4 % @ 0.1480 — **10.4× on unstranded**. ⚠ The existing `cap["solvable"]` flag agrees with it only **45.7 %** of the time, so it is not that predicate |
| — | **THEN** the re-solve + a ship candidate | the hyperprior fit itself is no longer blocked — see P2a |

---

> ### ⚠⚠ P1e IS PARTLY A DEBT — IT PRICES A BIAS AS A VARIANCE. DO NOT LET THIS BE FORGOTTEN
>
> P1e damps a message by the unexplained part of its conservation violation `δ = log(M/S)`. **On a large
> share of its firing mass `δ` is a systematic BIAS, not scatter** — measured `E[δ]` ≈ −0.5 to −1.5, with a
> bias share of **53–77 %** on graft × one-component messages and **98.9–99.2 %** at intergenic
> destinations. A variance prices random scatter; pricing a bias as a variance does not move the mode toward
> truth, it only weakens the message — **and it never shrinks**, which breaks the ledger §2.2 guarantee that
> a DerSimonian–Laird residual vanishes as the model improves.
>
> **It is landed anyway** because it is the only change measured to improve ACCURACY and honest PRECISION at
> the same time (suite 0.0870 → 0.0841, fit-substrate 0.0883 → **0.0845**, held-fixed `z2` ALL 8.53 → **1.74**,
> exon-single 3.57 → **0.88**), and because pass-0 is explicitly allowed to be weak-and-correctable. That is
> a pragmatic trade, not a derivation.
>
> **Four things that must not be overlooked:**
> 1. ⛔ **It was hoped the bias-dominated strata were inert** (intergenic nodes are `solvable = False`).
>    **Measured and REFUTED: 90–100 % of the damping mass lands on solvable destinations.** The bias is live.
> 2. **The magnitude is not what works.** Control arms: a flat pooled constant on the same firing set beats
>    the derived `b̂²` on 3 of 4 conditions, and `b̂² := δ²` with no null subtraction is identical. Permuting
>    `b̂²` within edge class FAILS and damping all grafts uniformly FAILS — so **`δ` identifies WHICH message
>    to distrust, and the calibrated magnitude adds nothing.** Exactly ω_graft's shape.
> 3. **`δ` is not the hard observable it was scoped as.** It is **M-free at every region node** (verified to
>    1.9e-16 — the message is built ∝ `M`, so `M` cancels). On a plain edge it is a composition disagreement
>    with the destination's own belief. The genuine content is in the ROUTING: **84.2 % of Σδ² is on peel and
>    graft edges.**
> 4. **The right fix for the bias half is a MODE fix, not a variance.** When the bias strata are diagnosed,
>    P1e must SHRINK — if it does not, the model has not improved.
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

## 0a. ⭐ THE OBJECTIVE (owner, 2026-07-26) — ACCURACY **and** honest precision, together

The Phase-2 re-solve **already resets** (`calibrate.py:419-420`: `belief = _init_belief()` before
`_sweep(gdna_hyperprior)`), so nothing from pass-0 is inherited but the fitted prior. That does **not**
demote precision — the owner's ruling is explicit:

> *"Honest precision is critical. We may adopt an iterative strategy. We may decide to use a
> precision-weighted gDNA hyperprior fit. If we push towards a combination of ACCURACY with honest
> precision, we will prevail."*

Three live consumers of honest precision: pass-0's own fused **mode** (hence the substrate's accuracy), the
**precision-weighted** hyperprior fit (`var_gdna` weights the NPMLE), and a possible iterative re-solve.

**The number that was missing:** ACCURACY restricted to the hyperprior's fit substrate
(`REGION & (single | gonly) & live` — exon-single, intron-single, and the exact intergenic anchors).
Suite-wide mwae is a proxy for it, not a substitute. **Track both from now on**, plus `z2`.

**Also standing (owner):** the code must reach Phase 2 *efficient, clear, concise, readable, maintainable —
not over-engineered, as simple as possible.* Simplification is a gate, not a nicety.

## 0b. ⭐ STATE OF PLAY (2026-07-26, end of session)

Suite **0.0855 (refit=0) / 0.0671 (refit=1)** against the session's starting HEAD of 0.0885 / 0.0678 —
**20 better / 5 worse / 7 flat** and **17 / 7 / 8**. And on the axis that decides whether pass-0 is finished:

| | session start | now |
|---|---|---|
| suite error | 12,344,845 | **11,928,101** |
| **confidently-wrong** | 1,777,658 (14.4 %) | **1,186,552 (9.9 %)** |
| exon AMBIG confidently-wrong | 693,557 (20.1 %) | **378,285 (11.8 %)** |

**But the mass is not the whole story, and this is the number that says pass-0 is NOT finished:** `z2` inside
the confident quartile — *is the confidence earned?* — is still **8.4 on single-strand exons, 92.1 on AMBIG
exons, 15.5 overall.** Introns (1.5) are the only honest class. A node at `z2 = 9` will be weighted nine times
too heavily by a hyperprior fit that trusts declared precision, which is precisely the poison mechanism.

**So: the error MASS is close to its prior-free ceiling; the CALIBRATION of that error is not.** Everything
still on the list below is chosen for its effect on `z2`, not on mwae. P1d alone is measured at
**exon-single `z2` 8.4 → 1.2** (`variance_ledger.md` §4).

## 0. Why the list is ordered this way — the goal is a SUBSTRATE, not a score

**On the error-MASS axis, pass-0 is already near its prior-free ceiling** on the arms that dominate.
`SESSION_2026_07_25_HANDOFF_10.md` §3 measured 67–75 % of the error in the largest scenarios as reachable
"by a better MODE only", premise-limited (neighbouring nodes genuinely differ in composition by
|Δf_g| = 0.24), with **×10 on every precision moving nothing**. Chasing mwae further has a low ceiling.

What is *not* done is the axis the hyperprior actually consumes. Pass-0's deliverable is a set of nodes the
hyperprior fit can trust — and it can down-weight a node that says it is unsure. It cannot do anything about
a node that is **confidently wrong**. Suite state at the head of this list:

> **1,516,367 of 12,490,634 error reads (12.1 %) sit on the most-confident quartile.**

Ranking scenarios by `ERR × errQ1conf` instead of by `ERR` produces a **completely different list**, and that
is the whole point of P0:

| scenario | ERR | errQ1conf | confidently-wrong | share |
|---|---|---|---|---|
| `gdna300 ss0.50 present capOFF` | 527,005 | **41.8 %** | 220,288 | 14.5 % |
| `gdna300 ss0.50 none capOFF` | 480,210 | **36.2 %** | 173,836 | 11.5 % |
| `none ss0.50 none capOFF` | 595,139 | 24.1 % | 143,428 | 9.5 % |
| `gdna300 ss0.99 present capON` | 627,063 | 16.2 % | 101,584 | 6.7 % |
| `gdna300 ss0.99 present capOFF` | 126,228 | **75.8 %** | 95,681 | 6.3 % |

The top three are all **capture-OFF**, and **none of them is in the error-mass top five** (which is entirely
unstranded × capture-ON — the hyperprior's problem, not pass-0's). Conversely `verystrong` puts 76–88 % of
its error on the *least*-confident quartile: badly wrong and honest about it, which is acceptable.

By node class:

| class | ERR | errQ1conf | confidently-wrong | share of all |
|---|---|---|---|---|
| exon AMBIG | 3,498,633 | 18.0 % | **629,754** | 41.5 % |
| exon single | 7,676,603 | 7.9 % | 606,452 | 40.0 % |
| intron (both) | 254,708 | **~91 %** | 231,716 | 15.2 % |
| boundary (both) | 1,060,690 | ~4.5 % | 46,703 | 3.1 % |

---

## P0 — change the loop's loss function

**Why.** The standing methodology (memory `pass0_debug_iteration_loop`) is *suite → worst scenario by error
MASS → worst nodes → root cause → fix*. That loss now points at the premise-limited arms every single time.
Swap it for **confidently-wrong mass** and the loop starts pointing at things pass-0 can actually fix.

**What.** `scripts/debug/pass0_error_table.py` already computes both columns; add the product, sort by it,
and make it the default view. Free.

**⚠ Do not read `errQ1conf` naively.** It is the share of a class's error sitting on the suite-wide
most-confident quartile. If a whole class is confident for a *legitimate* reason, its `errQ1conf` is high by
selection, not by over-confidence. The metric that separates the two is **`z2 = E[err²]/E[Var]`** (1.0 =
honest, > 1 = over-confident). P0 must report both, or P2 will chase a selection artifact.

## ✅ P0 — DONE (2026-07-26)

`pass0_error_table.py` now orders the scenario table by **confidently-wrong mass** and prints two new columns:
`CWRONG` (the absolute reads) and **`z2`** (the calibration, `E[(f_g−oracle)²]/E[Var(f_g)]`), plus `%nodeQ1`
beside `errQ1conf` so the selection confound is visible. It settled P2 on its first run.

### ⭐ The measurement that re-ordered this list

`z2` restricted to the CONFIDENT quartile — *is the confidence EARNED?* — cross-tabbed by regime and class.
**Confidently-wrong reads: 1,514,652 total.**

| regime | exon single | exon AMBIG | boundary | intron |
|---|---|---|---|---|
| stranded × capOFF | 8,565 [z2 0] | 125,029 [**280**] | 9,852 [7] | 113,973 [3] |
| stranded × capON | 121,001 [1] | 143,908 [**107**] | 32,441 [9] | 17,931 [3] |
| **unstranded × capOFF** | **378,815 [33]** | 133,405 [**343**] | 2,827 [9] | 87,504 [1] |
| unstranded × capON | 97,358 [6] | 140,988 [**96**] | 1,490 [8] | 12,332 [1] |
| unstranded × verystrong | 5 [1] | 87,228 [**717**] | — | — |
| **TOTAL** | **605,744 [10]** | **630,558 [183]** | 46,610 [8] | **231,740 [2]** |

Three conclusions, and they set the rest of the list:

1. **Introns are the BEST-calibrated class in the solver** — `z2|Q1` = 1–3 in every regime. P2 is dead.
2. **AMBIG exons are over-confident by two orders of magnitude in EVERY regime** (96–717), evenly spread
   (87 k–144 k per regime). This is **structural, not regime-specific** — the largest and most uniform defect
   in the solver, and it is 41.6 % of all confidently-wrong mass. → **P3, promoted.**
3. **Single-strand exon over-confidence is sharply LOCALIZED**: `z2|Q1` = 33 on unstranded × capOFF
   (378,815 reads, 25 % of all confidently-wrong) against 0–6 everywhere else. A different population from
   P3, and the most tractable regime on the list. → **P1, next.**

## ⛔ P2 — REFUTED (2026-07-26). Introns are not over-confident; the metric was.

The 91 % `errQ1conf` is a **selection artifact**, exactly as the P0 warning anticipated:

* **48.4 %** of intron NODES sit in the suite-wide confident quartile (against 12–29 % for every other
  class) — because the intron factory *legitimately* makes them confident;
* their calibration is **`z2` = 0.011 overall and 1.73 within the confident quartile** — honest to
  conservative, the best in the solver. Median |err| 0.0218 against a median declared sd of 0.0191.

This is consistent with `SESSION_2026_07_25_HANDOFF_10.md` §12.3 (`z2` = 0.25–1.22 at introns), which was
already on record; the `errQ1conf` framing simply obscured it. **There is no trust defect at introns and no
fix to make here.** 231,740 reads of "confidently wrong" intron error are a legitimately confident class
being slightly off.

⚠ **What survives:** the `intron_prior`-array vs `tau_lam` two-path inconsistency
(`SESSION_2026_07_25_HANDOFF_10.md` §6) is still a genuine code-level inconsistency — the same quantity
reaching ψ and `v_own` by two routes that disagree — and worth closing on correctness grounds. It is **not** a
trust defect and does **not** justify its own investigation. Fold it into whatever next touches that code.

## P2 (historical) — what the item said before it was refuted

**Why first among the fixes.** Introns are only 2.1 % of suite error mass, but they are **where the
hyperprior's gDNA density is measured**. This is the direct poison path into Phase 2, and it is the one item
with a *named, already-diagnosed* candidate cause rather than an open investigation.

**The candidate cause** (`SESSION_2026_07_25_HANDOFF_10.md` §6, recorded and never fixed):

> ψ weights the `intron_prior` ARRAY while `v_own` reads `tau_lam` — the same quantity down two paths that
> then disagree. Proved by a ×100 saturation test: scaling `tau_fac` alone saturates (intron error
> 6,872 → 5,753 at ×100) because the two paths do not move together.

**⚠ The premise must be checked first.** `SESSION_2026_07_25_HANDOFF_10.md` §12.3 measured intron
`z2 = 0.25–1.22` — i.e. *honest to conservative*. If ~91 % of intron NODES already sit in the confident
quartile (because the factory legitimately makes them confident), then `errQ1conf ≈ 91 %` is a selection
effect and there is nothing to fix. **Measure quartile membership and `z2` per class before changing
anything.**

## 🔬 P1 — ROOT-CAUSED (2026-07-26). The RNA measurement channel wears a borrowed mode.

**The population** (`scratchpad/p1_capoff.py`, paired against the same nodes at capture-ON — same genome,
same partition, so the node index is the same DNA). Single-strand exons, unstranded, capture-OFF, in the
suite-wide confident quartile: **216 nodes, 513,966 reads, z2 = 36.6**.

```
  stratum                                n       reads   |err|      sd      z2    c_tau     cm_g     cm_R
  capture OFF  all exon-single         487   1,616,671  0.2279  0.1115     5.7    10.36    63.92   929.09
  capture OFF  confident quartile      216     513,966  0.2815  0.0487    36.6     6.87   109.20   651.07
  capture ON   same nodes as above     216   1,663,473  0.0879  0.1396     0.5     4.10    72.61    61.36
```

**Where the posterior precision comes from, on exactly those nodes:**

| | `c_tau` (λ) | `cm_g` (gDNA meas.) | **`cm_R` (RNA meas.)** | `tau_own` |
|---|---|---|---|---|
| capture OFF | 6.87 (1.9 %) | 109.20 (23.2 %) | **651.07 (74.9 %)** | 0.000 |
| capture ON | 4.10 (14.9 %) | 72.61 (61.6 %) | 61.36 (23.5 %) | 0.000 |

The individual nodes are stark — oracle `f_g` = 0.016–0.069 (i.e. **93–98 % RNA**), self-solve 0.49–0.51
(ψ's uninformed reference, as expected with `τ_own = 0`), **solved 0.23–0.59** — barely moved off the
reference — at `sd` = 0.044–0.056, so **z = 4.5–10.1**. Enormous precision, almost no movement.

**Why it is capture-OFF specific.** Three damping terms exist and capture-OFF removes the last one:
* `τ_own = 0` on unstranded data ⇒ `v_own = ∞` ⇒ the DL mismatch term is **inert by design** (M7);
* `r ≈ 1` off-capture ⇒ M5's `σ²_transfer ≈ 0`;
* and **M8's `graft_frame_logvar = (log r)²` is identically 0 at `r = 1` by construction** — it is the ONLY
  thing damping the grafted spliced count `_spc = SP/(1 + SP·σ²)`, so off-capture the spliced count enters at
  its **raw Poisson precision** (measured `cm_R` up to 2,025 on a single node).

**⛔ Ablating the channel is NOT the fix** (arm `abl_r0`, `RIGEL_RNAMEAS_OFF=1`): **0.0895 → 0.1033, 4 better
/ 17 worse**, and catastrophic exactly where the channel is load-bearing — `gdna_none` 0.1063 → **0.1438**,
`none_ss0.50_nrna_none_capOFF` 0.3624 → **0.5049**. It is the only thing that lets a zero-gDNA library say
"my mass is all RNA". It did help the arms M11 regressed (stranded 0.0347 → 0.0331, 3 better / 1 worse), which
is consistent: the channel is *essential* where RNA is everything and *harmful* where gDNA is high.

### ⛔ CORRECTION (2026-07-26, later the same day) — the paragraph below is WRONG

The "borrowed mode" reading was inferred from precision bookkeeping, not measured. Dumping one node in
fragments and base pairs (`scratchpad/p1_worked.py`) refutes it: **the spliced-RNA claim is essentially
EXACT.** See §P1b. The passage is kept because the precision bookkeeping it describes is real and still
explains the *confidence*; it simply does not explain the *error*.

### The defect as first (mis)diagnosed — the mode/precision mismatch

`bp_solver` builds the RNA measurement ψ factor as

```
    rna_imp = (mo_p, mo_n), (cm_p, cm_n)      mo_p = log(cp·E_r/M)     cp = _fuse_v(ap, app, bp, bpp)
```

— the **mode comes from the FUSED density** (both directions, weighted by the full mode-fusion precisions)
while the **precision comes from the MEASUREMENT stream alone** (`mg`/`mp`/`mn`, which are *precisions with no
companion density*: the relay carries `mg[i] = mg_own[i] + tmg` and nothing else). So a spliced COUNT's
confidence is being attached to a mode that count does not determine. Where the fused mode happens to be right
(RNA-rich libraries) the high precision helps; where it is wrong (gDNA-rich) it pins the node hard and wrong.
That is exactly the observed split, and it is why "huge precision, almost no movement" can co-occur.

**The fix is a design call, not a patch:** the three-stream separation tracks three precisions but only one
density. Either give the measurement stream its own density (so its ψ factor is its own claim at its own
precision), or price the premise gap between the two — a spliced count measures the JUNCTION's flux exactly
but the destination exon's RNA density only under the imputation premise "my RNA all flows through this
junction", and that premise error is currently unpriced. **Do not start building before deciding which.**

## P1 (historical) — the item as written before the investigation

Flagged in `SESSION_2026_07_25_HANDOFF_10.md` §5 as "the more important Phase-2 risk" and **never
investigated**. It is also the cleanest regime available: **no capture cliff**, so `r ≈ 1` and the reframe is
near-exact — which means whatever is found is a composition/premise or precision defect, not a frame artifact.
Best signal-to-noise on the list. Start at `gdna300_ss_0.50_nrna_present_capture_off` (220 k confidently-wrong
reads, the largest single block in the suite).

## P3 — AMBIG exon over-confidence (630 k reads, 41.5 % of all confidently-wrong)

At `errQ1conf` 18.0 % they are 2.3× less honest than single-strand exons (7.9 %). `τ_own = 0` on every AMBIG
node **by construction**, so the DL mismatch term is inert (`v_own = ∞` ⇒ `b̂² = 0`) and nothing damps their
messages. The *mode* is Phase-2's — established in `SESSION_2026_07_25_HANDOFF_9.md` §0b, do not re-litigate
it — but the *confidence* is pass-0's: **a node with zero own composition evidence should not be landing in
the most-confident quartile.**

## P4 — FAR → a proper message, and M11's pre-DL precision

Two related correctness defects in what landed 2026-07-26:

* **The FAR level estimator is a LOOKUP, not a message.** It reads `op[far]` (the far node's message-free
  self-solve) and reframes it by the pre-pass `ρ_tot` ratio. It is charged the M5 hop variance but **never
  sees M7's `b̂²`**, and — worse — it breaks BP's cardinal rule: the message *into* `i` *from the left* now
  depends, through `w`, on the RIGHT neighbour's belief. The two messages then share information and the
  fused belief is over-confident. This project has already been bitten by exactly that failure mode
  (`SESSION_2026_07_25_HANDOFF_10.md` §4, the left/right gap as a DL second study: "both sides share the same
  `r`/graft bias, so they agree while both are wrong").
  **The fix is BP's own rule:** a message α may depend on everything *except* α's own source. So the LEFT
  message's share should use the **RIGHT message's** delivered RNA claim as its far level — which the combine
  already computes. Restructure `_transport` into two stages (form both raw claims → form each direction's
  share from the *other's* claim → apply → pin). The relay only ever has one direction, so it uses own ⊕ mass
  there; the combine, where the answer is formed, uses all three properly.
  ⚠ It is **load-bearing** (ablating it costs 0.0895 → 0.0905 / 0.0693 → 0.0703) — replace, do not delete.
* **M11 consumes `tpg` BEFORE the DL deflation** (`mismatch_deflate` runs at the end of `_transport`), so the
  gDNA claim it reads is missing `b̂²`. On unstranded data DL is inert anyway (`v_own = ∞`), but **on
  stranded data it is live — and stranded × capON is exactly where the M11 landing regressed.**

## P5 — `r` from the gDNA channel

`r_M/E_g` scored **0.808 → 0.330** against the true capture step and was shelved in
`SESSION_2026_07_25_HANDOFF_10.md` §4 as "**moot**, because `_pin_v` cancels `r` on 87.6 % of the error".
**That verdict is stale**: M11 is the first consumer in the solver for which `r` does *not* cancel. Re-test
scoped to M11's input only. The estimator already exists; cheap.

## P6 — fragment length as a FOURTH information source (GATED — scope before building)

Measured 2026-07-26 (`scratchpad/e1_enrichment.py`, recorded in HANDOFF_12 §5b): the gDNA and RNA
fragment-length distributions differ enough that a fragment's LENGTH is informative about which component it
came from. Per-fragment Fisher information about `f_g`: **`I_FL` = 1.19–3.98**, against the strand's 1.26 at
κ = 0.99 and **exactly 0** at κ = ½. In `τ_λ` terms on unstranded capture-ON: exons `τ_own = 0` → `τ_FL` ≈ 29,
boundaries `0` → 2.5, against the intron factory's anchoring `τ` of 0.254.

**Do not build until these two are answered:**
1. **Circularity.** `fl.py`'s gDNA/RNA pmfs must be estimated from anchors that do not depend on
   calibration's own answer — intergenic regions are pure gDNA and spliced fragments are pure RNA, so both
   marginals are directly observable, but this must be verified, not assumed.
2. **Realism.** TV distance in this suite is 0.46–0.997, i.e. the lengths nearly *label* each fragment. Real
   libraries size-select gDNA and RNA together and overlap far more, so the suite will over-sell this badly —
   the same family of artifact as `synthetic_suite_is_poisson_omega_zero`. Get a realistic separability
   estimate first.

Notes for when it is built: it does **not** violate count-zero-information (a fragment's *length* is a
per-fragment observable like its *strand*, not a count vote), and the minimal implementation is **one extra
float per node** in the C++ accumulator — the per-node mean fragment length is a sufficient statistic for a
two-component method-of-moments split, `E[ℓ] = f_g·μ_g + (1−f_g)·μ_r`. The accumulator's existing FL pools are
library-wide (3 region-types × 2 compartments), so they cannot serve this as they stand.


## ⭐ P1b — THE ACTUAL DEFECT: the conservation rescale punishes a measurement to accommodate a guess

One node, in fragments and base pairs (`scratchpad/p1_worked.py`, node 2651,
`gdna300_ss_0.50_nrna_present_capture_off`, 78,742 fragments — the largest confidently-wrong exon):

```
  the exon:      78,742 unspliced fragments over ~2,100 bp   ->  37.3 fragments/bp if it were 100 % RNA
  TRUTH:         98.4 % RNA  (gDNA density 0.58/bp,  RNA density 36.7/bp)
  the junction next door:  7,069 spliced fragments over 100 bp  ->  35.3 fragments/bp

  what the two channels claim, in the exon's own frame:
                        claimed      TRUTH    ratio
      gDNA density      27.6319     0.5827    47.4x        <-- an IMPUTATION, 47x too big
      RNA  density      36.6453    36.6734     1.0x        <-- a MEASUREMENT, essentially exact

  together they claim 138,000 fragments.  The exon observed 78,742.  ->  1.75x
```

**The RNA claim is right. The gDNA claim is 47× wrong. And the conservation rescale (`_pin_v`) divides BOTH
by 1.75** — so it takes 43 % off a direct measurement in order to make room for a guess. The result is
`f_gDNA = 0.322` against a truth of 0.016, stated with a std-dev of 0.049: **6.2 σ wrong**.

Had the whole correction gone to the gDNA claim (the uncertain one), the residual would be
78,742 − 77,395 = 1,347 gDNA fragments ⇒ **`f_gDNA` = 0.017 against a truth of 0.016.**

### Neither extreme is right — the correction must be APPORTIONED

Exon-single nodes, mass-weighted |f_g − oracle|, against two counterfactual rescales:

| condition | claim/obs | SHIPPED | RNA-claim exact | gDNA-claim exact |
|---|---|---|---|---|
| `gdna300 ss0.50 present capOFF` | **1.31×** | 0.2179 | **0.0647** | 0.3690 |
| `gdna300 ss0.50 none capOFF` | **1.33×** | 0.2096 | **0.0469** | 0.3713 |
| `gdna100 ss0.50 present capOFF` | **1.26×** | 0.1651 | **0.0434** | 0.2950 |
| `gdna1 ss0.50 present capON` | 1.03× | 0.0322 | **0.0047** | 0.0301 |
| `gdna300 ss0.99 present capOFF` | **1.01×** | **0.0056** | 0.0419 | 0.0440 |
| `none ss0.50 present capOFF` | 1.00× | **0.0007** | 0.0002 | 0.0001 |
| `gdna300 ss0.50 present capON` | 1.15× | **0.1524** | 0.2012 | 0.1814 |
| `gdna300 ss0.99 present capON` | 1.19× | **0.0477** | 0.2176 | 0.1203 |
| `gdna300 ss0.50 none capON` | 1.20× | **0.1834** | 0.3209 | 0.1893 |
| `gdna100 verystrong` | 1.04× | **0.1953** | 0.2238 | 0.2150 |

Two things fall out, and the second is the useful one:

1. **Neither extreme wins everywhere.** RNA-exact is 3.4–6.9× better on unstranded capture-OFF and on
   low-gDNA, and much worse on gDNA-rich capture-ON. So the fix is not "trust the measurement" — it is
   **apportion the conservation correction by how well each component is known**, which is the weighted pin
   derived earlier this session (`SESSION_2026_07_26_HANDOFF_12.md` working notes): with a common frame error
   `σ²_r` and independent composition errors, the correction splits as
   `a = δ(σ²_r + α·v_g)/(σ²_r + α²v_g + β²v_R)`, `b = δ(σ²_r + β·v_R)/(…)`. Today's blind common rescale is
   the `v → 0` limit of that formula (correct only when the *only* error is a common frame error), and
   "RNA-exact" is the `v_R → 0` limit. The general form is what is needed.
2. **`claim/obs` PREDICTS the error, with no oracle.** 1.00–1.01× ⇒ the shipped answer is excellent
   (0.0007, 0.0056). 1.26–1.33× ⇒ it is bad (0.165–0.218). This is a **directly observable, prior-free
   diagnostic of message error** — the owner's *"what is the likelihood my available counts account for this
   composition?"* — and nothing in the solver currently uses it as evidence. `_pin_v` merely erases it.

### Why the error is also CONFIDENT (the second half, and the first diagnosis survives here)

The size of the mistake is §P1b. The *confidence* attached to it is the precision bookkeeping: with
`τ_own = 0` (unstranded) the DL term is inert by design, `r ≈ 1` off-capture makes `σ²_transfer ≈ 0`, and
`graft_frame_logvar = (log r)²` is identically 0 at `r = 1` — so the grafted spliced count enters at its raw
Poisson precision (up to 2,025 on one node) as a claim about a whole exon. The owner's framing is exact: **a
spliced measurement taken over a ~100 bp junction window is being asked to speak for a ~2,100 bp exon — a 12–21×
extrapolation — and it is charged no transfer variance for that at all.**
