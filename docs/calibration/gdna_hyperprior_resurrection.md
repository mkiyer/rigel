# The gDNA hyperprior (Role B) — RESURRECTED at the rebuilt pass-0

**Status: the WIP is live again and re-measured. Date 2026-07-26, branch `calib-ambig-init-wip`, HEAD
`922abfcb`.** No production code changed — this is measurement only, on cached pass-0 beliefs.

This picks up the track that was **paused on 2026-07-21** and supersedes the archived docs' *numbers* (their
derivations stand). Read this first, then the archived pair below for the design reasoning.

---

## 0. Where the track was paused, and why

The hyperprior work stopped mid-flight at commit `243bd5ef` (2026-07-21 16:18); everything after it is pass-0
solver work. The pause is recorded verbatim in `archive/enrichment_sensitivity_worklog.md` §6 and
`archive/dna_prior_session_resume.md` §6:

> **⚑ CURRENT PRIORITY (owner's pivot): fix pass-0, not the projection.** We discovered the enriched
> under-call is a **pass-0 solver problem** on unstranded data […]. The asymmetric projection is a
> **Band-Aid** for a symptom; recovering gDNA at the source cascades to a better landscape AND sensitivity.

So the pause was **conditional**, and the condition has now been met: pass-0 was rebuilt (suite mwae ~0.15 →
**0.0841 / 0.0679**). This document answers the question the pause posed — *what did fixing pass-0 buy the
hyperprior?*

**The archived doc map** (all in `docs/calibration/archive/`, all still valid as *design*, superseded on
numbers by this file):

| doc | what it is |
|---|---|
| `dna_prior_session_resume.md` | **the resume guide** — the §0 starter prompt, the infra map, the open next steps |
| `enrichment_sensitivity_worklog.md` | **the definitive work log** — the under-call evidence (§1), the bracket concept (§4), the pass-0 dissection (§8) |
| `gdna_hyperprior_plan.md` | the authoritative *design*: projection = a **sampling likelihood**, prior = a **gentle anchor**, third line of defense |
| `gdna_crush_dissection_node1055.md` | the proof that the crush is **not** a solver bug — swap the prior's shape and the same node lands at 0.863 |
| `gdna_hyperprior_from_scratch.md` / `gdna_hyperprior_clean_slate.md` | the "gravity" principle; the KDE-vs-NPMLE strategic study |

---

## 1. Infrastructure — resurrected, with two data bugs fixed

The exploration harness (pure numpy on a cached substrate — **no `calibrate` re-runs**, ~seconds per sweep)
still imports and runs at HEAD unmodified. Two things were broken and are now fixed:

1. **The caches were gone.** `gdna_explore_lib._CACHE` pointed into a *per-session* `/private/tmp/claude-503/…`
   scratchpad, which is volatile — the 2026-07-21 caches died with that session. **Repointed to the repo's
   durable `scratchpad/`** and rebuilt at HEAD: `gdna_substrate_cache.pkl` (ambig, 32) and
   `…_quick.pkl` (quick, 16).

2. ⚠ **Two substring-matching bugs in `gdna_cache_build.py` corrupted the strata every metric is broken out
   by.** Both are fixed (`_group`, now a positional parse) and both invalidate some archived numbers:
   * **`"none" in cond` matched `nrna_none`** — so 6 real-gDNA conditions (`gdna100 … nrna_none …`) were
     bucketed as the **zero-gDNA** level. The `fabrication` specificity canary selects on exactly that
     bucket, so **it was scoring gDNA-bearing scenarios**. See §4.
   * **`capture_verystrong` fell through to "OFF"** — the *strongest* capture arm was labelled capture-off,
     which is why the archived capture-OFF rows look worse than they should.

**Run book** (`OMP_NUM_THREADS=1`, `conda activate rigel`):

```bash
python scripts/debug/gdna_cache_build.py --suite ~/Downloads/rigel_runs/ambig_dense_10mb \
       --out scratchpad/gdna_substrate_cache.pkl        # ~2 min, only after pass-0 changes
python scratchpad/gdna_undercall.py [ambig|quick]       # worklog §1 — the under-call table
python scratchpad/gdna_resurrect.py --plot              # GOAL 1 + GOAL 2 + the fit-vs-oracle figure
python scratchpad/gdna_enriched_mode.py                 # the enriched-mode census (mass + location)
```

Figure: `docs/calibration/figures/resurrect_fit_vs_oracle.png` (fit vs oracle, unstranded conditions).

---

## 2. What fixing pass-0 bought — the three re-measurements

### 2.1 The compression is FIXED; the enriched under-call is NOT

Single-strand REGION nodes, `obs − oracle` in decades (worklog §1's table, re-run):

| cap / ss | all single-strand (was) | **ENRICHED single-strand** (was) |
|---|---|---|
| OFF / 0.50 | **+0.046** (+0.415) | n/a (−0.135) |
| OFF / 0.99 | −0.017 (−0.026) | n/a (~0) |
| ON / 0.50 | **−0.076** (+0.097) | **−0.152** (−0.169) |
| ON / 0.99 | −0.035 (−0.048) | −0.028 (−0.012) |
| VSTRONG / 0.50 | +0.042 (—) | **−0.219** (—) |

**The depleted OVER-call — the other half of the "compression toward 0.5" — is gone** (+0.415 → +0.046, a 9×
cut). **The enriched UNDER-call is essentially unchanged** (−0.169 → −0.152; the `quick` suite reads
−0.169 → −0.187, i.e. flat-to-slightly-worse). Worst cells: `gdna1 ss0.50 verystrong` **−0.301**,
`gdna100 ss0.50 capON` −0.201, `gdna100 ss0.50 verystrong` −0.204.

This is exactly what the dissection predicted: worklog §8's verdict was that the enriched under-call **is the
identifiability ceiling** (κ=½ ⇒ zero strand Fisher information), not a fixable mode defect. The pass-0 work
fixed everything *around* it.

### 2.2 GOAL 1 — the landscape, and the enriched mode is NO LONGER BLIND

EMD-to-oracle: **0.267 → 0.250** cross-suite (ambig 0.209, quick 0.291). A modest move — but the EMD was never
the interesting number. The **enriched-mode census** is:

| stratum | enriched **mass** (fit/oracle) | enriched-mode **location** shift |
|---|---|---|
| **stranded** `ss_0.99` × capON | **0.99 – 1.09** | **−0.03 … +0.03** |
| **unstranded** `ss_0.50` × capON | **0.53 – 0.79** | **−0.12 … −0.23** |
| `gdna1`/`gdna5` × verystrong | 0.21 – 0.34 | −0.90 / +0.03 |

Mean recovery ratio **0.70**. At the pause this was recorded as *"the enriched mode has ~no mass"* — the
landscape was **enriched-blind**. **It is not any more**, and on stranded data the mode is now essentially
exact in both mass and position. The `docs/calibration/figures/resurrect_fit_vs_oracle.png` panels show it
directly: under capture the oracle is visibly bimodal and the fit reproduces both modes.

### 2.3 GOAL 2 — the projection still cannot pull up, and the reason is structural

`enr_recovery = mean(μ* − obs)` over truly-enriched nodes (want **positive**):

| projection | ambig | quick | (2026-07-21) |
|---|---|---|---|
| symmetric, `h=0.15` | **−0.039** | −0.075 | −0.050 |
| symmetric, `h=0.40` | −0.090 | −0.190 | — |
| symmetric, mode read-out | −0.033 | −0.084 | — |
| **asymmetric** (`hup=0.70`, the Band-Aid) | **+0.234** | +0.211 | +0.250 |

`enr_abs_err` improved across the board (0.226 → **0.188**) — the landscape *is* better — but the symmetric
sampling-likelihood projection is still **negative**: it pushes enriched nodes slightly *down*.

**Why, and this is the load-bearing finding.** Substituting the **ORACLE** landscape (same observations, a
perfect prior) recovers only:

| landscape | unstranded | stranded |
|---|---|---|
| fitted | −0.008 | −0.006 |
| **ORACLE** | **+0.025** | **+0.037** |

A *perfect* landscape buys **+0.03 … +0.08** against a **−0.15 … −0.22** under-call. So the projection's
failure is **not** the landscape's quality, and it is not fixable by improving the fit:

> **A symmetric sampling-likelihood projection cannot correct a BIAS. It is variance-shaped.** It asks
> "which DNA level was this observation drawn from?", which is the right question when the observation is
> *noisy* around the truth. Pass-0's enriched error is not noise — it is a systematic, one-directional
> −0.2-decade offset. An observation sitting inside a mode returns `μ* ≈ d` no matter how good the mode is;
> and because the landscape is fitted from those same biased observations, the mode is displaced **with** the
> observation, so the pair is self-consistently wrong.

This is the same structural failure the live docs already name for **P1e** (`variance_ledger.md` §6, "prices
a BIAS as a VARIANCE… a variance cannot move a mode toward truth"). The asymmetric projection works precisely
because it is the only variant that is **bias-shaped** (directional: the observation is a *lower* bracket,
worklog §4). That is why it is a Band-Aid and also why it is the only thing that works.

---

## 3. Two archived claims that the re-measurement REFUTES

1. **"The competitor is NASCENT everywhere"** (worklog §8b). Measured on matched pairs, nascent presence
   makes the unstranded enriched under-call **better**, not worse:

   | cap / dna / ss | nrna present | nrna none | Δ |
   |---|---|---|---|
   | ON / gdna100 / 0.50 | −0.129 | **−0.201** | −0.073 |
   | ON / gdna300 / 0.50 | −0.103 | **−0.166** | −0.063 |
   | ON / gdna100 / 0.99 | −0.075 | −0.005 | +0.070 |

   The under-call is **largest with no nascent RNA at all** on unstranded data, so nascent is not its driver.
   (On *stranded* data the sign flips and nascent does hurt — the archived claim was drawn from a stranded-
   adjacent dossier.) The competitor is unspliced-appearing RNA generally, not the nascent species.
   ⚠ Note the vocabulary rule: RNA is one species; only spliced vs unspliced is observable.

2. **The `fabrication` specificity crisis was half a labelling artifact.** Recorded as drifting −0.83 →
   −0.4…−0.6. Re-measured on the *corrected* zero-gDNA bucket it is **−1.79 (ambig) / −1.67 (quick)** — far
   safer. Running the old buggy bucket reproduces the drift (**−1.05**), confirming the cause: half of those
   "zero-gDNA" conditions (4 of 8) actually contained gDNA. **The deferred "specificity round" is much less
   urgent than the archive says.** The real fabrication risk is narrower and now localised: `none × ss_0.50 ×
   capON` in the quick suite puts a mode at **+0.73** where the oracle has none (mass ratio 4.4) — and its
   stranded twin is clean (0.11–0.17). Same root cause as everything else: unstranded.

---

## 4. The state of play in one paragraph

Fixing pass-0 delivered on **one** of its two promises. The landscape is materially better and the **enriched
mode is recovered** — perfectly on stranded data, at ~70 % mass and ~0.2 decades left-shifted on unstranded.
What survives is a **single directional bias**: on unstranded data pass-0 cannot separate gDNA from RNA inside
the unspliced mass (κ=½ ⇒ zero strand information), so enriched nodes under-call by ~0.15–0.22 decades. That
one bias is now measurable in **three** places and they agree numerically: the under-call table (§2.1), the
enriched mode's left-shift (§2.2), and the projection's inability to pull up (§2.3). **The remaining question
is not "how do we fit a better landscape" — a perfect landscape is worth +0.03. It is "how do we price a
directional bias".**

## 5. The decision the owner has to make (do not build unilaterally)

Two admissible routes, and they are not exclusive:

**(A) Correct the bias at the source — the fragment-length channel.** worklog §8's verdict: gDNA and RNA
fragment lengths are separable, per-fragment, *independently of strand* — the one signal that breaks the
κ=½ degeneracy. It is not a count-zero-information violation (lengths vote, exactly as the strand split
does). ⚠ Its own recorded caveat stands and is serious: **this suite generates fragments at delta-function
lengths (gDNA 100, RNA 200, `frag_std=0`), so the measured separation is a toy artifact** and real libraries
overlap heavily. Flagged at the time as *"a decision point for the owner, not a unilateral build"* — it still
is.

**(B) Price the bias directionally in the projection** — keep `project_asym`, and **derive** `hup`/`cap`
instead of tuning them. This is the archive's own #1 open item. The principled form is visible now in a way
it was not at the pause: the bracket is `obs ≤ truth ≤ log10(mass/E)` (worklog §4), and the *width* of that
bracket is a per-node observable, so `hup`/`cap` need not be constants at all. That would satisfy the
no-magic-numbers rule and remove the last tunables.

**Not admissible:** improving the landscape fit and expecting the projection to follow — measured at +0.03
(§2.3). And per `SESSION_2026_07_26_HANDOFF_15.md` §3 step 2, node **exclusion** rules should be treated
with suspicion: `archive/solve_gate_design.md` records exclusion derived, implemented and empirically
refuted (+0.010 / +0.025).

## 6. Relationship to the live Phase-2 plan

`SESSION_2026_07_26_HANDOFF_15.md` is the live handoff and its Phase-2 plan is **P2a** (feed the hyperprior's
λ-curvature into the DL null) then **P-sub** (the substrate research phase). This document is **upstream
context for P-sub, not a replacement for it** — P-sub asks *which nodes to fit from*, scored in TV against an
oracle-fitted prior. §2.3 here says something P-sub needs to know before it starts: **the substrate choice
cannot fix a directional bias**, so P-sub should be scored on the enriched mode's *mass and location*
(§2.2's census) alongside TV, or it will optimise a number that is already ~83 % saturated while the enriched
mode stays displaced.

Note also that production still ships the **δ-pin** (`DensityNPMLE.logprior`, `gdna_prior_additive=False`) —
the mechanism `gdna_crush_dissection_node1055.md` diagnosed and `gdna_hyperprior_plan.md` §3 proposed to
replace with the projection-plus-gentle-anchor. **That wiring was never done**; it is still open item #3 of
the resume guide's §6.
