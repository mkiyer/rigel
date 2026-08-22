# NEXT SESSION — the owner has a NEW IDEA for calibration. This is the ground it lands on.

    ⚠ **A DEV DOC. Nothing may cite it, and it is NOT the state.** It is a handoff: read it, do the
    work, MOVE what settles into the permanent docs, and DELETE this file in the same session
    (precedent: the three previous `NEXT_SESSION.md` files were each deleted per this instruction).

⛔⛔ **THIS IS DELIBERATELY NOT A QUEUE.** The owner is bringing a new idea for making calibration
work, so the useful handoff is not "do rung 4 next" — it is **what is measured, what is already
refused, and which constraints any new idea has to survive**. Read the sections below in order, then
hear the idea before proposing anything.

## ① WHERE CALIBRATION ACTUALLY IS — the numbers, not the narrative

* **The metric** is the calibration result against ORACLE calibration, per stratum, on the three
  IN-SCOPE strata (`docs/SUCCESS.md` Stage B; `CLAUDE.md`'s SCOPE section). The transcript table is a
  thermometer. ⛔ Never rank on the panel total: 64.5 % of it is the DEFERRED stratum.
* **The exon is the whole problem, and it is structural.** gDNA is directly measurable only where no
  mature transcript crosses — intergenic and intron REGIONs and the BOUNDARIES against them. An exon's
  unspliced mass is gDNA + RNA, which is the unknown itself, so an exon's gDNA level cannot be measured
  and must be imputed from neighbours. That imputation IS message passing
  (`TRAPS: we-keep-re-deriving-message-passing` — four sessions have re-derived this from scratch;
  read that rule before re-deriving it a fifth time).
* **Message propagation is net-HARMFUL on every in-scope contaminated stratum, for both policies**, and
  its value is concentrated where the local solve is BLIND (the deferred stratum, and the capture-ON
  zero controls). `SilentPolicy` wins 8 of 16 ladder rows. That is measured, repeatedly, and it is the
  central open tension in the phase.
* **Which exons are solved well is answered**: the STRAND BIT. Single-strand exons misplace 8.9 / 16.1 /
  15.8 per 1k against AMBIG's 37.5 / 63.2 / 39.4; AMBIG is 18–20 % of exon mass carrying ~50 % of exon
  misplacement (`exon_solvability.py`). The second predictor is the exon's own depth.
* **Strand carries exactly zero information about the gDNA fraction on unstranded data** — the fraction
  cancels from the mean. That is why unstranded × capture-ON is structurally blind and DEFERRED.
* **The measured-prior thread (`ROADMAP.md` §1 rank 3) is three rungs in**: the composition-free
  per-slot TOTAL is built and validated 32/32; two background consumers swapped; the
  `AbundanceLandscape` is fitted, censused on 40 conditions, and is now the QC panel's source. Rung 4
  (the two-atom pass-0 reference) is designed but NOT built.

## ② WHAT THE NEW IDEA MUST SURVIVE — hard constraints, each one measured

⭐ These are not preferences. Every one was paid for.

1. **AXIOM 0 — three populations, never four**: `gDNA`, `RNA+`, `RNA−`, with membership a function of
   two annotation bits. "Mature"/"nascent" is not a degree of freedom. If a design seems to need "a
   third component", re-ask it as *what is this channel's OPPORTUNITY for RNA at this object?*
2. **A TOTAL is composition-vacuous.** No function of a total density can state a gDNA/RNA split. A
   total is a PRE-solve instrument, valid in exactly three places (`DESIGN.md` §3.1a-i).
3. **No magic numbers.** Every divisor derived from the deposit rule and unit-tested against
   brute-force enumeration. A new constant needs its selection evidence or it does not land.
4. **`TRAPS: panel-before-src`** — five recorded times a toy-positive or single-condition-positive
   change was panel-negative. Price on the 16-condition ladder, per stratum, both zero controls, before
   writing a mechanism into `src/`.
5. **An upper bound is not an estimate** (`TRAPS: an-upper-bound-is-not-an-estimate`) and **a mean of
   ratios inherits the partition** — both refuted whole designs here.
6. **Both zero controls on every experiment.** ⚠ `g00` is ONE-SIDED: any RNA-favouring change reads as
   "right" there, so a win at `g00` alone is not evidence.
7. **A refused claim must lose its precision with its value**
   (`TRAPS: zero-the-precision-with-the-value`).

## ③ WHAT IS ALREADY REFUSED — do not rebuild these

⛔ `ROADMAP.md` §4.1 is **THE GRAVEYARD: eleven mechanisms priced and refused**, plus §4.2's seven from
the reference investigation and §4.3's truncation-free bank. Read it before proposing a mechanism — it is
a guardrail that points forward, and it will save a session. The ones that recur as ideas:

* a library-wide composition reference mean (2.65× worse at `g05` stranded × ON);
* BOTH single-mode per-object reference locations, `pooled` and flank-cancellation `local` (5–8× worse
  on stranded × ON, with identical `B exon|exon` damage);
* the soft-min per-transcript allocation rule (worse on every stratum, 1.26–1.60×) — ⭐ but its density
  identity and its opportunity weighting SURVIVE and are reusable;
* the mass-identity message rescale (`k = M/Σρ·E`: amplified a weak source claim 235,800×);
* the fragment-length COMPOSITION channel — ⛔ RETIRED until after 0.8.0; must not be proposed.

## ④ THE STANDING QUEUE, if the new idea is adjacent rather than a replacement

* **Rung 4** — the two-atom mixture through the widened location door (`PLAN_measured_prior.md` §3e),
  consuming `(rho_0, w, enrichment detector)` and ⛔ **NOT `span_R`**, which failed the grid sweep
  (`DESIGN.md` §3.1a-iii). Then rung 5's pricing.
* **Rank 1, the dissection loop** on the worst IN-SCOPE scenario — the owner's own protocol, and the
  item pointed most directly at the 0.8.0 metric.
* Outstanding measurements: `calibration_vs_oracle.py --background-abundance` on the two fl-gap arms;
  the contained-opportunity intron bias (±1.2 %, sign flips with gDNA length — its own ticket).

## ⑤ SESSION MECHANICS AND THREE OPEN CHORES

* FIRST: `python scripts/design/preflight.py` and `python -m pytest tests/ -q`.
  Baseline **0 failed / 3,652 passed / 9 xfail / 3,661 collected** (⚠ this file is not in that count —
  a `docs/dev/` file is `+1` for the jargon gate; re-derive rather than adjust).
* `source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel`; `OMP_NUM_THREADS=1`;
  `python -m pytest`, never bare pytest; export `RIGEL_SCRATCH` in every instrument invocation.
* Falsification first, verified failing; perturb the fixed code; surgical reverts. The owner drives
  commits.
* ⛔ **CHORE 1 — the tree is not `ruff format`-clean under the current ruff**: the documented
  `ruff format src/ tests/` reformats **87 files / ~900 lines** nobody touched. Two sessions have now
  reverted it rather than let it ride along with a measurement change. It wants its own commit, and it
  is an owner decision (`TRAPS: one-thing-varied`).
* ⚠ **CHORE 2** — `/tmp/rigel_ladder_ceiling` (73 GB, a prior session's self-test leak) is still on
  disk, awaiting the owner's word.
* ⚠ **CHORE 3** — `ladder_arm_ab.py --self-test` was interrupted mid-run and has not been seen green
  since; it is a pre-existing instrument nothing recent touched. `preflight.py` covers it — if it ✘,
  check that before assuming a regression.
