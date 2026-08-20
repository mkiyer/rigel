# NEXT SESSION — characterise the residual, then evolve the baseline

    ⚠ **A DEV DOC. Nothing may cite it, and it is NOT the state.** The state is `ROADMAP.md` §0/§1, the
    rulings `DESIGN.md` §0c, the derivations `EQUATIONS.md` §3.5f–h, the lessons `TRAPS.md` §D, the
    substrate `TESTING.md` §0a/§0b. ⛔ Delete this file when the plan it points at is executed.

    Rewritten 2026-08-20, when the third policy became the baseline.

## ⭐⭐⭐ FIRST COMMAND OF THE SESSION — prove you can run and regenerate everything

    python scripts/design/preflight.py          # ~2 min, or --fast to skip the instrument sweep
    python -m pytest tests/ -q                  # the baseline is in CLAUDE.md; ANY failure is a regression

⛔ **Do not start work on a ✘.** A missing DERIVED artifact is a command you have not run, not damage.

## WHERE THINGS STAND — read `ROADMAP.md` §0 for the numbers, this for the shape

⭐⭐⭐ **THERE IS A BASELINE THAT BEATS BOTH SHIPPED POLICIES ON EVERY IN-SCOPE CONTAMINATED STRATUM**
(owner, 2026-08-20: *"for the first time … we have a baseline that is beating all of the prior existing
implementations. This is a win … This is the beginning, not the end."*). It is
`src/rigel/calibration/messages/currency.py`, selected by `CalibrationConfig.message_policy`, with
`RelayPolicy` untouched beside it.

⛔ **It is a baseline to EVOLVE, not a finish.** The owner's framing for what comes next is
characterisation: understand the residual well enough to know **which exons are safe to train on**, so
the tool can be trusted on real data where there is no ground truth.

## ⭐⭐⭐ WHAT THE NEXT SESSION DOES

**① CHARACTERISE WHICH EXONS ARE SOLVED — by OBSERVABLE properties, never by truth.** This is the
training substrate for the gDNA landscape, and the whole bootstrap (`DESIGN.md` §0c.0d) turns on it: the
prior cannot be trained until enough exons are solved, and once it is trained it subsumes much of what
message propagation does. ⛔ Report a CURVE or a distribution, never a tolerance — a threshold here is a
magic number. Candidate observables, all available without truth: own composition evidence
(`has_own_composition_evidence`), hops from the nearest structurally pure-gDNA slot, flanked by a splice
site vs a terminus, one strand vs both, observed depth.

**② THEN EXPAND THE gDNA SPECTRUM — deliberately** (`ROADMAP.md` §1 rank 1b). Levels are informative
where behaviour CHANGES; interpolation is free. ⛔ Do not multiply the full cross.

**③ AND KEEP GROWING THE TEST CHROMOSOME** — it is becoming a gold-standard regression set (owner). Add
a case whenever a new stressing structure is conceived. Queued: a cassette (skipped) exon; a same-strand
intron-retention isoform; divergent and convergent gene pairs; a probed exon shared by isoforms whose
other exons are unprobed.

**The loop is ~1 min end to end** (`docs/TESTING.md` §0a has the commands; `<T> = ~/Downloads/rigel_runs/test_reference`).

## ⛔ THE MEASUREMENT DISCIPLINE THAT FOUND EVERY DEFECT THIS SESSION

**Split the error by whether the destination HAD ITS OWN composition evidence.** "The messages help" and
"the messages trample a measurement" are different findings, and a pooled number cannot tell them apart —
it is how an 8.09× degradation of the measured half of the chain was found hiding inside a modest total
(`TRAPS: an-imputation-must-cost-something-every-hop`).

⚠ And verify that a perturbation LANDED before believing a green result. Three separate times this
session an ablation silently did nothing: a `sed` that matched no line, a gate that read one slot past
the stretch it meant to isolate, and — the expensive one — a patched module binding that
`region_init` does not use, because it holds its own reference to the solver
(`TRAPS: an-ablation-that-never-ran`).

## STILL OPEN

* ⭐ **The capture-aware reference** (`ROADMAP.md` §1 rank 3 + the §0 reference row): a per-object mean
  is measured to take both zero controls to exactly 0 and to improve every capture-OFF row on the
  ladder, and to regress every stranded × capture-ON row by the documented anchor under-read. Rank 2 and
  rank 3 are ONE piece of work.
* `relay_pool_ab.py`'s docstring promises a `--table pipeline` that does not exist. Build it or remove
  the promise.
* Prose in ~18 instruments still describes the retired 36-condition ladder. Historical measurement
  stamps are PROVENANCE and must stay; what is worth fixing is prose presenting a retired panel as
  current.
