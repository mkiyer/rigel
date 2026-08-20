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

⛔⛔⛔ **THE THIRD POLICY IS BUILT AND THE PANEL SAYS IT IS THE WORST ARM ON EVERY IN-SCOPE STRATUM.**
`src/rigel/calibration/messages/currency.py` beat both shipped policies on all three in-scope strata **on
the test chromosome** and is the worst of three **on the 16-condition ladder**, where `SilentPolicy` wins
8 of 16 rows (`ROADMAP.md` §0 carries both tables). ⛔ The toy and the panel disagreed in RANK, not just
in magnitude — `TRAPS: a-toy-and-a-panel-can-disagree-in-rank`.

⭐ **What that does NOT invalidate**, and it is most of the session's value: the four defects it found
were real defects; the honest-precision result (an imputation must cost something every hop) is a
statement about any message layer; the model-free enrichment ratio is a correctness fix for the tool as
a whole; and **the reference finding was measured on the LADDER itself**, so it stands.

⛔ **So the next step is not to polish this policy.** It is the characterisation the owner asked for —
understand the residual well enough to know which exons are safe to train on — plus one specific
question this result raises: **message propagation is net-harmful on every in-scope contaminated stratum
for BOTH policies on this panel.** Its value is concentrated where the local solve is blind. That should
be measured deliberately rather than inferred.

## ⭐⭐⭐ WHAT THE NEXT SESSION DOES

**⓪ FIRST, AND IT IS CHEAP: does the policy help ONLY where the local solve is blind?** Split the
ladder's error by whether the destination had its own composition evidence and read the two halves
separately, per stratum (`TRAPS: an-imputation-must-cost-something-every-hop` has the recipe). If the
whole of message propagation's value is the evidence-free half, then the shippable arrangement may be a
policy that speaks ONLY there — and that is a measurement, not a switch.

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
