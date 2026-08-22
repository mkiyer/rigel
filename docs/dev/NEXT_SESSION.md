# NEXT SESSION — build stages 0–2 of the calibration first-pass redesign

    ⚠ **A DEV DOC, and it is a HANDOFF.** Read it, do the work, MOVE what settles into the permanent
    docs, and DELETE this file in the same session. (Precedent: four previous `NEXT_SESSION.md` files
    were each deleted per this instruction.)

## ① WHERE THE TREE IS

* **HEAD is `bbf4bb2b`** on `fragment-length-gold-standard`; the working tree is clean apart from the
  plan and this handoff.
* **Suite: 0 failed / 3,656 passed / 9 xfail.** `preflight.py` green — and it now runs in **~2 s**
  (the `--self-test` sweep is opt-in via `--full`, which costs ~1 hour and is dominated by one
  instrument). ⛔ Run both before anything else, as always.
* **BOTH PANELS WERE REBUILT WITH SPARSE NASCENT RNA and are certified 16/16** (ladder and test
  chromosome). ⛔ Every quantitative claim in `ROADMAP.md` predates that rebuild unless it says
  otherwise — `ROADMAP.md` now opens with a NUMBERS POLICY that says so and names the instrument for
  each claim instead of embedding a figure.

## ② THE WORK — `docs/dev/PLAN_first_pass_redesign.md`

⭐⭐⭐ **That file is the brief.** It is the owner's stepwise redesign of calibration's first pass,
stages 0–2, turned into a design: the substrate derived from structure (stage 0), the intergenic
initialization and background fit (stage 1), and the two-channel message-free intron deconvolution
(stage 2). The owner is writing the methodology and the later stages in parallel.

⭐ **What the plan already settles, so it is not re-derived:**
* **§A — the load-bearing assumption.** The substrate is a gDNA anchor set *only because nascent RNA is
  sparse*. Every slot it admits is one where MATURE RNA is structurally absent, which is what the
  annotation can prove; what remains is `gDNA + nascent`. This is why the panel rebuild had to come
  first, and why the design must fail LOUDLY when sparsity is violated rather than quietly.
* **§B.2 — the mechanism behind the owner's exon rule**, which a substrate needs or it is a guess: at a
  splice-site boundary with no contiguous exon crossing, the unspliced crossing has no mature term
  (an anchor), the spliced flux is certified RNA (the mature accounting), and the exon shares the
  neighbour's gDNA field. ⚠ I verified the owner's three worked examples against this reading.
* **§D.1 — the precision-weighted combine is LEGITIMATE here, with a derivation.** A slot's total count
  and its strand split are independent statistics (Poisson total, multinomial split given the total), so
  their informations add — which is what `TRAPS: two-gaussians-one-latent` forbids only for channels
  built from ONE latent. ⛔ The implementation must feed each channel its own raw sufficient statistic;
  re-parameterising either in terms of the other destroys the independence.

⛔ **FOUR OPEN DECISIONS are in §G and the owner should rule on them before much code is written**:
`the-one-hop` (does the exon's single structurally-terminated hop use the message framework, or is it
a different thing? — `TRAPS: one-hop-lifted-out-is-still-the-relay` says be careful), `substrate-home`
(index time and persisted, or derived at calibration init? — the plan recommends derived, and says why),
`clamp-semantics` (is the `ρ_bg ≥ ρ_total` clamp a location clamp or does it cap the strength too?), and
`combine-form` (lattice product with per-channel diagnostics, or two point estimates? — the plan
recommends the former with the latter published as diagnostics).

## ③ THE TWO RUNNING LOGS THE OWNER ASKED FOR

They live in the plan's **§F** and must be updated as each stage lands, not at the end: **§F.1** the code
this design makes OBSOLETE (with the evidence it is superseded and what must be re-priced before it can
go), **§F.2** the code that gets SIMPLER. The owner's stated goal is a simpler, cleaner codebase as an
outcome of this redesign, not a side effect.

⚠ **A grounding sweep was launched for exactly this and had not finished when the session ended.** It is
saved and re-runnable as one command — `Workflow({scriptPath: "<session-dir>/workflows/scripts/`
`ground-calibration-redesign-wf_40b9929c-ba0.js"})`, or just re-author it; the script's own prompts are the
specification. It maps, with file:line: the structural/statics layer and
whether the solvable-exon predicate is expressible from existing bits; the strand channel's per-slot
contribution and precision; `fit_intron_background`, the NegBinom intron factory and `total_abundance`;
psi's combine and the existing precision-weighted idiom; and the first draft of both logs. Re-run it, or
re-derive by hand, BEFORE writing code — the design's file-level claims are not yet verified.

## ④ WHAT IS ALREADY KNOWN TO BE CLOSE TO STAGES 1–2

⭐ **Stage 1 largely exists**: `g1_locked` already pins intergenic REGIONs and gene-edge BOUNDARIEs at
`f_g = 1`, variance 0, and `fit_intron_background` already pools intergenic REGIONs only (both live call
sites pass `include_introns=False`). ⛔⛔ **But the emission lock is scoped `~solvable ∧ REGION` where the
standing strict xfail says it should be `g1_locked ∧ REGION`, and rescoping it ALONE was priced and
REFUSED** — it is half of a recorded cancelling pair (`ROADMAP.md` §4.1). Making it structural is more
principled AND is not free; price it with whatever replaces the load the mis-scoped mask carries.

⭐ **Stage 2's density channel largely exists** as `calibrate.py`'s per-intron
`log NegBinom(f_g·C; ρ_bg·E_g, α_eff)` factor, applied to intron REGIONs only, whose curvature already
becomes a precision. Read it before building a second one.

## ⑤ THE VALIDATION BAR (plan §E) — in this order, nothing quoted before it

① a falsification test written FIRST and verified failing, then every gate fired by perturbing the fixed
code; ② stage 0 checked as a CONFUSION MATRIX against certified per-slot truth — it needs no solver, and a
"solvable" slot whose truth disagrees falsifies the predicate immediately; ③ per-stratum scoring
restricted to the claimed slots only, both zero controls included; ④ test chromosome first, then the
16-condition ladder (`TRAPS: panel-before-src`).

⛔ **JUDGE ONLY WHAT IS ACTIVELY SOLVED** (owner). While these stages are in flight, ignore the slots the
stage does not claim; a pooled number would hide both the win and the damage.

## ⑥ TWO OPEN CHORES AND ONE CAVEAT

* ⚠ The tree is not `ruff format`-clean under the current ruff (~87 files reformat). Three sessions have
  now declined to let it ride along with a measurement change. It wants its own commit — owner's call.
* ⚠ `/tmp/rigel_ladder_ceiling` (73 GB, an old self-test leak) is still on disk awaiting the owner's word.
* ⛔ **`object_composition.PURE_GDNA_STRATA` includes `R intron`**, which was exactly pure only while the
  panel held no nascent RNA and is now contaminated worst where gDNA is scarce. Recorded in the
  instrument and stamped in `ROADMAP.md`, deliberately NOT changed — the pool is a measurement-design
  decision, and the recorded exon-reference result came through it.
