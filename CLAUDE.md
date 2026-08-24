# CLAUDE.md

Guidance for Claude Code working in this repository.

## What Rigel is

A Bayesian RNA-seq transcript quantifier that separates **RNA from genomic-DNA contamination**. A
single-pass C++ BAM scanner tallies fragments, a **calibration** stage deconvolves the library into gDNA vs
RNA, and a per-locus EM solver assigns RNA to transcripts. PyPI package `rigel-rnaseq`; the import and CLI
are `rigel`.

⭐⭐⭐ **The release target is 0.8.0 and CALIBRATION is the focus.** The SCOPE section below is the frame
for every change in this phase. ⭐ Read it with `DESIGN.md` §0b's **NASCENT SCOPE RULING** (2026-08-22):
nascent RNA is sparse in real data and is modeled for ROBUSTNESS, so a number measured at the panel's
nascent share (20.2 % capture-OFF and 2.6 % capture-ON) is a stress reading and never a design driver.

⭐⭐ **THIS FILE IS DIRECTIONS, NOT A HISTORY.** It was torn down from 732 lines on 2026-08-22 because the
accumulated rules had started to inhibit design rather than protect it. ⛔ Keep it that way: a rule earns
its place by changing what the next session DOES, its evidence lives in the doc that owns it, and nothing
here is a changelog. Do not add a rule the owner did not ask for.

## ⛔⛔⛔ AXIOM 0 — RNA IS RNA. READ THIS BEFORE ANYTHING ELSE

**There are THREE populations and there is no fourth: `gDNA`, `RNA+`, `RNA−`.**

⛔ **"Mature" and "nascent" are NOT populations, NOT species, and NOT a degree of freedom.** RNA inside an
intron is RNA that has not spliced *at that position*. The only distinction that exists is whether a
fragment **is spliced** — certified RNA, gDNA cannot splice, needs no deconvolution — or **is not**, which
is the entire deconvolution problem.

⭐ **The set is a function of TWO BITS, which is what makes this structural rather than something to
remember:**

```
T(slot) = {gDNA}                                  # always — gDNA is genomically continuous
        ∪ {RNA+ if statics.free_pos[slot]}        # iff the annotation admits that strand here
        ∪ {RNA− if statics.free_neg[slot]}
```

so `|T| ∈ {1,2,3}`, always. **No expression in this codebase may produce a fourth population.**

⛔ **The tell that you have just violated this:** you have written a population set with more than three
members, or the words "mature"/"nascent"/"subspecies" while describing a *solver*, *composition* or
*population* question, or you have concluded that something needs "a third component". Stop, and re-ask the
question as **"what is this channel's OPPORTUNITY for RNA at this object?"** — the answer is a geometry, it
is derivable from the index, and it dissolves the objection every time.

⚠ The words survive as **simulator inputs and never as model concepts**: the simulator's
`nrna_abundance` knob and the toy harness's `--nrna` arm, which exist to put RNA inside introns so the
solver can be tested on it. ⭐ How much WEIGHT a nascent concern carries is the NASCENT SCOPE RULING's
question, not this axiom's — see the frame above.

## ⛔⛔⛔ THE 0.8.0 SCOPE — READ THIS BEFORE CHOOSING WHAT TO WORK ON

Owner ruling, 2026-08-14. The version on disk is **0.7.1** (`pyproject.toml`); the target is **0.8.0**.
This is the frame every doc, every experiment and every ranked list hangs off.

⭐⭐⭐ **THE FOCUS IS CALIBRATION, AND THE METRIC IS THE CALIBRATION RESULT SCORED AGAINST ORACLE
CALIBRATION — not the end-to-end transcript number.** The transcript number stays a thermometer
(`SUCCESS.md` says why); a change is judged by `solvability_audit.py` and `prior_vs_oracle.py`.

**THE FOUR STRATA — THREE ARE THE TARGET AND ONE IS NOT:**

| stratum | 0.8.0 |
|---|---|
| unstranded × capture-OFF | ⭐ **IN SCOPE** |
| stranded × capture-OFF | ⭐ **IN SCOPE** |
| stranded × capture-ON | ⭐ **IN SCOPE** |
| unstranded × capture-ON | ⛔ **DEFERRED** |

⛔⛔ **DEFERRED IS NOT DROPPED, AND THAT DISTINCTION IS THE WHOLE RULING.** unstranded × capture-ON
**remains in every benchmark and every measurement, and must keep being reported**. It is simply not a
development target until the other three are fully optimised. If it improves as a side effect of other
work, that is a free win.

⚠ **The deferred stratum is also where the error is** — 64.5 % of transcript error and 90 % of
gene-level error on the rebuilt ladder — because the tool emits a near-zero gDNA fraction there
regardless of truth, which looks acceptable at low gDNA by coincidence. ⛔ So the debug loop points there
every time: take the worst **IN-SCOPE** scenario instead. The algebra behind the blindness is that the
gDNA fraction cancels from the strand mean, so an unstranded AMBIG slot has no channel at all.

⛔⛔ **THE LENGTH CHANNEL IS RETIRED UNTIL AFTER 0.8.0 — DO NOT PROPOSE IT.** The fragment-length
likelihood **as a CALIBRATION composition channel** may not be listed or ranked. ⚠ It does not exist in
`src/`: it was A/B'd once and purged the same day, so this is a scope ruling and not a deletion.
⭐ **Three other things called "length" are NOT affected**: layer 2's `fl` / `effective_length` /
`capture_eff_length` (the OPPORTUNITY model, which calibration and the EM both need), `length_likelihood`
in `src/rigel/second_pass.py` (per-fragment assignment, same word, different thing), and the fl PMFs
priced by `length_ceiling.py`.

⭐⭐ **THE LADDER GIVES gDNA AND RNA EQUAL FRAGMENT LENGTHS BECAUSE THE EM ALREADY READS THE FL
DISTRIBUTION.** A large length gap lets the EM split the origins on LENGTH ALONE, bypassing calibration
and masking its bugs; equal lengths FORCE calibration to be exercised, leaving it strand, density and the
message layer. The fl-gap SIDE panel is what exerts the length mechanism deliberately.

⭐ **SCENARIOS ARE CACHED, and that is a scope requirement**: calibration must re-run in seconds off a
cached scan, so `panel.py cache` / `build_scan_cache.py` come before any experiment. ⭐ **The ranked list
of what to do next is `ROADMAP.md`** — this section says what counts, not what is next.

## ⭐ The docs — read them in this order

Six permanent docs. They do not overlap, and none of them is a changelog.

| doc | what it is |
|---|---|
| ⭐⭐ **`docs/SUCCESS.md`** | **START HERE.** How performance is measured: **Stage A** (the accumulator — faithful, unbiased, sufficient?) and **Stage B** (calibration against an oracle). ⭐ 0.8.0's metric is Stage B's |
| ⭐⭐ **`docs/ROADMAP.md`** | Where the tool is in measured numbers, and what to do next in priority order. ⛔ **Not a changelog — a completed item is DELETED, not stamped** |
| ⭐ **`docs/TRAPS.md`** | Mistakes already made. Lessons, not measurements — cite a rule by its NAME |
| **`docs/EQUATIONS.md`** | The derivations the code depends on |
| **`docs/DESIGN.md`** | What is built and the rulings behind it — settled, do not re-litigate. ⭐⭐ **§0 is the binding VOCABULARY** (REGION, BOUNDARY, step, structurally pure-gDNA object) and **§0b carries the 0.8.0 SCOPE and the NASCENT SCOPE RULING** |
| **`docs/TESTING.md`** | The panels, how to build them, the simulator's gates, and what the suite can and cannot judge. ⭐⭐ **§0b is the TOY HARNESS** |

Reference rather than design: `docs/MANUAL.md`, `docs/PUBLISHING.md`.

⛔ **THE MOVE RULE, which governs every doc and `docs/dev/` alike.** When a finding settles, MOVE it to its
one home — a current number to `ROADMAP.md` §0, a lesson to `TRAPS.md` as a NAMED rule, a ruling to
`DESIGN.md`, a derivation to `EQUATIONS.md` — and delete it from where it was, in the same edit. ⛔ **Move,
never copy**: copying is what creates two homes that then diverge. ⚠ Both halves of this have been paid
for — `ROADMAP.md` reached 787 lines under "never delete, only stamp", and two dev docs reached 1,181 lines
and were being read as THE STATE.

⭐⭐ **`docs/dev/` IS THE SANDBOX AND IT IS ENCOURAGED** — working notes, a half-finished argument, a
handoff. ⛔ Nothing there is authoritative and nothing may cite into it: not the source, not a test, not a
permanent doc. It is expected to be provisional and occasionally wrong.

⛔ **THE SOURCE DOES NOT CITE THE DOCS.** Docs rot — 73 % of the citations once in the source pointed at
deleted documents. A docstring may cite a TEST or the executable specification
(`tests/native/_accumulator_reference.py`, which wins over any document about the accumulator), never a doc.

## ⭐⭐⭐ WHERE DOES A CHANGE GO? — the calibration layering

**THE ONE RULE: an import may point DOWN a layer or SIDEWAYS within one, never UP.** A module reaching for
something one layer up is telling you the thing belongs lower — a TYPE almost always does.
`rigel/calibration/_layers.py` is authoritative and `tests/calibration/test_layering.py` enforces it.

| if the change is about… | it goes in |
|---|---|
| what a fragment tally MEANS | **1 · the payload view** — `splice_graph` `substrate` `region_arrays` |
| how many places a fragment COULD have sat | **2 · opportunity** — `effective_length` `capture_eff_length` `sj_opportunity` `gdna_opportunity` `fl` |
| one slot's own numbers, ψ, and its total | **3 · geometry + the per-slot solve** — `region_geometry` `simplex_logodds` `total_abundance` |
| which strand a fragment came from | **4 · strand** — `gdna_strand` `strand_deconv` `strand_balance` `strand_summary`, and `strand_likelihood` (a gated executable REFERENCE) |
| how dense a component is, and the priors | **5 · density and prior** — `density_model` `density_deconv` `landscape` `abundance_landscape` `run_fill` |
| what one neighbour tells another | **6 · the solve** — `sweep` (the backbone) + `messages/` (the policy) + `region_init` |
| turning the solve into a result | **7 · assemble** — `calibrate` `priors` `result` `derive` `diagnostics` `track` |

⭐ **Run `python scripts/design/module_census.py` rather than trusting this table** — it re-derives the
graph from the AST, flags any upward import, and flags docstrings naming a sibling with no import. ⚠ A layer
is not a claim that its modules are the right SIZE; layer 4 being five modules for one concept is open.

## ⛔ THE MESSAGE LAYER — where it stands, in one place

⭐ **Message propagation is ON** (`CalibrationConfig.message_propagation = True` since 2026-08-18), and the
policy is one config value, `message_policy`, default `relay`. Three policies sit on one gated backbone:
`SilentPolicy`, `RelayPolicy`, and `CurrencyPolicy` (`messages/currency.py`).

⛔ **`CurrencyPolicy` is a DEVELOPMENT BASELINE WITH A MEASURED DEFICIT, not an improvement** — it wins
every zero control and loses all three in-scope contaminated strata, where `SilentPolicy` still wins.
⛔ **On the test chromosome it wins those same strata, so every claim about it MUST name its substrate**
(`TRAPS: a-toy-and-a-panel-can-disagree-in-rank`). The tables are `ROADMAP.md` §0; the rulings are
`DESIGN.md` §0c and `EQUATIONS.md` §3.5f–h; the reference-first ordering — fix ψ's reference, THEN
re-contrast the three policies — is `ROADMAP.md` §1.

⛔ **`RelayPolicy` is being REBUILT rather than repaired bug by bug** (owner, 2026-08-18), so do not open a
per-bug dissection of it. Its known defects are kept as CONSTRAINTS on the replacement, chief among them
`TRAPS: zero-the-precision-with-the-value`.

⭐ **The vocabulary rename is in flight.** Landed: `graft` → SPLICE IN, the operator sense of the old
peel token → SPLICE OUT, the deconvolution verb → `deconvolve`, `mass_pin` → `mass_rescale`,
`HeadPolicy` → `RelayPolicy`, `lend` → `may_share_composition`, `s2t` → `hop_logvar`. ⛔ **Still to do and
DEFERRED for a measured reason: `arm` has THREE senses** — an experiment arm, a component arm (the sense
the owner ruled becomes `component`), and `__ARM_NEON` in the C++ scanner — across 1,358 sites. It needs
its own `rename_census.py --sense` pass, never a tail-end sweep. ⛔ Run `rename_census.py --sense <token>`
before renaming anything and `rename_identity.py --check` after each stage.

## ⛔⛔ CITE A RULE BY ITS NAME. NUMBERED LABELS ARE BANNED

Cite a trap as `TRAPS: off-grid-message-mode`, never as `A16`. The name IS the identifier, so a citation
says what it means without a lookup and stays one greppable string with one home. ⛔ The old numbers were
not merely opaque, they were AMBIGUOUS — `G1` meant "no magic numbers" AND "a structurally pure-gDNA
object" across 201 sites. `tests/test_no_jargon_labels.py` enforces this and carries the history; its
allowlist is scoped, never blanket.

## Working rules

- **No magic numbers.** Stop and discuss before adding any constant, heuristic or tunable. Every divisor
  must be derived from the deposit rule and unit-tested against brute-force enumeration.
- **A falsification test first, verified failing — then break the fixed code and watch each gate fire.**
  The second half is not optional; it has found holes in already-green gates repeatedly.
- ⭐⭐ **THE DEBUG LOOP IS THE DEFAULT METHOD**: run the panel → take the worst **IN-SCOPE** scenario
  (never the deferred stratum, which is otherwise worst every time) → dissect it to the highest-error
  objects → find the mechanism → fix → re-run the panel. It produced the 39 % win.
- ⭐ **A ceiling is sometimes the right instrument** (`calibration_truth_ab.py --ceiling`) and has
  re-ranked the project twice, but it prices something that may be unreachable, so it is not the default
  (owner, 2026-08-05). ⛔ Check which RULER a ceiling left installed: every arm patches `assemble_priors`
  while `pipeline.py` builds `effective_lengths_em` BEFORE calling it, so the effective-length shrinkage
  has never been inside any ceiling number.
- **One thing varied per experiment**, and a baseline re-recorded from the current tree in the same
  session.
- **Score against TRUTH, not against the previous run.** The simulator writes per-fragment ground truth
  into the oracle BAM's read names.
- **No legacy, no backwards compatibility, no speculative code.** Converge and delete. Code kept "for
  comparison with the old version" is a defect. ⛔ No version suffixes in file names — it is
  `accumulator.py`, never `accumulator_v5.py`.
- **No Greek letters in identifiers** (fine in maths write-ups).
- ⛔ **Real data is a TEST input, NEVER a DESIGN input.** The cfRNA on disk is one far end of the RNA-seq
  spectrum, not a sample of it. Sweep the plausible space, report the worst case, and bring the owner the
  domain call.
- **Profile on real data, never a small synthetic suite** — a toy ranks hotspots backwards. ⚠ For the
  0.8.0 performance work the target is HIGH-DEPTH real RNA-seq rather than cfRNA, which is too sparse to
  optimise against (owner, 2026-08-17).
- **The owner drives commits.** Do not commit unless asked.

## Build, test, lint

> **Every build, test and lint command runs inside the activated `rigel` conda environment** — it holds
> htslib and the compilers, and the C++ build finds htslib via `$CONDA_PREFIX`.

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel

pip install --no-build-isolation -e ".[dev]"   # rebuild after ANY src/rigel/native/ change
python -m pytest tests/ -q                     # ⛔ never bare `pytest` — the repo root leaves sys.path
python -m pytest tests/ --update-golden        # regenerate tests/golden/ after intended output changes
ruff check src/ tests/ scripts/ && ruff format src/ tests/   # ⚠ NEVER format scripts/
```

⭐ **THE STANDING BASELINE: 0 failed / 3,653 passed / 0 skipped / 8 xfail** (re-derived 2026-08-24,
sixth measurement that day. The `−29` over 3,682 is the FAN-OUT DELETION (`DESIGN.md` §6b.2 ✅):
`tests/calibration/test_fanout_policy.py` removed — 27 own cases plus its 2 per-file gate cases —
after the anchored-tree re-contrast measured `FanOutPolicy` dominated on every row; no shipped
behavior moved. Before that, `+2` over 3,680: two estimator-hardening gates in the EXISTING
`tests/calibration/test_rna_anchor.py`, no file moved. The `+16`
over 3,664 is the RNA-ANCHORED EVIDENCE FACTOR (`DESIGN.md` §6b.2): `src/rigel/calibration/rna_anchor.py` (+3: jargon, docs-boundary, layering)
and `tests/calibration/test_rna_anchor.py` (11 own cases +2 gate cases). Goldens did NOT move —
the golden toys sit below the factor's dispersion-estimator population minimum, where it is
deliberately flat. The prior `−64`: the REFERENCE-LOCATION DELETION (owner
refutation; `DESIGN.md` §6b.1): three gate files removed — `test_structural_reference.py`,
`test_measured_reference.py`, `test_reference_location.py` — 58 own cases plus 3×2 parametrized
gate cases. The `−1` xfail is the relay-transports-a-structural-claim strict xfail that lived in
the first of those files; the DEFECT stays recorded as a relay-rebuild constraint in `DESIGN.md`
§0c.2/§6b.1. A golden update landed in the same change — the deletion moves shipped output; diff
magnitudes and the truth-scored before/after live in `DESIGN.md` §6b.1 and `ROADMAP.md` §0.)
⛔ **ANY failure at all is a regression** — a stronger and
cheaper rule than counting the expected ones. ⚠ A commit that measures the suite updates this line, or the
next session reads a green run as a regression.

⛔ **RE-DERIVE A COUNT, NEVER ADJUST ONE** (`TRAPS: re-record-the-baseline`). Adding or retiring a file
moves the total, because several gates are PARAMETRISED over the files on disk. Account for the delta from
this table, then confirm it with `pytest --collect-only -q | grep <stem>` and EXACT bracket matching:

| adding one… | moves collected by | which cases |
|---|---|---|
| `src/rigel/calibration/` module | **+3** | jargon, docs-boundary, and layering *if declared in `_layers.py`* |
| `tests/calibration/` file | **+2** | jargon, docs-boundary (never +3 — `test_scripts_index` does not parametrise it) |
| top-level `tests/` file | **+3** | jargon, docs-boundary, scripts-index |
| `scripts/design/` (or `sim/`, `profiling/`) file | **+4** | imports, says-what-it-is-for, jargon, docs-boundary |
| `docs/dev/` file | **+1** | jargon only |

⭐ **A content-only sweep must move the collected total by ZERO, and that is the check** — edits to
comments, docstrings and prose add no cases, so a moved total after one means a file was added or removed.

⛔ **DERIVE THE FAILURE SET, NEVER EYEBALL THE TAIL** (`TRAPS: read-the-whole-failure-list`) — pytest
prints the last screen, and 21 reported "golden failures" once hid two real ones::

    python -m pytest tests/ -q 2>&1 | grep '^FAILED' | sed 's/::.*//' | sort | uniq -c

⛔ **A GOLDEN UPDATE IS WHERE A REGRESSION GETS LAUNDERED INTO "INTENDED".** Read the diff and record its
magnitude BEFORE `--update-golden`, and check the truth-scored instruments (panel, thermometer) first.

⛔ **RUN THE INSTRUMENTS, NOT ONLY THE SUITE, after a `src/` DELETION, a RENAME, or a SHIPPED-DEFAULT
FLIP** — `preflight.py --full` does it in one command. A green suite has hidden five dead instruments twice
(`TRAPS: a-green-suite-hid-five-dead-instruments`), because the tests install what the shipped default
does not.

⭐⭐ **THE 8 xfails ARE NOT ONE KIND OF THING.** Some are the recorded price of a config default; the rest
are executable records of proven defects whose fixes are panel-negative alone — for those, "fix the test"
is a category error, because the test is right and the code is wrong. ⛔ An xfail is closed by REPAIRING
the thing or by asserting the invariant STRUCTURALLY, never by widening a bound.

Always set `OMP_NUM_THREADS=1` when benchmarking or comparing runs.

## Tooling under `scripts/`

This table indexes `scripts/design/` plus four `sim/` rows and nothing else — ⛔ **an absence is NOT a
deletion.** The eight unindexed `design/` files are named in `tests/test_scripts_index.py`'s
`UNDOCUMENTED_DEBT` (an entry is a decision owed, and the list may only shrink); `scripts/profiling/` is
covered by the same gate and described in `scripts/README.md`. ⭐ **Each row carries only the QUESTION its
instrument answers** — verdicts and long cautions live in the instrument's own docstring, which a suite gate
requires. Groups are ordered by 0.8.0 priority; `docs/SUCCESS.md` has the run order.

| | |
|---|---|
| **⭐⭐⭐ START A SESSION HERE** | |
| `design/preflight.py` | ⭐⭐⭐ **CAN THIS SESSION RUN AND REGENERATE EVERYTHING? — one command, one verdict, before anything else.** Checks the toolchain (the `rigel` env, the native extension, the CLI), both references, both panels (scan caches, oracle caches with all five partitions, the certified `slot_truth`) and that every `scripts/design/` instrument IMPORTS. ⛔ It changes nothing and measures nothing — every check is a read or an import, and a ✘ prints the exact command that regenerates the missing artifact. ⭐ **The default is ~2 s**; `--full` adds every instrument's `--self-test` and costs **~1 hour measured** (7.5 CPU-hours, dominated by `ladder_arm_ab.py`) — run it after a deposit-rule change or a default flip, never every session. `--self-test` 8/8 |
| **⭐⭐⭐ 0.8.0'S METRIC — calibration against ORACLE CALIBRATION** | |
| `design/calibration_vs_oracle.py` | ⭐⭐⭐ **IS THE CALIBRATION RESULT ITSELF RIGHT, SCORED AGAINST AN ORACLE CALIBRATION? — 0.8.0's metric, and the only instrument that reaches the effective-length shrinkage.** `P = calibrate(...)` against the same payload with only the six deconvolved arrays swapped, per stratum, plus `U`, the no-enrichment null no other instrument carries. ⛔ Read `ruler_n_moved`, never the aggregate: the total can barely move while nearly every transcript is redistributed. No solver, no EM, no re-scan — ~5–12 s/condition. `--self-test` 21/21 |
| `design/object_composition.py` | ⭐⭐⭐ **MUST ψ's Beta REFERENCE BE ONE LIBRARY-WIDE NUMBER, OR CAN EACH OBJECT SUPPLY ITS OWN?** `m_i` per object from the two densities, scored as misplaced fragments against the shipped ½, per stratum. `--self-test` 25/25 |
| `design/abundance_landscape_census.py` | ⭐⭐⭐ **WHAT DOES THE TOTAL-DENSITY FIELD LOOK LIKE, PER CONDITION?** Fits `calibration.abundance_landscape` on the cached wall-exact totals — every mode's basin mass, `rho_0`, the anchor gap in nats, the per-class enrichment. `--self-test` 13/13 |
| `design/landscape_head_to_head.py` | ⭐⭐⭐ **HOW GOOD IS THE TOTAL-DENSITY LANDSCAPE, AND WHAT IS ITS BANDWIDTH REALLY?** The axis offset, a held-out predictive likelihood, and `--grid-sweep`. ⛔ `span_R` and the mode count are grid artefacts. `--self-test` 42/42 |
| `design/transfer_variance_audit.py` | ⭐⭐⭐ **WHERE DOES THE MESSAGE TRANSFER VARIANCE COME FROM, AND COULD THE `AbundanceLandscape` SUPPLY A BETTER ONE?** An audit: it proposes nothing and changes nothing. `--self-test` 20/20 |
| `design/calibration_oracle.py` | ⭐⭐⭐ **WHAT IS THE CERTIFIED PER-OBJECT TRUTH? — run this before debugging calibration against anything.** Every REGION and BOUNDARY's count, its realized `n_gdna`/`n_nrna`/`n_mrna` and `true_f_g`, at two certification levels: COMPOSITION (no opportunity model anywhere in it) and FIELD (densities too). ⛔ REFUSED unless its named gates pass — sum-to-full, partition-projects-exactly, gdna-field-uniformity, exact-zeros, nascent-in-annotation — because a merely plausible oracle is how a calibration bug and a truth bug survive each other. Writes `slot_truth.npz` beside each oracle cache; `--self-test` 11/11 |
| `design/total_abundance_audit.py` | ⭐⭐⭐ **IS THE MEASURED TOTAL A TRUE TOTAL?** Five arms against the origin partitions; read ⓔ START/END agreement first — the only field-free arm and the decisive test of the wall rule. `--self-test` 15/15 |
| `design/calibration_walk.py` | ⭐⭐⭐ **WHICH STAGE OF CALIBRATION INTRODUCES THE ERROR?** The solve as a ladder — init → strand → local → +messages → +refits → shipped — each rung scored per stratum against `calibration_oracle.py`, which it refuses to run without |
| `design/structural_claims_audit.py` | ⭐⭐⭐ **IS EVERY SLOT THE STAGE-0 SUBSTRATE ADMITS TRULY WHAT IT CLAIMS? — the confusion matrix against certified slot truth, no solver.** Each structural class scored on ITS OWN claim in fragments; the solvable-exon claim is tested at the licensing FLANK, and nascent inside an ss intron is not a violation. ⛔ REFUSED without `slot_truth.npz`. `--self-test` 8/8 |
| `design/pass0_claimed_ab.py` | ⭐⭐⭐ **HOW WELL DOES PASS-0 SOLVE THE SLOTS IT CLAIMS, PER POLICY?** silent/relay at the stage-0 substrate's two claimed populations (`ss_intron_boundary`, `solvable_exon`), misplaced gDNA fragments vs certified truth, DELIVER/REFUTE split and never pooled. ⛔ A whole-library number cannot judge pass-0 — that context is `calibration_vs_oracle.py`. ⚠ The `--dissect` survey died with `FanOutPolicy` (2026-08-24); its verdicts live in `DESIGN.md` §6b.2. `--self-test` 6/6 |
| `design/relay_pool_ab.py` | ⭐⭐ **WHAT DOES MESSAGE PROPAGATION DO, OFF vs ON, per condition and per pool?** Signed and misplaced-mass errors in fragments against origin-split truth, never collapsed, both arms in one process off one cached payload. `--self-test` 11/11 |
| `design/benchmark_report.py` | ⭐⭐ **WHAT DOES THE WHOLE BENCHMARK LOOK LIKE ON ONE HTML PAGE? — every scenario in counts, pooled only on the last row.** ⛔ It scores nothing: it renders `relay_pool_ab.py --out`. `--self-test` 10/10 |
| `design/hop_currency.py` | ⭐⭐⭐ **WHICH CURRENCY DOES EACH HOP TYPE CARRY — A LEVEL OR A COMPOSITION?** Every adjacent pair keyed by `object class × {sj, term}`, the source's true value transported both ways and scored against a Monte-Carlo noise floor. `--self-test` 36/36 |
| `design/solvability_audit.py` | ⭐⭐⭐ **WHICH OBJECTS ARE SOLVABLE, WHICH ARE SOLVED WRONG, AND WHICH ARE CONFIDENTLY WRONG? — where pass-0 and 0.8.0 are judged.** ⛔ Honest ignorance is excluded: `f_g ≈ ½` at zero precision with no own evidence is correct. `--suite` runs the panel |
| `design/prior_vs_oracle.py` | ⭐⭐⭐ **IS `LocusPriors` — the thing the EM actually reads — RIGHT?** Five arms separate calibration's own error from the assembler's, reporting the count, the composition claim and the scale apart, per stratum. ⛔ Undrained on every arm |
| `design/pass0_vs_oracle.py` | **HOW DOES PASS-0 COMPARE WITH THE ORIGIN-SPLIT PAYLOAD AND TWO LEVERED CEILINGS, per object and per class?** ⛔ Its mass-weighted headline is the wrong yardstick for pass-0 — honest ignorance reads as error there |
| `design/worst_objects.py` | ⭐⭐ **WHICH REGIONS AND BOUNDARIES CARRY ONE CONDITION'S ERROR MASS?** Read the concentration curve first — concentrated means a mechanism exists, diffuse a systematic bias; `fg_loc` vs `pred_fg` separates a bad local solve from bad messages |
| `design/calibration_truth_ab.py` | ⭐⭐ **HOW DOES THE DELIVERABLE SCORE AGAINST TRUTH, AND WHAT IS PERFECTING EACH fl PMF WORTH (`--ceiling`)?** ⛔ Read the ceiling caution in Working rules first: the effective-length shrinkage sits outside every arm's patch point |
| **⭐⭐⭐ the panel, and the caches that make calibration a seconds-long loop** | |
| `sim/configs/flgap_rna_long.yaml` · `flgap_rna_short.yaml` | ⭐⭐ **WHAT BREAKS WHEN gDNA AND RNA FRAGMENT LENGTHS DIFFER? — the fl-gap SIDE panel, two arms of opposite sign** at `g50`. ⛔ Never a ladder rung: its transcript-level number is not a calibration result, though everything before the EM is valid ⛔⛔ **NOT REGENERATED on 2026-08-22: both arms still carry the RETIRED UNIFORM nascent model (`mode: fragment_share`), while the ladder carries SPARSE.** Each panel's on-disk data matches its own config, so both are internally accurate and their recorded measurements stand — but a ladder-vs-side-panel comparison now varies TWO things, so no claim may be carried across them. Re-simulating is an open owner decision |
| `sim/configs/gdna_ladder.yaml` | ⭐ **THE STAGE-B PANEL, AND THE ONLY PANEL THE TOOL IS RANKED ON** — 16 conditions (`g00/g05/g50/g98` × ss `0.50/0.99` × capture off/on), gDNA 0 → 98 % at a fixed 10 M total. ⛔⛔ **EQUAL FRAGMENT LENGTHS, and that is a forcing function**: the EM already reads the fl distribution, so a length gap lets it split the origins on LENGTH ALONE and mask calibration bugs (owner, 2026-08-14). ⛔⛔ **Every row carries SPARSE nascent RNA (2026-08-22): `on_fraction 0.50` of gene SPANS on, level logU(1, 100) INDEPENDENT of the mature level — 20.2 % of RNA fragments, the retired uniform model's total distributed sparsely.** 0.50 is a DEVELOPMENT STRESS level and not real data (`DESIGN.md` §0b, THE NASCENT SCOPE RULING); realistic is 0.10 ⇒ 4.2 %. `docs/TESTING.md` §0 |
| `sim/panel.py` | ⭐⭐⭐ **HOW IS A PANEL BUILT, SIMULATED, CACHED AND SCORED? — one command per stage**: `status` / `build` / `simulate` / `cache` / `score` / `report`, every path derived from one panel YAML. ⭐ Run `status` FIRST — every stage is expensive and resumable, and it names the next one. ⛔ It adds no measurement code: each stage shells out to the instrument that already owns it. ⛔ `cache` builds BOTH caches, and the oracle one is the origin-split truth every scoring instrument reads. Gated by `tests/test_panel_workflow.py` |
| `design/build_scan_cache.py` | **SCAN ONCE, CALIBRATE MANY TIMES.** ⛔ The cache key hashes `accumulator.cpp`'s deposit rule and not `resolve.cpp`'s fragment construction, so for a change to which fragments are OFFERED use `--force` or delete the caches by hand |
| `design/rescan_panels.py` | ⭐⭐⭐ **DID A RE-SCAN CHANGE ONLY WHAT IT WAS SUPPOSED TO? — the first irreversible step, so it carries its own falsification.** Gates each condition on byte-identity against the stale cache read before the write. `--self-test` 14/14 |
| `sim/build_suite_reference.py` · `design_suite_probes.py` · `simulate_reads.py` | **HOW IS THE PANEL'S SUBSTRATE BUILT?** ⚠ `panel.py build` drives the last two; the reference carve needs the source genome/GTF, which a panel config does not name, so it stays manual |
| **⭐⭐ the prior assembler, and the end-to-end thermometer above it** | |
| `design/quant_accuracy.py` | ⭐⭐⭐ **HOW ACCURATE IS THE TOOL END TO END, AND WHAT IS A PERFECT PRIOR WORTH?** `--arm base` plus the oracle and per-field injection arms, scored count against count. ⚠ A THERMOMETER above 0.8.0's metric, never the target |
| `design/mass_prior_ab.py` | ⭐⭐⭐ **CAN THE PRIOR BE A CONSERVED FRAGMENT COUNT RATHER THAN ONE MANUFACTURED FROM A DENSITY?** Subsamples by qname hash so the whole and all three origin partitions stay consistent. ⛔ The subsample must reproduce the defect first |
| `design/transcript_truth.py` | ⭐⭐⭐ **WHAT IS THE TRUE PER-TRANSCRIPT COUNT, SPLIT BY SPLICEDNESS?** One pass over the oracle BAM, read names only. ⛔ Splicedness comes from spliced-transcript coordinates, NEVER the CIGAR, which misses every sj in the unsequenced inner gap |
| `design/transcript_weights.py` | ⛔⛔ **WHAT IS THE SOFT-MIN PER-TRANSCRIPT WEIGHT WORTH? — BUILT, PRICED AND REFUSED; a RECORD, not a proposal** (`ROADMAP.md` §4). ⚠ Standalone it runs only its own re-partition falsification; scoring belongs to `quant_accuracy.py` |
| `rigel.sim.net_flow` (a MODULE, not a script) | ⭐⭐ **WHERE DID EACH MISASSIGNED FRAGMENT GO?** The DIRECTION of transcript error, per transcript, split into gDNA-sourced and RNA-isoform-sourced flow — the one question an accuracy table cannot answer. Gate: `tests/test_net_flow.py` |
| **⭐⭐⭐ where to develop** | |
| `design/rename_identity.py` | ⭐⭐⭐ **IS THIS RENAME STAGE NUMERICALLY A NO-OP?** `--freeze` captures one reference, `--check` compares after every stage — on array CONTENT and the transcript table, never on names. ⚠ The reference is frozen, never rolling. `--self-test` 8/8 |
| `design/rename_census.py` | ⭐⭐⭐ **WHICH NAMES DOES A VOCABULARY RULING TOUCH, AND WHICH CARRY TWO SENSES?** Reports by kind — identifiers, C++, prose — and never renames; `--sense <token>` dumps every site with context. ⛔ Run it before renaming anything |
| `design/module_census.py` | ⭐⭐⭐ **WHERE DOES A CHANGE GO?** The calibration package re-derived from the AST: the layering with every upward import, each module's importers, docstrings naming a sibling with no import, dead public surface. ⛔ It reports; it does not judge |
| **⭐⭐⭐ the backbone** | |
| `design/arm_identity.py` | ⭐⭐⭐ **IS THIS ARM BYTE-IDENTICAL TO THAT ONE?** Compares every scored field of every row, where an aggregate hides a difference that cancels between two fields; the row-key sets must be EQUAL. ⛔ Falsified by a 1-ULP nudge |
| `design/backbone_parity.py` | ⭐⭐⭐ **WHAT DOES ONE MESSAGE OPERATOR DO, PER SLOT?** Two policies on one real chain in one process, every output array and diagnostic key compared element by element. ⭐ Strictly stronger than the panel per condition, so run it first |
| **⭐⭐ the toy harness** | |
| `design/toy_panel.py` | ⭐⭐ **HOW DOES ONE TOY SPEC BEHAVE ACROSS EVERY CACHED CONDITION AND AN RNA-DENSITY LADDER, scored per object?** It names which object carries the error and whether the messages helped it. ⚠ 13 s per condition — shard with `--conditions` |
| `design/verify_toy_substrate.py` | ⭐⭐⭐ **IS THE INPUT CORRECT? — no solver runs.** Every accumulator bank re-derived from per-fragment truth by an independent implementation, plus the splice combinatorics and the length marginal. ⛔ Run it on any new toy spec first |
| `design/verify_capture.py` | ⭐⭐ **WHAT DOES HYBRID CAPTURE DO ON IDENTICAL GEOMETRY, probes ON vs OFF?** The gDNA landscape, the length selection and the sj depletion, each gated on the direction the knobs predict |
| `design/toy_trace_error.py` | ⭐⭐⭐ **WHERE DOES ONE SCENARIO'S ERROR COME FROM, END TO END?** Every object ranked by error mass, the four init sources, the relay hop by hop with its licence state, then a ψ channel ablation |
| `design/zero_controls.py` | ⭐⭐⭐ **DOES THE TOOL HOLD AT ZERO RNA AND AT ZERO gDNA? — the owner requires both on every experiment.** The truth is a constant, so every deviation is a false positive. ⛔ Flags any EMPTY object: a degenerate zero arm tests nothing |
| `design/reframe_walk.py` | ⭐⭐⭐ **WHAT DOES THE REFRAME DO AT EVERY HOP, IN BOTH DIRECTIONS, on one two-exon transcript?** `r` as used, `r` as the predecessor's single total would give it, and what truth says — so each decision is a number |
| `design/toy_dissect.py` | ⭐⭐ **WHICH CHANNEL IS ALIVE AT EACH SLOT OF ONE SCENARIO?** The solver's own ladder beside `cm_g` / `c_tau` / `cg`, so a dead channel is visible as a zero rather than inferred |
| `design/certified_rna_audit.py` | ⭐⭐⭐ **IS THE CERTIFIED-RNA CHANNEL WIRED?** Audits whether the bank is populated, whether it has a divisor and whether a precision is emitted — any one failing looks identical. ⛔ Scores the MASS, never `f_g` |
| `design/certified_q_census.py` | ⭐⭐ **CAN A CERTIFIED COUNT SPEAK ABOUT THE UNSPLICED SPLIT? — the answer is NO, and this measures why**, straight off the origin-split oracle with no solver. ⭐ Its two extreme rungs are the two zero controls |
| `design/vertex_ceiling.py` | ⭐⭐ **WHAT IS PINNING ORACLE TRUTH AT ONE OBJECT CLASS WORTH ON THE REAL LADDER?** A `noop` arm must be byte-identical, and `--arm ref_c=A,B` drives ψ's two Beta reference exponents. ⛔ It prices missing information, not headroom |
| `design/toy_ceiling.py` | ⭐⭐ **WHAT IS A DIFFERENT OWN BELIEF AT ONE OBJECT CLASS WORTH, re-solving the whole chain?** Six arms share one simulation, and `--arms base noop` is its own falsification — those two must be byte-identical |
| `design/psi_channel_ablation.py` | ⭐⭐⭐ **WHICH ψ CHANNEL PUT THE ERROR THERE?** Every argument of the final ψ combine recorded, then re-solved with one imputed channel nulled at a time. ⛔ The as-is arm must reproduce the run bit-identically, write-back included |
| `design/arm_score.py` · `design/arm_sweep.py` | ⭐⭐ **HOW DO TWO `ladder_arm_ab` ARMS COMPARE, AND HOW DOES A WHOLE FAMILY SWEEP? — per stratum.** ⛔ Never pooled: the panel total hides a sign flip between strata. Both print the deliverable's error beside pass-0's |
| `design/ladder_arm_ab.py` | ⭐⭐ **DOES THE SAME OVERRIDE STILL HELP ON THE REAL LADDER?** ⛔⛔ Run it before writing a mechanism into `src/` — four times a toy-positive change has been panel-negative. `--messages {off,on}` is part of the arm. `--self-test` 25/25 |
| `design/length_ceiling.py` | ⭐ **WHAT IS A PERFECT LENGTH MODEL WORTH ON THE LADDER, ONE fl PMF AT A TIME?** ⚠ This is the OPPORTUNITY model's fl PMF, not the length-likelihood composition channel deferred past 0.8.0 |
| `design/toy_harness.py` | ⭐⭐ **HOW DOES A MINI CHROMOSOME YOU DEFINE CALIBRATE — in 0.1–5 s, with every object's answer beside its truth?** (`docs/TESTING.md` §0b) The priors a toy cannot fit are harvested from a real cached condition; `--list` for the ladder |
| **the substrate — are the panel and the index sound?** | |
| `design/simulator_gates.py` | **DOES THE SIMULATOR PASS ITS OWN GATES, scored on per-fragment truth?** ⛔ Run it before trusting the panel |
| `design/suite_resolves.py` | **CAN THE SUITE RESOLVE THE AXIS YOU ARE CHANGING?** ⛔ Run it before quoting any suite number |
| `design/index_census.py` | **WHAT IS ACTUALLY IN THIS INDEX?** ⛔ Re-derive the census; never quote a stored table |
| `design/verify_index_rebuild.py` | **DID AN INDEX REBUILD PRESERVE THE STRUCTURE?** Regions byte-identical, boundaries only in contiguous reach |
| **Stage A — the accumulator** | |
| `design/fl_anchor_gap.py` | **HOW DO THE ANCHOR AND BOTH LENGTH MODELS COMPARE WITH TRUTH?** `--drain` measures before and after the second pass |
| `design/gdna_pool_census.py` | ⭐ **DOES EACH OF THE FOUR gDNA POOLS AGREE WITH ITS OWN OPPORTUNITY, AND WITH TRUTH?** |
| `design/second_pass_accuracy.py` | **HOW ACCURATE IS THE SECOND PASS, PER FRAGMENT, against the oracle BAM's read names?** |
| `design/observable_efficiency.py` | **WHAT FRACTION OF THE LENGTH INFORMATION DOES A STORAGE CHOICE KEEP?** |
| `design/region_density_derivation.py` | **DOES THE RECIPROCAL-OPPORTUNITY THEOREM HOLD?** T0–T6, each perturbed |
| **diagnostics** | |
| `design/anchor_opportunity_census.py` | ⭐⭐ **IS A ZERO-COUNT ANCHOR'S DENSITY CLAIM TRUE OF ITS NEIGHBOURHOOD? — no solver runs.** The empty pure-gDNA population against what its own chain neighbours measure. ⚠ Its population is `g1_locked ∧ REGION`, not `struct_lock` |
| `design/composition_evidence_census.py` | **HOW MUCH LIBRARY MASS REACHES THE SOLVER WITH NO COMPOSITION EVIDENCE?** `--inject-kappa 0.5` is its falsification handle |
| `design/held_flux_census.py` | **HOW OFTEN DOES A HELD CANDIDATE HAVE ZERO FLUX EVIDENCE, AND BY WHAT CAUSE?** |
| `design/prior_units_check.py` | **IS THE EM PRIOR IN FRAGMENT UNITS, OR STILL THE OLD INCIDENCE SUM?** |
| `design/boundary_q_population.py` | ⭐⭐⭐ **IS THE PRIOR'S CROSSING→FRAGMENT CONVERSION POPULATION-BLIND?** Feeds a perfect `f_g` from the oracle so the `q` conversion is isolated, and scores it against the conserved truth. ⛔ Real but bounded: record the bound, build nothing |
| **plumbing** | |
| `design/native_parity_on_real_data.py` | **DOES NATIVE PARITY HOLD ON REAL cfRNA AT FULL SCALE?** |
| `design/scan_profile.py` | **HOW MANY ns PER FRAGMENT, regressed over several BAMs?** ⚠ A `profiling/` file shares this basename and is a different instrument |

## CLI

```bash
rigel index --fasta genome.fa --gtf annotation.gtf -o index/
rigel quant --bam sample.bam --index index/ -o results/
rigel sim --config scenario.yaml -o out/
rigel export results/ -f tsv
rigel report results/ -o report.html
```

Input BAM must be name-sorted with the `NH` tag.
