# CLAUDE.md

Directions for a session working in this repository. It says what Rigel is, what the current release is
about, where each kind of change goes, how to run the benchmarks and the suite, and which instrument
answers which question. It is not a history: the changelog is git, rulings live in `docs/DESIGN.md`, open
problems in `docs/ISSUES.md`, lessons in `docs/TRAPS.md`. A rule earns its place here by changing what the
next session does.

## What Rigel is

A Bayesian RNA-seq transcript quantifier that separates RNA from genomic-DNA contamination. A single-pass
C++ BAM scanner tallies fragments, a **calibration** stage deconvolves the library into gDNA vs RNA, and a
per-locus EM solver assigns RNA to transcripts. PyPI package `rigel-rnaseq`; the import and CLI are
`rigel`. The version on disk is 0.7.1 (`pyproject.toml`); the target is **0.8.0, a calibration release**.

## Axiom 0 — RNA is RNA

There are three populations and no fourth: `gDNA`, `RNA+`, `RNA−`. "Mature" and "nascent" are not
populations and not a degree of freedom: RNA inside an intron is RNA that has not spliced at that
position. The only distinction is whether a fragment **is spliced** (certified RNA, needing no
deconvolution) or **is not** (the whole deconvolution problem). The population set at a slot is a function
of two bits, which is what makes this structural:

```
T(slot) = {gDNA}                                  # always — gDNA is genomically continuous
        ∪ {RNA+ if statics.free_pos[slot]}        # iff the annotation admits that strand here
        ∪ {RNA− if statics.free_neg[slot]}
```

so `|T| ∈ {1,2,3}`. The tell that you have violated this: a population set with more than three members,
or the words "mature", "nascent" or "a third component" in a solver or composition question. Re-ask it as
"what is this channel's OPPORTUNITY for RNA at this object?" — the answer is a geometry, derivable from
the index. The words survive only as simulator inputs (`nrna_abundance`, the toy harness's `--nrna`).
How much weight a nascent concern carries is the nascent scope ruling's question (`docs/DESIGN.md` §0b):
nascent RNA is sparse in real data and is modelled for robustness, so a number measured at the panel's
nascent share is a stress reading, never a design driver.

## The 0.8.0 scope

The focus is CALIBRATION, and the metric is the calibration result scored against oracle calibration
(`calibration_vs_oracle.py`, `solvability_audit.py`, `prior_vs_oracle.py`), never the end-to-end
transcript number, which stays a thermometer (`docs/SUCCESS.md`). Three strata are in scope — unstranded ×
capture-OFF, stranded × capture-OFF, stranded × capture-ON — and **unstranded × capture-ON is DEFERRED**:
still reported on every benchmark, never a development target until the other three are optimised. It is
also where most of the error is (the gDNA fraction cancels from the strand mean, so an unstranded AMBIG
slot has no channel), so the debug loop must take the worst IN-SCOPE scenario, never the deferred one.

The fragment-length COMPOSITION channel is retired until after 0.8.0 and may not be proposed; it does not
exist in `src/`. Three other things called "length" are unaffected: layer 2's `fl` / `effective_length` /
`capture_eff_length` (the opportunity model), `length_likelihood` in `second_pass.py` (per-fragment
assignment), and the fl PMFs priced by `em_fl_ceiling.py`. The ladder gives gDNA and RNA EQUAL fragment
lengths on purpose: the EM already reads the fl distribution, so a gap would let it split origins on length
alone and mask calibration bugs. Scenarios are cached (`panel.py cache`) so calibration re-runs in seconds.
The full ruling is `docs/DESIGN.md` §0b; the ranked next steps are `docs/ROADMAP.md`.

## The docs

Nine permanent docs, none of them a changelog.

| doc | what it is |
|---|---|
| `docs/SUCCESS.md` | **Start here.** How performance is measured: Stage A (the accumulator) and Stage B (calibration against an oracle), with the instruments in run order |
| `docs/ROADMAP.md` | The short ranked view: the frame, one line per state claim naming its instrument, the ordered next steps. No numbers, no history |
| `docs/ISSUES.md` | The issue log: every open problem as a named entry (`ISSUES: kebab-name`) with priority and instrument, plus the append-only CLOSED / REFUSED record that keeps each refusal's killing number |
| `docs/TRAPS.md` | Mistakes already made, as rules that change what the next session does. Cite by name: `TRAPS: kebab-name` |
| `docs/EQUATIONS.md` | The derivations the code depends on, each named from the module that implements it |
| `docs/DESIGN.md` | What is built and the rulings behind it — settled, not re-litigated. §0 is the binding vocabulary; §0b the 0.8.0 scope and the nascent scope ruling |
| `docs/TESTING.md` | A manual: how to build each panel and reference, how to run each gate, what the suite can and cannot judge. §0a is the test chromosome, §0b the toy harness |
| `docs/MANUAL.md`, `docs/PUBLISHING.md` | The user's manual and the release procedure |

**The move rule.** When a finding settles, MOVE it to its one home (an open problem or refusal to
`ISSUES.md`, a lesson to `TRAPS.md` as a named rule, a ruling to `DESIGN.md`, a derivation to
`EQUATIONS.md`) and delete it where it was, in the same edit. Never copy: two homes diverge.

**`docs/dev/` is the sandbox** — working notes, half-finished arguments, handoffs. Nothing there is
authoritative and nothing may cite into it (`tests/test_docs_boundary.py`). **The source does not cite
the docs**: a docstring may cite a test or the executable specification
(`tests/native/_accumulator_reference.py`), never a doc.

## Where does a change go? — the calibration layering

An import may point DOWN a layer or SIDEWAYS within one, never UP. A module reaching for something a
layer up is telling you the thing belongs lower. `rigel/calibration/_layers.py` is authoritative,
`tests/calibration/test_layering.py` enforces it, and `python scripts/design/module_census.py` re-derives
the graph from the AST.

| if the change is about… | it goes in |
|---|---|
| what a fragment tally MEANS | **1 · the payload view** — `splice_graph` `substrate` `region_arrays` |
| how many places a fragment COULD have sat | **2 · opportunity** — `effective_length` `capture_eff_length` `sj_opportunity` `gdna_opportunity` `fl` |
| one slot's own numbers, ψ, and its total | **3 · geometry + the per-slot solve** — `region_geometry` `simplex_logodds` `total_abundance` |
| which strand a fragment came from | **4 · strand** — `gdna_strand` `strand_balance` `strand_summary`, and `strand_likelihood` (a gated executable reference) |
| how dense a component is, and the priors | **5 · density and prior** — `density_model` `density_deconv` `landscape` `abundance_landscape` |
| what one neighbour tells another | **6 · the solve** — `sweep` (the backbone) + `blocks` (one locus block cut out and put back) + `message_cache` + `messages/` (the policy) + `region_init` |
| turning the solve into a result | **7 · assemble** — `calibrate` `priors` `result` `derive` `diagnostics` `track` |

## The message layer

`CalibrationConfig.message_policy = "transfer"` ships (the default since 2026-09-09, `message_propagation
= True`). Two policies, selected by one config value (an unknown name raises), both on the two-phase
backbone (`docs/DESIGN.md` §6b.11–§6b.12): `prepare` (every node's own claim) → a forward pass and a
backward pass of `receive(source, destination)`, so every node ends with one message from each
neighbour it has → `solve(from_left, from_right)`, which hands ψ two row channels and nothing else.

| policy | |
|---|---|
| `silent` | the measured floor (`messages/silent.py`); the same policy `message_propagation = False` installs |
| `transfer` | the shipped default (`messages/transfer.py`; the face table in `messages/faces.py`, the level lanes in `messages/lanes.py`, pure row constructors in `messages/transfer_rows.py`). `prepare` is a table of contents, one named builder per message: `_claims`, `_splice_faces`, `_edge_level`, `_terminus_rules`, `_alternative_splice_site`, `lanes.gdna_lane`, `lanes.rna_lanes`. Every hop pays its pair's counting plus the disagreement beyond it |

Messages exist for the slots whose own solve has no composition channel — unstranded data and AMBIG
slots. **We do not expect to beat `silent`**: on strand-specific data a sighted exon's own solve is
excellent. The goal is to win on unstranded data while doing minimal harm on stranded data, and the two
halves are judged against different bars and never pooled (`policy_benchmark.py` prints them apart). The
certified flux is a message: spliced fragments are measured at boundaries, never solved, strictly one hop
(`docs/DESIGN.md` §6b.13). Two earlier policy campaigns were torn down with their numbers (`ISSUES:
the-message-policy-campaign`) — read them before re-proposing a mechanism.

## Running the benchmarks

Two substrates. Develop on the test chromosome; decide on the ladder.

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
python scripts/design/preflight.py                          # first: can this session run at all?

python scripts/design/policy_benchmark.py --panel test      # the development loop: 30 conditions, seconds
python scripts/design/policy_benchmark.py --panel ladder    # the shipping judgement: 16 conditions, minutes
python scripts/design/policy_benchmark.py --panel ladder --policies silent transfer --by-class
```

The test chromosome is one hand-edited YAML, `scripts/sim/test_reference/test_chr.yaml`; its GTFs,
abundances, probe panels and FASTA are rendered from it by `build_test_reference.py`, a suite gate refuses
a drifted render, and after editing it everything derived must be rebuilt (`docs/TESTING.md` §0a has the
recipe; `panel.py status` names the next stage). Read the two halves separately and never pool them:
unstranded rows are where a policy must win, stranded rows where it must do minimal harm. A toy and the
panel have inverted a ranking before (`TRAPS: a-toy-and-a-panel-can-disagree-in-rank`): confirm on the
ladder.

## Cite a rule by its name

Cite a trap as `TRAPS: off-grid-message-mode`, an issue as `ISSUES: two-sided-exon-row`, never by a
number. `tests/test_no_jargon_labels.py` enforces it: the old numbered labels were ambiguous (`G1` meant
both a process rule and a structurally pure-gDNA object).

## Working rules

- **DERIVE → DESIGN → PLAN → PROTOTYPE → A/B → only then `src/`.** No idea enters the production source
  before it has been derived on paper, prototyped outside `src/`, and A/B'd against the policies that
  already exist. One mechanism at a time: a change that cannot be A/B'd alone cannot be judged alone.
- **No magic numbers.** Stop and discuss before adding any constant, heuristic or tunable. Every divisor
  must be derived from the deposit rule and unit-tested against brute-force enumeration.
- **A falsification test first, verified failing — then break the fixed code and watch each gate fire.**
  The second half is not optional; it has found holes in already-green gates repeatedly.
- **The debug loop is the default method**: run the panel → take the worst IN-SCOPE scenario → dissect it
  to the highest-error objects (`worst_objects.py`, `calibration_walk.py`) → find the mechanism → fix →
  re-run the panel.
- **A ceiling is sometimes the right instrument** (`calibration_truth_ab.py --ceiling`) but prices
  something that may be unreachable, so it is not the default. Every ceiling arm patches
  `assemble_priors`, while `pipeline.py` builds `effective_lengths_em` before calling it, so the
  effective-length shrinkage has never been inside any ceiling number.
- **One thing varied per experiment**, a baseline re-recorded from the current tree in the same session,
  and **score against truth** (the oracle BAM's read names), never against the previous run.
- **No legacy, no backwards compatibility, no speculative code.** Converge and delete. No version
  suffixes in file names. No Greek letters in identifiers (fine in maths write-ups).
- **Real data is a test input, never a design input.** Sweep the plausible space, report the worst case.
  Profile on a deep real RNA-seq library, never a panel condition or a toy; the cfRNA libraries are smoke
  tests (`docs/TESTING.md` §7). Read a timing only from back-to-back A/B pairs — untouched stages drift
  25–40 % between runs taken at different times — and prove every speed-up a numeric no-op
  (`rename_identity.py --bam`, `profiling/sweep_replay.py`).
- **Renames**: run `rename_census.py --sense <token>` before renaming anything and `rename_identity.py
  --check` after each stage. `arm` still carries three senses (an experiment arm, a component arm, and
  `__ARM_NEON` in the scanner) and needs its own `--sense` pass, never a tail-end sweep.
- **The owner drives commits.** Do not commit unless asked.

## Build, test, lint

Every build, test and lint command runs inside the activated `rigel` conda environment — it holds htslib
and the compilers, and the C++ build finds htslib via `$CONDA_PREFIX`.

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel

pip install --no-build-isolation -e ".[dev]"   # rebuild after ANY src/rigel/native/ change
python -m pytest tests/ -q                     # never bare `pytest` — the repo root leaves sys.path
python -m pytest tests/ --update-golden        # regenerate tests/golden/ after intended output changes
ruff check src/ tests/ scripts/ && ruff format src/ tests/   # never format scripts/
```

**The standing baseline: 0 failed / 3,434 passed / 0 skipped / 2 xfail, 3,436 collected** (re-derived
2026-09-12 after the cleanup split, the received tables, the replay's tolerance report — one `tests/` file at +2
holding two gates, +4 — and the one ψ solver, which retired the float32 hoisting gate, −1: four `src/rigel/calibration/` modules added —
`blocks`, `message_cache`, `messages/faces`, `messages/lanes` — at +3 each, +12; one backbone test retired
with its premise (a kernel can no longer leave a hop unspoken) and one gate added on the face table, ±0). The 2 xfails are executable records of proven defects whose fixes are elsewhere
(`ISSUES: two-sided-exon-row`; the antisense prior-assembly casualty) — "fix the test" is a category
error, and an xfail is closed by repairing the thing or asserting the invariant structurally, never by
widening a bound. **Any failure at all is a regression.** A commit that measures the suite updates this
line.

**Re-derive a count, never adjust one** (`TRAPS: re-record-the-baseline`). Several gates are parametrised
over the files on disk, so adding or retiring a file moves the total; account for it from this table and
confirm with `pytest --collect-only -q | grep <stem>`:

| adding one… | moves collected by | which cases |
|---|---|---|
| `src/rigel/calibration/` module | **+3** | jargon, docs-boundary, and layering *if declared in `_layers.py`* |
| `tests/` file (any directory) | **+2** | jargon, docs-boundary |
| `scripts/design/` (or `sim/`, `profiling/`) file | **+4** | imports, says-what-it-is-for, jargon, docs-boundary |
| `docs/dev/` file | **+1** | jargon only |
| top-level `docs/` .md (an owner decision — the permanent-set gate pins the list) | **+2** | jargon, docs-boundary |

A content-only sweep moves the collected total by zero, and that is the check. Derive the failure set,
never eyeball the tail (`TRAPS: read-the-whole-failure-list`):

    python -m pytest tests/ -q 2>&1 | grep '^FAILED' | sed 's/::.*//' | sort | uniq -c

A golden update is where a regression gets laundered into "intended": read the diff and record its
magnitude before `--update-golden`, and check the truth-scored instruments first. After a `src/`
deletion, a rename, or a shipped-default flip, run the instruments and not only the suite —
`preflight.py --full` does it in one command (`TRAPS: a-green-suite-hid-five-dead-instruments`).
Always set `OMP_NUM_THREADS=1` when benchmarking or comparing runs.

## Tooling under `scripts/`

This table indexes `scripts/design/` plus four `sim/` rows; `tests/test_scripts_index.py` holds it
against the disk in both directions. `scripts/profiling/` is indexed in `scripts/README.md`. Each row carries only the
question its instrument answers; `docs/SUCCESS.md` has the run order.

| | |
|---|---|
| **⭐⭐⭐ START A SESSION HERE** | |
| `design/preflight.py` | ⭐⭐⭐ **CAN THIS SESSION RUN AND REGENERATE EVERYTHING? — one command, one verdict, before anything else.** Checks the toolchain (the `rigel` env, the native extension, the CLI), both references, both panels (scan caches, oracle caches with all five partitions, the certified `slot_truth`) and that every `scripts/design/` instrument IMPORTS. ⛔ It changes nothing and measures nothing — every check is a read or an import, and a ✘ prints the exact command that regenerates the missing artifact. ⭐ **The default is ~2 s**; `--full` adds every instrument's `--self-test` (minutes since the relay-era arm harness retired on 2026-09-09) — run it after a deposit-rule change or a default flip, never every session. `--self-test` 8/8 |
| **⭐⭐⭐ THE POLICY BENCHMARK — where a message-policy change is judged** | |
| `design/policy_prototype.py` | ⭐⭐⭐ **HOW DOES A PROTOTYPE MESSAGE POLICY SCORE, PER GENE TYPE AND PER SLOT, AGAINST CERTIFIED TRUTH?** — the harness every message rung is developed on before `src/`. Installs a class from `--module` in place of the shipped policy for the `transfer` arm; whole-library and per-type tables, `--by-class` (the error at each NODE CLASS — the view that judges a message at its destinations), `dissect` for one gene type slot by slot. ⛔ Compare src-vs-src across a landing (`TRAPS: a-harness-on-the-parent-class-dies-when-the-parent-gains-the-mechanism`). `--self-test` |
| `design/policy_benchmark.py` | ⭐⭐⭐ **HOW DOES EACH POLICY SCORE, PER CONDITION, AGAINST CERTIFIED TRUTH?** Whole-library gDNA error in fragments, per axis, one row per condition, for `silent` and `transfer`. ⭐ `--panel test` is the test chromosome (seconds — the development loop); `--panel ladder` is the 16-condition benchmark. ⭐ **`--by-class`: WHERE DOES A POLICY'S REMAINING ERROR SIT, BY NODE CLASS?** — per certified stratum, boundaries split by terminus flag, exons by reach (licensed face / edge only / walled); the instrument that ranks the rebuild's holes. ⛔ NEVER POOLED, and the two halves are judged against DIFFERENT bars: unstranded rows are where a policy must WIN, stranded rows are where it must do minimal HARM against silence |
| **⭐⭐⭐ 0.8.0'S METRIC — calibration against ORACLE CALIBRATION** | |
| `design/calibration_vs_oracle.py` | ⭐⭐⭐ **IS THE CALIBRATION RESULT ITSELF RIGHT, SCORED AGAINST AN ORACLE CALIBRATION? — 0.8.0's metric, and the only instrument that reaches the effective-length shrinkage.** `P = calibrate(...)` against the same payload with only the six deconvolved arrays swapped, per stratum, plus `U`, the no-enrichment null no other instrument carries; `--message-policy` prices a policy on this metric (the ship protocol's first item). ⛔ Read `ruler_n_moved`, never the aggregate: the total can barely move while nearly every transcript is redistributed. No solver, no EM, no re-scan — ~5–12 s/condition. `--self-test` 21/21 |
| `design/object_composition.py` | ⭐⭐⭐ **MUST ψ's Beta REFERENCE BE ONE LIBRARY-WIDE NUMBER, OR CAN EACH OBJECT SUPPLY ITS OWN?** `m_i` per object from the two densities, scored as misplaced fragments against the shipped ½, per stratum. `--self-test` 25/25 |
| `design/abundance_landscape_census.py` | ⭐⭐⭐ **WHAT DOES THE TOTAL-DENSITY FIELD LOOK LIKE, PER CONDITION?** Fits `calibration.abundance_landscape` on the cached wall-exact totals — every mode's basin mass, `rho_0`, the anchor gap in nats, the per-class enrichment. `--self-test` 13/13 |
| `design/calibration_oracle.py` | ⭐⭐⭐ **WHAT IS THE CERTIFIED PER-OBJECT TRUTH? — run this before debugging calibration against anything.** Every REGION and BOUNDARY's count, its realized `n_gdna`/`n_nrna`/`n_mrna` and `true_f_g`, at two certification levels: COMPOSITION (no opportunity model anywhere in it) and FIELD (densities too). ⛔ REFUSED unless its named gates pass — sum-to-full, partition-projects-exactly, gdna-field-uniformity, exact-zeros, nascent-in-annotation — because a merely plausible oracle is how a calibration bug and a truth bug survive each other. Writes `slot_truth.npz` beside each oracle cache; `--self-test` 11/11 |
| `design/total_abundance_audit.py` | ⭐⭐⭐ **IS THE MEASURED TOTAL A TRUE TOTAL?** Five arms against the origin partitions; read ⓔ START/END agreement first — the only field-free arm and the decisive test of the wall rule. `--self-test` 15/15 |
| `design/landscape_training_census.py` | ⭐⭐⭐ **WHICH SLOTS TRAIN THE gDNA LANDSCAPE PRIOR, WITH WHAT EVIDENCE, AND HOW MUCH OF THAT TRAINING IS FALSE?** Spies each refit's `fit_landscape` inputs and each sweep's held messages; per refit, per node class and per evidence class (anchor · locked · own:strand · own:factory · delivered:composition · delivered:bound · none), the weight the estimator summed, the gDNA trained and its share on certified-zero slots, with both zero controls from the final answer; `--estimator` re-fits the last population at the certified values. ⛔ Reads `Σw`, never the trained mass, on a `g00` row. `--self-test` 17/17 |
| `design/calibration_walk.py` | ⭐⭐⭐ **WHICH STAGE OF CALIBRATION INTRODUCES THE ERROR?** The solve as a ladder — init → strand → local → +messages → +refits → shipped — each rung scored per stratum against `calibration_oracle.py`, which it refuses to run without |
| `design/structural_claims_audit.py` | ⭐⭐⭐ **IS EVERY SLOT THE STAGE-0 SUBSTRATE ADMITS TRULY WHAT IT CLAIMS? — the confusion matrix against certified slot truth, no solver.** Each structural class scored on ITS OWN claim in fragments; the solvable-exon claim is tested at the licensing FLANK, and nascent inside an ss intron is not a violation. ⛔ REFUSED without `slot_truth.npz`. `--self-test` 8/8 |
| `design/transport_dispersion.py` | ⭐⭐⭐ **WHERE DOES THE FLANK-TRANSPORT DISPERSION COME FROM? — the decomposition against certified truth, no solver.** Pair disagreement vs common-mode center, each charged with counting (flank AND truth side), the length curve, structure and capture. ⛔ Fit nothing on shallow pairs; the truth count's own trigamma must be subtracted before quoting any certified scatter |
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
| `rigel.sim.net_flow` (a MODULE, not a script) | ⭐⭐ **WHERE DID EACH MISASSIGNED FRAGMENT GO?** The DIRECTION of transcript error, per transcript, split into gDNA-sourced and RNA-isoform-sourced flow — the one question an accuracy table cannot answer. Gate: `tests/test_net_flow.py` |
| **⭐⭐⭐ where to develop** | |
| `design/rename_identity.py` | ⭐⭐⭐ **IS THIS RENAME, REFACTOR OR SPEED-UP NUMERICALLY A NO-OP?** `--freeze` captures one reference, `--check` compares after every stage — on array CONTENT and the transcript table, never on names; `--bam` takes a real library instead of a panel condition. ⚠ The reference is frozen, never rolling. `--self-test` 8/8 |
| `design/rename_census.py` | ⭐⭐⭐ **WHICH NAMES DOES A VOCABULARY RULING TOUCH, AND WHICH CARRY TWO SENSES?** Reports by kind — identifiers, C++, prose — and never renames; `--sense <token>` dumps every site with context. ⛔ Run it before renaming anything |
| `design/module_census.py` | ⭐⭐⭐ **WHERE DOES A CHANGE GO?** The calibration package re-derived from the AST: the layering with every upward import, each module's importers, docstrings naming a sibling with no import, dead public surface. ⛔ It reports; it does not judge |
| **⭐⭐⭐ the backbone** | |
| `design/arm_identity.py` | ⭐⭐⭐ **IS THIS ARM BYTE-IDENTICAL TO THAT ONE?** Compares every scored field of every row, where an aggregate hides a difference that cancels between two fields; the row-key sets must be EQUAL. ⛔ Falsified by a 1-ULP nudge |
| `design/backbone_parity.py` | ⭐⭐⭐ **WHAT DOES ONE MESSAGE OPERATOR DO, PER SLOT?** Two policies on one real chain in one process, every output array and diagnostic key compared element by element. ⭐ Strictly stronger than the panel per condition, so run it first |
| **⭐⭐ the toy harness** | |
| `design/toy_panel.py` | ⭐⭐ **HOW DOES ONE TOY SPEC BEHAVE ACROSS EVERY CACHED CONDITION AND AN RNA-DENSITY LADDER, scored per object?** It names which object carries the error and whether the messages helped it. ⚠ 13 s per condition — shard with `--conditions` |
| `design/verify_toy_substrate.py` | ⭐⭐⭐ **IS THE INPUT CORRECT? — no solver runs.** Every accumulator bank re-derived from per-fragment truth by an independent implementation, plus the splice combinatorics and the length marginal. ⛔ Run it on any new toy spec first |
| `design/verify_capture.py` | ⭐⭐ **WHAT DOES HYBRID CAPTURE DO ON IDENTICAL GEOMETRY, probes ON vs OFF?** The gDNA landscape, the length selection and the sj depletion, each gated on the direction the knobs predict |
| `design/zero_controls.py` | ⭐⭐⭐ **DOES THE TOOL HOLD AT ZERO RNA AND AT ZERO gDNA? — the owner requires both on every experiment.** The truth is a constant, so every deviation is a false positive. ⛔ Flags any EMPTY object: a degenerate zero arm tests nothing |
| `design/vertex_ceiling.py` | ⭐⭐ **WHAT IS KNOWING THE TRUTH AT THE PARAMETER-VERTEX OBJECTS WORTH ON THE REAL LADDER?** The pin is the node's own claim plus its ψ row (re-pointed to the two-phase solve, 2026-09-09); the population is the PARAMETER vertex (silent genes, nascent-free introns, the zero rows — never the realized vertex, which priced chance), so it needs `--oracle-cache`. A `noop` arm must be byte-identical, and `--arm ref_c=A,B` drives ψ's two Beta reference exponents. ⛔ It prices missing information, not headroom |
| `design/toy_harness.py` | ⭐⭐ **HOW DOES A MINI CHROMOSOME YOU DEFINE CALIBRATE — in 0.1–5 s, with every object's answer beside its truth?** (`docs/TESTING.md` §0b) The priors a toy cannot fit are harvested from a real cached condition; `--list` for the ladder |
| **the substrate — are the panel and the index sound?** | |
| `design/simulator_gates.py` | **DOES THE SIMULATOR PASS ITS OWN GATES, scored on per-fragment truth?** ⛔ Run it before trusting the panel |
| `design/suite_resolves.py` | **CAN THE SUITE RESOLVE THE AXIS YOU ARE CHANGING?** ⛔ Run it before quoting any suite number |
| `design/index_census.py` | **WHAT IS ACTUALLY IN THIS INDEX?** ⛔ Re-derive the census; never quote a stored table |
| `design/verify_index_rebuild.py` | **DID AN INDEX REBUILD PRESERVE THE STRUCTURE?** Regions byte-identical, boundaries only in contiguous reach |
| **Stage A — the accumulator** | |
| `design/fl_pool_purity.py` | ⭐⭐⭐ **ARE THE FOUR gDNA LENGTH POOLS ACTUALLY PURE gDNA, AND WHAT DOES THE SHIPPED LENGTH MODEL SAY AGAINST TRUTH?** Per pool: the gDNA / nascent / mature counts and each component's mean length; then `TRUE` / `POOLED` / `SHIPPED`, so **contamination (`pool−true`) and the divisor+shrinkage (`ship−pool`) are attributed APART**. ⛔⛔ **Run it only where the two components' fragment lengths DIFFER** — the bias is `RNA_share × length gap`, and the ladder and test chromosome give them EQUAL lengths by design, so a 95 %-contaminated pool reads under a bp there. That is why the defect shipped |
| `design/em_fl_ceiling.py` | ⭐⭐⭐ **WHAT IS A PERFECT gDNA fl pmf WORTH END TO END, THROUGH THE EM? — the one fl question that stops at no earlier stage.** Every other fl instrument stops at `calibrate`, but `pipeline.py` also hands `gdna_pmf` to the fragment scorer, so a wrong length model is applied per fragment in the channel that separates origins. ⭐ Read `gdna_frac_est` against `gdna_frac_true` — the PRODUCT; ⛔ the transcript rows flip sign between the two fl-gap arms and are not the deliverable. Three gates: the injection counts its fires, `noop_fl` must be byte-identical, and `base_reseed` is the noise floor. ⛔⛔ Meaningless on an equal-length panel — run both sign arms AND the equal-length control |
| `design/gdna_pool_census.py` | ⭐ **DOES EACH OF THE FOUR gDNA POOLS AGREE WITH ITS OWN OPPORTUNITY, AND WITH TRUTH?** |
| **diagnostics** | |
| `design/prior_units_check.py` | **IS THE EM PRIOR IN FRAGMENT UNITS, OR STILL THE OLD INCIDENCE SUM?** |
| **plumbing** | |
| `design/native_parity_on_real_data.py` | **DOES NATIVE PARITY HOLD ON REAL cfRNA AT FULL SCALE?** |
| `design/accumulator_cost.py` | **HOW MANY ns PER FRAGMENT DOES THE ACCUMULATOR COST, regressed over several BAMs?** |


## CLI

```bash
rigel index --fasta genome.fa --gtf annotation.gtf -o index/
rigel quant --bam sample.bam --index index/ -o results/
rigel sim --config scenario.yaml -o out/
rigel export results/ -f tsv
rigel report results/ -o report.html
```

Input BAM must be name-sorted with the `NH` tag.
