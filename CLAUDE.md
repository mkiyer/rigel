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

0.8.0 is A RELEASE OF THE TOOL (owner, 2026-09-19; `docs/DESIGN.md` §0b's amendment). TWO numbers are primary
and they answer different questions: the transcript table against per-transcript truth is what the release
ships on (`quant_accuracy.py`, read per stratum and only above `--arm base_reseed`, its noise floor), and the
calibration result against an oracle calibration is what ranks a CALIBRATION mechanism
(`calibration_vs_oracle.py`, `solvability_audit.py`, `prior_vs_oracle.py`). Neither stands in for the other,
and ranking a calibration mechanism on the transcript number stays refused: the table is downstream of both
calibration and the EM, so one figure cannot say which moved. Three strata are in scope — unstranded ×
capture-OFF, stranded × capture-OFF, stranded × capture-ON — and **unstranded × capture-ON is DEFERRED**:
still reported on every benchmark, never a development target until the other three are optimised. It is
also where most of the error is (the gDNA fraction cancels from the strand mean, so an unstranded AMBIG
slot has no channel), so the debug loop must take the worst IN-SCOPE scenario, never the deferred one.

The fragment-length COMPOSITION channel is retired until after 0.8.0 and may not be proposed; it does not
exist in `src/`. Three other things called "length" are unaffected: layer 2's `fl` / `effective_length` /
`capture_eff_length` (the opportunity model), `length_likelihood` in `second_pass.py` (per-fragment
assignment), and the fl PMFs themselves (`calibration.fl.FLModels`). The ladder gives gDNA and RNA EQUAL fragment
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
| which strand a fragment came from | **4 · strand** — `gdna_strand` `strand_balance` `strand_summary` |
| how dense a component is, and the priors | **5 · density and prior** — `density_model` `density_deconv` `landscape` `abundance_landscape` |
| what one neighbour tells another | **6 · the solve** — `sweep` (the backbone) + `blocks` (the chain view's fields and the diagnostic capture) + `messages/` (the policy) + `region_init` |
| turning the solve into a result | **7 · assemble** — `calibrate` `priors` `result` `derive` `diagnostics` `track` |

## The message layer

`CalibrationConfig.message_policy = "transfer"` ships (the default since 2026-09-09; `"silent"` is the
floor, and the one field selects — the `message_propagation` switch retired 2026-09-13). Two policies, selected by one config value (an unknown name raises), both on the two-phase
backbone (`docs/DESIGN.md` §6b.11–§6b.12): every node's own claim → a forward pass and a backward pass in
which each recipient receives what its neighbour sends, so every node ends with one message from each
neighbour it has → the solve, which hands ψ two row channels and nothing else. THE SWEEP IS ONE NATIVE CALL
(`native.solve_blocks`, `native/solve_kernel.cpp`, since 2026-09-18): the chain's locus blocks are solved on a pool
of threads, one block at a time, end to end — the prior rows, the self-solve, the layer, the final solve, the
write-back — bit-identical at every thread count; `sweep.solve_chain` cuts the blocks, reduces the library, makes
the call and judges the assertions' counts.

| policy | |
|---|---|
| `silent` | the measured floor (`messages/silent.py`) |
| `transfer` | the shipped default (`messages/transfer.py`: the policy's name, its strand model and its LIBRARY — the three level lanes' coordinates and the strand witness's liveness, the only cross-block reduction a message may use). Everything else is the kernel's (`native/transfer_kernel.h`, run per block inside `native.solve_blocks`): the builders — a table of contents of named builders, one per message: `claims`, `splice_faces`, `edge_level`, `terminus_rules`, `alternative_splice_site`, `gdna_lane`, `rna_lane` — the two passes and the solve, on the row constructors of `native/transfer_rows.h` (bound for the gates as `native.transfer_rows`). Every hop pays its pair's counting plus the disagreement beyond it. The gates read the kernel's tables through `native.transfer_prepare` / `transfer_pass` / `transfer_solve` (`tests/calibration/_transfer_harness.py`) |

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
  before it has been derived on paper, prototyped OUTSIDE THE MAIN TREE, and A/B'd against what already
  ships, on the same conditions. The solve is native, so anything inside the block — a builder, a rule, a
  row constructor, ψ — is prototyped in C++ in a worktree and the two trees are scored against each other
  (`policy_benchmark.py --by-class`, `calibration_vs_oracle.py`); what is still Python (the library, the
  fits, the assembly) is prototyped in Python the same way. One mechanism at a time: a change that cannot
  be A/B'd alone cannot be judged alone.
- **No magic numbers.** Stop and discuss before adding any constant, heuristic or tunable. Every divisor
  must be derived from the deposit rule and unit-tested against brute-force enumeration.
- **A falsification test first, verified failing — then break the fixed code and watch each gate fire.**
  The second half is not optional; it has found holes in already-green gates repeatedly.
- **The debug loop is the default method**: run the panel → take the worst IN-SCOPE scenario → dissect it
  to the highest-error objects (`worst_objects.py`, `calibration_walk.py`) → find the mechanism → fix →
  re-run the panel.
- **A ceiling is sometimes the right instrument** (`quant_accuracy.py`'s injection arms)
  but prices something that may be unreachable, so it is not the default. Every prior-injection arm patches
  `assemble_priors`, while `pipeline.py` builds `effective_lengths_em` before calling it, so the
  effective-length shrinkage is inside no prior-injection ceiling; only `oracle_ruler`, which hands the EM
  the simulator's own capture-aware length, reaches it. Benchmark every arm with
  `--set em.assignment_mode=fractional` (owner, 2026-09-19): the shipped assignment is a sampled draw.
- **One thing varied per experiment**, a baseline re-recorded from the current tree in the same session,
  and **score against truth** (the oracle BAM's read names), never against the previous run.
- **No legacy, no backwards compatibility, no speculative code.** Converge and delete. No version
  suffixes in file names. No Greek letters in identifiers (fine in maths write-ups).
- **Real data is a test input, never a design input.** Sweep the plausible space, report the worst case.
  Profile on a deep real RNA-seq library, never a panel condition or a toy; the cfRNA libraries are smoke
  tests (`docs/TESTING.md` §6). Read a timing only from back-to-back A/B pairs — untouched stages drift
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
python -m pytest tests/ -q
python -m pytest tests/ --update-golden        # regenerate tests/golden/ after intended output changes
ruff check src/ tests/ scripts/ && ruff format src/ tests/   # never format scripts/
```

**The standing baseline: 0 failed / 3,440 passed / 0 skipped / 2 xfail, 3,442 collected** (re-derived
2026-09-19 after `quant_accuracy.py --markdown`, the per-scenario release report: +8 in
`tests/calibration/test_quant_accuracy.py`, every one a way the RENDERING can lie while the numbers
underneath are right — the gDNA column scored without its intergenic half, a ratio invented at a truth of
zero, the two RNA pools folded together, a ragged table, the deferred stratum unmarked, a column that
CANNOT FIRE reported as a score (`fn_mass` is identically 0 under fractional assignment), and the truth
table's ROW COUNT passed off as a transcript count (15,669 rows, 8,750 annotated; 9,385 gene rows, 2,466
real genes). Before that 3,432 / 3,434 after nascent RNA's share of the RNA prior was restored: +5 in `test_estimator.py` — five gates on
the restored allocation (two parametrised over four locus shapes, three over both EM modes) replacing the four
that pinned the eligibility test and the one that pinned its flag — and −1 in
`tests/native/test_grouped_prior_update.py`, where the identity's table loses its all-synthetic row and the
mask's byte-identity pair becomes one property gate over the space; THREE XFAILS CLOSED by the same landing
(5 → 2), which moves `passed` and not `collected` — `ISSUES: nascent-gets-no-rna-prior` and both rungs of
`ISSUES: nested-antisense-leak-under-the-sane-ruler`. ⚠ The goldens moved: 5 of 21 scenarios materially (the
ones holding a LIVE nascent entity — nascent mass up, `gdna_rate` down, e.g. `combo_moderate` 52.93 → 171.21
and 0.5672 → 0.3958), the other 16 at the denormal floor only. The calibration columns
`gdna_prior_count` / `rna_prior_count` also moved by 1e-16 to 1e-12, and that is NOT this change: regenerating
at HEAD moves them by byte-identical amounts, so the committed goldens carried that staleness already, under
`rtol=1e-6`. Before that 3,425 / 3,430 after the simulator's capture physics: +2 in `test_sim_capture.py` — gDNA and cDNA bind the same half
of a split probe alike, and a capture key the loader does not know is refused — and +1 for
`docs/dev/CAPTURE_WITHOUT_THE_PANEL.md` by the `docs/dev/` row; before that 3,422 / 3,427 after the cache key and
the ruler arm: +2 in `test_scan_cache.py`, where three gates — a thread count is
not the key, a tally setting is, the key is derived from the recorded settings — replaced the one that tampered a
stored digest; +4 in `test_quant_accuracy.py` — the ruler arm reaching the EM with its noop inert, its guard, the
capture truth's anchor, the refusal of mixed assignment modes; before that 3,416 / 3,421 when the seven SPENT
sandbox notes were retired — every plan whose work had landed and whose record
had moved to a permanent doc: −7 by the `docs/dev/` row; before that −1 for the performance plan itself, when
it was retired into `DESIGN.md` §6b.15 and `ISSUES: performance-memory-bounded-solve`; before that after its
phase 4: +2 gates in `test_fl.py` — the adjacent-pair table against a
hand-built three-reference walk, and the two boundary classes told apart on an asymmetric fixture; before that
after phase 3: +2 gates — the batched sj lookup's parity with the scalar rule in
`tests/native/test_accumulator_native_parity.py` and the two-reference id-and-motif gate in
`test_second_pass_scoring.py`; before that after phase 2: +1 the split's arithmetic gate in
`test_scan_order_independence.py`; before that after phase 1: +4 gates — the bisection's bit-equality against the old fixed
loop and the log-gamma table's against the direct form in `test_gdna_density.py`, the one-traversal spy and the
handed-in shares in `test_priors.py`; before that +1 for `docs/dev/PERFORMANCE_PLAN.md` by the `docs/dev/` row;
before that after `strand_likelihood.py` was converged into the gates' oracle module: −3 by the module row; before that after `policy_prototype.py` was retired: −4 by the `scripts/design/` row; before that after the message cache was deleted: −3 for `message_cache.py` by the module row, −6 gates — the five cache
gates of `test_sweep_backbone.py` and the served-injection gate of `test_landscape_training_population.py`; before that
after the block went into one native call: −10 for the four files deleted — `messages/faces.py` and
`messages/lanes.py` by the module row (−3 each), `native/psi_kernel.cpp` and `native/transfer_kernel.cpp` by the row
below (−2 each) — +6 for `native/solve_kernel.cpp`, `transfer_kernel.h` and `psi_kernel.h` (+2 each), −1 the `RowTable`
gate and −1 the no-copy gate (the Python tables they tested are gone), −2 the two backbone gates on those tables (a
lane reaching the solve in the same table, the block view's slicing), +1 the thread gate on `solve_chain`; before that
after ψ took its priors apart: +1 the bit-equality gate on the kernel's arm in `test_landscape.py`, +1 for
`docs/dev/BLOCK_NATIVE_PLAN.md` by the `docs/dev/` row; before that after the cache key moved to the factory's inputs: +2 the digest gates in `test_sweep_backbone.py`; before
that after ψ went threaded: +2 the thread-exactness and budget gates in `test_sweep.py`; before that after
the one-path convergence of the pass: −6 for `test_pass_kernel.py` (its four gates and the `tests/` row's 2),
+3 for the gates that moved or joined (the wiring and no-copy gates in `test_transfer_policy.py`, the trigamma gate in
`test_zero_count_is_a_measurement.py`), −3 for `messages/transfer_rows.py` by the module row; before that after the
tables' allocations: +1 the poison gate in `test_transfer_policy.py`, +1 `docs/dev/THREADS_PLAN.md`
by the `docs/dev/` row; before that after the solve went native and the transfer's two `.cpp` files became one `transfer_kernel.cpp`: −2 by
the row below; before that after ψ went native: +4 — +2 for `native/psi_kernel.cpp` and +2 for the gates' oracle module
`tests/calibration/_psi_reference.py` by the rows below, the rewritten ψ gates moving nothing; before that
after the one-path cleanup deleted the Python builders and their two-kernel gate: −7 — the five gates
of `test_prepare_kernel.py` and its `tests/` row's 2 — and +1 for `docs/dev/PSI_PORT_PLAN.md` by the `docs/dev/`
row; before that +11 after the native builders — those five gates
with the `tests/` row's +2, +2 for `native/transfer_kernel.cpp` and +2 for `native/transfer_rows.h` by the row
below; before that +2 after the
builders' layout — the `RowTable` gate in `test_transfer_faces.py` and the
tables-without-a-copy gate in `test_pass_kernel.py`; before that +7 after the native pass — the three gates of
`test_pass_kernel.py` with the `tests/` row's +2, and +2
for `native/pass_kernel.cpp` by the row below; before that +5 when the yield's floors were deleted — the
common-thinning invariance and the cannot-emit rule in `test_estimator.py`, each over both EM modes, and the
unfloored locus yield in `test_priors.py`; and +33 from 3,386 with the ruler's repair of 2026-09-16, file by file:
the members rule +3 — two selector gates in `test_abundance_landscape.py`, the result-schema gate for
`gdna_reference_members`; the truth instrument +5 — `ruler_vs_truth.py` by the `scripts/design/` row, its self-test
in `test_instrument_self_tests.py`; the expectation ruler +25 — `capture_efficiency.py` by the module row (+3) with
its seven gates in `test_capture_efficiency.py` (+2 by the `tests/` row, +7), the taper, crossing-share and
per-interval-length gates in `test_effective_length.py` (+18), `test_capture_eff_length.py` rewritten 15 → 12,
`test_priors.py` 29 → 27 with the floor's tests replaced by the object-form length tests; the goldens unchanged to
the bit). The 2 xfails are executable
records of proven defects whose fixes are elsewhere
(`ISSUES: two-sided-exon-row`; `ISSUES: the-lower-bound-noise-ratchet`, the encompassing locus's shallow
flank under an edge level), deferred by ruling to their threads — "fix the test" is a category error, and an xfail is closed
by repairing the thing or asserting the invariant structurally, never by widening a bound. **Any failure at all is a regression.** A commit that measures the suite updates this line.

**Re-derive a count, never adjust one** (`TRAPS: re-record-the-baseline`). Several gates are parametrised
over the files on disk, so adding or retiring a file moves the total; account for it from this table and
confirm with `pytest --collect-only -q | grep <stem>`:

| adding one… | moves collected by | which cases |
|---|---|---|
| `src/rigel/calibration/` module | **+3** | jargon, docs-boundary, and layering *if declared in `_layers.py`* |
| `tests/` file (any directory) | **+2** | jargon, docs-boundary |
| `src/rigel/` file outside `calibration/` (a `.py`, or a `native/*.cpp` or `*.h`) | **+2** | jargon, docs-boundary |
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
| `design/preflight.py` | ⭐⭐⭐ **CAN THIS SESSION RUN AND REGENERATE EVERYTHING? — one command, one verdict, before anything else.** Checks the toolchain (the `rigel` env, the native extension, the CLI), both references, both panels (scan caches, oracle caches with every part `panel.py` requires, the certified `slot_truth`) and that every `scripts/design/` instrument IMPORTS. ⛔ It changes nothing and measures nothing — every check is a read or an import, and a ✘ prints the exact command that regenerates the missing artifact. ⭐ **The default is ~2 s**; `--full` adds every instrument's `--self-test`, in seconds. `--self-test` 7/7 |
| **⭐⭐⭐ THE POLICY BENCHMARK — where a message-policy change is judged** | |
| `design/policy_benchmark.py` | ⭐⭐⭐ **HOW DOES EACH POLICY SCORE, PER CONDITION, AGAINST CERTIFIED TRUTH?** Whole-library gDNA error in fragments, per axis, one row per condition, for `silent` and `transfer`. ⭐ `--panel test` is the test chromosome (seconds — the development loop); `--panel ladder` is the 16-condition benchmark. ⭐ **`--by-class`: WHERE DOES A POLICY'S REMAINING ERROR SIT, BY NODE CLASS?** — per certified stratum, boundaries split by terminus flag, exons by reach (licensed face / edge only / walled); the instrument that ranks the rebuild's holes. `--set SECTION.FIELD=VALUE` applies a config value on top of every policy, the same spelling as `calibration_vs_oracle.py`. ⛔ NEVER POOLED, and the two halves are judged against DIFFERENT bars: unstranded rows are where a policy must WIN, stranded rows are where it must do minimal HARM against silence |
| **⭐⭐⭐ 0.8.0'S METRIC — calibration against ORACLE CALIBRATION** | |
| `design/calibration_vs_oracle.py` | ⭐⭐⭐ **IS THE CALIBRATION RESULT ITSELF RIGHT, SCORED AGAINST AN ORACLE CALIBRATION? — 0.8.0's metric, and it reaches the effective-length shrinkage, which no prior-injection arm does.** `P = calibrate(...)` against the same payload with only the six deconvolved arrays swapped, per stratum (the ruler's reference is the result's own, so at capture-OFF `O` is the no-enrichment null with no fitting); `--set SECTION.FIELD=VALUE` prices any config value on both arms — `--set calibration.message_policy=silent` is the ship protocol's first item — so a policy or a grid arm is a config value and nothing in `src/` moves to price it. ⛔ Read `ruler_n_moved`, never the aggregate: the total can barely move while nearly every transcript is redistributed. No solver, no EM, no re-scan — ~5–12 s/condition. `--self-test` 39/39 |
| `design/ruler_vs_truth.py` | ⭐⭐⭐ **IS THE RULER'S FORMULA RIGHT? — the EM's effective length under capture, per transcript, against the simulator's own capture-aware effective length** (`CaptureSampler.partition_array`, the truth the reads were drawn with), anchored on the fully probed transcripts and read per class (the probed fraction of the transcript's bases, off the sampler itself) and kind. Arms: `shipped`; `oracle_gdna` (the certified true gDNA counts through the shipped ruler — the ideal witness, so what remains under it is the ruler's or the panel geometry's); `--module` prototype rulers beside them, the harness every ruler mechanism is developed on before `src/`; `--set` prices a config value on every arm. `--panel-dir` takes the depth ladder (`scenarios_depth_d10` / `_d100` / `_full_lowg`, `docs/TESTING.md` §0c) for the operating curve of the reference. ⛔ Read the classes apart: the probed class is where the formula is exact when its witness is, the unprobed class is where a floor or a reference bites. `--self-test` 20/20 |
| `design/calibration_oracle.py` | ⭐⭐⭐ **WHAT IS THE CERTIFIED PER-OBJECT TRUTH? — run this before debugging calibration against anything.** Every REGION and BOUNDARY's count, its realized `n_gdna`/`n_nrna`/`n_mrna` and `true_f_g`, at two certification levels: COMPOSITION (no opportunity model anywhere in it) and FIELD (densities too). ⛔ REFUSED unless its named gates pass — sum-to-full, partition-projects-exactly, gdna-field-uniformity, exact-zeros, nascent-in-annotation, rna-strands-close — because a merely plausible oracle is how a calibration bug and a truth bug survive each other. Writes `slot_truth.npz` beside each oracle cache; `--self-test` 13/13 |
| `design/calibration_walk.py` | ⭐⭐⭐ **WHICH STAGE OF CALIBRATION INTRODUCES THE ERROR?** The solve as a ladder — init → strand → local → +messages → +refits → shipped — each rung scored per stratum against `calibration_oracle.py`, which it refuses to run without |
| `design/solvability_audit.py` | ⭐⭐⭐ **WHICH OBJECTS ARE SOLVABLE, WHICH ARE SOLVED WRONG, AND WHICH ARE CONFIDENTLY WRONG? — where pass-0 and 0.8.0 are judged.** ⛔ Honest ignorance is excluded: `f_g ≈ ½` at zero precision with no own evidence is correct. Omit `--condition` to run the panel |
| `design/prior_vs_oracle.py` | ⭐⭐⭐ **IS `LocusPriors` — the thing the EM actually reads — RIGHT?** Five arms separate calibration's own error from the assembler's, reporting the count, the composition claim and the scale apart, per stratum, in the drained frame |
| `design/pass0_vs_oracle.py` | **HOW DOES PASS-0 COMPARE WITH THE ORIGIN-SPLIT PAYLOAD AND TWO LEVERED CEILINGS, per object and per class?** ⛔ Its mass-weighted headline is the wrong yardstick for pass-0 — honest ignorance reads as error there |
| `design/worst_objects.py` | ⭐⭐ **WHICH REGIONS AND BOUNDARIES CARRY ONE CONDITION'S ERROR MASS?** Read the concentration curve first — concentrated means a mechanism exists, diffuse a systematic bias; `fg_loc` vs `pred_fg` separates a bad local solve from bad messages |
| **⭐⭐⭐ the panel, and the caches that make calibration a seconds-long loop** | |
| `sim/configs/flgap_rna_long.yaml` · `flgap_rna_short.yaml` | ⭐⭐ **WHAT BREAKS WHEN gDNA AND RNA FRAGMENT LENGTHS DIFFER? — the fl-gap SIDE panel, two arms of opposite sign** at `g50`. ⛔ Never a ladder rung: its transcript-level number is not a calibration result, though everything before the EM is valid ⛔⛔ **NOT REGENERATED on 2026-08-22: both arms still carry the RETIRED UNIFORM nascent model (`mode: fragment_share`), while the ladder carries SPARSE.** Each panel's on-disk data matches its own config, so both are internally accurate and their recorded measurements stand — but a ladder-vs-side-panel comparison now varies TWO things, so no claim may be carried across them. Re-simulating is an open owner decision |
| `sim/configs/gdna_ladder.yaml` | ⭐ **THE STAGE-B PANEL, AND THE ONLY PANEL THE TOOL IS RANKED ON** — 16 conditions (`g00/g05/g50/g98` × ss `0.50/0.99` × capture off/on), gDNA 0 → 98 % at a fixed 10 M total. ⛔⛔ **EQUAL FRAGMENT LENGTHS, and that is a forcing function**: the EM already reads the fl distribution, so a length gap lets it split the origins on LENGTH ALONE and mask calibration bugs (owner, 2026-08-14). ⛔⛔ **Every row carries SPARSE nascent RNA (2026-08-22): `on_fraction 0.50` of gene SPANS on, level logU(1, 100) INDEPENDENT of the mature level — 20.2 % of RNA fragments, the retired uniform model's total distributed sparsely.** 0.50 is a DEVELOPMENT STRESS level and not real data (`DESIGN.md` §0b, THE NASCENT SCOPE RULING); realistic is 0.10 ⇒ 4.2 %. `docs/TESTING.md` §0 |
| `sim/panel.py` | ⭐⭐⭐ **HOW IS A PANEL BUILT, SIMULATED, CACHED AND SCORED? — one command per stage**: `status` / `build` / `simulate` / `cache` / `score` / `report`, every path derived from one panel YAML. ⭐ Run `status` FIRST — every stage is expensive and resumable, and it names the next one. ⛔ It adds no measurement code: each stage shells out to the instrument that already owns it. ⛔ `cache` builds BOTH caches, and the oracle one is the origin-split truth every scoring instrument reads. Gated by `tests/test_panel_workflow.py` |
| `design/build_scan_cache.py` | **SCAN ONCE, CALIBRATE MANY TIMES.** ⛔ The cache key hashes `accumulator.cpp`'s deposit rule and not `resolve.cpp`'s fragment construction, so for a change to which fragments are OFFERED use `--force` or delete the caches by hand |
| `sim/build_suite_reference.py` · `design_suite_probes.py` · `simulate_reads.py` | **HOW IS THE PANEL'S SUBSTRATE BUILT?** ⚠ `panel.py build` drives the last two; the reference carve needs the source genome/GTF, which a panel config does not name, so it stays manual |
| **⭐⭐⭐ the prior assembler, and THE NUMBER THE RELEASE SHIPS ON** | |
| `design/quant_accuracy.py` | ⭐⭐⭐ **HOW ACCURATE IS THE TOOL END TO END, AND WHAT IS A PERFECT PRIOR WORTH?** `--arm base` plus the oracle and per-field injection arms, scored count against count. ⭐ Read per stratum, above `--arm base_reseed`, and under `--set em.assignment_mode=fractional` (owner, 2026-09-19): one of 0.8.0's TWO primary numbers beside the calibration metric, never a stand-in for it. ⭐ **`--report FILES… --markdown OUT`: THE PER-SCENARIO RELEASE REPORT** — every condition's three pools (gDNA / SYNTHETIC nascent / annotated) against truth in raw counts and per cent, then transcript and gene error inside the annotated pool alone, then the per-stratum roll-up. It renders arm jsonl and runs nothing |
| **⭐⭐⭐ where to develop** | |
| `design/rename_identity.py` | ⭐⭐⭐ **IS THIS RENAME, REFACTOR OR SPEED-UP NUMERICALLY A NO-OP?** `--freeze` captures one reference, `--check` compares after every stage — on array CONTENT and the transcript table, never on names; `--bam` takes a real library instead of a panel condition. ⚠ The reference is frozen, never rolling. `--self-test` 8/8 |
| `design/rename_census.py` | ⭐⭐⭐ **WHICH NAMES DOES A VOCABULARY RULING TOUCH, AND WHICH CARRY TWO SENSES?** Reports by kind — identifiers, C++, prose — and never renames; `--sense <token>` dumps every site with context. ⛔ Run it before renaming anything |
| `design/module_census.py` | ⭐⭐⭐ **WHERE DOES A CHANGE GO?** The calibration package re-derived from the AST: the layering with every upward import, each module's importers, docstrings naming a sibling with no import, dead public surface. ⛔ It reports; it does not judge |
| **⭐⭐ the toy harness** | |
| `design/zero_controls.py` | ⭐⭐⭐ **DOES THE TOOL HOLD AT ZERO RNA AND AT ZERO gDNA? — the owner requires both on every experiment.** The truth is a constant, so every deviation is a false positive. ⛔ Flags any EMPTY object: a degenerate zero arm tests nothing |
| `design/toy_harness.py` | ⭐⭐ **HOW DOES A MINI CHROMOSOME YOU DEFINE CALIBRATE — in 0.1–5 s, with every object's answer beside its truth?** (`docs/TESTING.md` §0b) The priors a toy cannot fit are harvested from a real cached condition; `--list` for the ladder |
| **the substrate — are the panel and the index sound?** | |
| `design/simulator_gates.py` | **DOES THE SIMULATOR PASS ITS OWN GATES, scored on per-fragment truth?** ⛔ Run it before trusting the panel |


## CLI

```bash
rigel index --fasta genome.fa --gtf annotation.gtf -o index/
rigel quant --bam sample.bam --index index/ -o results/
rigel sim --config scenario.yaml -o out/
rigel export results/ -f tsv
rigel report results/ -o report.html
```

Input BAM must be name-sorted with the `NH` tag.
