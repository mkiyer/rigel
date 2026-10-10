# CLAUDE.md

A session primer: what Rigel is, what 0.8.0 is about, where a change goes, the working rules, and which
instrument answers which question. It is not a history (git is) and not a manual (the nine docs below are).
A rule earns its place here by changing what the next session does.

## What Rigel is

A Bayesian RNA-seq transcript quantifier that separates RNA from genomic-DNA contamination. A single-pass
C++ BAM scanner tallies fragments, a **calibration** stage deconvolves the library into gDNA vs RNA, and a
per-locus EM assigns RNA to transcripts. PyPI package `rigel-rnaseq`; the import and CLI are `rigel`. The
version on disk is 0.7.1 (`pyproject.toml`); the target is **0.8.0**.

## Axiom 0 — RNA is RNA

There are three populations and no fourth: `gDNA`, `RNA+`, `RNA−`. "Mature" and "nascent" are not
populations and not a degree of freedom: RNA inside an intron is RNA that has not spliced at that
position. The only distinction is whether a fragment **is spliced** (certified RNA, needing no
deconvolution) or **is not** (the whole deconvolution problem). The population set at a slot is a function
of two bits:

```
T(slot) = {gDNA}                                  # always — gDNA is genomically continuous
        ∪ {RNA+ if statics.free_pos[slot]}        # iff the annotation admits that strand here
        ∪ {RNA− if statics.free_neg[slot]}
```

so `|T| ∈ {1,2,3}`. The tell that you have violated this: a population set with more than three members,
or the words "mature", "nascent" or "a third component" in a solver or composition question. Re-ask it as
"what is this channel's OPPORTUNITY for RNA at this object?" — a geometry, derivable from the index. The
words survive only as simulator inputs (`nrna_abundance`). Nascent RNA is sparse in real data and is
modelled for robustness, so a number measured at the panel's nascent share is a stress reading, never a
design driver (`docs/DESIGN.md` §0b, the nascent scope ruling).

## The 0.8.0 scope

0.8.0 is a release of the TOOL; the ruling is `docs/DESIGN.md` §0b. Two numbers are primary and neither
stands in for the other: the transcript table against per-transcript truth is what the release ships on
(`quant_accuracy.py`), and the calibration result against an oracle calibration ranks a CALIBRATION
mechanism (`calibration_vs_oracle.py`, `prior_vs_oracle.py`). Ranking a calibration mechanism on the
transcript number stays refused: the table is downstream of both calibration and the EM, so one figure
cannot say which moved.

**The pinned A/B.** An end-to-end A/B pair runs with the scan pinned (`--set scan.total_threads=1`;
`OMP_NUM_THREADS=1` does not pin it) and under `--set em.assignment_mode=fractional` (the shipped assignment
is a sampled draw), so it is exactly reproducible; an effect is judged by its size, with genes and pools read
beside the transcript table.

Three strata are in scope — unstranded × capture-OFF, stranded × capture-OFF, stranded × capture-ON — and
**unstranded × capture-ON is DEFERRED**: reported on every benchmark, not a 0.8.0 release requirement.
The owner explicitly authorizes investigating the factory-removal regression in this stratum for
robustness (2026-10-07); it does not replace the three in-scope release targets. The normal debug loop
takes the worst IN-SCOPE scenario. Every score is read per stratum, never pooled.

Fragment length is back in scope for 0.8.0 (owner, 2026-10-02): the composition channel, the capture-length
frame and every other length item that waited for the release (`docs/DESIGN.md` §0b;
`ISSUES: calibration-detects-capture-on-a-capture-off-library`, `ISSUES: the-scorer-reads-a-census-length-law`). The ladder still gives gDNA and RNA EQUAL fragment lengths on purpose:
a gap would let the EM split origins on length alone and mask calibration bugs. So a length mechanism is
judged on the gap panels beside it, in both directions, against the equal-length control.

## The docs

Nine permanent docs, none of them a changelog; each one's preamble says what it holds and what it does not.

| | |
|---|---|
| `docs/SUCCESS.md` | **Start here.** How performance is judged, and the instruments in run order |
| `docs/ROADMAP.md` | The ranked view: one line per state claim naming its instrument, and the ordered next steps |
| `docs/ISSUES.md` | Every open problem as a named entry, and the CLOSED / REFUSED record with each killing number |
| `docs/TRAPS.md` | Mistakes already made, as named rules |
| `docs/EQUATIONS.md` | The derivations the code depends on |
| `docs/DESIGN.md` | What is built and the rulings behind it: settled, not re-litigated. §0 is the binding vocabulary |
| `docs/TESTING.md` | The panels, the test chromosome (§0a), the gates, and what the suite can judge |
| `docs/MANUAL.md`, `docs/PUBLISHING.md` | The user's manual; the release procedure |

- **The move rule.** When a finding settles, MOVE it to its one home (an open problem or refusal to ISSUES, a
  lesson to TRAPS, a ruling to DESIGN, a derivation to EQUATIONS) and delete it where it was, in the same
  edit. Never copy: two homes diverge.
- **Cite a rule by its name** — `TRAPS: off-grid-message-mode`, `ISSUES: two-sided-exon-row` — never by a number.
- **`docs/dev/` is the sandbox**: nothing there is authoritative and nothing may cite into it. **The source
  cites no doc**: a docstring may cite a test or `tests/native/_accumulator_reference.py`, never a doc.
  `tests/test_docs_boundary.py` gates the naming, sandbox and source rules.

## Where does a change go? — the calibration layering

An import may point DOWN a layer or SIDEWAYS within one, never UP; a module reaching up is telling you the
thing belongs lower. `rigel/calibration/_layers.py` is authoritative, `tests/calibration/test_layering.py`
enforces it.

| if the change is about… | it goes in |
|---|---|
| what a fragment tally MEANS | **1 · the payload view** — `splice_graph` `substrate` `region_arrays` |
| how many places a fragment COULD have sat | **2 · opportunity** — `effective_length` `capture_eff_length` `sj_opportunity` `gdna_opportunity` `gdna_density` `fl` |
| one slot's own numbers and ψ | **3 · geometry + the per-slot solve** — `region_geometry` `simplex_logodds` |
| which strand a fragment came from | **4 · strand** — `strand_balance` |
| how dense a component is, and the priors | **5 · density and prior** — `density_model` `landscape` `capture_efficiency` |
| what one neighbour tells another | **6 · the solve** — `sweep` (the backbone) + `blocks` + `messages/` (the policy) + `region_init` |
| turning the solve into a result | **7 · assemble** — `calibrate` `priors` `result` `derive` `track` |

**The message layer**: `CalibrationConfig.message_policy = "transfer"` ships and `"silent"` is the measured
floor; the sweep is one native call (`native.solve_blocks`). The backbone, the policies and their rulings are
`docs/DESIGN.md` §6.1 and §6b; read `ISSUES: the-message-policy-campaign` before re-proposing a mechanism.

**The benchmarks**: develop on the test chromosome (`policy_benchmark.py --panel test`, seconds), decide on
the ladder (`--panel ladder`). Unstranded rows must WIN against silent, stranded rows must do minimal HARM;
the halves are never pooled. `docs/TESTING.md` §0a has the recipe and the bars, `docs/SUCCESS.md` the run order.

## Working rules

- **DERIVE → DESIGN → PLAN → PROTOTYPE → A/B → only then `src/`.** Nothing enters the production source
  before it has been derived on paper, prototyped OUTSIDE THE MAIN TREE and A/B'd against what ships, on
  the same conditions. Anything inside the native block — a builder, a rule, a row constructor, ψ — is
  prototyped in C++ in a worktree and the two trees are scored against each other
  (`policy_benchmark.py --by-class`, `calibration_vs_oracle.py`); what is still Python is prototyped in
  Python the same way. One mechanism at a time: a change that cannot be A/B'd alone cannot be judged alone.
- **No magic numbers.** Stop and discuss before adding any constant, heuristic or tunable. Every divisor is
  derived from the deposit rule and unit-tested against brute-force enumeration. A fixed number lives in
  `config.CONSTANTS`, documented with why it has its value (a run setting in `PipelineConfig` or `IndexConfig`);
  `tests/test_constants.py` refuses a bare module constant.
- **A falsification test first, verified failing — then break the fixed code and watch each gate fire.**
  The second half is not optional; it has found holes in already-green gates repeatedly.
- **The debug loop is the default method**: run the panel → take the worst IN-SCOPE scenario → dissect it
  to the highest-error node class (`policy_benchmark.py --by-class`) → find the mechanism → fix → re-run
  the panel. A ceiling (`quant_accuracy.py`'s injection arms) prices what may be unreachable, so it is not
  the default (`TRAPS: measure-the-ceiling-first`; which arm reaches which length is `docs/SUCCESS.md`'s).
- **One thing varied per experiment**, a baseline re-recorded from the current tree in the same session,
  and **scored against truth** (the oracle BAM's read names), never against the previous run.
- **No legacy, no backwards compatibility, no speculative code.** Converge and delete. No version suffixes
  in file names; no Greek letters in identifiers.
- **Real data is a test input, never a design input** (`TRAPS: real-data-is-a-test-input`). Profile on a
  deep real library, read a timing only from back-to-back A/B pairs, and prove every speed-up a numeric
  no-op (`docs/TESTING.md` §6).
- **Renames**: grep every sense of a token first, and run `rename_identity.py --check` after each stage.
  `arm` has three senses (an experiment arm, a component arm, `__ARM_NEON` in the scanner) and needs a
  pass per sense, never a tail-end sweep.
- **The owner drives commits.** Do not commit unless asked.

## Build, test, lint

Everything runs inside the `rigel` conda env, which holds htslib and the compilers (the C++ build finds
htslib via `$CONDA_PREFIX`).

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
export OMP_NUM_THREADS=1                        # always, when benchmarking or comparing runs

pip install --no-build-isolation -e ".[dev]"    # rebuild after ANY src/rigel/native/ change
python -m pytest tests/ -q
python -m pytest tests/ -q 2>&1 | grep '^FAILED' | sed 's/::.*//' | sort | uniq -c   # the failure set
python -m pytest tests/ --update-golden         # regenerate tests/golden/ after an intended output change
ruff check src/ tests/ scripts/ && ruff format src/ tests/   # never format scripts/
python scripts/design/preflight.py --full       # every instrument's --self-test
```

**The standing baseline: 1 failed / 2,804 passed / 0 skipped / 0 xfail, 2,805 collected** — re-derived
2026-10-10 after the capture reader landed (the goldens regenerated and read against truth). The one failure is
`test_multimap_counting.py::TestParalogMultimapping::test_gdna_sweep[gdna_100]`, red by ruling until the
accumulator's second pass assigns multimappers (`ISSUES: the-calibration-count-is-blind-to-multimappers`);
never rewrite it to green. **Any other failure is a regression.** Re-derive a count, never adjust one (`TRAPS: re-record-the-baseline`): every gate that scans
the files on disk is one case, so only adding or removing a test moves the total. Derive the failure set,
never eyeball the tail (`TRAPS: read-the-whole-failure-list`). A golden update is where a regression gets
laundered into "intended": read the diff, record its magnitude and check the truth-scored instruments
before `--update-golden`. After a `src/` deletion, a rename or a shipped-default flip, run the instruments
and not only the suite — `preflight.py --full` (`TRAPS: a-green-suite-hid-five-dead-instruments`).

## Tooling under `scripts/`

One row per instrument; its module docstring and `--help` are its manual, and `docs/SUCCESS.md` has the run
order. `tests/test_scripts_index.py` holds the `design/` rows against the disk in both directions;
`scripts/profiling/` is indexed in `scripts/README.md`.

| | the question it answers |
|---|---|
| `design/preflight.py` | START HERE: can this session run and regenerate everything? Reads and imports only; `--full` adds every `--self-test` |
| `design/policy_benchmark.py` | How does each message policy score per condition against certified truth? `--by-class`: where the error sits. ⛔ Never pooled: unstranded must win, stranded must do minimal harm |
| `design/calibration_vs_oracle.py` | Is the calibration result right against an oracle calibration? 0.8.0's calibration metric; `--set` prices a config value on both arms. ⛔ It cannot see the capture-contracted length (`O` keeps `P`'s efficiencies): read `ruler_vs_truth.py` beside every capture-ON arm |
| `design/ruler_vs_truth.py` | Is the EM's capture-contracted length right against the simulator's own yield? `--scale`: every class on one scale, and the within-gene spread; `--module` runs prototype rulers. ⛔ Read the probed and unprobed classes apart |
| `design/calibration_oracle.py` | What is the certified per-object truth? Refused unless its named gates pass; `--build` makes the oracle caches |
| `design/prior_vs_oracle.py` | Is `LocusPriors`' gDNA count right? `P − O` is calibration's error, `O − S` the pooled share's (`ISSUES: the-pooled-q-in-the-gdna-count`) |
| `design/quant_accuracy.py` | How accurate is the tool end to end, and what is a perfect prior worth? Read in a pinned A/B pair; `--report … --markdown` renders the per-scenario release report |
| `design/build_scan_cache.py` | Scan once, calibrate many times. ⛔ The key hashes the deposit rule, not `resolve.cpp`: after a change to which fragments are offered, `--force` |
| `design/rename_identity.py` | Is a rename, refactor or speed-up numerically a no-op? `--freeze` once, `--check` after every stage, `--bam` on a real library; the reference is frozen, never rolling |
| `sim/panel.py` | Build, simulate, cache, score and report a panel from one YAML; `status` names the next stage. ⛔ `cache` builds both caches; the oracle one is the truth every scorer reads |
| `sim/build_test_reference.py` | Renders the test chromosome from its one hand-edited YAML, `test_chr.yaml` (`docs/TESTING.md` §0a) |
| `sim/build_suite_reference.py` · `design_suite_probes.py` · `simulate_reads.py` | The panel's substrate; `panel.py build` drives the last two, and the reference carve is manual |
| `sim/configs/gdna_ladder.yaml` | The 16-condition panel the tool is ranked on (`docs/TESTING.md` §0) |
| `sim/configs/flgap_rna_long.yaml` · `flgap_rna_short.yaml` | The fl-gap side panel, re-simulated 2026-10-03 with the ladder's nascent block, so the length gap is its only difference from the ladder. ⛔ Never a ladder rung; read both arms, per stratum, beside the ladder |

## CLI

```bash
rigel index --fasta genome.fa --gtf annotation.gtf -o index/
rigel quant --bam sample.bam --index index/ -o results/
rigel sim --config scenario.yaml -o out/
rigel export results/ -f tsv
rigel report results/ -o report.html
```

Input BAM must be name-sorted with the `NH` tag.
