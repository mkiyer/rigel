# CLAUDE.md

Guidance for Claude Code working in this repository.

## What Rigel is

A Bayesian RNA-seq transcript quantifier that jointly models mRNA, nascent RNA and genomic DNA
contamination. A single-pass C++ BAM scanner tallies fragments, a **calibration** stage deconvolves the
library into gDNA vs RNA, and a per-locus EM solver assigns RNA to transcripts. PyPI package
`rigel-rnaseq`; the import and CLI are `rigel`.

## ⭐ The docs — read them in this order

There are six. They do not overlap, and none of them is a changelog.

| doc | what it is |
|---|---|
| ⭐⭐ **`docs/SUCCESS.md`** | **START HERE.** How performance is measured in this phase: **Stage A** (the accumulator — faithful, unbiased, sufficient?) and **Stage B** (calibration init + pass-0, against an oracle). Says what "done" means for each, and why the end-to-end number is a thermometer rather than a target |
| ⭐⭐ **`docs/ROADMAP.md`** | Where the tool is, in six measured numbers, and **what to do next in priority order** |
| ⭐ **`docs/TRAPS.md`** | **Mistakes this project has already made.** Read before designing anything. Lessons, not measurements — several of these were made twice |
| **`docs/EQUATIONS.md`** | The derivations the code depends on. Deposit rule, opportunity functions, the 2×2, strand, overdispersion, damping |
| **`docs/DESIGN.md`** | What is built and the rulings behind it. Settled — do not re-litigate. ⭐⭐ **§0 is the binding VOCABULARY** — NODE, EDGE, step, the mass pin, structurally pure-gDNA object — and it lists the banned synonyms. Read it before writing a comment or a docstring |
| **`docs/TESTING.md`** | The simulated panels, how to build them, the simulator's own gates, and what the suite can and cannot judge. ⭐⭐ **§0b is the TOY HARNESS** — how to define a test transcript and calibrate it in under a second; read it before hand-rolling any scenario |

Reference rather than design: `docs/MANUAL.md` (user-facing), `docs/PUBLISHING.md`.

⚠ **`docs/calibration/` holds EPHEMERAL working docs** — currently just `SESSION_HANDOFF.md` — used
to pass one task between the owner and a session. ⛔ A per-issue REGISTER used to live here too and was
deleted on 2026-08-05: it had become a second copy of `ROADMAP.md`'s ranked list, which is `TRAPS.md` A11's
"two homes for one predicate" in documentation form. One ranked list, in `ROADMAP.md`. They are **not** permanent fixtures and are **not** in
the six above. ⛔ Read the handoff if one is present and you are picking up its task; **delete it when that
task lands**, promoting anything worth keeping into `ROADMAP.md` (a current number), `TRAPS.md` (a lesson)
or `DESIGN.md` (a settled decision).

⛔ **THE SOURCE DOES NOT REFERENCE THE DOCS, AND MUST NOT START.** Docs evolve and rot — 73 % of the
citations that used to be in the source pointed at documents that had already been deleted. A docstring may
cite a **test** or the executable specification (`tests/native/_accumulator_reference.py`), because code
cannot rot silently; it may not cite a document.

⛔ **`tests/native/_accumulator_reference.py` is the executable specification** for the accumulator. The
C++ is gated on byte-identity to it; where it and a document disagree, it wins.

## Working rules

- **No magic numbers.** Stop and discuss before adding any constant, heuristic or tunable. Every divisor
  must be derived from the deposit rule and unit-tested against brute-force enumeration.
- **A falsification test first, verified failing — then break the fixed code and watch each gate fire.**
  The second half is not optional; it has found holes in already-green gates repeatedly.
- ⭐⭐ **Measure the CEILING before building the correction.** Hand the consumer the *exact* answer for one
  channel and see what perfecting it is worth. `calibration_truth_ab.py --ceiling`. This has re-ranked the
  project twice, and it is also how you learn a phase is finished.
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
- **Profile on real cfRNA, never a small synthetic suite** — a toy ranks hotspots backwards.
- **The owner drives commits.** Do not commit unless asked.

## Build, test, lint

> **Every build, test and lint command must run inside the activated `rigel` conda environment** — it holds
> htslib and the compilers, and the C++ build finds htslib via `$CONDA_PREFIX`.

⛔ **Run the suite as `python -m pytest`, never bare `pytest`.** Bare `pytest` does not put the repo root on
`sys.path`, so an intra-repo import raises and the suite reads one extra failure.

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel

pip install --no-build-isolation -e ".[dev]"   # rebuild after ANY src/rigel/native/ change
python -m pytest tests/ -q                     # baseline: 22 fail / 3 xfail — 21 goldens + the paralog row
python -m pytest tests/ --update-golden        # regenerate tests/golden/ after intended output changes
ruff check src/ tests/ scripts/ && ruff format src/ tests/   # ⚠ NEVER format scripts/
```

⭐ **"22 failures, 21 goldens and the paralog row" is the standing baseline.** A 23rd failure, or any other
non-golden name in the list, is a regression. See `docs/TESTING.md` §6.
⭐ **And 3 xfails, one of them STRICT and load-bearing**: `test_toy_harness`'s intron-composition gate is the
project's detector for the level defect that the splice-flux reframe un-masked. It must go green — not be
widened — when the joint arm lands (`SESSION_HANDOFF.md` §1, §4).

Always set `OMP_NUM_THREADS=1` when benchmarking or comparing runs.

## Tooling under `scripts/`

`docs/SUCCESS.md` lists these in the order you would run them. ⚠ Everything not listed here was deleted as
unrunnable; if a script is not in this table, treat it as not existing.
⛔ **That claim is currently FALSE for eight files** and the drift is recorded rather than silently
tolerated: `bam_spans.py`, `implicit_splice_census.py`, `inv_L_limits.py`, `length_sieve.py`,
`partition_test.py`, `reference_on_real_data.py`, `sigma_inv_L.py`, `spanning.py` are on disk and unlisted.
They pre-date 2026-08-05, none was run this campaign, and none should be trusted without re-running it —
either promote a row or delete the file, but do not assume the table is complete until that is done.

| | |
|---|---|
| **⭐⭐ the toy harness** | |
| `design/toy_panel.py` | ⭐⭐ **one toy spec × EVERY cached condition × an RNA-density ladder, scored per object** — what you run once a structure is the *target* rather than a probe. Prior-free pass-0 by default; the RNA density is a multiple of each donor's OWN gDNA density. ⛔ Reports which object carries the error, whether the messages helped or hurt it, and which objects are CONFIDENTLY wrong. ⚠ 13 s per condition; shard with `--conditions` |
| `design/verify_toy_substrate.py` | ⭐⭐⭐ **IS THE INPUT CORRECT? — no solver runs.** Every accumulator bank re-derived from the simulator's per-fragment TRUTH by an independent implementation, plus the splice combinatorics and the length marginal. ⛔ Run this on any new toy spec BEFORE reading a solver number off it |
| `design/verify_capture.py` | ⭐⭐ **hybrid capture, probes ON vs OFF on identical geometry** — the gDNA landscape, the length selection and the junction depletion, each gated on the direction the knobs predict |
| `design/toy_trace_error.py` | ⭐⭐⭐ **the full error trace for ONE scenario** — every node and edge ranked by error MASS, the four init sources, the relay hop by hop with the licence state, then a ψ channel ABLATION (fidelity-checked to 1e-6) and ⭐ a RELAY-level arm that withholds the composition licence in both twins at once. ⛔ This is what "dissect the error" means |
| `design/zero_controls.py` | ⭐⭐⭐ **THE TWO ZERO CONTROLS — and the owner requires them on every experiment.** ZERO RNA (silent genes, truth `f_g = 1.000` everywhere) and ZERO gDNA (the `g00` donor, truth `0.000`). The truth is a CONSTANT, so every deviation is a false positive with nothing to cancel it. ⭐ The zero-RNA arm is the **biologically dominant** case: most annotated transcripts are OFF. ⛔ Prints the three-rung ladder (strand → self-solve → final), what the MESSAGES delivered as an implied `f_g`, and flags any EMPTY object, because a zero arm can be degenerate and then it tests nothing (`TRAPS.md` A14) |
| `design/reframe_walk.py` | ⭐⭐⭐ **THE REFRAME, WALKED — every count, both flank totals, and every hop in BOTH directions, on one two-exon transcript.** Not a summary statistic: per hop it prints `r` as used, `r` as the predecessor's single total would have given it, what TRUTH says the same ratio is, and the true gDNA-density ratio beside it — so the include/exclude decision at each `(EDGE, side)` is a number. ⭐ Three gDNA arms (zero / high / very high), capture-OFF × unstranded, HIGH RNA so nothing is sparse. ⛔ The geometry is rebuilt and GATED against the solver's published frames to 1e-12 |
| `design/toy_dissect.py` | ⭐⭐ **one scenario opened to every slot and every channel** — the solver's own ladder (strand-only → self-solve → final) beside `cm_g` / `c_tau` / `cg`, so a dead channel is visible as a zero rather than inferred. ⛔ Step 3 of the debug loop, for a toy |
| `design/toy_ceiling.py` | ⭐⭐ **the RE-SOLVE ceiling** — hand one object class a different own belief and re-solve the whole chain, six arms sharing one simulation. ⛔ This is what `TRAPS.md` B17 demands instead of a substitution, and `--arms base noop` is its own falsification (must be byte-identical). Its docstring carries what it measured |
| `design/ladder_arm_ab.py` | ⭐⭐ **the same override on the REAL 36-condition ladder**, scored by `solvability_audit.py`. ⛔⛔ **Run this before writing a mechanism into `src/`** — twice now a toy-positive change has been panel-negative (`TRAPS.md` B18). 40 s per condition with the oracle cache |
| `design/length_ceiling.py` | ⭐ **what a PERFECT length model is worth, on the ladder, ONE pmf at a time.** ⛔ Pricing the pair together hid a 14× split between them (`TRAPS.md` B21). Reports the pass-0 solvable yardstick beside the mass-weighted one, because they disagree |
| `design/toy_harness.py` | ⭐⭐ **a mini chromosome you define, calibrated in 0.1–5 s** (`docs/TESTING.md` §0b), with every object's answer beside per-object truth. The library-level priors a toy cannot fit are harvested from a real cached condition and INJECTED; the gDNA depth is DERIVED to match that donor, never set by hand. `--list` for the spec ladder. ⭐ Reach for this FIRST when a mechanism needs isolating — it found C1's mechanism candidate in one sweep |
| **the substrate** | |
| `design/simulator_gates.py` | ⛔ the simulator's own gates G-S1…G-S6, scored on per-fragment truth. Run before trusting the panel |
| `design/suite_resolves.py` | ⛔ can the suite resolve the axis you are changing? Run before quoting any suite number |
| `design/index_census.py` | re-derive an index's census — never quote a stored table, run this |
| `design/verify_index_rebuild.py` | nodes byte-identical, edges only in contiguous reach |
| **Stage A — the accumulator** | |
| `design/fl_anchor_gap.py` | the anchor and both length models vs truth. `--drain` measures before and after the second pass |
| `design/gdna_pool_census.py` | ⭐ the four gDNA pools, each against its own opportunity and against truth |
| `design/second_pass_accuracy.py` | the second pass scored PER FRAGMENT against the oracle BAM's read names |
| `design/observable_efficiency.py` | what fraction of the length information a storage choice keeps |
| `design/node_density_derivation.py` | the reciprocal-opportunity theorem, T0–T6, each perturbed |
| **the scoping number** | |
| `design/calibration_truth_ab.py` | ⭐⭐ the deliverable against truth, and `--ceiling` — what perfecting each length model is *worth* |
| **Stage B — calibration** | |
| `sim/configs/gdna_ladder.yaml` | ⭐ **the Stage-B panel** — 36 conditions, EQUAL fragment lengths, gDNA 0 → 98 % at a fixed 10 M total. An ORACLE CACHE ships beside it at `ladder/oracle_cache/` and stays valid across every calibration change. `docs/TESTING.md` §0 explains both panels |
| `design/solvability_audit.py` | ⭐⭐ **START HERE FOR PASS-0.** Which objects are SOLVABLE, which of those are solved wrong, and — the one that matters — which are **confidently** wrong, since a tight-variance wrong value outvotes correct neighbours and anchors the prior. ⛔ Excludes the undetermined population: in pass-0 an object with no own evidence reporting `f_g ≈ ½` at zero precision is *correct*. Carries the ablation ladder (strand → local → final) and a calibration curve that needs no threshold. `--suite` alone gives the whole panel |
| `design/pass0_vs_oracle.py` | T (the origin-split payload) vs P (`calib_refit_iters=0`) vs two levered ceilings, per object and per class. ⛔⛔ **ITS HEADLINE MASS-WEIGHTED ERROR IS THE WRONG YARDSTICK FOR PASS-0** — it scores every object with mass, so honest ignorance reads as error and buried a 0.0456 answer inside a 0.3150 one. Use it for the T/C/P decomposition and the oracle plumbing; use `solvability_audit.py` to judge pass-0. ⛔ undrained on every arm. `--oracle-cache` makes repeat runs cheap: the oracle depends on the accumulator and index, never on calibration, so one cache serves a whole debugging campaign |
| `design/worst_objects.py` | ⭐⭐ **step 3 of the debug loop** — one condition dissected to individual nodes/edges, ranked by error **MASS**. Read the CONCENTRATION curve first (concentrated ⇒ a mechanism exists; diffuse ⇒ a systematic bias). `fg_loc` vs `pred_fg` separates a wrong local solve from wrong messages |
| **diagnostics** | |
| `design/composition_evidence_census.py` | how much library mass reaches the solver with NO composition evidence. `--inject-kappa 0.5` is its falsification handle |
| `design/held_flux_census.py` | how often a held candidate has ZERO flux evidence, by cause |
| `design/prior_units_check.py` | the EM prior in fragment units vs the old incidence sum |
| `design/length_likelihood_ab.py` | the A/B for the per-node length likelihood (currently gated OFF) |
| **plumbing** | |
| `design/build_scan_cache.py` | scan once, calibrate many times. ⚠ re-run after any accumulator change — the payload schema digest invalidates every cache, by design |
| `design/native_parity_on_real_data.py` | the native-parity gate on real cfRNA at full scale |
| `design/scan_profile.py` | ns/fragment, regressed over several BAMs |
| `sim/build_suite_reference.py` · `design_suite_probes.py` · `simulate_reads.py` | build the panel |
| `sim/evaluate_suite.py` | net fragment flow (`rigel.sim.analysis`) |

## CLI

```bash
rigel index --fasta genome.fa --gtf annotation.gtf -o index/
rigel quant --bam sample.bam --index index/ -o results/
rigel sim --config scenario.yaml -o out/
rigel export results/ -f tsv
rigel report results/ -o report.html
```

Input BAM must be name-sorted with the `NH` tag.
