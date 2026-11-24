# CLAUDE.md

A session primer: what Rigel is, where it stands, where a change goes, and how work is done here. Git holds the history
and the eight docs below hold the detail; open problems are GitHub issues.

## What Rigel is

Rigel is a Bayesian RNA-seq transcript quantifier that separates RNA from genomic-DNA contamination. The pipeline:
1. A single-pass C++ BAM scanner tallies fragments per genomic region.
2. A **calibration** stage splits each region's unspliced count into gDNA and RNA.
3. A per-locus EM assigns RNA to transcripts.

The PyPI package is `rigel-rnaseq`; the import and the CLI are `rigel`. The version on disk is 0.7.1
(`pyproject.toml`); the target is **0.8.0**.

## Axiom 0: RNA is RNA

There are three populations: `gDNA`, `RNA+` and `RNA−`. "Mature" and "nascent" are not populations. The only
distinction is whether a fragment **is spliced** (certified RNA) or **is not**, and the unspliced ones are the whole
deconvolution problem. The population set at a slot follows from the annotation:

```
T(slot) = {gDNA} ∪ {RNA+ if statics.free_pos[slot]} ∪ {RNA− if statics.free_neg[slot]}
```

A population set with more than three members, or "mature/nascent" in a solver question, violates the axiom. An output
reports RNA as **annotated** (reference transcripts) or **nascent** (Rigel's synthetic spans). Nascent RNA is modelled
for robustness, so a number measured at the panel's nascent share is a stress reading, never a design driver
(`docs/DESIGN.md` §0).

## Where things stand

- `src/`'s numbers are 1f97ff96's (the junction cap and the share over one tally) plus the **capture weight**, landed
  2026-11-15: the capture reading is raised to the share of the library's reads whose own gDNA would show. At 5 % gDNA
  and above no in-scope whole-genome row rises by more than about 1 %; the 0.1 % hole closes on genes; STAR's zero-gDNA
  rows fall, the switch row then about four times the no-capture line (`docs/ROADMAP.md`).
- The **scan phase** landed 2026-11-18: the splice blacklist is retired, and a multimapper's candidates are pruned
  within each place. The simulator's rows are bit-identical; on the STAR rung the stranded capture-ON rows move,
  through the one locus that holds most of the row (issue #78).
- The **intron abundance model** landed 2026-11-19: every intron's count is read against the library's own intergenic
  gDNA law, with a departure prior a stranded library measures from its strands. Transfer now beats silent on
  unstranded capture OFF; the stranded strata hold; VCaP is unchanged; sparse unstranded libraries book about 1 % of
  their gDNA as nascent (issue #82).
- The **unclamped capture reader** landed 2026-11-20: an object reads its own level, not its enriched basin's peak.
  Transcripts rise 4.6 % on the simulator and 9.7 % on STAR, between the isoforms of genes with many (issue #83).
- The simulator's **probes bind by their own GC** (2026-11-21), and the capture panels are regenerated with them.
- **A locked object reads its strands** (2026-11-21, issue #1): a stranded library's one-stranded fragments on an
  object the annotation admits no RNA at are no longer booked gDNA whole. Only stranded STAR rows move: transcript
  error falls on all six, most at no gDNA (431,143 → 287,469).
- **The cuts' prices** landed 2026-11-22 (issue #76): an RNA placement that leaves an exon across a junction is
  priced at the junction's crossing level, read from gDNA beside it, not by the half-intron gDNA object at the cut.
  Under stranded capture at 5 % gDNA and above the annotated over-call per locus falls 13–41 % on five of six rows;
  isoforms pay on rows whose capture reading is false or thin (STAR unstranded capture OFF 5 %, transcripts +9.8 %).
- **The index's one rule** landed 2026-11-23: nascent entities whose ends lie within `--end-tolerance` are one, and a
  single-exon transcript is its own nascent RNA (`docs/DESIGN.md` §2.2). The ladder's in-scope transcript error moves
  −4,117 fragments on the same reads. With it the E-step holds each distinct likelihood row once (issue #84). This
  tree is the baseline.
- The **design phase on fused loci and the EM's cap** (issues #79, #34, #78) closed for 0.8.0 on 2026-11-24: VBEM
  ships. MAP with removal by evidence is issue #86, annotated aliases issue #87; the receipts and their map are in
  `$RIGEL_RUNS/perf/2026-11-23_em_cap/`.
- **Every error is read per object** (the owner, 2026-11-21): Σ|estimate − truth| over objects, never a library net
  (`docs/SUCCESS.md`). Bar 5 holds on every row it covers, transcripts at most 0.72 of 0.7.1's and genes 0.55.
  Bar 1 holds on every row but one: stranded capture ON at 5 % gDNA, 9.5 % of the library's gDNA above chance.
- The validation of ec996bba beat 0.7.1 in every in-scope stratum and stated eight failures, mostly about capture
  (`docs/ROADMAP.md`).
- **No preview ships.** 0.8.0 waits for every `release-blocker` issue in the milestone; the milestone's other issues
  ship stated in the MANUAL.
- **The order of work** (`docs/ROADMAP.md`, Next): performance → code review → the release, with CI alongside. No session has a pass bar: each improves the current state, and what is left ships
  stated. One `src/` writer at a time, in the main checkout.
- The **capture reader** is developed together with the owner. Its working memory is `docs/dev/CAPTURE_IDEAS.md`.
- Everything retired at the 2026-11-03 reset is in `$RIGEL_RUNS/archive/cleanup_2026-11-03/`.

## The 0.8.0 scope

The ruling is in `docs/DESIGN.md` §0b.
- **Two primary numbers.** Neither stands in for the other:
  - the transcript table against per-transcript truth (`quant_accuracy.py`);
  - the calibration result against an oracle calibration (`calibration_vs_oracle.py`, `prior_vs_oracle.py`).
- **The pinned A/B.** Run with `--set scan.total_threads=1` and `--set em.assignment_mode=fractional`, so the pair is
  exactly reproducible. Read genes and pools beside the transcript table.
- **Strata.**
  - In scope: unstranded × capture-OFF, stranded × capture-OFF, stranded × capture-ON.
  - unstranded × capture-ON is **deferred**: reported on every benchmark, but not a release requirement.
  - Every score is read per stratum, never pooled.
- **Fragment length.** The EM's two length tables are one shared law unless the data show a difference capture cannot
  make. The ladder gives gDNA and RNA equal lengths on purpose; length mechanisms are judged on the fl-gap panels, in
  both directions.

## The docs

| | |
|---|---|
| `docs/SUCCESS.md` | **Start here.** How performance is judged, and the instruments in run order |
| `docs/ROADMAP.md` | The state, and the ordered next steps |
| `docs/TRAPS.md` | Mistakes already made, as named rules |
| `docs/EQUATIONS.md` | The derivations the code depends on |
| `docs/DESIGN.md` | What is built and the rulings behind it. §0 is the binding vocabulary |
| `docs/TESTING.md` | The panels, the test chromosome (§0a), and the gates |
| `docs/MANUAL.md`, `docs/PUBLISHING.md` | The user's manual and the release procedure |

Every open problem is a GitHub issue (https://github.com/mkiyer/rigel/issues). The `0.8.0` milestone
(https://github.com/mkiyer/rigel/milestone/1) is what the release must fix or state.

- **One home per finding.** When a finding settles, move it to its home and delete it where it was: an open problem to
  a GitHub issue, a lesson to TRAPS, a ruling to DESIGN, a derivation to EQUATIONS.
- **Cite an open problem as `issue #N`**, and a TRAPS rule by its name, e.g. `TRAPS: name`. A fix's commit says
  `Fixes #N`; a refusal closes as not planned, with the number that killed it.
- **`docs/dev/` is the sandbox.** Nothing there is authoritative and nothing may cite into it.
- **The source cites no doc and no issue.** `tests/test_docs_boundary.py` gates these rules.

## Where a change goes: the calibration layering

An import may point down a layer or sideways, never up. `src/rigel/calibration/_layers.py` is authoritative and
`tests/calibration/test_layering.py` enforces it.

| if the change is about… | it goes in |
|---|---|
| what a fragment tally means | **1 · payload view**: `splice_graph` `substrate` `region_arrays` |
| how many places a fragment could have sat | **2 · opportunity**: `effective_length` `capture_eff_length` `sj_opportunity` `gdna_opportunity` `gdna_density` `fl` |
| one slot's own numbers and ψ | **3 · geometry + the per-slot solve**: `region_geometry` `simplex_logodds` |
| which strand a fragment came from | **4 · strand**: `strand_balance` |
| how dense a component is, and the priors | **5 · density and prior**: `landscape` `capture_efficiency` |
| what one neighbour tells another | **6 · the solve**: `sweep` + `diagnostics` + `messages/` |
| turning the solve into a result | **7 · assemble**: `calibrate` `priors` `result` `track` |

The sweep is one native call (`native.solve_blocks`). `message_policy = "transfer"` ships; `"silent"` is the measured
floor.

## Working rules

- **Prototype outside `src/`, A/B against what ships, then land.** Derive first. Prototype outside the main tree (C++
  in a worktree, Python beside it). A/B against what ships on the same conditions. Change one mechanism at a time.
- **No magic numbers.** A fixed number lives in `config.CONSTANTS` with the reason for its value. A run setting lives
  in `PipelineConfig` or `IndexConfig`. `tests/test_constants.py` refuses a bare module constant.
- **A falsification test first**, verified failing. Then break the fixed code and watch each gate fire.
- **Score against truth.** Vary one thing per experiment. Re-record the baseline from the current tree in the same
  session, and score against the oracle BAM's read names, never against the previous run.
- **No legacy, no backwards compatibility, no speculative code.** No version suffixes in file names, and no Greek
  letters in identifiers.
- **Real data is a test input, never a design input.** Read a timing only from back-to-back A/B pairs.
- **A speed-up is held to a tolerance, not to the bit** (the owner, 2026-11-22): the payload, the calibration result
  and the tables' lengths and priors within `rename_identity.py --check --rtol 1e-9` of the tree it starts from, and
  the EM's counts unmoved beyond noise against truth on the ladder (`panel.py score`, A/B per stratum): the EM's split of
  near-tied isoforms moves with any last-bit change, main's own scan thread count included. Write the clean, fast
  formulation; no contortion to keep the last bit.
- **Renames and refactors:** run `rename_identity.py --check` after each stage.
- **The owner drives commits.** Do not commit unless asked.

## Build, test, lint

Everything runs in the `rigel` conda env, which holds htslib and the compilers.

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
export OMP_NUM_THREADS=1                        # always, when benchmarking or comparing runs

pip install --no-build-isolation -e ".[dev]"    # rebuild after ANY src/rigel/native/ change
python -m pytest tests/ -q
python -m pytest tests/ -q 2>&1 | grep '^FAILED' | sed 's/::.*//' | sort | uniq -c   # the failure set
ruff check src/ tests/ scripts/ && ruff format src/ tests/   # never format scripts/
python scripts/design/preflight.py --full       # every instrument's --self-test
```

**The baseline: 0 failed / 3,811 passed / 0 skipped, 3,811 collected** (re-derived 2026-11-24, the review's
stage 6a). Any failure is a regression.
- Re-derive the count; never edit it by hand. Only adding or removing a test moves it.
- A golden update (`--update-golden`) is where a regression gets laundered into "intended". Read the diff and check the
  truth-scored instruments first.
- After a deletion, a rename or a default flip, run `preflight.py --full` as well as the suite.

## Tooling under `scripts/`

Each instrument's module docstring and `--help` are its manual. `docs/SUCCESS.md` gives the run order.
`tests/test_scripts_index.py` holds the `design/` rows below against the disk, and `scripts/README.md` indexes
`scripts/profiling/` and the `_`-prefixed helper modules of `scripts/design/` and `scripts/sim/`.

| | the question it answers |
|---|---|
| `design/preflight.py` | Can this session run and regenerate everything? `--full` adds every `--self-test` |
| `design/quant_accuracy.py` | How accurate is the tool end to end, against truth? Read in a pinned A/B pair. `--arm no_capture` is the line a capture reader must beat; `--report … --markdown` renders the release report; `--bam` runs on an aligner's input (the STAR rung) |
| `design/paired_score.py` | Does one arm forecast each fragment's source better than another? The Brier score of P(gDNA) and its CORP miscalibration on identical read names, per stratum, closed on `quant_accuracy.py`'s pools. `--record` one arm per process (a second build under its own interpreter; `--arm no_capture` for the line a capture-ON row is read beside), then `--pair` the two records |
| `design/calibration_vs_oracle.py` | Is the calibration result right against an oracle calibration? `--set` prices a config value on both arms. It cannot see the capture-contracted length: read `ruler_vs_truth.py` beside every capture-ON arm |
| `design/ruler_vs_truth.py` | Is the EM's capture-contracted length right against the simulator's own yield? `--scale` puts every class on one scale; `--module` runs prototype rulers |
| `design/prior_vs_oracle.py` | Is the locus prior's gDNA count right? `P − O` is calibration's error, `O − S` the pooled share's |
| `design/policy_benchmark.py` | How does each message policy score per condition? `--by-class` shows where the error sits. Per stratum, never pooled |
| `design/calibration_oracle.py` | What is the certified per-object truth? `--build` makes the oracle caches |
| `design/build_scan_cache.py` | Scan once, calibrate many times. After a change to which fragments are offered, `--force` |
| `design/rename_identity.py` | Is a refactor or speed-up numerically the same? Bit-identical, or within `--rtol`, against frozen references in `$RIGEL_RUNS/arms/` |
| `sim/panel.py` | Build, simulate, cache, align (the STAR rung), score and report a panel from one YAML. `status` names the next stage |
| `sim/build_test_reference.py` | Renders the test chromosome from `test_chr.yaml` (`docs/TESTING.md` §0a) |
| `sim/configs/gdna_ladder.yaml` | The 16-condition panel the tool is ranked on |
| `sim/configs/flgap_rna_long.yaml`, `flgap_rna_short.yaml` | The fl-gap side panels; read both, per stratum, beside the ladder |
| `sim/configs/gdna_ladder_lowg.yaml` | The low-gDNA side panel: 0.1 % and 1 % gDNA, stranded, beside the ladder's `g00` and `g05` |

## CLI

```bash
rigel index --fasta genome.fa --gtf annotation.gtf -o index/
rigel quant --bam sample.bam --index index/ -o results/
rigel sim --config scenario.yaml -o out/
rigel export results/ -f tsv
rigel report results/ -o report.html
```

The input BAM must be paired-end and name-sorted or collated. A read pair's multiplicity is its usable alignments. Align
STAR input with `--sjdbScore 0 --outFilterMultimapScoreRange 2` (`docs/MANUAL.md`, input requirements).
