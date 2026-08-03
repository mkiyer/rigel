# CLAUDE.md

Guidance for Claude Code working in this repository.

## What Rigel is

A Bayesian RNA-seq transcript quantifier that jointly models mRNA, nascent RNA and genomic DNA
contamination. A single-pass C++ BAM scanner tallies fragments, a **calibration** stage deconvolves the
library into gDNA vs RNA, and a per-locus EM solver assigns RNA to transcripts. PyPI package
`rigel-rnaseq`; the import and CLI are `rigel`.

## ⭐ Current direction — read this first

⭐⭐ **FRAGMENT LENGTH IS DONE (2026-08-03), AND THE CRITICAL PATH IS NOW C3 — junction opportunity.**
Length was the blocker on this whole area. It is finished: **one** definition (C0–C2), an **accurate** one
(C2.6), and an **unbiased** one — the anchor's error against the simulator's own truth is **+0.00 % mean /
+0.02 % sd** on the zero-gDNA falsification condition, from −1.61 % / −1.48 %.

⛔ **BUT THE DELIVERABLE GOT WORSE, AND THAT IS WHY C3 IS RANK 0** (`LEDGER.md` **B4**). The library gDNA
fraction's mean |error| on the four contaminated pilot conditions went **0.0381 → 0.0472 (+23.9 %)** once
the side buffer drains. ~55 % of that is the drain being *right* (truth says **99.8 %** of held fragments
are RNA, so depositing them genuinely lowers the gDNA fraction, and the estimate was already too low);
~45 % is the fragment-length **pools** moving. ⭐ **`calibrate` fits from the two pure pools, not from the
anchor** — and `RNA_SPLICED` went **+2.4 % → +6.2 %** against truth, because the drain feeds long
junction-using fragments into a pool selected on "used an annotated junction", which is itself **+3.8 %**
longer than the library. C3 is the correction for exactly that.

**Read them in this order:**

| doc | what it is |
|---|---|
| **`docs/TODO.md`** | ⭐⭐ **START HERE — the one ranked list, and rank 0 is C3.** Everything landed is struck through with its measurement |
| **`docs/LEDGER.md`, entries C0 → B4** | ⭐ **THE MOST RECENT WORK, newest last.** The fragment-length track (C0–C2.6), the two-pass structure (S1, S2, S2.1), the second pass in full (P0–P4.2), and ⭐ **B4 — the composition baseline that sets the priority** |
| **`docs/JUNCTION_OPPORTUNITY.md`** | ⭐ **C3's formula, derived and proven** over 48,648 exhaustive configurations. ⛔ **EVERY number in §3 is STALE** — scored against the contaminated anchor and the pre-drain pool. §1's derivation is untouched |
| **`docs/SPEC_SECOND_PASS.md`** | ⭐ the second pass, in full: the score, the draw, the drain. **P0–P4 are all CLOSED**; §8's D-1/D-2/D-3/D-5/D-6 are all decided and **D-4 is the one left** |
| **`docs/SOLVER_OBSERVABLES_PLAN.md`** | ⚠ its P0/P1/P2 are **not** the second pass's. ⭐ **P2 — the fragment-length likelihood in the per-node solve — is built and gated OFF, blocked on the FL pools**, i.e. on C3. Mechanism proven: blind mass 100 % → 0 % |
| `docs/PLAN_TWO_PASS.md` | ⚠ **HISTORY** — why the gap-path and junction-opportunity problems are ONE problem. S1–S3 all landed; its §5 is superseded by `SPEC_SECOND_PASS.md` and its §2.4 numbers by B4 |
| `docs/FRAGMENT_LENGTH_AUDIT.md` | ⚠ **HISTORY** — how the critical path got here. THREE definitions of fragment length were live at once. C0/C1/C2/C2.6 all landed; its "C3 is next" is now correct again |
| `docs/SPEC_GAP_PATHS.md` | ⚠ **HISTORY** — the enumeration rule pass 1 implements, landed as S1. Supersedes `SPEC_GAP_INTRONS.md`, which is the C2.6 record |
| `docs/S5_DESIGN_LOG.md` | ⚠ **HISTORY** — S5 is finished. Kept for §1's accumulator derivations and §3's observable measurements, still the reference |
| `docs/IMPLEMENTATION_PLAN.md` §0 | the live handoff for everything not on the critical path |
| `docs/NODE_DENSITY_DERIVATION.md` | why the deposit weight is 1/opportunity, and what each stored channel buys |
| `docs/ACCUMULATOR_DESIGN.md` | the accumulator's design. ⚠ §8's purity claim for the two *splash* pools is **unverified** — see `TODO.md` |
| `docs/CARRY_FORWARD.md` | ⭐ **§3 traps then §2 equations** — the most-used reference in the project |
| `docs/BENCHMARK_SUITE.md` | the suite: how to build it, and **what it can and cannot judge** |
| `docs/LEDGER_ARCHIVE.md` | older ledger entries |

### ⭐ The four numbers that describe the tool today

| | measured on | |
|---|---|---|
| fragment length, the **anchor** | 8 pilot conditions vs truth | ⭐ **+0.00 % mean / +0.02 % sd** |
| the second pass's **per-fragment** accuracy | 171,534 held fragments vs per-fragment truth | ⭐ **90.5 % exact**, mean error **+0.12 bp** |
| fragments above the library's true longest molecule | all 8 conditions | ⭐ **0** |
| ⛔ the **deliverable** — library gDNA fraction | 4 contaminated conditions vs truth | ⛔ mean \|error\| **0.0472** (worst row **0.3704 against 0.5**) |

⚠ **The last row is the one to care about.** The first three are upstream plumbing; the fourth is the
product, and it is what C3 exists to move.

Reference rather than design: `BENCHMARKING.md` (how to evaluate — net fragment flow), `MANUAL.md`,
`PUBLISHING.md`, `docs/testing/testing_plan.md` (the owner's plan for the cached-substrate harness).

### What is settled, and must not be re-litigated

⭐ **ONE DEFINITION OF FRAGMENT LENGTH.** `L` = genomic span minus cut introns, proven (C0); the
accumulator bins every deposited fragment by it (C1); every consumer reads it (C2). The scanner's rival
histogram, `FragmentLengthModels` and the transcript-space definition are **deleted**, and
`build_fl_models` takes the payload and nothing else, so a mixed-frame call is unrepresentable.
⚠ `FragmentLengthModel` **singular** is the scorer and stays.

⭐ **A GAP INTRON IS CUT ON EVERY FRAGMENT** (C2.6), not only unspliced ones, with the gaps the CIGAR
already explained excluded by **exact `(start, end)` equality** — ⛔ overlap would let a *different* nearby
intron answer for one and make `L` too short.

⭐ **THE ACCUMULATOR ARBITRATES** (S1). A fragment arrives with its hypothesis *set*; exactly one survivor
deposits, two or more are held WHOLE in the side buffer. ⭐ **And the second pass drains it** (P0–P4.2):
score from pass-1 evidence alone, one multinomial draw, re-deposit through the **same** `deposit`.

⭐ **THE POOL IS KEYED ON DETERMINACY, NOT PROVENANCE** (D1, deleted at S1). A fragment enters a length
pool when exactly ONE hypothesis survived, however it got there. ⛔ **A purity filter on a length pool is a
length filter**: keying on provenance measured **−9.58 % / −22.46 %** against truth because the excluded
fragments are the LONG ones.

⭐ **THE SCORE'S COMBINATION RULE** (P4.1/P4.2). `score = rho x f(L) x s`, normalised within one fragment's
candidate set, factors applied in order of the evidence behind them and skipped when flat-zero among the
survivors. ⛔ **An all-zero factor is uninformative, not decisive** — it used to annihilate the other two
and collapse the record to a coin toss. ⭐ **gDNA's strand term is 0.5** (double-stranded, no sense
direction); a fitted mixture marginal was measured and **refuted** — it destroys 78 % of the orientation
signal, because the discrimination is `(1−p)/p` and any *constant* for ∅ cancels out of it.

⚠ **Three things the baseline says**, before anything is built on it:
1. ⛔ **A7's "11.0 % gDNA over-call" DID NOT SURVIVE MEASUREMENT.** The A/B is done (`LEDGER.md` S5.g-2):
   turning the contiguous-edge RNA taper on moves the library gDNA fraction by **≤ 0.0002**. The 11.0 %
   was a *bp-weighted* geometric mean; the estimator is *fragment*-weighted, and **89 % of edge mass sits
   on lines the taper does not touch**. `CARRY_FORWARD.md` §1 fact 6 is corrected.
2. ✅ **THE κ "MIRROR" WAS NOT A DEFECT — closed 2026-08-02, and it had been filed twice.** Two different
   quantities were both being called strand specificity: the simulator's is protocol **fidelity**
   (direction-agnostic), `rna_sense_frac` is the **directional** sense fraction. For an R1-antisense
   (dUTP) protocol — which is what the simulator emits, and what real cfRNA is — they are complements, so
   comparing them reads as a sign error. ⭐ `StrandModel.strand_specificity` is the matching quantity and
   already recovers the simulated knob: **1.00 → 1.0000, 0.75 → 0.7701, 0.50 → 0.5020**. So a fitted
   κ of 0.0101 on a "0.99 stranded" library is **correct**, which is also why forcing 0.99 measured 166×
   worse. ⛔ The cfRNA libraries are ordinary highly stranded **R1-antisense** ones and the scalar says so.
3. ✅ **ALL THREE ASSIGNMENT MODES NOW REPRODUCE EXACTLY with a fixed seed** (S2.1, 2026-08-01), so
   `sample` is a legitimate A/B arm. ⛔ It did not before, and the cause was **not** the seed: the
   scanner fills the fragment buffer in worker-completion order and the sampled assignment consumed its
   per-locus RNG in that order. `build_multi_loci` now orders each locus's units by `frag_id`.
   ⚠ Hold the mode fixed across both arms regardless: the three are different estimators
   (5441 / 6002 / 6277 on one scenario), not one answer at different precision.

⚠ **THE GOLDENS ARE STALE AND HAVE MOVED SIX TIMES: the suite is 1923 passed / 21 failed / 1 xfailed.**
All 21 are `test_golden_output` — P1 (the EM prior's units), C2 (the FL models), C2.6 (`L` itself), S1 (the
hold-out), P4 (the drain wired in), P4.2 (the combination rule). Every one expected.
⛔ **Regenerate ONCE, after C3, twice, and diff.** ⭐ "21 failed and all of them goldens" is the standing
baseline: a 22nd failure, or a non-golden name in the list, is a regression.

⭐ **`TODO.md` §7's `test_nrna_double_counting[g20_n0_s100]` now PASSES** (C2, 2026-08-01) — the silent
negative control reads **0 counts against a limit of 25** where it leaked ~30, stable across re-runs.
⚠ **Do not close §7 on that.** A negative control is one-sided (trap 19), this was not C2's target, and
§7's substance was always that `map` suppresses assignment most; re-read it against all three modes
before retiring it.

⛔ **The OLD benchmark baseline is void.** `r0 0.079005 / r3 0.046675` was the deleted `ambig_dense_10mb`
suite — do not quote it, compare against it, or try to reproduce it. The replacement suite is proven to
resolve 6 of its 8 requirements; ⛔ run `suite_resolves.py` before quoting any number from it.

**`tests/native/_accumulator_reference.py` is the executable specification** for the accumulator. The C++
is gated on byte-identity to it; where it and a document disagree, it wins.

⭐ **Every object stores THREE integer sums** (S5.a): `count` = Σ1, `inv_length_sum` = Σ round(2³²/placements),
`length_sum` = Σ L. ⚠ `inv_length_sum` is **not** called `density` on purpose — it is an exact model-free
density at an edge and is *not* one at a node. `length_sum` exists because the other two carry **zero**
information about the gDNA/RNA split when the two components share a mean length.

**The accumulator is being redesigned and replaced, not patched.** The accumulator is the tally built
during the BAM scan. Its central defect: it chops each fragment at partition borders and then decides
"this piece crosses its right border if another piece follows it" — which **cannot distinguish a
contiguous crossing from a splice jump**. Measured consequence: 18–22 % of "spliced" mass on real cfRNA
sits at positions with no annotated splice site.

The replacement, agreed with the owner:

> The genome is a graph. **Nodes** are intervals; **edges** connect them. A fragment is a **path**.
> Nodes count fragments *contained*; edges count fragments *crossing* (a 0-bp line, no width).
> Every object stores `(integer count, Σ 1/L)`.

- `L` is the fragment's own molecule length = genomic span **minus cut introns** (so a paired-end mate
  gap counts, an intron does not). Whatever counts toward `L` must also count as coverage for crossing,
  or the estimator is biased.
- **Two edge kinds.** *Contiguous* edges (the seam between genomically adjacent nodes) are bidirectional,
  carry gDNA + RNA, and their endpoints are **implicit** (edge `i` sits between node `i` and `i+1`).
  *Junction* edges are directed donor→acceptor, are **pure mature RNA** by construction, need explicit
  `(src, dst, strand)`, and drop the unspliced channel and the structural flags.
- A splice jump deposits on its **junction edge only** — never on the contiguous edges it splices over.
- `Σ 1/L` estimates density with **no fragment-length model**, but only at an edge and only where there
  is ≳1.5 fragment lengths of template on both sides. It is not model-free at a node, and near a
  transcript end it under-reads badly (0.11× at 20 bp). Both limits are *signal*, not defects.
- **ONE strand convention**: everything is stored by **genome** strand (`CHANNEL_PLUS`/`CHANNEL_MINUS`).
  *Sense*/*antisense* is the transcript-relative notion and is **derived, never stored**.
- **Two components only, gDNA and RNA.** "RNA is RNA" — no mature/nascent split in the accumulator.

**Phase order** (S1–S6 in the plan): index → **Python reference + spec matrix** → C++ → payload →
consumers → delete. The reference is written first and is the oracle for everything after it. Do not
graft new work onto the old accumulator.

⛔ **No version suffixes in file names.** It is `accumulator.py`, never `accumulator_v5.py`. Files are
rewritten in place and the old path is deleted, not kept for comparison.

## The index

`INDEX_FORMAT_VERSION 8`, shipped as `nodes.feather` + `edges.feather`, built and checked by
`calibration/splice_graph.py`.

✅ **2026-07-30: rebuilt from scratch at `~/Downloads/rigel_runs/refs/rigel_index`** after every index was
deleted. 69 s to build; the build is deterministic (two builds byte-identical in all seven artifacts).
⭐ **`manifest.json` now records the sources, their sha256 and the build flags**, so it never again has to
be inferred — the previous one held `format_version` and `rigel_version` and nothing else.

⚠ **The census below describes AN ANNOTATION, not the tool** — GENCODE v46 / Ensembl 112 with ERCC
controls. It re-derived exactly on 2026-07-30: **1,043,881 nodes** (median 151 bp, mean 2,970, 56.7 %
shorter than one 200 bp fragment), 1,043,595 contiguous edges, 404,168 junction edges. A rebuild from a
different GTF moves every one of those numbers, so **re-derive** — `scripts/design/index_census.py`.

- Nodes tile each reference, cut at **every exon endpoint** of every non-synthetic transcript, with no
  merging. 53.4 % of real transcript termini were invisible to the previous merged partition.
- Edges always run `src < dst`, so **genomic order is a topological order — there is no graph traversal
  anywhere**. The graph is a DAG; it is not a polytree (every junction edge closes one undirected loop),
  so junction edges must be *factors on their endpoint nodes*, never message channels.
- 8 structural flag bits per contiguous edge: TSS / TES / DONOR / ACCEPTOR × {+,−}, not mutually
  exclusive. Carry the raw bits to the consumer; do not pre-derive predicates in the plumbing.
- Validated by invariants I1–I13, two of which re-derive the answer by a **different algorithm**. That
  discipline caught two real bugs; a validator that calls the builder's own helper validates nothing.
- 1 bp nodes are legal (15,687 of them). Nothing may assume length > 1.

## Build, test, lint

> **Every build, test and lint command must run inside the activated `rigel` conda environment** — it
> holds htslib and the compilers, and the C++ build finds htslib via `$CONDA_PREFIX`.

⛔ **Run the suite as `python -m pytest`, never bare `pytest`.** Bare `pytest` does not put the repo root
on `sys.path`, so `tests/calibration/test_fl.py`'s `import tests.native._accumulator_reference` raises and
the suite reads one extra failure — which, under "a new failure is a regression", reads as one.

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel

pip install --no-build-isolation -e ".[dev]"   # rebuild after ANY src/rigel/native/ change
python -m pytest tests/ -q                     # 1923 pass / 21 fail — all 21 stale goldens (regen after C3)
python -m pytest tests/ --update-golden        # regenerate tests/golden/ after intended output changes
ruff check src/ tests/ scripts/ && ruff format src/ tests/   # ⚠ NEVER format scripts/
```

**Tooling that is current** (everything else under `scripts/` was deleted 2026-07-30 as unrunnable;
`scripts/sim/fl_estimation_stress.py` followed on 2026-08-01, having called a `calibrate()` signature that
had not existed for several milestones):

| | |
|---|---|
| `scripts/design/index_census.py` | re-derive an index's census — never quote the numbers, run this |
| `scripts/design/verify_index_rebuild.py` | nodes byte-identical, edges only in contiguous reach |
| `scripts/design/suite_resolves.py` | ⛔ **run before quoting any suite number** |
| `scripts/design/build_scan_cache.py` | scan once, calibrate many times. ⚠ **re-run after any accumulator change** — the payload schema digest invalidates every cache, by design |
| `scripts/design/composition_evidence_census.py` | ⭐ how much library mass reaches the solver with NO composition evidence. `--inject-kappa 0.5` is its falsification handle |
| `scripts/design/prior_units_check.py` | the EM prior in fragment units vs the old incidence sum, both arms from one calibration |
| `scripts/design/length_likelihood_ab.py` | the P2 A/B: `CalibrationConfig.length_likelihood` False vs True over the pilot |
| `scripts/design/native_parity_on_real_data.py` | the S3 gate on real cfRNA at full scale |
| `scripts/design/scan_profile.py` | ns/fragment, regressed over several BAMs |
| `scripts/design/observable_efficiency.py` | what fraction of the length information a storage choice keeps |
| `scripts/design/fl_anchor_gap.py` | ⭐ the zero-gDNA falsification: EB anchor vs RNA pool over the 8 pilot caches, scored on `truth_fragment_lengths.tsv` so no target is chosen. ⭐ **`--drain` measures every panel BEFORE and AFTER the second pass** |
| `scripts/design/calibration_truth_ab.py` | ⭐⭐ **THE DELIVERABLE, scored against truth** — the library gDNA fraction, undrained vs drained, against the simulator's own origin counts. ⛔ Run this before claiming any calibration improvement |
| `scripts/design/held_flux_census.py` | D-3: how often a held candidate has ZERO flux evidence, decomposed by cause and scored against `truth_abundances.tsv` |
| `scripts/design/node_density_derivation.py` | the reciprocal-opportunity theorem, T0–T6, each perturbed |
| `scripts/sim/build_suite_reference.py` · `design_suite_probes.py` · `simulate_reads.py` | build the suite |
| `scripts/sim/evaluate_suite.py` | net fragment flow (`rigel.sim.analysis`) |

Always set `OMP_NUM_THREADS=1` when benchmarking or comparing runs.

## Architecture

**Pipeline** (`pipeline.py`): three stages.

1. **Scan** (`scan_and_buffer`) — C++ htslib single-pass reader. Resolves fragments against the index,
   trains the strand and fragment-length models from unique mappers, buffers fragments into a columnar
   `FragmentBuffer`, and deposits per-object tallies into the C++ **accumulator** → `AccumulatorPayload`.
2. **Calibrate** (`calibration.calibrate`) — deconvolves each node's unspliced mass into gDNA vs RNA by
   a belief-propagation sweep over the node chain, and fits the library-level strand and density
   parameters. Output: `CalibrationResult`.
3. **Quantify** (`quant_from_buffer`) — scores fragments, builds loci by connected components, runs a
   per-locus EM with `n_transcripts + 1` components (every transcript row plus one gDNA component). The
   calibration prior enters as **two per-locus Dirichlet scalars**, never per-transcript.

**Python modules.** Top level: `cli` `pipeline` `config` `index` `scan` `scoring` `buffer`
`scan_payload` `_accumulator` `locus` `locus_partition` `scored_fragments` `estimator` `strand_model`
`frag_length_model` `splice` `splice_blacklist` `native` `gtf` `transcript` `annotate` `stats` `types`.

`calibration/` (33 modules): `calibrate` (orchestrator) · `splice_graph` (the v8 index) · `bp_solver`
(the solver) · `node_chain` `node_geometry` `node_init` (chain construction and per-node setup) ·
`substrate` `region_arrays` `signature` (payload views and geometry) · `effective_length`
`capture_eff_length` `fl` (length models) · `strand_likelihood` `gdna_strand` `strand_balance`
`strand_deconv` `strand_summary` (strand) · `density_deconv` `density_model` `gdna_landscape` `npmle`
`background_reference` (density and priors) · `simplex` `simplex_logodds` `enrichment_frame` `derive`
`run_fill` (solver internals) · `priors` `result` `errors` `diagnostics` `track`.

**C++** (`src/rigel/native/`, nanobind, C++17, `-O3`, LTO, OpenMP):

| Module | Source | Purpose |
|---|---|---|
| `_bam_impl` | `bam_scanner.cpp`, `calibration/accumulator.cpp` | BAM parsing, fragment grouping, model training, the accumulator |
| `_em_impl` | `em_solver.cpp` | per-locus EM, connected components (Kahan summation, `fast_exp.h` SIMD) |
| `_scoring_impl` | `scoring.cpp` | fragment likelihood scoring (`-ffast-math`, no SIMD) |
| `_resolve_impl` | `resolve.cpp` | fragment→transcript resolution via cgranges |
| `_cgranges_impl` | vendored | interval overlap |

The C++ accumulator must match `tests/native/_accumulator_reference.py` **byte for byte** — that Python
file is the specification, not the documentation.

**nRNA components** are not per-transcript shadows: unique nascent spans keyed by
`(ref, strand, start, end)` are shared across transcripts and materialized as ordinary transcript rows
in `index.t_df` (flagged `is_nrna` / `is_synthetic`). On a **non-synthetic** row `is_nrna` means
"single-exon, so mature ≡ nascent" — **not** "manufactured span". Using it as a realness filter deleted
52,104 real transcript termini. The real-transcript filter is `~is_synthetic`, alone.

## CLI

```bash
rigel index --fasta genome.fa --gtf annotation.gtf -o index/
rigel quant --bam sample.bam --index index/ -o results/
rigel sim --config scenario.yaml -o out/
rigel export results/ -f tsv
rigel report results/ -o report.html
```

Input BAM must be name-sorted with the `NH` tag.

## Working rules

- **No magic numbers.** Stop and discuss before adding any constant, heuristic or tunable. Every divisor
  must be derived from the deposit rule and unit-tested against brute-force enumeration.
- **No legacy, no backwards compatibility, no speculative code.** Converge and delete. Code kept "for
  comparison with the old version" is a defect.
- **No Greek letters in identifiers** (fine in maths write-ups).
- **One thing varied per experiment**, with a falsification test written *first* and verified failing,
  and a baseline re-recorded from the current tree in the same session.
- ⛔ **Real data is a TEST input, NEVER a DESIGN input.** The cfRNA on disk is one far end of the RNA-seq
  spectrum, not a sample of it. Sweep the plausible space, report the worst case, and bring the owner the
  domain call — do not adopt whatever the available data happens to say. Owner, 2026-07-30.
- **Profile and gate on real cfRNA, never a small synthetic suite** — the deleted 10 Mb one ranked
  hotspots backwards, was Poisson by construction, and had zero fragment-length variance. ⚠ The simulator
  is still Poisson (`sim/wgs_engine.py:473`), so the replacement inherits that unless it is built out.
- **The owner drives commits.** Do not commit unless asked.
