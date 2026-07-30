# CLAUDE.md

Guidance for Claude Code working in this repository.

## What Rigel is

A Bayesian RNA-seq transcript quantifier that jointly models mRNA, nascent RNA and genomic DNA
contamination. A single-pass C++ BAM scanner tallies fragments, a **calibration** stage deconvolves the
library into gDNA vs RNA, and a per-locus EM solver assigns RNA to transcripts. PyPI package
`rigel-rnaseq`; the import and CLI are `rigel`.

## ⭐ Current direction — read this first

**Read them in this order:**

| doc | what it is |
|---|---|
| **`docs/S5_DESIGN_LOG.md`** | ⭐ **START HERE for S5** — §0 status, §1 the accumulator changes, §2 the live plan. **It supersedes `IMPLEMENTATION_PLAN.md` §4/§5** |
| `docs/IMPLEMENTATION_PLAN.md` §0 | live state for everything else |
| `docs/NODE_DENSITY_DERIVATION.md` | why the deposit weight is what it is, and what each stored channel buys |
| `docs/TODO.md` | the one deferred-work list, ranked, each item with the reason it is deferred |
| `docs/ACCUMULATOR_DESIGN.md` | the design being implemented |
| `docs/LEDGER.md` | what has landed, its gates, and why. Older entries: `LEDGER_ARCHIVE.md` |
| `docs/CARRY_FORWARD.md` | ⭐ **§3 traps then §2 equations** — the most-used reference in the project |
| `docs/BENCHMARK_SUITE.md` | the suite: how to build it, and **what it can and cannot judge** |

Reference rather than design: `BENCHMARKING.md` (how to evaluate — net fragment flow), `MANUAL.md`,
`PUBLISHING.md`, `docs/testing/testing_plan.md` (the owner's plan for the cached-substrate harness).

⛔ **`calibrate()` DOES NOT RUN**, and it still gates the benchmark suite producing any number, the scan
cache's toy seed, and every future A/B. **S5.0/a/b/c/d have landed; S5.e is next**, then S5.f.

⚠ **The 291 failures are not one thing.** ~266 are the original consumer breakage; ~25 are tests of the
per-face geometry model that S5.c/S5.d deleted, and **those fail with a message naming S5.e or S5.f**. A
failure that says which step owns it is the normal state here.

⛔ **There is no benchmark baseline.** `r0 0.079005 / r3 0.046675` was the deleted `ambig_dense_10mb`
suite — do not quote it, compare against it, or try to reproduce it. The replacement suite exists and is
proven to resolve 6 of its 8 requirements, but cannot produce a calibration number until S5.

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

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel

pip install --no-build-isolation -e ".[dev]"   # rebuild after ANY src/rigel/native/ change
pytest tests/ -q                               # 1384 pass / 291 fail / 15 error — the failures ARE S5
pytest tests/ --update-golden                  # regenerate tests/golden/ after intended output changes
ruff check src/ tests/ scripts/ && ruff format src/ tests/   # ⚠ NEVER format scripts/
```

**Tooling that is current** (everything else under `scripts/` was deleted 2026-07-30 as unrunnable):

| | |
|---|---|
| `scripts/design/index_census.py` | re-derive an index's census — never quote the numbers, run this |
| `scripts/design/verify_index_rebuild.py` | nodes byte-identical, edges only in contiguous reach |
| `scripts/design/suite_resolves.py` | ⛔ **run before quoting any suite number** |
| `scripts/design/build_scan_cache.py` | scan once, calibrate many times |
| `scripts/design/native_parity_on_real_data.py` | the S3 gate on real cfRNA at full scale |
| `scripts/design/scan_profile.py` | ns/fragment, regressed over several BAMs |
| `scripts/design/observable_efficiency.py` | what fraction of the length information a storage choice keeps |
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
