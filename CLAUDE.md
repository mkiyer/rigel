# CLAUDE.md

Guidance for Claude Code working in this repository.

## What Rigel is

A Bayesian RNA-seq transcript quantifier that jointly models mRNA, nascent RNA and genomic DNA
contamination. A single-pass C++ BAM scanner tallies fragments, a **calibration** stage deconvolves the
library into gDNA vs RNA, and a per-locus EM solver assigns RNA to transcripts. PyPI package
`rigel-rnaseq`; the import and CLI are `rigel`.

## ⭐ Current direction — read this first

**Four documents are the whole story. Read them in this order:**

| doc | what it is |
|---|---|
| **`docs/IMPLEMENTATION_PLAN.md` §0** | ⭐ **START HERE** — live state, next actions, what is uncommitted |
| `docs/ACCUMULATOR_DESIGN.md` | the design being implemented |
| `docs/LEDGER.md` | what has landed, its gates, and why each thing is the way it is |
| `docs/CARRY_FORWARD.md` | 24 measured facts, 18 equations the code depends on, **27 traps**, 30 design ideas |

On 2026-07-29 the project's other 274 docs (74,823 lines) and 132 agent-memory files were deleted and
distilled into these. **Nothing else in `docs/` is a design document**, and a doc path not listed above
does not exist — several older references to one were dangling.

**`tests/native/_accumulator_reference.py` is the executable specification** for the accumulator. The C++
is gated on byte-identity to it; where it and a document disagree, it wins.

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

## The index (built, validated, current)

`INDEX_FORMAT_VERSION 8`, shipped as `nodes.feather` + `edges.feather`, built and checked by
`calibration/splice_graph.py`. Human scale: **1,043,881 nodes** (median 151 bp), 1,043,595 contiguous
edges, 404,168 junction edges.

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
pytest tests/ -q                               # 95 modules, 1298 tests
pytest tests/ --update-golden                  # regenerate tests/golden/ after intended output changes
ruff check src/ tests/ scripts/ && ruff format src/ tests/
```

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
- **Profile and gate on real cfRNA, never the 10 Mb synthetic suite** — it ranks hotspots backwards,
  is Poisson by construction, and has zero fragment-length variance.
- **The owner drives commits.** Do not commit unless asked.
