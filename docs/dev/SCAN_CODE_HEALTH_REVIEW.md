# Scan-path code-health review

2026-10-16 · base `41899ae4` · wave 1, step 1 (`PAIRED_END_WAVE1_PLAN.md` §5). Read-only: nothing in `src/`
changed. For the owner to select from.

## How it was done

- **Nine read-only slices, one checklist.** Each slice was reviewed by one agent, and every finding cites file and line.
  - The scan path in depth, six slices: the scanner's front half, its back half, the resolver, scoring and EM units,
    the calibration deposit with its second pass, and the Python scan plumbing.
  - The rest of `src/` swept, three slices: the calibration package; EM, index, CLI and config; and the simulator.
- **Known items.** Each slice started from `ISSUES: hygiene-ledger` and grepped ISSUES and TRAPS. Items already
  recorded are pointed to at the end, not repeated.
- **Verification.**
  - **Checked** means I re-read the cited code myself.
  - **Traced** means the slice traced the path and I did not repeat it. Two slices ran small numeric checks in pure
    Python; those are marked.
- **Evidence.** The nine full slice reports are in `scan_code_health_slices/`.
- **Baseline**, re-derived after a rebuild on a clean worktree at `41899ae4`: 0 failed / 2,959 passed / 0 skipped /
  0 xfail, in 88 s.

## The short version

About 90 distinct findings. Most are dead code, duplicated rules, stale comments and vocabulary: cleanup that should
move no number. The ones that matter, in order:

1. **The scan cache can return a stale tally.** Its key misses the splice blacklist and the BAM. Every wave-1 A/B reads
   caches, so this is fixed first.
2. **A BAM without NM is scored as if no read had a mismatch,** silently. STAR's default attributes omit NM, so the
   recommended STAR command must ask for it.
3. **Pairing without HI pairs mates by shared transcripts.** Every location inside one gene shares the gene's synthetic
   span, so all combinations pass. Step 3's bounded search must decide by geometry alone.
4. **The multimapper merge keeps each transcript's least-overhang hit,** not its most likely one. Overhang is counted
   twice, which corrects the wave-1 plan's §9.2 argument. This belongs to wave 2.
5. **Pruning uses three different references,** so gDNA can be cut while an equally likely span survives. Wave 2.
6. **Accounting and annotation are off in several places:**
   - `proper_pairs` counts hits and is shown as a share of records;
   - gated fragments are stamped "intergenic" in the annotated BAM;
   - a unique intergenic pair on two references is labelled chimeric;
   - about 20 scan counters are asserted by no test.

   All of this is step 6.
7. **The second pass decides "RNA" from the motif tag,** where the deposit uses the annotated intron, and it passes an
   AMBIGUOUS motif into its junction lookup. A candidate fourth correction in step 5.
8. **The resolver's strand table is installed only by the scan.** So every unit test runs without strands, and its
   ambiguous-strand branch is never exercised.

Proposed order:

1. **The cleanup stages.** Each is proven a numeric no-op.
2. **The defects inside the scan path.** Each one goes into its wave-1 step as its own arm.
3. **The rest.** ISSUES entries for later waves, which are the owner's call.

## 1. Fix inside wave 1

These touch the code wave 1 changes, or a measurement it relies on.

| name | what | where | slice · status | moves numbers |
|---|---|---|---|---|
| `scan-cache-key-misses-blacklist-and-bam` | The key hashes the partition, the sj CSR, the payload schema, the reach and the scan config. It does not hash the splice blacklist, which changes the tally and holds artifact fragments back, nor the BAM. `build_scan_cache` also skips an existing directory unless `--force`. | `scan_cache.py:35-45, 433-491`; `index.py:937-967` | Python scan · checked | no |
| `nm-absent-reads-as-zero` | A record without NM reads NM = 0, and nothing counts or warns. STAR's default `--outSAMattributes Standard` (NH HI AS nM) omits NM. That silently switches off the mismatch weighting between multimapper hits and the mismatch gate. | `bam_scanner.cpp:636-638` (writer `2756-2758`) | scanner front · checked | no |
| `no-hi-pairing-by-transcript-overlap` | Without HI, mates pair when their resolved transcript sets intersect. Every location in one gene shares the synthetic span, so all n1 × n2 combinations pass. The cross-pair step can only produce trans chimeras, so its comment is false. | `bam_scanner.cpp:938-951, 985-991` | scanner front · checked | yes, input without HI |
| `unindexed-references-share-id-minus-one` | Mates on references outside the index all get ref −1. `(-1) << 32` is undefined behaviour. Mates on two different unindexed contigs are merged as one reference. | `bam_scanner.cpp:658, 690-693` | scanner front · traced | labels only |
| `resolver-strands-installed-only-by-the-scan` | The index's resolver gets its transcript strands, K and maximum fragment length only in `scan_and_buffer`. Every unit test therefore runs with no strands (`ambig_strand` always 0), and the annotation writer depends on the scan's setters. | `pipeline.py:281-292`; `resolve_context.h:1322-1336` | resolver, Python scan · checked | no |
| `intergenic-append-is-dead-and-mislabels` | Unique fragments with no candidate are appended to the buffer for a consumer that no longer reads it. Its one effect is that an intergenic pair on two references is counted and annotated chimeric. | `bam_scanner.cpp:1881-1897`; `scoring.cpp:701-707, 738` | scanner back · checked | tags, a logged count |
| `proper-pairs-counts-hits` | `proper_pair` counts hit groups. summary.json files it under record-level counts and the report divides it by records, so a perfectly paired unique library shows 50%. | `bam_scanner.cpp:1797-1812`; `cli.py:549-555`; `report/model.py:179` | scanner back · checked | summary.json, report |
| `gated-fragments-stamped-intergenic` | A fragment the scorer gates gets no annotation row, and the writer reads "no row" as intergenic (ZF 0x25 or 0x21, ZW 1.0). "One row per frag_id" is assumed, not checked. | `bam_scanner.cpp:2415-2420, 2656-2670`; `scan.py:140-171` | scanner back, Python scan · traced | tags only |
| `annotation-pass-can-drift-a-group` | Pass 1 drops a filtered record before it compares read names; pass 2 flushes on it. On an uncollated BAM, every later read name could carry another fragment's tags, silently. | `bam_scanner.cpp:1461-1487` against `2714-2726` | scanner back · traced | no |
| `scan-counters-untested` | About 20 scan counters are asserted by no test. `proper_pairs` survived this way. | `bam_scanner.cpp:339-375, 2041-2080` | scanner back · traced | no |
| `no-indel-cigar-fixtures` | No fixture carries CIGAR D, I, =, X or H, so a coordinate regression would go unseen. | `bam_scanner.cpp:65-68` | scanner front · traced | no |
| `second-pass-restates-the-motif-rule` | The second pass passes an AMBIGUOUS motif into its sj lookup, where the deposit trusts no splice. It also decides "this placement is RNA" from the motif tag, where the deposit uses an annotated observed intron. The C++ comment saying AMBIGUOUS "never reaches here" is false. | `second_pass.py:278-280, 598-610`; `accumulator.cpp:427` | accumulator · checked in part | yes, rare held placements on tagged stranded libraries |
| `single-block-length-shortcut` | A single merged block gets its own length for every candidate. It skips the subtraction of the candidate's introns inside it, which `tx_frag_length` and DESIGN §3.4 apply. So one molecule gets two lengths, depending on whether its mates overlapped. | `resolve_context.h:960-965` | resolver · checked | yes, mostly candidates already deep in overhang |
| `multimapper-splice-type-drops-implicit` | A multimapper made only of IMPLICIT or ARTIFACT hits reports UNSPLICED in ZS. | `scoring.cpp:505-511` | scoring · traced | ZS only |

## 2. Wave 2: how a template's hits combine

| name | what | where | slice · status | moves numbers |
|---|---|---|---|---|
| `merge-keeps-least-overhang` | Per transcript, the merge keeps the least-overhang hit, then the most likely. But each hit's likelihood already contains `oh × ln 0.1`, so overhang counts twice. A hit with overhang 0 and NM 3 beats one with overhang 1 and NM 0, a choice 4.6 nats worse. At the mismatch gate the transcript is lost. The two rules agree for same-place alternatives with equal NM. | `scoring.cpp:569-582` | scoring · checked | yes |
| `pruning-uses-three-references` | Annotated candidates are pruned against the best annotated candidate, synthetic spans against the best span, and gDNA against the best of both. So gDNA can be cut while a span reading the same footprint, nearly as likely, survives. | `scoring.cpp:604-629, 843-868, 883-888` | scoring · checked | yes |
| `mismatch-gate-nan-in-multimapper-gdna` | At `--mismatch-alpha 0`, two gated hits make the gDNA log-sum NaN, so a later clean hit's gDNA reading is lost depending on hit order. | `scoring.cpp:57-73, 519-527, 646-653` | scoring · traced, numeric check | at the gate only |

## 3. For later waves: candidate ISSUES entries (the owner's call)

| name | what | where | slice · status | moves numbers |
|---|---|---|---|---|
| `coverage-weight-plateau-for-short-transcripts` | The plateau is min(f, L/2); enumeration gives min(f, L − f). They differ when L < 2f, and the gate restates the formula. | `scoring.cpp:93-96, 111`; `tests/test_pipeline_routing.py:1004` | scoring · traced, enumeration | yes, through the EM warm start |
| `two-own-composition-thresholds` | τ > 0 and τ > 1e-9 both mean "has its own composition". | `solve_kernel.cpp:420-421, 607`; `region_init.py:38, 58` | calibration · traced | slots in the window only |
| `chimera-span-over-exon-blocks-only` | The outermost span is measured over exon-touching blocks only, so an outer intronic or intergenic block is unbounded at this stage. The owner picks the rule. | `constants.h:405-410` | resolver · traced | yes, if changed |
| `capture-probes-on-dropped-duplicates` | The ladder's panel names 2 transcripts the index dropped as duplicates. Their probes are skipped at debug level and still counted as loaded. | `sim/capture/sampler.py:470-476, 514-517` | simulator · checked | ladder capture-ON, on re-simulation |
| `contiguous-reaches-unread` | Contiguous-boundary reaches are computed, stored, gated and hashed into the scan-cache key, but read by nothing. The taper A/B they wait for was never run. | `splice_graph.py:317`; `region_geometry.py:152-155` | calibration · traced | no (cache keys change) |
| `second-pass-geometry-in-python` | `placement_abundance` re-derives the deposit's segments, crossings and region in Python. They are equivalent today, but can drift. | `second_pass.py:374-461` | accumulator · traced | no |
| `draw-clamp-can-pick-a-zero` | The clamp keeps a draw inside its run, but it can land on a zero-score hypothesis, about 10⁻³ events per library. | `second_pass.py:767-771` | accumulator · traced, numeric check | no, in practice |
| `output-nascent-two-meanings` | "Nascent" names the synthetic spans in summary.json and the report, but annotated single-exon transcripts in `gene_quant`. The report row's description is false. Fix the row text now; renaming the outputs is a user contract (the MANUAL, the goldens) and needs an owner ruling. | `report/model.py:298-304`; `estimator.py:229-234, 564-575`; `cli.py:433-434` | EM sweep · checked | no |
| `loci-table-recomputes-mrna` | The loci table rebuilds mRNA by subtraction, behind a clamp that hides drift. | `estimator.py:861-863` | EM sweep · traced | ulps in `loci.feather` |
| `simulated-errors-redraw-the-base` | A "substitution" is drawn over all four bases, so a quarter change nothing and the realised rate is ¾ of `error_rate`. No test simulates errors. | `sim/wgs_engine.py:461-465` | simulator · checked | none today (every panel uses 0.0) |

## 4. Out of scope under the narrowed scope

These are reported only, per the owner's ruling of 2026-10-16.

- **`anchor-tolerance-absorbs-overrun`.** `--splicing-anchor-tolerance` (default 3) cuts an annotated intron whole
  even when up to K aligned bases lie inside it. By its own comment, it exists because "an aligner routinely runs a
  read a base or two past a donor". Location: `resolve_context.h:693-696, 846-857`. Found by the resolver and EM
  slices.
- **`deposit-merges-conflicting-introns`.** When mates disagree, the deposit cuts the union of their introns, a path
  neither mate aligned. Location: `accumulator.cpp:456-458`, repeated at `second_pass.py:378-384`. The stage
  disagreement this causes is already in `ISSUES: splicing-artifacts`.
- **The alignable splice blacklist.** Already recorded; wave 4 retires it once the comparison justifies it.

## 5. Cleanup stages

Each stage moves no number. It is proven with `rename_identity.py --check` against a reference frozen once before the
first stage, the full suite at 2,959, and `preflight.py --full` after any deletion or rename.

### 5.1 The words: comments, docstrings and vocabulary

No code changes.

- **The K tolerance is described as a label rule.** In fact it shapes every fragment's gap hypotheses, and so
  calibration and the scan-cache key.
  - Locations: `resolve_context.h:677-679`, `resolve.cpp:119-123`, `config.py:212-222`, `pipeline.py:285-288`,
    `cli.py:1452-1455`.
  - In the same stage, the setter throws on a negative value, and the comment states why the default is 3.
- **`--em-mode map`.** The help, and `docs/MANUAL.md:193`, describe a `max(0, n + a − 1)` update. The code adds the
  pseudocounts and normalises (`cli.py:1319-1321`; `em_solver.cpp:617-665, 774-791`).
- **Population words.** Descriptions that use "mRNA"/"nRNA" or "mature"/"nascent" as populations are reworded as
  "the winner's kind", or as RNA that has or has not spliced there:
  - the report's "Nascent RNA" row;
  - the ZF bit documentation (`annotate.py:108-110`);
  - comments at `accumulator.h:210-211`, `resolve_context.h:821-823`, `constants.h:232`, `estimator.py:460`, `cli.py:342`,
    `html.py:146`, `report.js:277-278` and `transcript.py:86-89`;
  - the calibration comments (`simplex_logodds.py:48`; `splice_graph.py:36, 49, 56, 849-876`; `fl.py:22, 481`;
    `sj_opportunity.py:153`; `region_geometry.py:78`).
- **`is_spliced`.** It means "gDNA is no candidate", not "spliced"; its docstring and that of `log_liks` are corrected
  (`scored_fragments.py:48-51`).
- **Stale comments**, file by file, as listed in the slices:
  - the scanner: header, the identity, the 48 B `Region`, the scratch;
  - the resolver: three comments;
  - the side buffer: three places;
  - `splice.py`;
  - calibration: seven places;
  - EM and index: nine places;
  - the simulator: about fifteen places.
- **Three comments describe consumers that do not exist:** the supporting-transcript comment, "AMBIGUOUS never reaches
  here", and the `n_sj_blacklisted` consumer.

### 5.2 Dead code

**Native:**
- `WorkerState::deposit_scratch` (`bam_scanner.cpp:424`) and `ParsedAlignment::sj_strand` (`192-200`);
- the unreachable junction-merge branch (`resolve_context.h:1227-1232`) and the `n != ti` guard (`1138-1150`);
- `SweepCapture`'s five unread fields, and the strand-only ψ solve that feeds `fg_strand` (`blocks.py:57-59`,
  `sweep.py:357-365`, `solve_kernel.cpp:597-603`);
- ψ's signal mask, which solves spliced-only slots and discards the result (`simplex_logodds.py:342-345`; prove the
  no-op);
- digamma's `x <= 0` branch (`em_solver.cpp:141-147`);
- `ScoredHypothesis::absorbed`, and `stat_chim_`, which is returned and discarded.

**Python:**
- the zero-regions guard (`pipeline.py:305-307`);
- in the index and estimator: `TranscriptIndex.g_to_strand_arr`, `AbundanceEstimator(is_nrna=…)` and
  `IntervalType.SJ`;
- `FragmentLengthModel`'s unreachable smoothing and `CONSTANTS.fragment_length.unseen_smoothing_ess`;
- in calibration: `pooled_log_rate`, `_RNA_SIGNATURE_BITS`, and `is_terminus` / `is_splice_site`;
- `build_fl_models`' optional divisors. Production always passes all four, and tests use gDNA laws production never
  builds. Test expectations may move, so this is its own stage.
- `MultiLocus.gdna_span`'s floor.

**Defaults:** native defaults that no caller relies on become required arguments, and an unknown sj spec throws
(`bam_scanner.cpp:1305-1308, 3312-3313, 3330-3333, 3418-3419`).

**The simulator:**
- `strand_kappa`;
- the PyYAML guards;
- the unreachable branches, and the second seed;
- the argparse CLI, which moves into its script.

### 5.3 One definition each

- **The qc dict keys.** They are written three times; a single helper replaces them (`bam_scanner.cpp:2298-2308,
  3102-3112, 3289-3299`).
- **Codes and versions shared with native code:**
  - the EM's mode codes are exported from native, and an unknown code is refused (`estimator.py:351, 372`;
    `em_solver.cpp:55-57, 880, 901, 906`);
  - the summary schema version is defined once (`cli.py:208`; `report/substrate.py:25`);
  - the `FRAG_*` codes are exported, with a parity test (`buffer.py:39-43`; `constants.h:59-63`).
- **"auto" is resolved once,** inside `scan_and_buffer`, and the annotation writer takes the resolved config. Today
  six copies resolve it (`pipeline.py:974-976`; `annotate.py:355-400`; five instruments).
- **The exonic strand diagnostic becomes four counters.** Today it is per-fragment labels (`bam_scanner.cpp:249-250,
  455-465, 1969-1982`).
- **Calibration helpers are reused:**
  - `overlapping_region_runs`, which `priors.py:95-96` re-implements;
  - `region_lengths_from_partition`, which `gdna_opportunity_from_index` re-implements;
  - one intergenic predicate, where there are three today.
- **The simulator:** the read-name encoding is owned by `read_name.py`, and one writer produces both truth tables.
- **One `hash_combine` helper** (`constants.h:131`).

### 5.4 Contracts: refusals and checks that change nothing on valid input

- **`Accumulator::merge_from`** checks the histogram width before it merges anything, with a falsification test
  (`accumulator.cpp:969-995`).
- **`MarshalledFragment`** reads its hypotheses once, and its copy constructor is deleted (`bam_scanner.cpp:1039-1066`).
- **`ResolverScratch`'s** move operations are defaulted, with RAII for its buffer (`resolve_context.h:368-408`).
- **Configuration is validated:**
  - `pruning_min_posterior` must be in [0, 1) (`scoring.py:129-130`);
  - an `IndexConfig` is actually built from the CLI, and the tolerance gets one name (`config.py:439-451`).
- **Every index file is checked up front** and raises `IndexLoadError`, instead of failing with a traceback
  (`index.py:1323-1326`).
- **The scan payload and cache:**
  - the region ledger identity is asserted, and the histogram dtypes too (`scan_payload.py:782-804`);
  - the strand-model cache fields are derived from the dataclasses (`scan_cache.py:328-330`).
- **The simulator's config:**
  - `shares` and `transcript_filter` are refused, and the unused aliases dropped;
  - `oracle_bam: false` with `emit_fastq: false` is refused;
  - the manifest records the index.

### 5.5 One re-cache

This stage changes the payload schema. It is measured, because the draw order may move.

- **The supporting-transcript lists** (`hypothesis_t`, `gap_supporting`) are exported, validated and hashed, but only
  tests read them. Three slices found this. They are part of the canonical record key, so removing them can reorder
  tied records and shift the seeded draw.
- **The same re-cache takes:**
  - `hygiene-ledger` (b)'s `region_span_count` and `region_end_count`;
  - the three unread QC counters (`scan_payload.py:696-701`);
  - the contiguous reaches, if deleted.

### 5.6 Renames

Grep every sense of each token first, and run `rename_identity.py --check` after each.

- **Synthetic, not "nRNA":**
  - `t_is_nrna` holds the synthetic flag, so it becomes `t_is_synthetic`;
  - `mrna_best` / `nrna_best` and `MergedMrna` become annotated / synthetic names (`scoring.py:132-136`;
    `scoring.cpp:181, 398, 533, 605-611, 793-800`).
- **Names that say the wrong thing:**
  - `chose_genomic`, which also counts multimapper choices;
  - `MultiLocus.gdna_span` → `span_bp`;
  - `SweepCapture.mature` (or delete it, per §5.2);
  - `exonic_truth`, which holds the annotated strand, not truth.
- **One word, two senses:**
  - the simulator's "shadow" against the index's;
  - `before_second_pass` → `before_annotated_bam`;
  - `FragmentRouter.scan`, which scores rather than scans.
- **Output column names** (`mature_count` / `nascent_count`, `mrna_*` / `nrna_*`) are a user contract, and need an
  owner ruling first.

### 5.7 Bare constants

These move into `CONSTANTS` or `constants.h`:
- the tag detection's `max_spliced_reads = 1000`, which is written twice;
- the QC tile thresholds and rounding, and the junction-depth cuts (`report/model.py`; `strand_model.py:160-161,
  488-489`; `index.py:152, 322-323`; `cli.py:353-355, 583`);
- the landscape's tile literals (8 and 16), with a gate that the tiling cannot move a number;
- the buffer's spill capacity (`buffer.py:291, 463`) and the blacklist's `min_count` defaults;
- the simulator's I/O literals;
- the resolver's −999 sentinel.

### 5.8 Test gaps outside §1

- the merge's replace branch: two hits on one transcript, at the default alpha and at the gate;
- pool separation in pruning;
- the coverage weight against enumeration;
- the mismatch gate on multimappers;
- an AMBIGUOUS motif in the sj-lookup gate;
- the landscape tiling;
- `write_truth_abundances` and `_introduce_errors_batch`.

## 6. What this changes in the wave-1 plan

These are proposed, pending the owner's selection.

1. **§9.2 is corrected.** It claimed that the least-overhang hit is the right likelihood for a transcript. That rule
   counts overhang twice (`merge-keeps-least-overhang`). Wave 1 keeps today's rule, because for same-place alternatives
   with equal NM both rules choose the same hit. Wave 2's combination rule decides.
2. **§6, locations:**
   - count the records without NM, and warn;
   - keep unindexed references' own ids;
   - check that NH and HI are integer-typed.
3. **§7.1, pairing:** the no-HI search decides by genomic geometry, never by transcript sets. A new fixture puts both
   copies inside one gene.
4. **§9.3:** a candidate fourth second-pass correction. A placement's annotated-sj use and strand come from one rule,
   the deposit's.
5. **§10, accounting:**
   - `proper_pairs` in the right unit;
   - gated fragments annotated as gated;
   - the intergenic append deleted;
   - the annotation pass's read-group count checked against the scan's;
   - every scan counter asserted exactly;
   - the multimapper ZS type.
6. **§12, measurement:** the scan-cache key covers the blacklist and the BAM before any wave-1 arm reads a cache.
7. **The index installs the resolver's strand table at load,** so tests exercise the ambiguous-strand branch.
8. **The resolver's single-block length shortcut** is removed in its own arm.
9. **The input contract (wave 4's manual):** the recommended STAR command includes NM, which STAR's default omits.

## 7. Already recorded

These are pointers by name only, not repeated here:
- `hygiene-ledger`, items (b) through (i);
- `latent-defects`;
- `splicing-artifacts`;
- `sj-strand-tag-chosen-from-the-first-reads`;
- `multimapper-intergenic-alignments`;
- `performance-memory-bounded-solve`;
- `binary-cuts-on-continuous-quantities`;
- `input-mistakes-pass-silently`;
- `the-em-prior-share-excludes-drawn-multimappers`;
- `scoring-penalties-are-underived-constants`;
- `overlapping-synthetic-shadows`;
- `an-unrecorded-splice-blacklist-is-dropped-silently`;
- `TRAPS: nrna-does-not-mean-synthetic`, `TRAPS: a-hash-that-misses-its-artifact` and
  `TRAPS: a-transcript-predicate-must-not-silently-drop-a-molecule`.

**Corrections to recorded entries:**
- `multimapper-intergenic-alignments` cites `bam_scanner.cpp:1808-1845`; the append is now at `1892-1897`.
- `scoring-penalties-are-underived-constants` says "a zero alpha gives NaN today" and points at a removed item.
  `mismatch-gate-nan-in-multimapper-gdna` is the live remainder.
- `splicing-artifacts` says the census counts uniquely mapped fragments; `latent-defects` says it counts placements.
- `hygiene-ledger` (i) should name `transfer_rows.h`'s uncommented `EPS` and `TINY`. `EPS` serves as a variance floor,
  a flatness threshold and a fraction clamp.

**Housekeeping:** 15 old prototype worktrees remain registered, for the owner to prune.
