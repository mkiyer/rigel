# Learning capture without a panel: implementation and validation plan

First written: 2026-10-04. Revised: 2026-10-05 after the updated `CAPTURE_CLEAN_SLATE.md`.

Status: executable development plan, not a description of shipped behavior. This is a sandbox
document. Nothing outside `docs/dev/` may cite it. Implement prototypes and measured A/B arms before
changing production defaults. Move accepted findings into their permanent documentation homes when
they land; do not make this file authoritative by citation.

**OWNER RULING, 2026-10-05 (evening): capture is a spectrum; a binary capture detector is absolutely wrong and must never be built.** Every library-level test of capture in this plan — the decision contract (§2), the conditional test, its envelope, tail and e-value, the gDNA admission gate and the power study (§4.3–§4.6), the statuses, Phase Detector and Phase Gate — is ABANDONED and deleted below. The evidence bank (§4.1–§4.2), the physical operator, the source laws, the field, the yields, the likelihood and the composition phases remain, now as the main line. Plasma panels enrich extremely scarce transcripts that stay a tiny fraction of the library, so a library-level statistic reads the uncaptured background; the shipped reader's mode-or-`None` is such a decision and its cost was measured in the Phase 0 sweep (the two sparse plasma libraries quantified as uncaptured).

## 0. Execution readiness (2026-10-05)

Verified on the day, in the campaign directory `~/Downloads/rigel_runs/prototypes/2026-10-05_capture_learning/`
(`README.md` there maps every file; `manifest.json` records commit `bdfd8709`, the dirty-diff and source-tree hashes,
the environment and every panel's manifest hash and cache counts; `protocol.json` holds the section 2 contract with
every unchosen value `null`):

| item | state |
|---|---|
| toolchain, references, panels | `preflight.py` passes: ladder 16/16 and test chromosome 30/30 scan and oracle caches, certified; the two suite gap arms (4 + 4) and the test chromosome's `fl_rna_long`, `fl_gdna_long`, `fl_equal200` and `probes_junction` panels (30 each) cached |
| the seams this plan names | every file, function and test file in sections 6, 8 and 10 exists at `bdfd8709`, and every instrument flag in section 10; the claims checked against the tree: index format 7 at the tag against 8 (`INDEX_FORMAT_VERSION`), `deposit_digest` hashes top-level arrays only, `SJStrandTable` counts once per fragment at its leftmost annotated junction, the tag's CLI has `--threads` (mapped to `scan.total_threads`), `--assignment-mode fractional` and `--seed` |
| `release_tag` | built from the worktree `~/proj/rigel-v071` (`pip wheel --no-build-isolation --no-deps`, unpacked into `release_tag/site`); `release_tag/rigel071.py` launches its CLI from that site and refuses the editable tree; its own index of the test chromosome built (format 7; `--collapse-duplicate-transcripts --no-mappability`, matching the tree's index manifests); one pinned quant: 5.9 s, 0.69 GB peak |
| the adapter | `phase0/score_release.py` scores any `rigel quant` output directory with `quant_accuracy.py`'s own `score_transcripts` / `score_genes` and rebuilds the library row from `summary.json`'s `quantification` block (`mrna_total`, `nrna_total`, `gdna_total`, `intergenic_total`, present at the tag and at HEAD); on today's CLI output it reproduces `quant_accuracy --arm base` on the same condition to the printed digit (transcript abs error 50,972.042, gene 6,563.284, gDNA fraction 0.503) |
| the first tagged number | test chromosome `gdna_g50_ss_0.99_nrna_file_capture_off`: 0.7.1 transcript abs error 57,845 / gene 8,070 / gDNA fraction 0.506 against today's 50,972 / 6,563 / 0.503 |
| Phase Audit task 6 | reproduced: 834/1709/3166 to 799/1744/3166; 206 fractions above 1 (max 5.3); 699 finite of 834 (`phase0/census_audit.py`) |
| Phase Audit task 7 | `phase0/pool_contrast.py` on 17 conditions (`results/pool_contrast_*.txt`): the pooled boundary count over its constant-intensity expectation is 0.98–1.00 on every capture-OFF row of every panel (the RNA-short defect row 0.99, the RNA-long arm 0.99, the test chromosome's fl arms 1.00 / 0.977 / 0.983), 10.2–19.8 on every captured row and 11.9 on `g00 ON` |
| not yet run | the tag on the ladder, the gap arms and the real libraries (tasks 2 and 9); task 8 (the crossing census on the VCaP transcriptome half); Phase Evidence's Python preflight |

### Readiness verdict after inspecting the driver

**The scientific Phase 0 work has a go. The full sweep is not safe to launch with the current driver.**
The completed pilot and adapter check remain useful. No numerical alpha, field model or native evidence
bank is needed for the remaining baselines. Four execution fixes below precede the long run.

`phase0/run_release_tag.sh` currently uses `INDEX_LABEL/CONDITION` as the output namespace, trusts
`summary.json` as completion, then appends scores unconditionally. The previous commands supplied
`suite` as the label for both gap arms and the ladder. Their four g50 condition names coincide:
24 requested suite quantifications map to only 16 output directories. The later panels would skip
quantification and score the first panel's output against different truth. Those commands are removed.
This collision does not affect the existing test-chromosome pilot.

Implement these changes in the outside-tree campaign scripts before running them:

1. **Separate reference identity from run identity.** Give the driver distinct `index_label` and
   `panel_id` arguments. Reuse `idx_suite` across the three suite panels, but write quantifications to
   `baseline/release_tag/<panel_id>/<condition>/`. Record panel ID, resolved BAM/condition paths and
   build identity in every score. A package version string is insufficient: both trees report 0.7.1.
2. **Make reuse verifiable and scoring idempotent.** Before reusing an index, validate a receipt of
   FASTA/GTF hashes, index options and tagged build. Before skipping quantification, validate a
   successful-exit receipt, BAM/index/build/config fingerprint and required output files. A summary
   alone is not completion. Record truth and scorer identity separately for scoring. Write each
   condition's three score rows atomically; rebuild aggregate JSONL deterministically with unique
   `(arm, panel_id, condition, axis)` keys. Rerunning must neither duplicate rows nor change their
   input association. Preserve the pilot and adopt its artifacts only after validating provenance.
3. **Complete the paired current arm.** Run today's instrument on all the same simulated conditions,
   with one output namespace per panel, explicit seed 0, one scan thread and fractional assignment.
   Freeze/record current code and native identity for the campaign. The single current test-chromosome
   command did not cover the ladder or either suite gap arm.
4. **Add a real-library path.** The existing tag runner hardcodes `sim_oracle.bam`; the adapter requires
   `truth_abundances.tsv` and `truth_summary.json`. Neither is a real-data runner. Add manifest-driven
   quantification for an explicit BAM path and separate real-data comparison, with no fabricated
   simulated truth. Build the tag's format-7 real index from the sources listed below and preserve
   the original real-library outputs. Run both tagged and current versions on the same BAMs.

Use this explicit simulated-run inventory, expanding conditions from each panel's manifest and
checking that each has the required BAM/truth files:

| `panel_id` | Input directory under `~/Downloads/rigel_runs/` | `index_label` | Expected conditions |
|---|---|---|---|
| `test_base` | `test_reference/scenarios` | `test` | 30 |
| `flgap_rna_short` | `suite/flgap_rna_short` | `suite` | 4 |
| `flgap_rna_long` | `suite/flgap_rna_long` | `suite` | 4 |
| `ladder` | `suite/ladder` | `suite` | 16 |

This first sweep is 54 simulated conditions per arm, 162 unique scored axis rows per arm. The other
cached test-chromosome variants are separate additions to the run inventory, not covered by this count.
Order the defect condition first within its panel, then the remaining gap rows and ladder. Keep one
genome-scale job active; record logs, exit codes and completion receipts and use resumable execution.

The required revised driver interface is a specification, **not yet implemented** in the inspected script:

```text
run_release_tag.sh PANEL_DIR REF_FASTA REF_GTF INDEX_LABEL PANEL_ID CONDITION...
```

Before the full simulated run, exercise the revised scripts with a lightweight stub runner/scorer:
two panels sharing the same condition names must execute distinct BAMs and produce distinct score
keys while reusing one matching index. A second invocation must add no duplicate rows; changed BAM,
reference or build identity must refuse reuse; a partial output containing a summary must not pass as
complete. Test scorer failure and interrupted output writing. Then rerun/reuse the verified pilot and
one gap condition, inspect their receipts, and release the serial queue. This is a driver correctness
gate, not another estimator experiment.

#### Repairs implemented and the gate's results (2026-10-05, afternoon)

The driver is now `phase0/driver.py` (Python, receipted), with `run_release_tag.sh` kept as the specified
interface `PANEL_DIR REF_FASTA REF_GTF INDEX_LABEL PANEL_ID CONDITION...` over it. What it does:

- **Identities apart.** `index_label` keys the index (`release_tag/idx_<label>/receipt.json`: FASTA and GTF
  hashes and sizes, options, build identity, format version, build command); `panel_id` keys the run
  (`baseline/<arm>/<panel_id>/<condition>/receipt.json`: BAM hash and size, the index receipt's hash, the build
  identity, the pinned config, the command, exit code, the hash of every required output, `ok`). The build
  identity is the wheel's hash for the tag and `commit:src_tree_sha:native_so_sha` for the current tree (the
  five editable `.so` files hashed), never a version string. A format-7 manifest records no sources, so an
  index without a receipt is never adopted: the pilot's index was moved aside and rebuilt (0.7 s).
- **Verifiable reuse.** Quantification is skipped only against a receipt with `ok`, an equal fingerprint and
  every output present at its recorded hash; a complete run with a different fingerprint is refused (exit 4)
  and left untouched; a partial output without a complete receipt is discarded and rerun; a failed run keeps
  a non-`ok` receipt and gets no score. Scores are written atomically beside the run, keyed
  `(arm, panel_id, condition, axis)`, with the run receipt's hash and the truth and scorer hashes, and are
  recomputed only when one of those changes; `collect` rebuilds `baseline/<arm>.jsonl` deterministically and
  refuses a duplicate key.
- **The paired current arm.** `--arm current_cli` runs today's `rigel quant` with the same pinning
  (`--threads 1` maps to `scan.total_threads`, `--assignment-mode fractional`, `--seed 0`) against the tree's own
  index, validated against the request's FASTA and GTF hashes through its manifest; it is scored by the same
  adapter, which reproduces `quant_accuracy --arm base` to the printed digit on the pilot. The instrument
  commands above remain valid as an independent check and are not the campaign's `current` rows.
- **The real path.** `driver.py real --bam --sample-id` quantifies an explicit BAM under `baseline/<arm>/real/<sample>/`
  with the same receipt and writes `qc.json` (the summary's quantification block, the calibration scalars, the
  output table's hash and count sum); it never touches the simulated-truth scorer. The tag's human index is
  built from `refs/genome_controls.fasta.bgz` and `refs/genes_controls.sorted.gtf` under `index_label` `real`.
- **The gate.** `phase0/test_driver.py` runs the stub runner and scorer through the real driver: 24 checks pass
  (distinct directories and keys for two panels sharing condition names; one index serving both; a rerun with
  no quantification and no byte changed; a changed BAM refused with the receipt untouched; a changed reference
  and a changed build identity refused before anything runs; a partial output with a summary rerun; a failed
  run's non-`ok` receipt, no score, then its rerun; a scorer failure writing no score and no temp file; a
  deterministic collect; the real path's receipt and `qc.json` with no score; a dry run creating nothing).
- **The pilot.** Preserved under `baseline/pilot_2026-10-05/`; rerun through the driver for both arms, every
  score identical to the pilot's, receipts in place (the tag 7.1 s, the current tree 7.2 s).
- **The queue.** `phase0/run_queue.py run sweep_jobs.json` runs the 16 jobs of `phase0/make_jobs.py` serially
  (both arms on the four simulated panels, 54 conditions per arm, the defect row first in its panel; then both
  arms on the four real libraries), with `phase0/mem_guard.py` beside it (kill any one process above 20 GB, the
  largest when together above 64 GB) and a state file that records every exit code; a failed job is recorded
  and the queue continues. The gap-condition check (`gap_check_jobs.json`: the defect row, both arms, the tag's
  suite index built on the way, 2.0 s) ran through the same queue before the sweep was released.
- **The gap check's result (the launch gate).** Both runs completed with `ok` receipts on the same BAM hash
  (`0bf25cc9…`), 261 s and 265 s, 8.9 GB peak; the first scoring attempt failed because the queue runner was
  named `queue.py` beside the driver and shadowed the standard library's `queue` for pyarrow — renamed
  `run_queue.py`, a 25th check (no file in `phase0/` shadows a stdlib module) added, and the two runs were
  re-scored from their receipts without re-quantifying. The defect row, RNA-short stranded × OFF: the tag reads
  **no enriched mode** (`capture.enriched: false`, peak-to-peak fold 1.0, on-target mass 0 over 20,557 nodes)
  and transcript error 6.31 % / gene 0.44 % / gDNA fraction 0.513 against the truth's 0.500; the current tree
  reproduces its defect, 42.30 % / 3.28 % with the 10.07/bp reference from 393 members. The design's prediction
  that the mass-weighted KDE would not read the false mode held. The sweep (`sweep_jobs.json`, 16 jobs) was
  released at 08:45 through `run_queue.py`; `run_queue.py status sweep_jobs.json` reads its progress.
- **The simulated half's table (2026-10-05, 11:37; 54 conditions per arm, every row receipted; `results/summary_simulated_all.txt`).** Transcript and gene error as Σ|Δ| over the stratum's conditions as a percentage of its true RNA fragments:

| panel | stratum | 0.7.1 tx % / genes % | current tx % / genes % | ahead |
|---|---|---|---|---|
| ladder | stranded × OFF | 3.95 / 0.44 | 2.02 / 0.25 | current |
| ladder | stranded × ON | 9.91 / 2.52 | 6.14 / 1.12 | current |
| ladder | unstranded × OFF | 10.29 / 1.49 | 2.29 / 0.31 | current |
| ladder | unstranded × ON (deferred) | 33.21 / 14.25 | 12.23 / 3.69 | current |
| ladder | g00 false gDNA (str OFF / str ON / unstr OFF / unstr ON) | 41,473 / 25,843 / 401,708 / 21,285 | 127 / 183 / 126 / 335 | current |
| RNA-short gap arm | stranded × OFF (the defect) | 6.31 / 0.44 | **42.30 / 3.28** | **0.7.1** |
| RNA-short gap arm | stranded × ON | 20.33 / 3.52 | 19.82 / 1.52 | current |
| RNA-short gap arm | unstranded × OFF | 6.00 / 0.61 | 3.78 / 0.27 | current |
| RNA-short gap arm | unstranded × ON (deferred) | 29.76 / 6.19 | 18.84 / 2.15 | current |
| RNA-long gap arm | stranded × OFF | 5.23 / 0.46 | 2.26 / 0.19 | current |
| RNA-long gap arm | stranded × ON | 7.38 / 1.62 | 6.17 / 0.62 | current |
| RNA-long gap arm | unstranded × OFF | 5.29 / 0.59 | 2.13 / 0.22 | current |
| RNA-long gap arm | unstranded × ON (deferred) | 28.99 / 5.84 | 8.24 / 0.98 | current |
| test chromosome | stranded × OFF | 8.20 / 1.20 | 7.39 / 1.08 | current |
| test chromosome | stranded × ON | 8.58 / 2.45 | **13.89 / 4.01** | **0.7.1** |
| test chromosome | ss 0.70 × OFF | 9.44 / 1.54 | 8.61 / 1.25 | current |
| test chromosome | ss 0.70 × ON | 10.28 / 3.32 | **15.09 / 5.75** | **0.7.1** |
| test chromosome | unstranded × OFF | 9.25 / 1.76 | 8.78 / 1.33 | current |
| test chromosome | unstranded × ON (deferred) | 18.94 / 10.31 | 18.15 / 7.80 | current |

  Read it this way. The current tree is behind the release in exactly two places: the RNA-short defect row, where the
  no-reference figure (3.53 %) would beat 0.7.1's 6.31 %; and the
  test chromosome's captured stranded rows (stranded × ON 13.89 against 8.58, ss 0.70 × ON 15.09 against 10.28), a loss
  the ladder does not show (6.14 against 9.91 the other way). Both trees detect capture on those rows, so that loss sits downstream of
  detection, in what each tree does with the reference, and it is the one place where the bar as stated (every required
  stratum of every panel) is not met by Phases 0–3 alone. The release's own weaknesses are elsewhere: it books 128k to
  586k unspliced-RNA fragments as gDNA on the capture-OFF rows of the genome-scale panels and reads 21k to 402k false gDNA
  fragments on the ladder's zero controls, which its transcript table does not see; the current tree's zero controls are
  two to three orders of magnitude cleaner. The receipted current-tree ladder rows reproduce the 2026-10-03 pinned record
  to the digit, so the CLI path and the instrument agree at genome scale as they did on the pilot.
- **The sweep completed (2026-10-05, 12:35): 16 of 16 jobs, none failed, 3.4 h of queue time, no guard kills.**
  The real half (`results/summary_real.txt`; no truth for plasma; `qc.json` per run): the tag's human index built in
  8 s (752,654 regions, format 7). Both arms agree on every library's total fragment count to within ~30 and on
  κ; they disagree on how much gDNA they read and on whether capture is detected:

  | library | 0.7.1: gDNA fraction, its KDE reading | current: gDNA fraction, its located mode |
  |---|---|---|
  | LBX0190 (plasma, 147k fragments) | 0.130; enriched, fold 14.7, on-target 0.62 | 0.083; no located mode |
  | LBX0588 (plasma, 1.29 M, gDNA-rich) | 0.966; enriched, fold 2,246, on-target 0.76 | 0.915; 0.297/bp from 12,154 members |
  | MO_3021 (plasma, 835k) | 0.225; enriched, fold 1,079, on-target 0.03 | 0.155; no located mode |
  | VCaP mix (18.5 M; exome DNA + transcriptome, per-fragment truth by read name) | 0.263; enriched, fold 254, on-target 0.72 | 0.236; 0.086/bp from 24,512 members |

  The tag declares capture on every real library, as the panels warrant; the current reader locates a reference
  on the two deep libraries only and quantifies the two sparse plasma libraries as uncaptured, reading less gDNA
  than the tag on all four. That is the binary reader's cost the owner's ruling names, and the regime the kernel is
  for; the VCaP origin labels give the one real truth for the gDNA fraction: by read-name prefix over the
  BAM's primary first-in-pair alignments (`results/vcap_origin_counts_r1.txt`; A00839 exome DNA 4,708,907,
  HWI-D00127 transcriptome 13,994,657) the gDNA fraction is 0.2518, which the tag reads +0.011 high (0.263) and
  the current tree −0.016 low (0.236). That domain is the aligner's, not the tools' fragment filter (18.70 M
  against 18.5 M fragments), so it is an approximate library-level truth; the matched-domain origin-pool
  comparison the real-input readiness section asks for is still to do.

Current-tree simulation commands, after recording the frozen build (these flags already exist):

```bash
CAPTURE_CAMPAIGN=~/Downloads/rigel_runs/prototypes/2026-10-05_capture_learning
CAPTURE_TEST=~/Downloads/rigel_runs/test_reference
CAPTURE_SUITE=~/Downloads/rigel_runs/suite
export OMP_NUM_THREADS=1
python scripts/design/quant_accuracy.py --arm base --em-seed 0 \
  --set scan.total_threads=1 --set em.assignment_mode=fractional --jobs 1 \
  --suite "$CAPTURE_TEST/scenarios" --index "$CAPTURE_TEST/idx" \
  --out "$CAPTURE_CAMPAIGN/baseline/current_test_base.jsonl"
for CAPTURE_PANEL in flgap_rna_short flgap_rna_long ladder; do
  python scripts/design/quant_accuracy.py --arm base --em-seed 0 \
    --set scan.total_threads=1 --set em.assignment_mode=fractional --jobs 1 \
    --suite "$CAPTURE_SUITE/$CAPTURE_PANEL" --index "$CAPTURE_SUITE/rigel_index" \
    --out "$CAPTURE_CAMPAIGN/baseline/current_${CAPTURE_PANEL}.jsonl" || exit "$?"
done
```

Attach `panel_id` and input/build receipts when ingesting these instrument rows into the comparison;
condition basenames are not unique across files. Do not append all panels into an unkeyed stream.

### Real-input readiness

Four local BAMs were found under `~/Downloads/rigel_runs/cfrna/<sample>/bam/star.srt.rmdup.collate.bam`:

| `sample` | Role |
|---|---|
| `mctp_LBX0190_SI_43883_HHHKGDRX7` | Plasma |
| `mctp_LBX0588_SI_43875_HHHKGDRX7` | Plasma |
| `mctp_MO_3021_SI_42465_HGL7VDRX7` | Plasma |
| `mctp_vcap_rna20m_dna05m` | VCaP mixture |

Their quick checks passed, and all 286 reference names/lengths match the FASTA index. The existing
format-8 real index is `~/Downloads/rigel_runs/refs/rigel_index`. Its source files are
`refs/genome_controls.fasta.bgz` and `refs/genes_controls.sorted.gtf` under that run root; their
recomputed hashes match its manifest. Build the tagged real index from those same sources and options.
Add these inputs, hashes, explicit sample IDs and output locations to the campaign manifest; they are
absent from the inspected manifest. No fourth plasma input has been inventoried.

For plasma, compare output/QC, count conservation and runtime; there is no transcript truth score.
For VCaP, origin labels permit an origin-pool comparison on a matched fragment domain. The RNA-only
and DNA-only BAMs exist at `prototypes/2026-09-27_strand_od/data/vcap_rna_only.bam` and
`vcap_dna_only.bam`; verify their provenance when using them. RNA-half transcript/gene estimates are
reference comparisons, not known transcript abundances. The VCaP crossing census remains a separate
Phase Audit task. The real branch can proceed once its manifest, index and quantification/QC path
pass the same receipt and resume checks; it must not be routed through the simulated-truth adapter.

### Decisions after Phase 0

The four 2026-10-05 choices are superseded: the decision rule and the decision-only deliverable are abandoned by the
owner's ruling; the bank (Python fixtures first, then the native bank) and the shared amplitude stand. What
replaces the veto as the 0.8.0 candidate is the owner's call, informed by the test chromosome's root cause in the
design document's §10.

## 1. Outcome, scope, and execution order

Build a capture estimator that uses ordinary sequencing observations, requires no probe panel, and
does not mistake a poorly determined gDNA density for enrichment. Separate three questions:

1. Is there admissible evidence of gDNA and of spatial enrichment?
2. Which capture-dependent placement probabilities and component yields are actually identified?
3. Does using those probabilities improve calibration and the delivered transcript counts?

Do not rebuild the clean-slate design's physical `beta` plus exon `probed_fraction` parameterization.
Those parameters are not jointly identified by its observations, and a fraction is not geometry.
A replacement for the ruler lands only after its evidence, source length laws, numerators, denominators and count
frame agree; there is no capture decision anywhere (owner ruling).

The deliverable (the decision-only deliverable is abandoned by the owner's ruling above):

- **Capture likelihood:** a learned placement kernel, identified yield functionals or honest bounds,
  and consistent EM numerators and normalizers.

The implementation must respect information limits. It must not return a precise RNA yield when
different fields fitting the available gDNA evidence give different RNA yields. Missing information is
an output of the estimator, not permission to impute a fully probed gene.

### Recommendations on the revised design's four questions

These are the execution recommendations for `CAPTURE_CLEAN_SLATE.md` section 10. They do not turn
unmeasured results or unspecified numerical policies into release approvals.

| Question | Recommendation | Evidence needed before promotion |
|---|---|---|
| When to build the bank | First prove the schema, ownership, observability and tiny-case statistics through the Python specification on synthetic offered events. Then build the minimal native banks in an isolated worktree before the real-BAM Gate experiment. | The Python accumulator is not a BAM scanner. A claimed Python-only end-to-end prototype needs a validated scanner-event replay; do not build a second BAM classifier merely to postpone the native bank. |
| One binding amplitude | Keep one shared amplitude as the first library model, conditional on its fit and transfer checks. VCaP can falsify it but raw per-context dispersion cannot decide the issue alone. | Simulated variable-binding/acceptance fixtures and held-out predictive/yield failures, then a locked VCaP diagnostic controlling geometry, source laws, depth and local DNA baseline. Add a population prior only after a separate-amplitude diagnostic shows a useful effect. |

The executable path is: **Audit (done) → Python evidence preflight → native Evidence → Geometry → Laws →
Field → Yields → Likelihood → Composition → Release**, one mechanism per phase. The richer spliced-law
observations and all-contributor deposit/path records are later dependencies; do not claim the first two banks
suffice for the complete likelihood.

### Start here

An executing agent should perform these steps in order:

1. Read `CLAUDE.md`, `docs/SUCCESS.md`, `docs/DESIGN.md` section 0 and its scope rulings,
   `docs/EQUATIONS.md` sections on count frames and conserved lengths, and `docs/TESTING.md`.
2. Read the 2026-10-05 revision of `CAPTURE_CLEAN_SLATE.md` and `FRAGMENT_LENGTH_POSTMORTEM.md` in this
   directory for context. Use this plan's explicit contracts where the revision compresses an
   assumption, omits required metadata, or makes an unmeasured exactness claim.
3. Record the current commit, dirty diff, environment, panel manifests, and source hashes. Preserve
   the existing dirty working tree. Create isolated prototype/worktree artifacts; do not reset it.
4. Execute the phases in section 8. Each phase names its inputs, outputs, tests, and stopping rule.
5. Maintain the artifact manifest and failure ledger in section 11. A failed prerequisite prevents
   dependent promotion, not independent measurement work.

The three release target strata remain stranded/OFF, unstranded/OFF, and stranded/ON. Report
unstranded/ON separately on every relevant experiment; do not use it to choose the next optimization.
Do not silently convert a captured unsupported library to a successful capture-OFF result. The Q3
captured-library collapse is a required fallback test even though its stratum is deferred.

### Constraints retained from the request and repository

- No production probe BED, simulator labels, truth origins, or known panel parameters.
- Calibration populations stay `gDNA`, `RNA+`, `RNA-`. Continuous RNA support is geometry, not a new
  composition population. Nascent simulation conditions are robustness tests, not a new solver class.
- A library without usable gDNA evidence cannot learn a map. Distinguish this operational convention
  from the mathematical fact that assumed-uniform continuous RNA can also show a boundary contrast.
- No member census, enriched-density mode, or deconvolved pseudo-count used as capture evidence.
- No unexplained density cutoff, minimum opportunity, fold-change cutoff, or minimum witness count.
- No unconditional exon/non-exon landscape split and no reintroduction of refused mechanisms.
- Real libraries are test inputs. Do not select a model or threshold to fit the real samples.
- One mechanism per A/B. Native prototypes belong in an isolated worktree. Do not commit unless asked.

## 2. (Deleted: the decision contract, by the owner's ruling of 2026-10-05)

## 3. Mathematical and data-frame contract

### 3.1 Placement, observation, and capture

A placement includes reference/template identity, start, molecular length, orientation, and the whole
molecular path. A pair of aligned read blocks is not automatically a uniquely determined molecule.
An unsequenced insert may admit different splice paths.

For a component `c` and placement `z`, define:

```text
source_pmf[c, length]        source placement law, before capture/opportunity selection
acceptance[c, z]             probability or indicator for the modeled observation process
overlap[c, z, part]          overlap with one contiguous probe part on c's molecular path
weight[c, z]                1 + max_part(binding[part] * overlap[c, z, part])
kernel[c, z]                source_pmf[c, length(z)] * acceptance[c, z] * weight[c, z]
yield[c]                    sum over the admitted placement domain of kernel[c, z]
```

The simulator's equal-binding model is the special case with one common binding amplitude. In a
variable-binding extension, use the maximum of the part responses, not a sum. Union coverage is a
descriptive label only. Adjacent genomic parts do not establish that a probe joins them across a
transcript junction.

Use the full interval-intersection operator even at a boundary. In coordinates relative to a left
exon edge, a crossing placement is `[x - length, x)` and a part is `[p, q)`, so overlap is
`max(0, min(x, q) - max(x - length, p))`. The revision's shorter formula subtracting only `p` fails
when the part extends past the fragment start: `length=100, x=50, p=-200, q=20` gives 70 bases,
not 220. Mirror and test the right edge, overhangs and intron controls with the same operator.

The physical operator belongs to production-owned code. Production must not import `rigel.sim`, a
prototype directory, a BED reader, or an oracle. Keep the simulator as a separately exercised oracle;
retain an independent molecule-centric enumerator so sharing a bug cannot certify the operator.

`rho_off` is a library-normalized source-placement intensity. It is neither absolute molecular input
abundance nor the current `gdna_density_global` QC statistic. Existing output fragment counts remain
post-capture counts. A source abundance, if later reported, is proportional to `count / yield` and
must be separately named. Do not silently change existing TPM semantics.

### 3.2 Incidence and conserved mass

Keep these measures distinct in types and function arguments:

```text
incidence_opportunity[o] = sum_z source_pmf[length(z)] * 1(z counted on object o)
incidence_capture[o]     = sum_z source_pmf[length(z)] * 1(z counted on o) * (weight[z] - 1)
conserved_opportunity[o] = sum_z source_pmf[length(z)] * deposit_share[o, z]
conserved_capture[o]     = sum_z source_pmf[length(z)] * deposit_share[o, z] * (weight[z] - 1)
```

Include observation acceptance in each sum when it is part of the domain. A boundary incidence can
count the same molecule at several boundaries. Conserved shares sum to one across the full deposit,
and need not sum to one within a particular locus. Never replace the last expression by a product of
an average share and an average capture weight.

For a locus with routing/deposit function `D_locus(z)`, its denominator is the integral containing
`D_locus(z)`. A whole-fragment overlap union is interchangeable only after proving that
`D_locus(z)` equals that union's indicator on every admitted placement. Shared objects, short external
neighbors, and straddlers require explicit checks.

### 3.3 Numerator/denominator identity

The EM estimates captured component counts. For observation `i`, its candidate kernel is the sum over
compatible latent placements/paths of `kernel[c, z] * P(i | c, z, accepted)`. This conditional
observation law must sum/integrate to one over admitted observations for each accepted placement;
`acceptance[c, z]` is already in the kernel and must not be multiplied in again. Its component
likelihood is `candidate_kernel[c, i] / yield[c]` on that same admitted observation domain.

Dropping capture from a numerator is valid only when the capture factor is common to every compatible
candidate placement for that observation. A denominator-only intervention is an experimental arm,
not a complete likelihood. Capture must enter before pruning, competition, and alternative-hit merging
when the final scorer is implemented.

### 3.4 Source laws

Give source and census laws different names and provenance. The current `gdna_pmf` is intended to be
uniform-frame; `gdna_realized_pmf` is census; the current RNA law does not establish capture-corrected
source sampling. Do not reinterpret these fields silently.

For each origin, observed lengths have probability proportional to the source probability multiplied
by the total capture-weighted opportunity at that length. Using that observed histogram as a source
law inside another capture integral selects for capture twice.

Retain a structural-support mask before numerical PMF floors. No admissible placement means exact
zero yield. Numeric floors cannot create physical support, and no component-specific rescaling or
clamping to plain span is allowed for raw physical yields.

## 4. The evidence bank: length-resolved, disjoint molecular events

The kernel needs new raw evidence; the old marginal boundary counts cannot reconstruct it. (This section once
also held a library-level capture test; it is deleted by the owner's ruling.)

### 4.1 Event ledger

Add a typed `CapturePlacements` bank to `AccumulatorPayload`. Initial fields:

```text
ref_id: int32
start: int64
end: int64
align_strand: uint8
read_layout_id: integer
acceptance_id: integer
count: uint64
```

Initially retain only unique, nonchimeric, pass-one determinate, unambiguously unspliced accepted
placements. The whole molecular length is `end - start`. Fill the bank using the existing native
deposit acceptance logic, not a separate Python interpretation of a BAM. Retain excluded counts by
reason: explicit splice, deferred/implicit path, multimap, chimeric, invalid length, undefined strand,
and unsupported acceptance geometry.

Canonical-sort and combine identical rows with integer multiplicity before export. Preserve a stable
independent molecule-split label if later cross-fitting is enabled: put the split identifier in the
aggregation key, not an origin-bearing read name. Do not deduplicate distinct molecules merely because
they share coordinates. PCR duplicate policy remains the scan's existing policy and is recorded.

Store large banks in bounded chunks using the established spill pattern; never allocate a dense
genome-by-length array. Maintain exact raw counters that reconcile offered, retained, and excluded
events. The capture bank must not change when a deferred gap is resolved using model-dependent drain
choices. Adding deferred evidence later requires marginalizing its hypotheses, not treating a chosen
path as raw truth.

Define `capture_offered` at a single documented scanner hook (`Accumulator::deposit` is that hook today: every
fragment passes it and `deposited + deferred + dropped_* == offered` already holds there, so the bank's partition
is taken at the same point) and classify every offered fragment into
exactly one retained/excluded outcome, using the scanner's existing rejection precedence. Then assert
`retained + excluded_all == capture_offered`. Separately assert
`retained + excluded_accepted_unspliced == pass1_accepted_unspliced` on the matching domain, and
reconcile to the existing deposited/deferred/dropped counters. Explicit splices, deferred paths and
invalid fragments cannot all be added to the accepted-unspliced side of one ledger, as in the revised
design's abbreviated identity.

Keep layout/acceptance IDs unless a single-class restriction has been proved and recorded. A key of
only `(ref, start, end, strand)` cannot support layout-specific exposures. Derive storage from the
finalized dtypes and report array bytes, sorting/merge peak, spill volume and cache size; do not rely
on the revision's approximate 25-byte row for a memory budget.

For the first field prototype retain all eligible unspliced placements, including contained exon
positions and off-target anchors, with streaming spill. A boundary-only bank is not sufficient evidence for the
field model.
If an annotation-only retention mask is needed, derive its support before counts are inspected,
apply it identically to exposure enumeration and record what additional scan the field will require.

The ordinary fragment buffer currently lacks explicit reference identity for this purpose; intergenic
records cannot recover reference from transcript candidates. This bank must carry `ref_id` directly.

Also emit `SplicedAdmissionCounts(stratum_id, n_sense, n_opposite)` with raw integer counts by the
predeclared orientation/read-layout/error strata used in section 4.5. Count each certified-spliced
fragment once, regardless of how many junctions it crosses, and retain any independent split label
in its aggregation key. Reconcile exclusions and prove disjointness from the unspliced bank.
The current junction-indexed `SJStrandTable` already counts each qualified fragment once at its
leftmost annotated junction (`strand_model.py` and `native/bam_scanner.cpp`). It lacks the required
layout metadata; that missing stratification, not duplicate junction counting, justifies the new bank.
If eligibility and error rate are genuinely common, its raw integer marginal is a valid initial
single-stratum input. Fitted
`StrandModel.p_r1_sense` and fractional `genuine_n_same` are not hypergeometric observations.
Do not defer this small bank until the richer spliced evidence in Phase Laws.

### 4.2 Annotation and ownership

Implement a pure geometry function that maps any hypothetical admitted placement to exactly one of:

- `boundary`: crosses at least one eligible boundary;
- `interior`: wholly contained in an eligible intron control and not assigned to boundary;
- `excluded(reason)`.

A boundary is eligible only if an annotated mature path cannot account for the supposedly continuous
observation. Reuse and test the scanner's path semantics; a coarse exon/intron type alone is not this
proof. Contiguous exon/exon boundaries are excluded from the clean contrast. Explicit splice junctions
are a different observation axis.

Use a fixed genomic ownership rule when several eligible boundaries are crossed: choose the smallest
`(ref_id, boundary_coordinate, boundary_id)`. This tie-breaker merely defines a partition. Its exact
choice cannot alter total retained molecular mass; exposure enumeration must use the identical rule.
Do not select a boundary by count, inferred enrichment, transcript expression, or fitted composition.

Partition the candidate start axis into maximal connected structural windows with the same complete
set of continuous RNA templates/spans able to contain the whole placement. Include strand and finite
template reach in that signature. The two strand-admissibility bits alone are insufficient for
overlapping genes. Group boundary/interior events within a window, exact molecular length, observed
alignment strand, read-layout class, and acceptance class. Do not pool disconnected genes with
different RNA amounts merely because their coarse signatures match. Strand-specific acceptance and
source rates must not be averaged before conditioning.

For each cell retain:

```text
WitnessCell:
  reference_id, length_bp, window_id, support_id, align_strand, read_layout_id, acceptance_id
  boundary_count, interior_count                   scalar integer counts
  boundary_exposure, interior_exposure
  boundary_start_moment, interior_start_moment
  min_start, max_start
  geometry_digest, eligibility_digest
```

Enumerate zero-count opportunities too. A cell with only one of boundary/interior opportunity has no
contrast information; record it as such. It does not receive an arbitrary pseudocount. Cells with one
observed molecule are legal. There is no minimum count or minimum opportunity threshold.

Implement a slow exact placement enumerator first, then an interval sweep over structural endpoints
and length-specific translated endpoints. Require exact agreement. Do not reuse the existing generic
`region_geometry.eff_rna` boundary divisor: it intentionally lacks the finite reach required here.

## 5. Source laws and learned capture information

### 5.1 Estimate source laws without double selection

Yield recovery needs source-law estimation. Implement an
explicit `SourceLengthModel` with PMF, structural support, estimation method, uncertainty, and data
provenance. Reuse valid parts of `calibration/fl.py`; do not create a competing undocumented law.

First construct an oracle-law control, strictly in the evaluation harness. Before a production RNA
law fit, add a `SplicedLawObservations` bank or a deterministic replay pass retaining each selected
molecule's splice blocks/junctions, candidate templates, starts, lengths, strand, layout, and
alternative-path observation likelihoods. Canonical aggregation is permitted only for identical
observation/candidate-support classes. The global `pool_lengths[RNA_SPLICED]` histogram and the
unspliced-only `CapturePlacements` bank cannot supply these observations.

With a frozen candidate capture model, use the joint raw-observation likelihood. For observation
class `j`, define its Poisson mean as

```text
mu[j] = sum_t abundance[t] * sum_(z compatible with j)
          source_rna[length(z)] * acceptance[t, z] * weight[t, z]
          * P(j | t, z, accepted)
```

Include all selected-observation classes, including zero-count opportunities, and the matching
normalizer. A conditional multinomial is valid when conditioning on the declared total. Use one
shared RNA source law initially; a transcript-specific extension requires additional identification
evidence. The pooled length histogram alone identifies only
`source_rna[length] * sum_t abundance[t] * opportunity[t, length]` up to normalization. Do not
silently plug baseline transcript abundances in as known or assume block updates remove this ambiguity.

Implementation sequence:

1. Enumerate observation-specific capture-weighted opportunities and candidate support at every
   supported length, including truncation and junction-selection probability.
2. Fit the joint likelihood above, or begin on a demonstrably identifiable homogeneous subset.
   Unsupported mixtures return unresolved source laws; retain source/abundance uncertainty.
3. For a single known homogeneous template, the unfloored source estimate reduces to observed count
   divided by its per-length opportunity, followed by normalization. Use this as an analytic test,
   not a pooled correction for a mixture of unknown transcript abundances.
4. In supported mixtures, fit source PMFs and abundance nuisance jointly or in likelihood-increasing
   block updates, keeping the capture model fixed for this arm. Retain the current EB prior only with
   its stated estimand and likelihood; an observed census is not an independent source prior.
5. Profile the source-law/capture tradeoff. If length selection and capture cannot be separated on
   the available support, return bounds/unresolved rather than a point law supplied by a PMF floor.

Start each fit from the uncaptured law, a capture-selected law, and the previous-stage fitted law;
require agreement of identified predictive quantities or retain multiple solutions. Record objective
values and convergence diagnostics. Numerical convergence tolerances must follow a stated numerical
error target, not a capture/no-capture cutoff.

The two-length fixture is mandatory: source probabilities `(1/2, 1/2)` and weighted opportunities
`(1, 3)` give observed probabilities `(1/4, 3/4)`. Treating that census as the source yields `2.5`
instead of the correct source-weighted opportunity `2.0`.

### 5.2 Establish what can be learned before fitting a geometric map

Use the retained position/length evidence, not one density per exon. The first transfer experiment is
an identified contiguous-placement functional. For common contiguous support `A`:

```text
gDNA intensity(z) = rho_off * source_gdna[length(z)] * acceptance_gdna[z] * weight[z]
RNA yield on A   = sum_(z in A) source_rna[length(z)] * acceptance_rna[z] * weight[z]
```

With truth gDNA labels in the diagnostic only, an unbiased estimator at known laws/intensity is:

```text
yield_hat = sum_DNA_in_A (source_rna[length] * acceptance_rna[z])
                       / (source_gdna[length] * acceptance_gdna[z]) / rho_off
```

Acceptance cancels only when it is identical
for both origins on this domain. Require positive `source_gdna * acceptance_gdna` wherever the RNA
target has positive `source_rna * acceptance_rna`; missing observation support remains unresolved.
Its Poisson variance has squared importance weights; use the multinomial covariance when conditioning
on library size. Estimated intensity and laws add uncertainty and must be fitted/profiled, not treated
as known.

This experiment transfers the actual best-part response without recovering probes. It must recover
the functional where RNA support lies within observed gDNA support, and become uncertain/unidentified
as that overlap disappears. It does not recover a spliced placement's capture.

Distinguish zero support from inadequate finite information. The current gap configurations give
both source laws bounds of 50–500 bases: different means do not establish disjoint support. Report
target probability outside gDNA support separately from importance-weight variance and effective
sample size where support is positive. Add a truly disjoint-support fixture. The direct estimator
is a diagnostic/control; it need not become the production kernel representation.

Replace truth labels with the original raw strand likelihood on eligible stranded supports. A gDNA
contribution supplies both orientations; RNA has its admitted orientation/error model. Fit RNA amounts
as nuisance using both raw strands. Opposite-strand contained fragments are not certified gDNA:
RNA strand leakage must remain in their likelihood. Never convert a posterior fractional gDNA count
into a new Poisson observation. In
unstranded RNA-admitting interiors, the functional may remain unidentified; record that result.

### 5.3 Bounded geometric inference prototype

Do not start by building an unrestricted inverse solver. Implement a finite correctness model:
zero or one unknown contiguous part in each small isolated genomic context, with one nonnegative
effective binding amplitude shared by the library. No hierarchical amplitude prior is needed for
this first model. Begin with a single context to certify the operator and then the joint fit.
Enumerate all integer endpoint pairs in the context plus the null. An interval
ending outside observable support needs a censored-end alternative, not an artificial clipped edge.
The input is the context and raw observations, never a target BED.

For each candidate interval:

1. Evaluate its whole-placement weight using the physical operator.
2. At each value of the shared amplitude, fit admitted RNA/background nuisance using the raw count/strand
   likelihood, or a conditional local placement likelihood where all origins share the same spatial
   factor. State which nuisance was conditioned out and which remains.
3. Include zero-count admissible placements in the normalizer. Never fit only observed start sites.
4. Retain alternative fits and label the result's uncertainty level: structural identified set,
   profile/sensitivity envelope, or calibrated confidence set as defined below. A point MLE is not
   evidence that the endpoints or amplitude are identified.
5. Minimize/maximize every requested yield over the full retained parameter set, including all
   admissible amplitudes, source PMFs and nuisance values for each geometry. Retaining only each
   geometry's MLE or several optimizer starts is insufficient. Evaluate the same field jointly for
   all components sharing the context. Censored endpoints and arbitrarily large admissible amplitudes
   must remain in the set; an infinite bound is an unresolved yield.

For genuinely disjoint contexts and conditionally fixed shared nuisances, the scalar outer profile is
`profile_loglik(binding) = sum_context max_(interval, local_nuisance) loglik_context`. Shared source
laws, global intensity constraints or overlapping observation domains can prevent this factorization;
profile them jointly and merge interacting contexts before claiming independence. Retain all relevant
interval/nuisance combinations at each amplitude, not only the individually best interval. Sparse
contexts borrow through the shared amplitude assumption, not through an all-exons-targeted rule.
A flat global profile or amplitude/width tradeoff is unresolved even if an optimizer converges.

For the conditional spatial likelihood, each length/support cell is a multinomial over its admissible
starts with probabilities proportional to the fitted spatial intensity times capture weight. Use it
only when the common-factor/source-support assumptions hold. This can identify within-cell shape
without identifying the relative scale of different cells; gDNA intensity evidence and clean
off-target support are required to put those scales into one library frame. A flat conditional shape
does not establish unit capture. To constrain those scales, use a declared joint Poisson likelihood
over admitted unspliced field-observation classes, combining DNA and RNA means of the form in section
5.1 and explicit
off-target anchor observations, or its correctly normalized fixed-total multinomial counterpart.
If only the conditional spatial likelihood is used, leave unconstrained laws and between-cell scales
free in the retained set. Intervals conditional on frozen estimated laws/intensities are diagnostic,
not unconditional confidence bounds.

Use three levels of uncertainty work, in order:

1. Exact noiseless observational-equivalence fixtures: retain every fitting geometry/amplitude and
   verify its yield range. These are structural identification statements, not sampling intervals.
2. Finite-data likelihood profiles and yield sensitivity envelopes. Expose the likelihood cutoff as
   an analysis axis; a cutoff or a few optimized fits is not a confidence guarantee. Report incomplete
   searches as incomplete and test repeated-draw coverage before any accuracy claim.
3. A calibrated retained set only when needed for promoted yield-confidence claims. The construction
   below is one available method, not a prerequisite for the first field profile.
   Exact geometry enumeration alone does not calibrate continuous-amplitude or nuisance uncertainty.

For an executable finite-sample confidence construction, randomly split independent molecules before
fitting. Fit `q`, a normalized joint predictive law for the entire held-out observation, on training
or independent data. Every fitted input to `q`, including source PMFs, strand error, selected contexts,
masks and model choices, must use those independent data. On the held-out split retain every parameter
tuple `theta` satisfying `likelihood_theta >= alpha_confidence * q`, using the confidence allocation
from `protocol.json`. Under a true model in the declared family, Markov's inequality gives coverage
at least `1 - alpha_confidence`, because `E_true[q / likelihood_true] <= 1` (equality if the true
support covers `q`). Include nuisance parameters in the set; profile by existence of an admissible
nuisance value. Conditioning on held-out split/cell totals is allowed only if both laws are normalized
conditional on exactly those totals. Missing multinomial normalizers invalidate the construction.
A balanced split maximizes the smaller fitting/validation sample size; record it before outcomes.

Bounds on multiple yields from one joint retained set inherit its simultaneous coverage. Separate
context sets need the declared familywise allocation. Certify extrema over continuous parameters;
multi-start agreement is not a certificate. Uncertified optimization returns `fit_incomplete`.
This construction can be conservative. A bootstrap profile alternative is diagnostic until its
coverage, including endpoint search and nuisance boundaries, has been demonstrated.

Never interpret this within-family coverage as proof that one interval is the correct model family.
Two-part, neighboring-part, nonlinear-binding and RNA-hotspot controls must test its misspecification.
If the simple geometry family is inadequate, extend first to a finite two-part family using
`1 + binding * max_j overlap(part[j])`, still with shared binding. A variable-binding diagnostic is a
separate mechanism. Keep the original one-part and null alternatives inside
the enlarged family and recalibrate model selection. The family claim remains conditional; passing
a lack-of-fit check cannot establish it. If held-out residuals trigger a new family or a refitted `q`,
use fresh validation observations. Alternatively, predeclare the larger confidence family and retain
a training-selected normalized `q`. Do not stack overlapping probes or replace them with their union.

Exact enumeration is the correctness prototype. Resource exhaustion returns `fit_incomplete` and
uncertified/unresolved bounds. An optimization may prune candidates only with a likelihood bound
proving they cannot enter the retained set. No arbitrary endpoint grid, maximum useful fragment
length, regularization constant, or gene-wide probe imputation may narrow that set silently.

For larger contexts, retain all source/probe interactions reachable by the molecular lengths being
integrated. Context boundaries are determined by placement support, including RNA paths, not a fixed
flank distance. Merge interacting contexts or fit them jointly. If exhaustive support is infeasible,
the general map is not yet promotable.

#### Source-law dependence and the shared-amplitude test

Phase Laws uses a known map only to isolate law estimation; the production fit must never depend on
that oracle. Where justified, fit field shape with a length-conditioned spatial likelihood and free
per-length background intensities. Fit the shared RNA source law/abundances from selected spliced
observations conditional on each retained field/topology, then propagate the union of the resulting
yield envelopes. If conditioning is invalid, use a joint likelihood of the admitted unspliced
observations with source laws among its nuisance parameters. Independent RNA nuisance calibration
is allowed, with its uncertainty; the existing gDNA uniform-frame PMF is not known source truth
under every new capture/acceptance regime either.

Preserve the refusal of a gene's junction RNA as its capture witness. Do not add the spliced-law
likelihood to an objective that selects probe geometry/topology, or discard a DNA-admissible topology
because it fits spliced RNA less well. If a conditional law fit fails, retain the field alternative
and report the law/yield unresolved. A fully joint RNA-to-map update would reopen that refusal and
is outside this plan. As a test, freeze strand and source-law inputs, vary only spliced
junction allocations, and require the admitted field/topology set to remain unchanged. Final
uncertainty must cover the modular field/law coupling; do not call two plug-in errors independent
or reuse an observation as a new independent likelihood term.

Keep one shared amplitude until its assumptions are tested. Predeclare variable-binding, saturation,
GC/acceptance and local-baseline fixtures; compare shared versus independent context-amplitude
diagnostics on identical geometry, laws and observation support. Read held-out count prediction and
target-yield error as well as fitted-amplitude dispersion. Wrong shared binding can produce biased
endpoints, an empty set, or a narrow wrong interval; uncertainty does not automatically widen under
model misspecification.

Use the VCaP exome half's known panel and origin labels only in a separately labeled evaluator. Check
whether independent input/uncaptured coverage exists; local DNA copy number, mapping/GC acceptance,
sampling depth, uncertain probe geometry and source laws can all cause apparent amplitude dispersion.
Profile these nuisance contributions or report that they cannot be separated. Estimate sampling-only
dispersion with matched synthetic draws under the frozen common-amplitude model. Use a held-out
VCaP partition for the declared predictive check; it cannot also select a new tolerance or prior.
VCaP may reject the common model, but is not the sole arbiter of it. If it fails, reproduce the
mechanism in the synthetic fixtures, measure independent context amplitudes as the next isolated
change, and add a population prior only if the measurement justifies that extra mechanism.

### 5.4 Identification and unresolved topology

Keep a fitted effective amplitude/interval representation, not a reported physical global binding
constant and per-exon fraction. With multiple geometries, deduplicate only models with identical
relevant placement kernels, not models with identical base coverage.

For a linearized field, test whether every observation-null direction is also null for a target yield
using the observation matrix and target kernel. For the actual nonlinear best-part model, directly
profile the target yield over the retained model set. A small optimizer standard error is not a
replacement for this structural check.

Separate genomic parts do not tell whether they are connected by a transcript-space probe. Preserve
all split/join alternatives consistent with the evidence. Integrate each candidate RNA path under
those alternatives. Return an unresolved yield if the resulting interval is material. Never use
`all joined`, `all split`, a mean exon fraction, or the old additive junction price as exact truth.

A joined part must actually be contiguous on the proposed RNA path; it cannot skip unprobed RNA
bases between an interval endpoint and a junction. With several junctions, part identities must be
consistent across all affected isoforms. The shipped additive price is not necessarily either
endpoint of the physical topology range; establish extrema using the operator. A central probe
invisible to boundary crossings may still be visible in stranded contained-position evidence. Report
the resolving channel, and call it fully unobserved only when all usable channels fail.

No-component mixing rule: all enabled components of a connected EM locus must use the same kernel,
source-law, support, and scale contract. Do not apply absolute new yields to selected transcripts while
keeping old reference-normalized lengths for their competitors. If the locus cannot be represented
coherently, the full replacement remains diagnostic for it. A production approximation requires an
explicit error bound and its own approved/scored policy; lack of a bound blocks promotion rather than
creating an undocumented fallback branch.

## 6. Production architecture and integration seams

The following new module names are proposals to implement after their prototypes pass. Register them
in `calibration/_layers.py`; imports may go down or sideways, never up.

| Location | Responsibility |
|---|---|
| Layer 0, proposed `capture_model.py` | Immutable kernel, source-law provenance, uncertainty and support contracts; no calibration dependencies. |
| Layer 1, payload/substrate | Typed raw placement bank and views; no fitted beliefs. |
| Layer 2, proposed `capture_geometry.py` | Ownership, opportunity, incidence and conserved integrals; pure placement operators. |
| Layer 5, proposed `capture_fit.py` | Raw likelihood and supported field/functionals; no test of capture. |
| Layer 7, `calibrate.py`, `result.py`, `priors.py` | Orchestration and frozen outputs handed to downstream consumers. |
| Scorer/native scorer | Per-candidate capture kernel and path aggregation using that same frozen model. |

Avoid an oversized module split before measurement: the prototypes can use these responsibilities in
fewer files. Types belong below their consumers. `sweep.py` remains a mechanism-agnostic backbone.

### Existing files/functions to change at the indicated phase

| Concern | Existing seam | Required change |
|---|---|---|
| Raw evidence | `native/calibration/accumulator.cpp`, `accumulator.h`, native scanner export | Emit determinate placements once using existing accepted-deposit semantics. |
| Payload/cache | `scan_payload.py`, `scan_cache.py`, `calibration/substrate.py` | Validate and persist the new bank, exclusions and support metadata. Update structural schema/deposit digests; incompatible caches require rescanning. |
| Laws | `pipeline.library_fl_models`, `calibration/fl.py`, scanner/payload | Add selected spliced observation/candidate evidence and explicit source-law contract/provenance without relabeling census fields. |
| Fit orchestration | `calibrate(...)`, pipeline setup | Construct one frozen capture/law result used by calibration, geometry and scoring; no capture decision anywhere. |
| Transcript yields | `pipeline._setup_geometry_and_estimator`, `capture_eff_length.transcript_capture_eff_lengths` | Introduce physical yield computation on the correct template; keep existing output plain lengths distinct. |
| gDNA prior/yield | `calibration.priors.assemble_priors` | Freeze `gdna_count` in denominator experiments; integrate the actual conserved/routed support; remove the old span clamp only in the tested physical path. |
| Numerators | `scoring.FragmentScorer.from_models`, `native/scoring.cpp` | Carry reference/path metadata, add log capture before pruning and both ordinary/MM aggregation paths. |
| EM | `native/em_solver.cpp` | Consume matched unfloored nonnegative yields. Zero support cannot emit. Update contracts and invariance tests. |
| Composition | `_fit_gdna_hyperprior`, `_sweep`, `simplex_logodds.py`, `sweep.solve_chain`, native prior API, `native/psi_kernel.h` | Carry a dedicated prior-offset array and shift prior support; do not alter message opportunities in the same arm. |
| Origin counts | Native deposit/export, payload or deterministic replay, calibration raw-likelihood API, `priors.py` | Retain all contributors' actual deposit shares and path alternatives; aggregate origin responsibilities without recycling final EM posteriors into their own priors. |
| Results/reports | `calibration/result.py`, `_log_summary`, `cli.py`, `report/model.py`, report assets/substrate | Publish statuses/evidence and remove enriched-mode language when its reader is retired. Coordinate summary schema handling. |

The `_FinalizedChunk` buffer currently carries genomic start/footprint but not a direct reference ID.
For the final numerator path, add explicit `ref_id` through scanner output, buffer validation,
spill/load, and native scoring arguments. Keep it for every alignment alternative, including gDNA-only
ones. Do not infer chromosome from coordinates or an arbitrary transcript candidate.

The existing multimapper scorer keeps a best RNA hit per transcript but averages eligible gDNA hits.
Before calling the final likelihood normalized, define which alternatives are distinct latent paths,
which are duplicate alignments of the same path, and their observation probabilities. Sum distinct
path kernels and deduplicate duplicate representations. If changing the current merger is necessary,
measure that merger change separately with capture disabled before combining it with capture weights.

Do not add a permanent compatibility branch retaining the old reader. Retaining it in named A/B arms
is experimental isolation, with a deletion dependency in phase Release. No current API should accept
physical weights greater than one through a validator for clipped efficiencies in `[0, 1]`.

## 7. Composition and conserved-count corrections

### 7.1 Map-conditioned density prior

Keep composition frozen during the denominator experiments. For the later composition
arm, use one pooled residual field:

```text
capture_offset[o] = log(captured_gdna_incidence_opportunity[o] / plain_gdna_incidence_opportunity[o])
residual_log_density[o] = log_gdna_density[o] - capture_offset[o]
```

Fit the landscape on the residual coordinate and shift it into each object's solve. Carry a dedicated
prior-offset array through `simplex_logodds.py`, `sweep.py`, the native prior interface and
`psi_kernel.h`; update the prior lattice window to include shifted support. Do not overwrite
`eff_global`: it also drives message opportunities, masks and search brackets. Handle zero
plain opportunity structurally, without a ratio. When capture is disabled, every offset is zero; no
extra exon-class population is fitted. Map uncertainty must widen/mix the prior rather than become
a deterministic pseudo-observation.

Use independent molecule splits for an initial clean experiment: learn the map, the residual prior
distribution, and every data-fitted nuisance entering that prior on one split; apply the resulting
prior to the other, then swap and score. Each held-out observation contributes its likelihood once.
Cross-fitting only the map while fitting the residual landscape on held-out observations does not
remove their reuse. A joint model may later recover the split's efficiency loss. If existing empirical
Bayes reuse is retained, label and measure that approximation; do not claim it is an exact joint model.

The correction for the reference prior entering twice remains a separate C++ arm. Test it with the
map-conditioned prior as its proposed partner. Preserve the previous refusals; the unconditional
class split is not reinstated.

### 7.2 gDNA conserved count

An additional count-frame defect must not be hidden by a new denominator. `priors.py` currently
converts gDNA boundary incidence using the observed mixture's `boundary_mass_per_crossing`. The
appropriate expected share for gDNA depends on its own law and capture:

```text
expected_gdna_share[o] = captured_conserved_opportunity_gdna[o] / captured_incidence_opportunity_gdna[o]
```

First audit this ratio against per-origin conserved truth. Multiplying it by an estimated incidence
is a moment approximation, not automatically a conserved realized origin count. Keep it as a named
diagnostic arm only.

Before implementing the coherent count arm, add an all-contributor deposit/path ledger or a
deterministic replay pass that retains actual object shares, locus routing and latent path alternatives
for every calibration contributor, including spliced, deferred and implicit-path observations.
The clean evidence bank excludes these observations and lacks object shares; aggregate composition
arrays cannot reconstruct them. Add a joint raw-likelihood responsibility API that marginalizes path
alternatives and returns origin/path responsibilities on this same deposit domain. Do not use final
EM posteriors to construct the priors consumed by that same EM fit.

For the coherent count implementation, aggregate per-molecule gDNA responsibilities times each
molecule's actual deposited shares, jointly with the corresponding RNA responsibilities. This
preserves total observed mass by construction and avoids using RNA's share for gDNA. Responsibilities
come from the raw observation likelihood; do not manufacture the count as `rho_off * model_yield`.
Retain the independent `P - O` calibration and `O - S` assembler comparisons to attribute the change.
Score this count correction separately from the map-conditioned prior and from the denominator change.

## 8. Phases and concrete completion gates

Use named phases below in code and reports; do not invent new numbered rule labels that violate the
repository's documentation gate. Each phase has a baseline arm, a candidate arm, and a frozen-input
manifest. A passing phase is not permission to silently change another mechanism.

### Phase Audit — repair and freeze the measurements

First establish four unambiguous artifact arms (labels for the prototype manifest, not additional
values accepted by the current CLI):

| Arm | Purpose |
|---|---|
| `release_tag` | Actual `v0.7.1`, the owner's end-to-end release benchmark. |
| `current` | The frozen current dirty tree, including the defect and its present reference reader. |
| `instrumented` | That tree plus evidence acquisition, with behavior disabled; must reproduce `current`. |

Every later mechanism is A/B'd against `instrumented` and read against `release_tag`. (A `decision_veto` arm was
planned here and is abandoned by the owner's ruling.)

Tasks:

1. Create an outside-tree prototype run directory and manifest. Copy the review audit out of any
   temporary location into this reproducible artifact set; do not depend on `/private/tmp` surviving.
2. Build `v0.7.1` in its own worktree, environment and native build. Its index format is 7, versus 8
   in the current tree; rebuild from the same frozen FASTA/GTF and semantically matching options.
   Do not share binary indices/caches or regenerate the BAMs/truth with the tag's simulator. The tag
   has neither today's scan-cache API nor today's `scripts/design` instruments. Run its own CLI
   with `--threads 1 --assignment-mode fractional --seed` using a recorded seed, after checking
   that tag's help. Record resolved configuration, loaded Python/native paths and build identities.
3. Export tagged `quant.feather`, `gene_quant.feather`, `nrna_quant.feather`, `summary.json` and the
   equivalent current outputs. Score neutral exports with one external adapter; never load the two
   native builds into one interpreter. Validate the adapter on current outputs against today's
   instruments. Join transcript identifiers/coordinates semantically; report unmatched estimate
   mass and synthetic entities rather than dropping them in a truth-left-join. Reconcile transcript,
   gene, synthetic-RNA and gDNA pools, including intergenic mass, before comparing release errors.
   Report unavailable historical calibration fields as unavailable, not converted pseudo-truth.
4. Correct the descriptive probe fraction to use union bases. Preserve separate probe-part identities
   for capture physics. Produce raw count/exposure tables, not just ratios.
5. Report, per geometry/class/condition, total eligible objects, finite ratios, positive/zero ratios,
   zero/zero ratios, excluded objects, and associated RNA/gDNA truth mass in the scorer only.
6. Reproduce the reviewed RNA-short ON audit: old class counts `834/1709/3166`, union-corrected
   `799/1744/3166`; 206 old fractions above one; only 699 finite truth-gDNA ratios in the original
   834-object fully-probed class. If source hashes changed, explain differences rather than forcing
   these historical numbers.
7. Read existing `kDnaIntronExon` versus `kDnaIntronic` genome-pooled length shapes as a cheap
   diagnostic. Those pools have different selection rules and lack structural-window labels;
   their histogram ratio is neither the proposed conditional test nor a capture decision.
8. On the VCaP transcriptome half, report continuous crossings of the proposed eligible boundaries
   by raw observation class/strand, using origin labels only in evaluation. Distinguish explicit
   splice, ambiguous insert path, annotated continuous RNA and unknown transcription; raw crossings
   are not proof that a mature transcript produced a continuous molecule.
9. Re-record current-tree calibration, ruler, priors and end-to-end metrics on the test panel, both
   gap arms and ladder; compare tagged and current end-to-end results on every required panel and
   real input. Use truth where available; plasma outputs alone do not certify quantification accuracy.
   Record B_match's source laws, support verdict and numerator approximation. Its historical gain is
   a qualified diagnostic, not a guaranteed attainable ceiling for the learned likelihood.

Artifacts: `manifest.json`, `census.tsv`, `coverage_audit.json`, `baseline/`, `comparison_contract.json`.

Exit: all count partitions reconcile, coverage lies in `[0, 1]`, zero denominators remain visible,
and every benchmark row has immutable provenance. No estimator is fitted yet.

### Phase Evidence — Python preflight, then native observation acquisition

Tasks:

1. Implement the raw placement ledger and once-per-fragment `SplicedAdmissionCounts` first in the
   Python accumulator specification, including matching strand strata and disjointness checks.
2. Implement the exact cell enumerator and ownership rule using synthetic
   offered-event streams. Stop for repair before native work if the known contrast fixture loses
   its information through conditioning. This is the inexpensive Python preflight.
3. Implement/export both banks in the native worktree, using the existing scanner's offered-fragment
   construction. `Accumulator.deposit` accepts already constructed extents, strands, introns and gap
   hypotheses; it is not a BAM reader. The current buffers cannot reconstruct all those raw inputs.
   If an offered-event trace is used for interim replay, verify scanner/spec input parity first.
   Simulator truth and selected drain paths cannot substitute for observed inputs.
4. Add schema validation, spill/reload and cache behavior coverage. `payload_schema_digest` can see
   nested field names, but today's `deposit_digest` hashes only top-level ndarray attributes. Add
   canonical recursive hashing or explicit bank exports and extend the fixed behavior fixture to
   exercise each new field/outcome. A deliberate new-bank-only behavior change must alter the
   digest. Certify layout/eligibility generation at the scanner seam too.
5. Produce the real-bank observability census and full scanner-to-cell contrast fixture from
   section 4.6. Measure memory/spill behavior; compare `instrumented` to `current` with the gate disabled.
6. Rebuild only affected scan caches in separate destinations; never overwrite a baseline cache with
   a new schema under the same claimed artifact identity.

Artifacts: placement banks, cell/exposure tables, exclusion tables, producer hashes, parity report.

Exit: native/spec parity; one event per admitted molecule; canonical integer banks identical across
workers; existing floating mass/reciprocal banks within their established derived tolerances.
At fixed worker count, the paired instrumentation change leaves existing fields unchanged. Cache
rejects missing/incompatible evidence; zero-count opportunities are
enumerated; the full raw-bank fixture retains its expected contrast. Stop if adding the bank changes
ordinary deposits or drain behavior. Record structurally lost information before the field is fitted.

### Phase Geometry — certify the physical operator and support

Tasks:

1. Implement part-based weights, candidate molecular paths, source-law injection, and incidence versus
   conserved integrals in the prototype.
2. Build the independent enumerator and all geometry fixtures in section 9.
3. Audit the actual EM routing/deposit support per locus. Compare union and conserved domains and
   retain both diagnostics; use the matched one for a claimed normalizer.
4. Define log-kernel and log-yield APIs on a common scale, including exact zero support.

Exit: every per-placement and per-object identity holds; the six-base covariance fixture gives
`3 + binding`, not `3 + 0.75 * binding`; two separate 40-base overlaps contribute 40, not 80;
reference/template endpoints and shared-locus counts are correct. Derive floating error tolerances
from accumulation/error analysis and test deliberate broken operators.

### Phase Laws — source-law estimation in isolation

Tasks:

1. Add the spliced observation/candidate bank or deterministic replay, source-law data structures,
   and fixed-map joint likelihood from section 5.1. Verify observation normalization and evidence
   sufficiency before fitting mixtures; the global length histogram is insufficient.
2. Compare oracle versus estimated source laws while holding the oracle map and all other inputs
   fixed. Separately preserve a measured-census-as-source negative control.
3. Validate origin/class-specific predicted observed length histograms and structural support.
4. Profile sparse/off-target-starved cases and carry uncertainty forward.

Exit: the selection fixture passes; uncertainty covers repeated draws; unsupported tails remain
unsupported; laws are not silently borrowed from the simulator in the production fit.

### Phase Field — infer supported response and yield functionals

Tasks:

1. Implement the direct contiguous-functional control using oracle origins, then the raw-strand
   likelihood version. Keep truth available only to evaluation.
2. Implement one interval per context with a shared amplitude and scalar outer profile from section
   5.3. Start with structural equivalence and finite-data profile envelopes. Calibrate a retained set
   before claiming confidence coverage; do not require the advanced construction for diagnostic runs.
3. Evaluate source-law uncertainty separately first, then use the modular production law/field
   procedure in section 5.3 with no oracle initialization requirement. Keep both raw strand channels
   and leakage; conditional spliced-law fits cannot select capture geometry or junction topology.
4. Add the two-part/neighbor extension only when the one-part model's falsifiers require it. Preserve
   censored endpoint and junction-topology alternatives; incomplete search cannot certify bounds.
5. Report per-component yield intervals and supported truth mass, not just apparent interval/coverage
   recovery. No unconditional borrowing across a gene's exons.
6. Test shared binding against the independent-context diagnostic on frozen simulated fixtures and
   the declared VCaP holdout. Keep geometry, baseline, law and acceptance uncertainty in the comparison;
   a population prior is a later separate experiment, not an automatic response to raw dispersion.

Exit: observationally equivalent models retain their different RNA yields; held-out count predictions
are calibrated; supported yields meet their accuracy contract. Unidentified components/loci are listed
explicitly. If coverage is inadequate for the release, this phase has found a remaining information
problem and must not be presented as a complete replacement.

### Phase Yields — denominator-only diagnostic

Tasks:

1. Freeze composition, count priors, source/scorer laws, candidate scores and routing.
2. Substitute both RNA-template and gDNA-locus denominators using one common scale and correct support.
3. Keep the full shipped pipeline as an external baseline. Inside the frozen-input comparison,
   compare the old denominator implementation evaluated under those same frozen laws (the law-only
   control), known-map matched denominators, and learned-map matched denominators. Separately score
   the law-only control against shipped behavior to attribute law injection. Compare the historical
   B_match only as an additional qualified diagnostic, not as an isolated denominator change.
   Clearly labeled point/profile candidates may run as diagnostics before uncertainty calibration;
   they do not acquire confidence or release claims by producing good denominator scores.
4. Read actual denominator vectors back from the estimator to prove the hooks fired. Compare class
   mean log length/yield ratios and within-gene spread beside the transcript/gene/pool table.

Exit: frozen-input hashes agree; denominator identities pass; learned/oracle differences are
attributable to field/law estimation rather than support or scale. Do not call this the final model.

### Phase Likelihood — matching numerators and observation paths

Tasks:

1. Carry explicit reference and full candidate-path metadata through the buffer/scorer interface.
2. Define and test alternative-hit aggregation with capture disabled first. Remove duplicate path
   representations without dropping distinct placements.
3. With known map and laws fixed, add candidate capture terms before pruning, gDNA competition and
   hit merging in both ordinary and multimapper scoring paths.
4. Integrate over any remaining implicit splice/insert paths with their own lengths and weights.
5. Repeat with the learned field/laws and shared uncertainty/support statuses.

Exit: candidate likelihoods normalize on their actual domains; the implicit-splice fixture reproduces
the physical posterior; duplicate alignment representations do not change a probability; unique
placements reduce to the single-path formula; the same kernel/law/support hashes reach scorer and
denominator builder. Measure changed pruning/routing separately; do not demand them stay unchanged
after a correct numerator is enabled.

### Phase Composition — map-conditioned prior and count frame, separate experiments

Tasks:

1. Implement and test the dedicated prior-offset API and shifted prior lattice support. Test the
   residual-density prior from section 7.1 with an oracle map, then with independently trained maps
   and residual prior distributions. Keep message policy and count conversion fixed.
2. Test the prior-accounting correction separately with that prior; retain the analytic gate.
3. Implement the all-contributor deposit/path ledger or replay and joint responsibility API required
   by section 7.2. Verify their raw mass/path reconciliation, then test conserved gDNA counts separately
   with the other components fixed.
4. Combine only individually attributed changes and remeasure their interaction.

Exit: calibration object and prior-count errors pass per-stratum bars; OFF offsets are identically
zero with one pooled population; no new g00 pseudo-mass; observed mass is conserved by origin-aware
deposits; no duplicate use of witness likelihood factors.

### Phase Release — integrate, delete and verify

Tasks:

1. Select the smallest complete passing scope against both the tagged release and current tree.
   State which mechanisms are proposed; do not claim a full likelihood replacement from a partial win.
2. Port passing prototype code into the designated layers, register modules, and remove experimental
   hooks from the production path. Source comments cite tests, not this document or prototype files.
3. Publish the new statuses, evidence provenance, support coverage and uncertainty in results/reports.
   Coordinate the summary schema identifier and its reader; incompatible artifacts receive an explicit
   diagnostic. Do not add a permanent old-reader fallback solely for compatibility.
4. For a complete likelihood replacement, retire `located_enriched_mode`, capture-only basin/member
   census consumers, posterior-mean efficiency machinery, additive junction price, and reference/member
   result fields. Search every caller first: retain any landscape functionality still used for ψ.
5. Remove the raw-density mixture role only after the residual prior has replaced it. Preserve the
   necessary conserved geometry, source laws, raw evidence and structural-support tests.
6. Run full checks, fresh pinned panels, and real-library test inputs. Promote conclusions to their
   permanent documentation homes only with the measured implementation. Do not commit automatically.

Exit: no unresolved required gate, no undeclared approximation, no oracle dependency, no incompatible
mixed yield scales, reproducible reports, and all applicable tests/instruments pass.

## 9. Mandatory falsification matrix

Build synthetic fixtures by changing one cause at a time. Cases marked negative controls must actually
fail a deliberately broken estimator; record the mutation and observed assertion failure. Synthetic
numbers below define fixtures, not claims about real libraries.

| Named case | Construction | Required read-out |
|---|---|---|
| Coverage union | Duplicate/overlap BED blocks in the evaluator. | Coverage bounded by one; physical part identities retained separately. |
| Release adapter | Score the same current outputs through the neutral adapter and current instrument; add an unmatched estimate entity. | Metrics/pools reconcile; unmatched mass is reported, never silently lost; tag and current use isolated builds and indices. |
| Evidence no-op | Current versus instrumented trees with the gate disabled, same inputs and workers. | Existing arrays and end-to-end outputs unchanged; only evidence/diagnostics added. |
| Bank-only cache mutation | Change a new bank's eligibility or layout field without changing legacy fields. | Behavior digest changes; old evidence cache is rejected. |
| Zero denominator | Positive boundary/zero intron and zero/zero cases. | Counts retained; no silent deletion, pseudocount or infinite ratio used as precision. |
| Repeated boundary | A molecule crossing several short regions. | One event; incidences remain diagnostic; unchanged by duplicate counting opportunities. |
| Raw-bank observability | One known structural window with both boundary and interior opportunities, scanned through the real ledger and cell builder. | Expected molecules/exposures survive length, layout and strand conditioning; no false success from testing tails alone. |
| Thread/spill identity | Same scan at existing worker-count fixtures and spill sizes. | Integer banks agree exactly; old floating banks obey established derived tolerances. |
| Reversed laws | Equal laws, RNA 75/gDNA 250, and reverse; zero and positive continuous RNA. | Conditional null calibration independent of the gap; source-law fit assessed separately. |
| Finite RNA reach | Short introns, gene termini, overlapping spans, same-strand and opposite-strand genes. | Correct support signatures/opportunities or explicit exclusion; never generic unbounded RNA reach. |
| Retained/missing annotation | Annotated and unannotated retained intron, readthrough and shadow transcripts. | Eligibility changes where annotated; unsupported assumptions/failure rates exposed where not annotated. |
| Nonlinear hotspot | Local RNA peak near a boundary, plus high depth. | Model limitation/misspecification exposed; no claim of universal gradient protection. |
| Shadow by gap | Fixed shadow contamination with both law gaps. | Old plug-in excess has the predicted sign reversal; conditional comparison avoids that rho/nu plug-in mechanism. |
| No full-coverage anchor | All targets partial, including infinite-depth expectations. | No falsely identified physical global binding/probed-fraction decomposition. |
| Invisible central probe | Probe farther from boundaries than gDNA support reaches, with and without usable contained-strand evidence. | Report which channel can see it; no inference from absent boundary contrast alone. |
| Same summary, different yield | Length-1000 template, gDNA 100, RNA 250; part `[110,230)` versus `[400,520)`. | Old summaries agree; overlap integrals for RNA are 20460 versus 30000. Summary-only confidence cannot select one. New positional evidence may distinguish them if it actually resolves the field. |
| Part competition | Two separate 40-base overlaps versus one connected 80-base part. | Best-part response distinguishes the layouts; no sum/union substitution. |
| Intronic overhang | Crossing `[-50,50)`, part `[-200,20)`, plus mirrored right edge and controls. | Overlap 70, never 220; the physical operator handles the complete placement. |
| Conserved covariance | Three 2-base pieces, fragment length 4, one-base probe at coordinate zero. | Yield `3 + binding`; incidence-average shortcut `3 + 0.75 * binding` fails. |
| Shared/neighbor locus | Shared objects, short outside neighbors, multi-block loci and reference ends. | Each count's routed/conserved domain matches its denominator. |
| Source/census law | Two-length selection example from section 5.1. | Selection applied once; source and observed laws remain distinct. |
| Structural zero | Transcript shorter than supported RNA lengths; zero counterfactual gDNA support. | Exact zero/unidentified support, no rescue by numeric PMF floors. |
| Support versus depth | Truly disjoint source support beside the current 50–500 gap laws and progressively thinner samples. | Structural nonidentification separated from rare-tail/high-variance transfer; different means do not imply disjoint support. |
| Implicit splice | Candidate weights 701 versus 501, with distinct molecular lengths. | Capture enters candidate odds; full posterior matches independent enumeration. |
| Multimap/path duplicate | Distinct alternative placements plus duplicate encodings of one path. | Alternatives summed under the stated measure; duplicate encodings do not increase likelihood. |
| Exon targeting alternates | Only alternating exons in a gene targeted. | No automatic all-exons-probed borrowing; uncertainty follows local evidence. |
| Binding/GC/nonlinearity | Variable amplitudes, saturation, GC-dependent acceptance; both gap directions. | Response and transfer errors attributed; unsupported shape assumptions do not masquerade as fractions. |
| Binding confounding | Shared binding with varying depth/local DNA baseline, then variable binding with matched depth. | Raw amplitude dispersion is not mistaken for binding heterogeneity; misspecified shared fits may be narrow/wrong and fail predictive/yield gates. |
| Sparse plasma | Thin to measured per-library usable gDNA depths, with RNA support appropriate to plasma tests, and a panel whose captured transcripts are a small minority of the library. | Per-locus field coverage and quantification cost over independent thinnings; no library-level reading switches capture off. Previously: separate admission power, contrast power, field coverage and quantification cost over independent thinnings. |
| Prior/count isolation | Truth map with current versus residual prior, then mixture-share versus origin-share counts. | Each change's calibration and pool effect measured independently. |
| Oracle isolation | Run fit process without truth arrays, probe path, or simulator parameters. | Identical estimator output; only the evaluator changes when truth is supplied. |
| RNA-to-map isolation | Freeze strand and source-law inputs; alter only spliced junction allocations. | Field/topology set unchanged; RNA law fitting cannot silently become a gene's own junction-capture witness. |

Measure every relevant row by strand specificity and capture condition, with g00 separate. Include
the ss 0.70 stress cases. Report per-geometry coverage and error as well as overall transcript, gene
and pool errors. A favorable total cannot conceal a class for which the model is unidentified.

## 10. Validation commands and integration tests

The commands below are for the current tree's `rigel` environment. The tagged baseline uses its own
environment/CLI and rebuilt index as described in Phase Audit; these instruments do not exist at the
tag. Run native candidate worktree commands in a separate environment as well; an editable install
must not replace the frozen baseline's imported build. Inspect each tool's current `--help` before a long
campaign; do not invent unsupported `--arm` values. Add prototype arms through an explicit isolated
wrapper or the established harness, and verify hook counts and read-back values.

```bash
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate rigel
export OMP_NUM_THREADS=1
python scripts/design/preflight.py --full
```

Use a new run directory for outputs. These are existing instrument examples for a single ladder row;
expand to all declared conditions and panel paths in the campaign manifest:

```bash
CAPTURE_RUN=$(mktemp -d /tmp/rigel-capture.XXXXXX)
python scripts/design/quant_accuracy.py --arm base \
  --conditions gdna_g50_ss_0.99_nrna_mid_capture_on \
  --set scan.total_threads=1 --set em.assignment_mode=fractional --jobs 1 \
  --out "$CAPTURE_RUN/baseline.jsonl"
python scripts/design/calibration_vs_oracle.py \
  --conditions gdna_g50_ss_0.99_nrna_mid_capture_on --jobs 1 \
  --json "$CAPTURE_RUN/calibration.json"
python scripts/design/prior_vs_oracle.py \
  --conditions gdna_g50_ss_0.99_nrna_mid_capture_on --jobs 1 \
  --json "$CAPTURE_RUN/priors.json"
python scripts/design/ruler_vs_truth.py --panel ladder \
  --condition gdna_g50_ss_0.99_nrna_mid_capture_on --scale
```

Order expensive runs: test chromosome, RNA-short defect and other gap rows, ladder, then real test
libraries. Re-record a same-session baseline for each change. Run one genome-scale job at a time.
Use detached resumable execution for long campaigns with logs, exit status, manifest and completion
records. A missing/failed row is incomplete, never an implicitly passing row.

Existing test groups to extend/run at the matching phases:

```bash
python -m pytest tests/native/test_accumulator_spec.py \
  tests/native/test_accumulator_native_parity.py \
  tests/native/test_accumulator_drain.py tests/native/test_gap_hypothesis_arbitration.py \
  tests/native/test_conserved_mass.py \
  tests/test_accumulator_payload.py tests/test_scan_cache.py \
  tests/test_scan_order_independence.py -q
python -m pytest tests/calibration/test_layering.py \
  tests/calibration/test_result_schema.py tests/calibration/test_fl.py -q
python -m pytest tests/calibration/test_effective_length.py \
  tests/calibration/test_capture_eff_length.py tests/calibration/test_priors.py \
  tests/calibration/test_prior_units.py -q
python -m pytest tests/test_buffer.py tests/test_second_pass_scoring.py \
  tests/test_estimator.py tests/test_em_pseudocounts.py -q
python -m pytest tests/test_summary_report.py tests/test_report.py \
  tests/test_docs_boundary.py -q
```

Add new focused tests for capture event ownership,
raw-likelihood fitting, source-law selection, part geometry, yield identification and candidate
normalization. These files do not exist merely because the plan names their responsibility. Pair
each important protection with a deliberate broken implementation in the prototype test report.

After native changes rebuild in the isolated environment, then run relevant native and Python tests.
For final integration:

```bash
pip install --no-build-isolation -e ".[dev]"
python -m pytest tests/ -q
ruff check src/ tests/ scripts/
python scripts/design/preflight.py --full
```

Re-derive the collected test count; the historical 2650-pass count is not permission to ignore new
tests or changed collection. Do not update goldens to make an unexplained behavioral change pass.
When a new instrument is promoted, update its script index, self-test registration and existing
instrument-contract tests. Do not add permanent citations to this sandbox plan.

## 11. Artifact, reproducibility, and stopping contract

Every phase writes machine-readable records beside its human-readable verdict. Minimum manifest:

```text
run identity: commit, dirty-diff hash, environment, native build identity
inputs: BAM/index/annotation hashes, scan schema/deposit hashes, panel generator identity
method: phase, arm, kernel/law/support/eligibility hashes, assumption set, protocol identity
randomization: independent seeds, fixed replicate counts, split assignment procedure
outputs: paths, hashes, hook counts, read-back vector hashes, timings, peak memory, exit codes
assessment: per-stratum results, unsupported classes/mass, gates passed/failed/incomplete
```

Keep fitting inputs physically separate from evaluation-only truth/panel inputs. Pass typed fit
inputs, not a general benchmark context containing the oracle. A test must prove oracle information
cannot reach the production estimator through configuration, globals, cache metadata or read names.

Suggested outside-tree prototype layout, with names chosen once for the campaign:

```text
capture_learning/
  protocol.json
  manifest.json
  evidence.py             raw-bank export/view and cell geometry prototype
  geometry.py             physical and conserved placement integrals
  fit.py                  supported law/field fitting and uncertainty
  checks.py               independent enumerator and deliberate breaks
  run_condition.py        explicit arm hooks and frozen-input assertions
  summarize.py            per-stratum reports; refuses incomplete/mixed contracts
  baseline/
  cases/
  results/
  failures.jsonl
```

For a failed gate record the exact condition, source/input hashes, expected and observed values,
failure direction, and the next discriminating measurement. Do not rerun a large panel without a
specific changed mechanism or a previously missing condition. For a resource failure preserve the
partial search/report as incomplete; never call its uncertainty bound certified.

The following block promotion:

- Any hidden oracle/panel dependence or inconsistent count/opportunity selection.
- A null error guarantee asserted outside its recorded assumptions.
- A capture test using repeated boundary incidences as independent molecules.
- Source/census law confusion or a numerical floor creating support.
- Different component scales or mismatched kernel/denominator support within an EM locus.
- An identified point yield claimed for observationally equivalent fields with materially different yields.
- A nuisance/fallback failure exceeding the protocol, including the measured captured-library collapse.
- A claimed full replacement with unmeasured numerator, composition, or conserved-count effects.

## 12. Completion checklist for the executing agent

The decision-only release candidate is abandoned (owner ruling, 2026-10-05); the items below are one milestone.

- [ ] Tagged 0.7.1 and current baselines built separately and scored through a reconciled adapter.
- [ ] Measurement labels/zero-denominator reporting corrected and reproduced.
- [ ] Raw molecule bank, spliced strand counts, disjoint ownership, full support signatures and
      cache invalidation tested.
- [ ] Full likelihood only: best-single-part geometry, count frames and routed support certified.
- [ ] Full likelihood only: source laws estimated with correct provenance and structural support.
- [ ] Full likelihood only: response/functionals assessed with held-out data and identified-set bounds.
- [ ] Full likelihood only: unresolved junction topology and missing support retained honestly.
- [ ] Full likelihood only: numerators, alternative paths and denominators share the kernel contract.
- [ ] Full likelihood only: map-conditioned prior, prior accounting and origin counts measured separately.
- [ ] Obsolete capture machinery deleted only after its consumers have complete replacements.
- [ ] Full tests/instruments and pinned panel/real-input checks pass; no source cites this plan.
- [ ] Final report names the completed deliverable and every remaining information/model limitation.

## 13. Evidence locations

These are provenance and diagnostic inputs, not runtime dependencies:

- `docs/dev/CAPTURE_CLEAN_SLATE.md` — the 2026-10-05 revised design, retained measurements and four
  questions answered by the recommendations in section 1. Its compressed claims are qualified here.
- `docs/dev/FRAGMENT_LENGTH_POSTMORTEM.md` — defect, refused mechanisms and campaign history.
- `~/Downloads/rigel_runs/prototypes/2026-10-03_fl_arms/VERDICT.md` — per-stratum Q3/Q4 results.
- `~/Downloads/rigel_runs/prototypes/2026-10-04_clean_slate/` — census/mass scripts and complete outputs.
- `~/Downloads/rigel_runs/prototypes/2026-10-02_stage3/` — exact operator fixtures, law controls,
  matched-denominator diagnostics and support audits. Its README includes historical generator notes;
  use each run's recorded generator identity rather than treating old notes as current panel status.
- `~/Downloads/rigel_runs/prototypes/2026-10-04_capture_compromise/bed_geometry.py` — incidence integral
  prototype; not a proof of conserved-share equality.
- `/private/tmp/rigel_capture_review_audit.py` — review census audit, if still present. Copy it into the
  campaign artifacts before use; if missing, reconstruct the union/finite-ratio audit from the named
  census inputs and verify the historical values in Phase Audit.

The postmortem's mass diagnosis remains supported: on the defect row the bulk contains 99.9% of
estimated gDNA mass while the false reference contains 0.024%. The implementation replaces the
member-census failure with raw-likelihood evidence; it does not elevate observed-mass medians into a
universal information estimator or use agreement with the shipped reference as truth.
