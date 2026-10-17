# Brief: can STAR report both readings of an ambiguous read pair?

2026-10-16. Written by the laptop session that owns Rigel's paired-end design, for a new Claude Code session on the
SLURM cluster: you. Read all of it before running anything; it is your only context. The owner is reachable in
this session. "Template" below means one read pair.

## 1. The question

Rigel (`rigel-rnaseq`) quantifies annotated transcripts from paired-end RNA-seq while separating genomic-DNA
contamination (gDNA) from RNA. A spliced alignment is the strongest evidence that a fragment is RNA, because DNA is
never spliced. So when the aligner reports a gDNA template as spliced, Rigel certifies it as RNA: the release-blocking
false positive called a **splice artifact**. The reverse also happens: a genuinely spliced RNA template is reported as
a contiguous alignment (to a processed pseudogene, or extended into its own intron) and is treated as possibly gDNA.

In both cases the template's bases usually fit two readings almost equally well, a **spliced reading** and a
**contiguous reading**, and STAR reports only one. Rigel's draft design compensates inside its BAM scanner: it
re-derives the other reading itself, moving bases across junctions, reconciling mates and proposing contiguous
alternatives. The owner's objection (2026-10-16): that makes Rigel an alignment-fixing tool. If STAR reported both
readings as alternative alignments of the same template (multimappers, with NH/HI and per-alignment NM/MD), Rigel's EM
would weigh them with the likelihood it already has, and the scanner would stay a reader.

**Your question.** Can STAR 2.7.11b be configured to report, for each template, every alignment whose score is close to
the best, including the contiguous reading of a spliced alignment and the spliced reading of a contiguous one? At what
cost? Where a single run cannot do it, find exactly why in STAR's source, and price the fallback: a second, unspliced
alignment of the same reads, merged by read name.

You are not asked to change Rigel. You are asked for evidence a design decision can rest on.

## 2. What Rigel does with an alignment

You need this to judge which outputs matter.

- **Input.** A name-collated paired-end BAM. One alignment of a template is a **hit**: both mates, with the same HI.
  NH is the number of hits. Rigel reads the flags, NH, HI, NM, and XS (from `--outSAMstrandField intronMotif`) for a
  junction's strand. Unmapped and QC-fail records are skipped. Chimeric alignments are not quantified.
- **Scoring one hit.** Every annotated transcript the hit overlaps is a candidate. Each base outside that transcript's
  exons (overhang) multiplies its likelihood by 0.1, and so does each NM mismatch. The fragment length is computed
  along the transcript. A gDNA candidate exists only if the hit is unspliced.
- **Pruning.** A candidate below 10⁻⁴ of the template's best likelihood is pruned. So an alternative that carries 4 or
  more extra mismatches contributes nothing (0.1⁴ = 10⁻⁴). This bounds the useful score range.
- **Several hits.** Each transcript keeps its best hit, and the gDNA term combines the unspliced hits. A same-locus pair
  of hits (spliced S, contiguous C) is therefore already scored sensibly: RNA takes the better reading per transcript,
  and gDNA takes C.
- **Unannotated junctions.** A hit whose junction is not annotated is withheld entirely. This is the owner's ruling:
  Rigel does not assemble transcripts. So a gDNA template that STAR reports only with an unannotated junction is lost.
  The same template with its contiguous reading also reported survives as an unspliced hit.
- **Soft clips.** Clipped bases are not penalized today. That is a known gap: an alternative that wins by clipping is
  not scored fairly, so count the alternatives that differ only by clipping.
- **Calibration.** Its first pass learns from unique templates (NH = 1), so every template that gains an alternative
  leaves that pass. Count those templates.

The situations the draft design would repair itself are listed below. Test whether STAR can cover each one.

1. **Short side.** A spliced alignment with only a few bases across the junction, whose contiguous reading fits as well
   or nearly as well.
2. **Mates that disagree.** One mate is spliced across a junction, and the other runs a few bases into the same intron,
   or is clipped there.
3. **A remote copy.** A gDNA template from a processed pseudogene is reported spliced at the parent gene. The RNA-side
   mirror image is a spliced template reported contiguous at the copy.
4. **An unannotated splice of gDNA.** Typically a sequencing error near a read end, plus a chance match of the end
   bases further downstream at a GT-AG or GC-AG motif.

**The owner's rulings that bear on this work.**

- Alternatives with more mismatches must still be reported, because sequencing errors are a common reason a fragment
  aligns wrongly.
- Every alternative is scored with its own mismatches.
- Rigel will not propose spliced readings for unspliced reads. That is the aligner's job, and part of what you are
  measuring.

## 3. What is known, and the hypotheses to test

Measured earlier, on the laptop:

- **The VCaP mix** (real data; STAR 2.7.11b; the production parameters in §4). The DNA library has 4.71 M templates,
  and 10,164 of them have a sequenced junction, every one an artifact.
  - Among the 2,944 (template, junction) units at DNA-majority junctions: 828 annotated, 692 of them with h = 0; and
    2,116 unannotated (72 %), mostly long-anchor. Here h = 0 means the reference at the spliced placement equals the
    reference at the contiguous placement over the moved bases, so the bases cannot tell the two readings apart.
  - Annotated artifacts with a short side of 1–10 bases: 96 % have h = 0, and 81 % have the short side on the
    acceptor.
  - Annotated artifacts with h ≥ 1 are mostly remote copies: 71 % fit contiguously near-perfectly, 996 of those 1,056 at
    another locus (minimap2 check).
  - The unannotated DNA artifacts are 62 % GT-AG and 31 % GC-AG.
- **minimap2 on the simulated panel** produced mostly unannotated long-anchor artifacts, and 60 % of them were
  multimappers whose secondary was often the true origin. The aligner matters a great deal, so nothing measured on
  minimap2 transfers to STAR.
- **The production parameters** include arriba's fusion settings, e.g. `--alignSJstitchMismatchNmax 5 -1 5 5`. STAR's
  default is `0 -1 0 0`, and the production value admits mismatches near non-canonical and GC-AG junctions.

Hypotheses, from a second-hand reading of STAR. Confirm or refute each from the source and the data:

- **Selection.** `--outFilterMultimapScoreRange` (default 1) reports alignments within 1 point of the best. A mismatch
  costs 2 points (+1 becomes −1). `--sjdbScore` (default 2) adds 2 to an alignment across an annotated junction.
  - So at an annotated junction with h = 0, if both readings span the same genome, the spliced reading wins by exactly
    2, and the default range never reports the contiguous one. The two spans are equal when the other mate lies beyond
    the junction.
  - The log2 genomic-length penalty (`--scoreGenomicLengthLog2scale`, −0.25) pulls the other way when the spliced
    reading spans more genome. The gap therefore depends on the intron length and on where the other mate lies.
- **De-duplication.** Within a window, STAR drops an alignment whose aligned bases are a subset of another alignment's.
  - A soft-clipped contiguous alternative (`97M3S`) is a subset of the spliced one (`97M300N3M`). So local alignment can
    delete the contiguous reading before selection ever sees it.
  - `--alignEndsType EndToEnd` would keep it, as `100M` with explicit mismatches.
- **Generation.** Some alternatives may never be built, because of seed or window limits. For those, no score range
  helps.
- **The cap.** A wider range raises NH, and a template with more than `--outFilterMultimapNmax` hits (production 50) is
  reported unmapped.

## 4. The cluster

Verify each item, and mark what you could not.

- **Scheduler.** SLURM, partition `standard`. Check the account with `my_accounts` or
  `sacctmgr show assoc user=$USER format=account%20,partition%20`. Run STAR and any heavy Python only through `sbatch`
  on compute nodes. The login node is for editing, small inspections and submitting. When waiting, check `squeue` and
  `sacct` at intervals matched to the job's expected length, never in a tight loop.
- **STAR 2.7.11b.** The VCaP alignment (July 2026) used
  `<cluster scratch>/shared_data/hulkrna/.snakemake/conda/56d06cf1324bdacee1a20d4400bbf302_/bin/STAR-avx2`.
  A second candidate is the env `alignable` under `<cluster home>/sw/miniforge3` (unverified). If neither works, create
  a NEW env, e.g. `star_params` with bioconda `star=2.7.11b samtools pysam numpy pandas`. Never modify an existing env.
  Confirm that `STAR --version` prints `2.7.11b`.
- **STAR's source.** `git clone --depth 1 --branch 2.7.11b https://github.com/alexdobin/STAR.git`. Its
  `source/parametersDefault` documents every parameter.
- **Where to work.**
  - Scripts, metrics and the report go in `<cluster home>/proj/rigel_star_params/`, where this brief is. It is small and
    is not purged.
  - BAMs and temporary files go in `<cluster scratch>/shared_data/rigel_star_params/work/`.
  - The inputs the owner pushed are in `<cluster scratch>/shared_data/rigel_star_params/inputs/`.

  Scratch purges files that go unused for a long time, so do not leave anything only there that you need later.
- **The production parameters**, as the VCaP BAM header records them. `--quantMode GeneCounts` and the read-group line
  are left out because they do not change alignments.

  ```
  --outSAMtype BAM Unsorted --outSAMunmapped Within KeepPairs
  --outSAMstrandField intronMotif
  --outSAMattributes NH HI AS NM MD MC
  --outFilterMultimapNmax 50
  --alignSJstitchMismatchNmax 5 -1 5 5
  --alignSplicedMateMapLminOverLmate 0.5
  --peOverlapNbasesMin 10
  --chimSegmentMin 10 --chimScoreDropMax 30 --chimScoreSeparation 1 --chimScoreJunctionNonGTAG 0
  --chimJunctionOverhangMin 10 --chimOutType WithinBAM HardClip --chimSegmentReadGapMax 3 --chimMultimapNmax 50
  ```

  Add `nM jM jI ch` to `--outSAMattributes` in your runs. `jM` marks an annotated junction with +20, and `jI` gives the
  intron coordinates. On a sample, show that the added attributes change no alignment.

## 5. Data and truth

**The simulated ladder** gives exact truth, and you rank configurations on it. The owner pushes two conditions:
`gdna_g50_ss_0.99_nrna_mid_capture_off` and `gdna_g50_ss_0.99_nrna_mid_capture_on`.

- **Files.** `inputs/ladder/<condition>/{sim_R1.fq.gz, sim_R2.fq.gz, sim_oracle.bam, truth_summary.json}` and
  `inputs/reference/{genome.fa, genome.fa.fai, genes.gtf}`.
- **Genome.** chr21, chr22 and 92 ERCC contigs, 97.6 Mbp, with the ladder's own annotation.
- **Size.** 10 M templates per condition, 2×100 bp. Fragments are 206 ± 98 bp, within [50, 500]. Capture OFF holds
  3.99 M mature RNA, 1.01 M nascent (unspliced) RNA and 5.00 M gDNA.
- **The reads are error-free** (simulator error rate 0.0). The sequencing-error mechanism the owner names is therefore
  absent; see the injection instrument in §6.
- **Truth is in the read name.**
  - `{transcript_id}:{start}-{end}:{strand}:{index}` is mature RNA.
  - `nrna_{transcript_id}:...` is nascent RNA.
  - `gdna:{ref}:{start}-{end}:{strand}:{index}` is gDNA.

  `sim_oracle.bam` holds every read's true alignment and is the truth for "which reading is right". Check it before
  relying on it: on error-free reads, every oracle record should match the reference exactly.
- **All-N reads.** About a fifth of the gDNA at capture OFF is all-N reads, because the simulator draws gDNA over the
  genome's N runs. Exclude them from every rate and report their count.
- **Remote copies.** The ladder genome holds two chromosomes, so remote copies are rare in it. Take the remote-copy class
  from VCaP, and say where the two inputs disagree.
- **The index.** The STAR index for the ladder is cheap to build: `--runMode genomeGenerate --genomeFastaFiles genome.fa
  --sjdbGTFfile genes.gtf --sjdbOverhang 100 --genomeSAindexNbases 12`.

**The VCaP mix** is real data, with truth only by library. Use it to test the ranking, never to make it.

- **Sources.** It is an in-silico blend of two real VCaP libraries. The instrument field of the read name names the
  source:
  - The DNA library's instrument prefix marks the exome DNA library (NovaSeq, 2×150). Every spliced alignment of it is an artifact.
  - The RNA library's instrument prefix marks the transcriptome library (HiSeq, 2×125, stranded). It holds some gDNA of its own, so report
    the RNA-side costs as upper bounds, with the strand-implied split (its wrong-strand spliced templates are its gDNA).
- **Inputs recorded in the BAM header**, relative to `<cluster scratch>/shared_data/hulkrna/`:
  - FASTQs: `runs/human/vcap_mix/fastq/k2_prefilter_{1,2}.dedup.fastq.gz`.
  - Index: `refs/human/star_index`. Read its `genomeParameters.txt` for the FASTA, the GTF and `sjdbOverhang`.
  - The original output: `runs/human/vcap_mix/star/`.

  These may have been purged. If they have, ask the owner: the laptop holds the aligned BAM, and
  `samtools fastq -F 0x900` regenerates the reads from it. That BAM has had duplicates removed.
- **A free check.** If the original `Aligned.out.bam` survives, your production run on the same reads must reproduce
  its primary alignments. That shows the setup is right.

## 6. The work, in order

**Phase 0: setup (login node, minutes).** Verify everything in §4–5. Write `RUNLOG.md`, which records every job: its
command, inputs, estimated and actual core-hours, MaxRSS and elapsed time (`sacct`). Then tell the owner what is present
and what is missing.

**Phase 1: the mechanism, from the source (no compute).** Trace one paired-end template through STAR 2.7.11b, citing
file:line:

- how alternative alignments of one template are BUILT: seeds, windows, stitching and extension;
- how they are de-duplicated within a window;
- how they are SCORED: matches, mismatches, indels, junction-motif penalties, sjdbScore, the genomic-length penalty and
  its rounding, clipping;
- how they are SELECTED: the score range, the filters, the multimapper cap, the primary;
- what is written for each one: NH, HI, AS, nM, NM, MD, jM, jI, and whether secondaries carry SEQ.

Likely files are `ReadAlign_mapOneRead.cpp`, `ReadAlign_stitchPieces.cpp`, `stitchWindowAligns.cpp`,
`stitchAlignToTranscript.cpp`, `extendAlign.cpp`, `ReadAlign_multMapSelect.cpp`, `ReadAlign_mappedFilter.cpp`,
`ReadAlign_outputAlignments.cpp` and `ReadAlign_peOverlapMergeMap.cpp`. These names are from memory, so correct them.

Then, for each situation in §2 and each hypothesis in §3, DERIVE three things under the production parameters before
you measure anything:

- the score of the spliced reading and of the contiguous reading;
- which step removes the alternative;
- which parameter would move it.

Write the predictions down. Phase 3 checks them.

**Phase 2: instruments.** Each one gets a self-test, and then a deliberate break that the self-test must catch.

- **Truth matcher.** A hit is the TRUE READING if every aligned base of both mates sits where the oracle puts it.
  Compare query position → reference position in the read's original orientation, e.g. with pysam
  `get_aligned_pairs(matches_only=True)`. Soft-clipped bases are allowed and counted. Otherwise classify the hit: same
  locus (it overlaps the oracle's span) or remote, and spliced or contiguous.
  - Self-test: the oracle matched against itself is 100 % true.
  - Break: shift one mate by 1 bp, and remove one intron from a CIGAR. Each must be rejected.
- **STAR's score, recomputed.** Implement the Phase 1 score.
  - Prove that it reproduces AS exactly on a large sample of reported hits. Explain every disagreement; overlapping
    mates merged by `--peOverlapNbasesMin` are a likely one.
  - Then score the oracle reading of every template STAR got wrong. The gap between STAR's best and the truth says which
    range would have reported it, if STAR built it. Comparing that with what STAR actually reports separates the
    failures of selection from the failures of generation.
  - Break: set the sjdb bonus to a wrong value, and agreement must fail on annotated-spliced hits only.
- **Error injection.** Write seeded copies of the ladder FASTQs, with uniform substitutions at two rates:
  - the per-base mismatch rate STAR reports for the VCaP DNA library in `Log.final.out` (an upper bound, since VCaP's
    own variants are in it);
  - 0.01, the owner's setting for the alignable store, as a stress.

  The names are kept, so the truth is kept. Self-test: injected substitutions per read equal the oracle reading's
  recomputed mismatches.
- **Census.** One row per template:
  - the truth class: gdna, nrna, mature-spliced or mature-unspliced, from the name and the oracle;
  - STAR's outcome: true reading primary; true reading secondary, with its score gap; or absent;
  - for the primary: spliced annotated, spliced unannotated or contiguous, and local or remote;
  - unmapped, with the reason, and chimeric;
  - NH;
  - for spliced hits: the short side's length, which side it is (donor or acceptor), and h, which you compute from the
    reference.

  VCaP gets the same census by library.

**Phase 3: the baseline, with the production parameters.**

- **Inputs.** A fixed, seeded subsample per ladder condition, the same template list in every arm, at three error levels
  (0, the VCaP rate and 0.01). Size it so that no class you must read is thin, and state each class's count beside its
  rate. On VCaP, the whole DNA library and a fixed RNA subsample. Run VCaP first, because it gives the injection rate.
- **Report** the census, and compare the observed gaps with the Phase 1 predictions.
- **Then stop and write `CHECKPOINT.md` for the owner.** It holds the mechanism, the baseline, where the predictions held
  or failed, and the sweep you propose. Continue only when the owner says so.

**Phase 4: one parameter at a time, starting from production.**

- Use the same inputs, re-record the baseline in the same batch, and score every arm against truth, never against
  another arm.
- One wide-range run screens every smaller range offline: keep the hits with AS ≥ best − r. Confirm on one direct run
  that this equals running at r, because NH, the cap and the primary all change.
- Candidates:
  - `--outFilterMultimapScoreRange`, from 1 to the end of Rigel's useful range. Derive that end from the 4-mismatch
    pruning bound, the 2 points per mismatch and the sjdb bonus.
  - `--sjdbScore` at 2, 1 and 0. At 0, a tie between equal bases is an exact tie.
  - `--alignEndsType`: Local, EndToEnd and Extend5pOfReads12.
  - `--scoreGenomicLengthLog2scale` at −0.25 and 0.
  - `--alignSJDBoverhangMin` (3 down to 1, and up) and `--alignSJoverhangMin` (5 and up): how much each situation hangs
    on short sides.
  - `--alignSJstitchMismatchNmax`: production `5 -1 5 5` against the default `0 -1 0 0`.
  - `--outFilterType BySJout`.
  - `--outFilterMultimapNmax`, varied only together with the range, to show what the cap removes.
  - Search breadth (`--winAnchorMultimapNmax`, `--seedSearchStartLmax`, `--alignTranscriptsPerWindowNmax`), only if
    Phase 3 finds alternatives that are never built.

**Phase 5: combinations, and the configurations users run.**

- Combine only parameters that helped on their own.
- Also run, as context for what users will hand Rigel: STAR's defaults; the ENCODE long-RNA pipeline's parameters;
  nf-core/rnaseq's; and `--twopassMode Basic`. Take each one from its source and cite it, never from memory.

**Phase 6: the double alignment.** If no single run reports the alternatives, price the fallback.

- An unspliced STAR run on the same reads: an index without annotation, `--alignIntronMax 1`, and `--alignMatesGapMax`
  taken from the fragment-length range.
- Merge it with the spliced run by read name: take the union of hits, de-duplicate identical pair alignments, recompute
  NH and HI, and choose one primary.
- Run the same census, and add the merge's cost and every assumption it makes about tags.

**Phase 7: confirm and hand over.**

- Run production and the recommended configuration(s) on both full ladder conditions and the full VCaP mix.
- Keep those BAMs, name-collated as Rigel reads them, for the later end-to-end Rigel A/B. That A/B is not this session's
  job: Rigel's pairing is under repair, and its numbers would confound the two.

## 7. Metrics

Give every number per truth class and per condition, never pooled across capture OFF and ON, with counts beside the
rates.

- **Contiguous truth (gdna, nrna):**
  - spliced-primary rate per million templates, annotated and unannotated, local and remote;
  - **recall**, the fraction of those templates whose true contiguous reading is among the hits;
  - the share lost to unmapped, the cap, or chimeric.
- **Spliced truth:**
  - true reading primary; true reading among the hits; true reading absent;
  - spurious contiguous alternatives, with their mismatches relative to the true reading (only those within 4 matter to
    Rigel).
- **Unspliced mature truth:** spurious spliced alternatives.
- **Mate disagreement:** one mate spliced at a junction while the other runs into that intron or is clipped there.
- **Clip-only alternatives:** alternatives that differ only by soft clipping.
- **Cost:**
  - the NH distribution;
  - the unique-template fraction, by class;
  - templates lost to the cap;
  - records per template and BAM size;
  - run time and peak memory.
- **VCaP:** the same metrics by library, with the RNA side given as an upper bound.

## 8. Rules

- **Budget.** Estimate the core-hours before each batch and record the actuals. Read the first job's MaxRSS and elapsed
  time before launching more. Ask the owner before any batch over 100 core-hours. Sweeps use subsamples; full inputs are
  for Phase 7 only.
- **Varying.** Vary one thing per arm.
- **Values.** Every parameter value you try is derived, or taken from a cited configuration, or an endpoint of a stated
  bracket. No number without its provenance.
- **Instruments.** No instrument is trusted before its self-test passes and its deliberate break fires.
- **Hands off.** Do not modify Rigel or any repository, and make no commits. Do not delete or overwrite the owner's
  existing files on scratch: the hulkrna tree, the aligned-ladder tree, the alignable store.
- **Reporting.** Report faithfully. A hypothesis from this brief that the source or the data refutes is a finding: say so
  plainly.
- **Asking.** Ask the owner when a choice is theirs: the account, the budget, a missing input.

## 9. Deliverables

Everything goes in `<cluster home>/proj/rigel_star_params/`:

- **`REPORT.md`.** It must stand alone for the laptop session and for an external reviewer. Write it for a reader who
  knows RNA-seq but not STAR's internals, and define each term at first use. It holds:
  - the mechanism, with source citations and one template traced through it;
  - the baseline;
  - each parameter's effect, with its cost;
  - the recommendation: configurations ranked, and what each one still misses and why;
  - the verdict on the double alignment;
  - what Rigel's input contract would need: tags, sorting, NH/HI semantics, clipping, the cap;
  - what Rigel must still do itself for BAMs from default STAR.
- **`CHECKPOINT.md`** (Phase 3) and **`RUNLOG.md`**.
- **`scripts/`**: the SLURM jobs, matcher, score, injection and census, each with its self-test.
- **`metrics/`**: one TSV or JSON per arm.
- **`bundle.tar.gz`**: all of the above, without BAMs, for the owner to carry back to the laptop.
