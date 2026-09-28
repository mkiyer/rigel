# Spliced-fragment artifacts — findings and plan, for a session on the cluster (owner, 2026-09-27)

The owner's highest-priority issue (`ISSUES: splicing-artifacts`), moved to the cluster: the data lives there,
`alignable` can be re-run there, and as many libraries as needed can be re-aligned and re-quantified. The
defect: an aligner writes a junction on a molecule that was not spliced there, and Rigel certifies the fragment
as RNA and trains its strand model and its RNA fragment-length law on it. The first session (2026-09-26, on the
owner's Mac) found that the defect is two-sided and that its largest artifact class is decided by the reference
sequence. This file holds those findings, the plan, and what the tree does today.

## The prompt

> Read `CLAUDE.md`, then this file, then the owner's four investigation notes in `docs/dev/`:
> `GDNA_SPLICE_ARTIFACTS_OVERVIEW.md` first, then `_RIGEL.md`, `_ALIGNABLE.md` and `_STAR.md`. They were written
> on 2026-08-03 against rigel 0.7.1: every file:line in them is a hypothesis, and "Corrections to the notes"
> below overrides them. Then read `ISSUES: splicing-artifacts` and `DESIGN.md` §3.1d.
>
> The task is detecting and handling spliced-fragment artifacts in Rigel, robustly to aligner settings, judged
> on BOTH errors together: artifacts certified as RNA, and genuine junction reads rejected. Work the plan below
> one phase at a time and stop to discuss at the end of each. Every mechanism here is a hypothesis until it is
> measured on the libraries. No magic numbers, no yes/no gates on continuous quantities; the owner drives
> commits.

## What the first session found (2026-09-26, from 0.7.1's production outputs)

All numbers below come from the 0.7.1 production run's own `annotated.bam` (production blacklist), primary
records. Phase 1 re-derives them on today's tree.

### The judge: the VCaP mix carries per-fragment truth

`mctp_vcap_rna20m_dna05m` blends two real libraries in silico, and the instrument field of the read name names
the source (owner-confirmed):

- `A00839:75:H7MFFDSXY` (NovaSeq, 150-bp reads): VCaP **exome DNA**. 4.71 M fragments, 10,164 with a sequenced
  junction, every one an alignment artifact.
- `HWI-D00127:189:C6EL5ANXX` (HiSeq, 125-bp reads): VCaP **transcriptome**, mostly RNA (it may hold a little
  gDNA). 14.0 M fragments, 9.31 M spliced.

Exome DNA is concentrated on exons, as the gDNA in a capture-enriched library is. One BAM therefore measures
artifacts caught and genuine reads lost against truth, not against a strand proxy. More such mixes can be built
(Phase 0).

### The defect is two-sided

| 0.7.1, production blacklist | DNA library (all artifacts) | RNA library |
|---|---|---|
| fragments certified spliced (`spliced_annot` + `spliced_unannot`) | 5,034 (escaped) | 9,072,988 |
| fragments tagged `splice_artifact` | 2,660 (caught) | 230,364 |
| uniquely mapped single-junction records rejected | 3,143 of 6,032 | 429,208 of 12.55 M (3.4 %) |
| … with both anchors ≥ 21 | 739 | 235,406 |
| multi-junction records rejected | 632 of 887 | 154,604 of 1.88 M (8.2 %) |

The RNA library's rejections are not its own gDNA. At the DNA library's rejection rate (0.073 % of records), even
a pure-gDNA library of that size would give about 20,000 rejected records, not 598,000. A rejected genuine read
becomes an unspliced, gDNA-eligible fragment whose footprint spans the intron (`DESIGN.md` §3.1d). The cost
therefore lands on the gDNA and synthetic-span pools of exactly the genes whose junctions are listed. The notes
measured only the first error: their 0.19–0.59 % is the site-pair rule's added loss, counted in junctions.

### The largest artifact class is decided by the reference sequence

Take a single-junction record with blocks `[s−a, s) | N | [e, e+b)`. Its unspliced alternative puts the short
block contiguous with the long one: the left block at `[e−a, e)`, or the right block at `[s, s+b)`. Let **h** be
the Hamming distance between the reference at the aligner's placement and at the alternative. At h = 0 the read's
bases fit both readings identically, whatever they are, so the aligner's choice is no evidence of splicing. At an
annotated junction STAR's `sjdbScore` bonus (+2) breaks the tie toward the spliced reading.

| DNA library (all artifacts), single-junction records | n | h = 0 |
|---|---|---|
| annotated junction, short side 1–10 bp | 5,297 | **5,092 (96 %)** |
| annotated, short side 11–20 | 197 | 22 |
| annotated, short side ≥ 21 | 1,538 | 46 (3 %) |
| unannotated, short side 1–10 | 282 | 172 |
| unannotated, short side 11–20 | 494 | 106 |
| unannotated, short side ≥ 21 | 3,302 | 147 (4 %) |

- **The short side is on the acceptor side** in 81 % of the annotated 1–10 bp records (4,272 of 5,297). An exon
  usually ends in `…AG`, like every intron, so an intron's last bases often equal the upstream exon's last bases.
  The splice-site consensus predicts this.
- **The blacklist catches this class only where the simulation sampled the junction.** It rejected 4,643 of
  these 5,092 h = 0 records; the 449 that escaped are the same physics. The reference gives h for every junction,
  with no sampling, no count threshold and no anchor ceiling.
- **RNA library, annotated junctions, short side ≤ 20 bp:** 202,770 records were rejected.
  - 40,557 of them have h = 0. They are ambiguous: genuine, or an unspliced molecule. Only abundance can split
    them.
  - The other 162,213 have h ≥ 1, so reading them as unspliced needs extra mismatches; 141,787 need two or more.
    Most are genuine reads lost.
  - Separately, 30,119 h = 0 records were kept and certified as spliced.

Method and tables: `~/Downloads/rigel_runs/prototypes/2026-09-26_splice_artifacts/` on the owner's Mac. It holds
`extract_single_junction.sh` and `overhang_ambiguity.py`, about 60 lines together, and re-derivable from the
definition above.

### Long-anchor artifacts and long-anchor rejections: remote origin (a hypothesis)

- **alignable records an artifact from every alignment** of a simulated fragment: primary, secondary and
  supplementary (`_scorer_fast.cpp`, `process_group`). It does not use only the alignment STAR ranks best.
- **Consequence.** Take a simulated read from a retrocopy or paralog whose best alignment is unspliced at its
  origin, with a secondary spliced hit on the parent gene. It lists the parent's junction with long anchors on
  both sides. Under the OR rule, an entry whose two stored anchors both reach m rejects every single-junction read
  there up to about 2m aligned bases.
- **This is consistent with the RNA library's data.** Its rejections at ≥ 21 bp are 98 % uniquely mapped
  (235,406 of 240,520).
- **Derived from STAR's default scoring, not yet observed on reads:** a read pair from a retrocopy identical to
  the parent over its span is reported as a unique spliced alignment on the parent when the intron is short
  (under ~5 kb at 300 aligned bases). The reason is that `sjdbScore` (+2) outweighs the genomic-length term
  (−0.25·log2 of the span). The origin then falls outside `outFilterMultimapScoreRange` (1). One mismatch makes it
  a multimapper instead.
- In the DNA library, the artifacts with a short side ≥ 21 bp total 4,840 records, 60 % of them uniquely mapped.

### Corrections to the notes

- **`_ALIGNABLE.md` §6, "defect 1 needs no re-simulation":** this holds only if the per-record data survives.
  - The table inside `alignable.zarr.zip` already has the half-blind aggregation applied, because
    `_aggregate_splice_artifacts` runs before `write_splice_blacklist`.
  - Each simulated record's own anchors and origin live in `<alignable output>/splice/chunk_*.feather`.
    `package_store` does not pack that directory; `scripts/reaggregate_splice.py` reads it in place.
  - The owner does not know whether it survives.
- **0.7.1 `summary.json` counts are the fragment-length models' `n_observations`,** not a census. For example,
  VCaP's `splice_artifact` is 63,187, against 233,024 tagged fragments.
  - The notes' ratios (0.491; 0.010–0.046) are in that unit.
  - Today's tree counts a census of the fragments offered to the accumulator (unique, resolved, non-chimeric).
    Re-derive before comparing.
- `_RIGEL.md` §1.1 and §6 (anchor definitions) and §5a (multimappers): see "What the tree does today".

## What the tree does today

Re-checked from the code at `062bc9ea`; since then only a CLI help string has changed. Line numbers drift, so
treat them as pointers.

- **Anchors are cumulative** (`bam_scanner.cpp` `parse_cigar`).
  - The left anchor is every reference-advancing base (M/D/=/X) before the junction; the right anchor is the
    rest. A later junction's anchors therefore span the earlier introns.
  - `alignable` counts anchors the same way (its `_scorer_fast.cpp`), so the two sides agree on the definition.
  - Blocks are cut at every `N` before any blacklist test.
- **The reject rule** (`filter_blacklisted_sjs`):
  - An exact `(ref, start, end)` match that ignores strand, then `left <= max_left || right <= max_right`. A stored
    0 is inert.
  - It runs per alignment record, before mates are joined, so a junction rejected on one mate and kept on the
    other survives.
- **The blacklist build** (`splice_blacklist.py`):
  - The minimum count (default 2) is applied per `(chrom, intron, strand, read length)` BEFORE aggregating over
    read lengths. A junction seen once at each of several read lengths therefore never enters.
  - The shipped feather holds only `ref, start, end, max_anchor_left, max_anchor_right`: no count, no site
    structure.
  - The index applies `splice_blacklist.feather` at quant time only when its `manifest.json` records the store
    (`sources.alignable_zarr`, written by `rigel index --alignable-zarr`; `index.py`). A feather dropped into an
    index built without a store is ignored (`sj_blacklist_loaded: false` in `summary.json`); in one built with a
    store, replacing the feather takes effect without a rebuild and removing it turns detection off.
- **Classification** (`resolve_context.h`): a surviving junction wins over ARTIFACT, which wins over IMPLICIT.
  IMPLICIT means more than one gap hypothesis survived.
- **An artifact fragment** (every junction rejected):
  - Calibration counts it in the census, then holds it out: it deposits nothing.
  - The EM lets gDNA explain it (`DESIGN.md` §3.1d), scored at its footprint, which still spans the rejected
    intron.
- **A mixed fragment** (one junction rejected, one surviving) is spliced and deposits. The rejected `N` becomes a
  gap:
  - re-implied as a junction where a candidate transcript has an annotated intron there;
  - otherwise kept inside the unspliced hypothesis's length;
  - or deferred to the second pass.

  No test covers this.
- **Only uniquely mapped fragments train the strand model and the RNA fragment-length law.**
  - The strand model trains on SPLICED_ANNOT at the leftmost annotated junction (`strand_model.py`).
  - The RNA law is the `RNA_SPLICED` pool (`fl.py`).
  - Multimappers never entered either, so the contamination is uniquely mapped artifacts at ANNOTATED junctions,
    which is the class the reference decides.
- **Multimappers** (`include_multimap` on) feed only the EM. gDNA may explain a multimapper if any of its hits is
  unspliced (`DESIGN.md` §3.1d), but the scanner drops a hit that overlaps no transcript
  (`ISSUES: multimapper-intergenic-alignments`).
- **Evidence Rigel reads:** NM only, as a whole-fragment penalty of `log(mismatch_alpha)` per mismatch
  (`mismatch_alpha` = 0.1, `scoring.py`, an underived constant). MD, sequence and qualities are never read; the
  scanner never reads the reference.
- **Counters:**
  - `sj_blacklisted` counts junction × record over every record, multimappers included.
  - The `splice.*` census counts uniquely mapped fragments.
- **Tests:** every blacklist fixture uses anchors of 500–10,000, longer than any read. Nothing pins the `<=`, the
  OR, the inert 0, cumulative anchors or the mixed fragment.

## Data

**On the cluster,** under `/scratch/mkiyer_root/mkiyer0/shared_data/`. `/scratch` is purged by policy.

- **The index case**, `mctp_LBX0077_SI_44263_HKWNYDRX7` (96.4 % gDNA, 138-bp mates), under
  `hulkrna_gdna_debug/`:
  - `runs/human/…/bam/star.srt.rmdup.collate.bam`;
  - `evidence/`: `spliced_only.bam`, the per-junction parquets, the four bug tables, `validation_targets.json`
    and the analysis scripts.
  - **GONE (owner, 2026-09-28):** this splice evidence did not survive `/scratch`. It must be regenerated from
    scratch by a new run, a separate task; nothing here can be copied off.
- **The production Rigel index,** `hulkrna/refs/human/rigel_index/`: its `splice_blacklist.feather` (5,833,092
  rows) and `manifest.json`. Before the next quant there, check the manifest's `sources.alignable_zarr` is
  present and not null: otherwise the current tree ignores the feather and detection is off.
- **alignable:**
  - `hulkrna/refs/human/alignable.zarr.zip` (26 GB). Inside, the per-read-length table has 34.9 M rows with
    `count`.
  - The unpackaged output with its `splice/` directory, if it survives.
- **The 344-library cohort's** `rigel/summary.json` and `star/SJ.out.tab.gz`.

**Only on the owner's Mac** (`~/Downloads/rigel_runs/cfrna/`), possibly the last copies — check the cluster
first:

- **The VCaP mix:** its BAM (4.5 GB), the 0.7.1 `rigel/annotated.bam` (5.1 GB) and `summary.json`.
- **LBX0190, LBX0588 and MO_3021** (cfRNA at 13 %, 97 % and 22 % gDNA): each BAM, 0.7.1 `annotated.bam`,
  `summary.json` and `SJ.out.tab.gz`.

## The plan

One phase at a time, each ending in a discussion with the owner. Measurements before mechanisms, one mechanism
per A/B.

### Phase 0 — the substrate on the cluster

1. **Regenerate the index case's evidence** with a new run (a separate task): the `/scratch` copy is gone.
2. **Build the working index** with today's `rigel index`, from the production GTF, FASTA and
   `alignable.zarr.zip`, after the index format bump lands (`ISSUES: the-format-changes-to-batch-before-release`
   (b)), so the cluster rebuilds once.
   - Confirm the checksums against the production index's `manifest.json`.
   - Keep a blacklist-free twin as the "no blacklist" arm.
3. **alignable:** find out whether `splice/` survives.
   - If it does not, plan a re-run that keeps per-record rows: anchors, origin, primary or secondary, NH, and the
     score of the best unspliced alignment. Run splice-only if that is possible (`_ALIGNABLE.md` §6).
   - Change none of its settings until Phase 2 says what the catalogue is for.
4. **Truth-labelled libraries.**
   - The VCaP mix (upload it if the cluster copy is gone).
   - New mixes of matched exome DNA and transcriptome RNA from the same patients or cell lines, spanning DNA
     fraction, read length and chemistry. Each read name carries its source.
5. **Real libraries re-run with `--notemp`.**
   - LBX0077; MI_1104 (98.6 % gDNA); the clean MI_1089 and MI_1083; LBX0054 (74 % gDNA, 104-bp mates, does not
     leak).
   - A spread from the cohort across gDNA fraction and read length, since the notes find leakage grows with read
     length.

### Phase 1 — today's baseline, both errors

- **Per library, on today's tree:**
  - strand specificity;
  - the census (spliced annotated, unannotated, implicit, artifact);
  - blacklisted junctions per fragment;
  - the artifact-signature fraction among fragments called RNA.
- **On the truth-labelled libraries:**
  - artifacts escaped and genuine reads lost, per fragment and per record;
  - broken down by short-side length, h, side, annotation, NH, and single- versus multi-junction;
  - where the lost genuine fragments land in the EM: gDNA, synthetic spans or transcripts.
- **Arms:** no blacklist, and the production blacklist.

### Phase 2 — the mechanism census, on truth

Split every artifact (a DNA-library spliced fragment) and every genuine loss (an RNA-library rejection) by
mechanism before designing anything:

- **Overhang ambiguity:** h = 0 or small at the short side, annotated and unannotated.
- **Remote origin** (retrocopy, paralog, repeat): the read's sequence occurs unspliced elsewhere. To measure it:
  - align each such read without splicing, end to end, against the genome;
  - use alignable's per-record origins where they exist;
  - build a junction mappability table (each annotated junction's spliced sequence searched unspliced in the
    genome — deterministic, no sampling).
- **Mismatch-tolerant stitches** (`alignSJstitchMismatchNmax 5 -1 5 5`): h ≥ 1 at a short side.
- **Genome structure:** deletions in the sample relative to the reference, such as VCaP's rearrangements. They
  show as novel junctions with long anchors on both sides, supported by DNA reads, with no remote origin.
- **Multi-junction records:**
  - which of their junctions is rejected;
  - cumulative against local anchors;
  - the precedence rule, under which one surviving junction shields the rejected ones.

The shares rank the mechanisms. The notes' candidates are ranked after the census, not before it.

### Phase 3 — one mechanism at a time: derive, prototype outside `src/`, A/B

The design direction to test is a hypothesis, not a ruling. An ambiguous or listed junction should ADD an
alternative reading for the EM instead of deleting the junction by a yes/no test: the fragment unspliced at the
alternative placement, eligible for gDNA and RNA. Each reading is weighted by how likely the observed alignment is
under it. This is the owner's framing ("the aligner's rate of writing that junction on genomic sequence at that
anchor length") made per fragment, and it treats both errors with one mechanism. The multimapper machinery already
scores a fragment with several readings.

**Overhang ambiguity first,** the largest measured class.

- The alternative's weight against the aligner's reading is the per-base error odds raised to the extra
  mismatches, about `(ε/3)^h`.
- ε comes from the library (NM, or MD). It must not be inherited from `mismatch_alpha`.
- h comes from the reference: at index time for annotated junctions, at scan time for novel ones.
- The alternative's footprint and length follow from the geometry. The owner's "edit the alignment" is its gDNA
  limit at h = 0.
- No constant, no blacklist, no aligner dependence.
- Derive it on a page first, including the multi-junction case: the alternative is defined per terminal block,
  which is local by construction.

**Remote origin second.** The alternative reading sits at the origin locus. Its sources:

- the aligner's own multimapper hits, which needs `ISSUES: multimapper-intergenic-alignments` so the scanner keeps
  the intergenic hit;
- for reads reported as unique, the junction mappability table or alignable's per-record origins.

**Then whatever the census ranks next.**

**A/B on three substrates together:**

- The truth-labelled mixes: both errors, per mechanism.
- The index case:
  - strand specificity, toward its siblings' 0.991–0.999;
  - the artifact-signature fraction.
- A clean control:
  - presumed-genuine junctions kept (annotated, overhang > 20, ≥ 5 unique crossing reads);
  - `mrna_fraction` stable.

### Phase 4 — the two training sets

Once each spliced fragment carries its probability of being genuine, train the strand model and the RNA
fragment-length law on it: by weight, or on the fragments whose alternative reading is negligible. Judge on the
index case's strand specificity, and on its RNA length law against a clean sibling's.

### Phase 5 — robustness to aligner settings, then the catalogue

- **Re-align the truth-labelled mixes** under several STAR settings:
  - `alignSJDBoverhangMin` 3 and 8;
  - `alignSJstitchMismatchNmax` `5 -1 5 5` and `0 -1 0 0`;
  - `outFilterType BySJout`;
  - another aligner if possible.

  The mechanisms' output should not move with the settings.
- **Only then decide what `alignable`'s catalogue is for, and rebuild it.** The candidate changes:
  - read lengths past production;
  - coverage;
  - the half-blind aggregation;
  - records from the best alignment only, or carrying their flags;
  - the count shipped, so any threshold is decided at run time.

### The owner's candidates and the smaller items, placed

| candidate | where it lands |
|---|---|
| 1. site-pair rule | Phase 2 decides whether any class still needs it after overhang ambiguity and remote origin |
| 2. stored 0 as "unknown"; an explicit "beyond the catalogue's range" outcome | moot if Phase 3 replaces the gate; otherwise Phase 5 |
| 3. clean training sets | Phase 4 |
| 4. explicit anchor floor | superseded if Phase 3 holds: a floor is a constant on a quantity the reference decides |
| 5. local anchors and the precedence rule | Phase 2 (multi-junction) and Phase 3 (the alternative is per terminal block) |
| 6. mismatch-near-junction and per-fragment strand evidence | Phase 2, then Phase 3; needs MD, which Rigel does not read |
| edit the alignment | the gDNA limit of the alternative reading |
| per-read-length minimum count; shipping the count | Phase 5 |
| the three fragment types calibration banks as unspliced but the EM treats as RNA only; the length rule applied differently (`ISSUES: splicing-artifacts`) | their own A/B, independent of the above |
| the unannotated-junction rule | the owner's call, informed by Phase 2's novel-junction census |
| the junction-strand tag chosen from the first 1,000 spliced reads | read it per record; small and independent |

## The rules that bite here

- **No magic numbers.** The notes' constants (minimum count 2, anchor floors of 4, 6 or 11,
  `alignSJDBoverhangMin 8`, the 5-bp mismatch window) are each derived or put to the owner, and so is
  `mismatch_alpha` wherever a new mechanism leans on it.
- **Real data is a test input.** Design from the mechanism; judge on the libraries; report the worst case.
- **Judge both errors together.** Never read the truth-labelled halves pooled, and never read one error alone.
- **The ladder cannot see any of this.** Its oracle BAM holds no artifacts; the owner regenerates the aligned
  ladder later.
- **Tests first.** Write the falsification tests for the current rule before changing it: the `<=`, the OR, the
  inert 0, cumulative anchors, the mixed fragment.
- **The owner drives commits.**
