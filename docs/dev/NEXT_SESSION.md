# NEXT SESSION — spliced-fragment artifacts on real data (owner, 2026-09-26)

The owner's highest-priority issue, and a session devoted to it. The in-scope accuracy work on the simulated
ladder is at diminishing returns, so the work turns to real data. The defect is genomic DNA fragments that the
aligner gives a spliced alignment and the splice blacklist lets through. Rigel then counts them as certified RNA.
The problem is multi-pronged: STAR's settings (the `hulkrna` pipeline), the blacklist's builder (`alignable`),
and Rigel. Rigel's defence is the most important, and considerable changes to its behaviour are expected.

## The prompt

> Read `CLAUDE.md`, then this file, then the owner's four investigation notes in `docs/dev/`:
> `GDNA_SPLICE_ARTIFACTS_OVERVIEW.md` first, then `_RIGEL.md`, `_ALIGNABLE.md` and `_STAR.md`. Then read
> `ISSUES: splicing-artifacts` and `DESIGN.md` §3.1d. The task is **spliced-fragment artifact detection and
> handling in Rigel**. Genomic fragments that the aligner splices, and the blacklist lets through, are counted as
> certified RNA. They also train the two laws the EM leans on for every unspliced fragment: the strand model and
> the RNA fragment-length law. Robustness to aligner settings is the design goal.
>
> The notes are the owner's evidence. Their code references are hypotheses: they were written on 2026-08-03
> against rigel 0.7.1 and index format 7, and the tree has changed since. "What the tree does today" below has
> the claims re-checked on 2026-09-26. Verify anything you build on against the current source, and measure it
> against the library before you believe it.
>
> Do, in order, and stop to discuss at the end of each step:
> 1. **The substrate.** The index case's BAM and evidence bundle are on the cluster's `/scratch`, and the local
>    human index has no splice blacklist. Agree with the owner which libraries, which index and which machine the
>    work runs on (see "Data" below).
> 2. **Today's baseline.** Reproduce the index case's readouts on today's tree: strand specificity, the artifact
>    share of spliced fragments, blacklisted junctions per fragment, and the artifact-signature fraction among
>    fragments called RNA. Do the same for the clean controls. The notes' numbers are from 0.7.1.
> 3. **One mechanism at a time.** Derive it on one page, prototype it outside `src/`, and A/B it against what
>    ships. Judge it on an artifact-rich library and a clean control together: the artifacts caught, and the
>    genuine junctions lost (presumed genuine = annotated, overhang > 20, ≥ 5 unique crossing reads). The owner's
>    candidates, in the notes' rank order, are listed below.
>
> The rules that bite here:
> - **No magic numbers.** The notes carry several constants: the blacklist's minimum count of 2, anchor floors of
>   4, 6 or 11, `alignSJDBoverhangMin 8`, a 5-bp mismatch window. Each is derived or put to the owner. The owner
>   dislikes yes/no gates on continuous quantities, and the blacklist is one: an artifact is the aligner's rate of
>   writing that junction on genomic sequence at that anchor length.
> - **Real data is a test input.** Design from the mechanism and judge on the libraries.
> - **The ladder cannot see any of this.** Its oracle BAM holds no artifacts, and the owner regenerates the
>   aligned ladder later.
> - **The owner drives commits.**

## What the tree does today

Re-checked from the code on 2026-09-26, at `062bc9ea`; line numbers drift, so treat them as pointers. The notes
still hold unless stated otherwise here.

- **Anchors are cumulative** (`bam_scanner.cpp` `parse_cigar`, ~456–505). The left anchor is every
  reference-advancing base (M/D/=/X) before the junction, and the right anchor is the rest, so a later junction's
  anchors span the earlier introns. Blocks are cut at every `N` before any blacklist test.
  - **Correction to `_RIGEL.md` §6:** `alignable` computes anchors the same cumulative way (its `scorer.py` and
    its test on `30M500N20M1000N50M`). The two sides agree on the definition; what differs is the reads each
    side measured.
- **The reject rule is unchanged** (`filter_blacklisted_sjs`, ~517–533).
  - It is an exact `(ref, start, end)` match that ignores strand, with the test
    `left <= max_left || right <= max_right`. A stored 0 is inert.
  - It runs per alignment record, before the mates are joined, so a junction rejected on one mate and kept on
    the other survives.
- **The blacklist build** (`splice_blacklist.py`, `aggregate_splice_blacklist`).
  - It drops rows under the minimum count (default 2) per `(chrom, intron, strand, read length)` BEFORE
    aggregating over read lengths, so a junction seen once at each of several read lengths never enters. This is
    a Rigel-side defect of its own.
  - The shipped feather holds only `ref, start, end, max_anchor_left, max_anchor_right`: no count, and no donor or
    acceptor site structure.
- **Classification is unchanged** (`resolve_context.h`, ~1338–1342): any surviving junction wins over ARTIFACT,
  which wins over IMPLICIT. IMPLICIT now means more than one gap hypothesis survived.
- **An artifact fragment** (every junction rejected):
  - calibration: counted in the census, then held out; it deposits nothing and enters no length pool;
  - the EM: gDNA may explain it (`DESIGN.md` §3.1d), scored at its footprint, which still spans the rejected
    intron. That is the length over-statement `ISSUES: splicing-artifacts` wants repaired by editing the
    alignment.
- **A mixed fragment** (one junction rejected, one surviving) is spliced, and it deposits. The rejected `N`
  becomes a gap:
  - re-implied as a junction where a candidate transcript has an annotated intron there;
  - otherwise kept inside the unspliced hypothesis's length, which inflates it or drops the fragment as too long;
  - possibly deferred to the second pass.

  No test covers this.
- **The strand model and the RNA fragment-length law take only uniquely mapped fragments.**
  - The strand model trains on SPLICED_ANNOT fragments at their leftmost annotated junction (`strand_model.py`,
    ~366–414).
  - The RNA law is the `RNA_SPLICED` pool: an annotated junction on the one surviving path, sequenced or implied,
    de-tilted (`fl.py`).
  - **Correction to `_RIGEL.md` §5a:** keeping multimapper-only novel junctions out of these training sets
    changes nothing, because multimappers never entered them. Their contamination is the uniquely mapped
    artifacts at ANNOTATED junctions: the short-anchor chance stitches, which only the anchor length separates.
- **Multimappers** (`include_multimap` on by default) feed only the EM. They get no census entry, no strand
  observation and no deposit, and the EM's gDNA term is the mean over a group's eligible hits.
- **Evidence Rigel reads:** NM, used only as a whole-fragment penalty in the scorer. MD, sequence and qualities are
  never read.
- **`--splicing-anchor-tolerance`** (default 3) is the gap-hypothesis slack on annotated introns, unrelated to the
  blacklist.
- **Counters:** `sj_blacklisted` counts junction × record, over primary, secondary and supplementary records,
  multimappers included. The `splice.*` census counts uniquely mapped fragments.
- **Tests:** every blacklist fixture uses anchors of 500, 1,000 or 10,000, larger than any read. Nothing pins the
  `<=`, the OR, the inert 0, cumulative anchors or the mixed fragment. Write the falsification tests before
  changing the rule.

## Data

- **Index case** (`mctp_LBX0077_SI_44263_HKWNYDRX7`, 96.4 % gDNA, 138-bp mates) and its evidence bundle are on
  the cluster at `/scratch/mkiyer_root/mkiyer0/shared_data/hulkrna_gdna_debug/`. That includes the BAM Rigel
  reads, `spliced_only.bam` (11 MB), the per-junction parquet tables, the four bug tables and
  `validation_targets.json`. **`/scratch` is purged by policy: copy it somewhere durable.**
- **Second artifact-rich library:** `mctp_MI_1104_SI_43016_HNMVYDMXY` (98.6 % gDNA, strand specificity 0.909).
- **Clean controls:** `mctp_MI_1089_SI_40085_HKF7JDRX5` (0.1 % gDNA) and `mctp_MI_1083_SI_40080_HKF7JDRX5`
  (0.3 %). **Intermediate, does not leak:** `mctp_LBX0054_SI_43886_HKWNYDRX7` (74 %, 104-bp mates).
- **The production blacklist** comes from `alignable.zarr.zip` (26 GB, cluster). The local human index
  (`~/Downloads/rigel_runs/refs/rigel_index`) was built without it (`alignable_zarr: null`, no
  `splice_blacklist.feather`), so a local run never exercises the blacklist. Local work needs a production index,
  or at least its `splice_blacklist.feather`.
- **Local real libraries** (`~/Downloads/rigel_runs/cfrna/`): LBX0190, LBX0588 and MO_3021 (cfRNA), and
  `mctp_vcap_rna20m_dna05m`, a real RNA + DNA mix. They are smoke tests, and their artifact load is unmeasured.
- **The other repositories:** `~/proj/alignable` (the catalogue) and `~/proj/hulkrna` (the pipeline and STAR's
  parameter file). The anchor definition must be agreed across Rigel and `alignable` before either changes.

## The owner's candidates for Rigel, in the notes' order

1. A **site-pair rule** beside the exact-junction match: a junction is a candidate when its donor AND its
   acceptor are both known artifact sites. The anchor test stays in force, and the two paths are reported apart.
   Measured on 0.7.1: novel-junction capture 36 % → 65 %, at a cost of 0.2–0.6 % of presumed-genuine junctions.
2. A **stored anchor of 0 read as "unknown"**, not as a threshold, and an explicit "beyond the catalogue's
   range" outcome, counted, instead of a silent pass.
3. **Clean training sets for the strand model and the RNA fragment-length law.** This is the second-order
   damage, and possibly the largest. The notes' multimapper filter misses it: the contamination is the uniquely
   mapped short-anchor stitches at annotated junctions (see "What the tree does today").
4. An **explicit anchor floor**, aligner-independent, with its true-positive cost measured first.
5. **Anchors near their own junction, and the precedence rule.** Both sides count anchors cumulatively. The
   question is whether the test should use the bases local to each junction, and whether one surviving junction
   should shield the rejected ones in the same fragment.
6. **Mismatch-near-junction and per-fragment strand evidence** for the novel multimapper population. This needs
   MD, which Rigel does not read today.

From `ISSUES: splicing-artifacts` itself, and the code check:
- The owner's first repair to consider: **edit the alignment**. Delete the rejected junction and keep the
  largest aligned block, so the fragment is unspliced in both stages at its true length.
- The blacklist build's per-read-length minimum count, and shipping the count so the threshold is decided at
  run time rather than at index time.
- The three fragment types calibration puts on its unspliced bank that the EM treats as RNA only.
- The unannotated-junction rule, which is the owner's call.

A related real-data item, small and clear: the junction-strand tag (`XS` or `ts`) is chosen from the first 1,000
spliced reads (`detect_sj_strand_tag`), so a library whose early spliced reads lack it loses the tag throughout;
it should be read per record.
