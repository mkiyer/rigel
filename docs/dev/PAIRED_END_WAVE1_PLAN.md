# Wave 1 plan: the scan stage

2026-10-17 · base `3941a0dd` · revision 3. Step 1 is complete: its review is delivered and every cleanup stage has
landed (§5). Nothing in steps 2–6 is implemented.

Revision 1 also scrutinized each pair's alignment and constructed alternative readings, using a genome stored in the
index. The owner then narrowed Rigel's scope (2026-10-16): Rigel reads alignments and does not repair or re-align
them. Revision 2 replaced revision 1 and was accepted on 2026-10-16. Revision 3 adds what has happened since: step 1
carried out (§5); the owner's rulings of 2026-10-16 on pairing, mates that disagree and the anchor tolerance (§7.1,
§7.3, §7.5); and those of 2026-10-17: the junction-strand tag read per record (§6), the anchor tolerance renamed
(§7.5) and the output names settled (§13). The waves are in `PAIRED_END_ROADMAP.md`. Part A is the prompt for the
reviewer; Part B is the plan.

# Part A — Prompt for the reviewer (copy everything below the line)

---

You are an independent technical reviewer with a background in short-read alignment (STAR), paired-end library
geometry, splice-junction biology, the likelihood models of RNA-seq quantification, and C++/Python code health. You
are reviewing the **design of one stage** of **Rigel**, a Bayesian RNA-seq quantifier that separates RNA from
genomic-DNA (gDNA) contamination.

Rigel reads a name-collated BAM, calibrates the library's gDNA and RNA content, and runs a per-locus EM over annotated
transcripts and one gDNA component. The design under review is **wave 1, the scan stage**: how alignment records
become complete, validated read-pair hits before any likelihood is computed. It covers:

- a code-health review of the scan path and its cleanup, each stage proven to change no number (carried out: §5);
- each alignment record as a mate location, its junction-strand tag read per record;
- hits from the aligner's pairs (STAR's HI groups) or, for aligners without HI, from a bounded R1 × R2 search, each
  validated as one molecule;
- one rule for unannotated junctions, applied hit by hit;
- the hand-off to calibration and the EM, with four corrections to calibration's second pass;
- the accounting of every template and every record.

**History.**

- A first design (`PAIRED_END_REPAIR.md`, revision 2) was reviewed (`PAIRED_END_REPAIR_INDEPENDENT_REVIEW.md`) and
  had four blockers.
- A first wave-1 plan followed, in which Rigel would also have constructed alternative readings: the contiguous
  reading of a short spliced side, and consistent readings of mates that disagree. It was reviewed
  (`PAIRED_END_WAVE1_INDEPENDENT_REVIEW.md`, seven blockers).
- The owner then narrowed Rigel's scope (2026-10-16):
  - Rigel reads alignments; it does not repair or re-align them.
  - Near-tied spliced and contiguous alignments are the aligner's to report. A parallel study is establishing STAR
    parameters that do so (`star_parameters/PROMPT.md`), and Rigel will recommend them.
  - The genome store is dropped for now.
- Revision 2 of this plan followed and was accepted. Its step 1, the code-health review and cleanup, has since been
  carried out (§5). This is revision 3; it has not been reviewed.

The owner's earlier decisions still hold:

- plausible RNA and gDNA alignments of one pair compete in the EM;
- calibration's two-pass approximation is accepted as an approximation;
- unannotated spliced alignments are withheld, because Rigel does not assemble transcripts;
- a template whose fragment length cannot be computed is withheld;
- the multimapper combination rules belong to wave 2;
- the downstream calibration defects belong to wave 3.

**Your job:**

1. Find the plan's defects, most severe first. Name every blocker: a defect that must be fixed before implementation,
   with the evidence that it is one.
2. Check the arguments the plan relies on:
   - the two pairing strategies (§7.1): HI groups for STAR, and a bounded R1 × R2 search for aligners without HI. Is
     anything lost, and is the search's bound complete?
   - passing mates that disagree on as reported, priced by the existing overhang term (§7.3; the owner agrees);
   - the existing EM combination is right for alternatives the aligner reports at one place (§9.2);
   - the unannotated rule, applied per hit (§8);
   - reading the junction-strand tag per record instead of choosing one mode from the first spliced reads (§6).
3. Check §14: does revision 3 answer the blockers of both earlier reviews, or only claim to?
4. Review the code-health step as carried out (§5): did the cleanup cover what steps 2–6 will change, and was each
   stage's proof sufficient? Three stages were measured rather than proven no-ops: the re-cache, the output names'
   relabelling, and the tile-invariant fits.
5. Answer the open questions (§13), and say what evidence must gate implementation and the end of wave 1.

**Materials.** All in the bundle, at base commit `3941a0dd` (the cleanup's commits are `41899ae4..3941a0dd`; their
diff is in the bundle):

- this plan, which you should read first and completely;
- `PAIRED_END_ROADMAP.md`;
- the STAR study's brief, `star_parameters/PROMPT.md`;
- the first design, revision 2 (`PAIRED_END_REPAIR.md`), with its evidence (`paired_end_repair_evidence/`);
- both earlier reviews and their evidence (`paired_end_repair_review_evidence/`, `paired_end_wave1_review_evidence/`);
- the owner's restart note (`paired_fragment_completion/RESTART_HANDOFF.txt`);
- the code-health review (`SCAN_CODE_HEALTH_REVIEW.md`, with its nine slice reports), the output-name convention
  (`OUTPUT_NAMES_PROPOSAL.md`, ruled 2026-10-17) and the constants rulings (`CONSTANTS_PROPOSAL.md`);
- the code, as it stands after the cleanup:
  - `src/rigel/native/bam_scanner.cpp`;
  - `resolve.cpp` and `resolve_context.h`;
  - `constants.h`;
  - `scoring.cpp` (the multimapper merge is `score_mm_alignment` and `flush_mm_group`);
  - `calibration/accumulator.*`;
  - `second_pass.py`;
  - `pipeline.py`;
  - `config.py`;
  - `locus.py`;
- `docs/ISSUES.md`, for the `hygiene-ledger`, `multimapper-intergenic-alignments`,
  `sj-strand-tag-chosen-from-the-first-reads` and `the-anchor-tolerance-applies-beyond-implicit-splicing` entries;
- the project's rules: `CLAUDE.md` and `docs/DESIGN.md` §0.

The BAMs are not included; the recorded measurements stand in for them.

**Rules the plan must keep** (the owner's). You may argue that one costs accuracy, but only with a measured or
derivable cost.

1. **Three populations:** gDNA, RNA+, RNA−.
2. **Both mates.** Every hit places both mates and has a computable fragment length. No single-mate scoring and no
   mate rescue.
3. **The emission model is kept:**
   - candidates by overlap (one base nominates);
   - an overhang penalty per base;
   - a mismatch penalty per mismatch;
   - per-candidate lengths;
   - no candidate deleted to make a path fit.
4. **The aligner proposes, the EM decides.** The scan stage never chooses between hits.
5. **Read, don't repair.** Rigel never moves bases, never forms an alignment the aligner did not report, and never
   re-aligns. Pairing mate alignments that an aligner reports per mate is reading, not repair.
6. **Annotated transcripts only.** No transcript assembly; chimeras are not quantified.
7. **No new constant.** A tuned cut on a continuous quantity counts as one.
8. **One mechanism per A/B.** Each A/B is pinned (`scan.total_threads=1`, `OMP_NUM_THREADS=1`, fractional
   assignment), scored against truth, and read per stratum (unstranded OFF, stranded OFF, stranded ON; unstranded ON
   reported separately), with both errors: artifacts counted as RNA, and genuine RNA lost.

**Deliver** a markdown report of about 3–6 pages:

1. **Verdict,** in at most 10 lines.
2. **Findings,** numbered. Each has:
   - a severity: blocker, fix before implementing, or consider later;
   - the section it concerns;
   - the reasoning, labelled *derived*, *measured* (cite the number) or *conjecture*;
   - a falsifiable test.
3. **The five arguments** above, each confirmed or refuted with a finite example.
4. **§14:** each earlier blocker, answered or not.
5. **The code-health step (§5):** what the cleanup missed that steps 2–6 depend on.
6. **Calibration (§9.3):** are the four corrections right? Is the all-held approximation acceptable now that the
   recommended STAR parameters will hold more templates, and how should it be measured?
7. **Tests (§11):** recompute at least four expected values; name the missing cases and break tests.
8. **The open questions (§13),** one paragraph each, with a recommendation.

Be specific and terse. A finding without a test is a conjecture: say so. Do not write production code.

---

# Part B — The plan

## 1. What wave 1 must achieve

For every template in the BAM, the scan stage produces either a set of validated hits, or one disposition saying why
there is none. A **hit** is one alignment of the template as the aligner reported it: both mates, checked to be one
molecule, with a computable fragment length.

Wave 1 delivers six things:

1. a code-health review of the scan path, and its cleanup, before anything else changes (carried out through the
   re-cache, §5);
2. one mate location per alignment record;
3. hits from the aligner's pairs, or from R1 × R2 pairing where the aligner leaves pairing to the reader, each
   validated;
4. the unannotated-junction rule, applied hit by hit;
5. the hand-off to calibration and the EM, with four corrections to calibration's second pass;
6. one disposition per template, with every record written once.

Wave 1 does **not**:

- construct, edit or re-align any alignment: no bases are moved, no reading is proposed, and no genome is stored
  (owner, 2026-10-16);
- change how a template's several hits combine (wave 2);
- change the gDNA length estimator or the capture ruler (wave 3);
- infer molecule length from clips (wave 3).

## 2. Background for a reader new to Rigel

**A template** is one sequenced molecule (fragment) and its two reads (mates), which sit at the molecule's two ends
and face inward. The aligner reports each alignment of the pair as two records. STAR groups the records of one pair
alignment by its HI tag, and NH is the number of alignments it reported.

**The EM** scores each hit's candidates:

- every annotated transcript the hit overlaps by at least one base;
- a gDNA candidate, if the hit is unspliced.

Each candidate's likelihood uses:

- the fragment's strand;
- its length along that candidate: along the transcript's path, or along the genome for gDNA;
- its **overhang**: aligned bases outside the transcript's exons, ×0.1 each;
- its mismatches: ×0.1 per NM.

A template with several hits is one EM unit (a multimapper).

**Calibration** estimates the library's gDNA and RNA content before the EM, in two passes. Pass 1 deposits
unambiguous templates (one hit). Pass 2 assigns each held template to one of its placements by a seeded draw,
weighted by pass-1 traffic. The owner treats this as a deliberate approximation to a circular problem: the EM needs
calibration, and calibration would need the EM.

**The bug that started this.** The scanner built "fragments" out of single mates. A lone secondary, a primary whose
mate was unmapped, or a half-filled HI group each became a read-length fragment, a false multimapper and a 100-bp
length deposit. Fixing it exposed two downstream defects that had compensated for it (the roadmap, wave 3).

**The scope (owner, 2026-10-16).** A gDNA molecule aligned as spliced is certified as RNA, because gDNA cannot
splice. Its bases usually fit a contiguous alignment almost as well, and STAR reports only one of the two.

- Rigel will not construct the other alignment itself.
- The remedy is the aligner's. A parallel study is establishing STAR parameters under which STAR reports both
  alignments when they nearly tie, and Rigel will recommend them. Rigel then scores what is reported.
- Wave 1 is what Rigel must get right whichever parameters produced the BAM.

## 3. Principles

1. **Read, don't repair.** A hit is an alignment the aligner reported, checked. Rigel never moves bases, never forms an
   alignment the aligner did not report, and never re-aligns. Pairing mate alignments that an aligner reports per mate
   (§7.1) is reading, not repair.
2. **Both mates, and a length.** A hit is usable only if it places both mates and its fragment length can be computed.
3. **The emission model is unchanged.**
4. **The aligner proposes, the EM decides.** The scan stage never chooses between hits.
5. **One disposition per template; every record written once.**
6. **Clean before changing.** Every cleanup is proven to change no number, and every behavior change is its own
   measured arm.
7. **No new constant.**

## 4. The scan stage in one picture

```
step 1  the code-health review, and its cleanup as proven no-ops, before anything below changes

the BAM records of one template
 │
 ├─ step 2  one mate location per record (CIGAR map, clips, skips, NM, motif); chains withheld
 ├─ step 3  pairs → hits: HI groups (STAR) or a bounded R1 × R2 search (no HI); each validated as one molecule
 │            one reference · opposite strands · inward · some path with 1 ≤ L ≤ max_frag_length
 ├─ step 4  hits with an unannotated junction are withheld, hit by hit
 ├─ step 5  usable hits → calibration and the EM; num_hits = the number of usable hits
 └─ step 6  one disposition per template; every record written once
```

## 5. Step 1 — the code-health review

**Why first.** The bug that started this lived for months in code no test exercised with real pairing, and two
downstream defects compensated for it. Changing that code safely needs it understood first:

- what each function assumes about its input;
- which paths are dead;
- which names mislead.

The owner asked (2026-10-16) that wave 1 include a general code-health review.

**Scope.**

- **In depth, the scan path**, about 14,000 lines from the BAM record to the EM unit and calibration's deposit:
  - native: `bam_scanner.cpp`, `resolve.cpp`, `resolve_context.h`, `scoring.cpp`, `constants.h`,
    `calibration/accumulator.*`;
  - Python: `buffer.py`, `scan.py` (`routing.py` since `fd012918`), `scan_payload.py`, `scan_cache.py`,
    `scored_fragments.py`, `scoring.py`, `splice.py`, `splice_blacklist.py`, `annotate.py`, `stats.py`,
    `second_pass.py`, `locus.py`, `locus_partition.py`, and the scan's part of `pipeline.py`.
- **A lighter sweep of the rest of `src/`**, about 29,000 lines, looking for:
  - dead code;
  - duplicated rules;
  - stale comments;
  - constants outside `config.CONSTANTS`;
  - code whose only purpose is to repair alignments, which the narrowed scope puts up for removal once a measured
    comparison justifies it.
- **What is already recorded.** `ISSUES: hygiene-ledger` holds what earlier reviews left. The review starts from it
  rather than re-finding it.

**What a finding records:**

- its file and line;
- its category: defect, dead code, duplicate, implicit contract, misleading name, stale comment, bare constant, test
  gap, or out of scope;
- the evidence: a test or an input that shows it;
- the proposed change;
- whether that change moves a number.

**Two kinds of change.**

- **A cleanup moves no number.** Cleanups land first, in small stages, on a fresh worktree. Each stage is proven a
  numeric no-op:
  - `rename_identity.py --check` passes against a reference frozen once before the first stage. The reference covers
    the ladder, the test chromosome and `--bam` on VCaP.
  - The full suite passes at the re-derived baseline.
  - `preflight.py --full` runs after any deletion or rename.

  A rename greps every sense of its token first.
- **A defect moves a number.** It gets a falsification test first, then its own arm in steps 2–6, or an ISSUES entry
  if it lies outside wave 1. It is never folded into a cleanup stage.

**Reviewed, not polished.** Steps 2–3 replace the record grouping and pairing: `group_records_by_hit`,
`pair_multimapper_reads`, `collect_placement` and the `num_hits` rule. That code is reviewed for its defects and its
contracts only.

**Baseline.** Re-derived on a clean worktree before the first stage (2,959 at `41899ae4`) and after every stage since:
3,051 at `3941a0dd`. On this disk, untracked review files in `docs/dev` trip the docs-boundary gate, so the baseline is
always derived on a clean checkout.

**Delivered** on 2026-10-16: `SCAN_CODE_HEALTH_REVIEW.md`, with the nine slice reports in `scan_code_health_slices/`.
The owner took all of it the same day:
- the cleanup stages run in the report's order, starting with the scan-cache key;
- the defects inside the scan path join steps 2–6 below, each as its own arm;
- the later-wave findings are ISSUES entries.

**Carried out,** 2026-10-16/17, each stage on a clean worktree, gated by the full suite, `rename_identity.py --check` on
three frozen references (the ladder's oracle, the minimap2 panel and VCaP) and `preflight.py --full`:

| stage | commit | what | the proof |
|---|---|---|---|
| the scan-cache key | `b6aec8a4` | the key covers the splice blacklist and the BAM's bytes | bit-identical |
| the words | `a7c5ebe3` | comments, docstrings and help text, 45 files | bit-identical |
| dead code | `26d37589` | what nothing reaches, 51 files | bit-identical |
| one definition each | `c1156e86` | native codes exported, `auto` resolved once, shared helpers, the read-name format | bit-identical |
| contracts | `1b4a46c6` | refusals checked before anything changes | bit-identical |
| one re-cache | `d5a08e46` | the fields nothing reads, deleted from the scan, the payload and the cache | measured: only the two deleted region banks gone; the transcript table and calibration identical |
| the oracle's provenance | `30357fb3` | an oracle part records the BAM it was split from | bit-identical |
| constants, part 1 | `bb3967cb` | six numbers given their reason or a name | bit-identical |
| the identity check | `8ebf71ea` | `rename_identity.py` digests the transcript table's content, not its column names | its self-test |
| internal names | `fd012918` | names that said the wrong thing, inside the code | bit-identical |
| the output names | `1fa1c143` | outputs report RNA as annotated or nascent (owner's convention); `summary.json` schema 4 | measured: every number unchanged on the ladder end to end, only the defined relabellings move |
| the public repository's scrub | `8401aeba` | no local path, cluster detail or lab identifier in a tracked file; the data root is `$RIGEL_RUNS` | bit-identical |
| index and scan names | `dc736e39` | index format 10; `--nascent-tolerance`, `--intron-overhang-tolerance`; the contiguous reaches no longer computed | bit-identical; every index rebuilt and checked file by file |
| the truth tables | `5c67c1f4` | the simulator's one truth writer, each quantity once | an A/B gives identical reads |
| constants, part 2 | `0c754afc` | the resource and display numbers given a home and a reason | bit-identical |
| tile-invariant fits | `4bcc72d9` | the row-tiled calibration fits give one answer whatever their tile height | measured: noise-level per stratum on the ladder (0.997–1.003 of HEAD) |

The suite went from 2,959 to 3,051 cases: gates added for each refusal and definition, cases deleted with the code
they tested. Every new native gate was broken in a rebuild and seen to fire. Every panel was re-cached twice (262
conditions re-certified each time: once for the payload's schema, once for index format 10), and the references were
re-frozen at each measured stage, each new freeze proven equal to the one before but for the measured change.
Preflight now checks each cache's key, not only its presence, and that check found the oracle defect.

**Moved out of step 1:** the tag detector's 1,000-read cutoff is repaired in step 2 (the junction-strand tag read per
record); and two renames wait for the next re-simulation, because every panel on disk carries the old names and must
stay readable: the synthetic spans' `RIGEL_NRNA_` ids and the truth columns `pre_capture_*` → `molar_*`.

## 6. Step 2 — records into mate locations

Each eligible, non-supplementary record becomes one **mate location**. It holds:

- reference and strand;
- the full CIGAR as a query-to-reference map;
- soft and hard clips at each end, kept apart;
- the raw CIGAR-N skips (before the blacklist) and the kept ones;
- NM, as the record reports it;
- the motif strand (XS or ts);
- provenance: the record, its HI and its flags.

No bases are kept, because nothing in wave 1 re-scores a sequence.

**Rules:**

- **One location per alignment.** Two records are one location when they map every aligned query base to the same
  genome base with the same clips. STAR repeats a mate's record in every pair alignment that shares it, and those
  records collapse to one location.
- **Conflicts are input defects.** Records of one location that disagree on NM are counted, and the location keeps the
  larger NM. Opposing motif strands make the motif ambiguous.
- **No M/=/X base, no location.**
- **Withheld, not absent.** A location that depends on a supplementary chain, or has a hard-clipped 5′ end, is
  withheld: counted, and not paired. A mate whose only locations are withheld is "unresolved", not "missing".
- **Reporting completeness.** Under STAR's semantics, a template whose NH exceeds its number of distinct HI groups
  stops the run, with the remedy `--outSAMmultNmax -1`. Without HI, completeness is recorded as unknown.
- **NM is required.** A record without NM is counted and warned about; its NM still reads 0. STAR's default attributes
  omit NM, so the recommended command asks for it (`nm-absent-reads-as-zero`).
- **Unindexed references keep their own identity.** Two mates on different references outside the index are a trans
  pair, not one reference (`unindexed-references-share-id-minus-one`).
- **Tag types are checked.** NH and HI must be integers. Another type is an input defect, and is counted.
- **The junction-strand tag is read per record** (owner, 2026-10-17: before 0.8.0). Today one mode (XS, ts, both or
  none) is chosen from the first 1,000 spliced reads, so a library whose first spliced reads lack the tag loses its
  strand (`ISSUES: sj-strand-tag-chosen-from-the-first-reads`). Each location reads its own record's tag by the rules
  the scanner applies today; an explicit `--sj-strand-tag` still restricts which tags are read, and a spliced record
  without one is counted. Its own arm: it moves numbers only on libraries whose tags are absent from their first
  spliced reads, or mixed.

**The index sets up the resolver.** It installs the transcript strands when it loads, and K and the maximum fragment
length come from the scan configuration. Today `scan_and_buffer` installs them on the shared resolver, so every unit
test runs without strands (`resolver-strands-installed-only-by-the-scan`).

## 7. Step 3 — hits from the aligner's pairs

### 7.1 Which pairs: two strategies (owner, 2026-10-16)

The strategy follows what the aligner writes. It is chosen per BAM, the way the junction-strand tag (XS or ts) is
detected today. A template that contradicts the BAM's mode is counted as an input defect.

- **With HI (STAR).** STAR groups the records of each pair alignment by HI. A hit starts from the R1 and R2 records of
  one HI group, and Rigel forms no other combination. Crossing HI groups would only duplicate reported hits, or revive
  pairs STAR did not report.
- **Without HI (minimap2).** minimap2 reports secondaries per mate, so pairing is the reader's job. Each R1 location
  is tried against the R2 locations that could be its partner, and each candidate pair goes through the validator
  (§7.2). The search is kept small:
  - only locations on the same reference, with opposite strands, are paired;
  - each mate's locations are sorted by position, so each R1 location walks its partners in order of distance;
  - the walk stops at a bound read off the index and the two locations: `max_frag_length`, plus both mates' observed
    skips, plus the longest total intron length of any annotated transcript on that reference. Beyond it, no path of
    §7.2 can bring L under `max_frag_length`, so the bound is exact, not a tuned cut;
  - the aligner caps the locations per mate (minimap2 keeps at most 5 secondaries by default), so a template has at
    most a few dozen candidate pairs;
  - the search decides by geometry alone. Today's pairing without HI pairs mates whose resolved transcript sets
    intersect, and every location inside one gene shares the gene's synthetic span, so every combination passes
    (`no-hi-pairing-by-transcript-overlap`).

Pairing mate alignments that the aligner reported per mate is reading, not repair.

### 7.2 One molecule

A candidate pair is a hit when two conditions hold.

1. **The mates are on one reference, with opposite strands.** Call the forward mate F and the reverse mate R; either
   may be R1.
2. **They can be one molecule along some path.**
   - The path P is the genome between the mates with every observed raw skip removed. An alternative P also removes
     one annotated transcript's introns where they lie wholly inside a gap that neither mate aligned. An aligned base
     inside a removed skip is off the path (§7.3).
   - Walking F's 5′ clip back along P gives the molecule's start a. Walking R's 5′ clip forward gives its end b.
   - Require a < b, both mates' aligned bases inside [a, b), and 1 ≤ L ≤ `max_frag_length`, where L is the number of
     path bases between a and b.
   - The path decides acceptance only. It never removes a candidate or changes a candidate's length.

**Outward and clipped pairs.**

- An unclipped pair whose 5′ ends point away from each other is rejected for geometry. Such pairs indicate circular RNA
  or a chimera, which Rigel does not quantify.
- Soft clips can hide a junction and make an inward molecule look outward.
  - Example: a 100-bp insert across TA2's junction, with each mate clipped at the junction (F `50S50M` at 2000, R
    `50M50S` at 1150). Walking the clips along TA2's path recovers it.
  - When no annotated path does, the reason is **unresolved geometry**: counted, not claimed impossible. VCaP has
    2,669 such pairs in 18.7 M templates.

### 7.3 Mates that disagree

Sometimes one mate is spliced across a junction while its partner runs a few bases into the same intron, because the
aligner would not splice so short an end. The pair is passed on **as reported**.

- The partner's bases inside the skip lie off the path. They do not fail containment, and the existing overhang term
  prices them.
- Today's resolver already does this. In the toy case (§11), the acceptor-side pair keeps one candidate, TA2, with 148
  exon bases, overhang 2 and length 148.
- With one candidate the penalty cancels. It matters only where the template has other hits.

**The owner agrees (2026-10-16).** One mate is spliced, so the hit is spliced: RNA only, with no gDNA candidate. The
intron overhang is penalized, but no other transcript is expected to compete for the template, so the penalty does not
matter. Revision 1 constructed both consistent readings instead; that is repair, so it is out, and no fallback is
planned.

### 7.4 What a hit carries

- **Contents.** Each hit carries:
  - both locations;
  - its raw and kept skips;
  - NM summed over the two mates.
- **Kind.** Hits of one template whose spans overlap on one reference are **same-place alternatives**: for example,
  the aligner's spliced and contiguous alignments of one molecule. They are counted apart from **multi-place** hits.
- **`num_hits`** is the number of usable hits after step 4. It replaces today's `max(NH, hits)`.

### 7.5 The anchor tolerance (owner, 2026-10-16)

`--intron-overhang-tolerance` (K, default 3; `--splicing-anchor-tolerance` until 2026-10-17) belongs to implicitly
spliced fragments: mates aligned to different exons with a whole intron between them, where up to K aligned bases may
sit inside the intron. The owner keeps it there, "but we must be judicious about deciding which cases apply it".

Today it is applied more widely (`ISSUES: the-anchor-tolerance-applies-beyond-implicit-splicing`):
- to every non-chimeric fragment's gap hypotheses, SPLICED_ANNOT included;
- to `transcript_contradicts_blocks`;
- to cutting an annotated intron whole when up to K sequenced bases lie inside it.

Wave 1 lists every use, keeps the implicit-splicing ones, and justifies or removes each other one in its own measured
arm. The mates that disagree in §7.3 are not a use of K: the overhang term prices those bases, as the owner ruled.

**Renamed** (owner, 2026-10-17): K is a valid tolerance for implicit splicing, and the renames stage gave it a name
that says what it is, `intron_overhang_tolerance` (`--intron-overhang-tolerance`, `dc736e39`).

## 8. Step 4 — unannotated junctions, hit by hit (owner, 2026-10-16)

**Rule.** A hit containing an observed junction that no annotated transcript contains is not usable:

- **Not as RNA:** no annotated transcript has its path, so its length on a transcript cannot be computed.
- **Not as gDNA:** it is spliced.

The hit is withheld. **The template's other hits stay.** Take a gDNA template that the aligner reports both as an
unannotated splice and as a contiguous alignment: it keeps the contiguous hit.

**Templates left with nothing.** A template with no usable hit left is withheld, with the reason **unannotated
splice**. It is reported in three classes, read off the genes whose span contains each side of the junction:

- within one gene;
- between genes (read-through);
- touching an intergenic region.

**No alternatives.** Rigel constructs no contiguous or partial reading of an unannotated junction. If the aligner
reports one, it is a hit like any other.

**What changes from today.** Today such hits are quantified through overlap candidates whose lengths keep the
junction's skip. The cost of withholding them is measured as an explicit policy cost:

- **on VCaP**, an earlier screen counted 114,519 RNA-library junction units (about 0.8% of spliced RNA units) against
  2,253 DNA-library artifacts;
- **on the minimap2 panel**, 90% of gDNA artifacts are unannotated long-anchor junctions.

Rigel does not assemble transcripts. Reads of unannotated isoforms, even ones that reuse annotated splice sites,
contribute to no transcript or gene.

**The blacklist** (`alignable`) is unchanged in wave 1. It removes a junction before step 4 sees it. It is itself a
repair of alignments, so the narrowed scope retires it once wave 4's matched comparison justifies that.

## 9. Step 5 — usable hits into calibration and the EM

### 9.1 Hits into units

- **Units.** Each usable hit is one EM hit. A template with more than one is a multimapper unit, and calibration holds
  it in pass 1, as with any multimapper today.
- **Intergenic hits are kept as hits**, so that `num_hits` and the accounting see them. Today the scanner drops a
  multimapper's hit that overlaps no transcript (`ISSUES: multimapper-intergenic-alignments`), while NH still counts
  it.
- **No consumer changes in wave 1.** The EM and calibration treat intergenic hits exactly as today until wave 2 sets
  their rule. The EM's gDNA term ignores a hit with no transcript candidate, and calibration's placements exclude
  intergenic hits.
- **For STAR input with complete pairs**, `num_hits` then equals NH, except where step 3 rejects a hit or step 4
  withholds one.

### 9.2 Why the existing EM combination is right for same-place alternatives

**RNA.** For each transcript, the EM keeps the hit where that transcript has the least overhang, then the highest
likelihood.

- Under transcript t, the molecule followed t's path. The right likelihood for t is the hit that best explains the read
  along t's path, counting its overhang and its mismatches once each.
- Today's rule takes the least overhang first, although overhang is already inside each hit's likelihood. So a hit
  with overhang 0 and NM 3 beats one with overhang 1 and NM 0 (`SCAN_CODE_HEALTH_REVIEW.md`,
  `merge-keeps-least-overhang`).
- For same-place alternatives with equal NM, both rules choose the same hit, so wave 1 keeps today's rule. Wave 2's
  combination rule decides.
- In the toy pair of §11 (a 3-base tie at TA1's donor, reported both spliced and contiguous), TA1 takes the spliced
  hit, with overhang 0, not 3, under either rule.
- TA2 and the gene's unspliced-RNA candidate appear only in the contiguous hit.

**gDNA.** Only unspliced hits are gDNA-capable, so the gDNA term comes from contiguous hits only. Today the term is
the average over gDNA-capable hits:

- **With one contiguous hit,** the average is that hit's own term. That is correct, and the spliced hit never enters
  it.
- **With two contiguous hits at different origins,** the average is half the correct sum. This is the wave-2 defect,
  bounded here by a factor of 2. It is counted in wave 1 and fixed in wave 2 by a sum over places (owner,
  2026-10-16).

**The prior.** The EM prior's count of calibrated fragments excludes multimapper units. Calibration holds these
templates too, so the two agree.

So wave 1 needs no EM change. The code-health review also found that pruning compares candidates against three
different references (`pruning-uses-three-references`). Like the merge rule, that belongs to wave 2.

### 9.3 Calibration's second pass: four corrections

Pass 2 draws one placement of a held template. The draw is weighted by:

- the pass-1 traffic at the objects that placement would deposit on: a junction for a spliced placement, an
  exon–intron boundary for a contiguous one;
- the placement's length likelihood;
- its strand.

**The all-held approximation (accepted, measured).**

- If every template that crosses a junction is held, the junction's pass-1 traffic is zero, and pass 2 deposits all of
  them as contiguous.
- The owner expects this to be rare: templates that map uniquely across the junction establish its traffic.
- The measurement is the number of templates whose spliced placement gets zero weight in pass 2, with their truth.
- The recommended STAR parameters will hold more templates, so this is re-measured on the STAR-aligned ladder.

**Four corrections, each its own arm.** Each also changes today's multimapper draws.

1. **The length law.** A placement whose observed junction is annotated is RNA, and should take the RNA law. Today the
   law is chosen from implied (mate-gap) introns only, so a spliced placement takes the global law.
2. **The density factor.** A placement with no contested region has no within-placement question, so its factor
   should be "not applicable". Today it reads 0: a hard zero against a placement that does have a gap hypothesis.
3. **The sequence evidence.** Pass 2 never reads mismatches, so a placement with 37 more mismatches is drawn as often
   as a perfect one.
   - Each held placement carries its NM into the deferred bank, as a new field in the bank's schema and its spill
     format.
   - Pass 2 multiplies the placement's weight by the EM's own mismatch term, α_mm^NM. Calibration and the EM then weigh
     hits by the same sequence evidence. This is the owner's rule of 2026-10-16: every hit is scored with its own
     mismatches.
   - At a tie the factor is equal and changes nothing.
4. **The motif rule** (`second-pass-restates-the-motif-rule`). Pass 2 restates the deposit's rule two other ways:
   - it passes an AMBIGUOUS motif into its sj lookup, which matches it at either strand, while the deposit trusts no
     splice on an AMBIGUOUS motif;
   - it decides that a placement is RNA from the motif tag, while the deposit and `A_p` decide it from an observed intron
     that resolves to an annotated sj. So an annotated splice whose tag is missing, or contradicts the annotation, gets
     gDNA's ½ in pass 2 while the deposit books it as spliced.

   Each placement's annotated-sj use and its strand are derived once, by the deposit's rule, and every factor reads
   them. It moves numbers only on rare held placements of tagged stranded libraries. The ladder cannot see it, because
   the simulator writes XS on every spliced read.

The first, second and fourth are local to `second_pass.py`. The third also touches the accumulator's deferred bank.

### 9.4 One length rule

The resolver gives a lone merged block its own length for every candidate. That skips the subtraction of the
candidate's introns inside the block, which `tx_frag_length` and DESIGN §3.4 apply (`single-block-length-shortcut`).
So one molecule gets two lengths, depending on whether its mates overlapped. The shortcut is removed in its own arm;
most of the candidates it affects already sit deep in overhang.

## 10. Step 6 — accounting for every template

**Outcomes.** Each template is exactly one of:

- **placed:** at least one usable hit;
- **reported chimera:** none usable, but a candidate pair resolves as chimeric;
- **unpaired:** none, with a reason;
- **excluded:** no eligible record, or the multimap option.

**Unpaired reasons.** Each candidate pair (an HI group, or an R1 × R2 combination) records how far it got, and the
template's reason is the furthest stage that any of them reached:

- missing mate;
- unresolved location;
- geometry;
- unresolved geometry;
- unannotated splice;
- length.

**A placed template records its route:**

- whether it reached the EM with candidates;
- whether it is intergenic only;
- in calibration, whether it was deposited or held, and its pass-2 draw;
- its hits by kind: same-place alternatives, multi-place;
- its records that were not used.

**The counts reconcile:**

> templates = placed + reported chimera + unpaired + excluded
>
> records = used + unused + records of unplaced templates

**Counters.**

- Today's `mate_unmapped` actually counts one-mate hit groups. It is renamed `one_mate_groups`, and a raw count of
  records flagged mate-unmapped sits beside it.
- New counters: `unpaired_alignments`, `nm_conflicts`, `withheld_chains`.
- Rates are reported against all templates and against eligible templates.

**The annotated BAM** writes every original record once. Tags carry the template's outcome and reason
(`ZC:unpaired` with the reason, `ZC:unused_alignment`), and stale tags are cleared.

**Fixes the code-health review found** (2026-10-16), each with its fixture:

- `proper_pairs` counts primary records that carry the proper-pair flag, so it matches its record-level label
  (`proper-pairs-counts-hits`).
- A gated fragment is annotated as gated, with its reason, not as intergenic, and the writer refuses duplicate
  frag_ids (`gated-fragments-stamped-intergenic`).
- The intergenic append is deleted, so an intergenic pair on two references is no longer counted and annotated
  chimeric (`intergenic-append-is-dead-and-mislabels`).
- The annotation pass's read-group count must equal the scan's, or the run stops (`annotation-pass-can-drift-a-group`).
- A multimapper's ZS ranks all five splice types (`multimapper-splice-type-drops-implicit`).
- Every scan counter is asserted exactly, on one hand-built BAM, at one worker and at four (`scan-counters-untested`).

## 11. Tests

One suite per step. Each fixture is one of two kinds:

- **invariant:** passes on the baseline and after;
- **defect demonstration:** fails on the baseline, passes after the fix, and fails again when the fix is broken on
  purpose.

**The toy reference.** These are the annotation's exons, 0-based and half-open:

| transcript | strand | exons |
|---|---|---|
| TA1 | + | [1000,1200) [1500,1700) [2000,2300) |
| TA2 | + | [1000,1200) [2000,2300) |
| TB1 | + | [5000,5150) |
| TC1 | − | [2200,2400) [2600,2800) |

**Expected hits.** The emission values come from the baseline resolver (`counterexample_baseline.py`). In the table,
"x/y/z" means exon bases / overhang / fragment length for that candidate.

| suite | case | records (F / R) | expected |
|---|---|---|---|
| pairing | same-place alternatives, as STAR would report them | HI 1: `97M300N3M` at 1103 / [1550,1650); HI 2: [1103,1203) / [1550,1650) | two hits, one unit. HI 1 gives TA1 200/0/247. HI 2 gives TA1 197/3/247, TA2 97/103/547 and unspliced A 200/0/547, with gDNA at 547 |
| pairing | mates disagree, acceptor side | F `50M800N50M` at 1150 / R [1998,2098) | one hit, as reported: TA2 148/2/148 |
| pairing | mates disagree, donor side | F [1103,1203) / R `50M300N50M` at 1150 | one hit, as reported: TA1 147/3/147 |
| pairing | clips hiding TA2's junction | F `50S50M` at 2000 / R `50M50S` at 1150 | placed along TA2's path, L 100 |
| pairing | outward, unclipped | F [1300,1400) / R [1100,1200) | unpaired: geometry |
| pairing | a lone secondary, no partner in its HI group | one mate only | an unpaired alignment, counted; no hit |
| pairing | a primary whose mate is unmapped | one mate only | unpaired: missing mate |
| pairing | two HI groups at different loci | both asserted | exactly the two asserted hits; no cross pair |
| pairing | per-mate secondaries, no HI (minimap2 style) | R1 at A and B, R2 at A′ and B′ | two hits, (A, A′) and (B, B′); the cross pairs fail validation |
| pairing | no HI, a partner inside the bound but beyond `max_frag_length` | the long-intron witness ([1050,1150) / [2100,2200)), its mates reported per mate | paired: the genomic span is 1,150 against the default 1,000, but L is 350 along TA2 and 550 along TA1 |
| pairing | no HI, a partner on another reference, on the same strand, or beyond the bound | R1 and R2 per mate | never tried; counted |
| pairing | NH above the number of HI groups | truncated output | the run stops, with the remedy |
| unannotated | a novel junction inside gene A, with a contiguous alternative | HI 1: junction [1250,1500); HI 2: contiguous | HI 1 withheld; HI 2 kept; the template is placed, unique |
| unannotated | a junction between genes A and B, alone | `50M2800N50M` at 2200 / [5100,5200) | withheld: unannotated splice, read-through |
| pairing | one-base overlap, both outer boundaries, long-intron witness | as in revision 2 of the first design, §12 | as there |

**Other cases:**

- **Locations:**
  - duplicates within and across HI groups;
  - `5S95M` against `95M5S`;
  - conflicting NM or motif;
  - missing NM;
  - hard-clipped carriers;
  - chains without SA;
  - adjacent H and S;
  - a record with no M/=/X base;
  - a BAM without NM: counted and warned;
  - mates on two unindexed references;
  - NH or HI typed as a string;
  - CIGAR D, I, =, X and H, before and after an N (`no-indel-cigar-fixtures`).
- **Pairing and the resolver:**
  - without HI, two copies inside one gene's introns: only the geometric pairs;
  - an antisense overlap with the strands installed by the index: `ambig_strand` is 1;
  - one fixture per use of the anchor tolerance (§7.5);
  - a merged block over an intron: the candidate's path length (§9.4).
- **Accounting:**
  - one fixture per reason;
  - a geometry failure beside a too-long consistent pair;
  - every record written once;
  - stale tags cleared;
  - the counts reconciling on a mixed BAM;
  - the counters BAM, every counter exact;
  - a gated fragment's annotation;
  - an annotation pass out of step with the scan: refused.
- **Calibration:**
  - an all-held junction;
  - the length-law correction;
  - the density-factor correction;
  - a held template whose placements differ only in NM, drawn in proportion to α_mm^NM;
  - NM surviving the deferred bank's spill and reload;
  - the motif rule: an AMBIGUOUS motif, a missing tag, and a tag that contradicts the annotation.
- **Code health:**
  - each cleanup stage passes `rename_identity.py --check`;
  - a deliberately moved output value in a stage makes `--check` fire.

**Break tests.** Each must fail its suite:

- restore the one-mate hit;
- restore `max(NH, hits)`;
- skip pair validation;
- form a pair across two HI groups;
- take only the primary pair on a BAM without HI;
- bound the no-HI walk by `max_frag_length` alone (the long-intron witness must fail);
- bound candidates by exon ends;
- prune by `L_t`;
- walk clips contiguously only;
- fail containment on a mate's bases inside its partner's skip;
- withhold the whole template when only one of its hits is unannotated;
- withhold a hit with only annotated junctions;
- drop intergenic hits from `num_hits`;
- drop NM from the deferred bank on spill;
- let a cleanup stage move a number;
- read a missing NM as 0 without counting it;
- pair without HI by transcript sets;
- file unindexed references under one id;
- leave the strands to the scan;
- stamp a gated fragment intergenic;
- accept an annotation pass whose read-group count differs;
- decide "RNA" in pass 2 from the motif tag;
- keep the single-block length shortcut.

## 12. Measurement and gates

**Arms,** one mechanism each, in order:

1. the cleanup stages, each a proven numeric no-op. The first is the scan-cache key: it must cover the splice
   blacklist and the BAM before any arm reads a cache (`scan-cache-key-misses-blacklist-and-bam`);
2. locations, a numeric no-op except for the enumerated input defects (with the NM count, the unindexed references'
   identity and the tag checks);
3. the resolver set up by the index;
4. pairs by HI, and the bounded R1 × R2 search without HI, with validation and the new `num_hits` (the one-mate
   fix);
5. the path witness for clipped and outward pairs;
6. chains withheld;
7. unannotated hits withheld, hit by hit;
8. each use of the anchor tolerance outside implicit splicing, one arm per use (§7.5);
9. the single length rule (§9.4);
10. the accounting fixes (tags and counters; no quantification number moves);
11. the length-law correction;
12. the density-factor correction;
13. mismatches in pass 2;
14. the motif rule in pass 2.

**Readings of each arm:**

- **the setup:** the gDNA law and the ruler frozen, with a hash of the arrays each consumer received, and also each
  arm's own; per stratum; scored against truth;
- **the errors:**
  - both errors (artifact RNA mass, and genuine RNA lost);
  - transcript and gene error;
  - calibration against the oracle;
  - the pools, including the unspliced-RNA over-call;
- **the census:** a per-template transition census with truth, and the hits counted by kind;
- **calibration:** the second pass's all-held count;
- **cost:** scan time and memory.

**Inputs:**

- **the oracle ladder**, as the identity check. Its alignments are perfect, so arms 2–7 must leave it unchanged apart
  from enumerated cases.
- **the local minimap2 panel** (five conditions, simulator truth). Its per-mate secondaries are paired by the R1 × R2
  search (§7.1).
- **VCaP** (STAR, production parameters), for transitions and library provenance.
- **the STAR-aligned ladder conditions** from the STAR study, once delivered, at both production and recommended
  parameters. This is the decisive input, because the recommended parameters are Rigel's supported input.

**The implementation gate:**

- the code-health report reviewed (done, 2026-10-16), and its cleanup landed;
- this plan reviewed and revised;
- every suite's fixtures written, and failing where they demonstrate a defect;
- the injection hooks proven to reach their consumers.

**The end-of-wave gate:**

- every suite and break test passes, and the full suite passes at the re-derived baseline;
- oracle identity holds, apart from enumerated cases;
- with the law and the ruler frozen, no in-scope stratum gets worse, except by the measured, accepted policy costs
  (unpaired templates, withheld unannotated hits);
- the owner reassesses before wave 2.

## 13. Open questions

1. **The no-HI search's bound** (§7.1). Is it complete: can any pair that §7.2 accepts lie beyond it?
2. **The per-hit unannotated rule** (§8), and the blacklist's interaction with it.
3. **The existing EM combination for same-place alternatives** (§9.2), and its factor-2 case.
4. **Calibration's all-held approximation and the four corrections** (§9.3), now that the recommended parameters will
   hold more templates.
5. **The code-health step as carried out** (§5). Is anything left that steps 2–6 depend on?
6. **The anchor tolerance** (§7.5). Which of K's uses are implicit splicing?
7. **The junction-strand tag per record** (§6). When a record carries both XS and ts and they disagree, which stands?
   Is a library-level warning still needed when most spliced records carry no tag?

The output names are settled (owner, 2026-10-17; `OUTPUT_NAMES_PROPOSAL.md`), and the renames stage applies them before
step 2.

## 14. How this answers the earlier reviews

| finding | here |
|---|---|
| tied alternatives meet best-hit RNA and mean gDNA | correct for same-place alternatives (§9.2); remote ones wait for wave 2's sum |
| preference is not identification | the scan stage never chooses between hits; the aligner proposes, the EM decides (§3) |
| one owner of the sequence evidence | the aligner's NM for each hit, in the EM and now in pass 2 (§9.3) |
| calibration can zero a held RNA placement | accepted by the owner and measured; four corrections (§9.3) |
| the clip walk is not exact | the path witness walks along annotated paths; unresolved geometry otherwise (§7.2) |
| location identity; the ledger | §6 and §10 |
| read-through membership | subsumed: unannotated hits are withheld; the read-through class is reporting only (§8) |
| revision 1's open questions: the generation tolerance, K where bases are missing, multi-junction enumeration | moot: wave 1 constructs no reading |

## Appendix A — what left wave 1, and where it went

- **Alignment scrutiny** (moving bases, Δm, the generation tolerance, reconciled readings, multi-junction enumeration):
  removed. Near-tied alternatives are the aligner's job. The STAR study (`star_parameters/PROMPT.md`) establishes the
  parameters Rigel will recommend.
- **The genome stored in the index:** dropped for now (owner, 2026-10-16). Nothing in the plan reads bases.
- **The splice-artifact machinery of 2026-10-14/15** (the local tie table and its second hit, the hold rule, the
  mismatch emission factor): superseded by the aligner's alternatives and their NM, pending the STAR report.
- **R1 × R2 enumeration:** kept only for aligners without HI, as a bounded search (§7.1); STAR input is paired by HI.
- **A dual-aligner input:** a spliced aligner and a genome aligner, merged by read name into one BAM for Rigel. This is
  the owner's eventual fallback if STAR cannot be parameterized, and the STAR study prices it.
- **Remote copies** (the index-time copy table): wave 2, and only if VCaP shows that STAR misses them.

## Appendix B — evidence pointers

- **The code-health review:** `SCAN_CODE_HEALTH_REVIEW.md` and `scan_code_health_slices/`; the cleanup is
  `41899ae4..3941a0dd`.

- **Diagnosis tables:** revision 2 of the first design, Appendix A, with `paired_end_repair_evidence/arms.json`.
- **The mate-disagreement census:** `census_*.json`.
  - End-extensions within 3 bases are 495 mRNA against 58 artifact templates on minimap2 OFF.
  - On VCaP they are 2,165 RNA-library against 5 DNA-library templates.
  - Recount after pair validation.
- **Emission values on the toy annotation:** `counterexample_baseline.py`, confirmed independently by the first
  review.
