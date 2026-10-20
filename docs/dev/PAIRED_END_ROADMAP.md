# The paired-end repair: roadmap to release

2026-10-18 · code base `07eba0fc` · sandbox roadmap; the fix order (owner, 2026-10-18) below supersedes the waves'
order, and its first two changes are implemented.

This file holds the waves, what each must deliver, and the gate between them.

- **The plan for wave 1:** `PAIRED_END_WAVE1_PLAN.md` (revision 4).
- **The diagnosis behind it:** `PAIRED_END_REPAIR.md` (revision 2, superseded for planning) and its independent
  review, `PAIRED_END_REPAIR_INDEPENDENT_REVIEW.md`.
- **The STAR study's brief:** `star_parameters/PROMPT.md`.

The owner's direction (2026-10-16) is to go slowly, one wave at a time: design it, review it, build it, measure it,
and reassess before the next.

## Where we are

- **The bug.** Rigel's BAM scanner built "fragments" from a single mate. A lone secondary alignment, a primary whose
  mate was unmapped, or a half-filled HI group each became a read-length fragment, a false multimapper and a 100-bp
  length deposit.
- **The first repair** proved that, but transcript error rose. The diagnosis found two downstream defects the bug had
  been hiding:
  - the gDNA fragment-length estimator is ill-conditioned on aligned input;
  - the capture ruler finds capture on uncaptured libraries.

  With those two held fixed in diagnostic runs, the pairing repair improves every in-scope condition that has gDNA.
- **Splice artifacts** remain the release-blocking false positive: gDNA aligned as spliced, and certified as RNA.
  Unannotated splicing is handled inconsistently today.
- **The scope is narrowed** (owner, 2026-10-16). Rigel reads alignments; it does not repair or re-align them.
  - Near-tied spliced and contiguous alignments are the aligner's to report.
  - A parallel study on the cluster has found STAR parameters that do so (owner, 2026-10-17), and Rigel will
    recommend them.

## Progress

Kept current until wave 4 closes (owner, 2026-10-16). ☑ done · ◐ in progress · ☐ not started.

**The fix order** (owner, 2026-10-18). Rigel is fixed in pipeline order: the scan, then calibration, then after
calibration, one small change at a time, built on the STAR panel at the recommended parameters. It folds wave 1's
steps 2–6 and waves 2 and 3 into one list; the waves below keep their record. Each change gets:
1. a failing test on a hand-built case;
2. the fix alone;
3. each gate broken on purpose and seen to fire;
4. the STAR panel against truth, per stratum, with the census by true origin; the oracle ladder for regressions; the
   full suite and preflight;
5. commit, push and record, then the next change.

*Phase A, the scan.* Its two goals come first (owner, 2026-10-18): keep competing spliced and unspliced alignments,
and give every accepted alignment a valid fragment length, so no unpaired reads. The rest is bundled after them because
the code is open; each is less critical.
1. ☑ The junction-strand tag read per record (2026-10-18; `ISSUES: sj-strand-tag-chosen-from-the-first-reads`,
  closed). Pinned and fractional, against baselines from the same code: the STAR panel at the recommended
  parameters 15/15 quant_accuracy rows and the oracle ladder 48/48 identical in every field; suite 3,050. The panel
  caches had recorded the old resolved tag, so the cache check now reads the scan settings (preflight and
  `build_scan_cache.py`'s skip did not), and every panel was re-scanned once: only float sums moved, by the scan's
  thread-order noise (two threaded scans of one BAM differ by up to 4e-15 relative; pinned scans are identical)
2. ☑ Transcript strands installed by the index (2026-10-18): a fragment resolved without a scan sees them, so an
  antisense overlap reads strand-ambiguous; pinned and fractional, the STAR panel 15/15 and the oracle ladder 48/48
  quant_accuracy rows identical in every field; suite 3,052
3. ☑ The pairing question: STAR's HI pairs against Rigel's own search, each scored against the molecule's truth
  (`ISSUES: star-hi-pairs-against-rigels-own-pairing`; owner, 2026-10-18). ☑ Template level (2026-10-18, 47.9 M
  templates, per-base truth): the validator rejects no correct pair and removes 1.40 M incorrect STAR pairs; HI misses
  the correct pair in 0.31 % of multi-hit templates, and ignoring HI adds 475 k pairs, 98 % incorrect. ☑ End to end
  (2026-10-19, four arms): the search wins every capture-OFF condition with the capture efficiencies at 1; the owner
  took it for every input, HI not read (`DESIGN.md` §1.1)
4. ☐ The gates' own proof (the plan's step 1b), before the first change whose identity claim leans on those gates
5. ◐ Location identity and aliases (D and N apart, M/=/X runs coalesced), which duplicate hits need. ☑ One place per
  alignment, its aliases written with its hit (2026-10-19, with the search). ☐ NM conflicts between aliases (item 12)
6. ☑ Complete pairs only: the search for every input (2026-10-19; HI not read, so no HI group, mixed template or
  duplicate group remains; a truncated output is counted against the records reported per mate). ☑ One-mate hits dropped, counted as `missing_mate` (2026-10-18): the STAR panel 15/15 and the ladder 48/48
  identical, as neither holds one; suite 3,063. ☑ Duplicate groups collapse, by construction in the search (held
  2026-10-18 alone: it tripped the EM's fragility at exact paralogs, g00's SIK1B 33,100 → 6,255 of 33,063)
7. ☑ `num_hits` counts the template's valid hits, not `max(NH, hits)` (2026-10-19): against the one-mate arm on the
  STAR panel, transcript error +0.08 % unstranded OFF, +0.16 % stranded OFF, −0.15 % stranded ON (gene −0.01 / +0.15
  / +0.07 %); the census moves by under 2.3 k per condition; the ladder identical; suite 3,067. Truncated outputs
  (NH above the HI groups) are counted: none on the STAR panel
8. ◐ Without HI (minimap2), every same-reference pair validated (`GEOMETRIC_PAIRING_DESIGN.md` revision 2, approved
  2026-10-19: the owner's rule, one transcript over both mates, in place of the span bound). ☑ Step 1, the witness
  rule (item 9). ☑ Step 2, the search, for every input (2026-10-19; item 3 decided it); suite 3,091
9. ◐ Each hit validated as one molecule along a path, with one witness transcript: K by the resolver's rule; outward
  and too-long pairs rejected; clips walked. ☑ The validator, clips walked along the genome (2026-10-18; the C++ counts
  match the prototype's on all five STAR conditions within their duplicate HI groups; suite 3,060). Against main on
  the STAR panel: transcript error −2.87 % stranded OFF, +1.18 % unstranded OFF, +0.29 % stranded ON (gene −3.33 /
  +3.27 / −0.62 %); the ladder identical. ☑ One transcript over both mates (2026-10-19, the owner's rule): the witness
  and a spliced hit's junction holder overlap both mates; transcript error −0.74 / −0.66 / +0.86 %, gene +0.06 /
  +3.82 / +0.31 % (g00's paralogs), true gDNA counted RNA −24 k capture-ON; the ladder identical; suite 3,086.
  ☐ Clip walks along a transcript path (minimap2 input); ☐ one-mate hits and
  mates with several records, which pass through unjudged
10. ◐ The annotated-path rule: a spliced hit is usable only if one annotated transcript holds all its junctions; its
  contiguous competitor stays. ☑ One junction set per hit, every record's (2026-10-19, 69d003ea): the STAR panel 15/15
  and the ladder 48/48 identical. ☑ The rule (2026-10-19, aff4be8e): transcript error −10.17 % unstranded capture-OFF,
  +4.55 % stranded capture-OFF, +1.28 % stranded capture-ON (gene −6.16 / +5.16 / +2.24 %), the ladder identical;
  landed regardless (owner), the losses in `ISSUES: the-annotated-path-rule-costs-the-stranded-strata`. ☑ The
  intergenic RNA pool deleted. ☑ `SPLICED_UNANNOT` deleted (2026-10-19): the resolver reads `junction_holders`
  and refuses a spliced fragment with no holder; the STAR panel and the ladder identical
11. ☑ One length rule for overlapping mates (2026-10-20, rule B, `ONE_LENGTH_RULE_DESIGN.md`): a candidate's length is
  the length of the path it implies, measured by the deposit's own `path_length`; the endpoint rule deleted. Against
  main on the STAR panel, pinned and fractional: transcript error −0.004 % unstranded OFF, −0.002 % stranded OFF,
  −0.058 % stranded ON (gene −0.003 / +0.001 / −0.030 %), the gDNA pool within 1 fragment; with every capture
  efficiency at 1, −0.043 / +0.057 / +0.005 %; calibration identical (every locus prior, the models, all 56 scan
  counters); the ladder within ±0.52 % per row; the pinned scan +1.3 %. The drop in gDNA counted as RNA the design
  expected did not appear: the candidates whose length changed already lose to their overhang penalty. Suite 3,092

*The bundle*
12. ☐ The rest of one location per alignment: NM conflicts and missing NM counted; unindexed references apart; tag
  types checked; count widths that cannot wrap
13. ☐ Supplementary chains and hard-clipped 5′ ends withheld, and counted
14. ☐ Each use of K outside implicit splicing, one at a time
15. ☐ Accounting: one outcome per template, each record's status from its template's outcome, the counters and the
  annotated-BAM tags; the annotation pass checked against the scan by identity

*Phase B, calibration*
0. ☐ First, the highest priority (owner, 2026-10-19): the capture reader reads false gDNA as capture on aligned input
  (`ISSUES: the-reader-reads-false-gdna-as-capture`); with every efficiency at 1, g00's transcript error falls 3.44 M →
  0.25 M and every capture-OFF condition's by 18–93 %, while capture-ON needs them (+232 %)
1. ☐ The second pass: an observed annotated splice takes the RNA length law
2. ☐ The second pass's weights: one placement factor times a unit-free share within the placement
3. ☐ The second pass reads NM, carried in the deferred bank
4. ☐ The second pass's motif rule, derived once as the deposit derives it
5. ☐ The gDNA fragment-length estimator on aligned input
6. ☐ The capture ruler's false capture on uncaptured libraries: item 0
7. ☐ The rule for intergenic placements
8. ☐ Observable opportunity, and the clip-length model

*Phase C, after calibration (the EM)*
1. ☐ The per-transcript merge chooses by full likelihood (today the least overhang first)
2. ☐ The gDNA term a sum over places (the owner's ruling; today the mean)
3. ☐ Intergenic hits enter the EM, release-blocking (owner, 2026-10-19): a fragment's valid intergenic and genic
  alignments compete (`ISSUES: multimapper-intergenic-alignments`)
4. ☐ The prior's denominator (today it excludes the multimappers calibration deposited)
5. ☐ Pruning against one reference

**Wave 1, the scan stage** (plan revision 2 accepted 2026-10-16; revision 3 reviewed 2026-10-17; revision 4 for a
targeted re-review)

- ☑ Step 1, the code-health review and its cleanup (2026-10-16/17)
  - ☑ the baseline re-derived on a clean checkout at `41899ae4` after a rebuild: 0 failed / 2,959 passed / 0 skipped
    / 0 xfail, 2,959 collected (2026-10-16; CLAUDE.md updated, uncommitted)
  - ☑ the review, in nine read-only slices: the scan path in depth (six), the rest of `src/` swept (three); the
    slice reports are in `scan_code_health_slices/` (2026-10-16)
  - ☑ the report, `SCAN_CODE_HEALTH_REVIEW.md` (2026-10-16)
  - ☑ the owner's selection: all of it (2026-10-16). The scan-path defects joined steps 2–6, the later-wave findings
    became ISSUES entries, and the anchor tolerance is kept for implicit splicing with every other use audited
  - ☑ the cleanup stages, each a proven no-op or measured, in the report's order (2026-10-16/17):
    - ☑ the scan-cache key (blacklist and BAM): committed `b6aec8a4` and pushed, 2026-10-16 (suite 2,968 / 0 failed;
      bit-identical on all three references; preflight clean). Every existing cache is refused until the re-cache
    - ☑ the words: comments, docstrings and help text, 45 files: committed `a7c5ebe3` and pushed (suite 2,968;
      bit-identical on all three references)
    - ☑ dead code, 51 files: committed `26d37589` and pushed (suite 2,964: five cases deleted with their code, one
      gate added for the raw fragment-length model's refusal; bit-identical on all three references)
    - ☑ one definition each, 49 files: committed `c1156e86` and pushed (suite 2,984: twenty gates added, the native
      ones broken one by one and each fired; bit-identical on all three references). The exonic strand vectors
      are gone, so the renames stage's `exonic_truth` item is moot
    - ☑ contracts, 23 files: committed `1b4a46c6` and pushed (suite 3,025: forty-one gates added, the four native
      ones broken in one rebuild and each fired; bit-identical on all three references)
    - ☑ one re-cache, 25 files: committed `d5a08e46` and pushed (suite 3,025; measured, since a canonical-key field
      went: on all three references exactly the two deleted region banks are gone, nothing appeared, and the
      transcript table and every calibration array are bit-identical). The contiguous reaches' deletion waits
      for an index-format change, so it joins the renames stage's index rebuild. Every panel's caches
      rebuilt (13 configs, 262 conditions certified) and the three identity references re-frozen at `d5a08e46`
      (proven equal to the old freeze but for the two deleted banks). Preflight's new key check then caught
      a defect from the first stage: an oracle part recorded its scratch split file as its BAM, which a
      same-named condition in another panel overwrites (12/16 ladder, 0/30 test chromosome refused). Fixed
      in the three writers (a part records the BAM it was split from) with a gate; every oracle part
      rebuilt, and preflight reads every cache valid (16/16 and 30/30, 262 conditions certified); committed
      `30357fb3` and pushed (suite 3,026)
    - ☑ stale benchmark data deleted for space, at the owner's request (107.6 GiB; text receipts kept)
    - ☑ renames, to the convention the owner ruled on 2026-10-17 (`OUTPUT_NAMES_PROPOSAL.md`, five steps):
      ☑ 1. the identity digest ignores column names (`8ebf71ea`), references re-frozen there; ☑ 2. the internal
      names (`fd012918`, bit-identical); ☑ 3. the outputs to the convention, `summary.json` schema 4 (`1fa1c143`:
      every output file compared on one ladder condition, g50 stranded capture-ON: each renamed field keeps its
      values bit for bit and only the defined relabellings move; suite 3,045);
      ☑ 4. the index and scan names, index format 10 (`dc736e39`: every index rebuilt and checked file by file
      against its format-9 build, every panel re-cached, bit-identical on all three references; the synthetic
      spans keep `RIGEL_NRNA_` ids, which the panels' read names carry); ☑ 5. the simulator's truth tables through
      one writer (each quantity once, each path under one key; reads and kept truth identical in an A/B; the
      renames that need a re-simulation are held)
    - ☑ the public repository scrubbed of local paths, cluster details and lab identifiers, owner's request
      (`8401aeba`, suite 3,048); docs/dev keeps markdown only (`20f187ab`); the data a commit of mine pushed was
      taken out of history by rewriting that one commit (owner-authorized)
    - ☑ constants, to the owner's rulings of 2026-10-17: the six that needed none (`bb3967cb`); the resource and
      display numbers (`0c754afc`, bit-identical); the row-tiled fits made independent of their tile height
      (`4bcc72d9`: measured on the ladder against truth, fix / before 0.998-1.001 per stratum, 0.996-1.007 per
      condition); the tag detector's
      1,000-read cutoff is repaired in step 2 instead (the junction-strand tag read per record)
- ☑ the readiness review of revision 3 (2026-10-17, `PAIRED_END_WAVE1_READINESS_REVIEW.md`): revise before
  implementing, eight findings; each verified against the source and accepted
- ☑ revision 4, answering the eight findings (2026-10-17)
- ☑ its targeted re-review (2026-10-18, `PAIRED_END_WAVE1_TARGETED_REVIEW.md`): two blockers, that one transcript
  must witness a path's observed and implied junctions alike, and that a mixed-HI template must keep its known pairs;
  two fixes before implementing, stronger evidence gates with the normalization contract, and each record's status
  from its template's final outcome; the six open questions answered with recommendations
- ☑ revision 5 (2026-10-18), answering the targeted review: one witness transcript per path; a mixed-HI template keeps
  its asserted pairs and invents none; unique hypotheses and log weights in pass 2; each record's status from its
  template's outcome; the six open questions adopted or put to the owner. ☐ The reviewer's check of it, optional,
  before the changes it governs
- ☐ step 1b, the gates' own proof: each gate broken on purpose and seen to fail (plan §12)
- ☑ the STAR panel (2026-10-18): the five conditions of the minimap2 panel aligned on the cluster by STAR 2.7.11b at
  production and recommended parameters; all twelve BAMs collated, complete and carrying NM. Error-injected copies not
  made (plan §13, open question 3)
- ☑ rustar approved for local alignment (owner, 2026-10-18): 99.9% of templates identical to STAR on a 100k-pair
  sample; wrapper in `$RIGEL_RUNS/rustar/`
- ☑ the first measurement (2026-10-18; `$RIGEL_RUNS/star_panel/first_measurement.md`), today's tree, pinned and
  fractional, scored against truth:
  - on STAR alignments 18–22% of the g50 conditions' gDNA is counted as RNA, mostly in nascent spans, against near zero
    on the oracle alignments; 99.5% of it is unspliced (the deconvolution, not splice artifacts);
  - the recommended parameters cut annotated-junction artifacts 25–55%, but raise unannotated-junction ones 12–60% and
    spliced intergenic ones about twofold, which the annotated-path rule (step 4) withholds;
  - with today's tree they raise gDNA counted as RNA 2–6% and transcript error 7–24% per stratum;
  - ☐ still to take: the census of the deferred defects (the merger, the gDNA mean, the prior's cohort)
- ☐ the fixtures, with the exhaustive pairing oracle
- ☐ Step 2, mate locations
- ☐ Step 3, pairing: HI groups, completeness per mate and duplicates; every R1 × R2 pair without HI; the anchor
  tolerance's uses audited
- ☐ Step 4, the annotated-path rule, hit by hit
- ☐ Step 5, the hand-off, the single length rule, and the four second-pass corrections (one arm each)
- ☐ Step 6, accounting
- ☐ the end-of-wave gate, and the owner's reassessment

**Wave 2, multimapping**

- ☐ the design: the combination rule, intergenic hits, pass-2 weights, fewer unique templates
- ☐ its review, the build, and the gate
- ☐ only if the STAR report shows the gap: remote copies; the clipped-bases rule

**Wave 3, calibration and downstream**

- ☐ the gDNA length estimator
- ☐ the capture ruler's false capture
- ☐ the aligned gDNA/RNA split
- ☐ the clip length model
- ☐ observable opportunity

**Wave 4, validation and release**

- ◐ the STAR study (cluster session): ☑ parameters that report both readings (owner, 2026-10-17); ☐ its report
- ☐ the STAR-aligned ladder, at the recommended and the production parameters
- ☐ independent artifact-containing libraries, and VCaP
- ☐ the fragment-length-gap panels and the zero-gDNA controls
- ☐ the blacklist-retirement comparison
- ☐ the manual's STAR command, and the BAM-header check
- ☐ the full suite, caches, persistence, docs and release steps

## The owner's rules that every wave keeps

- **Three populations:** gDNA, RNA+ and RNA−.
- **Both mates.** Every hit places both mates, and is scorable only if its fragment length can be computed. A template
  without such a hit is withheld and counted.
- **The emission model is kept:**
  - candidates by overlap;
  - an overhang penalty;
  - a mismatch penalty;
  - per-candidate lengths.
- **Plausible RNA and gDNA alignments of one read pair compete in the EM,** as the aligner reports them.
- **Rigel reads alignments and never repairs or re-aligns them.** It recommends STAR with tested parameters.
- **Annotated transcripts only.** Rigel does not assemble transcripts or quantify chimeras.
- **No new constant.**
- **Fix everything before release.**

## The waves

| wave | what it delivers | depends on |
|---|---|---|
| 1. The scan stage | a code-health review of the scan path, and its cleanup; pairs by HI (STAR) or by a bounded R1 × R2 search (aligners without HI), each validated; the unannotated rule, hit by hit; every template accounted for; four corrections to calibration's second pass | nothing |
| 2. Multimapping | how a template's hits combine, for RNA and gDNA alike, in the EM and in calibration | wave 1; the STAR study, which sets how many templates this touches |
| 3. Calibration and downstream | the gDNA length estimator; the capture ruler's false capture; the aligned gDNA/RNA split; the clip length model; observable opportunity | wave 1 (the order against wave 2 is reassessed after wave 1) |
| 4. Validation and release | the STAR-aligned ladder at the recommended parameters; independent artifact-containing libraries; the blacklist-retirement comparison; the manual's STAR recommendation; the full suite and release steps | waves 1–3; the STAR report |

### Wave 1: the scan stage

The scan stage reads the BAM faithfully, and turns every template into validated hits or a disposition. It has six
steps:

1. **Review the scan path's code health, and clean it up.** Each cleanup stage is proven to change no number.
2. **One mate location per alignment record.**
3. **Hits from pairs:** STAR's HI groups, or every R1 × R2 pair for aligners without HI (minimap2). Each is validated
   as one molecule.
4. **The annotated-path rule, hit by hit:** a spliced hit is usable only if one annotated transcript holds all its
   junctions. A template keeps its other hits.
5. **The hand-off to calibration and the EM,** with four corrections to calibration's second pass:
   - a spliced placement takes the RNA length law;
   - each candidate carries one placement factor and a share within its placement, so a placement with nothing
     contested is no longer zeroed, and the draw does not depend on units;
   - pass 2 weighs each hit's mismatches as the EM does;
   - pass 2 decides a placement's splice and strand by the deposit's rule.
6. **One disposition per template;** every record written once.

**Gate.**

- Every suite passes, and every deliberate break test fires.
- The full suite passes at the re-derived baseline.
- The oracle ladder is identical, except for enumerated cases.
- With the gDNA law and the ruler held fixed, no in-scope stratum gets worse, except by the explicitly accepted policy
  costs: the unpaired templates, the withheld unannotated hits and the contradictory mates. Each cost is measured
  and reported, on the STAR panel (simulated truth) and, where simulation cannot show it, on VCaP.

### Wave 2: multimapping

Today, when a template has several hits, the EM builds one unit from them in two ways:

- **RNA:** each transcript keeps the likelihood from the single hit where it fits best.
- **gDNA:** the term is the **average** of the per-hit gDNA likelihoods. Only hits that gDNA can explain and that
  overlap at least one transcript count.

A worked example:

- A template fits gene A perfectly (gDNA likelihood 1, in relative units) and a paralog B with two mismatches (gDNA
  likelihood 0.01).
- Today's gDNA term is (1 + 0.01) / 2 ≈ 0.5.
- gDNA could have come from either place, so its probability should be the **sum** over places, each weighted by that
  place's genomic prior: about 1.01 here, not 0.5.
- The RNA candidates of gene A, meanwhile, keep their full likelihood.

Adding a poor second place therefore halves gDNA but leaves RNA untouched, which biases multimappers toward RNA. A hit
in an intergenic region (no transcript) is dropped from the gDNA term altogether, so its gDNA hypothesis disappears.

**Why this wave is now critical.** With the recommended STAR parameters, every near-tied template becomes a
multimapper. This wave decides how all of them are scored. It derives:

- **the combination rule** (owner: a sum, not the mean):
  - gDNA sums over distinct places, each with its genomic prior, its capture efficiency and its observable
    opportunity;
  - alternatives at one place count once;
- **intergenic hits,** which keep their gDNA hypothesis;
- **the matching weights in calibration's second pass;**
- **the calibration consequence of fewer unique templates.**

Only if the STAR report shows the gap:

- **remote copies from the index.** At index build, the annotated transcripts are aligned to the genome unspliced
  (minimap2 via `mappy`), which proposes the remote copies a fragment could have come from. This applies only if STAR
  misses retrocopies. A 2026-10-14 measurement found such copies for 2,658 of 404,168 annotated junctions, and they
  recalled 95.5% of VCaP's retrocopy artifacts.
- **a rule for clipped bases,** so that an alignment cannot win by explaining fewer of the read's bases. This applies
  only if the recommended alignment keeps soft clipping.

**Gate.**

- Finite enumeration matches the combination rule.
- The per-stratum comparison of both errors on the STAR-aligned ladder.

### Wave 3: calibration and downstream

- **The gDNA length estimator.** Its two-pool contrast is ill-conditioned when the estimated pool purities are low.
  Its internal density is 7.5–72× too low on aligned input.
- **The capture ruler on uncaptured libraries.** At zero gDNA it carries 94% of the transcript error. At g50 capture
  OFF, about 2% of transcripts are off by more than 2×.
- **The aligned gDNA/RNA split** in the unstranded stratum.
- **The clip length model:** molecule ends inferred from 5′ soft clips, along each candidate's path.
- **Observable opportunity:** not only N bases. Two mates can span an N gap, and repeats can be unobservable.

Each is designed and measured on the repaired scan.

### Wave 4: validation and release

- **The STAR-aligned ladder** at the recommended parameters, which are Rigel's supported input; and at the production
  parameters, to price what a default BAM gets.
- **Independent artifact-containing libraries, and VCaP.**
- **Both fragment-length-gap panels, and the zero-gDNA controls.**
- **The blacklist-retirement comparison.** The blacklist repairs alignments, so the narrowed scope retires it once the
  comparison justifies it.
- **The manual,** with the recommended STAR command, and a check of the BAM header's `@PG` line that warns when a BAM
  was aligned otherwise.
- **The full suite, the cache and persistence checks, the docs, and the release steps.**

## Ideas held for evaluation

- **A dual-aligner input** (owner, 2026-10-16). Every read pair is aligned twice, by a spliced aligner and by a genome
  aligner, and the two are merged by read name into one BAM for Rigel. Rigel would then see the empirical gDNA
  alternatives directly.
  - It is the fallback if STAR cannot be parameterized to report them itself.
  - The STAR study prices it against a single STAR run.
- **Unannotated isoforms** (owner). Rescuing reads whose junctions reuse annotated splice sites would need implied
  transcripts, which is transcript assembly. That is outside Rigel's scope, so such reads are withheld and counted.

## Decisions

**Recorded (owner, 2026-10-15/16):**

- complete pairs, and a fragment length is required;
- plausible spliced and contiguous alignments both compete in the EM;
- unannotated spliced alignments are withheld;
- fix everything before release;
- go one wave at a time, and start no wave before its design is confident and reviewed.

**Recorded (owner, 2026-10-16, third round):**

- a multimapper's gDNA term is a sum over its places, not the mean;
- Rigel proposes no spliced reading for an unspliced read, because that is the aligner's job;
- mismatch differences are tolerated, and every hit is scored with its own mismatches;
- unannotated junctions are withheld completely;
- STAR's parameters are explored in a parallel session on the cluster (`star_parameters/PROMPT.md`).

**Recorded (owner, 2026-10-16, fourth round):**

- wave 1 is split: the scan stage goes first, without alignment scrutiny, and opens with a general code-health
  review;
- the genome store is dropped for now;
- Rigel's scope is narrowed: no alignment repair and no second aligner inside Rigel;
- Rigel documents and recommends STAR with the recommended parameters;
- a dual-aligner input is the eventual fallback.

**Recorded (owner, 2026-10-16, fifth round):**

- two pairing strategies: STAR by HI groups, and minimap2 by an R1 × R2 search made efficient (same reference, sorted
  by distance);
- a mate spliced at a junction beside a partner running into its intron is passed on as reported: the hit is spliced
  (RNA only), and the overhang penalty does not matter because no other transcript competes;
- the wave-1 plan (revision 2) is accepted; step 1, the code-health review, starts first;
- the roadmap is tracked to the completion of every wave and step (§Progress).

**Recorded (owner, 2026-10-16, sixth round):**

- the code-health review is taken whole: the cleanup stages in its order, starting with the scan-cache key; its
  scan-path defects in wave 1; its later-wave findings as ISSUES entries;
- the output column names must be consistent (a naming proposal first);
- the anchor tolerance stays in scope for implicitly spliced fragments, but every case that applies it is judged.

**Recorded (owner, 2026-10-17, seventh round):**

- output names: `count_<what>` word order; RNA is reported as annotated (every reference transcript, single-exon
  included) or nascent (Rigel's synthetic spans) — `DESIGN.md` §0; the gene table's mature / nascent split goes;
  the index-side names are renamed and every index rebuilt;
- the 1,000-read tag cutoff is repaired before 0.8.0 (the tag read per record, wave 1 step 2);
- K stays as the implicit-splicing tolerance, renamed;
- the report tiles, the summary's depth keys and rounding, and the resource constants as proposed;
- wave 1 is reviewed again, at revision 3, before step 2 is implemented.

**Recorded (owner, 2026-10-17, eighth round):**

- the STAR panel for wave 1's end-of-wave gate is the five conditions of the local minimap2 panel, at the production
  and the recommended parameters;
- splice-artifact behavior is measured on simulated data aligned with the recommended parameters, because its truth
  is exact;
- the owner has found STAR parameters that report both readings of a template whose alignments score alike; by
  default STAR's bonus for the spliced reading reports it alone;
- open questions or uncertainties after revision 4 go to a targeted re-review.

**Recorded (owner, 2026-10-18, ninth round):**

- rustar, the Rust port of STAR, is approved for local alignment; its residual difference from STAR (about 0.1% of
  templates on a 100k-pair sample, with the spliced-and-contiguous counts within one template per class) is accepted;
- the first measurement proceeds on the STAR panel aligned on the cluster.

**Recorded (owner, 2026-10-18, tenth round):**

- the recommended STAR parameters stay. They are more correct: they make ambiguous alignments visible instead of
  hiding them, and Rigel's job is to use that data. They are not evaluated against today's tool, which is known to
  be broken, so the first measurement's production-against-recommended comparison records today's tree and is no
  evidence about the parameters;
- Rigel is fixed in pipeline order, the scan, then calibration, then after calibration, one small change at a time
  (the fix order above);
- the caches are rebuilt once, when the scan's changes are done, and not after each one: the scan's changes are
  measured end to end, which reads no cache, and each one leaves the caches stale anyway;
- the production-parameter BAMs are deleted (their logs, junction tables and checks are kept).

**Recorded (owner, 2026-10-18, eleventh round):**

- the scan and buffer phase has two goals: keep competing spliced and unspliced alignments, so the aligner's choice no
  longer certifies a gDNA fragment as RNA; and a valid fragment length for every accepted alignment, so no unpaired
  reads. Every other scan defect is bundled with them, less critical;
- a gene-level count for templates that one gene holds but no single transcript is a future enhancement
  (`ISSUES: fragments-that-fit-one-gene-but-no-transcript`);
- whether STAR input is paired by its HI groups or by Rigel's own search, and whether that search is valid, is to be
  investigated (`ISSUES: star-hi-pairs-against-rigels-own-pairing`);
- stale or unneeded BAMs and results are removed once nothing needs them. Removed 2026-10-18, text receipts kept: the
  production-parameter BAMs, the parked threads' annotated BAMs and prototype wheels, the superseded pairing
  prototype's arrays and builds, the closed fl-tilt campaign's fragment dumps, and the parity and smoke BAMs.

**Recorded (owner, 2026-10-18, twelfth round):**

- the recommended STAR parameters are final: nothing runs on the production ones again, so the prior-denominator
  comparison uses the recommended inputs only;
- sequencing errors are shown with engineered cases, specific errors placed at specific junctions, not with whole
  error-injected copies of the panel;
- the richer-annotation simulation is not planned: the cost it would price is accepted by the annotated-only ruling;
- no header check of the aligner for HI: the per-group checks suffice.

**Open for wave 1:** the questions in `PAIRED_END_WAVE1_PLAN.md` §13.
