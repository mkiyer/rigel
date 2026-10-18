# Wave 1 plan: the scan stage

2026-10-18 · code base `01aba957` · revision 5. Revision 4 had a targeted re-review on 2026-10-18
(`PAIRED_END_WAVE1_TARGETED_REVIEW.md`): two blockers, two fixes before implementing, and a recommendation for each open
question. Revision 5 answers them (§14) and follows the owner's fix order of 2026-10-18 (`PAIRED_END_ROADMAP.md`): the
scan's changes one at a time, each measured end to end, and the panel caches rebuilt once, when the scan's changes are
done (§12). Two changes have landed: the junction-strand tag read per record (`01aba957`) and the transcript strands
installed by the index (`07eba0fc`). Revision 4 answered the readiness review of revision 3 (`PAIRED_END_WAVE1_READINESS_REVIEW.md`).

The history: revision 1 also scrutinized each pair's alignment and constructed alternative readings, using a genome
stored in the index. The owner then narrowed Rigel's scope (2026-10-16): Rigel reads alignments and does not repair or
re-align them. Revision 2 replaced revision 1 and was accepted on 2026-10-16. Revision 3 recorded step 1 carried out
(§5) and the owner's rulings of 2026-10-16 and 2026-10-17. The waves are in `PAIRED_END_ROADMAP.md`. Part A is the
prompt for the reviewer; Part B is the plan.

# Part A — Prompt for a check of revision 5 (copy everything below the line)

---

You reviewed revision 4 of this plan (`PAIRED_END_WAVE1_TARGETED_REVIEW.md`). Revision 5 answers your findings, and §14
maps each one to the section that answers it. Check only what changed:

1. **Your two blockers:** the one-witness rule (§7.2, §8) and the mixed-HI rule with its asserted pairs (§7.1). Does
   each falsifier you gave now hold, and does either rule admit anything it should refuse?
2. **Your two fixes:** the evidence gates and the normalization contract (§9.3, §11), and each record's status from its
   template's outcome (§10).
3. **The open questions (§13):** where revision 5 departs from your recommendation, is the reason sound?

**Deliver** at most one page: a verdict in at most three lines, then one line per finding of yours (answered, partly
or not) with the falsifier that decides it. Do not write production code.

---

# Part B — The plan

## 1. What wave 1 must achieve

For every template in the BAM, the scan stage produces either a set of validated hits, or one disposition saying why
there is none. A **hit** is one alignment of the template as the aligner reported it: both mates, checked to be one
molecule, with a computable fragment length.

**The two goals** (owner, 2026-10-18). The scan and buffer phase must:

- **keep competing spliced and unspliced alignments**, the potential splice artifacts that the recommended STAR
  parameters now report, so that a gDNA fragment is no longer certified as RNA by the aligner's choice alone;
- **give every accepted alignment a valid fragment length:** no more unpaired reads, because a length cannot be
  computed from one mate, and the EM scores by it.

Every other scan defect in this plan is bundled with them because the code is open anyway; each is less critical, and
the arms take the two goals first (§12). What the length rule loses, and one way to win part of it back, is
`ISSUES: fragments-that-fit-one-gene-but-no-transcript` (a future enhancement).

Wave 1 delivers six things:

1. a code-health review of the scan path, and its cleanup, before anything else changes (carried out, §5), and proof
   that the gates which judge the rest can fail (§12);
2. one mate location per alignment;
3. hits from the aligner's pairs, or from R1 × R2 pairing where the aligner leaves pairing to the reader, each
   validated;
4. the annotated-path rule, applied hit by hit;
5. the hand-off to calibration and the EM, with four corrections to calibration's second pass;
6. one disposition per template, with every record written once.

Wave 1 does **not**:

- construct, edit or re-align any alignment: no bases are moved, no reading is proposed, and no genome is stored
  (owner, 2026-10-16);
- change how a template's several hits combine (wave 2), although it measures how many templates that rule touches
  (§9.2);
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

**Calibration** estimates the library's gDNA and RNA content before the EM, in two passes. Pass 1 deposits each template
whose placement and path it can decide, and holds the rest: templates left with several placements, and templates whose
gap between the mates admits more than one path (`accumulator.cpp:1073-1163`). Pass 2 assigns each held template to one
of its placements by a seeded draw, weighted by pass-1 traffic. The owner treats this as a deliberate approximation to a
circular problem: the EM needs calibration, and calibration would need the EM.

**The bug that started this.** The scanner built "fragments" out of single mates. A lone secondary, a primary whose
mate was unmapped, or a half-filled HI group each became a read-length fragment, a false multimapper and a 100-bp
length deposit. Fixing it exposed two downstream defects that had compensated for it (the roadmap, wave 3).

**The scope (owner, 2026-10-16).** A gDNA molecule aligned as spliced is certified as RNA, because gDNA cannot
splice. Its bases usually fit a contiguous alignment almost as well, but by default STAR gives the spliced alignment a
bonus and reports it alone.

- Rigel will not construct the other alignment itself.
- The remedy is the aligner's. The owner has found STAR parameters under which STAR reports both alignments when they
  score alike (2026-10-17), and Rigel will recommend them. Rigel then scores what is reported.
- Wave 1 is what Rigel must get right whichever parameters produced the BAM.

## 3. Principles

1. **Read, don't repair.** A hit is an alignment the aligner reported, checked. Rigel never moves bases, never forms an
   alignment the aligner did not report, and never re-aligns. Pairing mate alignments that an aligner reports per mate
   (§7.1) is reading, not repair. The splice blacklist is the one declared exception until wave 4 (§8).
2. **Both mates, and a length.** A hit is usable only if it places both mates and its fragment length can be computed.
3. **The emission model is unchanged.**
4. **The aligner proposes, the EM decides.** The scan stage never chooses between hits.
5. **One disposition per template; every record written once.**
6. **Clean before changing.** Every cleanup is proven to change no number, and every behavior change is its own
   measured arm.
7. **No new constant.**
8. **Today's behavior is cited, never recalled.** Every statement of what the code does today names the file and line
   it was read from; every claim that a rule is exact or complete is checked by enumeration.

## 4. The scan stage in one picture

```
step 1  the code-health review and its cleanup as proven no-ops; then proof that the gates can fail

the BAM records of one template
 │
 ├─ step 2  one mate location per alignment (CIGAR, clips, skips, NM, motif); chains withheld
 ├─ step 3  pairs → hits: HI groups (STAR), or every R1 × R2 pair (no HI); each validated as one molecule
 │            one reference · opposite strands · inward · some path with 1 ≤ L ≤ max_frag_length
 ├─ step 4  a spliced hit is usable only if one annotated transcript holds all its junctions
 ├─ step 5  usable hits → the EM; offered placements → calibration; num_hits = distinct usable hits
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
| the output names | `1fa1c143` | outputs report RNA as annotated or nascent (owner's convention); `summary.json` schema 4 | measured on one ladder condition (g50, stranded, capture-ON), every output file compared: each renamed field keeps its values bit for bit, and only the defined relabellings move |
| the public repository's scrub | `8401aeba` | no local path, cluster detail or lab identifier in a tracked file; the data root is `$RIGEL_RUNS` | bit-identical |
| index and scan names | `dc736e39` | index format 10; `--nascent-tolerance`, `--intron-overhang-tolerance`; the contiguous reaches no longer computed | bit-identical; every index rebuilt and checked file by file |
| the truth tables | `5c67c1f4` | the simulator's one truth writer, each quantity once | an A/B gives identical reads |
| constants, part 2 | `0c754afc` | the resource and display numbers given a home and a reason | bit-identical; the simulator's reads identical, its BAM's `@PG` version now the package's |
| tile-invariant fits | `4bcc72d9` | the row-tiled calibration fits give one answer whatever their tile height | measured on the ladder, per stratum: 0.998–1.001 of the tree before (0.996–1.007 per condition) |

The suite went from 2,959 to 3,051 cases: gates added for each refusal and definition, cases deleted with the code
they tested. Every new native gate was broken in a rebuild and seen to fire. Every panel was re-cached twice (262
conditions re-certified each time: once for the payload's schema, once for index format 10), and the references were
re-frozen at each measured stage, each new freeze proven equal to the one before but for the measured change.
Preflight now checks each cache's key, not only its presence, and that check found the oracle defect.

**Moved out of step 1:** the tag detector's 1,000-read cutoff is repaired in step 2 (the junction-strand tag read per
record); and two renames wait for the next re-simulation, because every panel on disk carries the old names and must
stay readable: the synthetic spans' `RIGEL_NRNA_` ids and the truth columns `pre_capture_*` → `molar_*`.

**The readiness review's qualifications** (2026-10-17), accepted:

- the identity references cover one ladder condition, its minimap2 counterpart and VCaP. Re-certifying 262 caches
  proves that they are accepted, not that 262 conditions are unchanged;
- the tile fix is a measured numerical change, not a no-op, so every wave-1 baseline is recorded on the tree after it;
- the gates themselves must be shown to fail before the fixtures rely on them: step 1b, the gates' own proof (§12);
- no second general cleanup is needed.

## 6. Step 2 — records into mate locations

Each eligible, non-supplementary record becomes one **mate location**. It holds:

- its mate (R1 or R2), reference and strand;
- the full CIGAR as a query-to-reference map;
- soft and hard clips at each end, kept apart;
- the raw CIGAR-N skips (before the blacklist) and the kept ones;
- NM, as the record reports it;
- the motif strand (XS or ts, below);
- its aliases: every record and HI value that carries it.

No bases are kept, because nothing in wave 1 re-scores a sequence.

**Rules:**

- **One location per alignment.** A location is identified by its mate, reference, strand, leftmost aligned position
  and CIGAR, with M, = and X read as one operation, adjacent runs of them coalesced (`50=50X` and `100M` are one
  alignment), and soft and hard clips kept apart. So `50M3D50M` and `50M3N50M` at
  one start are two locations, a deletion and an intron, although they place the same aligned bases. Records with one
  identity are one location: STAR repeats a mate's record in every pair alignment that shares it, and each repeat is
  an alias. The identity is a property of each record alone, so record order cannot change it.
- **Conflicts are input defects.** Aliases of one location that disagree on NM are counted, and the location keeps the
  larger NM. Aliases that disagree on the motif make it ambiguous.
- **No M/=/X base, no location.**
- **Withheld, not absent.** A location that depends on a supplementary chain, or has a hard-clipped 5′ end, is
  withheld: counted, and not paired. A mate whose only locations are withheld is "unresolved", not "missing".
- **NM is required.** A record without NM is counted and warned about; its NM still reads 0, as today
  (`bam_scanner.cpp:662-664`). That is a disclosed approximation: it can favor a hit of unknown NM over one whose NM is
  known. STAR's default attributes omit NM, so the recommended command asks for it. STAR's `nM` is never read in its
  place: it counts the pair's mismatches, while NM counts one mate's edits, indels included.
- **Widths.** `num_hits` and NM are stored wide enough for any count the scan can produce, or refused explicitly when
  they would not fit. Today both are stored in 16 bits through unchecked casts, so 65,536 hits read as 0 and 65,537 as
  1 (`resolve_context.h:207,249,268`).
- **Unindexed references keep their own identity.** Two mates on different references outside the index are a trans
  pair, not one reference. Today every unindexed reference shares the id −1 (`bam_scanner.cpp:160-173`), so
  a pair on two unindexed references is merged into one fragment by coordinate (`bam_scanner.cpp:711-771`;
  `unindexed-references-share-id-minus-one`).
- **Tag types are checked.** NH and HI must be integers. Another type is an input defect, and is counted. A missing NH
  is counted too (§7.1).
- **The junction-strand tag is read per record** (owner, 2026-10-17: before 0.8.0). Today one mode (XS, ts, both or
  none) is chosen per BAM from its first 1,000 spliced records, a literal default (`bam_scanner.cpp:2877-2930`,
  resolved once per run at `pipeline.py:987`), so a library whose first spliced reads lack the tag
  loses its strand (`ISSUES: sj-strand-tag-chosen-from-the-first-reads`). Revision 4 reads each record's tags on their
  own:
  - XS is genomic. minimap2's ts is relative to the read, so it is flipped on a reverse record, as today
    (`bam_scanner.cpp:610-613`);
  - both present and agreeing give that strand; disagreeing, the motif is ambiguous and the record is counted. Today,
    with both tags, XS is read and ts only when XS is absent, so a conflict passes unseen (`bam_scanner.cpp:626-629`);
  - one present gives that strand; none, on a spliced record, is counted;
  - a value other than `+` or `−` is read as absent, and counted. Today an invalid XS also hides ts
    (`bam_scanner.cpp:594-613`);
  - an explicit `--sj-strand-tag` restricts which tags are read, and is recorded.

  There is no "mostly missing" threshold. The counts and fractions of missing and conflicting tags are reported, and the
  existing strand-model diagnostic warns when strand evidence is unavailable. This is its own arm (§12): it moves
  numbers only on libraries whose tags are absent from their first spliced reads, mixed, or conflicting.

**The index sets up the resolver.** It installs the transcript strands when it loads, and K and the maximum fragment
length come from the scan configuration. Today only the scan installs them (`resolve_context.h:613-615`, called from
`pipeline.py:296-300`), so every unit test runs without strands, and without them `ambig_strand` is always 0
(`resolve_context.h:1287-1303`; `resolver-strands-installed-only-by-the-scan`).

## 7. Step 3 — hits from the aligner's pairs

### 7.1 Which pairs: two strategies (owner, 2026-10-16)

**The mode is read per template**, from its own records, not from the BAM's first records. Whether STAR input
should be paired by HI at all, or by Rigel's own search with HI ignored, is open (owner, 2026-10-18,
`ISSUES: star-hi-pairs-against-rigels-own-pairing`): it is measured on the STAR panel, against each molecule's truth,
before the HI-group arm is built. The rules below hold if HI is kept.

- **With HI (STAR).** Every record of the template carries HI. STAR groups the records of each pair alignment by HI. A
  hit starts from the R1 and R2 locations of one HI group, and Rigel forms no other combination: crossing HI groups
  would only duplicate reported hits, or revive pairs STAR did not report.
- **Without HI (minimap2).** No record carries HI. minimap2 reports secondaries per mate, so pairing is the reader's
  job (below).
- **Mixed.** A template whose eligible records mix the two keeps STAR's pairs and invents none (targeted review,
  finding 2). It is read with HI: each HI group's asserted pair stands, and a record without HI is an input defect,
  counted and paired with nothing. Revision 4 paired a mixed template without HI, which crossed known groups: with
  forward locations at 100 and 120 under HI 1 and 2, and reverse ones at 250 under HI 1 and 270 untagged, it accepted
  all four lengths (250, 270, 230, 250), forward HI 2 with reverse HI 1 among them. Today a record without HI joins
  group 0 of its template (`bam_scanner.cpp:821-846`).

**What decides the mode.** The template's eligible records only: the mapped records the scan reads after its QC-fail
and duplicate filters. An excluded record never sets or switches the mode, so adding an unmapped untagged record to a
template changes nothing. HI is read as STAR documents it, the index of one pair alignment, carried by both of its
mates' records, and the completeness checks below test that reading template by template. Rigel does not name
producers: STAR and its port rustar write different program lines (`STAR`, `rustar-aligner`) with the same semantics,
and a list of trusted names would be one more thing to keep current. An aligner whose HI meant something else would
fail the per-group checks, counted, not paired across groups.

**An asserted pair** is an HI group with exactly one distinct location for each mate, after aliases collapse. A group
with two distinct locations for one mate is ambiguous: counted, and it yields no hit, because STAR asserts one pair per
group.

**Completeness.** Under STAR's semantics, a template's HI values run from 1 to NH, and every group holds both mates.

- **Truncated output:** the template's distinct HI values number fewer than NH. The run stops, with the remedy
  `--outSAMmultNmax -1`. No remedy recovers what STAR filtered before output: a template beyond
  `--outFilterMultimapNmax` is reported unmapped.
- **A group missing one mate**, checked group by group, yields no hit. Its record is an unpaired alignment, counted as
  incomplete mate reporting. The template-wide count above cannot see this: R1 with HI 1 and 2, R2 with HI 1 only, and
  NH 2, has two HI values against NH 2, and its group 2 holds R1 alone.
- A missing NH, or an NH that differs between a template's records, means completeness is unknown: counted, not
  refused. Today the scan takes the first NH above 1 it meets (`bam_scanner.cpp:1804-1807`).

**Duplicate hits.** Two groups whose R1 and R2 locations are identical are one hit, which keeps both HI values as
aliases. `num_hits` counts distinct hits.

**Without HI: every pair, the bound first.** Each R1 location is paired with each R2 location on the same reference
with the opposite strand, and every such pair goes to the validator (§7.2). The validator's first length test is the
span bound

    G ≤ max_frag_length + S1 + S2 + Tmax

where G is the two mates' aligned envelope, S1 and S2 their observed skip totals, and Tmax the largest total intron
length of any annotated transcript on that reference. A pair beyond it is rejected for length, which it would fail
along every path of §7.2:

- L counts the path bases in [a, b), and the removed bases inside [a, b) are at most S1 + S2 + Tmax, so
  G ≤ b − a ≤ L + S1 + S2 + Tmax;
- checked by enumeration on tiny references with clips, several introns and K from 0 to 3: 48,301 accepted pairs in
  200,000 random cases, with no violation; removing any one term of the bound produces violations (Appendix B).

The work is O(n1 · n2) per template, where n1 and n2 count the two mates' locations, with no cap. Under minimap2's
default of five secondaries it is at most 36 pairs. Revision 3 walked the partners by distance and stopped at the first
one beyond its own bound. That misses a farther partner whose longer skip raises its bound: the readiness review's
counterexample has a near partner at span 2,100 against a bound of 1,000, and a farther spliced partner at 3,100
against 3,850 that validates at L 250.

The search decides by geometry alone. Today, without HI, a unique mapper's primaries form one hit unconditionally
(`bam_scanner.cpp:892-899`). A multimapper's locations are paired first wherever their resolved transcript sets
intersect, then nearest first, then every remaining cross pair, and what is left becomes one-mate hits
(`bam_scanner.cpp:907-1043`). Every location inside one gene shares the gene's synthetic span, so every combination
there passes (`no-hi-pairing-by-transcript-overlap`).

Pairing mate alignments that the aligner reported per mate is reading, not repair.

### 7.2 One molecule

A candidate pair is a hit when two conditions hold.

1. **The mates are on one reference, with opposite strands.** Call the forward mate F and the reverse mate R; either
   may be R1.
2. **They can be one molecule along some path.**
   - The path P is the genome with every observed raw skip of either mate removed. An alternative P also removes one
     annotated transcript's introns that lie in a stretch neither mate aligned, between the mates or beyond a mate's 5′
     end where its clip is walked, by the rule the resolver applies (§7.5): an intron that overlaps the stretch and lies
     within K bases of its ends, on each side, is cut whole (`resolve_context.h:816-842`), and a transcript with an
     intron overlapping an aligned block by more than K bases offers no path (`resolve_context.h:682-704`). Aligned
     bases inside a cut intron are off the path, and the overhang term prices them.
   - **One witness** (targeted review, finding 1). The transcript whose introns an alternative P removes must hold
     every junction the hit observes, the junctions §8 checks, so the observed chain and the implied introns come
     from one transcript, clip walks included. For an unspliced hit any transcript qualifies. The base path, the
     observed skips alone, stands; the blacklist exception and the ordinary overhang are unchanged; and the witness
     limits neither the candidates nor their lengths. Revision 4 let the two come from different transcripts: with F
     `10M100N10M` at 100, R `10M` at 1700, `max_frag_length` 1,000 and K 0, T1 holding the observed intron [110,210)
     and T2 the intron [300,1300), the base path is 1,510 and fails, but T2's intron gave 510 while §8 passed on T1.
     Now the pair is rejected for length, and a transcript holding both introns admits it at 510.
   - Walking F's 5′ clip back along P gives the molecule's start a. Walking R's 5′ clip forward gives its end b. A walk
     that leaves the reference is unresolved geometry.
   - Require a < b, both mates' aligned bases inside [a, b), and 1 ≤ L ≤ `max_frag_length`, where L is the number of
     path bases between a and b.
   - The path decides acceptance only. It never removes a candidate or changes a candidate's length.

**One K rule with the resolver.** A strict reading of "wholly inside a gap" would reject pairs the resolver scores.
With intron [1200,2100), F [1103,1203), R [2150,2250) and K 3, the strict path gives L 1,147, beyond the default
1,000, while the resolver's rule removes the intron, keeps F's three intronic bases as overhang and gives L 247. The
validator uses the resolver's rule (checked, Appendix B: K 3 and 4 accept at 247, K 2 rejects).

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
- The penalty cancels only when every surviving explanation of the template carries it. A clean alternative elsewhere
  changes the odds, and a hard gate (overhang α = 0) removes the hit, which is then counted as gated.

**The owner agrees (2026-10-16).** One mate is spliced, so the hit is spliced: RNA only, with no gDNA candidate. The
intron overhang is penalized, but no other transcript is expected to compete for the template, so the penalty does not
matter. Revision 1 constructed both consistent readings instead; that is repair, so it is out, and no fallback is
planned.

**Mates whose junctions contradict.** Both mates spliced at overlapping but unequal junctions disagree about a donor or
an acceptor. No transcript holds two overlapping introns, so no path explains both, and the hit is withheld by §8 with
its own reason. Today calibration measures such a fragment along the union of the two introns, a path neither mate
reported (`accumulator.cpp:425-464`).

### 7.4 What a hit carries

- **Contents.** Each hit carries both locations, its raw and kept skips, its aliases, and NM summed over the two mates.
- **Kind.** Hits of one template whose spans overlap on one reference are **same-place alternatives**: for example, the
  aligner's spliced and contiguous alignments of one molecule. The kind is a reporting class only. Overlapping spans do
  not prove one physical origin, so no rule treats such hits as one (§9.2).
- **`num_hits`** is the number of distinct usable hits after step 4. It replaces today's `max(NH, hits)`
  (`bam_scanner.cpp:1857`). The multimap option (`include_multimap=false`) is applied after it is counted (§10).

### 7.5 The anchor tolerance (owner, 2026-10-16)

`--intron-overhang-tolerance` (K, default 3; `--splicing-anchor-tolerance` until 2026-10-17) belongs to implicitly
spliced fragments: mates aligned to different exons with a whole intron between them, where up to K aligned bases may
sit inside the intron. The owner keeps it there, "but we must be judicious about deciding which cases apply it".

**A use is classified by the gap it is applied to, not by the hit's splice label.** A hit with an observed junction can
also have an unsequenced gap between its mates, and an annotated intron in that gap is an implicit splice. So the
question at each use is whether it tests an implied intron in an unaligned gap.

Today K is read in one place, inside the gap enumeration, which runs for every non-chimeric fragment of any splice
type, SPLICED_ANNOT included (`resolve_context.h:780`, `1320-1322`;
`ISSUES: the-anchor-tolerance-applies-beyond-implicit-splicing`). It matters only where two blocks leave a gap that is
not exactly an observed intron (`resolve_context.h:736-784`), and from there it reaches three decisions:

- `transcript_contradicts_blocks` (`resolve_context.h:682-704`, called at `:790`): a transcript with an intron
  overlapping an aligned block by more than K bases contributes no hypothesis;
- `collect_transcript_introns_in_gap` (`resolve_context.h:816-842`, called at `:796`): an annotated intron that overlaps
  the gap and lies within K bases of its ends is cut whole;
- the number of hypotheses then decides the SPLICE_IMPLICIT label (`resolve_context.h:1334-1336`).

K is not used for fragment lengths, chimera detection, junction matching or scoring.

Wave 1 keeps the uses that test an implied intron in a gap, and justifies or removes each other one in its own measured
arm. Each use is tested at 0, K and K + 1, once with an implied gap and once beside an observed splice. No use may shift
an observed junction, delete alignment evidence or a candidate, or forgive two observed skips that conflict. The mates
that disagree in §7.3 are not a use of K: the overhang term prices those bases, as the owner ruled.

**Renamed** (owner, 2026-10-17): K is a valid tolerance for implicit splicing, and the renames stage gave it a name
that says what it is, `intron_overhang_tolerance` (`--intron-overhang-tolerance`, `dc736e39`).

## 8. Step 4 — the annotated-path rule, hit by hit (owner, 2026-10-16)

**Rule.** A spliced hit is usable only if **one annotated transcript holds every junction the hit observes**, both
mates' together. Otherwise it is not usable:

- **not as RNA:** no annotated transcript has its path, so its length on a transcript cannot be computed;
- **not as gDNA:** it is spliced.

The hit is withheld, with one of three reasons:

- **a novel junction:** some junction is in no annotated transcript;
- **an unsupported combination:** each junction is annotated, but no one transcript holds them all;
- **contradictory mates:** two observed junctions overlap without being equal (§7.3).

A junction matches an annotated one as the resolver matches it today: exactly, at the strand of a definite motif, and at
either strand when the motif is absent or ambiguous (`resolve_context.h:959-993`). The rule decides only whether a hit
is usable. Its candidates are still nominated by overlap, and the transcript that holds its junctions limits neither the
candidates nor their lengths. It is the transcript whose introns may fill the hit's unaligned stretches when §7.2
validates it (one witness).

**The template's other hits stay.** Take a gDNA template that the aligner reports both as an unannotated splice and as
a contiguous alignment: it keeps the contiguous hit.

**Templates left with nothing.** A template with no usable hit left, and with a hit withheld here, takes the outcome
**unannotated splice** (§10). It is reported in three classes, read off the genes whose span contains each side of the
junction:

- within one gene;
- between genes (read-through);
- touching an intergenic region.

**No alternatives.** Rigel constructs no contiguous or partial reading of an unannotated junction. If the aligner
reports one, it is a hit like any other.

**Not in wave 1:** counting a withheld template that one gene holds, but no single transcript, to that gene alone,
outside the EM (owner, 2026-10-18; `ISSUES: fragments-that-fit-one-gene-but-no-transcript`). The census of withheld
templates above, by class, is what it would be designed from.

**What changes from today.** Today such hits are quantified through overlap candidates, whose lengths count the novel
junction's skipped bases as if they were sequenced: a candidate's length subtracts only its own annotated introns
(`resolve_context.h:900-954`). A fragment with one annotated and one novel junction is labelled annotated-spliced
(`resolve_context.h:1188-1202`). A combination of annotated junctions that no transcript joins is scored like a
supported path, with overhang: the readiness review's same-gene fixture gives TJ1 150/50/500 and TJ2 150/50/250, as an
annotated splice. The cost of withholding them is an explicit policy cost:

- **on VCaP**, an earlier screen counted 114,519 RNA-library junction units, about 0.8% of spliced RNA units, against
  2,253 DNA-library artifacts. That is exposure by junction unit, not genuine RNA lost: it is recounted by validated
  template, under the three reasons;
- **on the minimap2 panel**, 90% of gDNA artifacts are unannotated long-anchor junctions.

Rigel does not assemble transcripts. Reads of unannotated isoforms, even ones that reuse annotated splice sites,
contribute to no transcript or gene. The simulator writes only annotated isoforms, so their cost on genuine RNA is
priced on VCaP's RNA library (§12).

**The blacklist: the one declared exception.** The blacklist (`alignable`) is unchanged in wave 1, as ruled. It removes
a junction before step 4 sees it: a record's N whose shorter anchor is within the junction's blacklisted anchor length
is erased from its junction list, which leaves an unaligned hole between two blocks
(`bam_scanner.cpp:567-583, 689-693`). With no other junction left, the fragment is labelled an artifact splice
(`resolve_context.h:1253-1257`), and calibration never receives it (`bam_scanner.cpp:1675`). So the rule reads the kept
junctions, and a hit whose only junction is blacklisted is a contiguous hit. That deletes alignment evidence. It is a
repair, and the one exception to "read, don't repair", until wave 4's matched comparison retires it. Until then a census
records every hit whose raw skips include a blacklisted junction, annotated or novel, with its disposition.

## 9. Step 5 — usable hits into calibration and the EM

### 9.1 Usable hits, offered placements, held templates

Three sets, told apart:

- **Usable hits** reach the EM. Each is one EM hit; a template with more than one is a multimapper unit.
- **Offered placements** reach calibration: every hit of the template, whatever its candidates, except an
  artifact-spliced hit (`bam_scanner.cpp:1675`) and a hit on an unindexed reference (`bam_scanner.cpp:1679-1683`).
  Today the scanner offers a multimapper's intergenic hit as a placement (`bam_scanner.cpp:1926-1929`), while the EM
  never sees it: only a unique mapper's intergenic fragment enters the buffer (`bam_scanner.cpp:1915-1920`). NH still
  counts it (`ISSUES: multimapper-intergenic-alignments`).
- **Held templates** are the ones calibration's pass 1 holds. Offering several placements does not by itself hold a
  template (`accumulator.cpp:1073-1163`):
  - identical placements collapse to one (`accumulator.cpp:1109-1115`);
  - one placement left takes the one-placement path, which deposits it, or holds it when the gap between its mates
    admits more than one path (`accumulator.cpp:1121-1126`, `549-553`);
  - with several left, every (placement, hypothesis) within the length limit survives; exactly one survivor is
    deposited, and otherwise the template is held whole (`accumulator.cpp:1128-1162`).

**No consumer changes in wave 1.** The EM and calibration keep these routes, including calibration's intergenic
placements, until wave 2 sets the rule for intergenic hits.

**For STAR input with complete pairs,** `num_hits` then equals NH, except where step 3 rejects a hit or collapses a
duplicate, or step 4 withholds one.

### 9.2 The EM's combination: kept, and measured

Wave 1 changes none of the EM's combination rules; wave 2 does (owner, 2026-10-16). Revision 3 said the existing rules
were right for same-place alternatives. In general they are not:

- **The RNA merge.** For each transcript, the EM keeps the hit with the least overhang, then the highest likelihood
  (`scoring.cpp:576-589`). Overhang is already inside each hit's likelihood, so the rule can keep a worse term:
  - (overhang, NM) = (0, 2) beats (1, 0): it keeps weight 0.01 instead of 0.1 (`merge-keeps-least-overhang`);
  - at equal NM, lengths 200 and 201 with overhang 0 and 1, where f(200)/f(201) = 10^-4, it keeps a term 1,000-fold
    worse (the readiness review's probe).

  The two rules agree when NM, length and strand agree, as in §11's same-place toy, where TA1 takes the spliced hit,
  with overhang 0, under either rule.
- **The gDNA term** is the mean of the gDNA-capable hits' likelihoods, not their sum (`scoring.cpp:654-661`). With one
  gDNA-capable hit it is that hit's own term, which is correct; in a spliced-and-contiguous pair only the contiguous hit
  is gDNA-capable. With k gDNA-capable hits at different origins it is 1/k of the sum the owner ruled for wave 2, not at
  most a factor 2, as revision 3 said.
- **The prior's denominator.** Calibration deposits each held template at its drawn placement before the EM prior is
  read, but `n_calibrated` counts unique fragments only (`pipeline.py:700-720`). In a locus holding multimappers the
  gDNA share is inflated, and a locus of multimappers alone gets no prior. The owner kept this denominator for 0.8.0
  on a measured A/B (`ISSUES: the-em-prior-share-excludes-drawn-multimappers`). Revision 3 said the two stages agree;
  they do not.
- **Pruning** compares candidates against three different references (`pruning-uses-three-references`).

**Their exposure grows with wave 1's input.** The recommended STAR parameters report a contiguous reading beside many
spliced ones, so more templates become two-hit units. The first measurement and every arm count it (§12). Open
questions 1 and 2 (§13) ask what would move a fix into wave 1.

### 9.3 Calibration's second pass: four corrections

Pass 2 draws one (placement, hypothesis) of each held template by a seeded draw (`second_pass.py:745`). Today the
draw's weight is, within the record (`second_pass.py:170-240`):

- the placement's traffic A_p, the unique pass-1 traffic at the objects the placement would deposit on, when the
  record has more than one placement (`second_pass.py:540-542`);
- the hypothesis's density ρ, in fragments per base: at the junction for a spliced path, at the exon–intron
  boundaries for the contiguous one, the scarcest object deciding (`second_pass.py:549, 567-593`);
- the length likelihood f(L), and the strand term s.

A factor that is zero for every surviving candidate is dropped as uninformative; a partial zero is a hard zero.

**The all-held approximation (accepted by the owner, and measured).** If every template that crosses a junction is
held, the junction's pass-1 traffic is zero. A spliced placement then gets weight 0 when its contiguous alternative has
positive traffic, and an even 1/2 against 1/2 when neither has any and the other factors tie. A contiguous draw is a
deposit in the unspliced tally, not a gDNA assignment.

The owner expects all-held junctions to be rare, because templates that map uniquely across a junction establish its
traffic. The recommended STAR parameters hold more templates, so it is measured on the STAR panel (§12):

- all-held, partly held, independently supported and unsupported junctions, across depths;
- the predicted probabilities and the actual draws, with truth;
- what support was available, and the RNA length and strand training;
- the deposits and calibration against the oracle, the priors, and the genuine RNA finally lost;
- matched seeds and further seeds.

A count of zero weights alone is not enough.

**Four corrections, each its own arm.** Each also changes today's multimapper draws.

1. **The length law.** A placement with an admissible annotated junction, observed or implied, is RNA and takes the RNA
   law. Today the law is chosen from the hypothesis's implied introns only (`second_pass.py:596`), so an observed
   spliced placement takes the global law.
2. **The weights.** Each candidate carries exactly one placement factor and one share within its placement:

   ```
   w(p, h) = A_p × ρ(h) / Σ ρ(h′ of p) × f(L) × s
   ```

   - the share is dimensionless, so the draw no longer depends on ρ's units;
   - a placement with one hypothesis has share 1, so a placement with no contested region is no longer zeroed against
     one that has (today its ρ reads 0, `second_pass.py:593`);
   - all-zero ρ within a placement gives an even share within it, and a partial zero stays a hard zero; the existing
     uninformative rule still applies to A_p, f and s;
   - on a one-placement record it equals today's weight;
   - the share is taken over **unique hypotheses**: two descriptions of one physical path are one hypothesis, merged
     before the share and kept merged through the deferred bank. Without it, duplicating one hypothesis moved a
     placement from 0.354839 to 0.411765 in the targeted review's probe. It is a contract on the input, not a claim
     that today's resolver emits duplicates;
   - the weights are formed in logs and normalized against the largest, so a small factor raised to a large power
     cannot underflow: 0.1^400 against 0.1^401 is odds of 10:1, where the direct products both read 0. A hard gate is an
     explicit −∞, kept apart from the zero-evidence fallback (an even share), which applies only where a factor is
     uninformative for every candidate.

   Revision 3 proposed ρ = 1, "not applicable", for such a placement. That gives some candidates two densities and
   others one, so a tenfold change of units moves an even pair to 1/11 against 10/11. The checks (Appendix B): 20,000
   random records unchanged under any rescaling of ρ and A_p; one-placement records equal to today; an uncontested
   placement never zeroed beside contested ones; and the raw-density break caught. The targeted review found them
   blind to the evidence itself: a uniform draw over every multi-placement record, with today's one-placement weight,
   also passes all four. So the gate adds independently computed cases (§11): two single-hypothesis placements with A_p
   1 and 3 and every other factor equal draw 1/4 and 3/4; a duplicated hypothesis keeps the physical paths'
   probabilities; a common NM offset keeps the odds.
3. **The sequence evidence.** Pass 2 never reads mismatches: the deferred bank carries none
   (`bam_scanner.cpp:1142-1158`), so a placement with 37 more mismatches is drawn as often as a perfect one.
   - Each held placement carries its NM into the deferred bank: a new field in the bank's schema, its spill format, the
     cache, and every copy `lift_choices` makes.
   - Pass 2 multiplies the weight by the EM's own mismatch term, α_mm^NM, with the configured α_mm, once.
   - α_mm = 0 is a hard gate, as in the EM: a placement with NM > 0 gets weight 0, and a record whose every placement
     is gated is not deposited, and is counted. α_mm = 1 is inert. At a tie the factor is equal and changes nothing.
4. **The motif rule** (`second-pass-restates-the-motif-rule`). Pass 2 restates the deposit's rule two other ways
   (`second_pass.py:282-304, 604-609`):
   - it passes an AMBIGUOUS motif into its sj lookup, which matches it at either strand, while the deposit trusts no
     splice on an AMBIGUOUS motif;
   - it decides that a placement is RNA from the motif tag, while the deposit and A_p decide it from an observed intron
     that resolves to an annotated sj. Upstream the motif is set only from an observed intron that matched an
     annotated junction (`resolve_context.h:1165-1167`), so an annotated splice whose tag is missing, or contradicts
     the annotation, gets gDNA's ½ in pass 2 while the deposit books it as spliced.

   Each placement's annotated-junction use and its strand are derived once, by the deposit's rule, and the length,
   placement and strand factors all read them. Fixtures cover a missing tag, a tag that contradicts the annotation, and
   two opposite-strand annotations sharing coordinates. It moves numbers only on rare held placements of tagged
   stranded libraries. The ladder cannot see it, because the simulator writes XS on every spliced read.

The first, second and fourth are local to `second_pass.py`. The third also touches the accumulator's deferred bank.

### 9.4 One length rule

The resolver gives a lone merged block its own length for every candidate (`resolve_context.h:935-940`). That skips the
subtraction of the candidate's introns inside the block, which `tx_frag_length` and DESIGN §3.4 apply
(`single-block-length-shortcut`). So one molecule gets two lengths, depending on whether its mates overlapped. The
shortcut is removed in its own arm; most of the candidates it affects already sit deep in overhang.

## 10. Step 6 — accounting for every template

**Outcomes.** Each template has exactly one, decided in this order:

1. **excluded:** no eligible record;
2. **excluded by the multimap option:** more than one usable hit under `include_multimap=false`, counted after step 4;
3. **placed:** at least one usable hit;
4. **unannotated splice:** none usable, and some hit was a valid molecule that §8 withheld, with its class (within one
   gene, read-through, touching an intergenic region). A chimera label on that hit is kept as a secondary reason,
   because a read-through is a withheld unannotated splice (owner, 2026-10-16);
5. **reported chimera:** none usable, none withheld by §8, and a candidate pair resolves as chimeric;
6. **unpaired:** none of the above, with the furthest reason any candidate pair reached: missing mate; unresolved
   location; geometry; unresolved geometry; length.

Every reason that any pair reached is recorded beside the outcome.

**Each record's status** follows from its template's final outcome (targeted review, finding 4):

- **used:** in a usable hit of a placed template;
- **unused:** in a placed template, but in no usable hit (a withheld hit's records, a lone secondary);
- **unplaced:** every record of a template with any other outcome, an excluded one included.

So `records = used + unused + unplaced` holds template by template. Revision 4 tied the status to usable hits, so a
template excluded by the multimap option left a record in none of the three: two usable HI groups and an orphan third,
under `include_multimap=false`, came out one short. The usable-hit counts taken before the option are reported apart,
and neither an option exclusion nor a hit §8 withheld reaches calibration or the EM.

**A placed template records its route:**

- whether it reached the EM with candidates;
- whether it is intergenic only;
- in calibration, whether it was deposited or held, and its pass-2 draw;
- its hits by kind: same-place alternatives, multi-place;
- its records that were not used.

**The counts reconcile, by identity as well as by total:**

> templates = excluded + placed + unannotated splice + reported chimera + unpaired
>
> records = used + unused + unplaced

The annotated BAM's records are matched to the scan's by template and record, not only counted.

**Counters.**

- Today's `mate_unmapped` actually counts one-mate hit groups (`bam_scanner.cpp:1824-1828`). It is renamed
  `one_mate_groups`, and a raw count of records flagged mate-unmapped sits beside it.
- New counters: `unpaired_alignments`, `incomplete_mate_reporting`, `duplicate_hits`, `mixed_pairing_mode`,
  `missing_nh`, `nm_conflicts`, `tag_conflicts`, `withheld_chains`, and the blacklist census of §8.
- Rates are reported against all templates and against eligible templates.

**The annotated BAM** writes every original record once, with only the intended tag changes. Tags carry the template's
outcome and reason and the record's status (`ZC:unpaired` with the reason, `ZC:unused_alignment`), and stale tags are
cleared.

**Fixes the code-health review found** (2026-10-16), each with its fixture:

- `proper_pairs` counts primary records that carry the proper-pair flag, so it matches its record-level label
  (`proper-pairs-counts-hits`). Today it counts two-mate hits by their first primary R1's flag
  (`bam_scanner.cpp:1829-1837`).
- A gated fragment is annotated as gated, with its reason, not as intergenic, and the writer refuses duplicate
  frag_ids (`gated-fragments-stamped-intergenic`; today `bam_scanner.cpp:2551-2566`).
- The intergenic append is deleted, so an intergenic pair on two references is no longer counted and annotated
  chimeric (`intergenic-append-is-dead-and-mislabels`).
- The annotation pass's groups must equal the scan's, by identity (each group's read name and record count) and not
  only by count, or the run stops (`annotation-pass-can-drift-a-group`). Today nothing compares the two passes
  (`annotate.py:441-449`), and a real drift exists: a filtered record between two records of one name splits the group
  in the annotation pass only, so every later frag_id is off by one (`bam_scanner.cpp:1494-1512` against
  `2730-2737`).
- A multimapper's ZS ranks all five splice types (`multimapper-splice-type-drops-implicit`).
- Every scan counter is asserted exactly, on one hand-built BAM, at one worker and at four (`scan-counters-untested`).

## 11. Tests

One suite per step. Each fixture is one of two kinds:

- **invariant:** passes on the baseline and after;
- **defect demonstration:** fails on the baseline, passes after the fix, and fails again when the fix is broken on
  purpose.

**The exhaustive pairing oracle.** A small, independent Python enumeration of every R1 × R2 pair of a template, along
every path of §7.2, stands beside the suites. On tiny references with clips, several introns and both orientations,
the scanner's hits and every rejection reason must equal the oracle's: the complete accepted and rejected sets, with
the clip and K cases among them, not lengths the validator returns. The span-bound check of Appendix B is necessary,
not complete: with every annotated path removed from the validator it still prints `ALL CHECKS PASS`, with 47,255
pairs accepted.

**The toy reference.** These are the annotation's exons, 0-based and half-open:

| transcript | strand | exons |
|---|---|---|
| TA1 | + | [1000,1200) [1500,1700) [2000,2300) |
| TA2 | + | [1000,1200) [2000,2300) |
| TB1 | + | [5000,5150) |
| TC1 | − | [2200,2400) [2600,2800) |

**Expected hits.** The emission values come from the baseline resolver (`counterexample_baseline.py`), and the
readiness review recomputed them independently. In the table, "x/y/z" means exon bases / overhang / fragment length
for that candidate. Cases marked † use their own tiny reference, named in the case.

| suite | case | records (F / R) | expected |
|---|---|---|---|
| locations | a deletion against an intron | `50M3D50M` and `50M3N50M` at one start, each with the same partner | two locations and two hits |
| locations | one pair under two HI values | HI 1 and HI 2 carry identical records | one hit with both HI values as aliases; `num_hits` 1 |
| pairing | same-place alternatives, as STAR would report them | HI 1: `97M300N3M` at 1103 / [1550,1650); HI 2: [1103,1203) / [1550,1650) | two hits, one unit. HI 1 gives TA1 200/0/247. HI 2 gives TA1 197/3/247, TA2 97/103/547 and unspliced A 200/0/547, with gDNA at 547 |
| pairing | mates disagree, acceptor side | F `50M800N50M` at 1150 / R [1998,2098) | one hit, as reported: TA2 148/2/148 |
| pairing | mates disagree, donor side | F [1103,1203) / R `50M300N50M` at 1150 | one hit, as reported: TA1 147/3/147 |
| pairing | clips hiding TA2's junction | F `50S50M` at 2000 / R `50M50S` at 1150 | placed along TA2's path, L 100 |
| pairing | outward, unclipped | F [1300,1400) / R [1100,1200) | unpaired: geometry |
| pairing | a lone secondary, no partner in its HI group | one mate only | an unpaired alignment, counted; no hit |
| pairing | one mate's record missing from an HI group | R1 HI 1 and 2, R2 HI 1 only, NH 2 | HI 1 is the hit; HI 2 is an unpaired alignment, counted as incomplete mate reporting |
| pairing | a primary whose mate is unmapped | one mate only | unpaired: missing mate |
| pairing | two HI groups at different loci | both asserted | exactly the two asserted hits; no cross pair |
| pairing | per-mate secondaries, no HI, far copies | R1 at A and B, R2 at A′ and B′, with A and B farther apart than the span bound | two hits, (A, A′) and (B, B′); the cross pairs fail the bound (length) |
| pairing | per-mate secondaries, no HI, near copies | 100-bp R1 forward at 100 and 120, R2 reverse at 250 and 270 | four hits, lengths 250, 270, 230 and 250; all kept |
| pairing | no HI, a partner inside the bound but beyond `max_frag_length` | the long-intron witness ([1050,1150) / [2100,2200)), its mates reported per mate | paired: the genomic span is 1,150 against the default 1,000, but L is 350 along TA2 and 550 along TA1 |
| pairing | † the variable bound (no annotation) | reverse R1 [10000,10100); forward partners `100M` at 8000 and `50M2850N50M` at 7000; `max_frag_length` 1,000 | the near partner fails the bound (span 2,100 against 1,000); the far one passes (3,100 against 3,850), validates at L 250, then §8 withholds it (novel junction) |
| pairing | no HI, a partner on another reference or on the same strand | R1 and R2 per mate | rejected for geometry; counted |
| pairing | NH above the number of HI groups | truncated output | the run stops, with the remedy |
| pairing | the tolerance K | intron [1200,2100) of a † one-transcript reference; F [1103,1203), R [2150,2250); `max_frag_length` 1,000 | K 3: one hit, L 247, overhang 3. K 2: rejected for length (1,147 along the genome) |
| annotated path | a novel junction inside gene A, with a contiguous alternative | HI 1: junction [1250,1500); HI 2: contiguous | HI 1 withheld; HI 2 kept; the template is placed, unique |
| annotated path | a junction between genes A and B, alone | `50M2800N50M` at 2200 / [5100,5200) | withheld: unannotated splice, read-through; the chimera label is a secondary reason |
| annotated path | the same junction, with a valid contiguous alternative | HI 1 as above; HI 2 contiguous | HI 2 placed; HI 1's reasons recorded |
| annotated path | mates whose junctions contradict | F `50M300N50M` at 1150 / R `50M800N50M` at 1150 (TA1's and TA2's overlapping junctions) | withheld: contradictory mates; calibration never measures their union |
| annotated path | † two annotated junctions that no transcript joins | one gene, TJ1 holding one junction and TJ2 the other | withheld: unsupported combination (today TJ1 150/50/500, TJ2 150/50/250) |
| annotated path | a known and a novel junction | one mate at TA1's [1200,1500), the other at a novel junction | withheld (novel) |
| annotated path | blacklisted junctions | one hit whose only junction is blacklisted; one with one of two junctions blacklisted | the kept junctions decide; the census records the raw ones |
| pairing | one-base overlap, both outer boundaries, long-intron witness | as in revision 2 of the first design, §12 | as there |
| pairing | † a hybrid witness (finding 1) | F `10M100N10M` at 100 / R `10M` at 1700; T1 holds [110,210), T2 [300,1300); `max_frag_length` 1,000; K 0 | rejected for length (1,510); with a T3 holding both introns, one hit at L 510; TA2's clip case still L 100 and the K case still 247 |
| pairing | a mixed-HI template (finding 2) | 100-bp forward at 100 (HI 1) and 120 (HI 2); reverse at 250 (HI 1) and 270 (no HI) | one hit, (100, 250); HI 2 an unpaired alignment; the untagged record a counted defect. The same with an excluded untagged record added, and with the records permuted |
| pairing | an HI group with two distinct locations for one mate | HI 1: R1 at 100 and at 400, R2 at 250 | ambiguous: counted; no hit |
| accounting | an excluded template's records (finding 4) | two usable HI groups and an orphan third; `include_multimap=false` | one outcome; five disjoint record statuses, all unplaced; conserved through the annotation and under permutation |

**Other cases:**

- **Locations:**
  - duplicates within and across HI groups, and an alias reused by several hits;
  - the same records in another order: identical locations and hits;
  - `5S95M` against `95M5S`;
  - conflicting NM or motif;
  - missing NM, and a BAM without NM: counted and warned;
  - hard-clipped carriers; chains without SA; adjacent H and S;
  - a record with no M/=/X base;
  - mates on two unindexed references;
  - NH or HI typed as a string; NH missing; NH conflicting between a template's records;
  - malformed mate flags; two templates sharing a read name;
  - CIGAR D, I, =, X and H, before and after an N (`no-indel-cigar-fixtures`);
  - 65,535, 65,536 and 65,537 hits, through persistence: exact counts or an explicit refusal.
- **The junction-strand tag:**
  - XS `+` with ts `−` on a reverse record: `+`; XS `+` with ts `+` on a reverse record: ambiguous, counted;
  - 1,000 untagged spliced records, then a tagged one: the tagged record's strand is read;
  - an explicit `--sj-strand-tag`, honored and recorded;
  - an invalid tag value: counted.
- **Pairing and the resolver:**
  - a template whose first records lack HI and whose later ones carry it; a template mixing the two;
  - without HI, two copies inside one gene's introns: only the geometric pairs;
  - an antisense overlap with the strands installed by the index: `ambig_strand` is 1;
  - each use of the anchor tolerance (§7.5) at 0, K and K + 1, once with an implied gap and once beside an observed
    splice;
  - a merged block over an intron: the candidate's path length (§9.4).
- **Accounting:**
  - one fixture per reason, and one per outcome in the order of §10;
  - a geometry failure beside a too-long consistent pair;
  - a group whose every hit the multimap option excludes;
  - every original record written once, with only the intended tag changes; stale tags cleared;
  - the counts reconciling on a mixed BAM, by identity as well as by total;
  - the counters BAM, every counter exact;
  - a gated fragment's annotation;
  - an annotation pass out of step with the scan, including equal-sized groups in another order: refused.
- **Calibration:**
  - the genic and intergenic pair (HI 1 [1010,1150), HI 2 [6000,6150)): one EM hit and two calibration placements,
    through a fresh scan, the cache and the drain, traced into the priors;
  - an all-held junction; a partly held one; one supported only by held templates; one supported by nothing;
  - the length-law correction;
  - the weights: every density scaled tenfold leaves every probability unchanged; a duplicated hypothesis description
    changes nothing; a placement with no contested region beside contested ones is not zeroed; two single-hypothesis
    placements with A_p 1 and 3, every other factor equal, draw 1/4 and 3/4; factors of 0.1^400 against 0.1^401 give
    odds of 10:1; a common NM offset keeps the odds;
  - placements that differ only in NM, 37 apart: odds 10^-37; α_mm 0 with every placement carrying NM: not deposited,
    counted; α_mm 1: inert;
  - NM surviving the deferred bank's spill and reload, the cache, and `lift_choices` across partitions;
  - the motif rule: an AMBIGUOUS motif, a missing tag, a tag that contradicts the annotation, and two opposite-strand
    annotations sharing coordinates.
- **Code health:**
  - each cleanup stage passes `rename_identity.py --check`;
  - a deliberately moved output value in a stage makes `--check` fire.

**Break tests.** Each must fail its suite:

- restore the one-mate hit;
- restore `max(NH, hits)`;
- skip pair validation;
- form a pair across two HI groups;
- take only the primary pair on a BAM without HI;
- stop the no-HI search at the first partner beyond its own bound (the variable bound must fail);
- bound the no-HI search by `max_frag_length` alone (the long-intron witness must fail);
- identify a location by its aligned bases alone (D against N must fail);
- count a pair reported under two HI values twice;
- choose the pairing mode from the first records of the BAM;
- check HI completeness over the template instead of per mate;
- store `num_hits` in 16 bits;
- bound candidates by exon ends;
- prune by `L_t`;
- walk clips contiguously only;
- fail containment on a mate's bases inside its partner's skip;
- validate with the strict gap rule instead of the resolver's K rule (the tolerance case must fail);
- withhold the whole template when only one of its hits is unannotated;
- withhold a hit whose junctions one transcript holds;
- check junctions one by one instead of together;
- merge contradictory mates' junctions into their union;
- let the chimera label outrank a withheld unannotated splice;
- drop intergenic hits from `num_hits`;
- read XS alone when both tags are present;
- drop NM from the deferred bank on spill;
- weight an uncontested placement by its raw density in pass 2;
- let a cleanup stage move a number;
- read a missing NM as 0 without counting it;
- pair without HI by transcript sets;
- file unindexed references under one id;
- leave the strands to the scan;
- stamp a gated fragment intergenic;
- accept an annotation pass whose groups differ, by count or by identity;
- decide "RNA" in pass 2 from the motif tag;
- keep the single-block length shortcut;
- borrow an implied intron from a transcript that does not hold the observed junctions (the hybrid witness must fail);
- remove every annotated path from the validator (the oracle's accepted and rejected sets must differ);
- pair a mixed-HI template without HI, and let an excluded record set the mode;
- draw every multi-placement record uniformly (the 1/4 against 3/4 case must fail);
- reverse each weight factor in turn, and make every NM gate fire;
- form the weights as direct products (the 0.1^400 case must fail);
- tie a record's status to usable hits rather than to its template's outcome (the excluded template must fail).

## 12. Measurement and gates

**Inputs:**

- **The STAR panel** (owner, 2026-10-17), the decisive input. It holds the five ladder conditions of the local minimap2
  panel:
  - `gdna_g00_ss_0.99_nrna_mid_capture_off`, the zero-gDNA control;
  - `gdna_g05_ss_0.50_nrna_mid_capture_off` and `gdna_g50_ss_0.50_nrna_mid_capture_off`;
  - `gdna_g50_ss_0.99_nrna_mid_capture_off` and `gdna_g50_ss_0.99_nrna_mid_capture_on`.

  They cover the three in-scope strata and compare condition for condition with minimap2. Each is aligned at the
  production parameters and at the recommended ones, under which STAR reports a near-tied contiguous reading beside
  the spliced one (the owner's result, 2026-10-17). Splice artifacts are measured here, on simulated reads, because
  their truth is exact (owner, 2026-10-17).
- **Sequencing errors (proposed; open question 3).** The ladder's reads are error-free, so the mechanism by which a
  sequencing error near a junction tips STAR to the spliced reading is absent from them. The proposal aligns one seeded
  copy of each condition beside them, with substitutions at the VCaP DNA library's mismatch rate (the STAR study's
  injection instrument). The names, and so the truth, are kept, and the two copies hold the same templates.
- **The local minimap2 panel** (five conditions, simulator truth). Its per-mate secondaries are paired by the search
  (§7.1).
- **The oracle ladder**, as the identity check. Its alignments are perfect, so the location and pairing arms must leave
  it unchanged apart from enumerated cases.
- **VCaP** (STAR, production parameters): a test input only. It covers transitions, library provenance, and the costs
  that simulation cannot show: the annotated-path rule's cost on genuine RNA, because the simulator writes only
  annotated isoforms, and artifacts from sample variants or library chemistry.

**The first measurement, before any arm** (owner, 2026-10-17; taken 2026-10-18, `PAIRED_END_ROADMAP.md`). Today's
tree on the STAR panel, production against recommended parameters, per stratum, scored against truth. By the owner's
ruling of 2026-10-18 it records today's tree and is no evidence about the parameters, which stay. It also takes the
census of the deferred defects (§9.2), still to do:

- templates by number of usable hits, and by kind;
- units in which the merger keeps a hit other than the most likely one, for some transcript;
- units with two or more gDNA-capable hits;
- per locus, the held templates deposited against `n_calibrated`;
- spliced placements given zero weight in pass 2, with their truth.

**Arms,** one mechanism each, in the owner's fix order of 2026-10-18 (`PAIRED_END_ROADMAP.md`), the two goals of §1
first:

1. the junction-strand tag read per record (§6): landed, `01aba957`, a numeric no-op on the STAR panel and the ladder;
2. the transcript strands installed by the index: landed, `07eba0fc`, a numeric no-op on both panels, as it is by
   construction (the scan installed them before every scan);
3. the pairing question, HI's pairs against Rigel's own search (`ISSUES: star-hi-pairs-against-rigels-own-pairing`):
   a measurement, which changes no number of the tool;
4. the gates' own proof (below), before the first arm whose identity claim leans on those gates;
5. location identity and aliases, which duplicate hits need (§6);
6. complete pairs only: HI groups, with duplicate hits, completeness per mate and the mixed-template rule (§7.1), or
   the search, as the pairing question decides;
7. `num_hits` counting distinct usable hits (the one-mate fix);
8. every pair validated without HI, the bound first;
9. the path witness for clipped and outward pairs, with one witness transcript (§7.2);
10. the annotated-path rule (§8), which withholds an unannotated spliced hit and keeps its contiguous competitor;
11. the single length rule (§9.4);

then the bundle, each its own arm:

12. the rest of §6: the NM count, the unindexed references, the tag types and the widths;
13. chains withheld;
14. each use of the anchor tolerance outside implicit splicing, one arm per use (§7.5);
15. the accounting fixes (tags and counters; no quantification number moves);
16. the length-law correction;
17. the weights in pass 2;
18. mismatches in pass 2;
19. the motif rule in pass 2.

**Every arm:**

- measures end to end on fresh scans: quant_accuracy's base arm reads each BAM and no cache, because the cache's key
  covers the deposit rule, not pairing or resolution. The panel caches are rebuilt once, when the scan's changes are
  done, before calibration's instruments read them (owner, 2026-10-18);
- records its inputs' content digests, the paths of the modules it imported, and the native build's hash;
- is compared with its parent and with the common baseline, recorded on the tree after the tile fix (§5);
- runs pinned (`scan.total_threads=1`, `OMP_NUM_THREADS=1`, fractional assignment), with pass-2 draws at matched seeds
  and at further seeds, so that sampling variation is told apart from the effect;
- is read per stratum, never pooled.

**Readings of each arm:**

- **the setup:** the gDNA law and the ruler frozen, with a hash of the arrays each consumer received, and also each
  arm's own; per stratum; scored against truth;
- **the errors:**
  - both errors, by truth origin: artifact RNA mass, and genuine RNA lost;
  - transcript and gene error;
  - calibration against the oracle;
  - the pools, including the unspliced-RNA over-call;
- **the census:** a per-template transition census with truth; the hits counted by kind; the deferred defects' census
  above;
- **calibration:** the second pass's predicted probabilities and actual draws, with truth;
- **cost:** scan time and memory.

**The gates' own proof** (step 1b; before any fixture relies on them). Each gate is broken on purpose and seen to fire:

- `stage_gate.sh` and every arm's gate exit non-zero when any check fails (today they record a failure and exit 0);
- the output comparison checks record counts, and each record's position, CIGAR and sequence, not only names, flags
  and tags (today it pairs records with `zip`, so a tail goes unseen);
- `rename_identity.py` checks fields by name where names are unchanged, the nested banks and the scalar QC, and captures
  the EM under fractional assignment;
- the mutants: a stale cache after a pairing-only change; an input changed under the same size and time; a failed
  suite; a BAM tail appended or deleted; a moved coordinate; two arrays swapped; a changed nested NM.

**The implementation gate:**

- this plan reviewed and revised: revision 5's check of the targeted review's findings, before the arms they govern
  (6, 8, 9, 15 and 17);
- the gates' own proof done;
- every suite's fixtures written, failing where they demonstrate a defect, with the exhaustive pairing oracle;
- the injection hooks proven to reach their consumers: the frozen calibration and pass-2 laws, the EM's realized laws,
  the priors where declared, and the transcript and gDNA rulers.

**The end-of-wave gate:**

- every fixture passes, and every deliberate reversal fails;
- the full suite passes at the re-derived baseline, and every annotated record is conserved;
- oracle identity holds, with every template-level exception enumerated;
- on the STAR panel, per stratum, with the law and the ruler frozen and also with each arm's own, no in-scope stratum
  gets worse, except by the measured, accepted policy costs (unpaired templates, withheld unannotated hits,
  contradictory mates);
- the readings above are reported for every arm;
- the owner reassesses before wave 2.

## 13. Open questions

The targeted review recommended an answer to each of revision 4's six (its finding 5). Revision 5 adopts them, except
where the owner must decide:

1. **The deferred defects' exposure** (§9.2). Adopted: the wave boundaries stand unless a separate causal arm shows that
   a deferred mechanism stops a stratum from meeting its gate; counts alone promote nothing. Likelihood losses, gDNA
   multiplicity and prior distortion are recorded beside the truth errors, and the validated-hit census is taken with
   the oracle, since today's `num_hits` is not it.
2. **The prior's denominator** (`ISSUES: the-em-prior-share-excludes-drawn-multimappers`). Adopted in substance: the
   three-denominator comparison is re-run before the wave completes, all-multimapper loci included, and the
   denominator is not replaced merely because its cohort mismatch is real. It runs on the recommended STAR inputs only:
   the parameters are final and nothing runs on the production ones again (owner, 2026-10-18).
3. **Sequencing errors** (§12). Not whole error-injected copies of the panel (owner, 2026-10-18): random errors over
   millions of templates spend a lot of compute to reach the handful that tip STAR at a junction. Instead, engineered
   cases: templates built at specific junctions with specific errors placed where they tip the aligner (a substitution
   in a short anchor, a mismatch just inside an exon end, an indel near a donor), aligned at the recommended parameters,
   to show that the near-tied contiguous reading is reported beside the spliced one. Needed for the end-of-wave gate.
4. **The all-held approximation** (§9.3). Adopted: reopened only if matched support, depth and seed experiments show a
   loss from zero support that stops an in-scope gate; changing support from zero to independently supported must move
   the probabilities and deposits as predicted. Exposure alone decides nothing.
5. **What simulation cannot show** (§12). Adopted for VCaP: withholding in its RNA library is reported as exposure, not
   as genuine RNA lost, since a withheld alignment can be genuine unannotated RNA or misaligned contamination, and VCaP
   corroborates. The review's instrument for the cost on genuine RNA, known-origin RNA simulated from a richer annotation
   and quantified against the normal one, is not planned: the cost it would price is accepted by the annotated-only
   ruling (owner, 2026-10-16), and VCaP's count bounds it (owner, 2026-10-18).
6. **The no-HI search's form** (§7.1). Adopted, as revision 4 had it: every pair validated, the bound applied to each
   pair; multiplicity tails and memory profiled, and optimized only if a cost asks, keeping exact hits and reasons.
   Joined by the owner's question whether STAR input is paired by HI at all (`ISSUES: star-hi-pairs-against-rigels-own-pairing`).

## 14. How this answers the reviews

**The targeted review of revision 4** (2026-10-18):

| finding | revision 5 |
|---|---|
| 1 (blocker). the observed and implied junctions can have different witnesses | §7.2: one witness, the transcript that holds the observed junctions supplying every implied intron, clip walks included; §8; the hybrid fixture and its break (§11) |
| 2 (blocker). the mixed-HI fallback discards known pairs | §7.1: a mixed template keeps its asserted pairs and invents none; the mode from eligible records only; one distinct location per mate in an asserted group; ambiguous groups counted. Producer naming declined, because the per-group checks test HI's reading (the owner agrees, 2026-10-18). Whether HI is kept at all is the owner's open question (§7.1) |
| 3. the evidence gates and the normalization contract | §9.3: unique hypotheses, weights in logs, hard gates apart from the zero-evidence fallback; §11: the oracle's complete accepted and rejected sets, the unequal-factor, duplicate and offset cases, and their mutants |
| 4. an excluded template leaves records without a status | §10: each status from the final outcome; the five-record fixture (§11) |
| 5. exposure, error and unidentifiable cost | §13: each recommendation adopted, or put to the owner |

The review judged revision 4's answers to the readiness review: 1, 3, 5 and 6 answered (1 with M/=/X runs to coalesce,
now §6); 2, 4, 7 and 8 partly, now completed by its findings 1, 3, 2 and 4.

**The readiness review of revision 3** (2026-10-17):

| finding | revision 4 |
|---|---|
| 1. location identity collapses D and N | §6: identity is the alignment itself, aliases kept; duplicate hits collapse (§7.1); order invariance tested |
| 2. junctions one by one are not an annotated path | §8: one transcript must hold every observed junction; contradictory mates withheld; the blacklist declared as an exception, with its census |
| 3. the search's stopping rule | §7.1: every pair validated, the span bound first; checked by enumeration, each term needed |
| 4. pass 2's normalization | §9.3: one placement factor and a share within the placement; checked by enumeration; zero and hard-zero behavior stated; NM through every path |
| 5. consumer descriptions and closure claims | §9.1–§9.2 restated from the source with citations; the deferred defects named with their ISSUES entries and measured (§12); this section says what is deferred |
| 6. cache validity and proof | §12: fresh scans, content digests and build hashes; the gates' own proof |
| 7. input and numeric contracts | §6, §7.1: the mode per template, completeness per mate, a missing NH, the widths, `nM` never read as NM |
| 8. outcomes and the tolerance | §10: the outcome order with every reason kept, record status, the multimap option after counting, groups checked by identity; §7.2, §7.5: K applied by gap, one rule with the resolver |

**The earlier reviews' blockers**, as the readiness review judged them, with revision 4's status:

| blocker | status |
|---|---|
| first design 1: best-hit RNA, mean gDNA | not answered: deferred to wave 2 (§9.2), with its exposure measured in wave 1 |
| first design 2: preference is not identification; one owner of the sequence evidence | removed: the aligner proposes and owns NM; a missing NM is counted and disclosed (§6) |
| first design 3: calibration can zero a held RNA placement | an accepted approximation, not resolved: measured, and the weights' units fixed (§9.3) |
| first design 4: the clip walk | design answered (§7.2; L 100 along TA2); record-level fixtures pending (§11) |
| wave-1 review 1: the local merger | not answered: wave 2 (§9.2); exposure measured |
| wave-1 review 2: the prior's denominator | not answered: kept for 0.8.0 by ruling; exposure measured; open question 2 |
| wave-1 review 3: the mismatch proposal cutoff | removed: Rigel proposes nothing; the aligner's recall is the STAR study's |
| wave-1 review 4: validation before reconstruction | removed: each reported alignment is validated as it stands |
| wave-1 review 5: incomplete reconstruction | removed |
| wave-1 review 6: unsupported paths; the blacklist bypass | answered by §8's whole-path rule; the blacklist is a declared exception until wave 4 |
| wave-1 review 7: sequence for recomputed NM | removed: nothing is recomputed; a missing NM is a disclosed approximation |

## Appendix A — what left wave 1, and where it went

- **Alignment scrutiny** (moving bases, Δm, the generation tolerance, reconciled readings, multi-junction enumeration):
  removed. Near-tied alternatives are the aligner's job. The STAR study (`star_parameters/PROMPT.md`) establishes the
  parameters Rigel will recommend.
- **The genome stored in the index:** dropped for now (owner, 2026-10-16). Nothing in the plan reads bases.
- **The splice-artifact machinery of 2026-10-14/15** (the local tie table and its second hit, the hold rule, the
  mismatch emission factor): superseded by the aligner's alternatives and their NM.
- **R1 × R2 enumeration:** kept only for aligners without HI (§7.1); STAR input is paired by HI.
- **A dual-aligner input:** a spliced aligner and a genome aligner, merged by read name into one BAM for Rigel. This is
  the owner's eventual fallback if STAR cannot be parameterized; the STAR study prices it.
- **Remote copies** (the index-time copy table): wave 2, and only if VCaP shows that STAR misses them.
- **The merger rule, the gDNA sum and the prior's frame:** wave 2 and 0.8.1, by ruling; wave 1 measures their
  exposure (§9.2, §12).

## Appendix B — evidence pointers

- **The targeted review of revision 4:** `PAIRED_END_WAVE1_TARGETED_REVIEW.md`, with its reproduction record (lab use
  only).
- **The readiness review of revision 3:** `PAIRED_END_WAVE1_READINESS_REVIEW.md`, with its reproduction record (lab
  use only).
- **Revision 4's checks:** `$RIGEL_RUNS/arms/2026-10-17_wave1_rev4_checks/`: the span bound (200,000 random cases,
  48,301 accepted pairs, no violation; each term's removal produces violations), the second pass's units (four checks,
  and the raw-density break caught), and the K example. Each script prints its receipt beside it.
- **The code-health review:** `SCAN_CODE_HEALTH_REVIEW.md` and `scan_code_health_slices/`; the cleanup is
  `41899ae4..3941a0dd`.
- **Diagnosis tables:** revision 2 of the first design, Appendix A, with `paired_end_repair_evidence/arms.json`.
- **The mate-disagreement census:** `census_*.json`.
  - End-extensions within 3 bases are 495 mRNA against 58 artifact templates on minimap2 OFF.
  - On VCaP they are 2,165 RNA-library against 5 DNA-library templates.
  - Recount after pair validation.
- **Emission values on the toy annotation:** `counterexample_baseline.py`, confirmed independently by the first
  review and again by the readiness review.
