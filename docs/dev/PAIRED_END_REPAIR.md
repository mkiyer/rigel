# Paired-end repair: the design, revision 2

2026-10-16 · base `41899ae4` · **superseded for planning** by `PAIRED_END_ROADMAP.md` and `PAIRED_END_WAVE1_PLAN.md`
(after the owner's 2026-10-16 decisions). Kept as the record of the response to the first review and of the diagnosis.

This revision incorporates the independent review (`PAIRED_END_REPAIR_INDEPENDENT_REVIEW.md`); §11 maps every finding
to its change. It replaces revision 1 of this file. Where they conflict, it supersedes `MATE_PAIRING_FIX.md`,
`MATE_PAIRING_IMPLEMENTATION.md` and the withdrawn continuation. The measurements are in
`paired_end_repair_evidence/` and the reviewer's checks in `paired_end_repair_review_evidence/`.

## 0. The design on one page

Rigel reads a name-collated BAM and must turn each read pair into evidence about one molecule. The repair has three
waves.

- **Wave 1: the scan stage.** It runs in six steps, and it adds no competing readings yet.
  1. **The reference genome in the index**, so the scanner can check bases. (Already approved, 2026-10-14.)
  2. **Read the records:** one *location* per alignment record, carrying its bases, clips, mismatch count and
     junctions.
  3. **Pair them:** every R1 location against every R2 location, through one molecule validator. A pair becomes a
     *placement* only if some path makes it one molecule of a possible length.
  4. **Account for every template:** each one is placed, a reported chimera, unpaired (with a reason) or excluded,
     and every record is written once.
  5. **Reconcile mates that disagree:** a mate overhanging into its mate's intron by up to K bases is corrected
     (moved across the junction). Every other disagreement is set aside, counted, until wave 2.
  6. **Classify read-through chimeras:** a junction whose two sides belong to different genes is a reported chimera,
     not quantified, whatever its span.
- **Wave 2: competing readings.** One pair can often be read two ways: spliced, or contiguous because the junction is
  an artifact. Three likelihood rules must be settled before any such alternative reaches the EM or calibration:
  - how readings and placements add up;
  - who owns the sequence evidence;
  - how calibration treats a reading whose support is entirely held.

  Only then are the alternatives switched on.
- **Wave 3: the downstream defects** the diagnosis found, which run independently of wave 2:
  - the gDNA fragment-length estimator;
  - false capture in the ruler;
  - the clip length model;
  - the rest of the release list (§9).

**What never changes:** Rigel's emission model. Any transcript overlapping the reads by even one base is a
candidate, and overhang pays its existing penalty. A placement's transcript candidates and lengths are computed
exactly as today. The repair decides which placements exist. It changes a placement's emission only where step 5
deliberately corrects an alignment.

## 1. Why: what went wrong, and what the review corrected

Rigel's scanner built "fragments" from a single mate: a lone secondary alignment, a primary whose mate was unmapped,
or an HI group holding one mate. Each became a read-length fragment, a false multimapper, and a 100-bp
fragment-length deposit.

Fixing that exposed two downstream defects that had compensated for it:

- **The gDNA length estimator is ill-conditioned on aligned input.** Its two-pool contrast divides by a purity
  difference that comes out near 0.014 on aligned input, so it amplifies small pool differences about 50–75 fold.
  The baseline's 26% spike at 100 bp happened to give the right mean. Without the spike the fit has mean 318 bp,
  against 217 on the oracle BAM of the same reads. That law moves calibration's capture efficiencies and hence the
  **capture ruler**: the EM's per-transcript, capture-adjusted effective length. The ruler moves reads between
  isoforms.
- **The ruler finds capture on uncaptured aligned libraries.** At zero gDNA, flattening it removes 93.8% of the
  transcript error. At the g50 capture-OFF conditions, about 2% of transcripts get a ruler factor off by more than
  2×. The oracle BAM shows 0.3–0.4%.

The evidence comes from diagnostic swaps: one input at a time (the gDNA law, or the ruler) taken from another run.
Swapping the baseline's ruler into stage 1 turns capture-ON 801k transcript error into 488k. With the law fixed,
each pairing stage improves capture ON and both capture-OFF strata.

**The review's corrections, accepted:**

- These swaps show the ruler and the law **mediate** the regression. They do not show that nothing else contributes.
- The law swap replaced only calibration's uniform-frame law. The EM's realized law and the second pass's fit were
  not touched.
- The capture-OFF ruler has outliers; it is not "flat".
- The zero-gDNA residual is *consistent* with the unpaired loss, not shown to be it.
- No arm crossed law and ruler together. §8 adds that experiment and requires every injection to hash the array its
  consumer actually received.

The tables are in Appendix A.

## 2. The words used here

| word | meaning |
|---|---|
| template | one read name: the two mates of one sequenced molecule and every alignment record of them |
| location | one alignment record of one mate, normalized (§4) |
| placement | one R1 location with one R2 location, accepted as one molecule (§5) |
| path | the genomic bases a molecule occupies between its ends, after removing introns |
| reading | one consistent interpretation of a placement's path; wave 2 lets one placement carry several |
| reconciliation | making two mates that disagree about a junction agree on one path |
| K | `--splicing-anchor-tolerance`, default 3 bp |
| the length law | a fitted fragment-length distribution; calibration fits one for gDNA and one for RNA |
| the ruler | the EM's per-transcript capture-adjusted effective length |
| disposition | what happened to a template, and to each of its placements |

## 3. The rules this design must keep (owner)

1. Three populations: gDNA, RNA+ and RNA−. Unpaired and chimeric are dispositions, not populations.
2. The emission model is preserved. One base of overlap nominates a transcript, overhang pays its per-base penalty,
   and no molecule must fit inside a transcript's exons. No candidate is deleted to make a path fit, and no target
   lends its length to another.
3. Every hypothesis places both mates. No single-mate likelihood and no mate rescue.
4. No new constant. K is the existing, approved tolerance.
5. Read-through chimeras are reported, not quantified.
6. A mate overhanging into its mate's intron is an aligner error to correct, neither penalized nor dropped.
7. Fix everything before release. The scan stage comes first; each fix has its own tests.

## 4. Step 1 and 2 of wave 1: the reference, then the records

### The reference genome in the index (first)

The index stores the genome it was built from: 2-bit packed, with N runs and other ambiguity codes in a separate
mask, mmapped at quant time. It also stores each contig's checksum, checked against the BAM header's `M5` where
present. A missing store stops the run with a rebuild instruction.

It comes first because everything that checks bases needs it: the mismatch counts below, reconciliation and wave 2.
Doing a correction first and adding its check later would create a confounded intermediate state (review).

### A location per record

Each eligible, non-supplementary record becomes one location holding:

- the full CIGAR as a query-to-reference mapping, not just its blocks, so that `5S95M` and `95M5S` at one start
  stay different;
- terminal soft and hard clips at each end, kept apart;
- the raw CIGAR-N skips, before the blacklist removes any, and the kept ones;
- the read's bases in query orientation. A secondary that omits them (minimap2 writes `*`) takes them from a full
  record of the same mate, reverse-complemented when the strands differ and sliced by original query position.
  Bases a donor lacks (its hard-clipped part) stay marked missing;
- **NM:** the record's value while the alignment is unchanged; recomputed against the stored reference when it is
  missing, conflicting, or the alignment is edited;
- motif strand, with provenance and record ids.

**Identity.**

- Two records are one location when they map every aligned query base to the same reference base, with the same
  clips. An insertion shifted inside a homopolymer is the same location only after left-normalization against the
  bases and the reference.
- Duplicates with different NM are resolved by recomputation, not by taking the smaller value. Opposing motif
  strands make the location's motif ambiguous.

**Withheld, not absent.** A location that needs a supplementary chain, or has a hard-clipped 5′ end, is *withheld*:
it exists, it is counted, and it is not paired (the approved chain policy, in its own arm). A withheld location is
never confused with a missing one, so "unresolved location" is a different reason from "missing mate".

**Reporting completeness.** Checked on raw records before any filter, and only under verified tag semantics
(STAR's NH/HI). A missing HI group, per mate or per template, stops the run with the read name and the remedy
(`--outSAMmultNmax -1`). Missing tags mean completeness is unknown, which is recorded.

## 5. Step 3: pairing

Every R1 location is tried with every R2 location. HI and the primary flag are provenance; they never forbid a
combination. Accepted pairs are deduplicated and ordered by geometry, never by record order.

**The validator.** A pair is a placement when all of the following hold:

1. **One reference, opposite orientations.** Call the forward mate F and the reverse mate R; either may be R1.
2. **The mates agree on one path** after step 5's reconciliation.
3. **The path witness.** Some path P makes the pair one molecule of a possible length:
   - P is either the **observed path** (both mates' bases, with every raw skip removed), or that path with **one
     annotated transcript's introns** removed where they lie wholly inside a gap neither mate aligned. Those
     introns get no K slack, so no aligned base is ever erased.
   - Walk F's 5′ clip back along P from F's first aligned base to get the molecule's start a, and R's 5′ clip
     forward along P to get its end b.
   - Require a < b, both mates' aligned bases inside [a, b), and 1 ≤ L ≤ `max_frag_length`, where L is the
     number of path bases in [a, b), counted as calibration counts them.

**Why the witness walks along annotated paths (review, blocker 4).** A walk along the observed path is not enough
when a clip crosses a gap neither mate aligned:

```
TA2               exon 1 =====|..........intron..........|===== exon 3
                             1200                       2000
true molecule        [1150====1200)                   [2000====2050)     100 bases, spliced on TA2
R  50M50S@1150       <====|ooooo                                         R's last 50 bases are clipped
F  50S50M@2000                                   ooooo|====>             F's first 50 bases are clipped
```

Walked contiguously, F's clip puts the start at 1950 and R's puts the end at 1250. That pair looks outward and
would be rejected. Walked along TA2, the start is 1150, the end is 2050 and the length is 100: a valid molecule.
When no path, annotated or observed, contains the pair, but a clip crosses an unobserved gap where an unannotated
intron could sit, the reason is **unresolved geometry**, not "impossible".

**What the witness is not.**

- It decides acceptance and nothing else. It never removes a candidate.
- It never changes a candidate's length `L_t`. TA1's 550-bp witness does not become TA2's length.
- It does not certify that every scoring interpretation fits the bound.

**Accepted placements and the count.** A pair is an accepted placement when it passes the validator and is neither
a baseline chimera nor a read-through (step 6). `num_hits` is the number of accepted placements. A fragment is a
multimapper when that exceeds 1, and `include_multimap=false` applies after assembly.

**What stays bit-identical.** A fragment keeps the baseline's buffer rows and per-hit emission exactly when all of
these are unchanged:

- its accepted placements are the same record pairs as the baseline's non-chimeric hits;
- its evidence (NM, motif) and splice tags;
- its hit count;
- its route.

The deliberate exceptions are:

- duplicate-NM resolution;
- withheld chains;
- motif-strand conflicts;
- geometry ordering;
- step 5's corrections.

Calibration's draw stays identical only where the held bank and its order are unchanged too.

## 6. Step 4: accounting for every template

Every template ends with exactly one **outcome**:

- **placed:** at least one accepted placement;
- **reported chimera:** none accepted, but a validated or aligner-asserted pair resolves as chimeric, read-through
  included;
- **unpaired:** neither;
- **excluded:** no eligible record, or dropped by the multimap option.

**Unpaired reasons.** Each R1 × R2 combination records how far it got through the validator; the template's reason
is the furthest stage any combination reached.

| reason | meaning |
|---|---|
| missing mate | one mate has no location at all |
| unresolved location | a mate's only locations are withheld |
| geometry | no combination has one reference, opposite orientations and containment on any path |
| unresolved geometry | containment fails only where a clip crosses an unobserved gap |
| contradiction | the mates disagree about a junction, and wave 1 cannot reconcile them (§7) |
| length | some path contains the pair, but none is within `max_frag_length` |

A template with an accepted placement is placed even if other combinations failed. Set-aside contradictions never
erase another accepted pair.

**A placed template records its route,** so that "placed" is never mistaken for "quantified". It says:

- whether it reached the EM with candidates;
- whether it is intergenic only (counted, not scored);
- whether it is spliced with no candidate (reported unsupported, never turned into gDNA);
- and in calibration, whether it was deposited, held, held as an artifact, or dropped.

The counts reconcile:

> templates = placed + reported chimera + unpaired + excluded
>
> records = placement records + unused records + records of unplaced templates

**Counters.**

- Today's `mate_unmapped` counts one-mate hit groups (review), not raw flags. It keeps that meaning under an explicit
  name, `one_mate_groups`, beside a new raw count of records flagged mate-unmapped.
- Rates are reported against all templates and against eligible templates, with units.

**The annotated BAM** writes every record once:

- an unpaired template's records get `ZC:unpaired` and the reason;
- unused records of a placed template get `ZC:unused_alignment`;
- excluded and chimeric records keep their own tags;
- stale assignment tags are cleared.

## 7. Step 5: reconciling mates that disagree

Two overlapping mates read the same molecule, so they must agree on its path. Aligners align mates separately and
sometimes disagree.

**The owner's case: end-extension.** X is spliced at [d, a); Y is not, and Y's aligned block runs k ≤ K bases past
the junction into X's intron, without covering the intron:

```
                    exon 1             intron             exon 2
genome          =========|.........................|=========
                         d                         a
X  (spliced)       ======|----------- N -----------|=====
Y  (as aligned)  ========>>>                        3 bases run past d into the intron
Y  (corrected)   ========|----------- N -----------|>>>   the same 3 bases, at [a, a+3)
```

**Rule.**

- Construct reading S: move Y's k bases across the junction, along X's path. No base is deleted. The mirror case at
  the acceptor moves them to [d − k, d).
- Recompute those k bases' mismatches against the stored reference at their new position. The existing mismatch
  term then scores them, so if they fit the far side worse, S pays for it automatically. There is no extra factor
  and no threshold.

**What it changes.** Only the overhang. On the toy annotation it falls 2 → 0 (TA2) and 3 → 0 (TA1). The review
recomputed both; the lengths and candidates are unchanged. The gain is k × 2.3 nats: 100× (k = 2) and 1,000× (k = 3)
more weight for the spliced transcripts. That can move abundance, compete with a remote placement and change pruning,
so it is measured in its own arm. (Revision 1 wrongly said it "costs nothing".)

**What it does not claim.** S is constructed, not proven. The contiguous reading C, in which X's few far-side bases
belong beside its near block, may also be valid. Weighing C against S is wave 2.

**Every other disagreement is a contradiction in wave 1:**

- an overhang beyond K;
- a block spanning the whole skip;
- a block lying inside the other mate's skip;
- two overlapping skips that differ.

The pair is set aside and counted, the template's other placements stand, and wave 2 recovers it. This is a
measured, temporary development state, not a release solution. In the census (Appendix A), the classes beyond K are
30–89% gDNA or nascent artifacts on minimap2. The exception is the k = 6–8 bump, where 93–94% of templates are of
mRNA origin. Origin is provenance, not proof that a correction would be right; the census is recounted after full
pair validation, with every class a template shows.

## 8. Step 6: read-through chimeras

The owner's rule: a pair joined across two genes by an unannotated splice is a read-through chimera, reported and not
quantified.

```
gene A (+)  ===exon===|.........                         .........|===exon===   gene B (+)
pair             =====|---------------- N ------------------|=====
                 left side: exons of gene A only      right side: exons of gene B only
```

**The predicate (review, finding 8).**

1. For each junction of a reconciled reading, let L be the genes with an exon under the aligned block that ends at
   its donor side. Let R be the genes with an exon under the block that starts at its acceptor side.
2. If both are non-empty and share no gene, the reading is an inter-gene chimera, whatever its span. It is called a
   read-through when the two genes share a strand.
3. If either side has no gene, it is not inter-gene: it is an unannotated splice of one gene, quantified as today.

Consequences:

- **An annotated junction is never chimeric:** the transcript that has it puts its gene on both sides. This covers
  GENCODE's read-through genes.
- **The common-gene test is strand-agnostic,** so a strand disagreement can never manufacture a chimera.
- **Unspliced pairs spanning two genes stay ordinary molecules.**
- **A chimeric reading is never a placement,** but the template's other placements stand.
- **The cost is measured three ways:** true read-throughs excluded by the policy, ordinary RNA misclassified, and
  artifacts removed. Long-span cases (3,544 on VCaP) were already excluded by the baseline, so only the short-span
  ones are new.

## 9. Wave 2: competing readings

A pair can often be read two ways at one place:

```
reported   X  ===========|-------- N --------|===     spliced, with a 3-base far anchor
           Y      ============>>>                     contiguous, 3 bases into the intron
reading S  X  ===========|-------- N --------|===     both take the junction (RNA only)
           Y      ==========|-------- N -------|>>>
reading C  X  ===========>>>                          both contiguous: X's 3 bases belong beside its near
           Y      ============>>>                     block (gDNA can explain this one)
```

Local contiguity (source A), the contradictions of §7, and remote copies all produce such alternatives. They are the
splice-artifact mechanism. The review showed they cannot go through today's multimapper path. Three rules come
first; each is a derivation checked against finite enumeration before any code.

### Rule A: how readings and placements add up (review, blocker 1)

Today a multimapper takes the **best** hit for each transcript and the **mean** gDNA term over hits. The reviewer
measured what that does:

- two distinct origins of one transcript gain nothing, where a sum would gain ln 2;
- adding a gDNA origin with relative likelihood 10⁻⁶ halves the gDNA term, where a sum would barely move it.

The proposal to derive and test distinguishes two kinds of "several":

| kind | what it is | how they combine |
|---|---|---|
| distinct physical locations | the molecule sat at one of several places | add |
| alternative readings at one location | which path the molecule took | under transcript t, t's own path is the reading; under gDNA, only a contiguous reading exists at a location |

The gDNA sum needs each location's opportunity and prior, intergenic support, and capture per placement. That is
the parked multimapper-likelihood design, so the alternatives depend on it.

### Rule B: who owns the sequence evidence (review, blocker 2)

**Proposal.** A reading's sequence evidence is the mismatch count of both its mate alignments against the stored
reference, over the same observed bases, scored by the **existing** mismatch term (NM × ln α_mm). It adds no factor
and no new estimand:

- the same final reading scores the same whatever parent it came from;
- nothing is counted twice;
- a tie (exactly equal mismatches) is the only "no evidence" case. Revision 1's "Δm ≤ 0" conflated ties with
  improvements.

The review's derived alternative replaces α_mm with the library's own `ε/[3(1−ε)]`. That is an emission-model change
with its own A/B; it is listed as a decision.

Bases that one reading clips and another aligns need an explicit treatment. A reading must not win by explaining
fewer bases.

### Rule C: calibration support for held readings (review, blocker 3)

Calibration's second pass weighs a held fragment's readings by pass-one deposited traffic. A zero for one reading is
decisive. If every fragment supporting a junction is held, S reads zero and loses to C by construction: the held
evidence manufactures a decision against its own RNA path.

The rule to derive: support that is **unavailable** because it is held is not zero. The second pass needs independent
support, or a coherent joint fit, and a test where an entirely held junction keeps its RNA reading. The second pass's
length law must also follow observed as well as implied introns, and calibration's spliced channel must agree with
the EM on wholly unannotated observed splices.

### Then, in order, one arm each

1. S and C for §7's contradictions, built jointly over both mates and every junction, in both anchoring directions.
   Partial RNA readings (keep a real junction, repair an artifact) are allowed. A gDNA reading must explain every
   junction contiguously.
2. Local contiguity alternatives (source A): ties first, then any mismatch difference.
3. Intergenic placements in the EM's gDNA term (the owner's "honor the gDNA hypothesis").
4. Remote copies from the index-time copy table (source B).

## 10. Wave 3 and the release list

These are independent of wave 2, each designed and measured on the repaired scan:

1. **The gDNA length estimator:** its conditioning, its deflated internal density, observability (not only N; two
   ends can span an N gap, and repeats can be unobservable), and its common-contaminant assumption.
2. **False capture in the ruler,** at zero gDNA and in the capture-OFF tails.
3. **The aligned gDNA/RNA split** in the unstranded stratum.
4. **The clip length model:** molecule ends walked along each candidate's own path.
5. **Observable opportunity, and the scope of sample-specific alleles (source C).**

The blacklist stays until a matched two-error comparison justifies removing it. Release waits for all of it (owner).

## 11. What the review changed

| finding | change |
|---|---|
| blocker 1: tied alternatives meet best-hit RNA and mean gDNA | alternatives leave wave 1; Rule A first |
| blocker 2: preference is not identification; Δm conventions; NM double-counting | S is a constructed reading, not a proof; one owner of sequence evidence (Rule B); ties mean exactly equal |
| blocker 3: calibration can zero a held RNA reading | Rule C before activation; wave 1 adds no alternatives |
| blocker 4: the observed clip walk is not exact | containment moved into the path witness, along annotated paths; new reason "unresolved geometry" |
| 5: location identity underspecified | full CIGAR mapping, explicit equivalence, recomputed NM, sequence-recovery rules |
| 6: the ledger is not sufficient | outcome, reason over all combinations, route for placed templates; withheld ≠ missing; `mate_unmapped` corrected |
| 7: diagnosis overstated | mediation, not "the whole"; capture-OFF ruler tails; crossed law × ruler arms; consumed arrays hashed |
| 8: read-through gene membership | flank gene sets, no common gene, strand-agnostic test |
| 9: enumeration cost | exhaustive reference first; no caps without a proven bound; overflow checks now |
| order | reference store first; alternatives after Rules A–C; policy costs explicit |

## 12. Tests

One suite per step. Each fixture is labelled as one of two kinds:

- **invariant:** passes on the baseline and after;
- **defect demonstration:** fails on the baseline, passes after the fix, and fails again when the fix is broken on
  purpose.

**The toy annotation** (0-based, half-open; `counterexample_baseline.py` builds it):

| transcript | strand | exons |
|---|---|---|
| TA1 | + | [1000,1200) [1500,1700) [2000,2300) |
| TA2 | + | [1000,1200) [2000,2300) |
| TB1 | + | [5000,5150) |
| TC1 | − | [2200,2400) [2600,2800) |

**Emission values the review independently recomputed** (exon_bp / overhang / L):

| case | expected | kind |
|---|---|---|
| one-base overlap, [700,800) / [901,1001) | TA1, TA2, nascent A: 1 / 199 / 301 | invariant |
| both outer boundaries | TB1 100 / 100 / 250 | invariant |
| one overhanging mate | TA1 150 / 50 / 200; TA2 50 / 150 / 500; nascent A 200 / 0 / 500 | invariant |
| long-intron witness | TA1 200 / 0 / 550; TA2 200 / 0 / 350; nascent A 200 / 0 / 1,150 | invariant |
| acceptor-side correction | TA2 148 / 2 / 148 → 148 / 0 / 148 | demonstration |
| donor-side correction | TA1 147 / 3 / 147 → 147 / 0 / 147 | demonstration |

**Cases per suite:**

- **Reference store:** round trip; N and ambiguity codes; soft-masked case; contig ends; reverse complement; a wrong
  contig checksum; a missing store.
- **Locations:**
  - duplicates within and across HI groups;
  - homopolymer-shifted indels;
  - `5S95M` against `95M5S` at one start;
  - conflicting NM and motif tags;
  - missing NM or bases;
  - a hard-clipped carrier;
  - a supplementary chain without SA;
  - adjacent H and S;
  - a record with no M/=/X base.
- **Pairing:**
  - R1-forward and R2-forward;
  - minus-strand mirrors;
  - outward and protruding pairs;
  - the clipped pair across TA2's junction (accept);
  - the clip across an observed skip (accept);
  - a contiguous-only walk, which must reject the TA2 case (the break test);
  - K and K + 1 edges; the length bound and bound + 1;
  - a raw blacklisted skip;
  - an RNA witness whose gDNA span is long;
  - a gap intron that would erase aligned bases (not a witness);
  - an orphan secondary;
  - cross-HI pairs;
  - hit-count overflow.
- **Accounting:**
  - one fixture per reason;
  - a geometry failure beside a too-long consistent pair;
  - an unresolved pair beside a valid alternative;
  - a placement with no transcript;
  - unused aliases;
  - every record written once, and stale tags cleared.
- **Reconciliation:**
  - both sides, with k = 1, K and K + 1;
  - the moved bases matching or mismatching the far side;
  - an indel at the repair boundary;
  - a spanning block;
  - a block inside the other mate's skip;
  - shifted skips with one or with two annotated placements.
- **Read-through:**
  - a novel junction inside one gene;
  - two nearby genes, at short and long span;
  - nested and antisense overlap;
  - an annotated read-through gene;
  - one side with no gene;
  - a template with a valid non-chimeric alternative.

**Wave 2 cases** come with Rules A–C:

- multi-junction S and C;
- two readings sharing a transcript with unequal overhang;
- a remote competitor;
- a duplicated report (no probability change);
- a negligible extra origin;
- an entirely held junction;
- canonical draw order;
- spill and reload.

**Break tests:**

- restore the one-mate hit;
- restore `max(NH, hits)`;
- skip asserted-pair validation;
- bound candidates by exon ends (the continuation's error);
- prune candidates by `L_t`;
- skip `detect_chimera`;
- merge one mate's HI records;
- walk clips contiguously;
- delete rather than move overhang bases;
- double-apply NM;
- make Δm depend on the parent;
- fabricate junction support;
- drop raw skips or sequence fields when serializing;
- omit filtered templates from the output;
- drop the read-through rule.

## 13. How each step is judged

**The A/B contract (from the review).** Each comparison varies one mechanism against its parent and the baseline. It
uses:

- `scan.total_threads=1`, `OMP_NUM_THREADS=1`, fractional assignment and recorded calibration seeds;
- the gDNA law and the ruler frozen in both arms, with a hash of the array each consumer actually received;
- scoring against truth.

It reports both errors — artifact RNA mass and genuine RNA lost — plus transcript, gene and calibration error,
separately for the three in-scope strata, with the deferred fourth beside them. Every arm also gets:

- a per-template transition census, joined on read id;
- the crossed law × ruler design, wherever attribution matters.

**Implementation gate:**

- the revised equations, checked by independent finite enumeration of paths and origins;
- the identity domain and its exceptions written down;
- complete location, candidate and disposition schemas;
- a failing test for every claimed defect;
- injection hooks proven to reach their consumers.

**Release gate:**

- every isolated repair passes its tests, and passes again when broken back in;
- the full suite and the cache and persistence checks pass;
- the whole system passes:
  - the truth-scored aligned and oracle ladders;
  - the zero-gDNA controls;
  - both length-gap panels;
  - STAR and minimap2;
  - VCaP transitions;
  - independent artifact-containing validation.

Policy costs are quantified and accepted explicitly, never hidden behind a tolerance. Examples: the templates the
complete-pair rule leaves unpaired, and the read-throughs excluded by policy.

## 14. Decisions

**Recorded (owner):**

- complete pairs;
- fix everything before release;
- scan stage first, with separate tests;
- read-through chimeras are reported, not quantified;
- mate overhangs are corrected using K;
- the genome is stored in the index.

**Open:**

1. **Approve wave 1** (steps 1–6) for prototyping. The reviewer approves the complete-pair architecture and the
   overlap-preserving witness in principle.
2. **Read-throughs at any span**, with §8's flank predicate. The reviewer recommends yes.
3. **Beyond K:** the reviewer recommends proposing and scoring S and C (wave 2), not correcting automatically, and
   not choosing a new cutoff. Agree?
4. **The owner of sequence evidence (Rule B):** the existing mismatch term with recomputed mismatch counts (my
   proposal), or a factor derived from the library's own error rate (the review's). Either way, the term counts the
   same evidence once.
5. **Order after wave 1:** wave 2's three rules first, or wave 3's length estimator and ruler first. They are
   independent; the length estimator is where the largest measured error sits.

## Appendix A: the measurements

Each table entry gives transcript / gene error, summed absolute counts, on the local aligned minimap2 panel. Every
run is pinned and fractional.

**Capture ON, stranded:**

| build | own law and ruler | gDNA law fixed | ruler fixed |
|---|---|---|---|
| baseline | 564,547 / 342,837 | 592,114 / 341,312 | 793,321 / 575,065 |
| stage 1 | 801,015 / 341,628 | 577,866 / 330,918 | 685,296 / 479,161 |
| stage 3 | 825,893 / 331,159 | 557,556 / 316,653 | 668,655 / 468,399 |
| oracle BAM | 216,127 / 103,401 | | |

Swaps: the baseline's ruler in stage 1 gives 488,421 / 275,533; stage 1's ruler in the baseline gives
907,999 / 442,235. The "ruler fixed" column uses the oracle run's ruler, whose scale differs, so read it only down the
rows.

**The contrast** (from the review's recomputation). Stage 1's inversion is 63.8·f0 − 62.8·f1, with a signed mean of
740 bp and negative mass 2.48 before projection. Projection, blending and smoothing then give 318. On the oracle,
purity near 1 makes it one-pool.

**Capture OFF:**

| | unstranded | stranded |
|---|---|---|
| baseline, own law | 260,773 / 116,794 | 309,282 / 112,797 |
| baseline, law fixed | 310,430 / 142,670 | 309,207 / 129,708 |
| stage 3, own law | 277,582 / 144,709 | 249,937 / 112,046 |
| stage 3, law fixed | 266,983 / 137,019 | 249,401 / 109,135 |

**Zero gDNA:**

| build | own ruler | flat ruler |
|---|---|---|
| baseline | 5,076,734 / 705,450 | 315,241 / 67,603 |
| stage 3 | 4,457,600 / 778,794 | 316,497 / 74,169 |

**The ruler's tails at capture OFF.** Transcripts off by more than 2×, and by more than 10×:

| input | > 2× | > 10× |
|---|---|---|
| aligned, baseline | 1.9–2.3% | 1.0–1.1% |
| aligned, stage 3 | 1.6–2.1% | 0.7–1.0% |
| oracle BAM | 0.3–0.4% | 0.3–0.4% |

**Mate disagreements** (primary pairs, first class per template; recount pending):

| class | minimap2 OFF: mRNA / artifact | minimap2 ON: mRNA / artifact | VCaP: RNA library / DNA library |
|---|---|---|---|
| end-extension, k ≤ 3 | 495 / 58 | 687 / 20 | 2,165 / 5 |
| end-extension, k > 3 | 996 / 423 | 946 / 1,023 | 247 / 90 |
| block inside the other's skip | 561 / 451 | 646 / 1,358 | 0 / 0 |
| block spanning the skip | 87 / 49 | 5 / 41 | 10 / 10 |
| all spliced pairs | 1,514,991 / 8,607 | 1,945,663 / 7,903 | 1,493,592 / 1,453 |
