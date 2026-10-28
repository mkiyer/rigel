# Splice artifacts — gDNA the aligner splices, and what Rigel does with it

Sandbox document for the session that takes `ISSUES: splicing-artifacts`: release-blocking, and step 3 of `ROADMAP.md`'s
Next, the scan stage, the first step in the tool's order that runs on this machine (step 1's remainder, the STAR rung's
own `alignable` store, and step 2 are cluster tasks). Nothing here is authoritative. The open record is that entry, with
every measurement; the scope is `DESIGN.md` §1.0 (Rigel reads alignments, never repairs them); the rulings this work
builds on are `DESIGN.md` §1.1 (a template's hits), §3.1d (which alignments gDNA may explain), §3.3b (κ from the
genuine junctions) and §3.1g with `EQUATIONS.md` §15 (a multimapper's places). This document replaces the issue's
design trail and the owner's comments of 2026-10-16, removed on 2026-10-28; section 9 says how to recover each.

## 1. The defect

An aligner writes a junction on a molecule that was never spliced there: a gDNA read whose short end happens to match
the far side of an annotated intron, or a read from a processed pseudogene that also fits its parent gene's junction.
Rigel takes a sequenced junction as certified RNA. gDNA cannot explain it (`DESIGN.md` §3.1d); calibration deposits it
as spliced RNA, so it anchors the RNA side of the deconvolution; it trains the strand model's κ and RNA's
fragment-length law; and in the EM it has RNA candidates only. One artifact is booked as RNA and also trains the models
that split every other fragment. On the index case, cfRNA-D (96.4 % gDNA), about 43 % of the RNA-called spliced
fragments were artifacts, and the library read a strand specificity of 0.786, last of 344.

The guard is the index's splice blacklist, from the `alignable` tool: junctions the same aligner wrote on simulated
genomic reads, rejected when a read's anchor is no longer than the longest seen there. It cuts both ways: on the VCaP
mix's RNA half it rejected 598,010 records where a pure-gDNA library predicts about 20,000, and a rejected genuine read
becomes unspliced and gDNA-eligible at a footprint that spans the intron. It also misses: on the index case 51.6 % of
the spliced fragments' junction rows had no entry, and another 10.8 % were listed but escaped.

Why now (owner, 2026-10-21, `ROADMAP.md` Next): "gDNA that aligns as spliced is the main error on real data". The
recommended STAR settings report the competing alignments and the multimapper handling weighs them; that handling
landed on 2026-10-28 (`DESIGN.md` §3.1g), so this issue is what the scan's repair was building toward. 0.8.0 waits for
the mechanism and for its validation on independent artifact-containing evidence (owner, 2026-10-14), which is the
cluster's.

## 2. What has been said, in order

- **2026-08-03, the first investigation's principle** (on the index case): "The system must be robust to aligner
  settings. STAR parameter tightening is a mitigation, not a fix." And: the repair "matters for a second-order reason:
  it de-contaminates the priors, not primarily the spliced counts". `DESIGN.md` §1.0 has since settled the input: STAR
  at the recommended settings, and anything else read as reported.
- **2026-09-27, the owner's plan for the cluster.** The work is "judged on BOTH errors together: artifacts certified as
  RNA, and genuine junction reads rejected". "Every mechanism here is a hypothesis until it is measured on the
  libraries. No magic numbers, no yes/no gates on continuous quantities." "Never read the truth-labelled halves
  pooled, and never read one error alone." "The shares rank the mechanisms."
- **2026-10-14, owner.** Spliced fragments are the most specific signal the tool has, so artifacts are removed by
  several prongs before any release; 0.8.0 waits for the mechanism and its validation.
- **2026-10-15, owner.** The VCaP mix's RNA library is real and holds some gDNA of its own, so a read from it is
  provenance, not proof of genuine splicing: read its losses as an upper bound, with the split the strand implies. And:
  "whenever a fragment aligns to an intergenic region we MUST honor the gDNA hypothesis".
- **2026-10-16, morning, owner** (`my_notes.md`). "Our method should admit BOTH "S" (the RNA hypothesis) and "C" (the
  gDNA hypothesis). Both hypotheses should be allowed to propagate through and ultimately compete in the EM." Then: "I
  agree that 'S' and 'C' need to compete in the EM. Why is this design parked? Rigel supports multimapping everywhere
  now". On evidence: "Currently the sequence evidence is 'overhang' and 'mismatches'. What else? These can be computed
  for every paired-end alignment possibility." Two ideas, both held: aligning every fragment twice, spliced and
  genome-only, merged by read name ("the "pure" way … this may be overkill and it has its own potential problems"); and
  aligning the annotation's own sequences back to the genome at index build, which "would supersede the alignable
  tool's splice blacklist".
- **2026-10-16, later, owner: the scope** (`DESIGN.md` §1.0). "The hope is that Rigel does not need to become an
  alignment repair utility or a second alignment tool. Rigel will document the optimal STAR parameters." The genome
  store was dropped; a dual-aligner input is the fallback; the blacklist is the one declared exception.
- **2026-10-17/18, owner.** Splice-artifact behaviour is measured on simulated reads aligned at the recommended
  settings, because there the truth is exact; those settings report both readings of a template whose alignments score
  alike, and they are final.
- **2026-10-19, owner.** Unannotated junctions are withheld (the annotated-path rule, `DESIGN.md` §1.1), and the
  withheld population must stay buildable for a gene-level salvage
  (`ISSUES: fragments-that-fit-one-gene-but-no-transcript`).
- **2026-10-27, owner** (the multimapper precedent, `DESIGN.md` §3.1g). gDNA called RNA is the error to cut, and RNA
  called gDNA the accepted risk at a measured exchange. That acceptance was for multimappers. No genuine-loss budget
  has been set for splice artifacts: a tradeoff here needs its own decision.

## 3. What changed since the trail, and what each change does here

| change | landed | what it does to splice artifacts |
|---|---|---|
| Rigel reads, never repairs (`DESIGN.md` §1.0) | 2026-10-16 | Rules out the trail's machinery that reads bases: the genome in the index, a tie table with a Rigel-built second hit, the hold rule decided on the bases, the sequence emission from SEQ, h computed from the reference, and an index-time copy table unless the VCaP mix shows STAR missing the copies. A contiguous reading must come from the aligner |
| The recommended STAR settings | 2026-10-17/18 | `--sjdbScore 0` removes the annotated-junction bonus and `--outFilterMultimapScoreRange 2` widens the reported range, so a spliced reading and a contiguous one that score alike are both reported. Every VCaP number in the entry was measured on a BAM aligned at STAR's defaults (its `Log.out` sets neither option), and the index case's on 0.7.1 at the production settings |
| The annotated-path rule (`DESIGN.md` §1.1) | 2026-10-19 | A spliced hit is usable only if one annotated transcript holds every junction it observes; a template left with none is withheld (`unannotated_splice`). The unannotated class, 72 % of VCaP's harmful artifacts, no longer reaches RNA, and its gDNA is no longer counted anywhere |
| The place model (`DESIGN.md` §3.1g) | 2026-10-28 | A template's hits compete in the EM: overlapping hits on one reference are one place, read at each component's best hit; distinct places sum; intergenic places are a background at the bulk rate. A retrocopy origin the aligner reports is a distinct place, and a contiguous hit at the same locus merges with the spliced one (section 4.3) |
| The STAR rung | 2026-10-21 | The ladder aligned by rustar at the recommended settings and scored against read names: a real aligner's artifacts with exact truth. Its suite index has no `alignable` store, so every run on the rung is the no-blacklist arm |
| The observable genome | 2026-10-22 → 27 | The GRCh38 index (format 11) now carries the store's mask and its blacklist, 5,833,092 junctions with counts, so a VCaP run on it is the production-blacklist arm, not the no-blacklist arm the trail described |

## 4. The code today (HEAD `84897128`)

### 4.1 The scan
- **Junctions.** `parse_cigar` (`bam_scanner.cpp:643-690`) makes each `N` a junction with two anchors; the left anchor
  counts every reference-advancing base before it, cumulative across earlier introns (`:671`), as `alignable` counts
  them. One junction strand per record, from `XS` or minimap2's `ts` (`read_sj_strand`, `:724-796`).
- **The blacklist.** Loaded only when the index records an `alignable` store (`index.py:1675-1709`); strand-blind, and
  the `count` column is never read (`resolve_context.h:607-622`). Applied per record, before the mates are joined: a
  junction is erased iff `anchor_left <= max_anchor_left || anchor_right <= max_anchor_right`
  (`bam_scanner.cpp:709-716`), so an entry whose one side is 0 tests the other side only, and the record's blocks stay
  split. A junction kept on either mate survives (`:891-905`). A hit becomes ARTIFACT only when none of its junctions
  survived (`resolve_context.h:1401-1403`), so one kept junction shields the rejected ones.
- **The annotated-path rule** reads the junctions the blacklist left (`hit_verdict`, `bam_scanner.cpp:1124-1141`), so
  a blacklisted novel junction leaves an ARTIFACT fragment, not a withheld one.
- **What a buffered row keeps** (`resolve_context.h:197-277`): splice type, strands, `num_hits`, `genomic_start`,
  `read_length` (aligned bases), `genomic_footprint` (introns and mate gap included), `nm` (summed over the hit's
  records, `bam_scanner.cpp:921`) and per candidate its length and exon bases; a multimapper's place masses. Not kept:
  junction ids, blocks, soft clips, MD, SEQ or qualities. No `bam_get_seq` exists in `src/`.
- **ARTIFACT** is censused, then never offered to calibration (`bam_scanner.cpp:1846-1864`). In the EM gDNA reads it at
  its intron-spanning footprint, and a transcript holding the rejected intron reads it at its spliced length.

### 4.2 Calibration
- A survivor with an annotated junction deposits as spliced (`spliced = !sj_ids.empty()`, `accumulator.cpp:587`):
  spliced counts at its boundaries, the junction's mass, and the `kRnaSpliced` length pool that becomes RNA's length
  law (`fl.py`'s own-count law, then `shared_length_law`). A multimapper's drawn placement deposits the same way
  (`DESIGN.md` §4.2), so the pool holds drawn spliced placements too. Where the shared law keeps one law, RNA's spliced
  census is both origins' table, so its artifacts reach gDNA's table as well (not measured).
- **κ** comes from unique, strand-qualified spliced fragments, one observation at the leftmost annotated junction
  (`bam_scanner.cpp:2126-2144`), fitted by `junction_fit` as three classes (genuine at κ, artifact at ½, reversed at
  1 − κ). Its per-junction responsibilities stay inside the fit and its class shares go to `summary.json` only
  (`cli.py:386-410`): nothing reads a per-junction artifact posterior.
- A multimapper's placements are held whole and drawn in the second pass (`DESIGN.md` §4.2), which carries no NM
  (`accumulator.h:251-260`).

### 4.3 The EM
- A SPLICED_ANNOT fragment has RNA candidates only (`gdna_can_explain`, `scoring.cpp:221-232`); unique with one
  candidate, it skips the E-step (`:796-836`).
- RNA's term is strand + f_R(length) + overhang · log α_oh + NM · log α_mm, plus log w on a multimapper; gDNA's is
  f_G(footprint) + log ½ + NM · log α_mm (`scoring.cpp:864-898, 942`). `overhang_alpha` and `mismatch_alpha` (0.1
  each) and the 1e-4 prune are run settings in `FragmentScoringConfig` (`config.py:110-157`), underived
  (`ISSUES: scoring-penalties-are-underived-constants`). Soft clips are priced nowhere.
- **S and C reported together at one locus** (two usable hits, a multimapper). A transcript holding the junction is a
  candidate of both, their extents overlap, so it is one place read at the better of its S term and its C term, and C
  pays overhang for its intronic bases. The nascent span and any intron-retaining isoform read C only, and gDNA reads C
  only (`merge_places`, `scoring.cpp:428-456`). The odds between the transcript's S and gDNA's C go as
  α_mm^(NM_C − NM_S) times the strand, length and place-weight terms. No gate pins this case yet.
- If the blacklist also rejects S's junction, S is ARTIFACT: held out of calibration, so C deposits as a unique
  mapper; in the EM S adds a gDNA term at its spanning footprint, merged with C's.

### 4.4 The index
`splice_blacklist.feather` holds ref, start, end, the two longest anchors and, since format 9, a count, one row per
junction, built from the store's rows per (chrom, intron, read length) after a `count >= 2` cut on each row
(`splice_blacklist.py:87-112`). The index stores no genome sequence (`index.py:132`).

## 5. What is measured

The entry holds the record: the index case, the catalogue's defects, the strand prong, h, the harmful class, the
blacklist both ways, remote origin, the copy table and the Rigel-built second hit's price. Three things to carry into
reading it:
- Every number predates the recommended settings: the index case was read on 0.7.1, and every VCaP number on the mix as
  STAR's defaults aligned it. With the bonus gone, ties that STAR wrote spliced alone may now come with their
  contiguous alternative, so the h = 0 class is the one most likely to have moved. Nothing has been measured since.
- The unannotated class, the largest harmful one, is withheld now: it no longer inflates RNA, and it no longer counts as
  gDNA either.
- The Rigel-built second hit's side effect returns with an aligner-reported one: a contiguous reading also fits an
  intron-retaining isoform or a longer exon, and it moved one ADAR isoform from 3,545 to 134 fragments.

## 6. What survives of the trail's design

- **The likelihood.** Every candidate is scored once as (θ_k / ℓ_k) · c(P) · f_k(L) · s_k(o) · e(x | P), and
  alternative origins add, never multiply. `DESIGN.md` §3.1g now does this for the hits the aligner reported
  (`EQUATIONS.md` §15.1). The aligner proposes paths; it does not testify: its score and MAPQ are never a probability.
- **Class evidence** (a junction's known ends, its h) enters only as a likelihood ratio, odds(class share) /
  odds(library share), with matched opportunity. Strand is used once, in the EM's s(o): a per-junction strand posterior
  would reuse the fragment's own strand (the review's finding 1). On an unstranded library a ratio of one means no
  evidence, not half artifacts. The wrong-strand excess, q = κ + (½ − κ) a + (1 − 2κ) b, cannot identify both the
  artifact share a and the reversed share b (finding 6).
- **The sequence ratio** between a contiguous and a spliced reading of the same bases is [ε / (3(1 − ε))]^(m_s − m_c),
  where m_s and m_c count the discriminating positions the read matches under each reading: exactly 1 at h = 0,
  whatever the read's errors. Rigel has no SEQ under `DESIGN.md` §1.0, but where the aligner reports both readings over
  the same bases (no clip, no shifted indel), NM_C − NM_S should equal m_s − m_c (an inference, not checked), which
  makes the derived per-mismatch factor ε / (3(1 − ε)) against today's α_mm = 0.1. ε, Σ(NM − I − D) / Σ(M, =, X) over
  uniquely mapped unspliced fragments, is alignment divergence rather than pure sequencing error, and it must come from
  the library, never from `mismatch_alpha`.
- **A tie is ambiguous, not disproved.** With B the length-and-opportunity ratio and q the prior's gDNA share, the MAP
  fixed point is z = B(z + q) / (B(z + q) + 2 − z − q): at B = 1, z = q; at q = 1 and B = 0.05 it leaves 0.947 on RNA.
  The mechanism claims nothing at a tie, and the prior does not decide a singleton.
- **The stage contract.** Share candidate geometry and observation evidence between calibration and the EM, never the
  realized assignment: fit from pass one, score once, draw or deposit once, then fit the final models. Training weights
  are the probability of RNA under a valid RNA geometry, not 1 − P(gDNA).
- **Multi-junction fragments**: "a DNA candidate must explain the whole fragment without RNA splicing … independent
  per-junction penalties do not establish that combined path's existence". Never derived; 32,019 of the VCaP RNA
  library's 82,001 h = 0 units are multi-junction.
- **Superseded**: the first proposal's weight r (a posterior share multiplied into the EM, review finding 1); the
  largest block kept, the short anchor clipped (finding 2); pruning on r alone (finding 3); every Rigel-built
  contiguous candidate (sources A and B, `DESIGN.md` §1.0); the hold rule decided on the bases; the genome in the
  index; a sample allele at the spliced path (source C, "not established").

## 7. The open questions

1. **What the aligner's alternatives and the place model already do**, before any mechanism. On the STAR rung, every
   gDNA template that reaches RNA through a spliced hit and every RNA template that reaches gDNA, by class (annotated
   short anchor at h = 0, annotated long anchor, mixed, implicit, ARTIFACT), by whether the aligner reported a
   contiguous alternative, and by where the EM put it, per stratum. This sizes everything after it.
2. **One place, two readings.** S and C merge at the better term for a transcript holding both (section 4.3), and
   alternatives at one place count once (`DESIGN.md` §3.1g, a ruling). Does a reported C pull genuine RNA into an
   intron-retaining isoform, the nascent span or gDNA, and how much?
3. **What decides a tie** (h = 0, equal NM): strand, length and abundance. Is the fragment's own strand term enough, or
   does the junction's strand record add evidence without reusing it? `DESIGN.md` §3.3b's limit, "The per-junction
   artifact posterior is the input to pruning", contradicts the review's finding 1 and needs the owner.
4. **A spliced fragment the aligner reported alone** stays certified RNA. How large is that class on the rung, and is it
   the aligner's to fix (`DESIGN.md` §1.0)?
5. **The blacklist**, the one declared exception: retire it, or replace its yes/no cut with a likelihood. Its catalogue
   was simulated at STAR's production settings, not the recommended ones, and carries its own defects (the entry's
   THE CATALOGUE); `IndexConfig` says a minimum count of 2 matches `alignable`'s own threshold, which defaults to 1
   (unchecked). The matched comparison needs both arms on one input: VCaP has the production arm, and the rung needs a
   store for the suite genome, built on the cluster. Its untested rules matter only if it stays.
6. **Calibration and the training sets.** A spliced fragment with a reported contiguous hit is held and drawn without
   NM; should the draw read the hit's evidence? κ and RNA's length law train on every spliced fragment, artifacts
   included, and the shared law can carry them into gDNA's table (section 4.2).
7. **The withheld unannotated class.** Its gDNA leaves the gDNA count; is that the accounting wanted, and would the
   gene-level salvage re-admit artifacts?
8. **The per-mismatch factor and soft clips**: derived ε / (3(1 − ε)) against α_mm = 0.1, and clipped bases that cost
   nothing (`ISSUES: scoring-penalties-are-underived-constants`).
9. **Multi-junction and implicit fragments**, never derived.

## 8. Data, instruments, constraints

- **The STAR rung**: `$RIGEL_RUNS/suite/ladder/<condition>/rustar_recommended.bam` with its `SJ.out.tab` and
  `.align.json`, the simulator's own `sim_oracle.bam` (no aligner artifacts) and `truth_abundances.tsv`; the suite index
  `$RIGEL_RUNS/suite/rigel_index` (format 11, the N runs as its mask, no blacklist). `panel.py score --bam
  rustar_recommended` and `quant_accuracy.py --bam` score it; per-fragment outcomes come from the annotated BAM, read
  against the read names. Harm is read on the oracle ladder, both fl-gap panels and the weak-strand flips
  (`$RIGEL_RUNS/prototypes/2026-10-22_census_law/flip_bams/`).
- **The VCaP mix**: `$RIGEL_RUNS/cfrna/vcap_mix/`, a BAM aligned by the lab pipeline at STAR's defaults, truth by read
  name (each library's instrument prefix), its 0.7.1 `summary.json`, and its RNA-only half; the GRCh38 index
  `$RIGEL_RUNS/refs/rigel_index` (format 11, the store's mask and blacklist). The mix at the recommended settings needs
  a re-alignment, and there is no local GRCh38 STAR index; the cluster has one. cfRNA-A, -B and -C carry their BAMs,
  `summary.json` and `SJ.out.tab.gz` (`$RIGEL_RUNS/cfrna/`); every 0.7.1 `annotated.bam` is gone.
- **The trail's measurement scripts**: `$RIGEL_RUNS/prototypes/2026-10-14_splice_design/` (`strand_prong.py`,
  `h_prong.py`, `seq_prong.py`, `junction_copies.py`, `retro_reads.py`, `unannot_ends.py`, `extract_spliced.py` and
  their outputs), `2026-09-26_splice_artifacts/` (the h table by short side), and the superseded Rigel-built second
  hit, `2026-10-15_splice_a/`. The index case's evidence bundle is gone.
- **The cluster**: its index must be rebuilt at format 11 naming its store before any cluster quant (ROADMAP step 2,
  `ISSUES: an-unrecorded-splice-blacklist-is-dropped-silently`); the lab pipeline still passes `--nrna-tolerance`, now
  `--nascent-tolerance`; `alignable`'s per-record rows, which a fix of its half-blind entries needs, survive only in an
  unpacked run (`splice/chunk_*.feather`), and the STAR index it simulated with is deleted. The aligned ladder's draft,
  `$RIGEL_RUNS/aligned_ladder/`, predates the STAR rung; `$RIGEL_RUNS/to_cluster/`.
- **Rules**: one whole-genome process at a time; every A/B pinned and fractional, read per stratum against truth, with
  every capture efficiency at 1 beside the shipped run (the capture reader amplifies any change,
  `ISSUES: the-reader-reads-false-gdna-as-capture`); a failing hand-built fixture first, then each gate broken on
  purpose; one mechanism per arm.

## 9. The trail, for recovery

`git show 84897128:docs/dev/<name>`, copies in `$RIGEL_RUNS/archive/docs_dev_2026-10-28/`:
- `GDNA_SPLICE_ARTIFACTS_OVERVIEW.md` (2026-08-03): the index case's tables (escape routes, the two artifact classes by
  anchor, mismatch and mapping, the cohort's correlations). `_PLAN.md` (2026-09-27, the owner's plan for the cluster):
  the VCaP two-sided table, the h table by short side, the phases. `_RIGEL.md`: the site-pair rule's measurements, the
  annotated-unlisted rows by anchor, the counters' units. `_ALIGNABLE.md`: the catalogue's facts and six defects, and
  two hypotheses refuted. `_STAR.md`: the effective-parameter audit and five recommendations, none adopted.
- `SPLICE_ARTIFACTS_DESIGN.md` (2026-10-14, the review packet; its §12 holds the verdicts), `_REVIEW_PROMPT.md`,
  `_REVIEW.md` (the external review: build the competing-placement architecture, not the first proposal),
  `_DERIVATION.md` (2026-10-14, superseded where it conflicts), `_IMPLEMENTATION.md` (2026-10-15, "not ready for
  production implementation") and `_PLAN_REVIEW.md` (2026-10-15, superseded in part).
- `my_notes.md`, the owner's comments of 2026-10-16.

Two reviews they cite: `git show f763270c:docs/dev/INTERGENIC_GDNA_REVIEW.md` (2026-10-15: the one-mate fragments,
since repaired, and the normalised model `DESIGN.md` §3.1g partly realizes) and `PARKED_WORK.md`.

## 10. Starting the session

> We are taking `ISSUES: splicing-artifacts` for Rigel 0.8.0, the next step in the tool's order. Read `CLAUDE.md`,
> then `docs/dev/SPLICE_ARTIFACTS.md`, then the ISSUES entry and `DESIGN.md` §1.0, §1.1, §3.1d and §3.1g. The scope is
> §1.0: Rigel reads alignments and never repairs them, so the contiguous reading of a spliced fragment must come from
> the aligner, and the EM weighs it as a place. Nothing has been measured since that scope, the annotated-path rule and
> the place model landed. First, without changing code, measure on the STAR rung, against the read names, every gDNA
> template that reaches RNA through a spliced hit and every RNA template that lands in gDNA: by class, by whether the
> aligner reported a contiguous alternative, and by where the EM put it, per stratum, pinned and fractional (question 1
> of the document's section 7). Then derive on one page, the multi-junction case included, what is left for a mechanism
> to do within §1.0, and bring me the measurement and the derivation before any prototype.
