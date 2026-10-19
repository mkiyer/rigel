# The annotated-path rule: design

2026-10-19 · code base `497d6f45` · approved by the owner and implemented; §10 records how the build departed from it. Every claim about today's code
cites a line re-read for this document.

## 1. The change in brief

A spliced hit is usable only if one annotated transcript holds every junction it observes (plan §8). The rule joins the
pair validator in the scan's one hit verdict, decided by one function. A hit it rejects is withheld with its reason; the
template keeps its other hits; a template left with none is withheld, as an unpaired one is today. After it, every
spliced fragment any stage sees has its junctions held by one transcript, so the code that exists for any other spliced
fragment is unreachable, and a later stage deletes it (§6).

**The scale.** On the STAR panel (recommended parameters, today's main) the resolver labels 0.02–0.06 % of each true
origin `spliced_unannot`: about 1,900–2,500 gDNA templates per g50 condition, and 1,000–4,400 RNA ones, which are
STAR's misreadings of annotated isoforms. The rule also withholds hits labelled annotated whose junctions no one
transcript holds, which the label hides. The panel's effect will be small; VCaP's RNA library carries more (plan §8:
about 0.8 % of spliced RNA junction units).

## 2. What it fixes

Today a spliced hit is labelled annotated as soon as one of its junctions matches, and its candidates come from
`merge_sets`, which falls back to the union of the junctions' transcript sets when no transcript holds them all
(`resolve_context.h:1306-1358`, `constants.h:291-338`). Only a hit with no matching junction is `SPLICED_UNANNOT`
(`resolve_context.h:1390-1402`). The stages treat those two kinds of hit four ways:

| stage | treatment | where |
|---|---|---|
| the EM | RNA only | `scoring.cpp:223-234` |
| calibration | an intron with no sj id is no splice: the fragment is banked unspliced, gDNA-eligible | `accumulator.cpp:641-653` |
| the outputs | with no candidate, its own intergenic RNA pool and ZF 0x21 | `bam_scanner.cpp:2045-2052, 2186-2197, 2738-2752`; `cli.py:350-382` |
| strand training | left out: it needs `SPLICED_ANNOT` | `resolve_context.h:109-115` |

`ISSUES: splicing-artifacts` names the three fragment types on which calibration and the EM disagree: an unannotated
sequenced junction; two annotated junctions on opposite motif strands; mates disagreeing on an acceptor. The rule
withholds all three (§3), so the disagreement ends with them. It also takes out of the RNA length and strand training
the hits whose candidate lengths count a novel junction's skipped bases as sequenced (plan §8).

## 3. The rule

**A hit's observed junctions** are the kept CIGAR junctions of every record of the hit, deduplicated by reference,
coordinates and motif strand: the set `build_fragment` hands the resolver (`bam_scanner.cpp:773-863`). Every record
includes the chimeric supplementary records that the recommended STAR parameters write (`chimOutType WithinBAM`). A hit
whose junctions were all blacklisted observes none and is unspliced, as today: the blacklist stays the one declared
exception (plan §8).

**A transcript holds a junction** when it has that intron exactly, on the motif's strand when the motif is definite and
on either strand when there is none: `sj_lookup_into` (`resolve_context.h:1116-1150`), unchanged.

**A spliced hit is usable** when one transcript holds every observed junction; those transcripts are its **holders**.
Otherwise the first reason that applies:

1. **contradictory:** the junctions cannot be one molecule's as observed: two overlap without being equal, or their
   motifs name both strands (a record whose XS and ts disagree, mates reading one junction on opposite motifs, two
   junctions on opposite strands);
2. **novel junction:** some junction is held by no transcript;
3. **unsupported combination:** each junction is held by some transcript, and none holds them all.

**The order of the checks.** Missing mate, reference, strand, outward, length, then the rule. So a hit the rule
withholds is always a molecule with a fragment length: the population the owner's gene-level salvage would take
(`ISSUES: fragments-that-fit-one-gene-but-no-transcript`). Until then it is withheld.

**What it leaves alone.** A usable hit's candidates and lengths are today's. With every junction held, `merge_sets`
returns the holders, and the candidates are the holders the hit's blocks overlap (`resolve_context.h:1344-1358`). The
plan's "the witness limits neither the candidates nor their lengths" is true of the validator's witness; the candidates
were already drawn from the holders. The template keeps its other hits: a gDNA template reported both as an unannotated
splice and as a contiguous alignment keeps the contiguous one. Rigel constructs no alternative reading.

## 4. The code: one junction set, one predicate, one verdict

- **`hit_junctions(hit, out)`** (`bam_scanner.cpp`): every record's kept junctions, sorted and unique. The validator,
  the rule and `build_fragment` all read it. Today there are two collections: `hit_verdict` reads the two primary
  records only (`bam_scanner.cpp:1007-1021`), `build_fragment` every record, into its own `std::set`. With two, a
  chimeric hit could take its length from a witness that holds the primaries' junctions and pass the rule on holders of
  them all: the hybrid that the targeted review's finding 1 refused.
- **`FragmentResolver::junction_holders(junctions, scratch, holders)`**: the holders, and `HIT_USABLE` or the reason
  there are none. It is the one definition of "one transcript holds every observed junction". The validator's witness
  search calls it in place of its own loop (`resolve_context.h:796-813`).
- **`hit_verdict`** runs today's checks unchanged (they skip the hits they do not judge: a mate on several records, a
  mate off the index), then the rule, on every hit. The verdict codes become the hit's: `HIT_USABLE`,
  `HIT_MISSING_MATE` … `HIT_LENGTH`, and `HIT_CONTRADICTORY`, `HIT_NOVEL_JUNCTION`, `HIT_UNSUPPORTED`, since a rejected
  hit is no longer only a pair that is not a molecule.
- **The annotation pass** replays `hit_verdict` (`bam_scanner.cpp:2706-2712`), so it applies the rule with no second
  copy.
- **The resolver is untouched.** It labels every hit the scan passes it `SPLICED_ANNOT`: every junction matches, and
  the block that ends at a held junction ends on the holder's exon, so the hit overlaps an exon
  (`resolve_context.h:1344-1358`).

Nothing new is stored on the buffered fragment.

## 5. Outcomes and counters

- **Per hit:** `hits_rejected`, by its eight reasons (`pairs_rejected` today).
- **Per template:** a template with no usable hit is withheld as **unannotated splice** when the rule withheld one of its
  hits, and as **unpaired** otherwise (the plan's §10 order). The two counts replace `withheld_no_valid_pair`.
- **The annotated BAM:** a withheld template's records keep ZF 0, and ZC carries the outcome, `unannotated_splice` or
  `unpaired`.
- **summary.json:** the renamed and new keys, with its schema version bumped.
- **The salvage census** (how many withheld molecules lie within one gene, between genes, or outside the genes) is
  computed by the measurement from the annotated BAM and the index's gene table, not by the scan (decision 2).

## 6. Stages

Each stage lands alone, after its gate.

0. **The rename:** `PAIR_*` to `HIT_*`, `pairs_rejected` to `hits_rejected`. No number moves (`rename_identity.py
   --check`).
1. **One junction set and `junction_holders`.** The validator reads every record's junctions. A verdict changes only for
   a hit with a spliced supplementary record whose genome path fails; the panel's verdict counts and buffer show
   whether any did.
2. **The rule**, its counters and the annotated BAM's outcome, with one more gate: the resolver's junction table
   (`index.py:1600-1619`) and calibration's (`splice_graph.build_sj_arrays`) are built along two paths, and no test
   holds them equal. Calibration finds an sj id for every observed junction of a usable hit only if they are.
   Measured (§8).
3. **Deletions**, each its own commit and a numeric no-op on the panel and the ladder:
   1. **The intergenic RNA pool.** A usable spliced hit always has candidates (§4), so `n_intergenic_spliced` is
      always 0. It goes, with `count_intergenic_rna`, `fraction_intergenic_rna`, `quant_accuracy.py`'s split
      (`:329-357`) and the report's segment. The annotation pass's intergenic-RNA stamp (ZF 0x21,
      `bam_scanner.cpp:2738-2752`) is then reached only by a gated spliced fragment, so it goes with the fix for
      `gated-fragments-stamped-intergenic`.
   2. **`SPLICED_UNANNOT`.** The scan never produces it. The splice type, its two count columns, its census key, its ZS
      value, `n_with_unannotated_sj` and the EM's arms for it (`scoring.cpp:229, 512-518`) go. The resolver's junction
      block reads `junction_holders`, and refuses a spliced fragment with no holder: a precondition the scan
      guarantees. Tests that hand the resolver an unannotated junction directly (`test_implicit_splice.py:259-269,
      411`) move to the scan. The splice-type codes are renumbered, which changes the buffer and the caches, so this
      lands before the scan phase's one re-cache.
   3. **Pass 2's motif branches.** An ambiguous motif on a placement's observed junctions, and a path contained in no
      object (`second_pass.py:435, 466`), are unreachable for a usable hit. They go with the plan's motif-rule
      correction (§9.3, item 4), which this rule simplifies.

**Not this change.** The deposit's own ambiguous-motif guard (`accumulator.cpp:641-643`) still serves implied introns:
when a fragment carries no motif, the deposit takes the chosen hypothesis's strand (`accumulator.cpp:619-620`), which is
ambiguous when opposite-strand transcripts imply one path (`resolve_context.h:1025`). Merging the two junction tables
into one is a refactor for later; after stage 2's gate they agree.

## 7. Tests

Each fails on today's code before the fix; each break below must fail at least one.

| fixture | expected |
|---|---|
| a novel junction in gene A, with a contiguous alternative under another HI | the spliced hit `novel_junction`; the contiguous hit kept; `num_hits` 1 |
| an annotated and a novel junction in one hit (labelled annotated today) | `novel_junction` |
| two annotated junctions that no transcript joins (labelled annotated today) | `unsupported` |
| mates' junctions overlapping without being equal | `contradictory` |
| mates reading one junction on opposite motifs; a record whose XS and ts disagree | `contradictory` |
| an annotated junction read on the opposite definite motif | `novel_junction` |
| an annotated junction with no motif tag | usable |
| a junction on a chimeric supplementary record that no transcript holds with the primaries' | rejected by the rule: every record is read |
| genome path too long; the primaries' junction held by T1, the supplementary's only by T2 | `length`: T1 is no witness |
| a mate on two records (the validator does not judge it), with a novel junction | `novel_junction` |
| a spliced hit outside every gene, and one on a reference off the index | `novel_junction`; the template `unannotated_splice`; nothing reaches the intergenic counters (rewrites `test_annotate.py:930-936`) |
| a hit both too long and novel | `length` |
| one hit `length`, another `novel_junction`; every hit `length` | `unannotated_splice`; `unpaired` |
| the only junction blacklisted; one of two blacklisted | unspliced and usable; decided by the kept one |
| a usable spliced hit | resolved exactly as today |
| the index | the resolver's and calibration's junction tables hold the same junctions |

Breaks: withhold the whole template when one hit fails; check each junction alone (the unsupported combination must
fail); ignore the motif strand; drop the overlap test; drop the motif test; read the primary records only; run the rule
before the length; invert the outcome order; let the validator alone read the primaries.

## 8. Measurement

- **The STAR panel** at the recommended parameters, scan pinned and fractional, against today's main (the num_hits
  arm): per stratum, transcript and gene error and the pools; the census by true origin and splice type; the rejected
  hits by reason and true origin; the withheld templates by outcome and true origin. The gDNA withheld is the gain. The
  RNA withheld is the cost: here STAR's misreadings of annotated isoforms, since the simulator writes no other.
- **The salvage census** from the annotated BAM: withheld molecules within one gene, between genes and outside the
  genes, by true origin.
- **The oracle ladder:** identical, since every simulated spliced read follows its transcript with its strand's XS. The
  new reasons must read zero there.
- **A regression** is dissected by gene before it is attributed to the rule (the paralog fragility,
  `ISSUES: the-em-stops-before-it-converges`).

## 9. Decisions for the owner

1. **Contradictory includes motifs that disagree** (§3). The plan's reason covers overlapping junctions only. Without
   the extension, a hit whose record carries disagreeing XS and ts stays usable: the EM certifies it RNA while
   calibration, reading its motif as ambiguous, banks it unspliced (`accumulator.cpp:641-653`), which is one of the
   disagreements the rule exists to end. Mates on opposite motifs and opposite-strand junctions are withheld either
   way; only their reason changes. On the STAR panel it should read zero: under `intronMotif`, STAR filters out
   alignments whose introns disagree on the strand. Recommended.
2. **The salvage census is measured, not built into the scan.** Plan §8 reports the withheld templates in three classes.
   Computing them in the scan needs each gene's span and each transcript's gene in the resolver, state nothing else
   there reads. The measurement computes them from the annotated BAM and the gene table, and the scan gains them with
   the salvage and the gene-level EM it needs. The measurement reads the class off the molecule's span (both mates
   inside one gene: the owner's criterion), not off each side of the junction as plan §8 put it; the two differ only
   for a molecule that reaches past a gene's end. Recommended.
3. **Stage 3 changes the output formats.** summary.json loses the intergenic RNA count and fraction and the
   unannotated-splice keys; the annotated BAM loses ZF 0x21 and ZS `spliced_unannot`. Recommended, as its own stage
   after the rule is measured.

## 10. Implementation notes (2026-10-19, after the owner's approval)

- **Stage 0** (`2b5bc0d1`) and **stage 1** (`69d003ea`) landed: the rename is bit-identical; one junction set per
  hit left the STAR panel (15/15 rows) and the oracle ladder (48/48) identical.
- **Stage 2** (`aff4be8e`), gated by 15 tests and 10 breaks that each fire, landed by the owner's ruling despite its
  stranded losses (transcript error −10.17 % unstranded capture-OFF, +4.55 % stranded capture-OFF, +1.28 % stranded
  capture-ON), which come through each condition's one giant EM component and calibration's gDNA prior:
  `ISSUES: the-annotated-path-rule-costs-the-stranded-strata`. Numbers and dissection in
  `$RIGEL_RUNS/arms/2026-10-19_annotated_path/` (`report_s2.md`, `dissect/`, `summary.md`).
- **Stage 3a** (`aaa0a93d`): the intergenic RNA pool deleted, a no-op on the STAR panel and the ladder.
- **A junction with no aligned base.** A record whose CIGAR has an N and no M (`5S40N5S`) carries a junction with no
  anchor. With that junction annotated, the hit is usable and can still resolve to no candidate, so
  `n_intergenic_spliced` is not zero by construction; the intergenic RNA pool's deletion books such a fragment as
  plain intergenic. STAR writes no such record.
- **Stage 3, revised.** The resolver still meets junctions no transcript holds: the pairing without HI
  (`pair_multimapper_reads`) resolves each mate's location to pair them before any hit verdict. So `SPLICED_UNANNOT`
  stays until the plan's geometric pairing replaces that path; the intergenic RNA pool can go now.
- **The test chromosome's shadow genes.** U5 and U6 are spliced and unannotated, so the scan withholds their spliced
  reads; their unspliced reads remain the control (`TESTING.md` §0a, updated). The scan caches show it after the
  scan phase's one re-cache.
- **The two junction tables agree** on the panel's index (13,482 junctions each), as the new gate requires of the toy
  indexes.
