# Multimapper intergenic alignments — the dedicated document

Sandbox document for a session devoted to it (owner, 2026-10-27). The open record is
`ISSUES: multimapper-intergenic-alignments` (release-blocking). It lands with `ISSUES: the-gdna-term-is-a-mean-over-places`
and `ISSUES: mismatch-gate-nan-in-multimapper-gdna`. Nothing is in `src/`. Receipts:
`$RIGEL_RUNS/prototypes/2026-10-27_mm_intergenic/` (`mm_census.py`, `assign_by_class.py`, their JSON).

## 1. The goal (owner, 2026-10-27)

**Mitigate the false positives where gDNA is called RNA.** gDNA bleeds into the RNA population, worst at extremely high
gDNA content, and worst of all where a gDNA fragment aligns as spliced. The opposite risk, RNA pulled into gDNA, is
accepted: today the error runs the other way.

## 2. What the data say

Each read pair's records are grouped into places by mate pairing. A place is genic if it overlaps an annotated
transcript's span, and spliced if it holds an `N`. A multimapper with a genic place and an intergenic unspliced one is
the owner's **case 5**. The table shows each class's fate today on the STAR rung (rustar at the recommended
parameters), pinned and fractional, against the read names.

**At g98 ss.99, capture OFF: the bleed.**

| true gDNA, by class | n | called gDNA | called RNA |
|---|---:|---:|---:|
| unique | 8,650,384 | 99.8 % | 16,341 |
| multimapper, all places genic | 165,794 | 91.4 % | 14,039 |
| multimapper, case 5 | 172,256 | 84.1 % | 27,369 |
| multimapper, all places intergenic | 671,492 | 96.9 % | 0 (3.1 % unresolved) |

Of the condition's 172,250 RNA calls, 57,751 are gDNA, a third of them; its true RNA is 199,975 templates. gDNA
multimappers carry 41,408 of the false calls (72 %), case 5 alone 27,369 (47 %), and almost all of them go to the nascent
spans. The opposite error also runs: 29 % of unique mature RNA is called gDNA at g98.

**At g50 ss.99, capture OFF.** gDNA multimappers read 65–68 % gDNA, case 5 (68.2 %) and all-genic (65.2 %) alike,
against 94.9 % for unique gDNA. That is about 57 k fragments, 1.2 % of the condition's gDNA. Mature RNA multimappers
read 97.9–98.6 % annotated.

**How common case 5 is.**

| | gDNA | mature RNA | nascent RNA |
|---|---:|---:|---:|
| STAR rung g50 ss.99 OFF | 88,274 (1.8 %) | 26,251 (0.7 %) | 6,088 (0.6 %) |
| of them, the genic place spliced only (certified RNA today) | 1,654 | 456 | 168 |
| STAR rung g50 ss.99 ON | 80,360 (1.6 %) | 23,243 (0.5 %) | 535 (0.4 %) |
| the VCaP mix | 5,631 DNA (0.1 %) | 31,364 RNA (0.2 %) | — |

On the real library case 5 is mostly RNA: a read inside an exon that also places on an intergenic copy, such as a
processed pseudogene. That is the direction of the accepted risk.

## 3. How the tool handles a multimapper today

- **The scan** (`bam_scanner.cpp`, the per-hit loop): every hit is a placement for calibration (`collect_placement`),
  an intergenic one included. A hit with no transcript candidate is appended to the EM's buffer only for a unique
  mapper, and the scorer then skips it. A multimapper's intergenic hits never reach the EM.
- **The accumulator's second pass** (`DESIGN.md` §4.2): a held multimapper is one fragment with several placements,
  and one is drawn per fragment, weighted by `A_p`, the pass-one UNIQUE traffic at the objects the placement deposits
  on. A placement without unique evidence beside one with some is a hard zero. That rule's named consequence, a silent
  processed pseudogene's gDNA going to its expressed parent, was accepted for now (owner, 2026-10-10), "the mixed class
  to be measured on real libraries before any uniform component is considered". §2 is that measurement.
- **The score** (`scoring.cpp`, `score_mm_alignment` and the unit's flush): the unit's gDNA term is the MEAN of
  `h_G(ℓ) · ½ · nm` over the hits gDNA can explain that have a transcript candidate (`n_cand > 0`).
- **The EM** (`em_solver.cpp`): per locus, normalised mixture weights with `log θ − log L̃` as each component's weight,
  one gDNA component per locus, and a grouped prior in the M-step. A multimapper's transcripts join one locus.

So calibration sees a case-5 fragment's intergenic place, through the draw, and the EM does not. With a spliced genic
place only, gDNA is no candidate and the fragment is certified RNA.

## 4. The owner's framing

- Calibration's accumulator deposits fragments into REGIONS and BOUNDARIES genome-wide.
- **Intergenic regions are not captured.** Under hybrid capture every intergenic region is off target.
- **An intergenic fragment does not belong to a locus.** A locus is an overlapping group of genes and transcripts, and
  intergenic regions are what separate loci.
- **Option 1, one intergenic meta-component**: every intergenic region of every chromosome combined into one, at one
  rate, the global sum of intergenic fragments over the global sum of observable intergenic effective lengths. Every
  intergenic alignment competes at that rate. Probably the simplest.
- **Option 2, every intergenic region modelled**: the accumulator deposits uniquely aligned fragments genome-wide,
  then multimappers in a second pass, so a multimapper's intergenic alignment is deposited exactly as a genic one is.

## 5. What the options share, and where they differ

**Both need:**
- **The EM must see the intergenic places**: the scan buffers a mixed multimapper's intergenic hits, keeping the
  fragment ledger exact.
- **Per-place weights and a sum.** A fragment is observed as its places `S`. It arises from one component at one place
  and the aligner reports the set, so `P(fragment | o) = Σ_{p ∈ S} κ_p · P_o(p, ℓ_p)`, with `κ_p` the place's capture
  weight.
  - For a unique fragment every candidate carries the same `κ_p`, which cancels (`DESIGN.md` §3.6b).
  - For a multimapper it does not. Every candidate at a place, RNA's and gDNA's, carries that place's weight, and gDNA's
    term is the sum over its places, not the mean.
  - An intergenic place has `κ = 1`, by the owner's premise.
  - Weighting gDNA's places alone would read an uncaptured intergenic place up to `κ` times too strong against a
    captured RNA place. The 2026-10-15 sum arm, with no place's efficiency, moved 33.6 k RNA multimappers to gDNA.
  - NaN-safe, at `--mismatch-alpha 0` too.
- **A rate no locus parameter carries.** An intergenic place is in no locus, so its gDNA term needs an EM element that
  does not scale with any locus's `θ_G`. Two forms:
  - an extra component whose weight the M-step does not update;
  - a per-unit background likelihood in the E-step that takes its share of the unit and leaves the rest to the locus.

  Either way its mass is reported as intergenic gDNA, outside the loci.

**Option 1.** The rate is calibration's own numbers over the intergenic regions: their deconvolved gDNA counts over
their observable gDNA effective lengths (`count_gdna_region`, `gdna_region_eff_len`). It is one number and matches the
premise. It is blind to how density varies between intergenic regions: partly mappable DNA, which loses a quarter of
its gDNA unplaced around the centromeres (`docs/dev/OBSERVABLE_OPPORTUNITY_DESIGN.md` §7), and repeats.

**Option 2.** Its calibration half largely exists: every intergenic region is a calibration object, and the second
pass already draws intergenic placements. Its EM half would read each intergenic place at its own region's calibrated
rate. It inherits the draw's unique-traffic weighting. A region all of whose fragments multimap has no unique traffic,
so it is drawn rarely or never, its calibrated rate falls toward zero, and its gDNA goes to the genic copy: the
pseudogene-to-parent rule (§3; the region form of `ISSUES: the-second-pass-corrections`' all-held approximation). So it
needs either that weighting repaired or a prior under each region's rate. Calibration's landscape is such a prior: with
intergenic `κ = 1` it puts each intergenic region near the bulk, which is option 1's rate. Option 2 then differs from
option 1 only where a region's own evidence is strong.

## 6. Questions for the session

1. Is "intergenic is uncaptured" a fixed `κ = 1`, or read from calibration? Today the capture reader reads every
   object's density against the bulk, intergenic regions included.
2. Which places are intergenic: the scanner's rule, a hit with no transcript candidate, which leaves intronic places
   genic?
3. Is the intergenic rate fixed from calibration, or re-estimated across loci after their EMs?
4. Does the second pass's draw read the same intergenic rate, so calibration and the EM agree on where a case-5
   fragment's gDNA can sit?
5. Spliced gDNA: case 5 with a spliced genic place only becomes a competition, its intergenic place the gDNA evidence.
   A unique spliced gDNA fragment with no other place is `ISSUES: splicing-artifacts`'.
6. The accounting: the intergenic term's mass joins the summary's intergenic gDNA (`count_gdna − count_gdna_em`), and the
   loci table holds none of it.

## 7. Plan

1. **Derive on paper**: the fixed-rate term in the EM (which of the two forms, its accounting, its interaction with the
   grouped prior), and the second pass's draw if question 4 says so.
2. **Tests first**, failing on main: the owner's five cases on hand-built multimappers; a mixed multimapper's
   intergenic hit reaching the EM with the ledger exact; the sum against an enumeration; a captured RNA place against
   an uncaptured intergenic one; a unique fragment bit-identical; the three hit orders at `--mismatch-alpha 0`.
3. **Prototype**: per-hit `log κ_p` in Python from calibration's efficiencies at the objects each footprint deposits
   on; the scan, the scorer and the EM term in C++ in a worktree.
4. **A/B**:
   - the STAR rung per stratum, with g98 read first (the bleed) and g50 beside it, pinned and fractional, with every
     capture efficiency at 1 beside the shipped run;
   - the per-class fates of §2 beside the pools;
   - the VCaP mix (gDNA against 0.2518, genes against its RNA-only half);
   - the oracle ladder, which holds no multimapper, bit-identical.

## 8. Instruments

- `mm_census.py BAM INDEX OUT.json [vcap]`: case 5's census by true origin.
- `assign_by_class.py BAM INDEX OUTDIR`: today's tree on a BAM with the annotated BAM, then each class's fate by true
  origin. When gDNA wins a multimapper no record carries `ZH = 1`, so the fate is read from any record's `ZF`.
- The STAR rung: each ladder condition's `rustar_recommended.bam`, scored by `quant_accuracy.py --bam`.

## 9. Starting the session

> We are taking `ISSUES: multimapper-intergenic-alignments` for Rigel 0.8.0, a release blocker. Read
> `docs/dev/MULTIMAPPER_INTERGENIC.md` first, then the ISSUES entry, `ISSUES: the-gdna-term-is-a-mean-over-places`,
> `DESIGN.md` §3.1d and §4.2, and the code it cites (`bam_scanner.cpp`'s per-hit loop, `scoring.cpp`'s multimapper
> path, `em_solver.cpp`'s E-step). The goal is to stop gDNA bleeding into RNA. At g98 a third of the RNA calls are gDNA
> and gDNA multimappers carry 72 % of them. RNA pulled into gDNA is the accepted risk. Explore my two options: one
> intergenic meta-component at a global rate, and every intergenic region modelled through the accumulator's second
> pass. Derive the EM's fixed-rate intergenic term and the per-place capture weights on paper first. Answer §6's
> questions with me, then write the failing tests. Prototype in a worktree and A/B on the STAR rung, g98 first, before
> anything enters `src/`.
