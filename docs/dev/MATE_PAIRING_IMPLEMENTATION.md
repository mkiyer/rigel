# Mate pairing: implementation map (2026-10-15)

Companion to `MATE_PAIRING_FIX.md`, which is the contract. This file maps the contract onto the code, records the
implementation-level choices made inside it, and tracks the stages. Work happens in the worktree
`<worktree>` (from 41899ae4); runs and evidence go in `$RIGEL_RUNS/prototypes/2026-10-15_pairing/`.
Sandbox; nothing here is in `src/`.

## 1. Gate 1: the baseline and the raw-input census (done)

- **The baseline** is a wheel built from the unmodified reviewed commit, the same way the prototype will be built,
  so compile flags cannot masquerade as numerical change. `MANIFEST_baseline.txt` holds the commit, the native
  modules' sha256, the compiler and `pip freeze`.
- **The census** (`census/input_census.py`, full files, units = templates; `census/CENSUS.md`):

| input | templates | primary pairs | one-mate primary | 5′ soft clip (pair) | SA chain | secondary | HI truncation |
|---|---|---|---|---|---|---|---|
| oracle `g50 ss.99 ON` | 10.0 M | 100 % | 0 | 0 | 0 | 0 | — (no HI) |
| minimap2 `g50 ss.99 OFF` | 10.0 M | 89.0 % | 1.05 % | 0.66 % | 0.005 % | 13.6 % | — (no HI) |
| VCaP (STAR 2.7.11b) | 18.7 M | 99.3 % | 0.002 % | 4.21 % | 0.36 % | 0.80 % | none |

  - Every input satisfies the reporting contract: no NH/HI truncation, no inconsistent NH, no ambiguous mate role.
    Primaries carry no 5′ hard clip anywhere.
  - VCaP's 5′ clip lengths: 78 % are 1–5 bases, 10 % 6–20, 12 % over 20.
  - **5′ provenance, declared:** VCaP's STAR command line shows de-duplicated, prefiltered FASTQs and no 5′
    trimming step. The analysis assumes no non-template 5′ sequence, to be confirmed by the owner. The simulated
    panels have no clips (oracle) or clips from simulated errors only (minimap2).

## 2. The representation

- **`MateLocation`** (new, `bam_scanner.cpp`): one mate's self-contained alignment.
  - reference and strand;
  - aligned blocks, split at N only, as `parse_cigar` does today;
  - **raw** skips, kept before the blacklist removes any, so a hidden contradiction stays visible;
  - the blacklist-surviving skips;
  - terminal clips in read orientation, soft and hard apart, with adjacent H/S consumed;
  - the SA tag's presence;
  - NM, the junction-strand tag, provenance (primary, secondary, HI) and the raw record ids it stands for.

  Supplementary records never become locations. Identical locations (same reference, strand, blocks, skips and
  clips) are de-duplicated into one location with several record aliases.
- **`FragmentAssembly`**:
  - the outcome: placed, unpaired, reported chimeric or excluded;
  - the rejection reasons;
  - the placements in canonical order, each an R1 location, an R2 location and its validated paths;
  - record → placement membership;
  - the reporting status.

  The scanner and the annotated-BAM writer both call one `assemble_fragment()`.

## 3. Validation of a candidate pair

Every R1 location is tried with every R2 location, HI or not (contract §"Normalize and enumerate"). A pair becomes a
placement only if it passes all of the following.

1. **One reference, opposite strands.** Call the forward mate F and the reverse mate R.
2. **Endpoints** (the conditional clip model), per candidate path P:
   - the molecule start is `u − cF` along P, where u is F's aligned 5′ (left) boundary;
   - the molecule end is `v + cR` along P, where v is R's aligned 5′ (right) boundary;
   - `cF` and `cR` are the 5′ soft-clip lengths;
   - a 5′ hard clip, or an SA chain, withholds the location.
3. **Containment.** On each surviving path, both mates' aligned blocks lie inside `[start_P, end_P)`. This refuses
   outward and protruding pairs and admits a genuine overlap.
4. **Mate-against-mate splice consistency**, on the raw skips. No aligned base of one mate lies inside the other
   mate's skip by more than K, and skips that overlap must be identical.
5. **Paths:**
   - **genomic:** L = end − start in genomic coordinates. It exists only when no blacklist-surviving skip remains,
     which keeps the legacy ARTIFACT reading.
   - **observed splices:** the observed skips removed, the mate gap contiguous. This is how an unannotated splice
     forms an RNA placement.
   - **one path per candidate transcript t:** the mates' observed skips plus t's introns in the unsequenced gap,
     with the clip extension walked along t's exons.
6. **The bound.** Each path must satisfy `1 ≤ L ≤ scan.max_frag_length`. A pair with no surviving path is no
   placement.

## 4. Implementation choices made inside the contract

Flag any of these and it changes before it is built.

1. **Overhang stays as it is.** A transcript whose intron overlaps one read's aligned bases by more than K remains an
   EM candidate with today's overhang penalty, and keeps today's length rule.
   - That is a per-read emission policy, not pairing. Calibration's gap hypotheses already exclude such a transcript.
   - The pairing contract's "remove invalid interpretations" is applied here to pair-level contradictions, exhausted
     paths and the length bound.
   - Changing overhang would move every library, the oracle included. It is its own mechanism and gets a parked
     entry.
2. **K for mate-against-mate consistency is the existing `scan.splicing_anchor_tolerance` (3 bp).** The resolver
   uses it for the same physical ambiguity, an aligner running a read a base or two past a donor. Strict equality
   would refuse pairs the resolver accepts.
3. **The clip extension along a transcript** walks t's exons from the aligned 5′ boundary when that boundary lies
   in an exon of t.
   - If the boundary is outside t's exons (an overhanging end), the extension continues contiguously, as the
     existing overhang convention measures it.
   - A walk past t's first or last base exhausts the path, and t is no candidate for this pair.
4. **Extents per hypothesis.** `GapHypothesis` and the deferred bank gain a start and an end per hypothesis.
   - A clip that crosses an annotated intron becomes an implied intron of that hypothesis, never an observed one.
   - Without clips, every hypothesis keeps the placement's extent, so deposits are bit-identical.
5. **The EM buffer:**
   - `frag_lengths` carries each candidate's validated L.
   - `genomic_footprint` carries the genomic path's L, or −1 when that path is invalid. The scorer then gives the
     unit no gDNA reading. Today a length ≤ 0 scores as neutral; the change goes through an explicit rule with a
     guard.
   - A new per-candidate molecule-start column feeds `coverage_weight`, so `genomic_start` keeps its observed
     meaning.
   - `num_hits` uses the existing checked 16-bit narrowing; today the cast is silent.
6. **Reporting status** is computed on the raw read-name group, before Rigel's filters. With HI tags, NH above the
   number of distinct HI values present means truncation, and the run stops with the read name, the counts and the
   STAR remedy.
7. **Supplementary records** follow the approved chain policy: a location whose record carries SA is
   chain-dependent and withheld. On VCaP that is about 0.36 % of templates, priced exactly at stage 1.

Made while building stages 2–3 (each is open to a flag):

8. **The witness paths are the resolver's own gap hypotheses**, plus the observed path (the mates' skips cut, the mate
   gap contiguous). Each is measured by `Accumulator::hypothesis_length`, now a public static member: the one
   definition of L that calibration deposits by. A pair is a placement if its shortest path is `1 ≤ L ≤ max`.
   - The resolver's candidates that do not contradict the blocks supply the transcript paths: their introns inside
     the unsequenced gaps, within ±K.
   - The observed path is never shorter than a transcript path, so it decides only when no transcript implies
     anything.
9. **A validated pair is resolved as one molecule**: `_resolve_core(..., one_molecule=true)` skips `detect_chimera`.
   - The only pairs this changes are mates joined by an unannotated splice between transcript-disjoint genes with
     an observed path within the limit. They were `CIS_STRAND_SAME` chimeras; now they are placements.
   - None exist on the oracle.
10. **A mate whose CIGAR aligns no reference base** (only clips, or an N between clips) forms no location. Its
    template is unpaired, reason missing mate. The baseline placed such a pair on the other mate's blocks alone.
11. **Duplicates.** Identical geometry is one location: the same reference, strand, blocks, raw skips and kept
    skips.
    - The representative is the report with the least (NM, motif strands), so the choice is order-free.
    - A disagreement is counted (`n_duplicate_conflicts`).
    - Locations sort canonically, and placements follow in (R1, R2) order.
12. **An asserted pair that fails** is a reported chimera if the resolver, run normally, calls it one; otherwise it is
    nothing. A constructed pair that fails is never reported.
13. **Supplementary records stay attached to their primary's location until stage 4.** A location whose records lie
    on two references or two strands fails as a chain; stage 4 replaces this with the SA withholding.
14. **Unpaired reasons**, summing to `n_unpaired`:
    - `n_unpaired_missing_mate`: no location for one mate;
    - `n_unpaired_too_long`: some pair passed the geometry checks but has no path within the limit;
    - `n_unpaired_geometry`: every other case.

    Constructed placements are counted as `n_placements_constructed`.

## 5. Stages (diagnostic arms, one worktree commit each)

| stage | what changes | expected on clip-free, secondary-free input (the oracle) |
|---|---|---|
| 1 | complete pairs only; accepted-placement count; the unpaired outcome; one assembler for scanner and writer; reporting status | bit-identical |
| 2 | exhaustive R1 × R2 construction with the full checks, for constructed pairs | bit-identical |
| 3 | the same checks on aligner-associated pairs (the primary pair, HI groups) | bit-identical (every oracle pair is proper, clip-free and consistent) |
| 4 | 5′ soft clips: path-specific endpoints, occupancy, lengths and coverage starts, compared with blanket withholding | bit-identical (no clips) |
| 5 | the hard bound on every retained path, in the EM as in calibration | **not** identical: about 1.6 % of fragments, implicit splices longer than 1,000 bp, lose their gDNA reading |

Then reporting, persistence and capacity: stats, `summary.json`, the report, `MANUAL.md`, the aligned scenarios'
accounting, the scan-cache assembly digest, and the 16-bit audit. Then the final same-input validation
(`MATE_PAIRING_FIX.md`).

## 6. Tests

- **Scanner-observable gates** (`tests/test_paired_placements.py`): the reviewer's six cases, the geometry probes
  and the counting, ledger and annotation cases. They fail on the baseline where the defect shows.
- **The independent reference** (`tests/native/_pairing_reference.py`) enumerates placements, paths and lengths from
  CIGARs and a tiny annotation, without the resolver. It is compared, case by case, with a native
  `assemble_fragment` exposed for testing.
- **Mutation gates:** each repair is broken back in, one at a time (contract §"Falsification and invariant tests"),
  with the failure saved as evidence.

## 7. Stage records

Worktree `<worktree>`. Cumulative patches are in `patches/` of the run directory:
- `stage1.patch`;
- `stage1b_seam.patch`: the test seam `BamScanner._assemble_bam` and the reference;
- `stage2.patch`;
- `stage3.patch`.

Census scripts are in `census/`, logs in `logs/`, and transition tables in `transitions/`.

**Gates.** The enumeration gate (`test_the_placements_are_exactly_the_specified_pairs`, four provenances × 33 cases)
compares every read name's accepted set, extent, shortest path, outcome and EM lengths with the reference.

| tree | pairing gates | enumeration (primary / hi / secondary / cross_hi) | full suite |
|---|---|---|---|
| stage 1b | the geometry and clip gates fail | all four fail | — |
| stage 2 | asserted geometry, clip and ledger gates fail | secondary passes; the other three fail on the aligner's own pair only | 2,971 passed, 13 failed (12 later-stage gates; one fixture, rewritten at stage 3) |
| stage 3 | only the five clip/chain gates (stage 4) fail | all pass | 2,978 passed, 7 failed (5 stage-4 gates; 2 fixtures, fixed below); with the fixtures fixed and the gates added, **2,982 passed, 5 failed** (the stage-4 gates) |

**Fixtures changed at stage 3,** each a pair the contract no longer places:
- `test_an_intron_on_another_reference_is_not_cut` built an R2 that aligns no base. It moved to the pairing gates as
  `test_a_mate_that_aligns_no_base_is_no_mate`. The adapter's other-reference intron filter is now unreachable from
  the scanner; deleting it is a follow-up.
- `test_scan_worker_exception`'s 70-kb fragment exceeds the default limit, so the test raises `max_frag_length` to
  keep its throw site.
- The annotate sentinel moved both mates to one start. A spliced donor then protrudes, so the pair now moves whole.

**Mutation sweep** (`logs/mutations/SUMMARY.md`): each of 11 repairs broken back in alone. Every one fires at
least one gate, and the control fires only the five stage-4 gates. The sweep found one hole and it was closed: a
chimeric resolve leaves `frag_lengths` empty, which the EM-length check zipped past. The seam now reports
`chimera_type`, and an accepted placement must carry no chimera type and one length per candidate.

**Oracle identity** (`rename_identity.py`, references frozen from the baseline wheel in `identity/`): stages 2 and 3
are BIT-IDENTICAL on `g00 ss.99 OFF`, and stage 3 on `g50 ss.99 ON`. That covers the transcript table and every
payload and calibration array.

**Transition census** (`census/transition_dump.py`, scan and drain, joined on frag_id).

| | minimap2 g50 ss.99 OFF (10.0 M templates) | VCaP (18.7 M) |
|---|---|---|
| baseline → stage 2: deposits at L = 100 | 243,990 → 23,456 | 722 → 719 |
| unpaired, missing mate (= the census's one-mate primaries) | 104,886 | 317 |
| multimapper groups | 975,299 → 627,801 | 136,188 → 125,279 |
| constructed placements (no-HI pairs; cross-HI combinations) | 2,581,196 | 11,207 |
| stage 2 → 3: newly unpaired, geometry / too long | 2,252 / 9,803 | 33,938 / 13,158 |
| stage 2 → 3: `CIS_STRAND_SAME` chimeras now placements | 12 | 3,544 |
| unique fragments whose footprint changed | 0 | 1 |

**By stage** (`transitions/*_base_to_proto.md`, `*_proto_to_s2.md`).
- **On minimap2, stage 1 is almost all of the length repair.** Deposits at L = 100 go 243,990 → 25,995 once a lone
  mate stops being a fragment. Stage 2 adds 25,995 → 23,456.
- **Stage 2's 1,847 new reported chimeras on minimap2** are fragments whose only placements were geometrically
  invalid combinations Rigel had built beside a chimeric aligner pair.
- **On VCaP, stage 2 turns 6,860 multimappers unique by de-duplication.** STAR writes one pair alignment twice
  (NH 2), differing only in where an insertion sits in a homopolymer (`109M2I39M` against `108M2I40M`): the same
  reference blocks and the same NM.
  - This confirms choice 11's key, the geometry the model reads rather than the CIGAR string.

**Who stage 3 unpairs** (`census/unpaired_breakdown.py`). Each template's primary pair is judged by the
specification; the tables are `transitions/*_s3_unpaired_breakdown.md`.

| VCaP, 46,875 templates | n | at stage 4 |
|---|---|---|
| primary carries SA (a supplementary chain) | 25,463 | stays unpaired: the approved chain policy |
| no path within 1,000 (spliced 7,928, unspliced 3,998) | 11,926 | stays unpaired: the hard bound |
| outward, unclipped (back-splice candidates among them) | 3,741 | stays unpaired |
| protruding only through 5′ clips; a contiguous clip walk contains it | 2,669 | restored by the clip model |
| splice conflict | 1,747 | stays unpaired |
| unclipped protrusion, same strand, other | ~1.3 k | stays unpaired |

On minimap2 every template is a true simulated molecule of at most 500 bp, so each of its 12,071 is a real molecule
the aligner placed impossibly. The examples inspected:
- one mate's primary on the wrong locus, 9–25 Mb away, with no reported alternative near the other mate;
- a spurious short-anchor junction (`5M1636N95M`) protruding past the other mate's 5′ end;
- a spurious 26-bp skip the other mate reads through;
- 51 bases placed in the true transcript's intron at NM 0, joined by an unannotated junction.

The baseline placed them on absurd lengths (one 25-Mb "fragment"). Recovering such molecules from alternative
alignments is the parked splice-artifact design's question, not this repair's.

**Carried forward:**
- **The old `mate_unmapped` unit.** It also counted the deleted pairing code's lone-mate leftovers: minimap2 g50
  595,634 → 516,346. Redefine it in the reporting stage; never change it silently.
- **Stage 5 must take each candidate's EM length from its validated path.** Today `tx_frag_length` ignores an
  unannotated observed skip. The pairs validated across two genes (VCaP 3,544) carry candidates 11 kb long, and on
  their observed path they are about 300 bp.
- **A placement whose only path in bound is an unannotated splice** keeps no EM candidate once stage 5 drops the
  unsupported interpretations. The owner has clarified its reporting/abundance scope; implement
  [the no-supported-transcript contract](MATE_PAIRING_FIX.md#spliced-placements-without-a-supported-transcript).
  Do not repair an incompatible target by assigning it the observed path's shorter length.
- **VCaP's 110 duplicate conflicts:** identical STAR alignments repeated with a different NM or strand tag. Not yet
  inspected.

## 8. Stage 4 design: the 5′ clip endpoint model, walked along each path

Written before building stage 4. It refines choice 3 above. The observed path, not only a transcript path, needs a
walk.

**The molecule on a path.** On path P, the molecule is `[u_P, v_P)`:
- `u_P` lies `cF` path bases before F's aligned 5′ boundary;
- `v_P` lies `cR` path bases after R's;
- a path base is one that P keeps.

Containment (both mates' aligned blocks inside `[u_P, v_P)`), L and the length bound are judged per path. The
observed alignment is untouched: blocks, NM and observed junctions stay as reported, and a clipped base is
neither coverage nor junction support.

**The walks.**
- **The observed path** keeps every base except the mates' kept skips, so its walk jumps a kept skip it meets.
  - Example: F is `5S95M` starting where R's skip `[1600, 2000)` ends.
  - F's five clipped bases are `[1595, 1600)`, so the pair is contained on this path, although a contiguous
    extension to 1,995 would put R's `[1597, 1600)` block outside the molecule.
- **A transcript path** keeps t's exons. From an anchor inside an exon of t, the walk follows t's exons across its
  introns, and each crossed intron joins that hypothesis's implied introns. An anchor outside t's exons (an
  overhang, choice 1) walks contiguously.
- **Exhaustion.** A walk past t's first or last base, or past either end of the reference, leaves that path with
  no molecule.
- **The genomic (gDNA) path** is the observed path of a pair with no kept skip: a contiguous walk.

**Where it is computed.** `_resolve_core` takes the two anchors and their clips. The resolver owns the transcript
exon index the walks need.
- A hypothesis becomes (implied introns, `u`, `v`); candidates implying the same triple are one hypothesis.
- `frag_lengths[t]` becomes the length of t's walked path (`tx_frag_length` over `[u_t, v_t)`). An overhang
  candidate keeps today's rule over the observed path's extent.
- With both clips zero, every walk is empty, every extent is the aligned extent, and every output is today's. That
  is the bit-identity the oracle checks.

**Where the extent goes.**
- `GapHypothesis` gains `start` and `end`. `hypothesis_length`, the deposit's occupancy, the deferred bank and
  `length_under` read them.
- In the EM buffer:
  - `genomic_footprint` becomes the genomic path's L;
  - a new per-candidate `molecule_start` column feeds `coverage_weight`;
  - `genomic_start` keeps its observed meaning.

**Containment before or after the resolver.**
- An unclipped pair's endpoints are the same on every path, so its containment stays the exact pre-resolve check.
- A clipped pair passes the clip-independent checks (one reference, opposite strands, splice consistency) before
  the resolver, and containment per path after it.
- A pair none of whose paths is contained is rejected for geometry; one whose contained paths are all too long is
  rejected as too long.

**Withheld, with reasons.**
- A location whose record carries SA is a chain.
- A location with a 5′ hard clip has an unknown endpoint.
- A supplementary record is never a location.
- Each is counted, and the unpaired reasons gain both.

**The comparison arm.** Blanket withholding of every 5′-soft-clipped location, on the same inputs. The census reports
exactly which fragments and paths return under the model.

**The reference gains clips.** `_pairing_reference.py` walks the same paths base by base, and the enumeration gains
clip cases: both orientations, crossed anchors, a walk across an annotated intron, an observed-path walk across a
kept skip, exhaustion at a transcript end, and the bound at `max` and `max + 1` after extension.
