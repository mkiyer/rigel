# Chimeric alignments and usable fragments: design (roadmap item 13)

2026-10-20 · code base `09f2aeb0` · STAR 2.7.11b source read for this document · for the owner's review. Every claim
about Rigel's code cites a line re-read for this document, and every claim about STAR cites its source. Every number is
measured; the receipts are in `$RIGEL_RUNS/arms/2026-10-20_chains/`. Approved 2026-10-20 (owner: chimeric alignments
withheld now, the overrun rule its own step, the tail rule deferred); steps 1 and 2 landed the same day, their
measurements the roadmap's item 13 (`PAIRED_END_ROADMAP.md`).

## 1. The recommendation in brief

The owner's criterion: a fragment is usable when its reads can be attributed to a contiguous stretch of the genome, or
to a reference transcript, within `max_frag_length`. Otherwise it cannot be used.

1. **Chimeric alignments are withheld (item 13).**
   - An alignment that is one piece of a chimeric chain is not usable: its record carries `SA`, or its 5′ end is
     hard-clipped.
   - A supplementary record never joins a fragment.
   - A template left with no usable alignment is withheld as `chimeric`, counted.
   - STAR replaces such a read's linear alignment with its chimeric pieces, so the BAM no longer holds a linear
     placement of the read.
2. **Mates that overrun each other are one molecule (a new step).** Two mates whose 5′ ends face each other are a
   molecule spanning both mates, even when one mate's aligned bases run past the other's 5′ end. Such reads are
   attributable to one contiguous stretch, which is the criterion. Mates whose 5′ ends face away from each other stay
   unusable.
3. **Deferred: a chain whose foreign piece is only a read's 3′ tail.** In cfRNA these are likely valid molecules, but
   the only data with truth says the rule would accept wrong ones, so it waits for a simulation that has the mechanism.

## 2. How STAR writes chimeric alignments

Read in STAR 2.7.11b's source:

- **When it looks.** Chimeric detection runs only when chimeric output is on (`--chimSegmentMin` > 0). It then runs
  only when the best linear alignment scores at most the read pair's length less `--chimNonchimScoreDropMin` (20): 180
  for 2×100 (`ReadAlign_chimericDetection.cpp:46-49`). It also runs on overlapping mates merged into one read (§2.1).
- **The linear alignment is dropped.** With `--chimOutType WithinBAM`, a read found chimeric is written as its chimeric
  alignment only; its linear alignment is not output (`ReadAlign_oneRead.cpp:99`). So a chained mate almost never has a
  self-contained alternative in the BAM: 23 of 2,692 in the panel's g50 ss.50 OFF.
- **Two kinds** (`ChimericAlign_chimericBAMoutput.cpp:13-26`):
  - **Spanning:** the junction lies inside one mate's read. The piece holding both mates is the *representative*,
    written as an ordinary pair with the proper-pair flag and a TLEN over its own span. The other piece is one
    supplementary record (`0x800`), hard-clipped with `HardClip`, SEQ trimmed to its aligned part, TLEN 0, its mate
    fields pointing at the other mate (`:43-56`; `ReadAlign_alignBAM.cpp:186-218, 548-577`).
    - `SA` sits only on the split mate's representative record and on the supplementary record, each naming the other
      (`:68-84`).
    - **The representative is chosen by mate structure, not by length.** It is whichever piece holds both mates, so the
      primary record can be the short piece.
  - **Encompassing:** the junction lies between the mates, each mate wholly in one place. Both mates are written as
    ordinary single-mate records with each other's mate fields: no `SA`, no supplementary, no proper-pair flag. Rigel's
    pair rule already rejects these as two references, one strand, outward or too long.
- **Multimapping chimeras** (`--chimMultimapNmax` > 0, the new scheme):
  - every chimeric alignment is written with NH = their number and its own HI;
  - only the best one is primary; the others' records, supplementary ones included, are secondary (`0x900`);
  - the representative's unsplit mate is repeated in each chimeric alignment, identical, which Rigel aliases.
- **Checked on the panel and cfRNA-A.** Every chimeric record has the flags, mate fields, TLEN and `SA` above. Every
  supplementary record is hard-clipped and carries `SA`; no primary or secondary record is hard-clipped.

### 2.1 Why the simulation's chains exist: false mate merges

The recommended parameters switch on mate merging (`--peOverlapNbasesMin 10`):

- `peMergeMates` calls two mates overlapping when R1's end matches the reverse complement of R2's start for 10 or more
  bases with at most 1 % mismatches. Below 100 bases that means an exact match. It judges the sequences alone, whatever
  the true fragment length (`ReadAlign_peOverlapMergeMap.cpp:80-105`, `SequenceFuns.cpp:317-339`).
- The merged read is mapped and converted back to a pair, and chimeric detection runs on it.
- The chimeric result is adopted whenever it is chimeric, whatever the original pair scored (`:72`).

On g50 ss.50 OFF:

- 2,430 of 2,435 chained templates have such an overlap (54 % at 10–14 bases), against a third of the rest;
- 84 % of them come from fragments longer than 200 bp, whose mates cannot overlap at all.

A chance match within a repeat makes STAR merge the mates. It then aligns the merged read in two pieces and puts the
short piece beside the other mate as the representative. The correct linear pair is discarded.

## 3. Today, in Rigel

- **The pieces are glued onto the primary.** `template_hits` attaches a mate's supplementary records to its primary
  location (`bam_scanner.cpp:963-966, 994-1005`), including the pieces of a secondary chain.
- **The pair is judged on the primary**, its 5′ soft clip walked back (`:1035-1047`; `resolve_context.h:840-842`).
- **The fragment is built from every piece**, with every piece's NM and junctions (`bam_scanner.cpp:848, 866-878`).
- **Clips: only the outermost CIGAR operation is read** (`:591-598`).
- **The pair rule demands both mates' aligned bases inside [a, b)**, a and b being the 5′ ends with their soft clips
  walked (`resolve_context.h:840-842`). Mates that overrun each other are rejected as `outward`.

## 4. What chains are, measured

| library | templates | with a chain | pieces on another chromosome |
|---|---:|---:|---:|
| STAR panel, 5 conditions | 10 M each | 0.014–0.024 % | 3 % (g50 ss.50 OFF) |
| `vcap_mix` | 18.6 M | 0.36 % | 63 % |
| cfRNA-B | 1.30 M | 0.9 % | — |
| cfRNA-C | 0.88 M | 6.7 % | — |
| cfRNA-A | 155 k | 8.6 % | 87 % |

- **In the simulation every chain is an artifact.** For gDNA, the supplementary piece is at the read's true position
  and the representative is not in 1,553 of 1,956 chains; the representative is at the true position in 9.
- **In real libraries** chains mix:
  - true chimeras and library chimeras;
  - false merges;
  - adapter tails that STAR aligned chimerically;
  - spliced reads placed partly on a processed pseudogene;
  - unannotated junctions: a local chain in splice order implies an unannotated junction in 460 of 516 cases in
    cfRNA-A and 7,355 of 7,851 in `vcap_mix`.
- **Today half of cfRNA-A's chains are already withheld** as `unpaired` (7,137) and 1,873 are called chimeric; 4,200
  are placed from stitched fragments.
- **Every chimeric fragment today comes from a chain**, except one misfire of the span rule (§8).

## 5. Rescue rules, measured

**Overrun.** Two mates whose 5′ ends face each other, with one mate's aligned bases running past the other's 5′ end.

| library | overrun pairs | |
|---|---:|---|
| STAR panel, g50 ss.50 OFF | 0.03 % | every gDNA one misplaced (1,621 of 1,621; 719 off by more than 100 bp) |
| `vcap_mix` | 0.006 % | |
| cfRNA-B | 0.1 % | |
| cfRNA-A | 2.0 % | distances peak at even values (2, 4, 6, 8) |

- The simulator truncates a read to its fragment (`wgs_engine.py:622`), so a correctly placed simulated pair can never
  overrun. A simulated overrun is a misplacement by construction.
- In real reads the mechanism is the protocol: a read whose 5′ bases were trimmed before alignment starts inside its
  molecule, and its mate reading through shows the true end. The rate follows the protocol (short, trimmed cfRNA
  fragments) rather than STAR's misplacement rate, which suggests most real overruns are such molecules. Real data has
  no truth to confirm it.
- The reads of an overrun pair are contiguous alignments facing each other, so they are attributable to the stretch
  spanning both mates: the owner's criterion. Today's containment test is stricter than the criterion: it assumes a
  read's 5′ end is its molecule's end, which trimming breaks.

**Tail.** A chained read whose foreign piece covers only its 3′ part, running past the other mate's 5′ end.

- In cfRNA-A it would recover 4.9 % of templates, which suggests read-through tails aligned chimerically.
- On the simulation it accepts 700 gDNA templates, every one wrong (570 off by more than 100 bp). These are the false
  merges, whose representative is the misplaced piece.
- The rule trusts STAR's representative pair, and the BAM cannot say whether STAR merged the mates correctly.

## 6. The plan

1. **Item 13: chimeric alignments withheld.**
   - A mate location is chimeric when its record carries `SA` (any of its aliases' records), or its 5′ end is
     hard-clipped.
   - A mate with a supplementary record named by no `SA` makes all of that mate's locations chimeric: the association
     is ambiguous, and the fix contract withholds then.
   - A hit with a chimeric mate is not usable (`HIT_CHIMERIC`, after `missing_mate` and `off_index`). A template left
     with no usable hit is withheld as `chimeric` (after `off_index`, before `unannotated_splice` and `unpaired`).
   - Supplementary records join no hit, so no block, junction or NM of theirs reaches a fragment. The annotated BAM
     writes each one with its mate's location, as an alias: a withheld template's records are all unresolved, `ZC`
     `chimeric`.
   - Counted as `hits_rejected.chimeric` and `withheld.chimeric`.
   - Cost:
     - the panel: 0.014–0.024 % of templates, about 40 % of them withheld today;
     - `vcap_mix`: 0.36 %;
     - cfRNA-A: 8.0 %, of which 5.9 % of the library is unquantified today; 2.6 % of the library is newly withheld,
       placed today from stitched fragments.
2. **The overrun rule.** Mates whose 5′ ends face each other (a < b) are one molecule from the lower of a and the
   reverse mate's first aligned base to the higher of b and the forward mate's last.
   - Its length, along the genome or a witness transcript, must lie in 1..`max_frag_length`.
   - Mates facing away (a ≥ b) stay `outward`.
   - Gain: 2.0 % of cfRNA-A's pairs, 0.84 % of cfRNA-C's, 0.1 % of cfRNA-B's. Cost: 0.03 % of the panel's pairs, each
     a misplaced mate: none of the gDNA ones within 3 bp of its true molecule, 44 % off by more than 100 bp.
   - It needs no constant, and it subsumes the 1–2-base overhang the owner named. Its own commit and A/B.
3. **Deferred: the tail rule.** It waits for a simulation with adapter read-through and 5′ trimming, which would give it
   truth. Until then those templates are withheld as `chimeric`.

## 7. Tests and breaks

**Failing on today's code.**

- Chimeric:
  - a spanning chimera, its representative pair and supplementary record as STAR writes them, is withheld `chimeric`,
    not chimeric-fragment, intergenic or placed;
  - a chimeric mate with a self-contained alternative keeps the alternative;
  - a secondary chimera's pieces no longer join the primary location;
  - a chain without `SA` withholds its mate's locations;
  - a hard-clipped 5′ end is withheld; a hard clip at the 3′ end is not;
  - the annotated BAM writes every record of a withheld template unresolved, `ZC` `chimeric`, supplementary included.
- Overrun:
  - 42-base mates offset by 4 at both ends are one 46-base molecule;
  - mates facing away stay `outward`;
  - an overrun whose span exceeds `max_frag_length` is `length`.

**Breaks, each seen to fire:** the pieces glued on again; `SA` ignored; the 5′ hard clip ignored; the outcome ignoring
`chimeric`; supplementary records left to the sweep; containment restored.

## 8. What item 13 lets us delete (separate steps)

- **`n_deposit_not_offered` and the deposit adapter's multi-reference refusal:** a hit's blocks then always lie on one
  reference.
- **The resolver's chimera classification.**
  - Every chimeric fragment today is a chain, except one misfire: a spliced pair 445 bp along the genome once its
    junction is removed, called chimeric over its 2,162-bp raw span.
  - With chains withheld the classification guards nothing the pair rule and the annotated-path rule do not.
  - Deleting it closes `ISSUES: chimera-span-over-exon-blocks-only`, and frees the word `chimeric` for the withheld
    outcome alone.

## 9. Measurement

- The STAR panel against main, pinned and fractional, per stratum, with the census and every capture efficiency at 1,
  for each step.
- The oracle ladder: no chains and no overruns, so it is identical by construction; to confirm.
- `vcap_mix` and cfRNA-A–C as test inputs: the counts against §4–5.
