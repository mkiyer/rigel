# Accounting: one outcome per template, one status per record: design (roadmap item 15)

2026-10-21 · code base `fef6b0ce` · for the owner's review. Every claim about Rigel's code cites a line re-read for this
document. Every number is measured; the receipts are in `$RIGEL_RUNS/arms/2026-10-21_accounting/`. Approved 2026-10-21
(owner: the `read_pairs` block, schema 5; ZC keeps `unpaired`, its reasons counted in `summary.json`). The steps land
one at a time, their measurements the roadmap's item 15 (`PAIRED_END_ROADMAP.md`).

## 1. The recommendation in brief

The scan already counts most of what it decides, and its counters reconcile: read groups = fragments + withheld on
every BAM measured, and genic fragments = the units routed on the STAR panel. Four things are missing or wrong:

1. **Read pairs the scan never sees.** A read name all of whose records are filtered (both mates unmapped, QC-failed,
   or skipped duplicates) never becomes a read group, and nothing counts it: 2.8–22 % of the cfRNA libraries' pairs,
   0.7 % of `vcap_mix`'s. The report's "Unmapped" read-fate segment is therefore always 0.
2. **Records have no status.** In a placed template, the records of a hit the pair rule rejected carry the winning
   fragment's tags, indistinguishable from a losing placement's: 7–46 k records in the cfRNA libraries, 0.2–0.7 M per
   10 M panel templates. A record in no candidate hit is stamped intergenic gDNA (8–1.5 k per real library), and a
   supplementary record is written as part of the hit beside it.
3. **Counters that mean something else.** `proper_pairs` sits among the record counts but counts pairs, exactly half
   of the records flagged on every BAM. `unique_reads` / `multimapping_reads` read the aligner's NH, not the usable
   alignments that are the template's multiplicity since 2026-10-19: 21–37 k templates per 10 M on the panel are
   multimapping by NH with one usable alignment, and `--no-include-multimap` drops them. A multimapping fragment's
   hits all carry one splice type (ZS), wrong for 0.4–5.6 % of its hits. The annotation pass's own counters count
   hits, and file the multimap option's drops as intergenic.
4. **Nothing compares the two passes.** The annotation pass re-judges every template with the scan's code, held to
   the scan only by a frag_id counter kept in step.

The design gives every template exactly one outcome and every record exactly one status, decided and counted by the
scan, stamped by the annotated BAM, and the annotation pass proves template by template that it saw what the scan saw,
or the run stops. No quantification number moves at the defaults; only `--no-include-multimap` changes, in its own
step.

## 2. What is ruled

- **The wave-1 plan, §10** (revision 5, answering the targeted review's finding 4; item 15 of the fix order of
  2026-10-18):
  - each template has exactly one outcome, decided in a fixed order;
  - each record's status follows from its template's final outcome: **used** (in a usable hit of a placed template),
    **unused** (in a placed template, in no usable hit), **unplaced** (every record of a template with any other
    outcome);
  - the multimap option excludes a template with more than one usable hit, counted after the annotated-path rule;
  - the counts reconcile by identity as well as by total; tags carry the outcome and the record's status
    (`ZC:unused_alignment`), with stale tags cleared;
  - the annotation pass's groups must equal the scan's by identity, not only by count, or the run stops.
- **What has landed since** shortens the plan's outcome list: HI is not read (2026-10-19), so its HI cases are gone; a
  chimeric alignment is withheld before resolution and the resolver's chimera call is deleted (2026-10-20), so
  "reported chimera" is gone and `chimeric` is a withheld reason; so is `off_index` (2026-10-20); and a template's
  multiplicity is its usable hits (2026-10-19).

## 3. Today

### 3.1 The scan

- **Records.** The reader counts each record's flags and drops the filtered ones before it reads the name: QC-failed,
  then unmapped, then a duplicate when duplicates are skipped (`bam_scanner.cpp:1584-1589`). The precedence makes them
  a partition, so a QC-failed unmapped record is counted QC-failed only and `mapped_reads = total − unmapped`
  (`cli.py:565`) calls it mapped. Secondary and supplementary flags are counted on kept records only (`:1598-1599`).
- **Templates.** A group is a run of kept records sharing a name, given the next frag_id (`:1601-1609`). A name with
  no kept record forms no group.
- **The outcome** of a group (`process_qname_group_threaded`):
  1. counted `unique` or `multimapping` by the aligner's NH or a secondary flag; with `--no-include-multimap` a
     multimapping group returns here, before any hit is judged, counted only as `multimapping` (`:1867-1873`);
  2. `mate_unmapped`, `proper_pair` or `improper_pair`, once from its first primary records (`:1881-1896`);
  3. every hit judged, each rejected hit counted by reason; a group with no usable hit is withheld, `off_index` >
     `chimeric` > `unannotated_splice` > `unpaired` (`judge_hits`, `:1101-1119`; `:1903-1910`);
  4. otherwise one fragment (`:1932`): intergenic when no usable hit resolves (`:2069`).
- **Two holes.** A group with no candidate hit, every kept record supplementary, is not withheld: `judge_hits` returns
  "not withheld" for an empty hit list (`:1114`), so it is counted a fragment and yields nothing. And the pair rule
  passes a mate with no aligned base unjudged (`:1065`), the scan skips a usable hit with no aligned block uncounted
  (`:1937`), which leaves `collect_placement`'s own no-block branch unreachable (`:1766`).

### 3.2 What is reported

- `summary.json` reports the records, the read groups (`read_groups`, `unique_reads`, `multimapping_reads`), the
  fragments, the withheld groups by reason and the rejected hits by reason (`cli.py:562-624`). `proper_pairs` sits
  among the record counts but is step 2's per-group count (`:569`); `mate_unmapped` and `improper_pair` are never
  reported.
- The report's read fate splits the read groups into unique, multimapping and "unmapped" `= read_groups − unique −
  multi` (`report/model.py:176`), 0 by construction; its "proper pairs" row divides a pair count by the record count
  (`:203`).

### 3.3 The annotation pass

`BamAnnotationWriter::write` (`bam_scanner.cpp:2433-2794`) re-reads the BAM, writes each filtered record through
untagged, groups the rest by name, and re-judges each group with the scan's `template_hits` and `judge_hits`
(8c6135f5). Its frag_id is a counter advanced once per group (`:2537`, `:2560`, `:2711`).
- **Withheld**: ZF 0, ZC its reason, ZH 0 (`:2541-2562`).
- **Dropped by the multimap option**: ZF 0x80, ZC `.`, and ZH 1 on every record; counted `intergenic` (`:2521-2539`).
- **Placed**: each hit's records, written once with the first hit holding them, get the fragment's annotation row,
  the winning placement ZH 1 and every other hit ZH 0, **a rejected hit's included** (`:2580-2664`). With no row, every
  hit is stamped intergenic gDNA, 0x25 (`:2665-2681`): the intergenic fast path, and also **a gated fragment**, which
  the scorer drops without a row (`scoring.cpp:738`, `:894`), at an alpha of 0 only. Records in no hit are swept up
  and stamped 0x25 (`:2687-2709`).
- **ZS** is the fragment's splice type; a multimapping fragment's is `spliced_annot` when any hit is, else `unspliced`
  (`scoring.cpp:507`), stamped on every hit.
- **The grouping differs from the scan's in one case.** The scan's group key moves with kept records only, so a
  filtered record of another name between two records of one name leaves the name one group; the pass closes the held
  group at a filtered record of another name (1cde0664) and splits it, and every later frag_id is one off. It takes an
  uncollated BAM; none of §4's inputs is one.
- **Docs.** `annotate.py`'s ZC list omits `off_index` and `chimeric` (`:66-73`); its bit-5 text and invariant 2 still
  allow a spliced intergenic fragment, 0x21 (`:120-122`, `:133-135`), as does the MANUAL's tally example
  (`MANUAL.md:851-852`), which the annotated-path rule has made impossible; the MANUAL defines ZF 0 as "scored, but
  assigned no pool" (`:656`), while a withheld pair is never scored.

## 4. Measured

**The instrument** (`audit_pass.py`; `audit/instrument.patch`, measurement only): the scan, then the annotation pass
on the same resolver, re-judging every template with the scan's code and counting each template's outcome and each
record's status. On all nine BAMs:
- every identity of §5.4 holds (72 checks), and the instrument's placed and withheld counts equal the scan's own;
- no read name reappears after another, so today's two passes group alike on these inputs;
- no template lacks a candidate hit, no usable hit lacks an aligned block, and no mapped record lacks an aligned
  base: §3.1's holes are unreachable here.

**Templates** (read names):

| BAM | read pairs | excluded | withheld: unpaired / unannotated splice / off index / chimeric | placed: unique / multimapping | NH multimapping, one usable hit |
|---|---:|---:|---:|---:|---:|
| g00 κ0.99 OFF | 10.00 M | 1.7 k (0.02 %) | 12 / 1.6 k / 0 / 1.5 k | 9.45 M / 544 k | 35 k |
| g05 κ0.50 OFF | 10.00 M | 102 k (1.02 %) | 12 / 1.6 k / 0 / 1.5 k | 9.34 M / 559 k | 33 k |
| g50 κ0.50 OFF | 10.00 M | 998 k (9.98 %) | 183 / 1.8 k / 0 / 2.4 k | 8.31 M / 687 k | 21 k |
| g50 κ0.99 OFF | 10.00 M | 997 k (9.97 %) | 174 / 2.0 k / 0 / 2.4 k | 8.31 M / 686 k | 21 k |
| g50 κ0.99 ON | 10.00 M | 36 k (0.36 %) | 22 / 1.9 k / 0 / 1.4 k | 9.33 M / 635 k | 37 k |
| cfRNA-A | 168 k | 12 k (7.41 %) | 1.4 k / 2.2 k / 0 / 13 k | 129 k / 9.9 k | 275 |
| cfRNA-B | 1.34 M | 38 k (2.84 %) | 1.1 k / 1.7 k / 0 / 12 k | 1.23 M / 56 k | 546 |
| cfRNA-C | 1.12 M | 248 k (22.06 %) | 5.8 k / 5.5 k / 0 / 57 k | 742 k / 66 k | 2.3 k |
| `vcap_mix` | 18.70 M | 135 k (0.72 %) | 32 k / 115 k / 0 / 67 k | 18.23 M / 128 k | 6.3 k |

**Records** (`unused` split by why; `proper_pairs`: the records flagged, against today's count):

| BAM | records | filtered | used | unused: rejected hit (supplementary) / no candidate hit | unplaced | `proper_pairs`: flagged / reported | multimapping hits stamped another hit's ZS (implicit) |
|---|---:|---:|---:|---:|---:|---:|---:|
| g00 κ0.99 OFF | 22.20 M | 3.5 k | 21.99 M | 189 k (4) / 0 | 9.1 k | 20.00 M / 10.00 M | 79 k of 1.52 M (48 k) |
| g05 κ0.50 OFF | 22.40 M | 204 k | 21.95 M | 234 k (45) / 0 | 9.2 k | 19.80 M / 9.90 M | 75 k of 1.60 M (45 k) |
| g50 κ0.50 OFF | 24.25 M | 2.00 M | 21.58 M | 658 k (15) / 0 | 23 k | 18.00 M / 9.00 M | 40 k of 2.26 M (24 k) |
| g50 κ0.99 OFF | 24.23 M | 1.99 M | 21.57 M | 649 k (18) / 0 | 22 k | 18.01 M / 9.00 M | 40 k of 2.25 M (24 k) |
| g50 κ0.99 ON | 22.87 M | 72 k | 22.58 M | 209 k (7) / 0 | 11 k | 19.93 M / 9.96 M | 61 k of 1.85 M (38 k) |
| cfRNA-A | 574 k | 25 k | 343 k | 6.8 k (2.1 k) / 638 | 199 k | 309 k / 154 k | 419 of 60 k (250) |
| cfRNA-B | 3.17 M | 76 k | 2.93 M | 8.3 k (341) / 254 | 156 k | 2.60 M / 1.30 M | 15 k of 270 k (14 k) |
| cfRNA-C | 3.50 M | 496 k | 2.07 M | 46 k (13 k) / 1.5 k | 892 k | 1.74 M / 872 k | 1.4 k of 328 k (365) |
| `vcap_mix` | 38.33 M | 271 k | 37.32 M | 17 k (226) / 8 | 722 k | 37.11 M / 18.56 M | 18 k of 427 k (2.6 k) |

- **The excluded pairs**, on the panel, are the simulator's gDNA drawn from the references' N gaps (996,761 of g50
  κ0.50's 997,628 unmapped pairs are gDNA, their source sequence all N); in the real libraries they are the pairs STAR
  did not align.
- **`unpaired` by its furthest check** (missing mate / reference / strand / outward / length): on the panel almost all
  `length`; cfRNA-C 0 / 2.8 k / 79 / 522 / 2.4 k; `vcap_mix` 317 / 7.9 k / 358 / 4.1 k / 19 k. No missing-mate pair
  has a mate its flags call mapped.
- **A rejected hit's records** are mostly secondary alignments beside the one Rigel kept; in cfRNA-A and cfRNA-C about
  three in ten are supplementary records. **Records in no candidate hit** occur only in the real libraries.

## 5. The design

### 5.1 One template, both passes

A template is a run of records sharing a read name, **filtered records included**: the scan's reader closes a group at
every change of name, as the annotation pass already does. The two passes then group alike by construction, and a
name whose every record is filtered is a template too: **excluded**.

### 5.2 One outcome per template

| outcome, in this order | when | its records |
|---|---|---|
| **excluded** | every record unmapped, QC-failed or a skipped duplicate | filtered |
| **withheld**, by reason | no usable hit: `off_index` > `chimeric` > `unannotated_splice` > `unpaired`; `unpaired` by the furthest check any hit reached, `missing_mate` < `reference` < `strand` < `outward` < `length`. No candidate hit at all is `unpaired`, `missing_mate` | unplaced |
| **excluded by the multimap option** | `--no-include-multimap` and more than one usable hit (step 4; until then today's NH test, counted) | unplaced |
| **placed**, `unique` or `multimapping` | one usable hit, or more | used or unused |

The aligner's NH stays read only for `nh_above_reported`.

### 5.3 One status per record

| status | record |
|---|---|
| **filtered** | unmapped, QC-failed, or a skipped duplicate |
| **used** | its alignment is a location of a usable hit (a secondary record repeating that alignment included) |
| **unused** | in a placed template, in no usable hit: a rejected hit's, one in no hit, a supplementary record |
| **unplaced** | every eligible record of a withheld or option-excluded template |

### 5.4 The counts and their identities

> read pairs = excluded + withheld + excluded by the option + placed unique + placed multimapping
>
> records = filtered + used + unused + unplaced

Both hold by construction, so the scan checks them and a failure stops the run as a bug. Each flag counter counts
every record carrying its flag, no longer a partition by precedence, so `mapped_reads = total − unmapped` holds;
`proper_pairs` counts the records flagged, as its family says. `mate_unmapped`, `improper_pair` and the NH-based
`unique` / `multimapping` go. Writing step 1's gates found one more: `nh_above_reported` counted a template none of
whose primary or secondary records is eligible (NH 1 against no record reported); it no longer does.

**`summary.json`, schema 5.** A `read_pairs` block, one outcome per read name and summing to its `total`: `excluded`,
`withheld` (the four reasons) with `unpaired_by_reason` (its five), `excluded_by_option`, `placed` (`unique`,
`multimapping`).
`alignment_stats` gains `records`: `filtered`, `used`, `unused`, `unplaced`, summing to `total_reads`. `read_groups`,
`unique_reads`, `multimapping_reads` and `fragment_stats.withheld` move into `read_pairs`. The report's read fate
becomes placed unique / placed multimapping / withheld / excluded over every read pair, so unmapped pairs are seen.

### 5.5 The annotated BAM stamps what the scan decided

| template, record | ZF | ZC | ZH | ZS |
|---|---|---|---|---|
| placed, used | the fragment's | its class | 1 on the winning placement, 0 on the others | **the hit's own splice type** |
| placed, unused | the fragment's | **`unused_alignment`** | 0 | `unknown` |
| withheld | 0 | its reason | 0 | `unknown` |
| option-excluded | 0x80 | `.` | **0** (today 1) | `unknown` |
| placed, gated at an alpha of 0 | **0**, from a row the scorer now writes | its class | 0 | the fragment's |
| filtered | written through untagged | | | |

ZT, ZG, ZR, ZI, ZJ, ZW, ZN and ZL keep their meaning, the fragment's on every record of a placed template. A placed
template with no annotation row is then exactly the intergenic fast path and keeps 0x25. The splice type of a used
record comes from one function that resolves a hit, shared by the scan and the pass. §3.3's docs residue is corrected
with it.

### 5.6 The two passes compared, template by template

The scan folds each template into a 64-bit digest: its frag_id, a hash of its name, its record count, its outcome and
its used and unused counts, combined with the existing `hash_combine` and summed over templates, so the digest does not
depend on which worker took which template. The annotation pass computes the same digest as it writes. At the end,
`write_annotated_bam` compares the digests and the outcome counts; a mismatch deletes the annotated BAM and stops the
run, naming both passes' counts. The tables are written before the pass, so no quantification is lost.

### 5.7 What goes

- `is_multimapper`, the NH test, once the option reads the usable hits (step 4).
- The pass's own counters (`n_annotated`, `n_intergenic`, …): it reports the scan's outcomes and statuses instead.
- A mate location with no aligned base becomes no location, as a supplementary record is: with it go the pair rule's
  pass-through (`bam_scanner.cpp:1065`), the scan's silent skip (`:1937`), `collect_placement`'s unreachable branch
  (`:1766`) and the pass's `exons.empty()` guards; a numeric no-op on every input of §4.
- The hygiene ledger's two entries 8c6135f5 fixed: the writer re-implementing `parse_bam_record`, and multimap
  detection written twice.

## 6. The plan

Each step its own commit, gated as every scan change is. No quantification number moves at the defaults, proven by
rename_identity on oracle g00, STAR g50 κ0.50 OFF, cfRNA-A and cfRNA-C, and each step's counts equal this audit's.

1. **Outcomes and statuses, counted** (§5.1–5.4): the grouping at every name, excluded pairs, the no-hit template
   withheld, `unpaired` by reason, the record statuses, the two identities, the flag counters, `summary.json` schema 5,
   the report's read fate, the MANUAL.
2. **The annotated BAM** (§5.5): unused records, the option's ZH, the gated rows, ZS per hit, the docs.
3. **The identity check** (§5.6).
4. **The multimap option on usable hits**: a behaviour change under `--no-include-multimap` only, so its own A/B on
   the STAR panel with the option on both arms, per stratum.
5. **The no-aligned-base location and the ledger** (§5.7).

## 7. Tests and breaks

- **A counters BAM**, hand-built, holding every outcome (each withheld reason and each `unpaired` reason, a name with
  every record filtered, a no-hit template, the multimap option's case, a placed unique and a placed multimapping
  template) and every status (a rejected hit's records, a location in no hit, a supplementary record, a secondary
  record repeating a used alignment): every counter asserted exactly, at one worker and at four, and with each name's
  records permuted.
- **The identities** on that BAM and on a simulated scenario.
- **The annotated BAM** on the counters BAM: each record's tags as in §5.5's table, every record written once, the
  digests equal. A BAM with a filtered record of another name inside a name is split alike by both passes (the digests
  equal); a pass run with another `max_frag_length` than the scan's, or given the BAM with two read names swapped,
  stops the run.
- **Breaks, each seen to fire:** a record's status tied to its hit rather than its template's outcome (the
  option-excluded template's records come out uncounted); the no-hit template counted a fragment; the scan grouping
  kept records only (the split name's digests differ, and the run stops); an unused record stamped 0x25; the option's
  ZH back at 1; ZS from the fragment; the digest without the name (the swapped names pass).

## 8. Measurement

- **Steps 1–3 and 5:** rename_identity bit-identical on the four inputs; the suite and preflight; every new counter
  equal to this audit's count on the nine BAMs of §4.
- **Step 4:** the STAR panel with `--no-include-multimap` on both arms, pinned and fractional, per stratum, with the
  census.

## 9. Questions for the owner

Both answered 2026-10-21, as recommended.

1. The `summary.json` layout of §5.4, a `read_pairs` block in place of the read-group counts (schema 5). Agree?
2. The plan stamps an `unpaired` template's reason into its records. Today ZC reads `unpaired`; stamping the reason
   means five ZC values (`missing_mate` … `length`) in its place, a change to a documented value. Recommended: keep ZC
   `unpaired` and count the reasons in `summary.json`. Agree, or stamp them?
