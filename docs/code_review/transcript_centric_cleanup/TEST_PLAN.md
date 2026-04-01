# Transcript-Centric Cleanup Test Plan

**Date**: 2026-03-31
**Status**: Draft for review

---

## Purpose

This plan defines the tests that should land before or alongside the cleanup work so the refactor is driven by behavior instead of by code motion.

---

## Test categories

## 1. Index and entity-construction tests

### 1.1 Same-strand cross-gene span support

Construct a locus where annotated transcripts on the same strand but with different gene labels contribute to one nRNA span.

Assertions:

- the nRNA entity exists
- it does not inherit a single `g_id` or `g_name`
- transcript membership is many-to-many or represented via links
- transcript outputs remain transcript-centric

### 1.2 Annotated equivalent without duplication

Construct a multi-exon transcript whose merged nascent span is fully covered by an annotated single-exon transcript.

Assertions:

- no duplicate synthetic transcript row is created
- the annotated single-exon transcript remains the only canonical row for that span
- any `nrna_t_index` or equivalent link resolves to that same row
- the EM never sees separate mature and nascent components for the same interval

### 1.3 Gene-table hygiene

Assertions:

- derived gene transcript counts exclude nRNA entities
- transcript counts reflect annotated transcripts only

---

## 2. Scoring and routing tests

### 2.1 Annotated-equivalent nRNA candidate participates in routing

Use an unspliced fragment compatible with:

- an annotated single-exon transcript
- the corresponding nRNA entity
- optional gDNA

Assertions:

- the annotated single-exon transcript is retained as a valid candidate
- no duplicate mature/nascent candidate pair is produced for the same single-exon interval
- routing metadata remains sufficient to distinguish annotated versus implied transcript rows

### 2.2 Cross-boundary nRNA entity does not require gene ownership

Assertions:

- candidate construction succeeds
- no code path assumes one owning gene for the nRNA entity

---

## 3. Locus and EM tests

### 3.1 No duplicate single-exon EM states

Build a small synthetic locus with:

- two annotated transcripts
- one annotated single-exon nascent-equivalent transcript
- one synthetic implied nascent transcript representing a genuinely distinct span

Assertions:

- the annotated single-exon span appears exactly once in the component layout
- the distinct implied nascent span appears once as its own transcript row
- no non-identifiable duplicate pair is introduced

### 3.2 Reporting views do not imply duplicate EM states

Assertions:

- transcript outputs are derived from the canonical transcript table
- any nascent reporting view is derived from the same underlying rows rather than from duplicate EM components

---

## 4. Output tests

### 4.1 Transcript output contains annotated transcripts only

Assertions:

- no synthetic nRNA rows appear in transcript output
- transcript `mrna` counts remain correct
- annotated single-exon nascent-equivalent transcripts remain present in transcript output

### 4.2 Nascent reporting does not duplicate model state

Assertions:

- any retained nascent report is a view over the same transcript rows
- annotated single-exon transcripts are not materialized as separate nascent components just for reporting
- synthetic implied nascent transcripts remain reportable as distinct transcript rows because they represent distinct modeled spans

### 4.3 Gene output is convenience-only and mature-only

Assertions:

- gene transcript counts exclude nRNA entities
- gene mRNA aggregates only annotated transcript mRNA
- no code path forces nRNA into one gene owner

---

## 5. API contract tests

### 5.1 `quant_from_buffer()` requires calibration

Assertions:

- calling without calibration fails immediately
- calling with calibration succeeds

### 5.2 Error messaging is explicit

Assertions:

- the raised error says that calibration is mandatory

---

## 6. Scheduler tests

### 6.1 Small-workload batching

Construct a case where total work is smaller than thread count.

Assertions:

- not every locus is classified as mega
- multiple loci are still batched together

### 6.2 Mixed workload balance

Construct a case with one large locus and many small loci.

Assertions:

- the large locus may run alone
- the tail of small loci is still balanced across batches

### 6.3 Stable scheduling under permutation

Assertions:

- shuffling locus order does not change batch quality materially
- any deterministic scheduler remains deterministic under fixed input

---

## 7. Golden-output refresh rules

Golden files should only be updated after the new ownership model is fully in place.

Expected intentional changes:

1. transcript outputs no longer contain synthetic nRNA rows
2. any retained nascent report becomes a view over the canonical transcript table rather than a second canonical state space
3. gene outputs stop inheriting transcript counts from synthetic implied nascent rows

Golden refresh should happen after the output contract is stable, not during the intermediate migration.

---

## Suggested landing order for tests

1. API contract tests for mandatory calibration
2. index/entity-construction tests
3. routing tests for annotated-equivalent entities
4. locus/EM layout tests
5. scheduler tests
6. golden-output refresh

---

## Acceptance checklist

The cleanup is ready to merge when all of these are true:

1. The new targeted tests pass.
2. Existing unit and scenario tests pass or are intentionally updated.
3. Golden outputs are regenerated only after the new output contract stabilizes.
4. No remaining test relies on synthetic nRNA rows living inside the transcript table.
5. No remaining test assumes an annotated single-exon transcript is duplicated into separate mature and nascent EM states.
