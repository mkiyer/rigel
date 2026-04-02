# Annotated BAM Tag Redesign (v4)

## Design Principles

1. **ZC is untouched.** It already captures alignment topology / pipeline
   disposition with full fidelity. Post-processing can collapse categories
   as needed — it is always better to preserve information in the BAM.

2. **The EM has exactly two component types**: transcript and gDNA. This is
   a binary classification — no need for a string tag.

3. **One compact integer bitfield (ZF) replaces ZP and adds transcript
   properties.** Assignment result + transcript flags in a single tag,
   eliminating the redundant ZP string entirely.

4. **Bit 0 is `is_resolved` — the most fundamental flag.** Analogous to
   SAM FLAG bit `0x4` (`is_unmapped`): technically derivable from other
   fields, but so fundamental that every downstream tool benefits from
   having it as a single bit test. All other ZF bits imply `is_resolved`.

5. **Transcript properties belong on the transcript, not the pool.** But
   for user convenience, `is_nrna` and `is_synthetic` are stamped directly
   on each record so users don't need an index join.

---

## Problem with the Current Schema

### Current ZP (5 string values — conflated)

| ZP value | What it encodes |
|----------|-----------------|
| `mRNA` | EM winner is a multi-exon transcript |
| `nRNA` | EM winner is a single-exon transcript |
| `gDNA` | EM winner is the gDNA component |
| `intergenic` | No transcript overlap (pipeline disposition) |
| `chimeric` | Chimeric alignment geometry (pipeline disposition) |

Problems:
- `mRNA` vs `nRNA` is not an EM distinction — the EM treats all transcripts
  identically. This is a post-hoc label derived from `is_nrna` on the winner.
- `intergenic` and `chimeric` are pipeline dispositions, not EM components.
  They duplicate information already in ZC.
- String encoding wastes space for what is fundamentally a 2-bit result
  (transcript / gDNA / unmodeled).

### Current ZC (6 string values — no changes needed)

| ZC value | Meaning |
|----------|---------|
| `unambig` | Single transcript candidate |
| `ambig_same_strand` | Multiple candidates, same strand |
| `ambig_opp_strand` | Multiple candidates, opposite strand |
| `multimapper` | Candidates across multiple loci |
| `chimeric` | Chimeric alignment geometry |
| `intergenic` | No transcript overlap |

ZC is correct as-is. The `ambig_same_strand` / `ambig_opp_strand` split
preserves scoring-relevant information. Collapsing is trivial in
post-processing if desired.

---

## Revised Schema

### ZC — Fragment Class (UNCHANGED)

No modifications. Retains all 6 current values.

### ZP — REMOVED

Replaced entirely by ZF. A string tag for a binary classification was
over-engineering.

### ZF — Assignment Flags (new integer bitfield)

A single `int32` tag encoding resolution status, EM assignment result, AND
transcript properties.

| Bit | Mask | Flag | Meaning |
|-----|------|------|---------|
| 0 | 0x1 | `is_resolved` | Fragment was scored and assigned to an EM component |
| 1 | 0x2 | `is_gdna` | Fragment assigned to the gDNA EM component |
| 2 | 0x4 | `is_nrna` | Assigned transcript is single-exon (nascent RNA candidate) |
| 3 | 0x8 | `is_synthetic` | Assigned transcript is a rigel-generated nRNA span |

#### Bit hierarchy (implication chain)

```
is_synthetic → is_nrna → is_resolved
is_gdna → is_resolved
is_nrna ⊕ is_gdna  (mutually exclusive)
```

All flags above bit 0 **imply** `is_resolved`. This means:
- Every valid non-zero ZF value has bit 0 set (is odd)
- `ZF = 0` unambiguously means "not resolved" — no ZC lookup needed
- `ZF & 1` is the fastest possible "was this fragment quantified?" check

#### Valid ZF values

| ZF | Binary | Interpretation |
|----|--------|----------------|
| 0 | `0000` | Not resolved (intergenic, chimeric, or filtered) |
| 1 | `0001` | Resolved → transcript (multi-exon, annotated) |
| 3 | `0011` | Resolved → gDNA component |
| 5 | `0101` | Resolved → transcript (single-exon nRNA, annotated) |
| 13 | `1101` | Resolved → transcript (single-exon nRNA, synthetic) |

Five valid values. Every resolved fragment is odd-numbered. Every unresolved
fragment is zero.

#### Invalid ZF values (structural constraints)

| ZF | Why invalid |
|----|-------------|
| 2 | `is_gdna` without `is_resolved` |
| 4 | `is_nrna` without `is_resolved` |
| 6 | `is_gdna` + `is_nrna` — mutually exclusive |
| 7 | `is_gdna` + `is_nrna` + `is_resolved` — contradictory |
| 8 | `is_synthetic` without `is_nrna` — impossible |
| 9 | `is_synthetic` + `is_resolved` without `is_nrna` — impossible |
| 10+ | Various contradictions |

Invalid values serve as free error detection — if observed, something is wrong
in the annotation pipeline.

#### Why `is_resolved` at bit 0?

1. **SAM precedent.** SAM FLAG puts the most fundamental classifications at
   the lowest bits (`0x1` = paired, `0x4` = unmapped). `is_resolved` is the
   most fundamental question about a rigel fragment: did it get quantified?

2. **ZF becomes fully self-interpreting.** Without `is_resolved`, `ZF=0`
   was ambiguous (transcript assigned? or unmodeled?). With it, `ZF=0`
   unambiguously means "not resolved." No cross-reference to ZC needed.

3. **The bit test is the most common query.** "How many fragments were
   quantified?" is the first thing every downstream analysis asks. `ZF & 1`
   is O(1), single-instruction. Compare to the alternative: string
   comparison on ZC against two values.

4. **Intentional redundancy.** Yes, `is_resolved` is derivable from ZC.
   So is SAM's `is_unmapped` from CIGAR. The redundancy is the point —
   the most fundamental classification deserves a direct bit test, not a
   derived computation.

5. **All other bits assume it.** `is_gdna`, `is_nrna`, `is_synthetic` are
   only meaningful for resolved fragments. Having `is_resolved` at bit 0
   establishes a clean hierarchy: check bit 0 first, then interpret the rest.

---

## Complete Tag Table

| Tag | Type | Description | Status |
|-----|------|-------------|--------|
| `ZT` | Z (string) | Assigned transcript ID (`.` if none) | unchanged |
| `ZG` | Z (string) | Assigned gene ID (`.` if none) | unchanged |
| `ZI` | i (int32) | Transcript index (`-1` if none) | unchanged |
| `ZJ` | i (int32) | Gene index (`-1` if none) | unchanged |
| `ZC` | Z (string) | Fragment class: `unambig`, `ambig_same_strand`, `ambig_opp_strand`, `multimapper`, `chimeric`, `intergenic` | **unchanged** |
| `ZF` | i (int32) | Assignment flags bitfield (see above) | **new (replaces ZP)** |
| `ZW` | f (float) | Posterior probability of assignment | unchanged |
| `ZH` | i (int32) | Primary-hit flag (1 = winning alignment, 0 = secondary) | unchanged |
| `ZN` | i (int32) | Number of candidate EM components | unchanged |
| `ZS` | Z (string) | Splice type: `spliced_annot`, `spliced_unannot`, `unspliced`, `unknown` | unchanged |
| `ZL` | i (int32) | Locus ID (`-1` if none) | unchanged |

11 tags total (same count as current). ZP removed, ZF added.

---

## Before / After Mapping

| Fragment scenario | Before | After |
|-------------------|--------|-------|
| EM assigns multi-exon transcript | `ZP=mRNA` | `ZF=1` |
| EM assigns single-exon transcript | `ZP=nRNA` | `ZF=5` |
| EM assigns synthetic nRNA span | `ZP=nRNA` | `ZF=13` |
| EM assigns gDNA | `ZP=gDNA` | `ZF=3` |
| No transcript overlap | `ZP=intergenic` | `ZF=0` |
| Chimeric geometry | `ZP=chimeric` | `ZF=0` |

Key improvements:
- ZP eliminated — no more redundancy with ZC
- `ZF=0` unambiguously means "not resolved" (no ZC check needed)
- nRNA/synthetic status encoded (previously lost in generic `nRNA` label)
- All resolved fragments have `ZF & 1` set (odd-numbered)

---

## Common User Queries

### samtools CLI

```bash
# All resolved (modeled) reads
samtools view -d ZF:1 -d ZF:3 -d ZF:5 -d ZF:13 sample.bam
# Or: all reads where ZF != 0 (awk)
samtools view sample.bam | awk '!/ZF:i:0/'

# All transcript-assigned reads (resolved, not gDNA)
samtools view sample.bam | awk '/ZF:i:1\b/ || /ZF:i:5\b/ || /ZF:i:13\b/'

# All gDNA reads
samtools view -d ZF:3 sample.bam

# All nRNA-assigned reads (ZF=5 or ZF=13)
samtools view sample.bam | awk '/ZF:i:5\b/ || /ZF:i:13\b/'

# All intergenic/chimeric (unresolved) reads
samtools view -d ZF:0 sample.bam
# Or by class:
samtools view -d ZC:intergenic sample.bam
```

### Python / pysam

```python
zf = read.get_tag("ZF")
is_resolved  = (zf & 0x1) != 0
is_gdna      = (zf & 0x2) != 0
is_nrna      = (zf & 0x4) != 0
is_synthetic = (zf & 0x8) != 0

# Assigned to a transcript? (resolved and not gDNA)
is_transcript = is_resolved and not is_gdna
```

### Numpy / vectorized (for batch analysis)

```python
zf = np.array([r.get_tag("ZF") for r in reads])
resolved_mask  = (zf & 1).astype(bool)
gdna_mask      = (zf & 2).astype(bool)
nrna_mask      = (zf & 4).astype(bool)
synthetic_mask = (zf & 8).astype(bool)
transcript_mask = resolved_mask & ~gdna_mask
```

---

## Internal Representation (Python)

```python
# ZF flag bits (written to BAM as ZF:i tag)
ZF_RESOLVED: int  = 0x1    # bit 0
ZF_GDNA: int      = 0x2    # bit 1
ZF_NRNA: int      = 0x4    # bit 2
ZF_SYNTHETIC: int  = 0x8   # bit 3

# Pre-computed valid ZF values for convenience
ZF_TRANSCRIPT: int     = ZF_RESOLVED                          # 1
ZF_GDNA_RESOLVED: int  = ZF_RESOLVED | ZF_GDNA               # 3
ZF_NRNA_RESOLVED: int  = ZF_RESOLVED | ZF_NRNA               # 5
ZF_SYNTH_RESOLVED: int = ZF_RESOLVED | ZF_NRNA | ZF_SYNTHETIC  # 13

# Internal pool codes (drive ZF computation, not written to BAM)
POOL_CODE_TRANSCRIPT: int = 0
POOL_CODE_GDNA: int = 1
POOL_CODE_NONE: int = 2   # intergenic, chimeric, unassigned
```

The `AnnotationTable` adds a `tx_flags: np.ndarray` (uint8) column. The C++
`stamp_and_write_hit()` receives the `tx_flags` value and writes it as the
`ZF:i` tag (replacing the old `ZP:Z` string tag).

---

## Impact Assessment

### Code changes required

| File | Change |
|------|--------|
| `annotate.py` | Remove 5 `POOL_*` string constants. Add `ZF_*` bit constants + `POOL_CODE_*` internal codes. Add `tx_flags` to `AnnotationTable`. Remove ZP from docstring, add ZF. |
| `bam_scanner.cpp` | Remove `pool_label()`. Remove ZP string tag from `stamp_and_write_hit()`. Add `tx_flags` param → write `ZF:i` tag. |
| `pipeline.py` | `_populate_em_annotations()`: transcripts → `ZF_TRANSCRIPT` or `ZF_NRNA_RESOLVED` or `ZF_SYNTH_RESOLVED`. gDNA → `ZF_GDNA_RESOLVED`. Unmodeled → `0`. |
| `scan.py` | Det-unambig: compute `tx_flags` from `is_nrna`/`is_synthetic`. Chimeric: `tx_flags=0`. |
| `test_annotate.py` | Remove ZP assertions, add ZF assertions. |
| `test_golden_output.py` | Regenerate golden files (`--update-golden`). |
| `MANUAL.md` | Update tag table: remove ZP, add ZF. |

### Files with NO changes needed

- `em_solver.cpp` — winner semantics unchanged
- `scoring.cpp` — internal scoring unchanged
- `resolve.cpp` — resolution logic unchanged
- `locus.py` — component layout unchanged
- `config.py` — no config impact
- `quant.feather` / output columns — unaffected
- `scripts/benchmarking/` — does not parse BAM tags
- `frag_class_label()` in `bam_scanner.cpp` — ZC unchanged

### Backward compatibility

**Breaking change**: ZP tag removed entirely. Users filtering on `ZP:Z:mRNA`,
`ZP:Z:nRNA`, or `ZP:Z:gDNA` must migrate to ZF bitfield queries. Document
in CHANGELOG.md.

---

## Summary

Two orthogonal tags, zero redundancy:

| Axis | Tag | Encoding | Question answered |
|------|-----|----------|-------------------|
| Alignment topology / disposition | **ZC** | string (6 values, unchanged) | What happened to this fragment? |
| Assignment + transcript properties | **ZF** | int32 bitfield | Was it modeled? What was assigned? What kind? |

ZF bit layout:

| Bit | Flag | Implication |
|-----|------|-------------|
| 0 | `is_resolved` | Most fundamental: was this fragment quantified? |
| 1 | `is_gdna` | Implies `is_resolved` |
| 2 | `is_nrna` | Implies `is_resolved`; mutually exclusive with `is_gdna` |
| 3 | `is_synthetic` | Implies `is_nrna` (and therefore `is_resolved`) |

ZP is eliminated. ZC is untouched. All information is preserved in fewer tags
with richer semantics. ZF is fully self-interpreting — `ZF=0` always means
"not resolved," no cross-reference needed.
