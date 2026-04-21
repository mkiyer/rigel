# ZF Tag Redesign — Unified Outcome Bitfield

**Status:** plan (not yet implemented)
**Type:** breaking change to annotated BAM schema (pre-release, no backward compatibility required)
**Scope:** `ZF` tag semantics + `ZC` tag cleanup

---

## 1. Motivation

The current annotated BAM has two partly-overlapping tags:

- `ZF:i` — 4-bit bitfield answering "what did the EM decide?" (only meaningful for EM-resolved records)
- `ZC:Z` — string enum mixing two orthogonal axes: input ambiguity (`unambig`, `ambig_same_strand`, `ambig_opp_strand`, `multimapper`) **and** terminal outcomes (`intergenic`, `chimeric`)

Consequences of the current design:

1. **"What happened to this fragment?" requires both tags.** A user who wants to separate intergenic (deterministic gDNA) from multimapper-dropped (unresolved) from chimeric (unresolved) has to inspect `ZC` because `ZF=0` collapses all three.
2. **Total gDNA is not a single mask.** `ZF & 0x2` counts only EM-gDNA, not intergenic fragments that were deterministically bucketed as gDNA.
3. **`ZC` does double duty.** "multimapper" means both "fragment had multiple alignments and went to EM" *and* "fragment was dropped by `--no-multimap`". These are different outcomes.
4. **Latent bug.** The `--no-multimap` skip path in `BamAnnotationWriter` currently stamps `ZC="intergenic"` on dropped multimappers. That's wrong regardless of schema.

Goals of this redesign:

- **ZF alone fully answers "what happened to this fragment?"** — pool assignment + reason if unresolved.
- **Total gDNA = one AND-mask on ZF.**
- **ZC becomes a pure input-ambiguity enum**, orthogonal to the outcome axis.
- **Every legitimate outcome is representable**, including currently-invisible cases (intergenic, chimeric, multimapper-drop).

---

## 2. Design

### 2.1 ZF bit layout

```
bit 0 (0x01)  is_resolved              rigel assigned this fragment to a pool
bit 1 (0x02)  is_mrna                  pool = mRNA
bit 2 (0x04)  is_gdna                  pool = gDNA (intergenic + EM-gDNA)
bit 3 (0x08)  is_nrna                  pool = nRNA (annotated or synthetic)
bit 4 (0x10)  is_synthetic             nRNA subtype: rigel-generated span
bit 5 (0x20)  is_intergenic            gDNA subtype: deterministic fast-path
bit 6 (0x40)  is_chimeric              unresolved reason: chimeric
bit 7 (0x80)  is_multimapper_dropped   unresolved reason: --no-multimap drop
```

### 2.2 Invariants

1. `is_resolved` XOR `(is_chimeric | is_multimapper_dropped)` — every stamped record is exactly one or the other.
2. When `is_resolved`, exactly one of `{is_mrna, is_gdna, is_nrna}` is set.
3. `is_synthetic` → `is_nrna` (synthetic spans are always nRNA).
4. `is_intergenic` → `is_resolved & is_gdna`.
5. `is_chimeric` and `is_multimapper_dropped` are mutually exclusive.
6. Pool subtype bits (`is_synthetic`, `is_intergenic`) are only set when the corresponding pool bit is set.

A unit test will assert invariants 1–6 over all produced ZF values.

### 2.3 Canonical constants

Defined in `src/rigel/annotate.py`:

| Name                     | Value (hex) | Bits                           | Meaning                           |
|--------------------------|-------------|--------------------------------|-----------------------------------|
| `AF_UNRESOLVED`          | `0x00`      | —                              | No assignment, no reason set      |
| `AF_MRNA`                | `0x03`      | resolved \| mrna               | mRNA (multi-exon transcript)      |
| `AF_NRNA`                | `0x09`      | resolved \| nrna               | Annotated single-exon nRNA        |
| `AF_NRNA_SYNTH`          | `0x19`      | resolved \| nrna \| synthetic  | Rigel-generated nRNA span         |
| `AF_GDNA_EM`             | `0x05`      | resolved \| gdna               | EM-resolved gDNA                  |
| `AF_GDNA_INTERGENIC`     | `0x25`      | resolved \| gdna \| intergenic | Deterministic intergenic gDNA     |
| `AF_CHIMERIC`            | `0x40`      | chimeric                       | Chimeric fragment, dropped        |
| `AF_MULTIMAPPER_DROP`    | `0x80`      | multimapper_dropped            | Multimapper dropped by config     |

Primitive bit constants (`AF_RESOLVED`, `AF_MRNA_BIT`, `AF_GDNA_BIT`, `AF_NRNA_BIT`, `AF_SYNTHETIC_BIT`, `AF_INTERGENIC_BIT`, `AF_CHIMERIC_BIT`, `AF_MULTIMAPPER_DROP_BIT`) are exposed for masking.

### 2.4 Derived queries (what users run downstream)

```awk
# Total mRNA
(ZF & 0x02) != 0
# Total gDNA (includes deterministic intergenic)
(ZF & 0x04) != 0
# Total nRNA
(ZF & 0x08) != 0
# Anything resolved (assigned to any pool)
(ZF & 0x01) != 0
# Dropped/unresolved
(ZF & 0xC0) != 0
# Reason for drop
(ZF & 0x40) != 0   # chimeric
(ZF & 0x80) != 0   # multimapper dropped
# Intergenic-only gDNA
(ZF & 0x20) != 0
# EM-resolved gDNA only
(ZF & 0x04) != 0 && (ZF & 0x20) == 0
# Synthetic nRNA only
(ZF & 0x10) != 0
```

### 2.5 ZC redesign

`ZC:Z` becomes a pure input-ambiguity enum reflecting the candidate-set structure **before** EM collapsed it:

| ZC                  | Meaning                                       |
|---------------------|-----------------------------------------------|
| `unambig`           | Single EM candidate                           |
| `ambig_same_strand` | >1 candidates, all same strand                |
| `ambig_opp_strand`  | >1 candidates spanning both strands           |
| `multimapper`       | Multiple genomic alignments                   |
| `.`                 | Fragment never entered EM (intergenic, chimeric, multimapper-dropped) |

Note: "multimapper" as a ZC value now means "was a multimapper **and** entered EM." Dropped multimappers get `ZC="."`. "intergenic" and "chimeric" are removed from ZC — they're ZF-encoded.

The C++ `FRAG_INTERGENIC` / `FRAG_CHIMERIC` enum values (if any) remain internal; the string mapping in `frag_class_label` drops `intergenic` from its output domain for the annotated-BAM stamp path (the internal frag_class code is separate from the ZC string).

### 2.6 ZF / ZC matrix (canonical expected values per outcome)

| Outcome                           | ZF    | ZC                   | ZW          | ZT/ZG/ZR | ZI/ZJ | ZN    |
|-----------------------------------|-------|----------------------|-------------|----------|-------|-------|
| Unambig mRNA                      | `0x03` | `unambig`           | 1.0         | filled   | ≥0    | 1     |
| Unambig nRNA (annotated)          | `0x09` | `unambig`           | 1.0         | filled   | ≥0    | 1     |
| Unambig nRNA (synthetic)          | `0x19` | `unambig`           | 1.0         | filled   | ≥0    | 1     |
| EM winner → mRNA                  | `0x03` | `ambig_*`/`multimapper` | post.    | filled   | ≥0    | ≥2    |
| EM winner → nRNA (annotated)      | `0x09` | `ambig_*`/`multimapper` | post.    | filled   | ≥0    | ≥2    |
| EM winner → nRNA (synthetic)      | `0x19` | `ambig_*`/`multimapper` | post.    | filled   | ≥0    | ≥2    |
| EM winner → gDNA                  | `0x05` | `ambig_*`/`multimapper` | post.    | `.`      | -1    | ≥2    |
| Intergenic (deterministic gDNA)   | `0x25` | `.`                 | 1.0         | `.`      | -1    | 0     |
| Chimeric (dropped)                | `0x40` | `.`                 | 0.0         | `.`      | -1    | 0     |
| Multimapper dropped (`--no-multimap`) | `0x80` | `.`             | 0.0         | `.`      | -1    | 0     |

---

## 3. Breaking changes

This redesign is **not backward compatible**. Consumers of the annotated BAM will need to update their logic. Summary of what changes from the current schema:

- `ZF=1` previously meant "multi-exon transcript assigned" (i.e. mRNA or annotated nRNA). Now it's an invalid value (it would mean "resolved but no pool"). mRNA is `ZF=0x03`.
- `ZF=0` previously conflated intergenic, chimeric, multimapper-dropped, and not-processed. Now intergenic gets `0x25`, chimeric gets `0x40`, multimapper-drop gets `0x80`. `ZF=0` becomes genuinely reserved.
- `is_gdna` (bit 1 in old scheme, value `0x2`) is now bit 2 (value `0x04`) and covers total gDNA, not just EM-gDNA.
- `ZC="intergenic"` and `ZC="chimeric"` no longer appear; those cases use `ZC="."` with ZF bits.

Acceptable because this is pre-release (no external consumers).

---

## 4. Implementation plan

### 4.1 `src/rigel/annotate.py`

- Replace current `AF_*` constants block with:
  ```python
  # Primitive bits
  AF_RESOLVED: int = 0x01
  AF_MRNA_BIT: int = 0x02
  AF_GDNA_BIT: int = 0x04
  AF_NRNA_BIT: int = 0x08
  AF_SYNTHETIC_BIT: int = 0x10
  AF_INTERGENIC_BIT: int = 0x20
  AF_CHIMERIC_BIT: int = 0x40
  AF_MULTIMAPPER_DROP_BIT: int = 0x80

  # Canonical composed values
  AF_UNRESOLVED: int = 0x00
  AF_MRNA: int = AF_RESOLVED | AF_MRNA_BIT              # 0x03
  AF_NRNA: int = AF_RESOLVED | AF_NRNA_BIT              # 0x09
  AF_NRNA_SYNTH: int = AF_NRNA | AF_SYNTHETIC_BIT       # 0x19
  AF_GDNA_EM: int = AF_RESOLVED | AF_GDNA_BIT           # 0x05
  AF_GDNA_INTERGENIC: int = AF_GDNA_EM | AF_INTERGENIC_BIT  # 0x25
  AF_CHIMERIC: int = AF_CHIMERIC_BIT                    # 0x40
  AF_MULTIMAPPER_DROP: int = AF_MULTIMAPPER_DROP_BIT    # 0x80
  ```
- Rewrite `winner_flag(is_nrna, is_synthetic) -> int`:
  - `is_synthetic` (implies nRNA) → `AF_NRNA_SYNTH`
  - `is_nrna` → `AF_NRNA`
  - otherwise → `AF_MRNA`
- Rewrite `winner_flags(is_nrna, is_synthetic)` vectorized to match.
- Update the BAM tag schema docstring (§BAM Tag Schema) and the `frag_class` list (remove `intergenic`, `chimeric`; add `.` to the enum domain).

### 4.2 `src/rigel/scan.py`

- Chimeric annotation site (`tx_flags=AF_UNRESOLVED`) → `tx_flags=AF_CHIMERIC`.
- Update import to include `AF_CHIMERIC`.

### 4.3 `src/rigel/pipeline.py`

- `_populate_em_annotations`: existing `AF_GDNA_RESOLVED` (=0x3 old) reference must be renamed to `AF_GDNA_EM` (=0x5 new). Same semantic (EM winner on gDNA component), new value.
- `AF_UNRESOLVED` init stays correct for pre-decision rows.
- Update imports.

### 4.4 `src/rigel/native/bam_scanner.cpp`

Four stamp sites need updated constants + ZC strings:

- **Multimap-skip (~L1716)**: currently `zf=0, zc="intergenic"`. Change to `zf=AF_MULTIMAPPER_DROP (0x80), zc="."`. Fixes the latent mis-labeling bug.
- **Main hit loop, has_ann branch (~L1840)**: unchanged — stamps `flags_val` from the annotation table (which now carries new-scheme values).
- **Main hit loop, no-ann branch (~L1855)**: currently `zf=0, zc="intergenic"`. Change to `zf=AF_GDNA_INTERGENIC (0x25), zc="."`. Also set `ZW=1.0` (deterministic).
- **Orphan sweep (~L1883)**: same change as no-ann branch. Represents a raw record left unclaimed by any hit — treat as intergenic (deterministic gDNA).

The constants can be defined at the top of `bam_scanner.cpp` next to `FRAG_*`:

```cpp
static constexpr int AF_UNRESOLVED        = 0x00;
static constexpr int AF_MRNA              = 0x03;
static constexpr int AF_NRNA              = 0x09;
static constexpr int AF_NRNA_SYNTH        = 0x19;
static constexpr int AF_GDNA_EM           = 0x05;
static constexpr int AF_GDNA_INTERGENIC   = 0x25;
static constexpr int AF_CHIMERIC          = 0x40;
static constexpr int AF_MULTIMAPPER_DROP  = 0x80;
```

`frag_class_label` keeps its full enum for internal use but the stamp sites that don't go through the annotation table pass literal `"."` for ZC.

### 4.5 Tests

- **Update** `tests/test_annotate.py`:
  - `test_zf_tag_schema_always_present` / existing tests: update expected ZF values from the old (0, 1, 3, 5, 13) set to the new set.
  - `test_annotated_bam_preserves_records_*`: no change (these test record-count parity, not tag values).
  - `test_zb_tag_*`: no change (orthogonal).
- **Add** `tests/test_annotate.py::test_zf_invariants` covering invariants 1–6 from §2.2: iterate over all records in a produced annotated BAM, assert each ZF value satisfies every invariant.
- **Add** `tests/test_annotate.py::test_zf_intergenic_marked_gdna`: produce a BAM with a read that lands clearly outside any locus; assert the read has `ZF & AF_GDNA_BIT` set and `ZF & AF_INTERGENIC_BIT` set.
- **Add** `tests/test_annotate.py::test_zf_chimeric_marked`: produce a chimeric pair (different refs); assert `ZF == AF_CHIMERIC` and `ZC == "."`.
- **Add** `tests/test_annotate.py::test_zf_multimapper_dropped`: run with `--no-multimap`, inject a multimapper; assert `ZF == AF_MULTIMAPPER_DROP` and `ZC == "."`.
- **Add** `tests/test_annotate.py::test_zc_no_longer_contains_intergenic_or_chimeric`: scan every record's ZC value in fixture outputs; assert the set of observed ZC values is a subset of `{unambig, ambig_same_strand, ambig_opp_strand, multimapper, .}`.
- Golden output regeneration: any golden BAM-derived tests will need `--update-golden`.

### 4.6 Documentation updates (apply after code + tests land)

- `src/rigel/annotate.py` docstring schema table (already required in §4.1).
- `docs/MANUAL.md` §"Annotated BAM" — rewrite ZF description + value table; update ZC enum list.
- `docs/METHODS.md` §12.5 — update tag table ZF entry; update ZC entry.
- `docs/parameters.md` — no tag-list change (tag names unchanged), but refresh the annotated-bam help text if it mentions specific ZF/ZC values.
- `CHANGELOG.md` — add a "Breaking" entry describing the annotated-BAM schema change.

---

## 5. Execution order

1. Edit `src/rigel/annotate.py` — new constants + rewritten `winner_flag`/`winner_flags` + docstring.
2. Edit `src/rigel/scan.py` — chimeric `AF_UNRESOLVED` → `AF_CHIMERIC`.
3. Edit `src/rigel/pipeline.py` — rename `AF_GDNA_RESOLVED` → `AF_GDNA_EM`.
4. Edit `src/rigel/native/bam_scanner.cpp` — new constants, update four stamp sites (multimap-skip, has_ann is unchanged, no-ann, orphan).
5. Recompile: `pip install --no-build-isolation -e .`
6. Update `tests/test_annotate.py` with new expected values + add 4 new tests.
7. Run `pytest tests/test_annotate.py -v` until green.
8. Full regression: `pytest tests/ -v`. If golden tests fail due to tag changes, regenerate with `--update-golden`.
9. Documentation updates (MANUAL, METHODS, parameters, annotate.py docstring, CHANGELOG).
10. Final `pytest tests/` sanity run.

---

## 6. Out of scope

- Chimeric fragments **remain unresolved** in this redesign. Adding them to the EM is a separate (large) feature.
- No change to ZT/ZG/ZR/ZI/ZJ/ZW/ZH/ZN/ZS/ZL/ZB tags.
- No change to the collated-in → collated-out record-parity contract.
- No change to Pass-1/Pass-2 splice-blacklist parity.
- No change to the overall two-pass architecture.
