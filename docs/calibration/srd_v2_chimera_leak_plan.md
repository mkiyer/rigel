# SRD v2 — Chimera Leak Into Calibration Pool: Diagnosis & Plan

> **STATUS: REJECTED / SUPERSEDED (2026-04-28).**
>
> Option B (introducing `CHIMERA_CIS_INTERGENIC` and `CHIMERA_CIS_NOVEL_SJ`
> types) was rejected. Rationale:
>
> - In short-read PE data with `max_frag_length ≤ 1 kb`, a fragment whose
>   two reads sit within that distance on the same contig is, by
>   definition, a contiguous piece of DNA. Calling it a "chimera"
>   because of an unannotated splice or a read-through into intergenic
>   space misdiagnoses the event — the simpler explanations are
>   annotation error, unannotated isoforms / TSSs / TTSs, or unannotated
>   splice junctions.
> - True chimeras (case D in the table below — disjoint exon sets across
>   blocks) remain handled by the existing `detect_chimera()` and are
>   not affected by this rejection.
> - Trying to root-cause every long-genomic-footprint fragment is a
>   losing game: we are at the whim of the gene annotation, and any
>   classification we add will misfire on the next annotation refresh.
>
> **Adopted policy (SRD v2 Phase 7):** drop fragments with
> `genomic_footprint > max_frag_length` from FL training for both the
> RNA and gDNA models. The drop count is reported as
> `CalibrationResult.n_pool_dropped_out_of_range` for telemetry. The
> `n_pool` denominator now reflects the in-range subset that actually
> fed the 1-D mixture EM. No new chimera enum values, no per-transcript
> max-intron lookup, no novel-SJ heuristic.
>
> Implementation: [src/rigel/calibration/_simple.py](../../src/rigel/calibration/_simple.py#L155-L175),
> [src/rigel/frag_length_model.py](../../src/rigel/frag_length_model.py#L149-L165) (already in effect for the per-category models).
>
> The original analysis below is preserved for historical context.

## Problem

In Phase 4 validation we observed ~5–11% of the SRD v2 pool fragments
being dropped as **out-of-range (OOR)** because their `genomic_footprint
> max_size (1000 bp)`. The FL histogram excludes them but they
contaminate `n_pool_categorized`, distorting acceptance metrics.

User assertion (correct): chimeras *should already be detected* and
routed away from calibration.  Why are they bleeding in?

## Root Cause

`detect_chimera()` in [src/rigel/native/constants.h](src/rigel/native/constants.h#L279-L352)
fires only when **≥2 alignment blocks each have a non-empty
`exon_t_set`** that are pairwise disjoint. This catches the canonical
trans-/cis-chimera (R1 in gene A exon, R2 in gene B exon) but **misses
several real-world failure modes** that produce huge `genomic_footprint`
values which then leak into the calibration pool:

| Case | R1 location | R2 location | exon_t_set R1 | exon_t_set R2 | detect_chimera fires? | Current fate |
|------|-------------|-------------|---------------|---------------|----------------------|--------------|
| A | gene A exon | **intergenic** | non-empty | empty | **No** (n=1) | EXON_INCOMPATIBLE → OOR drop |
| B | gene A exon | gene A intron (deep) | non-empty | empty (introns not in exon index) | No | INTRONIC pool — fine |
| C | gene A exon | gene A exon, separated by **un-annotated intron** | non-empty | non-empty, **same t_set** | No (sets not disjoint) | fl ≈ genomic_span → OOR drop |
| D | gene A exon | gene B exon | non-empty | non-empty disjoint | **Yes** ✅ | CHIMERIC, excluded |
| E | gene A exon | adjacent read-through into gene B intergenic gap | non-empty | empty | **No** (n=1) | EXON_INCOMPATIBLE → OOR drop |

Cases **A**, **C**, **E** are the OOR pollution. Together they are
biologically anomalous fragments (true chimeras, novel splicing, or
read-throughs to intergenic) that the user's intent says should be
routed to the chimera path, not into calibration.

## Plan

### Option B (recommended): extend `detect_chimera` to cover **one-sided** and **same-set + large-gap** cases

Add two complementary checks inside `detect_chimera` that operate
**after the existing pairwise-disjoint scan**:

1. **One-sided intergenic chimera** (cases A, E)
   If at least one block has a non-empty `exon_t_set` AND at least one
   other block has an empty `exon_t_set` AND the inter-block gap exceeds
   a threshold `MAX_INTERGENIC_GAP_BP` (default 1000 bp), emit
   `CHIMERA_CIS_INTERGENIC` (new enum, value 4) with `gap = max gap`.

2. **Same-isoform impossibly-large gap** (case C)
   If all non-empty `exon_t_set`s share at least one transcript AND the
   inter-block gap is larger than the **longest annotated intron of any
   shared transcript**, emit `CHIMERA_CIS_NOVEL_SJ` (new enum, value 5).
   Implementation note: we already have transcript spans in the
   resolver context; longest annotated intron per transcript can be
   precomputed at index build time, or approximated by transcript span
   minus spliced length.

Both new types should:
- Be tagged with `AF_CHIMERIC` in the alignment-flags output.
- Set `cr.chimera_type` so they take the existing chimeric branch
  (`FRAG_CHIMERIC`) and are skipped by the scorer (already implemented).
- Increment a new telemetry counter
  (`stats.n_chimera_cis_intergenic`, `stats.n_chimera_cis_novel_sj`)
  for visibility.

### Why not Option A (post-hoc threshold in calibration)

Option A would just drop large-`genomic_footprint` fragments before
the OOR check inside `_simple.py`. Two reasons against:
- The fragments still flow through the resolver, scorer, EM. Wasteful.
- Calibration would silently shoulder a classification responsibility
  that belongs in the chimera detector. Brittle.

### Why not Option C (cap by `1.5 × max_tx_length`)

Tempting but breaks down for true short transcripts where a legitimate
PE insert can exceed 1.5 × tx_length when one end overhangs (which is
already an EXON_INCOMPATIBLE case, not a chimera). Cleaner to keep the
intergenic / novel-SJ semantic distinction.

## Acceptance Criteria

After Option B lands and Phase 4 sweep is re-run:

- OOR drop rate `< 0.5%` of pool (down from 5–11%).
- New telemetry: `n_chimera_cis_intergenic` and `n_chimera_cis_novel_sj`
  account for the displaced fragments.
- gDNA π̂ headline metrics stay within ±0.5 percentage points of current
  Phase 4 numbers (the displaced fragments were anomalies, not signal).
- No regression in `tests/` (golden updates expected for any test that
  inspects chimera_type counts or pool sizes).

## Implementation Steps

1. Add `CHIMERA_CIS_INTERGENIC = 4`, `CHIMERA_CIS_NOVEL_SJ = 5` to the
   `ChimeraType` enum in [src/rigel/native/constants.h](src/rigel/native/constants.h).
2. Extend `detect_chimera` with the two new checks after the existing
   disjoint-set loop. Threshold constants:
   - `MAX_INTERGENIC_GAP_BP = 1000`
   - For novel-SJ: per-transcript max-intron length (loaded from index;
     fallback to `transcript_genomic_span - spliced_length` if not yet
     in the index format).
3. Add telemetry counters in `bam_scanner.cpp` PipelineStats.
4. Recompile, run pytest, update goldens.
5. Re-run Phase 4 sweep at
   `/scratch/mkiyer_root/mkiyer0/shared_data/rigel/srd_v2_vcap_mixture/default/`.
6. Update `docs/calibration/srd_v2_results.md` with corrected numbers.

## Out of Scope (explicitly)

- True novel splice-junction *discovery* (we don't try to find them, we
  just route them away from calibration).
- Multi-mapper chimera handling — multi-mappers already short-circuit
  the chimera path; this plan only changes unique-mapper behavior.
- Index format changes — the per-transcript max-intron lookup can land
  in a follow-up if the fallback (span − spliced_length) is empirically
  good enough.
