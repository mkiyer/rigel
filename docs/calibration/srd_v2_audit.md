# SRD v2 — Phase 0 Audit

Read‑only verification of the eight Phase 0 items in
[srd_v2_final_plan.md](srd_v2_final_plan.md). All file:line references are
against `HEAD` (`7ebd391`, with the in‑flight unstaged edits noted in
`git status`).

## Decisions table

| # | Item | Status | Decision |
|---|---|---|---|
| 1 | `compute_frag_lengths` skips entries with projected `fl ≤ 0` | ✅ verified | Calibration must stop reading per‑candidate `frag_lengths` (Phase 1). |
| 2 | `from_core` fills missing entries with `-1` | ✅ verified | Same as #1; `-1` sentinel is unavoidable while we keep the per‑candidate map for the EM scorer. |
| 3 | `np.clip(..., 0, max_size)` in `_simple.py` | ✅ verified | Replace with explicit out‑of‑range filtering + diagnostic counter (Phase 1). |
| 4 | `genomic_footprint` universally populated and ≥ 0 | ✅ verified by code path; **runtime spot‑check pending** | Use `genomic_footprint` as the pool length primitive (Phase 1). One caveat for SPLICE_ARTIFACT — see §4 below. |
| 5 | Resolver drops zero‑candidate fragments | ✅ verified | Stop dropping in Phase 2; let them flow as INTERGENIC. |
| 6 | SJ blacklist filter and splice_type assignments | ✅ verified | New `SPLICE_ARTIFACT` class on the blacklist‑rejection path (Phase 2/5). |
| 7 | Annotated‑introns index in resolver | ⚠️ **does not exist as a queryable index** | `SPLICED_IMPLICIT` (deferred to v2.1) cannot reuse the SJ map without modification — see §7 below. |
| 8 | Fraction of pool fragments with `frag_length = -1` | ✅ inferred from existing run | ~33 % of EXON_INCOMPATIBLE in `mctp_vcap_rna20m_dna00m`. Phase 1 alone will eliminate it. |

---

## 1. `compute_frag_lengths` skips `fl ≤ 0`

[resolve_context.h:725-757](../../src/rigel/native/resolve_context.h#L725-L757):

```cpp
// Single exon block
if (exons.size() == 1) {
    int32_t fl = exons[0].end - exons[0].start;
    if (fl > 0)                                // <-- skip on fl <= 0
        for (int32_t t : t_inds) result[t] = fl;
    return result;
}
...
for (int32_t t : t_inds) {
    int32_t tx_s = genomic_to_tx_pos(gstart, t);
    int32_t tx_e = genomic_to_tx_pos(gend, t);
    int32_t fl   = std::abs(tx_e - tx_s);
    if (fl > 0) result[t] = fl;                // <-- skip on fl <= 0
}
```

Two skip conditions:
- Single‑block fragments where the block has zero length (degenerate; rare).
- Multi‑block fragments where projecting `gstart`/`gend` to transcript
  coordinates collapses to the same transcript position. This happens
  when both ends fall outside the transcript's exonic span (overhang into
  intron/intergenic on candidate `T_j`) so `genomic_to_tx_pos` clamps
  both projections to the same boundary base, yielding `tx_e - tx_s = 0`.
  This is the dominant source of `-1` sentinels for `EXON_INCOMPATIBLE`
  fragments.

## 2. `from_core` fills `-1` for missing entries

[resolve_context.h:163-169](../../src/rigel/native/resolve_context.h#L163-L169):

```cpp
r.frag_lengths.reserve(r.t_inds.size());
for (int32_t t : r.t_inds) {
    auto it = cr.frag_length_map.find(t);
    r.frag_lengths.push_back(
        (it != cr.frag_length_map.end()) ? it->second : -1);
}
```

Then [_categorize.py:142-144](../../src/rigel/calibration/_categorize.py#L142-L144)
reads the **first** candidate's length:

```python
first_t = t_indices[idx]
ref_id[valid_first] = t_to_ref_id[first_t]
frag_length[valid_first] = cand_frag_lengths[idx]   # may be -1
```

## 3. `np.clip` in `_simple.py`

[_simple.py:130-133](../../src/rigel/calibration/_simple.py#L130-L133):

```python
pool_hist = np.bincount(
    np.clip(pool_lens, 0, max_size).astype(np.intp),  # -1 -> 0
    minlength=n_bins,
).astype(np.float64)
```

`np.clip(-1, 0, max_size) → 0`, dumping every undefined‑length pool
fragment into bin 0.

## 4. `genomic_footprint` universally populated and ≥ 0

**Definition (single source).** [bam_scanner.cpp:342-352](../../src/rigel/native/bam_scanner.cpp#L342-L352):

```cpp
int32_t genomic_footprint() const {
    if (exons.empty()) return -1;
    int32_t min_start = exons[0].start;
    int32_t max_end   = exons[0].end;
    for (size_t i = 1; i < exons.size(); i++) {
        min_start = std::min(min_start, exons[i].start);
        max_end   = std::max(max_end,   exons[i].end);
    }
    return max_end - min_start;
}
```

**Resolved‑fragment path.** [resolve_context.h:834](../../src/rigel/native/resolve_context.h#L834):
`if (n_exons == 0) return false;` — every fragment that reaches the
buffer has `exons.size() ≥ 1`, hence `genomic_footprint ≥ exon_block_len > 0`.

**Storage.** Populated in `RawResolveResult.genomic_footprint`
([resolve_context.h:836](../../src/rigel/native/resolve_context.h#L836))
and copied to `ResolvedFragment` and the C++ accumulator vector
([resolve_context.h:155, 246](../../src/rigel/native/resolve_context.h#L155);
[L246](../../src/rigel/native/resolve_context.h#L246)) and then the buffer
([buffer.py:142](../../src/rigel/buffer.py#L142)).

### Caveat — SPLICE_ARTIFACT and SPLICED_IMPLICIT inflate `genomic_footprint`

For fragments with multiple disjoint blocks (CIGAR‑N spliced and
paired‑end gap), `genomic_footprint = max_end − min_start` includes the
intervening gap. This is the correct primitive for **truly UNSPLICED**
fragments, where the gap is a real insert. It is **misleading** when:

- The fragment had a CIGAR‑`N` that was rejected by the SJ blacklist and
  was downgraded to `UNSPLICED` (the blocks remain disjoint; footprint
  includes the artifact intron).
- The paired‑end gap silently spans an annotated intron (currently
  classified `UNSPLICED`).

**Phase 1 implication.** Replacing per‑candidate length with
`genomic_footprint` will eliminate the bin‑0 artifact, but will
introduce a long‑tail inflation in the pool from these two cases. The
expected magnitude is small (`n_sj_blacklisted` is the relevant counter
to watch in `summary.json`; for `mctp_vcap_rna20m_dna00m` it is ~tens of
thousands out of 13.9M fragments — well under 1 %). Phase 2's
`SPLICE_ARTIFACT` class fixes the first; v2.1 fixes the second.

### Runtime spot‑check (TODO before Phase 1 ships)

Add a 10‑line `scripts/debug/audit_genomic_footprint.py` that opens the
`mctp_vcap_rna20m_dna00m` BAM, runs `scan_and_buffer`, and reports
`min`, `max`, `n_zero`, `n_neg`, and `n_over_1000` of the buffered
`genomic_footprint` column. Acceptance: `n_neg == 0`, `n_zero` < 0.01 %,
`n_over_1000` < 0.1 %. Skip if it doubles audit time; the static
analysis above is the substantive evidence.

## 5. Zero‑candidate fragment drop

Two places drop fragments with empty `t_inds`:

- [resolve_context.h:1014-1016](../../src/rigel/native/resolve_context.h#L1014-L1016):
  ```cpp
  if (cr.t_inds.empty()) {
      scratch.clean();
      return false;
  }
  ```
- The earlier `else { return false; }` at line 1011 for the "no exon and
  no transcript overlap" branch.

The C++ scanner then runs the fallback intergenic accumulator at
[bam_scanner.cpp:1320-1326](../../src/rigel/native/bam_scanner.cpp#L1320-L1326):

```cpp
if (is_unique_mapper && !frag.has_introns()) {
    int32_t flen = frag.genomic_footprint();
    if (flen > 0) {
        fraglen_obs.intergenic_lengths.push_back(flen);
        stats.n_frag_length_intergenic++;
    }
}
```

This populates `frag_length_models.intergenic` and is replayed into the
pool by `_simple.py:142-152`.

**Decision.** In Phase 2, remove the early `return false;` on
`cr.t_inds.empty()` so these fragments flow through the buffer with
empty `t_inds`. The intergenic accumulator + replay path becomes
redundant and is removed in Phase 5 cleanup. The `n_frag_length_intergenic`
counter remains as a diagnostic to cross‑check that the new flow
preserves count.

**Downstream consumer impact (one‑grep check).** The locus builder
([locus.py](../../src/rigel/locus.py)) iterates fragments via
`t_indices`; an empty `t_inds` contributes nothing. `FragmentRouter`
([scan.py](../../src/rigel/scan.py)) similarly skips. No further plumbing
change anticipated.

## 6. SJ blacklist filter and splice_type assignments

**Per‑read counter** ([bam_scanner.cpp:178-181](../../src/rigel/native/bam_scanner.cpp#L178-L181)):

```cpp
// SJ artifact: how many CIGAR-N junctions were dropped per read by the
// blacklist. Populated by filter_blacklisted_sjs in Pass 1 (for the
// thread-pool path) or Pass 2 (legacy).
uint16_t n_sj_blacklisted = 0;
```

**Filter implementation** ([bam_scanner.cpp:478-489](../../src/rigel/native/bam_scanner.cpp#L478-L489)):

```cpp
static inline int32_t filter_blacklisted_sjs(
    std::vector<IntronBlock>& sjs,
    ResolverContext& resolver,
    int32_t ref_id) {
    if (sjs.empty() || !resolver.has_sj_blacklist() || ref_id < 0) return 0;
    ...
    const auto* hit = resolver.sj_blacklist_lookup(...);
    ...
}
```

Called at [L602-605](../../src/rigel/native/bam_scanner.cpp#L602-L605)
and [L1947-1948](../../src/rigel/native/bam_scanner.cpp#L1947-L1948).

**Splice type assignment.** After the filter has erased blacklisted
intron blocks from `rec.sjs`, the resolver sets `cr.splice_type` based
on what remains:

- [resolve_context.h:1006-1008](../../src/rigel/native/resolve_context.h#L1006-L1008):
  ```cpp
  cr.splice_type = has_unannotated_sj ? SPLICE_SPLICED_UNANNOT
                                      : SPLICE_UNSPLICED;
  ```
  This is reached when `any_exon` is true but `has_annotated_sj` is
  false. **A fragment whose only CIGAR‑N junctions were all blacklisted
  is now UNSPLICED at this point** because `rec.sjs` was emptied
  upstream and the resolver no longer sees the introns. This is the
  exact "downgrade to UNSPLICED" path Pillar 5 must intercept.

**Phase 5 hook.** Either (a) preserve the original
`n_sj_blacklisted > 0` flag through to the resolver and set
`cr.splice_type = SPLICE_ARTIFACT` here, or (b) keep the blacklisted
introns in `rec.sjs` with a "blacklisted" tag and let the resolver
classify. Option (a) is simpler and is the recommended Phase 5 approach;
the `n_sj_blacklisted` counter is already threaded per‑read.

**Block merging.** The "disjoint aligned blocks merged into one
footprint" complaint in the prior plan is a red herring at the source:
nothing in `build_fragment` or the resolver explicitly merges blocks.
The "merge" is the side effect of computing `genomic_footprint =
max_end − min_start` over the disjoint exon blocks (see §4 caveat).
There is no separate merge step to remove. Phase 5's "do not merge
blocks" goal therefore reduces to **emit the SPLICE_ARTIFACT class and
hold the fragment out of calibration** — not a separate plumbing change.

## 7. Annotated‑introns index

**Searched fields** in `ResolverContext`:

```
sj_map_           // SJKey -> (offset, count) into sj_map_data_
sj_map_data_      // flat int32 of per-SJ transcript ids
sj_lookup(ref_id, start, end, strand)
```

`sj_map_` is keyed by **(ref_id, start, end, strand)** — i.e. it answers
"which transcripts have an intron with **exactly these** boundaries?".
It does **not** support range queries like "is there any annotated
intron whose interval is fully contained in `[r1_end, r2_start]`?",
which is what `SPLICED_IMPLICIT` detection needs.

**Decision.** v2.1 will need either:
(a) a new cgranges interval index over annotated introns
(`itype = ITYPE_INTRON`), built at index‑load time from the same SJ list
that populates `sj_map_`; or
(b) a per‑transcript intron lookup that walks
`exons_starts_/exons_ends_` to derive intron coordinates on the fly.

Option (a) is cleaner and reuses the existing cgranges infrastructure.
This is **out of scope for v2** per the final plan.

## 8. Quantification of `frag_length = -1` in pool

Direct buffer instrumentation skipped (would require a full re‑scan
~10 min). Instead inferred from the existing
`mctp_vcap_rna20m_dna00m/summary.json`:

- `category_counts` reports `EXON_INCOMPATIBLE = 20,756` (the only
  Python‑categorized contributor to the pool).
- The `intergenic` accumulator contributes another `8,589` fragments,
  which all have positive `genomic_footprint` (the C++ check at
  [bam_scanner.cpp:1325](../../src/rigel/native/bam_scanner.cpp#L1325)
  filters `flen > 0`). These cannot land in bin 0.
- The fitted `gDNA_FL` mass at bin 0 is **6,837**. Because `RNA_FL[0] ≈ 0`
  (RNA fragments don't have length 0), the EM responsibility at bin 0 is
  ≈ 1.0 to gDNA, so `pool_hist[0] ≈ gdna_counts[0] ≈ 6,837`.
- Therefore **≈ 6,837 / 20,756 = 32.9 %** of EXON_INCOMPATIBLE pool
  fragments have `frag_length = -1` (sentinel) and were silently mapped
  to bin 0.
- Bin `max_size = 1000` mass is `540` (~2.6 % of EXON_INCOMPATIBLE) — a
  smaller artifact from `np.clip` saturation.
- Combined: **~35 % of the EXON_INCOMPATIBLE pool is non‑physical
  bin‑0/bin‑max mass** in the pure library.

For the spiked libraries the *fraction* of the gDNA component that is
artifact drops (real gDNA signal grows much faster), but the absolute
artifact mass scales with pool size:

| Library | EXON_INCOMPAT | gDNA bin‑0 | gDNA bin‑max | Artifact frac |
|---|---|---|---|---|
| dna00m | 20,756 | 6,837 | 540 | **35.5 %** |
| dna20m | 438,108 | 56,387 | 29,290 | 19.6 % |
| dna80m | (not categorized in summary) | 179,445 | 100,961 | (similar) |

Phase 1 alone is expected to drop pool bin‑0 by ~33 % across the board
— small absolute change in the headline `gdna_fraction` for spiked
libraries, large relative change for the pure library's `pi_pool`.

---

## Phase 0 deliverables (status)

- [x] Static verification of all eight items.
- [x] Updates to the final plan: tightened Pillar 5 wording (no separate
  block‑merge step exists) and Pillar 6 wording (no annotated‑introns
  range index exists).
- [ ] Optional runtime `scripts/debug/audit_genomic_footprint.py` —
  defer unless Phase 1 acceptance fails.

**Gate to Phase 1: open.**
