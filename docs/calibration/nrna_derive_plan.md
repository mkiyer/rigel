# Drop synthetic nRNAs from cgranges; derive on demand from real-tx hits

## Verdict: **YES — implement.** This is a clean win.

## Quantified bloat (Gencode v44 human index)

| metric | count | share |
|---|---:|---:|
| total transcripts | 457,513 | — |
| synthetic nRNAs | 203,052 | **44 %** |
| total cgranges intervals | 2,362,845 | — |
| intervals belonging to synthetics (1 EXON + 1 TRANSCRIPT each) | 406,104 | **17.2 %** |

Synthetic intervals are **long** (entire pre-mRNA span, 5–500 kb), so
they dominate hit lists for every intronic / intragenic query. Removing
them shrinks the index, speeds up `cr_overlap`, and — most importantly
— **eliminates the structural source** of the INTRONIC=0 bug we just
fixed.

## Architectural property that makes this elegant

Every annotated multi-exon transcript carries a one-to-one
`nrna_t_index` field ([src/rigel/transcript.py#L44](src/rigel/transcript.py#L44),
populated in [src/rigel/index.py#L689-L709](src/rigel/index.py#L689-L709))
pointing to its parent nRNA entity (synthetic, annotated nascent-equiv,
or itself). The map is already there — the resolver just doesn't use it.

A nRNA span is *exactly* a single interval `[N.start, N.end)` ⊇ every
contributor mRNA's span. So bp-overlap of a fragment with N is trivial
to compute from the fragment's exon blocks and N's coordinates.

## Proposed change (elegant minimal version)

### Index build ([src/rigel/index.py](src/rigel/index.py))

In `_gen_transcript_intervals`, **skip synthetic transcripts entirely**
(no EXON, no TRANSCRIPT interval emitted). Annotated nascent-equiv
single-exon transcripts continue to emit their own intervals as today
(they're real annotation). The `t_df`, `nrna_t_index`, `is_nrna`,
`is_synthetic` columns are unchanged — synthetics still exist as
transcript rows for scoring/EM, they just no longer pollute cgranges.

One additional thing exposed to the C++ resolver: the
`nrna_t_index` int32 array (parallel to `t_strand`, `t_is_nrna`). Add
a `set_nrna_parent_index(...)` setter alongside the existing
`set_nrna_status(...)`.

### Resolver ([src/rigel/native/resolve_context.h](src/rigel/native/resolve_context.h) `_resolve_core`)

After the per-block exon/transcript hit collection loop, *before* the
classification logic, derive the nRNA hits in one tight pass:

```cpp
// Derive nRNA candidates from real-tx hits.  After this loop, exon_t_sets
// and transcript_t_sets contain only real-transcript indices; nrna_to_visit
// holds the unique nRNA parents whose spans this fragment overlaps.
std::unordered_set<int32_t> nrna_to_visit;
auto add_nrna_parent = [&](int32_t ti) {
    if (ti < 0 || ti >= (int32_t)nrna_parent_.size()) return;
    int32_t n = nrna_parent_[ti];
    if (n >= 0 && n != ti) nrna_to_visit.insert(n);
};
for (auto& s : exon_t_sets)        for (int32_t ti : s) add_nrna_parent(ti);
for (auto& s : transcript_t_sets)  for (int32_t ti : s) add_nrna_parent(ti);

// For each nRNA parent, compute fragment-bp inside its [start,end) span
// and accumulate into scratch.t_exon_bp[N] (single-exon nRNA: all overlap
// is "exonic" relative to N).  scratch.t_transcript_bp[N] gets the same
// value, so the later derivation t_intron_bp = transcript_bp - exon_bp
// correctly yields 0 for nRNAs.
for (int32_t n : nrna_to_visit) {
    int32_t ns = exon_starts_[exon_offsets_[n]];      // single-exon
    int32_t ne = exon_ends_  [exon_offsets_[n]];
    int32_t bp = 0;
    for (const auto& eb : exons) {
        int32_t lo = std::max(eb.start, ns);
        int32_t hi = std::min(eb.end,   ne);
        if (hi > lo) bp += (hi - lo);
    }
    if (bp > 0) {
        scratch.mark_dirty(n);
        scratch.t_exon_bp[n] += bp;
        scratch.t_transcript_bp[n] += bp;
        // Inject n into one of the per-block sets so downstream
        // merge_sets() picks it up.  Block 0 is fine — nRNA hits don't
        // need block-level chimera disambiguation.
        exon_t_sets[0].push_back(n);  // sort happens later via merge
    }
}
```

The `exon_bp_pos/_neg` and `tx_bp_pos/_neg` strand-aware accumulators no
longer need the `is_nrna`-skip guard added in the recent bug fix —
synthetics are simply not in cgranges. Remove that special-case logic
entirely.

### Things that stay unchanged

- `t_is_nrna_` mask (used by scorer to route mRNA vs nRNA scoring) ✓
- `nrna_t_index` column in `t_df` ✓
- nRNA-specific `frag_length` computation (single-exon: tx_pos = genomic_pos − N.start, trivially correct) ✓
- All downstream EM / quant code ✓
- Annotated nascent-equiv transcripts (real single-exon transcripts marked `is_nrna=True, is_synthetic=False`) — they remain in cgranges and behave as today ✓

## Behavioral diff (rare cases)

For a fragment landing in a "tolerance gap" inside a merged nRNA span
that is NOT covered by any contributor mRNA's transcript-span interval:

- **Today:** hits the synthetic EXON → `t_inds = [nrna]` only.
- **After:** no real hit → `t_inds = []` → INTERGENIC.

These cases are vanishingly rare (TSS/TES clustering tolerance is
typically tens of bp, far smaller than any contributor's span). For the
common case (fragment inside an mRNA's genomic span), the new code
yields **richer** candidate sets — the mRNA parents themselves enter
`t_inds` alongside the nRNA, giving the EM more signal.

## Implementation steps

1. `_gen_transcript_intervals`: skip `t.is_synthetic`.
2. C++: add `set_nrna_parent_index` setter on `FragmentResolver`; pass
   `index.t_df["nrna_t_index"].values` from `pipeline.scan_and_buffer`.
3. C++ `_resolve_core`: insert the ~25-line derivation block above.
4. Remove the `is_nrna` skip guards we just added in
   `resolve_context.h:920-961` (no longer needed).
5. Recompile, regen goldens, run pytest.
6. Re-run Phase 4 sweep — expect (a) faster scan, (b) further INTRONIC
   pool growth, (c) some shift in nRNA EM weights as more real-tx
   candidates reach the scorer.

## Acceptance criteria

- All 922 tests pass (after golden regen).
- Phase 4 wall-time decrease ≥ 5 % on dna80m (the largest lib) — bonus
  but expected.
- INTRONIC pool counts in the same ballpark as our just-fixed numbers
  (sanity check that the derivation produces equivalent semantics).
- nRNA quant `gene_quant.feather` numbers stay within ±2 % of post-bug-fix
  values (the architectural change should be score-equivalent up to
  the rare tolerance-gap cases above).

## Risk assessment

Low. The change is mechanically straightforward (~50 LoC), the
abstraction (`nrna_t_index`) already exists, and behavior diverges
only in edge cases that are themselves more correct after the change.
Rollback is trivial — revert one commit and rebuild the index.
