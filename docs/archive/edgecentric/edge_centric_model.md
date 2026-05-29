Short answer: the symmetry isn't what you think it is. The two slots are not duplicates — they store complementary portions of the same events, not the same data. But your *instinct* that the current layout is awkward is right, and there's a cleaner design.

## What the current store actually holds

For a fragment block spanning the r | r+1 junction at split position x, with extent [a, b]:

- `regions[r].boundary_right` gets mass ∝ `(x − a)` bp inside region r
- `regions[r+1].boundary_left` gets mass ∝ `(b − x)` bp inside region r+1

(Both are multiplied by the same `inv_total = w / fragment_total_bp` weight, and both fall into the boundary FL pool of their *own* coarse class — which may differ across the edge, e.g. exon|intron.) Verified in accumulator.cpp.

So per junction:
- **Event identity** is duplicated — every junction-crossing fragment touches both slots.
- **Mass values** are not. They partition the fragment, so generally `regions[r].boundary_right ≠ regions[r-1].boundary_left`.
- What *is* invariant: their sum (the fragment's mass) and, in aggregate, the count of crossing events.

The "fully spans" case (fragment longer than the whole region) is handled by depositing `w/2` into each of that region's left and right slots, which is a separate quirk — and is also more naturally expressed in an edge view as "fragment generated two distinct edge-crossing events".

## So is the storage redundant?

- In the *information* sense: no. You cannot recover `(boundary_right[r], boundary_left[r+1])` from any single per-edge aggregate, because the bp-split is a fragment-level quantity that was integrated out region-by-region.
- In the *conceptual* sense: yes. Both slots describe the same physical object (the edge between r and r+1), they share an FL-pool identity at the event level, and downstream prior assembly always wants "the flux across this junction" — never "boundary_left of region 7" in isolation.

## A cleaner architecture: edge-centric storage

Index boundary statistics by **edge** rather than by region-side. The edge set is well-defined:

- Enumerate adjacent region pairs within each reference (n_edges ≤ n_regions − n_refs).
- Each edge `e = (r, r+1)` carries its own channel/FL counters.

Per edge, store the *same information content* as today but with explicit semantics:

```
edge_counts[e, channel] = total boundary mass crossing edge e
edge_left_bp_share[e]   = mass-weighted bp share that landed in r
                          (so mass on r+1 side = total − r side)
```

Equivalently: two numbers per edge (`mass_left`, `mass_right`), which equals today's two slots but keyed on edges.

Benefits:
1. **One conceptual object per junction.** No "are these symmetric?" footgun. Summing boundary mass across regions can't accidentally double-count.
2. **Natural fit for prior assembly.** Phase 4 wants "gDNA flux across an intron/exon junction" — that's literally an edge property today, reconstructed from two regions.
3. **Class-pair FL pools become first-class.** The boundary FL pool currently picks "the receiving region's coarse class", which is asymmetric for an exon|intron junction (left-side mass lands in EXONIC_BOUNDARY, right-side in INTRONIC_BOUNDARY). An edge can carry a pair-tag (e.g., `EXON↔INTRON`) and split FL contributions deterministically.
4. **`tolerance_bp` slots in cleanly.** Tolerance becomes a per-edge split-rule (snap |b − x| ≤ k or |x − a| ≤ k entirely to one side). The two-number edge schema is exactly what asymmetric tolerance needs; today's region-keyed store works too, but the edge view makes the rule local.
5. **Memory roughly equal.** ~2 floats × n_edges vs. 2 × n_regions × n_strands × n_splice today — comparable, slightly smaller because we drop the per-region boundary channels.

Costs:
1. **C++ accumulator change.** The hot path addresses regions; you'd need an `edge_id` for each adjacent-region pair (or look up `r → right_edge_id` and `r → left_edge_id`). Cheap, but it's a touch in the inner loop.
2. **Cross-ref handling.** Edges must not span chromosomes; the enumeration is per-ref CSR, which you already have in `RegionArrays.ref_offsets`.
3. **Migration.** All downstream consumers (`PayloadArrays.boundary_*_unspliced_*`, exposure helpers, anything that projects boundary mass back to "this region") need a translation layer. The translation is mechanical: `regions[r].boundary_right = edges[right_edge_of_r].mass_left`, etc.
4. **Contained vs. fully-spans.** Single-region references have no edges, and "fully spans" generates two edge events (left edge + right edge); both work fine, but it's a small semantic shift to verify in tests.

## My recommendation

- Don't collapse to "store once". The current two slots are *not* duplicates; collapsing loses the bp-split.
- Do migrate to edge-centric storage *when* you start touching boundary semantics for tolerance or for Phase 4 prior assembly. Doing it now in isolation is a refactor with no behavioral payoff; doing it together with `tolerance_bp` or the boundary density estimator gets you a cleaner accumulator, cleaner invariants, and a natural place to hang the split-rule.
- Even before that migration, you can get most of the readability win in Python: expose a derived `EdgeView` over `PayloadArrays` that computes `(mass_left, mass_right, total)` per `(r, r+1)` pair. Zero-copy, no schema churn, and it sets up the API the downstream estimators will eventually want — so the C++ migration becomes a transparent backend swap.

If you want, I can sketch the `EdgeView` (pure-Python, derived from the existing sorted `RegionArrays` + `PayloadArrays`) as a no-risk first step.