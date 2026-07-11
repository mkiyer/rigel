# Measuring the nascent-RNA siphon — the one correct metric

> **TL;DR.** The nascent-RNA siphon is `estimator.em_counts.sum(axis=1)[index.t_df["is_synthetic"]].sum()`
> — the EM-assigned mass on the **synthetic nascent shadow spans**. Use `scripts/debug/_metrics.py::nascent_siphon`.
> Do **not** use `get_counts_df()["count"]`, do **not** filter by `is_nrna`, do **not** use `net_flow` ZF/ZT.

This document exists because, during the 2026-07 flagship-siphon investigation, the siphon was measured
three different ways and each gave a different wrong answer — flip-flopping the diagnosis between "the EM is
structurally biased" and "it's calibration error" until the metric was pinned down. The lesson: **fix the
measurement before interrogating the mechanism.**

## The model: three RNA component classes

The per-locus EM has `n_t + 1` components — one per transcript row plus one gDNA. Transcript rows fall into
three classes, and only one of them is the siphon:

| class | flags | is it siphon? |
|---|---|---|
| **multi-exon mature** | `is_nrna=False, is_synthetic=False` | no — real mature RNA |
| **annotated single-exon** | `is_nrna=True, is_synthetic=False` | **no** — a single-exon transcript is *both* mature and nascent (no intron to splice), so its mass is real |
| **synthetic nascent shadow** | `is_nrna=True, is_synthetic=True` | **YES** — the unspliced twin of a *multi-exon* transcript, added to model nascent RNA; it has no independent molecular existence, so on an `nrna_none` scenario any mass here is pure siphon |

The distinguishing flag is **`is_synthetic`**, not `is_nrna`. `is_nrna` is a superset that includes the
legitimate single-exon transcripts.

## The three traps

1. **`get_counts_df()["count"]` zeroes synthetics.** The display count is
   `np.where(is_synthetic, 0.0, t_total)`, so summing it over nascent rows returns **0** regardless of what
   the EM assigned. (This produced "the controlled toy never siphons".)
2. **Filtering by `is_nrna` conflates real single-exon transcripts with the shadows.** Their mass is real,
   not siphon; including it over-counts and manufactures a false "structural residual". (This produced the
   bogus "84% structural / 43% residual" split.)
3. **`net_flow` ZF-bit (0x08) and ZT hard labels are per-fragment hard MAP labels**, not the EM abundance —
   insensitive to the soft posterior and not comparable across runs (`rigel.sim.analysis`). They answer a
   different question ("how was each fragment hard-labelled") than "what abundance did the EM assign".

## The correct measurement

Read `estimator.em_counts` (indexed over **all** transcripts, including the synthetic rows that
`get_counts_df` drops), sum over channels, and select `is_synthetic`:

```python
from _metrics import nascent_siphon, rna_component_breakdown
siphon = nascent_siphon(estimator, index)                 # the one number
sy, single_exon, mature = rna_component_breakdown(estimator, index)  # the three buckets, never conflated
```

Validate on an `nrna_none` scenario (true nascent = 0), where the siphon is the whole nascent mass. To
separate **calibration error** from a **structural EM** contribution, compare the siphon under the fitted
calibration vs an **oracle calibration** (true per-node masses): `scripts/debug/oracle_lever_test.py`. In the
flagship this showed the siphon is ~100% calibration error (83,703 → 379 under oracle).

## Canonical diagnostic tools

- `scripts/debug/_metrics.py` — the metric definitions (import these; do not re-derive).
- `scripts/debug/oracle_lever_test.py` — siphon under fitted vs oracle calibration + FL/strand/eff-length levers.
- `scripts/debug/ambig_locus_em_proof.py` — controlled toy loci (oracle vs fitted).
- `scripts/debug/locus_component_audit.py` — per-component eff-lengths + masses for the worst locus.
- `scripts/debug/pass_trace.py` — calibration per-node vs oracle deconvolution.
