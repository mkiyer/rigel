# Redesign Phase 2 — the count field on strand-cleaned counts (dry-run plan, rev 4 — robust cleaning)

**Status:** dry-run for review, 2026-06-11 (rev 4). Rev-3 used a binary "clean iff `info>0`", which is
**unsafe** when the fitted κ lands just above ½ by sampling noise (the strand split there is the
near-meaningless prior median, and the binary trusts it fully). Rev 4 makes the cleaning **degrade
continuously to a no-op** via the agreed weight `w = I/(I+I₀)`, so the count model can trust the cleaned
counts at *every* strand specificity, κ=½±noise included. Built **alongside**; golden bit-identical until
Phase 3 wires it.

## Goal

Give the count model strand-**cleaned** gDNA counts so the imputed density at exon/AMBIG regions no longer
over-states gDNA by the RNA it can't otherwise see — **and guarantee the cleaning never injects garbage**
when the strand is uninformative. Output `g_count` (count-model gDNA fraction per region) for the Phase-3
combine. Validated by the toy AMBIG node `count_gdna_frac ≈ 0.120 → ~0` at stranded gDNA=0+nascent, and by
a κ-sweep showing graceful (no-cliff) degradation to the raw floor as κ→½.

## The hook (one place, confirmed in code)

`node_gdna_density` builds the whole field from three count arrays — `contained_gdna`, `left_gdna`,
`right_gdna` (raw `pos+neg`, density_model.py:221–223); everything downstream reads only those. One
optional param cleans the field, golden-identical by default:

```python
node_gdna_density(..., gdna_counts=None)   # None → raw (golden-identical) ; provided → strand-cleaned
```

## The cleaning — graceful, robust at every strand specificity

A new helper (the strand module's robust output) turns a Phase-1 `StrandSplit` + the raw count into a
cleaned gDNA count, **shrinking toward "no cleaning" as the strand information vanishes**:

```python
w = info / (info + I₀)                                   # info = N·(2κ−1)² (Phase 1) ; 0 at κ=½, →1 high
cleaned_gdna_count = ( w·g_strand + (1 − w)·1.0 ) · raw_count
```

- The cleaning fraction `w·g_strand + (1−w)·1` slides **continuously**: → **1.0** (all-gDNA = no RNA removed
  = the raw count, a **no-op**) as `info→0`; → `g_strand` (full strand clean) as `info→∞`.
- `g_strand` is `g_q` (the BB posterior read at the FP-quantile), **bounded to [0,1] by the Jeffreys grid**
  — never the wild unbounded linear-unmix MLE — so even when uninformative it is at worst the prior median,
  and `w≈0` discards it.
- `I₀` is the same scale as the Phase-3 combine weight — one knob, sensible default (≈1 effective
  discriminating fragment), validate-don't-theorize.

**Robustness assurance (the key point).** `w` is built from the *same* `info = N·(2κ−1)²` that governs
whether `g_strand` means anything, so garbage and trust are **coupled**: at `κ≈½` with ordinary depth
`info` is tiny ⇒ `g_strand` is the bounded prior median ⇒ `w≈0` ⇒ **no-op, raw passes through**; `g_strand`
is only *applied* (`w→1`) when `info` is genuinely large (high κ, or κ-off-½ with enough depth that the
strand truly resolves it). **There is no strand specificity at which a meaningless split is trusted** — by
construction, not a threshold. κ=0.502, N=100 ⇒ `info≈0.0016`, `w≈0.0016`, cleaned ≈ raw. This is *not* the
rev-2 precision-weighted *field* (the imputation/run-fill is untouched); it is the per-node cleaning
shrinking to a no-op via the agreed `w`.

Applied to all three views (contained for regions; the two boundary sides for the imputation anchors).

## The AMBIG carry-over (no-capture) = the run-fill on cleaned densities

An AMBIG region shares an exon bit with each neighbour, so it has no count-observable boundary and is filled
by `runfill_bidirectional` from the flanking **introns** — now **cleaned** (strand-observable), so it
inherits a nascent-free density. That is the no-capture AMBIG fix (the toy node 5000–6000 fills from the
cleaned introns). **Under capture this run-fill is not enough** — see capture below; that is Phase 3.

## Capture (concern carried from the design): single-strand robust now, AMBIG in Phase 3

- **Single-strand regions are already capture-robust** — the own strand gives the gDNA *fraction*, which is
  probe-enrichment-invariant. No imputation needed.
- **AMBIG regions under capture** need their gDNA imputed from their **own on-target boundary crossings**
  (captured gDNA spilling across the exon edge), via the **per-strand carry-over** (carry nascent⁺/nascent⁻
  from the flanking single-strand exons, subtract, solve for gDNA at the edge) — *not* the off-target intron
  the naive run-fill reaches. This is the boundary architecture's purpose and is capture-correct. **Phase 3
  deliverable**, validated *with* capture — committed, not deferred. Phase 2 validates the **no-capture**
  AMBIG fix first.

## Decisions

- **D1** `node_gdna_density(gdna_counts=None)` hook; `None` ⇒ raw (golden-identical).
- **D2** **Graceful** cleaning `cleaned = (w·g_strand + (1−w))·raw`, `w = info/(info+I₀)` — continuous no-op
  at `info→0`, full clean at `info→∞`. (Replaces rev-3's binary.) New helper in `strand_deconv`.
- **D3** No field precision-weighting; the strand confidence `I` enters the cleaning (above) and the Phase-3
  combine weight. The field imputation/run-fill is unchanged.
- **D4** Robustness is structural: `w` and `g_strand`-meaningfulness share `info`, so no κ trusts a garbage
  split. The Jeffreys grid bounds `g_strand` ∈ [0,1].
- **D5** No-capture AMBIG = run-fill on cleaned densities. Capture-aware AMBIG (per-strand carry-over at the
  exon's own on-target boundaries) = **Phase 3**, committed.
- **D6** Built alongside; live `node_gdna_density` untouched; golden bit-identical.

## Issues / risks

- **I1 — `I₀` is one scale.** Sensible default ~1; the bimodal-data argument says its exact value barely
  matters at the extremes; validate the κ-sweep (esp. the κ≈0.6–0.8 middle) and tune if a weakly-stranded
  condition is sensitive. A sharper `Iⁿ/(Iⁿ+I₀ⁿ)` is available if wanted.
- **I2 — capture+AMBIG is Phase 3, not Phase 2.** Phase 2's run-fill carries the off-target intron density
  into an on-target AMBIG exon — fine without capture, an under-estimate with it. Committed fix in Phase 3
  (per-strand carry-over at the exon's own boundaries). Don't pre-build; do validate no-capture first.
- **I3 — `count_gdna_frac_var`.** Recomputed on cleaned counts; its role in the new combine is a Phase-3
  question. Leave the `need_count_variance` gate as-is.
- **I4 — `g_q` (FP-quantile) flows through the cleaning** (`q≠½` shifts gDNA up/down). Intended.

## Files & changes

- `calibration/strand_deconv.py`: **new** `cleaned_gdna_count(split, raw_count, info_scale)` — the graceful
  shrinkage (+ export). Pure, trivially tested.
- `calibration/density_model.py`: `node_gdna_density(..., gdna_counts=None)`; when provided, set the three
  count arrays from it. **No other change.**
- New config scalar `gdna_strand_info_scale` (`I₀`, default ~1) on `CalibrationConfig`. Wired in Phase 3.
- **No** change to `calibrate.py` / `splice_junction.py`. Golden bit-identical.

## Test plan

- Unit (`cleaned_gdna_count`): `info→0 ⇒ cleaned→raw` (no-op); `info→∞ ⇒ cleaned→g_strand·raw`; continuous;
  `g_strand`∈[0,1]; a κ=0.502, N=100 case ⇒ cleaned ≈ raw (the robustness regression).
- Unit (`node_gdna_density(gdna_counts=None)`) == current (regression).
- **Toy scenario:** `count_gdna_frac[AMBIG]` 0.120→~0 at stranded gDNA=0+nascent; **κ-sweep
  {0.99, 0.8, 0.6, 0.54, 0.50}** showing smooth degradation to the raw floor with **no cliff** at 0.54/0.50
  (the robustness validation); no-nascent / unstranded sane.

## Exit

Full suite green + golden bit-identical (live path untouched); `cleaned_gdna_count` robust across the κ
sweep (no garbage at κ≈½); the toy AMBIG `count_gdna_frac` collapses at stranded gDNA=0+nascent and degrades
gracefully as κ→½; not yet consumed by `calibrate` (Phase 3).
