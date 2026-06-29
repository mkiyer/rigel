# Hybrid-Capture Calibration — Derivation

Companion to `docs/mydocs/calibration_fit_idea.md` (the working overview). This is the derivation we
converged on: how to deconvolve gDNA-vs-RNA **elegantly under hybrid capture**, in the log-density
solver, with the goal of shipping once capture behaves.

---

## 1. The factorization

> **`ρ_g(node) = floor × E(node)`**

- **`floor`** — the depleted, off-target gDNA density. Estimated as the mean ρ_g over **intergenic +
  intron** regions (the user's "global gDNA minimum, not mean"). Anchors `E = 1`.
- **`E(node)`** — the **enrichment ratio** (the probe-targeting field). `E = 1` off-target; `E ∈
  [10, 1000]×` on-target; a continuous spectrum (probes capture at different strengths).

This is the right model because hybrid capture targets a **subset** of transcripts (~20K of ~75K
genes): gDNA is ~uniform genomically at the floor, *except* on probe-targeted regions where it is
enriched. The same probe enriches RNA and gDNA from a targeted region by ~the same `E` — so gDNA, with
its uniform genomic baseline, is the clean readout of `E` (RNA's 10,000× dynamic range hides it).

## 2. Grounded structure (measured, `scripts/debug/enrichment_structure.py`)

Synthetic flagship (cap-ON ss99) — `E = ρ_g / floor` by class, floor = 0.0134/bp:

| class | E median | E p10–p90 | reading |
|---|---|---|---|
| intergenic | 1.0 | 0.9–1.1 | **tight floor** |
| intron | 1.0 | 0.8–1.4 | tight floor |
| exon single-strand | **1153×** | 0.8–2163 | **spans the spectrum** (E ≈ 1 → 7700) |
| exon AMBIG | **1212×** | 1.0–3002 | enriched, not strand-solvable |

Off-capture: **every class E ≈ 1** (unimodal, no enrichment). So `E` is **bimodal under capture,
collapses to ≡1 off-capture** — and the regime distinction lives in the *gDNA*.

## 3. Who observes E? (the training set = the user's "self-solving nodes")

| node | observes `E` via | precision |
|---|---|---|
| seeds (intergenic, intron, their exon boundaries) | structure ⇒ `E=1` | infinite |
| single-strand exons | strand deconvolution ⇒ `E = ρ_g/floor` | high (stranded) / zero (unstranded) |
| **AMBIG exons** | the enrichment transfer **`ê(z)`** trained on single-strand exons | `σ²_level(ê)` |

The single-strand exons **span the full E spectrum** (1 → 7700), so they are ideal training data for
`ê`. AMBIG exons sit at the same enriched level but cannot be strand-solved and **cannot be imputed
from their depleted (E≈1) seams** — so they require the enrichment prior.

## 4. The regime unification (resolves "capture: WIP")

Belief-propagation **messages impute `ρ_g` assuming spatial continuity**. That is valid **only where E
is smooth**:

- **Off-capture**: `E ≡ 1` everywhere ⇒ ρ_g smooth ⇒ messages impute it directly (and we measured the
  log-density solver *beats* the lattice off-capture — the model is right there).
- **On-capture**: `E` jumps ~1000× at the **enrichment seam** (exon body vs intron/intergenic). The
  message across a seam is a **density discontinuity** — it must **self-silence**, and the enriched
  body density comes from `ê(z)` instead.

No capture detector. The **disagreement** (source ρ_g vs destination ρ_g) drives the message variance;
`σ²_level` drives `ê`'s pin. One machinery, both regimes, continuous.

## 5. The open problem = message propagation (the suspected gap)

Diagnosed at the node level (`density_message_dissect.py`): under capture the gDNA messages drag
enriched exons DOWN, because the flanking seams impute their depleted density (~6× diluted below the
body). The **σ²_bio var~mean curve under-estimates the seam discontinuity** — it is a *single* curve
averaging a **bimodal edge population** (smooth within-regime + discontinuous across-seam), so it lands
in the middle and the seam message stays confidently wrong. In linear space the bounded fraction-pull
masked this; in log space the far-off mode drags hard.

**This is the gap to close**, and it is *message propagation*: the σ²_bio (message reliability) and
`ê(z)` (the enrichment transfer) are the two fits to get right.

## 6. Prototyping plan (on the REAL VCaP dataset — now unblocked by the log-density solver)

The logodds solver (O(m·K)) is what lets VCaP calibrate at all (the lattice OOM'd at >100 GB). We
extract the real training data pre-sweep (`scratchpad/vcap_extract.py`) and prototype offline:

1. **`ê(z)` enrichment transfer** — `vcap_ehat_points.tsv` (z = seam crossing density, ρ_g = body
   gDNA density, on single-strand exons spanning the spectrum). Candidates: the current monotone
   P-spline; an anchored Hill curve `floor + (ρ_max−floor)·zⁿ/(Kⁿ+zⁿ)` (smooth, monotone, collapses to
   `floor` off-cap); log-log linear; isotonic+smooth. Metric: real-data fit quality + the off-cap
   collapse + how tightly the seam z predicts the body (does the transfer even hold?).
2. **Message `σ²_bio`** — `vcap_gdna_edges.tsv` (adjacent gDNA density pairs by edge type). Confirm the
   bimodal `Δlog ρ` structure (within-regime vs seam) and prototype a **disagreement-aware**
   reliability so seam messages self-silence.

Goals: **simple, elegant, accurate.** Ship the log-density solver once capture reaches parity.

## 7. Checklist (current gaps → elegant capture)

1. **Floor, not mean** — anchor the global at the depleted floor; `−L` of the log-odds grid maps to it
   (the third-party critique #1: de-magics `L`, preserves resolution where data sits).
2. **`ê` anchored to the floor, predicting `E`** — not shrunk toward the enriched mean.
3. **`ê` strong enough on AMBIG exons** — honest `σ²_level` (measured too weak: pulls to 0.42 vs target
   0.68).
4. **Messages self-silence at E-discontinuities** — disagreement-aware `σ²_bio`.
