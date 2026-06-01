# PR 4c — Resurrect the gDNA fragment-length distribution

**Parent plan:** [`../00_implementation_plan.md`](../00_implementation_plan.md) §4 (effective length), §5.
**Builds on / fixes:** PR 4b ([`PR04b_ambig_sweep.md`](PR04b_ambig_sweep.md)) — replaces its **interim** gDNA FL.
**Type:** **C++** (accumulator FL-pool emission) **+ Python**. **Build required: yes.**
**Status:** **DESIGN FINAL — awaiting go-ahead.** Decisions ①–⑤ resolved (§6).
Audit of `archive/calibration_legacy_2026_05` + the burnt `fc96902` accumulator
below.

---

## 0. Problem & audit

**Problem.** PR 4b's effective lengths need the **gDNA** fragment-length (FL)
distribution. The pipeline currently passes `category_models[UNSPLICED].pmf` —
which is gDNA **+ nascent RNA**, an over-broad proxy. Determining the gDNA FL is
itself a primary calibration goal; the code that did it was burnt.

**Audit (what survives, what's burnt):**

| Component | Where | Verdict |
|---|---|---|
| `fl.py` `build_fl_models` — EB policy (good/weak/fallback), `FLModels` | `archive/calibration_legacy_2026_05` | **reusable** (slim — drop the FL-*scoring* surfaces, see §1) |
| `_fl_sources` — `global`/`rna`/`gdna` count-vector extractors | same archive | **reusable** (interface) |
| `fractional_evidence.gdna_fl_mass` — gDNA pool = INTERGENIC + INTRONIC, both compartments | `fc96902` (git) | **reusable** (the pool definition — exactly the intended set) |
| accumulator `fl_pool_mass[6][L]` emission (3 region-types × 2 compartments) | `fc96902` (git, burnt) | **rebuild in C++** (current PR 2.5 accumulator emits counts + boundary flux/mass, no FL pools) |
| SRD-v1 `_fl_mixture` (1-D mixture EM demix) | `5f4754f^` (git) | **optional** alternative (§4) |

So: the **Python surface is salvageable** (`fl.py` + `_fl_sources` +
`gdna_fl_mass`); the **missing piece is the C++ FL-pool histogram**.

---

## 1. Theory — the gDNA FL pool

gDNA FL is the empirical FL distribution of fragments drawn from **gDNA-dominated
territory**, which the annotation defines and the burnt code already encoded:

> **gDNA pool = INTERGENIC + INTRONIC regions, contained *and* boundary-crossing
> fragments.** Exonic pools are excluded (they carry mature mRNA). This is your
> set: intergenic regions, intronic regions, intergenic-exon crossings,
> intron-exon crossings.

The pool's FL histogram is finalised by **smooth empirical-Bayes shrinkage**
toward the global FL — **no quality cliff** (decision 5: "quality is a spectrum;
4999 fragments is still good"):

```
fl_pmf = (pool_counts + ρ_ess · global_pmf) / (pool_total + ρ_ess)
```

A large pool (`pool_total ≫ ρ_ess`) is essentially empirical; a sparse pool
relaxes continuously toward the global FL; an empty pool *is* the global FL — one
smooth curve, the single Dirichlet pseudo-count `ρ_ess` its only knob. No
`GOOD`/`WEAK`/`fallback` classification, no `≥5000` threshold, no adaptive cap.
Residual **nascent RNA** in introns is the one known contaminant — handled at
base by excluding exonic mass and (optionally) by the iterative refinement of §4.

**Both gDNA and RNA FL distributions are produced** (decision 4): the geometric
effective lengths are FL-agnostic (PR 4b), so each species uses its own FL —
gDNA FL for the gDNA exposure (consumed now), **RNA FL** (= the scanner's
SPLICED-ANNOT histogram, EB-finalised the same way) for the **RNA effective
length in PR 5's live count channel** (the fairness flag raised in PR 4b). PR 4c
delivers both as global FL-anchored distributions; only the gDNA FL has a
consumer today.

> **This is NOT re-introducing FL as a likelihood channel.** The archive README
> and `docs/accumulator/audit_phase1.md` #6 forbid FL *as a per-region
> calibration likelihood* (the burnt per-region FL-mixture, for substrate-memory
> + goal-directedness). PR 4c produces **global FL distributions** used only for
> the **geometric effective lengths** (PR 4b: boundary `μ_FL`, region
> `E_f[max(0,L−ℓ)]`) — no per-region FL histograms in the substrate, no FL
> likelihood term (decision 4: "for now we are not using FL as a likelihood
> term"). We therefore **drop `fl.py`'s scoring surfaces**
> (`rna_scoring`/`gdna_scoring`/contrast) — the downstream EM's FL-likelihood
> machinery, out of scope here.

---

## 2. The FL-pool scheme (C++ → payload)

The accumulator bins each **unspliced** fragment's FL mass (fractional per side,
as for region/boundary mass) into one of **6 pools** = {EXONIC, INTRONIC,
INTERGENIC} × {CONTAINED, BOUNDARY}, by the region signature it overlaps:

```
fl_pool_mass : float64[6][max_size + 1]      # emitted in the payload
gdna pool    = INTERGENIC_CONTAINED + INTERGENIC_BOUNDARY
             + INTRONIC_CONTAINED  + INTRONIC_BOUNDARY      # gdna_fl_mass(payload)
rna  FL      = scanner SPLICED_ANNOT histogram               # spliced ⇒ pure RNA
global FL    = scanner global histogram                      # EB anchor
```

Only 6 × (max_size+1) floats — pooled, not per-region, so the substrate-memory
concern that motivated dropping per-region FL does **not** apply. Spliced
fragments are excluded from the pools (their genomic span ≠ FL); the RNA FL comes
from the scanner's SPLICED_ANNOT channel, as before.

---

## 3. Implementation plan

### 3.1 C++ — accumulator FL pools (the rebuild)
- `accumulator.{h,cpp}`: add `fl_pool_mass[6][max_size+1]`. On each unspliced
  fragment deposit, add its fractional FL mass to the pool(s) for the region
  signature(s) it overlaps — CONTAINED for the contained portion, BOUNDARY for a
  crossing portion (reusing PR 2.5's per-side split). `merge_from` sums pools.
- `bam_scanner.cpp`: surface `fl_pool_mass` in the standalone Accumulator binding
  and `build_result`'s `calibration` dict; thread `max_frag_length`.
- `scan_payload.py`: `AccumulatorPayload` gains `fl_pool_mass: float64[6, L]`;
  `from_scan_result` reshapes + validates.

### 3.2 Python — derivation (resurrect, slim)
- `fl.py` (NEW, from the archive, **slimmed**): `build_fl_models(global, rna,
  gdna, max_size, prior_ess)` → a small `FLModels(global_, rna, gdna, n_*)`, each
  the **smooth-EB** `FragmentLengthModel` (§1; global_ is the un-shrunk anchor).
  Drop the scoring-surface machinery + the `Quality` classification (§1, §5).
- FL-pool constants + `gdna_fl_mass(payload)` (intergenic+intronic, both
  compartments): a small helper (in `fl.py` or `signature.py`).
- `_fl_sources` equivalents: `global`/`rna` from the scanner `FragmentLengthModels`,
  `gdna` from `gdna_fl_mass(payload)`.

### 3.3 Wire into calibrate
- `pipeline.py`: build the FL models via `build_fl_models(...)`, pass the **gDNA
  FL pmf** as `gdna_fl_pmf` to `calibrate` — replacing the interim `UNSPLICED.pmf`.
  The PR 4b effective-length/sweep code is unchanged (it already consumes a pmf).
  The **RNA FL pmf** is produced + logged now and threaded into `calibrate` when
  PR 5 adds the RNA effective length to the live count channel.

### 3.4 Tests
- C++ accumulator FL-pool spec test (deposit → expected pool masses; contained vs
  boundary split; spliced excluded) + native==reference cross-check.
- `build_gdna_fl`: good/weak/fallback branches; gDNA pool excludes exonic;
  empty → global; pmf normalised. Reuse archived test patterns where applicable.
- Pipeline: `gdna_fl_pmf` is the derived gDNA FL (mean distinct from UNSPLICED on
  a synthetic with nascent-rich introns); full-suite failure mode unchanged
  (only `quant_from_buffer` NotImplementedError).

---

## 4. Iterative refinement & the demix option

- **Iterative (your "part of the EM") — recommend PR 5 coupling.** Once the outer
  loop deconvolves regions/boundaries and confirms gDNA-high / RNA-low ones,
  re-derive the gDNA FL weighting the pool by the confirmed gDNA mass (anchored
  empirical Bayes; design docs recovered to `/tmp/srd_fl/anchored_*.md`). PR 4c
  delivers the **one-shot** base FL (a clean fixed point to iterate from); the
  re-derivation hooks into PR 5's outer loop.
- **Mixture-EM demix (SRD-v1 `_fl_mixture`) — optional, benchmark-decided.** A
  1-D EM `pool = π·gDNA + (1−π)·RNA` with RNA FL fixed (spliced) demixes residual
  nascent RNA explicitly. The v6 surface dropped it (mask-purity + EB worked); I'd
  keep it out unless intronic-nascent contamination shows up in benchmarks.

---

## 5. Constants (Q6)

**One** new tunable: `ρ_ess` (the global-FL Dirichlet pseudo-count for the smooth
EB shrinkage, §1). The archive's `GOOD_THRESHOLD = 5000` / `WEAK_THRESHOLD = 1`
quality cliffs and the adaptive cap are **removed** (decision 5 — smooth shrinkage
has no threshold); the scoring-surface constants are **dropped** with the scoring
surfaces (decision 4). `ρ_ess` is a Dirichlet strength, not a cliff: at
`pool_total ≫ ρ_ess` it vanishes; at `pool_total ≪ ρ_ess` the FL is the global
anchor. Default ~`1000` (revisit on real-data pool sizes — decision 5).

---

## 6. Decisions (resolved)

① C++ FL-pool emission ✓ · ② no demix (mask-pool + smooth EB; mixture-EM stays a
flagged §4 option) ✓ · ③ one-shot base here, iterative refinement → PR 5 ✓ · ④
produce **gDNA + RNA** FL distributions for effective lengths, **no FL likelihood
term**, scoring surfaces dropped ✓ · ⑤ **no `GOOD`/`WEAK` cliff** — smooth EB with
a single `ρ_ess` Dirichlet pseudo-count, revisit default on real data ✓.

Design final; ready to implement on your go-ahead.

## Rollback

Revert the sub-PR; `pipeline` reverts to the interim `UNSPLICED.pmf`; PR 4b's
sweep is otherwise unaffected. The C++ accumulator FL pools are additive (no
existing payload field changes).
