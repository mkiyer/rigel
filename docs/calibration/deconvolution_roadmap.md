# The three-signal deconvolution — implementation roadmap

**Status:** roadmap / north star, 2026-06-12. Operationalizes the two theory docs
([`deconvolution_generative_model.md`](deconvolution_generative_model.md) — the complete three-species
model; [`deconvolution_implementation.md`](deconvolution_implementation.md) — the capture-aware S0–S5
blueprint) into a phased, plan-then-implement build. This document is the **roadmap**; each phase gets
its own dry-run planning doc + execution, as in the prior redesign.

## 0. Why now — the forensic finding this fixes

The flagship leak (`gdna 3:1, capture-ON, ss=0.99`) is concentrated in **complex loci** — overlapping
opposite-strand genes (e.g. locus 21, `chr_syn:964416-1004165`: GENE0024(+) nested in GENE0023(−)). The
dissection (per-region, per-boundary, oracle-validated) found the failure is **not** a strand bug and
**not** a calibration-prior error in the simple regions. It is specific to **mixed AMBIG regions** — a
region carrying both a real mature transcript and gDNA, where:

- the region is **AMBIG** (overlapping `+`/`−`), so its own strand cannot deconvolve it (`w=0` today);
- under **capture**, the gDNA is enriched *inside* the exon and **depleted at the intron-facing edge**, so
  the unspliced boundary crossing (our gDNA probe) is depleted (region 226: 1,241 crossing vs 14,924
  contained, 12×), while the mature crosses by **splicing** — anchored at the junction, **not** depleted
  (1,944 crossing vs 11,542 contained, 6×);
- the current `splice_junction` collapses these into a **fraction** `ρ_g/(ρ_g+ρ_r)` — a depleted
  numerator over a non-depleted denominator → structural under-call of gDNA (region 226: imputed 0.40 vs
  oracle 0.56), which leaks gDNA → RNA (mostly to the nascent spans).

The fix is to stop collapsing the two boundary signals into one capture-fragile ratio and instead resolve
the region's **own** (capture-enriched) unspliced mass into the full three-term `{RNA+, RNA−, gDNA}`
using three complementary signals.

## 1. The three signals (all three are needed — none is sufficient alone)

A region's contained-unspliced mass is, per the generative model, a mixture: `gDNA (½/½ strand) +
nascent₊ + nascent₋ + mature₊ + mature₋`, where each RNA species reads its transcript's sense at rate κ
(`rna_sense_frac`) and antisense at `1−κ`. We resolve it from:

| # | signal | source | estimates | capture behaviour |
|---|---|---|---|---|
| **1** | unspliced boundary crossings | `density_model` (existing) | **gDNA** local density | depleted at exon edges → under-estimate, but a real, direct probe |
| **2** | **spliced** boundary crossings | **new** | **mature RNA**, *per strand* | junction-anchored → **not** depleted; an accurate absolute estimate |
| **3** | strand deconvolution of the residual | `strand_deconv` (existing, re-aimed) | splits residual **gDNA vs nascent**; informs AMBIG | capture-invariant (a fraction; enrichment cancels) |

**Signal #2 is strand-aware (the key correction).** A junction's splice motif identifies its strand
reliably. A `+` junction's spliced crossings estimate the `+` mature rate `m₊`. But that mature, when it
makes *unspliced* (within-exon) reads, distributes onto the **observed** strand counts by κ: **κ·m₊ to
RNA+, (1−κ)·m₊ to RNA−** (ss=0.9 → 90/10). So the predicted mature is subtracted **strand-aware** — κ·M
from the sense count, (1−κ)·M from the antisense count — preserving the residual's strand structure for
signal #3. A `−` junction contributes `m₋` with the mirrored split.

**Why all three.** #2 alone gives only RNA on the junction's strand; #1 alone gives only (depleted) gDNA;
neither resolves the opposite-strand RNA. Subtracting the **well-estimated** mature (#2) from the region's
**own enriched** total leaves a `gDNA + nascent` residual at the true interior density — already correcting
#1's depletion under-count — and #3 (plus #1 as a prior/lower-bound) splits that residual. They are
independent measurements of the same local latents → combined as **likelihoods**, not collapsed into a
ratio.

## 2. The per-node cascade (one feed-forward iteration — no loop)

Per region, after the library-global fits (κ, strand overdispersions, FL pmfs — unchanged):

1. **Mature rate, per strand** from the bounding eligible splice junctions:
   `m̂_s = (Σ_j S_{j,s}) / (n_{j,s} · fl_mean_RNA)` for `s ∈ {+,−}` (sum the eligible junctions' spliced
   crossings; divide by their count × the RNA crossing eff-length — see §3 geometry).
2. **Predicted contained-unspliced mature, per strand:** `M_s = m̂_s · E_mu`, `E_mu =
   region_eff_length(L, rna_pmf)` (exon-limited contained eff-length).
3. **Strand-aware subtraction** from the observed unspliced counts → residual `U'_pos, U'_neg`
   (`gDNA + nascent₊ + nascent₋`): remove κ·M_s from strand-s sense, (1−κ)·M_s from antisense.
4. **Resolve the residual** with the strand posterior (signal #3) + the gDNA crossing density (signal #1)
   as independent likelihoods on the local `ρ`:
   - **single silent strand** (`m̂` ≈ 0 ⇒ `ν` ≈ 0 — nascent tracks mature for the same gene): the region
     is *effectively single-strand* → strand identifies `ρ` vs the live strand's nascent. **The AMBIG
     annotation is overridden by the data.**
   - **both strands expressed:** residual under-determined by strand alone → the per-strand mature
     expression is the **directional prior** on nascent, blended (the existing `w=I/(I+I₀)` gradient) with
     the gDNA crossing estimate / neighbour propagation.
5. **gDNA mass** = residual − nascent; the spliced mass stays deterministic RNA.

The intertwining the design needs: **signal #2 supplies the per-strand expression that drives signal #3.**
The count and strand channels are no longer independent — mature expression (count side) becomes the prior
on which strand's nascent exists (strand side). Still a single pass: fits → mature subtract → residual
resolve → impute remainder. No EM loop (deferred; §6).

## 3. The geometry — validated, no new constant

The mature imputation mirrors the gDNA density imputation, on spliced crossings. The divisor is the count
of bounding eligible junctions × the RNA crossing eff-length:

```
m̂_s = (Σ eligible-junction spliced crossings on strand s) / (n_junctions_s · fl_mean_RNA)
M_s  = m̂_s · region_eff_length(L, rna_pmf)
```

**Validation (locus 21, oracle):** with the *single* `fl_mean` divisor the mature is over-predicted ~2×
(region 226: predicts 19,951 of 22,372 → residual 0.11 vs oracle 0.59); the back-solved divisor is
≈2·fl_mean — exactly the **two bounding junctions** of an internal exon. With `n_junctions·fl_mean`, the
residual fraction recovers the oracle on the failing regions (226/231/236/240) **and degrades to a no-op**
on the gDNA-dominated regions (sig=6: ~0 spliced crossings → ~0 subtraction → stays at 0.98). The geometry
is pure (reuses `region_eff_length` / `fl_mean`), introduces no tunable constant, and self-scales with
mature abundance.

Phase-A details to nail: terminal vs internal exons (`n_junctions` ∈ {1,2}); junctions whose other side is
a reference edge; and reconciling the **substrate's** fractional crossing-mass counts with true crossings
(substrate reported 3,146 vs oracle 1,944 at one junction — the accumulator's boundary-mass split must be
matched, not assumed).

## 4. Code health — reuse, new, retire (organization is the implementer's call)

This is a **re-base, not a rewrite** — the modules already have the shape. Treat it as an opportunity to
improve health: single-purpose modules, no double-counting, retire the superseded path.

**Reuse:** `strand_deconv.strand_posterior_gdna_frac` (works on any sense/antisense); the `w=I/(I+I₀)`
gradient blend; the boundary-anchored imputation *pattern* (`density_model`); `splice_junction_eligibility`
(already names each boundary's splicing strand); the RNA/gDNA eff-length functions; the substrate's
per-view `n_spliced_sense/antisense`.

**New:** a **mature-imputation module** (the RNA mirror of `density_model`) — per-region per-strand
predicted contained-unspliced mature; the strand-aware subtraction producing `U'_pos/U'_neg`; AMBIG
strand-observability driven by per-strand mature expression.

**Retire / re-base:** `splice_junction.region_splice_gdna_frac`'s fraction-at-boundary (subsumed by the
signal-#2 + residual resolution); `node_gdna_density`'s exon path re-aimed at the residual. Keep the
eligibility predicate.

## 5. Phased roadmap (each: dry-run plan → execute → validate → commit)

- **Phase A — mature imputation geometry (the crux).** Build the mature-imputation module + nail
  `E_mu / (n_junctions·fl_mean)`. Validate the predicted contained-unspliced mature against the locus-21 +
  toy-AMBIG oracle, in isolation, before touching `calibrate`. *This is the only thing that can be wrong;
  prove it first.*
- **Phase B — mature subtraction.** *(plan written, then **merged into C** — see below.)* Step 0 found the
  subtraction only benefits AMBIG nodes (a strand-observable node's strand already separates gDNA from all
  RNA — generative §4), so there is no separable subtraction-only phase. Step-0 record:
  [`phaseB_mature_subtraction_plan.md`](phaseB_mature_subtraction_plan.md).
- **Phase C (merged B+C) — the iterative propagation deconvolution** *(rev 3 — the unifying architecture).*
  The per-node prototype proved a region cannot be solved in isolation, so deconvolution is a **graph of
  region + boundary nodes**, each a 3-term `{RNA+, RNA−, gDNA}` partition, **solved by propagation from
  seeds** (intergenic / single-strand / clean-junction nodes) **to convergence**. Boundaries deconvolve
  their *unspliced* crossings into all three terms too (the key revelation), so information flows
  region↔boundary↔region. The AMBIG gDNA is the **residual of the region's own enriched count minus the
  propagated (non-depleted) RNA** → fixes the capture depletion. Consolidates `density_model` +
  `mature_density` + the strand blend + `run_fill` into one solver; retires the splice fraction, the `w`
  blend, and the rev-2 `ρ`. Plan: [`phaseC_ambig_resolution_plan.md`](phaseC_ambig_resolution_plan.md).
  Scope note: this **lifts the one-pass lock** — bounded iteration (a few sweeps to convergence) is now
  required, not optional.
- **Phase D — validate at scale.** Golden refresh; full gDNA suite + the flagship; confirm the leak drop and
  no regression on simple/unstranded/no-gDNA conditions. Update the authoritative theory doc + CLAUDE.md.

## 6. Scope decisions (locked)

- **One iteration, no loop** — the cascade above. The full joint-ML over a gDNA-homogeneous segment
  (generative §6, profile `m`/`ν`, ML `ρ`) is **deferred** until the one-pass version is correct and
  validated. It must work well before we consider iteration.
- **No capture-class gating** — the imputation is purely **local** (a region's own residual + its
  immediately-bounding boundary crossings, all within ~FL), so it is at the right local enrichment *by
  construction*; nothing global is modelled. We do **not** add an on/off-target gate — the mechanism is
  capture-agnostic. Locality (the nearest local anchor), not a gate, is what keeps an exon from inheriting
  a deep-intron density.
- **Protocol-agnostic strand, parametrized by `κ = ss ∈ [0,1]`** — the tool supports both library
  orientations: forward (read1 on the sense strand, `κ→1`) and reverse (read1 antisense, `κ→0`), plus
  unstranded (`κ→½`). `κ = rna_sense_frac = ss` is fit per-library and used consistently; the strand-aware
  mature split puts `κ` on the transcript's sense and `1−κ` on its antisense. Discriminability `(2κ−1)²`
  is **direction-symmetric**, so only `κ≈½` is uninformative. The design must never hardcode a direction
  (this dataset happens to be reverse, `κ≈0.01`).
- **No magic numbers** — the geometry and the blend are derived; the only knobs are the existing `I₀` and
  any prior strength, discussed before adding.

## 7. Open risks

1. The `E_mu/E_ms` geometry edge cases (§3) — the primary risk; pure geometry, oracle-testable.
2. Both-strands-expressed AMBIG — genuinely under-determined in one pass; relies on the directional prior +
   propagation. Where even that fails, **report low confidence** rather than manufacture a value.
3. Substrate crossing-count semantics vs. true crossings (§3) — must be reconciled in Phase A.
