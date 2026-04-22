# gDNA Locus-Prior Deep Diagnostic

**Date:** 2026-04-22
**Status:** Investigation in progress. No fix to the computation itself yet —
one display-level bug corrected (already staged in working tree).
**Run referenced:** `mctp_vcap_rna20m_dna80m` (dna80m = 80% gDNA contamination sim).
**Expected gdna_rate:** ~0.80. **Observed overall:** ~0.66. **Observed mega-locus:** 0.50.

---

## TL;DR

Three distinct issues are stacked on top of each other in the gDNA
locus-prior / effective-length machinery. In order of magnitude:

| # | Issue | Scope | Impact |
|---|---|---|---|
| D1 | `gdna_prior` column in `loci.feather` is always 0 | **display bug in `get_loci_df`** | zero — EM was unaffected, but led us down wrong diagnostic paths. Already fixed in uncommitted working tree. |
| P1 | Mega-locus: single gDNA state spans 569 Mb across chromosomes via multimapper-linked connected components | **EM model** | dominant (~20–30 pp `gdna_rate` deficit) |
| P2 | gDNA bias denominator uses `gdna_span` (intron-inclusive genomic total), whereas calibration's λ_G is a per-*mappable-bp* rate | **EM–calibration unit mismatch** | small on real data (~5% in mega-locus; median ≈0.1% elsewhere) |

P1 (locus-scale pathology) is the real enemy. P2 is a second-order
inconsistency worth cleaning up.

---

## D1. The `gdna_prior=0` display bug

### Symptom

All 25,396 loci in `loci.feather` report `gdna_prior = 0.0` exactly.

### Root cause

Uncommitted diff against `src/rigel/estimator.py:758` (working tree):

```diff
-                    "gdna_prior": r.get("gdna_init", 0.0),
+                    ...
+                    alpha_g = float(r.get("alpha_gdna", 0.0))
+                    alpha_r = float(r.get("alpha_rna", 0.0))
+                    gdna_prior = alpha_g / (alpha_g + alpha_r) if (alpha_g + alpha_r) > 0 else 0.0
+                    ...
+                    "gdna_prior": gdna_prior,
```

The per-locus result dict `r` (built by `pipeline._build_locus_meta`)
was renamed in an earlier refactor from `{"gdna_init": γ, ...}`
to `{"alpha_gdna": α_G, "alpha_rna": α_R, ...}`. `get_loci_df` was not
updated, so `r.get("gdna_init", 0.0)` returned the default **for every
locus**. The priors themselves were computed correctly and passed to
the C++ EM — only the reported value in the output TSV/Feather was 0.

### Verification

The `alpha_gdna` / `alpha_rna` numpy arrays are constructed by
[`_compute_priors`](../../src/rigel/pipeline.py) and passed to
`_run_locus_em_partitioned`, which uses them as the Dirichlet prior
and as the EM warm-start ratio (`θ_gDNA = α_G/α_R × θ_RNA`). This
path was never broken.

### Action

The working-tree fix is sufficient. After the next run, `loci.feather`
will show true γ values. No further code action required.

---

## P1. Mega-locus pathology

### Observation

| Locus | span (bp) | n_frags | γ_proxy\* | gdna_rate | expected |
|---|---:|---:|---:|---:|---:|
| **0 (mega)** | **569,789,373** | **10.5 M** | **0.279** | **0.500** | 0.80 |
| >10 Mb (n=1) | 569 M | 10.5 M | 0.279 | 0.500 | 0.80 |
| 1–10 Mb (n=31) | — | 0.7 M | 0.49–0.60 | 0.534 | 0.80 |
| 100 kb–1 Mb (n=3301) | — | 24.5 M | 0.19 (med) | 0.631 | 0.80 |
| 10–100 kb (n=9217) | — | 33.9 M | 0.08 (med) | 0.661 | 0.80 |
| <10 kb (n=12067) | — | 5.7 M | 0.18 (med) | 0.753 | 0.80 |

\* γ_proxy = λ_G · mappable_bp_locus / n_em_fragments_locus. **Upper
bound** on the true γ_ℓ (since the denominator `region_n_total` also
counts unique mappers outside the EM pool — see _Caveats_ below).

### Qualitative story

1. **Root of the pathology.** Rigel's locus is a connected component
   of transcripts linked by shared fragments (multimapper union-find).
   High-repeat content in the human genome plus NH-permissive
   alignment (minimap2 with `-N 20`) means chains of multimapper
   fragments collapse into a single giant "locus" spanning chunks of
   many chromosomes. In dna80m we get **one 569 Mb locus consuming
   10.5 M EM fragments** (41% of all EM fragments, 70.5k transcripts).

2. **Effective length blows up.** The C++ EM applies
   `log_liks[i] -= log(prof_len - frag_len + 1)` per candidate
   ([`em_solver.cpp:330`](../../src/rigel/native/em_solver.cpp#L330)).
   For the mega-locus gDNA component, `prof_len = 569,789,373 + gdna_flank ≈ 5.7e8`.
   For a typical RNA transcript in this locus, `prof_len` = spliced
   length ≈ 2 kb. The gDNA per-position probability is therefore
   ~2.9e5× lower than RNA's, or ~12.6 nats of log-likelihood deficit
   per fragment.

3. **The prior can't save it.** γ_mega ≈ 0.28 with `c_base = 5` gives
   `α_G/α_R ≈ 0.39`. That's a warm-start at 28% gDNA, not 80%.
   Even with γ=0.28 the per-fragment log-lik deficit dominates, and
   fragments that actually came from gDNA get pulled toward the
   thousands of low-effective-length RNA components — nRNA in
   particular, since it spans the full genomic range and looks gDNA-ish.

4. **Result.** Mega-locus gdna = 6.34 M of 12.67 M total (observed
   `gdna_rate = 0.500`), but nrna = 2.13 M has absorbed the remainder.
   In truth, for this simulation, nrna_abundance ≈ 0 — so nrna here is
   a pure siphon.

### The deeper issue: what is the "correct" gDNA denominator?

Each fragment has a **finite set of alignment positions** — its NH
hits on the genome. For gDNA, the appropriate "where could this
fragment have originated from" question is "across ALL gDNA positions
genome-wide that are compatible with this alignment," not "across
every bp of the multimapper-linked locus it happens to belong to."

The current model treats the locus's `gdna_span` as the gDNA support
set. This is only accurate when:

- the locus is a single contiguous genomic region (yes for ~small loci)
- the fragment's alignment positions are uniform over that region
  (approximately yes within one locus; false across 569 Mb)

Options (not implementing; for write-up discussion):

- **(A)** Drop the single-gDNA-per-locus state in favor of a
  genome-wide gDNA component parameterized by λ_G directly. This
  sidesteps locus-size dependence entirely. The downside is that it
  decouples gDNA from the isoform EM and complicates fragment
  re-normalization inside the locus.
- **(B)** Partition mega-loci. A locus with span > some threshold
  (e.g. 5 Mb) is split into per-chromosome or per-contiguous-interval
  sub-loci, each with its own gDNA state. Multimapper fragments then
  have one gDNA candidate per sub-locus they touch, proportionally.
- **(C)** Use per-fragment gDNA eff_len: for each fragment, compute
  the total genomic footprint of its NH alignments and use THAT
  (not the locus-wide `gdna_span`) as the gDNA denominator. This is
  correct per fragment but requires threading per-candidate length
  into the scoring/bias correction pipeline instead of the current
  one-size-fits-all-per-locus convention.
- **(D)** Post-hoc: detect mega-loci and re-assign their fragments to
  gDNA via λ_G outside the EM. Hacky; buys headline numbers but
  masks the underlying model issue.

---

## P2. Mappable-bp vs genomic-bp asymmetry (EM ↔ calibration)

### The mismatch

- **Calibration** estimates `λ_G` in units of fragments per
  *mappable bp* ([`_calibrate.py:91`](../../src/rigel/calibration/_calibrate.py#L91)).
  `region_e_gdna = λ_G × mappable_effective_length`.
- **EM bias correction** uses `gdna_span` (total genomic bp,
  intron-inclusive, mappability-unaware).

### Empirical magnitude

From [`scripts/debug/locus_mappable_vs_gdna_span.py`](../../scripts/debug/locus_mappable_vs_gdna_span.py) output:

| Size bin | median `mappable_bp / gdna_span` | Spearman(ratio, gdna_rate) |
|---|---:|---|
| <10 kb | 1.000 | |
| 10–100 kb | 0.999 | |
| 100 kb–1 Mb | 0.999 | |
| 1–10 Mb | 0.998 | overall +0.30 |
| >10 Mb (mega) | 0.946 | |

Impact on mega-locus: extra penalty = log(1/0.946) ≈ 0.055 nats —
dwarfed by the 12.6 nats from P1. Not a primary driver.

### Fix options (not implementing)

- **(B)** Plumb `mappable_bp_locus` (already computable from
  `region_e_gdna` aggregation in `compute_locus_priors`) into the
  C++ EM as the gDNA bias-correction denominator instead of
  `gdna_span`. Preserves λ_G's physical meaning. Requires one extra
  per-locus float array through the EM call path.

---

## Caveats & limitations of this analysis

1. **`γ_proxy` upper-bounds the real γ.** We use `n_em_fragments` as
   a substitute for `Σ region_n_total` within the locus. `region_n`
   counts ALL fragments in overlapping regions, including unique
   mappers that bypass the EM. The real `n_total` is larger, so real
   γ is smaller. Still, for the mega-locus: `n_em_frags = 10.5 M`
   vs global implied `total_n = 77.5 M` (`= E[gdna]/gdna_fraction`).
   If the mega-locus absorbs most of the genome, its `n_total` could
   be up to ~77 M → γ_mega ≈ `2.9M / 77M = 0.04`. Even the **lower**
   bound still leaves the gDNA massively under-weighted by the per-bp
   effective length penalty in P1.

2. **`gdna_flank ≈ 285` bp** is negligible vs the 569 Mb mega span.
   It matters only for small loci.

3. **`c_base = 5`** is the default pseudocount. Increasing it would
   strengthen the prior's pull on the warm start but would NOT
   overcome the per-fragment 12 nats log-lik deficit on the mega-locus.

4. The **non-mega-locus gdna_rate deficit** (0.75 in <10 kb bin vs
   expected 0.80) is not explained by P1 or P2 individually. The
   combined effect of prior pulling toward the locus-local γ
   (median 0.11 in the 10–100 kb bin) plus bias correction is
   consistent with these rates, but a proper decomposition requires
   running the EM with γ forced to its "true" value and comparing.

---

## Scripts

| Path | Purpose |
|---|---|
| [`scripts/debug/locus_mappable_vs_gdna_span.py`](../../scripts/debug/locus_mappable_vs_gdna_span.py) | Correlate `mappable_bp / gdna_span` with `gdna_rate` per locus |
| [`scripts/debug/estimate_locus_gdna_prior.py`](../../scripts/debug/estimate_locus_gdna_prior.py) | Estimate true per-locus γ using `λ_G × mappable_bp / n_frags` proxy |

Both write diagnostic tables under the run's `rigel/` directory.

---

## Next steps (suggested)

1. **Re-run a canonical condition** (e.g. `mctp_vcap_rna20m_dna80m`)
   with the current working-tree `get_loci_df` fix so `loci.feather`
   contains real γ values. Then repeat the correlation analysis
   with the actual priors rather than the `γ_proxy` estimate.
2. **Prototype P1 option (B)** — sub-locus partitioning of mega-loci.
   Hypothesis: splitting the 569 Mb locus into per-chromosome
   components (or smaller) will reduce the per-fragment effective
   length penalty for gDNA by 2–3 orders of magnitude and restore
   gDNA_rate on the mega-locus toward ~0.80.
3. **Decide P2 fate** after P1 is addressed. If P1 is fixed, P2 may
   become the dominant remaining asymmetry and warrant the cleanup
   for correctness, even if the magnitude is small.
