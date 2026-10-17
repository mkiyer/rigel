# The mode-free capture reader — external review packet (2026-10-13)

**Status.** This is a prototype measured on the full benchmark matrix; nothing is committed. The owner has chosen to
ship Rigel 0.8.0 as "Option A": the count foundation (commit `3e7e0fbc`) with the mode-free capture reader described
here. That integration waits on this review.

**What you are asked.**

1. Review the reader's theory, implementation and evidence, and find opportunities to improve it before it is
   integrated.
2. Assess the probe-footprint mechanism (§8), which is a candidate follow-up aimed at the reader's largest remaining
   loss.

The questions are in §10. The evidence bundle and the reproduction paths are in §11.

**After the review (2026-10-13, afternoon).** The external review (`RNA_RELEASE_0_8_0_MODEFREE_INDEPENDENT_REVIEW.md`)
found errors in this packet. Every correction below was re-checked. The sections it refers to are left as reviewed, so
that the review's citations still resolve. The new arms are pinned full-matrix runs:

| Arm | What it is | Results |
|---|---|---|
| `mf_cap` | Foundation + touched-max cap (`site_cap`) | `out/mf_cap` |
| `mf_tglbg` | tglb + cap at `S` (`site_tglbg`) | `out/mf_tglbg` |
| `mf_tglbm` | tglb + cap at `max(S, touched)` (`site_tglbm`) | `out/mf_tglbm` |

The comparison tables come from `cap_compare.py` and `arm_tally.py`.

- **§6.1 is withdrawn.** All 244 ladder "flanks" overlap no probe part. They are Poisson noise, selected on realized
  counts with a median of 3 true fragments. The mode-free reader reads them at exactly the bulk, which is right; the
  foundation reads them ×1.7.
  - By probe geometry (ladder g05 κ 0.99 ON), unprobed regions (truth: the bulk) have rms log error 0.93 under the
    mode-free reader against 2.96 under the foundation.
  - §8's premise falls with §6.1.
- **The in-scope loss was the junction cap of §4.3, not the reader.** All of it sits in junction-probed transcripts.
  - The touched-max cap alone, on the foundation's efficiencies: ladder κ 0.99 ON 3.87 → 4.22; RNA-short
    13.32 → 14.17; RNA-long 4.25 → 4.27.
  - Why it is wrong: the simulator binds a junction-spanning probe contiguously in RNA but as two parts in gDNA. A
    spliced fragment can therefore out-capture every gDNA object it touches, and the touched-max cap forbids that.
- **The corrected cap is `max(S, touched max)`**, where `S` is the most massive enriched basin's peak against the bulk
  (1 when there is none). In the foundation's unit this is exactly the old clip at 1.
  - Arm `tglbm`, κ 0.99 ON: ladder 3.87, RNA-short 13.50, RNA-long 5.17. Sparse 21.60 and junction 30.29 are kept.
  - `S` alone loses where enriched basins lie above `S`: sparse 25.85, junction 34.86.
  - Tally against the foundation: in-scope 7 better / 72 within ±0.3 / 17 worse; deferred 3 / 16 / 3.
- **The RNA-long residual (+0.92) is not object accuracy.** On well-measured captured objects the mode-free error
  against realized truth is median −0.02, sd 0.11; the foundation's is −0.19, sd 0.24. Long RNA fragments saturate
  probe overlap more uniformly than the short gDNA fragments that measure it, and the foundation's single clip happens
  to mimic that. This is a fragment-length frame problem.
- **§3.2, property 1.** "Every contaminated capture-OFF row matches" is false: lowg 6.97 → 7.51 and d10
  16.56 → 16.97. Both OFF landscapes have more than one basin.
- **§5.1.** Each gap panel holds four g50 conditions, not the ladder's grid.
- **§5.5, split by scope.** In-scope (96 strata): 10 better / 70 within ±0.3 / 16 worse. Deferred (22 strata):
  4 / 13 / 5.
- **§4.2.** The listing is a specification and is unbatched: one dense array is 4.37 GB at 2.1 M objects. The frozen
  implementation batches; VCaP peak RSS is 11.2 GB against the foundation's 10.7 GB.
- **§4.1.** The reference-reader equality covered 91,284 object evaluations across seven dumps. The 70,176 quoted
  there is one dump.
- **Appendix A.** The census drops an endpoint maximum whenever an interior peak exists. This is not live on any of
  the 120 saved landscapes.

---

## 1. Summary

Rigel corrects transcript and gDNA quantification for **hybrid capture**. Capture enriches the fragments that overlap
probes by orders of magnitude. gDNA is a uniform template before capture, so each genomic object's gDNA density
measures how strongly capture enriched it. The EM consumes those per-object **capture efficiencies** as contractions of
every component's effective length.

**The previous reader** (the foundation's) read efficiencies relative to a **located enriched mode** of calibration's
gDNA-density landscape. It returned no correction at all whenever no mode qualified. That gate fires on every
capture-OFF and zero-gDNA library, which is correct. It also fires on the single-probe panel, on the unstranded sparse
and junction-probe panels, and on **two of the three real cfRNA plasma libraries tested**. A small panel targeting scarce
transcripts against a dominant background, the owner's actual cfRNA use case, rarely forms a mode. Under the gate,
the probed target gene of the single-probe panel is over-called by +65% (κ 0.99) to +117% (κ 0.70) at 50% gDNA.

**The mode-free reader** keeps the landscape and drops the mode.

- **Prior.** The landscape's largest basin is collapsed to one level (the off-probe bulk: one uniform template). The
  rest of the landscape is the enriched spectrum, used as is.
- **Likelihood.** Each object's own Poisson count.
- **Readout.** The posterior mean of **log** density, relative to the bulk.
- **Saturation.** Each enriched basin is capped at its own peak.

There is no mode, no membership test, no width test, no gate and no new constant.

| Result | Foundation | Mode-free |
|---|---:|---:|
| Single-probe target gene, net error at g50, κ 0.99 / 0.70 | +64.7% / +116.7% | −9.6% / +0.1% |
| Sparse panel, stranded capture-ON (transcripts %) | 24.43 | 21.60 |
| Junction-probe panel, stranded capture-ON | 32.05 | 30.29 |
| Plasma objects read as enriched (cfRNA-A / cfRNA-C) | 0 / 0 | 611 / 376 |
| Strata within ±0.3 points of the foundation (all 16 panels) | — | 83 of 118 |
| Strata better by more than 0.3 / worse by more than 0.3 | — | 14 / 21 |

**Cost on the three in-scope stranded capture-ON strata of the deciding panels:**

- ladder 3.87 → 4.34;
- RNA-short gap panel 13.32 → 14.54;
- RNA-long gap panel 4.25 → 5.37.

On those strata the reader still beats the current `main` branch on two of three (6.05 / 14.24 / 6.09). The mechanism
is identified (§6.1), and §8 proposes a remedy.

---

## 2. Background the reviewer needs

### 2.1 Rigel and the rules this design must respect

Rigel is a Bayesian RNA-seq quantifier that separates RNA from gDNA contamination. Its stages:

1. A single-pass scan tallies fragments per genomic object.
2. **Calibration** deconvolves each object's unspliced mass into gDNA and RNA by strand and splicing evidence, using
   message passing along a chain of objects.
3. A per-locus EM assigns RNA to transcripts.

Calibration fits the **gDNA-density landscape**, a population prior on log gDNA density, and feeds it back into its
own solve.

The design rules that bind this work:

- **Axiom 0.** Every object holds at most three populations: gDNA, RNA on the + strand and RNA on the − strand.
  "Mature" and "nascent" are not populations.
- **Capture is a spectrum.** No detector and no library-level on/off gate. Capture is read per object, continuously and
  locally.
- **No new empirical constant.** Every number must be derived or already defined.
- **Rulings need benchmarks.** A design choice is judged on pinned, per-stratum A/B data, never on a derivation alone.

The three in-scope strata are unstranded × capture-OFF, stranded × capture-OFF and stranded × capture-ON.
**Unstranded × capture-ON is deferred**: it is reported but is not a release requirement.

### 2.2 What the reader consumes

**Objects.** The genome is partitioned into **regions**: intergenic, intron and exon pieces of the annotation. Between
consecutive regions are **boundaries**.

- A region holds the fragments **contained** in it.
- A boundary holds the fragments that **cross** it.

For each object `o`, calibration publishes a deconvolved gDNA count `k_o` (a real number ≥ 0) and its gDNA support
`S_o` (effective placements). `S_o = 0` means the object is too short to contain a fragment.

**The landscape** (`calibration/landscape.py`, unchanged) is a weighted sum of per-region Poisson kernels on a log-ρ
grid of `K = 260` points. Its construction:

- **Training set.** Regions only. Both-strand regions are excluded, and so are slots the solve cannot locate
  (Var(log f_g) > 1 nat²). Zero-count intergenic and intronic regions are included as anchors.
- **Weights.** Each kernel is weighted by a reliability factor.
- **Smoothing.** Each kernel is widened by its k-nearest-neighbour spacing, with `k = √n`.
- **Padding.** One uniform pseudo-region of Laplace smoothing is added. `logP` includes it; we call the density before
  it `logP_free`.

**The census** (`landscape._census`, `split_basins`, unchanged; Appendix A) splits the grid into basins at the interior
minima between local maxima. The **bulk** is the basin with the most rendered mass; `above` are the basins wholly above
it.

### 2.3 How efficiencies enter the EM, and why only ratios matter

**Transcripts.** A transcript's EM effective length is

    ℓ_t = fl_t · Σ_o a_to e_o / Σ_o a_to

- The sum runs over its pieces (at their region's efficiency) and its cuts.
- A contiguous cut takes its boundary's efficiency.
- A junction is priced by conservation of bases: `c_lo + c_hi − ½(c_intron,lo + c_intron,hi)`, clipped below at 0 and
  capped above.
- `a_to` are conserved placement shares.

**The gDNA component.** A locus's gDNA effective length is `Σ_r S_r e_r + Σ_b M_b e_b`.

**The EM.** The per-fragment weight of a component is `log θ − log ℓ`. Multiplying every efficiency in a locus by one
constant multiplies every length by it, so the EM's answer is unchanged; this was measured to 1e-15 on `main` when its
reader's unit changed. The gDNA prior pseudocount is a count and is not length-scaled.

**Hence a reader's unit is free.** Only these matter:

1. the ratios between objects within a locus;
2. any non-linear rule a consumer applies, such as the junction cap and the locus-length cap.

### 2.4 The foundation's reader, and what its mode does

    c_o = E[ min(ρ_o/ρ_ref, 1) | k_o, S_o ]

This is computed under the padded landscape, with a Poisson likelihood on the point count. `ρ_ref` is the peak of the
**located enriched mode**: the basin above the bulk holding the most located kernels (count ≥ 1). It qualifies only if
it has more than `√n_located` members and their median kNN width² is at most `located_var`. Otherwise
`ρ_ref = None` and every efficiency is exactly 1.

Since the unit is free, the mode does exactly two jobs.

1. **Saturation.** Every object at or above the captured level reads 1, so captured objects share one level instead of
   their Poisson scatter. This matters: removing only the clip, with everything else the foundation's (the `fu` arm),
   costs +0.65 / +0.37 points on the test chromosome's κ 0.99 / 0.70 capture-ON strata.
2. **The gate.**

Where the gate fires, measured from the foundation's published efficiencies (all exactly 1):

| Panel / library | Gate fires on |
|---|---|
| Every capture-OFF condition, every zero-gDNA (g00) condition | all — correct |
| Single-probe panel (one 125 bp probe) | all 6 conditions |
| Sparse probe panel (219 transcripts, 112 genes) | every unstranded (κ 0.50) capture-ON condition |
| Junction-spanning probe panel | every unstranded capture-ON condition |
| 1/10 depth, 0.1% gDNA, capture-ON | yes |
| cfRNA-A, cfRNA-C (real cfRNA plasma) | yes: no correction at all. The third plasma library, cfRNA-B, has a mode |

---

## 3. The mode-free reader: theory

### 3.1 The model

For object `o` with `S_o > 0`, use the foundation's own likelihood: the deconvolved count is treated as a Poisson
count.

    ℓ_o(x) = k_o (x + log S_o) − e^x S_o           x = log ρ

**The prior (two groups, Efron's construction on the fitted landscape).**

    G(x) = π_0 · δ(x − x_b)  +  Σ_{basins β ≠ bulk} p(x) · 1[x ∈ β]

| Symbol | Meaning |
|---|---|
| `x_b` | the bulk basin's peak |
| `π_0` | the bulk basin's mass under `p`, the normalized `logP_free` |
| `p` on other basins | the landscape's own density, at every level outside the bulk basin (above or below it) |

**Saturation.** For each grid level `x` inside an enriched basin β (a basin wholly above the bulk), `s(x) = peak(β)`.
Elsewhere `s(x) = +∞`.

**The efficiency.**

    log e_o = E_post[ min(x, s(x)) ] − x_b,      post ∝ G(x) · exp ℓ_o(x)
    e_o = 1                                       when S_o = 0

### 3.2 Properties, each with its reason

1. **A one-population landscape reads exactly 1 everywhere.** If the census finds one basin, `G` is a point mass at
   `x_b`, so every posterior is that point. Capture-OFF libraries come out uniform with no decision taken; the
   foundation's gate becomes a consequence of the formula. Measured: every contaminated capture-OFF row matches the
   foundation; the test chromosome's κ 0.99 g05 capture-OFF row reads exactly 1 on every region and boundary; and on
   its zero-gDNA capture-OFF row every genic region reads within 1.3% of 1. On that row the only exception is the
   unannotated contig itself (next property).

2. **A level supported by one object cannot lift other objects.** Take an object whose likelihood is flat over the
   support of `G`. Its log efficiency is the prior's mean saturated log level above the bulk:

       log e_o = Σ_β μ_β · Δ_β,   Δ_β = E[min(x,s) | β] − x_b

   A lone basin has `μ_β ≈ 1/n`, so the lift is `≈ e^{Δ_β/n}`.

   The linear posterior mean instead gives `1 + μ_β (e^{Δ_β} − 1)`, which scales with the level **ratio**. The case
   that decided it is the test chromosome's 1.5 Mb intergenic contig. At zero gDNA it carries 8,938 fragments of
   **unannotated mRNA**, which calibration must call gDNA under Axiom 0. That forms a basin of mass 0.00215 at
   10^−2.22 per bp, while the bulk sits at the grid floor, 10^−6.15.

   | Readout | Lift for short objects | Test κ 0.99 g00 capture-OFF |
   |---|---:|---:|
   | Linear posterior mean | about 4–20× | 6.99 → 32.61 |
   | Log posterior mean | e^{0.00215 · 9.05} = 1.02× | 6.95 |

3. **An informative object reads its own level wherever the landscape holds mass there.** As `k_o S_o` grows, the
   posterior concentrates at `log(k_o/S_o)` and `e_o → min(k_o/S_o, s)/ρ_b`. One probed exon's own kernel forms a small
   basin (single-probe panel: mass 0.0008 at 10^1.35, against the bulk at 10^−1.22), and the exon reads 340–490×.

4. **Saturation is per level.** Upward scatter within a level reads as the level, and distinct levels stay distinct.
   This matters because the spectrum is real.

   - On the sparse panel at g50 κ 0.99, the truth's captured objects sit at two levels: 10^0.26 and 10^0.6–0.75.
   - On the junction-probe panel, captured levels span 10^−0.66 to 10^0.63, because junction-spanning probes capture
     gDNA unevenly.

   The foundation's single clip at its mode flattens both.

5. **Unit invariance.** `e_o` is relative to `ρ_b = e^{x_b}`, and the consumers were made unit-free (§4.3).

### 3.3 Why each ingredient: the measured failures that forced it

| Ingredient | Without it | Measured |
|---|---|---|
| Bulk collapsed to one level | The smoothed bulk carries the training objects' own Poisson spread, so as a prior it under-shrinks. Unprobed objects read noise, which flips allocations in near-identical isoform clusters. | Single-probe g05 κ 0.99: 10.53 without the collapse, 7.70 with it, other ingredients held fixed (foundation 7.85) |
| Log readout | A lone basin above a floor bulk lifts every short object. | Test κ 0.99/0.70/0.50 g00 OFF: 6.99/6.35/5.67 → 32.61/31.85/32.72 |
| Saturation | Captured objects keep their scatter. | Test κ 0.99 ON: 8.33 without, 7.44 per level |
| Per level rather than one level | A single level flattens real spectra. | Junction κ 0.99 ON: 31.98 with one level, 30.29 per level, 24.36 with none |
| Point count rather than the solve's Var(log f_g) | Widening by the solve's variance erases capture at scarce gDNA, where the point estimate still carries the signal. | Widening against point counts, other ingredients held fixed: lowg g001 ON 38.30 vs 22.04; d10 g01 ON 40.38 vs 29.81 (foundation 19.22, 26.91) |

### 3.4 What the reader deliberately does not use

- no mode, membership count, width test or gate;
- no constant: the census and the landscape are the foundation's own;
- not the solve's Var(log f_g);
- not region types or locality: every object is read alone against one population prior.

---

## 4. Implementation, exactly

### 4.1 Where it is, and the identity check

- **Prototype worktree:** `<worktree>`, the foundation `3e7e0fbc` plus uncommitted edits behind an
  environment switch, `RIGEL_MODEFREE`. The shipped arm is `RIGEL_MODEFREE=tglb`.
- **Frozen build that produced every `tglb` number:**
  `$RIGEL_RUNS/prototypes/2026-10-13_modefree/site_tgl3`.
- **Its exact diff against the foundation:** `review/tglb_exact.diff` (425 lines, 6 files; it also carries other
  arms, selected by the switch).
- **Identity:** with the switch off, the build reproduced the foundation byte for byte on 108 of 108 result rows
  (single-probe and test panels, every axis).

### 4.2 The reader: standalone listing, verified equal to what ran

The listing below is `reference_reader.py`. It reproduces the prototype's `tglb` efficiencies to within 4.5e-14 in log
on seven rows: single, test (capture ON, OFF and zero-gDNA), sparse and ladder (70,176 objects). The prototype's
efficiencies in turn equal the end-to-end runs' published efficiencies exactly. The listing computes the bulk level's
likelihood in closed form; the prototype interpolates it from the grid, at a grid point, which is exact.

```python
def mode_free_efficiencies(log_rho, logP_free, bulk, above, k, S, rate_floor):
    x = np.asarray(log_rho, dtype=np.float64)
    k = np.maximum(np.asarray(k, dtype=np.float64), 0.0)
    S = np.maximum(np.asarray(S, dtype=np.float64), 0.0)
    tiny = np.finfo(float).tiny

    # THE TWO-GROUP PRIOR. The null: the bulk basin collapsed to its peak x_b, carrying the basin's whole mass pi0.
    # The alternative: the landscape's own (unpadded) density everywhere outside the bulk basin — any level, any mass.
    in_bulk = (x >= bulk.lo) & (x <= bulk.hi)
    p = np.exp(np.asarray(logP_free, dtype=np.float64) - np.max(logP_free))
    p /= p.sum()
    pi0 = float(p[in_bulk].sum())
    slab = np.where(in_bulk, 0.0, p)
    log_slab = np.where(slab > 0.0, np.log(np.maximum(slab, tiny)), -np.inf)
    xb = float(bulk.log_rho)

    # PER-LEVEL SATURATION: a grid level inside an enriched basin reads at most that basin's peak.
    cap = np.full(x.shape, np.inf)
    for m in above:
        cap[(x >= m.lo) & (x <= m.hi)] = float(m.log_rho)

    # THE OBJECT'S OWN POISSON LIKELIHOOD of its point count on its support, on the grid and at the bulk level.
    lam = np.exp(x)[None, :] * S[:, None]
    ll = np.where(S[:, None] > 0.0, k[:, None] * np.log(np.maximum(lam, rate_floor)) - lam, 0.0)
    lam_b = np.exp(xb) * S
    ll_b = np.where(S > 0.0, k * np.log(np.maximum(lam_b, rate_floor)) - lam_b, 0.0)

    # THE POSTERIOR over {the bulk level} ∪ {the grid outside the bulk basin}, in log space.
    a_null = (np.log(pi0) if pi0 > 0.0 else -np.inf) + ll_b
    a_alt = log_slab[None, :] + ll
    top = np.maximum(a_null, a_alt.max(1))
    top = np.where(np.isfinite(top), top, 0.0)
    w_null = np.exp(a_null - top)
    w_alt = np.exp(a_alt - top[:, None])
    den = np.maximum(w_null + w_alt.sum(1), tiny)

    # THE READOUT: the posterior mean of LOG (saturated) density, against the bulk level.
    mean_log = (w_null * xb + w_alt @ np.minimum(x, cap)) / den
    return np.where(S > 0.0, np.exp(mean_log - xb), 1.0)
```

**How it is called.**

- `bulk, above = split_basins(_census(landscape))` reads the census on the **padded** `logP`, as the foundation does.
- `logP_free` is the landscape's density before the pseudo-region, which `fit_landscape` now keeps (§4.3).
- `rate_floor` is `CONSTANTS.calibration.rate_log_floor`, the existing floor the foundation's reader uses.
- The reader is called once on the regions (contained counts and supports) and once on the boundaries (crossing counts
  and supports). It runs after calibration's solve and feeds nothing back into it.

### 4.3 Every other change in the shipped arm

All of these are in `review/tglb_exact.diff`; in the prototype each is gated on the switch.

| File | Change |
|---|---|
| `landscape.py` | `DensityLandscape.logP_free`: the rendered density before the uniform pseudo-region, stored beside `logP`. Calibration's own use of the landscape is unchanged. |
| `calibrate.py` | When a landscape exists, the efficiencies come from the reader above. `gdna_reference_density` publishes `ρ_b` and `gdna_reference_members` publishes the landscape's training count, both as metadata. With no landscape, efficiencies are 1, as before. |
| `capture_eff_length.py` | The junction cap changes from `1` to the **most captured touched object**: `max(c_region[A], c_region[B], c_boundary[lo], c_boundary[hi])`. This is the unit-free form of "a fragment across a junction is captured no more than a fully captured piece". It is `main`'s rule, measured there on the junction-probe layout. |
| `priors.py` | The locus gDNA length drops its `min(·, span)`. That was an ulp absorber while every efficiency was ≤ 1; with efficiencies > 1 it would bind and break unit invariance. |
| `result.py` | The `[0, 1]` validation of efficiencies becomes `≥ 0`, finite. |

**Units.** The bulk reads 1 and enriched objects read above 1. The reported `em_effective_length` column therefore
exceeds the fl-marginal length for captured transcripts. TPM uses the plain length by an earlier ruling and is
unaffected.

### 4.4 Cost

Measured on one thread, on 2.1 million objects (VCaP scale, resampled from the ladder's objects) with `K = 260`:

| Reader | Time |
|---|---:|
| Foundation's clipped posterior mean | 3.9 s |
| Prototype `tglb` path (per-row `np.interp` loop) | 5.7 s |
| Standalone listing (closed-form spike) | 4.3 s |

The listing equals the prototype to within 5e-15. No native port is needed. Real-library wall times are in §5.4 for
completeness; those runs were sequential under varying load, so they are indicative only.

---

## 5. Evaluation

### 5.1 Protocol

- **Tool runs.** End to end through `scripts/design/quant_accuracy.py --arm base`, pinned:
  `--set scan.total_threads=1 --set em.assignment_mode=fractional`, `OMP_NUM_THREADS=1`, EM seed 0.
- **Builds.** Every arm ran from its own frozen site, never the editable install.
- **Metric.** Transcripts: `Σ|est − truth| / Σ truth`. Genes: the same after summing by gene. Truth is the simulator's
  post-capture mRNA fragment count per transcript. Strata are ratios of sums over their conditions, with g00 shown
  apart. In the tables below a cell reads `transcripts % (genes %)`.
- **Deciding panels.**
  - **Ladder:** 16 conditions on a reference carved from the human genome: g00, g05, g50 and g98 gDNA × κ 0.99 and
    0.50 × capture ON/OFF. RNA and gDNA have equal fragment lengths by design.
  - **Gap panels:** RNA fragments shorter (RNA-short) or longer (RNA-long) than gDNA, otherwise the ladder's
    conditions.
- **Robustness panels**, on the 7.93 Mb test chromosome plus a 1.5 Mb unannotated contig:

  | Panel | Conditions |
  |---|---|
  | `test` | 5 gDNA levels × κ 0.50 / 0.70 / 0.99 × capture ON/OFF |
  | `single` | one 125 bp probe |
  | `sparse` | probes on 219 transcripts of 112 genes |
  | `junction` | junction-spanning probes |
  | `bind01`, `bind1`, `bind100` | binding strengths |
  | `odg05` | gDNA strand overdispersion 0.05 |
  | `fle200`, `flgl`, `flrl` | fragment-length variants |
  | `lowg` | 0.1% and 1% gDNA |
  | `d10` | 1/10 depth |

- **Real libraries.** VCaP: an exome-captured DNA plus transcriptome mix, with per-fragment origin truth by flowcell;
  the gDNA fraction truth ≈ 0.2518 is the primary-pair count, not the pipeline-filtered domain. VCaP
  transcriptome-only subset: a real zero-gDNA library. Three hybrid-captured plasma cfRNA libraries without truth
  (cfRNA-A, cfRNA-C, cfRNA-B).

### 5.2 Deciding panels

| Ladder stratum | Foundation | Mode-free | Main |
|---|---:|---:|---:|
| κ 0.99 OFF contaminated | 2.01 (0.26) | 2.02 (0.26) | 3.37 (0.33) |
| κ 0.99 OFF g00 | 1.93 (0.14) | 1.93 (0.14) | 13.61 (0.83) |
| κ 0.99 ON contaminated | 3.87 (1.38) | 4.34 (1.49) | 6.05 (1.37) |
| κ 0.99 ON g00 | 6.47 (2.36) | 6.47 (2.36) | 15.10 (2.61) |
| κ 0.50 OFF contaminated | 2.30 (0.31) | 2.36 (0.31) | 2.40 (0.32) |
| κ 0.50 OFF g00 | 1.59 (0.14) | 1.59 (0.14) | 1.59 (0.14) |
| κ 0.50 ON contaminated (deferred) | 7.63 (3.76) | 7.74 (3.86) | 62.96 (23.69) |
| κ 0.50 ON g00 (deferred) | 7.27 (3.22) | 7.27 (3.22) | 7.27 (3.22) |

| Gap-panel stratum | Foundation | Mode-free | Main |
|---|---:|---:|---:|
| RNA-short κ 0.99 OFF | 3.54 (0.23) | 3.55 (0.23) | 4.99 (0.36) |
| RNA-short κ 0.99 ON | 13.32 (1.42) | 14.54 (1.47) | 14.24 (1.50) |
| RNA-short κ 0.50 OFF | 3.80 (0.27) | 3.82 (0.27) | 3.95 (0.28) |
| RNA-short κ 0.50 ON (deferred) | 18.47 (2.32) | 20.10 (2.32) | 68.06 (11.46) |
| RNA-long κ 0.99 OFF | 2.26 (0.19) | 2.26 (0.19) | 3.52 (0.25) |
| RNA-long κ 0.99 ON | 4.25 (0.68) | 5.37 (0.64) | 6.09 (0.65) |
| RNA-long κ 0.50 OFF | 2.16 (0.23) | 2.17 (0.23) | 2.55 (0.23) |
| RNA-long κ 0.50 ON (deferred) | 43.22 (11.26) | 43.34 (11.53) | 34.11 (13.40) |

### 5.3 Where the foundation is blind

| Single-probe target gene gB1_cap (net error) | Foundation | Mode-free | Main |
|---|---:|---:|---:|
| g50 κ 0.99 | +64.68% | −9.55% | −10.99% |
| g50 κ 0.70 | +116.66% | +0.07% | −2.94% |
| g50 κ 0.50 (deferred) | +102.96% | −59.33% | −59.43% |
| g05 κ 0.99 | +5.74% | +0.43% | +0.41% |
| g05 κ 0.70 | +6.52% | −1.55% | −2.88% |
| g05 κ 0.50 (deferred) | +6.73% | +6.73% | +6.73% |

| Panel stratum, transcripts (genes) | Foundation | Mode-free | Main |
|---|---:|---:|---:|
| single κ 0.99 ON | 8.82 (1.52) | 8.23 (1.05) | 11.46 (1.08) |
| single κ 0.70 ON | 8.44 (2.04) | 7.57 (1.12) | 8.28 (1.20) |
| sparse κ 0.99 ON | 24.43 (3.10) | 21.60 (2.66) | 21.20 (2.66) |
| sparse κ 0.70 ON | 26.08 (4.39) | 25.86 (3.80) | 30.88 (4.76) |
| sparse κ 0.50 ON (deferred) | 94.42 (66.47) | 80.14 (52.61) | 86.56 (58.40) |
| junction κ 0.99 ON | 32.05 (1.12) | 30.29 (1.18) | 23.39 (1.49) |
| junction κ 0.70 ON | 32.37 (1.59) | 31.62 (1.72) | 21.22 (2.15) |

**The sparse panel's 112 probed genes**, transcripts / genes %:

| Condition | Foundation | Mode-free |
|---|---:|---:|
| g05 κ 0.99 | 21.47 / 1.63 | 18.52 / 1.10 |
| g25 κ 0.99 | 24.73 / 2.92 | 22.64 / 2.36 |
| g25 κ 0.70 | 23.02 / 4.25 | 21.46 / 2.94 |
| g50 κ 0.99 | 28.08 / 4.19 | 24.09 / 3.70 |
| g50 κ 0.70 | 27.55 / 5.73 | 27.16 / 4.89 |
| g98 κ 0.99 | 60.49 / 27.71 | 58.31 / 34.24 |

Unprobed genes worsen at high gDNA:

| Condition | Foundation | Mode-free |
|---|---:|---:|
| g98 κ 0.99 | 157.5 | 209.7 |
| g98 κ 0.70 | 161.9 | 251.6 |
| g25 κ 0.70 | 30.5 | 36.3 |

### 5.4 Real libraries

| Library | Reader | gDNA fraction | mRNA | nRNA | Objects read > e× bulk | Wall s | Peak RSS GB |
|---|---|---:|---:|---:|---:|---:|---:|
| cfRNA-A | foundation | 0.0850 | 126,512 | 7,699 | 0 | 10 | 3.9 |
| cfRNA-A | mode-free | 0.0889 | 125,048 | 8,590 | 611 | 18 | 6.2 |
| cfRNA-C | foundation | 0.1587 | 632,561 | 69,895 | 0 | 23 | 4.2 |
| cfRNA-C | mode-free | 0.1615 | 628,069 | 72,030 | 376 | 32 | 6.3 |
| cfRNA-B | foundation | 0.9265 | 46,596 | 48,333 | 0 | 32 | 6.4 |
| cfRNA-B | mode-free | 0.9259 | 47,942 | 47,717 | 50,243 | 58 | 6.2 |
| VCaP mix | foundation | 0.2395 | 13,336,410 | 755,354 | 0 | 206 | 10.7 |
| VCaP mix | mode-free | 0.2391 | 13,356,328 | 742,697 | 603,664 | 190 | 11.2 |
| VCaP transcriptome-only | foundation | 0.0080 | 13,419,032 | 343,733 | 0 | 130 | 8.8 |
| VCaP transcriptome-only | mode-free | 0.0077 | 13,410,440 | 355,694 | 1,824 | 160 | 8.8 |

The foundation's efficiencies are ≤ 1 by construction, so the "> e×" column is structurally 0 for it. VCaP DNA
recovery (mix minus transcriptome-only, against 4,708,907 primary-pair DNA fragments) is 91.9% under the foundation
and 91.8% under the mode-free reader. Neither tree has `main`'s multimapper second pass, which is what moves VCaP to
0.2528 on `main`; it is a separate component of Option A.

The plasma libraries have no truth. The reader reads hundreds of enriched objects where the foundation read none, and
moves their gDNA fraction by +4.6% and +1.8% relative.

### 5.5 Full tally against the foundation (118 strata, 16 panels)

**Better by more than 0.3 points (14):**

- sparse deferred −14.28; sparse κ 0.99 ON −2.83;
- junction deferred −1.90; junction κ 0.99 ON −1.75; junction κ 0.70 ON −0.75;
- fle200 κ 0.50 OFF g00 −1.19; fle200 κ 0.99 OFF g00 −0.41; fle200 κ 0.70 ON −0.38;
- bind01 deferred −0.91; bind1 κ 0.70 ON −0.48; bind100 κ 0.70 ON −0.37;
- single κ 0.70 ON −0.87; single κ 0.99 ON −0.58; single deferred −0.44.

**Worse by more than 0.3 points (21):**

| Stratum | Δ |
|---|---:|
| odg05 κ 0.70 ON | +3.11 |
| odg05 κ 0.70 OFF g00 | +2.82 |
| flgl κ 0.70 OFF g00 | +2.75 |
| flrl κ 0.99 OFF g00 | +2.44 |
| flrl κ 0.70 OFF g00 | +2.32 |
| d10 κ 0.99 ON | +1.68 |
| RNA-short deferred | +1.63 |
| lowg κ 0.99 ON | +1.55 |
| RNA-short κ 0.99 ON | +1.23 |
| RNA-long κ 0.99 ON | +1.12 |
| flrl κ 0.70 ON | +0.85 |
| test deferred | +0.61 |
| fle200 deferred | +0.59 |
| lowg κ 0.99 OFF | +0.55 |
| d10 κ 0.99 OFF g00 | +0.51 |
| ladder κ 0.99 ON | +0.47 |
| flgl κ 0.99 OFF g00 | +0.46 |
| d10 κ 0.99 OFF | +0.41 |
| odg05 deferred | +0.41 |
| odg05 κ 0.99 ON | +0.37 |
| flgl deferred | +0.31 |

Every stratum of every panel is in `review/tables_strata.md`; every condition is in `review/tables_per_condition.md`.

---

## 6. Where it loses, and why

### 6.1 Flank flattening: the in-scope stranded capture-ON cost

The bulk basin, collapsed to one level, also contains the **partially captured flanks**: introns and pieces within
about a fragment length of a probe, at 2–30× the bulk. They lie in the basin's upper shoulder, below the valley. With
no alternative mass at their level, a flank's posterior chooses between the bulk and the captured basin, and reads
near the bulk.

The table scores regions against the certified per-object truth (realized gDNA fragments over the object's support),
on the bulk scale. Cells are log error, median / rms. Both readers ran on the same calibration; their dumps have
identical counts.

| Truth class | Ladder g05 κ 0.99 ON: n | Foundation | Mode-free | Test g50 κ 0.99 ON: n | Foundation | Mode-free |
|---|---:|---:|---:|---:|---:|---:|
| bulk (< 2×) | 1017 | −0.03 / 0.21 | −0.11 / 0.35 | 639 | −0.01 / 0.17 | +0.01 / 0.31 |
| flank (2–30×) | 244 | **−0.46 / 0.71** | **−1.00 / 1.20** | 71 | −0.57 / 0.70 | −0.74 / 0.76 |
| captured (> 300×) | 3142 | −0.08 / 0.33 | −0.08 / 0.37 | 363 | −0.01 / 0.07 | −0.01 / 0.07 |

**Caveat.** The truth is *realized* counts. A reader that puts every bulk object exactly at the bulk level is scored
for the truth's own Poisson scatter, so the bulk-class rms overstates its error. The flank class's median shift is
the systematic effect.

**Attribution.** The closest linear variant measured (`tgtla`, §7: type-weighted, with a collapsed null) reaches 4.25
on the ladder against 4.34. So the log readout accounts for at most about 0.1 of the 0.47, and the collapsed bulk for
the rest.

### 6.2 False gDNA on micro-boundaries at zero gDNA

Four zero-gDNA capture-OFF strata lose +2.3 to +2.8: flrl κ 0.99 and 0.70, flgl κ 0.70, odg05 κ 0.70. In those
simulations calibration deposits false gDNA on boundaries between regions only 1–10 bp long:

- on the κ 0.70 rows, 42–145 fragments per boundary, a median gDNA share of 5–14%;
- on the κ 0.99 row, about one fragment;
- the solve's own Var(log f_g) for them is 19–32 nat².

The reader takes the point count as precise and reads those boundaries at 10^3× to 10^5×. The foundation's gate hid
them.

Widening by the solve's variance (`tglh`) fixes three of these rows: odg05 7.82, flgl 6.37 and flrl κ 0.99 5.00,
against the foundation's 8.01, 6.40 and 4.93. It reduces the fourth, flrl κ 0.70, to 8.21 against 6.64. It is still
refused (§3.3), because it erases capture at scarce gDNA. The defect is calibration's; the reader exposes it.

### 6.3 Strand noise splits the captured basin

On odg05 κ 0.70 ON (+3.11), noise in the deconvolved gDNA splits the captured level into several basins, and
per-level saturation keeps them as distinct levels. Saturating at one level (`tglc`) scores 13.31 against the
foundation's 12.87, but loses the sparse and junction gains.

### 6.4 Information-starved captured libraries

At 1/10 depth with 0.1% or 1% gDNA, the whole library holds about 117 or 1,170 gDNA fragments. Per-object reads are
then Poisson 0/1 flips: the d10 κ 0.99 ON stratum loses +1.68 and lowg κ 0.99 ON loses +1.55. At these rows the
foundation's refusal to read is genuinely the better answer.

---

## 7. Variants measured and refused

Every number is recomputed from the archived JSONL (`review/refusals.txt`).

| Variant | What failed | Measured |
|---|---|---|
| Linear posterior mean, no gate (`tgn`) | A lone basin lifts short objects at zero gDNA | test κ 0.99 g00 OFF 6.99 → 32.61 |
| Null ratio: own evidence ÷ the same measurement at the bulk, with Var widening (`tg`) | Under-reads short captured pieces | test κ 0.99 / 0.70 ON 7.32 / 7.76 → 10.79 / 13.06 |
| Log readout with Var(log f_g) widening (`tglh`) | Erases the capture signal at scarce gDNA | lowg g001 ON 19.22 → 38.30; d10 g01 ON 26.91 → 40.38 |
| Log readout, no saturation (`tgl`) | Wins junction, loses dense panels | junction κ 0.99 ON 24.36; test κ 0.99 ON 8.33 |
| Log readout, one saturation level (`tglc`) | Flattens real spectra | junction κ 0.99 ON 31.98; ladder ON 4.34 |
| Pooled bulk measurement + deepest-valley null (`tgp`) | No change on the deciding panels | ladder ON 4.35; test κ 0.99 g00 OFF 7.85 |
| Type-specific base rates, linear (`tgtla`) | Closes only 0.09 of the ladder gap; breaks the starved extreme | ladder ON 4.25; d10 g001 ON 14.75 → 27.28 |
| Type-specific maximum-likelihood weights, linear (`tgtm`) | Breaks the starved extreme | d10 g001 ON 39.36 |
| Deconvolved landscape (NPMLE) prior, offline | Less uniform off capture; worse on dense rows | capture-OFF max 0.79 nats; dense rms 0.52 vs 0.29 |
| Point null nested inside the bulk basin, ML split | The ML weight on the spread is ≈ 0 | ladder 0.000; test 0.012 |
| Foundation without its clip (`fu`, diagnostic) | Shows the saturation job is real | test κ 0.99 ON 7.32 → 7.97 |

---

## 8. Option (c): the probe-footprint mechanism (investigational)

**Hypothesis.** The in-scope loss (§6.1) is flank flattening. A flank is partially captured *because* it lies within
about one fragment length of a probed object. That is the probe physics: a probe enriches roughly one fragment length
around itself, and the footprint depends on probe length, both fragment-length laws, and whether the probe spans a
junction. So the information that separates a real flank from bulk noise is **local**. The population prior cannot
provide it, which is why the nested-null and NPMLE attempts put no weight there.

**Sketch.** Two passes, no new constant.

1. **First pass:** the reader above. For every object it gives `P_o`, the posterior probability of an enriched level
   (`1 − w_null` in the listing).
2. **Neighbours.** For each object, its neighbours `N(o)` are the objects within one fragment length on the region
   chain. The index already carries per-edge fragment-length reaches.
3. **Second pass**, for objects with an enriched neighbour. Let `q_o = max_{n ∈ N(o)} P_n`, or an fl-weighted
   combination. The prior opens the bulk basin's **shoulder**, the landscape's density between the bulk peak and the
   valley, with mass moved from the spike in proportion to `q_o`. An object beside a probed object can then read a
   flank level; an isolated object cannot, so capture-OFF libraries, small panels and the lone-level property are
   untouched.

**Expected ceiling.** At most the in-scope gap: 0.47 (ladder), 1.23 (RNA-short), 1.12 (RNA-long). A linear readout
with the same collapsed bulk recovers only 0.09 on the ladder, so the shoulder is the lever.

**Risks.**

- Local coupling has failed here before. Earlier attempts coupled **every** object's level (a nearest-neighbour Laplace
  prior over published posteriors; a local field on own likelihoods). This sketch only re-opens the shoulder beside
  enriched objects.
- It interacts with junction pricing, which already reads neighbours.
- It must be derived, gated and A/B'd **alone**.

**Proposed test.**

1. Offline on the saved calibration dumps: the flank class's median and rms (§6.1) on the ladder and test rows, the
   capture-OFF uniformity, and the single-probe reading.
2. Then the pinned A/B on the full matrix.

---

## 9. Integration plan for option (b), for the reviewer to stress-test

This is the owner's chosen path after the review. Each step is one mechanism at a time, with a pinned A/B or
identity check.

1. **Reader.** `capture_efficiency.py` becomes the listing in §4.2, with the bulk likelihood in closed form. The clipped
   posterior mean and `landscape.located_enriched_mode` / `LocatedMode` lose their only consumer and are deleted.
   `_census` and `split_basins` stay.
2. **Landscape.** Keep `logP_free` on `DensityLandscape`. The prototype also added `train_count`, `train_eff` and
   `train_type` for refused arms; those are dropped.
3. **Result and surfaces.** `gdna_reference_density` and `gdna_reference_members` change meaning; they are to be
   renamed or replaced, for example by the bulk density plus the number of objects read above the bulk. The affected
   surfaces:
   - `CalibrationResult` validation;
   - the `calibrate` log line;
   - the CLI summary (`cli.py`) and the HTML report's capture tile (`report/model.py`);
   - two design instruments that read the fields (`ruler_vs_truth.py`, `calibration_vs_oracle.py`).
4. **Consumers.** The junction cap becomes the touched-object maximum; the locus gDNA length loses `min(·, span)`. No
   other consumer assumes efficiencies ≤ 1; the tests named below are the check.
5. **Tests and goldens.**
   - Rewrite `tests/calibration/test_capture_efficiency.py` around the five properties of §3.2. Each gate must be shown
     to fire by breaking the code.
   - Update `test_capture_eff_length`, `test_priors`, `test_result_schema`, `test_landscape`,
     `test_landscape_training_population`, `test_prior_units`, `test_oracle`, `test_report` and
     `test_summary_report`.
   - Regenerate the goldens only after reading the diff against the truth-scored instruments.
6. **Option A's other components, separately.** Port the multimapper second pass and the ruler instrument fix from
   `main`, each A/B'd alone.
7. **Release gate.** Full matrix and real libraries on the integrated branch; suite green; `preflight.py --full`.
   Documentation: DESIGN (the reader's ruling), EQUATIONS (the two-group derivation), ISSUES (§6's losses as named
   entries), MANUAL (the efficiency unit and the summary fields).

**Open choices the review may settle.**

- Bulk level at the basin peak versus the pooled measurement: no deciding-panel difference was measured.
- Per-level versus single-level saturation.
- The log readout versus a linear readout made safe some other way.
- Whether to keep the reader in Python: 4.3 s at genome scale.

---

## 10. Questions for the reviewer

1. **The null.** Is collapsing the census's largest basin to its peak, with the basin's whole mass, the right
   two-group null for a population that physically has one off-probe level plus partial-capture flanks? Is there a
   principled null that keeps the flanks (§6.1) without reintroducing the bulk noise that hurt small panels (§3.3)?
2. **The readout.** The log readout is the Bayes estimate under squared loss in log efficiency. It is chosen for
   robustness: a level backed by one object lifts others only as `e^{Δ/n}`. The consumer, however, is linear in
   efficiency. Is there a better-founded readout that keeps that robustness and removes the Jensen under-read of
   mixed posteriors (about 0.1 point on the ladder)?
3. **Saturation.** Per-level saturation depends on the census's basins, and therefore on the smoothed landscape's
   resolution. Noise can split a level (§6.3) and nearby real levels can merge (junction). Is there a partition-free
   formulation of "capture saturates"?
4. **The solve's uncertainty.** Point counts let calibration's unresolved false gDNA through (§6.2). Symmetric
   log-normal widening by Var(log f_g) erased real scarce-gDNA capture (§3.3); one reason is that for a near-zero
   share the variance describes a one-sided lower tail. Is there a principled way to use the solve's uncertainty that
   discards §6.2's objects and keeps scarce-gDNA capture?
5. **The remaining losses.** Which of §6.1–§6.4 are structural for a per-object reader, and which are fixable within
   it?
6. **Option (c).** Is the footprint mechanism sound? Would you build it before 0.8.0, after, or never? If built, how
   should `q_o` and the shoulder mass be derived?
7. **Integration.** Do you see risks in §9: semantics of `em_effective_length`, consumers assuming efficiencies ≤ 1,
   or edge cases such as no landscape, empty basins, below-bulk basins, `S = 0` and one-basin libraries?
8. **The evidence.** Could anything in the protocol mislead: ratio-of-sums strata, realized-count truth in §6.1, the
   simulator's capture model, or differences in the comparison arms' code paths? Do the tables support the claims in
   §1?

---

## 11. Evidence and reproduction

**Portable bundle:** `.cache/rigel_runs/2026-10-13_modefree_review/rigel_modefree_review.zip`. `HANDOFF.md` inside
lists the contents and `HANDOFF_SHA256.json` the checksums. It contains:

- the foundation's source at `3e7e0fbc` (`src/`, `tests/calibration/`, `scripts/design/`);
- the shipped arm's exact diff and modified files;
- `reference_reader.py`;
- every analysis script;
- every arm's per-condition JSONL and receipts, and the comparison arms' JSONL;
- the real-library JSON files;
- the generated tables;
- the current permanent docs, `CLAUDE.md`, and this packet with its prompt.

It contains no BAMs, indexes, calibration dumps or binaries. Several scripts carry absolute paths from the machine
they ran on, so the bundle is not a portable end-to-end replay.

**On the machine they ran on:**

| Purpose | Path or command |
|---|---|
| Prototype worktree | `<worktree>` |
| Run directory | `$RIGEL_RUNS/prototypes/2026-10-13_modefree/` |
| One arm on one panel | `run_arm.sh`, with `MODEFREE=tglb SITE_DIR=$P/site_tgl3` |
| Every table | `review_tables.py`, `review_tally.py`, `review_refusals.py` |
| §6.1 | `review_flank.py` |
| Real libraries | `real_libs.sh` |
| The listing's equality check | `reference_reader.py` against `dumps_tglb/` |

---

## Appendix A: the census (foundation code, unchanged)

```python
def _census(landscape):                       # every interior local maximum with its basin
    x = landscape.log_rho
    p = np.exp(landscape.logP); p = p / max(p.sum(), _EPS)
    interior = np.where((p[1:-1] > p[:-2]) & (p[1:-1] >= p[2:]))[0] + 1
    if interior.size == 0:
        interior = np.array([int(np.argmax(p))])
    peaks = np.sort(interior)
    cuts = [0] + [int(a + np.argmin(p[a:b + 1])) for a, b in zip(peaks[:-1], peaks[1:])] + [p.size - 1]
    # basin i spans x[cuts[i]] .. x[cuts[i+1]], its mass the density inside, its peak x[peaks[i]]

def split_basins(modes):                      # the bulk = most rendered mass; `above` = basins wholly above it
    depleted = max(modes, key=lambda m: m.basin_mass)
    above = tuple(m for m in modes if m.lo >= depleted.hi - _EPS and m is not depleted)
    return depleted, above
```
