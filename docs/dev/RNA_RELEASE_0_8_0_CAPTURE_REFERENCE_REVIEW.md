# What capture efficiency is measured against — review packet (2026-10-13)

**Status.** Rigel 0.8.0 is being integrated on a release branch with the **mode-free capture reader**. It was reviewed
externally earlier today, then modified on that review's evidence: the junction cap is now
`max(saturation, best touched object)`.

The integrated build is verified:

- **The port is faithful.** It reproduces the measured prototype bit for bit across the full benchmark matrix.
- **One stress scenario in the test suite regresses badly.** The owner paused the decision there, to have the
  question underneath it reviewed: **what reference level should every capture efficiency be read against?**

Nothing in this packet is committed.

**What you are asked.**

1. Review the two references the project has used:
   - the **captured level**, the old reader;
   - the **uncaptured level, "the bulk"**, the mode-free reader.

   Cover their statistical properties and the failure modes measured here.
2. Judge the candidate change in §6, the **one-sided null**, and the alternatives in §7.
3. Recommend a path for 0.8.0, and what to investigate after it.

The questions are in §8. Figures are in `capture_reference_figures/`; evidence and reproduction are in §9.

---

## 1. Summary

**The problem in one paragraph.** Hybrid capture multiplies the sampling rate of fragments that overlap probes, by up to
~1,000× here. gDNA is a uniform template before capture, so each genomic object's gDNA density after capture measures
how much capture enriched it. Rigel turns those densities into per-object **capture efficiencies** and scales every
EM component's effective length by them. The EM is invariant to a factor common to every component in a locus. So only
**ratios** of efficiencies matter, and above all the ratio between captured and uncaptured objects. To read a ratio
you must place both ends, and whichever end you call "1" is the reference.

**Two references, two weaknesses.** The owner's framing, which the measurements below support:

| | Old reader: reference = the **captured** level | Mode-free reader: reference = the **uncaptured** level |
|---|---|---|
| Efficiency of a captured object | ≈ 1 | ≈ the enrichment, e.g. 850 |
| Efficiency of an uncaptured object | ≈ 1/850 | 1 |
| Must estimate | the captured level (a "mode" to locate) | the uncaptured level ("the bulk") |
| Breaks when | the captured population is too small to locate (small panels); a zero-fragment object reads 0, giving 0/0, so shrinkage is needed | the uncaptured gDNA is too sparse to place (very efficient capture, little gDNA); or another level outweighs it |

**What was measured today.**

- **Sparse uncaptured gDNA.** Where the whole library holds at least ~1,000 uncaptured gDNA fragments, the bulk is
  within 0.95–1.11× of the truth (92 of 92 calibrations). Below ~300 such fragments it collapses toward the bottom of
  the density grid, ending up 4–21× too low (Figure 3a, b). The libraries where that happens are those where the
  mode-free reader loses most against the old one (§5.2).
- **A wrong basin.** In a 12 kb, two-gene stress library with no capture, nascent RNA inflates the deconvolved gDNA
  inside the genes 5×. That inflated level outweighs the true one and becomes "the bulk", and the true gDNA level then
  reads 0.2. Transcript error against truth is 12.4% under the old reader and 36.8% under the mode-free one (Figure 3c).
- **No evidence.** An object with no evidence of its own reads its prior's mean in both readers. Under capture that is
  far above its truth: about ×190 in the old reader and up to ×7 in the new one, where the truth is ×1 (Figures 2d
  and 4).
- **The one-sided null fixes the wrong-basin case but is a wash on the matrix.** It makes every level at or below the
  bulk read 1. That fixes the stress library but loses as much elsewhere as it gains (§6), and it does not address
  sparsity.

**The open question.** Both references need both ends of the ratio. The old one fails when the captured end is
sparse; the new one fails when the uncaptured end is sparse or misidentified. Is there a formulation that is
well-conditioned in both regimes? If not, which failure is the better one to carry into 0.8.0?

---

## 2. Background

### 2.1 Rigel in brief

Rigel is a Bayesian RNA-seq quantifier that separates RNA from genomic-DNA (gDNA) contamination. It works in three
stages:

1. A native scanner tallies fragments onto genomic **objects**:
   - **regions**, the genome partitioned at annotation breakpoints, which hold contained fragments;
   - **boundaries**, the cuts between adjacent regions, which hold the fragments that cross them.
2. **Calibration** deconvolves every object's unspliced fragments into gDNA and RNA, by strand, composition and
   neighbouring evidence. It publishes each object's deconvolved gDNA count `k` on its gDNA support `S`, the expected
   number of fragment placements it offers.
3. A per-locus **EM** assigns RNA to transcripts against a gDNA component.

Every component's likelihood divides by its effective length. The project's binding rule is that there are three
populations per object: gDNA, RNA on +, and RNA on −.

### 2.2 What a capture efficiency is, and why only ratios matter

The efficiencies scale the effective lengths:

- **A transcript's** effective length is its plain length times the share-weighted mean efficiency of the objects its
  fragments land on: its exon pieces, the boundaries between them, and its junctions, priced from the gDNA objects
  beside them.
- **The locus gDNA component's** length is `Σ S_r e_r + Σ M_b e_b` over the locus's regions and boundaries.

Multiplying every efficiency in a locus by one constant multiplies every component's length by it, and the EM's
fractional assignments do not change. That invariance was measured to 1e-15. On the release branch, unit tests now pin
the half that the readers control: one factor on every efficiency is one factor on every length. So the unit is free, and **what matters is the ratio of efficiencies between objects**: chiefly captured against uncaptured,
and partly captured against fully captured.

The two references differ only in which end of that ratio is the unit, and in how the other end is estimated.

---

## 3. The gDNA landscape, and what "basin" and "bulk" mean

### 3.1 What the landscape is

Calibration fits one **landscape** per library: the population's density over `log ρ`, where `ρ = k / S` is a region's
deconvolved gDNA density in fragments per bp.

- **Grid.** A grid of `K = 260` points spanning the data's support.
- **Kernels.** Every training region contributes one Poisson kernel `P(k | ρ S)` in `log ρ`. Its width comes from
  the population's nearest-neighbour spacing (`k = √n`), and it is weighted by how reliably the solve located that
  region's gDNA.
- **Zero counts.** Regions with no gDNA fragment contribute kernels flat below their resolution wall `1/S`.
- **Padding.** One uniform pseudo-region of ignorance is added (Laplace smoothing). The solve (ψ) reads that padded
  density.
- **The reader's version.** The capture reader reads the population's **own** density, without the pseudo-region,
  called `logP_free`.

Read it as a smoothed histogram of every region's gDNA density.

### 3.2 The census: basins, and the bulk

The **census** finds every local maximum of the density, including a grid-edge cell that is higher than its one
neighbour, and cuts the grid at the minimum between adjacent maxima. Each segment is a **basin**, and the basins
partition the grid. Each basin has:

- a **peak**, its maximum;
- a range `[lo, hi]`, bounded by the minima on either side;
- a **mass**, the share of the population inside it.

The mode-free reader then names three kinds of basin:

- **the bulk:** the basin with the largest mass. The premise is that gDNA off the probes always outnumbers gDNA on
  them, since most of the genome is unprobed.
- **enriched basins:** basins lying wholly above the bulk. These are captured levels; partial capture and different
  probe strengths can make several.
- **below-bulk basins:** basins lying wholly below the bulk.

![Figure 1](capture_reference_figures/fig1_landscape_basins.png)

**Figure 1** shows three real landscapes, each on a linear scale (top) and a log scale (bottom).

- **Ladder, 5% gDNA, capture ON.** The bulk (68% of the population) sits at 10^−3.7 fragments/bp. One enriched basin
  (30%) sits at 10^−0.8, about 850× higher: the fully captured level. A 1.8% basin at the grid floor holds regions
  whose deconvolved gDNA is about zero.
- **Test chromosome, 50% gDNA, capture ON.** The bulk is 68%. There are three enriched basins: two partial-capture
  levels (2.9% and 0.5%) and the fully captured level (28%, about 1,084× the bulk).
- **Test chromosome, 5% gDNA, capture OFF.** One basin holds 100%: nothing is enriched.

---

## 4. The two readers

### 4.1 The old reader: reference = the captured level

The old reader is the release foundation's, used until today. It **locates an enriched mode**: the basin above the
largest one that holds the most kernels with at least one fragment. That basin must also pass two tests:

- more than `√n` such kernels;
- a resolution test on their spacing, the **gate**.

The basin's peak becomes the reference `ρ_ref`, the fully captured level. Each object then reads

    e_o = E[ min(ρ / ρ_ref, 1) | k_o, S_o ]

under the padded landscape as prior. That is the posterior mean of its density, clipped at full capture, so it lies
in `[0, 1]`. **No located mode, no correction:** every efficiency is 1.

Two properties follow:

- **Zero-fragment objects.** A plug-in reading, `k/S` against the reference, puts a zero-fragment object at exactly 0.
  A component whose objects all hold no fragment then gets length 0, and the EM meets 0/0. Earlier versions guarded
  this with a floor and then with the landscape posterior above, which reads such an object at the population's
  clipped mean.
- **Everything hangs on locating the captured level.** When capture is sparse, the gate refuses a mode and the
  library gets no correction at all.

### 4.2 The mode-free reader: reference = the uncaptured level (the bulk)

![Figure 2](capture_reference_figures/fig2_reader.png)

The mode-free reader has four parts.

- **The prior (Figure 2a), Efron's two groups.**
  - The **null** is the bulk basin collapsed onto its peak `x_b`, with the basin's whole mass `π0`. Off the probes
    there is one uniform template, and the basin's width is taken as measurement noise.
  - The **alternative** is the population's own density everywhere outside the bulk basin, as it is: the enriched
    spectrum, and any below-bulk basins.
  - No mode is located and nothing is gated. A one-basin landscape makes the prior a single point, so every object
    reads exactly 1.
- **The likelihood.** Each object's own count: `ℓ(x) = k·log(S eˣ) − S eˣ`.
- **The readout.**
  - `e_o = exp( E[ min(x, s(x)) ] − x_b )`: the posterior mean of the **log** density, against the bulk. The bulk
    reads exactly 1.
  - **Saturation per level (Figure 2b):** inside an enriched basin, a density above that basin's peak reads the peak.
    `s(x)` is that peak, or `+∞` outside enriched basins.
  - **The log readout limits the lift of an object without evidence.** A level of prior mass `μ` at `Δ` nats above the
    bulk lifts it by at most `e^{μΔ}`; a linear mean would give `1 + μ(e^Δ − 1)`.
- **The saturation level `S`.** This is the most massive enriched basin's peak over the bulk. A junction's price is
  bounded by `max(S, best touched object)`, which is exactly the old reader's bound of 1 in the old unit. An object
  with no gDNA support reads 1.

**Worked examples (Figure 2c–f):**

| Object | Reads |
|---|---|
| An unprobed intron with no fragment | 1: its zero count over 2 kb is strong evidence of the bulk |
| A tiny piece with no fragment (support 0.3) | The prior, about ×4.4 |
| A partially captured piece below its basin's peak | Its own level, ×13 |
| A fully captured exon above its basin's peak | The saturated level, ×1,084 |

### 4.3 Side by side

| | Old | Mode-free |
|---|---|---|
| Unit | the captured level (located mode) | the uncaptured level (largest basin's peak) |
| Prior | the padded landscape | two groups: bulk point + the population's own density elsewhere |
| Readout | linear mean, clipped at 1 | mean of log, clipped per enriched basin |
| Gate | `√n` located members + resolution test | none |
| Object with no evidence | the prior's arithmetic clipped mean | the prior's geometric mean |
| One-basin (capture-OFF) library | gated: all 1 | prior is a point: all 1 |

---

## 5. What each reference gets wrong, measured

The numbers come from the project's benchmark simulations, scored against the simulator's truth.

### 5.1 The old reference

- **No captured level, no correction.**
  - **Single-probe panel:** the target gene is over-called by +65% (κ 0.99) and +117% (κ 0.70) at 50% gDNA. The
    mode-free reader gives −10% and 0%.
  - **Real cfRNA:** two of three plasma libraries are gated, so their capture goes uncorrected.
- **Objects without evidence read far too high (Figure 4, left).** On the ladder, every unprobed region's true
  efficiency is the bulk's (log 0). The 16,032 unprobed regions holding no gDNA fragment read a median of +1.1 nats,
  rms 3.3. The smallest of them, with support below 5, read a median ×190: the population's clipped mean.
- **Objects with evidence spread widely (Figure 4, right).** Unprobed regions holding fragments spread about ±1 nat
  around the truth: rms 0.54 with 1–2 fragments, and 0.37 with 3 or more.

![Figure 4](capture_reference_figures/fig4_old_vs_new_unprobed.png)

### 5.2 The new reference

**(a) Sparse uncaptured gDNA: the owner's concern, confirmed.** For every saved calibration, the census's bulk is
compared with the true uncaptured gDNA density. That truth is computed as Σ true gDNA fragments over Σ true support,
on regions that overlap no probe according to the simulator's own probe geometry, or on every region when capture is
off.

| True uncaptured gDNA fragments in the library | Calibrations | Bulk ÷ truth |
|---|---:|---|
| ≥ 1,000 | 92 | 0.95 – 1.11 |
| 4 – 601 | 7 | 0.05, 0.06, 0.17, 0.25, 0.93, 1.15, 1.55 |

With little uncaptured gDNA, most regions hold no gDNA fragment at all. Their kernels pile up below their resolution
walls, so the bulk basin's peak slides toward the grid floor (Figure 3b: 22× below the truth). Every captured object's
efficiency is then too high by the same factor, relative to every uncaptured one. These are the libraries where the
mode-free reader loses most against the old one:

| Capture-ON library (test chromosome) | Uncaptured gDNA fragments | Bulk ÷ truth | Transcript error: old → mode-free |
|---|---:|---:|---|
| 10% depth, 0.1% gDNA | 4 | 1.55 | 14.75 → 17.94 |
| full depth, 0.1% gDNA | 27 | 0.25 | 19.22 → 21.36 |
| 10% depth, 1% gDNA | 39 | 0.17 | 26.91 → 30.66 |
| 10% depth, 5% gDNA | 140 | 0.05 | 14.15 → 15.10 |
| full depth, 1% gDNA | 267 | 0.93 (pooled 0.004) | 6.57 → 7.72 |
| 10% depth, 25% gDNA | 601 | 1.15 | 14.24 → 14.37 |
| 10% depth, 50% gDNA | 1,255 | 1.11 | 20.28 → 19.89 |

The association is strong, but it is not yet proof of cause. §7 gives the decisive test.

![Figure 3](capture_reference_figures/fig3_bulk_failures.png)

**(b) The largest basin can be the wrong level (Figure 3c).** The library:

- a 12 kb genome with two genes and no probes;
- heavy gDNA and heavy nascent RNA;
- weak strandedness, κ 0.65.

The truth is uniform gDNA at about 0.057 fragments/bp, which the RNA-free intergenic regions read. Inside the genes,
calibration cannot separate nascent RNA from gDNA at κ 0.65. It books the genic gDNA at 0.16–0.36 fragments/bp, up
to about 5× too high.

With only 13 regions, the landscape has two basins: the true level (44%) and the inflated genic level (56%). The
larger one becomes the bulk, so the true gDNA level reads 0.2.

The efficiencies then shrink lengths unevenly:

- Most of the transcripts' fragments sit on boundaries and junctions, which read 0.23–0.31. So t1 and t2 shrink to
  20% of their length.
- The locus's gDNA component also covers genic regions reading about 1, so it shrinks only to 35%.

The EM therefore sees the transcripts as denser than they are:

| | Old reader | Mode-free | Truth |
|---|---:|---:|---:|
| t1 | 34.8 | 43.1 | 33 |
| t2 | 14.6 | 17.9 | 11 |
| t3 | 1.9 | 1.7 | 3 |
| t4 | 6.8 | 5.4 | 7 |
| Transcript error | 12.4% | 36.8% | |

The old reader read 1 everywhere because its gate found nothing above the largest basin. The companion scenario,
lighter gDNA and κ 0.90, improves from 3.5% to 2.5%.

**Where the error comes from.** The mode-free reader measured faithfully. A deconvolution error became a reference
error, in a population too small to outvote it.

On genome-scale libraries the uncaptured population is far larger. Across the benchmark's contaminated capture-OFF
rows, which should all read 1, the cost is at most 0.7 points (lowg κ0.99 OFF, +0.7). Weakly stranded, nascent-heavy
capture-OFF libraries with few objects remain the risk.

The zero-gDNA capture-OFF rows lose up to 3.4 points through a different mechanism, tracked separately: false gDNA on
tiny boundaries forms a basin *above* the bulk.

**(c) Objects without evidence read the prior's geometric mean.** With 30% of the ladder's landscape enriched, a
zero-fragment unprobed region too small to carry evidence reads up to about ×7: the hump in Figure 4, left. That is
much closer to the truth (×1) than the old reader's ×190, but it is not ×1. Most such objects carry small shares, so
their weight in any length is small.

**(d) One physical level can split into two basins.** For example, odg05 κ0.70 g25, and the RNA-long gap panel's
captured level. Per-level saturation then caps fully captured objects at two levels 1.5× apart. This was found by the
earlier external review and is tracked separately.

### 5.3 The common structure

Both references need both ends of the captured-to-uncaptured ratio.

- **The old one** must locate the captured level; it fails when the captured population is small (small panels:
  exactly the owner's cfRNA use case).
- **The new one** must place the uncaptured level; it fails when that population is sparse (very efficient capture
  with little gDNA) or outweighed (§5.2b).

Which end is sparse depends on the library. A design that assumes either end is always well measured will fail on
some libraries.

---

## 6. The candidate change: the one-sided null

**Definition.** The alternative is restricted to levels **above** the bulk basin. Every level at or below the bulk
basin's top joins the null, so nothing reads below 1.

**Rationale.** Capture only ever adds to a fragment's chance of being kept: in the simulator its weight is the
off-target weight plus binding × probe overlap. So a gDNA density below the uncaptured level is not a capture effect.
It is copy loss, poor mappability or deconvolution error, and it should not shrink any transcript's length.

The same rule reproduces the old reader's safety in both wrong-bulk cases. In a capture-free library whose biased
level outweighs the true one (§5.2b), and in a library whose captured population outweighs its uncaptured one,
everything reads 1.

**Measured** on the full benchmark matrix, 16 panels and 265 conditions, against the integrated reader:

| Stratum | Old reader | Integrated mode-free | One-sided | Change |
|---|---:|---:|---:|---:|
| lowg κ0.99 capture-OFF, contaminated | 6.97 | 7.70 | 6.96 | −0.74 |
| fle200 κ0.70 capture-ON | 8.88 | 9.29 | 8.77 | −0.53 |
| bind100 κ0.50 capture-ON (deferred stratum) | 15.02 | 15.36 | 14.96 | −0.40 |
| Ladder κ0.99 capture-ON (deciding) | 3.87 | 3.87 | 3.90 | +0.02 |
| RNA-short / RNA-long gap κ0.99 capture-ON (deciding) | 13.32 / 4.25 | 13.50 / 5.17 | 13.50 / 5.17 | 0.00 |
| flgl κ0.99 / κ0.70 capture-ON | 9.07 / 10.63 | 9.05 / 10.49 | 9.37 / 10.85 | +0.32 / +0.35 |
| sparse κ0.50 capture-ON (deferred stratum) | 94.42 | 80.14 | 80.48 | +0.34 |
| bind100 κ0.99 capture-ON | 8.58 | 7.85 | 8.72 | +0.88 |
| combo_extreme stress library (§5.2b) | 12.44 | 36.77 | 12.44 | −24.33 |

In-scope tally against the old reader: integrated 7 better / 72 within ±0.3 / 17 worse; one-sided 6 / 73 / 17.

On bind100 the change moves 1,080 of 1,760 regions from just below 1 (median 0.96) to exactly 1. The cost lands on the
isoform split of one capture-cluster block of transcripts; the mechanism is not identified. The one-sided null does
not touch the sparse-null failure (§5.2a), where the bulk itself is misplaced.

---

## 7. Directions the reviewer may weigh

None of these is built or ranked; they are here to be judged.

1. **Estimate the uncaptured level by pooling, not by a peak.**
   - Σ k / Σ S over the null's members, zero-count regions included, stays at 0.8–1.6× the truth where the census
     peak collapses (Figure 3a, orange).
   - It fails badly on one library, full depth with 1% gDNA (0.004×). There the bulk basin holds only 41% of the
     population, so "which objects are uncaptured" is decided by the same census.
   - The membership rule is the crux.
   - An earlier variant used a pooled null level together with a deepest-valley null and was flat on the deciding
     panels.
2. **Estimate the ratio, not an end.** Treat both levels as parameters with uncertainty, for example in a hierarchical
   model over the regions' counts, and integrate over them. Neither end then has to be "located".
3. **Use the better-measured end as the unit.** The unit is free (§2.2), so a library could read its efficiencies
   against whichever end its population measures more precisely.
4. **Known probe support.** Users of capture panels usually have the probe BED. Taken as an optional input, it would
   define the uncaptured population geometrically and remove the census from the reference. This conflicts with the
   design rule that capture is read from the data, per object, and is the owner's call; it is listed for completeness.
5. **Fix the input, not the reader, for §5.2b.** The stress library's failure begins as a deconvolution error (nascent
   RNA booked as gDNA at κ 0.65). Robustness could also come from calibration, upstream.
6. **The decisive experiment for §5.2a.** An oracle arm that gives the reader the true uncaptured density, and nothing
   else, on the sparse libraries. If the losses vanish, the bulk's placement is the cause and priced. A matching
   oracle for the old reader's captured level would price its failure the same way.

---

## 8. Questions for the reviewer

1. **The null's identification.** Is "the largest basin is the uncaptured level" defensible, and under what
   conditions? What is a better identification of the null?
2. **Sparse null.** When uncaptured gDNA is sparse (tens to hundreds of fragments library-wide), how should its level
   be estimated: census peak, pooled rate, or a model over the null's members (e.g. zero-inflated)? How can members be
   chosen without circularity?
3. **Which unit is better conditioned?** Given the EM's invariance to a common factor, and the two sparsity regimes of
   §5.3, is there a formulation that needs neither end located? If not, which end should be the unit, and should that
   be decided per library?
4. **Below-bulk levels.** They are physically impossible as capture effects but measurable (copy loss, mappability,
   deconvolution error). Should the reader allow them (two-sided) or not (one-sided), given the trade measured in §6?
5. **Objects without evidence.** What should such an object read? Both readers return the prior's mean, which lifts
   unprobed objects above their truth: ×190 (arithmetic) or up to ×7 (geometric), against ×1.
6. **The stress library.** Should the reader be robust to calibration's deconvolution errors (for example, by requiring
   the null to be the lowest well-populated level), or should the fix be upstream?
7. **0.8.0.** Ship the integrated mode-free reader with §5.2 as documented limits, adopt the one-sided null, return to
   the captured-level reference, or something else? What would you need to see first?

---

## 9. Evidence and reproduction

| What | Where |
|---|---|
| Release branch | `release-0.8.0` from foundation `3e7e0fbc`, worktree `<worktree>`, uncommitted |
| Release reader | `src/rigel/calibration/capture_efficiency.py`; census `landscape.py`; junction bound `capture_eff_length.py` |
| Old reader | foundation `3e7e0fbc`: `capture_efficiency.py`, `landscape.located_enriched_mode` |
| Prototype directory | `$RIGEL_RUNS/prototypes/2026-10-13_modefree/` |
| Bulk accuracy (§5.2a, Figure 3a) | `bulk_vs_truth.py` → `review/bulk_vs_truth.tsv` (99 calibrations) |
| Figures | `review_figures.py`, run under the release site; data from `dumps_tglb/*.npz` |
| Stress library (§5.2b) | `tests/test_golden_output.py::TestCombinedStress`, scored by a truth script in the session scratchpad |
| One-sided arm (§6) | site `site_tglbm1`; results `out/mf_tglbm1`; comparison `onesided_compare.py`, `arm_tally.py` |
| Integrated arm | site `site_tglbm`; results `out/mf_tglbm`; the release build's own run `out/rel_identity` is identical bit for bit |
| Earlier review of the reader | `docs/dev/RNA_RELEASE_0_8_0_MODEFREE_REVIEW.md` and `RNA_RELEASE_0_8_0_MODEFREE_INDEPENDENT_REVIEW.md` |

**Protocol for every matrix number:**

- scan pinned (`scan.total_threads=1`);
- fractional assignment;
- `OMP_NUM_THREADS=1`;
- seed 0.

Transcript error is `Σ|estimate − truth| / Σ truth`, per stratum as a ratio of sums.
