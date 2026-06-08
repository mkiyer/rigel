# The count-clue density sweep: audit, failure analysis, and redesign

**Status:** root-cause analysis + design proposal (no code changed). 2026-06-08.
**Scope:** the Phase-1 "count clue" in `src/rigel/calibration/density_model.py` (`node_gdna_density`),
its role in the acyclic calibrator, why it produces a near-uniform global gDNA density (and thereby
the in-locus gDNA→RNA leak), and a proposed replacement.

---

## 1. The calibration algorithm (where the sweep sits)

Rigel's calibrator (`calibration/calibrate.py`) is a **single acyclic feed-forward pass** — no EM
loop. It deconvolves the library, per *node*, into gDNA vs RNA. A **node** is a region's *contained*
mass plus its two *boundary sides* (the `CalibrationSubstrate` 3-view). The flow:

```
payload (fractional accumulator: per-region contained counts + per-boundary side counts)
   │  CalibrationSubstrate.from_payload → contained / left-side / right-side views
   │
   ├─ effective lengths (FL geometry): region_eff_len, boundary_side_eff_len, fl_mean
   ├─ fit_strand_balance ............................ rna_sense_frac (κ)         [strand MEAN]
   ├─ strand-clean gdna_frac per region (closed form from κ)
   │
   ├─ node_gdna_density  ◀───────────────────────────  PHASE 1: THE COUNT CLUE / "SWEEP"
   │      → per-region gDNA *density* (frag/eff-bp) + count_evidence (= density·eff_len)
   │
   ├─ fit gDNA / RNA strand Beta-Binomial overdispersions
   │
   ├─ deconv_regions / deconv_sides  ◀───────────────  PHASE 3: JOINT per-node posterior
   │      posterior(gdna_frac) ∝ Beta(count prior FROM PHASE 1) × BB_strand(sense,anti | gdna_frac)
   │      → per-node gDNA / RNA mass
   │
   └─ derive ........................................ gdna_density_global + gdna_geom_len
```

Downstream, `priors.assemble_priors` projects the deconvolved per-node mass to per-locus Dirichlet
scalars and runs `_transport_boundary_flux` (re-attributes deconvolved boundary **mass** to its two
regions ∝ local density ratio × boundary capacity).

**Two clues feed the joint deconvolution (Phase 3):**
- **count clue** — Phase 1's per-region density → a `Beta(a_c, b_c)` prior on `gdna_frac` with
  mean `count_gdna_frac = density·eff_len / unspliced_mass` and **concentration `count_evidence = density·eff_len`** (the expected gDNA count). `joint_deconv.py:94-97`.
- **strand clue** — the Beta-Binomial strand likelihood (gDNA symmetric ½, RNA skewed to κ).

> **Correcting a common misconception.** The sweep is **Phase 1, *before* the joint deconvolution** —
> it *builds the count prior*. It is **not** a post-deconvolution step that imputes nodes that "could
> not be deconvolved." (The only thing that runs after deconvolution is `_transport_boundary_flux`,
> which moves already-deconvolved *mass*; that mechanism works correctly.) Non-strand-observable
> nodes (AMBIG) fall back to *count-only* inside the joint — i.e. they lean on this same Phase-1
> density. So Phase 1 is load-bearing for exactly the hardest nodes.

---

## 2. The sweep algorithm (`node_gdna_density`)

**Count-observability** (`density_model.count_observable_masks`) — where unspliced fragments are
gDNA *by construction* (no mature RNA can be present):
- a **region** is observable ⇔ it carries **no exon bit** (intron-only or intergenic); an exonic
  region's contained unspliced mass is contaminated by mature RNA, so it is **not** observable.
- a **boundary** is observable ⇔ its two regions **share no exon bit** (an exon–intron or
  exon–intergenic *seam*); no single-exon transcript continues across it, so no unspliced *mature*
  RNA crosses — the crossing unspliced mass is gDNA (+ nascent).

**Direct observations** that seed the sweep:
- observable region `r`: `reg_mass[r] = gdna_frac[r] · contained.mass_unspliced[r]`, length `reg_len[r] = region_eff_len[r]`.
- observable boundary `r|r+1`: `bnd_mass[r] = right.mass_unspliced[r] + left.mass_unspliced[r+1]` (raw crossing unspliced mass, treated as gDNA), length `bnd_len[r] = fl_mean`.

**The sweep itself** (`density_model.py:128-165`) — an *alternating region↔boundary propagation*,
per reference, in two passes (forward `from_left`, reverse `from_right`):

```
run_mass = run_len = 0
for region i, left → right:
    if i not first:
        w = weight[i-1]                              # conduit weight of the boundary to the left
        run_mass = w · (run_mass + bnd_mass[i-1])    # carry accumulated evidence across the boundary
        run_len  = w · (run_len  + bnd_len[i-1])
    from_left_mass[i] = run_mass;  from_left_len[i] = run_len
    run_mass += reg_mass[i];  run_len += reg_len[i]  # add this region's own observation
# (reverse pass symmetric → from_right_*)

density[i] = (reg_mass[i] + from_left_mass[i] + from_right_mass[i])
           / (reg_len[i]  + from_left_len[i]  + from_right_len[i])      # global fallback if 0
count_evidence[i] = density[i] · region_eff_len[i]
```

The **conduit weight** `weight[b] = cross_flux[b] / (cross_flux[b] + 1)` (`density_model.py:126`)
is meant to attenuate propagation across low-traffic boundaries: a boundary with no crossing
fragments should not conduct density.

**Intent:** impute the gDNA density at non-observable nodes (exons, AMBIG) by conducting the density
of nearby observable nodes across boundaries, with attenuation, so each node's density reflects its
*local* neighborhood.

---

## 3. The bug

Two structural flaws make the sweep collapse to a single global average instead of a local estimate.

### 3a. The conduit weight does not attenuate at real coverage
`weight = flux / (flux + 1)`. This was tuned for `flux ≈ 0…1`. At realistic depth the per-boundary
crossing flux is **thousands** of fragments, so `weight ≈ 1.0000` at *every* boundary. With `w ≈ 1`
the running accumulators never decay:

```
from_left[i]  →  Σ_{j<i} (reg_mass[j] + bnd_mass[j])     (entire left side of the locus, undamped)
from_right[i] →  Σ_{j>i} (reg_mass[j] + bnd_mass[j])     (entire right side, undamped)
density[i]    →  (Σ_all observable gDNA mass) / (Σ_all observable length)  =  the GLOBAL average
```

So **every** node's density converges to one number. There is **no locality**: an exon's density is
influenced no more by its two flanking boundaries than by a node 40 kb away.

### 3b. Length normalization dilutes the boundary signal
Even with attenuation, the numerator/denominator are mismatched. The boundary-crossing gDNA mass
(the *enrichment* signal — large at exon edges under capture) goes in the **numerator**, but it is
divided by accumulated lengths dominated by the **long, depleted intron** `region_eff_len`
(5–8 kb each). The local enrichment present at the boundaries is length-normalized away to the
global mean.

### 3c. Consequence for the count prior
A non-observable **exon** therefore gets `density ≈ global` (≈ 0.62 frag/bp in the worked example).
Its count prior then has:
- **mean** `count_gdna_frac = density·eff_len / contained_mass ≈ 0.02` (global density × short exon
  eff-len ÷ the *huge* capture-enriched exon mass), and
- **concentration** `count_evidence = density·eff_len ≈ 700–960` — a **fabricated** confidence (it is
  not a count of locally observed gDNA; it is `global_density × eff_len`).

This high-confidence, near-zero prior **fights the strand clue** — which is *correct*
(`strand_gf ≈ true ≈ 0.8–0.94`). The joint posterior compromises at `gdna_frac ≈ 0.4`, so the exon's
gDNA is under-deconvolved ~2×. That under-deconvolved gDNA mass becomes the per-locus prior the EM
receives, so the EM then leaks the corresponding unspliced fragments gDNA→RNA. (The downstream
boundary *transport* works — it depletes introns and fills exons — but it can only redistribute the
already-too-low total.)

---

## 4. Worked example — GENE0037 (locus 36), gdna400 / ss0.99 / capture-on

Tooling: `scripts/debug/gene37_sweep_walk.py` (node-by-node), `gene37_region_boundary_autopsy.py`
(region/boundary vs oracle). All accumulator data validated against oracle truth (template spans):
intron contained gDNA, boundary crossing gDNA, and exon contained mass all match truth exactly — the
**accumulator is correct**; the failure is purely in the Phase-1 sweep.

Node-by-node (abridged; `unspl(+,-)` = contained unspliced counts):

```
R368 intergenic obs=Y eff=19341  mass=0             DIRECT gdna=0       SWEPT dens=0.618 swept_gdna=11958
 └B368 interg→EXON obs=Y flux=10296  conduit_w=0.9999  DIRECT bnd=3246
R369 ex+        obs=n eff= 1141  unspl(12077,19127) mass=31204  (impute) SWEPT dens=0.618 swept_gdna=  705
 └B369 EXON→intron obs=Y flux=12681  conduit_w=0.9999  DIRECT bnd=4206
R370 in+        obs=Y eff= 1667  unspl(9,17) mass=26   DIRECT gdna=18    SWEPT dens=0.618 swept_gdna= 1031
 └B372 intron→EXON obs=Y flux=14412  conduit_w=0.9999  DIRECT bnd=4795
R372 in+        obs=Y eff= 4959  mass=79              DIRECT gdna=60     SWEPT dens=0.618 swept_gdna= 3066
   …every node: conduit_w = 0.9999, SWEPT density = 0.618 ± 0.0001…
```

Per-region density vs truth (exons):

| region | true gDNA density (frag/eff-bp) | **swept density** | count_gf | strand_gf | joint cal_gf | **true gf** |
|---|---|---|---|---|---|---|
| 369 | 21.1 | **0.618** | 0.02 | 0.77 | 0.41 | 0.77 |
| 371 | 20.7 | **0.618** | 0.03 | 0.93 | 0.44 | 0.94 |
| 378 | 20.9 | **0.618** | 0.02 | 0.77 | 0.39 | 0.77 |

What the data shows:
- **conduit_w = 0.9999 at every boundary** (flux 8k–16k) → undamped → uniform 0.618 everywhere (§3a).
- The **boundary crossing gDNA is large and correct** (DIRECT bnd 3k–5.5k per seam, matching oracle) —
  the enrichment signal *is present at the boundaries* — but the sweep dilutes it (§3b).
- The sweep **globalizes both ways**: the empty intergenic flank R368 (`mass=0`) is *hallucinated*
  up to `swept_gdna=11958`; the enriched exon R369 (`mass=31204`, ~94% gDNA) is *crushed* down to
  `swept_gdna=705`. **It cannot distinguish an enriched exon from an empty intergenic region.**
- The error is **systematic and uniform across all exons** (not a subset); introns/intergenic
  deconvolve correctly (`gf=1.0`). The leak lives entirely in the non-observable exon regions.

---

## 5. What the ideal algorithm should do

The exon is **not directly count-observable** (its contained unspliced mass is gDNA + mature RNA),
so its gDNA density must be *imputed*. The design question: **how should an exon obtain its gDNA
density from its immediately-flanking exon–intron boundaries, normalized to its own length, rather
than from an undamped global accumulation?**

Requirements:
1. **Locality.** A region's density must come from its *nearest* observable evidence — its two
   bracketing boundaries (and adjacent observable regions) — with genuine distance attenuation, so
   distant nodes do not dominate.
2. **Correct, count-based normalization.** Boundary crossing gDNA is a *rate* observation: the count
   of fragments crossing a seam = (local gDNA density) × (boundary count effective length). Convert
   it to a density by dividing by the **boundary** effective length (≈ `fl_mean`), **never** by a
   region's length. Then a region's density is the (length-weighted) blend of its bracketing seam
   densities.
3. **Per-side attribution.** The accumulator already splits each crossing fragment into the two
   region sides (`substrate.left[r]`, `substrate.right[r]`). The **exon side** of each flanking seam
   samples the *exon's* edge density directly; under capture the probe coverage is ~uniform across the
   exon, so the edge density ≈ the interior density. Use the exon-side mass for the exon's estimate
   (and the intron-side for the intron — though introns, being observable, use their own contained
   density).
4. **Honest concentration.** For a non-observable region the count prior is an *imputation*. Its
   Dirichlet concentration must reflect the **uncertainty of that imputation** (e.g. the number of
   crossing fragments actually observed at the bracketing seams), **not** `density × full eff_len`.
   A fabricated high concentration is what lets the count prior overrule the (accurate) strand clue.
5. **Graceful endpoints.** Zero gDNA → density 0 (no hallucination onto empty intergenic regions);
   an isolated non-observable node with no observable neighbor → fall back to the global density, but
   that should be the *rare* exception, not the universal outcome.

**Empirical feasibility (already verified).** The flanking exon-side crossing mass alone recovers
**6–8 frag/bp** for GENE0037's exons — vs the true ~21 and the broken 0.618. So the local boundaries
carry the correct *order of magnitude* (~10–30× the global). The signal is recoverable; what is
missing is (a) using it *locally* and (b) an unbiased count→density normalization.

---

## 6. Proposed solution

Replace the undamped global sweep with a **local, boundary-anchored density imputation**, plus an
**honest prior concentration**. Concretely, three changes:

### 6a. Observable regions: keep their own direct density (unchanged)
`density[r] = gdna_frac[r] · contained.mass_unspliced[r] / region_eff_len[r]` for intron/intergenic
regions. These are accurate (the depleted baseline under capture).

### 6b. Non-observable regions: impute from the *bracketing seams only*
For an exon (or AMBIG) region `r`, anchor the density on the **two boundaries that bracket `r`**, using
the **exon-side** crossing gDNA converted to a density via the boundary effective length:

```
seam_density(boundary b on r's side) = side_gdna_mass_on_r / boundary_side_eff_len[r]
density[r] = mass-weighted mean of r's left-seam and right-seam densities
```

No locus-wide accumulation; no decay chain. A region sees *only* its two neighbors. (For a region
with one non-seam side — e.g. an internal exon–exon junction, which is not observable — use the one
available seam, then the adjacent observable region, then the global fallback, in that order.)

> The exact count→density constant needs derivation so the estimate is **unbiased** (the naive
> per-side-mass / `boundary_side_eff_len` form lands ~0.3× of truth because edge-crossing
> under-samples the interior; a crossing-*count* / `fl_mean` form over-samples). The right form is a
> short FL-geometry derivation — analogous to the existing `effective_length.py` count-effective-length
> derivation — and must be validated on the benchmark so an unstranded/clean locus reproduces its
> known density. This is the one open quantity; §4's tooling makes it directly checkable.

### 6c. Make the imputed prior's concentration honest
For a non-observable region, set the count-prior concentration to the **count of crossing fragments
actually observed** at its bracketing seams (the real evidence), not `density × region_eff_len`. An
exon with, say, ~15k crossing fragments still gets a real (large) concentration — but one *grounded in
observations*, and the mean is now in the right place (§6b), so the count prior and the strand clue
**agree** instead of fighting. Where evidence is genuinely thin, the concentration is small and the
strand clue rightly governs (this is what the Phase-1 comment *claims* happens today but does not).

### Why this fixes the leak
With §6a–c, each exon's count prior has the correct order-of-magnitude mean (`count_gf` ≈ the true
~0.8 rather than 0.02) and an honest concentration, so the joint posterior lands near the truth
(`gdna_frac ≈ 0.8`) instead of 0.4. The per-locus gDNA prior handed to the EM is then ~correct, and
the in-locus gDNA→RNA leak collapses. The downstream `_transport_boundary_flux` (already correct)
then redistributes the *correct* total gDNA across exon vs intron.

### What is removed
- The decay-chain accumulation in `node_gdna_density` (the `from_left`/`from_right` running sums and
  the `weight = flux/(flux+1)` conduit). It is replaced by the local two-seam estimate.
- The fabricated `count_evidence = density · region_eff_len` for non-observable regions, replaced by
  the observed crossing-fragment count.

### Validation plan
1. **Unit:** on a synthetic locus with known uniform gDNA density, the imputed exon density equals
   the true density (within tolerance) for both stranded and unstranded inputs.
2. **GENE0037 (this doc):** re-run `gene37_sweep_walk.py` — exon swept density should jump from 0.618
   to ~20, `count_gf` from 0.02 to ~0.8; intergenic flanks stay near 0 (no hallucination).
3. **Benchmark:** `evaluate_suite.py` on `gdna_benchmark_5mb` — net gDNA→RNA leak in the
   capture-on / high-gDNA conditions should drop sharply, **without** inflating the gdna_none false
   gDNA (no new RNA→gDNA siphon) and without disturbing the capture-off conditions. Golden tests
   regenerated; full suite green.

### Open questions / for discussion
- The exact unbiased count→density normalization constant (§6b note) — derive from FL geometry, or
  calibrate empirically against the benchmark?
- Should an exon also borrow from its adjacent *introns'* observed baseline (a two-component "baseline
  + capture-excess" model), or is the bracketing-seam estimate sufficient on its own?
- How to treat AMBIG regions (overlapping opposite-strand genes), which are non-observable for a
  *different* reason — the same local-seam imputation should apply, but worth a dedicated test.

---

## Appendix — reproduce
```bash
SUITE=/Users/mkiyer/Downloads/rigel_runs/gdna_benchmark_5mb
COND=gdna_gdna400_ss_0.99_nrna_none_capture_on
python scripts/debug/gene37_sweep_walk.py --index $SUITE/rigel_index \
  --bam $SUITE/$COND/sim_oracle.bam --ref chr_syn --start 1659774 --end 1701056
python scripts/debug/gene37_region_boundary_autopsy.py --index $SUITE/rigel_index \
  --bam $SUITE/$COND/sim_oracle.bam --ref chr_syn --start 1659774 --end 1701056 \
  --gtf $SUITE/reference/genes.gtf
```
Key code: `calibration/density_model.py` (sweep), `calibration/joint_deconv.py` (count prior use),
`calibration/effective_length.py` (FL geometry), `calibration/priors.py` (`_transport_boundary_flux`).
