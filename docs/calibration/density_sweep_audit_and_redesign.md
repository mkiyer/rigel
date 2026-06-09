# The count-clue density sweep: audit, failure analysis, and redesign (v2)

**Status:** root-cause analysis + design proposal (no calibration code changed). 2026-06-08.
**Scope:** the Phase-1 "count clue" in `src/rigel/calibration/density_model.py` (`node_gdna_density`) —
why it collapses to a global density (causing the in-locus gDNA→RNA leak), and a redesign that
handles **runs of consecutive non-observable regions**, the **full region-signature case space**, and
the **effective-length contraction**, which is the dominant error.

> **v2 changes (from v1):** added the consecutive-non-observable-run problem and the directional
> inward sweep; the full 16-signature / region-pair observability case table; a second imputation
> method (DNA *fraction* vs density); the empirical finding that **intergenic contained mass is zero**
> (so a global average is a deflated, bad fallback); and a dedicated treatment of the **IPR
> effective-length contraction** and its uncertainty.

---

## 1. The calibration flow (where the sweep sits)

Single acyclic feed-forward pass (`calibration/calibrate.py`). A **node** = a region's *contained*
mass + its two *boundary sides* (`CalibrationSubstrate` 3-view).

```
payload (fractional accumulator) → substrate (contained / left-side / right-side)
  ├─ effective lengths (FL geometry): region_eff_len, boundary_side_eff_len, fl_mean
  ├─ fit_strand_balance → rna_sense_frac (κ)
  ├─ strand-clean gdna_frac per region (closed form from κ)
  ├─ node_gdna_density  ◀── PHASE 1: the count clue / SWEEP → per-region gDNA density + count_evidence
  ├─ fit gDNA / RNA strand overdispersions
  ├─ deconv_regions / deconv_sides ◀── PHASE 3: joint posterior(gdna_frac) ∝ Beta(count prior) × BB_strand
  └─ derive → gdna_density_global + gdna_geom_len      (later: priors.assemble_priors → IPR eff-len + transport)
```

The sweep is **Phase 1, before deconvolution** — it builds the **count prior** (a `Beta` on
`gdna_frac` with mean `density·eff_len / unspliced_mass` and concentration `count_evidence = density·eff_len`,
`joint_deconv.py:94-97`). It is *not* post-deconvolution imputation. The per-region gDNA mass it
implies also feeds the **IPR effective-length** in `assemble_priors` — so the sweep drives **two**
downstream quantities: the prior *count* and the gDNA *effective length*.

---

## 2. Observability — the full case space

**Region** is count-observable ⇔ it has **no exon bit** (intergenic / intron-only): its unspliced
mass is gDNA by construction. **Boundary** is count-observable ⇔ its two regions **share no exon bit**.
The 4 signature bits are `{intron+ 0x8, intron- 0x4, exon+ 0x2, exon- 0x1}`.

Region observability (16 signatures, condensed to the realistic ones):

| signature | class | region observable? |
|---|---|---|
| `intergenic` (0000) | intergenic | **Y** (direct) |
| `in+` / `in-` / `in+\|in-` | intron | **Y** (direct) |
| `ex+` / `ex-` | exon | n (impute) |
| `ex+\|ex-` | exon, AMBIG | n (impute) |
| `ex+\|in+`, `ex-\|in-` (alt-splice, same strand) | exon | n (impute) |
| `ex+\|in-`, `ex-\|in+` (alt / antisense) | exon, AMBIG | n (impute) |

Boundary observability (computed from the real code; `NON` = shared exon bit ⇒ not observable):

```
 L \ R        interg  in*   ex+   ex-   ex+|ex-   ex+|in-   ex-|in+   ex+|in+   ex-|in-
 interg/in*    obs    obs   obs   obs    obs       obs       obs       obs       obs
 ex+           obs    obs   NON   obs    NON       NON       obs       NON       obs
 ex-           obs    obs   obs   NON    NON       obs       NON       obs       NON
 ex+|ex-       obs    obs   NON   NON    NON       NON       NON       NON       NON
```

**Consequences (the crux of v2):**
- Intergenic & intron regions are **always** observable; every boundary touching them is observable.
- The **only** non-observable boundaries are exon–exon seams **sharing a strand** (`ex+|…` next to `ex+|…`).
- `ex+ | ex-` (a sense exon abutting an antisense exon) **is observable** — no single transcript's
  mature RNA crosses it.
- A **run** of non-observable regions forms wherever consecutive regions share a strand's exon bit
  (alt-splice on one strand; or an antisense overlap `ex+ | ex+|ex- | ex-`). Inside a run, the interior
  regions have **both** boundaries non-observable.

---

## 3. The bug in the current sweep

`node_gdna_density` (`density_model.py:128-165`) propagates observable density across boundaries with a
running accumulator decayed by a **conduit weight** `w = cross_flux / (cross_flux + 1)`.

1. **`w` never attenuates.** At real depth the per-boundary crossing flux is thousands, so
   `w ≈ 1.0000` everywhere. The running sums never decay → every node's density →
   `Σ observable gDNA mass / Σ observable length` = a single **global** value. No locality.
2. **Length mismatch.** The boundary-crossing gDNA (the enrichment signal) is divided by lengths
   dominated by the long, depleted **intron** eff-lengths → diluted to the global mean.
3. **Intergenic blindness (new, measured).** Genome-wide on the benchmark, **all 101 intergenic
   regions have exactly 0 contained mass** (3.19 Mb / ~64% of the genome) — the accumulator does not
   deposit deep-intergenic gDNA into region contained mass, yet truth has ~7.5k gDNA fragments per
   intergenic region. So the global density is **deflated** (huge empty length in the denominator), and
   the gDNA *baseline* is **not observable from intergenic contained mass** — only at gene-flanking
   boundaries and (depleted) introns.

Result: a non-observable exon gets `density ≈ global`, so its count prior is
`mean ≈ 0.02` with a **fabricated concentration** `count_evidence ≈ 700–960` that overrules the
*correct* strand clue → joint `gdna_frac ≈ 0.4` (truth ~0.8) → the under-deconvolved gDNA mass becomes
a too-low per-locus prior → the EM leaks the matching unspliced fragments gDNA→RNA.

---

## 4. Worked examples

Tooling: `scripts/debug/gene37_sweep_walk.py`, `gene37_region_boundary_autopsy.py`,
`build_stress_overlap_scenario.py`. Accumulator validated against oracle (template spans): intron
contained, boundary crossing, and exon contained masses all match truth — the failure is purely the
Phase-1 sweep.

### 4a. GENE0037 (capture-on, 4:1 gDNA) — uniform global, all exons
Every node's swept density ≈ 0.618 (`conduit_w = 0.9999`); exons true ~21 frag/bp → count_gf ≈ 0.02;
strand_gf ≈ true ≈ 0.8–0.94; joint cal_gf ≈ 0.4. Empty intergenic flanks hallucinated up
(`swept_gdna ≈ 12k` from `mass = 0`); enriched exons crushed down. Systematic across all exons.

### 4b. Antisense-overlap stress — a run with a truly-interior region
`ex+ (R7) | EXON→EXON NON-obs | ex+|ex- (R8) | EXON→EXON NON-obs | ex- (R9)`:

```
R7 ex+      obs=n   left boundary B6 (interg→EXON) OBSERVABLE  → reachable from LEFT
R8 ex+|ex-  obs=n   BOTH boundaries (B7,B8) NON-observable     → reachable from NEITHER side
R9 ex-      obs=n   right boundary B9 (EXON→interg) OBSERVABLE → reachable from RIGHT
```

R8 can only be imputed by sweeping **inward** from R7 (anchored on B6) and R9 (anchored on B9). The
current sweep again returns the global (≈ 2.456) for all of R6–R10, hallucinating onto the empty
intergenic flanks. The observable edge boundaries B6/B9 carry the real local density (~14 frag/bp,
matching the true uniform density) — the signal the redesign must route inward.

---

## 5. What the ideal algorithm must do

The exon's gDNA density must be **imputed** (its contained mass is gDNA + mature RNA). Requirements:

1. **Locality with real attenuation.** A region's estimate comes from its *nearest* observable
   evidence; influence decays with distance. (The current `flux/(flux+1)` does not attenuate.)
2. **Count-based normalization.** Boundary crossing gDNA is a *rate*: convert to density by dividing by
   the **boundary** effective length (≈ `fl_mean`), never by a region's length.
3. **Directional inward sweep over runs.** Within a run of consecutive non-observable regions, impute
   from the **leftmost** region (whose left boundary side is observable) rightward, and from the
   **rightmost** region (right side observable) leftward; combine the two directions at each interior
   region. Interior regions whose both *original* boundaries are non-observable inherit the estimate
   carried in from the run's observable edges (this is the *legitimate* core of a "sweep" — local and
   directional, not a global accumulation).
4. **Per-side attribution.** Use the **region-facing side** of each anchoring boundary
   (`substrate.left[r]` / `substrate.right[r]`), which samples that region's edge density.
5. **Honest concentration.** A non-observable region's count prior is an imputation; its Dirichlet
   concentration must reflect the **observed crossing-fragment count** that supports it, not
   `density × full eff_len`. Thin evidence ⇒ weak prior ⇒ the strand clue rightly governs.
6. **No global-average fallback.** Empirically the global is deflated by 3.19 Mb of empty intergenic
   length and mixes depleted with enriched; it is a *bad* default. Fall back, in order, to: the nearest
   observable boundary/region density → an intron baseline → and only as a last resort a conservative
   prior length with near-zero concentration (add little to the prior rather than wrong mass).

**Empirical feasibility.** The flanking observable boundaries already recover the right order of
magnitude (GENE0037 exons: 6–8 frag/bp from boundary sides vs true ~21 and broken 0.618; antisense
run edges: ~14 vs true ~14). The signal is local and present; what is missing is using it locally and
normalizing it correctly.

---

## 6. Proposed redesign

### 6a. Observable regions — direct (unchanged)
`density[r] = gdna_frac[r] · contained_unspliced[r] / region_eff_len[r]` for intron/intergenic.
(Intergenic contributes 0 today — see §3.3; treat its density as the **intron-baseline**, not 0, and
flag the intergenic-deposit gap for a separate fix.)

### 6b. Non-observable regions — local, directional, two-anchor imputation
For each **run** of consecutive non-observable regions bracketed by observable boundaries:
- **Left-to-right pass:** seed at the run's left edge with the observable boundary's region-facing
  side density `d = side_gdna_mass / boundary_side_eff_len`; carry it inward, updating each region from
  the boundary on its left (using that boundary's region-facing gDNA, when observable, else the carried
  estimate) with distance attenuation.
- **Right-to-left pass:** symmetric from the run's right edge.
- **Combine** the two directional estimates per region (e.g. evidence-weighted mean by the crossing
  counts each direction carries). Interior regions with no observable boundary of their own are thus
  filled from both edges of the run.

Two ways to express the per-anchor estimate (both worth trying):
- **(A) Density:** `boundary region-side gDNA mass / boundary_side_eff_len` → a frag/eff-bp density,
  applied as the region's gDNA density (then `count_gdna_frac = density·eff_len / mass`).
- **(B) DNA fraction:** at an intron→exon (or intergenic→exon) seam the crossing fragments are
  *unspliced* (gDNA, after strand-cleaning) **and** *spliced* (mature RNA). Their ratio
  `gdna_unspliced / (unspliced + spliced)` is a **gDNA fraction** that can be applied directly to the
  region's *own total* contained counts: `gdna_mass ≈ region_total · boundary_dna_fraction`. This
  sidesteps the density/length normalization entirely and may be more robust where lengths are
  uncertain. (Tradeoff: assumes the seam's DNA fraction transfers to the region interior; needs
  testing, especially when the two flanking seams disagree.)

### 6c. Honest concentration
Set the non-observable region's count-prior concentration to the **crossing-fragment count** observed
at its anchoring boundaries (real evidence), not `density · region_eff_len`. The mean is now ~correct
(6b) and the concentration is grounded, so count and strand **agree** instead of fighting.

### 6d. The effective-length contraction (the make-or-break)
The IPR (`priors.assemble_priors`) contracts the gDNA component's effective length from the per-region
gDNA *mass* concentration. This length is a **divisor** in the per-fragment EM weight `θ/L`; when gDNA
(~350 bp, spread) competes with RNA (~spliced length), a wrong `L` is a first-order bias — worse than a
wrong prior count, which the EM can outvote.

Key points for the redesign:
- **It is coupled to 6b.** The IPR is built from the same per-region gDNA mass the sweep imputes.
  Getting the density imputation right *automatically* contracts the eff-len correctly (gDNA mass
  concentrates on the captured exons → short eff-len). So 6b is the primary lever for the eff-len too.
- **Model the uncertainty, don't hide it.** For imputed regions the gДНК mass is uncertain, so the IPR
  length is uncertain. Rather than dividing by a single hard `L`, the imputation uncertainty should
  enter as **prior strength** (6c) — a region we cannot resolve contributes *little* to both the prior
  count and the eff-len, deferring to the EM + strand, instead of injecting a confident wrong `L`.
- **Principled fallback length.** When a region is genuinely unidentifiable, do **not** contract toward
  a global; regularize the eff-len toward a sensible local prior — the **geometric span** is the
  uninformative null (no contraction; the existing Laplace `+1` already does this for `G→0`), and under
  *suspected* capture the **exon/spliced footprint** is the informed prior. Choosing between them
  should be driven by the local evidence strength, not a global constant.

### 6e. Two error modes — design stance
- **Prior count.** Adding *no* mass where we have no evidence is acceptable (a weaker prior; the EM
  still solves it). Adding *wrong* mass is harmful. So bias toward honest (low) concentration over
  confident imputation when evidence is thin.
- **Effective length.** This is the dangerous one: an unmodeled, confident divisor. Treat its
  uncertainty explicitly (6d). "We don't know" must translate to *less contraction confidence*, not a
  fabricated short or global length.

### What is removed
The undamped `from_left`/`from_right` accumulation and `w = flux/(flux+1)`; the fabricated
`count_evidence = density · region_eff_len` for non-observable regions.

---

## 7. Empirical validation plan

Generate stress genomes covering **every** region-pair signature combination and exercise the
imputation against oracle truth (`scripts/debug/build_stress_overlap_scenario.py` is the seed; extend
to a generator):
- **Single neighbor pairs:** for each (left-sig, right-sig) in §2, a 2-region locus — check the
  imputed density/fraction and eff-len vs truth, observable and non-observable cases.
- **Runs:** 3–5 consecutive non-observable regions (alt-splice on one strand; antisense overlap;
  mixed), with an interior region reachable only inward — check both directional passes and their
  combination, including conflicting left/right anchors.
- **Capture vs no-capture × stranded vs unstranded:** confirm enriched exons get high density,
  depleted introns/intergenic get baseline, and unstranded/AMBIG degrade gracefully (no hallucination).
- **Truly-unidentifiable:** overlapping opposite-strand exons at low strand-specificity — confirm the
  fallback adds little (honest weak prior) rather than a confident wrong value.
- **Regression:** `evaluate_suite.py` on `gdna_benchmark_5mb` — net gDNA→RNA leak drops in capture-on /
  high-gDNA without inflating gdna_none false gDNA; golden regenerated; full suite green.

Targets on the worked examples: GENE0037 exon swept density 0.618 → ~20, `count_gf` 0.02 → ~0.8;
antisense-run interior region R8 imputed to ~true from both edges; intergenic flanks stay near baseline
(no hallucination).

## 8. Open questions
- The exact unbiased count→density constant (per-side mass / `boundary_side_eff_len` ran ~0.3× of truth
  on GENE0037; a crossing-count / `fl_mean` form over-samples) — derive from FL geometry and validate.
- Density (6b-A) vs DNA-fraction (6b-B): which is more robust, and how to combine conflicting left/right
  anchors over a run.
- The **intergenic contained-mass = 0** behavior (§3.3): intended (genic-proximal-only calibration) or
  a deposit gap? It deflates the global density and removes the baseline signal — fix or formally model.
- AMBIG regions (opposite-strand overlap) at low strand-specificity: the genuinely-unidentifiable
  floor — quantify the best achievable and ensure the fallback hits it.

---

## Appendix — reproduce
```bash
SUITE=/Users/mkiyer/Downloads/rigel_runs/gdna_benchmark_5mb
COND=gdna_gdna400_ss_0.99_nrna_none_capture_on
python scripts/debug/gene37_sweep_walk.py --index $SUITE/rigel_index \
  --bam $SUITE/$COND/sim_oracle.bam --ref chr_syn --start 1659774 --end 1701056
python scripts/debug/build_stress_overlap_scenario.py
python scripts/debug/gene37_sweep_walk.py --index /tmp/rigel_sim/antisense_run/index \
  --bam /tmp/rigel_sim/antisense_run/antisense_run_oracle.bam --ref antisense_run --start 4500 --end 8500
```
Key code: `density_model.py` (sweep + observability masks), `joint_deconv.py` (count-prior use),
`effective_length.py` (FL geometry), `priors.py` (IPR eff-len + `_transport_boundary_flux`),
`signature.py` (4-bit signatures).
