# Phase B — strand-aware mature subtraction + residual deconvolution: dry-run plan

**Status:** dry-run plan, 2026-06-12. Second phase of [`deconvolution_roadmap.md`](deconvolution_roadmap.md),
building on Phase A ([`phaseA_mature_imputation_plan.md`](phaseA_mature_imputation_plan.md), shipped
`83b3610`: `mature_density.py` produces per-strand `M₊`/`M₋`). Phase B subtracts that mature from the
region's contained unspliced counts — **strand-aware**, by κ — and feeds the residual into the
deconvolution. This doc walks the implementation and, critically, surfaces a **phasing coupling** the
dry-run revealed and the **unit-consistency** risk, then recommends the scope.

## 1. The subtraction math

Per region, observed contained unspliced counts `U_pos`, `U_neg` (genomic orientation). A mature transcript
on strand `s` reads its **sense** (transcript-strand-matching genomic strand) at rate κ and antisense at
`1−κ` — the same convention `deconv_regions` already uses (`sense = neg if NEG else pos`). So:

```
mat_pos = κ·M₊ + (1−κ)·M₋        # + mature reads genomic+ at κ; − mature reads genomic+ at 1−κ
mat_neg = (1−κ)·M₊ + κ·M₋
U'_pos  = max(U_pos − mat_pos, 0)
U'_neg  = max(U_neg − mat_neg, 0)
```

`U' = U'_pos + U'_neg` is the residual unspliced mass (`gDNA + nascent₊ + nascent₋`). The `max(·,0)` clamp
handles the zero-mass slivers (locus-21 regions 227/232: `U≈0`, tiny run-filled `M` → `U'=0`) and any
over-subtraction; **κ is protocol-agnostic** (forward κ→1, reverse κ→0 — this dataset κ≈0.01).

**Contained view only.** Boundary crossing-*unspliced* is already mature-free (a mature fragment crossing a
junction is *spliced*, not unspliced — `deconvolution_implementation.md` §1). So `deconv_sides` is
**unchanged**; only the contained view is mature-subtracted.

## 2. Mass bookkeeping — keep the blend consistent

The deconvolution splits the contained mass into gDNA/RNA via `frac = w·g_strand + (1−w)·g_count`, then
`gdna = frac·mass`. The mature subtraction must not break the blend's single denominator. Two equivalent
options; **recommend (a)** (smaller diff):

- **(a) fraction-of-`U`, strand reads the residual.** Keep `mass = U` and `gdna = frac·U`. The *strand*
  module estimates the gDNA fraction of the **residual** `g'_strand` (from `U'_pos/U'_neg`), converted to a
  fraction of `U`: `g_strand_U = g'_strand · (U'/U)`. `g_count` stays a fraction of `U`. Both fractions of
  `U` → the blend and the `frac·U` mass split are unchanged; mature stays implicitly in the RNA remainder
  `(1−frac)·U`. **Mass is conserved automatically** (no separate mature term).
- **(b) split the residual explicitly.** `mass = U'`; `gdna = frac·U'`; `rna = (1−frac)·U' + M + spliced`.
  Conserves (`U' + M = U`), but reworks the mass plumbing in `_deconv_per_node`. More invasive.

## 3. THE PHASING COUPLING (the dry-run's key finding) — read before scoping

The roadmap's Phase B said "retire the splice_junction fraction; failing regions recover." The dry-run shows
**these cannot both happen in Phase B without a regression**, because:

- The splice fraction is the *current* (flawed) AMBIG recovery: it lifts region 226 to **0.37**.
- The *new* AMBIG recovery is mature-subtraction **+** strand-on-the-residual **+** the per-strand-mature
  prior. The count module's gDNA estimate is the **capture-depleted boundary density** — unchanged by the
  subtraction (subtracting mature rescales the fraction but not the depleted gDNA *mass*). So an AMBIG
  region left **count-only** after retiring the splice fraction falls to **0.26** (the raw depleted
  density) — **worse than today's 0.37**.
- What actually recovers AMBIG region 226 to ~0.54 is the **strand on the residual**: a no-nascent AMBIG
  residual is *symmetric* → strand reads it as gDNA → `g'≈1`. **But** a symmetric residual can equally be
  *balanced nascent on both strands* (the `test_ambig_scenario`: gDNA=0+nascent, both genes expressed).
  Naively reading symmetric→gDNA there **fails that test**. Distinguishing the two needs the **per-strand
  mature expression** (both strands have mature ⇒ both likely have nascent ⇒ symmetric = balanced nascent,
  not gDNA) — which is **Phase C** (silent-strand override + directional prior).

**Therefore the splice fraction can only be retired once the strand-on-residual + per-strand-mature prior
(Phase C) is in place.** B and C are more coupled than the roadmap assumed.

### Recommended scope (decision needed)

**Phase B (recommended): the subtraction, made correct, wired to feed the *strand* channel — no splice
retirement, no AMBIG change.**
- Compute `M₊/M₋`, subtract strand-aware → `U'_pos/U'_neg` (clamped).
- Feed the residual to the **strand** module (`g'_strand` on `U'`, option 2a). Single-strand exons
  (`POS`/`NEG`, `w→1`) get a *cleaner* gDNA/nascent split (mature contamination removed). AMBIG (`w=0`)
  is **unchanged** — still `g_count` = the splice fraction (no regression).
- Validate: subtraction accuracy (§5), single-strand exons hold/improve, **no regression**, golden refresh,
  `test_ambig_scenario` still passes.
- **Phase C** then: make AMBIG strand-informative on the residual using the per-strand mature
  (silent-strand override, directional prior, neighbour carry-over) and **retire the splice fraction** —
  this is where the flagship AMBIG leak actually drops.

**Alternative: merge B+C** into one "residual deconvolution" phase (subtraction + strand-on-residual for all
regions incl. AMBIG with the per-strand-mature prior + retire splice). Coherent, no intermediate state, but
a larger single phase and it folds the `test_ambig`-critical prior in at once. *I lean to keeping them
separate (above) so the subtraction is validated in isolation before the AMBIG prior logic lands.*

## 4. Dry-run walkthrough (recommended scope; `calibrate.py`)

Insert after the FL/strand fits, before/within the region deconvolution:
1. `md = mature_density(substrate, region_arrays, region_eff_len_rna, rna_fl_mean)` (Phase A).
2. strand-aware subtract (§1) → `up_pos`, `up_neg` (and `u_resid = up_pos+up_neg`).
3. `deconv_regions(...)` gains the residual strand counts: the **strand** module orients `up_pos/up_neg`
   (not the raw `c.n_unspliced_pos/neg`); the count path (`g_count`) is unchanged. Inside
   `_deconv_per_node`, `g_strand` is computed on the residual and rescaled by `U'/U` (option 2a) so the
   blend stays a fraction of `U`. `deconv_sides` unchanged.
4. `region_splice_gdna_frac` stays (retired in Phase C).

`mature_density` is called once; `region_eff_len_rna` and `rna_fl_mean` already exist in `calibrate`.

## 5. Step 0 (do FIRST) — unit-consistency of `M` vs `U`

The subtraction `U − M` is only valid if `M` (built from the **spliced boundary flux** × `E_mu/fl_mean`)
and `U` (the **contained** count) are in the **same accumulator units**. Both are accumulator fragment
counts, so they *should* match — but the per-region check shows scatter:

- AMBIG regions match well: 226 residual 0.54 vs oracle 0.59; 231 0.57 vs 0.57.
- A gene-edge single-strand exon over-subtracts: region 222 (first exon of GENE0023, **one** bounding
  junction, the other side intergenic) residual 0.49 vs oracle 0.63 — `M₋` ~40% high in substrate units.
- The substrate `U` is itself ~20–34% below my **read1** oracle totals (the oracle counts whole read1
  fragments by searchsorted start; the accumulator deposits fractional contained mass) — so **validate `M`
  against the substrate-consistent contained mature, NOT the read1 oracle.**

**Action:** before wiring, a focused check (extend `scripts/debug/phaseA_issueD_mature_density.py`) that
`M` recovers the *substrate's* contained mature per region, especially for **single-junction / gene-edge
exons** (one anchor, no averaging) and regions adjacent to the AMBIG zone. If a systematic per-region bias
appears (e.g. single-junction exons), correct the geometry (it is pure — likely the `n_junctions`
edge-handling or `E_mu` at gene-edge exons) before relying on the subtraction.

## 5.5. Step 0 OUTCOME (2026-06-12) — units consistent; and B should MERGE into C

`scripts/debug/phaseB_step0_unit_consistency.py` scanned a **mature-only** BAM (the accumulator's
contained-unspliced of mrna-only reads = substrate mature, exact units) and compared to `M`:

- **Units consistent, geometry accurate where covered.** On the locus-21 leak regions, `M /
  substrate_mature` = 226:1.00, 231:0.96, 236:0.98, 222:1.16, 244:0.94, 246:1.11. Region 222's earlier
  "40% over" was a **read1-oracle artifact** — it is 16% in substrate units. The subtraction `U − M` is
  valid on exactly the AMBIG leak regions.
- **Global 29% under = a coverage gap** (`total M / total substrate mature = 0.708`, p10 of the per-region
  ratio = 0). Regions with **no eligible bounding junction** get `M=0` — predominantly **single-exon
  transcripts**, which signal #1 *fundamentally cannot see* (no junction → no spliced anchor). Across
  classes: POS ratio 0.56, NEG 0.83, AMBIG 0.58 (4 AMBIG regions under-covered, 26k mature missed → a
  **Phase-C input**: those need the gDNA crossing / neighbour / a floor).

**The decisive finding: mature subtraction has NO benefit for POS/NEG regions.** A strand-observable
region's strand posterior on the **raw** counts already separates symmetric gDNA from **all** stranded RNA
(mature *and* nascent) — generative-model §4 ("the gDNA prior is robust to the RNA composition; mature
cancels"). So removing mature before the strand is redundant there, and with the ~10–16% scatter it only
*risks harm* via over-subtraction (region 222: 0.54→0.49, away from oracle 0.63). The subtraction only
helps where the strand **cannot** separate gDNA from RNA — **AMBIG regions**.

**⇒ Phase B (subtraction feeding the single-strand channel) provides no leak benefit and is not separable
from the AMBIG resolution. RECOMMEND merging B into C** as one "AMBIG residual resolution" phase: subtract
mature from AMBIG regions, run the strand on the AMBIG residual with the per-strand-mature prior, retire
the splice fraction. (`mature_density` from Phase A stands; `scripts/debug/phaseB_step0_unit_consistency.py`
is the unit-consistency record.) **The §1–§4 subtraction math + mass bookkeeping below carry into the
merged phase, applied to AMBIG regions only.**

## 6. Risks

- **Unit consistency / per-region scatter** (§5) — the #1 risk; resolve in Step 0. The strand-on-residual
  must be **robust to imperfect subtraction** (it re-estimates from the residual's symmetry, not the raw
  `U'/U`, which softens subtraction error).
- **Single-strand improvement is not guaranteed** — region 222's over-subtraction could make it *worse*
  (0.54→0.49) if the subtraction bias isn't fixed first. Hence Step 0 gates Phase B.
- **Clamp distortion** — `max(U'_pos,0)` per strand can skew the residual's strand balance when one strand
  clamps; in practice it binds only at zero-mass slivers and gross over-prediction. Log where it binds.
- **Mass conservation** — option (2a) conserves by construction (`gdna+rna = frac·U + (1−frac)·U + spliced
  = U + spliced`); add an invariant check.
- **Golden change** — Phase B wires into `calibrate`, so calibration outputs shift; regenerate golden
  (`pytest --update-golden`) and eyeball the diff (single-strand exon gDNA up, AMBIG ~unchanged).
- **`test_ambig_scenario`** — must still pass (AMBIG unchanged in the recommended scope). If it moves, the
  subtraction is leaking into the AMBIG count path — investigate.

## 7. Validation

1. **Subtraction accuracy** (Step 0) — `M` vs substrate contained mature, per region, on the clean
   multi-exon scenario + locus 21; the gene-edge / single-junction cases specifically.
2. **Per-region gDNA fraction** on locus 21 — single-strand exons (222/244/246) should hold or improve
   toward oracle; AMBIG (226/231/…) unchanged (Phase C territory).
3. **No regression** — full calibration unit suite + `test_ambig_scenario`; golden refreshed and reviewed.
4. **Debug script** `scripts/debug/phaseB_residual_deconv.py` (committed): per-region raw-vs-residual strand
   split + the gDNA fraction, joined to oracle.

## 8. Out of scope for Phase B (→ Phase C)

- AMBIG strand-observability on the residual (silent-strand override, directional nascent prior from
  per-strand mature, neighbour carry-over).
- **Retiring `region_splice_gdna_frac`** — only after the Phase-C AMBIG recovery replaces it.
- The flagship AMBIG-dominated leak drop (the big number) — lands in Phase C.
