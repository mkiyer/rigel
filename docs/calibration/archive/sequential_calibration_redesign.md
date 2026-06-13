# Sequential strand→count calibration — phased implementation plan

**Status:** authoritative redesign plan, 2026-06-11. Supersedes the decoupled *blend* (`w·g_strand +
(1−w)·g_count` on raw counts) and the `phase3_implementation_plan.md` direct-3-term. Grounded in the
strand-prior analysis (`strand_deconvolution_explained.md`) and validated against the AMBIG toy scenario
(`scripts/debug/phase3_ambig_scenario.py`). **Most of the implementation is intact — this is a
restructure, not a rewrite.**

## 0. The idea in one paragraph

The latent quantity is a per-node **gDNA fraction**, observed through two channels. The **strand module**
deconvolves each node's *unspliced* fragments into gDNA vs unspliced-RNA, **standalone**, and — crucially
— emits the **information** (precision) of that split. The **count module** consumes the *strand-cleaned*
counts: it builds a gDNA **density field** (the spatial prior), imputes it where strand is silent, and
combines the unspliced split with the deterministic **spliced** RNA (the 3-term). The strand module can
stand alone precisely *because* it communicates its uncertainty; the count module *applies* that precision.
This fixes the red flag — the count module currently imputes from **raw** counts, so nascent leaks into the
gDNA density (the AMBIG node reads 0.120 gDNA at 0 gDNA + nascent; its single-strand neighbours, which the
strand module *did* clean, read ~0).

## 1. The architecture (end state)

```
fit:   FL models · strand balance (ss) · gDNA/RNA strand overdispersions          [INTACT]

STRAND MODULE  (standalone, per node = each region-contained view + each boundary side):
   BB likelihood of the gDNA fraction g from (pos, neg | ss, overdispersion).
   - point split  g_q   = the likelihood read at the FP-quantile q   [the knob lives HERE — benefit #1]
   - information  I     = N·(2·ss−1)²   (the per-fragment Fisher info × N; → 0 at ss=½ and at AMBIG)
   - emit per node: unspliced_gdna = g_q·N,  unspliced_rna = (1−g_q)·N,  I.        [3-term inputs ready]

COUNT MODULE  (on the CLEANED counts):
   - DENSITY FIELD = the gDNA density prior, built from the strand-cleaned gDNA counts at OBSERVABLE
     nodes, precision-weighted by I; imputed across regions (runfill); eff-length / FL-consistency
     corrected.  Where I=0 everywhere (unstranded) it falls back to the RAW density (today's behaviour).
   - 3-TERM per node: gDNA fraction = unspliced_gdna / (unspliced_gdna + unspliced_rna + spliced),    [#2]
     density→count via density_frac_to_count_frac.
   - POSTERIOR per node = field-PRIOR ⊗ strand-LIKELIHOOD, precision-weighted:
        g = (g_field·P_field + g_q·I) / (P_field + I)
     → field governs where strand is silent (AMBIG/unstranded); strand governs where confident.

derive · priors · per-locus EM                                                     [INTACT]
```

### 1.1 Why the precision is the *likelihood* information, not the posterior variance

The Jeffreys `Beta(½,½)` prior gives the strand module a **proper** posterior even at `ss=½` — variance
`0.125`, *finite*. If the count module weighted nodes by that posterior precision, a fully **uninformative**
node (`ss=½`) would carry weight `1/0.125 = 8` and inject its prior-driven median (½) into the field. Wrong.

The quantity that correctly → 0 when the data say nothing is the **likelihood Fisher information**
`I = N·(2·ss−1)²` (derived in `strand_deconvolution_explained.md` §2,6). So the strand module emits a
**likelihood** `(g_q, I)`, and the **density field is the single prior**. This is what resolves the
circularity that stalled the old design: there is exactly one prior (the spatial density field), the strand
contributes a likelihood, and `posterior = field-prior ⊗ strand-likelihood` is a clean Bayesian update with
no EB-in-the-blind and no double-counted prior. At `ss=½`, `I=0` ⇒ the field/prior governs; at high
specificity, `I` is large ⇒ strand governs. (The Jeffreys prior is retained only to *regularise the
standalone split* `g_q` — keep it in [0,1], tame the `1/(½−ss)` blow-up — not to set the handoff weight.)

### 1.2 The combine weight — the one real risk (capture bias)

`posterior = field-prior ⊗ strand-likelihood` by precision is the *optimal* (BLUE) update **iff both are
unbiased**. The strand likelihood is unbiased. The density **field is capture-biased** (gDNA depleted at
exon edges vs interior — the junction-vs-body residual `~0.17` flagged in `fl_consistency_diagnostic.md`),
so naive precision-weighting can *launder* that bias into the unbiased strand at low/moderate SS — the exact
failure the old `w=(2ss−1)²` blend was built to avoid (and that the Phase-2 Fisher-blend diagnostic
confirmed). **Decision (validate, don't assume):** start with the principled precision combine **and**
debias the field as far as the eff-length machinery allows; if the suites show laundering at moderate SS,
fall back to capping the field's effective precision (a bias-robust ceiling, e.g. the `(2ss−1)²` effect-size
weight as a *cap* on `I/(P_field+I)`). The toy scenario + the gDNA/nascent suites are the arbiter — §6.

### 1.3 What this buys (the user's two benefits, realised)

1. **The FP-quantile becomes the strand deconvolution's own knob.** `q` reads the *strand* likelihood
   (`g_q`), the allowable false-positive rate of the **strand** split — clean and interpretable, no longer
   a fuzzy shift on a strand+count *center*. It is uncertainty-aware (wide likelihood ⇒ larger shift) and
   bites only where strand has information (`I>0`); AMBIG/unstranded FP is the field's business (the floor).
2. **Strand cleans first → count gets a pre-cleaned 3-term.** The count module receives, per region *and*
   boundary side, `(unspliced_gdna, unspliced_rna)` from strand, ready to impute and combine with the
   deterministic `spliced` RNA — the 3-term is then a *consumer* of the strand output, not a strand-split
   re-derived inside `splice_junction`.

## 2. Intact vs restructured (most code survives)

| Component | Fate |
|---|---|
| FL models, strand balance (`ss`), gDNA/RNA overdispersion fits | **intact** |
| `strand_posterior_gdna_frac` (BB posterior median+var, vectorised) | **extend** → emit `g_q` (quantile) + `I = N·(2ss−1)²` |
| `_grid_posterior_median` | **generalise** → `_grid_posterior_quantile(q)` |
| `runfill_bidirectional`, `same_ref_left_right` | **intact** (the field imputation) |
| `density_frac_to_count_frac` (FL-consistency) | **intact** (the density→count step of the 3-term) |
| `boundary_gdna_fraction` (3-term arithmetic) | **intact** (the count module's 3-term core) |
| `node_gdna_density` | **restructure** → density field from *cleaned* counts, precision-weighted by `I` |
| `deconv_regions` / `deconv_sides` / `_deconv_per_node` (the blend) | **split** into (a) strand deconvolve, (b) field-⊗-likelihood combine |
| `region_splice_gdna_frac` (the direct strand-split-inside) | **simplify** → consume the strand module's split |
| FP-quantile (`gdna_deconv_quantile`) | **move** from the blend to the strand likelihood |
| `derive`, `priors`, `result`, EM | **intact** |
| AMBIG toy scenario harness | **intact** → promoted to a locked test in Phase 4 |

## 3. Phases

Built **alongside** the live pipeline so each phase is independently testable; the switch happens in
Phase 3. Golden stays green until then.

### Phase 1 — Strand module: the standalone likelihood emitter
**Goal.** A pure strand deconvolver producing, per node, `(g_q, unspliced_gdna, unspliced_rna, I)`.
**Changes.**
- `strand_deconv.strand_posterior_gdna_frac` → also return `I = (sense+antisense)·(2·ss−1)²` and read the
  posterior at the **FP-quantile** `q` (generalise `_grid_posterior_median` → `_grid_posterior_quantile`).
- New `strand_deconvolve(substrate, region_arrays, ss, od_gdna, od_rna, q, n_grid)` → per-region and
  per-boundary-side `(g_q, unspliced_gdna, unspliced_rna, I)`. Strand-unobservable (AMBIG) and `ss=½`
  nodes ⇒ `I=0`, no split (passed through whole).
**Tests.** confident node splits; AMBIG/`ss=½` ⇒ `I=0`; `g_q` monotone in `q`; total invariant
(`gdna+rna = N`); reproduces `strand_posterior_gdna_frac` median at `q=½`.
**Exit.** strand module standalone + unit-tested; live pipeline unchanged.

### Phase 2 — Count module: density field on cleaned counts (the red-flag fix)
**Goal.** The gDNA density field is built from strand-cleaned gDNA counts, precision-weighted by `I`.
**Changes.**
- `node_gdna_density` → take the strand module's per-node `unspliced_gdna` + `I` (not raw counts); the
  observable-node density uses the cleaned gDNA count; the imputation precision-weights anchors by `I`;
  the unstranded fallback (`I=0` everywhere) reproduces today's raw-count density.
**Tests / validation.** the AMBIG toy scenario, calibration-level: AMBIG region gDNA fraction at
gDNA=0+nascent drops **0.120 → ~0** and at gDNA=100+nascent **0.343 → ~0.297**; the single-strand
neighbours stay correct; unstranded conditions unchanged.
**Exit.** the field is strand-cleaned; the AMBIG fix demonstrated in isolation.

### Phase 3 — Sequential wiring + 3-term + combine (the switch)
**Goal.** `calibrate.py` runs the new flow; the blend is gone.
**Changes.**
- `calibrate.py`: `fit → strand_deconvolve → count field (cleaned) → 3-term → posterior (field ⊗
  likelihood, precision-weighted, §1.2) → derive`.
- `region_splice_gdna_frac` → consume the strand module's `(unspliced_gdna, unspliced_rna)` for the 3-term
  (drop the strand-split-inside); keep `density_frac_to_count_frac` + `runfill` carry-over.
- Implement the §1.2 combine; settle precision vs bias-robust cap by measurement.
- Retire the old `_deconv_per_node` blend and the direct-3-term.
**Validation.** toy scenario across the full gDNA×nascent sweep (AMBIG solved, no harm at gDNA=0/no-nascent);
then the suites (§6). Regenerate golden (intended change).
**Exit.** new flow live; toy + suites pass.

### Phase 4 — Lock + cleanup
**Goal.** Make it permanent and tidy.
**Changes.** promote the AMBIG scenario to `tests/scenarios/test_phase3_ambig_carryover.py` (+ a small
golden over the sweep); full 3-pool suite validation isolating FL-fix vs redesign; delete dead blend code;
update `calibration_theory.md` + `decoupled_calibration_design.md` (note the evolution to sequential).
**Exit.** locked tests, green suite, docs current.

## 4. The uncertainty math (precise)

Per node `i` with unspliced (`pos_i`, `neg_i`), `N_i = pos_i+neg_i`, spliced `S_i`:
- **Likelihood** of the gDNA fraction `g`: BB(`pos_i`; `N_i`, `p(g)=ss+(½−ss)g`, overdispersion). Read the
  point `g_q,i` at quantile `q`; information `I_i = N_i·(2·ss−1)²` (0 at `ss=½`/AMBIG).
- **Cleaned counts:** `gdna_i = g_q,i·N_i`, `rna_unspl_i = (1−g_q,i)·N_i` (total invariant).
- **Density field (prior):** observable-node gDNA density `ρ_i = gdna_i / eff_len_i`; aggregate/impute with
  weights `I_i` (precision) → per-region `(g_field, P_field)` (P_field from propagated `I` + imputation
  spread; reuses the Phase-1 variance machinery). Unstranded ⇒ all `I_i=0` ⇒ raw-density fallback.
- **3-term fraction:** `g_3 = (gdna/eff_g) / (gdna/eff_g + (rna_unspl+S)/eff_r)`, then
  `density_frac_to_count_frac` to the region count fraction.
- **Posterior:** `g = (g_field·P_field + g_q·I)/(P_field + I)`; `P = P_field + I` (subject to §1.2).

## 5. Acyclicity

`ss/overdispersion fit → strand deconvolve → count field → 3-term → combine → derive`. Forward only: the
strand module needs no count input (Jeffreys-regularised); the field is built from strand-confident
observable nodes and is the prior for the rest; the combine consumes both. Assert no field→strand feedback.

## 6. Validation plan

- **Toy AMBIG scenario** (`phase3_ambig_scenario.py`) — the primary instrument: AMBIG gDNA fraction vs
  oracle across gDNA ∈ {0,20,100} × nascent ∈ {0,30}; require est→~0 at gDNA=0+nascent, est→oracle
  elsewhere, and **no regression** at gDNA=0/no-nascent.
- **gDNA suites** (`gdna_benchmark_5mb`, `gdna_shortfl_5mb`) — gDNA-vs-RNA separation, FP rate; **isolate**
  the FL-fix and the redesign as separate A/Bs (the lesson from the conflated quick_3to1 A/B).
- **Nascent + antisense** (`quick_3to1_5mb`) — 3-pool net flow; the AMBIG+nascent regions should improve,
  the rest stay flat (the §1.2 combine decision rides on this). Always rerun the full suite after any
  gDNA-shifting change before trusting a downstream A/B.

## 7. Open questions / risks

1. **The combine weight (§1.2)** — pure precision vs bias-robust cap. The single biggest decision; measured,
   not assumed.
2. **Field bias** — the junction-vs-body capture residual; how much eff-length correction is enough before
   precision-weighting is safe.
3. **`P_field` estimation** — propagating `I` + imputation spread into a usable field precision (reuse the
   Phase-1 LOESS/Poisson variance, now on cleaned counts).
4. **Boundary-side AMBIG carry geometry** — does the per-strand split need transcript-extent carry beyond
   the run-fill (design §5.1), or does the cleaned density field already deliver it? The toy scenario tells
   us.
