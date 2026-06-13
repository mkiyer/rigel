# Calibration deconvolution — the authoritative plan (current state → production)

**Status:** SINGLE SOURCE OF TRUTH, 2026-06-13. Consolidates and **supersedes** the scattered
propagation/simplex planning docs (see §7). Where this doc and any older note disagree, this doc wins (and
`calibration_theory.md` remains the authoritative *theory* reference). Written after a full code+docs audit
and an outside-reviewer reconciliation.

The calibration stage deconvolves each genomic region's unspliced fragment mass into **gDNA vs RNA** (a
per-region gDNA fraction `f_g`) that becomes the per-locus EM prior. This doc states exactly where the code
is, the one target architecture we are building toward, every resolved design decision, the phased path,
and the consolidation actions.

> **§0 — rev-2 (2026-06-13, post-review refinement).** Two simplifications, now baked into §2–§5: **(a) the
> coupling flips to per-strand RNA-density propagation; gDNA is the residual and does not propagate** — so
> all capture-enrichment machinery (`γ_ij`, enrichment-class chains) is **dropped** (the gDNA *fraction* is
> capture-invariant; gDNA *density* — the only discontinuous quantity — is never propagated). **(b)
> `σ²_RNA_local` is SOLVED** by a splice-junction-pair `var~mean` LOESS (the RNA mirror of the gDNA curve),
> retiring the last open modelling piece. Boundaries are confirmed **mass-bearing nodes** (the accumulator
> gives them their own crossing mass), not edge factors.

---

## 1. WHERE WE ARE (audited 2026-06-13)

**The one method that actually runs** is a per-node **inverse-variance fusion** of two independent
estimators of the gDNA fraction:

```
f_g = w·g_strand + (1−w)·g_count ,   w = I/(I+β) ,   I = N·(2κ−1)²
  g_strand : Beta-Binomial strand posterior (capture-invariant; defined on strand-observable POS/NEG nodes)
  g_count  : splice-upgraded count gDNA-fraction (count magnitude; AMBIG/intergenic use this alone, w=0)
  β        : the count-trust penalty (= the old hard-coded I₀=10; now the explicit knob count_trust_beta)
```

This is the production path (`strand_deconv.deconv_regions` / `deconv_sides`, `calibrate.py:244-266`). It is
log-likelihood fusion: strand precision `I` vanishes at κ=½ (strand silent → count governs) and grows with
depth×discriminability (good deep strand → governs). `β` caps the count's precision so a confidently-biased
count cannot beat a good strand.

**Status of the "simplex" work:**
- `simplex_propagate.deconv_regions_simplex` (the `use_propagation=True` branch) is **byte-identical** to
  production at β=I₀ — it is the same fusion, β-named. It adds nothing yet; it is the entry point the sweep
  will extend.
- The per-node grid-MAP (`simplex.solve_node`) was tried as the combine and **regressed +8.7 pt net leak**
  (the strand posterior's overdispersion skew biases `f_g` low); the inverse-variance fusion was restored to
  close that gap (`5acb9c8`). **Any future grid solver must beat that regression bar.**
- The **bidirectional sweep is NOT wired** — `propagate_simplex`/`_rts_smooth` (a *scalar* Kalman/RTS on the
  gDNA density only) are built and order-independence-tested but unused. The first sweep attempt **smeared**
  (propagating the raw count density), which is why it is off.

**Dead / superseded code** (audit-confirmed, 0 production callers): `propagation.py` (the entire
count-cascade: `propagate`, `propagate_regions`, `_propagate_strand`, `_solve_ambig`, …); the scalar-RTS
sweep scaffold (`propagate_simplex`, `_rts_smooth`, `_smooth_density`, `_coupling_process_var`,
`_local_gdna_density`, `_solve_view`, `_side_allow`) is built-but-unwired and will be **replaced** by the
grid sum-product (§2), not extended.

**Reusable primitives** (keep): `simplex._simplex_lattice` (the 2-simplex triangular lattice),
`simplex._mixture_strand_loglik` (the exact 3-component gDNA/RNA₊/RNA₋ strand likelihood on the lattice),
`strand_likelihood.strand_loglik`, `strand_deconv.strand_posterior_gdna_frac`,
`density_model.density_variance_curve` (the `var~mean` LOESS), the `substrate` 3-view layout,
`region_arrays`/`region_adjacency` (the genomic chain), `priors._transport_boundary_flux`.

**Bottom line:** there is effectively **one per-node method** (the fusion) duplicated behind a flag, plus
two dead scaffolds. The real, unbuilt innovation is the **spatial sweep** that resolves AMBIG/thin nodes.

---

## 2. THE TARGET ARCHITECTURE (reconciled — reviewer + our docs + our corrections)

**Belief propagation on a per-locus chain of REGION nodes.** A locus is a 1-D genomic chain → a tree → BP
is exact in two sweeps (forward, backward), order-independent.

- **State:** the 3-term pie `θ_i = (f_rna₊, f_rna₋, f_g)_i` on the 2-simplex (mass-conserving `Σ=1`;
  zero-injection safe). The calibration output is the marginal `E[f_g]`.
- **Message representation: a `(P,)` log-likelihood vector over the triangular lattice** (resolution K≈20 ⇒
  P=(K+1)(K+2)/2=231 — **not** the production `n_grid=200` lattice). This handles the simplex edges and the
  one-sided spliced inequality **exactly** (loglik = −∞ off-support). The **scalar inverse-variance message
  is rejected for the sweep** — it cannot carry the inequality or the simplex truncation, and post-hoc
  clipping re-introduces order-dependence (the reviewer is right; our own theory doc said the same).
- **Local evidence `ψ_i`** (a `(P,)` loglik) **= `solve_node` MINUS the count term**: the 3-component strand
  mixture (`_mixture_strand_loglik`) + the **sided spliced lower bound** (`f_c·U ≥ S_c` → −∞ on violation,
  fed from `substrate.contained.n_spliced_sense/antisense`) + a weak gDNA prior. The count is **not** local
  evidence — it arrives via the coupling (below). Seeds emerge from `ψ_i`: an intergenic seed is a sharp
  spike at the `f_g=1` vertex (no special-case code).
- **Boundaries are first-class MASS-BEARING NODES (not edge factors).** In our accumulator a boundary
  *owns* the fragment mass that overlaps it (it is not added to the regions), so the chain interleaves
  region–boundary–region; every node deconvolves its own mass into a pie. (The reviewer's "collapse
  boundaries into `φ_ij`" is **rejected** — it would silently drop the boundary mass + the flux transport
  that consumes it.) A boundary's **spliced flux** is the RNA-propagation signal; its **unspliced crossing**
  is gDNA+nascent.
- **Coupling = per-strand RNA-density propagation (PRIMARY); gDNA is the residual and does NOT propagate.**
  This is the rev-2 correction (see §0): the flip from "gDNA carries far + RNA local" to "RNA propagates,
  gDNA = `1−f₊−f₋`".
  - **`Q_RNA`** (per strand): finite where the RNA density is continuous — *within a contiguous same-strand
    exon stretch* (e.g. the non-overlap part of an exon → its AMBIG-overlap part: same exon, same capture
    enrichment, uniform RNA) — `∞` (decouple) at exon↔intron transitions (RNA density jumps mature→nascent)
    and where a strand is silent. `σ²_RNA_local` from the **splice-junction-pair `var~mean`** (§3 Q2).
    The AMBIG node inherits `f₊` from its `+`-neighbour and `f₋` from its `−`-neighbour → `f_g` = residual.
  - **gDNA does not propagate.** It is resolved per node by the strand (capture-invariant) or as the
    RNA-residual (AMBIG); junction-less nodes (single-exon / intergenic) fall back to the count + a **global
    off-target gDNA-density baseline** prior (the clean deep-intron/intergenic anchors fit this baseline).
- **Capture needs NO special handling — no enrichment classes, no `γ_ij`.** The gDNA *fraction* is
  capture-invariant (probe enrichment scales gDNA and RNA on an exon by the same factor → it cancels in
  `f_g`); the strand fraction is capture-invariant; the RNA propagation follows same-exon/same-transcript
  edges (same enrichment by construction). gDNA **density** (the only capture-discontinuous quantity) is
  never propagated, so the discontinuity never bites. `γ_ij` was unbuildable anyway (no per-region
  enrichment exists); class-chains were a flawed stand-in (they sever exon↔intron edges that, on inspection,
  carry no usable cross-enrichment gDNA signal). Both are dropped.
- **Inference:** forward then backward grid sum-product — `m_{i→j} = logsumexp(ψ_i + m_{prev→i} + M_ij, axis
  over θ_i)`, `M_ij` the `(P,P)` transition `log φ_ij`; belief `b_i = ψ_i + m_fwd + m_bwd`; `softmax →
  E[f_g] → CalibrationResult`. Exact in 2 sweeps; **per-locus chunked + parallel** (mandatory at scale, §5).
- **Boundary flux transport** (`priors._transport_boundary_flux`) stays **post-sweep, decoupled** (it
  translocates mass across junctions; the sweep only re-proportions a node's pie).

---

## 3. RESOLVED DESIGN DECISIONS (issues A–H, Q1–Q3)

| # | decision |
|---|---|
| **A** state vs coupling | State = **fractions** (pie, mass-conserving). The **gDNA fraction is capture-invariant** (enrichment cancels in the ratio), so gDNA needs **no enrichment normalization and is not propagated**. `γ_ij` and per-enrichment-class chains are both **dropped** (rev-2, §0/§2). |
| **B** boundary 1 node or 2 | **Boundaries are mass-bearing NODES** (the accumulator gives them their own crossing mass) — region–boundary–region chain. The reviewer's edge-factor collapse is **rejected** (it drops the boundary mass). |
| **C** spliced flux ownership | The **spliced flux is the RNA-propagation signal** + a one-sided RNA lower bound on its node; sided (sense/antisense), never double-counted. |
| **D / Q1** message rep | **2-simplex grid `(P,)` message + `(P,P)` transitions** (K≈20). Scalar/Gaussian rejected for the sweep; defer Gaussian/EP to the perf phase. |
| **F** var~mean circularity | **1-pass** fit from anchor seeds (intergenic + deep single-strand) pre-sweep (non-circular); optional outer IRLS/EM refit. |
| **G** boundary transport | **Decoupled, post-sweep** (`_transport_boundary_flux`). |
| **H** seeds | **No special-case code** — intergenic seed = a sharp `ψ_i` delta at the `f_g=1` vertex; "signal dies into a seed" emerges from the math. |
| **Q2** RNA variance model | **SOLVED (rev-2):** `σ²_RNA_local` from a **splice-junction-pair `var~mean` LOESS** — a region flanked by two same-strand splice boundaries with nonzero spliced counts is two measurements of the same RNA (`¼(d_L−d_R)²`), reusing `density_variance_curve`. Works at `ss=0.5` (splice is motif-stranded). Start with adjacent pairs; defer transcript-set dedup. |
| **Q3** per-strand chains | RNA propagates along same-strand genomic edges, `Q_RNA=∞` at exon↔intron RNA-discontinuities + silent-strand transitions; the full splice-graph (intron-skipping exon→exon) is a deferred extension. |
| count trust `β` | The integration lever; the explicit successor to `I₀`. Phased: single (now) → 2-level by count-observability → continuous → derived. |

---

## 4. THE PHASED PATH (current → production)

Each phase is gated on the **net-flow benchmark** (the primary metric), reference OFF baseline: capture-off
~1%, capture-on ss0.99 ~8.7%, ss0.5 ~20.3%. A phase ships only if it does **not regress** (and ideally
improves) — and **must beat the prior grid-MAP's +8.7 pt regression**.

- **Phase 1 — per-node fusion (DONE).** `f_g = w·g_strand+(1−w)·g_count`, `w=I/(I+β)`. Byte-identical to
  production; `β` explicit. This is the local resolve the sweep builds on; the no-propagation baseline.
- **Phase 2a — the splice-junction RNA `var~mean` model (NEXT, standalone).** Per the Q2 construction: find
  regions flanked by two same-strand splice-junction boundaries with nonzero spliced counts, take
  `¼(d_L−d_R)²` at `½(d_L+d_R)`, fit the LOESS via `density_variance_curve`. Validate the fit on the toy +
  flagship (does it produce a sensible `σ²_RNA(μ)`?). No regression risk (standalone).
- **Phase 2b — the RNA-propagation grid sum-product sweep (the core fix).** The `(P,P)` grid sum-product
  over per-locus region–boundary chains: `ψ_i` = `solve_node` count-stripped (3-component strand + sided
  spliced bound + weak gDNA prior); edge `φ_ij` = **per-strand RNA-density continuity** (`Q_RNA` from 2a,
  finite within same-strand exon stretches, `∞` at RNA-discontinuities/silent strands); **gDNA = residual**.
  Resolves AMBIG (`f₊` from the `+`-neighbour, `f₋` from the `−`-neighbour) and benefits `ss=0.5` (splice is
  motif-stranded). Validate **≥ parity and ideally an AMBIG/capture-on improvement** on the 16-condition net
  flow + flagship (must beat the prior grid-MAP's +8.7 pt regression).
- **Phase 3 — refinements.** Full splice-graph RNA propagation (intron-skipping exon→exon); transcript-set
  dedup for the variance model; `β` sophistication (2-level → continuous → derived); a global off-target
  gDNA baseline prior for junction-less nodes.
- **Phase 4 — productionize + teardown + perf.** Flip the default; **delete the old path + dead code**
  (collapse to one path); per-locus parallelism; Gaussian/EP messages only if perf demands.
- **β sophistication (parallel track):** single → 2-level (count-observable introns MAE 0.005 vs imputed
  exons 0.47) → continuous → derived; retires the magic number.

---

## 5. RISKS / OPEN MODELLING (read before building)

1. **`σ²_RNA_local`: SOLVED in principle** (Q2 splice-pair `var~mean`), but the *fit quality* is the
   Phase-2a validation gate — sparse same-strand splice-pairs, and the multi-transcript "count each
   boundary-set once" dedup (deferred). If the curve is unstable, fall back to a Poisson/constant `Q_RNA`.
2. **RNA-propagation soundness (the main modelling risk now).** RNA density is uniform *within an exon* but
   jumps at exon↔intron and across genes; the `Q_RNA=∞` gating (RNA-discontinuity + silent-strand) must be
   correct or RNA smears. The validated win is the *within-exon* AMBIG case (overlap part ← non-overlap
   part); intron-skipping exon→exon is the deferred splice-graph extension. Measure on the toy AMBIG sweep.
3. **Perf: the `(P,P)` sum-product cannot run whole-genome as one chain.** ~1e6 regions × P²(≈53k) is
   ~5e10 flops/sweep and a single `(n,P,P)` tensor is unmaterializable. **Per-locus chunking is mandatory**
   (`region_adjacency` + the per-reference run skeleton; loci are <~100 nodes, embarrassingly parallel). Use
   **K≈20 (P=231)** for messages — *not* the production `n_grid=200` (P=20301).
4. **The +8.7 pt regression bar** — the prior grid-MAP-as-combine failed it (strand-posterior overdispersion
   skew). The grid sum-product must verifiably **beat production net flow**, not just match pool fraction
   (the earlier "parity" that was a uniform leak increase).
5. **FL-unit consistency** — `Q` (from `var~mean`) and the state `ρ=f·U/L` must use the **same per-node L**;
   boundary densities use `1/fl_mean`, contained uses `/region_eff_len` — correct today but fragile in a
   rewrite.
6. **Boundary nodes own mass** — the region–boundary–region chain must conserve mass (boundary crossing
   mass deconvolved, not dropped); the flux transport still consumes the boundary gDNA post-sweep.
7. **gDNA = residual** assumes `f₊+f₋` are well-estimated; a junction-less node (single-exon / no spliced)
   has no RNA propagation → falls back to per-node strand + count + the global off-target gDNA baseline (the
   honest "unknown" floored by the baseline). Confirm this degrades gracefully (no false gDNA injection).

---

## 6. CODE CONSOLIDATION ACTIONS

- **DELETE NOW (safe, audit-confirmed 0 production callers):** `propagation.py` (the count-cascade) + its
  test `tests/calibration/test_propagation.py` + the debug script `scripts/debug/phaseC_validate_propagation.py`.
- **KEEP:** the production fusion (`deconv_regions`/`deconv_sides`); the simplex primitives
  (`_simplex_lattice`, `_mixture_strand_loglik`, `solve_node` as the `ψ_i` compiler — call with
  `count_precision=0` to strip the count); `density_variance_curve`; the strand machinery; the substrate;
  `priors`/`derive`.
- **REPLACE (Phase 2):** the scalar-RTS sweep (`propagate_simplex`, `_rts_smooth`, `_smooth_density`,
  `_local_gdna_density`, `_solve_view`, `_side_allow`) → the grid sum-product. Keep `_coupling_process_var`
  (the `var~mean` → `Q_g`).
- **COLLAPSE (Phase 4 teardown):** one deconvolution path; remove `use_propagation` + the redundant
  byte-identical `deconv_regions_simplex`; `solve_node`'s grid-MAP combine stays retired.
- **CONFIG:** `count_trust_beta` (= `β`, documented placeholder, successor to `gdna_strand_info_scale`=`I₀`).

---

## 7. DOCS CONSOLIDATION

**This file is the single authoritative plan.** It **supersedes** (→ archived to `archive/`):
`simplex_roadmap.md`, `count_trust_design.md`, `signal_integration_investigation.md`,
`propagation_message_passing.md`, `propagation_simplex_plan.md`, `propagation_sweep_plan.md`,
`propagation_increment4_plan.md`, `propagation_implementation_plan.md`, `deconvolution_roadmap.md`,
`sequential_calibration_redesign.md`, `decoupled_calibration_design.md`, `remaining_phases.md`,
`phaseA_mature_imputation_plan.md`, `phaseB_mature_subtraction_plan.md`, `phaseC_ambig_resolution_plan.md`,
`phase0_phase1_implementation_plan.md`, `phase2_design.md`, `phase3_implementation_plan.md`,
`redesign_phase1_plan.md`, `redesign_phase2_plan.md`, `redesign_phase3_plan.md`, `strand_first_plan.md`,
`count_channel_capture_design.md`, `count_local_dispersion_design.md`, `count_mean_bias_design.md`,
`count_posterior_design.md`, `density_phase1_local_imputation_design.md`,
`density_phase2_dna_fraction_design.md`.

**KEPT as theory/reference** (current, not superseded):
- `calibration_theory.md` — authoritative theory (the acyclic pipeline + the per-node model). Lightly note
  this plan is the forward spec.
- `deconvolution_generative_model.md` — the three-species (gDNA/nascent/mature) generative theory.
- `deconvolution_implementation.md` — the boundaries-as-atoms blueprint (background for §2).
- `strand_deconvolution_explained.md` — the ground-up strand-posterior explainer.
- `calibration_theory.md` + the FL/capture diagnostics (`fl_consistency_diagnostic.md`,
  `capture_effective_length_design.md`, `accumulator_fragment_span_redesign.md`) — accuracy records.
- `CALIBRATION_TODO.md` — the live tracker. `my_notes.md` — user design notes.

The `archive/` docs remain in git for history; this plan is the truth.
