# Calibration — current status & resume point (2026-07-19)

**Purpose:** a snapshot of where the calibration effort stands, written to resume cleanly in a new session.
Companion: [`CALIBRATION_MASTER.md`](CALIBRATION_MASTER.md) (design + phased plan),
[`CALIBRATION_ARCHITECTURE.md`](CALIBRATION_ARCHITECTURE.md) (the count-zero-information invariant, authoritative).

---

## TL;DR — the one thing that matters

**The root cause of our recent struggles is upstream and now identified:** the **pass-0 solve does not subtract
mature RNA at single-strand exons**, so it never leaves clean gDNA. The mechanism that should do this — the
**"nascent factory," `ρ_nascent = ρ_RNA − ρ_mature`** — is a *documented, un-landed gap*
([`bp_solver.py:442-446`](../../src/rigel/calibration/bp_solver.py#L442-L446): *"an exon's ~95%-mature unspliced
payload leaks in as nascent, since we do not yet subtract mature from the RNA-total factor … until the nascent
factory land."*). Everything else we tried this session (a pass-0 background anchor, strand-deferral) was
**patching this defect downstream**, and has been **reverted**. The next work is to **finish the mature
subtraction through the message system** and re-measure.

**Evidence** (`scratchpad/ss_exon_pass0.py`): single-strand exon pass-0 f_g vs oracle correlates **0.90 when
stranded** (strand carries f_g) but **collapses to 0.62 unstranded** (regresses to the 0.5 reference — the
missing subtraction). The *arithmetic* spliced-derived gDNA already hits **corr 0.895 at ss=0.5** (memory
`unstranded_spliced_derived_rho`), so completing the subtraction should recover most of that gap.

---

## The problem (unchanged)

Deconvolve each genomic node's **unspliced** fragment mass into `(f_pos, f_neg, f_g)` — sense-RNA / antisense-RNA
/ gDNA. **Invariant:** a fragment *count* carries ZERO composition information (count-zero-info); the only
intrinsic composition signals are the **strand tilt** (Fisher info `∝ N(2κ−1)²`, **= 0 at κ=½ / unstranded**)
and the **spliced count** (gDNA cannot splice ⇒ spliced reads are pure RNA — the one composition-informative
count, and the key to the unstranded regime). The **hard case** is **unstranded + hybrid-capture**, where strand
is silent and gDNA-vs-RNA must come from the spliced signal + structure + messages + the prior.

## The pipeline (as it stands, after the revert)

```
 1. ENRICHMENT NPMLE (Role A)  fit DensityNPMLE on TOTAL density → σ²_transfer (message precision) ONLY.
 2. PASS-0 (PRIOR-FREE)        strand + Beta(½,½) reference + belief-free FB messages + σ²_transfer.
                               gdna_prior=None. ← reverted to this; NO pass-0 DNA prior.
 3. FIT gDNA HYPERPRIOR (RoleB)  DensityNPMLE on the DECONVOLVED gDNA (ρ_g=f_g·M/E) at SELECTED nodes,
                               with the ρ_bg background floor. (`_fit_gdna_hyperprior`.)
 4. REFIT                       re-solve with the hyperprior as the gDNA arm. (`calib_refit_iters`, default 1.)
 5. PROJECT                     per-region / per-boundary-side gDNA-RNA mass → per-locus Dirichlet prior.
```

The whole chain is bottlenecked at **step 2**: pass-0 is where the mature subtraction should happen, and it
doesn't (see TL;DR). A bad pass-0 ⇒ a contaminated hyperprior ⇒ downstream compensation.

---

## What LANDED this session (kept)

- **Additive-KDE representation for the gDNA hyperprior (Role B).** `DensityNPMLE.fit(additive=True)` +
  `_fit_gdna_hyperprior(additive=…)` + `config.gdna_prior_additive` (**default OFF, byte-identical**). Motivation
  ([`kde_vs_npmle_enriched_mode`](../../) memory; [`gdna_kde_restore_plan.md`](gdna_kde_restore_plan.md)): the EM
  Poisson-mixture NPMLE *silences* the capture-enriched mode (a competitive mixture + a heavy `n_regions`
  in-mixture background). The additive path (occupancy-weighted, fixed-bandwidth kernels + a weak 1-pseudo-obs
  floor) keeps it — the shipped v0.7.1 KDE's representation. Tests in `tests/calibration/test_npmle.py`
  (`test_additive_*`). **Status: valid, but its value depends on a clean pass-0 (step 2) — re-evaluate after the
  nascent factory.**
- **OI-1 resolved (the RNA-parsimony reference).** The `½·log(1−f_g)` Beta(½,½) reference is **correct and
  unchanged**; the shipped `−log(1−f_g)` was an improper f_g→1 *promoter*, not a parsimony (a red herring). The
  real pull toward RNA at unstranded nodes comes from the (mature-subtracted) RNA content, not the reference.
  Reviewer doc for a statistician: [`rna_parsimony_reviewer_question.md`](rna_parsimony_reviewer_question.md)
  (its open §5 — two-sided vs one-sided ρ_bg — is **moot** if the mature subtraction fixes pass-0).
- **Diagnosis tooling:** `scripts/debug/anchor_regression_dissect.py` (node-level regression dissection, by
  class/type, nascent-stratified), `scripts/debug/confident_fp_trace.py` (init→strand→messages→final trace),
  `scratchpad/ss_exon_pass0.py` (the pass-0 mature-subtraction measurement), `scratchpad/bench36.py`
  (32-scenario benchmark).

## What was EXPLORED and CLOSED / REVERTED

- **Pass-0 background DNA anchor — REVERTED (owner, 2026-07-19).** A weak two-sided ρ_bg Gaussian at pass-0 +
  strand-deferral. It fixed the zero-DNA FP but **crushed capture-enriched exons** (the anchor references the
  intergenic floor, the wrong reference for a captured exon) — a downstream patch for the missing mature
  subtraction. Pass-0 is now prior-free again (`_BackgroundAnchor` + `pass0_bg_anchor` flags removed).
- **Fitted RNA arm (`logP_r` / S4) — REJECTED** (full-solve A/B failed; improper promoter). Not an option.
- **Signature-class stratification of the gDNA fit — REJECTED** (owner directive). Not pursued.

---

## THE NEXT WORK — the nascent factory (mature subtraction)

**Goal:** make pass-0 subtract the spliced-measured mature RNA density from each node's RNA-total *before* the
nascent imputation, so the residual left for gDNA is clean. `ρ_nascent = ρ_RNA − ρ_mature`, intron-baselined.
The spliced density is already carried in the message MODE (`rho_pos += SPs/esp`, `bp_solver.py:547`) and credits
PRECISION (the I_spliced work); what's missing is treating the node's **own** mature-unspliced RNA as a
*measurement* to subtract, rather than lumping it into the fuzzy nascent imputation.

**Acceptance:** single-strand-exon pass-0 corr(f_g, oracle) **0.62 → ~0.90** on the unstranded conditions; then
re-fit the hyperprior on the (now clean) gDNA and re-check whether any pass-0 DNA prior / anchor is needed at
all (likely not).

**Owner's open choice (from the last exchange):** (1) prototype the nascent factory directly behind a flag, vs
(2) write a short design/derivation note first (the mature-vs-nascent split touches the RNA-total factor, the
boundary spliced-absorption, single-exon genes with no junctions, and the intron-baseline for nascent — subtle,
and the message layer has bitten us before). Recommendation was **(2) a tight design note, then implement**.

## Other open items (lower priority, after the nascent factory)

- **Additive-KDE default-on decision** — re-run the recovery/FP/stranded gates once pass-0 is clean; flip
  `gdna_prior_additive` on if it holds. Data-driven bandwidth (the fixed `npmle_bandwidth=0.15`) is a later knob;
  final value needs **real data**.
- **The message-absorption correctness item** ([`message_absorption_fix.md`](message_absorption_fix.md)) — the
  splice-junction RNA-imputation false-positive; related to the same message layer as the nascent factory,
  likely addressed together.
- **Real-data validation** (LBX0190 cfRNA, MO_3005 in `~/Downloads/rigel_runs`) — the sim under-represents the
  intergenic flood and has unrealistically abundant nascent; final tuning of the background weight, bandwidth,
  and intron-inclusion must be settled on real libraries.
- **Golden regen** — deferred to finalization (there is a pre-existing known-red test
  `test_gdna_sweep_zero_gdna_pin_and_monotone`, value ~0.68, unrelated to this work).

## Key memory notes (context)

`pass0_mature_subtraction_gap` (★★★ the root cause), `kde_vs_npmle_enriched_mode` (★★★ why KDE keeps the enriched
mode), `unstranded_spliced_derived_rho` (the 0.895 arithmetic result), `gdna_projection_snap_fix`,
`calibration_architecture_count_zero_info`, `nascent_rna_identifiability_intron_required`.
