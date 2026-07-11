# AMBIG / boundary deconvolution and the spliced-mass message — reference guide

**Status:** the density-prior application is **SOLVED and SHIPPED** (see the RESOLUTION banner below);
the remaining residuals (unstranded-with-gDNA, boundary class) still point at the spliced-mass anchor
and are scoped for a dedicated session. This document is the guide for that session.

---

> ## RESOLUTION (shipped): the generative two-density prior
>
> The density-prior integration challenge in §2.5 / §3 was solved by an ultracode design workflow and
> shipped in `bp_solver._kde_logprior`. The per-node posterior over f_g is
>
> ```
>   post(f_g) ∝ L_strand(f_g) · P_KDE(log(f_g·d_tot)) · 1/(1−f_g)
> ```
>
> — the strand likelihood landscape × the empirical gDNA-density KDE (evaluated with **real Gaussian
> tails**, `GdnaDensityPrior.logpdf_kernel`, NOT the clamped interpolation) × a **Jeffreys RNA density
> prior** `1/(1−f_g)`. The Jeffreys term is the crux: it removes the gDNA-priored/RNA-free asymmetry
> (§2.5) that caused both the cliff and the false-positives, with **no tuned constants** (Jeffreys
> exponent 1, native-log KDE coordinate). Applied in pass-2 on top of the weak stability floor; the
> per-node posterior is read out as the existing median (the floor + messages regularise it, so the
> harness's median-cliff does not appear in the full solver).
>
> **Validated:** contained gDNA calibration error −59,360 → −2,998 (95%); AMBIG-exon class −56,218 →
> −6,574 (88%); no gdna_none false-positive regression; end-to-end target net leak 235k (floor) → 171k
> (−27%) with nRNA hallucination halved (79k → 42k); 1084 tests pass, goldens regenerated. Rejected
> alternatives (proximity self-pin, M/E anchor, global-amount, raw-grid-multiply) are recorded in §5.
>
> **Still open (the spliced-mass session):** (1) unstranded-WITH-gDNA single-strand exons over-call
> (Jeffreys diverges as f_g→1 with no strand envelope; needs the motif-stranded spliced anchor, §7.1);
> (2) the boundary class regressed slightly (−20.9k → −35.7k) because the region-trained KDE penalises
> boundary crossing-densities in its valley — needs boundary-appropriate density handling. Both are the
> spliced-mass work below.

---
**Audience:** anyone, no prior code knowledge assumed.
**Scope:** the calibration stage only (per-region gDNA-vs-RNA deconvolution). The EM is untouched.
**Convention:** everything below is in **counts and densities (log-space)**, never fractions.
`f_g` (the gDNA fraction of a node) appears only as a *label* for oracle comparison — it is never a
working quantity or a prior target.

---

## 0. TL;DR

Calibration under-calls genomic-DNA (gDNA) badly on regions annotated with genes on **both strands**
(**AMBIG** nodes) and on **boundary** nodes — 95% of the calibration error on our benchmark is AMBIG
exons. The strand of the reads cannot separate "unstranded gDNA" from "two equal-and-opposite RNAs"
there, so these nodes depend on the *prior* + *messages*. Three things were found:

1. The **strand tilt** already provides the right *envelope* (it fixes the majority-strand RNA and caps
   gDNA at `2·min(N₊,N₋)`); the job of the prior is only to resolve the balanced **ambiguous core**.
2. The density prior, as applied, has a **cliff**: any weak-strand node whose total density is below a
   threshold is crushed to the intergenic (depleted) density mode, i.e. declared almost pure RNA. This
   is an artifact of an **asymmetric prior** (gDNA density is priored; RNA density is free) meeting a
   tall, low depleted KDE mode.
3. The evidence that *should* resolve these nodes — the **spliced (junction) fragment mass**, which is
   pure RNA and strand-resolved even in unstranded libraries — is propagated through the message system,
   but its handling has a known divergence and a live A/B toggle. **Getting the spliced-mass message
   right is the fix** and is the subject of the dedicated session.

Dead ends (do not revisit — see §5): the proximity self-pin prior, the M/E total-density anchor, and the
global-amount term.

---

## 1. The deconvolution problem and the node model (background)

A library mixes RNA fragments and gDNA-contamination fragments; they are often sequence-identical, so
individual fragments cannot be labelled. Calibration estimates, **per region**, the gDNA vs RNA
composition of the fragment mass, splitting RNA by strand: each region/boundary node carries
`(gDNA, RNA₊, RNA₋)` masses.

**The only intrinsic gDNA/RNA signal is strand** (`CALIBRATION_ARCHITECTURE.md`, the
count-zero-information principle): a fragment *count* carries no gDNA/RNA information; only the strand
*pattern* does, as overdispersed (Beta-Binomial) Fisher information. gDNA is **unstranded** (≈50/50
across strands); RNA is **stranded** (a + transcript reads "+" with probability κ).

**Node classes** come from the annotation signature:
- **POS / NEG** — single-strand; strand is informative; pass-1 solves these well without any prior.
- **AMBIG** — both strands annotated (e.g. overlapping opposite-strand genes); strand alone cannot
  separate gDNA from two-strand RNA.
- **NONE** — intergenic (structurally pure gDNA).

**Densities.** A node's **total density** is `d_tot = M/E` (fragment mass ÷ effective length). Its gDNA
density is `ρ_g = f_g·d_tot`; RNA density `ρ_r = (1−f_g)·d_tot`. The density prior lives in `log ρ_g`.

**The two-pass sweep** (`calibrate.calibrate` → `bp_solver.node_sweep`):
1. **Pass 1** — solve every node from: its **strand likelihood**, cross-node **imputation messages**,
   and a weak **global floor**. (Intended as "prior-free"; see §4 — it is not, quite.)
2. **Phase 2** — fit a nonparametric density prior `P(log ρ_g)` (a weighted KDE) on the *solved
   single-strand + structural* nodes (**AMBIG excluded → non-circular training**; ≥ `_MIN_KDE_TEACHERS`
   = 10 nodes, else skip).
3. **Pass 2** — re-solve every node with the density prior added.

**Spliced fragments live only at boundaries.** A region's interior has only unspliced fragments;
spliced (junction-spanning) fragments deposit at boundary nodes and are strand-resolved by the splice
**motif** — pure RNA, stranded *even in an unstranded library*.

---

## 2. Diagnosis (all numbers validated on `gdna300 / ss0.99 / capture-on`)

### 2.1 The error is one node class
Per-node oracle comparison (`scripts/debug/dissect_regions.py`), gDNA **counts**:

| class × type | n | oracle gDNA | production gDNA | error |
|---|---|---|---|---|
| **AMBIG × exon** | 92 | 198,429 | 142,211 | **−56,218** (95% of total) |
| POS × exon | 263 | 666,302 | 668,364 | +2,062 |
| NEG × exon | 294 | 739,205 | 736,186 | −3,019 |
| NONE × intergenic | 102 | 42,448 | 42,448 | 0 |

Total contained-gDNA error −59,360; boundary-side gDNA error a further −20,874.

### 2.2 The strand tilt is an envelope, and the oracle sits at its ceiling
With unspliced strand counts `N₊, N₋` (`N = N₊+N₋`), because gDNA is symmetric (`G/2` per strand):
`N₊ = G/2 + RNA₊`, `N₋ = G/2 + RNA₋`. Hence, as basic arithmetic (the Beta-Binomial softens it):
- **majority-strand excess `|N₊−N₋|` is unambiguously RNA** (gDNA cannot make an imbalance) — determined
  at high precision;
- **gDNA is capped: `G ≤ 2·min(N₊,N₋)`**;
- the balanced remainder `2·min(N₊,N₋)` is the **ambiguous core** = gDNA + minority-RNA — only this
  needs the prior.

Validated — the oracle gDNA tracks the tilt ceiling `2·min/N` in every balance bucket:

| balance (min/N) | n | tilt ceiling | oracle | production |
|---|---|---|---|---|
| 0.20–0.35 (skewed) | 6 | 0.548 | **0.538** | 0.457 |
| 0.45–0.51 (balanced) | 23 | 0.969 | **0.970** | 0.656 |

So the AMBIG exons are **gDNA-saturated up to the strand ceiling** (under capture the balanced high-density
reads really are unstranded gDNA), production **respects** the ceiling (never exceeds it) but sits far
**below** it, and **69% of the error is in the balanced bucket** — the hardest case, where the whole mass
is the ambiguous core and the prior must do all the work. These are the predicted
**nested opposite-strand transcript** loci (region 242: a 27 kb +strand intron with 7 −strand
transcripts nested inside; region 1083: 16 nested pairs) and their region interiors carry **zero spliced
fragments** (0/92) — nothing local disambiguates them.

### 2.3 The strand likelihood is a broad, multi-root curve — not a point
Reconstructing a node's strand posterior over `f_g` (`scripts/debug/ambig_posterior_shape.py`): balanced
nodes are **unimodal but heavily skewed with a long low tail** (median 0.87, mean 0.75, mass spread from
0→1); skewed nodes are **genuinely multi-modal** (3–5 modes). A single `(mode, precision)` summary is a
poor representation; the long tail is what lets a downward message/prior slide the estimate.

### 2.4 Root cause of the crush (measured on region 242, oracle 0.976)
```
strand likelihood alone   ≈ 0.87   accurate — balanced parsimoniously reads as gDNA
+ global floor             → 0.62   the floor drags it down
+ imputation messages      → 0.30   depleted-intron neighbours impute low gDNA density
+ pass-2 density prior      → 0.30   confirms the crash (proximity self-pin, see §5)
```
Two mechanisms:
- **The messages contaminate.** Ablation: removing the imputation messages lifts the AMBIG-exon
  aggregate 0.543 → 0.678. Under capture, gDNA density is exon-concentrated, but the gDNA-density
  imputation flows **depleted intron density into captured exons** (genomically adjacent, opposite
  capture status). The disagreement gate that should stop this is anchored on the *floor-crushed* local
  belief, so a low message "agrees" with the crush and keeps high precision.
- **The pass-2 prior cements it** (see §5).

### 2.5 The density-prior CLIFF (the crux open problem)
Fit the production KDE and sweep a node's total density under a **flat (weak) strand**
(`scripts/debug/` cliff probe). The KDE has a depleted mode at `log ρ_g = −4.39` (ρ≈0.012 — intergenic
gDNA) that is **4.5× taller** than the enriched mode at `log ρ_g = 2.78` (ρ≈16 — captured gDNA); 90% of
training nodes sit at the depleted mode.

```
 log d_tot   d_tot   solved log ρ_g   result
   0.0       1.0      −4.5            depleted  (crushed)
   1.0       2.7      −4.1            depleted  (crushed)
   2.0       7.4      −4.9            depleted  (crushed)   ← boundaries live here
   2.5      12.2      −4.4            depleted  (crushed)
   3.0      20.1      +2.77           enriched  (snaps)     ← hard cliff
```

**A hard cliff at `d_tot ≈ 12–20.`** Below it, any weak-strand node is crushed to the intergenic gDNA
density (≈ all RNA); above it, it snaps to enriched. This is exactly why the raw grid-multiply prior
*fixes AMBIG exons* (d_tot ≈ 15–30, above the cliff) but *crushes pure-gDNA boundaries* (d_tot ≈ 7,
below it): raw grid-multiply gave AMBIG-exon 0.543→0.612 but boundary error −20,874→**−78,504**.

**Why the cliff exists (the asymmetry).** gDNA density carries a prior `P(log ρ_g)`; RNA density carries
**none** (its 10,000× dynamic range makes any RNA level free). So for a flat-strand node the cheapest
explanation is always "dump the mass into free RNA and park gDNA at the tall, cheap depleted mode,"
reachable by lowering `f_g` until `f_g·d_tot` hits ρ≈0.012 — possible whenever `d_tot ≲ 20`. The argmax
(or median of the bimodal posterior) then pins there. **The prior is being used to invent RNA and pull
the density down, which is backwards.** This is the "brittle, cliff-like, argmax single-solution
pinning" — a KDE that looks correct as a plot but destroys a node when its density is mapped onto it.

---

## 3. The fix direction

### 3.1 Density prior: the raw grid-multiply (settled)
Apply the single-strand-trained prior by **multiplying it into the strand likelihood at every candidate
`f_g` and marginalising**:

```
log posterior(f_g) = ψ_strand(N₊,N₋ | f_g)  +  log P( log(f_g · d_tot) )
```

- **Non-circular:** trained on single-strand nodes, applied to AMBIG; the application never reads back
  the node's own solved `f_g`.
- **Robust to RNA dynamic range:** `d_tot` only converts each candidate `f_g` into its implied gDNA
  density at that grid point; it is *not* a fixed anchor. A highly-expressed node (`d_tot = 10⁴`) can
  reach a plausible gDNA density (ρ ≲ 16) only at tiny `f_g`, so the prior correctly forces it to RNA.
- **Strand precision is automatic and honest:** `ψ_strand` is the Beta-Binomial log-likelihood *in
  counts*, so its sharpness *is* the strand Fisher information. A confident stranded node has a peaked
  `ψ_strand` that dominates the sum → the prior is barely heard; a balanced/AMBIG node has a `ψ_strand`
  flat along the tilt-degenerate direction → the prior decides. No separate precision knob is needed for
  the prior (honest precision remains only for the inter-node *messages*).

**But the raw grid-multiply alone has the §2.5 cliff.** It is necessary, not sufficient. The cliff is
removed by the spliced-mass anchor (§3.2), which supplies the missing RNA evidence so the prior is no
longer guessing.

### 3.2 The real fix: spliced mass as the RNA anchor (dedicated session)
The evidence that breaks the flat-strand degeneracy is the **spliced (junction) fragment mass** — pure
RNA, strand-resolved by the splice motif, present at boundaries. Correctly propagated, it *measures* the
RNA at/adjacent to a node, so the gDNA is the residual and the density prior no longer has to guess (and
cannot crush):
- a boundary/region with real RNA gets its RNA **measured** from the spliced mass → gDNA is the residual;
- a pure-gDNA boundary has ~no spliced mass and a balanced unspliced crossing → gDNA, not crushed;
- crucially, spliced strand comes from the **motif**, so this works for **unstranded** libraries too,
  where it is the *only* RNA signal.

This is why the spliced-mass message is the fix and deserves its own careful session.

---

## 4. Current spliced-mass architecture — audit (starting point for the session)

All in `bp_solver.py`. Spliced mass is one-sided and strand-routed to its exon flank by the junction
motif (`build_node_geometry`, `spliced_pos_left/right`, `spliced_neg_left/right`; regions carry zero
spliced). Per-face RNA density (`node_densities`) already sums nascent + a one-sided mature/spliced
floor with the correct half-triangle eff-length `E[min²/2ℓ]` (`spliced_efflen_not_2x_…`).

The message propagation (`node_sweep._scan`, RNA-pos/RNA-neg branches) implements **three** spliced
terms:

```
ρ_message = (source nascent  fbp·M_src/E_r)
          + (source-face mature  SP_src/E_spl_src)     ADDED  only B→exon: MEASURES the exon's mature
          − (dest-face mature   SP_dst/E_spl_dst)      ABSORBED only exon→B: mature crosses as spliced,
                                                       so only NASCENT is imputed into the intron side
```

Mapping to the intended directions:
- **SJ boundary → region** (impute spliced boundary mass into the region's RNA): **present** as the
  ADDED source-face mature term. ✅
- **region → SJ boundary** (impute region unspliced RNA onto the junction by junction strand): the code
  instead **subtracts/absorbs** the dest-face mature so only nascent crosses. ⚠️ **This diverges from
  the stated intent** and is gated by a live A/B toggle `RIGEL_DBG_MAT` (`''|off|no_absorb|no_meas`) —
  a prime audit target.
- Precision: all RNA messages (mature MEASUREMENT and nascent IMPUTATION alike) use the
  disagreement-aware `σ²_edge` anchored on the dest's message-free local belief.

Known caveat already in the comments: junction-spanning spliced reads are only **partially captured**,
so `n_mat/E_spl` can *under*-estimate the exon's RNA — the reason the mature message was made
disagreement-gated rather than full-count-precision. This capture bias of the spliced eff-length is
directly relevant to getting the anchor's magnitude right.

---

## 5. Rejected alternatives (do NOT revisit — with reasons)

1. **Proximity self-pin prior** (added earlier this session; the current live `gdna_prior` path). It
   computes `x_obs = log(f_g·M/E)` from the node's **own current f_g** and pins the node there with
   precision `1/σ_node²`; for a high-count node `σ_node→0` so the precision explodes (6199 on region
   242), cementing the pass-1 crash. A prior conditioned on the datum, applied to that datum. **Remove.**
2. **M/E total-density anchor** (anchor the prior at the fixed `log(M/E)`). Fatal on real data: total
   density is RNA-dominated with >10,000× dynamic range, so two nodes with identical gDNA density but
   1000× different expression get opposite anchors. Only looked good because the simulator's RNA range is
   compressed.
3. **Global-amount term** (default a flat-strand node toward the library gDNA fraction). Not viable under
   hybrid capture: the library averages enriched and depleted regions, and the global under-estimates the
   **enriched** density by >30×. At most a hyper-hyperprior, never the per-node resolver.

The through-line: you cannot resolve a flat-strand node by a *fraction* prior or a *fixed-density*
anchor. You resolve it by (a) the strand envelope where strand exists, and (b) **measured RNA (spliced
mass)** where it does not.

---

## 6. Plan for the dedicated session

Prototype on the cached payload (`scripts/debug/dissect_regions.py` cache), measured node-by-node in
counts/log-density against the by-origin oracle. One node class at a time.

1. **Audit + verify the spliced-mass message** end to end: quantify, per boundary and per adjacent
   region, how much spliced (RNA) mass exists, how much is imputed, and compare the resulting RNA to the
   oracle RNA. Resolve the exon→B **absorb-vs-impute** divergence (`RIGEL_DBG_MAT`) against first
   principles and the oracle.
2. **Correct the spliced eff-length capture bias** so the mature MEASUREMENT is unbiased (`n_mat/E_spl`
   should equal the true exon RNA density under uniform capture).
3. **Confirm the anchor removes the cliff**: with RNA measured from spliced mass, re-run the §2.5 flat-
   strand sweep and the boundary node class — the crush should disappear without a global/M-E hack.
4. **Then** switch the density prior to the raw grid-multiply (§3.1) and re-check AMBIG exons, POS/NEG,
   and boundaries together.
5. **Make pass-1 truly prior-free** (move the floor out of pass-1) and verify no regression on zero-count
   nodes.
6. **Validate on stranded AND unstranded** conditions (§7 concern 1), then the full 16-condition
   benchmark; regenerate goldens only after net-flow is confirmed non-regressing.

### 6.1 Diagnostic tools (already built this session)
- `scripts/debug/dissect_regions.py` — per-region production-vs-oracle gDNA (splits BAM by origin).
- `scripts/debug/ambig_exon_probe.py` — per-node solve decomposition (strand / floor / messages / prior)
  with a monkeypatch that injects `_capture` into `node_sweep`.
- `scripts/debug/pass1_ablation.py` — ablate floor / messages / prior on the AMBIG-exon class.
- `scripts/debug/ambig_posterior_shape.py` — reconstruct a node's strand posterior over `f_g`.

### 6.2 Temp scaffolding to remove before shipping
Uncommitted diagnostic toggles from this investigation (all default-off, but must be cleaned):
`RIGEL_DBG_MAT`, `RIGEL_DBG_GLOBAL_OFF`, `RIGEL_DBG_MSG_OFF`, `RIGEL_DBG_FLOOR_OFF`, `RIGEL_DBG_RAW_KDE`,
`RIGEL_DBG_ME_ANCHOR` (all in `bp_solver.py`); the `_capture` debug additions in `node_sweep`;
`RIGEL_DBG_PRIOR_WEIGHT`, `RIGEL_NO_CAPTURE_EFFLEN` (`pipeline.py`); `RIGEL_DBG_GDNA_LLR`
(`estimator.py`); `RIGEL_DBG_GDNA_EFFLEN_SCALE`, `RIGEL_DBG_GDNA_EFFLEN_CONTAINED` and the
`c_sq_over_len`/`c_span` projection entries (`priors.py`).

---

## 7. Concerns / open questions

1. **Unstranded data is the real test.** The whole fix leans on: strand envelope where strand exists,
   spliced mass where it does not. In an unstranded library the strand envelope is empty everywhere, so
   the **spliced motif strand is the only RNA signal** — the spliced-mass message *must* be correct or
   unstranded+capture will over-call gDNA on expressed exons. Validate on `ss0.50` explicitly; never ship
   on stranded results alone.
2. **Genuinely-expressed overlapping opposite-strand genes.** A balanced, high-density AMBIG node that is
   *actually* two equally-expressed opposite-strand RNAs is intrinsically ambiguous from strand; only the
   spliced mass (two-strand spliced) can reveal it. Confirm the spliced anchor distinguishes this from
   unstranded gDNA. Test with an overlapping-gene stress scenario, and never test AMBIG without a robust
   single-strand substrate (`ambig_needs_robust_training_subset` memory).
3. **Spliced eff-length capture bias** (partial capture of junction-spanning reads) biases `n_mat/E_spl`
   low — the anchor's magnitude depends on fixing this (session step 2).
4. **exon→B absorb vs impute** — the current absorption may be right (mature leaves as spliced, only
   nascent crosses) or may need to become an imputation onto the junction as stated; resolve against the
   oracle, not intuition.
5. **KDE valley / bimodality** — the depleted mode is 4.5× taller and 90% of training nodes; even with
   the raw grid-multiply and the spliced anchor, confirm intermediate-density gDNA nodes are represented
   (or that the strand+spliced evidence makes the KDE shape moot for them).
6. **Benchmark is hard-label.** Use the per-node oracle dissection as the primary metric; the
   16-condition net-flow (`calibration-benchmark` skill) is secondary and partly blind to soft shifts.

---

## 8. Pointers
- Architecture / count-zero-info: `CALIBRATION_ARCHITECTURE.md`
- Spliced deposit + junction strand: `spliced_fragment_length_deposit.md`, `spliced_mature_nascent_message.md`
- Disagreement-aware message precision: `dispersion_aware_message_precision.md`
- Density prior history (superseded proximity design): `KDE_boundary_prior_review.md` (§10 — add a
  correction banner pointing here)
- Log-density solver: `log_density_1d_solver_design.md`
