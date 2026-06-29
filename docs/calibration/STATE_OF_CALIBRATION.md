# State of the Calibration Algorithm

**Last updated:** 2026-06-28. Branch `calib-spliced-junction-strand` (uncommitted). This is a living
reference for where the gDNA-vs-RNA calibration stands — the architecture, the two solver backends, the
log-density overhaul, the benchmark, and the open problem with its root cause + fix direction.

---

## 1. What calibration does

Calibration deconvolves each genomic **node's unspliced fragment mass** into three components —
**sense-RNA, antisense-RNA, and gDNA** — and fits the library hyperparameters that the per-locus EM
consumes as a prior. It runs over a **bipartite region↔boundary node chain** by belief propagation
(forward-backward). Each node sees three evidence sources (`CALIBRATION_ARCHITECTURE.md`, the
count-zero-info principle):

1. **Strand likelihood** — the only INTRINSIC signal. A Beta-Binomial tilt of the per-strand counts;
   its precision is the Fisher information `N·(2κ−1)²`. At κ=½ (unstranded) it carries *zero* gDNA/RNA
   information — the count alone cannot separate balanced gDNA from balanced RNA.
2. **Cross-node imputation (messages)** — neighbours impute each other's component **densities** along
   the chain (gDNA flows genomically; per-strand RNA only where that strand is continuous).
3. **Global / enrichment prior** — the population gDNA density `ρ_global`, shrunk toward an
   enrichment-aware transfer `ê(z)` on capture-enriched exons.

> **Counts and densities, not fractions.** A node's *own* solve deconvolves its counts into
> (RNA₊, RNA₋, gDNA) **counts**; the composition fraction is only that solve's internal latent.
> **All cross-node communication is in DENSITY** (ρ = count / eff-length). Fractions are *not*
> comparable across nodes — an intron can be f_g≈1 at trivial density, an mRNA-rich exon f_g≈0 at high
> gDNA density — so they must never appear in message passing.

---

## 2. Two solver backends

`CalibrationConfig.calibration_solver` (env-overridable via `RIGEL_CALIBRATION_SOLVER`):

| | **`lattice`** (default, shipped) | **`logodds`** (the overhaul) |
|---|---|---|
| per-node latent | 2-simplex grid `(f₊,f₋,f_g)`, `P=(K+1)(K+2)/2` points | 1-D log-odds `λ=logit(f_g)`; AMBIG marginalizes the RNA tilt τ |
| memory @ 1.5M nodes | `O(m·K²)` ≈ **>100 GB → OOMs at genome scale** | `O(m·K)` ≈ **~3 GB** |
| model space | linear-fraction priors/messages | **log-density** (multiplicative noise, the right model for counts) |
| status | production reference | implemented, **regresses under capture** (this doc, §5) |

The `logodds` backend is the strategic direction: it is the only one that runs at genome scale, and
log-density is the correct noise model for HTS counts (PCR doubling → log-normal). It was verified
**O(1/K)-equivalent** to the lattice on the same model (Phases 0–2). The lattice remains the default +
A/B reference until logodds reaches parity, then is retired. Plan:
`docs/calibration/log_density_overhaul_plan.md`; design: `log_density_1d_solver_design.md`.

---

## 3. The log-density overhaul (logodds)

Migrates the model into log space, refit consistently (the §3.3 "fit-in-log, solve-in-log" fix):

- **Messages** are log-density Gaussians: a source imputes its component **density** `ρ_c_src`; the
  receiver expresses it in its own log-f_c solve frame via its own density base. No `[0,1]` mode clip.
- **Precision** composes in log space: `prec = 1/(Var(log ρ_c)_src + σ²_bio + pois)`. `pois = 1/count`
  (the source's log-density sampling variance) — **no `+1` floor**: a non-detection (count→0) has ~∞
  variance → ~0 precision → sends no message. ("Zero density is not a measurement.")
- **σ²_bio** (the between-node density spread) refit directly in log space (`var(log ρ)` vs the level).
- **Global** is a Gaussian on `log ρ_g`; M-independent precision `N_log = 1/(var_mean_log + σ²_g_log)`,
  `var_mean_log = 1/(1+G)` (the delta-method 1/count — transparent, no `trigamma`).
- **Grid**: fixed log-odds window `λ∈[−L,L]` (L=10; `log_expit` keeps it exact, no underflow cap).
- **Reference measure** unchanged (Jeffreys + the uniform-f Jacobian); only the *empirical* priors moved
  to log-density. The spliced/mature floor is node-local; messages carry **nascent** only.

---

## 4. Benchmark

Suite: `quick_3to1_5mb` (16 conditions = gDNA{300,none} × ss{0.99,0.50} × nascent{none,rnd} ×
capture{off,on}). Two views:

### 4a. Calibration accuracy — per-region gDNA mass vs oracle (the prior the EM consumes)

Signed error `(Σ solved gDNA mass − Σ oracle gDNA mass) / oracle`. Measured directly (dissect cache):

| condition | lattice | logodds (current) | note |
|---|---|---|---|
| gDNA300 cap-ON ss99 (flagship) | −0.4% | **−7.5%** | AMBIG exons under-called |
| gDNA300 cap-OFF ss99 | −0.9% | −6.2% | |
| gDNA300 cap-ON ss50 (unstranded) | −7.4% | **−28.2%** | worst — no strand, messages dominate |
| gDNA300 cap-OFF ss50 | −12.6% | **−0.9%** | **logodds BEATS lattice** |
| gDNA300 cap-ON ss50 +nascent | −0.2% | −23.4% | |
| zero-gDNA cap-OFF ss50 | 26K phantom | 64K phantom | phantom worse |

**Read:** off-capture, logodds is excellent (beats the lattice) — the log-density framework works where
gDNA is genomically smooth. **Under capture it under-calls gDNA**, badly when unstranded. The
non-detection precision fix (`pois=1/count`) already recovered off-cap ss50 (−17.5% → −0.9%) and halved
the unstranded cases — proof the direction is right.

### 4b. End-to-end net fragment flow (calibration → EM → assignment)

Net gDNA→RNA leak after the per-locus EM (`net_to_rna`; + = gDNA mis-assigned to RNA, − = RNA siphoned
to gDNA). Both backends, fresh quants over all 16 conditions:

| gdna | cap | ss | true gDNA | lattice net (leak%) | logodds net (leak%) | |
|---|---|---|---|---|---|---|
| 300 | off | 0.5 | 3.0M | 126,725 (4.2%) | **5,345 (0.2%)** | **LO ≫ better** |
| 300 | off | 0.99 | 3.0M | −15,228 (−0.5%) | 15,239 (0.5%) | comparable |
| 300 | **on** | 0.5 | 3.0M | 158,936 (5.3%) | **718,652 (24.0%)** | **LO ≫ worse** |
| 300 | **on** | 0.99 | 3.0M | 91,853 (3.1%) | 222,388 (7.4%) | LO worse (flagship) |
| none | off | 0.5 | 0 | −171,863 | −262,984 | LO worse (phantom siphon) |
| none | off | 0.99 | 0 | −10,213 | −5,117 | LO better |
| none | on | 0.5 | 0 | −116,864 | −277,627 | LO worse |
| none | on | 0.99 | 0 | −19,993 | −5,514 | LO better |
| **Σ\|net leak\|** | | | | **711,675** | **1,512,866 (+113%)** | |

**Three clear reads:**

1. **The EM does NOT absorb the calibration error — it propagates ~1:1.** The calibration prior *is* the
   EM's gDNA-vs-RNA split, so an under-called gDNA prior becomes a gDNA→RNA leak directly (cap-on ss0.5:
   calibration −28% → EM 24% leak). Fixing this must happen in calibration; the EM cannot save it.
2. **Off-capture, logodds is better or comparable** — and dramatically so at ss0.5 (4.2% → 0.2%). Where
   gDNA is genomically smooth, the log-density framework is the *better* model. This is the proof the
   direction is right.
3. **On-capture, logodds is 2× worse overall**, entirely from the seam density discontinuity (§5):
   under-called enriched exons → gDNA→RNA leak, worst when unstranded (24%). The zero-gDNA unstranded
   phantom is also worse (the κ=½ siphon). Stranded zero-gDNA is *better* (the honest precision helped).

**Bottom line:** logodds has a *better* model off-capture and a single, well-localized failure
on-capture (the seam). It is not at parity yet (Σ leak 2× the lattice), so the lattice stays the
default; closing the seam gap (§6) is the gating work before logodds can ship and the lattice retire.

---

## 5. The open problem: gDNA density discontinuity at capture seams

Diagnosed at the node level, in **density** terms (unstranded cap-ON):

- Exon **bodies** are at **uniform ρ_g ≈ 15–16** (the probe enriches every exon similarly).
- Enrichment **seams** (exon↔intron/intergenic boundaries) are at **ρ_g ≈ 0–2.5** (the crossing reaches
  into depleted territory).
- The gDNA message asserts density **continuity** (`ρ_g_dst ≈ ρ_g_src`), so a flanking seam imputes its
  *depleted* density onto an enriched exon — biased low by **6× to billions×**.
- **σ²_bio (≈8.7) under-estimates the discontinuity.** The true exon↔seam jump is `(Δlog ρ)² ≈ 64`, but
  σ²_bio is a *single* var~mean curve averaging a **bimodal edge population** (smooth within-the-exon
  regime + discontinuous across-seam) → it lands in the middle → the seam message is *confidently
  wrong* (the truth sits ~2.8σ outside its claim). In log space the far-off mode drags hard; the
  lattice's *bounded* fraction-pull masked the same bias.

This is why the regression is capture-specific: off-capture gDNA *is* genomically smooth (messages
correct), but under capture gDNA density is set by **probe targeting** (a bimodal/discontinuous
structure), and a smooth-continuity message model is the wrong prior.

---

## 6. Recommendations

1. **Disagreement-aware message precision (the lever).** A gDNA message across a large density
   disagreement must self-silence — the communication variance should grow with the observed source↔dst
   density gap (a robust / Student-t "W-regime" weight), so seam messages go quiet *without a capture
   detector*. This is the honest σ²_bio the bimodal edge population demands.
2. **Let `ê(z)` carry the enriched body.** Under capture the exon-body gDNA density is a function of the
   enrichment (the probe), which `ê(z)` already models (it learns the seam-crossing→body transfer). The
   enriched body should come from the global, not from boundary messages that impute the wrong location.
   Verify the `ê`-global is strong enough on AMBIG exons to hold them at their enrichment target.
3. **Validate stranded + unstranded together** (and off-cap, which already works) — the fix must be
   robust across regimes, not tuned to one.
4. **Genome-scale + real data** remain the end goal (the whole reason for logodds): VCaP slice → full
   genome, then the (K,L) sweep and the real `ê(z)` extraction. Production-readiness items surfaced on
   real data: stale index format, GTF duplicate-transcript collapse, mappability (`--no-mappability`),
   and the calibrate-stage memory (which logodds fixes).

---

## 7. Pointers

- Solver: `src/rigel/calibration/simplex_logodds.py`; driver + messages/global: `bp_solver.py`
  (`node_sweep`, `_scan`, `_global_logprior`, `_edge_varmean`); config: `config.py`.
- Plan / design: `log_density_overhaul_plan.md`, `log_density_1d_solver_design.md`.
- Diagnostics: `scripts/debug/{density_message_dissect, flagship_message_dissect, logodds_node_diagnosis,
  phase0_logdens_census, phase1_logodds_equiv}.py`.
