# The NPMLE crush — a node-level dissection of the unstranded+capture f_g collapse

**Status:** diagnosis / workup (2026-07-20). Resumes the KDE-restore WIP
([`gdna_kde_restore_plan.md`](gdna_kde_restore_plan.md), [`kde_vs_npmle_enriched_mode`], archived
[`npmle_struggles.md`](archive/npmle_struggles.md)). The shipped production release uses the **additive KDE** for
the gDNA composition prior (Role B); this branch's default is the **EM-mixture NPMLE**, which still exhibits the
crush documented here. This note traces **one node** all the way through — initialization → pass-0 (strand +
messages) → the exact NPMLE `logprior` arithmetic — so the collapse is followed number-by-number.

---

## 0. The behavior

On **captured + unstranded** scenarios (`ambig_dense_10mb`, `ss_0.50 … capture_on/verystrong`, mid/high gDNA),
the **Phase-2 hyperprior refit makes the gDNA fit worse than not refitting at all** (EMD-to-oracle e.g.
`gdna100_ss0.50_none_capture_on` 0.67 → 2.09). The refit's fitted gDNA-density NPMLE grows a spurious peak at
`log10 ρ_g ≈ −4` — which is simply **every unstranded node's f_g crushed to ~0**, piling their deconvolved
density `ρ_g = f_g·M/E` up near zero. Aggregate on that scenario: captured single-strand exons go pass-0 f_g
**0.435 → refit 0.006** (oracle 0.498); introns **0.511 → 0.008** (oracle 0.937). The stranded (`ss_0.99`)
scenarios are unaffected — the refit tracks the oracle there.

---

## 1. The worked example — node 1055 (`gdna100_ss_0.50_nrna_none_capture_on`)

A **captured, single-strand (POS) exon**:

```
    mass  M = 24235      eff-length E = 2011      total density  ρ_tot = M/E = 12.05   (log10 = +1.08)
    ORACLE: G = 21860 gDNA fragments,  R = 2375 RNA  →  true f_g = 0.902   (true ρ_g = 10.87, log10 +1.04)
```

This node is genuinely **90% gDNA** — 21,860 real gDNA fragments, enriched by capture to a high density
(`log10 ρ_g = +1.04`). Here is what happens to it:

### 1.1 Initialization → pass-0

| stage | f_g | note |
|---|---|---|
| init (signature-binary) | single-strand seed | it *admits* a POS RNA strand, so it is solvable (not a G1 lock) |
| **strand-only** (Beta-Binomial tilt alone) | **0.510** | ss_0.50 is **unstranded** ⇒ the strand likelihood is ~flat ⇒ no composition information |
| local (strand + Jeffreys ½ reference) | 0.510 | the reference is symmetric; still uninformative |
| gDNA message received (mode, prec) | mode **−7.483** (log-f_g), prec **0.022** | a neighbour imputes "f_g ≈ e^−7.5 ≈ 0.0006" — but at **precision 0.022** (σ²_transfer correctly damps the cross-enrichment message) |
| **PASS-0 final f_g** | **0.375** | var_gdna = 3.56 (honestly **imprecise**) |

Pass-0 already leans low (0.375 vs true 0.902) — the unstranded node has no intrinsic signal and takes a weak
low-f_g message — but it stays *honestly uncertain* (var 3.56) and nowhere near collapsed. **The strand channel
carries zero composition information here (unstranded), so whatever prior the refit supplies will decide f_g
almost unopposed.** That is the vulnerability the refit then exploits.

### 1.2 The refit consults the fitted gDNA hyperprior

The Phase-2 refit fits `DensityNPMLE` on pass-0's deconvolved gDNA and re-solves with it as the composition
(gDNA) arm of ψ. The fitted density `P(log10 ρ_g)`:

```
    top mixture modes:  [−2.96, w=0.20]  [−2.93, w=0.08]  [−1.88, w=0.13]  [−1.85, w=0.17]  [−1.83, w=0.13]
    dominant mode at log10 ρ_g ≈ −2.96  (the DEPLETED level)      TRUE ρ_g for this node = +1.04
```

**The fitted density has essentially no mass at the node's true enriched density (+1.04).** Its mass sits at the
depleted level (−2.96) and a mild "enriched" bump at −1.85 — both far below where this captured node's gDNA
actually lives. §2 shows *why* the fit silences the true enriched mode.

### 1.3 THE CRUSH — the exact `logprior` arithmetic

`DensityNPMLE.logprior(f_g; M, E)` evaluates, for each candidate `f_g`,

```
    log ρ_g(f_g)  =  log(f_g) + log(M/E)  =  log(f_g) + 2.489          # 2.489 = ln ρ_tot
    logP_prior    =  interp( log ρ_g(f_g) ,  self.log_rho ,  self.logP )
```

i.e. as `f_g` sweeps, `ρ_g = f_g·ρ_tot` sweeps the density, and the prior scores each. Evaluated on node 1055:

```
      f_g     log10 ρ_g     logP (prior)
    0.001       −1.92        −3.42        ← prior PREFERS this (ρ_g near the −1.85/−2.96 modes)
    0.005       −1.22        −8.97
    0.100       +0.08        −5.84
    0.300       +0.56        −7.55
    0.838       +1.00        −7.15        ← TRUE f_g — sits in the density's clamped LOW TAIL
    0.950       +1.06        −7.34
```

The prior term is **maximised as f_g → 0**: to place `ρ_g` at the fitted density's depleted mode (−2.96), the
node — whose `ρ_tot` is high (+1.08) — must set

```
    f_g*  =  ρ_g^mode / ρ_tot  =  10^(−2.96 − 1.08)  =  10^−4.04  ≈  0.00009  ≈  0.
```

The **true** f_g (0.838–0.902) forces `ρ_g` up to +1.0, which lands in the density's disfavoured high tail
(`logP ≈ −7.1`). With the strand flat (unstranded) and messages weak, nothing resists the prior:

```
    SUMMARY   true f_g 0.902  |  pass-0 0.375  →  REFIT 0.0001        (21,860 gDNA fragments → called ~0)
```

The refit **snapped f_g to `depleted_mode / ρ_tot ≈ 0`**. That is the −4 peak.

---

## 2. Into the NPMLE code — why the fit has no enriched mode, and how it collapses the node

Two code paths matter: **the fit** (`DensityNPMLE.fit`, `additive=False`) that builds the depleted density, and
**the projection** (`logprior`) that applies it. File: `src/rigel/calibration/npmle.py`.

### 2.1 The fit silences the true enriched mode — three compounding mechanisms

The refit's training set is a **depleted-majority** population (this scenario: ~722 depleted intron/intergenic
nodes vs ~470 exons), and every mechanism in the EM fit tilts toward that majority:

1. **Per-node τ-width smears the enriched nodes** (`_cell_loglik`, npmle.py:79). Each cell's likelihood over the
   grid uses that cell's belief width `τ = √Var(log f_g)` as the kernel width. The capture-enriched exons are the
   **low-count, unstranded, imprecise** ones (node 1055: var_gdna = 3.56 ⇒ τ ≈ 1.9 decades) → their likelihood
   bump is **broad and short**, so even before any competition their mode is barely a ripple.

2. **Competitive EM starves the minority mode** (`_em_weights`, npmle.py:113). The component weights are fit by a
   *competitive* mixture EM: cells compete for responsibility. The many sharp, tall depleted cells win the
   responsibility; the few broad enriched cells are down-weighted toward zero. An enriched minority can be
   **competed away entirely** — the exact pathology the additive KDE was designed to avoid
   (`kde_vs_npmle_enriched_mode`).

3. **The aggregate background cell is an `n_regions` tower** (npmle.py:264–276). When a background is measured, a
   single pooled cell is appended at the depleted floor `ρ_floor` with weight **`w_c = n_regions`** (the whole
   intergenic population) and eff-length `ΣE` (genome-scale ⇒ a razor-sharp low mode). This tower dominates the
   EM and pins the density's mass at the depleted level.

Net: `dens_grid = kk @ w_full` (npmle.py:283) is a density whose mass is at −2.96 / −1.85 and whose top,
`hi = max(log ρ̂_g) + 2h` (npmle.py:217), is set by pass-0's *own* deconvolved gDNA — for these under-called
nodes that top is low, so the true enriched level (+1.0) is **above the support and clamps to the low tail**
(`logP[-1]`, npmle.py:296). The enriched mode is not merely down-weighted — it is **off the grid**.

### 2.2 The projection turns "depleted density" into "f_g → 0"

`logprior` (npmle.py:308) is a δ-pin: `log ρ_g = log f_g + log(M/E)`, interpolated on `(log_rho, logP)`. This is
the load-bearing coupling — **it conditions the composition (f_g) on the node's observed total density.** A prior
that says "gDNA density is ≈ ρ_depleted" therefore says, at a node of total density `ρ_tot`, "f_g ≈
ρ_depleted / ρ_tot". The **larger the enrichment `ρ_tot`, the smaller the implied f_g** — so the crush is
*hardest exactly where capture enriches most*. Add a flat strand (unstranded), and ψ ≈ logprior alone ⇒ f_g
snaps to that near-zero value.

---

## 3. Why only unstranded + capture

- **Stranded (ss_0.99) — safe.** The strand Beta-Binomial pins f_g by the intrinsic signal; the weak prior cannot
  override it. The refit tracks the oracle.
- **Unstranded (ss_0.50) — vulnerable.** The strand is flat ⇒ the prior is the only informative term ⇒ f_g snaps
  to the prior's implied value.
- **Capture — the trigger.** `ρ_tot` is high, so `f_g* = ρ_depleted/ρ_tot` is driven hardest toward 0. And
  capture is *exactly* the condition that creates a real enriched-gDNA mode the depleted-fit erases. Uncaptured
  unstranded barely regresses (smaller `ρ_tot`, closer to the fit's mode).
- **Self-reinforcing.** The crushed f_g feeds the next refit's fit, which is then even more depleted-dominated.

The message channel shows the *same* enrichment-blindness in miniature — node 1055's pass-0 gDNA message had mode
−7.483 (a depleted neighbour imputing "≈0 gDNA") — **but σ²_transfer correctly damped it to precision 0.022**, so
it barely moved the node. The **prior has no such damping**: it is applied full-strength at every node regardless
of the enrichment gap between the (depleted-majority) fit and the (enriched) node. The bug is an *un-damped
enrichment-blind prior*.

---

## 4. Root cause & the fix direction

**Root cause.** The gDNA composition prior is a distribution over the **observed** gDNA density `ρ_g = e(x)·d_g`,
pooled across all enrichment levels and — via the EM's competition, τ-width, and `n_regions` tower — dominated by
the depleted majority. Applied through the δ-pin `f_g = ρ_g/ρ_tot`, it imposes a near-**uniform (depleted) gDNA
density** on every node. But hybrid capture makes gDNA density **non-uniform** (enriched at exons), so the prior
is wrong precisely at the enriched nodes, and — where the strand cannot resist (unstranded) — it crushes them.

**The fix (already designed, partially landed):** the **additive occupancy-weighted KDE** for Role B
(`gdna_kde_restore_plan.md`), the representation the shipped production release uses. It removes all three
silencing mechanisms: (1) additive kernels — **no EM competition**, so a minority enriched mode is never competed
away; (2) a **common fixed bandwidth** — occupancy equals height, so the low-count enriched nodes are no longer
τ-discounted; (3) a **weak 1-pseudo-observation floor** at `ρ_floor` — not an `n_regions` tower. With the enriched
mode preserved, node 1055's prior would carry mass at `log10 ρ_g ≈ +1`, and its f_g would no longer be crushed.

**Status / open.** The additive KDE is landed behind `config.gdna_prior_additive` (default OFF, byte-identical).
Its earlier blocker — it faithfully represented a pass-0 *over-call* (`gdna_none` FP) that the EM's heavy
background used to mask — was traced upstream to pass-0; the pass-0 τ-precision + honesty work since then
(`CALIBRATION_STATUS.md`; the pass-0 error is now 80.5% honest-high-variance, not confidently-wrong) removes that
blocker's premise. **Resume point:** re-run the KDE A/B gates (recovery / `gdna_none` FP / stranded
no-regression) on the current honest pass-0, then the default-on flip. The deeper modeling question this
dissection sharpens — whether the prior should condition on the node's **enrichment** (put the density model on
the intrinsic rate `d_g = ρ_g/e(x)`, or stratify by enrichment) rather than pool the observed `ρ_g` — is the
principled generalization of "keep the enriched mode," and worth weighing against the KDE restore.

---

## Appendix — the one-line mechanism

`f_g_refit ≈ ρ_g^{prior-mode} / ρ_tot`. The EM NPMLE's mode is the **depleted** level (competition + τ-width +
`n_regions` tower erase the enriched mode); dividing it by a **capture-enriched** `ρ_tot` gives `f_g ≈ 0`; an
**unstranded** node has no strand to resist. True 0.902 → refit 0.0001. The additive KDE keeps the enriched mode,
so the division no longer collapses.
