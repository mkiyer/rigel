# Message-precision collapse: capture-ON regression & fix plan

**Scope:** `src/rigel/calibration/bp_solver.py`
**Regressed in:** `1863ef57` ("ship disagreement-variance message precision as the single production path")
**Baseline:** v0.6.4 (`2a5d0935`)
**Status:** ⏸️ **PAUSED 2026-07-09** — real solution exists, but masked by a bigger error (the nascent
siphon). Resume later. See "Status & findings" below.
**Published doc:** https://claude.ai/code/artifact/6751163a-dfd2-4758-b7c2-7c52ddd5850b
**WIP code preserved on branch:** `calib-msgprec-eb-maxform-wip` (the EB + env-gated max-form implementation).

## Status & findings (2026-07-09) — PAUSED, this is the resume bookmark

We implemented and benchmarked two candidate fixes on the flagship suite (`quick_3to1_5mb`, 16 conditions).
Both are **empirically confirmed non-ships on this data**, for informative reasons:

1. **Empirical-Bayes shrinkage (Change 1, §6) — `w = 0` on the flagship.** The data-fit weight is zero on
   BOTH the total-density and the per-component gDNA-resolved bases (`Var(log resid²) ≤ π²/2`). With one
   residual per edge, and capture seams being *structural outliers* rather than exchangeable samples, the
   EB decomposition honestly concludes "can't distinguish edges → use the population average," so the
   precision collapses to today's global scalar. Symmetric per-edge shrinkage **cannot fire** with one
   residual per edge. (It also proved that a single agreeing residual isn't reliable evidence of extra
   precision — so the asymmetry below is justified.)

2. **Max form (per-edge self-silencing, belief-free base) — a real lever, but a stranded win / unstranded
   loss.** `σ²_edge = max(resid² − var_loc, σ²_imp + 1/n)`, env-gated `RIGEL_MSG_MAXFORM=1`. Aggregate abs
   error vs the global-scalar baseline, by bucket:

   | bucket | Σ\|gDNA err\| | Σ\|nascent err\| |
   |---|--:|--:|
   | ss99 (stranded) ON  | −5% ✓  | −27% ✓ |
   | ss99 (stranded) off | −36% ✓ | −18% ✓ |
   | ss50 (unstranded) ON  | +22% ✗ | +3% ✗ |
   | ss50 (unstranded) off | +34% ✗ | −31% ✓ |
   | **ALL 16** | **+16% ✗** | **−16% ✓** |

   Self-silencing helps where strand carries the solve (stranded) and hurts where the relay is load-bearing
   (unstranded) — the user predicted this exactly. It halves the worst nascent leaks (`gNONE none ON`
   195k→81k) but by *bluntly* muting the mature-into-intron relay, which also costs unstranded gDNA.

**Why paused:** the dominant, gDNA-independent **nascent siphon** (mature→nascent, e.g. `gNONE ss50 none
ON` = 171,861 mature fragments called nascent vs a true zero) is so large it masks any true message-precision
win. Message precision cannot fix a gDNA-independent error. Fix the siphon first (reference-density /
eff-length + mature-gate), then resume here with a clean signal.

**How to resume:** (a) `git checkout calib-msgprec-eb-maxform-wip` restores the EB + max-form code; (b) the
two viable directions are the max form gated/weighted by *strand information content* (self-silence only
where strand carries the solve), or a per-component basis with more than one residual per edge (aggregate
edges so the EB weight can rise). Validate on the soft 3-pool surplus, NOT `net_flow`.

## At a glance

- **What broke:** the BP-sweep message precision changed from a **per-edge, disagreement-aware** term that
  self-silences at density seams to a **global count-only scalar** with no seam term.
- **Why capture:** probe capture creates the library's largest legitimate adjacent-density jumps (enriched
  on-target exon beside depleted off-target flank). The old form muted messages there; the new one lets the
  enriched node impose its density on its depleted neighbour → gDNA/RNA separation bleeds across the seam.
- **Fix:** restore per-edge, **symmetric** precision via **empirical-Bayes shrinkage** — an edge that
  agrees better than the population earns *above-average* precision, a seam edge self-silences — with NO
  tunable weight. The shrinkage weight is the measured signal-to-noise ratio, not a constant. Self-silencing
  returns; the runaway cannot.
- **Sequencing:** independent of the eff-length inversion fix; run both as separate A/B arms.

## 1. What the message precision controls

The forward–backward sweep passes **density messages** between adjacent region↔boundary nodes; each
message carries a **precision** (inverse variance) that sets how hard the source pulls the destination's
gDNA-vs-RNA composition toward its own. The decisive case is a **density seam** — two genomically adjacent
nodes with sharply different true gDNA densities. Hybrid capture manufactures these (on-target exon ~100×
enriched beside a depleted off-target flank). A message from the enriched node must **not** impose its
density on the depleted neighbour; the precision has to **self-silence** across a seam.

## 2. v0.6.4 — per-edge, disagreement-aware, self-silencing

`node_sweep._scan`, gDNA channel (bp_solver.py:1064, v0.6.4); identical form for RNA±:

```python
pois     = 1.0 / max(n_src, _EPS)                                     # source log-density sampling var
base_var = vbg[lsrc] + pois                                          # sampling + source-belief log-var
s2_edge  = max((mo - lfg_loc[i])**2 - (base_var + vg_loc[i]), 0.0)   # per-edge surprise
pr       = 1.0 / max(base_var + s2_edge, _EPS)
```

`resid = mo − lfg_loc[i]` = message minus the destination's **message-free** local belief (non-circular
anchor). Agreeing edge → `s2_edge ≈ 0` → full precision. Seam → `resid` large → `pr ≈ 1/resid²` → self-mutes.

**Flaw it carried:** `base_var = vbg[lsrc] + pois` used the source's own *belief* variance `vbg`, which
collapses toward 0 as beliefs sharpen → precision runaway (~1e9). That runaway — not the self-silencing —
is what the next change targeted (`a1d0f17a`).

## 3. Current HEAD — the global count-only scalar

`node_sweep._scan`, gDNA channel (bp_solver.py:724); identical for RNA±:

```python
# σ²_msg = σ²_imp + 1/n_src  ⇒  pr = n_src/(n_src·σ²_imp + 1)
pr = n_src / (n_src * sig_imp + 1.0)
```

`sig_imp` (σ²_imp) is one number, fit once by `adjacent_disagreement_variance` (bp_solver.py:474) over
adjacent boundary↔region edges. **No `resid` term** — an edge's precision no longer depends on whether the
message agrees with its destination. Cures the runaway (no belief-derived variance) but deletes the seam
term. The v2 design was a *shrinkage* (blend residual with this floor); what shipped kept only the floor.

**Secondary change, same commit:** intron/intergenic **floor nodes** gained a new *uncapped* density
likelihood `1/(s2_bg + 1/N)` + per-node Jeffreys (`_global_logprior`, bp_solver.py:339-362) that reads
density-excess-over-background as nascent RNA, and were **excluded** from the pass-2 KDE (bp_solver.py:638).
v0.6.4 gave them only a capped (≤ 1 pseudo-obs) floor hyperprior and kept them in the KDE.

## 4. Why it regresses capture-ON specifically

Off-capture, adjacent densities are similar so the missing seam term costs little (clean scenarios never
regressed). Under capture:

- **Enrichment bleeds across the boundary.** An enriched on-target exon (high count → high precision under
  the new rule) sends a confident message into its depleted neighbour, dragging its gDNA density up. The
  old per-edge form muted exactly this message.
- **One global σ²_imp is the wrong knob.** Fit once over the whole chain, it is inflated by capture's
  heterogeneity and applied uniformly — over-trusting relay in quiet regions, under-muting it at seams.
- **The floor-node change reinforces it.** Under capture, on-target introns sit far above the
  off-target-estimated background, so their enriched gDNA is mis-attributed to nascent — the same leak
  direction the weakened relay encourages.

**Why it slipped through:** validated on the hard-label `net_flow` metric, which is insensitive to a
soft-prior shift (a real change moves the soft 3-pool counts by tens of thousands of fragments while
`net_flow` stays byte-identical). Validate the fix on the soft 3-pool surplus, never `net_flow`.

## 5. Evidence — flagship 3-pool net flow

`quick_3to1_5mb`, ss 0.99, 3:1 gDNA. Net surplus = assigned − true fragments per pool.

| condition | gDNA net | nascent net | mature net | gDNA→nascent |
|---|--:|--:|--:|--:|
| capture-OFF · nascent | +23,705 | +46,866 | −70,571 | −19,099 |
| capture-ON · nascent | **−66,939** | **+57,300** | +9,639 | **+36,112** |
| capture-OFF · no-nascent | −41 | +72,524 | −72,483 | 0 |
| capture-ON · no-nascent | **−110,491** | **+176,660** | −66,169 | **+78,648** |

Capture flips gDNA→nascent from −19k (off) to +36k (on); the no-nascent arm calls 176,660 fragments
nascent against a true zero. This leak has **two** contributors — the eff-length inversion (fixed
separately) and this message-precision regression. The §7 A/B attributes the split.

## 6. The fix — empirical-Bayes disagreement shrinkage (parameter-free)

Restore per-edge, **symmetric** precision: an edge that agrees better than the population earns precision
*above* the average; a seam edge self-silences *below* it. The obstacle was the shrinkage weight, which
looks like a tunable. It is not — estimating a variance from residuals has a *known* amount of noise, so
the weight is the empirical-Bayes signal-to-noise ratio, measured once.

Work in **log-variance** space (variances are positive; their sampling noise is multiplicative). For each
adjacent edge, `d_i = log(resid_i²)` is a one-sample estimate of that edge's log-variance, where
`resid_i = message − dst message-free local belief` (the non-circular anchor). A single squared residual is
a draw from `σ²_i·χ²₁`, so

```
d_i = log σ²_i + log χ²₁ ,   Var(log χ²₁) = ψ′(½) = π²/2 ≈ 4.93     # CONSTANT: the irreducible noise
                                                                    # of a variance from ONE residual
```

The population spread decomposes, `Var(d_i) = tau2_between + π²/2`, so the weight falls out with nothing to
choose (fit ONCE over all adjacent edges, pre-sweep):

```python
PI2_2        = math.pi**2 / 2.0                       # ψ'(1/2): within-edge 1-sample log-variance noise
tau2_between = max(var(d) - PI2_2, 0.0)               # between-edge dispersion of true log-variances
w            = tau2_between / (tau2_between + PI2_2)   # signal fraction in [0,1) — MEASURED, not tuned
mu           = mean(d)                                 # population mean log(resid^2) (bias-corrected below)
# per edge, during the scan:
d_i          = log(max(resid_i**2, samp_i))            # floor at Poisson sampling var (not a knob)
log_s2_edge  = mu + w * (d_i - mu)                     # shrink each edge toward the population mean
pr_i         = exp(-(log_s2_edge - E_LOGCHI2))         # E[log chi2_1] = -1.27 bias correction
```

- **Ceiling/floor framing.** The population-average variance `σ²_imp` (≈ `exp(mu − E_LOGCHI2)`) is the
  shrinkage *target*, no longer a hard ceiling. Genuine agreement (`d_i < mu`) drives `σ²_edge < σ²_imp` →
  precision *above* average; a seam (`d_i > mu`) drives it below. The random-non-adjacent-pair disagreement
  variance is the natural upper bound on `σ²_edge` (a message no better than random is dead) — an optional
  clip, deferred to Change 2.
- **Generalizes both extremes; the data picks.** `w = 0` (all spread is sampling noise) ⇒ every edge
  collapses to `σ²_imp` — exactly today's global scalar. `w → 1` (edges genuinely differ) ⇒ each edge is
  trusted near its own observation. Under capture the field is heterogeneous, so `w` sits high — more local
  trust exactly where it is needed.
- **Cannot run away.** `π²/2 > 0` ⇒ `w < 1` strictly, so one lucky-agreement residual is always shrunk
  back. High precision comes only from high count + genuine agreement (non-circular: the residual is vs the
  *message-free* belief), never the circular belief-sharpening that caused the v0.6.4 runaway.

**Cost, honestly:** log-space + a sampling floor on `resid²` (floored at the Poisson variance `1/n` — you
cannot observe agreement finer than the count allows; a floor, not a knob) + the `E[log χ²₁] ≈ −1.27` bias
correction. One modeling assumption — **one effective residual observation per edge** — is where `π²/2`
comes from; still count-derived, never tuned.

## 7. Implementation plan (Change 1 = the EB form)

1. **Reinstate destination local-belief anchors** — the residual needs each destination's *message-free*
   local belief (`lfg_loc`/`vg_loc` + RNA± equivalents). Confirm they survive the `node_geometry`
   refactor; recompute per node if the collapse removed them. Anchor on the message-free belief (avoid the
   echo chamber).
2. **Fit the EB weight once, pre-sweep** — extend `adjacent_disagreement_variance` to also return the
   between-edge log-variance dispersion `tau2_between = var(log resid²) − π²/2` and hence `w` and `mu`,
   from the same adjacent-edge residual population it already builds. No new inputs.
3. **Swap the precision in `_scan`** — replace `pr = n_src/(n_src·sig_imp+1)` with the log-space EB
   shrinkage above, in all three channels (gDNA, RNA⁺, RNA⁻), using the per-edge `resid_i` and the
   pre-fit `(w, mu)`. Floor `resid²` at the Poisson sampling variance.
4. **Decide the floor-node change on its own** — evaluate reverting the uncapped intron/intergenic density
   likelihood to the v0.6.4 capped floor (+ KDE re-inclusion) as a *separate* arm; may be partly redundant
   once the relay self-silences.
5. **A/B on the soft metric** — in-process A/B over the four flagship conditions, measuring the 3-pool net
   surplus (not `net_flow`). Target: recover the v0.6.4 capture-ON separation. Report each arm's share
   (message-precision vs floor-node vs eff-length).
6. **Unit tests** — (a) on a *homogeneous* synthetic field, `w → 0` and the EB precision reproduces the
   global scalar `1/σ²_imp` (the limiting case); (b) an above-average-agreement edge gets precision *above*
   `1/σ²_imp`; (c) a high-`resid` seam edge self-silences. Regenerate the calibration goldens.

Queued behind the effective-length inversion fix (landed 2026-07-09), which removes the other capture-ON
leak contributor. The two are independent; treat the v0.6.4 capture-ON 3-pool separation as the target.
