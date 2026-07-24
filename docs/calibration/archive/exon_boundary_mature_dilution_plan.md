# Exon ↔ IE-boundary message: the MATURE-DILUTION rule — derivation, evidence, implementation plan

**Status:** derived + empirically validated on the grounded full-transcript toy; **not yet implemented**.
**Scope:** pass-0 message propagation ONLY (`calib_refit_iters=0`). No Phase-2 hyperprior work.
**Date:** 2026-07-22.

---

## 1. The problem this solves

Pass-0 message propagation fails across the hybrid-capture **enrichment cliff**. The audit
(`scripts/debug/msg_audit.py`, grounded toy, population priors injected) traced the failure precisely:

| edge | mode | imputed f_g | oracle | prec_g | σ²_T |
|---|---|---|---|---|---|
| `intgc → B` | (structural) | 0.93 ✓ | — | 7.85 | 0.12 |
| **`B → exon`** | density | **0.0016** ✗ | **0.747** | 0.03 | **39** |
| `exon → B` | geo-mean | — | — | 0.01 | 39 |
| `B → intron` | shift | 0.37 | — | 0.01 | 0.12 |

The intergenic ρ_bg anchor reaches the first boundary correctly, then **dies at the exon**: the density mode
computes `f_g ≈ ρ_g^src·E_g^dst / md_dst`, dividing a *depleted* source gDNA density by the *enriched* exon's
total density, under-imputing gDNA by the enrichment factor e(x) (~600× in the toy). ρ_bg therefore never
enters the transcript body, every interior node falls to the no-evidence ~0.5, and the interior introns
under-call by −0.41 (0.51 solved vs 0.92 oracle) — the largest single error in the transcript.

The **intron → IE-boundary** half is already fixed (composition invariance / shift — landed, 260 calibration
tests pass): an intron and its boundary crossing have identical active components (gDNA + nascent), so the
fraction scales across the cliff. What remains is the **exon ↔ boundary** direction, where the exon carries a
third component (mature) that the boundary crossing does not.

---

## 2. Derivation

Capture is nucleic-acid-agnostic: every species at position *x* scales by the same enrichment `e(x)`.
Let `ρ_bg` = gDNA (uniform, genomic), `ν` = nascent (continuous along the gene body), `μ` = mature (exonic).

| node | unenriched composition | observable |
|---|---|---|
| intron | ρ_bg + ν | D_I |
| **boundary, unspliced (crossing)** | ρ_bg + ν | **D_B** |
| **boundary, spliced (junction)** | μ | **S_B** |
| exon, unspliced | ρ_bg + ν + μ | D_E |

A **mature fragment cannot cross an intron-exon junction contiguously** — it is spliced, so it jumps the
intron. Hence the boundary's unspliced crossing is mature-FREE. (Verified: §3.1.)

The boundary's **total** (unspliced + spliced) therefore has exactly the same component set as the exon:

```
D_B + S_B = (ρ_bg + ν + μ)·e_B        D_E = (ρ_bg + ν + μ)·e_E
```

so `e_E/e_B = D_E/(D_B + S_B)` — the enrichment ratio is determined. **But it is not even needed**, because it
cancels in the composition:

```
f_g^B = ρ_bg/(ρ_bg + ν)          f_g^E = ρ_bg/(ρ_bg + ν + μ)

   f_g^B / f_g^E  =  (ρ_bg + ν + μ)/(ρ_bg + ν)  =  (D_B + S_B) / D_B
```

### The rule

```
    f_g^boundary = f_g^exon · (D_B + S_B)/D_B          (exon → boundary: remove mature)
    f_g^exon     = f_g^boundary · D_B/(D_B + S_B)      (boundary → exon: restore mature)
```

Message modes are already `log f`, so this is a single **additive** term. Define, per boundary node *b*:

```
    c_b  =  log( (D_B + S_B) / D_B )        the MATURE-DILUTION term  (≥ 0)
```

and the exon↔boundary message mode becomes the ordinary composition shift ± `c_b`:

```
    mode = log f_g^src + [ log(E_g^dst/E_g^src) − log(E_r^dst/E_r^src) ]  ±  c_b
                         └──────── the existing shift ────────┘      + for exon→B, − for B→exon
```

**Properties.** Enrichment-invariant (e cancels — the cliff height is irrelevant); **zero free constants**;
the mature correction uses only the boundary's OWN measured spliced/unspliced split, so it does **not** require
extrapolating within-exon mature via eff-length ratios — the step that failed before (§9, "300× incomplete").

---

## 3. Empirical validation (grounded toy)

Harness: `scripts/debug/mature_dilution_check.py` on the full 3-exon transcript with intergenic ends
(`toy_inject.build_toy`), population priors injected from the `ambig_dense_10mb` cache. Oracle truth from the
simulator BAM read names (`gdna` / `nrna_*` / mature), boundary composition measured by **contiguous** crossing
(`get_blocks`), exon composition by contained mate+TLEN span.

### 3.1 Premise: mature does not cross contiguously — **CONFIRMED**

Oracle mature-crossing fraction at every splice-junction boundary: **0.000**. The boundary crossing is
gDNA + nascent only, exactly as the derivation requires.

### 3.2 The identity — holds to ~2%

Base condition (capture ON, gDNA=3.0, nascent=200, mature=100), per junction:

| junction | LHS = f_g^B/f_g^E (oracle) | RHS = (D_B+S_B)/D_B (measured) | rel. error |
|---|---|---|---|
| 7000 | 1.0882 | 1.1067 | +1.7% |
| 11000 | 1.0616 | 1.0876 | +2.4% |
| 12000 | 1.0641 | 1.0910 | +2.5% |
| 16000 | 1.0692 | 1.0844 | +1.4% |

**~2% agreement, zero fitted constants.**

### 3.3 The FRAME matters — use summed-eff, not per-side

Two conventions were tested for forming the densities:

| convention | RHS | verdict |
|---|---|---|
| **summed-eff**: `D_B=(m_l+m_r)/(e_l+e_r)`, `S_B=(s_l+s_r)/(es_l+es_r)` | 1.084–1.107 | **correct (~2%)** |
| per-side density sum: `S_B = s_l/es_l + s_r/es_r` | 1.169–1.213 | over-corrects ~2× |

Spliced mass lands on **one face only** (the exon side: `spl_l=196, spl_r=0`), so the per-side sum divides by a
single side's eff-length and double-counts. **The implementation must use the summed-eff frame.**

### 3.4 Probe placement is LOAD-BEARING — a real risk, discovered here

Splice-junction reads live at exon **edges**. A centred partial probe excludes them:

| capture / probe layout | mature reads | **junction-spanning (spliced)** |
|---|---|---|
| capture OFF | 4592 | **331** (7.2%) |
| capture ON, probes = exon **centre 50%** | 13310 | **1** (0.008%) |
| capture ON, probes = **full exon** | 13000 | **685** (5.3%) |

With centred probes the junction channel is annihilated ~300×, `S_B → 0`, and the correction silently becomes a
no-op (`c_b → 0`) — the message then transfers the exon's mature-contaminated composition **at full
confidence**. This is the dominant failure mode of the rule and **must** be handled by the precision model
(§5), not ignored. It also confirms the owner's standing point that probe placement varies per junction and
cannot be assumed.

### 3.5 Grid + probe-layout + FL sweep (22 cells, `n_rna=4000`)

**`mature-cross = 0.000` in ALL 22 cells** — the premise is universal, not condition-specific.

capture × gDNA × nascent × mature (relative error of RHS vs oracle LHS):

| capture | gDNA | nascent=20,mat=20 | n=20,m=200 | n=200,m=20 | n=200,m=200 |
|---|---|---|---|---|---|
| off | 0.05 | *degenerate* | *degenerate* | −30.6% | *degenerate* |
| off | 3.0 | −13.6% | −17.3% | **+0.4%** | −13.6% |
| on | 0.05 | −13.1% | +24.3% | −15.5% | −13.1% |
| **on** | **3.0** | **+6.1%** | **+11.0%** | **+1.8%** | **+6.1%** |

*degenerate* = at gDNA=0.05 with this depth the boundary crossing receives **zero** gDNA fragments, so the
oracle `f_g^B = 0` and the ratio is undefined — a measurement floor, not a model failure. It marks exactly the
regime §5's precision term must suppress.

The **production-like hard case (capture ON, gDNA 3:1)** is the best-behaved row: **+1.8% … +11%**.

**Residual is counting noise, not model bias.** Same base cell at two depths: **~2% at `n_rna=20000`**
(§3.2) vs **~6% at `n_rna=4000`**. Errors change sign across cells (+24% … −31%) rather than drifting
systematically — the signature of variance, not bias. The boundary crossing is intrinsically tiny (eff-len
~60 bp vs the exon's ~1000 bp), so `D_B`/`S_B` are the noisiest quantities in the chain. **This is why §5 is
not optional.**

Probe layout (capture ON) — the owner's point that probes may sit *over* the junction:

| probe layout | rel. error | note |
|---|---|---|
| full exon | **+3.4%** | baseline |
| centre 50% | −2.1% | only 3/4 junctions survive (§3.4 depletion) |
| **over the junction** | **−18.6%** | boundary enriched *above* the exon; rule under-corrects |

The junction-probe case inverts the enrichment ratio (boundary > exon). The rule survives directionally — `e`
does cancel — but under-corrects ~19%, the largest systematic deviation found. **Flagged as a known
limitation**, not a blocker.

FL sweep (the direct test of the eff-length frame, since `E_g` vs `E_r` is what the shift corrects):

| RNA_FL / gDNA_FL | 200/100 | 300/100 | 150/150 |
|---|---|---|---|
| rel. error | +3.4% | +11.2% | −3.2% |

Frame is robust across a 2× FL spread and across `E_g = E_r`. No blow-up.

**Verdict: proceed to implementation.** The premise is exact, the identity is unbiased, and the residual is
depth-driven variance that the precision model converts into (correct) low confidence.

---

## 4. Implementation

All changes in `src/rigel/calibration/bp_solver.py`, `node_sweep` / `_scan`.

### 4.1 Precompute the per-boundary dilution term (vectorized, before the scans)

Alongside the existing `rho_g_cross` precompute:

```python
# MATURE-DILUTION term per boundary (docs/calibration/exon_boundary_mature_dilution_plan.md §2).
# D_B = unspliced crossing density, S_B = spliced (mature) density — SUMMED-EFF frame (§3.3).
_msum = MS[0] + MS[1]
_egsum = np.maximum(EG[0] + EG[1], _EPS)
_ssum = SP[0] + SN[0] + SP[1] + SN[1]
_essum = np.maximum(ESP[0] + ESP[1], _EPS)
D_B = _msum / _egsum
S_B = _ssum / _essum
mature_dilution = np.where(~_is_reg_arr & (D_B > _EPS), np.log1p(S_B / np.maximum(D_B, _EPS)), 0.0)
```

### 4.2 Apply it on exon↔boundary edges, with sign by direction

In `_scan`, on the gDNA message mode:

- `exon → boundary` (`_is_bnd[i] and _ex_s`): `mo += mature_dilution[i]`
- `boundary → exon` (`_ex_d and _is_bnd[lsrc]`): `mo -= mature_dilution[lsrc]`

and **enable `use_shift` for both**, since with the mature removed the components match:

```python
use_shift = (not _ex_d and not _ex_s and (not _is_bnd[i] or is_intron_node[lsrc]))   # existing
            or (_is_bnd[i] and _ex_s)          # exon → boundary  (mature removed by +c_b)
            or (_ex_d and _is_bnd[lsrc])       # boundary → exon  (mature restored by −c_b)
```

### 4.3 Retire the geo-mean crossing

`rho_g_cross` (the unweighted `√(ρ_g^exon·ρ_g^intron)` pre-scan hack) is **superseded**: it stood in for the
missing exon→boundary imputation and bypassed the precision machinery entirely. With both directions carrying
real, precision-bearing messages, the boundary's two flanks are reconciled by the existing α⊗β integrator —
which is already a **precision-weighted geometric mean** (`_comb`: precisions add, modes are log-fractions).
Remove `rho_g_cross` and its `_ex_s` RNA-suppression branch **only after** §6 gates pass; keep them until then.

---

## 5. Precision model (the part that makes it safe)

`c_b` is estimated from counts and must carry its own uncertainty, or §3.4 becomes a silent corruption.
With `n_s` = raw spliced fragment count and `n_d` = raw unspliced crossing count at the boundary, and
`r = S_B/D_B`:

```
    Var(c_b) ≈ [ r/(1+r) ]² · ( 1/(n_s+1) + 1/(n_d+1) )
```

(Poisson counting noise on both channels, propagated through `c = log(1+r)`; the Gamma-posterior `n+1`
keeps it finite at `n=0`.) This enters the existing precision denominator:

```python
pr = _pred_precision(sm, v_logfg + var_c, s2t)     # var_c = Var(c_b) on exon↔boundary edges, else 0
```

so a count-starved junction (§3.4) yields a **low-precision** message that the α⊗β integrator automatically
down-weights against the clean intron→boundary message. **No new constants, no gate.**

Measured on the base cell: `r ≈ 0.107`, `n_s ≈ 196`, `n_d ≈ 1431` → `Var(c_b) ≈ 5.4e-5` — negligible where the
junction is well sampled, exactly as intended. In the *degenerate* grid cells (§3.5) it is `n_d`, not `n_s`,
that starves, and the existing `1/M_src` term already handles that.

**Known residual risks — documented, not silently handled:**

1. **Probe-depleted junctions (§3.4).** When the junction channel is depleted by *probe placement* rather than
   merely sparse, `n_s ≈ 0` makes `r ≈ 0`, so both `c_b ≈ 0` **and** `Var(c_b) ≈ 0` — an uncorrected message
   at full confidence. The counting model cannot distinguish "no mature" from "mature unobservable here". A
   diagnostic (compare the boundary's spliced count against the adjacent exon's mature evidence; flag gross
   inconsistency) is proposed as **follow-up**, not part of this change.
2. **Probes over the junction (§3.5).** ~19% systematic under-correction when the boundary is enriched above
   the exon. Direction is still right; magnitude is not. Accept for now; revisit if the audit shows it
   dominating.

---

## 6. Test plan / acceptance gates

1. **Identity regression test** (new, `tests/calibration/`): assert on the toy that `(D_B+S_B)/D_B` reproduces
   the oracle `f_g^B/f_g^E` within tolerance, and that oracle mature-crossing is ~0. Guards the premise and
   the frame (§3.1–3.3) — this is the class of bug that the summed-eff-vs-per-side error would have shipped.
2. **`msg_audit.py`**: `boundary→exon` imputed f_g must move from 0.0016 toward the 0.747 oracle; the interior
   introns must rise off the 0.51 no-evidence floor toward 0.92.
3. **Grid + FL sweep** re-run (§3.5) after implementation — mode error, not just the identity.
4. **`gdna_none` phantom guard** — hard gate, must not regress.
5. **Golden + full suite.** This is a behavioural change: goldens will move. Regenerate **only** after the
   benchmark A/B below is reviewed.
6. **Benchmark A/B** (`scripts/debug/benchmark_ab_report.py`) — the intron→boundary shift already in place
   plus this change together alter §14 geo-mean behaviour; validate on the real suite before committing.

## 7. Sequencing

1. Land §4.1 + §4.2 + §5 behind the existing structure (geo-mean retained) → run gates 1–4 on the toy.
2. If gates pass, retire the geo-mean (§4.3) → re-run gates 1–4.
3. Benchmark A/B (gate 6), review, then goldens (gate 5) and commit.

**Open question for the owner:** the TSS/TES boundaries (intergenic↔exon) are a distinct case — the intergenic
flank carries gDNA only (no nascent), so neither composition invariance nor the mature-dilution rule applies
unmodified. They are currently left on the density path and are **out of scope** for this change; worth a
separate derivation if the audit shows them mattering.
