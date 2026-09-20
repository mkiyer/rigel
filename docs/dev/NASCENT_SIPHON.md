# The nascent RNA siphon — what causes it, what does not, and what to build

**Status.** Diagnosis complete and priced; no repair landed. Written 2026-09-19/20 against `f77d43e1`.
Nothing in `src/` changed while this was measured. The settled parts have moved to their permanent
homes — the derivation to `EQUATIONS.md` §9b, the record and the numbers to
`ISSUES: nascent-siphons-gdna-under-capture` and `ISSUES: per-transcript-prior-lane`, the lesson to
`TRAPS: an-oracle-column-that-omits-a-population`. This file is the argument joining them up, and it is
`docs/dev/`: a working document, authoritative for nothing.

---

## 1. The defect in one table

`quant_accuracy.py --arm base`, fractional assignment, est − true in fragments:

| condition | nascent Δ | gDNA Δ | annotated Δ |
|---|---:|---:|---:|
| `g05 ss.99 ON` | +69,268 | −66,516 | −2,752 |
| `g50 ss.99 ON` | +541,216 | −534,656 | −6,560 |
| `g98 ss.99 ON` | +590,406 | −611,173 | +20,767 |
| `g50 ss.99 OFF` | −66,752 | +70,572 | −3,820 |

The synthetic nascent entities and gDNA exchange fragments almost one for one while the annotated pool
barely moves. The worst IN-SCOPE condition is `g98 ss.99 ON`, where nascent reads **99.6× its truth**.

---

## 2. What is actually happening

### 2.1 The mass is on entities that hold nothing

At `g50 ss.99 ON`, of the +541,116 fragments of nascent error:

| | n | true | estimate | Δ |
|---|---:|---:|---:|---:|
| LIVE shadows (true nascent > 0) | 1,379 | 150,405 | 238,667 | +88,262 |
| **SILENT shadows (true nascent = 0)** | **5,540** | **0** | **452,854** | **+452,854** |

**84 % of the siphon sits on spans that emitted no fragment at all**, and **97.5 % of that sits on
shadows whose own GENE is expressed** (4,805 entities, 91.8 fragments each) rather than on shadows of
silent genes (681 entities, 16.9 each). It is diffuse — the top 100 of 6,919 carry 32 %, the top 1,000
carry 90 % — so it is a systematic bias and not a handful of pathological loci.

### 2.2 The mass comes from gDNA, and the control is exact

At **every locus with zero certified gDNA the silent-shadow mass is EXACTLY 0** — 33 such loci at
`g50 ON`, 44 at `g50 OFF`, 28 and 14 at `g98`. The leak is gDNA-fed, without exception.

### 2.3 The root cause: `θ_n = 0` is an unstable fixed point

The EM gives a whole **MultiLocus** — a connected component of transcripts linked by shared fragments —
ONE gDNA component with ONE opportunity `L_g` covering the entire component. Every synthetic nascent
entity inside it carries only its own gene's span `L_n`. A connected component is a union of gene spans,
so **`L_g > L_n` structurally**, and the gap grows with the number of genes it holds.

In the E-step a component's weight is `θ_c / L_c`. Writing `a = θ_g/L_g` and `b = θ_n/L_n`, near `b = 0`

```
    (b/a)  ←  (b/a) · L_g/L_n      per iteration
```

so `L_g/L_n > 1` makes `θ_n = 0` **unstable**: a shadow holding nothing climbs off zero and settles
where only the strand channel and the gDNA pseudocount stop it. On a pool of `N` fragments that only
gDNA and the shadow can explain, split evenly by genome strand, the fixed point in `r = b/a` is

```
    r · L_n/L_g = [ ss·r/(ss·r + ½) + (1−ss)·r/((1−ss)·r + ½) ]
                / [    ½/(ss·r + ½) +        ½/((1−ss)·r + ½) ]
```

whose only root at `L_n ≥ L_g` is `r = 0` — the correct answer, and stable.

**Put to the shipped solver** (20,000 fragments, all gDNA in truth, `ss 0.99`, no pipeline):

| `L_g/L_n` | 0.5 | 1.0 | 2.0 | 4.33 | 6.22 | 20 |
|---|---:|---:|---:|---:|---:|---:|
| no gDNA prior | 0 | **0** | 35.2 % | 48.1 % | 52.7 % | 82.0 % |
| `gdna_prior = N/2` | 0 | 0 | 17.2 % | 39.0 % | 44.4 % | 58.2 % |

4.33 and 6.22 are the ladder's own paired medians at `g50 ss.99` OFF and ON. The closed form and the
native solver agree to two decimals, so this is the model's ML answer and not an EM artefact. **The
threshold is exactly 1 — derived, not chosen.** Gated in `tests/test_estimator.py`.

Remove the strand term and the contest is fully degenerate: the shorter component takes **everything**
(1.0 at every ratio above 1, exactly 0.5 at 1, 0 below). Strand does not open this channel; it only
bounds it.

### 2.4 The ratio is set by how many genes share a locus

At `g50 ss.99 ON`, the silent-shadow false positives by the gene count of their locus:

| genes in the locus | shadows | FP mass | share | median `L_g/L_n` |
|---|---:|---:|---:|---:|
| 1 | 1,372 | 64,317 | 14.2 % | 2.50 |
| 2–3 | 1,911 | 134,558 | 29.7 % | 8.11 |
| 4–9 | 2,013 | 222,722 | 49.2 % | 15.75 |
| 10+ | 182 | 31,251 | 6.9 % | 17.25 |

**85.8 % of the leak is in multi-gene loci.** A single gene's own shadow is nearly harmless; the damage
is the pooling.

### 2.5 End-to-end confirmation, and why a length knob is not the repair

Scaling **only** the shadows' EM effective length by `λ` at `g50 ss.99 ON` — a falsification probe, not
a mechanism:

| `λ` | nascent (true 150,432) | gDNA (true 5,000,000) | transcript Σ\|Δ\| |
|---|---:|---:|---:|
| 1 (shipped) | 691,709 | 4,465,354 | 252,376 |
| 2 | 18,023 | 4,903,115 | 386,614 |
| 6.22 | 738 | 4,916,794 | 390,454 |

Doubling `L_n` alone removes **97 %** of the siphon and returns 437,761 fragments to gDNA — a cliff, the
signature of a threshold rather than a bias. It also **destroys the true nascent signal** (the
transcript table worsens 252,376 → 386,614): it kills the live entities with the dead ones. At
capture-OFF the same knob moves the same mass (`λ = 2`: nascent 946,781 → 449,852), so the channel is
open in both regimes.

---

## 3. Capture does not reverse the arbitration

The old framing — "capture flips the sign" — is wrong. The same channel is open in both regimes at the
same order of magnitude, and the library-level sign flip is arithmetic:

| `g50 ss.99` | capture OFF | capture ON |
|---|---:|---:|
| false positive on SILENT shadows | +257,002 | +452,854 |
| Δ on LIVE shadows | −323,623 | +88,262 |
| **net nascent Δ** | **−66,621** | **+541,116** |
| TRUE nascent pool | 1,013,400 | 150,405 |
| shadow-exclusive gDNA pool | 2,212,935 | 2,142,952 |
| leak rate on that pool | 11.6 % | 21.1 % |
| mass-weighted `L_g/L_n` | 3.57 | 9.70 |

Capture changes two things and neither is a direction:

1. it raises `L_g/L_n` **2.7×** — the shadow contracts harder than the pooled component opportunity —
   which roughly doubles the leak rate;
2. it collapses the TRUE nascent pool **6.7×** (nascent RNA is intron-heavy and unprobed), so the
   compensating under-call on live shadows that was masking the false positive off capture disappears.

⚠ The shadow-exclusive pool itself barely changes (2.21 M → 2.14 M). Under capture it is made of
**crossing** fragments — 4,336,797 of them against 385,866 off capture, of which 2,072,664 straddle an
exon|intron edge, which no mature isoform can produce.

---

## 4. What is ruled out, with the killing number

**The ruler's witness geometry is not the driver.** The shipped ruler ALREADY reproduces the 13.5×
annotated-vs-synthetic capture gap: contraction factors at `g50 ss.99 ON` are 0.4491 (annotated median)
and 0.0347 (synthetic), ratio **0.0773**, against the simulator's own anchored 41.0/555.7 = **0.0738** —
5 % agreement. Over **coincident footprints** (single-gene single-shadow loci) the gDNA component's mean
efficiency and the shadow's agree to a median **1.031** ON and exactly **1.000** OFF, so the two
contractions are not mis-weighted against each other either. The junction-probe blindness is real and is
`ISSUES: ruler-witness-geometry-on-transcript-panels`; it is worth 13–28 % here and is a different
defect.

**The intronic pool is not the source under capture.** Only 2.76 % of gDNA's *contained* fragments sit in
intron-only regions at `g50 ss.99 ON` (70,288 against a leak of 452,854). The pool is the crossing
fragments — see §3.

**A per-locus prior cannot reach it, structurally.** The RNA prior enters every RNA component as the
same factor (`EQUATIONS.md` §9b), so it moves the gDNA:RNA split and not a within-RNA one.
`--arm oracle` reads 541,216 → 541,762.

**A stronger gDNA pseudocount is not a route.** Sweeping its strength on the shipped solver, closing the
channel needs a pseudocount of **1× the data at ratio 2, 5–10× at 6.2 and 50× at 20**. A prior many
times the data is not a prior.

**⛔ And one instrument was lying.** `quant_accuracy.truth_weights` read `observed_mrna_fragments`,
identically 0 on all 6,919 SYNTHETIC rows (the nascent truth is `observed_nrna_fragments`), so
`--arm oracle_alloc*` handed every shadow a weight of **zero** — the retired `alpha = 0` rule wearing an
oracle's name, tracking the pre-restoration baseline to within 15 % on every in-scope condition. Fixed
and gated (`d651515d`). Re-measured with the true weights, a perfect per-transcript allocation removes
67 % at `g50 ss.99 ON` and 52 % at `g98 ss.99 ON` — a large lever, not the cause, since 181,136 survives
a perfect allocation.

---

## 5. Where the information is: calibration already has it

This is the part that matters for the repair.

### 5.1 Calibration is right

The per-locus prior the EM receives, against the certified truth through the **same** assembler:

| condition | `gdna_prior` cal / true | err | `rna_prior` cal / true | err |
|---|---|---:|---|---:|
| `g05 ss99 ON` | 505,889 / 501,158 | +0.9 % | 5,381,719 / 5,384,115 | **−0.0 %** |
| `g50 ss99 OFF` | 2,380,225 / 2,404,082 | −1.0 % | 3,358,599 / 3,333,710 | +0.7 % |
| `g50 ss99 ON` | 4,954,189 / 4,937,801 | +0.3 % | 2,793,710 / 2,805,267 | −0.4 % |
| `g98 ss99 ON` | 9,551,684 / 9,614,034 | −0.6 % | 180,806 / 109,915 | **+64.5 %** |

Per object it is as good. At `g98 ss.99 ON`, by region class (contained fragments):

| class | regions | true gDNA | true mRNA | true nRNA | CAL RNA |
|---|---:|---:|---:|---:|---:|
| exon-bearing | 24,018 | 4,676,099 | 73,716 | 2,641 | 94,789 |
| **intron-only** | 9,805 | 137,888 | 0 | **92** | **5,951** |
| intergenic | 1,312 | 180,276 | 0 | 0 | 0 |

At the intron-only regions — the only place nascent RNA can sit — **calibration says 5,951 fragments
against a true 92, while the EM's shadows hold 476,270.**

### 5.2 And the EM discards it

On the 20 loci carrying the most shadow false positives:

| condition | gDNA: calibration / certified | the EM outputs | discarded | silent shadows hold |
|---|---|---:|---:|---:|
| `g98 ss99 ON` | 1,101,626 / 1,110,372 (−0.8 %) | 914,155 | **196,217** | 150,579 |
| `g50 ss99 ON` | 674,692 / 670,985 (+0.6 %) | 509,508 | **161,476** | 152,777 |

One locus, walked end to end (`g98 ss.99 ON`, locus 1049 — 9 genes, 180 regions):

* calibration is exact per region — e.g. region 27950, true gDNA 1,251 / mRNA 218 / nRNA 2, calibration
  1,263 / 208;
* the prior the EM receives is gDNA 31,152 (true 31,838) and RNA 1,722 (true 983);
* the one annotated transcript with real RNA reads **2,483 against a truth of 2,381**;
* and **eight shadow spans hold 11,187 fragments against a truth of zero.**

The EM's answer for the real transcript is excellent. The shadows are the entire error.

### 5.3 The lane that would carry it has no producer

`rna_prior_weight` is plumbed end to end — `pipeline.py` → `estimator.py` → `em_solver.cpp` — and
**nothing in `src/` fills it.** There is a parameter and no producer. So `apply_grouped_prior_update`
always takes the fallback `w_i = raw[i]`: the prior's allocation echoes the EM's own current belief and
therefore *cannot contradict it*. `nascent-gets-no-rna-prior` (CLOSED) left the lane free for exactly
this.

---

## 6. What a measured prior is worth, measured

### 6.1 The producer

Calibration publishes an RNA mass per region with gDNA already removed. `assemble_priors` then SUMS
those into one `rna_prior_count` per locus, which destroys the per-object resolution that says "no RNA
in this intron". Deconvolving them back onto the transcripts under the **same** opportunity model the EM
uses —

```
    r_o  ≈  Σ_t θ_t · a[t,o] / L_t ,    a[t,o] = t's taper-weighted bases in object o
```

— is a Poisson mixture with a closed EM and no free parameter. What it says about the shadow pool:

| condition | deconvolution | TRUE nascent | ratio | the EM | ratio |
|---|---:|---:|---:|---:|---:|
| `g05 ss99 ON` | 103,024 | 286,337 | 0.36× | 355,617 | 1.24× |
| `g50 ss99 OFF` | 997,068 | 1,013,538 | **0.98×** | 946,781 | 0.93× |
| `g50 ss99 ON` | 56,489 | 150,432 | 0.38× | 691,515 | **4.60×** |
| `g98 ss99 ON` | 9,744 | 5,989 | 1.63× | 596,171 | **99.5×** |

Calibration's own statement is within a factor of ~2.6 everywhere; the EM is out by up to 100×.

### 6.2 Filling the lane wholesale does not work

| arm | `g50 ss.99 ON` nascent Δ | siphon left | transcript Σ\|Δ\| |
|---|---:|---:|---:|
| base | +541,216 | 100 % | 5.21 % |
| deconvolution as the whole weight | +113,081 | 21 % | **53.37 %** |
| + certified spliced mass | +79,211 | 15 % | **52.86 %** |
| + `em.warm_start=prior` | +71,841 | 13 % | **53.90 %** |

The pool improves and the transcript table is destroyed. **The lane is a single static array, so filling
it reallocates the WHOLE RNA pseudocount** — 2,793,710 fragments at `g50 ss.99 ON`, as large as the
unspliced RNA itself — and a coverage-derived weight is a far worse isoform allocator than the EM's own
per-fragment likelihood. Adding the certified spliced mass (13,482 junctions, 45,609 (sj, transcript)
pairs, **0 synthetic holders**) does not fix it, so the damage is the isoform split and not the
shadow/annotated balance. The warm-start half contributes almost nothing: the allocation does the work.

### 6.3 The half that works, priced alone

A component that reaches an object **no other component's structure reaches** has an independently
measurable mass; one whose opportunity is wholly shared has none. That test selects **97.4 % of shadow
spans and 0 % of annotated transcripts**, measured on all four conditions. Correcting only the tested
components and leaving every other weight alone:

| `ss 0.99` | nascent Δ base | tested-only | siphon left | tx Σ\|Δ\| base | tested-only |
|---|---:|---:|---:|---:|---:|
| **`g50 ON`** | +541,216 | **+22,187** | **4 %** | 5.21 % | **5.85 %** |
| `g50 OFF` | −66,752 | −45,507 | 68 % | 2.50 % | **2.53 %** |
| `g98 ON` | +590,406 | +278,411 | 47 % | 32.99 % | 65.69 % |

**96 % of the siphon at the worst well-calibrated condition, for 0.64 points of transcript error, and no
harm off capture.**

⛔ Two caveats, both load-bearing:

* **This arm is a DIAGNOSTIC and cannot ship.** The untested half reads the base run's own counts, which
  is circular. A shippable version needs `raw[i]`, which lives in the kernel, not in a static lane.
* **`g98` is a cancelling defect pair.** Its RNA pseudocount (180,806) is nearly the whole true RNA
  (194,011) and over-states by 64.5 %, so the shadows had been acting as its **sink**. Removing the sink
  without fixing the over-call simply moves the error onto the transcript table (33.0 → 65.7 %).

---

## 7. The two candidate repairs

### Candidate A — a per-gene gDNA opportunity ⭐ recommended

**What.** Give the gDNA component a per-gene opportunity and pseudocount rather than one pooled over the
whole connected component, so `L_g/L_n` stays near 1 and the existing prior suffices.

**Why it should work, with the number.** On the shipped solver, at the ratio a per-gene opportunity
gives — **1.25**, since single-gene loci measure `span_g/fl_n` = 1.21 and the two contractions agree to
1.03 — the shadow takes **0.00 %** with the gDNA pseudocount the EM **already receives**, against 45–49 %
at the panel's mass-weighted 9.7:

| `L_g/L_n` | 1.00 | **1.25** | 1.50 | 2.00 | 2.50 | 9.70 |
|---|---:|---:|---:|---:|---:|---:|
| no prior | 0 | 17.4 % | 26.2 % | 35.2 % | 40.0 % | 59.3 % |
| `gdna_prior = 0.5× N` | 0 | **0.00 %** | **0.00 %** | 17.2 % | 26.0 % | 49.6 % |
| `gdna_prior = 1× N` | 0 | **0.00 %** | **0.00 %** | 0.00 % | 12.9 % | 45.1 % |

**Why I prefer it.** It needs **no new information and no new estimator**. It does not touch the isoform
allocation, so the transcript table is not put at risk by a weaker allocator. It removes the instability
rather than counterweighting it, so it does not depend on `rna_prior_count`'s accuracy — which is the
thing that breaks candidate B at `g98`. And it explains the measured gradient (§2.4) directly: 85.8 % of
the leak is in exactly the loci where the pooling is wrong.

**What it costs.** It changes `LocusPriors`, so calibration's own consumers move with it and all three
controls — `calibration_vs_oracle.py`, `zero_controls.py`, `policy_benchmark.py` — become live rather
than being the "must not move" guardrails they were this session. It is the expensive option. It also
needs a ruling on what "per gene" means where genes overlap, which is the same geometry question
`assemble_priors` already answers for regions by share.

**How I would build it.** Derive on paper first: the locus's objects already carry per-object gDNA mass
and support, and `_region_locus_shares` already apportions an object across loci by share — the same
apportionment across the genes within a locus is the whole change. Prototype in a worktree, A/B with
`quant_accuracy.py` per stratum above `--arm base_reseed` and the three controls beside it, on
`g50 ss.99 ON` and `g50 ss.99 OFF` first, then the full ladder.

### Candidate B — the measured per-transcript prior on tested components

**What.** Fill `rna_prior_weight` from the §6.1 deconvolution, but only for components with exclusive
objects; every other component keeps `raw[i]`.

**Worth.** 96 % of the siphon at `g50 ss.99 ON` for 0.64 points (§6.3).

**What it costs.** It cannot be done through the static lane — expressing "keep `raw[i]` where the data
cannot speak" needs `raw[i]`, so the shippable form is a small `em_solver.cpp` change plus a Python
producer for the deconvolution. And it does not help `g98` until `rna_prior_count`'s +64.5 % over-call
there is fixed.

**Why it is second.** It counterweights the instability instead of removing it, its benefit is bounded
by how well calibration localises RNA (0.36–1.63× across conditions), and it adds a new estimator to the
pre-EM chain.

### What I would NOT build

* **A length knob on the shadows** — it kills the live entities with the dead ones (§2.5).
* **A stronger gDNA pseudocount** — it needs 50× the data at the ratios that matter (§4).
* **An annotation-based null on synthetic entities** — that is the retired `alpha = 0` rule
  (`ISSUES: nascent-gets-no-rna-prior`, CLOSED by the owner as "a hack"), and it was masking this defect
  rather than fixing it. Any rule here must be about EVIDENCE, never about whether a span is
  manufactured (Axiom 0).

---

## 8. What is still unexplained

1. **The damping factor.** The bare two-component contest predicts 44–53 % of the shadow-exclusive pool
   at the measured geometry; the panel leaks 21.1 % (ON) and 11.6 % (OFF). The mature isoforms competing
   at exonic positions — pinned by spliced fragments the shadow can never claim — and the gDNA
   pseudocount are the two damping forces, and neither is closed quantitatively.
2. **The live shadows' under-call off capture** (−323,623 at `g50 ss.99 OFF`), which is what masks the
   false positive there. Not the same defect and not yet its own entry.
3. **Why a perfect allocation makes `g98` capture-ON markedly worse** (33.0 → 49.6 %, and 102.0 → 321.9 %
   on the deferred stratum). Both the broken and the repaired arm show it. §6.3's cancelling-pair reading
   is a hypothesis, not a measurement.
4. **`rna_prior_count`'s +64.5 % over-call at `g98`**, which is the reason candidate B stalls there. It
   is calibration-side and has its own instrument (`calibration_vs_oracle.py`, `prior_vs_oracle.py`).

---

## 9. Size it before paying for it

⚠ The whole defect is a competition between gDNA and a nascent pool the ladder runs at `on_fraction`
**0.50**, a DEVELOPMENT STRESS level (`DESIGN.md` §0b). Realistic is 0.10. `ISSUES:
nascent-stress-sensitivity` is therefore the **gating question, not a footnote**: re-simulate
`g50 ss.99 ON` at the realistic level and re-read the siphon before paying for a mechanism. A repair
worth 541,216 fragments at stress may be worth a fifth of that in the expected case.

---

## 10. How to reproduce every number here

* **The instruments**: `quant_accuracy.py` (the pool rows and `nrna_est`, per stratum above
  `--arm base_reseed`, `--set em.assignment_mode=fractional`), `calibration_oracle.py` (the certified
  per-object truth), `ruler_vs_truth.py`, and `calibration_vs_oracle.py` / `zero_controls.py` /
  `policy_benchmark.py` as the controls — all three confirmed unmoved across this session.
* **The gate that pins the threshold**: `tests/test_estimator.py`, the shadow-vs-gDNA rows.
* **The session's own harnesses**, kept at
  `~/Downloads/rigel_runs/prototypes/2026-09-19_siphon_root_cause/`:
  `investigation/` — the per-locus dissection (`dissect2.py`, eight `d2_*.npz`), the closed-form fixed
  point (`derive.py`), the solver toy (`shadow_toy.py`), the λ probe (`lambda_probe.py`), the corrected
  allocation arm (`alloc_fixed.py`);
  `investigation2/` — the per-region/boundary dissection (`dissect3.py`, four `d3_*.npz`), the object
  interrogation (`probe.py`), the locus walk (`walk.py`), the weight studies (`weights.py`,
  `weights2.py`, `deconv.py`, `lift.py`, `gate.py`), the sj map (`sjmap.py`), the prototype arms
  (`proto_weight.py`, `proto_weight2.py`, `proto_split.py`) and the anchor sweep (`anchor.py`).
