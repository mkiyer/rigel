# PLAN — the gDNA strand-overdispersion fit (derived 2026-08-29, owner-approved)

⚠ A dev doc: a plan and its working record, not a ruling. Findings that settle move to their home
(`TRAPS.md` as a named rule, `DESIGN.md`, `ROADMAP.md`) and are deleted here.

## 1. The defect, in one paragraph (measured on the test chromosome, cross-verified)

`gdna_strand.fit_gdna_strand_from_substrate` pools intron regions and exon|intron boundaries as seeds
with gDNA weight **w ≡ 1** (regions: `count_gdna_frac` is algebraically 1 for a count-observable
region; boundaries: `np.ones`, "gDNA by signature"). Nascent RNA in those seeds is strand-skewed on
stranded libraries and is read as gDNA strand overdispersion: od_g fits to **0.10–0.20 (the ceiling)
against a simulator truth of 0.0**; at κ = ½ the contamination is invisible (mixture mean is ½ either
way) so unstranded rows fit 0. The od term is `(n·f_ref)²·od/4` of the frozen strand variance — 95 %
of it at a 5 %-gDNA exon — inflating σ_f 8×; the wide likelihood is then pulled by ψ's reference
(`σ_f²/(2f)`) and truncated at the f = 1 vertex, giving the universal gDNA under-call on stranded rows.
κ̂ is unbiased. The λ grid is not involved. The fl model is not the pass-0 mechanism (bit-identical
under true pmfs). The pure-gDNA oracle partition has MoM od = 0.00000: the 0.113 is an artifact.

## 2. Derivation of the replacement (no constants introduced)

**Seeds.** A *structural* selection cannot truncate the observed overdispersion — it changes the
population, not the sample — so: (a) structurally pure objects (intergenic regions; exon|intergenic
gene-edge boundaries — certified zero RNA on every condition, and capture-ENRICHED because they abut
the probed exon) enter at w = 1; (b) every other count-observable seed enters at a **strand-free,
count-derived** weight — the density deconvolution's own posterior `P(g | C) ∝ NegBinom(g; ρ_bg·E_g,
α_eff)·1[g ≤ C]`, `w = E[g | C] / C`, against the pooled INTERGENIC background. This restores the
count ⊥ strand independence the estimator's docstring already invokes; the identity that made w ≡ 1
was using the seed's OWN density as its reference. ⛔ A data-driven filter on the strand split itself
("drop skewed seeds") conditions on the estimand and caps od — never write it.

**Shrinkage target.** Beta(14,14) ⇒ 0.0345 is abolished (conjured). The only defensible fallback is a
measurement from the same library: **od_r**, fitted from the per-SJ strand table (certified pure RNA),
weighted by od_r's OWN information — i.e. the prior "od_g is distributed like our estimate of od_r".
With clean seeds the prior is nearly inert (seed information ≫ any prior); it binds only where there
is ~no gDNA, where od_g has nothing to act on. The physical support [0, 0.2) is kept as a BOUND.

**Why a binomial (od := 0) is not the safe choice on real data.** od moves only the variance; with
genuine gDNA strand ICC a binomial is over-confident by `(n_g − 1)·od`, and at gDNA-rich slots skews
toward the RNA strand read as RNA while skews away clip at f = 1 → systematic gDNA under-call, at
inflated confidence. Beta-Binomial with od fitted well dominates in every world.

## 3. The plan — one mechanism per A/B, each on the test chromosome first, then the ladder

1. **Falsification test first** (`tests/calibration/test_gdna_strand_fit.py` when it lands; prototype
   copy beside the prototype): synthetic pure seeds (binomial ½) + contaminated intron-like seeds
   (sense 1 %, w ≡ 1 today) ⇒ shipped fitter returns od ≫ 0 (documents the defect, test FAILS on
   shipped code); new fitter returns ≈ 0 within its own null sd `1/√I`; PERTURB: force w ≡ 1 → fires;
   GUARD against the threshold trap: pure Beta-Binomial seeds at od = 0.05 ⇒ recovers 0.05 ± sd.
2. **A/B 1 — seeds + weights + od_r target** (prototype `scratchpad/odfit/proto_gdna_strand.py`,
   injected via `InjectedCalibrationPriors.gdna_strand_overdispersion`): silent and relay, all 30
   conditions, per stratum, both zero controls, vs shipped. Ceiling already priced (od := 0):
   silent ss0.99×ON 0.57×, ss0.99×OFF 0.86×, ss0.70×ON 0.51×, unstranded 1.00×.
3. **A/B 2 — functional form**: exact Beta-Binomial log-likelihood in place of the frozen-variance
   Gaussian in ψ's strand term, od_g held at its FITTED value (runtime patch of
   `simplex_logodds._mixture_strand_loglik` in the harness, never `src/`). Decides whether a CORRECT
   od > 0 on real data would still reproduce today's bias through the readout.
4. **True-od arm**: `test_reference_odg05.yaml` (gDNA strand ICC 0.05, same 30-condition cross) —
   binomial vs BB(fitted) vs BB(true), so the real-data claims in §2 are measured, not argued.
5. Only then `src/`: the fit lands in `gdna_strand` with weights computed in `calibrate` via
   `density_deconv` (layer 7 → 5 → 4, no upward import); the test moves to `tests/calibration/`.
6. **Next mechanism, same method:** the gDNA fl pmf's two intronic pools weighted the same way, and
   its EB anchor moved from the global MIXTURE pmf to the pure-gDNA pooled pmf (the fl side arms are
   the falsification substrate — the field gate fails on both gap arms today).

## 4. Open questions for the owner (each decided by an A/B above, not here)

* od_r as target with od_r's full information (a genuine blend when both are well measured) vs
  od_r as a FALLBACK only — A/B 1 reports raw pooled MoM beside the shrunk value per condition.
* Whether the 0.2 ceiling (Beta(2,2)) stays: it is a support bound, not a location, but it is chosen.
* The exact-zero branch sensitivity at the unstranded g00 control (8,873 → 10,573 FP) — a code-path
  question surfaced by the injection, to be located when the fit lands.

## 5. Working record — A/B 1 (2026-08-29, test chromosome, 30 conditions × silent/relay)

**Prototype v1** (intergenic background as the reference for EVERY impure seed) FAILED under capture:
a capture-enriched exon|intron boundary gets w ≈ 0.0005 against an 800×-too-low reference; a near-zero
weight does not drop the seed — it sets its mixture mean to κ and puts the whole residual in the
numerator over a zero denominator → od pinned at the ceiling, ss0.99×ON 8.1× WORSE. Lesson (a rule
candidate): **a seed with no trustworthy reference must be EXCLUDED (structural), never weighted ≈ 0.**

**Prototype v2** — region seeds referenced to the intergenic background; boundary seeds referenced to the
SAME LOCUS's gene-edge boundaries (Gamma posterior, Jeffreys ½; same probe chemistry); exclusion where a
locus has no pure boundary; od_r's own BetaBinom excess removed from the numerator. Fits od ≈ 0.000–0.005
on every truth-0 row (max raw 0.045 at g05 ss0.99 ON). Per stratum vs shipped (silent): ss0.99×ON 0.59×,
ss0.99×OFF 0.86×, ss0.70×ON 0.51×, ss0.70×OFF 0.98×, unstranded 1.00×, deferred 0.99×; relay 0.70× /
0.86× / 0.94×; relay's stranded capture-ON g00 inventions 2,999 → 49. No exclusions were needed on this
annotation (every locus has two gene-edge boundaries) — real annotations will exercise the exclusion.

**Shrinkage — DECIDED BY THE TRUE-od ARM (`test_reference_odg05.yaml`, gDNA strand ICC 0.05):**
od_r-information shrinkage OVER-shrinks where gDNA information is scarce (capture-ON g05/g25: 0.0005 vs
truth 0.05) and hides real overdispersion. **The raw pooled MoM clipped to the support, with od_r ONLY as
the no-seed fallback, tracks truth on both arms**: od05 arm ss0.99×ON raw 15,992 vs true-od 16,049
(shrunk 21,607; shipped 20,621; binomial 30,898); truth-0 arm ss0.99×ON raw 6,868 vs true 6,135 (shrunk
6,307). No location constant survives. Cost of no shrinkage: ~10 % noise at low-gDNA capture-ON rows.

**The binomial question, measured (od05 arm, silent):** od := 0 is 1.8× worse than the true od at
ss0.99×OFF (9,147 vs 5,189) and 1.9× at ss0.99×ON (30,898 vs 16,049) — over-confidence bites exactly
where derived. ⚠ BUT the TRUE od plugged into the frozen-variance Gaussian is WORSE than binomial at
ss0.70×ON (109,935 vs 72,175) and on the deferred stratum (320k vs 248k): an honestly wide likelihood is
pulled by ψ's reference measure (`σ_f²/(2f)`) and truncated at the vertex. **A/B 2 (functional form) is
therefore REQUIRED before the fit ships**: on real data a CORRECT od will still under-call wherever the
strand information is moderate. A/B 2 arms (runtime patches of `simplex_logodds._mixture_strand_loglik`
and the reference, never `src/`): (i) live-composition variance instead of the frozen reference variance;
(ii) a flat (uniform-in-f) reference against the Jeffreys f^{-½}; (iii) an exact Beta-Binomial marginal at
the sighted slots — each alone, at TRUE od on the od05 arm and at od ≈ 0 on the base arm.

**The g00-unstranded exact-zero sensitivity is NOT this fix's:** with κ injected at exactly ½ every od
gives 8,863 FP; with κ̂ = 0.5016 od = 0 gives 10,573. It is the fitted residue κ̂ − ½ carrying a spurious
slope whose weight od scales (`TRAPS: a-threshold-on-a-fitted-residue`) — the κ ≈ ½ zero-control defect.

## 6. Working record — A/B 2 and A/B 3 (2026-08-29): the functional form

Harness: runtime patches of `simplex_logodds._mixture_strand_loglik` / `_JEFFREYS_REF` (never `src/`),
silent, shipped K=3 unless stated, both arms (truth od 0 at its fitted od; od05 at od := 0.05).

* **Reference exponent** (½ shipped): 1 (uniform in f) is catastrophic everywhere (unstranded 6×, g00
  50–60k FP); 0 zeroes every g00 control and fixes the vertex but over-calls mixed exons at K=3 (+3,101
  at ss0.99×ON). Neither is a candidate. ⚠ The current λ-grid + ½ reference IS the arcsine-uniform
  (Jeffreys) prior; an arcsine grid would change only the vertex SUPPORT (f = 0, 1 as nodes), not the pull.
* **Exact BetaBinom⊕Binom marginal** (pass-0, four focus rows): 0.92–1.29× base, ≈1.00× od05; the vertex
  under-call barely moves (−1,340 → −1,174), mixed exons swing to over-call. Not the mechanism.
* **Dropping the −½ log var term**: blows up the zero controls (38k FP at g00 ss0.50/0.70). Out.
* ⭐ **LIVE-composition variance (variance evaluated at the solved f, log-var term KEPT)** — A/B 3 with the
  v2 od fit installed: truth-0 arm NEUTRAL (4,982→4,973; 6,868→7,107 = +3 %; deferred +6 %) because the od
  term vanishes; od05 arm LARGE WINS at K=3: deferred 320k→81k, ss0.70×ON 120k→34k, ss0.70×OFF 6,657→5,188,
  ss0.99×OFF 5,261→4,687, ss0.99×ON 15,992→14,283; g00 ss0.50 controls +27 % (9,436→11,955). Per-slot: the
  frozen `f_ref` never leaves ≈½ where the solve has no strand information, so the frozen variance is
  `(n/2)²·od/4` — 12× too wide at a 14 %-gDNA exon — and the shipped form COLLAPSES both a pure-gDNA exon
  (1,447 → 1.3) and a mixed exon (1,502 → 1.4) on the deferred stratum at true od 0.05.
* ⭐⭐ **The κ ≈ ½ win is NOT the κ̂ residue slope**: injecting κ := ½ exactly leaves it (88,264 → 15,135;
  mirrored residue 15,416). With the mean flat in f the live term is `−½r²/var(f) − ½ log var(f)`, a
  SECOND-MOMENT composition channel: gDNA's strand split wanders (od_g), RNA's is pinned at κ (od_r ≈ 0),
  so the strand residual's SIZE measures gDNA. The "count-zero-information freeze" suppressed exactly this
  term as "a manufactured preference toward the variance-minimum" — which it WAS while od_g was mis-fitted
  from contaminated seeds; with an honest od it is information (one χ²₁ draw per slot; noisy; Jensen-
  biased for small f; the only composition information the deferred stratum has). ⚠ Its validity on real
  data rests on od_g being fitted honestly and od_r being small — the v2 fit is the precondition.

**Decision proposed (owner's call):** land MECHANISM 1 (the v2 od fit: structural purity + strand-free
weights + raw clipped MoM, od_r fallback) in `src/` now — its A/B is clean on both arms. Treat the LIVE
variance as MECHANISM 2 with its own falsification test and ladder A/B; it is neutral where od ≈ 0 and
decisive where od > 0, and it re-opens a channel the freeze closed, so it deserves its own ruling.

## 7. Working record — the LADDER gate on mechanism 1 (2026-08-29)

`TRAPS: panel-before-src` earned its keep. Prototype **v2 on the 16-condition ladder**: stranded × OFF lands
on the truth ceiling exactly (silent 306,679 → 296,331 vs 296,310 at od := 0; relay 441,878 → 413,444 vs
413,312) and the stranded g00 controls the shipped fit poisons are repaired (silent 61,383 → 30,503;
relay capture-ON 43,983 → 20,625); unstranded rows untouched. ⛔ **But stranded × capture-ON REGRESSED
badly at g50/g98**: v2 fitted od ≈ 0.11–0.13 where the shipped fit happened to land ≈ 0 (silent 249k →
598k at g50, 276k → 1,032k at g98). The twin block — every locus one gene with two clean gene edges —
could not expose it (`TRAPS: a-toy-and-a-panel-can-disagree-in-rank`).

**Decomposition with certified truth per seed (g50/g98 ss0.99 ON):** 99.9 % of the excess variance sits
on the LOCUS-REFERENCED exon|intron boundary seeds, which truth says are 97–99.9 % gDNA — class "od" 1.28,
impossible for an ICC, i.e. a wrong mixture MEAN from wrong weights. On a real annotation under capture a
merged locus mixes probed and unprobed exons, so gene-edge and interior boundaries do not share an
enrichment class and the locus-local reference is untrustworthy. The pure classes alone — 1,177 intergenic
regions + 1,671 certified gene-edge boundaries (176k crossings; my structural mask matched the certified
`B gene edge` stratum 100 %) — fit od 0.0000–0.0001; the intron-REGION weighting is sound (class od 0.0005
at 29 % true RNA, capture-OFF). **v3 = pure seeds (w = 1) + count-weighted intron REGIONS against the
intergenic background; exon|intron boundaries EXCLUDED from the od fit** — the v1 lesson applied
consistently: a seed with no trustworthy reference is excluded, never weighted. v3 on the test chromosome
is indistinguishable from v2 (raw = truth ceiling at ss0.99×ON, 6,135 vs 6,135; od05 arm 15,556 vs
16,049). Ladder v3 verdict: pending at the time of writing (§8 when it lands).

## 8. The LADDER verdict on v3 (2026-08-29) — mechanism 1 clears both gates

v3 fits od 0.0000–0.0006 on every ladder row (truth 0; shipped fit 0.2000 at g00/g05 ss0.99 OFF, 0.1551 at
g05 ss0.99 ON). Per stratum, shipped → v3 (od := 0 ceiling in brackets): stranded × OFF silent 306,679 →
296,329 (296,310), relay 441,878 → 413,418 (413,312); stranded × ON silent 638,844 → 611,107 (620,710),
relay 961,484 → 910,909 (901,130); unstranded and deferred strata unchanged to <0.01 %. Zero controls: the
stranded g00 rows the shipped fit poisons are repaired (silent 61,383 → 30,503 and 34,538 → 20,810; relay
capture-ON 43,983 → 20,625); unstranded g00 within +0.2 % (the κ̂ − ½ residue slope at od exactly 0).
The largest single-row gain is g05 ss0.99 ON (114,173 → 85,724). ⚠ The ladder's improvement is MODEST
(3–5 % per in-scope stranded stratum) because its deep pure-gDNA seeds already held the shipped fit near
0 at g50/g98; the test chromosome's stress nascent made the same defect 2–3× there. No regression anywhere.

**Final shape of mechanism 1 (for the owner's go):** seeds = structurally pure objects at w = 1
(intergenic regions; exon|intergenic gene-edge boundaries) + count-observable intron REGIONS at the
density-deconvolution posterior weight against the intergenic background; exon|intron boundaries excluded;
pooled MoM with od_r's own excess removed; od = clip(raw MoM, 0, ceiling); od_r as the fallback when no
gDNA seed carries pairs. No location prior; no new constant. Lands in `gdna_strand` with the weights
computed in `calibrate` via `density_deconv` (layer 7 → 5 → 4); falsification test to `tests/calibration/`
with a structural case asserting no exon|intron boundary enters the seed set.

## 9. Mechanism 2 on the LADDER (2026-08-29) — the live-composition variance, with the v3 od fit installed

Truth od is 0 on the ladder, so only a COST could show — and none did: the live variance is a small WIN on
both in-scope stranded strata (silent, K=3: stranded × OFF 296,329 → 291,860 = 0.985×; stranded × ON
611,107 → 573,686 = 0.939×, the exon|intron boundary class −158,575 → −124,948), neutral on unstranded and
deferred (< 0.01 %), zero controls: g00 ss0.99 OFF +3.6 % (30,503 → 31,589), ON +16 % (20,810 → 24,166),
unstranded g00 ≤ +0.3 %. ⭐ Why there is a gain even at od ≈ 0: the freeze also evaluates the BINOMIAL
term at f_ref — at κ = 0.01 and f_ref ≈ ½ that is `n·0.255·0.745 ≈ 0.19n` against `0.033n` at f = 0.05,
5.8× too wide — so the frozen strand channel is under-weighted on stranded data independently of od.
Combined with the true-od arm (§6) the case for mechanism 2 is: no measured cost on any substrate beyond a
zero-control tick of +4–16 % on stranded g00 rows, wins on every in-scope stranded stratum, and a large
win wherever gDNA strand ICC is real. It still reverses a documented design choice (the
"count-zero-information freeze") whose derivation no longer exists in the docs — an owner ruling.

## 10. LANDED (2026-08-29, owner-approved) — mechanism 1 in `src/`

* `density_deconv.gdna_share_posterior` (layer 5): the strand-free seed weight, `E[g | C]/C` under the
  intron factory's own NegBinom against the intergenic background.
* `gdna_strand.fit_gdna_strand_overdispersion`: raw pooled MoM clipped to the support, the RNA
  component's excess removed, `fallback_overdispersion` handed in; `_region_seeds` takes the weight from
  the caller; `fit_gdna_strand_from_substrate(..., region_gdna_weight, rna_strand_overdispersion,
  fallback_overdispersion)`. The Beta(14,14) target survives ONLY in the RNA fit (a separate decision).
* `strand_deconv.boundary_seeds`: gene-edge (a flank intergenic) boundaries only.
* `calibrate`: RNA fit first; contained intergenic background → weights for count-observable
  non-intergenic regions → the gDNA fit with od_r as correction and fallback.
* Dead surface deleted: `RegionDensity.count_gdna_frac` (no consumer left) and its test.
* Gate `tests/calibration/test_gdna_strand_fit.py` (10 cases): watched FAILING on the shipped code (5/10 —
  the two selector gates on assertion, three on the missing contract); PERTURBED after landing (drop the
  gene-edge rule → 2 fire; force w ≡ 1 → 1 fires); restored → 10/10. Old gates rewritten to the new
  contract (fallback = handed od_r; no shrinkage; gene-edge fixtures).
* Suite 3,656 / 8 xfail / 0 failed; goldens updated with diffs ≤ 1.3e-3 relative recorded first.
  Shipped-path reproduction: test panel silent ss0.99×ON 6,136 (prototype 6,135), ×OFF 3,795 (3,794);
  ladder rows identical to §8. ⚠ Accepted, documented consequence: at the unstranded g00 control the
  fit now returns the od_r fallback (0.0000) and silent's FP read 10,573 (was 8,873) — the κ̂ − ½ residue
  slope (§5), a separate defect.

## 11. Mechanism 2 split — 2A is NOTHING, all of it is 2B (2026-08-29)

Owner question: is the second-moment channel real, or over-engineering? A hybrid arm — the BINOMIAL variance
evaluated at the solved composition, the ρ (ICC) terms still frozen at f_ref ("2A") — against frozen and
fully live, both test arms: **2A is neutral everywhere** (truth-0 arm ss0.99×ON 6,135 → 6,126; od05 arm
16,049 → 16,327, 109,935 → 109,222; every other stratum within ±1 %). So the ladder's 0.985× / 0.939× and
every od05 gain belong to the ρ term evaluated live — at ladder depths even ρ ≈ 10⁻⁴ frozen at f_ref = ½
dominates the binomial variance ((n/2)²·ρ/4 vs n/4). There is no premise-free half to land. Whether 2B is
information or a trap depends entirely on whether real gDNA strand ICC exists and is honestly estimated:
on the simulator it was put there by fiat (od05) or absent (truth 0). ⛔ Not adopted. The evidence to gather
is on the real libraries in `~/Downloads/rigel_runs/cfrna` (four, incl. the vcap RNA/gDNA spike-in): the
landed fit's ρ_g and ρ_r with their information, and the fit-free check — at structurally pure objects, does
the strand split's excess variance over Binomial(½) scale as n(n−1), the ICC signature?

## 12. REAL DATA (2026-08-29) — four libraries in `~/Downloads/rigel_runs/cfrna`, whole-genome index built from `refs/`

Pure seeds = intergenic regions + gene-edge boundaries (structurally pure on the annotation). od by MoM ± null sd.

| library | frags | κ̂ | od_r | LANDED od_g | OLD raw | pure, all | pure < 20 | pure ≥ 20 | 1st-moment RNA share, ≥ 50 |
|---|---|---|---|---|---|---|---|---|---|
| vcap RNA20M+gDNA5M | 18.0 M | 0.0002 | 0.016 | **0.2000** | 0.114 | 0.612 | **0.020 ± 0.001** | 0.689 | **+0.40** |
| LBX0588 cfRNA (heavily gDNA) | 1.23 M | 0.064 | 0.070 | **0.174** | 0.017 | 0.018 | **0.006 ± 0.002** | 0.023 | −0.06 |
| MO_3021 cfRNA | 0.74 M | 0.0023 | 0.008 | **0.2000** | 0.652 | 0.166 | **0.076 ± 0.003** (flat 0.076/0.076/0.075 over 2–20) | 0.213 | 0.00 (gene-edge mean sense 0.34) |
| LBX0190 cfRNA (sparse) | 0.13 M | 0.0025 | 0.007 | 0.2000 | 0.931 | 0.894 | 0.83 ± 0.02 | — | — |

**Three findings, each changing a recommendation.**

1. ⛔ **The weighted-intron seeds FAIL on real data.** Real introns are 75–90 % RNA by the count background
   (weight medians 0.11–0.24), and at small w the mixture mean is wrong (intronic RNA's sense fraction is
   not the spliced κ; the intergenic rate is not the intron's): the intron class reads od 8–351 and drives
   the landed fit to the ceiling on 3/4 libraries — and to 0.174 on LBX0588 where the OLD fit read 0.017.
   On the simulator the weighting was inert (§7: class od 0.0005) because nascent sat exactly at κ over
   exactly-uniform gDNA. **Proposed: pure seeds only** — identical on the simulator, strictly better on real
   data (LBX0588 0.018). `TRAPS: real-data-is-a-test-input`, paid.
2. ⛔ **Even structurally pure seeds are contaminated at DEPTH on real data** — unannotated transcription:
   vcap pure seeds ≥ 50 fragments read od 0.81 with a first-moment RNA share of +40 %; MO_3021's gene-edge
   boundaries have mean sense 0.34 (RNA-skewed, 1,168 seeds). Annotation incompleteness is the real-data
   contaminant, exactly as the owner predicted. Shallow pure seeds are clean; the honest gDNA strand ICC is
   **0.006 (LBX0588), 0.020 (vcap), 0.076 (MO_3021)** — MO_3021's is flat across depth bins, the ICC
   signature, so real gDNA strand ICC EXISTS and varies by library. Neither the old nor the landed fit is
   trustworthy on 3/4 libraries; both sit at the ceiling. A robust estimator against contaminated deep seeds
   is the next derived mechanism (candidates: the first-moment contamination share per depth stratum as the
   mixture weight — it only sees NET skew; an information-weighted median over depth strata; never a depth
   constant). The test chromosome needs an UNANNOTATED expressed transcript to reproduce this failure mode.
3. **2B (the second-moment channel): not adopted.** Real ρ_g is 0.006–0.08 and od_r 0.007–0.07 — comparable,
   with the RNA term muted only by κ(1−κ). What the fit reads at depth on real data is unannotated RNA, not
   ICC: the freeze's caution was right on real data. 2A alone was neutral (§11). No part of mechanism 2 ships.

## 13. DERIVED — the AWAY-HALF estimator: overdispersion without purity (2026-08-29, owner direction)

**Why purity cannot be the premise (owner):** "pure" is a property of an annotation AND a sample. Pervasive
transcription is real; the intergenic space is whatever the user's GTF leaves over (a protein-coding-only
reference makes most of the genome "intergenic"); references grow; and most genes are OFF in any one sample
— perfectly good gDNA seeds — but nobody knows which. An estimator that trusts a structural class is an
estimator that changes with the GTF.

**Model.** Any count- and strand-observable object of a GENE (intron region, exon|intron boundary, gene-edge
boundary; sense = the gene's strand) with N fragments and K on the sense strand holds an unknown gDNA share
w ∈ [0,1]: K ~ BetaBinom(N, μ(w), ·), μ(w) = ½·w + κ·(1 − w). The estimand is ρ_g at w = 1; the distribution
of w over seeds (which genes are on, nascent, the GTF) is nuisance.

**Lemma (one-sidedness).** On a strand-specific library RNA of the seed's own gene can only pull its sense
fraction TOWARD κ. Orient d = K − N/2 so that RNA pulls d negative. Under w = 1, d is symmetric about 0 and
the moment excess (d² − N/4) is an EVEN function of d. Hence over the AWAY half — seeds with d > 0, ties at
weight ½ — the expected excess is exactly half the pure population's, and

    ρ̂_g = Σ_away (d² − N/4) / Σ_away N(N−1)/4

is unbiased for ρ_g for EVERY distribution of w: a contaminated seed reaches the away side only by noise,
with a small d, biasing ρ̂ DOWN, never up. No weights, no purity class, no threshold; half the information.
The seed set is every genic object in whatever GTF is supplied — it GROWS with annotation. Intergenic seeds
have no gene strand and get no protection: they leave the od fit (they remain the density reference).

**Evidence.** Synthetic (3,000 pure BB seeds, ρ = 0.05): pure 0.0513; + 30 % same-strand RNA w~U(0,1):
full MoM 0.113, away 0.0509; + 30 % nascent-dominated (w = 0.1): 0.220 vs 0.0513; + 60 % contaminated
(a MAJORITY): 0.161 vs 0.0496; + 10 % ANTISENSE: 0.075 vs 0.0965 (the known blind spot). Simulator, ALL
genic seeds, no purity: truth 0 → 0.0001 / 0.0014 / −0.0001 (full MoM 0.176 / 0.017 / 0.017); truth 0.05 →
0.037 / 0.040 / 0.056. Real data (pooled genic): vcap 0.114 → **0.023** (shallow pure seeds read 0.020),
LBX0588 0.017 → 0.019 (flat across depth), MO_3021 0.65 → 0.48, LBX0190 0.93 → 0.96.

**What the coherence test says about the two libraries that still read ~0.5–0.9.** At N = 2 a binomial puts
25 % of seeds at K = N and 25 % at K = 0. LBX0588: 0.25/0.26 (binomial; ρ ≈ 0.01). vcap and MO_3021 genic:
K = 0 at 0.32/0.39 vs K = N at 0.21/0.20 — ONE-sided, same-strand RNA, removed by the away half; their
residual at depth (vcap ≥ 100 frags: 0.16–0.32; MO_3021 ≥ 20: 0.35) is away-side coherence in a few deep
seeds — antisense / unannotated overlapping transcription, the lemma's blind spot, amplified by the MoM's
N(N−1) leverage. LBX0190: 0.37/0.40 on BOTH sides, genic AND intergenic alike — residual PCR-duplicate
families: a GENUINE technical strand ICC near 0.9 that the Beta-Binomial must carry (a binomial there would be
catastrophically over-confident). Real gDNA strand ICC therefore ranges 0.01 → 0.9 across libraries.

**Two residual weaknesses, each measured, neither a constant:** (i) mixed shallow seeds crossing to the away
side bias ρ̂ DOWN (vcap N 2–5: −0.018; od05 capture-OFF 0.037 vs 0.05) — conservative, and de-weighted by the
information weighting; (ii) antisense / unannotated transcription inflates the away side at depth, and the
pooled MoM's quadratic leverage lets a few deep seeds dominate (MO_3021). Candidate remedies to derive next:
per-depth-stratum estimates combined by null information with each stratum clamped to the model's OWN support
(no new constant), or representing unannotated antisense as the RNA− population it is (AXIOM 0) rather than
as noise. ⛔ Not a depth cutoff.

**Proposed correction to mechanism 1 (owner's call):** replace the seed selection + count weights in `src/`
with the away-half estimator over all genic count/strand-observable objects; intergenic seeds out of the od
fit; od_r stays the no-pair fallback and the numerator correction. Identical to the od := truth ceiling on
the simulator; reference-agnostic; strictly better than both the old and the landed fit on real data.

## 14. BUILT — the BLANK CHROMOSOME and SHADOW transcripts (owner design, 2026-08-29; my authority for the layout)

* **Simulator mechanism**: `WholeGenomeSimConfig.shadow_gtf` — a supplemental GTF the simulator draws
  fragments from and the index NEVER sees. `whole_genome.merge_shadow_transcripts` appends them after the
  index's rows (no nascent; a shadow id the index knows is REFUSED); their fragments are named like any RNA
  fragment, so the oracle files them as `mrna` with no read-name change. Manifest records `shadow_gtf` and
  `n_shadow_transcripts`. Gate: `tests/test_sim_shadow_transcripts.py`.
* **Reference**: `build_test_reference.py` now writes TWO contigs into `test_chr.fa` — `test_chr` (annotated)
  and `test_blank` (1 Mb, own seed, NO annotation) — and injects splice motifs for shadow introns on the blank
  genome; `test_shadow.gtf` is the fourth hand-edited file (checked: shadows on `test_blank` only, ids unique
  across both GTFs, abundance rows over the union). Self-test 23/23. The index built from `test_chr.gtf`
  gives `test_blank` exactly ONE intergenic region [0, 1,000,000) — region 151 — so the control is not vacuous.
* **Layout** (`test_shadow.gtf`): 8 shadows at 100 kb spacing, strands alternating: SE 2 kb (250, 75, 3, 1
  molar), SE 10 kb (25), SE 500 bp (250 — short and deep), two 3-exon (25, 8); ≈ 5 % of the library's RNA
  fragments — pervasive transcription at a realistic dose with two DEEP shadows reproducing the real-data
  failure mode. No probes (an annotation-designed panel cannot target what it does not know). `gdna.genomic_refs`
  includes `test_blank` in all five panel configs; all five panels regenerated against the new index.
* **What the control asserts**: on region 151 the annotation says pure gDNA and the truth says gDNA + shadow
  RNA. Any od estimator that trusts intergenic purity (the OLD premise; PURE-only) must move; the AWAY-HALF
  over genic seeds (no intergenic, no purity) must be UNMOVED. Driver: `scratchpad/control/blank_control.py`.

## 15. THE CONTROL PASSES (2026-08-29) — `scratchpad/control/blank_control.py`, main panel, 30 conditions, truth od 0

Region 151 (`test_blank`, the annotation's "pure gDNA") received the shadows: 8.9–9.2k RNA fragments at g00
(100 % RNA), 63 % RNA at g05 OFF, 21 % at g25, 8 % at g50, 0.2 % at g98; ~20 fragments under capture (unprobed).
| estimator | g00 ss0.99 OFF | g05 ss0.99 OFF | g25 ss0.99 OFF | g25 ss0.70 OFF |
|---|---|---|---|---|
| OLD premise (all seeds, w ≡ 1, intergenic in) | 0.2000 | 0.2000 | 0.0914 | 0.0153 |
| LANDED fit (pure seeds + weighted introns) | 0.2000 | 0.1130 | 0.0113 | 0.0021 |
| PURE-only (intergenic + gene edge) | 0.2000 | 0.1130 | 0.0113 | 0.0022 |
| **AWAY-HALF, genic seeds, no intergenic** | (no gDNA) | **−0.008 ± 0.017** | **−0.003 ± 0.003** | **−0.002 ± 0.003** |
Every estimator that trusts intergenic purity moves; the away-half is UNMOVED on every contaminated row
(|ρ̂| ≤ 0.008, ≤ ~1 sd). Exactly the owner's point: "pure" is a property of the annotation, and unannotated
transcription is real. ⭐ Cleared to land as the correction to mechanism 1.

## 16. LANDED — the AWAY-HALF estimator replaces mechanism 1's seed weights (2026-08-29)

* `gdna_strand.away_half_moment(sense, total, kappa) → (od_mom, information)`: the core; `d = (K − N/2)·sign(½ − κ)`,
  weight `1[d>0] + ½·1[d=0]`, information = ½ × pair-count null information of the entering seeds.
* `fit_gdna_strand_overdispersion(sense, total, rna_sense_frac, *, fallback_overdispersion)`: raw clipped away-half
  moment; NO weights, NO od_r numerator correction (nothing to correct — RNA cannot reach the away side);
  fallback = the handed od_r when no pair is in the away half.
* Seeds = every count- and strand-observable GENIC object: `_region_seeds` → intron regions of TS_POS/TS_NEG
  genes (intergenic OUT, AMBIG OUT); `strand_deconv.boundary_seeds` → every strand-observable contiguous boundary,
  exon|intron AND gene-edge (returns `(sense, total)`, the weight column is gone).
* `calibrate`: the intergenic-background fit for weights and `density_deconv.gdna_share_posterior` DELETED (no
  consumer). `density_model` docstring corrected.
* Gates `tests/calibration/test_gdna_strand_fit.py` rewritten (13 cases: the lemma at 30 %/60 % contamination and
  nascent-dominated w = 0.1; two-sided moment as the PERTURBATION; truth-0 under contamination; antisense as the
  recorded limit; κ ↔ 1−κ symmetry; exact N = 2 tie weight; fallback; the three selectors). PERTURBED after
  landing — two-sided moment → 6 fire; intergenic readmitted → 1; tie at full OR zero weight → 1 (a gate was
  ADDED for this: the first tie perturbation fired NOTHING); orientation ignoring κ → 4. `test_gdna_strand.py`
  moved to the API (uniform contamination never INFLATES; wrapper fixture intron+ instead of intergenic);
  `test_gdna_strand_integration.py` fixture intron+ (intergenic is not a seed).
* Goldens: two files, `em_effective_length` only, ≤ 2.4e-5 relative, recorded then updated.
* Suite: 3,672 collected = 3,664 + 8 xfail (shadow test file +5, fit gates +3) — pending the confirming run.

**TEST PANEL A/B (main panel WITH shadows; `scratchpad/awayhalf/ab_inject.py`; ship = away-half; landed = the §10
fit's od per condition injected; true = od 0):**
| stratum | silent ship / landed / true | relay ship / landed / true |
|---|---|---|
| ss 0.50 × off | 23,362 / 23,365 / 23,362 | 24,255 / 24,248 / 24,256 |
| ss 0.50 × on (DEFERRED) | **296,779** / 324,961 / 296,759 | 24,728 / 24,681 / 24,681 |
| ss 0.70 × off | 23,697 / 23,695 / 23,697 | 25,317 / 25,330 / 25,317 |
| ss 0.70 × on | 13,161 / 12,991 / 12,991 | 17,944 / 17,838 / 17,838 |
| ss 0.99 × off | **22,727** / 22,842 / 22,711 | **24,559** / 24,715 / 24,484 |
| ss 0.99 × on | 5,919 / 5,820 / 5,800 | 11,510 / 11,416 / 11,370 |
Ship ≈ true on every stratum; where the landed fit was moved by the shadows (ss 0.99 × OFF: od 0.113 at g05,
0.0113 at g25) ship is better; zero controls: ship = true, landed poisoned relay at g00 ss0.99 ON (3,757 vs 134).
⚠ HONEST COST: under capture on this small panel (130 seeds, few pairs) the away-half is NOISIER — od 0.008 ± 0.010
at g05 ON — costing +1.3 % (ss 0.70 × ON) and +1.7 % (ss 0.99 × ON) against the landed fit. Half the pairs is the
price of needing no purity; on genome-scale data (ladder: 5–14 k seeds) the variance is negligible.

**LADDER (16 rows, `ab_inject.py`, landed = §8's v3 values):** shipped od 0.0000 on every contaminated row
(0.0013 at g00 ss0.50 ON, 0.0003 at g05 ss0.50 ON — noise on 5.5–12.7 k seeds). Per stratum silent ship / landed /
true: ss0.50×off 371,017 / 371,017 / 371,017; deferred 18,612,500 / 18,610,397 / 18,610,397; ss0.99×off 296,396 /
296,310 / 296,310; ss0.99×on 621,247 / 622,279 / 620,710. Relay likewise (ss0.99×on 904,881 / 909,879 / 901,130).
No regression; at 12–27 k seeds the half-pairs variance is invisible.

**REAL DATA (shipped path, `realdata_od.py`):** vcap 0.114 (old) → **0.0201**; LBX0588 0.174 (landed) → **0.0165**
(its pair-level truth 0.016); LBX0190 → 0.2000 (ceiling), MO_3021 → 0.2000 (ceiling) — the two κ̂ ≈ 0.002 libraries
whose pure seeds already showed ICC 0.55–0.9 / 0.08–0.45 at the SHALLOWEST depths. Raw-value check below.

**THE TWO CEILING LIBRARIES ARE NOT CONTAMINATION (`realdata_raw.py`, raw unclipped values):**
| library | away-half (all genic) | full moment | intergenic (excluded) | by depth [2,5) / [5,20) / [20,100) / 100+ |
|---|---|---|---|---|
| LBX0190 (130 k frags) | 0.974 ± 0.007 | 0.938 | 0.899 | 0.69 / 0.89 / 0.97 / 0.98 |
| MO_3021 | 0.675 ± 0.003 | 0.763 | 0.213 | **0.03** / 0.21 / 0.46 / 0.87 |
| vcap | 0.0201 ± 0.0003 | 0.093 | 0.182 | −0.02 / −0.003 / 0.01 / 0.17 |
LBX0190 is a REAL near-1 ICC on every class and every depth (mean sense 0.496 — symmetric, not RNA-pulled): each
seed's fragments sit on one strand — the duplicate-family library. The 0.2 ceiling clips a real 0.97; whether
`Beta(2,2)` is the right ceiling is an OWNER question (it is an asserted constant; a real library exceeds it 5×).
MO_3021 is DEPTH-GRADED: ≈0 at shallow seeds (where the full moment still reads 0.19 — contamination removed, as
designed) rising to 0.87 on the 8 deepest seeds; the pooled moment weights a seed by its PAIRS (∝ n²), so a
handful of deep seeds carry it to the ceiling. Two readings remain possible for the deep seeds — duplicate
clumping that grows with depth, or unannotated ANTISENSE at expressed loci (the recorded limit) — and this
instrument cannot separate them. ⚠ Open, not fixed: no mechanism proposed; brought to the owner.
vcap's 0.0201 is carried by its 394 deepest seeds (0.17) over ~350 k shallow ones at ≈ 0 — same weighting note.

**THE 0.2 CEILING IS LOAD-BEARING ON REAL DATA (`scratchpad/awayhalf/ceiling_price.txt`, gDNA mass on chain slots):**
| library | od 0.0 | od 0.2 (the ceiling = what SHIPS) | od = its own raw away-half |
|---|---|---|---|
| LBX0190 (raw 0.974) | 8,236 = 6.3 % | **7,130 = 5.5 %** | 4,572 = **3.5 %** |
| MO_3021 (raw 0.675) | 79,561 = 10.7 % | **74,286 = 10.0 %** | 69,893 = 9.4 % |
Direction confirms the original defect analysis: a HIGHER od_g widens the gDNA strand variance, the `−½·log var`
term penalises gDNA, and LESS gDNA is called. So on LBX0190 the ceiling changes the whole-library gDNA call by
**1.56×** (5.5 % vs 3.5 %) — `_MAX_OVERDISPERSION` is an ASSERTED constant that is NOT inert on real data (unlike
the RNA fit's shrinkage weight, measured inert). ⛔ There is no truth for a real library, so this prices the
constant; it does not adjudicate it. OWNER QUESTION: does the `Beta(2,2)` ceiling stay, and is a duplicate-family
library's true ICC ≈ 0.97 a thing the strand model should be allowed to believe?

**INSTRUMENTS AFTER THE `src/` DELETION:** `preflight.py --full` — every `scripts/design/` instrument imports
(73 files) and every `--self-test` passes (17/17); both references, both panels, all five oracle partitions and
all 46 certified `slot_truth` present. No ✘. (Required by CLAUDE.md after a `src/` deletion — `gdna_share_posterior`
and the `boundary_seeds` weight column went — because a green suite has hidden dead instruments before.)

## 17. ADVERSARIAL REVIEW OF THE LANDING, AND THE THREE REPAIRS (2026-08-30)

A 92-agent review of the §16 landing; 67 agents died on the session limit, so only findings marked
`[0/3 refuted]` were actually verified — everything `[0/0]` is UNVERIFIED and is NOT cleared.

**VERIFIED AND REPAIRED (each re-derived by hand before acting):**
1. ⛔ **`away_half_moment`'s information halved TWICE** — `0.5 * _null_information(n[weight > 0])` halves a
   set that is already the away half. Correct: `I = P/2` over the TOTAL pair count (`Var(e_s)|₀ = n(n−1)/8`,
   `E[a_s] = ½` ⇒ `Var(num) = P/8`, `E[den] = P/4`, `Var(od_mom) = 2/P`). Measured 1.41× too wide against a
   Monte-Carlo null. ⚠ Production DISCARDS the value, so no shipped number was ever wrong — but it is
   exported and documented, and every `± sd` quoted in §15/§16 is ~1.4× too wide (conclusions unchanged: the
   away-half readings there are ≤ 1 sd either way).
2. ⛔ **`κ = ½` EXACTLY collapsed the fit to a hard `od = 0`** with `fallback_used` False. Reachable:
   `κ = (n_same+1)/(n_obs+2)` is exactly ½ whenever `2·n_same = n_obs` — the modal outcome on an unstranded
   library. Repaired with the FULL two-sided moment at that κ (RNA is symmetric there, so the one-sided
   guarantee survives); no new constant.
3. ⛔ **`_region_seeds`' NEG-strand column flip was load-bearing and UNPINNED** — the whole 1,268-case
   calibration suite passed with it deleted. Without it a minus-strand gene's residual is oriented backwards
   and contamination INFLATES od (measured: 0.000 → 0.153 on a contaminated fixture).
Plus three docstrings the change had made false (`strand_deconv`'s 3-tuple; `gdna_strand`'s "symmetric by
design" / "prior overdispersion" prose; `density_model`'s claim that its density feeds the fit).

**GATES: 13 → 16.** The tie gate had asserted the implementation's own info expression — so it certified
defect 1 and fired on the repair; it is now a property gate (closed form + Monte-Carlo null se). New:
`test_the_information_is_HALF_THE_TOTAL_PAIR_COUNT_not_half_the_away_halfs`,
`test_an_unstranded_library_uses_the_FULL_moment_and_is_not_forced_to_zero`,
`test_a_NEG_strand_genes_seed_orients_to_the_NEG_column`. **Seven perturbations, all firing, each new gate on
exactly its own defect**: revert the info fix → 1; delete the κ=½ branch → 1; delete the NEG flip → 1;
two-sided moment → 8; readmit intergenic → 1; ties at full weight → 1; orientation ignores κ → 5.

**Suite 3,667 / 8 xfail / 0 failed** (re-derived). ⭐ **The test panel is BYTE-IDENTICAL to the pre-repair
run** — all three were latent (κ ≠ ½ exactly on every panel condition; the information is discarded in
production), so this is a correctness repair with no measured movement.

**UNVERIFIED, REPORTED NOT ACTED ON:** `RegionDensity.density` (and `run_fill`, its only caller) has had NO
production consumer since the weight deletion — grep finds reads only in `tests/calibration/test_density_model.py`.
Deleting them is a separate mechanism and an owner call; the docstring now says so instead of claiming a
consumer that does not exist.

## 18. THE DEAD DENSITY MACHINERY IS DELETED (owner ruling, 2026-08-30)

`RegionDensity`, `region_gdna_density` (~65 lines of flank anchoring → run-fill → global baseline) and
`run_fill.py` existed to WEIGHT the strand fit's seeds. The away-half moment needs no weight, so they had no
production consumer; only `count_observable_masks` survives, and `calibrate` now calls it directly — the fit
takes the two boolean masks, not a density object (`fit_gdna_strand_from_substrate(..., region_count_observable=,
boundary_count_observable=)`, `boundary_seeds(substrate, region_arrays, boundary_count_observable)`).
`test_density_model.py` was rewritten to gate the surviving function (4 cases: the region rule, the boundary
SHARING rule, exon+|exon− sharing no bit, and the per-reference axes); `test_run_fill.py` deleted with its module.
Layer-5 tables in `CLAUDE.md` and `DESIGN.md` updated. **Suite 3,655** (3,667 −3 module −8 test-file −1 case,
exact); ⭐ test panel BYTE-IDENTICAL; module census: no upward imports, no dead public surface; preflight 17/17.

## 19. THE OUTLIER MEASUREMENT AND THE INFLUENCE WEIGHTING (owner request, 2026-08-30)

### 19a. The concentration is EXTREME — one seed can be the whole estimate

`scratchpad/awayhalf/outlier_census.py`, shipped seed set, real libraries. "leave-out" is a DIAGNOSTIC, never a proposal.

| library | seeds | od | top-1 seed's share of the numerator | seeds for 90 % | leave top-100 out |
|---|---|---|---|---|---|
| MO_3021 | 51,199 | 0.6750 | **77.75 %** | **7** (0.014 %) | 0.0338 |
| LBX0190 | 3,367 | 0.9736 | **56.92 %** | 10 (0.297 %) | 0.3369 |
| vcap | 415,773 | 0.0201 | 7.26 % | 73 (0.018 %) | 0.0011 |
| LBX0588 | 86,383 | 0.0165 | 10.33 % | 150 (0.174 %) | 0.0042 |
⚠ On three libraries the top 1 % of seeds carries MORE than 100 % of the numerator (112 %, 151 %, 184 %) — the
rest sums NEGATIVE. The estimate is a small difference between a handful of huge positive contributions and a
mass of small negative ones. ⭐ **od̂ RISES MONOTONICALLY WITH DEPTH on every library** (MO_3021 0.029 → 0.214
→ 0.457 → 0.868; vcap −0.024 → −0.003 → 0.010 → 0.173). A genuine ICC is depth-INVARIANT, so the single-ρ
Beta-Binomial is MISSPECIFIED on real data — that is the finding behind the ceiling, not just an outlier story.

### 19b. The weighting — derived, no constant, no selection

`Var(od̂_s | ρ) ≈ 2/(n(n−1)) + V∞(ρ)`, `V∞(ρ) = 2ρ²(1−ρ)/(1+2ρ)` (MC-verified to 4 s.f. at ρ = 0.01/0.05/0.2/0.5).
Minimum-variance pooling ⇒ `w_s = 1/(½ + c_s·V∞(ρ))`, fixed point from ρ = 0. Weights depend ONLY on `n_s` and ρ,
never on a seed's own data, so nothing is selected, trimmed or discarded and the ratio stays consistent.
⭐ **At ρ = 0 it IS the current estimator** (w ≡ 2, cancels).

**Unbiased at known ρ** (`synth_validate.py`, 40 reps): ρ = 0/0.01/0.05/0.2 recovered as 0.0001/0.0101/0.0502/0.1974
on U[2,400] depths. ⭐ On "4,000 shallow + 8 deep n=2000" at ρ = 0.2 the UNWEIGHTED moment is biased LOW (0.1467)
and the weighted one is right (0.1967) — the weighting is more ACCURATE, not merely more robust.

**Real data:** vcap 0.0201 → **0.0085**; LBX0588 0.0165 → **0.0147**; MO_3021 0.2000 (clamped, raw 0.675) →
**0.1327** — off the ceiling, so it MEASURES instead of clamping; LBX0190 stays clamped (raw 0.974). Top-8 share
of influence: MO_3021 71.8 % → **1.7 %**, LBX0190 87.1 % → **15.8 %**, vcap 3.4 % → 0.3 %.
⚠ My §18 prediction of ≈0.05 for MO_3021 was WRONG (0.1327): the plateau caps a deep seed at ~1/V∞ ≈ 41
shallow-equivalents, not at 1, so the mid-depth bins — which carry genuinely higher od̂ — still dominate.

**Test panel (30 conditions):** every stratum within 3 fragments of the shipped fit (≤ 0.03 %), two strata
marginally better. No harm, as ρ ≈ 0 requires.

### 19c. Antisense: consistent with, NOT proven

`asymmetry.py`. Fully-ANTISENSE seeds OUTNUMBER fully-SENSE at every depth (MO_3021 ratio 0.31–0.67, vcap
0.11–0.62) — that is ordinary own-gene RNA, exactly what the away half is built to ignore, and confirmation the
premise matches the data. The OUTLIERS run the other way: MO_3021's single dominant seed is n = 426 at sense
fraction **0.9930**, and 8 of its top 10 are at 1.0000 — deep, extremely one-sided, and on the SENSE side, where
own-gene RNA cannot put them. ⛔ **This is consistent with unannotated antisense transcription but does not
prove it**: the influence ranking is taken on the away half, which is sense-heavy BY CONSTRUCTION, so the
direction of the top seeds is not independent evidence. A depth-correlated technical artefact (duplicate
families, repeat pile-ups) predicts the same monotone depth gradient. The instrument cannot separate them.

**Ladder (16 conditions):** per stratum, shipped → weighted: ss0.50×off 371,017 → 371,017; deferred
18,612,500 → 18,612,500; ss0.99×off 296,396 → 296,398; ss0.99×on 621,247 → 621,247 (relay likewise, max drift
9 fragments). ⭐ NO HARM on either panel. ⚠ And neither panel can show a BENEFIT: their top-8 influence share is
2.8–57 % against a real library's 71.8–87.1 %, because at truth ρ = 0 with 12–27 k seeds no seed can dominate.
The concentration pathology is a REAL-DATA phenomenon the simulator does not reproduce — the same asymmetry
the blank chromosome was built for, and the reason the census is the evidence rather than the panels.

### 19d. SOLVE IT BY BISECTION, NOT BY ITERATION — and no iteration constant is needed

`converge.py`: plain iteration reaches 6 decimals on every library but takes 13–60 steps to hit a 1e-12
tolerance, which would force an arbitrary `max_iters`. Unnecessary — the root is **bracketed by construction**:
`g(ρ) = F(ρ) − ρ` has `g(0) = F(0) ≥ 0` (the unweighted clipped estimate) and `g(ceiling) = F(ceiling) − ceiling ≤ 0`
(F is clipped to the ceiling), so a root always exists in `[0, ceiling]`. Measured: exactly **1 sign change** on an
81-point grid for all three libraries, and bisection agrees with plain iteration to 6 decimals
(vcap 0.008530, MO_3021 0.132694, LBX0588 0.014665). ⭐ **Bisect to machine precision in a fixed number of
halvings — deterministic, always terminates, introduces no constant.**

### 19e. What the weighting is WORTH on the deliverable — small, and that is the honest framing

`weighted_price.txt`: MO_3021 gDNA call 74,286 (10.0 %) at the clamped 0.2000 → 75,444 (10.2 %) at the weighted
0.1327; vcap 7,065,388 (39.2 %) at 0.0201 → 7,141,903 (39.6 %) at 0.0085. **≈1–2 % relative.** ⭐ So the case for
the weighting is NOT a big correction to the answer — it is that the answer stops depending on ONE seed
(MO_3021: 77.75 % of the numerator) and stops being a clamp. Integrity of the estimate, at ~1–2 % on the call
and 0 % on both panels.

## 20. LANDED — the clamp is visible and the concentration is reported (owner ruling, 2026-08-30)

⭐ `Beta(2,2)` = 0.2 STAYS as the ceiling (owner). `GdnaStrandModel` gains three QC facts, none of which changes
a number: `raw_overdispersion` (the moment before the clip), `clamped_at_ceiling`, and `effective_seeds` —
`gdna_strand.seed_participation`, the participation ratio `(Σ|x_s|)²/Σx_s²` of the seeds' contributions to the
numerator. ⭐ **Threshold-free**: it is the seed count when all seeds contribute alike and → 1 when one seed IS
the estimate, so reporting concentration needs no `k` and no cutoff. `calibrate`'s QC line prints the effective
count and, on a clamp, `CLAMPED at the ceiling from a raw <x> - NOT a measurement`.

| library | od shipped | raw | seeds | **effective** | clamped |
|---|---|---|---|---|---|
| vcap | 0.0201 | 0.0201 | 415,773 | 393.7 | no |
| LBX0588 | 0.0165 | 0.0165 | 86,383 | 935.4 | no |
| LBX0190 | 0.2000 | **0.9736** | 3,367 | **2.9** | ⛔ yes |
| MO_3021 | 0.2000 | **0.6750** | 51,199 | **1.9** | ⛔ yes |
⭐ The two clamped libraries are exactly the two resting on fewer than three effective seeds — concentration and
clamping travel together, and both are now visible instead of inferred.

Gates 16 → 19; ⛔ **five perturbations, and THREE of them initially fired NOTHING** (substituting the seed count
for the ratio; dropping the `|·|`; dropping the away-half weight) — the first two gates tested the function and
the model separately and pinned neither the wiring nor the magnitude. Holes closed, all five now fire.
Suite **3,658**, test panel BYTE-IDENTICAL, preflight 17/17.

## 21. LANDED — the influence weighting (owner approval, 2026-08-30)

`gdna_strand`: `_away_half_parts` (the split), `between_seed_variance` (V∞, derived), `influence_weights`
(`w_s = 1/(½ + c_s·V∞(ρ))`), `away_half_moment(..., overdispersion=0.0)` — the default is the pair-count moment
— and `fit_gdna_strand_overdispersion` solving `g(ρ) = clip(moment(ρ)) − ρ` by **bisection on a bracket that
exists by construction**, halving until the interval is one ULP wide, so **no iteration limit is asserted**.
`seed_participation` now reports concentration AT the fit's ρ, so it describes the estimate actually made.

**Gates 19 → 28.** New: unbiasedness at ρ = 0/0.01/0.05/0.2; `V∞` against the Beta's own fourth moment by
Monte Carlo; the ρ = 0 identity to the pair-count estimator; **the mechanism's own gate** — 4,000 shallow seeds
and 8 deep ones at true ρ = 0.2, where the pair-count moment reads ~0.147 and the weighted fit 0.197, so the
weighting is more ACCURATE and not merely more robust; de-concentration; and the bracket/termination.
⛔ **Six perturbations, and one fired NOTHING**: dropping `(1−ρ)/(1+2ρ)` from `V∞` left every behavioural gate
green — because a wrong weight costs efficiency and never bias, exactly as the derivation claims. The
derivation is therefore pinned DIRECTLY against the Beta's moments; both `V∞` mutations now fire.

**Real data:** vcap 0.0201 → **0.0085** (effective seeds 394 → 3,668); LBX0588 0.0165 → **0.0147** (935 →
2,041); MO_3021 0.2000-clamped → **0.1327, no longer clamped** (1.9 → **1,241**); LBX0190 still clamped from a
raw 0.9736 (2.9 → 86). **Test panel:** 1 fragment on 2 of 30 rows (bisection vs the old exact clip path).
Suite **3,667**, preflight 17/17.

## 22. CLOSE-OUT — what moved to a permanent home, and what stays open

MOVED (and deleted from here in the same edit is not possible for a working record, so this doc is now the
JOURNEY and the permanent docs are the STATE): derivations → `EQUATIONS.md` §6a (away half) and §6b (weighting);
rulings → `DESIGN.md` §3.3a; lessons → `TRAPS.md` as `purity-is-a-property-of-the-annotation`,
`pair-count-weighting-lets-one-seed-decide`, `a-gate-that-restates-the-implementation`; current state →
`ROADMAP.md` §0; the two open questions → `ROADMAP.md` §2.

⛔ **STILL OPEN, and deliberately not done in this session** — the `Beta(14,14)` shrinkage in the *RNA* fit.
Measured 2026-08-30: inert at depth (identical to 6 d.p. on 20 k deep sj) but it BINDS on sparse data (0.0259
shrunk vs 0.0000 unshrunk on 500 sj at depth 3). Deleting it makes a sparse library assert PERFECTLY BINOMIAL
RNA, the most confident strand likelihood assertable. The constant-free candidate — shrink toward the CEILING,
reusing the one constant that exists, so both components coincide when neither is measured — is a behaviour
change needing its own A/B. `ROADMAP.md` §2 question 2.

**LADDER, LANDED CODE (`ab_weighted.py`, 16 conditions):** `od ship == od wtd` on every row — the shipped fit
IS the weighted estimator — and every stratum reproduces the PROTOTYPE's `wtd` column exactly: ss0.50×off
371,017; deferred 18,612,500; ss0.99×off silent 296,398 / relay 413,678; ss0.99×on silent 621,247 / relay
904,882. Drift from the pre-weighting ladder is +2 fragments of 296,398 and +9 of 413,678 (0.002 %). ⭐ Prototype
→ `src/` carried across with no surprise, which is what the DERIVE → PROTOTYPE → A/B → `src/` order is for.

## 23. THE DEPTH GRADIENT IS DIAGNOSED — under-annotated OPPOSITE-STRAND transcription (2026-08-30)

⭐ **PCR duplication is EXCLUDED, twice over**: the scanner skips `BAM_FDUP` by default
(`bam_scanner.cpp`, `skip_duplicates_ = true`), and the cfRNA BAMs are `star.srt.rmdup.collate.bam` — already
de-duplicated. So the gradient is not a technical clumping artefact.

`scratchpad/awayhalf/locate_outliers.py` puts MO_3021's dominant seeds on the genome. All are INTRONS whose
fragments sit ~100 % on the strand that ONLY AN OPPOSITE-STRAND TRANSCRIPT CAN PRODUCE (κ = 0.0023, so
own-gene RNA lands on the other column — the away half selects exactly the seeds own-gene RNA cannot reach):

| share | locus | annotated context | sense frac |
|---|---|---|---|
| **77.75 %** | chr16:3,841,010-3,850,296 | CREBBP intron (−). ⭐ `DPPA3P6`, a **+ strand** processed pseudogene, ends at **3,841,010 — exactly the seed's start** | 0.9930 |
| 4.73 % | chr10:97,678,025-97,679,618 | AVPI1 intron (−). ⭐ `PI4K2A` (**+**) ends 1.6 kb upstream at 97,676,434 | 0.9906 |
| 3.31 % | chr19:46,388,731-46,390,050 | PPP5C intron (+). ⭐ `ENSG00000291298`, a **− strand lncRNA**, starts 464 bp past the seed's end | 1.0000 |

⭐⭐ **Three of three carry an ADJACENT opposite-strand annotated feature whose annotated extent stops at or
near the seed** — the signature of read-through / an under-annotated 3′ extent, not of a novel phenomenon.
⛔ **So the gradient is an ANNOTATION defect, and it is the away half's one structural blind spot.** It also
means a user with a better GTF gets a better fit automatically, which is a property of the design and not a
defect. ⛔ Do NOT filter on it: an "adjacent opposite-strand feature" rule re-imports the purity assumption the
away half exists to remove, and an "impossible under the model" test needs a multiplicity threshold — a new
constant. The influence weighting already caps the effect (MO_3021 top-8 influence 71.8 % → 1.7 %), and what
remains should be REPORTED, not corrected.

## 24. LANDED — the two components reconcile against each other; `Beta(14,14)` is GONE (owner, 2026-08-30)

⛔ **Both conjured constants deleted**: `_PRIOR_ALPHA_BETA` / `_PRIOR_OVERDISPERSION` (the 0.0345 target) and
`_PRIOR_INFORMATION` / `_prior_information()` (the ≈909 weight), plus the shared `_fit_overdispersion` core
that existed to apply them. `test_sweep.py`'s two tests on those constants collapse into one that asserts the
CEILING is the only asserted constant left and that none of the four names can come back.

⭐ **The owner's question — "does shrinking toward the ceiling need a shrinkage constant?" — answered NO, it
needs a different one**: the weight was derived as the max-entropy distribution on `[0, ceiling]` with mean
`od₀`, and at `od₀ = ceiling` that degenerates to a point mass (variance 0, `W = ∞`, the fit would always
return the ceiling). So the ceiling cannot be the target. The reference is the OTHER COMPONENT'S MEASURED
VALUE instead, weighted by its own precision.

**Three pieces:**
1. `between_seed_variance(ρ, mean)` generalised — `V∞(ρ,μ) = 3ρ²[2ρ + μ(1−μ)(1−7ρ)]/[μ(1−μ)(1+ρ)(1+2ρ)] − ρ²`,
   which reduces algebraically to the μ=½ form and is MC-verified on a 3×4 grid. Plus `binomial_scale`,
   `b_s = (2n·pq + 1 − 6pq)/(n−1)`, exactly ½ at μ=½ for every depth.
2. ⛔ **THE RNA FIT IS NOW INFLUENCE-WEIGHTED TOO, and this was REQUIRED, not a bonus.** The first cut compared
   the two components on their NULL informations and it was badly wrong: at κ = 0.0023 a seed's `V∞` at ρ = 0.05
   is **0.285** against gDNA's 0.0043, so pair counts over-credit whichever component sits nearer ½ by orders
   of magnitude. Both fits now report `Σ w_s c_s` — the precision of the estimate actually made — which is the
   only comparable currency.
3. `reconcile_overdispersions`: the weaker borrows its DEFICIT, `(I_s − I_w)/I_s`. ⛔ Found by perturbing my
   own gate: the naive `(I_w·od_w + I_s·od_s)/(I_w+I_s)` blend drags two EQUALLY well-measured components'
   weaker half to the midpoint while the stronger does not move — asymmetric, and a claim neither measurement
   supports.

**⭐ TRUTH-SCORED RESULT — it IMPROVES the panel.** Test panel: silent 460,771 → 460,574 (−0.043 %), relay
168,648 → 168,257 (−0.232 %), and every material row moves TOWARD truth. On the capture-ON rows it lands on
the oracle EXACTLY: `g05 ss0.70 ON` silent 1,281 → **1,169 = true**, `g05 ss0.99 ON` 610 → **573 = true**.
That is the mechanism doing what it is for — where hybrid capture depletes the genomic seeds and the gDNA fit
is thin, it borrows from the component that is not blind.

**Real data (own → shipped):** vcap od_g 0.0085 → 0.0159; LBX0588 od_r 0.0705 → 0.0550 (gDNA better informed
there, so the borrow runs the other way); MO_3021 od_g 0.1327 → **0.0084**; LBX0190 od_g 0.2000-clamped →
**0.0073**. ⚠ **HONEST READING: on 3 of 4 libraries the gDNA dispersion is now dominated by the RNA
measurement.** That is defensible — the spliced seeds are CERTIFIED pure and the genomic ones are contaminated
by the unannotated antisense transcription of §23 — and it dissolves the ceiling clamping. ⚠ But it also means
the away-half's own value is largely superseded on real libraries, and a component that measures HIGH
dispersion is mechanically judged less informative (V∞ grows as ρ²), so the reconciliation carries a
systematic pull toward the LOWER of the two values. Statistically that is just inverse-variance weighting with
a mean-dependent variance; it is recorded here because it is the thing to re-examine first if a real-data
verdict ever looks wrong.

Suite **3,671**, preflight 17/17. Goldens updated twice, magnitudes recorded first (≤1.9 % on one tiny
`loci.gdna`, ≤0.8 % on one `transcript.count`, rest ≤4e-4 relative). Gates 28 → 33; six perturbations on the
reconciliation, all firing.

### 24a. ⛔ THE LADDER'S RELAY MOVEMENT IS A KNIFE-EDGE IN THE RELAY PATH, NOT THE RECONCILIATION

The ladder showed silent improving on 3 of 4 strata (unchanged on the 4th) but relay's `ss 0.50 × off`
worsening 476,246 → 481,195 (+1.04 %), all of it one row: `g98 ss0.50 capture-OFF`, 212,581 → 217,531.
⛔ **It is not the new mechanism.** On that row `od_g` and `od_r` are BOTH 0 before and after, and the change
appears in ALL THREE arms including `true`, which injects `od_g = 0` and bypasses reconciliation entirely.

Bisecting the one parameter left — `od_r`, injected on an otherwise fixed condition:

| od_r | 0 | 1e−12 | 1e−9 | 1e−7 | 1e−6 | 2e−6 | 5e−6 | **1e−5** |
|---|---|---|---|---|---|---|---|---|
| relay error | 217,531 | 217,531 | 217,531 | 217,532 | 217,543 | 217,558 | 217,589 | **212,581** |

⭐ **A 5e−6 change in a nuisance parameter moves relay by 2.3 %** — a threshold in the relay/anchor path, not a
response. The OLD run sat on the far side of it because the deleted `Beta(14,14)` shrinkage lifted that row's
NEGATIVE raw moment (−2.4e−5, truth 0) to a tiny positive **1.1e−5**; the honest fit clips it to exactly 0.
Reconstructing the old estimator exactly reproduces 1.1e−5, and injecting 1e−5 reproduces 212,581 to the
fragment. ⛔ **So the old number was better by ACCIDENT, and no comparison across this edge is attributable.**
Same family as the recorded κ̂ − ½ residue slope at unstranded g00. Not chased here: it is in the relay's
anchor, not the estimator, and it is its own mechanism.

**LADDER VERDICT, silent (no such edge on these rows):** ss0.50×off unchanged at truth; ss0.50×on (deferred)
18,612,500 → 18,606,642; ss0.99×off 296,398 → 296,390; ss0.99×on 621,247 → **616,817** (−4,430). In-scope
strata −0.34 %. ⚠ Relay's ladder movement is NOT attributable until the discontinuity is fixed.
