# Ablation re-test — does the locus-EM per-fragment `log 0.5` bias abundances *under overdispersion*?

**Status: COMPLETE (2026-06-07) → VERDICT UPHELD: do NOT build the EM strand term; REVERT the
ablation.** The Phase-0 conclusion holds even with overdispersion now present in both the
simulator and calibration. Re-opened it (the runbook below) because a load-bearing premise was
gone; the re-test (§9 RESULTS) re-confirms it. Design / runbook retained below.

## 9. RESULTS (2026-06-07) — Q1 negative; verdict upheld

Sweep suite: gDNA {none, high} × strand-overdispersion {od00, od05, od14} × ss {0.99, 0.50},
nRNA none, 25 genes / 500 kb / 80k RNA fragments (`simulate_suite … --gdna-strand-overdispersions
0.0,0.05,0.14`); baseline quant via `bench_calibration --run`; gating via `em_strand_gating`.

- **The hard-label leak exists and GROWS with overdispersion** (theory's mechanism confirmed):
  gDNA→RNA leak at ss99 rises **2.9% → 3.0% → 3.4%** across od00 → od05 → od14 (higher, ~4.5%, at
  ss50). The per-fragment ρ_d = 0 penalty bites harder as the data overdisperses.
- **But soft abundances are NOT biased (Q1 negative).** Total mRNA `count_em` is **never
  inflated** — ≈ **−1.7% (under, FP-safe), flat across od00/05/14** at gdna_high; the zero-gDNA
  control is *more* under (−2 to −15%, the gDNA-independent baseline). Per-locus, the largest
  Δmrna at gDNA-heavy loci are **negative** (−686, −454 at od14/ss99); the few positives are tiny
  and do not track gDNA load (Δmrna/gdna ≈ 0.04). Leak and siphon cancel in the soft counts.

**Conclusion.** With overdispersion now in the data and fit by calibration, the per-fragment
`log 0.5` leak remains a **hard-label artifact** (it grows with od) that does **not** bias the soft
per-transcript abundance Rigel reports. Pillar 1 of the Phase-0 verdict survives the removal of
pillar 2. Q2 (fix efficacy) is moot — there is no soft-abundance bias to fix. **Action: revert the
held-back `strand_neutral_split` ablation (§7); close the em_strand line.**

---

**Status:** design / runbook (2026-06-07). Re-opens the Phase-0 verdict of
[02_mstep_implementation_plan.md](02_mstep_implementation_plan.md) because a load-bearing premise
of that verdict no longer holds.

## TL;DR

The locus EM scores every gDNA fragment's strand with a constant `log 0.5` — the **ρ_d = 0**
(rigid 50/50) limit of the correct Beta-Binomial(½, ρ_d) ([01](01_theory_and_fix.md) §2). Phase 0
concluded this is harmless for *quantification* (the co-stranded gDNA→RNA leak is a hard-label
artifact; **soft** abundances are never inflated) and that the clean latent-κ fix is "provably
inert because ρ_d = 0." **That conclusion was reached on data where ρ_d really was 0** — the old
simulator generated gDNA strand rigidly 50/50, and the calibration decode used the Binomial limit.

Both are now false:
- The simulator **generates** gDNA *and* RNA strand overdispersion (`GDNASimConfig.strand_overdispersion`, the `--gdna-strand-overdispersions` sweep axis).
- Calibration **fits and applies** `gdna_strand_overdispersion` / `rna_strand_overdispersion` (shipped: docs/em_strand/03–05).

So there is now a genuine **mismatch**: the data has ρ_d > 0, calibration knows it, but the locus
EM per-fragment factor still assumes ρ_d = 0. [01](01_theory_and_fix.md) §2 predicts this is when
the `log 0.5` over-confidence penalty bites *hardest* (it grows as the local split deviates from
50/50, and overdispersion is exactly extra deviation). **The re-test asks whether that mismatch
biases soft abundances.** If yes → build the fix (productionize the held-back ablation). If no →
the verdict stands, now confirmed *with* overdispersion → revert the ablation and close the line.

## 1. Background (what was decided, and on what)

- **The defect** ([01](01_theory_and_fix.md) §2): per-fragment `log 0.5` is the overconfident
  ρ_d = 0 gDNA strand model; at a stranded library it penalizes co-stranded gDNA → that gDNA
  leaks into RNA. Measured signature: co-stranded leak ≫ antisense leak at ss99, symmetric at ss50.
- **The fix sketch** ([02](02_mstep_implementation_plan.md) §2): a per-locus latent gDNA strand
  rate `κ_g ~ Beta(a_d, a_d)` (a_d from calibration's ρ_d), replacing the constant `log 0.5` with
  a digamma orientation factor (closed-form VBEM); reduces to `log 0.5` as ρ_d → 0. "B + C" =
  aggregate count×strand split in the M-step + neutralize the per-fragment split factor
  (`strand_neutral_split`).
- **Phase-0 verdict** ([02](02_mstep_implementation_plan.md) §4): **do NOT build it.** Two pillars:
  1. **Soft abundances unbiased** — `em_strand_gating.py` showed total mRNA (`count_em`) is never
     *inflated* (always slightly under; a gDNA-independent −3% baseline dominates). The leak is a
     hard-label QC artifact (leak and siphon cancel in soft counts).
  2. **ρ_d = 0 by construction** — the decode used the Binomial limit *and the sim had no
     overdispersion*, so the latent-κ_g model is provably inert (≡ `log 0.5`).
- **A theoretical caveat that survives even if we re-open** ([02](02_mstep_implementation_plan.md)
  §3.1): a Beta(½, ρ_d) rate does **not** anchor "sense count = antisense count"; with large ρ_d it
  *learns the observed asymmetry* (the wrong direction). So even a correct latent-κ_g model may not
  reduce the leak. The re-test must check fix *efficacy*, not just leak existence.

## 2. Why re-test now — pillar 2 is gone

Pillar 1 (soft unbiased) was measured at ρ_d = 0. With the simulator and calibration now carrying
overdispersion, the per-fragment ρ_d = 0 penalty operates against data it no longer matches. The
open question is whether the resulting hard-label leak now **also** biases the soft abundances —
i.e. whether pillar 1 still holds when pillar 2 is removed. Nothing in the Phase-0 data answers
this; it must be re-measured on overdispersed simulations.

## 3. The re-test questions

- **Q1 (the decision gate).** On simulations *with* gDNA strand overdispersion, does the locus-EM
  per-fragment `log 0.5` bias **soft** per-transcript abundance (`count_em`) beyond the ±2.2% pool
  tolerance — specifically, does mRNA over-estimation now track the per-locus gDNA load and grow
  with the simulated overdispersion? (Re-run the Phase-0 gating measurement across an
  overdispersion sweep.)
- **Q2 (only if Q1 is positive — fix efficacy).** Does neutralizing the per-fragment split factor
  (`strand_neutral_split`, both scoring paths) — as a *stand-in for* the full latent-κ_g M-step —
  reduce that abundance bias, or does it make the leak worse (the [01](01_theory_and_fix.md) §4.2
  finding: C-alone regressed the antisense catch) / move in the wrong direction
  ([02](02_mstep_implementation_plan.md) §3.1)? This decides whether the *full* B+C build is worth
  it before we spend on the native VBEM update.

## 4. Prerequisites

1. **Build with the ablation present.** The `strand_neutral_split` instrumentation is in the
   working tree (uncommitted): `scoring.cpp` (both scoring paths — region 542 multimapper *and*
   ~983 unique-mapper, per [01](01_theory_and_fix.md) §4.2), `scoring.py`, `config.py`
   (`FragmentScoringConfig.strand_neutral_split`), `pipeline.py` (`RIGEL_STRAND_NEUTRAL_SPLIT`
   env-var toggle). Recompile so the toggle is live:
   `pip install --no-build-isolation -e .` (in the `rigel` conda env).
2. **Overdispersion-capable simulator** — shipped (`--gdna-strand-overdispersions`).
3. **Diagnostic scripts** (working tree, untracked): `scripts/debug/em_strand_gating.py` (soft
   abundance vs truth, per-locus vs gDNA load), `em_strand_leak_stratify.py` (co/anti hard-label
   leak), `em_strand_fate.py` (ZF fate decode). Verify each still parses the current
   `quant`/`loci`/`summary` schema; fix if drifted (they predate the calibration rebuild).

## 5. Exact steps

> All commands inside the activated `rigel` env:
> `source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel`.
> Use a scratch dir, e.g. `SUITE=$(mktemp -d)/od_retest`.

**Step 0 — rebuild** with the ablation toggle (Prereq 1).

**Step 1 — generate an overdispersion sweep suite.** One reference, sweeping gDNA strand
overdispersion at a stranded library (where the penalty bites) and unstranded (the null), at
gDNA-heavy load (where bias, if any, is largest). Example:
```bash
python scripts/sim/simulate_suite.py --outdir "$SUITE" --profile full \
    --gdna-rates 0.0,1.0 --gdna-labels none,high \
    --strand-specificities 0.99,0.50 \
    --gdna-strand-overdispersions 0.0,0.05,0.14 \
    --gdna-strand-overdispersion-labels od00,od05,od14
```
(`od00` is the Phase-0 baseline; `od14` ≈ the calibration prior. The sweep is the key new axis.)

**Step 2 — quant every condition (baseline, ablation OFF).**
```bash
python scripts/sim/bench_calibration.py --sim-base "$SUITE" --run --force
```

**Step 3 — Q1 gating.** Soft `count_em` vs truth, per condition and per locus:
```bash
python scripts/debug/em_strand_gating.py --sim-base "$SUITE"   # adapt to its actual CLI
```
Read out, **as a function of simulated overdispersion**: (a) total mRNA inflation vs the
gdna_none −3% baseline; (b) does per-locus mRNA over-estimation track locus gDNA load *and grow
with od*? Cross-check the hard-label leak trend with `em_strand_leak_stratify.py`.

**Step 4 — Q2 (only if Q1 positive): ablation ON.** Re-quant with the per-fragment split
neutralized, into a parallel dir, and re-run the gating:
```bash
RIGEL_STRAND_NEUTRAL_SPLIT=1 python scripts/sim/bench_calibration.py \
    --sim-base "$SUITE" --run --force --out-suffix _neutral   # or a second SUITE
python scripts/debug/em_strand_gating.py --sim-base "$SUITE_neutral"
```
Compare baseline vs neutralized abundance bias (and the co/anti leak), at each od level. Watch for
the §4.2 regression (antisense catch lost) and the §3.1 wrong-direction effect.

**Step 5 — decide** per §6 and record the numbers in this doc.

## 6. Decision tree

- **Q1 negative** — soft abundances stay unbiased even as simulated overdispersion grows (mRNA
  never inflated beyond the gDNA-independent baseline; no od-dependent per-locus trend). → **The
  Phase-0 verdict holds, now confirmed under overdispersion.** **REVERT the ablation:** delete
  `strand_neutral_split` from `scoring.cpp` / `scoring.py` / `config.py` and the
  `RIGEL_STRAND_NEUTRAL_SPLIT` block in `pipeline.py`; keep the debug scripts (or move under
  `docs/em_strand/` provenance); mark this line CLOSED. Rebuild, full suite green, commit.
- **Q1 positive + Q2 shows neutralization reduces the bias** → the fix is warranted. **Productionize:**
  implement the latent-κ_g closed-form VBEM M-step ([02](02_mstep_implementation_plan.md) §2/§5)
  (aggregate count×strand split using `(N_sense, N_anti)` + calibration's ρ_d, reduces to `log 0.5`
  as ρ_d→0) **and** enable `strand_neutral_split` (both paths) so the per-fragment factor stops
  double-counting the split. Validate against [01](01_theory_and_fix.md) §6.1 targets (co-leak ↓
  without losing the antisense catch; ss50 unchanged; gdna_none FP unchanged; pool within ±2.2%).
- **Q1 positive + Q2 shows neutralization does NOT help (or worsens it)** — consistent with the
  §3.1 wrong-direction concern. → Document that the per-fragment ablation alone is insufficient and
  the latent-κ_g model would likely "learn the asymmetry"; **revert the ablation** and record the
  residual as a known hard-label artifact (not a quantification bug). Only revisit with a model that
  *anchors* sense=antisense (not a plain Beta(½) rate).

## 7. What "revert" vs "productionize" touches (code inventory)

**Held-back ablation surface** (uncommitted working tree): `src/rigel/native/scoring.cpp`
(`strand_neutral_split_` flag + both gDNA-strand-term branches), `src/rigel/scoring.py`
(`strand_neutral_split` param), `src/rigel/config.py` (`FragmentScoringConfig.strand_neutral_split`),
`src/rigel/pipeline.py` (`import os` + `RIGEL_STRAND_NEUTRAL_SPLIT` block + the kwarg pass-through),
plus `scripts/debug/em_strand_*.py` and `scripts/profiling/profiler.py` (interim).

- **Revert** = remove all of the above (`git checkout`/edit out), rebuild, confirm full suite green.
- **Productionize** = keep + commit `strand_neutral_split`, and add the native latent-κ_g M-step in
  `em_solver.cpp` (`apply_grouped_prior_update`) per [02](02_mstep_implementation_plan.md) §5, with
  tests (reduces to `log 0.5` at ρ_d→0; recovers abundance under overdispersion).

## 8. Risk / honesty notes

- The §3.1 caveat is real: a symmetric Beta(½, ρ_d) latent rate may not pull toward sense; the
  re-test's Q2 (cheap neutralization stand-in) is precisely to avoid building the expensive VBEM
  update before knowing it helps.
- The simulator's overdispersion is generated per *exon-derived region* (matching the calibrator's
  seed partition); confirm the locus EM's `(N_sense, N_anti)` actually sees that overdispersion
  (a single small locus may not span enough regions). If the per-locus signal is too weak even at
  high simulated od, that itself is an answer (the effect is sub-threshold for quantification).
