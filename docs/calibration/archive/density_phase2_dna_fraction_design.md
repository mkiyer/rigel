# Density Phase 2 — DNA-fraction estimator at splice seams (design + implementation plan)

**Status:** design / for review. 2026-06-09. Implements Phase 2 of
`density_sweep_implementation_plan.md`. **Not yet implemented** — this plan surfaces a measurement
that reframes Phase 2 vs Phase 3 and asks for a decision before coding (§2, §8).

---

## 0. TL;DR

Phase 1 fixed the count-prior **density** (0.62 → ~11 on GENE0037 exons) but the deconv gДНК
fraction still lands at ~0.45 vs oracle ~0.77. A benchmark measurement (`gdna400/ss0.99/capture-on`,
oracle) shows **why**, and reshapes the plan:

- **The strand clue alone already recovers the true fraction (~0.77 to ~1–2%) at ss=0.99.** The deconv
  fails only because the count-prior **concentration** (`density·eff_len` ~12k) drowns it. ⇒ **Phase 3
  (honest concentration) closes the high-SS leak with no new estimator.**
- **The seam DNA-fraction also recovers ~0.77.** Its *distinct* value is the **low-SS / unstranded**
  regime: when nascent RNA is low it is an **SS-independent** gДНК:RNA ratio, so it works where the
  strand clue is uninformative — and it repairs the degenerate empty exon (eff_len=0) that the strand
  clue mis-calls.

**Recommendation:** treat Phase 3 (concentration) as the primary high-SS fix and Phase 2 (DNA-fraction)
as the low-SS complement; build **both estimators + a comparison harness** and let the data pick the
per-condition selection rule (the plan already called for this). Decision requested in §8.

Measured (GENE0037 exons, ss0.99/capture-on; `oracle gf | strand-only gf | seam DNA-frac`):
`369 0.770|0.770|0.762  371 0.937|0.933|0.904  376 0.732|0.726|0.773  378 0.772|0.769|0.770
382 0.719|0.704|0.776` (r374, an eff_len=0 empty exon: oracle 0, strand mis-calls 1.0, DNA-frac 0.762).

---

## 1. The DNA-fraction estimator

At a **valid splice seam** — a boundary where one side carries the **exon** bit and the other the
**intron** bit of the **same strand** (`ex+|in+`, `in+|ex+`, `ex−|in−`, `in−|ex−`) — spliced reads
genuinely flow across the junction. There, for the exon region's side:

```
gdna_cross  = strand_clean_gdna_frac(side_sense, side_unspliced_total, κ) · side_unspliced_total
spliced     = side_spliced_sense + side_spliced_antisense        # mature multi-exon RNA (SS-independent)
dna_fraction = gdna_cross / (gdna_cross + spliced)               # require spliced > 0 (a real seam)
```

Averaged over the exon's qualifying seam sides (left boundary `(r−1,r)`, right boundary `(r,r+1)`),
this is the exon's count-prior **mean** gДНК fraction, applied to the whole region's contained mass.

**Gating (your earlier guidance).** Compute it **only** at same-strand exon↔intron seams with
`spliced > 0`. NOT intergenic↔exon (TSS/TES — no splicing), NOT `ex+|ex−` (antisense — no intron bit):
at a non-splice boundary we would see 0 spliced + some unspliced and wrongly call 100 % DNA. The gate
falls out of the signature transition (one side exon-S, other intron-S) plus the `spliced > 0` guard;
where it fails, fall back to the Phase-1 density-derived fraction.

**Why it can beat the density (and the empirical confirmation).** It is a *same-position* DNA/RNA
ratio: the local probe enrichment multiplies both the crossing gДНК and the spliced RNA, so it
largely cancels (the worry that the crossing gДНК straddles the depleted intron and is *under*-captured
relative to the fully-exonic spliced read turned out empirically small — see §0, DNA-frac ≈ oracle —
consistent with near-binary capture). And it is **SS-independent** when nascent is low: the spliced
count is RNA regardless of strand specificity, so the ratio holds at ss=0.5 where the strand clue is
blind.

---

## 2. The measurement that reframes Phase 2 vs 3 (READ FIRST)

The deconv's count prior is `Beta(mean, concentration)` with `mean = density·eff_len/mass` and
`concentration = count_evidence = density·eff_len`. On a captured exon the concentration is ~12k, so
the count prior is near-deterministic at its (Phase-1) mean ~0.40 and **overrides** the strand
likelihood — even though the strand clue's own estimate is ~0.77 (≈ oracle). So:

- The high-SS residual is **not** a missing estimator; it is **count over-confidence**. **Phase 3**
  (concentration = the *observed* crossing/seed count, not `density·eff_len`) lets the strand clue
  govern and should recover ~0.77 at ss≥~0.9 on its own.
- **Phase 2** earns its keep at **low SS** (strand uninformative) and on **degenerate/empty** exons
  (no contained signal), where the DNA-fraction's SS-independence is the only handle.

This argues for **measuring Phase 3's effect first** (cheap: it is a concentration change) and scoping
Phase 2 to the regimes where the comparison harness shows it adds value — rather than implementing a
new estimator on the assumption it is needed everywhere. **Open decision in §8.**

---

## 3. Integration + the count-clue→deconv interface cleanup

Today `deconv_regions` re-derives the count-prior mean from `node_density.density`
(`density·eff_len/mass`). To inject a *fraction* (Phase 2) and a separate *concentration* (Phase 3)
cleanly, the count clue should own and expose both explicitly:

- **`NodeDensity` gains `count_gdna_frac` (per-region count-prior MEAN).** `node_gdna_density` sets it:
  - observable region: `clip(density·eff_len/mass, 0,1)` (= the strand-cleaned contained fraction);
  - exon with a qualifying seam: the **seam DNA-fraction** (Phase 2);
  - exon without a seam: `clip(density·eff_len/mass)` (Phase-1 behaviour).
- **`deconv_regions` consumes `count_gdna_frac` directly** (drop its `density·eff_len/mass` line).
- `density` stays (the `deconv_sides` fallback + diagnostics); `count_evidence` stays the concentration
  (Phase 3 revises). **Bit-identical** for observable + no-seam exons (only seam exons change), so the
  refactor lands behind a golden diff isolated to seam exons.

This is the natural seam for both phases and removes the "density re-interpreted as a fraction" smell.
(`node_gdna_density` could then be renamed — it returns more than a density — but that is optional
churn; noted, not bundled.)

---

## 4. Codebase cleanup opportunities (found while designing)

**Bundle with Phase 2 (low-risk, on-path):**
1. **Explicit `count_gdna_frac`** in `NodeDensity` (§3) — clarifies the count-clue→deconv contract.
2. **Remove the unused `NodeDensity.gdna_mass`** field — nothing consumes it (`derive` uses the
   *deconv* gdna_mass); it is dead weight.
3. **Delete obsolete debug scripts** that target the removed sweep / old signature:
   `scripts/debug/gene37_sweep_walk.py` (describes the deleted sweep),
   `gene37_region_boundary_autopsy.py`, `impute_prototype.py`, `locus_calibration_dive.py` all call
   `node_gdna_density(..., gdna_frac=…)` (removed signature) → currently broken. Delete, or port the
   one or two still useful to the new API. They are TEMP debug artifacts.

**Note, do NOT bundle (avoid scope creep):**
4. `boundary_eff_length(fl_pmf)` is just `fl_mean(fl_pmf)` — collapse the two-name duplication later.
5. The seam test belongs in `signature.py` as a small `is_splice_seam(sig_a, sig_b)` helper next to
   `count_observable_masks`'s logic (single home for signature predicates).
6. CLAUDE.md still references the long-removed `density.py`/`estep.py`/sweep wording — refresh after
   this line of work settles.

---

## 5. Open questions (must resolve before / during implementation)

1. **Phase 2 vs Phase 3 ordering (the big one).** Given §2, do we (a) implement Phase 3 first and
   re-measure the residual before deciding if Phase 2 is needed, or (b) implement both + the
   comparison harness now and pick the selection rule from data? **Recommend (a)→(b): a quick Phase-3
   concentration experiment first, then build Phase 2 scoped to where it helps.** Needs your call.
2. **Nascent RNA in the DNA-fraction numerator.** The benchmark is `nrna_none`. With nascent present,
   the seam's unspliced crossing = gДНК **+ nascent**, so the numerator needs `strand_clean` to
   subtract nascent — which is **SS-dependent**, eroding the low-SS SS-independence that is Phase 2's
   whole point. Must measure on an `nrna>0` condition before claiming low-SS robustness.
3. **Denominator population mismatch.** The DNA-fraction's RNA term is *spliced multi-exon* mature RNA,
   while the exon's contained-unspliced RNA is *single-exon mature + nascent*. They are proportional
   for a well-expressed gene (same molecules, different positions), but verify on multi-isoform /
   alt-spliced loci where it may break.
4. **Degenerate / empty seams.** Require `spliced > 0` and `(gdna_cross + spliced) > 0`; otherwise the
   gate fails → density fallback. The eff_len=0 exon (r374) shows the DNA-fraction *helps* here, but
   the rule must be explicit.
5. **Selection rule.** Per-condition: seam DNA-fraction where the seam qualifies AND (low SS OR no
   strand identifiability), else strand-governed (Phase 3); density elsewhere. Decide the exact rule
   from the comparison harness, not a priori.
6. **Multiple seams / disagreement.** If the left and right seam DNA-fractions disagree (alt-splice,
   partial enrichment), average vs min vs evidence-weight? Decide from data.

---

## 6. The comparison harness (the deliverable that answers §5)

Extend the benchmark tooling into a per-exon table over the condition matrix
(`gdna×ss×capture`, ± nrna): for each exon region report `oracle gf`, `strand-only gf`,
`Phase-1 density gf`, `seam DNA-frac`, `deconv gf` (current), and `deconv gf` under a
reduced (Phase-3) concentration. Summaries: per-condition leak (net gДНК→RNA) under each estimator;
where DNA-fraction beats strand and vice-versa. This is the empirical basis for the §5 selection rule
and the Phase-2-vs-3 decision — and reuses the oracle-vs-calibration scaffolding already written for
the Phase-1 validation.

---

## 7. Tests & acceptance

- **Unit:** `is_splice_seam` truth table (ex+|in+ ✓, in+|ex+ ✓, ex−|in− ✓, ex+|ex− ✗, interg|ex+ ✗,
  ex+|in− ✗ cross-strand); DNA-fraction on a hand-built seam (known gДНК/spliced → known fraction);
  gate falls back to density when `spliced=0`.
- **Bit-identical:** the `count_gdna_frac` refactor leaves observable + no-seam-exon deconv unchanged
  on golden (isolate the seam-exon diff).
- **Worked (oracle):** GENE0037 exon deconv gf → ~0.77 (from ~0.45) under the chosen path; introns
  stay ~depleted; `gdna_none` shows no false gДНК rise.
- **Suite + golden** regenerate (seam-exon scenarios shift); full suite green.
- **Benchmark:** net gДНК→RNA leak drops at capture-on/high-gДНК across ss, with the comparison harness
  documenting which estimator carried each condition.

## 8. Decision requested

1. **Ordering:** Phase 3 concentration experiment first (recommended), or both-now + harness?
2. **Cleanup scope:** OK to bundle items 1–3 of §4 (explicit `count_gdna_frac`, drop dead `gdna_mass`,
   delete obsolete debug scripts) with the implementation?
3. **DNA-fraction selection rule:** decide a priori, or strictly from the comparison harness?

## 9. Implementation phases (once §8 is decided)

- **P2.0 — Comparison harness + Phase-3 concentration experiment** (measure-first; no production
  behaviour change beyond a concentration knob to probe).
- **P2.1 — Interface cleanup:** explicit `count_gdna_frac`; drop `gdna_mass`; `is_splice_seam` helper;
  delete obsolete debug scripts. Bit-identical golden.
- **P2.2 — DNA-fraction estimator** at gated seams feeding `count_gdna_frac`, with the selection rule
  from P2.0. Golden regen; benchmark.
- (Phase 3 — honest concentration — sequenced per the §8.1 decision.)
