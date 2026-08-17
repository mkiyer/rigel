# Handoff — 0.8.0 is scoped, and the next session is CALIBRATION

## ⭐ WHERE THE SETTLED PARTS WENT — this file is a HANDOFF, not the state

⛔ It was written to hand one session to the next and is now superseded: the ranked list is
`ROADMAP.md` §1, the state is `ROADMAP.md` §0, and the reference work it points at has moved to
`EQUATIONS.md` §9c / `DESIGN.md` §6b / `ROADMAP.md` §4.2. ⚠ **Nothing may cite it**, and any number
here that disagrees with `ROADMAP.md` is stale by definition.


    ⚠ **A DEV DOC.** Provisional; nothing may cite it. When something settles, MOVE it to its
    permanent home and delete it here in the same edit.

    ⭐ Written 2026-08-14. The permanent docs were rewritten to this scope in the same session —
    `ROADMAP.md` is the ranked list, this file is only the *next steps* and the traps that would
    cost you a morning.

---

## 0. ⭐⭐⭐ THE SCOPE, IN FOUR LINES (owner, 2026-08-14)

1. **Optimise three strata**: unstranded × capture-OFF, stranded × capture-OFF, stranded × capture-ON.
2. ⛔ **DEFER unstranded × capture-ON.** It stays in every benchmark and every table — it is not a
   development target until the other three are done. If it improves as a side effect, free win.
3. ⛔⛔ **The length channel is RETIRED until after 0.8.0.** Do not propose it. Do not rank it.
4. ⭐ **The focus is CALIBRATION**, and the metric is the **calibration result against ORACLE
   calibration** — not the end-to-end transcript number.

---

## 1. ⛔ WHY EQUAL FRAGMENT LENGTHS — the reason is not the one that used to be written down

The old rationale was "neutralise length so residual error is attributable to density and strand."
The real one is stronger and it is why the length channel is retired rather than merely unbuilt:

⭐⭐ **The EM ALREADY USES THE FL DISTRIBUTION.** Give gDNA and RNA a large length difference and the
EM assigns fragments on length alone — **bypassing calibration entirely and masking every bug in
it.** Equal lengths FORCE the calibration phase to be exercised. Calibration has strand and density
(plus message propagation across objects, currently off), and those are the things under test.

⚠ So an equal-length panel is not a limitation of the ladder. It is the instrument working.

---

## 2. WHERE IT STANDS — measured 2026-08-13/14, noise floor 0.996–1.013

⭐ Panel: the rebuilt 16-condition ladder, committed at `59a341e0`,
`{g00, g05, g50, g98} × ss {0.50, 0.99} × capture {off, on}`, 10 M fragments each.
Gates: `simulator_gates.py` **6/6**, `suite_resolves.py` **11/12** (only `(c)`, deferred 2026-07-30).

| | |
|---|---|
| error concentration | unstranded × capture-ON carries **64.5 %** of transcript error and **90 %** of gene error — ⛔ and that is the DEFERRED stratum |
| the three in-scope strata | carry the remaining ~36 % / ~10 %, and are where 0.8.0 is won |
| perfect **per-locus** prior | takes unstranded × capture-ON to **0.099** gene-level; **neutral to worse** (0.79–1.21) on the three in-scope strata |
| perfect **per-transcript** prior | takes the three in-scope strata to **0.37–0.58** gene-level and cuts false-positive mass **100–200×**; does NOT rescue the deferred stratum (0.641) |
| ⛔ the per-transcript lane | is built end to end and **never passed in production** — `pipeline.py:841` omits `rna_prior_weight`. Without it the EM's default rule carries ZERO per-transcript information |

⭐ **Read that table twice before choosing work.** The per-locus prior is the lever on the deferred
stratum; the per-transcript prior is the lever on the three that matter for 0.8.0.

---

## 3. ⭐⭐⭐ NEXT STEPS, in order

### 3a. ~~Make the ruler measurable~~ — ✅ **DONE 2026-08-14. MOVED OUT OF THIS FILE.**

The build is `quant_accuracy.py --arm oracle_ruler` / `oracle_ruler_noop` and
`scripts/design/calibration_vs_oracle.py`; what it measured is `ROADMAP.md` §0 and its rank 1, and the
session write-up is `session_2026-08-14_the_ruler.md` beside this file.

⛔ **Nothing about it is repeated here on purpose** — a dev doc that keeps its own copy of a settled
finding is the second home that then diverges, which is the one rule this sandbox has.

### 3b. Run the scenarios, look at the whole error spectrum, pick ONE target

⭐ Likely the **zero-gDNA** scenarios first, and there is a measured reason:

| condition | ρ_ref shipped → truth | shrink factor shipped → **truth** | transcripts changed |
|---|---|---|---|
| g00 unstranded capture-OFF | 1.600 → **None** | 0.345 → **1.000** | 13,673 / 15,669 |
| g00 stranded capture-ON | 0.957 → **None** | 0.350 → **1.000** | 11,467 / 15,669 |
| g50 unstranded capture-ON | 0.00209 → 2.054 (**981×**) | 0.834 → **0.401** | 13,662 / 15,669 |

⛔ **At zero gDNA the shrinkage has no signal, so it fabricates one from false-positive gDNA and
contracts every transcript to a mean factor of 0.345 when the correct answer is exactly 1.000** —
including on capture-OFF, where the module's own contract says the factor must be 1 by construction.
⭐ The true value is DERIVABLE with no fitting (no gDNA ⇒ no reference density ⇒ factor ≡ 1), which
is the property `ROADMAP.md` §4.1's eleven refused mechanisms all lacked.

⚠ **Real libraries are never at exactly zero gDNA** (owner), so `g00` is a CONTROL, not a target
condition. But it is the only place with a provable answer, and the code fails it by 2.9×.

### 3c. Fix composition where a channel exists — the shrinkage then follows for free

⭐⭐ **The effective-length defect is a SYMPTOM, and this is measured rather than argued.**
`priors.py:341` does `from .capture_eff_length import _global_reference_density` — one split, two
consumers, **one function**. Substituting ONLY the composition arrays and re-running the **shipped,
unmodified** shrinkage function gives the correct factor in both directions (table above).
⛔ So do not "fix the effective length". It is not broken. Its input is.

⚠ **Before choosing a mechanism, establish that a channel EXISTS at the target.** On unstranded ×
capture-ON both channels are dead — strand carries exactly zero at κ=½, and a perfect gDNA landscape
was measured to move exon recovery only 0.0021 → 0.0153 against a truth of ~1.0. That is one reason
that stratum is deferred rather than merely hard.

---

## 4. ⭐⭐ THE SETUP THE NEXT SESSION NEEDS

### Caching, so calibration re-runs in seconds

| cache | state | cost |
|---|---|---|
| `ladder/scan_cache/` | ✅ **16/16** | re-calibrating one condition from cache is **~6 s**, no BAM re-scan |
| `ladder/oracle_cache/` | 12/16 complete; the four `g00` have `gdna`/`mrna`/`nrna` but no `_main` | the origin-split truth |

⭐ **The `g00` `_main` gap is not a blocker and two designs wrongly called it one.** `_main` is the
UNDRAINED FULL PAYLOAD, i.e. the same scan as the plain scan cache — so read it from
`ladder/scan_cache/<cond>` when `oracle_cache/<cond>/_main` is absent. A two-line fallback; verified
working.

⚠ **They are NOT byte-identical, and the difference is the known one rather than a defect.** Measured
on `g05 ss0.99 capture_off`: `manifest.json` and `strand.npz` identical, and of `payload.npz`'s 31
arrays **25 are exactly equal and 6 differ by 3.6e-14 … 6.4e-12** — precisely the six float64 banks
(`boundary_spliced_mass`, `boundary_unspliced_mass`, `boundary_unspliced_inv_length_sum`,
`region_contained_inv_opportunity_sum`, `sj_mass`, `sj_inv_length_sum`). They are two independent
scans of one BAM and float addition is not associative across worker threads — the same fact that
makes `rescan_panels.py`'s gate a derived budget instead of exact equality. ⭐ Every INTEGER bank is
identical, so the substitution is sound for measurement. ⛔ It would not survive a byte-identity gate,
so do not use it as one.

⚠ `panel.py status` will still print `oracle cache 12/16` and a ✘, because it counts directories.

⛔ **`pass0_vs_oracle.py` HOLDS OUT every zero-gDNA condition** ("N zero-gDNA row(s) held out as
false-positive checks"), which is why `g00` has no `_main`. That is deliberate on its part; it is not
a reason `g00` cannot be measured.

### The metric

✅ **BUILT 2026-08-14 as `scripts/design/calibration_vs_oracle.py`** — measured 5–12 s per condition,
the whole ladder in ~2 min, no solver. ⛔ Its numbers live in `ROADMAP.md` §0 and are not copied here.

⛔ **Score per stratum, never pooled.** On 8 of 12 conditions the shipped composition is already
within ~1 % of truth, so a pooled number reads any calibration experiment as a small effect.

---

## 5. ⛔ TRAPS THIS SESSION FOUND — each would cost you a morning

**An arm name can silently disable the thing it names.** `quant_accuracy.run_condition` keys on three
separate prefix tests — `seeded()` at `:420`, the oracle test at `:590`, and
`install_computed_weights` at `:596`, the last two both `arm.startswith(("alloc_", "allocg_"))`. An
arm called `oracle_cal_alloc_harmonic` passes the oracle test and **fails the other two**: it gets the
shipped warm start and **never installs a weight vector**, then scores as a prior-substitution arm and
reads as "the weighting function with perfect input". `TRAPS: an-ablation-that-never-ran`, caused by
the name alone.

**`quant_accuracy.py --report` refuses arms with different row sets**, which is correct and will catch
you the first time you compare a 12-condition oracle arm against a 16-condition base arm. Subset the
base arm; rows are scored independently per condition, so a subset is exactly equivalent.

**The zero-gDNA control cannot test the shrinkage even in principle** — shrinkage needs a gDNA signal
and there is none (owner). It tests the FALSE-POSITIVE behaviour, which is a different and also
valuable thing.

**Aggregate accuracy is the wrong statistic for a soft-min.** The per-transcript weighting function
is a lower-tail statistic: at `g05` the composition arrays are 99.6 % accurate in aggregate and the
algorithm's output still moves **28 %** under perfect input, because the support changes (9,250
regions and 16,234 boundaries have truth RNA exactly 0 while calibration hands them mass). ⛔ At
`g00` the support is **bit-identical** under perfect input — so the condition the weighting function
was refused on is the one condition where better input provably cannot help.

---

## 6. What is NOT next

⛔ **Not the length channel.** Retired until after 0.8.0 by owner ruling; `EQUATIONS.md`'s derivations
are kept as a record and carry a deferral banner.

⛔ **Not unstranded × capture-ON.** Keep measuring it; do not optimise it.

⛔ **Not a new per-transcript weighting mechanism yet.** The lane is not wired in production
(`pipeline.py:841`), the ceiling that would price it was measured with a wrong ruler, and
`ROADMAP.md` §4.1's graveyard is eleven mechanisms that died on the zero control. Wire, then price,
then propose.
