# The preconditions, run — and what they change

    ⚠ **A DEV DOC.** Provisional; nothing may cite it. When something settles, MOVE it to its
    permanent home and delete it here in the same edit.

    ⭐ Measured 2026-08-13 at `59a341e0`, on the rebuilt 16-condition ladder. Suite green:
    3,402 passed / 10 xfail / 0 skipped. `messages OFF`, `length_likelihood OFF`, EM seed pinned
    at 20260807.

---

## 0. What was run, and what it cost

| precondition | state |
|---|---|
| panel committed, reproducible from the tree | ✅ `59a341e0` |
| `simulator_gates.py` | ✅ **6/6** |
| `suite_resolves.py` | ✅ **11/12** — only `(c)` replicate-pairs, deferred by owner 2026-07-30 |
| `panel.py status` | ✅ build / simulate / scan-cache complete |
| `quant_accuracy` arms on the panel that exists | ✅ `base`, `base_reseed`, `oracle`, `oracle_gdna`, `oracle_rna`, `oracle_efflen` |

⛔ **The oracle arms cover 12 of 16 conditions, and the four missing ones are all `g00`.**
`pass0_vs_oracle.py` — which populates the oracle cache — HOLDS OUT every zero-gDNA condition
("N zero-gDNA row(s) held out as false-positive checks"). This is structural, not a half-built
cache. ⭐ It happens to be exactly the set the strata are computed over (`ALL` excludes `g00`), so
every stratum row below is like-for-like. ⛔ **The one thing it costs is the `g00` row of every
ceiling table, which reads `nan`** — so "a perfect prior scores 0.95 at `g00`" is currently
UNMEASURABLE on this panel, and that number should not be carried forward.

⚠ `quant_accuracy.py --report` REFUSES to aggregate arms with different row sets — the
`arm_identity.py` lesson, enforced in the reporter. The `base` arm is therefore subset to the same
12 conditions for every comparison here; rows are scored independently per condition, so a subset
is exactly equivalent to having run the subset.

⭐ **THE NOISE FLOOR, re-recorded this session** (`base` vs `base_reseed`, seed+1, all 16
conditions): **0.996–1.013**, mostly 0.998–1.000. Every ratio below is far outside it except where
noted.

---

## 1. ⭐⭐⭐ THE THREE POOLS — gDNA vs SYNTHETIC NASCENT vs ANNOTATED

⛔ "nascent" is a SIMULATOR input and an INDEX construct, never a population the solver has
(AXIOM 0). This scores ASSIGNMENT. On an `nrna_none` panel its truth is **0** everywhere, so that
column is pure false positive and a relative error against it is undefined, not large.

⚠ gDNA TOTAL folds the intergenic attribution back in (`gdna_est + n_intergenic`): the arm file
splits the gDNA pool into the part the EM resolved at annotated loci and the part attributed to
intergenic space, and both are gDNA.

### 1a. By stratum — `base`

| stratum | pool | est | true | Δ | Δ% |
|---|---|---:|---:|---:|---:|
| stranded × OFF | gDNA | 15,248,662 | 15,300,004 | −51,342 | −0.3 % |
| | nascent | 868 | 0 | +868 | n/a |
| | annotated | 14,744,019 | 14,699,996 | +44,023 | +0.3 % |
| stranded × ON | gDNA | 14,518,528 | 15,300,004 | −781,476 | −5.1 % |
| | nascent | 42,290 | 0 | +42,290 | n/a |
| | annotated | 15,381,871 | 14,699,996 | +681,875 | +4.6 % |
| unstranded × OFF | gDNA | 15,223,957 | 15,300,004 | −76,047 | −0.5 % |
| | nascent | 585 | 0 | +585 | n/a |
| | annotated | 14,769,116 | 14,699,996 | +69,120 | +0.5 % |
| ⛔ **unstranded × ON** | **gDNA** | **3,270,675** | 15,300,004 | **−12,029,329** | **−78.6 %** |
| | **nascent** | **1,735,497** | 0 | **+1,735,497** | n/a |
| | **annotated** | **24,936,781** | 14,699,996 | **+10,236,785** | **+69.6 %** |
| ⛔ g00 control (4 cond) | gDNA | 1,500,175 | 0 | +1,500,175 | n/a |
| | nascent | 563,643 | 0 | +563,643 | n/a |
| | annotated | 37,936,182 | 40,000,000 | −2,063,818 | −5.2 % |

⭐ **The failure is a mass-conservation story and it is ONE stratum.** 12.0 M gDNA fragments go
missing; 10.2 M reappear as annotated RNA and 1.7 M as synthetic nascent. Three strata are within
−0.3 % … −5.1 %.

⭐ **The sub-split names the mechanism.** Off capture, intergenic carries ~8.0 M (52 % of all
gDNA); under capture it collapses to ~0.33 M (2 %), so essentially all gDNA must be identified by
the EM at annotated loci. With strand the EM does it (14.19 M); without strand it recovers 2.94 M.
The failure needs BOTH — capture removes the intergenic reference and κ=½ removes the strand
channel.

### 1b. By stratum — `oracle` (perfect `LocusPriors` injected)

| stratum | pool | est | true | Δ% |
|---|---|---:|---:|---:|
| stranded × OFF | gDNA | 15,337,084 | 15,300,004 | +0.2 % |
| | nascent | 805 | 0 | n/a |
| | annotated | 14,655,660 | 14,699,996 | −0.3 % |
| stranded × ON | gDNA | 14,793,162 | 15,300,004 | −3.3 % |
| | nascent | 51,448 | 0 | n/a |
| | annotated | 15,098,079 | 14,699,996 | +2.7 % |
| unstranded × OFF | gDNA | 15,344,082 | 15,300,004 | +0.3 % |
| | nascent | 680 | 0 | n/a |
| | annotated | 14,648,896 | 14,699,996 | −0.3 % |
| **unstranded × ON** | **gDNA** | **14,782,077** | 15,300,004 | **−3.4 %** |
| | **nascent** | **157,840** | 0 | n/a |
| | **annotated** | **15,003,036** | 14,699,996 | **+2.1 %** |
| ⛔ g00 control | — | — | — | not measurable (held out) |

⭐ A perfect prior takes the blind stratum from **−78.6 % → −3.4 %** on gDNA and **+69.6 % →
+2.1 %** on annotated, i.e. into line with the stranded twin.

### 1c. ⚠ The residual after a perfect prior is CONCENTRATED AT `g05`, not spread

| condition | gDNA Δ% | nascent est |
|---|---:|---:|
| **g05 ss0.50 on** | **−56.4 %** | **149,220** |
| g50 ss0.50 on | −2.2 % | 3,800 |
| g98 ss0.50 on | −1.3 % | 4,820 |
| g05 ss0.99 on | −1.5 % | 9,094 |

⛔ **A perfect prior fixes `g50` and `g98` on the blind stratum and does NOT fix `g05`.** The
largest nascent false positive anywhere in the `oracle` arm is `g05 ss0.50 on` at 149,220 — 6×
the next. ⚠ This is the rung that only exists because `suite_resolves.py` requirement (f) demanded
it; on a `{g00, g50, g98}` panel this residual would have been invisible.

### 1d. ⚠ The nascent leak INVERTS between the zero control and the contaminated rungs

| condition | nascent est (truth 0) |
|---|---:|
| g00 ss0.99 **off** | 224,125 |
| g00 ss0.99 **on** | 328,333 |
| g00 ss0.50 off / on | 5,626 / 5,559 |
| g50 ss0.50 **on** | 585,122 |
| g98 ss0.50 **on** | 1,117,305 |
| g50 / g98 ss0.99 on | 15,323 / 18,147 |

At `g00` the leak is on the **STRANDED** conditions and ~40× the unstranded ones. At `g50`/`g98`
it is the **UNSTRANDED × capture-ON** cells. ⭐ These are plausibly TWO distinct absorption paths,
and only the second is explained by the `rho_ref` collapse. A repair aimed at the second should
not be expected to close the first.

---

## 2. ⭐⭐⭐ THE CEILING, ON THE PANEL THAT EXISTS

⛔ Every ceiling number ranking this project was measured on the 36-condition ladder DELETED on
2026-08-13. These replace them.

### 2a. TRANSCRIPT level — Σ|count_est − count_true|

| stratum | base | oracle | ratio | documented (retired ladder) |
|---|---:|---:|---:|---:|
| ALL (g00 excl.) | 20,755,019 | 10,827,911 | **0.522** | ~0.678 |
| stranded × OFF | 2,333,871 | 2,384,316 | 1.022 | 1.031 |
| stranded × ON | 2,907,697 | 2,722,525 | 0.936 | 0.981 |
| unstranded × OFF | 2,132,084 | 2,143,886 | 1.006 | 1.013 |
| **unstranded × ON** | 13,381,367 | 3,577,184 | **0.267** | **0.405** |

### 2b. GENE level

| stratum | base | oracle | ratio |
|---|---:|---:|---:|
| ALL (g00 excl.) | 12,008,407 | 2,123,691 | **0.177** |
| stranded × OFF | 144,015 | 173,870 | 1.207 |
| stranded × ON | 907,895 | 721,029 | 0.794 |
| unstranded × OFF | 151,498 | 158,202 | 1.044 |
| **unstranded × ON** | 10,804,999 | 1,070,590 | **0.099** |

⭐ **The documented PATTERN survives exactly** — neutral-to-worse on three strata, transformative
on one — **and the MAGNITUDE is larger than recorded**: 0.267 not 0.405 at transcript level, and
at gene level a perfect prior removes **90 %** of the blind stratum's error.

⚠ Two things a perfect prior does NOT buy: false-NEGATIVE mass gets **worse** (1.231 overall,
1.423 on stranded × OFF), and TPM error only reaches 0.756. It is not a free win everywhere.

⛔⛔ **THIS CEILING IS A LOWER BOUND, AND THE REASON IS STRUCTURAL.** `install_arm`
(`quant_accuracy.py:199, :223`) wraps **`PRIORS.assemble_priors` ONLY** and never touches
`CAL.calibrate` — unlike the `alloc_*` arms, which patch both (`:496`). So `pipeline.py:521`
builds `effective_lengths_em` from the **SHIPPED, UNMODIFIED** `CalibrationResult` in every oracle
arm. Whatever is wrong with that ruler is still installed in the 0.267 / 0.099 above.

---

## 3. ⭐⭐⭐ THE PRIOR IS ALL-OR-NOTHING — no single field is a lever

`_ARM_FIELDS` (`quant_accuracy.py:146-154`) injects, per arm, which `LocusPriors` fields come from
truth. GENE level, **unstranded × capture-ON**:

| arm | injects | error | ratio |
|---|---|---:|---:|
| base | — | 10,804,999 | 1.000 |
| `oracle_gdna` | `gdna_prior_count` | 9,636,829 | 0.892 |
| `oracle_rna` | `rna_prior_count` | 9,409,202 | 0.871 |
| `oracle_efflen` | `gdna_eff_len` | 9,552,460 | 0.884 |
| **`oracle`** | **all three** | **1,070,590** | **0.099** |

⭐⭐ **Every field alone buys 11–13 %. All three together buy 90 %.** The interaction is the whole
story; the noise floor is 0.996–1.013, so none of this is sampling.

⛔ **And a PARTIALLY-correct prior is actively HARMFUL on the healthy strata.** `oracle_efflen`
alone, gene level: **1.389** (stranded × OFF), 1.066 (stranded × ON), **1.375** (unstranded ×
OFF). Injecting one truth value beside two wrong ones is worse than three consistently-wrong ones.

⭐⭐⭐ **THIS REFRAMES THE GRAVEYARD.** `ROADMAP.md` §4.1's eleven refused mechanisms were all
single-lever injections into `assemble_priors` or `build_node_init` — one lever, eleven times.
This measurement says a single lever **cannot** work on this stratum however good it is, because
the three prior quantities must be MUTUALLY CONSISTENT. That is one structural explanation for
eleven failures rather than eleven separate stories, and it is a falsifiable prediction: any new
mechanism that moves one of the three numbers in isolation should score ≈0.88 on this stratum.

---

## 4. ⭐⭐⭐ CALIBRATION'S OWN ENDPOINT — `prior_vs_oracle.py`, 16 conditions, 1,385 s

`LocusPriors` is what the EM actually reads. **P** = shipped, **O** = the same assembler fed the
origin-split truth masses. ⚠ Undrained (`--no-drained-arm`), so the drain caveat is unpriced here.

### 4a. ⛔ The gDNA prior is PINNED on exactly two cells

| condition | true f_g | O_g | P_g | **P/O** |
|---|---:|---:|---:|---:|
| g05 ss0.50 off | 0.0516 | 240,819 | 217,993 | 0.905 |
| **g05 ss0.50 on** | 0.0511 | 492,785 | 322,499 | **0.654** |
| g05 ss0.99 on | 0.0511 | 492,745 | 468,811 | 0.951 |
| g50 ss0.50 off | 0.5085 | 2,404,230 | 2,323,520 | 0.966 |
| ⛔ **g50 ss0.50 on** | 0.5060 | 4,890,926 | 223,939 | **0.046** |
| g50 ss0.99 on | 0.5059 | 4,890,735 | 4,779,609 | 0.977 |
| g98 ss0.50 off | 0.9807 | 4,698,859 | 4,572,464 | 0.973 |
| ⛔ **g98 ss0.50 on** | 0.9805 | 9,548,775 | 434,347 | **0.045** |
| g98 ss0.99 on | 0.9805 | 9,548,137 | 9,384,967 | 0.983 |

⭐ Every other condition sits at **0.905–0.983**. The two catastrophic cells sit at **0.045–0.046**
— *pinned*, and at essentially the same value despite a 2× difference in true `f_g`. This
reproduces the documented "`P/O = 0.040`" pinning, now localised to `ss0.50 × capture_on` only.
⚠ `g05 ss0.50 on` at 0.654 is degraded but NOT pinned — the same cell where a perfect prior still
leaves −56.4 % (§1c). Three regimes, not two.

### 4b. ⛔ `gdna_eff_len` is 4–8× WRONG on the blind stratum

| stratum | weighted rel err | median rel |
|---|---:|---:|
| stranded × OFF | 0.0082 | 0.0081 |
| stranded × ON | 0.0183 | 0.0062 |
| unstranded × OFF | 0.0132 | 0.0069 |
| ⛔ **unstranded × ON** | **8.2098** | **4.3258** |
| g00 control | nan | 0.9664 |

⭐ **This is the `rho_ref` defect measured on calibration's own endpoint**, independent of the EM:
`gdna_eff_len` (`priors.py:347`, computed from `_global_reference_density`) is off by a MEDIAN
FACTOR of 4.3 on unstranded × capture-ON and by <2 % everywhere else.

### 4c. The composition claim `phi = a_g/(a_g+a_r)`

| stratum | mwae_phi | median |
|---|---:|---:|
| stranded × OFF | 0.0110 | 0.0220 |
| stranded × ON | 0.0127 | 0.0249 |
| unstranded × OFF | 0.0143 | 0.0302 |
| ⛔ **unstranded × ON** | **0.5683** | **0.4647** |
| g00 control | 0.0750 | 0.0581 |

⭐ The SCALE is exact on every stratum (`log10 P/O = +0.000`) — the assembler's arithmetic is
right and it is the COMPOSITION that is wrong, on one stratum, by ~0.57.

### 4d. ⚠ A SEPARATE, UNIFORM SCALE ERROR IN `gdna_eff_len` — present everywhere, including g00

Table ⑩ asks whether `gdna_eff_len`'s clamp is a genomic extent or an INCIDENCE sum:

| stratum | support/genomic | Σ support | Σ genomic |
|---|---:|---:|---:|
| stranded × OFF | **1.15** | 148,942,857 | 139,628,517 |
| stranded × ON | **1.17** | 150,860,199 | 139,291,592 |
| unstranded × OFF | **1.15** | 148,946,348 | 139,632,903 |
| unstranded × ON | **1.18** | 150,827,646 | 139,269,875 |
| g00 control | **1.12** | 184,280,257 | 172,534,934 |

⛔ The instrument's own note: *"`support/genomic` well above 1 means every interior boundary is
adding ~mu_g − 1 to the locus's clamp. The EM divides the gDNA component's abundance by
`gdna_eff_len`, so this is a direct scale error on a shipped number."* ⭐ **This one is NOT
stratum-specific** — it is 12–18 % on every stratum including the zero control, so it is a
different defect from §4b and is invisible to any blind-stratum experiment.

---

## 5. What this says about the next priority

⭐ **The target is a COHERENT prior on unstranded × capture-ON.** It is worth 90 % of that
stratum's gene-level error, the stratum carries 90 % of the panel's gene-level error, and the
noise floor is 0.4 %.

⭐⭐ **AND THE CALIBRATION ENDPOINT NOW NAMES THREE SEPARABLE DEFECTS RATHER THAN ONE**, which is
the thing §3's all-or-nothing result was missing. They are independent and only one is
stratum-specific:

| defect | where | measured |
|---|---|---|
| **`pinned-gdna-prior`** | `g50`/`g98` × unstranded × capture-ON ONLY | `P/O` 0.045–0.046 against 0.905–0.983 everywhere else (§4a) |
| **`eff-len-on-the-blind-stratum`** | unstranded × capture-ON | weighted rel err **8.21**, median **4.33**, against <0.02 elsewhere (§4b) |
| **`clamp-is-an-incidence-sum`** | **EVERYWHERE, including g00** | `support/genomic` **1.12–1.18** on every stratum (§4d) |

⛔ **`clamp-is-an-incidence-sum` is the one nothing has ever been able to see**, because it is uniform: every blind-stratum
experiment, every ceiling, and every zero control carries it identically, so it cancels out of
every contrast that has ever been run. It is a direct 12–18 % scale error on a shipped number the
EM divides by. ⭐ It also satisfies the two admission criteria the graveyard's eleven all failed:
its correct value is DERIVABLE from the definition (a genomic extent is not an incidence sum), and
the repair REMOVES a fabricated contribution rather than adding a prior.

⚠ **`pinned-gdna-prior` and `eff-len-on-the-blind-stratum` are plausibly the same defect seen twice** — both are confined to unstranded ×
capture-ON, and `_global_reference_density` feeds both `gdna_eff_len` (`priors.py:347`) and
the `effective_lengths_em` ruler. ⛔ But that is a HYPOTHESIS: `oracle_efflen` fixes only the
prior-side half and buys 0.884, so if the two were one defect with one cause, fixing that half
should have bought more. Testing it needs the arm that does not exist (§below).

⭐ **`g05 ss0.50 on` is a THIRD regime and the panel only has it by accident.** `P/O` = 0.654 (not
pinned, not healthy), and it is the one cell a PERFECT prior fails to fix (−56.4 %, §1c). It
exists only because `suite_resolves.py` requirement (f) forced the rung in.

⛔ **Per-transcript weighting (the work that was in flight) cannot move this at first order** —
`apply_grouped_prior_update` skips `gdna_index`, so no per-transcript weight touches the library
gDNA fraction.

⚠ **`capture_eff_length` / `_global_reference_density` is NOT ranked by this measurement, and that
is a gap rather than an acquittal.** `oracle_efflen` prices only the PRIOR-side use of `rho_ref`
(`priors.py:347`) and buys 0.884. The OTHER consumer — `capture_eff_length.py:272` →
`effective_lengths_em` — is wrong in EVERY arm here including `oracle`, and **no instrument that
exists can price it.** Building that arm is the cheapest way to close the question.

⭐ **The simplex closure defect is real, is 36–52× concentrated on this exact stratum, and is not
the cause of the blindness.** Mass-weighted deficit: 8.76 % (REGION) / 9.61 % (BOUNDARY) on
unstranded × capture-ON against 0.168 % / 0.267 % on stranded × capture-ON. But it is ~9 % at ALL
FOUR gDNA levels while the blindness only appears at `g50`/`g98`; r ≈ 0.53. Shared upstream cause,
not the same defect. ⛔ Its root cause is **`f_g` is a posterior MEDIAN while `f_pos`/`f_neg` are
posterior MEANS** (`simplex_logodds.py:460-463, :615-625`), so `s = median(f_g) + 1 − mean(f_g)`,
which closes only for a symmetric posterior — this CORRECTS `ROADMAP.md` §0b, which attributes it
to `sweep.py`'s independent clip and itself notes the clip cannot be the whole cause.

⚠ **`ROADMAP.md` §0b's headline is measuring EMPTINESS.** "74.72 % of objects close" — every
zero-mass slot closes trivially; among objects carrying library mass the closure rate is **0 % at
g00** and 5.8–8.5 % elsewhere. Only the mass-weighted view is informative.
