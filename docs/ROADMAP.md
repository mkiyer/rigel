# ROADMAP — where the tool is, and what to do next

Reading order for a new session: **`SUCCESS.md`** (how performance is measured) → **this file** (current
state and priority) → **`TRAPS.md`** (mistakes not to repeat) → `DESIGN.md` / `EQUATIONS.md` as needed.

⛔⛔ **THIS FILE POINTS FORWARD. IT IS NOT A CHANGELOG.** A closed item earns **one line** in §3, and only
because that line is what stops it being rebuilt. The derivation goes to `EQUATIONS.md`, the lesson to
`TRAPS.md`, the ruling to `DESIGN.md`, and the investigation stays in git. ⚠ This file reached 1,030 lines
before its first prune (2026-08-05) by accumulating a section per campaign; if a section here reads like a
report on work already done, delete it.

⛔ **Every number below was measured on the current tree and the current panel.** Re-derive rather than
trust — a number that has moved is a result, not a documentation bug.

---

## §0 THE STATE

⛔ **Re-derive rather than trust.** A number that has moved is a result, not a documentation bug.

⚠ This was "six numbers" and is now ten, because the 2026-08-07 validation campaign measured four
things nothing had measured before. If it grows again, prune it rather than letting it become a report.

| | | |
|---|---|---|
| **Stage A — the accumulator** | ✅ **DONE**, and that is a measurement | perfecting BOTH fragment-length models is worth **2.6 %** of the deliverable, down from 22.2 % |
| **calibration, 3 of 4 strata** | ✅ median library `f_gdna` error **0.005–0.012**, and ⭐ the PRIOR the EM reads is within **2.5–4.6 %** of a perfect one | stranded × on/off and unstranded × capture-OFF |
| ⛔ **calibration, unstranded × capture-ON** | ⛔ **BLIND** — reports **0.033–0.058** while truth spans **0.00 → 0.98**, and hands the EM a gDNA prior **94.4 %** short | not noisy; a flat line. This is the whole open problem |
| ⛔ **the prior ASSEMBLER, under capture** | ⛔ **+15.1 % over-call** (`rel` 0.179) with PERFECT masses in; **0.005** off capture | ⭐ NEW 2026-08-07. Larger than calibration's own error on the three strata it handles. §1.1 |
| **message propagation** | ⛔ **OFF** since 2026-08-07 (`config.message_propagation`) | net better on 3 of 4 strata (−58 / −44 / −32 %); +155 % on the fourth, which carries 73 % of panel error |
| **the price of that** | ⚠ zero-gDNA golden scenarios go **0.029 → 89.93** and **0.005 → 9.58** | both AMBIG loci — the stratum above |
| **end to end — the LIBRARY figure** | ✅ mean `\|f_gdna − truth\|` **0.1060**, and a perfect prior takes it to **0.0097** | calibration is the whole bottleneck here. §1.3 ① |
| ⛔ **end to end — TRANSCRIPT assignment** | ⛔ **31.1 %** of fragments misassigned, and a perfect prior removes only **32 %** of that | 67.5 % survives, and on capture-OFF a perfect prior is 3–4 % WORSE. §1.3 ② |
| ⛔ **the NASCENT channel** | ⛔ **20.2 M fragments** parked on entities whose truth is exactly **0** | 9.2 % of all true RNA, invisible in every transcript table. §1.2 |
| **reproducibility** | ⛔ **the tool does not reproduce itself**: `EMConfig.seed` defaults to `None` | `TRAPS: the-deliverable-is-not-reproducible-by-default`. Counts still conserve exactly (`mrna+nrna+gdna == n_unambig+n_em`) |

⭐⭐ **THE ONE SENTENCE, AND IT NOW HAS TWO HALVES BECAUSE THE TOOL HAS TWO DELIVERABLES.**

**On the LIBRARY gDNA fraction, the tool is accurate everywhere except unstranded × capture-ON, where it
cannot see gDNA at all** — `κ = ½` makes the strand λ-term identically 0 and no other channel reaches an
AMBIG slot. That is the whole of the calibration problem, and §1.3 ① proves it is the *whole* problem for
this deliverable: a perfect prior fixes it completely. Everything in §2 is about giving that slot its own
evidence.

⛔ **On TRANSCRIPT-LEVEL assignment, that sentence is FALSE and was never tested until 2026-08-07.** The
tool misassigns 15.6–20.6 % of fragments even on the three strata calibration handles perfectly, and
42.9 % at `g00` where there is no gDNA at all to get wrong. A perfect prior removes ~32 % of the panel
total and makes two strata *worse*. **These are two different problems in two different files, and
conflating them is how "is calibration the bottleneck?" got asked for months without an answer.**

⭐ **And the blindness now has a second, sharper statement, in the EM's own units.** The gDNA prior at
unstranded × capture-ON is not merely wrong, it is **PINNED**: `P/O = 0.040` at every one of
`g50 / g75 / g90 / g98`, with `P_g` sitting at 221 k → 428 k while the truth `O_g` runs 5.5 M → 10.8 M.
A flat line in the deliverable was already known; this is the same flat line one stage earlier.

⛔ **Read `mwae_all` / `Σ|err|`, never `solv%` / `mwae` / `conf-wrong` / `calib`** — those four share a
denominator the solver moves by declining to answer. And quote the SHIPPED column, not pass-0: a −37.2 %
pass-0 win was −3.9 % shipped (TRAPS: the-intermediate-is-not-the-deliverable).

## §1 ⭐⭐⭐ WHAT TO DO NEXT — **VALIDATE THE WHOLE TOOL BEFORE ADDING ANYTHING**

⛔⛔ **OWNER RULING 2026-08-07: no new features until the tool is characterised end to end.** All of the
work below runs with **message propagation OFF** and **`length_likelihood` OFF** — the two switches stay
down so that every number is attributable to the tool as it stands, not to a feature landing mid-campaign.

⭐ The four items are ONE campaign and they share a substrate. Do them in order: each answers a question the
next one needs.

| | the item | the question it answers | state |
|---|---|---|---|
| ⭐⭐⭐ **1** | **calibration-prior-vs-oracle** | Calibration's endpoint is a PRIOR for the EM — three per-locus scalars. **How wrong are they?** | ✅ **DONE** — `prior_vs_oracle.py` (+14 gates, 10/10 perturbations). **§1.1** |
| ⭐⭐⭐ **2** | **tool-absolute-accuracy** | The panels are simulated with ground truth. **What is the tool's transcript-level accuracy, end to end?** | ✅ **DONE** — `quant_accuracy.py --arm base` (+11 gates, 8/8 perturbations). **§1.2** |
| ⭐⭐⭐ **3** | **error-downstream-of-calibration** | Inject the ORACLE prior and re-quantify. **How much of the error is calibration's, and how much is the EM's?** | ✅ **DONE** — `--arm oracle`, `noop` byte-identical, noise floor −0.01 %. **§1.3** |
| ⭐⭐ **4** | **performance** | It has been a while, and new code means slowdown. | ⏸ **NOT STARTED** — ⛔ profile on real cfRNA, not on this panel. §1.0 |

⛔⛔ **ITEMS 1–3 ARE DONE (2026-08-07). THE ANSWER IS BOTH, AND THE QUESTION HAD A HIDDEN ASSUMPTION.**
The framing was "calibration or the EM?" and it presumed one deliverable. There are two, and they
answer differently:

* **the library gDNA/RNA separation** — calibration is the whole bottleneck. A perfect prior takes it
  from 0.1060 to **0.0097** and both zero controls to exactly right (§1.3 ①). `length_likelihood` is
  therefore confirmed as the right next build *for this deliverable*.
* **transcript-level assignment** — calibration is **not** the bottleneck. 67.5 % of the error survives
  a perfect prior, and on the two capture-OFF strata a perfect prior is **worse** by 3–4 %, at 57–1,129×
  the noise floor (§1.3 ②). That is a defect in how the prior is CONSUMED and it lives in the EM.

⭐ So §2's `length_likelihood` work is unblocked and correctly aimed — but it should be scored on the
library figure, and a second, separate line of work exists in the EM that nothing was previously
pointing at.

### §1.1 ⭐⭐⭐ ITEM 1 IS DONE — CALIBRATION'S ENDPOINT, MEASURED AGAINST THE ORACLE

`scripts/design/prior_vs_oracle.py`, all 36 ladder conditions, messages OFF, `length_likelihood` OFF,
noop gate byte-identical on all 36 × 3 arrays. `LocusPriors` is what the EM reads; **P** is the shipped
prior, **O** is the same assembler fed the origin-split truth masses, **F** is the direct per-locus
fragment count from `node_start_count` (EXACT on the gDNA arm — gDNA does not splice).

**① `P − O`, calibration's own error, gDNA arm, as a fraction of the true prior:**

| stratum | `rel` | `mwae_phi` | | stratum | `rel` | `mwae_phi` |
|---|---|---|---|---|---|---|
| stranded × capture OFF | 0.027 | 0.0099 | | unstranded × capture OFF | 0.046 | 0.0166 |
| stranded × capture ON | 0.025 | 0.0139 | | ⛔ unstranded × capture **ON** | **0.944** | **0.5270** |

36.8 M of the panel's 39.0 M misplaced prior fragments — **94 %** — are that one stratum, and `mwae_phi`
there climbs monotonically with the gDNA level: 0.0345 (`g01`) → 0.3083 (`g25`) → 0.5890 (`g50`) →
0.9486 (`g98`).

⭐⭐ **THE SCALE IS RIGHT EVERYWHERE AND ONLY THE SPLIT IS WRONG, WHICH IS NEW AND IS ACTIONABLE.**
`log10(scale_P / scale_O)` is **+0.000 to +0.003** on every stratum including the broken one: the prior's
total pseudocount `a_g + a_r` is correct to a fraction of a percent, so the EM is being pulled with the
right STRENGTH in the wrong DIRECTION. Nothing here is a regularisation or a prior-weight problem, and a
mechanism that changes how hard the prior pulls is aimed at a number that is already right.

**④ `O − F`, the ASSEMBLER's own error** — perfect masses in, true fragment counts as the reference, so
this prices the mass → density → fragment-count conversion, the overlap projection and the
support-weighted pooling alone:

| | capture OFF | capture ON |
|---|---|---|
| gDNA arm `rel` | **0.005** (both strata) | ⛔ **0.179**, net **+5.1 M on 33.9 M (+15.1 %)** |

⛔ **That is bigger than calibration's own error on the three strata calibration handles (0.025–0.046).**
So on those strata the assembler is now the dominant prior error, and it appears only under capture —
consistent with the second-order residual `assemble_priors`' own docstring flags, since `Σm/ΣS` is the
support-weighted mean density and capture puts a strong gradient *inside* a locus. ⚠ Diagnosed, not fixed.

⛔ **The `g00` ZERO-gDNA control fails at the prior, and its worst stratum is NOT the blind one.** The
shipped prior claims **2,067,637 gDNA fragments in libraries containing none**:

| `g00` condition | false gDNA prior | `mwae_phi` |
|---|---|---|
| ⛔ ss0.50 capture **OFF** | **1,707,321** | 0.2524 |
| ss0.50 capture ON | 312,200 | 0.0452 |
| ss0.99 capture OFF | 28,843 | 0.0043 |
| ss0.99 capture ON | 19,274 | 0.0028 |

Unstranded × capture-OFF is worst by **5.5×** — and that is a stratum that reads healthy (`rel` 0.046) on
every contaminated row. A zero control was the only thing that could have found it.

⭐ **Both `flgap` panels reproduce the ladder exactly**, so a ±40 % gDNA−RNA fragment-length gap changes
nothing while `length_likelihood` is off: three strata at `rel` 0.024–0.036, unstranded × capture-ON at
**0.966** (`flgap_short`) and **0.837** (`flgap_long`).

⚠ **Undrained on both sides, and priced rather than waved**: the drain moves the shipped prior by
**0.153 %** (gDNA) and **0.462 %** (RNA). A drained oracle is inadmissible on this panel —
`TRAPS: an-equal-length-panel-defeats-the-lift`.

### §1.2 ⭐⭐⭐ ITEM 2 IS DONE — THE TOOL'S ABSOLUTE ACCURACY, END TO END

`scripts/design/quant_accuracy.py --arm base`, all 36 ladder conditions, seed pinned. Scored **count
against count** against each condition's OWN `truth_abundances.tsv` (the realised observed fragment
count). ⚠ The suite-level `truth_abundances_nrna_none.tsv` is the *pre-capture molar* abundance — a
different quantity, log-space correlation **0.72** with the realised one under capture — and scoring
against it would charge the tool for hybrid capture it never claims to invert.

| stratum | true fragments | misassigned `Σ\|Δ\|` | **rate** | false-positive mass | Spearman | MARD |
|---|---|---|---|---|---|---|
| ALL (g00 excluded) | 178,399,996 | 55,488,787 | **0.311** | 16,162,221 (0.091) | 0.825 | 0.393 |
| stranded × capture OFF | 44,599,999 | 6,936,270 | 0.156 | 1,610,144 (0.036) | 0.910 | 0.286 |
| stranded × capture ON | 44,599,999 | 8,811,983 | 0.198 | 1,454,490 (0.033) | 0.778 | 0.448 |
| unstranded × capture OFF | 44,599,999 | 9,205,612 | 0.206 | 1,776,020 (0.040) | 0.896 | 0.315 |
| ⛔ unstranded × capture ON | 44,599,999 | 30,534,922 | **0.685** | 11,321,567 (0.254) | 0.716 | 0.523 |
| ⛔ **`g00`, ZERO gDNA** | 40,000,000 | 17,149,596 | **0.429** | 2,610,770 (0.065) | 0.816 | 0.374 |

⛔⛔ **READ THE `g00` ROW FIRST, AND IT IS THE MOST IMPORTANT NUMBER ON THIS PAGE.** Those four
conditions contain **no gDNA at all** — the library is 100 % RNA — and the tool still misassigns
**42.9 %** of its fragments, worse than three of the four contaminated strata. At the two
`g00 ss0.99` rows the gDNA channel is *essentially silent and essentially correct*
(`gdna_est` = 7,683 and 10,754 against a truth of 0, i.e. ~0.1 % of the library) and **38–49 % of
fragments are still misassigned.** ⭐ That error cannot be calibration's, because at those conditions
calibration is already right. It is transcript assignment.

⭐⭐ **AND THERE IS A THIRD FALSE-POSITIVE CHANNEL THAT NOTHING WAS WATCHING.** The panel is
`nrna_none`: the true nascent count is **exactly 0** at every condition. The EM parks **20,151,758
fragments — 9.2 % of all true RNA — on synthetic nascent entities.** ⛔ `get_counts_df` DROPS the
synthetics, so this is invisible in every transcript-level table ever printed here; it shows up only
because `quant_accuracy.py` reports it on the `library` row.

| stratum | `nrna_est` (truth 0) | as a share of true RNA |
|---|---|---|
| stranded × capture OFF | 263,818 | 0.006 |
| unstranded × capture OFF | 816,933 | 0.018 |
| stranded × capture ON | 1,398,112 | 0.031 |
| ⛔ unstranded × capture ON | **13,373,410** | **0.300** |
| ⛔ `g00` (no gDNA, no nascent) | 4,299,485 | 0.107 |

⭐ **The blind stratum's missing gDNA is going into the NASCENT channel, and that is mechanically
exactly what should happen.** gDNA is genomically continuous and unspliced; so is nascent RNA. When the
prior cannot say "gDNA" (§1.1: `P/O = 0.040` there), the only remaining unspliced-compatible hypothesis
is the nascent entity, and 13.4 M fragments take it. ⚠ This is an EM assignment channel, not a
calibration population — AXIOM 0 is about the deconvolution, and the EM's mature/nascent split is
downstream of it.

⚠ **A constant, negligible floor**: 5 transcripts in each condition's truth table are absent from the
index (the transcript filter drops them), carrying 2–327 fragments. Identical in every arm.

### §1.3 ⭐⭐⭐ ITEM 3 IS DONE — AND THE ANSWER IS **BOTH, ON TWO DIFFERENT DELIVERABLES**

`scripts/design/quant_accuracy.py`, `base` vs `oracle` (the §1.1 prior injected in place of
calibration's, rebuilt in-process on the run's own loci), all 36 conditions, seed pinned.

⛔ **The gate passed first**: `arm_identity.py qa_base qa_noop` is **byte-identical on all 1,296 scored
fields of all 72 rows** — the wrapper builds O in full and discards it, so the whole injection path is
proven inert. And the noise floor is **−0.01 %** (`base_reseed`, seed + 1, nothing else varied), so
every effect below is 54–4,700× the size of a different draw.

**① THE LIBRARY-LEVEL SEPARATION IS CALIBRATION-LIMITED, AND A PERFECT PRIOR SOLVES IT OUTRIGHT.**
Mean `|f_gdna − truth|`, intergenic included (the `cli.py` denominator):

| stratum | base | oracle |
|---|---|---|
| stranded × capture OFF | 0.0036 | 0.0028 |
| stranded × capture ON | 0.0242 | 0.0158 |
| unstranded × capture OFF | 0.0075 | 0.0018 |
| ⛔ unstranded × capture ON | **0.3885** | **0.0183** |
| ⛔ `g00` zero control | 0.0374 | **0.0000** |
| **ALL (g00 excluded)** | **0.1060** | **0.0097** |

**10.9× better, and the zero control goes to exactly right** — `g00 ss0.50 capture_off` falls from
**1,263,494 phantom gDNA fragments to 42.** Nothing survives there. For the number this tool exists to
produce, calibration is confirmed as the whole bottleneck.

**② TRANSCRIPT ASSIGNMENT IS NOT, AND ON TWO STRATA A PERFECT PRIOR MAKES IT WORSE.** Misassigned
fragments `Σ|count_est − count_true|`:

| stratum | base | oracle | effect | noise floor |
|---|---|---|---|---|
| stranded × capture OFF | 6,936,270 | 7,161,996 | ⛔ **+3.25 %** | +0.00 % |
| stranded × capture ON | 8,811,983 | 8,544,966 | −3.03 % | −0.06 % |
| unstranded × capture OFF | 9,205,612 | 9,606,938 | ⛔ **+4.36 %** | +0.08 % |
| unstranded × capture ON | 30,534,922 | 12,165,799 | **−60.16 %** | −0.02 % |
| `g00` zero control | 17,149,596 | 16,320,735 | −4.83 % | −0.01 % |
| **ALL (g00 excluded)** | **55,488,787** | **37,479,699** | **−32.46 %** | −0.01 % |

⭐⭐ **THE TWO POSITIVE ROWS ARE THE RESULT, AND THEY ARE NOT NOISE — THEY ARE 1,129× AND 57× THE
FLOOR.** A *more accurate* prior costs the two capture-OFF strata 3–4 % of their transcript accuracy,
and the mechanism is visible in the directional split: false-NEGATIVE mass rises **+26.0 %** and
**+23.6 %** on exactly those two strata while false-positive mass is flat. **The prior is correct at
the LOCUS and the EM spends it indiscriminately WITHIN the locus** — told (truthfully) that a locus
holds more gDNA, the EM removes it from whichever transcripts its likelihood ranks lowest, and off
capture those are partly real RNA. ⛔ That is a defect in how the prior is *consumed*, not in the
prior, and no improvement to calibration can reach it.

⭐ **67.5 % of the transcript error survives a perfect prior**, and on the three strata where
calibration already works the prior is worth ±4 % — i.e. nothing.

**③ AT GENE LEVEL — ISOFORM AMBIGUITY SUMMED AWAY — MOST OF THE RESIDUE DISAPPEARS, AND THE PRIOR
MATTERS TWICE AS MUCH.** Summing a gene's isoforms collapses exactly the error that comes from not
knowing WHICH isoform a fragment came from; what survives is error in deciding whether the fragment
was RNA from this gene at all, which is the question Rigel is for.

| stratum | transcript | **gene** | gene ÷ tx | gene, as a share of true RNA | gene, oracle arm |
|---|---|---|---|---|---|
| stranded × capture OFF | 6,936,270 | **311,100** | 4.5 % | **0.7 %** | 532,652 (⛔ 1.71×) |
| unstranded × capture OFF | 9,205,612 | 1,033,478 | 11.2 % | 2.3 % | 1,389,836 (⛔ 1.35×) |
| stranded × capture ON | 8,811,983 | 2,636,387 | 29.9 % | 5.9 % | 2,373,626 (0.90×) |
| ⛔ unstranded × capture ON | 30,534,922 | 22,951,782 | 75.2 % | **51.5 %** | 4,766,115 (**0.21×**) |
| **ALL (g00 excluded)** | 55,488,787 | **26,932,747** | 48.5 % | 15.1 % | **9,062,229 (0.34×)** |
| ⛔ `g00` zero control | 17,149,596 | 7,053,026 | 41.1 % | 17.6 % | 5,969,001 (0.85×) |

⭐⭐ **95.5 % of the transcript error on stranded × capture-OFF is isoform ambiguity** — 6.94 M falls to
311 k, which is **0.7 %** of the true fragments. That is an ordinary quantification result and there is
no Rigel-specific defect there to chase. ⛔ **And a perfect prior removes 66 % of the gene-level error
against 32 % at transcript level**, so calibration matters roughly twice as much as the transcript
number suggested. ⚠ The two capture-OFF strata still get **WORSE** with a perfect prior, and more so
at gene level (1.71× and 1.35×) — §1.3 ②'s finding survives the harder test.

⛔⛔ **`g00` gene-level is 7.05 M with NO gDNA present, and a perfect (all-zero) prior leaves 5.97 M of
it.** Nothing about calibration can reach that. It is the nascent channel below.

⭐⭐ **THE ONE PART THAT IS CLEARLY RIGEL'S AND IS UNTOUCHED BY CALIBRATION: THE NASCENT CHANNEL.**
§1.2's 20.2 M fragments on zero-truth nascent entities falls to 9.8 M under the oracle prior — but at
`g00`, where there is neither gDNA nor nascent RNA, it goes **UP**: 4,299,485 → 4,636,624. The
fragments the shipped prior was calling gDNA become *nascent* rather than mature. That channel is
downstream of the prior and a perfect prior cannot fix it.

### §1.0 ⚠ THE PERFORMANCE SUBSTRATE IS A TRAP, AND THE NUMBERS BELOW SHOW WHY

Measured 2026-08-07 on one 10 M-fragment ladder condition (a **35,135-node** chr22 index):

| stage | time | share |
|---|---|---|
| per-locus EM | **15.9 s** | **47 %** |
| native scan | 6.5 s | 19 % |
| calibration | 6.5 s | 19 % |
| second-pass drain | 3.5 s | 10 % |
| **total** | **33.5 s** | |

⚠ **And message propagation costs nothing measurable** — 33.7 s ON vs 33.5 s OFF, calibration 6.6 vs 6.5 s.
A hypothesis that the relay was the expensive part is REFUTED.

⛔⛔ **BUT THIS PROFILE IS UPSIDE DOWN FROM THE REAL ONE, AND BOTH ARE CORRECT.** Calibration cost is
**depth-INDEPENDENT** — every node in the index is solved regardless of read depth — so it scales with the
INDEX while the EM scales with the DATA. On real cfRNA at genome scale (~1.5 M nodes) calibration has
measured ~66 s against the EM's ~24 s: **the exact reverse of the table above.** ⛔ So optimising against
this panel would tune the wrong stage — `TRAPS: toys-rank-hotspots-backwards`, which cost a whole analysis
once already. **Profile on real cfRNA.** Both are on disk: `~/Downloads/rigel_runs/cfrna/` (4 samples,
including `mctp_vcap_rna20m_dna05m`) and the genome-scale index at `~/Downloads/rigel_runs/refs/rigel_index`.

## §2 ⛔ AFTER THE CAMPAIGN — NOT NOW, kept because the briefs are complete and should not be re-derived

⛔⛔ **BOTH ITEMS BELOW ARE ON HOLD BY OWNER RULING (2026-08-07) until §1 is finished.** They are written up
in full because the analysis is done and re-deriving it would be waste — not because they are next. ⚠ In
particular `length_likelihood` STAYS OFF: its panels and oracle caches are built and ready, and turning it
on mid-campaign would make every number in §1 unattributable.

### §2.1 THE BRIEF FOR **the-cancelling-pair** — everything that experiment needs, in one place

⭐ Promoted from a working doc when it was deleted (2026-08-07). This is the next build.

**THE TWO HALVES, AND WHY NEITHER HAS AN HONEST PRICE ALONE.** The certified-RNA channel is a LOWER BOUND
(`ρ_R(exon) ≥ ρ_ν(B) + ρ_μ(B)` — the exon may hold molecules that never touch that seam) and ψ delivers it
as a two-sided Gaussian, penalising a destination holding MORE RNA exactly as hard as one holding less.
Making it one-sided — `−½·p·max(0, mo − log f)²`, no new constant — is **the only mechanism the zero-gDNA
control has ever ENDORSED: −81.9 %, 8/8.** The panel is nonetheless **+5.2 %**, and the reason is
`TRAPS: a-cancelling-defect-pair`: the two-sided term was doing TWO jobs, the bound *and* — by its upward
side — a de-facto gDNA LEVEL channel. Removing the accident exposes the real gap, +7.7 % on exactly the
stratum with no working level channel and −9.2 % where the prior already suffices.

⛔⛔ **AND THE LEVEL CHANNEL IS STRUCTURALLY DISCONNECTED, WHICH IS A THEOREM AND NOT A MEASUREMENT.**
Receivers of a one-slot-step channel at `g00`, by slot type:

| channel | EDGE | intergenic NODE | intron NODE | exon NODE |
|---|---|---|---|---|
| gDNA level | 2,592 | 0 | 0 | **0** |
| certified RNA | **0** | 0 | 0 | 10,493 |

The chain is strictly `N E N E … N`, so a one-slot-step channel is **BIPARTITE**: a NODE emitter can only
reach an EDGE. The only licensed originators of a gDNA level are structurally pure-gDNA objects, which are
NODEs — so **no NODE can ever receive a gDNA level.** It is not weak, it is disconnected from every object
carrying mass. ⛔ Two repairs are already refused: letting an EDGE originate a level (a patch on a
symptom), and letting the level cross two steps (the chain-fused level is dominated by the 1,312 intergenic
anchors and hands every exon the off-probe floor, against a true neighbouring density **346×** higher).

⚠ **The measurement to run it against**: at `g01 ss0.50 capture_on`, ψ-boundary ablation with the identity
exact, every single channel ablation is small (+2.3 / +6.8 / +4.9 / −0.1 %) and the joint one is +60.7 % —
`TRAPS: all-small-singly-large-jointly`, so the four share an upstream quantity. And at HEAD's top-12 error
slots the self-solve WITH the fitted prior is nearly correct while the messages destroy it; muting the
certified-RNA channel alone recovers the truth at 8 of the 12.

### §2.2 ⛔⛔ THE LENGTH CHANNEL CANNOT BE PRICED ON THE LADDER — the substrate is BUILT AND CACHED

`length_likelihood` defaults **False** and is the **only** channel that can give an unstranded slot its OWN
composition evidence: it is θ-independent, so the Schur complement that zeroes the strand term at an AMBIG
node does not apply, and it enters `tau_lam` ungated. ⛔ But `TRAPS: equal-lengths-carry-no-composition` —
at equal component mean lengths the 2×2 is identified only through `μ_g − μ_r`, and the ladder was built
with equal configured lengths: the realised post-capture gap is **+1.5 %**. Enabling it there would measure
approximately nothing and that would be read as "the feature does not work".

⭐⭐ **THE PANELS ARE BUILT AND THE ORACLE CACHES EXIST (2026-08-07)**: `suite/flgap_short` (realised gap
**−41.0 %**) and `suite/flgap_long` (**+40.4 %**), 4 conditions each, 32 MB of cache each, verified
byte-identical on a cached re-run (2m39s cold → **21 s** warm). Nothing blocks the experiment but the
campaign in §1.

⭐ **And both panels reproduce the ladder's failure exactly, which is what makes them a clean A/B** — with
messages off and the length channel off, only unstranded × capture-ON is broken:

| | truth | reported | | | truth | reported |
|---|---|---|---|---|---|---|
| `flgap_short` ss0.50 OFF | 0.4854 | 0.4771 ✅ | | `flgap_long` ss0.50 OFF | 0.4933 | 0.4843 ✅ |
| ⛔ ss0.50 **ON** | 0.5288 | **0.0324** | | ⛔ ss0.50 **ON** | 0.5895 | **0.1048** |
| ss0.99 OFF | 0.4854 | 0.4785 ✅ | | ss0.99 OFF | 0.4932 | 0.4859 ✅ |
| ss0.99 ON | 0.5290 | 0.5147 ✅ | | ss0.99 ON | 0.5894 | 0.5764 ✅ |

A fix should move exactly those two rows and leave the other six alone.

## §3 ⛔ THE VERTEX PROBLEM — CLOSED, kept as one paragraph so it is not rebuilt

ψ lands 5–8 % short of `f_g = 1` on unexpressed genes, and that was ranked #1 for a campaign. It is
**not a build**: `vertex_ceiling.py` prices it and evidence-free objects sit at median `|Δ|/sd(f_g)` =
**z = 0.5–0.6**, inside their own 1σ — the 24.4 % it measures is *the value of missing information, not
headroom*. No prior-free solve can beat it: every proper prior on [0,1] has a median strictly inside
(0,1). ⭐ `test_certified_rna_licence.py`'s **zero-spliced-count** gate independently closed the one channel that might have
supplied vertex evidence — a zero certified count is consistent with `f_g = 1` too. ⚠ Its companion
hypothesis, the **phantom-gDNA floor**, turned out to be a DIFFERENT and much larger bug and is now
fixed: TRAPS: a-zero-count-is-a-measurement/TRAPS: a-ratio-cannot-carry-zero, §0's 39 %.

## §4 ⛔ WHAT IS DELIBERATELY NOT NEXT — one line each, so it is not rebuilt

| | closed by | verdict |
|---|---|---|
| **the gDNA scale rule** · **the mass pin** · **TSS/TES as the population licence** | landed 2026-08-04 | ✅ `EQUATIONS.md` §3.5/§3.5b/§3.5c, gates in `test_gdna_scale_rule.py`, `test_relay_mass_pin.py`, `test_terminus_population_licence.py`. ⚠ The ceiling says the mass pin cost the panel **nothing** (+0.0002 to delete it outright); it landed on the derivation and on being free |
| **face (I) of the `intron\|exon` EDGE** | re-solve ceiling + panel arm | ⛔ **DO NOT BUILD.** The derivation (`EQUATIONS.md` §3.6) is re-verified and is not what failed: handing both EDGEs the ORACLE truth and re-solving is worth **−0.000** off capture, and the ladder prototype is **negative** (mwae 0.0413 → 0.0426, confidently-wrong +10.7 %). TRAPS: panel-before-src |
| **a LEVEL transfer from the intron** | toy + panel | ⛔ **REFUTED**, +0.207 on capture-ON × unstranded — capture inverts which side is well-counted (TRAPS: capture-inverts-the-counted-side) |
| **the RNA fragment-length model** | `length_ceiling.py`, one pmf at a time | ⛔ **−0.02 %** at pass-0, **+0.21 % (worse)** over all objects. Root cause exact (`pi(w)` scores junction *crossing*, the pool requires the splice to be *seen*). ⭐ Its value is the BOUND: the whole fragment-length-model cluster costs ≤0.43 % of the shipped solve. TRAPS: price-the-halves-separately |
| **TRAPS: pure-and-length-censored's κ residue, as an ACCURACY fix** | κ injected at exactly ½, all 36 conditions | ⛔ **−0.2 %** unstranded, worse on the shipped solve. ⭐ But the *general* defect — a boolean licence flipped by a small residue — is **the-capture-level-residual**, and the destruction control taught TRAPS: honesty-metrics-reward-ignorance |
| **a nascent-bearing ladder condition** | toy, 36 conditions × 7 rungs | ⚠ **−5 %**, and the wrong way on one stratum. Keep it as a harness arm (`--nrna 60`); it no longer justifies re-simulating the panel |
| **the gDNA prior's BIMODAL CAPACITY, and "give the prior more signal"** | a read of `gdna_landscape.py` + the production refit on real conditions | ⛔ **BOTH BRANCHES CLOSED.** The prior already renders the landscape correctly — **2.98 decades** of mode separation at `g75 ss0.99 capture_ON`, 30× more enriched mass ON than OFF, a single pile at the wall at `g00`. And a prior fitted from ORACLE truth is the same prior (0.04 dec). Not capacity, not signal, not location. §3 below |
| **§2's Jeffreys MEAN density location** | `--arm eta`, the `g00` zero control | ⛔ **REFUTED at +96,299 %.** It cannot say ZERO (`node_init.rho_g` is an exact 0 at 60,544/70,176 slots — the statement earning the −98 % at `g00`), and the TRAPS: a-ratio-cannot-carry-zero benefit it was credited with belongs to **§4**. ⭐ If revisited the derived form is the Gamma **MODE** `max(a−½,0)/E`, which is exactly 0 at a zero count |
| **a threshold anywhere in the licence family** | TRAPS: a-threshold-on-a-fitted-residue implemented and refuted one | ⛔ τ is continuous across the region, so any floor is a tuned constant (TRAPS: a-threshold-on-a-fitted-residue, TRAPS: a-licence-with-no-floor, TRAPS: a-multiplication-gated-by-a-trace — refused three times) |

### §4.1 ⛔⛔⛔ THE GRAVEYARD — ELEVEN MECHANISMS PRICED, ELEVEN REFUSED. DO NOT REBUILD THESE.

⭐ Promoted from a working doc when it was deleted (2026-08-07). Every row is a
real build that was measured and refused. ⛔ **`g00` is the owner-required ZERO-gDNA control: its truth is
exactly 0, so every fragment there is a false positive with nothing to cancel it** — which is why a
mechanism can look good on its target and still be inadmissible.

| candidate | what it did | `g00` | its target | why it died |
|---|---|---|---|---|
| `zc_jeffreys_mean` | `ρ_g = ½/E_g` at zero mass | ⛔ +7,269 % | −13.9 % | moves the mode UP |
| `zc_logmean` | `ρ_g = e^{ψ₀(½)}/E_g` | ⛔ +6,264 % | −11.3 % | moves the mode UP |
| `zc_anchor_mute` | no `prec_g` at empty locked slots | ⛔ +5,554 % | −7.7 % | kills the zero-gDNA win |
| `zc_struct_lock_g1` | scope `struct_lock` to `g1_locked ∧ NODE` | ⛔ +3,207 % | −1.2 % | ⭐ the MIS-SCOPED mask is load-bearing |
| `zc_reference_var` | `Var(f_g) = ⅛` where `τ = 0` | ✅ +0.0 % | −0.3 % | ⭐ passes the control and is INERT |
| `zc_discrepancy` | `+½ log D` shift, `(log D)²/12` | ⛔ +982 % | panel +4.5 % | moves the mode UP |
| `zc_disc_var` | the variance alone, mode untouched | ⛔ +255 % | panel +0.9 % | damping cannot bite |
| `zc_ref_prior` | own belief = ψ's reference, `τ + 1/π²` | ⛔ +3,792 % | −14.9 % | moves the mode UP |
| `zc_ref_prior_damp` | the two above, PAIRED | ⛔ +3,809 % | −15.5 % | ditto |
| the `eta` rebuild | a clean frame-free re-derivation | ⛔ unbounded | +85–103 % | see `DESIGN.md` §6.1 |
| §2's mean-location as a structural floor | the same idea in the LEVEL channel | ⛔ +96,299 % | — | it cannot say ZERO |

⭐⭐ **THE PATTERN, AND IT IS THE REAL RESULT: every one of the eleven was a rule for how to resolve DOUBT,
and at `g00` the doubt must resolve to NO gDNA.** A rule that lifts an evidence-free slot off zero is
inadmissible there however well it scores elsewhere. ⛔ The only candidate the control has ever ENDORSED is
the one-sided certified-RNA bound (−81.9 %, 8/8), and that one is panel-negative alone — it is half of
`TRAPS: a-cancelling-defect-pair` and is §1's **the-cancelling-pair**.

⚠ Four more `zc_*` arms exist as decomposition REVERTS used to attribute the 39 % win, not as proposals:
`zc_own_count`, `zc_live_count`, `zc_total_n` (inert) and `zc_transfer`, which reproduces the pre-fix tree.

---

---

## §5 THE OPEN QUESTIONS, ranked by how much they would change a decision

1. ⭐⭐⭐ ~~**WHY IS PRIOR FIDELITY ANTI-CORRELATED WITH DELIVERABLE QUALITY?**~~ ⭐⭐ **LIKELY ANSWERED,
   AND THE ANSWER IS THAT IT IS NOT THE PRIOR — IT IS THE MESSAGES** (2026-08-06, ψ-boundary ablation on
   HEAD at `g01 ss0.50 capture_on`, TRAPS: byte-identity-gate exact). At every one of HEAD's 12 worst slots the **self-solve WITH
   the fitted prior is nearly correct** (0.0043–0.0344 against truths of 0.0007–0.0136) and **the message
   layer then destroys it** (0.05–0.83). Muting the certified-RNA channel alone recovers the truth at 8 of
   the 12. So a better prior cannot show up in the deliverable: it is already right at the objects that
   carry the error, and something downstream overwrites it. ⛔ That also re-reads the earlier evidence
   rather than contradicting it — "the prior is degraded on unstranded where the solve improves" is what a
   prior looks like when it is not the binding constraint. ⚠ Confirm on a second stratum before closing.
   §1.1 above, TRAPS: the-intermediate-is-not-the-deliverable.
1b. ⭐⭐ ~~**Is the phantom-gDNA floor the same bug as the vertex problem?**~~ ⛔ **ANSWERED, and NO** — it was
   a different and much larger bug (§0's 39 %), and its residual is the CAPTURE LEVEL defect, not the vertex.
   TRAPS: a-zero-count-is-a-measurement/TRAPS: a-ratio-cannot-carry-zero/TRAPS: the-divergence-was-a-barrier.
2. ⭐⭐ **What is the joint price of the splice-flux reframe and the level defect?** **the-cancelling-pair**. Neither has an
   honest price alone, and one of them is already in `src/`.
3. ⭐ **Do we solve ALTERNATIVE SPLICING correctly?** The `alt_splice` rung exists and is **unverified** — see
   the handoff. Cheap, and it is the only structure where several junctions share an EDGE.
4. ⚠ **The two capture-ON pilot rows disagree about the SIGN of every length correction**, across two
   independently built panels. Unexplained. ⛔ Do not average them; find which one is lying.
5. ⚠ **Does the reframe's own `σ²_transfer` correctly price a ratio built on 19 counts?** `EQUATIONS.md`
   §3.5d says it is the right medicine for the noise and the wrong medicine for the share-weighting.
   Unmeasured.

---

## §6 THE RULES THAT MADE THESE NUMBERS TRUSTWORTHY

⭐ These are `TRAPS.md`'s operational core, repeated here because they are what a new session must do on day
one, not read about later.

- **Measure the CEILING before building the correction** (TRAPS: measure-the-ceiling-first). It has re-ranked this project
  **five** times, and it is also how you learn a phase is finished — §0 is that story for Stage A.
- **Run the PANEL arm before writing a mechanism into `src/`** (TRAPS: panel-before-src). Three toy-positive changes have been
  panel-negative.
- **Zero-gDNA and zero-RNA controls on every experiment** (owner, 2026-08-05). The truth is a constant, so
  every deviation is a false positive with nothing to cancel it — and **re-anchor-on-the-deliverable** was found this way.
  ⚠ But check the arm could have fired at all (TRAPS: could-the-arm-have-fired).
- **Quote `mwae` over ALL objects and the raw Σ|err|**, never the honesty metrics alone (TRAPS: honesty-metrics-reward-ignorance).
- **One thing varied**, with the baseline re-recorded from the current tree in the same session (TRAPS: re-record-the-baseline).
- **A falsification test first, verified failing — then break the fixed code and watch each gate fire**
  (TRAPS: perturb-every-gate). And name the observable for *each place* the change was made (TRAPS: name-the-observable-per-site).
