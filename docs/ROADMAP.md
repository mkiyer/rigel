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
| ⛔ **calibration, unstranded × capture-ON** | ⛔ **BLIND** — reports **0.033–0.058** while truth spans **0.00 → 0.98**, and hands the EM a gDNA prior **94.4 %** short | not noisy; a flat line. Still the whole open problem. ⛔ `length_likelihood` does NOT fix it — see the row below |
| ⛔⛔ **`length_likelihood`** | ⛔ **MEASURED 2026-08-10 AND IT IS INADMISSIBLE.** On the g00 ZERO-gDNA control, truth exactly 0, it reports **0.539 / 0.573** gDNA (unstranded). Mean `|err|` over the five ladder conditions tested **0.228 → 0.393 (+72 %)** | ⚠ It looks like an **87 %** win on the flgap pair — but that panel is ALL g50 and the channel returns a near-CONSTANT ~0.5 (g00 0.539, g50 0.522, g98 0.287 unstranded), so g50's truth of ~0.5 flatters it. §1.4 |
| ✅ **the prior ASSEMBLER and its POPULATION** | ✅ `rel` **0.0019–0.0027** with perfect masses in and **4.9e-4** with perfect per-component shares; the composition claim `a_g:a_r` is exact against the unspliced pool (`Δphi` **≤ 5e-4**) | ⭐ the assembler was **0.179**: the conserved-count rewrite took it to 0.0202 and the YARDSTICK took it to 0.0027 (`TRAPS: score-the-consumers-own-count`). ⛔ Two entries here once claimed a 72 % residual and a +0.07 tilt; both were the reference, not the tool |
| ✅ **the fragment-length models** | ✅ accurate — `pi(w)`'s de-tilt reads **211.77** against a true **212.20**, and the gDNA pmf is exact to **0.02 bp** off capture | ⛔ a claimed +10.7 % bias was a truth-parser bug (`TRAPS: a-truth-table-of-aggregates`). ⚠ §4 retires the RNA-FL row for the DIVISOR consumer only; a pricing is per consumer |
| ⚠ **the gDNA pmf under capture** | ⛔ the last unfixed length defect: **+13.6 bp** at a 330 bp gDNA mean and **+3.5 bp** at 120, drained per-bin `fit/true` **1.22 … 4.18** in the tail | ⭐ EXACT off capture (+0.1 / +0.0 bp) and untouched by the drain (Δ ≤ 0.1 bp), so it is a PLACEMENT problem — `EQUATIONS.md` §4.4 — and must not be attacked by editing `gdna_opportunity` |
| **message propagation** | ⛔ **OFF**, and it stays off until the tool is optimised end to end across all scenarios (owner, 2026-08-10) | net better on 3 of 4 strata (−58 / −44 / −32 %); +155 % on the fourth. ⭐ That fourth stratum is the one `length_likelihood` just fixed, so the relay's price there must be RE-priced, never inherited |
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

## §1 ⭐⭐⭐ THE VALIDATION CAMPAIGN — items 1–3 DONE, item 4 open

⛔⛔ **OWNER RULING 2026-08-07, STILL IN FORCE: no new features until the tool is characterised end to
end**, with message propagation OFF and (until 2026-08-10) `length_likelihood` OFF, so every number is
attributable to the tool as it stands rather than to a feature landing mid-campaign.

| | the item | the question | state |
|---|---|---|---|
| **1** | calibration-prior-vs-oracle | calibration ships a PRIOR — how wrong is it? | ✅ `prior_vs_oracle.py` (+23 gates) |
| **2** | tool-absolute-accuracy | transcript-level accuracy, end to end? | ✅ `quant_accuracy.py --arm base` |
| **3** | error-downstream-of-calibration | inject the oracle prior — whose error is it? | ✅ `--arm oracle` |
| ⏸ **4** | performance | it has been a while, and new code means slowdown | **NOT STARTED** — §1.0 |

⭐⭐ **THE ANSWER TO 1–3, AND THE QUESTION HAD A HIDDEN ASSUMPTION.** "Calibration or the EM?" presumed
ONE deliverable. There are two and they answer differently:

* **the library gDNA/RNA separation** — calibration is the whole bottleneck. A perfect prior takes mean
  `|f_gdna − truth|` from **0.1060 to 0.0097**.
* **transcript assignment** — it is NOT calibration. **31.1 %** of fragments are misassigned and a
  perfect prior removes only **32 %** of that; on the two capture-OFF strata a perfect prior is **3–4 %
  WORSE**. 95.5 % of the stranded × capture-OFF error is ordinary isoform ambiguity (6.94 M → 311 k when
  isoforms are collapsed), which is not a Rigel-specific defect.

⛔ **Two things a perfect prior cannot touch, and they are §0's open rows**: the nascent channel (20.2 M
fragments on zero-truth entities, and at `g00` a perfect prior makes it WORSE, 4.30 M → 4.64 M), and
`g00`'s 7.05 M gene-level misassignment with no gDNA present at all.

⛔ **THE `g00` ZERO-gDNA CONTROL FAILS AT THE PRIOR, AND ITS WORST STRATUM IS NOT THE BLIND ONE.** The
shipped prior claims **2,067,637 gDNA fragments in libraries containing none**, and **1,707,321** of
them are at unstranded × capture-**OFF** — a stratum that reads healthy (`rel` 0.046) on every
contaminated row, and worse by **5.5×** than the blind one. Only a zero control could have found it.

⭐ **Where the live numbers are.** They are not repeated here: `prior_vs_oracle.py` regenerates the
`P−O` / `O−Fo` / `S−Fo` / `Fo−F` tables and the composition and strength diagnostics on demand
(~48 s/condition), and §0 carries the handful that rank the work. A closed item earns a line, not a
report.

### §1.4 ⛔⛔ `length_likelihood` — MEASURED FOR THE FIRST TIME 2026-08-10, AND IT IS INADMISSIBLE

`scripts/design/length_likelihood_ab.py`, one thing varied (the config flag) on the same cached scan,
scored against each condition's TRUE `f_gdna` from the origin-split caches.

**On the flgap pair (8 conditions) it looks like a triumph:** better on 7 of 8, mean
`|f_gdna − truth|` **0.133 → 0.0175 (−87 %)**, and the blind stratum apparently resolved — 0.0324 → 0.5222
against a truth of 0.507. Mass with no own composition evidence falls **97–99 % → 5–7 %**, which is
exactly the hole the channel was built to fill.

⛔⛔ **AND THE ZERO-gDNA CONTROL DESTROYS IT.** Truth is exactly 0 at `g00`:

| ladder condition (truth 0.000) | OFF | **ON** |
|---|---|---|
| g00 ss_0.50 capture OFF | 0.1776 | ⛔ **0.5394** |
| g00 ss_0.50 capture ON | 0.0350 | ⛔ **0.5729** |
| g00 ss_0.99 capture OFF | 0.0030 | 0.1059 |
| g00 ss_0.99 capture ON | 0.0022 | 0.0539 |
| g98 ss_0.50 capture ON (truth 0.980) | 0.0576 | 0.2865 |

Mean `|err|` over those five: **0.228 → 0.393 (+72 %)**. A library with no gDNA is reported as 54–57 %
gDNA.

⭐⭐ **THE MECHANISM, AND IT EXPLAINS BOTH TABLES AT ONCE: on unstranded data the channel returns a
near-CONSTANT answer of about 0.5 rather than a truth-tracking one** — 0.539 at g00, 0.522 at g50, 0.287
at g98, while truth spans 0.00 → 0.98. **The whole flgap panel is g50, where that constant is right.**
On stranded data the effect is much smaller (0.054–0.106 at g00) because the strand channel still
supplies real evidence and dilutes it.

⭐ It is consistent with the channel census: at `g00` the two fitted pmfs are **1.2 bp** apart (there is no
gDNA to fit one on, so it falls back to the anchor), the channel's rows are near-flat, and its argmax
carries bias **+0.66** with median `|Δ| = 1.0000` at a precision pinned to `1/Var(λ grid) = 0.029`
(`EQUATIONS.md` §3d). Enabling it lets that speak at **100 %** of slots.

⛔ **So the smooth shrinkage (§3d) is NECESSARY AND NOT SUFFICIENT**: it corrects the PRECISION, and this
failure is in the MODE — the near-flat row's argmax, added to ψ. Fixing the precision alone would not
have prevented it.

⚠ **THE LESSON FOR THE PANEL, NOT JUST THE CHANNEL: an 8-condition panel at ONE gDNA level cannot
validate a composition estimator.** The flgap pair varies fragment length and capture and strand — and
holds `g50` fixed. `TRAPS: a-single-level-panel-cannot-see-a-constant`.

### §1.0 ⚠ ITEM 4 — AND THE PERFORMANCE SUBSTRATE IS A TRAP

Measured 2026-08-07 on one 10 M-fragment ladder condition (a 35,135-node chr22 index): per-locus EM
**15.9 s (47 %)**, native scan 6.5 s, calibration 6.5 s, second-pass drain 3.5 s, total 33.5 s.
⚠ Message propagation costs nothing measurable (33.7 s ON vs 33.5 s OFF) — the "the relay is the
expensive part" hypothesis is REFUTED.

⛔⛔ **BUT THAT PROFILE IS UPSIDE DOWN FROM THE REAL ONE AND BOTH ARE CORRECT.** Calibration is
depth-INDEPENDENT — every node in the index is solved regardless of read depth — so it scales with the
INDEX while the EM scales with the DATA. On real cfRNA at genome scale (~1.5 M nodes) calibration has
measured ~66 s against the EM's ~24 s: the exact reverse. **Profile on real cfRNA**
(`~/Downloads/rigel_runs/cfrna/`, genome index at `~/Downloads/rigel_runs/refs/rigel_index`), never on
this panel — `TRAPS: toys-rank-hotspots-backwards`, which cost a whole analysis once.

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
