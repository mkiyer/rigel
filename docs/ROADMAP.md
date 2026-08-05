# ROADMAP — where the tool is, and what to do next

Reading order for a new session: **`SUCCESS.md`** (how performance is measured) → **this file** (current
state and priority) → **`TRAPS.md`** (mistakes not to repeat) → `DESIGN.md` / `EQUATIONS.md` as needed.

⛔ **Every number in this file was measured on the current tree and the current panel.** Re-derive rather
than trust: `SUCCESS.md` §"The instruments" runs the whole set in ~15 minutes. A number that has moved is a
result, not a documentation bug.

---

## §0 THE STATE, in six numbers

Measured on the 8-condition pilot (chr21 + chr22 + 92 ERCC, 10 M RNA fragments), post-drain, against the
simulator's own per-fragment truth.

| | scored against | |
|---|---|---|
| fragment length, the **anchor** | 4 zero-gDNA conditions, where anchor and pool are ONE population | ✅ **−0.00 % mean / ±0.01 % sd** |
| the **RNA** length model | the same | ✅ **−0.21 … +1.12 %** |
| the **gDNA** length model | 4 contaminated conditions | ✅ **−0.01 %** off capture · ⚠ **+7.9 %** under it |
| the second pass, **per fragment** | 339,765 held fragments with an unambiguous true length | ✅ **90.6 % exact**, mean error **−0.02 bp** |
| ⛔ **the deliverable** — library gDNA fraction | 4 contaminated conditions vs the simulator's origin counts | ⛔ mean \|error\| **0.0263** (all 8: **0.0245**) |
| ⛔ **the worst row** | `gdna100 ss0.50 capture_on` | ⛔ **0.4607 against 0.5** |

⭐ **And the PER-OBJECT answer, which is the one Stage B is about** — `pass0_vs_oracle.py`, 4 contaminated
conditions, mass-weighted `|Δf_g|` against the origin-split oracle. ⚠ **Undrained**, so not comparable to
the drained deliverable above; see `SUCCESS.md` Stage B for why it cannot be otherwise.

| | node axis | edge axis |
|---|---|---|
| ⛔ **pass-0** (prior-free) | **0.1192** | **0.1350** |
| the shipped final solve (3 refits) | 0.0586 | 0.0633 |
| ⭐⭐ what perfecting BOTH length models is worth | **−0.7 %** (pass-0) · −2.5 % (final) | +1.1 % · −1.5 % |

---

## §1 ⭐⭐ THE LENGTH PHASE IS DONE, AND THAT IS A MEASUREMENT

`SUCCESS.md` defines a length model as accurate enough **when handing calibration the exact distribution
changes nothing**. That is directly measurable — `calibration_truth_ab.py --ceiling` swaps the fitted pmf
for the simulator's own — and it has now collapsed:

| the FL pmf handed to `calibrate` is the only thing varied | before the four-pool model | ⭐ **now** |
|---|---|---|
| shipped | 0.0329 | **0.0263** |
| the EXACT **RNA** length distribution | 0.0325 (−1.1 %) | **0.0261 (−0.6 %)** |
| the EXACT **gDNA** length distribution | 0.0258 (−21.5 %) | **0.0258 (−1.7 %)** |
| ⭐⭐ **BOTH exact — the ceiling on the whole length phase** | 0.0256 (**−22.2 %**) | ⭐⭐ **0.0256 (−2.6 %)** |

⭐⭐ **Perfecting both fragment-length models would now buy 2.6 % of the deliverable, down from 22.2 %.**
The four-pool gDNA model harvested ~92 % of the headroom that existed, and **97 % of the remaining error is
not in the accumulator's length models.**

⚠ **The gDNA pool is still 7.9 % off under capture and it no longer matters much**, which is the whole
point of measuring a ceiling rather than a residual. The reason: what the solver consumes is the **gap**
`μ_g − μ_r`, and correct hybrid capture narrows the *true* gap to **−8.0 bp on a ~230 bp mean** — a 3.5 %
separation, where the 2×2 deconvolution is barely identified at all (`EQUATIONS.md` §3.1, `TRAPS.md` F3).
The fitted gap is now **+7.4 bp**: right magnitude, wrong sign, and worth 1.7 % because the channel is
nearly uninformative in this regime either way.

⛔ **So do NOT chase the +7.9 %.** Its cause is known and is not a divisor: both opportunity functions
assume gDNA is placed **uniformly** along the genome, and under capture it is not. Closing it needs a
capture-aware placement model, which is a large build for ~1.7 % of a number that is 97 % dominated by
something else.

### The four gDNA pools, each against its own opportunity

`scripts/design/gdna_pool_census.py`, on `gdna100 ss0.50`. ⭐ **The de-tilt is validated by four pools with
opposite raw tilts all landing within 0.3 % of truth off capture** — that is what says the opportunity
functions are right, and it is a much stronger statement than the combined number alone.

| pool | n (off / on) | raw vs truth | ⭐ de-tilted vs truth |
|---|---|---|---|
| `DNA_INTERGENIC` (contained) | 5.21 M / 237 k | −0.12 % / −14.9 % | **−0.00 %** / −14.8 % |
| `DNA_INTRONIC` (contained) | 4.01 M / 183 k | −1.06 % / −15.7 % | **−0.00 %** / −14.8 % |
| `DNA_INTRON_EXON` (crossing) | 206 k / 1.10 M | **+9.9 %** / +12.5 % | **−0.06 %** / +7.4 % |
| `DNA_INTERGENIC_EXON` (crossing) | 33.7 k / 199 k | **+18.1 %** / +28.7 % | **−0.28 %** / +20.1 % |
| shipped before — the contained pair | 9.22 M / 420 k | −0.53 % / −15.3 % | −0.00 % / −14.8 % |
| ⛔ the four **pooled raw** | 9.46 M / 1.72 M | −0.23 % / **+7.6 %** | — (there is no single divisor) |
| ⭐ the four, **each de-tilted** | 9.46 M / 1.72 M | — | ⭐ **−0.01 %** / **+7.9 %** |

⭐⭐ **Note the pool sizes: under capture the crossing pair holds 1.30 M fragments and the contained pair
420 k.** Off capture it is the reverse, 9.22 M to 240 k. That is the whole argument for fitting from four
pools rather than two — and it is why the contained pair alone cannot be rescued by a better divisor.

⚠ **The +7.9 % residual is entirely in the two crossing pools** (+7.4 % and +20.1 %), and it is the
non-uniform placement, not the divisor: their opportunity assumes gDNA is spread evenly along the line, and
under capture it is concentrated at the probe.

---

## §2 ⭐⭐ WHAT TO DO NEXT — Stage B, and its first step is DONE

`SUCCESS.md` gates Stage B behind "Stage A condition 2": the accumulator must stop being the limiting
factor. **At a 2.6 % ceiling it has stopped.** The priority therefore moved, for the first time in several
milestones, out of the accumulator and into calibration.

### ✅ What has landed — the pass-0 oracle comparison, and what it measured

`scripts/design/pass0_vs_oracle.py`, gated by `tests/calibration/test_pass0_vs_oracle.py` (nine gates,
nine perturbations, 9/9 fired). This closed the *previous* §2 steps 1 and 2 — build the comparison, and
localise the residual per object class. ⚠ **The ranked list further down is a NEW list**; its step 1 is
not the old step 1. On the four contaminated conditions:

| | |
|---|---|
| ⛔⛔ **THIS ROW WAS THE WRONG QUESTION AND ITS NUMBER MUST NOT BE QUOTED** — "the relay carries 93.1 % of pass-0's error" | It counted every object with mass, so **honest ignorance was scored as error**. Pass-0 is the PRIOR-FREE solve; an object with no own evidence reporting `f_g ≈ ½` at **zero precision** is *correct* — it is saying "I cannot be solved without a prior", which is true. Re-measured with the undetermined population excluded, pass-0's error on the objects it CAN solve is **0.0456**, not 0.3150, and 99.5 % of what the old number called error was the undetermined class. See `solvability_audit.py` |
| ⭐⭐ on an **unstranded** library the **density model carries 100 % of the own-evidence budget**, and it is wired to ONE slot population | strand contributes exactly 0 at κ = ½; the intron factory reaches **100 % of intron-node mass, both-stranded included**. It is applied to `(NODE, INTRON)` only, so the relay-only set is *exon nodes + edges*: 27.6 % + 20.6 % of mass off capture, 53.3 % + 44.6 % on it. ⛔ **100.0 % of that mass has a count AND a gDNA opportunity** — the information is present and unread |
| ⭐ **structurally-locked objects carry 0.0 % of the error**, every row | `DESIGN.md` §7's figure reproduced exactly. They are also 90 % EMPTY (12,466 nodes, 1,260 with mass) — a scorer that turned 0/0 into a number would report them as solved |
| ⭐⭐ **perfecting both length models is worth ≤ 2.5 % of the per-object error, and is NEGATIVE on 4 of 8 arms** | the per-object confirmation of §1. `C_input − P` is a rounding error; the length inputs are not the limiter at the object level either |
| ⚠ the refit loop **halves** the per-object error (0.1192 → 0.0586 node, 0.1350 → 0.0633 edge) | ⛔ and that **contradicts** the recorded "refit iterations 1→5 went 0.0840 → 0.1056, monotonically worse" (`TRAPS.md` D3). Different axis (0→3 vs 1→5) and a different metric, so it is not a refutation — but the D3 number must be re-derived before it is cited again |

⭐⭐ **THE RIGHT QUESTION, AND THE ONE THE OLD CALIBRATION MODULE ALREADY ASKED.** Pass-0's job is not
accuracy — it is to be a **substrate the gDNA hyperprior can be fitted against**. So the measurement is
a three-way partition, not an average (owner, 2026-08-03; the deleted `node_error_attribution.py` and
`confident_fp_trace.py` had this and it was lost in the refactor):

1. **UNDETERMINED** — no own evidence. ⛔ **Excluded from the error denominator.** Its only failure
   mode is the opposite one: claiming precision it has not earned.
2. **SOLVABLE and right.**
3. **SOLVABLE and wrong**, split by confidence. ⭐⭐ **Confidently wrong is the defect** — a wrong value
   with a tight variance outvotes correct neighbours and *anchors the prior*, so it propagates.

⭐⭐ **MEASURED ON THE gDNA LADDER** — 36 conditions, equal fragment lengths, gDNA 0 → 98 % at a fixed
10 M total, all 6 simulator gates passing (`scripts/sim/configs/gdna_ladder.yaml`). Reported **per
stratum, never pooled** (`TRAPS.md` B10): the strata behave so differently that any average over them
describes nothing.

Re-recorded 2026-08-04 with the corrected instrument. ⭐ `weak%` is the new column and it is the one to
read first: the share of the **scored** error sitting on objects whose own evidence has sd(λ) ≥ 10 nats.

| stratum | solv % | ⭐ weak % | mwae | calib ratio | relay hurts |
|---|---|---|---|---|---|
| ss0.50 capture_**off** | 2.0 – 96.2 % | **0.0 %** except ⛔ **79.1 %** at g75 | 0.0102 – 0.0852 | 0.58 – 3.78 | 5/8 |
| ss0.50 capture_**on** | **0.1 – 6.3 %** | 0.0 % | **0.0455 – 0.2825** | 0.62 – 0.86 | 5/8 |
| ss0.99 capture_off | 76.9 – 99.1 % | 0.0 – 0.2 % | 0.0055 – 0.0142 | ⛔ **3.81 – 8.90** | **8/8** |
| ss0.99 capture_on | 81.8 – 81.9 % | 0.2 – **16.2 %** (g98) | 0.0029 – 0.0231 | ⛔ 1.59 – 3.94 | **8/8** |

Panel aggregates, unchanged: **the relay hurts the solvable set on 26/32 contaminated conditions**
(mean Δ +0.0027); **the declared precision is not earned on 17/32**.

⭐⭐ **`weak%` ISOLATES THE τ ≈ 0 CONTAMINATION TO ESSENTIALLY ONE ROW — AND THAT ROW IS THE ONE THAT WAS
ALREADY UNEXPLAINED.** It is 0.0 % on 33 of 36 conditions and **79.1 % on
`gdna_g75_ss_0.50_nrna_none_capture_off`** — which is exactly the row the previous roadmap flagged as
breaking the pattern on four axes at once (223,898 confidently-wrong, mwae 0.0495 against 0.0133 and
0.0109 at the rungs either side, the only unstranded row with calib > 1, the largest negative relay
delta). **All four anomalies are one thing:** at that rung κ̂ landed 1.9σ from ½, τ crossed 1e-9 on 79 %
of the error mass, and a population of unsolvable exons was promoted into the scored set. ⭐ That closes
the "anomalous row" question outright.

⛔⛔ **AND THE `solv %` COLUMN OF THE UNSTRANDED ROWS IS AN ARTIFACT, AS IS THE "DIFFERENT REGIME"
CONCLUSION DRAWN FROM IT.** The old cut called a slot solvable at `tau_lam > 1e-9`. The strand
arm's information is `I(f_g) ∝ (2κ−1)²` — exactly zero only at κ = **½** — and κ is *fitted*:

| unstranded row | fitted κ | (2κ−1)² | strand τ on `edge/intron\|exon` | sd(λ) it earns | old class |
|---|---|---|---|---|---|
| `capture_off` | 0.500689 | 1.90e-6 | **5.3e-7** | **1,377 nats** | ⛔ "SOLVABLE" |
| `capture_on` | 0.499972 | 3.23e-9 | below 1e-9 | ∞ | "undetermined" |

**Same physics, opposite classification, decided by the 4th decimal of a nuisance parameter.** The
solver's λ grid spans ±10 nats, so an object declaring sd(λ) = 1,377 has no answer of its own — yet
0.00 % of that class fell below the grid span and 100 % of it was scored. ⭐ Nothing upstream is broken:
the derived deadband correctly *admits* κ̂, because κ̂−½ is a real **1.9σ** fluctuation at 10 M
fragments. The signal is statistically genuine and physically nil; only a **binary cut** cannot tell
those apart (`TRAPS.md` B11).

⭐⭐ **REPLACED BY A CURVE, and the curve says the thing directly.** `solvability_audit.py` now reports
`sd(λ) = 1/√τ` per decade — same device as its confidence curve, no cut. On
`g75 ss0.50 capture_off`, node axis:

| own-evidence sd(λ) | objects | what they ARE | mass | Σ\|err\| | err share | pred f_g | true f_g |
|---|---|---|---|---|---|---|---|
| CERTAIN (G1) | 1,252 | intergenic | 3,904,109 | 0 | 0.0 % | — | — |
| 1–10 nats | 8,942 | ⭐ **100 % introns**, factory 100 % | 2,994,792 | 82,473 | 20.9 % | 0.9725 | 1.0000 |
| ⛔ 10–100 | **56** | **100 % exons**, median **2,210 bp** | 343,550 | 125,708 | **31.9 %** | 0.4022 | 0.0363 |
| ⛔ 100–1000 | 2,716 | **100 % exons**, median 614 bp | 675,898 | 176,062 | **44.6 %** | 0.3985 | 0.2282 |
| >= 1000 | 6,476 | exons + edges, shallow | 46,475 | 10,113 | 2.6 % | 0.4315 | 0.4968 |

⭐⭐ **79.1 % of the "solvable" error sits on objects whose own evidence is wider than half the
solver's representable range.** The 20.9 % that is genuinely earned is *entirely* `node/intron` at
2.75 % low, and it is the only band with a **factory**.

⭐⭐ **AND THE BANDS ARE STRUCTURAL CLASSES IN DISGUISE — so "weak evidence" is the wrong reading.**
`τ = N·disc·[…]`, and on this row `disc = 1.36e-6` for *every* slot. So on exons τ ∝ **depth**, and the
sd(λ) ordering across the three exon bands is just exon size (2,210 bp → 614 bp → shallow). The 56
objects are simply **the 56 deepest exons**; their evidence is exactly as nil per fragment as the
2,716's. ⛔ There is no gradient of evidence quality here — there are **two** populations, the intron
factory and depth-scaled strand noise, and only the first is evidence.

⭐⭐ **What the 56 actually expose is a RELAY defect, and it INVERTS on a minority.** Every one has
`fg_loc` in 0.45–0.55 — the uninformative ½, i.e. no own answer — so the relay decides them entirely.
Split by what it did (exact, all 56):

| the relay… | n | mass | Σ\|err\| | share of the band | mean pred | mean true |
|---|---|---|---|---|---|---|
| ⛔ **pushed UP** | **12** | 83,424 | **69,460** | **55.3 %** | ⛔ **0.8713** | 0.0507 |
| left alone at ½ | 8 | 61,892 | 29,758 | 23.7 % | 0.5073 | 0.0187 |
| ⭐ pulled DOWN toward truth | **36** | 198,234 | 26,490 | 21.1 % | 0.1765 | 0.0353 |

⭐⭐ **So the relay is RIGHT on 36 of 56 and catastrophically inverted on 12** — and those 12 carry
55 % of the band's error. Tighter still: **10 objects with `pred ≥ 0.85` carry 63,001 fragments of
error, which is 16 % of pass-0's entire scored error on this condition.**

| node | bp | mass | fg_loc | pred | true | Σ\|err\| |
|---|---|---|---|---|---|---|
| 6,700 | 3,954 | 12,615 | 0.5098 | ⛔ **0.8888** | 0.0240 | 10,909 |
| 32,140 | 14,779 | 12,213 | 0.5098 | ⛔ **0.8888** | 0.0906 | 9,748 |
| 2,761 | 3,503 | 9,894 | 0.5489 | ⛔ **0.9162** | 0.0255 | 8,813 |
| 3,005 | 3,233 | 18,327 | 0.4511 | 0.4511 (untouched) | 0.0133 | 8,024 |
| 14,154 | 13,492 | 29,332 | 0.5489 | ⭐ 0.1277 (rescued) | 0.0365 | 2,673 |

**The deepest, most RNA-rich exons in the library are being told by their neighbours that they are
~87 % gDNA.** ⭐ That the *same* machinery rescues 36 of them means a fix has to explain the **sign**,
not damp the channel.

### ⛔⛔ AND THE DEFECT IS GENERAL — IT WAS 0.0 % REPORTED ON EVERY RUNG BUT ONE

The 56 were visible only because `g75 ss0.50 capture_off` accidentally scored them. The error does not
depend on being scored. Selecting **every** exon node with no own answer (`fg_loc ∈ [0.4, 0.6]`) and
≥ 1,000 fragments, without reference to any solvable/undetermined class:

| condition | κ̂ | n | mass | pushed UP | Σ\|err\| on those | mean pred | mean true | **in scored set** |
|---|---|---|---|---|---|---|---|---|
| `g25 ss0.50 off` | 0.500063 | 605 | 2,869,843 | 87 | ⛔ **395,251** | 0.8291 | 0.0089 | ⛔ **0.0 %** |
| `g50 ss0.50 off` | 0.499690 | 471 | 1,829,672 | 79 | ⛔ **272,599** | 0.8552 | 0.0264 | ⛔ **0.0 %** |
| `g75 ss0.50 off` | 0.500689 | 268 | 801,757 | 52 | 121,313 | 0.8476 | 0.0848 | 79.1 % |
| `g90 ss0.50 off` | 0.501823 | 104 | 222,352 | 23 | 31,205 | 0.8639 | 0.2550 | 0.0 % |
| `g98 ss0.50 off` | 0.501537 | 22 | 35,124 | 14 | 4,746 | 0.8207 | 0.6744 | 0.0 % |
| `g75 ss0.99 off` | 0.009767 | 3 | 5,047 | 1 | 28 | 0.6062 | 0.6209 | 100.0 % |

⭐⭐ **The defect is 3× LARGER at `g25` than on the row it was found on, and 0.0 % of it is reported.**
`g25 ss0.50 capture_off` publishes **mwae 0.0170**; its excluded class carries **1,056,019 fragments**
of error, of which **6,425 objects (2,028,011 fragments) sit a mean 0.414 away from ½ with 100 %
declaring a finite precision.** ⭐ It is stranded-specific in the right direction too: with strand live
the same selection finds **3** objects, not 605.

⛔⛔ **This is a hole in the INSTRUMENT as much as in the solver.** Excluding the undetermined class is
correct, but an exclusion is a promise to check that class another way — `SUCCESS.md` even named the
check ("claiming a precision it has not earned") and nothing implemented it (`TRAPS.md` B13).
`solvability_audit.py` now reports the class bucketed by `|f_pred − ½|` with its error and its
precision claim, which is where the numbers above came from. ⚠ Note the lowest band is *innocent*: 1,716
objects at a mean 0.021 off ½ carry 309,155 fragments of error because truth there is near 0 or 1 and ½
is the honest answer. The **0.30–0.50** band is the defect.

⚠ Zero-gDNA rows have an **empty** solvable set and cannot exercise this metric: with no gDNA,
`n_gdna_obs = 0` gates the strand seed off and the intron background is uninformative, so nothing has
own evidence. Honest, but they are false-positive checks only.

⚠ **A second, smaller correction in the same place:** `struct_lock`/`locked` was `~solvable & is_node`
in three instruments, so every structurally-locked **edge** was filed as honest ignorance rather than as
certain — 42,459 fragments of exactly-right answers on the row above. The definition now has one home
(`node_geometry.g1_locked`) and is deliberately **not** the same mask as `node_init`'s same-named,
node-only one (`TRAPS.md` E14).

⛔ **REFUTED, DO NOT RE-PROPOSE: "the relay damage and the overconfidence are one bug."** The idea was
that an overconfident node outvotes correct neighbours through the relay. The pooled correlation
between calibration ratio and relay delta is **+0.219**, which looks weakly supportive and is entirely
an artifact of the strand split; stratified it is **−0.839** (unstranded) and **−0.049** (stranded),
both non-positive. They are two independent defects. ⭐ And the sharpest test refutes it directly:
`g75_ss0.50_capture_off` is the one unstranded row that IS overconfident (3.78) and it has the largest
*negative* relay delta (−0.0153) — the relay rescues it.

### ⛔⛔ THE JUNCTION LEAK: MECHANISM CONFIRMED, CEILING MEASURED, AND IT IS **NOT** THE THING TO BUILD

The previous entry here reported `edge/intron|exon` as "the first badly broken link", 21–28 % low with
the junction flux isolated as the cause, and called for a derivation. **The mechanism is real and is
now pinned exactly. Its ceiling is ~0, and the framing was wrong in two ways.**

⛔ **The control was not a control.** `edge/intergenic|exon` predicting exactly 1.000 is **G1** — it is
never solved at all and keeps its pinned `{0,0,1}` init (solv 0.0 %). It cannot be wrong, so its
correctness said nothing (`TRAPS.md` B12). ⭐ **The control that works** is the same class split on the
variable, junction flux present vs absent:

| condition | flux = 0 | flux > 0 | inflation on the flux subset | flux's share of the class error |
|---|---|---|---|---|
| `g75 ss0.50 capture_off` | 0.770 | 0.699 | 11.5× | 26 % |
| `g75 ss0.99 capture_off` | 0.818 | 0.810 | 11.7× | **4 %** |
| `g75 ss0.50 capture_on` | 0.630 | 0.591 | 2.9× | 10 % |
| `g75 ss0.99 capture_on` | 0.920 | 0.918 | 3.0× | **2.5 %** |

(truth is 1.000 on every row). **74–97.5 % of the class's error is present where no junction attaches
at all**, and on the two strata where these objects have own evidence the flux explains 2.5–4 %.

⭐⭐ **THE MECHANISM, stated exactly — it is a COMPONENT-SET MISMATCH in the reframe, not "the flux is
read as RNA".** `node_total_density` puts the junction channel into `ρ_tot` at an EDGE, and a NODE slot
structurally cannot have one (an exon's mature molecules are inside its *contained* count, not in a
junction bank). So every `node → line` hop is scaled by

    r = [ρ_contig(line) + ρ_J(line)] / ρ_contig(node)  =  (1 + ρ_J/ρ_contig) · r_true

— measured median **11.5×** on the flux-bearing subset off capture, and the delivered gDNA over-claim
`φ = c_g·E_g/M` tracks it linearly across five decades (inflation 1.04 → φ 0.93; 50.9 → 23.3).
⭐ **But it is COMMON-MODE**, so it barely moves the composition: `λ = log(f_g/f_R)` is invariant to a
shared scale, and a **50×** inflation costs only 0.12 of `f_g` (0.768 → 0.651). That is the design
working as intended, and it is why the mechanism is large and the consequence is small.

⭐⭐ **THE CEILING (`scratchpad`, 4 conditions × 4 arms, one thing varied, in-process):**

| arm | what it varies | `edge/intron\|exon` | verdict |
|---|---|---|---|
| **A** baseline | — | 0.1865 / 0.0809 / 0.2677 / 0.2631 | |
| **B** matched-component reframe | the J term out of `ρ_tot` on both sides | **worse on 2 of 4**, best case **−6.8 %** | ⛔ not worth a derivation |
| **C** junction channel deleted | frame + graft + peel together | ⛔ **worse overall**; library f_gdna 0.774 → **0.841** (truth 0.75) on the unstranded row, `node/exon` 0.285 → **0.611** | ⭐ the channel is **load-bearing** |
| **D** `edge_spliced` added to `ρ_tot` | the other asymmetry (S5.a2) | ±0.0005, uniformly inert | not a lever |

⭐ **Arm C is the useful result: it confirms the owner's refusal to gate the message off, on
measurement rather than principle.** Removing the junction channel is much worse, so it is carrying
real information — the graft is what lets an exon's measured mature flux set its RNA level.
⚠ Arm C also improves `node/exon` by **25 %** under capture (0.0431 → 0.0317) while wrecking it off
capture — a capture-dependent **graft frame** problem, which is the M8/P1d term the solver already
flags ("under capture the graft OVER-states, median φ = 2.45"). That, not the intron|exon line, is
where the junction channel actually costs something.

⛔ **So: do not spend the derivation on the intron|exon seam.** ⚠ And the panel cannot license a
stronger statement either way — it is `nrna_none`, so `edge/intron|exon` truth is 1.0 *because there is
no nascent RNA to cross*, and `TRAPS.md` F9 only says *mature* RNA does not cross such a seam
contiguously. A nascent-bearing condition would move that truth off 1.0 and is the substrate this
question actually needs.

⭐ **The anomalous row, still worth dissecting after that:** `g75_ss0.50_capture_off` breaks the pattern on four
axes at once — **223,898** confidently-wrong fragments (3.5× the next worst), mwae 0.0495 against
0.0133 and 0.0109 at the rungs either side, the only unstranded row with calib > 1, and the largest
negative relay delta. A ladder that is otherwise smooth going non-monotone at exactly one rung usually
has one mechanism behind it.

⭐⭐ **TWO SOLVER DEFECTS, BOTH SYSTEMATIC:**

1. **The relay degrades the solvable set on 26/32 contaminated conditions — and on 8/8 of BOTH
   stranded strata.** The message-free local solve beats the final solve wherever strand is live. This
   is the class the retired tooling named `P1_OVERRULED`, and it matters out of proportion to its size:
   these are exactly the objects the gDNA hyperprior is fitted on, so the damage propagates. ⚠ On
   unstranded rows it is only 5/8 and sometimes *helps*, so the defect tracks **strand being live**,
   not gDNA level and not confidence.
2. ⛔ **The declared precision is NOT earned where the strand channel is live** — the calibration ratio
   (realised RMS log error ÷ claimed sd) is **3.81–8.90** stranded/off-capture and 1.59–3.94
   stranded/on-capture, against **0.58–0.86** on every unstranded row but one. The solver is up to
   **8.9× overconfident** exactly where strand speaks. ⚠ Off capture it *worsens* as gDNA rises
   (3.81 → 8.90); under capture it improves (3.94 → 1.59). ⚠ A single-condition read of this metric
   said "every decile ≤ 1, not overconfident" — that was one unstranded row and it does not generalise
   (`TRAPS.md` B10).

⛔ **On that 36.2 % the accumulator kept the information and pass-0 does not read it.** The channel that
would read it is the per-node fragment-length likelihood, which is gated OFF — and ⛔ **this is not
permission to turn it on.** C_info was computed with the simulator's OWN pmfs.

⭐ **But the reason it is gated off is narrower than "the fitted gap has the wrong sign".** Re-measured
post-drain this session, the fitted gap is **right off capture and flipped only under it**:

| post-drain | gDNA model | RNA model | fitted gap | true gap |
|---|---|---|---|---|
| capture **off** | −0.01 % | −0.21 % | **−16.7 bp** | −17.1 bp ✅ |
| capture **on** | **+7.85 %** | +1.11 % | **+7.4 bp** | −7.9 bp ⛔ |

⛔ **The sign flip is NOT a sign error anywhere in the code — it is the gDNA pool's +7.9 % capture bias
(≈ +18 bp) swamping a true gap of −8 bp.** The estimator assumes no ordering, and that is the point: the
bias is one-directional (it inflates gDNA relative to RNA), so on a library where gDNA is genuinely the
LONGER component the same bias would *exaggerate* the gap instead of flipping it. **The failure mode
depends on the true ordering, which is exactly what may never be assumed** (`TRAPS.md` F7).

⚠ **This re-prices the +7.9 %.** §4 rules it not worth chasing at 1.7 % of the library deliverable — a
judgement made against the library yardstick. It is now the single thing standing between the solver and
a channel that could address 36.2 % of pass-0's per-object error, and it is capture-specific.

### ⭐ THE NEW RANKING

⭐⭐ **The per-issue detail lives in `docs/calibration/ISSUES.md`** — measured cost, how sure each
diagnosis is, and the reasoning behind the order. This table stays the single ranked list.

⭐ **Re-ranked 2026-08-03 by measurement, not by argument.** The previous #0 (the junction leak) was
measured and dropped: mechanism confirmed, ceiling ~0, and the channel is load-bearing. This is the
third time the ceiling discipline has re-ranked the project (`TRAPS.md` B1).

⭐⭐ **Steps 1, 1= and 4b LANDED 2026-08-04 and are struck below**; the scale rule's numbers are in §2b and
the message licence's in §2c. What §2c exposed is the new step 1.

| | step | shape | why now |
|---|---|---|---|
| ✅ **1** | ~~**the gDNA scale rule**~~ **LANDED** — `EQUATIONS.md` §3.5, `TRAPS.md` D4c, gates in `tests/calibration/test_gdna_scale_rule.py`. §2b has the ladder numbers | — | — |
| ✅ **1=** | ~~**the mass pin** (C11)~~ + ~~**TSS/TES, plumbed and unconsumed** (C4, was 4b)~~ **BOTH LANDED** — one licence with two conjuncts. `EQUATIONS.md` §3.5b/§3.5c, `TRAPS.md` D4d/D4e/D4f + B14, gates in `test_relay_mass_pin.py` + `test_terminus_population_licence.py`. §2c has the ladder and toy numbers, and the CEILING that says C11 cost the panel nothing | — | — |
| **1** | ⭐⭐ **PROPAGATE `Var(κ̂)` INTO THE STRAND ARM'S PRECISION** (`ISSUES.md` C2 / C7, `TRAPS.md` B11 + D4f) | a derivation, then one term | ⭐ **What the message licence exposed, and it is now a SOLVER defect rather than a reporting one.** κ = ½ means *exactly* zero strand information, but κ̂ = 0.500689 leaves `I(f_g) ∝ (2κ−1)²` at 1e-6…1e-4 — and the composition licence is a **precision** test with no floor, so those precisions now grant the imputation licence at every evidence-free exon. Measured cost: on `nested_exons` the licensed pin still fires through the gene interior, 0.2264 where the pin-off arm reaches 0.0760. ⛔ **Not a floor** — B11 implemented and refuted that (τ is continuous across the region, so any floor is a tuned constant). The honest fix is that the strand arm's precision must carry κ̂'s own variance, which drives it to ~0 by derivation on an unstranded library. ⭐ It also pays twice: C7's "declared precision not earned, up to 8.9× overconfident" is the same term |
| ⭐⭐⭐ **1=** | **THE TWO FACES OF AN `intron\|exon` EDGE** (`ISSUES.md` C12, derivation `EQUATIONS.md` §3.6) | ⛔ a **RE-SOLVE** ceiling, then build | ⭐⭐ **Owner's target, derived and already VERIFIED against oracle truth — this is a build, not a question.** One EDGE carries one gDNA density and two composition statements: the INTRON face matches its unspliced population `{gDNA, nascent}`, the EXON face adds the measured junction to match `{gDNA, nascent, mature}`. That makes the EDGE *solvable* (two equations, two unknowns, the intron's composition prior-free from the factory) and the exon's composition **arithmetic on measured terms**. ⭐ Face (I)'s job is to carry the intron's **349** counts to the edge rather than trust the edge's **8**: measured at m=1 the exon goes 0.625 → ~0.499 against a truth of 0.458. ⭐ Face (II) closes it two ways — the junction (tight at high RNA, useless at low) and the exon's own mass identity (the reverse) — which FUSE by inverse variance. ⛔ Its substitution ceiling understates it (`TRAPS.md` B17); run a real arm |
| **2** | ⚠ **Add a 1–10 % gDNA rate to `pilot.yaml`** | one config line + a re-run | `suite_resolves.py` requirement (f). Unchanged in substance: every §1/§3 number was measured at 0 %/100 %. ⚠ The *ladder* already covers 0.01–0.98, so this now blocks only the Stage-A (length) numbers, not Stage B |
| **3** | ⭐ **A nascent-bearing ladder condition** | one config axis | `nrna_none` is why `edge/intron|exon` truth is exactly 1.0, so the panel structurally **cannot** distinguish "no RNA crosses this seam" from "no *mature* RNA crosses it" (`TRAPS.md` F9). Every conclusion about intron↔exon seams is conditional on that, including the one dropped above |
| **4** | ⭐ **The graft's frame under capture** (M8 / P1d) | design, and a ceiling first | Arm C says deleting the junction channel *improves* `node/exon` by 25 % under capture and wrecks it off capture. So the graft is right off capture and over-states under it — the solver's own note already measures median φ = 2.45. This is where the junction channel actually costs something |
| **4b** | ⭐ **`graft_premise_logvar` PER STRUCTURAL CLASS — the DEBT the terminus bit was waiting for** | one term, re-derived | ⭐ The TSS/TES bits now have a consumer (§2c), so the half of C4 that remains is the one it was documented for: `graft_premise_logvar` is "**A DEBT, not a model** … stands in for a quantity that splits ≥30× on whether the boundary carries a transcript TERMINUS — a bit the region map does not have". ⚠ It pairs with 4, and ⛔ its scope is **DONOR/ACCEPTOR**, which §2c deliberately left alone: a splice site changes the population too, but its flux is *measured* and the graft and the peel route it, so extending the licence there is a separate experiment against C3's contradictory arms |
| **5** | ⭐ **Make the fitted length gap's SIGN right, or prove it cannot be** | a measurement first | ⚠ **Demoted, and its supporting number needs re-deriving.** The "36.2 % identified-but-unread" figure was computed over a population defined by the `τ > 1e-9` cut, so it inherits B11 |
| **6** | ⭐ **The relay, panel-wide** | design, not measurement | ⚠ Step 1 is this defect's sharpest instance; do that first. ⚠ And the two numbers that used to rank this — "63.8 % on genuinely undetermined objects" and "the relay hurts 8/8 stranded" — both partition on the retired cut. Re-derive against the sd(λ) curve before designing. ⭐ Note the relay is not uniformly harmful even on the 56: it rescues node 14,154 from ½ to 0.128 (truth 0.037) while wrecking node 6,700 to 0.889 (truth 0.024). A fix must explain the SIGN, not damp the channel |
| **7** | ⚠ **The multimapper hole** | `SUCCESS.md` A3 | Unchanged in rank. Its cost is still **unmeasured**: measure it as a ceiling first |

⛔ **Steps 5 and 6 are the two that were previously #2 and #3, and they are demoted for a reason that is
not about their importance:** both were ranked on shares of an error total taken over the population the
`τ > 1e-9` cut defined, and that population was 79 % objects with no usable own evidence. The *defects*
are almost certainly real — the overconfidence and the relay damage were measured on the **stranded**
strata, where solvability is genuine — but the numbers that ranked them are not, and a ranking is a
number.

---

## §2b ✅ THE gDNA SCALE RULE — landed 2026-08-04, and what it cost

Derivation `EQUATIONS.md` §3.5 · lesson `TRAPS.md` D4c · gates
`tests/calibration/test_gdna_scale_rule.py` (6, perturbation-verified). **One line of arithmetic**: the
reframe is a composition imputation, so it needs the λ-emission gate's own **licence**, and where the
source cannot lend a composition the gDNA **level crosses unscaled**. No branch on capture anywhere.

### The five acceptance gates the owner set, each with its number

| | gate | result |
|---|---|---|
| 1 | the exon's `f_g` **TRACKS** truth across the density sweep instead of sitting flat | ⭐⭐ on `TA_single_exon` it was **flat at 0.896–0.910 across a 10,000×** RNA sweep; now **0.833 → 0.001** against a truth of **0.914 → 0.001**. \|err\| at 3× / 10× / 100× / 1000× RNA: **0.677 / 0.812 / 0.895 / 0.895 → 0.079 / 0.020 / 0.000 / 0.000** |
| 2 | every pure-gDNA object stays exact | ⭐ all four toy objects read **1.000**; on the panel `node/intergenic` mwae is **0.0000** on every condition, both arms |
| 3 | capture OFF and ON both work, **no branch** | ⭐ one expression; the unit gate is *parametrized* over a flat landscape and 20×/200× steps, and all 8 `ss0.50 × capture_on` ladder rows improve (−0.0007 … −0.0011) |
| 4 | the ladder's `weak%` curve **and** the undetermined-overreach table both improve | ⭐ `weak%` mean **3.18 → 2.91**, and on the row where it is large, **79.1 → 71.3**. Overreach Σ\|err\| at `g25 ss0.50 off` — the **1,056,019** fragments `TRAPS.md` B13 named as the panel's largest error at 0.0 % reported — is now **710,128 (−33 %)**, with the worst band (\|f−½\| 0.30–0.50) **450,429 → 212,909 (−53 %)**. Capture-ON: **2,029,961 → 1,433,967 (−29 %)**, worst band **−72 %**. Stranded rows are *inert*, correctly |
| 5 | ⛔ D4 not re-introduced — the delivered level independent of the destination's own total | ⭐⭐ on the unit fixture the delivered level is **exactly the laid-down field to 1e-9** and its spread over a 700,000× RNA sweep is **< 1.001×**, against a pre-fix **982.7×** |

### The panel, and the one stratum that got worse

| | baseline | fix |
|---|---|---|
| `node/exon` mwae, `g25 ss0.50 off` | 0.3002 | ⭐ **0.2019** (−33 %) |
| `node/exon` mwae, `g25 ss0.50 **on**` | 0.3850 | ⭐ **0.2686** (−30 %) |
| `node/exon` mwae, `g50 ss0.50 off` | 0.3084 | ⭐ **0.2150** (−30 %) |
| mwae, 32 contaminated rows | 0.0419 | **0.0415** — 15 better / 7 worse / 14 flat |
| confidently-wrong Σ\|err\|, mean | 23,759 | ⭐ **21,154 (−11 %)**; the largest row **223,898 → 114,036 (−49 %)** |
| relay Δ (the relay's own harm), mean | +0.0027 | ⭐ **+0.0023** |
| declared-precision ratio, mean | 2.802 | 2.787 |

⛔⛔ **AND THE HONEST COST: 7 of 32 rows are worse, and they are ALL `ss0.99 × capture_on`** — mwae
+0.0001 … +0.0021, confidently-wrong +2 … +15 %. ⭐ **The mechanism is known and it is the residual the
derivation declares**, not a surprise: the level an exon receives is the one measured at its flanking
EDGE, and a fragment spanning an EDGE at a gene end lies only PARTLY under the capture probe while one
contained in the exon lies wholly under it (measured 2.7× apart), so under capture the delivered level is a
**lower bound**. On an unstranded library the exon has no answer of its
own and a partly-low message beats a wildly-high one — every `ss0.50 capture_on` row improves. On a
**stranded** library the exon *does* have a correct own answer, and a low message now drags it. ⛔ Closing
it needs a probe-geometry model the tool does not have, and the affine extrapolation that would close it
is the simulator's own retention law (`TRAPS.md` D4c).

⭐ **A third result, unlooked-for:** the intron-mediated face of the same defect is fixed too. The toy
harness gate that *pinned* it (an exon inheriting its neighbouring intron's composition, 31× on the ladder
donor) **failed because the defect was gone** — 0.0107 vs 0.0112, a factor of 1.04 — and has been rewritten
from an ORDERING assertion to an INDEPENDENCE one, which is the stronger statement.

---

## §2c ✅ THE MESSAGE LICENCE — landed 2026-08-04, and it is where the licence became TWO conjuncts

Derivation `EQUATIONS.md` §3.5b/§3.5c · lessons `TRAPS.md` D4d/D4e/D4f + B14 · gates
`tests/calibration/test_relay_mass_pin.py` (7 + 1 strict xfail) and
`tests/calibration/test_terminus_population_licence.py` (11), both perturbation-verified with the matrix
recorded in each file's docstring. Two changes, measured separately against a baseline re-recorded from the
same tree in the same session.

**A — the MASS PIN is licensed** (this closes `ISSUES.md` C11). It fires iff no BELIEF can reach its
budget: the composition was supplied, **or** the destination is a structurally pure-gDNA object.
**B — the LICENCE gains the POPULATION conjunct** (this closes C4, and it is the first consumer of the
TSS/TES bits): `T(EDGE) = T(left) ∩ T(right)`, so a transcript terminus makes one flank's population
larger and the imputation across that pair is withheld.

### The ladder, all 36 conditions, `solvability_audit.py` — ⛔ `weak%` read before `mwae`

| | baseline | A: pin licensed | B: + population |
|---|---|---|---|
| solvable mwae, 32 contaminated rows | 0.0415 | 0.0414 | ⭐ **0.0413** |
| ⛔ `weak%` of the scored error | 0.0291 | 0.0291 | ⭐ **0.0280** |
| ⭐⭐ **confidently-wrong Σ\|err\|**, mean | 21,154 | 21,096 | ⭐⭐ **20,173 (−4.6 %)** |
| declared-precision ratio, mean | 2.788 | 2.786 | ⭐ **2.688** |
| relay Δ (the relay's own harm), mean | +0.0027 → +0.0023 | +0.0023 | ⭐ **+0.0021** |
| rows better / worse | — | 22 / 10, **worst +0.0000** | 20 / 12, worst +0.0008 |

⭐ Every aggregate moves the right way and `weak%` **falls**, so the gain is not error being moved onto
objects with no answer of their own. Best row `g75 ss0.50 off` **−0.0039**; worst `g98 ss0.99 on` +0.0008.

### ⭐⭐ The toy, where the mechanism is visible — `nested_exons`, mass-weighted `|Δf_g|` over the gene

| | + strand | − strand mirror |
|---|---|---|
| baseline | 0.2618 | 0.2633 |
| A: pin licensed | 0.2264 | 0.2145 |
| B: + population | ⭐⭐ **0.0541** | ⭐⭐ **0.0535** |
| *(the pin deleted outright — the CEILING run)* | *0.0760* | *0.0756* |

⛔⛔ **AND THE CEILING RUN IS THE RESULT THAT RE-RANKED THIS WORK, so read it before quoting the toy.**
Deleting the mass pin entirely moves the **ladder** by **+0.0002** (worse; better on 18/32, worse on 14,
worst row +0.0108) while making confidently-wrong **+13 %** and overconfidence **2.79 → 3.14**. So C11 had
**no measurable cost on the panel at all** — its −71 % on the toy lives almost entirely in the population
the audit deliberately excludes. ⭐ The licence is landed because it is *derived* and because it is free
(A regresses no row); it is **not** landed on the strength of the toy number.
⭐ And the ceiling run earned its keep twice: the pin-off arm **breaks the capture landscape** (the
off-probe floor leaks into every exon), which is what forced the structural second branch — `TRAPS.md` D4e.

### ⚠ What is still short, and it is now a DIFFERENT issue

The delivered level is exact on the synthetic uniform-field chain (1.000 at every slot, was up to 1.96×),
but on the real toy the pin still fires through the gene interior. Cause, measured: `lend` is a
**precision** test with no floor, and κ̂ = 0.500689 on a nominally unstranded library leaves the strand
arm's precision at 1e-6…1e-4 rather than 0 — so every evidence-free exon "supplies" a composition whose own
statement is 10²–10³ nats wide. ⛔ **That is `ISSUES.md` C2 / `TRAPS.md` B11 reaching the solver, not
C11**, and the repair is to propagate `Var(κ̂)`, not to add a floor (B11 already refuted the floor).
⭐ Two residuals are now separated and each has its own gate: the **level** is exact, and the RNA-free
interior exon still reads `f_g` **0.914** against 1.000 — ψ's uninformative reference at the vertex, pinned
as a strict xfail against the truth.

---

## §3 THE TWO ROWS THAT DISAGREE — the sharpest open question

The two capture-on conditions disagree about the *sign* of every length correction, and they have done so
across two independently built panels:

| condition | shipped | both length models exact | verdict |
|---|---|---|---|
| `gdna100 ss0.50 capture_off` | 0.0053 | 0.0053 | length-limited **+1 %** — byte-identical to 4 dp |
| `gdna100 ss0.99 capture_off` | 0.0364 | 0.0364 | **+0.1 %** — likewise byte-identical |
| `gdna100 ss0.50 capture_on` | 0.0393 | 0.0363 | **+7.6 %** |
| ⛔ `gdna100 ss0.99 capture_on` | 0.0240 | 0.0243 | ⛔ **−1.2 %** — a *perfect* length model makes it worse |

⭐ **Off capture the length channel is worth exactly nothing** — both rows are byte-identical across all
four ceiling arms. So whatever dominates those rows (0.0053 and 0.0364, a 7× spread between two conditions
differing only in strandedness) is **not** length and never was.

⛔ **And on `ss0.99 capture_on` a perfect length model is still slightly worse.** The effect shrank with the
four-pool model (it was −20 %) but it did not change sign, and the deleted panel showed the same sign on its
own worst row. **Something on the capture-on rows is partly cancelled by the length error.**

⭐ **The oracle has now answered half of that.** Scored PER OBJECT on the same row, an exact length model
is **better**, not worse: node `Σ|err|` −5.9 %, edge +0.5 % at the shipped depth. So the library-level
−1.2 % is a *cancellation* artefact, not a real degradation — the exact model moves per-object errors that
happened to be offsetting each other in the net. ⚠ That does not make the row solved; it relocates the
question from "why does truth hurt?" to "why do these objects' errors cancel?", and `B4`'s 274× says the
net was never the right place to ask it.

⚠ **The unstranded zero-gDNA row is the largest false positive by 3×**: `f_gdna = 0.0646` against a truth
of 0, where the stranded rows read 0.0021–0.0026. With no strand information composition comes from length
alone, and at a 3.5 % length separation there is nothing to come from.

---

## §4 WHAT IS DELIBERATELY NOT NEXT

| | why not |
|---|---|
| **chasing the gDNA pool's +7.9 %** | §1. Measured worth: 1.7 %. Needs a capture-aware placement model |
| **turning on the per-node fragment-length likelihood** | It is a *comparator* — it reads `μ_g − μ_r` as composition. Under correct capture the true gap is 3.5 % of the mean, so the channel is nearly uninformative, and the fitted gap currently has the **wrong sign**. ⛔ Turning it on would feed it a manufactured discriminant on exactly the conditions it exists to rescue. Its A/B harness exists (`length_likelihood_ab.py`) |
| **iterating the second pass** | Answered by measurement: the anchor is exact, the assignment is right to −0.02 bp, and the pool's residual is a **selection** effect — iteration cannot fix one |
| **regenerating the goldens** | ⛔ Not until the accumulator work is finished, then **twice, and diff**. They run under the default sampling mode, so a flaky expectation baked in now is permanent — and regenerating is not validating |
| **the traffic term's Poisson likelihood, and D-4** | `DESIGN.md` §4.2. Both deferred on a 15-fragment measurement taken on the deleted panel; ⚠ that measurement has **not** been re-derived, so re-measure before citing it as the reason |
| **overdispersion in the simulator** | `suite_resolves.py` requirement (c). Needs a mechanism decision *and* replicate conditions; nothing dispersion-dependent can be validated until then |
| **the paralog EM degeneracy** | `TRAPS.md` D9. A real unidentifiability, and an EM design call rather than accumulator work. It is the 22nd test failure and it is expected |
| **`POOL_EB_PRIOR_ESS = 1000`** | A magic number, and changing it is the owner's call. ⚠ Its previously-recorded effect is void: the pool it was said to dominate now holds 419,799 fragments under capture and 9.2 M off it. Re-measure before treating it as a live defect |

---

## §5 THE RULES THAT MADE THESE NUMBERS TRUSTWORTHY

Not process decoration — each caught something real, and each is expanded in `TRAPS.md`.

1. **A falsification test written first, verified failing, then the fixed code broken to watch each gate
   fire.** The four-pool model's 8 perturbations found one blind gate that review had not.
2. **No magic numbers.** The four-pool combination's weights are inverse-variance and its divisors are
   enumerated against an oracle; nothing in it was chosen.
3. **Score against truth, not against the previous run.**
4. ⭐⭐ **Measure the ceiling before building the correction** — and again afterwards, because that is how
   you learn a phase is over rather than assuming it is not.
5. **Prove the substrate before the code.**
6. **Real data is a TEST input, never a DESIGN input.**
