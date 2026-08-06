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

**Stage A — the accumulator — is DONE, and that is a measurement.** `SUCCESS.md` defines a length model as
accurate enough when handing calibration the exact distribution changes nothing, and
`calibration_truth_ab.py --ceiling` measures exactly that:

| | |
|---|---|
| perfecting **BOTH** fragment-length models is worth | ⭐ **2.6 %** of the deliverable, down from 22.2 % before the four-pool model |
| fragment length, the **anchor** | ✅ −0.00 % mean / ±0.01 % sd |
| the **RNA** / **gDNA** length models | ✅ −0.21…+1.12 % · ✅ −0.01 % off capture, ⚠ +7.9 % under it |
| the second pass, per fragment | ✅ 90.6 % exact, mean error −0.02 bp |

⛔ **Do NOT chase the gDNA model's +7.9 % under capture.** Its cause is known and is not a divisor: both
opportunity functions assume gDNA is placed uniformly and capture does not. It needs a capture-aware
placement model, for ~1.7 % of a number that is 97 % dominated by something else. And what the solver
consumes is the *gap* `μ_g − μ_r`, which correct capture narrows to −8.0 bp on a ~230 bp mean — a regime
where the 2×2 is barely identified at all (`EQUATIONS.md` §3.1, `TRAPS.md` F3).

**Stage B — calibration — is where the error is.** Per-object, mass-weighted `|Δf_g|` against the
origin-split oracle, ⚠ undrained:

| | node axis | edge axis |
|---|---|---|
| ⛔ **pass-0** (prior-free — the substrate the gDNA hyperprior is fitted against) | **0.1192** | **0.1350** |
| the shipped final solve (3 refits) | 0.0586 | 0.0633 |
| the deliverable, library gDNA fraction, 4 contaminated pilot conditions | ⛔ mean \|error\| **0.0263** | — |

---

## §1 ⭐⭐⭐ WHAT TO DO NEXT

> ⛔⛔ **RE-RANKED 2026-08-05, and two entries below are STALE — read this box first.**
>
> **Step 1 (the vertex) is CLOSED as a build.** `vertex_ceiling.py` prices it and its own docstring leads
> with the verdict: the 24.4 % is **the value of missing information, not headroom** — evidence-free
> objects sit at median `|Δ|/sd(f_g)` = **z = 0.5–0.6**, inside their own 1σ, and no prior-free solve can
> do better. `test_certified_rna_licence.py` `C2b` independently closed the one channel that might have
> supplied vertex evidence. §2 below is kept as the *measurement*; it is not a build.
>
> **Step 3 (`Var(κ̂)`) is REPRICED — it is a REPORTING fix, not an accuracy fix, and it is now urgent for
> a different reason.** Measured 2026-08-05: the deadband opens on **3 of 18** unstranded conditions, and
> at `g75 ss0.50 capture_off` forcing κ = ½ takes `solv%` **77–90 % → 0.0 %** while Σ|err| moves
> −0.6 %/−5.3 %/**+4.6 %** across the three big classes (a wash, reproducing §3's −0.2 %). ⛔ **So the
> panel's `solv%` / `weak%` / `conf-wrong` / `calib` columns are inflated by the defect at those three
> conditions, and any ranking that reads them is reading the bug** (`TRAPS.md` A12b). Fix it to make the
> instrument trustworthy, not to move the error.
>
> ⭐ **What the same dissection says the ERROR actually is** — the classes, which are ablation-invariant,
> at `g75 ss0.50 capture_off`:
>
> | class | objects | mass | true `f_g` | pred `f_g` | Σ\|err\| |
> |---|---|---|---|---|---|
> | `edge/exon\|exon` | 12,811 | 1,742,598 | **0.1177** | **0.2291** | **320,297** |
> | `node/exon` | 11,304 | 1,383,829 | 0.1725 | 0.2683 | 242,309 |
> | `edge/intron\|exon` | 19,610 | 319,314 | 1.0000 | 0.7462 | 81,032 |
>
> **88 % of the error is the two exon classes, both biased the same way — gDNA over-called ~2×** — because
> on an unstranded off-capture library an exon object has **no working evidence channel at all**: strand is
> structurally silent, `length_likelihood` is OFF by ruling, and the intron factory is intron-NODE-only. It
> falls back to ψ's uninformative reference (~0.5) against a true `f_g` of 0.12–0.17. ⭐ **That, not the
> vertex, is the modal defect on the contaminated panel.**

| | step | shape | why now |
|---|---|---|---|
| ⛔ ~~**1**~~ | ~~**ψ CANNOT REACH A SIMPLEX VERTEX**~~ | **CLOSED — see the box above** | §2 below. The zero-RNA control reads **0.92–0.99 against a truth of exactly 1.000** on every rung, 0.65–2.51 % of mass. ⚠ Kept because the measurement is sound and re-deriving it is expensive; ⛔ the *conclusion* that it is a build is withdrawn |
| ⭐⭐ **2** | **PRICE THE SPLICE-FLUX REFRAME JOINTLY WITH THE LEVEL DEFECT** (`EQUATIONS.md` §3.6c + §4) | nothing to design: one is in `src/`, the other is `toy_trace_error.py`'s relay-level arm | The reframe is correct by derivation and now DEFAULT (owner's ruling, 2026-08-05), but alone it is panel-negative on the node axis. `TRAPS.md` D4j says that is not a verdict on it: it is half of a **cancelling pair**. Neither defect has an honest price until the joint arm runs |
| ⭐⭐ **3** | **PROPAGATE `Var(κ̂)` INTO THE STRAND ARM'S PRECISION** (`ROADMAP.md` §1 step 3) | a derivation, then one term | κ = ½ means *exactly* zero strand information, but κ̂ = 0.500689 leaves `I(f_g) ∝ (2κ−1)²` at 1e-6…1e-4 — and a boolean licence on a precision is flipped as hard by 1e-6 as by 1e+6 (`TRAPS.md` D4f). ⛔ **Not a floor** — B11 implemented and refuted that. Pays for C7's "declared precision not earned, up to 8.9× overconfident" at the same time |
| **4** | **The graft's frame under capture** (`ROADMAP.md` §1 step 4) + **`graft_premise_logvar` per structural class** (C4b) | design, and a ceiling first | Deleting the junction channel *improves* `node/exon` 25 % under capture and wrecks it off capture, so the graft is right off capture and over-states under it (median φ = 2.45). C4b is the DEBT the terminus bit was waiting for; ⛔ its scope is DONOR/ACCEPTOR, which the terminus licence deliberately left alone |
| ⭐ **4=** | **THE `intergenic\|exon` EDGE ↔ intergenic NODE PAIR AS A gDNA REFERENCE** | small build, ceiling first | ⭐ Owner: *"we must use this as a measure of gDNA density. The intergenic NODE ↔ intergenic-exon EDGE can be solved. They can be assumed to be pure gDNA and have the same composition."* — so it is a **solvable pair**, not an anchor. Mechanism identified: `struct_lock` is NODE-only, so a G1 EDGE gets `Var = ∞` and can only RELAY a level, never ORIGINATE one. ⚠ Under capture that matters: the EDGE carries **8.6×** more mass per object (139.7 vs 16.2) while the intergenic NODEs it depends on are depleted ~24×, and on the capture-ON toy the EDGE→exon message is **dead** (`cm_g = 0`) at low RNA. ⭐ On the toy the EDGE has 23 counts and the NODE beside it 60 — neither alone is the answer, the inverse-variance fuse of the two is. ⛔ `strand_evidence`'s docstring argues for the NODE-only scope and that argument must be re-examined for a G1 EDGE specifically: RNA cannot cross a gene boundary, so fragments spanning it are **not** RNA-contaminated |
| **5** | **The relay, panel-wide** | design, not measurement | Steps 1–2 are this defect's sharpest instances; do them first. ⚠ The two numbers that used to rank this both partition on the retired `τ > 1e-9` cut (`TRAPS.md` B11) — re-derive before designing. ⭐ A fix must explain the SIGN: the relay rescues some objects and wrecks others |
| **6** | ⚠ **Add a 1–10 % gDNA rate to `pilot.yaml`** | one config line + a re-run | `suite_resolves.py` requirement (f); every §0 length number was measured at 0 %/100 %. The *ladder* already covers 0.01–0.98, so this blocks only the Stage-A numbers |
| **7** | ⚠ **The multimapper hole** | `SUCCESS.md` A3 | Cost still **unmeasured**; measure it as a ceiling first |

---

## §2 ⛔⛔⛔ STEP 1 IN FULL — THE VERTEX PROBLEM

⭐ **Why it ranks first, in one sentence:** the annotation has >50,000 genes and perhaps ~10,000 are
expressed in any one sample, so **most annotated transcripts are OFF** — their objects are pure gDNA, they
are the majority of the genome's objects, and there is nothing to deconvolve, yet the solver is 5–8 % wrong
on them. Instrument: `scripts/design/zero_controls.py`.

| spec | G1 objects (ψ BYPASSED) | gene-body objects (ψ SOLVES) | worst Δ | error share of mass |
|---|---|---|---|---|
| `silent` | ✅ **1.0000 exact** | ⛔ 0.946 – 0.989 | −0.0540 | 0.65 % |
| `TA_single_exon` | ✅ **1.0000 exact** | ⛔ 0.922 | −0.0779 | 2.51 % |
| `spliced_exons` | ✅ **1.0000 exact** | ⛔ 0.928 – 0.984 | −0.0725 | 2.01 % |

**It is not the message layer, and three independent lines say so:**

1. **What the messages deliver is right.** The relay hands each gene-body slot a gDNA density whose implied
   composition `cg·E_g/M` is **0.886 – 1.148** — `f_g ≈ 1` to within the Poisson noise of the 14–21 counts at
   the source EDGEs — and on 4 of 10 objects it already implies **≥ 1**. ψ then returns 0.93–0.96.
2. **It does not shrink with counts.** Silent-exon length swept 1 kb → 30 kb at fixed gDNA density (57 →
   2,266 counts, 40×): Δ = −0.084, −0.078, −0.063, **−0.111**. Not a sampling error.
3. **It moves with the SOLVE's own parameterisation.** `sweep_logodds_window` 10 → 20 → 40 gives Δ = −0.073,
   −0.093, **−0.115** — ⛔ **widening the grid makes it WORSE**, which is `TRAPS.md` D5 (a grid supplies a
   prior whether you ask for one or not, and a wider window puts more mass away from the vertex).
   `sweep_n_grid` converges by 120, so it is not a resolution problem.

⭐ **The contrast that settles it:** the G1 objects are structurally PINNED, so ψ is bypassed — and they are
exact at 1.0000, while their own incoming messages imply 0.77–1.33. **Where ψ solves, the answer lands 5–8 %
short of the vertex; where ψ is bypassed, it is exact.**

⚠ This is the residual already pinned as a strict xfail (`ROADMAP.md` §2, "ψ's uninformative reference at
the vertex"). What is new is that it is measured on the dominant population, localised away from the
messages, and shown insensitive to depth and sensitive to the grid.

### ⭐⭐ Two things this already decides for you

* **A perfect per-component `r_g` (`EQUATIONS.md` §3.5d) is worth ZERO on the unexpressed population.** On a
  silent gene the true gDNA density is uniform, so the true `r_g` is exactly 1.0 on every hop — and that is
  already what the solver uses, because a G1 source has no RNA precision, `lend` is False, and the level
  crosses unscaled. This arm therefore **already has a perfect `r_g`** and still lands 5–8 % short. ⛔ Do not
  build the `r_g` fuse before the vertex is settled.
* ⚠ **The phantom-gDNA floor may be the SAME bug at the other vertex, and that is a HYPOTHESIS not a
  finding.** The zero-gDNA arm reads `f_g` = 0.0218 against exactly 0.0000; the zero-RNA arm reads 0.946
  against 1.000. Both ~2–7 %, both at a vertex. ⭐ Establishing whether one fix covers both is step 1's first
  measurement, because it doubles or halves the value of everything after it.

### The intron row, a second and smaller thing in the same table

The intron has real own evidence (`tau` 0.22–0.25) and its factory gets it to within 1 % (`fg_loc`
0.986–0.991). ⛔ The messages then make it **worse** — gap closed −17 %, delivered message implies
`f_g` = 0.68. Small (0.011–0.016), but it is the message layer degrading an object that already had the
right answer, on the simplest scenario that exists.

---

## §3 ⛔ WHAT IS DELIBERATELY NOT NEXT — one line each, so it is not rebuilt

| | closed by | verdict |
|---|---|---|
| **the gDNA scale rule** · **the mass pin** · **TSS/TES as the population licence** | landed 2026-08-04 | ✅ `EQUATIONS.md` §3.5/§3.5b/§3.5c, gates in `test_gdna_scale_rule.py`, `test_relay_mass_pin.py`, `test_terminus_population_licence.py`. ⚠ The ceiling says the mass pin cost the panel **nothing** (+0.0002 to delete it outright); it landed on the derivation and on being free |
| **face (I) of the `intron\|exon` EDGE** | re-solve ceiling + panel arm | ⛔ **DO NOT BUILD.** The derivation (`EQUATIONS.md` §3.6) is re-verified and is not what failed: handing both EDGEs the ORACLE truth and re-solving is worth **−0.000** off capture, and the ladder prototype is **negative** (mwae 0.0413 → 0.0426, confidently-wrong +10.7 %). `TRAPS.md` B18 |
| **a LEVEL transfer from the intron** | toy + panel | ⛔ **REFUTED**, +0.207 on capture-ON × unstranded — capture inverts which side is well-counted (`TRAPS.md` B19) |
| **the RNA fragment-length model** | `length_ceiling.py`, one pmf at a time | ⛔ **−0.02 %** at pass-0, **+0.21 % (worse)** over all objects. Root cause exact (`pi(w)` scores junction *crossing*, the pool requires the splice to be *seen*). ⭐ Its value is the BOUND: the whole C13/C14 cluster costs ≤0.43 % of the shipped solve. `TRAPS.md` B21 |
| **C2's κ residue, as an ACCURACY fix** | κ injected at exactly ½, all 36 conditions | ⛔ **−0.2 %** unstranded, worse on the shipped solve. ⭐ But the *general* defect — a boolean licence flipped by a small residue — is step 3, and the destruction control taught `TRAPS.md` A12 |
| **a nascent-bearing ladder condition** | toy, 36 conditions × 7 rungs | ⚠ **−5 %**, and the wrong way on one stratum. Keep it as a harness arm (`--nrna 60`); it no longer justifies re-simulating the panel |
| **a threshold anywhere in the licence family** | B11 implemented and refuted one | ⛔ τ is continuous across the region, so any floor is a tuned constant (`TRAPS.md` B11, D4f, D4g — refused three times) |

---

## §4 THE OPEN QUESTIONS, ranked by how much they would change a decision

1. ⭐⭐ **Is the phantom-gDNA floor the same bug as the vertex problem?** §2. Doubles or halves the value of
   step 1.
2. ⭐⭐ **What is the joint price of the splice-flux reframe and the level defect?** Step 2. Neither has an
   honest price alone, and one of them is already in `src/`.
3. ⭐ **Do we solve ALTERNATIVE SPLICING correctly?** The `alt_splice` rung exists and is **unverified** — see
   the handoff. Cheap, and it is the only structure where several junctions share an EDGE.
4. ⚠ **The two capture-ON pilot rows disagree about the SIGN of every length correction**, across two
   independently built panels. Unexplained. ⛔ Do not average them; find which one is lying.
5. ⚠ **Does the reframe's own `σ²_transfer` correctly price a ratio built on 19 counts?** `EQUATIONS.md`
   §3.5d says it is the right medicine for the noise and the wrong medicine for the share-weighting.
   Unmeasured.

---

## §5 THE RULES THAT MADE THESE NUMBERS TRUSTWORTHY

⭐ These are `TRAPS.md`'s operational core, repeated here because they are what a new session must do on day
one, not read about later.

- **Measure the CEILING before building the correction** (`TRAPS.md` B1). It has re-ranked this project
  **five** times, and it is also how you learn a phase is finished — §0 is that story for Stage A.
- **Run the PANEL arm before writing a mechanism into `src/`** (B18). Three toy-positive changes have been
  panel-negative.
- **Zero-gDNA and zero-RNA controls on every experiment** (owner, 2026-08-05). The truth is a constant, so
  every deviation is a false positive with nothing to cancel it — and step 1 was found this way.
  ⚠ But check the arm could have fired at all (A14).
- **Quote `mwae` over ALL objects and the raw Σ|err|**, never the honesty metrics alone (A12).
- **One thing varied**, with the baseline re-recorded from the current tree in the same session (B8).
- **A falsification test first, verified failing — then break the fixed code and watch each gate fire**
  (A2). And name the observable for *each place* the change was made (B14).
