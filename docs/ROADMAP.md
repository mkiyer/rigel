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

**Stage B — calibration.** ⭐⭐ **A 39 % PANEL-WIDE ERROR REDUCTION LANDED 2026-08-06** (`80254008`,
`0ed8eceb`). Σ|err| over every LIVE object, pass-0, all 36 ladder conditions, against origin-split truth:

| | Σ\|err\| before | after | |
|---|---|---|---|
| **THE PANEL** | 24,682,284 | **15,048,672** | **−39.0 %** · 23/36 better |
| unstranded | 22,288,865 | 12,442,326 | −44.2 % · 16/18 |
| capture_OFF | 7,060,421 | 2,599,789 | −63.2 % · 12/18 |
| capture_ON | 17,621,863 | 12,448,883 | −29.4 % · 11/18 |
| ⛔ **stranded** | 2,393,419 | **2,606,346** | **+8.9 %** · 7/18 — the residual, `SESSION_HANDOFF.md` |

⭐ **The zero-gDNA controls are now essentially exact** — truth is 0 there, so every fragment was a false
positive: `mwae_all` **0.3776 → 0.0076** (`g00 ss0.50 capture_on`), **0.3411 → 0.0080** (capture_off),
**0.0172 → 0.0003** and **0.0128 → 0.0002** on the stranded pair. And the win runs far past `g00` along
the whole low-gDNA unstranded axis: `g01` −85 %, `g05` −72/−77 %, `g10` −58/−73 %, `g25` −28/−62 %.

⚠ **The per-axis `mwae` figures this table used to carry (node 0.1192 / edge 0.1350 at pass-0, 0.0586 /
0.0633 shipped) are STALE and are not restated here** — they pre-date the fix and nobody has re-derived
them. Re-run `solvability_audit.py --suite` rather than quoting a number from memory.

⛔ **And read `mwae_all` / `Σ|err|`, never `solv%` / `mwae` / `conf-wrong` / `calib`.** Those four share a
denominator the solver moves by declining to answer, and it flips on fitting noise at three conditions
(`TRAPS.md` A12b). The two fixed-denominator columns landed in `a3d45031` for exactly this reason.

## §1 ⭐⭐⭐ WHAT TO DO NEXT

| | step | shape | why now |
|---|---|---|---|
| ⭐⭐⭐ **1** | **THE STRANDED × CAPTURE-ON RESIDUAL** | dissect → localize → fix, the standard loop | `SESSION_HANDOFF.md` has it in full. Every `ss_0.99 capture_on` row at `g10`+ is **20–34 % WORSE** after the 39 % win — one coherent stratum, not noise. Hypothesis (UNPROVEN): under capture an empty intergenic node means *no probe here*, not *no gDNA here*, so the new zero-count anchor asserts `ρ_g ≈ 0` where it should assert nothing. Worst row `g75 ss0.99 capture_on`, **+100,599** |
| ⭐⭐ **2** | **THE CAPTURE LEVEL RESIDUAL** (`EQUATIONS.md` §3.5c) | a build, ceiling-free dissection first | Still the largest untouched axis: an exon's own total density equals its true gDNA density to **0.2 %** (923× vs 921×) while the relay delivers **270×**, because it is measured at the flanking EDGE whose fragments straddle the probe boundary. ⚠ Plausibly the SAME root as step 1 — both are "capture broke the opportunity model" — so do step 1 first and check whether it subsumes this |
| ⭐⭐ **3** | **PRICE THE SPLICE-FLUX REFRAME JOINTLY WITH THE LEVEL DEFECT** (`EQUATIONS.md` §3.6c) | nothing to design: one is in `src/`, the other is `toy_trace_error.py`'s relay-level arm | The reframe is correct by derivation and DEFAULT, but alone it is panel-negative on the node axis. `TRAPS.md` D4j: half of a **cancelling pair**, so neither has an honest price alone. ⚠ Its numbers pre-date the 39 % win and must be re-measured |
| ⭐ **4** | **THE `intergenic\|exon` EDGE ↔ intergenic NODE PAIR AS A gDNA REFERENCE** | small build, ceiling first | Owner-requested. Mechanism identified: `struct_lock` is NODE-only, so a G1 EDGE gets `Var = ∞` and can only RELAY a level, never ORIGINATE one. ⚠ **Re-check first** — the zero-count fix may already have moved this |
| **5** | ⚠ **Add a 1–10 % gDNA rate to `pilot.yaml`** | one config line + a re-run | `suite_resolves.py` requirement (f); the pilot is only 0 %/100 %. Blocks the Stage-A numbers only |
| **6** | ⚠ **The multimapper hole** | `SUCCESS.md` A3 | Cost still **unmeasured**; measure it as a ceiling first |

⛔ **`Var(κ̂)` / the strand deadband is CLOSED as an accuracy item** — measured a wash (−0.6 / −5.3 / +4.6 %
across the three big classes) and it is a REPORTING defect, fixed instrument-side in `a3d45031`
(`TRAPS.md` A12b). Do not re-open it as a solver change; three designs were refuted.

## §2 ⛔ THE VERTEX PROBLEM — CLOSED, kept as one paragraph so it is not rebuilt

ψ lands 5–8 % short of `f_g = 1` on unexpressed genes, and that was ranked #1 for a campaign. It is
**not a build**: `vertex_ceiling.py` prices it and evidence-free objects sit at median `|Δ|/sd(f_g)` =
**z = 0.5–0.6**, inside their own 1σ — the 24.4 % it measures is *the value of missing information, not
headroom*. No prior-free solve can beat it: every proper prior on [0,1] has a median strictly inside
(0,1). ⭐ `test_certified_rna_licence.py` `C2b` independently closed the one channel that might have
supplied vertex evidence — a zero certified count is consistent with `f_g = 1` too. ⚠ Its companion
hypothesis, the **phantom-gDNA floor**, turned out to be a DIFFERENT and much larger bug and is now
fixed: `TRAPS.md` C0c/C0d, §0's 39 %.

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
