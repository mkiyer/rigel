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
| ⛔ **stranded** | 2,393,419 | **2,606,346** | **+8.9 %** · 7/18 — the residual, now **ATTRIBUTED**, below |

⭐⭐⭐ **THE STRANDED × CAPTURE-ON RESIDUAL IS ATTRIBUTED, AND IT IS THE CAPTURE LEVEL RESIDUAL** (2026-08-06,
second pass). **100 %** of it is owned by ONE term — `composition_logvar`'s counting term, i.e. by
`σ²_transfer` no longer being `∞` at a zero-mass slot. Reverting only that term reproduces the "before"
column **to the fragment** (295,453 / 304,815 on the node axis), confirmed independently against a from-git
pre-fix tree to 0.17 % (`TRAPS.md` A5 applied and passed). ⭐ The mechanism: **`1/n = ∞` made a zero-mass
slot a relay BARRIER, and `count_logvar` made it a CONDUIT** — under capture the off-probe gaps between probe
islands *are* the zero-mass slots, so a gDNA level now crosses between capture strata unscaled. ⛔ It is
PROPAGATION, not origination: muting every zero-mass emitter recovers **4–27 %** of the regression, and
restoring the barrier recovers **100 %**. **Six candidate fixes priced, six refused** — five by the `g00` zero control (+3,207 % to
+7,269 %) and one as inert. `TRAPS.md` C0e/A11b/D4k, and `SESSION_HANDOFF.md` §4 carries every number so
none is re-tried.

⛔⛔⛔ **AND THE 39 % IS A PASS-0 NUMBER THAT IS ONLY WORTH −3.9 % ON THE DELIVERABLE** (measured 2026-08-06,
both arms, one session, `abs_err_all_final` — the shipped solve after the gDNA prior is fitted and it
re-solves). The column was in the data all along and had never been read:

| | pass-0 | **final (shipped)** |
|---|---|---|
| **PANEL** | −37.2 % | ⭐ **−3.9 %** |
| unstranded | −42.9 % | −9.9 % |
| `g00` zero control | −98.0 % | −97.0 % |
| ⛔ **stranded** | +14.9 % | ⛔ **+30.0 %** |
| ⛔ stranded × capture-ON | +28.5 % | ⛔ **+38.7 %** |

⭐⭐ **The prior refit was already doing most of the work the pass-0 fix does**, so the win attenuates 10×.
⛔⛔ **And the regression AMPLIFIES through the refit** — the "pass-0 is the training set" mechanism running
backwards: a pass-0 that under-calls gDNA on high-gDNA stranded conditions fits a prior that under-calls it
there, and the refit compounds the bias instead of correcting it. **A prior cannot rescue a substrate biased
in its own direction.**
⛔ **CONSEQUENCE FOR HOW THIS PROJECT IS MEASURED: pass-0 `Σ|err|` is only weakly coupled to what ships, and
nobody had checked the coupling.** Quote both columns from now on. `SESSION_HANDOFF.md` §3 and `TRAPS.md` A15.

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
| ⭐⭐⭐ **0** | ⭐⭐⭐ **CONSOLIDATE ONTO A BACKBONE — and the measurement that sizes it is DONE** | a restructure gated on BYTE-IDENTITY, then one arm per operator | ⛔⛔ **The whole message layer at pass-0 is worth +0.2 % on the deliverable** (`--arm msgfree_p0`, 36 conditions, `g00` BYTE-IDENTICAL) while pass-0's own Σ\|err\| moves 77.5 % — so the arm demonstrably fired and the shipped answer did not care (A14 in its strongest form). ⭐⭐ And muted EVERYWHERE (`--arm msgfree_all`) the messages are worth −50 % overall but are a **net HARM on three of the four strata**: stranded × capture-ON **−58.3 % (16/16)**, stranded × capture-OFF **−43.7 % (16/16)**, unstranded × capture-OFF **−32.1 %**. Their entire value is **unstranded × capture-ON**, +154.8 %, 0/16 — the stratum where κ = ½ makes the strand λ-term exactly 0 and messages are the only composition evidence there is. ⛔⛔ **So the standing stranded × capture-ON residual is not a missing mechanism: deleting the message layer on that stratum is −58.3 %.** Eleven candidates were tuned on a panel where 3 of 4 strata are harmed by the thing being tuned. `docs/calibration/CONSOLIDATION.md` |
| ⭐⭐⭐ **1** | ⭐⭐ **RE-ANCHOR THE PROJECT ON THE DELIVERABLE, NOT ON PASS-0** | measurement + a scoring ruling, no build | §0: a **−37.2 %** pass-0 win is **−3.9 %** shipped, and the stranded regression **grows** +14.9 % → +30.0 % through the refit. ⛔ Until this coupling is understood, every pass-0 ranking in this file is suspect. ⭐ Cheapest first step: re-read `abs_err_all_final` for the campaigns already on disk — it costs nothing and may re-rank several closed items. Then decide whether `solvability_audit`'s headline should be the final solve |
| ⭐⭐⭐ **2** | ⭐⭐ **THE `η` REBUILD'S `g00` MECHANISM IS FOUND AND FIXED — PRICE THE CERTIFIED-RNA CHANNEL'S *COVERAGE* NEXT** | one panel arm | ⭐ **The mechanism was a UNIT ERROR in the tilt channel** — a raw log-odds delivered into ψ's `theta_imp`, whose grid is the ANGLE `arcsin(τ)` spanning exactly `[−π/2, π/2]`; the modes were ±4.6, so the tilt pinned at the boundary and destroyed the AMBIG Schur protection. **74 % of `g00`'s error, in the one channel the prototype did not publish into `_capture`** (`TRAPS.md` A16). Fixed, gated (33 tests, 8/8 firing perturbations): `g00 ss0.99 cap_ON` Σ\|err\| **1,516,976 → 403,073**, panel **+103 % → +85.1 %**. ⛔ Still negative, but the residual is now MONOTONE IN gDNA CONTENT on unstranded data — −76 % at `g98`, +3,806 % at `g01` — i.e. a systematic upward `f_g` bias on the stratum where κ = ½ silences strand and the MEASUREMENT channels carry everything. ⭐⭐ **The suspect is WHICH STREAMS ARE IN THE RELAYED STATE**: HEAD relays ten arrays including the three MEASUREMENT precisions, so certified-RNA evidence travels the chain; the rebuild relays four, and delivers the certified-RNA measurement one step to exon NODEs only (EDGEs receive it 0.0 % of the time). ⛔ **NOT a forward-backward-vs-Gauss-Seidel difference — an earlier draft of this row said so and was WRONG** (`TRAPS.md` A17); both sweeps are one forward pass plus one backward pass and both combine from the neighbour's state. ⭐ Unlike the eleven in the graveyard this is a claim about coverage, not about how to resolve doubt. `SESSION_HANDOFF.md` §1/§3 |
| ⭐⭐ **3** | **THE CAPTURE LEVEL / STRANDED × CAPTURE-ON RESIDUAL** (`EQUATIONS.md` §3.5c) | ⛔ **NOT a pass-0 patch — that route is now closed by measurement** | ⭐⭐ Steps 1 and 2 of the old list were ONE item, as that table guessed. ⛔⛔ **ELEVEN candidates have now been priced and all eleven died on the `g00` zero control or were inert** — including the `η` rebuild itself and `variance_model_notes.md`'s own §2 (`variance_model_notes.md` §6, `SESSION_HANDOFF.md` §6). The reason is structural: the error mass sits on objects with **zero own precision**, so only the MODE can move it — and `D` is ~10⁶ at `g00` where the floor is right and ~10 at `g75` where a lift is right, so **no fixed function of the observables sets it**. The separating quantity is the library's gDNA content, which pass-0 excludes by construction. ⭐ So the work is at the REFIT, and its precondition is an unbiased pass-0 mode that does not yet exist. ⚠ The 0.2 %/270 % numbers pre-date this campaign and must be re-derived |
| ⭐⭐ **4** | **PRICE THE SPLICE-FLUX REFRAME JOINTLY WITH THE LEVEL DEFECT** (`EQUATIONS.md` §3.6c) | nothing to design: one is in `src/`, the other is `toy_trace_error.py`'s relay-level arm | The reframe is correct by derivation and DEFAULT, but alone it is panel-negative on the node axis. `TRAPS.md` D4j: half of a **cancelling pair**, so neither has an honest price alone. ⚠ Its numbers pre-date the 39 % win and must be re-measured |
| ⭐ **5** | **THE `intergenic\|exon` EDGE ↔ intergenic NODE PAIR AS A gDNA REFERENCE** | small build, ceiling first | Owner-requested. ⚠ **Do NOT approach it as a pass-0 mechanism** — that route is closed (§0, `variance_model_notes.md` §6). A G1 EDGE has `tau_lam = 0` so it can only RELAY a level, and the mass pin's case (ii) is what gives it one. ⭐ Under the `η` rebuild the pin is DELETED, so this item is subsumed by it: do it inside the rebuild, not separately |
| **6** | ⚠ **Add a 1–10 % gDNA rate to `pilot.yaml`** | one config line + a re-run | `suite_resolves.py` requirement (f); the pilot is only 0 %/100 %. Blocks the Stage-A numbers only |
| **7** | ⚠ **The multimapper hole** | `SUCCESS.md` A3 | Cost still **unmeasured**; measure it as a ceiling first |

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
| **the gDNA prior's BIMODAL CAPACITY, and "give the prior more signal"** | a read of `gdna_landscape.py` + the production refit on real conditions | ⛔ **BOTH BRANCHES CLOSED.** The prior already renders the landscape correctly — **2.98 decades** of mode separation at `g75 ss0.99 capture_ON`, 30× more enriched mass ON than OFF, a single pile at the wall at `g00`. And a prior fitted from ORACLE truth is the same prior (0.04 dec). Not capacity, not signal, not location. `SESSION_HANDOFF.md` §2/§3 |
| **§2's Jeffreys MEAN density location** | `--arm eta`, the `g00` zero control | ⛔ **REFUTED at +96,299 %.** It cannot say ZERO (`node_init.rho_g` is an exact 0 at 60,544/70,176 slots — the statement earning the −98 % at `g00`), and the C0d benefit it was credited with belongs to **§4**. ⭐ If revisited the derived form is the Gamma **MODE** `max(a−½,0)/E`, which is exactly 0 at a zero count |
| **a threshold anywhere in the licence family** | B11 implemented and refuted one | ⛔ τ is continuous across the region, so any floor is a tuned constant (`TRAPS.md` B11, D4f, D4g — refused three times) |

---

## §4 THE OPEN QUESTIONS, ranked by how much they would change a decision

1. ⭐⭐⭐ ~~**WHY IS PRIOR FIDELITY ANTI-CORRELATED WITH DELIVERABLE QUALITY?**~~ ⭐⭐ **LIKELY ANSWERED,
   AND THE ANSWER IS THAT IT IS NOT THE PRIOR — IT IS THE MESSAGES** (2026-08-06, ψ-boundary ablation on
   HEAD at `g01 ss0.50 capture_on`, A5 exact). At every one of HEAD's 12 worst slots the **self-solve WITH
   the fitted prior is nearly correct** (0.0043–0.0344 against truths of 0.0007–0.0136) and **the message
   layer then destroys it** (0.05–0.83). Muting the certified-RNA channel alone recovers the truth at 8 of
   the 12. So a better prior cannot show up in the deliverable: it is already right at the objects that
   carry the error, and something downstream overwrites it. ⛔ That also re-reads the earlier evidence
   rather than contradicting it — "the prior is degraded on unstranded where the solve improves" is what a
   prior looks like when it is not the binding constraint. ⚠ Confirm on a second stratum before closing.
   `SESSION_HANDOFF.md` §3.2b, `TRAPS.md` A15.
1b. ⭐⭐ ~~**Is the phantom-gDNA floor the same bug as the vertex problem?**~~ ⛔ **ANSWERED, and NO** — it was
   a different and much larger bug (§0's 39 %), and its residual is the CAPTURE LEVEL defect, not the vertex.
   `TRAPS.md` C0c/C0d/C0e.
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
