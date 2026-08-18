# THE PLAN — FIX THE RELAY, THEN GIVE THE EXON A REFERENCE. START HERE.

    ⚠ **A DEV DOC. Nothing may cite it, and it is NOT the state.** The state is `ROADMAP.md` §0, the
    rulings are `DESIGN.md`, the derivations are `EQUATIONS.md`, the lessons are `TRAPS.md`. Everything
    settled has already MOVED there. ⛔ When this plan is executed, DELETE this file in the same commit.

    ⭐ Written 2026-08-17 as the last act of a session, at the owner's instruction: *"I wanna spend the
    last few turns of this session planning out exactly what we should be building... an implementation
    plan that we totally agree on, we believe in, we think it's gonna be successful."*

## ⛔⛔⛔ READ THIS PARAGRAPH BEFORE ANYTHING ELSE

**You are not going to invent a mechanism. FOUR sessions have already invented it, independently, and
each one spent most of its turns getting there.** The mechanism is MESSAGE PASSING, it is already built
in `src/rigel/calibration/messages/`, and it is switched off because of ONE NAMED BUG. Your job is to fix
that bug and give the exon a reference prior — not to design a way around the relay.

⭐ If you find yourself reasoning about anchors, ladders, pooled densities, local imputation, run-fill,
throttles or bounds **without having said the words "message passing"**, stop and read
`TRAPS: we-keep-re-deriving-message-passing`, `DESIGN.md`'s message-passing section, and `EQUATIONS.md`'s
capture-reference section. They contain the whole derivation. It cost four sessions; do not pay again.

## §1 WHAT IS ALREADY TRUE — do not re-measure these

| | |
|---|---|
| ✅ **the λ bracket is LANDED and unconditional** | `L = max(reference floor, the fitted prior's own demand)`, derived as the landscape's own log-dynamic range. All four strata improved on the 16-condition panel. No flag |
| ✅ **capture-OFF is SOLVED by a per-object reference** | `m_i = ρ_g,i·E_g,i / M_i` with ρ_g locally imputed. Measured **0.4037–0.4977** vs base on the contaminated capture-OFF conditions, and the `g00` control to 0 |
| ✅ **off-target objects are SOLVED** | intergenic and intron REGIONs read `m_i` 0.9878–0.9916 against a truth of 1.0000 |
| ⛔ **capture-ON at EXONS is the ONLY gap** | every pooled or local gDNA anchor is 2.6–3.6× low there, because probes enrich exons and every anchor sits outside or at the edge of the footprint |
| ⛔ **the relay's composition licence is BUGGY** | it checks transcript TERMINI only; `mrna_active` flipping is a population change too and is unchecked. Strict xfail in `tests/calibration/test_structural_reference.py` |

⛔ **AND THE THINGS THAT ARE REFUTED — do not retry them.** A pooled scalar reference (3.90× worse on
stranded × capture-ON), a nearest-rung one (1.27×), a locally-imputed one (1.50×), and a
variance-throttle on any of them (`w` never fell below 0.635; the defect under capture is BIAS, and
widening a prior cannot repair bias). `ROADMAP.md` §4/§4.1/§4.2 carry eighteen more. Read them.

## §2 THE BUILD, IN ORDER. Each step has a gate, and step N+1 is not started until step N's gate is green.

### STEP 0 — SOLVE THE EXONS WE CAN, AT PASS-0, ON A TIGHTLY SCOPED POPULATION (owner, 2026-08-17)

⭐ **This is the bridgehead, and it comes BEFORE turning the relay on.** Rebuild message passing from the
ground up, confined to exactly one structure:

    intron REGION  <->  exon|intron BOUNDARY  <->  exon REGION

**THE CRITERIA — an exon REGION is solvable at pass-0 iff:**
1. it is **SINGLE-STRANDED** — `free_pos XOR free_neg`, never AMBIG (an AMBIG exon has two RNA
   populations and its composition is 2-DOF);
2. it has **at least one adjacent `exon|intron` BOUNDARY**;
3. that boundary **qualifies**: single-stranded, **NOT a TSS/TES** (a terminus is a population change —
   `terminus_flank_gain` already knows this), and **carries splice-junction flux** (`sj_count > 0`).

⭐ **MEASURED REACH, 2026-08-17** — the qualifying share of EXON mass:

    g50 ss0.99 capture_off   7,120 exons   570,196 of 2,451,080   23.3 %
    g50 ss0.50 capture_off   7,133 exons   569,893 of 2,450,555   23.3 %
    g50 ss0.99 capture_ON    5,373 exons   996,109 of 5,375,588   18.5 %
    g98 ss0.99 capture_ON    2,916 exons   530,100 of 5,558,330    9.5 %

⛔ **THE ONE CAUTION TO CARRY: the `sj_count > 0` clause self-limits exactly where gDNA dominates.**
`g98` capture-ON falls to 9.5 % because there is little RNA left to splice. Those exons are near-pure
gDNA anyway, so it may not matter — but it means the bridgehead is THINNEST in the regime that needs it
most, and that must be checked rather than assumed.
⚠ Criterion 1 is ANNOTATION-derived, not data-derived: the `ss_0.50` row reads the same 23.3 % as its
stranded twin, which confirms it.

**Why it does not need to cover every exon:** the point is a SUBSTRATE. Solved exons train the gDNA
landscape; the landscape generalises to the rest. ⭐ So the gate on this step is not coverage — it is
whether the solved exons SPAN the enrichment range well enough for the landscape to fit both modes.

**Gate.** Both zero controls; per stratum; and the solved exons' `f_g` scored against the origin-split
oracle *on the qualifying population only*, so a partial win is not diluted by the exons it never
claimed.

### STEP 1 — FIX THE COMPOSITION LICENCE (the blocker for turning the RELAY on)

**What.** `HeadPolicy`'s composition licence must refuse to impute a composition across a hop where
`mrna_active` FLIPS, exactly as it already refuses across a transcript terminus (`terminus_flank_gain`).
The predicate exists and is the solver's own; the licence simply does not consult it.

**Why first.** With the bug live, an intron's correct "pure gDNA" claim is relayed into the adjacent exon
and drives it to a confident wrong vertex. Every relay measurement taken while that is true — including
the recorded **+154.8 %** — is a measurement of the bug, not of the relay.

**Gate.** The strict xfail in `tests/calibration/test_structural_reference.py` goes green **by the code
being fixed**, never by widening a bound. ⛔ Then break the fix and watch it fail again
(`falsification_needs_perturbation`).

**Watch for.** This is half of `TRAPS: a-cancelling-defect-pair`'s family — check `ROADMAP.md` §4.1's
`zc_struct_lock_g1` row before assuming the licence change is independently priceable.

### STEP 2 — RE-PRICE `message_propagation = True`, WITH THE LICENCE FIXED

**What.** `ladder_arm_ab.py --arm base --messages off` against `--messages on`, `--jobs 8`, all 16
conditions, per stratum, both zero controls. Then `calibration_vs_oracle.py`.

**Why.** The owner's 2026-08-10 ruling (off until the tool is optimised end to end) is not re-opened by a
measurement — but the number that justified it is stale in two independent ways: the licence bug was live
and the panel it was measured on was retired 2026-08-13. **RE-PRICE, NEVER INHERIT.**

**Gate.** Per stratum, never pooled. The `g00` controls are owner-required. ⚠ The recorded price on the
zero-gDNA golden scenarios was 0.029 → 89.93 and 0.005 → 9.58 — if that survives the licence fix, the
relay is still not shippable and step 3 must carry the exon on its own.

### STEP 3 — THE EXON REFERENCE: SPIKE-AND-SLAB

**What.** `EQUATIONS.md`'s capture-reference section is the derivation; implement it as the value of
`CompositionPriors.location`'s successor. The term's SHAPE changes (a mixture, not a point location), so
this touches `simplex_logodds._location_term`.

    p(log ρ_g,i) = π·N(log ρ_0, σ_0²) + (1−π)·Unif[ log ρ_anchor,i , log ρ_max,i ]
    ρ_max,i      = ( M_i − (ρ_r·E_r,i − S) ) / E_g,i

⛔ **`ρ_max,i` is the sj-flux bound in density coordinates, NOT the trivial cap `M_i/E_g,i`** (which is
just `f_g ≤ 1`). It degenerates to that cap where there is no flux to subtract. A draft of this plan wrote
the cap and the bound as if they were the same thing — they are not, and `EQUATIONS.md` §9d.4 carries the
corrected form. ⚠ And the atom is the median only when it carries `π ≥ ½`; at `g00` the slab still spans
to 0.6039 in `f_g`, so a MINORITY spike does not reach the vertex.

**Where it goes.** The mixture's evaluation is layer 3 (`simplex_logodds`). Its PARAMETERS come from
neighbours and therefore from layer 6 (`messages/`) — which is the whole point, and the reason step 1
comes first.

**Gate.** ⛔ It must reduce EXACTLY to the shipped form at capture-OFF (measured anchor ratio 0.98 ⇒ the
slab collapses onto the spike) — a `noop`-equivalent byte-identity check on the capture-OFF conditions.
Both zero controls. All 16, per stratum. And a SHUFFLE control: a deliberately permuted `ρ_anchor` must
NOT win (`TRAPS: attribution-must-survive-a-shuffle`).

**⭐ The prize that makes this worth doing, and it is a theorem rather than a hope.** A Beta prior cannot
put its median at `f_g = 0` — `EQUATIONS.md` §9a. A spike-and-slab has an ATOM and can. So this is the
only proposed mechanism that can reach the vertex at all, and the "≈¾ of the vertex shortfall is
irreducible" framing is an artifact of the prior FAMILY.

### STEP 4 — THE LANDSCAPE TRAINS ON A REAL SUBSTRATE

**What.** Nothing to build — this is the payoff, and it is why the exon matters at all. With exons
bracketed rather than pinned, `_fit_gdna_hyperprior` trains on objects whose `f_g` reflects data, and it
can fit the bimodal `P(log ρ_g)` that hybrid capture actually produces.

**Gate.** The landscape must show mode SEPARATION on a capture-ON condition and a single pile at the wall
at `g00`. ⭐ Both are already recorded as achievable (`ROADMAP.md` §4: 2.98 decades at `g75 ss0.99 ON`,
30× more enriched mass ON than OFF) — so this is a reproduction, not a discovery.

**⚠ The circularity to watch.** `_fit_gdna_hyperprior` trains on `belief.f_g * mass` over expressed
REGIONs including exons. If step 3's reference is too strong the landscape learns the prior back. The
measurement that says whether this is happening: at `g50`/`g98` capture-ON, **~80 % of the training mass
has its own composition evidence** — so the landscape has real data to learn from *provided the reference
does not overwrite it*. Measure the no-evidence share of training mass before and after.

### STEP 5 — PERFORMANCE, AND IT IS MANDATORY BEFORE 0.8.0 (owner)

**What.** The grid solve in memory-bounded parallel chunks, with an advanced CLI flag spanning
one-object-at-a-time through many-objects-in-parallel.

⛔ **OPTIMISE ON HIGH-DEPTH REAL RNA-SEQ, NOT cfRNA** (owner, 2026-08-17) — cfRNA is too sparse and small
to be the target. ⚠ This REVERSES the standing "profile on real cfRNA" instruction, which still appears
in several places; the sites are listed in the scripts-track report and need correcting.

⚠ `scripts/profiling/` is the natural home and is currently covered by NO gate (`ALL_SCRIPTS` is
`design + sim`); 2 of its 3 files were found dead this session. Its resurrection did not run — its agent
hit an API error — so it is still open.

## §3 THE STANDING DISCIPLINE — the things that actually caught defects this session

* ⛔ **`TRAPS: panel-before-src`.** Five times now a toy- or single-condition-positive change has been
  panel-negative. Price on all 16, per stratum, before `src/`.
* ⛔ **A monkeypatch arm must ASSERT ITS TARGET EXISTS and COUNT ITS FIRINGS before any number is read.**
  `rigel/calibration/__init__.py` does `from .calibrate import calibrate`, so
  `import rigel.calibration.calibrate as CAL` binds the FUNCTION — patching `CAL.solve_chain` is inert
  and a 314-second run reported "all arms byte-identical", which read as a real negative.
* ⛔ **Derive the failure set, never eyeball the tail.**
  `python -m pytest tests/ -q 2>&1 | grep '^FAILED' | sed 's/::.*//' | sort | uniq -c`
* ⛔ **Re-derive the suite count, never adjust it.** `pytest --collect-only -q | grep <stem>`.
* ⛔ **Both zero controls on every arm**, and report in FRAGMENTS, per stratum, never pooled.
* ⚠ **`g00` can be right for the wrong reason.** When every anchor density is 0, an arm can assert
  `f_g = 0` unconditionally through a degenerate clip. The answer is right; the control is testing
  nothing. Check the code path before quoting it.

## §4 THE TREE IS CLEAN — everything below is COMMITTED

    76702b91  docs: re-derive the standing baseline on a clean tree
    4cfd35d9  scripts: retire the profiling tree's legacy, and let the gate cover it
    c2599979  docs: lock down message passing so it is not re-derived a fifth time
    3425e341  tests: restore the fixture ids and gate labels the bulk rename overwrote
    21bf6aca  scripts: fix toy_dissect's leaked loop variable — a printed density moves 2.25x
    46bb8ee1  scripts: resurrect ten dead instruments, and delete two that cannot run
    648cf9a1  calibration: derive ψ's λ bracket from the fitted prior's own support

**Baseline: 3,501 passed / 10 xfail / 0 skipped / 0 failed, 3,511 collected.** `ruff` clean.
⛔ Re-derive it rather than trusting this line (`TRAPS: re-record-the-baseline`).

⚠ **Every reference experiment from that session was PROBE-ONLY and is NOT in the tree.** Nothing is
half-built in `src/`, and there is nothing to revert. The probes live in the session scratchpad and are
gone; the MEASUREMENTS they produced are in `EQUATIONS.md` §9d, `DESIGN.md` §0c and `ROADMAP.md` §0.

⚠ **One judgement to know about, because it is a live disagreement resolved against my first instinct.**
Mid-session I restored `scripts/profiling/pyspy_driver.py` and 8 config YAMLs, reading their deletion as
over-pruning ahead of the mandatory performance work. `scripts/README.md` documents the opposite and is
right: the configs name a `--config` flag no instrument has and keys absent from `src/rigel/config.py`,
and `py-spy record -- python scripts/profiling/profiler.py` replaces the driver exactly. The tree was
never deleted — **the GATE was extended to it**, which is the thing that was actually missing. If you
want the driver back, that is an owner call and `4cfd35d9^` has it.
