# SESSION HANDOFF — ⭐⭐⭐ THE TWO FACES OF AN `intron|exon` EDGE: derived, verified, NOT BUILT

    Written 2026-08-04. ⚠ **WORKING DOC, NOT A PERMANENT FIXTURE.** The six permanent docs are
      `CLAUDE.md`, `docs/SUCCESS.md`, `docs/ROADMAP.md`, `docs/TRAPS.md`, `docs/EQUATIONS.md`,
      `docs/DESIGN.md`, `docs/TESTING.md`.
    ⛔ **DELETE THIS FILE when its steps land**, promoting anything worth keeping into `ROADMAP.md`
      (a current number), `TRAPS.md` (a lesson), `EQUATIONS.md` (a derivation) or `DESIGN.md` (a ruling).

## Read in this order

| | |
|---|---|
| 1 | `CLAUDE.md` — doc map, working rules, the script table |
| 2 | ⭐⭐ **`docs/DESIGN.md` §0 — THE BINDING VOCABULARY.** NODE, EDGE, step, the mass pin, structurally pure-gDNA object. Read it before writing a comment |
| 3 | ⭐⭐⭐ **`docs/EQUATIONS.md` §3.6 — THE DERIVATION THIS HANDOFF IS ABOUT.** It is verified against oracle truth; you are building it, not deciding it |
| 4 | ⭐⭐⭐ **`docs/TESTING.md` §0b** — the toy harness, the `spliced_exons` rung, `toy_panel.py`, and the capture-ON genome-length rule. **Everything you need to run the bench is there** |
| 5 | `docs/EQUATIONS.md` §3.5 / §3.5b / §3.5c — the composition licence the last two steps landed; §3.6 extends the same idea to a matched component set |
| 6 | `docs/calibration/ISSUES.md` **C12, C13, C14** (this work) and **C2, C6** (what it collides with) |
| 7 | `docs/TRAPS.md` **B15, B16, B17** — three mistakes made *this* session on the bench, all avoidable |

---

## 1. State of the tree

**Branch `fragment-length-gold-standard`, and it is COMMITTED** — a change from the last handoff. Four
commits, newest last:

| | |
|---|---|
| `c6b2ea89` | the message licence: the mass pin gated (C11) + the population conjunct (C4) |
| `b6a50e65` | the `spliced_exons` rung + `toy_panel.py` |
| `34da3d71` | probes tile PER EXON; capture-ON edges are lifted by LENGTHENING |
| *(this one)* | the docs |

**Test baseline: 22 failed / 2214 passed / 2 xfailed.** The 22 are the standing set (21
`test_golden_output` + `test_paralogs::test_gdna_sweep[gdna_100]`, `TRAPS.md` D9). A 23rd failure or any
other non-golden name is a regression. `ruff check src/ tests/ scripts/` clean.
⚠ The oracle cache at `~/Downloads/rigel_runs/suite/ladder/oracle_cache` is **VALID — do not delete.**

---

## 2. ⭐⭐⭐ WHAT IS ALREADY ESTABLISHED — do not re-derive it

`EQUATIONS.md` §3.6 has the full derivation. The one-paragraph version, because everything below turns on
it: at one line the accumulator stores three populations whose **component sets differ** —
`unspliced_count` is `{gDNA, nascent}` (mature cannot cross an exon↔intron seam contiguously),
`junction_count` is certified mature. So one EDGE carries one gDNA density and **two** composition
statements, differing only in what is in the denominator:

    (I)  INTRON face:  rho_g / (rho_g + rho_nas)             == phi_g(INTRON)
    (II) EXON   face:  rho_g / (rho_g + rho_nas + rho_mat)   == phi_g(EXON)     rho_mat = J/E_J, MEASURED

**Verified against oracle truth** on `spliced_exons` × `g50 ss0.50 capture_off`: face (I) residual
+0.000000 (`nrna_none`) and +0.0266 / −0.0198 on the nascent control; face (II) −0.0011 / −0.0031 at
m=100. ⚠ The `nrna_none` arm tests (I) only trivially — both sides are 1.000 — so **the nascent control is
not optional**.

⭐⭐ **And the reason it is worth building** — measured on the exon where pass-0 is worst (m=1):

| route to the exon's `f_g` | value | vs truth 0.458 |
|---|---|---|
| the shipped pass-0 answer | 0.625 | +0.167 |
| the EDGE's own `rho_g` (**8** counts) | 0.350 | −0.108 |
| ⭐ the INTRON's `rho_g` (**349** counts), carried by face (I) | **0.499** | **+0.041** |

---

## 3. ⛔ THE THREE BENCH MISTAKES I MADE, so you do not repeat them

1. **I called capture-ON depletion "starvation" and wrote off half the panel.** It is the signal *moving*
   to the EDGEs abutting the exon. `--genome-length 120000` turns the two `intron|exon` EDGEs from 2 and 5
   counts into 20 and 36, at the exon's own capture stratum. `TRAPS.md` B15, table in `TESTING.md` §0b.
2. **Probes tiled across the transcript put a junction-straddling (split) probe over every internal
   junction**, and `gdna_split_penalty` then suppressed exactly the edge-crossing gDNA. Fixed to tile per
   exon. `TRAPS.md` B16.
3. **`toy_panel`'s per-object ceiling is a SUBSTITUTION and understates a message SOURCE.** It priced the
   two EDGEs at 9.1 % of the gene's error while the exon they feed carries 82.7 %. `TRAPS.md` B17.
   ⛔ **This is why step 2 below demands a re-solve arm, not a substitution.**

---

## 4. ⭐⭐⭐ THE PLAN — owner-approved, in this order

### STEP 0 — RE-EARN THE CAPTURE-ON NUMBERS (do this first; everything else is scoped by it)

```bash
python scripts/design/toy_panel.py --spec spliced_exons --genome-length 120000 \
    --conditions $(ls ~/Downloads/rigel_runs/suite/ladder | grep capture_on) --out capon.jsonl
```

⛔ **Every capture-ON number anywhere in the git history of this campaign was measured at 12 kb and is an
empty chromosome.** That includes my claim that capture-ON × unstranded is the worst stratum (0.3044) —
treat it as unmeasured. Confirm the edges hold ~20–40 counts across the whole gDNA ladder, not just at
g75. ~8 min sharded. ⭐ Also worth a rung of its own: the **split-probe** geometry, which is real in a real
panel and now deliberately absent from the toy.

### STEP 1 — THE NASCENT CONTROL, as a panel arm

`--nrna 60`. Every cached condition is `nrna_none`, so an `intron|exon` EDGE's truth is exactly 1.000 and
face (I) is untestable non-trivially. This is also `ROADMAP.md` step 3 and `ISSUES.md` C10 — the panel may
be measuring a robustness corner rather than the real case.

### STEP 2 — ⛔ THE RE-SOLVE CEILING FOR FACE (I)

Pin each `intron|exon` EDGE's gDNA density to the **intron's** and RE-SOLVE (not substitute). That is the
upper bound on step 3 and it is the number that decides whether to build. ⚠ Do it on capture-OFF *and*
capture-ON (post step 0), because the two regimes reference different objects.

### STEP 3 — BUILD FACE (I): the intron↔EDGE composition transfer

Same shape as the population conjunct that just landed (`EQUATIONS.md` §3.5b): the intron and the EDGE's
*unspliced* population have the **same component set**, so the composition may cross — and the transport
must run from the **well-counted** side. ⚠ Mirror it in BOTH twins (`_relay` and `_transport`) and gate the
relay side directly — `TRAPS.md` B14 records that a combine-only gate let a relay-only deletion pass the
entire suite.

### STEP 4 — BUILD FACE (II) as an INVERSE-VARIANCE FUSE

Two estimators of the exon's RNA density, strong in opposite regimes: `rho_mat = J/E_J` (tight at high RNA
— the junction sees only ~11 % of the mature fragments, so it is 3 counts at m=1) and the exon's own mass
identity closed with the EDGE's `rho_g` (strong at low RNA, a small difference of large numbers at high).
⛔ **Fuse, do not rank** — a precedence silenced the RNA channel before (`bp_solver`'s peel-level note).

### STEP 5 — the two small ones, both measured

* **C13**: `E_J` vs the exon's `E_r` disagree by 5–10 % — the whole residual of face (II).
* **C14**: `bp_solver`'s P1d says the graft's premise is a LOWER bound; measured it is an UPPER bound by
  5–10 %. The variance term built on it is sized from the wrong direction.

### ⚠ Two open questions to settle on the way

* **The e1/e2 asymmetry.** At m=100 the second exon is given ~10× too much gDNA (pred 0.0721 vs truth
  0.0073) while the first is nearly exact, in both strand regimes. ⚠ **Possibly not a separate mechanism**:
  the two flanking EDGEs had 13 and 8 gDNA counts and e2 is the one behind the 8-count edge, so it may just
  be the variance problem step 3 fixes. Confirm before treating it as its own defect.
* **Does the prior already hide any of this?** Everything above is `--refit-iters 0`. The shipped solve runs
  3 refits and has never been compared on this rung. One run tells you whether you are fixing something the
  hyperprior already covers.

---

## 5. ⭐⭐ THE BENCH — every command

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
export OMP_NUM_THREADS=1
SUITE=~/Downloads/rigel_runs/suite/ladder
CACHE=$SUITE/oracle_cache            # ⚠ VALID — do not delete

# ⭐ ONE condition, every object beside per-object truth — 13 s including the harvest
python scripts/design/toy_harness.py --spec spliced_exons \
    --donor gdna_g50_ss_0.50_nrna_none_capture_off

# ⭐⭐ THE PANEL: 36 conditions x 7 RNA rungs, PRIOR-FREE pass-0, per object
python scripts/design/toy_panel.py --spec spliced_exons --out rows.jsonl
python scripts/design/toy_panel.py --report rows.jsonl          # re-aggregate, no re-measurement
#   ⛔ capture-ON: add --genome-length 120000.   ⭐ the nascent control: --nrna 60

# the regression suite for the licence that already landed
python -m pytest tests/calibration/test_relay_mass_pin.py \
    tests/calibration/test_terminus_population_licence.py \
    tests/calibration/test_gdna_scale_rule.py -q
python -m pytest tests/ -q      # must read 22 failed / 2214 passed / 2 xfailed
ruff check src/ tests/ scripts/ && ruff format src/ tests/   # ⚠ NEVER format scripts/

# ⛔ the acceptance test for any SOLVER change — single-process, ~25 min, read weak% BEFORE mwae
python scripts/design/solvability_audit.py --suite $SUITE --oracle-cache $CACHE
```

⚠ **`solvability_audit.py` is single-process and ~25 min. Do NOT edit `src/` while it runs.** ⭐ It shards
cleanly by condition if you need it faster — one process per subset, then concatenate the rows.

⭐ **Re-record the baseline from the current tree in the same session as any arm** (`TRAPS.md` B8), and
⚠ **a falsification test first, verified failing — then break the fixed code and watch each gate fire.**
That second half found a real hole in each of the last two steps.

---

## 6. What this session measured, for the record

| | |
|---|---|
| ⭐⭐⭐ **the message licence landed** | the mass pin gated by "no BELIEF may enter its budget" (C11) + the terminus POPULATION conjunct (C4). Ladder: confidently-wrong 21,154 → 20,173, overconfidence 2.788 → 2.688, `weak%` 2.91 → 2.80. `ROADMAP.md` §2c |
| ⛔⛔ **the C11 CEILING said it cost the panel nothing** | deleting the pin outright moves the ladder **+0.0002 (worse)**. It landed on the derivation and on being free. Fourth time the ceiling discipline has re-ranked something |
| ⭐⭐ **the two-face derivation verified** | `EQUATIONS.md` §3.6, both faces, with and without nascent |
| ⭐ **the `spliced_exons` rung + `toy_panel.py`** | the owner's two-exon transcript and the instrument that sweeps it over the whole cache |
| ⛔ **three bench errors found and fixed** | `TRAPS.md` B15/B16/B17 — capture depletion ≠ starvation, split probes, substitution ceilings |
| ⚠ **C2 promoted to a solver defect** | a boolean licence is flipped at full strength by κ̂ = 0.500689. `TRAPS.md` D4f. It is `ROADMAP.md` step 1 and it caps what step 3 can deliver |
