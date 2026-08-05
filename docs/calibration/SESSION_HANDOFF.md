# SESSION HANDOFF — ⭐⭐⭐ ψ CANNOT SOLVE TO A SIMPLEX VERTEX, AND ONE VERTEX IS MOST OF THE ANNOTATION

    Written 2026-08-05. ⚠ **WORKING DOC, NOT A PERMANENT FIXTURE.** The six permanent docs are
      `CLAUDE.md`, `docs/SUCCESS.md`, `docs/ROADMAP.md`, `docs/TRAPS.md`, `docs/EQUATIONS.md`,
      `docs/DESIGN.md`, `docs/TESTING.md`.
    ⛔ **DELETE THIS FILE when its task lands**, promoting anything worth keeping into `ROADMAP.md`
      (a current number), `TRAPS.md` (a lesson), `EQUATIONS.md` (a derivation) or `DESIGN.md` (a ruling).

## Read in this order

| | |
|---|---|
| 1 | `CLAUDE.md` — doc map, working rules, the script table |
| 2 | ⭐⭐⭐ **`docs/ROADMAP.md` §2** — the task, measured. §1 is the ranked list, §3 is what must NOT be rebuilt |
| 3 | ⭐⭐ **`docs/TRAPS.md` D5** (a grid supplies a prior whether you ask or not), **D1** (a variance cannot fix a bias), **A14** (a "no change" arm is not a control until you count the opportunities), **A12** (never the honesty metrics alone) |
| 4 | ⭐ **`docs/EQUATIONS.md` §9** (priors on a grid) and **§3.5d** (why `r_tot` is not a component's ratio — and why a perfect `r_g` is worth zero here) |
| 5 | `docs/TESTING.md` §0b — the toy harness and its spec ladder |

---

## 1. State of the tree

**Branch `fragment-length-gold-standard`.** Everything from the 2026-08-05 sittings is committed. `src/`
carries the splice-flux reframe as the **default, ungated** path (owner's ruling — it is theoretically
correct; the old junction-inclusive total was catastrophically wrong).

**Test baseline: 22 failed / 2214 passed / 3 xfailed.** The 22 are the standing set — 21
`test_golden_output` + the paralog row (`TRAPS.md` D9). ⭐ The third xfail is new and is
`test_toy_harness::test_the_harness_REPRODUCES_the_intron_composition_dependence`, `strict=True`: it is the
project's detector for the level defect that the reframe's correction un-masked, and **it must go green
again when step 2's joint arm lands**. ⛔ Its 2.0 bound was deliberately NOT widened.
`ruff check src/ tests/ scripts/` clean. ⚠ **Never `ruff format scripts/`.**
⚠ The oracle cache at `~/Downloads/rigel_runs/suite/ladder/oracle_cache` is **VALID — do not delete.**

---

## 2. ⭐ WARM-UP: FINISH `alt_splice` — 20 minutes, and it closes an owner question

The owner asked whether alternative splicing solves correctly and the rung exists but is **UNVERIFIED**. No
solver number may be read off it until the substrate is verified (`TRAPS.md` A1).

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
export OMP_NUM_THREADS=1

python scripts/design/verify_toy_substrate.py --spec alt_splice --no-gdna --n-rna 200000
python scripts/design/verify_toy_substrate.py --spec alt_splice --n-rna 200000 --nrna 60
for p in pos_blocks transcript_order single_geom drop_junction; do
    python scripts/design/verify_toy_substrate.py --spec alt_splice --n-rna 200000 --perturb $p; done

python scripts/design/reframe_walk.py --arms high --n-rna 200000 --spec alt_splice
python scripts/design/zero_controls.py --specs alt_splice
```

**`alt_splice`** is TA+ (1,000–2,000)(5,000–6,000)(9,000–10,000) and TB+ (1,000–2,000)(9,000–10,000) —
inclusion and skipping isoforms. Three junctions over two shared sites: 2,000→5,000, 6,000→9,000 and the
skipping jump 2,000→9,000.

| what to check | expected, from `EQUATIONS.md` §3.6c |
|---|---|
| EDGE @2,000 | genomic-LOW end of **TWO** junctions ⇒ both fluxes in `junction_count_lo`, `_hi` = 0, pooled `Σcount/ΣE` |
| EDGE @9,000 | genomic-HIGH end of **TWO** ⇒ both in `_hi` |
| EDGE @5,000 / @6,000 | one junction each, `_hi` / `_lo` respectively |
| node [5,000, 6,000) | ⭐ exon AND intron **on the same strand** — TA's middle exon inside TB's intron. No strand bit separates them and `coarse_type_array` calls it `exon` |

⚠ **The prediction I would most like falsified:** the skipping junction (2,000→9,000) passes *over* node
[5,000,6,000), which is TA's exon. So that node holds TA's mature RNA contained in it **and** TB's mature
passing over it via a junction it does not attach to. The flank rule says nothing about that population, and
it is where I expect a surprise.

⚠ What `alt_splice` does **not** cover: an EDGE that is one junction's LOW end and another's HIGH end at
once. That needs one transcript's intron to END where another's BEGINS; it is gated synthetically in
`tests/calibration/test_splice_flux_reframe.py::test_ONE_EDGE_can_be_the_LOW_end_of_one_junction_and_the_HIGH_end_of_another`.

---

## 3. ⭐⭐⭐ THE TASK — ψ CANNOT REACH A VERTEX

`ROADMAP.md` §2 has the full measurement. The one-paragraph version:

A silent gene — every transcript `abundance = 0`, so every object is pure gDNA and the truth is `f_g = 1.000`
exactly — solves to **0.92 – 0.99**, on every rung, for 0.65–2.51 % of mass. The objects where ψ is BYPASSED
(structurally pure-gDNA G1 objects) are **exact**. The objects where ψ SOLVES land 5–8 % short of the vertex.
It is not the messages: what they deliver already implies `f_g` = 0.886–1.148, it does not shrink over a 40×
count sweep, and it gets **worse** as the log-odds window widens.

⭐ **Why this is step 1 and not the message layer.** Most annotated transcripts are unexpressed, so this is
the modal object in a real library, not a corner. And because the zero-RNA arm has a *perfect* `r_g` and
*correct* delivered messages and still lands 5–8 % short, the vertex residual is a **floor under every
message-layer arm** — which is very likely why panel A/Bs keep returning ±0.5 %.

### ⛔ STEP 1 — FIRST MEASUREMENT, BEFORE ANY DESIGN: is the phantom-gDNA floor the SAME bug?

```bash
python scripts/design/zero_controls.py --specs silent TA_single_exon spliced_exons nested_exons
```

Both arms, one command. The zero-**gDNA** arm reads `f_g` = 0.0218 against exactly 0.0000; the zero-**RNA**
arm reads 0.946 against 1.000. Both are ~2–7 %, both sit at a vertex of the composition simplex
(`f_g = 0` and `f_g = 1`). ⭐ **If one mechanism covers both vertices, one fix closes two of the owner's
stated concerns and the value of step 1 doubles. If not, they are separate builds.** Establish which
before designing anything — this is the ceiling discipline applied to a mechanism rather than a channel.

The falsification handle is already measured and works: sweep `sweep_logodds_window` and watch. If the
zero-gDNA arm's phantom moves the same way the zero-RNA arm's does, it is one bug.

### ⚠ What is NOT yet known, and must not be assumed

* **Whether the fix is the grid, the reference prior, or the parameterisation.** Three candidates, and the
  measurement above does not choose between them:
  1. the log-odds lattice cannot represent a vertex (`log(f_R/f_g) → −∞`), so no posterior median reaches it;
  2. the reference prior on that lattice is effectively Haldane-ish and amplifies toward the vertices
     (`TRAPS.md` D5 measured a 0.045-vs-0.525 spread between Jeffreys and Haldane) — ⚠ note the *sign* here
     is the puzzle: widening made it worse, which a pure "cannot reach the vertex" story does not predict;
  3. ψ's uninformative reference is being fused as if it were evidence when the object has none.
* ⛔ **No magic numbers.** A "snap to the vertex when close" rule is exactly the tuned constant this project
  has refused three times (`TRAPS.md` B11, D4f, D4g). The fix must be a derivation.
* ⛔ **Measure the CEILING first** (`TRAPS.md` B1 — it has re-ranked this project five times). Hand every
  evidence-free pure-gDNA object the exact answer and re-solve: what does perfecting the vertex buy on the
  36-condition ladder? That number decides whether this is a build or a note.
* ⚠ **Then the panel arm before `src/`** (`TRAPS.md` B18 — three toy-positive changes have been
  panel-negative). `ladder_arm_ab.py`, ~40 s/condition with the oracle cache; shard with `--conditions` and
  say which rows carried it.

---

## 4. ⭐⭐ THE OTHER LIVE THREAD — the joint arm (`ROADMAP.md` §1 step 2)

Do not start this before step 1 unless step 1's ceiling comes back ~0. Both halves already exist:

* **half A** is in `src/` and is the default — the splice-flux reframe (`EQUATIONS.md` §3.6c);
* **half B** is `toy_trace_error.py`'s relay-level arm — withholding the composition licence on the hops out
  of the intron, measured at **34 %** of that gene's error on its own.

⛔ Alone, half A is panel-negative on the node axis (+0.5 % mwae, +36.9 % confidently-wrong) and positive on
the edge axis (−0.3 %, −7.7 %, 6/6 better on the shipped solve). `TRAPS.md` D4j says that is not a verdict on
it: it is half of a **cancelling pair**, and the strict xfail in §1 is the detector. ⭐ The joint arm is the
only measurement that prices either.

---

## 5. ⭐ THE INSTRUMENTS BUILT THIS CAMPAIGN, and what each is for

All are in `CLAUDE.md`'s script table with their verdicts in their own docstrings.

| | |
|---|---|
| `zero_controls.py` | ⭐⭐⭐ **the two zero arms — the owner requires them on every experiment.** How step 1 was found |
| `reframe_walk.py` | ⭐⭐⭐ every count, both flank totals, every hop in both directions, with `r` as used / as the predecessor gave it / as TRUTH says. §3b decomposes `r_tot` into its per-component ratios |
| `verify_toy_substrate.py` | ⭐⭐⭐ is the INPUT correct — every bank re-derived from per-fragment truth. Any number of transcripts, either strand. ⛔ Run before reading a solver number off a new rung |
| `toy_trace_error.py` · `toy_dissect.py` · `toy_ceiling.py` · `ladder_arm_ab.py` · `length_ceiling.py` · `verify_capture.py` | the error trace, the channel dissection, the re-solve ceiling, the panel arm, the length ceiling, the capture A/B |

---

## 6. ⛔ WHAT WAS CLEANED UP, so the next session is not confused by its absence

* **`docs/calibration/ISSUES.md` — DELETED.** Its register had become a second copy of `ROADMAP.md` §1/§3
  (`TRAPS.md` A11's shape: two homes for one thing). Every live item was folded into `ROADMAP.md` §1, every
  closed item into §3 as one line, and every `ISSUES.md C-number` citation in `docs/`, `scripts/` and
  `tests/` was repointed. There are zero dangling references.
* **`docs/calibration/solver_derivation.md` — DELETED.** Superseded by `EQUATIONS.md` §3.5d/§3.6c.
* **`docs/calibration/splice_flux_reframe_src.patch` — DELETED.** ⛔ It was a "revert convenience" copy of a
  `git diff`, i.e. exactly the code-kept-for-comparison defect `TRAPS.md` G3 bans. Git is the revert
  mechanism.
* **`ROADMAP.md` 1,030 → ~175 lines**, and it now points forward: state, the ranked list, one line per
  closed item, the open questions, the rules. ⛔ It had grown a section per campaign, which is a changelog.
* **`TRAPS.md` 741 → ~514 lines with every label intact.** The `D4` family's ten entries are now one lesson
  with a ten-row instance table, and the longest B/A entries are compressed to the mistake–tell–rule form.
  ⛔ **Labels are cited from source and tests, so a label is never renumbered or deleted** — a merged label
  survives as a table row. The file's header now carries that rule and the pruning discipline.
