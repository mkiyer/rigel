# Session handoff — 2026-08-12

    ⚠ **A DEV DOC.** Provisional; nothing may cite it. When an item settles, MOVE it to its permanent
    home and delete it here in the same edit.

    ⭐ Deliberately SHORT. The previous session's handoff and a 675-line design doc were deleted when
    their content moved to `EQUATIONS.md` §9b, `TRAPS.md` and `ROADMAP.md` §0. ⛔ If this file starts
    to grow, that is the signal to move things out, not to add a section.

---

## 0. WHERE THE STATE LIVES — read these, not this

| question | answer |
|---|---|
| where is the tool? | `ROADMAP.md` §0 |
| what do I do next? | ⭐ **§1 below**, then `ROADMAP.md` §2 |
| how does the RNA prior work today? | `EQUATIONS.md` §9b |
| what must I not repeat? | `TRAPS.md` |
| how do I run anything? | `CLAUDE.md`; the panel loop is `scripts/sim/panel.py` |

**Suite: 3,299 passed / 0 skipped / 10 xfail / 0 failed.** `ruff` clean. C++ built. Tree committed.

---

## 1. ⭐⭐⭐ THE TASK: ALLOCATE THE RNA PRIOR PER TRANSCRIPT

⛔⛔ **`docs/dev/rna_prior_alloc.md` IS THE OWNER'S OWN WRITE-UP AND IT IS THE SPEC. READ IT FIRST AND
IN FULL.** It carries the theorem, the harmonic-mean proposal, the zero-handling options and the KEY
REALIZATION about unspliced mass. This section does not restate it — it records only what was already
verified against the code, so the next session does not re-derive it.

### 1.1 What is verified TRUE about today's prior

* **The RNA prior is UNSPLICED fragments only.** `priors.py` subtracts `mass_rna_spliced_edge`, and the
  junction flux is deliberately never added: *"a locus whose RNA is fully spliced SHOULD get a
  near-zero `rna_prior_count`"* (owner ruling, 2026-07-30). ⭐ **The owner's geometry argument is
  therefore already accepted policy at LOCUS level** — this work extends it one granularity down.
* **It was never applied uniformly.** It is shared in proportion to each component's CURRENT COUNT — a
  common multiplicative factor over the eligible components, so before this session it was provably
  neutral on the within-RNA split. `EQUATIONS.md` §9b has the arithmetic.
* **The warm start is `unambig_totals` + coverage-weighted shares, then projected through the SAME
  prior update** the M-step uses. So "the prior multiplies the warm start" is correct, and there is one
  projection with two call sites.
* **Strength is exactly 1 pseudo-fragment per real unspliced fragment, with no knob.**
* ⚠ But the EM's `rna_count` it competes against ALSO contains spliced units, so the prior is weakest
  at spliced-heavy loci and **strongest exactly at the unspliced-dominated loci** where a wrong
  allocation does the most damage.

### 1.2 What already exists and must not be rebuilt

* ⭐⭐ **`_transcript_node_incidence` (`capture_eff_length.py`) is the transcript → node/edge/junction
  map the write-up asks for**, and returns exactly its three requirements: regions, interior boundaries
  (contiguously crossed edges), and splice junctions with both flanks. Annotation-only,
  **0.36 s** on the suite index. Do not write a second one.
* ⭐ **`rna_node_eff_len` / `rna_edge_eff_len` on `CalibrationResult` are the RNA population's OWN
  opportunity per object** — the divisor a density argument needs. They have NO consumer in `src/`
  today, and that is because the prior was rewritten as a conserved fragment count that divides by
  nothing (`d045d820`, `5591cc01`). ⛔ **They are ORPHANED, not dead** — audited 2026-08-12, three
  stale docstrings claiming a vanished consumer were corrected. They are PROJECTED off
  `NodeGeometry.eff_rna`, never recomputed, so they are byte-identically what the SOLVER divided by.
* `transcript_capture_eff_lengths` already performs "the same operation over a transcript's own
  objects", including the capture contraction. Read it before designing the effective-length half.

### 1.3 ⛔ THE FOUR THINGS THAT WILL BITE

1. ⛔⛔ **A GEOMETRY-PROPORTIONAL ALLOCATION IS AN IMPLICIT UNIFORM-ABUNDANCE PRIOR AT THE STRENGTH OF
   THE DATA, AND IT MOVES THE FIXED POINT.** `apply_grouped_prior_update` runs inside EVERY M-step, not
   just the warm start. Allocating in proportion to unspliced opportunity asserts equal ABUNDANCE
   across a locus's transcripts and would shrink isoform ratios toward equal usage — worst at
   unspliced-dominated loci, and it would damp the paralog vertex `TRAPS: identical-paralogs-are-bimodal`
   says to localise rather than smooth. **Any allocation must be checked against this.**
2. ⛔ **INTRONIC MASS HAS NO HOME.** A `Locus` is the merged span of its ACTIVE transcripts, so
   `assemble_priors` already sums intronic nodes — and the RNA prior is specifically unspliced RNA,
   which lives disproportionately in introns. But `_transcript_node_incidence` puts a multi-exon
   transcript on its EXON nodes only; the sole entity covering an intronic node is the synthetic
   nascent shadow, which currently receives zero prior. **Drop that mass and the locus total breaks
   (moving the gDNA:RNA split); keep it and only the nascent shadow can hold it.** This is the crux.
3. ⚠ **VBEM HAS AN UNDECLARED ACTIVATION THRESHOLD.** With per-row max subtraction, `EXP_CUTOFF = -708`
   and `digamma(a) ≈ −1/a`, a per-transcript pseudocount must exceed roughly **1/708 ≈ 0.0014** in
   alpha units before its component receives ANY responsibility. A per-transcript scheme inherits that
   as a magic number unless it is handled explicitly. VBEM is the shipped default (`mode: "vbem"`).
4. ⚠ **`_region_locus_shares` drops the nodes of UNEXPRESSED genes** (sample-dependent), so "the prior
   is a conserved fragment count of the library" does not survive contact with real data. Scope any
   conservation claim to the loci that exist.

### 1.4 Numbers already measured — do not re-measure

| | |
|---|---|
| nodes contained by exactly ONE transcript | **6.9 %** (median 6 transcripts/node, max 391, suite index) |
| `P(fragment unspliced)` per annotated transcript @212 bp | p10 **0.000**, median 0.463, p90 **0.953** |
| annotated transcripts that should get ~ZERO unspliced prior | **16.7 %** |
| annotated transcripts that should take ~ALL of it | **10.1 %** |
| `P(fragment unspliced)` for SYNTHETIC nascent entities | **1.000 by construction** |

⭐ The first row is why the SPLIT rule does nearly all the work. ⛔ The last row is why a pure-geometry
split is dangerous: a nascent entity has maximum unspliced opportunity and typically minimum abundance.

### 1.5 Plumbing and gates

* A per-LOCUS array rides `array[ids]` (`pipeline.py`); a per-TRANSCRIPT array rides the
  `t_eff_lens` / `t_is_synthetic` lane (flat `float64[n_transcripts]`, remapped in `em_solver.cpp`).
  **Mixing the two conventions is the easy plumbing error.**
* ⛔ **The conservation identity has NO TEST and cannot have one from Python** —
  `apply_grouped_prior_update` is `static`. Exposing it behind a test-only binding is a PREREQUISITE:
  a per-transcript vector makes that identity harder to hold, not easier. See the comment block at the
  end of `tests/test_batch_em_impl.py`.
* ⛔ `gdna_eff_len` is NOT a pseudocount — it is base pairs entering as the gDNA component's log
  effective length. An arm that injects "the prior" and means all three `LocusPriors` fields is
  measuring two mechanisms at once.
* No instrument scores a PER-TRANSCRIPT prior; `prior_vs_oracle.py` is per-locus. **Build the scoring
  arm before the mechanism** (`TRAPS: measure-the-ceiling-first`).

---

## 2. WHAT THIS SESSION CHANGED

| | |
|---|---|
| ⭐ **synthetic nascent entities get ZERO RNA prior** | in scope: nascent **6,238,406 → 691,796**, gene Σ\|err\| **0.608**, `g00` gene **0.476**. `EQUATIONS.md` §9b |
| the end-to-end baseline | re-derived on all 36 ladder conditions, BOTH halves. `ROADMAP.md` §0 |
| ⭐ `scripts/sim/panel.py` | the sim + benchmark workflow, one command per stage. `TESTING.md` §2 |
| a per-script IMPORT GATE | a green suite had hidden **five dead instruments**; 4 repaired |
| deleted | `scripts/benchmarking/` (4,036 lines) and `rigel/sim/analysis.py` (1,589 → 618 kept as `net_flow.py`) |

---

## 3. ⛔⛔ PROCESS — THE SAME CLASS OF ERROR TWICE MORE

⭐ **CHECK `git status` BEFORE BELIEVING ANY REPORTED STATE, INCLUDING YOUR OWN.** Last session a
subagent edited production calibration code during an "interrogation". This session **a workflow agent
wrote an unrequested file into `docs/dev/`** — invisible except as a +1 in the parametrised doc-test
count, which is the only reason it was caught. ⛔ Re-derive the test count after any change; never
adjust it.

⚠ **And three separate times a PER-M-STEP algebraic identity was stated as an END-TO-END invariant**
("the library gDNA fraction cannot move"). It can, and it did (+4.21 M). One of those became a test
whose premise was false and which failed on the shipped configuration. **When a quantity is conserved by
one function, say which function — the EM iterates around it.**

⚠ A subagent reported `rna_node_eff_len` / `rna_edge_eff_len` as "dead surface"; that was relayed
without checking and the owner caught it. They are orphaned by a design change and are **exactly the
divisor the next task needs**. Verify a subagent's negative before repeating it.
