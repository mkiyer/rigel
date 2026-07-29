# Accumulator v5 / P1g — SESSION HANDOFF, 2026-07-30

**Read this file, then `accumulator_ledger.md` (bottom four entries), then `06_implementation_plan.md`.**
Everything below was measured or run in the session that wrote it. Nothing is quoted from an older
document without being re-verified — where an older number turned out to be wrong, the correction and
its measurement are given.

The previous handoff (`SESSION_2026_07_29_HANDOFF.md`) is **superseded by this one**. It is kept only
because its §1 (the design in one page) and §4 (things not to re-derive) are still accurate.

---

## 0. STATE IN TEN LINES

| | |
|---|---|
| branch / HEAD | `main` @ `3c293038`. ⛔ **Nothing committed — the owner drives commits.** ~562 modified/untracked paths. |
| test suite | **1280 pass**, ruff clean on `src/ tests/ scripts/`, `ruff format` parity, warning-clean under `-W error::RuntimeWarning` |
| baseline (32-cond `ambig_dense_10mb`, mass-weighted mwae) | **r0 `0.079005` / r3 `0.046675`** — v8 partition, reproduced 32/32 exactly at the end of every arm below |
| phases done | **W0** · **W0b** · **W1a** · **W1b** · **W1b-clean** · **W1c** · **W1-review** |
| next phase | ⭐ **W2a (C3)** — the structural `accept_l`/`accept_r`. Fully specified in §4 below; it needs no derivation, only implementation. |
| indexes | **8/8 at INDEX_FORMAT_VERSION 8.** `regions.feather`/`boundaries.feather` deleted from all of them. |
| C++ | **not rebuilt this session and not touched.** No `src/rigel/native/` edit is pending. |
| goldens | **did NOT move** at W1b (measured: zero mergeable adjacencies in all 20 golden scenarios). One regeneration still expected, at W7. |
| cfRNA fast loop | **restored on v8** — 3 caches rebuilt and `partition_hash`-verified. |
| scratchpad | **338 → 33 files.** Survivors all run against the current API. |

---

## 1. WHAT TO READ, IN ORDER

| # | doc | why |
|---|---|---|
| 1 | **this file** | state, the remaining implementation, the open questions |
| 2 | `accumulator_ledger.md` — the **last four** entries (W1b, W1b-clean, W1c, W1-REVIEW) | every measured delta, and the two real bugs found |
| 3 | `../calibration/P1G_SCOPE.md` §7 (C1/C2/C3) and §8 (the A/B protocol) | the W2 contract. ⚠ §9's biggest stated risk is now **RETIRED** — see §3 item 6 |
| 4 | `06_implementation_plan.md` | ⚠ read its **amendment banner first**; F1's two-filter rule and W4's shadow mode are both struck |
| 5 | `05_accumulator_v5.md` | THE SPEC for W3–W7. Still current. Read through the plan's F-findings |
| 6 | `../index/00_splice_graph_design.md` | the graph. Its amendment banner records which gates died |

⛔ **Do NOT read as guidance:** `02_redesign_derivation.md`, `03_path_accumulator.md`,
`04_accumulator_v3.md`, `../calibration/PATH_MARGINALIZATION.md`, `../calibration/archive/**`.
`dag_design.md` is the owner's concept draft — context, not spec.

---

## 2. WHAT LANDED THIS SESSION (four arms, all bit-identical on the bench)

| arm | what it did | gate |
|---|---|---|
| **W1b** | `build_node_partition_arrays` feeds the scanner the v8 node cut array; `RegionArrays.from_index` follows it | 32/32 exact |
| **W1b-clean** | **deleted the v7 partition entirely** — `calibration/regions.py` whole, both feathers, `index.region_df`/`.boundary_df`, `coarsen_nodes_to_regions`. P3 → **I3b** | 32/32 exact |
| **W1c** | the 8 structural flags reach `NodeStatics.boundary_flags` (raw `uint16`). **Nothing reads them yet — W2a is the first consumer** | 32/32 exact |
| **W1-review** | 5 adversarial lenses / 64 agents / 59 findings / **42 confirmed and fixed** | 32/32 exact |

**Net `src/` change across the four: −849 lines.**

### The two real bugs, both found by tests rather than by reasoning

1. ⛔ **The flags filter deleted 26,475 real transcripts' termini.** Plan F1 specified
   `~is_synthetic & ~is_nrna` for flags/reaches. On a **non-synthetic** row `is_nrna` means *this real
   transcript is single-exon, so mature ≡ nascent* — **not** "manufactured span". Measured: all 26,475
   such human rows have `n_exons == 1` and none is a `RIGEL_NRNA_*` row, while every manufactured span
   is already `is_synthetic`. **52,104 distinct terminus positions** were being dropped. Now **ONE
   filter, `~is_synthetic`** (`splice_graph._is_real`).
2. ⛔ **I13 was supplying false assurance.** It validated the flags by calling *the builder's own event
   emitter*, so it checked placement and never the event definition. Proven by mutation: delete the
   NEG-strand TSS/TES swap and the **entire 1,289-test suite passed and `validate_graph` accepted the
   graph**. Now `_events_independently` re-derives all eight event sets by a different algorithm
   (per-transcript `min`/`max` vs the builder's cumulative-exonic-offset arithmetic), and it rejects
   the swap. A second hole closed with it: bits whose event class is empty on a reference were never
   compared at all.

---

## 3. THE MEASURED FACTS — DO NOT RE-DERIVE ANY OF THESE

1. ⭐ **The 32-condition bench CANNOT see a partition change.** `ambig_dense_10mb`'s v8 `nodes_df` is
   **row-for-row identical** to its old merged `region_df` (1,698 == 1,698 on
   `ref_name/start/end/length/signature`). That suite has no alternative TSS/TES. Its
   `partition_hash` did not change, so its 32 cached payloads survived W1b untouched.
2. ⭐ **THE PARTITION EFFECT, measured on `quick_3to1_5mb`** (16 oracle conditions, +25.0 % nodes):
   **+0.0069 at refit=0, −0.0014 at refit=3.** It splits structurally — with gDNA present v8 **wins**
   (up to −0.0152); with **zero** gDNA v8 **loses** (up to +0.0836), because thinner nodes drift off
   `f_g = 0` and there is nothing to win. **Pass-0 pays; the refit recovers it.**
   ⚠ The Jensen confound was bounded first and is **0.998× on that suite** (vs 0.739× on human), so
   the sign is interpretable. **Do not carry 0.998 to another suite.**
3. ⭐ **Scan cost of the finer partition, real cfRNA, +38.7 % nodes:** LBX0190 **+9.7 %**, MO_3021
   **+6.9 %**. This is graph-doc §8's required acceptance number; it is measurable and small.
4. ⭐ **53.4 %, not 59.5 %**, of real human transcript termini are invisible to a merged partition
   (232,451 of 435,291). The 59.5 % figure was itself computed under the buggy `~is_nrna` filter.
5. ⭐ **Positions that are BOTH terminus and splice site: 0.99 %** of human contiguous edges (10,337 of
   1,043,595), **not "the majority"** as the graph doc claimed. Terminus-only 40.70 %, splice-only
   58.31 %, **neither 0** — GENCODE contains no bookended exons.
6. ⭐⭐ **P1G_SCOPE §9's biggest stated risk is RETIRED.** It warned: *"an alternative TSS inside an
   exon does not even produce a region boundary — so a share of real termini will be invisible to a
   boundary-keyed bit no matter how it is built"*. **Under v8 every exon endpoint is a cut**, so every
   real terminus is a boundary and carries its flag — guaranteed by I13, not merely by construction.
   That is precisely what the partition change bought. **P1G_SCOPE's step 3 ("measure the bit on the
   real cfRNA index — the one that could invalidate the plan") is CLOSED**; the census is item 7.
7. ⭐ **Boundary-slot census** (`build_boundary_flags_array`, all slots incl. terminals):

   | index | slots | terminus-only | splice-only | both | neither |
   |---|---|---|---|---|---|
   | `ambig_dense_10mb` | 1,699 | 40.7 % | 47.4 % | **11.8 %** | 0.1 % |
   | `quick_3to1_5mb` | 1,527 | 33.7 % | 66.2 % | **0.0 %** | 0.1 % |
   | **human (cfRNA)** | 1,044,167 | 40.7 % | 58.3 % | **1.0 %** | 0.1 % |

   ⚠ **The toy over-represents the both-bits case 12× and `quick_3to1_5mb` cannot exercise it at
   all.** Any W2 branch keyed on "both" is validated by the toy against a population 12× too large.
   "neither" is the two per-reference terminals.
8. **`partition_hash` covers `nodes.feather` ONLY** — deliberately, because it is a cache key for a
   *scan* and the scan sees the cut array. `edges.feather` can change without it moving, **and did**:
   the flag fix rewrote every edge file and left every node file byte-identical. **Anything caching an
   edge-derived artifact must carry its own provenance.** `calib_cache` therefore derives
   `boundary_flags` from the index rather than caching them.
9. **I3b, I4, I11 and I13 run at BUILD only** — they are gated on `transcripts`, which `load()` never
   passes (reconstructing them costs ~3 s at human scale). A `signature` or `flags` column that has
   drifted from the annotation **loads clean**. Load-time validation is the graph-internal half
   (I1/I2/I5–I9/I12) and costs **0.23 s**.

---

## 4. ⭐⭐ THE NEXT ARM — W2a (C3): STRUCTURAL `accept_l` / `accept_r`

**This section is written so it can be implemented without deriving anything.**

### 4.1 What `accept_l` / `accept_r` mean today

`src/rigel/calibration/bp_solver.py:383-386`:

```python
accept_l = (SP[0] + SN[0]) > _EPS   # the LEFT face carries the spliced (acceptor) => WITH-spliced rho_tot
accept_r = (SP[1] + SN[1]) > _EPS
```

where `SP = (geometry.spliced_pos_left, geometry.spliced_pos_right)` and `SN` the neg twin
(`bp_solver.py:186-189`). They are **observational**: *"spliced mass was observed on this face"*.

Their only consumer is `_rho_faces` at `bp_solver.py:537-546`:

```python
def _rho_faces(fgc):
    ru, rw = node_total_density(chain, geometry, fgc)
    rs = rw - ru                                   # the one-sided spliced density
    return (ru,
            ru + np.where(accept_l, rs, 0.0),      # left face
            ru + np.where(accept_r, rs, 0.0))      # right face
```

`_rho_faces` is the **weld** that makes C1/C2/C3 one function (plan F9): it produces `rho_l0`/`rho_r0`
(`:549`), which feed `_seam_pair` (`:1042`, C2's fit population), `_relay` (`:1046-1047`) and
`_transport` (`:1330-1333`) — the reframe C1 gates — and `_flank_dom` (`:993`, C2's fit).

### 4.2 ⭐ The structural predicate — DERIVED AND MEASURED, just implement it

A boundary's LEFT face lies inside the left node. The `DONOR_s` flag means *"this position is the
`intron_start` of a strand-`s` intron"* — the **genomically low** end — so the exon abutting it is on
the **left**, and a spliced fragment's exonic block lands on the **left face**. `ACCEPTOR_s`
(`intron_end`, genomically high) puts the exon on the **right**. Therefore:

```python
accept_l = is_splice_donor(flags)      # DONOR_pos | DONOR_neg     -> exon on the LEFT
accept_r = is_splice_acceptor(flags)   # ACCEPTOR_pos | ACCEPTOR_neg -> exon on the RIGHT
```

⚠ **This is simpler than plan F10's phrasing and supersedes it.** F10 correctly refuted the graph
doc's `rna_crosses_contiguously = exon_s(src) AND exon_s(dst) AND NOT terminus` (only 13.7 % of
junction seams have `exon_s` on both flanks, so a literal swap drops the mature term at ~86 % of
them). F10 then proposed *"`DONOR_s or ACCEPTOR_s` on the flank carrying `exon_s`"*. The `exon_s`
qualifier is unnecessary: `DONOR` already implies an exon on the left and `ACCEPTOR` an exon on the
right, and invariant I9 asserts that both endpoints of every junction carry that strand's exon bit.
**Use the two-line form above.**

⚠ **`splice_graph.is_splice_site` ORs DONOR and ACCEPTOR together — do not use it here.** W2a needs
them **separately**, per face. Add `is_splice_donor` / `is_splice_acceptor` beside it (three lines
each, mirroring `is_terminus`), or read the bits inline.

### 4.3 ⭐ MEASURED, on the real cached payloads — this is W2a's prize and its risk

Observational vs structural on `ambig_dense_10mb`, per condition (boundary nodes only):

| condition | obs `accept_l` | **struct** | agree | obs `accept_r` | **struct** | agree |
|---|---|---|---|---|---|---|
| `gdna100 ss0.50 capture_off` | 526 | **503** | 95.8 % | 528 | **503** | 95.7 % |
| `gdna100 ss0.50 capture_on` | 392 | **503** | 89.5 % | 393 | **503** | 89.3 % |
| `gdna100 ss0.50 capture_verystrong` | **263** | **503** | 82.9 % | **264** | **503** | 82.9 % |
| `gdna100 ss0.99 capture_off` | 524 | **503** | 95.6 % | 527 | **503** | 95.6 % |

**Read this carefully — it contains both directions of the change:**

* ⭐ **The structural count is CONSTANT at 503.** It is a property of the annotation. That is the whole
  point of C3: coverage-independence.
* ⭐ **The observational count collapses with coverage** — 526 → 392 → **263**. At capture-verystrong
  the observational predicate is blind to **48 %** of the structurally-real spliced faces, and it is
  blindest exactly where coverage is thinnest, which is the wrong direction.
* ⚠ **At capture-OFF the observational count EXCEEDS the structural one** (526 vs 503). Those ~23
  faces carry observed spliced mass at a position with **no annotated junction** — unannotated
  junctions / `SPLICED_IMPLICIT`. **C3 will REMOVE the spliced term from `rho_tot` at those faces.**
  That is a real behaviour change in the opposite direction to the main effect and it must be in the
  pre-registration.

### 4.4 Implementation steps for W2a

1. **Re-record the baseline** from the current tree, both refits, into `/tmp/w0_baseline_{r0,r3}.tsv`
   as arm `w2a_pre`, and confirm 32/32 exact against arm `w1c_rev`. (Standing rule; every stored
   number goes stale the moment the tree moves.)
2. Add `is_splice_donor` / `is_splice_acceptor` to `splice_graph.py` beside `is_terminus`, and export
   them. Keep `is_splice_site` (it is the OR, and C4 may want it).
3. Thread the flags to the point of use. `NodeStatics.boundary_flags` is already on the chain
   (`uint16`, 0 on region nodes and terminals) — `bp_solver` receives `statics` and reads
   `statics.u_pos` etc. at `:208`, so the array is already in scope at `:383`.
4. Replace `bp_solver.py:383-386` with the structural form. ⚠ **Keep the observational arrays
   available under a different name** if a diagnostic wants them, but do **not** leave a fallback: a
   silent "structural if available else observational" is exactly the class of defect the review
   removed everywhere else.
5. ⚠ **`node_total_density` returns the two per-face spliced densities SUMMED**
   (`node_geometry.py:370-372`). That is safe today only because `accept_*` is observational, so a
   face with no spliced mass contributes zero anyway. **With a structural predicate that is no longer
   true**: a face can be structurally accepting while carrying no observed mass, and `rs` would then
   deliver the *other* face's density to it. **Split the return into `(rho_unspliced, rho_spliced_left,
   rho_spliced_right)` and have `_rho_faces` use the matching face.** This is not optional and it is
   the part most likely to be got wrong.
6. Run the gates in §6.

### 4.5 ⭐ PRE-REGISTRATION for W2a — predictions and the falsification test

| # | prediction | falsified by |
|---|---|---|
| P1 | **capture-ON and capture-verystrong improve most** — that is where the observational predicate loses 22–48 % of its faces | flat or worse under capture |
| P2 | **capture-OFF moves least**, and may move *slightly negative* from the ~23 unannotated-junction faces losing their spliced term | a large capture-OFF move in either direction (means the predicate is not doing what §4.2 says) |
| P3 | the number of accepting faces becomes **identical across all 32 conditions** (503 on this index) | any per-condition variation — the predicate is then still reading data |
| P4 | ⭐ **FALSIFICATION: a boundary with NO annotated junction on either side must have `accept_l == accept_r == False` in every condition.** Under the observational rule ~23 of them are True at capture-OFF | any such boundary still accepting ⇒ the flags are not being read, or the wrong bit is |
| P5 | `gdna_none` conditions unmoved (no gDNA to redistribute) | a material move there |

⚠ **Write a test that FAILS against the pre-W2a tree first** — e.g. a fixture where a face carries
observed spliced mass at a position with no annotated junction, asserting the structural predicate
excludes it. That test is the arm's own proof it did something.

### 4.6 Then W2b (C1) and W2c (C2) — one arm each, baseline re-recorded between

* **W2b — C1, the reframe at terminus seams.** The prize: capture-OFF **−0.0178 / −0.0124**
  (P1G_SCOPE §7). ⛔ **Do NOT divide `r` by `f_g(dst)`** — the closed form says that is exactly the
  missing factor and it is the destination's own belief, i.e. the pin bug rebuilt. The narrowest fix
  consistent with the architecture: **a pure-gDNA source's density level must not be reframed by a
  total-density ratio.** Gate the **gDNA leg only** of the reframe (`_relay` `:1046-1047`,
  `_transport` `:1330-1333`) on `is_terminus(flags)`.
  ⭐ **Falsification: junction-only edges must be UNMOVED** — they are already exact at 1.0×, so any
  movement there means the bit is applied too widely.
  ⚠ **capture-ON must not regress** — the `r_g ≡ 1` sizing arm failed exactly there, at **+0.1728**.
* **W2c — C2, `ω_graft` per structural class.** `ω̂` is 1.7–1.9 at termini vs 0.04–0.06 at
  junction-only seams. The partial-pooling block in
  `calibration/enrichment_frame.graft_premise_logvar` is the plug-in point.
  ⚠ **Watch the fit-population bias:** exons with one live seam never enter the fit and are the
  worse-behaved half, and terminus-flanked exons are disproportionately one-seam. **Decide the
  minimum per-class pair count before falling back to the pooled scalar, and state it.**
  ⚠ Real data says this matters more than the toy: `ω_graft` spans **15×** across four cfRNA samples.

---

## 5. AFTER W2 — W3 THROUGH W7 (unchanged from the plan except where noted)

* **W3 — Python reference accumulator + tests, written to target and FAILING.** A1–A14, E1–E4, S1–S2.
  Lift `scratchpad/acc_proto_e.py` whole as S1 and the A3/A4 fixture; `acc_proto_g.py` as S2;
  `acc_proto_d.py`'s enumerators as the E1–E4 oracle. Run it on a **real cfRNA payload** through a shim
  before any C++ exists.
  ⚠ **Freeze three contract items in the header BEFORE W4**: the exact rounding of `R = round(2^32/L)`
  (`llround` vs banker's differ at ties; byte-identity is undefined without it); whether `L` is pre- or
  post-clip (observable in A14); and the per-object-class channel API.
* **W4 — C++ accumulator. ⛔ SHADOW MODE IS CANCELLED** (owner directive, 2026-07-29). Replace the C++
  accumulator outright, gated on the Python reference and the oracle bench. `set_graph(...)` as a
  **second** method, not an overload (the one-shot guard at `bam_scanner.cpp:1118-1122` throws today).
  ⚠ **Invert the W0.7 float-determinism test here** — that inversion IS the §6/A9 deliverable.
* **W5 — consumers, one observable per arm.** `W5a` contained counts · `W5b` contiguous-edge flux ·
  `W5c` junction edges + `E_J` **(and the reach A/B — reach earns its keep here or gets deleted)** ·
  `W5d` the §4 effective lengths — ⚠ **structural-zero masks BEFORE removing `_EPS`** (plan F6: node
  divisors collapse to exactly 0.0 on 12.4 % of v8 nodes; §8's "delete the floors" is right for edges
  and wrong for nodes) · `W5e` the `(count, Σ1/L)` moment into `node_init`, **gated on real cfRNA**.
  ⚠ Highest-risk mechanical change in the rework: the `_relay`/`_transport` and
  `_peel_share`/`_peel_share_scalar` twin pairs are deliberate, forbidden from merging (a measured
  15.7×/op), and the face index threads through six functions. **Every face-removal edit must be
  hand-mirrored into both arms with nothing enforcing it.** Re-measure on `cfrna:LBX0190`, never the toy.
* **W6 — optimization (budgeted).** Object count becomes 1.04 M nodes + 1.04 M contiguous + 404 K
  junction ≈ **2.49 M vs today's ~2.09 M chain nodes**.
* **W7 — deletions, then goldens.** `ruff check src/ tests/ scripts/` and treat undefined-name failures
  as the authoritative delete list. Plan §4 has the census with a **⛔ DO NOT DELETE** list of seven
  items that were each proposed and are each wrong.

---

## 6. STANDING METHODOLOGY — NON-NEGOTIABLE

* **No magic numbers.** Pause and discuss before ANY new constant or heuristic.
* **One thing varied per arm**, with a **pre-registered falsification test written FIRST and verified
  failing** against the current tree.
* **Re-record the baseline from the current tree in the same session**, both refits. If HEAD-vs-baseline
  is not 32/32, the BASELINE is broken, not your change.
* **Append to `accumulator_ledger.md` at EVERY arm**, as it lands, never retroactively.
* **ruff + `ruff format` + full suite on every arm.** Goldens at W7 only.
* ⚠ **Profile and gate on REAL cfRNA, never the 10 Mb toy.** It overstates contained efficiency 4.6×,
  understates multi-crossing 3.8×, has no alternative TSS/TES, is Poisson by construction, and — new
  this session — **over-represents the both-bits seam case 12×**.
* ⛔ **The owner drives commits. Do not commit unless asked.**
* ⛔ **No legacy, no backwards compatibility, no speculative code.** Anything kept only "for comparison
  with the old version" is a defect. Converge and delete.

---

## 7. ENVIRONMENT AND TOOLS

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
# ALWAYS OMP_NUM_THREADS=1.  After any src/rigel/native/ edit:
pip install --no-build-isolation -e .
```

| tool | use |
|---|---|
| `scripts/debug/pass0_oracle_bench.py --arm NAME` | THE A/B. `P0_REFIT=0\|3`, `P0_BENCH_OUT=...`, `--suite NAME`, `--report --vs A --new B` |
| `scripts/debug/graph_human_gate.py [INDEX]` | human-scale graph gate — I1–I13 incl. I3b + the I11 walk, census, §8 budgets |
| `scripts/debug/calib_cache.py` | cfRNA fast loop. `load(path, index=...)` verifies `partition_hash`; pass `index=` to `calibrate_from_cache` for the flags |
| `scripts/debug/z2.py` | THE `z2` denominator — `lin_var`, `z2`. Import it; never re-implement |
| `scratchpad/perf_scale.py`, `perf_identity.py` | the cfRNA perf loop |
| suites | `~/Downloads/rigel_runs/ambig_dense_10mb` (32 cond, oracle — **THE bench**) · `.../quick_3to1_5mb` (16 cond, oracle, +25 % partition — **the only suite that shows the partition effect**) · `.../cfrna/_calib_cache` (3 real, v8) · `.../refs/rigel_index_v7` (human) |

⚠ `/tmp/rigel_selfsolve` and `_selfsolve_cache` are namespaced by `partition_hash`; a stale payload
cannot silently load.

---

## 8. OPEN QUESTIONS AND CONCERNS FOR THE OWNER

1. ⚠ **cfRNA `f_g` is not what the 2026-07-14 memory note says.** Measured on v8: LBX0190 **0.0562**,
   MO_3021 **0.1264**, LBX0588 **0.7718**. The note records LBX0190 as "~15 % gDNA". There is **no v7
   comparator** (that index no longer exists) and the note predates D6, the hyperprior work and much
   else, so this is recorded as the v8 baseline rather than as a regression. **It is worth an
   explicit look before W2's real-data gating.** LBX0588 at 77 % gDNA in particular deserves a sanity
   check — is that sample genuinely that contaminated?
2. ⚠ **`rna_sense_frac` (κ) is 0.002–0.012 on all three cfRNA samples** — i.e. almost fully antisense.
   The strand model is untouched by this session's work, so this is pre-existing, but it is extreme
   enough that someone should confirm it is the expected read-orientation convention and not a bug.
3. ⚠ **`calibrate` wall-clock on v8 human is 137–198 s** (refit=3, and measured with ten review agents
   on the machine). `PERFORMANCE.md` records **67 s**, but at **refit=1** — the two are **not
   comparable**. Someone should re-measure like-for-like on a quiet machine before W6 is scoped. If it
   really is ~2× for +38.7 % objects, that is superlinear and worth understanding.
4. **`CalibrationResult`'s dtype gate accepts integers for a producer that does not exist until W5a.**
   The review confirmed this as speculative code. It was kept because the owner named it explicitly as
   a W1b item. **Revert it, or leave it?**
5. **I3b/I4/I11/I13 do not run at load** (§3 item 9), so a hand-edited or drifted `signature`/`flags`
   column loads clean. Reconstructing transcripts at load costs ~3 s at human scale. **Accept, or pay
   the 3 s?** A cheaper middle option exists: store a content hash of `nodes`+`edges` in
   `manifest.json` at build and check it at load — but that is a stored derived value, which is the
   thing `partition_hash` was deliberately designed to avoid.
6. **The `both terminus and splice site` case is 1.0 % on human and 11.8 % on the toy.** If W2 grows a
   branch for it, the toy will validate it against a 12×-inflated population. Worth deciding up front
   whether that case gets its own handling at all.
7. **Two W2 questions P1G_SCOPE leaves genuinely open** (it says so): *what is the right operator at a
   terminus* (three candidates, one ruled out), and *per-strand or per-boundary* (the measurement used
   a per-boundary binary; the per-strand refinement is untested and now cheap, since the flags are
   per-strand).
8. **Nothing is committed.** 562 modified/untracked paths, including the deletion of
   `calibration/regions.py` and 305 scratchpad files. A commit before W2 would make the next arm's
   `git diff` meaningful — currently it is not.

---

## 9. WHAT I WOULD DO FIRST IN THE NEW SESSION

1. `conda activate rigel`; `pytest tests/ -q` → expect **1280 pass**; `ruff check src/ tests/ scripts/`.
2. Re-record the baseline, both refits, arm `w2a_pre`; confirm **32/32 exact** vs arm `w1c_rev`.
   If it is not 32/32, **stop** — the baseline or the tree is broken, and nothing after that is valid.
3. Read §4 of this file, then `P1G_SCOPE.md` §7 C3 and §8.
4. Write W2a's falsification test (§4.5 P4) and **verify it fails**.
5. Implement §4.4, being careful about step 5 (the per-face split of `node_total_density`).
6. Measure, record in the ledger, and only then move to W2b.
