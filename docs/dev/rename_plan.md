# The bulk rename — audit and plan

    ⚠ **A DEV DOC.** Provisional; nothing may cite it. When a stage lands, MOVE what settled — the
    vocabulary to `DESIGN.md` §0, a lesson to `TRAPS.md` — and delete it here in the same edit.

    ⭐ Written 2026-08-13 against `4576c3be`, the checkpoint commit taken *before* any of this.
    Re-derive the numbers with `python scripts/design/rename_census.py`; do not trust the tables below.

---

## 0. THE VOCABULARY (owner, 2026-08-12, sharpened 2026-08-13)

| term | what it IS | replaces |
|---|---|---|
| **REGION** | a genomic **interval** | `node` |
| **BOUNDARY** | a single genomic **position between two regions** | `edge`, `cut`, `line`, `seam` |
| **SPLICE JUNCTION** / **sj** | a connection **between two boundaries**, non-contiguous. ⭐ Both spellings allowed; ⛔ bare `junction` is not | `junction` |
| — | ⛔ **splice donor / splice acceptor are BANNED** as names for a genomically-ordered pair | ✅ done, stage 1 |

⭐ **"Between two regions" is load-bearing and it decides the plan** — see §2a.

---

## 1. THE MEASURED LANDSCAPE — ~10,900 occurrences

| where | distinct | occurrences |
|---|---:|---:|
| Python identifiers (class / function / argument / attribute / variable) | 845 | **5,192** |
| C++ identifiers | 86 | **978** |
| Python **docstrings** | — | **3,653** |
| `docs/*.md` + `CLAUDE.md` | — | **979** |
| **memory** files (`~/.claude/.../memory/`) | — | **101** |
| module **filenames** (an import-graph change) | 16 | — |
| **payload wire schema** names ⛔ refuses every cache | 12 | — |

⛔⛔ **PROSE IS 43 % OF THE JOB — 4,733 of ~10,900.** A code-only rename leaves nearly half the
vocabulary stale, which is the exact rot this exercise exists to remove. Per token, prose dominates
everywhere: `node` 1,264 docstrings against ~2,700 code sites; `junction` 646; `edge` 783.

⚠ **The memory directory is in scope** and has 101 (`node` 62, `junction` 19, `edge` 10).

---

## 2. ⛔⛔⛔ FOUR AMBIGUITIES — a global replace is WRONG for each, and each needs a ruling

### 2a. `cut` — ⭐⭐⭐ ANSWERED FROM THE CODE: **terminal boundaries serve no purpose. There is no stage 0.**

The question asked was the right one: *do terminal boundaries make the code simpler or more elegant?*
Measured, the answer is **no**, and it was already litigated — `node_chain.py` says so in its own words:

> *"It starts and ends with a NODE, and there are no terminal slots. That is the whole shape change from
> the predecessor, which ran `B R B R … R B` with **k+1 boundary slots per reference — the two outermost
> carrying no data and existing only so every region had an object on each side**. A contiguous edge is
> the line BETWEEN two adjacent nodes; there is no such line before the first or after the last, so **the
> object does not exist rather than existing empty**. An edge therefore always has a node on both sides,
> **an invariant the old shape could not state.**"*

⛔⛔ **The `N+1` anchored layout IS the predecessor, and it was deliberately deleted.** Proposing it here
was on the way to rebuilding a refused structure — the `ROADMAP.md` §4 failure mode, in structural form.

**The ledger, measured rather than argued:**

| | |
|---|---|
| special cases they would DELETE | **3, all trivial** — one `flatnonzero(right >= 0)` in `edge_node_indices`, two `edges_[line - 1]` in the C++ deposit |
| special cases they would ADD | a `-1` in `boundary → regions` (a terminal boundary has no region on one side) |
| slots they would ADD to the SOLVE | **188**, all permanently empty — `NodeChain` builds the solve from this axis, so the quant digest would move |
| invariants they would BREAK | *"an edge always has a node on both sides"*, and the chain's *"starts and ends with a NODE"* — which `sweep.py`'s propagation-sink logic rests on |
| modelling need | **none.** Nothing can cross a chromosome end, so a terminal boundary has no crossing population, permanently |

⭐ **And the vocabulary resolves with no structural change at all, because `cut` and `boundary` were
never the same thing:**

| name | count | what it is |
|---|---:|---|
| **REGION** | 35,135 | a genomic interval, with `start`/`end` |
| **BOUNDARY** | **35,041** | the position BETWEEN two adjacent regions. ⭐ It does not exist before the first region or after the last — *the object does not exist rather than existing empty* |
| `cut_positions` → **`region_bounds`** | 35,229 | ⛔ NOT a boundary axis. The packed region ENDPOINTS — a coordinate grid, an input to the partition, not a derived object |

⚠ **A flip-flop, stated plainly so the record is honest.** `region_bound` was the original
recommendation; it was retracted after reading the respecification as *"boundaries include the
chromosome ends"*. Re-reading, the respecification says **terminal anchors (intergenic↔exon)** — which
are already ordinary INTERIOR boundaries and always were — and says explicitly that boundaries at `0`
and `L` are *not needed*. The retraction was mine and was wrong; the code settles it in the same
direction as the original reading.

### 2b. `junction` — bare form banned, two allowed spellings

> **RECOMMENDATION — `sj` in identifiers, "splice junction" in prose.** `sj` is already the accepted
> short form and keeps names readable: `n_junctions → n_sj`, `JunctionEdge → SpliceJunction`,
> `junction_opportunity.py → sj_opportunity.py`, `junctions_df → sj_df`.

⭐ **One place where bundling SAVES a change**: `ROADMAP.md` §2.5 already wants `mass_rna_junction →
count_rna_junction`, because the field is an incidence and not a mass. The wire schema is being
rewritten here anyway, so land both at once as **`count_rna_sj`** rather than touching it twice.
⚠ State it in the commit — it is two renames in one field, both provably numeric no-ops.

### 2c. `line` — 120 boundary sites vs **60 ordinary lines of text**

`for line in text.splitlines()` is 60 sites. A global replace corrupts file-reading code in the
instruments. ⇒ **per-site**, no automation.

### 2d. `node` — and `ast` nodes, including in the census tool itself

10 sites are `ast.NodeVisitor` / `generic_visit` / `node` locals in `module_census.py` and
`rename_census.py`. ⇒ exempt by file, not by token.

### 2e. ✅ `donor` / `acceptor` — DONE (stage 1, in `4576c3be`)

117 identifiers renamed to `boundary_left` / `boundary_right`. **Two exemptions are load-bearing** and
must not be undone:
* `donor` already had a second sense — the toy harness's **source condition** (`donor_dir`, `donor_on`,
  `_donor_sim_params`). Renaming those would have corrupted the simulator.
* `sim/splice_motif.py::splice_donor_acceptor` takes the **strand** and returns the GT..AG
  dinucleotides. It is the one place the biology term is used correctly, and it stays.

---

## 3. ⭐⭐⭐ THE GATE — every stage PROVEN BIT-IDENTICAL, not merely suite-green

**`scripts/design/rename_identity.py`.** A rename is a numeric no-op; this makes that claim falsifiable.

### 3a. The keystone: bit-identity is achievable, but only pinned

⛔ The shipped scan is **not** bit-reproducible — two scans of the same BAM by the same binary differ on
six banks by ~3.5e-14, which is what left `rescan_panels.py` unsatisfiable for three days.
⭐⭐ **Measured 2026-08-13: `total_threads=1` removes it entirely — 0 of 62 banks differ.** So the harness
pins `total_threads`, `bgzf_threads`, `OMP_NUM_THREADS` and the EM seed, and only then is "bit-identical"
a claim anyone can make. ⚠ A `--check` therefore says nothing about the shipped thread configuration; it
is a statement about the CODE, which is exactly what a rename needs.

### 3b. It compares CONTENT, never names — because the names are what is changing

A gate keyed on bank names would fail on every stage by construction and prove nothing. Two
rename-invariant views, and they are complementary:

| | view | catches | blind to |
|---|---|---|---|
| ① | **content multiset** — sorted `(dtype, shape, sha256)` of every array, names discarded | any VALUE moving, any array appearing/disappearing/reshaping | two arrays having their names SWAPPED |
| ② | **the tool's output** — the transcript table, which carries no region/boundary vocabulary at all | a swap, and anything else that moves an answer | nothing relevant |

⛔ **No name map, deliberately.** It would close ①'s hole, and it is the compatibility hack that was
refused — it would have to be maintained across nine stages while being exactly the thing under test.

### 3c. ⭐ ONE NUMBER TO WATCH, FROM NOW UNTIL STAGE 9

    quant digest   ac0b82b4e71ae431…      8,750 rows
    payload        18 arrays              calibration  19 arrays
    frozen at `4576c3be` + the audit, on `g00 ss0.99 capture_off`

⭐⭐ **The quant digest must not change at ANY stage, and now there is no exception.** With stage 0
struck (§2a), every remaining stage is a PURE RENAME — so ① and ② must BOTH be unchanged, every time.
The reference is frozen once and never re-frozen.

### 3d. FROZEN, not rolling — this is the answer to "the renames compound"

Every stage compares to the SAME capture, never to the stage before it. ⛔ A rolling baseline lets a
defect introduced at stage 2 become the accepted truth for stages 3-9 — which is precisely the
compounding failure. This is the one place `TRAPS: re-record-the-baseline` is INVERTED, and deliberately:
the claim under test is *"nothing has moved since the freeze"*. `--freeze` refuses to overwrite an
existing reference for that reason, and with stage 0 struck it is **never re-frozen at all**.

### 3e. Falsified, on the comparator AND on real data

`--self-test` is 8/8: a pure rename is silent; a changed value, a shape change, a moved output, a
DISAPPEARING array and an EMPTY capture are all caught; and ①'s name-swap hole is pinned as a hole so
nobody mistakes it for coverage. ⭐ On real data a **one-ULP** nudge to a single element of a
13,482-element array was injected and the gate fired, then read clean again on restore.

---

## 3f. THE STAGING — each stage is one commit, and lands only on a green `--check`

⭐ Ordered by **risk of collision**, ascending, so the method is proven on the safe stages first.
⛔ **There is no stage 0** — §2a struck it. Every stage below is a PURE RENAME, so every one of
them must read ✅ on `--check` with no exception and no re-freeze.

| # | stage | scope | gate beyond suite + `--check` |
|---|---|---|---|
| ✅ 1 | `donor`/`acceptor` → `boundary_left`/`right` | 117 sites, 17 files | done in `4576c3be` (suite-gated; predates the harness) |
| 2 | **`node` → `region`** | ~2,700 code + 1,264 docstrings; exempt the 10 `ast` sites | C++/spec parity |
| 3 | **`edge`/`seam` → `boundary`** | ~1,800 code + 877 docstrings | C++/spec parity |
| 4 | **`line` → `boundary`**, per site | 120 of 180; leave the 60 text lines | — |
| 5 | **`cut` → `region_bound`** (§2a — it is the region ENDPOINT grid, not a boundary) | ~690 code + 214 docstrings | — |
| 6 | **`junction` → `sj`**, incl. §2.5's `mass`→`count` | ~900 code + 646 docstrings | — |
| 7 | **module filenames** | remaining `edge_*`/`junction_*` files + imports | ⚠ ONLY the ones whose stage did not already move them — see below |
| **8a** | ⭐ **ON-DISK ARTIFACT names** — `nodes.feather`, the `node_id` column, `gdna_density_nodes.feather` | 24 strings | ⛔ needs an INDEX REBUILD, so it lands with stage 9 |
| 8b | **docs, `CLAUDE.md`, `DESIGN.md` §0, memory** | 1,080 | `test_docs_boundary`, `test_no_jargon_labels`; ⭐ `--check` is a NO-OP here and that is the point: prose cannot move a number |
| 9 | **the caches** — one regeneration | 4 panels | `rescan_panels.py` |

### The per-stage protocol, in order — no step is optional

1. `python scripts/design/rename_census.py` — the leftovers list is this stage's checklist.
2. Apply the rename **longest-identifier-first, word-boundary-anchored**. `sj_acceptor_cut_` must not be
   half-eaten by `acceptor_cut`; stage 1 used exactly this and it held.
3. Rebuild if `src/rigel/native/` moved: `pip install --no-build-isolation -e ".[dev]"`.
4. `ruff check src/ tests/ scripts/` — ⚠ never `ruff format scripts/`.
5. `python -m pytest tests/ -q` — green, with the count delta **re-derived** via
   `--collect-only | grep <new file stem>`, never adjusted.
6. ⭐⭐ `python scripts/design/rename_identity.py --check --stage <n>` — **must be ✅**.
7. `rename_census.py` again: the stage's token now appears only in exempt or prose-correct positions.
8. Commit. One stage, one commit, so any stage can be reverted alone.

⛔ **If step 6 fails, STOP and revert the stage.** A rename that moves a number is not a rename, and the
next stage would build on it.

## 3g. ⛔⛔ THREE THINGS STAGE 2 TAUGHT, ALL OF WHICH CHANGE THE REMAINING STAGES

**① A MODULE NAME IS AN IDENTIFIER, so stage 7 cannot be deferred per-module.** Renaming
`node_geometry` → `region_geometry` rewrites `import` statements, and the suite then fails to COLLECT
until the file moves. ⇒ **a module's file moves in the same stage as its contents**; stage 7 keeps only
whatever a stage did not already move.

**② PERSISTED ARTIFACT NAMES ARE NOT IDENTIFIERS AND MUST NOT MOVE WITH THEM.** `nodes.feather`, the
`node_id` COLUMN inside it, `gdna_density_nodes.feather`. Renaming those invalidated the shipped index
mid-stage and made the bit-identity gate unrunnable — so they were reverted and given stage 8a, which
needs an index rebuild and therefore belongs with the caches. ⭐ The right shape is
`REGIONS_FEATHER = "nodes.feather"`: the CONSTANT renamed, the on-disk VALUE not.

**③ ⭐⭐ THE NEW VOCABULARY IS THE PREDECESSOR'S VOCABULARY, and two anti-regression gates were keyed on
that.** `test_the_old_BOUNDARY_vocabulary_is_GONE` banned the attributes `REGION` and `BOUNDARY`;
`test_the_per_face_fields_are_gone` banned the fields `gdna_region_eff_len` and `n_regions`. All four are
names the ruling RE-ADOPTS, so the gates began banning live code.
⛔ **The words were never the defect — the `B R B R` SHAPE was.** Both gates are re-pointed at the shape:
`k` regions give `2k − 1` slots, the chain starts and ends with a REGION, and every BOUNDARY has a region
on both sides — the invariant the predecessor could not state. That is strictly stronger than the name
check it replaces (`TRAPS: a-gate-that-reconstructs` — prefer an invariance).
⚠ **Expect the same collision in stages 3–6**, since `boundary` is the other re-adopted word.

**⚠ And one false positive the census did not predict:** `request.node` is **pytest's API**. It was
renamed to `request.region` and two scenario tests failed. Nothing else in the tree hit a third-party
`node` API, but the class of error is real — a token can belong to a LIBRARY, not just to two of our own
senses.

## 4. THE METHOD, and the two rules it follows

1. **Identifier-by-identifier, word-boundary-anchored, longest-first** — `sj_acceptor_cut_` must not be
   half-eaten by `acceptor_cut`. Stage 1 used exactly this and it worked.
2. **Re-run `rename_census.py` after every stage.** The leftovers list is the checklist; a stage is done
   when its token appears only in exempt or prose-correct positions.
3. ⚠ **Do not run `ruff format` over `scripts/`** — `CLAUDE.md` forbids it.
4. ⛔ Rebuild the C++ after any stage touching `src/rigel/native/`.

---

## 5. WHAT THIS PLAN DELIBERATELY DOES NOT DO

* ⛔ **The structural half of `region_boundary_sj_design.md` §5** — the boundary-event bitmask, the two
  sj CSRs (`sj_at_left`/`sj_at_right`), collapsing five sj structures into two tables. Those change
  STRUCTURES, not names, and `one-thing-varied` says they follow the rename rather than ride it.
* ⛔ **Re-keying the sj CSR from cut index to boundary index** — same reason. ⭐ Note that stage 0's
  unification makes this nearly free afterwards: once the boundary axis IS the cut axis, the CSR is
  already keyed on a boundary index and only the NAME changes.
* ⛔ **`SJStrandTable`'s 1.576× disagreement with the accumulator** — blocking for the SJ output, not
  for the rename.
