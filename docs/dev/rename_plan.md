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

### 2a. `cut` — ⭐⭐ RESPECIFIED BY THE OWNER 2026-08-13, and the answer inverted

**The ruling:** a BOUNDARY **includes the terminal anchors** — an intergenic↔exon terminus is a boundary.
Chromosome-end boundaries at `0` and `L` are *not required* but **may be added if they help
algorithmically**, and then a chromosome with `N` regions has **`N+1` boundaries**.

⛔ **MY EARLIER RECOMMENDATION (`cut` → `region_bound`) IS RETRACTED.** It was built on reading BOUNDARY
as interior-only, which the respecification overturns. Measured on the shipped index, the arithmetic is
exact and it lands the other way:

| array | count | what it is |
|---|---:|---|
| `cut_positions` | **35,229** | `N + one per chromosome` — ⭐ **already exactly the anchored `N+1` boundary set** |
| regions `N` | 35,135 | |
| today's `edge` axis | **35,041** | `N − one per chromosome` — the **interior-only** subset, anchors excluded |

⭐ So `cut` → **`boundary`** after all, as `region_boundary_sj_design.md` §5 originally said. What §5
missed is that the codebase carries the boundary set **twice, at two different extents**, differing by
exactly 2 per reference (188 slots).

> **RECOMMENDATION — UNIFY THE TWO AXES *BEFORE* RENAMING, NOT DURING.**
> Grow the deposit axis from `N−R` to `N+R` so there is ONE boundary concept at ONE length. The two
> terminal slots per reference simply carry zero crossings — nothing crosses a chromosome end — so it is
> additive, and 188 zero slots out of 35,229 is free.
>
> ⭐⭐ **It is not cosmetic: it makes §5's stated arithmetic actually TRUE.** §5 already claims
> `region r → boundaries (r, r+1)` and `boundary b → regions (b, b+1)` are ARITHMETIC. Today that is
> false at every reference end, because the edge axis skips the anchors; with the anchored axis it holds
> everywhere with no special case. It also **deletes a population** rather than renaming one.
>
> ⛔ **But it is a STRUCTURAL change, not a rename** — the payload's boundary axis changes length — so it
> gets its own commit and its own gate (`one-thing-varied`). Doing it FIRST is what converts the hardest
> naming ambiguity into a non-issue, instead of carrying `two-masks-one-name` through 690 rename sites.

⚠ **If the owner prefers not to unify**, the fallback is to keep both and name the extents apart —
`boundary_positions` (`N+R`, anchored) versus the interior deposit axis. That works, but it preserves the
ambiguity the rename exists to remove, so it is the worse option and is recorded only as the alternative.

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

⭐⭐ **The quant digest must not change at ANY stage — stage 0 included.** Stage 0 alters the payload's
SHAPE, so ① legitimately moves there; the OUTPUT must not, because adding all-zero slots changes no
number. That makes stage 0's gate sharper than the renames', not looser.

### 3d. FROZEN, not rolling — this is the answer to "the renames compound"

Every stage compares to the SAME capture, never to the stage before it. ⛔ A rolling baseline lets a
defect introduced at stage 2 become the accepted truth for stages 3-9 — which is precisely the
compounding failure. This is the one place `TRAPS: re-record-the-baseline` is INVERTED, and deliberately:
the claim under test is *"nothing has moved since the freeze"*. `--freeze` refuses to overwrite an
existing reference for that reason; it is re-frozen exactly once, after stage 0.

### 3e. Falsified, on the comparator AND on real data

`--self-test` is 8/8: a pure rename is silent; a changed value, a shape change, a moved output, a
DISAPPEARING array and an EMPTY capture are all caught; and ①'s name-swap hole is pinned as a hole so
nobody mistakes it for coverage. ⭐ On real data a **one-ULP** nudge to a single element of a
13,482-element array was injected and the gate fired, then read clean again on restore.

---

## 3f. THE STAGING — each stage is one commit, and lands only on a green `--check`

⭐ Ordered by **risk of collision**, ascending, so the method is proven on the safe stages first.

| # | stage | scope | gate beyond suite + `--check` |
|---|---|---|---|
| **0** | ⭐⭐ **STRUCTURAL: unify the boundary axis at `N+R`** (§2a) | payload axis 35,041 → 35,229; 188 new slots | ⛔ the new slots must be **all zero**; the quant digest must NOT move; then **re-freeze** |
| ✅ 1 | `donor`/`acceptor` → `boundary_left`/`right` | 117 sites, 17 files | done in `4576c3be` (suite-gated; predates the harness) |
| 2 | **`node` → `region`** | ~2,700 code + 1,264 docstrings; exempt the 10 `ast` sites | C++/spec parity |
| 3 | **`edge`/`seam` → `boundary`** | ~1,800 code + 877 docstrings | C++/spec parity |
| 4 | **`line` → `boundary`**, per site | 120 of 180; leave the 60 text lines | — |
| 5 | **`cut` → `boundary`** | ~690 code + 214 docstrings | — |
| 6 | **`junction` → `sj`**, incl. §2.5's `mass`→`count` | ~900 code + 646 docstrings | — |
| 7 | **module filenames** | 16 files + every import | — |
| 8 | **docs, `CLAUDE.md`, `DESIGN.md` §0, memory** | 1,080 | `test_docs_boundary`, `test_no_jargon_labels`; ⭐ `--check` is a NO-OP here and that is the point: prose cannot move a number |
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
