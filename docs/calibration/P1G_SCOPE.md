# P1g — TSS/TES in the region/boundary map: the problem, and the scope

**Written 2026-07-27 for a fresh session.** Read `gdna_reframe_terminus.md` first — it is the measurement
this brief acts on. Nothing here is implemented; `src/` is bit-identical 32/32 to the P-2 baseline.

---

## 1. THE PROBLEM IN ONE PARAGRAPH

The calibration region/boundary map records, per genomic interval, a **4-bit signature** — `{intron±,
exon±}` — and nothing else. A **boundary** row is `[boundary_id, ref_name, position]`. So the solver knows
that a seam separates an exon from an intron, and it does **not** know whether that seam is a **splice
junction** (RNA crosses it, spliced) or a **transcript terminus** (RNA stops there). Those two are
physically opposite and the solver treats them identically.

This was already on record as the structural debt behind `ω_graft`, which splits **≥30×** on exactly that
bit (1.7–1.9 at termini vs 0.04–0.06 at junction-only seams). **It is now measured to do a second, larger
and independent kind of damage: it corrupts the reframe, as a MODE error of up to 190×.**

## 2. THE EVIDENCE (2026-07-27, `gdna_reframe_terminus.md`)

The gDNA density a message delivers into an exon is **5.4× too large** at capture-OFF. Decomposing it —
an identity, residual 4.4e-16, over 1,289 exon destinations:

| term | mass-wt mean |
|---|---|
| **TOTAL** `log10(delivered/true)` | **+1.564 dec (37×)** |
| the source's own gDNA error | +0.110 |
| **the REFRAME `log10 r`** | **+1.508 — 96 %** |
| true spatial difference | −0.054 (gDNA is uniform ✓) |

split by the structural bit:

| source boundary | n | median `log10 r` | median delivered gDNA error |
|---|---|---|---|
| **TERMINUS** | 426 | **+0.836** | **7.0× too big** |
| **junction-only** | 732 | **+0.021** | **1.0× — exact** |

**33 % of edges carry 66–68 % of the error mass.** Reproduced on `gdna300` and on the stranded twin.

### 2.1 ⭐ The closed form — and it predicts the measurement to 3 %

`r = ρ_tot(dst)/ρ_tot(src)`. At a terminus **no RNA crosses**, so the boundary face is pure gDNA,
`ρ_tot(src) = ρ_g(src)`, and with `e` the true capture step (`ρ_g(dst) = e·ρ_g(src)`):

```
    r  =  ρ_tot(dst) / ρ_g(src)  =  e / f_g(dst)
    tg =  ρ_g(src) · r           =  ρ_tot(dst)  =  ρ_g(dst) / f_g(dst)
```

> **The delivered gDNA density is too large by EXACTLY `1/f_g(dst)` — the reciprocal of the destination's
> own composition, which is the very quantity being estimated.**

| node | true `f_g` | predicted `1/f_g` | **measured** |
|---|---|---|---|
| 2651 | 0.0054 | 185.2× | **190.1×** |
| 3087 | 0.0080 | 125.0× | **127.6×** |

Verified independently: on **63–66 %** of terminus edges `tg == ρ_tot(dst)` to within **1e-9**, and the
composition it implies is `f_g = 1.18–1.27` — the message tells an exon that is 99.5 % RNA that it is more
than 100 % gDNA.

### 2.2 Why this is the same shape as the pin bug

`_pin_v` manufactured a message from the destination's own **belief**; this manufactures one from the
destination's own **total density**. Both produce a claim carrying zero source information. That is why P-2
exposed it — the pin was scaling it back down by ×0.648 and hiding it.

**It also sets the constraint on any fix:** the error factor is `1/f_g(dst)`, so *any* correction that
divides it back out using the destination's belief re-creates the pin bug. **The fix must be structural —
do not form that ratio — not corrective.**

## 3. WHAT THE BIT ACTUALLY IS

"Terminus" is shorthand. The quantity the solver needs is:

> **does transcribed RNA continue across this seam, on this strand, on this side?**

Three refinements the design must decide on, in increasing order of how much they matter:

1. **Per strand.** A seam can be a terminus for `+` and a junction for `−`. The reframe faces are already
   per-side; the strand split is not represented at all.
2. **Per side.** `_rho_faces` builds a left face and a right face. A terminus is asymmetric by nature (RNA
   exists on one side and not the other), so a single per-boundary flag under-specifies it.
3. **Not binary.** 73 positions in this annotation are **both** a terminus (for one transcript) and a
   junction (for another). The honest quantity is a *continuing share*, which is exactly M10's peel share
   `w = ρ_ν/(ρ_ν + ρ_μ)` — already derived and implemented, but estimated from **data**, not structure.

⚠ Today the solver answers this question from **observation**: `accept_l/accept_r = (SP+SN) > 0`
(`bp_solver.py:381`) — "did we see spliced mass on this face". That conflates *structurally no RNA crosses*
with *no RNA observed here*, and it is silently wrong on any low-coverage junction.

## 4. ⭐ P1g IS A RESTORATION, NOT A NEW BUILD

**Index format v6 already carried these columns and v7 removed them.** `index.py` version history:

```
#   6 — three-channel calibration: … boundaries.feather carried per-boundary annotation flags
#        (is_tss/is_tes/is_splice_junction/genomic_sj_strand).
#   7 — those v6 precompute columns removed (the message-precision collapse retired the
#        mature/nascent overlay that consumed them …)
```

Removed by commit **`1863ef57`** (2026-07-09, "ship disagreement-variance message precision"), which also
dropped `build_boundary_partition`'s `transcripts` argument. **The builder is recoverable verbatim:**

```bash
git show 1863ef57^:src/rigel/calibration/regions.py    # ~25 lines, the sj/tss/tes event sweep
```

Its logic, for reference — per transcript, skipping synthetic rows:

```python
for intron_start, intron_end in tx.introns():      # every intron endpoint is a splice junction
    sj.add(intron_start); sj.add(intron_end)
first, last = tx.exons[0].start, tx.exons[-1].end
tss_p = first if tx.strand == POS else last        # 5' terminus, orientation by strand
tes_p = last  if tx.strand == POS else first
```

`BOUNDARY_COLUMNS` and `BOUNDARY_COLUMN_DTYPES` both existed, including
`genomic_sj_strand: int8  # 0 none / 1 POS / 2 NEG / 3 both`.

**This is the single biggest scoping fact: the annotation half of P1g is a revert with a version bump.**

## 5. ⛔ CAN WE SKIP THE INDEX AND DERIVE IT FROM THE SIGNATURE? — MEASURED: NO

Tempting, because the signature ships today and the calibration path already has it. A boundary would be a
terminus for strand `s` if `s` has coverage (exon **or** intron) on one side and not the other. Measured
against the annotation-derived truth on all 1,699 boundaries:

| | |
|---|---|
| agreement | **83.3 %** (1415/1699) |
| precision / recall | 0.695 / 0.868 |
| **false negatives** | **73** |
| **false positives** | **211** |

Both error classes are explained **exactly**, and the explanations are the verdict:

* **FP = 211 = 100 % nRNA-span edges.** The signature includes synthetic nRNA spans; the annotation
  classifier skips them. A *definitional* disagreement, not an error — fixable either way.
* ⛔ **FN = 73 = precisely the set of positions that are BOTH a terminus and a junction** (verified as a
  subset, exactly). Where one transcript ends at a point another spans, transcript coverage is continuous,
  so the **union** signature cannot see the terminus at all. **The signature is structurally blind on
  exactly the mixed seams that §3.3 says matter most.**

**→ Restore the exact annotation bit. Do not derive it from the signature.** It is cheaper (the code exists)
and it is exact.

## 6. THE PLUMBING — what changes, and what does NOT

### Does not change

* ⭐ **No C++.** The accumulator keys on `boundary_id` and consumes
  `build_region_partition_arrays(index) -> (boundary_positions, ref_pos_offsets, region_types)`. Adding
  columns to `boundary_df` does not touch that ABI. **No `_bam_impl` rebuild, no accumulator spec change,
  no `tests/native/` change.**
* No change to `region_df`, the signature encoding, or `RegionArrays`.
* No change to the payload schema (`AccumulatorPayload`).

### Does change

| file | change | size |
|---|---|---|
| `calibration/regions.py` | restore the 4 columns + the event sweep; `build_boundary_partition` regains its `transcripts` arg | ~30 lines (revert) |
| `index.py` | `INDEX_FORMAT_VERSION` 7 → **8**; version-history note; pass `transcripts` at the build call site | ~10 lines |
| `index.py` loaders / `validate_boundaries_against_regions` | accept the new columns | small |
| **⚠ NEW: a path from `boundary_df` to `calibrate`** | **the calibration package never reads `boundary_df` today** — verified by grep. `calibrate()` takes `(payload, region_arrays, …)`. The bit must reach `node_geometry` / `bp_solver`. | **the real work** |
| `node_geometry.py` | carry the per-boundary, per-strand bit into `NodeGeometry` (it already carries ~20 per-face arrays, so this is idiomatic) | ~20 lines |
| `bp_solver.py` | the consumer (§7) | design-dependent |

**The plumbing question is the one genuine unknown in the scope.** Two shapes, pick one early:

* **(A) widen `RegionArrays`** into a `PartitionArrays` that carries both region and boundary annotation.
  Touches every `RegionArrays.from_region_df` call site (many, including every debug script).
* **(B) a separate small `BoundaryAnnotation` array passed alongside**, threaded
  `pipeline → calibrate → build_node_geometry`. Fewer call sites; one more argument.
  **Recommended** — it matches how `region_arrays` is already threaded and keeps the blast radius small.

⚠ **Every existing index must be rebuilt.** `~/Downloads/rigel_runs/*/rigel_index`, the cfRNA index, and
the test fixtures. Budget for that: it is the slow part, not the code.

## 7. THE CONSUMERS — what to DO with the bit, ranked

This is the design work; the sections above are mechanical. **Do not implement a consumer before running
the A/B protocol in §8 on it.**

### C1 ⭐ — the reframe on the gDNA component (the new finding, and the big one)

The reframe's premise is *"my neighbour and I share a composition"*. **At a terminus that premise is void by
construction** — one side has RNA and the other cannot. So:

* ⛔ **Do NOT "correct" `r` by `f_g(dst)`.** The closed form (§2.1) says that is exactly the missing factor,
  and it is the destination's own belief — the pin bug rebuilt. This is the trap.
* ⭐ **The precedent already in the code is the λ-EMISSION GATE.** `bp_solver.py` gates the composition
  message on the source supplying both components: *"A source that carries only ONE component has no such
  claim to make — λ is not 'large' for it, it is UNDEFINED"*, and its own comment continues:
  **"a pure-gDNA source's authority is a DENSITY LEVEL, and that travels on the measurement stream
  (`tmg`)"**. That is precisely right — **and the measurement stream is still multiplied by `r`.** The
  three-stream separation is correct in design and the reframe is applied uniformly across all three, which
  is the defect. **The narrowest fix consistent with the existing architecture: a pure-gDNA source's
  density level must not be reframed by a total-density ratio.**
* **A gDNA-specific frame.** `r_g = ρ_g^init(dst)/ρ_g^init(src)` instead of the total ratio. At
  `_RHO_ITERS = 1` the frames come from `_init_belief()`, which is belief-free (P4b), so this is at no worse
  BP standing than today's `r`. ⚠ Untested, and on unstranded data the init composition is the
  uninformative ½, so it may just rescale the same error. **Measure before building.**

**Sizing (already done):** `r_g ≡ 1` everywhere takes unstranded × capOFF × gDNA **0.0495 → 0.0146** at r0
and **0.0313 → 0.0077** at r1 (6/6 better), capture-OFF overall **−0.0178 / −0.0124** (9 b / 0 w) — while
costing unstranded × capON **+0.1728**. So the prize is real and the constraint is sharp: **any operator
must reduce to `r` where the seam is a junction and a capture step is genuinely being crossed.**

### C2 — `ω_graft` per structural class (already specified)

The ROADMAP's standing instruction: *"same estimating equation, one scalar per structural class, each
fitted over a homogeneous population; the partial-pooling block already in `graft_premise_logvar` is the
plug-in point."* `ω̂` is measured at 1.7–1.9 at termini vs 0.04–0.06 at junction-only seams.
⚠ **Real data says this matters more than the toy did:** `ω_graft` spans **15×** across four real cfRNA
samples (0.25 → 3.83), with the two strands agreeing closely *within* each sample.

### C3 — replace the observational `accept_l/accept_r` with the structural bit

`(SP+SN) > 0` currently means "spliced mass was observed on this face". With the bit available this becomes
"RNA structurally crosses here", which is the intended semantics and is coverage-independent.
⚠ Behaviour-changing on every library; needs its own A/B. Lowest priority of the three, but it is the one
that makes the other two coherent.

### C4 — re-examine P1e's scope

`variance_ledger.md` §6 says P1e must SHRINK once the bias strata are diagnosed. Terminus edges are a named
bias stratum now. Re-measure `δ`'s bias share by structural class after C1 lands.

## 8. THE A/B PROTOCOL AND PRE-REGISTERED GATES

Non-negotiable, from the standing methodology:

1. **Re-record the baseline from the current tree in the same session**, both refits. If HEAD-vs-baseline
   is not 32/32, the *baseline* is broken.

   ⚠ **The reference quoted here was STALE and is corrected (2026-07-29).** It read
   *"r0 0.078786 / **r1** 0.052470"*, but HEAD has run `calib_refit_iters = 3` since `daa32a13` — the
   refit=1 arm is retired, and the r0 figure predates the P1e/landscape work. **Current reference,
   re-recorded from the working tree at `3c293038` this session and reproduced 32/32 exactly at both
   settings: `refit=0` **0.079005** / `refit=3` **0.046675**.** See
   `../accumulator/accumulator_ledger.md` W0.1. Re-record again before any arm — every stored number goes
   stale the moment the tree moves.
2. **Vary one thing per arm.** The index change and each consumer are separate arms. The index change alone
   must be **bit-identical 32/32** — it adds data nobody reads yet. That is the wiring gate.
3. **Pre-register predictions including a falsification test** before measuring.

**Suggested pre-registration for C1:**

* terminus edges improve, junction-only edges are **unmoved** — the falsification test, because junction
  edges are already exact (1.0×) and any movement there means the bit is being applied too widely;
* unstranded × capOFF × gDNA improves materially (the sizing arm says up to −0.0349 at r0);
* **capture-ON must not regress** — this is the gate the `r_g ≡ 1` arm fails at +0.1728;
* `gdna_none` unmoved or better (the sizing arm gave −0.0031);
* stranded capture-OFF improves slightly (its terminus split is the same: 0.796 vs 0.017).

**Also gate on:** ruff, the full test suite, held-fixed `z2` (do not regress), and **goldens LAST** — an
index version bump plus a solver change will move them, so regenerate once at the end.

## 9. RISKS AND OPEN QUESTIONS

* ⚠ **The synthetic suite's annotation is simple.** Nested truncations plus exon skips; no alternative
  TSS/TES *inside* exons, no mutually exclusive exons, no retained introns. Real annotations have all of
  these, and an alternative TSS inside an exon does not even produce a region boundary — **so a share of
  real termini will be invisible to a boundary-keyed bit no matter how it is built.** Size this on the real
  cfRNA index before trusting the toy's 33 %/68 % concentration.
* ⚠ **Index rebuild cost and coordination.** Every cached index and every scenario suite. The
  `_selfsolve_cache` payloads are keyed to an index; check whether they must be regenerated too.
* ⚠ **`/tmp/rigel_selfsolve` is shared and non-namespaced** — two sessions doing scenario work corrupt each
  other. This has cost two full rebuilds already.
* **Open: what is the right operator at a terminus?** §7/C1 gives three candidates and rules one out. The
  honest state is that the *diagnosis* is closed and the *operator* is not.
* **Open: per-strand or per-boundary?** The measurement used a per-boundary binary and got a clean 2×
  concentration. Whether the per-strand refinement pays has not been measured, and it is cheap to test
  once the columns exist.
* **Open: does this also explain the capture-ON residual?** At capture-ON the terminus concentration
  *inverts* (25 % of error mass on terminus edges, 69 % on junction). Something else dominates there and it
  has not been diagnosed.

## 10. RECOMMENDED ORDER

```
1. restore the v6 columns + version bump  ─ gate: bit-identical 32/32          (mechanical, ~1 session)
2. thread them to node_geometry (option B) ─ gate: bit-identical 32/32          (the real plumbing)
3. measure the bit on the REAL cfRNA index — how many real termini does it see? (before trusting the toy)
4. C1: the reframe at terminus seams       ─ the pre-registered A/B of §8       (the prize: −0.0178 capOFF)
5. C2: omega_graft per class               ─ re-fit, per structural class
6. C3 / C4, then goldens LAST
```

Steps 1–2 are low-risk and independently verifiable. **Step 3 is the one that could invalidate the plan**,
and it is cheap — do it before step 4.
