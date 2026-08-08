# Accumulator + prior — design and implementation plan

    ⚠ **A DEV DOC.** Provisional; nothing may cite it. When an item settles, MOVE it to its permanent
    home and delete it here in the same edit.

    Opened 2026-08-07. ⛔ **Supersedes `prior_derivation.md` entirely** — that document accumulated
    superseded proposals and two retracted claims and should be deleted once this one is accepted.
    Everything below that came from it has been re-derived against the code, not copied.

    **Status: DESIGN SETTLED, NOTHING IMPLEMENTED.** All three open questions were answered by the
    owner on 2026-08-07 and are recorded in §8 with the evidence. Implementation is §7.

---

## 1. Goals

| # | goal | what it needs from the accumulator |
|---|---|---|
| 1 | **strand model** | DISCRETE integer counts, on BOTH genome-strand columns |
| 2 | **length likelihood** (built, gated off) | a per-slot length moment pair — `inv_length_sum`, `length_sum` |
| 3 | **spliced fragments deposited correctly** | correct arithmetic for multi-junction / multi-node paths |
| 4 | **density deconvolution** | a per-NODE density (intergenic vs intronic, to peel RNA off introns) |
| 5 | **the prior** (this campaign) | a FRACTIONAL mass that sums to 1 per fragment → conservation |

⭐ **The organising principle, and it is what the current code violates: a COUNT and a MASS are two
different deposits and one number cannot be both.** A count is extensive and discrete (goal 1 needs
integers for a Beta-Binomial). A conserved mass is fractional and must sum to 1 per fragment (goal 5).
The current accumulator deposits only counts and densities, and `assemble_priors` then manufactures a
count from a density — which is the defect measured in §3.

---

## 2. ⭐⭐ THE PRODUCTION ACCUMULATOR ALREADY SOLVED THIS — read it before designing

`v0.7.1`, `src/rigel/native/calibration/accumulator.h`:

```cpp
struct Region   { uint32_t contained[kNChannels]; };                    // 16 B
struct Boundary { float    mass_left[4],  mass_right[4];                // fractional, conserved
                  uint32_t flux_left[4],  flux_right[4]; };             // discrete, for strand
```

⭐ **It carried the dual.** `flux` = discrete integer per side, NOT split (the strand model's
Beta-Binomial input). `mass` = float, split so the fragment sums to 1. Same object, two deposits,
two consumers. **This is the design to restore.**

The deposit rule (`accumulator.cpp`, crossing path):

```cpp
share = (slice.end - slice.start) * inv_L / n_cross;   // n_cross = crosses_left + crosses_right
if (crosses_right) { b[r+1].mass_left  += share; b[r+1].flux_left  += 1u; }
if (crosses_left)  { b[r  ].mass_right += share; b[r  ].flux_right += 1u; }
```

* **coverage-weighted** — a slice's mass is its share of the fragment's bases, `slice_len / L`
* **deposited at BOUNDARIES, never into the region** — edges keep ownership of crossing fragments
* **conserved** — `Σ slices (slice_len / L) = L / L = 1`, exactly
* a fully-traversed interior region splits 50/50 across its two boundaries (`n_cross == 2`)

⚠ **It is strictly better than the `1/K` share I proposed earlier**: `1/K` gives every crossed line an
equal share regardless of geometry; this weights by the bases actually adjacent to that line. Both
conserve; only this one is unbiased about *where* the fragment sat.

⚠ Production also carried a **left/right** split per boundary, which the new accumulator dropped. Its
purpose is recoverable from the new `edge_unspliced` / `edge_spliced` bank split — see §5.3 — but the
decision must be explicit, not inherited.

---

## 3. What is wrong today, measured

`assemble_priors` computes `rho_c = Σ mass_c / Σ S_c` then `a_c = rho_c · span_bp`. With a **perfect**
deconvolution (oracle masses in), `g50 ss0.99`:

| | capture OFF | capture ON |
|---|---|---|
| node-contained `rho` vs truth | 1.001 | **0.633** |
| edge-crossing `rho` vs truth | 0.977 | **7.621** |
| pooled | 0.999 | **1.247** |
| edge share of pooled mass / support | 7.4 % / 7.5 % | **53.7 % / 8.8 %** |

Capture puts gDNA on exons; a fragment on a short exon cannot be *contained*, so it deposits as a
*crossing*, K times for K lines. `Σm/ΣS` pools a K-inflated numerator against a flat per-line divisor.
The two halves nearly cancel, which is why it has read as roughly right. **+24.7 % genome-wide is the
+15.1 % `prior_vs_oracle` measures through the locus projection.**

⭐ Supporting fact: **40.2 % of gDNA fragments are crossings under capture, against 4.8 % off it.**

⛔ **The effective-length contraction is NOT the defect and must not be removed.** Capture changes the
sampling probability — it depletes introns and intergenic regions — so a locus's gDNA normalised by its
full span instead of its probed footprint understates gDNA and tilts the EM toward RNA. That is
`gdna_eff_len`'s job and it stays. What goes is the use of that machinery to *manufacture a count*.

---

## 4. ⭐⭐⭐ A FINDING THAT CHANGES GOAL 3 — verified against the deposit code

**A typical spliced fragment deposits on NO node and NO edge.** `contained` requires `not sj_ids`;
crossings are lines *strictly* inside a block, and a spliced block begins/ends exactly at an exon edge,
which is a cut — `searchsorted(..., 'right')` / `(..., 'left')` exclude it. Verified:

| path | lines crossed |
|---|---|
| spliced, 1 junction, blocks inside their exons | **none** |
| spliced, 2 junctions, 3 blocks | **none** |
| unspliced across one boundary | 1 |
| unspliced spanning a whole node | 2 |

⭐ **This is CORRECT, and it is why the design closes.** The prior arbitrates only the **unspliced**
competition — spliced fragments are certified RNA (gDNA cannot splice) and the EM assigns them
directly. So the conservation target is the unspliced population, and for an unspliced fragment
*contained ∪ crossing* is **exhaustive**: a contiguous interval either stays inside one node or crosses
≥ 1 boundary.

⚠ It also explains two audit results that otherwise look like bugs: `edge_spliced_*`'s density channels
are dead because that bank only ever catches the rare spliced block that spans an interior boundary,
and `CalibrationResult` already documents `mass_rna_spliced_edge` as "routinely two orders of magnitude
smaller than the junction flux at the same place".

⛔⛔ **BUT DO NOT OVER-GENERALISE THIS — AN EARLIER DRAFT DID.** The table above was produced on a toy
with a SINGLE isoform and no staggered exon boundary, which *by construction* cannot produce a spliced
contiguous crossing. With staggered isoforms a spliced fragment routinely crosses a contiguous line
(§5.1), and on the panel that population is **larger than the unspliced one at mRNA edges**.
⚠ `test_pass0_vs_oracle`'s fixture already documents the stagger as load-bearing for exactly this bank.

⭐ So goal 3 is: *a spliced fragment's junctions are accounted on the junction axis, its contiguous
crossings on the SPLICED edge bank, and neither enters the unspliced competition.* All three channels
already route that way; the new `mass` channel must route identically.

---

## 5. The design

### 5.1 Deposits to ADD

    edge_unspliced_mass[l] += share       // fixed-point; share = (slice_len / w) / n_cross
    edge_spliced_mass[l]   += share       // SAME rule, routed by the same `spliced` flag

⭐⭐ **TWO new banks, and no left/right split** (§8 Q2). Nodes need nothing new:
`node_contained_count` is already 1 per contained fragment, i.e. already the conserved node mass.

⛔⛔ **`edge_spliced_mass` IS REQUIRED — an earlier draft of this plan dropped it and was wrong**
(owner, 2026-08-07). A spliced fragment CAN cross a contiguous line, whenever another isoform's exon
boundary falls INSIDE one of its blocks:

    TA+ exons (1000,2000) (9000,10000)      a TA+ fragment splicing 2000->9000 and running to 9200
    TB+ exons (1000,2000) (9050,10000)      crosses the line at 9050 CONTIGUOUSLY

⭐ **And it is not a corner case — measured on the panel it is the DOMINANT edge population for mRNA:**

| | gDNA | mRNA |
|---|---|---|
| `Σ edge_spliced`, capture OFF | **0** | **1,593,441** on 5,992 lines (unspliced: 1,485,580) |
| `Σ edge_spliced`, capture ON | **0** | 1,146,842 on 3,959 lines (unspliced: 1,270,161) |

Exactly zero on gDNA and larger than the unspliced crossings on mRNA — the certified-RNA signature.
**These are proven pure RNA at that line and must never enter the unspliced competition.**

⭐ **THE STRUCTURAL ARGUMENT, WHICH IS THE STRONGEST ONE.** At deposit time every edge channel is
selected by ONE tuple — `(edge_count, edge_inv_length, edge_length) = spliced ? spliced : unspliced`.
Adding `mass` to that tuple routes it correctly for free. **Omitting it would make `mass` the only
channel that does not honour the spliced/unspliced split**, which is precisely the latent bug this
guards against: a spliced fragment's mass silently landing in the unspliced bank.

⚠ **THE TWO BANKS HAVE DIFFERENT CONTRACTS AND A CONSUMER MUST NOT CONFUSE THEM:**

| bank | sums to, per fragment | is it a conservation ledger? |
|---|---|---|
| `edge_unspliced_mass` | **1** over the lines an unspliced fragment crosses | ✅ **yes** — see the proof below |
| `edge_spliced_mass` | `crossed_block_len / w` — a PARTIAL | ⛔ **no**, and it is not meant to be |

A spliced fragment's blocks that contain no interior boundary deposit nothing (their accounting is on
the junction axis), so `edge_spliced_mass` is a partial by construction. ⭐ **That is correct, because
it is a per-LINE certified-RNA term, not a per-FRAGMENT count** — and per line it is exactly
commensurate with the unspliced mass, both being "the share of this fragment's bases adjacent to this
line", which is what makes the two safe to compare at one edge. ⛔ Do not read it as "the number of
spliced fragments at this line".

**Conservation, proved on the population that matters.** An unspliced path is ONE segment. If it stays
inside a node it is `contained` (1 deposit). Otherwise it has `n ≥ 2` slices: slice 0 crosses right only
(`n_cross = 1`), slice `n−1` crosses left only (`n_cross = 1`), every interior slice crosses both
(`n_cross = 2`). **Every slice has `n_cross ≥ 1`**, so `Σ share = Σ slice_len / w = 1` exactly, with
nothing dropped by the `n_cross == 0` guard.

⚠ **Fixed-point, not float.** `TRAPS: integer-channels-reproduce` — float channels are not bit-identical
across worker counts, and the whole parity gate against `_accumulator_reference.py` rests on integer
merges. Use the existing `INV_LENGTH_SCALE` convention.

### 5.2 Deposits to REMOVE — no production consumer

Traced through `CalibrationSubstrate` → `PopulationView` including the dynamic `getattr` in
`node_geometry._channel`. Each appears only in the schema, the C++, the substrate wrapper, and parity
fixtures.

| bank | verdict |
|---|---|
| `edge_spliced_inv_length_sum` | ⭐ REMOVE — `length_likelihood` reads `node_contained` + `edge_unspliced` only |
| `edge_spliced_length_sum` | ⭐ REMOVE — redundant with `pool_lengths` |
| `sj_length_sum` | ⭐ REMOVE |
| `node_spanning_count` / `_inv_length_sum` / `_length_sum` | ⭐ **REMOVE — decided on evidence, §8 Q3.** No consumer; the mass is already at the edges the fragment crosses; and the one argument for keeping it (spanning is the only per-node observable at short exons where `contained` is empty by geometry) **fails when measured** |

### 5.3 KEEP, with the reason

| bank | why |
|---|---|
| `node_contained_count` | goals 1, 4, 5 — the node population the solve deconvolves, the strand model's node counts (both columns), AND the conserved node mass |
| `node_contained_inv_length_sum`, `node_contained_length_sum` | goal 2 — dormant by design; the designed input to the only channel that can give an AMBIG slot its own composition evidence |
| `edge_unspliced_count` | goals 1, 4 — `strand_deconv` reads both columns |
| `edge_unspliced_inv_length_sum` | ⭐ **LIVE** in `second_pass.py:461` (the drain scores held fragments) + goal 2 |
| `edge_unspliced_length_sum` | goal 2 |
| `edge_spliced_count` | `calibrate` (`mass_rna_spliced_edge`), `node_geometry`, `sweep` |
| `sj_count` | `calibrate` (`mass_rna_junction`, flux), `node_geometry` |
| `sj_inv_length_sum` | ⭐ **LIVE** in `second_pass.py:451` |
| `node_start_count` | the conservation invariant `Σ == qc.deposited` — the only externally checkable statement the accumulator makes |
| `pool_lengths`, `deposited_lengths` | `fl.py` — the five pure pools and the EB anchor |

⭐ **Net: two banks added, six removed** — `edge_spliced_inv_length_sum`, `edge_spliced_length_sum`,
`sj_length_sum`, and the three `node_spanning_*`.

### 5.4 Collapse the unused strand splits (owner agreed)

`_build_length_loglik` calls `.sum(axis=1)` on all three inputs and `length_likelihood` is explicitly
strand-agnostic — which strand a read aligned to says nothing about whether the molecule was gDNA or
RNA. **The `[n,2]` shape on every surviving `*_inv_length_sum` / `*_length_sum` is half wasted by
construction**; collapsing to `[n]` touches no consumer.
⛔ **The `count` banks keep both columns — that is goal 1.** `strand_deconv` reads
`edge_unspliced.count` per column and `gdna_strand` reads `node_contained.count` per column, and a
Beta-Binomial needs discrete integers on both strands.

Cost at rest today: **84 B/node** (40 contained + 40 spanning + 4 start), **80 B/edge**. Scale by the
INDEX (~1.5 M nodes genome-wide), not the panel — calibration's cost is depth-independent.

### 5.5 The prior

`mass_gdna_node[r]` is already `f_g(r) · contained_count[r]`, so the node term needs no arithmetic. The
edge term rescales from the K-inflated count onto the conserved mass:

    a_g(locus) = Σ_{r∈locus} share(r,locus) · mass_gdna_node[r]
               + Σ_{l∈locus} share(l,locus) · mass_gdna_edge_unspliced[l] · (edge_unspliced_mass[l] / edge_unspliced_count[l])

    a_r(locus) = the same on the RNA masses (spliced withheld, as today)

* `a_g + a_r` = the locus's unspliced fragments. Exactly the population the EM arbitrates.
* Each object's composition applies to its OWN population — a line's composition is never used for a
  node, which matters because **mature RNA never crosses an exon↔intron seam (0 of 1,146)**, so the two
  are different populations by construction.
* `gdna_eff_len` unchanged.

---

## 6. Consumer impact

| consumer | reads | impact |
|---|---|---|
| `gdna_strand.py:375` | `node_contained.count` (2 col) | none |
| `strand_deconv.py:86` | `edge_unspliced.count` (2 col) | none |
| `density_deconv.py:132` | `node_contained.count` | none — goal 4 is node-only, already supported |
| `background_reference.py:92` | `node_contained.count` | none |
| `sweep.py:634,668,669` | contained + both edge counts | none |
| `node_geometry.py:219-224,298` | contained, edge_unspliced, edge_spliced, junction | none |
| `calibrate.py:652,659,685` | `edge_spliced.count`, `junction.count` | none |
| `second_pass.py:451,461` | `sj_inv_length_sum`, `edge_unspliced_inv_length_sum` | none |
| `fl.py` | `pool_lengths`, `deposited_lengths` | none |
| **`priors.py`** | the four eff-len arrays + masses | ⭐ **rewritten** — §5.5 |

⭐ **Every calibration consumer is untouched.** The change is additive in the accumulator and confined
to `assemble_priors` on the consuming side. That is what makes it a tractable commit.

---

## 7. Implementation order — each step independently gated

1. ⭐⭐ **PROTOTYPE THE WHOLE CHANGE IN PYTHON FIRST — no production edit.** Compute
   `edge_unspliced_mass` from `_accumulator_reference.py` on one capture-ON and one capture-OFF
   condition, apply §5.5, and score with `prior_vs_oracle.py` against `F`. **This prices the entire
   change before a line of C++.** Target: the **+15.1 %** under capture → ~0, and capture-OFF's 0.999
   unchanged. ⛔ `TRAPS: measure-the-ceiling-first` — if it does not land, nothing downstream is worth
   building. It also answers §8 Q4 for free.
2. **Add the bank to `_accumulator_reference.py`** (the executable spec) + its gates.
3. **Add it to the C++** and restore byte-identity parity with the reference.
4. **Bump `payload_schema_digest`** — every scan cache invalidates by design. ⚠ The 341 M ladder oracle
   cache and both flgap caches must be rebuilt; budget for it.
5. **Rewrite `assemble_priors`** per §5.5. Gate with `prior_vs_oracle.py` (O vs F) and
   `arm_identity.py`.
6. **Remove the dead banks** (§5.2) — separate commit, after the above is green, so a regression is
   attributable to one change.
7. **Re-run the campaign**: `prior_vs_oracle.py`, then `quant_accuracy.py --arm base/noop/oracle`.
   ⚠ ROADMAP §1.3's ceiling was measured with the OLD assembler and must be re-recorded
   (`TRAPS: re-record-the-baseline`).

---

## 8. ✅ DECISIONS — all three answered (owner, 2026-08-07)

**Q1 — Is `contained ∪ crossing` exhaustive for UNSPLICED fragments?** ✅ **Sufficient for our
purposes.** Unannotated spliced fragments are a tiny fraction and are DEFERRED, consistent with the
other classes already held out of the accumulator (chimeras, multimappers). ⚠ Record the deferral
rather than the assumption: conservation is exact over *deposited, unspliced, annotated* fragments, and
the excluded classes are named.

**Q2 — Does an edge need a left/right split?** ✅ **NO.** Production carried `mass_left`/`mass_right`
and `flux_left`/`flux_right` and **never used the split — every consumer summed the two sides.** One
number per edge. This simplifies both the deposit and §5.5's projection, and it is why the new
accumulator's single-number-per-line shape is right rather than a regression.

**Q3 — Keep or delete `node_spanning_*`?** ✅ **DELETE — decided on evidence, not argument.** The mass
is not lost (a spanning fragment crosses ≥ 2 lines and is already deposited there), and the only
substantive case for keeping it — that `spanning` is the sole per-node observable at short exons, where
`contained` is empty by geometry, which is exactly ROADMAP's evidence-starved AMBIG slot — **fails when
measured**:

| in-locus nodes with NO contained evidence | capture OFF | capture ON |
|---|---|---|
| starved nodes | 13,544 | 19,209 |
| …reached by `spanning` | 13,352 | 9,774 |
| …reached by a bounding **edge** | 13,495 | **13,583** |
| ⭐ **reachable ONLY by `spanning`** | **0** | **141 nodes / 822 fragments** |

⭐ Under capture the **edges reach strictly more starved nodes than spanning does**, and the exclusive
case is 822 fragments of ~10 M — **0.008 %**. Spanning carries no evidence the edges do not.
⚠ Reversal is not free: re-adding a deposit is mechanical, but it costs another
`payload_schema_digest` bump and a full cache rebuild. Deleting now, while this change is already
paying that cost, is the cheap moment.

⚠ **Q4, still open and lower stakes** — is `f_g(l)` unbiased for a crossing fragment's gDNA
probability? Same population by construction, so it should be, but it is an assumption. Step 1
measures it against the oracle at no cost, before any C++.

---

## 9. Risks

* ⛔ **The schema digest invalidates every cached scan.** ~400 MB of oracle cache across three panels.
  This is by design and is the single largest time cost of the change.
* ⛔ **The C++/reference byte-identity gate** is the project's strongest invariant. Steps 3 and 4 must
  land together or the suite is red in between.
* ⚠ **`mass_rna_edge` is spliced-INCLUSIVE** by an existing convention (`mass_gdna_edge + mass_rna_edge
  == unspliced + spliced` per edge). §5.5 must re-read that against the new mass, or the RNA arm
  double-counts the certified fraction.
* ⚠ **Float vs fixed point** — see §5.1. Production used `float` and got away with it; this tree gates
  on bit-identical worker merges and must not.
