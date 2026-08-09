# Accumulator + prior — WHERE THIS CAMPAIGN IS, AND WHAT TO DO NEXT

    ⚠ **A DEV DOC.** Provisional; nothing may cite it. When an item settles, MOVE it to its permanent
    home and delete it here in the same edit.

    Opened 2026-08-07. Rewritten 2026-08-08 as the RESUME POINT, twice. The design that used to live
    here is implemented and its content has moved out: derivations to `EQUATIONS.md` §3b/§3c, rulings to
    `DESIGN.md` §3.1 and §3.1b, lessons to `TRAPS.md`. What is left is state, sequence and open
    questions. ⛔ If this file starts growing design again, move it out.

---

## 0. ⭐⭐⭐ START HERE — THE NEXT SESSION IN FIVE LINES

1. Tree is **GREEN at `411999d5`**: 0 failed / 3,253 passed / 2 skipped / 9 xfail, lint clean.
2. **`python scripts/design/flgap_study_cache.py --list`** — four conditions cached, ~1 s to load.
   `priors.py` is outside the cache key, so assembler changes are a one-second loop.
3. ⛔⛔ **DO NOT chase the assembler's residual until §3 step 1 lands.** The yardstick it was measured
   against (`F`) is disqualified, so `O−F` and `S−F` are currently uninterpretable.
4. ⭐ **The assembler was never the main term.** `P−F` is **0.82–0.97** under capture against `O−F`
   **0.008–0.04**. The error is upstream, in calibration.
5. Two length biases remain open and are named in §3 steps 3–4.

---

## 1. ⭐⭐⭐ THE TREE, EXACTLY

| | |
|---|---|
| HEAD | **`411999d5`** — "an edge owns its own crossings" |
| suite | ✅ **0 failed / 3,253 passed** / 2 skipped / 9 xfail · lint clean |
| `payload_schema_digest` | **`19ee4ba867ff0441`** — unchanged by the edge-ownership work, and FINAL unless §4 says otherwise |
| caches | ✅ 176 oracle caches + `pilot/scan_cache` on that digest · ✅ the 4-condition flgap **study cache** |
| switches | ⛔ `message_propagation = False`, `length_likelihood = False` — both deliberate, both study-only |

⚠ **The 9 xfails are untouched.** All five strict ones still fail as intended.

---

## 2. WHAT LANDED

**Accumulator schema — DONE, and FINAL unless §4 says otherwise.**

* ⭐ Two conserved-mass banks (`edge_unspliced_mass`, `edge_spliced_mass`), one column each.
* ⛔ Six dead banks removed; five length moments collapsed `[n,2] → [n]`.
* Structs: `Node` 80 → **24 B**, `ContiguousEdge` 80 → **48 B**, `JunctionEdge` 40 → **16 B**.
* Byte-identity with the C++ restored; **16 new gates** in `tests/native/test_conserved_mass.py`.
* ⭐ Falsified: 6 injected defects, all caught. **An exact `1/K` deposit passes every conservation
  gate** and fails only the per-base re-derivation — conservation cannot pin this rule.
* ⭐ `sj_count` keeps both strand columns, now as a stated ruling (aligner-artifact detection).

**`assemble_priors` — REWRITTEN and GREEN.** Sums conserved counts; no density, no `span_bp`
integration, no support-weighted pooling. Measured `O/F` **1.149 → 1.019** on three ladder gDNA levels;
capture-OFF `Σ|Δ|` −61 %.

**⭐⭐ THE EDGE-OWNERSHIP RESTORE — `411999d5`, and it is a RULING, not a tuning.** A node owns contained
fragments, an edge owns crossings, nothing is re-attributed; a locus collects the edges touching its
nodes. `DESIGN.md` §3.1b is the ruling, `TRAPS: a-fold-grows-a-heuristic` the lesson.

* ⛔ DELETED: `edge_owner_nodes`, the intergenic re-key, `_component_node_arrays`,
  `_mass_where_there_is_opportunity`, `_left_keyed_edge_arrays`.
* ⭐ The zero-opportunity guard is now **structural** — `min(m/ρ_ref, S)` per object means `S = 0`
  contributes exactly 0, with no test and no `1e-9` floor.
* ⭐ The transcript path emits an EDGE index, so the locus path and `transcript_capture_eff_lengths`
  are the same operation over different object sets.
* ⭐ **Numerically identical on real data** (`S−F`/`O−F` unchanged on all four flgap conditions);
  2 of 21 goldens moved by ≤ 4.8e-4; `payload_schema_digest` unchanged.
* ⭐ **8/8 injected defects caught** — after the first pass exposed TWO gate holes: `node_right_edge[r]
  == r` on a single reference made every fixture blind to the edge-vs-node index distinction, and the
  per-object contraction clamp had no gate at all. Both closed.
* ⛔ **Contended edges are NOT impossible**: 20–34 per flgap condition at ~0.01 % of the mass.
  `contended_edges` reports them; an assert would have died on real data.

**Phase-0 instrumentation — DONE.** `prior_vs_oracle.py` gained arms **S** (`S_vs_F`, `O_vs_S`) and the
`gdna_eff_len` clamp diagnostic; `OracleTruth.component_shares()` measures the true per-component share
off the origin split; the two mass banks joined `_BANKS` so the split is validated on them.

**What phase 0 measured** (gDNA arm, `rel`) — ⚠ **each row is a MEAN over the panel's 4 conditions, which
the original entry did not say, and the per-stratum spread it hides is the interesting part.
Re-run and decomposed 2026-08-08:**

| panel | arm | ④ `O−F` | ⑧ `S−F` | ⑨ `O−S` | the share's portion |
|---|---|---|---|---|---|
| `flgap_long` (gDNA +40 %) | capture **OFF** | 0.0146 | 0.0084 | 0.0061 | **42.6 %** |
| | capture **ON** | **0.0403** | **0.0302** | 0.0102 | 25.0 % |
| | *recorded (4-cond mean)* | *0.0319* | *0.0229* | *0.0088* | *28 %* |
| `flgap_short` (gDNA −41 %) | capture **OFF** | 0.0041 | 0.0035 | 0.0009 | **14.7 %** |
| | capture **ON** | 0.0080 | 0.0058 | 0.0036 | 28.0 % |
| | *recorded (4-cond mean)* | *0.0067* | *0.0050* | *0.0027* | *25 %* |

⭐ **The recorded row reproduces** — my two-condition means run 0.83–0.93× of it, the residual being the
`ss_0.99` arms I did not run. ⛔ **But "the share is 25–28 %" is an artefact of the averaging: per
stratum it is 14.7 % … 42.6 %**, so the open residual is **57–85 %**, not a flat 72 %. And the stratum
where the share matters MOST (`flgap_long` capture OFF, 42.6 %) is *not* the one carrying the most error
(`flgap_long` capture ON, `O−F` 0.0403) — a single panel row cannot express that and should not be quoted.

⭐⭐ **THE DRAIN DOES NOT MOVE ANY OF IT.** Re-run with both sides drained through the shipped route,
every arm shifts by **≤ 0.0003** in `rel` and the share portion by **≤ 0.7 pp**. So the ranking survives
the caveat that broke the FL measurement — ⚠ which was worth checking precisely because it was not
predictable from the FL result: there the drain was worth +4.5 bp.

⚠ **Admissibility, stated not assumed.** Sum-to-full is EXACT on every drained partition. The gDNA
spliced leak is **1** record on `flgap_short` against **1,010 / 5,827** on `flgap_long`, so flgap_short is
the admissible panel and flgap_long's drained arm carries real contamination — it agrees to 3 decimals
anyway, which is why the conclusion holds on both.

⚠ **`O−F` against the pre-rewrite artifact on disk** (`arms/prior_vs_oracle_*.json`, 2026-08-07 16:55,
the retired density rule — and it has **no S arm at all**): capture-ON went **0.0884 → 0.0080**
(flgap_short) and **0.2639 → 0.0403** (flgap_long), an 11× and 6.5× win. ⛔ **But capture-OFF went the
other way: 0.0029 → 0.0041 and 0.0082 → 0.0146**, a 1.4–1.8× regression that the "1.149 → 1.019 /
−61 %" summary above does not mention. ⚠ Confounded — that artifact also predates the accumulator schema
change — so it is a flag to re-derive cleanly, not a verdict.

⭐ Two earlier claims were **corrected by measurement**: the pooled share's capture contribution is
**0.36 %**, not the modelled −2.49 % (the uniform-placement model fails under capture); and
`gdna_eff_len`'s clamp is **1.03–1.25×**, not the predicted ~3× — with `nodes only` at 0.70–0.84×, so
the edge term is partly compensating a real deficit rather than purely inflating.

---

## 3. ⭐⭐⭐ WHAT TO DO NEXT, IN ORDER

⛔ **Gate everything on the flgap PAIR, both capture arms — never the ladder alone.** The ladder's
realised gDNA/RNA length gap is +1.5–2.1 % and it is structurally blind to this axis. Use
`scripts/design/flgap_study_cache.py`.

---

### 1 · ⛔⛔ FIX THE YARDSTICK — nothing downstream is interpretable until this lands

**`prior_vs_oracle`'s `F` arm counts a gDNA fragment in whichever locus holds its FIRST BASE.** The EM
counts every fragment that is a *candidate* — any overlap — once, in the one multi-locus it belongs to
(`DESIGN.md` §3.1b). So `F` drops exactly the fragments that overlap a locus but start outside it, and
`O−F` / `S−F` measure against a target the EM does not use.

⭐ **The hypothesis to test, and it is cheap**: a straddling fragment deposits its conserved mass on
lines that all touch the locus's nodes, so `S` may already give it **exactly 1.0** — in which case the
assembler is correct and *the "72 % residual" does not exist*. ⚠ Not asserted; measured.

* Build `F_overlap`: each gDNA fragment counted once in the multi-locus whose objects it touches.
  ⛔ The owner's ruling that overlapping loci MERGE settles ambiguity *between* loci; it does not
  supply the start-base-vs-overlap correction, which is what this is.
* Re-score `O` and `S` against it on all four flgap conditions, drained.
* ⛔ Correct `prior_vs_oracle.py`'s docstring, which still asserts `F` is exact truth for the gDNA arm.

**Done when**: the assembler's residual is a number measured against the EM's own semantics — or is
demonstrated to be zero.

---

### 2 · ⭐⭐ GO AT CALIBRATION, NOT THE ASSEMBLER — this is where 95 % of the error is

Measured 2026-08-08 on the flgap pair, gDNA arm, `rel`:

| | `S−F` | `O−F` | **`P−F`** |
|---|---|---|---|
| long / capture ON | 0.0302 | 0.0406 | **0.8154** |
| short / capture ON | 0.0058 | 0.0081 | **0.9656** |
| long / capture OFF | 0.0084 | 0.0146 | 0.0272 |
| short / capture OFF | 0.0035 | 0.0041 | 0.0325 |

⭐⭐ **Under capture, `P−F` is 20–100× `O−F`.** Every prior-assembly change this campaign made — the
conserved count and the edge-ownership restore — moved the smaller term. ⚠ The `S`/`O`/`F` columns
inherit step 1's caveat; **`P−O` does not**, and it is the one that matters here.

⭐ **First arm, and it is one config flag**: `message_propagation = True`. It is OFF as a study
configuration, and `CLAUDE.md` records the price as **+154.8 % on unstranded × capture-ON** — which is
exactly the stratum above. Establish how much of `P−O` is the muted relay before designing anything.

---

### 3 · ⚠ DEFECT 2b — the gDNA pmf under capture, the surviving length bias

The fitted gDNA pmf reads **+13.6 bp long** at μ_g 330 and **+3.5 bp** at μ_g 120, with drained per-bin
`fit/true` running **1.22 … 4.18** in the tail. **Exact off capture** (+0.1 / +0.0 bp), and
**untouched by the drain** (Δ ≤ 0.1 bp).

⛔ **It is a PLACEMENT problem, not a divisor problem** — `EQUATIONS.md` §4.4: the opportunity
correction assumes uniform placement within a flank and capture concentrates gDNA at the probes. ⛔ Do
not attack it by editing `gdna_opportunity`; that divisor is measurably exact where its assumption holds.

---

### 4 · THE PER-COMPONENT SHARE — the other length bias

One `q` for two components tilts the g:r split. `O−S` measures it at **0.0009–0.0102**, i.e.
**14.7 %–42.6 %** of the assembler's residual depending on stratum (§2).

⛔ **It must act PER LINE, inside `assemble_priors`' crossing term, before the contained term is added.**
`EQUATIONS.md` §3b: the cancellation is exact per line and dissolves at locus granularity, so a
per-locus factor built from `share_r/share_g` is wrong by the contained fraction (53 % on the unit
fixture). ⚠ Prefer the hybrid — empirical scale from the accumulator, analytic ratio — and note that the
analytic ratio consumes the pmfs, so **step 3 gates it**.

---

### 5 · `length_likelihood` ON, scored on the library figure

---

## 4. DOES THE ACCUMULATOR CHANGE AGAIN?

⚠ **Possibly once more, and it would partly UNDO a removal above.** The candidate is
**`edge_spliced_inv_length_sum`** (8 B/edge, ~8.3 MB genome-wide). `edge_spliced_count /
edge_spliced_inv_length_sum` is a **local, model-free mean length of the certified-RNA population** —
gDNA cannot splice — which breaks the dependence on the globally-fitted, splice-censored `rna_pmf`.

⛔ **Do not land it on that argument.** Measure first, off the existing origin-split caches: (a) is the
*spliced* population's length representative of the *unspliced* RNA the deconvolution arbitrates? and
(b) how many lines carry enough spliced crossings for the ratio to be stable?

⛔⛔ **AND IF IT LANDS, BATCH IT.** Each schema change costs a digest bump plus a ~2 h rebuild of 176
caches plus the pilot. That has been spent once. There must be exactly **one** more, carrying every
remaining schema change at once.

⚠ The original removal was not wrong on its evidence — nothing read the bank. What changed is that a
consumer may now exist. Record it as a reversal with its reason, not as a mistake.

---

## 5. OPEN — do not treat these as settled

1. ⛔ **The ladder's `O−S` is 0.0037, larger than flgap_short's 0.0027 at a 20× smaller gap.**
   ⚠ **And every `O−F`/`S−F` figure in this doc now also inherits §3 step 1's yardstick caveat.** Partly
   explained (realised gap +1.5–2.1 %, and `share_c` is a censored functional sensitive to *shape*), but
   the ordering is not. **Do not quote the absolute `O−S` values until it is.** The ratios between
   panels are the trustworthy part. ⭐ The arm itself was checked and is sound: at the 11,341 lines where
   gDNA never crossed, its truth mass is exactly 0, so the `1.0` share default is inert on that arm.
2. ⚠ **The flgap panels vary the standard deviation as well as the mean** (1.2× and 2.0×). Every
   analysis so far is mean-only, and `share_c` is variance-sensitive even at equal means. The "±40 %
   gap" shorthand is not what the panels test.
3. ⛔⛔ **MEASURED 2026-08-08, AND IT IS NOT A CAVEAT — IT WAS THE DEFECT.** Every cached payload is
   `drain: null` while production always drains. On the FL models the drain is worth **+4.5 bp on `μ_r`
   off capture** and it removed ~90 % of what had been attributed to the junction opportunity (§3 step 2).
   The gDNA side moved **≤ 0.1 bp**. ⛔ **So the caveat is NOT uniformly small and must not be waived by
   its old bound**: it is large and RNA-only. ⚠ Every *other* FL and share measurement in this campaign
   is still undrained, including the phase-0 `O−F` / `S−F` decomposition — **re-run the ones that carry a
   decision before quoting them again.**
4. ✅ **RESOLVED 2026-08-08, in the negative.** A junction inside an unsequenced mate gap is **HELD in the
   side buffer, not silently filed as contiguous** — that is what the deferred bank is for, and draining
   deposits it. The RNA pool's long-tail deficit that looked like this (per-bin `fit/true` 0.907 across
   w = 200) drains to **0.990**. ⛔ Nothing needs checking in `build_fragment`. ⚠ It is only fully
   resolved where the deferral fires; a fragment whose gap admits exactly ONE hypothesis deposits in pass
   one, so this says the *bulk* mechanism is handled, not that the edge case is impossible.
5. ✅ **Debt — CLEARED 2026-08-08.** `_density_times_span` deleted; `assemble_priors`' docstring now
   teaches the conserved count and says why `ρ·span_bp` is retired. ⭐ Two more found while clearing it:
   the folded gDNA mass out of `_component_node_arrays` was **stored and never read** (it was the retired
   rule's numerator — now discarded explicitly, with the reason), and the docstring taught a
   **Laplace-smoothed `(2G+1)/span` IPR that is not in the code** — the mechanism is `min(m/ρ_ref, S)`
   per object plus the `w = C/(C+1)` contained-evidence shrinkage.

---

## 6. WHAT NOT TO RE-LITIGATE

* **The `S` arm is a DIAGNOSTIC CEILING, not a candidate implementation.** It feeds the truth
  per-component share, which production cannot have. `O−S` prices the pooled share; it is not a design.
* **A node owns contained fragments, an edge owns crossings, nothing is re-attributed.** `DESIGN.md`
  §3.1b is the ruling and `TRAPS: a-fold-grows-a-heuristic` the lesson. ⛔ The fold is not a hack this
  branch introduced — shipped v0.7.1 has it too, for the gDNA arm only.
* **`contended_edges` REPORTS and must never renormalise.** It fires on real data (20–34 edges,
  ~0.01 % of mass); making it an assertion would kill a production run.

* The conserved-mass rule is **coverage-weighted, not `1/K`** — both conserve; only one is expressible
  per base, and that is the only gate that separates them. `EQUATIONS.md` §3b.
* The prior's target is the **conserved fragment count**, not `ρ·span`. `ρ·span` was the approximation
  the density conversion happened to produce.
* The **counts** keep both strand columns; the **moments and mass** keep one. `DESIGN.md` §3.1.
* `sj_count` keeps both columns for **aligner-artifact detection**, not for the strand model.
* `edge_mass_per_crossing` is **geometry, not a deconvolved mass** — it must never join
  `prior_vs_oracle.OVERRIDE_FIELDS`.
* You may **not** make the tool gap-robust by shrinking the estimated gap: `μ_g − μ_r` is the only
  θ-independent composition evidence an AMBIG slot can get. `EQUATIONS.md` §3c.
