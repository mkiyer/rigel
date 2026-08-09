# Accumulator + prior — WHERE THIS CAMPAIGN IS, AND WHAT TO DO NEXT

    ⚠ **A DEV DOC.** Provisional; nothing may cite it. When an item settles, MOVE it to its permanent
    home and delete it here in the same edit.

    Opened 2026-08-07. **Rewritten 2026-08-08 as the RESUME POINT** — the design sections that used to
    live here are implemented, and their content moved out: derivations to `EQUATIONS.md` §3b/§3c,
    rulings to `DESIGN.md` §3.1, lessons to `TRAPS.md`. What is left is state, sequence and open
    questions. ⛔ If this file starts growing design again, move it out.

---

## 1. ⭐⭐⭐ THE TREE, EXACTLY

**GREEN as of 2026-08-08. §3 step 1 is done; the red state below is history, kept only so the diff reads.**

| | |
|---|---|
| suite | ✅ **0 failed / 3,240 passed** / 2 skipped / 9 xfail — the standing baseline again |
| was | ⛔ 40 failed / 3,197 passed: `test_golden_output` 21, `test_prior_units` 13, `test_priors` 6 |
| net | +3 tests: `test_prior_units` 16 → 17, `test_priors` 19 → 21 (two gates the rewrite had left uncovered) |
| `payload_schema_digest` | **`19ee4ba867ff0441`** — FINAL for this campaign |
| caches | ✅ all **176** oracle caches (ladder 36 × 4, both flgap 4 × 4) + `pilot/scan_cache` rebuilt on that digest, verified 0 stale |
| lint | ✅ clean |

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

**`assemble_priors` — REWRITTEN, not yet green.** Sums conserved counts; no density, no `span_bp`
integration, no support-weighted pooling. Measured `O/F` **1.149 → 1.019** on three ladder gDNA levels;
capture-OFF `Σ|Δ|` −61 %.

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

**1 · GET TO GREEN.** ✅ **DONE 2026-08-08.** What it turned up, since none of it was a value re-base:

* ⛔ **Per-component partition invariance is BROKEN, and the earlier note here saying it "already
  passes" was reading the RNA side of a gate that fails first on the gDNA side.** Only the TOTAL is
  invariant — exactly, on all four tilings. Re-tiling 1200 bp from one node to twelve moves gDNA
  **−19.9 %** and RNA **+13.7 %** with the library unchanged. That is the pooled share, arriving in the
  unit tests at full size, and the tests now pin the biased value to 1e-9 rather than tolerate it.
* ⚠ **`EQUATIONS.md` §3b's mixture-independence was over-stated and is now scoped there.** The
  cancellation is exact PER LINE; contained mass never passes through the share at all, so at locus
  granularity a 900× mixture swing moves the distortion **0.578 → 0.767**. ⭐ **A correction derived from
  `share_r/share_g` and applied per locus would be wrong by the contained fraction** — this changes
  where step 4 has to act, not whether.
* ⛔ **`TRAPS: a-guard-outlives-its-divisor` — a new one, and it was live in this tree.** Removing the
  division left `_mass_where_there_is_opportunity` inert on the prior path while still load-bearing on
  the eff-length, where **nothing gated it**. Its perturbation test passed throughout and its assertion
  had silently reversed. Two gates now, naming each other.
* ⭐ Falsified: **7/7 injected defects caught** (`raw-incidence-sum`, `rescale-the-contained-term-too`,
  `share-clamped-to-one`, `retired-density-rule`, `drop-zero-opportunity-mass-from-the-COUNT`,
  `floor-the-eff-len-divisor`, `no-intergenic-rekey`). ⚠ `floor-the-eff-len-divisor` fires **exactly one**
  test — the one added here; before it, that defect was ungated.
* Goldens regenerated last, as planned. The prior's reach into the deliverable, measured across all 21
  scenarios: transcript `count` ≤ **3.0 %**, `tpm` ≤ **2.8 %**, `effective_length` ≤ **0.006 %**, and
  **`gdna_em_count` ≤ 33 %** — concentrated exactly where the prior arbitrates, which is the signature
  to expect.

**2 · FIX THE pmf ESTIMATOR.** ⛔ **Prerequisite for everything downstream.** Both `length_likelihood`
and any analytic share consume these pmfs; shipping either on top is the net-negative trap in
`fragment_length_bias.md` §0.

⛔⛔ **THE PRESCRIPTION THIS ENTRY USED TO CARRY — "pool de-tilt by membership probability" — IS ALREADY
IMPLEMENTED.** `detilt_pool` divides by `pi(w) = A(w)/T(w)` and is wired for both pools
(`fl.py:288-294`); `EQUATIONS.md` §4.1 is the derivation. ⚠ Anyone reading the old line would have
re-written working code. **Measured 2026-08-08** on the flgap pair, both capture arms, against the
origin-split payloads' own deposited truth:

| | fitted − true | without the de-tilt |
|---|---|---|
| phantom gap, 8/8 conditions | **+4.7 … +14.9 bp**, fixed sign | −9.4 … −28.4 bp |
| gDNA, capture OFF | **+0.1 / +0.0 bp** — exact to 5 s.f. | −1.0 / −0.2 |
| gDNA, capture ON | **+13.7** (μ_g 330) · **+3.5** (μ_g 120) | −19.3 / −12.1 |
| RNA, OFF / ON | **−5.0** / −1.2 | +9.2 / +9.0 |

⭐ **So the audit's magnitude was right and its mechanism was one defect where there are two**, with
different substrates and different repairs:

⛔⛔ **AND THEN THE DRAIN DISSOLVED HALF OF IT.** Every number above is on `drain: null` payloads while
production always drains, so both sides were re-measured DRAINED through the shipped route
(`_drain_side_buffer(_lift=)` → `lift_choices` → `second_pass.drain` on each origin partition, so the
fit and the truth move together and nothing is compared across a drain boundary):

| | undrained | **DRAINED** |
|---|---|---|
| **RNA** mean err, capture OFF | −4.9 / −5.1 | **−0.5 / −0.6** |
| **RNA** mean err, capture ON | −1.2 / −1.3 | **+2.7 / +2.5** (sign flips) |
| **RNA** shape `fit/true` hi÷lo across w = 200 | 0.907 / 0.909 | **0.990 / 0.992** (OFF) · 1.050 (ON) |
| **gDNA** mean err, ON | +13.7 / +3.5 | **+13.6 / +3.5** — Δ ≤ 0.1 bp |

* **2a · MOSTLY A DRAIN ARTIFACT, AND THE MECHANISM WAS ALREADY HANDLED.** ~90 % of the RNA error off
  capture was the undrained payload. ⭐ It is the same mechanism `_drain_side_buffer`'s docstring
  states — the held fragments are *"systematically the **long** ones, because a longer gap admits more
  hypotheses"* — so undrained, the RNA pool is missing its own long tail and reads short. ⛔ **The
  "sequenced-window opportunity" hypothesis was therefore a re-derivation of a problem the SECOND PASS
  already solves**, and the uniform-θ term it needed was explaining a plateau (1.049) that drains to
  1.006. Neither was built. ⚠ What survives is **+2.5 bp / 5 % under capture only**, sign-flipped by the
  drain — real, small, and probably not separable from 2b.
* **2b · THE gDNA DIVISOR IS CORRECT AND THE PLACEMENT IS NOT — and it is now the whole defect.** Exact
  off capture on both panels; the entire residual is under capture, scales with the gDNA mean, and is
  **untouched by the drain**. Its shape signature is far larger than anything on the RNA side: drained
  per-bin `fit/true` runs **1.22 … 4.18** in the tail under capture against 0.999–1.001 off it. That is
  `EQUATIONS.md` §4.4 verbatim — *"a placement model, not a better divisor"* — so ⛔ it must not be
  attacked by editing `gdna_opportunity`.

⛔ **De-ranked: `POOL_EB_PRIOR_ESS`.** Ablating 1000 → 0 moves every fitted mean by **≤ 0.1 bp**; pool
totals are 0.70–4.75 M so the shrink weight is ≤ 0.14 %. It is inert here. ⚠ It is *not* thereby
harmless on sparse real data — the anchor is the MIXTURE, so it shrinks the very gap `EQUATIONS.md` §3c
says is the only θ-independent evidence — but that is a separate claim needing a separate measurement.

⚠ **Two things the drained arm exposed about the ORACLE, which cost nothing to report and would have
cost a lot to assume:** the lift's ambiguity is **1.24–2.17 %**, not the ~0 that
`TRAPS: an-equal-length-panel-defeats-the-lift` predicts for a distinct-span substrate — the panels are
distinct BETWEEN components, not within one. And flgap_long's drained gDNA partition leaks
**1,491 (OFF) / 8,641 (ON)** spliced records where flgap_short leaks **1**; gDNA cannot splice, so
**flgap_short is the admissible panel for any drained measurement** and flgap_long is a cross-check
carrying a stated error bar.

**3 · FIND THE RESIDUAL — 57–85 %, NOT A FLAT 72 %.** ⭐ **The largest unexplained term, and the highest
-value question.** Perfect per-component shares remove **14.7 % … 42.6 %** of the assembler's residual
depending on stratum (§2), so the rest is 57–85 % and the "72 %" was a panel average of four numbers
that differ by 3×. ⭐⭐ **Attack `flgap_long` capture ON**: it carries the largest error by far
(`O−F` **0.0403**, `S−F` **0.0302**) and the share explains the *least* of it there (25.0 %).

⭐ **The "5× worse when gDNA is LONGER" asymmetry is CONFIRMED and is not the share**: at capture ON,
`S−F` is **0.0302** (long) against **0.0058** (short) = **5.2×**, measured after perfect per-component
shares are already applied. ⚠ And it is a capture-ON statement — off capture the same ratio is 2.4×.

### ⭐⭐⭐ FOUND, 2026-08-08 — `edge_owner_nodes` CREDITS A WHOLE LINE TO ONE FLANK, AND HALF OF IT IS OUTSIDE

**A contiguous line's conserved mass belongs to its TWO flanks**, in the ratio the deposit rule's own
opportunity function already states (`EQUATIONS.md` §3b):

    A_mass(w; a, b) = [ min(w−1, a) + min(w−1, b) ] / 2      ⇒  far share = E[min(w−1,a_far)]
                                                                            ────────────────────────────
                                                                            E[min(w−1,a_far)] + E[min(w−1,a_own)]

⛔ **`edge_owner_nodes` assigns the ENTIRE line to ONE node.** At a locus's outer lines the far flank is
intergenic — **median 8,066 bp against the locus flank's 211 bp**, longer on **95 %** of them — so
**~61 %** of that line's mass geometrically sits outside the locus and *all* of it is credited inside.
`_project_regions_to_loci` then drops nothing, and the locus's gDNA prior is over-called.

**The evidence, on all four flgap conditions:**

| | long/ON | long/OFF | short/ON | short/OFF |
|---|---|---|---|---|
| net `Σ(S−F)` | +141,783 | +19,914 | +27,762 | +8,303 |
| net ÷ abs | 0.998 | **1.000** | 0.998 | **1.000** |
| net ÷ (½ · boundary mass) | 1.296 | 1.244 | 1.352 | 1.136 |
| ⭐ net ÷ **exact `E[min(w−1,a)]` split** | 1.131 | **1.002** | 1.184 | **1.000** |
| corr(err, exact split) | +0.84 | +0.73 | +0.78 | +0.71 |
| corr(err, *interior*-line mass) — the null | +0.20 | −0.04 | +0.25 | −0.05 |
| corr(err, locus size bp) — the null | +0.02 | +0.01 | +0.08 | −0.02 |

⭐⭐ **Off capture the exact split predicts the whole over-call to 1.002 / 1.000, with NO free
parameter** — the flank lengths are `region_size_bp` and the expectation is over the TRUE gDNA histogram
off the origin-split payload, never the fitted pmf. ⭐ It also explains the **5× length asymmetry**: the
boundary mass scales with `E[w−1]`, and long/short boundary mass is 5.3× at capture ON against a measured
`S−F` ratio of 5.1×. ⛔ And `net ÷ abs ≈ 1.000` means nothing cancels — this is a pure systematic
over-call, not boundary redistribution.

⚠ **Under capture it closes to 1.13 / 1.18 and no further.** The residual is the uniform-placement
assumption inside the flank — `EQUATIONS.md` §4.4, the *same* cause as defect 2b, reached independently.

⭐ **THE REPAIR, and it deletes a heuristic rather than adding one.** Split each line's conserved mass
between its two flank nodes by `min(w−1,a) : min(w−1,b)` instead of assigning it whole. The intergenic
share is then dropped by the projection **correctly**, and ⛔ `edge_owner_nodes`' intergenic RE-KEY
becomes unnecessary — it exists because keying left "silently loses that line's crossing gDNA", and under
the split rule the locus keeps its own share by construction. One change removes a heuristic and a bias
together. ⚠ `_component_node_arrays` shares the same owner map, so the eff-len path moves with it.

⛔ **NOT YET IMPLEMENTED, and one question is genuinely open first:** `F` counts a fragment where its
FIRST BASE lands and `S` splits it by mass, so at a boundary they are two *definitions*, not a right and
a wrong. The defect above is real either way — the assembler uses conserved-mass semantics everywhere
else and then attributes a line whole, which is internally inconsistent — but **which target the EM
actually wants is an owner call.** Gate on the flgap pair, both capture arms.

⛔ **REFUTED: this is NOT defect 2b.** Arm S's `gdna_prior_count` is byte-identical under (a) a 3×
perturbation of every eff-length array — the pmf's only path into `assemble_priors` — and (b) a full
re-`calibrate` on a pmf tilted +15 bp. `gdna_eff_len` moved in both, so the perturbations ran. The two
defects are independent and must be repaired separately.

⚠ Superseded suspect: the contained/crossing **routing** was the prime candidate. The crossing term IS
where the error lives (relative error runs +1.7 % → +55 % across crossing-fraction quintiles), but the
cause is the ATTRIBUTION of a line to a node, not which bank a fragment routes into.

**4 · THE PER-COMPONENT SHARE** — re-ranked DOWN by phase 0. Prefer the hybrid: empirical scale `M(e)`
from the accumulator, analytic ratio `r = share_g/share_r`. Gate on the flgap pair.
⛔ **It must act PER LINE, inside `_conserved_count`, before the contained term is added** — step 1
measured the identity dissolving at locus granularity (`EQUATIONS.md` §3b), so a per-locus factor built
from the same ratio is wrong by the contained fraction, which on the unit fixture is 53 %.

**5 · `length_likelihood` ON**, scored on the library figure.

⛔ **Gate everything on the flgap PAIR, both capture arms — never the ladder alone.** The ladder's
realised gap is only +1.5–2.1 % and it is structurally blind to this whole axis.

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

1. ⛔ **The ladder's `O−S` is 0.0037, larger than flgap_short's 0.0027 at a 20× smaller gap.** Partly
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
