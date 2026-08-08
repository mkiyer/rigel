# Accumulator + prior — WHERE THIS CAMPAIGN IS, AND WHAT TO DO NEXT

    ⚠ **A DEV DOC.** Provisional; nothing may cite it. When an item settles, MOVE it to its permanent
    home and delete it here in the same edit.

    Opened 2026-08-07. **Rewritten 2026-08-08 as the RESUME POINT** — the design sections that used to
    live here are implemented, and their content moved out: derivations to `EQUATIONS.md` §3b/§3c,
    rulings to `DESIGN.md` §3.1, lessons to `TRAPS.md`. What is left is state, sequence and open
    questions. ⛔ If this file starts growing design again, move it out.

---

## 1. ⭐⭐⭐ THE TREE, EXACTLY

**Committed mid-development, deliberately, with the suite RED. That is not the standing baseline.**

| | |
|---|---|
| suite | ⛔ **40 failed / 3,197 passed** / 2 skipped / 9 xfail |
| the 40 | `test_golden_output` **21**, `test_prior_units` **13**, `test_priors` **6** — *all* consequences of the `assemble_priors` rewrite, nothing else |
| `payload_schema_digest` | **`19ee4ba867ff0441`** — FINAL for this campaign |
| caches | ✅ all **176** oracle caches (ladder 36 × 4, both flgap 4 × 4) + `pilot/scan_cache` rebuilt on that digest, verified 0 stale |
| lint | ✅ clean |

⛔ **The 40 are a known, understood, deliberate state — not a mystery to debug.** See §3 step 1.

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

**What phase 0 measured** (gDNA arm, `rel`):

| panel | ④ `O−F` | ⑧ `S−F` | ⑨ `O−S` | the share's portion |
|---|---|---|---|---|
| `flgap_long` (gDNA +40 %) | 0.0319 | **0.0229** | 0.0088 | **28 %** |
| `flgap_short` (gDNA −41 %) | 0.0067 | **0.0050** | 0.0027 | **25 %** |

⭐ Two earlier claims were **corrected by measurement**: the pooled share's capture contribution is
**0.36 %**, not the modelled −2.49 % (the uniform-placement model fails under capture); and
`gdna_eff_len`'s clamp is **1.03–1.25×**, not the predicted ~3× — with `nodes only` at 0.70–0.84×, so
the edge term is partly compensating a real deficit rather than purely inflating.

---

## 3. ⭐⭐⭐ WHAT TO DO NEXT, IN ORDER

**1 · GET TO GREEN.** The only tree-blocking task.

* `test_prior_units` (13) + `test_priors` (6): their target is the OLD rule's `ρ·span`. The new rule
  produces the **conserved fragment count**, which on a finite reference of span `S` is `ρ·(S−w+1)` —
  the number of fragments that FIT. `ρ·span` counted `w−1` start positions no fragment can occupy.
  ⭐ `test_prior_units`'s fixture is already made self-consistent (masses AND share enumerated through
  the specification) and **partition invariance already passes** — 1,001 on every tiling. What remains
  is re-basing the *target* and the ratio-flatness sweep.
* `test_golden_output` (21): the prior genuinely changed. ⛔ Regenerate **last**, after step 2, because
  step 4 may move it again — regenerating twice is wasted.

**2 · FIX THE pmf ESTIMATOR.** ⛔ **Prerequisite for everything downstream.** The gDNA pools are
geometry-limited and the RNA pool splice-limited — *opposite* censorings — so the fitted `μ_g − μ_r`
carries a phantom **+5 to +18 bp** gap that is not in the library. Both `length_likelihood` and any
analytic share consume these pmfs; shipping either on top is the net-negative trap in
`fragment_length_bias.md` §0. Work: pool de-tilt by **membership probability**
(`TRAPS: divide-by-a-probability`), and the shared EB anchor (`POOL_EB_PRIOR_ESS = 1000.0`,
`fl.py:86`, sits directly on the composition determinant).

**3 · FIND THE 72 %.** ⭐ **The largest unexplained term, and the highest-value question.** Perfect
per-component shares remove only **25–28 %** of the assembler's residual; the rest survives, and is **5×
worse when gDNA is LONGER** (`S−F` 0.030 vs 0.006 at capture-ON). That asymmetry is not the share.
⭐ Prime suspect: the **contained/crossing ROUTING** — a longer component deposits less as contained and
more as crossings, so *which bank* it lands in is itself gap-dependent, and the fraction of its prior
that passes through the share at all moves with the gap.

**4 · THE PER-COMPONENT SHARE** — re-ranked DOWN by phase 0. Prefer the hybrid: empirical scale `M(e)`
from the accumulator, analytic ratio `r = share_g/share_r`. Gate on the flgap pair.

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
3. ⚠ **Every cached payload is `drain: null` while production always drains.** The bound on `μ_r` is
   **[−3.96 %, −0.60 %]** — comparable to the phantom gap itself. Until that is a number, every FL
   measurement in this campaign carries a caveat its own size.
4. ⚠ **Unverified:** a junction inside an unsequenced mate gap may be filed as a *contiguous* unspliced
   crossing, tilting the unspliced RNA population long exactly beside junctions. Needs checking in
   `build_fragment`. If true it is a length-driven composition bias upstream of all of the above.
5. ⚠ **Debt:** `_density_times_span` in `priors.py` is now **dead** (no call sites), and
   `assemble_priors`' docstring still teaches the retired `ρ_c = Σm/ΣS; prior = ρ·span_bp` rule at
   `priors.py:293-297`. Clear both in step 1.

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
