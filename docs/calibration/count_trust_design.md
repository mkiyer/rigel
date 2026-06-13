# Count-trust (β) — design proposal + phased implementation plan

**Status:** design + plan, 2026-06-13. Builds on
[`signal_integration_investigation.md`](signal_integration_investigation.md). The integration of the strand
and count signals happens inside `simplex.solve_node` by **adding log-likelihoods**; the only thing we
*design* is how much to **trust the count** relative to the strand. This doc proposes a **per-node count
penalty `β`** (a tempered/power likelihood) and a phased plan from a single hard-coded value to a derived,
calibrated function — explicitly retiring the old hard-coded `I₀` magic number along the way.

## 1. The lever and the key insight

`solve_node` maximises `log L = log L_strand + β · log L_count + log prior` over the pie. Two precisions
govern the strand-vs-count balance:

- **Strand precision** `I = N·(2κ−1)²` — *intrinsic* to the Beta-Binomial strand likelihood (Fisher
  information for `f_g`). It **vanishes at κ=½** (unstranded) and is large at κ=0.99 deep nodes. This is
  automatic; nothing to set.
- **Count penalty `β`** — how much of the count log-likelihood's curvature we keep. `β=1` full trust;
  `β→0` ignore. This is the knob.

**The required behaviour** (the user's key insight): the count must contribute **much less** than the
strand when the strand is excellent (κ=0.99) and **much more** (take over) when the strand is silent
(κ=½). Because `I → 0` as `κ → ½`, the *relative* trust already swings correctly with strand quality — at
κ=½ *any* positive count precision dominates the vanishing strand, and at κ=0.99 the large `I` dominates a
modestly-trusted count. `β` sets the **absolute** count trust so that "modest" is calibrated: low enough
that a good strand wins, high enough that the count cleanly governs where the strand is gone.

This is exactly what the old `w = I/(I+I₀)` did — Gaussian fusion with **count precision hard-capped at
`I₀=10`**. `I₀` was a **hard-coded magic number** we accepted "temporarily, to derive later." `β`
generalises it (`I₀` is the special case `β` chosen so the count contributes a fixed `I₀` of precision) and
gives us the path to *derive* it.

## 2. Why `β` must eventually be per-node (the measurement)

The grounding experiment (flagship, ss 0.99, capture-on; strand vs count error vs oracle) showed the
count's reliability is **not uniform**:

| node class | count MAE | reliable? |
|---|---|---|
| count-observable intron / intergenic | **0.005** | yes — trust it (`β` high) |
| imputed exon (across the capture jump) | **0.47** (bias −0.46) | no — distrust (`β` low) |
| AMBIG (count only) | 0.38 | no — distrust; defer to propagation |

A single `β` (like the old single `I₀`) under-trusts the excellent introns and over-trusts the biased
exons. The endgame is a `β_node` that is **high where the count is a direct measurement and low where it is
imputed across the enrichment discontinuity** — using signals we already compute (count-observability, the
`var~mean` curve, capture-class crossing). The count error is **bias**, so `var~mean` (variance) alone is
insufficient; observability/class-crossing carry the bias information.

## 3. Phased plan (simple → sophisticated)

Each phase is gated on the **flagship benchmark** (the old net leak: **~1%** capture-off, **8.7%** ss 0.99
capture-on, **20.3%** ss 0.50 capture-on) — a phase ships only if it does **not regress** and ideally
improves.

- **Phase 1 — single hard-coded `β` (this PR).** One global count trust, wired as a per-node array filled
  with one value. Reproduces the old `I₀`-cap robustness inside the simplex pie. **Acceptance:** matches the
  old flagship leak (no regression) and the key-insight check holds — at κ=0.99 the strand governs
  single-strand nodes, at κ=½ the count governs. `β` is the new acknowledged magic number (documented,
  config-exposed), the successor to `I₀`.
- **Phase 2 — 2-level `β` by count-observability.** `β_high` for count-observable nodes (direct
  measurement: introns, intergenic), `β_low` for imputed nodes (exons, AMBIG). Maps directly onto the
  measured introns-good / exons-bad split. **Acceptance:** capture-on leak drops vs phase 1 (the biased
  exon count is down-weighted).
- **Phase 3 — continuous `β(observability, var~mean, capture-class-crossing)`.** A smooth reliability
  function; `var~mean` modulates within a class, the capture-class-crossing flag carries the bias. No cliff.
- **Phase 4 — derived / calibrated `β`.** Fit the `β`-mapping on the benchmark so the *fused* estimate's
  MAE provably beats either signal alone (the §2 measurement becomes the calibration target). This also
  **derives the `I₀`/`β` scale** — fulfilling the original "derive later" promise; no magic number remains.
- **(Parallel track) Propagation of strand quality to count-only nodes.** Orthogonal to `β`: the RTS carries
  single-strand neighbours' *strand-derived* density (bias −0.15) into AMBIG nodes, where a low `β` lets it
  govern over the local biased count (−0.38). Built (simplex_propagate); the exclude-self/no-double-count
  rule is the remaining care. Layered after phase 2.

## 4. Phase-1 implementation

**Where `β` enters.** `solve_node`'s count term is `−½ · count_precision · (f_g − count_gdna_frac)²`. The
**effective count precision = `β`** (phase-1: a fixed value, the `I₀`-analog) — i.e. the count contributes a
fixed, capped amount of evidence, exactly like the old `I₀`, while the strand term keeps its intrinsic
`N·(2κ−1)²` curvature. So the per-node fusion is `≈ w·g_strand + (1−w)·g_count`, `w = I/(I+β)` — the old
behaviour, now inside the pie (with no-over-subtraction + the gDNA prior).

**The count value** fed to `solve_node` is the **splice-upgraded** `region_count_frac` (the count-mean-bias
clue), optionally RTS-smoothed along enrichment-class chains (the propagation architecture, kept on). The
**strand** is the node's own `u_pos/u_neg`. AMBIG nodes (no strand) ⇒ the count governs (β trust), as in
the old path.

**Wiring.** `calibrate.use_propagation` → `propagate_simplex`, passing `count_gdna_frac=region_count_frac`
and a config `count_trust_beta`. Output interface unchanged (`derive`/`priors`/`result` untouched).

**Config.** Add `CalibrationConfig.count_trust_beta: float` — documented as the phase-1 placeholder
(successor to `gdna_strand_info_scale=I₀`), to be replaced by the derived per-node `β` in phases 2–4.

**No-magic-numbers.** `β` is one new, *documented*, config-exposed constant that **replaces** the role of
the existing hard-coded `I₀` — not an addition on top. Its value is set by the key-insight constraint
(count ≪ strand at κ=0.99) and tuned on the flagship; the phases derive it.

## 5. Validation

1. **Unit:** at κ=½ the count governs `solve_node` (strand flat); at κ=0.99 a deep single-strand node
   follows the strand, not the count (the key-insight check) — add to `test_simplex`.
2. **Flagship A/B** (`--use-propagation` on vs off): pool gDNA fraction and net leak must **not regress**
   vs the old path; capture-off ~1% preserved.
3. **Full gDNA suite** before flipping the default.

## 5b. Phase-1 result (2026-06-13) — parity achieved

Implemented: `simplex_propagate.deconv_regions_simplex` — `solve_node` (own strand at `N·(2κ−1)²` + the
splice-upgraded count at fixed `β=count_trust_beta`) for **strand-observable** POS/NEG regions; **count-only
for AMBIG / intergenic** (no valid sense — exactly production's `w=0` rule). The gDNA log-prior is dropped
(`gdna_prior_count=0`): the count signal carries the gDNA default (its global-fallback baseline), and the
log-prior blew up at `f_g→0`, manufacturing false gDNA. Boundary sides keep the production `deconv_sides`.

Flagship A/B (pool gDNA fraction, true ≈ 0.75; production OFF vs simplex-β ON):

| condition | OFF | ON | Δ |
|---|---|---|---|
| gdna300 ss 0.99 capture-on | 0.6859 | 0.6793 | −0.0066 |
| gdna300 ss 0.50 capture-on | 0.5991 | 0.5925 | −0.0065 |
| gdna300 ss 0.99 capture-off | 0.7424 | 0.7398 | −0.0026 |
| gdna_none ss 0.99 capture-on (FP) | 0.0000 | 0.0038 | +0.0038 |

**Parity** (within ~0.7%) — the catastrophic gut regression (0.686→0.628) is gone, the elegant
`solve_node` + pie is restored, and the toy AMBIG-nascent test (gDNA=0) reads RNA (<0.2). 232 calibration
tests pass. The small consistent −0.005 (the β combine trusts the count marginally more than the old
`w=I/(I+I₀)`) and the tiny `gdna_none` FP (0.004) are within phase-1's no-regression bar; phases 2–4
(per-node β; strand→AMBIG propagation) are the improvement path, not regression fixes.

## 6. Open questions

- The exact phase-1 `β` scale (solve_node's count-term precision units vs the strand BB curvature) — tune
  on the flagship so the key-insight check holds; document the value.
- Whether to keep the RTS count-density smoothing on in phase 1 or defer it (it is a minor enhancement; the
  count trust `β` is the phase-1 subject).
- Phase 4 calibration target: minimise fused MAE, or directly minimise net leak on the suite?
