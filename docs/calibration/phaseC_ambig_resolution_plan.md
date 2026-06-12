# Phase C (merged B+C) — the per-node likelihood deconvolution: implementation plan

**Status:** dry-run plan, rev 2 (2026-06-12). Third phase of
[`deconvolution_roadmap.md`](deconvolution_roadmap.md), merging former Phase B. **Rev 2 replaces the
heuristic blend + hard mature-subtraction with a proper per-node likelihood combination over the 3-simplex
`{RNA+, RNA−, gDNA}`** — every piece of evidence is a likelihood that pushes the partition; no estimate is
ever "fully trusted." Builds on Phase A (`mature_density.py`, `83b3610`) and the Step-0 unit-consistency
record ([`phaseB_mature_subtraction_plan.md`](phaseB_mature_subtraction_plan.md) §5.5).

## 1. The model — one pie, three slices, evidence as likelihoods

Each node's total unspliced mass `U` is partitioned into `(f₊, f₋, f_g)` on the 2-simplex
(`f₊+f₋+f_g = 1`): the `+`-strand RNA, the `−`-strand RNA (each = nascent + contained-mature), and gDNA.
The **region's `strand_class` sets the simplex constraints**, so one model spans every region type:

| class | constraint | free params |
|---|---|---|
| `NONE` (intergenic) | `f₊=f₋=0` | none → `f_g=1` (all gDNA) |
| `POS` | `f₋=0` | 1 |
| `NEG` | `f₊=0` | 1 |
| `AMBIG` | — | 2 |

The posterior is the **product of independent likelihoods** (× a weak prior); the answer is its MAP (or
mean). The calibration reads off **`f_g`** (gDNA fraction); `f₊/f₋` are secondary (they let the strand and
mature constrain `f_g` and feed the EM downstream).

## 2. The key reduction — what each evidence actually constrains

Reparameterize by **RNA fraction `r = f₊+f₋ = 1−f_g`** and **tilt `t = f₊−f₋`** (so `f₊=(r+t)/2`,
`f₋=(r−t)/2`, `|t| ≤ r`). Then the genome-`+` read fraction is

```
p_pos = ½·f_g + κ·f₊ + (1−κ)·f₋  =  ½ + (κ−½)·t
```

**The strand sees ONLY the tilt `t` — not `f_g`.** gDNA and *balanced* RNA both contribute ½ to pos, so
they are indistinguishable by direction (the flat direction, made exact). This is the whole identifiability
story in one line, and it tells us how to combine the evidence:

- **Strand likelihood** `L_strand(t)` = `BetaBinom(U_pos | U, ½+(κ−½)t, od)` — pins **t**. Sharp when κ is
  far from ½; **flat when κ≈½** (then `p_pos≡½`, t unconstrained by strand — correct).
- **Mature likelihoods** `L_mat₊(f₊)`, `L_mat₋(f₋)` — the bounding **spliced** crossings (motif-oriented,
  **strand-independent**) say each strand's RNA is **at least** its mature: a soft **lower bound** on `f₊`,
  `f₋` (hence on `r` and `t`), width from the flux Poisson **+ the geometry scatter** (§5).
- **gDNA likelihood** `L_gdna(f_g)` — the bounding **unspliced** crossings give a capture-**depleted**
  estimate: a soft **lower bound** on `f_g`, width from the crossing Poisson (+ a depletion allowance).

So `t` is pinned by strand (and/or the mature tilt); `f_g` is squeezed between the **gDNA lower bound**
(pushes `f_g` up) and the **mature-implied RNA** (pushes `f_g` down). The MAP balances them by precision —
**the silent-strand override, the balanced-nascent floor, and the single-strand result all emerge from
this, with no `ρ`, no gate, no subtraction.**

## 3. How every regime falls out (the unification, verified)

- **Single-strand `POS`** (`f₋=0 ⇒ t=r`): `L_strand` pins `r` directly (`p_pos = ½+(κ−½)r`) ⇒ `f_g=1−r` —
  **exactly today's strand module**. The mature lower bound only *refines* (corrects a strand RNA
  under-estimate); with a realistic width (§5) it cannot override a sharp, reliable strand. Subsumed.
- **AMBIG, one strand silent** (`M₋≈0`): the `−` mature lower bound is ~0, so nothing opposes pushing the
  symmetric mass into gDNA; `f_g → 1−|t|` (max). Strand pins `t`. Recovers gDNA. *(flagship 226.)*
- **AMBIG, both expressed** (`M₊,M₋` large): both mature lower bounds push `r` up ⇒ `f_g` down toward the
  gDNA lower bound — symmetric mass read as balanced RNA, not gDNA. *(passes `test_ambig`: gDNA=0+nascent ⇒
  gDNA lower bound ≈0 ⇒ `f_g≈0`.)*
- **Unstranded (κ≈½)**: `L_strand` flat. `t` is then pinned by the **mature tilt** (`M₊−M₋`, motif-based,
  still available); `f_g` squeezed between the gDNA lower bound and the mature-implied RNA. The
  nascent-vs-gDNA part that neither mature nor gDNA-crossing constrains stays **genuinely uncertain** — the
  posterior is *wide* there rather than inventing a value. **Optimal without strand**: every identifiable
  quantity (mature by motif, gDNA lower bound) is used; only the truly-confounded part is left open.
- **Intergenic / pure gDNA**: `r=0 ⇒ f_g=1`. Trivial.

## 4. Implementation

A per-node MAP over `(r, t)` (1 free param for single-strand, 2 for AMBIG) — **independent per node, one
feed-forward pass, no genome loop.** Reuses the grid machinery the strand module already has.

- **`strand_likelihood` / the strand posterior**: generalize from the 1-D gDNA-fraction grid to the
  `(r, t)` grid (the existing `strand_posterior_gdna_frac` is the `f₋=0` slice). `L_strand` depends only on
  `t`; the BB overdispersion is the composition-weighted mix of `od_gdna` / `od_rna` (detail §6).
- **count + mature → likelihoods, not points**: `density_model` already carries a posterior **variance**
  (`count_gdna_frac_var`); expose `g_count` + its width as `L_gdna`. `mature_density` (Phase A) gains a
  **width** per strand (flux Poisson + geometry scatter) → `L_mat₊/₋`. Both are **soft one-sided** (lower
  bounds): `RNA_s ≥ mature_s`, `gDNA ≥ depleted estimate`.
- **combine + solve**: a new `node_posterior` that builds the log-posterior on the `(r,t)` grid (respecting
  the `strand_class` constraint), adds a weak Dirichlet prior, and returns the MAP `f_g` (+ optionally the
  mean and a width for the FP-quantile). `deconv_regions` calls it; **`deconv_sides` unchanged** (sides are
  mature-free; their 2-component strand posterior is the `f₋=0`/`f₊=0` slice — already handled, or routed
  through the same code with the side's class).
- **retire**: the heuristic `w = I/(I+I₀)` blend in `_deconv_per_node` **and** `region_splice_gdna_frac` —
  both are subsumed by the likelihood combination. Keep `splice_junction_eligibility` (mature anchors) and
  the BB posterior core (`strand_loglik`).

Output unchanged downstream: `f_g·U` = gDNA mass, `(1−f_g)·U + spliced` = RNA mass; `derive`/`priors` as-is.

## 5. The one real derivation — the mature (and gDNA) likelihood WIDTH

A lower-bound likelihood that is **too sharp** lets a biased `M` override a reliable strand (e.g. region
222: `M₋` over-predicts ~16% in substrate units; a sharp bound would wrongly pull `f_g` 0.55→0.49). **Too
wide** and it can't resolve the AMBIG flat direction. So the width must reflect `M`'s *true* uncertainty:

```
Var(M_s) = (flux Poisson: 1/S_s)  +  (geometry scatter: the empirical M/substrate-mature spread)
```

Step 0 already measured the geometry scatter (median ratio 0.93, p10–p90 spread) — so the width is
**derived from data, not a tuned constant**. With it, a sharp strand (single-strand) dominates a wide
mature bound (no harm), while at an AMBIG node — where the strand cannot constrain `f_g` at all — even a
wide mature bound is the governing evidence. This is the Bayesian self-balancing that replaces `ρ`. The
exact functional form of the one-sided term (clipped-Gaussian on `f_s·U − M_s`, or a profiled
Poisson-with-nascent-nuisance) is to be finalized + calibrated on the toy sweep — **no new constant**.

The gDNA likelihood is the analogous one-sided lower bound at `g_count` with its `count_gdna_frac_var`
width, plus a one-sided **depletion allowance** (true gDNA ≥ the depleted boundary estimate) — off-capture
(intron/intergenic) it tightens to two-sided.

## 6. Open modeling choices (settle in execution; none introduce a magic number)

1. **Mature/gDNA likelihood width + one-sided form** (§5) — the primary derivation; calibrate the width
   against the Step-0 scatter + the toy sweep.
2. **3-component BB overdispersion** — the strand sees only `t`; the overdispersion around `½+(κ−½)t` is a
   composition-weighted mix of `od_gdna`/`od_rna`. Confirm the mix (or that a single od suffices).
3. **MAP vs posterior-mean**, and a cheap solve: profile `t` (≈ closed-form from the BB) then 1-D in `f_g`,
   vs a small 2-D grid. (Perf is deferred; correctness first — but the profiled 1-D is likely both.)
4. **Boundary-imputation uncertainty** — the mature/gDNA come from *imputed* boundary crossings; their
   widths should include the imputation variance (the run-fill distance), not just the local flux Poisson.

## 7. Scope (locked) & what is deferred

- **One feed-forward pass.** The per-node MAP is a *local* solve, not iteration. The **genome-level**
  strand⇄count⇄boundary propagation (re-deconvolving boundaries from refined nodes and re-imputing) is the
  deferred "iteration" — add only after one-pass validates. If added, it is a fixed-point over the same
  likelihoods (the curvature is the precision — no new knobs).
- **Protocol-agnostic** (κ∈[0,1]); **no magic numbers** (the only shape choice is the §5 width, derived
  from data — pause & discuss if a bare constant threatens to appear); **capture-agnostic** (local
  evidence only).

## 8. Validation

1. **`test_ambig_scenario`** — both cases pass (§3).
2. **Toy AMBIG sweep** — gDNA 0→high × nascent {none, both, one-strand} × κ {0.99, 0.5}: silent-strand
   recovers gDNA; both-expressed → floor; unstranded → mature+gDNA carry it, nascent-vs-gDNA stays wide;
   **calibrate the §5 width here.**
3. **Single-strand / simple / no-gDNA / unstranded conditions** — **must not regress** (the model reduces to
   today's strand result there); golden refreshed and the diff reviewed for byte-stability where expected.
4. **Flagship locus 21** — AMBIG 226/231/236 recover to ~0.54–0.57; net leak drops.
5. **Full gDNA suite** — net `gdna→rna` leak falls on capture/stranded/complex; no regression elsewhere.
6. New debug script `scripts/debug/phaseC_node_posterior.py` (committed) — per-node `(r,t)` posterior,
   each likelihood term, MAP `f_g`, joined to oracle.

## 9. Risks

- **Scope** — this replaces the deconvolution *core* (strand+count+blend) with the unified posterior. Big,
  but principled and subsuming; the single-strand reduction (§3) is the guard against regression — verify it
  byte-reduces to the current strand module when the mature/gDNA likelihoods are flat.
- **Mature width** (§5) — too sharp harms single-strand, too wide forfeits AMBIG; the derived width is the
  fix, calibrated on real scatter.
- **2-D MAP cost** — per-node, deferred to the perf phase; profile-`t`→1-D keeps it ~the current cost.
- **Golden/suite drift** — review; confirm simple conditions are stable.
