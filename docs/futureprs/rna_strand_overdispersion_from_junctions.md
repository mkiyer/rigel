# Future PR — Learn the (technical) strand overdispersion from per-junction counts → RNA fit + gDNA cross-channel prior

**Status:** Backlog / future work. Not urgent — PR 9's robustness floor (RNA) and the
MLE-plus-`0.01`-fallback (gDNA) are adequate until real-data evidence shows the strand
overdispersion matters. Recorded so the open modelling question isn't lost.

## What "strand overdispersion" is — and that it is purely TECHNICAL

A region's / junction's observed **sense fraction** deviates from its mean (`κ_rna`
for RNA, `κ_d = 0.5` for gDNA) for two reasons: binomial **sampling**, and
**overdispersion** — extra region-to-region spread beyond sampling. That
overdispersion is **entirely technical on both channels**:

- **RNA is single-stranded.** A transcript's reads are all one strand — there is **no
  biological junction-to-junction strand variation**. The residual *antisense* reads
  are **imperfect strand-specific library prep** (a technical artifact), and how that
  residual fraction varies across junctions is driven by **sequence accessibility,
  open/closed chromatin, GC content, secondary structure, and local mappability**.
- **gDNA is double-stranded** — reads come from both strands, so `κ_d = 0.5`; its
  deviation from 0.5 is driven by the **same technical factors** (which strand maps
  better in a region).

So RNA and gDNA strand overdispersion are **the same technical phenomenon**, exposed
to almost the same factors — not two different quantities. Call the shared parameter
**`ρ_od`** (the technical strand overdispersion). Because the Beta-Binomial
overdispersion is an **intra-class correlation** (κ-independent in
`Var = nκ(1−κ)[1+(n−1)ρ]`), `ρ_od` transfers between the channels' different means
(0.5 vs `κ_rna`) **without rescaling**: `ρ_r ≈ ρ_d`. (An earlier draft wrongly posited
an RNA "biological" excess giving `ρ_r ≥ ρ_d` — there is none; the transfer is direct.)

## The limitation today

Neither channel actually *learns* `ρ_od`:

- **RNA (PR 9):** `ρ_r_bb = 1/(n_obs+3)` is the **estimation uncertainty** of the
  pooled `κ_rna` (→ 0 with data) — not the region-to-region spread. Pooled 2×2 counts
  cannot measure overdispersion (the pooling destroys the per-junction scatter).
- **gDNA:** `ρ_d_bb` **is** learned — a 1-D MLE on the symmetric BB NLL over the
  per-region soft-allocated gDNA counts (`fit_rho_d_bb`, contained regions only →
  independent, no double-count). But with `< 2` gDNA regions it falls back to a magic
  constant `_RHO_D_BB_FALLBACK = 0.01`, with no principled basis.

## Goal

Learn the single shared `ρ_od` from the **well-powered** channel — the RNA junctions
(10⁵–10⁶ of them) — and use it both as the RNA `ρ_r_bb` and as the **principled,
data-driven prior/anchor for the sparse gDNA `ρ_d_bb`**, replacing the `0.01` magic
number.

## Design — a per-junction accumulator in the C++ scanner

Add a **junction table** to the accumulator, keyed by the junction's unique coordinate:

```
(ref_id, intron_start, intron_end, sj_strand)  →  (sense_count, antisense_count)
```

- We **already parse** these: each spliced read carries its `sjs` (start, end, motif
  strand) and is tested against the splice-artifact blacklist (`filter_blacklisted_sjs`).
- **sense vs antisense** = `align_strand == sj_strand` → sense, else antisense (after
  the R2 flip; the self-consistent R1-orientation frame from the PR 8 diagnosis).
- **Summing** the columns reproduces the pooled `(n_same, n_opp)` (a cross-check on the
  `StrandModel` 2×2); the **per-junction rows** are the new signal for the `ρ_od` fit.
- Tabulate **post-blacklist** (clean RNA junctions only).

## The double-count (fragments spanning multiple junctions)

A spliced fragment can cross >1 junction; counting it at each correlates those
junctions (the per-junction analog of PR 10's boundary double-count) and **biases
`ρ_od` (toward underestimation** — shared fragments make junctions look alike). The
fit needs **independent** per-junction observations. Options:

1. **Fractional `1/k`** — conserves the total but the junctions still *share* the
   fragment → not independent, and needs rounding / a weighted likelihood. Good for
   *mass* conservation, not the variance fit.
2. **Random single junction** — independent, but **non-deterministic** → breaks the
   accumulator's byte-for-byte reproducibility (`tests/native/_accumulator_reference.py`).
   Avoid raw randomness.
3. **Deterministic single junction** (canonical — e.g. leftmost `(start, end)`).
   Independent **and** deterministic. *Recommended if multi-junction fragments are
   worth keeping.*
4. **Single-junction fragments only** — tabulate only fragments crossing exactly one
   junction. Independent, deterministic, simplest; drops the multi-junction minority.
   *Recommended starting point* (a global `ρ_od` has ample single-junction support).

**Recommendation:** for the `ρ_od` fit prioritise independence + determinism — option 4
to start, or option 3 to keep all fragments. Fractional accumulation is the wrong tool
here (it's for mass conservation).

## The fit

Hierarchical Beta-Binomial over the per-junction rows `(k_j sense, n_j total)`:

```
κ_j   ~ Beta(κ_rna, ρ_od)            # technical region-to-region spread
k_j | κ_j ~ Binomial(n_j, κ_j)
⇒ k_j ~ BetaBinomial(n_j, κ_rna, ρ_od)   (marginal)
```

Fit `(κ_rna, ρ_od)` by MLE (MoM warm-start → 1-D/2-D optimise). With 10⁵–10⁶ junctions
this is well-determined — this *learns* `ρ_od`, vs PR 9 which assumes the spread is 0.

## gDNA cross-channel companion — the principled prior for the sparse `ρ_d_bb`

This is the second half, and the answer to "how do we set the gDNA fallback
principledly?" — **use the RNA-learned `ρ_od`.** Because the two channels' strand
overdispersion is the *same technical phenomenon* (above), the well-powered RNA
estimate is a direct, data-driven estimate of the gDNA overdispersion too:

- **Replace** `_RHO_D_BB_FALLBACK = 0.01`: when gDNA is too sparse to fit `ρ_d_bb`
  (`< 2` regions), use `ρ_od` — a data-driven, principled value, not a magic number,
  and one that's *moderate* by construction (neither the Binomial knife-edge nor the
  uniform collapse).
- **When gDNA is well-sampled**, keep the existing MLE but as a **MAP** with `ρ_od` as
  a weakly-informative prior (strength ≈ unit information), so the gDNA channel can
  still depart from `ρ_od` if its own data demands — but is anchored when sparse.
- **Why this is the right source (esp. for hybrid capture):** intergenic gDNA is
  *depleted* under capture, so any gDNA-only prior would itself be sparse; the RNA
  junctions live in the *captured* regions and stay well-powered. And there is no
  implemented pure-gDNA per-region strand sample today (the `intergenic` StrandModel in
  the docstring was never built), so RNA is the only well-powered strand-spread source.

## Integration

```
ρ_r_bb = ρ_od                                  # RNA channel uses the learned spread
ρ_d_bb = MAP(gDNA per-region MLE | prior = ρ_od)   # gDNA: learn when it can, anchor to ρ_od when sparse
```

## Constraints / notes

- **Determinism:** the accumulator must match the Python reference byte-for-byte
  (CLAUDE.md) → deterministic multi-junction assignment, no raw RNG; update
  `tests/native/_accumulator_reference.py` in lock-step.
- **Memory:** a hash map over ~10⁵–10⁶ junctions × 2 counts is small.
- Independent of the remaining `nrna_dc g20` count/EM failures and the edge-case
  test-pool PR.

## Open questions

1. Double-count strategy — single-junction-only (simplest) vs deterministic-single?
   Decide on real data.
2. Fit method — MLE vs MoM vs a light Bayesian posterior on `ρ_od` (the latter gives
   *its* uncertainty, to set the gDNA prior strength cleanly).
3. gDNA prior strength — unit information, or scaled by how confidently `ρ_od` is
   measured?
4. Validation — a real stranded RNA-seq + gDNA library where `ρ_od` is present and
   measurable; confirm fewer deep-region mis-calls on both channels.

## Priority / effort

Low priority (PR 9's RNA floor + the gDNA MLE/fallback suffice until real data shows
otherwise); moderate effort (C++ junction accumulator + Python reference + the fit +
wiring into both `fit_strand_balance` and `fit_rho_d_bb`). A clean, self-contained
enhancement that closes the overdispersion story on both channels at once.
