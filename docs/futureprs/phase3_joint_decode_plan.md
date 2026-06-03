# Phase 3 — Joint Decode: Theoretical & Implementation Plan (v2)

Redesigned after review. **Goal:** independently deconvolve every **node** (each region and
each boundary, kept separate) into gDNA mass and RNA mass, combining the count clue (Phase 1)
and the strand clue (Phase 2) as **two orthogonal estimates**, single-pass and acyclic, with
one continuous confidence knob, graceful at zero-gDNA and zero-RNA, and numerically stable at
all inputs.

## 0. Decisions carried in from review

1. **Orthogonal clues; raw density.** We do **not** strand-correct the density or use
   leave-one-out (deferred enhancement). The count density is the raw Phase-1 estimate and is
   **honestly upward-biased** by unspliced nascent-RNA fragments in introns / exon-intron
   boundaries. On toy scenarios with heavy nascent RNA we will *over-call gDNA*; on real data
   nascent RNA is real but sparse, so the global density bias is small, and the **downstream
   per-locus EM** (with the strand model) deconvolves most nascent RNA anyway.
2. **Single pass first.** The iterative multi-pass (the unstranded expressed/unexpressed
   latent, §10) is a planned enhancement, validated only after single-pass works.
3. **Fractional mass for the amount; discrete fragment count for the statistical power.** The
   accumulator's fractional mass is conserved (no double-count) and gives the *amount*; the
   discrete supporting fragment count gives the *uncertainty/precision*. The two are never
   conflated.
4. **Strand likelihood is a Beta-Binomial** with a **data-driven** overdispersion and a
   data-driven sparse floor (no magic numbers); the very-sparse limit should land near a
   Beta(3,3)-shaped split (the conservative high-overdispersion regime).
5. **Nodes are independent and regions/boundaries stay separate.** Each node is deconvolved on
   its own; masses are never merged across nodes; the deconvolved per-node masses are combined
   into loci only **after** calibration. (The density *rate* is still shared via the Phase-1
   sweep — that borrows a rate, never mass.)
6. **Graceful at both extremes:** zero (or near-zero) gDNA and zero (or near-zero) RNA must be
   handled with no cliffs and no magic numbers.
7. **Efficient and provably numerically stable** at all inputs.

---

## 1. The node and the estimand

A **node** is a single region **or** a single boundary. For each node, independently:

- `M` — fractional mass (conserved by the accumulator; the *amount*).
- `(S, A)` — discrete oriented sense / antisense flux, `N = S + A` (the *power*).
- `eff_len` — `region_eff_len` for a region, `μ_FL` for a boundary.
- signature → count-decodability (Phase 1) and strand-decodability (POS/NEG vs AMBIG/NONE).

Latent: the gDNA fraction `π_g ∈ [0, 1]`. Output: `gdna_mass = π_g · M`, `rna_mass = (1−π_g)·M`.
**Regions and boundaries are decoded separately and combined only after calibration.**

## 2. The two orthogonal clues, per node

- **Count clue (Phase 1, raw):** the swept density `d` predicts an expected gDNA mass
  `μ_G = d · eff_len`, hence a count-predicted fraction
  `π_g^count = clip(μ_G / M, 0, 1)`. Its **strength** is the discrete gDNA count evidence
  behind `d` (the sweep's accumulated supporting count), `κ_c`. *Raw → upward-biased by
  nascent RNA (§0.1), accepted.*
- **Strand clue (Phase 2, Beta-Binomial):** the discrete `(S, A)` split gives a likelihood
  over `π_g`. Its strength is the discrete `N`, capped by the data-driven overdispersion (§6).

The two are orthogonal observables — *count = how much mass / how many fragments*, *strand =
their orientation* — so they multiply without double-use, even on the same node.

## 3. The joint, per node

Posterior over the gDNA fraction:

```
posterior(π_g) ∝ P_count(π_g)  ·  L_strand(π_g)
```

- `P_count(π_g) = Beta(a_c, b_c)` with mean `π_g^count` and concentration `κ_c` (the discrete
  count evidence): `a_c = κ_c·π_g^count + ½`, `b_c = κ_c·(1−π_g^count) + ½`. The `½` is the
  **Jeffreys** floor (principled, not magic): at zero count evidence the prior is the diffuse
  Beta(½,½) and the strand decides.
- `L_strand(π_g)` is the §6 Beta-Binomial mixture likelihood of `(S, A)` given `π_g`.

The reported gDNA fraction is the posterior **`z`-quantile** (§9; `z=0` = mean). Then
`gdna_mass = π_g · M` (fractional amount), `rna_mass = (1−π_g)·M`. Precision-weighting between
the clues falls out of the product, with the correct limits (§4).

## 4. Graceful limits (no cliffs)

| input | behaviour |
|---|---|
| **zero gDNA** (`d → 0`) | `π_g^count → 0` → `a_c → ½` → posterior → strand; if strand also says RNA, `π_g → 0` |
| **zero RNA** (pure gDNA) | strand unstranded + `μ_G ≥ M` → `π_g^count = 1` → `π_g → 1` |
| **silent thin node** | `μ_G > M` ⇒ `π_g^count = 1`; strand diffuse ⇒ posterior → 1 (all gDNA) |
| **expressed exon** | sharp strand ⇒ posterior → mostly RNA |
| **AMBIG** | no strand likelihood ⇒ posterior = `P_count` ⇒ `π_g = π_g^count` |
| **N = 0** (no flux) | strand flat ⇒ posterior = `P_count` |
| **M = 0** (no mass) | skip — `gdna = rna = 0` |

Both extremes (§0.6) degrade continuously: zero gDNA → `π_g→0`, zero RNA → `π_g→1`, with no
thresholds.

## 5. Statistical power — fractional mass vs discrete count

The **fraction** `π_g` is inferred from the **discrete** signals — the count evidence `κ_c`
(prior concentration) and the strand `(S, A)` (likelihood) — so a node with large fractional
mass but few supporting fragments is correctly *uncertain*. The **mass** is then `π_g · M`
using the conserved fractional `M`. This is the §0.3 separation: discrete counts set every
precision; fractional mass sets only the final amount.

## 6. The strand likelihood (Beta-Binomial, data-driven)

Of the `N` discrete fragments, a fraction `π_g` are gDNA (oriented-sense rate ½) and `1−π_g`
are RNA (rate `κ_rna`). The per-fragment sense rate is **Beta-distributed** (overdispersion),
so the count is Beta-Binomial, not Binomial. We use the moment form: mean and variance of `S`
given `π_g`,

```
μ_S(π_g) = N·[½·π_g + κ_rna·(1−π_g)]
σ²_S(π_g) = N·p(1−p)·[1 + ρ·(N−1)]          p = ½·π_g + κ_rna·(1−π_g)
L_strand(π_g) ≈ Normal(S | μ_S, σ²_S)        (exact BB-mixture for very small N — §8)
```

`ρ` is the **strand overdispersion**, estimated from the data (the across-node excess
variance of the sense fraction beyond Binomial — the existing `ρ_d_bb` machinery), with a
**data-driven floor**: as the supporting evidence vanishes, `ρ` rises to its conservative
limit so a node with `N=1–2` cannot be confidently split. Target sanity check: the very-sparse
posterior over `π_g` should look ≈ **Beta(3,3)** (concentration ~6, the reviewer's intuition)
— we validate against that rather than hard-coding it. `ρ = 0` recovers the Binomial; the
ceiling is the Binomial limit (well-sampled data).

## 7. Numerical stability (at all inputs)

Claim: the posterior over `π_g ∈ [ε, 1−ε]` is proper and bounded for every admissible input.

- `P_count` is `Beta(a_c, b_c)` with `a_c, b_c ≥ ½ > 0` (Jeffreys floor) ⇒ always a proper,
  bounded density; `κ_c → 0` ⇒ Beta(½,½) (diffuse, finite).
- `L_strand` is a Normal in `S` with `σ²_S ≥ N·p(1−p) > 0` for `π_g ∈ (0,1)` and `κ_rna ∈ (0,1)`
  away from the singular endpoints; the `ε` grid guard keeps `p ∈ (0,1)`. For `N=0`, `L_strand`
  is constant ⇒ posterior = `P_count`.
- `μ_G/M` is computed only when `M>0` (else the node is skipped) and clipped to `[0,1]`, so
  `π_g^count` is always defined; `d → ∞` ⇒ `π_g^count = 1` (finite).
- No logs of zero, no divisions by zero, no unbounded terms. The product of a bounded prior
  and a bounded likelihood over a compact interval is integrable ⇒ the normalized posterior,
  its mean, and any quantile exist and are finite. ∎ (sketch; to be unit-tested at the
  extremes: `M=0`, `N=0`, `d=0`, `d=∞`, `κ_rna→0.5`).

## 8. Computational efficiency

Both factors are ~Gaussian in `π_g` for moderate `N`: `P_count` is Beta (≈Gaussian when
concentrated), `L_strand` is Normal in `S` (≈Gaussian in `π_g`). So the **default path is a
closed-form Gaussian product** (precision-weighted mean + variance) — O(R), no grid. A small
`π_g`-grid (or the exact BB-mixture) is used only for **low-`N` / near-boundary** nodes where
the Gaussian approximation is poor. This keeps 10⁶-node genomes O(R) with a bounded grid
fraction.

## 9. Confidence knob `z`

`z` selects a quantile of the per-node posterior of `π_g`: `z=0` = mean (unbiased, default);
`z>0` = upper quantile → more gDNA → purer RNA (specificity); `z<0` → RNA sensitivity.
Continuous, applied uniformly (strand-decodable and AMBIG count-only). It propagates to the
reported gDNA mass `= π_g(z)·M`.

## 10. Single-pass flow, and the deferred iterative enhancement

**Single pass (this version):** Phase-1 raw swept density (per node, regions + boundaries) +
Phase-2 BB strand → §3 per-node joint → per-node `(gdna, rna)`.

**Deferred — iterative classification for unstranded data.** When strand is unavailable, one
pass gives the seed density from the count-decodable "seed" nodes (intergenic/intronic regions,
exon-intron / exon-intergenic boundaries). Then classify the unknown nodes (exons, AMBIG) as
**expressed vs unexpressed** against that density: a node with **spliced** fragments is
*definitively expressed*; an **unspliced-only** node may be *unexpressed* and, if its density is
consistent with the seed, folds into the gDNA pool → re-estimate density → iterate, holding the
fixed seed pool and fixed expressed pool, only the unknown pool migrating. Cyclic but anchored
by the fixed pools; deferred until single-pass is validated.

## 11. Implementation

**New — `src/rigel/calibration/joint_decode.py`:**
- `JointDecode` dataclass: per-node `gdna_mass`, `rna_mass`, `pi_g`, `pi_g_var`, `is_boundary`.
- `joint_decode(nodes, *, kappa_rna, strand_overdispersion, confidence=0.0)` — the §3 posterior
  per node (closed-form Gaussian product default; grid for low-`N`), returning the `z`-quantile.

**Modifications:**
- `density_model.py` (Phase 1): emit **per-node** (region **and** boundary) raw density + the
  **discrete count evidence** `κ_c` (already in the sweep `(α, β)` — expose, and carry a
  discrete-count companion to `α`).
- `strand_decode.py` (Phase 2): upgrade to the Beta-Binomial likelihood (§6) with the
  data-driven `ρ` + floor; expose `L_strand(π_g)` as a function reusable by the joint.
- A small **node abstraction** (regions ∪ boundaries) so the joint treats them uniformly while
  keeping their masses separate.

**Torn down:** nothing yet (additive; Phase 6 removes the old `estep`/`mstep`/`exposure`).

## 12. Validation

1. **nrna_dc toy:** RNA recovery high; gDNA *somewhat over* (honest nascent bias, §0.1) but
   bounded; `ρ₀` bounded.
2. **Overlapping locus:** AMBIG interior count-only, edges via strand, consistent.
3. **Zero-gDNA scenario:** `π_g → 0` everywhere, RNA ≈ total, no NaNs.
4. **Zero-RNA scenario:** `π_g → 1`, gDNA ≈ total, no NaNs.
5. **Stability unit tests:** `M=0`, `N=0`, `d=0`, `d=∞`, `κ_rna→0.5` — all finite (§7).
6. **3-regime harness** + **`z` sweep** (monotone purity/recall).
7. **Sparse strand sanity:** very-low-`N` `π_g` posterior ≈ Beta(3,3) (§6).

## 13. Remaining open items

- The **data-driven `ρ` floor** (§6): the exact estimator + that its sparse limit matches the
  Beta(3,3) intuition — to be fit and validated, not hard-coded.
- **Boundary-node deconvolution is new**; confirm per-node mass conservation (Σ decoded gDNA +
  RNA over all nodes = total fragment mass) and that keeping boundaries separate doesn't
  drop/duplicate mass.
- The **raw-density nascent bias** is accepted for v1; quantify it on toys so we know the size
  of the effect the downstream EM must absorb, and revisit the strand-corrected density as the
  documented enhancement if it proves too large.
- **Single-pass self-consistency** for stranded data — verify empirically before considering
  the §10 iterative path.
