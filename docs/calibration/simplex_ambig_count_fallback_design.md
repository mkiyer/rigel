# Design + implementation plan: gate the simplex-sweep Jeffreys prior at strand-uninformative nodes

**Status:** design + implementation plan, ready for review → implement. **Date:** 2026-06-14.
**Scope:** fixes the #1 failure family of the 2-D simplex sweep (`use_propagation=True`) exposed by the
complex-loci benchmark — the **no-anchor / anchor-starved AMBIG over-call**. The multi-transcript-per-strand
family (#2) is a *separate* follow-up doc. Prereq context: the sweep's per-node solver was missing the
`Beta(½,½)` Jeffreys prior (`simplex_sweep_missing_jeffreys_prior` memory); adding it (calibrated,
`_PRIOR_EPS=1e-3`, K≈60) made the sweep ≥ production on single-strand loci and 2× better on simple AMBIG —
but introduced this regression at AMBIG nodes that have no resolving evidence.

---

## 1. TL;DR

The Jeffreys prior the sweep adds to its per-node `ψ_i` is **unconditional**. At a single-strand node
that is the correct, essential tie-breaker (the strand likelihood is informative, the prior nudges the
balanced-gDNA node to the `f_g=1` vertex). But at an **AMBIG node whose `f_g` is under-determined**
(strand uninformative, no spliced anchor, no propagated-odds inflow) the *flat* likelihood lets the
Jeffreys **vertex-push** the node toward `f_g=1` — over-calling gDNA — and it **overrides the count
clue** (`g_count`) that the 1-D production solver correctly defers to there. **Fix:** scale the Jeffreys
contribution by the per-node strand-information weight `w = I/(I+I₀)`, `I = N·(2κ−1)²`, **with `I=0` at
AMBIG/strand-unobservable nodes** — exactly the gate the 1-D uses (`w=0 ⇒ count governs`). The Jeffreys
then only breaks the tie where there *is* strand evidence to break; where there isn't, the count term (and
the propagation's odds inflow, when present) governs.

---

## 2. The problem (evidence)

Complex-loci benchmark (`scripts/debug/complex_loci_benchmark.py`, 49 loci, gDNA=120/nRNA=25/ss=0.99,
sweep K=60). Mass-weighted AMBIG |f_g − oracle|, 1-D (production fusion) vs 2-D (sweep):

| locus | structure | 1-D err | 2-D err |
|---|---|---|---|
| `full_single` | + and − identical single exons (no anchor, no SJ) | 52 | **968** |
| `anchor_sj_neither` | AMBIG with SJ on neither strand | 116 | **589** |
| `nest_anchor_starved_core` | wide AMBIG core, no interior single-strand anchor, no interior SJ | 99 | **502** |

These three alone contribute ~+1800 of the sweep's total regression (battery total 2D/1D = 1.33). The
common structure: an AMBIG region with **no path to resolution** — strand uninformative (both strands),
no spliced reads to lower-bound an RNA strand, and no strand-informative neighbour for the propagation to
inherit odds from. The 1-D solver handles these correctly (it defers to `g_count`); the 2-D over-calls.

(Note: `anchor_intron_moat` 696/751 and `anchor_exon_vs_intron_flank` 367/391 are hard for *both* solvers —
near-unidentifiable — and are out of scope here.)

---

## 3. Root cause (the algebra)

`simplex_sweep._local_loglik` builds `ψ_i` over the 2-simplex lattice as

```
ψ_i(θ) = mixture_strand_loglik(θ)            # 3-component gDNA/RNA₊/RNA₋ Gaussian strand likelihood
       − ½·β·(f_g − g_count)²                 # count clue (β = count_trust_beta = 10)
       + (½−1)·[log f_g + log f₊ + log f₋]    # Dirichlet(½,½,½) Jeffreys prior  ← UNCONDITIONAL
       − spliced lower bound                  # 0 when no spliced reads
```

At an unanchored AMBIG node the strand likelihood is **flat along a ridge** (many `(f_g,f₊,f₋)` give the
same observed sense fraction) and the spliced bound is 0, so only the count term and the Jeffreys shape
the posterior. Evaluate the two at the `f_g=1` vertex vs the count's preferred `f_g≈g_count` (take
`g_count≈0.3`, the oracle):

- **count term:** `−½·β·(f_g−0.3)²` penalises `f_g=1` by `½·10·0.49 = 4.9` relative to `f_g=0.3`.
- **Jeffreys (3-component):** at the `f_g=1` vertex `f₊=f₋→0`, so the term `≈ (−½)(log 1 + 2·log _PRIOR_EPS)
  = (−½)(2·log 1e-3) = +6.9`; at `f_g=0.3, f₊=f₋=0.35` it is `(−½)(log .3 + 2 log .35) ≈ +1.65`. The
  Jeffreys thus **rewards `f_g=1` by +5.25** over `f_g=0.3`.

`+5.25 (Jeffreys) > 4.9 (count)` ⇒ the posterior is pushed to `f_g=1`. The **3-component** Jeffreys is
*stronger* than the single-strand one because it rewards *both* `f₊→0` and `f₋→0` at the gDNA vertex (two
log terms, not one) — precisely at the nodes where the strand cannot justify it.

Why the 1-D solver is right here: `_deconv_per_node` uses `center = w·g_strand + (1−w)·g_count`,
`w = I/(I+I₀)`, and AMBIG nodes are **strand-unobservable ⇒ w=0 ⇒ `g_count` governs**. The 2-D sweep has
no equivalent gate; its Jeffreys fires regardless of whether the strand can resolve the node.

---

## 4. The fix

**Scale the Jeffreys prior by the per-node strand-information weight `w`** (the same quantity the 1-D
uses), so it is active only where the strand is informative:

```
w_i = I_i / (I_i + I₀),   I_i = N_i·(2κ−1)²   for strand-OBSERVABLE nodes (TS_POS / TS_NEG),
      0                                        for AMBIG / intergenic (strand-unobservable),
ψ_i += w_i · (½−1)·[log f_g + log f₊ + log f₋]     # Jeffreys, gated
```

- `I₀ = config.gdna_strand_info_scale` (= 10, the same half-trust scale the 1-D uses).
- **Single-strand nodes** (`w→1` at ss=0.99): full Jeffreys — the balanced-gDNA tie-breaker is preserved,
  so the validated single-strand / simple-AMBIG wins are unchanged.
- **Unanchored AMBIG** (`w=0`): no Jeffreys ⇒ the count term `−½β(f_g−g_count)²` governs ⇒ matches the
  1-D's `g_count` fallback ⇒ the `full_single` / `anchor_sj_neither` / `nest_anchor_starved_core`
  over-call is removed.
- **Propagation-anchored AMBIG** (the 2-D's wins — promoter / interleave / telescope families): unchanged,
  because those nodes are resolved by the **odds inflow over the edges**, not by the local Jeffreys; the
  edge messages still apply. The fix only removes the *local* vertex-push that had nothing to justify it.

### Why this is principled, not a patch
The Jeffreys prior's *job* is to break the binomial-proportion identifiability tie **for the strand
posterior**. Where there is no strand evidence (AMBIG), there is no such tie to break — the resolving
evidence is the count (locally) or the propagated odds (spatially), exactly as in the 1-D. Gating the
prior by the strand information makes the sweep degrade to the count fallback at unanchored nodes
*by construction*, which is the correct behaviour the 1-D already embodies. No new constant (`I₀` is
reused).

---

## 5. Implementation plan

All changes in `src/rigel/calibration/simplex_sweep.py`; no C++; Python-only (editable install).

1. **Compute `w` per node** in `deconv_regions_sweep` (it already has `ts = region_arrays.strand_class`,
   `U = u_pos + u_neg`, `kappa`):
   ```python
   info = U * (2.0 * kappa - 1.0) ** 2
   strand_obs = (ts == TS_POS) | (ts == TS_NEG)         # AMBIG / NONE → unobservable
   w = np.where(strand_obs, info / (info + info_scale), 0.0)   # info_scale = gdna_strand_info_scale
   ```
   Add `info_scale` as a `deconv_regions_sweep` parameter (default `10.0`), threaded from
   `config.gdna_strand_info_scale` in `calibrate.py` (the sweep call already passes the other params).

2. **Pass `w` into `_local_loglik`** and gate the Jeffreys:
   ```python
   psi = psi + (w[:, None]) * (_STRAND_PRIOR - 1.0) * (
       np.log(np.clip(f_g, _PRIOR_EPS, 1.0))[None, :]
       + np.log(np.clip(f_pos, _PRIOR_EPS, 1.0))[None, :]
       + np.log(np.clip(f_neg, _PRIOR_EPS, 1.0))[None, :]
   )
   ```
   (The existing optional `gdna_prior_count·log f_g` stays as-is, default 0.)

3. **`calibrate.py`**: pass `info_scale=config.gdna_strand_info_scale` to the `deconv_regions_sweep` call.

4. **No change** to the strand mixture, the count term, the spliced bound, the edge/propagation, or the
   median extraction.

Estimated change: ~10 lines.

---

## 5b. Prototype validation (done — 2026-06-14)

The fix (§5) was prototyped (`simplex_sweep._local_loglik` `jeffreys_w` gate + the `w` computation in
`deconv_regions_sweep`, `info_scale=10`) and run on the 49-locus battery (K=60). It is a clear net win:

- **Battery total 2D/1D: 1.33 → 0.92** (the 2-D goes from worse-overall to better-overall).
- **Family-1 recovered substantially:** `nest_anchor_starved_core` 502→187, `full_single` 968→652,
  `anchor_sj_neither` 589→413.
- **2-D wins preserved AND improved** (no regression): `interleave_picket` 298→199,
  `nest_telescope_4deep` 462→339, all 6 promoter loci still win.

**Residual the prototype exposed (refines Open Question #2):** `full_single` is only *partially* fixed
(652 vs the 1-D's 52). With the Jeffreys gated off, the **strand-mixture likelihood term itself still
pulls f_g** at the fully-unanchored AMBIG node, whereas the 1-D at `w=0` discards the strand entirely and
uses *pure* `g_count`. The mixture has a ridge (f_g free with f₊=f₋ balancing); the count term competes
with the mixture's residual pull instead of governing outright. **To fully match the 1-D count fallback,
also down-weight the f_g-component of the strand mixture by `w`** (or scale β up as `w→0`), keeping the
f₊-vs-f₋ discrimination. This is the open design choice to settle in review (see §7.1/§7.2).

## 6. Validation plan (gates)

1. **Complex-loci battery** (`complex_loci_benchmark.py`, 49 loci): expect Family-1 loci (`full_single`,
   `anchor_sj_neither`, `nest_anchor_starved_core`) to drop to ≈ the 1-D error, **with no regression** on
   the 2-D wins (promoter, interleave, telescope, asymmetry families). Target: battery total 2D/1D < 1.0.
2. **Single-node controls** (`/tmp/node_mechanism.py` style): a single-strand pure-gDNA node must still
   give `f_g≈1` (Jeffreys still active there, `w→1`); an unanchored AMBIG node must give `≈ g_count`.
3. **Flagship + net-flow** (the standard suite, `use_propagation=True`): no regression vs the
   calibrated-Jeffreys sweep (the flagship has few AMBIG nodes, so this should be ≈ unchanged).
4. **Unit tests + goldens**: unchanged (the sweep is opt-in; production default `use_propagation=False`).

---

## 7. Open questions for review

1. **Gate the count term / strand mixture too? (NOW THE KEY QUESTION — the prototype residual.)** The
   prototype (§5b) shows gating *only* the Jeffreys leaves `full_single` partially over-called, because the
   strand-**mixture** term still pulls `f_g` at the unanchored AMBIG node. The 1-D at `w=0` uses *pure*
   `g_count` (no strand at all). Two ways to close it: (a) scale the **f_g-component of the mixture** by
   `w` too (keep f₊-vs-f₋ discrimination, drop the f_g pull where strand can't resolve it); or (b) scale
   the count precision up as `w→0` (`β_eff = β/(w+ε)` or `(1−w)·β_max`) so the count dominates the
   mixture's residual. (a) is more principled (it makes the sweep degrade to pure-count at `w=0` exactly
   like the 1-D); recommend (a), validate both.
2. **AMBIG with partial strand info.** A `κ`-stranded AMBIG node's strand carries info about `f₊` vs `f₋`
   (not `f_g`). Whichever option in #1, keep the f₊/f₋ discrimination — only the `f_g`-resolving part of
   the strand should be `w`-gated. Confirm the decomposition is clean on the lattice.
3. **Interaction with the propagation.** For an AMBIG node *with* odds inflow, gating the local Jeffreys
   should not weaken the propagation (the edge messages are separate). Verify on the promoter/interleave
   wins that they are unchanged. (Validation gate #1 covers this.)
4. **Does this generalize the per-node-β idea** the plan deferred? Effectively yes — this *is* a per-node
   strand-trust gate, applied to the prior rather than to a count blend.

---

## 8. Out of scope (follow-ups)

- **Family 2 — multi-transcript-per-strand / double-coverage** (`multistrand_*`, `deep_nest`,
  `nest_double_cover_inside`): the per-node mixture + odds-propagation assume one transcript per strand.
  Separate design doc.
- **Near-unidentifiable loci** (`anchor_intron_moat`, `anchor_exon_vs_intron_flank`): both solvers fail;
  these are a graceful-degradation question, not a 2-D regression.
