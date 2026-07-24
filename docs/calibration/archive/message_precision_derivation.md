# Message Precision in Rigel's Prior-Free Belief-Propagation Pass — v2

**A reviewer-ready derivation, diagnosis, and repair — the Compositional Evidence Compiler**

*Subject: whether the node-to-node message precision in `bp_solver._scan` manufactures
composition confidence on a composition-vacuous chain, and what the correct precision rule
is. v2 supersedes v1: it replaces the strand-only gate (v1 "Candidate B") with a
multi-channel reference-free evidence compiler, corrects the two adversarially-refuted term
definitions (`I_struct`, `I_motif`), states the identifiability ceiling honestly, and adds
the message↔hyperprior division of labor.*

> **Status.** The v1 diagnosis survives intact (§1–§2). The v1 fix does **not** ship as
> written: the owner and an external math reviewer corrected it, and successive adversarial
> verifications (reproduced with real numbers on `ambig_dense_10mb`) then corrected the
> reviewer, the owner's proposals, and each other. The net design is the **Compositional
> Evidence Compiler** of §3, with the caveats and open questions of §7–§8. Two structural
> reframes land in v2.1: **(i)** `I_motif` is a boundary→**exon EMISSION** message (§3.3) — the
> junction is a *source* that injects measured pure-RNA precision into its transcript-continuous
> exon flank, per-strand and one-directional — **not** a pin on the boundary's own `f_g` (the
> reviewer's "`pr_motif` into the boundary's own `τ`" is refuted: it re-commits `S ⊥ f_g`);
> **(ii)** the same junction is an RNA **SINK** in the reverse direction (§3.7) — as a
> *destination* it absorbs the exon's mature and relays nascent-only onward. The former
> `I_motif`-**mode** ship-blocker stays **RESOLVED** (§7.1, now reaffirmed against Derivation A's
> density-mode *revival*, which was adversarially refuted on the short-exon majority): `I_motif`
> ships **precision-only**, on the short-flank-safe `boundary_side_eff_length`. The
> phantom-collapse core and the `I_motif` precision seed **are** ship-ready. What is **not** yet
> ship-ready, flagged honestly: the onward-composition-precision `τ`-gate does **not** currently
> reach 0 on a fitted-`κ` unstranded chain (§7.8), so the sink's promised onward silence is a
> **goal with an unmet precondition**, not a delivered result. §10 records the reproduced
> end-to-end touchstone numbers; §11 is the Phase-B implementation plan.

---

## 0. Reading guide and headline

This document analyzes one sub-step of Rigel's calibration: the **prior-free
forward-backward belief-propagation (BP) pass** that deconvolves each genomic node's
unspliced fragment mass into a 2-simplex composition
`(f_pos, f_neg, f_g) = sense-RNA / antisense-RNA / gDNA`.

**Headline diagnosis (v1, survives adversarial review).** On an *unstranded* library the
pass manufactures **confidence** — collapsing per-node posterior variance by up to 35× to a
sharp answer — even though the chain carries **zero intrinsic composition information**
(strand Fisher info `∝ N(2κ−1)²`, identically 0 at `κ=½`). Message precision is sourced from
the running **belief** variance (which contains the shared reference prior and all pooled
inflow) instead of a reference-free **evidence** quantity. This corrupts the downstream
precision-weighted fit of the gDNA hyperprior (`1/var_g` weighting → circularity).

**Headline repair (v2).** Gate every message's composition precision on a **reference-free,
non-negative evidence count `τ`**, maintained separately from the belief and folded from
**all** legitimate composition channels — not strand alone. There are exactly three:

- **(a) STRUCTURE** — a signature-locked node (intergenic ⇒ pure gDNA; a strand switched off
  by signature) is a *composition-certain* measurement.
- **(b) SPLICING** — spliced reads at a boundary are *measured pure RNA*, a genuine
  composition signal even when unstranded.
- **(c) STRAND** — the per-strand tilt where `κ ≠ ½`.

`τ = 0` (⇒ emitted `pr = 0`) **only** when all three are silent. This subsumes v1's Candidate
B (strand-only), which the owner correctly rejected as *discarding* the load-bearing splice
channel (768 boundaries / 322k reads, §5).

**Headline honesty (the identifiability ceiling).** The first pass is a **phantom-confidence
gate + local-anchor propagator**, *not* a gDNA detector. On the decisive `gdna300` unstranded
case the messages currently move the gDNA-rich captured exons **away** from truth (mass-wtd
`f_g` local `0.509` → final `0.495`, oracle `0.853`) while fake-collapsing variance
`2.81 → 1.61`. Gating therefore **loses nothing** there. Enriched-island gDNA is
first-pass-**unidentifiable** and is grounded by the **Phase-2 stratified gDNA hyperprior**
(the lineage of released v0.7.1's background prior). §4 makes this division of labor explicit;
it is the single most important correction to the reviewer's proposal (their acceptance test
"detect real gDNA in the first pass" is unachievable and is reframed in §9).

**What resolved (the former ship-blocker).** The `I_motif` **precision** seed is derived and
sound and ships as the **whole** of `I_motif` — **precision-only**. The former open **mode**
question (how a measured mature-RNA density maps to an `f_g` floor) is now **closed by
rejection** (§7.1): a spliced count `S` carries **zero** honest information about the
*unspliced* composition the sweep solves (`S ⊥ f_g` under the nascent invariant), so **no**
composition mode-floor is in scope. Both the boundary self-floor `f_g ≤ N/(N+S)` and the
one-hop exon ceiling `f_g ≤ 1 − ρ_m/ρ_tot` are refuted with real numbers (the self-floor is
violated by **36.6 %** of `gdna300` boundary mass; the exon ceiling hard-clamps real gDNA
islands to `0` on short flanks). The residual gDNA magnitude belongs to **strand + the Phase-2
hyperprior**, never to `S`. **v2.1 sharpens the *framing* without reopening the mode:** the honest
content `S` carries is an **emission to the exon** (measured pure RNA is present and continuous with
the flank), delivered as *precision* on the exon's RNA channel (§3.3), not as a floor on anyone's
`f_g`; and the reverse direction — the junction as RNA **sink** — is added in §3.7.

---

## 1. The system and the governing invariant

### 1.1 Nodes, chain, latent parameterization

Each node accumulates `N` unspliced fragments split `s` (genome sense strand) / `a`
(antisense), `N = s + a`. We deconvolve the gDNA-vs-RNA split; per node the latent is

```
  λ = logit(f_g),   f_g = σ(λ),   f_r = 1 − f_g.
```

Nodes form a **linear bipartite chain** region ↔ boundary ↔ region ↔ …. The running BP belief
on `λ` is Gaussian `N(μ_λ, σ²_λ)`. The reported per-node posterior on `log f_g` has variance
`var_g` (a genuine grid second central moment, §2.3).

This is the **very first** step of calibration. It is *prior-free*: the gDNA composition
hyperprior (the expected DNA level, fit later as an enrichment shape over this pass's output)
is **not yet available**. The only structural prior present is the derived reference.

### 1.2 The reference prior and its exact moments (reproduced)

The reference is `Beta(1/2, 1/2)` (Jeffreys) on `f_g`, pushed through `λ = logit f_g` and
truncated to `λ ∈ [−10, 10]`:

```
  E_ref(f_g)        = 0.5000
  Var_ref(λ)        = 8.6743      ⇒  reference precision on λ:       a  = 1/8.6743 = 0.1153
  Var_ref(log f_g)  = 2.8037      ⇒  reference precision on log f_g: Λ_ref^ℓ = 0.3567
```

`Var_ref(log f_g) = 2.804` is the message-free "I don't know" belief width.

> **Note on the released code.** v0.7.1 did **not** use this `Beta(½,½)` reference in the
> message pass; it used the **background gDNA prior only**, with its precision **capped at 1
> pseudo-observation**. That cap is the crude ancestor of `τ`: `τ` *is* an honest evidence
> pseudo-observation count, 0 for an unpinned node and growing only with real strand / motif
> / structural Fisher information. The v2 fix generalizes the hard cap-at-1 into a
> per-channel evidence count.

### 1.3 The count-zero-information invariant (authoritative)

> A fragment **count** (or density) carries **zero intrinsic information** about a node's
> gDNA/RNA **composition**. 100 unspliced fragments split 50/50 across the two genome strands
> is equally consistent with pure gDNA (strand-symmetric) and pure RNA in an unstranded
> library. The count may enter only as **precision** (statistical power), never as a
> composition vote.

The **only** intrinsic per-node *unspliced-strand* composition likelihood is the **strand
tilt**: RNA tilts the genome-sense rate toward the library sense fraction `κ`; gDNA stays at
`½`. The sense count is `Binomial(N, p)`, `p(f_g) = ½·f_g + κ·(1−f_g) = κ + f_g·(½−κ)`, with
Fisher information for `λ`

```
  I_strand(λ) = N · (2κ − 1)² · [f_g(1 − f_g)]² / (4 p(1−p))   ≥ 0,               (1.1)
```

**identically 0 at `κ = ½`** (unstranded), for every count `N` and working point `f_g`.
Under Beta-Binomial overdispersion `N → N_eff = N/(1+(N−1)ω) → 1/ω`.

**Consequence (composition-vacuity).** In an unstranded library `Σ_i I_strand,i = 0`. A node
is then determined only by (a) **structure/signature**, (b) cross-node **imputation
messages**, and (c) the gDNA **hyperprior** (absent this pass). Any node ending the pass with
nonzero composition precision, absent a legitimate structure/splice/strand channel, is
transporting something that is **not evidence**.

The v1 update was to gate on **(c) strand only** — the owner corrects this: (a) structure and
(b) splicing are *also* legitimate, direct composition measurements and must not be discarded.

### 1.4 Why this matters downstream (the motivation, and the circularity)

After this pass, the gDNA hyperprior is fit as an **enrichment shape over per-node deconvolved
gDNA densities, precision-weighted by `1/var_g`**. Nodes ending confident (small `var_g`)
receive high weight. On the zero-gDNA diagnostic library the fitted hyperprior has L1 distance
≈ 2.0 to the oracle (near maximal) — corrupted by confident false positives. **The pass feeds
its own manufactured confidence into the hyperprior it is meant to inform.** Honest precision
out of the first pass is a prerequisite for a clean hyperprior; the gate is not an alternative
to the hyperprior but a **precondition** for it (§4).

---

## 2. The worked example (all numbers reproduced to the digit)

### 2.1 Condition

`gdna_none_ss_0.50_nrna_present_capture_on`: unstranded (`κ = 0.50117`, `(2κ−1)² = 5.5e-6`),
hybrid-capture, **zero-gDNA**. Ground truth `f_g = 0` at every node (2,677,093 RNA fragments;
max oracle `f_g = 0` over 1691 live nodes).

### 2.2 The current message-precision rule (`bp_solver._scan`, line ~446)

For a directed edge source `s` → destination `d`, the gDNA density message is a Gaussian
factor on `ℓ = log f_g` at the destination, with mode `mo` and precision

```
  pr = sm / ( sm·(v_logfg + s2t) + 1 )  ==  1 / ( v_logfg + s2t + 1/sm ),        (2.1)

    sm      = source's facing unspliced mass (a count; 1/sm → 0 for a big exon)
    v_logfg = (1 − f_g,s)² · σ²_λ,s        [delta Var(log f_g) from the RUNNING BELIEF]
    s2t     = var_proj[d] + (μ_proj[d] − μ_proj[s])²   [belief-free enrichment-transfer var]
```

Per-strand RNA messages are analogous with `v_logfr = f_g,s² · σ²_λ,s`. The load-bearing
defect: `v_logfg` is built from `σ²_λ,s` — the source's running **belief** variance
(line ~418), which contains the reference and all pooled inflow.

### 2.3 The chain, genomic order (oracle `f_g = 0` throughout)

```
 node  class   oracle | local_fg local_var | recv gDNA-msg pr | σ²_λ(fwd,bwd)   | final_fg final_var  sharpen
 2929  single   0.000 |  0.490    2.836     | 0.08             | (5.331, 5.160)  | 0.529    0.872      3.3×
 2930  bndry    0.000 |  0.490    2.812     | 0.16             | (9.329, 5.507)  | 0.757    0.675      4.2×
 2931  single   0.000 |  0.490    2.818     | 4.41             | (6.308, 0.521)  | 0.606    0.080     35.1×
 2932  bndry    0.000 |  0.510    2.778     | 2.76             | (1.976, 0.663)  | 0.510    0.136     20.5×
 2933  AMBIG    0.000 |  0.542    2.827     | 2.27             | (1.137, 0.988)  | 0.376    0.216     13.1×
 2934  bndry    0.000 |  0.510    2.772     | 1.74             | (0.876, 2.145)  | 0.451    0.293      9.4×
```

`var_g` is a true grid posterior second central moment of `log f_g`
(`simplex_logodds.py:356`), in the same space as `Var_ref(log f_g)=2.804`; grid-swept
40→240 and fold 17→65 leave it fixed — **not a quadrature artifact**.

### 2.4 Node 2931: a 35× sharpening to a confident-WRONG answer

Local: `f_g=0.490, var=2.818`. After the FB pass: `f_g=0.606, var=0.080` — a **35.1×**
collapse to a confident false positive (node 2931 has 119,244 RNA fragments and 0 gDNA),
built by messages alone.

### 2.5 Precision budget and cascade

`local_share = var_final/var_local = 0.080/2.818 = 0.028` ⇒ **97.2%** of the final precision
is message-supplied. Across all 1619 FP nodes the median message share is 0.726 (boundaries
0.781). Node 2931's dominant message (`pr = 4.315 = 1/(0.0836 + 0.1481 + 0.00009)`) comes
from a source (boundary 2932) whose own `σ²_λ=0.66` is 13× below reference and 95%
message-driven. Of 210 FP nodes receiving `pr>1`, **100%** have a source with running `σ²_λ`
below reference (median 1.99). **The confidence is built upstream and compounds.**

### 2.6 The GMRF pooling proof (the model defect, reproduced)

Trace the reference `π` through the sum-product messages. Each cavity contains `φ = π·L`; in
the vacuous case `L ≡ const`, so the belief is built **entirely from transported copies of
`π`**, stacking like independent pseudo-observations `λ ≈ 0`. Exact GMRF marginal
`Var(λ) = 1/√(a(a+4b))`, `a = 0.1153`, `b = 1/σ_e²`:

```
   σ_e²   3-node Var(λ_2)   sharpen    401-node interior Var(λ)   sharpen
   1.00      3.105           2.79×          1.452                  6.0×
   0.50      3.001           2.89×          1.034                  8.4×
   0.15      2.925           2.97×          0.569                 15.2×   ← code regime
```

**Decisive control:** `a → 0` (drop the per-node prior precision, keep smoothness) makes the
precision matrix **singular** — the marginal is improper. Every bit of manufactured confidence
is pooled **reference precision**; none is evidence. And this pooled precision is exactly the
`1/var_g` currency weighting the downstream hyperprior — the **circularity**.

### 2.7 Representativeness

1619 of 1691 live nodes (81.2% of mass) are false positives; 290 sharpen >10×; median 3.54×.
The same direction/magnitude appears in the nonzero-gDNA capture library (single 4.3×, AMBIG
4.2×, bndry 6.1×): the mechanism is the precision model itself, independent of true gDNA.

---

## 3. The Compositional Evidence Compiler (the v2 fix)

Maintain, **separately** from the belief `N(μ_λ, σ²_λ)` used for point estimates, a
**reference-free composition-evidence precision** on **both** latent axes — `τ^λ` on
`λ = logit f_g` and `τ^θ` on the strand-split — folding **only** genuinely
composition-informative Fisher information from all three legitimate channels.

### 3.0 The compiler

```
  τ^λ_i = I_strand,i  +  I_struct,i  +  I_motif→λ,i  +  Σ_{k∈N(i)} pr^λ_{k→i}      (3.0a)
  τ^θ_i = I_strand^θ,i +                                Σ_{k∈N(i)} pr^θ_{k→i}      (3.0b)
```

The belief precision is `τ + a` — **the reference is added once, locally, and never
transmitted.** Messages are emitted from the **cavity** evidence (exclude the return edge,
exclude the reference), with the convention `τ = 0 ⇒ v_evid = ∞ ⇒ pr = 0`. Each term below is
given its count-zero-info-clean definition; the two axes combine at emission per §3.5.

### 3.1 I_strand — the strand tilt (channel c)

The only intrinsic *unspliced-strand* composition Fisher information, eq. (1.1):

```
  I_strand,i = N_eff (2κ − 1)² [f_g(1 − f_g)]² / (4 p(1−p)),   p = κ + f_g(½−κ).
```

The count `N` enters only as multiplicative power. `≡ 0` at `κ = ½`. Verified: `κ = 0.9,
N = 100 ⇒ I_strand = 4.76` at `f_g=½`; `κ=½ ⇒ 0.00000` for `N` up to `1e5`. The parallel
`θ`-axis term `I_strand^θ` is the same tilt resolved onto the per-strand split; it carries the
`κ≠½` information into `τ^θ` so the strand channel is not laundered back through `Var(θ)`
(without it the per-strand RNA message precision climbs `0.35 → 7.69` on a vacuous chain, the
θ-analogue of the λ pooling — §5).

### 3.2 I_struct — structural certainty (channel a), CORRECTED to a count-free gate

A signature-locked node (intergenic ⇒ pure gDNA `f_g=1`; a strand switched off by signature)
has its composition fixed **by structure, not by count**. The reviewer proposed
`I_struct = N_eff`. **Adversarial verification (reproduced on `gdna300`) refutes the `N_eff`
scaling as provably inert, and shows the clean form is a boolean gate:**

The emitted precision is `pr = 1/(v_evid + s2t + 1/M_src)`, `v_evid = (1−f_g,src)²/τ_cav`. A
true lock has `f_g,src = 1` **exactly** (`init_beliefs`: `f_g=1`, `var_g=0` for G1), so
`(1−f_g,src)² = 0` and `v_evid = 0` **regardless of `τ`**. The magnitude of `τ` — `N_eff`, 1,
or the reference — **cannot enter**.

> **Vertex-annihilation proof (real numbers, `gdna300_ss0.50_capture_on`,
> `(2κ−1)²=2.5e-6 ⇒ I_strand≡0`).** Over **611 forward + 635 backward** edges whose source is
> a `gonly` lock (median `f_g,src = 1.0000`): median `v_logfg = (1−f_g,src)²·vls = 0.000000`.
> Holding `f_g, s2t, M_src` fixed and sweeping the source variance across its full range, the
> emitted `pr_g` is **byte-identical**:
>
> ```
>   dir   n_edges   pr@vls=ref(8.674)   pr@vls=0   pr_actual   med s2t   med 1/M_src
>   fwd     611          2.084           2.084       2.084       0.123      0.203
>   bwd     635          2.067           2.067       2.067       0.123      0.217
> ```
>
> `I_struct`'s magnitude is **annihilated**. Its only function is to flip the `0/0` at the
> vertex from the convention `τ=0 ⇒ pr=0` to `composition-certain ⇒ v_evid=0 ⇒
> pr = 1/(s2t + 1/M_src) > 0`. **A boolean bit, not a count.** The lock emits at `pr ≈ 2`
> governed **entirely** by the density message's honest `1/M_src` (median `M_src ≈ 5`, so
> `1/M_src ≈ 0.20` dominates) and `s2t` — count-zero-info-compliant by construction.

**Corrected definition.**

```
  I_struct,i  ≡  "composition-certain" gate: τ^λ_i = +∞  (equivalently var_composition = 0)
                 when node i is signature-locked (gonly / a switched-off strand),  else 0.
                 It injects NO count-scaled composition precision. The count precision flows
                 ONLY through the density message's 1/M_src + s2t.
```

Off-vertex (a hypothetical `lock` at `f_g < 1`) the `N_eff` form would inject
`(1−f_g)²/N_eff` — count-scaled composition precision unmediated by any density-smoothness
term, a direct §0 violation. The boolean gate gives `v_evid = 0` there too (a *certain*
composition, whatever its value, emits at full density-message precision, zero composition
precision injected) — clean in both cases. **Drop "must scale with `N_eff`".**

**Over-imputation.** Because the vertex annihilates `τ`, `I_struct` adds no new
over-imputation risk. The residual density-message risk is `pr → 1/s2t` as `M_src → ∞`,
bounded solely by `s2t` (the enrichment-crossing transfer variance, orthogonal to
`I_struct`). Data: locks are low-count (median `M_src≈5`, `1/M_src≈0.20` caps `pr≈2`) and
directionally **safe** toward enriched exons (a lock carries *background* gDNA density; a
capture-enriched exon has *locally elevated* gDNA, oracle `f_g=0.853`, so a lock if anything
**under**-imputes there).

### 3.3 I_motif — the boundary→exon spliced EMISSION (channel b): RNA-channel-only, additive, per-strand, one-directional

A spliced count at a boundary is a **direct, pure-RNA composition measurement** (a spliced
read cannot be gDNA). It is the load-bearing unstranded anchor: **764 / 749 boundaries** carry
spliced mass (**767,718 / 763,520 reads**; **322,761 / 321,485 mass**) in `gdna300` / `gdna5`
respectively — library-invariant and abundant, vs the sparse `I_struct` substrate (§4). Gating
on strand alone (v1) would discard it.

**The emission frame (the owner's authoritative correction).** `I_motif` is **not** a statement a
junction makes about *its own* `f_g`. It is a **message the boundary EMITS to its
transcript-continuous EXON flank**: "there is a definitely-measured, guaranteed-pure-RNA
population continuous with you, at a known count level, that does not itself need to be solved."
The junction is the **SOURCE**; the exon is the **DESTINATION**. This directionality is what keeps
the emission honest under `S ⊥ f_g`:

- **It does NOT pin the boundary's own unspliced `f_g`.** The reviewer's proposal to fold
  `pr_motif` into *the boundary's own* `τ^λ` is **refuted**: the boundary's unspliced pool is
  `nascent + gDNA` (mature has spliced out into the *disjoint* spliced pool), so `S` carries
  **zero** information about the boundary's own unspliced composition (`S ⊥ f_g`). A self-`τ`
  there would flag the boundary high-precision and **block Phase-2** from grounding node 2930 from
  `~0.5` to its true `~0.007`. The boundary's own `f_g` is set by **strand + λ-imputation +
  Phase-2 only** — never by its own `S`.
- **It is per-strand and one-directional.** A junction is single-stranded
  (`boundary_substrate.junction_strand ∈ {+, −}`, the observed genomic motif) and its spliced mass
  sits on **one face** — the same-strand exon flank (`node_geometry._spliced_faces`,
  `node_geometry.py:170`: `spliced_pos_*` fires iff `js==+` **and** that flank carries
  `BIT_EXON_POS`). So `pr_motif` enters **only** `τ^{f_pos}` (or `τ^{f_neg}`) of the **exon on that
  face** — never the opposite strand, never the intron flank, never symmetrically. In the
  forward–backward sweep an *acceptor* face (mass on the right, `_scan sf=1`) fires in the
  **forward** pass only and a *donor* face (mass on the left, `sf=0`) in the **backward** pass
  only, so each motif face contributes to exactly one direction and the FB `_comb` never
  double-counts it. Verified: node 2930 `js=+`, `spliced_pos_right=7501.8`, all other faces 0,
  `spliced_n_pos_right=15753 ≈ 2·mass`; donor node 132 `spliced_pos_left=2709`; neg junction node
  274 `js=−`, `spliced_neg_left=3233`.

The routing above is **already fully realized by the shipped geometry** (`_spliced_faces`); Phase-B
changes only the *precision consumption* in `_scan`, not the geometry. **Three adversarial
corrections to the reviewer's `I_motif = S_eff·h(f_r)`:**

**(1) The Jacobian `h(f_r) = f_g²` is correct — as an RNA-channel-only term.** A spliced count
`S` is a Poisson measurement of the boundary's mature-RNA log-density `log ρ_m = log(S/E_spl)`,
Fisher info `S_eff`. On the RNA-total factor (a factor on `log f_r`; Jacobian
`d log f_r/dλ = −f_g`), setting `v_logfr = f_g²/τ = 1/S_eff` yields
`h(f_r) = f_g² = (1−f_r)²`. Verified emitted RNA precision `= S_eff = 7502`.

**(2) It must be FORBIDDEN from the gDNA channel.** A single `τ^λ` also feeds the gDNA factor
as `v_logfg = f_r²/τ ⇒ pr_g = f_g²/f_r²·S_eff`. Reproduced: **`pr_g = 86,965`** (gdna300) /
`85,329` (gdna5) — a phantom gDNA precision from a **pure-RNA** junction. A spliced read
measures RNA presence, never gDNA background (both coexist). So `I_motif` is an **RNA-only**
term `τ^{f_r}` and one scalar `τ^λ` cannot express that.

**(3) It must be a separate ADDITIVE factor, and BYPASS `s2t`.** In the current code the same
7502-read pure-RNA measurement is **zero-counted**, not double-counted: absorption drives the
incoming residual negative, the honest clamp zeros `n_eff`, and RNA precision **collapses
`119,007 → 0`** (reproduced at node 2930; its own 7502-read measurement yields zero
precision, mis-solving `f_g` to `0.771` vs oracle `0.007`/`0.305`). The fix is **precision
addition** of two disjoint factors:

```
  pr_RNA = pr_imputation + pr_motif,
    pr_imputation : current honest-clamp / absorption path (mode + source-count; legitimately zeroable)
    pr_motif      : f_g² · S_eff, seeded from the RAW spliced count, NOT from rho_pos.
```

A zero in `pr_imputation` can never subtract `pr_motif` — the **algebraic can-never-cancel
guarantee**; node 2930 restores `pr_motif ≈ 7502`. `pr_motif` also **bypasses `s2t`**: the
boundary's own spliced is a **measurement AT the node** (owner's channel b), not an imputation
across an edge, so no crossing-variance applies (at node 2930 `s2t ≈ 3.9–29.6` would otherwise
crush a correctly-seeded `I_motif` to `pr ≈ 0.03–0.3`).

**`S_eff` choice.** `S_eff = Kish n_eff = (Σw)²/Σw²`. The geometry plumbs the mass `Σw`
(`spliced_* = 7502`) and integer count `N` (`spliced_n_* = 15753 ≈ 2·mass`, half-triangle
partial deposits) but not `Σw²`. Use `S_eff = MASS Σw` (**conservative lower bound**, never
over-confident); `N` is the 2×-over-confident upper bound. Exact Kish needs one extra
accumulator reduction `Σw²`. Because `pr_motif` ships **precision-only**, `S_eff` is the raw
mass `Σw` directly — **no eff-length divides it** (an eff-length enters only when a *mode* is
placed, and none is; §3.3(4)). If a mode-scale were ever wanted it must use the short-flank-safe
`boundary_side_eff_length = E[min(ℓ,R)]/2`, never the half-triangle `spliced_side_eff_length`
(2–199× low on short flanks — `effective_length.py` docstrings).

**Double-counting is orthogonal (precision vs mode), only in the target design.** One
measurement `S` contributes both moments to a **single mature-RNA factor** (mode `log ρ_m`,
precision `S_eff` — a Gaussian factor legitimately owns both). The absorption is a **mode**
operation on a **disjoint** factor (the incoming nascent imputation), subtracting the
already-measured mature mean so `nascent = incoming − mature` is not re-added — and it must
touch **no** precision. In the current code the mature self-factor does not exist, so its
precision has nowhere to live but the incoming message, where the clamp zeros it — the
`119007 → 0` coupling.

**(4) The COMPOSITION mode is precision-only — the emission carries NO `f_g` floor (resolves the
§7.1 ship-blocker, reaffirmed against Derivation A).** `pr_motif` contributes **precision** to the
**exon's** RNA channel; the composition **point** stays the ordinary RNA imputation mode (the
exon flank's own belief), **never** a data-derived `ρ_m/ρ_tot` ceiling on `f_g`. A spliced count is
`S ⊥ f_g` about the *unspliced* composition the sweep solves
(nascent-cannot-be-inferred-from-abundance: `π(ρ_g,ρ_r,ρ_m)` factorizes, so the `ρ_m`-integral
cancels in `p(f_g|·)`), so any floor `f_g ≤ g(S)` is a smuggled nascent-from-abundance prior.

Derivation A attempted to **revive a mode** — a *soft, one-sided* density-mode lower bound
`m_r = log(ρ_m/ρ_tot)` on the **corrected** `ρ_m = S/boundary_side_eff_length`, penalising only
`f_r < ρ_m/ρ_tot`. It is **honest on the long-exon touchstone** (exon 2931, 1273 bp, a genuine
73,639-fragment *unspliced-mature* population: `ρ_m = 75.0 ≈ ρ_mat_uns = 68.6`, delivering the
honest ceiling `f_g ≤ 0.32` from the uninformed `0.49`). But it was **adversarially refuted on the
short-exon majority** and **must not ship as a mode**: on `gdna300`, **52/517 junction faces
(25.9 % of all spliced mass)** assert `ρ_m` *greater than* the exon's true unspliced-RNA density;
**31 land on ~pure-gDNA short exons** (114–178 bp, oracle `f_g=1.0`, `mat_uns=0`) yet inject
`ρ_m=39–85` → a dishonest lower bound on the very unspliced axis the sweep deconvolves. The reason
is decisive: the fraction of a transcript's mature coverage appearing as *unspliced* fragments in
an exon body is `P(unspliced | R) → 0` for `R ≪ FL`, so on a short exon **all** the mature is in
the disjoint spliced pool (`mat_uns=0`) and `S` carries *zero* info about the exon's *unspliced*
`f_g` — exactly `S ⊥ f_g`, now at the exon. The honest mode `log(ρ_m·P(unspliced|R)/ρ_tot)` needs
the length factor `P(unspliced|R)`, which reintroduces a length threshold (a magic number). **Net:
the mode has no regime that is simultaneously honest AND load-bearing** — long exons: honest but
*redundant* with the ordinary unspliced imputation, which already sees the 73,639 `mat_uns`; short
exons: dishonest.

So `I_motif` ships **precision-only**: the emission *raises the exon flank's RNA precision* (the
exon genuinely has continuous mature RNA), riding the ordinary imputation mode, and the residual
exon gDNA/RNA split is left to **strand + the Phase-2 hyperprior**. `pr_motif` adds a **magnitude**
on the RNA/`f_r` axis with **no mode of its own**; `ρ_m` is used only (if at all) to *scale*
`S_eff`, never to place a ceiling.

### 3.4 Cavity emission

The outgoing message uses the cavity evidence (exclude the destination edge, exclude the
reference), non-negativity guaranteed since `τ` is a sum of non-negatives:

```
  τ^{λ,cav(j)}_i = I_strand,i + I_struct,i(gate) + I_motif→λ,i + Σ_{k≠j} pr^λ_{k→i},        (3.4)
  (analogously τ^{θ,cav(j)}_i for the strand-split axis.)
```

### 3.5 The two-axis v_evid mapping at emission

The per-strand RNA message precision draws from **both** axes; the destination folds `pr` into
**both** the belief (point estimate) **and** `τ` (evidence channel):

```
  v_evid,pos = (1 − f_g,s)² · (1/τ^{λ,cav}_s)  +  (cos θ /(1 + sin θ))² · (1/τ^{θ,cav}_s)   (3.5)

  pr_{i→j} = 1 / ( v_evid + s2t + 1/sm ).
```

**If both `τ^λ` and `τ^θ` are 0** (unstranded, no structural/motif anchor) ⇒ `v_evid → ∞` ⇒
`pr → 0` on **both** channels. This is the phantom-collapse: node 2931's source (unstranded
boundary 2932, no lock, no motif) has `τ^λ = τ^θ = 0 ⇒ pr_{2932→2931} = 0` (vs current 4.315).
The `I_motif` and structural terms keep the honest anchors alive.

For hygiene, replace the delta factor `(1−f_g)²/τ` by the exact
`Var[log f_g | N(μ, 1/τ^{cav})]` (v1 §3.5; a ≤7% correction at mid-range, up to 2.5× near the
vertices — adopt as the variance transform, do not gate on it).

### 3.6 Limiting checks

- **(i) Vacuous unstranded chain → 0.** `I_strand = I_strand^θ = 0`, no lock, no motif ⇒ every
  `τ = 0`, every `pr = 0` by induction. Node 2931 → `pr = 0`, `f_g ≈ 0.49`, `var ≈ 2.82`; the
  35× sharpening eliminated at the root.
- **(ii) Genuine strand anchor propagates.** `κ ≠ ½ ⇒ I_strand > 0` seeds `τ`; the anchor's
  own marginal variance is exactly `1/I_strand` (gains nothing spurious), attenuating `+s2t`
  per hop.
- **(iii) Structural anchor propagates one hop where substrate exists.** A lock emits
  `pr = 1/(s2t + 1/M_src) > 0` to a directly-adjacent open node (governed by count+transfer,
  not `τ` magnitude — §3.2).
- **(iv) Motif anchor propagates.** A spliced boundary emits `pr_motif = f_g²·S_eff > 0`,
  pinning adjacent exon flanks to RNA; can-never-cancel by a negative absorption residual.
- **(v) Non-negativity / soundness.** `τ ≥ 0`, `v_evid ≥ 0`, `pr ≥ 0`. No clamp.

### 3.7 The reverse direction — the junction as RNA SINK (absorption + onward-nascent)

The same junction that EMITS to its exon (§3.3, source direction) is an RNA **SINK** in the reverse
direction (destination direction): it receives the exon's RNA message, absorbs the mature it has
itself measured, and relays **nascent-only** onward to the next region.

**Physics.** Every population crossing a junction seam partitions cleanly: **mature** reads
*splice* across (counted as spliced `S`, one-sided, motif strand — they never cross *unspliced*);
**nascent** (pre-mRNA) crosses *unspliced*; **gDNA** crosses *unspliced*, strand-symmetric. So the
junction's **unspliced** pool is `nascent + gDNA` only. The exon, however, sends its *total*
unspliced-RNA density, which mixes the exon-internal **mat_uns** (mature fragments that merely
never crossed a junction) with **nas_uns**. The sink's job is to **remove the mature part** so what
relays onward is nascent-only — this is the `_RNA_ABSORB` subtraction, and it is physically
required.

**Decouple MODE from PRECISION (the one coherent mechanism).** The current `_scan` couples them:
`rho_nascent = ρ_r^incoming + SPs/esp − SPd/ESPd` (mode) and then derives the onward precision from
that same residual via `n_eff` + the belief-variance `v_logfp`. Split into two orthogonal gates:

- **MODE gate = the absorption (kept in spirit, fixed in frame).** `rho_nascent =
  ρ_r^incoming − ρ_mature^own`, floored at `1/erd`. `residual ≤ 0 ⇒ nothing nascent continues ⇒
  drop the RNA message`. Purely a *mode* decision ("is there any nascent left?"), carrying **no**
  precision claim. **Frame fix (required):** the subtraction must be on the boundary's **own
  crossing eff-length** — `boundary_side_eff_length = E[min(ℓ,R)]/2` — the **same** frame as the
  nascent density it subtracts from, **not** the half-triangle `spliced_side_eff_length` the code
  uses as `ESPd` (`bp_solver.py:538`). On a long flank the two agree (node 2930: 100.0 vs 99.99);
  on a short flank the half-triangle is **2–199× low** (177 short faces on `ambig_dense_10mb`),
  inflating the measured mature density and manufacturing a **spurious** negative residual → a
  spurious "no nascent" → the onward nascent wrongly killed. Like-for-like density is the invariant.
- **PRECISION gate = `τ` (the compiler).** The onward *composition* precision is
  `pr_onward = 1/(f_·²/τ^{cav} + s2t + 1/M_src)`, `= 0` when `τ^{cav}=0`. The absorbed mature count
  `S` contributes **exactly zero** to it: `S ⊥ (f_g, f_nascent)` under the nascent invariant, so `S`
  carries no Fisher information about the onward `λ`. The onward cavity evidence is therefore
  `τ^{λ,cav} = I_strand + Σ_{k≠dst} pr` — **no** `I_motif` term (the intron is not this junction's
  exon flank, §3.3), **no** `I_struct` (a spliced junction is never signature-locked:
  `S>0 ⇒ free_s=True ⇒` solvable, not a `{0,0,1}` vertex). On a vacuous unstranded chain
  (`I_strand=0`, no inflow) the target is `τ^{λ,cav}=0 ⇒ pr_onward=0`: the boundary relays its
  point-estimate **mode** as a chain anchor but injects **zero** composition precision, deferring
  its onward gDNA-vs-nascent split — and its own — to Phase-2.

**The no-compensating-gDNA-vote rule (closes the sink-laundering phantom; reconciles §6 / v1
open-Q#8).** Today a *zeroed* RNA residual sets `n_eff=0 ⇒ pr_RNA=0`, which pie-complementarity
(`f_g = 1 − f_r`) then launders into a **confident gDNA vote** — the `f_g ≈ 0.77`/`0.97` phantom at
node 2930. The rule: `residual ≤ 0 ⇒ emit NO RNA message AND add NO gDNA composition precision from
the RNA absence` (`τ += 0` on both). A fully-absorbed junction emits `≈0` nascent **and** `≈0` gDNA
composition precision *coherently*; the honest clamp's `n_eff→0` and the "no gDNA vote" are the
**same** `τ=0` — one weak-zero, not two.

**Honest caveat (the onward `τ`-gate does not yet reach 0 — §7.8).** The onward-silence target above
is a **goal with an unmet precondition**. Reproduced on the touchstone with the shipped
`_TAU_PRECISION` machinery: the fitted `κ = 0.500975 ≠ ½` makes `I_strand ≠ 0`, which *clears* the
hard `τ > 1e-9` emit gate (`bp_solver.py:491`); the cavity accumulation `tau_lam[i] += Σpr`
(`bp_solver.py:585`) is then a positive-feedback loop (`τ↑→1/τ↓→v_logf↓→pr↑→τ↑`) that self-inflates
`τ` ~6 orders of magnitude to `τ≈0.84` at the exon — so the onward `pr_g≈0.041` is **identical**
with `τ` on vs off, and node 2930 stays at the `0.76` phantom. **Under capture-off the same formula
gives `pr_g≈4.8`** — the full cascade; the `s2t`/capture term is currently the *only* thing holding
the onward precision down. So the sink's **mode-op + frame-fix + no-vote rule** are ship-ready; the
**onward-composition-precision `τ`-gate is NOT** until §7.8 is resolved.

---

## 4. Division of labor — the message pass vs the stratified gDNA hyperprior

This section is the most important v2 addition and the reframe the reviewer's proposal was
missing. **The first-pass message layer is a phantom-confidence gate + local-anchor
propagator, not a gDNA detector.** Enriched-island gDNA is deferred to the Phase-2 stratified
hyperprior *by design*.

### 4.1 The identifiability ceiling (honest, reproduced)

On `gdna300_ss0.50_capture_on` (unstranded, true live gDNA frac **0.691**), the gDNA-rich open
exons (oracle `f_g > 0.5`) are **59.3% of live mass** (n=1165, mass-wtd oracle **0.853**), yet
the current messages move them **away** from truth:

```
  cut        n     massfrac  oracle  local_fg  final_fg  |loc−or|  |fin−or|  var_loc  var_final
  fo>0.5   1165     0.593     0.853    0.509     0.495     0.344     0.361      2.811     1.611
  fo>0.7   1016     0.486     0.906    0.510     0.505     0.396     0.401      2.811     1.607
```

The messages push toward **RNA** (down from local 0.509) while **fake-collapsing** variance
`2.81 → 1.61`. This is first-pass **unidentifiable** by count-zero-information: strand is
silent (`I_strand = 0`), the count carries zero composition info, and the only reachable local
structure is a spliced boundary saying "RNA". A captured exon at true `f_g = 0.85` is
indistinguishable from one at 0 using strand + local structure alone. **Gating loses nothing
here** — the current "detection" is negative.

### 4.2 Why `I_struct` cannot rescue it (substrate-limited, reproduced)

`gdna300` has **216** live intergenic `gonly` locks (**4.1% of mass**), **zero** 1-hop
adjacent to any gDNA-rich exon. BFS: 206/216 reach one at **median 2 hops** (max 4). But even
at 2 hops the lock supplies **background** gDNA density, not the exon's **elevated
capture-gDNA**; and 288 of those 2-hop paths cross an intervening spliced boundary whose
`I_motif` pushes the exon toward RNA — actively **opposing** the lock. Structural imputation
from sparse background locks provably cannot detect capture-enriched gDNA in the first pass.

### 4.3 `τ` as the arbitration key

The first pass produces two decoupled per-node outputs: a point estimate (mode) and an
**honest evidence precision `τ`**. `τ` is the mature descendant of v0.7.1's cap-at-1
pseudo-obs — 0 for unpinned nodes, growing only with real Fisher information. It is the key
that arbitrates the two stages:

- **high-`τ` nodes** are message-decided (the hyperprior must not override them);
- **low-`τ` nodes** (the 59%-mass unpinned enriched exons) are **hyperprior-decided**.

The Phase-2 stratified hyperprior owns exactly what the messages cannot: the aggregate DNA
background scalar `ρ_bg` (pooled from pure intron/intergenic — precise even at true zero) ×
the enrichment **shape** grounds the low-`τ` enriched islands (`ρ_g = ρ_bg · enrichment`).
Without the gate this arbitration is **corrupt**: the pass hands the hyperprior a
fake-confident `1/var_g` weight (pooled reference precision, hyperprior L1 ≈ 2.0 to oracle) —
the circularity of §2.6. **The gate is a prerequisite for the hyperprior, not an alternative
to it.** Legitimate enriched-gDNA "detection" lives in Phase 2, fed the clean `τ` substrate.

---

## 5. Real-number honesty on the refined design

- **`I_strand` — zero unstranded.** `κ = 0.4992` (gdna300) / `0.5010` (gdna5) /
  `0.50117` (gdna_none) ⇒ `(2κ−1)² = 2.5e-6 / 3.8e-6 / 5.5e-6` ⇒ `I_strand ≈ 0` **everywhere**
  in the diagnostic suite. The strand channel is silent; the fix stands or falls on structure
  + splicing.
- **`I_struct` — sparse, thin reach.** 216 locks, 4% mass, 0 one-hop gDNA-rich adjacency,
  median 2 hops, wrong density level. Right for grounding the intron/intergenic neighbourhood;
  **substrate-limited** on enriched captured exons. Will not, and should not be asked to, solve
  `gdna300`.
- **`I_motif` — abundant, load-bearing.** 749–764 boundaries, 322k mass, 763–768k reads,
  library-invariant. The dominant unstranded anchor in every regime. This is why gating on
  strand alone (v1) is wrong.
- **Both regimes.** `gdna_none`/`gdna5`: messages *improve* the mass-weighted mean
  (`0.508 → 0.386` toward oracle 0) but fake-collapse variance and corrupt the precision
  channel; gating loses the mean-nudge, and the background hyperprior (low `ρ_bg`) then grounds
  those nodes honestly. `gdna300`: messages move enriched exons *away* from truth (§4.1);
  gating strictly improves them (`|err| 0.361 → 0.344`, honest `var ≈ 2.8`). In **both** the
  first pass carries honest evidence only and the stratified hyperprior does the grounding.

---

## 6. Answering the reviewer's closing question — `I_struct` (structural precision) vs the Phase-1b boundary-subtraction bug

> *"Given that structural evidence must flow into adjacent exons, how do we keep the
> downstream boundary subtraction (Phase 1b) from accidentally cancelling true structural
> evidence with local negative splicing calculations?"*

**Answer: the two never touch — they are orthogonal channels on orthogonal quantities. Three
rules make it airtight, and they reconcile v1 open-Q#8.**

1. **`I_struct` is a PRECISION on `τ^λ` — a sum of non-negatives.** No splicing subtraction can
   drive it negative. The absorption term
   (`rho_pos = f_s·sm/E_r + SP_s/E_spl − absorb_p`) acts only on the **RNA-channel density
   MODE**, never on `τ`. Structural certainty survives even when a local RNA residual momentarily
   goes negative.

2. **Structural gDNA rides the gDNA `λ`-face, which has NO subtraction.** In `_scan` the gDNA
   message (`rho = f_g·sm/E_g`, `pr = sm/(sm(v_logfg + s2t)+1)`) has no absorption term; the
   RNA-residual subtraction is RNA-channel-only. Route `I_struct`'s precision onto the gDNA `τ`
   axis; **never net it against an RNA density residual.** Structural gDNA evidence is
   structurally immune to the Phase-1b bug.

3. **The two fixes are the SAME quantity (reconciling v1 open-Q#8).** The absorption
   honest-clamp already sets `n_eff` to the residual's **own** count (a weak zero when a
   boundary fully absorbs its incoming RNA into measured spliced). That `n_eff` **is** the `τ`
   contribution of that RNA edge: a fully-absorbed junction emits `≈0` residual **and** `≈0`
   imputation precision **coherently** — one mechanism, they cannot double-suppress. The bug is
   only that today a *cancelled* residual is laundered by pie complementarity into a confident
   gDNA vote (`f_g → 1`, the `f_g≈0.97` cascade). **Rule: `residual ≤ 0 ⇒ emit NO RNA message
   AND NO compensating gDNA vote` (`τ += 0`), and do the spliced subtraction on the boundary's
   OWN eff-len (like-for-like density).** With the additive `pr_motif` seed (§3.3) the
   boundary's raw 7502-read measurement is a *separate* standing factor, so even a fully-absorbed
   junction keeps `pr_motif ≈ 7502` — the structural/motif evidence can never be cancelled by
   the local subtraction.

---

## 7. Open questions and caveats for the reviewers

§7.1–§7.2 are now **resolved** (the former `I_motif`-mode ship-blocker); §7.3–§7.7 we do
**not** consider closed.

### 7.1 `I_motif` MODE — RESOLVED (was the ship-blocker): no composition mode-floor, precision-only

The precision seed `pr_motif = f_g²·S_eff` is sound and ships as the **whole** of `I_motif`.
The open **mode** question is closed **by rejection**: there is no honest composition
mode-floor from a spliced count.

**Boundary-only (the reviewer's region confusion).** Only **boundary splice-junction nodes**
carry a spliced count `S` — a spliced read requires a motif-stranded donor–acceptor junction to
align across. **Region nodes, introns included, have `S = 0`.** The reviewer's worked example
applied the floor to an **intron** (a region), where it degenerates to `N/(N+0) = 1` — no
constraint at all. The essence is boundary-only, and even there it does not do what the
reviewer claimed.

**The floor is on the wrong pool (decisive numbers).** The sweep solves the **unspliced**
composition `f_g = N_g/N`; the reviewer's `f_g ≤ N/(N+S)` is a bound on the **combined**
`(N+S)` composition. Both readings are dead:

- **Combined:** `f_g,comb = N_g/(N+S) ≤ N/(N+S)` is **hard but trivial** — it is the counting
  identity `N_g ≤ N` in disguise (**0 violations** in `gdna_none/gdna5/gdna300` **by
  construction**), exactly achievable, and therefore never binding.
- **Unspliced:** applied to the `f_g` the sweep actually solves, the floor is **violated by
  36.6 % of `gdna300` boundary mass** (394 boundaries) — fast-splicing, gDNA-rich junctions
  whose unspliced crossing is nearly all gDNA *despite* abundant `S`. This is the
  nascent-from-abundance violation: `S ⊥ f_g` under the nascent invariant
  (`π(ρ_g,ρ_r,ρ_m) = π(ρ_g)π(ρ_r)π(ρ_m)` ⇒ the `ρ_m`-integral cancels in `p(f_g | ·)`), so any
  `f_g ≤ g(S)` smuggles a `π(ρ_m,ρ_r)` cross-term. A hard unspliced cap pushes real gDNA
  **below truth** and **blocks Phase-2 recovery** on exactly the mass Rigel most needs.

**The one-hop exon ceiling `f_g ≤ 1 − ρ_m/ρ_tot` also fails (adversarially reproduced).**
Displacing the floor one hop to the transcript-continuous exon flank (`ρ_m = S/E_spl`) is
invariant-clean in *form* but **not parameter-free in practice**: `E_spl` is the half-triangle
`spliced_side_eff_length`, which `effective_length.py` flags **2–199× low on short internal
exons**, so `ρ_m` over-states the local RNA density on short flanks (`ρ_m /` true mature
density p90 ≈ 5e11). Result: **32 `gdna5` exon faces hard-clamped to `f_g = 0`** by a *negative*
barrier, **193 `gdna300` faces ceiled below oracle**, and **binding on real gDNA islands**
(boundary 640 → exon 639, `S_eff = 2233`, caps `f_g` at `0.645` against oracle `1.000` — a
pure-gDNA capture island; Derivation 1's `(f_g=1, S>0)` fast-splice counterexample realized at
the exon). `ambig_dense_10mb`'s long exons (median 440–495 bp) hide this at ~0.1 % mass; real
transcriptomes (median internal exon ≈ 120 bp ≪ FL) make short flanks the **majority**. The
`ρ_m/ρ_tot` mode is `E_spl`-geometry-dependent and re-commits the over-cap one hop out — **do
not ship it.**

**Derivation A's revival (one-sided + corrected eff-length) also fails — the short-exon
`P(unspliced|R)→0` wall.** A later attempt hardened the one-hop ceiling into a *soft, one-sided*
density-mode lower bound `m_r = log(ρ_m/ρ_tot)` on the **corrected** short-flank-safe
`ρ_m = S/boundary_side_eff_length`, penalising only `f_r < ρ_m/ρ_tot`. This removes the
half-triangle bias and the hard clamp, and is **honest on the long-exon touchstone** (exon 2931,
`ρ_m=75.0 ≈ ρ_mat_uns=68.6`, ceiling `f_g ≤ 0.32`). It is still **refuted on the short-exon
majority**: on `gdna300`, 52/517 faces (25.9 % of spliced mass) assert `ρ_m` above the exon's true
unspliced-RNA density, 31 on ~pure-gDNA short exons (`mat_uns=0`, oracle `f_g=1.0`). The reason is
`S ⊥ f_g` *at the exon*: for `R ≪ FL` the mature coverage deposits *no* unspliced fragments in the
exon body (`P(unspliced|R)→0 ⇒ mat_uns=0`), so `S` cannot bound the exon's *unspliced* `f_g`. The
honest mode `log(ρ_m·P(unspliced|R)/ρ_tot)` needs a length factor = a magic number. **The corrected
`boundary_side_eff_length` is adopted — but for the absorption *frame* (§3.7) and as an
`S_eff` mode-scale only, never to license a composition mode.** Precision-only stands (§3.3(4)).

**The phantom is precision, not mode.** The boundary over-call (`f_g ≈ 0.62–0.97` at truth `0`)
is manufactured confidence, orthogonal to any mode-floor: at node 2930 the current solved `f_g`
(`0.77`, or `0.97` with the honest clamp off) sits **below** the combined floor (median `0.43`),
so even the combined floor does not catch it. It is cured by removing the dst-spliced absorption
subtraction (`_RNA_ABSORB → 0`) + the τ-gate — not by a floor.

**Structural guarantee — `S>0 ⟹ already off-vertex`.** A node sits at the confident all-gDNA
`{0,0,1}` vertex **only if** it is G1-locked (no free strand). `S>0` requires a motif-stranded
junction ⇒ that strand is transcript-continuous ⇒ `free_s = True` ⇒ the node is **solvable**
(G2/G3), initialized at the strand-solve reference midpoint `f_g ≈ 0.5`, **not** the vertex. So
the reviewer's "move it off the vertex" motivation has **no work**; the off-vertex state is
emergent from precision, never a floor.

**Resolution.** `I_motif` ships **precision-only** — `pr_motif = f_g²·S_eff` (`S_eff = spliced
mass Σw`), RNA-channel only, additive / can-never-cancel, `s2t`-bypassed, attached to the
**ordinary RNA imputation mode** — with **no** data-derived composition mode. The residual gDNA
fraction of the boundary's own unspliced pool is left to **strand** (`I_strand`, silent at
`κ=½`) + the **Phase-2 capture-safe intron hyperprior** (`ρ_g = ρ_bg·enrichment`), which grounds
node 2930 to `≈ 0.01` (oracle 0.007). **What `S` must NOT do:** pin/cap the unspliced `f_g`
(`S ⊥ f_g`); enter the gDNA channel (a shared `τ` manufactures `pr_g ≈ 8.7e4` from a pure-RNA
junction); assert an *equality* on the exon's RNA density (a two-sided pin re-injects
nascent-from-abundance); subtract from a neighbour's RNA density (the absorption bug); or infer
the nascent amount from the mature amount. **If a reviewer insists on an explicit statement,**
the only parameter-free honest form is the *combined-pool soundness check* `f_g,comb ≤ N/(N+S)`
— provably inert (`N_g ≤ N`), a soundness assertion **never** applied to the unspliced/exon axis
and **never** a force.

### 7.2 The mature→nascent damping is MOOT under the precision-only resolution

The damping was needed only to soften a `ρ_m → unspliced-f_g` **mode** mapping. With no
composition mode (§7.1), there is no mature→nascent mode to damp: `S` never enters the unspliced
`f_g`, so no transfer variance `s2t` is applied to it (`pr_motif` bypasses `s2t` because it is a
measurement of the boundary's *own* mature density, not an imputation across a crossing). The
only spliced-side sub-item that remains open is `S_eff` exactness (§7.3).

### 7.3 `S_eff` exactness

We ship `S_eff = mass Σw` (conservative). Plumbing `Σw²` for exact Kish
`n_eff = (Σw)²/Σw²` is one extra accumulator reduction; confirm it is worth the cost vs the
conservative lower bound. **Caveat (adversarial):** `S_eff = Σw` is the *statistical* Fisher
information of `S` about `log ρ_m` only; it does **not** absorb the *systematic* `E_spl`
short-flank geometry error. This is safe **because** `pr_motif` carries no composition mode
(§7.1) — the precision rides the ordinary imputation mode, so a biased `ρ_m` cannot manufacture
a confidently-wrong composition. Were a `ρ_m/ρ_tot` mode ever revived, `S_eff` would have to be
down-weighted by the `E_spl` geometry error; the precision-only design sidesteps this.

### 7.4 `I_struct` boolean vs a soft signature certainty

We treat a lock as composition-certain (`var_composition = 0`). If the G1/G2/G3 signature init
is itself uncertain, a soft certainty may be warranted — but that reintroduces a magnitude and
risks a magic number. The vertex-annihilation proof says the *emission* is insensitive to it
anyway; confirm no downstream consumer needs the soft form.

### 7.5 The `θ`-axis evidence-only principle

`τ^θ` must gate the per-strand RNA precision exactly as `τ^λ` gates the total — else the pooled
prior launders through `Var(θ)` (`0.35 → 7.69`). Confirm the two-axis mapping (3.5) fully closes
this.

### 7.6 Delta vs exact `Var(log f_g)`

≤7% mid-range, up to 2.5× near vertices. Confirm it does not perturb genuinely stranded
near-vertex nodes unintentionally.

### 7.7 Point-estimate vs variance harm

The defect is in the variance/precision channel; confirm the downstream fit's dependence is
genuinely through `1/var_g` (it is, per code) and that suppressing messages does not mask a
compensating point-estimate regression.

### 7.8 The onward-composition-precision `τ`-gate does NOT yet reach 0 (ship-blocker for the sink's silence)

The compiler's central promise — `τ^{λ,cav}=0 ⇒ pr=0` on a vacuous unstranded chain — is **not
delivered** by the shipped `_TAU_PRECISION` machinery, and this is the honest open item behind
§3.7. Reproduced on the touchstone (`gdna5`/`gdna_none`, `_TAU_PRECISION=True`):

- The strand balance is *fitted*, not exactly ½: `κ = 0.500975`, so `I_strand = N_eff·(2κ−1)²·[…]
  = O(1e-6·N) ≠ 0`. This clears the hard `τ > 1e-9` emit gate (`bp_solver.py:491`).
- The cavity accumulation `tau_lam[i] += math.fsum(pr)` (`bp_solver.py:585`) is a **positive-feedback
  loop**: a nonzero seed `τ` lowers `v_logf = f²/τ`, which *raises* the emitted `pr`, which *raises*
  the destination `τ`. From the `~1e-6` seed it self-inflates to `τ ≈ 0.84` at exon 2931 (~6 orders
  of magnitude), so `pr_g` is **byte-identical** with `τ` on vs off, and node 2930 stays at the
  `0.76` phantom.
- **The capture `s2t` term is the only thing holding `pr` down.** Strip capture (`s2t→0`) and the
  onward `pr_g → 4.8` — the full manufactured cascade the gate was meant to prevent.

**What Phase-B must prove before certifying the sink's onward silence:** (1) force `I_strand ≡ 0`
when `|2κ−1|` is within its own fitted sampling noise — a *strand-power* threshold from `SE(κ̂)`
(parameter-free), replacing the current `τ > 1e-9` value-gate — so the seed cannot clear the gate;
(2) show the cavity accumulation cannot self-bootstrap — exclude sub-threshold-seeded precision from
the running `τ`, or prove the feedback is contractive; (3) re-run the touchstone and demonstrate
onward `pr_g→0` **and** node 2930 `f_g→~0.01` **and** the capture-off regression guard holding at
`~0` (currently 4.8). Until then the sink ships its **mode-op + frame-fix + no-vote rule** (all
`τ`-independent) but **not** the τ-gated onward precision (§11.4–§11.5).

---

## 8. Proposed change to `bp_solver._scan` (formula-level)

**Current** (line ~446, both channels):

```
  v_logfg = (1 − f_g,s)² · var_lam[lsrc]          # belief variance  ← the defect
  v_logfr = f_g,s²       · var_lam[lsrc]
  pr_g = 1 / (v_logfg + s2t + 1/sm)
  pr_r = 1 / (v_logfr + s2t + 1/sm)
```

**Proposed.**

```
  # ---- seeded before the sweep (statics), per node i ----
  tau0_lam[i]    = I_strand(N_i, κ, f_g,i)                   # (3.1); 0 unstranded
  struct_lock[i] = signature-locked?  (G1 gonly / switched-off strand)   # boolean gate (3.2)
  S_eff[i]       = spliced MASS Σw over both faces (spliced_pos_* + spliced_neg_*)   # (3.3)
  tau0_th[i]     = I_strand^θ(...)                           # (3.1) θ-axis

  # ---- during the forward (resp. backward) sweep, source s → dest d ----
  tau_lam_cav = tau_running_lam[s]        # left-context-only; excludes the d-edge (cavity)
  tau_th_cav  = tau_running_th[s]
  if struct_lock[s]:  v_lam = 0.0                           # composition-certain: v_evid=0
  else:               v_lam = (1 − f_g,s)² / tau_lam_cav    # ∞ (⇒ pr 0) when tau_lam_cav = 0
  v_th  = (cos θ_s/(1+sin θ_s))² / tau_th_cav

  # gDNA channel — structural λ-face only, NO absorption:
  v_logfg = v_lam                                           # (exact VarLogFg optional)
  pr_g = 1 / (v_logfg + s2t + 1/sm)                         # = 0 when unpinned

  # RNA channel — ADDITIVE: imputation (zeroable) + motif (standing, s2t-bypassed):
  v_evid_r = f_g,s²·v_lam + v_th                            # two-axis (3.5)
  pr_imputation = 1 / (v_evid_r + s2t + 1/sm)               # honest-clamp path unchanged
  pr_motif      = (f_g,s² · S_eff[s]) if S_eff[s] > 0 else 0   # RNA-ONLY, bypasses s2t
  pr_r = pr_imputation + pr_motif                           # can-never-cancel

  # Phase-1b: residual ≤ 0  ⇒  pr_imputation = 0 AND no compensating gDNA vote (§6)

  # destination folds pr into BOTH belief (as now) AND tau:
  #   tau_running_lam[d] += pr_g            (pr_motif routed to the RNA/f_r τ, NOT gDNA)
  #   tau_running_th[d]  += pr_th
  # belief keeps its reference: belief_prec = tau + a, applied once at the marginal.
```

`s2t` and `1/sm` are unchanged for the imputation/structural channels. `pr_motif` is the only
term that bypasses `s2t` (a measurement, not a crossing). **`I_motif` carries NO composition
mode of its own (§7.1, resolved): `pr_motif` rides the ordinary RNA imputation mode. The
precision plumbing above is the whole of `I_motif` and is safe to land; there is deliberately no
`ρ_m/ρ_tot` (or `N/(N+S)`) floor on `f_g`.**

This section is the *source-direction* skeleton; §11 is the full **Phase-B plan** — it adds the
per-strand one-directional routing that makes `pr_motif` an **emission to the exon** (§3.3), the
**reverse/sink** direction (absorption frame-fix + no-vote rule, §3.7), and the `τ`-gate
prerequisite (§7.8) that gates the onward composition precision. Land B1 (emission seed + routing +
sink frame-fix + no-vote) first; B2 (onward `τ`-precision) is blocked on §11.5.

---

## 9. Acceptance tests (corrected — the reviewer's "detect in pass 1" is reframed)

- **T1 — PHANTOM COLLAPSE (zero/low gDNA).** `gdna_none` & `gdna5`: every live node's final
  `var_g` within ~10% of the message-free local (~2.8); median live sharpening → ~1.0× (from
  3.54×); node 2931 `pr → 0`, `f_g → ~0.49`, `var → ~2.82`.
- **T2 — CLEAN PRECISION SUBSTRATE.** `gdna300`: the 1165 gDNA-rich open exons (59% mass, all
  ≥2 hops from any lock) end at `τ ≈ 0 ⇒ var ≈ 2.8`, **not** the current fake `var 1.61`; the
  downstream `1/var_g` weight is no longer inflated (hyperprior L1-to-oracle drops from ≈2.0).
- **T3 — REFRAMED: NO DEGRADATION on `gdna300` (delete "detection").** Gated final `f_g` on the
  enriched exons ≈ local `0.509` (`|err| → 0.344`, **improved** from 0.361) with honest
  `var ≈ 2.8`, not a confident-wrong `0.495 / var 1.61`. *(The reviewer's original T3 —
  "detect real gDNA in the first pass" — is unachievable; see T8.)*
- **T4 — MOTIF ANCHOR PROPAGATES AS PRECISION, NEVER A COMPOSITION PIN (load-bearing
  positive).** Both libraries: the 764/749 spliced-carrying boundaries still emit `pr_r > 0`
  (node 2930 `pr_motif ≈ 7502`, not 0), raising adjacent exon-flank **RNA precision**; no drift
  to gDNA at heavy-spliced junctions. **Touchstone node 2930** (`N=962, S=7502, floor
  N/(N+S)=0.114`): its own unspliced `f_g` must land **neither** at the confident FP `~0.97`
  (the absorption-clamp phantom — killed by removing the dst-spliced subtraction + the τ-gate)
  **nor** hard-pinned to `0` (no `N/(N+S)` or `ρ_m/ρ_tot` floor touches it). Pass 1 leaves it
  **unpinned** at `f_g ≈ 0.5`, `var ≈ 2.8`; the Phase-2 capture-safe intron hyperprior then
  grounds it to `≈ 0.01` (oracle 0.007). *(This is what strand-only gating would wrongly kill,
  and what a hard mode-floor would wrongly pin.)*
- **T5 — STRUCTURAL LOCK PROPAGATES 1 HOP where substrate exists.** A locked intergenic node
  emits `pr = 1/(s2t + 1/M_src) > 0` to a directly-adjacent open node. **Explicit note:**
  `gdna300` has 0 lock-adjacent gDNA-rich exons, so run this on the intron/intergenic
  neighbourhood or a constructed lock-adjacency toy — **never** as a `gdna300` enriched-exon
  test.
- **T6 — STRANDED PROPAGATION intact.** `κ ≠ ½` library: a genuine strand anchor pins its node
  at `~1/I_strand`, attenuates `~+s2t` per hop; no regression.
- **T7 — BOTH LATENTS.** Vacuous chain: the `θ` (per-strand RNA) precision no longer climbs
  `0.35 → 7.69`; `τ^θ` gating collapses it.
- **T8 — HYPERPRIOR RECOVERS ENRICHED gDNA (Phase-2, the deferred half of old-T3).** Fed the
  clean `τ` substrate from T2, the stratified `ρ_bg × enrichment` hyperprior moves the
  `gdna300` 59%-mass exons from ~0.5 toward oracle 0.85 while `gdna_none` stays ~0. This is
  where legitimate detection lives.
- **T9 — NON-NEGATIVITY / DETERMINISM.** `τ ≥ 0`, `pr ≥ 0`, no clamp firing; cross-process
  determinism unchanged.
- **T10 — NO SPLICED COMPOSITION MODE-FLOOR (the resolved §7.1).** No first-pass factor caps
  the unspliced `f_g` (nor an exon flank's `f_g`) from a spliced count: no `f_g ≤ N/(N+S)`, no
  `f_g ≤ 1 − ρ_m/ρ_tot`. Regression guard on `ambig_dense_10mb`: **zero** hard-clamps to
  `f_g = 0` from a negative barrier (the 32 `gdna5` / 50 `gdna_none` short-flank faces the
  `ρ_m/ρ_tot` ceiling would clamp) and **zero** gDNA-island exons ceiled below oracle (the 193
  `gdna300` faces / boundary 640 → exon 639 case). The only permitted explicit
  spliced-composition statement is the provably-inert combined-pool soundness check
  `f_g,comb ≤ N/(N+S)`, never applied to the unspliced/exon axis and never as a force.

---

## 10. The boundary emission ↔ sink, end to end (touchstone 2930→2931, reproduced)

The neighbourhood `2929(INTRON⁺) — 2930(BND, splice-acceptor, js=+) — 2931(EXON⁺)`,
`ambig_dense_10mb / gdna5`, pass-0:

| node | oracle unspliced (gDNA / mat_uns / nas) | spliced | oracle `f_g` | local `f_g` | current solved `f_g` |
|---|---|---|---|---|---|
| 2929 INTRON⁺ | 0 / 0 / 78 | 0 | 0.000 | 0.510 | 0.510 |
| 2930 BND | 7 / 0 / 955 | **7502** (right/exon face) | **0.0073** | 0.510 | **0.771** (phantom) |
| 2931 EXON⁺ | 493 / 73639 / 44889 | 0 | 0.0041 | 0.490 | 0.625 |

- **Emission 2930→2931 (§3.3):** `S_eff = 7502` (mass `Σw`), RNA-`pos` channel only,
  `s2t`-bypassed, additive/can-never-cancel. It **raises exon 2931's `f_pos` precision** with
  **no mode** — the exon's own 73,639 `mat_uns` already carry the composition point. The
  density-mode lower bound is *honest here* (long exon) but **refuted on short exons** and does not
  ship (§3.3(4) / §7.1).
- **Boundary's own `f_g` (§3.3):** untouched by its own `S` (`S ⊥ f_g`). Target: strand
  (`I_strand≈0`) + λ-imputation leave it unpinned `≈0.5`; Phase-2 grounds it to `≈0.01` (oracle
  0.007). The current `0.771` is the sink-laundering phantom (§3.7) — **not** yet fixed on `main`.
- **Sink 2931→2930→2929 (§3.7):** the exon's incoming RNA density (55.6) is *below* the junction's
  own measured mature (75.0) ⇒ `residual < 0` ⇒ no nascent relays ⇒ (rule) no RNA message **and**
  no gDNA vote. Onward composition-precision **target** 0.
- **Honest verdict:** the emission precision-only + routing (§3.3, Derivation C) and the sink
  mode-op + frame-fix + no-vote rule (§3.7) are **ship-ready**; the onward `τ`-gate that would drive
  `pr_g→0` and collapse the `0.771` phantom is **NOT** (§7.8: the `τ` machinery self-bootstraps;
  capture `s2t` masks a `pr_g→4.8` leak).

---

## 11. Phase-B implementation plan (formula → code)

The plan lands in **two tranches**: **B1 (ship-ready, `τ`-independent)** — the emission precision
seed + routing, the sink absorption frame-fix, and the no-vote rule; **B2 (blocked on §7.8)** — the
τ-gated onward composition precision. Do **not** land B2 before §11.5.

### 11.1 Effective-length / accumulator

- `effective_length.py`: **no new function.** `boundary_side_eff_length` (`E[min(ℓ,R)]/2`, line 131)
  is the short-flank-safe divisor; `spliced_side_eff_length` (the half-triangle, line 165) stays for
  the emission's genuinely *one-sided crossing* geometry but is **removed as the ABSORPTION divisor**
  (§11.4).
- Accumulator: **none** to ship the conservative `S_eff = mass Σw` (already plumbed as `spliced_*`).
  *Optional* exact Kish needs a new reduction `Σw²` (`spliced_w2_{pos,neg}_{left,right}`) through the
  same `_spliced_faces` gate — defer unless the conservative bound proves too loose (§7.3).

### 11.2 `node_geometry`

- Routing: **no change** — `_spliced_faces` (`node_geometry.py:170`) already routes per-strand,
  one-sided, motif-gated (Derivation C, verified).
- Expose the boundary's exon-flank `boundary_side_eff_length` (`b_eff_r_*`, already computed at
  `node_geometry.py:147`) as the SINK's absorption divisor, replacing the half-triangle
  `b_eff_spl_*` (line 149) at the `ESPd` consumption site.

### 11.3 `bp_solver._scan` — EMISSION (source boundary → dest exon) — **B1**

After the existing RNA-pos `pr` (`bp_solver.py:547`; symmetric for neg at 559):

```
S_eff_p    = SP[sf][lsrc]                                   # source-face spliced-pos MASS; 0 off-strand / non-junction
pr_motif_p = (fg_s*fg_s) * S_eff_p  if S_eff_p > _EPS else 0.0   # RNA-pos ONLY; bypasses s2t & 1/sm
app[i]     = pr + pr_motif_p                                # ADDITIVE: imputation (zeroable) + standing motif
```

- Route `pr_motif_p` to `app` / `apn` (RNA) **only** — never `apg` (gDNA); a shared `τ^λ` would
  manufacture `pr_g = f_g²/f_r²·S_eff = 86,965` (verified phantom).
- Gate on `S_eff_p > 0` **alone**, not `emit_p` — `pr_motif` is a *standing* factor; coupling it to
  the imputation gate risks zeroing the anchor when the imputation path legitimately closes
  (Derivation C §C.7).
- Mode: **unchanged** (the existing imputation `mo`, `bp_solver.py:541`). No `ρ_m/ρ_tot`.

### 11.4 `bp_solver._scan` — SINK (dest boundary absorbs; onward relay)

- **Frame fix — B1 (ship now):** change the absorption divisor from `ESPd[i]` (half-triangle) to the
  boundary's `boundary_side_eff_length` at `bp_solver.py:538` (pos) / `:554` (neg), so
  `absorb_p = SPd[i]/b_eff_r_dst[i]` — like-for-like with the nascent density it subtracts from.
- **No-compensating-vote — B1 (ship now):** when `rho_pos ≤ 1/erd` (residual absorbed), emit no RNA
  message **and** add no gDNA precision from the RNA absence (`τ += 0` on both). `τ`-independent.
- **Onward composition precision — B2 (BLOCKED on §11.5):** replace `v_logfp = f_g²·σ²_λ` with
  `f_g²/τ^{λ,cav}` and the belief-variance `n_eff` precision with
  `pr_onward = 1/(f²/τ^{cav} + s2t + 1/M_src)`. Keep behind a toggle until §11.5 passes.

### 11.5 The `τ`-gate prerequisite (§7.8) — must precede §11.4's B2

1. Replace the `tau_lam[lsrc] > _EPS` value-gate (`bp_solver.py:491`) with a **strand-power** gate:
   `I_strand ≡ 0` when `|2κ−1| ≤ z·SE(κ̂)` (the fitted-`κ` sampling noise from the strand-balance
   fit — parameter-free, not a magic threshold).
2. Prove/enforce that the cavity `tau_lam[i] += Σpr` (`bp_solver.py:585`) is non-self-bootstrapping:
   exclude precision seeded below the strand-power threshold from the running `τ`, or demonstrate
   contractivity.
3. Acceptance: onward `pr_g→0` on `gdna5`/`gdna_none`; node 2930 `f_g→~0.01`; the **capture-off**
   regression guard (`s2t=0`) holding at `~0` (currently 4.8).

### 11.6 Acceptance tests

- **PB-1 — emission anchor.** Node 2930 emits `pr_motif_p ≈ 7502` into exon 2931's `app`; exon RNA
  precision rises; `prec_g` at the exon is **byte-identical** with `pr_motif` on/off (no gDNA leak —
  guards the `86,965` phantom).
- **PB-2 — boundary own `f_g`.** Node 2930's own `μ_λ, σ²_λ` byte-identical with `pr_motif` on/off;
  pass-1 leaves it `f_g≈0.5, var≈2.8`; Phase-2 grounds it to `≈0.01` (oracle 0.007).
- **PB-3 — routing.** `pr_motif` lands on `app` iff `js=+`, `apn` iff `js=−`; intron-flank neighbour
  gets 0; each face fires in one sweep direction (no 2× in `_comb`). (Node 274 `js=−`: `apn` only.)
- **PB-4 — `S=0` degrade.** Every `S=0` face emits `pr_motif=0`; determinism unchanged.
- **PB-5 — sink frame.** Switching the absorption divisor to `boundary_side_eff_length` removes the
  spurious negative residuals the 2–199× half-triangle manufactures on the 177 short faces (zero
  short-flank faces flipped to a spurious "no nascent + gDNA vote").
- **PB-6 — no-vote.** With `rho_pos ≤ 0`, the gDNA precision added at the boundary is 0; node 2930
  does not launder to `f_g→1`.
- **PB-7 — τ-gate (blocks B2).** Onward `pr_g→0` on the vacuous chain **and** the capture-off guard
  at `~0` (regression against the reproduced 4.8). Until PB-7 passes, §11.4's onward precision stays
  behind a toggle.

---

## Appendix: reference numbers (all reproduced this session)

```
  Var_ref(λ)          = 8.6743     ref precision a          = 0.1153
  Var_ref(log f_g)    = 2.8037     ref precision on log f_g = 0.3567
  E_ref(f_g)          = 0.5000
  I_strand(κ=½)       = 0.00000    (every N, every f_g);  κ suite = 0.4992/0.5010/0.50117

  --- v1 worked example (gdna_none, oracle f_g=0) ---
  node 2931 current pr = 4.315 = 1/(0.0836 + 0.1481 + 0.00009);  local→final var 2.818→0.080 (35.1×)
  GMRF interior 401-node vacuous: Var(λ) 8.674→0.569 (15.2×) at σ_e²=0.15;  a→0 ⇒ improper (all prior)

  --- v2 vertex annihilation (gdna300, I_struct=N_eff refuted) ---
  lock-sourced edges: 611 fwd / 635 bwd, med f_g,src=1.0000, med v_logfg=0.000000
  pr@vls=ref(8.674) == pr@vls=0 == pr_actual = 2.084 (fwd) / 2.067 (bwd)   ← magnitude INERT
  med s2t=0.123, med 1/M_src=0.203/0.217 (M_src≈5)   ⇒ boolean gate, not N_eff

  --- v2 I_struct reach (gdna300) ---
  216 gonly locks, 4.1% mass, 0 one-hop gDNA-rich adjacency; 206/216 reach at 2/2/4 hops

  --- v2 I_motif (node 2930) ---
  own spliced: mass 7645.7/7501.8, count 16031/15753 (gdna300/gdna5) ≈ 2×mass
  absorb ON:  rho_pos −7.5/−19.4 ⇒ n_eff=0 ⇒ RNA pr COLLAPSES 148812/119007 → 0
  I_motif RNA pr = S_eff(mass) = 7645.7/7501.8;  phantom gDNA pr (shared τ, REFUTED) = 86965/85329
  substrate: 764/749 boundaries, 767718/763520 reads, 322761/321485 mass (library-invariant)

  --- v2 identifiability ceiling (gdna300, true live gDNA 0.691) ---
  gDNA-rich open exons fo>0.5: n=1165 (59.3% mass), oracle 0.853, local_fg 0.509 → final 0.495
    |loc−or| 0.344 → |fin−or| 0.361 (messages move AWAY);  var 2.811 → 1.611 (fake collapse)
    ⇒ gating loses nothing; enriched-island gDNA is Phase-2's job
```

Core code: `src/rigel/calibration/bp_solver.py` (`node_sweep`, `_scan` ~L405–503; L418
`vls = var_lam[lsrc]`, L446 the `pr` formula), `node_geometry.py` (signature locks;
`spliced_*` mass / `spliced_n_*` count faces; `init_beliefs` G1 `f_g=1, var_g=0`),
`strand_likelihood.py` (`I_strand` moments), `simplex_logodds.py:356` (the `var_g` grid
moment). Invariants: `docs/calibration/CALIBRATION_ARCHITECTURE.md` (count-zero-information),
`docs/calibration/CALIBRATION_MASTER.md` (the 3 sources, the two NPMLE roles, the stratified
hyperprior). Verification tools (scratchpad): `vertex_probe.py`, `motif_verify.py`,
`struct_substrate.py`, `local_vs_final.py`, `chain_accum.py`, `prec_audit.py`,
`delta_check.py`.
