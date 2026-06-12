# Phase C (merged B+C) — AMBIG residual resolution: dry-run plan

**Status:** dry-run plan, 2026-06-12. Third phase of [`deconvolution_roadmap.md`](deconvolution_roadmap.md),
**merging the former Phase B** (mature subtraction) into C per the Step-0 finding
([`phaseB_mature_subtraction_plan.md`](phaseB_mature_subtraction_plan.md) §5.5: the subtraction only earns
its keep on AMBIG nodes — a strand-observable node's strand already separates gDNA from *all* RNA). This is
where the flagship leak drops. Builds on Phase A (`mature_density.py`, shipped `83b3610`).

## 1. The problem, exactly

An **AMBIG** node (overlapping `+`/`−` transcripts) carries contained unspliced mass that is a mixture of
**three** components we must allocate — `RNA+`, `RNA−`, `gDNA` — and its own strand **cannot** separate them
(both strands present). With observed `U_pos, U_neg` (κ = `rna_sense_frac`) and the three components
(gDNA symmetric ½/½; `RNA_s` = nascent+mature on strand `s`, reading its sense at κ):

```
U_pos = ½·gDNA + κ·RNA₊ + (1−κ)·RNA₋
U_neg = ½·gDNA + (1−κ)·RNA₊ + κ·RNA₋
```

Two equations, three unknowns → **under-determined by the node's own counts**. Resolution needs outside
information, drawn from the bounding boundaries (spliced → mature, unspliced → gDNA) and the strand tilt.

## 2. Three estimators, then blend (your framework, mapped)

**(1) Strand deconvolution — independent, naive (the "tilt").** Run the strand posterior on the node's own
`U_pos/U_neg` with no prior on which strand is expressed. It yields the *symmetric vs tilted* split: the
difference pins the **net RNA tilt** `RNA₊ − RNA₋ = (U_pos − U_neg)/(2κ−1)`; the sum is `gDNA + RNA₊ + RNA₋`.
Symmetric ⇒ "looks like gDNA *or* balanced two-strand RNA" (the irreducible ambiguity); tilted ⇒ net RNA on
one strand. Information, not a solution.

**(2) Count imputation — independent, no strand.** From the bounding boundaries:
- **mature** per strand, `M₊/M₋`, from the **spliced** crossings (Phase A `mature_density`) — accurate
  (Step 0: 0.96–1.00 on the leak regions);
- **gDNA**, `g_count`, from the **unspliced** crossings (`density_model`) — capture-depleted, a **lower
  bound** on the true gDNA.

These need **not** sum to `U` (nascent is unaccounted, gDNA is depleted) — so they enter as
**fractions/anchors** that help allocate `U`, not as absolute truth.

**(3) Blend.** Combine (1) and (2). We reuse the established precision-weighted combine `g = w·g_strand +
(1−w)·g_count`, with `w` extended to encode *how identifiable the gDNA is at this AMBIG node* (§4). The
richer MLE/iterative form (strand ⇄ count) is **deferred** (§7) — one feed-forward pass first.

## 3. The cascade for an AMBIG node (one iteration)

1. **Subtract mature, strand-aware** (former Phase B math, AMBIG-only): `U'_pos = max(U_pos − κ·M₊ −
   (1−κ)·M₋, 0)`, mirror for neg (clamp ≥0 handles zero-mass slivers). Residual `U' = gDNA + nascent₊ +
   nascent₋`. Mature `M = M₊+M₋` joins the deterministic RNA (like spliced).
2. **Orient the residual** by the **dominant expressed strand** (`M₊` vs `M₋`): sense = that strand's sense.
   (The per-strand mature is the directional prior — generative §6 "supply the missing equation".)
3. **Strand on the residual** → `g'_strand` = gDNA fraction of `U'` (symmetric ⇒ gDNA; tilt ⇒ the tilted
   strand's nascent is RNA). Reuse `strand_posterior_gdna_frac` on the oriented residual; convert to a
   fraction of `U` by `×U'/U` (the Phase-B mass bookkeeping, option 2a).
4. **Blend** with `g_count` (the depleted gDNA lower bound) by `w` (§4).
5. `gdna = g·U`; `rna = (1−g)·U + spliced` (mature `M` is inside `(1−g)·U`). Mass conserved.

## 4. The strand-trust weight `w` at an AMBIG node — the key design

The strand cleanly separates gDNA from RNA **only when the RNA sits on one strand** (the other strand
RNA-free). The error in "gDNA = symmetric part" is exactly the **balanced (two-strand) nascent**, whose size
tracks the **weaker strand's** RNA. The mature is our proxy for nascent (same gene), so the confusion scales
with **`min(M₊, M₋)`**:

- one strand silent (`min ≈ 0`): the symmetric residual **is** gDNA → trust the strand (`w → 1`). This is
  the **silent-strand override** — emergent from the weighting, not a hard gate.
- both strands expressed (`min` large, balanced): symmetric ⇄ balanced-nascent ambiguity → **don't** trust
  the strand (`w → 0`) → defer to the count (depleted) gDNA + the documented floor.
- **no mature either strand** (`M₊+M₋ ≈ 0`): no nascent ⇒ symmetric is gDNA ⇒ **trust** the strand
  (`w → 1`). (Note: this is *RNA-free-ness of the weaker strand*, **not** the imbalance ratio
  `|M₊−M₋|/(M₊+M₋)`, which mis-handles the no-mature case — a correction from the Step-0 analysis.)

So extend the existing `w = I/(I+I₀)` with an AMBIG **identifiability factor**: `I_ambig = N'·(2κ−1)²·ρ²`,
where `N'` = residual count, and `ρ ∈ [0,1]` is the weaker-strand RNA-free-ness — first cut `ρ = 1 −
min(M₊,M₋)/max(M₊,M₋)` (with `ρ = 1` when `M₊+M₋ ≈ 0`). For a single-strand node `ρ = 1` and this reduces to
today's `I = N·(2κ−1)²` — **the same formula**, AMBIG just carries the extra factor. The exact form of `ρ`
(and whether to scale `min(M)` by the residual or use a likelihood-curvature derivation) is to be
**calibrated against the toy scenarios in execution** — it is the one place a shape choice enters; keep it
derived, not tuned (no magic number — pause & discuss if one appears).

**Why this passes the AMBIG regression tests:**
- `test_ambig_no_false_gdna_from_nascent` (gDNA=0, nascent, **both** transcripts expressed): both strands
  have mature ⇒ `ρ→0` ⇒ `w→0` ⇒ `g_count` governs; gDNA=0 ⇒ boundary gDNA density ≈ 0 ⇒ `g≈0`. ✓
- `test_ambig_reads_gdna_when_present` (gDNA, no nascent, both transcripts): both have mature ⇒ `ρ→0` ⇒
  `g_count` governs ⇒ reads the (present) boundary gDNA density. ✓
- **flagship region 226** (− mature dominant, ~0 `+` mature): `ρ→1` ⇒ `w→1` ⇒ `g'_strand` governs; the
  no-nascent residual is symmetric ⇒ `g' ≈ U'/U = 0.54` (oracle 0.59). ✓ — the leak drops.

## 5. Retiring the splice-junction fraction

`region_splice_gdna_frac` is **removed** (the call in `calibrate.py`; keep `splice_junction_eligibility`,
used by `mature_density`). Its job — debias the AMBIG exon gDNA fraction — is now done by the mature
subtraction + residual strand. Safe because: high-imbalance AMBIG (where the splice fraction tried hardest
and failed, e.g. 226 at 0.37) is now governed by `g'_strand` (0.54); low-imbalance AMBIG falls to `g_count`
(the same depleted density the splice fraction started from, minus its flawed upgrade) — acceptable, and the
`test_ambig` cases confirm no regression.

## 6. Implementation (the changes)

- `calibrate.py`: compute `md = mature_density(...)`; pass `M₊/M₋` into the **region** deconvolution; remove
  the `region_splice_gdna_frac` call and its `region_node_density` replace.
- `strand_deconv.deconv_regions` / `_deconv_per_node`: for **AMBIG** regions, (a) strand-aware subtract
  mature → residual sense/antisense (oriented by the dominant strand), (b) make AMBIG `use_strand=True` on
  the residual, (c) `w` carries the `ρ²` identifiability factor (regions only), (d) rescale `g'_strand` by
  `U'/U`. **Single-strand (POS/NEG) and the boundary sides are unchanged** (Step-0: no benefit there; raw
  strand already separates gDNA from all RNA). `deconv_sides` untouched (sides are mature-free).
- Keep the change surgical: AMBIG-only new path; the single-strand path is byte-for-byte the same.

## 7. Edge cases & deferred

- **No-junction AMBIG regions** (Step 0: 4 regions, 26k mature, `M=0` — single-exon overlapping genes): no
  mature anchor ⇒ no orientation. `ρ=1` (no mature ⇒ presumed no nascent) ⇒ trust the strand on the **raw**
  AMBIG counts (symmetric ⇒ gDNA). If that is unsafe (genuine balanced nascent with no mature, e.g.
  single-exon nascent), fall back to `g_count` / **neighbour carry-over** (the nearest informative
  same-class boundary). Decide in execution against the toy scenario; log where it triggers.
- **Both-strands-expressed with real gDNA** (balanced nascent *and* gDNA): `ρ→0` ⇒ `g_count` (depleted) ⇒
  gDNA under-called — the documented one-pass floor. The neighbour carry-over / a richer prior is the lever;
  revisit only if the suite shows it dominates the residual leak.
- **Iteration / MLE (deferred, your step 3 tail).** The gDNA estimate refines the nascent/mature split which
  refines the tilt which refines gDNA — a fixed point. The one-pass cascade is the first iterate. Defer the
  loop until the one-pass is validated; if added later, profile-ML the per-strand nascent/mature and ML the
  local gDNA over the node + its boundaries (generative §6) — the curvature *is* the precision, replacing
  the hand-chosen `ρ`.

## 8. Validation

1. **AMBIG regression tests** (`tests/calibration/test_ambig_scenario.py`) — both must still pass (§4).
2. **Toy AMBIG scenario** — sweep gDNA 0→high × nascent {none, both, one-strand} and confirm: silent-strand
   → strand governs (recovers gDNA); both-expressed → count floor; the `ρ` weighting behaves monotonically.
   Calibrate the `ρ` form here.
3. **Flagship locus 21** — AMBIG regions 226/231/236 recover to ~0.54–0.57 (the residual fraction Phase A
   validated); single-strand exons unchanged; net leak drops.
4. **Full gDNA suite** (20 conditions) — net `gdna→rna` leak falls on the capture/stranded/complex
   conditions; **no regression** on the simple / unstranded / no-gDNA conditions (the `ρ` and the
   subtraction must be no-ops there). Golden refresh + review.
5. New debug script `scripts/debug/phaseC_ambig_resolution.py` (committed).

## 9. Risks

- **The `ρ` form** (§4) — the one shape choice; calibrate on the toy sweep, keep derived. Over-trusting the
  strand at balanced expression re-introduces the `test_ambig` failure; under-trusting forfeits the leak fix.
- **Orientation at near-balanced expression** — the dominant-strand pick flickers when `M₊≈M₋`; but there
  `ρ≈0` so the (mis)oriented strand is down-weighted anyway. Confirm the transition is smooth.
- **No-junction AMBIG** (§7) — the fallback path; validate it doesn't manufacture gDNA where balanced
  single-exon nascent exists.
- **Golden + suite drift** — AMBIG outputs shift; regenerate and review; confirm single-strand/simple
  conditions are byte-stable (the surgical-change guarantee).
