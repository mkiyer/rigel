# NEXT SESSION — the second feeder of the zero-control gap: an EMPTY slot's own zero-count claim into the mass pin

    ⚠ **A DEV DOC. Nothing may cite it, and it is NOT the state.** The state is `ROADMAP.md` §0/§1, the
    rulings `DESIGN.md`, the lessons `TRAPS.md`. ⛔ When this plan is executed, DELETE this file in the
    same commit.

    Written 2026-08-18 (session 2), at the end of the session that walked `g00 ss0.50 nrna_mid capture_on`
    (62 k → 228 k under the per-strand licence) to its mechanism, fixed one of its two feeders, and priced
    the other feeder's two candidate repairs on all 32 conditions. The lesson is
    `TRAPS: zero-the-precision-with-the-value`; the state is `ROADMAP.md` §1 rank 1.

## WHERE THE TREE IS

* `message_propagation = True`; the per-strand licence (`rna_level_scale`) is landed; and as of this
  session **the licence's zeroing runs AFTER the SPLICE IN block in both twins**, so a refused RNA arm delivers
  neither value nor precision (`test_gdna_scale_rule.py::test_the_splice_in_delivers_no_precision_on_a_refused_arm_*`,
  each fired by breaking its own twin).
* Panel (all 32, certified truth, Σ|f_g − true|·mass over ALL live slots): in-scope contaminated
  **0.931 / 0.999 / 1.000** vs the licence tree; `g00` **430,871 → 324,377**, all 8 rows better;
  deferred 1.000 / 0.827. Suite 3,523 / 9 xfail / 0 failed / 3,532 collected.
* ⛔ **THE BAR STILL OPEN: `g00` 324 k against the pre-licence relay's 154 k; 193 k of it is the one
  condition above, and all of that is FEEDER ②.**

## THE MECHANISM (so it is not re-derived)

One shape: **a zero RNA claim with LIVE PRECISION, relayed one hop, into the mass pin.** The pin's licence
`_lend = pop & pg[s] > 0 & (pp[s]+pn[s]) > 0` reads the arm as SUPPLIED, restores `Σρ_c E_c = M` with a
common factor `k = M/(tg·E_g + 0·E_r)`, and delivers `tg = M/E_g` — "all your mass is gDNA" — at the gDNA
arm's unchanged precision. Walk it yourself: `slot 3565 (TES−) → 3566 (EMPTY exon, M = 0, E_r = 0.018) →
3567 (ACC−, M = 2,160) → 3568 (exon, M = 26,296, true 0.000 → 0.9037)`.

* Feeder ① (FIXED): the SPLICE IN added the sj COUNT's precision back onto a zeroed arm (`n = 0 @ pn = 280`).
* Feeder ② (OPEN): `struct_lock = ~solvable ∧ REGION` (the strict xfail in `test_region_init.py`) makes an
  EMPTY exon composition-CERTAIN, so `own_precision` gives it `1/trigamma(½) = 0.2026` on EVERY arm at
  density 0. Before the licence the neighbours' relayed RNA level overwrote that zero (12.6 ≫ 0.2); at a
  terminus the licence now isolates it, and it becomes the source's whole RNA claim.

## THE TWO CANDIDATES PRICED, AND WHY NEITHER LANDED (all 32, per stratum, arm ÷ shipped-with-①)

| arm | in-scope contaminated (unstr×OFF / str×OFF / str×ON) | `g00` (8) | deferred | the loss, NAMED |
|---|---|---|---|---|
| `g1` — `struct_lock = g1_locked ∧ REGION`, own precisions recomputed | 1.025 / 1.021 / 1.005 | **0.71×** (232 k) | 1.17 | the four `nrna_none` zero controls go **1.8–15×** worse: at AMBIG `exon\|exon` boundaries in an RNA-only library ψ was held at `f_g ≈ 0` by the empty exons' "gDNA = 0 @ 0.2026" (a zero count over ~0.05 placements — right at `g00` for ANY reason). Without it, an RNA+ level claim (0.56 @ 72) with NO RNA− claim drifts ψ to **0.38** at slot 44476 (`g00 ss0.50 nrna_none off`). ψ-side question |
| `pinM` — `_lend` additionally requires `M[src] > 0` (a source with no mass has no composition) | 0.999 / 1.004 / **1.032** | **0.41×** (134 k — 0.31× the licence tree, better than the pre-licence relay's 154 k) | 1.137 | capture-ON only, six of six stranded × ON rows worse (1.02–1.13): under capture the empty slots between probe-covered stretches are the CONDUITS a RELAYED composition travels through (`TRAPS: the-divergence-was-a-barrier`) and the pin at those hops is load-bearing |

⭐ **The sharp predicate is neither.** Refuse the pin for a source's OWN zero-count RNA artefact (feeder ②)
and KEEP it for a relayed composition passing through an empty slot. Structurally that is what `g1` does
(an empty slot then has no own precision, so whatever it carries is relayed) — so `g1`'s remaining cost is
the thing to open next, and it is on the ψ side, not in the relay.

## THE NEXT DISSECTION, STEP BY STEP

1. Reproduce: run the shipped tree and the `g1` arm on `gdna_g00_ss_0.50_nrna_none_capture_off` (arm plumbing
   below), rank slots by err(g1) − err(shipped) — the top slots are AMBIG `B exon|exon` (44476, 44504, 18449 …).
2. The hop walk at 44476 (backward, 6 hops) shows the delivered packet is IDENTICAL between the arms except
   `cm_g` (0.0055 vs 0.136). So: **why does ψ book the remainder as gDNA rather than RNA− when it receives
   RNA+ = 0.56 @ 72, RNA− @ 0, and no gDNA claim?** Open `simplex_logodds` at that slot with the captured
   `global_lp` / `solve_grid` (the `psi_channel_ablation.py` replay is the tool). Suspects: the tilt channel
   `theta_prec = cm_p + cm_n` on an AMBIG boundary at κ = ½ (rank 10's f32 ridge is nearby), or the RNA
   channel's mode/precision mismatch when one strand's precision is 0.
3. If ψ is right and the boundary genuinely needs a gDNA claim: the honest source is the empty exon's
   zero-count gDNA MEASUREMENT — `a-zero-count-is-a-measurement` — which `g1` silences because
   `own_precision` needs a finite COMPOSITION variance. Consider: at count = 0 every component's density is
   bounded by the total's, so the composition variance is irrelevant to a zero — a zero-count slot could keep
   `1/count_logvar(0)` on the gDNA arm ONLY IF the same argument does not hand it back to the RNA arms
   (which is feeder ②). Derive it; do not tune it (`no-magic-numbers`).
4. Price the PAIR on all 32 with a bar WRITTEN BEFORE the panel returns (this session's: accept iff `g00`
   total improves with no single `g00` row > 5 % worse AND every in-scope stratum ≤ 1.02× — reject a `g00`
   win bought on the contaminated in-scope strata), never a half.

## ARM PLUMBING — rebuild it from the repo in minutes; the scratchpad is SESSION-SCOPED and will be gone

Everything used this session is a thin composition of three things already in the tree:
* the COUNTING arm — `relay_pool_ab.py::arm` (assert `_uni` present, assert the policy was built): patch
  `rigel.calibration.calibrate.HeadPolicy` (the MODULE attribute, reached via `importlib.import_module` —
  `rigel.calibration.calibrate` the NAME is the FUNCTION) with a factory that counts constructions and returns
  `HeadPolicy(HeadSwitches(...))`; `pre` = `rna_level_scale=False`, `nopin` = `mass_pin=False`;
* the `g1` arm — `ladder_arm_ab.py::_install_zc_struct_lock_g1`, verbatim (wraps `build_region_init` in BOTH
  `region_init` and `sweep`, rescopes `struct_lock`, recomputes the three own precisions, counts firings);
* the SCORER — `calibration_walk.py`'s: `slot_truth.npz` (certified) → `Σ|f_g − true|·mass` over `live`, per
  `stratum`; per condition, per arm, one JSON row with the tree stamp (`git rev-parse` + dirty flag), sharded
  by condition across processes; report per stratum (`ss_0.50`/`ss_0.99` × `capture_off`/`on`, contaminated
  rungs vs `g00`) and per condition, arm ÷ base. `pinM` was a TEMPORARY src edit (`and M_l[s] > 0.0` on
  `_lend` in `scan`, `& (M[src] > 0.0)` on `lend` in `_transport`), never landed.
The hop walk reads `_uni_static`'s `fwd_*`/`bwd_*` (the raw running state, published by the backbone),
`_pin[-2:]` (the final sweep's two transports) and `_uni[-1]` (the delivered packet), with the per-strand
licences recomputed from `boundary_flags` exactly as `head.py`'s `__init__` does.

## THE QUEUE BEHIND IT (unchanged)

* `transfer_var` off read 0.0010 at the first crime scene — understand before touching.
* the 2.5 % boundary-crossing deficit (owner: priority, not highest).
* `strand_population` keep-or-cut (measured inert; owner call).
* `npmle.py` retirement scope (still LIVE at `calibrate.py:591`; owner sign-off).
* the VOCABULARY RENAME (`splice_in`/`peel`/`pin`; dedicated effort; do not rename piecemeal).

## THE DISCIPLINE THAT PAID (keep doing exactly this)

* The loop: worst condition → rank by error MASS → concentration → flags → channel bisect → HOP WALK →
  READ the code line. The hop walk is what named the pin; the flags alone pointed at the SPLICE IN.
* Score the mechanism where it LANDS, not where it is written — the signature census at the delivered
  destination read 2.6 % of the error (`TRAPS: a-scorer-scoped-to-the-mechanisms-targets`, second form).
* `TRAPS: panel-before-src` fired a sixth time: `pinM` closes the zero-control gap outright and costs
  stranded × ON six of six. Set the bar in writing BEFORE the panel returns.
* Baseline arms on the COMMITTED tree, candidate arms after the edit, ONE scorer, tree stamped on every row.
