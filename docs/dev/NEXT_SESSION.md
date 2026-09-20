# NEXT SESSION — the isoform allocation, in three gated steps (2026-09-20, end of day)

This file is only how to begin. The ranked view is `docs/ROADMAP.md`, the open problems are `docs/ISSUES.md`,
the rulings and the record are `docs/DESIGN.md`, what "done" means is `docs/SUCCESS.md`, the lessons are
`docs/TRAPS.md` cited by name.

## Where the tool is

`main` is at `c52c9b93`, pushed. The nascent siphon under capture is repaired: the gDNA component's opportunity
converts its crossing support by the count's own `q` (`EQUATIONS.md` §11), `g50 ss.99 ON` reads a siphon of
+32,905 where it read +541,216, and the three calibration instruments are identical to before it. The suite is
**0 failed / 3,452 passed / 0 skipped / 2 xfail, 3,454 collected**. The report at the owner's link renders the
committed tree (`arms/qa_ladder_base_q.jsonl`; the `*_q` arms are the current ones, the unsuffixed arms the
2026-09-19 tree's).

## Three rulings of 2026-09-20 to carry

* Rigel stays PANEL-AGNOSTIC and takes no panel input; the junction witness is post-release
  (`ISSUES: ruler-witness-geometry-on-transcript-panels`, deferred).
* The pre-EM prior chain is deferred: a perfect prior moves nothing in scope on the committed tree.
* `g98`'s RNA prior floor is a composition-solve estimator item at the stress rung, not a quick fix
  (`ISSUES: rna-prior-floor-at-pure-gdna-loci`), post-release.

## The job: the isoform allocation — `ISSUES: per-transcript-prior-lane`

The ceiling, on the committed tree (`quant_accuracy.py --arm oracle_alloc_seed`, transcript Σ|Δ| as a share of the
true annotated RNA at `g00` / `g05` / `g50` / `g98`): stranded OFF 2.01 / 1.59 / 2.49 / 15.70 → 0.43 / 0.43 / 0.70 /
6.38 %; unstranded OFF 1.73 / 1.93 / 2.55 / 18.35 → 0.42 / 0.45 / 0.81 / 8.35 %; stranded ON 6.53 / 3.49 / 5.50 /
31.03 → 2.91 / 1.08 / 2.07 / 18.00 %. That arm hands over the TRUE support; it is a ceiling, never headroom.

1. **The support probe as an ARM, nothing in `src/`.** An isoform is supported iff its exclusive spliced
   evidence says so — its exclusive junctions carry reads at a depth where its siblings' do — and every other
   weight stays the EM's own. Build it beside `oracle_alloc_seed` in `quant_accuracy.py` (the arm machinery
   already feeds a per-transcript array through `rna_prior_weight`), score it per stratum above the
   `base_reseed` floor against `base` and the ceiling, keep `g98` apart, never pool. Derive the rule on paper
   first: what "exclusive", "reads" and "depth where its siblings' do" mean is the whole design, and there are
   no magic numbers — a threshold is a ratio to a sibling's evidence, never a count.
2. **Only if ① recovers a meaningful fraction at `g05` and `g50`** (a third of the ceiling would already be the
   largest in-scope gain available): the producer for `rna_prior_weight` and the kernel's own-count fallback for
   components whose structure cannot speak (`em_solver.cpp`, `apply_grouped_prior_update`; `EQUATIONS.md` §9b.1
   says why a flat share is refused and why `raw[i] = 0` must stay absorbing). ⛔ The lane is one static array
   and filling it wholesale reallocates the WHOLE RNA pseudocount (5.21 → 53.37 % once, 2026-09-19). A
   falsification test first, verified failing, then the fixed code broken and every gate watched to fire; A/B
   against `base` on `g50 ss.99 ON` and `OFF` first, then the ladder, then the report with the `ladder-report`
   skill; snapshots for the owner's go.
3. **If ① recovers little:** the census of which genes carry the error — the near-tied multi-isoform genes
   (`quant_accuracy.py` per gene, `worst_objects.py`) — and stop to discuss before designing anything larger.

## Where everything is

* The session's harnesses, traces with true origins, the truth E-steps and the A/Bs of 2026-09-20:
  `~/Downloads/rigel_runs/prototypes/2026-09-20_siphon_mechanism/`.
* The 0.10 condition: `~/Downloads/rigel_runs/suite/ladder_nrna_lo/` (one condition, cached, certified).
* The teaching page: the Artifact "The Twice-Counted Crossing" (the siphon's mechanism and the ruler's).
