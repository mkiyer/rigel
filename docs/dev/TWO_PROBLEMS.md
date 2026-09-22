# The two problems (owner, 2026-09-21; restated 2026-09-22)

The work on the tool, after the cleanup, is two problems. Stated plainly:

1. **Estimate a prior for each transcript.**
2. **Correctly estimate the capture-contracted effective length for each transcript.**

That is all. If a derivation of either needs more than a page, it is off track.

Six words, meaning exactly this:

| word | meaning |
|---|---|
| **capture efficiency** | of a region or boundary: the fraction of fully captured that gDNA measures there (`gdna_capture_efficiency_*`) |
| **capture-contracted length** | of a component: the number of fragments it yields per unit abundance under capture — the EM divides its count by it |
| **the gDNA component** | the one hypothesis per locus that a fragment came from genomic DNA |
| **the synthetic spans** | one unspliced template per gene, the stand-in for RNA that has not spliced |
| **the pseudocount prior** | the two numbers per locus the EM adds to its counts today, `P_g` to gDNA and `P_R` to RNA |
| **the per-fragment gDNA prior** | a proposal for problem 1's gDNA side: each unspliced fragment carries its own probability of being gDNA, read from calibration at the objects it sits on |

## What production does today, and where it is measured to be wrong

An unspliced fragment at a locus is assigned among the gDNA component, the isoforms and the synthetic spans by the
E-step weight `count_c / length_c × (fragment terms)`. Two inputs come from calibration: each component's
capture-contracted length, and the pseudocount prior. Three things are measured wrong, each an `ISSUES` entry with
its numbers:

1. **The pseudocount prior is biased toward gDNA** — `P_R` is spread over the whole RNA pool, spliced fragments
   included, so the odds the EM starts from carry a factor `(1 + P_g/G) / (1 + P_R/R)` that is not 1
   (`ISSUES: the-pseudocount-prior-is-biased-toward-gdna`).
2. **The gDNA component's length is computed by a different rule than the transcripts'** and the two disagree by
   14–16 % on capture; on capture this error and the first cancel, so neither can be repaired alone
   (`ISSUES: the-gdna-component-length-rule-differs-from-the-transcripts`).
3. **Spliced and unspliced templates disagree by 3–5 % on capture** under any per-region efficiency — the junction
   problem (`ISSUES: ruler-witness-geometry-on-transcript-panels`).

## What is known that bears on the two problems (measured, with the instrument)

* Calibration's per-object gDNA fraction is a calibrated probability at the panel's depth: on every in-scope ladder
  row the realised gDNA share of the fragments deposited on objects publishing a fraction `p` reads `0.93–1.12 p`
  in every populated bin, the whole-library gDNA bias within ±2.7 %. Its shallow objects (five fragments or fewer)
  understate gDNA by 15–40 %, on unstranded rows as much as stranded ones. At `g98` it is right to 0.5–0.8 % of the
  gDNA, which is 32–78 % of the RNA; on the cfRNA libraries its per-object precision is near zero and 15–44 % of the
  unspliced incidence sits on shallow objects. (Measured 2026-09-22 with a since-retired instrument; recorded here
  because it sizes what problem 1's gDNA side can lean on: calibration at depth, and nothing at `g98` or on sparse
  real data without pooled evidence.)
* The per-transcript prior lane (`rna_prior_weight`) is plumbed end to end and nothing fills it
  (`ISSUES: per-transcript-prior-lane`); the allocation ceiling recovers most of the in-scope transcript error
  (`quant_accuracy.py --arm oracle_alloc_seed`).
* The one-scale read-out for problem 2 exists: `ruler_vs_truth.py --scale`, read within each probed class.

## The order

The cleanup first (`docs/dev/CLEANUP.md`): a production-ready tree with the shipped infrastructure and no new
mechanism. Then problem 2 (the length, one rule for every component, judged per class against the simulator's
yield with no EM), then problem 1 (the prior), each derived on one page and A/B'd against what ships.
