# Nascent / mature EM effective-length: bug investigation & the capture identifiability limit

**Status:** 2026-06-14. Records the dissection that led to a confirmed eff-len bug fix
(`capture_eff_length.py`, the exon→region incidence) and the decisive finding that the *residual*
nascent leak under hybrid capture is an identifiability limit, not a code bug.

## 1. The dissection (what started it)

Ranking the flagship (`gdna300/ss0.99/capture_on`) loci by total transcript error (3-pool:
gDNA / nRNA / mRNA, `scripts/debug/dissect_loci.py`) showed the dominant artifact is a **nascent
sink**: in a library with *zero* true nascent, the synthetic nRNA spans absorb large mass from both
gDNA and mature. It is present in every capture-on condition (and even capture-off), and it appears
even in a **pure-RNA** (`gdna_none`) library — so it is an EM/geometry artifact, not a gDNA-calibration
error. A prior-vs-truth split (`scripts/debug/prior_vs_truth.py`) confirmed the calibration prior is already
~97% accurate at ss0.99 — the flagship leak is **EM-side**.

## 2. The confirmed bug — exon→region incidence off-by-one (FIXED)

A nascent span's footprint (exons + introns) must always be ≥ its mature analog's (exons only), so its
EM effective length must be ≥ the mature's. It was not: `nascent em_eff = 2928 < mature 3703` at
flagship locus 2. An adversarial 3-agent audit traced the root cause to
`_exon_region_incidence` ([capture_eff_length.py](../../src/rigel/calibration/capture_eff_length.py)):
it mapped an exon `[a, b)` to regions via `searchsorted(region_starts, a, side="left")`, which **skips
the region that contains `a`** whenever `a` falls in a region's interior. Region boundaries do **not**
always coincide with exon edges — `build_region_partition` merges adjacent *same-signature* segments,
so a shorter alternative-transcript exon starts interior to a merged region (only ~90% of exon starts
coincide with a region start). Consequences on the flagship index: **46 single-exon nRNA spans dropped
entirely** (zero incidence → no contraction) and **165/432 mature transcripts got `len_t < exonic`**
(geometrically impossible).

**Fix:** read the lower bound from region *ends* — `lo = searchsorted(region_ends, a, side="right")`
(containment semantics, matching `tests/native/_accumulator_reference.py`). Verified: 0 dropped,
`len_t/exonic ≥ 1.0` everywhere. Regression test: `tests/calibration/test_capture_eff_length.py`
(fails on the old code, passes on the fix). All unit + golden tests pass with zero drift.

**Effect** (incidence fix on the baseline IPR contraction; 16-condition soft-mass sweep):
flagship no-nascent phantom nascent **313K → ~88K**; real nascent **recovered** (`nrna_rnd`
ss0.99 cap_on 168K/190.6K = 88%, cap_off 178K/190.6K = 93%); gDNA recovery improved
(−9.6% → −8.3%); capture-off clean. Better than the committed (buggy) baseline on *both* the phantom
sink and real-nascent recovery.

## 3. The residual is an identifiability limit, not a bug

The only evidence separating nascent from gDNA is **intron coverage** (gDNA is unstranded-uniform;
mature has no introns; spliced reads are mature-only). Measured unspliced read density in intron-class
regions, `gdna_none` so intron reads = nascent-only (`scripts/debug/intron_signal.py` logic):

| condition | intron density | intron/exon |
|---|---|---|
| real-nascent, capture-**off** | 0.098 | 0.050 |
| real-nascent, capture-**on**  | 0.0013 | 0.0006 |

Hybrid capture **crushes the intron signal ~75×** (probes target exons; a nascent transcript's reads
pile up at exons, indistinguishable from mature/gDNA exonic reads). So under capture, nascent vs gDNA
vs unspliced-mature for an exonic read is **fundamentally unidentifiable** — no eff-len or prior
recovers what the probes deleted. The EM effective length only chooses *which way to be wrong*:
contract nascent → phantom sink; lengthen it → nascent→gDNA siphon + real nascent lost (verified: the
"don't contract" extreme makes real-nascent recovery *worse* and re-opens the siphon). The corrected
(incidence-fixed) IPR is a reasonable operating point: it recovers real nascent where identifiable
(capture-off) and the residual ~88K capture-on phantom is the identifiability floor.

## 4. Rejected: the global enrichment-ratio contraction

An `eff = fl · ρ_t/ρ_ref` enrichment-ratio model was explored as a monotone replacement for the
(non-monotone) IPR contraction factor. Every global reference tried regressed: `ρ_global` (simple
mean) over-contracts fully-exonic mature ~12.8× (gDNA recovery 91%→71%); the mass-weighted mean
over-contracts everything (phantom sink returns); the exon-class density is **annotation-derived**
(probes need not sit on exons — using exon/intron labels to define "enriched" is not data-driven and
was rejected). A correct contraction must key on a component's *internal* density heterogeneity
(which the IPR already does) and be monotone in the footprint; that refinement is deferred.

## 5. Open fronts (separate from this fix)

- **gDNA→mature under capture** — with the phantom nascent removed, the freed exonic-unspliced mass
  partly lands on mature (mRNA overshoots +11-33% capture-on). This is the EM gDNA-vs-mature
  identifiability front (`em_gdna_strand_likelihood_fix.md`).
- **ss=0.5 capture-on** — unstranded + capture is the worst identifiability and also carries the
  count-mean calibration bias (the calibration prior under-calls gDNA ~14% there); the simplex-sweep
  / count-channel work in `CALIBRATION_PLAN.md` targets this.
- **Monotone, data-driven contraction model** — a per-component internal-concentration contraction
  that restores the nascent⊇mature invariant structurally without an annotation reference.
