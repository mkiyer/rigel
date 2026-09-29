# ISSUES — the issue log

This file is the issue log: one entry per open problem, question, decision or risk, and an append-only record
of what was measured and turned down. An OPEN entry is a `### kebab-name` heading with a `priority:` line (now
/ next / later / parked), the question in a sentence or two, the numbers a ranking turns on, and the
instrument that re-derives them. The OPEN section runs now, then next, later and parked, and within a tier it
follows `ROADMAP.md`'s order. A CLOSED / REFUSED entry keeps its name, its verdict, the killing number and
the date, so a refused mechanism is not rebuilt. Cite an entry as `ISSUES: <name>`; names are the only
identifiers (`tests/test_docs_boundary.py`). What does not belong here: the ranked view (`ROADMAP.md`),
rulings and derivations (`DESIGN.md`, `EQUATIONS.md`), lessons (`TRAPS.md`), and any record of what was done —
the changelog is git.

---

## OPEN

Ordered by priority. An entry says what is open and the number a ranking turns on; what was done is git,
what was ruled is `DESIGN.md`.

### splicing-artifacts
`priority: now — the cluster track, the owner's highest priority (2026-09-26), a bigger problem than previously thought (owner, 2026-09-24) · kind: defect · 2026-09-24`
A gDNA read that the aligner splices across an annotated junction reads as certified RNA unless the index's
blacklist — junctions the same aligner wrote on simulated genomic reads, kept at two or more, with the longest
anchor seen — rejects it. Today a fully rejected fragment is labelled artifact: the EM treats it as unspliced
(`DESIGN.md` §3.1d) but its blocks stay split at the rejected junction, so its footprint still spans the intron
and over-states gDNA's length when the short anchor is the outer block; calibration holds it out entirely. The
owner's first repair to consider: edit the alignment — delete the rejected junction and keep only the largest
aligned block — so the fragment is unspliced in both stages at its true length. Beyond that: the blacklist is a
yes/no cut (count ≥ 2, anchor ≤ the longest seen) on a continuous quantity, the aligner's rate of writing that
junction on genomic sequence at that anchor length. Measurable only on aligned data: the aligned ladder is drafted
for the cluster (`~/Downloads/rigel_runs/aligned_ladder/`), and on the production store 23.3 % of the ladder's
annotated junctions are blacklisted. Also here, from
`ISSUES: calibration-and-the-em-disagree-on-what-can-be-gdna` (closed): three fragment types calibration still puts on
its unspliced bank, gDNA-eligible, that the EM treats as RNA only — an unannotated sequenced junction; two annotated
junctions on opposite motif strands; mates disagreeing on an acceptor, which calibration merges into an unannotated
intron and the resolver calls annotated — and one length rule the stages apply differently: calibration drops an
unspliced or artifact reading past the maximum fragment length, the EM keeps gDNA while its likelihood competes.
Those three fragment types are a local A/B, independent of the cluster work.
- BOTH WAYS. On the VCaP RNA half the production blacklist rejected 598,010 records where a pure-gDNA library's rate
  predicts about 20,000 (uniquely mapped single-junction records 429,208 of 12.55 M, 3.4 %; multi-junction 154,604 of
  1.88 M, 8.2 %). A rejected genuine read becomes an unspliced, gDNA-eligible fragment whose footprint spans the
  intron, so the cost lands on the gDNA and synthetic-span pools of exactly the listed genes. A repair is judged on
  both errors together.
- THE LARGEST ARTIFACT CLASS IS DECIDED BY THE REFERENCE, through h, the reference mismatches between the spliced and
  unspliced placements. On the DNA library at annotated junctions with a 1–10 bp short side, 96 % of artifacts have
  h = 0 (5,092 of 5,297) and 81 % sit on the acceptor side (4,272), where STAR's sjdbScore breaks the tie toward
  splicing; the blacklist catches them only where alignable sampled the junction (4,643 rejected, 449 h = 0 records
  escaped). On the RNA library 40,557 rejected records are h = 0 and 162,213 need h ≥ 1. The hypothesis to derive on
  one page, the multi-junction case included: an ambiguous junction ADDS an unspliced reading weighted about
  (ε/3)^h, ε from the library, h from the reference (at index time for annotated junctions, at scan time for novel
  ones); the owner's "edit the alignment" is its gDNA limit at h = 0. Rigel reads NM only — never MD, the sequence,
  the qualities or the reference — so this needs new scanner and index infrastructure.
- THE LONG-ANCHOR ARTIFACTS ARE NOT h = 0 (short side ≥ 21 bp: 4,840 DNA records, 60 % uniquely mapped; 98 % of the RNA
  library's rejections, 235,406 of 240,520, uniquely mapped). The hypothesis is remote origin, from retrocopies and
  paralogs: alignable records an artifact from every alignment, secondary and supplementary included, so a parent
  gene's junction is listed with long anchors on both sides and, under the OR rule, rejects every single-junction read
  up to about 2m aligned bases. Measuring it needs an unspliced end-to-end realignment, alignable's per-record origins
  and a junction mappability table; the EM side needs `ISSUES: multimapper-intergenic-alignments`.
- THE BUILDER applies `min_count` (2, underived) per (chrom, intron, strand, read_length) BEFORE aggregating over read
  lengths, so a junction seen once at each of several read lengths never enters, and the shipped feather keeps no
  count, so no threshold can be decided at run time. Fixed with the catalogue rebuild.
- UNTESTED RULES, to pin with local falsification tests before any rule changes: the `<=`, the OR of the anchors, the
  inert stored 0, cumulative anchors across earlier introns, and the mixed fragment (one junction rejected, one kept,
  deposited as spliced with the rejected N as a gap); every fixture's anchors (500–10,000 bp) are longer than any
  read. Undecided: the reject runs per record before the mates are joined, so a junction rejected on one mate and
  kept on the other survives, and any surviving junction beats ARTIFACT, so one kept junction shields the rejected
  ones.
- THE TRAINING SETS. The spliced strand table that sets κ and supplies the od's RNA junction seeds holds artifacts
  whose sense fraction is ½ (the VCaP exome-DNA half: 3,223 spliced observations on 2,255 junctions, κ 0.5008, feeding
  the RNA od with information 5,701), and κ pools per-junction rates that vary with depth
  (`ISSUES: strand-overdispersion-one-shared-value`, step 4), so on real libraries the RNA strand MEAN is itself off:
  LBX0588's 0.9365 against MO_3021's 0.998 is plausibly artifacts (unmeasured). The repair trains the strand model and
  the RNA length law weighted by each fragment's probability of being genuine.
- THE COUNTERS ARE IN THREE UNITS: `sj_blacklisted` counts junction × record over every record, multimappers
  included; the `splice.*` census counts uniquely mapped fragments; 0.7.1's `summary.json` splice counts are the
  length models' observations. So the 344-library cohort's 0.7.1 summaries do not compare with today's census.
- THE WORK runs on the cluster in phases: the substrate (per-record splice rows, a working index checksummed against
  production and built after the index format bump, `ISSUES: the-format-changes-to-batch-before-release` (b),
  `--notemp` re-runs, more truth mixes); both errors measured on today's tree; the mechanism census; the h-weighted
  reading (`ISSUES: scoring-penalties-are-underived-constants`); the training sets; aligner-settings robustness and the
  catalogue rebuild. The index case's splice evidence on the cluster's `/scratch` is gone (owner, 2026-09-28): it is
  regenerated from scratch by a new run, a separate task. The VCaP mix (its BAM, 0.7.1 `annotated.bam` and
  `summary.json`) and LBX0190 / LBX0588 / MO_3021 (the same, with `SJ.out.tab.gz`) may exist only in
  `~/Downloads/rigel_runs/cfrna/`: confirm cluster copies before the substrate phase and upload any that are missing.

### an-unrecorded-splice-blacklist-is-dropped-silently
`priority: now — the cluster track: check the manifest before the next cluster quant (owner, 2026-09-28) · kind: risk · 2026-09-28`
`TranscriptIndex.load` applies `splice_blacklist.feather` only when `manifest.json` records `sources.alignable_zarr`;
otherwise it sets `sj_blacklist_size = 0` with no log line. The cluster's production index
(`hulkrna/refs/human/rigel_index/`) holds a 5.8 M-row feather: if it was placed by hand, the next quant there runs
with splice-artifact detection off and no error — format 8 (eeba4b0a) predates the sources block (8c673324), with no
version bump. The only trace is `sj_blacklist_loaded: false` in `summary.json`. Before any cluster quant, check the
manifest or rebuild with `--alignable-zarr`. The one-line warning for a feather without a recorded store is in
`ISSUES: latent-defects`.
Instrument: the manifest; `summary.json`'s `sj_blacklist_loaded`.

### strand-overdispersion-one-shared-value
`priority: now — Tier 0, the owner's first (2026-09-28): steps 0–3, then the fl.py fix (ISSUES: the-realized-gdna-length-law-reads-rna-counts), then the field collapse; steps 4–6 are research · kind: design · 2026-09-27`
Calibration's strand overdispersion (od, `EQUATIONS.md` §6) reaches the solve through `reconcile_overdispersions`, fed
two mismatched value/precision pairs: gDNA's ρ = 0 pair-count moment beside the fitted estimate's precision (the
influence-weighted fit is computed, then discarded), and RNA's unweighted moment beside its null information, which
credits it about 700,000× too much evidence on a real library (MO_3021 looks calm only because RNA is credited 1.0e7).
The rule keeps the two apart unless one is far better measured — the weaker is pulled 16 / 50 / 159 of its own SEs at
S = 2k / 20k / 200k — which is not what the owner intended. On `odg05` (planted gDNA 0.05, RNA 0) the shipped gDNA
value is ≤ 0.0007 on all 12 g05–g50 rows against 0.038–0.058 for the gDNA fit alone; at g98 the RNA value reaches
0.049 against a truth of 0; within one 1 % draw the pair sits up to 0.19 apart.

OWNER RULINGS (2026-09-28):
1. gDNA and RNA share ONE overdispersion.
2. First prototype the gDNA-evidence fix on top of the shared estimator, with a new panel condition that plants
   gDNA alongside RNA on the opposite strand of its seeds; then land the two in sequence, each A/B'd alone.
3. The no-evidence fallback becomes 0 (binomial), once it is justified theoretically and, ideally, empirically.
4. The shared estimating equation can have two roots at low gDNA; it needs a clear, simple rule.
5. No hard-coded upper limit, but no estimate may be able to sabotage the tool.

Accepted with the triage (owner, 2026-09-28):
6. Now: steps 0–3, the fl fix and the field collapse. Steps 4–6 are research, after the splice training sets or with
   them.
7. Ruling 1 stands at its measured price (below), reopened only if the RNA-od arm shows harm.
8. Step 1 lands before any evidence fix, under the hold rule.
9. Ruling 3's 0 rests on theory, with forced-fallback injections as its empirical half; `summary.json` gains an
   information and evidence label.
10. Ruling 4's rule is the least root, with no cap on its passes; the pass count and each seed set's own root are
    logged as diagnostics.
11. The 0.2 clip stays, labelled as a clamp, until steps 4–6.
12. Step 1 skips the fit where the strand channel is dead and reports the od as not measured.
13. The simulator plants a matched RNA od (a per-transcript Beta sense rate) before step 1. At ρ_r = 0 it draws
    nothing, so every existing panel stays bit-identical.
14. The panel: an exon-edge shadow class placed so that no seed spans two rates; capture ON priced at g25 only
    (λ ≈ 3.5, +3.1 % RNA); one contaminated fraction for step 1 and a second before step 5; ss 0.70 kept; gene-edge
    seeds clean; the panel stays outside the tree until step 1's verdict.
15. The robustness rule reads the gDNA-total spread across seeds at 10 % depth. The vetoes are
    `calibration_vs_oracle.py` and `quant_accuracy.py` per in-scope stratum; the od spread is a diagnostic.

THE RECORD IS FALSE TODAY. DESIGN §3.3a and EQUATIONS §6b–§6c describe a weighted RNA fit and fitted-precision pairs
that never shipped, and DESIGN's "lands on the oracle value exactly" was measured on mixed code. Three stated reasons
are wrong:
- The −½·log var contrast given for the reconcile and the 0.2 fallback is vacuous since the strand variance is frozen
  at the reference composition: od sets a width only. This is ruling 3's theoretical half.
- EQUATIONS §6a says both "unbiased under any distribution of RNA content" and "so it biases down, never up". Seeds
  mixing gDNA with sense RNA do bias the away-half moment LOW, because away-half selection makes the cross pairs
  negative: a planted 0.03 reads 0.0264 ± 0.0021 against 0.0303 on truly pure seeds, weighting by the true pair
  fraction still leaves −5 to −8 %, and mixed seeds are unobservable, so no estimator repair exists.
- "A wrong weight costs efficiency, never correctness" (`gdna_strand.py`) is false under contamination: the weights
  move odg05's joint root 0.0406 → 0.0278, while the gDNA side alone moves 0.0446 → 0.0429.

THE LANDING. Every step re-records the baseline on the current tree; reads `calibration_vs_oracle.py` per stratum,
beside `ruler_vs_truth.py` on capture-ON arms, and `quant_accuracy.py --set em.assignment_mode=fractional`; verifies
each gate failing, then breaks the fix and watches it fire; reads the golden diff's size before regenerating.

STEP 0, THE INSTRUMENT (designed on paper; no file yet). The contaminated-seed panel condition: 15 hosts (gB2 / gB3 /
gB4 × 5, 9 of them probed), each with one opposite-strand single-exon shadow; 90 of 1,782 count-observable seeds
contaminated (5.1 %), shadows at 1.23 % of RNA; 40 conditions at 1.17 M fragments on their own `RIGEL_SCRATCH`.
- Sizing, 100 replicates at capture OFF: od00, the joint fit reads 0.118 against 0 at g05 (z 17), 0.036 at g25 (z 12),
  0.0099 at g50 (z 15); od05, G reads 0.164 against 0.049 at g05 (z 16). Capture ON resolves only at g05 (z 15; g25
  z 4.0; g50 null).
- Gates: each of the five shadow refusals broken on purpose; `calibration_oracle.py` accepts every row; each host's π
  within binomial noise (per host under capture); the od05 capture-OFF target G in 0.048–0.051; ss 0.50 must not
  rise; FASTA, index and cache keys unchanged. Contaminated slots are scored apart.
- Its first reading tests ruling 2's premise — seeds with more than ¾ of their fragments away carry 49–79 % of the
  pair weight on the real libraries, against 97 / 34 / 1.8 % at the panel's g05 / g25 / g50 — and whether the lower
  root is the target on the capture-ON g05 rows. The numpy-only MO_3021 check (step 5) can run here.
- The RNA-od arm (ruling 13): today every panel scores the shared value at an RNA od of 0, so odg05 has no single
  right answer.
- The planted gDNA od is drawn per simulator region, so a seed straddling regions has a lower ICC than planted and
  the oracle arm is not an ideal seed-level target: the fits read 0.0446 ± 0.0059 at g05 ss.99 ON and 0.0516 ± 0.0047
  at g25 ss.99 OFF against 0.05. Size unmeasured; fix it with the RNA arm if the first reading needs it.

STEP 1, THE JOINT FIT UNDER THE CLIP — NOT a numeric no-op. One influence-weighted fit over the gDNA seeds and the RNA
junctions replaces the reconcile, still clipped at 0.2, with CLAMPED keyed on the fitted value. MEASURED (prototype,
2026-09-27): odg05 calibration −10.1 % / −18.7 % and transcripts −6,370 / −28,829, −3 % / −7 % (stranded OFF / ON);
ladder and test panel within noise (a fitted feed moves 3 of 16 ladder rows by about 1e-6). On real low-gDNA stranded
libraries the gDNA seeds are mostly opposite-strand RNA (raw moments 0.4–0.86; VCaP transcriptome: 118,622 away-half
seeds, moment 0.47), so the joint fit
equals the 0.2 value on 13 of 14 real runs with identical `rigel quant` output (Σ|Δ| 2,280.8 in both); against
shipped it moves calibration gDNA −15.4 % (VCaP RNA) and −5.1 % (MO_3021). HOLD if, on the panel's contaminated od00
rows, the joint fit costs more than shipped in a pinned A/B pair, judged by size.
- With it: key the log's "own-evidence od" labels and CLAMPED on the shipped value, not the raw moment (MO_3021: raw
  0.675, fitted 0.1327); key `clamped_at_ceiling` and `effective_seeds` (set from the discarded fit, read nowhere) the
  same way, then drop ROADMAP's "read neither until step 1" caveat; rewrite the ~15 docstrings in `gdna_strand.py` and
  `calibrate._fit_strand` that describe no weighting, a null information, prior shrinkage and a closed form, and the
  "per-sj SJ strand table" stutter, also in a runtime error message; correct EQUATIONS §6a–§6c and DESIGN §3.3a (the
  reconcile ruling is superseded, the ceiling ruling is not) and the owner's auto-memory note on the
  strand-overdispersion estimator, which describes the influence-weighted design that never shipped; append the
  reconcile's refusal to CLOSED with the numbers above; check the pipeline's "reads as unstranded" warning, which
  re-derives κ instead of reading calibration's decision; log the chosen root, the pass count and each seed set's own
  root as diagnostics only.
- Gates: the prototype's 26 reduction cases; two disagreeing sides return one value; zero gDNA; zero RNA; unequal
  depth (the raw feed gives 0.1667, the fitted 0.0139); the harness's exact-variance check (Var(e) against
  enumeration, 2.6e-12), the gDNA-only reduction (bit-exact 26/26), the RNA-only fit against brentq (2.8e-17), the tie
  term and the between-seed quadrature; a pin on `GdnaStrandModel.information` when seeds exist. Today's recovery test
  builds equal od in both components, so the reconcile is a no-op in it; restate or retire the HANDFUL gate
  (`ISSUES: hygiene-ledger` (e)).
- The dead-channel skip (ruling 12): at κ = ½ the strand term is flat and `strand_evidence` is multiplied by disc = 0;
  on the ladder's unstranded OFF rows the 0.2 value moves 0.014 fragments, and on odg05 g05 unstranded the joint fit
  reads 0.0003 against a gDNA fit of 0.042.

STEP 2, THE EMPTY-EVIDENCE PREDICATE — a numeric no-op. The empty case moves inside the shipped predicate (den > 0 and
clip(num/den) − ρ > 0), so an empty seed set returns exactly 0.0 and every fit with evidence stays bit-identical (437
of 437 seed sets). The theory, on the minimal ψ relative to the true od: 0 is at most ×1.07 worse (best in 32 of 48
cells), ρ_eff at R = 0.2 ×2.56, 0.2 ×10. A wrong value costs: 0.2 where the truth is 0, +255 % (test panel ON),
+199.8 % (ladder ON, 436,286 → 1,307,989), +18.1 % (stranded OFF); 0 where the truth is 0.05, up to +19.8 % (odg05
stranded OFF), and at g98 an under-call of 9.6k / 10.3k. The branch has never fired: 1,834 fits over 136 substrates
(125 both sides, 10 RNA only, 1 gDNA only, 0 neither). Unrun: odg05 at d10 and d100 with od injected at {0, 0.05,
0.2}; the per-slot condition β = J₀b² ≥ nρ; the bridge "fallback firing implies no gDNA-dense slot", which holds only
for today's away-half rule (a mostly-gDNA seed lands on the RNA side with probability 0.84 at f = 0.8, n = 200) and is
re-derived after step 5; a junction-poor deep library, which does not exist (birthday thresholds m from 10 to 183).
The information and evidence label are in `ISSUES: the-format-changes-to-batch-before-release`.

STEP 3, THE LEAST ROOT on [0, ρ_max]: the shipped bisection plus a certificate, no constant. At low gDNA the joint
equation has two stable roots and the bisection's search path takes the upper one (≈0.0002 and 0.036 at odg05 g05
ss.99 capture OFF; 0 and 0.040 ON). 59 of 61 seed files match bit for bit; on the two two-root rows the upper root
moves transcripts −1,697 (within the floors) and −8,622 (unresolved: summed floor 2,398, worst single row 12,668), and
genes +1,525 against floors ≤ 15. The search takes a median 119 passes (max 314, about 0.03 s) against the bisection's
1,072 at root 0, growing as ε^−½ near a fold — about 10⁵ passes (about 20 s) at 1e-6 from an edge, P(passes > N) ∝ N⁻²
— accepted uncapped (ruling 10). Gates: the band-edge cases, the one-root reduction cases, the clamp at the upper end.
A/B on odg05 g05 ss.99 and the panel's capture-ON g05 od00 rows (two roots in 30 of 30 replicates, the lower the
target).

THEN, BIT-IDENTICALLY, ONE FIELD. Collapse the gDNA and RNA fields in the `_Strand` triple, `CalibrationResult`, the
cli summary, `sweep`'s od_g / od_r and `policy_od_*`, and the native Grid's `pol_od_*`. After it: re-read the strand-input
drift in `ISSUES: the-gdna-landscape-collapses-at-low-depth` and `ISSUES: strand-likelihood-over-confident-beyond-od`,
and refresh the ladder report.

STEP 4, THE JUNCTION WITNESS — NOT DERIVED. A junction's antisense rate depends on its depth (VCaP 22κ at 2–3
fragments, about κ/4 at ≥ 128; LBX0588 0.14 at 2–3, 0–0.01 at ≥ 16, against κ 0.064), so a single-rate fit reads a
wrong mean as dispersion: the misfit alone gives M_r(1) = 1.27–2.17 against ≤ 0.17 for a single rate, and deep VCaP
junctions (≥ 128) read −0.0038 ± 0.00027, about −1/(n−1). Alone, the junction set has the two roots {0, 1} (VCaP RNA at
10 %, all three seeds; MO_3021 at 1 %), chosen by the sign of a near-zero moment (M_r(0) −0.005 to +0.002 against a
bootstrap SD of 0.0014–0.0044): a hidden yes/no gate. It swings with depth (VCaP RNA 0.0 / 0.2 / 0.2 at 1 % / 10 % /
full; resampled at full depth, 0.55 on 76–80 % of draws) and reads above a planted 0 (d100 g00 OFF: 0.065). Model the
misfit, never threshold it. The junction seeds and κ carry splice artifacts, which the training-sets work of
`ISSUES: splicing-artifacts` retrains. Gate: a bootstrap distribution on every real substrate. Falsifier: M_r(1) > 1
on VCaP after the repair.

STEP 5, A gDNA PURITY WEIGHT (after steps 3–4 and ruling 5): weight each seed by h_s = E[p_s² | counts]. Five arms —
the joint fit, the posterior E[p²], the true pair fraction, truly-pure seeds only, shipped — ranked on the panel's G
readout, never on "the clip stops binding": the Jeffreys E[p²] ≈ 0.38·√μ_a acts as a depth gate (LBX0588's gDNA weight
cut 14× at full depth, 40× at 10 %), and taking weight off the gDNA side hands the fit to the unstable junction witness.
- Failing today: the own count against λ_off under capture (truly pure intron seeds at odg05 g05 ON read π̂² 0.10
  against 0.98, `ISSUES: intron-seeds-near-probes-are-capture-enriched`); transport from a neighbour (77 % / 82 % of
  antisense exon-flank boundaries have an empty neighbour, VCaP RNA 10 % / LBX0190).
- Open, with its numbers in `04_evidence.md`: where the weight gets capture information (a refit after the solve from
  the located enriched mode, unweighted seeds, or fl's transport without the ½); weight or coefficient; a
  κ-continuous weight (p² + (1−p)² is not derived); where λ_off lives (fl's one-sided rate is QC-only and would break
  `fit_gdna_strand_from_substrate`'s no-density contract). The per-object efficiencies may not be used: the solve
  would certify the contaminated seeds as pure.
- The premise first: LBX0588's moment rises with seed depth even among seeds ≥ 90 % gDNA (0.006 / 0.012 / 0.020 at n
  in [2,5) / [5,20) / [20,100), median away fraction 0.50), and only 26 % of its excess comes from seeds with an away
  fraction ≥ 0.9, against 87–100 % on MO_3021, VCaP RNA and LBX0190. Stratifying MO_3021's census by λ_off·E (95.5 %
  of its excess at an away fraction ≥ 0.9) is a numpy-only test of whether antisense explains it. Unpriced: real gDNA
  together with contaminated seeds. Gates: the evidence page's four, and a rerun of the blank-chromosome shadow
  control.

STEP 6, REPLACE THE CLIP, then delete `_MAX_OVERDISPERSION` and `_CEIL_ALPHA_BETA`. On [0, 1] the joint fit reads
0.512 (VCaP RNA full), 0.53–0.86 (VCaP at 10 % / 1 %), 0.23–0.75 (MO_3021) and 0.820 (LBX0190); one pair hits the
ceiling with probability ⅓ at ρ = 0 given an away-side pair; the clip binds on 14 of 61 surveyed rows, and the other 47
are bit-identical under it. Candidates: an agreement rule on the repaired witnesses, or a two-width strand mixture (a
native change and a learned weight).
- Checks: the bootstrap spread inside the fixed-value floor; a Monte Carlo of the measured misfits; the harness on the
  ladder, odg05 and the new panel; exact Godambe weights, since §6b's per-seed variance drops a sampling term and two
  between-seed terms and is validated only for ρ ≤ 0.1 (the sizes: the harness's `DERIVATION.md`, `02_roots.md`).
- Its risk is the low-information regime, where one or a few pairs make the od a two-point law: LBX0588 at 1 %
  (150–160 spliced reads) reads 0–0.098 shipped and 0–0.069 joint, a gDNA-total cv of 5.2 % shipped and 4.3 % joint
  against a 0.8–2.3 % fixed-value floor, and the od spread predicts 3.3 % of the observed 4.9 %. Run the fallback
  page's §7.3 (withhold all but m ∈ {1, 2, 5} pairs).
- Deleting the clip supersedes DESIGN §3.3a's ceiling ruling, the `test_sweep.py` pin, `config.py`'s constant and the
  CLAMPED flag. Falsifier: with ρ ≤ 0.1 planted and clean seeds, the fit on [0, 1] exceeds 0.2.

RULING 1'S PRICE, watched at every step:
- MO_3021 rejects a common value: at full depth gDNA reads 0.133 (information 19,436) against RNA's 0.008 (honest
  information 14,106), z ≈ 11, and the joint 0.121 comes 87 % from the gDNA seeds. LBX0588 is consistent (z ≈ −0.6).
- Where the truths differ, one od is a weight-dependent compromise: an agreement rule reads ≈ 0 on odg05, giving up
  step 1's gain.
- The gDNA–RNA pair mean is 0 only under independence; a shared local bias would make it positive. Unverified.
- Sharing widens RNA's strand term: first-pass Σ|err| on odg05 stranded OFF 75,873 → 94,819 (+25 %, the oracle +3 %),
  stranded ON +7 % (oracle +4 %), ss 0.70 +71 %; N_eff = N/(1 + (N−1)·od_r) falls 111 → 34 at N = 1,000 as od_r goes
  0.008 → 0.0287. Each landing reports the first-pass cost, the boundary axis (+14.8 %), the gene level at capture ON
  (+473) and the gDNA pool.

THE ROBUSTNESS RULE'S FOUR READINGS (ruling 15 takes the third): the od spread across seeds at one depth (every
estimator fails at 1 %); the od across depth; the gDNA-total spread across seeds (the joint fit passes at 10 %: 1.3 %
against a 0.7–1.5 % floor); the gDNA total across depth (calibration's own drift, unusable). Read literally, "swings"
fails even the planted truth (odg05: gene stranded ON +145, gDNA pool +2,606).

RUNS THE LANDING A/B STILL NEEDS: `policy_benchmark.py --by-class`, `ruler_vs_truth.py` and `prior_vs_oracle.py` under
the arms; AMBIG slots scored apart (at n = 1,000, u = 380 the tilt peak moves 0.73 → 0.70 → 0.85 as od goes 0 → 0.05 →
0.2); the `strand_evidence` consumer priced apart from the row width; a like-for-like fixed-0.05 arm on the ladder
(LBX0588 at 10 % reads 0.47 at ρ = 0 against 0.14 at 0.2); one harness run to give the two-witness refusal a measured
number (its chord prediction: a cv of 2.2 % against a 0.78–1.28 % band). Falsifiers beyond each step's: a
forced-fallback injection puts 0 beyond the oracle's od-nil floor in an in-scope stratum; an fl landing moves an od
value.

Instrument: the harness, substrates and report in `~/Downloads/rigel_runs/prototypes/2026-09-27_strand_od/`; one
derivation page per ruling in `~/Downloads/rigel_runs/prototypes/2026-09-28_od_design/` (`00_synthesis.md` the landing
order; `01_fallback`, `02_roots`, `03_bound`, `04_evidence`, `05_panel` the rulings' pages). Refused designs:
`ISSUES: the-overdispersion-design-refusals`.

### the-realized-gdna-length-law-reads-rna-counts
`priority: now — Tier 0: step 1, a numeric no-op, can land first; step 2 after strand-overdispersion step 3, never in its A/B window (owner, 2026-09-28) · kind: defect · 2026-09-28`
`_realized_gdna_counts` (`fl.py`) normalises the gDNA law but not its RNA input: the caller passes the de-tilted
spliced census (`rna_fl_mass`) as `rna_pmf`, and `contained_opportunity` is linear, so every boundary's RNA rate is
divided by the spliced count N_s, R_b → 0, and every exon-flanking boundary reads as pure gDNA. It reaches the realized
law, the uniform law (through `_couple_estimands`) and the E-step, which scores with the realized law
(`ISSUES: the-scorer-reads-a-census-length-law`). The tests pass a normalised pmf and cannot see it; a fixture reads
a2 = 0.026 / 0.55 / 0.963 / 1.000 at an input scale of 1 / 30 / 10³ / 10⁶.
MEASURED on the ladder's no-EM census, realized mean in bp, shipped → fixed (truth): g05 OFF 219.95 → 216.99 (216.90;
on-target share 0.171 → 0.018); g50 OFF 216.93 → 216.74 (216.67); g05 ON 236.22 → 244.97 (240.97); g50 ON 234.52 →
239.11 (240.84); g98 ON 234.87 → 238.58 (240.87). At g05 ON the error flips −4.7 → +4.0 bp: the bug was partly
cancelling the unmeasurable-region over-pricing, which remains
(`ISSUES: the-fl-boundary-inversion-reads-missing-evidence-as-a-value`). The uniform law at g05 ON stays 3.6 bp high
after the fix (221.76 → 220.46 against 216.90; g50 and g98 close within 0.35 bp): one realization, unattributed —
re-measure after step 2 before giving it a home.
THE LANDING, four steps, each alone (owner, 2026-09-28: the EB pmf, one refresh, step 2 in its own A/B window):
1. Delete the seven guards that never bind, bit-identical (`rename_identity.py --check`): the four `max(·, 1e-30)`
   and `max(μ_g − 1, 1e-9)` floors in the boundary odds, `not np.isfinite(mu_next)`, and the two
   `max(m_C + m_B, 1e-30)` after the early return (a fitted ρ_off ≥ min(Σ_{n>0} E/ΣE², 1/E_max) > 0).
   `max(μ_r − 1, 1e-9)` binds at N_s = 0 and goes with step 2. Add unit fixtures at zero gDNA and at N_s = 0, which no
   frozen reference reaches.
2. The fix: normalise inside the function; pass `FLModels`' EB pmf, built once (today it is built after this
   function, in `_fl_models_from_histograms`); write E_b = ρ_off(μ_g − 1) + ρ_adj(μ_r − 1). The EB pmf carries
   `POOL_EB_PRIOR_ESS` into the boundary odds, and plain normalisation invents a law at N_s = 0
   (`ISSUES: eb-shrinkage-magic-ess`). Gates, each failing today: the output is invariant to c·rna; on an
   off-capture expected-count fixture the on-target share is 0 and the realized law equals the uniform; the array
   passed is `FLModels.rna_pmf`. Re-take the census in the same session (the 2026-09-26 one predates the purity-tie
   fix); judge the off-capture closure, then the capture-ON rows against truth, then `calibration_vs_oracle.py`, then
   `quant_accuracy.py` fractional. Falsifiers: g05 OFF realized − uniform outside +0.3 to +0.6 bp; g50 OFF above its
   +0.04 to +0.05 bp rectification floor by more than the replicate spread; g50 and g98 ON not rising about 4 bp.
   Settle `GdnaContrast` and `GdnaRealized` here — computed and discarded, read only by `test_fl.py`, their docstring
   still calling them a surfaced record: surface them in the diagnostics or call them test surface
   (`ISSUES: hygiene-ledger` (c)).
3. An empty crossing pool is structural. `den` counts every crossing but n2 and n3 only single crossers, so den > 0
   with an empty pool is reachable (n3 = 0 at g00), and `_normalized`'s uniform 0–1000 law then enters r̂ (n2 = 0) or
   becomes g_B and pushes μ_g toward 500 (n2 = n3 = 0). No pool means no boundary stratum; one empty pool skips the
   inversion. Gates, both failing today: den > 0 and n2 = 0 gives g_B = the normalised f3; n2 = n3 = 0 leaves μ_g
   unmoved and m_B = 0. It moves numbers: on the capture-ON ×30 toy the on-target share goes 0.042 → 0.029.
4. One unconditional refresh replaces the 0.25 bp stopping test, as straight-line code with f2, f3, n2, n3 and den out
   of the loop. The test fires 98 / 95 / 75 % of the time at 260 / 2.6k / 26k crossings; the secant gives q ≈ 0.03, so
   one refresh lands within 0.011 bp of the fixed point (g05 ss.99 ON: 244.569 / 244.965 / 244.954 at zero refreshes /
   one / the fixed point), 30× below the 0.32 bp replicate spread. Unpriced: off capture the refresh swaps the
   contained pool's mean for the noisier crossing pool's.
It shares no code with the od work and never shares its A/B window: both move `calibration_vs_oracle.py`.
Instrument: `~/Downloads/rigel_runs/prototypes/2026-09-28_od_design/06_fl.md` and its `scratch/` (census, toys);
`calibration_vs_oracle.py`, `quant_accuracy.py --set em.assignment_mode=fractional`.

### the-gdna-landscape-collapses-at-low-depth
`priority: next — Tier 1, after the strand-overdispersion and fl.py landings, A/B'd apart from them · kind: defect · 2026-09-28`
Calibration's gDNA estimate falls with sequencing depth. LBX0588 (stranded, capture ON) reads 0.11 / 0.45 / 0.83 gDNA
per deposited fragment at 1 % / 10 % / full depth. The EM follows it (gDNA fraction 0.16 / 0.62 / 0.91), so at 1 %
the RNA pool is about 9× too large. The VCaP mix, which has truth by read name, reads P/O 0.125 / 0.52 / 0.92.

The mechanism is the refits' gDNA landscape.
- A slot trains as located only when Var(log f_g) ≤ 1 nat², which in practice means it holds a gDNA fragment.
- At 1 % depth an object at VCaP's bulk density needs about 1.26 Mb to expect one fragment. So 0.23 % of the ~323k
  kernels locate, and those few are selected above the bulk.
- The rest are zero-count anchors, flat down to the grid floor. The refit loop converges to one mode at the floor
  (log ρ ≈ −17) in all three real libraries.
- Under that prior, sparse unspliced fragments go to RNA.

How much each part owns:
- The refits own 80 % / 69 % of LBX0588's drift (ln P_full/P_d, from 1 % / 10 %) and 86 % / 63 % of VCaP's.
- Pass 0 drifts less. It nets an under-call on no-RNA objects against an over-call on no-gDNA objects.
- The capture reference owns none of it, because it is computed after the solve, but it goes missing: its members
  fall to 0 at 1 % for VCaP, LBX0588 and MO_3021 (0 at every depth for MO_3021), which sets every efficiency to 1. On
  the test chromosome at capture ON the members are g001 0 / 0 / 57, g01 0 / 73 / 329, g05 20 / 313 / 352 (d100 /
  d10 / full), and a reference that does form at low depth reads above proportional (×1.97 at g05 d100, ×1.48 at g01
  d10). What that costs the EM is unmeasured, and `calibration_vs_oracle.py` cannot see it. The reference is a yes/no
  switch; its continuous form, an efficiency whose uncertainty grows as members fall, is a direction to price.
- MO_3021 barely drifts (0.82 of full at 1 %): 65–80 % of its gDNA sits on intergenic regions and gene edges, which
  read their count at any depth.
- The level survives in aggregate. VCaP's intergenic regions at 1 % hold 0.0097 of full depth's true gDNA.
- The drift is not the od's: calibration's gDNA share falls with depth under every od arm. LBX0588 per deposited
  fragment at 1 % → 10 % → full, beside shipped's above: binomial 0.113 → 0.471 → 0.835, the 0.2 value 0.098 → 0.145 →
  0.679. The VCaP transcriptome half drifts too (gDNA fraction 0.50–0.63 % at full depth against 0.34–0.41 % at 10 %
  and 1 %, every arm); MO_3021 (0.075 → 0.077 → 0.091) and the VCaP DNA half (0.045 → 0.044 → 0.054) barely move. The
  0.2 value deepens LBX0588's collapse at 10 % (0.145 against 0.450), so the od landing and this repair are A/B'd
  apart, in that order.
- Strand inputs, an untested cause: two inputs fitted before pass 0 move with depth. MO_3021's gDNA od reads
  0.171 ± 0.042 at 1 % against 0.0084 at full depth; LBX0588's strand specificity rests on about 150 spliced reads
  (0.913 / 0.919 / 0.934 at 1 % against 0.937), and its pass 0 follows the same order (0.5465 / 0.5555 / 0.5641).
  MO_3021's pass 0 rises as depth falls (0.199 / 0.164 / 0.139), the opposite of the others. On the test chromosome,
  planting the true strand model moves the error by at most 6.3 points. Check once
  `ISSUES: strand-overdispersion-one-shared-value` step 1 lands.
- The EM does not always follow calibration. On VCaP it recovers much of what calibration loses (EM recall 0.614 /
  0.747 / 0.853 against calibration P/O 0.125 / 0.518 / 0.922); LBX0588's EM follows calibration down (0.164 / 0.620 /
  0.906). Why is untested: the libraries differ on three axes at once — strand specificity 0.9998 against 0.91–0.94,
  spliced share 51 % against 1.3 %, fragment length (VCaP DNA 307 bp and RNA 244 bp against LBX0588's 104 bp).

On the test chromosome's stranded depth family:
- Under a depth-invariant true landscape, the landscape owns 69–94 % of the capture-ON low-gDNA drift at d100 and
  d10.
- The capture-OFF drift is per-object.
- Injecting the full-depth landscape, shifted by log f, cures g01 and g05 at d100 capture ON (−81.1 → −0.8 and
  −52.6 → −8.3 points of O). It only halves g001, does nothing at capture OFF, and worsens g25 and g50 at d100.
- True values on every structural region at depth d still fail at g01 and g001 d100. So fitting count/E kernels on
  sparse counts is itself part of the defect.
- The same fragile regime at full depth: at g001 capture ON, removing the one 19-fragment `test_blank` slot from the
  training set moves the annotated chromosome's net gDNA error +3 → −630 (Σ|Δ| 983.1 → 1,090.9, n_train 754 → 722).
  The depleted mode rests on very few located kernels at low gDNA.
- Under a correct, depth-invariant landscape the capture-ON residual is misplacement, not loss: g01 d100 no-gDNA
  +62.95, mixed −45.91, no-RNA −23.12 (net −6.08); g001 d100 +83.34 / −55.79 / −24.81 (net +2.75); g001 full +20.63,
  cancelled by the landscape term −16.23. A repair is scored BY CLASS (no-gDNA, mixed, no-RNA), never on the net.
- The collapsed landscape beats the true one at g25 d100 ON (shipped −6.8 against −10.9), and the shifted full-depth
  landscape moves g25 −6.8 → −10.6 and g50 −4.2 → −8.0: two compensating errors, and a repair must not regress these
  rows.
- Not separated: whether the landscape or the density fit removes pass 0's mirror floor (RNA-only objects 1,894.5
  against a predicted 1,806.9, and 65.3 after the refits). No arm isolates the message layer: shipped − silent is
  +33.65 (g05 d100 ON), +79.39 (g01 d10 ON), +34.88 (g25 d100 ON).
- The depth family is ss 0.99 only, so the drift on in-scope unstranded × capture OFF is unmeasured; TESTING §0c still
  describes these configs as the ruler's substrate only. RULED (owner, 2026-09-28): add ss 0.50 rows at d100 for g01
  and g05 only.
- An oracle landscape trained on the library's own true counts leaks the answer (it gave the landscape 29–63 % of the
  capture-OFF drift, against about 0 under the depth-invariant true landscape): use the depth-invariant one.

First, alone: a failed refit overwrites the last landscape with None after the belief was solved under it, losing the
reference, the efficiencies and the debug record. It fires in sparse fits.

Direction: at its maximum, the exact Poisson mixture keeps the prior's implied total equal to the pooled count
(Σ E_i·E[ρ | c_i] = Σ c_i). The shipped fit loses that, and a fit at depth f should equal the full-depth fit shifted by
log f wherever the data resolve it.
- The location floor and the anchor E-step that brought the zero controls down
  (`ISSUES: gdna-landscape-trains-on-false-positives`) are what collapse here; a repair must not regress them.
- `estep_all` and `row_kernel` stay refused (`ISSUES: the-landscape-training-population-arms`).

Not this entry: VCaP's 7.8 % full-depth deficit, already present in pass 0
(`ISSUES: vcap-dna-only-objects-lose-gdna-at-full-depth`).

Still unmeasured:
- The real-data arm: the full-depth VCaP landscape, shifted by log f, on the 1 % and 10 % caches; check that it fires
  (`TRAPS: could-the-arm-have-fired`). The VCaP 10 % seeds put the mode at −14.70 / −14.56 / −17.05 with P/O 0.514 /
  0.526 / 0.516 against a bulk at −11.74, so P is insensitive to where the collapsed mode sits. Run it through
  `v2/scripts/cal_nodebug.py` under `one_at_a_time.sh` (4.4–4.8 GiB peak): one real-genome job, never `_debug`.
- What the drift costs the transcript table: `quant_accuracy.py` and `prior_vs_oracle.py` never ran on the depth
  family, nor `prior_vs_oracle.py` on `full_lowg`. The real-library EM numbers are sampled draws (seed 0), and the
  read-name scorer needs hard labels, against the fractional-benchmark rule. If full depth is roughly right, LBX0588's
  RNA pool is 8.9× too large at 1 % and 4.0× at 10 %; VCaP's DNA sent to annotated transcripts is 11.8 / 5.9 / 4.4 %
  of the DNA.
- Full-depth cfRNA has no truth: MO_3021's refits remove 34 % of pass 0's gDNA (0.139 → 0.091 per deposited fragment;
  exons 21,552 → 4,156; exon|exon 11,363 → 367). Only a truth-labelled mix can settle it (the truth mixes of
  `ISSUES: splicing-artifacts`); a √n depth extrapolation fails on VCaP (8.2 / 2.6 against a true 7.0).

Instrument: `~/Downloads/rigel_runs/prototypes/2026-09-28_depth_lowg/`: `pages/explanation_depth.md`; the depth-drift
chain `rerun/scripts/depth_harness.py`; the arms `shifted_landscape.py`; the tables `rerun/scripts/tables.py` (test
chromosome) and `v2/scripts/tab_v2.py` (real libraries); the substrates `test_reference/scenarios_depth_d100`, `_d10`
and full depth. RULED (owner, 2026-09-28): the harness stays outside the tree until a repair is A/B'd; the one tree
addition is the gDNA-only row of `ISSUES: the-background-dispersion-assumes-a-pure-intergenic-pool`.

### gdna-landscape-trains-on-false-positives
`priority: next — Tier 1, a constraint on the landscape repair: its zero controls must not regress · kind: question · 2026-09-02`
What is still open in the landscape estimator, each with its number: (a) under capture a short dark region's
mass is split over both modes of the previous fit (deferred 1.01–1.02×, stranded ON 1.004×); (b) the zero-RNA
controls move ±4 % (`g05 ss.99 ON` +4.4 %, `g98 ss.50 OFF` −1.4 %); (c) the Beta(½,½) reference still decides
a blind slot. CLOSED here 2026-09-14, (d): a slot whose solve is wider than one nat² in `log f_g` no longer
trains (`DESIGN.md` §7.1 rule 4) — the vertex-solved exons that trained 3,771 false fragments at the first
refit on `g00 ss.99 OFF` are out, and the four ladder zero controls read 282 / 194 / 265 / 172. Refused here:
excluding κ-dead exons (`g50 ss.50 ON` 2,691 → 56,422), AMBIG in the final fit (worse 25/32), and the two
other readings of the floor (`DESIGN.md` §7.1 rule 4). Its instrument, `landscape_training_census.py`, was retired 2026-09-14 (in git).

### the-background-dispersion-assumes-a-pure-intergenic-pool
`priority: next — Tier 1, its own A/B after the landscape fit · kind: defect · 2026-09-28`
`fit_intron_background` (`density_deconv.py`) fits the gDNA background's mean and its dispersion α by moments on the
intergenic regions. It assumes they hold only gDNA, and its own comment says unannotated transcription breaks that
(`TRAPS: purity-is-a-property-of-the-annotation`).

On the test chromosome's `full_lowg` (stranded, capture OFF):
- The unannotated transcription on `test_blank` raises the fitted mean 12.7×. It drops α from ∞ (no extra spread) to
  0.142 at g001 and 0.405 at g01.
- The lower α widens the prior on intron gDNA. Changed alone, it moves the annotated chromosome's gDNA by −197 / −687
  fragments; the shipped figures are −189 / −643.
- It owns 27 % / 38 % of the annotated chromosome's gDNA under-call (51 of 192 and 235 of 624 fragments).

Real intergenic regions hold strand-coherent RNA too, about 10–11 % on VCaP, so the route is live on real data. Its
size there is unmeasured.

Read it with the gDNA-only objects apart:
- On `full_lowg` capture OFF, `test_blank` itself carries 8,844 of the 9,194 Σ|Δ| at g001: the designed control,
  locked to gDNA by structure, not a defect. `calibration_vs_oracle.py` reads 2.31× = (9,825.1 + 20,009.7) / (1,170 +
  11,764), driven by it; without it the reading is 0.92× (0.923× shipped, 0.950× with a gDNA-only background). Every
  od arm reads the same 2.31×; the 0.2 value lifts the capture-ON P/O 1.095 → 1.307.
- The annotated remainder still carries the shadow, through the background's α, so a clean reading needs the gDNA-only
  background fit beside it.
- `calibration_vs_oracle.py` does not report structurally gDNA-only objects on their own row. RULED (owner,
  2026-09-28): add the row, the only tree addition for Tier 1. Until it lands, the debug loop must not pick
  `full_lowg` capture OFF as the worst in-scope case on this number.

A repair should read the half of the statistic the contaminant cannot reach, as its own A/B after the landscape's.
Whether it needs a new constant is for the derivation.
- fl's one-sided off-target rate has the same "intergenic = gDNA" contamination in a different estimator
  (`ISSUES: the-fl-boundary-inversion-reads-missing-evidence-as-a-value`).
- Under capture, intron regions next to probes are not off target either
  (`ISSUES: intron-seeds-near-probes-are-capture-enriched`).

Instrument: `~/Downloads/rigel_runs/prototypes/2026-09-28_depth_lowg/review_alt/bg_split.py`. It fits the background
on the gDNA fragments alone, beside the shipped fit. See also `pages/explanation_lowg.md`.

### capture-on-overcalls-gdna-at-low-gdna
`priority: next — Tier 1, read after the landscape repair · kind: defect · 2026-09-28`
On the test chromosome's `full_lowg` with the `test_blank` control excluded, capture ON is about twice as wrong as
capture OFF, and over-calls: annotated-chromosome Σ|Δ| 983.1 ON against 349.5 OFF at g001 and 2,846.3 against 1,506.7
at g01; the net is positive ON and negative OFF (`blank_arm.txt` +1,089 / −832; a second table +1,193 / −816). It
looks fine overall because capture removes the unprobed shadow. The mechanism is the strand resolution limit: capture
concentrates the low gDNA onto probed exons holding about 1,500 RNA fragments each (at g001, exon true gDNA 70 → 717
under capture), where the summed 0.42·√(n/I) is 1,135 against 1,107 true gDNA, and it errs both ways (g01 ON
single-strand exons +406). The located part is the AMBIG exons at g01 ON — 84 slots holding 524 true gDNA fragments,
read as 1,517 — at +993 shipped, +939 under the oracle landscape and background, +1,151 silent, +2,870 at pass 0. At
g001 the silent policy reads AMBIG +8,398 (capture OFF, P/O 15.58) and +5,029 (ON) against shipped −3 / +51: only
transfer's final solve removes it. The capture-OFF AMBIG over-call is +9 to +22. Untested candidate:
`ISSUES: the-atom-at-an-unwitnessed-both-strand-slot`. Whether any in-scope mechanism beats the resolution limit is
open.
Instrument: `~/Downloads/rigel_runs/prototypes/2026-09-28_depth_lowg/`: `results/blank_arm.txt`,
`review_alt/arm_classes.txt`, `results/floor_test.txt`.

### strand-plug-in-bias-on-sparse-libraries
`priority: next — Tier 1: the shallow-object loss below, not the plug-in alone, owns the capture-OFF depth drift; invisible on the ladder, no instrument in the tree · kind: defect · 2026-09-21`
Calibration's per-object gDNA fraction is a plug-in that conditions on the object's own fragments, this one
included; the exact statement is the leave-one-out posterior, and the plug-in understates gDNA by −21 / −22 /
−19 / −11 / −2 % at 1 / 2 / 5 / 10 / 50 fragments. The pool-level bias is (live objects) / (unspliced incidence):
0.2–0.4 % on the ladder, where it is guaranteed invisible, and MEASURED 2026-09-21 on real data 0.216 / 0.120 /
0.308 on the cfRNA libraries LBX0190 / LBX0588 / MO_3021 (72–91 % of their live objects hold five or fewer
fragments and carry 15–44 % of the incidence) and 0.055 on the deep VCaP library. Any per-fragment use of the
fraction must take the leave-one-out posterior and must not apply the fragment's own strand term a second time.
The plug-in is not the whole cause: shallow gDNA-only objects lose far more gDNA than it or the closed zero-RNA floor
predicts, even under the oracle landscape and background. On `full_lowg` g001 capture OFF 214 objects (O 304) lose
193 shipped, 155 with a gDNA-only background and 149 under the oracle landscape and background, against bounds of 107
(strand only) and 52–70 (the closed floor); at g01 OFF 556 objects (O 2,158) lose 502 / 338 / 301 against a bound of
429. The ladder's unstranded rows, where no strand plug-in exists, show the same staircase: shallow objects understate
gDNA by 15–40 % (real / predicted 1.04–1.40 at 1–2 fragments, 1.00 at depth). The depth family's [1,3) bin loses
49–99 % at d100 against 4–45 % at full depth; VCaP at 1 % reads its pass-0 introns at 169 against 1,915. The candidate
is the median read-out of a skewed shallow posterior; splitting it from the plug-in needs the kernel's posterior. This
loss owns the capture-OFF depth drift of `ISSUES: the-gdna-landscape-collapses-at-low-depth`.

### nascent-stress-sensitivity
`priority: later — Tier 2, first: a cheap re-measure of Tier 0's verdicts at the realistic nascent level · kind: question · 2026-08-22`
Does any in-scope verdict depend on the nascent stress level? The ladder runs `on_fraction 0.50`; realistic is ~0.10
(`DESIGN.md` §0b). MEASURED 2026-09-20 for the siphon: `g50 ss.99 ON` re-simulated at 0.10
(`~/Downloads/rigel_runs/suite/ladder_nrna_lo`) reads a siphon of +522,205 against +541,216 at 0.50, so that verdict
does not depend on the stress level — the repair is worth the full amount in the expected case. Re-simulate the worst
in-scope scenario at the realistic level and check whether any rank moves; a verdict that holds only at stress is a
robustness finding. The verdict to re-measure: at the realistic level (`g50 ss.99 ON`, 2026-09-21) the synthetic
nascent pool read 4.2× its truth — baseline 5.50 (nascent 4.2×), pseudocount-corrected 4.98 (6.3×), oracle 1.85 (2.2×)
— before the one shared length rule and the count-form pseudocount. `sim/panel.py`, `quant_accuracy.py`,
`policy_benchmark.py`.

### psi-reads-kappa-where-the-strand-channel-is-dead
`priority: later — Tier 2, the first mechanism: in scope, unstranded × capture OFF · kind: defect · 2026-09-28`
ψ reads the estimated κ̂ where the protocol's Bayes factor calls the strand channel dead: the discriminability gates
only `strand_evidence`, and `solve_kernel.cpp` builds ψ's grid with κ̂ unconditionally. At od = 0, κ̂'s sampling error
gives the strand term a slope n(½ − κ̂), so real overdispersion reads as composition along it, at a weight of about
n/m for m spliced reads (SD(κ̂) ≈ 0.08 at about 40 spliced reads). On odg05 ss 0.50 g98 (true od 0.05, κ̂ = 0.515),
assuming od 0 costs +1,178 in scope at capture OFF and +29,106 on the deferred capture-ON row, and nothing against a
true 0 (test panel, κ̂ 0.485). Not dissected. The candidate — ψ reads κ = ½ where the channel is dead — is its own
A/B, beside the od step-1 skip of the fit there (`ISSUES: strand-overdispersion-one-shared-value`).
Instrument: `~/Downloads/rigel_runs/prototypes/2026-09-28_od_design/01_fallback.md` (its last open question);
`policy_benchmark.py --by-class`, unstranded rows.

### the-scorer-reads-a-census-length-law
`priority: later — Tier 2, re-measured after fl step 2 (ISSUES: the-realized-gdna-length-law-reads-rna-counts) · kind: defect · 2026-09-24`
The E-step scores an unspliced fragment's length with gDNA's library census (`gdna_realized_pmf`,
`pipeline.py`) against RNA's spliced census. Under capture the length selection depends on where a fragment sits,
and every hypothesis at one footprint shares that footprint's capture — so the ratio wants the two origins' laws in
ONE frame, and two censuses that average capture over different positions tilt short exon-contained gDNA fragments
toward RNA even when both are estimated exactly. Confirmed free of the count and length confounds (g50 ss.99 ON, exact
count, one-scale lengths): the same law for both origins gives annotated −2.7k, spans 156.1k (truth 150.4k),
transcripts 1.70 %; the true census gives +12.7k, 139.2k, 1.79 %; shipped +5.2k, 158.6k, 1.75 %. At `g98 ss.99 ON`
(oracle arm) the scorer's law is worth 22.2k gDNA and 9.8k transcripts, 18.4k / 4.6k of it the estimator
(`ISSUES: the-gdna-length-law-falls-back-at-identical-purities`). Nil off capture. The same law for both is a
diagnostic, not a candidate: real libraries with different gDNA and RNA chemistry need the channel. The repair is a
derivation of the per-fragment length term at a shared footprint. 2026-09-26: the estimator's half is closed
(`ISSUES: the-gdna-length-law-falls-back-at-identical-purities`); at `g98 ss.99 ON` the realized law now reads
234.9 against a true 240.9, the same offset as every other capture-ON row. The realized law it reads is itself fed
RNA counts (`ISSUES: the-realized-gdna-length-law-reads-rna-counts`): re-measure after that fix, and land
`ISSUES: the-pooled-q-in-the-gdna-count` with the repair.

### the-pooled-q-in-the-gdna-count
`priority: later — Tier 2, lands with ISSUES: the-scorer-reads-a-census-length-law · kind: defect · 2026-09-24`
Calibration's per-locus gDNA count (`priors.assemble_priors`) converts each boundary's gDNA mass to fragments by
`boundary_mass_per_crossing`, pooled over gDNA and RNA, so where RNA dominates a boundary gDNA is over-stated even
under a perfect calibration (`prior_vs_oracle.py` O − S: +33.7k at `g50 ss.99 ON`, +2.9k at g98 ON, +1.9k off
capture). The EM passes a g50 count change through at only 0.37–0.50, so its size in the result is about 12k: with
gDNA's own share the g50 gDNA pool moves −18.2k → −30.4k, spans 142.0k → 152.1k (truth 150.4k), transcripts
unchanged. It cancels part of the capture likelihood's lean toward RNA today, so it lands with those repairs. The
gDNA LENGTH already uses gDNA's own share (`ISSUES: the-pooled-q-in-the-gdna-length`, replaced for the length only).
`prior_vs_oracle.py` (P, O, S) is its one instrument and retires when it lands (owner, 2026-09-28).

### the-fl-boundary-inversion-reads-missing-evidence-as-a-value
`priority: later — Tier 2, the fl second wave, after the realized-law fix · kind: defect · 2026-09-28`
Missing evidence must drop out, never be replaced by an invented value; after
`ISSUES: the-realized-gdna-length-law-reads-rna-counts`, four places in `fl.py` still break this.
(a) UNMEASURABLE INTRON REGIONS. A region too short to hold a fragment has e_r ≈ 0, so its RNA density reads 0 rather
than unknown, a_b → 1, and nascent exon–intron–exon fragments are priced as capture. At g05 OFF after the fix the
signed w(ε − 1) reads +0.66 against a rectified +0.70 (bias); at g50 −0.006 against +0.072 (rectification only).
Exon|intron pairs with ρ_adj = 0 hold 10 % of boundary counts at a count-weighted ε ≈ 20–21.6; pairs with ε > 5 border
introns of median 94 bp (31–185). The likely cause of the +4 bp g05 ON overshoot. OWNER: carry the variance (the
weight goes to 0) or borrow from the gene's measurable introns? One derivation could also serve the purity weight's
empty-neighbour failure (`ISSUES: strand-overdispersion-one-shared-value`, step 5) and
`ISSUES: the-efficiency-posterior-floor-on-empty-pieces`.
(b) ZERO gDNA. The one-sided rate declines only at total_n = 0, so at g00 ρ_off is a geometric floor (2.1e-6, against
5.3e-3 at g05): ε reads about 2,300·n_b wherever ρ_adj = 0, and the realized law turns 99 % "capture excess" (+28.9 bp
over uniform after the fix, +19.1 before: the fix makes it worse). `if den[0] <= 0 or den[1] <= 0: break` drops the
whole boundary stratum at g00 (den = (185k, 0)), since simulated RNA never crosses an exon|intergenic edge. The
one-class limit would admit a2·n2 RNA crossings against m_C = 80; derive it once for both den and n. OWNER: how should
ε fade? It is a floor, not noise, so a variance term may not cure it.
(c) THE CONTAINED PAIR still declines on an empty pool in the four g00 rows, falling back to the capture-selected
four-pool census (243–244 bp under capture; at g00 ON the realized law reads 263 bp against 237, RNA 227). Left open
by `ISSUES: the-gdna-length-law-falls-back-at-identical-purities`.
(d) INTERGENIC RNA. The off-target rate is fitted on intergenic and intron counts, so RNA spread over them lifts it:
toy fixed point 202 bp against 217 (a3 0.63 against a true 0.37); noise-free off-capture steps −10.0 / −15.2 bp at
ig_rna 0.5 / 2. The fl counterpart of `ISSUES: the-background-dispersion-assumes-a-pure-intergenic-pool`, in another
estimator (`ISSUES: capture-degeneracy-standing-risk`). Ladder RNA never reaches intergenic space, so its real-data
size is unmeasured.
Instrument: `06_fl.md` and `scratch/fl_loop2.py` under `~/Downloads/rigel_runs/prototypes/2026-09-28_od_design/`.

### the-fl-boundary-inversion-has-underived-pieces
`priority: later — Tier 2, the fl second wave, with the missing-evidence entry · kind: derivation · 2026-09-28`
Four pieces of `fl.py`'s boundary inversion are borrowed, not derived.
(a) The one-sided excess clip (cnt − ρ_off·e_g)₊ makes a pure-gDNA region's expected share read below 1 at every
finite depth: E[a_b] = 0.914 / 0.779 / 0.853 / 0.910 / 0.965 at λ = 0.1 / 0.5 / 1 / 10 / 100; real introns average
0.910 (LBX0588), 0.940 (MO_3021), 0.978 (LBX0588 at 10 %). a2, a3 and m_B carry the bias, and the flatness premise
a3 = 1 never holds on a real library; in the toy it moves the realized mean by at most 0.03 bp. Gate: a Poisson
pure-gDNA fixture at λ ≈ 0.5–1 bounding a2 against the exact E[a_b].
(b) The boundary pair's standard error, se_b = √(a2²/n2 + a3²/n3), is borrowed from the contained pair and was never
derived for the incidence-weighted a2 = Σa_b·n_b / Σn_b; it sets the inversion's resolution weight.
(c) The boundary's RNA side is priced at the spliced census mean, not at its own RNA length law; under a capture that
prefers some lengths E_g[s]/E_r[s] need not cancel. Unmeasured, and distinct from
`ISSUES: the-scorer-reads-a-census-length-law`.
(d) UNCONFIRMED: the inversion mixes by incidence shares (n2, n3) while the de-tilted pools mix by start density. A
first draft read the fixed point's |Φ′| falling 0.026 → 0.005 on switching, and no results file exists. Measure first.
Instrument: `06_fl.md` under `~/Downloads/rigel_runs/prototypes/2026-09-28_od_design/`.

### capture-blind-gdna-divisor
`priority: later — Tier 2, the fl second wave, after fl step 2 · kind: defect · 2026-08-31`
`gdna_opportunity_from_index` is computed from the index alone, so under capture it removes ~6 bp of a ~30 bp
length selection — the gDNA control moved +6.0 % on all six capture-ON rows (gDNA has no introns to miss), and with
`ISSUES: eb-shrinkage-magic-ess` it owns the −5.90 % capture-ON length ceiling. `capture_eff_length` already models
the panel; it also blocks `ISSUES: crossing-pool-contrast`. ⛔ Not priced on the deliverable: no `quant_accuracy`
arm reaches the locus gDNA length (the oracle overrides only calibration's count arrays; the `oracle_efflen` reading
is withdrawn in `ISSUES: end-to-end-error-unattributed`). A price needs the simulator's per-locus gDNA opportunity
as the truth.

### eb-shrinkage-magic-ess
`priority: later — Tier 2, the fl second wave: replaced on the re-simulated fl-gap panels · kind: defect · 2026-08-31`
`POOL_EB_PRIOR_ESS = 1000.0` shrinks the gDNA pmf toward `global_pmf` (mostly RNA whenever gDNA is a minority) at a
magic ESS: inert on the ladder (0.01 bp), dominant on the fl-gap arm at `g05` capture-ON (`ship−pool` −23.7 of −31.7
bp). Replacement: reconcile the pools by their precision (`EQUATIONS.md` §6c). Its instrument, `fl_pool_purity.py`,
was retired 2026-09-14 (in git). The ladder's deliverable cannot see it; only the fl-gap arm can. The realized-law fix
(`ISSUES: the-realized-gdna-length-law-reads-rna-counts`) passes the EB-smoothed RNA pmf (owner, 2026-09-28), so the
ESS also reaches the boundary composition, shifting μ_r at small N_s by ess/(N_s + ess)·(μ_anchor − μ_RNA). Plain
normalisation is no alternative: at N_s = 0 it invents a uniform law (a2 0.109 against a truth of 0.143). In the toy
a2 moves only 0.1472 → 0.1482 as N_s goes 0 → 1e6, and on the ladder EB and plain agree to 0.002 bp; N_s ≈ 1e2–1e4 is
unmeasured. The replacement is measured on the re-simulated fl-gap panels
(`ISSUES: flgap-panels-stale-nascent-model`).

### intron-seeds-near-probes-are-capture-enriched
`priority: later — Tier 2, the od research tail; answered before any purity weight or background repair that treats intron regions as off target · kind: question · 2026-09-28`
Probes overhang exon edges and bind the bases there, so the intron regions next to probes are capture-enriched and
"every region seed is off target under capture" is false. On odg05 g05 capture ON truly pure intron seeds read
n/(λ_off·E) at a median 1.89, 90th percentile 8.2, 99th 44 (1.00 / 1.23 / 1.84 at capture OFF): a pooled 1.9× λ_off.
Three places treat introns as off target: `density_deconv`'s background (introns at the intergenic depletion); fl's
one-sided rate (intergenic and intron regions together); the proposed purity weight
(`ISSUES: strand-overdispersion-one-shared-value`, step 5). Proposal: measure enrichment against distance to the
nearest probe, then treat near-probe introns as enriched or price the loss.
Instrument: `~/Downloads/rigel_runs/prototypes/2026-09-28_od_design/04_evidence.md`.

### strand-likelihood-over-confident-beyond-od
`priority: later — Tier 2, the od research tail; re-read after strand-overdispersion step 1 · kind: question · 2026-09-28`
A nonzero RNA overdispersion beats the simulator's planted RNA od of 0 on several strata, and the 0.2 value lowers
false gDNA on the zero controls, so the RNA strand likelihood looks over-confident even where its od is truly 0. Its
variance is frozen per slot at the reference composition (`psi_kernel.h`'s slot cube, `transfer_rows.h`): a
misspecification in the strand term, outside the od estimator. MEASURED (the od harness, 2026-09-27): on odg05
stranded ON the joint fit beats the planted pair on calibration (172,329 against 179,181) and transcripts (−28,829
against −11,248); stranded OFF the 0.2 value beats it on calibration (50,971 against 52,004); a planted RNA od of 0
costs +57 % at g98 ss 0.70 ON (184,431 against 117,753); on the ladder's g00 the 0.2 value cuts false gDNA 1,313 →
1,185 (−9.8 %). No mechanism identified; re-read after the joint fit lands.
Instrument: `~/Downloads/rigel_runs/prototypes/2026-09-27_strand_od/`.

### the-pseudocount-strength-is-not-derived
`priority: later — Tier 2, the exact VB update is a ready A/B; the odds are fixed (DESIGN.md §3.1c) · kind: open question · 2026-09-24`
The EM's two pseudocounts carry `P_g + P_R = U`, one pseudo-fragment per unit with a gDNA candidate — the total
the replaced rule carried, kept so the count form changed only the odds. Nothing derives it. The rule that
counts each fragment once is refused as buildable (`ISSUES: a-count-once-density-prior-for-the-strength`): under
capture calibration holds no estimate from outside the locus. Under MAP the strength cancels at the neutral
point; under the shipped VBEM it does not, and the grouped update lets it move the split WITHIN RNA — a
6-fragment isoform reads 4.0 / 5.0 / 5.8 at `P_R` = 0 / S / 10S. The exact VB update for the same prior holds it
at 4.0 at every strength and needs no constant; it is its own A/B
(`~/Downloads/rigel_runs/prototypes/2026-09-24_spliced_rule/toy_vb_nested.py`). Why a prior is needed at all:
with none the EM's gDNA error is −12.1k / −119.7k / −731.8k at `g05` / `g50` / `g98 ss.99 ON` and −0.4k off
capture, and −35.5k remains at `g98 ss.99 ON` under an exact calibration — a lean of the capture likelihood
toward RNA that no strength may be chosen to cover. The owner's framing (2026-09-24): calibration gives the
locus's gDNA FRACTION, and a multiplier converts it to pseudocounts; that multiplier depends on the number of gDNA
candidates and their likelihoods. Today it is the count of units whose gDNA candidate survives pruning (`U`), and
the pruning (`DESIGN.md` §3.1d) is what stops a unit with a vanishing gDNA reading from counting.

### per-transcript-prior-lane
`priority: later — Tier 2, research: the largest in-scope lever, no candidate mechanism yet; the owner's problem after the junction price (2026-09-22, 2026-09-23) · kind: build · 2026-08-31`
`rna_prior_weight` is PLUMBED end to end — `pipeline.py` → `estimator.py` → `em_solver.cpp` — and ⛔ NOTHING IN
`src/` FILLS IT, so the solver takes the fallback `w_i = raw[i]`, which echoes the EM's own belief and cannot
contradict it. ⚠ The lane is ONE static array, so filling it reallocates the WHOLE RNA pseudocount (2,793,710
fragments at `g50 ss.99 ON`, as large as the unspliced RNA itself), and a coverage-derived weight is a far worse
isoform allocator than the EM's own likelihood (5.21 → 53.37 %); a weight that keeps `raw[i]` where the data cannot
speak needs `raw[i]`, which lives in the kernel.
THE CEILING (`quant_accuracy.py --arm oracle_alloc_seed`, the true weights mature and nascent, the committed tree
`c52c9b93`, fractional, transcript Σ|Δ| as a share of the true RNA at `g00` / `g05` / `g50` / `g98`): stranded OFF
2.01 / 1.59 / 2.49 / 15.70 → 0.43 / 0.43 / 0.70 / 6.38 %, unstranded OFF 1.73 / 1.93 / 2.55 / 18.35 → 0.42 / 0.45 /
0.81 / 8.35 %, stranded ON 6.53 / 3.49 / 5.50 / 31.03 → 2.91 / 1.08 / 2.07 / 18.00 %; the deferred stratum 7.32 /
11.00 / 10.11 / 90.63 → 3.41 / 2.67 / 4.03 / 221.15 %. The largest lever on every in-scope stratum. A capability
proof, never headroom: it hands over the true support, and a zero weight is absorbing.
WHERE ITS GAIN SITS (2026-09-20, the four `g05`/`g50` `ss.99` rows): not in support — an exclusive-junction support
rule recovers ~0 % of it (`ISSUES: exclusive-evidence-support-rule`, refused) and the TRUE support only 3–28 % —
but in the PROPORTIONS among expressed isoforms that share every object: 49–61 % of the gain is on isoforms with no
exclusive sj, region or boundary (94 % of them mosaics), 65–76 % in genes with 11+ isoforms. Two weightings are
refused (`ISSUES: refused-transcript-weights`, `ISSUES: refused-soft-min-path-weighting`). The next candidate moves
proportions among mosaic isoforms or is a sparsity argument on silent mosaics. On the one shared length rule the
synthetic pool's error at `g50 ss.99 ON` is the pseudocount's
(`ISSUES: the-pseudocount-prior-is-biased-toward-gdna`), not a siphon's residual
(`ISSUES: nascent-siphons-gdna-under-capture`).
RULED (owner, 2026-09-28): the lane is kept, and with it `warm_start='prior'`, plumbed alongside it, whose only
producer is `quant_accuracy.py`'s oracle-allocation arm; neither is dead code for a sweep to delete.
Instrument: `quant_accuracy.py`, per stratum.

### message-layer-open-cases
`priority: later — Tier 2, the message layer on unstranded capture OFF · kind: question · 2026-09-09`
Four residual cases, none a hole (`DESIGN.md` §6b.4–§6b.14): (a) an exon with both faces speaking — the
arrivals are summed and the price where both carry the same intron's claim is unmeasured; (b) the intron
factory on a region with both exon and intron bits, under capture; (c) the chain of termini — an empty
outside piece (median 12 bp) whose far face is another terminus, half the ladder's terminus-boundary error,
every upper side refused (`ISSUES: the-edge-upper-side`); (d) substrate `nest` (the prior already serves it,
0.666 vs 0.630) and the antisense's nascent variant (`docs/TESTING.md` §0a; `div` was built 2026-09-14).
`policy_benchmark.py --by-class`.

### refit-vs-message-arbitration
`priority: later — Tier 2, the message layer on unstranded capture OFF · kind: design · 2026-08`
At the unstranded × capture-OFF exon cell the refitted gDNA prior and the message both impute one slot with
nothing arbitrating them; the message is the accurate voice there and the refit displaces it. Re-read under the
E-step (with an instrument since retired): the prior does the unstranded rows and the messages the stranded
capture-ON ones. `policy_benchmark.py --by-class`, `calibration_vs_oracle.py`. Belongs with `ISSUES: gdna-landscape-trains-on-false-positives`.

### the-intron-own-solve-on-unstranded-capture-off
`priority: later — Tier 2, the message layer on unstranded capture OFF (the plan of 2026-09-28; parked 2026-09-05); re-measure before ranking · kind: question · 2026-09-28`
On unstranded data off capture an intron region's own solve has no strand channel, so its gDNA fraction is whatever
the landscape prior and the message layer deliver. Measured 2026-09-05 at 43–45 % of the in-scope off-capture error
under every policy; that predates the location floor, the E-step and the one shared rule, so re-measure before
ranking. Related: `ISSUES: refit-vs-message-arbitration`, `ISSUES: two-sided-exon-row`.
Instrument: `policy_benchmark.py --panel ladder --by-class`, the unstranded capture-OFF rows, intron class.

### the-efficiency-posterior-floor-on-empty-pieces
`priority: later — Tier 2: the unprobed class's scale at low gDNA, its EM cost unmeasured · kind: defect · 2026-09-23`
An object's capture efficiency is the posterior mean of its clipped gDNA density under the landscape, from its own
count (`capture_efficiency`). On an object whose support holds far less than one expected gDNA fragment the count
cannot tell a depleted object from a partly captured one, so the posterior reads above the truth. It is what lifts
the unprobed multi-exon transcripts at `g05 ss.99 ON` to L / Y 2.84 (sd of log 0.69, 706 transcripts) against the
synthetic spans' 1.148 and the gDNA footprints' 1.135 in the same probed class (1.56 against 1.005 / 1.013 before
the one shared rule; at `g50` 1.214 against 1.104 / 1.096). MEASURED 2026-09-23: it is the transcripts' PIECES, not
the junction sum — with the pieces at the truth the class reads 1.05, with the junctions at the truth it is
unchanged. The unprobed transcripts' exon pieces are short (median contained support 4.2 starts, ~0.001 expected
gDNA fragments each at `g05`) and read 2.57× their truth from the simulator's noise-free counts and 4.87× from
calibration's (support-weighted), where the spans' long pieces (support ~97) read 1.07× / 1.05×; at `g50` the same
exon pieces read 1.16× / 1.18×. The short-intron fallback of `ISSUES: the-junction-price-is-noisy-within-a-gene` is
the same posterior on the same kind of object. What it costs through the EM is unmeasured: the class competes with
the gDNA component in unprobed genes, where every fragment is off target. A floor or a threshold on support is a
constant; the open question is what an object holding under one expected fragment should read — the question of
`ISSUES: the-gdna-landscape-collapses-at-low-depth` on capture-ON transcripts, and of fl's unmeasurable regions
(`ISSUES: the-fl-boundary-inversion-reads-missing-evidence-as-a-value`). `ruler_vs_truth.py
--scale` (the unprobed rows); the census in `~/Downloads/rigel_runs/prototypes/2026-09-23_junction_price/`.

### scoring-penalties-are-underived-constants
`priority: later — Tier 3, with the splice work's h-weighted reading on the cluster · kind: question · 2026-09-28`
Two fragment-score penalties are 0.1 and underived (`scoring.py`): log(`mismatch_alpha`) per NM mismatch, and
log(`overhang_alpha`), RNA only. Both also set where the pruning lattice cuts. The splice-artifact work adds more to
derive or put to the owner: the minimum count of 2, the anchor floors 4 / 6 / 11, `alignSJDBoverhangMin` 8, the 5-bp
mismatch window. The ε of the proposed h-weighted reading (`ISSUES: splicing-artifacts`) must come from the library,
never from `mismatch_alpha`. A zero alpha gives NaN today (`ISSUES: latent-defects`).
Instrument: none yet; the splice truth mixes.

### multimapper-intergenic-alignments
`priority: later — Tier 3, on the aligned ladder after the splice phases; a separate feature after the end-to-end work (owner, 2026-09-24) · kind: defect · 2026-09-24`
A multimapper whose alignments all lie outside genes is intergenic and counted as gDNA — correct. One with an
alignment in a gene is solved by the EM, and the owner's rule is that it is compatible with gDNA if ANY alignment
is unspliced. The scanner never buffers a multimapper's alignment that overlaps no transcript
(`bam_scanner.cpp:1808-1845`: only unique mappers are appended), so a gDNA read from an unannotated processed
pseudogene whose other alignment splices across the parent gene reaches the EM with its spliced alignment alone
and is certified RNA on the parent. The owner expects multimappers to matter most. The repair buffers those
alignments, lets the scorer take a gDNA term from an alignment with no transcript candidate
(`scoring.cpp`'s `n_cand > 0`) while an all-intergenic group still drops silently, keeps the fragment ledger
exact, and decides which locus's gDNA component owns a gDNA reading at an intergenic position, outside that
locus's gDNA length. Unmeasurable on the ladder (`NH` 1); the aligned ladder
(`~/Downloads/rigel_runs/aligned_ladder/`) and `tests/scenarios_aligned/test_multimap_counting.py` are where it is
seen. Taken with `ISSUES: multimapper-blind-support`. Also open for multimappers: the group's gDNA term is the MEAN
over its eligible alignments (`scoring.cpp`, pinned by `tests/test_pipeline_routing.py`), where uniform gDNA
placement derives the sum.

### multimapper-blind-support
`priority: later — Tier 3, a deferred session with ISSUES: multimapper-intergenic-alignments (owner, 2026-09-24) · kind: defect · 2026-09-16`
The accumulator drops every fragment with `NH > 1`, so an object whose gDNA fragments are multimappers holds no
count the calibration can see, and its capture efficiency — read from its own gDNA count against the fully captured
level — reads it at the DEPLETED level: the transcript is contracted as if unprobed. Measured read-only on the two
captured real libraries on the per-base rule (the session's `multimapper_check.py`: per transcript with ≥ 20
fragments over its pieces, the share of multimapping reads among them from the BAM's `NH` tag, then the length
factor's quantiles per share bin): LBX0588 (reference 2.977e-01/bp, 11,579 kernels) median factor 0.076 / 0.003 /
0.013 / 0.001 at a share below 1 % / 1–10 % / 10–50 % / ≥ 50 % (n 36,446 / 6,161 / 1,429 / 928), the share below
1/10 rising 55 → 90 %; VCaP (8.593e-02/bp, 24,496 kernels) 0.231 / 0.149 / 0.028 / 0.002 (n 280,227 / 35,494 /
2,237 / 1,103), the share below 1/10 36 → 95 %. A hundredfold gap between the unique-mapping and the repeat-rich
transcripts, monotone in the share, on both libraries; the two libraries without a reference read every transcript
at 1 and are uninformative. The confound is real and unmeasured here — a repeat-rich transcript may also be
unprobed by design — so the number is an upper bound on the blindness, not its size; the floor `C/(C+1)` the repair
deleted had hidden it behind a +3.4 nat bias on every unprobed transcript
(`ISSUES: ruler-multimapper-floor-caps-the-correction`, CLOSED). The repair is in the support, not the posterior:
an object's support (a region's contained, a boundary's crossing) should count only the starts a fragment could be
uniquely placed at — a mappability of the support from the index, so that a wholly repeated object has no support,
no evidence, and reads the population's clipped mean rather than the depleted level. Instrument:
`ruler_vs_truth.py` cannot see it (the simulator's reads are unique); the truth on a real library is the
multimapper share against the factor as above, and the repair's gate is that the four bins read alike.
RULED (owner, 2026-09-24): "add multimapping deposit to the accumulator and let calibration include multimapping fragments, which is correct behavior and something we have intended to add".

### unwitnessed-loci-and-multimappers-at-the-em
`priority: later — Tier 3, on the aligned ladder with the multimapper pair; the ladder cannot show either case, and no instrument is in the tree (the 2026-09-21 scratch census) · kind: defect · 2026-09-21`
MEASURED 2026-09-21: every real library carries MultiLoci with a gDNA candidate and no witnessed gDNA (`P_g = 0`:
257 / 237 / 288 / 194 loci holding 374 / 723 / 941 / 2,182 gDNA-bearing units on LBX0190 / LBX0588 / MO_3021 /
VCaP; `L_g = 0` nowhere), and a multimapper share of the gDNA-bearing units of 13.3 / 4.0 / 13.2 / 1.0 %,
concentrated in the sparse loci (per-locus median 0, q90 1–5 % over loci with ≥ 20 such units). A kernel that
pins gDNA at calibration's count sends every unspliced unit of an unwitnessed locus to RNA and under-calls
gDNA by the multimapper share one for one; the ladder has neither case (0 of 1,213 loci). A per-fragment prior
prices each multimapper placement at its own objects, which is the cleaner rule.

### sj-strand-tag-chosen-from-the-first-reads
`priority: later — Tier 3, local, small and independent of the splice phases · kind: defect · 2026-09-28`
The scanner picks the junction-strand tag mode from the first 1,000 spliced reads (`detect_sj_strand_tag_native`) and
returns "none" if none carries XS or ts, so a library whose first reads lack the tag loses its strand. An unknown
`--sj-strand-tag` spec silently becomes XS_TS (`parse_sj_tag_spec`'s default), and the resolved mode is never logged
or written to `summary.json` or `config.yaml`, so a run cannot be reproduced from its outputs. The repair reads the tag
per record; it moves numbers only on libraries whose first 1,000 spliced reads lack it.
Instrument: a fixture with the tag absent from the first reads.

### vcap-dna-only-objects-lose-gdna-at-full-depth
`priority: later — Tier 3, local real data with truth by read name · kind: defect · 2026-09-28`
At full depth the VCaP mix reads P/O 0.922, and the deficit sits where no truth label is needed: P − O = −365,700,
+17,267 on structural objects and −382,967 elsewhere — −345,473 on no-RNA objects (−15.0 % of their O, 7.35 points of
the library), −58,503 on mixed, +21,009 on no-gDNA. It is already in pass 0 (P/O 0.930). By object size: [1,3)
−26.5 %, [10,30) −12.2 %, [30,100) −11.4 %, [100,1k) −5.0 %; the c/√n floor at c = 0.39 predicts 5.5 % at n = 50 and
2.3 % at n = 300, so objects of 30–1,000 fragments lose about twice that. Mechanism not isolated;
`ISSUES: the-gdna-landscape-collapses-at-low-depth` excludes it.
Instrument: `~/Downloads/rigel_runs/prototypes/2026-09-28_depth_lowg/pages/explanation_depth.md`; the VCaP-mix scorer
under `v2/scripts/`.

### em-gdna-exceeds-calibration-on-the-vcap-transcriptome-half
`priority: later — Tier 3, local real data with no truth · kind: question · 2026-09-28`
On the VCaP transcriptome half (no added gDNA) the EM assigns more gDNA than calibration called, under every od arm,
and about 2.5× the gDNA prior it was handed. At full depth: calibration 86,650 (od 0), 81,612 (shipped), 69,004 (the
joint fit or the 0.2 value); the EM 96,362–97,856 (0.695–0.705 %); the EM's gDNA prior sum 33.1k–39.8k. Calibration's
`library_gdna_fragments` includes intergenic regions and the EM's total excludes `intergenic_total` = 9,521, so like
for like the EM plus intergenic exceeds calibration by about 24–53 %. Across arms the EM shrinks calibration's 20 %
od-driven swing to 1.5 % (transcripts move 1.8–2.8k of 13.4M). The half is not certified gDNA-free, so neither number
can be scored.
Instrument: `~/Downloads/rigel_runs/prototypes/2026-09-27_strand_od/results/robust_no-gdna/`.

### latent-defects
`priority: later — Tier 4, each its own commit with a falsification test verified failing · kind: defect · 2026-09-28`
Correctness defects from the 2026-09-27 review and the cleanup. A test's expectations may move; no panel number moves
unless marked.
MOVES A NUMBER (each alone):
- `index.py`'s `s_tx.nrna_n_contributors = …` should be `+=`: a single-exon transcript covering two merged spans keeps
  only the last count (the nascent table's `n_contributing_transcripts`).
- `Scenario.build_oracle`'s fixed-total branch never reads `gdna_fraction`, so six fixtures simulate NO gDNA, among
  them `test_scan_order_independence`, whose docstring says it puts gDNA in introns. Raise on the combination, then
  repair; the expectations move.
- `lift_choices` keys records on (ref, start, end, align_strand, sj_strand) while the native key adds the introns and
  the hypotheses, so tied records that differ in introns pool: oracle truth can replay the wrong record, or `drain`
  can raise. Frequency unmeasured; only oracle truth moves.
- The simulator's fragment-length sampler: `.astype(int)` truncates (mean frag_mean − 0.5 while `fl_pmf` is evaluated
  at integers), and an empty [frag_min, frag_max] window loops forever. Land with the next re-simulation.
- Antisense at p_sense = ½: the deterministic path books a mismatched fragment as sense, the EM as antisense. Goes with
  the scoring unification (`ISSUES: hygiene-ledger` (h)).
- The leading intron: a first intron starting at the fragment start makes `build_segments` emit no leading segment
  and `left_sj` mislabel every junction end. Fix the executable specification first; it re-caches.
LATENT (no current path reaches it):
- `run_parallel`: a throw from fn(0) leaves `task_` racing, and a worker throw terminates; the ψ and EM callers do not
  catch.
- nanobind: `nb::cast<A>(h).data()` in `em_solver.cpp` and `solve_kernel.cpp` converts into a temporary on a dtype or
  order mismatch and dangles (GIL released); outputs without `.noconvert()` write into a discarded copy; the scatter
  bindings downcast float64 silently; `StreamingScorer` has no `nb::keep_alive`.
- `--overhang-alpha 0` / `--mismatch-alpha 0` give 0 × −inf = NaN: require alpha in (0, 1] or carry the gate as a bool.
- The scan cache: `deposit_digest`'s fixture offers `hypotheses=()`, blind to gap arbitration and deferral;
  `_payload_from_parts` skips `from_dict`'s conservation checks; the two `.npz` files and the manifest are written in
  place, with no temp-and-rename. One validating `AccumulatorPayload` constructor would serve the scan, the cache and
  the drain.
- GTF parsing: a final attribute without `;` is dropped; a missing `transcript_id` escapes `warn-skip`; an exon on
  another chromosome or strand is merged silently.
- Simulator config: a capture entry without `enabled:` and without probes yields an empty `CaptureConfig`; gene and
  transcript keys are unchecked (a typo runs on the default of 100); an empty section fails with a TypeError;
  `sim_command` creates the output directory before validating.
- `total_abundance`: IndexError when the last reference owns no regions.
- `splice_graph._ref_slices`: unsorted input drops exons silently.
- The transfer RNA coordinate's fallback keys on opportunity, not count, so all-zero-count single-strand exons build no
  RNA lane.
- An arm that injects `rna_sense_frac` without `n_rna_obs` silently kills the strand channel.
- `fast_exp` on an all −inf row: `static_cast<int64_t>(NaN)` is undefined behaviour.
- `connected_components_native` returns a 2-tuple early and a 5-tuple normally.
- `FragmentAccumulator::finalize` drops `t_offsets_`'s leading 0.
- `_resolve_core`'s non-empty-exons precondition is stated only in its doc comment (all four callers guard it): an
  assert.
- `genomic_span` is 0 for a transcript at position 0 (a truthiness test).
- `set_sj` narrows to int32 unchecked, justified only by a census.
- The buffer's finalizer, run by garbage collection on a thread holding a logging lock, deadlocks joining the writer.
- The report's meta line puts the sample and index names into `innerHTML` unescaped. (Its capture note is in
  `ISSUES: the-format-changes-to-batch-before-release`.)
- An unknown quant YAML key raises without listing the known keys or hinting at an older build, so a 0.7.1
  `config.yaml` is refused unhelpfully; MANUAL's "rerun the exact analysis" needs a caveat.
- A splice blacklist present without a recorded store is dropped with no log line: a one-line warning
  (`ISSUES: an-unrecorded-splice-blacklist-is-dropped-silently`).
- `transfer_rows.h`'s `level_bound_row` `std::max(v, 1e-12)` never binds on a reachable input: delete it together
  with an assert v > 0 in the binding (today v = 0 reads 1e-12; without the floor it would be NaN).

### performance-memory-bounded-solve
`priority: later — Tier 4: the numeric no-ops, the first two unparked (owner, 2026-09-28); the rest of the thread parked (owner, 2026-09-19: the method is the focus) · kind: build · 2026-08-17`
A deep run must be fast enough to iterate on, and memory-bounded. The port (the sweep as one native call) and the
work outside calibration are done (`DESIGN.md` §6b.15). THE BASELINE (VCaP, 18,568,456 fragments, `--threads 8`,
interleaved pairs, 2026-09-19): 125 s — calibrate 39.5 (four sweeps 24.0, landscape fits 4.1, its Python 10.3),
quant 33.7 (locus EM 14.4, capture effective lengths 5.2, scoring 5.7), the scan 28.4, the second pass 13.0, a
second fragment-length fit 2.8, the index load 6.4; peak RSS 10.3 GB, in quant.
OPEN, ranked (a cProfile share ranks, `TRAPS: a-profile-share-is-a-ranking`): ① the capture effective lengths
(quant's 5.2 s) — the short-template taper table that held most of it is deleted with the per-base rule, and the
one shared rule (`capture_eff_length.transcript_objects`, a closed form per object) is re-timed on the deep library
before it is ranked; ② the second pass's
factor combiner, 355,180 calls on 2–4-element arrays, hoistable to segment operations over the hypothesis CSR
(exactness to check); ③ the landscape fits, ~5 s, unexamined; ④ the sweeps' 24.0 s, where the blur's loop
interchange is priced at −1.5 s for a summation-order change and is the owner's call; ⑤ past 100 M fragments.
REFUSED as targets: the scan (at its floor for eight threads) and the index load (mostly native already).
⑥ THE REVIEW'S NUMERIC NO-OPS, each proven with `rename_identity.py --bam` from back-to-back pairs; UNPARKED (owner,
2026-09-28) for the first two, the rest batched:
- `get_loci_df` builds `locus_ids == lid` per locus, O(loci × transcripts): 0.48 ms per locus at 450k transcripts,
  about 19 s at 40k loci, in the CLI writer; `np.bincount` fixes it.
- `iter_chunks_consuming` holds the previous chunk while loading the next, and `_SpillWriter._run` keeps its last
  chunk while blocked: with chunks of up to 1 M fragments against a 2 GiB budget, a `del` fixes both.
- `build_region_partition_arrays` runs 3+ times per run, `build_sj_arrays` twice.
- `--annotated-bam` calls `add` per fragment where a batch API exists.
- The EM clears `local_map` up to the largest global index per locus, and allocates per E-step (heap rows for
  k > 512; the task list rebuilt).
- The scan path copies `gap_introns` per hit, one malloc each.
- The capture sampler's memo evicts across the mRNA and gDNA spaces, and its overlap lookup is O(T) per probe.
- The solve builds whole cubes for slots it then discards, and `pool_of` rebuilds the pool when the task count
  changes.

Instrument: `profiling/profiler.py`, `profiling/sweep_replay.py`.

### hygiene-ledger
`priority: later — Tier 4, batched between A/B windows; a numeric no-op lands between any two · kind: hygiene · 2026-08-31; the review of 2026-09-22`
What the reviews left, each its own commit; content-only changes keep the collected count unchanged.
(a) ROTTEN BUT LIVE, moving an instrument's numbers when repaired: `quant_accuracy`'s oracle arms are undrained
(documented there).
(b) DEAD, DEFERRED: `region_span_count`, tallied per fragment by the accumulator and carried through the payload,
the substrate and the caches, read by nothing (the retired length channel's). Deleting it changes the payload schema
and re-caches both panels; step 1b re-cached every panel without taking it, so it now rides the leading-intron fix's
re-cache (`ISSUES: latent-defects`; owner, 2026-09-28). Until then `accumulator.h`'s "each population stores only the
channels something READS" is false. The native `FragmentAccumulator` also still reserves, fills and exports
`sj_strand` and `merge_criteria` per fragment, which the buffer drops at `from_raw`.
(c) PRODUCTION-DEAD, KEPT AS TEST OR INSTRUMENT SURFACE: the resolver's `intron_bp`, the tested half of its overlap
profile; `rigel.sim.benchmark` and `Scenario.build_oracle`, test tooling inside the package;
`region_init.has_own_composition_evidence`, which tests use as a check; `priors._project_regions_to_loci`, kept for
`prior_vs_oracle.py` (`assemble_priors` calls `_region_locus_shares` directly); `fl`'s `GdnaContrast` and
`GdnaRealized`, computed and discarded (decided at `ISSUES: the-realized-gdna-length-law-reads-rna-counts` step 2);
`AbundanceLandscape`'s census fields (`modes`, `depleted`, `enriched`, `n_train`, `AbundanceMode.width`), read only by
tests until the format batch gives `depleted` a reader; the vendored `_cgranges_impl` extension, still compiled,
optimised, LTO'd and shipped for a test-only `query()`, while `index.py`'s docstrings present it as a quantification
facility.
(d) CLAIMS NOT RE-DERIVED on the current tree, left standing: that most in-scope error sits at the simplex
vertices (`simplex_logodds`, a relay-era measurement); `sweep`'s refused deferral of UNIDENTIFIED slots to the
prior (priced 2026-07, with no refusal entry here); `region_geometry`'s "no per-region spliced floor" A/B
(relay-era); and `fl`'s crossing pools called "gDNA by structure" because mature RNA never crosses an
exon|intron boundary, while RNA that has not spliced there does.
(e) GATES THAT CHECK LESS THAN THEIR NAME:
- FIRST, since any failure is a regression: `test_THE_FIXTURE_REALLY_DOES_REORDER_THE_BUFFER` asserts that a thread
  race fired; it failed in two full runs and passed on rerun.
- `test_sweep.py::test_gdna_sweep_zero_gdna_pin_and_monotone` asserts `f_g < ½` on slot 3, the AMBIG|intron−
  boundary, while the AMBIG region it means ends at 0.949 under the silent policy, and nothing in it checks
  monotonicity; the mature-exon chain's tests give the same `f_g` with the junction spliced or not, so none of them
  exercises the junction reads.
- `test_conserved_mass.py::test_the_mass_is_the_PER_BASE_attribution` allows `count` whole fragments where its
  docstring claims `count` half-ulps; `test_gdna_strand_fit.py` and `test_region_geometry.py` still build float64
  banks as uint64 fixed point in their fixtures.
- The HANDFUL gate (`test_gdna_strand_fit.py`) calls the ρ = 0 pair-count moment "biased LOW" and asserts mean(raw) <
  0.17 on 8 fixed seeds; in a 300-draw Monte Carlo of its own fixture the raw mean is 0.2007 (sd 0.117) against a
  fitted sd of 0.0044, and only 19 % of random 8-seed batches pass. Restate it as a variance claim, or retire it with
  the raw field at `ISSUES: strand-overdispersion-one-shared-value` step 1.
- The report test's synthesized v3 summary omits `sj_blacklist_size` and `sj_blacklist_loaded`, so it renders
  "detection off", a shape no real run produces; no test asserts the splice-artifact note.
- The buffer's finalizer test asserts nothing about "still writes what it holds", and its writer-error test calls
  `cleanup()` before consuming, the reverse of production.
- The scan-worker exception gate runs one worker, so deleting `output_queue.abort()` stays green.
- The annotated-BAM ZS check is set membership only, so a label swap passes.
- The fragment-length proof's L check is still guarded by `if contained is not None`; two spec tests carry leftovers.
- The docs-citation gate never checks that it scanned anything (an empty glob or a `.hpp` suffix passes); a URL
  containing `/docs/` fails with no remedy; the TRAPS heading regex is written three times.
- The log-odds window's single source is ungated: restoring a default stays green, and about 23 literal 10.0s sit in
  the tests.
- `test_a_toy_with_NO_anchor_pool_falls_back_to_the_largest_basin` holds by construction on its unimodal fixture; a
  bimodal one (0.731 / 0.269) would bite.
- The deleted float64 pass-through gate was never retargeted at the live boundary masses; `test_substrate.py` still
  says `sj_mass` arrives per strand "for artifact detection".
- `test_a_blacklist_built_from_a_store_loads_into_the_resolver` asserts a size that Python sets, and stays green with
  the native map stubbed.
- The step-1c instrument fixes have no gates: `policy_benchmark`'s apart grouping, `calibration_vs_oracle`'s APART
  strata and its `--json --jobs` merge, `prior_vs_oracle`'s strata, the `find_spec` guard.
(f) SMALL: the `--no-mappability` flag names a store whose mappability the index no longer reads (it opts out of
the splice blacklist); three unread scalars in `solve_kernel.cpp`'s positional aggregates; `pgzip` is an undeclared
optional import in the simulator; eight test sites import `rigel._resolve_impl` for names `rigel.native` re-exports.
- `_solve_impl` builds with `-g` only. Step 1b tried `-march=native` and IPO and reverted them: bit-identical on arm64,
  but on x86-64 they turn on FMA contraction in the ψ, transfer and solve kernels (clang: 0 FMA instructions at -O3,
  159 with -march=haswell) and move numbers (`TRAPS: a-kernel-edit-moves-ulps-under-fma-contraction`). RULED (owner,
  2026-09-28): they stay off; correct the two CMakeLists comments that are false for it. The 3bb35bdd message's
  "both proven bit-identical" is wrong for `_solve_impl` and stands corrected here.
- The numpy bit-matching machinery: the native solve still carries numpy's pairwise-sum order (`PW_BLOCKSIZE = 128`),
  the `lgamma` binding and named temporaries, with no Python twin left. Retiring it moves results by ulps and needs the
  goldens and the sweep replay re-recorded. RULED (owner, 2026-09-28): retire it after the strand-overdispersion and
  fl landings, in one golden re-record window.
Kept as coverage GAPS, not dead code: the five CLI command bodies, the silent policy through `calibrate`, the
simulator's sharded writers and its whole-genome grid, and the zarr splice blacklist.
(g) COMMENT AND DOC RESIDUE, one content-only sweep:
- Native: `bam_scanner.cpp`'s references to the deleted Python Fragment, "Fractional accumulator deposit", an
  orphaned reference, a lost rationale, the intron filter justified twice, `std::mutex` without `<mutex>`;
  `resolve_context.h` and `constants.h`'s history narration, "for compatibility" notes, a misplaced banner, four AF
  bits unused in C++; `scoring.cpp`'s `-Wreorder-ctor` warning and its gDNA per-hit comment omitting `n_cand > 0`;
  `fast_exp.h` listing three "call-site properties", one of them the function's own.
- Calibration: the sweep docstring's "counted in the kernel"; the test named for owned slots; "bound for the gates",
  although production reads `transfer_rows`; `region_geometry.py`'s "five populations" (there are four), the Axiom-0
  tell "nascent-RNA-active" (also in `signature.py`) and mature/phantom history; `result.py` citing the deleted
  `simplex_logodds._compose`; `calibrate.py`'s "two supports", a docstring that lost its conclusion, and the test
  helper `_zero_rna_opportunity` zeroing a field already 0; `strand_model.py`'s spliced model's qualification stated
  two ways, an undefined ε_CI, `posterior_95ci` not marked QC-only; the four texts stating "a quarter of the RNA side"
  as a general property (a deleted measurement on a rebuilt ladder).
- The EM, locus and CLI: `estimator.py`'s "every unambiguous fragment is deterministic" and a 102-character line
  hard-coding (N_t, 10); `locus.py`'s false (ref_id, start) sort claim (it sorts by categorical code) and its pointer
  to `calibration.priors` for the EM pseudocounts; `quant_command` omitting the `--annotated-bam` second pass;
  `report/substrate.py` ("omits the panels") and `report/__init__.py` (the extra is needed only for the charts).
- Tests: `test_estimator.py`'s "no EM at 0" (SQUAREM runs at least one step); `test_resolution.py`'s "(gap=0)";
  `test_buffer.py`'s "None for intergenic"; over-long test docstrings; "Boundary cases" headings from the
  edge→boundary rename.
- The simulator: `whole_genome.py`'s "…and filtering" banner; `splice_motif.py`'s over-long lines.
- Docs: DESIGN names the solver's `mrna_active`, which no longer exists, cites a test by its pre-rename name, and its
  compiler-flags row omits `-fno-finite-math-only` (as does `.github/copilot-instructions.md`); EQUATIONS §3.5g, DESIGN
  §3.5's "four ways", `accumulator.h` and `test_conserved_mass.py` still describe the deleted
  `region_contained_inv_opportunity_sum`, and one measurement table can no longer be reproduced; the executable
  specification and its C++ twin word one mass rule differently, and two payload comments point at inexact text
  (UNCERTAIN: check when touched); `test_chr.yaml` says the shadows are "≈ 5 % of RNA fragments" while the render
  realizes 0.76 %; the headers and TESTING list the YAML's extra keys without `blank:`; TESTING says `rigel sim`
  accepts `index:`, which the CLI now rejects; TESTING has no section for the VCaP read-name truth mix (A00839 exome
  DNA, HWI-D00127 transcriptome; its halves differ in fragment length, 307 / 244 bp, the confound the ladder forbids;
  read-name truth labels the RNA library's own gDNA as RNA; calibration can be scored only as slot-mass overlap).
- RELEASE-GATING: the MANUAL, README, report and cli texts still describe the deleted capture census or the old
  density-file semantics ("Capture on-target enrichment"; MANUAL's density row contradicting its correct neighbour);
  the MANUAL does not say that `sj_blacklisted` counts junction × record while the splice census counts fragments.
(h) SIMPLIFICATIONS, each proven a numeric no-op:
- `solve_one_block` and `transfer_prepare` each wire the chain, so the gates exercise the copy: one shared builder.
- Scoring is written twice, for multimappers and singles; unifying it also removes the p = ½ antisense disagreement
  (`ISSUES: latent-defects`).
- SQUAREM's VBEM and MAP branches (about 150 lines) differ only in the weight, the floor and the normalisation.
- `BamAnnotationWriter::write` re-implements `parse_bam_record`, and multimap detection is written twice.
- `TranscriptIndex.load` round-trips through frozensets: build the CSR directly.
- `intervals.feather`'s exons are read three ways: one walk would serve all three.
- Two functions named `contained_opportunity` sit in layer 2, and one ramp is written twice.
- `FragmentScorer` copies 13 fields into native: return the native scorer and store alphas.
- `_apply_scan_stats` reads 33 hand-listed keys with `.get(key, 0)`: iterate the dataclass and index strictly.
- The simulator's truth is spread over nine `ground_truth_*` methods and three file passes.
- Five instruments import `tests/calibration/_oracle.py` through `sys.path`; four scripts redefine `_shared.RUNS`.

### instrument-ledger
`priority: later — Tier 4, batched between A/B windows · kind: instrument · 2026-09-28`
What the instruments do wrong or leave out; each fix is its own commit, with its `--self-test` or a gate.
(a) `slot_truth.npz` counts boundary gDNA per crossing, in incidence units, so under capture its sums exceed the
library's true gDNA (slot-frame O 1,281 / 12,485 at g001 / g01 ON against 1,170 / 11,700 fragments; P/O 1.015 in the
slot frame against 1.055 in fragments), while `policy_benchmark.py` calls its Σ|Δ| "in fragments". Document the frame.
(b) No instrument output records the code revision, the `--set` overrides or the prototype arm
(`TRAPS: an-ablation-that-never-ran`).
(c) `sweep_replay.py` pickles slots dataclasses; on Python 3.12 an old capture loads by position, shifted silently, so
every capture taken before 1747b282 must be retaken.
(d) `build_test_reference.py` puts `src` first on `sys.path` unconditionally; fix it before the contaminated panel's
in-tree stage.
(e) The ladder-report skill hard-codes the four ladder rungs (StopIteration on g25 or g001).
(f) `panel.py --index` reaches the cache and oracle stages but not `simulate`.
(g) `calibration_vs_oracle.py` has no arm that carries the true capture efficiencies, so no calibration A/B sees the
capture-contracted length move; `quant_accuracy.py --arm oracle_ruler` is the model.

### debug-capture-memory-is-unbounded
`priority: later — Tier 4: the fix is in src; the one-real-genome-job rule stands meanwhile · kind: instrument · 2026-09-28`
`calibrate(_debug=…)` builds a `SweepCapture` for every sweep with no size or region bound; on the whole human genome
it reached 25 GB in 20 s and nearly crashed the machine (2026-09-28). `ruler_vs_truth.py` always passes `_debug` and
does not refuse a real index, and no capture is region-restricted or streaming.
The per-session memory guards miss Python from other environments, processes under 1 GB and growth between 3-s polls.
For scale: a normal quant peaks at about 3 GB, and `rename_identity.py --check` on one ladder condition holds 8.9–10.7
GB. RULED (owner, 2026-09-28): `_debug` refuses an index whose manifest does not mark a panel, and the capture is
bounded to named regions. The one-real-genome-job rule stands meanwhile.
Instrument: none in the tree; the session's guards (`mem_guard_env.sh`, `mem_guard_fp.sh`).

### the-format-changes-to-batch-before-release
`priority: later — Tier 4, before 0.8.0 while summary.json schema 3 is unreleased; the index bump before the cluster's working-index build (owner, 2026-09-28) · kind: decision · 2026-09-28`
Changes to an output or on-disk format, batched so each format changes once. RULED (owner, 2026-09-28): the whole
batch, before 0.8.0.
(a) `summary.json` (schema 3, unreleased):
- one od field instead of two, with its information and an evidence label, so a no-evidence 0 reads apart from a
  measured 0 (`ISSUES: strand-overdispersion-one-shared-value`);
- `gdna_fraction` without the spliced intergenic fragments: `cli.py` adds `stats.n_intergenic`, which includes them,
  and a spliced fragment is certified RNA (Axiom 0). The count is so far inferred (27 at g001 OFF, 29 at g01 OFF, from
  the three-exon shadows), not read from `n_intergenic_spliced`; `quant_accuracy.py` copies it into the release
  report's gDNA pool;
- the depleted gDNA density, the mode `split_basins` computes inside `located_enriched_mode` and discards: it restores
  the report's on-target / off-target fold (dividing by `gdna_density_global` instead is biased, 25× against about
  1,000× from the modes) and gives `AbundanceLandscape`'s `depleted` census a reader; it touches about six test
  construction sites and two instruments;
- `gdna_density_kde.feather` and `gdna_density_regions.feather` renamed: they hold the TOTAL-density landscape;
- the report's capture note, which says "no calibration" for a 0.7.1 or schema-2 summary that has a calibration block.
(b) The index (`INDEX_FORMAT_VERSION` 8): drop the written-never-read columns (`transcripts.feather`'s `abundance`,
`nrna_abundance` and `n_exons`; `sj.feather`'s `interval_type`), and add the duplicate map as an alias map
`dropped_t_id → kept_t_id`. Every index is rebuilt, the cluster's included, so land it before the cluster's
working-index build and the cluster rebuilds once. `ISSUES: overlapping-synthetic-shadows`' index merge is a separate
decision.

### flgap-panels-stale-nascent-model
`priority: later — Tier 4, the panels: re-simulate, then delete fragment_share (owner, 2026-09-28) · kind: decision · 2026-08-22`
The two fl-gap side panels were not regenerated in the sparse-nascent rebuild and carry the retired uniform
nascent model, so a claim spanning the ladder and a side panel varies two things. RULED (owner, 2026-09-28):
re-simulate them, then delete the simulator's retired `fragment_share` nascent mode, whose only live users are these
two configs (the test-chromosome configs carry it only in `nrna` blocks that are dead under `abundance.mode: file`,
so odg05 and the depth family are unaffected). These panels are where `ISSUES: eb-shrinkage-magic-ess` is visible,
and the real VCaP halves differ in fragment length by 63 bp.

### expand-the-gdna-spectrum
`priority: later — Tier 4, the panels · kind: decision · 2026-08`
Fill the gDNA spectrum (1, 5, 10, 25 % up past 90) without multiplying benchmarks: a level is justified by a measured
transition and crosses a reduced set of the other axes until an interaction is shown; each condition costs a simulate,
two caches and a certification (`sim/panel.py`). The junction-probed test-chromosome twin predates the capture-physics
change and is stale; RULED (owner, 2026-09-28): retire it, since the capture-length campaign that needed it is closed
for now. See `ISSUES: flgap-panels-stale-nascent-model`; the depth family's unstranded rows are ruled in
`ISSUES: the-gdna-landscape-collapses-at-low-depth`.

### a-pure-gdna-library-reads-as-nascent-rna
`priority: parked — Tier 5, the deferred stratum; kept as an open challenge (owner, 2026-09-28) · kind: defect · 2026-09-28`
Rigel reads a library of pure genomic DNA as 93 % RNA. The case is the exome-DNA half of the VCaP mix, 4,671,916
fragments, separated by read name. It is unstranded (strand specificity 0.5008) and capture-enriched (calibration's
gDNA reference density is 180× its genome-wide density).

Where its fragments went:
- synthetic nascent RNA, over 53,962 spans: 3,781,970 (81 %);
- annotated multi-exon mRNA: 452,338 (9.7 %), only 4,213 of which come from spliced fragments;
- annotated single-exon transcripts: 109,003 (2.3 %);
- gDNA: 313,624 (6.7 %).

The library has 3,225 fragments spliced at an annotated junction (0.07 %), about what alignment artifacts alone
would give.

This is the deferred unstranded × capture-ON stratum at its extreme. The gDNA fraction cancels from the strand mean,
and capture depletes the intergenic space that would otherwise measure the gDNA density, so a gene span's unspliced
coverage fits nascent RNA as well as gDNA.

What could still tell them apart is the SPLICED FRACTION. Mature RNA at a multi-exon transcript produces
junction-crossing fragments at a rate its geometry and the fragment-length law set (layer 2's `sj_opportunity`).
Here the multi-exon mRNA assignment carries under 1 % spliced support where real RNA carries tens of per cent. And by
the nascent scope ruling, nascent RNA is sparse beside mature RNA, which it cannot outweigh 7:1 as it does here.
Whether either becomes a likelihood term, and how a pure-gDNA library should be recognised, is open.

Its plausible in-mix face: inside the VCaP mix the EM sends the true DNA fragments to RNA at 1 % / 10 % / full depth
(sampled draws, seed 0, not fractional) as synthetic spans 0.254 / 0.185 / 0.099, other synthetic RNA 0.014 / 0.008 /
0.004 and mRNA 0.118 / 0.059 / 0.044: 39 % / 25 % / 14.7 % in all.

Measured with the strand-overdispersion harness, `rigel quant` under every arm (the arms change nothing here):
`~/Downloads/rigel_runs/prototypes/2026-09-27_strand_od/results/robust_no-rna/`.

### unannotated-transcription-is-booked-as-gdna
`priority: parked — Tier 5: a DESIGN ruling after 0.8.0, recorded now (owner, 2026-09-28) · kind: design · 2026-09-28`
By Axiom 0 an unannotated slot admits only gDNA, so real unannotated transcription is booked as gDNA twice — by
calibration (structural objects give P = count) and by the output's intergenic rule — and it biases any pooled
structural rate a landscape repair might read. DESIGN §0b says the tool does not model it and calls that safe for a
pooled level; `ISSUES: the-background-dispersion-assumes-a-pure-intergenic-pool` covers one route only. MEASURED:
`test_blank` books 8,844 of 8,844 RNA fragments as gDNA at g001 capture OFF and 8,889 of 8,889 at g01, its eight
shadow transcripts all of the false-negative mass (12.9 % / 11.5 % of transcript Σ|Δ|); VCaP's structural P − O at
full depth is +17,267 — intergenic +8,654 (12.4 %), gene edges +8,613 (14.1 %) — at a median minor-strand share of
0.000 (74.5 % of 208 objects at ≤ 5 %); on the test chromosome's g00 the false gDNA is 100 % shadow transcription.
Admitting RNA at unannotated slots on stranded data needs a ruling and has costs: the closed floor would point toward
RNA on truly pure intergenic objects, and that RNA would have no transcript to go to. RULED (owner, 2026-09-28): for
0.8.0 only the output's spliced intergenic count changes (`ISSUES: the-format-changes-to-batch-before-release`).
Instrument: `~/Downloads/rigel_runs/prototypes/2026-09-28_depth_lowg/results/blank_arm.txt`,
`pages/explanation_depth.md`.

### ruler-witness-geometry-on-transcript-panels
`priority: parked — Tier 5, deferred past 0.8.0 (owner, 2026-09-20): Rigel stays panel-agnostic and takes no panel input · kind: defect · 2026-09-15`
Only a transcript holding the junction a probe spans binds that probe whole, so the extra capture of junction-spanning
fragments is ISOFORM-SPECIFIC, and the EM splits a gene's shared fragments by the ratio of its isoforms' capture-aware
lengths. gDNA holds the probe's parts apart and binds the better one, and at zero gDNA there is no witness at all; the
shipped junction price reads it from gDNA at the junction's low and high boundaries by conservation of bases, with no
panel input — 0.959 of the simulator's own junction capture over 395 junctions of 60 probed transcripts in aggregate,
but per junction too noisily (`ISSUES: the-junction-price-is-noisy-within-a-gene`). MEASURED on the rebuilt ladder
(the simulator binds every fragment through one contiguous part of a probe since 2026-09-19; earlier numbers measured
an asymmetric half-match the owner ruled unphysical), stranded × capture-ON transcript Σ|Δ| as a share of the true RNA
at `g00` / `g05` / `g50` / `g98`, fractional: the per-base rule 6.53 / 3.49 / 5.50 / 31.03 %, the simulator's own
lengths (`quant_accuracy.py --arm oracle_ruler`) 1.34 / 1.71 / 2.90 / 24.66 %. The error sits in highly expressed
multi-isoform genes whose TOTALS are right, and a ruler is judged by the WITHIN-GENE spread of its error
(`TRAPS: judge-a-ruler-by-its-within-gene-spread`). The junction-probed test-chromosome twin is stale; its retirement
is `ISSUES: expand-the-gdna-spectrum`'s. THE ONE SHARED RULE'S RESIDUAL (2026-09-23, 60 probed transcripts against the
simulator): exact (1.000) on the objects where a transcript's fragments are gDNA's — contained pieces and cuts away
from junctions, exon edges and ends, 41 % of a probed transcript's yield — and its junction cuts, 40 % of a probed
transcript's yield, are what the junction price covers. Its cuts within a fragment of an exon edge are captured 1.10×
more than gDNA prices them (19 % of the yield), and the simulator binds a fragment through its best single probe part,
which captures a little less than the sum when both exons carry separate probes: that is the ~2 % the classes stay
apart with the true gDNA on every object (`ISSUES: the-gdna-component-length-rule-differs-from-the-transcripts`).
THE DECISION IT WAITS ON: the probe design is the one observable that sees isoform-specific capture directly, and
Rigel reads no panel (`DESIGN.md` §7.2). Real panels often span junctions and the design file is usually
unavailable (owner, 2026-09-19), so the candidates are data-derived: the junction price from gDNA (shipped, its
precision open); a capture field fitted from the coverage shape around each probe and the gDNA footprint; a
per-kit capture profile learned across a cohort; and a sparsity prior on isoform support as the safety net. The
spliced reads at each junction are refused (`ISSUES: the-spliced-read-junction-price`). `quant_accuracy.py --arm
oracle_ruler`, `ruler_vs_truth.py`.

### overlapping-synthetic-shadows
`priority: parked — Tier 5, an owner decision on the index; not a 0.8.0 number · kind: design · 2026-09-20`
`create_nrna_transcripts` clusters TSS/TES within a tolerance per strand and manufactures one synthetic
nascent entity per merged span, so a gene with several isoform clusters carries several near-coincident
shadows: 6,366 of the ladder index's 6,919 overlap a same-strand twin, 5,879 at ≥ 90 %, degrees up to 240.
The twins are indistinguishable to the EM, which puts a family's mass on an arbitrary holder (the live
entity is the holder in 78 % of live families off capture, 42 % on), so every PER-ENTITY nascent number is
a coin flip and the family is the unit that can be scored. They do not cause the capture-ON siphon —
collapsing every family to one component in the shipped EM reads 843k against 692k nascent at `g50 ss.99
ON` — and they cost 7 % of live nascent off capture (live families 942,567 against 1,013,400 at
`g50 ss.99 OFF`). The pool rows and the transcript table do not see them (synthetic rows are dropped);
per-entity diagnostics do. The decision is whether the index should merge them
(one span per gene per strand) or the scoring should. `quant_accuracy.py`, the per-family scoring in
`~/Downloads/rigel_runs/prototypes/2026-09-20_siphon_mechanism/`.

### yield-variance-beside-the-count
`priority: parked — Tier 5, with the per-transcript prior lane; the release ships the contraction as it stands (owner, 2026-09-17) · kind: build · 2026-09-17`
The capture-contracted yield is a posterior mean and carries no variance, so a count on a small yield reads as a
large abundance with the same error bar as any other. Measured on VCaP on the per-base rule by drawing every piece
efficiency from its posterior on the landscape grid and re-running the EM (eight draws, the EM's seed fixed; the
session's `yield_draws.py`, scored by `draws_census.py`): the counts themselves are stable — the dominant isoform
under the draws' mean equals the point estimate's in 98.5 % of multi-isoform genes with ≥ 20 fragments and is
identical in every draw for 90.6 % — while the yield's uncertainty is a CV of 0.006 / 0.061 / 0.317 at the 10th /
50th / 90th percentile of transcripts with ≥ 20 fragments, above the Poisson CV (median 0.100) for 36.4 % of them.
So the variance need not be propagated through the EM; it is a per-transcript number to publish beside the count:
the yield's posterior sd from its objects' posterior variances under the landscape, so that a user's abundance
carries `CV² = 1/k + Var[eff_t]/eff_t²`. On the one shared rule the form runs over the objects — a piece at its
contained share, a cut at its conserved share, a junction's price a signed sum of four objects' posteriors
(`EQUATIONS.md` §11); the objects' posteriors are what `capture_efficiency` already computes, and the second moment
is one more `w @ clipped²`. Derive → gate against the draws' spread on VCaP → `src/`: a result field per object and
an `em_effective_length_sd` column. Its consumer is the allocation of the pooled RNA prior across a locus's
transcripts (`ISSUES: per-transcript-prior-lane`): today the transcripts duel for the ambiguous fragments with no
per-transcript prior, and a variance per transcript is what an allocation better than equal shares would read
(owner, 2026-09-17). Instrument: the draws' spread is the truth the analytic form is gated against.

### capture-premise-untested-on-cdna
`priority: parked — Tier 5, watch: no library in hand can test it · kind: risk · 2026-09-17`
The contracted length reads the panel from gDNA and applies the same efficiencies to cDNA: the premise that
off-target cDNA is depleted as off-target gDNA is, which the simulator satisfies by construction and which
cross-hybridisation of cDNA to paralog probes or nonspecific binding could make milder by an order of magnitude.
Its exposure is not the scale (TPM is on the plain length by ruling and counts see only ratios within a locus) but
the isoform split: on VCaP 4,677 of 12,039 multi-isoform genes with ≥ 20 fragments change dominant isoform between
the contracted and the plain yield on the per-base rule (`DESIGN.md` §7.2, the yield's two consumers), two thirds
of the isoform-level fragment mass with them, stable under the yield's posterior — so the flips are the model's
answer, right if the premise holds. The one experiment that settles it: a sample sequenced both ways, panel and
whole transcriptome — the count ratio on unprobed transcripts is the cDNA depletion directly, and the capture
efficiencies give the gDNA depletion on the same objects. The lab has no such pair (owner, 2026-09-17); the
session's `lever_census.py` is the read-only smoke test to re-run on any new captured library, and `count_unambig`
beside `count` is what tells a user which isoform assignments rest on shared fragments alone.

### the-atom-at-an-unwitnessed-both-strand-slot
`priority: parked — Tier 5, accepted as a limit of the information (owner, 2026-09-14) · kind: known limit · 2026-09-14`
The tilt atom (`DESIGN.md` §6b.15.13, `EQUATIONS.md` §9f) admits "all of this slot's RNA is on strand s" unless
a held RNA level on the other strand rules it out. On stranded data at a both-strand slot whose minor strand
has no witness — no junction certifies it and no single-strand piece of its own exists anywhere, the case of
a single-exon gene wholly inside another gene's exon — the strand split cannot tell "pure s with gDNA" from
"both strands without": the pure hypothesis fits the split with no parameter and the atom pushes toward it.
Measured on the golden toy `antisense_contained` (1,000 fragments): the exon∩exon slot's truth is 0 gDNA;
prior-free the continuum reads 0.244 and the atom 0.483; with the toy's landscape of 8 training regions
0.005 and 0.242 — the antisense transcript's count 81 → 0 and 177.6 false gDNA fragments in a gDNA-free locus.
The mono shared-exon stress (no junction) doubles its false gDNA at 2–20 % minor. On the ladder, where the
landscape trains on ~15k regions, the same slots cost a few fragments each (the (0.2, 0.5] band, `g05 ss.99 ON`
4,942 → 7,005, `g00 ss.99 OFF` 35 → 73) against the strand-pure band's −11k / −13k. THE STANCE: this is not
enough information, and it is accepted as such — no presence witness (a locus with RNA elsewhere does not
imply this slot is expressed; a per-transcript active set would add bookkeeping and no measurement, since the
failing structure has no node where the transcript can be measured alone). On a real library the landscape
prior is the deciding voice; its location floor (`DESIGN.md` §7.1 rule 4, 2026-09-14) is the lever that was
pulled, and the stranded zero controls read 265 / 172 with it — while the 1,000-fragment golden, whose prior
now fits from about four anchors, reads 200.8 false gDNA (from 177.6): the toy's limit. On unstranded data the
atom is inert
(κ = ½ makes the three hypotheses equal) and such a slot is the prior's entirely. Watch it on the census
(`ISSUES: the-tilt-census-as-an-instrument`) and the mono stress; a real library that shows it would be a
nested single-exon antisense gene on a shallow stranded run.

### two-sided-exon-row
`priority: parked — Tier 5, deferred by the owner (2026-09-13), with the enrichment witness · kind: problem · 2026-09-04`
On unstranded data an exon's held row is the intron's composition through the face map, whose upper side is
the map's plateau above its ceiling — a lower bound on gDNA — so at pass zero an unstranded licensed exon
reads ~9× its true gDNA (+2,000 % at `g25 ss.50 OFF`) and forwarding it compounds the bias. A toy gate
(retired 2026-09-28 with its harness) read a factor of 88 (exon |Δf_g| 0.762 beside a pure-gDNA intron, 0.0086
beside a nascent-bearing one); on the ladder the class is 4–6 % of the unstranded error with
`transfer` at parity there, on the test chromosome 20–32 % and transfer worse (2,647 → 4,056 at
`g50 ss.50 OFF`). Every cap, two-sided row or wall is refused where junction probes enrich the flux more than
the crossing (`ISSUES: the-certified-flux-row-as-a-level`, `ISSUES: two-sided-exon-row-forms`,
`ISSUES: the-discrepancy-priced-cap`, `ISSUES: the-two-sided-level-lane`,
`ISSUES: the-wall-above-the-face-map-ceiling`). What would license a two-sided level is whether this library's
gDNA is enriched, witnessed on unstranded data by silent genes' exons against the intergenic density — the
landscape prior's job. The first step when taken up is that witness's derivation; judge at `R exon
(licensed)` and the walled classes, halves apart. `policy_benchmark.py --by-class`.

### the-lower-bound-noise-ratchet
`priority: parked — Tier 5, with the enrichment witness · kind: defect · 2026-09-05`
A level from an RNA-rich node's own strand profile has a mode that is noise around zero; its lower side bounds
its neighbours and the tightest noisy neighbour wins (ladder `g05 ss.99 OFF` 44,714 → 45,076; test chromosome
`g00 ss.99 ON` 64 → 216, 187 of it at `capcluster_ab`'s inner termini). A two-sided own-profile level keeps
fewer fragments overall, so lower-only stays. The gDNA lane's EDGE level does it too (2026-09-14, on a toy
of the encompassing locus): a shallow single-strand flank (404 fragments, truth 0.530, local solve 0.546)
reads 0.596 under the intergenic neighbour's Poisson level — a lower bound at that neighbour's sampled
density, 0.298/bp against the flank's realised 0.27, a 1.6σ excursion the hop's price blurs but does not
move. The cure is the enrichment witness `ISSUES: two-sided-exon-row` waits for. The ladder's `g00` rows
(`calibration_vs_oracle.py`).

### flux-floor-dispersion
`priority: parked — Tier 5, with the transport-dispersion decomposition · kind: question · 2026-09-08`
The certified flux at a junction is a lower-sided estimate of the exon's strand RNA level priced by the pair
(`DESIGN.md` §6b.13); the route rate scatters beyond counting (median −3 %, 5–9 % at depth; 0–40 % over on
nine block readings) and at a pure-RNA exon the price cannot see it, so a lucky over-read is a sharp floor a
few points too high. Its decomposition instrument, `transport_dispersion.py`, was retired 2026-09-14 (in git).

### splice-out-premise-bias-uncorrected
`priority: parked — Tier 5, the owner's call · kind: decision · 2026-09-02`
The splice-out message assumes spliced and unspliced fragments at one face share capture affinity; measured,
the premise fails as a bias under capture (`log a` ≈ 0 off, +0.28 under exon probes, +0.78 under junction
probes), so only a subtracted level — a cross-locale fudge, the owner's call — would correct it. No instrument
yet.

### the-tilt-census-as-an-instrument
`priority: parked — Tier 5 · kind: build · 2026-09-13`
"Where does the strand tilt matter, and how does the tool do there" is answered only by a scratch script: per
stratum, the AMBIG slots by the RNA on each strand (bands of the minor strand's share), their depth, the gDNA
error, the TILT read-out's error (`|Δτ|·R/2`, RNA fragments on the wrong strand — measured by nothing else) and the
predicted θ peak width. `policy_benchmark.py --by-class` ranks by node class and not by strand split or tilt
error, so the tilt-error column is the new part; a census of what it overlaps comes before it becomes an
instrument (the instrument ruling: few instruments, kept current).

### transfer-variance-premise
`priority: parked — Tier 5 · kind: question · 2026-08`
Does a hop's transfer variance price a ratio built on a handful of counts? The policy prices every hop by both
witnesses' counting plus the pair's disagreement (``hop_price``, `native/transfer_rows.h`); whether that is right where a
pair agrees by coincidence is the open half (the landscape is no substitute: ~10× over-stated). `EQUATIONS.md`
§3.5h.

### drain-contaminates-certified-rna
`priority: parked — Tier 5 (owner, 2026-09-01) · kind: defect · 2026-08-31`
The second pass deposits some true-gDNA fragments into the certified-RNA banks: 233 records at
`g50 ss.99 OFF`, 1,482 at `g98 ss.99 ON` (1.9 % of that in-scope channel). The leak is exactly posterior
sampling and a no-leak counterfactual moves the 0.8.0 metric −0.60 %/−0.05 % at worst, so the harm is the
certainty claim and no in-solve correction is licensed (`ISSUES: drain-provenance-split`). To build: `DrainQC`
records `Σ q_null`; the middle-bin posterior bias repaired as its own A/B. `calibration_vs_oracle.py`.

### crossing-pool-contrast
`priority: parked — Tier 5, blocked on ISSUES: capture-blind-gdna-divisor · kind: question · 2026-08-31`
A second gDNA length contrast on the crossing pools: with oracle weights it beats the contained one under
capture (TV 0.076 vs 0.136; 0.078 vs 0.182) and starves off capture (pool 3 at 29–630 fragments). Blocked: the
weight estimator does not transfer under capture (`ISSUES: capture-blind-gdna-divisor`), and no shadow
transcript overlaps a gene edge, so pool 3 reads exactly 1.0000 pure
(`TRAPS: purity-is-a-property-of-the-annotation`).

### capture-degeneracy-standing-risk
`priority: parked — Tier 5, watch · kind: risk · 2026-08-31`
The gDNA two-pool contrast survives capture by a degeneracy: the shared-contaminant assumption is false under
capture (TV 0.95 vs 0.06–0.14 off) and it is safe only because the intergenic pool is depleted-not-impure, so
`a_0` clips to 1 and the algebra collapses to `g = f_0`. A probe panel that put RNA into intergenic space would
break it silently; `_deconvolved_gdna_counts` carries the derivation. No panel can fire it today. fl's one-sided
off-target rate absorbs evenly spread intergenic RNA, the same break in another estimator
(`ISSUES: the-fl-boundary-inversion-reads-missing-evidence-as-a-value`).

### pure-rna-mirror-asymmetry
`priority: parked — Tier 5 · kind: defect · 2026-08`
Two exact per-fragment mirrors of a pure-RNA library deconvolve differently in `count_gdna_region` by a few
percent, neither boundary-only nor monotone in strandedness. An R1-sense library is simulable; no instrument
yet.

### binary-cuts-on-continuous-quantities
`priority: parked — Tier 5, tracking only (owner, 2026-09-28: one entry; the inventory stays outside the tree) · kind: tracking · 2026-09-28`
The owner's rule is no yes/no gates on continuous quantities; the 2026-09-24 inventory ranks 55 such cuts by harm.
Homed elsewhere: whole counts and the fixed seed; the closed ranks; the od ceiling
(`ISSUES: strand-overdispersion-one-shared-value`); the splice and multimapper cuts. Untracked, with measured harm:
- the reference ρ_ref is the landscape's argmax cell: one cell moves 91–100 % of transcripts' EM length by more than
  1 %;
- `calib_refit_iters = 3` stops unconverged: an L1 change of 85–34k fragments at a fourth refit;
- kernels with count < 1 get no location (the low-depth regime of `ISSUES: the-gdna-landscape-collapses-at-low-depth`);
- pruning on a k·ln10 lattice;
- the zero-density veto, which hits 27–43 % of held records;
- one unit merges two loci;
- `TAIL_DECAY`;
- `max_frag_length = 1000`, in three roles.
Each is raised when its mechanism is next touched, never as a sweep.
Instrument: `~/Downloads/rigel_runs/prototypes/2026-09-24_cut_inventory/INVENTORY.md`.

---

## CLOSED / REFUSED — do not rebuild these; append-only

Every entry keeps its stamped measurement exactly as recorded: a graveyard row without its number is an
invitation to rebuild. A row measured on "all 36 conditions" or quoting `g01`/`g10`/`g25`/`g75`/`g90` predates
the ladder retired 2026-08-13 — the verdict stands as a record, and re-opening one means re-running it on the
current panel. Where a mechanism's only target was unstranded × capture-ON the row is moot as a 0.8.0
candidate on top of being refused; the `g00` zero-control column is never moot.

### benchmark-noise-floors-unmeasured
CLOSED 2026-09-28 (owner): not a project. A last-bit change moves the genes by at most 15 fragments and the gDNA pool
by 113, while the per-transcript Σ|Δ| amplifies it through the EM's start-dependent isoform split (thousands of
fragments: 7,968 on a 4-thread rerun at `g25 ss.50 ON`; `ISSUES: the-em-answer-depends-on-where-it-starts`).
THE RULE: an A/B pair runs with the scan pinned (`--set scan.total_threads=1`), so it is exactly reproducible, and an
effect is judged by its size, with genes and pools read beside the transcript table.

### calibration-vs-oracle-cannot-see-a-ruler-move
CLOSED 2026-09-28: the output that read 0 by construction is deleted — `calibration_vs_oracle.py`'s ruler section and
`prior_vs_oracle.py`'s `gdna_eff_len` score, since the `O` arm keeps `P`'s efficiencies, reference density and gDNA
region lengths, everything the length reads. The length's truth is `ruler_vs_truth.py`; the true-efficiency `O` arm is
`ISSUES: instrument-ledger` (g).

### rna-prior-floor-at-pure-gdna-loci
CLOSED 2026-09-26 (owner) for 0.8.0, to be reopened only with a fundamentally different algorithm. The prize is real —
a perfect calibration prior takes `g98 ss.99 ON` from 25.26 to 19.99 % of the transcripts' true count — but neither
mechanism built for it improves an in-scope stratum without costing another: a no-RNA state read out beside the solve
(`ISSUES: a-no-rna-state-read-beside-the-solve`, refused as a design) and end states inside the solve
(`ISSUES: no-rna-and-no-gdna-end-states-in-psi`, refused at the A/B). What a new approach has to answer, from both:
the weight of "no RNA" at an object must come from evidence whose precision grows with depth, and on unstranded data
an object's strand term carries none (`DESIGN.md` §6b.1); it must not pull objects holding 1–50 % RNA to "no RNA",
where both overshot; and what it changes reaches every contracted length through the capture efficiencies, so it is
judged on the capture-ON rows as well as the pools. The entry as it stood
(`priority: NEXT — the owner's residual pass (2026-09-24): the largest owner of g98's gDNA error · kind: defect · 2026-09-20`):
DISSECTED AGAIN 2026-09-24 on the count form, against per-fragment truth (four lenses, a synthesis and two refuters;
`~/Downloads/rigel_runs/prototypes/2026-09-24_g98_dissection/`). The read-out is the posterior MEDIAN of λ
(`psi_kernel.h`, `DESIGN.md` §6c), not a mean: under the Jeffreys reference with no atom at zero RNA, an object whose
true RNA is below its resolution `1/√(n·I)` reads about 0.42·√(n/I) fragments of RNA, always toward RNA (measured
0.21–0.28·√n stranded, 0.32–0.37·√n unstranded; the stranded/unstranded ratio 0.721 against the Fisher-information
ratio 0.714). It is per OBJECT, and 78–85 % of calibration's count error sits on objects with no RNA, mostly INSIDE
expressed loci; pure-gDNA loci hold only −3.9k / −5.5k / −9.7k. The EM passes the count through at g98 (per-locus
slope 0.90–1.00), and off capture it lands in synthetic spans at introns. Size (the truth restored on no-RNA objects
alone): gDNA −36.7k → −4.7k at `g98 ss.99 OFF`, −61.6k → −12.3k at `g98 ss.50 OFF`; about 47k at `g98 ss.99 ON` once
`ISSUES: the-gdna-length-law-falls-back-at-identical-purities` is separated out. Its transcript cost is small (1.3k
of 1.8k at ss.99 OFF, none at ss.50 OFF). The repair: an atom at zero RNA whose weight is the population's own
no-RNA share, fitted by marginal likelihood across objects (no constant), and a CONTINUOUS read-out — a median of a
posterior carrying an atom is a yes/no cut — which amends the §6c median ruling: the owner's decision. A test bar
must be stated per `n`: at an atom weight of 0.75 the floor falls about 74 % at n = 10 and 89 % at n = 1,000.
The history below measured the pre-count-form RNA pseudocount; its "posterior MEAN" is corrected above.
`rna_prior_count` over-states by +64.5 % at `g98 ss.99 ON` (180,806 against 109,915) and +32.4 % at
`g98 ss.99 OFF`, and it is a diffuse positive floor, not a few loci: 1,089 of 1,149 loci read high and 55k
of the 73k excess sits in loci that are over 99 % gDNA where the true RNA is 544 fragments over 867 loci.
DISSECTED per object (2026-09-20, `g98 ss.99 ON`): on the 14,280 regions holding gDNA and NO RNA the
calibration places 19,538 RNA fragments, on 97 % of them, and the floor is per OBJECT and grows roughly as
the square root of the object's gDNA — a median 0.14 fragments on regions of 1–10 gDNA (11.6 % of their
gDNA), 1.8 at 100–1,000 (1.2 %), 26 at 10,000+ (0.19 %) — the signature of a posterior MEAN of a composition
whose likelihood sits at the boundary `f_r = 0` with a width `~1/√n` (`TRAPS: zero-target-guards-are-one-sided`).
Most of the mass is on BOUNDARIES: 231,410 RNA incidences against 63,972 true, +85,112 on boundaries with
no RNA at all (1.6 % of their gDNA incidences) and +82,326 on mixed ones; regions carry +24,291. ⛔ NOT a
constant fraction and NOT a cheap fix: it is the composition solve's estimator at a pure-gDNA object — a
posterior mean cannot read zero — so the repair is an atom at `f_r = 0` in ψ's hypothesis space (the tilt
already carries {pure +, pure −, mixed}; the composition does not) or a different estimand there, a solver
design item, and it lives at the stress rung (at `g50` the same floor is +4,309 on 519 near-pure loci,
0.15 % of the RNA, invisible in the pool). `prior_vs_oracle.py`, `calibration_vs_oracle.py`; the confidently-wrong
class's z-band table retired with `solvability_audit.py` (2026-09-28). It is why the per-transcript allocation made `g98`
capture-ON worse before the gDNA opportunity was corrected (`ISSUES: per-transcript-prior-lane`).

### a-no-rna-state-read-beside-the-solve
REFUSED 2026-09-26 by the owner, as a design: a second solver on top of the first. The fix belongs in the solve, which
should be able to say "no RNA" itself, and it must leave the tool simpler. Kept as the measurement of the prize of
`ISSUES: rna-prior-floor-at-pure-gdna-loci` (closed). The mechanism
(`~/Downloads/rigel_runs/prototypes/2026-09-24_zero_rna_atom/`, measured through the EM in
`~/Downloads/rigel_runs/prototypes/2026-09-26_zero_atom_em/`): a no-RNA state per object, its weight fitted per
annotation class by marginal likelihood, read as `P0 + (1 − P0)·median`, computed on the pipeline's own final sweep
and handed ONLY to the EM's prior gDNA count; calibration, its landscape and training, and every length as shipped.
All 16 rows, fractional, each against a fresh shipped run and a replicate that gives the row's own floor (0–59
transcript fragments). Transcripts % shipped → with it (perfect calibration prior): stranded OFF 1.988 → 1.979
(1.970), stranded ON 6.203 → 6.150 (6.149), unstranded OFF 2.010 → 2.017 (1.987); genes 0.207 → 0.202, 1.594 → 1.547,
0.238 → 0.249; the deferred stratum 10.28 → 9.80 (9.26). Per row: `g98 ss.99 ON` 25.26 → 20.14 % (19.99), genes
15.03 → 11.06; `g98 ss.99 OFF` 15.14 → 14.72 (14.25); `g50 ss.99 ON` 6.13 → 6.06 (6.07); `g50 ss.99 OFF` −1,262
fragments; every `g00` row's prior unchanged (the fitted weight is 0 there) and its table within 52 transcript
fragments, run-to-run noise that one replicate under-reads; `g05 ss.99 ON` +760 transcript fragments (genes −172);
`g98 ss.50 OFF` WORSE, 16.88 → 17.62 %, genes 5.66 → 6.95 %. The price: it only ever raises gDNA and it overshoots on
objects holding 1–50 % RNA (in fragments, `g50 ss.99 ON` 1–10 % band |error| 40.7k → 50.6k, 10–50 % 38.9k → 44.5k;
`g98 ss.99 ON` signed +4.5k → +26.1k and +1.1k → +3.8k), which lands on the synthetic spans (est ÷ true
`g98 ss.99 OFF` 1.88 → 0.36, `g98 ss.50 OFF` 2.52 → 0.27, `g50 ss.99 ON` 0.91 → 0.77). On unstranded data the strand
term cancels and the weight is read from the density terms alone, which is where it loses (capture OFF; the deferred
unstranded capture-ON stratum gains). Unmeasured: realistic nascent (on_fraction 0.10) at `g98`.

### no-rna-and-no-gdna-end-states-in-psi
REFUSED 2026-09-26 at the A/B; the owner closed the thread for 0.8.0. The owner's fix at the source: ψ's gDNA-fraction
axis given its two end states — "no RNA" (`f_g = 1`) and "no gDNA" (`f_g = 0`) — beside its continuum, each at the
continuum's reference mass, as the strand axis carries "all RNA on +" and "all RNA on −" (`DESIGN.md` §6b.15.13); in
the final solve of every sweep only; read out as each state's answer weighted by its posterior, the §6c median inside
the continuum. A C++ prototype in a worktree, switched off bit-identical to shipped on all 30 calibration arrays
through three refits at `g05 ss.50 OFF` and `g98 ss.99 ON`
(`~/Downloads/rigel_runs/prototypes/2026-09-26_presence_states/`: `DERIVATION.md`, `README.md`, and
`presence_states.patch`, the diff against `062bc9ea`).
Calibration against the oracle (region + boundary |gDNA error|, `g05`–`g98`): stranded OFF −14.8 %, stranded ON
−10.1 %, unstranded OFF −12.1 %, but `g05 ss.99 ON` +6.3 % (80.7k → 85.8k), `g50 ss.99 ON` +1.4 % and `g50 ss.50 OFF`
+1.2 % (the deferred `g05 ss.50 ON` +94.7 %); the `g00` controls stranded 353 → 80 (OFF) and 329 → 44 (ON), unstranded
366 → 552 (OFF) and 266 → 299 (the deferred ON); the "no RNA" state alone takes those two unstranded `g00` rows to
74,442 and 70,536.
Through the EM (fractional; transcripts / genes as a % of truth per stratum, shipped → this, the perfect calibration
prior in brackets): stranded OFF 1.988 → 1.975 (1.970) / 0.207 → 0.202; stranded ON 6.203 → 7.026 (6.149) / 1.594 →
1.778; unstranded OFF 2.010 → 2.013 (1.987) / 0.238 → 0.245; the deferred stratum 10.28 → 18.88. Per row:
`g98 ss.99 ON` 25.26 → 21.73 % (19.99), genes 15.03 → 11.53; `g05 ss.99 ON` 5.56 → 7.78 %, genes 0.59 → 1.16 (+204k
transcript fragments against a replicate floor of 25); the deferred `g05 ss.50 ON` 11.03 → 34.19 %.
Why it fails. (1) The transcript table's harm runs through the capture efficiencies, and so through every contracted
length: with the efficiencies and the capture reference held at the switched-off values, `g05 ss.99 ON`'s +204k
transcript fragments fall to −13 (floor 25), the deferred `g05 ss.50 ON`'s +2.13M to −230, and `g98 ss.99 ON` reads
21.42 %; the pools do not all return — at the deferred `g05 ss.50 ON` the counts alone move gDNA into the synthetic
spans (gDNA 391k → 251k, truth 500k; spans 380k → 520k, truth 285k). At `g05 ss.99 ON` 1,635 objects' efficiencies
move by more than 0.05 (up to 0.75), and objects holding 10–50 % RNA (median count 33) take P(no RNA) > ½ 31 % of the
time. (2) Three states at equal weight give each end even prior odds against the whole continuum, above the no-RNA
share the record fitted by marginal likelihood (`ISSUES: a-no-rna-state-read-beside-the-solve`) on every class at
`g00` and `g05` and at `g50` OFF, and on all but the intron classes at `g50 ss.99 ON`: 0 at `g00`, ≤ 0.024 at `g05`
OFF, 0–0.497 at `g05 ss.99 ON`, ≤ 0.43 at `g50` OFF. (3) On unstranded data the strand term is ½ in every state, so
the states are decided by terms whose precision does not grow with depth — the landscape's density and one-sided level
rows — and the fixed weight is the answer at every depth, the pattern `DESIGN.md` §6b.1 refuses. (4) A one-sided row
read at its held end favours the end state it does not bound, and the "no gDNA" state reads the landscape's density at
its floor cell as the mass of a point (`TRAPS: a-ratio-cannot-carry-zero`). (5) The published `f_g` mixes the states
while its `Var(log f_g)` stays the continuum's, so the landscape trains on a pair that no longer matches. Not built:
the weight learned per class in the refit loop — the fitted shares are near 0 where the states help least and reach
~0.5 on the intron boundaries at `g05 ss.99 ON`, where the harm is.

### the-gdna-length-law-falls-back-at-identical-purities
CLOSED 2026-09-26 by the fix, in the working tree for the owner's commit (`calibration/fl.py`,
`_deconvolved_gdna_counts`): at a purity separation of exactly 0 the contained pair answers with its own de-tilted
mixture — the limit its resolution-weighted fade approaches, and what the boundary pair already returned at its own
tie — instead of declining to the capture-selected four-pool census. MEASURED 2026-09-26 (no EM,
`~/Downloads/rigel_runs/prototypes/2026-09-26_gdna_law_tie/`): only `g98 ss.99 ON` ties on the ladder (both purities
1.000); there the uniform-frame law reads 245.1 → 217.2 bp (truth 216.7) and the realized law 254.3 → 234.9 (truth
240.9, where every other capture-ON row reads 234.5–236.5 against ~241: the census's own offset,
`ISSUES: the-scorer-reads-a-census-length-law`); the other 15 rows' laws are bit-identical, so nothing else moves.
Through the EM (fractional, `g98 ss.99 ON`), transcripts / genes / gDNA est − true / synthetic spans est (truth 6.0k)
/ annotated est − true: 26.18 → 25.26 % / 17.60 → 15.03 % / −281.0k → −246.5k / 79.7k → 53.7k / +26.9k → +18.2k; under
the oracle prior 21.34 → 20.04 % / 12.69 → 9.42 %; `--arm base_reseed` equals `base`. The landed pipeline reproduces
the prototype to 0.24 of 49,013 transcript fragments. The 2026-09-24 reading that transcripts worsened (25.95 →
26.46 %) does not reproduce on this tree. Gates (`tests/calibration/test_fl.py`), verified failing on the shipped
code: the tie at purity 1, the uniform law continuous through the tie end to end, and a tie below purity 1;
restoring the decline fires all three, deleting the tie branch fires the two `errstate` guards (bin `L = 0` is
empty in both pools, so `0/0` makes the inverted law's sum NaN and the old arithmetic reached the mixture by
accident). Left as it was: the contained pair still declines on an empty pool — the `g00` rows, which then read the
four-pool census (243–244 bp under capture) — the same kind of cut, a separate mechanism. The entry as it was opened:
`calibration/fl.py`'s `build_fl_models` estimates gDNA's fragment-length law by contrasting two pools of different
gDNA purity; when the separation is EXACTLY 0.0 ("purities identical", as at `g98 ss.99 ON`, both pools at gDNA
share 1.0) it declines and falls back to the capture-selected four-pool census: `gdna_pmf` reads 245.1 bp against a
true 216.7. Where the contrast is applied the estimator is right (216.7–217.0 at `g50 ss.99 ON` and `g98 ss.50 ON`,
separations 0.021 and −0.0018). A yes/no fallback on a continuous quantity. At `g98 ss.99 ON` this one root owns the
scorer's gDNA law reading +13 bp long (254.3 against 240.9), the gDNA component's offset from the RNA lengths'
scale (gDNA/spans 1.014, gDNA/isoforms 1.041, which fall to g50's levels with the true law), and 13.6k of
calibration's count deficit. Fixed alone: gDNA −100.5k → −67.7k, spans 80.0k → 54.1k, the oracle arm −35.5k →
−21.7k; but transcripts do not improve (25.95 → 26.46 %): correcting the law the contracted lengths use exposes an
opposite error in the length rule, so the repair is judged with its effect on the lengths. Class-mean L/Y moving
toward one scale is NOT a valid gate for it: that moved toward one scale while the EM got worse (−4.8k to −5.9k
gDNA, +6.2k to +9.8k transcripts). Scripts: `analysis/refute-mechanism/fl_root.py`, `scale_arm.py` in `~/Downloads/rigel_runs/prototypes/2026-09-24_g98_dissection/`.

### the-capture-length-owns-stranded-capture-on
CLOSED 2026-09-26 (owner), for now: the tool is at diminishing returns, a change must be simple and improve it, and
the shipped length stays. The frame as it was opened:
revisited with a new argument and a measurement · kind: defect · 2026-09-26`
Stranded × capture-ON is the worst in-scope stratum: transcript Σ|Δ| 6.21 % of truth against 1.99 % stranded OFF
and 2.01 % unstranded OFF (the ladder, fractional, arms of 2026-09-25). Handing the EM the simulator's own
capture-aware lengths for every transcript and synthetic span (`quant_accuracy.py --arm oracle_ruler`; the gDNA
component's length and the priors stay as shipped) takes it to 2.02 %, while a perfect calibration prior
(`--arm oracle`) leaves it at 6.16 %: the length the EM divides each component's reads by under capture
(`DESIGN.md` §7.2, `EQUATIONS.md` §11) owns about two-thirds of the stratum's transcript error. Per condition,
transcripts (genes) as a share of truth, shipped → true lengths: `g00` 6.47 → 1.42 % (2.36 → 0.19 %), `g05`
5.56 → 1.72 % (0.59 → 0.52 %), `g50` 6.13 → 2.84 % (1.43 → 1.37 %), `g98` 26.18 → 25.90 % (17.60 → 18.65 %). At
`g05` and `g50` the length decides the isoform split inside genes; at `g00`, with no gDNA to witness capture, it is
nearly the whole error and moves the split between annotated transcripts and synthetic spans as well; at `g98` it
is not the owner. The capture-OFF rows do not move, as they must not, and the deferred stratum moves the same way
(`g05 ss.50 ON` 11.03 → 1.46 %). A capability proof, never headroom: it hands over the true lengths.
Where it sits, as far as measured: at `g50 ss.99 ON`, 161k of the 168k within-gene error that lengths on their
yield remove sits in genes with a junction-probed isoform (`ISSUES: the-junction-price-is-noisy-within-a-gene`,
parked); the short exon pieces of unprobed transcripts read above their capture truth
(`ISSUES: the-efficiency-posterior-floor-on-empty-pieces`); the reference needs gDNA (`DESIGN.md` §7.2). Not
measured: the transcript cost of each class's length error alone. The per-class length error against the truth is
`ruler_vs_truth.py --scale`; the truth is the simulator's `CaptureSampler.partition_array`; the deliverable is
`quant_accuracy.py --arm base` beside `--arm oracle_ruler`, per stratum, under
`--set em.assignment_mode=fractional`.
WHAT THE CAMPAIGN MEASURED (2026-09-26, `~/Downloads/rigel_runs/prototypes/2026-09-26_capture_length/README.md`;
fractional, stranded × capture ON, transcripts / genes as a % of the true annotated RNA at `g00` / `g05` / `g50`):
* The frame reproduced exactly (base 6.47 / 5.56 / 6.13 / 26.18 with `g98`, `--arm oracle_ruler` 1.42 / 1.72 /
  2.84 / 25.90); under fractional assignment `--arm base_reseed` equals `base` to 0.01 point.
* `--arm oracle_ruler` is not a clean ceiling: it keeps the gDNA component's shipped length, so at `g50` the
  synthetic spans read 3.25× their truth and the gDNA pool −459k. With the gDNA component at its own true yield on
  the same scale the ceiling reads 1.44 / 1.45 / 2.08 (genes 0.19 / 0.23 / 0.58) with clean pools.
* Where it sits: at `g05` / `g50` about 80 % of what the true lengths remove is the isoform split inside probed
  multi-exon genes (within-gene error in genes with a junction-probed isoform 4.61 → 0.96 / 4.43 → 1.19 points);
  at `g00`, with no reference and nothing contracted, half of it is the spans against mRNA (gene-level 2.36 → 0.19).
  Probed genes without a junction probe, unprobed genes and the classes' mean offsets cost almost nothing, and the
  shipped spans already sit on the gDNA component's scale.
* It is the junction price's STRUCTURE, not its noise: with every other length true, the shipped sum on the
  simulator's noise-free gDNA weights alone reads 5.11 / 5.12 / 7.10 — the whole shipped error. On the simulator's
  own capture landscape (4,264 probed multi-exon transcripts) fragments across a junction carry 46 % of a probed
  transcript's yield, lie wholly in exons and are captured nearly uniformly (weight median 1,051 of a full probe's
  1,201; 1,092 where a probe spans the junction, 906 where none does; ~1.19× the reference), while a gDNA fragment
  across a junction's boundary lies half in the intron: noise-free the sum is off by a log sd of 0.64 per junction.
* The isoform split needs within-gene length ratios near 1 %: with every other length true, one true price per gene
  leaves 2.57 / 2.45 / 3.23, one price for every junction 3.23 / 2.98 / 4.05, two prices by whether a probe spans
  the junction 2.65 / 2.38 / 3.41.
* Swapping one class's true lengths into an otherwise true set overstates that class's cost — a gene's span and its
  isoforms are priced from the same noisy objects and their errors cancel (each class alone cost about the whole
  base error at `g50`) — so a class is also attributed from the shipped side, with every gene's own scale kept.
Refused on the way: `ISSUES: a-junction-price-from-its-genes-own-rna`,
`ISSUES: a-junction-price-at-the-best-object-its-fragments-reach`. What stays open beside it:
`ISSUES: the-efficiency-posterior-floor-on-empty-pieces`, `ISSUES: ruler-witness-geometry-on-transcript-panels`.

### the-junction-price-is-noisy-within-a-gene
CLOSED 2026-09-26 (owner), for now, with `ISSUES: the-capture-length-owns-stranded-capture-on`: the sum stays. That
campaign found the within-gene error to be the rule's structure — a fragment across a junction lies wholly in exons
and the gDNA beside it half in the intron — and refused the two prices tried against it. The record as it stood
(PARKED 2026-09-23; the pseudocount's odds fixed 2026-09-24):
SIZED 2026-09-24 (the g98 dissection's g50 contrast, `~/Downloads/rigel_runs/prototypes/2026-09-24_g98_dissection/`): with every length on its
yield at its locus's gDNA scale, the within-gene transcript error at `g50 ss.99 ON` falls 231.6k → 64.0k
(transcripts 6.14 → 1.74 %), and 161k of the 168k sits in genes with a junction-probed isoform; at `g98 ss.99 ON`
it is about 3.5k. Still parked: nothing derivable moves the noise-free ceiling.
The one shared length rule (`DESIGN.md` §7.2, `EQUATIONS.md` §11) prices a junction by conservation of bases,
`c_junction = c_lo + c_hi − ½(c_intron,lo + c_intron,hi)` (`capture_eff_length._cut_efficiencies`). It is right on
average — it is what puts the transcripts on the gDNA component's scale — but it widens the within-gene spread of
the transcripts' lengths, and the EM splits isoforms by that spread. The within-gene sd of log(L / Y) about each
gene's mean (annotated multi-exon transcripts with ≥ 20 fragments, no EM) beside the class mean L / Y ×10⁻³,
`g05` / `g50 ss.99 ON`:

| the junction priced by | within-gene sd | class mean |
|---|---|---|
| the per-base rule, before the port (`ISSUES: the-per-base-rule-for-every-component`) | 0.052 / 0.040 | 0.944 / 0.957 |
| its adjacent pieces (`ISSUES: the-junction-price-at-its-neighbours-mean`) | 0.057 / 0.041 | 0.936 / 0.965 |
| the sum, SHIPPED | 0.072 / 0.057 | 1.076 / 1.112 |

On the shipped tree (`ruler_vs_truth.py --scale`, width step 1) the gDNA component and the synthetic spans read
1.108 / 1.118 at `g05` and 1.131 / 1.141 at `g50`, 3.9 % / 2.6 % from one scale with the isoforms. Giving back the
scale is not free: the synthetic pool's elasticity to its own length is −21 at `g50 ON` (2026-09-21).
MEASURED 2026-09-23 — THE SPREAD IS THE RULE'S, NOT THE COUNTS'. The same rule fed four count sources, within-gene
sd (the shipped row is `--scale`'s, reproduced exactly):

| every region's and boundary's gDNA count | `g05` | `g50` |
|---|---|---|
| the simulator's expected count, as a plug-in — no noise, no prior | 0.049 | 0.053 |
| the same counts through the posterior | 0.054 | 0.053 |
| the true realized counts (`slot_truth.npz`) | 0.059 | 0.054 |
| calibration's, SHIPPED | 0.072 | 0.057 |

Swapping one term at a time (on the transcripts' own routed shares, where the shipped rule reads 0.074 / 0.058): the
junction at the transcript's own truth 0.048 / 0.030; all four junction terms noise-free 0.065 / 0.056 — the ceiling
of any price that only removes noise; the intron terms alone noise-free 0.075 / 0.060, no gain. Per junction the
noise-free rule is 0.27 / 0.31 rms from the truth on prices near 1, and the noise adds 0.23 / 0.13 in quadrature.
The rule's own error is the capture physics, which no gDNA object sees: the simulator binds a fragment through its best
single probe part, so a junction a probe spans (one part in the transcript, two in gDNA) is priced exactly by the
sum, while a junction with separate probes on its two sides binds one of them and the sum prices both — and the two
read the same to every region and boundary (`ISSUES: ruler-witness-geometry-on-transcript-panels`). Where both sides
carry capture the sum over-prices by 18 % (truth / sum 0.847); where one side dominates it under-prices by 25 % (each
boundary is clipped at 1, a junction reaches 1.3); the per-junction median error is 19 % whatever the exon lengths.
Second, the geometry: exons beside junctions are short (median 128 bp, 93 % below the longest fragment), so a
fragment across a junction holds bases of up to four exons where the gDNA crossings at its boundaries hold exon ends
and intron; under additive capture the rule is exact where both exons exceed every fragment (median error 3.6 %)
and 15 % off elsewhere. The count noise is Poisson on ~37 / ~340 expected crossings at a fully captured boundary,
plus calibration's bias at mid-density boundaries (+0.10 in the 0.3–0.7 band at `g05`); the posterior neither adds
nor removes it (posterior and plug-in on the same counts agree to 0.005 rms in every band). The unprobed class's
inflation at `g05` is not the junction's (`ISSUES: the-efficiency-posterior-floor-on-empty-pieces`). Refused on the
way: `ISSUES: a-junction-price-read-from-a-fitted-capture-map`, `ISSUES: the-junction-sum-on-unclipped-posteriors`.
The census, the attribution and the derivation: `~/Downloads/rigel_runs/prototypes/2026-09-23_junction_price/`.
THROUGH THE EM (`panel.py score` on the shipped tree, 16 ladder conditions, fractional assignment, stranded ×
capture ON): transcript and gene Σ|Δ| as a % of the true annotated RNA, and the pools as estimate − truth — gDNA
(intergenic included) / synthetic spans est/true / annotated. The floor is |base − base_reseed|; the third column
adds the population-corrected RNA pseudocount of `ISSUES: the-pseudocount-prior-is-biased-toward-gdna` (not
shipped; the prototype's run of 2026-09-23, which the shipped column reproduces within 0.03 points without it):

| | before the port | the shipped rule | + the corrected pseudocount |
|---|---|---|---|
| `g05` transcripts (floor 0.009) | 3.46 | 5.59 | 5.55 |
| `g50` transcripts (floor 0.003) | 5.50 | 6.18 | 6.20 |
| `g98` transcripts (floor 0.008) | 31.05 | 26.24 | 28.00 |
| `g05` / `g50` / `g98` genes | 0.44 / 1.91 / 16.87 | 0.64 / 2.12 / 14.96 | 0.62 / 1.45 / 18.54 |
| `g05` pools | +22,538 / 0.97 / −15,185 | +76,880 / 0.76 / −8,503 | +3,187 / 0.96 / +7,001 |
| `g50` pools | +35,572 / 1.22 / −68,501 | +172,932 / 0.28 / −64,052 | −25,216 / 0.98 / +28,209 |
| `g98` pools | −122,373 / 22.41 / −5,839 | −54,229 / 8.64 / +8,481 | −117,852 / 15.97 / +28,218 |

At `g05` the transcript error rises 2.13 points against 0.20 for the genes on the shipped rule (2.09 / 0.18 with
the corrected pseudocount): it is isoform allocation. Capture OFF the transcripts and genes move within ±0.03
points except at `g98` (unstranded 18.36 → 18.23, stranded 15.70 → 15.53); the deferred stratum's `g98`
transcripts read 90.63 → 115.07 (116.50 with the corrected pseudocount). The pools improve under the corrected
pseudocount while the isoforms regress either way (`TRAPS: judge-a-ruler-by-its-within-gene-spread`). The ceiling
it is priced against: the simulator's own lengths for every hypothesis with the corrected pseudocount read 1.89 %
and the synthetic pool 0.96× at `g50 ss.99 ON` (2026-09-21).

STILL OPEN FROM THE PORT'S REVIEW (2026-09-23, `g05 ss.99 ON` unless stated; recorded, not fixed):
* THE SHORT-INTRON FALLBACK keys on `S > 0` exactly (`S = gdna_region_eff_len`, the intron piece's gDNA contained
  support), so a piece with 0 < S < 1 (less than one expected start) holds essentially no information yet reads its
  own posterior (≈ the population's mean, 0.244, sd 0.058). Beside junctions 3,057 / 2,724 intron pieces (low /
  high side) have S = 0 and 2,759 / 2,819 have 0 < S < 1; treating S < 1 as no count moves 4,007 of 45,609
  junction prices (1,284 by more than 0.2, at most 0.79). A real library's shorter fragment tail shrinks the S = 0
  set (7,275 → 4,075 regions with a tail to 20 bp). No threshold is proposed: a threshold is a constant. MEASURED
  2026-09-23: beside 15–20 % of scored junctions the intron piece holds under 10 expected starts and the posterior
  moves the price 0.07–0.13 from its noise-free value, yet the intron terms made exactly noise-free leave the
  within-gene sd where it was (above) — a per-junction defect with no isoform cost on the ladder.
* THE CLIP AT 0: 137 transcripts carry more than half their share on junctions priced exactly 0; their contraction
  factor has median 0.0033 (minimum 0.00069), 30–130× below the adjacent-piece price's.
* ABOVE fl AND ABOVE THE PARENT: 2,593 of 8,750 annotated transcripts have `eff_em > fl` (at most 1.80×), and 23 of
  the 7,709 with a parent synthetic span read longer than it (median 1.010, at most 1.150) — possible with junction
  probes, no longer excluded by construction.
* gDNA'S REACH AT A REFERENCE END: gDNA's conserved share takes `UNBOUNDED_REACH` at every boundary, including one
  within a fragment of a reference end, where the accumulator clips (a boundary 10 bp from a reference start: true
  share 10.0, shipped 129.5) — gDNA's reach ruling (`region_geometry`: gDNA takes no reach argument), consequential
  only on short contigs.

THE CONSTRAINTS: the owner's rulings in `DESIGN.md` §7.2 (capture is local; junctions never pooled,
`ISSUES: pooling-junctions`; capture stays outside the EM) and §0b (robustness); no panel input (2026-09-20).
Refused on the way: `ISSUES: the-spliced-read-junction-price`, `ISSUES: a-junction-price-clipped-at-one`,
`ISSUES: a-junction-price-read-from-a-fitted-capture-map`, `ISSUES: the-junction-sum-on-unclipped-posteriors`. A
candidate is read on both columns of the first table at once — the within-gene sd toward the adjacent pieces', the
class mean on the gDNA component's — and only then through the EM; a candidate that only removes noise cannot pass
0.065 / 0.056. What would move it is an observable of the probe physics. `ruler_vs_truth.py --scale` (its
within-gene sd line and the junction-probed rows); `quant_accuracy.py` per stratum under
`--set em.assignment_mode=fractional`.

### a-junction-price-from-its-genes-own-rna
REFUSED 2026-09-26 by the owner, before it was built: one junction price per gene, read from the gene's own spliced
against its unspliced mature RNA. Its noise-free ceiling (the true per-gene price, every other length true) read
2.57 / 2.45 / 3.23 % transcripts at `g00` / `g05` / `g50 ss.99 ON` against the sum's 5.11 / 5.12 / 7.10. Why not: a
gene's isoforms start, end, include and skip a junction independently, so reading capture along a transcript needs
its abundance, which is the EM's output; capture efficiency stays outside the EM (`DESIGN.md` §7.2), and fitting the
two at once is not a narrow change. It would also relax the locality and no-pooling rulings. Do not re-propose a
junction price read from the RNA.

### a-junction-price-at-the-best-object-its-fragments-reach
REFUSED 2026-09-26 at the A/B. The owner's narrow candidate: a fragment across a junction binds its best single
probe part, so the junction takes the largest published efficiency among the objects its fragments touch in the
transcript's own coordinates — the junction's two boundaries, the exon pieces beside it and, where an exon is
shorter than a fragment, the next exons' boundaries — averaged over its placements with the deposit rule's weights;
no constant (`junction_reach/` in `~/Downloads/rigel_runs/prototypes/2026-09-26_capture_length/`). No EM, class-mean
L / Y ×10⁻³ and within-gene sd, shipped → this: `g05 ss.99 ON` isoforms 1.076 → 0.953 (spans / isoforms 1.039 →
1.173), within-gene 0.072 → 0.081; `g50` 1.112 → 0.989 (1.026 → 1.153), 0.057 → 0.049; the unprobed isoforms 2.84 →
59.7 / 1.21 → 21.8. Through the EM (fractional, stranded ON), transcripts / genes / spans est ÷ true: `g05` 5.56 →
4.16 / 0.59 → 1.13 / 0.98 → 0.99; `g50` 6.13 → 5.32 / 1.43 → 2.35 / 0.91 → 0.69; `g98` 26.18 → 27.29 / 17.60 → 19.48 /
13.30 → 12.88; `g00` and capture OFF unchanged (no reference). The isoform split improves while genes, pools and `g98`
regress, through the scale: every published efficiency is clipped at 1 and a fragment across a junction is captured
~1.19× the reference, so the best neighbour reads low and the isoforms fall 11–17 % below the spans and the gDNA
component. Refused on the way, noise-free: each edge read as an exon level (`2·c_edge − c_intron`) under the max
overshoots (class mean 1.20–1.31, within-gene sd 0.085–0.114 against the sum's 0.060); the boundary crossings split
by which side holds the fragment's midpoint price a junction no better than the sum (per-junction log sd 0.58
against 0.64), so no accumulator change is worth it.

### spike-in-references-read-as-depleted
CLOSED 2026-09-26 (owner), for now. A spike-in (ERCC) reference holds no gDNA template, so under capture its objects
read the landscape's depleted level and a probed spike-in's published `em_effective_length` is ~800× short
(ERCC-00131: 0.7 bp at `g05 ss.99 ON`). Its COUNT is right — Σ|Δ| 14.8 / 16.9 fragments of 24,747 / 12,982 at `g05` /
`g50 ss.99 ON`, the same with the true lengths — because the spike-in and its reference's gDNA component are
contracted alike. "No gDNA on the reference" cannot recognise one without a threshold: calibration places gDNA on
the ERCC contigs, up to 5,007 fragments (1,366 on one contig) at `g50 ss.50 ON` and 11.9 on one contig at `g05 ss.99
ON` (`ercc/` in `~/Downloads/rigel_runs/prototypes/2026-09-26_capture_length/`). If reopened: a reference covered end
to end by one single-exon transcript is recognisable from the index alone, with no constant; the gain is the
published yield only.

### the-em-answer-depends-on-where-it-starts
RULED 2026-09-26 (owner): accepted as a known limit and documented (`DESIGN.md` §3.1e; `MANUAL.md`, the FAQ on a
component that dies). The shipped start stays the coverage-weighted one, the most accurate measured.
What is left after the SQUAREM repair (`ISSUES: the-squarem-clamp-decided-which-components-live`, CLOSED): VBEM's
own start-dependence. With the clamp gone MAP reaches one answer from any start (204 fragments apart at
`g00 ss.99 ON`, all in loci still at the iteration cap), while VBEM still ends 8,096 fragments apart there (136
loci, 132 of them winner forks — a candidate dead in one answer and holding fragments in the other). Its E-step
weight ψ(α), with no per-component prior, penalises a small component by about −1/α, so components sharing fragments
race and the start picks the winner; the two answers' likelihoods differ (the uniform start's is higher in 121 of 136
loci), so it is not a flat valley. Candidate repairs, each measured or excluded below: a warm start that is itself
start-free (the MAP optimum, which is unique), or a per-component prior. Instrument: `em_start_lab.py` in
`~/Downloads/rigel_runs/prototypes/2026-09-25_em_start/` (every EM setting re-solved from one pre-EM state in one
process, which repeats itself to 0.0 fragments; `analyze.py` classifies each moved locus). The minus-strand
coverage-weight fix, which seeds only the start, landed on its own (`ISSUES:
the-minus-strand-coverage-weight-sat-one-fragment-length-off`, CLOSED); what it moved is this entry speaking.
MEASURED 2026-09-25 on `g00 ss.99 ON` (the landed solver's prototype): the 136 loci VBEM's two starts disagree on carry
EQUAL truth error — 385,299 (coverage start) against 385,088 (uniform), MAP's single answer 385,376; the coverage start
is closer in 50 loci, the uniform in 64, 22 tie — so what is left is reproducibility, not accuracy: the forks pick
among near-equivalent explanations of genuinely ambiguous fragments.
TWO START-FREE ARMS MEASURED 2026-09-25 (`~/Downloads/rigel_runs/prototypes/2026-09-25_vbem_start/`; the ladder
against the landed build, all 16, fractional, `ladder/compare_vbem.txt`): (1) VBEM STARTED FROM THE MAP OPTIMUM
(prototype `~/proj/rigel-em`, `RIGEL_PROTO_VBEM_START=map`: MAP first, then VBEM from the counts one E-step at θ_MAP
implies) — nearly start-free (the two starts 294 fragments apart at `g00 ss.99 ON` against 8,096; truth error there
936.6k–936.9k against 937.8k–938.0k), in-scope transcript Σ|Δ| −4.7 / +0.2 / −3.8 % (stranded OFF / stranded ON /
unstranded OFF), but genes +0.6 / +1.0 / +4.0 % and the gDNA pool's Σ|error| +3.3 / +2.1 / +3.2 %: it inherits MAP's
signature, more mass on annotated transcripts and less on the synthetic spans and on gDNA (`g50 ss.99 ON` −13.4k →
−15.5k), and it costs 30–50 % more pipeline time on the heavy conditions (measured beside the uniform arm under the
same load). (2) VBEM FROM THE UNIFORM START (a config value, `em.warm_start=uniform`): worse everywhere, transcripts
+2.0 / +0.1 / +0.5 %, genes +1.6 / +0.4 / +1.8 %. So which basin VBEM starts in is not neutral — the coverage start's
favour genes and gDNA, the MAP optimum's favour the transcripts off capture — and neither arm is a clean win. SIZED on an
unstranded row (`g05 ss.50 OFF`, the same lab): the landed VBEM's two starts end 115,159 fragments apart (134 loci, 58
by more than 100), and there the coverage start is the better one — truth error 1,371.7k against 1,376.9k — though the
uniform start is closer in more of the forked loci (76 against 44). There VBEM from the MAP optimum is neither
start-free nor better: its two starts end 19,404 fragments apart, and its truth error is 1,413.8k / 1,410.6k against
the landed 1,371.7k (+3 %) — it moves mass off the synthetic spans and gDNA onto the annotated transcripts, which is
why the transcript table alone read better. The coverage start has been the most accurate start measured on every
condition. Not
tried: choosing among VBEM's optima by VBEM's own objective (the evidence lower bound), which needs its derivation
under the grouped prior. A per-component reference prior is excluded by `EQUATIONS.md` §9b.1 (any lift above the
~0.16–0.47 activation threshold revives every shadow entity; the equal share was refused 2026-09-19).

### whole-counts-optimised-assignment
REFUSED 2026-09-26 (owner): the draw is kept. An assignment chosen for the most reads on their true origin under the
same counts wins reads and loses where they came from. The owner's framing had been: the draw a WARM START, the
repair a FEASIBILITY step (count error zero, `EQUATIONS.md` 13.3–13.4), and what is left an OPTIMISATION over
feasible assignments, wanted with a proof of what it reaches and a stated cost. The exact optimum was derived, built
and priced (`EQUATIONS.md` §13.5: prices, then shortest augmenting paths; exact by construction, certified by its own
dual). Prototype in
`~/Downloads/rigel_runs/prototypes/2026-09-26_whole_counts/` (`opt.c`, `exact.py`, `fidelity.py`): 9,000 random
loci equal brute-force enumeration, and dropping or corrupting the price update, or moving the wrong fragment, fails
that gate. MEASURED on two full libraries, the counts identical under every method (`g05 ss.99`, capture OFF / ON,
9.2 M fragments each):

| | reads on their true origin | reads misplaced by class | time |
|---|---|---|---|
| the fractional posterior (the model's own floor) | 66.53 / 72.41 % (expected) | 8.13 / 4.86 % | — |
| the draw (ships) | 66.53 / 72.43 % | 8.45 / 5.01 % | 0.4 / 0.2 s |
| the optimum, most reads expected correct (`Σ p`) | 68.09 / 74.07 % | 24.20 / 18.38 % | 130 / 7 s |
| the optimum, most probable (`Σ log p`) | 67.44 / 73.14 % | 16.09 / 11.29 % | 106 / 5 s |

"Misplaced by class" is `½ Σ |assigned − true|` over (candidate set, component) cells: whether each component's reads
come from where its reads really came from. An optimum is a vertex of the transportation polytope, so it hands groups
of look-alike fragments to one component instead of sharing them out — its gain and its cost are the same act. The
draw is the posterior sampled with the counts held, within 0.15–0.32 points of the model's own floor. The 130 s is
one locus (737 k fragments over 521 components, a median 47 effective candidates each: ten Gauss–Seidel sweeps still
left 334 k excess, then 110 s of shortest paths); a two-scale solve (the candidate-set classes first, then the
fragments) is designed and not built. REFUSED on the way: the auction algorithm with ε-scaling, which re-bid all
9.2 M fragments at each of six ε stages while near-identical fragments raised prices by ε per bid; it was stopped
after 25 minutes without finishing its first library. Greedy shortcuts do not substitute either: a greedy quota fill
strands 52–95 k fragments per library, and taking the largest re-weighted candidate instead of drawing strands
23–28 k. The draw's shortfall per read is the price of calibration, and an annotated BAM is read for where each
transcript's reads are.

### whole-counts-by-rounding-then-assigning
FIXED 2026-09-25 (`DESIGN.md` §3.1f; `em_solver.cpp`'s `assign_posteriors`; gates `tests/test_estimator.py`'s
`TestWholeCounts`, 14 cases). Whole-count mode gave each fragment to one transcript, drawn from its posterior after every
candidate under 1 % of THAT fragment was zeroed, to keep noise isoforms from collecting whole counts — but on one
fragment a real minority and noise look alike, so the floor took a minor isoform, a gene's unspliced RNA or low gDNA off
every fragment it shared (off capture +24 to +35 % transcript and 2.3–3.4× gene error against fractional). Now COUNT
FIRST: every component's fractional count, the locus's gDNA included, is rounded within its EM locus by largest
remainder; each fragment is drawn from its own posterior re-weighted by (count still owed ÷ posterior mass still to
come); a fragment the draw leaves with every candidate full is repaired onto the counts by a breadth-first chain of
moves, each fragment only to another of its own candidates. The whole counts ARE the rounded fractional counts,
whatever the seed, whenever that rounding is reachable; where it is not (Hall's condition fails — two components one
fragment can reach, both rounded up), a fallback of the same augmenting paths settles every count within one of its
fractional count, which is always reachable. The derivation — termination, the success condition, the cost, and the
optimality it gives up — is `EQUATIONS.md` §13. The fallback was added after a search found the exact-target repair
alone leaving a count two from its expectation; on 30,000 draws over 6,000 small loci built around an unreachable
core (exact targets unreachable in 9,310), no count left its floor or ceiling and no locus total moved.
`assignment_mode="map"` and `assignment_min_posterior` are deleted (owner). The gates hold the counts
to an independent largest-remainder reference on four loci (one ordered so the draw strands and the repair must work;
one where every remainder is under a half, so rounding each count alone would lose a fragment), a minority under 1 % on
every fragment keeping its count, a transcript expecting under half a fragment getting none, seed-independent counts
with seed-dependent fragments, every fragment on one of its own candidates, and — for identical fragments — the first
and the last drawn with the same odds (the draw is an urn). Each of six perturbations fires its own gates: the old
floor, no repair, rounding each count alone, no re-weighting, no randomness, and a lost per-read output. LADDER, the
shipped whole counts against fractional (all 16, `~/Downloads/rigel_runs/prototypes/2026-09-26_whole_counts/ladder/`):
transcript Σ|Δ| −0.02 / +0.04 / −0.01 % (stranded OFF / stranded ON / unstranded OFF), genes −0.04 / −0.00 / −0.13 %,
every condition within ±0.07 % on transcripts and ±0.32 % on genes, the gDNA pool within 51 fragments; transcripts
reported for zero-truth isoforms 3,707–5,849 → 123–537. Fractional output is bit-identical (the goldens regenerate
unchanged). The design was measured first (`~/Downloads/rigel_runs/prototypes/2026-09-26_whole_counts/`, every EM fragment's posterior dumped
and joined to its true origin): the draw alone — before the repair — landed transcripts within 0.3 % of fractional;
a greedy quota fill stranded 52–95k fragments per library, and the exact optimum, priced in
`ISSUES: whole-counts-optimised-assignment`, trades region fidelity for reads assigned correctly.
The earlier record: Whole-count mode (`assignment_mode="sample"`) gives each fragment to one transcript, drawn from
its posterior after every candidate under 1 % (`assignment_min_posterior`) is zeroed, to keep noise isoforms from
collecting whole counts. On one fragment a real minority and noise look the same — a small share — so the floor takes
a minor isoform, a gene's unspliced RNA or low gDNA off every fragment it shares with a dominant candidate: off
capture the synthetic pool loses 16–29k fragments, transcript error rises 15–23 % and gene error 2.0–2.7× (fractional
mode is untouched; `~/Downloads/rigel_runs/prototypes/2026-09-24_cut_inventory/`). A cutoff relative to the best
candidate (the owner's first idea) spares a fragment with many near-equal candidates but not this. The design (owner:
"I like it in theory"): COUNT FIRST — each component's EM expected count, rounded within its locus so the locus total
is exact (a transcript expecting 0.3 fragments gets 0, one expecting 300 gets 300) — THEN ASSIGN each transcript that
many fragments, the ones it most plausibly produced (a transport problem over the unit posteriors, per locus). Every
fragment still goes to exactly one transcript, with no seed. The rounding half is measured
(`~/Downloads/rigel_runs/prototypes/2026-09-25_floor_study/`): transcripts reported for zero-truth isoforms 4,186 →
196 at `g05 ss.99 OFF` and 4,704 → 403 at `g50 ss.99 ON`, with transcript and gene error unchanged. The assignment
half needs a design and a C++ prototype, then an A/B against today's floor and a relative cutoff. False-positive MASS
(90–96 % of it in a few dozen to a few hundred transcripts that each take a small share of many fragments) is the EM's
own and no assignment rule separates it from a real minority. THE ASSIGNMENT HALF, PROTOTYPED 2026-09-25
(`~/Downloads/rigel_runs/prototypes/2026-09-26_whole_counts/`, `summary.txt`): every EM fragment's posterior row
dumped from the solver (it re-sums to `count_em` to 7e-8) and joined to its true origin by read name, on `g05` and
`g50` of the three in-scope strata. Quotas are each component's fractional count (transcripts and the locus's gDNA)
rounded within the locus. The COUNT-EXACT DRAW — one pass, each fragment drawn from its posterior re-weighted by
(quota left ÷ posterior mass left) per candidate, a fragment whose candidates are all full keeping its best — lands
transcripts within +0.19 to +0.26 % of fractional off capture and −0.03 to 0.00 % on it, genes +2.0 to +4.1 % off
capture and −0.3 % on it (the 916–2,634 fragments per library it could not place), reports 190–400 zero-truth
transcripts against fractional's 4,186–4,704, gives each read a draw from its own posterior (reads correct = the
posterior's own expectation, 65.2–78.4 %), and runs in 0.3–0.5 s per library; visiting each locus's fragments in a
random order changes nothing (the simulator's BAM is grouped by origin). Against it: today's floor +24.5 to +35.5 %
transcripts and 2.3–3.4× genes off capture (+1.6 to +3.2 % / +6 to +45 % on); the relative cutoff +16 to +23 % /
1.9–2.9×; a plain draw, with or without dust removed, +5.6 to +7.6 % / +13 to +16 %; a greedy quota fill strands
52–95k fragments per library (their candidates already full) and costs +40 to +53 %; taking the largest re-weighted
candidate instead of drawing strands 23–28k and costs +7.6 to +13.9 %; the argmax assigns 72–85 % of reads correctly
but multiplies transcript error 5–20×. The optimum under the same counts is priced in
`ISSUES: whole-counts-optimised-assignment`.

### the-scan-fraction-banks-are-not-reproducible
REFUSED 2026-09-25 (owner): bit-identical results are not a goal for now; the tally keeps its float sums and the
scan's thread noise is accepted. The facts it was refused with: the scan's workers take batches from one queue as
they finish and each sums its own float64 fraction banks, so the last bits change from run to run at any scan thread
count above one; at `g50 ss.99 OFF` (two runs, EM on one thread, fractional) the six banks differ at ~1e-15, 64 loci's
gDNA prior at ~3e-16, and the EM's forks carry that to 317 transcripts by up to 0.88 fragments and nascent parents by
up to 42; with the scan on one thread all 214 stage arrays repeat bit for bit
(`~/Downloads/rigel_runs/prototypes/2026-09-26_repro/FINDINGS.md`). The repair priced: exact sums (every deposit a
whole number of 2^-94 units in an unsigned 128-bit integer, rounded once on export) made the whole pipeline repeat at
the default thread budget (213 of 214 arrays; the last is the buffer's row order), moved the answer once by about one
run-to-run spread, and cost ~50 MB per scan worker at human scale plus an exact-summing specification
(`exact_sum.patch` beside the findings). `--threads 1` gives bit-identical output (`MANUAL.md`).

### the-minus-strand-coverage-weight-sat-one-fragment-length-off
FIXED 2026-09-25 (`scoring.cpp`'s `coverage_weight`; gate `tests/test_pipeline_routing.py`: every fragment's weight
against the trapezoid coverage model integrated exactly — 29 placements on one exon, a short transcript and a spliced
one, both strands — beside a strand-mirror test and a multimapper test over 14 offsets each). The EM's warm start
weights each candidate by where the fragment sits on the transcript: the trapezoid `min(x, w, L − x)`,
`w = min(f, L/2)`, over the fragment's bases. On a minus-strand transcript the old rule flipped the genomic start
into transcript space and read the fragment forward from it, but there the genomic start is the fragment's 3′ end,
so the window sat one fragment length off: `[L − s, L − s + f)` in place of `[L − s − f, L − s)`. The rule now
measures the start as the resolver does (a start in an intron back from the next exon, `ISSUES:
a-start-in-an-intron-was-measured-from-the-wrong-exon`), reverses the window on the minus strand, and cuts it at both
transcript ends, so an overhang weighs nothing. Each perturbation fails the gate: the old flip 50 cases, an intron
start snapped to the previous exon 4, the overhang shifted 6, no cut 20; the flip alone, without the other two, 8.
It seeds only the start, so it reaches the answer only through VBEM's start-dependence (`ISSUES:
the-em-answer-depends-on-where-it-starts`). LADDER against the SQUAREM-landed arm (all 16, fractional;
`~/Downloads/rigel_runs/prototypes/2026-09-26_coverage_landed/ab_table.txt`): in-scope transcript Σ|Δ| 390.0k →
391.6k (+0.41 %, stranded OFF), 1,489.8k → 1,487.8k (−0.13 %, stranded ON), 395.5k → 396.1k (+0.16 %, unstranded
OFF); genes +0.32 / −0.02 / +0.27 %; the gDNA pool's Σ|error| within 142 fragments. Rows move by −730 to +1,159
transcript fragments (`g50 ss.99 OFF` +1.24 %), above the run-to-run spread (same-build reseed pairs differ by 0–345,
median 7), so the moves are real: the fix changes only the start, and VBEM's answer depends on its start. (Corrected
2026-09-25: first recorded as inside the floor, read against a reseed arm from an older build — a cross-build
difference, not a floor.) Goldens moved by at most 2.3e-8 relative.

### the-squarem-clamp-decided-which-components-live
FIXED 2026-09-25 (`DESIGN.md` §3.1e; `em_solver.cpp`'s `backtracked_squarem_step`, VBEM and MAP; gate
`tests/test_em_start_independence.py`: three real loci the clamp forked by 7–72 fragments between the coverage and the
uniform start, every case failing on the old solver, and each mode's cases alone failing when its backtrack is
removed). SQUAREM's extrapolation overshot a shrinking component below the floor and CLAMPED it there, and a component
at the floor takes no responsibility and no share of the evidence-proportional prior, so it never came back: which of
the components sharing fragments survived depended on the warm start. `MANUAL.md`'s FAQ said a component with genuine
read support recovers — true only for one holding deterministic fragments of its own. It surfaced through the
minus-strand coverage-weight fix, which changes only the warm start yet moved the ladder ±1–2 %, and at
`em.iterations=10000`, `em.convergence_delta=1e-9` the two trees still differed by −1.8k to +5.8k transcript
fragments. Now the step is halved toward the plain double step until no component that step keeps alive is carried
below the floor. MEASURED on `g00 ss.99 ON`, every setting re-solved from one pre-EM state in one process
(`~/Downloads/rigel_runs/prototypes/2026-09-25_em_start/`): the start moved 41,343 fragments under VBEM and 22,447
under MAP, every moved locus clamped and none a flat valley; with backtracking MAP moves 204 (loci at the iteration
cap only) and VBEM 8,096 (its own, `ISSUES: the-em-answer-depends-on-where-it-starts`); from one start the repaired
MAP beats the clamp's likelihood in 116 of 120 loci; false gDNA there (truth 0) VBEM 4,287 → 177, MAP 1,187 → 158.
LADDER, the landed build against the code it replaces (all 16, fractional; `landed/landed_vs_committed.txt`):
in-scope transcript Σ|Δ| 407.7k → 390.0k (−4.3 %, stranded OFF), 1,496.1k → 1,489.8k (−0.4 %, stranded ON),
409.6k → 395.5k (−3.4 %, unstranded OFF); genes −6.1 / −0.3 / −7.6 %; the gDNA pool's Σ|error| −2.3 / −5.8 / −2.3 %,
and every `g00` row's false gDNA ~4.5k → ~150. It unmasks an under-call at `g05`, which the clamp's false gDNA had
partly offset: stranded OFF −1.3k → −4.1k, unstranded OFF −4.0k → −6.6k. MAP, with or without the repair, is worse
than VBEM on the gDNA pool (`g50 ss.99 OFF` −15.2k against −9.5k), so VBEM stays. Goldens moved by at most 3.3e-8
relative (the plain step is now copied at step 1, not recomputed). COST, MEASURED 2026-09-25 on VCaP at 8 threads
(`~/Downloads/rigel_runs/prototypes/2026-09-26_squarem_timing/`): counted, so independent of machine load, +8.8 % SQUAREM
iterations and +13.7 % E-step work, loci at the 333-cycle cap 103 → 123 (the largest locus, 12,890 components, now
among them: 299 → 333); the halvings themselves are ~1e8 component checks, negligible. Three interleaved wall-clock
pairs, taken beside another workload, put the locus-EM stage at +12 % on average (−3 to +27 % per pair, ~+3 s of
~23 s) and the whole run unchanged within ±7 %.

### a-start-in-an-intron-was-measured-from-the-wrong-exon
FIXED 2026-09-25 in the resolver (`resolve_context.h`'s `tx_frag_length`, the one length definition of `DESIGN.md`
§3.4; gate `tests/test_transcript_space_fl.py`'s `TestLengthRuleByEnumeration`, a base-by-base count over every
two-mate fragment with its ends at or beside an exon boundary). The EM's per-candidate fragment length projected both
ends into transcript space and measured an end lying in an intron forward from the PREVIOUS exon's end — right for the
fragment's end, wrong for its start, which landed one intron length downstream: the length read |true − intron| (a
200-bp fragment starting 1 bp before a 2.8-kb intron's end read 2,610), and a length of exactly 0 was stored as
missing (−1), which the scorer reads as log-likelihood 0, the best fit there is. An end whose last base was an
intron's last base read as ending at the next exon's start and dropped the whole intron. Two hand-written tests
asserted the defective values (398 and 2,006; now 3,397 and 5,005), and neither boundary error is caught by any
hand-written test: perturbing either one fails only the enumeration. Before the fix
(`~/Downloads/rigel_runs/prototypes/2026-09-24_cut_inventory/`), 1.35–10.7 % of multi-block candidate entries per
condition had a wrong length, median error 1–1.2 kb, and 1.9–2.6k candidates read −1. MEASURED on all 16 conditions
(fractional, the worktree build against the shipped `base` arm;
`~/Downloads/rigel_runs/prototypes/2026-09-25_start_endpoint/ab_table.txt`): transcript Σ|Δ| summed per in-scope
stratum −560 (stranded OFF), −4,277 (stranded ON), −887 (unstranded OFF); the largest row `g50 ss.99 ON` 301.5k →
298.0k; two rows worse beyond the run-to-run spread (`TRAPS: the-deliverable-is-not-reproducible-by-default`),
`g05 ss.99 OFF` +496 (+0.41 %) and `g00 ss.50 OFF` +168; genes within ±1.2 %; the gDNA pool within ±1k everywhere.
The inventory's estimate that the defect pushed ~1.1k gDNA fragments per `g98` capture-ON condition toward RNA did
not survive: the pool moved by 22 there.

### calibration-and-the-em-disagree-on-what-can-be-gdna
FIXED 2026-09-24 in the EM (`DESIGN.md` §3.1d; `scoring.cpp`'s `gdna_can_explain` and `gdna_competes`; gates in
`tests/test_pipeline_routing.py`); the leftovers moved to `ISSUES: splicing-artifacts` and
`ISSUES: multimapper-intergenic-alignments`. Calibration and the EM decided by different rules whether a fragment
can be gDNA: calibration offered every fragment its unbroken reading and drew among the survivors, while the EM
denied gDNA to every label but "unspliced" and removed it from a multimapper's every hit when one hit was spliced.
Now an unspliced, implicit or artifact alignment may be gDNA, a multimapper may if any alignment may, and gDNA is
pruned exactly as an RNA candidate is. The implicit splice needs no new likelihood: gDNA's term is the unspliced term
at the footprint and the fragment length decides (derived three ways; `DERIVATION.md` in
`~/Downloads/rigel_runs/prototypes/2026-09-24_spliced_rule/`). MEASURED on a one-gene toy on the shipped solver:
denying implicit fragments gDNA put the error on the synthetic span and the annotated transcripts under every
prior — +6 % at a 60-bp intron and gDNA 50 %, +39–62 % at 90 %, +13 % at 150 bp and 90 %. On the ladder, all 16
conditions (fractional, VBEM, on the count form; `~/Downloads/rigel_runs/prototypes/2026-09-24_implicit_gate/ab/`):
gDNA error at stranded ON `g50` −28.0k → −14.1k, `g98` −112.8k → −100.5k, `g05` +1.3k → +2.6k, off capture within
±1.0k; transcript error `g98 ss.99 ON` 54.2k → 50.3k, elsewhere within ±1.7k; gene error lower in 12 of 16 rows.
The gDNA gained at `g50 ON` came from the synthetic spans (1.00× → 0.91×), not the annotated excess (+28.3k →
+27.6k). A first form that limited implicit splices to the maximum fragment length admitted 18k–55k units per
condition whose gDNA reading is about −708 nats, each counting in the prior's strength, and raised `g00 ss.99
OFF`'s false gDNA to +7.0k against +4.4k under the pruning.

### the-pseudocount-prior-is-biased-toward-gdna
FIXED 2026-09-24 by the count form (`DESIGN.md` §3.1c, `EQUATIONS.md` §9b.3, `pipeline.em_pseudocounts`; gate
`tests/test_em_pseudocounts.py`). The EM's two pseudocounts were calibration's gDNA and UNSPLICED RNA counts, but
θ is each component's share of every fragment in the locus — the deterministic spliced fragments enter every
M-step — so the unspliced split stated as the whole pool's leaned toward gDNA in proportion to the spliced
share: its centre read 63.9 % against 50.0 % at `g50 ss.99 ON`, and it over-called an exact toy by 70–944 of 900
as the spliced share ran 11–60 %. Now `P_g = U·min(G_c/N_c, 1)`, `P_R = U − P_g`, the share taken over the
`N_c` fragments calibration counts from (on the ladder `N_c = N`). MEASURED on all 16 conditions
(fractional, VBEM with a MAP pair, the odds changed and the strength held): gDNA error at `g05` / `g50`
unstranded OFF +59.2k → −3.1k and +225.7k → −9.7k, stranded OFF +24.6k → −1.2k and +128.5k → −8.7k, stranded ON
+76.9k → +1.3k and +172.9k → −29.4k (the spans 0.28× → 1.01×); gene error at `g50` 18.3k → 14.4k, 13.6k →
11.4k, 102.6k → 70.2k; the transcript table within 2 %; `g00` unmoved. `g98` worse in every stratum (−18.7k →
−62.8k, −3.1k → −37.8k, −54.2k → −118.3k; genes ON 29.0k → 36.2k): the lean was offsetting
`ISSUES: rna-prior-floor-at-pure-gdna-loci` and the capture likelihood's lean toward RNA, the owner's next
targets. The deferred unstranded × capture-ON stratum moves the other way (`g50` −125.6k → −589.5k). MAP moves
the same way as VBEM, 2–8k further toward RNA. The per-locus A/B, census and tables:
`~/Downloads/rigel_runs/prototypes/2026-09-23_pseudocount/` (`DERIVATION.md`, `ab16/table16.py`).

### a-count-once-density-prior-for-the-strength
REFUSED 2026-09-24 as a buildable rule for the pseudocount's strength; the principle stands. A prior restating
calibration's full answer counts the locus's fragments twice, because the E-step reads them again; the rule that
counts them once puts a Gamma prior on the locus's gDNA density at calibration's estimate from OUTSIDE the locus,
and its exact M-step is `out_g = (raw_g + a − 1) / (1 + a/μ)` with RNA untouched (matched on the shipped MAP solver
to 0.01 fragment). Killing numbers: under capture, where the prior carries the answer, no outside estimate exists —
each object's capture efficiency is read from its own count, so the "outside" centre equals calibration's own gDNA
count to 0.1–1.3 % on the test chromosome, and the locus's flanks are off target, 75–440× below its gDNA. Off
capture it carries three choices that move it at first order: how over-dispersion correlates within a locus (the
test chromosome's `g50 ss.99 OFF` per-locus error 2,375 or 3,349 against the shipped strength's 2,461), the
moment estimate's `α = ∞` boundary (three of four ladder conditions), and MAP against mean against VB (30 % apart
at `a` = 3.3 on the unstranded ridge). At the ladder's own reading it is worse than the shipped strength on all
three capture-OFF test conditions (3,349 / 4,754 / 1,542 against 2,461 / 4,262 / 1,189), and on the shipped
grouped VBEM it erases the split within RNA (a 6-fragment isoform 4.75 → 7.37). Scripts: `toy_gamma_form.py`,
`refute_kappa_*.py` in `~/Downloads/rigel_runs/prototypes/2026-09-24_spliced_rule/`.

### a-junction-price-read-from-a-fitted-capture-map
REFUSED 2026-09-23 at the no-EM gate. The general form of conservation of bases: capture additive over bases, every
region a flat per-base capture read from gDNA's own objects (non-negative least squares over every region's
contained and every boundary's crossing efficiency, each row weighted by its support, per reference), and each
junction priced by routing the transcript's own placements through that map — the shipped sum is its limit where
the pieces beside a junction exceed every fragment. Within-gene sd at `g05` / `g50 ss.99 ON`, on the transcripts'
own routed shares where the sum reads 0.050 / 0.054 on the simulator's noise-free counts and 0.074 / 0.058 on
calibration's: the map 0.060 / 0.064 and 0.089 / 0.073, the unprobed class inflated further (1.52 / 1.50
noise-free against 1.00 / 1.01). With the TRUE per-region probe coverage under additive capture it fixes the
geometry (junction rms 0.58 → 0.19, 1.013 on average); against the simulator's capture it prices junctions 2.1×
high, because probes designed on different isoforms overlap, and overlapping probes add under additive capture but
not where a fragment binds its best single part (`ISSUES: the-junction-price-is-noisy-within-a-gene`).

### the-junction-sum-on-unclipped-posteriors
REFUSED 2026-09-23 (`ruler_vs_truth.py --module`, the prototype's `junction_rulers.py`). The sum by conservation of
bases is linear in density, so it was taken on each of its four terms' unclipped posterior `E[ρ / ρ_ref]` in place
of the published `E[min(ρ / ρ_ref, 1)]`, every other object unchanged. The class medians move toward 0 (partial
mRNA +0.035 → +0.024 at `g05`, +0.019 → −0.020 at `g50`) and the within-gene sd widens, 0.074 → 0.080 / 0.060 →
0.066 (on the simulator's noise-free counts 0.050 → 0.052 / 0.054 → 0.060). The junction keeps the published
clipped efficiencies.

### the-gdna-component-length-rule-differs-from-the-transcripts
RESOLVED 2026-09-23 by the one shared rule (`DESIGN.md` §7.2, `EQUATIONS.md` §11): every EM component — the locus
gDNA component, every synthetic span, every annotated transcript — has as its capture-contracted length its
conserved share of every object its fragments deposit on times that object's capture efficiency. The defect:
transcripts and spans took the per-base rule (their own bases at their pieces' efficiencies with the fragment-end
taper) while the gDNA component took region and boundary objects with the crossing supports converted by the
pooled `q`, and the two rules disagreed by 14–16 % on capture. Killing numbers (`ruler_vs_truth.py --scale`, L / Y
×10⁻³ against the simulator's exact yields, gDNA component / spans / isoforms, and the largest gap): `g50 ss.99 ON`
1.108 / 0.933 / 0.957 — 18.7 % → 1.126 / 1.136 / 1.111 — 2.3 %; `g05 ss.99 ON` 1.107 / 0.932 / 0.944 — 18.7 % →
1.087 / 1.105 / 1.075 — 2.8 %, in the prototype (calibration's counts on the count frame as clipped plug-ins). The
shipped tree, with the published posterior efficiencies (width step 1), reads 1.131 / 1.141 / 1.112 — 2.6 % and
1.108 / 1.118 / 1.076 — 3.9 % at the two rows. Within each locus, at the two rows, in the prototype:
gDNA/isoforms 1.149 / 1.153 → 0.983 / 0.986, gDNA/spans 1.202 / 1.208 → 0.984 / 0.981, spans/isoforms 0.978 /
0.977 → 1.014 / 1.014. The true gDNA on every object reads the same ~2 %, the rule's own residual
(`ISSUES: ruler-witness-geometry-on-transcript-panels`). The port reproduces the prototype A/B'd to 9e-15 on
transcripts except the 113 at the fl floor (shorter than every fragment), and on the gDNA component except the loci
sharing a region with another locus (34 of 1,216 at `g50 ON`, 3.7 % of the gDNA prior; `P_g`-weighted totals
0.999998), where the length now takes the count's own overlap weights. What it cost through the EM is two open
entries: the pseudocount's bias no longer cancelled (`ISSUES: the-pseudocount-prior-is-biased-toward-gdna`) and
the isoforms split by the junction price's noise (`ISSUES: the-junction-price-is-noisy-within-a-gene`). Refused on
the way: `ISSUES: synthetic-spans-on-the-object-rule`, `ISSUES: the-per-base-rule-for-every-component`,
`ISSUES: the-pooled-q-in-the-gdna-length`.

### pooling-junctions
REFUSED 2026-09-23 by the owner, before it was built: a posterior pulling each junction's price toward its
neighbours' price times a ratio learned across all junctions, proposed to cut the sum's noise
(`ISSUES: the-junction-price-is-noisy-within-a-gene`). Some junctions are probed across the junction itself and
some are not, so a ratio pooled over junctions is theoretically invalid: a gain on the ladder would be chance, and
another probe design could read much worse (`DESIGN.md` §7.2). Do not re-propose a junction price that pools
across junctions.

### a-junction-price-clipped-at-one
REFUSED 2026-09-23. The sum by conservation of bases is floored at 0 and not clipped at 1 — each boundary is an
unclipped half-exon level — and at `g50 ss.99 ON` 18,308 of 45,609 junctions price above 1 (max 1.999) and 4,440
below 0 (`g05 ss.99 ON`: 17,167 above, 7,409 below). A clip at 1 shortened 25 % of transcripts by a median 13 %
(Σ 0.965 of the validated rule): a second clip on pieces already clipped inside their own posteriors.

### the-spliced-read-junction-price
REFUSED 2026-09-23: the owner's proposal to impute a junction's capture from the transcripts' own reads left and
right of it (derived with the correction ÷ `f_spliced`, not ×) — exact in principle, and starved. The pieces
beside junctions have a median length of 107 bp; 58 % of single-isoform junctions hold no unspliced read within a
fragment; and "exonic on the strand" neighbourhoods reach retained-intron-only regions and inflate the price
(per-transcript L up to ×60).

### the-crossing-apportionment
DELETED 2026-09-23 with the one shared rule (`capture_efficiency.capture_efficiencies`, no `shares` argument):
a region's efficiency had read its own contained count plus the crossings at every boundary within a fragment's
reach, apportioned to the pieces its fragments cover (`ISSUES: ruler-multimapper-floor-caps-the-correction`). Each
object now reads only its own count — a region its gDNA contained count on its contained support, a boundary its
gDNA crossing count on its crossing support — and the length prices every crossing at the boundary that holds it,
so a region never borrows the crossings around it. The apportionment had read a tiny piece's efficiency off its
edge crossings, which the rule now prices at those boundaries directly; a piece too short to contain a fragment
has no share, and its efficiency (the population's) multiplies nothing. A deletion by structure, not by a number.

### the-pooled-q-in-the-gdna-length
REPLACED 2026-09-22/23 by gDNA's own conserved share at each boundary (`CalibrationResult.gdna_boundary_conserved_len`,
`EQUATIONS.md` §11). The gDNA component's length had converted a boundary's crossing support by
`q = boundary_mass_per_crossing`, a mass per crossing pooled over gDNA and RNA and so RNA's where RNA dominates a
boundary: against the oracle's gDNA-only payload it overstated gDNA's mass by 1.4 % overall at `g50 ss.99 ON`,
3.4 % at exon|exon boundaries, 8.4 % where RNA is 50–90 % of the crossings, 5 % at `g05 ss.99 ON`, and it read
transcripts 9 % high off capture. Storing the unspliced boundary mass per strand column would give gDNA's own mass
exactly on stranded data (a 2×2 solve per boundary; unstranded data cannot split it) — not built: the length
needs neither. The two pseudocounts still convert a boundary's mass by `q`.

### the-junction-price-at-its-neighbours-mean
REFUSED 2026-09-22 as the scale. Against the simulator's own junction capture (60 probed transcripts, 395
junctions) the junction's two boundaries' mean prices it at 1/1.68, the adjacent pieces 1/1.73, the nearest
measured object along the transcript 1/1.74 and the transcript's own measured mean 1/1.24 (a transcript-level
anchor, barred by the owner's locality ruling, `DESIGN.md` §7.2), where the sum by conservation of bases reads
0.959. On the class scale, with the junction at its two boundaries or its adjacent pieces the gDNA component and
the spans agree to 1 % (0.991 / 0.993 with the true gDNA, 1.137 / 1.147 with calibration's) while the transcripts
sit 10–20 % below (0.910 / 0.970 with calibration's efficiencies). KEPT AS DATA: priced at its adjacent pieces a
junction is as precise within a gene as the per-base rule (within-gene sd 0.057 / 0.041 at `g05` / `g50 ss.99
ON`) and brings the class gap back (`ISSUES: the-junction-price-is-noisy-within-a-gene`).

### the-per-base-rule-for-every-component
REFUSED 2026-09-22 by the owner's ruling (`DESIGN.md` §7.2): a length over each component's own bases at its
pieces' capture efficiencies with the fragment-end taper omits the boundaries that own the crossing mass — any
fragment that crosses a boundary is deposited into it, and a piece shorter than a fragment is seen only through its
boundaries. Its no-EM reading was one scale for gDNA and isoforms (0.956 / 0.957 ×10⁻³ at `g50 ss.99 ON`) with
the spans 2.5 % low (0.933), and only because its clipped per-piece efficiency under-priced gDNA's short probed
pieces by about as much as it under-priced the transcripts' junctions (with the locus gDNA's taper on the outside
of its footprint, 3.0 / 2.0 / 0.7 % from one scale at `g50 ss.99 ON` / `g05 ss.99 ON` / `g50 ON` at `on_fraction`
0.10, 2026-09-21). Within a probed class its spans read 0.950 / 0.949 / 0.964 of the isoforms (`partial ≤ ½` /
`partial > ½` / `probed ≥ 0.9`) and the locus footprints 0.978 / 0.917 / 0.521, and no region-constant per-piece
efficiency closed it (the sampler's own covering truth 1.271 / 1.155 / 1.124, its contained truth 1.150 / 1.007 /
1.059, the probed fraction 1.122 / 1.022 / 0.979; the clipped estimator was the closest; 2026-09-21).
Through the EM with the corrected pseudocount it read 5.49 → 5.14 % transcript error at `g50 ss.99 ON` with the
synthetic pool 1.22× → 1.65×. Re-measured 2026-09-23 beside the one shared rule (fractional, stranded × capture
ON, with the corrected pseudocount on both): transcripts `g05` / `g50` / `g98` 3.47 / 5.15 / 30.43 % against 5.55 /
6.20 / 28.00 % — refused by structure (owner's ruling), not by this number (`DESIGN.md` §0b). Its form for the gDNA
component alone — the locus's bases at the regions' efficiencies — was refused 2026-09-20: it drops the boundary
objects the count keeps; test chromosome `g50 ss.99 ON` gene-level Σ|Δ| 25,633 → 38,174 against 23,967 with the
object form, the gDNA pool +5.5 % against −0.1 %. The taper and the per-base frame are deleted.

### background-abundance-pair-unruled
DELETED 2026-09-22 (the cleanup; owner: "I don't know what measured_total is for"). `CalibrationConfig.background_abundance`
chose which (counts, exposure) pair the pooled gDNA background took — the shipped `"contained"` pair or the START/END banks
over the region's own length. An alternative nobody ruled on is a tunable, and it went with its refusal path, its
`counts_exposure` parameter and its three tests; the contained pair is the only one. `rename_identity.py --check`
bit-identical.

### the-fixed-gdna-level
REFUSED 2026-09-21 as the shipped kernel (33 converged arms on the ladder, one EM thread, fractional). Holding the
gDNA component's count at calibration's locus count every iteration and dropping the RNA pseudocount repairs the
capture-OFF pools exactly as the corrected pseudocount does and halves the `g00` false positives (+3,913 → +1,902),
but at `g98 ss.99 ON` it collapses on every length rule — 46.1 % (object) and 45.4 % (per-base) transcript error
against 28.6 / 30.4 % for the corrected pseudocount, the synthetic pool at 128× / 61× — and 36.6 % even with the
simulator's own lengths for every hypothesis (the corrected pseudocount 23.5 %): at 98 % gDNA a 1 % error in the
level is half the RNA, and the pseudocount's damping is what absorbs it. On real data it would also need a rule
for the 194–288 unwitnessed loci per library and a `(1 − M)` multimapper term (`ISSUES:
unwitnessed-loci-and-multimappers-at-the-em`). Not rebuilt without both.

### the-unclipped-capture-efficiency
REFUSED 2026-09-21. Publishing `E[ρ / ρ_ref]` unclipped in place of `E[min(ρ / ρ_ref, 1)]` — so that pieces above
the reference keep their ratios — moves the per-base family from 3.0 % to 13.2 % from one scale at `g50 ss.99 ON`
(spans/isoforms 0.897, gDNA/isoforms 0.883): it inflates the isoforms' partial exons (1.35) more than the spans'
(1.21), because against the sampler's covering truth the unclipped estimator over-separates the exon classes
(fully probed 1.45, `½ < p < 0.9` 1.31, `p ≤ ½` 1.28) where the clip under-separates them (1.01 / 1.12 / 1.22). The
clip stays.
NOTE 2026-09-23: measured on the per-base rule, since refused (`ISSUES: the-per-base-rule-for-every-component`);
the one shared rule keeps the clip per object, and the unclipped form is re-priced on it before this is quoted
against it.

RE-EXAMINED 2026-09-24 and not re-opened: a reading that captured entities' contracted lengths follow only part
of the true capture range ("elasticity" 0.64 multi-exon at g50) came from weighting transcripts by the EM's own
gDNA exchange, which depends on the length ratio itself; with weights from outside the EM it reads 0.91–1.06
(`analysis/refute-mechanism/elasticity_exog.py` in `~/Downloads/rigel_runs/prototypes/2026-09-24_g98_dissection/`).
### nascent-siphons-gdna-under-capture
REPAIRED 2026-09-20 (`c52c9b93`). THE MECHANISM (per fragment against the read names' truth): the locus gDNA
component's opportunity counted a crossing start at EVERY boundary its fragment crossed, while its pseudocount
counted the fragment once — `assemble_priors` converted a boundary's incidence count by the accumulator's `q` and
left its incidence SUPPORT unconverted (`EQUATIONS.md` §11). Under capture the crossing support is 64 % of a
locus's opportunity (15 % off), so gDNA was priced at `ρ q̄`, handed its sense fragments at the probed exons to the
RNA hypotheses, and the unpinned synthetic shadows took them (95 % of their gDNA was exon-contained sense
fragments; the leak 17.1 % of a locus's gDNA at `q̄` 0.45–0.55, 0.5 % where crossings cross one boundary). One
E-step from the TRUE counts drifted +24 % per step on capture and was a fixed point off it. THE REPAIR converts the
crossing support by the same `q` (gates in `tests/calibration/test_priors.py`, watched to fire): `g50 ss.99 ON`
nascent +541,216 → +32,908, gDNA −626,550 → −56,319, transcript Σ|Δ| 5.21 → 5.50 % (the isoform ruler's own
over-statement unmasked, `ISSUES: ruler-witness-geometry-on-transcript-panels`); the deferred stratum's `g50`
siphon +1,567,550 → +569,674; `calibration_vs_oracle.py` and `policy_benchmark.py` unmoved. RULED OUT, with the
killing number: the ruler's witness geometry (shipped annotated/synthetic contraction ratio 0.0773 against the
simulator's 0.0738); the shadow-exclusive intronic pool (2.76 % of contained gDNA under capture); the per-locus
prior (structurally a gDNA:RNA lever, `--arm oracle` 541,216 → 541,762); nascent RNA itself (+522,205 at
`on_fraction` 0.10). A length knob on the shadows (doubling `L_n` removes 97 % and destroys the true nascent
signal) is a falsification probe, never a repair; the per-gene gDNA opportunity is refused
(`ISSUES: per-gene-gdna-opportunity`). The residual rides with `ISSUES: per-transcript-prior-lane`.
NOTE 2026-09-23: the length's `q` conversion is replaced by gDNA's own conserved share, which holds each start once
as this repair required (`ISSUES: the-pooled-q-in-the-gdna-length`).

### exclusive-evidence-support-rule
REFUSED 2026-09-20: an isoform supported iff its EXCLUSIVE spliced evidence carries reads (the per-transcript
prior's step ①), fed as a static weight through `rna_prior_weight` (the control reproducing base to 0.00 points),
recovers ~0 % of the allocation ceiling — transcript Σ|Δ| base → rule → true support → ceiling: `g50 OFF` 2.49 →
2.49 → 2.23 → 0.70, `g50 ON` 5.49 → 5.47 → 4.54 → 2.07, `g05 OFF` 1.59 → 1.59 → 1.55 → 0.43, `g05 ON` 3.49 → 3.47 →
2.83 → 1.08 (`ss 0.99`, fractional, `71037f6c`). The rule separates the testable set perfectly, but the EM's
likelihood already switches those isoforms off (122 false fragments of 9,950 at `g50 OFF`). Do not re-propose an
exclusive-evidence rule on junctions, regions or boundaries: the likelihood already uses exclusive evidence.

### a-per-object-gdna-rate-in-the-e-step
REFUSED 2026-09-20 (the calibrated-likelihood campaign's step 2, re-verified by the second adversarial review on a
corrected harness). Replacing the locus-pooled gDNA rate `count_g / length_g` in the E-step by calibration's
per-object gDNA density at the fragment's object gave no gain off capture and was worse on capture. One E-step from
the truth, the synthetic pool as a multiple of its truth: off capture every form reads 1.00–1.02; at `g50 ss.99 ON`
the locus-pooled fixed rate 1.14× (1.16× with the certified locus count), the per-object gDNA density 4.21×
(`local_gdna`), 1.18× even with the CERTIFIED per-object density (`local_true`), and 2.72× with a per-object RNA
rate beside it (`local_both`, 87 % of it a harness hard zero; the strike stands without it); at `g05 ss.99 ON`
1.04 / 1.25 / 0.95 / 1.18 in the same order.
A per-object RATE re-enters the gDNA-versus-span contest in calibration's density units against the ruler's; the
per-fragment gDNA PROBABILITY of `ISSUES: the-pseudocount-prior-is-biased-toward-gdna` is dimensionless and removes
that contest, and has not been measured.

### synthetic-spans-on-the-object-rule
REFUSED 2026-09-20. Pricing the synthetic spans by the gDNA component's object rule while the isoforms keep the
ruler puts the spans 18.7 % above the isoforms against the simulator's yield (`g50 ss.99 ON`) and the converged EM
collapses their pool to 0.55× its truth; a mixed state is the knife edge, never a repair.
NOTE 2026-09-23: the spans now share one rule over regions and boundaries with the isoforms and the gDNA
component (`ISSUES: the-gdna-component-length-rule-differs-from-the-transcripts`); what stays refused is a mixed state.

### per-gene-gdna-opportunity
REFUSED 2026-09-20, before it was built, on the measurement its premise fails. The candidate (2026-09-19)
rested on a toy in which every fragment of the locus sits inside the shadow's footprint, where the
shadow's growth factor is `L_g/L_n`; in a locus whose gDNA is uniform the footprint holds its share and
the factor is the density ratio, so pooling a locus's genes under one opportunity destabilises nothing by
itself. Measured on `g50 ss.99`: off capture `L_g/L_fam` reaches 5 with families UNDER-calling 7 %; on
capture families over-claim 3.5× where their own footprint's gDNA density per unit of the longest twin's
opportunity is exactly the locus's (D 0.8–1.2), and single-twin families leak 4.3× against 3.7× for
families of ten. A per-gene split of the shipped opportunity would inherit the twice-counted crossing per
gene (`ISSUES: nascent-siphons-gdna-under-capture`). Do not rebuild it.

### em-overturns-the-calibrated-gdna-split
CLOSED 2026-09-19, in two halves and by two different things. THE CAPTURE-OFF HALF CLOSED BY LANDING
`nascent-gets-no-rna-prior`: the EM's gDNA over-call there was the RNA prior's eligibility rule and not
the EM's arbitration — `g50 ss.50 OFF` reads a library gDNA fraction of 0.5127 against 0.50 where it
read 0.574, `g05 ss.50 OFF` 0.0538 where it read 0.104, and the nascent shortfall that mirrored it
closed with it (`g00 ss.50 OFF` nascent 1,901,090 → 2,018,540 against 2,024,341 true). THE CAPTURE-ON
HALF IS SUPERSEDED by `ISSUES: nascent-siphons-gdna-under-capture`, which names it correctly: the
missing gDNA does not go to the isoforms of probed genes, it goes to the SYNTHETIC NASCENT entities,
one fragment for one. The entry's own leading hypothesis — one gDNA rate spread uniformly along a locus
while capture concentrates gDNA at probed exons — survives intact and is carried there as a candidate;
what did not survive is the claim that the mass lands on annotated isoforms, which the pool rows refute
(annotated Δ is −2,752 / −6,560 / +20,767 at `g05` / `g50` / `g98` `ss.99 ON` against a nascent Δ of
+69,268 / +541,216 / +590,406). Do not re-open this name; the open thread is the siphon's.

### nascent-gets-no-rna-prior
CLOSED by landing 2026-09-19 (owner: "it's a hack; restore nascent RNA fairness"). The EM's RNA pseudocount
now goes to EVERY RNA component in proportion to the evidence it already carries, with none singled out for
zero; the eligibility test, the per-component `component_is_synthetic` flag, the `t_is_synthetic` lane and the
`index` parameter it was the only reader of are deleted (`EQUATIONS.md` §9b–§9b.2, `DESIGN.md` §0b).
THE WEIGHT — the one open choice — is `w_i = raw[i]`, the shipped weights, admitting every component at them
rather than an equal share (owner, 2026-09-19). Three reasons, each with its number: it keeps §9b.1's ABSORBING
STATE at zero evidence for free, since the state is a property of the weights and not of the eligibility test;
an equal share would be a strong informative prior rather than a neutral one, because `rna_prior_count` is a
conserved FRAGMENT COUNT (tens to thousands on an expressed locus), so `P/n` clears the ~0.16–0.47-alpha-unit
activation threshold by one to two orders of magnitude at every component; and it leaves
`component_rna_prior_weight` free for `ISSUES: per-transcript-prior-lane`, whose whole point is to fill that
lane with a MEASURED weight. THE PRICE, stated: the prior no longer helps a shadow entity decay — the rate
falls from `kappa/(1 + P/R)` per M-step to `kappa = w_N/w_T`, still strictly below 1 for free.
MEASURED. The defect was far larger than the entry recorded: on a tied two-component locus with one component
synthetic, a prior of 500 drove it from 100.0 fragments to 2.79e-298 and handed all 200 to the other. The
strict xfail this entry owned reads 70 → **0** fragments of leak onto the unexpressed antisense `t2`
(`test_nrna_multiexon_t2_low_ss`; the "72" in this entry and in `NEXT_SESSION.md` was stale). The gDNA:RNA
split does not move: the per-M-step identity `Σ out over RNA = rna_count + rna_prior` is gated unconditionally
in `tests/native/test_grouped_prior_update.py` and watched to fail under a perturbation that leaks the prior
across the gDNA boundary. A locus with no synthetic component is BIT-IDENTICAL — the operation order is
preserved for exactly that reason. Gates in `tests/test_estimator.py` (the prior moves no component's share of
the RNA pool, exact under MAP and within the derived digamma residual under VBEM; a shadow span still loses to
the transcript it shadows; a supported component keeps its mass; a zero-evidence component is not revived),
each watched to fail under its perturbation. ⚠ It also closed `nested-antisense-leak-under-the-sane-ruler`.

### nested-antisense-leak-under-the-sane-ruler
CLOSED by landing 2026-09-19, as a consequence of `nascent-gets-no-rna-prior` rather than by its own thread.
Both strict xfail rungs pass: the leak onto the unexpressed antisense `t2` falls 124 → **14** at SS 0.65 and
24 → **2** at SS 0.9, against the tight bounds of 20 and 5 the tests were written with. The entry blamed "the
EM's assignment at a nested transcript nothing witnesses" and ranked its lever as the per-transcript prior
lane; the diagnosis was one layer too deep. The leaked fragments were the HOST's nascent entity's, denied
their share of the locus's RNA pseudocount and landing on the one annotated transcript that could also explain
them — so the allocation rule was the mechanism, and the ruler was never implicated. ⚠ The residue at SS 0.65
is real: 14 of 2,000 fragments still reach a transcript whose truth is 0, at the strand specificity where the
separating channel is weakest. It is inside the bound and no longer an xfail; if it grows, it is the
unwitnessed-nested-transcript effect and `ISSUES: the-atom-at-an-unwitnessed-both-strand-slot` is its twin.
`tests/scenarios/test_antisense_intronic.py`, `quant_accuracy.py`.

### oracle-cache-key-hashes-a-thread-count
CLOSED by landing 2026-09-19 (owner): the scan settings' part of a cache's key is DERIVED at read time from the
settings the manifest records, and leaves out the two thread counts (`scan_cache._WORK_ONLY_FIELDS`), which divide
the work and not the tally (`test_scan_order_independence.py`, one worker to eight). No digest string is stored
any more, so the key's definition can move without stranding a cache: every oracle cache written 2026-08-22 at
`bgzf_threads` 4 loads again under the default config, with no re-scan. The five other resource settings
(`log_every`, `fragments_per_chunk`, `read_name_batch_size`, `buffer_size_bytes`, `spill_dir`) stay in the key:
their invariance is not gated. Gates in `test_scan_cache.py` (a thread count is not the key; a tally setting is; the
key is derived and no stored string decides), each watched to fail under its perturbation.

### oracle-ruler-arm-cannot-reach-the-ruler
CLOSED by landing 2026-09-19 (owner): `quant_accuracy.py --arm oracle_ruler` now hands the EM `fl × factor`, the
factor being the simulator's own capture-aware length (`ruler_vs_truth.load_truth`), anchored on the fully probed
class, cached per capture label beside the oracle caches; `oracle_ruler_noop` builds the same and hands back the
shipped lengths. The arm it replaced swapped calibration's count arrays, which the ruler stopped reading at
`c44fc306`, and refused every condition. Gates in `test_quant_accuracy.py` (the lengths reach the EM's published
`em_effective_length` and move the split, the noop reproduces `base`, the guard, the anchor), each watched to fail
under its perturbation — the anchor's first fixture could not tell the two pools apart and was repaired.

### end-to-end-error-unattributed
CLOSED by measurement 2026-09-19: the residual is the EM's and the assignment's. `quant_accuracy.py` on the
16-condition ladder at `cffca248`, seed 20260807, arms `base`, `base_reseed`, `noop`, `oracle`, `oracle_gdna`,
`oracle_rna`, `oracle_efflen` and `oracle_alloc_seed` (`~/Downloads/rigel_runs/arms/2026-09-19_e2e_baseline/`,
`decomposition.txt`). Per stratum, the `g05`–`g98` rows summed, in the order unstranded OFF / stranded OFF /
stranded ON and then the deferred one: transcript-level Σ|Δ| 520,117 / 453,135 / 1,833,492 / 4,487,120 fragments,
4.4 / 3.9 / 12.8 / 31.3 % of the true annotated RNA, against floors |base − base_reseed| of 3,478 / 5,802 / 3,532 /
5,079; gene level 246,441 / 210,165 / 601,452 / 2,640,220 (2.1 / 1.8 / 4.2 / 18.4 %), floors 558 / 282 / 2,228 /
3,605. The gDNA-free `g00` rows read 3.4 / 3.6 / 7.5 / 9.1 % at transcript level and `g05` 3.7 / 3.2 / 10.5 /
16.6 %, so most in-scope error exists with no gDNA at all. A PERFECT PRIOR (`oracle`) recovers 26,401 / 12,837 /
75,687 fragments at transcript level (5.1 / 2.8 / 4.1 %) and 24,913 / 11,497 / 74,687 at gene level (10.1 / 5.5 /
12.4 %); `oracle_rna` alone carries 21,697 / 14,074 / 73,855 of it, `oracle_gdna` 1,805 / 2,433 / 3,523, and
`oracle_efflen` −50 / +22 / −1,104 — ⛔ WITHDRAWN 2026-09-21: that arm overrides count arrays only and the
length reads none, so it cannot fire and those deltas are noise (`ISSUES: capture-blind-gdna-divisor`). By rung it is `g98`'s: at `g05` the move is below the floor
in all three strata (−1,694 against 2,431; +2,013 against 3,722; −1,373 against 1,942), at `g50` it is 2.0 / 2.8 /
1.2 %, at `g98` 37 / 25 / 35 %. Deferred: 40.9 % transcript, 70.4 % gene. So 95–97 % of the in-scope transcript
error and 88–95 % of the gene-level error survive a perfect prior, and what survives splits in two — the
per-transcript allocation (`ISSUES: per-transcript-prior-lane`: truth as the weights removes 49 / 46 / 61 %) and
the capture-OFF gDNA over-call (`ISSUES: em-overturns-the-calibrated-gdna-split`), which neither arm
moves. THE OLD CLAIM — in scope a perfect prior no longer improved the transcript number (2026-09-10:
`oracle`/`base` 1.019 / 1.026 / 0.978 transcript, 1.087 / 1.169 / 0.905 gene) — no longer holds as worded: today
0.949 / 0.972 / 0.959 and 0.899 / 0.945 / 0.876, each above its floor except at `g05`; which of the commits
between the two readings dissolved it is not attributed. Its substance holds and is now measured: the prior is
worth little in scope and the EM owns the rest. Two instruments were broken on the way in: `oracle_ruler` could
not fire (`ISSUES: oracle-ruler-arm-cannot-reach-the-ruler`), so the effective-length shrinkage sits inside no
ceiling here, and the oracle arms ran through a key-only wrapper (`ISSUES: oracle-cache-key-hashes-a-thread-count`).
The panel is not bit-reproducible at a pinned seed: `noop` matched `base` to ≤ 9 fragments per condition and two
direct `base` runs of one condition differed by 67, all below the floor.

### prior-fidelity-vs-deliverable
CLOSED by measurement 2026-09-19 (`ISSUES: end-to-end-error-unattributed`): the anti-correlation is gone on this
tree. A perfect prior now improves the deliverable in every in-scope stratum at both axes — `oracle`/`base` 0.949 /
0.972 / 0.959 at transcript level and 0.899 / 0.945 / 0.876 at gene level on unstranded OFF / stranded OFF /
stranded ON — where on 2026-09-10 it made both capture-OFF strata worse (1.019 / 1.026 transcript, 1.087 / 1.169
gene). The one in-scope row where it still reads worse, `g05 ss.99 OFF` (+2,013), is below its floor (3,722). The
leading answer (the messages destroying a nearly right self-solve at the worst slots) was never confirmed on a
second stratum and is no longer needed. `quant_accuracy.py`.

### scan-thread-split-starves-the-workers
CLOSED by landing 2026-09-19: the budget is split by the measured RATIO — one BGZF decompression thread keeps
about eight scan workers fed — instead of reserving a fixed four, so `bgzf_threads` defaults to deriving
`total // 8` and an explicit `--scan-bgzf-threads` still overrides. The old rule picked the worst measured cell
at every budget. Scan seconds on the 18.6M-fragment library, two interleaved rounds each, by (budget,
decompression threads): 4 — (0) 32.0, (1) 41.2; 8 — (1) 25.8, (2) 26.2, (0) 28.3, (4) 34.5; 16 — (2) 17.6,
(1) 19.6, (4) 20.1; the 2026-09-11 table it reproduces is (3,1) 113.5, (2,2) 59.5, (1,3) 42.4, (0,4) 33.7 at
total 4 and (4,4) 34.4, (2,6) 26.4, (1,7) 24.8 at total 8. BIT-IDENTICAL on all three identity references
including the real-library one that runs the scan, which is structural rather than lucky: every accumulator
bank is a sum of integers, and `test_scan_order_independence.py` holds the tally identical at 1, 2, 4 and 8
workers. `profiling/profiler.py --scan-only`.

### ruler-multimapper-floor-caps-the-correction
CLOSED by landing 2026-09-16 (`DESIGN.md` §7.2, `EQUATIONS.md` §11): the expectation ruler on the per-base
length — a transcript's own bases at their pieces' capture efficiencies with the fragment-length end taper, no
boundary or junction object in the length; a piece's efficiency the posterior mean of its clipped gDNA density
under the fitted landscape from its own contained count and the crossings at every boundary within a fragment's
reach, each crossing apportioned to the pieces its fragments cover by geometry times their own-count densities,
one pass — in `capture_efficiency` (computed by `calibrate`, published on the result), `capture_eff_length` and
`assemble_priors`; the multimapper floor `w = C/(C+1)`, the splice-junction objects and the flank imputation
deleted. The defect: the floor pulled every factor toward 1 by `1/(C+1)`, a +3.4 to +3.7 nat bias on the
unprobed class at every gDNA level; a piece shorter than a fragment had no contained support and read 0 then
the floor; forty junction objects imputed from unsupported flanks swamped four hundred bases of measurement on
the ladder. Killing numbers (`ruler_vs_truth.py`, the unprobed class's median log error; the test chromosome's
`g05 ss.99 / g25 ss.99 / g50 ss.50 / g50 ss.99 ON` rows): shipped +3.41 / +3.32 / +3.74 / +3.76 → landed +0.20 /
−0.03 / −0.15 / −0.03, the probed class 97–98 % → 99 % within ±0.1; the tiny-exon block's `capmixed` −0.63 →
+0.07 and `mixed` +5.9 → +0.2 (the unprobed `tiny` +7.0 → +0.7 / +1.9 at 5 % gDNA, where its edges hold 0.04
expected fragments and the posterior is the prior's valley); the depth ladder's probed class 75 / 93 / 97 % →
86 / 99 / 99 % at a tenth of the depth and 63 / 62 % → 91 / 91 % at a hundredth; the ladder's `g05 ss.99 ON`
unprobed class +4.40 → +0.45. Refused with numbers, each open question tried each way: the solve's own
posterior (+1.25 / +4.96 nat unprobed), the clip outside the expectation (indistinguishable), the own count
alone (+4.5), the plug-in (+3.36), the apportionment iterated to convergence (identical on the test
chromosome, +0.44 against +0.46 on the ladder), a joint update of neighbours (never converges). The
falsification tests were verified failing on the shipped ruler through its own fixture: an unprobed exon with
no gDNA fragment read the floor at +3.19 nat over the depleted level, a transcript of 40 bp exons read 0 and
then the floor. The multimapper blindness the floor stood for is `ISSUES: multimapper-blind-support`.
NOTE 2026-09-23: the per-base length and the crossing apportionment this landing built on are replaced
(`ISSUES: the-per-base-rule-for-every-component`, `ISSUES: the-crossing-apportionment`); the floor's deletion and
the posterior-mean efficiency stand. The own count alone, refused here at +4.5 nat because a piece shorter than a
fragment had bases to price and no count, now ships: under the one shared rule such a piece has no share, so its
efficiency multiplies nothing (`ISSUES: the-crossing-apportionment`).

### the-ruler-reference-on-sparse-real-libraries
CLOSED by landing 2026-09-16 (`DESIGN.md` §7.2, `EQUATIONS.md` §11): a basin's members are the kernels with a
location (count ≥ 1, `DensityLandscape.located`), the enriched candidate is the basin above the depleted one
with the most located members, and it is a mode iff more than `k = √n_located` members resolve it at the
population's k (the median width to the k-th nearest member ≤ 1 nat); the regime is on the result
(`gdna_reference_members`). The defect: `located_enriched_mode` tested membership on every kernel centre, and
a zero-count anchor's centre is its resolution wall `1/E` — a human index trains ~250–325 k anchors whose walls
span every decade, so on a sparse library a basin above the bulk was "located" by walls: LBX0190 (3,434 gDNA
fragments on regions, 1,118 located kernels) chose 10^-1.57/bp from 16,931 members, 16,918 of them anchors and
10 located, while its probed level (599 fragments in 3 regions near 10^-0.5) is no population; MO_3021
(65,874 / 15,088) chose 10^-1.35/bp from 20,053 anchors and 16 located members. Killing numbers, landed: both
panels' 46 rows and the depth ladder's 26 capture-ON rows unchanged to the reference (the enriched modes hold
257–3,606 located members); LBX0190 and MO_3021 read `None` (no basin above the bulk holds more than 11 / 16
located kernels against k = 33 / 123); LBX0588 (33,943 located, 11,579 members at 10^-0.53) and the VCaP
library (74,968, 24,496 members at 10^-1.07) unchanged; two capture-OFF rows of the depth ladder at a tenth of
the depth (`g00`, `g001`) that had read a reference from one located kernel among 37–43 walls read `None`,
restoring the capture-OFF contract there. The choice rule (largest rendered mass, most located members, largest
located weight, highest located basin) agrees on every row measured once the members are located, so the
most-located-members rule ships by derivation, not by number. Falsification test verified failing on the
shipped reader (a fitted basin of 600 short anchors' walls around twelve located kernels read a located mode of
612 members) and the perturbations fired: walls admitted as members, the k rule dropped, the members read at
their own √n_members. The three `review_identity_*` references re-frozen (LBX0190's reference is `None` now,
by design). The measured operating boundary stands: the reference is `None` below about 120 gDNA fragments on
the test chromosome, and about √n_train located probed pieces at one fragment or more are needed for a mode.

### capture-on-strand-pure-ambig-undercall
CLOSED by landing 2026-09-14 (`DESIGN.md` §6b.15.13, `EQUATIONS.md` §9f): the tilt atom — the AMBIG tilt's
hypothesis space {pure +, pure −, mixed} at equal reference weight, two columns at `τ = ±1` in ψ's cube beside
the continuum (`−log π` on its trapezoid weights), a held RNA level on a strand ruling the other strand's atom
out. The mechanism: the continuum's median sat below the strand cap at a strand-pure slot (prior-free, a truth
of 0.50 read 0.31–0.37; 0.46–0.50 now). Killing numbers, the L3 tree → landed: ladder stranded ON
470,862 → 427,046 (`g50 ss.99 ON` 217,636 → 199,409, `g98 ss.99 ON` 176,468 → 149,603), every stratum
better, 15 of 16 contaminated rows better, `g05 ss.99 ON` +1.7 %; stranded zero controls 497 → 550 / 224 → 231,
unstranded within a fragment; test chromosome −0.1 / +0.2 / −0.06 / +0.09 %; the strand-pure census band
−28 % / −38 % on `g50` / `g98 ss.99 ON`, its tilt error 40–80 % lower on every stranded row; the spliced stress
identical on every both-strand row; the encompassing region between TA+'s exons 0.366 → 0.511 against 0.544.
The source reproduces the prototype exactly on the test chromosome and to ≤ 0.007 % on the ladder. Gates in
`test_vertex_reference` (four perturbations fired). Not implemented: the structural witness (within 2 % of the
delivered-level witness). The cost where no witness exists is `the-atom-at-an-unwitnessed-both-strand-slot`.

### deadband-gates-a-gdna-free-library
CLOSED by landing 2026-09-14 (`DESIGN.md` §6b.15.12, `EQUATIONS.md` §5.2b): the strand channel is live iff the
protocol preserves strand — the Bayes factor on the spliced 2×2 of a free κ (the fit's own Beta(1, 1)) against
κ = ½ exactly, closed form, no constant (`region_init.strand_discriminability(κ̂, N)`); `disc = 4(κ̂−½)²` where
live; `n_gdna_obs` deleted throughout. Killing numbers: the replaced form (an unbiased estimate of `(κ−½)²`
floored at zero) is positive on 32 % of unstranded libraries — `g98 ss.50 OFF` shipped LIVE at z = 1.22, and
`g00 ss.50 OFF` was dead only by the gDNA term, without which 499 → 21,484; all four ladder `g00` rows had
`N_gdna = 0` and the channel dead at κ = 0.0099, so every g00-donor toy of the θ thread had run without a strand
channel. Landed: the ladder identical everywhere but unstranded OFF −0.22 % (`g98 ss.50 OFF` 123,657 →
122,981); the test chromosome identical on 30 rows; two contaminated rows bit-identical; the stranded zero
controls 405 → 497 / 194 → 224, located at the refit rung (the landscape's vertex bias,
`gdna-landscape-trains-on-false-positives` (d)); the encompassing exon∩exon slots on a gDNA-free donor
0.054 / 0.046 → 0.002 / 0.005; 17 goldens regenerated (≤ 2e-3 relative, `antisense_contained`'s false gDNA
78.7 → 5.6). REFUSED with it: an RNA level read from a slot's belief — the prior's answer re-emitted as data,
the relay (`TRAPS: one-hop-lifted-out-is-still-the-relay`). Gates in `test_region_init` (five perturbations
fired) and `test_encompassing_locus` on a gDNA-free donor (retired 2026-09-28).

### measured-prior-rung-4
SUPERSEDED 2026-09-14: ψ's composition reference fitted from composition-free observables (`rho_0`, the
enrichment responsibility, the detector as a boolean) was the measured-prior thread's fourth rung
(2026-08-26). The gDNA landscape hyperprior (`DESIGN.md` §7.1, trained on the previous solve's deconvolved
gDNA, the zero controls at a few hundred fragments) is the population prior ψ carries, and ψ has no reference
location (`DESIGN.md` §6b.1); the total-density landscape it would have consumed is QC-only. The requirements
it listed (span both lattice ends, exact at `g00`, one pseudo-fragment) are met or moot. Re-open only with a
measured gap the landscape prior cannot close.

### reference-prior-refuted-at-concept-level
RECORD (2026-08-24; the ruling is `DESIGN.md` §6b.1): ψ's reference tilt was refuted at the concept level — on
λ the data's information `I ∝ N_eff·disc·[f_g(1−f_g)]²` is zero at κ = ½ while a tilt holds fixed nats there,
overturned at N = 3 where the strand channel is alive and never at κ = ½ (0.7471 from N = 10 to 10⁶); 82–95 %
of unstranded pass-0 error sat on slots it decided. Refused: another constant (0.75 optimal on none of 16),
`σ(L)` (3.3×/10.8× worse at the zero controls), re-weighting (up to 2.8× worse), an information-weighted tilt,
the arcsine coordinate (`ISSUES: arcsine-magnitude-coordinate`). The reference location is deleted.

### ambig-node-as-a-gdna-source
REFUSED as built 2026-09-08: a both-stranded node's own counts plus its held RNA levels, read as a gDNA level
and emitted lower-sided — `g05 ss.70 OFF` +2.6 % on all three test panels, +4…+38 % on the weak-κ zero
controls, ladder `g05 ss.50 OFF` 1.745× / `g05 ss.99 ON` 1.465× (at low gDNA the level's mode is noise and
travels as a floor); its target, the walled host exon of a `span` locus, gains ~100 fragments. Re-open only
with a gate on the emitted level's own width.

### psi-lambda-bracket-unshipped
CLOSED 2026-09-14, landed with the one-lattice ruling (W5, `DESIGN.md` §6b.15.10): the λ bracket follows the
landscape prior's derived demand (`DensityLandscape.required_logodds_window`, read by `calibrate._sweep`) and
the point count follows the bracket at the fixed step; nothing ships off.

### the-cancelling-pair
REFUSED twice (2026-08-26 the last, with the measured intron reference as load): `struct_lock` rescoped to
`g1_locked ∧ REGION` and the `intergenic|exon` boundary claiming its RNA-contaminated crossing mass as gDNA —
neither half prices alone (`TRAPS: a-cancelling-defect-pair`), worse on two of three in-scope strata with the
wins confined to `g00`. The one-sided certified-RNA bound is the only mechanism the zero control endorsed on
every row (−81.9 %, 8/8) and is panel-negative alone. Revive only with messages on at `g05 ss0.50 capture_on`;
the confusion matrix is re-derivable from `build_structural_claims` (deleted unread 2026-09-22; in git)
against `slot_truth.npz`.

### parked-capture-pilot-sign
MOOT 2026-09-14: two capture-ON pilot rows disagreed about the sign of every length correction (2026-08-13);
both panels were deleted and the correction lives in the length composition channel, retired until after
0.8.0. If the channel returns, find which row lied rather than averaging.

### strand-marginal-volume-factor
REFUSED 2026-09-13, both forms, with their numbers (the derivation and the arms: the sandbox's θ note §11). The
finding stands as a property, not a defect: the exact θ-marginal of the strand term at an interior tilt is
∝ σ_τ(λ) ∝ 1/(1 − f_g), an Occam factor of order √n toward gDNA at a balanced both-strand slot — under a
proper prior on the tilt it is intrinsic (any normalised prior gives it; "uniform in p" is uniform in τ and
keeps it), and it is what the strand channel genuinely says: balanced strands are explained by gDNA with no
parameter. Two ways to remove it were priced on the ladder (Σ|Δ gDNA| both axes, g00 excluded; reference
249,703 / 471,202 / 308,529 / 3,622,383 on stranded OFF / stranded ON / unstranded OFF / deferred):
* the tilt's Jeffreys volume, `+ log(1 − f_g)` on the AMBIG cube (the AMBIG reference Beta(½, 3/2)): 267,373 /
  761,711 / 331,064 / 5,263,361 — +7 / +62 / +7 / +45 %; the g00 rows 2–4 % better; the test chromosome +1–3 %
  on every stratum. Where: the gDNA-rich AMBIG slots under capture (g98 ss.99 ON strand-pure band 34,056 →
  120,845, g50 38,823 → 81,151) — where `a(λ)` is small the strand term does not constrain the tilt and the
  term is a prior shift toward RNA on slots that are mostly gDNA;
* the PROFILE form, the strand term entering the λ posterior through `max_τ L` (the cap, flat inside it), the
  tilt integrated only for the rows and the read-out: 253,117 / 687,783 / 308,557 / 3,653,341 — +1.4 / +46 /
  0 / +0.9 %; g00 marginally better; the test chromosome +0.2 / +0.8 / 0 / 0 %, the same g98 capture-ON bands.
On the mono shared-exon toy (no junction, nothing guards the slot) both forms do what they claim — a balanced
500k exon's own solve 0.997 → 0.49 — and on the spliced toy neither moves a number: the factor bites only
where nothing else informs the slot, and that slot's honest answer under a flat strand marginal is the
reference median under the cap, not zero. What survives: the width factor at gDNA-rich AMBIG slots is
evidence the panels reward, so the measure stays the arcsine and the marginal stays exact. The saturating
variant (normalising by the truncated volume) was refused on paper: it inverts the strand-pure tail. The
exposure — a deep, balanced, junction-free AMBIG exon in a gDNA-free library reads mostly gDNA — is recorded
here with the mono toy as its instrument (`deep_stress.py --mono`), and is the RNA level lanes' business (a
junction on either strand pins it), not ψ's.

### theta-quadrature-at-zero-gdna
CLOSED by landing 2026-09-13 (`DESIGN.md` §6b.15.11; the derivation `EQUATIONS.md` §9e): the θ nodes follow the
strand term's peak (`simplex_logodds._tilt_window`), the node count derived (24 = 2T/π + 1, T = −log ε₆₄). The
recorded mechanism was wrong: the K_t 30 failure was ONE slot (`g00 ss.99 ON`, slot 37345, n 25,242, 28 % of
its RNA on the minor strand, a 0.006 rad peak) — the fixed lattice's sum is a comb across λ on any deep
interior-tilt slot, and 60 held the control by landing a node on that one peak; the strand-purity quartic was
never it, and the trapezoid endpoint weights alone (REFUSED, ≤ 0.7 fragments) removed nothing because the
resolution was the whole term. Killing numbers: the marginal's λ-shape error 90–130 nats at 500k under the
lattice at 60 (0.1–0.4 at 2,400 nodes), 2·10⁻⁶ under the window at 24; the shared-exon stress at 500k, lattice
→ window: `g50` balanced 102,076 → 19 false fragments, 20 %-minor 11,448 → 11, `g00` 10,994 → 632 and
9,907 → 25, tilt error down 6–70×; the ladder a numeric near-no-op (stranded OFF / unstranded unchanged,
stranded ON +0.4 %, `g00` rows identical, the both-strand tilt error at 24 nodes what the lattice reached at
120), the test chromosome within ±3 fragments. Gates: `test_vertex_reference` — the marginal against adaptive
quadrature at 500 and 500k fragments, the derived count converged, a delivered row evaluated at the nodes to
the bit, κ = ½ the whole domain; five perturbations fired. The second step the same day removed the last θ
lattice: the lanes deliver a row's ingredients (`CubeRow`) and `sweep_n_tilt` is deleted — no tilt count
exists in the tool. What it exposed is `ISSUES: strand-marginal-volume-factor`.

### f32-strand-tilt-at-half
CLOSED by landing 2026-09-12: the AMBIG cube is float64 like the rest of ψ — the float32 cube was a memory
choice the tiling made moot, and one solver (`simplex_logodds._solve_logodds`, a single-strand slot the
``K_t = 1`` case) has one precision. At κ = ½ the strand mean is ½ identically in float64 and `w_pos` reads ½.

### landscape-trains-on-real-substrate
CLOSED 2026-09-10, superseded: the no-evidence share of the prior's training mass is now zero by construction
(`RegionBelief.has_composition`, gated); `landscape_training_census.py` reports the population per evidence class.

### the-landscape-training-population-arms
A/B'd on the test chromosome (30) and the ladder (16), halves apart, both zero controls on every row
(2026-09-10): two landed, six refused. Do not rebuild a discount that reaches the delivered rows, nor an
E-step on counted kernels. Whole-library |gDNA − truth|, ratio to the previous default; `noop` byte-identical.

| arm | what | ladder zero controls (ss.50 OFF / ON / ss.99 OFF / ON) | in scope | deferred | verdict |
|---|---|---|---|---|---|
| `no_echo` | a slot with no evidence, or a bound only, does not train | 0.742 / 0.753 / 0.962 / 0.966 | within 0.1 % on every row | 0.98–1.01 | ⭐ LANDED (`DESIGN.md` §7.1) |
| `no_echo_1sided` | also a one-sided delivered composition row | 0.512 / 0.660 / 0.748 / 0.910 | ≤ 1 % either way | ⛔ 1.90 / 1.28 / 1.20 | REFUSED: the probed exons' one-sided rows are the enriched mode's witness |
| `own_var` | the own-evidence variance (`tau_lam` + the row's precision) as `_reliability`'s `v` | 0.737 / 0.742 / 0.956 / 0.959 | ⛔ `g05 ss.50 OFF` 1.058; test chromosome `g05 ss.70 ON` 1.24, `g50 ss.70 ON` 1.12 | ⛔ 2.20 / 3.27 / 3.57 | REFUSED |
| `sweep1_no_echo` | `no_echo` at the first fit only | 0.890 / 0.892 / 0.994 / 0.995 | 1.000 | 1.00 | REFUSED: weaker than the rule at every refit |
| `shuffle` | `no_echo`'s excluded count on random non-anchor slots | 0.914 / 0.789 / 0.978 / 0.968 | 1.000–1.002 | 1.01–1.07 | the control: wins the zero rows (every non-anchor slot is false there) and nothing in scope |
| `refits6` | `calib_refit_iters` 3 → 6, no patch | 0.822 / 0.805 / 0.987 / 0.991 | 1.000; test chromosome unstranded ON 0.85–0.91 | 0.96–0.98 | NOT TAKEN: 0.725 with the rule against 0.742 without, at twice the refit cost (`ISSUES: performance-memory-bounded-solve`) |
| `row_kernel` (unit weight / reliability weight) | a delivered-row slot's own λ-likelihood, mapped to the density grid, as its kernel | test chromosome only: `g00 ss.50 ON` 437× / 0.95 | `g25 ss.50 OFF` 5.40 / 1.57; `g25 ss.50 ON` 2.10 / 0.95 | `g50 ss.50 ON` 1.62 / 1.00 | REFUTED: the sum of per-slot likelihoods is the flat-prior posterior, not a density; a bound spreads uniform mass under itself |
| `oracle_values` (weights as trained / var 0) | the CEILING: the population trained at certified gDNA | 0.686 / — / 0.956 / — (var 0: 0.749 / — / 1.113 / —) | `g05 ss.50 OFF` 1.016 / 1.061 | `g50 ss.50 ON` 0.964 / 0.753 | a price, not a mechanism: the values fix recovers most of it |
| `estep_zero` (ratios to the landed population fix) | the E-step on the kernels with no location (count < 1): kernel × the previous refit's landscape | **0.006 / 0.002 / 0.038 / 0.018** | unstranded OFF 0.952 / 0.973 / 0.988; stranded OFF 0.986 / 0.994 / 0.997; stranded ON 0.980 / 1.004 / 1.004 | 0.960 / 1.020 / 1.011 | ⭐ LANDED (`DESIGN.md` §7.1); the test chromosome wins or ties all 30 rows |
| `estep_all` | the E-step on every kernel | 0.006 / 0.002 / 0.038 / 0.018 | `g05 ss.99 OFF` 1.000, `g50 ss.99 ON` 1.018, `g98 ss.99 ON` 1.047; test chromosome `g25 ss.50 OFF` 1.091 | 0.956 / 1.045 / 1.012 | REFUSED: a counted minority competed away by the bulk, the δ-pin EM predecessor's failure |
| `estep_mirror` | the control: the previous landscape reversed on the grid | test chromosome capture-ON zero rows 1.48 / 12.7 / 8.9 | neutral on capture-OFF rows (the mirrored bulk lands in the depleted zone by the grid's asymmetry) | — | the control discriminates where it can: the placement is what acts |

### data-derived-reference-location
REFUSED 2026-08-24: the intron-reference estimator's selector widened to `solvable_exon` wins all 8
capture-OFF rungs and all four `g00` rows exactly and loses all 6 contaminated capture-ON rungs. Priced on all
16, deleted.

### flux-source-skipped-at-an-empty-exon-piece
CLOSED by landing 2026-09-09: `lanes.rna_lanes` builds a junction's flux level at an empty exon piece too,
priced by `hop_price` on the piece's zero count; ladder through the pipeline every in-scope row within 0.5 %
(worst `g98 ss.99 ON` 1.0048×), stranded zero controls 0.977×/0.958×. Only the ladder can judge it.

### the-empty-flux-source-at-the-junctions-counting-alone
REFUSED 2026-09-09: the sharper price for the empty-piece source (the junction's counting alone) wins the zero
controls more (0.944×/0.970×) and costs `g98 ss.99 ON` 1.0179× (+3,144 at the AMBIG exon|exon boundaries) and
`g50 ss.99 ON` 1.0091× — a sharp RNA floor reaching a node whose price reads the column count
(`ISSUES: flux-price-witness-units`).

### flux-witness-in-strand-units
REFUSED 2026-09-09: the column split's asymmetry as the flux price's witness reads worse at pass zero on the
test chromosome (`g05 ss.50 OFF` +12.6 %, `g05 ss.99 ON` +20 %) — it reads zero at an equal-abundance overlap
exon.

### relay-od-r-discontinuity
CLOSED with the deleted relay policy, 2026-09-09: `g98 ss0.50 capture-OFF` read 217,531 at `od_r ≤ 1e−7` and
212,581 at 1e−5, a threshold not a response (`TRAPS: a-constant-parked-a-value-off-a-knife-edge`). Standing
half: a 1e−5 nudge moved 1/30 rows more than 0.5 % (worst 1.65 %), so a single-row policy difference below
~2 % needs the noise floor re-recorded in the same session.

### levels-always-travel-for-the-gdna-lane
REFUSED 2026-09-08 (three forms): the gDNA lane on every directed face wins every capture-ON ladder row
(unstranded × ON −15 to −17 %) and loses `g05 ss.99 OFF` +7.1 % (45,076 → 48,295) and `g05 ss.50 OFF` +7.8 %
(50,948 → 54,926) — a one-sided floor at a node with no channel of its own is a tilt; with phase 2's ceiling
worse (`g25 ss.50 OFF` 1.59×, the weak-κ zero control 1.8–3.4×). Re-opened only by
`ISSUES: two-sided-exon-row`.

### the-edge-upper-side
REFUSED 2026-09-04, three upper sides on the intergenic|exon edge's level: undampened (`g98 ss.99 ON` +75 %,
`g98 ss.70 ON` +360 %); dampened both sides (sparse panel `g98 ss.99 ON` 36,645 vs 6,981); dampened above only
(no in-scope win, stranded 0/8, deferred 1.065–1.146×). No local witness prices the interior's enrichment.

### the-abundance-discrepancy-map
LANDED 2026-09-02, REPLACED by the level rule 2026-09-04: a fitted step between hypotheses is a pooled
premise; removing it improved `g50 ss.50 OFF` 12,328 → 11,044 and `g25 ss.50 ON` 28,414 → 19,854. Also
refused: a Gaussian summary of the held profile (12,328 → 13,911); the form built from what the boundary holds
(11,503). A level is made from the sender's measurement, never from what it holds; removing the dampening
costs 47 % on `g50 ss.99 ON`.

### the-certified-flux-row-as-a-level
REFUTED by probe placement 2026-09-04: a two-sided flux-implied composition row at certified faces wins the
exon-probed test chromosome (`g50 ss.50 OFF` 11,242 → 8,876) and fails where the probes move — junction-probed
`g50 ss.99 ON` 6,502 → 26,497 (spliced fragments enriched 20–35× more than contained ones), ladder zero
controls 3.3× (`g00 ss.50 OFF` 299,380 → 995,880) and 6.9×. A hard s ≥ 1 truncation refused alone
(5,712 vs 1,803).

### two-sided-exon-row-forms
REFUSED 2026-09-03 on the test chromosome: a two-sided Poisson row fixes capture-OFF (`g50 ss.50 OFF`
11,242 → 8,945) and is catastrophic capture-ON (`g50 ss.99 ON` 6,367 → 15,558); an abundance-bounded row still
fails `g50`/`g05` ON (12,333 / 3,090 against 6,367 / 2,293); a flux cap claims 19 % gDNA at the zero control
(16,819 → 121,263).

### the-discrepancy-priced-cap
REFUSED 2026-09-04: a cap on the face map's plateau priced by the pair's discrepancies fixes the first pass
(−20 to −24 %) and loses the junction-probed panel's weakly stranded rows (`g25 ss.70 ON` 1.122×).

### the-two-sided-level-lane
REFUSED 2026-09-05: with every node reached, a two-sided level lane reads −26 % at pass zero on
`g50 ss.50 OFF` and +33 % on `g50 ss.99 ON` — an exon complex's minimum total density bounds its gDNA from
above only off capture.

### the-wall-above-the-face-map-ceiling
REFUSED 2026-09-09: a wall above the face map's ceiling priced by the junction–exon pair closes the toy
harness gate (|Δf_g| 0.123/0.126 against 0.848) and wins the ladder's target rows at pass zero
(`g05 ss.50 OFF` 0.921×, `g50 ss.50 OFF` 0.847×), and through the pipeline loses the stranded half 0/6
(`g50 ss.99 ON` 1.142×), the deferred rows 1.15–1.43×, the junction and sparse panels' capture-ON stranded
rows 1.4–2.5×.

### the-pooled-hop-step
REFUSED by the owner 2026-09-03: a per-hop-kind pooled disagreement premise. The per-pair rule alone stands
within 0.25 % of it on every ladder row and ahead on four of six stranded rows (`g05 ON` −1.34 vs −1.38 %,
`g98 ON` −2.70 vs −2.49 %); where pooling wins (junction panel `g05` +3.5…+5.6 % vs +8…+9 %) is a systematic
offset the owner declines to extrapolate from other pairs.

### the-flux-factor-hop-premise
REFUSED 2026-09-02: fitting `1/delta` on the flux at the alternative splice site's E hop. Alone it fits
2.07 ± 0.26 on `g50 ss.99 ON` and harms (boundaries 231 → 280; sparse panel 394 → 562); jointly with the step
the two are degenerate at n ≈ 14 (sparse k = 2.59 ± 0.13, χ² = 47.7 for 14 pairs; boundaries 382 → 526).

### arcsine-magnitude-coordinate
REFUSED 2026-09-01: `f_g = sin²φ` in place of λ = logit(f_g). The premise is false (a posterior median cannot
reach a vertex in any coordinate; logit is finer there — top f-gap 1.83e-5 vs 1.37e-3 at K = 60), the leverage
is bounded at ~2.5 % of the vertex-band defect, and the ladder A/B loses every in-scope unstranded × OFF row
(1.107/1.071/1.014), winning only stranded × OFF (1.008/0.939/0.898). Commit `48be3d0a` holds the prototype.
Vertex information re-priced under the shipped policy (2026-09-09, `vertex_ceiling.py`): ≤ 1 % on stranded
in-scope rows, 2–7 % unstranded OFF, 25–35 % deferred.

### refused-transcript-weights
REFUSED (all 36 conditions, seed pinned): a soft-min per-transcript allocation over exclusive objects with a
per-object Jeffreys half — worse on every stratum and the zero control, transcript Σ|err| 57.5 M → 81.6 M
(1.42×; length-proportional 2.10×), the gDNA:RNA split +0.2 % so the allocation alone was priced. Exclusivity
hard-zeroed 38.7 % of transcripts; a density from < 200 bp amplified variance up to 6,534×; the +½ revived the
silent half (false-positive mass 18.6 M → 41.6 M). Deleted; the lane kept
(`ISSUES: per-transcript-prior-lane`).

### refused-soft-min-path-weighting
REFUSED: the owner's path theorem, 12 arms on `g00 ss0.99 capture_off` and 3 on `g50 ss0.50 capture_on` —
worse than base at transcript level on every rung (1.317–1.604× at `g00`, 1.262–1.331× blind), false-positive
mass 1.76–2.20× worse; the gene-level 0.395–0.527× at `g00` collapses to 1.006–1.041× blind. Mechanism
`TRAPS: an-upper-bound-is-not-an-estimate`: 3,644 of 4,839 silent transcripts (75.3 %) share an object with an
expressed one. Kept: the dial is monotone, 0.0 % of expressed transcripts zeroed. The instrument was retired
with the refusal (2026-09-10); git carries it.

### the-rna-length-law-fix
REFUTED 2026-08-31: `rna_pmf` is sound as shipped. The −4 bp sj-observability diagnosis was measured on the
undrained pool; on the drained payload the residual vs mature truth is −0.24/−0.11 bp off capture, −1.06 bp
under it, and dividing by `pi·o(w)` overcorrects +17 bp. The −17,754 transcript win from
`calibrate(rna_fl_pmf)` at `g05 ss.99 ON` is compensation (shifting the shipped shape +10 bp wins −171,342,
far past truth) on the EM's short-exon isoform misassignment; on the 0.8.0 metric the true pmf moves
±0.3–2.9 %, mixed sign. `r_hat` is accurate off capture (−0.04/+0.19 bp), broken under it (−11/−20 bp). Never
rank a calibration input on the thermometer.

### the-general-verdict-table
Mixed verdicts of 2026-08, each with its instrument. "The RNA fragment-length model" is the accumulator's fl
geometry; the length-channel retirement is of a calibration composition channel.

| | closed by | verdict |
|---|---|---|
| **the gDNA scale rule** · **the mass rescale** · **TSS/TES as the population licence** | landed 2026-08-04 | ✅ `EQUATIONS.md` §3.5/§3.5b/§3.5c (the gates named here retired with the relay policy, 2026-09-09). ⚠ The ceiling says the mass rescale cost the panel **nothing** (+0.0002 to delete it outright); it landed on the derivation and on being free |
| **face (I) of the `intron\|exon` BOUNDARY** | re-solve ceiling + panel arm | ⛔ **DO NOT BUILD.** The derivation (`EQUATIONS.md` §3.6) is re-verified and is not what failed: handing both BOUNDARIES the ORACLE truth and re-solving is worth **−0.000** off capture, and the ladder prototype is **negative** (mwae 0.0413 → 0.0426, confidently-wrong +10.7 %). TRAPS: panel-before-src |
| **a LEVEL transfer from the intron** | toy + panel | ⛔ **REFUTED**, +0.207 on capture-ON × unstranded — capture inverts which side is well-counted (TRAPS: capture-inverts-the-counted-side) |
| **the RNA fragment-length model** | a per-pmf ceiling arm (`length_ceiling.py`, retired 2026-09-10) | ⛔ **−0.02 %** at pass-0, **+0.21 % (worse)** over all objects. Root cause exact (`pi(w)` scores sj *crossing*, the pool requires the splice to be *seen*). ⭐ Its value is the BOUND: the whole fragment-length-model cluster costs ≤0.43 % of the shipped solve. TRAPS: price-the-halves-separately |
| **TRAPS: pure-and-length-censored's κ residue, as an ACCURACY fix** | κ injected at exactly ½, all 36 conditions | ⛔ **−0.2 %** unstranded, worse on the shipped solve. ⭐ But the *general* defect — a boolean licence flipped by a small residue — is **the-capture-level-residual**, and the destruction control taught TRAPS: honesty-metrics-reward-ignorance |
| **a nascent-bearing ladder condition** | toy, 36 conditions × 7 rungs | ⚠ **−5 %**, and the wrong way on one stratum. It was a toy harness arm (`--nrna 60`, retired 2026-09-28); it no longer justifies re-simulating the panel |
| **the gDNA prior's BIMODAL CAPACITY, and "give the prior more signal"** | a read of `gdna_landscape.py` + the production refit on real conditions | ⛔ **BOTH BRANCHES CLOSED.** The prior already renders the landscape correctly — **2.98 decades** of mode separation at `g75 ss0.99 capture_ON`, 30× more enriched mass ON than OFF, a single pile at the wall at `g00`. And a prior fitted from ORACLE truth is the same prior (0.04 dec). Not capacity, not signal, not location. ⭐ Why an evidence-free object cannot reach the vertex at all — and why that is the value of missing information rather than headroom — is `EQUATIONS.md` §9a |
| **the Jeffreys MEAN density location** | `--arm eta`, the `g00` zero control | ⛔ **REFUTED at +96,299 %.** It cannot say ZERO (`region_init.rho_g` is an exact 0 at 60,544/70,176 slots — the statement earning the −98 % at `g00`), and the TRAPS: a-ratio-cannot-carry-zero benefit it was credited with belongs to that fix, not to this arm. ⭐ If revisited the derived form is the Gamma **MODE** `max(a−½,0)/E`, which is exactly 0 at a zero count |
| **a threshold anywhere in the licence family** | TRAPS: a-threshold-on-a-fitted-residue implemented and refuted one | ⛔ τ is continuous across the region, so any floor is a tuned constant (TRAPS: a-threshold-on-a-fitted-residue, TRAPS: a-licence-with-no-floor, TRAPS: a-multiplication-gated-by-a-trace — refused three times) |
| **simulator captures a pre-mRNA through every probe its genomic blocks span** | NASCENT SCOPE RULING, 2026-08-22 | ✅ **WON'T-FIX** — nascent × capture fidelity is out of scope (`DESIGN.md` §0b); the sparse rebuild shrank the residual it explains (capture depletes nascent ~an order of magnitude), and no verdict depends on it. Re-open only if a real library forces it |

### the-reference-mean-family
Attempts to give ψ's reference a mean, plus two out of shipping it, each measured on the panel with both zero
controls and refused (2026-08-15/16); the surviving form is `DESIGN.md` §6b.1. Outside the table: giving `τ_λ`
the location term's curvature — bit-identical on the deliverable on all 32 rows, the 3,227× fall in `τ_λ` at a
pinned slot being ~98 % the `[f(1−f)]²` Jacobian (`TRAPS: a-priors-curvature-is-not-the-datas-information`); a
per-object one-pseudo-fragment floor `m_i = E[g]_i/(E[g]_i+1)` — worse on every stratum (0.609 / 1.045 / 0.580
/ 1.000 against 0.381 / 0.659 / 0.363 / 0.800), replaced by the lattice's top point `σ(L)` (`EQUATIONS.md`
§9c.1). A caution: `f_g ≤ 1 − S/M` from certified RNA is false (the truth violated it by 302) —
`boundary_spliced` is a separate bank, and the identity is `ρ_r·E_r = unspliced_RNA + S`.

| # | mechanism | why it was refused |
|---|---|---|
| 1 | **a fitted RNA density `logP_r`**, the mirror of the gDNA landscape | ⛔ the only non-circular form (fit from the solver's own belief) reads **0.988 / 0.997 / 1.037** — nothing, then worse: feeding ψ a density fitted from ψ's own belief tells it what it already believes. ⭐ And the ORACLE version's gain did not survive a shuffle — a shape that is wrong on purpose BEAT the true one at `g98` (0.786 vs 0.854), so the attribution was never established (`TRAPS: attribution-must-survive-a-shuffle`) |
| 2 | **a library-wide Beta mean, `a = f_lib`** | ⛔ `g05` regresses **1.43×**; `f_lib` is calibration's own output so the loop has positive feedback with both vertices attracting; and moving `a`/`b` sets the TAILS as well as the location, so `b = 0.03` leaves **57 %** of the prior outside `L = 10` |
| 3 | **the OBJECT-weighted mean instead of `f_lib`** | ⛔ **the two split by STRAND and neither wins everywhere**: object-weighted 0.584 / 0.452 on the two stranded strata but **5.570×** on unstranded × capture-OFF. ⚠ The sweep that motivated it was `ss_0.99` on three of its four rows — `TRAPS: never-pool-the-strata`, met on a sweep rather than a panel |
| 4 | **a stratified ASSERTION** — pure-gDNA strata claim `f_g = 1`, reweighted by stratum size | ⛔ reads 1.000 at every condition **including `g00`**: an assertion cannot see a library with no gDNA in it. ⭐ What replaced it is a per-object DENSITY, which needs no reweighting at all — the strata select the training set, not the answer |
| 5 | **a pooled RNA density from sj flux** | ⛔ RNA spans six decades with no genomic autocorrelation, so a pooled flux is not a population parameter (owner). It scored well only by sitting on the mass-weighted centre — `TRAPS: a-mean-hits-the-mass-weighted-centre-by-luck`. ⭐ Replaced by RNA-as-residual, which predicts no RNA at all |

### the-truncation-free-region-bank
REFUSED 2026-08-20: `region_start_count / ell` behind `CalibrationConfig.region_abundance_bank` improves the
pass-0 exon solve (no-evidence mass coverage 45.8 % → 66.1 % at `g98 ss0.99 ON`) and regresses the deliverable
— four zero controls 2.18×, deferred 1.84×, stranded × ON 1.20×, the only win 0.843× on stranded × OFF; under
the shipped bank 58.6 % of live exon hops carry a REGION bank of exactly 0, an accidental mute. Survives: the
truncation algebra (`EQUATIONS.md` §2) and the fill gate; the implementation is one commit before `a2b81b34`.
2026-09-27: the fill gate was deleted with the field it filled, `RegionGeometry.inv_abundance`.

### the-message-policy-campaign
Six mechanisms built, measured and refuted; closed 2026-08-27, the code deleted the same day. The bar it
missed is the one the transfer policy was later judged by — not to beat `SilentPolicy`, but to win on
unstranded data at minimal harm on stranded data, halves never pooled. The one measurement that survived:
propagation is net-harmful wherever the local solve has its own evidence, and its value concentrates where
that solve is blind.

| mechanism | what happened |
|---|---|
| **A composition-transporting policy** (`CurrencyPolicy`) | Best zero controls ever measured, but lost every in-scope contaminated stratum to silence. Deleted. |
| **The gDNA-continuity rule** (an unsupplied source's gDNA level crosses unscaled) | Built THREE ways — a static per-slot licence (a value RATCHET, gDNA densities to 3.9e+32, from breaking the knob's telescoping cancellation), a running-state licence (killed the ratchet, still lost capture-ON), and a fuse-based pure-gDNA re-anchor lattice (too weak — a fuse negotiates where the relay's mass rescale overwrites). Halved one row, lost others. Only safe beside a scan-time mass rescale. ⭐ The transfer policy's LEVEL LANE is this rule rebuilt as a one-sided profile with a priced hop (`DESIGN.md` §6b.12) — a different mechanism, judged separately. |
| **The premise's exon-end scoping** | The dispersion decomposition proves intron-end hops carry no COMPOSITION cost, yet scoping the charge to exon-end hops REGRESSED the panel: freeing intron chains before a measured LEVEL charge exists releases un-priced level drift. Restoring the pooled charge recovered `g98 ss.50 ON` 4.37 M → 2.44 M. |
| **A class-keyed method-of-moments fit on the observed log-ratio** (as the runtime law for the transport variance) | Tracks truth at intron and plain classes, REFUTED at sj classes: the route-summed flux cancels the visible step exactly where the true error is largest (0.08 observed vs 4.4 true). |
| **A totals-form pair fit** (for the same variance) | Refuted by construction — the knob consumes the totals, so transported totals agree by the (1−w) algebra and carry almost no information. |
| **`FanOutPolicy`** | Measured dominated once the certified-flux anchor gave its destinations own evidence. Deleted 2026-08-24. |

### the-doubt-graveyard
Eleven mechanisms priced and refused (2026-08-07). `g00` is the owner-required zero-gDNA control: its truth is
exactly 0, so every fragment there is a false positive with nothing to cancel it. The pattern is the result:
every one was a rule for resolving doubt, and at `g00` the doubt must resolve to no gDNA. The only candidate
the control has ever endorsed is the one-sided certified-RNA bound (−81.9 %, 8/8), panel-negative alone
(`ISSUES: the-cancelling-pair`); `zc_struct_lock_g1` is the one row still live, as half of that pair.

| candidate | what it did | `g00` | its target | why it died |
|---|---|---|---|---|
| `zc_jeffreys_mean` | `ρ_g = ½/E_g` at zero mass | ⛔ +7,269 % | −13.9 % | moves the mode UP |
| `zc_logmean` | `ρ_g = e^{ψ₀(½)}/E_g` | ⛔ +6,264 % | −11.3 % | moves the mode UP |
| `zc_anchor_mute` | no `prec_g` at empty locked slots | ⛔ +5,554 % | −7.7 % | kills the zero-gDNA win |
| `zc_struct_lock_g1` | scope `struct_lock` to `g1_locked ∧ REGION` | ⛔ +3,207 % | −1.2 % | ⭐ the MIS-SCOPED mask is load-bearing |
| `zc_reference_var` | `Var(f_g) = ⅛` where `τ = 0` | ✅ +0.0 % | −0.3 % | ⭐ passes the control and is INERT |
| `zc_discrepancy` | `+½ log D` shift, `(log D)²/12` | ⛔ +982 % | panel +4.5 % | moves the mode UP |
| `zc_disc_var` | the variance alone, mode untouched | ⛔ +255 % | panel +0.9 % | damping cannot bite |
| `zc_ref_prior` | own belief = ψ's reference, `τ + 1/π²` | ⛔ +3,792 % | −14.9 % | moves the mode UP |
| `zc_ref_prior_damp` | the two above, PAIRED | ⛔ +3,809 % | −15.5 % | ditto |
| the `eta` rebuild | a clean frame-free re-derivation | ⛔ unbounded | +85–103 % | see `DESIGN.md` §6.1 |
| the mean-location as a structural floor | the same idea in the LEVEL channel | ⛔ +96,299 % | — | it cannot say ZERO |
| `struct_lock = g1_locked ∧ REGION` **re-priced** | the standing strict xfail's own fix, on top of the SPLICE IN precision repair | ⭐ **0.71×** | in-scope **+2.5 / +2.1 / +0.5 %**, deferred +17 % | ⛔ **2026-08-18**, and the sign on `g00` FLIPPED from the 2026-08-11 row above once the relay was fixed underneath it — but the verdict did not: the four `nrna_none` zero controls go **1.8–15× WORSE**. The empty exons' "gDNA = 0 @ 0.2026" is LOAD-BEARING at AMBIG `exon\|exon` boundaries in an RNA-only library (an RNA+ level claim with no RNA− claim drifts ψ to 0.38 without it). ⭐ Still half of a pair |
| the mass rescale refuses a ZERO-MASS source (`pinM`) | `may_share_composition` additionally requires `M[src] > 0` | ⭐⭐ **0.31×** (134 k — better than the pre-licence relay's 154 k) | unstr × OFF 0.999, str × OFF 1.004, ⛔ **str × ON 1.032, six of six worse** | ⛔ **2026-08-18.** Under capture the empty slots between probe-covered stretches are the CONDUITS a relayed composition travels through (`TRAPS: the-divergence-was-a-barrier`), so the rescale at those hops is load-bearing. ⭐ The sharp predicate is NEITHER this nor the row above: refuse a source's OWN zero-count artefact, KEEP a relayed composition passing through an empty slot |

Four more `zc_*` arms exist as decomposition reverts used to attribute the 39 % win, not as proposals:
`zc_own_count`, `zc_live_count`, `zc_total_n` (inert) and `zc_transfer`, which reproduces the pre-fix tree.

### drain-provenance-split
REFUSED by arithmetic 2026-08-31: splitting the certified channel by drain provenance evicts 134,850 correct
drained records to remove 75 contaminants at `g50 ss.99 OFF`, resurrecting the −4 bp spliced-pool bias the
drain exists to repair (`ISSUES: drain-contaminates-certified-rna`).

### the-second-lambda-grid-and-its-regrid
CLOSED by landing 2026-09-13 (`DESIGN.md` §6b.15.10, W5 the grid study): one λ lattice for every consumer,
`sweep_logodds_step` 0.2; the fine single-strand grid, `_regrid_global` and `_scaled_grid` deleted. Priced on the
way and REFUSED, each with its number: a CUBIC regrid on λ (`g00` 892×, deferred 1.15 — the spline overshoots the
landscape prior's floor wall); a SMOOTH-RECONSTRUCTION read-out (a cubic or a monotone cubic of the log-posterior
on a 16× sub-grid: exact on Gaussians, 0.8 / 0.5 of a step on a balanced bimodal posterior against the histogram
quantile's 0.17, a wall overshoot on the cubic); a λ-axis LINEAR regrid (0.994–1.000 on the ladder: the axis was
never the harm, the interpolation was); a step DERIVED from the sharpest posterior (200–930 points where the metric
plateaus by 138, because an unresolved heavy slot costs ≤ n·f(1−f)·dλ/4 fragments and a fragment-budget rule needs
a tolerance); K_t below 60 (`ISSUES: theta-quadrature-at-zero-gdna`).

### rename-the-drain
RULED 2026-09-13 (owner): the term stays. Prepared and not pursued — the census counts (141 src / 146 scripts /
298 tests / 50 docs sites) and the candidates (`settle`, `decide`; `resolve` and `assign` collide) are in the
plan's §7 should the ruling ever be revisited.

### rename-row-and-face
RULED 2026-09-13 (owner): both terms stay. `row` (a slot's max-normalised log-profile over the solve grid;
1,236 src sites, two senses) and `face` (one directed side of a boundary, the `(destination, side)` pair a
rule is keyed by; 303 src sites) keep their names; the census and the candidates are in the plan's §7.

### g00-shrinkage-upstream-repair
CLOSED by landing 2026-09-14 (`DESIGN.md` §7.2, `EQUATIONS.md` §11): the ruler's reference is the located
enriched mode of the fitted gDNA landscape, published as `CalibrationResult.gdna_reference_density` and read
by `capture_eff_length` and `priors`; the kernel-density detector and its two constants are deleted. On both
panels (`calibration_vs_oracle.py` ③): the zero controls' factor 0.154 / 0.141 → 1.000 with nothing moved
(51,436 / 5,108 transcripts had moved); both in-scope capture-OFF strata P = O = 1.000 with nothing moved
(from P 0.957 / 0.970 and O 0.923 / 0.926 on the ladder, P 0.946 / 0.949 on the test chromosome — the
"exactly 1.000 off capture" contract line was not stale, the per-object clip of Poisson noise was the
defect); stranded capture-ON P/O 1.013 / 1.015 against 1.011 / 1.011, the reference within 4 % of the
truth's mass-weighted median on every capture-ON row measured; the deferred stratum 0.941 / 1.080 against
0.990 / 1.000 (one reference for both arms now, where two detectors had agreed by luck). The composition
was fixed first and the factor did not follow: 178–189 false fragments on 35,135 regions still made a
reference. The solve is untouched (`policy_benchmark.py` identical on both panels).

### u-ruler-arm
CLOSED 2026-09-14 with `g00-shrinkage-upstream-repair`: the "~2× loss on the two capture-OFF in-scope strata"
was the per-object clip, and with the reference a property of the solve the U ruler reads 1.000 by
construction (no reference) or `ρ̄/ρ_ref` against P's reference (0.06 on capture-ON, a number about
nothing). Its question — what a noise-free uniform field leaves — is answered structurally: nothing, the O
arm at capture-OFF reads 1.000 with no fitting. The arm is deleted from `calibration_vs_oracle.py`.

### flux-price-witness-units
CLOSED by landing 2026-09-14 (`DESIGN.md` §6b.13, `EQUATIONS.md` §12): the exon's witness of the junction→exon
price is its column count on the PROTOCOL'S SHARE of its RNA opportunity, `(c_s, κ_read · a_r)` — a column is
that much opportunity for the strand's RNA to be counted on it — so the junction's whole-strand rate and the
exon's density are in one unit and the pair's agreement is priced as counting alone at every κ (gated on a
hand-built chain at κ 0.99 / 0.7 / ½). The record: the golden `strand_ss65_multi_iso`'s gDNA-free exon reads
0.008 gDNA against 0.235. Both panels: the ladder in scope within 0.01 % on the metric and 0.01 % on the benchmark (stranded ON −0.01 %), the
deferred stratum +0.3 % / +0.07 %, the four zero rows identical (366 / 258 / 353 / 233); the test chromosome in scope within 0.2 %
(`policy_benchmark.py` −0.17 / 0 / −0.02 %), the deferred stratum −10 % on the benchmark and −8 % on the
metric, the zero controls identical. Three other forms priced and refused (`EQUATIONS.md` §12): the total
unspliced count (right unstranded, wrong at a both-stranded exon on stranded data — the ladder's
`g00 ss.99 ON` 233 → 435 on two AMBIG exons), the total as a one-sided bound (stranded capture-ON −8–10 %),
and the split's strand count at its own precision (the golden's ceiling 0.323, the refused
`flux-witness-in-strand-units` again).

### the-overdispersion-design-refusals
REFUSED 2026-09-28, in the derivations for `ISSUES: strand-overdispersion-one-shared-value`. Do not rebuild any of
these without the number that killed it.
ROOT AND WITNESS RULES:
- The two-witness random-effects rule (each seed set's smallest root, combined by a random-effects mean): the sign of
  a near-zero moment decides the junction witness, a hidden yes/no gate; a lone witness passes straight through (0.87
  once the gDNA-evidence fix empties VCaP's gDNA set); it costs efficiency when one shared ρ is clean (SD 0.0016 →
  0.0058 on LBX0588's depths at ρ = 0.05); it reads ≤ 0.003 on odg05 ss 0.99, where od 0 cost +19.8 %.
- min(x_g, x_r): 0.507 on VCaP at full depth, where both witnesses are inflated; the binomial value on odg05.
- The largest root: with the junctions alone it returns 0.2 on the 6 of 23 real runs whose least root is 0, and 0.2 is
  the most harmful value measured where the truth is 0.
- Maximum Q: Q is not a quasi-log-likelihood, and ΔQ on the two-root rows (−2.06, −2.60) sits inside its ±36 spread
  across draws.
- The walk from 0: its first step is the refused pair-count moment (one deep seed gives 0.83 against a fit of 0.017),
  and it overshoots in 4 of 1,200 low-count draws (0.2 where the least root is 0.031).
- The prototype's bisection is not a rule: it takes the upper root only because the unstable root (0.0063) lies below
  its midpoint.
- A posterior-mean od needs Q = ∫U to be a log-likelihood, which the data-dependent a_s rules out.
- Midpoint values for the fallback: ρ_eff at R = 0.2 is ×2.56, and at R = 1 (a uniform prior on the ICC) ×3.71,
  against ×1.07 for 0 on the minimal ψ.
- The sign-of-Ψ empty-evidence predicate moved 22 of 61 joint and 16 of 61 gDNA-only fits by 2–8 ulps, where the
  shipped-predicate form moves none.
PURITY WEIGHTS:
- A weight from the seed's own split ((2·min(A,T)/n)²): −66 % at ρ = 0.05; sample-splitting each seed, −28 %.
- The capture-blind own-count purity: truly pure intron seeds at odg05 g05 ON read π̂² 0.10 against 0.98.
- Transport with a Jeffreys ½: the ½ was an underived prior, and it decided every real-library row.
fl:
- Bisection to fl's fixed point: it finds a root, and on an empty pool its bracket is [0, 1000]; one refresh lands
  within 0.011 bp of the fixed point.
Sources: `~/Downloads/rigel_runs/prototypes/2026-09-28_od_design/` (`02_roots.md` §3, `03_bound.md` §3,
`04_evidence.md` §2–§4, `01_fallback.md` §2 and §5, `06_fl.md` §3).
