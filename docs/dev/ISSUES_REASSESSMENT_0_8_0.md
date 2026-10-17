# ISSUES reassessment for the 0.8.0 cut (2026-10-14, sandbox, not authoritative)

All 63 open entries in `docs/ISSUES.md` were checked against `main` at 55fe212e: the count foundation, the mode-free
reader, the multimapper pass and the EM scheduling fix. The check was read-only: code, git history, and archived run
outputs. The headline findings are verified by hand at the cited lines.

Almost every panel number quoted inside these entries predates the foundation (3e7e0fbc). Re-derive a number before
ranking on it.

---

## 1. Decisions only the owner can make

**Decision 1. The reader reads calibration's false gDNA as capture on uncaptured libraries (new; not in ISSUES).** On a library
with almost no gDNA, the bulk sits at the landscape's floor, so a few false-gDNA fragments on an object read as
enormous enrichment.

Measured on the real VCaP transcriptome-only library (no gDNA spike, no capture), the shipped reader against the
foundation's:

| | Foundation | Shipped |
|---|---|---|
| Objects read above e× the bulk | 0 | 1,824 (1,786 above 1,000×; max 1.04e8×) |
| Saturation level the log prints, which `summary.json` publishes | — | **2.9e7×** |
| Transcript counts moved | — | **3.56 % of the library** (Σ\|Δ\| 478k of 13.4 M), inside 376 loci |
| Gene counts moved | — | 0.13 % |
| EM lengths of 2,730 annotated transcripts | — | inflated 1e5–1.6e7× |

In simulation the same family costs +2.3 to +2.8 transcript points on the zero-gDNA capture-OFF robustness rows, which
the owner accepted on 2026-10-13. The real-library magnitude was not in the packet, which reported only "1,824
objects > e×".

There is no truth for that library, but an efficiency of 1e8 on an uncaptured library is unphysical. Ordinary
gDNA-poor RNA-seq is the commonest input.

The likely upstream cause is already an open entry. ψ adds the fitted gDNA landscape on top of both Jeffreys halves
(`psi_kernel.h:129,181`; `ISSUES: the-gdna-prior-enters-psi-twice`). The only A/B of its fix ("Part 1"), run before
the foundation, cut the zero-gDNA false gDNA 366 → 13 fragments, but cost the capture-ON strata (stranded ON 6.14 →
8.98 %).

- *Recommendation:* one pinned A/B of Part 1, ported to the current tree, on the ladder, both gap panels, the four
  zero-gDNA capture-OFF robustness rows and VCaP transcriptome-only. That is M, about a day. If it removes the false
  capture without the old capture-ON cost, it goes in 0.8.0. Otherwise ship with the limitation documented, and do
  the repair in 0.8.1.
- *Alternatives:* ship now and document; or hold the tag for a reader-side guard. No candidate guard exists, and a
  detector or gate is refused.

**Decision 2. Whether the two long campaigns gate the release.** ROADMAP still lists, as release-critical, Tier 0's strand
pruning and the splicing-artifacts cluster track. Both are L and need cluster infrastructure.

- *Recommendation:* move both after 0.8.0, and keep only their cheap pieces in the cut: the `summary.json` strand
  diagnostics, the blacklist count column in the index bump, and the blacklist warning.

**Decision 3. The census length law.** CLAUDE.md, DESIGN §0b and ROADMAP:62 put fragment length in 0.8.0, while ROADMAP:212
defers `ISSUES: the-scorer-reads-a-census-length-law`. It has no candidate: its prize was measured on the deleted
ruler, and the real length gaps lose.

- *Recommendation:* 0.8.1, with the docs reconciled.

**Decision 4. Overlapping synthetic shadows.** `nrna_quant.feather` has one row per synthetic span, and 6,366 of 6,919 spans
overlap a same-strand twin. `quant.feather`'s `nrna_parent_count` can read about 0 while a sibling twin holds the
family's mass. The MANUAL lists seven `nrna_quant` columns the code never writes, and no test calls
`get_nrna_counts_df`.

- *Recommendation for 0.8.0:* fix the docs, add a column gate, and document that nascent counts are not identifiable
  among overlapping spans.
- *After 0.8.0:* the index merge, which is L. It changes the EM's component set and orphans the panels' nascent truth
  IDs.

**Decision 5. The message layer misses its own bar on one in-scope stratum (new).** On the ladder, transfer now loses to silent
on all three contaminated unstranded × capture-OFF rows: 1.02×, 1.07× and 1.31× (g05, g50, g98;
`policy_benchmark.py --panel ladder`). It still wins unstranded × ON by about 4× and wins every contaminated stranded
row. This is most likely the foundation's rebuilt strand claims, and DESIGN still says every contaminated row favours
transfer.

- *Recommendation:* record it now, and investigate in 0.8.1 beside
  `ISSUES: the-intron-own-solve-on-unstranded-capture-off`. Its transcript-level cost is unmeasured.

---

## 2. Proposed 0.8.0 cut

Each item lands in its own commit, with a falsification test verified failing first. "Moves numbers" marks the items
that need a pinned A/B.

**A. Silent failures and wrong output on user paths (S each).**
1. `--overhang-alpha 0`, which the CLI and the MANUAL document as a hard gate, gives 0 × −∞ = NaN for every
   zero-overhang candidate (`scoring.cpp:553,815`). The prune tests then drop every EM-routed fragment silently. The
   mismatch term is already guarded (`:492,766`).
2. A failed later refit overwrites the landscape the final solve ran under with `None` (`calibrate.py:368-379`). The
   capture reader is skipped and the result says no landscape was fitted. The fix keeps the last good landscape, and
   is a numeric no-op wherever no refit fails.
3. **Moves numbers (M).** The EM prior's gDNA share divides a count that includes the multimappers drawn since the port
   by `n_calibrated`, which excludes every multimapper unit (`pipeline.py:709-719`). The docstring and
   `tests/test_em_pseudocounts.py` still say "calibration deposits no multimapper". Validate on VCaP truth and the
   aligned multimapper scenarios.
4. `index.py:487` writes `nrna_n_contributors =` where it should be `+=`. It feeds the user column
   `n_contributing_transcripts`, and rides the index bump (B).
5. GTF parsing:
   - a final attribute without its `;` is dropped (`gtf.py:42`);
   - a missing `transcript_id` raises outside the warn-and-skip path;
   - an exon on another chromosome or strand merges silently.
6. The splice-junction strand tag:
   - validate the spec (an unknown one silently becomes XS_TS);
   - log the resolved mode and write it to `summary.json` and `config.yaml`, which today record "auto";
   - warn loudly when "auto" resolves to "none" on a library with spliced reads.

   A BAM with no XS/ts tag loses its strand entirely. The repair, taking the annotated junction's strand, is M and
   goes in 0.8.1.
7. A dropped splice blacklist is silent in both directions (`index.py:1503-1521`), and the tests pin the silence.
8. An unknown `config.yaml` key errors with no hint, and every 0.7.1 config carries four removed keys.
9. The report's meta line is written through `innerHTML` (`report.js:301`).
10. `rigel sim` hangs forever on an empty fragment-length window (`sim/sampling.py:27-37`).

**B. The format batch, in one window (M; every benchmark and cluster index rebuilt once).**
- `summary.json`:
  - one od field, labelled binomial by policy;
  - κ pooled and clean, the artifact and reversed shares, and a strand-live flag;
  - the resolved SJ tag;
  - `gdna_fraction` without the spliced intergenic fragments (`n_intergenic_spliced` is already counted);
  - the report's capture note for older summaries.
- Index format 9:
  - drop `abundance`, `nrna_abundance`, `n_exons` and `sj.feather`'s `interval_type`;
  - add the duplicate-transcript alias map, fix 4, and the blacklist count column;
  - update the CHANGELOG's "version 8".
- The annotated BAM labels spliced intergenic fragments as gDNA (ZF 0x25), the same Axiom-0 slip as `gdna_fraction`.

**C. Two numeric no-ops (S each; `rename_identity.py`).**
- `get_loci_df` builds a mask per locus (`estimator.py:821-823`).
- The buffer holds one chunk too many (`buffer.py:328-350,506-515`).

**D. The release docs, made true.**
- CHANGELOG:24 advertises a strand-overdispersion estimator, but od has been 0 since 2026-10-01.
- README:46 advertises a `map` assignment mode the CLI rejects.
- The MANUAL says STAR writes XS, which it does only with `--outSAMstrandField intronMotif`.
- MANUAL:296 says `pip install alignable`, while `mamba_env.yaml` says it is not on PyPI.
- Decide the `--no-mappability` rename now, since it is a public flag.
- Document `sj_blacklisted`'s unit.
- Fix the report subtitle "Capture on-target enrichment".
- DESIGN §7.2's later paragraphs still give the old clipped efficiency and "never above 1" (my miss on 10-14).
- ISSUES' grid-point entry wrongly calls a bulk step a common factor: bulk-level objects read 1 wherever the bulk
  sits, while enriched objects move with it (my error).
- `simplex_logodds.py:9-11` contradicts `:89-91`.
- DESIGN's message-policy sentence (decision 5).
- A MANUAL known-limitation line for decision 1, unless the A/B fixes it.

**E. The list itself.** Apply the closures (§3), condense the two campaign logs, and rewrite ROADMAP's "Next — the
order", which predates the reboot.

**Then the release checklist** (ROADMAP item 7):
- the deliverable measured per stratum;
- the suite at its standing count;
- `preflight.py --full`;
- CI run once by hand;
- a real-data smoke run;
- the standing risks re-read.

---

## 3. Close, merge or condense

| Entry | Action | Why |
|---|---|---|
| calibration-detects-capture-on-a-capture-off-library | close the headline; refusals to REFUSED; the live remainder becomes decision 1's new entry | RNA-short stranded OFF 42.30 → 3.55 % (3e7e0fbc's face maps plus e74da234) |
| the-background-dispersion-assumes-a-pure-intergenic-pool | close | `density_deconv` and the intron factory are deleted (3e7e0fbc); its gDNA-only object row moves to `instrument-ledger` |
| gdna-landscape-trains-on-false-positives | merge into the collapse entry and the ψ-twice entry | (c) is the ψ double prior; (d) is closed; (a) is negligible |
| multimapper-blind-support | close | fixed by the multimapper pass |
| strand-likelihood-over-confident-beyond-od | close | its trigger never landed; od = 0 by policy |
| em-gdna-exceeds-calibration-on-the-vcap-transcriptome-half | close | no truth; its od arms are moot; it predates the port |
| splice-out-premise-bias-uncorrected | close | DESIGN §6b.4 already rules it recorded and not corrected |
| the-tilt-census-as-an-instrument | close | never built; its only consumer is an accepted limit |
| flux-floor-dispersion, transfer-variance-premise | merge into `message-layer-open-cases`; transfer-variance takes the TRANSPORT and SPLICE_OUT width terms | no instrument since 2026-09-14 |
| crossing-pool-contrast | merge into `capture-blind-gdna-divisor` | the census law already ships a crossing contrast |
| message-layer-open-cases (b), instrument-ledger (c), binary-cuts ρ_ref item | close the sub-items | the factory, the pickles and ρ_ref are deleted |
| the-gdna-prior-enters-psi-twice | condense 2,054 → about 12 lines | about half of it is fixed in 3e7e0fbc; its density-reader half is superseded by DESIGN §7.2 |
| the-edge-density-floor-under-capture, the-scorer-reads-a-census-length-law | condense to about 15 and 25 lines; drop "now" | §4 |
| pure-rna-mirror-asymmetry | re-measure once and close, or close as stale | August numbers on a rebuilt calibration |

New entries:
- decision 1, the false capture on uncaptured libraries;
- fix 3, the multimapper denominator;
- decision 5, the unstranded × OFF policy bar;
- the deferred RNA-long regression after the factory removal (6.55 → 43.22 → 43.34), which the owner authorized
  investigating for robustness.

---

## 4. Every open entry: now → proposed

Status codes: **open** (as described); **partly**; **fixed**; **stale** (its code was deleted or replaced);
**superseded**.

| Entry | Tier now | Proposed | Status | Size | User-facing | Note |
|---|---|---|---|---|---|---|
| strand-overdispersion-one-shared-value | now | diagnostics 0.8.0 (B); pruning later | partly | S–M / L | diagnostics | od = 0 and κ from genuine junctions have landed |
| splicing-artifacts | now | campaign after 0.8.0 (decision 2); cheap pieces 0.8.0 | open | L | yes, real data | reject rule, builder threshold and count column unchanged |
| the-format-changes-to-batch-before-release | now | 0.8.0 (B) | 1 of 7 done | M | yes | 0.7.1 shipped index 7, so users rebuild anyway |
| an-unrecorded-splice-blacklist-is-dropped-silently | now | 0.8.0 (fix 7) | open | S | yes | |
| calibration-detects-capture-on-a-capture-off-library | now | close; decision 1 | partly | — | — | §3 |
| the-edge-density-floor-under-capture | now | later | open | L | maybe, under capture | its withdrawal was not selected |
| the-scorer-reads-a-census-length-law | now | 0.8.1 (decision 3) | open | L | under capture | |
| the-gdna-prior-enters-psi-twice | next | 0.8.1, or 0.8.0 via decision 1's A/B; condense now | partly; core open | S + M | indirect (decision 1) | gate 1 re-run on this tree: 0.85 nats |
| the-gdna-landscape-collapses-at-low-depth | next | refit bug 0.8.0 (fix 2); re-measure, then 0.8.1 | partly | S / M / L | yes, shallow | its depth numbers predate the foundation |
| gdna-landscape-trains-on-false-positives | next | merge | superseded | — | no | |
| the-background-dispersion-assumes-a-pure-intergenic-pool | next | close | stale | — | no | |
| capture-on-overcalls-gdna-at-low-gdna | next | later | open, stale numbers | L | modest | the owner rules over-calling gDNA preferable (10-09) |
| strand-plug-in-bias-on-sparse-libraries | next | 0.8.1 with the low-depth work | open | L | yes, shallow | ψ still reads out a median |
| latent-defects | next | 7 items 0.8.0 (fixes 4–10); rest 0.8.1 or later | 23 of 25 present | S each | 7 on user paths | the zero alpha is live, not latent |
| sj-strand-tag-chosen-from-the-first-reads | next | 0.8.0 minimum (fix 6); repair 0.8.1 | open | S + M | yes | wider than stated: a BAM with no tag loses its strand |
| hygiene-ledger | next | release-gating docs 0.8.0 (D); rest later | mostly open | S–M | docs | the reorder gate is closed, but its replacement has no gDNA |
| performance-memory-bounded-solve | next | 2 no-ops 0.8.0 (C); rest later | open | S | writer time | ①–⑤ need a re-time before ranking |
| the-em-stops-before-it-converges | next | 0.8.1 | open | L | indirect | written 10-14 |
| nascent-stress-sensitivity | later | later | open, stale | M | no | |
| psi-reads-kappa-where-the-strand-channel-is-dead | later | 0.8.1 | open | S + M A/B | minor | ψ builds κ̂ unconditionally (`solve_kernel.cpp:701`) |
| the-fl-boundary-inversion-reads-missing-evidence-as-a-value | later | later; (e) 0.8.1 | open | L / S–M | indirect | |
| the-fl-boundary-inversion-has-underived-pieces | later | later; (e) 0.8.1 | open | M–L / S | indirect | (e) deletes a tunable |
| capture-blind-gdna-divisor | later | later | open | L | indirect | |
| eb-shrinkage-magic-ess | later | later | open | M–L | indirect | `pool_prior_ess = 1000`, not derived |
| the-capture-reference-is-read-at-a-grid-point | later | later; correct the text now | open | L | indirect | §2 D |
| intron-seeds-near-probes-are-capture-enriched | later | later | partly moot | M | indirect | |
| strand-likelihood-over-confident-beyond-od | later | close | superseded | — | no | |
| the-pseudocount-strength-is-not-derived | later | 0.8.1, with the EM work | open | M | indirect | |
| per-transcript-prior-lane | later | later | open by ruling | L | no | |
| message-layer-open-cases | later | later; close (b) | partly moot | M | no | |
| refit-vs-message-arbitration | later | later; merge with the next row | open | M–L | indirect | |
| the-intron-own-solve-on-unstranded-capture-off | later | 0.8.1, with decision 5 | open, stale | M | indirect | |
| the-efficiency-posterior-floor-on-empty-pieces | later | 0.8.1, with the sparse losses | new form | M–L | indirect | a scant object reads the prior; the read is discontinuous as S → 0⁺ |
| scoring-penalties-are-underived-constants | later | NaN 0.8.0 (fix 1); constants later | partly | S / L | yes (NaN) | |
| multimapper-intergenic-alignments | later | later | EM side open | L | yes, real data | |
| multimapper-blind-support | later | close | fixed | — | — | |
| unwitnessed-loci-and-multimappers-at-the-em | later | denominator 0.8.0 (fix 3); rest later | partly | S + M | yes | |
| vcap-dna-only-objects-lose-gdna-at-full-depth | later | later; re-measure | open, stale | M / L | internal | |
| em-gdna-exceeds-calibration-on-the-vcap-transcriptome-half | later | close | stale | — | no | |
| instrument-ledger | later | later; (e) before any spectrum expansion; close (c) | 5 of 7 open | S each | no | gains the gDNA-only object row |
| debug-capture-memory-is-unbounded | later | later | open | S–M | no | |
| expand-the-gdna-spectrum | later | later | open | L | no | |
| a-pure-gdna-library-reads-as-nascent-rna | parked | parked; a MANUAL limitation line | open, stale | L | deferred stratum | |
| the-pooled-q-in-the-gdna-count | parked | 0.8.1, with the census law | open | M | mildly | |
| unannotated-transcription-is-booked-as-gdna | parked | its slice 0.8.0 (B); design later | open by design | S / L | yes (`gdna_fraction`) | |
| ruler-witness-geometry-on-transcript-panels | parked | later | partly | L | indirect | the junction bound covers part of it |
| overlapping-synthetic-shadows | parked | 0.8.0 docs and gate (decision 4); merge later | open | S / L | **yes** | |
| yield-variance-beside-the-count | parked | later | not built | M | a feature | |
| capture-premise-untested-on-cdna | parked | parked; optional MANUAL line | risk | — | potentially | |
| the-atom-at-an-unwitnessed-both-strand-slot | parked | parked | open, stale | L | narrow | |
| two-sided-exon-row | parked | parked | open, stale | L | indirect | |
| the-lower-bound-noise-ratchet | parked | parked | open, stale | L | indirect | |
| flux-floor-dispersion | parked | merge | no instrument | — | no | |
| splice-out-premise-bias-uncorrected | parked | close | superseded | — | small | |
| the-tilt-census-as-an-instrument | parked | close | never built | — | no | |
| transfer-variance-premise | parked | parked; gains the width terms | open | M–L | no | |
| drain-contaminates-certified-rna | parked | later; the QC field 0.8.1 | open | S / M | no | |
| crossing-pool-contrast | parked | merge | partly moot | — | no | |
| capture-degeneracy-standing-risk | parked | parked (watch) | risk | — | not today | |
| pure-rna-mirror-asymmetry | parked | re-measure or close | stale | S | possibly | |
| binary-cuts-on-continuous-quantities | parked | parked; ρ_ref item closed | partly | M–L | mixed | `max_frag_length` = 1000 and TAIL_DECAY are on user paths |
| the-sparse-capture-losses, the-bulk-can-be-calibrations-nascent-rna | later | unchanged | open | — | — | written 10-13 |
