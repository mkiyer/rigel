# Output names: the convention, after the owner's rulings

Revision 2, 2026-10-17. The owner ruled on revision 1's four questions on 2026-10-17: word order `count_<what>`;
annotated against nascent as the two RNA categories, with single-exon transcripts annotated; "annotated
transcripts" and "nascent RNA" as the words a user reads; the index-side names renamed, with every index rebuilt.
The binding vocabulary is `DESIGN.md` §0 (annotated / nascent, in outputs only). This document is the convention
that implements it and the stage that applies it.

## The convention

1. **Two RNA categories, everywhere a user reads.**
   - **annotated**: RNA on every reference transcript, single-exon included. A single-exon transcript does not
     splice, so its nascent and mature sequence are one, and it is reported with the annotation (owner).
   - **nascent**: RNA on the synthetic spans Rigel adds to model the unspliced RNA of multi-exon transcripts.
   - beside them, **gdna**, **rna** (all RNA), and **intergenic** as a place.
   - `mature`, `mrna` and `nrna` leave the outputs. Inside the solver a span stays a synthetic transcript (AXIOM 0);
     only what a user reads calls its RNA nascent. The simulator keeps `mrna` / `nrna` for the templates it draws.
2. **A structural fact is named as one:** `is_single_exon`, never `is_nrna`.
3. **Units in the name, one spelling each:** `count` is fragments (fractional where the EM splits one); `fraction`
   a share; `effective_length` spelled out; `abundance` molar only. No `mass`, `rate`, `frac`, `_total`, or a
   bare population name for a count.
4. **Word order:** `count_<what>` and `fraction_<what>`, as the transcript table already reads (`count_unambig`).

## What changes

| output | today → the convention |
|---|---|
| `quant` | `is_nrna` → `is_single_exon`; `nrna_id` → `nascent_id`, the transcript's synthetic span (`.` when none, including when an annotated single-exon transcript covers the span); `nrna_parent_count` → `count_nascent`, that span's count (0 when none) |
| `gene_quant` | `mature_count`, `nascent_count` deleted: a gene's RNA is its annotated `count`, and the nascent spans are gene-neutral |
| `nrna_quant` → `nascent_quant` | synthetic spans only (an annotated single-exon coverer is a row of `quant`); `nrna_id` → `nascent_id`; `is_synthetic` deleted (every row is); `n_mrna`, `mrna_count` → `n_annotated`, `count_annotated` (the transcripts the span models); `tpm` → `tpm_total_rna` |
| `loci` | `mrna`, `nrna`, `gdna`, `total` → `count_annotated`, `count_nascent`, `count_gdna`, `count`; `n_transcripts` → annotated only, `n_nascent` the spans, `n_annotated_transcripts` and `n_nrna_entities` deleted; `gdna_rate` → `fraction_gdna`; `gdna_eff_len_em` → `gdna_em_effective_length`; `enable_gdna` 0/1 → `has_gdna` bool |
| `summary.json` → schema 4 | `quantification`: `count_annotated`, `count_nascent`, `count_gdna` (EM + unspliced intergenic: the numerator of `fraction_gdna`), `count_gdna_em`, `count_intergenic_rna`, and the matching `fraction_*`; `n_transcripts` annotated only, plus `n_nascent`; `intergenic_total` deleted (it is `fragment_stats.intergenic`); `strand_model.sj.n_sj_depth_ge_100` / `_ge_1000` deleted (the fit they served is gone: constants ruling 8); the rounding digits named beside `SUMMARY_DECIMALS` (ruling 9) |
| `calibration_track` | `gdna_mass`, `rna_mass`, `gdna_frac` → `count_gdna`, `count_rna`, `fraction_gdna` |
| `fragment_lengths` | categories `pool_dna_*` → `pool_gdna_*` |
| `locus_stats` | `n_units` → `n_em_fragments`; `gdna_eff_len` → `gdna_em_effective_length` |
| annotated BAM `ZF` | the winner-kind bits name the categories: 0x02 `is_annotated` (any reference transcript, single-exon included), 0x04 `is_gdna`, 0x10 `is_nascent` (a synthetic span); 0x08 is retired, so no old value is reused with a new meaning. A single-exon winner goes from 0x09 to 0x03, a synthetic one from 0x19 to 0x11; the other values stand |
| report | "Annotated transcripts", "Nascent RNA", "gDNA (EM)", "Intergenic RNA", "Intergenic gDNA"; the gene table loses its mature / nascent columns; "Transcripts" counts annotated transcripts; the strandedness tile reads `strand_live` (ruling 7) |
| the index | `--nrna-tolerance`, `IndexConfig.nrna_merge_tolerance` and the manifest's `nrna_tolerance` → one name, `nascent_tolerance`; `INDEX_FORMAT_VERSION` bumped and every index rebuilt, which also deletes the contiguous boundaries' unread reaches |
| the scan | `splicing_anchor_tolerance` / `--splicing-anchor-tolerance` (K) → `intron_overhang_tolerance` / `--intron-overhang-tolerance`: how many bases an alignment may run into an intron and an intron overhang a mate gap (ruling 6) |
| simulator truth | per-condition table: the `*_abundance` and `post_capture_*` copies deleted, `observed_*_fragments` kept, `pre_capture_*` → `molar_*`; the suite-level molar file keeps `mrna_abundance` / `nrna_abundance`, the simulator's input names; the manifest and `truth_summary.json` keep one key per path |

Four columns change VALUE as well as name, because the convention moves single-exon transcripts into annotated:
`quant.count_nascent` and `nascent_id` for a transcript whose span is an annotated single-exon transcript (it has no
nascent span), `nascent_quant`'s rows, and the `ZF` bits. Every count of the tool stays: these are relabellings of
which category a fragment is reported under.

## The stage

In order, each gated (suite on a clean worktree, the identity check on three references, preflight):

1. `rename_identity.py` digests the transcript table WITHOUT its header (today it hashes the header, against its own
   docstring), and the references are re-frozen at the commit before the stage.
2. The internal renames, which change no output: `t_is_nrna` (it holds `is_synthetic`) → `t_is_synthetic`, and the other
   review §5.6 items.
   Bit-identical.
3. The outputs to the convention, `summary.json` schema 4, the report, the MANUAL and the goldens. The diff of every
   golden is read before `--update-golden`, and the four value changes above are checked against their definition.
4. The index and scan names, `INDEX_FORMAT_VERSION` bumped, every index rebuilt (with the contiguous reaches gone)
   and every panel re-cached. The index's `is_nrna` column becomes `is_single_exon`, and with it the descriptions that
   call it "covers a multi-exon span" (`transcript.py`, `scoring.py`, `calibration/splice_graph.py`,
   `calibration/sj_opportunity.py`); the drain's `chose_genomic` / `chose_spliced`, held back from step 2 because the
   cache's schema digest covers them, become `chose_no_implied_intron` / `chose_implied_intron`.
5. The simulator's truth tables, through one writer, and the instruments that read them (`quant_accuracy.py`'s
   fallback for older tables deleted: every panel on disk carries `observed_*`). Held for the next re-simulation,
   because every panel on disk carries the old names and its readers must keep reading it: `pre_capture_*` →
   `molar_*`, the synthetic spans' `RIGEL_NRNA_` id prefix, the manifest key `post_capture_abundances` and the file
   name `truth_abundances.tsv`.
