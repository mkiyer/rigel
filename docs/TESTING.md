# TESTING — the panels, the gates, and what the suite can judge

The manual for Rigel's test substrates: the panels (§0), the test chromosome (§0a), the length-truth and STAR panels
(§0c, §0d), the simulator (§1), the build commands (§2), the test suite (§5), and profiling and identity checks
(§6). Not here: the scope and its rulings (`DESIGN.md` §0b), how performance is judged and the
instruments' run order (`SUCCESS.md`), the standing pass count (`CLAUDE.md`), lessons (`TRAPS.md`).

Two rules govern the file: when a simulated axis is the axis being judged, gate the simulator on it
(`TRAPS: prove-the-substrate`); and prove a suite resolves the axis being changed before quoting a number from it
(`TRAPS: can-the-benchmark-resolve-it`).

## 0. The panels — one ladder and its side panels

`$RIGEL_RUNS` is the benchmark data root, outside the repository. Everything simulated lives under `$RIGEL_RUNS/suite/`
(the ladder, the fl-gap panels, the carved reference and its index) and `$RIGEL_RUNS/test_reference/` (§0a). No panel
keeps FASTQs: every scorer reads the simulator's BAM, and `panel.py align` rebuilds the reads from it.

**The ladder** (`scripts/sim/configs/gdna_ladder.yaml`) is the only panel the tool is ranked on.

| | |
|---|---|
| conditions | 16: gDNA `g00 / g05 / g50 / g98` × strand specificity `0.50 / 0.99` × capture off / on |
| depth | 10 M fragments per condition, fixed; the gDNA rate decides only the split |
| fragment lengths | identical for both origins (206 ± 98, truncated to [50, 500], reads of 100) — a forcing function, `DESIGN.md` §0b |
| sequence | known bases only (no fragment covers an `N` run); substitution errors at 0.002 per base, each oracle record's NM their count |
| nascent RNA | `mode: sparse`, `on_fraction 0.50` of spans on, level logU(1, 100) independent of the annotated level — a stress level, not real data |

| rung | g00 | g05 | g50 | g98 |
|---|---|---|---|---|
| rate (`f_gdna = rate/(1+rate)`) | 0.0 | 0.052632 | 1.0 | 49.0 |
| f_gdna | 0 | 0.05 | 0.50 | 0.98 |
| n_rna (the total is held) | 10.0 M | 9.50 M | 5.00 M | 0.20 M |

`g00` is the zero control on every stratum, never tuned on; `g05` is the 1–10 % corner real libraries live in. Three
levels is a floor (`TRAPS: a-single-level-panel-cannot-see-a-constant`). Report per stratum, never pooled.

**The fl-gap side panels** are never a ladder rung: the ladder equalises lengths, so it is a null for any mechanism
that reads a gDNA-vs-RNA length difference. The two arms have opposite sign (RNA is not reliably longer than gDNA), so
a real repair moves them in opposite directions; they share everything else with the ladder.

| | `suite/flgap_rna_long/` | `suite/flgap_rna_short/` |
|---|---|---|
| config | `scripts/sim/configs/flgap_rna_long.yaml` | `scripts/sim/configs/flgap_rna_short.yaml` |
| gDNA / RNA fl configured | 75 ± 20 / 250 ± 60 | 250 ± 60 / 75 ± 20 |
| conditions | `g50` × ss {0.50, 0.99} × capture {off, on} | the same 4 |

Read the realised gap count-weighted off each condition's `truth_fragment_lengths.tsv`, never the config
(`TRAPS: configured-lengths-are-not-realised`). A calibration mechanism is ranked here only on what stops before the
EM, because a length gap hands the EM part of the answer (`TRAPS: a-length-gap-bypasses-calibration`).

**The test chromosome's fl arms** (`scripts/sim/configs/test_reference_fl_{rna_long,gdna_long,equal200}.yaml`, the 30
conditions of §0a each) set gDNA / RNA to 100 ± 50 / 250 ± 150, 250 ± 150 / 100 ± 50, and 200 ± 100 for both. The
equal arm is the control: an fl-gap result that also moves there is an artefact.

### Every scenario is cached

Nothing in the loop *change calibration → re-run calibration → score against oracle calibration* re-scans a BAM.
`panel.py cache` builds both caches, and neither is invalidated by a calibration change:

| cache | holds | invalidated by |
|---|---|---|
| scan (`build_scan_cache.py`) | the accumulator payload | an accumulator change (`payload_schema_digest`), the index (`graph_hash`, `reach_digest`), a tally setting of the scan (`sj_strand_tag` among them, never a thread count) |
| oracle (`<panel>/oracle_cache/`) | five origin partitions per condition (`gdna`, `mrna`, `nrna`, `rna_pos`, `rna_neg`), checked sum-to-full against the scan cache, and the certified `slot_truth.npz` | the same, and nothing else |

`panel.py status` names what is missing, and an instrument fed a stale cache refuses it. ⛔ The scan cache's key
hashes the deposit rule, not the resolver (`native/resolve/`): after a change to which fragments are offered, rebuild with `--force`.

## 0a. The test chromosome — where the message policy is developed

A synthetic chromosome designed one structure at a time in one YAML file, swept over 30 conditions and scored in
seconds, so a policy change is scored on every condition between edits; the ladder makes the shipping judgement.

| axis | values | where it lives |
|---|---|---|
| gDNA fraction | `g00` `g05` `g25` `g50` `g98` | `scripts/sim/configs/test_reference.yaml` |
| strand specificity | `0.50` · `0.70` (the transition) · `0.99` | the config |
| capture | off · on; about half the genes probed, half with no probe | the config, `test_probes.bed` |
| nascent RNA | ten transcripts (the twin block's `nasc`, `capnasc`), hand-authored below, near and above the annotated level | the `nrna_abundance` column |

Depth is 1.17 M fragments per condition. Five more configs differ from `test_reference.yaml` in one thing each: two
adversarial probe panels (`test_reference_probes_{junction,sparse}.yaml`) and the three fl arms (§0); three more are the
gDNA-depth ladder (§0c).

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
T=$RIGEL_RUNS/test_reference
# 1. render the GTFs, abundances and probe BEDs (repo and runs dir) and the FASTA; `--check` only compares
python scripts/sim/build_test_reference.py
rigel index --fasta $T/test_chr.fa --gtf $T/test_chr.gtf --collapse-duplicate-transcripts --no-mappability --no-tsv -o $T/idx
# 2. per panel. ⛔ one panel's `cache` at a time, or a RIGEL_SCRATCH each: the panels share condition names
CFG=scripts/sim/configs/test_reference.yaml
python scripts/sim/panel.py simulate --config $CFG
python scripts/sim/panel.py cache    --config $CFG     # both caches, then certifies slot_truth.npz
# 3. score every policy on every condition — seconds
python scripts/design/policy_benchmark.py --panel test
```

A failed FIELD gate at certification is reported, not fatal: the row is then COMPOSITION-certified. ⛔ **Rebuild after
every edit to the YAML**: the FASTA, index, reads, caches and `slot_truth.npz` all describe the annotation they were
built from, and a benchmark on stale caches answers a different question.

**How to read the benchmark — two halves, two bars.** `silent` is the measured floor. On strand-specific data a sighted
exon's own strand solve is excellent and a message can mostly only disturb it; message passing exists for unstranded
data, where the strand channel is dead and the local answer is a default.

| rows | the bar |
|---|---|
| unstranded (`ss 0.50`) | the policy must win — this is what the message layer is for |
| stranded (`ss 0.99`) | the policy must do as little harm as possible — near 1.00× of silence is a pass |

`ss 0.70` rows count against neither bar. Never pool the halves, nor conditions within them; confirm on the ladder,
since a toy and the panel have inverted a ranking before (`TRAPS: panel-before-src`).

### The substrate — one hand-edited file, thirteen blocks

`scripts/sim/test_reference/test_chr.yaml` is the `rigel sim` scenario schema
(`genes → {gene_id, strand, transcripts: [{t_id, exons, abundance, nrna_abundance}]}`, exons 0-based half-open) plus
`probed` per gene and a `shadow_genes` list. Its header is the design record of every block.

| block | what it stresses | gene types |
|---|---|---|
| anchored twin | the message layer's five controls on one shape (exon 1 kb · intron 7 kb · exon 1 kb · intron 7 kb · exon 1 kb), each varying one bit — probed, nascent, silent — × five abundance blocks | `clean` `nasc` `cap` `capnasc` `silent` |
| mono | single-exon: no junction, the gene edge the only message source | `mono` `capmono` + two silent controls |
| isoform | a host and a second isoform, replicated across A ≫ B, A ≪ B, A ≈ B | `altstart` `altss` |
| walled | exon pieces with no licensed face | `chain` `tssalt` `tandem` `altlast` |
| terminus-cluster | ten transcript ends 126–147 bp into a shared last exon (after MIR99AHG) | `cluster` `capcluster` |
| both-stranded | a host and its antisense in one locus (after TTC28-AS1, PPM1F-AS1, a convergent pair) — the AMBIG substrate | `asin` `asinrev` `span` `conv` |
| sj+terminus | one boundary carrying a junction and a same-strand terminus (after RUNX1, LARGE1) | `sjterm` `capsjterm` |
| encompassing | a 20 kb single-exon antisense over a whole twin-shape host (after ENSG00000280007 over TUBA8) | `enc` `capenc` |
| in-exon | a 500 bp single-exon antisense inside the host's last exon (after ENSG00000273300 in UFD1's 3' UTR) | `inexon` `capinexon` |
| shared-exon | two spliced genes sharing one exact 10 kb last exon (after SMARCB1 × DERL3): one walled AMBIG exon | `shared` `capshared` |
| in-intron | a 600 bp single-exon gene in the opposite strand's first intron: walled, its only source its own termini | `inintron` `capinintron` |
| head-to-head | `conv` with the strands swapped: two 5' ends overlap (after TRMT2A ⟷ RANBP1) | `div` `capdiv` |
| tiny-exon | exons shorter than a fragment, seen only through their boundaries | `tiny` `captiny` `mixed` `capmixed` |

273 genes on a 7.93 Mb chromosome, equal + / − representation (a sign error is invisible on one strand; the builder
refuses an imbalance), abundances in half-decade molar ladders, annotated up the blocks and nascent down. The
`shadow_genes` on the blank contig `test_blank` are transcription the simulator draws and the index never sees: the
control for anything that trusts the annotation's notion of "pure gDNA". The scan withholds the spliced reads of the
two spliced shadows (U5, U6); their unspliced reads, and the single-exon shadows', are the control.

The GTFs, `test_abundances.tsv` and the probe BEDs are rendered by `build_test_reference.py` and versioned beside the
YAML (`tests/test_test_reference_renders.py` refuses a drifted render): the benign panel tiles each probed gene's exon
union 8 × 125 bp, the sparse panel centres one probe per exon, the junction panel places a two-block probe per
junction. The FASTA is generated from a fixed seed with a GT..AG motif at every declared intron, so the aligner, the
simulator and the annotation agree; the builder reports every problem at once, and `--self-test` perturbs each check.

⚠ κ is fitted from spliced reads, so a transcriptome of only single-exon transcripts calibrates as unstranded,
silently (`tests/calibration/test_calibrate.py::test_an_empty_spliced_census_calibrates_as_unstranded`): a stranded
number needs a multi-exon transcript with real depth. `rigel index` makes the synthetic nascent spans itself.

## 0c. The length's truth instrument and the gDNA-depth ladder

`scripts/design/ruler_vs_truth.py` scores the EM's capture-contracted length
(`capture_eff_length.transcript_capture_eff_lengths`) against the simulator's own: `CaptureSampler.partition_array` drew
the reads, so the truth is `Σ_w f_pre(w) · partition_t(w)`. Ratios are compared in the log, anchored on the fully
probed transcripts, per probed class and per kind (annotated multi-exon, single-exon, synthetic span).

```bash
python scripts/design/ruler_vs_truth.py --panel test --condition gdna_g05_ss_0.99_nrna_file_capture_on
python scripts/design/ruler_vs_truth.py --panel ladder                     # one line per capture-ON row
python scripts/design/ruler_vs_truth.py --panel ladder --condition gdna_g50_ss_0.99_nrna_mid_capture_on --scale
python scripts/design/ruler_vs_truth.py --panel test --condition C --module proto.py --out table.tsv
```

`shipped` is the length in `src/`; `oracle_gdna` feeds the same code the certified true gDNA per object, so the gap is
calibration's and the rest the length rule's or the geometry's. `--module` runs prototype rulers beside the shipped
one on the same calibration, and `--set` prices a config value on every arm. Read the classes apart: the probed class
is where the formula is exact when its witness is, the unprobed where the bulk bites, the partial where the junction
price lives. `--scale` is the one-scale read-out of `SUCCESS.md`.

**The gDNA-depth ladder**: `test_reference_depth_d10.yaml` (117 k fragments) and `test_reference_depth_d100.yaml`
(11.7 k), each at `g00 g001 g01 g05 g25 g50`, and `test_reference_depth_full_lowg.yaml` (1.17 M at `g001 g01`),
stranded only, otherwise `test_reference.yaml`. They locate where the length's reference holds as gDNA fragments per
probed piece fall; `ruler_vs_truth.py --panel-dir` scores each set. Its transcript number is not a calibration result.

## 0d. The STAR rung — aligned input, where a scan change is decided

Every ladder condition is also read from `<condition>/rustar_recommended.bam`: its reads aligned by `rustar` (STAR's
Rust port) with the MANUAL's recommended settings (`scripts/sim/configs/star_recommended.params.txt`). `panel.py align`
writes it once; a condition counts as aligned when its `rustar_recommended.align.json` records the aligner, its version
and the parameters, after the BAM was checked to hold every template drawn. `panel.py score --bam rustar_recommended`
scores the base arm against the same truth, adding what the aligner left unplaced by true origin (`unaligned_*`) and
the gDNA fraction's truth over what it placed (`gdna_frac_true_aligned`).

A scan change is decided here, pinned, fractional and per stratum, because these settings report the ambiguity the
simulator's own alignments hide; the simulator's alignments stay the regression check, which a scan change that should
not touch them leaves identical. The rung holds one realisation per condition.

## 1. The simulated backbone

`chr21` + `chr22` + the 92 ERCC spike-ins, carved from the GRCh38 + GENCODE v46 sources of the production index. Two
genomic chromosomes and the spike-ins make the reference-id space non-trivial on the RNA and gDNA paths
(`TRAPS: one-reference-hides-refid-bugs`); every config states `gdna.genomic_refs: [chr21, chr22]`, since the engine
infers nothing (`TRAPS: annotated-is-not-genomic`). The lengths, 206 ± 98, are a real library's RNA pool, pre-capture;
capture selects for length (`f_post(w) ∝ f_pre(w) · total_eff(w)`), so ⛔ score against `truth_fragment_lengths.tsv`.

`rigel sim` builds (or is given, `index:`) a rigel index and simulates its transcript list, so what is simulated is
what `rigel quant` reads: the annotated transcripts plus the synthetic spans `index.create_nrna_transcripts` makes.

| | |
|---|---|
| nascent RNA is a transcript | a span's molecules are its `nrna_abundance`, on its own template, reads tagged `nrna_`. Under `sparse` a span is off with probability `1 − on_fraction`, else log-uniform over `abundance_ranges` |
| one multinomial over every RNA row | `prob ∝ abundance × capture-aware effective length`, annotated and span rows together |
| capture binds by genomic overlap | a fragment binds through ONE contiguous part of a probe, never the sum: a transcript holding a probed junction holds the probe whole, while gDNA, a span and an isoform without the junction bind the better part; gDNA and nascent RNA enrich alike |

Capture depletes nascent RNA, and the post-capture lengths and abundances are the truth (`truth_summary.json`). The
simulator does not do: dispersion beyond Poisson (it can plant a gDNA strand overdispersion, which no panel uses); R1-sense panels
(`ReadSimConfig.r1_sense` exists, gated in `tests/calibration/test_strand_sense_convention.py`, but
`run_condition_grid` in `scripts/sim/_sim_grid.py` does not expose it); gDNA on the spike-ins, so `taper_g = 1` is untested.

## 2. Building it — `panel.py`, one command per stage

`scripts/sim/panel.py` derives every path from one panel YAML; each stage is resumable and shells out to the
instrument that owns it.

```bash
export OMP_NUM_THREADS=1 RIGEL_SCRATCH=<scratch>   # else the instruments write tens of GB under /tmp
CFG=scripts/sim/configs/gdna_ladder.yaml
python scripts/sim/panel.py status   --config $CFG              # always first
python scripts/sim/panel.py build    --config $CFG              # index + capture probes
python scripts/sim/panel.py simulate --config $CFG --jobs 8
python scripts/sim/panel.py cache    --config $CFG --jobs 8     # BOTH caches — mandatory, §0
python scripts/sim/panel.py align    --config $CFG --aligner <rustar> --aligner-index <its index>
python scripts/sim/panel.py score    --config $CFG --jobs 2 --arms base oracle noop
python scripts/sim/panel.py report   --config $CFG --arms base oracle
```

`build` cannot carve the reference (it needs the source genome and GTF), so it prints the command and stops:
`build_suite_reference.py --fasta $REFS/genome_controls.fasta.bgz --gtf $REFS/genes_controls.sorted.gtf --refs chr21
chr22 --ercc -o $SUITE/reference`. ⛔ **`--force` does not reach `simulate`**: the simulator skips a condition whose
oracle BAM exists, so delete the condition directories before re-simulating. `status` counts an oracle condition
complete only with all five partitions; one row by hand is
`calibration_oracle.py --suite $SUITE/ladder --index $SUITE/rigel_index --build --condition <name>`.

## 5. The test suite

```bash
python -m pytest tests/ -q
python -m pytest tests/ --update-golden        # regenerate tests/golden/ after an intended output change
```

`CLAUDE.md` holds the standing pass count; any failure is a regression. The goldens run fractional with the EM on one
thread and a pinned seed (`tests/test_golden_output.py`), but the scan's float sums vary in the last bits at its
default thread budget (`TRAPS: the-deliverable-is-not-reproducible-by-default`): regenerate twice and diff.
`tests/scenarios_aligned/test_multimap_counting.py` holds identical paralogs to the total and the shape
(`assert_identical_paralogs`: even, or a vertex, never between) with the EM run to convergence (`CONVERGED_EM`); never
fix it by moving a seed (`TRAPS: identical-paralogs-are-bimodal`).

**The accumulator's identities**, a closed stage kept as a regression check after any accumulator or native change:

| identity | holds | gate |
|---|---|---|
| C++ against the executable specification | byte-identical | `tests/native/test_accumulator_spec.py` |
| one BAM at 1/2/4/8 workers | count banks bit-identical, float64 mass banks within n·ε | `tests/test_scan_order_independence.py` |
| `Σ region_start_count == deposited`; `deposited + deferred + dropped_* == offered` | exact | the same |
| the origin partitions sum to the full payload | exact on counts, within the representation on float64 | `scripts/design/_oracle.py` |

**A multimapping fragment is one fragment with several placements** (`DESIGN.md` §4), the rule stated once in the
specification's `Accumulator.offer` and held three ways: `tests/native/test_multimapper_offer.py` case by case,
`tests/native/test_multimapper_native_parity.py` array by array against `AccumulatorSet::offer`, and
`tests/test_second_pass_placements.py` for the placement term and the draw.

**The capture-contracted length is gated per object against the deposit rule**: every placement of every fragment
width runs through the reference accumulator (`tests/native/_accumulator_reference.py`), holding each piece's contained
share and each boundary's and junction's conserved share to 1e-12 (`test_effective_length.py` and
`test_capture_eff_length.py`, `test_each_share_is_what_the_deposit_rule_gives_the_object` and its gDNA twin) — per
object, never on the total (`TRAPS: conservation-misses-mis-attribution`), on pieces shorter than a fragment.

**The EM's two length tables** are gated on the design's argument (`DESIGN.md` §8.3) in
`tests/calibration/test_shared_length_law.py` — equal chemistry under capture is one table, a real gap in either
direction is never partly shrunk, the bound is the worst case over admissible capture responses, the weight is
continuous — with a deliberate break behind each gate; `tests/calibration/test_fl.py` holds each pmf's variance.

## 6. Profiling, timing and numeric identity

**Profile on a deep, real, high-complexity library**, never a panel condition or a toy; the cfRNA libraries on disk
are sparse smoke tests (`TRAPS: a-profile-share-is-a-ranking`). **Read a timing only from a back-to-back A/B pair**:
untouched stages drift between sittings, so alternate the arms in one sitting and check the untouched stages read
1.00. `scripts/profiling/profiler.py` runs the pipeline as a tree of named stages with per-stage peak and held RSS
(`--set`, `--scan-only`, `--compare`). A native port inside the sweep is held to the transfer gates
(`tests/calibration/_transfer_harness.py`, against `_solve_impl.testing`'s `transfer_prepare` and `transfer_pass`).

**Every speed-up, rename or refactor is a proven numeric no-op** (`scripts/design/rename_identity.py`): `--freeze` a
reference once from the tree a stage starts on, then `--check` each step against it, on a ladder condition and on a
real library (`--bam`). The default `--reference` does not exist and prints ⛔ rather than failing, so pass it and read
the log. The current references, frozen on the second pass's draw:

    P=$RIGEL_RUNS/arms/2026-10-29_identity S=$RIGEL_RUNS/suite
    python scripts/design/rename_identity.py --check --suite $S/ladder --index $S/rigel_index \
      --condition gdna_g05_ss_0.99_nrna_mid_capture_off --reference $P/identity_head_g05off.json
    python scripts/design/rename_identity.py --check --index $RIGEL_RUNS/refs/rigel_index \
      --bam $RIGEL_RUNS/cfrna/vcap_mix/bam/star.srt.rmdup.collate.bam --reference $P/identity_head_vcap.json
