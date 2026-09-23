# TESTING — the panels, the gates, and what the suite can judge

**What this file is.** The manual for Rigel's test substrates: the simulated panels and what each can
judge (§0), how the test chromosome is built, swept and read (§0a), the toy harness (§0b), the
simulator's backbone and its gaps (§1), the build commands (§2), the simulator's gates (§3), how results
are evaluated (§4), the test suite (§5) and the profiling substrate (§6). Not here: the 0.8.0 scope, the equal-fragment-length ruling and the nascent scope ruling
(`DESIGN.md` §0b); how performance is judged and the instruments' run order (`SUCCESS.md`); the suite's
standing pass count (`CLAUDE.md`); lessons (`TRAPS.md`, cited by name); open problems and refusals
(`ISSUES.md`). The panels' history is git.

Two rules govern the whole file (`TRAPS: prove-the-substrate`, `TRAPS: can-the-benchmark-resolve-it`):

> ⛔ Prove the substrate before you prove the code: when a simulated axis is the axis you are judging,
> gate the simulator on it. And prove the suite can resolve the axis you are changing before quoting a
> number from it.

---

## 0. The panels — one ladder, two fl-gap side panels, and a toy harness

Everything simulated lives under `~/Downloads/rigel_runs/suite/` (the ladder, the two side panels, the
carved reference and the index) and `~/Downloads/rigel_runs/test_reference/` (the test chromosome, §0a).
The panel says how much error there is and where; a toy (§0b) says why.

**The ladder** (`scripts/sim/configs/gdna_ladder.yaml`, `ladder/`) is the only panel the tool is ranked on.

| | |
|---|---|
| judges | Stage B, the calibration solver (`SUCCESS.md`) |
| conditions | 16: gDNA `g00 / g05 / g50 / g98` × strand specificity `0.50 / 0.99` × capture off / on |
| depth | 10 M fragments total per condition, fixed; the gDNA rate decides only the split |
| fragment lengths | identical for both origins (206 ± 98, truncated to [50, 500], read 100) — a forcing function, `DESIGN.md` §0b |
| nascent RNA | sparse on every row: `mode: sparse`, `on_fraction 0.50` of gene spans on, level logU(1, 100) independent of the mature level. 0.50 is a development stress level, not real data (`DESIGN.md` §0b) |

The 0.8.0 scope — three strata in scope, unstranded × capture-ON deferred but reported on every table,
the length composition channel retired, and why the ladder gives gDNA and RNA equal fragment lengths —
is ruled in `DESIGN.md` §0b. Report per stratum, never pooled (`TRAPS: never-pool-the-strata`); the
`g00` rung is the false-positive check on every stratum, and a control is never tuned on.

### The gDNA rungs

`f_gdna = rate/(1+rate)`, and the RNA share thins as gDNA rises because the total is held:

| rung | g00 | g05 | g50 | g98 |
|---|---|---|---|---|
| rate | 0.0 | 0.052632 | 1.0 | 49.0 |
| f_gdna | 0 | 0.05 | 0.50 | 0.98 |
| n_rna | 10.0 M | 9.50 M | 5.00 M | 0.20 M |

`g00` is the required zero-gDNA control, `g98` the top of the range, and `g05` exists because
real libraries live at 1–10 % gDNA and that corner needs a rung with capture on (`0 < rate ≤ 0.10`). Three
levels is a floor (`TRAPS: a-single-level-panel-cannot-see-a-constant`). Real libraries run from almost
zero gDNA to over 98 %, and the RNA-side accuracy that thins at the top rungs is a property of such
libraries, not an artefact.

### The fl-gap side panels — two arms of opposite sign

A side panel, never a ladder rung: the ladder equalises fragment lengths, so it is a null for anything
whose mechanism is a gDNA-vs-RNA length difference, and these two arms exert exactly that. RNA is not
reliably longer than gDNA (true for cfRNA, false elsewhere), and `E_r/E_g − 1` flips sign between the
arms, so a real repair must move them in opposite directions.

| | `flgap_rna_long/` | `flgap_rna_short/` |
|---|---|---|
| config | `scripts/sim/configs/flgap_rna_long.yaml` | `scripts/sim/configs/flgap_rna_short.yaml` |
| gDNA / RNA fl configured | 75 ± 20 / 250 ± 60 | 250 ± 60 / 75 ± 20 |
| gDNA / RNA fl measured | 78.58 / 247.62 (gap +169 bp) | 249.59 / 78.43 (gap −171 bp) |
| conditions | `g50` × ss {0.50, 0.99} × capture {off, on} | the same 4 |
| nascent RNA | `mode: fragment_share`, `shares: [0.20]` | the same |

The configured lengths are a configuration and the gap is a measurement (the sampler is a rejection draw
truncated to `[frag_min, frag_max]`): read it off the payload's deposit histograms. What may be read off these arms is everything that stops
before the EM (`calibration_vs_oracle.py`) and the library gDNA fraction; the
transcript-level number is not a calibration result here, because a length gap hands the EM the answer.

⛔ **Both side panels carry a different nascent model from the ladder** — `fragment_share` at a flat
0.20 on every expressed multi-exon span against the ladder's `sparse`. Each panel's data matches its own
config, so every measurement taken on either arm stands, but a ladder-vs-side-panel comparison varies
two things and no claim may be carried across them (`ISSUES: flgap-panels-stale-nascent-model`).

### The third fl panel — on the test chromosome, on the current nascent model

| | `scenarios_fl_rna_long/` | `scenarios_fl_gdna_long/` | `scenarios_fl_equal200/` |
|---|---|---|---|
| config (`scripts/sim/configs/`) | `test_reference_fl_rna_long.yaml` | `test_reference_fl_gdna_long.yaml` | `test_reference_fl_equal200.yaml` |
| gDNA / RNA fl configured | 100 / 250 | 250 / 100 | 200 / 200 — the control |
| conditions | the full 30 of §0a | the same 30 | the same 30 |

These carry the current nascent model structurally: the test reference is `abundance.mode: file`, so
each config's `nrna:` block is dead by design and nascent comes from `test_abundances.tsv`. The suite
arms' stale nascent model therefore does not block work on the fl gap. The equal-200 arm is the control:
an fl-gap result that also moves on the equal-length arm is an artefact.

### Every scenario must be cached — a requirement of the 0.8.0 loop

The loop is *change calibration → re-run calibration → score it against oracle calibration*, and nothing
in it may re-scan a BAM. Two caches carry it and `panel.py cache` builds both:

| cache | holds | invalidated by |
|---|---|---|
| scan (`build_scan_cache.py`) | the accumulator payload — scan once, calibrate many times | any accumulator change (`payload_schema_digest`), the index (`graph_hash`, `reach_digest`) |
| oracle (`<panel>/oracle_cache/`) | the origin-split truth every scorer reads: five partitions per condition (`gdna` / `mrna` / `nrna` and the per-strand `rna_pos` / `rna_neg`), the undrained `_main` payload, and the certified `slot_truth.npz` | the accumulator or the index, and nothing else |

Neither is invalidated by a calibration change. A scenario without both caches is not usable for
development: `panel.py status` names what is missing, and an instrument fed a stale cache refuses it
(the oracle cache is keyed by the scan cache's own key). The toy harness's donor bundle is the deliberate
exception (§0b): it is a function of the calibration code that fit it, so caching it would serve a stale
answer. How `status` counts a `g00` row is a §2 gotcha.

---

## 0a. The method-development test chromosome — where the message policy is developed

A synthetic chromosome the owner designs one structure at a time, described in one YAML file, swept over
30 conditions, cached, and scored in seconds — so a policy change is scored on every condition between
edits, and the ladder is kept for the shipping judgement.

### The sweep — 30 conditions

`scripts/sim/configs/test_reference.yaml` sweeps three axes; two more requirements are carried by the
substrate rather than the config, so that no condition is uniform in either:

| axis | values | where it lives |
|---|---|---|
| gDNA fraction | `g00` `g05` `g25` `g50` `g98` | the config (`rate = f/(1−f)`) |
| strand specificity | `0.50` unstranded · `0.70` · `0.99` strand-specific | the config |
| capture | off · on | the config |
| which transcripts are captured | about half probed, half with no probe | `test_probes.bed` |
| nascent RNA | sparse — up to about half the transcripts carry it, the rest exactly none, at levels below, near and above mature | the `nrna_abundance` column |

`g00` is the zero control (every reported gDNA fragment is a false positive) and `ss 0.70` is the
transition rung between the two regimes the bars below are written for. The nascent pattern is
hand-authored, not drawn: on a chromosome this small a draw is mostly sampling noise, and the point of a
debug substrate is knowing which structure carries nascent RNA before reading a number. Depth is
`n_total_fragments: 1170000` per condition, raised with each block so per-transcript depth holds. The
config's `nrna:` block is dead by design (`abundance.mode: file` makes the rendered TSV the one source).

### The commands — from the YAML to a scored benchmark

Everything derived is rebuilt from the one YAML. Run in order; each stage is resumable and
`panel.py status` names the next one. The `cache` stage builds both caches, pre-warms the `g00` rows
(held out of the oracle sweep) and copies their `_main` from the scan cache, then certifies
(`calibration_oracle.py` writes `slot_truth.npz`; a failed FIELD gate is reported, not fatal — the row is
COMPOSITION-certified). One row by hand is a §2 gotcha.

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
T=~/Downloads/rigel_runs/test_reference

# 1. renders (GTFs, abundances, probe panels — into the repo AND the runs dir) + the two-contig FASTA;
#    `--check` alone says whether the checked-in renders still match the YAML
python scripts/sim/build_test_reference.py
rigel index --fasta $T/test_chr.fa --gtf $T/test_chr.gtf --collapse-duplicate-transcripts --no-mappability --no-tsv -o $T/idx

# 2. per panel (seven configs: the benign panel, the two adversarial probe panels, od05, three fl arms)
#    ⛔ ONE panel's `cache` at a time, or give each its own RIGEL_SCRATCH: the origin split writes
#    $RIGEL_SCRATCH/rigel_oracle_build/w_<condition>/<condition>.<origin>.bam, and six of the seven configs share their
#    condition names, so two caching side by side corrupt each other's partition (sum-to-full fails).
CFG=scripts/sim/configs/test_reference.yaml
python scripts/sim/panel.py simulate --config $CFG
python scripts/sim/panel.py cache    --config $CFG

# 3. score every policy on every condition — seconds
python scripts/design/policy_benchmark.py --panel test
```

⛔ **Rebuild after every edit to the YAML.** The FASTA, index, reads, caches and `slot_truth.npz` all
describe the annotation they were built from; a benchmark against stale caches answers a different
question (`TRAPS: a-green-suite-hid-five-dead-instruments` in another dress). Move a superseded derived
set out of the way the same day.

### How to read the benchmark — two halves, two bars

The bar is not "beat silent". `SilentPolicy` is the measured floor, and on strand-specific data a
sighted exon's own strand solve is excellent, so a message can mostly only disturb it. Message
propagation exists for unstranded data, where the strand channel is dead and the local answer is a
default rather than a measurement.

| rows | the bar |
|---|---|
| unstranded (`ss 0.50`) | the policy must win — this is what the message layer is for |
| stranded (`ss 0.99`) | the policy must do as little harm as possible — near 1.00× of silence is a pass |

⛔ Never pool the two halves, and never pool conditions within them: a total hides a sign flip between
strata, and the halves are judged against different bars. A toy and the panel can disagree in rank
(`TRAPS: a-toy-and-a-panel-can-disagree-in-rank`): develop on the test chromosome, confirm on the ladder.

### The substrate — one hand-edited file, twelve blocks

`scripts/sim/test_reference/test_chr.yaml` is the one hand-edited file: the `rigel sim` scenario schema
(`genes → {gene_id, strand, transcripts: [{t_id, exons, abundance, nrna_abundance}]}`, exons 0-based
half-open) plus `probed` per gene and a `shadow_genes` list. Its header is the design record of every
block; this table says only what each one stresses.

| block | what it stresses | gene types |
|---|---|---|
| anchored twin | the message layer's five controls on one shape (exon 1 kb · intron 7 kb · exon 1 kb · intron 7 kb · exon 1 kb), each type varying one bit — probed, nascent, silent — × five abundance blocks; plus 8 shadow transcripts the index never sees | `clean` `nasc` `cap` `capnasc` `silent` |
| mono | single-exon: no sj, the gene edge is the only message source | `mono` `capmono` + two silent controls |
| isoform | host + one second isoform, grown one structure at a time, every structure replicated across A ≫ B, A ≪ B, A ≈ B | `altstart` `altss` (`nest` queued, `ISSUES: message-layer-open-cases`) |
| walled | exon pieces with no licensed face — the scan's stress test, the four classes the ladder's walled-exon census found | `chain` `tssalt` `tandem` `altlast` |
| terminus-cluster | ten transcript ends 126–147 bp into a shared last exon (mirrored from MIR99AHG) — the empty exon pieces the level lane crosses | `cluster` `capcluster` |
| both-stranded | two genes per locus on opposite strands, the host and its antisense (mirrored from TTC28-AS1, PPM1F-AS1, and a convergent pair) — the AMBIG nodes' substrate | `asin` `asinrev` `span` `conv` |
| sj+terminus | one boundary carrying a junction and a terminus of the same strand (mirrored from RUNX1 and LARGE1) | `sjterm` `capsjterm` |
| encompassing | a 20 kb single-exon antisense over the whole twin-shape host (mirrored from ENSG00000280007 over TUBA8): every host slot admits both strands, every boundary carries only the host's bits, the antisense's level lives in its two 1.5 kb single-strand flanks — the owner's encompassing locus on the panel | `enc` `capenc` |
| in-exon | a 500 bp single-exon antisense wholly inside the host's 3 kb last exon (mirrored from ENSG00000273300 in UFD1's 3' UTR): no junction and no single-strand piece of its own, the tilt atom's accepted limit | `inexon` `capinexon` |
| shared-exon | two spliced genes sharing one exact 10 kb last exon (the W12 deep stress; the tail-to-tail class of SMARCB1 × DERL3): a level on each strand into one walled AMBIG exon whose faces each pair a junction with a terminus | `shared` `capshared` |
| in-intron | a 600 bp single-exon gene centred in the opposite strand's first intron (mirrored from the ladder's 203 intronic pseudogenes and lncRNAs): the majority AMBIG class, walled, its only message source its own termini | `inintron` `capinintron` |
| head-to-head | `conv` with the strands swapped and nothing else: two 5' ends overlap in a 2 kb exon∩exon piece (mirrored from TRMT2A ⟷ RANBP1), a TSS and a donor of different strands on one boundary | `div` `capdiv` |
| tiny-exon | exons shorter than a fragment — ten 40 bp exons at 1,040 bp pitch, and the same run between two 1 kb exons: no piece holds a contained fragment, so a piece is seen only through its boundaries: the crossing counts at its edges price the conserved shares of the cuts around it (`EQUATIONS.md` §11) | `tiny` `captiny` `mixed` `capmixed` |

273 genes on a 7.930 Mb chromosome (`genome_length` in the YAML). Every gene carries an explicit strand
and the chromosome keeps equal + / − representation — a sign error is invisible on one strand, and the
builder refuses an imbalance. Abundances are molar ladders in half-decade steps with mature up the
blocks and nascent down (10/30/100/300/1000 against 100/30/10/3/1), independent levels, as the ladder
draws them.

`shadow_genes` sit on the blank contig `test_blank`: unannotated transcription the simulator draws from
and the index never sees — the control for anything that trusts the annotation's notion of "pure gDNA".
The GTFs, `test_abundances.tsv` and the three probe BEDs are rendered by `build_test_reference.py` and
versioned beside the YAML (`--check` and `tests/test_test_reference_renders.py` refuse a drifted render);
the benign panel tiles each probed gene's exon union 8 × 125 bp, never spanning an sj (a union piece shorter
than one probe gets a single probe centred on it, spilling into its flanks as a real panel's does), the sparse panel
centres one probe per exon, the junction panel places one two-block probe per sj. The FASTA
(`test_chr.fa`, both contigs), the index, the reads and the caches are derived, never hand-edited: a
spliced transcript needs a GT..AG at every intron or the aligner and the simulator disagree with the
annotation, and generating the chromosome from a fixed seed with a motif injected at every declared
intron makes that impossible to get wrong. The builder refuses a reference it cannot simulate and reports
every problem at once (an exon off the end, overlapping exons, an intron under 4 bp, a duplicate id, an
abundance row without a transcript or the reverse, the wrong reference name); `--self-test` perturbs each.

⚠ κ is fitted from spliced reads, so a transcriptome of only single-exon transcripts cannot calibrate
(`CalibrationStrandError`): at least one multi-exon transcript with real depth is needed before a number
means anything. Nascent RNA needs no declaration: `rigel index` creates a single-exon nascent entity
over each multi-exon transcript, fed from the contributor's `nrna_abundance` column.

## 0c. The ruler's truth instrument and the gDNA-depth ladder

`scripts/design/ruler_vs_truth.py` scores the EM's effective length under capture — the ruler,
`capture_eff_length.transcript_capture_eff_lengths`, each object a transcript's fragments deposit on at its
capture efficiency (`EQUATIONS.md` §11) — against the simulator's own capture-aware effective
length, per transcript: `CaptureSampler.partition_array` is what drew the reads, so the truth is
`Σ_w f_pre(w) · partition_t(w)` over the pre-capture fragment-length pmf and the plain length is
`Σ_w f_pre(w) · off_target_weight · (L_t − w + 1)+`; their ratio is the truth factor, an arm's factor is
its effective length over the same plain lengths, and the two are compared in the log anchored on the
fully probed transcripts (a global scale is free) per class — the probed fraction of the transcript's
bases, read off the sampler's partition at a one-base fragment — and per kind (mRNA / annotated single-exon
/ synthetic nascent entity). It reads the same caches as every other instrument and calibrates in
seconds; `--condition` prints the per-class table, without it every capture-ON condition prints one line
(the probed class's share within ±0.1 nat, the unprobed and partial classes' median error, per arm).

```bash
python scripts/design/ruler_vs_truth.py --panel test --condition gdna_g05_ss_0.99_nrna_file_capture_on
python scripts/design/ruler_vs_truth.py --panel ladder                       # one line per capture-ON row
python scripts/design/ruler_vs_truth.py --panel test --condition C --module proto.py --out table.tsv
python scripts/design/ruler_vs_truth.py --panel test --panel-dir ~/Downloads/rigel_runs/test_reference/scenarios_depth_d10
python scripts/design/ruler_vs_truth.py --panel ladder --condition gdna_g50_ss_0.99_nrna_mid_capture_on --scale
```

The arms: `shipped` is the ruler in `src/` on the shipped calibration; `oracle_gdna` feeds the same ruler
the certified true gDNA counts per object (`slot_truth.npz`) with the reference read off a landscape
fitted on the truth — the ideal witness, so the gap between the two is the calibration's and what remains
under `oracle_gdna` is the ruler's own or the panel's geometry; `--module proto.py` names a file defining
`ARMS = {name: ruler}` with `ruler(calibration, region_arrays, index, fl, **inputs)` — `inputs` is
`calibrate`'s debug bundle (the last refit's `DensityLandscape` under `gdna_hyperprior`, the prior an
expectation arm reads; the chain and the belief) plus the two fragment-length pmfs — run beside the shipped
one on the same calibration
— DERIVE → PROTOTYPE happens there, and nothing in `src/` moves to price an arm; `--set SECTION.FIELD=VALUE`
prices a config value on every arm. Read the classes apart: the probed class is where the formula is exact
when its witness is; the unprobed class is where a floor or the reference bites; the partial classes are
where the junction price lives. A probe spanning a junction or centred on an exon shorter than a fragment
is captured on gDNA at a fraction of the cDNA's overlap. The junction price reads that extra capture from
the gDNA objects beside the junction by conservation of bases — close on average, too noisy per junction to
split isoforms (`ISSUES: the-junction-price-is-noisy-within-a-gene`) — and what no gDNA object sees is
declared, not repaired (`ISSUES: ruler-witness-geometry-on-transcript-panels`). `--scale` reads the
shipped lengths: every hypothesis class's L / Y unanchored and the within-gene spread of the
multi-exon transcripts, so a prototype length is scored on it from a worktree, and a class-scale repair is
judged by the spread as well as the means (`TRAPS: judge-a-ruler-by-its-within-gene-spread`).

**The gDNA-depth ladder** is three side configs of the test chromosome, everything but the fragment
budget and the gDNA rungs identical to `test_reference.yaml`, stranded only, built the same way (§0a's
commands with `--config` pointing at each; one `cache` at a time or a `RIGEL_SCRATCH` each):

| config | depth | gDNA rungs | what it answers |
|---|---|---|---|
| `test_reference_depth_d10.yaml` | 117 k fragments | `g00 g001 g01 g05 g25 g50` | where the reference is `None` and where the probed class reads within ±0.1 nat, at a tenth of the depth |
| `test_reference_depth_d100.yaml` | 11.7 k | the same six | the same curve at a hundredth — the axis is gDNA fragments per probed piece, not the fraction |
| `test_reference_depth_full_lowg.yaml` | 1.17 M | `g001 g01` | the two rungs the panel lacks at full depth |

`ruler_vs_truth.py --panel-dir` scores each set; the curve is one line per condition. It is a side panel
in the sense of §0 — its transcript-level number is not a calibration result — and never a ladder rung.

---

## 1. The simulated panel (shared backbone)

A real human backbone, not a generated mini-genome: `chr21` + `chr22` + the 92 ERCC spike-in references,
carved from the same GRCh38 + GENCODE v46 sources the production index is built from. Neither is filler:
a single-reference synthetic index once hid a reference-id-space mismatch that silently dropped nearly
every fragment inside `deposit()` while every golden test passed
(`TRAPS: one-reference-hides-refid-bugs`); the spike-ins make the id space non-trivial for RNA, two
genomic chromosomes make it non-trivial on the gDNA path too.

Every fragment-length parameter is measured — the count-weighted mean and sd of a real cfRNA library's
own RNA pool, 206.1 ± 98.3, handed to both origins (`DESIGN.md` §0b). Those are pre-capture parameters
and the truth files are post-capture: hybrid capture selects for length, so the simulator draws the
length marginal as `f_post(w) ∝ f_pre(w) · total_eff(w)`. ⛔ Score against `truth_fragment_lengths.tsv`,
never against `frag_mean` (`TRAPS: capture-selects-for-length`).

### The simulator's transcriptome is a Rigel index — including the nascent entities

`rigel sim` builds (or is given, config key `index:`) a rigel index and simulates its transcript list, so
what is simulated is exactly what `rigel quant` reads: annotated transcripts plus the synthetic nascent
entities `index.create_nrna_transcripts` makes (one single-exon transcript over each multi-exon span,
TSS/TES clustered within `NRNA_MERGE_TOLERANCE`, an annotated single-exon transcript adopted where one
covers the span). `transcript_filter` is refused — filter the GTF before building the index.

| | |
|---|---|
| nascent RNA is a transcript, not a parallel space | its molecules are `entity.nrna_abundance`, sampled on its own template, its reads keep the `nrna_` origin tag. Under `sparse` the level is drawn per entity — off with probability `1 − on_fraction`, else log-uniform over `abundance_ranges`, independent of the mature level; under `additive_ratio` and `fragment_share` it is `Σ abundance × nrna_ratio` over contributors |
| one multinomial over every RNA row | mature and entity rows together, `prob ∝ abundance × capture-aware effective length`; the mature/nascent split follows from molecules and lengths, so a nascent-on condition and its nascent-off twin do not share a bit-identical mature stream |
| capture binds by genomic overlap | every probe → genomic blocks → gDNA, projected onto every transcript whose exons it touches, either strand (ds-cDNA at capture). A fragment binds through ONE contiguous part of a probe: a transcript holding the junction a probe spans holds it whole, while gDNA, a nascent span and an isoform without the junction hold its parts apart and bind the better one — never the sum, and with no other penalty, since the part a molecule holds binds alike whatever the molecule (owner, 2026-09-19; it replaced a `gdna_split_penalty` of 0.2 that bound a gDNA half-match at a fifth of the identical cDNA one). gDNA and nascent are enriched at the same rate under one probe, which is the physics |

Two consequences ruled correct (owner, 2026-08-19): capture depletes nascent RNA about 8× (roughly 20 %
of RNA fragments off capture, under 3 % on — read a condition's `truth_summary.json`), because probes
tile exons and a pre-mRNA is mostly intron, and adding nascent after the capture simulation is wrong —
molecules exist first; and the post-capture fl distribution and abundances are the ground truth
(`truth_kind: post_capture_empirical`), which is why gate G-S5 (§3) is stated per RNA population.

### What the simulator still does not do

| | consequence |
|---|---|
| counts are Poisson by construction | nothing dispersion-dependent validates here |
| the panel is all R1-antisense | the engine can emit either (`ReadSimConfig.r1_sense`, gated in `test_strand_sense_convention.py`) but `orchestrator.run_condition_grid` does not expose it. Real cfRNA is dUTP, so this is not urgent |
| the tool's gDNA reach assumption (`taper_g = 1`) is untested | latent: gDNA is not simulated on the spike-ins |
| each population is written as one contiguous block of read names | a per-fragment truth join checked by "does an impossible label appear?" is nearly blind. ⛔ Gate such a join on a count identity against the scanner's own `stats.total` / `stats.n_read_names` (`_oracle.check_walk_alignment`); `tests/calibration/test_prior_vs_oracle.py` pins both halves |

---

## 2. Building it — `panel.py`, one command per stage

`scripts/sim/panel.py` is the whole loop: every path derived from one panel YAML, every stage resumable,
each stage shelling out to the instrument that owns it. It adds no measurement code.

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
export OMP_NUM_THREADS=1
CFG=scripts/sim/configs/gdna_ladder.yaml       # the panel the tool is ranked on

python scripts/sim/panel.py status   --config $CFG              # always first
python scripts/sim/panel.py build    --config $CFG              # index + capture probes
python scripts/sim/panel.py simulate --config $CFG --jobs 8     # ~21 min, ~16 GB
python scripts/sim/panel.py cache    --config $CFG --jobs 8     # BOTH caches — mandatory, §0
python scripts/sim/panel.py score    --config $CFG --jobs 8 --arms base oracle noop
python scripts/sim/panel.py report   --config $CFG --arms base oracle
```

`status` prints what exists, what is missing, and which stage to run next. `build` cannot carve the
reference — that needs the source genome/GTF, which the config does not name — so it prints the exact
command and stops:

```bash
python scripts/sim/build_suite_reference.py \
    --fasta $REFS/genome_controls.fasta.bgz --gtf $REFS/genes_controls.sorted.gtf \
    --refs chr21 chr22 --ercc -o $SUITE/reference
```

### Two workflow gotchas

* ⛔ **`--force` does not reach the `simulate` stage.** `cmd_cache` passes it to `build_scan_cache.py`
  and `cmd_build` honours it for the index and probes, but `cmd_simulate` shells out to
  `simulate_reads.py` with no such flag, and the simulator skips a condition whose oracle BAM already
  exists (`skip_existing`). So editing a config and re-running `panel.py simulate --force` reports
  success and reproduces the old reads. Delete the condition directories first.
* `calibration_oracle.py --build` builds every row's oracle cache alike — the zero-gDNA rows too, there is no
  hold-out — and `status` counts an oracle condition complete only when `gdna`, `mrna`, `nrna`, the two strand
  partitions and `_main` are all present. One row by hand:

  ```bash
  python scripts/design/calibration_oracle.py --suite $SUITE/ladder --index $SUITE/rigel_index --build \
      --condition gdna_g00_ss_0.50_nrna_mid_capture_off
  ```
* Every panel config states `gdna.genomic_refs: [chr21, chr22]` explicitly; the engine does not infer
  which references carry genomic DNA (`TRAPS: annotated-is-not-genomic`).
* Export `RIGEL_SCRATCH` before a sweep, or the instruments write their per-condition work under `/tmp`
  (tens of GB after a ladder rebuild).

---

## 4. How results are evaluated

The 0.8.0 metric is the calibration result scored against oracle calibration, not the end-to-end
transcript number (`DESIGN.md` §0b; the judging rules are `SUCCESS.md`). Three questions, three
instruments:

| question | instrument |
|---|---|
| how wrong is calibration, against oracle calibration? — the 0.8.0 metric | `calibration_vs_oracle.py` · `prior_vs_oracle.py` (the `LocusPriors` the EM reads) · `solvability_audit.py` (pass-0, per object) |
| how wrong is the end-to-end answer? — per transcript, per pool, against per-fragment truth | `panel.py score` / `report`, i.e. `quant_accuracy.py` — the only end-to-end scorer, a thermometer rather than the target |

⛔ A ceiling only prices what its arm can reach: the effective-length shrinkage is built before
`assemble_priors`, which every measurement arm patches (`SUCCESS.md`, the ruler). Say which call your
arm patches and check it sits downstream of everything you mean to price.

Hard per-fragment label recovery is the wrong target: an unspliced RNA fragment and a gDNA fragment from
the same locus can be sequence-identical. A net-flow reduction, `net(a→b) = flow[a][b] − flow[b][a]`, cancels
the unrecoverable part and keeps only systematic bias; the module that computed it per locus was retired
2026-09-13 for want of an entry point (in git), and the absolute per-transcript error stays the report. Hard-label
metrics are nearly blind to a
calibration-prior change (`TRAPS: hard-labels-miss-soft-change`); the soft 3-pool surplus is not built.

---

## 5. The test suite

```bash
python -m pytest tests/ -q
python -m pytest tests/ --update-golden        # regenerate tests/golden/ after intended output changes
```

`CLAUDE.md` is the home of the standing pass count and accounts every delta: re-derive it, never adjust
it, because several gates are parametrised over the files on disk. Any failure is a regression.
`tests/scenarios_aligned/test_multimap_counting.py::TestParalogMultimapping::test_gdna_sweep[gdna_100]`
is a real EM unidentifiability (`TRAPS: identical-paralogs-are-bimodal`); if it fails again, do not fix
it by moving a seed.

The goldens run under the default sampling mode (`EMConfig.seed = None`, `assignment_mode = "sample"`),
and two runs of the identical pipeline on the identical BAM return different transcript counts
(`TRAPS: the-deliverable-is-not-reproducible-by-default`). Regenerate the goldens twice and diff, and pin
`EMConfig.seed` in any instrument that compares two end-to-end runs (`quant_accuracy.py` does, and prints
a `base_reseed` noise floor beside the effect).

**The capture-contracted length is gated per object against the deposit rule.** Three gates run every
placement of every fragment width through the reference accumulator (`tests/native/_accumulator_reference.py`)
and hold each piece's contained share and each boundary's and junction's conserved share to 1e-12:
`tests/calibration/test_effective_length.py::test_each_share_is_what_the_deposit_rule_gives_the_object` (a
spliced template cut into 1–3 bp pieces, one fragment crossing up to six cuts),
`::test_the_gdna_shares_are_what_the_deposit_rule_gives_each_boundary` (gDNA, every start on the chromosome)
and `tests/calibration/test_capture_eff_length.py::test_each_share_is_what_the_deposit_rule_gives_the_object`
(the index's own partition). Per object, never on the total (`TRAPS: conservation-misses-mis-attribution`),
and on pieces shorter than a fragment (`TRAPS: perturb-every-gate`). Two more stand beside them:
`test_capture_eff_length.py::test_a_junction_reads_the_objects_beside_it_on_a_multi_region_intron` puts
every region and boundary at its own efficiency and finds the junction's objects by coordinate, never by the
module's index arithmetic, and `::test_a_zero_length_fragment_places_nowhere` holds a length model's mass at
`w = 0` out of every share. `test_priors.py` holds the gDNA component's length to gDNA's own conserved share
and never the count's `q`; `test_capture_efficiency.py` holds each region and boundary to its own count.

---

## 6. Development discipline for test substrates

Develop on controlled toys, validate on real data — both, in that order: a big suite has confounds that
hide mechanisms, and a toy ranks hotspots backwards (`TRAPS: toys-rank-hotspots-backwards`).

**Profiling is judged on a deep, real, high-complexity RNA-seq library** (owner, 2026-08-17) — never a
panel condition and never a toy. The cfRNA libraries on disk are sparse and small, so they under-fill
exactly the structures whose cost the performance work is about (the grid solve, the per-slot arrays, the
fragment buffer); they are smoke tests that find a defect cheaply, and the deep library decides. This
leaves `TRAPS: real-data-is-a-test-input` and the accuracy panel (the ladder) untouched.

**Read a timing only from back-to-back A/B pairs.** Stages nobody touched drift 25–40 % between runs
taken at different times — page cache, concurrent work — so a before-and-after taken hours apart credits
the machine's mood to the change. Alternate the two arms in one sitting (stash `src/` only) and check
that the untouched stages read 1.00. **Every speed-up is proven a numeric no-op** before it is believed:
`design/rename_identity.py --bam` end to end on a real library and on both ladder references, and
`profiling/sweep_replay.py` for anything inside the sweep.

The instruments are `scripts/profiling/profiler.py` (the whole pipeline as a tree of named stages, with
per-stage peak and held RSS; `--set` for any config field, `--scan-only` for the scan alone, `--compare`
for two reports) and `scripts/profiling/sweep_replay.py` (one calibration sweep, replayed and compared
bit for bit; `--block-slots N|none` replays it at another locus-block size, which must move nothing —
the chunk-exactness of the whole sweep on real data); set `OMP_NUM_THREADS` deliberately. The frozen
references are `~/Downloads/rigel_runs/arms/review_identity_<condition>.json` — two ladder conditions and the
LBX0190 library, last re-frozen at `71037f6c` with the reason in `arms/refreeze_2026-09-21.log`. The instrument's
default `--reference` does not exist and it prints ⛔ rather than failing, so pass each one and read the log, never
the exit code:

    A=~/Downloads/rigel_runs/arms
    for c in gdna_g05_ss_0.99_nrna_mid_capture_on gdna_g05_ss_0.50_nrna_mid_capture_off; do
      python scripts/design/rename_identity.py --check --condition $c --reference $A/review_identity_$c.json
    done
    python scripts/design/rename_identity.py --check --index ~/Downloads/rigel_runs/refs/rigel_index \
      --bam ~/Downloads/rigel_runs/cfrna/mctp_LBX0190_SI_43883_HHHKGDRX7/bam/star.srt.rmdup.collate.bam \
      --reference $A/review_identity_LBX0190.json

The captured sweeps are `~/Downloads/rigel_runs/perf/sweeps_VCaP_step19` (the deep library, VCaP) — the replay
target; a step that moves numbers re-captures and deletes the superseded one (`DESIGN.md` §6b.15). A port is
held to the replay's `--tolerance` budget and to the transfer gates, which hold the kernel's tables — returned by
`native.transfer_prepare` — to independent recomputes and drive single hops of the native pass through
`native.transfer_pass` (`tests/calibration/_transfer_harness.py`, where the tables' containers live; the row constructors
and flag predicates the gates recompute with are the native ones, `native.transfer_rows`). The sweep itself is one native
call (`native.solve_blocks`): `sweep_replay.py replay --threads N` and `--block-slots N` hold it to a capture at any thread
count and block size. The deep library's timing baseline is
`perf/plan_final_2026-09-19/pair{1,2}_post.json` (the landed tree at 8 threads on VCaP, the post arm of two
interleaved pairs; the pre arm beside it is the tree before the work outside calibration, and the drift between two
runs of one arm, 0.98–1.07 per stage, is the noise floor a pair is read against). ⛔ Read a stage row from a PAIR,
never a wall from two sittings, and never a profiler's share as a saving (`TRAPS: a-profile-share-is-a-ranking`). Every
captured sweep replays the whole message layer at its own bracket (the refit sweeps at K = 233 against the first
sweep's 101). `ISSUES: performance-memory-bounded-solve` carries the work.

A both-strand stress test needs ample single-stranded regions (the population prior trains on them).
How to A/B honestly: in-process, opposite extremes, never on a saturated condition, one thing varied,
both arms sharing their random input (`TRAPS: perturb-every-gate`) — a byte-identical hard-label result
is no evidence.
