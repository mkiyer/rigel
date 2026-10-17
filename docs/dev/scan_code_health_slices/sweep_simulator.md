## Scope read
- Every file under `src/rigel/sim/`, read in full: `__init__.py`, `annotation.py`, `bam.py`, `benchmark.py`, `genome.py`, `intervals.py`, `manifest.py`, `orchestrator.py`, `read_name.py`, `reads.py`, `sampling.py`, `scenario.py`, `splice_motif.py`, `truth.py`, `wgs_config.py`, `wgs_engine.py` (1–1343), `whole_genome.py` (1–1252), and `capture/{__init__,config,design,sampler}.py`.
- Context:
  - `docs/ISSUES.md`: `hygiene-ledger`, `latent-defects`, `flgap-panels-stale-nascent-model` and `overlapping-synthetic-shadows`, plus greps. `docs/TRAPS.md`: greps.
  - Source: `config.py:763–781`, `tests/test_constants.py`, `cli.py:735–830` and `:1465–1510`.
  - Scripts: `scripts/sim/panel.py:155–200`, `design_suite_probes.py`, `simulate_reads.py`, and excerpts of `quant_accuracy.py` and `ruler_vs_truth.py`.
  - Tests: excerpts of `tests/calibration/_oracle.py`, `test_quant_accuracy.py`, `test_sim_truth.py` and `test_sim_gdna_overdispersion.py`.
  - Data: the ladder's `capture_panel.tsv`, plus `transcripts.tsv` and `transcript_aliases.tsv` from its index. I read these with `python3 -I` and csv only. Nothing from rigel was imported, built or run.

## Findings

**1. A probe in transcript coordinates that names a transcript the index dropped is skipped silently, and this happens on the ladder**
- **Location:** `capture/sampler.py:514–517`
  ```
          t_idx = self._tx_id_to_index.get(transcript_id)
          if t_idx is None:
              logger.debug("Skipping probe for unknown transcript %s", transcript_id)
              return
  ```
- **Category / severity:** contract · medium
- **What is wrong:**
  - The ladder's panel is designed on the raw GTF (`design_suite_probes.py:37`, `load_transcripts(args.gtf)`). It is simulated against an index built with `--collapse-duplicate-transcripts` (`panel.py:173`).
  - A probe placed on a dropped duplicate is discarded at DEBUG level. `_load_transcript_probes` still counts it in its INFO line "Loaded N probes", because `n_loaded += 1` (`:476`) runs after the silent return.
- **Evidence:**
  - Of 13,824 rows in `$RIGEL_RUNS/suite/reference/capture_panel.tsv` (the ladder uses `probe_format: transcript`), 2 name t_ids that are not in the index: `ENST00000397293.6` [1841,1961) and `ENST00000703226.1` [277,397).
  - Both appear in `transcript_aliases.tsv` as dropped duplicates. Their kept twins carry different probes ([847,967) and [145,265)), because the designer masks probes already placed. So the loss is real coverage, not a redundant copy.
  - No test covers the skip.
  - `ruler_vs_truth.py` and `quant_accuracy.py` skip the same two probes, so they agree with the simulator. Only the designed panel differs.
- **Change:** refuse an unknown t_id, or resolve it through `transcript_aliases`. Design the panel on `load_transcripts_from_index`, the reason that function's own docstring gives for not using the GTF loader.
- **Moves numbers:** yes (ladder capture-ON conditions; re-simulation).
- **Confidence:** confirmed.

**2. The deleted `fragment_share` mode's `shares` key is still on the allow-list, so a config using it runs a nascent-free scenario without any error**
- **Location:** `whole_genome.py:130–132`
  ```
      "nrna": frozenset(
          {"mode", "shares", "ratios", "abundance_ranges", "ratio_labels", "on_fraction", "seed"}
      ),
  ```
- **Category / severity:** defect · medium
- **What is wrong:**
  - Nothing reads `shares`, and `_MODE_FIELDS` (`:320–323`) does not list it.
  - So `nrna: {shares: [0.2], ratio_labels: [s20]}` with `mode` omitted parses cleanly and runs `additive_ratio` at ratio 0.0, labelled `s20`. That is exactly the failure the comment at `:316–319` says the parser refuses.
  - The CLOSED entry `flgap-panels-stale-nascent-model` records `shares` as deleted.
  - Same kind of residue: `transcript_filter` (`:101`, `:259–263`) is accepted silently when set to `all`. The `format`, `name` and `capture.scenarios` aliases (`:78–90`, `:186`, `:202`, `:214`) are used by no config; only `format:` appears, in `tests/test_sim_capture.py`.
- **Evidence:** `rg -n shares src/rigel/sim` finds only this set and the stale comment at `wgs_config.py:105`. A test that parses such a YAML and expects ValueError would fail today.
- **Change:** drop `shares` and `transcript_filter` from the allow-lists, so the unknown-key refusal names them; drop the unused aliases.
- **Moves numbers:** no.
- **Confidence:** confirmed by reading. Not executed, because importing `rigel.sim` loads the native modules.

**3. Two truth tables with two writers; one column name carries two units, and each count is written three times**
- **Location:** `truth.py:210–215`
  ```
                  "post_capture_mrna_fragments": mrna_count,
                  "post_capture_nrna_fragments": nrna_count,
                  "post_capture_total_rna_fragments": total_count,
                  "mrna_abundance": float(mrna_count),
                  "nrna_abundance": float(nrna_count),
                  "total_rna": float(total_count),
  ```
- **Category / severity:** duplicate · medium
- **What is wrong:**
  - Each condition's `truth_abundances.tsv` holds every realised count under three names: `mrna_abundance`, `post_capture_*` and `observed_*`. The docstring (`:291–294`) keeps `mrna_abundance` "compatible with existing analysis code".
  - The suite-level `truth_abundances_nrna_<label>.tsv`, written by a second writer (`whole_genome.write_truth_abundances`, `:926–947`), uses the same names for **molar** abundance.
  - `_detect_abundance_format` (`whole_genome.py:633`) accepts either file as molar input.
- **Evidence:**
  - `test_quant_accuracy.py:229–237` exists only to catch scoring against the wrong file.
  - `quant_accuracy.py:385–388` keeps a fallback pair for truth tables that "predate the observed columns" (see the test at `:940–942`).
  - `post_capture_*` is read only by `tests/test_sim_truth.py`.
  - Commit 23e22cb8 fixed the same genomic-span bug in both writers but gated only `truth.py`'s. `write_truth_abundances` has no test at all, so reverting its fix would stay green.
  - The same doubling appears in the manifest entry (`orchestrator.py:259–264`: three paths, each under two names) and in `truth_summary.json` (`truth.py:304–318`, repeated under `files`). No instrument reads either; the instruments list the condition directories on disk.
- **Change:** one writer in `truth.py` for both tables, with column names that state the unit (molecules versus `observed_*` fragments), no aliases, and one key per path. Retire the instrument's fallback. This sits next to `hygiene-ledger` (h), which covers the readers.
- **Moves numbers:** no.
- **Confidence:** confirmed.

**4. The read-name truth encoding has no single owner**
- **Location:** `wgs_engine.py:613`, `:668`
  ```
          qname_prefix = f"nrna_{t_id}" if is_nrna else t_id
                  f"{qname_prefix}:{int(frag_starts[i])}-{int(frag_starts[i]) + frag_len}:{strand_char}:{base_n + i}"
  ```
- **Category / severity:** duplicate · medium
- **What is wrong:**
  - The names are built by inline f-strings in the engine (`:613`, `:668`, gDNA at `:889–891`, strand chars at `:291–293`) and parsed in `read_name.py`.
  - The docstrings that state the format (`read_name.py:3–4`, the module that owns it; `reads.py:11–12`; `scenario.py:191`) all omit the `nrna_` prefix and the `f`/`r` strand alphabet.
  - That alphabet is enforced only in test tooling (`tests/calibration/_oracle.py:108–113`). `parse_origin` accepts any strand, and `tests/test_sim_shadow_transcripts.py:73` passes it `+`, which the engine never writes.
  - `truth.py:20` re-exports `parse_origin`, giving a second import path.
- **Evidence:** `rg 'nrna_|gdna:'` finds the prefixes as literals in `wgs_engine.py`, `read_name.py` and `tests/test_oracle_bam.py:240–242`.
- **Change:** `read_name.py` owns `format_rna_name`/`format_gdna_name`, the two prefixes and the strand alphabet, which `parse_origin` validates. The engine calls these functions. Delete the docstring copies and the re-export.
- **Moves numbers:** no (the read names stay bit-identical).
- **Confidence:** confirmed.

**5. "Shadow" means two opposite things**
- **Location:** `wgs_config.py:124–125`
  ```
      #: Shadow transcripts: a supplemental GTF the simulator draws fragments from but the index never
      #: sees — unannotated transcription, simulated. Their ids must be unknown to the index (a shadow
  ```
- **Category / severity:** naming · medium
- **What is wrong:**
  - In the simulator (`shadow_gtf`, `merge_shadow_transcripts`, `n_shadow_transcripts`), a shadow is transcription the index does not know about.
  - In the EM, the index and the instruments, a shadow is the synthetic nascent span the index itself manufactures: `ISSUES: overlapping-synthetic-shadows`, `em_solver.cpp:593–606`, and `quant_accuracy.py:380–383` ("a synthetic shadow entity emits no mature RNA").
  - The truth files the two under different origins (`mrna` versus `nrna`).
- **Evidence:** `rg -i shadow` across src, scripts, tests and docs shows both senses; ISSUES:1161–1163 uses the simulator's.
- **Change:** rename one sense, grepping both first and checking with `rename_identity.py`.
- **Moves numbers:** no.
- **Confidence:** confirmed.

**6. `oracle_bam: false` together with `emit_fastq: false` gives a panel with no reads, which is re-simulated on every run**
- **Location:** `orchestrator.py:274`, `:322–325`
  ```
                      _exists_key = cond_dir / ("sim_oracle.bam" if oracle_bam else "sim_R1.fq.gz")
  ```
  ```
                      if not emit_fastq:
                          for _fq in (cond_dir / "sim_R1.fq.gz", cond_dir / "sim_R2.fq.gz"):
                              if _fq.exists():
                                  _fq.unlink()
  ```
- **Category / severity:** contract · low
- **What is wrong:**
  - Nothing refuses this combination: not `whole_genome.py:381–386`, and not `main()` after `--no-oracle`.
  - Each condition is simulated and its truth written. Then the FASTQs are deleted, and they are both the only reads and the resume key.
  - The run exits 0 with no reads. The next run repeats everything. `panel.py:81–85`, which lists only directories holding `sim_oracle.bam`, never sees these conditions.
  - The manifest still names `fastq_r1/2` after they are deleted (`:266–267`).
- **Evidence:** the code path above.
- **Change:** refuse the combination at parse time and after the CLI overrides; record FASTQ paths only when the files are kept.
- **Moves numbers:** no.
- **Confidence:** confirmed.

**7. The manifest omits the index the transcriptome comes from**
- **Location:** `manifest.py:88–93`
  ```
      for key in ("genome", "gtf", "shadow_gtf", "strand_specificities"):
          if hasattr(config, key):
              value = getattr(config, key)
              if key in {"genome", "gtf", "shadow_gtf"} and value:
                  value = str(Path(value).resolve())
              manifest[key] = _jsonable(value)
  ```
- **Category / severity:** contract · low
- **What is wrong:** the simulated transcriptome is the index's (`WholeGenomeSimConfig.index`), yet the manifest records the GTF instead. Nothing ties `cfg.gtf` to the index: `run_simulation` requires the GTF to exist even when `index:` is set (`:1047–1048`).
- **Evidence:** instruments pass `--index` again and re-check that the rows line up by hand (`ruler_vs_truth.py:185–189`).
- **Change:** record the resolved index path and its source digests; refuse a mismatched `gtf:`.
- **Moves numbers:** no.
- **Confidence:** confirmed.

**8. `GDNAConfig.strand_kappa` is a legacy setting for a consumer that no longer exists**
- **Location:** `reads.py:79`, `:92`
  ```
      strand_kappa: float | None = None
              self.strand_kappa = None if od == 0.0 else (1.0 - od) / od
  ```
- **Category / severity:** dead · low
- **What is wrong:** od is converted to kappa here and back to od in `scenario.py:82–88`. Nothing else reads kappa.
- **Evidence:** the test docstring states the reason it exists: "legacy strand_kappa … so locus_sweep's strand_kappa still works" (`tests/test_sim_gdna_overdispersion.py:34–36`). `rg locus_sweep` finds only that docstring.
- **Change:** delete the field, the round trip, and the kappa asserts in that test (`:19–24`, `:42`).
- **Moves numbers:** no.
- **Confidence:** confirmed.

**9. Error injection: a "substitution" can redraw the original base, and nothing exercises the path**
- **Location:** `wgs_engine.py:461–465`
  ```
          mask = self._rng.random((count, read_len)) < self.sim_params.error_rate
          n_errors = int(mask.sum())
          if n_errors > 0:
              bases = np.array([65, 67, 71, 84], dtype=np.uint8)  # A, C, G, T
              seqs[mask] = bases[self._rng.integers(4, size=n_errors)]
  ```
- **Category / severity:** defect / test-gap · low (latent)
- **What is wrong:** the replacement base is drawn uniformly over ACGT, so a quarter of draws leave the base unchanged. The realised substitution rate is ¾ of the `error_rate` that `ReadSimConfig` documents.
- **Evidence:** every panel config and test sets `error_rate` to 0.0, and no test calls `_introduce_errors_batch`.
- **Change:** draw from the three other bases and add a test on the realised rate, or delete the setting.
- **Moves numbers:** no.
- **Confidence:** confirmed.

**10. Leftover from the `line → boundary` rename (ffa05dd2), and a replace that rewrites the wrong part of a condition name**
- **Location:** `benchmark.py:165`, `annotation.py:283`, `manifest.py:39`
  ```
          boundaries = [
              f"({sum(len(t.exons) for t in transcripts)} exon boundaries) → {gtf_path}"
          name = f"{name}_od_{gdna_strand_overdispersion:g}".replace(".", "p")
  ```
- **Category / severity:** naming · low
- **What is wrong:**
  - The rename turned two text-line senses into "boundaries": the summary's lines and the count of exon records in a log line.
  - Separately, `.replace(".", "p")` applies to the whole name prefix, not just the od value. An unlabelled rate 0.5 gives `gdna_r0.5_…` at od 0 but `gdna_r0p5_od_0p1_ss_1.00_…` when od > 0.
- **Evidence:** `git show ffa05dd2 -- src/rigel/sim/`.
- **Change:** go back to `lines` and "exon lines"; apply the replace to the od value only. Seeds are keyed at od 0, so they do not change.
- **Moves numbers:** no.
- **Confidence:** confirmed.

**11. Comments that describe old behaviour**
- **Representative quote:** `wgs_config.py:78`
  ```
      Three modes are supported, and ``mode`` selects exactly one; a field belonging to another mode
  ```
- **Category / severity:** stale-comment · low · confirmed. Each item below was checked against the code:
  - Retired nascent modes:
    - `wgs_config.py:78`: only two modes remain.
    - `wgs_config.py:105`: still says "shares".
    - `whole_genome.py:785` and `:888–889`: still say "the ratio modes", plural.
  - `whole_genome.py` docstrings:
    - `:6–7` says it writes `truth_abundances.tsv`; the file actually written is `truth_abundances_nrna_<label>.tsv`. `manifest.py:69` defaults to the name that is never written.
    - `:551–553` says `run_simulation` applies the nascent ratio; it is the orchestrator.
    - `:1032` says "Load genome + GTF"; it loads the index.
    - `:609` mentions "sim.py".
  - `scenario.py` and `__init__.py`:
    - `scenario.py:332–337`, `:446–450` and `__init__.py:15` put the index build last; it runs before the reads.
    - `scenario.py:336` says "Collate BAM"; the code runs `sort -n`.
    - `scenario.py:116` and `:443` say the oracle BAM is name-sorted; it is written `SO:unsorted GO:query` (`wgs_engine.py:366`).
    - `scenario.py:177` says only gDNA reads are excluded; nascent reads are excluded too.
  - Other:
    - `truth.py:144` and `whole_genome.py:939` end mid-sentence: "as start and end are".
    - `sampler.py:718` names a `project_genomic_block(s)` variant that does not exist.
    - `config.py:766–767` documents `ceil(needed × ratio)`, but `sampling.py:37` floors with `int()`.
    - `benchmark.py:126–128` says synthetic counts are "not tracked via n_nrna_pipeline". They are: `n_synthetic_observed` (`:266–268`) recomputes exactly `estimator.nrna_em_count`.
- **Moves numbers:** no.

**12. Inline I/O and format literals**
- **Location:** `wgs_engine.py:127–131`
  ```
  def _open_gzip_write(path: Path, threads: int = 4):
      """Open a gzip file for text writing, using pgzip if available."""
      if pgzip is not None:
          return pgzip.open(path, "wt", thread=threads, compresslevel=3)
      return gzip.open(path, "wt", compresslevel=3)
  ```
- **Category / severity:** bare-constant · low
- **What is wrong:**
  - The thread count and compression level sit inline, although the file already reads `CONSTANTS.resources` for similar settings.
  - The quality character `"I"` is written twice (`wgs_engine.py:408`, `bam.py:118`).
  - The `@PG` `VN: "1.0"` (`:375`) is not rigel's version.
  - `pgzip` is not declared in `pyproject.toml`, so the compressed bytes depend on the environment.
- **Change:** move these to `CONSTANTS.resources` and a single quality constant; stamp `rigel.__version__`.
- **Moves numbers:** no.
- **Confidence:** confirmed.

**13. Dead branches and tooling living in the package**
- **Representative quote:** `whole_genome.py:58–61`
  ```
  try:
      import yaml
  except ImportError:
      yaml = None  # type: ignore[assignment]
  ```
- **Category / severity:** dead / out-of-scope · low · confirmed
  - PyYAML is a declared dependency (`pyproject.toml:29`), so this guard and the one at `:244–245` cannot fire.
  - `whole_genome.py:965–966` (`lo == hi == 0.0`) is unreachable, because both `:343–344` and `:809–812` refuse `lo ≤ 0`.
  - `sampler.py:489`, `:535–537`: the BED12 strand is parsed and then deleted.
  - `wgs_engine.py:263`, `:275` add a second seed beside `sim_seed`. Only `tests/test_sim_genomic_refs.py:94` passes it, with the same value.
  - `whole_genome.py:1153–1252` is an argparse CLI reached only through `scripts/sim/simulate_reads.py`. `whole_genome.write_manifest` (`:1012–1021`) only wraps `manifest.write_manifest`.
- **Change:** delete the dead code; move `main()` into the script.
- **Moves numbers:** no.

## Already recorded
- `hygiene-ledger` (c): `rigel.sim.benchmark` and `Scenario.build_oracle` are test tooling inside the package.
- `hygiene-ledger` (f): the sharded writers and the whole-genome grid are coverage gaps.
- `hygiene-ledger` (g): `splice_motif.py`'s over-long lines.
- `hygiene-ledger` (h): nine `ground_truth_*` methods and three file passes (`scenario.py:151–228`).
- `hygiene-ledger` (i): the `rigel sim` CLI restates `Scenario`'s defaults.
- `latent-defects`:
  - `build_oracle`'s fixed-total branch ignores `gdna_fraction`, and so the reorder gate's fixture has no gDNA. The entry does not cover the reverse case: the fixed-RNA branch (`scenario.py:498`) drops `gdna_fraction` when no `GDNAConfig` is given. Every fixed-RNA fixture therefore passes `GDNAConfig(abundance=0.0)` purely to switch the branch on, and its `abundance` is never read.
  - The `.astype(int)` truncation in the fragment-length sampler (`sampling.py:40`). It also makes the claim at `whole_genome.py:749` that `fl_pmf` "matches … in distribution" untrue.
  - A capture entry with neither `enabled:` nor probes gives an empty `CaptureConfig` (`whole_genome.py:177–181`).
- ISSUES' WON'T-FIX table row "simulator captures a pre-mRNA through every probe its genomic blocks span" (nascent scope ruling).
- From your brief: gDNA drawn over the genome's N runs (`wgs_engine.py:861`).
- Context only, not records of my findings: `flgap-panels-stale-nascent-model` (CLOSED; claims `shares` was deleted, see 2) and `overlapping-synthetic-shadows` (the index's sense of "shadow", see 5).
