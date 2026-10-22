---
name: ladder-report
description: Build and publish Rigel's report card for the 16-condition gDNA ladder — per-scenario pool error (gDNA / synthetic nascent RNA / annotated RNA), transcript- and gene-level error, each on the simulator's own alignments AND on STAR alignments of the same reads (the STAR rung), beside the tagged 0.7.1 release; the fragment-length gap panels, calibration against the oracle and the VCaP mix — as one Artifact page that updates in place, and a markdown file. Use after any change to Rigel that could move a number, when the owner asks how the tool is doing end to end, or when the release report needs refreshing.
---

# The ladder report card

This is how progress toward a production release is assessed. It renders files the instruments already wrote;
**it runs nothing and measures nothing.**

⛔⛔ **THE REPORT IS A RENDERING, NOT A MEASUREMENT.** If Rigel changed since its inputs were written, the card will
render stale numbers and look perfectly current. Re-score first — step 1 is not optional after a code change.

## 1. Bring every input up to date

Export `RIGEL_RUNS` (the data root, `RUNS` in `scripts/design/_shared.py`), activate the `rigel` env, set
`OMP_NUM_THREADS=1`, and preserve each panel's previous arms first (`cp -r arms arms_<what-it-was>_<date>`): a
before-and-after is what makes the next delta attributable. `panel.py score` passes
`--set em.assignment_mode=fractional` and `--set scan.total_threads=1` itself.

| input | command | after a change to Rigel? |
|---|---|---|
| the ladder on the simulator's alignments | `python scripts/sim/panel.py score --config scripts/sim/configs/gdna_ladder.yaml --jobs 8 --arms base oracle` | always |
| **the STAR rung's alignments** | `python scripts/sim/panel.py align --config scripts/sim/configs/gdna_ladder.yaml --aligner <rustar or STAR> --aligner-index <STAR index of the suite genome> --threads 16` | never: once per panel (about 50 s and 1 GB per condition); the BAMs are kept, the FASTQs are not (`align` rebuilds a pair from the simulator's BAM) |
| **the STAR rung** | `python scripts/sim/panel.py score --config scripts/sim/configs/gdna_ladder.yaml --jobs 8 --arms base --bam rustar_recommended` | always |
| the two gap panels | `python scripts/sim/panel.py score --config scripts/sim/configs/flgap_rna_long.yaml --jobs 4 --arms base`, and `flgap_rna_short.yaml` | always |
| calibration against the oracle | `python scripts/design/calibration_vs_oracle.py --suite $RIGEL_RUNS/suite/<panel> --index $RIGEL_RUNS/suite/rigel_index --jobs 8 --json <dir>/cvo_<panel>.json`, for `ladder`, `flgap_rna_long`, `flgap_rna_short` | always; it reads the cached scan, so re-cache after a scan change |
| the VCaP mix | `rigel quant --bam <the VCaP mix's BAM> --index $RIGEL_RUNS/refs/rigel_index -o <dir>/vcap_rc --threads 8`, alone (one whole-genome job at a time) | always |

`panel.py status` shows the STAR rung as `align rustar_recommended 16/16`. The rung's aligner is rustar (the Rust
port of STAR, which reads STAR's index and options and matches STAR on 99.9 % of templates); its options are
`scripts/sim/configs/star_recommended.params.txt`, the MANUAL's recommended STAR settings.

## 2. Build the page and the markdown

```bash
python .claude/skills/ladder-report/build_report.py \
    --html <scratchpad>/ladder_report.html \
    --cvo <dir> --vcap <dir>/vcap_rc/summary.json --vcap-truth 0.2518 --vcap-v071 0.263117 \
    --v071 $RIGEL_RUNS/prototypes/2026-10-05_capture_learning/baseline/release_tag.jsonl \
    --baseline "v0.7.1, scored 2026-10-05 on the panel before its re-simulation — a ballpark" \
    --markdown $RIGEL_RUNS/reports/ladder_accuracy_<DATE>.md
```

The 0.7.1 rows were scored on the panel before the re-simulation of 2026-10-21 and stay on the card as a ballpark
(owner, 2026-10-21): their reads drew gDNA over the `N` runs and carried no errors, so a delta against them crosses
simulations. Re-score the tagged release on the current panel before reading one.

`--vcap-truth` is the mix's DNA-origin share by read name, and `--vcap-v071` 0.7.1's reading of it; `--head` and
`--scored` default to HEAD and the base arm's date. The markdown comes from `quant_accuracy.py --markdown`, which the
suite gates (`tests/calibration/test_quant_accuracy.py`); the page is `report_template.html` with the payload in place
of `__DATA__`.

## 3. Publish the page — to the SAME artifact

**`https://claude.ai/artifact/Cek2wmKtitgbDyfM5gNqyj`**

Pass that as `url` to the Artifact tool so the owner's link keeps working and the history accumulates:

```
Artifact(action="publish", url="https://claude.ai/artifact/Cek2wmKtitgbDyfM5gNqyj",
         file_path="<scratchpad>/ladder_report.html")
```

⛔ Read the artifact first if this conversation has not published it. Publishing **without** `url`
creates a second report and the owner is left with two links and no way to tell which is current.

## What the report must never get wrong

The numbers underneath can be right while the report lies. Each of these is gated and each has been
watched to fail under its own perturbation — if you change the rendering, keep them true.

| the lie | the truth |
|---|---|
| gDNA scored EM-only | `gdna_est` **excludes** intergenic fragments; the truth counts them. The comparable estimate is `gdna_est + n_intergenic`, and off capture the intergenic half is the larger one |
| a percent error at a truth of zero | The `g00` rung has no gDNA at all. Print `n/a`; the raw count is the whole answer (`TRAPS: a-ratio-cannot-carry-zero`) |
| the two RNA pools folded together | The split is `is_synthetic`, **never** `is_single_exon`. A single-exon transcript is **annotated** RNA (`TRAPS: nrna-does-not-mean-synthetic`) |
| a false-negative column | `fn_mass` needs an estimate of *exactly* zero, which a fractional posterior never is — it reads 0 on all 16 conditions and looks like a perfect score for something unmeasured (`TRAPS: could-the-arm-have-fired`). Report `count_under` |
| the truth table's row count as a transcript count | The table carries a row per synthetic entity too (6,919 of 15,669 on the current index), zero on **both** sides. Only 2,466 of 9,385 gene rows are real genes. `expressed` and `detected` are the scored sets |
| the deferred stratum unmarked | Unstranded × capture-ON carries most of the error; unmarked, a reader takes a pooled total for the tool's accuracy (`TRAPS: never-pool-the-strata`) |
| the STAR rung's gDNA gap read as Rigel's | The aligner leaves a little of the simulated library unplaced: since the re-simulation of 2026-10-21, 1.4 % of the gDNA off capture and 0.07 % on (it was a fifth off capture, drawn over the suite genome's N runs). The rung is scored against the same truth as the simulator's alignments, so the two compare directly; `unaligned` counts what the aligner left out by origin, and the gDNA fraction's truth is read over what it placed (`gdna_frac_true_aligned`) |
| a rung that mixes aligners | Every condition comes from one `panel.py align` run: its `<name>.align.json` records the aligner, its version and the parameter file, and `quant_accuracy.py --bam` stamps the input on every row |
| a ragged markdown table | It renders as a *wrong* table, not an error. The escaped pipes in a `Σ\|Δ\|` header are cell content, not delimiters |

Two labels that must stay explicit: **MARD is symmetric** (`|est−true|/(|est|+|true|)`, bounded [0,1] —
0.5 is a 3× error, not 50 %), and the pool rows carry one irreducible **asymmetry** — the estimate
splits on the entity, the truth on the simulator's template kind, so a few annotated single-exon
transcripts contribute ~300 fragments to nascent truth while their tool-side counts land in annotated.

## How the report is read

- **Per stratum, never pooled.** Three strata are in scope; unstranded × capture-ON is deferred.
- **Two inputs, one truth.** The simulator's alignments are the regression check; the STAR rung is what a user's
  library looks like to Rigel, and the difference between the two is the alignment's.
- **By size.** An A/B pair runs with the scan pinned (`--set scan.total_threads=1`) and is exactly
  reproducible; judge an effect by its size, genes and pools beside the transcript table.
- **The capture-OFF magnitudes are a stress reading.** The panel runs 20.2 % nascent fragments against
  a realistic ~4.2 % (`DESIGN.md` §0b's nascent scope ruling), so quote the share with the number.

## Files

| | |
|---|---|
| `build_report.py` | reads every input above, derives the truth table's composition and the library depth, writes the page and (optionally) the markdown |
| `report_template.html` | the page: tokens for both themes, the headline tiles, the per-stratum and STAR-rung tables, per-scenario pool cards, the metric tables, the gap panels, calibration and VCaP. `__DATA__` is the payload placeholder |

The palette is the data-viz reference instance's categorical slots 1–3 in fixed order (gDNA blue,
nascent orange, annotated aqua), both themes. If you restyle, keep the assignment stable — colour
follows the pool, never its rank.
