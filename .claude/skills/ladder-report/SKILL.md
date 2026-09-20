---
name: ladder-report
description: Build and publish Rigel's end-to-end accuracy report for the 16-condition gDNA ladder — per-scenario pool error (gDNA / synthetic nascent RNA / annotated RNA), then transcript- and gene-level error inside the annotated pool — as a markdown file and as an Artifact page that updates in place. Use after any change to Rigel that could move a number, when the owner asks how the tool is doing end to end, or when the release report needs refreshing.
---

# The ladder accuracy report

This is how progress toward a production release is assessed. It renders the arm files
`quant_accuracy.py` writes; **it runs nothing and measures nothing.**

⛔⛔ **THE REPORT IS A RENDERING, NOT A MEASUREMENT.** If Rigel changed since the arms were scored, the
report will render stale numbers and look perfectly current. Re-score first — step 1 is not optional
after a code change.

## 1. Are the arms current?

The arms live in `~/Downloads/rigel_runs/suite/ladder/arms/`. They are current only if they were
written *after* the last change that could move a number.

```bash
ls -la ~/Downloads/rigel_runs/suite/ladder/arms/*.jsonl
git log -1 --format='%h %ad %s' --date=iso
```

If the tool has moved, re-score — ~25 minutes for four arms, and it rebuilds nothing:

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
export OMP_NUM_THREADS=1
python scripts/sim/panel.py score --config scripts/sim/configs/gdna_ladder.yaml --jobs 8 \
    --arms base base_reseed oracle oracle_ruler
```

⚠ `panel.py score` passes `--set em.assignment_mode=fractional` itself. Preserve the previous arms
first (`cp -r arms arms_<what-it-was>_<date>`) — a before-and-after is what makes the next delta
attributable, and the report prints the seed floor from `base_reseed`.

## 2. Build both renderings

```bash
python .claude/skills/ladder-report/build_report.py \
    --html /tmp/ladder_report.html \
    --markdown ~/Downloads/rigel_runs/reports/ladder_accuracy_<DATE>.md
```

The markdown comes from `quant_accuracy.py --markdown`, which the suite gates
(`tests/calibration/test_quant_accuracy.py`); the page is the same payload through
`report_template.html`. The builder imports the instrument rather than re-implementing its label
rules, so the two cannot disagree about what a field means.

## 3. Publish the page — to the SAME artifact

**`https://claude.ai/artifact/Cek2wmKtitgbDyfM5gNqyj`**

Pass that as `url` to the Artifact tool so the owner's link keeps working and the history accumulates:

```
Artifact(action="publish", url="https://claude.ai/artifact/Cek2wmKtitgbDyfM5gNqyj",
         file_path="/tmp/ladder_report.html")
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
| the two RNA pools folded together | The split is `is_synthetic`, **never** `is_nrna`. A single-exon **annotated** transcript carries `is_nrna` and is annotated RNA (`TRAPS: nrna-does-not-mean-synthetic`) |
| a false-negative column | `fn_mass` needs an estimate of *exactly* zero, which a fractional posterior never is — it reads 0 on all 16 conditions and looks like a perfect score for something unmeasured (`TRAPS: could-the-arm-have-fired`). Report `count_under` |
| the truth table's row count as a transcript count | The table carries a row per synthetic entity too (6,919 of 15,669 on the current index), zero on **both** sides. Only 2,466 of 9,385 gene rows are real genes. `expressed` and `detected` are the scored sets |
| the deferred stratum unmarked | Unstranded × capture-ON carries most of the error; unmarked, a reader takes a pooled total for the tool's accuracy (`TRAPS: never-pool-the-strata`) |
| a ragged markdown table | It renders as a *wrong* table, not an error. The escaped pipes in a `Σ\|Δ\|` header are cell content, not delimiters |

Two labels that must stay explicit: **MARD is symmetric** (`|est−true|/(|est|+|true|)`, bounded [0,1] —
0.5 is a 3× error, not 50 %), and the pool rows carry one irreducible **asymmetry** — the estimate
splits on the entity, the truth on the simulator's template kind, so a few annotated single-exon
transcripts contribute ~300 fragments to nascent truth while their tool-side counts land in annotated.

## How the report is read

- **Per stratum, never pooled.** Three strata are in scope; unstranded × capture-ON is deferred.
- **Above the floor.** Re-running the identical command moves these figures by about the reseed
  floor's own size, and it is not the seed (`TRAPS: the-deliverable-is-not-reproducible-by-default`).
  Never act on a smaller difference.
- **The capture-OFF magnitudes are a stress reading.** The panel runs 20.2 % nascent fragments against
  a realistic ~4.2 % (`DESIGN.md` §0b's nascent scope ruling), so quote the share with the number.

## Files

| | |
|---|---|
| `build_report.py` | reads the arm jsonl, derives the truth table's composition and the library depth, writes the page and (optionally) the markdown |
| `report_template.html` | the page: tokens for both themes, per-scenario cards with paired truth/estimate bars, the metric tables. `__DATA__` is the payload placeholder |

The palette is the data-viz reference instance's categorical slots 1–3 in fixed order (gDNA blue,
nascent orange, annotated aqua), both themes. If you restyle, keep the assignment stable — colour
follows the pool, never its rank.
