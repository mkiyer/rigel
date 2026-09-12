# NEXT SESSION — start here (2026-09-11, after commit `38755cb1`)

The whole picture, the reasoning and the ordered plan are in `docs/dev/CALIBRATION_PERFORMANCE_PLAN.md`
— read it after `CLAUDE.md`. This file is only how to begin.

## Before anything

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
python scripts/design/preflight.py                 # ~2 s: can this session run?
python -m pytest tests/ -q                         # CLAUDE.md's baseline line is the count to reproduce
```

The bit-identity baseline is the chunk-exact tree: `~/Downloads/rigel_runs/perf/sweeps_MO_3021_step2`
(four captured sweeps; `sweep_replay.py replay --dir … --call 0..3`, `--block-slots N|none`) and the
three references `~/Downloads/rigel_runs/arms/locus_identity_*.json` (`rename_identity.py --check`; the
real BAMs are under `~/Downloads/rigel_runs/cfrna/mctp_<lib>_*/bam/star.srt.rmdup.collate.bam`). The
older `sweeps_MO_3021`, `cleanup_identity_*` and `perf_identity_LBX0190.json` differ from the tree by the
priced ≤ 1e-14 and are superseded.

## The first two tasks, in order (plan §3 A and B)

1. **Names and cleanup** — DONE 2026-09-11/12 (plan §3 A and §5): the names ruled and applied, `sweep.py`
   split into backbone / `blocks` / `message_cache`, the lanes and the face table out of `transfer.py`, the
   cache key over the context's fields, the small things; the tiling folded into D. Every step
   bit-identical on the replay, the three references and the suite.
2. **The received messages as tables** — the design is written in plan §3 B with its names; the owner
   agreed to it and deferred the implementation to the session's judgement: clear, concise, efficient.
   Falsification tests first, verified failing; then the code; then break it and watch the gates fire.

Do not start the port. Steps C (the tolerance gate) and D (memory) come before it, and the port begins only
when the owner is satisfied with the Python design and implementation.

## Decisions on record

* The refit count stays; the scan's thread split is the owner's (`ISSUES: scan-thread-split-starves-the-workers`).
* The thread count reuses `--threads`; the block size is a measured tunable (`sweep_block_slots`).
* ψ's chunk-exact read-out moved numbers by ≤ 1e-14: accepted as identical to a tolerance.
* Parallelism waits for the port; threads are refuted for the Python passes.
* The shared message layer's on/off switch is an open ruling (plan §3 E).
* Patch `calibrate` through `importlib.import_module("rigel.calibration.calibrate")`, never `import … as`.
