# NEXT SESSION — the state after the 2026-08-27 tear-down

    ⚠ **A DEV DOC, and it is a HANDOFF.** It describes where things stand; it does not tell you
    what to build. MOVE anything that settles into the permanent docs and DELETE this file.

## What this session did: it removed code, and left one skeleton

A long message-policy campaign was **torn down**. Deleted: `CurrencyPolicy`, a unified bridge
(`FrameAwarePropagation`, `AllocationSolve` and a stack of mechanisms built on them), their test
files, and the config values that reached them. Git carries all of it.

What remains in the message layer, and it is deliberately small:

| | |
|---|---|
| `messages/foundation.py` | the ratified spec — one `Message` with provenance lanes, the propagate/solve timepoints, and the laws the skeleton ENFORCES (spliced lanes never relay; a single-source transform may only lower a precision; lanes never mix) |
| `messages/policy.py` | `MessagePolicy` — the runner that plugs a `(PropagationModel, SolveModel)` pair into the backbone. With trivial models it is **byte-identical to `SilentPolicy`**, gated and confirmed on the panel |
| `messages/silent.py` | `SilentPolicy` — the measured floor. **Frozen** |
| `messages/relay.py` | `RelayPolicy` — the shipped default. **Frozen** |

`CalibrationConfig.message_policy` selects between them (`"relay"` / `"silent"` / `"message"`); an
unknown name raises rather than silently falling back to the relay.

## The test chromosome is EMPTY, and it is the owner's to design

`scripts/sim/test_reference/` — the GTF, the abundances TSV and the probes BED were cleared to zero
transcripts on 2026-08-27. The owner adds transcripts one at a time; the FASTA, index, reads and
caches are all DERIVED and must be rebuilt after every edit.

⛔ The 42-transcript derived artifacts were moved to
`~/Downloads/rigel_runs/test_reference_STALE_42tx_2026-08-27/` rather than left in place, because a
benchmark run against stale caches answers a different question and says nothing about it.

**`docs/TESTING.md` §0a has the four commands** (build → simulate → cache → score) and the
constraints that bite: κ is fitted from spliced reads, so at least one multi-exon transcript with
real depth is needed before any number means anything.

## How the work is judged

⛔ **The bar is NOT "beat `SilentPolicy`."** Silence is the floor, and on strand-specific data a
sighted exon's own strand solve is excellent — a message can mostly only disturb it. Message
propagation exists for **unstranded data and AMBIG slots**, where the strand channel is dead.

* **unstranded rows** — the policy must **WIN**.
* **stranded rows** — the policy must do **as little harm as possible** (near 1.00× of silence).

⛔ Never pool the two halves. `python scripts/design/policy_benchmark.py --panel test` prints them
apart, in seconds; `--panel ladder` is the shipping judgement.

## Before believing anything

* `python scripts/design/preflight.py` — one command, one verdict, before anything else.
* `python -m pytest tests/ -q` — ANY failure is a regression; the count lives in `CLAUDE.md`.
* A claim must name its SUBSTRATE (`TRAPS: a-toy-and-a-panel-can-disagree-in-rank`).
* `HONEST_PRECISION.md` in this directory is the record of what the deleted campaign tried,
  measured and REFUTED. Read it before re-proposing a mechanism, so a dead end is not re-run.
