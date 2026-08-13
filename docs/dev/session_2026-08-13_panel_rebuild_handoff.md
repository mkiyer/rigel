# Handoff — rebuild the panels from scratch, smaller

    ⚠ **A DEV DOC.** Provisional; nothing may cite it. When something settles, MOVE it to its permanent
    home and delete it here in the same edit.

    ⭐ Written at `2490393b`. The tree is committed and green: **0 fail / 3,401 pass / 0 skip / 10 xfail**.

---

## 0. ⛔⛔ READ THIS FIRST — deleting the ladder invalidates the rename's gate

`scripts/design/rename_identity.py`'s frozen reference is pinned to **one condition**:

    condition      gdna_g00_ss_0.99_nrna_none_capture_off      (ladder)
    quant digest   ac0b82b4e71ae431…      8,750 rows
    reference      ~/Downloads/rigel_runs/arms/rename_identity_reference.json

⭐ **That is acceptable, and here is why it is not a loss.** Stages 1–8a are all committed, each
verified ✅ BIT-IDENTICAL against that reference *before* the panels were touched. The gate has already
discharged its claim; it is not carrying anything forward.

⛔ **But do not let it fail silently.** On the new panel, either:
* **re-freeze** — delete the reference and run `--freeze` on a new condition, if any rename work remains
  (only stage 9 does, and §1 says it is superseded); or
* **retire it** — delete the reference file and say so, so nobody reads a stale digest as live.

⚠ Whichever you pick, `--check` against a deleted BAM fails with a confusing path error rather than a
useful one. Decide before the first `rm`.

---

## 1. ⭐ WHERE THE RENAME STANDS — 8 of 9, and stage 9 is SUPERSEDED

| stage | | |
|---|---|---|
| 1–8a | donor/acceptor · node→REGION · edge/seam→BOUNDARY · line→BOUNDARY · cut→REGION_BOUND · junction→SJ · EDGE_KIND_SJ · docs+memory · on-disk artifacts | ✅ every one proven bit-identical |
| 9 | the caches | ⛔ **SUPERSEDED by the rebuild-from-scratch decision** |

⛔ **`rescan_panels.py`'s byte-identity gate is INAPPLICABLE to a from-scratch rebuild, and that is worth
stating rather than discovering.** Its whole method is *old cache vs fresh scan, every non-target bank
identical*. With the BAMs re-simulated there is no old cache and no baseline — nothing to compare. Do not
run it expecting a verdict; it will either refuse or report every bank as new.
⭐ What replaces it: the simulator's own gates (`simulator_gates.py`), `verify_toy_substrate.py`, and
`suite_resolves.py` — *can this panel resolve the axis you are changing?* — run **before** quoting any
number off the new panel.

---

## 2. WHAT IS ON DISK, AND WHAT MUST SURVIVE

| | size | |
|---|---:|---|
| `suite/ladder/` | **27 G** | 36 conditions — delete |
| `suite/pilot/` | **8.5 G** | 8 conditions — delete |
| `suite/flgap_long/` | 3.2 G | 4 conditions — delete |
| `suite/flgap_short/` | 2.9 G | 4 conditions — delete |
| ⛔ `suite/reference/` | 145 M | **KEEP** — `genome.fa` + `genes.gtf` are the source of everything |
| ⛔ `suite/rigel_index/` | 12 M | **KEEP or rebuild** — see below |
| `~/Downloads/rigel_runs/arms/` | — | arm outputs + the frozen reference; stale once the panels change |

⭐ **The index was rebuilt at `2490393b` and is already on the new vocabulary** (`regions.feather`,
`region_id`). Rebuilding it again is free and deterministic — **measured: 6/6 data files byte-identical**,
only `manifest.json`'s `build_flags` differs. The exact command:

    rigel index --fasta <suite>/reference/genome.fa --gtf <suite>/reference/genes.gtf \
                --no-mappability --collapse-duplicate-transcripts -o <suite>/rigel_index

⚠ Both flags are required and neither is the default: without `--no-mappability` it demands an
`--alignable-zarr`; without `--collapse-duplicate-transcripts` it refuses the GTF's 5 duplicate groups.

---

## 3. ⭐⭐ THE WORKFLOW — one command per stage, and `status` first

`scripts/sim/panel.py` is the loop, driven by one panel YAML:

    status → build → simulate → cache → score → report

⛔ **`status` FIRST.** Every stage is expensive and resumable, and it names the next one.
⛔ **`cache` builds BOTH caches** — the scan cache and the ORACLE cache. The oracle one is the
origin-split truth every scoring instrument reads; a panel without it is one every scorer refuses.
⚠ A cached oracle condition needs all FOUR parts (`gdna`/`mrna`/`nrna`/`_main`).

---

## 4. ⭐ CHOOSING THE SMALLER PANEL — what may shrink, and what may not

The ladder is `9 gDNA levels × 2 strandedness × 2 capture = 36` at 10 M fragments each, which is the 27 G.

⛔ **THREE THINGS MUST SURVIVE ANY SHRINK, each for a measured reason:**

1. **`g00` — the ZERO-gDNA control.** Owner-required on every experiment. Its truth is exactly 0, so
   every deviation is a false positive with nothing to cancel it. `ROADMAP.md` §4.1's graveyard is
   eleven mechanisms that looked good on their target and died here.
2. **unstranded × capture-ON — the BLIND stratum.** It carries ~73 % of panel error and is the only
   open problem. A panel without it cannot see the thing the tool is worst at.
3. **At least THREE gDNA levels.** `TRAPS: a-single-level-panel-cannot-see-a-constant` — a channel that
   returns a near-CONSTANT scores well on a single-level panel and is refuted the moment a second level
   exists. That trap cost a whole analysis.

⭐ **A defensible 12-condition panel**: `{g00, g50, g98} × {ss 0.50, ss 0.99} × {capture off, on}`.
Keeps both controls, both strata, and three levels spanning the range. **~9 G at 10 M fragments.**
⚠ Depth is the other lever and is cheaper than conditions: 10 M → 4 M cuts the footprint 2.5× and keeps
every axis. ⛔ But check `suite_resolves.py` at the chosen depth before trusting a number from it — a
thinner panel resolves less, and that is exactly what that instrument measures.

⚠ `flgap_short` / `flgap_long` are a *pair* and only work as one: the short one is the admissible
drained panel (1 leaked spliced gDNA record) and the long one the cross-check (8,641). Rebuild both or
neither.

---

## 5. WHAT THIS SESSION LEFT BEHIND, briefly

* ⭐ The per-transcript RNA prior: foundation PROVEN, weighting function BUILT AND REFUSED
  (`ROADMAP.md` §4, `TRAPS: an-upper-bound-is-not-an-estimate`). Next candidate is a SPARSITY
  mechanism, not a better mean.
* ⭐ `sj_mass[2]`, the two dead banks, and `CalibrationResult.junction_conserved_mass` — all landed.
* ⭐ `rescan_panels.py`'s gate was UNSATISFIABLE for three days and is repaired
  (`TRAPS: a-stale-gate-accuses-the-newest-change`). It still matters for the *next* schema change,
  just not for this rebuild.
* ⚠ `docs/dev/rename_plan.md` holds the ten lessons the rename produced. `rename_census.py` re-derives
  the landscape; it is deliberately exempt from every rename stage because it searches for the old
  tokens.
