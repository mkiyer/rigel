# NEXT SESSION — start here (2026-09-14, after W12, W13 and the lanes worklist L1/L2/L4)

The whole picture, the reasoning and the ordered plan are in `docs/dev/CALIBRATION_PERFORMANCE_PLAN.md`
(THE LANES WORKLIST table is the section before §C's neighbour "Not on the list"; §F is the port). The θ
quadrature's derivation, the tilt study and the encompassing-locus audit are `docs/dev/THETA_QUADRATURE.md`
§§1–14. This file is only how to begin.

## The prompt for the next session (paste as the first message)

> Start from `main` (clean, at `602d5ab8` or later). Read `CLAUDE.md`, then `docs/dev/NEXT_SESSION.md`,
> `docs/dev/CALIBRATION_PERFORMANCE_PLAN.md` (THE LANES WORKLIST and §F) and `docs/dev/THETA_QUADRATURE.md`
> §§12–14. Run `scripts/design/preflight.py` and the suite before touching anything; `CLAUDE.md`'s baseline
> line is the count to reproduce (3,434 passed / 4 xfail / 3,438 collected).
>
> THE STANDING RULINGS are unchanged: elegant, simple, efficient, clear, concise, maintainable code, judged
> on the oracle metric (`calibration_vs_oracle.py`, per stratum, BOTH zero controls), the panel, the suite,
> the encompassing-locus test and the shared-exon stress; a restructure is proven bit-identical; a change
> that moves numbers is judged with its magnitudes read first; one mechanism at a time; no magic numbers —
> every constant derived; a falsification test first, verified failing, then break the fix and watch each
> gate fire; each item its own commit, committed on my go; never `ruff format scripts/`; patch `calibrate`
> through `importlib.import_module`; any config value is a `--set` arm.
>
> TWO ITEMS, IN ORDER, BEFORE THE PORT:
>
> L3 — the g00 case, `ISSUES: deadband-gates-a-gdna-free-library`: a METICULOUSLY FOCUSED derivation,
> design and fix. A library whose fitted gDNA count is exactly zero — the modal real case — sends no
> strand-derived RNA level, because the deadband's `1/N_gdna` term switches its strand channel off; the
> term has no derivation (gDNA's strand mean is ½ by symmetry), but deleting it alone is REFUSED with its
> number (the unstranded zero control 499 → 21,484): it was killing the unstranded phantom by accident.
> Derive the floor that kills the phantom on its own merits and the RNA level read from a slot's belief
> rather than only its strand claim; prototype outside `src/`; A/B on both panels with both zero controls,
> `test_encompassing_locus.py` and `deep_stress.py` on a `g00` donor; then `src/`, the derivation to
> `EQUATIONS.md` §5.2b, the ruling to `DESIGN.md`, the xfail flipped, the issue closed with its numbers.
>
> L5 — the witnessed atom, `ISSUES: capture-on-strand-pure-ambig-undercall`, approved: the AMBIG tilt's
> hypothesis space {pure +, pure −, mixed} at equal reference weight, a delivered RNA level on a strand
> ruling the other strand's pure hypothesis out. The prototype is `tilt_atom.py` in the session scratchpad
> (arm `atom_w`); the numbers to reproduce are in `THETA_QUADRATURE.md` §12 (ladder stranded ON 471,202 →
> 427,069, every stratum better, `g00 ss.99 OFF` +28). Gates first (the smoke-test truths, the encompassing
> exon∩exon xfail flipping), then `src/`, then the metric, the census bands and both stresses; the ruling to
> `DESIGN.md` §6b.15, the derivation to `EQUATIONS.md`, the issue closed.
>
> Then the deep-library baseline (two back-to-back profiler pairs at 8 threads), fresh `port_identity_*`
> references and `sweep_replay.py capture` on the landed tree, and THE PORT (plan §F).

## Before anything

```bash
source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel
python scripts/design/preflight.py                 # ~2 s: can this session run?
python -m pytest tests/ -q                         # CLAUDE.md's baseline line is the count to reproduce
```

## Where things stand

Pushed to `origin/main` as four commits after the W11 handoff:

* `23a431a9` — W12 the θ quadrature, both steps: the nodes follow the strand term's peak (`_tilt_window`,
  `_TILT_NODES` = 24 derived), the lanes deliver a row's ingredients (`CubeRow`), `sweep_n_tilt` and every
  θ lattice deleted. W13 — `sweeps_MO_3021_step5` and `port_identity_*` — was done on this commit and its
  docs rode along with L1.
* `cfb18e70` — L1: the RNA lanes exist whenever their own coordinate does (not the factory's, not the gDNA
  lane's); the rung-0 identity gate retired.
* `9d1375bf` — L2: one RNA coordinate for both lanes, a junction's certified flux always a source.
* `602d5ab8` — L4: `tests/calibration/test_encompassing_locus.py` (the owner's locus at four regimes, the
  exon∩exon solve gate an xfail until L5) and the `encompassing` rung of the toy ladder; L3 refused as a
  one-liner and filed.

W13's captures and references (`sweeps_MO_3021_step5`, `port_identity_*`) describe `23a431a9`; L2 moved
numbers, so they must be re-taken on the landed tree before the port (the last step above).

## The session scratchpad (persists across sessions; cited from the θ note, nothing in the tree cites it)

`/private/tmp/claude-503/-Users-mkiyer-proj-rigel/25f7f3df-7c89-49bd-826e-4c2efd9184b7/scratchpad/w12/`:
`tilt_census.py` (where the tilt matters, per stratum — the owner's instrument; arms `reference` /
`atom` / `atom_w` / `atom_s` / `jeffreys` / `profile`), `deep_stress.py` (the shared-exon toy, spliced and
mono), `audit_encompass.py` (the owner's locus, slot by slot), `tilt_atom.py` (L5's prototype),
`tilt_measure.py` (the two refused measures), `theta_window.py` (W12's prototype), `quadrature_check.py`,
`oracle_summary.py`, and every result json under `oracle/`, `census/`, `measure/`, `atom/`, `study/`,
`dissect/`, `deep*/`. `study/crosstab.py` is the presence-against-witness cross-tab.

## L3 — what is known, so the derivation starts where the last session stopped

* The floor: `σ²_d = ¼(1/N_rna + od_r) + ¼(1/N_gdna + od_g)`; `disc = 4·max(0, (κ−½)² − σ²_d)`
  (`region_init.strand_discriminability`). Binary readers of `disc > 0`: `sweep.solve_chain`'s
  `strand_live` → `ChainView.strand_live` → the library's `split_live` (the lanes' witness column) and
  `region_init`'s strand precision `tau_lam` (hence `has_own_composition`, the landscape's training
  population, the own claims the RNA lanes read as levels).
* On the ladder's unstranded g00 row: κ̂ = 0.500298, od_g = od_r = 0 (Poisson simulator), N_rna ≈ 1.5M, so
  the RNA half of the floor at 1σ is ~3e−4 against |κ̂−½| = 3e−4: a coin toss. With `1/N_gdna` gone,
  `disc` is positive there and the zero control reads 21,484 false fragments (499 with the term).
* On real data od_r > 0 widens the floor, but no multiple of σ is a derivation. Candidates: κ̂'s posterior
  width from `fit_strand_balance` as the floor; a continuous use of `disc` in place of the binary gates.
* The levels: `rna_lanes` builds a single-strand exon's own level from `own[x]` (the strand claim) only;
  the gDNA lane builds a full node's level from its own profile through its total. An RNA level read from
  the slot's belief (`belief_fg` always exists) would flow on a gDNA-free library whatever the deadband says.
* Instruments: the ladder's four g00 rows (both zero controls), `test_encompassing_locus.py` with a g00
  donor variant, `deep_stress.py` and `audit_encompass.py` on `gdna_g00_ss_0.99_nrna_mid_capture_off`.

## Decisions on record

* Float64 for the whole of ψ; ONE solver; ONE λ lattice at a dimensionless step with a stated guarantee; no
  θ lattice anywhere (the tilt count is derived).
* The tilt measure stays the arcsine and the marginal exact (both flattenings REFUSED on the ladder).
* The tilt is not a message; presence per strand is what the atom adds and the lanes witness.
* `drain`, `row` and `face` stay; the arcsine coordinate stays refused; the vertex atom is parked.
* The message cache's on/off switch (plan §E) is the owner's; the refit count stays; the scan's thread
  split is the owner's; parallelism waits for the port.
