# Holdout / reset investigation — how to reproduce

⚠ **A concurrent agent was editing `src/rigel/calibration/bp_solver.py` during this session and reverted my
edit mid-run.** Two bench arms in `/tmp/pass0_oracle_bench.tsv` are therefore CONTAMINATED and must be
ignored: `hoamb_r0`, `hoexcl_r0`, `hobnd_r0`, `howb_r0`, `hoamb_r1`, `hoexcl_r1` (the `ho*` prefix). Every
number reported was re-measured against a **frozen source snapshot**; those arms carry the `S_` prefix.

## The frozen snapshot

    /tmp/ho_snap/rigel/                 copy of src/rigel at HEAD 0bf11a37 + the two flags below
    /tmp/ho_snap/rigel/*.so             symlinks to the compiled extensions in site-packages/rigel
    /tmp/ho_snap/sitecustomize.py       strips scikit-build's ScikitBuildRedirectingFinder from
                                        sys.meta_path, which otherwise beats PYTHONPATH
    /tmp/ho_snap/debug/                 copy of scripts/debug (paths rewritten)
    /tmp/ho_snap/sp/                    copy of the scratchpad scripts (paths rewritten)

Run anything against it with:

    OMP_NUM_THREADS=1 PYTHONPATH=/tmp/ho_snap python /tmp/ho_snap/sp/<script>.py

## The flags (both DEFAULT OFF, both verified bit-identical unset — 32/32 at refit=0 and refit=1)

* `RIGEL_HOLDOUT=amb|bnd|excl|ambany` (`bp_solver.node_sweep`) — withhold a node set's OWN EVIDENCE
  (`pg_own`/`pp_own`/`pn_own`/`tau_own`/`mg_own` → 0). The node still relays.
* `RIGEL_HOLDOUT_WB=1` (`bp_solver`) — withhold the excluded half's WRITE-BACK only.
* `RIGEL_NORESET=1` (`calibrate`) — carry the pass-0 belief into the Phase-2 re-solve instead of resetting.

Patch of the `bp_solver` half against HEAD: `scratchpad/ho_holdout.patch`.

## Scripts

| script | what it measures |
|---|---|
| `ho_taint.py` | task 1 — per-substrate-node share of composition precision traceable to an excluded node |
| `ho_resetcheck.py` | tasks 2+3 — the exact reset / write-back-inertness proofs |
| `ho_solvable.py` | task 5 — the solvability predicates, suite-wide (writes `/tmp/ho_solvable.npz`) |
| `ho_subz2.py` | task 4 — substrate vs excluded accuracy + `z2\|Q1`, from `subacc_dump.py` npz files |
| `ho_summary.py` | arm rollup of `/tmp/pass0_oracle_bench.tsv` by stratum |
| `ho_ident.py` | bit-identity between two arms |
