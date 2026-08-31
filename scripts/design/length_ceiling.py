#!/usr/bin/env python
"""⛔⛔ MEASURE THE CEILING BEFORE BUILDING THE CORRECTION — what is a PERFECT length model worth here?

The RNA length model is 8.4 bp low on every condition of this panel and the cause is localised. Before
fixing it, the project's own rule says price it: hand the solver the simulator's exact length
distributions and re-score.

⭐ **The arm already exists.** `pass0_vs_oracle.measure_condition` builds `c_input_pass0` — the same
prior-free solve with the simulator's own post-capture pmfs substituted for the fitted ones — on every
run, and the previous session's ladder work was throwing it away.

Two yardsticks, because they answer different questions and the project has been burned by quoting one
for the other:

* ``mwae_all`` — `pass0_vs_oracle`'s own mass-weighted error over every object with mass. Comparable
  with `ROADMAP.md`'s Stage-A line (what perfecting BOTH length models is worth) (−0.7 % on the pilot). ⛔ It
  scores honest ignorance as error, so it is not how pass-0 should be judged.
* ``mwae_solvable`` — the same error restricted to objects that HAVE own evidence, with the class
  partition **held fixed at the fitted-pmf run's** so the two arms are scored on the same population
  (`solvability_audit`'s own rule: a class that moves between arms cannot compare them).

⭐⭐ **AND IT PRICES EACH pmf SEPARATELY, WHICH IS THE WHOLE POINT** (TRAPS: price-the-halves-separately). The both-exact arm
alone hid a 14x split: measured on the ladder, the shipped solve's **-6.31 %** is **-5.90 % gDNA and
-0.43 % RNA**, all of it capture-ON (capture-OFF: -0.0000), while at pass-0 the RNA model is worth
**-0.02 %** on the solvable set and **+0.21 % (worse)** on all objects. A ceiling on an INPUT is not a
ceiling on its components.

⚠ Only the error is swapped, never the declared variance — `audit` reads `var_g` from the pass-0
capture, so a "confidently wrong" count across arms would mix one arm's answer with the other's
precision. This reports the error only, which is clean.
"""

from __future__ import annotations

import argparse
import copy
import importlib.util
import json
import sys
import time
from pathlib import Path

import numpy as np

DESIGN = Path(__file__).resolve().parent
sys.path.insert(0, str(DESIGN))


def _sib(name):
    key = name[:-3]
    if key not in sys.modules:
        sp = importlib.util.spec_from_file_location(key, DESIGN / name)
        m = importlib.util.module_from_spec(sp)
        sys.modules[key] = m
        sp.loader.exec_module(m)
    return sys.modules[key]


SA = _sib("solvability_audit.py")
P0 = _sib("pass0_vs_oracle.py")

from rigel.config import CalibrationConfig, PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402


def solvable_mwae(m, arm_name, axis, config):
    """`solvability_audit`'s error, on the class partition the FITTED run defined."""
    mm = copy.copy(m)
    mm.arms = dict(m.arms)
    mm.arms["pass0"] = m.arms[arm_name]
    a = SA.audit(mm, axis=axis, config=config)
    det = a["determined"]
    mass = float(a["total"][det].sum())
    return float(np.abs(a["err"])[det].sum()) / max(mass, 1.0)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--conditions", nargs="*", default=None)
    ap.add_argument("--suite", type=Path, default=Path.home() / "Downloads/rigel_runs/suite/ladder")
    ap.add_argument("--index", type=Path, default=P0.DEFAULT_INDEX)
    ap.add_argument("--oracle-cache", type=Path, default=None)
    ap.add_argument("--work-dir", type=Path, default=Path("/tmp/rigel_length_ceiling"))
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    index = TranscriptIndex.load(str(args.index))
    config = CalibrationConfig()
    names = args.conditions or sorted(
        p.name for p in args.suite.iterdir() if (p / "sim_oracle.bam").is_file()
    )
    with args.out.open("w") as fh:
        for name in names:
            t0 = time.perf_counter()
            cond = args.suite / name
            m = P0.measure_condition(
                bam=str(cond / "sim_oracle.bam"), index=index, pipeline_config=PipelineConfig(),
                calibration_config=config, work_dir=args.work_dir / "rigel_pass0_oracle", tag=name,
                truth_pmfs=lambda size, d=cond: (
                    P0.truth_length_pmf(d, "gdna", size), P0.truth_length_pmf(d, "rna", size)
                ),
                oracle_cache=args.oracle_cache,
            )
            row = {"condition": name, "f_gdna": P0.truth_f_gdna(cond) or 0.0}
            # ⚠ A ZERO-gDNA condition's truth file has no `gdna` row, so `truth_length_pmf` returns
            # None and `measure_condition` builds no exact-input arm. That is correct — there is no
            # gDNA length distribution to be exact about — and those rows are excluded from every
            # aggregate anyway. Record the absence rather than crashing on it.
            row["has_ceiling"] = "c_input_pass0" in m.scores
            arms = ("pass0", "c_input_pass0") if row["has_ceiling"] else ("pass0",)
            for axis in ("region", "boundary"):
                for arm in arms:
                    row[f"mwae_all_{arm}_{axis}"] = float(m.scores[arm][axis]["ALL"].mwae)
                    row[f"mwae_solvable_{arm}_{axis}"] = solvable_mwae(m, arm, axis, config)
            row["mwae_all_final"] = float(m.scores["final"]["region"]["ALL"].mwae)
            if row["has_ceiling"]:
                row["mwae_all_c_input_final"] = float(m.scores["c_input_final"]["region"]["ALL"].mwae)
                # ⭐ SPLIT THE BOTH-EXACT ARM. `measure_condition` only builds both-at-once, and
                # "which of the two models earns it" is the whole question here: the gDNA model is
                # already exact off capture and +6 % under it, while the RNA model is 8.4 bp low
                # everywhere. One variable each, same payload, same everything else.
                size = int(m.payload.max_length)
                g_t = P0.truth_length_pmf(cond, "gdna", size)
                r_t = P0.truth_length_pmf(cond, "rna", size)
                from dataclasses import replace as _replace

                p0cfg = _replace(config, calib_refit_iters=0)
                for tag_, cfg_, gp, rp in (
                    ("rna_only_pass0", p0cfg, None, r_t),
                    ("gdna_only_pass0", p0cfg, g_t, None),
                    ("rna_only_final", config, None, r_t),
                    ("gdna_only_final", config, g_t, None),
                ):
                    arm = P0.calibrate_arm(
                        m.payload, m.calibrate_kwargs, cfg_, gdna_pmf=gp, rna_pmf=rp
                    )
                    s = P0.score_axis(
                        arm.mass_gdna_region, arm.mass_rna_region,
                        m.truth.mass_gdna_region, m.truth.mass_rna_region,
                    )
                    row[f"mwae_all_{tag_}"] = float(s.mwae)
                    if tag_.endswith("pass0"):
                        mm = copy.copy(m)
                        mm.arms = dict(m.arms)
                        mm.arms["pass0"] = arm
                        a = SA.audit(mm, axis="region", config=config)
                        det = a["determined"]
                        row[f"mwae_solvable_{tag_}"] = float(
                            np.abs(a["err"])[det].sum() / max(a["total"][det].sum(), 1.0)
                        )
            fh.write(json.dumps(row) + "\n")
            fh.flush()
            print(f"  {name} {time.perf_counter() - t0:.0f}s", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
