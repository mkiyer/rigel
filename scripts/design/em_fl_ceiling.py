#!/usr/bin/env python
"""⭐⭐⭐ **WHAT IS A PERFECT gDNA FRAGMENT-LENGTH pmf WORTH END TO END, THROUGH THE EM? — the one fl
question that stops at no earlier stage.**

⛔ **Every other fl instrument stops at `calibrate`.** `length_ceiling.py` and
`calibration_truth_ab.py --ceiling` both price the pmf as an OPPORTUNITY input and go no further, and
`quant_accuracy.py`'s arms all inject `LocusPriors` fields, never an fl pmf. But `pipeline.py` also turns
`fl_models.gdna_pmf` into a `FragmentLengthModel` and hands it to the fragment scorer, so a wrong gDNA
length model is applied as a **per-fragment length likelihood in exactly the channel that separates
origins**. That path had no instrument; this is it.

⚠ **This is the fl pmf inside the OPPORTUNITY and SCORING models — NOT the fragment-length COMPOSITION
channel, which is retired until after 0.8.0 and must not be proposed.** They share a word and nothing else.

⭐⭐ **READ `gdna_frac_est`, WHICH IS THE PRODUCT** — `cli.py`'s own `gdna_fraction`, intergenic included —
and read it against `gdna_frac_true` from the simulator's `origin_counts`. ⛔ **The TRANSCRIPT rows are not
the deliverable here and have been measured flipping sign between the two fl-gap arms**, so a one-arm
transcript reading is how this instrument would be quoted backwards.

⛔⛔ **RUN IT ON A PANEL WHERE THE TWO COMPONENTS' FRAGMENT LENGTHS DIFFER, AND RUN BOTH SIGNS PLUS THE
EQUAL-LENGTH CONTROL.** The ladder and the main test chromosome give both origins equal lengths by design,
which is a forcing function for calibration and a NULL for this question. An effect that also appears on
the equal-length arm is an artefact rather than the length gap.

Three gates, and they are the reason a number from here can be believed:

* the injection **counts its fires** and raises if it never ran (`TRAPS: an-ablation-that-never-ran`);
* ``noop_fl`` replaces the pmf **with itself** and must be byte-identical — this fired on its first run,
  catching a renormalisation that moved the pmf by 8.7e-19 (`TRAPS: perturb-every-gate`);
* ``base_reseed`` is run as the **noise floor**, because the shipped pipeline is not reproducible run to
  run and no `quant_accuracy` delta below that floor is attributable.

⛔ It re-implements no scorer: the arms are run and scored by `quant_accuracy.run_condition`
(`TRAPS: score-the-consumers-own-count`).

    python scripts/design/em_fl_ceiling.py --panel <scenarios dir> --index <index dir> \\
        --conditions <cond> [<cond> ...]
"""

from __future__ import annotations

import argparse
import dataclasses
import importlib.util
import sys
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


QA = _sib("quant_accuracy.py")
P0 = _sib("pass0_vs_oracle.py")

from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

#: the metrics printed per axis. ⭐ `count_abs_err` is in FRAGMENTS, the unit the noise floor is quoted in.
_TX_METRICS = ("count_abs_err", "fp_mass", "mard")


def install_gdna_pmf(pmf: np.ndarray | None):
    """Swap ``FLModels.gdna_pmf`` for ``pmf``. ``None`` replaces it with ITSELF — the noop falsification.

    ``pipeline`` imports ``build_fl_models`` function-locally in both call sites, so patching the module
    attribute reaches production without touching the pipeline.
    """
    import rigel.calibration.fl as FL

    inner = FL.build_fl_models
    fired = {"n": 0, "max_abs_delta": 0.0}

    def wrapper(*a, **kw):
        models = inner(*a, **kw)
        cur = np.asarray(models.gdna_pmf, dtype=np.float64)
        if pmf is None:
            # ⛔ The noop must be EXACT — no renormalisation. Normalising here moved the pmf by one
            # float64 round-trip and the inertness gate caught it, which is the gate's whole point.
            new = cur.copy()
        else:
            new = np.asarray(pmf, dtype=np.float64)[: cur.size].copy()
            s = float(new.sum())
            if s > 0.0:
                new /= s
        fired["n"] += 1
        fired["max_abs_delta"] = max(fired["max_abs_delta"], float(np.abs(new - cur).max()))
        return dataclasses.replace(models, gdna_pmf=new)

    FL.build_fl_models = wrapper

    def restore():
        FL.build_fl_models = inner

    return restore, fired


def run_arm(arm, panel, index, cond, pipeline_config, truth_pmf):
    """One arm's scored rows, keyed by axis. ``base``/``base_reseed`` run unpatched."""
    if arm in ("base", "base_reseed"):
        return {r["axis"]: r for r in QA.run_condition(arm, panel, index, cond, pipeline_config, None)}
    restore, fired = install_gdna_pmf(None if arm == "noop_fl" else truth_pmf)
    try:
        out = QA.run_condition("base", panel, index, cond, pipeline_config, None)
    finally:
        restore()
    if fired["n"] == 0:
        raise RuntimeError(f"{cond} [{arm}]: the fl injection never ran. This is not a measured zero.")
    if arm == "noop_fl" and fired["max_abs_delta"] != 0.0:
        raise RuntimeError(
            f"{cond} [noop_fl]: replacing the pmf with itself moved it by {fired['max_abs_delta']:.3e}. "
            "Replacing an array with itself must be inert."
        )
    if arm == "true_gdna_fl" and fired["max_abs_delta"] == 0.0:
        raise RuntimeError(
            f"{cond} [true_gdna_fl]: the substitution did not move the pmf by one ULP. An arm that "
            "cannot move the quantity it names has not measured it."
        )
    return {r["axis"]: r for r in out}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--panel", type=Path, required=True, help="the scenarios directory")
    ap.add_argument("--index", type=Path, required=True)
    ap.add_argument("--conditions", nargs="+", required=True)
    ap.add_argument("--max-size", type=int, default=1000)
    args = ap.parse_args()

    index = TranscriptIndex.load(str(args.index))
    pc = PipelineConfig()
    for cond in args.conditions:
        truth_pmf = P0.truth_length_pmf(args.panel / cond, "gdna", args.max_size)
        arms = ["base", "base_reseed", "noop_fl"]
        if truth_pmf is None:
            print(f"\n══ {cond}: zero gDNA — no gDNA length distribution to be exact about, skipped")
            continue
        arms.append("true_gdna_fl")
        rows = {a: run_arm(a, args.panel, index, cond, pc, truth_pmf) for a in arms}

        base, floor_arm, exact = rows["base"], rows["base_reseed"], rows["true_gdna_fl"]
        print(f"\n══ {cond}")
        print(f"   {'axis':11s} {'metric':16s} {'base':>13s} {'exact gDNA fl':>14s} "
              f"{'delta':>12s} {'reseed floor':>13s}")
        for axis in ("transcript", "gene"):
            for key in _TX_METRICS:
                if key not in base[axis]:
                    continue
                b, e = float(base[axis][key]), float(exact[axis][key])
                fl = abs(float(floor_arm[axis][key]) - b)
                mark = "" if abs(e - b) > fl else "  <- INSIDE the floor"
                print(f"   {axis:11s} {key:16s} {b:13.4f} {e:14.4f} {e - b:+12.4f} {fl:13.4f}{mark}")
        lib, lib_e, lib_f = base["library"], exact["library"], floor_arm["library"]
        t = float(lib.get("gdna_frac_true", float("nan")))
        b, e = float(lib["gdna_frac_est"]), float(lib_e["gdna_frac_est"])
        fl = abs(float(lib_f["gdna_frac_est"]) - b)
        removed = 100.0 * (abs(b - t) - abs(e - t)) / abs(b - t) if abs(b - t) > 0 else float("nan")
        print(f"   {'library':11s} {'gdna_frac_est':16s} {b:13.6f} {e:14.6f} {e - b:+12.6f} {fl:13.6f}")
        print(f"   ⭐ THE PRODUCT: true {t:.6f} — a perfect gDNA fl pmf removes {removed:+.1f} % of the "
              f"standing bias, at {abs(e - b) / fl if fl > 0 else float('inf'):.0f}x the noise floor")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
