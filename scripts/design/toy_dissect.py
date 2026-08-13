#!/usr/bin/env python
"""⭐⭐ ONE SCENARIO, DISSECTED TO EVERY SLOT AND EVERY CHANNEL — where does the invented gDNA come from?

The survey says the two-exon transcript is solved to <1 % on 26 of 36 scenarios and OVER-CALLS gDNA on
the rest. This opens one scenario and follows the number through the three rungs the solver itself
computes, which are already in ``_debug['capture']`` and need no re-solve:

    fg_strand   the strand likelihood ALONE — no prior, no factory, no messages
    fg_loc      the message-free SELF-SOLVE — strand + intron factory + psi's own reference
    f_g         the FINAL answer, after the forward-backward relay

⭐ Reading the ladder localises the defect without any guessing: if ``fg_loc`` is already wrong the
fault is in the per-object initialisation, and if ``fg_loc`` is right and ``f_g`` is not, it is the
messages. Then the message CHANNELS are printed beside it — ``cm_g`` is the gDNA measurement precision
psi receives, ``c_tau`` the composition precision, ``cg`` the fused gDNA density — so a dead channel is
visible as a zero rather than inferred.

⭐⭐ **WHAT IT FOUND, 2026-08-05, on `gdna_g00_ss_0.50_nrna_none_capture_off`** (zero gDNA, so truth is
exactly 0) — the chain, measured:

* the exon has NO own evidence (``tau = 0``, ``fg_loc = fg_strand = 0.4902`` — psi's reference);
* every gDNA channel into it is DEAD (``cg``, ``cm_g``, ``c_tau`` all 0);
* so the only channel left is the certified-RNA measurement grafted from the junction, and taken at
  face value it is RIGHT: it implies ``f_g = -0.003``;
* ⛔ but its precision is capped by ``spl_prec = J/(1 + J·(log r)²)`` with
  ``r = rho_tot(exon)/rho_tot(intron|exon EDGE)``. Measured r = 0.914 where it should be 1, so 15,639
  junction reads deliver **122** effective observations and the phantom gDNA is FLAT at ~2 % over a
  64x depth range. Handing the solver the toy's own realised RNA length pmf moves r to 0.964,
  ``spl_prec`` to 714, and halves the error.
"""

from __future__ import annotations

import argparse
import dataclasses
import importlib.util
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402

DESIGN = Path(__file__).resolve().parent
_s = importlib.util.spec_from_file_location("toy_ceiling", DESIGN / "toy_ceiling.py")
TC = importlib.util.module_from_spec(_s)
sys.modules["toy_ceiling"] = TC
_s.loader.exec_module(TC)
TH = TC.TH

from rigel.calibration.calibrate import calibrate  # noqa: E402
from rigel.calibration.region_chain import REGION  # noqa: E402
from rigel.calibration.signature import coarse_type_array  # noqa: E402
from rigel.config import CalibrationConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402
from rigel.scan_cache import index_derived_inputs  # noqa: E402

SUITE = Path.home() / "Downloads/rigel_runs/suite/ladder"
INDEX = Path.home() / "Downloads/rigel_runs/suite/rigel_index"
TYPES = {0: "intergenic", 1: "intron", 2: "exon"}


def dissect(cond: str, *, n_rna: int, genome_length: int, work_dir: Path):
    index = TranscriptIndex.load(str(INDEX))
    cfg = dataclasses.replace(CalibrationConfig(), calib_refit_iters=0)
    donor = TH.harvest(SUITE / cond, index, config=cfg)
    spec = TH.SPECS["spliced_exons"]
    if genome_length:
        spec = dataclasses.replace(spec, genome_length=genome_length)
    spec = dataclasses.replace(spec, n_rna_fragments=n_rna, name=f"dis_{cond[5:12]}")
    sub = TC.simulate(spec, donor, work_dir)
    dbg: dict = {}
    calibrate(payload=sub.payload, strand_model=sub.strand_model, gdna_fl_pmf=donor.gdna_fl_pmf,
              rna_fl_pmf=donor.rna_fl_pmf, config=cfg, injected_priors=donor.priors, _debug=dbg,
              **index_derived_inputs(sub.index))
    cap, chain = dbg["capture"], dbg["chain"]
    uni = cap["_uni"][-1]
    st = cap["_uni_static"]
    ra = sub.region_arrays
    rtype = coarse_type_array(np.asarray(ra.signature)).astype(np.int64)
    starts, sizes = np.asarray(ra.start, np.int64), np.asarray(ra.region_size_bp, np.int64)
    kind, obj = np.asarray(chain.kind), np.asarray(chain.obj_idx, np.int64)

    ov = sub.truth.override_masses(ra)
    tg = {"region": np.asarray(ov["mass_gdna_region"], float),
          "edge": np.asarray(ov["mass_gdna_edge"], float)}
    tr = {"region": np.asarray(ov["mass_rna_region"], float),
          "edge": np.asarray(ov["mass_rna_edge"], float)}

    print("=" * 132)
    print(f"⭐⭐ DISSECTION — {cond}")
    print(f"   toy: {spec.genome_length:,} bp · {n_rna:,} RNA fragments · "
          f"{sub.n_gdna_target:,} gDNA fragments (donor rate {donor.gdna_rate_per_base:.6g}/bp)")
    print(f"   kappa = {donor.priors.rna_sense_frac:.6f}   strand specificity "
          f"{donor.strand_specificity:g}   capture {'ON' if donor.capture_on else 'off'}")
    print("=" * 132)
    print("\n── THE SOLVER'S OWN LADDER, per slot.  strand-only → self-solve → final ───────────────")
    print(f"\n   {'slot':<26} {'n':>7} {'TRUTH':>7} {'fg_strand':>10} {'fg_loc':>8} {'f_g':>8} "
          f"{'Δ final':>8} {'tau_own':>9} {'sd(f_g)':>8}")
    print("   " + "-" * 118)
    rows = []
    for s in range(int(chain.n_slots)):
        i = int(obj[s])
        ax = "region" if kind[s] == REGION else "edge"
        if ax == "region":
            label = f"{TYPES[int(rtype[i])]} [{starts[i]:,},{starts[i] + sizes[i]:,})"
        else:
            hi = s + 1
            b = int(rtype[obj[hi]]) if hi < int(chain.n_slots) and kind[hi] == REGION else -1
            a = int(rtype[obj[s - 1]]) if s > 0 and kind[s - 1] == REGION else -1
            label = "|".join(TYPES.get(x, "?") for x in sorted((a, b)) if x >= 0) + \
                    (f" @{starts[obj[hi]]:,}" if b >= 0 else "")
        tot = tg[ax][i] + tr[ax][i]
        truth = tg[ax][i] / tot if tot > 0 else float("nan")
        n = float(np.asarray(cap["count"], float).sum(axis=1)[s])
        rows.append((s, label, n, truth))
        print(f"   {label:<26} {n:>7.0f} {truth:>7.4f} {cap['fg_strand'][s]:>10.4f} "
              f"{cap['fg_loc'][s]:>8.4f} {cap['f_g'][s]:>8.4f} "
              f"{cap['f_g'][s] - truth:>+8.4f} {cap['_tau0_lam'][s]:>9.3g} "
              f"{np.sqrt(max(cap['var_g'][s], 0)):>8.3f}")

    print("\n── THE MESSAGE CHANNELS INTO EACH SLOT ────────────────────────────────────────────────")
    print("   cm_g  = the gDNA MEASUREMENT precision psi receives (0 ⇒ the level arrives with no weight)")
    print("   c_tau = the COMPOSITION precision (the single-λ message)")
    print("   cg    = the fused gDNA DENSITY the messages deliver;  own = this slot's own self-solve")
    print(f"\n   {'slot':<26} {'own rho_g':>10} {'msg cg':>10} {'cm_g':>9} {'c_tau':>9} "
          f"{'lam_msg':>9} {'M/E_g':>9} {'E_g':>8} {'E_r':>8}")
    print("   " + "-" * 118)
    for s, label, n, truth in rows:
        M, Eg, Er = st["M"][s], st["E_g"][s], st["E_r"][s]
        print(f"   {label:<26} {st['og'][s]:>10.4g} {uni['cg'][s]:>10.4g} {uni['cm_g'][s]:>9.3g} "
              f"{uni['c_tau'][s]:>9.3g} {uni['lam_msg'][s]:>9.3g} "
              f"{(M / Eg if Eg > 0 else float('nan')):>9.4g} {Eg:>8.1f} {Er:>8.1f}")

    print("\n── WHERE THE gDNA LEVEL AT THE EXONS COMES FROM ───────────────────────────────────────")
    for s, label, n, truth in rows:
        if not label.startswith("exon"):
            continue
        M, Eg = st["M"][s], st["E_g"][s]
        print(f"\n   {label}   n = {n:,.0f}   truth f_g = {truth:.4f}   final f_g = "
              f"{cap['f_g'][s]:.4f}")
        print(f"      the slot's own total density M/E_g            {M / Eg:>12.5g}")
        print(f"      the gDNA density the messages delivered (cg)  {uni['cg'][s]:>12.5g}")
        print(f"      implied gDNA fraction  cg·E_g/M               "
              f"{uni['cg'][s] * Eg / max(M, 1e-9):>12.5g}")
        print(f"      forward relay level  fwd_g / backward  bwd_g  "
              f"{st['fwd_g'][s]:>12.5g} {st['bwd_g'][s]:>12.5g}")
        print(f"      psi inputs: mo_g {uni['mo_g'][s]:>8.4f} at cm_g {uni['cm_g'][s]:>8.4g}   "
              f"lam {uni['lam_msg'][s]:>8.4f} at c_tau {uni['c_tau'][s]:>8.4g}")
        # ⭐⭐ THE RNA MEASUREMENT CHANNEL — the graft. A mature-RNA count at the flanking junction is
        # a CERTIFIED-RNA measurement of this exon's RNA density (gDNA cannot be spliced), and on a
        # slot with no gDNA evidence and no strand it is the ONLY thing that can move f_g off psi's
        # uninformative reference. It is what the gDNA channels being dead leaves behind.
        cmp_, cmn_ = uni["cm_p"][s], uni["cm_n"][s]
        rho_r = (np.exp(uni["mo_p"][s]) + np.exp(uni["mo_n"][s])) * M / max(Er, 1e-9)
        print(f"      ⭐ RNA measurement: mo_p {uni['mo_p'][s]:>8.4f} at cm_p {cmp_:>9.4g}   "
              f"mo_n {uni['mo_n'][s]:>8.4f} at cm_n {cmn_:>9.4g}")
        print(f"         ⇒ delivered RNA density {rho_r:>12.5g}   against the slot's own total "
              f"M/E_r = {M / max(Er, 1e-9):>10.5g}")
        print(f"         ⇒ implied f_g = 1 − rho_R·E_r/M = "
              f"{1.0 - rho_r * Er / max(M, 1e-9):>10.5g}   (final {cap['f_g'][s]:.4f})")
        print(f"      junction flux at the flanking lines: "
              f"{float(np.asarray(cap['mature'], float)[s - 1]):,.0f} / "
              f"{float(np.asarray(cap['mature'], float)[s + 1]):,.0f}"
              f"   E_J-frame density = J/E_J")
    return cap, st, uni, rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--conditions", nargs="+", required=True)
    ap.add_argument("--n-rna", type=int, default=20000)
    ap.add_argument("--genome-length", type=int, default=0)
    ap.add_argument("--work-dir", type=Path, default=Path("/tmp/rigel_dissect"))
    args = ap.parse_args()
    for c in args.conditions:
        dissect(c, n_rna=args.n_rna, genome_length=args.genome_length, work_dir=args.work_dir)
        print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
