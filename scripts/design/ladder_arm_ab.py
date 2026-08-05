#!/usr/bin/env python
"""⭐⭐ PROTOTYPE A SOLVER MECHANISM AS A PER-OBJECT OVERRIDE AND SCORE IT ON THE REAL LADDER — before
writing a line of it into ``src/``.

⛔ **THE WORKFLOW THIS EXISTS FOR.** A toy isolates a mechanism; it cannot rank one (`TESTING.md` §0b).
`toy_ceiling.py` re-solves the owner's two-exon toy under a set of arms and says what each is worth
THERE. This applies the SAME override to the 36-condition gDNA ladder and scores it with the project's
own acceptance instrument (`solvability_audit.py`), so the two arms differ by exactly the mechanism and
the answer is in the currency the project actually ships.

⭐ **The override goes in at `node_init.build_node_init`** — the per-object message-free self-solve —
which is the one place a mechanism can be expressed without touching either relay twin, so a prototype
cannot accidentally be gated in one twin and not the other (`TRAPS.md` B14). `region_arrays` reaches the
wrapper by way of a `node_sweep` wrapper, because `build_node_init` is not handed it.

⚠ **A prototype here is an UPPER bound on the built form, not a model of it**: it overwrites the
object's own belief unconditionally, with no licence and no inverse-variance fuse. If the upper bound
loses, the built form cannot win; if it wins, the built form still has to earn it.

⛔⛔ **THE ARM IT CARRIES, AND ITS VERDICT (2026-08-04).** ``intron_phi`` is `EQUATIONS.md` §3.6 face (I):
an `intron|exon` EDGE takes the flanking INTRON's COMPOSITION — never its level — and closes the level
with its own observed mass. On the toy it was worth −0.009 panel-wide and −0.027 on capture-ON ×
unstranded. On the ladder, 32 contaminated conditions, node axis:

======================  ==========  =============  ===================
metric                  base        intron_phi     better/worse/flat
======================  ==========  =============  ===================
solvable mwae           0.0413      **0.0426**     14 / 18 / 0
confidently-wrong       20,173      **22,336**      6 / 26 / 0
declared-precision      2.688       **2.821**       6 / 26 / 0
relay Δ                 +0.0021     **+0.0034**    14 / 18 / 0
======================  ==========  =============  ===================

⛔ **REFUTED.** And on the EDGE axis it is worse still (mwae 0.0216 → 0.1449) because it ADMITS those
EDGEs into the scored population (``solv%`` 43.1 → 48.2) carrying a wrong composition — `TRAPS.md` B13
inverted. ⭐ The one real signal in it: the split is by STRAND, not by capture. Every unstranded row is
worse and every stranded capture-OFF row is better, which is the DL mismatch test having nothing to
check an imported composition against when the destination has no self-solve.

⭐ C11 is the precedent that makes this run mandatory rather than optional: its toy delta was LARGER
(−0.035) and the ladder said it cost the panel nothing. `TRAPS.md` B1.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np  # noqa: E402

DESIGN = Path(__file__).resolve().parent
sys.path.insert(0, str(DESIGN))


def _sibling(name: str):
    key = name[:-3]
    if key not in sys.modules:
        spec = importlib.util.spec_from_file_location(key, DESIGN / name)
        mod = importlib.util.module_from_spec(spec)
        sys.modules[key] = mod
        spec.loader.exec_module(mod)
    return sys.modules[key]


SA = _sibling("solvability_audit.py")
P0 = _sibling("pass0_vs_oracle.py")

import rigel.calibration.strand_balance  # noqa: E402,F401  (registers the module for patching)
from rigel.calibration import bp_solver, node_init as NI  # noqa: E402
from rigel.calibration.node_chain import NODE  # noqa: E402
from rigel.calibration.node_geometry import node_global_geometry  # noqa: E402
from rigel.calibration.signature import coarse_type_array  # noqa: E402
from rigel.config import CalibrationConfig, PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

CAL = sys.modules["rigel.calibration.calibrate"]

_EPS = 1.0e-9
INTRON, EXON = 1, 2
#: set by the `node_sweep` wrapper one call before `build_node_init` needs it.
_RA: dict = {}


def _install_kappa_half():
    """⭐⭐ THE ARM `ROADMAP.md` §1 step 3 ASKS FOR IN ITS OWN WORDS: "inject κ = 0.5 exactly and diff".

    On a genuinely unstranded library the strand channel carries **exactly** zero information about
    composition — the Fisher information is ``∝ (2κ−1)²`` and κ is ½. But κ is FITTED, so it lands on
    0.500689, and the tool reads that 0.000689 as evidence: ``disc = 4·max(0, (κ̂−½)² − σ²_d)`` survives
    the derived noise floor about one unstranded library in three, and a **boolean** composition licence
    is then flipped at full strength by it.

    Forcing κ to exactly ½ makes ``disc`` identically 0, which is what propagating ``Var(κ̂)`` would
    achieve by derivation — so on the UNSTRANDED conditions this is the ceiling on that fix.
    ⛔ On the STRANDED conditions it is not a fix at all, it is a DESTRUCTION TEST: κ̂ is ~97 standard
    errors from ½ there and the channel is carrying real signal, so those rows must get much worse. If
    they do not, the lever is not connected and no reading of the unstranded rows means anything.

    ⭐ Only κ moves. ``n_observations`` is carried through unchanged, so the noise floor ``σ²_d`` — which
    depends on the SAMPLE SIZE, not on κ — is identical in both arms; overriding it too would conflate
    "κ is exactly ½" with "there are no spliced reads"."""
    import dataclasses as _dc

    SB = sys.modules["rigel.calibration.strand_balance"]
    orig = CAL.fit_strand_balance

    def wrapper(strand_model):
        return _dc.replace(orig(strand_model), rna_sense_frac=0.5)

    CAL.fit_strand_balance = wrapper
    SB.fit_strand_balance = wrapper


def _targets(chain, region_arrays):
    """Vectorised: every `intron|exon` EDGE slot and the slot of its flanking INTRON NODE.

    ⭐ Genomic order IS slot order and the chain strictly alternates NODE/EDGE, so an EDGE's flanks are
    ``i-1`` and ``i+1`` and no traversal is needed. Breaks at a reference terminal are handled by
    `chain.left`/`chain.right` being −1 there."""
    kind = np.asarray(chain.kind)
    obj = np.asarray(chain.obj_idx, np.int64)
    rtype = coarse_type_array(np.asarray(region_arrays.signature)).astype(np.int64)
    is_node = kind == NODE
    ntype = np.where(is_node, rtype[np.clip(obj, 0, rtype.shape[0] - 1)], -1)
    left = np.asarray(chain.left, np.int64)
    right = np.asarray(chain.right, np.int64)
    ok = (~is_node) & (left >= 0) & (right >= 0)
    lt = np.where(ok, ntype[np.clip(left, 0, ntype.size - 1)], -1)
    rt = np.where(ok, ntype[np.clip(right, 0, ntype.size - 1)], -1)
    hit = ok & (((lt == INTRON) & (rt == EXON)) | ((lt == EXON) & (rt == INTRON)))
    edges = np.flatnonzero(hit)
    srcs = np.where(lt[edges] == INTRON, left[edges], right[edges])
    return edges, srcs


def _wrap_node_sweep():
    """Stash `region_arrays` — `node_sweep` receives it and calls `build_node_init` after."""
    orig = CAL.node_sweep

    def wrapper(chain, statics, geometry, belief, region_arrays, *a, **kw):
        _RA["region_arrays"] = region_arrays
        return orig(chain, statics, geometry, belief, region_arrays, *a, **kw)

    CAL.node_sweep = wrapper
    return orig


def _install_face_one():
    """⭐ FACE (I) AS DERIVED (`EQUATIONS.md` §3.6): the `intron|exon` EDGE takes the flanking INTRON's
    COMPOSITION — never its level — and closes the level with its OWN observed mass.

    ``T(INTRON) = {gDNA, nascent} = T(EDGE, unspliced only)``, so the two measure the same population and
    the share may cross. ⛔ The LEVEL may not: measured on the toy, transferring the intron's gDNA
    DENSITY instead is **+0.017 panel-wide and +0.207 on capture-ON × unstranded**, because capture
    depletes the intron ~1000× while enriching the EDGE."""
    orig = NI.build_node_init

    def wrapper(chain, statics, geometry, **kw):
        ni = orig(chain, statics, geometry, **kw)
        ra = _RA.get("region_arrays")
        if ra is None:
            return ni
        edges, srcs = _targets(chain, ra)
        if edges.size == 0:
            return ni
        f_g = np.array(ni.f_g, np.float64)
        f_pos = np.array(ni.f_pos, np.float64)
        f_neg = np.array(ni.f_neg, np.float64)
        tau = np.array(ni.tau_lam, np.float64)
        lock = np.asarray(ni.struct_lock, bool)
        M, E_g = node_global_geometry(geometry)
        M = np.asarray(M, np.float64)
        E_g = np.asarray(E_g, np.float64)
        E_r = np.asarray(geometry.eff_rna, np.float64)
        n_node = np.asarray(geometry.unspliced_count, np.float64).sum(axis=1)

        new_fg = f_g[srcs]
        rna = np.maximum(0.0, 1.0 - new_fg)
        tot = f_pos[edges] + f_neg[edges]
        fp_ok = np.asarray(statics.free_pos, bool)[edges]
        fn_ok = np.asarray(statics.free_neg, bool)[edges]
        k = fp_ok.astype(np.float64) + fn_ok.astype(np.float64)
        share_p = np.where(tot > _EPS, f_pos[edges] / np.maximum(tot, _EPS),
                           np.where(k > 0, fp_ok / np.maximum(k, 1.0), 0.0))
        share_n = np.where(tot > _EPS, f_neg[edges] / np.maximum(tot, _EPS),
                           np.where(k > 0, fn_ok / np.maximum(k, 1.0), 0.0))
        f_pos[edges] = rna * share_p
        f_neg[edges] = rna * share_n
        f_g[edges] = new_fg
        tau[edges] = tau[srcs]

        v_fg, v_fr = NI.own_composition_logvar(f_g, tau, lock)
        rho_g = np.maximum(
            np.where((M > _EPS) & (E_g > _EPS), f_g * M / np.maximum(E_g, _EPS), 0.0), 0.0
        )
        prec_g = NI.own_precision(n_node, v_fg, rho_g > _EPS)
        touched = np.zeros(f_g.shape[0], bool)
        touched[edges] = True

        def _rna(frac, free_s, rho_old):
            raw = np.where(
                (M > _EPS) & (E_r > _EPS) & np.asarray(free_s, bool),
                frac * M / np.maximum(E_r, _EPS), 0.0,
            )
            live = (n_node > 0.0) & (raw > _EPS) & ((rho_old > 0.0) | touched)
            rho = np.where(live, raw, 0.0)
            return rho, NI.own_precision(n_node, v_fr, rho > _EPS)

        rho_pos, prec_pos = _rna(f_pos, statics.free_pos, np.asarray(ni.rho_pos, np.float64))
        rho_neg, prec_neg = _rna(f_neg, statics.free_neg, np.asarray(ni.rho_neg, np.float64))
        return NI.NodeInit(
            f_g=f_g, f_pos=f_pos, f_neg=f_neg,
            rho_g=rho_g, rho_pos=rho_pos, rho_neg=rho_neg,
            prec_g=prec_g, prec_pos=prec_pos, prec_neg=prec_neg,
            struct_lock=lock, tau_lam=tau,
        )

    bp_solver.build_node_init = wrapper


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--arm", choices=("base", "intron_phi", "kappa_half"), required=True)
    ap.add_argument("--conditions", nargs="*", default=None)
    ap.add_argument("--suite", type=Path, default=P0.DEFAULT_SUITE)
    ap.add_argument("--index", type=Path, default=P0.DEFAULT_INDEX)
    ap.add_argument("--oracle-cache", type=Path, default=None)
    ap.add_argument("--work-dir", type=Path, default=Path("/tmp/rigel_ladder_ceiling"))
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    _wrap_node_sweep()
    if args.arm == "intron_phi":
        _install_face_one()
    elif args.arm == "kappa_half":
        _install_kappa_half()

    index = TranscriptIndex.load(str(args.index))
    config = CalibrationConfig()
    names = args.conditions or sorted(
        p.name for p in args.suite.iterdir() if (p / "sim_oracle.bam").is_file()
    )
    with args.out.open("w") as fh:
        for name in names:
            t0 = time.perf_counter()
            cond = args.suite / name
            truth = P0.truth_f_gdna(cond) or 0.0
            m = P0.measure_condition(
                bam=str(cond / "sim_oracle.bam"), index=index, pipeline_config=PipelineConfig(),
                calibration_config=config, work_dir=args.work_dir / "rigel_pass0_oracle", tag=name,
                truth_pmfs=lambda size, d=cond: (
                    P0.truth_length_pmf(d, "gdna", size), P0.truth_length_pmf(d, "rna", size)
                ),
                oracle_cache=args.oracle_cache,
            )
            for axis in ("node", "edge"):
                s = SA.summarise(SA.audit(m, axis=axis, config=config))
                # ⛔⛔ FIXED-DENOMINATOR COMPANIONS. `solvable_mwae`'s denominator is the SOLVABLE set,
                # and an arm that changes what counts as solvable changes its own denominator — so a
                # fall in it can be "solved better" or "stopped scoring the hard ones" and the number
                # cannot tell you which (`TRAPS.md` B13/B20). These three cannot move that way:
                #   * `mwae_all`  — every object with mass, so the denominator is the library.
                #   * `abs_err_all` — the raw error MASS, no denominator at all.
                #   * `library_f_gdna` — the deliverable itself, one number per condition.
                sc = m.scores["pass0"][axis]["ALL"]
                s["mwae_all"] = float(sc.mwae)
                s["abs_err_all"] = float(sc.abs_err)
                s["mass_all"] = float(sc.mass)
                s["net_err_all"] = float(sc.net_err)
                fin = m.scores["final"][axis]["ALL"]
                s["mwae_all_final"] = float(fin.mwae)
                s["abs_err_all_final"] = float(fin.abs_err)
                s["library_f_gdna_pass0"] = float(m.library_f_gdna.get("pass0", float("nan")))
                s["library_f_gdna_final"] = float(m.library_f_gdna.get("final", float("nan")))
                s["library_f_gdna_truth"] = float(m.library_f_gdna.get("T", float("nan")))
                fh.write(json.dumps({"arm": args.arm, "condition": name, "axis": axis,
                                     "f_gdna": truth, **s}) + "\n")
                fh.flush()
            print(f"  {name} {time.perf_counter() - t0:.0f}s", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
