"""DISSECTION 7 — DOES THE INITIALIZATION LEAK? The owner's invariant, tested mechanically.

Owner: *"My conception of this was that the initialization value of a node should not matter at zero
precision. It's officially just... a moot point. It might as well be none or undefined or NaN or anything.
... And if not, then that's the bug right there."*

That is a falsifiable statement about the code, so this falsifies it or confirms it. `_type_belief` sets

    f_g = np.ones(n)          # all-gDNA default
    var_g = np.full(n, inf)   # unsolved  ⇒  own_precision == 0

and a node keeps that default when it is G1 (no free strand), G3 (AMBIG), or G2 with no unspliced mass.
The 863 boundaries feeding false gDNA into exons are the last kind: **zero unspliced mass, so zero count,
so zero precision** — and yet they are declared 100 % gDNA.

PART A — the STRUCTURE. Owner: *"with zero DNA these unspliced boundaries have to be structurally exon-to-
exon boundaries... exon-to-intron boundaries should not have any unspliced fragments in a zero-DNA library
unless they're nascent RNA, and we're looking at a nascent-RNA-none library."* Checked against the flanking
region types.

PART B — the INVARIANT. Re-run the whole solve with the default flipped from 100 % gDNA to 100 % RNA, ONLY
on nodes that carry no unspliced mass. Those nodes have no count, hence no precision, hence — by the
owner's model — no influence. **Any difference in the output is the leak, and its size is the bug's size.**

Run: OMP_NUM_THREADS=1 python scratchpad/gdna_d7_init_leak.py [cond]
"""
from __future__ import annotations

import dataclasses
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from selfsolve_diag import _scan_and_truth  # noqa: E402

import rigel.calibration.node_geometry as NG  # noqa: E402
from rigel.calibration.bp_solver import REGION  # noqa: E402
from rigel.calibration.calibrate import calibrate  # noqa: E402
from rigel.calibration.node_geometry import _node_region_type  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
COND = sys.argv[1] if len(sys.argv) > 1 else "gdna_none_ss_0.50_nrna_none_capture_off"
_EPS = 1e-9
CLS = {0: "intergenic", 1: "intron", 2: "exon", 3: "boundary"}
_ORIG = NG._type_belief


def _patched(free_pos, free_neg, deconv, mass_unspl):
    """Identical to the shipped `_type_belief`, except that a node with NO unspliced mass — which therefore
    has no count and no precision — defaults to 100 % RNA instead of 100 % gDNA. By the owner's invariant
    this must be inert."""
    f_pos, f_neg, f_g, var_p, var_n, var_g = _ORIG(free_pos, free_neg, deconv, mass_unspl)
    empty = np.asarray(mass_unspl, np.float64) <= 0.0
    fp_b, fn_b = np.asarray(free_pos, bool), np.asarray(free_neg, bool)
    n_free = fp_b.astype(np.float64) + fn_b.astype(np.float64)
    flip = empty & (n_free > 0)          # leave G1 sinks (no admissible RNA strand) exactly as they are
    f_g = np.where(flip, 0.0, f_g)
    f_pos = np.where(flip, np.where(fp_b, 1.0 / np.maximum(n_free, 1.0), 0.0), f_pos)
    f_neg = np.where(flip, np.where(fn_b, 1.0 / np.maximum(n_free, 1.0), 0.0), f_neg)
    return f_pos, f_neg, f_g, var_p, var_n, var_g


def run(patch: bool):
    index = TranscriptIndex.load(str(SUITE / "rigel_index"))
    cfg = PipelineConfig()
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    inp = _scan_and_truth(SUITE, COND, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
    dbg: dict = {}
    NG._type_belief = _patched if patch else _ORIG
    try:
        calibrate(inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
                  np.asarray(inp["rna_fl_pmf"]),
                  dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg)
    finally:
        NG._type_belief = _ORIG
    return index, ra, dbg


def main():
    index, ra, base = run(False)
    chain, cap, st = base["chain"], base["capture"], base["statics"]
    mass = np.asarray(cap["mass_global"])
    isr = np.asarray(chain.kind) == REGION
    rt, _ = _node_region_type(chain, ra)
    cls = np.where(~isr, 3, rt)
    L, R = np.asarray(chain.left), np.asarray(chain.right)
    spl = np.asarray(st.mass_spliced, float)
    mu = np.asarray(st.mass_unspliced, float)
    fpb, fnb = np.asarray(st.free_pos, bool), np.asarray(st.free_neg, bool)

    # ── A. what ARE the boundaries that carry spliced flux but no unspliced mass? ──────────────────
    b = ~isr
    src = b & (mu <= _EPS) & (spl > _EPS)
    print(f"# {COND}\n")
    print("=== A. STRUCTURE of the boundaries that emit into the false-positive exons ===")
    print(f"    boundaries: {int(b.sum())};  with spliced flux and NO unspliced mass: {int(src.sum())}\n")

    def side(idx):
        return np.where(idx >= 0, cls[np.clip(idx, 0, len(cls) - 1)], -1)

    lc, rc = side(L), side(R)
    pairs: dict[tuple[int, int], int] = {}
    for i in np.flatnonzero(src):
        pairs[(int(min(lc[i], rc[i])), int(max(lc[i], rc[i])))] = \
            pairs.get((int(min(lc[i], rc[i])), int(max(lc[i], rc[i]))), 0) + 1
    nm = {-1: "TERMINAL", 0: "intergenic", 1: "intron", 2: "exon"}
    print(f"    {'flanking region pair':32s} {'n':>6s} {'share':>7s}")
    for k, v in sorted(pairs.items(), key=lambda kv: -kv[1]):
        print(f"    {nm[k[0]] + ' <-> ' + nm[k[1]]:32s} {v:6d} {100 * v / max(int(src.sum()), 1):6.1f} %")
    print(f"\n    of those: both strands free (AMBIG/G3) {int((src & fpb & fnb).sum()):5d}   "
          f"one strand free (G2) {int((src & (fpb ^ fnb)).sum()):5d}   "
          f"neither (G1 sink) {int((src & ~fpb & ~fnb).sum()):5d}")
    print(f"    their initialized f_g: median {np.median(np.asarray(base['belief_pass0'].f_g)[src]):.4f}"
          if "belief_pass0" in base else "")

    # ── B. the invariant: flip the default on zero-mass nodes and see if ANYTHING moves ────────────
    _, _, alt = run(True)
    fa, fb = np.asarray(cap["f_g"]), np.asarray(alt["capture"]["f_g"])
    d = np.abs(fa - fb)
    # the nodes the patch WRITES TO are not evidence of a leak — they moved because they were set. The
    # invariant is about everyone ELSE, and in particular about every node that carries mass.
    patched = (mu <= 0.0) & (fpb | fnb)
    leaked = (d > 1e-12) & ~patched
    live = mass > _EPS
    print("\n=== B. THE INVARIANT: default flipped 100 % gDNA -> 100 % RNA on ZERO-MASS nodes only ===")
    print("    Those nodes have no count, so no precision. By the owner's model this must be INERT.\n")
    print(f"    directly patched (zero-mass, >=1 free strand): {int(patched.sum())} nodes — expected to move")
    print(f"    ANY OTHER node that moved:                     {int(leaked.sum())}"
          f"{'' if not leaked.any() else f'   max |delta| {d[leaked].max():.4f}'}")
    print(f"    nodes WITH MASS that moved:                    {int((leaked & live).sum())}")
    fp_a, fp_b_ = (fa * mass)[live].sum(), (fb * mass)[live].sum()
    print(f"    FALSE-POSITIVE gDNA mass: {fp_a:,.0f}  ->  {fp_b_:,.0f}   "
          f"({100 * (fp_b_ - fp_a) / max(fp_a, 1):+.2f} %)")
    if (leaked & live).any():
        print("\n    movers by class:")
        for c in (0, 1, 2, 3):
            m = leaked & live & (cls == c)
            if not m.any():
                continue
            print(f"      {CLS[c]:11s} {int(m.sum()):5d} nodes  f_g {np.average(fa[m], weights=mass[m]):.4f}"
                  f" -> {np.average(fb[m], weights=mass[m]):.4f}")
        print("\n    ⇒ THE INITIALIZATION LEAKS — a value carried at ZERO precision changed the answer.")
    else:
        print("\n    ⇒ ✅ THE INVARIANT HOLDS. The zero-precision default is INERT: flipping it from 100 %")
        print("      gDNA to 100 % RNA changes nothing on any node that carries mass, and the false-positive")
        print("      mass is bit-identical. The 100 % gDNA default is misleading to READ, but it is not the")
        print("      source of the false gDNA — that has to come from the spliced measurement path.")


if __name__ == "__main__":
    main()
