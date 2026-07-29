"""Is the gDNA FACTORY severed at the intergenic-exon SEAM? (owner, 2026-07-23)

Owner's model: intergenic regions bookend every locus. They are the ANCHORS — they measure pure gDNA with
COUNT precision, they SOURCE gDNA into the locus, and they SINK messages (no RNA crosses them, so a message
that reaches one dies). Splice junctions are the mirror: they source and sink RNA, not gDNA. Those are the
only two places measured signal enters the chain.

The intergenic-exon SEAM boundary is exclusively gDNA too — no RNA strand is continuous across a TSS/TES — so
by the same argument it is composition-CERTAIN and should relay the factory's gDNA at counting precision.

But `bp_solver._compile_strand_evidence` sets

    struct_lock = locked & is_region                       <-- BOUNDARIES EXCLUDED

so a seam falls through to the tau path. A strand-free node has I_strand = 0 (the deadband, and
`(f_g(1-f_g))^2 = 0` at f_g=1), so tau = 0 => `v_logfg = inf` => `_pred_precision` returns 0 => the seam
emits NOTHING. The factory reaches the seam and stops one hop short of the locus.

This measures, per edge family: the message gDNA precision, the source's lock status, and the source count.
Reference scale: a message above 1/Var_Beta(1/2,1/2)(log f) = 0.304 is stronger than the destination's prior.

Usage:  OMP_NUM_THREADS=1 python scripts/debug/intergenic_factory_probe.py
"""
from __future__ import annotations

import dataclasses
import importlib
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from selfsolve_diag import _scan_and_truth  # noqa: E402

from rigel.calibration.bp_solver import REGION  # noqa: E402
from rigel.calibration.node_geometry import _node_region_type  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")

_EPS = 1e-9
REF_PREC = 0.3040
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
INTERGENIC, INTRON, EXON = 0, 1, 2

index = TranscriptIndex.load(str(SUITE / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_index(index)
conds = sorted(d.name for d in SUITE.iterdir() if (d / "sim_oracle.bam").exists())

fam: dict[tuple, list] = defaultdict(list)

for cond in conds:
    capst = "ON " if "capture_on" in cond else ("VSTR" if "verystrong" in cond else "OFF")
    inp = _scan_and_truth(SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
    dbg: dict = {}
    cc = dataclasses.replace(cfg.calibration, calib_refit_iters=0)
    calmod.calibrate(
        inp["payload"], ra, inp["strand_model"],
        np.asarray(inp["gdna_fl_pmf"]), np.asarray(inp["rna_fl_pmf"]), cc, _debug=dbg,
    )
    chain, cap, st = dbg["chain"], dbg["capture"], dbg["statics"]
    kind = np.asarray(chain.kind)
    rtype, _ = _node_region_type(chain, ra)
    fp = np.asarray(st.free_pos, bool)
    fn = np.asarray(st.free_neg, bool)
    no_rna = ~(fp | fn)  # structurally gDNA-only: intergenic region OR a no-RNA-crossing seam

    for e in cap.get("_edge_modes", []):
        s, d = e["src"], e["dst"]
        if kind[s] == REGION and rtype[s] == INTERGENIC and kind[d] != REGION:
            name = "intergenic-REGION -> seam"
        elif kind[s] != REGION and no_rna[s] and kind[d] == REGION and rtype[d] == EXON:
            name = "SEAM -> exon  (the factory hop)"
        elif kind[s] != REGION and no_rna[s] and kind[d] == REGION and rtype[d] == INTRON:
            name = "SEAM -> intron"
        elif kind[s] == REGION and rtype[s] == INTRON and kind[d] != REGION:
            name = "intron-REGION -> junction"
        else:
            continue
        fam[(capst, name)].append((e["prec_g"], e["sm"], float(no_rna[s])))

print("\nmessage gDNA precision by edge family. Reference scale 1/Var_Beta(1/2,1/2)(log f) = 0.304")
print("-- a source that is composition-CERTAIN should emit at COUNT precision n/(n*s2t+1), NOT zero.\n")
hdr = (
    f"{'capture':>7} {'edge family':>30} | {'edges':>7} | {'prec_g=0':>9} {'med prec':>9} {'>ref':>6} |"
    f" {'med src n':>10}"
)
print(hdr)
print("-" * len(hdr))
for capst in ("OFF", "ON ", "VSTR"):
    for k in sorted(fam):
        if k[0] != capst:
            continue
        v = np.array(fam[k], dtype=float)
        p, sm = v[:, 0], v[:, 1]
        live = sm > _EPS  # the source actually had mass to speak with
        if live.sum() < 3:
            continue
        pl, sml = p[live], sm[live]
        print(
            f"{capst:>7} {k[1]:>30} | {live.sum():>7} | {np.mean(pl <= _EPS):>8.1%}"
            f" {np.median(pl):>9.4f} {np.mean(pl > REF_PREC):>6.1%} | {np.median(sml):>10.2f}"
        )
    print()
