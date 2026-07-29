"""Is pass-0 subtracting the spliced (mature) RNA to leave gDNA, at SINGLE-STRAND EXON nodes?
For single-strand exon region nodes, compare pass-0 f_g to oracle f_g and stratify by the oracle MATURE-RNA
content (what the spliced signal should be subtracting). If the subtraction works, nodes with lots of mature RNA
should have their RNA removed → f_g near oracle. Systematic OVER-call (pass0 f_g > oracle) on mature-rich nodes =
the subtraction is NOT happening.
"""
from __future__ import annotations
import dataclasses
import os
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np
import sys

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from selfsolve_diag import _scan_and_truth  # noqa: E402
from flagship_interrogate import _oracle_per_node  # noqa: E402
from rigel.calibration import calibrate  # noqa: E402
from rigel.calibration.bp_solver import REGION  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.signature import BIT_EXON_NEG, BIT_EXON_POS  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

_EXON = BIT_EXON_POS | BIT_EXON_NEG
suite = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
index = TranscriptIndex.load(str(suite / "rigel_index"))
cfg = PipelineConfig()
ra = RegionArrays.from_index(index)


def pass0(cond):
    inp = _scan_and_truth(suite, cond, index, cfg, Path("/tmp/rigel_selfsolve"), suite / "_selfsolve_cache")
    cc = dataclasses.replace(cfg.calibration, calib_refit_iters=0)  # PASS-0 only
    dbg = {}
    calibrate(inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
              np.asarray(inp["rna_fl_pmf"]), cc, _debug=dbg)
    chain, cap = dbg["chain"], dbg["capture"]
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
    rp = inp["region_pools"]
    mat = np.asarray(rp["mat_uns_pos"], float) + np.asarray(rp["mat_uns_neg"], float)
    nas = np.asarray(rp["nas_uns_pos"], float) + np.asarray(rp["nas_uns_neg"], float)
    spl = np.asarray(rp["mat_spl_pos"], float) + np.asarray(rp["mat_spl_neg"], float) if "mat_spl_pos" in rp else None
    ridx = np.asarray(chain.ref_idx, np.int64)
    isr = np.asarray(chain.kind) == REGION
    mass = np.asarray(cap["mass_global"], float); eff = np.asarray(cap["eff_global"], float)
    fp = np.asarray(cap["free_pos"], bool); fn = np.asarray(cap["free_neg"], bool)
    live = (eff > 1e-9) & (mass > 1e-12)
    sig = np.asarray(ra.signature)
    sg = np.where((ridx >= 0) & (ridx < len(sig)), sig[np.clip(ridx, 0, len(sig) - 1)], 0)
    G, R = Gp + Gn, Rp + Rn
    fo = np.where((G + R) > 0, G / np.maximum(G + R, 1e-12), np.nan)
    fg0 = np.asarray(dbg["belief"].f_g)
    matn = np.where(isr, mat[np.clip(ridx, 0, len(mat) - 1)], 0.0)
    nasn = np.where(isr, nas[np.clip(ridx, 0, len(nas) - 1)], 0.0)
    # single-strand EXON region nodes with data
    ss_exon = live & isr & (fp ^ fn) & (sg & _EXON > 0)
    return ss_exon, fo, fg0, matn, nasn, G, R, mass


print(f"{'condition':40s} {'n':>4s} {'oracle':>6s} {'pass0':>6s} {'|err|':>6s} {'over%':>6s} {'corr':>5s}")
for cond in ["gdna_gdna100_ss_0.99_nrna_present_capture_on",
             "gdna_gdna100_ss_0.50_nrna_present_capture_on",
             "gdna_gdna300_ss_0.50_nrna_present_capture_on",
             "gdna_gdna5_ss_0.50_nrna_present_capture_on"]:
    sel, fo, fg0, matn, nasn, G, R, mass = pass0(cond)
    s = sel & np.isfinite(fo)
    oc, p0 = fo[s], fg0[s]
    err = np.abs(p0 - oc)
    over = np.mean(p0 > oc + 0.1)
    corr = np.corrcoef(oc, p0)[0, 1] if s.sum() > 2 else np.nan
    print(f"{cond.replace('gdna_',''):40s} {int(s.sum()):>4d} {oc.mean():6.3f} {p0.mean():6.3f} {err.mean():6.3f} {over*100:5.0f}% {corr:5.2f}")
    # stratify by mature-RNA fraction of the node's unspliced RNA
    matfrac = matn[s] / np.maximum(matn[s] + nasn[s], 1e-9)
    hi = matfrac > 0.6  # mature-dominated RNA (spliced should subtract)
    if hi.sum() > 2:
        print(f"     mature-RNA-rich nodes (n={int(hi.sum())}): oracle {oc[hi].mean():.3f} pass0 {p0[hi].mean():.3f}  (spliced should pull pass0 DOWN toward oracle)")
