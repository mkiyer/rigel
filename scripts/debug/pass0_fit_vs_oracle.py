"""Figure: the pass-0 gDNA-density prior fit vs the ORACLE, for the crushing scenarios — 'how badly do we miss
the mode based on pass-0 data?'. For each scenario: the ORACLE gDNA-density landscape P(log10 rho_g) (what the
prior SHOULD be) vs the density FIT ON PASS-0's deconvolved gDNA (what we actually get), on the RELEASED v0.7.1
training substrate (single-strand + structural region nodes, incl zero-count). node 1055's true / pass-0 / max
rho_g are marked so the miss is visible: the oracle has an ENRICHED mode; the pass-0 fit does NOT, because
pass-0 under-called every unstranded enriched node (the circularity). Writes docs/figures/pass0_npmle_vs_oracle.png."""
from __future__ import annotations

import dataclasses
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from selfsolve_diag import _scan_and_truth  # noqa: E402
from flagship_interrogate import _oracle_per_node  # noqa: E402
from rigel.calibration import calibrate  # noqa: E402
from rigel.calibration.bp_solver import REGION  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.signature import coarse_type_array  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

_EPS = 1e-12
GRID = np.linspace(-5.0, 2.5, 320)  # log10 rho_g
H = 0.15  # decades
suite = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
index = TranscriptIndex.load(str(suite / "rigel_index"))
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
cfg = PipelineConfig()
work = Path("/tmp/rigel_selfsolve")
cache = suite / "_selfsolve_cache"

PANELS = [
    ("gdna_gdna100_ss_0.50_nrna_none_capture_on", "gdna100 ss0.50 capON  (the crush)"),
    ("gdna_gdna300_ss_0.50_nrna_none_capture_on", "gdna300 ss0.50 capON"),
    ("gdna_none_ss_0.50_nrna_none_capture_on", "none ss0.50 capON  (zero gDNA)"),
    ("gdna_gdna100_ss_0.99_nrna_none_capture_on", "gdna100 ss0.99 capON  (stranded control)"),
]


def render(a_cells):
    """UNIT mass per node, fixed bandwidth H (decades)."""
    z = (GRID[:, None] - a_cells[None, :]) / H
    d = (np.exp(-0.5 * z * z) / (H * np.sqrt(2 * np.pi))).sum(axis=1)
    return d / max(float(d.sum()), _EPS)


def pass0(inp):
    cc = dataclasses.replace(cfg.calibration, calib_refit_iters=0)
    dbg = {}
    calibrate(inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
              np.asarray(inp["rna_fl_pmf"]), cc, _debug=dbg)
    return dbg


fig, axes = plt.subplots(2, 2, figsize=(13, 8))
axes = axes.ravel()
for ax, (cond, title) in zip(axes, PANELS):
    inp = _scan_and_truth(suite, cond, index, cfg, work, cache)
    dbg = pass0(inp)
    chain = dbg["chain"]
    cap = dbg["capture"]
    b0 = dbg["belief"]
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
    G = Gp + Gn
    Rt = Rp + Rn
    mass = np.asarray(cap["mass_global"], float)
    eff = np.asarray(cap["eff_global"], float)
    fp = np.asarray(cap["free_pos"], bool)
    fn = np.asarray(cap["free_neg"], bool)
    isr = np.asarray(chain.kind) == REGION
    idx = np.asarray(chain.ref_idx, np.int64)
    rtype = coarse_type_array(np.asarray(ra.signature)).astype(np.int64)
    ntype = np.where(isr, rtype[np.clip(idx, 0, len(rtype) - 1)], -1)
    single = fp ^ fn
    gonly = ~fp & ~fn
    live = (eff > 1e-9) & (mass > _EPS)
    sel = live & isr & (single | gonly)  # RELEASED v0.7.1 substrate: single + structural
    f0 = np.asarray(b0.f_g, float)
    ghat0 = f0[sel] * mass[sel]
    a_pass0 = np.log10(np.maximum(ghat0, 1.0)) - np.log10(eff[sel])
    a_orc = np.log10(np.maximum(G[sel], 1.0)) - np.log10(eff[sel])
    dO = render(a_orc)
    dP = render(a_pass0)
    ax.fill_between(GRID, dO, color="0.80", label="ORACLE landscape (target)")
    ax.plot(GRID, dO, "k", lw=1.6)
    ax.plot(GRID, dP, "tab:red", lw=1.8, label="fit on PASS-0 gDNA (what we get)")

    # node 1055 markers (captured single exon, highest true f_g)
    exsel = live & isr & (ntype == 2) & single & (G / np.maximum(G + Rt, _EPS) > 0.6)
    if exsel.any():
        i = int(np.argmax(G * exsel))
        rho_tot = mass[i] / eff[i]
        fo = G[i] / max(G[i] + Rt[i], _EPS)
        true_lrg = np.log10(max(fo * rho_tot, 1e-9))
        p0_lrg = np.log10(max(f0[i] * rho_tot, 1e-9))
        acc = np.log10(rho_tot)
        ax.axvline(true_lrg, color="tab:green", lw=2.0, ls="-",
                   label=f"node {i} TRUE rho_g ({true_lrg:+.2f}, f_g={fo:.2f})")
        ax.axvline(p0_lrg, color="tab:orange", lw=1.5, ls="--",
                   label=f"node {i} PASS-0 rho_g ({p0_lrg:+.2f}, f_g={f0[i]:.2f})")
        ax.axvline(acc, color="0.5", lw=1.0, ls=":", label=f"accessible max log(M/E) ({acc:+.2f})")
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("log10 gDNA density  rho_g", fontsize=8)
    ax.set_xlim(-5, 2.5)
    ax.legend(fontsize=6.5, loc="upper right")
    ax.tick_params(labelsize=7)
fig.suptitle("Pass-0 gDNA prior vs the ORACLE landscape (released v0.7.1 substrate). The oracle has an "
             "ENRICHED mode; the pass-0 fit misses it — the circularity.", fontsize=12)
fig.tight_layout(rect=(0, 0, 1, 0.98))
out = Path("/Users/mkiyer/proj/rigel/docs/figures/pass0_npmle_vs_oracle.png")
fig.savefig(out, dpi=120)
print(f"wrote {out}")
