"""Build the ENRICHED BOUNDARY-MODE STUDY PANEL (owner's setup, 2026-07-22).

A curated panel of splice-junction exon-intron boundary TRIPLETS (exon region → exon-intron boundary → intron
region) selected as HIGH-ERROR in the worst ambig_dense_10mb scenarios, then tracked across ALL scenarios as
controls. For each panel triplet, across every scenario, we store the full message observables on BOTH incoming
edges (exon→boundary, intron→boundary) so a candidate transfer MODE can be studied/derived OFFLINE (no re-run):

  triplet oracles/solves:  bnd_fo/bnd_fg, exon_fo/exon_fg, intron_fo/intron_fg, bnd_mass, bnd_spl
  exon→boundary edge:      ex_rho_g, ex_rho_r, ex_rna_src, ex_mat_abs, ex_egd, ex_erd, ex_md, ex_sm, ex_mode_g, ex_prec_g
  intron→boundary edge:    in_rho_g, in_rho_r, in_egd, in_erd, in_md, in_sm, in_mode_g, in_prec_g

The dominant error is the MATURE-CONTAMINATION under-call: at high-gDNA / no-nascent junctions the crossing is
pure gDNA (bnd_fo≈1) but the exon's within-exon mature is imputed as the boundary's nascent, crushing bnd_fg.
Output: `<scratch>/boundary_panel.npz`."""
from __future__ import annotations
import dataclasses
import importlib
import sys
from pathlib import Path
import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from selfsolve_diag import _scan_and_truth
from flagship_interrogate import _oracle_per_node
calmod = importlib.import_module("rigel.calibration.calibrate")
from rigel.calibration.bp_solver import REGION
from rigel.calibration.node_geometry import _node_region_type
from rigel.calibration.region_arrays import RegionArrays
from rigel.config import PipelineConfig
from rigel.index import TranscriptIndex

_EPS = 1e-9
K_PER_SCENARIO = 12
N_WORST = 3
SCRATCH = Path("/private/tmp/claude-503/-Users-mkiyer-proj-rigel/4f7a248b-0c78-4b40-9030-462373aefb19/scratchpad")
suite = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
index = TranscriptIndex.load(str(suite / "rigel_index")); cfg = PipelineConfig()
ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
conds = sorted(d.name for d in suite.iterdir()
               if (d / "sim_oracle.bam").exists() and d.name.startswith("gdna_"))

# ---- structural: chain, node types, exon-intron boundary set (scenario-independent) ----
inp0 = _scan_and_truth(suite, conds[0], index, cfg, Path("/tmp/rigel_selfsolve"), suite / "_selfsolve_cache")
dbg0 = {}; calmod.calibrate(inp0["payload"], ra, inp0["strand_model"], np.asarray(inp0["gdna_fl_pmf"]),
                            np.asarray(inp0["rna_fl_pmf"]), dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg0)
chain = dbg0["chain"]
node_type, _ = _node_region_type(chain, ra)
kind = np.asarray(chain.kind); left = np.asarray(chain.left); right = np.asarray(chain.right)
lt = np.where(left >= 0, node_type[np.clip(left, 0, None)], -1)
rt = np.where(right >= 0, node_type[np.clip(right, 0, None)], -1)
exon_intron = (kind != REGION) & (((lt == 2) & (rt == 1)) | ((lt == 1) & (rt == 2)))
ei_nodes = np.where(exon_intron)[0]
exon_flank = {int(b): (int(left[b]) if lt[b] == 2 else int(right[b])) for b in ei_nodes}
intron_flank = {int(b): (int(left[b]) if lt[b] == 1 else int(right[b])) for b in ei_nodes}

_EDGE_KEYS = ("rho_g", "rho_pos", "rho_neg", "rna_src", "mat_abs", "egd", "erd", "md", "sm", "mode_g", "prec_g")

# ---- per-scenario: run calibrate, collect node arrays + the ei-boundary incoming-edge observables ----
per = {}
for cond in conds:
    inp = _scan_and_truth(suite, cond, index, cfg, Path("/tmp/rigel_selfsolve"), suite / "_selfsolve_cache")
    dbg = {}; calmod.calibrate(inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
                               np.asarray(inp["rna_fl_pmf"]), dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg)
    cap = dbg["capture"]
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain); G, R = Gp + Gn, Rp + Rn
    fo = np.where(G + R > _EPS, G / np.maximum(G + R, _EPS), np.nan)
    # index the exon/intron incoming edge per ei-boundary from _edge_modes (dst == boundary, src == the flank)
    ex_edge = {}; in_edge = {}
    for e in cap.get("_edge_modes", []):
        d = e["dst"]
        if d in exon_flank:
            if e["src"] == exon_flank[d]:
                ex_edge[d] = e
            elif e["src"] == intron_flank[d]:
                in_edge[d] = e
    per[cond] = dict(fg=np.asarray(cap["f_g"]), fo=fo, mass=np.asarray(cap["mass_global"]),
                     spl=np.asarray(cap["spl_l"]) + np.asarray(cap["spl_r"]), ex_edge=ex_edge, in_edge=in_edge)

# ---- SJ set + worst scenarios + panel ----
has_spl = np.zeros(len(kind), bool)
for c in conds:
    has_spl |= per[c]["spl"] > _EPS
sj_ei = exon_intron & has_spl
print(f"{len(conds)} scenarios | {int(sj_ei.sum())} splice-junction exon-intron boundary nodes\n")

def _scen_err(c):
    p = per[c]; m = sj_ei & np.isfinite(p["fo"]) & (p["mass"] > _EPS)
    return float(np.average(np.abs(p["fg"][m] - p["fo"][m]), weights=p["mass"][m])) if m.sum() else 0.0
scen_err = sorted(((_scen_err(c), c) for c in conds), reverse=True)
print("worst scenarios by SJ exon-intron boundary error:")
for e, c in scen_err[:6]:
    print(f"  {e:.3f}  {c}")
worst = [c for _, c in scen_err[:N_WORST]]
panel = set()
for c in worst:
    p = per[c]; m = sj_ei & np.isfinite(p["fo"]) & (p["mass"] > _EPS)
    idx = np.where(m)[0]
    panel.update(int(i) for i in idx[np.argsort(-np.abs(p["fg"][idx] - p["fo"][idx]))[:K_PER_SCENARIO]])
panel = sorted(panel)
exon = [exon_flank[b] for b in panel]; intron = [intron_flank[b] for b in panel]
print(f"\nPANEL: {len(panel)} test boundaries\n")

# ---- assemble arrays [n_panel, n_scen] ----
def _node(nodes, key):
    return np.array([[per[c][key][n] for c in conds] for n in nodes], float)
def _edge(side, key):
    out = np.full((len(panel), len(conds)), np.nan)
    for bi, b in enumerate(panel):
        for ci, c in enumerate(conds):
            e = per[c][side].get(b)
            if e is not None:
                out[bi, ci] = e[key]
    return out
save = dict(conds=np.array(conds), worst=np.array(worst),
            bnd=np.array(panel), exon=np.array(exon), intron=np.array(intron),
            bnd_fg=_node(panel, "fg"), bnd_fo=_node(panel, "fo"), bnd_mass=_node(panel, "mass"), bnd_spl=_node(panel, "spl"),
            exon_fg=_node(exon, "fg"), exon_fo=_node(exon, "fo"),
            intron_fg=_node(intron, "fg"), intron_fo=_node(intron, "fo"))
for k in _EDGE_KEYS:
    save["ex_" + k] = _edge("ex_edge", k)
    save["in_" + k] = _edge("in_edge", k)
save["ex_rho_r"] = save["ex_rho_pos"] + save["ex_rho_neg"]
save["in_rho_r"] = save["in_rho_pos"] + save["in_rho_neg"]
SCRATCH.mkdir(parents=True, exist_ok=True)
np.savez(SCRATCH / "boundary_panel.npz", **save)
print(f"saved → {SCRATCH / 'boundary_panel.npz'}  ({len(panel)} triplets × {len(conds)} scenarios)")

# ---- readout: the mature-contamination in the worst scenario ----
wi = list(conds).index(worst[0])
print(f"\nworst scenario ({worst[0]}) — exon→boundary edge (mature contamination):")
print(f"  {'bnd':>5} | {'bnd_fo':>6} {'bnd_fg':>6} | {'ex_rho_g':>8} {'ex_rna_src':>10} {'ex_mat_abs':>10} | {'md':>7} {'ex_egd':>7}")
for bi, b in enumerate(panel[:10]):
    print(f"  {b:>5} | {save['bnd_fo'][bi,wi]:>6.2f} {save['bnd_fg'][bi,wi]:>6.2f} | "
          f"{save['ex_rho_g'][bi,wi]:>8.3f} {save['ex_rna_src'][bi,wi]:>10.3f} {save['ex_mat_abs'][bi,wi]:>10.3f} | "
          f"{save['ex_md'][bi,wi]:>7.1f} {save['ex_egd'][bi,wi]:>7.1f}")
