"""Shared offline evaluator for the boundary-mode derivation (loads the enriched panel from boundary_panel.py).

A candidate transfer MODE is a function ``fn(P) -> implied_fg[n_panel, n_scen]`` of the panel observables ``P``
(a dict of arrays). ``score(fn)`` reports, vs the boundary oracle ``bnd_fo``:
  * UNDER-fix : mean implied f_g on panel×worst where bnd_fo>0.8 (must be HIGH — fix the mature-contamination under-call)
  * OVER-guard: mean implied f_g on ALL scenarios where bnd_fo<0.2 (must stay LOW — no gDNA over-call / intron bleed)
  * panel_err : mean |implied − bnd_fo| on panel×worst  (must DROP vs baseline)
  * ctrl_err  : mean |implied − bnd_fo| over ALL scenarios where bnd_fo finite (must not rise much)
The offline evaluator scores a candidate FORMULA against the oracle (derivation), not the full fold; the winning
formula is then implemented in bp_solver and verified with real calibrate (panel + antisense + goldens).

Observables per [n_panel, n_scen]: bnd_fo/bnd_fg, exon_fo/exon_fg, intron_fo/intron_fg, bnd_mass, bnd_spl, and
per incoming edge ex_/in_: rho_g, rho_pos, rho_neg, rho_r, rna_src, mat_abs, egd, erd, md, sm, mode_g, prec_g."""
from __future__ import annotations
from pathlib import Path
import numpy as np

SCRATCH = Path("/private/tmp/claude-503/-Users-mkiyer-proj-rigel/4f7a248b-0c78-4b40-9030-462373aefb19/scratchpad")
_EPS = 1e-9


def load_panel(path=SCRATCH / "boundary_panel.npz"):
    z = np.load(path, allow_pickle=True)
    P = {k: z[k] for k in z.files}
    P["conds"] = [str(c) for c in P["conds"]]
    P["worst"] = [str(c) for c in P["worst"]]
    P["worst_mask"] = np.array([c in set(P["worst"]) for c in P["conds"]])  # [n_scen]
    return P


def clip01(x):
    return np.clip(x, 0.0, 1.0)


def score(fn, P=None, name=""):
    if P is None:
        P = load_panel()
    fg = fn(P)
    fo = P["bnd_fo"]
    fin = np.isfinite(fo) & np.isfinite(fg)
    wm = P["worst_mask"][None, :] & fin
    hi = wm & (fo > 0.8)            # panel×worst, truly-gDNA → the under-call cases
    lo = fin & (fo < 0.2)          # all scenarios, truly-RNA → the over-call / intron-bleed guard
    under_fix = float(np.nanmean(fg[hi])) if hi.any() else float("nan")
    over_guard = float(np.nanmean(fg[lo])) if lo.any() else float("nan")
    panel_err = float(np.nanmean(np.abs(fg[wm] - fo[wm]))) if wm.any() else float("nan")
    ctrl_err = float(np.nanmean(np.abs(fg[fin] - fo[fin]))) if fin.any() else float("nan")
    print(f"{name:>22} | UNDER-fix(→1) {under_fix:>5.2f} | OVER-guard(→0) {over_guard:>5.2f} | "
          f"panel_err {panel_err:>5.2f} | ctrl_err {ctrl_err:>5.2f}")
    return dict(under_fix=under_fix, over_guard=over_guard, panel_err=panel_err, ctrl_err=ctrl_err)


# ---- built-in candidate modes (baseline + the leading physical candidates) ----
def m_current(P):
    """The ACTUAL solved f_g (folds strand+prior+both messages) — the baseline to beat."""
    return P["bnd_fg"]

def m_density_exon(P):
    """Exon gDNA message alone, nascent=RESIDUAL: f_g = ρ_g^exon·E_g^dst / md  (suppress the mature-contaminated exon RNA)."""
    return clip01(P["ex_rho_g"] * P["ex_egd"] / np.maximum(P["ex_md"], _EPS))

def m_intron_gdna(P):
    """Crossing gDNA from the DEPLETED intron flank (background), nascent=residual: f_g = ρ_g^intron·E_g^dst / md."""
    return clip01(P["in_rho_g"] * P["in_egd"] / np.maximum(P["in_md"], _EPS))

def m_max_flank_gdna(P):
    """gDNA density = MAX of exon/intron flank (the crossing sees whichever side carries gDNA), ÷ own md."""
    rg = np.fmax(np.nan_to_num(P["ex_rho_g"], nan=0.0), np.nan_to_num(P["in_rho_g"], nan=0.0))
    egd = np.where(np.isfinite(P["ex_egd"]), P["ex_egd"], P["in_egd"])
    md = np.where(np.isfinite(P["ex_md"]), P["ex_md"], P["in_md"])
    return clip01(rg * egd / np.maximum(md, _EPS))


def m_geomean(P):
    """The SYNTHESIS: geometric-mean cliff-interpolated crossing gDNA (partial-capture crossing sits at the
    log-midpoint of its two flank gDNA densities), density frame in the mature-free crossing (nascent=residual).
    Count-LEGAL: per-node source densities; md only as the mature-free total being split. ZERO free constants.
    Single-flank fallback where one neighbour carries no gDNA (== today's density mode)."""
    ex, inn = np.nan_to_num(P["ex_rho_g"], nan=0.0), np.nan_to_num(P["in_rho_g"], nan=0.0)
    gm = np.sqrt(np.maximum(ex, 0.0) * np.maximum(inn, 0.0))
    rg = np.where(gm > 0.0, gm, np.fmax(ex, inn))  # fallback to the single available flank
    egd = np.where(np.isfinite(P["ex_egd"]), P["ex_egd"], P["in_egd"])
    md = np.where(np.isfinite(P["ex_md"]), P["ex_md"], P["in_md"])
    return clip01(rg * egd / np.maximum(md, _EPS))


CANDIDATES = {"current": m_current, "density_exon(residual)": m_density_exon,
              "intron_gdna": m_intron_gdna, "max_flank_gdna": m_max_flank_gdna,
              "GEOMEAN(synthesis)": m_geomean}

if __name__ == "__main__":
    P = load_panel()
    print(f"panel: {P['bnd'].shape[0]} triplets × {len(P['conds'])} scenarios; worst={P['worst']}\n")
    for nm, fn in CANDIDATES.items():
        score(fn, P, nm)
    # count-zero-info audit: does the own-crossing MASS md vote composition? (the synthesis' rejection of the
    # top scorers). corr(md/egd, oracle f_g) strongly negative ⇒ using md as an anchor is a count-vote.
    fo = P["bnd_fo"]; dcross = P["ex_md"] / np.maximum(P["ex_egd"], _EPS)
    m = np.isfinite(fo) & np.isfinite(dcross)
    print(f"\ncount-vote audit: corr(own-crossing density md/egd, oracle f_g) = "
          f"{np.corrcoef(dcross[m], fo[m])[0,1]:+.2f}  (strongly −ve ⇒ md is a count-vote; geomean avoids it)")
