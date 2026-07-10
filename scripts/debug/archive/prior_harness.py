"""Shared harness for prototyping the density-prior INTEGRATION (the "two likelihood landscapes" problem).

A candidate is a pure function:

    integrate(strand_ll, fg_grid, d_tot, kde) -> f_g estimate (scalar in [0,1])

  strand_ll : np.ndarray[K]  — the node's strand LOG-likelihood landscape over fg_grid (tau-marginalised,
              tilt-constrained). This is landscape #1. NOT normalised.
  fg_grid   : np.ndarray[K]  — the f_g grid (linear in (0,1)).
  d_tot     : float          — the node's total density M / E_gdna (so gDNA density at f_g is rho_g = f_g*d_tot).
  kde       : GdnaDensityPrior with .logpdf(log_rho)->log P(log rho_g) and .train_x/.bandwidth/.modes.
              This is landscape #2 (the population prior over log gDNA density).

`evaluate(integrate)` runs the candidate on every region node of two conditions and returns metrics:
  - gdna300 AMBIG-exon accuracy (mass-weighted f_g vs oracle 0.757; the CRUSH test)
  - gdna300 single-strand exon/intron error (must not regress)
  - gdna_none FALSE-POSITIVE gDNA (must stay ~0; the self-pin/FP test) — stranded AND unstranded
  - cliff sweep: flat-strand f_g vs d_tot (must be smooth/monotone, no crush of mid-density nodes)

Build once (slow: fits KDE per condition), then evaluate() is fast.

    from prior_harness import load, evaluate
    def integrate(strand_ll, fg_grid, d_tot, kde): ...
    print(evaluate(integrate))
"""
from __future__ import annotations
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys, pickle
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

_PKL = Path("/tmp/dissect_cache/prior_harness.pkl")
_CONDS = {
    "g300_s99": "gdna_gdna300_ss_0.99_nrna_none_capture_on",
    "g300_s50": "gdna_gdna300_ss_0.50_nrna_none_capture_on",
    "none_s99": "gdna_none_ss_0.99_nrna_none_capture_on",
    "none_s50": "gdna_none_ss_0.50_nrna_none_capture_on",
}
_KGRID = 161  # f_g grid resolution
_NTILT = 61


def _strand_landscape(upos, uneg, fg_grid, kappa, od_g, od_r, free_pos, free_neg):
    """tau-marginal strand LOG-likelihood over fg_grid for one node (landscape #1)."""
    from rigel.calibration.simplex import _mixture_strand_loglik
    N = upos + uneg
    K = fg_grid.shape[0]
    if free_pos and free_neg:                     # AMBIG: marginalise tau in [-1,1]
        tau = np.linspace(-1.0, 1.0, _NTILT)
        fgc = fg_grid[:, None]
        fpos = (1.0 - fgc) * (1.0 + tau[None, :]) / 2.0
        fneg = (1.0 - fgc) * (1.0 - tau[None, :]) / 2.0
        psi = _mixture_strand_loglik(np.array([[upos]]), N, fg_grid[None, :, None],
                                     fpos[None], fneg[None], kappa, od_g, od_r)[0]  # (K, Kt)
        m = psi.max(axis=1)
        return m + np.log(np.exp(psi - m[:, None]).sum(axis=1))
    # single-strand: all RNA on the free strand
    fpos = (1.0 - fg_grid) if free_pos else np.zeros(K)
    fneg = (1.0 - fg_grid) if free_neg else np.zeros(K)
    return _mixture_strand_loglik(np.array([[upos]]), N, fg_grid[None, :],
                                  fpos[None], fneg[None], kappa, od_g, od_r)[0]


def build():
    from dissect_regions import build_or_load_cache
    from rigel.calibration.region_arrays import RegionArrays
    from rigel.calibration.substrate import CalibrationSubstrate
    from rigel.calibration.signature import coarse_strand_from_signature, coarse_type_from_signature
    import rigel.calibration.gdna_density_prior as gdp
    fg_grid = np.linspace(1e-3, 1 - 1e-3, _KGRID)
    out = {"fg_grid": fg_grid, "conds": {}}
    for key, cond in _CONDS.items():
        index, blob = build_or_load_cache(cond, False)
        ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
        # capture the fitted KDE + geometry (mass_global/eff_global) from a real calibrate run
        CAP = {}
        _ofit = gdp.GdnaDensityPrior.fit.__func__
        def fit(cls, *a, **k):
            p = _ofit(cls, *a, **k); CAP["kde"] = p; return p
        gdp.GdnaDensityPrior.fit = classmethod(fit)
        import rigel.calibration.bp_solver as bp
        _ons = bp.node_sweep
        def w(*a, **k):
            k["_capture"] = CAP; return _ons(*a, **k)
        bp.node_sweep = w
        sys.modules["rigel.calibration.calibrate"].node_sweep = w
        from rigel.calibration.calibrate import calibrate
        from rigel.config import CalibrationConfig
        cal = calibrate(payload=blob["payload_full"], region_arrays=ra, strand_model=blob["strand_full"],
                        gdna_fl_pmf=blob["gdna_pmf"], rna_fl_pmf=blob["rna_pmf"], config=CalibrationConfig())
        bp.node_sweep = _ons; gdp.GdnaDensityPrior.fit = classmethod(_ofit)
        kde = CAP["kde"]
        # region-node arrays (the chain interleaves; _capture arrays are node-indexed — pull REGION nodes)
        from rigel.calibration.node_chain import REGION
        # rebuild chain to map region->node (same as calibrate)
        from rigel.calibration.node_chain import build_node_chain
        chain = build_node_chain(blob["payload_full"].ref_region_offsets, blob["payload_full"].ref_boundary_offsets)
        kind = np.asarray(chain.kind); rref = np.asarray(chain.ref_idx)
        r2n = {int(rref[i]): i for i in range(len(kind)) if kind[i] == REGION}
        mass_g = np.asarray(CAP["mass_global"]); eff_g = np.asarray(CAP["eff_global"])
        fp = np.asarray(CAP["free_pos"], bool); fn = np.asarray(CAP["free_neg"], bool)
        sf = CalibrationSubstrate.from_payload(blob["payload_full"], ra)
        sg = CalibrationSubstrate.from_payload(blob["payload_gdna"], ra)
        sr = CalibrationSubstrate.from_payload(blob["payload_rna"], ra)
        upos = np.asarray(sf.contained.n_unspliced_pos, float); uneg = np.asarray(sf.contained.n_unspliced_neg, float)
        go = np.asarray(sg.contained.mass_unspliced, float); ro = np.asarray(sr.contained.mass_unspliced, float)
        sig = index.region_df["signature"].to_numpy()
        scls = np.array([coarse_strand_from_signature(int(s)) for s in sig])
        tcls = np.array([coarse_type_from_signature(int(s)) for s in sig])
        kappa = float(cal.rna_sense_frac); od_g = float(cal.gdna_strand_overdispersion); od_r = float(cal.rna_strand_overdispersion)
        R = len(sig)
        rows = []
        for r in range(R):
            n = r2n.get(r)
            if n is None: continue
            m = mass_g[n]; e = max(eff_g[n], 1e-9)
            if m <= 0: continue
            ll = _strand_landscape(upos[r], uneg[r], fg_grid, kappa, od_g, od_r, bool(fp[n]), bool(fn[n]))
            oracle_fg = go[r] / (go[r] + ro[r] + 1e-9) if (go[r] + ro[r]) > 0 else np.nan
            rows.append(dict(reg=r, node=n, upos=upos[r], uneg=uneg[r], mass=m, eff=e, d_tot=m / e,
                             strand_ll=ll.astype(np.float32), oracle_fg=oracle_fg,
                             scls=int(scls[r]), tcls=int(tcls[r]),
                             ambig=bool(fp[n] and fn[n]), mass_gdna_or=go[r], mass_rna_or=ro[r]))
        out["conds"][key] = dict(rows=rows, kappa=kappa, od_g=od_g, od_r=od_r,
                                 kde=dict(x_grid=kde.x_grid, logP_grid=kde.logP_grid, bandwidth=kde.bandwidth,
                                          train_x=kde.train_x, modes=kde.modes))
    _PKL.parent.mkdir(parents=True, exist_ok=True)
    with open(_PKL, "wb") as f:
        pickle.dump(out, f)
    return out


class _KDE:
    """Lightweight KDE view (logpdf + fields) reconstructed from the pickle."""
    def __init__(self, d):
        self.x_grid = d["x_grid"]; self.logP_grid = d["logP_grid"]
        self.bandwidth = float(d["bandwidth"]); self.train_x = d["train_x"]; self.modes = d["modes"]
    def logpdf(self, log_rho):
        x = np.asarray(log_rho, float)
        return np.interp(x, self.x_grid, self.logP_grid, left=float(self.logP_grid[0]), right=float(self.logP_grid[-1]))


def load():
    if not _PKL.exists():
        build()
    with open(_PKL, "rb") as f:
        return pickle.load(f)


def evaluate(integrate, data=None, verbose=True):
    """Run a candidate integrate() over all conditions; return a metrics dict + (optional) printout."""
    if data is None:
        data = load()
    fg_grid = data["fg_grid"]
    res = {}
    for key, cd in data["conds"].items():
        kde = _KDE(cd["kde"])
        rows = cd["rows"]
        fg_hat = np.array([float(np.clip(integrate(r["strand_ll"].astype(np.float64), fg_grid, r["d_tot"], kde), 0, 1))
                           for r in rows])
        mass = np.array([r["mass"] for r in rows])
        ambig = np.array([r["ambig"] for r in rows])
        tcls = np.array([r["tcls"] for r in rows])
        go = np.array([r["mass_gdna_or"] for r in rows]); ro = np.array([r["mass_rna_or"] for r in rows])
        gpred = fg_hat * mass                       # predicted gDNA count per node
        # metrics
        ae = ambig & (tcls == 2)
        ss_exon = (~ambig) & (tcls == 2)
        m = {}
        if "g300" in key:
            M = mass[ae]
            m["ambig_exon_mw_fg"] = float((M * fg_hat[ae]).sum() / M.sum()) if M.sum() else float("nan")
            m["ambig_exon_oracle"] = float(go[ae].sum() / (go[ae] + ro[ae]).sum()) if (go[ae]+ro[ae]).sum() else float("nan")
            m["contained_gdna_err"] = float((gpred - go).sum())
            m["ss_exon_err"] = float((gpred[ss_exon] - go[ss_exon]).sum())
        else:  # gdna_none: FP gDNA (oracle gDNA = 0)
            m["fp_gdna"] = float(gpred.sum())
        res[key] = m
    # cliff sweep: flat strand, sweep d_tot on the g300_s99 KDE
    kde = _KDE(data["conds"]["g300_s99"]["kde"])
    flat = np.zeros_like(fg_grid)
    sweep = {}
    for ld in [-1, 0, 1, 1.5, 2, 2.5, 3, 3.5]:
        sweep[ld] = float(np.clip(integrate(flat, fg_grid, float(np.exp(ld)), kde), 0, 1))
    res["cliff_sweep_fg_vs_logdtot"] = sweep
    if verbose:
        import json
        print(json.dumps(res, indent=2, default=float))
    return res


if __name__ == "__main__":
    if "--build" in sys.argv:
        build(); print("built", _PKL)
    else:
        # sanity: the raw grid-multiply candidate (known to fix AMBIG, crush boundaries/mid-density)
        def raw(strand_ll, fg_grid, d_tot, kde):
            log_rho = np.log(np.maximum(fg_grid * d_tot, 1e-9))
            post = strand_ll + kde.logpdf(log_rho)
            post = np.exp(post - post.max())
            return (fg_grid * post).sum() / post.sum()   # posterior mean
        evaluate(raw)
