"""Pass-0 MESSAGE-PROPAGATION audit on the grounded full-transcript toy (population priors INJECTED, so κ /
NPMLE / ρ_bg are trustworthy; calib_refit_iters=0 — pass-0 ONLY, no Phase-2).

Goal (owner): find the ARITHMETIC (mode-determination) bugs in message propagation under HYBRID CAPTURE, where
total counts swing hard across depleted↔enriched (intron↔exon) transitions.

Reports:
  A. PER-NODE: oracle f_g, local (message-free) f_g, final f_g — does the message move it TOWARD or AWAY from
     oracle? (the sign of the message's net effect)
  B. PER-EDGE across the cliff: src→dst (type, density), the imputed gDNA message mode, its precision, and the
     σ²_transfer damping — localizing WHICH edge's mode arithmetic is wrong.
"""
from __future__ import annotations
import argparse
import contextlib
import dataclasses
import io
import numpy as np

import toy_inject as ti
from rigel.calibration.calibrate import calibrate
from rigel.calibration.node_chain import REGION
from rigel.calibration.signature import coarse_type_array
from rigel.config import CalibrationConfig

_EPS = 1e-9
_TN = {0: "intgc", 1: "intron", 2: "exon"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cond", default="gdna_gdna300_ss_0.50_nrna_present_capture_on")
    ap.add_argument("--n-rna", type=int, default=20000)
    ap.add_argument("--gdna", type=float, default=3.0)
    ap.add_argument("--nascent", type=float, default=200.0)
    ap.add_argument("--mature", type=float, default=100.0)
    args = ap.parse_args()
    capture_on = "capture_on" in args.cond

    priors = ti.extract_priors(args.cond)
    exons = [(6000, 7000), (11000, 12000), (16000, 17000)]
    with contextlib.redirect_stdout(io.StringIO()):
        toy = ti.build_toy("audit", exons=exons, gdna_fraction=args.gdna, nascent=args.nascent,
                           mature=args.mature, capture_on=capture_on, n_rna=args.n_rna, genome_length=30000)

    cfg = dataclasses.replace(CalibrationConfig(), calib_refit_iters=0)  # PASS-0 ONLY
    dbg = {}
    calibrate(toy["payload"], toy["ra"], toy["strand_model"], toy["gdna_fl_pmf"], toy["rna_fl_pmf"],
              cfg, _debug=dbg, injected_priors=priors)
    cap = dbg["capture"]
    chain = toy["chain"]; ra = toy["ra"]; truth = toy["truth"]
    kind = np.asarray(chain.kind); ref_idx = np.asarray(chain.ref_idx)
    ctype = coarse_type_array(np.asarray(ra.signature)).astype(int)
    fg = np.asarray(cap["f_g"]); fgl = np.asarray(cap["fg_loc"])
    mass = np.asarray(cap["mass_global"]); eff = np.asarray(cap["eff_global"])
    _us = cap["_uni_static"]
    mu_p = np.log(np.maximum(np.asarray(_us["rho_node0"]), 1e-12)); var_p = np.asarray(_us["logvar_tot"])
    dens = np.where(eff > _EPS, mass / np.maximum(eff, _EPS), 0.0)

    def ntype(n):
        return _TN.get(int(ctype[ref_idx[n]]), "?") if kind[n] == REGION else "BND"

    def orac(n):
        return truth[int(ref_idx[n])] if kind[n] == REGION else np.nan

    print(f"=== {args.cond}  (pass-0 only, injected priors, capture={'ON' if capture_on else 'off'}) ===")
    print(f"    κ={priors.rna_sense_frac:.3f}  NPMLE_cells={priors.enrichment_prior.n_cells}\n")

    print("A. PER-NODE — does the message move f_g toward oracle?")
    print(f"  {'node':>4} {'type':>6} | {'dens':>9} {'muProj':>7} | {'oracle':>6} {'local':>6} {'final':>6} "
          f"| {'msgΔ':>7} {'toward?':>8}")
    for n in range(chain.n_nodes):
        o = orac(n)
        loc, fin = fgl[n], fg[n]
        md = fin - loc
        toward = ""
        if not np.isnan(o):
            toward = "TOWARD" if abs(fin - o) < abs(loc - o) - 1e-4 else ("away" if abs(fin - o) > abs(loc - o) + 1e-4 else "-")
        os_ = "  nan" if np.isnan(o) else f"{o:.3f}"
        print(f"  {n:>4} {ntype(n):>6} | {dens[n]:>9.4f} {mu_p[n]:>7.2f} | {os_:>6} {loc:>6.3f} {fin:>6.3f} "
              f"| {md:>+7.3f} {toward:>8}")

    # B. PER-EDGE across the cliff
    em = cap.get("_edge_modes", [])
    print("\nB. PER-EDGE gDNA message (src→dst): the imputed mode + precision + σ²_transfer across the cliff")
    print(f"  {'src':>4}{'':>7}->{'dst':>4}{'':>7} | {'ρg_src':>8} {'mode_g':>7} {'prec_g':>7} | "
          f"{'σ²T':>7} {'shift?':>6} | {'src_or':>6} {'dst_or':>6}")
    for e in em:
        s, d = e["src"], e["dst"]
        s2t = var_p[d] + var_p[s]  # M5: sigma^2_transfer = Var(log r), the derived law
        so, do = orac(s), orac(d)
        print(f"  {s:>4}{('(' + ntype(s) + ')'):>7}->{d:>4}{('(' + ntype(d) + ')'):>7} | {e['rho_g']:>8.4f} "
              f"{e['mode_g']:>7.2f} {e['prec_g']:>7.2f} | {s2t:>7.3f} {str(e['use_shift']):>6} | "
              f"{('nan' if np.isnan(so) else f'{so:.3f}'):>6} {('nan' if np.isnan(do) else f'{do:.3f}'):>6}")


if __name__ == "__main__":
    main()
