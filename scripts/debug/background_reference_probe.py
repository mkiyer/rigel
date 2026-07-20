"""BACKGROUND-REFERENCE PROBE — validate the two claims under the aggregate-DNA-background principle:

 (1) DISCRIMINATOR: the pure-DNA (intergenic + RNA-free intron) AGGREGATE density ρ_bg = Σg/ΣE cleanly
     separates the library DNA level (gdna_none ≈ 0 ; gdna100 ; gdna300 high) — strand-independently — where the
     per-node NPMLE (fit on total density) cannot.
 (2) RESOLUTION: a SINGLE pure-DNA region cannot measure the faint background (its count is 0/1, info ρE<1);
     only the pooled aggregate (Σg over the giant support ΣE) resolves it. We show the per-region count is
     mostly 0 while Σg is large enough to pin ρ_bg precisely (rel. Poisson error 1/√Σg).

Uses the ORACLE to define "pure DNA" (R≈0) so the validation is clean; the SOLVER would use the signature
(intergenic/intron) to pick the same regions without truth.

    OMP_NUM_THREADS=1 python background_reference_probe.py
"""

from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")

import numpy as np
import sys

sys.path.insert(0, str(Path("/Users/mkiyer/proj/rigel/scripts/debug")))
from selfsolve_diag import _scan_and_truth  # noqa: E402
import refit_loop_study as R  # noqa: E402
from rigel.calibration.background_reference import measure_background  # noqa: E402
from rigel.calibration.effective_length import region_eff_length  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.substrate import CalibrationSubstrate  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

_EPS = 1e-12


def _production_rho_bg(inp, index):
    """The PRODUCTION `measure_background` scalar (intergenic-only, signature-based) — to compare against the
    oracle-based ρ_bg the probe computes from R<1 nodes."""
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    sub = CalibrationSubstrate.from_payload(inp["payload"], ra)
    reg_el = region_eff_length(ra.region_size_bp, np.asarray(inp["gdna_fl_pmf"]))
    bg = measure_background(sub, ra, reg_el)
    rho = 0.0 if not np.isfinite(bg.log_rho_bg) else float(np.exp(bg.log_rho_bg))
    return rho, bg.n_counts


def main():
    suite = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
    index = TranscriptIndex.load(str(suite / "rigel_index"))
    cfg = PipelineConfig()
    work = Path(os.environ.get("RIGEL_SCRATCH", "/tmp")) / "rigel_selfsolve"
    cache = suite / "_selfsolve_cache"
    conds = sorted(p.name for p in suite.iterdir() if p.is_dir() and p.name.startswith("gdna_"))
    print(f"{'condition':40s} {'ρ_bg(oracle)':>12} {'ρ_bg(PROD)':>11} {'Σg_prod':>9} "
          f"{'ρ_tot':>9} {'%reg=0':>7} {'relerr':>7}")
    for c in conds:
        inp = _scan_and_truth(suite, c, index, cfg, work, cache)
        s = R._setup(inp, index, cfg)
        M, E, G, Rn = s["mass_g"], s["eff_g"], s["G"], s["R"]
        live = (E > 1e-9 * 1.001) & (M > _EPS)
        pure = live & (Rn < 1.0)  # oracle-pure DNA regions (no RNA)
        # (1) aggregate DNA background over the pooled pure-DNA support
        Sg = float(G[pure].sum())
        SE = float(E[pure].sum())
        rho_bg = Sg / max(SE, _EPS)
        rho_tot = float(M[live].sum()) / max(float(E[live].sum()), _EPS)
        bg_fg = rho_bg / max(rho_tot, _EPS)  # background DNA fraction the prior should default to
        # (2) resolution: per pure region, is the count resolvable? Σg precision.
        gpure = G[pure]
        pct0 = float((gpure < 0.5).mean()) * 100.0 if gpure.size else 100.0
        relerr = 1.0 / np.sqrt(max(Sg, 1.0))  # Poisson rel. error of the aggregate scalar
        # the PRODUCTION measure_background scalar (intergenic-only) vs the oracle-based ρ_bg above
        prod_rho, prod_sg = _production_rho_bg(inp, index)
        print(f"{c.replace('gdna_',''):40s} {rho_bg:>12.2e} {prod_rho:>11.2e} {prod_sg:>9.0f} "
              f"{rho_tot:>9.2e} {pct0:>6.1f}% {relerr:>7.3f}", flush=True)


if __name__ == "__main__":
    main()
