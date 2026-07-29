"""Precompute the per-node SUBSTRATE for the whole scenario battery (pass-0 once) and pickle it, so landscape
exploration (workflow agents) is pure numpy with NO calibrate re-runs. Includes BOTH region and BOUNDARY nodes.
Run once:  OMP_NUM_THREADS=1 python scripts/debug/gdna_cache_build.py"""
from __future__ import annotations

import dataclasses
import os
import pickle
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np  # noqa: E402

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from selfsolve_diag import _scan_and_truth  # noqa: E402
from flagship_interrogate import _oracle_per_node  # noqa: E402
from rigel.calibration import calibrate  # noqa: E402
from rigel.calibration.calibrate import _fit_gdna_hyperprior  # noqa: E402
from rigel.calibration.bp_solver import REGION  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.signature import RegionType, coarse_type_array  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

import argparse  # noqa: E402

_SCR = "/Users/mkiyer/proj/rigel/scratchpad"  # durable; NOT a session scratchpad (see gdna_explore_lib._CACHE)
_ap = argparse.ArgumentParser()
_ap.add_argument("--suite", default="/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
_ap.add_argument("--out", default=f"{_SCR}/gdna_substrate_cache.pkl")
_args = _ap.parse_args()
OUT = Path(_args.out)
suite = Path(_args.suite)
index = TranscriptIndex.load(str(suite / "rigel_index"))
ra = RegionArrays.from_index(index)
cfg = PipelineConfig()
work = Path("/tmp/rigel_selfsolve")
cache = suite / "_selfsolve_cache"
conds = sorted(d.name for d in suite.iterdir()
               if (d / "sim_oracle.bam").exists() and d.name.startswith("gdna_"))
rtype_all = coarse_type_array(np.asarray(ra.signature)).astype(np.int64)


def _group(cond: str) -> tuple[str, str, str, str]:
    """(capture, gDNA level, strand-specificity, nascent) parsed POSITIONALLY from
    `gdna_<dna>_ss_<ss>_nrna_<nrna>_capture_<cap>`.

    ⚠ Two substring-matching bugs lived here until 2026-07-26, and both corrupted the strata every
    downstream metric is broken out by:
      * `"none" in cond` matched `nrna_none`, so 6 real-gDNA conditions (e.g. gdna100 … nrna_none …) were
        bucketed as the ZERO-gDNA level. The `fabrication` specificity canary selects on exactly that
        bucket, so it was scoring gDNA-bearing scenarios.
      * `capture_verystrong` fell through `"capture_on" in cond` to "OFF" — the STRONGEST capture arm was
        labelled capture-off. Now its own level, so it never silently joins either.
    """
    p = cond.split("_")           # gdna | <dna> | ss | <ss> | nrna | <nrna> | capture | <cap>
    dna, ss, nrna, cap = p[1], p[3], p[5], p[7]
    return ({"on": "ON", "off": "OFF", "verystrong": "VSTRONG"}.get(cap, cap.upper()), dna, ss,
            "none" if nrna == "none" else "nrna")


_GRID = np.linspace(-5.0, 2.5, 260)          # log10 rho_g — must match gdna_explore_lib.GRID
_LN10 = np.log(10.0)


def _production_prior(dbg, chain, mass_global, eff_global):
    """Run the PRODUCTION hyperprior fit on the pass-0 belief — exactly what `calibrate` does at the top of
    its Phase-2 refit loop — and render it on the shared log10 grid so it can be plotted against the oracle.

    Returns (density_on_GRID, n_train, selected_node_mask) or (None, 0, mask) if the fit declines.

    ⚠ Updated 2026-07-27 for **W4**: the production hyperprior is now `GdnaLandscape` and
    `_fit_gdna_hyperprior` is substrate selection only — its signature lost `background`/`bandwidth`
    /`additive` (the last was deleted from the config) and gained the `strength` temperature. The substrate
    predicate below mirrors the shipped one: REGION nodes, AMBIG excluded, boundaries excluded, plus the
    zero-count structural anchor on non-exon regions."""
    prior = _fit_gdna_hyperprior(
        chain, dbg["belief"], dbg["statics"], ra, mass_global, eff_global,
        strength=cfg.calibration.gdna_prior_strength,
    )
    # reproduce the substrate mask the fit selected, for the plots (same predicate, kept in sync by eye)
    isr = np.asarray(chain.kind) == REGION
    fp = np.asarray(dbg["statics"].free_pos, bool)
    fn = np.asarray(dbg["statics"].free_neg, bool)
    rtype = coarse_type_array(np.asarray(ra.signature))
    ridx = np.clip(np.asarray(chain.ref_idx, dtype=np.int64), 0, rtype.shape[0] - 1)
    expressed = isr & (eff_global > 1.0e-9) & (mass_global > 1.0e-12)
    anchor = isr & (eff_global > 1.0e-9) & (mass_global <= 1.0e-12) & (rtype[ridx] != RegionType.EXON)
    sel = (expressed & ((fp ^ fn) | (~fp & ~fn))) | anchor
    if prior is None:
        return None, 0, sel
    # the landscape carries logP on its OWN natural-log grid; interpolate onto the shared log10 grid.
    d = np.exp(np.interp(_GRID * _LN10, prior.log_rho, prior.logP,
                         left=prior.logP[0], right=prior.logP[-1]))
    tot = float(d.sum())
    return (d / tot if tot > 0 else None), int(prior.n_train), sel


scen = []
for cond in conds:
    inp = _scan_and_truth(suite, cond, index, cfg, work, cache)
    dbg = {}
    calibrate(inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
              np.asarray(inp["rna_fl_pmf"]),
              dataclasses.replace(cfg.calibration, calib_refit_iters=0), _debug=dbg)
    chain = dbg["chain"]
    cap = dbg["capture"]
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
    isr = np.asarray(chain.kind) == REGION
    idx = np.asarray(chain.ref_idx, np.int64)
    ntype = np.where(isr, rtype_all[np.clip(idx, 0, len(rtype_all) - 1)], 3)  # 0/1/2 region-type, 3 boundary
    f0 = np.asarray(dbg["belief"].f_g, float)
    mg = np.asarray(cap["mass_global"], np.float64)
    eg = np.asarray(cap["eff_global"], np.float64)
    prod_P, prod_cells, prod_sel = _production_prior(dbg, chain, mg, eg)
    scen.append(dict(
        cond=cond,
        group=_group(cond),
        is_region=isr,
        ntype=ntype.astype(np.int8),
        fp=np.asarray(cap["free_pos"], bool),
        fn=np.asarray(cap["free_neg"], bool),
        mass=mg,
        eff=eg,
        f0=f0,
        g_hat=(f0 * mg),
        var=np.asarray(dbg["belief"].var_gdna, np.float64),
        G=(Gp + Gn).astype(np.float64),
        R=(Rp + Rn).astype(np.float64),
        # the PRODUCTION hyperprior fit (DensityNPMLE), rendered on the shared log10 grid + its substrate
        prod_P=prod_P,
        prod_cells=prod_cells,
        prod_sel=prod_sel,
    ))
    print(f"  cached {cond}: {isr.sum()} region + {(~isr).sum()} boundary nodes; "
          f"production prior {'OK' if prod_P is not None else 'DECLINED'} "
          f"({prod_cells} cells, {int(prod_sel.sum())} training nodes)")

OUT.write_bytes(pickle.dumps(scen))
print(f"\nwrote {len(scen)} scenarios -> {OUT}  ({OUT.stat().st_size / 1e6:.1f} MB)")
