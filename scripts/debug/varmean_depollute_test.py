"""Does dropping degenerate short-exon contained densities (E_cont < gDNA fragment length) from the gDNA
var~mean DEFLATE σ²_g and restore imputation precision N for the legitimate (long-exon) edges?

Captures the seed arrays _fit_seed_varmean trains on, then refits σ²_g three ways:
  (a) ALL seed edges (production)
  (b) drop edges where either endpoint is a REGION with eff < frag_len (degenerate short-exon contained M/E)
  (c) drop, AND for kept region endpoints below frag_len cap nothing — pure exclusion
and reports σ²_g(μ) + implied N=μ²/σ²_g at μ∈{2,5,8} (the enriched-exon range). If (b) ≪ (a), the short-exon
degeneracy is the imputation-precision killer and the fix is to exclude sub-fragment-length contained densities.

    OMP_NUM_THREADS=1 python scripts/debug/varmean_depollute_test.py [condition]
"""
from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dissect_regions import build_or_load_cache  # noqa: E402

import rigel.calibration.bp_solver as bp  # noqa: E402
from rigel.calibration.calibrate import calibrate  # noqa: E402
from rigel.calibration.node_chain import REGION  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.variance_model import MonotoneVarMean  # noqa: E402
from rigel.config import CalibrationConfig  # noqa: E402

_EPS = 1e-9
COND = sys.argv[1] if len(sys.argv) > 1 else "gdna_gdna300_ss_0.99_nrna_none_capture_on"
_CAP = {}


def main():
    index, blob = build_or_load_cache(COND, False)
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    pmf = np.asarray(blob["gdna_pmf"], float)
    flen = float(np.sum(np.arange(pmf.size) * pmf) / max(pmf.sum(), _EPS))

    orig = bp._fit_seed_varmean

    def cap(chain, dens, eff, is_seed, seed_w):
        _CAP.update(chain=chain, dens=np.asarray(dens, float), eff=np.asarray(eff, float),
                    is_seed=np.asarray(is_seed, bool), seed_w=np.asarray(seed_w, float))
        return orig(chain, dens, eff, is_seed, seed_w)

    bp._fit_seed_varmean = cap
    calibrate(payload=blob["payload_full"], region_arrays=ra, strand_model=blob["strand_full"],
              gdna_fl_pmf=blob["gdna_pmf"], rna_fl_pmf=blob["rna_pmf"], config=CalibrationConfig())
    bp._fit_seed_varmean = orig

    ch = _CAP["chain"]
    dens, eff, is_seed, seed_w = _CAP["dens"], _CAP["eff"], _CAP["is_seed"], _CAP["seed_w"]
    is_reg = np.asarray(ch.kind) == REGION
    # a node carries a DEGENERATE contained density if it is a region whose gDNA eff < frag length
    degenerate = is_reg & (eff < flen)

    def fit(drop_degenerate):
        left, right = np.asarray(ch.left), np.asarray(ch.right)
        m, r, o, w = [], [], [], []
        for nbr in (left, right):
            idx = np.where((nbr >= 0) & is_seed)[0]
            s = nbr[idx]
            keep = is_seed[s]
            idx, s = idx[keep], s[keep]
            dr, sr, de, se = dens[idx], dens[s], eff[idx], eff[s]
            ok = (dr > 0) & (sr > 0)
            if drop_degenerate:
                ok &= ~degenerate[idx] & ~degenerate[s]
            m.append(0.5 * (dr[ok] + sr[ok]))
            r.append((dr[ok] - sr[ok]) ** 2)
            o.append(dr[ok] / de[ok] + sr[ok] / se[ok])
            w.append(np.minimum(seed_w[idx][ok], seed_w[s][ok]))
        c = lambda p: np.concatenate(p) if p else np.zeros(0)  # noqa: E731
        return MonotoneVarMean.fit_offset(c(m), c(r), c(o), c(w)), int(c(m).size)

    f_all, n_all = fit(False)
    f_clean, n_clean = fit(True)
    n_deg = int((degenerate & is_seed).sum())
    print(f"=== {COND}: var~mean de-pollution test (gDNA frag len ≈ {flen:.0f} bp) ===")
    print(f"seed nodes with degenerate contained density (region, eff<{flen:.0f}): {n_deg}")
    print(f"var~mean edges: all={n_all}  clean(no degenerate endpoint)={n_clean} "
          f"(dropped {n_all-n_clean})\n")
    print(f"{'μ':>5} | {'σ²_g ALL':>10} {'N_all':>7} | {'σ²_g CLEAN':>11} {'N_clean':>8} | {'N gain':>7}")
    print("-" * 60)
    for m in (1.0, 2.0, 5.0, 8.0, 12.0):
        sa = float(f_all.predict(np.array([m]))[0])
        sc = float(f_clean.predict(np.array([m]))[0])
        na, nc = m * m / max(sa, _EPS), m * m / max(sc, _EPS)
        print(f"{m:>5.1f} | {sa:>10.3f} {na:>7.2f} | {sc:>11.3f} {nc:>8.2f} | {nc/max(na,_EPS):>6.1f}×")


if __name__ == "__main__":
    main()
