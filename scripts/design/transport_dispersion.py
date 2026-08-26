"""⭐⭐⭐ WHERE DOES THE FLANK-TRANSPORT DISPERSION COME FROM? — the decomposition, certified truth,
no solver runs.

Per two-complete-flank exon this instrument measures the pair disagreement D = log(rate_L/rate_R)
and the common-mode center c = log(geomean rate / certified contained rate), and charges each
candidate source with its measurable signature: per-route trigamma COUNTING (flank side), the
TRUTH count's own trigamma (⛔ the instrument artifact that inflated every prior "certified
scatter" reading ~sd 0.44 — median deep-pair n_mrna is 11-19), the LENGTH-CURVE center by exon
decile (geometry — both flanks move together), termini strictly INSIDE the exon body (structure),
the credit-leftmost deposit rule (a SIDE-SPLIT center would betray it), a strand-matched
re-pooling arm, and the capture OFF/ON contrast. What survives every charge is the honest
transport premise. Also reports the single-complete-flank population's relative center.

⛔ Fit nothing on shallow pairs: below ~100 flank flux the disagreement is 100 % counting.
⚠ The strand-matched arm keys the substrate's count columns naively per sj strand — measured to
be the WRONG key (columns are genome-strand); its reading is retained as the falsification of
that key, not as a candidate.
"""

import sys
from pathlib import Path

import numpy as np
from scipy.special import polygamma

REPO = Path("/Users/mkiyer/proj/rigel")
sys.path.insert(0, str(REPO / "scripts" / "design"))

from rigel.calibration.fl import build_fl_models  # noqa: E402
from rigel.calibration.gdna_opportunity import gdna_opportunity_from_index  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.region_chain import REGION, build_region_chain  # noqa: E402
from rigel.calibration.region_geometry import (  # noqa: E402
    build_region_geometry,
    build_region_statics,
)
from rigel.calibration.rna_anchor import build_route_table  # noqa: E402
from rigel.calibration.sj_opportunity import crossing_probability_from_index  # noqa: E402
from rigel.calibration.splice_graph import (  # noqa: E402
    build_boundary_flags_array,
    build_sj_geometry_arrays,
)
from rigel.calibration.structural_claims import build_structural_claims  # noqa: E402
from rigel.calibration.substrate import CalibrationSubstrate  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402
from rigel.scan_cache import read_scan_cache  # noqa: E402

RUNS = Path.home() / "Downloads" / "rigel_runs"
SUITE = RUNS / "suite" / "ladder"

DEFAULT_CONDITIONS = [
    "gdna_g00_ss_0.99_nrna_mid_capture_off",
    "gdna_g00_ss_0.99_nrna_mid_capture_on",
    "gdna_g05_ss_0.99_nrna_mid_capture_off",
    "gdna_g50_ss_0.99_nrna_mid_capture_off",
    "gdna_g50_ss_0.99_nrna_mid_capture_on",
    "gdna_g50_ss_0.50_nrna_mid_capture_off",
]

_EPS = 1e-12


def flank(routes_j, flux, opp, col, matched):
    """Route-summed rate + delta-method log-variance for one flank; (rate, vlog, flux_total, k)."""
    js = np.asarray(routes_j, np.int64)
    if js.size == 0:
        return 0.0, np.inf, 0.0, 0
    f = flux[js, col[js]] if matched else flux[js].sum(axis=1)
    a = np.maximum(opp[js], _EPS)
    k = js.size
    r_j = (f + 0.5 / k) / a
    rate = float(r_j.sum())
    # Var(log rate) by delta method: each route's log-rate variance is trigamma(f + 1/2)
    v = float(np.sum(r_j**2 * polygamma(1, f + 0.5)) / max(rate, _EPS) ** 2)
    return rate, v, float(f.sum()), k


def main() -> int:
    import argparse

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--suite", type=Path, default=SUITE)
    ap.add_argument("--index", type=Path, default=RUNS / "suite" / "rigel_index")
    ap.add_argument("--conditions", nargs="*", default=None)
    ap.add_argument("--out-dir", type=Path, default=Path("/tmp/rigel_transport_decomposition"))
    args = ap.parse_args()
    conditions = args.conditions or DEFAULT_CONDITIONS
    args.out_dir.mkdir(parents=True, exist_ok=True)

    index = TranscriptIndex.load(str(args.index))
    ra = RegionArrays.from_frame(index.regions_df, index.ref_name_to_id)
    bflags = build_boundary_flags_array(index)
    sj = build_sj_geometry_arrays(index)

    # mature transcript termini per ref, for the termini-inside classification
    t = index.t_df
    mature = t[~t["is_nrna"].astype(bool) & ~t["is_synthetic"].astype(bool)]
    term_by_ref = {}
    for ref, grp in mature.groupby("ref"):
        rid = index.ref_name_to_id.get(ref)
        if rid is None:
            continue
        term_by_ref[rid] = np.sort(
            np.concatenate([grp["start"].to_numpy(np.int64), grp["end"].to_numpy(np.int64)])
        )

    from rigel.calibration.splice_graph import Strand

    sj_col = np.where(np.asarray(sj.strand) == np.int8(Strand.POS), 0, 1)

    for condition in conditions:
        cache = read_scan_cache(args.suite / "oracle_cache" / condition / "_main", index)
        payload = cache.payload
        fl = build_fl_models(
            payload,
            sj_opportunity=crossing_probability_from_index(index, int(payload.max_length)),
            gdna_opportunity=gdna_opportunity_from_index(index, int(payload.max_length)),
        )
        substrate = CalibrationSubstrate.from_payload(payload, ra)
        chain = build_region_chain(payload.ref_region_offsets, payload.ref_boundary_offsets)
        geometry = build_region_geometry(chain, substrate, ra, sj, fl.gdna_pmf, fl.rna_pmf, None)
        statics = build_region_statics(chain, ra, bflags)
        claims = build_structural_claims(chain, statics)
        routes = build_route_table(sj, substrate, fl.rna_pmf)
        truth = dict(np.load(args.suite / "oracle_cache" / condition / "slot_truth.npz"))

        kind = np.asarray(chain.kind)
        obj = np.asarray(chain.obj_idx, np.int64)
        is_k = np.asarray(truth["kind"]) == REGION
        tobj = np.asarray(truth["obj"], np.int64)[is_k]
        tm = np.zeros(int(tobj.max()) + 1)
        tm[tobj] = np.asarray(truth["n_mrna"], np.float64)[is_k]

        comp_l = np.asarray(claims.exon_flank_left_complete, bool)
        comp_r = np.asarray(claims.exon_flank_right_complete, bool)
        eff_r = np.asarray(geometry.eff_rna, np.float64)
        flux_2col = np.asarray(substrate.sj.count, np.float64)
        opp = np.asarray(routes.opportunity, np.float64)

        rows = []
        for s_ in np.flatnonzero((kind == REGION) & (comp_l | comp_r)):
            r = int(obj[s_])
            jl = routes.into.get(r, [])
            jr = routes.outof.get(r, [])
            has_l = comp_l[s_] and len(jl) > 0
            has_r = comp_r[s_] and len(jr) > 0
            if not (has_l or has_r):
                continue
            rl, vl, fl_tot, kl = flank(jl, flux_2col, opp, sj_col, False) if has_l else (
                np.nan, np.nan, 0.0, 0)
            rr, vr, fr_tot, kr = flank(jr, flux_2col, opp, sj_col, False) if has_r else (
                np.nan, np.nan, 0.0, 0)
            rlm, _, _, _ = flank(jl, flux_2col, opp, sj_col, True) if has_l else (
                np.nan, 0, 0, 0)
            rrm, _, _, _ = flank(jr, flux_2col, opp, sj_col, True) if has_r else (
                np.nan, 0, 0, 0)
            length = int(ra.end[r] - ra.start[r])
            rid = int(ra.ref_id[r])
            terms = term_by_ref.get(rid)
            n_inside = 0
            if terms is not None:
                lo = np.searchsorted(terms, int(ra.start[r]) + 1, side="left")
                hi = np.searchsorted(terms, int(ra.end[r]) - 1, side="right")
                n_inside = int(hi - lo)
            e_r = float(eff_r[s_])
            rho_true = tm[r] / max(e_r, _EPS) if e_r > 0 else np.nan
            rows.append((
                r, length, kl, kr, fl_tot, fr_tot, rl, rr, rlm, rrm, vl, vr,
                rho_true, tm[r], n_inside, int(has_l), int(has_r),
            ))

        arr = np.array(rows, np.float64)
        np.savez_compressed(args.out_dir / f"transport_{condition}.npz", rows=arr)
        (
            r_id, length, kl, kr, flt, frt, rl, rr, rlm, rrm, vl, vr,
            rho_t, n_m, n_in, hl, hr,
        ) = arr.T

        pair = (hl > 0) & (hr > 0) & (rl > 0) & (rr > 0)
        deep = pair & (flt + frt >= 100)
        print(f"\n== {condition}")
        print(f"   flanked exons {arr.shape[0]:,}  pairs {int(pair.sum()):,}  "
              f"deep pairs (flux>=100) {int(deep.sum()):,}")

        def dstats(mask, rL, rR, label):
            d = np.log(rL[mask] / rR[mask])
            vc = vl[mask] + vr[mask]
            var_tot = float(np.var(d))
            var_cnt = float(np.mean(vc))
            print(f"   {label:<26} n {int(mask.sum()):>6,}   Var(D) {var_tot:7.4f}   "
                  f"counting {var_cnt:7.4f}   EXCESS {max(var_tot - var_cnt, 0.0):7.4f}   "
                  f"p50|D| {np.median(np.abs(d)):6.3f}  p90|D| {np.percentile(np.abs(d), 90):6.3f}")
            return var_tot, var_cnt

        dstats(pair, rl, rr, "ALL pairs (summed cols)")
        dstats(deep, rl, rr, "deep pairs (summed cols)")
        dstats(deep, rlm, rrm, "deep pairs STRAND-MATCHED")
        clean = deep & (n_in == 0)
        dirty = deep & (n_in > 0)
        if clean.sum() >= 20:
            dstats(clean, rl, rr, "deep, NO terminus inside")
        if dirty.sum() >= 20:
            dstats(dirty, rl, rr, "deep, terminus(es) INSIDE")

        # the side-split center by exon-length decile — the deposit-rule / opportunity signature
        live = pair & (rho_t > 0)
        if live.sum() >= 50:
            q = np.quantile(length[live], np.linspace(0, 1, 6))
            print(f"   {'length bin':<18}{'n':>6}{'med cL':>9}{'med cR':>9}{'med D':>9}"
                  f"{'MAD(D)':>9}")
            for i in range(5):
                m = live & (length >= q[i]) & (length <= q[i + 1])
                if m.sum() < 10:
                    continue
                cl = np.median(np.log(rl[m] / rho_t[m]))
                cr = np.median(np.log(rr[m] / rho_t[m]))
                d = np.log(rl[m] / rr[m])
                print(f"   [{q[i]:6.0f},{q[i+1]:7.0f}] {int(m.sum()):>6,}{cl:>9.3f}{cr:>9.3f}"
                      f"{np.median(d):>9.3f}{np.median(np.abs(d - np.median(d))):>9.3f}")

        # the COMMON-MODE split: total, length-curve, truth-counting, and the left-over —
        # the half the pair estimator is structurally blind to
        if deep.sum() >= 50:
            m = deep & (rho_t > 0)
            c = np.log(np.sqrt(rl[m] * rr[m]) / rho_t[m])
            v_truth = float(np.mean(polygamma(1, n_m[m] + 0.5)))
            v_cnt_c = float(np.mean(vl[m] + vr[m]) / 4.0)
            L10 = length[m]
            qq = np.quantile(L10, np.linspace(0, 1, 11))
            fit10 = np.zeros_like(c)
            for i in range(10):
                mm = (L10 >= qq[i]) & (L10 <= qq[i + 1])
                if mm.sum() >= 5:
                    fit10[mm] = np.median(c[mm])
            resid = float(np.var(c - fit10))
            print(f"   COMMON-MODE Var(c) {np.var(c):7.4f} = length-curve {np.var(fit10):7.4f}"
                  f" + truth-counting {v_truth:7.4f} + flank-counting {v_cnt_c:7.4f}"
                  f" + LEFT OVER {max(resid - v_truth - v_cnt_c, 0.0):7.4f}")
            hi = m.copy()
            hi[m] = n_m[m] >= 100
            if hi.sum() >= 50:
                ch = np.log(np.sqrt(rl[hi] * rr[hi]) / rho_t[hi])
                vt = float(np.mean(polygamma(1, n_m[hi] + 0.5)))
                Lh = length[hi]
                qh = np.quantile(Lh, np.linspace(0, 1, 6))
                fh = np.zeros_like(ch)
                for i in range(5):
                    mm = (Lh >= qh[i]) & (Lh <= qh[i + 1])
                    if mm.sum() >= 5:
                        fh[mm] = np.median(ch[mm])
                print(f"   n_mrna>=100 subset (n {int(hi.sum()):,}): within-bin common Var "
                      f"{np.var(ch - fh):7.4f}, truth-counting {vt:7.4f}, LEFT OVER "
                      f"{max(float(np.var(ch - fh)) - vt, 0.0):7.4f}")

        # the single-flank population's center (the recorded residual candidate)
        sl = (hl > 0) & (hr == 0) & (rl > 0) & (rho_t > 0)
        sr = (hr > 0) & (hl == 0) & (rr > 0) & (rho_t > 0)
        if sl.sum() + sr.sum() >= 20:
            c = np.concatenate([np.log(rl[sl] / rho_t[sl]), np.log(rr[sr] / rho_t[sr])])
            print(f"   single-flank exons: n {int(sl.sum() + sr.sum()):,}  med center "
                  f"{np.median(c):+.3f}  (pairs' pooled med "
                  f"{np.median(np.log(np.sqrt(rl[live] * rr[live]) / rho_t[live])):+.3f})")
        sys.stdout.flush()
    return 0


if __name__ == "__main__":
    sys.exit(main())
