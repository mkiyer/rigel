"""P1e — the shared substrate: run one condition, and lay the solver's PRE-PIN message state flat.

Every row is ONE MESSAGE (one direction into one destination node) at the FINAL rho-iteration, with:

    delta   = log(M / S)                       the conservation VIOLATION  (S = the claim's asserted mass)
    aSa     = alpha^T Sigma alpha              the violation the claim DECLARED it could produce, with
                                               Sigma = sigma_cm^2 * 11^T + diag(w)   (conservation_rescale's model)
    z2      = delta^2 / (aSa + 1/n_dst)        the SURPRISE
    fg_msg  = the composition the message delivers   (= sigmoid of the lambda it emits)
    fo      = the destination node's ORACLE f_g

Nothing here is a model change; it re-reads `_capture["_pin"]` / `_capture["_uni_static"]`.
"""

from __future__ import annotations

import dataclasses
import importlib
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from flagship_interrogate import _oracle_per_node  # noqa: E402
from selfsolve_diag import _scan_and_truth  # noqa: E402

from rigel.calibration.bp_solver import REGION  # noqa: E402
from rigel.calibration.node_geometry import _node_region_type  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

calmod = importlib.import_module("rigel.calibration.calibrate")
SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
_EPS = 1e-9

CONDS = {
    "gdna300_ss0.99_present_capOFF": "gdna_gdna300_ss_0.99_nrna_present_capture_off",
    "gdna300_ss0.50_present_capOFF": "gdna_gdna300_ss_0.50_nrna_present_capture_off",
    "gdna300_ss0.99_present_capON": "gdna_gdna300_ss_0.99_nrna_present_capture_on",
    "gdna100_ss0.50_present_VERYSTRONG": "gdna_gdna100_ss_0.50_nrna_present_capture_verystrong",
    "none_ss0.50_present_capOFF": "gdna_none_ss_0.50_nrna_present_capture_off",
    # the P1d residual-regression stratum (task item 5)
    "gdna100_ss0.50_none_capOFF": "gdna_gdna100_ss_0.50_nrna_none_capture_off",
    "gdna300_ss0.50_none_capOFF": "gdna_gdna300_ss_0.50_nrna_none_capture_off",
}

_index = None
_ra = None
_cfg = None


def _boot():
    global _index, _ra, _cfg
    if _index is None:
        _index = TranscriptIndex.load(str(SUITE / "rigel_index"))
        _cfg = PipelineConfig()
        _ra = RegionArrays.from_region_df(_index.region_df, _index.ref_name_to_id)
    return _index, _ra, _cfg


def solve(cond: str, refit: int = 0):
    """Run calibrate on one condition and return (inp, dbg)."""
    index, ra, cfg = _boot()
    inp = _scan_and_truth(
        SUITE, cond, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache"
    )
    dbg: dict = {}
    calmod.calibrate(
        inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
        np.asarray(inp["rna_fl_pmf"]),
        dataclasses.replace(cfg.calibration, calib_refit_iters=refit), _debug=dbg,
    )
    return inp, dbg


def node_frame(inp, dbg):
    """Per-NODE arrays: class, mass, oracle f_g, solved f_g, own densities."""
    index, ra, _ = _boot()
    chain, cap = dbg["chain"], dbg["capture"]
    us = cap["_uni_static"]
    kind = np.asarray(chain.kind)
    n = kind.shape[0]
    rt, _ = _node_region_type(chain, ra)
    CLS = {0: "intergenic", 1: "intron", 2: "exon"}
    cls = np.array(
        [CLS.get(int(rt[j]), "?") if kind[j] == REGION else "boundary" for j in range(n)]
    )
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
    T = Gp + Gn + Rp + Rn
    fo = np.where(T > _EPS, (Gp + Gn) / np.maximum(T, _EPS), np.nan)
    return {
        "cls": cls,
        "kind": kind,
        "M": np.asarray(us["M"], float),
        "E_g": np.asarray(us["E_g"], float),
        "E_r": np.asarray(us["E_r"], float),
        "og": np.asarray(us["og"], float),
        "op": np.asarray(us["op"], float),
        "on": np.asarray(us["on"], float),
        "fo": fo,
        "fg": np.asarray(cap["f_g"], float),
        "var_g": np.asarray(cap["var_g"], float),
        "tau_own": np.asarray(us["tau_own"], float),
        "n_unspl_l": np.asarray(us["n_unspl_l"], float),
        "n_unspl_r": np.asarray(us["n_unspl_r"], float),
        "is_amb": np.asarray(us["fp"], bool) & np.asarray(us["fn"], bool),
        "mass_oracle": T,
    }


def message_table(inp, dbg, last_iter_only=True):
    """One row per (direction, destination node). Returns a dict of flat arrays."""
    cap = dbg["capture"]
    nf = node_frame(inp, dbg)
    pins = cap["_pin"]
    # `_pin` is appended twice per rho-iteration (left msg then right msg). Keep the FINAL iteration.
    sel = pins[-2:] if last_iter_only else pins
    M, E_g, E_r = nf["M"], nf["E_g"], nf["E_r"]
    og, op, on = nf["og"], nf["op"], nf["on"]
    is_reg = nf["kind"] == REGION
    n_dst = np.where(is_reg, nf["n_unspl_l"], nf["n_unspl_l"] + nf["n_unspl_r"])

    out = {k: [] for k in (
        "df", "dst", "src", "cls", "graft", "delta", "aSa", "z2", "bhat2",
        "alpha_g", "alpha_p", "alpha_n", "s_g", "s_p", "s_n", "w_g", "w_r", "s2cm",
        "fg_msg", "fo", "M", "S", "lam_emit", "tau_own", "is_amb", "cls_src",
        "dv_g", "dv_r", "n_dst", "n_src", "v_g", "v_r", "prec_g", "prec_r", "spl_prec",
        "nsup", "sup_g", "sup_r",
    )}
    for p in sel:
        valid = np.asarray(p["valid"], bool)
        src = np.asarray(p["src"], np.int64)
        tg, tp, tn = p["tg"], p["tp"], p["tn"]
        pg, pp, pn = p["tpg"], p["tpp"], p["tpn"]
        s2cm = np.asarray(p["s2t"], float) + 1.0 / np.maximum(np.asarray(p["n_src"], float), _EPS)
        supplied = np.stack([pg > 0.0, pp > 0.0, pn > 0.0], axis=-1)
        rho = np.stack([tg, tp, tn], axis=-1)
        own = np.stack([og, op, on], axis=-1)
        eff = np.stack([E_g, E_r, E_r], axis=-1)
        m = np.where(supplied, rho, own) * eff
        S = m.sum(axis=-1)
        ok = valid & (S > _EPS) & (M > _EPS)
        alpha = m / np.maximum(S, _EPS)[..., None]
        v = np.where(supplied, 1.0 / np.maximum(np.stack([pg, pp, pn], axis=-1), _EPS), np.inf)
        w = np.maximum(np.where(supplied, v, 0.0) - s2cm[..., None], 0.0)
        s = np.where(supplied, s2cm[..., None] + alpha * w, 0.0)
        aSa = np.sum(alpha * s, axis=-1)
        delta = np.where(ok, np.log(np.maximum(M, _EPS) / np.maximum(S, _EPS)), np.nan)
        denom = aSa + 1.0 / np.maximum(n_dst, _EPS)
        z2 = delta * delta / np.maximum(denom, _EPS)
        bhat2 = np.maximum(delta * delta - denom, 0.0)
        # the DERIVED per-component attribution: Sigma += bhat2 * s s^T / (aSa)^2  ⇒  dv_c = bhat2*(s_c/aSa)^2
        dv = bhat2[..., None] * (s / np.maximum(aSa, _EPS)[..., None]) ** 2
        tR = tp + tn
        lam_emit = (tg > _EPS) & (tR > _EPS)
        fg_msg = np.where(
            (tg * E_g + tR * E_r) > _EPS,
            tg * E_g / np.maximum(tg * E_g + tR * E_r, _EPS),
            np.nan,
        )
        k = np.flatnonzero(ok & np.isfinite(nf["fo"]))
        out["df"].append(np.full(k.shape, p["df"]))
        out["dst"].append(k)
        out["src"].append(src[k])
        out["cls"].append(nf["cls"][k])
        out["cls_src"].append(nf["cls"][src[k]])
        out["graft"].append(np.asarray(p["graft"], bool)[k])
        out["delta"].append(delta[k])
        out["aSa"].append(aSa[k])
        out["z2"].append(z2[k])
        out["bhat2"].append(bhat2[k])
        out["alpha_g"].append(alpha[k, 0])
        out["alpha_p"].append(alpha[k, 1])
        out["alpha_n"].append(alpha[k, 2])
        out["s_g"].append(s[k, 0])
        out["s_p"].append(s[k, 1])
        out["s_n"].append(s[k, 2])
        out["w_g"].append(w[k, 0])
        out["w_r"].append(w[k, 1] + w[k, 2])
        out["s2cm"].append(s2cm[k])
        out["fg_msg"].append(fg_msg[k])
        out["fo"].append(nf["fo"][k])
        out["M"].append(M[k])
        out["S"].append(S[k])
        out["lam_emit"].append(lam_emit[k])
        out["tau_own"].append(nf["tau_own"][k])
        out["is_amb"].append(nf["is_amb"][k])
        out["dv_g"].append(dv[k, 0])
        out["dv_r"].append(dv[k, 1] + dv[k, 2])
        out["n_dst"].append(n_dst[k])
        out["n_src"].append(np.asarray(p["n_src"], float)[k])
        out["v_g"].append(v[k, 0])
        out["v_r"].append(np.minimum(v[k, 1], v[k, 2]))
        out["prec_g"].append(pg[k])
        out["prec_r"].append((pp + pn)[k])
        out["spl_prec"].append(np.asarray(p["spl_prec"], float)[k])
        out["nsup"].append(supplied[k].sum(axis=-1))
        out["sup_g"].append(supplied[k, 0])
        out["sup_r"].append(supplied[k, 1] | supplied[k, 2])
    return {kk: np.concatenate(vv) for kk, vv in out.items()}, nf


def q(x, w=None):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return dict.fromkeys(("n", "med", "p25", "p75", "p90", "p99"), float("nan")) | {"n": 0}
    return {
        "n": int(x.size),
        "med": float(np.median(x)),
        "p25": float(np.percentile(x, 25)),
        "p75": float(np.percentile(x, 75)),
        "p90": float(np.percentile(x, 90)),
        "p99": float(np.percentile(x, 99)),
    }


def spearman(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 5:
        return float("nan")
    ra = np.argsort(np.argsort(a[m])).astype(float)
    rb = np.argsort(np.argsort(b[m])).astype(float)
    ra -= ra.mean()
    rb -= rb.mean()
    d = np.sqrt((ra * ra).sum() * (rb * rb).sum())
    return float((ra * rb).sum() / d) if d > 0 else float("nan")
