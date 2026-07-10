"""Flagship (gdna300, capture-on, ss0.99, zero-nascent) per-edge MESSAGE dissection for the logodds
under-call. For the worst AMBIG-exon region nodes, recompute each flanking message from the captured
local belief + geometry, splitting it into MODE (the imputed density/fraction — is it BIASED?) and the
THREE precision terms Var_own / σ²_bio / pois (is the precision OVERINFLATED?). Tests the user's two
hypotheses directly.

    OMP_NUM_THREADS=1 python scripts/debug/flagship_message_dissect.py
"""
from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
import importlib
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dissect_regions import build_or_load_cache  # noqa: E402

from rigel.calibration import calibrate  # noqa: E402
from rigel.calibration.node_chain import REGION, build_node_chain  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.signature import (  # noqa: E402
    coarse_strand_from_signature,
    coarse_type_from_signature,
)
from rigel.calibration.substrate import CalibrationSubstrate  # noqa: E402
from rigel.config import CalibrationConfig  # noqa: E402

_EPS = 1e-9
COND = "gdna_gdna300_ss_0.99_nrna_none_capture_on"


def main():
    index, blob = build_or_load_cache(COND, False)
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    payload = blob["payload_full"]
    cm = importlib.import_module("rigel.calibration.calibrate")
    orig = cm.node_sweep
    cap = {}
    cm.node_sweep = lambda *a, **k: orig(*a, _capture=cap, **k)
    try:
        calibrate(payload=payload, region_arrays=ra, strand_model=blob["strand_full"],
                  gdna_fl_pmf=blob["gdna_pmf"], rna_fl_pmf=blob["rna_pmf"],
                  config=CalibrationConfig(calibration_solver="logodds"))
    finally:
        cm.node_sweep = orig

    ch = build_node_chain(payload.ref_region_offsets, payload.ref_boundary_offsets)
    left, right = np.asarray(ch.left), np.asarray(ch.right)
    kind = np.asarray(ch.kind)
    refidx = np.asarray(ch.ref_idx, np.int64)
    is_reg = kind == REGION
    R = len(index.region_df)
    reg_node = {int(refidx[i]): int(i) for i in np.where(is_reg)[0]}

    sub_g = CalibrationSubstrate.from_payload(blob["payload_gdna"], ra)
    sub_r = CalibrationSubstrate.from_payload(blob["payload_rna"], ra)
    g_or = np.asarray(sub_g.contained.mass_unspliced, float)
    r_or = np.asarray(sub_r.contained.mass_unspliced, float)
    tot = g_or + r_or
    fg_true = np.where(tot > _EPS, g_or / np.maximum(tot, _EPS), np.nan)

    sig = np.asarray(ra.signature)
    scl = np.array([coarse_strand_from_signature(int(s)) for s in sig])
    tcl = np.array([coarse_type_from_signature(int(s)) for s in sig])

    fg_loc = np.asarray(cap["fg_loc"], float)
    vg_loc = np.asarray(cap["vg_loc"], float)
    f_fin = np.asarray(cap["f_g"], float)
    EGl, EGr = np.asarray(cap["eff_gdna_l"], float), np.asarray(cap["eff_gdna_r"], float)
    MSl, MSr = np.asarray(cap["mass_l"], float), np.asarray(cap["mass_r"], float)
    gvm = cap["gdna_vm_obj"] if "gdna_vm_obj" in cap else cap["gdna_vm"]

    # per-node, per-face gDNA density (the message currency) = fg_loc · M_face / E_face
    dGl = fg_loc * MSl / np.maximum(EGl, _EPS)
    dGr = fg_loc * MSr / np.maximum(EGr, _EPS)

    def msg_from(src, dst, src_face, dst_face):
        """One directed gDNA message src→dst. src_face/dst_face: 0=left,1=right. Returns dict."""
        d_src = dGr[src] if src_face == 1 else dGl[src]          # source facing density ρ_g_src
        eg_src = (EGr if src_face == 1 else EGl)[src]
        sm = (MSr if src_face == 1 else MSl)[src]
        md = (MSl if dst_face == 0 else MSr)[dst]
        egd = (EGl if dst_face == 0 else EGr)[dst]
        d_dst = (dGl if dst_face == 0 else dGr)[dst]
        n_src = fg_loc[src] * sm                                  # gDNA count at source
        mode = np.log(max(d_src * egd / max(md, _EPS), 1.0 / max(md, _EPS)))
        s2 = float(gvm.predict(np.array([max(0.5 * (d_src + d_dst), _EPS)]))[0])
        pois = 1.0 / max(n_src, _EPS)
        vown = vg_loc[src]
        prec = 1.0 / max(vown + s2 + pois, _EPS)
        return dict(src=src, n_src=n_src, d_src=d_src, mode=mode, impl=float(np.exp(min(mode, 0.0))),
                    vown=vown, s2=s2, pois=pois, prec=prec, fg_src=fg_loc[src])

    amb = np.array([r for r in range(R) if (scl[r] == 3) & (tcl[r] == 2) and tot[r] > 200])
    nodes = np.array([reg_node[int(r)] for r in amb])
    err = np.abs(f_fin[nodes] - fg_true[amb]) * tot[amb]
    worst = np.argsort(-err)[:5]
    print(f"================= {COND} : AMBIG-exon message dissection =================")
    print("Each AMBIG-exon region: its local f_g, final f_g, ORACLE, and the two flanking-boundary gDNA")
    print("messages (mode→impl_f, + precision split Var_own / σ²_bio / pois → prec).\n")
    for j in worst:
        r = int(amb[j]); n = nodes[j]
        lb, rb = int(left[n]), int(right[n])
        print(f"REGION {r} (node {n}) mass={tot[r]:,.0f}  ORACLE f_g={fg_true[r]:.3f}  "
              f"local={fg_loc[n]:.3f}  final={f_fin[n]:.3f}")
        for side, b, sf, df in (("LEFT-bnd", lb, 1, 0), ("RIGHT-bnd", rb, 0, 1)):
            if b < 0:
                print(f"    {side}: (ref end)")
                continue
            m = msg_from(b, n, sf, df)
            print(f"    {side} node {b}: src f_g={m['fg_src']:.3f} src_dens={m['d_src']:.4f} "
                  f"n_src(gDNAct)={m['n_src']:.1f}")
            print(f"        message mode={m['mode']:.2f} (impl f_g={m['impl']:.4f})  "
                  f"prec={m['prec']:.3f} = 1/(Var_own {m['vown']:.3f} + σ²_bio {m['s2']:.3f} + pois {m['pois']:.3f})")
            print(f"        ⇒ MODE biased? impl_f {m['impl']:.4f} vs oracle {fg_true[r]:.3f}  |  "
                  f"dominant prec term: {max((m['vown'],'Var_own'),(m['s2'],'σ²_bio'),(m['pois'],'pois'))[1]}")
        print()


if __name__ == "__main__":
    main()
