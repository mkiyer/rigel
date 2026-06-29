"""DENSITY-space message dissection (no fractions). For the worst under-called exon nodes, compares the
node's TRUE gDNA density ρ_g (oracle counts / eff-len) to the gDNA densities its flanking messages
IMPUTE (ρ_g_src) — testing whether the gDNA message asserts a FALSE density continuity across the
capture enrichment seam, and whether σ²_bio (the between-node density log-spread) is honest about it.

    OMP_NUM_THREADS=1 python scripts/debug/density_message_dissect.py [cond]
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
from rigel.calibration.effective_length import region_eff_length  # noqa: E402
from rigel.calibration.node_chain import REGION, build_node_chain  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.signature import (  # noqa: E402
    coarse_strand_from_signature,
    coarse_type_from_signature,
)
from rigel.calibration.substrate import CalibrationSubstrate  # noqa: E402
from rigel.config import CalibrationConfig  # noqa: E402

_EPS = 1e-9
DEFAULT = "gdna_gdna300_ss_0.50_nrna_none_capture_on"
SCN = {1: "POS", 2: "NEG", 3: "AMBIG"}


def main():
    cond = sys.argv[1] if len(sys.argv) > 1 else DEFAULT
    index, blob = build_or_load_cache(cond, False)
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
    reg_eff_g = region_eff_length(np.asarray(ra.region_size_bp, float), blob["gdna_pmf"])
    rho_g_true = np.where(g_or > _EPS, g_or / np.maximum(reg_eff_g, _EPS), 0.0)  # TRUE gDNA density
    fg_true = np.where(tot > _EPS, g_or / np.maximum(tot, _EPS), np.nan)

    sig = np.asarray(ra.signature)
    scl = np.array([coarse_strand_from_signature(int(s)) for s in sig])
    tcl = np.array([coarse_type_from_signature(int(s)) for s in sig])

    fg_loc = np.asarray(cap["fg_loc"], float)
    f_fin = np.asarray(cap["f_g"], float)
    vg_loc = np.asarray(cap["vg_loc"], float)
    EGl, EGr = np.asarray(cap["eff_gdna_l"], float), np.asarray(cap["eff_gdna_r"], float)
    MSl, MSr = np.asarray(cap["mass_l"], float), np.asarray(cap["mass_r"], float)
    gvm = cap["gdna_vm"]
    dGl = fg_loc * MSl / np.maximum(EGl, _EPS)   # solved gDNA density, left face
    dGr = fg_loc * MSr / np.maximum(EGr, _EPS)

    def src_density(src, face):  # the message's imputed gDNA DENSITY ρ_g_src (this is the message content)
        return (dGr if face == 1 else dGl)[src]

    print(f"\n========= {cond} : DENSITY-space message dissection =========")
    print("Per exon node: TRUE ρ_g vs the gDNA densities its flanking boundaries IMPUTE (the message says")
    print("dst ρ_g ≈ src ρ_g). σ²_bio = between-node density log-spread; ratio = true ρ_g / imputed ρ_g.\n")
    exo = np.array([r for r in range(R) if (tcl[r] == 2) and tot[r] > 1000])
    nodes = np.array([reg_node[int(r)] for r in exo])
    err = np.abs(f_fin[nodes] - fg_true[exo]) * tot[exo]
    for j in np.argsort(-err)[:7]:
        r = int(exo[j]); n = nodes[j]
        print(f"R{r} [{SCN.get(int(scl[r]),'?')}] mass={tot[r]:,.0f}  TRUE ρ_g={rho_g_true[r]:.3f}  "
              f"solved ρ_g(loc)={dGl[n]:.3f}  (oracle f_g={fg_true[r]:.3f} → final {f_fin[n]:.3f})")
        for side, b, sf in (("L", int(left[n]), 1), ("R", int(right[n]), 0)):
            if b < 0:
                continue
            ds = src_density(b, sf)
            mean = 0.5 * (ds + dGl[n])
            s2 = float(gvm.predict(np.array([max(mean, _EPS)]))[0])
            ratio = rho_g_true[r] / max(ds, _EPS)
            print(f"    {side}-bnd {b}: imputes ρ_g={ds:.3f}  vs TRUE {rho_g_true[r]:.3f}  "
                  f"(true/imputed = {ratio:.1f}×)  σ²_bio={s2:.2f}  Var_own={vg_loc[b]:.2f}")
        print()


if __name__ == "__main__":
    main()
