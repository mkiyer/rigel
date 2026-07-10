"""Node-by-node FB dissection on the flagship: where is FB's error injected — refit/local, or propagation?

For each chain REGION node compares three f_g values:
  - true  : the oracle gDNA fraction (payload_gdna vs payload_rna contained — the split BAMs).
  - local : the message-FREE belief (FB step A: strand + global prior + Jeffreys). Wrong here ⇒ refit/global fault.
  - final : after the FB messages (α⊗β). Wrong here while local was right ⇒ PROPAGATION fault.

Decomposes the gDNA-mass error into local_err = (local−true)·M and prop_err = (final−local)·M, aggregates by
region type × strand class, ranks the worst nodes, and for each dumps the forward α and backward β gDNA messages
(mode, prec) + the neighbours, so we can see WHICH message dragged it and WHY (a cliff = a low-gDNA neighbour
relayed into a high-gDNA exon).

    OMP_NUM_THREADS=1 python scripts/debug/fb_node_dissect.py [cond] [--top N]
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
from rigel.calibration.node_chain import BOUNDARY, REGION, build_node_chain  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.signature import (  # noqa: E402
    coarse_strand_from_signature,
    coarse_type_from_signature,
)
from rigel.calibration.substrate import CalibrationSubstrate  # noqa: E402
from rigel.config import CalibrationConfig  # noqa: E402

DEFAULT = "gdna_gdna300_ss_0.50_nrna_none_capture_on"
SC = {0: "NONE", 1: "POS", 2: "NEG", 3: "AMBIG"}
TC = {0: "intgc", 1: "intron", 2: "exon"}
_EPS = 1e-9


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    cond = args[0] if args else DEFAULT
    top = int(sys.argv[sys.argv.index("--top") + 1]) if "--top" in sys.argv else 15
    index, blob = build_or_load_cache(cond, False)
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    payload = blob["payload_full"]

    # capture the FB internals via a node_sweep wrapper
    calmod = importlib.import_module("rigel.calibration.calibrate")
    orig = calmod.node_sweep
    cap = {}

    def wrap(*a, **k):
        return orig(*a, _capture=cap, **k)

    calmod.node_sweep = wrap
    try:
        calibrate(payload=payload, region_arrays=ra, strand_model=blob["strand_full"],
                  gdna_fl_pmf=blob["gdna_pmf"], rna_fl_pmf=blob["rna_pmf"], config=CalibrationConfig())
    finally:
        calmod.node_sweep = orig

    # rebuild the SAME chain calibrate used (deterministic from payload offsets)
    ch = build_node_chain(payload.ref_region_offsets, payload.ref_boundary_offsets)
    kind = np.asarray(ch.kind)
    refidx = np.asarray(ch.ref_idx, np.int64)

    is_reg = kind == REGION

    # oracle per region (split BAMs); f_g_true on the UNSPLICED mass (matches the solver's f_g)
    g_or = np.asarray(CalibrationSubstrate.from_payload(blob["payload_gdna"], ra).contained.mass_unspliced, float)
    r_or = np.asarray(CalibrationSubstrate.from_payload(blob["payload_rna"], ra).contained.mass_unspliced, float)
    M = g_or + r_or
    fg_true_reg = g_or / np.maximum(M, _EPS)

    sig = np.asarray(ra.signature)
    scl = np.array([coarse_strand_from_signature(int(s)) for s in sig])
    tcl = np.array([coarse_type_from_signature(int(s)) for s in sig])

    fg_loc, fg_fin, fg_str = cap["fg_loc"], cap["f_g"], cap["fg_strand"]
    amg, apg = cap["a_fwd"][0], cap["a_fwd"][1]   # forward gDNA message (mode, prec) into each node
    bmg, bpg = cap["b_bwd"][0], cap["b_bwd"][1]   # backward gDNA message
    solv = np.asarray(cap["solvable"], bool)

    # per REGION node row
    rows = []
    R = M.shape[0]
    for i in np.where(is_reg & solv)[0]:
        rgi = int(refidx[i])
        if rgi >= R or M[rgi] <= _EPS:
            continue
        t, st, lo, fi = fg_true_reg[rgi], float(fg_str[i]), float(fg_loc[i]), float(fg_fin[i])
        rows.append(dict(
            node=i, reg=rgi, sc=int(scl[rgi]), tc=int(tcl[rgi]), mass=float(M[rgi]),
            true=t, strand=st, local=lo, final=fi,
            strand_err=(st - t) * M[rgi], global_err=(lo - st) * M[rgi],
            local_err=(lo - t) * M[rgi], prop_err=(fi - lo) * M[rgi], tot_err=(fi - t) * M[rgi],
            amode=float(amg[i]), aprec=float(apg[i]), bmode=float(bmg[i]), bprec=float(bpg[i]),
        ))

    print(f"=== {cond} : FB node-by-node (REGION nodes; f_g on unspliced; gDNA-mass-weighted error) ===")
    def tot(key):
        return sum(abs(r[key]) for r in rows), sum(r[key] for r in rows)
    a_str, s_str = tot("strand_err")
    a_glo, s_glo = tot("global_err")
    a_prop, s_prop = tot("prop_err")
    a_tot, s_tot = tot("tot_err")
    print("\nERROR DECOMPOSITION (gDNA-mass-weighted, Σ over solvable region nodes; local = strand + global):")
    print(f"  strand likelihood   |err|={a_str:>11,.0f}   signed {s_str:>+12,.0f}")
    print(f"  global prior  ê(z)   |err|={a_glo:>11,.0f}   signed {s_glo:>+12,.0f}")
    print(f"  FB propagation       |err|={a_prop:>11,.0f}   signed {s_prop:>+12,.0f}")
    print(f"  ── TOTAL (final−true) |err|={a_tot:>11,.0f}   signed {s_tot:>+12,.0f}")
    denom = max(a_str + a_glo + a_prop, 1)
    print(f"  → shares of |err|:  strand {100*a_str/denom:.0f}%  global {100*a_glo/denom:.0f}%  "
          f"propagation {100*a_prop/denom:.0f}%")

    print("\nBY region-type × strand-class  (signed gDNA-mass error: local / prop):")
    print(f"  {'class':>6} {'type':>7} {'nreg':>5} {'Σlocal_err':>12} {'Σprop_err':>12} {'Σtot_err':>12}")
    for tc in (2, 1, 0):
        for s in (3, 1, 2, 0):
            sub = [r for r in rows if r["tc"] == tc and r["sc"] == s]
            if not sub:
                continue
            print(f"  {SC[s]:>6} {TC[tc]:>7} {len(sub):>5} {sum(r['local_err'] for r in sub):>+12,.0f} "
                  f"{sum(r['prop_err'] for r in sub):>+12,.0f} {sum(r['tot_err'] for r in sub):>+12,.0f}")

    print(f"\nTOP {top} OFFENDING nodes by |total err|  (true→strand→local→final; a=fwd α, b=bwd β gDNA msg):")
    print(f"  {'node':>5} {'cls':>5} {'type':>6} {'mass':>8} {'true':>5} {'strnd':>5} {'local':>5} {'final':>5} "
          f"{'aMode':>6} {'aPrec':>6} {'bMode':>6} {'bPrec':>6}  diagnosis")
    for r in sorted(rows, key=lambda r: -abs(r["tot_err"]))[:top]:
        tag = []
        if abs(r["strand"] - r["true"]) > 0.15:
            tag.append("STRAND" + ("↑" if r["strand"] > r["true"] else "↓"))
        if abs(r["local"] - r["strand"]) > 0.10:
            tag.append("GLOBAL" + ("↑" if r["local"] > r["strand"] else "↓"))
        if abs(r["final"] - r["local"]) > 0.10:
            who = "fwd" if r["aprec"] > r["bprec"] else "bwd"
            tag.append(f"PROP{'↑' if r['final'] > r['local'] else '↓'}({who})")
        if not tag:
            tag.append("·")
        print(f"  {r['node']:>5} {SC[r['sc']]:>5} {TC[r['tc']]:>6} {r['mass']:>8,.0f} {r['true']:>5.2f} "
              f"{r['strand']:>5.2f} {r['local']:>5.2f} {r['final']:>5.2f} {r['amode']:>6.2f} {r['aprec']:>6.1f} "
              f"{r['bmode']:>6.2f} {r['bprec']:>6.1f}  {' '.join(tag)}")

    # ── BOUNDARY EMISSION: which boundaries emit a gDNA message, which don't, why ──
    left, right = np.asarray(ch.left), np.asarray(ch.right)
    is_bnd = kind == BOUNDARY
    ml, mr = cap["mass_l"], cap["mass_r"]          # unspliced crossing mass per face
    sl, sr = cap["spl_l"], cap["spl_r"]            # spliced (mature) mass per face
    rt_node = np.where(is_reg, tcl[np.clip(refidx, 0, R - 1)], -1)  # region type per chain node (−1 on bnd)

    def flank(b):  # (left flank type, right flank type) of boundary b, via its chain neighbours
        lt = rt_node[left[b]] if left[b] >= 0 else -1
        rtt = rt_node[right[b]] if right[b] >= 0 else -1
        return int(lt), int(rtt)

    emits_g = (ml > _EPS) | (mr > _EPS)  # gDNA: facing unspliced crossing mass (strand-agnostic structural gate)
    has_spl = (sl > _EPS) | (sr > _EPS)
    print("\nBOUNDARY EMISSION (a boundary emits a gDNA msg iff solvable & has UNSPLICED crossing mass):")
    cats = {}
    for b in np.where(is_bnd)[0]:
        lt, rtt = flank(b)
        key = "-".join(sorted([TC.get(lt, "end"), TC.get(rtt, "end")]))
        d = cats.setdefault(key, dict(n=0, g=0, spl_only=0, neither=0, unspl_mass=0.0, spl_mass=0.0))
        d["n"] += 1
        d["unspl_mass"] += float(ml[b] + mr[b])
        d["spl_mass"] += float(sl[b] + sr[b])
        if emits_g[b]:
            d["g"] += 1
        elif has_spl[b]:
            d["spl_only"] += 1            # emits the mature/RNA floor but NO gDNA
        else:
            d["neither"] += 1
    print(f"  {'flank-pair':>16} {'nbnd':>5} {'emit_gDNA':>9} {'spliced-only':>12} {'neither':>8} "
          f"{'Σunspl':>10} {'Σspliced':>10}")
    for k, d in sorted(cats.items(), key=lambda kv: -kv[1]["n"]):
        print(f"  {k:>16} {d['n']:>5} {d['g']:>9} {d['spl_only']:>12} {d['neither']:>8} "
              f"{d['unspl_mass']:>10,.0f} {d['spl_mass']:>10,.0f}")

    # exons starved of any gDNA message (BOTH flanking boundaries emit no gDNA) — stuck on local
    starved_mass = starved_n = 0
    for rrow in rows:
        i = rrow["node"]
        lb, rb = left[i], right[i]
        lg = lb >= 0 and emits_g[lb]
        rg = rb >= 0 and emits_g[rb]
        if rrow["tc"] == 2 and not lg and not rg:  # an exon with no gDNA message either side
            starved_n += 1
            starved_mass += rrow["mass"]
    print(f"\n  EXON region nodes with NO gDNA message (both flanks gDNA-silent) = {starved_n} "
          f"holding {starved_mass:,.0f} unspliced mass → gDNA solved on LOCAL only (no cross-exon relay).")


if __name__ == "__main__":
    main()
