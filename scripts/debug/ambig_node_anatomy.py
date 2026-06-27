"""Granular anatomy of an AMBIG region node: transcript structure, exact region/boundary values, the strand
solve (the 'tilt'), the global gDNA density ê(z), the propagation messages, and where the error originates.

For the worst-mispredicted AMBIG exon nodes on a condition (or explicit region ids), prints:
  - the overlapping transcripts (the +/− overlap that makes it AMBIG),
  - the region's raw accumulator values (per-strand unspliced counts, mass, spliced, eff-lengths),
  - the two flanking boundaries (neighbour type/strand, per-side gDNA-crossing mass),
  - the SOLVE TRAJECTORY: strand-only f_g (the tilt) → +global ê(z) (local) → +messages (final) → vs ORACLE,
  - so the error's origin (strand / global / propagation) is visible with exact numbers.

    OMP_NUM_THREADS=1 python scripts/debug/ambig_node_anatomy.py [cond] [--regions r1,r2,...] [--top N]
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
from rigel.calibration.node_chain import BOUNDARY, REGION, build_node_chain  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.signature import (  # noqa: E402
    coarse_strand_from_signature,
    coarse_type_from_signature,
)
from rigel.calibration.substrate import CalibrationSubstrate  # noqa: E402
from rigel.config import CalibrationConfig  # noqa: E402

DEFAULT = "gdna_gdna300_ss_0.99_nrna_none_capture_on"
SC = {0: "NONE", 1: "POS", 2: "NEG", 3: "AMBIG"}
TC = {0: "intergenic", 1: "intron", 2: "exon"}
_EPS = 1e-9


def main():
    argv = [a for a in sys.argv[1:] if not a.startswith("--")]
    cond = argv[0] if argv else DEFAULT
    want = None
    if "--regions" in sys.argv:
        want = [int(x) for x in sys.argv[sys.argv.index("--regions") + 1].split(",")]
    top = int(sys.argv[sys.argv.index("--top") + 1]) if "--top" in sys.argv else 3

    index, blob = build_or_load_cache(cond, False)
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    payload = blob["payload_full"]

    calmod = importlib.import_module("rigel.calibration.calibrate")
    orig = calmod.node_sweep
    cap = {}
    calmod.node_sweep = lambda *a, **k: orig(*a, _capture=cap, **k)
    try:
        cal = calibrate(payload=payload, region_arrays=ra, strand_model=blob["strand_full"],
                        gdna_fl_pmf=blob["gdna_pmf"], rna_fl_pmf=blob["rna_pmf"], config=CalibrationConfig())
    finally:
        calmod.node_sweep = orig

    ch = build_node_chain(payload.ref_region_offsets, payload.ref_boundary_offsets)
    kind = np.asarray(ch.kind)
    refidx = np.asarray(ch.ref_idx, np.int64)
    left, right = np.asarray(ch.left), np.asarray(ch.right)
    is_reg = kind == REGION
    R = len(index.region_df)
    reg_node = {int(refidx[i]): int(i) for i in np.where(is_reg)[0]}  # region id → chain node

    # oracle composition (split BAMs)
    sub_g = CalibrationSubstrate.from_payload(blob["payload_gdna"], ra)
    sub_r = CalibrationSubstrate.from_payload(blob["payload_rna"], ra)
    sub_f = CalibrationSubstrate.from_payload(blob["payload_full"], ra)
    g_or = np.asarray(sub_g.contained.mass_unspliced, float)
    r_or = np.asarray(sub_r.contained.mass_unspliced, float)
    fg_true = g_or / np.maximum(g_or + r_or, _EPS)
    reg_eff_g = region_eff_length(np.asarray(ra.region_size_bp, float), blob["gdna_pmf"])

    upos = np.asarray(sub_f.contained.n_unspliced_pos, float)
    uneg = np.asarray(sub_f.contained.n_unspliced_neg, float)
    nspl = np.asarray(sub_f.contained.n_spliced_sense, float)
    massu = np.asarray(sub_f.contained.mass_unspliced, float)

    sig = np.asarray(ra.signature)
    scl = np.array([coarse_strand_from_signature(int(s)) for s in sig])
    tcl = np.array([coarse_type_from_signature(int(s)) for s in sig])
    start = index.region_df["start"].to_numpy()
    end = index.region_df["end"].to_numpy()
    refname = index.region_df["ref_name"].to_numpy()

    fg_loc, fg_str, fg_fin = cap["fg_loc"], cap["fg_strand"], cap["f_g"]
    amg, apg = cap["a_fwd"][0], cap["a_fwd"][1]
    bmg, bpg = cap["b_bwd"][0], cap["b_bwd"][1]
    pp, pn = cap["prec_p"], cap["prec_n"]
    ehat, z, eff_glob, mass_glob = cap["ehat"], cap["z_enrich"], cap["eff_global"], cap["mass_global"]
    rho_global = float(cap["rho_global"])
    kappa = float(cal.rna_sense_frac)
    massl, massr = cap["mass_l"], cap["mass_r"]
    rtype_node = np.where(is_reg, tcl[np.clip(refidx, 0, R - 1)], -1)
    sc_node = np.where(is_reg, scl[np.clip(refidx, 0, R - 1)], -1)

    # ── SS-CENSUS: same table for SINGLE-STRAND (POS|NEG) exon nodes (the candidate larger error source) ──
    if "--ss-census" in sys.argv:
        ssx = np.where(((scl == 1) | (scl == 2)) & (tcl == 2))[0]
        cols = ("rid node cls mass true strand local final z ehat_rho mu_g "
                "strand_err global_err prop_err tot_err").split()
        print("\t".join(cols))
        for rid in ssx:
            n = reg_node[rid]
            m = g_or[rid] + r_or[rid]
            if m <= _EPS:
                continue
            zi = float(z[n]) if np.isfinite(z[n]) else float("nan")
            rho_hat = (float(max(ehat.predict(np.array([zi]))[0], 0.0)) if np.isfinite(zi) else rho_global)
            mu = float(np.clip(rho_hat * eff_glob[n] / max(mass_glob[n], _EPS), 0.0, 1.0))
            t_, st, lo, fi = fg_true[rid], fg_str[n], fg_loc[n], fg_fin[n]
            print(f"{rid}\t{n}\t{SC[int(scl[rid])]}\t{m:.0f}\t{t_:.3f}\t{st:.3f}\t{lo:.3f}\t{fi:.3f}\t{zi:.3f}\t"
                  f"{rho_hat:.3f}\t{mu:.3f}\t{(st-t_)*m:.0f}\t{(lo-st)*m:.0f}\t{(fi-lo)*m:.0f}\t{(fi-t_)*m:.0f}")
        return

    # ── CENSUS: a TSV of EVERY AMBIG-exon node (the full population), one row per region ──
    if "--census" in sys.argv:
        amb = np.where((scl == 3) & (tcl == 2))[0]
        cols = ("rid node mass u_pos u_neg true strand local final z ehat_rho mu_g "
                "strand_err global_err prop_err tot_err").split()
        print("\t".join(cols))
        for rid in amb:
            n = reg_node[rid]
            m = g_or[rid] + r_or[rid]
            if m <= _EPS:
                continue
            zi = float(z[n]) if np.isfinite(z[n]) else float("nan")
            rho_hat = (float(max(ehat.predict(np.array([zi]))[0], 0.0)) if np.isfinite(zi) else rho_global)
            mu = float(np.clip(rho_hat * eff_glob[n] / max(mass_glob[n], _EPS), 0.0, 1.0))
            t_, st, lo, fi = fg_true[rid], fg_str[n], fg_loc[n], fg_fin[n]
            print(f"{rid}\t{n}\t{m:.0f}\t{upos[rid]:.0f}\t{uneg[rid]:.0f}\t{t_:.3f}\t{st:.3f}\t{lo:.3f}\t"
                  f"{fi:.3f}\t{zi:.3f}\t{rho_hat:.3f}\t{mu:.3f}\t{(st-t_)*m:.0f}\t{(lo-st)*m:.0f}\t"
                  f"{(fi-lo)*m:.0f}\t{(fi-t_)*m:.0f}")
        return

    # pick targets: explicit regions, else the worst AMBIG-exon nodes by |final−true|·mass
    tdf = index.t_df
    if want is None:
        amb = np.where((scl == 3) & (tcl == 2) & (g_or + r_or > 100))[0]
        err = np.abs(fg_fin[[reg_node[int(r)] for r in amb]] - fg_true[amb]) * (g_or + r_or)[amb]
        want = [int(amb[i]) for i in np.argsort(-err)[:top]]

    print(f"================= {cond} : AMBIG node anatomy =================")
    print(f"library: κ(rna_sense_frac)={kappa:.3f}  od_gdna={cal.gdna_strand_overdispersion:.3f}  "
          f"od_rna={cal.rna_strand_overdispersion:.3f}  ρ_global(gDNA density)={rho_global:.3f}")

    for rid in want:
        n = reg_node[rid]
        s, e, rf = int(start[rid]), int(end[rid]), str(refname[rid])
        print("\n" + "═" * 100)
        print(f"REGION {rid} (chain node {n})  {rf}:{s}-{e}  len={e-s}  sig={int(sig[rid])} "
              f"class={SC[int(scl[rid])]}/{TC[int(tcl[rid])]}")

        # ---- transcript structure (the +/− overlap that makes it AMBIG) ----
        ov = tdf[(tdf.ref == rf) & (tdf.start < e) & (tdf.end > s)]
        print(f"\n  TRANSCRIPT STRUCTURE — {len(ov)} transcript(s) overlap this region:")
        for _, t in ov.iterrows():
            strand = "+" if int(t.strand) == 1 else "-"
            kindt = "nRNA" if bool(t.is_nrna) else "mRNA"
            print(f"    {t.t_id:<14} {t.g_name:<8} strand={strand} {kindt}  span={int(t.start)}-{int(t.end)} "
                  f"n_exons={int(t.n_exons)}")
        strands = {("+" if int(x) == 1 else "-") for x in ov.strand}
        print(f"    ⇒ strands present: {sorted(strands)}  "
              f"{'(opposite-strand overlap ⇒ AMBIG: own strand undefined)' if strands == {'+', '-'} else ''}")

        # ---- region accumulator (the data the node sees) ----
        print("\n  REGION DATA (what the node measures):")
        print(f"    unspliced strand counts:  u_pos={upos[rid]:,.0f}  u_neg={uneg[rid]:,.0f}  "
              f"(ratio +:{upos[rid]/max(upos[rid]+uneg[rid],_EPS):.2f})   mass_unspliced={massu[rid]:,.0f}")
        print(f"    spliced (mature) sense counts={nspl[rid]:,.0f}")
        print(f"    eff_len gdna(contained)={reg_eff_g[rid]:,.0f}")

        # ---- flanking boundaries ----
        print("\n  FLANKING BOUNDARIES (gDNA crossing each seam):")
        for side, b in (("LEFT ", int(left[n])), ("RIGHT", int(right[n]))):
            if b < 0:
                print(f"    {side} boundary: (reference end)")
                continue
            # the boundary's OTHER neighbour (the region on the far side of the seam)
            far = int(left[b]) if int(right[b]) == n else int(right[b])
            far_t = TC.get(int(rtype_node[far]), "end") if far >= 0 else "end"
            far_s = SC.get(int(sc_node[far]), "") if far >= 0 else ""
            far_r = int(refidx[far]) if far >= 0 else -1
            print(f"    {side} boundary node {b}: far neighbour R{far_r} ({far_s}/{far_t})  "
                  f"gDNA-cross mass L={massl[b]:,.0f} R={massr[b]:,.0f}")

        # ---- THE STRAND SOLVE (the tilt) ----
        print("\n  ① STRAND SOLVE (the TILT — essential, but cannot solve f_g alone at AMBIG):")
        print(f"     strand-only f_g = {fg_str[n]:.3f}   "
              f"(at AMBIG, balanced ± counts read as 'could be gDNA OR balanced RNA' — the count-zero-info wall)")

        # ---- GLOBAL gDNA density ê(z) ----
        zi = float(z[n]) if np.isfinite(z[n]) else float("nan")
        applied = bool(cap["enrich_apply"][n])
        rho_hat = float(max(ehat.predict(np.array([zi]))[0], 0.0)) if np.isfinite(zi) else rho_global
        mu = float(np.clip(rho_hat * eff_glob[n] / max(mass_glob[n], _EPS), 0.0, 1.0))
        print("\n  ② GLOBAL gDNA DENSITY  ê(z)  (the population gDNA-density prior AMBIG nodes lean on):")
        print(f"     z(enrichment predictor)={zi:.3f}  ê(z)=ρ̂_g={rho_hat:.3f}  "
              f"{'[ê applied]' if applied else '[flat ρ_global=%.3f]' % rho_global}")
        print(f"     ⇒ implied prior μ_g = clip(ρ̂_g·E/M) = {mu:.3f}   "
              f"(oracle gDNA density ρ_g_true={g_or[rid]/max(reg_eff_g[rid],_EPS):.3f})")
        print(f"     local f_g (strand ⊗ global) = {fg_loc[n]:.3f}")

        # ---- PROPAGATION ----
        print("\n  ③ PROPAGATION (neighbour messages — the other essential AMBIG input):")
        print(f"     fwd α gDNA: mode={amg[n]:.3f} prec={apg[n]:.2f}   "
              f"bwd β gDNA: mode={bmg[n]:.3f} prec={bpg[n]:.2f}")
        print(f"     incoming RNA precision: +={pp[n]:.2f}  −={pn[n]:.2f}")

        # ---- trajectory + error origin ----
        t_, lo, fi = fg_true[rid], fg_loc[n], fg_fin[n]
        print("\n  SOLVE TRAJECTORY → ORACLE:")
        print(f"     strand-only {fg_str[n]:.3f}  →  +global (local) {lo:.3f}  →  +messages (final) {fi:.3f}"
              f"   ‖  ORACLE true f_g = {t_:.3f}")
        de_str = (fg_str[n] - t_)
        de_glo = (lo - fg_str[n])
        de_prop = (fi - lo)
        m = g_or[rid] + r_or[rid]
        print(f"     error decomposition (×mass {m:,.0f}):  strand {de_str:+.3f}  global {de_glo:+.3f}  "
              f"propagation {de_prop:+.3f}   TOTAL {fi-t_:+.3f}  ({(fi-t_)*m:+,.0f} frags)")
        worst = max((abs(de_str), "strand"), (abs(de_glo), "global"), (abs(de_prop), "propagation"))[1]
        print(f"     ⇒ dominant error source: {worst.upper()}")


if __name__ == "__main__":
    main()
