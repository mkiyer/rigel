"""DISSECTION 5 — ZERO gDNA: WHERE DOES THE FALSE gDNA COME FROM?

Owner: *"If there's zero DNA there's no [gDNA]. There's no seeding at the ends and the boundaries. Where is
this DNA coming from? ... The truth is we don't have any DNA in the intergenic space or the intronic space.
Clean as can be, so it's gonna be strange when we find what the false positive's coming from."*

`gdna_node_dissect` on the worst condition showed every top node with the same signature — `spliced = 0`,
`tau0_lam = 0`, neighbours with zero mass, and the strand channel ALONE landing at f_g ~ 0.51. That is the
uninformative reference, not a leak. But its replay reported a max fidelity error of 2.2e-01, so nothing
from it can be trusted until the mismatch is localized.

This does three things, in order:
  A  per-node REPLAY FIDELITY, so the attribution below is known-good rather than assumed;
  B  the FP mass censused by WHAT INFORMATION THE NODE HAD — the decisive split between "no evidence, so
     the reference answered 1/2" (a designed weakness the population prior exists to fix) and "evidence
     existed and was not used, or was used and was wrong" (a pass-0 defect);
  C  for the nodes that DID receive an RNA imputation message, what the message CLAIMED against what the
     truth is — because on a zero-gDNA library every node's true RNA fraction is exactly 1.

Run: OMP_NUM_THREADS=1 python scratchpad/gdna_d5_zero_source.py [cond]
"""
from __future__ import annotations

import dataclasses
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from flagship_interrogate import _oracle_per_node  # noqa: E402
from selfsolve_diag import _scan_and_truth  # noqa: E402

from rigel.calibration.calibrate import calibrate  # noqa: E402
from rigel.calibration.node_geometry import _node_region_type  # noqa: E402
from rigel.calibration.region_arrays import RegionArrays  # noqa: E402
from rigel.calibration.simplex_logodds import _solve_nodes_logodds_all  # noqa: E402
from rigel.config import PipelineConfig  # noqa: E402
from rigel.index import TranscriptIndex  # noqa: E402

SUITE = Path("/Users/mkiyer/Downloads/rigel_runs/ambig_dense_10mb")
_EPS = 1e-9
COND = sys.argv[1] if len(sys.argv) > 1 else "gdna_none_ss_0.50_nrna_none_capture_off"
CLS = {0: "intergenic", 1: "intron", 2: "exon", 3: "boundary"}


def load():
    index = TranscriptIndex.load(str(SUITE / "rigel_index"))
    cfg = PipelineConfig()
    ra = RegionArrays.from_region_df(index.region_df, index.ref_name_to_id)
    inp = _scan_and_truth(SUITE, COND, index, cfg, Path("/tmp/rigel_selfsolve"), SUITE / "_selfsolve_cache")
    dbg: dict = {}
    cc = dataclasses.replace(cfg.calibration, calib_refit_iters=0)
    calibrate(inp["payload"], ra, inp["strand_model"], np.asarray(inp["gdna_fl_pmf"]),
              np.asarray(inp["rna_fl_pmf"]), cc, _debug=dbg)
    return index, ra, inp, dbg, cc


def main():
    index, ra, inp, dbg, cc = load()
    chain, cap, st = dbg["chain"], dbg["capture"], dbg["statics"]
    Gp, Gn, Rp, Rn = _oracle_per_node(inp, chain)
    G = Gp + Gn
    fg = np.asarray(cap["f_g"])
    mass = np.asarray(cap["mass_global"])
    rt, _ = _node_region_type(chain, ra)
    cls = np.where(np.asarray(chain.kind) != 0, 3, rt)
    live = mass > _EPS
    fp = fg * mass                                   # on a zero-gDNA library this IS the error, exactly

    print(f"# {COND}\n# total unspliced mass {mass[live].sum():,.0f}   "
          f"true gDNA {G.sum():,.0f}   FALSE-POSITIVE gDNA {fp[live].sum():,.0f} "
          f"({100 * fp[live].sum() / mass[live].sum():.1f} % of mass)\n")

    # ── A. replay fidelity, per node ───────────────────────────────────────────────────────────────
    prs = dbg["calibration_priors"]
    kw = dict(kappa=float(dbg["rna_sense_frac"]),
              od_g=float(prs.gdna_strand_overdispersion), od_r=float(prs.rna_strand_overdispersion),
              n_grid=cc.sweep_n_grid, L=cc.sweep_logodds_window, n_tilt=cc.sweep_n_tilt,
              n_grid_ss=cc.sweep_n_grid_single_strand,
              fg_ref=cap["fg_init"], fpos_ref=cap["fpos_init"], fneg_ref=cap["fneg_init"])
    pos = (st.u_pos, st.u_neg, st.free_pos, st.free_neg, st.mass_unspliced, st.mass_spliced)

    def solve(**extra):
        return np.asarray(_solve_nodes_logodds_all(*pos, **{**kw, **extra}).gdna_frac)

    msgs = dict(global_logprior=cap["global_lp"], lam_logprior=cap["intron_prior"],
                gdna_imp_mode=cap["mode_g"], gdna_imp_prec=cap["prec_g"],
                rna_imp_mode=(cap["mode_p"], cap["mode_n"]),
                rna_imp_prec=(cap["prec_p"], cap["prec_n"]))
    f_all, f_strand = solve(**msgs), solve()
    err = np.abs(f_all - fg)
    print("=== A. REPLAY FIDELITY |replay(ALL) - shipped| ===")
    print(f"    live nodes {int(live.sum())}   median {np.median(err[live]):.2e}   "
          f"p95 {np.percentile(err[live], 95):.2e}   max {err[live].max():.2e}")
    good = err <= 1e-6
    print(f"    faithful (<=1e-6): {int((good & live).sum())}/{int(live.sum())} nodes, "
          f"carrying {100 * mass[good & live].sum() / mass[live].sum():.1f} % of the mass and "
          f"{100 * fp[good & live].sum() / max(fp[live].sum(), 1):.1f} % of the FP\n")

    # ── B. the FP mass by WHAT THE NODE HAD TO GO ON ───────────────────────────────────────────────
    has_spl = np.asarray(st.mass_spliced) > _EPS
    pg = np.asarray(cap["prec_g"], float)
    pp, pn = np.asarray(cap["prec_p"], float), np.asarray(cap["prec_n"], float)
    has_g, has_r = pg > 0, (pp > 0) | (pn > 0)
    fp_pos = np.asarray(st.free_pos, bool)
    fp_neg = np.asarray(st.free_neg, bool)
    strat = [
        ("NO evidence at all (no spliced, no msg)", live & ~has_spl & ~has_g & ~has_r),
        ("RNA message only", live & ~has_spl & ~has_g & has_r),
        ("gDNA message only", live & ~has_spl & has_g & ~has_r),
        ("both messages", live & ~has_spl & has_g & has_r),
        ("has OWN spliced mass", live & has_spl),
    ]
    print("=== B. WHERE THE FALSE gDNA IS, BY THE EVIDENCE THE NODE ACTUALLY HAD ===")
    print(f"{'stratum':42s} {'nodes':>6s} {'mass':>12s} {'FP mass':>12s} {'FP %':>6s} "
          f"{'share of FP':>12s} {'f_g strand-only':>16s}")
    for name, m in strat:
        if not m.any():
            continue
        print(f"{name:42s} {int(m.sum()):6d} {mass[m].sum():12,.0f} {fp[m].sum():12,.0f} "
              f"{100 * fp[m].sum() / max(mass[m].sum(), 1):6.1f} "
              f"{100 * fp[m].sum() / max(fp[live].sum(), 1):11.1f} % "
              f"{np.average(f_strand[m], weights=mass[m]):16.3f}")
    print("\n    by node class:")
    for c in (0, 1, 2, 3):
        m = live & (cls == c)
        if not m.any():
            continue
        amb = m & fp_pos & fp_neg
        print(f"    {CLS[c]:12s} {int(m.sum()):5d} nodes  mass {mass[m].sum():11,.0f}  "
              f"FP {fp[m].sum():11,.0f} ({100 * fp[m].sum() / max(mass[m].sum(), 1):5.1f} %)  "
              f"AMBIG share of FP {100 * fp[amb].sum() / max(fp[m].sum(), 1):5.1f} %")

    # ── C. what the RNA message CLAIMED, when there was one ────────────────────────────────────────
    print("\n=== C. THE RNA IMPUTATION MESSAGE — its claim vs the truth (true RNA fraction is 1.0) ===")
    print("    `rna_imp_mode` is a mode on log(f_active); the DEAD strand carries prec 0 and a sentinel")
    print("    mode, so only the live strand's may be read. exp(mode) is the claimed RNA fraction.\n")
    mp, mn = np.asarray(cap["mode_p"], float), np.asarray(cap["mode_n"], float)
    live_mode = np.where(pp > 0, mp, np.where(pn > 0, mn, np.nan))
    live_prec = np.where(pp > 0, pp, np.where(pn > 0, pn, np.nan))
    claim = np.exp(np.clip(live_mode, -50, 0.0))
    print(f"{'stratum':24s} {'n':>6s} {'mass':>12s} {'claimed f_rna':>14s} {'implied f_g':>12s} "
          f"{'actual f_g':>11s} {'med prec':>9s}")
    for name, m in (("exon, RNA msg", live & (cls == 2) & has_r),
                    ("boundary, RNA msg", live & (cls == 3) & has_r)):
        if not m.any():
            continue
        w = mass[m]
        c = float(np.average(claim[m], weights=w))
        print(f"{name:24s} {int(m.sum()):6d} {mass[m].sum():12,.0f} {c:14.3f} {1 - c:12.3f} "
              f"{np.average(fg[m], weights=w):11.3f} {np.nanmedian(live_prec[m]):9.2f}")
    print("\n    quintiles of the message's own precision (exons with an RNA message):")
    m = live & (cls == 2) & has_r
    q = np.nanpercentile(live_prec[m], [0, 20, 40, 60, 80, 100])
    print(f"    {'precision bin':22s} {'n':>5s} {'mass':>12s} {'claimed f_rna':>14s} {'f_g':>8s} "
          f"{'FP %':>7s}")
    for i in range(5):
        b = m & (live_prec >= q[i]) & (live_prec <= q[i + 1] if i == 4 else live_prec < q[i + 1])
        if not b.any():
            continue
        w = mass[b]
        print(f"    [{q[i]:8.2f},{q[i + 1]:8.2f}) {int(b.sum()):5d} {mass[b].sum():12,.0f} "
              f"{np.average(claim[b], weights=w):14.3f} {np.average(fg[b], weights=w):8.3f} "
              f"{100 * fp[b].sum() / max(mass[b].sum(), 1):7.1f}")
    print("\n    On a zero-gDNA library the truth is f_rna = 1.000 at every node. A message claiming less")
    print("    than that leaves the remainder to gDNA BY CONSTRUCTION — so if the claim is short, the")
    print("    defect is the message's MODE, upstream of every prior, and it is a pass-0 fix.")


if __name__ == "__main__":
    main()
