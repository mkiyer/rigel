"""§7.1 counterfactual: current σ²_bio(μ) message precision vs the proposed disagreement-aware σ²_edge,
computed per gDNA edge on the CONVERGED belief. NO solver change, NO golden regen.

For each directed gDNA edge src→dst (the message _scan would build):
  mo        = log( ρ_g(src, facing dst) re-expressed in dst's log-f_g frame )   [message mode]
  base_var  = vg_loc[src] + pois,   pois = 1/gdna_count_src                      [sampling + source belief]
  CURRENT : pr = 1/( base_var + σ²_bio(μ) ),   μ = ½(ρ_src+ρ_dst)                [the monotone curve]
  PROPOSED: resid    = mo − lfg_loc[dst]            (vs dst's MESSAGE-FREE local anchor — non-circular)
            expected = base_var + vg_loc[dst]       (variance sampling+belief already explain)
            σ²_edge  = max(resid² − expected, 0)
            pr = 1/( base_var + σ²_edge )

Criteria:
  (i)   confident (single-strand) dst + large |resid| (a depleted-seam message) → trust ≪ current (silenced)
  (ii)  any dst + small |resid| (within-regime relay)                          → trust ≈ current (relay survives)
  (iii) uncertain (AMBIG) dst + large |resid| (a correct enriched lift)         → trust ≈ current (lift survives)

Usage:  python precision_counterfactual.py [condition_subdir]
"""
import sys
from pathlib import Path

import numpy as np

from rigel.calibration import bp_solver as bp
from rigel.calibration.fl import build_fl_models, gdna_fl_mass
from rigel.calibration.node_chain import build_node_chain
from rigel.calibration.region_arrays import RegionArrays
from rigel.calibration.strand_balance import fit_strand_balance
from rigel.calibration.substrate import BoundarySubstrate, CalibrationSubstrate
from rigel.config import BamScanConfig, CalibrationConfig
from rigel.index import TranscriptIndex
from rigel.pipeline import scan_and_buffer
from rigel.splice import SpliceType

_EPS = 1.0e-9
ROOT = Path.home() / "Downloads/rigel_runs/quick_3to1_5mb"
COND = sys.argv[1] if len(sys.argv) > 1 else "gdna_gdna300_ss_0.99_nrna_none_capture_on"
BAM = ROOT / COND / "sim_oracle.bam"

idx = TranscriptIndex.load(str(ROOT / "rigel_index"))
st_, sm, flm, buf, pl = scan_and_buffer(str(BAM), idx, BamScanConfig())
fl = build_fl_models(global_counts=flm.global_model.counts,
                     rna_counts=flm.category_models[SpliceType.SPLICED_ANNOT].counts,
                     gdna_counts=gdna_fl_mass(pl), max_size=flm.max_size)
ra = RegionArrays.from_region_df(idx.region_df, idx.ref_name_to_id)
sub = CalibrationSubstrate.from_payload(pl, ra)
bsub = BoundarySubstrate.from_payload(pl)
chain = build_node_chain(pl.ref_region_offsets, pl.ref_boundary_offsets)
geom = bp.build_node_geometry(chain, sub, bsub, ra, fl.gdna_pmf, fl.rna_pmf)
statics = bp.build_node_statics(chain, sub, bsub, ra)
kappa = float(fit_strand_balance(sm).rna_sense_frac)
belief0 = bp.init_beliefs(chain, sub, bsub, ra, rna_sense_frac=kappa, n_grid=60, statics=statics)
cap = {}
belief, _ = bp.node_sweep(chain, statics, geom, belief0, ra, bsub, rna_sense_frac=kappa,
                          n_grid=60, max_outer=CalibrationConfig().sweep_max_outer, _capture=cap)

# converged belief (message content/relay) + message-free local anchor (phase A of the last outer iter)
fg = np.asarray(belief.f_g)
fg_loc = np.asarray(cap["fg_loc"]); vg_loc = np.asarray(cap["vg_loc"])
lfg_loc = np.log(np.maximum(fg_loc, _EPS))
gdna_vm = cap["gdna_vm"]
fp = np.asarray(cap["free_pos"], bool); fn = np.asarray(cap["free_neg"], bool)
ambig = fp & fn                          # both strands live  → strand carries ~0 gDNA info (uncertain dst)
single = fp ^ fn
EG = (np.asarray(geom.eff_gdna_left), np.asarray(geom.eff_gdna_right))
MS = (np.asarray(geom.mass_left), np.asarray(geom.mass_right))
left = np.asarray(chain.left); right = np.asarray(chain.right)
node_rtype, _ = bp._node_region_type(chain, ra)  # 0 intergenic,1 intron,2 exon,-1 boundary
node_rtype = np.asarray(node_rtype)
# local densities (for the σ²_bio μ query, matching _scan) + converged densities (message content)
dGloc = (fg_loc * MS[0] / np.maximum(EG[0], _EPS), fg_loc * MS[1] / np.maximum(EG[1], _EPS))
dGcon = (fg * MS[0] / np.maximum(EG[0], _EPS), fg * MS[1] / np.maximum(EG[1], _EPS))

rows = []  # (resid, |resid|, pr_cur, pr_new, ratio, s2bio, s2edge, dst_ambig, dst_single, src_rtype)
for nbr, sf, df in ((left, 1, 0), (right, 0, 1)):   # forward (src=left), backward (src=right)
    for i in np.where(nbr >= 0)[0]:
        s = int(nbr[i])
        sm_ = MS[sf][s]
        if sm_ <= _EPS:
            continue
        eg = max(EG[sf][s], _EPS); egd = max(EG[df][i], _EPS); md = max(MS[df][i], _EPS)
        n_src = fg[s] * sm_                        # source gDNA count (converged)
        if n_src <= _EPS:
            continue
        rho = n_src / eg
        mo = np.log(max(rho, 1.0 / egd) / (md / egd))
        pois = 1.0 / max(n_src, _EPS)
        base_var = vg_loc[s] + pois
        mu = max(0.5 * (dGloc[sf][s] + dGloc[df][i]), _EPS)
        s2bio = max(float(gdna_vm.predict(np.array([mu]))[0]), 0.0)
        pr_cur = 1.0 / max(base_var + s2bio, _EPS)
        resid = mo - lfg_loc[i]
        expected = base_var + vg_loc[i]
        s2edge = max(resid * resid - expected, 0.0)
        pr_new = 1.0 / max(base_var + s2edge, _EPS)
        rows.append((resid, abs(resid), pr_cur, pr_new, pr_new / max(pr_cur, _EPS),
                     s2bio, s2edge, bool(ambig[i]), bool(single[i]), int(node_rtype[s])))

R = np.array([r[:7] for r in rows])
A = np.array([r[7] for r in rows], bool); S = np.array([r[8] for r in rows], bool)
SRT = np.array([r[9] for r in rows])  # source region coarse type (0 intergenic,1 intron,2 exon,-1 bnd)
resid, ares, pr_cur, pr_new, ratio, s2bio, s2edge = (R[:, k] for k in range(7))
drag = resid < 0.0   # message says LESS gDNA than dst's own anchor (the capture-leak direction)
lift = resid > 0.0   # message says MORE gDNA (a possible enriched lift)
src_dep = (SRT == 0) | (SRT == 1)  # source is depleted-by-structure (intergenic/intron)
print(f"\n===== {COND}  (kappa={kappa:.3f}) =====")
print(f"gDNA directed edges: {R.shape[0]:,}   (AMBIG-dst={int(A.sum()):,}  single-strand-dst={int(S.sum()):,})")
print(f"current σ²_bio range: [{s2bio.min():.3f}, {s2bio.max():.3f}]   median |resid|={np.median(ares):.3f}")
vgl_amb = float(np.median(vg_loc[ambig])) if ambig.any() else float("nan")
vgl_ss = float(np.median(vg_loc[single])) if single.any() else float("nan")
print(f"dst belief log-var vg_loc median:  AMBIG={vgl_amb:.2f}   single-strand={vgl_ss:.2f}   "
      f"(expected absorbs resid² up to ~vg_loc)")


def summary(label, mask):
    if mask.sum() == 0:
        print(f"  {label:<40} n=0"); return
    rr = ratio[mask]
    print(f"  {label:<40} n={int(mask.sum()):>6}  med|resid|={np.median(ares[mask]):5.2f}  "
          f"med ratio={np.median(rr):7.3f}  silenced(<0.1)={np.mean(rr < 0.1):4.2f}  kept(>0.5)={np.mean(rr > 0.5):4.2f}")


print("\n  dst-class × DIRECTION × magnitude  →  trust ratio (proposed / current):")
for cls, m in (("single-strand (CONFIDENT) dst", S), ("AMBIG (UNCERTAIN) dst", A)):
    print(f" [{cls}]")
    summary("   |resid|<0.5  (agree / relay)", m & (ares < 0.5))
    summary("   DRAG  resid≤−1.5 (depleted-seam pull-DOWN)", m & drag & (ares >= 1.5))
    summary("     └ of which src=depleted (intron/intergenic)", m & drag & (ares >= 1.5) & src_dep)
    summary("   LIFT  resid≥+1.5 (push-UP)", m & lift & (ares >= 1.5))
    summary("     └ of which src=depleted (cross-seam, silence OK)", m & lift & (ares >= 1.5) & src_dep)
    summary("     └ of which src=exon (within-regime lift, keep)", m & lift & (ares >= 1.5) & (SRT == 2))

relay = ares < 0.5
print("\n  Headline:")
print(f"  relay (|resid|<0.5)            median trust ×{np.median(ratio[relay]):.2f}  (kept>0.5: {np.mean(ratio[relay] > 0.5):.2f})")
print(f"  DRAG seam (resid≤−1.5)         median trust ×{np.median(ratio[drag & (ares >= 1.5)]):.3f}  (silenced<0.1: {np.mean(ratio[drag & (ares >= 1.5)] < 0.1):.2f})")
print(f"  → relay/drag trust discrimination = {np.median(ratio[relay]) / max(np.median(ratio[drag & (ares >= 1.5)]), 1e-9):.0f}×")
