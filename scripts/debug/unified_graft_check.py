"""GRAFT MAGNITUDE CHECK — is the boundary's measured spliced density the exon's mature density?

The unified solver's `boundary → exon` GRAFT (`unified_solver_design.md` §6) adds the junction's measured
mature density ``ρ_μ = spliced_mass / E_spl`` to the RNA the exon receives, and the combine then re-forms a
mass as ``ρ_μ · E_r(exon)``. That is only right if ``E_spl`` and ``E_r`` are commensurate frames — and there
is a standing warning in `bp_solver.py` that they are not (*"a density in the SPLICED eff-length frame
while rho_pos is in the CONTAINED RNA frame, and that spliced frame is known-wrong on a face where the exon
continues"*).

This probe settles it against the oracle: for every EXON node, the grafted mature MASS predicted from its
flanking boundaries vs the oracle's own contained **unspliced mature** (``mat_uns_pos + mat_uns_neg``).
A ratio far from 1 names the frame error and its size.

    RIGEL_UNIFIED=1 OMP_NUM_THREADS=1 python scripts/debug/unified_graft_check.py [cond]
"""

from __future__ import annotations

import sys

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scripts/debug")
from unified_message_audit import run  # noqa: E402

from rigel.calibration.bp_solver import REGION  # noqa: E402

_EPS = 1e-12


def _oracle_mature_nascent(inp, chain):
    """Per chain-node oracle contained-unspliced MATURE and NASCENT mass (both strands)."""
    kind = np.asarray(chain.kind)
    idx = np.asarray(chain.ref_idx, np.int64)
    isr = kind == REGION

    def pick(p, keys):
        return sum(np.asarray(p[k], float) for k in keys)

    mat_k = ("mat_uns_pos", "mat_uns_neg")
    nas_k = ("nas_uns_pos", "nas_uns_neg")
    out = []
    for keys in (mat_k, nas_k):
        a, b = pick(inp["region_pools"], keys), pick(inp["boundary_pools"], keys)
        ri = np.clip(idx, 0, a.shape[0] - 1)
        bi = np.clip(idx, 0, b.shape[0] - 1)
        out.append(np.where(isr, a[ri], b[bi]))
    return out  # mature, nascent


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    cond = args[0] if args else "gdna_gdna300_ss_0.50_nrna_present_capture_on"
    d = run(cond)
    cap, cls, chain, inp = d["cap"], d["cls"], d["chain"], d["inp"]
    st = cap["_uni_static"]
    M, E_r = st["M"], st["E_r"]
    spl = st["spl_p"] + st["spl_n"]  # node-pooled mature DENSITY at a boundary (0 on regions)
    left, right = st["left"], st["right"]
    mat, nas = _oracle_mature_nascent(inp, chain)

    n = M.size
    sl, sr = np.clip(left, 0, n - 1), np.clip(right, 0, n - 1)
    # the mature density an exon would receive: whichever flanking boundary carries the junction
    rho_mu = np.where(left >= 0, spl[sl], 0.0) + np.where(right >= 0, spl[sr], 0.0)
    pred_mature = rho_mu * E_r  # the mass the combine re-forms in the exon's contained frame

    ex = (cls == 2) & (M > 1e-9) & (mat + nas > 1e-9)
    print(f"# cond={cond}   exon nodes with RNA: {int(ex.sum())}")
    print(f"# oracle mature fraction of exon RNA: "
          f"{np.average((mat / np.maximum(mat + nas, _EPS))[ex], weights=M[ex]):.3f}")

    r = pred_mature[ex] / np.maximum(mat[ex], _EPS)
    fin = np.isfinite(r) & (mat[ex] > 1e-9) & (pred_mature[ex] > 0)
    rr = r[fin]
    print(f"\n=== grafted mature MASS / oracle mature MASS, at exons (n={rr.size}) ===")
    for q in (10, 25, 50, 75, 90):
        print(f"  p{q:<3} {np.percentile(rr, q):>10.3f}")
    print(f"  geomean {np.exp(np.mean(np.log(np.maximum(rr, 1e-12)))):>8.3f}"
          f"   mass-wtd mean {np.average(rr, weights=M[ex][fin]):>8.3f}")

    # the same comparison as a FRACTION of the exon's total mass — what the message actually asserts
    fm_pred = (pred_mature / np.maximum(M, _EPS))[ex]
    fm_true = (mat / np.maximum(M, _EPS))[ex]
    print("\n=== as a fraction of the exon's unspliced mass (mass-weighted) ===")
    print(f"  grafted f_mature = {np.average(fm_pred, weights=M[ex]):.3f}"
          f"    oracle f_mature = {np.average(fm_true, weights=M[ex]):.3f}")

    # the frame ratio itself
    bnd = cls == 3
    print(f"\n# median E_spl(boundary) = {np.median(cap['_uni_static'].get('E_r', E_r)[bnd]):.1f} (E_r proxy)")
    print(f"# median E_r(exon) = {np.median(E_r[ex]):.1f}   median M(exon) = {np.median(M[ex]):.1f}")


if __name__ == "__main__":
    main()
