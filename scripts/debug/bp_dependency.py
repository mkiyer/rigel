"""The interdependency (pie-chart) problem + the Trojan-horse failure, and the correct fix.

Components sum to the observed total density D (a simplex / pie chart), so they are DEPENDENT: you cannot move
one without the others moving. Degrees of freedom: single-strand (2 active comps) = 1 DoF; AMBIG (3) = 2 DoF.
D is observed, NOT a DoF.

TROJAN HORSE (the failure): if you update components INDEPENDENTLY and re-impose the sum, a message aimed at a
WEAK component moves it, and the sum then FORCES a CONFIDENT component to change — even though nothing informed
the confident one. A confident gDNA can be flipped by a message about RNA.

FIX: never update components independently. Solve JOINTLY on the degrees of freedom, with EVERY component's
information (its local belief precision AND its incoming message) entered as a simultaneous constraint. The
confident component then pins its share; a message to a weak component is absorbed by the OTHER weak components,
never by the confident one.

Run:  python scripts/debug/bp_dependency.py
"""
from __future__ import annotations
import math
import numpy as np

EPS = 1e-9
LOG = lambda x: math.log(max(x, EPS))  # noqa: E731


def joint_solve(D, comps, n=600):
    """JOINT solve on the simplex: components ρ_c ≥ 0 with Σρ_c = D, minimizing Σ_c p_c·(log ρ_c − m_c)².
    comps = [(name, m_c, p_c)]; p_c may come from the local belief AND/OR a message (already combined per
    component as independent log-density Gaussians — that part is fine; the DEPENDENCY is handled here by the
    single constrained solve). 2 comps ⇒ 1-D grid; 3 comps ⇒ 2-D simplex grid."""
    names = [c[0] for c in comps]; m = np.array([c[1] for c in comps]); p = np.array([c[2] for c in comps])
    k = len(comps)
    if k == 2:
        f0 = np.linspace(EPS, 1 - EPS, n)
        rho = np.stack([f0 * D, (1 - f0) * D])
        cost = p[0] * (np.log(rho[0]) - m[0]) ** 2 + p[1] * (np.log(rho[1]) - m[1]) ** 2
        i = int(np.argmin(cost)); return names, rho[:, i]
    # k == 3: grid the 2-simplex (f0, f1), f2 = 1−f0−f1 > 0
    g = np.linspace(EPS, 1 - EPS, n)
    F0, F1 = np.meshgrid(g, g); F2 = 1 - F0 - F1
    ok = F2 > EPS
    r0, r1, r2 = F0 * D, F1 * D, F2 * D
    cost = np.where(ok,
                    p[0] * (np.log(np.where(ok, r0, 1)) - m[0]) ** 2
                    + p[1] * (np.log(np.where(ok, r1, 1)) - m[1]) ** 2
                    + p[2] * (np.log(np.where(ok, r2, 1)) - m[2]) ** 2, np.inf)
    j = np.unravel_index(int(np.argmin(cost)), cost.shape)
    return names, np.array([r0[j], r1[j], r2[j]])


def naive_update(D, conf_name, conf_rho, weak_name, weak_rho, weak_prec, msg_m, msg_prec):
    """The WRONG (Trojan-horse) update for a single-strand node: combine the WEAK component with its message
    independently, then make the CONFIDENT component the residual (D − weak) — ignoring the confident
    component's own precision. This is what 'update the addressed component + re-impose the sum' does."""
    lw = (weak_prec * LOG(weak_rho) + msg_prec * msg_m) / (weak_prec + msg_prec)
    weak_new = math.exp(lw)
    conf_new = D - weak_new  # residual — the confident component is dragged
    return conf_new, weak_new


def show(names, rho, note=""):
    D = rho.sum()
    print("      " + "  ".join(f"{n}={r:6.2f}" for n, r in zip(names, rho)) + f"   (Σ={D:.1f}) {note}")


print("=" * 92)
print("SINGLE-STRAND (2 comps, 1 DoF): confident gDNA=30 (prec 50), weak RNA+=3 (prec 1). D=33.")
print("A message arrives about RNA+ saying ρ_RNA+ ≈ 10 (prec 3).\n")
print("  NAIVE (update weak RNA+ + residual gDNA) — the TROJAN HORSE:")
cg, cr = naive_update(33.0, "gDNA", 30.0, "RNA+", 3.0, 1.0, LOG(10.0), 3.0)
show(["gDNA", "RNA+"], np.array([cg, cr]), "← confident gDNA got DRAGGED from 30")
print("\n  JOINT solve on the DoF (gDNA local prec 50 + RNA+ local 1 combined with the msg prec 3):")
# RNA+ local(3,prec1) ⊗ msg(10,prec3): combined target
rp_m = (1 * LOG(3.0) + 3 * LOG(10.0)) / (1 + 3); rp_p = 1 + 3
names, rho = joint_solve(33.0, [("gDNA", LOG(30.0), 50.0), ("RNA+", rp_m, rp_p)])
show(names, rho, "← confident gDNA PROTECTED; the weak RNA+ barely moves")
print("  ⇒ the confident component pins its share; a weak message cannot flip it. The DoF is solved ONCE with")
print("    all precisions, not by updating one component and residual-ing the other.")

print("\n" + "=" * 92)
print("AMBIG (3 comps, 2 DoF): confident gDNA=25 (prec 50), weak RNA+=4 (prec 1), weak RNA-=4 (prec 1). D=33.")
print("A message arrives about RNA+ saying ρ_RNA+ ≈ 12 (prec 4).  Which component absorbs it?\n")
names, rho0 = joint_solve(33.0, [("gDNA", LOG(25.0), 50.0), ("RNA+", LOG(4.0), 1.0), ("RNA-", LOG(4.0), 1.0)])
show(names, rho0, "before the message")
rp_m = (1 * LOG(4.0) + 4 * LOG(12.0)) / (1 + 4); rp_p = 1 + 4
names, rho1 = joint_solve(33.0, [("gDNA", LOG(25.0), 50.0), ("RNA+", rp_m, rp_p), ("RNA-", LOG(4.0), 1.0)])
show(names, rho1, "after the RNA+ message")
print(f"  ⇒ gDNA stays ~{rho1[0]:.0f} (confident, PROTECTED); RNA+ rises ({rho0[1]:.1f}→{rho1[1]:.1f}); the")
print(f"    OTHER weak component RNA- absorbs the trade ({rho0[2]:.1f}→{rho1[2]:.1f}). The change flows to the")
print("    least-confident components, never the confident one. This is the correct simplex behaviour.")
