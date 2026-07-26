"""audit_3c — the FALSIFIABLE minimal test: does a 1000x reframe r leave the message PRECISION too high?

Replicate the _transport precision arithmetic in isolation and sweep r over 1e-3 .. 1e3.
If r leaks into precision, tpp/cm_p change with r. If not, they are invariant (the mode moves, precision does not).
"""
import numpy as np
from rigel.calibration.enrichment_frame import transfer_logvar, composition_logvar

EPS = 1e-9


def dv(psrc, s2t):
    return 1.0 / (1.0 / max(psrc, EPS) + s2t) if psrc > 0.0 else 0.0


# a source boundary: spliced-pos count SP=17 fragments, own measurement stream mp_src=9.
SP = 17.0
mp_src = 9.0
# dst is a high-mass exon (graft): s2t=0 by the graft exemption; also test the peel/plain s2t.
logvar_dst = 1.589e-5   # node 1909
logvar_src = 0.05168    # boundary 1910

print("edge = GRAFT (boundary->exon): transfer_logvar returns 0 regardless of r")
s2t_graft = float(transfer_logvar(np.array(logvar_dst), np.array(logvar_src), np.array(True)))
print(f"  s2t_graft = {s2t_graft}")
print(f"{'r':>10} {'mode=rho*r':>14} {'graft_spc':>12} {'prop=_dv(mp)':>14} {'tmp(=cm_p)':>12}")
rho_src = 0.03
for r in (1e-3, 1.0, 1e2, 1e3, 1e6):
    mode = rho_src * r
    spc = SP / (1.0 + SP * s2t_graft)
    prop = dv(mp_src, s2t_graft)
    tmp = spc + prop
    print(f"{r:>10.0e} {mode:>14.4g} {spc:>12.4f} {prop:>14.4f} {tmp:>12.4f}")

print("\nedge = PLAIN/PEEL: s2t = logvar_dst+logvar_src (Var(log r)), ALSO independent of r's magnitude")
s2t_plain = float(transfer_logvar(np.array(logvar_dst), np.array(logvar_src), np.array(False)))
print(f"  s2t_plain = {s2t_plain:.4g}  (= Var(log r), depends on COUNTS/composition, NOT on r itself)")
for r in (1e-3, 1.0, 1e2, 1e3, 1e6):
    mode = rho_src * r
    prop = dv(mp_src, s2t_plain)
    print(f"  r={r:.0e}  mode={mode:.4g}  _dv(mp_src)={prop:.4f}")

print("\nKEY POINT: composition_logvar (=Var(log rho_tot), the ingredient of s2t=Var(log r)) is a function of")
print("counts n and eff-lengths, NOT of r. For a well-counted node it is ~1e-5 even when r=1000x:")
for n in (10, 100, 1000, 66544):
    lv = float(composition_logvar(np.array(0.5), np.array(300.0), np.array(300.0), np.array(0.0), np.array(n)))
    print(f"  n={n:>6}: Var(log rho_tot) = 1/n contribution = {lv:.3g}  ->  a 1000x r transported at precision ~n")
