"""Verify the (count, 1/L) estimators by direct simulation. No rigel imports -- pure first principles.

Ground truth: molecules of length L ~ f(L) are deposited on a template with START density rho
(molecules per template position). We ask what three observables estimate:

  (a) COUNT of molecules COVERING a point p
  (b) SUM 1/L over molecules COVERING p
  (c) SUM 1/L over molecules SPANNING a seam (>=1 base on each side)
  (d) COUNT of molecules CONTAINED in a region of length ell
"""

import numpy as np

rng = np.random.default_rng(7)

T = 2_000_000  # template length
RHO = 0.02  # molecules per position
MU, SD = 200.0, 60.0  # fragment length distribution


def draw(n):
    L = rng.normal(MU, SD, n).round().astype(np.int64)
    return np.clip(L, 30, 600)


n_mol = rng.poisson(RHO * T)
starts = rng.integers(0, T, n_mol)
lens = draw(n_mol)
ends = starts + lens  # half-open [start, end)

fl = np.bincount(lens, minlength=601)[:601].astype(float)
fl /= fl.sum()
Lv = np.arange(601)
EL = float((fl * Lv).sum())
EinvL = float((fl[1:] / Lv[1:]).sum())
print(f"truth: rho={RHO}  E[L]={EL:.3f}  E[1/L]={EinvL:.6f}  E[L]*E[1/L]={EL * EinvL:.4f}")

# probe points well inside the template
probes = rng.integers(1000, T - 1000, 4000)
order = np.argsort(starts)
s_s, e_s, l_s = starts[order], ends[order], lens[order]

cov_count = np.empty(probes.size)
cov_dens = np.empty(probes.size)
span_dens = np.empty(probes.size)
lo_all = np.searchsorted(s_s, probes - 600, side="left")
hi_all = np.searchsorted(s_s, probes, side="right")
for i, p in enumerate(probes):
    a, b = lo_all[i], hi_all[i]
    s, e, L = s_s[a:b], e_s[a:b], l_s[a:b]
    cov = e > p  # covers position p  (start <= p < end)
    cov_count[i] = cov.sum()
    cov_dens[i] = (1.0 / L[cov]).sum()
    span = (s <= p - 1) & (e > p)  # >=1 base each side of the seam between p-1 and p
    span_dens[i] = (1.0 / L[span]).sum()

print()
print("(a) COUNT covering a point")
print(f"    mean {cov_count.mean():.4f}   predicted rho*E[L] = {RHO * EL:.4f}")
print("(b) SUM 1/L covering a point")
print(
    f"    mean {cov_dens.mean():.6f}   predicted rho       = {RHO:.6f}   "
    f"ratio {cov_dens.mean() / RHO:.5f}"
)
print(f"    var  {cov_dens.var():.3e}    predicted rho*E[1/L]= {RHO * EinvL:.3e}")
print("(c) SUM 1/L SPANNING a seam")
print(
    f"    mean {span_dens.mean():.6f}   predicted rho*(1-E[1/L]) = {RHO * (1 - EinvL):.6f}   "
    f"ratio {span_dens.mean() / RHO:.5f}  <- {100 * EinvL:.3f}% LOW"
)
print("    the fix: sum 1/(L-1) instead -> predicted exactly rho")

span_dens_fix = np.empty(probes.size)
for i, p in enumerate(probes):
    a, b = lo_all[i], hi_all[i]
    s, e, L = s_s[a:b], e_s[a:b], l_s[a:b]
    span = (s <= p - 1) & (e > p)
    span_dens_fix[i] = (1.0 / (L[span] - 1.0)).sum()
print(f"    sum 1/(L-1) mean {span_dens_fix.mean():.6f}  ratio {span_dens_fix.mean() / RHO:.5f}")

print()
print("(d) CONTAINED in a region of length ell -- does 1/L rescue it?")
for ell in (100, 300, 1000, 5000):
    lo = rng.integers(1000, T - ell - 1000, 3000)
    cnt = np.empty(lo.size)
    den = np.empty(lo.size)
    a_ = np.searchsorted(s_s, lo, side="left")
    b_ = np.searchsorted(s_s, lo + ell, side="right")
    for i in range(lo.size):
        s, e, L = s_s[a_[i] : b_[i]], e_s[a_[i] : b_[i]], l_s[a_[i] : b_[i]]
        m = (s >= lo[i]) & (e <= lo[i] + ell)
        cnt[i] = m.sum()
        den[i] = (1.0 / L[m]).sum()
    eff = float((fl * np.maximum(ell - Lv + 1, 0)).sum())
    print(
        f"    ell={ell:>5}: count {cnt.mean():>8.3f}  rho*E[(ell-L+1)+] = {RHO * eff:>8.3f}"
        f"   |  sum1/L {den.mean():>7.4f}  vs rho*ell = {RHO * ell:>7.4f}"
        f"  ratio {den.mean() / (RHO * ell):.4f}"
    )
