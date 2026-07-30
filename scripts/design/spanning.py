"""Verify a node's THREE populations (contained, spanning, neither) against brute-force simulation."""

import numpy as np

rng = np.random.default_rng(19)
T, RHO, MU, SD = 300_000, 0.05, 200.0, 60.0
L = np.arange(1, 1201)
w_ = np.exp(-0.5 * ((L - MU) / SD) ** 2)
w_ /= w_.sum()

n = rng.poisson(RHO * T)
lens = np.clip(rng.normal(MU, SD, n).round().astype(np.int64), 30, 600)
starts = rng.integers(0, T - lens)
ends = starts + lens

print(
    f"{'node len':>9} {'contained obs':>14} {'pred':>9} | {'encompass obs':>14} {'pred':>9} | {'neither':>8}"
)
for ell in (25, 50, 100, 150, 200, 300, 600, 1200):
    a = rng.integers(1000, T - ell - 1000, 3000)
    b = a + ell
    c = e = 0
    for lo, hi in zip(a, b):
        m = (starts < hi) & (ends > lo)
        s, en = starts[m], ends[m]
        c += int(((s >= lo) & (en <= hi)).sum())
        e += int(((s <= lo - 1) & (en - 1 >= hi)).sum())
    c /= a.size
    e /= a.size
    pc = RHO * float((w_ * np.maximum(ell - L + 1, 0)).sum())
    pe = RHO * float((w_ * np.maximum(L - ell - 1, 0)).sum())
    print(f"{ell:>9} {c:>14.4f} {pc:>9.4f} | {e:>14.4f} {pe:>9.4f} | {'w=l+1 only':>8}")
