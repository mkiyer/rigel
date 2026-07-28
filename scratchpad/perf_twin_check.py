"""Bitwise equality of the scalar twins against their array originals — the standing check for Target A.

Far wider than the unit test (which is trimmed to keep pytest fast): a full product grid over every
argument corner PLUS a large randomized battery over the solver's operating ranges.
"""
from __future__ import annotations

import itertools
import math

import numpy as np

from rigel.calibration.enrichment_frame import (
    peel_continue_share,
    peel_continue_share_scalar,
    residual_level,
    residual_level_scalar,
)


def same(a, b):
    a, b = float(a), float(b)
    return np.float64(a).view(np.int64) == np.float64(b).view(np.int64) or (
        math.isnan(a) and math.isnan(b)
    )


CORNERS = (
    [0.0, 1e-13, 1e-12, 1.0, 12.0, 1e6, np.inf, np.nan],  # mass
    [0.0, 1e-13, 1.0, 33.0, 1e6, np.inf, np.nan],  # n_mass
    [-1.0, -0.0, 0.0, 1e-15, 0.4, 1e3, np.inf, np.nan],  # rho_g
    [0.0, 1e-9, 210.0, 1e5],  # E_g
    [0.0, 1e-12, 1e-13, 260.0],  # E_r
    [0.0, 1e-30, 1.0, 4.0, 1e6, 1e18, 1e300, np.inf, np.nan, -1.0],  # v_g
)

bad = 0
n = 0
for a in itertools.product(*CORNERS):
    arr, sca = residual_level(*a), residual_level_scalar(*a)
    n += 1
    for c in range(3):
        if not same(arr[c], sca[c]):
            bad += 1
            if bad < 12:
                print(f"MISMATCH slot{c} args={a} array={float(arr[c])!r} scalar={sca[c]!r}")
print(f"residual_level corner grid: {n} cases, {bad} mismatches")

rng = np.random.default_rng(7)
N = 200_000
draws = (
    np.exp(rng.uniform(-14, 14, N)) * (rng.random(N) > 0.05),
    np.floor(np.exp(rng.uniform(-1, 12, N))) * (rng.random(N) > 0.05),
    np.exp(rng.uniform(-16, 8, N)) * (rng.random(N) > 0.1),
    np.exp(rng.uniform(-2, 12, N)),
    np.exp(rng.uniform(-2, 12, N)),
    np.where(rng.random(N) > 0.1, np.exp(rng.uniform(-30, 6, N)), np.inf),
)
arr = residual_level(*draws)
bad2 = 0
for i in range(N):
    sca = residual_level_scalar(*(float(d[i]) for d in draws))
    for c in range(3):
        if not same(arr[c][i], sca[c]):
            bad2 += 1
            if bad2 < 12:
                print(f"RANDOM slot{c} i={i} args={[float(d[i]) for d in draws]}")
print(f"residual_level random battery: {N} cases, {bad2} mismatches")

pv = [0.0, -0.0, 1e-13, 1e-12, 1e-9, 0.5, 1.0, 1e6, np.inf, -np.inf, np.nan, -1.0]
bad3 = sum(
    0 if same(peel_continue_share(x, y), peel_continue_share_scalar(x, y)) else 1
    for x in pv
    for y in pv
)
u, v = np.exp(rng.uniform(-20, 12, N)), np.exp(rng.uniform(-20, 12, N))
ref = peel_continue_share(u, v)
bad3 += sum(
    0 if same(ref[i], peel_continue_share_scalar(float(u[i]), float(v[i]))) else 1 for i in range(N)
)
print(f"peel_continue_share: {len(pv) ** 2 + N} cases, {bad3} mismatches")

print("TOTAL MISMATCHES:", bad + bad2 + bad3)
raise SystemExit(1 if (bad + bad2 + bad3) else 0)
