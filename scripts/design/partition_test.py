"""Should a fragment's 1/L be PARTITIONED across the boundaries it crosses?

Three candidate rules, each evaluated at every boundary against the KNOWN true density rho:

  R1  no partition      each crossed boundary gets  1/L
  R2  equal partition   each crossed boundary gets  (1/L) * (1/K)     <- the current tool's rule
  R3  overlap-weighted  boundary between region a,b gets (1/L) * (overlap_a + overlap_b)/(2L)

Ground truth: molecules of length L land with START density rho per template position, uniform.
Everything is exact -- no approximation anywhere -- so any departure from rho is the rule's own bias.
"""

import numpy as np

rng = np.random.default_rng(3)
T = 400_000
RHO = 0.05
MU, SD = 200.0, 60.0
N_REP = 40


def run(spacing):
    cuts = np.arange(0, T + 1, spacing)  # region boundaries; boundaries sit at cuts[1:-1]
    n_boundary = cuts.size - 2
    acc = {k: np.zeros(n_boundary) for k in ("R1", "R2", "R3")}
    cnt = np.zeros(n_boundary)
    for _ in range(N_REP):
        n = rng.poisson(RHO * T)
        lens = np.clip(rng.normal(MU, SD, n).round().astype(np.int64), 30, 600)
        starts = rng.integers(0, T - lens)
        ends = starts + lens  # half-open
        # boundaries strictly inside the fragment: cuts[j] with start < cuts[j] < end
        lo = np.searchsorted(cuts, starts, side="right")  # first cut > start
        hi = np.searchsorted(cuts, ends, side="left")  # first cut >= end
        K = np.maximum(hi - lo, 0)
        for s, e, L, a, b, k in zip(starts, ends, lens, lo, hi, K):
            if k == 0:
                continue
            idx = np.arange(a, b) - 1  # boundary index (cut j -> boundary j-1)
            ok = (idx >= 0) & (idx < n_boundary)
            idx = idx[ok]
            if idx.size == 0:
                continue
            cnt[idx] += 1.0
            acc["R1"][idx] += 1.0 / L
            acc["R2"][idx] += (1.0 / L) / k
            # overlap-weighted: the two region-pieces flanking this boundary, over 2L
            cj = cuts[a:b][ok]
            left_ov = cj - np.maximum(s, cj - spacing)
            right_ov = np.minimum(e, cj + spacing) - cj
            acc["R3"][idx] += (1.0 / L) * (left_ov + right_ov) / (2.0 * L)
    for k in acc:
        acc[k] /= N_REP
    cnt /= N_REP
    interior = slice(5, n_boundary - 5)
    return {k: float(np.mean(v[interior])) / RHO for k, v in acc.items()}, float(
        np.mean(cnt[interior])
    )


print(f"truth rho = {RHO};  fragment length ~ N({MU:.0f},{SD:.0f})")
print(f"{'region spacing':>13} {'mean K':>8} | {'R1 no-part':>11} {'R2 1/K':>9} {'R3 overlap':>11}")
print("-" * 62)
for s in (50, 100, 200, 400, 1000, 4000):
    r, c = run(s)
    meanK = c / (RHO * MU) * (MU / s) if False else None
    print(
        f"{s:>13} {'~' + str(round(MU / s, 2)):>8} | "
        f"{r['R1']:>11.4f} {r['R2']:>9.4f} {r['R3']:>11.4f}"
    )
print("\n(values are estimate / truth; 1.0000 = unbiased)")
print("R2 and R3 depend on the REGION SPACING -- a property of the annotation, not the sample.")
