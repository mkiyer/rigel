"""ADJ — landmine check for the post-pin graft: does the mature density escape the fp_a/fn_a strand gate?

In `rc1_replay`'s `post_pin` arm the graft is added AFTER the `np.where(fp_a, tp, 0.0)` structural gate, so a
strand that the destination does not admit can still receive mature density. Its precision is gated (so ψ's
rna_imp is inert), but `cR = cp + cn` — hence `mo_R` and the single-λ message — is NOT. Count the exposure.
"""

from __future__ import annotations

import sys

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
import rc1_replay as R  # noqa: E402

_EPS = 1e-9
CONDS = [
    "gdna_gdna300_ss_0.99_nrna_none_capture_on",
    "gdna_gdna300_ss_0.50_nrna_none_capture_on",
    "gdna_gdna100_ss_0.99_nrna_present_capture_on",
    "gdna_none_ss_0.50_nrna_none_capture_off",
]


def main():
    for cond in CONDS:
        ctx = R.Ctx(cond)
        n_bad = 0
        bad_mass = 0.0
        n_graft = 0
        for src, valid, sf in ((ctx.sl, ctx.vl, 1), (ctx.sr, ctx.vr, 0)):
            graft = ctx.ex_a & ctx.is_bnd[src] & valid
            gp = np.where(graft, ctx.spl_p_f[sf][src], 0.0)
            gn = np.where(graft, ctx.spl_n_f[sf][src], 0.0)
            bad = ((gp > _EPS) & ~ctx.fp_a) | ((gn > _EPS) & ~ctx.fn_a)
            n_graft += int((graft & ((gp > _EPS) | (gn > _EPS))).sum())
            n_bad += int(bad.sum())
            bad_mass += float(ctx.mass[bad].sum())
        print(f"{cond:<46} live grafts={n_graft:5d}  strand-DEAD graft density={n_bad:4d}"
              f"  dst mass={bad_mass:12,.0f}")


if __name__ == "__main__":
    main()
