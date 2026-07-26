"""ADJ — does the RC1 'post-pin graft' arm create NEW f_g=1 vertex pins?

Under `graft_mode='post_pin'` the graft's DENSITY leaves the `_pin_v` normalizer but its PRECISION
(`_spc`) still enters `tpp`, so the normalizer's presence test `np.where(pp_ > 0, p, op)` sees an RNA leg
of *exactly zero* whenever the relayed RNA is 0 — i.e. it renormalizes the gDNA leg onto the f_g = 1
vertex (RC3 defect #4's shape) before the graft is added back. Count it, and measure the delivered mo_g.
"""

from __future__ import annotations

import sys

import numpy as np

sys.path.insert(0, "/Users/mkiyer/proj/rigel/scratchpad")
import rc1_replay as R  # noqa: E402

_EPS = 1e-9
CONDS = [
    "gdna_gdna300_ss_0.99_nrna_none_capture_on",
    "gdna_gdna300_ss_0.99_nrna_present_capture_on",
    "gdna_none_ss_0.99_nrna_none_capture_on",
]


def probe(ctx, arrs, src, valid, df, sf, dst_face_v, src_face_v):
    """Recompute the pre-pin legs exactly as `_transport` does, for one direction."""
    rg, rp, rn, pg, pp, pn, mg, mp, mn, tau = arrs
    from rigel.calibration.enrichment_frame import transfer_logvar

    sfv = src_face_v[src]
    framed = valid & (sfv > _EPS) & (dst_face_v > _EPS)
    r = np.where(framed, dst_face_v / np.maximum(sfv, _EPS), np.where(valid, 1.0, 0.0))
    graft = ctx.ex_a & ctx.is_bnd[src] & valid
    gp = np.where(graft, ctx.spl_p_f[sf][src], 0.0)
    gn = np.where(graft, ctx.spl_n_f[sf][src], 0.0)
    s2t = transfer_logvar(ctx.logvar_tot, ctx.logvar_tot[src], graft)

    def _dv(p):
        return np.where(valid & (p[src] > 0.0), 1.0 / (1.0 / np.maximum(p[src], _EPS) + s2t), 0.0)

    tpp, tpn, tpg = _dv(pp), _dv(pn), _dv(pg)
    _sp = np.where(graft, ctx.SP[sf][src], 0.0)
    _sn = np.where(graft, ctx.SN[sf][src], 0.0)
    _s2 = np.where(np.isfinite(s2t), s2t, 0.0)
    tpp = tpp + np.where(_sp > _EPS, _sp / (1.0 + _sp * _s2), 0.0)
    tpn = tpn + np.where(_sn > _EPS, _sn / (1.0 + _sn * _s2), 0.0)
    # relayed-only RNA legs (post_pin arm) vs shipped (graft folded in before the reframe)
    rel_p, rel_n = rp[src] * r, rn[src] * r
    ship_p, ship_n = (rp[src] + gp) * r, (rn[src] + gn) * r
    rel_p = np.where(ctx.fp_a, rel_p, 0.0)
    rel_n = np.where(ctx.fn_a, rel_n, 0.0)
    ship_p = np.where(ctx.fp_a, ship_p, 0.0)
    ship_n = np.where(ctx.fn_a, ship_n, 0.0)
    tpp = np.where(ctx.fp_a, tpp, 0.0)
    tpn = np.where(ctx.fn_a, tpn, 0.0)
    tpg = np.where(valid, tpg, 0.0)
    # the pin's presence test, post_pin arm
    dead_rna = (tpp + tpn > 0.0) & (rel_p + rel_n <= _EPS)
    live_g = tpg > 0.0
    return graft, dead_rna & live_g, ship_p + ship_n, rel_p + rel_n


def main():
    for cond in CONDS:
        ctx = R.Ctx(cond)
        fwd, bwd = ctx.build_relays(graft_mode="dst")
        _, rho_lf, rho_rf = ctx.rho_faces(ctx.fg_init)
        tot_g = tot_bad = 0
        bad_mass = 0.0
        for arrs, src, valid, df, sf, dfv, sfv in (
            (fwd, ctx.sl, ctx.vl, 0, 1, rho_lf, rho_rf),
            (bwd, ctx.sr, ctx.vr, 1, 0, rho_rf, rho_lf),
        ):
            graft, bad, _, _ = probe(ctx, arrs, src, valid, df, sf, dfv, sfv)
            tot_g += int(graft.sum())
            tot_bad += int((bad & graft).sum())
            bad_mass += float(ctx.mass[bad & graft].sum())
        print(f"{cond:<46} graft msgs={tot_g:6d}  post-pin RNA-leg-dead={tot_bad:5d} "
              f"({100 * tot_bad / max(tot_g, 1):5.1f}%)  dst mass={bad_mass:12,.0f}")


if __name__ == "__main__":
    main()
