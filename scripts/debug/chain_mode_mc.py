"""Monte-Carlo validation of the message-propagation MODE arithmetic on the

        intron region  <->  intron|exon boundary  <->  exon region

chain, with gDNA + nascent RNA (unspliced, gene-body) + mature RNA (spliced, exon).
Capture enrichment E on any fragment overlapping the exon (positional, NA-agnostic).

The point: derive/validate the DENSITY-currency message math. Densities ρ_c = C_c/E_c
(mass over the per-component eff-length) are the frame-invariant currency; FRACTIONS are
frame-dependent. We validate four claims (see the derivation doc):

  C1  intron -> boundary  : the FL-frame log-odds shift predicts the boundary's UNSPLICED f_g
                            (mature-free cliff; capture cancels).
  C2  frame-invariance    : ρ_mature at the boundary (spliced) == ρ_mature contained in the exon;
                            ρ_g / ρ_nasc are shared (no capture cliff) between boundary & exon.
  C3  boundary -> exon    : with the boundary's spliced added (ρ_r = ρ_nasc + ρ_mature) the
                            composition shift predicts the exon's f_g (which INCLUDES mature).
  C4  exon -> boundary    : the ABSORPTION — subtract the boundary's own spliced ρ_mature from the
                            exon RNA message to recover ρ_nasc, then predict the boundary's f_g.

Run inside the `rigel` env.  OMP_NUM_THREADS=1.
"""

from __future__ import annotations

import numpy as np

from rigel.calibration.effective_length import (
    boundary_eff_length,  # crossing count eff-len for the NODE = fl_mean
    region_eff_length,
    spliced_side_eff_length,
)


# ----------------------------------------------------------------------------- FL pmfs
def make_pmf(kind: str, m: float, s: float, L: int = 1200) -> np.ndarray:
    x = np.arange(1, L)
    if kind == "gauss":
        p = np.exp(-0.5 * ((x - m) / s) ** 2)
    elif kind == "gamma":
        k = (m / s) ** 2
        th = s * s / m
        p = x ** (k - 1) * np.exp(-x / th)
    elif kind == "bimodal":
        p = np.exp(-0.5 * ((x - m) / s) ** 2) + 0.7 * np.exp(-0.5 * ((x - (m + 200)) / (s * 0.6)) ** 2)
    elif kind == "unif":
        p = ((x >= m - s) & (x <= m + s)).astype(float)
    else:
        raise ValueError(kind)
    p = np.maximum(p, 0.0)
    pmf = np.zeros(L)
    pmf[1:] = p / p.sum()
    return pmf


# ----------------------------------------------------------------------------- the simulation
def simulate(
    gpmf,
    rpmf,
    *,
    Li=3000.0,   # intron length
    Le=3000.0,   # exon length
    d_g=0.40,    # gDNA molecular density (per bp)
    d_nasc=0.50, # nascent RNA molecular density over the gene body
    d_mat=0.90,  # mature RNA molecular density over the exon (0 => no splice junction)
    E=1.0,       # capture enrichment on exon-overlapping fragments
    capture_mode="binary",  # "binary" (per-node scalar) | "proportional" (per-fragment overlap weight)
    Nbase=8_000_000,
    seed=1,
):
    """Deposit fragments and return per-node per-component MASSES (weighted counts).

    Geometry (genomic axis): intron [0, Li), exon [Li, Li+Le), boundary point at x=Li.
    A previous exon sits off-stage to the LEFT (across the intron) so mature fragments can
    span the acceptor at Li (=> one-sided spliced deposit on the exon flank).
    """
    rng = np.random.default_rng(seed)
    Ltot = Li + Le
    maxfl = max(len(gpmf), len(rpmf))

    def capture_w(s, end):
        """Capture weight for a fragment [s,end).

        binary      — E if it overlaps the exon [Li, Ltot) at all, else 1 (a per-NODE scalar that cancels in
                      density ratios ⇒ the boundary⟷exon "no cliff").
        proportional — 1 + (E−1)·(exon-overlap bp / ℓ): a per-FRAGMENT weight. A boundary-crossing fragment lies
                      only partly in the exon ⇒ under-enriched vs a fully-contained exon fragment, so the cliff
                      no longer fully cancels (the differential mature-vs-unspliced cliff — see the doc §6a)."""
        if capture_mode == "binary":
            return np.where((end > Li) & (s < Ltot), E, 1.0)
        overlap = np.clip(np.minimum(end, Ltot) - np.maximum(s, Li), 0.0, None)
        ell = np.maximum(end - s, 1e-9)
        return 1.0 + (E - 1.0) * overlap / ell

    # ---- UNSPLICED components (gDNA over all genome; nascent over gene body) ---------------
    def run_unspliced(pmf, dens, lo, hi):
        n = int(Nbase * dens)
        s = rng.uniform(lo - maxfl, hi, n)
        ell = (rng.choice(len(pmf), size=n, p=pmf)).astype(float)
        end = s + ell
        w = capture_w(s, end)
        cont_i = (s >= 0) & (end <= Li)          # contained in intron
        cont_e = (s >= Li) & (end <= Ltot)       # contained in exon
        cross = (s < Li) & (end > Li)            # crosses the boundary point Li (UNSPLICED)
        return {
            "cont_i": float((cont_i * w).sum()),
            "cont_e": float((cont_e * w).sum()),
            "cross": float((cross * w).sum()),
        }

    g = run_unspliced(gpmf, d_g, 0.0, Ltot)          # gDNA everywhere
    nasc = run_unspliced(rpmf, d_nasc, 0.0, Ltot)    # nascent over the gene body [0, Ltot)

    # ---- MATURE (spliced) RNA: fragments on the mature transcript, exon + acceptor junction --
    # Mature transcript coords: [prev_exon(len Pe) | this_exon(len Le) | ...]; the prev|this
    # junction maps to genomic Li. A mature fragment either (a) lies within this_exon => genomic
    # contained in exon (adds to exon RNA), or (b) spans the acceptor => one-sided spliced deposit
    # on the exon flank of the boundary at Li (deposit = exon-overlap slice share a/ℓ).
    Pe = 3000.0  # prev-exon length (large flank)
    mat = {"cont_e": 0.0, "spliced": 0.0}
    if d_mat > 0.0:
        n = int(Nbase * d_mat)
        # start uniformly over [prev_exon .. this_exon], mature coords 0..(Pe+Le)
        ms = rng.uniform(-maxfl, Pe + Le, n)
        ell = (rng.choice(len(rpmf), size=n, p=rpmf)).astype(float)
        me = ms + ell
        # this_exon occupies mature coords [Pe, Pe+Le); the acceptor junction is at mature coord Pe.
        cont_e = (ms >= Pe) & (me <= Pe + Le)                      # wholly within this exon
        spans = (ms < Pe) & (me > Pe) & (me <= Pe + Le)           # spans acceptor, ends inside exon
        # exon-side overlap length a for a spanning fragment = me - Pe  (bp in this exon)
        a = np.clip(me - Pe, 0.0, Le)
        # capture: contained mature is fully in the exon (weight E in both modes); a spanning fragment is
        # captured by its THIS-exon overlap only (proportional => under-captured ~a/ℓ; binary => E).
        if capture_mode == "binary":
            w_cont = E
            w_spl = np.where(spans, E, 0.0)
        else:
            w_cont = E  # contained ⇒ overlap == ℓ ⇒ 1+(E-1)·1 == E
            w_spl = np.where(spans, 1.0 + (E - 1.0) * (a / ell), 0.0)
        mat["cont_e"] = float((cont_e).sum()) * w_cont
        # one-sided spliced deposit = Σ (a/ℓ)·weight  (END slice, intron flank never credited)
        dep = np.where(spans, (a / ell) * w_spl, 0.0)
        mat["spliced"] = float(dep.sum())

    return {"g": g, "nasc": nasc, "mat": mat, "Li": Li, "Le": Le, "Pe": Pe}


# ----------------------------------------------------------------------------- densities & truth
def node_densities(sim, gpmf, rpmf):
    """Per-node per-component DENSITIES ρ = mass / eff-length, and the TRUE (measured) fractions."""
    Li, Le = sim["Li"], sim["Le"]
    Eg_reg_i = region_eff_length(np.array([Li]), gpmf)[0]
    Er_reg_i = region_eff_length(np.array([Li]), rpmf)[0]
    Eg_reg_e = region_eff_length(np.array([Le]), gpmf)[0]
    Er_reg_e = region_eff_length(np.array([Le]), rpmf)[0]
    Eg_bnd = boundary_eff_length(gpmf)   # = fl_mean(gpmf)
    Er_bnd = boundary_eff_length(rpmf)
    Espl = spliced_side_eff_length(rpmf, np.array([Le]))[0]  # one-sided, exon flank R=Le

    g, nasc, mat = sim["g"], sim["nasc"], sim["mat"]

    # intron region (contained): gDNA + nascent
    rho_g_i = g["cont_i"] / Eg_reg_i
    rho_r_i = nasc["cont_i"] / Er_reg_i
    fg_i = rho_g_i * Eg_reg_i / (rho_g_i * Eg_reg_i + rho_r_i * Er_reg_i)  # == count fraction

    # boundary (unspliced crossing): gDNA + nascent ; plus one-sided spliced mature
    rho_g_b = g["cross"] / Eg_bnd
    rho_nasc_b = nasc["cross"] / Er_bnd
    rho_mat_b = (mat["spliced"] / Espl) if mat["spliced"] > 0 else 0.0
    fg_b = rho_g_b * Eg_bnd / (rho_g_b * Eg_bnd + rho_nasc_b * Er_bnd)  # UNSPLICED f_g (mature excluded)

    # exon region (contained): gDNA + (nascent + mature)
    rho_g_e = g["cont_e"] / Eg_reg_e
    rho_nasc_e = nasc["cont_e"] / Er_reg_e
    rho_mat_e = (mat["cont_e"] / Er_reg_e) if mat["cont_e"] > 0 else 0.0
    rho_r_e = rho_nasc_e + rho_mat_e
    C_g_e = g["cont_e"]
    C_r_e = nasc["cont_e"] + mat["cont_e"]
    fg_e = C_g_e / (C_g_e + C_r_e)  # exon f_g INCLUDES mature in the RNA

    return dict(
        Eg_reg_i=Eg_reg_i, Er_reg_i=Er_reg_i, Eg_reg_e=Eg_reg_e, Er_reg_e=Er_reg_e,
        Eg_bnd=Eg_bnd, Er_bnd=Er_bnd, Espl=Espl,
        rho_g_i=rho_g_i, rho_r_i=rho_r_i, fg_i=fg_i,
        rho_g_b=rho_g_b, rho_nasc_b=rho_nasc_b, rho_mat_b=rho_mat_b, fg_b=fg_b,
        rho_g_e=rho_g_e, rho_nasc_e=rho_nasc_e, rho_mat_e=rho_mat_e, rho_r_e=rho_r_e, fg_e=fg_e,
    )


def _fg_from_shift(fg_src, log_shift):
    odds = fg_src / (1.0 - fg_src)
    return 1.0 / (1.0 + 1.0 / (odds * np.exp(log_shift)))


def _fg_from_densities(rho_g, rho_r, Eg, Er):
    """recipient-side composition shift: f_g = ρ_g·E_g / (ρ_g·E_g + ρ_r·E_r)."""
    Mg = rho_g * Eg
    Mr = rho_r * Er
    return Mg / (Mg + Mr)


# ----------------------------------------------------------------------------- the four claims
def check(name, gk, rk, **kw):
    gpmf = make_pmf(*gk)
    rpmf = make_pmf(*rk)
    sim = simulate(gpmf, rpmf, **kw)
    d = node_densities(sim, gpmf, rpmf)

    # C1: intron -> boundary FL-frame shift predicts the boundary UNSPLICED f_g
    shift = np.log(d["Eg_bnd"] / d["Eg_reg_i"]) - np.log(d["Er_bnd"] / d["Er_reg_i"])
    c1 = _fg_from_shift(d["fg_i"], shift)

    # C2: frame-invariance of ρ_mature (spliced-frame vs exon-contained-frame) & shared g/nasc scale
    inv_mat = (d["rho_mat_b"], d["rho_mat_e"])          # should be equal
    inv_g = (d["rho_g_b"], d["rho_g_e"])                # should be equal (no cliff boundary<->exon)
    inv_nasc = (d["rho_nasc_b"], d["rho_nasc_e"])       # should be equal

    # C3: boundary -> exon, add the boundary's spliced ρ_mat into the RNA density
    rho_r_msg = d["rho_nasc_b"] + d["rho_mat_b"]
    c3 = _fg_from_densities(d["rho_g_b"], rho_r_msg, d["Eg_reg_e"], d["Er_reg_e"])

    # C4: exon -> boundary ABSORPTION. exon sends (ρ_g_e, ρ_r_e = nasc+mat); boundary subtracts
    #     its own spliced ρ_mat to recover ρ_nasc, then shifts into the crossing frame.
    rho_nasc_abs = max(d["rho_r_e"] - d["rho_mat_b"], 0.0)
    c4 = _fg_from_densities(d["rho_g_e"], rho_nasc_abs, d["Eg_bnd"], d["Er_bnd"])

    print(f"\n=== {name}   (E={kw.get('E',1.0)}, d_mat={kw.get('d_mat',0.9)}) ===")
    print(f"  eff-len: Eg_reg_i={d['Eg_reg_i']:.1f} Er_reg_i={d['Er_reg_i']:.1f} "
          f"Eg_bnd={d['Eg_bnd']:.1f} Er_bnd={d['Er_bnd']:.1f} Espl={d['Espl']:.1f}")
    print(f"  TRUTH:  fg_intron={d['fg_i']:.4f}  fg_boundary(unspl)={d['fg_b']:.4f}  fg_exon={d['fg_e']:.4f}")
    print(f"  C1 intron->bnd shift: pred={c1:.4f}  truth={d['fg_b']:.4f}  |Δ|={abs(c1-d['fg_b']):.4f}")
    print(f"  C2 ρ_mat  bnd={inv_mat[0]:.4f} exon={inv_mat[1]:.4f}  |Δ|={abs(inv_mat[0]-inv_mat[1]):.4f}")
    print(f"     ρ_g    bnd={inv_g[0]:.4f} exon={inv_g[1]:.4f}  |Δ|={abs(inv_g[0]-inv_g[1]):.4f}")
    print(f"     ρ_nasc bnd={inv_nasc[0]:.4f} exon={inv_nasc[1]:.4f}  |Δ|={abs(inv_nasc[0]-inv_nasc[1]):.4f}")
    print(f"  C3 bnd->exon (+spliced): pred={c3:.4f}  truth={d['fg_e']:.4f}  |Δ|={abs(c3-d['fg_e']):.4f}")
    print(f"  C4 exon->bnd (absorb):   pred={c4:.4f}  truth={d['fg_b']:.4f}  |Δ|={abs(c4-d['fg_b']):.4f}")
    return d


def main():
    base = dict(Li=3000.0, Le=3000.0, d_g=0.40, d_nasc=0.50, d_mat=0.90)
    # vary FL_g, FL_r
    check("gDNA=RNA gauss200", ("gauss", 200, 40), ("gauss", 200, 40), **base)
    check("gDNA300>RNA150", ("gauss", 300, 50), ("gauss", 150, 30), **base)
    check("gDNA120<RNA350", ("gauss", 120, 25), ("gauss", 350, 60), **base)
    check("gamma FLs", ("gamma", 300, 90), ("gamma", 180, 60), **base)
    # vary capture (the cliff)
    check("gDNA300>RNA150 E=30", ("gauss", 300, 50), ("gauss", 150, 30), **{**base, "E": 30.0})
    check("gDNA300>RNA150 E=300", ("gauss", 300, 50), ("gauss", 150, 30), **{**base, "E": 300.0})
    # vary densities
    check("gDNA-rich", ("gauss", 250, 45), ("gauss", 180, 35), **{**base, "d_g": 2.0, "d_nasc": 0.3, "d_mat": 0.5})
    check("RNA-rich", ("gauss", 250, 45), ("gauss", 180, 35), **{**base, "d_g": 0.2, "d_nasc": 1.0, "d_mat": 3.0})
    # no splice junction (d_mat=0): C3/C4 collapse to the pure nascent case
    check("NO splice (d_mat=0)", ("gauss", 300, 50), ("gauss", 150, 30), **{**base, "d_mat": 0.0})
    check("NO splice E=30", ("gauss", 300, 50), ("gauss", 150, 30), **{**base, "d_mat": 0.0, "E": 30.0})

    # PROPORTIONAL capture (doc §6a): the boundary⟷exon cliff no longer fully cancels. The composition transfer
    # survives (C1 robust; C3/C4 exact at d_mat=0) but the mature add/subtract acquires a DIFFERENTIAL-cliff bias
    # ∝ mature fraction (C3/C4 degrade; worst when mature is present).
    prop = {**base, "E": 30.0, "capture_mode": "proportional"}
    check("PROP-capture E=30 (mature)", ("gauss", 300, 50), ("gauss", 150, 30), **prop)
    check("PROP-capture E=30 NO mature", ("gauss", 300, 50), ("gauss", 150, 30), **{**prop, "d_mat": 0.0})


if __name__ == "__main__":
    main()
