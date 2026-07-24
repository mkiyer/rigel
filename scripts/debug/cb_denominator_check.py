"""A3 — which denominator makes  c_b = log1p(S_B/D_B)  recover the TRUE mature dilution?

c_b converts an EXON's mature-INCLUSIVE composition into the BOUNDARY's mature-FREE crossing composition:

    exon unspliced pool     = {gDNA, nascent, mature}      (within-exon mature spans no junction)
    boundary unspliced pool = {gDNA, nascent}              (mature is SPLICED -> the spliced channel)

    f_g^B = f_g^E * (D_B + S_B)/D_B      =>   c_b = log1p(S_B/D_B)
    D_B = the boundary's unspliced crossing density,  S_B = its spliced (mature) density.

TRUTH being targeted:   r_true = rho_mu / (rho_g + rho_nu).

The deposit rule (accumulator/00_design.md 4.3, verified in tests/native/_accumulator_reference.py) gives the
expected deposited mass per unit density in CLOSED FORM, so this needs no simulation:

    unspliced crossing, per FACE f :  rho * E_fl[min(l, R_f)] / 2        (mass lands on BOTH faces)
    spliced crossing,   ONE face   :  rho_mu * eff_spl(R_exon)           (mass lands on ONE face only)
        eff_spl = E_rna[min(l,R)]/2   if the transcript CONTINUES past the flank  (A1/A2)
                = E_rna[min(l,R)^2/(2l)]  (half-triangle) if it TERMINATES there

Candidate denominators for S_B:
    (a) SUMMED-EFF   spl_mass / (eff_spl_l + eff_spl_r)     <- what production does today
    (b) PER-FACE     spl_mass / eff_spl_face                <- the face that can actually receive the mass

and D_B = (m_l + m_r) / (EG_l + EG_r), where EG uses the **gDNA** FL even though the mass it divides also
contains nascent RNA (a different FL) -- a separate frame question this script also quantifies.

Usage:  OMP_NUM_THREADS=1 python scripts/debug/cb_denominator_check.py
"""
from __future__ import annotations

import numpy as np

from rigel.calibration.effective_length import boundary_side_eff_length, spliced_side_eff_length


def gauss_pmf(mu, sd, n=1200):
    x = np.arange(float(n))
    p = np.exp(-0.5 * ((x - mu) / sd) ** 2)
    return p / p.sum()


def analyse(R_e, R_i, rho_g, rho_nu, rho_mu, gpmf, rpmf, continues=True):
    """Exact expected deposits, then each candidate's implied r = S_B/D_B."""
    Re, Ri = np.array([float(R_e)]), np.array([float(R_i)])
    # per-face two-sided crossing opportunity, per component FL
    Eg_e, Eg_i = boundary_side_eff_length(gpmf, Re)[0], boundary_side_eff_length(gpmf, Ri)[0]
    Er_e, Er_i = boundary_side_eff_length(rpmf, Re)[0], boundary_side_eff_length(rpmf, Ri)[0]
    # one-sided spliced opportunity on the EXON face
    esp_e = (boundary_side_eff_length(rpmf, Re) if continues else spliced_side_eff_length(rpmf, Re))[0]
    # the INTRON face carries no spliced mass, but production still adds its opportunity to the divisor
    esp_i = (boundary_side_eff_length(rpmf, Ri) if continues else spliced_side_eff_length(rpmf, Ri))[0]

    # deposited unspliced mass per face = gDNA (gDNA FL) + nascent (RNA FL)
    m_e = rho_g * Eg_e + rho_nu * Er_e
    m_i = rho_g * Eg_i + rho_nu * Er_i
    spl = rho_mu * esp_e  # one-sided

    D_B = (m_e + m_i) / (Eg_e + Eg_i)  # production: gDNA-FL frame for a mixed-FL mass
    r_true = rho_mu / (rho_g + rho_nu)
    r_sum = (spl / (esp_e + esp_i)) / D_B
    r_face = (spl / esp_e) / D_B
    # what D_B would be with a composition-correct divisor (isolates the mixed-FL frame error)
    D_ideal = rho_g + rho_nu
    r_face_ideal = (spl / esp_e) / D_ideal
    return r_true, r_sum, r_face, r_face_ideal, D_B / D_ideal


gp, rp = gauss_pmf(300, 60), gauss_pmf(200, 50)
print("gDNA FL ~ N(300,60);  RNA FL ~ N(200,50);  rho_g=0.5  rho_nu=0.3  rho_mu=1.0  (r_true = 1.25)")
print("'c_b' columns are log1p(r) — the actual additive correction applied to the log-fraction mode.\n")
hdr = (
    f"{'R_exon':>7} {'R_intr':>7} {'cont':>5} | {'r_true':>7} {'r_SUMMED':>9} {'r_PERFACE':>10} |"
    f" {'c_b true':>8} {'c_b sum':>8} {'c_b face':>8} | {'D_B/D_true':>10}"
)
print(hdr)
print("-" * len(hdr))
for cont in (True, False):
    for R_e in (100, 300, 1000, 3000):
        for R_i in (2000,):
            t, s, f, fi, dratio = analyse(R_e, R_i, 0.5, 0.3, 1.0, gp, rp, continues=cont)
            print(
                f"{R_e:>7} {R_i:>7} {str(cont):>5} | {t:>7.3f} {s:>9.3f} {f:>10.3f} |"
                f" {np.log1p(t):>8.4f} {np.log1p(s):>8.4f} {np.log1p(f):>8.4f} | {dratio:>10.4f}"
            )
    print()

print("Isolating the two error sources at R_exon=1000, R_intron=2000, continues=True:")
t, s, f, fi, dratio = analyse(1000, 2000, 0.5, 0.3, 1.0, gp, rp, True)
print(f"   truth                                r = {t:.4f}   c_b = {np.log1p(t):.4f}")
print(f"   PER-FACE S_B, IDEAL D_B (= rho_g+rho_nu)  r = {fi:.4f}   c_b = {np.log1p(fi):.4f}   <- frame-error-free")
print(f"   PER-FACE S_B, production D_B          r = {f:.4f}   c_b = {np.log1p(f):.4f}")
print(f"   SUMMED   S_B, production D_B (TODAY)  r = {s:.4f}   c_b = {np.log1p(s):.4f}")
print(f"\n   D_B / (rho_g+rho_nu) = {dratio:.4f}  <- the MIXED-FL frame error: D_B divides a mass containing")
print("      nascent (RNA FL) by the gDNA-FL opportunity, so it does not equal rho_g+rho_nu.")
