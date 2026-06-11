"""FL-consistency diagnostic (Phase-3 prep, accuracy) — expose the f_b·M_region bias.

The splice-junction partition currently sets a region's gDNA fraction to the boundary fraction
`f_b` directly (`splice_junction.region_splice_gdna_frac`). But `f_b` is a *molecular density*
fraction (the boundary divides crossing counts by the **crossing** eff-length `fl_mean`, which
cancels to `d_g/(d_g+d_r)`), while it is then applied to the region's contained **count**
`M_region`. The count fraction of the contained mass is

    g_true = d_g·E^g_region(L) / ( d_g·E^g_region(L) + d_r·E^r_region(L) )

with `E^{g,r}_region(L) = region_eff_length(L, fl_pmf_{g,r})`. Writing `r = E^g_region/E^r_region`
and `f_b = d_g/(d_g+d_r)`:

    g_corr = (f_b·E^g_region) / (f_b·E^g_region + (1−f_b)·E^r_region)   ==  g_true   (exact identity)
    g_bare = f_b

So the bare form is correct iff `r = 1` (equal gDNA/RNA region eff-lengths). They diverge when the
gDNA and RNA FL distributions differ AND the exon is short (short L weights the FL-containability gap;
for L ≫ FL the eff-lengths both → L − fl_mean and r → 1). This probe sweeps exon length × FL gap with
the REAL eff-length functions and reports the bare-form bias `g_bare − g_corr`, with no sampling noise
(isolating the *structural* error). Stage (b) — a realized read-deposit scenario — follows only if this
shows a material bias.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from rigel.calibration.effective_length import boundary_eff_length, region_eff_length

MAXLEN = 1000


def gaussian_pmf(mean, sd, floor=50, maxlen=MAXLEN):
    """Discrete fragment-length pmf ~ N(mean, sd) truncated to [floor, maxlen]."""
    ell = np.arange(maxlen + 1, dtype=np.float64)
    p = np.exp(-0.5 * ((ell - mean) / sd) ** 2)
    p[: int(floor)] = 0.0  # no fragments shorter than the floor (~2·readlen)
    return p / p.sum()


def g_corr(f_b, eg, er):
    """The eff-length-consistent region count fraction (== g_true under uniform density)."""
    num = f_b * eg
    den = num + (1.0 - f_b) * er
    return np.where(den > 0, num / den, np.nan)


# Scenarios use the EXACT benchmark FL parameters (RNA frag_mean=250/std=50/min=80; gDNA per config)
# so the bias maps directly onto the suites. The standard suites set gDNA 350/100 (LONGER than RNA);
# gdna_shortfl_5mb sets gDNA 100/20 (much SHORTER). Both gaps are real and exercised.
RNA = gaussian_pmf(250, 50, floor=80)
SCENARIOS = [
    ("control: gDNA == RNA (250/50)", gaussian_pmf(250, 50, floor=80)),
    ("BENCHMARK std: gDNA 350/100 (longer)", gaussian_pmf(350, 100, floor=100)),
    ("gdna_shortfl: gDNA 100/20 (much shorter)", gaussian_pmf(100, 20, floor=50)),
]
LENGTHS = np.array([60, 80, 100, 120, 150, 200, 250, 300, 400, 600, 800, 1000], dtype=float)
F_B_GRID = [0.2, 0.5, 0.8]

er_region = region_eff_length(LENGTHS, RNA)
er_boundary = boundary_eff_length(RNA)
print(f"RNA: fl_mean(boundary)={er_boundary:.1f};  human median exon ≈ 120 bp\n")

fig, ax = plt.subplots(1, 2, figsize=(13, 5))
for label, gdna in SCENARIOS:
    eg_region = region_eff_length(LENGTHS, gdna)
    ratio = np.where(er_region > 0, eg_region / np.maximum(er_region, 1e-9), np.nan)
    print(f"=== {label}   fl_mean_g={boundary_eff_length(gdna):.1f} ===")
    print(f"{'L':>6} {'E^g_reg':>9} {'E^r_reg':>9} {'r=Eg/Er':>9} "
          f"{'g_bare':>7} {'g_corr':>7} {'bias@.5':>8}  {'|bias| @ f_b=.2/.5/.8':>22}")
    for i, L in enumerate(LENGTHS):
        eg, er = eg_region[i], er_region[i]
        biases = [abs(fb - g_corr(fb, eg, er)) for fb in F_B_GRID]
        gc5 = g_corr(0.5, eg, er)
        bias5 = 0.5 - gc5 if np.isfinite(gc5) else np.nan
        flag = "  <-- exon-scale" if 80 <= L <= 250 else ""
        print(f"{L:>6.0f} {eg:>9.1f} {er:>9.1f} {ratio[i]:>9.2f} "
              f"{0.5:>7.2f} {gc5:>7.3f} {bias5:>+8.3f}  "
              f"{biases[0]:>6.3f}/{biases[1]:.3f}/{biases[2]:.3f}{flag}")
    # bias vs L at f_b=0.5 (signed: + = bare OVER-estimates gDNA)
    bias_curve = np.array([0.5 - g_corr(0.5, eg_region[i], er_region[i]) for i in range(len(LENGTHS))])
    ax[0].plot(LENGTHS, bias_curve, marker="o", label=label)
    ax[1].plot(LENGTHS, ratio, marker="o", label=label)
    print()

for a, (ttl, yl) in zip(ax, [("(A) bare-form gDNA-fraction bias (f_b=0.5)\n+ = bare over-calls gDNA",
                              "g_bare − g_corr"),
                             ("(B) region eff-length ratio r = E^g/E^r\n(=1 ⇒ no bias)", "r")]):
    a.axhline(0.0 if "bias" in ttl else 1.0, color="k", lw=0.6, ls=":")
    a.axvspan(80, 250, color="0.9", label="typical exon scale")
    a.set_xlabel("exon (region) length L (bp)")
    a.set_ylabel(yl)
    a.set_title(ttl)
    a.set_xscale("log")
    a.legend(fontsize=7)
fig.tight_layout()
fig.savefig("/tmp/fl_consistency.png", dpi=120)
print("saved /tmp/fl_consistency.png")
