"""Fix-1 prototype: does a VALLEY-FILLING mixture gDNA-density prior remove the collapse?

The AMBIG local objective (strand is ~flat for balanced AMBIG) reduces to:
    argmax_{f_g in (0,1)}  prior(log rho_g = log(f_g*d))  +  Jeffreys(-log(1-f_g))
We compare the current KDE prior against valley-free mixture priors, and plot the
"collapse curve": argmax f_g as a function of the node's TOTAL density d. The baseline
KDE should show a CLIFF (argmax f_g -> ~0) for d in the between-modes valley; a correct
mixture prior should NOT.

    OMP_NUM_THREADS=1 python mixture_prior_prototype.py <pass_trace.npz>
"""
import sys
import numpy as np

z = np.load(sys.argv[1])
tx = np.asarray(z["kde_train_x"], float)          # training log rho_g (solved non-AMBIG nodes)
tk = np.asarray(z["kde_train_kind"])
bw = float(z["kde_bw"])                            # 0.19
modes = z["kde_modes"]
mu_dep, mu_enr = float(modes[0, 0]), float(modes[1, 0])   # -4.39, +2.78 (log rho)
rho_dep, rho_enr = np.exp(mu_dep), np.exp(mu_enr)         # 0.012, 16.09
print(f"training n={tx.size}  bw={bw:.3f}  modes: dep log={mu_dep:.2f}(rho {rho_dep:.3f})  enr log={mu_enr:.2f}(rho {rho_enr:.2f})")

# ---------- prior variants: each returns log-prior(log_rho) ----------
W = tx.size
def kde_baseline(lr):                              # the current narrow-bandwidth KDE
    lr = np.atleast_1d(lr)
    Z = -0.5 * ((lr[:, None] - tx[None, :]) / bw) ** 2
    m = Z.max(1)
    return np.log(np.exp(Z - m[:, None]).sum(1)) + m - np.log(W) - np.log(bw) - 0.5 * np.log(2 * np.pi)

# 2-component GMM (manual EM in log space) — wide components fill the valley via their tails
def fit_gmm2(x, iters=200):
    mu = np.array([mu_dep, mu_enr]); s = np.array([0.5, 1.0]); pi = np.array([0.6, 0.4])
    for _ in range(iters):
        r = np.stack([pi[k] * np.exp(-0.5 * ((x - mu[k]) / s[k]) ** 2) / s[k] for k in range(2)])
        r /= r.sum(0) + 1e-300
        Nk = r.sum(1) + 1e-9
        pi = Nk / Nk.sum()
        mu = (r * x).sum(1) / Nk
        s = np.sqrt((r * (x - mu[:, None]) ** 2).sum(1) / Nk + 1e-6)
    return pi, mu, s
pi_g, mu_g, s_g = fit_gmm2(tx)
print(f"GMM2: pi={pi_g.round(3)}  mu={mu_g.round(2)}  sigma={s_g.round(2)}")
def gmm2(lr):
    lr = np.atleast_1d(lr)
    return np.log(sum(pi_g[k] * np.exp(-0.5 * ((lr - mu_g[k]) / s_g[k]) ** 2) / (s_g[k] * np.sqrt(2 * np.pi)) for k in range(2)) + 1e-300)

# explicit enriched-fraction mixture: rho_g = m*rho_enr + (1-m)*rho_dep, m ~ Uniform[0,1]
# pushforward P(rho_g) = Uniform[rho_dep, rho_enr] (linear); log-space adds Jacobian +log rho_g.
def mmix_uniform_linear(lr):
    lr = np.atleast_1d(lr); rho = np.exp(lr)
    inside = (rho >= rho_dep) & (rho <= rho_enr)
    lp = np.where(inside, lr - np.log(rho_enr - rho_dep), -50.0)   # +log rho Jacobian; hard floor outside
    return lp

# RECOMMENDED: valley-free mixture in LOG space with SOFT tails (no hard rho_enr cap).
# = a wide, near-flat plateau across [mu_dep, mu_enr] (the "you are some mixture" region) with Gaussian
#   shoulders of width = the true between-node spread (use the enriched cluster's own spread, not bw).
sig_shoulder = max(float(np.std(tx[tx > 0.5*(mu_dep+mu_enr)])), 0.5) if (tx > 0.5*(mu_dep+mu_enr)).sum() > 2 else 1.0
def mmix_softlog(lr, tilt=0.0):
    """flat across [mu_dep,mu_enr], Gaussian shoulders (sigma=sig_shoulder). tilt<0 adds a gentle
    downward slope across the plateau (a weak 'depleted is a priori more common' preference)."""
    lr = np.atleast_1d(lr)
    d = np.where(lr < mu_dep, -0.5*((lr-mu_dep)/sig_shoulder)**2,
        np.where(lr > mu_enr, -0.5*((lr-mu_enr)/sig_shoulder)**2, 0.0))
    d = d + tilt*(lr - mu_dep)                                    # optional gentle tilt on the plateau
    return d - np.log((mu_enr-mu_dep) + sig_shoulder*np.sqrt(2*np.pi))
# U-shaped m prior (Beta(0.5,0.5)) pushforward — favors extremes but KEEPS middle support (fills valley)
def mmix_ushape(lr, a=0.5):
    lr = np.atleast_1d(lr); rho = np.clip(np.exp(lr), rho_dep+1e-9, rho_enr-1e-9)
    m = (rho - rho_dep)/(rho_enr - rho_dep)
    # Beta(a,a) density in m, times |dm/drho_g| (const) times |drho_g/dlog rho| = rho
    logbeta = (a-1)*np.log(m) + (a-1)*np.log1p(-m)
    return np.where((np.exp(lr)>rho_dep)&(np.exp(lr)<rho_enr), logbeta + lr, -50.0)
# flat-in-log mixture: log rho_g ~ Uniform[log rho_dep, log rho_enr] with soft Gaussian tails (bw)
def mmix_flatlog(lr):
    lr = np.atleast_1d(lr)
    lo, hi = mu_dep, mu_enr
    d = np.zeros_like(lr)
    d = np.where(lr < lo, -0.5 * ((lr - lo) / bw) ** 2, d)
    d = np.where(lr > hi, -0.5 * ((lr - hi) / bw) ** 2, d)
    return d - np.log((hi - lo) + bw * np.sqrt(2 * np.pi))         # normalization (approx)

# RECOMMENDED MINIMAL FIX: KDE + a uniform "mixture bridge" across [mu_dep,mu_enr] at weight eps.
# Floors the valley (no collapse) while PRESERVING the KDE's real Gaussian tails outside the mode
# range (so the high-rho_g false-positive suppression is unchanged). eps = the a-priori mixture-node
# mass (one knob; to be derived, not magic).
def kde_bridge(lr, eps=0.05):
    lr = np.atleast_1d(lr)
    uni = np.where((lr >= mu_dep) & (lr <= mu_enr), -np.log(mu_enr - mu_dep), -50.0)  # uniform on [dep,enr]
    return np.logaddexp(np.log1p(-eps) + kde_baseline(lr), np.log(eps) + uni)

PRIORS = {"KDE(base)": kde_baseline, "kde+bridge.01": lambda lr: kde_bridge(lr, 0.01),
          "kde+bridge.05": lambda lr: kde_bridge(lr, 0.05), "kde+bridge.20": lambda lr: kde_bridge(lr, 0.20),
          "softlog": mmix_softlog}

# ---------- objective argmax over f_g, as a function of total density d ----------
fg = np.linspace(1e-3, 1 - 1e-3, 400)
jeff = -np.log1p(-fg)
def argmax_fg(prior, d):
    lr = np.log(fg * d)
    obj = prior(lr) + jeff
    return float(fg[np.argmax(obj)])

print("\nCollapse curve — argmax f_g vs total density d (AMBIG, prior+Jeffreys only):")
ds = [0.02, 0.1, 0.5, 1.0, 2.0, 5.17, 8.0, 12.0, 18.9, 30.0]
hdr = "  d      " + "".join(f"{name[:16]:>18s}" for name in PRIORS)
print(hdr)
for d in ds:
    row = f"  {d:6.2f}  " + "".join(f"{argmax_fg(p, d):18.3f}" for p in PRIORS.values())
    tag = ""
    if d == 5.17: tag = "  <- BND 2170 (want HIGH: mostly gDNA)"
    if d == 18.9: tag = "  <- EXON 1085 (want ~0.9)"
    if d == 0.02: tag = "  <- depleted intron (want ~1.0, pure gDNA)"
    print(row + tag)

print("\nInterpretation:")
print("  * baseline KDE should COLLAPSE (argmax f_g -> ~0) across the valley band (d ~ 0.5 .. 12).")
print("  * a valley-free mixture prior should keep argmax f_g monotone/high (no cliff), so BND 2170 does")
print("    not collapse and never emits the pathologic RNA message.")

# sanity: where is the baseline KDE valley + how deep, at the two key nodes
print(f"\n(softlog shoulder sigma={sig_shoulder:.2f})")

# ---------- OVER-CALL RISK: are there real mid-density AMBIG nodes that are truly RNA (low true_gf)? ----------
# If common, a valley-free prior + Jeffreys would over-call gDNA on them. Measure from the actual data.
kind = np.asarray(z["chain_kind"]); fp = np.asarray(z["free_pos"], bool); fn = np.asarray(z["free_neg"], bool)
mg = np.asarray(z["mass_global"], float); eg = np.asarray(z["eff_global"], float)
dtot = mg/np.maximum(eg, 1e-9)
ambig = fp & fn
# region truth join (per-region true_gf via reg_node)
reg_node = np.asarray(z["reg_node"]); true_g = np.asarray(z["reg_true_g"]); reg_raw = np.asarray(z["reg_raw"])
true_gf_by_node = np.full(kind.shape[0], np.nan)
ok = reg_node >= 0
true_gf_by_node[reg_node[ok]] = (true_g[ok]/np.maximum(reg_raw[ok], 1e-9))
valley = (dtot > 0.3) & (dtot < 12) & ambig & np.isfinite(true_gf_by_node)
print(f"\nOVER-CALL RISK check — AMBIG REGION nodes in the collapse band (0.3<d<12), n={int(valley.sum())}:")
if valley.sum():
    tg = true_gf_by_node[valley]
    print(f"  true gDNA fraction: min={tg.min():.2f} p25={np.percentile(tg,25):.2f} med={np.median(tg):.2f} p75={np.percentile(tg,75):.2f} max={tg.max():.2f}")
    print(f"  fraction truly-mostly-RNA (true_gf<0.3): {(tg<0.3).mean():.2f}   truly-mostly-gDNA (true_gf>0.7): {(tg>0.7).mean():.2f}")
    print("  => if mostly-gDNA dominates, a valley-free prior (which leans f_g up via Jeffreys) is a net WIN here;")
    print("     if mostly-RNA is common, the over-call risk is real and messages/strand must arbitrate.")
