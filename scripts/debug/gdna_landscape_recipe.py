def recipe(s):
    import numpy as np
    from scipy.special import gammaln
    import gdna_explore_lib as L
    LN10 = np.log(10.0); GRID = L.GRID; lnat = GRID * LN10; EPS = 1e-12
    # kernel/grid resolution floor = the library's own fit_kde bandwidth (0.15 decades) as a
    # log-rate variance. A fixed reference SCALE on the rho axis, not a per-sample percentile.
    S0 = (0.15 * LN10) ** 2
    mk = L.masks(s)
    # --- substrate: region base + zero-DNA structural anchor, PLUS boundary nodes ---
    region = (mk['base'] & (mk['single'] | mk['gonly'])) | mk['struct_zero']
    boundary = mk['boundary'] & (s['eff'] > 1e-9) & (s['mass'] > EPS) & (mk['single'] | mk['gonly'])
    # --- boundary admission = DERIVED information criterion. tau_prior = per-node log10-density
    # spread of the CLEAN intergenic reference (no genes -> no nascent leak). Belief-tau is
    # bimodal (confident cluster vs an uninformative "give-up" cluster at f0~0.5) with tau_prior
    # in the empty gap; a boundary joins only if its belief pins density TIGHTER than the
    # population background varies. ---
    im = s['is_region'] & (s['ntype'] == 0) & (s['eff'] > 1e-9) & (~mk['ambig'])
    g_i, E_i = s['g_hat'][im], s['eff'][im]
    if im.sum() > 3 and g_i.sum() > 0:
        ln = np.log10(np.maximum(g_i, 1.0)) - np.log10(E_i)
        tau_prior = float(np.clip(np.std(ln), 0.15, 0.7))
    else:
        tau_prior = 0.3
    tau_n = np.sqrt(np.maximum(s['var'], 0.0)) / LN10
    info = np.isfinite(s['var']) & (tau_n < tau_prior)
    sel = region | (boundary & info)
    if not sel.any():
        return None
    g = s['g_hat'][sel]; E = s['eff'][sel]
    v = np.nan_to_num(s['var'][sel], nan=1e9, posinf=1e9); v = np.maximum(v, 0.0)
    # --- per-node RELIABILITY as MASS, DERIVED reference = sum of the two IRREDUCIBLE variance
    # sources of the log gDNA rate: vf = 1/max(g,1) (own Poisson counting floor, self-referential)
    # and S0 (kernel/grid resolution floor). w = ref/(var+ref) = irreducible share vs deconvolution
    # AMBIGUITY. Flooring by S0 keeps confident high-count nodes (var~S0 -> w~0.5) so REAL enriched
    # modes survive, while give-up nodes (var>>ref) collapse to ~0 mass, killing the nascent-leak
    # fabrication. ---
    ref = 1.0 / np.maximum(g, 1.0) + S0
    w = ref / (v + ref)
    w = np.where(mk['struct_zero'][sel], 1.0, w)   # zero-count anchor = trusted zero-DNA
    # --- zero-native additive Poisson landscape: the count sets each node's own width ---
    lam = np.exp(lnat)[None, :] * E[:, None]
    ll = g[:, None] * np.log(np.maximum(lam, EPS)) - lam - gammaln(g[:, None] + 1.0)
    ll -= ll.max(1, keepdims=True)
    pn = np.exp(ll); pn /= np.maximum(pn.sum(1, keepdims=True), EPS)
    d = (w[:, None] * pn).sum(0); ss = float(d.sum())
    return d / ss if ss > 0 and np.isfinite(ss) else None