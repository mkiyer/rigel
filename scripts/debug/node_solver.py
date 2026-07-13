"""General node solver prototype — the boundary/region unification (theory sandbox, numpy only).

ONE solver for regions AND boundaries. The difference is DATA, not a branch:
  - regions have no spliced counts;   boundaries carry motif-stranded spliced (mature RNA);
  - eff-lengths encode contained vs crossing geometry.

Solves the UNSPLICED count split {RNA+, RNA-, gDNA} (sum = N_unspliced) from:
  (1) structural constraints (free_pos/free_neg — which RNA strands can exist here);
  (2) strand balance (informative ONLY when strand-specific; precision -> 0 at kappa=0.5);
  (3) the SPLICED anchor (boundaries): pure mature RNA of known motif strand, held OUT of the
      solve variable but folded INTO belief + precision, and used to gate the nascent strand;
  (4) OPTIONAL gDNA prior (OFF at init);
  (5) OPTIONAL incoming messages (OFF at init).

Native currency = COUNTS. density(count, E) is a derived helper. Spliced counts are STATIC
(never modified). Message-free self-solve = solve_node(gdna_prior=None, msg=None).

VALUE vs PRECISION (the crux):
  - VALUE  (f_g): at INIT the unspliced defaults to gDNA (nascent~0 — the safe prior-free
    assumption, revisable by later passes). STRAND refines it in-solve (stranded peels nascent);
    unstranded it stays at gDNA. ns_prec is the init nascent~0 default strength, NOT a tuned prior.
  - PRECISION: the honest COMPOSITION precision 1/(1/N_u + 1/N_s + 2phi) from the two independent
    count pools (unspliced gDNA-candidate + spliced RNA), plus the strand's marginal info. This is
    the 'solved' state — non-zero even unstranded, and it is the threshold a later message must beat
    to revise the value ('until proven otherwise', S8).
"""
import numpy as np
from dataclasses import dataclass

PHI = 0.02          # overdispersion floor on log-density (honest precision cap 1/PHI)


# ----------------------------------------------------------------------------- helpers
def density(count, eff_len):
    return count / max(eff_len, 1e-9)


def _logsumexp(a):
    m = np.max(a)
    return m + np.log(np.sum(np.exp(a - m)))


# --------------------------------------------------------------- the intrinsic signals
def strand_loglik(n_pos, n_neg, f_g, f_np, f_nn, kappa, od_g, od_r):
    """Three-component gDNA/RNA strand Beta-Binomial (Normal-moment form, matches simplex.py).
    The COUNT enters ONLY as Fisher information (count-zero-info principle)."""
    N = n_pos + n_neg
    if N <= 0:
        return np.zeros_like(f_g)
    p = 0.5 * f_g + kappa * f_np + (1.0 - kappa) * f_nn
    mean = N * p
    rscale = kappa * (1.0 - kappa)
    var = (mean * (1.0 - p)
           + (N * f_g) ** 2 * 0.25 * od_g
           + (N * f_np) ** 2 * rscale * od_r
           + (N * f_nn) ** 2 * rscale * od_r)
    var = np.maximum(var, 1e-9)
    return -0.5 * (n_pos - mean) ** 2 / var - 0.5 * np.log(var)


def nascent_sparse_prior(f_np, f_nn, ns_prec):
    """WEAK structural prior: nascent RNA is a small fraction of the unspliced pool (nascent-sparse).
    Sets only the VALUE (breaks the unstranded f_g degeneracy toward gDNA); it is deliberately weak
    so it does NOT drive the precision. ns_prec=0 -> off; large -> assert nascent=0 (nrna_none, exact).
    Pulls the nascent fraction (f_np+f_nn) toward 0 in a mild quadratic."""
    if ns_prec <= 0:
        return np.zeros_like(f_np)
    return -0.5 * ns_prec * (f_np + f_nn) ** 2


def composition_precision(N_gdna_pool, N_rna_pool, phi):
    """The HONEST 'solved' precision (theory sec.3): the gDNA/RNA composition log-odds combines two
    INDEPENDENT count pools, so Var(lambda) = 1/N_gdna + 1/N_rna + 2*phi. Needs BOTH a gDNA-candidate
    count (unspliced) AND an independent RNA count (spliced) — that is why a spliced boundary is
    'solved' and a spliced-free unstranded region is not. Returns 0 if either pool is empty."""
    if N_gdna_pool <= 0 or N_rna_pool <= 0:
        return 0.0
    return 1.0 / (1.0 / N_gdna_pool + 1.0 / N_rna_pool + 2.0 * phi)


@dataclass
class NodeSolution:
    # native currency: COUNTS (unspliced split sums to N_unspliced; spliced is static)
    n_gdna: float
    n_rna_pos: float          # nascent (unspliced) RNA+
    n_rna_neg: float          # nascent (unspliced) RNA-
    n_spliced_pos: float      # mature (spliced) RNA+  — static passthrough
    n_spliced_neg: float
    # precision state = 1/Var(log f_c). prec_gdna = the HONEST emitted precision.
    prec_gdna: float
    prec_rna_pos: float
    prec_rna_neg: float
    # diagnostics
    f_g: float
    prec_composition: float   # honest spliced-anchor composition precision (1/N_u+1/N_s+2phi)
    strand_marginal: float    # strand's OWN marginal info above the reference floor (-> 0 at kappa=0.5)


def solve_node(n_u_pos, n_u_neg, n_s_pos=0.0, n_s_neg=0.0,
               free_pos=True, free_neg=True, motif_strand=None,
               E_g=100.0, E_rna=200.0, E_spl=100.0,
               kappa=0.5, od_g=0.0, od_r=0.0, phi=PHI,
               ns_prec=3.0,
               gdna_prior=None, msg=None, n_grid=161):
    """gdna_prior: None (OFF) or (mode_logdensity, prec) on log gDNA density.
       msg: None (OFF) or dict with 'gdna'=(mode_logdensity,prec) / 'rna_pos' / 'rna_neg'.
       motif_strand: 'pos'/'neg'/None — a boundary's spliced motif; gates the nascent strand
                     (nascent shares the transcript strand of the mature RNA)."""
    N_u = n_u_pos + n_u_neg
    N_s = n_s_pos + n_s_neg

    # nascent strand gate: on a splice junction, nascent RNA shares the mature motif strand.
    nas_pos = free_pos and (motif_strand in (None, 'pos'))
    nas_neg = free_neg and (motif_strand in (None, 'neg'))

    fg = np.linspace(1e-3, 1 - 1e-3, n_grid)
    tau = np.linspace(-1.0, 1.0, n_grid)
    FG, T = np.meshgrid(fg, tau, indexing='ij')
    f_act = 1.0 - FG
    if nas_pos and nas_neg:
        f_np = f_act * (1.0 + T) / 2.0
        f_nn = f_act * (1.0 - T) / 2.0
    elif nas_pos:
        f_np, f_nn = f_act, np.zeros_like(f_act)
    elif nas_neg:
        f_np, f_nn = np.zeros_like(f_act), f_act
    else:                                   # no RNA admissible -> all gDNA
        f_np = f_nn = np.zeros_like(f_act)

    psi = np.zeros_like(FG)
    psi = psi + strand_loglik(n_u_pos, n_u_neg, FG, f_np, f_nn, kappa, od_g, od_r)
    psi = psi + nascent_sparse_prior(f_np, f_nn, ns_prec)

    # (4) OPTIONAL gDNA prior on the gDNA log-DENSITY (rho_g = f_g*N_u/E_g)
    if gdna_prior is not None:
        mode, prec = gdna_prior
        log_rho_g = np.log(np.maximum(FG * N_u / E_g, 1e-9))
        psi = psi - 0.5 * prec * (log_rho_g - mode) ** 2

    # (5) OPTIONAL incoming messages on per-component log-DENSITY
    if msg is not None:
        if 'gdna' in msg and msg['gdna'] is not None:
            mode, prec = msg['gdna']
            log_rho_g = np.log(np.maximum(FG * N_u / E_g, 1e-9))
            psi = psi - 0.5 * prec * (log_rho_g - mode) ** 2
        for key, fc in (('rna_pos', f_np), ('rna_neg', f_nn)):
            if key in msg and msg[key] is not None:
                mode, prec = msg[key]
                log_rho = np.log(np.maximum(fc * N_u / E_rna, 1e-9))
                psi = psi - 0.5 * prec * (log_rho - mode) ** 2

    # reference measure (neutral log-odds Jacobian, matches production)
    psi = psi + np.log(FG) + np.log(1.0 - FG)

    post = np.exp(psi - _logsumexp(psi))
    post_fg = post.sum(axis=1)              # marginal over tau
    post_fg = post_fg / post_fg.sum()

    f_g_hat = float(np.sum(post_fg * fg))
    # per-component fractions (posterior means over the full grid)
    f_np_hat = float(np.sum(post * f_np))
    f_nn_hat = float(np.sum(post * f_nn))

    def _var_log(vals_1d, w):
        lv = np.log(np.maximum(vals_1d, 1e-9))
        m = np.sum(w * lv)
        return max(float(np.sum(w * lv * lv) - m * m), 1e-9)

    # --- honest EMITTED precision = strand marginal info + spliced composition, capped at 1/phi ---
    # (1) the strand's OWN marginal information above the reference-measure floor.
    def _grid_prec(psi_terms):
        p = np.exp(psi_terms - _logsumexp(psi_terms))
        pf = p.sum(axis=1); pf /= pf.sum()
        return 1.0 / (_var_log(fg, pf) + phi), pf
    jac = np.log(FG) + np.log(1.0 - FG)
    prec_jac_only, _ = _grid_prec(jac)                                   # reference floor (no strand)
    prec_strand_grid, _ = _grid_prec(
        strand_loglik(n_u_pos, n_u_neg, FG, f_np, f_nn, kappa, od_g, od_r) + jac)
    strand_marginal = max(prec_strand_grid - prec_jac_only, 0.0)          # -> 0 at kappa=0.5
    # (2) the spliced composition precision — the strand-free 'solved' handle.
    prec_composition = composition_precision(N_u, N_s, phi)
    # (3) honest total (independent info adds), capped by the overdispersion ceiling.
    prec_emit = min(strand_marginal + prec_composition, 1.0 / phi)

    prec_rp = prec_emit if nas_pos else 0.0
    prec_rn = prec_emit if nas_neg else 0.0

    return NodeSolution(
        n_gdna=f_g_hat * N_u, n_rna_pos=f_np_hat * N_u, n_rna_neg=f_nn_hat * N_u,
        n_spliced_pos=n_s_pos, n_spliced_neg=n_s_neg,
        prec_gdna=prec_emit, prec_rna_pos=prec_rp, prec_rna_neg=prec_rn,
        f_g=f_g_hat, prec_composition=prec_composition, strand_marginal=strand_marginal)


# ============================================================================ scenarios
def show(tag, s):
    N_u = s.n_gdna + s.n_rna_pos + s.n_rna_neg
    rho_g = density(s.n_gdna, 100.0)
    print(f"{tag}")
    print(f"    unspliced split: gDNA={s.n_gdna:6.1f}  RNA+={s.n_rna_pos:6.1f}  RNA-={s.n_rna_neg:6.1f}"
          f"   (N_u={N_u:.0f})   spliced: +={s.n_spliced_pos:.0f} -={s.n_spliced_neg:.0f}")
    print(f"    f_g={s.f_g:.3f}   PREC_EMIT={s.prec_gdna:6.2f}  (= strand_marginal {s.strand_marginal:5.2f}"
          f" + composition {s.prec_composition:5.2f})   -> emit rho_gdna={rho_g:.3f}")


print("=" * 96)
print("A general node solver — canonical cases. Native currency = counts; spliced static; prior/msg OFF unless noted.")
print("PREC_EMIT = honest strand-marginal + spliced-composition precision (capped 1/phi=50).")
print("=" * 96)

print("\n[S1] UNSTRANDED nrna_none boundary (kappa=0.5): unspliced 300 sym, spliced 100 (motif +). TRUE f_g=1, rho_g=3.0.")
print("     THE CRUX: the composition PRECISION (18.75) is honest, but the f_g VALUE needs the nascent attribution.")
print("     Sweeping ns_prec (nascent-sparse strength) shows value depends on it; precision (composition) does not:")
for nsp in (0.0, 3.0, 30.0, 300.0):
    s = solve_node(150, 150, n_s_pos=100, motif_strand='pos', kappa=0.5, ns_prec=nsp)
    print(f"     ns_prec={nsp:6.1f}:  f_g={s.f_g:.3f}  emit rho_g={density(s.n_gdna,100):.3f}"
          f"   composition={s.prec_composition:.2f}  strand_marginal={s.strand_marginal:.2f}")
print("     -> value reaches the true f_g~1 (rho_g~3.0) only under a STRONG nascent-sparse assumption (exact for nrna_none).")

print("\n[S2] STRANDED nrna_present boundary (kappa=0.9, ns_prec=0): unspliced (220,80) + excess, spliced 100 (+).")
print("     EXPECT: SYNERGY — strand splits unspliced into gDNA + nascent+ (NOT 100% gDNA); strand adds precision.")
show("     ", solve_node(220, 80, n_s_pos=100, motif_strand='pos', kappa=0.9, od_g=0.05, od_r=0.05, ns_prec=0.0))

print("\n[S3] UNSTRANDED region, NO spliced (kappa=0.5): unspliced 300 sym, spliced 0.")
print("     EXPECT: PREC_EMIT ~ 0 — no independent RNA count -> not self-solved; needs prior/messages.")
show("     ", solve_node(150, 150, n_s_pos=0, motif_strand=None, kappa=0.5))

print("\n[S4] Same region + gDNA PRIOR ON (enriched mode, log rho_g ~ log 3.0, prec 8):")
print("     EXPECT: prior resolves the value -> f_g high. Toggle demonstrates prior inclusion.")
show("     ", solve_node(150, 150, n_s_pos=0, motif_strand=None, kappa=0.5, gdna_prior=(np.log(3.0), 8.0)))

print("\n[S5] gdna_none boundary self-gating (kappa=0.5): unspliced 20 sym, spliced 400 (+).")
print("     EXPECT: emitted gDNA DENSITY low (small N_u) — does NOT hallucinate gDNA even though f_g(unspliced) high.")
show("     ", solve_node(10, 10, n_s_pos=400, motif_strand='pos', kappa=0.5))

print("\n[S6] ss=0.5 PRECISION diagnostic across kappa (boundary S1 counts):")
print("     EXPECT: strand_marginal -> ~0 at kappa=0.5, rising with kappa; composition (spliced) keeps the node solved.")
for k in (0.5, 0.6, 0.75, 0.9, 0.99):
    s = solve_node(150, 150, n_s_pos=100, motif_strand='pos', kappa=k, od_g=0.05, od_r=0.05)
    print(f"     kappa={k:.2f}:  strand_marginal={s.strand_marginal:6.2f}   composition={s.prec_composition:6.2f}"
          f"   PREC_EMIT={s.prec_gdna:6.2f}   f_g={s.f_g:.3f}")

print("\n[S7] init vs sweep TOGGLE on the SAME node (S1 counts) — self-defense: honest precision LETS prior/msg adjust.")
base = dict(n_u_pos=150, n_u_neg=150, n_s_pos=100, motif_strand='pos', kappa=0.5)
s0 = solve_node(**base)                                                    # init: both OFF
s1 = solve_node(**base, gdna_prior=(np.log(3.0), 8.0))                     # + prior (prec 8 < 18.8)
s2 = solve_node(**base, msg={'gdna': (np.log(2.0), 40.0)})                 # + strong disagreeing msg (prec 40)
print(f"     init (prior OFF, msg OFF)     : f_g={s0.f_g:.3f}  emit rho_g={density(s0.n_gdna,100):.3f}  prec={s0.prec_gdna:.2f}")
print(f"     + gDNA prior (rho 3.0, p 8)   : f_g={s1.f_g:.3f}  emit rho_g={density(s1.n_gdna,100):.3f}  (barely moves: weaker than node)")
print(f"     + strong msg (rho 2.0, p 40)  : f_g={s2.f_g:.3f}  emit rho_g={density(s2.n_gdna,100):.3f}  (moves: msg stronger than node)")


# ---------------------------------------------------------------- init -> revision
def _logit(p): return np.log(p / (1.0 - p))
def _sig(x): return 1.0 / (1.0 + np.exp(-x))

print("\n[S8] nascent~0 is an INIT default, REVISABLE by later passes (user's principle).")
print("     Init belief (unstranded boundary): value f_g~1 (nascent~0 default), precision = composition 18.75.")
print("     A later pass sends an RNA-evidence message (f_g lower); it revises ONLY if it out-precises the node:")
lam_init, prec_init = _logit(0.98), 18.75            # nascent~0 default @ honest composition precision
for lam_p, pr, tag in [(_logit(0.50), 8.0,  "weak   RNA msg (p 8  < 18.75)"),
                       (_logit(0.50), 40.0, "strong RNA msg (p 40 > 18.75)")]:
    lam_post = (prec_init * lam_init + pr * lam_p) / (prec_init + pr)
    print(f"     + {tag}:  f_g  {_sig(lam_init):.3f} -> {_sig(lam_post):.3f}   (prec {prec_init:.1f} -> {prec_init+pr:.1f})")
print("     -> the assumption holds against weak evidence, yields to strong evidence. Exactly 'until proven otherwise'.")
