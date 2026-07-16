"""A standalone belief-propagation SANDBOX for the calibration initialization + message theory.

Pure math/stats — NO biology, NO rigel imports. A node has an OBSERVED total density D and a hidden
COMPOSITION over components {RNA+, RNA-, gDNA} that sum to 1. We infer the composition. The point of the
sandbox is to develop the initialization + first-message theory cleanly:

  1) the honest DENSITY precision (count / FL-length, Poisson with overdispersion) — the precision derivation;
  2) node SELF-SOLVE (structure / strand tilt / boundary spliced-vs-unspliced) → (composition, precision);
  3) a first belief-propagation pass with NO gDNA hyperprior — messages carry COMPOSITION (log-odds), not
     absolute density, so a low-density boundary correctly imputes the composition of a high-density region.

Everything is in COMPOSITION log-odds. For a single-strand '+' node the composition is 1-D:
    λ = log(f_g / f_+),   f_g = σ(λ),   f_+ = 1 − σ(λ),   gDNA density = f_g · D.
(The 2-D {RNA+,RNA-,gDNA} AMBIG case adds a second log-odds; the 1-D case is the flagship single-strand exon.)

Run:  python scripts/debug/bp_sandbox.py
"""
from __future__ import annotations
import math
from dataclasses import dataclass

import numpy as np

# ------------------------------------------------------------------------------------------------
# 1. THE HONEST DENSITY-PRECISION MODEL
# ------------------------------------------------------------------------------------------------
# A component's density is estimated as   ρ = n / L_eff ,  where n = fragment count and L_eff is the
# FL-based effective length (for a boundary CROSSING, L_eff = mean fragment length μ_FL; for a contained
# REGION, L_eff = E[max(0, L−ℓ)]). We need the precision of log ρ.
#
# Leading order (Poisson). Under a Poisson molecular field of rate ρ (per bp), the number of fragments
# crossing a point is n ~ Poisson(ρ·μ_FL) — longer fragments cross more, exactly ∝ length, so μ_FL (the
# MEAN FL) is the correct normalizer and ρ̂ = n/μ is unbiased. Then Var(log ρ̂) = Var(n)/E[n]² = 1/n.
#   ⇒ leading log-precision = n (the COUNT). The FL SPREAD is ancillary here: given a global FL model, the
#     per-fragment lengths do not change the estimate of ρ (they are ancillary to n), so "1 count from a
#     100bp vs a 300bp fragment" is a NAIVE 1/ℓ estimate; the correct n/μ estimate has variance 1/n only.
#     The FL spread instead shapes L_eff (esp. for short regions) and feeds the overdispersion below.
#
# Overdispersion (Poisson is too generous). Real counts are over-dispersed vs Poisson (GC/mappability,
# clustering, FL-model error, local density heterogeneity). Model a relative-variance floor φ (a
# negative-binomial-style dispersion): Var(log ρ̂) = 1/n + φ. Then
#     precision(log ρ) = 1 / (1/n + φ) = n / (1 + n·φ)
# which SATURATES at 1/φ — a confident (large-n) measurement cannot exceed the overdispersion ceiling.
def log_density_precision(n: float, phi: float) -> float:
    """Honest precision of log(density) from a Poisson-with-overdispersion count. n=0 ⇒ 0 (no measurement)."""
    if n <= 0:
        return 0.0
    return 1.0 / (1.0 / n + phi)


def density(n: float, L_eff: float) -> float:
    return n / max(L_eff, 1e-9)


SIG = lambda x: 1.0 / (1.0 + math.exp(-x))  # noqa: E731  logistic


# ------------------------------------------------------------------------------------------------
# 2. THE NODE + its SELF-SOLVE (composition log-odds λ = log(f_g/f_+) and its precision)
# ------------------------------------------------------------------------------------------------
INF = 1e12  # "infinite" precision (locked / unmovable)


@dataclass
class Node:
    name: str
    D: float                    # OBSERVED total density (high precision; a direct count/L_eff)
    free_pos: bool = True       # is the + RNA strand admissible here?
    free_neg: bool = True       # is the − RNA strand admissible here?
    lam: float = 0.0            # composition log-odds log(f_g / f_+)  (1-D single-strand model)
    prec: float = 0.0           # precision of λ  (0 ⇒ zero-information, will move)
    locked: bool = False        # structural lock (intergenic / off-strand sink): blocks messages

    def f_g(self) -> float:
        return SIG(self.lam)

    def rho_gdna(self) -> float:
        return self.f_g() * self.D


def init_default(name, D) -> Node:
    """Default init: 100% gDNA at ZERO precision. The value is irrelevant (prec 0 ⇒ it will move)."""
    return Node(name, D, lam=+8.0, prec=0.0)  # λ large ⇒ f_g≈1, but prec 0 so it does not matter


def init_intergenic(name, D) -> Node:
    """Structure init — intergenic / TSS / TES: pure gDNA, INFINITE precision. RNA sink, gDNA source,
    unmovable, and BLOCKS messages (nothing can change it, and it emits only its locked gDNA)."""
    return Node(name, D, free_pos=False, free_neg=False, lam=+INF, prec=INF, locked=True)


def selfsolve_strand(node: Node, n_pos: float, n_neg: float, kappa: float, phi: float) -> None:
    """Strand-tilt self-solve for a single-strand region. The strand model gives a tilt whose PRECISION
    scales with strand-specificity via w(κ) = (2κ−1)²: 0 at unstranded (κ=½) → the tilt carries NO
    information; →1 at fully stranded. (In a real BB model the strand info also scales with count; here
    we keep the count effect in the density precision and let w(κ) gate the strand's contribution.)"""
    w = (2.0 * kappa - 1.0) ** 2
    n = n_pos + n_neg
    # a stranded node's own counts say f_g via the strand imbalance; unstranded ⇒ w→0 ⇒ prec→0.
    node.prec = w * log_density_precision(n, phi)
    # (λ mean from strand omitted in the 1-D toy — the point is the PRECISION vanishes at κ=½.)


def selfsolve_boundary(name, D, n_spl_pos, n_spl_neg, n_uns, mu_rna, mu_gdna, phi,
                       kappa=0.5, free_pos=True, free_neg=True) -> Node:
    """Boundary self-solve (the NEW init). A boundary directly observes:
       * spliced fragments — KNOWN mature RNA, and KNOWN strand (the splice motif GT/AG is single-stranded),
         so they go straight into the correct RNA strand density: ρ_+ = n_spl_pos/μ_RNA, ρ_- = n_spl_neg/μ_RNA.
       * unspliced fragments — gDNA + (sparse) nascent; at init assume gDNA: ρ_g = n_uns/μ_gDNA.
    The composition log-odds (for a '+' single-strand boundary) is
         λ = log(ρ_g / ρ_+) = log(n_uns/n_spl_pos) + log(μ_RNA/μ_gDNA)          [FL means set the offset]
    and its precision comes from BOTH counts (independent log-densities add variances):
         Var(λ) = Var(log ρ_g) + Var(log ρ_+) = (1/n_uns + φ) + (1/n_spl + φ).
    CRUCIAL: this precision is NONZERO even for UNSTRANDED data — it comes from the spliced-vs-unspliced
    STRUCTURE (known RNA vs assumed gDNA), not from strand. Strand data would only ADD precision (deconvolving
    nascent from gDNA inside the unspliced count)."""
    rho_pos = density(n_spl_pos, mu_rna)
    rho_neg = density(n_spl_neg, mu_rna)
    rho_g = density(n_uns, mu_gdna)
    # 1-D '+' composition: gDNA vs RNA+ (assume neg strand off / handled separately)
    eps = 1e-9
    lam = math.log(max(rho_g, eps)) - math.log(max(rho_pos, eps))
    var = (1.0 / max(n_uns, eps) + phi) + (1.0 / max(n_spl_pos, eps) + phi)
    prec = 0.0 if (n_uns <= 0 and n_spl_pos <= 0) else 1.0 / var
    nd = Node(name, D, free_pos=free_pos, free_neg=free_neg, lam=lam, prec=prec)
    nd._diag = dict(rho_pos=rho_pos, rho_neg=rho_neg, rho_g=rho_g, f_g=SIG(lam))  # type: ignore[attr-defined]
    return nd


# ------------------------------------------------------------------------------------------------
# 3. BELIEF PROPAGATION on a chain (composition messages, NO gDNA hyperprior)
# ------------------------------------------------------------------------------------------------
def combine(lam_a, p_a, lam_b, p_b):
    """Precision-weighted combine of two Gaussian beliefs on λ. p=0 contributes nothing (zero info)."""
    p = p_a + p_b
    if p <= 0:
        return lam_a, 0.0
    return (p_a * lam_a + p_b * lam_b) / p, p


def bp_pass(nodes, transfer="composition"):
    """One forward + one backward message pass on a linear chain. Each node sends its CURRENT belief
    (λ, prec) to its neighbour; a locked node blocks (emits its lock, absorbs nothing). Returns the solved
    nodes. `transfer='composition'` sends λ (the fraction, scale-free — correct); `transfer='density'`
    sends the absolute gDNA density (the current-tool behaviour — wrong when densities differ)."""
    lam = np.array([n.lam for n in nodes], float)
    prec = np.array([n.prec for n in nodes], float)
    D = np.array([n.D for n in nodes], float)
    locked = np.array([n.locked for n in nodes], bool)

    def msg(src, dst):
        if locked[src] or locked[dst]:
            return  # intergenic/TSS/TES are BARRIERS: unmovable (dst) + do not impose composition on
            # genic neighbours (src). Their gDNA-SOURCE role is to anchor the depleted-background mode of
            # the gDNA prior (a population estimate), NOT to inject f_g=1 into an adjacent expressed exon.
        if transfer == "composition":
            m_lam, m_prec = lam[src], prec[src]           # scale-free: fraction transfers
        else:  # density transfer (the current tool): impute the SOURCE's gDNA density onto dst
            rho_g_src = SIG(lam[src]) * D[src]
            f_g_implied = min(rho_g_src / max(D[dst], 1e-9), 1 - 1e-9)
            m_lam = math.log(max(f_g_implied, 1e-9) / max(1 - f_g_implied, 1e-9))
            m_prec = prec[src]
        lam[dst], prec[dst] = combine(lam[dst], prec[dst], m_lam, m_prec)

    order = list(range(len(nodes)))
    for i in order[:-1]:      # forward L→R: i sends to i+1
        msg(i, i + 1)
    for i in order[:0:-1]:    # backward R→L: i sends to i-1
        msg(i, i - 1)
    for i, nd in enumerate(nodes):
        nd.lam, nd.prec = float(lam[i]), float(prec[i])
    return nodes


def show(nodes, title):
    print(f"\n--- {title} ---")
    print(f"  {'node':14} {'D':>7} {'prec':>10} {'f_g':>7} {'ρ_gDNA':>8}")
    for n in nodes:
        p = "inf" if n.prec >= INF else f"{n.prec:.2f}"
        print(f"  {n.name:14} {n.D:>7.1f} {p:>10} {n.f_g():>7.3f} {n.rho_gdna():>8.2f}")


# ------------------------------------------------------------------------------------------------
# 4. THE TOY EXAMPLE + scenarios
# ------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    MU_RNA, MU_GDNA, PHI = 200.0, 100.0, 0.02

    print("=" * 90)
    print("TOY: '+' splice-junction boundary B  →  adjacent exon region R (unstranded, no gDNA prior)")
    print(f"  FL: RNA μ={MU_RNA}, gDNA μ={MU_GDNA};  overdispersion φ={PHI}")
    print("  B: 100 spliced(+)  → ρ_RNA+ = 100/200 = 0.5 ;  300 unspliced → ρ_gDNA = 300/100 = 3.0")
    print("  R: 10000 unspliced over eff-len ~ (500 − FL)  → total density ≈ 33  (composition UNKNOWN)")
    print("  TRUTH: R is mostly gDNA (B's composition ⇒ f_g ≈ 0.857) ⇒ R gDNA density ≈ 0.857·33 ≈ 28")

    def fresh():
        B = selfsolve_boundary("B(bnd,+)", D=3.5, n_spl_pos=100, n_spl_neg=0, n_uns=300,
                               mu_rna=MU_RNA, mu_gdna=MU_GDNA, phi=PHI, free_neg=False)
        R = init_default("R(exon,+)", D=33.0); R.free_neg = False
        return [B, R]

    nodes = fresh()
    show(nodes, "INIT (B self-solved from spliced/unspliced; R at zero precision)")
    print(f"    B self-solve: ρ_RNA+={nodes[0]._diag['rho_pos']:.2f} ρ_gDNA={nodes[0]._diag['rho_g']:.2f} "
          f"f_g={nodes[0].f_g():.3f} prec(λ)={nodes[0].prec:.2f}  ← NONZERO precision with NO strand data")

    bp_pass(fresh_nodes := fresh(), transfer="composition")
    show(fresh_nodes, "AFTER BP  (COMPOSITION transfer — correct)")
    bp_pass(fresh_nodes2 := fresh(), transfer="density")
    show(fresh_nodes2, "AFTER BP  (DENSITY transfer — the current tool's behaviour, WRONG)")
    print("  ⇒ composition-transfer gives R f_g≈0.857 (ρ_gDNA≈28, correct); density-transfer imputes B's"
          "\n    absolute ρ_gDNA=3.0 onto R ⇒ f_g≈0.09 (ρ_gDNA≈3) — the undersampling collapse.")

    print("\n" + "=" * 90)
    print("SCENARIO gdna_none: same boundary but 0 unspliced (no gDNA), only spliced RNA")
    Bn = selfsolve_boundary("B(no gDNA)", D=0.5, n_spl_pos=100, n_spl_neg=0, n_uns=0,
                            mu_rna=MU_RNA, mu_gdna=MU_GDNA, phi=PHI, free_neg=False)
    Rn = init_default("R(exon,+)", 33.0); Rn.free_neg = False
    bp_pass([Bn, Rn], transfer="composition"); show([Bn, Rn], "AFTER BP")
    print("  ⇒ SELF-GATING: no unspliced ⇒ ρ_gDNA=0 ⇒ f_g→0 ⇒ R imputed as RNA. No false gDNA, no prior needed.")

    print("\n" + "=" * 90)
    print("SCENARIO stranded: strand data ADDS precision to R's own self-solve (does not need the message)")
    for kappa in (0.5, 0.99):
        Rk = init_default("R", 33.0); Rk.free_neg = False
        selfsolve_strand(Rk, n_pos=9000, n_neg=1000, kappa=kappa, phi=PHI)
        print(f"  κ={kappa}: R strand self-solve precision = {Rk.prec:.3f}  "
              f"(w(κ)=(2κ−1)²={ (2*kappa-1)**2:.3f})  ← 0 at unstranded, high when stranded")

    print("\n" + "=" * 90)
    print("SCENARIO intergenic blocks: [Intergenic(locked gDNA) | R(exon) | B(RNA boundary)]")
    IG = init_intergenic("IG(locked)", 0.05)
    Rb = init_default("R(exon,+)", 33.0); Rb.free_neg = False
    Bb = selfsolve_boundary("B(RNA,+)", 0.5, n_spl_pos=200, n_spl_neg=0, n_uns=2,
                            mu_rna=MU_RNA, mu_gdna=MU_GDNA, phi=PHI, free_neg=False)
    bp_pass([IG, Rb, Bb], transfer="composition"); show([IG, Rb, Bb], "AFTER BP")
    print("  ⇒ IG stays locked (prec inf, unmovable); R takes the RNA boundary's composition (f_g→0)."
          "\n    The intergenic node does not inject its gDNA into R (R's evidence is the RNA boundary).")
