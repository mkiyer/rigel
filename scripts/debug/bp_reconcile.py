"""How do messages reconcile DIFFERENCES IN TOTAL DENSITY? (the message-currency derivation)

Neither pure composition nor pure absolute-density transfer is right:
 * composition transfer breaks when the ACTIVE COMPONENT SET changes (intergenic gDNA-only → intron→exon
   boundary where RNA+ turns on): composition jumps discontinuously, so you cannot copy it across.
 * absolute-density transfer leaves the receiver's OWN observed total density unallocated (a boundary with
   total 10 cannot tell a region with total 100 how to split its 100).

Correct model — PER-COMPONENT density messages + the receiver's observed total as a HARD SUM CONSTRAINT:
 * each component is a DENSITY FIELD with its own continuity: gDNA is active everywhere (a smooth field
   modulo capture enrichment); each RNA strand is a coverage field active only along its transcript, 0 where
   off. Messages carry per-component DENSITY (log ρ_c, honest precision), gated by that component's continuity.
 * the receiver observes its TOTAL density D at high precision. Its components MUST sum to D. The solve finds
     ρ_c ≥ 0 with Σ ρ_c = D   minimizing   Σ_c prec_c · (log ρ_c − m_c)² .
 * this reconciles total-density differences: the messages set the per-component targets; the sum constraint
   scales/allocates them to the receiver's own D; the EXCESS (D − Σ imputed) flows to the LOWEST-precision
   active component. The split therefore interpolates between "preserve the composition (scale to D)" when the
   component precisions are balanced, and "pin the confident component, residual the rest" when one dominates.
 * THE COLLAPSE is exactly: gDNA imputed at LOW precision (depleted/undersampled) + RNA left as the free
   (0-precision) residual ⇒ the real excess gDNA is dumped into RNA.

Run:  python scripts/debug/bp_reconcile.py
"""
from __future__ import annotations
import math
import numpy as np

EPS = 1e-9


def reconcile_2(D, m_g, p_g, m_r, p_r, n=4001):
    """Two active components (gDNA + one RNA strand). Solve ρ_g+ρ_r=D minimizing
    p_g(log ρ_g − m_g)² + p_r(log ρ_r − m_r)² over the split f_g∈(0,1). m_* are log-density targets;
    p_*=0 means that component is active but UNINFORMED (it absorbs the residual)."""
    fg = np.linspace(EPS, 1 - EPS, n)
    rg = fg * D
    rr = (1 - fg) * D
    cost = p_g * (np.log(rg) - m_g) ** 2 + p_r * (np.log(rr) - m_r) ** 2
    i = int(np.argmin(cost))
    return rg[i], rr[i], fg[i]


def show(label, D, rows, rg, rr, fg):
    print(f"  {label}")
    for nm, m, p in rows:
        tgt = "—" if p <= 0 else f"{math.exp(m):.2f}"
        print(f"      msg {nm:5} target ρ={tgt:>6}  prec={p:>5.1f}")
    print(f"      → SPLIT: ρ_gDNA={rg:6.2f}  ρ_RNA={rr:6.2f}  (f_g={fg:.3f})   [total {D:.0f} fully allocated]\n")


print("=" * 92)
print("Q: HOW DO MESSAGES RECONCILE TOTAL-DENSITY DIFFERENCES?")
print("Setup: boundary B (total 10) sends to region R (total 100).")
print("  B composition {RNA+ 90% → ρ_RNA=9, gDNA 10% → ρ_gDNA=1}.  R must allocate its OWN total of 100.")
print("  Messages carry B's per-component densities as log-targets: m_gDNA=log(1), m_RNA=log(9).\n")

D = 100.0
mg, mr = math.log(1.0), math.log(9.0)
# regime 1: gDNA imputed confidently, RNA uninformed  → gDNA pinned near 1, RNA = residual (COLLAPSE shape)
rg, rr, fg = reconcile_2(D, mg, 50.0, mr, 1.0)
show("regime A — gDNA-confident, RNA weak:", D, [("gDNA", mg, 50), ("RNA", mr, 1)], rg, rr, fg)
# regime 2: RNA imputed confidently (e.g. spliced coverage), gDNA uninformed → RNA pinned near 9, gDNA=residual
rg, rr, fg = reconcile_2(D, mg, 1.0, mr, 50.0)
show("regime B — RNA-confident (spliced), gDNA weak:", D, [("gDNA", mg, 1), ("RNA", mr, 50)], rg, rr, fg)
# regime 3: balanced precisions → does NOT preserve composition; the excess flows to the LARGER target.
rg, rr, fg = reconcile_2(D, mg, 10.0, mr, 10.0)
show("regime C — balanced precisions (excess → larger target):", D, [("gDNA", mg, 10), ("RNA", mr, 10)], rg, rr, fg)
print("  ⇒ the TOTAL is always fully allocated (sum constraint). But the SPLIT is NOT composition-preserving:")
print("    matching independent per-component ABSOLUTE densities in log-space PINS the smaller/confident target")
print("    and dumps the excess (D − Σtargets) into the LARGER / weaker-precision component (here RNA, 9>1).")
print("    So the deconvolution reduces to: impute at least ONE component confidently AND CORRECTLY; the other")
print("    is the residual. The gDNA→RNA COLLAPSE = gDNA target imputed too LOW (depleted/undersampled) while")
print("    RNA is the large free residual ⇒ the real excess gDNA is dumped into RNA. Fixing it REQUIRES a")
print("    correct high-precision gDNA-field target (the enriched mode) — total-density reconciliation cannot")
print("    save us if no component is confidently+correctly imputed.")

print("\n" + "=" * 92)
print("ACTIVE-SET CHANGE: intron I (gDNA only) → intron↔exon boundary B (RNA+ turns ON).")
print("  Composition is NOT transferable (I is 100% gDNA; B is not). But the gDNA DENSITY FIELD transfers.")
print("  I: total 5, all gDNA → gDNA density 5.   B: total 12; RNA+ activates, informed by B's OWN spliced.\n")
# gDNA message from I is an ABSOLUTE density (continuous field): m_g = log(5).  RNA+ from B's spliced coverage.
DB = 12.0
mg_field = math.log(5.0)                    # gDNA field density transfers from the intron (continuous)
mr_spliced = math.log(6.0)                  # RNA+ coverage measured at B (its spliced fragments)
rg, rr, fg = reconcile_2(DB, mg_field, 8.0, mr_spliced, 8.0)
show("B solve (gDNA-field + own-spliced RNA, both informed):", DB,
     [("gDNA", mg_field, 8), ("RNA+", mr_spliced, 8)], rg, rr, fg)
print("  ⇒ gDNA transfers as a DENSITY (5, from the intron field) — NOT as I's composition (which was 100%).")
print("    RNA+ is a NEW active component, informed by B's own spliced coverage, not by the intron.")
print("    Composition emerges (≈ {gDNA 42%, RNA+ 58%}) from per-component density fields + the sum-to-D solve.")

print("\n" + "=" * 92)
print("FLAGSHIP RECOVERY: gDNA-DOMINANT boundary (unspliced ≫ spliced) → adjacent ENRICHED exon.")
print("  B: 300 unspliced → gDNA density 3.0 ;  100 spliced(+) → RNA+ 0.5.   R (exon) total density 33.")
print("  gDNA is the LARGER target (3.0 > 0.5), so the excess (33 − 3.5) flows to gDNA — even though B's")
print("  crossing gDNA density (3.0) UNDERSAMPLES R's interior (~30). Direction transfers; magnitude = R's own D.\n")
rg, rr, fg = reconcile_2(33.0, math.log(3.0), 10.0, math.log(0.5), 10.0)
show("R solve:", 33.0, [("gDNA", math.log(3.0), 10), ("RNA+", math.log(0.5), 10)], rg, rr, fg)
print("  ⇒ R recovers to f_g ≈ 0.97 (ρ_gDNA ≈ 32), NOT the collapsed 0.002 — and it is ROBUST to the boundary")
print("    undersampling, because the spliced:unspliced RATIO sets the DIRECTION and R's own total sets the")
print("    magnitude. The collapse only happens when gDNA is imputed from a DEPLETED off-target neighbour")
print("    (tiny target) while RNA is the larger residual — i.e. from the WRONG source, not the own junction.")
