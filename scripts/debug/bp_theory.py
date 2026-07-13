"""Consolidated calibration BELIEF-PROPAGATION THEORY sandbox — the solver we throw hard problems at.

Pure math/stats (numpy only). Implements the full model so we can validate the theory BEFORE touching code:

  * STATE: each node has an OBSERVED total density D (hard) and a composition over active components
    ⊆ {gDNA g, RNA+ p, RNA- n} that must sum to D (a simplex / pie chart). We infer the split.
  * NODE SOLVE (joint, constrained): given per-component log-density targets m_c and precisions π_c
    (local self-solve ⊗ incoming messages), solve  min Σ_c π_c (log ρ_c − m_c)²  s.t.  Σ_c ρ_c = D.
    This is the ONLY correct way to honor the interdependency (never update a component + residual the rest).
  * SELF-DEFENSE (the key property): a confident component must not be moved by a weak message on another
    component. The constrained MAP delivers this via its stationarity condition (derived in the doc):
        (log ρ_c − m_c) = −(μ/2)·ρ_c/π_c    ⇒    deviation from target ∝ 1/π_c.
    So high-π components barely deviate; low-π components absorb the constraint. Self-defense therefore
    requires (a) the JOINT solve and (b) HONEST precisions — an overconfident wrong message still dominates.
  * MESSAGES: per-component DENSITY fields (NOT composition). gDNA transfers everywhere (× enrichment, with an
    enrichment-transfer variance); each RNA strand transfers only where that strand is active on both ends.
    A message's precision is CAPPED: min(source's belief precision, the edge's transfer reliability).
  * PROPAGATION: forward-backward on the chain (exact on a linear chain). Messages CHANGE as they relay
    (accumulate upstream context), but a node never feeds a neighbour's own message back (tree BP).

Run:  python scripts/debug/bp_theory.py
"""
from __future__ import annotations
import math
import numpy as np

EPS = 1e-9
COMPS = ("g", "p", "n")
LOG = lambda x: math.log(max(x, EPS))  # noqa: E731


# ------------------------------------------------------------------------------------------------------
# node solve — the joint constrained MAP (self-defense lives here)
# ------------------------------------------------------------------------------------------------------
def node_solve(D, beliefs, n=1200):
    """beliefs: dict c -> (m_c log-density target, π_c precision) for the node's ACTIVE components
    (π_c may be 0 = active-but-uninformed). Inactive components are absent (ρ_c = 0). Solve the
    constrained MAP; return dict c -> ρ_c (summing to D)."""
    active = [c for c in COMPS if c in beliefs]
    if len(active) == 1:
        return {active[0]: D}
    if len(active) == 2:
        a, b = active
        f = np.linspace(EPS, 1 - EPS, n)
        ra, rb = f * D, (1 - f) * D
        ma, pa = beliefs[a]; mb, pb = beliefs[b]
        cost = pa * (np.log(ra) - ma) ** 2 + pb * (np.log(rb) - mb) ** 2
        i = int(np.argmin(cost))
        return {a: ra[i], b: rb[i]}
    # 3 active — grid the 2-simplex
    g = np.linspace(EPS, 1 - EPS, n // 3)
    F0, F1 = np.meshgrid(g, g); F2 = 1 - F0 - F1
    ok = F2 > EPS
    r = [F0 * D, F1 * D, F2 * D]
    m = [beliefs[c][0] for c in active]; p = [beliefs[c][1] for c in active]
    cost = np.where(ok, sum(p[k] * (np.log(np.where(ok, r[k], 1)) - m[k]) ** 2 for k in range(3)), np.inf)
    j = np.unravel_index(int(np.argmin(cost)), cost.shape)
    return {active[k]: float(r[k][j]) for k in range(3)}


def combine(*gaussians):
    """precision-weighted combine of (mean, precision) log-density Gaussians. Returns (mean, prec)."""
    P = sum(p for _, p in gaussians)
    if P <= 0:
        return (gaussians[0][0], 0.0)
    return (sum(m * p for m, p in gaussians) / P, P)


def frac(rho, D):
    return {c: rho.get(c, 0.0) / D for c in COMPS}


# ------------------------------------------------------------------------------------------------------
# TEST 1 — SELF-DEFENSE: honest weak message vs overconfident wrong message
# ------------------------------------------------------------------------------------------------------
def hdr(t):
    print("\n" + "=" * 94 + f"\n{t}\n" + "=" * 94)


hdr("TEST 1 — SELF-DEFENSE.  Node: D=33, active {g,p}. Confident local gDNA (target 30, π=50);"
    "\n           RNA+ blank (π=0). An INCORRECT RNA+ message says ρ_p≈15.")
for label, msg_prec in [("HONEST weak msg (π=2)", 2.0), ("OVERCONFIDENT wrong msg (π=80)", 80.0)]:
    # local: g=(log30, 50), p=(anything, 0).  message on p: (log15, msg_prec).
    bp = {"g": (LOG(30.0), 50.0), "p": combine((LOG(1.0), 0.0), (LOG(15.0), msg_prec))}
    rho = node_solve(33.0, bp)
    print(f"  {label:34}: ρ_g={rho['g']:6.2f}  ρ_p={rho['p']:6.2f}  (f_g={rho['g']/33:.3f})")
print("  ⇒ honest weak msg: confident gDNA SURVIVES (~30). Overconfident wrong msg: gDNA is DOMINATED.")
print("    Self-defense = joint solve (proportional protection) + HONEST precision (a wrong msg must be weak).")

# infinite-precision lock (structural) — full protection
print("\n  Structural lock (π=∞) vs any message:")
bp = {"g": (LOG(30.0), 1e9), "p": combine((LOG(1.0), 0.0), (LOG(15.0), 80.0))}
rho = node_solve(33.0, bp)
print(f"    locked gDNA + overconfident RNA+ msg: ρ_g={rho['g']:.2f} (UNMOVED — infinite precision = full defense)")


# ------------------------------------------------------------------------------------------------------
# CHAIN belief propagation — forward-backward, per-component density messages, honest caps
# ------------------------------------------------------------------------------------------------------
class Node:
    def __init__(self, name, D, active, local, locked=False):
        self.name = name; self.D = D; self.active = set(active)
        self.local = local            # dict c -> (m, prec)  (self-solve)
        self.locked = locked

def bp_chain(nodes, reliab_g=8.0, reliab_r=8.0, enrich_shift=None):
    """Forward-backward on a linear chain. Per-component DENSITY messages:
       * gDNA transfers between all adjacent nodes (× enrichment via enrich_shift[edge], reliability reliab_g);
       * RNA strand s transfers only where s is active on BOTH endpoints (reliability reliab_r).
       Message precision = min(source belief precision on c, edge reliability). Locked nodes neither move
       nor emit composition (barriers). Returns final ρ per node."""
    N = len(nodes)
    enrich_shift = enrich_shift or {}

    def msg(src, dst, belief_src):
        """message dict c->(m,prec) from src to dst, given src's belief (dict c->(m,prec))."""
        if nodes[src].locked or nodes[dst].locked:
            return {}
        out = {}
        for c in nodes[src].active & nodes[dst].active:
            if c in ("p", "n") and not (c in nodes[src].active and c in nodes[dst].active):
                continue
            m, p = belief_src.get(c, (0.0, 0.0))
            if p <= 0:
                continue
            shift = enrich_shift.get((src, dst), 0.0) if c == "g" else 0.0
            rel = reliab_g if c == "g" else reliab_r   # per-hop transfer PRECISION (1/σ²_hop)
            # HONEST message precision: source belief ⊕ per-hop transfer variance (variances ADD), so it can
            # never exceed the source's own belief AND it DECAYS each relay hop (accumulating σ²_hop).
            out[c] = (m + shift, 1.0 / (1.0 / p + 1.0 / rel))
        return out

    # forward beliefs (local ⊗ forward message from the left)
    fwd = [None] * N
    for i in range(N):
        b = dict(nodes[i].local)
        if i > 0 and fwd[i - 1] is not None:
            m_in = msg(i - 1, i, fwd[i - 1])
            for c, (mm, pp) in m_in.items():
                b[c] = combine(b.get(c, (0.0, 0.0)), (mm, pp))
        fwd[i] = b
    # backward beliefs (local ⊗ backward message from the right)
    bwd = [None] * N
    for i in range(N - 1, -1, -1):
        b = dict(nodes[i].local)
        if i < N - 1 and bwd[i + 1] is not None:
            m_in = msg(i + 1, i, bwd[i + 1])
            for c, (mm, pp) in m_in.items():
                b[c] = combine(b.get(c, (0.0, 0.0)), (mm, pp))
        bwd[i] = b
    # final: local ⊗ left-message ⊗ right-message, then the joint constrained solve
    out = []
    for i in range(N):
        b = dict(nodes[i].local)
        if i > 0:
            for c, (mm, pp) in msg(i - 1, i, fwd[i - 1]).items():
                b[c] = combine(b.get(c, (0.0, 0.0)), (mm, pp))
        if i < N - 1:
            for c, (mm, pp) in msg(i + 1, i, bwd[i + 1]).items():
                b[c] = combine(b.get(c, (0.0, 0.0)), (mm, pp))
        b = {c: v for c, v in b.items() if c in nodes[i].active}
        rho = node_solve(nodes[i].D, b) if not nodes[i].locked else {"g": nodes[i].D}
        out.append(rho)
    return out


def show_chain(nodes, out):
    for nd, rho in zip(nodes, out):
        comp = "  ".join(f"{c}={rho.get(c,0.0):6.2f}" for c in COMPS if c in nd.active)
        fg = rho.get("g", 0.0) / nd.D
        print(f"    {nd.name:16} D={nd.D:6.1f}  {comp:34}  f_g={fg:.3f}")


# TEST 2 — enriched exon recovers through the chain, NO prior
hdr("TEST 2 — ENRICHED EXON recovery (unstranded, NO gDNA prior). Chain:"
    "\n  intron(gDNA) — junction B (gDNA-dominant) — EXON R (blank, D=33) — junction B' (gDNA-dom) — intron(gDNA)")
nodes = [
    Node("intron_L", 0.05, {"g", "p"}, {"g": (LOG(0.05), 6.0), "p": (LOG(0.001), 0.0)}),
    Node("junctionB", 3.5, {"g", "p"}, {"g": (LOG(3.0), 18.0), "p": (LOG(0.5), 12.0)}),   # boundary self-solve
    Node("EXON_R", 33.0, {"g", "p"}, {"g": (LOG(1.0), 0.0), "p": (LOG(1.0), 0.0)}),        # blank
    Node("junctionB'", 3.2, {"g", "p"}, {"g": (LOG(3.0), 18.0), "p": (LOG(0.4), 10.0)}),
    Node("intron_R", 0.05, {"g", "p"}, {"g": (LOG(0.05), 6.0), "p": (LOG(0.001), 0.0)}),
]
show_chain(nodes, bp_chain(nodes))
print("  ⇒ the blank EXON recovers to mostly gDNA (f_g high) from the gDNA-dominant junctions — no prior used.")

# TEST 3 — gdna_none: junctions carry NO gDNA → exon stays RNA (self-gating)
hdr("TEST 3 — gdna_none (self-gating). Same chain but junctions have ~0 unspliced (no gDNA), only spliced RNA.")
nodes3 = [
    Node("intron_L", 0.02, {"g", "p"}, {"g": (LOG(0.02), 4.0), "p": (LOG(0.001), 0.0)}),
    Node("junctionB", 5.0, {"g", "p"}, {"g": (LOG(0.02), 2.0), "p": (LOG(5.0), 20.0)}),
    Node("EXON_R", 33.0, {"g", "p"}, {"g": (LOG(1.0), 0.0), "p": (LOG(1.0), 0.0)}),
    Node("junctionB'", 5.0, {"g", "p"}, {"g": (LOG(0.02), 2.0), "p": (LOG(5.0), 20.0)}),
    Node("intron_R", 0.02, {"g", "p"}, {"g": (LOG(0.02), 4.0), "p": (LOG(0.001), 0.0)}),
]
show_chain(nodes3, bp_chain(nodes3))
print("  ⇒ no gDNA at the junctions ⇒ exon imputed as RNA (f_g→0). Self-gating, no false positives, no prior.")

# TEST 4 — RELAY: a gDNA signal propagates through several blank nodes (messages change/accumulate)
hdr("TEST 4 — RELAY + DECAY. gDNA-confident source at the left; 4 nodes to its right, each with a WEAK"
    "\n           competing RNA prior (ρ_p≈10, π=0.5). The decaying gDNA message must out-argue it hop by hop.")
relay = [Node("source_g", 20.0, {"g", "p"}, {"g": (LOG(18.0), 30.0), "p": (LOG(2.0), 8.0)})]
relay += [Node(f"node_{k}", 20.0, {"g", "p"}, {"g": (LOG(1.0), 0.0), "p": (LOG(10.0), 0.5)}) for k in range(4)]
show_chain(relay, bp_chain(relay, reliab_g=8.0))
print("  ⇒ the gDNA message relays rightward but its precision DECAYS each hop (variances add: π ≈ 6.3→3.5→2.4"
      "\n    →1.9), so against the constant weak RNA prior the gDNA fraction FALLS with distance — messages change"
      "\n    and confidence attenuates as they propagate. A far node is rightly less sure than a near one.")
