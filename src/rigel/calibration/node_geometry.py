"""rigel.calibration.node_geometry — per-slot chain geometry, beliefs, densities, statics, and init.

The lower layer beneath the belief-propagation solver (`bp_solver`): everything that describes *what a
chain slot is* before any message passing. Pure functions of the accumulator substrate + the node
geometry + the junction axis + the FL pmfs — no sweep state, no global prior.

⭐ **EVERY FACE CONCEPT IS GONE (S5.e).** The predecessor built **18 per-face arrays** — a
``_left``/``_right`` pair for the unspliced mass, its integer flux, three effective lengths and four
spliced channels — because a boundary's two sides lay in differently-sized flanks and therefore had
**different divisors**. A contiguous edge is a 0-bp line: one set of numbers, seen identically from both
sides. Dissolved with the pairs:

* the **junction-strand routing** and the **exon-bit flank gating**, which existed to *guess* which flank
  a spliced deposit belonged to — the old accumulator attributed a splice to the node's edge rather than
  to the junction's own coordinate (`CARRY_FORWARD.md` §3 trap 6). The v8 index states ``(src, dst,
  strand)`` explicitly, so the guess is replaced by the chain's own adjacency;
* the ``_continues`` / ``_eff_spl_face`` far-boundary machinery, which chose between two per-face spliced
  divisors. S5.c left ONE crossing formula, so there is nothing to choose;
* ``mass`` versus ``n``. The old accumulator split one fragment's mass across objects, so the density
  numerator was fractional and a separate integer flux was carried for the Poisson variance. The new one
  deposits ``+1`` on every object the fragment touched, so ``count`` is both and ``Var(log rho) = 1/n`` is
  honest against it.

Contents:

* `NodeGeometry` / `build_node_geometry` — the per-slot static geometry: the unspliced count being
  deconvolved, its two per-component divisors, and the mature (junction) flux with its own.
* `NodeBelief` — the per-slot pie ``(f_pos, f_neg, f_g)``; the per-component message densities
  ``rho = f·M/E`` are computed inline in ``bp_solver.node_sweep``.
* `NodeStatics` / `build_node_statics` — the static per-slot solver inputs (per-strand counts, masks).
* `init_beliefs` — the signature-binary G1/G2/G3 initial belief.
* `_node_region_type` — per-slot coarse node type (intergenic/intron/exon), shared with the prior
  subsystem and `gdna_density_prior`.

Layering: imports only `node_chain`, `signature`, `effective_length`, `simplex_logodds`, `strand_deconv`
(all lower layers) — never `bp_solver` or `gdna_density_prior`, so it sits cleanly below both.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..types import Strand
from .effective_length import UNBOUNDED_REACH, contained_eff_length, crossing_eff_length
from .node_chain import EDGE, NODE, NodeChain
from .signature import (
    TS_AMBIG,
    TS_NEG,
    TS_POS,
    coarse_type_array,
    mrna_active_strands,
    nrna_active_strands,
)
from .simplex_logodds import _solve_nodes_logodds_all

__all__ = [
    "NodeGeometry",
    "build_node_geometry",
    "node_global_geometry",
    "node_total_density",
    "NodeBelief",
    "NodeStatics",
    "build_node_statics",
    "init_beliefs",
]

_EPS = 1.0e-9


def _rate(numerator: np.ndarray, divisor: np.ndarray) -> np.ndarray:
    """``numerator / divisor``, and **0 where the divisor is 0** — never a floored division.

    ⛔ `CARRY_FORWARD.md` §3 trap 23 and `effective_length`'s own contract: an object with no opportunity
    for a component must emit nothing, at zero precision. The predecessor floored every divisor to
    ``_EPS`` instead, which is what produced densities of ~1e9 on the **12.4 %** of fine-partition nodes
    where the contained effective length collapses to exactly 0 — turning "no data" into a confident
    wrong answer, and seeding false gDNA into the neighbouring exons.
    """
    d = np.asarray(divisor, np.float64)
    live = d > 0.0
    return np.where(live, np.asarray(numerator, np.float64) / np.where(live, d, 1.0), 0.0)


@dataclass(frozen=True, slots=True)
class NodeGeometry:
    """Per-chain-slot static geometry (length ``n_slots``). **ONE set of numbers per slot.**

    ⭐ **THREE POPULATIONS, NAMED FOR THE ACCUMULATOR'S OWN THREE BANKS.** ``unspliced`` / ``spliced`` /
    ``junction`` is what ``edge_unspliced``, ``edge_spliced`` and ``sj_count`` are called in the
    executable specification, and a consumer that renames a bank on the way in is how one quantity comes
    to have two names. At one line, of the molecules that touched it::

        unspliced_count   crossed it CONTIGUOUSLY, spliced nowhere       a gDNA/RNA MIXTURE
        spliced_count     crossed it CONTIGUOUSLY, spliced elsewhere     certified RNA
        junction_count    never crossed it -- it JUMPED from here        certified RNA

    ⚠ **The last two are different molecules and routinely differ by two orders of magnitude at the same
    line.** At a donor seam the junction flux is the gene's whole mature output while the spliced
    crossing is the handful of molecules that read through without splicing; at a cut inside an exon
    there is no junction at all, so ``junction_count`` is 0 while ``spliced_count`` carries everything.
    ⛔ The word *mature* fits both and therefore distinguishes neither, which is why it is not used.

    ⭐ **TWO STRAND AXES, AND THEY ARE NOT THE SAME AXIS.** ``unspliced_count`` and ``spliced_count``
    are keyed by **genome** strand — where the read aligned, the accumulator's one storage convention.
    ``junction_count``/``eff_junction`` are keyed by **transcript** strand, derived from each junction's
    own annotated strand. That derivation is precisely why the accumulator stores no sense/antisense
    column; putting the two axes in one array under one name is the two-conventions-in-one-schema defect
    the redesign removes.
    """

    n_slots: int

    #: float64[n_slots, 2] — ``edge_unspliced`` at an EDGE, ``node_contained`` at a NODE. By GENOME
    #: strand. **The mixture being deconvolved**, and the only population that is one.
    #: ⭐ It is BOTH the density numerator and the Poisson ``n``: the accumulator deposits ``+1`` on
    #: every object the fragment touched, so there is no fractional mass to carry separately.
    unspliced_count: np.ndarray
    #: float64[n_slots] — the gDNA divisor. Contained placements at a NODE, crossing placements at an EDGE.
    eff_gdna: np.ndarray
    #: float64[n_slots] — the RNA divisor, the same frames on the RNA length pmf.
    #: ⚠ Both are **0 where there is no opportunity**, never floored — see :func:`_rate`.
    eff_rna: np.ndarray
    #: float64[n_slots, 2] — ``edge_spliced``: molecules that crossed this line CONTIGUOUSLY having
    #: spliced somewhere else in the same molecule. By GENOME strand. **Certified RNA — gDNA cannot
    #: splice** — which is what makes it a floor on the RNA inside this line's own population.
    #: 0 on every NODE slot, structurally: ``node_contained`` is credited only when the fragment used no
    #: junction at all (`_accumulator_reference.deposit`).
    #: ⚠ Its divisor is :attr:`eff_rna` — a contiguous crossing is a contiguous crossing whatever the
    #: molecule did elsewhere — so it needs no effective length of its own.
    spliced_count: np.ndarray
    #: float64[n_slots, 2] — ``sj_count``, gathered onto the two lines each junction leaves and enters,
    #: by TRANSCRIPT strand. The flux that **changes template** here: what the graft hands an exon, and
    #: what the peel measures the continuing share against.
    junction_count: np.ndarray
    #: float64[n_slots, 2] — the SUMMED junction divisor, same keying. Several junctions on one line are
    #: several estimates of one rate, so the pooled statement is ``Σcount / ΣE`` — the ratio of sums, never
    #: the mean of ratios (`CARRY_FORWARD.md` §2, ``ρ_bg = Σg/ΣE``).
    eff_junction: np.ndarray


def build_node_geometry(
    chain: NodeChain,
    substrate,
    region_arrays,
    junctions,
    gdna_fl_pmf: np.ndarray,
    rna_fl_pmf: np.ndarray,
) -> NodeGeometry:
    """Assemble the per-slot geometry from the substrate's five populations onto the chain.

    **The divisors, and A7** (`S5_DESIGN_LOG.md` §1 A7, ruled 2026-07-30):

    ================  ==========================================================================
    NODE, both        ``contained_eff_length(node_len, pmf)`` — no reach argument exists
    EDGE, gDNA        ``UNBOUNDED_REACH`` both sides ⇒ ``mu_g − 1``. gDNA's template is the
                      chromosome, so ``taper_g = 1``: settled by physics, not by the ruling
    EDGE, RNA         ⭐ ``UNBOUNDED_REACH`` **by the ruling** ⇒ ``mu_r − 1``. An unspliced crossing
                      is a MIXTURE whose RNA half alone is bounded, so reach there is per
                      COMPONENT; turning it on is S5.g, where it can be A/B'd against S5.f's
                      first baseline. ⚠ Carries a known 11.0 % genome-wide gDNA over-call
    junction, RNA     ⭐ the **real exonic** per-strand reach. A junction is used only by a molecule
                      that spliced across it, and it is a brand-new population with no predecessor
                      divisor, so wiring it regresses nothing
    ================  ==========================================================================

    **Where a junction attaches.** Its donor is the line to the RIGHT of ``src_node`` and its acceptor
    the line to the LEFT of ``dst_node``; molecules leave the template at the first and arrive at the
    second, so both lines genuinely saw the flux. ⭐ Both are read off the chain's own adjacency
    (``right``/``left``) rather than re-derived from per-reference offsets — which is also what makes
    cross-reference leakage impossible, since a reference terminal links to ``-1``.
    """
    kind = np.asarray(chain.kind)
    obj = np.asarray(chain.obj_idx, dtype=np.int64)
    is_node = kind == NODE
    is_edge = kind == EDGE
    n = chain.n_slots

    # ── the two CONTIGUOUS populations: the mixture, and the certified-RNA floor beside it ───────
    unspliced_count = np.zeros((n, 2), dtype=np.float64)
    unspliced_count[is_node] = np.asarray(substrate.node_contained.count, np.float64)[obj[is_node]]
    unspliced_count[is_edge] = np.asarray(substrate.edge_unspliced.count, np.float64)[obj[is_edge]]
    # ⚠ NODE slots stay 0, and that is structural rather than a shortcut: a contained fragment used no
    # junction, so a node's contained population cannot hold a spliced molecule.
    spliced_count = np.zeros((n, 2), dtype=np.float64)
    spliced_count[is_edge] = np.asarray(substrate.edge_spliced.count, np.float64)[obj[is_edge]]

    # ── the two per-component divisors ───────────────────────────────────────────────────────────
    node_len = np.asarray(region_arrays.region_size_bp, dtype=np.float64)
    n_nodes = node_len.shape[0]
    node_idx = np.clip(obj, 0, max(n_nodes - 1, 0))
    unbounded = np.full(1, UNBOUNDED_REACH)

    def divisor(pmf: np.ndarray) -> np.ndarray:
        contained = contained_eff_length(node_len, pmf) if n_nodes else np.zeros(0)
        crossing = float(crossing_eff_length(pmf, unbounded, unbounded)[0])
        return np.where(
            is_node, contained[node_idx] if n_nodes else 0.0, np.where(is_edge, crossing, 0.0)
        )

    eff_gdna = divisor(gdna_fl_pmf)
    eff_rna = divisor(rna_fl_pmf)

    # ── the JUMPING population: a junction edge is a FACTOR on the lines it leaves and enters ───
    junction_count = np.zeros((n, 2), dtype=np.float64)
    eff_junction = np.zeros((n, 2), dtype=np.float64)
    if junctions.n_junctions:
        slot_of_node = np.zeros(int(chain.n_nodes_total), dtype=np.int64)
        slot_of_node[obj[is_node]] = np.flatnonzero(is_node)
        donor = np.asarray(chain.right)[slot_of_node[np.asarray(junctions.src_node, np.int64)]]
        acceptor = np.asarray(chain.left)[slot_of_node[np.asarray(junctions.dst_node, np.int64)]]
        if np.any(donor < 0) or np.any(acceptor < 0):
            raise ValueError(
                "a junction edge attaches to a reference terminal, which has no line beside it. "
                "Both of an intron's endpoints are interior interfaces of the same reference by "
                "construction (splice_graph I5), so this is a junction axis addressing a different "
                "partition than the payload was scanned on."
            )
        # ⭐ The flux is the junction's own count summed over the GENOME-strand columns, filed under the
        # junction's TRANSCRIPT strand. That join is the whole of "sense is derived, never stored".
        flux = np.asarray(substrate.junction.count, np.float64).sum(axis=1)
        eff = crossing_eff_length(rna_fl_pmf, junctions.reach_lo, junctions.reach_hi)
        column = np.where(np.asarray(junctions.strand) == np.int8(Strand.POS), 0, 1)
        for line in (donor, acceptor):
            np.add.at(junction_count, (line, column), flux)
            np.add.at(eff_junction, (line, column), eff)

    return NodeGeometry(
        n_slots=int(n),
        unspliced_count=unspliced_count,
        eff_gdna=eff_gdna,
        eff_rna=eff_rna,
        spliced_count=spliced_count,
        junction_count=junction_count,
        eff_junction=eff_junction,
    )


def node_global_geometry(geometry: NodeGeometry):
    """Per-slot 'global' gDNA support ``(mass, eff)`` — the basis the enrichment NPMLE (`DensityNPMLE`) is
    fit on and projected onto, shared by :func:`bp_solver.node_sweep` and ``calibrate`` so the fit and the
    projection use one definition.

    ⭐ **It no longer sums two faces.** The predecessor returned ``mass_l + mass_r`` over ``E_l + E_r`` at
    a boundary, with a long note about a ``½`` here that was silently cancelling a ``½`` missing from the
    per-face length — which is why that frame read the correct ρ while every per-face MESSAGE read ρ/2.
    With one set of numbers there is nothing to sum and nothing to cancel.
    """
    return np.asarray(geometry.unspliced_count, np.float64).sum(axis=1), np.asarray(
        geometry.eff_gdna, np.float64
    )


def node_total_density(geometry: NodeGeometry, f_g):
    """The LAZY, composition-aware total density: the SUM of component densities, each in its OWN FL
    frame, from the current belief ``f_g``::

        ρ_unspliced = f_g · (M/E_g)  +  (1−f_g) · (M/E_r)        gDNA-FL for gDNA, RNA-FL for RNA
        ρ_junction  = Σ_s  junction_count_s / eff_junction_s      per TRANSCRIPT strand, edges only

    Returns ``(rho_unspliced, rho_with_junction)`` per slot. This is NEVER a pure-gDNA precompute —
    ``f_g`` is the best current composition; gDNA-FL alone (``f_g = 1``) is only the fallback where
    composition is genuinely unknown.

    ⚠ **``spliced_count`` is deliberately NOT in either total.** It is a contiguous crossing, so it does
    belong in the level in principle — but the predecessor's ``mass_spliced`` entered only the strand
    solve, and folding it into ρ_tot here would be a modelling change smuggled into a rename. It is
    S5.a2's question ("how the new channels enter the solve"), recorded in `S5_DESIGN_LOG.md`.

    ⭐ **The per-face acceptor test dissolves with the faces.** The predecessor returned a triple (node,
    left face, right face) because only the acceptor FACE carried a junction; a line either carries
    junction flux or it does not, so there is one answer per slot and ``junction_count == 0`` says so.
    """
    mass, eff_g = node_global_geometry(geometry)
    fg = np.clip(np.asarray(f_g, dtype=np.float64), 0.0, 1.0)
    rho_unspl = mass * (
        _rate(fg, eff_g) + _rate(1.0 - fg, np.asarray(geometry.eff_rna, np.float64))
    )
    junction = np.asarray(geometry.junction_count, np.float64)
    eff_j = np.asarray(geometry.eff_junction, np.float64)
    rho_junction = _rate(junction[:, 0], eff_j[:, 0]) + _rate(junction[:, 1], eff_j[:, 1])
    return rho_unspl, rho_unspl + rho_junction


@dataclass(frozen=True, slots=True)
class NodeBelief:
    """Per-node solved state on the chain: the composition pie `(f_pos, f_neg, f_g)` over the node's UNSPLICED
    mass + its per-component posterior variance in LOG-FRACTION space, `(var_pos, var_neg, var_gdna)` =
    **`Var(log f_c)`, NOT `Var(f_c)`**. All length ``n_nodes``.

    ⚠ **The variances are log-space** — grid moments of `log f_c` over the λ lattice
    (`simplex_logodds._solve_nodes_logodds`), matching the log-density message currency. They are therefore
    **not bounded by ¼** and routinely exceed it; a consumer needing the LINEAR `Var(f_c)` must convert
    (delta method `Var(f_c) ≈ f_c²·Var(log f_c)`, as `bp_solver.node_sweep` does for `composition_logvar`).

    The variance is the **precision state** (`docs/CARRY_FORWARD.md`): `Var(log f_c)=0` ⇒
    locked/certain (e.g. a forbidden strand), `=∞` ⇒ no information (unsolved). It feeds the honest message
    send — a source's outgoing precision is degraded from its own `Var_own` by the
    communication noise, so an unsure node speaks quietly (Phase 2). The composition is stored as a FRACTION
    (the face-invariant quantity — a boundary has two faces but one composition); density `ρ=f·M_face/E_face`
    is the message currency (computed inline in `bp_solver.node_sweep`), mass `m=f·M_face` (`NodeDeconv`) the
    output."""

    f_pos: np.ndarray
    f_neg: np.ndarray
    f_g: np.ndarray
    var_pos: np.ndarray
    var_neg: np.ndarray
    var_gdna: np.ndarray


# ---------------------------------------------------------------------------
# Initialization — the signature-binary G1/G2/G3 belief on the chain.
# ---------------------------------------------------------------------------
#
# A strand axis is hard-LOCKED (a forbidden strand, an intergenic sink) by the per-node ``allow_pos`` /
# ``allow_neg`` forbid mask in the solve. The init ALSO sets the per-component precision state ``var(f_c)``
# (`docs/CARRY_FORWARD.md`): ``0`` = locked/certain (a forbidden strand, an intergenic gDNA sink), ``inf`` =
# no information (an admissible-but-unsolved axis — it will listen to messages, and emits none until solved).
# A solved single-strand (G2) node takes the strand-solve posterior variance.


def _type_belief(free_pos, free_neg, deconv, mass_unspl):
    """Build the per-node composition ``(f_pos, f_neg, f_g)`` for ONE node type (regions OR boundaries) from its
    signature-binary classification + its strand-only solve.

    ``free_pos``/``free_neg`` are the per-node booleans for whether each strand's RNA axis is admissible (a
    region's own ±transcript bits; a boundary's ±strand CONTINUITY across the seam). ``deconv`` is the
    strand-only :class:`NodeDeconv` (no global, no imputation). The signature-binary default is all-gDNA
    ``{0,0,1}`` (`ARCHITECTURE §3`). The class overrides:

    * **G1** (neither strand free — intergenic region / no-RNA-crossing boundary): a LOCKED gDNA sink — keep
      ``{0,0,1}``.
    * **G2** (exactly one strand free, with data): the STRAND DECONVOLUTION alone resolves the pie (a
      single-strand node is 1-D: ``f_active = 1 − f_g``).
    * **G3** (both strands free — AMBIG): unresolvable by strand → keep the ``{0,0,1}`` default at MAX
      (``inf``) variance; the sweep resolves it from neighbour messages + the global prior.

    Returns the six per-node arrays ``(f_pos, f_neg, f_g, var_pos, var_neg, var_gdna)`` — the composition + the
    precision state (`docs/CARRY_FORWARD.md`): ``var=0`` locked, ``inf`` no-information, else the strand-solve
    posterior variance. The variances are ``Var(log f_c)`` (log-space, unbounded above), never ``Var(f_c)``.
    """
    n = free_pos.shape[0]
    f_pos = np.zeros(n)
    f_neg = np.zeros(n)
    f_g = np.ones(n)  # all-gDNA default (count plays NO role; ARCHITECTURE §3)
    # precision state: gDNA unsolved (inf); a strand axis is locked (0) iff forbidden, else unsolved (inf).
    var_g = np.full(n, np.inf)
    var_p = np.where(free_pos, np.inf, 0.0)
    var_n = np.where(free_neg, np.inf, 0.0)

    g1 = ~free_pos & ~free_neg
    g2 = free_pos ^ free_neg
    g2_active = g2 & (np.asarray(mass_unspl, dtype=np.float64) > 0.0)

    # G2-active: take the strand-only solve (median f_g, mean f±, and the posterior variances). G1 sinks + G3
    # AMBIG keep the {0,0,1} default at MAX variance.
    fgv = np.asarray(deconv.gdna_frac_var, dtype=np.float64)
    fpv = np.asarray(deconv.rna_pos_frac_var, dtype=np.float64)
    fnv = np.asarray(deconv.rna_neg_frac_var, dtype=np.float64)
    f_g[g2_active] = np.asarray(deconv.gdna_frac, dtype=np.float64)[g2_active]
    f_pos[g2_active] = np.asarray(deconv.rna_pos_frac, dtype=np.float64)[g2_active]
    f_neg[g2_active] = np.asarray(deconv.rna_neg_frac, dtype=np.float64)[g2_active]
    var_g[g2_active] = fgv[g2_active]
    var_p[g2_active & free_pos] = fpv[g2_active & free_pos]
    var_n[g2_active & free_neg] = fnv[g2_active & free_neg]

    # G1 sink: lock the gDNA axis (the fractions are already the {0,0,1} default).
    var_g[g1] = 0.0
    return f_pos, f_neg, f_g, var_p, var_n, var_g


@dataclass(frozen=True, slots=True)
class NodeStatics:
    """Per-slot STATIC structural masks (length ``n_slots``). The sweep mutates only the dynamic
    :class:`NodeBelief`; these never change.

    ⭐ **NO COUNT LIVES HERE — this class is structure only.** ``u_pos``/``u_neg``/``mass_unspliced``/
    ``mass_spliced`` used to be carried alongside :class:`NodeGeometry`'s ``mass_*``/``n_unspl_*``
    because the old accumulator's mass was fractional while the strand likelihood needed an integer.
    They are the same number now, so **all three populations sit together on the geometry** — where
    their difference is visible — and a consumer slices them there. One quantity, one place
    (`CARRY_FORWARD.md` §3 trap 27).

    ``free_pos``/``free_neg`` are the nascent-RNA-active axes (a node's own ±transcript bits; an edge's
    ±continuity — the RNA-crossing gate); ``mrna_active_pos``/``mrna_active_neg`` are the mature-RNA-active
    axes (a node's ±exon bits; an edge's ±contiguous exon) that select the per-node solver prior.

    ``edge_flags`` carries the splice graph's 8 structural bits (``TSS_s``/``TES_s``/``DONOR_s``/
    ``ACCEPTOR_s``) at each EDGE slot, ``0`` on NODE slots.
    ⭐ **Raw bits, not pre-derived predicates.** Every consumer wants a different combination of them, and
    P1G_SCOPE's own specified predicate was measured to be nearly the COMPLEMENT of what it was meant to
    replace (plan F10). Compose with :func:`~rigel.calibration.splice_graph.is_terminus` /
    :func:`~rigel.calibration.splice_graph.is_splice_site`. It is ``0`` when no graph was supplied.
    """

    n_slots: int
    free_pos: np.ndarray  # bool — nascent-RNA-active (transcript continuity); the RNA-crossing gate
    free_neg: np.ndarray  # bool
    mrna_active_pos: (
        np.ndarray
    )  # bool — mature-RNA-active (contiguous exon); selects the node prior
    mrna_active_neg: np.ndarray  # bool
    edge_flags: np.ndarray  # uint16 — graph structural bits; 0 on NODE slots


def build_node_statics(
    chain: NodeChain,
    region_arrays,
    edge_flags: np.ndarray | None = None,
) -> NodeStatics:
    """Gather the structural masks onto the chain, in ONE slot-keyed pass.

    ⚠ **It no longer takes a substrate**: with the counts moved onto :class:`NodeGeometry` this reads
    nothing but the signatures and the chain's own adjacency, so an unused parameter would only invite a
    caller to think the two are coupled.

    ⭐ **The region/boundary twin pair collapses.** ``_region_strand_stats`` and
    ``_boundary_strand_stats`` computed the same predicates against two different keyings, and the
    boundary one carried a ``max(left, right)`` over its two sides — a de-duplication of the straddle
    that only existed because the old accumulator deposited a crossing fragment on both sides. An edge
    has one count, so the ``max`` goes with the sides.

    ⭐ **And the flank lookup is now the chain's own adjacency.** ``BoundarySubstrate`` carried explicit
    ``left_region``/``right_region`` arrays with ``-1`` at a reference terminal. Edge endpoints are
    implicit — edge ``i`` lies between node ``i`` and node ``i+1`` — and **an edge always has a node on
    both sides**, so ``chain.left``/``chain.right`` answer it and the ``-1`` branch has no cases.

    The allow mask is the transcript-structure CONTINUITY gate: a strand-``s`` unspliced crossing is
    nascent RNA only where strand ``s`` is present on BOTH flanks. This blocks RNA at a TSS/TES
    (intergenic↔exon → neither strand continuous → a gDNA sink) and at a mixed exon↔AMBIG seam.
    ``mrna_active_s`` is the tighter **mature**-crossing gate (contiguous exon on both flanks).

    ``edge_flags`` is the per-contiguous-edge ``uint16[E]`` from
    :func:`~rigel.calibration.splice_graph.build_edge_flags_array`; ``None`` leaves the field zero,
    which every current consumer treats as "no structural information".
    """
    kind = np.asarray(chain.kind)
    obj = np.asarray(chain.obj_idx, dtype=np.int64)
    is_node = kind == NODE
    is_edge = kind == EDGE
    n = int(chain.n_slots)
    flags = _check_edge_flags(edge_flags, int(chain.n_edges_total))

    sig = np.asarray(region_arrays.signature).astype(np.int64)
    n_nodes = sig.shape[0]
    node_idx = np.clip(obj, 0, max(n_nodes - 1, 0))

    # the signature at each NODE slot, then read through the chain's adjacency at each EDGE slot
    slot_sig = np.where(is_node, sig[node_idx] if n_nodes else 0, 0)
    left = np.clip(np.asarray(chain.left), 0, max(n - 1, 0))
    right = np.clip(np.asarray(chain.right), 0, max(n - 1, 0))
    sig_l = np.where(is_edge, slot_sig[left], 0)
    sig_r = np.where(is_edge, slot_sig[right], 0)

    ts = np.where(is_node, np.asarray(region_arrays.strand_class)[node_idx] if n_nodes else 0, -1)
    nrp_l, nrn_l = nrna_active_strands(sig_l)
    nrp_r, nrn_r = nrna_active_strands(sig_r)
    mrp_l, mrn_l = mrna_active_strands(sig_l)
    mrp_r, mrn_r = mrna_active_strands(sig_r)
    mr_self_p, mr_self_n = mrna_active_strands(slot_sig)

    free_pos = np.where(is_node, (ts == TS_POS) | (ts == TS_AMBIG), nrp_l & nrp_r)
    free_neg = np.where(is_node, (ts == TS_NEG) | (ts == TS_AMBIG), nrn_l & nrn_r)

    return NodeStatics(
        n_slots=n,
        # No per-node spliced FLOOR: spliced (mature) handling is OWNED by the message system (the
        # edge→exon MEASUREMENT source + the exon→edge absorption in `_scan`). A node-local floor would
        # double-count it AND inflate an edge's UNSPLICED f_pos with mature → phantom nascent into
        # introns (matrix-confirmed: removing it is ≥ keeping it in every κ × capture × ±gDNA regime).
        free_pos=free_pos,
        free_neg=free_neg,
        mrna_active_pos=np.where(is_node, mr_self_p, mrp_l & mrp_r),
        mrna_active_neg=np.where(is_node, mr_self_n, mrn_l & mrn_r),
        edge_flags=np.where(is_edge, flags[np.clip(obj, 0, max(flags.shape[0] - 1, 0))], 0).astype(
            np.uint16
        ),
    )


def _check_edge_flags(edge_flags, n_edges: int) -> np.ndarray:
    """Validate the per-contiguous-edge flags against the chain, BEFORE anything else is computed.

    A mis-sized array would shift every flag by one seam — a defect invisible in aggregate and
    undetectable by a bit-identity gate while nothing reads the flags. Refuse it at the door.
    """
    if edge_flags is None:
        return np.zeros(max(n_edges, 1), dtype=np.uint16)
    flags = np.asarray(edge_flags, dtype=np.uint16)
    if flags.shape != (n_edges,):
        raise ValueError(
            f"edge_flags has shape {flags.shape}; expected ({n_edges},), one per contiguous edge. "
            f"Build it with splice_graph.build_edge_flags_array(index) against the SAME index the "
            f"payload was scanned on. ⚠ There are no terminal slots: a reference with k nodes owns "
            f"k-1 lines, not k+1."
        )
    return flags if flags.size else np.zeros(1, dtype=np.uint16)


def init_beliefs(
    chain: NodeChain,
    geometry: NodeGeometry,
    statics: NodeStatics,
    *,
    rna_sense_frac: float,
    gdna_strand_overdispersion: float = 0.0,
    rna_strand_overdispersion: float = 0.0,
    n_grid: int,
    n_grid_ss: int | None = None,
    logodds_window: float = 10.0,
) -> NodeBelief:
    """The signature-binary G1/G2/G3 initial :class:`NodeBelief` on the unified chain.

    All slots are strand-solved by the log-density 1-D/2-D log-odds solver (:mod:`simplex_logodds`,
    ``O(m·K)``; the bare strand likelihood + the Jeffreys reference at single-strand nodes; NO global prior,
    NO imputation — those enter the sweep, P2/P3). The signature-binary class overrides (:func:`_type_belief`)
    then set the G1/G2/G3 belief. Single-strand introns resolve to ``f_g≈0`` from the BB tilt alone (the
    zero-gDNA gate); intergenic / TSS sinks lock at ``{0,0,1}``; AMBIG nodes hold ``{0,0,1}`` at MAX variance
    for the sweep."""
    st = statics
    count = np.asarray(geometry.unspliced_count, np.float64)
    deconv = _solve_nodes_logodds_all(
        count[:, 0],
        count[:, 1],
        st.free_pos,
        st.free_neg,
        count.sum(axis=1),
        # ⚠ the strand solve's certified-RNA floor is strand-agnostic, so the two GENOME-strand columns
        # are summed here rather than stored pre-summed — the geometry keeps the axis it was deposited on.
        np.asarray(geometry.spliced_count, np.float64).sum(axis=1),
        kappa=float(rna_sense_frac),
        od_g=gdna_strand_overdispersion,
        od_r=rna_strand_overdispersion,
        n_grid=n_grid,
        n_grid_ss=n_grid_ss,
        L=logodds_window,
    )
    f_pos, f_neg, f_g, var_p, var_n, var_g = _type_belief(
        st.free_pos, st.free_neg, deconv, count.sum(axis=1)
    )
    return NodeBelief(
        f_pos=f_pos, f_neg=f_neg, f_g=f_g, var_pos=var_p, var_neg=var_n, var_gdna=var_g
    )


# ---------------------------------------------------------------------------
# Node-type helper (the sweep itself lives in bp_solver.node_sweep).
# ---------------------------------------------------------------------------


def _node_region_type(chain, region_arrays):
    """Per-slot coarse node type at NODE slots (0=intergenic, 1=intron, 2=exon; exon>intron), −1 at EDGE
    slots; plus the per-node type array. Single source of truth: :func:`signature.coarse_type_array`."""
    kind = np.asarray(chain.kind)
    idx = np.asarray(chain.obj_idx, dtype=np.int64)
    rtype = coarse_type_array(np.asarray(region_arrays.signature)).astype(np.int64)  # per node
    ri_ = np.clip(idx, 0, max(rtype.shape[0] - 1, 0))
    return np.where(kind == NODE, rtype[ri_], -1), rtype
