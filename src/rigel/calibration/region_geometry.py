"""rigel.calibration.region_geometry — per-slot chain geometry, beliefs, densities, statics, and init.

The lower layer beneath the belief-propagation backbone (`sweep`): everything that describes *what a
chain slot is* before any message passing. Pure functions of the accumulator substrate + the region
geometry + the sj axis + the FL pmfs — no sweep state, no global prior.

⭐ **EVERY FACE CONCEPT IS GONE (S5.e).** The predecessor built **18 per-face arrays** — a
``_left``/``_right`` pair for the unspliced mass, its integer flux, three effective lengths and four
spliced channels — because a boundary's two sides lay in differently-sized flanks and therefore had
**different divisors**. A contiguous boundary is a 0-bp boundary: one set of numbers, seen identically from both
sides. Dissolved with the pairs:

* the **sj-strand routing** and the **exon-bit flank gating**, which existed to *guess* which flank
  a spliced deposit belonged to — the old accumulator attributed a splice to the region's boundary rather than
  to the sj's own coordinate. The v8 index states ``(src, dst,
  strand)`` explicitly, so the guess is replaced by the chain's own adjacency;
* the ``_continues`` / ``_eff_spl_face`` far-boundary machinery, which chose between two per-face spliced
  divisors. S5.c left ONE crossing formula, so there is nothing to choose;
* ``mass`` versus ``n``. The old accumulator split one fragment's mass across objects, so the density
  numerator was fractional and a separate integer flux was carried for the Poisson variance. The new one
  deposits ``+1`` on every object the fragment touched, so ``count`` is both and ``Var(log rho) = 1/n`` is
  honest against it.

Contents:

* `RegionGeometry` / `build_region_geometry` — the per-slot static geometry: the unspliced count being
  deconvolved, its two per-component divisors, and the mature (sj) flux with its own.
* `RegionBelief` — the per-slot pie ``(f_pos, f_neg, f_g)``; the per-component message densities
  ``rho = f·M/E`` are computed inline in ``sweep.solve_chain``.
* `RegionStatics` / `build_region_statics` — the static per-slot solver inputs (per-strand counts, masks).
* `init_beliefs` — the signature-binary G1/G2/G3 initial belief.
* `_region_region_type` — per-slot coarse region type (intergenic/intron/exon), a thin re-keying of
  `signature.coarse_type_array`. ⛔ **It has NO caller today — inside this package or outside it** (grep,
  2026-08-17). This line claimed it was *"shared with the prior subsystem and ``gdna_density_prior``"*;
  the sharing is `coarse_type_array`'s, and `gdna_density_prior` is a module DELETED at `9b0f7419` whose
  successor is `landscape`. Retiring the helper is a source change and is left to the owner.

Layering: LAYER 3. It imports DOWN to `region_chain` and `signature` (0), `splice_graph` (1) and
`effective_length` (2), and SIDEWAYS to `simplex_logodds` (3) — never `sweep` (layer 6) or `landscape`
(layer 5), so it sits cleanly below both. ⚠ This line read "(all lower layers)" until 2026-08-17, which
is wrong about `simplex_logodds`: `_layers` puts it in layer 3 beside this module, so that import is
sideways, which its own rule (down or sideways, never up) allows. ⚠ `splice_graph` is imported for the
four TERMINUS FLAG BITS alone, which
:func:`terminus_flank_gain` has to know the meaning of; the module already carried the flags array
through `RegionStatics` without knowing what any bit meant.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..types import Strand
from .effective_length import UNBOUNDED_REACH, contained_eff_length, crossing_eff_length
from .region_chain import BOUNDARY, REGION, RegionChain
from .signature import (
    TS_AMBIG,
    TS_NEG,
    TS_POS,
    coarse_type_array,
    mrna_active_strands,
    nrna_active_strands,
)
from .simplex_logodds import _solve_regions_logodds_all
from .splice_graph import FLAG_TES_NEG, FLAG_TES_POS, FLAG_TSS_NEG, FLAG_TSS_POS

__all__ = [
    "RegionGeometry",
    "build_region_geometry",
    "region_gdna_geometry",
    "region_total_density",
    "RegionBelief",
    "RegionStatics",
    "build_region_statics",
    "init_beliefs",
    "g1_locked",
    "terminus_flank_gain",
]

_EPS = 1.0e-9


def _rate(numerator: np.ndarray, divisor: np.ndarray) -> np.ndarray:
    """``numerator / divisor``, and **0 where the divisor is 0** — never a floored division.

    ⛔ and `effective_length`'s own contract: an object with no opportunity
    for a component must emit nothing, at zero precision. The predecessor floored every divisor to
    ``_EPS`` instead, which is what produced densities of ~1e9 on the **12.4 %** of fine-partition regions
    where the contained effective length collapses to exactly 0 — turning "no data" into a confident
    wrong answer, and seeding false gDNA into the neighbouring exons.
    """
    d = np.asarray(divisor, np.float64)
    live = d > 0.0
    return np.where(live, np.asarray(numerator, np.float64) / np.where(live, d, 1.0), 0.0)


@dataclass(frozen=True, slots=True)
class RegionGeometry:
    """Per-chain-slot static geometry (length ``n_slots``). **ONE set of numbers per slot.**

    ⭐ **THREE POPULATIONS, NAMED FOR THE ACCUMULATOR'S OWN THREE BANKS.** ``unspliced`` / ``spliced`` /
    ``sj`` is what ``boundary_unspliced``, ``boundary_spliced`` and ``sj_count`` are called in the
    executable specification, and a consumer that renames a bank on the way in is how one quantity comes
    to have two names. At one boundary, of the molecules that touched it::

        unspliced_count   crossed it CONTIGUOUSLY, spliced nowhere       a gDNA/RNA MIXTURE
        spliced_count     crossed it CONTIGUOUSLY, spliced elsewhere     certified RNA
        sj_count    never crossed it -- it JUMPED from here        certified RNA

    ⚠ **The last two are different molecules and routinely differ by two orders of magnitude at the same
    boundary.** At a donor boundary the sj flux is the gene's whole mature output while the spliced
    crossing is the handful of molecules that read through without splicing; at a region_bound inside an exon
    there is no sj at all, so ``sj_count`` is 0 while ``spliced_count`` carries everything.
    ⛔ The word *mature* fits both and therefore distinguishes neither, which is why it is not used.

    ⭐ **TWO STRAND AXES, AND THEY ARE NOT THE SAME AXIS.** ``unspliced_count`` and ``spliced_count``
    are keyed by **genome** strand — where the read aligned, the accumulator's one storage convention.
    ``sj_count``/``eff_sj`` are keyed by **transcript** strand, derived from each sj's
    own annotated strand. That derivation is precisely why the accumulator stores no sense/antisense
    column; putting the two axes in one array under one name is the two-conventions-in-one-schema defect
    the redesign removes.
    """

    n_slots: int

    #: float64[n_slots, 2] — ``boundary_unspliced`` at an BOUNDARY, ``region_contained`` at a REGION. By GENOME
    #: strand. **The mixture being deconvolved**, and the only population that is one.
    #: ⭐ It is BOTH the density numerator and the Poisson ``n``: the accumulator deposits ``+1`` on
    #: every object the fragment touched, so there is no fractional mass to carry separately.
    unspliced_count: np.ndarray
    #: float64[n_slots] — the gDNA divisor. Contained placements at a REGION, crossing placements at an BOUNDARY.
    eff_gdna: np.ndarray
    #: float64[n_slots] — the RNA divisor, the same frames on the RNA length pmf.
    #: ⚠ Both are **0 where there is no opportunity**, never floored — see :func:`_rate`.
    eff_rna: np.ndarray
    #: float64[n_slots, 2] — ``boundary_spliced``: molecules that crossed this boundary CONTIGUOUSLY having
    #: spliced somewhere else in the same molecule. By GENOME strand. **Certified RNA — gDNA cannot
    #: splice** — which is what makes it a floor on the RNA inside this boundary's own population.
    #: 0 on every REGION slot, structurally: ``region_contained`` is credited only when the fragment used no
    #: sj at all (`_accumulator_reference.deposit`).
    #: ⚠ Its divisor is :attr:`eff_rna` — a contiguous crossing is a contiguous crossing whatever the
    #: molecule did elsewhere — so it needs no effective length of its own.
    spliced_count: np.ndarray
    #: float64[n_slots, 2] — ``sj_count``, gathered onto the two boundaries each sj leaves and enters,
    #: by TRANSCRIPT strand. The flux that **changes template** here: what the graft hands an exon, and
    #: what the peel measures the continuing share against.
    sj_count: np.ndarray
    #: float64[n_slots, 2] — the SUMMED sj divisor, same keying. Several sj on one boundary are
    #: several estimates of one rate, so the pooled statement is ``Σcount / ΣE`` — the ratio of sums, never
    # the mean of ratios (``ρ_bg = Σg/ΣE``).
    eff_sj: np.ndarray
    #: ⭐⭐ float64[n_slots, 2] — the SAME flux, split by WHICH GENOMIC END of its sj this BOUNDARY is.
    #: ``_lo`` is the flux of sj whose genomic-LOW end is here; ``_hi`` the genomic-HIGH end's.
    #: ⛔ **This split is what makes the reframe's total well defined**, and the reason is that the two
    #: halves belong to DIFFERENT FLANKS of the same BOUNDARY. A molecule that splices at this position has
    #: its body in the exon on ONE side of it — the low side if this is the sj's low end, the high
    #: side if this is its high end — and it never enters the other flank at all. So the total this BOUNDARY
    #: presents to its low neighbour must count ``_lo`` and not ``_hi``, and vice versa
    #: (:func:`region_total_density`). Summing them, as ``sj_count`` above does, is right for the
    #: GRAFT — which is about the whole flux leaving this boundary — and wrong for the reframe.
    #: ⚠ Written in GENOMIC terms, never donor/acceptor: the index's ``FLAG_DONOR_s`` bit marks the
    #: genomic-LOW end of an ``s``-strand intron on BOTH strands, so on ``−`` it sits at the transcript's
    #: biological ACCEPTOR. Naming these ``_lo``/``_hi`` is what stops that being a sign error
    #: (`test_splice_flux_reframe`).
    sj_count_lo: np.ndarray
    sj_count_hi: np.ndarray
    #: float64[n_slots, 2] — the matching divisors, split the same way.
    eff_sj_lo: np.ndarray
    eff_sj_hi: np.ndarray


def build_region_geometry(
    chain: RegionChain,
    substrate,
    region_arrays,
    sj,
    gdna_fl_pmf: np.ndarray,
    rna_fl_pmf: np.ndarray,
    boundary_rna_reach=None,
) -> RegionGeometry:
    """Assemble the per-slot geometry from the substrate's five populations onto the chain.

    **The divisors, and TRAPS: prove-the-substrate** (ruled 2026-07-30):

    ================  ==========================================================================
    REGION, both        ``contained_eff_length(region_len, pmf)`` — no reach argument exists
    BOUNDARY, gDNA        ``UNBOUNDED_REACH`` both sides ⇒ ``mu_g − 1``. gDNA's template is the
                      chromosome, so ``taper_g = 1``: settled by physics, not by the ruling
    BOUNDARY, RNA         ⭐ ``UNBOUNDED_REACH`` **by the ruling** ⇒ ``mu_r − 1``. An unspliced crossing
                      is a MIXTURE whose RNA half alone is bounded, so reach there is per
                      COMPONENT; turning it on is S5.g, where it can be A/B'd against S5.f's
                      first baseline. ⚠ Carries a known 11.0 % genome-wide gDNA over-call
    sj, RNA     ⭐ the **real exonic** per-strand reach. A sj is used only by a molecule
                      that spliced across it, and it is a brand-new population with no predecessor
                      divisor, so wiring it regresses nothing
    ================  ==========================================================================

    ``boundary_rna_reach`` is the **TRAPS: prove-the-substrate switch**: ``None`` (the default) keeps ``UNBOUNDED_REACH`` at
    contiguous boundaries and is byte-identical to the S5.f path; a ``(reach_lo, reach_hi)`` pair per
    contiguous boundary (:func:`~rigel.calibration.splice_graph.build_contiguous_boundary_reach_arrays`) turns
    the taper on. ⚠ It is ONE argument so that an A/B varies one thing and shares every boundary of code.

    **Where a sj attaches.** Its donor is the boundary to the RIGHT of ``src_region`` and its acceptor
    the boundary to the LEFT of ``dst_region``; molecules leave the template at the first and arrive at the
    second, so both boundaries genuinely saw the flux. ⭐ Both are read off the chain's own adjacency
    (``right``/``left``) rather than re-derived from per-reference offsets — which is also what makes
    cross-reference leakage impossible, since a reference terminal links to ``-1``.
    """
    kind = np.asarray(chain.kind)
    obj = np.asarray(chain.obj_idx, dtype=np.int64)
    is_region = kind == REGION
    is_boundary = kind == BOUNDARY
    n = chain.n_slots

    # ── the two CONTIGUOUS populations: the mixture, and the certified-RNA floor beside it ───────
    unspliced_count = np.zeros((n, 2), dtype=np.float64)
    unspliced_count[is_region] = np.asarray(substrate.region_contained.count, np.float64)[obj[is_region]]
    unspliced_count[is_boundary] = np.asarray(substrate.boundary_unspliced.count, np.float64)[obj[is_boundary]]
    # ⚠ REGION slots stay 0, and that is structural rather than a shortcut: a contained fragment used no
    # sj, so a region's contained population cannot hold a spliced molecule.
    spliced_count = np.zeros((n, 2), dtype=np.float64)
    spliced_count[is_boundary] = np.asarray(substrate.boundary_spliced.count, np.float64)[obj[is_boundary]]

    # ── the two per-component divisors ───────────────────────────────────────────────────────────
    region_len = np.asarray(region_arrays.region_size_bp, dtype=np.float64)
    n_regions = region_len.shape[0]
    unbounded = np.full(1, UNBOUNDED_REACH)

    def divisor(pmf: np.ndarray, boundary_reach=None) -> np.ndarray:
        """Per-slot effective length: contained at a REGION, crossing at an BOUNDARY.

        ``boundary_reach`` is ``(reach_lo, reach_hi)`` per contiguous boundary, or ``None`` for
        :data:`UNBOUNDED_REACH` — the TRAPS: prove-the-substrate switch. ``None`` is byte-identical to the pre-S5.g path, so
        the two arms of the A/B differ in ONE argument and share every boundary of code.
        """
        contained = contained_eff_length(region_len, pmf) if n_regions else np.zeros(0)
        n_boundaries = max(int(chain.n_boundaries_total), 1)
        if boundary_reach is None:
            crossing = np.full(n_boundaries, float(crossing_eff_length(pmf, unbounded, unbounded)[0]))
        else:
            crossing = crossing_eff_length(pmf, boundary_reach[0], boundary_reach[1])
        out = np.zeros(n, dtype=np.float64)
        if n_regions:
            out[is_region] = contained[obj[is_region]]
        if is_boundary.any():
            out[is_boundary] = crossing[np.clip(obj[is_boundary], 0, crossing.shape[0] - 1)]
        return out

    # ⭐ gDNA takes NO reach argument, ever: its template is the chromosome, so ``taper_g = 1``. That is
    # physics, not the TRAPS: prove-the-substrate ruling — the ruling is only about the RNA component.
    eff_gdna = divisor(gdna_fl_pmf)
    eff_rna = divisor(rna_fl_pmf, boundary_rna_reach)

    # ── the JUMPING population: a sj boundary is a FACTOR on the boundaries it leaves and enters ───
    sj_count = np.zeros((n, 2), dtype=np.float64)
    eff_sj = np.zeros((n, 2), dtype=np.float64)
    #: the same flux kept apart by which genomic END of its sj the boundary is — see the dataclass.
    jc_lo = np.zeros((n, 2), dtype=np.float64)
    jc_hi = np.zeros((n, 2), dtype=np.float64)
    ej_lo = np.zeros((n, 2), dtype=np.float64)
    ej_hi = np.zeros((n, 2), dtype=np.float64)
    if sj.n_sj:
        slot_of_region = np.zeros(int(chain.n_regions_total), dtype=np.int64)
        slot_of_region[obj[is_region]] = np.flatnonzero(is_region)
        donor = np.asarray(chain.right)[slot_of_region[np.asarray(sj.src_region, np.int64)]]
        acceptor = np.asarray(chain.left)[slot_of_region[np.asarray(sj.dst_region, np.int64)]]
        if np.any(donor < 0) or np.any(acceptor < 0):
            raise ValueError(
                "a sj boundary attaches to a reference terminal, which has no boundary beside it. "
                "Both of an intron's endpoints are interior interfaces of the same reference by "
                "construction (splice_graph I5), so this is a sj axis addressing a different "
                "partition than the payload was scanned on."
            )
        # ⭐ The flux is the sj's own count summed over the GENOME-strand columns, filed under the
        # sj's TRANSCRIPT strand. That join is the whole of "sense is derived, never stored".
        flux = np.asarray(substrate.sj.count, np.float64).sum(axis=1)
        eff = crossing_eff_length(rna_fl_pmf, sj.reach_lo, sj.reach_hi)
        column = np.where(np.asarray(sj.strand) == np.int8(Strand.POS), 0, 1)
        # ⚠ ``donor`` is the boundary at the sj's genomic-LOW end and ``acceptor`` the genomic-HIGH
        # one — for BOTH strands, because ``chain.right``/``chain.left`` are genomic and boundaries run
        # ``src < dst`` (`splice_graph`, DESIGN §2). The names are the index's; the meaning is genomic.
        for boundary, jc, ej in ((donor, jc_lo, ej_lo), (acceptor, jc_hi, ej_hi)):
            np.add.at(sj_count, (boundary, column), flux)
            np.add.at(eff_sj, (boundary, column), eff)
            np.add.at(jc, (boundary, column), flux)
            np.add.at(ej, (boundary, column), eff)

    return RegionGeometry(
        n_slots=int(n),
        unspliced_count=unspliced_count,
        eff_gdna=eff_gdna,
        eff_rna=eff_rna,
        spliced_count=spliced_count,
        sj_count=sj_count,
        eff_sj=eff_sj,
        sj_count_lo=jc_lo,
        sj_count_hi=jc_hi,
        eff_sj_lo=ej_lo,
        eff_sj_hi=ej_hi,
    )


def region_gdna_geometry(geometry: RegionGeometry):
    """Per-slot gDNA support ``(unspliced count, eff_gdna)`` — the basis the enrichment NPMLE
    (`DensityNPMLE`) is fit on and projected onto, shared by :func:`sweep.solve_chain` and ``calibrate``
    so the fit and the projection use one definition.

    ⛔ **It was called ``region_global_geometry`` and bound as ``eff_global`` until 2026-08-10, and the
    name was a lie in the one place it could do most harm.** It returns ``eff_gdna`` and nothing else,
    so ``rho_g = f_g·M/E_g`` — which is exactly what makes ``sum_c rho_c·E_c = M`` hold and therefore
    what makes ``f_g`` a COUNT share rather than a density share. A reader auditing whether the EM prior
    mixes the two had to open two files to learn that "global" meant "gDNA"
    (`TRAPS: two-masks-one-name`).

    ⭐ **It no longer sums two faces.** The predecessor returned ``mass_l + mass_r`` over ``E_l + E_r`` at
    a boundary, with a long note about a ``½`` here that was silently cancelling a ``½`` missing from the
    per-face length — which is why that frame read the correct ρ while every per-face MESSAGE read ρ/2.
    With one set of numbers there is nothing to sum and nothing to cancel.
    """
    return np.asarray(geometry.unspliced_count, np.float64).sum(axis=1), np.asarray(
        geometry.eff_gdna, np.float64
    )


def region_rna_geometry(geometry: RegionGeometry):
    """Per-slot RNA support ``(unspliced count, eff_rna)`` — the exact mirror of
    :func:`region_gdna_geometry`, for the OTHER component.

    ⭐ **The count is the SAME object total ``M``, and only the divisor differs.** ψ splits one
    unspliced population, so ``rho_r = (1-f_g)*M/E_r`` exactly as ``rho_g = f_g*M/E_g`` — which is what
    makes ``sum_c rho_c*E_c = M`` hold and ``f_g`` a COUNT share rather than a density share.
    :func:`region_total_density` already pairs ``1-f_g`` with ``eff_rna`` for the same reason.

    ⚠ Returning the pair rather than just ``eff_rna`` is deliberate: a caller that took the count from
    one helper and the divisor from another could silently mix bases, and the two helpers exist so the
    fit and the projection use ONE definition per component.
    """
    return np.asarray(geometry.unspliced_count, np.float64).sum(axis=1), np.asarray(
        geometry.eff_rna, np.float64
    )


def region_total_density(geometry: RegionGeometry, f_g):
    """⭐⭐⭐ The LAZY, composition-aware total density — **as a PAIR, one per FLANK**::

        ρ_unspliced = f_g · (M/E_g)  +  (1−f_g) · (M/E_r)        gDNA-FL for gDNA, RNA-FL for RNA
        ρ_lo        = ρ_unspliced + Σ_s sj_count_lo_s / eff_sj_lo_s
        ρ_hi        = ρ_unspliced + Σ_s sj_count_hi_s / eff_sj_hi_s

    Returns ``(rho_lo, rho_hi)``: the total to use when this slot is compared against its genomic-LOW
    neighbour, and the one for its genomic-HIGH neighbour. ⚠ **Equal at every REGION** — a REGION stores
    only CONTAINED fragments and a contained fragment used no sj, so both banks are 0 there and the
    pair degenerates to ``ρ_unspliced``. The distinction exists only at an BOUNDARY.

    ⛔⛔ **WHY IT IS A PAIR, AND WHY ONE NUMBER PER SLOT CANNOT BE RIGHT.** The reframe
    ``r = ρ_tot(dst)/ρ_tot(src)`` is a COMPOSITION imputation, so its numerator and denominator must be
    totals over the SAME component set — the intersection of what the two slots can carry. At an BOUNDARY the
    sj flux is the density of molecules that SPLICE at this position, and such a molecule's body
    lies in the exon on exactly ONE side of it: the low side if this BOUNDARY is the sj's genomic-LOW
    end, the high side if it is the genomic-HIGH end. It never enters the other flank. So:

        against the LOW neighbour  → count the flux of sj that START here    (ρ_lo)
        against the HIGH neighbour → count the flux of sj that END here      (ρ_hi)

    Using one sj-inclusive total on both sides — which is what the predecessor's second return
    value did, on every hop, in both directions and both twins — inflates the side facing the INTRON by
    exactly ``ρ_J/ρ_unspliced``: measured **1.28×** and **1.43×** at the two ``intron|exon`` BOUNDARIES of a
    two-exon toy, against a TRUTH ratio the split reproduces to 3 %. ⚠ It inflates each hop of a
    two-hop pair in opposite directions and therefore CANCELS in a compounded ratio, which is why no
    endpoint or aggregate check saw it. Using the unspliced-only total on both sides is the other
    mistake, and it was measured worse: at an BOUNDARY→EXON step the exon genuinely contains the spliced
    population and the term belongs.

    ⭐ This is `EQUATIONS.md` §3.6's two faces, made per-step and per-sj rather than per-object: the
    INTRON face is whichever flank is not the exonic side of the sj, and the EXON face is the other.

    ⛔ Written in GENOMIC terms, never donor/acceptor or TSS/TES. The index's ``FLAG_DONOR_s`` marks the
    genomic-LOW end of an ``s``-strand intron on BOTH strands, so on ``−`` it sits at the transcript's
    biological ACCEPTOR — a predicate phrased biologically flips sign with the strand and a genomic one
    does not. Gated on a ``−``-strand sj specifically (`test_splice_flux_reframe`).

    This is NEVER a pure-gDNA precompute — ``f_g`` is the best current composition; gDNA-FL alone
    (``f_g = 1``) is only the fallback where composition is genuinely unknown.

    ⚠ **``spliced_count`` is deliberately NOT in either total.** It is a contiguous crossing, so it does
    belong in the level in principle — but the predecessor's ``mass_spliced`` entered only the strand
    solve, and folding it into ρ_tot here would be a modelling change smuggled into a rename.
    """
    mass, eff_g = region_gdna_geometry(geometry)
    fg = np.clip(np.asarray(f_g, dtype=np.float64), 0.0, 1.0)
    rho_unspl = mass * (
        _rate(fg, eff_g) + _rate(1.0 - fg, np.asarray(geometry.eff_rna, np.float64))
    )

    def _flux(count, eff):
        c, e = np.asarray(count, np.float64), np.asarray(eff, np.float64)
        return _rate(c[:, 0], e[:, 0]) + _rate(c[:, 1], e[:, 1])

    return (
        rho_unspl + _flux(geometry.sj_count_lo, geometry.eff_sj_lo),
        rho_unspl + _flux(geometry.sj_count_hi, geometry.eff_sj_hi),
    )


@dataclass(frozen=True, slots=True)
class RegionBelief:
    """Per-region solved state on the chain: the composition pie `(f_pos, f_neg, f_g)` over the region's UNSPLICED
    mass + its per-component posterior variance in LOG-FRACTION space, `(var_pos, var_neg, var_gdna)` =
    **`Var(log f_c)`, NOT `Var(f_c)`**. All length ``n_slots``.

    ⚠ **The first axis is the unified region+boundary CHAIN, not the region axis.** This line said
    ``n_regions`` until 2026-08-17, while its own first sentence already said "on the chain" and its last
    one already said a BOUNDARY has one composition: :func:`init_beliefs` builds every array from
    ``geometry.unspliced_count``, which is ``float64[n_slots, 2]``. The same defect
    `region_init.RegionInit` and `landscape.DensityLandscape.logprior` each carried and had corrected.

    ⚠ **The variances are log-space** — grid moments of `log f_c` over the λ lattice
    (`simplex_logodds._solve_regions_logodds`), matching the log-density message currency. They are therefore
    **not bounded by ¼** and routinely exceed it; a consumer needing the LINEAR `Var(f_c)` must convert
    (delta method `Var(f_c) ≈ f_c²·Var(log f_c)`, as `sweep.solve_chain` does for `composition_logvar`).

    The variance is the **precision state**: `Var(log f_c)=0` ⇒
    locked/certain (e.g. a forbidden strand), `=∞` ⇒ no information (unsolved). It feeds the honest message
    send — a source's outgoing precision is degraded from its own `Var_own` by the
    communication noise, so an unsure region speaks quietly (Phase 2). The composition is stored as a FRACTION
    (the face-invariant quantity — a boundary has two faces but one composition); density `ρ=f·M_face/E_face`
    is the message currency (computed inline in `sweep.solve_chain`), mass `m=f·M_face` (`RegionDeconv`) the
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
# A strand axis is hard-LOCKED (a forbidden strand, an intergenic sink) by the per-region ``allow_pos`` /
# ``allow_neg`` forbid mask in the solve. The init ALSO sets the per-component precision state ``var(f_c)``
# ``0`` = locked/certain (a forbidden strand, an intergenic gDNA sink), ``inf`` =
# no information (an admissible-but-unsolved axis — it will listen to messages, and emits none until solved).
# A solved single-strand (TRAPS: one-thing-varied) region takes the strand-solve posterior variance.


def g1_locked(free_pos, free_neg) -> np.ndarray:
    """The **TRAPS: no-magic-numbers** class: neither RNA strand admissible, so the composition is structurally CERTAIN.

    ⭐ This is the predicate :func:`_type_belief` pins ``{0,0,1}`` at ``Var(log f_g) = 0`` on, and it
    applies to **both axes** — an intergenic region and an intergenic↔exon boundary are both TRAPS: no-magic-numbers, because
    RNA cannot cross a gene boundary any more than it can occupy intergenic space.

    ⚠⚠ **DO NOT CONFUSE THIS WITH ``region_init.strand_evidence``'s ``struct_lock``, which is
    deliberately REGION-ONLY.** They are two different quantities that happen to share a word:

    * *this* one answers "is the belief pinned and certain?" — a question about the belief, so both axes;
    * ``struct_lock`` answers "may this slot EMIT composition certainty into its messages?" and
      excludes G1 boundaries on purpose: a boundary is structurally gDNA but sits between RNA-carrying exons, so
      its crossing mass is RNA-contaminated and a certainty there compounds into a phantom-gDNA
      emitter. ``strand_evidence``'s own docstring carries that reasoning.

    It lives here, beside the code that applies it, so the instruments that classify objects by it read
    ONE definition rather than each re-deriving it — two homes for one predicate is how a region-only
    variant survived in two scripts and a test at the same time.
    """
    return ~np.asarray(free_pos, bool) & ~np.asarray(free_neg, bool)


#: A transcript's genomic **LOW** end. On ``+`` that terminus is the TSS, on ``−`` it is the TES — and
#: pairing the bits by which genomic end they mark, rather than by TSS/TES, is the whole point.
_RNA_LOW_END = np.uint16(FLAG_TSS_POS | FLAG_TES_NEG)
#: A transcript's genomic **HIGH** end: the TES on ``+``, the TSS on ``−``.
_RNA_HIGH_END = np.uint16(FLAG_TES_POS | FLAG_TSS_NEG)


def terminus_flank_gain(boundary_flags) -> tuple[np.ndarray, np.ndarray]:
    """``(right_gains, left_gains)`` — which FLANK of each BOUNDARY carries RNA the BOUNDARY itself cannot see.

    An BOUNDARY is a single genomic position and it counts the fragments spanning it CONTIGUOUSLY, so the
    transcripts it can see are exactly the ones continuous across it::

        T(BOUNDARY)  =  T(REGION_left)  ∩  T(REGION_right)

    A transcript whose body BEGINS at the BOUNDARY is in the right flank and not in the BOUNDARY, so
    ``T(BOUNDARY) = T(right)`` fails; one that ENDS there breaks the left equality the same way. Those are
    equalities and not containments: a composition imputation between two objects needs them to be
    measuring the same population, and ``phi_R`` too high or too low both corrupt ``phi_g``.

    ⛔⛔ **WRITTEN IN GENOMIC TERMS, AND THAT IS NOT A STYLE CHOICE — TSS/TES DOES NOT DETERMINE THE
    SIDE, BECAUSE THE STRAND FLIPS IT.** A ``+`` transcript's body extends toward higher coordinates from
    its TSS; a ``−`` transcript's extends that way from its TES. So the bits are paired by which genomic
    END they mark and the arrays are named for the FLANK. Gated on two mirror-image annotations —
    identical geometry, opposite strands — where the flank answer is the same and the TSS bit alone points
    at the opposite BOUNDARIES (`test_terminus_population_licence`).

    ⚠ **OR over both strands, deliberately.** A composition is a claim about the whole pair
    {gDNA, RNA+, RNA−}, so a population break on either strand breaks it.

    ⚠ **The two masks are INDEPENDENT, not complements** — a transcript can end where another begins, and
    then neither flank matches. ``0`` (no graph supplied) means no terminus and so both flanks match.

    ⛔ **TERMINI ONLY: DONOR/ACCEPTOR are excluded on purpose.** A splice site also changes the population
    — RNA splices out or in — but there the flux is MEASURED (``sj_count``) and the graft and the
    peel exist to route it. A terminus has no flux to measure: a transcript simply begins. That is the
    boundary between the two treatments.

    ``boundary_flags`` may be on either axis — the per-contiguous-boundary array from
    :func:`~rigel.calibration.splice_graph.build_boundary_flags_array`, or :class:`RegionStatics`'s per-slot
    copy of it (``0`` at REGION slots, which reads as "no terminus" and is correct: a REGION is not a
    position and breaks no population).
    """
    flags = np.asarray(boundary_flags, dtype=np.uint16)
    return (flags & _RNA_LOW_END) != 0, (flags & _RNA_HIGH_END) != 0


def _type_belief(free_pos, free_neg, deconv, mass_unspl):
    """Build the per-region composition ``(f_pos, f_neg, f_g)`` for ONE region type (regions OR boundaries) from its
    signature-binary classification + its strand-only solve.

    ``free_pos``/``free_neg`` are the per-region booleans for whether each strand's RNA axis is admissible (a
    region's own ±transcript bits; a boundary's ±strand CONTINUITY across the boundary). ``deconv`` is the
    strand-only :class:`RegionDeconv` (no global, no imputation). The signature-binary default is all-gDNA
    ``{0,0,1}`` (`ARCHITECTURE §3`). The class overrides:

    * **TRAPS: no-magic-numbers** (neither strand free — intergenic region / no-RNA-crossing boundary): a LOCKED gDNA sink — keep
      ``{0,0,1}``.
    * **TRAPS: one-thing-varied** (exactly one strand free, with data): the STRAND DECONVOLUTION alone resolves the pie (a
      single-strand region is 1-D: ``f_active = 1 − f_g``).
    * **TRAPS: converge-and-delete** (both strands free — AMBIG): unresolvable by strand → keep the ``{0,0,1}`` default at MAX
      (``inf``) variance; the sweep resolves it from neighbour messages + the global prior.

    Returns the six per-region arrays ``(f_pos, f_neg, f_g, var_pos, var_neg, var_gdna)`` — the composition + the
    precision state: ``var=0`` locked, ``inf`` no-information, else the strand-solve
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

    g1 = g1_locked(free_pos, free_neg)
    g2 = free_pos ^ free_neg
    g2_active = g2 & (np.asarray(mass_unspl, dtype=np.float64) > 0.0)

    # G2-active: take the strand-only solve (median f_g, mean f±, and the posterior variances). G1 sinks + TRAPS: converge-and-delete
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
class RegionStatics:
    """Per-slot STATIC structural masks (length ``n_slots``). The sweep mutates only the dynamic
    :class:`RegionBelief`; these never change.

    ⭐ **NO COUNT LIVES HERE — this class is structure only.** ``u_pos``/``u_neg``/``mass_unspliced``/
    ``mass_spliced`` used to be carried alongside :class:`RegionGeometry`'s ``mass_*``/``n_unspl_*``
    because the old accumulator's mass was fractional while the strand likelihood needed an integer.
    They are the same number now, so **all three populations sit together on the geometry** — where
    their difference is visible — and a consumer slices them there. One quantity, one place


    ``free_pos``/``free_neg`` are the nascent-RNA-active axes (a region's own ±transcript bits; an boundary's
    ±continuity — the RNA-crossing gate); ``mrna_active_pos``/``mrna_active_neg`` are the mature-RNA-active
    axes (a region's ±exon bits; an boundary's ±contiguous exon) that select the per-region solver prior.

    ``boundary_flags`` carries the splice graph's 8 structural bits (``TSS_s``/``TES_s``/``DONOR_s``/
    ``ACCEPTOR_s``) at each BOUNDARY slot, ``0`` on REGION slots.
    ⭐ **Raw bits, not pre-derived predicates.** Every consumer wants a different combination of them, and
    P1G_SCOPE's own specified predicate was measured to be nearly the COMPLEMENT of what it was meant to
    replace (plan TRAPS: a-boundary-with-rna-is-not-a-sj). Compose with :func:`~rigel.calibration.splice_graph.is_terminus` /
    :func:`~rigel.calibration.splice_graph.is_splice_site`. It is ``0`` when no graph was supplied.
    """

    n_slots: int
    free_pos: np.ndarray  # bool — nascent-RNA-active (transcript continuity); the RNA-crossing gate
    free_neg: np.ndarray  # bool
    mrna_active_pos: (
        np.ndarray
    )  # bool — mature-RNA-active (contiguous exon); selects the region prior
    mrna_active_neg: np.ndarray  # bool
    boundary_flags: np.ndarray  # uint16 — graph structural bits; 0 on REGION slots


def build_region_statics(
    chain: RegionChain,
    region_arrays,
    boundary_flags: np.ndarray | None = None,
) -> RegionStatics:
    """Gather the structural masks onto the chain, in ONE slot-keyed pass.

    ⚠ **It no longer takes a substrate**: with the counts moved onto :class:`RegionGeometry` this reads
    nothing but the signatures and the chain's own adjacency, so an unused parameter would only invite a
    caller to think the two are coupled.

    ⭐ **The region/boundary twin pair collapses.** ``_region_strand_stats`` and
    ``_boundary_strand_stats`` computed the same predicates against two different keyings, and the
    boundary one carried a ``max(left, right)`` over its two sides — a de-duplication of the straddle
    that only existed because the old accumulator deposited a crossing fragment on both sides. An boundary
    has one count, so the ``max`` goes with the sides.

    ⭐ **And the flank lookup is now the chain's own adjacency.** ``BoundarySubstrate`` carried explicit
    ``left_region``/``right_region`` arrays with ``-1`` at a reference terminal. Boundary endpoints are
    implicit — boundary ``i`` lies between region ``i`` and region ``i+1`` — and **an boundary always has a region on
    both sides**, so ``chain.left``/``chain.right`` answer it and the ``-1`` branch has no cases.

    The allow mask is the transcript-structure CONTINUITY gate: a strand-``s`` unspliced crossing is
    nascent RNA only where strand ``s`` is present on BOTH flanks. This blocks RNA at a TSS/TES
    (intergenic↔exon → neither strand continuous → a gDNA sink) and at a mixed exon↔AMBIG boundary.
    ``mrna_active_s`` is the tighter **mature**-crossing gate (contiguous exon on both flanks).

    ``boundary_flags`` is the per-contiguous-boundary ``uint16[E]`` from
    :func:`~rigel.calibration.splice_graph.build_boundary_flags_array`; ``None`` leaves the field zero,
    which every current consumer treats as "no structural information".
    """
    kind = np.asarray(chain.kind)
    obj = np.asarray(chain.obj_idx, dtype=np.int64)
    is_region = kind == REGION
    is_boundary = kind == BOUNDARY
    n = int(chain.n_slots)
    flags = _check_boundary_flags(boundary_flags, int(chain.n_boundaries_total))

    sig = np.asarray(region_arrays.signature).astype(np.int64)
    n_regions = sig.shape[0]
    region_idx = np.clip(obj, 0, max(n_regions - 1, 0))

    # the signature at each REGION slot, then read through the chain's adjacency at each BOUNDARY slot
    slot_sig = np.where(is_region, sig[region_idx] if n_regions else 0, 0)
    left = np.clip(np.asarray(chain.left), 0, max(n - 1, 0))
    right = np.clip(np.asarray(chain.right), 0, max(n - 1, 0))
    sig_l = np.where(is_boundary, slot_sig[left], 0)
    sig_r = np.where(is_boundary, slot_sig[right], 0)

    ts = np.where(is_region, np.asarray(region_arrays.strand_class)[region_idx] if n_regions else 0, -1)
    nrp_l, nrn_l = nrna_active_strands(sig_l)
    nrp_r, nrn_r = nrna_active_strands(sig_r)
    mrp_l, mrn_l = mrna_active_strands(sig_l)
    mrp_r, mrn_r = mrna_active_strands(sig_r)
    mr_self_p, mr_self_n = mrna_active_strands(slot_sig)

    free_pos = np.where(is_region, (ts == TS_POS) | (ts == TS_AMBIG), nrp_l & nrp_r)
    free_neg = np.where(is_region, (ts == TS_NEG) | (ts == TS_AMBIG), nrn_l & nrn_r)

    return RegionStatics(
        n_slots=n,
        # No per-region spliced FLOOR: spliced (mature) handling is OWNED by the message system (the
        # boundary→exon MEASUREMENT source + the exon→boundary absorption in `_scan`). A region-local floor would
        # double-count it AND inflate an boundary's UNSPLICED f_pos with mature → phantom nascent into
        # introns (matrix-confirmed: removing it is ≥ keeping it in every κ × capture × ±gDNA regime).
        free_pos=free_pos,
        free_neg=free_neg,
        mrna_active_pos=np.where(is_region, mr_self_p, mrp_l & mrp_r),
        mrna_active_neg=np.where(is_region, mr_self_n, mrn_l & mrn_r),
        boundary_flags=np.where(is_boundary, flags[np.clip(obj, 0, max(flags.shape[0] - 1, 0))], 0).astype(
            np.uint16
        ),
    )


def _check_boundary_flags(boundary_flags, n_boundaries: int) -> np.ndarray:
    """Validate the per-contiguous-boundary flags against the chain, BEFORE anything else is computed.

    A mis-sized array would shift every flag by one boundary — a defect invisible in aggregate and
    undetectable by a bit-identity gate while nothing reads the flags. Refuse it at the door.
    """
    if boundary_flags is None:
        return np.zeros(max(n_boundaries, 1), dtype=np.uint16)
    flags = np.asarray(boundary_flags, dtype=np.uint16)
    if flags.shape != (n_boundaries,):
        raise ValueError(
            f"boundary_flags has shape {flags.shape}; expected ({n_boundaries},), one per contiguous boundary. "
            f"Build it with splice_graph.build_boundary_flags_array(index) against the SAME index the "
            f"payload was scanned on. ⚠ There are no terminal slots: a reference with k regions owns "
            f"k-1 boundaries, not k+1."
        )
    return flags if flags.size else np.zeros(1, dtype=np.uint16)


def init_beliefs(
    chain: RegionChain,
    geometry: RegionGeometry,
    statics: RegionStatics,
    *,
    rna_sense_frac: float,
    gdna_strand_overdispersion: float = 0.0,
    rna_strand_overdispersion: float = 0.0,
    n_grid: int,
    n_grid_ss: int | None = None,
    logodds_window: float = 10.0,
) -> RegionBelief:
    """The signature-binary G1/G2/G3 initial :class:`RegionBelief` on the unified chain.

    All slots are strand-solved by the log-density 1-D/2-D log-odds solver (:mod:`simplex_logodds`,
    ``O(m·K)``; the bare strand likelihood + the Jeffreys reference at single-strand regions; NO global prior,
    NO imputation — those enter the sweep, P2/P3). The signature-binary class overrides (:func:`_type_belief`)
    then set the G1/G2/G3 belief. Single-strand introns resolve to ``f_g≈0`` from the BB tilt alone (the
    zero-gDNA gate); intergenic / TSS sinks lock at ``{0,0,1}``; AMBIG regions hold ``{0,0,1}`` at MAX variance
    for the sweep."""
    st = statics
    count = np.asarray(geometry.unspliced_count, np.float64)
    deconv = _solve_regions_logodds_all(
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
    return RegionBelief(
        f_pos=f_pos, f_neg=f_neg, f_g=f_g, var_pos=var_p, var_neg=var_n, var_gdna=var_g
    )


# ---------------------------------------------------------------------------
# Region-type helper (the sweep itself lives in sweep.solve_chain).
# ---------------------------------------------------------------------------


def _region_region_type(chain, region_arrays):
    """Per-slot coarse region type at REGION slots (0=intergenic, 1=intron, 2=exon; exon>intron), −1 at BOUNDARY
    slots; plus the per-region type array. Single source of truth: :func:`signature.coarse_type_array`."""
    kind = np.asarray(chain.kind)
    idx = np.asarray(chain.obj_idx, dtype=np.int64)
    rtype = coarse_type_array(np.asarray(region_arrays.signature)).astype(np.int64)  # per region
    ri_ = np.clip(idx, 0, max(rtype.shape[0] - 1, 0))
    return np.where(kind == REGION, rtype[ri_], -1), rtype


# ---------------------------------------------------------------------------
# ⛔ THE gDNA CAPTURE-CLASS LANDSCAPE WAS BUILT HERE AND DELETED — it is INERT. Do not rebuild it.
# ---------------------------------------------------------------------------
#
# The reasoning was sound and the conclusion is the interesting part. gDNA is uniform BEFORE capture, so a
# gDNA density claim transports between adjacent objects by the ratio of their CAPTURE EFFICIENCIES and
# nothing else; that efficiency is fixed by probe geometry, probes are designed from the annotation, so it
# is a property of the object's structural class; and it is directly observable on any class with
# structurally-pure-gDNA members (intergenic REGIONS, ``intergenic|exon`` BOUNDARIES). Five classes were
# implemented — REGION/BOUNDARY x off-probe / half-covered / fully-covered — with the pooled rate
# ``Σcount/ΣE`` per class supplying ``r_g = rate[class(dst)]/rate[class(src)]``.
#
# ⭐⭐ **MEASURED INERT, AND THE REASON IS AN OPERATOR THAT ALREADY EXISTS.** The relay's MASS PIN
# (the scan's ``Σ_c rho_c·E_c = M``) rescales the running level to each object's OWN observed total,
# and at a structurally pure-gDNA object that total IS its gDNA density — measured at its own capture
# stratum. So the landscape is already carried, per object and locally, by the pure-gDNA population's own
# measurements; a pooled class ratio only re-derives it, worse. A/B on the ladder: **byte-identical off
# capture**, and 1.2 % of one class on one capture-ON condition (``region/exon`` mwae 0.2719 → 0.2686).
# That does not pay for five constants, a helper and two gates.
#
# ⭐ The rule that remains is therefore one boundary and is in `messages.head`: a gDNA LEVEL crosses a
# composition-unlicensed hop **unscaled**. Capture-OFF is not a case — it is the same expression.
