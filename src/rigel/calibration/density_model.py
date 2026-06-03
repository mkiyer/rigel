"""Phase 1 — the density model ("count clue"): per-node gDNA density from OBSERVED counts.

Acyclic by construction: the gDNA density is read **directly** from count-decodable nodes
(where fragments are gDNA by construction) and **swept** to every other node. It never
consults the global ``ρ_0·ω·L`` — so there is no ``ρ_0 → decode → ρ_0`` feedback loop.

Count-decodability is a property of the region **signature** (4-bit exon/intron ± flags):

* **region** is decodable ⇔ it has **no exon bit** (intergenic or intron-only). Its
  contained unspliced mass is gDNA (+ nascent RNA — an upper bias the Phase-3 strand clue
  removes); an exonic region's contained mass is contaminated by mature RNA.
* **boundary** is decodable ⇔ **no exon bit is shared** across its two sides → no single
  exon-strand continues across it → no *unspliced mature RNA* crosses (a single-exon
  transcript spanning the seam would put unspliced mature RNA there). Its crossing
  unspliced mass is then gDNA(+nascent).

Everything else — exonic regions, exon-spanning boundaries, AMBIG — carries no direct
gDNA observation and is **imputed by the alternating region↔boundary sweep** (decodable
node → conduct across boundaries → impute → iterate). This sweep is mandatory: in an
overlapping locus (a single-exon ``+`` transcript over a multi-exon ``−`` one) the only
decodable nodes can be the two locus-edge boundaries, and the whole interior is filled
inward from them.

Counts → density via the gDNA FL effective lengths: contained mass ÷ ``region_eff_len``,
crossing mass ÷ ``μ_FL``. For uniform gDNA at density ``ρ`` both recover ``ρ``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .signature import BIT_EXON_NEG, BIT_EXON_POS

_EXON_BITS = BIT_EXON_POS | BIT_EXON_NEG
# Unit pseudo-fragment for the boundary conduit weight (weakly-informative; matches the
# sweep.py convention — a boundary with little crossing traffic conducts density weakly).
_TRAFFIC_PSEUDOCOUNT = 1.0


@dataclass(frozen=True, slots=True)
class NodeDensity:
    """Per-region gDNA density (count clue) after the sweep."""

    density: np.ndarray  # float64[R] — local gDNA density (fragments per effective bp)
    gdna_mass: np.ndarray  # float64[R] — density × region_eff_len (count-clue gDNA mass)
    count_evidence: np.ndarray  # float64[R] — swept DISCRETE gDNA-rep flux (count-prior precision)
    region_decodable: np.ndarray  # bool[R] — count-decodable region (non-exonic)
    boundary_decodable: np.ndarray  # bool[R] — decodable boundary to the right of region r
    n_region_dec: int
    n_boundary_dec: int


def count_decodable_masks(
    signature: np.ndarray, ref_id: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Signature-based count-decodability for regions and (right-) boundaries.

    Returns ``(region_decodable, boundary_decodable)``, both ``bool[R]``. ``boundary_decodable[r]``
    describes the internal boundary between region ``r`` and ``r+1`` (defined iff same ref).
    """
    sig = np.asarray(signature).astype(np.int64)
    ref = np.asarray(ref_id)
    r = sig.shape[0]
    region_dec = (sig & _EXON_BITS) == 0
    bnd_dec = np.zeros(r, dtype=bool)
    if r > 1:
        same = ref[:-1] == ref[1:]
        shared_exon = (sig[:-1] & sig[1:] & _EXON_BITS) != 0
        bnd_dec[:-1] = same & ~shared_exon
    return region_dec, bnd_dec


def node_gdna_density(
    substrate,
    region_arrays,
    region_eff_len: np.ndarray,
    mu_fl: float,
) -> NodeDensity:
    """Per-region gDNA density from the count clue + the region↔boundary density sweep."""
    sig = np.asarray(region_arrays.signature)
    ref_id = np.asarray(region_arrays.ref_id)
    ref_offsets = np.asarray(region_arrays.ref_offsets, dtype=np.int64)
    region_eff_len = np.asarray(region_eff_len, dtype=np.float64)
    r = sig.shape[0]
    region_dec, bnd_dec = count_decodable_masks(sig, ref_id)

    # --- direct observations ---
    # region node: contained gDNA-representative mass (decodable regions only).
    a_reg = np.where(region_dec, substrate.contained.mass_unspliced, 0.0)
    b_reg = np.where(region_dec, region_eff_len, 0.0)
    # boundary node (indexed by left region r): the unspliced crossing mass/flux at the
    # r/r+1 seam = region r's RIGHT view + region r+1's LEFT view (the two sides).
    same_right = np.zeros(r, dtype=bool)
    if r > 1:
        same_right[:-1] = ref_id[:-1] == ref_id[1:]
    cross_mass = np.zeros(r, dtype=np.float64)
    cross_flux = np.zeros(r, dtype=np.float64)
    if r > 1:
        cross_mass[:-1] = substrate.right.mass_unspliced[:-1] + substrate.left.mass_unspliced[1:]
        cross_flux[:-1] = substrate.right.n_unspliced[:-1].astype(
            np.float64
        ) + substrate.left.n_unspliced[1:].astype(np.float64)
    cross_mass = np.where(same_right, cross_mass, 0.0)
    cross_flux = np.where(same_right, cross_flux, 0.0)
    a_bnd = np.where(bnd_dec, cross_mass, 0.0)  # gDNA crossing mass (decodable boundaries)
    b_bnd = np.where(bnd_dec, mu_fl, 0.0)
    # conduit reliability: a boundary with more crossing traffic propagates density better.
    weight = np.where(same_right, cross_flux / (cross_flux + _TRAFFIC_PSEUDOCOUNT), 0.0)
    # DISCRETE gDNA-rep flux (statistical power; swept alongside the fractional mass so the
    # count prior's precision is a fragment count, never the fractional mass — review §3).
    c_reg = np.where(region_dec, substrate.contained.n_unspliced.astype(np.float64), 0.0)
    c_bnd = np.where(bnd_dec, cross_flux, 0.0)

    # --- alternating region↔boundary sweep (per reference; cf. sweep.py) ---
    # Three parallel tracks per node: a = gDNA mass, b = effective length, c = DISCRETE flux.
    fl_a = np.zeros(r)
    fl_b = np.zeros(r)
    fl_c = np.zeros(r)
    fr_a = np.zeros(r)
    fr_b = np.zeros(r)
    fr_c = np.zeros(r)
    for f in range(ref_offsets.shape[0] - 1):
        s, e = int(ref_offsets[f]), int(ref_offsets[f + 1])
        if e <= s:
            continue
        run_a = run_b = run_c = 0.0  # forward: left-side evidence, decayed per boundary
        for i in range(s, e):
            if i > s:
                w = weight[i - 1]
                run_a = w * (run_a + a_bnd[i - 1])
                run_b = w * (run_b + b_bnd[i - 1])
                run_c = w * (run_c + c_bnd[i - 1])
            fl_a[i] = run_a
            fl_b[i] = run_b
            fl_c[i] = run_c
            run_a += a_reg[i]
            run_b += b_reg[i]
            run_c += c_reg[i]
        run_a = run_b = run_c = 0.0  # reverse
        for i in range(e - 1, s - 1, -1):
            if i < e - 1:
                w = weight[i]
                run_a = w * (run_a + a_bnd[i])
                run_b = w * (run_b + b_bnd[i])
                run_c = w * (run_c + c_bnd[i])
            fr_a[i] = run_a
            fr_b[i] = run_b
            fr_c[i] = run_c
            run_a += a_reg[i]
            run_b += b_reg[i]
            run_c += c_reg[i]

    # density = own + swept-neighbour evidence (α = gDNA mass, β = effective length).
    alpha = a_reg + fl_a + fr_a
    beta = b_reg + fl_b + fr_b
    count_evidence = c_reg + fl_c + fr_c  # discrete fragment count behind the density
    # global fallback for nodes the sweep never reached (no decodable evidence in the ref).
    tot_b = float(b_reg.sum() + b_bnd.sum())
    seed_density = float(a_reg.sum() + a_bnd.sum()) / tot_b if tot_b > 0.0 else 0.0
    density = np.where(beta > 0.0, alpha / np.maximum(beta, 1e-12), seed_density)
    gdna_mass = density * region_eff_len
    return NodeDensity(
        density=density,
        gdna_mass=gdna_mass,
        count_evidence=count_evidence,
        region_decodable=region_dec,
        boundary_decodable=bnd_dec,
        n_region_dec=int(region_dec.sum()),
        n_boundary_dec=int(bnd_dec.sum()),
    )


__all__ = ["NodeDensity", "count_decodable_masks", "node_gdna_density"]
