"""Phase 1 — the density model ("count clue"): per-node gDNA density from OBSERVED counts.

Acyclic by construction: the gDNA density is read **directly** from count-observable nodes
(where fragments are gDNA by construction) and **swept** to every other node. It never
consults a global ``gdna_density_global * L`` product, so there is no
density->deconv->density feedback loop.

Count-observability is a property of the region **signature** (4-bit exon/intron ± flags):

* **region** is observable ⇔ it has **no exon bit** (intergenic or intron-only). Its
  contained unspliced mass is gDNA (+ nascent RNA — an upper bias the Phase-3 strand clue
  removes); an exonic region's contained mass is contaminated by mature RNA.
* **boundary** is observable ⇔ **no exon bit is shared** across its two sides → no single
  exon-strand continues across it → no *unspliced mature RNA* crosses (a single-exon
  transcript spanning the seam would put unspliced mature RNA there). Its crossing
  unspliced mass is then gDNA(+nascent).

Everything else — exonic regions, exon-spanning boundaries, AMBIG — carries no direct
gDNA observation and is **imputed by the alternating region↔boundary sweep** (observable
node → conduct across boundaries → impute → iterate). This sweep is mandatory: in an
overlapping locus (a single-exon ``+`` transcript over a multi-exon ``−`` one) the only
observable nodes can be the two locus-edge boundaries, and the whole interior is filled
inward from them.

Counts → density via the gDNA FL effective lengths: contained mass ÷ ``region_eff_len``,
crossing mass ÷ ``fl_mean``. For uniform gDNA at a given density both recover that density.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .signature import BIT_EXON_NEG, BIT_EXON_POS

_EXON_BITS = BIT_EXON_POS | BIT_EXON_NEG
# Unit pseudo-fragment for the boundary conduit weight (weakly-informative — a boundary with
# little crossing traffic conducts density weakly).
_TRAFFIC_PSEUDOCOUNT = 1.0


@dataclass(frozen=True, slots=True)
class NodeDensity:
    """Per-region gDNA density (count clue) after the sweep."""

    density: np.ndarray  # float64[R] — local gDNA density (fragments per effective bp)
    gdna_mass: np.ndarray  # float64[R] — density × region_eff_len (count-clue gDNA mass)
    count_evidence: (
        np.ndarray
    )  # float64[R] — density·eff_len: expected gDNA count (count-prior precision)
    region_count_observable: np.ndarray  # bool[R] — count-observable region (non-exonic)
    boundary_count_observable: np.ndarray  # bool[R] — count-observable boundary right of region r
    n_region_count_observable: int
    n_boundary_count_observable: int


def count_observable_masks(
    signature: np.ndarray, ref_id: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Signature-based count-observability for regions and (right-) boundaries.

    Returns ``(region_count_observable, boundary_count_observable)``, both ``bool[R]``. ``boundary_count_observable[r]``
    describes the internal boundary between region ``r`` and ``r+1`` (defined iff same ref).
    """
    sig = np.asarray(signature).astype(np.int64)
    ref = np.asarray(ref_id)
    r = sig.shape[0]
    region_count_observable = (sig & _EXON_BITS) == 0
    boundary_count_observable = np.zeros(r, dtype=bool)
    if r > 1:
        same = ref[:-1] == ref[1:]
        shared_exon = (sig[:-1] & sig[1:] & _EXON_BITS) != 0
        boundary_count_observable[:-1] = same & ~shared_exon
    return region_count_observable, boundary_count_observable


def node_gdna_density(
    substrate,
    region_arrays,
    region_eff_len: np.ndarray,
    fl_mean: float,
    gdna_frac: np.ndarray | None = None,
) -> NodeDensity:
    """Per-region gDNA density from the count clue + the region↔boundary density sweep.

    ``gdna_frac`` is the per-region strand-deconvolved gDNA fraction of the contained mass
    (Phase-2 strand clue). Supplying it cleans the nascent-RNA upper bias from the
    count-observable nodes **before** gdna_density_global is read, so the swept density is clean gDNA, not
    gDNA+nascent — the empirical-Bayes estimate of the global density hyperparameter. ``None``
    ⇒ all 1.0 (the raw, pre-clean behaviour).
    """
    sig = np.asarray(region_arrays.signature)
    ref_id = np.asarray(region_arrays.ref_id)
    ref_offsets = np.asarray(region_arrays.ref_offsets, dtype=np.int64)
    region_eff_len = np.asarray(region_eff_len, dtype=np.float64)
    r = sig.shape[0]
    region_count_observable, boundary_count_observable = count_observable_masks(sig, ref_id)
    gdna_frac = np.ones(r) if gdna_frac is None else np.asarray(gdna_frac, dtype=np.float64)

    # --- direct observations ---
    # region node: strand-cleaned contained gDNA mass (observable regions only).
    reg_mass = np.where(
        region_count_observable, gdna_frac * substrate.contained.mass_unspliced, 0.0
    )
    reg_len = np.where(region_count_observable, region_eff_len, 0.0)
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
    bnd_mass = np.where(
        boundary_count_observable, cross_mass, 0.0
    )  # gDNA crossing mass (observable boundaries)
    bnd_len = np.where(boundary_count_observable, fl_mean, 0.0)
    # conduit reliability: a boundary with more crossing traffic propagates density better.
    weight = np.where(same_right, cross_flux / (cross_flux + _TRAFFIC_PSEUDOCOUNT), 0.0)

    # --- alternating region↔boundary density sweep (per reference) ---
    # Two parallel tracks per node: gDNA mass and effective length.
    from_left_mass = np.zeros(r)
    from_left_len = np.zeros(r)
    from_right_mass = np.zeros(r)
    from_right_len = np.zeros(r)
    for f in range(ref_offsets.shape[0] - 1):
        s, e = int(ref_offsets[f]), int(ref_offsets[f + 1])
        if e <= s:
            continue
        run_mass = run_len = 0.0  # forward: left-side evidence, decayed per boundary
        for i in range(s, e):
            if i > s:
                w = weight[i - 1]
                run_mass = w * (run_mass + bnd_mass[i - 1])
                run_len = w * (run_len + bnd_len[i - 1])
            from_left_mass[i] = run_mass
            from_left_len[i] = run_len
            run_mass += reg_mass[i]
            run_len += reg_len[i]
        run_mass = run_len = 0.0  # reverse
        for i in range(e - 1, s - 1, -1):
            if i < e - 1:
                w = weight[i]
                run_mass = w * (run_mass + bnd_mass[i])
                run_len = w * (run_len + bnd_len[i])
            from_right_mass[i] = run_mass
            from_right_len[i] = run_len
            run_mass += reg_mass[i]
            run_len += reg_len[i]

    # density = own + swept-neighbour evidence: swept gDNA mass / swept effective length.
    swept_mass = reg_mass + from_left_mass + from_right_mass
    swept_len = reg_len + from_left_len + from_right_len
    # global fallback for nodes the sweep never reached (no count-observable evidence in the ref).
    total_len = float(reg_len.sum() + bnd_len.sum())
    seed_density = float(reg_mass.sum() + bnd_mass.sum()) / total_len if total_len > 0.0 else 0.0
    density = np.where(swept_len > 0.0, swept_mass / np.maximum(swept_len, 1e-12), seed_density)
    gdna_mass = density * region_eff_len
    # count-prior precision = the expected gDNA count (density·eff_len, which equals gdna_mass):
    # the count clue is only as confident as the number of gDNA events it expects to see, so it
    # defers to the strand clue where RNA dominates (imputed-low density => small expected count =>
    # weak prior, Jeffreys floor in joint_deconv). At RNA-rich exons it is small and strand governs.
    count_evidence = gdna_mass
    return NodeDensity(
        density=density,
        gdna_mass=gdna_mass,
        count_evidence=count_evidence,
        region_count_observable=region_count_observable,
        boundary_count_observable=boundary_count_observable,
        n_region_count_observable=int(region_count_observable.sum()),
        n_boundary_count_observable=int(boundary_count_observable.sum()),
    )


__all__ = ["NodeDensity", "count_observable_masks", "node_gdna_density"]
