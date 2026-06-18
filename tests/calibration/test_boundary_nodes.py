"""Boundary-node solve (`calibration.boundary_nodes`) — the bipartite R↔B↔R graph, Phase 2.

A boundary node is solved by the SAME per-node core as a region (`_solve_nodes`), with boundary-specific
sufficient statistics. Pins:
  1. an intergenic-flanked boundary → the (0,0,1) gDNA vertex (the forbid mask);
  2. a single-strand {POS,POS} seam is resolved by its strand (sense-tilted → RNA, balanced → gDNA);
  3. the count `N` is ONE side, NOT the sum — a boundary with `flux_left == flux_right` has the same posterior
     as a region carrying that many distinct fragments (no √2 over-sharpening);
  4. the one-sided motif-stranded spliced floor (§4.16): the sense spliced on a strand's exon flank pins that
     strand's fraction, read from ONE side only.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from rigel.calibration.boundary_nodes import (
    _spliced_floor,
    boundary_node_class,
    deconv_boundaries_sweep,
)
from rigel.calibration.signature import (
    BIT_EXON_POS,
    BIT_INTRON_POS,
    TS_NEG,
    TS_NONE,
    TS_POS,
)
from rigel.calibration.simplex_sweep import deconv_regions_sweep


def _view(pos, neg, *, spl_sense=None, spl_anti=None):
    pos = np.asarray(pos, dtype=np.float64)
    neg = np.asarray(neg, dtype=np.float64)
    n = pos.size
    ss = np.zeros(n) if spl_sense is None else np.asarray(spl_sense, dtype=np.float64)
    sa = np.zeros(n) if spl_anti is None else np.asarray(spl_anti, dtype=np.float64)
    return SimpleNamespace(
        n_unspliced_pos=pos, n_unspliced_neg=neg, n_spliced_sense=ss, n_spliced_antisense=sa,
        mass_unspliced=pos + neg, mass_spliced=ss + sa,
    )


def _bsub(left_region, right_region, left, right):
    return SimpleNamespace(
        left_region=np.asarray(left_region), right_region=np.asarray(right_region),
        left=left, right=right,
    )


def _ra(strand_class, signature):
    return SimpleNamespace(
        strand_class=np.asarray(strand_class), signature=np.asarray(signature, dtype=np.int64),
    )


def _solve_one_boundary(ts_l, ts_r, sig_l, sig_r, lpos, lneg, rpos, rneg, *, spl_sense=(0.0, 0.0),
                        kappa=0.99):
    """One boundary b0 between region 0 (left) and region 1 (right). Returns its f_g."""
    ra = _ra([ts_l, ts_r], [sig_l, sig_r])
    bs = _bsub([0], [1],
               _view([lpos], [lneg], spl_sense=[spl_sense[0]]),
               _view([rpos], [rneg], spl_sense=[spl_sense[1]]))
    return deconv_boundaries_sweep(bs, ra, rna_sense_frac=kappa, n_grid=40)


def test_node_class_flank_or_and_consistent_sense():
    tl = np.array([TS_POS, TS_POS, TS_POS, TS_NONE, TS_NEG])
    tr = np.array([TS_POS, TS_NONE, TS_NEG, TS_NONE, TS_NEG])
    ap, an, so = boundary_node_class(tl, tr)
    # {POS,POS}, {POS,NONE}: allow_pos, strand-observable. {POS,NEG}: both allowed, NOT observable.
    assert list(ap) == [True, True, True, False, False]
    assert list(an) == [False, False, True, False, True]
    assert list(so) == [True, True, False, False, True]  # {NONE,NONE} not observable; {NEG,NEG} is


def test_spliced_floor_reads_one_side_only():
    # + exon (left) → + intron (right): the spliced is on the LEFT (exon) side; never the sum.
    sig_l = np.array([BIT_EXON_POS])
    sig_r = np.array([BIT_INTRON_POS])
    pos = _spliced_floor(sig_l, sig_r, np.array([7.0]), np.array([3.0]), BIT_EXON_POS, BIT_INTRON_POS)
    assert pos[0] == 7.0  # the LEFT (exon) side only — not 7+3
    # the − strand has no junction here → 0
    from rigel.calibration.signature import BIT_EXON_NEG, BIT_INTRON_NEG
    neg = _spliced_floor(sig_l, sig_r, np.array([7.0]), np.array([3.0]), BIT_EXON_NEG, BIT_INTRON_NEG)
    assert neg[0] == 0.0


def test_intergenic_boundary_locks_to_gdna_vertex():
    # both flanks intergenic ⇒ allow_pos = allow_neg = False ⇒ the forbid mask leaves only (0,0,1).
    r = _solve_one_boundary(TS_NONE, TS_NONE, 0, 0, 50.0, 50.0, 50.0, 50.0)
    assert r.gdna_frac[0] > 0.99


def test_single_strand_seam_resolves_by_strand():
    # {POS,POS} exon↔intron seam, ~all-sense (99/1) ⇒ pure +RNA ⇒ f_g ≈ 0.
    rna = _solve_one_boundary(TS_POS, TS_POS, BIT_EXON_POS, BIT_INTRON_POS, 99.0, 1.0, 99.0, 1.0)
    assert rna.gdna_frac[0] < 0.10
    # balanced 50/50 ⇒ impossible for pure +RNA at κ=0.99 ⇒ gDNA.
    gd = _solve_one_boundary(TS_POS, TS_POS, BIT_EXON_POS, BIT_INTRON_POS, 50.0, 50.0, 50.0, 50.0)
    assert gd.gdna_frac[0] > 0.5


def test_count_is_one_side_not_summed():
    # A boundary with flux_left == flux_right == (50,50) has N = 100 (one side), NOT 200. Its posterior must
    # equal a single-strand-+ REGION carrying (50,50) contained — same N ⇒ same f_g and f_g_var (no √2
    # over-sharpening from doubling).
    bnd = _solve_one_boundary(TS_POS, TS_POS, BIT_EXON_POS, BIT_INTRON_POS, 50.0, 50.0, 50.0, 50.0)
    contained = _view([50.0], [50.0])
    reg = deconv_regions_sweep(
        SimpleNamespace(contained=contained), _ra([TS_POS], [BIT_EXON_POS]),
        rna_sense_frac=0.99, n_grid=40,
    )
    assert np.isclose(bnd.gdna_frac[0], reg.gdna_frac[0], atol=1e-9)
    assert np.isclose(bnd.gdna_frac_var[0], reg.gdna_frac_var[0], atol=1e-9)


def test_spliced_floor_pins_the_strand_fraction():
    # A + exon↔intron boundary whose crossing reads ~balanced (would read gDNA) but with a strong + spliced
    # mature floor on the exon side ⇒ f_pos is pinned up ⇒ f_g pulled below the no-floor case.
    floored = _solve_one_boundary(TS_POS, TS_POS, BIT_EXON_POS, BIT_INTRON_POS, 50.0, 50.0, 50.0, 50.0,
                                  spl_sense=(45.0, 0.0))
    bare = _solve_one_boundary(TS_POS, TS_POS, BIT_EXON_POS, BIT_INTRON_POS, 50.0, 50.0, 50.0, 50.0)
    assert floored.gdna_frac[0] < bare.gdna_frac[0]  # the mature floor forces RNA ⇒ less gDNA


def test_terminal_boundary_uses_present_flank():
    # left_region = -1 (ref-start terminal): the off-edge side is empty; the present (right) flank determines
    # the node class and no -1 index is dereferenced.
    ra = _ra([TS_POS], [BIT_EXON_POS])
    bs = _bsub([-1], [0], _view([0.0], [0.0]), _view([80.0], [80.0]))
    r = deconv_boundaries_sweep(bs, ra, rna_sense_frac=0.99, n_grid=40)
    assert np.isfinite(r.gdna_frac[0]) and 0.0 <= r.gdna_frac[0] <= 1.0
