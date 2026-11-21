"""A locked object reads its strands.

An object whose annotation admits no RNA strand is locked as gDNA. Where a stranded library puts one-stranded fragments
on it (a transcript the annotation ends too soon, a read placed at the wrong copy of a duplicated sequence), its two
strands say they are not gDNA: gDNA lands on either strand evenly. The lock then holds what its strands can share, and
the exon beside a locked boundary is not floored at the boundary's count. Where the strands say nothing, and where
every locked object is balanced, the lock stands as it was.

Each test runs the whole solve (`calibrate._sweep_inputs` and `_solve`) on a library of one-exon genes between
intergenic regions, so every exon has a locked boundary on each side.
"""

from __future__ import annotations

import importlib
from types import SimpleNamespace

import numpy as np
import pytest
from _synthetic import delta_pmf, make_chain_parts

from rigel.calibration.region_chain import BOUNDARY, REGION
from rigel.calibration.signature import BIT_EXON_NEG, BIT_EXON_POS
from rigel.calibration.strand_balance import StrandBalance
from rigel.config import CONSTANTS, CalibrationConfig

CAL = importlib.import_module("rigel.calibration.calibrate")

KAPPA = 0.99
LIVE = StrandBalance(KAPPA, 10_000, False)
DEAD = StrandBalance(0.5, 10_000, False)
N_GENES = 400
INTERGENIC, EXON, FRAGMENT = 4000, 1000, 200
BACKGROUND = 5e-3  # gDNA per placement
RNA = 2000  # unspliced RNA fragments per exon
OVERHANG = 40  # one-stranded fragments on a marked gene's right boundary
BEYOND = 60  # … and in the intergenic region beyond it
MARKED = np.arange(N_GENES) % 20 == 7  # one gene in twenty


def _library(seed: int, *, overhang: bool, strand: StrandBalance):
    """``N_GENES`` one-exon genes of alternating strand, an intergenic region before each and after the last. gDNA
    lands at ``BACKGROUND`` per placement everywhere, its strands balanced; each exon holds ``RNA`` fragments on its
    strand at κ. With ``overhang`` a marked gene's right boundary and the region beyond it also hold fragments on the
    gene's own strand, as a transcript running past its annotated end would leave them."""
    rng = np.random.default_rng(seed)
    kappa = float(strand.rna_sense_frac)
    sig, size, plus, minus, b_plus, b_minus = [], [], [], [], [], []

    def gdna(placements):
        g = rng.poisson(BACKGROUND * placements)
        p = rng.binomial(g, 0.5)
        return p, g - p

    def rna(n, on_plus):
        s = rng.binomial(n, kappa)
        return (s, n - s) if on_plus else (n - s, s)

    def add(target_plus, target_minus, counts):
        target_plus.append(counts[0])
        target_minus.append(counts[1])

    def both(a, b):
        return a[0] + b[0], a[1] + b[1]

    for k in range(N_GENES):
        on_plus = k % 2 == 0
        marked_before = overhang and k > 0 and bool(MARKED[k - 1])
        before = gdna(INTERGENIC - FRAGMENT + 1)
        if marked_before:
            before = both(before, rna(BEYOND, (k - 1) % 2 == 0))
        sig.append(0)
        size.append(INTERGENIC)
        add(plus, minus, before)
        add(b_plus, b_minus, gdna(FRAGMENT - 1))  # the gene's left boundary
        sig.append(BIT_EXON_POS if on_plus else BIT_EXON_NEG)
        size.append(EXON)
        add(plus, minus, both(gdna(EXON - FRAGMENT + 1), rna(RNA, on_plus)))
        right = gdna(FRAGMENT - 1)
        if overhang and MARKED[k]:
            right = both(right, rna(OVERHANG, on_plus))
        add(b_plus, b_minus, right)
    last = gdna(INTERGENIC - FRAGMENT + 1)
    if overhang and MARKED[N_GENES - 1]:
        last = both(last, rna(BEYOND, (N_GENES - 1) % 2 == 0))
    sig.append(0)
    size.append(INTERGENIC)
    add(plus, minus, last)
    return make_chain_parts(
        sig,
        region_size_bp=np.asarray(size, float),
        region_pos=np.asarray(plus, float),
        region_neg=np.asarray(minus, float),
        boundary_pos=np.asarray(b_plus, float),
        boundary_neg=np.asarray(b_minus, float),
        gdna_fl=delta_pmf(FRAGMENT, 400),
        rna_fl=delta_pmf(150, 400),
    )


def _solve(monkeypatch, parts, strand: StrandBalance):
    """The whole solve on ``parts``: calibration's own inputs and its two phases, the chain and its geometry handed in
    where `calibrate` would build them from a payload."""
    monkeypatch.setattr(CAL, "build_region_chain", lambda *a, **k: parts.chain)
    monkeypatch.setattr(CAL, "build_region_geometry", lambda *a, **k: parts.geometry)
    monkeypatch.setattr(CAL, "build_region_statics", lambda *a, **k: parts.statics)
    monkeypatch.setattr(CAL, "fit_strand_balance", lambda *a, **k: strand)
    payload = SimpleNamespace(ref_region_offsets=None, ref_boundary_offsets=None)
    inputs = CAL._sweep_inputs(
        payload,
        parts.substrate,
        parts.region_arrays,
        None,
        delta_pmf(FRAGMENT, 400),
        delta_pmf(150, 400),
        None,
        None,
        CalibrationConfig(),
    )
    belief, _ = CAL._solve(inputs, None)
    return belief, inputs


def _slots(parts):
    """The library's slots by role: each gene's exon and its right boundary, and every locked slot."""
    kind = np.asarray(parts.chain.kind)
    obj = np.asarray(parts.chain.obj_idx)
    region_slot = np.flatnonzero(kind == REGION)[np.argsort(obj[kind == REGION])]
    boundary_slot = np.flatnonzero(kind == BOUNDARY)[np.argsort(obj[kind == BOUNDARY])]
    genes = np.arange(N_GENES)
    locked = ~(np.asarray(parts.statics.free_pos, bool) | np.asarray(parts.statics.free_neg, bool))
    return SimpleNamespace(
        exon=region_slot[2 * genes + 1],
        right=boundary_slot[2 * genes + 1],
        beyond=region_slot[2 * genes + 2],
        locked=locked,
        count=np.asarray(parts.geometry.unspliced_count, np.float64).sum(axis=1),
    )


def test_one_stranded_fragments_on_a_locked_boundary_are_not_gdna(monkeypatch):
    """A marked gene's right boundary holds about forty fragments on the gene's strand and one of gDNA. Its strands
    cannot share them, so it books a few fragments, not forty; and the region beyond, with sixty more on one strand
    over nineteen of gDNA, books about its gDNA."""
    parts = _library(11, overhang=True, strand=LIVE)
    belief, _ = _solve(monkeypatch, parts, LIVE)
    s = _slots(parts)
    booked = np.asarray(belief.f_g) * s.count
    assert booked[s.right[MARKED]].mean() < 0.25 * s.count[s.right[MARKED]].mean(), (
        booked[s.right[MARKED]].mean(),
        s.count[s.right[MARKED]].mean(),
    )
    expected = BACKGROUND * (INTERGENIC - FRAGMENT + 1)
    assert booked[s.beyond[MARKED]].mean() < 2.0 * expected, (
        booked[s.beyond[MARKED]].mean(),
        expected,
    )


def test_the_remainder_is_rna_on_the_genes_strand_and_the_lock_has_no_location(monkeypatch):
    """What a marked boundary's strands cannot share is booked as RNA of the strand it sits on, the gene's, so the
    simplex closes. Where its other strand holds one fragment or none, its share is a fragment or two give or take as
    much: wider than the landscape's location floor, so the boundary does not train the landscape."""
    parts = _library(11, overhang=True, strand=LIVE)
    belief, _ = _solve(monkeypatch, parts, LIVE)
    s = _slots(parts)
    f_pos, f_neg, f_g = (np.asarray(x) for x in (belief.f_pos, belief.f_neg, belief.f_g))
    np.testing.assert_allclose((f_pos + f_neg + f_g)[s.locked], 1.0, rtol=0, atol=1e-12)
    plus_gene = np.arange(N_GENES) % 2 == 0
    on_plus, on_minus = s.right[MARKED & plus_gene], s.right[MARKED & ~plus_gene]
    assert np.all(f_pos[on_plus] > 0.5) and np.all(f_neg[on_plus] < 0.05)
    assert np.all(f_neg[on_minus] > 0.5) and np.all(f_pos[on_minus] < 0.05)
    var = np.asarray(belief.var_gdna)
    marked = s.right[MARKED]
    bare = np.asarray(parts.geometry.unspliced_count, np.float64)[marked].min(axis=1) <= 1.0
    assert bare.sum() >= 10 and np.all(var[marked[bare]] > CONSTANTS.landscape.located_var), (
        bare.sum(),
        var[marked[bare]].min(),
    )
    assert np.all(np.isfinite(var[s.locked]))


def test_the_exon_beside_it_is_not_floored_at_its_count(monkeypatch):
    """The edge rule asks an exon for its locked boundary's density. With the overhang counted as gDNA a marked gene's
    exon is asked for 160 fragments and books tens against its own strands; with the boundary reading its strands the
    exon books what the unmarked genes' exons book, within a factor of two."""
    parts = _library(11, overhang=True, strand=LIVE)
    belief, _ = _solve(monkeypatch, parts, LIVE)
    s = _slots(parts)
    booked = np.asarray(belief.f_g) * s.count
    marked, plain = booked[s.exon[MARKED]].mean(), booked[s.exon[~MARKED]].mean()
    assert marked < 2.0 * plain, (marked, plain)


def test_a_lock_holding_only_gdna_keeps_it(monkeypatch):
    """gDNA until proven otherwise. In the same library the locked objects that hold only gDNA keep it: a boundary's
    one or two fragments say nothing against the lock, and an intergenic region's nineteen are one-sided only by
    chance. The departed twentieth costs the boundaries under a fragment in a hundred and the regions under two; a lock
    holding one fragment keeps it."""
    parts = _library(11, overhang=True, strand=LIVE)
    belief, _ = _solve(monkeypatch, parts, LIVE)
    s = _slots(parts)
    plain = s.locked.copy()
    plain[s.right[MARKED]] = False
    plain[s.beyond[MARKED]] = False
    booked = np.asarray(belief.f_g) * s.count
    is_region = np.asarray(parts.chain.kind) == REGION
    for name, sel, kept in (
        ("regions", plain & is_region, 0.98),
        ("boundaries", plain & ~is_region, 0.99),
    ):
        assert booked[sel].sum() > kept * s.count[sel].sum(), (
            name,
            booked[sel].sum(),
            s.count[sel].sum(),
        )
    single = plain & (s.count == 1.0)
    assert single.sum() > 100 and np.all(np.asarray(belief.f_g)[single] == 1.0)


@pytest.mark.parametrize("seed", [3, 11, 29])
def test_a_library_with_no_overhang_is_the_lock(monkeypatch, seed):
    """Every locked object balanced: no strand refutes the lock, and every locked slot is all gDNA at no variance."""
    parts = _library(seed, overhang=False, strand=LIVE)
    belief, _ = _solve(monkeypatch, parts, LIVE)
    s = _slots(parts)
    held = s.locked & (s.count > 0)
    assert np.all(np.asarray(belief.f_g)[held] > 0.999), np.asarray(belief.f_g)[held].min()


@pytest.mark.parametrize("strand", [DEAD, StrandBalance(0.6, 4, False)])
def test_a_dead_strand_channel_is_the_lock(monkeypatch, strand):
    """The overhang library read without a strand protocol: an unstranded one, and one whose κ of 0.6 from four
    spliced reads is a fluctuation (`strand_discriminability`). The solve reads no strand there, and neither does the
    lock: every locked slot is exactly all gDNA at no variance, and the sweep is handed the geometry itself."""
    parts = _library(11, overhang=True, strand=strand)
    belief, inputs = _solve(monkeypatch, parts, strand)
    s = _slots(parts)
    assert np.all(np.asarray(belief.f_g)[s.locked] == 1.0)
    assert np.all(np.asarray(belief.var_gdna)[s.locked] == 0.0)
    assert inputs.lock.departure == 0.0 and inputs.message_geometry is inputs.geometry


def test_a_departed_lock_floors_nothing_and_still_books_what_its_strands_share(monkeypatch):
    """The edge rule's floor is a structurally pure count's. A boundary whose strands proved it departed is no longer
    pure: what it shares is an estimate with a width, as any single-strand slot's gDNA is, and no such slot floors its
    neighbour. So the sweep is handed the part of its count that is still a lock, under half a fragment, while the
    boundary itself books its shared part, a fragment or more. A lock that stands hands the sweep its whole count."""
    parts = _library(11, overhang=True, strand=LIVE)
    belief, inputs = _solve(monkeypatch, parts, LIVE)
    s = _slots(parts)
    handed = np.asarray(inputs.message_geometry.unspliced_count, np.float64).sum(axis=1)
    booked = np.asarray(belief.f_g) * s.count
    marked = s.right[MARKED]
    assert handed[marked].mean() < 0.5, handed[marked].mean()
    assert booked[marked].mean() > 1.0, booked[marked].mean()
    stands = s.locked & (np.asarray(belief.f_g) == 1.0)
    assert stands.sum() > 500 and np.array_equal(handed[stands], s.count[stands])
    assert np.array_equal(handed[~s.locked], s.count[~s.locked])


def test_the_backgrounds_law_reads_what_the_locks_share(monkeypatch):
    """The intergenic background is the locked regions' gDNA over their opportunity. With the overhang counted the
    pooled rate reads a sixth high; read at each region's share it is the library's rate within a twentieth."""
    parts = _library(11, overhang=True, strand=LIVE)
    _, inputs = _solve(monkeypatch, parts, LIVE)
    assert abs(inputs.background.pooled_rate / BACKGROUND - 1.0) < 0.05, (
        inputs.background.pooled_rate
    )
    s = _slots(parts)
    regions = s.locked & (np.asarray(parts.chain.kind) == REGION)
    counted = s.count[regions].sum() / (regions.sum() * (INTERGENIC - FRAGMENT + 1))
    assert counted / BACKGROUND > 1.1, counted
