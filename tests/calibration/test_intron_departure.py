"""The introns' departure prior ε, estimated per library from their antisense counts.

A stranded library measures how often an intron's gDNA departs from the intergenic background: on a single-strand
intron gDNA lands on either strand evenly while RNA lands on the antisense strand at 1 − κ, so the antisense count,
given the intron's count, reads its gDNA. An unstranded library cannot, and keeps the prior's mode, the run setting.

`calibrate` fits the intergenic pool at the run setting and hands the introns' rows the library's own ε.
"""

from __future__ import annotations

import dataclasses
import importlib

import numpy as np
import pytest
from _synthetic import (
    delta_pmf,
    make_chain_parts,
    make_strand_models,
    make_synthetic_payload,
    make_synthetic_sj,
    spike_pmf,
)

from rigel.calibration import calibrate
from rigel.calibration.diagnostics import CalibrationDiagnostics
from rigel.calibration.gdna_density import BackgroundLaw
from rigel.calibration.intron_departure import DepartureFit, fit_departure_prior
from rigel.calibration.messages.silent import SilentPolicy
from rigel.calibration.region_chain import REGION
from rigel.calibration.region_geometry import init_beliefs
from rigel.calibration.signature import (
    BIT_EXON_NEG,
    BIT_EXON_POS,
    BIT_INTRON_NEG,
    BIT_INTRON_POS,
)
from rigel.calibration.strand_balance import StrandBalance
from rigel.calibration.sweep import solve_chain
from rigel.config import CalibrationConfig

CAL = importlib.import_module("rigel.calibration.calibrate")
M0 = CalibrationConfig().intron_departure_prior
RATE = 0.01  # the background's gDNA per unit of opportunity


def _law(n_regions=200):
    """A Poisson background at ``RATE``, as the intergenic pool's fit would hand it over."""
    return BackgroundLaw(float(np.log(RATE)), np.inf, n_regions, 1e6, M0, RATE, float(n_regions))


def _truncated(rng, mu, top):
    """A draw from the background (Poisson at ``mu``) truncated to ``0…top``, by its inverse distribution."""
    from scipy.stats import poisson

    u = rng.random(mu.shape) * poisson.cdf(top, mu)
    return np.minimum(poisson.ppf(u, mu), top)


def _model_library(rng, *, n=2000, share=0.25, kappa=0.999):
    """``n`` single-strand introns drawn from the two-part law itself: each intron's count is its background plus RNA at
    0, 1 or 5 times it; given the count, its gDNA is the background truncated to the count, or, for a ``share`` of
    introns, uniform on 0…count; its antisense count is the strand likelihood's draw."""
    exposure = rng.uniform(500.0, 5000.0, n)
    mu = RATE * exposure
    count = rng.poisson(mu * (1.0 + rng.choice(np.asarray((0.0, 1.0, 5.0)), n))).astype(float)
    departs = rng.random(n) < share
    g = np.where(departs, np.floor(rng.random(n) * (count + 1.0)), _truncated(rng, mu, count))
    f = np.divide(g, count, out=np.zeros(n), where=count > 0)
    antisense = rng.binomial(count.astype(int), 0.5 * f + (1.0 - kappa) * (1.0 - f)).astype(float)
    return count - antisense, antisense, exposure, float(departs[count > 0].mean())


def _library(rng, *, n=2000, share=0.25, kappa=0.999, rna=(0.0, 1.0, 5.0)):
    """``n`` single-strand introns as a capture makes them: a ``share`` departing, their gDNA 5 to 50 times the
    background's (log-uniform), the rest's a draw from the background. Each carries RNA at a multiple of its background
    drawn from ``rna``. gDNA splits evenly between the strands, RNA lands antisense at ``1 − kappa``."""
    exposure = rng.uniform(500.0, 5000.0, n)
    mu = RATE * exposure
    departs = rng.random(n) < share
    fold = np.where(departs, np.exp(rng.uniform(np.log(5.0), np.log(50.0), n)), 1.0)
    g = rng.poisson(mu * fold)
    r = rng.poisson(mu * rng.choice(np.asarray(rna, float), n))
    antisense = rng.binomial(g, 0.5) + rng.binomial(r, 1.0 - kappa)
    return (
        (g + r - antisense).astype(float),
        antisense.astype(float),
        exposure,
        float(departs.mean()),
    )


def test_a_stranded_library_recovers_its_departing_share():
    """A quarter of the introns depart, drawn from the two-part law itself: at κ = 0.999 the estimate finds that share,
    decided by the library's own data rather than the prior."""
    rng = np.random.default_rng(1)
    sense, antisense, exposure, share = _model_library(rng)
    fit = fit_departure_prior(sense, antisense, exposure, _law(), 0.999, M0)
    assert fit.departure == pytest.approx(share, abs=0.03), (fit.departure, share)
    assert fit.evidence > 100.0, fit  # the strands prefer the estimate to the setting decisively


def test_capture_like_departures_are_not_underestimated():
    """Departures as a capture makes them, 5 to 50 times the background's gDNA, sit at one end of the departed law's
    span, so the fitted weight of that law, which must reach them, errs high. High is the safe side: a larger ε lets an
    intron's strand evidence decide sooner."""
    rng = np.random.default_rng(1)
    sense, antisense, exposure, share = _library(rng)
    fit = fit_departure_prior(sense, antisense, exposure, _law(), 0.999, M0)
    assert fit.departure >= share, (fit.departure, share)


def test_an_unstranded_library_keeps_the_prior_mode():
    """The same library at κ = ½: gDNA and RNA split alike, so the antisense count says nothing about departure, and
    the estimate is the prior's mode exactly."""
    rng = np.random.default_rng(1)
    sense, antisense, exposure, _ = _library(rng, kappa=0.5)
    fit = fit_departure_prior(sense, antisense, exposure, _law(), 0.5, M0)
    assert fit.departure == M0
    assert fit.evidence == 0.0


def test_rna_far_above_the_background_does_not_raise_the_estimate():
    """No intron departs, but every intron holds ten times its background in RNA, so its count sits far above the
    background. The estimate reads each intron's strand given its count, so the count says nothing, and the data, which
    show no departure, pull the estimate below the prior's mode."""
    rng = np.random.default_rng(2)
    sense, antisense, exposure, share = _library(rng, share=0.0, rna=(10.0,))
    assert share == 0.0
    fit = fit_departure_prior(sense, antisense, exposure, _law(), 0.999, M0)
    assert fit.departure < M0, fit
    assert fit.evidence > 10.0, fit


def test_counts_are_whole_fragments():
    """The estimate reads an intron's antisense count as a binomial draw, so a fractional count is refused rather than
    read as impossible at the background, which would count it as a certain departure."""
    with pytest.raises(ValueError, match="whole fragments"):
        fit_departure_prior([1.5], [1.0], [1000.0], _law(), 0.999, M0)


# --- the library's ε in calibrate -------------------------------------------------------------------------------------

LIVE = StrandBalance(
    0.99, 10_000, False
)  # a stranded library: read 1 on its transcript's strand 99 % of the time
LENGTHS = {"intergenic": 4000, "exon": 400, "intron": 2000}
SPARSE = 1e-4  # gDNA per placement of a 200 bp fragment: an intron's background expectation is 0.18 fragments
DEEP = 2.8e-3  # an intron's background expectation is 5 fragments


def _genes(rng, *, background, fold=1.0, rna=0.0, both=0, n_genes=1000, kappa=0.99):
    """``n_genes`` one-intron genes of alternating strand, each after an intergenic region, every count a whole
    fragment. Intergenic regions and exons hold gDNA at ``background`` per placement, the exons ten times as much RNA;
    each intron holds ``fold`` times the background's gDNA and ``rna`` times it in RNA. ``both`` more introns admit both
    strands and hold RNA on each at ten times the background. gDNA splits evenly between the strands, RNA lands on its
    transcript's strand at ``kappa``."""
    sig, size, pos, neg = [], [], [], []

    def region(signature, kind, gdna, rna_fold, strands):
        mean = background * (LENGTHS[kind] - 199)
        g = rng.poisson(gdna * mean)
        p = rng.binomial(g, 0.5)
        n = g - p
        for plus in strands:
            r = rng.poisson(rna_fold * mean)
            s = rng.binomial(r, kappa)
            p, n = (p + s, n + r - s) if plus else (p + r - s, n + s)
        sig.append(signature)
        size.append(LENGTHS[kind])
        pos.append(p)
        neg.append(n)

    for k in range(n_genes):
        plus = k % 2 == 0
        exon, intron = (BIT_EXON_POS, BIT_INTRON_POS) if plus else (BIT_EXON_NEG, BIT_INTRON_NEG)
        region(0, "intergenic", 1.0, 0.0, ())
        region(exon, "exon", 1.0, 10.0, (plus,))
        region(intron, "intron", fold, rna, (plus,))
        region(exon, "exon", 1.0, 10.0, (plus,))
    for _ in range(both):
        region(0, "intergenic", 1.0, 0.0, ())
        region(BIT_INTRON_POS | BIT_INTRON_NEG, "intron", 1.0, 10.0, (True, False))
    region(0, "intergenic", 1.0, 0.0, ())
    return make_chain_parts(
        sig,
        region_size_bp=np.asarray(size, float),
        region_pos=np.asarray(pos, float),
        region_neg=np.asarray(neg, float),
        gdna_fl=delta_pmf(200, 400),
        rna_fl=delta_pmf(150, 400),
    )


def _departure(parts, strand=LIVE):
    law = CAL._background(parts.chain, parts.geometry, parts.region_arrays, M0)
    fit = CAL._departure(
        parts.chain, parts.statics, parts.geometry, parts.region_arrays, strand, law, M0
    )
    return fit, law


def test_a_stranded_library_reads_its_own_departure():
    """Every intron holds ten times its background's gDNA, 50 fragments, and no RNA, its strands balanced. At κ = 0.99
    each intron's strands say its gDNA is not at the background, and a thousand introns decide the library's ε."""
    fit, _ = _departure(_genes(np.random.default_rng(3), background=DEEP, fold=10.0))
    assert fit.departure > 0.5, fit
    assert fit.evidence > 100.0, fit


def test_a_sparse_library_whose_introns_hold_gdna_reads_a_high_departure():
    """The introns hold their background's gDNA exactly, no departure, but the background expects 0.18 fragments, so
    every fragment an intron holds is an excess, which the row, charging nothing for RNA, reads as RNA. The strands say
    the fragments are gDNA, so the library's ε is high: what the estimate measures is how often the background misreads
    an intron, whatever the reason. It rises toward 1 as the library's introns outweigh the prior."""
    fit, _ = _departure(_genes(np.random.default_rng(3), background=SPARSE, n_genes=3000))
    assert fit.departure > 0.3, fit
    assert fit.evidence > 100.0, fit


def test_rna_on_either_strand_is_not_read_as_a_departure():
    """No intron departs. Each holds ten times its background in RNA on its gene's strand, half the genes on each
    strand, and a hundred more introns admit both strands and hold RNA on each. Read against its own annotated strand an
    intron's RNA is sense and its gDNA balanced, so the library shows no departure; an intron that admits both strands
    has no antisense to read, and is left out."""
    fit, _ = _departure(_genes(np.random.default_rng(4), background=DEEP / 2.5, rna=10.0, both=100))
    assert fit.departure < M0 / 2.0, fit  # moved by the data: the prior alone stays at M0


def test_a_dead_strand_channel_keeps_the_run_setting():
    """The departing library, but its spliced reads show no protocol: κ 0.6 from four observations is a fluctuation
    (`strand_discriminability`), so the solve reads no strand, and neither does the estimate. An unstranded library
    keeps the run setting too."""
    parts = _genes(np.random.default_rng(3), background=DEEP, fold=10.0)
    for strand in (StrandBalance(0.6, 4, False), StrandBalance(0.5, 10_000, False)):
        fit, _ = _departure(parts, strand)
        assert fit.departure == M0 and fit.evidence == 0.0, (strand, fit)


def test_the_introns_read_the_librarys_departure_and_the_pool_the_setting(monkeypatch):
    """`calibrate` fits the intergenic pool at the run setting, since intergenic regions depart at their own rate,
    which no strand measures, and then hands the introns' rows the library's own ε: the law the sweep reads carries the
    estimate and is otherwise the pool's fit at the setting."""
    estimate = DepartureFit(0.37, M0, 1, 12.0)
    monkeypatch.setattr(CAL, "_departure", lambda *args, **kwargs: estimate)
    payload, region_arrays = make_synthetic_payload()
    diagnostics = CalibrationDiagnostics()
    calibrate(
        payload=payload,
        region_arrays=region_arrays,
        strand_model=make_strand_models(0.95, 40),
        gdna_fl_pmf=spike_pmf(),
        rna_fl_pmf=spike_pmf(),
        config=CalibrationConfig(),
        sj=make_synthetic_sj(),
        boundary_flags=np.zeros(int(payload.n_boundaries), np.uint16),
        diagnostics=diagnostics,
    )
    assert diagnostics.departure is estimate
    assert diagnostics.background.departure == 0.37
    pool = CAL._background(diagnostics.chain, diagnostics.geometry, diagnostics.region_arrays, M0)
    assert dataclasses.replace(diagnostics.background, departure=M0) == pool


def _gdna_share(parts, law, strand=LIVE):
    """The single-strand introns' gDNA share, fragment-weighted, from a silent solve: each read on its own evidence."""
    kappa = float(strand.rna_sense_frac)
    belief = init_beliefs(
        parts.geometry, parts.statics, rna_sense_frac=kappa, n_grid=101, logodds_window=10.0
    )
    out = solve_chain(
        parts.chain,
        parts.statics,
        parts.geometry,
        belief,
        parts.region_arrays,
        rna_sense_frac=kappa,
        n_rna_obs=float(strand.n_observations),
        n_grid=101,
        logodds_window=10.0,
        background=law,
        policy=SilentPolicy(),
    )
    sig = np.asarray(parts.region_arrays.signature)
    single = np.isin(sig, (BIT_INTRON_POS, BIT_INTRON_NEG))
    obj = np.asarray(parts.chain.obj_idx)
    slots = np.flatnonzero((np.asarray(parts.chain.kind) == REGION) & single[obj])
    count = np.asarray(parts.geometry.unspliced_count, np.float64)[slots].sum(axis=1)
    return float((np.asarray(out.f_g)[slots] * count).sum() / count.sum())


def test_a_sparse_librarys_gdna_comes_back_once_its_strands_set_the_departure():
    """The sparse library: every intron holds its background's gDNA and no RNA. At the run setting the row reads each
    fragment as RNA and outweighs its few fragments' strands, so the introns read three-quarters RNA. At the library's
    own ε the row steps aside for the strands, and most of the gDNA comes back; a lone fragment's strand is weak
    evidence, so not all of it."""
    parts = _genes(np.random.default_rng(3), background=SPARSE, n_genes=3000)
    fit, law = _departure(parts)
    at_setting = _gdna_share(parts, law)
    own = _gdna_share(parts, dataclasses.replace(law, departure=fit.departure))
    assert at_setting < 0.3 and own > at_setting + 0.15, (at_setting, own, fit)
