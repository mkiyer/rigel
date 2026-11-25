"""The intron abundance model: an intron's unspliced count read against the intergenic background.

An intergenic region admits no RNA strand, so the intergenic regions are the library's one structurally pure gDNA pool,
and `gdna_density.background_law` reads their rate and spread. An intron's gDNA is one more draw from that law at the
intron's own gDNA opportunity, so on an unstranded library, where the strand channel is dead, the background is the
intron's only composition evidence: its count's excess over the background's expectation is RNA. With prior
probability ε (the run setting ``calibration.intron_departure_prior``) a region is not at the background (a probe's
reach, a copy-number change, a mappability gap, unannotated transcription), so the row's tail is capped and an intron's
own decisive strand evidence outweighs a background it contradicts, and a region that departs does not set the law.
Under the transfer policy the row is the intron's own claim.
"""

from __future__ import annotations

import dataclasses
import sys

import numpy as np
import pytest

from _synthetic import make_chain_parts
from rigel._solve_impl import testing as R
from rigel.calibration.diagnostics import SweepDiagnostics
from rigel.calibration.effective_length import contained_eff_length
from rigel.calibration.gdna_density import background_law
from rigel.calibration.messages.silent import SilentPolicy
from rigel.calibration.messages.transfer import TransferPolicy
from rigel.calibration.region_geometry import init_beliefs
from rigel.calibration.signature import BIT_EXON_POS, BIT_INTRON_POS
from rigel.calibration.strand_balance import StrandBalance
from rigel.calibration.sweep import solve_chain
from rigel.config import CalibrationConfig

CAL = sys.modules["rigel.calibration.calibrate"]
EPS = CalibrationConfig().intron_departure_prior

L = 2000.0  # every region's length
N_INTERGENIC = 6  # intergenic regions on each side of the gene
# the chain: intergenic × 6 | exon+ | intron+ | exon+ | intergenic × 6
SIGNATURES = [0] * N_INTERGENIC + [BIT_EXON_POS, BIT_INTRON_POS, BIT_EXON_POS] + [0] * N_INTERGENIC
INTRON = N_INTERGENIC + 1  # the intron's region index; its chain slot is twice that
INTRON_SLOT = 2 * INTRON


def _delta_pmf(length):
    p = np.zeros(length + 1)
    p[length] = 1.0
    return p


def _gene(*, rho_g, intron_gdna, intron_rna, kappa, boundary_rna=0.0):
    """The chain with every intergenic region exactly at gDNA density ``rho_g`` and the intron holding
    ``intron_gdna + intron_rna`` fragments: gDNA split evenly between the strands, RNA on its sense strand at ``kappa``
    (the + strand's sense share; ½ is unstranded). The exons hold background gDNA and ten times that in RNA."""
    gdna_fl, rna_fl = _delta_pmf(200), _delta_pmf(150)
    eg = float(contained_eff_length(np.full(1, L), gdna_fl)[0])
    g = rho_g * eg
    pos = np.full(len(SIGNATURES), g / 2.0)
    neg = np.full(len(SIGNATURES), g / 2.0)
    for e in (INTRON - 1, INTRON + 1):
        pos[e] += 10.0 * g * kappa
        neg[e] += 10.0 * g * (1.0 - kappa)
    pos[INTRON] = intron_gdna / 2.0 + intron_rna * kappa
    neg[INTRON] = intron_gdna / 2.0 + intron_rna * (1.0 - kappa)
    parts = make_chain_parts(
        SIGNATURES,
        region_size_bp=L,
        region_pos=pos,
        region_neg=neg,
        boundary_pos=1.0 + boundary_rna * kappa,
        boundary_neg=1.0 + boundary_rna * (1.0 - kappa),
        gdna_fl=gdna_fl,
        rna_fl=rna_fl,
    )
    return parts, eg


def _solve(parts, *, kappa, background, policy=None):
    belief = init_beliefs(
        parts.geometry, parts.statics, rna_sense_frac=kappa, n_grid=101, logodds_window=10.0
    )
    cap = SweepDiagnostics()
    out = solve_chain(
        parts.chain,
        parts.statics,
        parts.geometry,
        belief,
        parts.region_arrays,
        strand=StrandBalance(kappa, 10_000, False),
        n_grid=101,
        logodds_window=10.0,
        background=background,
        policy=SilentPolicy() if policy is None else policy,
        diagnostics=cap,
    )
    return out, cap


def _law(parts, departure=EPS):
    return CAL._background(parts.chain, parts.statics, parts.geometry, departure)


def _scattered_pool(rng, *, n, exposure, rate, alpha):
    """``n`` regions of one exposure whose gDNA rate scatters as ``Gamma(α, rate/α)``: a negative binomial pool."""
    lam = rng.gamma(alpha, rate * exposure / alpha, size=n)
    return rng.poisson(lam).astype(float), np.full(n, exposure)


def test_a_clean_pool_gives_its_rate_and_its_spread():
    """A Poisson pool reads its rate and a size that is the rate's own posterior width alone; an over-dispersed pool's
    spread joins it; a pool with no opportunity is not informative."""
    rng = np.random.default_rng(7)
    e = np.full(400, 1800.0)
    poisson = background_law(rng.poisson(90.0, size=400).astype(float), e, EPS)
    assert np.exp(poisson.log_rate) == pytest.approx(0.05, rel=0.01)
    assert poisson.size > 1_000.0  # no excess spread: the posterior's width alone
    counts, e = _scattered_pool(rng, n=400, exposure=1800.0, rate=0.05, alpha=20.0)
    spread = background_law(counts, e, EPS)
    assert np.exp(spread.log_rate) == pytest.approx(0.05, rel=0.03)
    assert 10.0 < spread.size < 40.0, spread.size  # α = 20 recovered
    assert not background_law(np.zeros(3), np.zeros(3), EPS).informative


def test_one_region_that_departs_does_not_set_the_background():
    """220 regions at the background and one long region that also holds unannotated transcription: its RNA nearly
    doubles its count. Pooled, it would raise the rate by two thirds and widen the law to a size of a few; read so that a
    region the law cannot explain does not set it, the background is the clean regions' own."""
    rng = np.random.default_rng(11)
    rate = 0.006
    e = np.full(220, 2000.0)
    n = rng.poisson(rate * e).astype(float)
    big = 1.5e6
    n = np.append(n, rng.poisson(rate * big) + 8000.0)
    e = np.append(e, big)
    law = background_law(n, e, EPS)
    assert np.exp(law.log_rate) == pytest.approx(rate, rel=0.02)
    assert law.size > 1_000.0, law.size
    assert (
        law.pooled_rate > 1.5 * rate
    )  # the naive rate, reported beside it, is what the region would have set
    # the clean regions keep nearly all their weight (at a dozen fragments the flat departure is not negligible), the
    # departed one none
    assert 210.0 < law.kept < 220.5, law.kept


def test_an_unstranded_intron_reads_its_gdna_against_the_background():
    """Unstranded, the intron's strand channel is dead; its count is the background's gDNA plus four times as much RNA,
    so it is one-fifth gDNA. Read against the background it says so, and its precision is its own evidence; without the
    background it has no evidence and stays at the reference."""
    parts, eg = _gene(
        rho_g=0.05, intron_gdna=0.05 * 1801.0, intron_rna=4 * 0.05 * 1801.0, kappa=0.5
    )
    law = _law(parts)
    assert law.informative and law.n_objects == 2 * N_INTERGENIC
    with_model, cap = _solve(parts, kappa=0.5, background=law)
    without, _ = _solve(parts, kappa=0.5, background=None)
    fg = float(with_model.f_g[INTRON_SLOT])
    assert fg == pytest.approx(0.2, abs=0.03), fg
    assert float(np.asarray(cap.tau_lam)[INTRON_SLOT]) > float(R.OWN_EVIDENCE_EPS)
    assert abs(float(without.f_g[INTRON_SLOT]) - 0.5) < 0.1, without.f_g[INTRON_SLOT]


def test_a_peak_between_lattice_nodes_is_kept():
    """A deep unstranded intron, 30,000 fragments over a background expecting 6,420, so 0.214 gDNA. The background's
    peak is narrower than a lattice cell and falls between two nodes, both far enough down its flanks to sit under the
    departure floor; read at the nodes alone, the row would be flat and the intron would read ½. Each node stands for
    its cell, so the background is read at the cell's best point and the intron reads a node beside the truth."""
    mu, count = 6420.0, 30000.0
    eg = float(contained_eff_length(np.full(1, L), _delta_pmf(200))[0])
    parts, _ = _gene(rho_g=mu / eg, intron_gdna=mu, intron_rna=count - mu, kappa=0.5)
    out, _ = _solve(parts, kappa=0.5, background=_law(parts))
    step = 20.0 / 100  # the lattice of `_solve`: 101 nodes over ±10
    truth = np.log(mu / (count - mu))
    lo, hi = 1.0 / (1.0 + np.exp(-(truth - step))), 1.0 / (1.0 + np.exp(-(truth + step)))
    assert lo <= float(out.f_g[INTRON_SLOT]) <= hi, (out.f_g[INTRON_SLOT], lo, hi)


def test_decisive_strand_evidence_outweighs_a_background_it_contradicts():
    """A stranded intron holding 25 times the background's gDNA and no RNA (its strands balanced): the background alone
    calls most of it RNA, but the intron's own strand evidence is decisive, and the departure prior caps how far the
    background can argue against it. With no room for departure the background's Poisson tail would win."""
    rho_g = 2.3 / 1801.0
    parts, _ = _gene(rho_g=rho_g, intron_gdna=25 * 2.3, intron_rna=0.0, kappa=0.99)
    out, _ = _solve(parts, kappa=0.99, background=_law(parts))
    assert float(out.f_g[INTRON_SLOT]) > 0.9, out.f_g[INTRON_SLOT]


def test_no_departure_room_lets_the_background_tail_win():
    """The departure prior is what lets the strand evidence decide: with it near zero, the same intron is read mostly
    RNA. A guard on the honest-tail test, so the cap cannot pass by accident."""
    rho_g = 2.3 / 1801.0
    parts, _ = _gene(rho_g=rho_g, intron_gdna=25 * 2.3, intron_rna=0.0, kappa=0.99)
    law = dataclasses.replace(_law(parts), departure=1e-300)
    out, _ = _solve(parts, kappa=0.99, background=law)
    assert float(out.f_g[INTRON_SLOT]) < 0.5, out.f_g[INTRON_SLOT]


def test_a_row_the_transfer_reads_as_flat_is_no_evidence_to_the_solve():
    """One flatness: a row whose spread is within the transfer's flat threshold sends no claim, so it carries no
    composition evidence to the solve either (the precision read a thousand-fold narrower threshold before); a row
    with a shape does."""
    lam = np.linspace(-5.0, 5.0, 11)
    flat = float(R.FLAT_ROW_SPREAD)
    assert R.row_precision(0.4 * flat * lam / 5.0, lam) == 0.0
    assert R.row_precision(-(lam**2), lam) > 0.0


def test_the_solve_refuses_a_background_with_no_finite_size():
    """An informative law's size is finite (``background_law`` joins the rate's own posterior width to its
    over-dispersion), so the solve reads no Poisson limit: an infinite size is refused, never read as one."""
    parts, _ = _gene(rho_g=0.05, intron_gdna=0.05 * 1801.0, intron_rna=0.0, kappa=0.5)
    law = dataclasses.replace(_law(parts), size=np.inf)
    with pytest.raises(ValueError, match="size is finite"):
        _solve(parts, kappa=0.5, background=law)


def test_the_intron_claims_its_row_under_the_transfer_policy():
    """Under the transfer policy the intron's row is its own claim, so its exon|intron boundaries, which share its
    densities and carry nothing of their own on an unstranded library, read the intron's composition."""
    parts, _ = _gene(
        rho_g=0.05,
        intron_gdna=0.05 * 1801.0,
        intron_rna=4 * 0.05 * 1801.0,
        kappa=0.5,
        boundary_rna=0.0,
    )
    law = _law(parts)
    policy = TransferPolicy(kappa=0.5)
    with_model, _ = _solve(parts, kappa=0.5, background=law, policy=policy)
    without, _ = _solve(parts, kappa=0.5, background=None, policy=policy)
    flanks = [INTRON_SLOT - 1, INTRON_SLOT + 1]
    moved = np.abs(np.asarray(with_model.f_g)[flanks] - np.asarray(without.f_g)[flanks])
    assert np.all(moved > 0.05), (with_model.f_g[flanks], without.f_g[flanks])


def test_the_departure_prior_is_a_run_setting():
    """``calibration.intron_departure_prior`` is the law's departure prior, and only a probability strictly between 0
    and 1 is accepted: at 0 the tail is gone, at 1 the background says nothing."""
    parts, _ = _gene(rho_g=0.05, intron_gdna=90.0, intron_rna=0.0, kappa=0.5)
    assert _law(parts, 0.25).departure == 0.25
    for bad in (0.0, 1.0, -0.1, 1.5):
        with pytest.raises(ValueError, match="intron_departure_prior"):
            CalibrationConfig(intron_departure_prior=bad)
