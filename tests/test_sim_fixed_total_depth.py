"""Falsification gates for the FIXED-TOTAL depth mode (``simulation.n_total_fragments``).

⭐ **WHY THE MODE EXISTS.** The panel's gDNA axis is a *rate*: ``n_gdna = rate x n_rna``, on top of a
fixed RNA depth. That cannot reach the high-gDNA end of the real spectrum — a 98 % gDNA library is
``rate = 49``, i.e. 490 M fragments at a 10 M RNA depth. Real libraries do not work that way either: a
sequencing run has a fixed **total** budget and the gDNA fraction decides how it is *split*, so holding
the total and varying the split is both feasible and the more faithful model. Owner ruling, 2026-08-03:
10 M total per condition is the ceiling, and the RNA-side accuracy loss at high gDNA is accepted.

⛔ Each gate below carries its own perturbation, because a gate that has never been watched to fire has
not been written yet.
"""

from __future__ import annotations

import pytest

from rigel.sim.orchestrator import resolve_depths
from rigel.sim.wgs_config import SimulationParams

#: The ladder the panel uses. ⚠ Rates, not fractions: ``f_gdna = rate / (1 + rate)``.
LADDER = [0.0, 0.010101, 0.052632, 0.111111, 0.333333, 1.0, 3.0, 9.0, 49.0]
TOTAL = 10_000_000


def _params(**kw) -> SimulationParams:
    return SimulationParams(n_rna_fragments=10_000_000, **kw)


def test_the_total_is_CONSERVED_EXACTLY_at_every_rung():
    """``n_rna + n_gdna`` must equal the requested total on every rung, exactly — not to within
    rounding drift. A ladder whose rungs have different depths cannot be compared across rungs, which
    is the entire purpose of fixing the total.

    PERTURBATION: the legacy path (``n_total_fragments=None``) must NOT conserve it — otherwise this
    gate is passing on a property the old code already had and proves nothing about the new mode.
    """
    sim = _params(n_total_fragments=TOTAL)
    for rate in LADDER:
        d = resolve_depths(sim, gdna_rate=rate, nrna_ratio=None, nrna_mode="file")
        assert d.n_rna + d.n_gdna == TOTAL, f"rate {rate}: {d.n_rna} + {d.n_gdna} != {TOTAL}"

    legacy = _params(n_total_fragments=None)
    totals = {
        resolve_depths(legacy, gdna_rate=r, nrna_ratio=None, nrna_mode="file").total for r in LADDER
    }
    assert len(totals) > 1, "the legacy path already fixes the total; the new mode adds nothing"


def test_the_realised_gDNA_FRACTION_tracks_the_rung():
    """The point of the ladder is to sweep ``f_gdna``, so the realised fraction must be the requested
    one. ⚠ Checked against ``rate/(1+rate)``, which is what the rate MEANS — not against the label.

    PERTURBATION: a ladder that produced a constant fraction would still conserve the total, so the
    gate above cannot catch it. This one requires the fractions to span the range they claim.
    """
    sim = _params(n_total_fragments=TOTAL)
    fracs = []
    for rate in LADDER:
        d = resolve_depths(sim, gdna_rate=rate, nrna_ratio=None, nrna_mode="file")
        want = rate / (1.0 + rate)
        got = d.n_gdna / d.total
        assert abs(got - want) <= 1.0 / TOTAL, f"rate {rate}: f_gdna {got:.6f} != {want:.6f}"
        fracs.append(got)
    assert fracs[0] == 0.0, "the zero rung must be EXACTLY zero — it is the false-positive control"
    assert fracs[-1] > 0.97, "the top rung must actually reach the high-gDNA regime"
    assert fracs == sorted(fracs), "the ladder must be monotone"


def test_the_nascent_split_is_INSIDE_the_RNA_budget_not_on_top():
    """With nascent enabled the total must still hold: nRNA comes out of the RNA share, it is not
    added on top. Otherwise the nascent axis silently changes the depth and every cross-condition
    comparison acquires a confound.

    PERTURBATION: the same call with ``nrna_ratio=0`` must give a strictly larger ``n_mrna`` — if the
    ratio is being ignored, this gate is vacuous.
    """
    sim = _params(n_total_fragments=TOTAL)
    d = resolve_depths(sim, gdna_rate=1.0, nrna_ratio=0.25, nrna_mode="additive_ratio")
    assert d.n_mrna + d.n_nrna == d.n_rna
    assert d.n_rna + d.n_gdna == TOTAL
    assert d.n_nrna > 0

    none = resolve_depths(sim, gdna_rate=1.0, nrna_ratio=0.0, nrna_mode="additive_ratio")
    assert none.n_nrna == 0 and none.n_mrna > d.n_mrna
    assert none.n_rna + none.n_gdna == TOTAL


def test_omitting_the_total_leaves_the_LEGACY_path_byte_identical():
    """⛔ The mode must be opt-in. Every existing config omits ``n_total_fragments``, and those panels
    must simulate exactly what they simulated before — a depth change would silently invalidate every
    stored number measured against them.

    PERTURBATION: switching the mode on must move the numbers, or "byte-identical" is vacuous.
    """
    legacy = _params(n_total_fragments=None)
    for rate in LADDER:
        d = resolve_depths(legacy, gdna_rate=rate, nrna_ratio=None, nrna_mode="file")
        assert d.n_rna == 10_000_000
        assert d.n_gdna == round(rate * 10_000_000)

    on = resolve_depths(
        _params(n_total_fragments=TOTAL), gdna_rate=1.0, nrna_ratio=None, nrna_mode="file"
    )
    off = resolve_depths(legacy, gdna_rate=1.0, nrna_ratio=None, nrna_mode="file")
    assert (on.n_rna, on.n_gdna) != (off.n_rna, off.n_gdna)


def test_a_total_too_small_for_the_rung_RAISES_rather_than_rounding_to_zero():
    """⛔ At rate 49 a 100-fragment total leaves 2 RNA fragments; at some point the RNA side rounds to
    zero and the condition silently stops being an RNA-seq library at all. That must fail loudly —
    a mis-stated depth must raise, never quietly produce nothing.

    PERTURBATION: a total that IS large enough must not raise.
    """
    with pytest.raises(ValueError, match="RNA"):
        resolve_depths(
            _params(n_total_fragments=10), gdna_rate=49.0, nrna_ratio=None, nrna_mode="file"
        )
    ok = resolve_depths(
        _params(n_total_fragments=TOTAL), gdna_rate=49.0, nrna_ratio=None, nrna_mode="file"
    )
    assert ok.n_rna == 200_000
