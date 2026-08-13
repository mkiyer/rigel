"""The four gDNA length pools' opportunity functions, against a BRUTE-FORCE ENUMERATING ORACLE.

    `rigel.calibration.gdna_opportunity`

⭐ **The oracle enumerates; the code computes.** Every closed form here is checked against a loop that
walks every start position on a small partition and asks the *deposit rule's own question* — "is this
fragment contained in exactly one region?", "does it cross exactly one boundary?" — so the divisor is verified
against the selection it is meant to invert, not against a re-derivation of itself.
`docs/TRAPS.md` trap 1: a validator that calls the builder's own helper validates nothing.

⚠ **The rule the oracle implements is `_accumulator_reference.FragmentPool`'s, verbatim:**

* contained in exactly one region -> typed by that region's coarse type (intergenic / intronic only)
* crossing exactly one boundary -> typed by the sorted pair of flanking region types
* anything else — an exonic contained fragment, a multi-boundary crossing — enters **no** pool
"""

from __future__ import annotations

import numpy as np
import pytest

from rigel.calibration.gdna_opportunity import (
    GdnaOpportunity,
    contained_opportunity,
    crossing_opportunity,
    total_opportunity,
)

TYPE_INTERGENIC, TYPE_INTRON, TYPE_EXON = 0, 1, 2


# ── the oracle ──────────────────────────────────────────────────────────────────────────────────


def enumerate_pools(cuts: list[int], types: list[int], width: int) -> dict:
    """Walk every start position on one reference and classify it exactly as the deposit rule does.

    ``cuts`` is the reference's cut axis, ascending, first 0 and last ``L_ref``. ``types`` is one coarse
    type per region, so ``len(types) == len(cuts) - 1``.
    """
    reference_length = cuts[-1]
    boundaries = cuts[1:-1]  # interior boundaries only
    counts = {
        "contained": {TYPE_INTERGENIC: 0, TYPE_INTRON: 0, TYPE_EXON: 0},
        "crossing": {},
        "total": 0,
    }
    for start in range(0, reference_length - width + 1):
        end = start + width
        counts["total"] += 1
        crossed = [i for i, p in enumerate(boundaries) if start < p < end]
        if not crossed:
            # Contained: find the single region holding [start, end).
            region = max(i for i, c in enumerate(cuts[:-1]) if c <= start)
            counts["contained"][types[region]] += 1
        elif len(crossed) == 1:
            boundary = crossed[0]
            pair = tuple(sorted((types[boundary], types[boundary + 1])))
            counts["crossing"][pair] = counts["crossing"].get(pair, 0) + 1
    return counts


def region_lengths_of(cuts: list[int]) -> np.ndarray:
    return np.diff(np.asarray(cuts, dtype=np.int64))


# ── the partitions the oracle is run over ───────────────────────────────────────────────────────

#: Deliberately awkward: a 1 bp region (legal, and 15,687 exist in the human index), a region shorter than
#: the widest fragment, adjacent same-type regions, and both flank pairs present.
PARTITIONS = [
    ([0, 40, 41, 90, 200], [TYPE_INTERGENIC, TYPE_EXON, TYPE_INTRON, TYPE_EXON]),
    ([0, 10, 60, 61, 62, 130], [TYPE_INTRON, TYPE_EXON, TYPE_INTRON, TYPE_EXON, TYPE_INTERGENIC]),
    ([0, 100], [TYPE_INTERGENIC]),
    (
        [0, 5, 10, 15, 20, 25, 30],
        [TYPE_INTERGENIC, TYPE_EXON, TYPE_INTRON, TYPE_EXON, TYPE_INTRON, TYPE_EXON],
    ),
    ([0, 3, 250], [TYPE_EXON, TYPE_INTRON]),
]
WIDTHS = [1, 2, 3, 7, 20, 41, 55, 99, 130]


class TestAgainstTheEnumeratingOracle:
    @pytest.mark.parametrize("cuts,types", PARTITIONS)
    @pytest.mark.parametrize("width", WIDTHS)
    def test_contained_opportunity_matches_enumeration(self, cuts, types, width):
        oracle = enumerate_pools(cuts, types, width)
        lengths = region_lengths_of(cuts)
        max_width = max(WIDTHS)
        for region_type in (TYPE_INTERGENIC, TYPE_INTRON, TYPE_EXON):
            selected = lengths[np.asarray(types) == region_type]
            computed = contained_opportunity(selected, max_width)[width]
            assert computed == pytest.approx(oracle["contained"][region_type]), (
                f"type {region_type}, width {width}, cuts {cuts}"
            )

    @pytest.mark.parametrize("cuts,types", PARTITIONS)
    @pytest.mark.parametrize("width", WIDTHS)
    def test_crossing_opportunity_matches_enumeration(self, cuts, types, width):
        oracle = enumerate_pools(cuts, types, width)
        lengths = region_lengths_of(cuts)
        max_width = max(WIDTHS)
        if len(lengths) < 2:
            return
        left, right = lengths[:-1], lengths[1:]
        pairs = np.sort(np.stack([np.asarray(types)[:-1], np.asarray(types)[1:]]), axis=0)
        for pair in {tuple(sorted(p)) for p in zip(types, types[1:])}:
            mask = (pairs[0] == pair[0]) & (pairs[1] == pair[1])
            computed = crossing_opportunity(left[mask], right[mask], max_width)[width]
            assert computed == pytest.approx(oracle["crossing"].get(pair, 0)), (
                f"pair {pair}, width {width}, cuts {cuts}"
            )

    @pytest.mark.parametrize("cuts,types", PARTITIONS)
    @pytest.mark.parametrize("width", WIDTHS)
    def test_total_opportunity_matches_enumeration(self, cuts, types, width):
        oracle = enumerate_pools(cuts, types, width)
        computed = total_opportunity(np.asarray([cuts[-1]]), max(WIDTHS))[width]
        assert computed == pytest.approx(oracle["total"])

    @pytest.mark.parametrize("cuts,types", PARTITIONS)
    @pytest.mark.parametrize("width", WIDTHS)
    def test_every_start_is_accounted_for(self, cuts, types, width):
        """⭐ The pools plus the unpooled remainder must exhaust every admissible start.

        This is the one assertion that would catch an opportunity function that is individually
        plausible and collectively wrong — double-counting a start, or losing one at a boundary.
        """
        oracle = enumerate_pools(cuts, types, width)
        pooled = sum(oracle["contained"].values()) + sum(oracle["crossing"].values())
        assert pooled <= oracle["total"]
        multi_boundary = oracle["total"] - pooled
        assert multi_boundary >= 0
        # And the closed forms sum to the same pooled figure.
        lengths = region_lengths_of(cuts)
        max_width = max(WIDTHS)
        computed = sum(
            contained_opportunity(lengths[np.asarray(types) == t], max_width)[width]
            for t in (TYPE_INTERGENIC, TYPE_INTRON, TYPE_EXON)
        )
        if len(lengths) >= 2:
            left, right = lengths[:-1], lengths[1:]
            pairs = np.sort(np.stack([np.asarray(types)[:-1], np.asarray(types)[1:]]), axis=0)
            for pair in {tuple(sorted(p)) for p in zip(types, types[1:])}:
                mask = (pairs[0] == pair[0]) & (pairs[1] == pair[1])
                computed += crossing_opportunity(left[mask], right[mask], max_width)[width]
        assert computed == pytest.approx(pooled)


class TestTheTiltsAreOpposite:
    """⛔ The reason the four pools cannot be pooled raw, asserted rather than asserted-in-prose."""

    def test_contained_opportunity_falls_with_length_and_crossing_rises(self):
        lengths = np.array([500, 1200, 80, 3000], dtype=np.int64)
        contained = contained_opportunity(lengths, 400)
        # Strictly decreasing over the live range: every extra base of fragment removes one start.
        assert np.all(np.diff(contained[1:301]) < 0)

        left = np.array([500, 700, 900], dtype=np.int64)
        right = np.array([600, 800, 1000], dtype=np.int64)
        crossing = crossing_opportunity(left, right, 400)
        # Strictly increasing: a longer fragment covers more of the boundary's neighbourhood.
        assert np.all(np.diff(crossing[2:301]) > 0)

    def test_a_short_region_caps_the_crossing_opportunity(self):
        """⭐ Where a flank is shorter than the fragment, the crossing count STOPS growing — that is the
        `(w-1-a)+` term, and it is why `(w-1)` alone is not the opportunity."""
        crossing = crossing_opportunity(
            np.array([30], dtype=np.int64), np.array([40], dtype=np.int64), 300
        )
        # Beyond a + b = 70 the fragment always crosses a neighbour too, so it leaves the pool.
        assert crossing[200] == 0.0
        assert crossing[60] > 0.0


class TestTheCombination:
    def test_combined_probability_is_the_opportunity_weighted_average(self):
        """⭐ `(sum count)/(sum A)` IS `sum(A_p f_p)/sum(A_p)` — the identity the design rests on."""
        rng = np.random.default_rng(0)
        max_width = 60
        pools = [rng.random(max_width + 1) * 1000 + 1 for _ in range(4)]
        total = np.full(max_width + 1, 5000.0)
        opportunity = GdnaOpportunity(*pools, total)
        counts = [rng.random(max_width + 1) * 50 for _ in range(4)]

        summed = np.sum(counts, axis=0) / opportunity.combined_probability()
        per_pool = [c / (a / total) for c, a in zip(counts, pools)]
        weighted = np.sum([a * f for a, f in zip(pools, per_pool)], axis=0) / np.sum(pools, axis=0)
        assert summed == pytest.approx(weighted)

    def test_an_empty_annotation_yields_an_inert_correction(self):
        zero = np.zeros(11)
        opportunity = GdnaOpportunity(zero, zero, zero, zero, zero)
        assert np.all(opportunity.combined_probability() == 0.0)


class TestWiredIntoTheModel:
    def test_the_fallback_is_the_CONTAINED_pair_not_the_raw_four(self):
        """⛔ Without a divisor the gDNA pool must fall back to the two contained pools.

        Pooling the four raw is measurably WORSE than either the contained pair or the de-tilted four,
        so "no annotation offered" must not silently pick it.
        """
        from types import SimpleNamespace

        from rigel.calibration.fl import build_fl_models, gdna_contained_fl_mass, gdna_fl_mass

        max_size = 20
        pools = np.zeros((5, max_size + 1))
        pools[0, 5] = 100.0  # intergenic contained
        pools[1, 6] = 50.0  # intronic contained
        pools[2, 18] = 400.0  # intron-exon crossing — much longer, and much larger
        pools[3, 19] = 300.0  # intergenic-exon crossing
        pools[4, 10] = 70.0  # RNA
        payload = SimpleNamespace(
            pool_lengths=pools,
            deposited_lengths=pools.sum(axis=0),
            max_length=max_size,
        )

        assert float(gdna_contained_fl_mass(payload).sum()) == 150.0
        assert float(gdna_fl_mass(payload).sum()) == 850.0

        models = build_fl_models(payload, prior_ess=0.0)
        # The fallback must be the contained pair — 150 fragments at 5 and 6 bp, not 850 dragged long.
        assert models.n_gdna == pytest.approx(150.0)
        assert float((models.gdna_pmf * np.arange(max_size + 1)).sum()) < 7.0

    def test_the_divisor_moves_the_gdna_pool_and_preserves_its_evidence_weight(self):
        """⭐ The four-pool arm uses all 850 fragments, and the EB weight reflects that.

        ⚠ **The RNA pool here spans SEVERAL lengths on purpose.** `detilt_pool` renormalises back to the
        pool's own total, so a **single-bin** histogram is invariant under any divisor whatsoever — and a
        one-bin RNA pool made this test blind to the gDNA divisor leaking into the RNA one. Found by
        perturbation, not review.
        """
        from types import SimpleNamespace

        from rigel.calibration.fl import build_fl_models
        from rigel.calibration.gdna_opportunity import GdnaOpportunity

        max_size = 20
        pools = np.zeros((5, max_size + 1))
        pools[0, 5] = 100.0
        pools[1, 6] = 50.0
        pools[2, 18] = 400.0
        pools[3, 19] = 300.0
        pools[4, 8] = 30.0
        pools[4, 12] = 25.0
        pools[4, 17] = 15.0
        payload = SimpleNamespace(
            pool_lengths=pools, deposited_lengths=pools.sum(axis=0), max_length=max_size
        )
        width = np.arange(max_size + 1, dtype=np.float64)
        opportunity = GdnaOpportunity(
            intergenic_contained=np.maximum(1000.0 - width, 0.0),
            intronic_contained=np.maximum(800.0 - width, 0.0),
            intron_exon_crossing=np.maximum(width - 1.0, 0.0) * 10.0,
            intergenic_exon_crossing=np.maximum(width - 1.0, 0.0) * 5.0,
            total=np.full(max_size + 1, 5000.0),
        )

        plain = build_fl_models(payload, prior_ess=0.0)
        detilted = build_fl_models(payload, gdna_opportunity=opportunity, prior_ess=0.0)

        assert detilted.n_gdna == pytest.approx(850.0), (
            "the crossing pools must contribute evidence"
        )
        assert not np.allclose(plain.gdna_pmf, detilted.gdna_pmf)
        # The RNA pool is untouched by the gDNA divisor — the two selections are different objects.
        assert plain.rna_pmf == pytest.approx(detilted.rna_pmf)

    def test_EVERY_production_caller_of_build_fl_models_passes_THE_GDNA_DIVISOR(self):
        """⛔ Same gate as the junction divisor's, for the same reason: optional means silently absent.

        ⚠ Source-level on purpose — a runtime check would need a full pipeline run per call site, and
        the failure this guards against is somebody adding a fourth caller.
        """
        import ast
        import inspect

        from rigel import pipeline, scan_cache

        for module in (pipeline, scan_cache):
            tree = ast.parse(inspect.getsource(module))
            calls = [
                region
                for region in ast.walk(tree)
                if isinstance(region, ast.Call)
                and isinstance(region.func, ast.Name)
                and region.func.id == "build_fl_models"
            ]
            assert calls, f"{module.__name__} no longer calls build_fl_models — retarget this test"
            for call in calls:
                assert any(k.arg == "gdna_opportunity" for k in call.keywords), (
                    f"{module.__name__}:{call.lineno} builds FL models without the gDNA divisor, so "
                    "the gDNA pool falls back to the contained pair and reads ~15 % short under capture"
                )


class TestFromTheIndex:
    def test_the_four_pools_are_built_off_the_index_partition(self, mini_index):
        """⭐ The index path, gated on the one invariant that must hold on any real partition:

        the four pools are disjoint selections of the same start positions, so their sum can never
        exceed the total number of admissible starts. A divisor that violates that is over-counting.
        """
        from rigel.calibration.gdna_opportunity import gdna_opportunity_from_index

        opportunity = gdna_opportunity_from_index(mini_index, 200)
        for pool in opportunity.pools:
            assert pool.shape == (201,)
            assert np.all(pool >= 0)
        assert np.all(np.sum(opportunity.pools, axis=0) <= opportunity.total + 1e-9)
        assert opportunity.total[0] > 0, "the fixture partition must host at least one placement"
        probability = opportunity.combined_probability()
        assert np.all((probability >= 0) & (probability <= 1 + 1e-12))
