"""Tests for rigel.stats — PipelineStats."""

from rigel.stats import PipelineStats


class TestPipelineStats:
    def test_defaults_are_zero(self):
        ps = PipelineStats()
        assert ps.total == 0
        assert ps.n_fragments == 0
        assert ps.deterministic_unambig_units == 0
        assert ps.em_routed_ambig_same_strand_units == 0
        assert ps.em_routed_ambig_opp_strand_units == 0

    def test_mutability(self):
        ps = PipelineStats()
        ps.total = 100
        ps.n_fragments = 42
        assert ps.total == 100
        assert ps.n_fragments == 42

    def test_to_dict_returns_all_fields(self):
        ps = PipelineStats()
        d = ps.to_dict()
        assert isinstance(d, dict)
        # Spot-check a few keys
        assert "total" in d
        assert "n_fragments" in d
        assert "deterministic_unambig_units" in d
        assert "em_routed_ambig_same_strand_units" in d
        assert "em_routed_ambig_opp_strand_units" in d
        # All values should be zero by default
        assert all(v == 0 for v in d.values())

    def test_to_dict_reflects_updates(self):
        ps = PipelineStats(total=50, n_fragments=10)
        d = ps.to_dict()
        assert d["total"] == 50
        assert d["n_fragments"] == 10
