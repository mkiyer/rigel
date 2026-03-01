"""Tests for hulkrna.stats — PipelineStats and _BamStatsProxy."""

import pytest

from hulkrna.stats import PipelineStats, _BamStatsProxy


class TestPipelineStats:
    def test_defaults_are_zero(self):
        ps = PipelineStats()
        assert ps.total == 0
        assert ps.n_fragments == 0
        assert ps.deterministic_unique_units == 0
        assert ps.em_routed_isoform_ambig_units == 0
        assert ps.em_routed_gene_ambig_units == 0

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
        assert "deterministic_unique_units" in d
        assert "em_routed_isoform_ambig_units" in d
        assert "em_routed_gene_ambig_units" in d
        # All values should be zero by default
        assert all(v == 0 for v in d.values())

    def test_to_dict_reflects_updates(self):
        ps = PipelineStats(total=50, n_fragments=10)
        d = ps.to_dict()
        assert d["total"] == 50
        assert d["n_fragments"] == 10


class TestBamStatsProxy:
    def test_setitem_propagates(self):
        ps = PipelineStats()
        proxy = ps.as_bam_stats_dict()
        proxy["total"] = 999
        assert ps.total == 999

    def test_getitem_reads_from_dataclass(self):
        ps = PipelineStats(total=42)
        proxy = ps.as_bam_stats_dict()
        assert proxy["total"] == 42

    def test_setdefault(self):
        ps = PipelineStats()
        proxy = ps.as_bam_stats_dict()
        proxy.setdefault("total", 5)  # already exists (0)
        assert proxy["total"] == 0  # should not overwrite
        proxy.setdefault("new_key", 99)  # not a BAM key
        assert proxy["new_key"] == 99

    def test_non_bam_keys_allowed(self):
        ps = PipelineStats()
        proxy = ps.as_bam_stats_dict()
        proxy["custom_stat"] = 123
        assert proxy["custom_stat"] == 123
        # But it should NOT set a field on PipelineStats
        assert not hasattr(ps, "custom_stat") or getattr(ps, "custom_stat", None) != 123

    def test_isinstance_dict(self):
        ps = PipelineStats()
        proxy = ps.as_bam_stats_dict()
        assert isinstance(proxy, dict)

    def test_bam_keys_are_present(self):
        ps = PipelineStats()
        proxy = ps.as_bam_stats_dict()
        for key in ("total", "qc_fail", "unmapped", "secondary",
                     "supplementary", "duplicate", "n_read_names",
                     "unique", "multimapping", "proper_pair",
                     "improper_pair", "mate_unmapped"):
            assert key in proxy
