"""Preflight reads a panel's expected conditions from its manifest, never from the caches it is checking, and checks an
index by loading it: a missing cache is reported missing, a stray one named, an index of another format refused at the
reference check rather than at every cache check downstream."""

from __future__ import annotations

import json

import preflight


def _panel(tmp_path, names, cached):
    root = tmp_path / "panel"
    (root / "scan_cache").mkdir(parents=True)
    (root / "manifest.json").write_text(json.dumps({"conditions": [{"name": n} for n in names]}))
    for n in cached:
        (root / "scan_cache" / n).mkdir()
    return root


def test_a_condition_without_a_cache_stays_expected(tmp_path):
    root = _panel(tmp_path, ["c1", "c2", "c3"], cached=["c1", "c3"])
    assert preflight._conditions(root) == ["c1", "c2", "c3"]
    assert preflight._cache_refusal(root / "scan_cache" / "c2", None) == "missing"


def test_a_cache_the_manifest_does_not_name_is_reported(tmp_path):
    root = _panel(tmp_path, ["c1"], cached=["c1", "old_name"])
    assert preflight._unexpected_caches(root, preflight._conditions(root)) == [
        ("old_name", "not a condition of the panel's manifest")
    ]


def test_an_unsimulated_panel_expects_nothing(tmp_path):
    assert preflight._conditions(tmp_path / "nothing_here") == []


def test_an_index_of_another_format_is_refused_where_it_is_checked(tmp_path):
    index = tmp_path / "index"
    index.mkdir()
    (index / "manifest.json").write_text(
        json.dumps({"format_version": 11, "rigel_version": "0.7.1"})
    )
    loaded, why = preflight._load_index(index)
    assert loaded is None and "format_version=11" in why
