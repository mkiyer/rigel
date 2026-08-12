"""⭐⭐ THE SIMULATION + BENCHMARKING WORKFLOW HAS ONE ENTRY POINT, AND ITS PREREQUISITES ARE GATED.

`scripts/sim/panel.py` sequences five expensive, resumable stages. The failure it exists to prevent is
running stage 4 for twenty minutes on a panel where stage 3 never happened — so every stage names its
prerequisite and REFUSES, and these tests are what keep the refusals real.

⛔ **No stage is executed here.** Each one costs minutes to hours and needs a 27 GB panel; what is
testable without that is the part that rots: path derivation, condition discovery, the completeness rule
for a cached condition, and every refusal. `TRAPS: a-gate-that-reconstructs` — a test that reads as more
coverage than it has is worse than none, so this file says plainly what it does not cover.
"""

from __future__ import annotations

import importlib.util
import pathlib
import sys

import pytest
import yaml

ROOT = pathlib.Path(__file__).resolve().parents[1]
PANEL_PY = ROOT / "scripts" / "sim" / "panel.py"


def _load():
    """Import `panel.py` the way running it imports it (its own dir on the path)."""
    spec = importlib.util.spec_from_file_location("_panel_mod", PANEL_PY)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["_panel_mod"] = mod
    sys.path.insert(0, str(PANEL_PY.parent))
    try:
        spec.loader.exec_module(mod)
    finally:
        sys.path.remove(str(PANEL_PY.parent))
        sys.modules.pop("_panel_mod", None)
    return mod


PANEL = _load()


def _config(tmp_path, **over):
    """A minimal panel config — the three keys the workflow derives every path from."""
    ref = tmp_path / "suite" / "reference"
    ref.mkdir(parents=True)
    body = {
        "genome": str(ref / "genome.fa"),
        "gtf": str(ref / "genes.gtf"),
        "outdir": str(tmp_path / "suite" / "mypanel"),
        **over,
    }
    p = tmp_path / "cfg.yaml"
    p.write_text(yaml.safe_dump(body))
    return p


# ── path derivation ──────────────────────────────────────────────────────────────────────────────


def test_every_path_derives_from_the_config(tmp_path):
    """⭐ One config in, every path out — the property that makes the workflow one command."""
    p = PANEL.Panel(_config(tmp_path))
    assert p.dir == tmp_path / "suite" / "mypanel"
    assert p.reference == tmp_path / "suite" / "reference"
    assert p.scan_cache == p.dir / "scan_cache"
    assert p.oracle_cache == p.dir / "oracle_cache"
    # ⚠ the index is the ONE derived-by-convention path: a sibling of the reference directory.
    assert p.index == tmp_path / "suite" / "rigel_index"


def test_the_index_can_be_overridden(tmp_path):
    """⛔ Because the index is a CONVENTION, not a config key, it must be overridable — otherwise a
    panel built against a non-default index is silently scored against the wrong one."""
    p = PANEL.Panel(_config(tmp_path), index=tmp_path / "elsewhere")
    assert p.index == tmp_path / "elsewhere"


@pytest.mark.parametrize("missing", ["genome", "gtf", "outdir"])
def test_a_config_missing_a_path_key_is_REFUSED(tmp_path, missing):
    """⛔ Not defaulted, not guessed. A panel that does not say where it lives cannot be driven, and
    inventing a path here is how a run writes 27 GB into the wrong directory."""
    cfg = _config(tmp_path)
    body = yaml.safe_load(cfg.read_text())
    del body[missing]
    cfg.write_text(yaml.safe_dump(body))
    with pytest.raises(SystemExit, match=missing):
        PANEL.Panel(cfg)


def test_a_config_that_does_not_exist_is_REFUSED(tmp_path):
    with pytest.raises(SystemExit, match="no such config"):
        PANEL.Panel(tmp_path / "nope.yaml")


# ── discovery ────────────────────────────────────────────────────────────────────────────────────


def test_a_condition_is_one_with_an_ORACLE_BAM(tmp_path):
    """⛔ The marker is `sim_oracle.bam`, not the directory. A condition whose simulation died leaves
    the directory behind, and counting directories would report it as simulated."""
    p = PANEL.Panel(_config(tmp_path))
    (p.dir / "cond_a").mkdir(parents=True)
    (p.dir / "cond_a" / "sim_oracle.bam").touch()
    (p.dir / "cond_b").mkdir()  # directory only — a died-halfway condition
    assert p.conditions == ["cond_a"]


def test_no_panel_directory_is_zero_conditions_not_a_crash(tmp_path):
    assert PANEL.Panel(_config(tmp_path)).conditions == []


# ── the refusals ─────────────────────────────────────────────────────────────────────────────────


def test_need_names_the_fix(tmp_path):
    """⭐ A refusal that does not say what to run next is a traceback with better manners."""
    with pytest.raises(SystemExit) as e:
        PANEL.need(False, "the oracle cache", "panel.py cache")
    assert "the oracle cache" in str(e.value) and "panel.py cache" in str(e.value)
    PANEL.need(True, "satisfied", "never printed")  # must not raise


def test_score_REFUSES_without_the_oracle_cache(tmp_path):
    """⛔⛔ THE REFUSAL THIS WORKFLOW WAS BUILT FOR. Every truth-scoring instrument needs the
    origin-split oracle cache, and until 2026-08-11 building it was a SIDE EFFECT of an unrelated
    instrument — so the documented recipe produced a panel that every scorer rejected."""
    p = PANEL.Panel(_config(tmp_path))
    args = type("A", (), {"jobs": 1, "arms": ["base"], "conditions": None})()
    with pytest.raises(SystemExit, match="oracle cache"):
        PANEL.cmd_score(p, args)


def test_cache_REFUSES_without_simulated_conditions(tmp_path):
    p = PANEL.Panel(_config(tmp_path))
    args = type("A", (), {"jobs": 1, "conditions": None})()
    with pytest.raises(SystemExit, match="simulated conditions"):
        PANEL.cmd_cache(p, args)


def test_report_REFUSES_and_names_the_missing_arm(tmp_path):
    p = PANEL.Panel(_config(tmp_path))
    args = type("A", (), {"arms": ["base", "oracle"]})()
    with pytest.raises(SystemExit) as e:
        PANEL.cmd_report(p, args)
    assert "qa_mypanel_base.jsonl" in str(e.value), "the refusal must name the file it looked for"


def test_a_failing_stage_STOPS_the_workflow():
    """⛔ A stage that failed must not look like a stage that was skipped — that is how a partial
    panel gets scored as a complete one."""
    with pytest.raises(SystemExit, match="FAILED"):
        PANEL.run([sys.executable, "-c", "raise SystemExit(3)"], what="a stage that fails")


# ── the completeness rule ────────────────────────────────────────────────────────────────────────


def test_an_oracle_condition_needs_ALL_FOUR_PARTS(tmp_path, capsys):
    """⛔ `status` counts a condition cached only when `gdna`, `mrna`, `nrna` AND the undrained
    `_main` payload are all present. Counting directories would call a half-written condition done,
    and the next stage would fail deep inside an instrument instead of here."""
    p = PANEL.Panel(_config(tmp_path))
    (p.dir / "c1").mkdir(parents=True)
    (p.dir / "c1" / "sim_oracle.bam").touch()
    for part in ("gdna", "mrna", "nrna"):  # three of four — deliberately incomplete
        (p.oracle_cache / "c1" / part).mkdir(parents=True)
        (p.oracle_cache / "c1" / part / "payload.npz").touch()
    PANEL.cmd_status(p, None)
    assert "oracle cache 0/1" in capsys.readouterr().out

    (p.oracle_cache / "c1" / "_main").mkdir(parents=True)
    (p.oracle_cache / "c1" / "_main" / "payload.npz").touch()
    PANEL.cmd_status(p, None)
    assert "oracle cache 1/1" in capsys.readouterr().out


def test_status_names_the_next_stage(tmp_path, capsys):
    """⭐ The whole point of `status`: not a dump, an instruction."""
    PANEL.cmd_status(PANEL.Panel(_config(tmp_path)), None)
    assert "next: `panel.py build`" in capsys.readouterr().out


def test_the_shipped_panel_configs_all_load(tmp_path):
    """⛔ A config the workflow cannot parse is a panel nobody can rebuild. ⚠ `example_*.yaml` are
    documentation templates and are deliberately excluded."""
    cfgs = [c for c in (ROOT / "scripts" / "sim" / "configs").glob("*.yaml")
            if not c.name.startswith("example_")]
    assert cfgs, "no panel configs found — this test would pass vacuously"
    for c in cfgs:
        p = PANEL.Panel(c)
        assert p.dir.name and p.index.name, f"{c.name} produced an empty path"
