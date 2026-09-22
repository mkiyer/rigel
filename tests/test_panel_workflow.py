"""The simulation and benchmarking workflow has one entry point, and its prerequisites are gated.

`scripts/sim/panel.py` sequences five expensive, resumable stages. The failure it exists to prevent
is running a stage for twenty minutes on a panel where the previous one never happened, so every
stage names its prerequisite and REFUSES, and these tests are what keep the refusals real.

No stage is executed here: each one costs minutes to hours and needs a panel tens of gigabytes
across. What is testable without that is the part that rots — path derivation, condition discovery,
the completeness rule for a cached condition, and every refusal — so this file says plainly what it
does not cover rather than reading as more coverage than it has
(`TRAPS: a-gate-that-reconstructs`).
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
    """One config in, every path out — the property that makes the workflow one command."""
    p = PANEL.Panel(_config(tmp_path))
    assert p.dir == tmp_path / "suite" / "mypanel"
    assert p.reference == tmp_path / "suite" / "reference"
    assert p.scan_cache == p.dir / "scan_cache"
    assert p.oracle_cache == p.dir / "oracle_cache"
    # with no `index:` the index is the reference directory's sibling, by convention
    assert p.index == tmp_path / "suite" / "rigel_index"


def test_the_index_is_the_configs_and_can_be_overridden(tmp_path):
    """The simulator reads `index:`, so the workflow must read the same one — a panel scored against an
    index other than the one it was simulated from is silently wrong. Without the key the reference's
    sibling is the convention (above), and `--index` overrides either."""
    cfg = _config(tmp_path, index=str(tmp_path / "idx"))
    assert PANEL.Panel(cfg).index == tmp_path / "idx"
    assert PANEL.Panel(cfg, index=tmp_path / "elsewhere").index == tmp_path / "elsewhere"


def test_the_probe_panel_is_the_one_the_capture_config_names(tmp_path):
    """A substrate whose probes are rendered (the test chromosome's BED) names them in its capture arm;
    the ladder's designed `capture_panel.tsv` is named the same way. `status` must check the file the
    config reads, or a fully built panel reads as unbuilt and `build` is named next."""
    probes = tmp_path / "rendered" / "probes.bed"
    capture = {"configs": [{"label": "on", "probes": str(probes)}]}
    p = PANEL.Panel(_config(tmp_path, capture=capture))
    assert p.probes == probes


@pytest.mark.parametrize("missing", ["genome", "gtf", "outdir"])
def test_a_config_missing_a_path_key_is_REFUSED(tmp_path, missing):
    """Not defaulted, not guessed. A panel that does not say where it lives cannot be driven, and
    inventing a path here is how a run writes tens of gigabytes into the wrong directory."""
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
    """The marker is `sim_oracle.bam`, not the directory. A condition whose simulation died leaves
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
    """A refusal that does not say what to run next is a traceback with better manners."""
    with pytest.raises(SystemExit) as e:
        PANEL.need(False, "the oracle cache", "panel.py cache")
    assert "the oracle cache" in str(e.value) and "panel.py cache" in str(e.value)
    PANEL.need(True, "satisfied", "never printed")  # must not raise


def test_score_REFUSES_without_the_oracle_cache(tmp_path):
    """The refusal this workflow was built for. Every truth-scoring instrument needs the
    origin-split oracle cache, so a recipe that builds it only as a side effect of some other
    instrument produces a panel every scorer rejects."""
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
    """A stage that failed must not look like a stage that was skipped — that is how a partial
    panel gets scored as a complete one."""
    with pytest.raises(SystemExit, match="FAILED"):
        PANEL.run([sys.executable, "-c", "raise SystemExit(3)"], what="a stage that fails")


# ── the completeness rule ────────────────────────────────────────────────────────────────────────


def test_an_oracle_condition_needs_every_part(tmp_path, capsys):
    """`status` counts a condition cached only when the three origin partitions, the per-strand RNA
    pair the certifier requires, and the undrained `_main` payload are all present. Counting
    directories would call a half-written condition done, and the next stage would fail deep inside an
    instrument instead of here."""
    p = PANEL.Panel(_config(tmp_path))
    (p.dir / "c1").mkdir(parents=True)
    (p.dir / "c1" / "sim_oracle.bam").touch()

    def write(part):
        (p.oracle_cache / "c1" / part).mkdir(parents=True)
        (p.oracle_cache / "c1" / part / "payload.npz").touch()

    for part in ("gdna", "mrna", "nrna", "_main"):  # no strand pair — deliberately incomplete
        write(part)
    PANEL.cmd_status(p, None)
    assert "oracle cache 0/1" in capsys.readouterr().out

    for part in ("rna_pos", "rna_neg"):
        write(part)
    PANEL.cmd_status(p, None)
    assert "oracle cache 1/1" in capsys.readouterr().out


def test_status_names_the_next_stage(tmp_path, capsys):
    """The whole point of `status`: not a dump, an instruction."""
    PANEL.cmd_status(PANEL.Panel(_config(tmp_path)), None)
    assert "next: `panel.py build`" in capsys.readouterr().out


def test_the_shipped_panel_configs_all_load(tmp_path):
    """A config the workflow cannot parse is a panel nobody can rebuild. `example_*.yaml` are
    documentation templates and are deliberately excluded."""
    cfgs = [
        c
        for c in (ROOT / "scripts" / "sim" / "configs").glob("*.yaml")
        if not c.name.startswith("example_")
    ]
    assert cfgs, "no panel configs found — this test would pass vacuously"
    for c in cfgs:
        p = PANEL.Panel(c)
        assert p.dir.name and p.index.name, f"{c.name} produced an empty path"


# ── the cache stage ──────────────────────────────────────────────────────────────────────────────


def test_cache_builds_the_scan_cache_then_builds_and_certifies_the_oracle_in_one_run(
    tmp_path, monkeypatch
):
    """Two commands, in order: the scan cache (a stage that must STOP the workflow on failure), then
    `calibration_oracle.py --build`, which builds every row's origin-split cache alike — the zero-gDNA
    rows too, there is no hold-out — and certifies `slot_truth.npz` in the same run. Its exit code is
    reported, never fatal: a failed FIELD gate still writes the COMPOSITION table."""
    p = PANEL.Panel(_config(tmp_path))
    for c in ("gdna_g00_ss_0.50_x", "gdna_g50_ss_0.50_x"):
        (p.dir / c).mkdir(parents=True)
        (p.dir / c / "sim_oracle.bam").touch()
    issued = []
    monkeypatch.setattr(
        PANEL, "run", lambda cmd, *, what: issued.append((what, [str(c) for c in cmd]))
    )
    launched = []

    class _Done:
        returncode = 1  # a certification gate "failed" — must be reported, not fatal

    monkeypatch.setattr(
        PANEL.subprocess, "run", lambda cmd, **k: (launched.append([str(c) for c in cmd]), _Done())[1]
    )
    args = type("A", (), {"jobs": 3, "conditions": ["gdna_g00_ss_0.50_x"], "force": False})()
    assert PANEL.cmd_cache(p, args) == 0
    assert [w for w, _ in issued] == [f"scan cache -> {p.scan_cache}"]
    assert len(launched) == 1
    cmd = launched[0]
    assert cmd[1].endswith("calibration_oracle.py")
    assert "--build" in cmd and cmd[cmd.index("--jobs") + 1] == "3"
    assert cmd[cmd.index("--condition") + 1] == "gdna_g00_ss_0.50_x"
    assert not any("prewarm" in c for c in cmd), "there is no zero-gDNA hold-out to pre-warm"

