"""The panel workflow's one entry point, `scripts/sim/panel.py`, and the refusals that gate its stages.

`panel.py` sequences five expensive, resumable stages, and each names its prerequisite and refuses without it, so no
stage runs for minutes on a panel whose previous stage never happened. No stage is executed here (each costs minutes
to hours on a panel tens of gigabytes across): these tests cover path derivation, condition discovery, the
completeness rule for a cached condition and every refusal, and nothing more.
"""

from __future__ import annotations

import json
import pathlib
import sys

import panel
import pytest
import yaml

ROOT = pathlib.Path(__file__).resolve().parents[1]


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


def _args(cfg, stage, *flags):
    """What `panel.py <stage> --config <cfg> <flags>` parses to, from the script's own parser."""
    return panel.build_parser().parse_args([stage, "--config", str(cfg), *flags])


# ── path derivation ──────────────────────────────────────────────────────────────────────────────


def test_every_path_derives_from_the_config(tmp_path):
    """One config in, every path out — the property that makes the workflow one command."""
    p = panel.Panel(_config(tmp_path))
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
    assert panel.Panel(cfg).index == tmp_path / "idx"
    assert panel.Panel(cfg, index=tmp_path / "elsewhere").index == tmp_path / "elsewhere"


def test_the_probe_panel_is_the_one_the_capture_config_names(tmp_path):
    """A substrate whose probes are rendered (the test chromosome's BED) names them in its capture arm;
    the ladder's designed `capture_panel.tsv` is named the same way. `status` must check the file the
    config reads, or a fully built panel reads as unbuilt and `build` is named next."""
    probes = tmp_path / "rendered" / "probes.bed"
    capture = {"configs": [{"label": "on", "probes": str(probes)}]}
    p = panel.Panel(_config(tmp_path, capture=capture))
    assert p.probes == probes


def test_every_path_is_expanded_as_the_simulator_expands_it(tmp_path, monkeypatch):
    """The in-repo configs name their paths under `$RIGEL_RUNS`. The workflow must derive the paths the
    simulator reads, or `status` checks one directory while `simulate` writes another. Perturbation:
    any one of the five without the expansion is a directory literally named `$RIGEL_TEST_ROOT`."""
    root = tmp_path / "root"
    monkeypatch.setenv("RIGEL_TEST_ROOT", str(root))
    paths = {
        "genome": "$RIGEL_TEST_ROOT/ref/genome.fa",
        "gtf": "$RIGEL_TEST_ROOT/ref/genes.gtf",
        "outdir": "${RIGEL_TEST_ROOT}/panel",
        "index": "$RIGEL_TEST_ROOT/idx",
        "capture": {"configs": [{"label": "on", "probes": "$RIGEL_TEST_ROOT/probes.bed"}]},
    }
    p = panel.Panel(_config(tmp_path, **paths))
    assert (p.genome, p.gtf, p.dir, p.index, p.probes) == (
        root / "ref" / "genome.fa",
        root / "ref" / "genes.gtf",
        root / "panel",
        root / "idx",
        root / "probes.bed",
    )


@pytest.mark.parametrize("missing", ["genome", "gtf", "outdir"])
def test_a_config_missing_a_path_key_is_refused(tmp_path, missing):
    """Not defaulted, not guessed. A panel that does not say where it lives cannot be driven, and
    inventing a path here is how a run writes tens of gigabytes into the wrong directory."""
    cfg = _config(tmp_path)
    body = yaml.safe_load(cfg.read_text())
    del body[missing]
    cfg.write_text(yaml.safe_dump(body))
    with pytest.raises(SystemExit, match=missing):
        panel.Panel(cfg)


def test_a_config_that_does_not_exist_is_refused(tmp_path):
    with pytest.raises(SystemExit, match="no such config"):
        panel.Panel(tmp_path / "nope.yaml")


# ── discovery ────────────────────────────────────────────────────────────────────────────────────


def test_a_condition_is_one_with_an_oracle_bam(tmp_path):
    """The marker is `sim_oracle.bam`, not the directory. A condition whose simulation died leaves
    the directory behind, and counting directories would report it as simulated."""
    p = panel.Panel(_config(tmp_path))
    (p.dir / "cond_a").mkdir(parents=True)
    (p.dir / "cond_a" / "sim_oracle.bam").touch()
    (p.dir / "cond_b").mkdir()  # directory only — a died-halfway condition
    assert p.conditions == ["cond_a"]


def test_no_panel_directory_is_zero_conditions_not_a_crash(tmp_path):
    assert panel.Panel(_config(tmp_path)).conditions == []


# ── the refusals ─────────────────────────────────────────────────────────────────────────────────


def test_need_names_what_to_run_next():
    """A refusal that does not say what to run next is a traceback with better manners."""
    with pytest.raises(SystemExit) as e:
        panel.need(False, "the oracle cache", "panel.py cache")
    assert "the oracle cache" in str(e.value) and "panel.py cache" in str(e.value)
    panel.need(True, "satisfied", "never printed")  # must not raise


def test_score_refuses_without_the_oracle_cache(tmp_path):
    """Every truth-scoring instrument needs the origin-split oracle cache, so `score` refuses without it: a
    recipe that builds it only as a side effect of some other instrument produces a panel every scorer rejects."""
    cfg = _config(tmp_path)
    with pytest.raises(SystemExit, match="oracle cache"):
        panel.cmd_score(panel.Panel(cfg), _args(cfg, "score", "--jobs", "1"))


def test_cache_refuses_without_simulated_conditions(tmp_path):
    cfg = _config(tmp_path)
    with pytest.raises(SystemExit, match="simulated conditions"):
        panel.cmd_cache(panel.Panel(cfg), _args(cfg, "cache", "--jobs", "1"))


def test_report_refuses_and_names_the_missing_arm(tmp_path):
    cfg = _config(tmp_path)
    with pytest.raises(SystemExit) as e:
        panel.cmd_report(panel.Panel(cfg), _args(cfg, "report", "--arms", "base", "oracle"))
    assert "qa_mypanel_base.jsonl" in str(e.value), "the refusal must name the file it looked for"


def test_a_failing_stage_stops_the_workflow():
    """A stage that failed must not look like a stage that was skipped — that is how a partial
    panel gets scored as a complete one."""
    with pytest.raises(SystemExit, match="FAILED"):
        panel.run([sys.executable, "-c", "raise SystemExit(3)"], what="a stage that fails")


# ── the completeness rule ────────────────────────────────────────────────────────────────────────


def test_an_oracle_condition_needs_every_part(tmp_path, capsys):
    """`status` counts a condition cached only when the three origin partitions and the per-strand RNA
    pair the certifier requires are all present. Counting directories would call a half-written
    condition done, and the next stage would fail deep inside an instrument instead of here."""
    p = panel.Panel(_config(tmp_path))
    (p.dir / "c1").mkdir(parents=True)
    (p.dir / "c1" / "sim_oracle.bam").touch()

    def write(part):
        (p.oracle_cache / "c1" / part).mkdir(parents=True)
        (p.oracle_cache / "c1" / part / "payload.npz").touch()

    for part in ("gdna", "annotated", "nascent"):  # no strand pair — deliberately incomplete
        write(part)
    panel.cmd_status(p, None)
    assert "oracle cache 0/1" in capsys.readouterr().out

    for part in ("rna_pos", "rna_neg"):
        write(part)
    panel.cmd_status(p, None)
    assert "oracle cache 1/1" in capsys.readouterr().out


def test_status_names_the_next_stage(tmp_path, capsys):
    """The whole point of `status`: not a dump, an instruction."""
    panel.cmd_status(panel.Panel(_config(tmp_path)), None)
    assert "next: `panel.py build`" in capsys.readouterr().out


def test_the_shipped_panel_configs_all_load(tmp_path):
    """A config the workflow cannot parse is a panel nobody can rebuild."""
    cfgs = sorted((ROOT / "scripts" / "sim" / "configs").glob("*.yaml"))
    assert cfgs, "no panel configs found — this test would pass vacuously"
    for c in cfgs:
        p = panel.Panel(c)
        assert p.dir.name and p.index.name, f"{c.name} produced an empty path"


# ── the build stage ──────────────────────────────────────────────────────────────────────────────


def _built_reference(tmp_path, **over):
    cfg = _config(tmp_path, **over)
    p = panel.Panel(cfg)
    p.genome.touch()
    p.gtf.touch()
    return cfg, p


def test_build_force_never_redesigns_a_rendered_bed_panel(tmp_path, monkeypatch):
    """The test chromosome's probes are a BED its renderer draws to a design (`build_test_reference.py`); the probe
    designer writes a TSV with its BED beside it, so handed a `.bed` target it would replace the rendered panel with a
    random one, through two handles on one file. `build --force` rebuilds the index and leaves the BED alone."""
    probes = tmp_path / "rendered" / "test_probes.bed"
    probes.parent.mkdir()
    probes.write_text("rendered\n")
    cfg, p = _built_reference(
        tmp_path, capture={"configs": [{"label": "on", "probes": str(probes)}]}
    )
    issued = []
    monkeypatch.setattr(panel, "run", lambda cmd, *, what: issued.append([str(c) for c in cmd]))
    assert panel.cmd_build(p, _args(cfg, "build", "--force")) == 0
    assert not any("design_suite_probes" in c for cmd in issued for c in cmd), issued
    assert probes.read_text() == "rendered\n"


def test_build_refuses_a_missing_rendered_bed_naming_its_renderer(tmp_path, monkeypatch):
    probes = tmp_path / "rendered" / "test_probes.bed"
    cfg, p = _built_reference(
        tmp_path, capture={"configs": [{"label": "on", "probes": str(probes)}]}
    )
    monkeypatch.setattr(panel, "run", lambda cmd, *, what: None)
    with pytest.raises(SystemExit, match="build_test_reference"):
        panel.cmd_build(p, _args(cfg, "build", "--force"))


def test_build_force_redesigns_a_designed_tsv_panel(tmp_path, monkeypatch):
    """The ladder's `capture_panel.tsv` is the designer's own output, so `--force` redesigns it."""
    cfg, p = _built_reference(tmp_path)
    issued = []
    monkeypatch.setattr(panel, "run", lambda cmd, *, what: issued.append([str(c) for c in cmd]))
    assert panel.cmd_build(p, _args(cfg, "build", "--force")) == 0
    (design,) = [cmd for cmd in issued if any("design_suite_probes" in c for c in cmd)]
    assert design[design.index("-o") + 1] == str(p.probes) and p.probes.suffix == ".tsv"


# ── the cache stage ──────────────────────────────────────────────────────────────────────────────


def test_cache_builds_the_scan_cache_then_builds_and_certifies_the_oracle_in_one_run(
    tmp_path, monkeypatch
):
    """Two commands, in order: the scan cache (a stage that must stop the workflow on failure), then
    `calibration_oracle.py --build`, which builds every row's origin-split cache alike, the zero-gDNA
    rows too, and certifies `slot_truth.npz` in the same run. Its exit code is reported, never fatal:
    a failed field gate still writes the composition table."""
    cfg = _config(tmp_path)
    p = panel.Panel(cfg)
    for c in ("gdna_g00_ss_0.50_x", "gdna_g50_ss_0.50_x"):
        (p.dir / c).mkdir(parents=True)
        (p.dir / c / "sim_oracle.bam").touch()
    issued = []
    monkeypatch.setattr(
        panel, "run", lambda cmd, *, what: issued.append((what, [str(c) for c in cmd]))
    )
    launched = []

    class _Done:
        returncode = 1  # a certification gate "failed" — must be reported, not fatal

    monkeypatch.setattr(
        panel.subprocess,
        "run",
        lambda cmd, **k: (launched.append([str(c) for c in cmd]), _Done())[1],
    )
    args = _args(cfg, "cache", "--jobs", "3", "--conditions", "gdna_g00_ss_0.50_x")
    assert panel.cmd_cache(p, args) == 0
    assert [w for w, _ in issued] == [f"scan cache -> {p.scan_cache}"]
    assert len(launched) == 1
    cmd = launched[0]
    assert cmd[1].endswith("calibration_oracle.py")
    assert "--build" in cmd and cmd[cmd.index("--jobs") + 1] == "3"
    assert cmd[cmd.index("--condition") + 1] == "gdna_g00_ss_0.50_x"
    assert not any("prewarm" in c for c in cmd), "there is no zero-gDNA hold-out to pre-warm"
    assert "--force" not in cmd and not any("--force" in c for _, c in issued)

    # --force reaches both caches: the scan cache and the origin-split partitions, keyed alike and alike stale
    issued.clear()
    launched.clear()
    args.force = True
    assert panel.cmd_cache(p, args) == 0
    assert "--force" in issued[0][1], "--force did not reach the scan cache"
    assert "--force" in launched[0], "--force did not reach the oracle's partitions"


# ── an aligned input: the ladder's STAR rung ─────────────────────────────────────────────────────


def test_the_parameter_file_is_read_as_the_aligner_reads_it(tmp_path):
    """`align` passes each line of the parameter file as `--<option> <values>`, as STAR's own
    `--parametersFiles` reads it (rustar takes none): a comment or a blank line is no option, and an
    option keeps every value on its line."""
    f = tmp_path / "p.txt"
    f.write_text(
        "# a comment\n\noutSAMtype BAM Unsorted\n  # indented comment\nalignSJstitchMismatchNmax 5 -1 5 5\n"
    )
    assert panel.star_options(f) == [
        "--outSAMtype",
        "BAM",
        "Unsorted",
        "--alignSJstitchMismatchNmax",
        "5",
        "-1",
        "5",
        "5",
    ]


def test_the_rung_aligns_with_the_manuals_recommended_star_settings():
    """The STAR rung is the input the user's manual tells a user to make: its parameter file carries the two
    settings the manual recommends, so the benchmark cannot drift from the advice."""
    opts = " ".join(panel.star_options(panel.STAR_PARAMS))
    manual = (ROOT / "docs" / "MANUAL.md").read_text()
    for setting in ("--sjdbScore 0", "--outFilterMultimapScoreRange 2"):
        assert setting in opts, f"the rung's parameters lack {setting}"
        assert setting in manual, f"the MANUAL no longer recommends {setting}"


def test_align_refuses_without_the_aligner_or_its_index(tmp_path):
    cfg = _config(tmp_path)
    p = panel.Panel(cfg)
    (p.dir / "c1").mkdir(parents=True)
    (p.dir / "c1" / "sim_oracle.bam").write_bytes(b"")
    args = _args(cfg, "align", "--threads", "1")
    with pytest.raises(SystemExit, match="a STAR or rustar binary"):
        panel.cmd_align(p, args)
    args.aligner = sys.executable  # a file that exists
    args.aligner_index = tmp_path  # with no STAR index in it
    with pytest.raises(SystemExit, match="genome index"):
        panel.cmd_align(p, args)


def _fake_aligner(tmp_path):
    """An executable that answers `--version`, writes the pair it is handed (read through `--readFilesCommand`) to
    `<prefix>reads_seen`, and writes the condition's simulator BAM as its alignment, every template, as `align`
    checks."""
    script = tmp_path / "fake_aligner"
    script.write_text(
        f"#!{sys.executable}\n"
        "import pathlib, shutil, subprocess, sys\n"
        "a = sys.argv[1:]\n"
        "if a == ['--version']:\n"
        "    print('fake 1.0')\n"
        "    sys.exit(0)\n"
        "i = a.index('--readFilesIn')\n"
        "cmd = a[a.index('--readFilesCommand') + 1]\n"
        "prefix = a[a.index('--outFileNamePrefix') + 1]\n"
        "seen = b''.join(subprocess.run([cmd, r], capture_output=True, check=True).stdout for r in a[i + 1 : i + 3])\n"
        "pathlib.Path(prefix + 'reads_seen').write_bytes(seen)\n"
        "shutil.copy(pathlib.Path(prefix).parent / 'sim_oracle.bam', prefix + 'Aligned.out.bam')\n"
    )
    script.chmod(0o755)
    return script


def _simulator_bam(path, pairs):
    """A simulator BAM as the engine writes one: each pair's R1 forward, then its R2 on the reverse strand, whose
    stored sequence is the read's reverse complement."""
    import pysam

    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 10_000}]}
    )
    comp = str.maketrans("ACGT", "TGCA")
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for name, r1, r2 in pairs:
            for stored, flag, pos, mate in (
                (r1, 99, 100, 300),
                (r2[::-1].translate(comp), 147, 300, 100),
            ):
                rec = pysam.AlignedSegment(header)
                rec.query_name, rec.query_sequence, rec.flag = name, stored, flag
                rec.reference_id, rec.reference_start, rec.cigar = 0, pos, [(0, len(stored))]
                rec.next_reference_id, rec.next_reference_start = 0, mate
                rec.query_qualities = pysam.qualitystring_to_array("I" * len(stored))
                out.write(rec)


def test_align_rebuilds_a_missing_pair_from_the_simulators_bam(tmp_path):
    """The panels keep no FASTQs: the simulator's BAM holds every read as it was written, so `align` rebuilds a
    missing pair with `samtools fastq -N`, byte for byte the FASTQs the simulator wrote, aligns it and deletes it.
    With neither the reads nor the BAM it refuses, naming both.

    Perturbation: no rebuild refuses the condition; a rebuild without `-N` drops the `/1`, `/2` the simulator's
    names carry; a rebuilt pair left behind is what the panels shed."""
    cfg = _config(tmp_path)
    p = panel.Panel(cfg)
    d = p.dir / "c1"
    d.mkdir(parents=True)
    pairs = [
        ("rna_T1:10-210:r:0", "ACGTTGCA", "GGATCCAT"),
        ("gdna_chr1:400-640:f:1", "TTGACCAG", "CAGTACGT"),
    ]
    _simulator_bam(d / "sim_oracle.bam", pairs)
    (d / "truth_summary.json").write_text(
        json.dumps({"origin_counts": {"annotated": 1, "gdna": 1}})
    )
    index = tmp_path / "index"
    index.mkdir()
    (index / "SA").write_bytes(b"")
    args = _args(
        cfg,
        "align",
        "--conditions",
        "c1",
        "--aligner",
        str(_fake_aligner(tmp_path)),
        "--aligner-index",
        str(index),
        "--threads",
        "1",
    )
    assert panel.cmd_align(p, args) == 0
    expected = "".join(f"@{n}/1\n{r1}\n+\n{'I' * len(r1)}\n" for n, r1, _ in pairs) + "".join(
        f"@{n}/2\n{r2}\n+\n{'I' * len(r2)}\n" for n, _, r2 in pairs
    )
    assert (d / "rustar_recommended.reads_seen").read_text() == expected
    assert not (d / "sim_R1.fq.gz").exists() and not (d / "sim_R2.fq.gz").exists(), (
        "the rebuilt pair was left behind"
    )
    assert (d / "rustar_recommended.align.json").is_file()
    (d / "sim_oracle.bam").unlink()
    args.force = True
    with pytest.raises(SystemExit, match="simulator's BAM"):
        panel.cmd_align(p, args)


def test_an_aligned_input_counts_only_with_its_record(tmp_path):
    """`status` counts a condition aligned only when `align` finished it: the BAM and the record it writes
    last, after the BAM's templates were checked against the simulator's. A BAM alone is an interrupted run."""
    p = panel.Panel(_config(tmp_path))
    for c in ("c1", "c2", "c3"):
        (p.dir / c).mkdir(parents=True)
        (p.dir / c / "sim_oracle.bam").write_bytes(b"")
    for c in ("c1", "c2"):
        (p.dir / c / "rustar_recommended.bam").write_bytes(b"")
    for c in ("c1", "c3"):  # c3's record outlived its BAM
        (p.dir / c / "rustar_recommended.align.json").write_text("{}")
    assert panel.aligned_inputs(p) == {"rustar_recommended": 1}


def test_an_aligned_input_scores_into_an_arm_file_of_its_own(tmp_path):
    """The simulator's alignments and an aligner's are two inputs to one arm, and each keeps its own
    file, so neither overwrites the other and `report` names the one it lacks."""
    cfg = _config(tmp_path)
    p = panel.Panel(cfg)
    assert panel.arm_file(p, "base", None).name == "qa_mypanel_base.jsonl"
    assert (
        panel.arm_file(p, "base", "rustar_recommended").name
        == "qa_mypanel_base__rustar_recommended.jsonl"
    )
    with pytest.raises(SystemExit) as e:
        panel.cmd_report(p, _args(cfg, "report", "--bam", "rustar_recommended"))
    assert "qa_mypanel_base__rustar_recommended.jsonl" in str(e.value)
