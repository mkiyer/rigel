"""The whole-genome simulator's configuration, a run's own refusals, and its nascent-abundance helpers.

The YAML parser and its refusals; a run's refusals of a GTF its index was not built from and of a run
that keeps no reads, and the manifest's record of that index; the fragment-length pmf the engine draws
from; and the sparse nascent model: nascent lives on the index's nRNA entities, the on/off draw is per
gene SPAN rather than per isoform, the level is log-uniform over the configured range, and it is drawn
independently of the mature level. Most of these properties are only observable statistically, so
each gate states the band it uses and where that band comes from, and several exist because a
perturbation to the mechanism passed every other gate in the file.
"""

import hashlib
import importlib.util
import json
import os
import re
import sys
from pathlib import Path

import numpy as np
import pytest

from rigel.sim.whole_genome import apply_sparse_nrna, parse_yaml_config, run_simulation
from rigel.transcript import Transcript
from rigel.types import Interval, Strand

ROOT = Path(__file__).resolve().parents[1]


def _transcript(
    transcript_id: str,
    abundance: float,
    exons: list[tuple[int, int]],
) -> Transcript:
    return Transcript(
        ref="chr1",
        strand=Strand.POS,
        exons=[Interval(start, end) for start, end in exons],
        t_id=transcript_id,
        g_id=f"gene_{transcript_id}",
        abundance=abundance,
    )


def _entity(transcript_id: str, start: int, end: int, t_index: int) -> Transcript:
    """An nRNA entity as the index makes it: single-exon over the clustered span, synthetic."""
    t = Transcript(
        ref="chr1",
        strand=Strand.POS,
        exons=[Interval(start, end)],
        t_id=transcript_id,
        g_id=transcript_id,
        t_index=t_index,
        is_single_exon=True,
        is_synthetic=True,
        abundance=0.0,
    )
    t.length = t.compute_length()
    return t


def test_sparse_nascent_pools_onto_entities_and_the_SPAN_is_the_unit():
    """Nascent molecules live on the index's nRNA ENTITIES, and under the sparse model the entity is
    also the unit of the on/off draw — one gene span is transcribed or it is not, however many
    isoforms share it. A span whose only contributor is SILENT gets nothing, single-exon transcripts
    contribute nothing, and the contributors themselves end at ``nrna_abundance = 0``. Drawing per
    CONTRIBUTOR instead compounds the probability over the isoforms, so the intron slots calibration
    reads would be several times less sparse than configured."""
    multi_a = _transcript("multi_a", 100.0, [(0, 100), (200, 300)])
    multi_b = _transcript("multi_b", 200.0, [(400, 500), (700, 800)])
    single = _transcript("single", 300.0, [(900, 1200)])
    off = _transcript("off", 0.0, [(1300, 1400), (1500, 1600)])
    shared = _entity("N_shared", 0, 800, 4)  # multi_a and multi_b cluster onto one span
    lone = _entity("N_off", 1300, 1600, 5)
    for i, t in enumerate((multi_a, multi_b, single, off)):
        t.t_index = i
    multi_a.span_t_index = multi_b.span_t_index = 4
    off.span_t_index = 5
    transcripts = [multi_a, multi_b, single, off, shared, lone]

    realized_ratio = apply_sparse_nrna(transcripts, (150.0, 150.0), on_fraction=1.0, seed=17)

    # the two contributors SHARE one span, so the span is ONE draw — not one per isoform
    assert shared.nrna_abundance == pytest.approx(150.0)
    assert lone.nrna_abundance == 0.0, "a span whose only contributor is silent is not transcribed"
    assert all(t.nrna_abundance == 0.0 for t in (multi_a, multi_b, single, off))
    assert realized_ratio == pytest.approx(150.0 / 600.0)


def test_a_multi_exon_contributor_with_no_entity_is_a_defect_not_a_skip():
    """Every multi-exon transcript of a rigel index links to an entity; a transcript list that does
    not is not the index's, and pooling onto a missing entity raises rather than dropping molecules."""
    from rigel.sim.whole_genome import apply_nrna_ratio

    orphan = _transcript("orphan", 100.0, [(0, 100), (200, 300)])
    orphan.t_index = 0
    with pytest.raises(ValueError, match="nascent entity"):
        apply_nrna_ratio([orphan], 0.25)


def test_the_fl_pmf_is_the_one_the_engine_draws_from():
    """The share solve integrates over the fragment-length pmf. If that pmf were not the engine's,
    the solved ratio would be right about a distribution nothing samples.

    PERTURBATION: the analytic pmf must match a large SAMPLE from `truncated_normal_frag_lengths`,
    and must NOT match one drawn with a different sd.
    """
    import numpy as np

    from rigel.sim.sampling import truncated_normal_frag_lengths
    from rigel.sim.wgs_config import SimulationParams
    from rigel.sim.whole_genome import fl_pmf

    sim = SimulationParams(frag_mean=206, frag_std=98, frag_min=50, frag_max=500)
    widths, p = fl_pmf(sim)
    rng = np.random.default_rng(0)
    draws = truncated_normal_frag_lengths(
        rng, 400_000, sim.frag_mean, sim.frag_std, sim.frag_min, sim.frag_max
    )
    emp = np.bincount(draws - sim.frag_min, minlength=len(widths))[: len(widths)] / len(draws)
    assert np.abs(emp - p).max() < 2e-3, "the analytic pmf must be the sampler's"
    wrong = fl_pmf(SimulationParams(frag_mean=206, frag_std=40, frag_min=50, frag_max=500))[1]
    assert np.abs(emp - wrong).max() > 5e-3, "a different sd must be distinguishable"


def test_an_empty_fragment_length_window_is_refused_not_drawn_from_forever():
    """No length lies in [400, 300], so the rejection loop never fills: `rigel sim` hung on a scenario whose
    frag_min exceeded its frag_max. The generator here raises after a bounded number of draws, so without the
    refusal this gate fails instead of hanging. A one-length window is not empty."""
    from rigel.sim.sampling import truncated_normal_frag_lengths

    class Bounded:
        def __init__(self):
            self.rng, self.draws = np.random.default_rng(0), 0

        def normal(self, *args):
            self.draws += 1
            if self.draws > 100:
                raise RuntimeError("the sampler kept drawing from an empty window")
            return self.rng.normal(*args)

    with pytest.raises(ValueError, match=r"\[400, 300\] holds no length"):
        truncated_normal_frag_lengths(Bounded(), 5, 250.0, 50.0, 400, 300)
    assert truncated_normal_frag_lengths(Bounded(), 5, 250.0, 50.0, 250, 250).tolist() == [250] * 5


def test_parse_sparse_nrna_config(tmp_path):
    config_path = tmp_path / "sim.yaml"
    config_path.write_text(
        "genome: /tmp/genome.fa\n"
        "gtf: /tmp/genes.gtf\n"
        "nrna:\n"
        "  mode: sparse\n"
        "  abundance_ranges:\n"
        "    - [1.0, 10.0]\n"
        "    - [10.0, 1000.0]\n"
        "  ratio_labels: [low, mixed]\n"
        "  on_fraction: 0.5\n"
        "  seed: 123\n"
    )

    config = parse_yaml_config(config_path)

    assert config.nrna.mode == "sparse"
    assert config.nrna.abundance_ranges == [(1.0, 10.0), (10.0, 1000.0)]
    assert config.nrna.ratio_labels == ["low", "mixed"]
    assert config.nrna.on_fraction == 0.5
    assert config.nrna.seed == 123


def test_sparse_requires_abundance_ranges(tmp_path):
    config_path = tmp_path / "sim.yaml"
    config_path.write_text("genome: /tmp/genome.fa\ngtf: /tmp/genes.gtf\nnrna:\n  mode: sparse\n")

    with pytest.raises(ValueError, match="abundance_ranges"):
        parse_yaml_config(config_path)


def test_a_log_uniform_range_may_not_touch_zero(tmp_path):
    """A log-uniform draw has no zero end, so ``[0, x]`` is not "no nascent" — it is undefined.
    Absence is expressed by ``on_fraction``, which is the model's own switch."""
    config_path = tmp_path / "sim.yaml"
    config_path.write_text(
        "genome: /tmp/genome.fa\ngtf: /tmp/genes.gtf\nnrna:\n  mode: sparse\n"
        "  abundance_ranges:\n    - [0.0, 10.0]\n  ratio_labels: [bad]\n"
    )
    with pytest.raises(ValueError, match="0 < min"):
        parse_yaml_config(config_path)


def test_an_unknown_key_is_refused_at_the_top_level_and_in_each_section(tmp_path):
    """A misspelt key must not run on the default of the key it meant."""
    config_path = tmp_path / "sim.yaml"
    for body, key in [
        ("strand_specificity: [0.5]\n", "strand_specificity"),
        ("simulation:\n  n_total_fragment: 1000\n", "n_total_fragment"),
        ("abundance:\n  frac_expresed: 0.5\n", "frac_expresed"),
        ("nrna:\n  ratio: [0.1]\n", "ratio"),
        ("gdna:\n  rate: [0.5]\n", "rate"),
    ]:
        config_path.write_text("genome: /tmp/genome.fa\ngtf: /tmp/genes.gtf\n" + body)
        with pytest.raises(ValueError, match=rf"unknown key\(s\) \['{key}'\]"):
            parse_yaml_config(config_path)


@pytest.mark.parametrize(
    ("body", "key"),
    [
        ("nrna:\n  shares: [0.2]\n  ratio_labels: [s20]\n", "shares"),
        ("transcript_filter: all\n", "transcript_filter"),
    ],
    ids=["shares", "transcript_filter"],
)
def test_a_retired_key_is_refused_not_run_on_a_default(tmp_path, body, key):
    """`shares` belonged to the deleted fragment-share nascent mode. Still allowed and never read,
    `nrna: {shares: [0.2], ratio_labels: [s20]}` ran `additive_ratio` at ratio 0.0: a nascent-free
    condition named `s20`. `transcript_filter` has nothing to filter once the transcriptome is the
    index's, and `all` was accepted silently. Both are unknown keys. PERTURBATION: putting either back
    on its allow-list fails this gate."""
    config_path = tmp_path / "sim.yaml"
    config_path.write_text("genome: /tmp/genome.fa\ngtf: /tmp/genes.gtf\n" + body)
    with pytest.raises(ValueError, match=rf"unknown key\(s\) \['{key}'\]"):
        parse_yaml_config(config_path)


def test_no_oracle_bam_and_no_fastq_is_refused_when_the_config_is_read(tmp_path):
    """With no oracle BAM written and the FASTQs deleted once the truth is, a condition keeps no read
    and loses its resume key, so the run exits 0 with nothing to score and the next run simulates it
    all again. Either output alone is a panel. PERTURBATION: `or` for `and` refuses the fl-gap panels'
    `emit_fastq: false`; deleting the parse-time check passes the last config."""
    config_path = tmp_path / "sim.yaml"
    base = "genome: /tmp/genome.fa\ngtf: /tmp/genes.gtf\n"
    for body in ("oracle_bam: false\n", "emit_fastq: false\n"):
        config_path.write_text(base + body)
        parse_yaml_config(config_path)
    config_path.write_text(base + "oracle_bam: false\nemit_fastq: false\n")
    with pytest.raises(ValueError, match="oracle_bam and emit_fastq are both false"):
        parse_yaml_config(config_path)


def _paths_named(cfg) -> list[str]:
    """Every path a parsed config names: the five top-level paths, the abundance file, the probes."""
    top = [cfg.genome, cfg.gtf, cfg.shadow_gtf, cfg.index, cfg.outdir, cfg.abundance.file]
    probes = [cfg.capture.probes, *(scenario.config.probes for scenario in cfg.capture_configs)]
    return [path for path in [*top, *probes] if path]


def test_a_config_path_expands_its_variables_and_the_home_directory(tmp_path, monkeypatch):
    """The in-repo configs name their data root as `$RIGEL_RUNS`, so every path a config names goes
    through one expansion, `$NAME`, `${NAME}` and `~`, in each of the seven places a path is read.
    PERTURBATIONS: dropping `expandvars` fails every `$` path, dropping `expanduser` the two `~`
    paths, and a call site that skips the expansion fails its own key."""
    root, home = tmp_path / "root", tmp_path / "home"
    monkeypatch.setenv("RIGEL_TEST_ROOT", str(root))
    monkeypatch.setenv("HOME", str(home))
    config_path = tmp_path / "sim.yaml"
    config_path.write_text(
        "genome: $RIGEL_TEST_ROOT/genome.fa\n"
        "gtf: ${RIGEL_TEST_ROOT}/genes.gtf\n"
        "shadow_gtf: ~/shadow.gtf\n"
        "index: $RIGEL_TEST_ROOT/idx\n"
        "outdir: ~/out\n"
        "abundance:\n  mode: file\n  file: $RIGEL_TEST_ROOT/abundances.tsv\n"
        "capture:\n  configs:\n"
        '    - {label: "off", enabled: false}\n'
        '    - {label: "on", probes: $RIGEL_TEST_ROOT/probes.bed}\n'
    )

    cfg = parse_yaml_config(config_path)

    assert cfg.genome == str(root / "genome.fa")
    assert cfg.gtf == str(root / "genes.gtf")
    assert cfg.shadow_gtf == str(home / "shadow.gtf")
    assert cfg.index == str(root / "idx")
    assert cfg.outdir == str(home / "out")
    assert cfg.abundance.file == str(root / "abundances.tsv")
    assert cfg.capture_configs[1].config.probes == str(root / "probes.bed")


def test_a_config_path_naming_an_unset_or_empty_variable_is_refused_by_name(tmp_path, monkeypatch):
    """`os.path.expandvars` leaves an unset variable as written, so an unset `RIGEL_RUNS` would read and
    write a directory literally named `$RIGEL_RUNS`, and an empty one would move the path to `/`. Both
    are refused, naming the key and the variable. PERTURBATIONS: deleting the refusal fails every row;
    testing that the name is in the environment rather than that its value is not empty fails the
    empty one."""
    monkeypatch.delenv("RIGEL_TEST_UNSET", raising=False)
    monkeypatch.setenv("RIGEL_TEST_EMPTY", "")
    config_path = tmp_path / "sim.yaml"
    genome, gtf = "genome: /tmp/genome.fa\n", "gtf: /tmp/genes.gtf\n"
    probes = 'capture:\n  configs:\n    - {label: "on", probes: $RIGEL_TEST_UNSET/p.bed}\n'
    for body, key, variable in [
        ("genome: $RIGEL_TEST_UNSET/genome.fa\n" + gtf, "genome", "RIGEL_TEST_UNSET"),
        (genome + "gtf: ${RIGEL_TEST_UNSET}/genes.gtf\n", "gtf", "RIGEL_TEST_UNSET"),
        (genome + gtf + "outdir: $RIGEL_TEST_EMPTY/out\n", "outdir", "RIGEL_TEST_EMPTY"),
        (genome + gtf + probes, "capture config 'on': probes", "RIGEL_TEST_UNSET"),
    ]:
        config_path.write_text(body)
        with pytest.raises(ValueError, match=rf"^{re.escape(key)}: .*\(\${variable}\)"):
            parse_yaml_config(config_path)


def test_every_in_repo_sim_config_parses(monkeypatch):
    """A key the parser refuses would stop a panel from being re-simulated. The configs name their
    paths under `$RIGEL_RUNS`, set here as `panel.py` and `simulate_reads.py` set it, to the data root
    the instruments default to (`RUNS` in `scripts/design/_shared.py`), and every path such a config
    names lands under it. PERTURBATIONS: without the variable every such config is refused; a path in
    one written under another root fails the last assertion."""
    monkeypatch.syspath_prepend(str(ROOT / "scripts" / "design"))
    runs = importlib.import_module("_shared").RUNS
    monkeypatch.setenv("RIGEL_RUNS", str(runs))
    # the template's `outdir` names `$USER`; set it for a shell that has none
    monkeypatch.setenv("USER", os.environ.get("USER") or "rigel")
    configs = sorted((ROOT / "scripts" / "sim" / "configs").glob("*.yaml"))
    assert configs, "no sim configs found — this test would pass vacuously"
    for path in configs:
        cfg = parse_yaml_config(path)
        if "$RIGEL_RUNS" in path.read_text():
            outside = [p for p in _paths_named(cfg) if not Path(p).is_relative_to(runs)]
            assert not outside, f"{path.name} names a path outside $RIGEL_RUNS: {outside}"


def _sparse_population(n: int, *, seed: int = 3) -> tuple[list[Transcript], list[Transcript]]:
    """``n`` expressed multi-exon transcripts, one entity each, so an entity's ``nrna_abundance`` reads
    back exactly one draw. Mature abundances span the ladder's own four decades (log-uniform 1 →
    10,000), which is what lets the independence gate below see a correlation if there is one."""
    rng = np.random.default_rng(seed)
    contributors, entities = [], []
    for i in range(n):
        base = i * 10_000
        t = _transcript(
            f"t{i}",
            float(10.0 ** rng.uniform(0, 4)),
            [(base, base + 100), (base + 200, base + 300)],
        )
        t.t_index = i
        t.span_t_index = n + i
        contributors.append(t)
        entities.append(_entity(f"N{i}", base, base + 300, n + i))
    return contributors, entities


def test_the_on_fraction_leaves_most_gene_spans_at_EXACTLY_zero():
    """Sparsity is the model, and it is only observable below full eligibility: at
    ``on_fraction = 1.0`` the Bernoulli draw cannot be seen at all.

    Nascent RNA is absent from most gene spans and present in a minority, so ``on_fraction = 0.1``
    must leave about 90 % of eligible spans at nascent EXACTLY zero — not small, zero — while the
    rest carry a real level. The band is the binomial's own five sigma, derived rather than tuned.
    """
    n, frac = 400, 0.1
    contributors, entities = _sparse_population(n)
    apply_sparse_nrna(contributors + entities, (1.0, 1000.0), on_fraction=frac, seed=11)

    on = [e for e in entities if e.nrna_abundance > 0.0]
    off = [e for e in entities if e.nrna_abundance == 0.0]
    assert len(on) + len(off) == n, "a span is either on or exactly zero"
    assert on and off, "zero-inflation is not observable — the Bernoulli draw did nothing"
    sd = (n * frac * (1 - frac)) ** 0.5
    lo, hi = n * frac - 5 * sd, n * frac + 5 * sd
    assert lo <= len(on) <= hi, (
        f"{len(on)}/{n} spans on at on_fraction={frac}: outside the binomial's own "
        f"[{lo:.0f}, {hi:.0f}] five-sigma band, so the draw is not Bernoulli({frac})"
    )


def test_on_fraction_zero_transcribes_NOTHING_and_one_transcribes_EVERY_ELIGIBLE_span():
    """The two closed ends, exactly — a fraction is a fraction only if its endpoints are honoured."""
    contributors, entities = _sparse_population(40)
    assert apply_sparse_nrna(contributors + entities, (1.0, 1000.0), on_fraction=0.0, seed=5) == 0.0
    assert all(e.nrna_abundance == 0.0 for e in entities)

    contributors, entities = _sparse_population(40)
    apply_sparse_nrna(contributors + entities, (1.0, 1000.0), on_fraction=1.0, seed=5)
    assert all(e.nrna_abundance > 0.0 for e in entities)


def test_the_level_is_LOG_uniform_over_the_range_not_linear():
    """Where nascent is present its level spans decades — very low in some spans, high in others.
    A LINEAR draw on (1, 1000) puts 90 % of its mass in the top decade and cannot express "very low";
    a log-uniform draw spreads the ORDERS OF MAGNITUDE evenly, so the decade occupancies are
    near-equal. That is the property gated here, and it is what distinguishes the two draws.
    """
    lo, hi = 1.0, 1000.0
    contributors, entities = _sparse_population(600)
    apply_sparse_nrna(contributors + entities, (lo, hi), on_fraction=1.0, seed=7)
    levels = np.array([e.nrna_abundance for e in entities])

    assert levels.min() >= lo and levels.max() <= hi, "levels escaped the range"
    # three decades, so a log-uniform draw puts ~1/3 of the mass in each; a linear draw puts ~0.9 in
    # the last one. The band is generous — it only has to separate 1/3 from 9/10.
    per_decade = [np.mean((levels >= 10.0**k) & (levels < 10.0 ** (k + 1))) for k in range(3)]
    assert max(per_decade) < 0.5, (
        f"decade occupancies {[round(x, 3) for x in per_decade]} — the level is concentrated in one "
        f"decade, which is a linear draw, not a log-uniform one"
    )
    assert min(per_decade) > 0.15, f"a decade is nearly empty: {[round(x, 3) for x in per_decade]}"


def test_the_nascent_LEVEL_IS_INDEPENDENT_OF_THE_MATURE_LEVEL():
    """The decision this mode exists for: the nascent level is drawn INDEPENDENTLY of the mature
    level, so ``nascent > mature`` is a real case the tool has to survive. A ratio model setting
    ``nascent = mature x ratio`` makes the two perfectly rank-correlated and puts a ceiling under
    exactly that case.

    Two things are asserted: the rank correlation between the mature and nascent levels is ~0 (a
    ratio model reads 1.0), and spans where nascent EXCEEDS mature actually occur.
    """
    contributors, entities = _sparse_population(500, seed=9)
    apply_sparse_nrna(contributors + entities, (1.0, 1000.0), on_fraction=1.0, seed=23)

    mature = np.array([c.abundance for c in contributors])
    nascent = np.array([e.nrna_abundance for e in entities])
    rank_m = np.argsort(np.argsort(mature)).astype(float)
    rank_n = np.argsort(np.argsort(nascent)).astype(float)
    rho = float(np.corrcoef(rank_m, rank_n)[0, 1])
    # 500 independent pairs: sd of Spearman's rho under the null is ~1/sqrt(n-1) = 0.045, so five
    # sigma is 0.22 — derived from the sample size, not tuned
    assert abs(rho) < 5.0 / (len(mature) - 1) ** 0.5, (
        f"rank correlation between mature and nascent levels is {rho:.3f} — the nascent level is "
        f"tracking the mature one, which is the ratio model this mode replaced"
    )
    assert (nascent > mature).mean() > 0.1, (
        f"only {(nascent > mature).mean():.1%} of spans carry more nascent than mature — the "
        f"robustness case the owner asked for is not being exercised"
    )


def test_the_on_off_DRAW_IS_PER_SPAN_NOT_PER_ISOFORM():
    """The unit of sparsity, gated where it is actually observable. With ``k`` isoforms sharing one
    gene span, a per-CONTRIBUTOR draw gives that span a ``1 - (1 - p)^k`` chance of carrying nascent
    — at ``p = 0.1`` and ``k = 5`` that is 41 %, four times the configured sparsity — while a
    per-SPAN draw gives exactly ``p``. The intron slots calibration reads belong to the span, so the
    per-span reading is the one that makes the configured number mean what it says.

    PERTURBATION: switching to per-contributor drawing passes every other test in this file, because
    the shared-span test runs at ``on_fraction = 1.0`` where the two are indistinguishable
    (`TRAPS: perturb-every-gate`).
    """
    n_spans, k, p = 300, 5, 0.1
    rng = np.random.default_rng(4)
    rows: list[Transcript] = []
    entities: list[Transcript] = []
    for sp in range(n_spans):
        base = sp * 100_000
        entity = _entity(f"N{sp}", base, base + 50_000, 10_000 + sp)
        entities.append(entity)
        for j in range(k):
            off = base + j * 1_000
            t = _transcript(
                f"t{sp}_{j}",
                float(10.0 ** rng.uniform(0, 4)),
                [(off, off + 100), (off + 200, off + 300)],
            )
            t.t_index = sp * k + j
            t.span_t_index = 10_000 + sp
            rows.append(t)
    apply_sparse_nrna(rows + entities, (1.0, 1000.0), on_fraction=p, seed=31)

    on_frac = float(np.mean([e.nrna_abundance > 0.0 for e in entities]))
    per_isoform = 1.0 - (1.0 - p) ** k  # 0.41 — what a per-contributor draw would give
    sd = (p * (1 - p) / n_spans) ** 0.5
    assert abs(on_frac - p) < 5 * sd, (
        f"{on_frac:.3f} of spans on at on_fraction={p} with {k} isoforms each: the per-span value is "
        f"{p}, a per-isoform draw would give {per_isoform:.3f} — the draw is not per span"
    )


def test_the_SEED_actually_drives_the_draw():
    """PERTURBATION: every other gate here is statistical, so replacing `default_rng(seed)` with a
    constant seed passes all of them. The orchestrator adds `seed + nrna_index` precisely so a
    multi-range sweep gets DIFFERENT on-sets, which would otherwise collapse to one pattern silently
    (`TRAPS: perturb-every-gate`).

    Two seeds must disagree, and one seed must reproduce itself exactly.
    """

    def draw(seed):
        contributors, entities = _sparse_population(200)
        apply_sparse_nrna(contributors + entities, (1.0, 100.0), on_fraction=0.5, seed=seed)
        return np.array([e.nrna_abundance for e in entities])

    a, b, a_again = draw(1), draw(2), draw(1)
    assert np.array_equal(a, a_again), "the same seed did not reproduce its own draw"
    assert not np.array_equal(a, b), (
        "two different seeds produced identical draws — the seed is inert"
    )
    assert (a > 0).sum() and (b > 0).sum(), "a draw with nothing on cannot distinguish anything"


def test_a_multi_exon_contributor_with_no_ENTITY_raises_in_the_SPARSE_path_too():
    """PERTURBATION: with the rule gated only through `apply_nrna_ratio`, replacing the sparse path's
    own raise with `continue` passes every test. A transcript list without entities is not a rigel
    index's, and dropping its molecules silently is the failure the rule exists to stop."""
    orphan = _transcript("orphan", 100.0, [(0, 100), (200, 300)])
    orphan.t_index = 0  # span_t_index stays unset/-1: no entity
    with pytest.raises(ValueError, match="nascent entity"):
        apply_sparse_nrna([orphan], (1.0, 100.0), on_fraction=1.0, seed=1)


def test_the_FUNCTION_validates_its_own_range_and_fraction_not_just_the_PARSER():
    """The parser is not the only door. `apply_sparse_nrna` is public, and a caller that builds its
    arguments directly never passes through `parse_yaml_config`, so the function's own guards are all
    that stands between a range of `(0, 100)` and `log(0) = -inf` inside the draw. PERTURBATION:
    deleting either guard passes every other gate in this file."""
    rows = sum(_sparse_population(4), [])
    for bad, why in [
        ((0.0, 100.0), "zero lower end has no log"),
        ((-1.0, 10.0), "negative lower end"),
        ((100.0, 1.0), "inverted endpoints"),
    ]:
        with pytest.raises(ValueError, match="0 < lo <= hi|abundance_range"):
            apply_sparse_nrna(list(rows), bad, on_fraction=0.5, seed=1)
    for bad in (1.5, -0.1):
        with pytest.raises(ValueError, match="on_fraction"):
            apply_sparse_nrna(list(rows), (1.0, 10.0), on_fraction=bad, seed=1)


# ── a run: the index the transcriptome comes from, and the outputs it keeps ─────────────────────────

_TINY_GTF = (
    'chr1\tsrc\texon\t101\t400\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
    'chr1\tsrc\texon\t501\t800\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
    'chr1\tsrc\texon\t1001\t1800\t.\t-\t.\tgene_id "G2"; transcript_id "T2";\n'
)


@pytest.fixture(scope="module")
def tiny_reference(tmp_path_factory):
    """A 2 kb genome, a two-gene GTF and the rigel index built from them: enough for `run_simulation`
    to simulate one small condition."""
    import pysam

    from rigel.index import TranscriptIndex

    base = tmp_path_factory.mktemp("tiny_reference")
    genome = base / "genome.fa"
    rng = np.random.default_rng(0)
    genome.write_text(">chr1\n" + "".join(rng.choice(list("ACGT"), 2000)) + "\n")
    pysam.faidx(str(genome))
    gtf = base / "genes.gtf"
    gtf.write_text(_TINY_GTF)
    index = base / "index"
    TranscriptIndex.build(str(genome), str(gtf), str(index), write_tsv=False)
    return genome, gtf, index


def _tiny_config(path, genome, gtf, index, outdir, extra=""):
    """One condition of 40 RNA fragments, every transcript expressed."""
    path.write_text(
        f"genome: {genome}\ngtf: {gtf}\nindex: {index}\noutdir: {outdir}\n"
        "simulation:\n  n_rna_fragments: 40\n  frag_mean: 150\n  frag_std: 20\n"
        "  frag_min: 100\n  frag_max: 200\n  read_length: 50\n"
        "abundance:\n  frac_expressed: 1.0\n  min: 1.0\n  max: 10.0\n" + extra
    )
    return path


def test_the_manifest_records_the_index_the_transcriptome_comes_from(tiny_reference, tmp_path):
    """The simulated transcriptome is the index's, so the manifest names the index: its resolved path and
    its own manifest, which carries the sha256 of the FASTA and the GTF it was built from. The config's
    `gtf:` here is a byte-identical COPY at another path; the check is on content, so the run proceeds
    and the record still names the file the index was built from. PERTURBATIONS: comparing paths instead
    of digests refuses this run; dropping the record fails the equality."""
    genome, gtf, index = tiny_reference
    copy = tmp_path / "copy_of_genes.gtf"
    copy.write_bytes(gtf.read_bytes())
    outdir = tmp_path / "out"
    cfg = parse_yaml_config(_tiny_config(tmp_path / "sim.yaml", genome, copy, index, outdir))

    run_simulation(cfg)

    recorded = json.loads((outdir / "manifest.json").read_text())["index"]
    assert recorded == {
        "path": str(index.resolve()),
        "manifest": json.loads((index / "manifest.json").read_text()),
    }
    sources = recorded["manifest"]["sources"]
    assert sources["gtf"]["path"] == str(gtf.resolve())
    assert sources["gtf"]["sha256"] == hashlib.sha256(gtf.read_bytes()).hexdigest()
    assert sources["fasta"]["sha256"] == hashlib.sha256(genome.read_bytes()).hexdigest()


def test_a_gtf_the_index_was_not_built_from_is_refused(tiny_reference, tmp_path):
    """`run_simulation` reads its transcriptome from `index:` and never parses `gtf:`, so a config naming
    another GTF simulated the index's transcripts while its manifest named the other file. A one-byte
    edit is refused, before any condition is simulated. PERTURBATION: deleting the digest comparison
    simulates this config."""
    genome, gtf, index = tiny_reference
    edited = tmp_path / "edited.gtf"
    edited.write_text(gtf.read_text().replace("\t1800\t", "\t1799\t"))
    outdir = tmp_path / "out"
    cfg = parse_yaml_config(_tiny_config(tmp_path / "sim.yaml", genome, edited, index, outdir))

    with pytest.raises(ValueError, match="is not the GTF the index"):
        run_simulation(cfg)
    assert not list(outdir.glob("gdna_*")) and not (outdir / "manifest.json").exists()


def test_no_oracle_on_a_config_without_fastqs_is_refused_before_anything_is_simulated(
    tiny_reference, tmp_path, monkeypatch
):
    """`simulate_reads.py --no-oracle` turns the oracle BAM off AFTER the config is read, so the parser
    cannot see it; `run_simulation` refuses it before writing anything. The config alone is valid: it is
    the fl-gap panels' `emit_fastq: false`. PERTURBATION: without the run-time check this run exits 0
    having simulated the condition and deleted its only reads."""
    genome, gtf, index = tiny_reference
    outdir = tmp_path / "out"
    config = _tiny_config(tmp_path / "sim.yaml", genome, gtf, index, outdir, "emit_fastq: false\n")
    spec = importlib.util.spec_from_file_location(
        "_simulate_reads", ROOT / "scripts" / "sim" / "simulate_reads.py"
    )
    script = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(script)
    monkeypatch.setattr(sys, "argv", ["simulate_reads.py", "--config", str(config), "--no-oracle"])

    with pytest.raises(ValueError, match="oracle_bam and emit_fastq are both false"):
        script.main()
    assert not outdir.exists()
