"""G-S3 / G-S4 / G-S5 — hybrid capture SELECTS FOR LENGTH, and it must do so in the marginal.

    `docs/TESTING.md` §1, §3

⭐ **The principle, and it is the owner's.** Hybrid capture works by hybridisation of probes to
sequence. A short fragment presents less sequence to hybridise, binds less well, and is captured less
efficiently. **Capture therefore selects for longer fragments**, for gDNA and RNA alike. Two
consequences that are not optional:

1. the post-capture distributions are the true baseline — the pre-capture ones describe a library
   that was never sequenced;
2. capture **narrows** the gDNA<->RNA length gap whenever gDNA is the shorter component, because the
   short tail it removes is disproportionately gDNA.

⛔ **The defect these gates exist for.** The engine drew the length marginal FIRST, capture-blind,
then computed a capture-aware per-template opportunity ``total_eff(w)`` and discarded it by
normalising within each length. So capture could move only *where* a fragment landed, never *whether*
its length survived — and the panel's post-capture truth summaries were byte-identical to the
pre-capture ones (gDNA mean 195.57 on both arms). The fix is to keep the term:
``f_post(w) ∝ f_pre(w) · total_eff(w)``.

⚠ **G-S3 and G-S5 are DIRECTIONAL, with no threshold.** Choosing a magnitude would be inventing the
capture efficiency curve, which is `docs/TRAPS.md` no-magic-numbers. The sign is physics; the
magnitude is whatever ``binding_per_base`` and the probe length imply.

⚠ **G-S5 is conditional on gDNA being the shorter, broader component**, which is what this fixture
and the pilot panel both are (gDNA 157 +/- 125, RNA 206 +/- 98, measured off a real cfRNA library).
Never assume that ordering in general — for other library types it reverses.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pysam
import pytest

from rigel.sim.capture import CaptureConfig, CaptureSampler
from rigel.sim.genome import random_dna_array
from rigel.sim.wgs_config import GDNASimConfig, SimulationParams
from rigel.sim.wgs_engine import WholeGenomeSimulator
from rigel.transcript import Transcript
from rigel.types import Interval, Strand


GENOMIC_REFS = ["chrA", "chrB"]

#: The pilot panel's own fragment-length parameters — the count-weighted mean and sd of a real cfRNA
#: library's pools (MO_3021), not chosen round numbers.
RNA_LENGTH = dict(frag_mean=206.0, frag_std=98.0, frag_min=50, frag_max=800)
GDNA_LENGTH = dict(frag_mean=157.0, frag_std=125.0, frag_min=40, frag_max=1000)


def _write_fasta(tmp_path: Path, sequences: dict[str, str]) -> Path:
    path = tmp_path / "genome.fa"
    with open(path, "w") as handle:
        for name, sequence in sequences.items():
            handle.write(f">{name}\n")
            for i in range(0, len(sequence), 60):
                handle.write(sequence[i : i + 60] + "\n")
    pysam.faidx(str(path))
    return path


def _transcript(t_id, ref, exons, abundance) -> Transcript:
    transcript = Transcript(
        ref=ref,
        strand=Strand.POS,
        exons=[Interval(start, end) for start, end in exons],
        t_id=t_id,
        g_id=f"g_{t_id}",
        g_name=f"g_{t_id}",
        abundance=abundance,
    )
    transcript.compute_length()
    return transcript


@pytest.fixture(scope="module")
def panel(tmp_path_factory):
    """Two 200 kb genomic references carrying 24 transcripts, half of them on a probe panel.

    Sized so the geometry is the pilot's rather than a toy's: a probe run is a small neighbourhood
    inside a large template, which is the regime in which capture's length preference is a real
    reweighting instead of a boundary effect.
    """
    tmp_path = tmp_path_factory.mktemp("capture_length")
    rng = np.random.default_rng(11)
    lengths = {"chrA": 200_000, "chrB": 200_000}
    sequences = {name: "".join(random_dna_array(n, rng)) for name, n in lengths.items()}
    fasta = _write_fasta(tmp_path, sequences)

    transcripts = []
    for ref_index, ref in enumerate(GENOMIC_REFS):
        for i in range(12):
            start = 5_000 + i * 15_000
            # Two exons, so the transcript is 4,000 nt of spliced sequence over a 9,000 bp span.
            exons = [(start, start + 2_000), (start + 7_000, start + 9_000)]
            transcripts.append(
                _transcript(f"T{ref_index}_{i}", ref, exons, abundance=100.0 * (i + 1))
            )

    # Probes on every other transcript, one 120 bp probe well inside the spliced sequence.
    probes = tmp_path / "probes.tsv"
    rows = ["transcript_id\tstart\tend"]
    for index, transcript in enumerate(transcripts):
        if index % 2 == 0:
            rows.append(f"{transcript.t_id}\t1500\t1620")
    probes.write_text("\n".join(rows) + "\n")
    return fasta, transcripts, str(probes)


def _simulator(fasta, transcripts, probes, *, seed=23):
    return WholeGenomeSimulator(
        fasta,
        transcripts,
        SimulationParams(sim_seed=seed, read_length=100, **RNA_LENGTH),
        GDNASimConfig(**GDNA_LENGTH),
        seed=seed,
        genomic_refs=GENOMIC_REFS,
        capture_config=(
            CaptureConfig(probes=probes, probe_format="transcript", binding_per_base=10.0)
            if probes is not None
            else None
        ),
    )


def _mean(length_counts: dict[int, int]) -> float:
    n = sum(length_counts.values())
    return sum(length * count for length, count in length_counts.items()) / n


def _gdna_lengths(simulator, n: int) -> dict[int, int]:
    counts: dict[int, int] = {}
    for (_ref_index, length), count in simulator._accumulate_gdna_counts(n).items():
        counts[length] = counts.get(length, 0) + count
    return counts


def _mrna_lengths(simulator, n: int) -> dict[int, int]:
    counts: dict[int, int] = {}
    mrna, _nrna = simulator._accumulate_rna_counts(n)
    for per_length in mrna.values():
        for length, count in per_length.items():
            counts[length] = counts.get(length, 0) + count
    return counts


@pytest.fixture(scope="module")
def arms(panel):
    """``{"off"|"on": {"gdna": mean, "rna": mean}}`` — capture is the ONLY thing varied.

    ⚠ **Each pool gets its own freshly seeded simulator, and that is what makes this an A/B.** A
    simulator's first use of its RNG is the pre-capture length draw, so a fresh instance per pool
    means both arms draw *the same* pre-capture lengths and every difference downstream is
    attributable to capture. Measuring both pools from one instance does not: the gDNA pool consumes
    a different amount of RNG in each arm, so the RNA draw that follows it differs between the arms
    by accident. That was measured — with the reweighting deliberately removed from the RNA pool, the
    RNA gate still read a rise and said nothing.
    """
    fasta, transcripts, probes = panel
    n = 20_000
    result = {}
    for label, probe_file in (("off", None), ("on", probes)):
        gdna_sim = _simulator(fasta, transcripts, probe_file)
        try:
            gdna = _gdna_lengths(gdna_sim, n)
        finally:
            gdna_sim.close()
        rna_sim = _simulator(fasta, transcripts, probe_file)
        try:
            rna = _mrna_lengths(rna_sim, n)
        finally:
            rna_sim.close()
        result[label] = {
            "gdna": _mean(gdna),
            "rna": _mean(rna),
            "n_gdna": sum(gdna.values()),
            "n_rna": sum(rna.values()),
        }
    return result


class TestCaptureMovesTheLengthMarginal:
    def test_the_arms_are_the_same_size(self, arms):
        """Every requested fragment is still allocated — the reweighting is not a filter."""
        for label, values in arms.items():
            assert values["n_gdna"] == 20_000, label
            assert values["n_rna"] == 20_000, label

    def test_g_s3_gdna_mean_length_rises_under_capture(self, arms):
        """⭐ G-S3 — strictly greater. Byte-identical arms are the falsification."""
        off, on = arms["off"]["gdna"], arms["on"]["gdna"]
        assert on > off, f"gDNA mean length {off:.2f} (off) -> {on:.2f} (on)"

    def test_capture_also_lengthens_the_rna_pool(self, arms):
        """The principle applies to gDNA and RNA alike — a probe hybridises to either."""
        off, on = arms["off"]["rna"], arms["on"]["rna"]
        assert on > off, f"RNA mean length {off:.2f} (off) -> {on:.2f} (on)"

    def test_g_s5_the_gdna_to_rna_length_gap_narrows_under_capture(self, arms):
        """⭐ G-S5 — the quantity the solver actually consumes. `mu_g - mu_r` is the ONLY thing that
        identifies the fragment-length channel: at equal component means it carries exactly zero
        information about composition at any depth."""
        gap_off = abs(arms["off"]["gdna"] - arms["off"]["rna"])
        gap_on = abs(arms["on"]["gdna"] - arms["on"]["rna"])
        assert gap_on < gap_off, f"|mu_g - mu_r| {gap_off:.2f} (off) -> {gap_on:.2f} (on)"

    def test_the_short_tail_is_what_capture_removes(self, panel):
        """The mechanism, not just its sign: the deficit is concentrated at the short end."""
        fasta, transcripts, probes = panel
        n = 20_000
        histograms = {}
        for label, probe_file in (("off", None), ("on", probes)):
            simulator = _simulator(fasta, transcripts, probe_file)
            try:
                histograms[label] = _gdna_lengths(simulator, n)
            finally:
                simulator.close()

        def share_below(counts, limit):
            total = sum(counts.values())
            return sum(c for length, c in counts.items() if length < limit) / total

        assert share_below(histograms["on"], 100) < share_below(histograms["off"], 100)


class TestTheConditionalWasAlreadyRight:
    """⚠ **A REGRESSION GUARD, NOT A FALSIFICATION — and the distinction is the finding.**

    "On-target gDNA is longer than off-target" already held BEFORE capture was made length-selective,
    because the engine's *conditional* was correct all along: at a fixed width, the probability of landing
    on a probe rises with the width. What was missing was the *marginal* — the engine normalised
    ``total_eff(w)`` away, so the length distribution never moved. ⭐ So the marginal gate above is the
    falsification and this one is a guard on the piece that already worked (`docs/TRAPS.md` a-gate-that-already-passed, TRAPS: right-conditional-wrong-marginal).

    ⛔ **And an on-target population defined by the START's territory cannot be satisfied by any capture
    model of this form.** Measured: on-target (probe-overlapping) 197.30 vs off-target 155.05 — correct —
    while by the territory the START lands in, exonic reads 190.93 against intronic 237.63. That inversion
    is geometry: conditioned on being captured, a fragment whose start is in the intron is one that was
    long enough to *reach* the probe, so an intronic start carries weight ~ w^2/2 while an exonic start
    carries ~ p^2/2, flat in w. `docs/TRAPS.md` on-target-by-start-is-geometry, and the sibling of TRAPS: a-purity-filter-is-a-length-filter and TRAPS: pure-and-length-censored — ask what a pool's
    selection rule correlates with, not only what the selected fragments are.
    """

    def test_the_on_probe_share_of_capture_mass_rises_with_fragment_length(self, tmp_path):
        """⭐ Checked against the closed form, not against a recorded number.

        For a probe ``[b, e)`` sitting far enough inside a template that no start is truncated,
        ``sum_s overlap(s, w) = (e - b) * w`` exactly — every probe base is covered by exactly ``w``
        of the admissible starts. So the capture-weighted opportunity is
        ``off_target * (L - w + 1) + binding_per_base * p * w``, whose probe-attributable share is
        strictly increasing in ``w``. That is the whole of "capture selects for length".
        """
        transcript_length, probe_start, probe_length = 20_000, 8_000, 120
        transcript = _transcript("T", "chrA", [(0, transcript_length)], 100.0)
        probes = tmp_path / "probes.tsv"
        probes.write_text(
            f"transcript_id\tstart\tend\nT\t{probe_start}\t{probe_start + probe_length}\n"
        )
        config = CaptureConfig(
            probes=str(probes),
            probe_format="transcript",
            off_target_weight=1.0,
            binding_per_base=10.0,
        )
        sampler = CaptureSampler.from_config(config, [transcript], {"chrA": transcript_length})

        shares = []
        for width in (60, 120, 240, 480, 800):
            total = sampler.partition("mrna", 0, transcript_length, width)
            baseline = config.off_target_weight * (transcript_length - width + 1)
            expected = baseline + config.binding_per_base * probe_length * width
            assert total == pytest.approx(expected), f"width {width}"
            shares.append((total - baseline) / total)

        assert shares == sorted(shares)
        assert shares[0] < shares[-1]
