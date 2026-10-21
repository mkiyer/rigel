"""The simulator draws no fragment over sequence the reference does not know.

A reference's `N` runs are assembly gaps: no molecule reads as `N`, and an aligner places no read there. The
simulator once drew gDNA over them anyway — about a fifth of the ladder's gDNA off capture — which no aligner could
place and the oracle BAM placed by truth, so the oracle ladder agreed with calibration's region lengths, which count
those bases too, and hid that defect. The capture sampler now masks every unknown base per template: a fragment of
gDNA or RNA that would overlap one is never drawn and adds no opportunity, in the allocation across templates and in
the placement within one, with capture and without. These gates hold the masks to the templates, the masked totals to
brute-force enumeration over every start, the masked draws to the masked weights, and the simulated reads to known
sequence.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
import pysam
import pytest

from rigel.sim.capture import CaptureConfig, CaptureSampler
from rigel.sim.capture.sampler import ProbeInterval, unknown_runs
from rigel.sim.genome import MutableGenome
from rigel.sim.read_name import parse_origin
from rigel.sim.wgs_config import GDNASimConfig, SimulationParams
from rigel.sim.wgs_engine import WholeGenomeSimulator
from rigel.transcript import Transcript
from rigel.types import Interval, Strand


def _transcript(t_id, exons, strand=Strand.POS, ref="chr1", **kw):
    t = Transcript(
        ref=ref,
        strand=strand,
        exons=[Interval(s, e) for s, e in exons],
        t_id=t_id,
        g_id=f"g_{t_id}",
        **kw,
    )
    t.length = t.compute_length()
    return t


def _runs_of(seq: str) -> list[tuple[int, int]]:
    """The `N` runs of a string, written out the slow way."""
    runs, start = [], None
    for i, c in enumerate(seq + "A"):
        if c == "N" and start is None:
            start = i
        elif c != "N" and start is not None:
            runs.append((start, i))
            start = None
    return runs


def _revcomp(seq: str) -> str:
    return seq[::-1].translate(str.maketrans("ACGTN", "TGCAN"))


# ── the masks ────────────────────────────────────────────────────────────────────────────────────


def test_unknown_runs_are_the_N_runs_of_a_sequence():
    seq = "ACNNNGTNACGTNN"
    got = unknown_runs(np.frombuffer(seq.encode(), dtype=np.uint8))
    assert got.tolist() == [list(r) for r in _runs_of(seq)]
    assert unknown_runs(np.frombuffer(b"ACGT", dtype=np.uint8)).shape == (0, 2)


def test_every_template_is_masked_at_its_own_unknown_bases(tmp_path):
    """The mask sits in each template's own coordinates: a gDNA reference's, a spliced transcript's (an N run in
    an intron is no base of it, one in an exon is), and the reverse strand's, read 5′ to 3′ as the engine builds it.

    Perturbation: masking in genomic coordinates, or ignoring the strand, misplaces every mRNA run."""
    genome = MutableGenome(3000, seed=1, name="chr1")
    genome.edit(150, "N" * 20)  # intergenic
    genome.edit(620, "N" * 7)  # in t_pos's first exon
    genome.edit(1300, "N" * 30)  # in an intron of both, and inside the nascent span
    genome.edit(2210, "N" * 5)  # in t_neg's last exon
    fasta_path = genome.write_fasta(tmp_path)
    t_pos = _transcript("t_pos", [(500, 800), (1500, 1700)])
    t_neg = _transcript("t_neg", [(1000, 1200), (2000, 2300)], strand=Strand.NEG)
    span = _transcript("span", [(500, 1700)], is_synthetic=True)
    sampler = CaptureSampler.disabled([t_pos, t_neg, span], {"chr1": len(genome)})
    with pysam.FastaFile(str(fasta_path)) as fa:
        sampler.mask_unknown_sequence(fa, ["chr1"])
        seq = fa.fetch("chr1")
    assert sampler.unknown("gdna", "chr1").tolist() == [list(r) for r in _runs_of(seq)]
    templates = {
        0: seq[500:800] + seq[1500:1700],
        1: _revcomp(seq[1000:1200] + seq[2000:2300]),
        2: seq[500:1700],
    }
    for k, template in templates.items():
        assert sampler.unknown("mrna", k).tolist() == [list(r) for r in _runs_of(template)], k


# ── the masked totals and draws, against enumeration ───────────────────────────────────────────────

SEQ_LEN = 1000
MASKS = {
    "none": [],
    "mid": [(400, 430)],
    "at the start": [(0, 15)],
    "at the end": [(980, 1000)],
    "two, closer than a fragment": [(300, 310), (350, 360)],
    "two, far apart": [(100, 120), (700, 760)],
}
LAYOUTS = {
    "no probes": [],
    "one probe": [ProbeInterval(380, 500)],
    "two probes": [ProbeInterval(50, 170), ProbeInterval(640, 780)],
}
WIDTHS = [1, 50, 120, 200, 999, 1000]


def _sampler(intervals, enabled=True, off=1.0, binding=10.0):
    config = CaptureConfig(
        probes="unused", probe_format="transcript", off_target_weight=off, binding_per_base=binding
    )
    sampler = CaptureSampler(config, transcripts=(), ref_lengths={}, enabled=enabled)
    sampler._mrna_intervals = {0: list(intervals)} if intervals else {}
    return sampler


def _weights(intervals, runs, width, enabled, off=1.0, binding=10.0):
    """Every start's weight, enumerated: zero where the fragment overlaps an unknown base."""
    out = np.zeros(max(0, SEQ_LEN - width + 1))
    for s in range(len(out)):
        if any(s < hi and s + width > lo for lo, hi in runs):
            continue
        best = 0
        if enabled:
            for p in intervals:
                best = max(best, min(s + width, p.end) - max(s, p.start))
        out[s] = off + binding * max(best, 0)
    return out


@pytest.mark.parametrize("enabled", [True, False])
@pytest.mark.parametrize("mask", list(MASKS))
@pytest.mark.parametrize("layout", list(LAYOUTS))
@pytest.mark.parametrize("width", WIDTHS)
def test_the_masked_partition_is_the_enumerated_sum_over_known_starts(enabled, mask, layout, width):
    """With capture or without, a template's opportunity is the sum of its starts' weights over the starts whose
    fragment covers no unknown base; the admissible count is their number.

    Perturbation: a mask that blocks only the run's own starts (not the `width − 1` before it), or one the binding
    term ignores, departs from the enumeration."""
    sampler = _sampler(LAYOUTS[layout], enabled=enabled)
    sampler.mask_unknown("mrna", 0, MASKS[mask])
    lengths = np.array([SEQ_LEN], dtype=np.int64)
    w = _weights(LAYOUTS[layout], MASKS[mask], width, enabled)
    assert sampler.partition_array("mrna", [0], lengths, width)[0] == pytest.approx(
        w.sum(), rel=1e-12, abs=1e-9
    )
    assert sampler.partition("mrna", 0, SEQ_LEN, width) == pytest.approx(
        w.sum(), rel=1e-12, abs=1e-9
    )
    assert int(sampler.admissible_counts("mrna", [0], lengths, width)[0]) == int((w > 0).sum())


@pytest.mark.parametrize("enabled", [True, False])
@pytest.mark.parametrize("mask", ["mid", "two, closer than a fragment", "at the start"])
def test_the_draws_avoid_unknown_bases_and_follow_the_masked_weights(enabled, mask):
    """No drawn fragment covers an unknown base, and the starts follow the masked weights: a chi-square over every
    admissible start, on a draw large enough to see a weight the mask got wrong.

    Perturbation: drawing the off-target share uniformly over every start, or reading the binding landscape
    unmasked, puts mass where the weights say none."""
    layout = LAYOUTS["two probes"]
    sampler = _sampler(layout, enabled=enabled)
    sampler.mask_unknown("mrna", 0, MASKS[mask])
    width, n = 120, 400_000
    starts = sampler.sample_starts("mrna", 0, SEQ_LEN, width, n, np.random.default_rng(5))
    w = _weights(layout, MASKS[mask], width, enabled)
    assert (w[starts] > 0).all(), "a draw covered an unknown base"
    observed = np.bincount(starts, minlength=len(w))
    expected = n * w / w.sum()
    live = expected > 0
    chi2 = float((((observed - expected) ** 2)[live] / expected[live]).sum())
    dof = int(live.sum()) - 1
    assert abs(chi2 - dof) < 6 * np.sqrt(2 * dof), (chi2, dof)
    for s in range(len(w)):  # every start, the blocked ones too, which no draw reaches
        assert sampler.fragment_weight("mrna", 0, SEQ_LEN, s, width) == pytest.approx(w[s]), s


def test_a_mask_is_kept_sorted_and_merged_whatever_order_it_arrives_in():
    """The blocked starts are read off the runs in order, so a mask given out of order, or overlapping, is put in
    order on install: given unsorted, one run once hid the other entirely. An empty or negative run is refused.

    Perturbation: storing the runs as given."""
    sampler = _sampler([])
    sampler.mask_unknown("mrna", 0, [(700, 760), (100, 120), (110, 130)])
    assert sampler.unknown("mrna", 0).tolist() == [[100, 130], [700, 760]]
    w = _weights([], [(100, 130), (700, 760)], 50, False)
    assert sampler.partition_array("mrna", [0], np.array([SEQ_LEN]), 50)[0] == pytest.approx(
        w.sum()
    )
    for bad in ([(5, 5)], [(-3, 4)]):
        with pytest.raises(ValueError):
            sampler.mask_unknown("mrna", 0, bad)


def test_fragments_asked_of_a_template_with_no_known_start_are_refused():
    """A template whose every start covers an unknown base has no opportunity, so the allocation gives it no
    fragment; asked for some anyway, the sampler refuses and names why rather than drawing over the gap."""
    sampler = _sampler([], enabled=False)
    sampler.mask_unknown("mrna", 0, [(100, 900)])
    with pytest.raises(ValueError, match="unknown base"):
        sampler.sample_starts("mrna", 0, SEQ_LEN, 950, 5, np.random.default_rng(0))


def test_an_unchanged_mask_keeps_the_memo_and_a_changed_one_clears_it():
    """A sampler is shared by every condition of a panel and each one installs the same masks: reinstalling an
    identical mask must keep the per-width memo (it is minutes of work), and a different one must not reuse it."""
    sampler = _sampler(LAYOUTS["one probe"])
    lengths = np.array([SEQ_LEN], dtype=np.int64)
    sampler.mask_unknown("mrna", 0, MASKS["mid"])
    before = sampler.partition_array("mrna", [0], lengths, 120)[0]
    sampler.mask_unknown("mrna", 0, MASKS["mid"])
    assert sampler._partition_memo, "an identical mask cleared the memo"
    sampler.mask_unknown("mrna", 0, MASKS["none"])
    after = sampler.partition_array("mrna", [0], lengths, 120)[0]
    assert after > before


# ── the reads ────────────────────────────────────────────────────────────────────────────────────


def test_no_simulated_read_covers_an_unknown_base_and_gdna_goes_where_sequence_is_known(tmp_path):
    """End to end: a reference half of which is unknown sequence, beside one with none. No read holds an `N`, no
    fragment's interval (its read name) overlaps an unknown run, and gDNA splits between the two references by their
    known opportunity, not their length.

    Perturbation: drawing gDNA over the whole reference writes all-`N` reads and gives the gapped reference twice
    its share."""
    a = MutableGenome(20_000, seed=2, name="chrA")
    a.edit(5_000, "N" * 10_000)
    a.edit(2_000, "N" * 3)  # inside the nascent span below
    b = MutableGenome(10_000, seed=3, name="chrB")
    fasta = tmp_path / "genome.fa"
    fasta.write_text(f">chrA\n{a.seq}\n>chrB\n{b.seq}\n")
    pysam.faidx(str(fasta))
    span = _transcript(
        "span", [(1_000, 3_000)], ref="chrA", is_synthetic=True, nrna_abundance=100.0
    )
    tx = _transcript("tx", [(16_000, 16_500), (17_000, 17_500)], ref="chrA", abundance=100.0)
    sim = WholeGenomeSimulator(
        fasta,
        [tx, span],
        SimulationParams(
            sim_seed=9, frag_mean=200, frag_std=20, frag_min=150, frag_max=250, read_length=50
        ),
        GDNASimConfig(frag_mean=200, frag_std=20, frag_min=150, frag_max=250),
        genomic_refs=["chrA", "chrB"],
    )
    try:
        r1, r2, _bam = sim.simulate_and_write(
            tmp_path / "out", n_rna=4_000, n_gdna=20_000, oracle_bam=True, n_workers=1
        )
    finally:
        sim.close()
    unknown = {"chrA": [(2_000, 2_003), (5_000, 15_000)], "chrB": []}
    span_start = 1_000
    per_ref = {"chrA": 0, "chrB": 0}
    for path in (r1, r2):
        lines = gzip.decompress(Path(path).read_bytes()).decode().splitlines()
        assert not any("N" in s for s in lines[1::4]), "a read holds an unknown base"
    for name in gzip.decompress(Path(r1).read_bytes()).decode().splitlines()[0::4]:
        o = parse_origin(name)
        if o.kind == "gdna":
            per_ref[o.ref] += 1
            lo, hi = o.start, o.end
            assert not any(lo < e and hi > s for s, e in unknown[o.ref]), name
        elif o.kind == "nrna":
            lo, hi = span_start + o.start, span_start + o.end
            assert not (lo < 2_003 and hi > 2_000), name
    # known opportunity: chrA's two known stretches against chrB's whole length, at the mean width
    known_a = (2_000 - 200 + 1) + (5_000 - 2_003 - 200 + 1) + (20_000 - 15_000 - 200 + 1)
    known_b = 10_000 - 200 + 1
    p = known_a / (known_a + known_b)
    n = per_ref["chrA"] + per_ref["chrB"]
    assert abs(per_ref["chrA"] - n * p) < 5 * np.sqrt(n * p * (1 - p)), (per_ref, p)
