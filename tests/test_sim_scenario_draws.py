"""What a scenario's rows ask for is what the simulator draws.

A transcript's own nascent abundance (``add_gene``'s ``nrna_abundance``) is independent of its mature one, so a
transcript with no mature molecules still has nascent ones. Each gate carries the control it is read against.
"""

from __future__ import annotations

from collections import Counter

import pysam

from rigel.sim import ReadSimConfig, Scenario
from rigel.sim.read_name import parse_origin

_SIM = ReadSimConfig(frag_mean=150, frag_std=20, frag_min=80, frag_max=250, read_length=50, seed=3)


def _origins(bam_path) -> Counter:
    with pysam.AlignmentFile(str(bam_path), "rb", check_sq=False) as bam:
        return Counter(parse_origin(r.query_name).kind for r in bam if r.is_read1)


def _nascent_by_span(result) -> Counter:
    """The nascent reads by the start of the span that emitted them (T1's span starts at 500, T2's at 3,500)."""
    start = dict(zip(result.index.t_df["t_id"], result.index.t_df["start"]))
    with pysam.AlignmentFile(str(result.bam_path), "rb", check_sq=False) as bam:
        origins = (parse_origin(r.query_name) for r in bam if r.is_read1)
        return Counter(int(start[o.transcript_id]) for o in origins if o.kind == "nrna")


def _scenario(
    tmp_path, *, abundance: float, nrna_abundance: float, other: float = 100.0
) -> Scenario:
    """A two-exon transcript on + at the abundances given, and one on − at ``other`` with no nascent abundance of
    its own, so under the scenario's own abundances every nascent read is the first transcript's."""
    sc = Scenario("draws", genome_length=6000, seed=5, work_dir=tmp_path / "sim")
    sc.add_gene(
        "g1",
        "+",
        [
            {
                "t_id": "T1",
                "exons": [(500, 1000), (2000, 2500)],
                "abundance": abundance,
                "nrna_abundance": nrna_abundance,
            }
        ],
    )
    sc.add_gene(
        "g2", "-", [{"t_id": "T2", "exons": [(3500, 4000), (4800, 5300)], "abundance": other}]
    )
    return sc


def test_a_transcript_with_no_mature_molecules_keeps_its_own_nascent_ones(tmp_path):
    """Abundance 0 with a nascent abundance of 50: about half the 400 RNA fragments are T1's nascent reads (50
    molecules over its 2,000 bp span against 100 over T2's 1,000 spliced bases), so a gate on the mature abundance,
    which drew none, cannot pass by chance. The control: the build's own override (``nrna_abundance=`` on the build)
    names the nascent abundance of every expressed multi-exon transcript, so T2 gets it and the silent T1 none."""
    sc = _scenario(tmp_path / "own", abundance=0.0, nrna_abundance=50.0)
    result = sc.build_oracle(n_rna_fragments=400, sim_config=_SIM)
    assert sum(_origins(result.bam_path).values()) == 400
    by_span = _nascent_by_span(result)
    assert set(by_span) == {500} and by_span[500] > 100, by_span

    sc = _scenario(tmp_path / "override", abundance=0.0, nrna_abundance=50.0)
    result = sc.build_oracle(n_rna_fragments=400, sim_config=_SIM, nrna_abundance=30.0)
    by_span = _nascent_by_span(result)
    assert set(by_span) == {3500}, by_span
