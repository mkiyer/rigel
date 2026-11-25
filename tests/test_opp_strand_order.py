"""Two transcripts on opposite strands, overlapping, swept across gDNA and nascent settings.

With neither gDNA nor nascent RNA the two transcripts are separate loci; nascent spans or gDNA candidates make them
share equivalence classes and so one locus. Each of the four settings is re-run with the fragments permuted, and the
answer must not depend on the buffer's order: an opposite-strand overlap is where a unit's candidate set stops being
obvious, so it is where an order dependence would first show.
"""

import numpy as np
import pysam
import pytest

from rigel.config import BamScanConfig, EMConfig, PipelineConfig
from rigel.pipeline import run_pipeline
from rigel.sim import GDNAConfig, ReadSimConfig, Scenario

N_FRAGMENTS = 2000
SIM_SEED = 42
PIPELINE_SEED = 42
GENOME_LENGTH = 20_000


def _sim_config(strand_specificity: float) -> ReadSimConfig:
    return ReadSimConfig(
        frag_mean=200,
        frag_std=30,
        frag_min=80,
        frag_max=450,
        read_length=100,
        strand_specificity=strand_specificity,
        seed=SIM_SEED,
    )


def _gdna_config(abundance: float) -> GDNAConfig | None:
    if abundance == 0:
        return None
    return GDNAConfig(abundance=abundance, frag_mean=350, frag_std=100, frag_min=100, frag_max=1000)


def _build_scenario(tmp_path, name: str) -> Scenario:
    """The 20 kb scenario: T1 on +, T2 on −, their spans overlapping and their exons not."""
    sc = Scenario(name, genome_length=GENOME_LENGTH, seed=SIM_SEED, work_dir=tmp_path / name)
    sc.add_gene(
        "G1",
        "+",
        [{"t_id": "T1", "exons": [(1000, 1500), (3000, 3500), (8000, 10000)], "abundance": 100}],
    )
    sc.add_gene(
        "G2",
        "-",
        [{"t_id": "T2", "exons": [(6000, 7000), (12000, 13000), (15000, 15500)], "abundance": 100}],
    )
    return sc


def _shuffle_bam(input_bam, output_bam, seed: int = 99) -> None:
    """Write *input_bam*'s read-name groups to *output_bam* in a permuted order, each group kept whole."""
    groups: dict[str, list] = {}
    with pysam.AlignmentFile(str(input_bam), "rb") as inf:
        header = inf.header.to_dict()
        for rec in inf:
            groups.setdefault(rec.query_name, []).append(rec)
    names = list(groups)
    np.random.default_rng(seed).shuffle(names)
    header["HD"]["SO"] = "queryname"
    with pysam.AlignmentFile(str(output_bam), "wb", header=header) as out:
        for name in names:
            for rec in groups[name]:
                out.write(rec)


def _answer(pr, result) -> dict:
    """The run's per-transcript counts, its gDNA and its loci table."""
    est = pr.estimator
    answer = {"gdna_em": float(est.gdna_em_count), "n_loci": len(est.locus_results)}
    for t in result.transcripts:
        answer[f"{t.t_id}_em"] = float(est.em_counts[t.t_index].sum())
        answer[f"{t.t_id}_unambig"] = float(est.unambig_counts[t.t_index].sum())
    for i, row in enumerate(est.get_loci_df(result.index).itertuples(index=False)):
        answer[f"locus_{i}"] = (
            row.n_transcripts,
            row.n_nascent,
            row.n_em_fragments,
            row.count_gdna_prior,
        )
    return answer


@pytest.mark.parametrize(
    ("gdna", "nascent", "strand_specificity", "n_loci"),
    [(0, 0, 1.0, 2), (0, 50, 0.9, 1), (50, 0, 0.9, 1), (50, 50, 0.9, 1)],
    ids=["rna_only", "nascent", "gdna", "nascent_and_gdna"],
)
def test_the_loci_and_the_answer_do_not_depend_on_fragment_order(
    tmp_path, gdna, nascent, strand_specificity, n_loci
):
    """The two transcripts are separate loci until nascent spans or gDNA candidates join them, and the run on the
    permuted BAM gives the same answer to the bit."""
    sc = _build_scenario(tmp_path, "opp_strand")
    result = sc.build_oracle(
        n_fragments=N_FRAGMENTS,
        sim_config=_sim_config(strand_specificity),
        gdna_config=_gdna_config(gdna),
        nascent_abundance=nascent,
    )
    config = PipelineConfig(
        em=EMConfig(seed=PIPELINE_SEED), scan=BamScanConfig(sj_strand_tag="auto")
    )
    original = _answer(run_pipeline(result.bam_path, result.index, config=config), result)
    assert original["n_loci"] == n_loci, original

    shuffled_bam = tmp_path / "shuffled.bam"
    _shuffle_bam(result.bam_path, shuffled_bam)
    shuffled = _answer(run_pipeline(shuffled_bam, result.index, config=config), result)
    assert shuffled == original
    sc.cleanup()
