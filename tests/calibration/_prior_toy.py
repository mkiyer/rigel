"""The prior-vs-oracle toy — the one contaminated three-gene scenario the ``LocusPriors`` gates and the
message-layer gates (``_transfer_harness.capture_sweep_inputs``) both build on, so that one structure is defined
once. A ``_``-prefixed helper: pytest collects nothing here.
"""

from __future__ import annotations

from rigel.sim import GDNAConfig, ReadSimConfig, Scenario


def _scenario(name: str, seed: int, work_dir) -> Scenario:
    """Three genes, staggered isoforms, and one region shorter than the shortest fragment.

    The stagger is load-bearing for the ``boundary_spliced`` bank (a contiguous crossing by a molecule
    that spliced elsewhere can only land where a region_bound falls inside another transcript's exon), and the
    short region is what gives the toy genuinely EMPTY objects — the population the NaN-not-zero gate is
    about. Both are the same structures ``test_oracle_arms`` relies on, for the same reasons.
    """
    sc = Scenario(name, genome_length=9000, seed=seed, work_dir=work_dir)
    sc.add_gene(
        "g1",
        "+",
        [
            {"t_id": "t1", "exons": [(600, 1100), (1800, 2300)], "abundance": 60},
            {"t_id": "t1b", "exons": [(1000, 1100), (1800, 1900)], "abundance": 30},
        ],
    )
    sc.add_gene("g2", "-", [{"t_id": "t2", "exons": [(4000, 4500), (5200, 5700)], "abundance": 40}])
    sc.add_gene(
        "g3",
        "+",
        [
            {"t_id": "t3", "exons": [(6600, 7100), (7800, 8300)], "abundance": 25},
            {"t_id": "t3b", "exons": [(6600, 7100), (7400, 7440), (7800, 8300)], "abundance": 12},
        ],
    )
    return sc


_SIM = ReadSimConfig(
    frag_mean=180,
    frag_std=30,
    frag_min=80,
    frag_max=400,
    read_length=90,
    strand_specificity=0.99,
    seed=17,
)


def build_toy(tmp_path_factory):
    """Contaminated: gDNA at 1.0x the RNA fragment count, plus nascent RNA inside the introns."""
    sc = _scenario("pv", 17, tmp_path_factory.mktemp("pv_sim"))
    return sc.build_oracle(
        n_rna_fragments=4000,
        gdna_fraction=1.0,
        nrna_abundance=15.0,
        sim_config=_SIM,
        gdna_config=GDNAConfig(abundance=0.0, frag_mean=240, frag_std=45),
    )


def build_toy_zero_gdna(tmp_path_factory):
    """The zero-gDNA control, at the instrument level. Truth is exactly 0 at every locus, so there
    is nothing for a false positive to cancel against."""
    sc = _scenario("pv0", 23, tmp_path_factory.mktemp("pv0_sim"))
    return sc.build_oracle(
        n_rna_fragments=4000,
        gdna_fraction=0.0,
        nrna_abundance=15.0,
        sim_config=_SIM,
        gdna_config=GDNAConfig(abundance=0.0, frag_mean=240, frag_std=45),
    )
