"""PROVE (or refute) the isoform-multiplicity EM siphon — pure EM, FAIR competition.

One locus. M = G+R ambiguous unspliced-exonic fragments, each equally likely (flat log-lik) under gDNA AND
every one of N identical single-exon transcripts (molecularly indistinguishable). We give BOTH sides a real
prior — gDNA prior = G, RNA prior = R — so at N=1 the correct split is G/R. Then we vary ONLY N and watch
the gDNA. If the split is N-invariant, the prior handles multiplicity correctly. If gDNA erodes as N grows
(same total RNA prior, just split among more transcripts + their nascent shadows), the isoform-multiplicity
siphon is a real structural EM bias: 1 gDNA component vs N RNA components.
"""
import sys
sys.path.insert(0, "/Users/mkiyer/proj/rigel/tests")
import numpy as np
from rigel.estimator import AbundanceEstimator
from rigel.config import EMConfig
from rigel.locus_partition import partition_and_free
from conftest import _make_locus_em_data


def run(N, G, R, rna_prior, em_iterations=200):
    M = G + R
    rc = AbundanceEstimator(N, em_config=EMConfig(seed=42))
    em_data, loci, gdna_prior_arr, index = _make_locus_em_data(
        [list(range(N))] * M, num_transcripts=N, include_gdna=True,
        gdna_prior_count=float(G), gdna_log_lik=0.0)
    partitions = partition_and_free(em_data, loci)
    parts = [partitions[i] for i in range(len(loci))]
    ptuples = [(p.offsets, p.t_indices, p.log_liks, p.coverage_weights, p.count_cols,
                p.is_spliced, p.gdna_log_liks, p.locus_t_indices, p.locus_count_cols) for p in parts]
    total_gdna, _, _ = rc.run_batch_locus_em_partitioned(
        ptuples, [loc.transcript_indices for loc in loci], gdna_prior_arr, index,
        rna_prior_count=np.array([float(rna_prior)]), em_iterations=em_iterations)
    return {"gdna": float(total_gdna), "mrna": float(rc.em_counts.sum()), "nrna": float(rc.nrna_em_count)}


for label, G, R, rna_prior in [("FAIR: gDNA prior=G=1000, RNA prior=R=1000", 1000, 1000, 1000),
                               ("gDNA minority: G=1000 R=3000, RNA prior=3000", 1000, 3000, 3000)]:
    print(f"\n=== {label} ===")
    print(f"{'N_tx':>5} {'gdna':>10} {'mrna':>10} {'nrna':>10} {'gdna_frac':>10} {'Δgdna_vs_N1':>12}")
    base = None
    for N in [1, 2, 5, 10, 30, 100, 300]:
        p = run(N, G, R, rna_prior)
        tot = p["gdna"] + p["mrna"] + p["nrna"]
        if base is None:
            base = p["gdna"]
        print(f"{N:>5} {p['gdna']:>10.1f} {p['mrna']:>10.1f} {p['nrna']:>10.1f} "
              f"{p['gdna']/tot:>10.4f} {p['gdna']-base:>+12.1f}")
