#!/usr/bin/env python
"""Quick test of the simulator zero-abundance fix."""
from rigel.sim.reads import SimConfig, ReadSimulator
from rigel.sim.genome import MutableGenome
from rigel.transcript import Transcript
from rigel.types import Strand, Interval

genome = MutableGenome(50000, seed=42)

cfg = SimConfig(frag_mean=250, frag_std=50, frag_min=80, frag_max=600,
                read_length=150, seed=42)

# Test 1: all-zero abundance
t_zero = Transcript(
    t_id="TA1", g_id="GA", ref="chr1", strand=Strand.POS,
    exons=[Interval(1000, 1020), Interval(5000, 5500), Interval(12000, 13000)],
    abundance=0.0, nrna_abundance=0.0,
)
sim = ReadSimulator(genome, [t_zero], config=cfg)
n_mrna, n_nrna = sim.compute_rna_split(10000)
print("Test 1 - all zero abundance:")
print("  compute_rna_split(10000) = ({}, {})".format(n_mrna, n_nrna))
assert n_mrna == 0, "FAIL: expected 0 mRNA"
assert n_nrna == 0, "FAIL: expected 0 nRNA"
print("  PASS")

n_m, n_n, n_g = sim._compute_pool_split(10000)
print("  _compute_pool_split(10000) = ({}, {}, {})".format(n_m, n_n, n_g))
assert n_m == 0 and n_n == 0 and n_g == 0, "FAIL: expected all zeros"
print("  PASS")

# Test 2: non-zero abundance
t_pos = Transcript(
    t_id="TA1", g_id="GA", ref="chr1", strand=Strand.POS,
    exons=[Interval(1000, 1020), Interval(5000, 5500), Interval(12000, 13000)],
    abundance=100.0, nrna_abundance=50.0,
)
sim2 = ReadSimulator(genome, [t_pos], config=cfg)
n_mrna2, n_nrna2 = sim2.compute_rna_split(10000)
print("\nTest 2 - non-zero abundance:")
print("  compute_rna_split(10000) = ({}, {})".format(n_mrna2, n_nrna2))
assert n_mrna2 + n_nrna2 == 10000, "FAIL: total should be 10000"
assert n_mrna2 > 0, "FAIL: expected some mRNA"
assert n_nrna2 > 0, "FAIL: expected some nRNA"
print("  PASS (total = {})".format(n_mrna2 + n_nrna2))

# Test 3: mRNA only (nrna_abundance=0)
t_mrna_only = Transcript(
    t_id="TA1", g_id="GA", ref="chr1", strand=Strand.POS,
    exons=[Interval(1000, 1020), Interval(5000, 5500), Interval(12000, 13000)],
    abundance=100.0, nrna_abundance=0.0,
)
sim3 = ReadSimulator(genome, [t_mrna_only], config=cfg)
n_mrna3, n_nrna3 = sim3.compute_rna_split(10000)
print("\nTest 3 - mRNA only:")
print("  compute_rna_split(10000) = ({}, {})".format(n_mrna3, n_nrna3))
assert n_mrna3 == 10000, "FAIL: expected all mRNA"
assert n_nrna3 == 0, "FAIL: expected 0 nRNA"
print("  PASS")

# Test 4: nRNA only (abundance=0)
t_nrna_only = Transcript(
    t_id="TA1", g_id="GA", ref="chr1", strand=Strand.POS,
    exons=[Interval(1000, 1020), Interval(5000, 5500), Interval(12000, 13000)],
    abundance=0.0, nrna_abundance=100.0,
)
sim4 = ReadSimulator(genome, [t_nrna_only], config=cfg)
n_mrna4, n_nrna4 = sim4.compute_rna_split(10000)
print("\nTest 4 - nRNA only:")
print("  compute_rna_split(10000) = ({}, {})".format(n_mrna4, n_nrna4))
assert n_mrna4 == 0, "FAIL: expected 0 mRNA"
assert n_nrna4 == 10000, "FAIL: expected all nRNA"
print("  PASS")

print("\nAll tests passed!")
