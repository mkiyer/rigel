"""The per-node grid solve (`calibration.simplex_sweep`).

Step-1 precision rebuild: each region node is solved by its own local evidence ψ_i — the strand likelihood
(the only count→precision path) + the node-class prior (Jeffreys at single-strand, the global gDNA prior at
AMBIG / intergenic). There is NO cross-node propagation (the `q_rna` odds-coupling was removed; cross-node
imputation returns in Step 2) and NO count/RNA prior. These pin the three per-node behaviours:

  1. intergenic → the (0,0,1) gDNA vertex (the forbid mask), regardless of counts;
  2. a single-strand node is resolved by its strand (sense-tilted → RNA, balanced → gDNA);
  3. the global foundation breaks the AMBIG strand degeneracy (50/50 is equally pure-gDNA or
     balanced-two-strand-RNA) — ρ_global ≈ 0 ⇒ no phantom gDNA; high ρ_global ⇒ gDNA. This is what now
     resolves AMBIG nodes in production (the iterated calibrate() drives ρ_global down); the end-to-end
     AMBIG behaviour is locked by test_ambig_scenario.py.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from rigel.calibration.signature import TS_AMBIG, TS_NONE, TS_POS
from rigel.calibration.simplex_sweep import deconv_regions_sweep


def _sweep(ts, pos, neg, *, kappa=0.99, rho_global=0.0, region_eff_len=None, global_tau=None):
    pos = np.asarray(pos, dtype=np.float64)
    neg = np.asarray(neg, dtype=np.float64)
    n = pos.size
    contained = SimpleNamespace(
        n_unspliced_pos=pos, n_unspliced_neg=neg,
        mass_unspliced=pos + neg, mass_spliced=np.zeros(n),
        n_spliced_sense=np.zeros(n),
    )
    sub = SimpleNamespace(contained=contained)
    ra = SimpleNamespace(strand_class=np.asarray(ts))
    return deconv_regions_sweep(
        sub, ra, rna_sense_frac=kappa, gdna_strand_overdispersion=0.0,
        rna_strand_overdispersion=0.0, n_grid=40, rho_global=rho_global,
        region_eff_len=region_eff_len, global_tau=global_tau,
    )


def test_intergenic_locks_to_gdna_vertex():
    # intergenic (TS_NONE): allow_pos = allow_neg = False ⇒ the forbid mask leaves only the (0,0,1) vertex,
    # so f_g = 1 regardless of the per-strand counts (gDNA by signature).
    r = _sweep([TS_NONE], [50.0], [50.0])
    assert r.gdna_frac[0] > 0.99


def test_single_strand_sense_tilt_reads_rna():
    # a + region reading ~all-sense (99 / 1) at κ=0.99 is pure +RNA → f_g ≈ 0 (the strand resolves it).
    r = _sweep([TS_POS], [99.0], [1.0])
    assert r.gdna_frac[0] < 0.10


def test_single_strand_balanced_reads_gdna():
    # a + region reading 50/50 — half antisense is impossible for pure +RNA at κ=0.99, so the strand
    # likelihood favours gDNA. The Jeffreys reference concentrates f_g at the likelihood-favoured vertex.
    r = _sweep([TS_POS], [50.0], [50.0])
    assert r.gdna_frac[0] > 0.5


def test_global_foundation_breaks_ambig_degeneracy():
    # AMBIG 50/50: the strand likelihood is degenerate (pure-gDNA f_g=1 and balanced-two-strand-RNA
    # f_±=½ both predict 50/50). The global foundation breaks it: ρ_global ≈ 0 (pure-RNA library) pulls
    # f_g → 0 (no phantom gDNA); a high ρ_global pulls it toward gDNA. (Propagation's replacement.)
    eff = np.array([300.0])
    tau = np.array([20.0])
    lo = _sweep([TS_AMBIG], [50.0], [50.0], rho_global=0.0, region_eff_len=eff, global_tau=tau)
    hi = _sweep([TS_AMBIG], [50.0], [50.0], rho_global=1.0, region_eff_len=eff, global_tau=tau)
    assert lo.gdna_frac[0] < hi.gdna_frac[0]
    assert lo.gdna_frac[0] < 0.2  # the global says ~0 ⇒ no phantom gDNA on the degenerate ridge
    assert hi.gdna_frac[0] > 0.4  # a high baseline pulls the degenerate ridge toward gDNA
