**Bottom Line**

Your concern is valid. The current calibration path is directionally sensible, but I do not think the locoregional gDNA prior is mathematically clean yet. The biggest issues are not just `mean_FL` versus `mean_FL - 1`; they are how local counts are paired with local exposure, how boundary-crossing fragments are reused, and how the estimated physical gDNA mass is compressed into a fixed-strength EM prior.

**What The Code Does Now**

The scanner records one `obs_mask` per unique, non-chimeric molecule, increments per-region counts, and separately increments integer `u_left/u_right` boundary counters for unspliced fragments crossing exon boundaries in accumulator.cpp.

Global densities are:

- `INTERGENIC`: intergenic-only count divided by contained FL-weighted effective length.
- `INTRON`: intron-only count divided by contained FL-weighted effective length.
- `EXON-INTRON`: boundary event count divided by `eligible_boundary_sides * gdna_fl.mean` in density_global.py.

For each locus, the code clips region lengths to the transcript-defined locus, shrinks local densities toward global densities, predicts `n_gdna`, converts that to `pi_gdna = n_gdna / n_obs`, then sets `alpha_gdna = c_base * pi_gdna` and `alpha_rna = c_base * (1 - pi_gdna)` in locus_prior.py and locus_prior.py. Native EM then appends a gDNA component for unspliced units and uses `alpha_gdna` as the gDNA Dirichlet prior in em_solver.cpp and em_solver.cpp.

**Main Findings**

1. The locoregional numerator and denominator can be mismatched.

In `estimate_locus_gdna()`, region lengths are clipped to the locus before computing `leff`, but `_shrink_one_type()` sums the entire region’s count with `count_col[rids].sum()` in locus_prior.py. If a locus overlaps only part of a large intron or exon-derived region, the numerator is whole-region while the denominator is clipped-region. That can inflate local density.

The same pattern exists for boundary flux: `u_left/u_right` and boundary flags are taken for the whole exon region, while `leff_full` is clipped to the locus in locus_prior.py. This is probably the highest-priority modeling issue.

2. The boundary denominator should be boundary exposure, not `mean_FL`.

For one boundary and fixed fragment length `l`, the number of start positions that cross the boundary is `l - 1`, not `l`. So the count-based MLE should be approximately:

```text
rho_boundary = N_crossing / (n_boundaries * E[L - 1])
```

Current code uses `n_boundaries * E[L]`, so it is conceptually off and slightly density-deflating. For typical `L ~ 300`, this is a small numerical error, but it is the wrong primitive.

Your proposed per-fragment quantity `1 / (l_f - 1)` is theoretically defensible as a Horvitz-Thompson style estimator:

```text
rho_boundary_HT = sum_i 1 / (l_i - 1) / n_boundaries
```

That is valid if you store the boundary flux as a float and do not divide by `E[L - 1]` again. I would either keep integer counts and use `E[L - 1]`, or store weighted float flux and use exposure `n_boundaries`. Mixing the two would double-normalize.

3. Current EXON-INTRON projection estimates exon-contained gDNA, not boundary-crossing gDNA.

The code estimates boundary density, then projects it onto contained exonic effective length via `leff_ei` and `n_gdna_ei = rho_loco_ei * leff_ei` in locus_prior.py. That imputes EXON_ONLY gDNA, but it does not explicitly add the observed EXON-INTRON fragments themselves.

Those boundary fragments still enter EM as unspliced units with a gDNA candidate, so they can be assigned to gDNA by likelihood. But the prior fraction `pi_gdna` uses `n_obs` over EM units, and boundary units are part of that denominator. If the numerator only contains imputed exon-contained gDNA, the prior is biased downward for loci where boundary fragments are a large share of the ambiguous units.

A cleaner decomposition would be:

```text
n_gdna_prior =
    n_gdna_intron_contained
  + n_gdna_intergenic_local_or_flank
  + n_gdna_boundary_observed_or_expected
  + rho_boundary_loco * exon_contained_exposure
```

Use boundary flux twice, but in two different roles: observed boundary gDNA mass, plus density anchor for EXON_ONLY imputation.

4. EXON_ONLY raw counts should not directly estimate gDNA density.

EXON_ONLY unspliced fragments are a mixture of mature RNA fragments contained within one exon, nascent RNA, and gDNA. They should stay latent in the EM unless you do a second-stage or iterative calibration using EM posterior expected gDNA counts. Raw EXON_ONLY counts are not a valid gDNA density estimator.

Also, EXON-INTRON is mature-mRNA-free, but it is not necessarily RNA-free. Nascent RNA can produce exon-intron unspliced fragments. So boundary flux is a strong gDNA anchor, but in nRNA-rich loci it can still overstate gDNA unless the model lets posterior nRNA evidence push back.

5. Intergenic evidence is mostly global, not locoregional.

The code has an intergenic branch, and global intergenic density is computed in density_global.py. But transcript-defined `Locus` intervals normally do not include neighboring intergenic regions, so the locoregional intergenic branch in locus_prior.py is usually zero. Deterministic intergenic fragments are counted as gDNA in stats/output, but they do not enter the EM.

Adding neighboring intergenic windows would likely improve local density estimates, especially for capture or regional contamination effects, but it needs a clear windowing/mappability policy.

6. The EM prior loses local evidence strength.

After estimating physical `n_gdna`, the code keeps only `pi_gdna` and scales it by a fixed `c_base` in locus_prior.py. A locus with 2 boundary events and a locus with 200 boundary events can get the same prior strength if their `pi_gdna` is the same. Kappa affects the density estimate, but not the Dirichlet evidence strength. That may be too blunt.

**Recommended Direction**

I would refactor the model around separate exposure primitives:

```text
B_eff = sum_boundary sum_l h(l) * max(l - 1, 0)
C_exon = sum_exon sum_l h(l) * max(exon_len - l + 1, 0)
```

Then estimate:

```text
rho_boundary = N_boundary_events / B_eff
n_exon_only_gdna = rho_boundary * C_exon
```

For the EM prior, include observed/expected boundary gDNA separately from imputed EXON_ONLY gDNA. Also fix locoregional clipping by either prorating counts to clipped exposure, expanding loci to whole evidence windows, or collecting finer positional bins so local numerators and denominators describe the same territory.

So, yes: boundary flux is the right signal, but the current implementation is best viewed as a first-order approximation. The next correctness step is to make the local numerator/exposure geometry explicit and separate “observed boundary gDNA” from “exon-only gDNA imputed from boundary density.”