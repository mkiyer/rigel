

# imputation variance modeling

node -> node imputation requires a variance model

currently we traverse pairs of adjacent nodes in the bipartite graph

for each node pair src, dst:
- extract src unspliced fractions { fpos_src, fneg_src, fgdna_src }
- extract dst unspliced fractions { fpos_dst, fneg_dst, fgdna_dst }

## gDNA model fit

- convert unspliced fractions to unspliced mass
- convert unspliced mass to unspliced density (mass/effective length)

- predictor: mean(density_gdna_src, density_gdna_dst)
- response: variance(density_gdna_src, density_gdna_dst)

Fit should be in log space

## RNA model fit

- Convert unspliced fractions to unspliced mass (fraction * count)
- mass_rnapos = fpos_src * (unspliced counts)
- mass_rnaneg = fneg_src * (unspliced counts)
- mass_gdna = fgdna * (unspliced counts)

If src is a BOUNDARY node, it may have spliced and unspliced mass. Add the spliced mass:
- mass_rnapos += mass_spliced_pos
- mass_rnaneg += mass_spliced_neg

Convert mass to density (divide by effective length)

RNA model fits EACH strand.

Positive strand:
- predictor: mean(density_rnapos_src, density_rnapos_dst)
- response: variance(density_rnapos_src, density_rnapos_dst)

Negative strand:
- predictor: mean(density_rnaneg_src, density_rnaneg_dst)
- response: variance(density_rnaneg_src, density_rnapos_dst)

Fit should be in log space


## Direct fit?

Is there a role for modeling `density_dst ~ density_src`?
Direct prediction of dst node density by src node density?
This is something I haven't considered


## Count variance consideration

Count variance is missing from our imputation fit

When we model using density, we BURY two important contributors to statistical power and variance:
1) the count itself
2) the effective length

Example:
- 10 counts / 2bp = 5 density
- 1000 counts / 200bp = 5 density

If our var ~ mean models densities, it fails to model the variance associated with random sampling of counts

My initial solution would be to seed the model with multiple poisson random samples from the observed counts.

Imputation variance should add the variance of the src and dst nodes

Nodes with low counts will have high variance and densities will fluctuate wildly

Nodes with high counts will have more stable density

So the variance measurement is comparing independent samples (poisson model is simplest sampling) with independent length (base pairs)

Modeling densities buries the discrete count nature of the data. 

We still want to build the imputation model

How can we do this?


draw poisson from 10 observed counts -> 10, 9, 11, 6, 10, 10, 12, 14, 15, etc..

We could do this for src and dst nodes, convert to densities, and add these as samples 


We can model var ~ mean

