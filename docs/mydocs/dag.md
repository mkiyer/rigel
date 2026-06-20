
# Transition from Linear to Directed Acyclic Bipartite Graph

## Premise

Rigel's index step partitions the genome into a linear bipartite graph consisting of nodes.


## Node design


# Fitting variance mean model

We have source node and destination node.
During imputation, destination is the unknown, source is the known.

During training, we fit our model on both source and destination.

Currently, I believe we fit using the mean of source and destination. We have variance response variable given the mean of two nodes as the predictor variable.

Then during imputation, we query the variance-mean fit using the current value of the destination node? 

That seems counterintuitive. Why would we query the variance the current state of the node we are trying to impute? Doesn't the variance depend on both nodes and not just one? We fit the model using the mean of two nodes. Should we query the variance based on the mean of two nodes?

I want to get this part right. 



## Pseudocount stabilizer

We measure abundance as discrete counts. The true relative molecular abundance becomes highly uncertain when the number of counts approaches zero. When we have a count of zero, we have high uncertainty.

Imagin a 1kb genomic region that we sample from. If we randomly sample from a very low ground truth absolute molecular density, we will need extremely sparse counts. We could often see 0, sometimes 1, rarely 2, etc. 

When counts are very low, we cannot precisely estimate the true molecular abundance. Our estimates are notoriously imprecise.

Our variance model needs to reflect this. Currently, we fit a model that returns ZERO variance when we give it a density of 0. This is mathematically correct but biologically egregious because it doesn't respect the imprecision associated with sampling discretely from an extremely rich molecular pool.

Our variance ~ mean model should be imprecise at low counts.

We did incorporate poisson sampling variance into the model fit.

However, I believe we still fit with count = 0 and density = 0.

A few ideas:
- Add a pseudocount of 1 to all observations to reflect the molecular uncertainty that is associated with low counts.
- density used for fitting is computed from (count + 1)/effective_length instead of (count)/effective_length
- OR, we add pseudocount to poisson so that we never have ZERO in the model.

Other observations:
- Typically we fit these models in log space. They tend to behave much better that way. We deviated from log-space fitting because of the complexity introduced when we added the poisson variance to the biological variance.
