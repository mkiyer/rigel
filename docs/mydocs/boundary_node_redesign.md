# Linear Graph to Splice Graph

We need to complete a substantial architectural redesign starting from the Rigel index phase, cascading through to the EM phase. The change will touch many modules and have relatively high risk.

We partition the transcriptome into Regions and Boundaries. 
- Regions: contiguous genomic intervals
- Boundaries: single genomic point between regions

The partitioning creates a bipartite graph. Regions are connected by boundaries. Starting at any region, one must visit a boundary to reach another region.

The current bipartite graph is linear along the genome. This structure models unspliced fragments correctly. We stored a vector of regions and a vector of boundaries. At index 'i', boundary[i] corresponded to the left boundary of region[i], and boundary[i+1] corresponded to the right boundary of region[i]. This was a convenient and efficient indexing system.


We have a few realizations that create an impetus to redesign the boundary nodes:

1) boundaries should be modeled as independent first-class nodes that track fragments crossing a single point, and own the fragment mass on either side. therefore, a boundary has the effective length of the fragments crossing it.

2) spliced fragments are non-linear paths through the graph. our current design captures spliced fragments by splitting boundaries into two 'sides' (left/right) and tracking spliced fragments on only one side.





I think we have the opportunity to transition from the 2-sided boundary node design to a single boundary node design (1 sided, no split).

Boundary density calculations should be more robust if we combined the 2-sides into one. We should be able to do this without a larger refactor of the fractional accumulator code.

The effective length contraction code was able to do this correctly, combining the two sides of a boundary into one. You might be able to derive much of the implementation there.

The boundary object was originally developed with the concept of boundaries as flux measurements, not as independent entities. Now, it is more reasonable for boundaries to own their crossing fragments.

So how do we combine the two sides of a boundary into one?

We need access to both sides. A simple 




Also, it would be helpful for me to see plots of the node-pair model fit. It may be especially helpful to visualize the fit (and the decrease in error in fit) over iterations. We have a folder for intermediate results at "~/Downloads/rigel_runs" if you can share plots there and point me to them, it would be helpful.