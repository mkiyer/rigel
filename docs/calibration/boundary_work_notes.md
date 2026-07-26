

# precision of enrichment ratios

if we compute an enrichment ratio over a very sparse denominator signal, it will be unstable.

for example if intron is depleted (near zero) then enrichment ratio is boundary density / (intron density ~ 0) --> goes to infinity.

so we shouldnt be computing enrichment ratios with imprecise density in the denominator. the enrichment ratio will be extremely imprecise.

that's a source of variance/precision for our composition/reframe method. to scale composition precisely, shouldnt the most precise observation be in the denominator?

i think this means, "the node with the highest total unspliced counts should be in the denominator".. that might not be correct, because length could be tiny. But can we derive the precision of the ratio?

for example:
- exon node (1000 counts / 1000 bp) = 1.0 density
- boundary node (50 counts / 100 bp) = 0.5 density
- intron node (1 count / 10000 bp) = 0.0001 density


we compute ratio of total density: 
- (boundary / intron) = 0.5 / 0.0001 = 5000

but intron is just one count that we are trusting... is it trustworthy? if the intron counts = 2, then ratio is 2500. a huge difference because of one count difference.

a single count scales the enrichment ratio by 2X

what if we reverse the ratio?

(intron / boundary) = 0.0001 / 0.5 = 0.0002

a single count still scales the ratio by 2X, but doesn't having the ratio in the numerator help us? 


# density is related to fragment length, and enrichment ratio (scaling) is dependent on density

lets work an example, where:
- DNA FL = 50
- RNA FL = 200

take a splice junction boundary with 100 spliced counts and 10 unspliced counts

for the spliced counts:
- spliced density = 100 / 200 = 0.5

for the unspliced counts:
- if 100% gdna, then density = 10 / 50 = 0.2
- if 100% rna, then density = 10 / 200 = 0.05

so the difference in fragment length can have a dramatic effect on total density.

boundary density = (0.5 spliced) + (0.05 - 0.2 unspliced)
- so boundary density upper and lower bound is known
- upper bound is 0.5 spliced + 0.2 = 0.7 (unspliced is gdna)
- lower bound is 0.5 spliced + 0.05 = 0.55 (unspliced is rna)


## how does this affect enrichment ratio?

enrichment ratio = boundary density / exon density

so just from fragment length alone, the enrichment ratio can vary. can we can compute the upper and lower bound if RNA and gDNA are free parameters.

*so enrichment ratio can vary by a known range*

but this neglects a key assumption:
** gDNA is (relatively, usually) uniform. uniform is a reasonable approximation **

so it should not be possible for boundary to be 100% dna and exon to be 100% rna... it would require some kind of rare genomic event (gdna deletion at exactly that location).. extremely unlikely. safer and more correct that gdna is uniform (has variance though, but relatively uniform)

# so can we derive this?

**the key missing piece is capture efficiency**
capture effiency is an unknown and makes a huge difference. capture efficiency is also not constant because it depends on probe position.

for this example, let's say capture effiency = 10X at locations with perfect probe overlap.

let's explore derivation.

given a sample with absolute gdna density 10 counts / 1000 bp = 0.01

then apply hybrid capture 10x:
- enriched gdna density = 10 x 0.01 = 0.1
- for simplicity we won't model depletion. just say depleted density remains 0.01

===========================

let's assume we have an exon with unknown density:
- 1kb exon with 2000 counts, unknown composition
- if 100% gdna, then exon density = 2000/(1000-50) = 2.10
- if 100% rna, then exon density = 2000/200 = 2.5

ENRICHED 1kb exon with 2000 counts

- gdna counts = 1kb x 0.1 counts/bp = 100
- rna counts = 1900
- density = (1900/(1000-200)) + (100/50) = 2.375 + 0.1 = 2.475


now the boundary:
- say the boundary is partially enriched, just 5X
- gdna density 0.05 x 50 = 2.5 counts (round to 2)
- we said we had 100 spliced counts and 10 unspliced counts

- rna density = (100 / 200) + (8 / 200) = 0.54
- dna density = (2 / 50) = 0.04
- total density = 0.58

so total density ratio is 2.475 / 0.58 = 4.26X

# enrichment ratio error is affected by fragment length difference

- above it is clear that difference in dna vs rna fragment length distribution affects enrichment ratio

enrichment ratio variance includes:
- rna count variance
- dna count variance
- rna FL distribution variability (RNA can have many different fragment lengths)
- dna FL distribution variability (DNA can have many fragment lengths)

All of this can vary for each node, and we take the ratio of two nodes.

# precision (belief) of a node affects the precision of the enrichment ratio

- an imprecise node can vary between RNA and gDNA, so the fragment length differential -> total density variation -> translates to larger variation of the enrichment ratio
- a PRECISE node has hardened its RNA vs DNA fraction -> fixed composition -> fixed enrichment ratio.


# open question: can we derive the "most likely" composition?

- does minimizing the "disagreement" in total density between two nodes lead to a "correct" solution?

- i don't think so.. but we may be onto a helpful computation or arithmetic that helps our solver.



=======



This term psi "φ = (the RNA at the seam, unspliced + spliced) / (the exon's actual RNA density), both in a common frame"

I think the literature has referred to this as percent spliced in or percent spliced out.

My intuition is that we might be able to try to derive this at least for Poisson counting without needing to model it in the data. The problem is the data will be over dispersed, and the Poisson model will be overconfident, but it probably adds what we're missing now. So let me talk you through what I'm envisioning.

We have a boundary with measured splice fragments, and so we have a defined density at the boundary. It's an RNA density. It's measured. It's fixed. It's invariant, and so we have that. And we can assume Poisson for the count process, but we really will have over dispersion in real data. It would be best to model with real data and not derive this from theoretical standpoint, but I think starting from theoretical underpinnings will help us to develop this better, and then we can figure out how to add over dispersion or change the model to a nonparametric model and of some sort.

So the derivation that needs to be done is assume we have a unknown RNA density. The unknown ground truth RNA density is being estimated. The exon itself is estimating the RNA density. If we assume we have ZERO gDNA, then we can model the variance of the imputation.

If we have zero DNA, then the boundary density is pure RNA regardless of whether it is spliced or unspliced, and we have a boundary density. That boundary density is imprecise because the boundary is short. In other words, it's a short, uh, fragment length that gives us our estimate. So whatever the counts are are measured over a very short length and adjusting the counts can cause significant changes to the density. Of course, this obscures away the problem of DNA, which, as you know, creates another source of variance because the density changes as the composition of DNA versus RNA changes because the fragment lengths change, and so that creates additional variation.

but in the most simple DNA free case, we just have RNA density at the boundary. And then in the exon, we have pure RNA density as well that's measured over a much longer exon length usually. Sometimes exons can be tiny and empty. In fact, if the length of an exon is too small, then the exon will have zero fragments. And so this all collapses in those cases. The exon length has to be large enough to accommodate many fragments. And so if the exon length is much larger than a fragment length, all of this is reasonably accurate. So we do have a case where the exon length becomes too small, in which case we can't really model the disagreement aware variance.

So we have two estimates of the same thing. One is coming from a short fragment length boundary with some fairly imprecise estimate. The other is coming from the exon itself, which may be more precise. Both of these are trying to estimate the same unknown ground truth RNA density.

So from a theoretical standpoint, we can understand this transfer variance concept. If we say we're measuring Poisson estimate of RNA at the boundary and Poisson estimate of RNA at the exon, we can derive what the error will be between the two estimates. and the error should be dependent on the region length because the exon region could be very short, in which case the error... exon will be very imprecise. and the axon could be very long, in which case it will be more precise. Either way, that's a baseline for what the disagreement will be and what the error will be based on Poisson counting principles.

So my question is, are we modeling at the very least this type of error, which is the disagreement between two Poisson processes estimating the same latent RNA density that is unknown. 

I think we can build up a correct count based error model this way. Then we need to add in over dispersion, which is difficult to measure again for the same reasons. We need some sort of observable measure of over dispersion, which we don't really have, and we have the circularity issues. And then we need to account for DNA, which is ALSO unknown. when we account for DNA, we have a similar error where the DNA is uniform between the boundary and the exon, and we have two estimations of the unknown DNA level, one at the exon, one at the boundary, and we can model the disagreement of those estimates the same way. So we can have a DNA transfer error imputation error, uh, and an RNA imputation error.

So this is my concept or mental model of what we're trying to model here. And I could be completely off base or off target, but I'd like you to take this, integrate it with what your thoughts are, and see if we can improve.