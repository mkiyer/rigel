

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

so just from fragment length alone, the enrichment ratio:
- lower bound = (boundary 100% rna = 1.0) / (exon 100% gdna = 40) = 0.025
- upper bound = (boundary 100% dna = 2.5) / (exon 100% rna = 10) = 0.25

*so enrichment ratio can vary by factor of 10 maximum*

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
- gdna density 0.05 x 50 = 2.5 counts
- we said we had 1000 
rna density = 














# error is bounded by fragment length difference


# idea: minimize absolute density difference?

here is an concept probably worth exploring:

"the composition that minimizes the absolute difference between two nodes is most likely to be correct"

If the gDNA FL == RNA FL, then this does not matter.

When gDNA FL != RNA FL, then differences in the composition affect the total density of the node.

The total density