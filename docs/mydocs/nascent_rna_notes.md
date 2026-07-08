
# thoughts on 3-component vs 2-component model

three component model: { mature RNA (+, -), nascent RNA (+, -), gDNA }
two component model: { RNA, gDNA }

it is quite possible that a 3-component model is over-engineered, and a simpler 2-component model is actually correct.


## problem that motivated 3 component model

this example shows the problem.
one transcript

here is the scenario:

TA+ exons (1000, 5000), (10000, 14000)
-region (1000,5000) unspliced density = 10 (rna vs gdna)
-boundary 5000 density spliced=1, unspliced=1
-region (5000,10000) intron unspliced density (variable, change to unveil pathology)
-boundary 10000 density spliced=1, unspliced=1 (rna vs gdna)
-region (10000, 14000) unspliced density = 10 (assume RNA=9, gDNA=1)
-intergenic gdna density (variable, sweep to unveil pathology)

the problem:
- often in hybrid capture, density is affected by placement of capture probes. this can lead to a discrepancy between region density and boundary density. total density region = 10, boundary total density = 2
- boundary spliced density acts as an RNA SINK, absorbing RNA, BUT leftover RNA could also BLEED into the intron.
- this bleeding behavior created phantom gDNA and phantom nascent RNA
- indirectly, the phantom gDNA + nascent RNA prevented intron nodes from collapse (effective length contraction). so the EM readily REASSIGNED fragments to mature RNA. Counterintuitive! so fragments were actually correctly assigned by the EM, BECAUSE of the phantom bleed gDNA + nascent RNA.
- cleaning up the phantom gdna + nascent RNA problem led to intronic node effective length collapse, making gDNA and nascent RNA more competitive, and then causing mature RNA to be outcompeted (!?!?!) this is the pathology.

