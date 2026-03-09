We recently updated the code to decouple mature RNA from nascent RNA. Now the model is that nascent RNA "parents" can be associated with multiple mature RNA "children" depending on the splicing that happens between nascent RNA -> mature RNA. Our test scenario lack cases that examine and stress the new code. By developing a test scenario, we will likely improve the robustness of the code and even expose areas that can be improved. The test scenarios are being run using a simulation framework in the 'sim' module. We create a synthetic DNA genome (random nucleotides), then add the transcript structures and edit the genome to add the appropriate splicing motifs. This works wonderfully for testing. We need to add a scenario that stresses the new nascent RNA code. Here's the model: Start with a 50kb (50000) genome. Then lets create a first group of transcripts 'TA' on the positive strand: TA1 with exons [(1000, 2000), (5000,5500), (7000,7500), (9000,10000)], TA2 [(1000, 2000), (9000,10000)], TA3 [(1000, 2000), (5000,5500), (9000,10000)], TA4[(4500,5500), (9000,10000)]. Then we will create negative control transcripts, TB - strand exons [(20000, 24000)], and TC - exons [(30000,31000), (35000,36000)]. For now this should suffice. Now, we need to be able to set the abundances for each transcript (TA1, TA2, TA3, TA4), the nascent RNA for the 'TA' group, and the genomic DNA 'contamination' abundance. The negative control transcripts TB and TC will be ZERO (these will detect false positives). We need to run the scenario with many different abundance settings. Strand specificity can be set to 1.0 (perfect) to start. The first set of tests should be "one-hot" abundances for each TA transcript, a.k.a Test-TA1-1hot (TA1=100, TA2=0, TA3=0, TA4=3, nRNA-TA=0, gDNA=0), Test-TA2-1hot (TA1=0, TA2=100, TA3=0, TA4=3, nRNA-TA=0, gDNA=0) and so forth setting each transcript one at a time. Then, we should repeat the "one-hot" test but add a constant nascent RNA abundance e.g. Test-TA1-1hot-nrna100 (TA1=100, TA2=0, TA3=0, TA4=3, nRNA-TA=100, gDNA=0), Test-TA2-1hot-nrna100 (TA1=0, TA2=100, TA3=0, TA4=3, nRNA-TA=100, gDNA=0) and so forth. We need to see how accurately the tool can estimate the correct transcript and nascent RNA abundances. Then, for the third set of tests, we should repeat the 'one-hot' pattern with BOTH nascent RNA and genomic DNA. e.g. Test-TA1-1hot-nrna100-gdna100 (TA1=100, TA2=0, TA3=0, TA4=3, nRNA-TA=100, gDNA=100), Test-TA2-1hot-nrna100 (TA1=0, TA2=100, TA3=0, TA4=3, nRNA-TA=100, gDNA=100) and so forth. This will be a good start and give us some helpful information. It will initially confirm that the abundance estimation works under "pristine conditions" (nRNA=0, gdna=0, ss=1.0). Then examine how well we can deconvolve nascent RNA from mature RNA in an "easy" scenario (gDNA=0, ss=1.0) and gDNA-contaminated "easy" scenario (gDNA=100, ss=1.0). We should use the results to determine whether Rigel is working well. The results may help us to tune our priors. Once we are satisfied, we will expand the scenario to include more tests (sweep transcript abundance from 0 to 1000, sweep nRNA abundance, sweep gDNA abundance, sweep strand specificity). Lets start with a narrow scope. To facilitate this kind of testing, please create a new script in 'scripts' called 'synthetic_sim_sweep.py'. This script will look a lot like the test scenarios in the 'tests/scenarios' folder except that we should support a CLI. The script could accept a GTF file that defines the transcript structures and a reference size (reference genome length) representing the random synthetic genome size (single chromosome reference). That way the script will be extensible and can be used to stress the tool under a variety of conditions. Then the script should accept a configuration file YAML. The YAML should allow specification of transcript abundances, nascent RNA abundance, genomic DNA abundance, and strand specificity. For each entity/structure, we can specify a constant (ex. TB: 0), a list of values (TB: [0, 1, 10, 100]), or a 'sweep' of values in the form of (start, end, step). Then the YAML will be used to generate a combinatorial test of all values of all  transcript/nascentRNA/gDNA structures (can be a huge number of simulations). The results of each simulation run should capture estimated transcript counts, nascent RNA counts, gDNA counts. This can be compared versus truth. The total number of fragments (reads) that we generate for each test can be another parameter (constant, list, sweep). This is just the first pass idea for this, feel free to refine yourself.


# section where we assign labels to nascent RNAs
nrnas:
  NTA: (ref, +, 1000, 10000)

then in 'sweep' or 'patterns' we can reference this nrna

```
sweep:
  NTA: [0, 100]
```

or in patterns:

```
patterns:
- {TA1: 100, TA2: 0, TA3: 0, TA4: 3, NTA: 100}
```

sometimes 'sweep' and 'patterns' will conflict.. (can't both sweep and pattern of an entity..)

the YAML will fail validation if the nRNA 'key' doesn't match a transcript genomic span (ref, strand, start, end)




