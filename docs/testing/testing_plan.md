
# Goals

1) Build a new, better set of realistic benchmarks for Rigel
2) Support the ability to rapidly stress test Rigel on specific scenarios/problems


# Problem

Suppose we want to see how calibration performs on a small toy scenario. Currently we would need to generate realistic simulation data for 100s of transcripts in order to have reasonable fragment length distributions, strand model, and accumulator state. We need this to run our solver. This would cost a large amount of recomputation churn each time we want to test a small change.

Instead, we build the master scenarios ONCE and cache them. It gives us the substrate to rapidly iterate over small tests.


## Simulated != Real

Rigel has an excellent simulation system. It has the ability to generate a synthetic genome and synthetic transcripts. However, the synthetic transcripts do not mirror properties of the human transcriptome

## Need for rapid development

Rigel needs a large data substrate to train its models. If we need to run benchmarks with millions of reads, we will crush our development cycle.


# Solution

1) Simulate from the human genome and transcriptome

- Generate simulation data from the human genome and transcriptome. The data can be realistic (ex. 10M fragments). It can be oracle BAM data to start. Eventually, we would generate FASTQ and process it with an RNA-seq pipeline (alignment, etc).

2) Run Rigel's BAM scanner / accumulator and cache the result.

- The BAM scanner processes the fragments and deposits them into the accumulator

- We need to determine which objects calibration requires. Fragment length distributions. Strand model. Accumulator state.

- Cache everything (npz files / pickle / etc)


3) The cached results become the substrate for further testing

- Now, we can build a battery of small toy scenarios that can be tested against the background of the genome-scale cached state.

4) Limitations

- The toy scenarios must be run under the same configuration as the genome-scale state!

Constraints (fixed hyperparameters):
- strand specificity must be the same
- fragment length distributions must be the same
- overall gDNA levels must be the same (we model global gDNA from the genome-scale data)
- Hybrid capture modeling must be the same

Free parameters:
- Our toy genome (added as a synthetic chromosome)
- Transcript structures themselves (the GTF file)
- RNA abundances (including nascent) are free


# Objectives

1) Build the caching infrastructure

2) Develop a 'small' benchmark suite. Use a single human chromosome (e.g. chr22, relatively small) as the training substrate (should be big enough to train models, but small relative to full genome)

3) Variables

- total RNA fragments (5M)
- gDNA as a percentage of RNA: (zero, 5%, 50%, 100%, 200%)
- strand specificity (unstranded, 0.99)
- nascent RNA (none, present)
- hybrid capture (off, on, verystrong)

4) Run simulation and cache results

Now we have our SMALL benchmark suite based on one chromosome from human genome + transcriptome. As realistic as it gets.

5) Add "toy" scenarios on top

Next, we need the ability to simulate specific toy scenarios. We would use the simulator to generate a random mini-chromosome with transcripts. Run the oracle BAM simulator and Rigel's BAM scan pipeline.

- We then load the appropriate cached scenario and extract its models. We incorporate the toy scenario on top.

- Then run calibration. Voila!





