
# Goals

1) Build a new, better set of realistic benchmarks for Rigel
2) Support the ability to rapidly stress test Rigel on specific scenarios/problems


# Problem

## Simulated != Real

Rigel has an excellent simulation system. It has the ability to generate a synthetic genome and synthetic transcripts. However, the synthetic transcripts do not mirror properties of the human transcriptome

## Need for rapid development

Rigel needs a large data substrate to train its models. If we need to run benchmarks with millions of reads, we will crush our development cycle.


# Solution

1) Simulate from the human genome and transcriptome

- Generate simulation data from the human genome and transcriptome. The data can be realistic (ex. 10M fragments). It can be oracle BAM data to start. Eventually, we would generate FASTQ and process it with an RNA-seq pipeline (alignment, etc).

2) Run Rigel's BAM scanner / accumulator and cache the result.

- 



Use real human genes as a backbone. This means that we can seed calibration with a backbone of real human genes (giving us RNA and DNA FL distributions, strand model, etc, all of the required items for calibration). 4a) Cache this intermediate (everything that we build from BAM scan -- the FL models, the strand model, the accumulated counts). 4b) Then piggyback on the cached objects, 
 'fake' chromosome with toy examples that stress test the system and allow us to test exact scenarios. 