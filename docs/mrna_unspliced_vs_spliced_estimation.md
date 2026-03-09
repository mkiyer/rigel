This is an incredibly powerful idea. You are describing what we can call a **Geometric Splicing Expectation** model.

By utilizing the known exon/intron geometry of the transcript and the fragment length distribution, you can deterministically predict the ratio of spliced to unspliced reads for mature mRNA. If you know that ratio, and you observe the spliced reads (which are *exclusively* mRNA), you can perfectly back-calculate the "phantom" unspliced mRNA reads and remove them from the pool before estimating gDNA and nRNA.

This elegantly solves the exact pathology you uncovered in your diagnostic sweep (where highly expressed mRNA was inflating the unspliced pool and causing the gDNA EB shrinkage to hallucinate).

Here is how we can develop the math and implement this efficiently.

### 1. The Math: Geometric Splicing Expectation

You mentioned sliding a fragment across the transcript. In RNA-seq, this is the classic definition of **Effective Length** ($\tilde{L}$).

Let $\mu_F$ be the mean fragment length.
For a mature transcript $t$ composed of exons with lengths $E_1, E_2, \dots, E_k$:

**Total Effective Length:**
The total number of valid starting positions for a fragment on the mature transcript is:


$$\tilde{L}_{total} = \left( \sum_{i=1}^k E_i \right) - \mu_F + 1$$

**Unspliced Effective Length (Exon-Internal):**
To be "unspliced," the fragment must fall entirely within a single exon. The number of valid starting positions inside exon $i$ is $E_i - \mu_F + 1$. Therefore:


$$\tilde{L}_{unspliced} = \sum_{i=1}^k \max(0, E_i - \mu_F + 1)$$

**Spliced Effective Length (Junction-Spanning):**
By definition, any valid fragment that is not entirely within one exon must cross at least one splice junction:


$$\tilde{L}_{spliced} = \tilde{L}_{total} - \tilde{L}_{unspliced}$$

**The Magic Ratio:**
Assuming uniform coverage (which is the standard, safe assumption for Bayesian priors), the expected ratio of unspliced to spliced reads for transcript $t$ is an immutable geometric property:


$$R_t = \frac{\tilde{L}_{unspliced}}{\tilde{L}_{spliced}}$$

### 2. The Strategy: 3-Step Deconvolution

With $R_t$ pre-computed, we can completely decouple the unspliced pool.

#### Step A: Estimate the mRNA "Unspliced Burden"

During initialization, count the number of **unambiguously spliced** reads assigned to transcript $t$ (let's call this $C_{spliced, t}$).
Because these reads crossed a junction, we are 100% certain they belong to mature mRNA.
We can now calculate exactly how many unspliced reads this transcript *must* have left in the pool:


$$\text{Expected\_mRNA\_Unspliced}_t = C_{spliced, t} \times R_t$$

Sum this across all transcripts in the locus to get the total mRNA burden:


$$B_{mRNA} = \sum_{t} \left( C_{spliced, t} \times R_t \right)$$

#### Step B: Isolate the True nRNA/gDNA Pool

Observe the total number of sense-strand unspliced reads in the locus ($U_{sense}$).
Subtract the mRNA burden to find the reads that *actually* belong to nascent RNA and genomic DNA:


$$U_{remainder} = \max(0, U_{sense} - B_{mRNA})$$


*(Note: We use `max(0, ...)` because real coverage isn't perfectly uniform, and we don't want negative counts).*

#### Step C: Strand-Specific Deconvolution

Now we apply your strand-specificity math, but strictly to the isolated remainder!
Let $A$ be the number of antisense unspliced reads.
Assuming a highly stranded library ($SS \approx 1.0$), all antisense reads are gDNA, and gDNA is symmetric.

* **gDNA Initial Estimate:** $2 \times A$
* **nRNA Initial Estimate:** $\max(0, U_{remainder} - A)$

### 3. Why This is a Massive Upgrade

1. **It solves H1 perfectly:** In your diagnostic sweep, when you introduced 10,000 mRNA reads, the old model applied the $2A$ rate to *all* of them, inflating gDNA. This new model will see the spliced reads, calculate $B_{mRNA} \approx 10,000$, subtract it from $U_{sense}$, and leave $U_{remainder} \approx 0$. gDNA initialization will correctly remain $0$.
2. **It dynamically scales with fragment length:** A dataset with 400bp fragments will span more junctions than a dataset with 150bp fragments. Because $\mu_F$ is in the $\tilde{L}$ equations, $R_t$ will automatically adjust, keeping the priors perfectly calibrated to the specific sequencing run.
3. **It operates entirely outside the EM:** All of this math can be executed during the `locus.py` initialization phase to generate the `alpha / beta` priors. It adds zero overhead to the C++ EM solver loop.

### 4. Implementation Thoughts

**Where to compute $\tilde{L}_{unspliced}$?**
You can actually compute this vector entirely in Python before the EM starts.
When building the `AbundanceEstimator`, you already have the `transcript_exonic_lengths` and the fragment length distribution ($\mu_F$).
You just need to pass the list of individual exon lengths per transcript from the `Index` to compute the $\max(0, E_i - \mu_F + 1)$ sum.

How does this logic feel to you? If you agree with the geometric derivation, we can draft the exact numpy operations to fold this into the `estimator.py` initialization pipeline.




I am considering another improvement related to gdna initialization.  

By utilizing the known exon/intron geometry of the transcript and the fragment length distribution, we can deterministically predict the ratio of spliced to unspliced reads for mature mRNA. If you know that ratio, and you observe the spliced reads (which are *exclusively* mRNA), you can perfectly back-calculate the "phantom" unspliced mRNA reads and remove them from the pool before estimating gDNA and nRNA.

This elegantly solves the exact pathology you uncovered in your diagnostic sweep (where highly expressed mRNA was inflating the unspliced pool and causing the gDNA EB shrinkage to hallucinate).

Here is how we can develop the math and implement this efficiently.

### 1. The Math: Geometric Splicing Expectation

You mentioned sliding a fragment across the transcript. In RNA-seq, this is the classic definition of **Effective Length** ($\tilde{L}$).

Let $\mu_F$ be the mean fragment length.
For a mature transcript $t$ composed of exons with lengths $E_1, E_2, \dots, E_k$:

**Total Effective Length:**
The total number of valid starting positions for a fragment on the mature transcript is:


$$\tilde{L}_{total} = \left( \sum_{i=1}^k E_i \right) - \mu_F + 1$$

**Unspliced Effective Length (Exon-Internal):**
To be "unspliced," the fragment must fall entirely within a single exon. The number of valid starting positions inside exon $i$ is $E_i - \mu_F + 1$. Therefore:


$$\tilde{L}_{unspliced} = \sum_{i=1}^k \max(0, E_i - \mu_F + 1)$$

**Spliced Effective Length (Junction-Spanning):**
By definition, any valid fragment that is not entirely within one exon must cross at least one splice junction:


$$\tilde{L}_{spliced} = \tilde{L}_{total} - \tilde{L}_{unspliced}$$

**The Magic Ratio:**
Assuming uniform coverage (which is the standard, safe assumption for Bayesian priors), the expected ratio of unspliced to spliced reads for transcript $t$ is an immutable geometric property:


$$R_t = \frac{\tilde{L}_{unspliced}}{\tilde{L}_{spliced}}$$

### 2. The Strategy: 3-Step Deconvolution

With $R_t$ pre-computed, we can completely decouple the unspliced pool.

#### Step A: Estimate the mRNA "Unspliced Burden"

During initialization, count the number of **unambiguously spliced** reads assigned to transcript $t$ (let's call this $C_{spliced, t}$).
Because these reads crossed a junction, we are 100% certain they belong to mature mRNA.
We can now calculate exactly how many unspliced reads this transcript *must* have left in the pool:


$$\text{Expected\_mRNA\_Unspliced}_t = C_{spliced, t} \times R_t$$

Sum this across all transcripts in the locus to get the total mRNA burden:


$$B_{mRNA} = \sum_{t} \left( C_{spliced, t} \times R_t \right)$$

#### Step B: Isolate the True nRNA/gDNA Pool

Observe the total number of sense-strand unspliced reads in the locus ($U_{sense}$).
Subtract the mRNA burden to find the reads that *actually* belong to nascent RNA and genomic DNA:


$$U_{remainder} = \max(0, U_{sense} - B_{mRNA})$$


*(Note: We use `max(0, ...)` because real coverage isn't perfectly uniform, and we don't want negative counts).*

#### Step C: Strand-Specific Deconvolution

Now we apply your strand-specificity math, but strictly to the isolated remainder!
Let $A$ be the number of antisense unspliced reads.
Assuming a highly stranded library ($SS \approx 1.0$), all antisense reads are gDNA, and gDNA is symmetric.

* **gDNA Initial Estimate:** $2 \times A$
* **nRNA Initial Estimate:** $\max(0, U_{remainder} - A)$

### 3. Why This is a Massive Upgrade

In your diagnostic sweep, when you introduced 10,000 mRNA reads, the old model applied the $2A$ rate to *all* of them, inflating gDNA. This new model will see the spliced reads, calculate $B_{mRNA} \approx 10,000$, subtract it from $U_{sense}$, and leave $U_{remainder} \approx 0$. gDNA initialization will correctly remain $0$.
2. **It dynamically scales with fragment length:** A dataset with 400bp fragments will span more junctions than a dataset with 150bp fragments. Because $\mu_F$ is in the $\tilde{L}$ equations, $R_t$ will automatically adjust, keeping the priors perfectly calibrated to the specific sequencing run.
3. **It operates entirely outside the EM:** All of this math can be executed during the `locus.py` initialization phase to generate the `alpha / beta` priors. It adds zero overhead to the C++ EM solver loop.

### 4. Implementation Thoughts

**Where to compute $\tilde{L}_{unspliced}$?**
You can actually compute this vector entirely in Python before the EM starts.
When building the `AbundanceEstimator`, you already have the `transcript_exonic_lengths` and the fragment length distribution ($\mu_F$).
You just need to pass the list of individual exon lengths per transcript from the `Index` to compute the $\max(0, E_i - \mu_F + 1)$ sum.

How does this logic feel to you? If you agree with the geometric derivation, we can draft the exact numpy operations to fold this into the `estimator.py` initialization pipeline.. We should be able to estimate the number of estimate mature RNA fragments based on the fragment length distribution and the transcript structures. Geometric Splicing Prediction: The DerivationYour idea to predict $P(\text{unspliced} \mid \text{transcript}, \text{frag\_len})$ is completely viable. In RNA-seq, sliding a fragment across a transcript is formalized as computing its Effective Length ($\tilde{L}$).Let $\mu_F$ be the mean fragment length. For a mature transcript $t$ composed of exons with lengths $E_1, E_2, \dots, E_k$:1. Total Valid Starting Positions (Total Effective Length):$$\tilde{L}_{total} = \left( \sum_{i=1}^k E_i \right) - \mu_F + 1$$2. Unspliced Starting Positions:For a read to be strictly "unspliced," it must fall entirely within the boundaries of a single exon.$$\tilde{L}_{unspliced} = \sum_{i=1}^k \max(0, E_i - \mu_F + 1)$$3. The Splicing Ratio:Assuming uniform transcript coverage, the probability that a fragment drawn from mature mRNA $t$ appears unspliced is an immutable geometric property:$$P(\text{unspliced} \mid t) = \frac{\tilde{L}_{unspliced}}{\tilde{L}_{total}}$$4. Executing the Deconvolution:During initialization, count the number of unambiguously spliced reads assigned to transcript $t$ ($C_{spliced}$). Because they cross a junction, we know with 100% certainty they are mature mRNA.We can calculate the "Phantom Unspliced Burden" that this mRNA must have left in the locus:$$Burden_{mRNA} = \sum_{t} \left( C_{spliced, t} \times \frac{P(\text{unspliced} \mid t)}{1 - P(\text{unspliced} \mid t)} \right)$$Subtract $Burden_{mRNA}$ from the total sense-strand unspliced pool before you ever calculate the gDNA or nRNA priors!  -- Can you do the following 1) Review and Explain the current gDNA initialization process. 2) Explain if and how gDNA can be overestimated when mRNA levels are very high. 3) Validate the accuracy of this theory, 4) After the BAM scanning pass, we could compute the P(unspliced | spliced, fragment length distribution) for each transcript. Ideally, we would initialize genomic DNA more conservatively based on 