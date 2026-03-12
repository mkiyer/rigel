# Bridge: Current System to Target Model

This note connects the current Rigel EM machinery to the redesign target.

Its purpose is to make the migration explicit rather than jumping directly from
the audit to the new theory.

## 1. Core Mapping

| Current mechanism | Target counterpart | Main gap |
|---|---|---|
| targeted excess penalty | Beta prior on gDNA strand split | current code is one-sided and approximate |
| `gdna_init` eligibility gate | explicit gDNA abundance prior or eligibility rule | binary gate versus continuous prior |
| `nrna_init` eligibility gate | explicit nRNA occupancy prior | binary gate versus continuous prior |
| intergenic density anchor | regional purity-weighted calibration | single anchor versus regional evidence model |
| intergenic FL training | purity-weighted gDNA FL calibration | sparse intergenic dependence |
| OVR warm start plus prior | simplex prior plus initialization | current dual role is conflated |

## 2. Current Targeted-Excess Penalty

The native solver currently protects a symmetric gDNA baseline and discounts
only asymmetric excess.

Conceptually it decomposes gDNA EM mass into:

1. a symmetric baseline `2 * min(sense, antisense)`
2. an asymmetric excess `abs(sense - antisense)`

The excess is discounted by a balance-dependent factor whose kernel matches a
symmetric Beta density.

This means the current implementation is already approximating the same prior
family proposed in the redesign. The gap is not biological intent. The gap is:

- the current implementation is implicit and one-sided
- the redesign target is explicit and parameterized in terms of the gDNA pair

## 3. Why This Bridge Matters

This bridge is important for two reasons.

1. It explains why the redesign is an evolution of the current system rather
   than a rejection of it.
2. It identifies the one-sided protection in the current penalty as a behavior
   that may need to be preserved explicitly because of the gDNA/nRNA
   identifiability problem.

## 4. Immediate Migration Consequences

Before the main EM refactor, the project should answer three concrete questions.

1. How closely does the current targeted-excess penalty match the intended Beta
   prior MAP update on synthetic data?
2. How much of the current gDNA/nRNA robustness comes from the one-sided nature
   of that penalty rather than from the prior family itself?
3. What asymmetric constraint on the `g_pos/g_neg` pair best preserves that
   robustness in the new solver?

This bridge analysis should be completed before the full solver migration is
treated as statistically frozen.