# THE REFERENCE PRIOR — a self-contained brief for external review

    ⚠ **A DEV DOC, written to be read by someone with no Rigel context.** It states the problem, the
    code as it is, the measurements, and the open question. It proposes nothing. Every number is
    reproducible with the command named beside it.

---

## 1. WHAT THE TOOL IS DOING, IN ONE PARAGRAPH

Rigel quantifies RNA-seq transcripts from libraries contaminated with genomic DNA. Before it can assign
RNA to transcripts it must **deconvolve each genomic location into gDNA versus RNA**. The genome is cut
into an alternating chain of **REGIONs** (intervals: intergenic, intronic, exonic) and **BOUNDARYs**
(the points between them). At each slot the accumulator supplies a count of fragments, and calibration
must split that count into gDNA and RNA. That per-slot split is `f_g`, the gDNA fraction, and solving
for it everywhere is "calibration".

There are exactly three populations at any slot — `gDNA`, `RNA+`, `RNA−` — and which are admissible is
fixed by the annotation. Where only one RNA strand is admissible the problem is one-dimensional and is
solved on `f_g` alone; where both are, a nuisance "tilt" parameter appears.

---

## 2. THE SOLVER

### 2.1 The coordinate

The solve runs on the **log-odds** axis

    lambda = log( f_g / (1 - f_g) )        f_g = sigmoid(lambda)

discretised to a uniform grid of `K` points on `[-L, +L]`, shipped `L = 10`, so `f_g` ranges over
`[4.54e-5, 0.99995]`. The posterior is evaluated at every grid point and the answer is read off as a
posterior summary. **The vertices `f_g = 0` and `f_g = 1` are not representable** — they sit at
`lambda = -inf` and `+inf`. This matters: several true compositions on the benchmark ARE exactly 0 or
exactly 1, so the lattice can only ever approach them.

### 2.2 What is added up at each grid point

`simplex_logodds.py:776` — the whole posterior, in log space:

    psi = <data likelihoods and messages>
        + _gdna_arm(lam, priors.gdna)          # +1/2 * log f_g      [+ fitted gDNA prior, if any]
        + _rna_arm (lam, priors.rna)           # +1/2 * log(1 - f_g) [+ fitted RNA prior, if any]
        + _location_term(lam, priors.location) # the LOCATION TILT — the subject of this document

The data terms are: a **strand likelihood** (a Beta-Binomial on the observed +/- read split, whose
information is `I ∝ N_eff · disc · [f_g(1-f_g)]^2 / (4p(1-p))`), an optional **density factor** (a
NegBinom comparing the slot's count against a fitted background rate), and optional **messages** from
neighbouring slots.

### 2.3 THE REFERENCE HAS TWO SEPARABLE PARTS, AND CONFLATING THEM IS THE ROOT OF THE CONFUSION

**Part A — the STABILIZER.** `+1/2*log f_g + 1/2*log(1-f_g)` is a **Jeffreys Beta(1/2, 1/2)** on the
composition. It is written **unconditionally**, by `_gdna_arm` and `_rna_arm`, and it is not
configurable. It is symmetric: it asserts nothing about which way the composition leans.

**Part B — the LOCATION TILT.** `_location_term` adds

    -log[ (1 - m) * f_g + m * (1 - f_g) ]

which moves the reference's mean from 1/2 to `m` while leaving Part A's tails untouched. Its range over
the grid is `|logit m|` nats. **It returns EXACTLY 0.0 when `m = 1/2`**, and it is fully switchable off
by two existing config flags (`structural_reference`, `measured_intron_reference`).

⭐ **Part A is the numerical stabilizer. Part B is the only prior ASSERTION in the system.**

### 2.4 WHY A STABILIZER IS GENUINELY REQUIRED (the Haldane question, settled)

`simplex_logodds.py:18-20` states it, and it is a property of the coordinate rather than a modelling
choice:

> *"Omitting a component's term is not 'no prior' — the grid's own measure supplies one. A bare
> uniform-`lambda` grid IS Haldane per component ⇒ Beta(0,0) on the composition: improper at BOTH
> vertices, a vertex amplifier. There is no third option. (This is what shipped before; it hid while
> psi stayed symmetric.)"*

A uniform grid in `lambda` is **not** a uniform prior on `f_g`; it is Beta(0,0), which has infinite mass
piled at both vertices. So *something* must bound the two ends. Jeffreys Beta(1/2, 1/2) is the weakest
standard choice that is proper, and it is what is written.

⛔ **So there are two different "weaken the reference" knobs and they have opposite consequences:**

| knob | what it does | is it proper? |
|---|---|---|
| remove **Part B** (the tilt, `m = 1/2`) | drops the only assertion; keeps Beta(1/2,1/2) | ✅ **proper, stable** |
| remove **Part A** (the Jeffreys halves) | leaves the bare grid measure = Beta(0,0) | ⛔ **Haldane, improper, vertex amplifier** |

**Removing the assertion entirely does NOT give Haldane.** The weakest proper configuration is
reachable today with two existing flags, and it is Jeffreys.

⚠ *Verified empirically*: a single-strand slot at `kappa = 1/2` with the tilt removed returns exactly
`f_g = 0.5000` — finite, centred, no pole.

---

## 3. HOW THE LOCATION IS SET TODAY

Two estimators compose, in `calibrate.py`'s sweep:

**(a) `structural_reference_location`** — from the ANNOTATION alone:

    m = (a+1)/(a+b+1) = 0.75   wherever NOT mrna_active   (a = b = 1/2, i.e. "one pseudo-observation
                                                            of gDNA on the Jeffreys exponents")
    m = 1/2                    wherever mrna_active        (⇒ the term is identically zero)

`mrna_active` means an annotated mature transcript runs continuously across the slot. So **intergenic
regions, intron regions and intron-flanking boundaries get `m = 0.75`; every exon gets NOTHING.**

**(b) `measured_reference_location`** — from the DATA, overriding (a) at single-stranded intron REGIONs:

    m_i = rho_bg * E_g,i / M_i         rho_bg = (sum of intergenic counts)/(sum of intergenic gDNA
                                                 effective lengths)

i.e. "the gDNA I expect here from the global background rate, over what I actually observe here",
clipped into the lattice window.

⚠ **(b) is data-derived, and that is the objection this document exists to record.** A prior is what is
known before the data; (b) is a statistic of the data being deconvolved.

---

## 4. THE CONFLICT, STATED PRECISELY

### 4.1 The nominal strength claim is honest as a pseudo-count and false as leverage

The term is documented as worth "one pseudo-fragment" (`a + b = 1`). Measured directly on the solver —
one slot, true `f_g = 0` (pure RNA), reference at `m = 0.75`, noiseless strand split, so the only
contest is reference versus likelihood:

| strandedness `kappa` | N=10 | N=100 | N=1,000 | N=10,000 | N=10^6 | overturned at |
|---|---|---|---|---|---|---|
| 0.99 | 0.120 | 0.030 | 0.009 | 0.003 | 0.0004 | **N = 3** |
| 0.90 | 0.184 | 0.039 | 0.012 | 0.004 | 0.0005 | **N = 3** |
| 0.75 | 0.441 | 0.073 | 0.020 | 0.006 | 0.0008 | **N = 10** |
| **0.50** | **0.747** | **0.747** | **0.747** | **0.747** | **0.747** | ⛔ **NEVER** |

**Where the strand channel carries information the term is worth about one fragment, exactly as
claimed.** Where it does not, the term is the entire answer at any depth.

**The mechanism is not mysterious.** Leverage is the reference's nats against the data's own Fisher
information *on the lambda axis*, and that information is

    I_strand  ∝  N_eff * disc * [ f_g (1 - f_g) ]^2 ,      disc ∝ (kappa - 1/2)^2

which is **identically zero at `kappa = 1/2`** (unstranded) and which **collapses as `f_g -> 0 or 1`**
even when stranded. The tilt's nats do not collapse. So the ratio "prior : data" diverges exactly where
the composition is extreme or the library is unstranded — which is most of the benchmark.

⭐ This is the same coordinate asymmetry that makes imputed messages over-confident elsewhere in the
system: one mechanism, several symptoms.

### 4.2 The score is largely a measurement of the reference, not of the algorithm

Pass-0 scores misplaced gDNA fragments over the slots it CLAIMS to solve. But a claimed slot with no own
evidence is not solved — its answer is whatever the tilt asserts. Partitioning claimed slots by whether
their posterior moves when the tilt is neutralised (a derived partition; no threshold needed to report
it):

| condition | policy | claimed live slots the reference materially moves | share of the SCORED error |
|---|---|---|---|
| `g50` stranded, capture-ON | silent | 6,633 / 15,261 | **25.8 %** |
| `g50` stranded, capture-ON | fanout | 9,451 / 15,261 | **40.3 %** |
| `g98` stranded, capture-ON | fanout | 11,117 / 17,216 | **28.7 %** |
| `g50` unstranded, capture-OFF | silent | 17,380 / 23,013 | **81.9 %** |
| `g50` unstranded, capture-ON | fanout | 14,540 / 15,311 | **95.2 %** |

⛔ **At the unstranded conditions, four fifths to nineteen twentieths of what we score is the tilt.** An
algorithmic change is then being graded against a number it barely controls, and a tilt that happens to
point the right way reads as a good algorithm.

⚠ **Caveat on this table**: neutralising the tilt also changes a fitted prior that is refit three times,
so "moved" mixes the tilt's direct effect with that indirect one.

⚠ VOCABULARY (owner, 2026-08-24): in the codebase "tilt" means the RNA+ vs RNA− nuisance
axis; this document's "location tilt" is a DIFFERENT object and is called the **location term**
below to avoid the collision.

**RESOLVED 2026-08-24 — the partition re-run at `calib_refit_iters = 0`, per population, exact-move
criterion, with the movement's MAGNITUDE beside it** (`Σ|Δest|` = total per-slot movement between the
shipped location term and `w = 0`, in fragments, against the shipped arm's scored error):

| condition | policy | pop | claimed | moved | err (loc) | err (w=0) | Σ\|Δ\| |
|---|---|---|---|---|---|---|---|
| `g50` ss.99 ON | silent | B | 17,518 | 11,733 | 38,650 | 43,533 | 6,824 |
| `g50` ss.99 ON | silent | E | 10,162 | **0** | 17,359 | 17,359 | **0** |
| `g50` ss.99 ON | fanout | B | 17,518 | 11,733 | 41,416 | 54,971 | 20,629 |
| `g50` ss.99 ON | fanout | E | 10,162 | 2,498 | 19,106 | 20,618 | 2,670 |
| `g98` ss.99 ON | fanout | B | 17,518 | 13,785 | 54,571 | 82,888 | 28,561 |
| `g98` ss.99 ON | fanout | E | 10,162 | 2,695 | 23,910 | 28,543 | 6,266 |
| `g50` ss.50 OFF | silent | B | 17,518 | 17,518 | 68,577 | 85,742 | 67,435 |
| `g50` ss.50 OFF | silent | E | 10,162 | **0** | 268,498 | 268,498 | **0** |
| `g50` ss.50 ON | fanout | B | 17,518 | 11,786 | 187,485 | 448,749 | 270,884 |
| `g50` ss.50 ON | fanout | E | 10,162 | 2,529 | 270,741 | 319,825 | 66,813 |

Three clean separations the shipped-refits table could not make:

1. **On the boundary population the location term's influence is DIRECT and total**: every live claimed
   boundary slot moves when it is neutralised (the unmoved remainder carries exactly zero error),
   and at the unstranded conditions the movement magnitude is the size of the entire scored error
   (67,435 vs 68,577) or larger (270,884 vs 187,485). At the stranded rows it is 18–52 % of the
   error in magnitude — present but no longer the whole answer.
2. **On exons the location term has ZERO direct effect** — under silence not one claimed exon moves, because
   `structural_reference_location` gives every `mrna_active` slot the neutral ½ and the measured
   reference is intron-only. Exons hear the location term ONLY through delivered messages (fanout rows:
   2,498–2,695 slots, 59–85 % of the exon error mass sits on them).
3. **The refit indirection the caveat worried about is real but secondary**: at zero refits the
   location term still helps every contaminated row (err loc < err w=0 throughout), so the earlier
   conclusion stands on the direct effect alone.

### 4.3 And the tilt does happen to point the right way, which is the trap

Certified truth for the mass-weighted gDNA fraction at intron-flanking boundaries — the population the
`m = 0.75` tilt speaks to — is not a constant:

| condition | certified `f_g` at intron boundaries |
|---|---|
| `g05` capture-OFF | **0.115** |
| `g05` capture-ON | 0.662 |
| `g50` capture-OFF | 0.712 |
| `g50` capture-ON | 0.974 |
| `g98` capture-OFF | 0.992 |
| `g98` capture-ON | 0.999 |

It spans from 0.115 to 0.999 — **the required tilt changes sign across the benchmark.** A sweep of the
constant confirms the optimum tracks that truth and no single value is best at more than a couple of
conditions; the shipped 0.75 is optimal at none of them.

**But removing the tilt makes the benchmark worse**, because on this panel it more often points the
right way than not. Sweeping the tilt's weight (1 = shipped, 0 = Jeffreys only), misplaced gDNA
fragments at claimed intron boundaries:

| condition | policy | w=1 | w=0.3 | w=0.1 | w=0 |
|---|---|---|---|---|---|
| `g00` unstranded OFF (zero control) | silent | 100,297 | 66,840 | 54,478 | **48,346** |
| `g05` stranded ON | silent | **22,714** | 25,289 | 26,454 | 27,131 |
| `g50` stranded OFF | fanout | **20,886** | 24,504 | 25,245 | 25,505 |
| `g50` unstranded ON | fanout | **251,503** | 530,581 | 656,091 | 712,437 |
| `g98` unstranded OFF | fanout | **20,020** | 34,544 | 38,057 | 39,574 |

⚠ **`w` is a scale on the location term's log-density**, not on the Jeffreys stabilizer, and not a
member of the derived family at intermediate values — `w = 1` and `w = 0` are exact, the middle is an
interpolation.

**So: weakening the assertion helps at the zero controls and hurts everywhere the library actually has
gDNA — up to 2.8x.** The tilt is load-bearing not because it is right but because it is the only voice
at slots where no channel speaks.

---

## 5. THE POSITION THIS DOCUMENT IS RECORDING

The tool's author's position, which the measurements above support:

1. **There is no prior here.** Pass-0's purpose is to *learn* the prior (a fitted gDNA-abundance
   landscape) for a later full solve. Manufacturing one in order to compute it is circular, and
   estimating it from the data being deconvolved (§3b) is circular twice over.
2. **What is genuinely known a priori is STRUCTURE only** — which slots are intergenic, intronic,
   exonic, and the corresponding boundaries — plus, arguably, one global scalar: the intergenic gDNA
   background rate, which imposes a floor everywhere but does not account for hybrid capture.
3. **Therefore the reference should be as weak as it can be while keeping the solver proper**, and
   everything else should come from data and from information propagated between slots.
4. **The benchmark should not grade slots the tool never solved.** Pass-0's deliverable is a confident
   SUBSTRATE for training the landscape prior; slots outside it should be excluded from its score, and
   the end-to-end judgement belongs after the full solve.

## 6. THE QUESTIONS FOR THE REVIEWER

1. **Is Jeffreys Beta(1/2,1/2) the right stabilizer**, given the solve is on a bounded log-odds lattice
   and several true compositions are exactly 0 or 1? Is there a weaker proper choice, or a
   reparameterisation that needs none?
2. **Is there a defensible prior at all** for "this interval is unannotated, therefore presumed gDNA
   until nascent-RNA evidence appears" — or is that a structural *constraint* that belongs in the
   likelihood rather than a prior?
3. **Given `I ∝ [f_g(1-f_g)]^2` vanishes at the vertices,** is the log-odds axis the right coordinate
   at all for compositions that are genuinely near 0 or 1? Would a different parameterisation make the
   prior's leverage uniform?
4. **Is the data-derived location (§3b) admissible** as an empirical-Bayes step, or is it the circularity
   the author suspects?
5. **What is the right scoring population** for a pass whose stated purpose is to produce a
   high-confidence subset rather than a complete answer?

---

## 7. HOW TO REPRODUCE EVERY NUMBER

* the leverage table (§4.1): a standalone solve of one slot through
  `rigel.calibration.simplex_logodds._solve_regions_logodds_all` with `CompositionPriors(location=...)`
* the weight sweep (§4.3) and the who-decides partition (§4.2): monkeypatch
  `simplex_logodds._location_term`, or set `CalibrationConfig(structural_reference=False,
  measured_intron_reference=False)` for the exact `w = 0` arm
* the certified truth (§4.3): `oracle_cache/<condition>/slot_truth.npz`, key `true_f_g`, produced and
  gated by `scripts/design/calibration_oracle.py`
* the claimed-slot scores: `scripts/design/pass0_claimed_ab.py [--dissect]`

⛔ One harness trap, paid three times: `rigel.calibration.__init__` rebinds the name `calibrate` to the
FUNCTION, so `import rigel.calibration.calibrate as m` yields a function and patching it is a silent
no-op; and `calibrate.py` binds its imports into its own namespace, so patching the defining module
misses too. Reach the module through `sys.modules` and assert the patch took.
