# Shared capture field: follow-up review and executable design (2026-10-11)

Sandbox, not authoritative. Baseline `624dfade4afce37d795890ea8ab537086ea20a93`.
The owner asked this reviewer to finish the scientific and implementation work, not to request another
reviewer's permission. No commits. The submitted design is preserved in `SHARED_CAPTURE_FIELD_DESIGN.md`.
Measurements and isolated prototypes for this follow-up go in `<scratch>/rigel-shared-field-20261011`.

## Before implementation: the contracts to test

The promising idea is to infer one capture field jointly with allocation. Four assertions in the submitted
design need correction before that idea can become an implementation.

1. The conserved deposit is a **vector**, not necessarily one object. A contained fragment deposits on its
   region. A crossing fragment's slices distribute their bases among the cuts bounding them. The sum of those
   shares is one, but several shares can be nonzero. This is the executable accumulator specification and the
   definition used by `conserved_cut_shares`. Moreover, implicit splice hypotheses and multiple alignments can
   give different paths to different candidate components. The scorer currently collapses those alignments.
2. The EM bound with a responsibility independent of the field does not have the integrated responsibility
   as its maximizing update. Mixing that bound's field update with an integrated allocation update has no
   stated ascent guarantee. This is a mathematical question, to test on exact enumerations first.
3. The closed-form covariance in the design is `Cov(L_k,L_m)`, not `Cov(log L_k,log L_m)`. Applying the delta
   method to the former adds another approximation. Even with exact log-length moments a second-order softmax
   expansion need not be positive or accurate over a broad posterior.
4. Shared exons do not in general identify abundances separately from a freely varying field. The toy's
   captured-exon success may depend on its chosen prior and finite grid. An equal-likelihood ridge is a
   sharper test than another replicate average under that same prior.

## A kernel that conserves the actual shares

Let `d_fkh,o` be the conserved deposit vector for fragment `f`, candidate `k` and alignment/path `h`,
with `sum_o d_fkh,o = 1`. Let `phi_fkh` include the fragment and alignment likelihood on that path. Define

    C_fkh(c) = sum_o d_fkh,o c_o
    K_fk(c) = sum_h phi_fkh C_fkh(c)
    A_ko = sum_(admissible f,h) phi_fkh d_fkh,o
    L_k(c) = sum_o A_ko c_o.

Then `sum_f K_fk(c) = L_k(c)` by interchanging two finite sums. A fragment contributes exactly one mixture
factor `sum_k a_k K_fk(c)/L_k(c)`, even when it has many nonzero deposits. A numerator cancels only when
every compatible candidate has the same field response. One object per unit, or cancellation before the
alignment/path mixture, is not that model. This is a linear capture model on the conserved-share objects;
its agreement with the simulator's best-probe capture law remains a measurement, not an identity.

Two gates in the submitted design can therefore conflict: arbitrary-field equality to the current scorer
and normalized candidate-specific capture. The appropriate regression contract is equality at a flat field
for all candidates, and arbitrary-field equality only when the candidate field responses coincide. Before
any production change, quantify how often they differ on the test chromosome and real annotation.

## Two internally consistent objectives

For the conditional fragment-mixture model, a structured variational posterior retains assignments
conditional on the field:

    F(a,q) = E_q[sum_f log sum_k a_k K_fk(c)/L_k(c)] - KL(q || p).
    r_fk(c) = a_k K_fk(c)/L_k(c) / sum_j a_j K_fj(c)/L_j(c).

The abundance update uses `E_q[r_fk(c)]`. A mean-field coordinate update for `q_o` uses the conditional
expectation of the **whole log-sum likelihood**, not just averaged component totals. Exact enumeration is
the reference. It is legitimate to approximate this objective numerically, but it must be the same
objective in both blocks. SQUAREM may accelerate a fixed objective; changing the field during an abundance
extrapolation requires a new acceptance check.

Alternatively factorize the assignment/object latent variables from the field. With an auxiliary deposit
object, the complete-data probability is proportional to `a_k phi_fkh d_fkh,o c_o / L_k(c)`. The valid
coordinate updates then are

    q_f(k,h,o) proportional to a_k phi_fkh d_fkh,o exp(E_q[log c_o] - E_q[log L_k])
    log q_o(c_o) = log p_o(c_o) + m_o log c_o - sum_k N_k E_(q_-o)[log L_k(c)] + constant.

Here `m_o` is an **expected latent assignment count**, not a fixed observed count, and `N_k` is the same
variational assignment's component total. This is a valid, cheaper lower bound. It is not integrated
responsibility inference and must be priced against the structured reference, especially on correlation
and ambiguous-annotation tests. No choice between these two is justified by naming both updates “EM”.

## Decisions to earn with measurements

- Objective: test the submitted hybrid against the two consistent objectives on exact small models.
- Expansion: test it using exact log-length moments before allowing an additional delta approximation.
- Prior: changing a unit must change neither likelihood nor prior probabilities. Centering on a background
  mode changes the model if it changes relative support; it is not a harmless choice of unit. No library
  reference will be introduced. Examine sensitivity along unidentifiable directions explicitly.
- Junction: retain the current physical relation for the first comparable arm. A free junction field is
  another mechanism and can confound splicing with capture; it needs its own truth experiment.
- Reader: keep production intact until a complete replacement passes. Its removal is an outcome of a
  successful A/B, not a prerequisite for an isolated prototype.
- Identification: measure the likelihood's rank/ridges on the enumerable locus, test chromosome and real
  annotation. A prior-selected solution must be labeled as such.

## Run ledger

Initial exact-model audit: budget under one minute, one process, under 1 GB. No native build and no timing
claim. The submitted toy's constants are retained as diagnostic fixtures, not proposed production knobs.
Results and subsequent budgets will be appended before advancing to a larger instrument.

## Repair of the costly field step: exact augmentation, not the inconsistent bound
+
+There is a third implementation route that retains the integrated allocation objective while using
+component/object sufficient statistics in the field step. It augments the model rather than approximating
+its log-sum. This derivation precedes the prototype.
+
+Given the latent component, path and deposit-object assignments, let their integer totals be `N_k` and
+`m_o`. The field-dependent complete likelihood is
+
+    p(c) product_o c_o^m_o product_k L_k(c)^(-N_k).
+
+For every `N_k > 0`, the gamma integral gives the exact identity
+
+    L_k^(-N_k) = [1/Gamma(N_k)] integral_0^infinity u_k^(N_k-1) exp(-u_k L_k) du_k.
+
+The three Gibbs blocks are therefore:
+
+1. Assign each fragment once to `(k,h,o)` under its normalized `a_k phi_fkh d_fkh,o c_o/L_k` weights.
+   Identical fragments can use a multinomial draw of their multiplicity. Keep integer `N_k`, `m_o`.
+2. Draw `u_k ~ Gamma(shape=N_k, rate=L_k(c))`; when `N_k=0`, use `u_k=0` exactly (no epsilon or fake count).
+3. For linear free-object fields, draw all objects conditionally independently on the existing grid with
+
+       log p(c_o | rest) = log p_o(c_o) + m_o log c_o - c_o sum_k u_k A_ko + constant.
+
+This is an exact posterior sampler for fixed abundances, including cross-object posterior dependence.
+Average the conditional responsibilities to estimate the integrated allocation (Rao–Blackwellization),
+rather than reporting the internal categorical draws. Internal Monte Carlo is distinct from the requested
+fractional final assignment. Existing seed semantics can make the inference reproducible.
+
+If junction weights remain deterministic functions of the field, use the same augmentation but update
+one free object on its grid using `sum_j m_j log c_j(c) - sum_k u_k L_k(c)`. The cheap independent update
+then becomes a local coordinate update; it still never revisits fragments. This retains the current
+junction relation for the first comparable arm instead of introducing a second mechanism.
+
+The cost is not automatically 8 or 16 E-steps: mixing and Monte Carlo error must be measured. This route
+removes the incorrect second-order expansion and the factorized-posterior assumption, but does not resolve
+identifiability or missing candidate-path data. Test its stationary distribution against exact enumeration,
+the common/independent/opposed identity, zero-count components and grouping/repetition before a panel arm.
+Budget for that sampler audit: under two minutes, one process, under 1 GB; no performance claim.

### Sampler screen: exact conditionals do not imply useful mixing
+
+The gamma-augmented sampler targets the correct distribution but failed the first practical gate.
+Four starts, 32,768 sweeps each, gave component-one counts 133.056, 127.978, 130.589 and 109.787 versus
+the exactly enumerated 150.836. The first chain's batch standard error was only 0.00057 counts: a chain
+trapped in a mode can report misleading precision. No panel run is justified with this sampler.
+
+The same posterior admits an exact **collapsed** field update: integrate the gamma variables back out,
+retain the sampled integer assignments, and draw one field coordinate from
+
+    log q_o(g | c_-o,z) = log p_o(g) + m_o log c_g - sum_k N_k log L_k(c_-o,c_g).
+
+This is not the submitted hybrid: assignments are sampled jointly with the field, not fixed to an
+averaged responsibility. It preserves all cross-object dependence. Its cost is objects × grid ×
+components, with no fragment loop inside the field step. Add an exact common-scale orbit update:
+among all grid states proportional to the current field, sample with probability proportional to
+their prior masses. The likelihood is identical along that orbit; this introduces no proposal scale
+or acceptance constant and avoids trapping the arbitrary field unit. For a deterministic junction
+relation homogeneous in the field, the same invariance holds.
+
+Next numerical budget: under two minutes, one process, under 1 GB. Compare collapsed updates, with and
+without the orbit update, to the same exact posterior from all four starts. This changes a numerical
+kernel, not the statistical model; no panel-based parameter selection is involved.

## Identification on annotation, measured independently of the fitted answer
+
+For a path-complete kernel tensor `B_fko`, let `K_fk = sum_o B_fko c_o`, `A_ko = sum_f B_fko`,
+and `p_f = sum_k a_k K_fk/L_k`. At a flat field, the derivative in a log-field coordinate is
+
+    J_fo = sum_k a_k [B_fko/L_k - K_fk A_ko/L_k^2].
+
+For a two-component abundance logit its derivative is `a_1 a_2 (K_f1/L_1 - K_f2/L_2)`.
+Project that derivative onto the field Jacobian's columns, weighting rows by `1/sqrt(p_f)`.
+Zero residual means a local abundance change is indistinguishable from changing the capture field;
+a nonzero residual identifies that direction in this model. This is a rank diagnostic, not an
+accuracy score. Enumerating fully observed paths gives more information than partially sequenced
+fragments and is therefore a favorable identification test. Selection is by annotation structure,
+never by a quantification error or a panel outcome. No-gDNA, same-strand RNA is the deliberate test
+of the zero-gDNA identification claim.

## Verdicts on the six questions
+
+The follow-up review is complete. The reference implementation below is built and checked outside the main
+tree. It is **not a production pipeline implementation or a release candidate**. The original plan cannot be
+ported literally: its numerator/normalizer and unit-representation contracts disagree with the source.
+
+| Question | Verdict and measurement |
+|---|---|
+| 1. EM bound or exact log-sum? | Reject the hybrid updates. On the submitted toy with noiseless 50/50 truth, exact integration gives 50.101/49.899%, consistent structured mean-field 50.045/49.955%, and the hybrid 40.614/59.386%. The hybrid reduces its stated objective on 12/16 seeded fixtures, worst −0.051942. Use the exact-assignment **collapsed field Gibbs** update derived above as the reference. It needs integer component/object totals in the field step and preserves cross-object dependence. |
+| 2. May the second-order expansion ship? | Not this expansion under this gate. Even supplying **exact** log-length moments gives responsibility error 0.036033 for lengths {1,9}, 0.154981 for {1,30}, and a negative responsibility −0.345210 for {1,1000}, whose exact answer is 0.650934. Clipping would be another approximation. The proposed S=16 reference is not precise either: the {1,9} fixture has independent-sample standard error 0.062098. A covariance formula for lengths does not give log-length covariance. |
+| 3. Which prior? | Retain the existing landscape's probability masses for the first isolated model comparison. Use its native log-density grid with a numerical translation applied to **both** kernel and lengths; a median or mode used only as a unit cannot change the inference. Do not anchor a new prior to a background reference. The C^-2 proposal was authorized as an isolated experiment, not selected for production. The prior is substantive: the same toy data gives T1=0.501010 under the submitted prior versus 0.405347 under uniform grid masses. Its role must be reported, especially where the likelihood is deficient. |
+| 4. Free junction or hard cap? | The hard cap is **not a universal physical constraint**. On the existing simulator's split-probe fixture, noiseless junction average 44.559322 exceeds the largest touched gDNA-object average 27.779661; the current capped rule underprices it 37.6569%. A full RNA fragment weighs 61 versus a maximum gDNA fragment's 41. A free junction is a legitimate latent object, but it adds a confounding direction and must be compared separately after the joint-field arm. Keep the current relation only as a first-arm comparator, not as a truth falsifier. Replace the cap gate by simulator-yield and allocation truth. |
+| 5. Delete the reader now? | No. It remains the production comparator. Deletion follows a successful full replacement A/B and coverage of its published surfaces, including objects outside EM loci. Nothing in this session changes production inference. |
+| 6. Does isoform identification survive annotation? | Sometimes; it is not guaranteed. The toy has a continuous equal-likelihood abundance/field ridge (maximum probability difference 1.1e-16). In six structurally selected test-chromosome pairs, five abundance directions are confounded with the field to numerical precision; the sixth has a 0.249 projection residual fraction. In six selected real-annotation pairs, one is confounded (5.2e-16), while five have residual fractions 0.798–1.000. These are local-rank results at a flat field using **fully observed paths**, an optimistic experiment, not VCaP accuracy or a population prevalence estimate. |
+
+The initial gamma augmentation is mathematically valid but practically refused: its four chains miss the exact
+component count by 17.78–41.05 fragments at 32,768 sweeps. Collapsing those gamma variables works on this gate:
+four starts at 16,384 sweeps recover 150.835980 expected fragments within 0.065497, without the orbit update.
+The orbit update also passes, but is not needed for this fixture. This is a numerical implementation result,
+not a whole-genome cost forecast. The exhaustive four-state transition matrix preserves the exact posterior
+to 5.6e-17; stale-coordinate and omitted-normalizer defects miss it by 0.067807 and 0.499019 respectively.
+
+### Production representation and normalization findings
+
+One object per fragment is not the conserved rule. Before draining the test chromosome's deferred records,
+object-count incidence / deposited-fragment count is 1.034261 stranded OFF, 1.069435 stranded ON, 1.033909
+unstranded OFF and 1.069722 deferred ON. The corresponding **conserved mass ratio is exactly one**. There are
+1,580–2,246 held placements with multiple paths on these four g50 rows. Counting their object incidences as
+independent observations would reuse evidence; picking just one object would change the kernel.
+
+After the production drain, all 24 ladder/gap conditions still conserve their full fragment mass exactly.
+The incidence ratio ranges below describe structure, not quantification errors:
+
+| panel | stranded OFF | stranded ON | unstranded OFF | deferred ON |
+|---|---:|---:|---:|---:|
+| ladder, range across rungs | 1.042–1.668 | 1.388–1.824 | 1.042–1.667 | 1.388–1.824 |
+| RNA-long | 1.410 | 1.501 | 1.410 | 1.501 |
+| RNA-short | 1.085 | 1.286 | 1.084 | 1.286 |
+
+The implemented `path_response` conserves these shares while giving each fragment one likelihood factor.
+Exhaustive placement sums agree with `TranscriptObjects`' length contraction to 2.95e-16 relative across
+all 918 test transcripts (19,336,758 placements at three diagnostic widths). On 16 deterministic structural
+quantiles of VCaP's real annotation, the error is 4.40e-16 (1,321,546 placements). The full annotation has
+457,371 transcripts, 10,558,249 pieces and 10,100,878 cut rows; 424,441 transcripts have multiple cuts.
+Its per-transcript piece-count median/p90/p99/max is 14/53/128/988. These counts also qualify the design's
+memory/cost assumptions. No sequencing-data accuracy was measured in this annotation audit.
+
+The scorer deliberately uses `gdna_realized_pmf`, whereas the opportunities in the gDNA length use
+`gdna_pmf`. On captured stranded g50 rows, geometry/scorer means in bp are:
+
+| panel | geometry | scorer | total variation between PMFs |
+|---|---:|---:|---:|
+| ladder | 216.697715 | 239.112613 | 0.112461 |
+| RNA-long | 78.759388 | 80.284799 | 0.041831 |
+| RNA-short | 249.803091 | 258.088516 | 0.062509 |
+
+For a 100-bp overlap window, directly summing the proposed scorer kernel over its admissible starts and
+dividing by the proposed geometry length gives 1.071001, 1.008581 and 1.023754, respectively, instead of one.
+This tests the **literal new generative identity**, not the old ruler's raw L/Y ratio. It does not establish
+that replacing a production PMF or length alone improves allocation. That swap was already measured badly
+in earlier work. A production joint model must specify a common measure for its kernels and normalizers;
+it cannot assert that the current two-frame construction is already that normalized generative model.
+
+There is another retained evidence path: `em_pseudocounts` uses calibration's estimated gDNA share, with
+one pseudo-observation per gDNA-eligible fragment. The proposed landscape and this grouped prior are both
+empirical Bayes on the library. Saying “one likelihood factor per fragment” is correct for the new kernel;
+saying it makes every part of the unchanged pipeline free of evidence reuse would not be correct. Keep
+the grouped prior fixed in the first comparison, and disclose it rather than silently replacing it.
+
+## The single mechanism I would take forward
+
+**The shared field with candidate-path kernels and collapsed posterior sampling**, not another independent
+weight readout. Its executable core is `shared_field.py`; `path_geometry.py` evaluates the conserved kernel.
+`fit_exact` is the allocation reference for enumerable loci. The core deliberately has no detector, new
+reference, empirical shrinkage coefficient, panel fit, production switch or reader dependency.
+
+The implementation contract is now concrete:
+
+1. Preserve a candidate's path until scoring has applied capture. RNA candidates must retain the path chosen
+   by the scorer; gDNA's multiple-hit mixture needs its eligible paths, not a single averaged scalar. Store
+   sparse object/share vectors and multiplicities. Retain deterministic RNA fragments' paths too. An `int32`
+   per unit is insufficient. Group only identical candidate/path likelihoods and deposit vectors.
+2. Build `A_ko` from the same kernel measure used by `B_fko`. The existing RNA conserved-share generator is
+   now gated by actual path enumeration. Resolve the gDNA observed-law/opportunity-law contract explicitly
+   before connecting the prototype; do not hide a PMF swap in a field update.
+3. For fixed abundance weights, draw each unit's component/path/object once. Update a field coordinate on
+   its existing grid from `prior + m_o log c_o - sum_k N_k log L_k(c)`. These sampled sufficient statistics
+   retain field dependence; replacing them with mean counts is the rejected hybrid. Average conditional
+   responsibilities for fractional output. A coupled junction relation uses the generalized conditional
+   already derived above; it does not use the free-object shortcut.
+4. Retain the existing abundance/prior block for the initial comparison, freezing its inputs during each
+   field update. A posterior sampler tested at fixed abundances is not by itself a tested stochastic EM.
+   The production convergence rule must account for integration error and different starting chains;
+   a fixed seed, a small within-chain standard error or an entropy cut is not that proof.
+5. Compare the first coherent pipeline candidate to shipping source on the test chromosome, then every
+   ladder rung and both gap directions, then VCaP. Only then separate the free-junction arm. Read transcript,
+   gene, false-positive and population errors per stratum; keep calibration-vs-oracle beside them. Counts
+   from this mathematical reference are not substitutes for that A/B.
+
+Falsifiers: omitted normalizer, candidate numerator cancelled when paths differ, lost cut shares, duplicate
+observations, stale-coordinate field updates, independent length averages, prior not transformed with the
+unit, sensitivity along a likelihood ridge mislabeled as identification, a trapped chain with tiny reported
+error, and any captured contaminated regression on either gap panel. A negative responsibility immediately
+rejects an approximation. A junction exceeding adjacent gDNA averages does **not** reject the model.
+
+## What is refused, and what is still unfinished
+
+Refused: the submitted hybrid update; the second-order production shortcut on broad posteriors; one object
+per unit; universal numerator cancellation; S=8/16 treated as a validated integration budget; a hard adjacent
+object cap presented as physical truth; background-mode anchoring; reader deletion before a replacement
+passes; and any 0.8.0 readiness claim from these mathematical or annotation tests.
+
+Unfinished: candidate-path data in the production scorer, a common kernel/normalizer measure for gDNA, the
+whole-locus abundance/field solver with measured numerical convergence and cost, and its full accuracy A/B.
+Consequently **no new end-to-end ladder/gap/VCaP allocation arm was launched or selected**. Running an arm
+with the demonstrated representation/normalization defects would not test the derived mechanism. The
+standing release regression and the standing 2,837-test suite have not been relabeled as fixed or rerun.
+
+The reviewed design resolves the six questions and supplies a working mathematical reference, not a claim
+that the remaining implementation is already finished. The viable part of the idea has survived; the
+submitted route to production has not. Production work should start at the representation/measure boundary
+above, rather than implementing the known-invalid E-step or deleting its comparator.
+
+## Execution receipts and verification
+
+All code in this review was prototyped outside the main tree. Only sandbox documentation is changed in the
+main checkout, and the owner-committed ruler fix at `624dfade` is untouched. The original design is archived
+unchanged beside the reference code. The exact-model audit took under a second; the first sampler audit
+about 10 seconds; the collapsed sampler audit about 30 seconds; test/real geometry audits about 0.52/8.95
+seconds; annotation identification about 16 seconds; all 24 input-frame audits 134 seconds (two-minute
+budget extended to four while running, observed RSS about 1.6 GB). These are budget receipts, **not idle-machine
+performance claims**. The model core has 13 passing focused gates, plus exhaustive and source-mutated
+sampler checks. Repository documentation checks are recorded in the evidence archive.

A further qualification to the submitted “blindnesses”: zero observations at an object with positive
opportunity are negative evidence through the component normalizer. Only a field coordinate that actually
cancels from the likelihood retains its prior unchanged. The added silent-object gate distinguishes these
cases. Likewise the deferred stratum's outcome cannot be declared unchanged just because its strand channel
is absent: it still has length and path information. It remains deferred, and its outcome needs measurement.
