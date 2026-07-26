export const meta = {
  name: 'message-variance-verify',
  description: 'Independently verify the pass-0 message-variance derivation + adjudicate the two-message combine fork',
  phases: [
    { title: 'Derive', detail: '3 independent re-derivations from substrate, each MC-validated' },
    { title: 'Adversary', detail: '2 adversarial critics attack the laws + the M6 combine recommendation' },
    { title: 'Adjudicate', detail: 'synthesize into a settled verdict' },
  ],
}

const ENV = 'source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel && OMP_NUM_THREADS=1'
const REPO = '/Users/mkiyer/proj/rigel'

const SUBSTRATE = `
Rigel calibration pass-0 MESSAGE-VARIANCE model. Read these first (do NOT read docs/calibration/archive/):
- ${REPO}/docs/calibration/variance_model_handoff.md   (§1-4: the message ground truth, per-component framework, σ²_transfer)
- ${REPO}/docs/calibration/variance_foundation_proposal.md  (the settled τ_λ FOUNDATION: Var(log f_g)=(1−f_g)²/τ_λ, Var(log f_r)=f_g²/τ_λ)
- ${REPO}/scripts/debug/message_precision_mc.py   (the MC VALIDATION TEMPLATE — k-mode T1/T2/T5, peel T3/T4)
- ${REPO}/src/rigel/calibration/enrichment_frame.py   (the pure transport arithmetic incl. composition_logvar)
- ${REPO}/src/rigel/calibration/bp_solver.py   (_unified_solve: the relay _relay / combine _transport / _fuse / ÷M_dst that consume the message precisions)
- ${REPO}/src/rigel/calibration/simplex_logodds.py   (_local_loglik_logodds: how a message enters ψ as −½·p_c·(log f_c − mo_c)²)

THE SETUP. The unified solver transports PER-COMPONENT densities ρ_c (c ∈ {gDNA g, RNA-total R}; RNA splits +/−
by tilt). ψ consumes, per component, a Gaussian on log f_c^dst with precision p_c, where the message mode is
    mo_c = log( ρ_c^src · r · E_c / M_dst ),   r = ρ_tot(dst)/ρ_tot(src).
Source component log-variances ("transport seed"): gDNA v_g = Var(log f_g)+1/n; RNA-continue v_ν = Var(log f_R)+1/n
(g and ν SHARE the node count M); RNA-splice v_μ = 1/n_spl (measured → composition-CERTAIN). n = node unspliced count.
Two message directions: GRAFT (boundary→exon) the RNA is a SUM ρ_R=ρ_ν+ρ_μ; PEEL (exon→boundary) the RNA-continue
is a DIFFERENCE ρ_ν = ρ_R(x)/r − ρ_μ. σ²_transfer = Var(log r) = Var(log ρ_tot^dst)+Var(log ρ_tot^src).
Counts are POISSON (ω=0); no overdispersion term. NO magic numbers.

METHOD (mandatory): derive → MC-VALIDATE before concluding. Write a pure-numpy MC (mirror message_precision_mc.py:
draw f_g~Beta with the composition variance, M~lognormal with Var(log)=1/n, spliced count ~lognormal 1/n_spl,
and the enrichment r~lognormal with Var(log r)=σ²_transfer; hold the TRUTH fixed = the imputation premise) and
compare predicted vs empirical to <8% rel. Run it with:  ${ENV} python <your_script>.py
`

phase('Derive')
const DERIVE_SCHEMA = {
  type: 'object',
  additionalProperties: false,
  required: ['mc_script_path', 'mc_ran', 'laws', 'combine_finding', 'combine_recommendation', 'combine_rationale', 'unstated_assumptions', 'summary'],
  properties: {
    mc_script_path: { type: 'string', description: 'absolute path to the MC script you wrote' },
    mc_ran: { type: 'boolean', description: 'did your MC actually run to completion in the rigel env' },
    laws: {
      type: 'array', description: 'one entry per derived law',
      items: {
        type: 'object', additionalProperties: false,
        required: ['name', 'formula', 'mc_rel_error', 'confirmed'],
        properties: {
          name: { type: 'string', description: 'e.g. graft-gDNA / graft-RNA-sum / peel-diff / transfer-var / anchor-limit / divM-conversion' },
          formula: { type: 'string' },
          mc_rel_error: { type: 'string', description: 'the empirical vs predicted rel error you measured' },
          confirmed: { type: 'boolean' },
        },
      },
    },
    combine_finding: { type: 'string', description: 'ψ combines gDNA (log f_g) and RNA (log f_R) messages from ONE source as INDEPENDENT Gaussians on the same λ DOF. MC-test: does that reproduce the correct destination composition posterior? Report the ratio (two-message var / joint-truth var).' },
    combine_recommendation: { type: 'string', enum: ['single-lambda-message', 'two-message-with-correction', 'accept-and-defer-to-AB', 'other'] },
    combine_rationale: { type: 'string' },
    unstated_assumptions: { type: 'array', items: { type: 'string' }, description: 'assumptions the model makes that are not obviously stated/justified' },
    summary: { type: 'string' },
  },
}

const framings = [
  'Frame it as the DELTA METHOD on log-densities (sum vs difference), and separately verify the ÷M_dst→Var(log f_c) conversion by checking M_dst cancels in the mode error.',
  'Frame it as a BAYESIAN measurement model: each source component density is a noisy observation; derive the posterior on the destination composition λ and check what precision ψ must receive per message.',
  'Frame it from COUNT-ZERO-INFORMATION first principles: identify every common-mode factor (M_src, M_dst, enrichment r) that cancels in the composition ratio, and derive which noise sources actually survive into the composition — then verify by MC.',
]
const derivations = await parallel(
  framings.map((f, i) => () =>
    agent(
      `${SUBSTRATE}\n\nYou are independent derivation agent #${i + 1}. ${f}\n\n`
      + `Derive the COMPLETE per-component message-variance model: (1) the transport seed per component (the 1/n split: `
      + `composition-only vs composition⊕sampling); (2) the GRAFT (RNA sum, share-weighted) and PEEL (RNA difference, `
      + `u-weighted) laws; (3) σ²_transfer=Var(log r) and its DIRECTION-dependence (does it cancel on the graft? is it `
      + `load-bearing on the peel/anchor?); (4) the pure-gDNA ANCHOR limit (is it finite where the ratio k=ρ_g/ρ_R is `
      + `singular?); (5) the ÷M_dst → precision conversion. THEN critically examine the COMBINE: ψ receives a gDNA message `
      + `on log f_g AND an RNA message on log f_R, both from one source's belief, treated as INDEPENDENT — MC-test whether `
      + `that reproduces the exact joint destination-composition variance, and if not, quantify and recommend. `
      + `Write your OWN MC (unique filename ${REPO}/scratchpad/mc_derive_${i + 1}.py), RUN it, and report measured rel errors. `
      + `Do not read message_variance_derivation.md or message_variance_mc.py — derive independently.`,
      { label: `derive#${i + 1}`, phase: 'Derive', schema: DERIVE_SCHEMA, effort: 'high' }
    )
  )
)

phase('Adversary')
const ADVERSARY_SCHEMA = {
  type: 'object', additionalProperties: false,
  required: ['refutations', 'law_verdicts', 'combine_verdict', 'combine_recommendation', 'mc_script_path', 'summary'],
  properties: {
    refutations: {
      type: 'array', description: 'concrete attacks you MC-tested',
      items: {
        type: 'object', additionalProperties: false,
        required: ['target', 'attack', 'mc_result', 'survives'],
        properties: {
          target: { type: 'string' },
          attack: { type: 'string' },
          mc_result: { type: 'string' },
          survives: { type: 'boolean', description: 'true if the claim SURVIVED the attack (attack failed to refute)' },
        },
      },
    },
    law_verdicts: {
      type: 'array',
      items: {
        type: 'object', additionalProperties: false,
        required: ['law', 'verdict'],
        properties: { law: { type: 'string' }, verdict: { type: 'string', enum: ['CONFIRMED', 'WRONG', 'INCOMPLETE'] } },
      },
    },
    combine_verdict: { type: 'string', description: 'is the two-message combine really over-confident? by how much? does resolution (1) single-λ actually fix it without breaking the anchor?' },
    combine_recommendation: { type: 'string', enum: ['single-lambda-message', 'two-message-with-correction', 'accept-and-defer-to-AB', 'other'] },
    mc_script_path: { type: 'string' },
    summary: { type: 'string' },
  },
}

const myDoc = `${REPO}/docs/calibration/message_variance_derivation.md`
const myMc = `${REPO}/scratchpad/message_variance_mc.py`
const derivationsJson = JSON.stringify(derivations.filter(Boolean), null, 1)
const adversaryAngles = [
  'Attack the LAWS (M1-M5): find a regime (extreme w_μ, tiny counts, f_g→0 or →1, huge enrichment, u≫1 peel) where a law is WRONG or its linearization silently under-states variance (over-confidence). MC every attack.',
  'Attack the COMBINE resolution: pressure-test whether the "single-λ message" fix (resolution 1) actually reproduces the joint AND stays finite at the pure-gDNA anchor AND handles a partial (one-component) message correctly. Try to show it either does not generalize or introduces a hidden singularity/magic number.',
]
const adversaries = await parallel(
  adversaryAngles.map((angle, i) => () =>
    agent(
      `${SUBSTRATE}\n\nYou are ADVERSARIAL verifier #${i + 1}. Your job is to REFUTE, not agree. ${angle}\n\n`
      + `The proposed derivation is in ${myDoc} (read it) and its MC arbiter is ${myMc} (read + run it: ${ENV} python ${myMc}). `
      + `Three INDEPENDENT re-derivations returned this JSON:\n${derivationsJson}\n\n`
      + `Find every disagreement between them and the proposed doc. For each law and for the M6 combine fork, MC-test the `
      + `disputed point yourself (write ${REPO}/scratchpad/mc_adv_${i + 1}.py, RUN it). Default to WRONG/INCOMPLETE if `
      + `uncertain. Give a final per-law verdict and a combine recommendation grounded in YOUR MC numbers.`,
      { label: `adversary#${i + 1}`, phase: 'Adversary', schema: ADVERSARY_SCHEMA, effort: 'high' }
    )
  )
)

phase('Adjudicate')
const VERDICT_SCHEMA = {
  type: 'object', additionalProperties: false,
  required: ['laws_settled', 'combine_verdict', 'combine_decision', 'combine_decision_rationale', 'corrections_needed', 'implementation_go', 'final_summary'],
  properties: {
    laws_settled: {
      type: 'array',
      items: {
        type: 'object', additionalProperties: false,
        required: ['law', 'status', 'note'],
        properties: { law: { type: 'string' }, status: { type: 'string', enum: ['SETTLED', 'NEEDS-FIX'] }, note: { type: 'string' } },
      },
    },
    combine_verdict: { type: 'string', description: 'the consensus on whether/how much the two-message combine is over-confident, with the MC numbers all agents produced' },
    combine_decision: { type: 'string', enum: ['single-lambda-message', 'two-message-with-correction', 'accept-and-defer-to-AB', 'other'] },
    combine_decision_rationale: { type: 'string' },
    corrections_needed: { type: 'array', items: { type: 'string' }, description: 'any fix the proposed derivation doc needs before implementation' },
    implementation_go: { type: 'boolean', description: 'is the derivation solid enough to begin the Step-2 implementation (M1-M5 + 1/n split + struct_lock override)?' },
    final_summary: { type: 'string' },
  },
}
const verdict = await agent(
  `${SUBSTRATE}\n\nYou are the ADJUDICATOR. Synthesize the proposed derivation (${myDoc}), the 3 independent `
  + `re-derivations, and the 2 adversarial critiques into ONE settled verdict.\n\n`
  + `INDEPENDENT DERIVATIONS:\n${derivationsJson}\n\nADVERSARIAL CRITIQUES:\n${JSON.stringify(adversaries.filter(Boolean), null, 1)}\n\n`
  + `Decide: (a) which laws (transport-seed / graft-sum / peel-diff / transfer-var / anchor-limit / divM-conversion) are `
  + `SETTLED vs NEED-FIX; (b) the definitive verdict on the two-message COMBINE over-confidence and what to do about it `
  + `(single-λ message vs correction vs defer-to-A/B) — weigh count-zero-info (never confidently wrong), the anchor `
  + `regularity, no-magic-numbers, and minimal-change; (c) whether to GREENLIGHT the Step-2 implementation. Read the docs `
  + `and re-run any MC you need (${ENV} python ...) to resolve disagreements — do not just average opinions.`,
  { label: 'adjudicate', phase: 'Adjudicate', schema: VERDICT_SCHEMA, effort: 'high' }
)

return { derivations: derivations.filter(Boolean), adversaries: adversaries.filter(Boolean), verdict }
