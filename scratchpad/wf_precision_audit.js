export const meta = {
  name: 'precision-honesty-audit',
  description: 'Root-cause the honest-precision violation: a tiny-count source becomes a high-precision message',
  phases: [
    { title: 'Trace', detail: '4 independent instrument-and-trace audits of the precision flow' },
    { title: 'Verify', detail: '2 adversaries stress-test the leading root-cause hypotheses' },
    { title: 'Adjudicate', detail: 'synthesize the precise root cause + the exact code locus + fix direction' },
  ],
}

const ENV = 'source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate rigel && OMP_NUM_THREADS=1'
const REPO = '/Users/mkiyer/proj/rigel'

const SUBSTRATE = `
Rigel calibration pass-0 — the belief-propagation solver (bp_solver._unified_solve) deconvolves each genomic
node's unspliced fragment mass into (f_pos, f_neg, f_g) by a forward-backward message-passing sweep over a
region↔boundary chain. Files (do NOT read docs/calibration/archive/):
- ${REPO}/src/rigel/calibration/bp_solver.py — THE SWEEP: the relay (_relay), the combine (_transport), the
  reframe r = ρ_tot(dst)/ρ_tot(src), σ²_transfer (transfer_logvar/logvar_tot), ÷M_dst, and the THREE precision
  streams (full pg/pp/pn → mode fusion; measurement mg/mp/mn → gdna_imp/rna_imp; composition τ → the λ-message).
- ${REPO}/src/rigel/calibration/node_init.py — the OWN per-node precisions: own_precision(n,v_log,live) =
  n/(n·Var(log f)+1) = 1/(Var(log f)+1/n); tau_lam (the Schur composition precision). The count n is the
  discrete fragment count.
- ${REPO}/src/rigel/calibration/enrichment_frame.py — composition_logvar (= 1/n + comp), message_precision,
  reframe_density, density_mode_logfrac. The message-VARIANCE laws (transport_seed_logvar/graft_rna_logvar/
  peel_rna_logvar/transfer_logvar). The MC arbiter: ${REPO}/scripts/debug/message_variance_mc.py.
- ${REPO}/src/rigel/calibration/node_geometry.py — n_unspl_left/right (the discrete COUNTS), the eff-lengths
  (discrete base pairs), node_total_density, node_global_geometry.
- ${REPO}/docs/calibration/SESSION_2026_07_24_HANDOFF_4.md §6 — the FINDING + the node-1909 dissection.
- ${REPO}/docs/calibration/message_variance_derivation.md — the settled variance laws (M1-M5).

THE OWNER'S INVARIANT (the lens for this audit). **Precision is fundamentally constrained by (1) discrete
COUNTS over (2) a discrete LENGTH in base pairs.** Enrichment ratios (up to 1000×) scale the MODES so
compositions are comparable across a capture cliff — but **enrichment must NEVER scale the PRECISION.** A node
with tiny counts is fundamentally imprecise; its message must stay imprecise no matter how the reframe stretches
its mode. THE BUG: somewhere a tiny-count source becomes a HIGH-precision message.

THE ANCHOR CASE (already dissected). Scenario gdna_gdna300_ss_0.99_nrna_present_capture_on, node 1909 (an exon,
oracle f_g=0.985, mostly gDNA, mass 66,544). Its own message-free solve fg_loc=0.953 is CORRECT but WEAK
(τ=1.6). It receives a STRONG WRONG RNA+ measurement cm_p=26.45, mode f_pos=0.718 (claims 72% RNA+ for a
~1.5%-RNA node) → f_g collapses to 0.51. Its neighbours are tiny boundaries (mass 6, 14). Where does cm_p=26.45
come from, and is it justified by the DISCRETE COUNTS at its origin?

TOOLS (rigel env, ${ENV} python <script>):
- ${REPO}/scratchpad/dump_node.py — dumps a node's own belief + all message modes/precisions + neighbours
  (edit it to instrument more; the solver's _capture dict already publishes cm_g/cm_p/cm_n/c_tau/mo_*).
- ${REPO}/scratchpad/trace_stranded.py — runs calibrate with σ²_transfer ON vs OFF (env RIGEL_S2T_OFF) and
  ranks nodes by regression.
- RIGEL_S2T_OFF=1 env var zeros σ²_transfer inside the solver (a diagnostic hook already in bp_solver).
- ${REPO}/scripts/debug/pass0_node_dissect.py, pass0_oracle_bench.py.

METHOD (mandatory): INSTRUMENT and TRACE to NUMBERS — do not theorise. The invariant is FALSIFIABLE: for a given
message, compare precision_CLAIMED vs the discrete COUNT available at its origin (Poisson: precision ≤ n). Find
exactly WHERE claimed > available, and WHY (which line inflates it). Write your instrumentation to
${REPO}/scratchpad/audit_<n>.py, RUN it, and report measured provenance.
`

phase('Trace')
const AUDIT_SCHEMA = {
  type: 'object', additionalProperties: false,
  required: ['ran', 'root_cause', 'invariant_violation', 'code_locus', 'evidence', 'fix_direction', 'confidence', 'summary'],
  properties: {
    ran: { type: 'boolean', description: 'did your instrumentation actually run in the rigel env' },
    root_cause: { type: 'string', description: 'the precise mechanism by which a tiny-count source becomes a high-precision message (or a statement that you could not confirm one)' },
    invariant_violation: { type: 'string', description: 'the SPECIFIC place where precision_claimed exceeds the originating discrete count, with the measured numbers (e.g. cm_p=26.45 at node 1909 traced to N spliced counts = M)' },
    code_locus: { type: 'string', description: 'file:line(s) where the inflation happens (or the missing floor)' },
    evidence: { type: 'array', items: { type: 'string' }, description: 'measured numbers from your instrumentation supporting the root cause' },
    fix_direction: { type: 'string', description: 'the fix that would restore honest precision (respect originating counts/length; enrichment scales modes not precisions)' },
    confidence: { type: 'string', enum: ['high', 'medium', 'low'] },
    summary: { type: 'string' },
  },
}
const angles = [
  'PROVENANCE: instrument the relay to record, per accumulated message precision at node 1909, the origin (which node, which discrete spliced/count) of every contribution and how σ²_transfer damped it hop-by-hop. Does cm_p=26.45 trace to ≥26 real spliced counts at 1909\'s own junctions, or is it summed/propagated from distant nodes whose counts do not belong to 1909? Quantify.',
  'INVARIANT SCAN: across ALL nodes of the anchor scenario, compute per-message precision_claimed (cm_g, cm_p, cm_n, c_tau, and the fused ψ precisions) vs the honest Poisson bound (the discrete count n at the ORIGIN, over the eff-length). Rank the worst violations (claimed ≫ count). Is the violation systematic (tiny-mass sources → high precision)? Quantify the distribution.',
  'ENRICHMENT LEAK: audit every place a precision is FORMED or TRANSPORTED (own_precision; the graft SP/(1+SP·s2t) term; _damp; _dv; ÷M / message_precision; the additive measurement fuse mg[i]=mg_own[i]+tmg). Prove whether the reframe r or the ÷M or the mode inflation leaks into the PRECISION (it should scale only the MODE). Construct a minimal numeric example where a 1000× reframe leaves the precision too high.',
  'PROPAGATION SEMANTICS: a spliced MEASUREMENT at junction J measures J\'s RNA. Trace whether the measurement stream (mp) PROPAGATES J\'s measurement to distant nodes and ADDS it (additive fuse) as if it measured them — inflating precision at nodes far from J. Is the measurement precision decaying with distance as honest imprecision requires, or accumulating? Test with σ²_transfer ON vs OFF.',
]
const audits = await parallel(
  angles.map((a, i) => () =>
    agent(
      `${SUBSTRATE}\n\nYou are audit agent #${i + 1}. FOCUS: ${a}\n\n`
      + `Instrument (${REPO}/scratchpad/audit_${i + 1}.py), RUN it in the rigel env, and report the MEASURED root cause `
      + `of the honest-precision violation with numbers. Default to "not confirmed" if your instrumentation does not `
      + `actually demonstrate precision_claimed > originating_count. Name the exact file:line.`,
      { label: `audit#${i + 1}`, phase: 'Trace', schema: AUDIT_SCHEMA, effort: 'high' }
    )
  )
)

phase('Verify')
const VERIFY_SCHEMA = {
  type: 'object', additionalProperties: false,
  required: ['ran', 'confirmed_root_cause', 'refuted', 'best_code_locus', 'residual_doubts', 'summary'],
  properties: {
    ran: { type: 'boolean' },
    confirmed_root_cause: { type: 'string', description: 'the root cause you could REPRODUCE with your own instrumentation' },
    refuted: { type: 'array', items: { type: 'string' }, description: 'any audit claim your instrumentation refuted' },
    best_code_locus: { type: 'string' },
    residual_doubts: { type: 'array', items: { type: 'string' } },
    summary: { type: 'string' },
  },
}
const auditsJson = JSON.stringify(audits.filter(Boolean), null, 1)
const verifiers = await parallel([
  'Reproduce the leading root cause INDEPENDENTLY (your own instrumentation). Then adversarially test it: construct a case that SHOULD trigger it and one that should NOT, and confirm the precision behaves as the root cause predicts. Try to find a simpler/deeper cause the trace agents missed.',
  'Stress the honest-precision INVARIANT itself: is the claim "precision ≤ originating count" exactly right for a MESSAGE (vs the own belief)? Consider the graft (spliced measurement — a real count), the composition τ (count enters only via the strand Fisher), and the fusion of independent sources (precisions legitimately add). Pin down precisely which precisions are honestly-bounded-by-count and which the audits wrongly flagged, so the fix does not over-correct.',
].map((angle, i) => () =>
  agent(
    `${SUBSTRATE}\n\nYou are VERIFIER #${i + 1}. ${angle}\n\nThe trace audits returned:\n${auditsJson}\n\n`
    + `Instrument (${REPO}/scratchpad/audit_verify_${i + 1}.py), RUN it, and report what you could REPRODUCE vs REFUTE, `
    + `with numbers. Name the best code locus.`,
    { label: `verify#${i + 1}`, phase: 'Verify', schema: VERIFY_SCHEMA, effort: 'high' }
  )
))

phase('Adjudicate')
const VERDICT_SCHEMA = {
  type: 'object', additionalProperties: false,
  required: ['root_cause', 'code_locus', 'invariant_statement', 'reproduction', 'fix_direction', 'what_NOT_to_touch', 'open_questions', 'ready_for_design'],
  properties: {
    root_cause: { type: 'string', description: 'the single precise mechanism, agreed + reproduced' },
    code_locus: { type: 'string', description: 'the exact file:line(s)' },
    invariant_statement: { type: 'string', description: 'the honest-precision invariant, stated precisely enough to implement + test (which precisions are count-bounded, and how enrichment must be kept out)' },
    reproduction: { type: 'string', description: 'the minimal command/script that demonstrates the violation with numbers' },
    fix_direction: { type: 'string', description: 'the fix that restores honest precision (respecting original counts/length; enrichment scales modes not precisions)' },
    what_NOT_to_touch: { type: 'array', items: { type: 'string' }, description: 'precisions that ARE honest and must not be over-corrected' },
    open_questions: { type: 'array', items: { type: 'string' } },
    ready_for_design: { type: 'boolean', description: 'is the root cause precise enough to hand to the derivation/design workflow?' },
  },
}
const verdict = await agent(
  `${SUBSTRATE}\n\nYou are the ADJUDICATOR. Synthesize the 4 trace audits and the 2 verifiers into ONE precise, `
  + `reproduced root cause of the honest-precision violation — the exact mechanism + code locus by which a tiny-count `
  + `source becomes a high-precision message. State the honest-precision INVARIANT precisely enough to implement and `
  + `test (which precisions are count-bounded; how enrichment must be kept out of precision). Give the minimal `
  + `reproduction and the fix DIRECTION (respecting original counts/length). List what must NOT be over-corrected `
  + `(precisions that are honest). Re-run any instrumentation you need to resolve disagreements — do not average opinions.\n\n`
  + `TRACE AUDITS:\n${auditsJson}\n\nVERIFIERS:\n${JSON.stringify(verifiers.filter(Boolean), null, 1)}`,
  { label: 'adjudicate', phase: 'Adjudicate', schema: VERDICT_SCHEMA, effort: 'high' }
)

return { audits: audits.filter(Boolean), verifiers: verifiers.filter(Boolean), verdict }
