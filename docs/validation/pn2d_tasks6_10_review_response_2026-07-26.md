# PN2D Tasks 6-10 review response and stop ledger

Date: 2026-07-26

## Scope completed

- Task 6: constrained-obtuse RED and opt-in area-conservation patch completed; the corrected true-relative source Jacobian gate failed.
- Task 7: exact imported-state residual, paired controls, carrier-only/coupled first updates, topology/provenance, and independent v2 verification completed.
- Task 8: not entered; Task 6 Jacobian and Task 7 cross-topology causality gates failed.
- Task 9: not entered because Task 8 was not authorized.
- Task 10: decision ledger and independent review responses completed; no default-change proposal.

## Independent review result

Scientific and code reviews independently agreed on the primary blocker: the original full-continuity finite-difference test did not isolate avalanche derivatives. The scientific review additionally required paired controls, independently derived authorization, exact provenance, honest zero-vector direction, first causal node/term, and limited claims for obtuse support distribution. The code review required failed gates to fail verification and complete topology/configuration sealing.

## Responses

| Finding | Response |
|---|---|
| full transport masked the `1e-8` source gate | focused test now compares analytic impact-only difference with central FD of independent source replay |
| full forward-FD matrix subtraction polluted real-state blocks | diagnostic source/recombination blocks now central-difference isolated term diagnostics and remove gauge-activation differences |
| double-subtracted baseline made real-state evidence transport-dominated | isolated pairs are passed directly; source blocks are no longer transport-dominated |
| verifier preselected negative outcome | v2 verifier independently derives first-material directions, Jacobian gates, authorization, and outcome |
| `rel_diff` used an absolute floor of one for small source norms | denominator is now the larger analytic/FD norm; fresh maximum is `0.9640948767506723`, so v2 verifier correctly returns `pass: false` |
| controls existed only for opt-in | production and opt-in avalanche-off/SRH-off controls are paired at every topology/bias/node/carrier |
| topology/state/config provenance incomplete | verifier checks both 117-config roots, complete input-hash keys, sealed path bindings, exact coarse `psi/eQFP/hQFP`, mesh nodes, contacts, and element permutations |
| zero-reference update direction overstated | recorded `undefined_zero_reference`; absolute delta is only an error magnitude |
| first causal node/term absent | v2 manifest records the first maximum avalanche departure for each topology |
| row scale was analyzer-defined | renamed `diagnostic_incident_term_scale`; actual continuity unit scale recorded separately |
| obtuse local partition lacks native oracle | claim limited to nonnegative exact-area closure; helper stays opt-in |

## Final typed outcomes

- Task 6 authorization: `jacobian_gate_failed` (maximum nonzero true-relative `0.9640948767506723`; no tolerance change).
- Task 7: `operator_improvement_without_qfp_causality`.
- Tasks 8-9: `not_entered_task6_jacobian_and_task7_causality_gates_failed`.
- Task 10 production decision: `keep_opt_in_diagnostic_no_default_change`.

No generated simulation root is committed. The user-owned untracked `tmp/` remains untouched.
