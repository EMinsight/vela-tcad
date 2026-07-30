# PN2D WP6 branch-resolved source-Jacobian closure

Date: 2026-07-30

Outcome: `source_jacobian_dependency_identified_and_closed`.

## Fixed evidence

The audit reuses the frozen coarse7x3 `-20 V` imported state:

| Input | SHA-256 |
|---|---|
| state | `895cb75415c40058dbfab5f06710120401d54f79f58c6b4d6b16f940e585d3c2` |
| mesh | `c9aaf5f3130f2e1e78e399d155390ed8f19a306ff9ab5af4904230b5e328bc7e` |
| doping | `714bb5c461d0acba49b1f9211318cc120a2e3891f92367776b636eda4b7fd155` |
| materials | `212b896ed86ace76d0d02e86b90f6ced6a2851417b2ec321d0dd702cfbfa6524` |

The state has eight exactly flat hole-QFP edges. The smallest nonzero
hole-QFP edge difference is `4.227729278e-13 V`.

## Compared derivative paths

1. Shared Real source evaluation and solver-used source/scatter records.
2. Production nine-slot Tri3 local forward AD.
3. Ordinary double cell-local symmetric finite difference.
4. Audit-only 50-decimal-digit cell-local symmetric finite difference.

All four paths reuse
`elementEdgeGssLauxAvalancheSourceIntegralsLocal`; the independent numerical
references change only scalar precision and perturbation. The high-precision
reference keeps nonzero branches fixed by selecting a step below the measured
margin. At exact-zero abs/norm branches, symmetric differentiation selects
zero, matching the production semismooth active-set choice.

No smoothing, empirical factor, physics formula, default, or nonlinear solver
behavior changed.

## Step convergence

The current source-Jacobian norm is `3.668331191164069e-8`.

| Reference | Relative step | Reference norm | Difference norm | True relative difference |
|---|---:|---:|---:|---:|
| ordinary double | `1e-8` | `3.6675206398499735e-8` | `3.896966580910351e-10` | `1.0623268123382636e-2` |
| ordinary double | `1e-10` | `3.668330318573974e-8` | `6.746090352222108e-14` | `1.839007985012764e-6` |
| ordinary double | `3e-11` | `3.6683282286866345e-8` | `1.1807262982132756e-13` | `3.2187014658253814e-6` |
| 50-digit branch-resolved | `1e-14` | `3.668331191164067e-8` | `4.1580808147659477e-22` | `1.1335074719484275e-14` |
| 50-digit branch-resolved | `3e-15` | `3.668331191164068e-8` | `4.3981554769767924e-23` | `1.1989526702416226e-15` |
| 50-digit branch-resolved | `1e-15` | `3.668331191164068e-8` | `2.375942778418628e-23` | `6.476903678003707e-16` |

The nonzero `1e-8` relative gate passes at all three branch-resolved steps.
The focused exact-zero test uses the `1e-12` absolute gate.

Configuration fingerprint `ab2dbf93089c7fe3` and active-branch fingerprint
`8e39422ae0ff24ba` are identical across all six probes. Existing independent
source-component, source-record, nodal-scatter, and assembled-block closure
tests retain their `1e-12` gate.

## Classification

The former `1.839007985e-6` mismatch is caused by ordinary double
source-subtraction cancellation and inability to resolve the smallest
nonzero abs branch margin. It is not:

- an incomplete analytic derivative;
- residual/Jacobian formula drift;
- a configuration mismatch;
- a missing scatter contribution; or
- justification for smoothing a physical branch.

Generated evidence is under
`build-release/pn2d-wp6-branch-resolved-jacobian-20260729` and is not
committed.
