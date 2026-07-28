# PN2D Task 12 branch audit and Task 14 same-support closure

Date: 2026-07-27

## Outcomes

- Task 12: `nonsmooth_branch_derivative`.
- Task 14: `exact_high_bias_oracle_available`.
- Nonlinear authorization: stopped. Tasks 15-19, coarse self-consistent BV,
  and fine PN2D remain unauthorized because Task 12 did not return
  `source_jacobian_dependency_identified_and_closed`.
- Production decision remains `keep_opt_in_diagnostic_no_default_change`.

## Task 12 direct source-only audit

The block audit no longer obtains the canonical element-edge avalanche
Jacobian by subtracting two full continuity Jacobians. The assembler now
provides direct APIs for:

- the shared Real impact-ionization source residual;
- the nine-slot Tri3 local-forward-AD source Jacobian; and
- an independent cell-local symmetric-FD source Jacobian that scatters each
  cell contribution into the global carrier rows.

The analytic and FD paths use the same canonical configuration fingerprint:

`3ff38a69c6e14b011668bcbd5984ba2f7034d276338ff34b1d69f98bc893b2bd`

The frozen coarse7x3 `-20 V` provenance is:

| Input | SHA-256 |
|---|---|
| state | `895cb75415c40058dbfab5f06710120401d54f79f58c6b4d6b16f940e585d3c2` |
| mesh | `c9aaf5f3130f2e1e78e399d155390ed8f19a306ff9ab5af4904230b5e328bc7e` |
| doping | `714bb5c461d0acba49b1f9211318cc120a2e3891f92367776b636eda4b7fd155` |
| materials | `212b896ed86ace76d0d02e86b90f6ced6a2851417b2ec321d0dd702cfbfa6524` |

Generated evidence is under
`build-release/pn2d-task12-direct-source-20260727` and is not committed.

### Frozen results

| State | Relative step | Analytic norm | FD norm | Difference norm | True relative |
|---|---:|---:|---:|---:|---:|
| coarse7x3 `-10 V` | `3e-11` | `8.853246323e-15` | `3.424605916e-13` | `3.423828847e-13` | near-zero absolute pass |
| coarse7x3 `-20 V` | `1e-8` | `3.668331191e-12` | `3.667520640e-12` | `3.896966581e-14` | `1.062326812e-2` |
| coarse7x3 `-20 V` | `1e-10` | `3.668331191e-12` | `3.668330319e-12` | `6.746090352e-18` | `1.839007985e-6` |
| coarse7x3 `-20 V` | `3e-11` | `3.668331191e-12` | `3.668328229e-12` | `1.180726298e-17` | `3.218701466e-6` |

The `-10 V` near-zero absolute gate closes at `3e-11`. The nonzero `-20 V`
gate does not reach `1e-8`; its best recorded result is `1.84e-6`.

### Branch classification

At `-20 V`, eight mesh edges have exactly flat hole QFP. The smallest nonzero
hole-QFP edge difference is `4.227729278e-13 V`. Perturbations large enough to
remain numerically resolvable cross the `abs(QFP difference)` branch used by
QFP-gradient mobility; perturbations small enough to remain on one side enter
double-precision subtraction noise. The error decreases with the outer step,
then reverses below the optimum. This is a branch-resolution limit, not a
missing matrix scatter, scaling factor, or relaxed-tolerance opportunity.

The frozen Task 12 contract requires this case to be typed
`nonsmooth_branch_derivative`; it must not be made GREEN by changing the
`1e-8` gate.

## Task 14 same-support closure

The accepted `9 branches x 8 biases x 2 roots = 144` matrix was reused without
rerunning Sentaurus. The new analyzer verifies that the two implicit-default
runtime record streams are identical, maps native vertex and element records
to the same coarse7x3 Tri3 cells, and emits:

- `same_cell_process.csv`;
- `support_class_summary.csv`;
- `fixed_hotspot_process_chain.csv`; and
- `same_support_summary.json`.

Generated outputs are under
`build-release/pn2d-task14-same-support-20260727` and are not committed.

The fixed `-20 V` generation hotspot is interior cell 13 with nodes
`11;10;8`. On that same support, the `-20 V / -19 V` ratios are:

| Stage quantity | Ratio |
|---|---:|
| mean electron density | `3.77935` |
| electric field | `1.02776` |
| electron QFP gradient | `1.02908` |
| electron mobility | `0.97261` |
| electron current density | `5.69013` |
| electron alpha | `1.14342` |
| total generation | `4.70196` |

Using the frozen material-change threshold of `1.10`, the first material stage
in dependency order is density. Field/QFP-gradient and mobility are not the
first departure. All generation-active cells at `-19 V` and `-20 V` are in the
interior class; contact-adjacent cells remain outside the active set.

This closes the Task 14 same-support, contact/interior, active-region, and
hotspot-coincident exit gate. It does not override the independent Task 12
nonlinear stop.

## Verification

- `test_newton_solver`: 728 assertions in 63 cases passed.
- `test_impact_ionization`: 523 assertions in 42 cases passed.
- Task 11/14 Python regression selection: 14 tests passed.
- `ascii_sources`: passed.
- `git diff --check`: passed.