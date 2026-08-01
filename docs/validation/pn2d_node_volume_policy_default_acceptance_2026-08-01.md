# PN2D node-volume policy default acceptance — 2026-08-01

## Outcome

The prospective contract completed with typed outcome:

`ready_for_independent_default_policy_reviews`

All executable scientific, geometry, compatibility, forward-IV, determinism,
and Release regression gates passed.  The contract deliberately does not
authorize a production default change.  Independent scientific and code
reviews of the eventual atomic template patch remain required.

Contract SHA-256:
`3a2d65879d1d7446a259afb6f81af9c17e7da2686536e9f6832eb5b5810f6c22`.

## Geometry and compatibility

The focused geometry tests establish:

- mixed-Voronoi areas are opt-in, positive, and conservative on the boundary
  refinement patch;
- an obtuse triangle receives the declared `1/2, 1/4, 1/4` vertex split;
- reversing triangle winding leaves every vertex volume unchanged;
- omitting `mesh_geometry`, providing an empty object, and explicitly selecting
  `barycentric` all preserve the legacy parser result;
- `mixed_voronoi` is selected only explicitly, and an unknown policy fails
  closed.

This acceptance scope is compatible with an atomic PN2D BV template change. It
does not change the low-level parser default, so existing configurations that
omit the field remain barycentric.

## M0 and M2 reverse-bias controls

Both grids were freshly rerun after the unified contract was frozen.  Each
branch was executed twice on the 29-point exact lattice.  IIC and SG/Laux-on
were entered only after the preceding control gate passed.

| Grid | Branch | Two complete runs | IV/state deterministic | Maximum continuity closure ratio |
| --- | --- | --- | --- | ---: |
| M0 | off | 29/29 each | yes | 2.200e-3 |
| M0 | IIC | 29/29 each; identical to off | yes | 2.200e-3 |
| M0 | SG/Laux-on | 29/29 each | yes | 2.200e-3 |
| M2 | off | 29/29 each | yes | 9.138e-5 |
| M2 | IIC | 29/29 each; identical to off | yes | 9.138e-5 |
| M2 | SG/Laux-on | 29/29 each | yes | 1.917e-6 |

The off curves retain the same-grid Sentaurus golden baseline:

| Grid | 28-point log-current RMSE | Maximum error |
| --- | ---: | ---: |
| M0 | 8.876e-5 dex | 1.562e-4 dex |
| M2 | 3.868e-5 dex | 9.223e-5 dex |

The self-consistent SG/Laux-on gates also pass:

| Metric | M0 | M2 |
| --- | ---: | ---: |
| 28-point log-current RMSE | 0.001766 dex | 0.001923 dex |
| 28-point maximum error | 0.003749 dex | 0.004647 dex |
| 11-point knee RMSE | 0.002810 dex | 0.003051 dex |
| V_break, Vela | -19.653 V | -19.390 V |
| V_break, Sentaurus | -19.654 V | -19.391 V |
| absolute V_break error | 0.001 V | 0.001 V |
| V_slope outcome | both cross; difference 0.00082 V | neither crosses |
| non-monotonic exact-lattice intervals | none | none |

The exact IIC/off state and IV equality on both grids excludes collision-
ionization feedback leakage from the control solution.

## 201-point forward IV

Explicit barycentric was run once and mixed-Voronoi twice from the same sealed
forward-IV configuration.  All three runs converged at `201/201` biases.
The two mixed IV files have identical SHA-256 hashes.

At the predeclared `1, 2, 5, 10, 15, 20 V` Sentaurus anchors:

- median mixed-Voronoi relative error: `0.2700%`;
- maximum mixed-Voronoi relative error: `0.4066%`;
- maximum error degradation relative to barycentric: `2.42e-14` absolute.

The node-volume selection is therefore numerically neutral for this forward
IV branch at the reported precision and does not consume the `0.05%`
degradation allowance.

## Release regression

The complete Release suite passed twice during this task.  The persisted run
contains `509/509` passing tests and zero failures.  It explicitly contains
passing entries for:

- mixed-Voronoi conservation;
- obtuse split and winding invariance;
- omitted-policy legacy fallback;
- DCSweep mixed-Voronoi selection.

Persisted log SHA-256:
`302425d807bc0edd0107b1d4f7aca0e5658f2bc9c3278cbc3cd4892ea3e5d64a`.

## Review boundary

The evidence supports control-volume geometry as the dominant prior M2
discretization mismatch: mixed-Voronoi improves both off and on agreement on
M0/M2 without changing the postprocess-only IIC state or forward IV.

Two limitations remain explicit:

1. Sentaurus box-node volumes have not been exported directly.  The geometry
   conclusion is a controlled inference from same-grid output parity, not a
   byte-level internal-volume comparison.
2. No production template has been modified in this task.  Independent
   reviewers must assess the scientific evidence and the eventual minimal
   template/default-render patch before authorization.

## Evidence

- contract:
  `docs/validation/contracts/pn2d_node_volume_policy_default_acceptance_v1.json`;
- aggregate machine decision:
  `build-release/pn2d-node-volume-default-acceptance-v1-20260801/acceptance.json`;
- M0/M2 control reports under the same build root;
- forward report:
  `build-release/pn2d-node-volume-default-acceptance-v1-20260801/forward-iv/acceptance.json`;
- persisted CTest log:
  `build-release/pn2d-node-volume-default-acceptance-v1-20260801/release-ctest.log`.

## Independent review packet

The scientific reviewer should verify the same-grid input hashes, off/IIC
controls, on-curve and knee estimators, forward-IV guardrail, and the stated
inference boundary around unexported Sentaurus box volumes.

The code reviewer should inspect only an atomic PN2D BV template/default-render
patch, confirm that parser omission still means barycentric, require explicit
barycentric rollback coverage, and confirm that no general solver, forward-IV,
or legacy configuration default changes are bundled with it.
