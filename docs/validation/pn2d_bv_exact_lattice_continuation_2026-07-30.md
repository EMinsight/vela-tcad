# PN2D BV exact-lattice continuation qualification

Date: 2026-07-30

Outcome: `complete_exact_lattice_process_manifest`.

## Scope

This qualification addressed the sealed Vela avalanche-on failure at the
transition
`-19.692187499999644 -> -19.693749999999643 V`. It did not change a physics
model, production default, residual, Jacobian, damping rule, or continuation
algorithm.

The baseline configuration used `max_iter=40`. Its rejected transition ended
after 40 Newton iterations with residual `1.0690707411832896e-7` and a
monotonically decreasing late-iteration trace. The qualification control
changed only the explicit run-local Newton budget to `max_iter=80`, including
`solver.handoff.newton_max_iter=80`.

## Exact-lattice result

The lattice was read directly from the accepted Sentaurus manifest:

`0, -1, ..., -18, -18.5, -19, -19.25, -19.5, -19.7, -19.8, -19.85,
-19.9, -19.95, -20 V`.

All three Vela branches completed all 29 requested points:

| Branch | Exact points | Attempts | Rejected attempts | Max accepted Newton iterations |
|---|---:|---:|---:|---:|
| avalanche-off | 29/29 | 401 | 0 | 3 |
| IIC/postprocess-only | 29/29 | 401 | 0 | 3 |
| avalanche-on | 29/29 | 3,243 | 638 | 59 |

The avalanche-on branch reached `-20 V` without a nearest-bias substitution.
The six required knee states from `-19.7 V` through `-20 V` are all present.
The 29 avalanche-off and IIC state files have identical SHA-256 values at
every bias.

The run-local controls therefore classify the old blocker as a Newton
iteration-budget exhaustion, not a physical divergence. The repository
default remains 40 and the sealed failure regression remains valid.

## Process manifest

The manifest builder:

- rejects any missing or inexact bias row;
- verifies byte-identical off/IIC states before using the IIC fixed-state
  avalanche probe for the off observation branch;
- emits physical-node state/density, same-cell drive/mobility/current,
  reconstructed node alpha/generation, element-local-vertex qG, terminal
  current, and exact-at-target Newton attempt records;
- records canonical units and support/provenance declarations; and
- validates every artifact and normalized-output hash with
  `vela.pn2d_bv_process_run.v1`.

Final manifest:

`build-release/pn2d-vela-exact-lattice-maxiter80-on-20260730/manifest.json`

SHA-256:
`b882ece81a9cd1e7633e5685adbdd1a9ffde8b4adf2d14dea4fbc2286d6ddf6d`.

It contains 68,034 field records, 783 aggregate records, and 87 exact-target
Newton-attempt records.

## WP7 rerun

The paired analyzer was rerun under:

`build-release/pn2d-wp7-process-chain-after-continuation-20260730`.

The run now has complete Sentaurus and Vela state-through-terminal process
records on the same exact lattice. All 203 independently derived
source-reintegration and carrier/total terminal-current closure rows pass.

This initial rerun returned `insufficient_observation` because the accepted
Sentaurus manifest had no spatial `residual_jacobian` or `newton_update`
records. A subsequent fixed-transition Newton follow-up supplied those
stages for both simulators at all six knee targets.

The final WP7 report is now under
`build-release/pn2d-wp7-process-chain-newton-complete-20260730`. It has no
missing stage observations and accepts `state` for Sentaurus IIC versus
avalanche-on at adjacent `-19.7/-19.8 V` biases. Its typed outcome is
`density_qfp_feedback_cause`.

The Vela `-20 V` avalanche-on solver-used integrated source is
`4.719448196205013e-5 A/um`; its independently scattered/reintegrated value
closes within the frozen `1e-12` relative gate.

## Code and test coverage

- `scripts/run_pn2d_bv_exact_lattice_process.py`
- `scripts/build_pn2d_bv_exact_lattice_manifest.py`
- provenance-priority, requested-bias normalization, knee-only adjacency, and
  explicit missing-stage reporting in
  `scripts/analyze_pn2d_bv_process_chain.py`
- focused regression coverage in
  `tests/regression/test_pn2d_bv_exact_lattice_process.py` and
  `tests/regression/test_pn2d_bv_process_chain.py`
- complete Release CTest: 503/503 passed.
