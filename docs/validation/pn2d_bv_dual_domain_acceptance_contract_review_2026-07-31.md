# PN2D BV dual-domain prospective acceptance-contract review

Date: 2026-07-31

Reviewed contract:

`docs/validation/contracts/pn2d_bv_dual_domain_acceptance_v1.json`

Typed decision:

`bv_model_consistent_low_current_precision_floor_open`

## Purpose

This review freezes a prospective reporting and authorization contract that
separates:

1. BV-active impact-ionization model consistency; and
2. low-current nonlinear-solver precision.

It does not edit or replace the historical Task 7
`tradeoff_without_parity` score. It does not add a current floor to that
historical score, change a production default, alter the SG/Laux candidate, or
change any physical parameter.

## Independent-review method

The contract was authored from previously frozen gate definitions rather than
new thresholds fitted to the SG/Laux result:

- curve and knee thresholds come from
  `scripts/analyze_pn2d_avalanche_on_bv_parity.py`;
- fixed-state current/source thresholds come from
  `scripts/score_pn2d_task7_frozen_sg_candidate.py`;
- self-consistent state and determinism gates come from
  `scripts/analyze_pn2d_bv_task7_candidate.py`; and
- low-current classification evidence comes from
  `scripts/audit_pn2d_task7_low_current_nonmonotonicity.py`.

The new evaluator is independent of those scorers: it reads their sealed
outputs, verifies that the two bias domains are disjoint, applies the frozen
thresholds, and writes a new review artifact. It never rewrites an input
scorecard.

## Frozen domains

| Domain | Exact bias set | Purpose |
|---|---|---|
| BV model consistency | `-15, -16, -17, -18, -18.5, -19, -19.25, -19.5, -19.7, -19.8, -19.85, -19.9, -19.95, -20 V` | Curve, gain, knee, internal current/source, state feedback, closure, and determinism |
| Low-current solver precision | `-3, -4, -5, -6, -7 V` | Classify residual-floor behavior without using it as a BV-model gate |

The low-current classification is eligible only if all Vela/Sentaurus on/off
currents at the exact low-current biases remain below
`1e-15 A/um`. Consequently, a low-bias high-current runaway such as the sealed
triangle baseline cannot be waived as solver-floor behavior.

## BV-active replay result

All BV-model gates pass.

| Metric | Observed | Frozen threshold | Result |
|---|---:|---:|---|
| effective-curve median error | `0.0116262 dex` | `<= 0.05 dex` | pass |
| effective-curve P95 error | `0.0217208 dex` | `<= 0.10 dex` | pass |
| effective-curve maximum error | `0.0220198 dex` | `<= 0.15 dex` | pass |
| effective-gain median error | `0.0325488 dex` | `<= 0.05 dex` | pass |
| effective-gain maximum error | `0.0994154 dex` | `<= 0.10 dex` | pass |
| knee median error | `0.0170771 dex` | `<= 0.05 dex` | pass |
| knee maximum error | `0.0220198 dex` | `<= 0.10 dex` | pass |
| `V_break` error | `0.0210 V` | `<= 0.10 V` | pass |
| `V_slope` error | `0.0164329 V` | `<= 0.10 V` | pass |
| adjacent-slope RMSE | `0.0416135 dex/V` | `<= 0.20 dex/V` | pass |

The following previously sealed gate groups also pass:

- fixed-state integrated source, matching electron/hole current, active local
  source, vector direction, duplicate determinism, and the sign-only negative
  control;
- self-consistent QFP and density state improvement;
- WP7 process closure;
- duplicate exact-lattice determinism; and
- Task 6 feedback-state improvement and determinism.

## Low-current replay result

All classification gates pass:

- maximum exact-domain current is `3.1399296735e-17 A/um`, below the
  `1e-15 A/um` eligibility limit;
- reverse intervals are shared with avalanche-off/IIC;
- off/IIC states and duplicate avalanche-on states satisfy the sealed hash
  checks;
- raw avalanche source is monotonic;
- SG and residual terminal currents agree;
- there are no continuation retries;
- the contact QFP drop is at the precision floor;
- drift/diffusion cancellation exceeds `1e5`; and
- tightening Newton tolerances changes the interval pattern and terminates at
  `stall_residual_floor`.

Therefore the low-current behavior is classified as a solver-precision floor,
not an impact-ionization model inconsistency.

## Review checks

| Check | Result |
|---|---|
| BV and low-current exact bias domains are disjoint | pass |
| historical `tradeoff_without_parity` score is preserved | pass |
| raw low-current monotonicity is excluded from the BV gate | pass |
| high-current runaway cannot use the low-current waiver | pass |
| production-default change remains unauthorized | pass |

## Authorization boundary

This review authorizes only:

`opt_in_bv_model_validation_only`

It does not authorize:

- changing the production default;
- adding a minimum-field or current threshold;
- fitting Van Overstraeten coefficients or a source scale;
- changing the continuation schedule to hide a failed branch; or
- declaring the low-current nonlinear precision task complete.

A production-default change still requires a separate prospective deployment
contract and its required regression, mesh, forward-IV, and review gates.

## Verification

Focused regression:

```text
python -m unittest tests.regression.test_pn2d_bv_dual_domain_contract -v
```

Result: 4/4 passed.

Replay evaluator:

```text
python scripts/review_pn2d_bv_dual_domain_contract.py \
  --contract docs/validation/contracts/pn2d_bv_dual_domain_acceptance_v1.json \
  --curve-acceptance build-release/pn2d-task7-sg-laux-selfconsistent-curve-20260730/acceptance.json \
  --candidate-scorecard build-release/pn2d-task7-sg-laux-selfconsistent-scorecard-20260730/acceptance.json \
  --frozen-state-score build-release/pn2d-task7-frozen-sg-candidate-score-20260730/result.json \
  --low-current-audit build-release/pn2d-task7-low-current-nonmonotonicity-audit-20260730/acceptance.json \
  --output-root build-release/pn2d-bv-dual-domain-contract-replay-20260731
```

Generated replay artifacts:

- `build-release/pn2d-bv-dual-domain-contract-replay-20260731/acceptance.json`;
- `build-release/pn2d-bv-dual-domain-contract-replay-20260731/gates.csv`.
