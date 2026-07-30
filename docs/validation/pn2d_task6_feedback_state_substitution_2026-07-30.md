# PN2D Task 6 density/QFP feedback-state substitution

Date: 2026-07-30

Outcome: `continuation_only_cause`.

Localization: `coupled_poisson_qfp_cross_block_reversal`.

## Controlled intervention

The diagnostic starts from the converged Vela avalanche-on exact state and
uses the exact Sentaurus avalanche-on state as the replacement source at
`-19.7 V` and `-19.8 V`. Every variant uses:

- one production baseline Jacobian;
- identical scaling, constraints, contact rows, and update caps;
- frozen residual-operator substitutions only; and
- no production configuration or default change.

The matrix contains baseline, electron/hole/combined density-only,
electron/hole/combined QFP-only, and combined density+QFP interventions.
Contact values remain on the Vela baseline. Both the full coupled Newton
update and the carrier-only block update are evaluated from the same residual.

## Cross-bias result

The isolated density and QFP interventions do not pass the coupled-update
causal gate at either bias.

| Variant/metric | `-19.7 V` | `-19.8 V` |
|---|---:|---:|
| density-only combined QFP-error improvement | `2.000%` | `2.330%` |
| hole-density-only combined improvement | `9.208%` | `9.191%` |
| hole-density-only electron-carrier change | `-14.842%` | `-14.490%` |
| hole-density-only hole-carrier improvement | `12.739%` | `12.770%` |
| QFP-only full-coupled QFP-error change | `-7.023%` | `-7.285%` |
| QFP-only full-coupled update-direction cosine | `-0.1994` | `-0.2132` |
| QFP-only carrier-only QFP-error improvement | `13.128%` | `12.942%` |
| QFP-only carrier-only update-direction cosine | `0.6412` | `0.6367` |

Thus QFP substitution moves both carriers toward the Sentaurus branch when
Poisson is frozen, but the complete coupled solve reverses that direction.
The same reversal occurs at both adjacent biases. Density-only substitutions
also show carrier antagonism: improving one carrier worsens the other.

This evidence rejects a direct density/QFP state correction and localizes the
remaining difference to the coupled Poisson-QFP/continuation path. Under the
companion plan, `continuation_only_cause` is not one of the three Task 6
outcomes that may authorize a physical/operator correction.

## Integrity gates

- residual decomposition maximum absolute closure error:
  `4.547473508864641e-13`, below `1e-12`;
- six contact nodes per bias: zero byte-level residual-row mismatches;
- duplicate output hashes:
  - `-19.7 V`:
    `453de2a850890f5da6b66591c14be914154aeed5079b1a7346d7dc19e3bba9bb`;
  - `-19.8 V`:
    `c6e346a681ce3a20a4f934431d62dbdcd1e25f4e465fe4164090ca48021ac923`;
- the pre-existing production `-19.7 V` Newton probe remained byte-identical
  at SHA-256
  `9c132f5c8f53f5c4714fddec2b704bb5781b78b4d6a25d09489d960d9b5ec8a3`.

The machine-readable scorecard is:

`build-release/pn2d-task6-feedback-substitution-final-scorecard-20260730/acceptance.json`

Its SHA-256 is:

`19a66d504adab81efa9ad521ad45ccb3bd0083cee2b74b4e1c712aad3668f60b`.

The final diagnostic executable SHA-256 is:

`12ffc3facd0cc320f7f855dc38e7516b5237edc6e49c9b1779554da33834c6c3`.

## Decision

Task 8 remains prohibited. The next authorized scientific slice is Task 7's
predeclared continuation-schedule branch-invariance control with unchanged
physics. It must still satisfy the complete curve/internal-causality gates
before any implementation task can proceed.
