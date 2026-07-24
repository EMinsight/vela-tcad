# PN2D Minimal6 mobility velocity-unit fix validation

Date: 2026-07-24

Status: schema-valid, deterministic, 40-state fixed-state validation passed.

## Resolution

The `unit_scaling` mobility path now has an explicit velocity dimension.
Physical SI defaults are converted from m/s to the active internal cm/s unit:

| Carrier | Physical default (m/s) | Old internal value (cm/s) | Fixed internal value (cm/s) |
|---|---:|---:|---:|
| electron | 1.07e5 | 1.07e5 | 1.07e7 |
| hole | 8.37e4 | 8.37e4 | 8.37e6 |

Explicit JSON saturation velocities are routed through the same velocity API.
As with the other public `unit_scaling` inputs, their numeric values use TCAD
units, so saturation velocity is supplied in cm/s and remains numeric
internally. Legacy-SI JSON inputs continue to use m/s.

The Masetti low-field expression and the high-field limiter formula were not
changed.

## RED and GREEN evidence

Before the implementation, the physical-parity test failed because the
electron unit-scaled default was `107000` instead of `10700000`.

After the implementation:

- electron and hole legacy-SI versus unit-scaled high-field mobility agree at
  relative tolerance `1e-13`;
- explicit unit-scaled JSON values `1.23e7` and `7.89e6` remain unchanged in
  active cm/s units;
- the complete mobility suite passes 102 assertions in 19 cases;
- the complete scaling suite passes 49 assertions in 12 cases.

## Direct mobility closure

The direct C++ triangle local-edge mobility was compared with two diagnostic
branches that differ only by a factor of 100 in saturation velocity.

| Support | Carrier | Branch | N | Median abs error (dex) | P95 (dex) | Maximum (dex) |
|---|---|---|---:|---:|---:|---:|
| triangle local edge | electron | legacy velocity | 480 | 1.874802 | 1.963048 | 1.968603 |
| triangle local edge | electron | correct velocity | 480 | 9.643e-17 | 1.929e-16 | 1.929e-16 |
| triangle local edge | hole | legacy velocity | 480 | 1.789944 | 1.939787 | 1.949354 |
| triangle local edge | hole | correct velocity | 480 | 9.643e-17 | 1.929e-16 | 1.929e-16 |
| Sentaurus native element | electron | correct, cell-average doping | 160 | 0.052688 | 0.207159 | 0.312971 |
| Sentaurus native element | hole | correct, cell-average doping | 160 | 0.047814 | 0.121378 | 0.184062 |

The direct production operator therefore selects the correctly converted
branch to floating-point precision. The remaining approximately 0.05 dex
native-element median difference is a support/model residual, not a remaining
factor-of-100 velocity error.

## Forty-state current and avalanche-source closure

The exact contract is 40 states: mirror/sketch topologies, reverse biases
-1 through -20 V, with 6 nodes, 9 global edges, and 4 cells per state.

Reference current rows are now aligned by the unordered node pair instead of
exporter-specific `edge_id`. This removes an independent support-order
ambiguity.

| Metric | Before velocity fix | After velocity fix | Interpretation |
|---|---:|---:|---|
| electron inferred mobility median abs error (dex) | 1.926183 | 0.034809 | endpoint current proxy |
| hole inferred mobility median abs error (dex) | 1.880891 | 0.032083 | endpoint current proxy |
| electron SG current median abs error (dex) | 2.004884 | 0.397171 | node-pair aligned after fix |
| hole SG current median abs error (dex) | 2.081550 | 0.474486 | node-pair aligned after fix |
| electron current sign agreement | 100% | 100% | nonzero proxy edges |
| hole current sign agreement | 50% | 100% | also benefits from node-pair alignment |
| self-consistent integrated generation median abs error (dex) | 2.273491 | 0.347282 | improves by 1.926209 dex |

The self-consistent generation ratio relative to Sentaurus changes from a
median `0.0053273` to `0.4494879`. The repaired result is still low by a
median factor of approximately `2.225`, but the prior approximately 188-fold
underestimate is removed.

The baseline Vela-state generation error increases from `3.222279` to
`5.132780` dex because the baseline Vela QFP state already overpredicts the
Sentaurus source and the corrected mobility raises its current magnitude.
This does not contradict the fixed-state replacement result: it shows that
mobility scaling and state/QFP agreement are separate error sources.

## Determinism and evidence

Two independent self-consistent roots have byte-identical node, edge,
generation, summary, report, and manifest files:

- `build-release/pn2d-minimal6-self-consistent-replacement-velocityfix-20260724-a`
- `build-release/pn2d-minimal6-self-consistent-replacement-velocityfix-20260724-b`

Two independent direct mobility roots also have byte-identical local-edge,
native-element, summary, report, and manifest files:

- `build-release/pn2d-minimal6-mobility-velocityfix-validation-20260724-d`
- `build-release/pn2d-minimal6-mobility-velocityfix-validation-20260724-e`

Key SHA-256 values:

| Evidence file | SHA-256 |
|---|---|
| self-consistent edge samples | `e92ecb8c961a21e6153fbc6d44b94d8d60f880d963496c3760054b679f3eb7b4` |
| self-consistent generation samples | `08fc568514be1652f169c5f4f0340d1831acc8f99c5b67c16d2a2e307c7ea9ac` |
| self-consistent summary | `d69a48255ee775c5e67166f321a52bd2eee0fdbbfb493f3c7c00b224e5266cbc` |
| direct mobility local-edge samples | `2e3c508c0e069d3e44672e15aed44fc8452330ed6e0c3c4ae8b4d46435d388eb` |
| direct mobility native-element samples | `5d2ded369bcb51c11360bd48ae4e54245b08f35fe27eaf05739fdc77795a2d1e` |
| direct mobility summary | `c8ad7745328b0a6ca0fb12910c188854f835e6baecb3a933d21e7f834da5b7c6` |

## Decision

The saturation-velocity unit defect is fixed and independently verified.
No further factor-of-100 mobility correction is justified. Remaining work
should focus on the approximately 0.35 dex integrated-source residual, the
endpoint-current proxy versus native directed-edge support, and the remaining
native-element mobility interpolation/temperature dependence.
