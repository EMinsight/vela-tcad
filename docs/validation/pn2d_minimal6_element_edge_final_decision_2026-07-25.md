# PN2D Minimal6 element-edge final decision

Date: 2026-07-25

Status: complete; production default unchanged.

## Decision

Keep `element_edge_sg_gss_laux` as an explicit diagnostic operator. Do not
change the production avalanche-current default.

The fixed-state evidence establishes that the opt-in operator reproduces the
documented Sentaurus box-method avalanche support on the four-triangle
Minimal6 mesh. The self-consistent evidence does not establish parity:
electron QFP is the first failed dependency metric, followed by carrier
density, mobility, directed current, terminal current, and avalanche source.

## Dependency-chain result

| Stage | Best current evidence | Decision |
| --- | --- | --- |
| electrostatic potential | self-consistent maximum `2.36565e-11 V` | keep Poisson formula |
| imported QFP to density | maximum density replay error `4.42618e-6 dex` | keep Old-Slotboom/BGN statistics |
| self-consistent QFP | electron/hole medians `0.0514769/0.0591186 V` | typed `model_difference`; no formula patch |
| self-consistent density | electron/hole medians `0.864770/0.993144 dex` | downstream of QFP mismatch |
| fixed-state mobility/current | cell-vector medians `0.0337634/0.0318710 dex` | keep mobility and SG formulas |
| self-consistent directed current | electron/hole medians `0.823152/0.939636 dex` | unresolved upstream state/support difference |
| avalanche coefficient | fixed-state maximum `1.55e-9/2.56e-9 dex` | keep Van Overstraeten formula; use electric-field driver for this Sentaurus target |
| fixed-state source integral | electron/hole medians `0.0548177/0.0517363 dex` | element-edge diagnostic support validated |
| self-consistent source integral | full-lattice median `1.23508 dex`; `-16..-20 V` median `0.154333 dex` | improved but misses frozen target |

## Formula and implementation ledger

1. Poisson and electrostatic-field formulas are unchanged because the
   self-consistent electrostatic potential already passes by five orders of
   magnitude relative to the frozen tolerance.
2. Carrier statistics are unchanged because imported Sentaurus potential and
   QFP reproduce the Sentaurus carrier state to below `5e-6 dex`.
3. Masetti/high-field mobility and variable-intrinsic-density SG formulas are
   unchanged. Documented mobility parameters match after unit conversion;
   remaining native-element and element-to-edge differences are typed model
   and support differences.
4. The Van Overstraeten coefficient formula is unchanged. The fixed-state
   coefficient matches Sentaurus to numerical precision when driven by the
   element electric-field magnitude.
5. The new `element_edge_sg_gss_laux` current support and
   `element_vertex_box_measure` mapping remain opt-in. They resolve the
   fixed-state avalanche-source support discrepancy but do not repair the
   earlier self-consistent QFP/current mismatch.
6. The final implementation audit corrected one diagnostic-path defect:
   cell-local residual and Jacobian replay now honor
   `MobilityModelConfig.highFieldDrivingForce` instead of unconditionally
   using the QFP edge gradient. Regenerated 40-state metrics are unchanged at
   report precision.

## Geometry conclusion

Minimal6 consists of right triangles. The hypotenuse box partial volume is
zero, so perturbing that edge current does not change the GSS/Laux cell
vector. The equality between GSS/Laux and the two-positive-edge control is
therefore expected.

This is not a general three-edge deletion rule. On an acute scalene triangle,
all three partial volumes are positive, perturbing any one edge changes the
reconstructed vector, forward and reverse cell orientation recover the same
constant vector, and the weights sum to the exact cell area.

## Evidence and verification

- Remote Sentaurus regeneration: 40 exact `mirror/sketch x -1..-20 V`
  states, zero failures.
- Fixed-state verifier: 480 element-vertex rows, 40 states, 160 exact
  zero-diagonal partial volumes, source identity maximum relative error
  `4.10e-16`.
- Phase F verifier: `status: passed`, `deterministic_pair: true`,
  `outcome: model_difference`.
- A/B sweep comparison: all 40 paired state CSV SHA-256 hashes are identical.
- Full Release CTest: 475 of 475 passed.
- Final focused geometry/mobility test: 59 assertions in 4 test cases passed.
- ASCII source contract: passed.
- Independent scientific review: APPROVE, no blocking findings; its
  determinism, diagonal-support, general-mesh, and conclusion-scope findings
  were incorporated.
- Independent code review: APPROVE, no blocking findings; configuration
  propagation, residual/Jacobian consistency, test coverage, and unchanged
  production defaults were confirmed.

## Remaining scientific boundary

The result proves operator closure for the documented Minimal6 fixed-state
Sentaurus reconstruction, not general-device or self-consistent Sentaurus
parity. Further production work should begin at the imported-state
continuity residual and first QFP update on a larger non-right-triangle mesh.
It should not tune mobility, SG, or impact coefficients to absorb the
remaining QFP-state difference.
