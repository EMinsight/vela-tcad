# PN2D BV M2 mixed-Voronoi self-consistent control — 2026-08-01

## Decision

The independently contracted execution order completed without opening a
later stage early:

1. `avalanche_off` passed;
2. `iic_postprocess` passed and was exactly state-equivalent to off;
3. only then did `avalanche_on` run, and it also completed.

The typed outcome is `completed_all_stages`.  This is observation-only evidence
and does not change the production node-volume or avalanche defaults.

## Frozen inputs

- node-volume candidate: `mesh_geometry.node_volume_policy=mixed_voronoi`;
- M2 mesh and nodal doping: unchanged from the paired Sentaurus/Vela M2 input;
- avalanche current/source candidate: complete SG/GSS-Laux edge current with
  element-vertex box source mapping;
- exact lattice: 29 biases from 0 V through -20 V, including the declared knee
  points;
- independent repetitions: two per branch;
- contract SHA-256: `1e85163812df6951ccbda36a60c7f24030a41b971e7bac974788989f01e9213e`;
- generated mixed config SHA-256:
  `1f373e1b944839044c1f019f3f189a542c9cff7024cac2beb22e0a7e73cae21b`.

## Prospective gate results

| Stage | Exact lattice | Run A/B IV | Run A/B states | Maximum continuity closure ratio | Gate |
| --- | ---: | --- | --- | ---: | --- |
| avalanche-off | 29/29 in both runs | identical SHA-256 | 29/29 identical hashes | 9.138e-5 | pass |
| IIC postprocess | 29/29 in both runs | identical SHA-256 | 29/29 identical hashes | 9.138e-5 | pass |
| SG/Laux-on | 29/29 in both runs | identical SHA-256 | 29/29 identical hashes | 1.917e-6 | pass |

The IIC run-A IV SHA-256 is exactly equal to mixed-off run A, and every one of
its 29 state hashes is equal to the corresponding off state.  The collision
ionization calculation therefore remained postprocess-only and did not feed
back into the carrier equations.

The mixed-off comparison against the Sentaurus native total terminal current
uses all 28 nonzero exact-lattice biases.  Its log-current RMSE is
`3.868e-5 dex`, and its maximum error is `9.223e-5 dex`, far inside the frozen
`0.01/0.015 dex` gates.

## Avalanche-on observation

These metrics were computed only after both controls passed and therefore did
not participate in the decision to enter SG/Laux-on.

| Metric | mixed-Voronoi Vela vs Sentaurus |
| --- | ---: |
| all 28 nonzero points, log-current RMSE | 0.001923 dex |
| all 28 nonzero points, maximum log-current error | 0.004647 dex |
| 11 knee points, log-current RMSE | 0.003051 dex |
| 11 knee points, maximum log-current error | 0.004647 dex |
| 11 knee points, maximum avalanche-gain log error | 0.004686 dex |
| fitted V_break | -19.390 V vs -19.391 V |
| absolute V_break error | 0.001 V |
| V_slope crossing | absent in both simulators |

Relative to the previously sealed barycentric M2 SG/Laux result, the
mixed-Voronoi on curve changes materially (maximum `0.08347 dex`).  The change
is toward the Sentaurus golden curve.  The earlier barycentric M2 V_break error
was about `0.014 V`; mixed-Voronoi reduces it to `0.001 V` under the same
estimator.  The off branch also improves rather than being sacrificed: its
all-point RMSE changes from `0.01008 dex` to `3.868e-5 dex`.

## Interpretation and boundary

The frozen first-step result showed that node-volume policy primarily moves the
Poisson stationary state while leaving the direct frozen SG/Laux source and
carrier soft-mode direction almost unchanged.  The self-consistent control now
shows the expected downstream consequence: after the Poisson state is rebuilt
with mixed-Voronoi volumes, both avalanche-off and avalanche-on move toward the
Sentaurus result.  The evidence therefore points to control-volume geometry as
the dominant remaining M2 discretization mismatch, rather than an IIC feedback
leak or a new SG/Laux current/source defect.

This result does not authorize a production default change.  A default-policy
proposal must separately bind M0, forward IV, boundary/obtuse-triangle
mixed-Voronoi behavior, the full Release test suite, and legacy configuration
fallback.

## Evidence

- contract: `docs/validation/contracts/pn2d_bv_m2_mixed_voronoi_self_consistent_control_v1.json`;
- machine report:
  `build-release/pn2d-bv-m2-mixed-voronoi-self-consistent-control-20260801/gate_report.json`;
- gate runner:
  `scripts/run_pn2d_bv_m2_mixed_voronoi_self_consistent_control.py`.
