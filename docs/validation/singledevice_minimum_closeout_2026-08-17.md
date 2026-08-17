# SingleDevice minimum closeout (2026-08-17)

## Outcome

The minimum closeout path is complete.  The original Sentaurus Save/Load
dependency graph, both self-consistent Id-Vg branches, the hybrid low-current
rule, derived device metrics, and three exact-bias field/KCL checks all pass.
No physics model or device equation was added.

## Low-current rule

The ordinary point gate remains a 10% relative-current limit.  Only when the
Sentaurus reference current is below `1e-13 A/um` may a point use the fallback:

- logarithmic error at most `0.2 dex`; and
- absolute error at most `1e-14 A/um`.

This rule is fail-closed on both guards.  The linear `Vg=-0.5 V` point uses the
fallback: `10.9097%`, `0.0501694 dex`, and `1.36067e-15 A/um`.  All other
points pass the ordinary relative rule.

## Save/Load workflow

The automated graph is:

```text
common equilibrium at Vg=-0.5 V, Vd=0 V -> Save
  -> load -> Vd 0...0.1 V -> linear Id-Vg (21 exact points)
  -> load -> Vd 0...1.1 V -> saturation Id-Vg (41 points, including 21 exact reference points)
```

The first four stages converged during the clean current-source build run.  A
branch-specific solver audit then retained the established 0.025 V quasi-Fermi
update limit for saturation and used 0.1 V for cold equilibrium and the linear
gate sweep.  Hash-checked resume reused only unchanged, passed stages and
reran the changed saturation stage.  All five final stage statuses are pass.

The distinction is numerical rather than physical: using 0.1 V on saturation
caused a hole-continuity line-search failure near `Vg=-0.3 V`; restoring the
proven 0.025 V setting completed all 41 points.  Mesh, material data, mobility,
SRH, BGN, density-gradient physics, contacts, and target biases were unchanged.

## Curve and derived metrics

| Metric | Linear, Vds=0.1 V | Saturation, Vds=1.1 V |
| --- | ---: | ---: |
| exact Sentaurus comparison points | 21 | 21 |
| median absolute log-current error | 0.00339653 dex | 0.00259518 dex |
| P95 absolute log-current error | 0.0329827 dex | 0.0316721 dex |
| Ion relative error | 0.654136% | 0.259824% |
| constant-current Vth delta | 1.96427 mV | 1.92168 mV |
| subthreshold-swing relative error | 0.168043% | 0.183335% |

DIBL absolute error is `4.25899e-05 V/V`.  Every acceptance check in the
reference manifest passes.

## Three-bias field and KCL audit

All three linear-branch checkpoints use exact Sentaurus TDR states on the
original 3584-node mesh; no field interpolation is used.

| Vg | psi median / P95 | electron density median / P95 | eQP median / P95 | KCL residual / drain |
| ---: | ---: | ---: | ---: | ---: |
| -0.5 V | 0.009852 / 0.022238 V | 0.021800 / 0.332307 dex | 0.000264 / 0.008321 V | 0.182605% |
| 0.31 V | 0.009852 / 0.022436 V | 0.017656 / 0.331806 dex | 0.000277 / 0.008040 V | 0.003686% |
| 2.2 V | 0.009852 / 0.028201 V | 0.017737 / 0.332572 dex | 0.000342 / 0.011029 V | 0.021408% |

KCL passes the hybrid rule `absolute residual <= 1e-14 A/um OR residual/drain
<= 1%` at every point.  Field aggregates are diagnostic: they document where
local tails remain wider than terminal quantities, but do not independently
authorize implementation of another model.

## Reproduction and evidence

Reusable scripts:

- `scripts/run_singledevice_workflow.py` materializes, executes, hashes, and
  safely resumes the five-stage graph;
- `scripts/close_singledevice_validation.py` applies the point and derived
  metric gates;
- `scripts/analyze_singledevice_closeout_fields.py` compares the three exact
  checkpoints and audits terminal KCL.

Compact durable results are in
`reference_tcad/singledevice_sentaurus2018/singledevice_validation_20260817.json`.
Large states, diagnostics, and comparison CSVs remain ignored under
`build/singledevice_closeout/workflow_current_run1`.
