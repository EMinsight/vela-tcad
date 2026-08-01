# PN2D BV M2 node-volume-policy frozen first-step comparison

Date: 2026-08-01

Status: complete; observation-only; production defaults unchanged.

Typed outcome:

```text
material_node_volume_policy_sensitivity
```

## Decision

Switching the frozen M2 equations from barycentric to mixed-Voronoi node
volumes materially changes the full first Newton update, but it does not change
the frozen SG/Laux integrated impact source and barely changes the dominant
joint-QFP carrier soft modes.  The measured sensitivity is dominated by
Poisson inconsistency: a state converged with barycentric control volumes is
not a stationary state of the mixed-Voronoi Poisson equation.

Therefore node-volume policy can alter self-consistent state formation, but
this experiment does not support it as the direct cause of the previously
isolated carrier-block near-null modes.  Nodal doping redistribution remains
unjustified.

## Frozen contract

The prospective contract was written before the simulation:

```text
docs/validation/contracts/pn2d_bv_m2_node_volume_policy_first_step_v1.json
SHA-256 352c9f0fc87027a30b8009ab6ab7438558e775dade81263e4d1e8495170072e6
```

It fixes:

- biases `-18`, `-19.5`, `-19.7`, and `-20 V`;
- Vela baseline and joint Sentaurus-QFP frozen states;
- barycentric and mixed-Voronoi policies;
- the same M2 mesh, nodal doping, SG/Laux settings, and physical parameters;
- no state advancement and no doping redistribution; and
- two byte-identical independent runs.

The predeclared classification is material when any decision metric changes
by at least `5%`.  The raw classification is material, but the physical cause
must be read from the block-resolved observations below.

## Joint Sentaurus-QFP result

The joint-QFP state is the discriminating carrier-soft-mode state.

| Bias (V) | full-step change | full-step cosine | carrier-step change | carrier-step cosine | carrier residual change | source change |
|---:|---:|---:|---:|---:|---:|---:|
| -18.0 | `18.23%` | `0.41981` | `1.34%` | `0.999972` | `0.00021%` | `0` |
| -19.5 | `11.99%` | `0.60041` | `1.80%` | `0.999994` | `0.01544%` | `0` |
| -19.7 | `12.73%` | `0.59876` | `1.74%` | `0.999994` | `0.00957%` | `0` |
| -20.0 | `14.50%` | `0.59392` | `1.65%` | `0.999995` | `0.00352%` | `0` |

The L2-equilibrated carrier-block condition-number ratio remains within
`1.000000009-1.000000020`.  At `-20 V`, the dominant-mode normalized terms
are also effectively unchanged:

| Policy | transport / sigma | avalanche diagonal / sigma | avalanche cross / sigma |
|---|---:|---:|---:|
| barycentric | `9.1821249` | `-3.2702188` | `-4.9119090` |
| mixed Voronoi | `9.1821265` | `-3.2702187` | `-4.9119107` |

Thus the transport-avalanche stiffness cancellation survives the policy
change almost exactly.  The carrier update changes slightly in magnitude but
not direction.  The larger full-step change comes from the Poisson block.

## Poisson-state consistency

For the barycentric Vela baseline state, initial Poisson residuals are about
`1.8e-8` to `2.2e-8`.  Re-evaluating exactly the same state with mixed-Voronoi
volumes raises them to about `25.98-26.13`.  Consequently:

- full first-step norm changes by `277-286x`;
- full-step direction cosine is approximately zero;
- carrier-only step changes by `3.3-4.7x`, although its absolute magnitude
  remains only about `0.0087-0.011 V`; and
- the initial carrier residual changes by less than `0.7%`.

This is expected for a frozen state solved under a different Poisson
control-volume discretization.  It proves material state-equation sensitivity,
not superiority of either volume policy.

## Frozen SG/Laux source

The integrated electron, hole, and combined impact sources are exactly equal
between policies for every bias and state.  The maximum relative combined
source change is `0`.

This separates the roles cleanly:

- node-volume policy changes Poisson and the subsequent coupled update path;
- it does not change SG/Laux read-only source evaluation on the same state;
- it does not remove or materially rotate the dominant carrier soft modes.

## Determinism and guards

All 16 policy/bias/state artifact groups are byte-identical across two
independent runs (`32` run cases).  The experiment changed only the opt-in
`mesh_geometry.node_volume_policy` field.  It did not alter:

- mesh nodes or triangles;
- nodal donors or acceptors;
- mobility doping basis;
- SG/Laux current, alpha, or source mapping;
- continuation settings;
- production template defaults; or
- acceptance thresholds.

## Interpretation

The prior `1.0-1.5` mixed/barycentric volume ratio on soft-mode support is a
real sensitivity, but its dominant immediate effect is the Poisson fixed-charge
balance.  In the carrier block that matters for the M2 knee discrepancy:

```text
condition number          essentially unchanged
dominant modal balance    essentially unchanged
carrier-step direction    cosine >= 0.99997
frozen impact source      exactly unchanged
```

Therefore the direct carrier-soft-mode root cause remains the signed
transport-avalanche modal cancellation.  Control volumes may still move the
self-consistent branch indirectly through Poisson, so the next experiment
must evaluate both avalanche-off and avalanche-on branches rather than only an
avalanche-on fit.

## Recommended next experiment

Do not redistribute doping.  If Sentaurus box volumes cannot be exported,
run an opt-in self-consistent mixed-Voronoi control with the following order:

1. avalanche-off at the four frozen knee biases, requiring that the existing
   same-grid golden agreement is not degraded;
2. IIC/postprocess-only, requiring state/hash equivalence to the corresponding
   avalanche-off policy branch;
3. SG/Laux-on only if the off/IIC controls pass; and
4. compare QFP, carrier density, modal balance, source, terminal current, and
   knee against both barycentric Vela and Sentaurus.

A mixed-Voronoi candidate should be rejected as a BV correction if it worsens
the already excellent avalanche-off agreement, even if it happens to improve
the avalanche-on knee.

## Reproduction

```powershell
python scripts\run_pn2d_bv_m2_node_volume_policy_first_step.py `
  --runner build-release\vela_example_runner.exe `
  --base-config build-release\pn2d-bv-m2-qfp-carrier-jacobian-verification-20260731\self_consistent_sg_laux_probe.json `
  --prior-root build-release\pn2d-bv-m2-qfp-carrier-jacobian-verification-20260731 `
  --contract docs\validation\contracts\pn2d_bv_m2_node_volume_policy_first_step_v1.json `
  --output-root build-release\pn2d-bv-m2-node-volume-policy-first-step-20260801 `
  --repeats 2
```

Generated outputs remain under `build-release/` and are not committed.
