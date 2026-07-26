# PN2D general-Tri3 fixed-state avalanche audit

Date: 2026-07-26

Plan:
`docs/superpowers/plans/2026-07-25-pn2d-avalanche-current-support-qfp-bv-followup.md`

Scope: Tasks 2-5 only. No production default or production physics formula
was changed.

## Outcome

The upstream imported-state transformations pass on the coarse and acute
skewed scientific-oracle meshes. The general current and source gates do not
pass:

- native general-Tri3 element alpha is unavailable and remains typed
  `insufficient_native_observation`;
- the documented Sentaurus output is vertex current density. The attempted
  vector `/Element` interpretation is undocumented and remains typed
  `insufficient_native_observation_undocumented_element_vector`;
- matching-support cell-current medians pass the `0.05 dex` gate, but P95,
  carrier-sign, terminal-current, and KCL gates fail;
- the coarse fixed-state source closes, but the acute skewed source does not;
  and
- the constrained obtuse mesh exposes nonconservative Vela negative-weight
  truncation and is diagnostic-only.

The typed Tasks 3-5 result is:

`support_limited_current_difference` plus
`insufficient_native_observation`.

The Van Overstraeten coefficient formula and the QFP-gradient driver remain
unchanged. Task 6 may use the reproducible obtuse-area defect as a RED case,
but it must keep the element-edge implementation opt-in.

## Deterministic oracle coverage

Three meshes were generated at exact biases `-1`, `-10`, and `-20 V`, with
eight declared Sentaurus variants and two independent roots:

| Mesh | Scientific role | Elements per state | A/B result |
|---|---|---:|---|
| coarse7x3 | device physics oracle | 32 | exact |
| skewed_tri3 | device physics oracle, acute scalene | 12 | exact |
| skewed_tri3_constrained | diagnostic only, obtuse | 8 | exact |

This is `3 meshes x 8 variants x 3 biases x 2 roots = 144` exact generated
states. The raw roots passed the independent Task 2 verifier. Parsed physics
outputs for each A/B pair have identical SHA-256 hashes.

## Task 3 - imported-state upstream replay

### Density and P1 vectors

| Mesh | max electron density error (dex) | max hole density error (dex) | max E relative error | max electron-QFP-gradient relative error | max hole-QFP-gradient relative error |
|---|---:|---:|---:|---:|---:|
| coarse7x3 | 4.70229e-6 | 4.72900e-6 | 7.62594e-15 | 1.66936e-13 | 1.06412e-15 |
| skewed_tri3 | 4.65443e-6 | 4.69293e-6 | 7.98131e-13 | 4.71706e-14 | 1.69898e-15 |
| constrained obtuse | 4.73051e-6 | 4.75249e-6 | 1.69198e-13 | 2.67437e-14 | 5.65988e-16 |

All listed values pass the fixed-state density and vector gates. The values
prove the Vela imported-state statistics and P1 gradient calculations on
these states; they do not prove self-consistent QFP parity.

### Mobility support

| Mesh | electron cell mobility median/P95/max (dex) | hole cell mobility median/P95/max (dex) |
|---|---:|---:|
| coarse7x3 | 1.31e-12 / 0.0692 / 0.2289 | 1.22e-12 / 0.0680 / 0.1506 |
| skewed_tri3 | 0.0144 / 0.1233 / 0.2460 | 0.0136 / 0.0834 / 0.1562 |
| constrained obtuse | 0.0363 / 0.3161 / 0.4169 | 0.0254 / 0.1853 / 0.2328 |

Endpoint-averaged edge mobility has substantially larger tails than native
element mobility. It is therefore a support difference, not a constant
mobility coefficient error.

Native general-Tri3 element alpha was not exported by a documented supported
Sentaurus observation. Task 3 alpha is therefore not claimed as a direct
native-element pass. The coarse Task 5 replacement provides indirect source
closure evidence only.

## Task 4 - current support and cell vectors

### Matching-support GSS/Laux cell-vector magnitude

| Mesh | electron median/P95 (dex) | hole median/P95 (dex) |
|---|---:|---:|
| coarse7x3 | 0.00585 / 1.27481 | 0.00285 / 0.92365 |
| skewed_tri3 | 0.03674 / 0.84526 | 0.01599 / 0.61346 |
| constrained obtuse | 0.33447 / 0.66846 | 0.14593 / 0.55743 |

The medians pass on the two scientific meshes. The P95 and 100% sign gates do
not pass. Near-zero and support-sensitive tails are retained rather than
converted into a fitted scale.

### Geometry, terminal current, and KCL

| Mesh | Sentaurus area max relative error | Vela truncated area max relative error | Sentaurus/Vela total-terminal max relative error | Sentaurus/Vela internal-total KCL max relative residual |
|---|---:|---:|---:|---:|
| coarse7x3 | 0 | 0 | 82.08 / 0.0458 | 41.47 / 0.372 |
| skewed_tri3 | 3.85e-16 | 4.34e-16 | 83.78 / 84.28 | 58.80 / 108.86 |
| constrained obtuse | 3.08e-16 | 3.2037 | 0.0522 / 6.19 | 0.0522 / 7.08 |

The element-local mobility plus ReadCoefficient replay is not the native
global Sentaurus continuity operator. In particular, element mobility can
differ on the two sides of an internal edge, so independently reconstructed
element-local fluxes need not cancel globally. The terminal and KCL failures
must not be used to fit an approximately `1.4258e6` conversion factor.

The constrained obtuse result identifies a separate Vela geometry defect:
clipping negative circumcentric coefficients destroys exact area
conservation. This mesh is diagnostic-only because it was deliberately
constrained to expose that policy.

## Task 5 - source support and driver selection

ReadMeasure integration and carrier-split CurrentPlot source integrals close
to:

| Mesh | maximum relative error |
|---|---:|
| coarse7x3 | 2.57318e-15 |
| skewed_tri3 | 3.45547e-15 |

The following table reports active-range median source-integral error for the
matching QFP-gradient alpha. `Sent box` replaces the current vector with the
Sentaurus box-operator reconstruction. `Vela` uses the Vela-recomputed
element-local SG current.

| Mesh/method | electron Sent box (dex) | hole Sent box (dex) | electron Vela (dex) | hole Vela (dex) |
|---|---:|---:|---:|---:|
| coarse GSS/Laux | 2.08e-6 | 0.000351 | 0.00681 | 0.00607 |
| coarse Charon HCurl | 0.000168 | 0.000494 | 0.06343 | 0.06110 |
| coarse Genius LS | 0.000126 | 0.000458 | 0.04944 | 0.04749 |
| coarse active-edge | 2.08e-6 | 0.000351 | 0.00681 | 0.00607 |
| skewed GSS/Laux | 0.46872 | 0.36769 | 1.31010 | 0.96076 |
| skewed Charon HCurl | 0.36347 | 0.27886 | 1.15358 | 0.81651 |
| skewed Genius LS | 0.46150 | 0.36498 | 1.28677 | 0.94140 |
| skewed active-edge | 0.46872 | 0.36770 | 1.31007 | 0.96071 |

For coarse GSS/Laux, the matching reconstructed-Sentaurus maximum errors are
`3.48e-6 dex` electron and `6.99e-4 dex` hole, and the Vela-current maximum
errors are `0.00940 dex` and `0.00902 dex`. The coarse source gates pass.

For skewed, the best reconstructed-Sentaurus medians are Charon HCurl:
`0.36347 dex` electron and `0.27886 dex` hole. Its maxima are `0.40740 dex`
and `0.32446 dex`; the general source gates fail. Changing the vector
reconstruction does not close the skewed discrepancy.

The contact-fallback and global-QFP candidates are identical in the active
integrals for these cases. A forced global-electric-field driver is worse.
This supports the documented global QFP-gradient default and does not support
changing alpha to an electric-field default.

Node alpha averages are consistently worse than the Vela matching-driver
alpha plus reconstructed box current. A vertex alpha average is not a native
element-alpha substitute.

## Formula decision ledger after Tasks 2-5

| Item | Decision | Evidence |
|---|---|---|
| P1 electric field | validated_unchanged | max relative error <= 7.99e-13 |
| electron/hole QFP gradient | validated_unchanged | max active relative error <= 1.67e-13 |
| Old-Slotboom/BGN density | validated_unchanged | max error <= 4.76e-6 dex |
| low/high-field mobility | diagnostic_only | medians small; element/edge tails remain |
| SG edge current | insufficient_data | no documented native directed-edge observation |
| current-vector reconstruction | diagnostic_only | median passes, P95/global closure fail |
| Van Overstraeten alpha | validated_unchanged | coarse replacement closure and driver controls |
| element-vertex measure | validated_unchanged | ReadMeasure/CurrentPlot <= 3.46e-15 relative |
| source mapping | diagnostic_only | coarse passes, skewed fails |
| obtuse negative-weight truncation | diagnostic_only | reproducible area nonconservation |
| production default | validated_unchanged | later parity gates not met |

## Independent verification and next gate

`scripts/verify_pn2d_general_tri3_fixed_state_pipeline.py` does not import the
generator or analyzer computations. It independently checks schemas, status,
exact biases, source hashes, output hashes, fixed-state gates, scientific
roles, and A/B output identity. It passes for:

- coarse7x3 A/B;
- skewed_tri3 A/B; and
- constrained-obtuse A/B as diagnostic-only.

Task 6 should begin with the constrained-obtuse area-conservation RED test and
then determine whether a nonnegative conservative support exists without
regressing acute/right triangles. It must not reinterpret reconstructed
Sentaurus edge current as native, must not fit a geometry/current scale, and
must keep the candidate opt-in.
