# PN2D BV M2 hotspot-edge transport Jacobian verification

Date: 2026-08-01

## Outcome

The read-only hotspot-edge audit passed with typed outcome:

`transport_edge_decomposition_verified`

The first material change between the Vela baseline and the frozen joint
Sentaurus-QFP state is the carrier-population part of the SG transport
derivative.  Mobility response, continuity-row scaling, and contact-row
replacement are not the first source of the M2 knee-state difference.

No SG/Laux parameter, physical model, production default, continuation
schedule, or acceptance threshold was changed.

## Frozen contract

- Mesh: shared M2 mesh and doping.
- Biases: `-18`, `-19.5`, `-19.7`, and `-20 V`.
- States: Vela baseline and joint Sentaurus electron/hole QFP replacement.
- Hotspot node: the prior interior transport-residual peak for each bias and
  carrier.
- Hotspot edge: the incident edge with maximum absolute baseline-to-joint-QFP
  flux change; ties select the smallest edge ID.
- Derivative columns: both same-carrier QFP endpoints.
- Derivative rows: both carrier-continuity endpoint rows.
- Finite-difference step: `1e-7 V`, central symmetric.
- Repeats: two complete independent runs; all 16 raw CSV files must be
  byte-identical.

## Operator decomposition

For the electron edge flux, the frozen-mobility production form is

\[
F_n = C_n\left[B(-\eta_n)n_0-B(\eta_n)n_1\right],
\qquad
C_n=\mu_n V_T\frac{A_e}{h_e},
\]

with the corresponding sign-reversed hole population form.  For a
same-carrier QFP column:

\[
\frac{\partial B(\pm\eta)}{\partial QFP}=0,
\]

because \(\eta\) depends on electrostatic potential and the effective-ni ratio,
not on QFP.  The production derivative therefore contains the exponential
carrier-population derivative while holding high-field mobility fixed.

The diagnostic-only live-mobility counterfactual adds

\[
\left.\frac{\partial F}{\partial QFP}\right|_{\mu}
+ \frac{F}{\mu}\frac{\partial\mu}{\partial QFP}.
\]

This counterfactual does not change the production contract, whose
`mobility.jacobian_field_derivatives` remains `false`.

## Numerical verification

| Check | Maximum error | Result |
|---|---:|---|
| Production analytic vs frozen-mobility FD, normalized by dominant derivative on the same edge | `5.12344e-10` | pass |
| Mobility product `(F/mu) dmu/dQFP` closure | `2.11374e-16` relative | pass |
| Live total derivative vs production plus mobility response | `1.01736e-7` edge-scaled | pass |
| Bernoulli/GSS QFP derivative | exactly `0` | pass |
| Row-scaling multiplication closure | exactly `0` | pass |
| Contact-eliminated physical edge derivative on constrained rows | exactly `0` | pass |
| Contact identity-row error | exactly `0` | pass |
| Independent-run determinism | all 16 raw CSVs byte-identical | pass |

Some upwind-suppressed endpoint derivatives are tens of orders of magnitude
smaller than the dominant derivative on the same edge.  Their elementwise
relative FD error can therefore be one even when the absolute error is
irrelevant.  The raw elementwise values are retained, while the formal check
uses the dominant same-edge Jacobian scale.  This is a numerical-resolution
classification, not a relaxed physics threshold.

## Knee-region findings

The following ratios compare the frozen joint Sentaurus-QFP state with the
Vela baseline on the selected hotspot edge.

| Bias | Carrier | Edge | QFP-drive ratio | Mobility ratio | Dominant transport-derivative ratio |
|---:|---|---:|---:|---:|---:|
| -19.5 V | electron | 256 | `0.7271` | `1.1899` | `0.08115` |
| -19.5 V | hole | 122 | `0.7329` | `1.1316` | `0.06818` |
| -19.7 V | electron | 256 | `0.7678` | `1.1691` | `0.08404` |
| -19.7 V | hole | 122 | `0.7684` | `1.1240` | `0.06966` |
| -20.0 V | electron | 257 | `0.8532` | `1.1008` | `0.17202` |
| -20.0 V | hole | 122 | `0.8100` | `1.1122` | `0.07335` |

Across `-19.5` to `-20 V`:

- The Bernoulli/GSS coefficients are exactly unchanged because psi and
  effective ni are held fixed.
- QFP drive decreases by about `14.7%` to `27.3%`.
- High-field mobility increases by about `10.1%` to `19.0%`.
- The mobility-response derivative is only about `3.75%` to `7.00%` of the
  dominant transport derivative.
- The dominant production transport derivative falls to `6.82%` to `17.20%`
  of baseline, a reduction of roughly `6x` to `14.7x`.

The large derivative reduction must therefore come from the QFP-controlled
carrier-population factor, not from Bernoulli/GSS, mobility, or QFP-drive
sign handling.

At `-18 V`, both states are in the low-current/pre-knee regime.  QFP gradients
are near the numerical floor and their ratios are not used to infer high-field
physics.  The dominant transport derivatives remain equal within
`5.2e-12` relative, consistent with the earlier low-current precision finding.

## Row scaling and contact handling

In the knee region, hotspot row weights span roughly `3.8e9` to `5.0e15`, but
their joint-QFP/baseline ratios remain between `0.999998` and `1.000899`.
Thus row scaling strongly equilibrates the carrier system but does not create
the state-to-state change.

All constrained edge-row records have zero physical edge contribution after
contact replacement and the correct unit identity entry.  Fourteen
unconstrained neighbor rows retain derivatives with respect to constrained
contact columns, confirming the implemented semantics: replace constrained
rows only; do not eliminate constrained columns.  No sign or scale defect was
found in this path.

## Interpretation

This experiment narrows the first self-consistent discrepancy to how the
carrier block responds to QFP-controlled exponential populations.  It does
not show that the local SG derivative is wrong: the production derivative
passes finite differences.  The remaining question is how those locally
correct, strongly anisotropic derivatives combine in the complete carrier
linear solve and then couple back to Poisson and avalanche.

## Next step

Keep SG/Laux, mobility settings, row scaling, contact handling, and acceptance
thresholds unchanged.  Decompose the complete electron/hole carrier-block
linear solve on the same hotspot support:

1. identify the dominant QFP columns and left/right singular directions;
2. compare diagonal dominance and variable scaling between baseline and joint
   QFP states;
3. quantify electron-hole coupling through recombination and common avalanche
   rows;
4. project the carrier-only Newton update onto the local population-derivative
   directions.

Do not propose a production correction until that solve-level decomposition
identifies a wrong sign, missing derivative, or inconsistent variable scale.

## Evidence

- `build-release/pn2d-bv-m2-transport-edge-jacobian-verification-20260801/result.json`
- `build-release/pn2d-bv-m2-transport-edge-jacobian-verification-20260801/hotspot_selection.csv`
- `build-release/pn2d-bv-m2-transport-edge-jacobian-verification-20260801/hotspot_decomposition.csv`
- `build-release/pn2d-bv-m2-transport-edge-jacobian-verification-20260801/hotspot_state_change.csv`
- `build-release/pn2d-bv-m2-transport-edge-jacobian-verification-20260801/contact_row_audit.csv`
- `build-release/pn2d-bv-m2-transport-edge-jacobian-verification-20260801/determinism.csv`
