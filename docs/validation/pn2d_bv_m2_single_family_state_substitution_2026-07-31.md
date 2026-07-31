# PN2D M2 single-family state substitution and first Newton update

## Technical summary

The predeclared four-bias experiment returns
`qfp_dominant__density_feedback_moves_qfp_away_from_sentaurus`. On the shared
M2 mesh, with SG/Laux and all physics settings unchanged, replacing only the
electron and hole quasi-Fermi potentials recovers 73.81% of the full
Sentaurus-state integrated-source error removal at -20 V and wins three of the
four bias comparisons. Replacing only `n/p` does not change the frozen SG/Laux
source under the current Masetti plus field-dependent mobility configuration.

The first coupled update does not preserve the golden QFP direction. A frozen
Sentaurus density substitution moves QFP slightly away from the Sentaurus
target at every bias; at -20 V its target projection is -0.0025435. Direct QFP
feedback produces a much larger negative projection. The first discrepancy is
therefore localized to the coupled residual/Jacobian response to QFP and its
state consistency, not to read-only SG/Laux evaluation on a fixed state.

## Fixed-source single-family substitutions

The fraction column is relative to the error removed by substituting the full
Sentaurus state. Negative values mean the one-family substitution makes the
integrated source less accurate than the Vela baseline.

| Bias (V) | Vela ratio | psi-only ratio / recovery | QFP-only ratio / recovery | n/p-only ratio / recovery | Full Sentaurus-state ratio |
|---:|---:|---:|---:|---:|---:|
| -18.0 | 0.936135 | 1.212593 / -1.9957 | 1.123844 / -0.7991 | 0.936135 / 0.0000 | 1.002481 |
| -19.5 | 0.881532 | 1.212453 / -0.5379 | 1.070904 / 0.4655 | 0.881532 / 0.0000 | 1.002381 |
| -19.7 | 0.865341 | 1.198302 / -0.2550 | 1.064073 / 0.5802 | 0.865341 / 0.0000 | 1.002404 |
| -20.0 | 0.825003 | 1.154008 / 0.2586 | 1.053517 / 0.7381 | 0.825003 / 0.0000 | 1.002371 |

QFP is not uniformly beneficial over the entire range: at -18 V, replacing it
alone overshoots the golden source. Its explanatory power grows through the
knee region and is dominant by -20 V. The result therefore supports a
bias-dependent coupled-state mechanism rather than a constant source scaling.

The density-neutral result is an algorithmic property of this configuration,
not evidence that the imported densities are identical. At -20 V the maximum
nodewise Vela/Sentaurus relative differences are about 62.0% for electrons and
50.4% for holes. The canonical SG edge flux is reconstructed from `psi`, QFP,
and intrinsic density. In this case, Masetti mobility uses doping and its
high-field correction uses the QFP-gradient drive, so explicit `n/p` does not
change the frozen current/source path. It still affects Poisson charge,
recombination, and other residual paths.

## First coupled Newton update

Two intervention types are intentionally separated:

- `psi` and QFP are independent Newton-state substitutions.
- `n/p` and QFP feedback substitutions are frozen residual-operator inputs
  evaluated with the same baseline production Jacobian. Density enters the
  Poisson charge, mobility, recombination, and avalanche paths; it is not a
  fourth independent Newton unknown.

| Bias (V) | QFP-only trial residual / initial | QFP-only combined target-distance ratio | Density-feedback QFP projection | QFP-feedback QFP projection |
|---:|---:|---:|---:|---:|
| -18.0 | 0.19657 | 1.61405 | -0.004072 | -1.32404 |
| -19.5 | 4.04726 | 1.86310 | -0.002975 | -1.13269 |
| -19.7 | 3.96492 | 1.85140 | -0.002807 | -1.14825 |
| -20.0 | 3.90581 | 1.82348 | -0.002543 | -1.18400 |

Starting from the Sentaurus QFP while retaining Vela `psi` increases the
combined distance from the full Sentaurus target at all four biases. In the
knee region (-19.5 to -20 V), the first production update also increases the
combined residual by 3.91-4.05 times. In the operator-feedback probe, the
Sentaurus QFP replacement produces a negative QFP target projection of roughly
-1.13 to -1.32, so the update opposes rather than preserves the imported QFP.

The density-only feedback experiment is also inconsistent with the baseline
Vela state: its initial combined residual is dominated by the Poisson block
(4.25 at -18 V, rising to 5.36 at -20 V), while the first update reduces it by
only about 2% and moves QFP 0.25-0.41% farther from Sentaurus. This does not make
density the frozen-source cause; it shows that independently replacing density
violates the thermodynamic/coupled state relation represented by Vela's three
unknowns.

## Contract and robustness

- Mesh: M2, 115 common physical nodes and 191 triangles.
- Bias lattice: exactly -18, -19.5, -19.7, and -20 V.
- Golden reference: Sentaurus on the identical mesh and input parameters.
- Vela operator: unchanged production element-edge SG/GSS/Laux source.
- Frozen source mode: `postprocess_only`; no continuity feedback or state
  advancement.
- Newton probe: first production update only; no continuation step is accepted.
- Two independent runs produced 120 byte-identical node, edge, triangle,
  element, process, and Newton-feedback artifacts.
- No physics model, default, continuation schedule, or acceptance threshold was
  changed.

## Conclusion and next localization step

The source-side observable first becomes QFP-dominant near the knee, while the
dynamic probe shows that Vela's coupled equations reject the golden QFP
direction. The strongest supported localization is the carrier-QFP residual
and its cross-coupling with Poisson/density consistency. It is not yet evidence
for changing SG/Laux or for a particular sign error.

The next read-only experiment should split electron and hole QFP, decompose the
first residual into transport, recombination, avalanche, and boundary terms,
and finite-difference the corresponding Poisson-QFP and carrier-QFP Jacobian
blocks on both the Vela baseline and the mixed Sentaurus-QFP state. Production
defaults remain unchanged until that operator-level derivative comparison is
complete.

## Reproduction

Run `scripts/run_pn2d_bv_m2_single_family_state_substitution.py` with the frozen
M2 inputs, then verify `result.json`, `source_substitution.csv`,
`newton_first_update.csv`, and `determinism.csv` using
`scripts/verify_pn2d_bv_m2_single_family_state_substitution.py`.
