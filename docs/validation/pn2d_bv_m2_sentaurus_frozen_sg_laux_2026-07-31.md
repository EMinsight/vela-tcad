# PN2D M2 Sentaurus-state SG/Laux frozen replay

## Technical summary

The four-point discriminating experiment classifies the M2 knee discrepancy as `state_feedback_dominant`. When Sentaurus electrostatic potential, electron/hole quasi-Fermi potentials, and electron/hole densities are frozen and evaluated by Vela SG/Laux in `postprocess_only` mode, the mean total-source error is 0.001045 dex (about 0.24%). The self-consistent Vela mean error on the same biases is 0.057445 dex. Freezing the golden state therefore removes 0.056400 dex of source error on average.

This result localizes the dominant discrepancy to the coupled formation and feedback of the self-consistent state. It does not authorize a production-default change or identify one state equation as the initiating cause.

## Key findings

| Bias (V) | Sentaurus total source (A/um) | Vela self-consistent (A/um) | Vela on frozen Sentaurus state (A/um) | Self / Sentaurus | Frozen / Sentaurus | Error removed (dex) |
|---:|---:|---:|---:|---:|---:|---:|
| -18.0 | 2.434559e-16 | 2.279076e-16 | 2.440599e-16 | 0.936135 | 1.002481 | 0.027585 |
| -19.5 | 6.263344e-16 | 5.521336e-16 | 6.278259e-16 | 0.881532 | 1.002381 | 0.053729 |
| -19.7 | 7.629276e-16 | 6.601925e-16 | 7.647615e-16 | 0.865341 | 1.002404 | 0.061770 |
| -20.0 | 1.117532e-15 | 9.219670e-16 | 1.120182e-15 | 0.825003 | 1.002371 | 0.082516 |

The frozen total-source ratio remains nearly constant at 1.00237–1.00248. Electron and hole ratios independently remain within 1.00224–1.00264. In contrast, the self-consistent total-source deficit grows monotonically from 6.39% at -18 V to 17.50% at -20 V.

The combined frozen mobility-current-alpha-source chain therefore closes against Sentaurus. A large error in the SG/Laux frozen operator, van Overstraeten coefficient, or element-vertex source measure is incompatible with this result at the integrated level.

## Scope and state mapping

- Grid: common M2 mesh, 115 Vela physical nodes and 191 triangles.
- Biases: -18, -19.5, -19.7, and -20 V, avalanche-on branch.
- Imported state: `psi`, `phin`, `phip`, `n`, and `p`.
- Sentaurus reference: native total/electron/hole integrated avalanche source.
- Operator: van Overstraeten, quasi-Fermi-gradient drive, current-density generation, complete element-edge SG/GSS/Laux current vector, and element-vertex box-measure mapping.
- Coupling: forced to `postprocess_only`; no continuity solve and no advancement to another voltage.

Sentaurus exports 122 physical-node records because seven contact-support vertices duplicate the coordinates of existing physical contact nodes. The replay imports the common node IDs 0–114 only. Across all four biases, the maximum coordinate mismatch is 1.11e-16 um.

## Methodology

For each bias, the script writes the five Sentaurus state fields to Vela's fixed-state CSV, forces `impact_ionization.coupling_mode=postprocess_only`, and evaluates the production SG/Laux operator. It sums carrier-resolved `qG_contribution_A_per_m` over element-vertex records and converts A/m to A/um. The result is compared with both the native Sentaurus source and the self-consistent Vela `solver_used` source on the exact lattice.

The typed decision rule was declared in the runner:

- `state_feedback_dominant`: mean frozen error <= 0.02 dex and mean error reduction >= 0.02 dex.
- `mixed_state_and_operator`: mean reduction >= 0.02 dex but the frozen error remains > 0.02 dex.
- `operator_or_support_dominant`: the frozen mean is no more than 0.01 dex below the self-consistent mean.

The observed 0.001045 dex frozen mean and 0.056400 dex mean reduction satisfy the first outcome with substantial margin.

## Robustness checks

- All 115 imported nodes round-trip with zero relative error.
- All process records have `solver_coupled=0`.
- No electron or hole residual-feedback contribution is nonzero.
- Carrier source integral and `qG` close within 2.83e-15 relative or better.
- Two independent runs are byte-identical for all 20 node, edge, triangle, element, and process artifacts.
- No physics model, production default, continuation schedule, or acceptance threshold was modified.

## Limitations

This frozen replay proves the combined Vela operator chain is accurate on the Sentaurus state, but it does not order the causal loop among potential, quasi-Fermi potentials, and densities. It also cannot exclude small compensating local errors hidden by an accurate integrated source. The earlier spatial and carrier-resolved decomposition remains relevant for judging local fidelity.

## Recommended next step

Keep SG/Laux and the acceptance contract unchanged. Use the same four frozen states for one-family-at-a-time substitutions (`psi`, QFP, and `n/p`) and then inspect the first coupled Newton update and carrier-row residual. The goal is to determine which state family recovers most of the 0.082516 dex improvement at -20 V and whether the first coupled update moves it in the same direction as the final self-consistent deficit.

## Further questions

1. Is density replacement alone sufficient, or is Poisson-QFP cross-coupling required?
2. Does the first coupled update reproduce the bias-growing source deficit before continuation can influence the branch?
3. After integrated closure, do the local hotspot source and carrier partitions remain within the existing spatial contract?
