# TransportModels DD deep-off fixed-state audit

Bias: `Vg=-1 V`, `Vd=1.1 V`. Sentaurus 2022 and Vela use the exact imported 3315-node topology.

## First result

| Quantity | Result |
|---|---:|
| Sentaurus exported SRH integral | -1.624301551e-15 A/um |
| Vela with the uncorrected imported config on the Sentaurus state | -5.232841500e-18 A/um |
| Vela SRH integral on its self-consistent state | -1.063505058e-15 A/um |
| Vela corrected fixed state with 1e16 cm^-3 internal reference | -1.083904487e-15 A/um |
| Fixed-state Vela/Sentaurus magnitude ratio | 0.00322159 |
| Self-consistent Vela/Sentaurus magnitude ratio | 0.654746 |
| Uncorrected imported config with Sentaurus n,p substituted | -5.138133518e-18 A/um |

## Diagnosis

The Sentaurus importer defect has been fixed: for a v1 `unit_scaling` deck it now writes the intended internal `1e16 cm^-3` value instead of the SI literal `1e22`. The corrected fixed-state path restores 66.7305% of the Sentaurus generation integral without changing the state.

Substituting the exact Sentaurus electron and hole densities changes the corrected Vela source by 1.8099%, so density reconstruction is not the leading cause. Generalized Fermi-SRH/ni-BGN semantics and source quadrature remain candidates for the fixed-state gap pending separate A/B tests.

Self-consistent diagnostic converged: `True`; drain current: `1.0664533466122335e-15` A/um.
Against the Sentaurus terminal current `1.634684064e-15` A/um, the corrected self-consistent relative error is `34.7609%`; current magnitude increases by `195.081x` over the uncorrected Vela run.

## Self-consistent state differences over silicon nodes

| Field | median | p95 | maximum |
|---|---:|---:|---:|
| psi_mV | 13.3656 | 111.634 | 113.787 |
| phin_mV | 1.22125e-11 | 130.589 | 173.862 |
| phip_mV | 1.45696e-09 | 207.657 | 261.059 |
| electron_density_dex | 0.0245836 | 2.23492 | 2.88013 |
| hole_density_dex | 0.644788 | 2.81649 | 5.38071 |

The fixed-state carrier-term and SG-edge CSV files are direct outputs from Vela production operators; no Python reimplementation is used for those terms.

Raw artifact directory: `D:\code-repo\vela-tcad\build-release\reference_tcad\transportmodels_sentaurus2022\reports\idvg_dd_deep_off_fixed_state_20260823`
