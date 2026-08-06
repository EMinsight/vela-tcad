# BVmethods NMOS avalanche source-mapping closure (2026-08-05)

## Scope

The accepted Fermi-Dirac SG transport state and the
`transport_cell_vector` high-field mobility discretization were held fixed.
Only the avalanche-generation interpolation, cell integration, and edge/node
source support were varied.

## 6.4 V source ledger

The Sentaurus node fields were integrated over the exact semiconductor P1
measure, `sum_cell(A/3)`, and compared with the Vela edge diagnostic after
mapping both endpoint contributions back to nodes.

| Quantity | Sentaurus P1 | Vela legacy | Vela / Sentaurus |
|---|---:|---:|---:|
| electron source (m^-1 s^-1) | 4.39356486458694e19 | 4.31422003110557e19 | 0.981941 |
| hole source (m^-1 s^-1) | 1.15408022494721e19 | 1.15792613281322e19 | 1.003332 |
| total source (m^-1 s^-1) | 5.54764508953415e19 | 5.47214616391880e19 | 0.986391 |

The corresponding current equivalents are 8.888307336e-6 A/um for the
Sentaurus exported nodal field and 8.767344722e-6 A/um for Vela. The raw
Sentaurus plot gives the more authoritative internal
`IntegrSemiconductor AvalancheGeneration=5.55281372910591e25`, or
8.896588410e-6 A/um. Thus the TDR P1 reconstruction is within 0.0931% of the
Sentaurus internal integration, and Vela is 98.5473% of the internal integral
at 6.4 V. The source-integral discrepancy is therefore 1.45%, not the
previously reported order-8% gap.

The previously used 9.5354995e-6 A/um reference is the integrated
`Band2BandGeneration` result. It must not be labelled or used as an avalanche
generation reference.

## Mapping controls

Two opt-in diagnostic controls were added without changing the default:

- `source_volume_policy=genius_conservative`: normalize the three truncated
  edge-box pieces in each transport triangle to the exact cell area and exclude
  insulator cells.
- `source_mapping_mode=nodal_eparallel_p1`: reconstruct nodal field/current,
  evaluate carrier-specific Eparallel/alpha at nodes, and use exact P1 node
  measures. This mode is post-processing only.

At 6.4 V the total-source ratios to the Sentaurus nodal P1 integral are:

| Vela mode | Vela / Sentaurus |
|---|---:|
| legacy `genius_truncated` plus edge/node mapping | 0.986391 |
| `genius_conservative` | 0.968380 |
| `nodal_eparallel_p1` plus conservative support | 0.957269 |

The legacy mapping is retained for the BVmethods target because it is the
closest of the tested mappings. The alternatives remain diagnostics and do not
change global or template defaults.

## Fixed-mapping high-voltage continuation

With the legacy source mapping frozen, the post-processed current-density
criterion gives:

| drain bias (V) | drain current (A/um) | avalanche current (A/um) | Iava / Id |
|---:|---:|---:|---:|
| 6.8 | 1.399271523e-5 | 1.368804531e-5 | 0.978227 |
| 6.9 | 1.529831934e-5 | 1.523881308e-5 | 0.996110 |
| 7.0 | 1.670013713e-5 | 1.693706847e-5 | 1.014187 |

Linear interpolation of `Iava - Id` locates the Vela crossing at
6.920073790 V.

The raw Sentaurus high-voltage plot provides the same current-source observable
at matching points:

| drain bias (V) | Vela / Sentaurus source | Vela / Sentaurus Id | Sentaurus Iava / Id |
|---:|---:|---:|---:|
| 6.4 | 0.985473 | 1.016351 | 0.937893 |
| 6.8 | 0.977372 | 1.011348 | 1.012233 |
| 6.9 | 0.975348 | 1.010309 | 1.031815 |
| 7.0 | 0.973416 | 1.009350 | 1.051626 |

Interpolating the Sentaurus `Iava - Id` values at 6.7 and 6.8 V gives
6.734425890 V. The like-for-like current-source crossing difference is therefore
0.185647899 V (2.76%), not the difference to 6.377494 V.

Correction after checking the O-2018.06-SP2 manual and the original `n4` plot:
the previously reported 6.377494 V is **not** the
`BreakAtIonIntegral(3 1.)` threshold. It is a linear interpolation of
`Iava-Id` between the sparse adaptive points 5.939778 V and 7.107808 V in the
original ABA-coupled plot. The dense fixed-bias replay resolves the curved
observable and gives the 6.734425890 V crossing above. `BreakAtIonIntegral`
only terminates the original quasistationary run at 10.448267 V when its three
largest path integrals exceed one. Because the 6.4 V local source integral is
already within 1.45%, the remaining current-source crossing difference is
dominated by the bias evolution of source and terminal current, not a constant
cell/source-support factor.

## Reproducible artifacts

- Legacy nodal ledger:
  `build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela_validation/btbt_e2_iic_qf_vector_fixed6p4_rerun_20260805/nodal_source_audit`
- Conservative mapping replay:
  `build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela_validation/btbt_e2_iic_qf_vector_conservative_fixed6p4_20260805`
- Nodal P1 replay:
  `build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela_validation/btbt_e2_iic_qf_vector_nodal_p1_fixed6p4_20260805`
- High-voltage accepted states and fixed points:
  `build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela_validation/btbt_e2_iic_qf_vector_branch_6p5_7p1_20260805`,
  `btbt_e2_iic_qf_vector_fixed6p8_20260805`,
  `btbt_e2_iic_qf_vector_fixed6p9_20260805`, and
  `btbt_e2_iic_qf_vector_fixed7p0_20260805`.

## Verification

`test_impact_ionization.exe` passes all 50 test cases and 668 assertions,
including conservative-area closure, nodal-source conservation, mode
validation, carrier-specific Eparallel, and Fermi-Dirac SG diagnostics.
