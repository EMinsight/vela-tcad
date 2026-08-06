# BVmethods NMOS nodal colocated P1 avalanche audit (2026-08-05)

## Scope

This audit closes the location mismatch that remained after introducing the
`nodal_vertex_star` electric-field recovery.  The frozen-state
`nodal_eparallel_p1` path now evaluates, at every transport node:

1. the electric-field vector recovered from the complete semiconductor vertex
   star;
2. the conventional electron and hole current vectors recovered from the same
   node star;
3. carrier-specific `Eparallel` and ionization coefficients;
4. carrier generation multiplied by the exact semiconductor P1 measure,
   `sum(incident triangle area / 3)`.

The mode is rejected unless the configuration uses `postprocess_only`,
`eparallel`, `current_density`, `nodal_vector_current_reconstructed`, and
`eparallel_field_recovery=nodal_vertex_star`.  It therefore cannot silently
fall back to the earlier edge-adjacent field stencil.

The SG avalanche diagnostic CSV now exports the electron and hole source
partitions for both edge endpoints.  Summing those partitions reconstructs the
exact nodal P1 source ledger; no symmetric edge-to-node projection is needed.

## Verification

- `test_impact_ionization`: 51 test cases, 693 assertions passed.
- `test_dc_sweep`: 80 test cases, 2982 assertions passed.
- Fixed-state replay: accepted 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7.0, 7.1,
  and 7.2 V states.
- The generated configurations explicitly contain
  `eparallel_field_recovery=nodal_vertex_star` and
  `source_mapping_mode=nodal_eparallel_p1`.

## Current-source diagnostic

On the unchanged converged DD states, the colocated P1 source gives:

| Bias (V) | Iava / Id |
|---:|---:|
| 6.4 | 0.882550 |
| 6.7 | 0.930542 |
| 7.0 | 0.979571 |
| 7.1 | 0.996169 |
| 7.2 | 1.012896 |

Linear interpolation brackets the current-source crossing at
`7.121442639 V`.  Relative to the previous node-star-field/edge-source result,
the exact P1 source is 0.8415% lower at 6.4 V and 1.1042% lower at 7.0 V; the
crossing moves upward by 68.168 mV from `7.053274758 V`.

This current-source crossing is a diagnostic and is not the official
Sentaurus `BreakAtIonIntegral(3 1.)` path criterion at `6.377494 V`.

## Sentaurus spatial comparison

Using the exact carrier-specific endpoint ledger and the common semiconductor
P1 measures:

- integrated electron avalanche source, Vela/Sentaurus, falls from 0.96526 at
  6.4 V to 0.94542 at 7.0 V;
- peak electron generation, Vela/Sentaurus, is 0.90377 at 6.4 V and 0.89965 at
  7.0 V;
- the region above 10% of the peak grows by 1.05087x in Vela versus 1.17522x
  in Sentaurus between 6.4 and 7.0 V;
- the largest remaining growth deficit is in the high-field corridor around
  x=0.079--0.097 um and y=0.004--0.014 um.

The location mismatch is therefore removed, but it does not explain the
remaining BV discrepancy.  The remaining evidence points to insufficient
bias growth of the high-field generation shoulder and to the distinction
between the global current-source crossing and Sentaurus path ionization
integrals.

## Artifacts

- Fixed-state branch:
  `build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela_validation/btbt_e2_iic_qf_vector_nodal_vertex_star_p1_colocated_branch_6p4_7p1_20260805`
- Current-source summary:
  `analysis/branch_closure/summary.json`
- Exact P1 spatial support summary:
  `analysis/generation_support_exact_p1/summary.json`
