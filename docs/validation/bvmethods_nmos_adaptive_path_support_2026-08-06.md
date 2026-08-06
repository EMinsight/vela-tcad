# BVmethods NMOS adaptive path-support audit

## Implemented policy

`sentaurus_eparallel_adaptive` keeps the established nodal-local-maximum seed
generation, numbered peak groups, best-vertex alpha sampling, and arithmetic
WriteAll mean. It changes only the trajectory direction selected for each
peak:

- interior transport peaks follow the majority-carrier reconstructed SG
  current, with the existing electric-field reliability fallback;
- peaks on the transport boundary, including their two-ring P1 aliases,
  follow the minority-carrier quasi-Fermi gradient;
- net doping selects electron minority paths in p-type material and hole
  minority paths in n-type material.

The selection contains no seed IDs, bias thresholds, avalanche parameter
changes, or target-integral fitting.

## Common-state result at 10.4482667308 V

| rank | Vela electron | Sentaurus electron | Vela hole | Sentaurus hole | Vela mean | Sentaurus mean | mean error |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 1.515512276 | 1.71545 | 2.060742713 | 1.90119 | 1.788127494 | 1.80832 | -1.12% |
| 2 | 1.318392593 | 1.39435 | 1.766196767 | 1.76272 | 1.542294680 | 1.578535 | -2.30% |
| 3 | 1.278886169 | 1.28569 | 1.675239385 | 1.65037 | 1.477062777 | 1.46803 | +0.62% |

Rank 2 is shortened from 0.262254 um to 0.214739 um. Rank 3 is extended from
0.117040 um to 1.348466 um. The original rank-3 hole deficit
`1.314998 -> 1.65037` becomes `1.675239 -> 1.65037`, a +1.51% residual.

## Full-branch result

All 16 frozen states from 7.0 V to 10.4482667308 V replay successfully. The
terminal common-state agreement is substantially improved, but the full
BreakAtIonIntegral branch is not closed:

- Vela adaptive rank 3 is 0.968844 at 7.0 V and exceeds one at 7.2 V;
- the Sentaurus 7.0 V third distinct MeanIonIntegral plateau is approximately
  0.001759;
- Vela therefore splits and numbers the surface high-field corridor too early,
  even though the final three path amplitudes are now close.

`corridor_deduplicated` removes the exact rank-2/rank-3 alias at low bias, but
still leaves the p-type surface path at 0.813352 at 7.0 V. Fixed geometric
deduplication is therefore insufficient. The remaining issue is the
bias-dependent local-peak prominence and retention rule, not alpha, mobility,
Mean ordering, or carrier path integration.

## Artifacts

- final-state run:
  `build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela_validation/adaptive_minority_qf_10p448_20260806`;
- final-state exact WriteAll comparison:
  `build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela_validation/adaptive_minority_qf_mean_compare_20260806`;
- 16-state replay:
  `build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela_validation/adaptive_minority_qf_branch_20260806`;
- fixed-corridor control:
  `build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela_validation/adaptive_corridor_audit_20260806`.

## Next root-cause target

Export the seed field and saddle field for each surface maximum at 7.0, 8.0,
9.0, 10.0, and 10.448 V. Compare peak prominence against the Sentaurus path
inventory, then merge maxima connected above the validated saddle fraction.
The criterion must be fixed across bias and mesh states before it is allowed to
control `BreakAtIonIntegral(3 1.)`.
