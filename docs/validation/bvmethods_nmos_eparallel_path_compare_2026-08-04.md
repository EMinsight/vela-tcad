# BVmethods NMOS Eparallel path comparison, 2026-08-04

## Outcome

The Vela postprocessed `Avalanche(Eparallel)` branch was continued from the
validated 7 V checkpoint to the official Sentaurus terminal bias
`10.4482667308 V`.  All 27 checkpoints in the final 7.9--10.448 V segment
converged.  Vela did not meet `BreakAtIonIntegral(3 1.0)`: its third mean path
integral was only `0.374585` at the Sentaurus stopping voltage.

The official `n4_des.tdr` contains nine distinct positive
`MeanIonIntegral` plateaus.  The three largest plateaus were exported as the
raw Sentaurus path-support coordinates and compared against Vela's three
ordered edge paths at the same terminal bias.

## Reproducible artifacts

- Vela converged sweep and accepted states:
  `build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela_validation/eparallel_extend_linear_7p9_10p45_20260804/postprocess_only`
- Imported Sentaurus TDR geometry and path fields:
  `build-release/reference_tcad/bvmethods_sentaurus2018/run01/sentaurus_path_export_20260804/neutral`
- Same-bias comparison products:
  `build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela_validation/eparallel_path_compare_10p448_20260804`
- Comparison implementation:
  `scripts/compare_bvmethods_nmos_path_ionization.py`

## Rank-wise physical comparison

All field and alpha maxima below are evaluated on each solver's own rank-N
path support.  This table intentionally preserves rank semantics instead of
forcing a spatial assignment.

| Rank | Vela max Eparallel (V/m) | Sent max E (V/m) | E ratio | Vela max alpha e/h (1/m) | Sent max alpha e/h (1/m) | Vela integral e/h/mean | Sent integral e/h/mean | Mean ratio |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 1.6603e8 | 1.8871e8 | 0.880 | 0 / 2.4204e7 | 3.6128e7 / 2.7367e7 | 0 / 0.75109 / 0.37554 | 1.71545 / 1.90119 / 1.79456 | 0.209 |
| 2 | 1.5718e8 | 2.1677e8 | 0.725 | ~0 / 2.2854e7 | 3.8414e7 / 3.0538e7 | ~0 / 0.74935 / 0.37468 | 1.39435 / 1.76272 / 1.54419 | 0.243 |
| 3 | 1.7570e8 | 3.6249e8 | 0.485 | 0 / 2.5601e7 | 4.8299e7 / 3.9411e7 | 0 / 0.74917 / 0.37459 | 1.28569 / 1.65037 / 1.42877 | 0.262 |

The hole alpha deficit is moderate on ranks 1--3: Vela maxima are about
`88%`, `75%`, and `65%` of Sentaurus.  The decisive discrepancy is the
electron branch: Vela's selected Eparallel paths carry essentially zero
electron alpha, whereas Sentaurus has `3.61e7--4.83e7 1/m` and electron
injection integrals `1.29--1.72`.

## Geometry comparison

The Sentaurus coordinates below are the raw node support of each distinct
ion-integral plateau.  They are authoritative TDR data, but they should not be
misrepresented as a unique interpolated centerline.

| Rank | Vela x extent (um) | Vela y extent (um) | Sent x extent (um) | Sent y extent (um) | Vela length (um) | Sent support nodes |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 0.2119--0.2825 | 0--1.000 | 0.1413--0.1766 | 0.0140--0.2112 | 1.0633 | 19 |
| 2 | 0.2472--0.4238 | 0--1.000 | -0.0177--0.1942 | 0.00005--0.2488 | 1.1126 | 73 |
| 3 | 0.1766--0.2825 | 0--1.000 | 0--0.1236 | 0.00005--0.1361 | 1.1126 | 55 |

Vela traces every monotone-potential path to a mesh boundary.  Sentaurus path
support remains localized near the drain junction and terminates by roughly
`y=0.25 um`.  The minimum-distance geometry assignment is Vela 1 -> Sentaurus
2, Vela 2 -> Sentaurus 1, and Vela 3 -> Sentaurus 3, so rank ordering also does
not identify the same spatial paths.  Even under that best assignment the
symmetric mean support distances are `0.152--0.195 um` and Hausdorff distances
are `0.768--0.873 um`.

## Root-cause localization

1. **Path termination semantics are different.** Vela currently follows a
   monotone graph to the device boundary.  The official TDR shows that
   Sentaurus stops the useful breakdown-path support in the high-field drain
   junction.  A Sentaurus-compatible stop-field/path-search termination rule
   is missing in Vela.
2. **The Eparallel carrier branches are not represented independently.** The
   Vela paths selected by the shared graph are hole dominated and suppress the
   electron alpha branch to zero.  Sentaurus reports nonzero electron and hole
   injection integrals on every top-three path.  Vela needs carrier-specific
   Eparallel support and a deliberate rule for combining/ranking those paths.
3. **Seed and spatial interpolation semantics remain different.** Sentaurus
   reports use of the best element vertex for the impact-ionization model,
   whereas Vela seeds paths from local edge-field maxima and uses
   piecewise-constant edge alpha.  This can explain the spatial rank mismatch
   and part of the remaining hole-alpha deficit.
4. **Peak field alone is not the closure blocker.** Rank-1 Vela path field is
   already 88% of Sentaurus and its hole-alpha maximum is 88%, yet its mean
   integral is only 21%.  Geometry/termination and the missing electron branch
   dominate before a coefficient refit should be considered.

## Code and diagnostics added

`sweep.diagnostics.path_ionization_integrals` now supports
`segments_csv_file`.  The ordered per-segment export includes coordinates,
endpoint potential, Eparallel, electron/hole alpha, alpha-ds, cumulative
alpha-ds, prefix injection integrals, and full-path injection integrals.

The comparison script exports:

- `sentaurus_top3_path_support_nodes.csv`: raw coordinates and physics fields;
- `sentaurus_top3_path_summary.csv`: plateau integrals and support maxima;
- `rankwise_physics_compare.csv`: direct rank-wise field/alpha/integral table;
- `vela_sentaurus_matched_segments.csv`: nearest-support local comparison;
- `vela_sentaurus_path_summary_compare.csv`: best geometry assignment;
- `summary.json`: machine-readable provenance and conclusions.

## Verification

- `test_path_ionization_integral`: 17 assertions in 4 test cases passed.
- `test_dc_sweep "[path_ionization]"`: 18 assertions in 2 test cases passed,
  including the new ordered segment CSV end-to-end test.
- `python -m py_compile scripts/compare_bvmethods_nmos_path_ionization.py`
  passed.
- `git diff --check` passed.

## Next implementation milestone

The next closure step should implement Sentaurus-compatible path-search
termination and independent electron/hole Eparallel path construction before
changing avalanche coefficients.  The same 10.448 V state can then be rerun
and evaluated against the exported top-three path support without another
Sentaurus simulation.
