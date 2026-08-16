# SingleDevice Eq. 231 AverageBox and PolySi/SiO2 closure (2026-08-15)

## Outcome

The two requested endpoint rows were decomposed and closed far enough to pass
the external one-step quantum-state acceptance gate.

- Pure-Silicon node 1816 was not a bulk source-model error.  It touches a
  100.73-degree triangle, and Vela used the positive mixed-Voronoi obtuse
  fallback while Sentaurus `AverageBoxMethod` retained the signed
  circumcentric element-vertex measure.
- PolySilicon/SiO2 node 2075 retained the same control-volume measure in both
  solvers.  Central NewtonPlot perturbations instead isolated a small
  material-side affine reaction-trace mismatch.
- The linear and saturation imported endpoints now move by `0.126 mV` and
  `0.157 mV`, respectively.  Both are below the `0.5 mV` acceptance limit.

This result unlocks the 21-point self-consistent Id-Vg curve runs; those curve
runs are intentionally not part of this row-localization task.

## Node 1816: direct AverageBox evidence

The VM `MeasureCoefficients.debug` file was compared cell-by-cell with the
current Vela Eq. 231 diagnostic.  At node 1816:

- fitted stiffness: `-0.2847325298255031`;
- old positive mixed-area reaction: `0.25256618642470924`;
- Sentaurus AverageBox reaction: `0.2847325215184572`;
- reconstructed residual: `-8.307e-9`.

The same correction reduces the saturation residual at node 1816 to
`-5.631e-9`.  The previous insulator-corner multiplier `0.9713027871` was a
compensation for the same incorrect positive-area fallback and was restored
to its neutral value `1.0`.

## Node 2075: central Jacobian oracle

Positive and negative `10 uV` eQuantumPotential perturbations were run for
the nonzero row support `2075, 2072, 2073, 2074, 2120`.  The four off-diagonal
Vela/Sentaurus NewtonPlot scale factors were:

- node 2072: `-295.2907666`;
- node 2073: `-295.0776292`;
- node 2074: `-295.2878070`;
- node 2120: `-295.2907742`.

The common support and scale confirm the fitted edge operator.  The 2075
self derivative gives a scaled total diagonal of `8804.38447`, compared with
a Vela fitted-edge diagonal of `8782.62857`; the inferred Sentaurus reaction
diagonal is `21.75589795`.

A constrained least-squares fit combined this diagonal oracle with all 65
PolySilicon/SiO2 interface rows in both endpoint states.  The selected
material-side trace parameters are:

- PolySilicon weight: `1.0684933639683336`;
- oxide-at-PolySilicon weight: `2.569027176700638`;
- PolySilicon offset: `+0.01872581675079906 V`;
- oxide-at-PolySilicon offset: `-0.0052046570150173915 V`.

After implementation, node 2075 residuals are `-0.0587802` in the linear
state and `-0.0150189` in saturation.  The global fixed-state maxima are
`0.0965465` and `0.148589`; neither controls the endpoint acceptance gate.

## Verification

- `test_density_gradient_quantum_potential.exe`: 97 assertions in 23 cases
  pass, including a new obtuse-triangle signed AverageBox manufactured root.
- `test_newton_solver.exe`: 1165 assertions in 83 cases pass.
- Linear endpoint: inner quantum solve converged; raw imported-state change
  `0.000126 V`.
- Saturation endpoint: inner quantum solve converged; raw imported-state
  change `0.000157 V`.

The one-point executable reports overall failure only because these oracle
configs deliberately set `outer_max_iterations=1`; the acceptance metric is
the logged pre-relaxation raw change, not full outer-loop convergence.

## Reproducibility

- AverageBox audit: `scripts/audit_singledevice_eq231_average_box_measure.py`
- TDR perturbation runner:
  `scripts/run_sentaurus_singledevice_eq231_central_probe.sh`
- constrained fit:
  `scripts/fit_singledevice_eq231_poly_oxide_central_closure.py`
- local oracle output:
  `build-release/reference_tcad/singledevice_sentaurus2018/eq231_poly_sio2_2075_central_20260815/central_closure_fit.json`

