# SingleDevice fixed-state Eq. 231 residual decomposition (2026-08-13)

## Scope and result

The global all-material Eq. 231 assembly now exports the initial residual at
cell-vertex support. The export records the stiffness, squared-gradient, and
reaction terms separately, and aggregates the same contributions by node,
region, and material-interface neighborhood. Both endpoints use the imported
Sentaurus checkpoint for `psi`, `phin`, `Lambda`, and the reconstructed
potential-like restart state; the diagnostic does not use the subsequent
Vela drift-diffusion update to define the audited Eq. 231 drive.

The failure is spatially and algebraically localized. It is not a general
silicon-channel residual and it is not led by the `2*Lambda/C` reaction term.
The dominant contribution is `theta*|grad(w)|^2` in the narrow gate-oxide
strip next to `R.Polygate`:

| endpoint | max free-node residual | squared-gradient share of component L1 | gate-oxide share of cell-total L1 | max node |
|---|---:|---:|---:|---:|
| linear, `Vg=2.2 V`, `Vd=0.1 V` | `1.05199e6` | `96.763%` | `96.648%` | `2074` |
| saturation, `Vg=2.2 V`, `Vd=1.1 V` | `1.05961e6` | `96.735%` | `96.486%` | `2074` |

The same mesh node and the same gate-oxide cell family dominate both biases.
Cell `412` is the largest local contributor in both states. It is adjacent to
the `R.Gateox|R.Polygate` interface and has centroid approximately
`(-0.00119622, -0.0192708)` in mesh coordinates. Immediately adjacent
gate-oxide cells that are not themselves interface cells have nearly the same
residual. Consequently the interface-neighborhood share is about `51.4%`, but
the effect propagates through the thin oxide column and is not confined to the
single row of triangles touching the material edge.

On the cell-support L1 metric, insulator cells contribute `99.932%` at the
linear endpoint and `99.927%` at saturation. The largest transport-only
aggregated node residual is approximately `140`, compared with the all-domain
peak near `1.05e6`. This is an attribution of the currently assembled
equation, not a claim that Sentaurus has the same insulator residual.

## Term-level evidence

| endpoint | stiffness L1 | squared-gradient L1 | reaction L1 | cell-total L1 |
|---|---:|---:|---:|---:|
| linear | `1.49665e6` | `4.50455e7` | `1.01272e4` | `4.50617e7` |
| saturation | `1.49003e6` | `4.44944e7` | `1.18968e4` | `4.45102e7` |

The largest linear gate-oxide cell has
`|grad(w)|^2 = 1.309e26 1/m^2`; its squared-gradient contribution is about
`3.57e5` on each local row. The local `w` values span roughly `35` to `160`
over an oxide element only about `0.0113 nm` wide. Saturation produces the
same pattern and scale. The reaction contribution at the dominant node is
zero and remains below one in the largest oxide cells.

## Interpretation

This closes the broad hypothesis "the remaining failure is somewhere in the
multi-material nonlinear/interface discretization." The immediate fault is
more specific: the current all-material transformation feeds a very steep
potential-like drive into the insulator `theta*|grad(w)|^2` term. Because the
dominant oxide response is approximately bias invariant and many neighboring
non-interface triangles reproduce it, additional damping or restart fields
cannot remove it. The next implementation task is to rederive and implement
the `xi=eta=0` insulator-side Eq. 231 transformed variable and its interface
trace, then add a two-material manufactured solution that verifies cancellation
of the electrostatic contribution before rerunning the two endpoint oracles.

This diagnostic establishes localization and term dominance. It does not by
itself establish the final corrected Sentaurus-equivalent formula because the
commercial solver's internal integration-point residual is not exported.

## Reproduction and artifacts

The two endpoint configurations enable `residual_diagnostic_prefix` and
`residual_diagnostic_use_initial_state`. Run the configurations with
`build-release/vela_example_runner.exe`, then create the compact report with:

```powershell
python scripts/summarize_singledevice_eq231_residual.py `
  --root build-release/reference_tcad/singledevice_sentaurus2018/vela_import_fixedmaterials/vela/reports/eq231_fixed_state_residual_20260813
```

The generated report directory contains `lin_*` and `sat_*` cell, node,
region, and scalar summaries plus `summary.json` and `summary.md`.

Verification:

- density-gradient unit tests: 36 assertions in 10 cases passed;
- Newton density-gradient configuration tests: 33 assertions in 2 cases
  passed;
- both fixed endpoint probes emitted all residual artifacts and reproduced the
  existing expected nonconvergence of the experimental all-material solve.
