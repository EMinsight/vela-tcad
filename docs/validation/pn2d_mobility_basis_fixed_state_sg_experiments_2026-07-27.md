# PN2D mobility-basis and fixed-state SG experiments

Date: 2026-07-27

## Purpose

These experiments test two hypotheses behind the remaining PN2D forward-bias
differences from Sentaurus:

1. Mobility in compensated junction nodes should use total impurity rather than
   signed net doping.
2. After fixing the Sentaurus state, the Vela SG operator may still have a
   material discretization or geometric residual.

The production default remains `net_doping`. The two new opt-in candidates are:

- `total_impurity`
- `cell_reconstructed_total_impurity`

## Experiment 1: self-consistent forward IV

All three candidates converged at all 201 bias points from 0 V to 20 V.

| Mobility doping basis | Electron mobility RMSE, nodes 7-15 (cm2/V/s) | Hole mobility RMSE, nodes 7-15 (cm2/V/s) | Median absolute current error at 1/2/5/10/15/20 V | Current error at 20 V | Potential RMS error at 20 V |
|---|---:|---:|---:|---:|---:|
| `net_doping` | 452.364 | 103.394 | 3.911% | +3.785% | 0.285801 V |
| `total_impurity` | 48.220 | 14.190 | 0.299% | -0.289% | 0.026305 V |
| `cell_reconstructed_total_impurity` | **8.695** | **1.872** | **0.271%** | **-0.262%** | **0.022758 V** |

Relative to `net_doping`, the cell-reconstructed candidate reduces the
20 V current-error magnitude by 93.1%, the potential RMS error by 92.0%, and
the junction electron-mobility RMSE by 98.1%.

At 20 V, the reconstructed electron mobilities on nodes 7-15 are within about
9 cm2/V/s of the imported Sentaurus values. The corresponding hole-mobility
error is below about 2.3 cm2/V/s.

## Experiment 2: Sentaurus fixed-state SG audit

Sentaurus `psi`, electron/hole quasi-Fermi potential, electron density, and
hole density were imported at 1, 2, 5, 10, 15, and 20 V. Vela then evaluated
edge SG currents and nodal continuity residuals without solving or modifying
the state.

The 20 V results are representative of all six bias points:

| Mobility doping basis | Junction electron edge-current ratio | Junction hole edge-current ratio | Global normalized L1 residual, e/h | Junction normalized L1 residual, e/h | Maximum nodal residual, e/h |
|---|---:|---:|---:|---:|---:|
| `net_doping` | 1.1564 | 1.0926 | 5.368% / 3.186% | 11.030% / 6.867% | 18.395% / 11.357% |
| `total_impurity` | 0.9931 | 0.9955 | 1.011% / 0.626% | 2.364% / 1.457% | 5.868% / 3.777% |
| `cell_reconstructed_total_impurity` | **0.9914** | **0.9959** | **0.440% / 0.257%** | **1.029% / 0.599%** | **2.092% / 1.335%** |

The reconstructed candidate therefore removes about 91% of the junction
continuity mismatch caused by the net-doping mobility basis. The remaining
fixed-state discrepancy is sub-percent globally and about 1% in the junction,
which bounds the residual attributable to SG/control-volume discretization and
state-import differences.

Sixteen diagonal edges have zero box coupling in this orthogonal dual mesh.
Their raw SG current density is non-zero, but their continuity support is
exactly zero. This is not by itself an error: a zero dual-face measure can be
geometrically valid. The result means these edges should remain an explicit
diagnostic category rather than be treated as the primary root cause.

## Implementation and recommendation

- The mobility configuration now accepts `doping_concentration_basis`.
- The selected basis is applied consistently to SG transport, VTK mobility
  diagnostics, triangle avalanche Real/AD paths, and element-edge local AD.
- The fixed-state audit supports an opt-in general Tri3 scope and preserves
  donor/acceptor fields and effective intrinsic density/BGN.
- The legacy/default `net_doping` behavior is unchanged.

For PN2D Sentaurus parity, `cell_reconstructed_total_impurity` is the best
candidate. It should not become the global default until reverse-bias,
breakdown, and avalanche regression matrices have been rerun because the
mobility field also enters the impact-ionization source and its Jacobian.

## Reproduction

```powershell
python scripts/run_pn2d_forward_mobility_doping_basis_experiment.py
python scripts/run_pn2d_sentaurus_fixed_state_sg_audit.py
ctest --test-dir build-release --output-on-failure
```

Primary result files:

- `build-release/pn2d-forward-mobility-doping-basis-20260727/candidate_summary.csv`
- `build-release/pn2d-forward-mobility-doping-basis-20260727/exact_anchor_current_comparison.csv`
- `build-release/pn2d-forward-mobility-doping-basis-20260727/junction_mobility_comparison_20V.csv`
- `build-release/pn2d-sentaurus-fixed-state-sg-audit-20260727/fixed_state_summary.csv`
- `build-release/pn2d-sentaurus-fixed-state-sg-audit-20260727/fixed_state_edge_comparison.csv`

Validation: all 485 CTest cases passed.
