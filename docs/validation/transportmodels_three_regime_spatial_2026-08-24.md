# TransportModels DD/DG three-regime spatial-field comparison

Status: **complete**. The comparison uses exact shared node IDs and only strict interior nodes of `R.Substrate`.

Operating points: deep off `Vg=-1.00 V`, threshold/transition `Vg=0.12 V`, and on state `Vg=0.92 V`; all use `Vd=1.10 V`.

Vector and positive-magnitude fields use active-node symmetric percentage error. Signed/cross-zero fields use absolute error normalized by the Sentaurus p95 absolute field scale. Ec/Ev are compared after removing the one arbitrary global Sentaurus energy origin.

| Mode | Regime | Priority | Field | P95 error (%) | Max error (%) | Direction P95 (deg) |
|---|---|---|---|---:|---:|---:|
| DD | deep_off | P0 | Electron current density | 199.5 | 200 | 119 |
| DD | deep_off | P0 | Hole current density | 199.6 | 200 | 34.2 |
| DD | deep_off | P0 | Total current density | 199.6 | 200 | 101 |
| DD | deep_off | P1 | Electric field | 13.71 | 69.5 | 5.54 |
| DD | deep_off | P1 | Electron GradQuasiFermi | 192.3 | 198.4 | 162 |
| DD | deep_off | P1 | Hole GradQuasiFermi | 16.5 | 86.59 | 8.04 |
| DD | deep_off | P1 | Electron Eparallel | 179.4 | 200 | - |
| DD | deep_off | P1 | Hole Eparallel | 160.5 | 199.8 | - |
| DD | deep_off | P1 | Electron Enormal | 102.7 | 195.6 | - |
| DD | deep_off | P1 | Hole Enormal | 102.7 | 195.6 | - |
| DD | deep_off | P1 | Electron mobility | 78.85 | 184.9 | - |
| DD | deep_off | P1 | Hole mobility | 46.55 | 111.9 | - |
| DD | deep_off | P1 | SRH recombination | 1.734 | 6030 | - |
| DD | deep_off | P1 | Space charge | 5.88 | 30.75 | - |
| DD | deep_off | P1 | Band gap | 1.058 | 1.058 | - |
| DD | deep_off | P1 | Bandgap narrowing | 75.82 | 76.36 | - |
| DD | deep_off | P1 | Electron affinity | 2.35 | 2.381 | - |
| DD | deep_off | P1 | Conduction band | 1.6 | 1.628 | - |
| DD | deep_off | P1 | Valence band | 1.685 | 1.71 | - |
| DD | on | P0 | Electron current density | 31.03 | 171.6 | 9.32 |
| DD | on | P0 | Hole current density | 190.4 | 200 | 19.7 |
| DD | on | P0 | Total current density | 31.03 | 171.6 | 9.32 |
| DD | on | P1 | Electric field | 10.21 | 115.1 | 3.95 |
| DD | on | P1 | Electron GradQuasiFermi | 23.16 | 197.6 | 9.87 |
| DD | on | P1 | Hole GradQuasiFermi | 15.59 | 195.7 | 6.11 |
| DD | on | P1 | Electron Eparallel | 162 | 197.8 | - |
| DD | on | P1 | Hole Eparallel | 46.94 | 199.3 | - |
| DD | on | P1 | Electron Enormal | 100 | 193 | - |
| DD | on | P1 | Hole Enormal | 100 | 193 | - |
| DD | on | P1 | Electron mobility | 73.18 | 184.3 | - |
| DD | on | P1 | Hole mobility | 44.22 | 106.1 | - |
| DD | on | P1 | SRH recombination | 1.661 | 23.68 | - |
| DD | on | P1 | Space charge | 5.729 | 85.42 | - |
| DD | on | P1 | Band gap | 1.058 | 1.058 | - |
| DD | on | P1 | Bandgap narrowing | 76.25 | 90.38 | - |
| DD | on | P1 | Electron affinity | 2.35 | 2.381 | - |
| DD | on | P1 | Conduction band | 1.595 | 1.621 | - |
| DD | on | P1 | Valence band | 1.685 | 1.71 | - |
| DD | threshold | P0 | Electron current density | 58.43 | 199.8 | 14 |
| DD | threshold | P0 | Hole current density | 192.7 | 200 | 17.7 |
| DD | threshold | P0 | Total current density | 58.43 | 199.8 | 14 |
| DD | threshold | P1 | Electric field | 10.86 | 63.63 | 3.69 |
| DD | threshold | P1 | Electron GradQuasiFermi | 186.7 | 198.2 | 167 |
| DD | threshold | P1 | Hole GradQuasiFermi | 13.61 | 172.9 | 6.7 |
| DD | threshold | P1 | Electron Eparallel | 165.2 | 199.3 | - |
| DD | threshold | P1 | Hole Eparallel | 57.09 | 196 | - |
| DD | threshold | P1 | Electron Enormal | 99.86 | 196.6 | - |
| DD | threshold | P1 | Hole Enormal | 99.86 | 196.6 | - |
| DD | threshold | P1 | Electron mobility | 75.3 | 184.7 | - |
| DD | threshold | P1 | Hole mobility | 44.3 | 106.2 | - |
| DD | threshold | P1 | SRH recombination | 1.616 | 23.67 | - |
| DD | threshold | P1 | Space charge | 5.262 | 44.03 | - |
| DD | threshold | P1 | Band gap | 1.058 | 1.058 | - |
| DD | threshold | P1 | Bandgap narrowing | 75.82 | 76.36 | - |
| DD | threshold | P1 | Electron affinity | 2.35 | 2.381 | - |
| DD | threshold | P1 | Conduction band | 1.601 | 1.627 | - |
| DD | threshold | P1 | Valence band | 1.683 | 1.706 | - |
| DG | deep_off | P0 | Electron current density | 199.3 | 200 | 126 |
| DG | deep_off | P0 | Hole current density | 199.5 | 200 | 31.7 |
| DG | deep_off | P0 | Total current density | 199.4 | 200 | 111 |
| DG | deep_off | P1 | Electric field | 37.97 | 185.3 | 9.39 |
| DG | deep_off | P1 | Electron GradQuasiFermi | 192.6 | 198.4 | 167 |
| DG | deep_off | P1 | Hole GradQuasiFermi | 18.35 | 97.22 | 10.4 |
| DG | deep_off | P1 | Electron Eparallel | 186.4 | 200 | - |
| DG | deep_off | P1 | Hole Eparallel | 154.1 | 199.8 | - |
| DG | deep_off | P1 | Electron Enormal | 120.7 | 194.9 | - |
| DG | deep_off | P1 | Hole Enormal | 120.7 | 194.9 | - |
| DG | deep_off | P1 | Electron mobility | 94.56 | 184.9 | - |
| DG | deep_off | P1 | Hole mobility | 47.37 | 110.4 | - |
| DG | deep_off | P1 | SRH recombination | 7.373 | 4489 | - |
| DG | deep_off | P1 | Space charge | 5.576 | 11.48 | - |
| DG | deep_off | P1 | Band gap | 1.058 | 1.058 | - |
| DG | deep_off | P1 | Bandgap narrowing | 75.82 | 76.36 | - |
| DG | deep_off | P1 | Electron affinity | 2.35 | 2.381 | - |
| DG | deep_off | P1 | Conduction band | 1.596 | 1.644 | - |
| DG | deep_off | P1 | Valence band | 1.665 | 1.85 | - |
| DG | on | P0 | Electron current density | 38.59 | 153.6 | 11 |
| DG | on | P0 | Hole current density | 190 | 200 | 19.2 |
| DG | on | P0 | Total current density | 38.59 | 153.6 | 11 |
| DG | on | P1 | Electric field | 24.41 | 138.4 | 3.41 |
| DG | on | P1 | Electron GradQuasiFermi | 51.69 | 197.6 | 13.6 |
| DG | on | P1 | Hole GradQuasiFermi | 21.76 | 199.4 | 7.87 |
| DG | on | P1 | Electron Eparallel | 165 | 199.6 | - |
| DG | on | P1 | Hole Eparallel | 71.56 | 195.4 | - |
| DG | on | P1 | Electron Enormal | 107.7 | 192.6 | - |
| DG | on | P1 | Hole Enormal | 107.7 | 192.6 | - |
| DG | on | P1 | Electron mobility | 78.88 | 184.3 | - |
| DG | on | P1 | Hole mobility | 45.81 | 120.6 | - |
| DG | on | P1 | SRH recombination | 7.498 | 231.4 | - |
| DG | on | P1 | Space charge | 7.696 | 16.32 | - |
| DG | on | P1 | Band gap | 1.058 | 1.058 | - |
| DG | on | P1 | Bandgap narrowing | 75.82 | 76.36 | - |
| DG | on | P1 | Electron affinity | 2.35 | 2.381 | - |
| DG | on | P1 | Conduction band | 1.596 | 1.65 | - |
| DG | on | P1 | Valence band | 1.665 | 1.848 | - |
| DG | threshold | P0 | Electron current density | 72.48 | 190.4 | 14.2 |
| DG | threshold | P0 | Hole current density | 192.1 | 200 | 17.6 |
| DG | threshold | P0 | Total current density | 72.48 | 190.4 | 14.2 |
| DG | threshold | P1 | Electric field | 30.37 | 182.3 | 4.48 |
| DG | threshold | P1 | Electron GradQuasiFermi | 186.2 | 198.1 | 168 |
| DG | threshold | P1 | Hole GradQuasiFermi | 16.56 | 167.6 | 9.24 |
| DG | threshold | P1 | Electron Eparallel | 170 | 199.1 | - |
| DG | threshold | P1 | Hole Eparallel | 86.8 | 196.2 | - |
| DG | threshold | P1 | Electron Enormal | 113.3 | 191.5 | - |
| DG | threshold | P1 | Hole Enormal | 113.3 | 191.5 | - |
| DG | threshold | P1 | Electron mobility | 78.96 | 184.7 | - |
| DG | threshold | P1 | Hole mobility | 45.62 | 106.1 | - |
| DG | threshold | P1 | SRH recombination | 6.242 | 333.1 | - |
| DG | threshold | P1 | Space charge | 6.624 | 13.6 | - |
| DG | threshold | P1 | Band gap | 1.058 | 1.058 | - |
| DG | threshold | P1 | Bandgap narrowing | 75.82 | 76.36 | - |
| DG | threshold | P1 | Electron affinity | 2.35 | 2.381 | - |
| DG | threshold | P1 | Conduction band | 1.596 | 1.646 | - |
| DG | threshold | P1 | Valence band | 1.663 | 1.846 | - |

## Interpretation limits

- Vela current vectors are node-reconstructed diagnostics from the solved state; terminal-current/KCL acceptance remains the conservative current criterion.
- Local magnitude percentages are evaluated only where the Sentaurus magnitude exceeds `1e-3` of its case maximum; values below that active-field floor are retained in CSV but excluded from the primary percentage statistic.
- SRH and SpaceCharge cross zero, so their reported percentages are full-field scale-normalized errors rather than pointwise relative errors.

## Artifacts

- Total-current vector figure: `D:\code-repo\vela-tcad\build-release\reference_tcad\transportmodels_sentaurus2022\reports\transportmodels_three_regime_spatial_20260824\total_current_density_vectors.png`
- Transport-chain heatmap: `D:\code-repo\vela-tcad\build-release\reference_tcad\transportmodels_sentaurus2022\reports\transportmodels_three_regime_spatial_20260824\transport_chain_error_heatmap.png`
- SRH/charge/band heatmap: `D:\code-repo\vela-tcad\build-release\reference_tcad\transportmodels_sentaurus2022\reports\transportmodels_three_regime_spatial_20260824\source_charge_band_error_heatmap.png`
- Summary CSV: `D:\code-repo\vela-tcad\build-release\reference_tcad\transportmodels_sentaurus2022\reports\transportmodels_three_regime_spatial_20260824\field_error_summary.csv`
- Long-form comparison CSV: `D:\code-repo\vela-tcad\build-release\reference_tcad\transportmodels_sentaurus2022\reports\transportmodels_three_regime_spatial_20260824\spatial_field_comparison.csv`
