# Slot-LDMOS Sentaurus BVDS staged reference

Date: 2026-08-20

Sentaurus release: T-2022.03-SP2

## Result

All seven staged decks completed with exit code 0 on the Sentaurus virtual machine. The two-dimensional current/resistance convention is calibrated, the avalanche-off, IIC and self-consistent avalanche controls are available, and the final external-resistor run reached its requested drain-current stop.

| Stage | Final outer drain voltage | Final inner drain voltage | Final drain current | Result |
|---|---:|---:|---:|---|
| Equilibrium | 0 V | 0 V | -5.9567e-24 A/um | Converged |
| 1 Tohm, 1 V calibration | 1 V | 0.999409604 V | 5.90255e-16 A/um | Converged |
| Direct 1 V control | 1 V | 1 V | 5.90523e-16 A/um | Converged |
| Avalanche off | 60 V | 59.9934120 V | 6.58796e-15 A/um | Converged |
| IIC/AvalPostProcessing | 60 V | 59.9934120 V | 6.58796e-15 A/um | Converged; IV matches avalanche-off |
| Self-consistent avalanche | 60 V | 38.5048827 V | 2.14951e-11 A/um | Converged; internal voltage is avalanche-clamped |
| Final external resistor | 109076.963 V | 38.5215816 V | 1.09038e-7 A/um | Stopped normally after exceeding 1e-7 A/um |

Linear interpolation of the final curve at exactly `1e-7 A/um` gives an outer source voltage of 100038.521 V and an internal drain voltage of 38.5209012 V. The report therefore uses `BVDS = 38.5209 V` for this current criterion. The approximately 100 kV outer voltage is the drop required across the 1 Tohm two-dimensional series resistor and is not the device breakdown voltage.

## Two-dimensional unit calibration

The 1 Tohm branch at 1 V gives:

- `Vouter - Vinner = 5.90395780855e-4 V`
- `Idrain * R = 5.90254985148e-4 V`
- load-line residual = `1.40795707e-7 V`
- direct/resistor drain-current relative difference = `4.53032e-4` (0.0453%)

Consequently, the SDevice two-dimensional deck uses terminal current in `A/um` and contact resistance in `ohm*um`. Vela must transcribe the original `Resistor=1e12` as `resistance_ohm_um=1e12`, with no additional width multiplication for the reference comparison.

## Sentaurus treatment of obtuse and non-Delaunay elements

The Sentaurus Device T-2022.03 User Guide distinguishes a geometrically obtuse element from a non-Delaunay element. The imported mesh has many obtuse triangles under a simple angle audit, but SDevice reports only two non-Delaunay elements out of 27,275 cells, no interface non-Delaunay elements, and total `DeltaVolume=1.5e-5%`.

| Manual mechanism | Applicability to this run | Decision |
|---|---|---|
| AverageBoxMethod obtuse-element truncation | Default SDevice box construction handles obtuse control volumes | Active through the normal SDevice discretization |
| `ElementVolumeAvalanche` | Renormalizes truncated element-vertex avalanche volumes so their sum equals the geometric element volume, avoiding exaggerated avalanche generation and premature breakdown | Already present in the original BVDS deck; retained |
| `AvalDensGradQF` | Uses the quasi-Fermi-gradient current approximation for avalanche, recommended for improved power-device stability | Already present in the original BVDS deck; retained |
| `AvalFlatElementExclusion=1..2` | Intended only for nearly flat avalanche-producing semiconductor elements | Not enabled: the worst Silicon triangle has maximum angle 139.74 degrees and minimum angle 3.08 degrees; the 177.07-degree worst cell is PolySilicon, where avalanche is not evaluated |
| `WeightedVoronoiBox` | Requires `DelVorWeight` arrays written by SProcess with `StoreDelaunayWeight` | Not enabled: the frozen process did not store these arrays, and the observed non-Delaunay volume error is already negligible |
| `BM_StableCalculation` | Optional accuracy enhancement for challenging sliver/non-Delaunay meshes | Keep as a sensitivity option; not required to reproduce the original deck or pass the present baseline |

The IIC best path reports a maximum field of 2.3515e6 V/cm at `(0.0060533, 4.42555) um`, with electron/hole ionization integrals 0.12874/0.0983688. The Silicon triangles nearest that point have approximately 35/55/90-degree angles, so the IIC hotspot is not located in the nearly flat or worst obtuse cells.

This Sentaurus handling does not change the separate Vela qualification rule: Vela's PN2D `element_edge_sg_gss_laux` bundle remains fail-closed on the imported mesh. The current Vela exact-topology baseline therefore continues to use the complete `legacy_cell_reconstructed` profile.

## Reproducible artifacts

The local result root is `build-release/reference_tcad/slot_ldmos_sentaurus2022/run01/sentaurus_bvds_result`. It contains:

- SHA-256-sealed archives for logs/curves and TDR states;
- extracted `.plt`, SDevice logs, exit codes, final TDR files and bias checkpoints;
- normalized per-stage CSV curves;
- `analysis/slot_ldmos_bvds_reference_summary.json` and the generated Markdown summary.

The tracked analyzer is `scripts/analyze_slot_ldmos_bvds_reference.py`, with regression coverage in `tests/regression/test_analyze_slot_ldmos_bvds_reference.py`.
