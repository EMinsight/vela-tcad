# PN2D Minimal6 element-edge avalanche replay

Date: 2026-07-25

Status: `valid_diagnostic_replay_with_closure`

## Scope

This audit resolves the discrete-current support used by Sentaurus
O-2018.06-SP2 for the default avalanche model on the Minimal6 mesh.

The selected lattice is:

- topologies: `mirror`, `sketch`;
- reverse biases: `-1 V`, `-10 V`, `-20 V`;
- 6 physical semiconductor vertices, 4 triangles, 9 global edges, and
  12 element-local edges per state;
- Sentaurus branches: implicit default, `AvalDensGradQF`, and
  `ElementVolumeAvalanche`.

The experiment is post-processing and remote-oracle diagnostics only. No
production formula in `include/` or `src/` was changed.

## Runtime support boundary

The CurrentPlot Tcl interface returned:

- vertex potential, carrier density, electron and hole QFP, mobility,
  avalanche coefficient, and carrier-split avalanche generation;
- element mobility, electric field, QFP gradients, and current-density
  vectors;
- element-local-edge `ReadCoefficient`;
- element-local-vertex `ReadMeasure`.

An isolated runtime probe attempted to read
`eAvalancheGeneration` on `ElementVertex-RegionWise` support. Sentaurus
terminated with:

```text
Tried to read undefined ElementVertex-RegionWise eAvalancheGeneration !
```

Therefore native element-vertex avalanche generation is unavailable through
this interface. Element-vertex source contributions in this audit are typed
reconstructions made from native vertex generation and exact `ReadMeasure`.

## Geometry

The four triangles return local box-coefficient permutations of:

```text
(1, 0, 0.25)
(0, 0.25, 1)
(1, 0.25, 0)
(0.25, 0, 1)
```

The two diagonal hypotenuse edges have zero box coefficient. Each triangle
returns element-vertex measures that are permutations of:

```text
0.0625, 0.125, 0.0625 um^2
```

and sum to the exact `0.25 um^2` triangle area.

The zero diagonal coefficient applies to box-flux assembly. It does not
remove the physical SG current density on that triangle edge before
element-vector reconstruction.

## Element-edge SG replay

For every element-local edge from `i` to `j`, the runtime probe recomputed
the carrier edge current density before multiplication by the box
coefficient. The same current was also multiplied by the exact box support
to reproduce the conservative continuity flux.

The apparent native element-vector current is approximately six decades
larger than the SG edge current on this mesh. It must not be substituted for
the default avalanche-current support.

## Reconstruction candidates

Five candidates were evaluated:

1. `gss_laux_edge_volume_weighted`: pairwise edge-vector constructions and
   GSS/Laux dual-volume weighting;
2. `charon_whitney_hcurl_cell_average`: lowest-order Whitney/HCurl mapping
   corresponding to the Charon SGCVFEM design;
3. `genius_least_squares_tangent`: all-three-edge tangent least squares,
   matching the Genius triangle reconstruction;
4. `box_active_edge_exact`: exact reconstruction from the two positive-box
   edges of each right triangle;
5. `native_element_vector_control`: direct native element vector, retained
   only as a rejected control.

For each candidate and carrier:

```text
G_n = alpha_n * |J_n| / q
G_p = alpha_p * |J_p| / q
```

was distributed with the exact element-local-vertex measures and compared
with native carrier-split Sentaurus generation.

| Candidate | Samples | Median integral error, dex | Maximum, dex |
|---|---:|---:|---:|
| GSS/Laux | 12 | `5.111719e-5` | `4.292374e-4` |
| Two active box edges | 12 | `5.111719e-5` | `4.292374e-4` |
| Charon Whitney/HCurl | 12 | `1.260280e-3` | `2.855924e-3` |
| Genius tangent least squares | 12 | `1.543705e-3` | `3.651102e-3` |
| Native element-vector control | 12 | about `6.15218` | rejected |

For this right-triangle mesh, GSS/Laux and the two-active-edge
reconstruction are numerically identical. This identifies the current
support used by the default Sentaurus avalanche source to within less than
`5e-4 dex` over both topologies, both carriers, and all three biases.

The native nodal generation comparison for GSS/Laux has a maximum node error
of about `0.0187 dex` at `-1 V`; it falls below `0.006 dex` at `-10 V` and
below `0.0032 dex` at `-20 V`. The much tighter integrated closure shows
that the remaining node difference is source redistribution, not current
magnitude.

## Existing Vela triangle proxy

After converting both sources to the same two-dimensional integral unit,
the existing Vela triangle proxy differs from native Sentaurus by:

| Bias | Absolute source error, dex |
|---:|---:|
| `-1 V` | `1.108073` |
| `-10 V` | `10.609656` |
| `-20 V` | `11.950516` |

This confirms that the production triangle QFP-gradient proxy is not the
Sentaurus default element-edge SG avalanche-current support.

## Conservation and integration closure

| Gate | Maximum result |
|---|---:|
| ReadMeasure integral vs CurrentPlot Integrate | `1.265365e-15` relative |
| carrier-specific terminal current | `2.318564e-3` relative |
| total terminal current | `1.137692e-7` relative |
| internal total-current KCL | `3.163255e-9` relative |

The carrier-specific maximum is the small cathode hole current at `-20 V`.
Electron-plus-hole total current gives the physically relevant terminal
gate.

## Sentaurus control branches

`ElementVolumeAvalanche` is bitwise identical to the default branch for
both topologies and all selected biases, including CurrentPlot files.

`AvalDensGradQF` changes the source as follows:

| Bias | Total source log10 ratio to default, dex | Maximum terminal relative change |
|---:|---:|---:|
| `-1 V` | `-0.356851` | `0` |
| `-10 V` | `+0.232399` | `3.22697e-9` |
| `-20 V` | `+0.262428` | `2.17036e-4` |

`mirror` and `sketch` produce the same control metrics within floating-point
roundoff.

The control confirms that the default Sentaurus target must remain the
element-edge SG branch. `AvalDensGradQF` is useful for sensitivity analysis
but is not the production parity target.

## Decision

The current-support cause is now identified:

1. Sentaurus default avalanche does not use the exported native element
   `CurrentDensity` vector.
2. It is reproduced by element-local SG edge currents followed by a
   GSS/Laux-equivalent reconstruction on this right-triangle mesh.
3. Zero-box-coefficient diagonal edges do not contribute to continuity
   flux, but they remain defined edge-current observations during the
   general three-edge vector reconstruction.
4. For the current Minimal6 geometry, the GSS/Laux result reduces exactly to
   reconstruction from the two positive-box edges.
5. The remaining Vela avalanche mismatch is therefore an implementation
   support mismatch, not evidence for changing the Van Overstraeten alpha
   law, mobility law, or QFP sign convention.

This evidence supports implementing a separately named diagnostic
element-edge SG/GSS-Laux avalanche-current mapping in Vela. It does not
authorize silently changing the existing production mapping.

## Evidence

Deterministic evidence root:

`build-release/pn2d-minimal6-element-avalanche-replay-20260725`

Key files:

- `analysis/element_edges.csv`
- `analysis/element_reconstructions.csv`
- `analysis/node_generation_replay.csv`
- `analysis/state_source_summary_corrected.csv`
- `analysis/box_current_closure.csv`
- `analysis/source_integral_closure.csv`
- `controls/control_vs_default.csv`
- `independent_verification.json`

The independent verifier passed all 240 element-reconstruction rows, source
formulas, hashes, GSS/Laux equivalence, source-integral closure,
terminal/KCL gates, control comparisons, and the typed unsupported
element-vertex observation with zero failures.
