# Sentaurus box measure direct export and obtuse-mesh audit

Date: 2026-08-01
Sentaurus Device: O-2018.06-SP2
Status: read-only physics audit; no production default changed

## Question and outcome

This audit directly exports the Sentaurus Device element-vertex box measure and
tests the statement that the obtuse vertex of an obtuse triangle contributes a
negative control-volume area.

The result has two parts:

1. On the actual balanced M0 and M2 BV grids, the exported Sentaurus
   `Measure[k][j]` agrees with Vela `mixed_voronoi` to floating-point roundoff.
   Neither grid contains an obtuse triangle. This replaces the earlier indirect
   inference with a direct same-grid measurement.
2. On a controlled obtuse mesh, Sentaurus reports a negative signed
   non-Delaunay intersection coefficient (`CoeffIntersection = -0.4`), but its
   exported element-vertex `Measure[k][j]` values remain non-negative. The
   default `AverageBoxMethod` represents the problem as overlapping control
   volumes and violates total-volume conservation by 2%. `MixAverageBoxMethod`
   truncates those control volumes and restores exact total-volume conservation.

Therefore, “the obtuse contribution is negative” is accurate for the raw
circumcentric construction or the signed non-Delaunay intersection diagnostic,
but not for Sentaurus O-2018.06-SP2 `AverageBoxMethod`'s exported `Measure` array.

## Official export path

The Sentaurus Device User Guide O-2018.06, Chapter 36, pages 1008–1015 defines:

- `Measure[k][j]` as the control volume associated with local vertex `j` of
  element `k`, in 2-D units of μm²;
- `Coefficients[k][j]` as the box discretization coefficient;
- `BoxMeasureFromFile(GrdNumbering)` and
  `BoxCoefficientsFromFile(GrdNumbering)` as the grid-numbered debug-file path;
- creation of `MeasureCoefficients.debug` when the file does not already exist;
- `BM_CoeffIntersectionNonDelaunayElements` and related fields as plot variables
  for box-method statistics.

The guide also states that `AverageBoxMethod` coefficients are positive on a
non-Delaunay mesh, but total-volume conservation can fail. `MixAverageBoxMethod`
uses the AverageBox coefficients and truncates obtuse elements when computing
control volumes.

The probe adds the following to `Math` without changing device physics:

```text
AverageBoxMethod                       # or MixAverageBoxMethod
BoxMeasureFromFile(GrdNumbering)
```

It performs only one Poisson initialization and writes the BM diagnostic fields.

## Actual M0 and M2 results

| Quantity | M0 | M2 |
|---|---:|---:|
| Vertices | 27 | 115 |
| Triangles | 32 | 191 |
| Maximum triangle angle | 90° | 90° |
| Obtuse triangles | 0 | 0 |
| Geometry area | 1.0 μm² | 1.0 μm² |
| Sentaurus `Measure` sum | 1.0 μm² | 1.0 μm² |
| Negative Sentaurus `Measure` entries | 0 | 0 |
| Max local difference, Sentaurus vs Vela mixed | 1.041e-17 μm² | 1.041e-17 μm² |
| Max assembled node-volume difference | 1.388e-17 μm² | 1.388e-17 μm² |
| L1 assembled node-volume difference | 9.021e-17 μm² | 2.559e-16 μm² |
| L1 difference, Sentaurus vs barycentric | 2.083e-2 μm² | 9.668e-2 μm² |

Input TDR hashes:

- M0: `999057B81108361FDA6FFDF6C0D8CC40CDFC44E00B7F7B7F2B5B6D6D0A32BAC0`
- M2: `5B52F9D16454DB2FC9E34C44185A1C8BD8DA468631AD7AA214C769AD4A87C889`

Interpretation:

- The old barycentric node volumes are directly proven not to match Sentaurus.
- The opt-in Vela `mixed_voronoi` node volumes are directly proven to match
  Sentaurus on the actual M0/M2 grids.
- Consequently, any residual M2 SG/Laux-on mismatch after enabling
  `mixed_voronoi` should not be attributed to the node box measure itself. The
  next search should remain in transport/source mapping and the coupled feedback
  path.

## Controlled obtuse experiment

The synthetic mesh preserves the M0 topology and all node doping arrays. Only
vertex 1 is moved from `(0.25, 0)` to `(0.05, 0)` μm. This creates two obtuse,
non-Delaunay triangles (elements 1 and 3), with maximum angle 128.6598°.

| Quantity | AverageBoxMethod | MixAverageBoxMethod |
|---|---:|---:|
| Geometry area | 1.000 μm² | 1.000 μm² |
| Sum of exported `Measure` | 1.020 μm² | 1.000 μm² |
| Total-volume error | +2.0% | approximately 0 |
| Negative exported `Measure` entries | 0 | 0 |

Additional diagnostics:

- Sentaurus run log signed `CoeffIntersection`: `-0.4`.
- Exported `BM_CoeffIntersectionNonDelaunayElements`: magnitude `0.4` on the
  affected elements.
- Raw local circumcentric formula: 3 negative element-vertex areas.
- `MixAverageBoxMethod` versus Vela's current per-triangle
  half/quarter/quarter obtuse rule:
  - maximum local difference: `6.034e-3 μm²`;
  - L1 assembled node-volume difference: `9.567e-3 μm²`.

The last comparison is important: Vela `mixed_voronoi` is equivalent to the
Sentaurus measure on the current non-obtuse M0/M2 grids, but it is not a general
implementation of Sentaurus `MixAverageBoxMethod` for arbitrary non-Delaunay
meshes. Sentaurus's truncation is mesh-neighborhood aware; Vela currently applies
a local per-triangle positive split.

## Reproduction and artifacts

Tracked tools:

- `scripts/generate_sentaurus_box_measure_probe.py`
- `scripts/create_obtuse_dfise_probe.py`
- `scripts/audit_sentaurus_box_measure_probe.py`
- `scripts/compare_sentaurus_box_measure_to_vela.py`

Generated outputs remain under the ignored directory:

```text
build-release/pn2d-box-measure-probe-20260801/
```

Key machine-readable outputs:

```text
audit/m0_direct_compare.json
audit/m2_direct_compare.json
audit/box_measure_audit.json
audit/element_local_box_measure_compare.csv
```

The controlled obtuse TDR hash is
`1A6C6423C3C16FBF429913D2176A5674B5EE1A3C4825DD693E05C30F334930B5`.

## Recommendation

1. Keep the production default unchanged until the existing default-policy
   acceptance process finishes.
2. Treat `mixed_voronoi` as validated against Sentaurus box measures for the
   present M0/M2 grid family, which contains no obtuse elements.
3. Add a separate future mesh-policy task if arbitrary obtuse/non-Delaunay grids
   must match Sentaurus `MixAverageBoxMethod`. Do not claim that the current
   half/quarter/quarter fallback is equivalent.
4. For the current M2 BV discrepancy, stop treating node control-volume geometry
   as an open root cause once `mixed_voronoi` is active; continue with the
   transport/source/coupled-feedback diagnostics already planned.
