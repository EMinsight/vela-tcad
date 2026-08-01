# PN2D node-volume policy independent scientific review

Date: 2026-08-01

Review mode: independent, read-only scientific review. The reviewer did not
modify the repository and did not read or depend on the code-review verdict.

Verdict:

```text
APPROVE_WITH_CONDITIONS
```

The scientific evidence authorizes a narrowly scoped proposal: the PN2D BV
template may explicitly render `mesh_geometry.node_volume_policy =
mixed_voronoi` for the validated, non-obtuse M0/M2 PN2D mesh family. It does
not authorize changing the global parser default or claiming equivalence to
Sentaurus `MixAverageBoxMethod` on arbitrary obtuse/non-Delaunay meshes.

## Blocking findings

### S1 - an unrestricted arbitrary-mesh default is not supported

Priority: P1 for an unrestricted template.

The PN2D BV template accepts a replaceable `mesh_file`, so an unconditional
template default is not intrinsically limited to M0/M2. The controlled obtuse
experiment found:

- maximum local Sentaurus MixAverage versus Vela mixed difference:
  `6.0337e-3 um^2`;
- assembled node-volume L1 difference: `9.5674e-3 um^2`.

The current local half/quarter/quarter obtuse rule is therefore not a general
implementation of Sentaurus MixAverage truncation. This finding becomes
nonblocking only if the authorization is explicitly limited to qualified
non-obtuse PN2D grids, or if the final template path includes a fail-closed
mesh qualification gate.

### S2 - the exact accepted physics bundle must be retained

Priority: P1.

The frozen avalanche-on candidate jointly fixes:

```text
mesh_geometry.node_volume_policy = mixed_voronoi
impact_ionization.current_approximation = element_edge_sg_gss_laux
impact_ionization.source_mapping_mode = element_vertex_box_measure
impact_ionization.cell_reconstructed_midpoint_density = bernoulli
```

The current PN2D BV template still defaults to `legacy_cell_reconstructed`.
Changing only node volume would create the unvalidated combination
`mixed_voronoi + legacy_cell_reconstructed`; the existing on/knee evidence
cannot authorize that combination. The eventual default render must match the
full accepted bundle, or the exact new combination must be rerun.

## Evidence findings

### Prospective ordering and binding

The reviewer verified the local artifact sequence:

| Artifact | UTC timestamp |
|---|---:|
| frozen contract | 2026-08-01 12:02:46 |
| M0 gate | 2026-08-01 12:06:08 |
| M2 gate | 2026-08-01 12:21:49 |
| forward-IV report | 2026-08-01 12:32:30 |
| Release log | 2026-08-01 12:39:42 |
| aggregate acceptance | 2026-08-01 12:42:18 |

The reports bind contract SHA-256
`3a2d65879d1d7446a259afb6f81af9c17e7da2686536e9f6832eb5b5810f6c22`.
File timestamps are not an immutable provenance system, but the ordering and
embedded hashes are adequate for this local prospective audit.

M0 and M2 use directly converted Sentaurus TDR inputs with complete global
vertex mapping:

- M0 TDR SHA-256 `999057b8...32bac0`, 27 vertices and 32 triangles;
- M2 TDR SHA-256 `5b52f9d1...7c889`, 115 vertices and 191 triangles.

### Reverse-bias controls

Every off, IIC, and SG/Laux-on branch completed `29/29` exact-lattice points in
two runs. IV and per-bias state hashes are deterministic. IIC is identical to
off in IV and state hashes, excluding IIC source feedback leakage.

| Grid | off RMSE | off maximum | on RMSE | on maximum | knee RMSE | abs V_break error |
|---|---:|---:|---:|---:|---:|---:|
| M0 | 8.876e-5 dex | 1.562e-4 dex | 0.001766 dex | 0.003749 dex | 0.002810 dex | 0.001 V |
| M2 | 3.868e-5 dex | 9.223e-5 dex | 0.001923 dex | 0.004647 dex | 0.003051 dex | 0.001 V |

Both exact lattices are monotonic. M0 has two slope crossings separated by
`0.00082 V`; M2 has the predeclared shared no-crossing typed outcome.

Maximum global continuity closure ratios remain below the `0.01` contract
limit: approximately `2.200e-3` on M0, `9.138e-5` on M2 off/IIC, and
`1.917e-6` on M2 on.

### Forward and regression controls

- Forward IV converged at `201/201` points.
- The two mixed runs have identical IV hashes.
- Six-anchor Sentaurus median/maximum errors are `0.2700%/0.4066%`.
- Maximum degradation over barycentric is `2.42e-14`.
- The frozen Release acceptance recorded `509/509` tests passing.

### Direct box-measure evidence

Both actual meshes have maximum angle 90 degrees and zero obtuse elements.

| Grid | max local difference | max assembled node difference | area sum |
|---|---:|---:|---:|
| M0 | 1.041e-17 um^2 | 1.388e-17 um^2 | both 1.0 um^2 |
| M2 | 1.041e-17 um^2 | 1.388e-17 um^2 | both 1.0 um^2 |

This changes the M0/M2 control-volume statement from an output-based inference
to a direct measurement.

## Mandatory scope limits

The final proposal must:

1. change only the PN2D BV template path;
2. preserve the global omitted-policy barycentric behavior;
3. retain explicit legacy/barycentric rollback;
4. render the complete accepted SG/Laux plus mixed-volume bundle atomically;
5. exclude or fail closed on unqualified obtuse/non-Delaunay meshes;
6. avoid claiming general Sentaurus MixAverage equivalence;
7. avoid bundled solver, forward-IV, physics-parameter, or threshold changes.

## Evidence reviewed

- `docs/validation/contracts/pn2d_node_volume_policy_default_acceptance_v1.json`
- `docs/validation/pn2d_node_volume_policy_default_acceptance_2026-08-01.md`
- `docs/validation/sentaurus_box_measure_direct_export_2026-08-01.md`
- `build-release/pn2d-node-volume-default-acceptance-v1-20260801/acceptance.json`
- `build-release/pn2d-node-volume-default-acceptance-v1-20260801/M0/gate_report.json`
- `build-release/pn2d-node-volume-default-acceptance-v1-20260801/M2/gate_report.json`
- `build-release/pn2d-box-measure-probe-20260801/audit/m0_direct_compare.json`
- `build-release/pn2d-box-measure-probe-20260801/audit/m2_direct_compare.json`
- `build-release/pn2d-box-measure-probe-20260801/audit/box_measure_audit.json`
