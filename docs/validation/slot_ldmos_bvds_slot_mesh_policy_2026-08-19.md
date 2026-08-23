# Slot-LDMOS BVDS SLOT boundary and mesh policy

Date: 2026-08-19

## Decision

The imported SProcess topology is preserved exactly for the first Vela BVDS baseline. The `SLOT` contact is represented as a fixed-potential electrostatic conductor boundary (`metal_gate`), not as an ohmic carrier contact. The imported process mesh is accepted only with the complete `legacy_cell_reconstructed` avalanche profile; `element_edge_sg_gss_laux` remains fail-closed because this mesh does not satisfy its non-obtuse qualification.

## SLOT boundary contract

| Property | Resolved value | Reason/evidence |
|---|---|---|
| Imported owner region | `Oxide_1` | The final `!Gas` TDR collapses the Tungsten body and retains the contact edge set. |
| Contact nodes | 62 | Preserved without remapping from the neutral Sentaurus export. |
| Incident material cells | 93 `SiO2`, 14 `Nitride` | No incident Silicon or PolySilicon transport cells were found. |
| Vela contact type | `metal_gate` | Pins electrostatic potential while intentionally applying no electron/hole Dirichlet boundary condition. |
| BVDS voltage | 0 V | Matches the source SDevice electrode definition. |
| Zero-bias check | max absolute SLOT potential = 0 V | Verified from all 62 SLOT nodes in the emitted VTK result. |

Using an ohmic or generic Dirichlet carrier contact for `SLOT` would incorrectly inject or remove carriers through an insulator boundary. The same electrostatic-only classification is used for `gate`; `source`, `drain`, and `substrate` remain ohmic semiconductor contacts.

## Coordinate and topology contract

The neutral export labels its columns `x_um` and `y_um`, but the values from this SProcess TDR are in centimetres. The preparation workflow therefore applies an explicit scale of `1e4` before writing Vela's micrometre mesh:

| Quantity | Raw TDR/export bounds | Vela bounds |
|---|---:|---:|
| x | -4.163438e-5 to 8.0e-4 cm | -0.4163438 to 8.0 um |
| y | 0 to 4.69e-4 cm | 0 to 4.69 um |

The sealed conversion preserves 13,881 vertices, 27,275 triangles, all five contact node sets, material assignments, and `NetActive` doping values.

## Mesh policy

| Item | Selected policy | Qualification result |
|---|---|---|
| Topology | Exact TDR topology | Preserves the Sentaurus oracle geometry and contacts. |
| Control volumes | `barycentric` | Positive on all nondegenerate triangles. |
| Box geometry | Negative-cotangent fallback enabled | 883 negative cotangent contributions use the positive `area/(3*edge.length)` fallback. |
| Avalanche current | `cell_reconstructed` | Part of the legacy atomic profile. |
| Avalanche source mapping | `triangle_gss_gradqf_truncated` | Part of the legacy atomic profile. |
| Midpoint density | `gss_logistic` | Part of the legacy atomic profile. |
| Non-obtuse requirement | `false` | Required for the exact imported topology. |
| `element_edge_sg_gss_laux` | Forbidden | 883 strictly obtuse cells; maximum angle 177.072209764 degrees. |

This is a baseline compatibility policy, not a claim that the mesh is high quality. Its maximum edge aspect ratio is 1291.58 and minimum angle is 0.04434 degrees. Any later use of the PN2D `element_edge_sg_gss_laux` bundle requires a new mesh that passes its complete non-obtuse qualification; the SG/GSS-Laux bundle must not be partially enabled on this mesh.

## Executed validation

| Check | Result |
|---|---|
| Preparation regression (`python -m unittest tests.regression.test_prepare_slot_ldmos_bvds`) | Passed, 2 tests. |
| CTest registration (`slot_ldmos_bvds_preparation`) | Passed, 1/1 in independent `build-slot-prep`. |
| Real export preparation | Passed; generated mesh, doping, materials, boundary-check and BVDS configurations. |
| SG/GSS-Laux selection on the real mesh | Rejected by preparation with an explicit mesh-qualification error. |
| Poisson SLOT boundary check | Converged; 27,275 cells, zero degenerate cells, 883 fallbacks; SLOT and gate potentials exactly 0 V. |
| Full legacy BVDS configuration at 0 V | Converged; one DC point, zero degenerate cells, 883 fallbacks. |

## Artifacts and remaining scope

The reproducible preparation entry point is `scripts/prepare_slot_ldmos_bvds.py`; its behavior is protected by `tests/regression/test_prepare_slot_ldmos_bvds.py`. The generated production bundle is under `build-release/reference_tcad/slot_ldmos_sentaurus2022/run01/vela_ready`.

SLOT semantics and the initial exact-topology mesh strategy are resolved. They do not yet establish quantitative BVDS equivalence. Remaining P0 work is the Sentaurus two-dimensional resistor/current unit calibration and staged BVDS reference run. The provisional PolySilicon/Nitride material entries and the IALMob model delta must also be tracked before making a full-physics equivalence claim.
