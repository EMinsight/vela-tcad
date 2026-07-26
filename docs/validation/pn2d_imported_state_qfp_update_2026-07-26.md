# PN2D imported-state continuity residual and first-QFP-update audit

Date: 2026-07-26

Status: Task 7 complete; Task 8 not authorized.

## Decision

The typed outcome is `operator_improvement_without_qfp_causality`.

The opt-in element-edge source strongly reduces the high-field residual and first update at `-20 V` on both Minimal6 and coarse7x3, and it strongly improves coarse7x3 at `-10 V`. However, at the first non-negligible Minimal6 source difference (`-10 V`) it worsens both carrier residual norms by about `2.7%` and all carrier-only/coupled QFP-update norms by about `2.4%`, while coarse7x3 improves by orders of magnitude. Therefore the source-operator change does not improve the same residual and first-update QFP error across both topology classes at the first departure. Task 7 does not authorize Task 8.

Tasks 8 and 9 are not entered. The production default remains unchanged; `element_edge_sg_gss_laux` remains opt-in. No field, mobility, alpha, or geometry scale was fitted. Van Overstraeten, QFP-gradient driving force, and Masetti plus QFP-gradient high-field mobility remain unchanged.

## Frozen inputs and topology gate

The audit uses exact imported `psi`, electron QFP, and hole QFP at `-1`, `-10`, and `-20 V` for:

- Minimal6 `mirror` and `sketch`;
- the regenerated coarse7x3 semiconductor mesh from the frozen explicit-GradQF Sentaurus root.

The coarse probe contains 33 raw vertices, of which 27 are semiconductor mesh vertices. All 27 Vela node IDs and coordinates match exactly (`0` maximum coordinate error), and all 32 element IDs and local Tri3 vertex permutations match exactly. The six non-semiconductor probe vertices are excluded. This avoids the older 21-node/24-cell coarse export, which is not topology-identical to the frozen 27-node/32-cell oracle.

The Task 6 area-conservation normalization is causal on the constrained-obtuse RED mesh. Minimal6 and the regenerated coarse mesh are right-triangle controls with no pre-fix truncated-area closure defect, so that normalization itself cannot explain their QFP updates. Task 7 compares the current canonical triangle source with the full opt-in element-edge source; it does not call box reconstruction a native directed edge current.

## Configurations

The production-triangle branch is the canonical configuration:

- `model: van_overstraeten`;
- `driving_force: quasi_fermi_gradient`;
- `generation: current_density`;
- `current_approximation: cell_reconstructed`;
- `cell_reconstructed_midpoint_density: gss_logistic`;
- `source_mapping_mode: triangle_gss_gradqf_truncated`.

The candidate is explicitly opt-in:

- `current_approximation: element_edge_sg_gss_laux`;
- `source_mapping_mode: element_vertex_box_measure`;
- `quasi_fermi_gradient_discretization: cell_gradient`.

Both branches use unchanged Vela Masetti mobility with QFP-gradient high-field driving. Separate avalanche-off and SRH-off probes change only the requested term family. Native C++ SG edge diagnostics provide the box-flux reconstruction and incident row scale.

## Sealed evidence

Two independently generated roots have byte-identical hashes for all five sealed CSVs:

- `build-release/pn2d-imported-state-qfp-update-20260726-a`;
- `build-release/pn2d-imported-state-qfp-update-20260726-b`.

| Evidence | Rows |
|---|---:|
| topology gate | 9 |
| residual decomposition | 468 |
| carrier-only and coupled first QFP updates | 600 |
| analytic/FD Jacobian blocks | 90 |
| avalanche-off and SRH-off controls | 468 |

The independent verifier does not import analyzer calculations. It confirms:

| Gate | Result |
|---|---:|
| sealed A/B hashes | identical |
| residual decomposition maximum relative error | `1.6654229802532582e-16` |
| boundary branch maximum absolute difference | `0` |
| topology maximum coordinate error | `0 um` |
| exact biases present | `-1, -10, -20 V` |
| update modes present | carrier-only, coupled |
| non-improved topology/bias/mode/carrier groups | present in Minimal6 and coarse7x3 |
| Task 8 authorized | `false` |

Normalized residual terms are direct C++ carrier-term diagnostics. Physical residual and row-scale columns multiply the unit-scaled residuals by `C0*D0`, with the same unit-scaling constants used by the assembler. The residual rows record SG divergence, SRH, avalanche, gauge, boundary/Dirichlet, incident row scale, normalized residual, physical residual, and final assembled residual.

## Residual and update causality

The table gives candidate/production L2 ratios. Values below one improve.

| Topology | Bias | electron residual | hole residual | carrier-only electron update | carrier-only hole update | coupled electron update | coupled hole update |
|---|---:|---:|---:|---:|---:|---:|---:|
| Minimal6 mirror/sketch | -1 | `1.000000` | `1.000000` | `1.000000` | `1.000000` | `1.000000` | `1.000000` |
| Minimal6 mirror/sketch | -10 | `1.027156` | `1.027976` | `1.023601` | `1.024308` | `1.023601` | `1.024308` |
| Minimal6 mirror/sketch | -20 | `3.409014e-5` | `3.373837e-5` | `0.0297334` | `0.0294365` | `0.0297335` | `0.0294366` |
| coarse7x3 | -1 | `1.000000` | `1.000000` | `1.000000` | `1.000000` | `1.000000` | `1.000000` |
| coarse7x3 | -10 | `1.099286e-7` | `1.099887e-7` | `0.00158892` | `0.0140316` | `0.00106766` | `0.00924058` |
| coarse7x3 | -20 | `0.0180704` | `0.0181181` | `0.0120429` | `0.00857828` | `0.00891144` | `0.00639328` |

Because each initial state is the exact Sentaurus QFP state, the pre-update QFP error is zero and the post-update error is `abs(delta_qfp)`. Thus a smaller first update is the direct first-step error metric; no target fitting is involved.

The first material source-support departure is `-10 V`. Minimal6 and coarse7x3 move in opposite causal directions there. Avalanche-off controls remove the branch difference, while SRH-off controls retain the source-support difference, locating the discriminator in avalanche current/source support rather than SRH, boundary conditions, or mobility configuration. Boundary rows remain bitwise identical between source branches.

## Jacobian evidence

Task 6's focused central-FD test remains the release gate and passes at `<=1e-8` nonzero relative error, with a frozen absolute gate for near-zero sources. Its diagnostic replay and assembled residual use identical mobility and driver configuration.

The broader imported-state matrix probe exposes an additional unresolved diagnostic: for non-negligible `sg_avalanche` blocks its maximum reported relative difference is `0.004877997535809977`, above `1e-8`. The largest row is coarse7x3, `-10 V`, opt-in. This is reported as a failed imported-state diagnostic gate, not hidden as a near-zero case and not used to weaken the Task 6 threshold. It must be resolved before any future authorization, but the cross-topology QFP causality failure already blocks Task 8 independently.

## Task 10 decision ledger

The applicable decision-matrix row is: general fixed-state operator evidence exists, but the first QFP update does not improve consistently across topologies. Allowed action: keep the operator opt-in as a diagnostic and investigate the QFP/source-support branch. No default-change proposal is made.

Tasks 8 and 9 are recorded as `not_entered_task7_gate_failed`; no Minimal6 curve is labeled a physical BV curve, and no fine PN2D knee claim is made.

## Reproduction

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
$env:PYTHONPATH = "."
python scripts\diagnose_pn2d_imported_state_qfp_update.py `
  --runner build-release\vela_example_runner.exe `
  --minimal-phase-e-root build-release\pn2d-minimal6-phase-e-continuity-sourcefix2-20260724-a `
  --coarse-log build-release\pn2d-general-tri3-sentaurus-avalanche-controls-20260725-a\coarse7x3\explicit_grad_qf\fetched\run_explicit_grad_qf.out `
  --coarse-mesh-root build-release\pn2d-general-tri3-task7-imported-mesh-20260726\vela `
  --output-root build-release\pn2d-imported-state-qfp-update-20260726-a
python scripts\verify_pn2d_imported_state_qfp_update.py `
  --root-a build-release\pn2d-imported-state-qfp-update-20260726-a `
  --root-b build-release\pn2d-imported-state-qfp-update-20260726-b `
  --output build-release\pn2d-imported-state-qfp-update-20260726-a\independent_verification.json
```

Generated `build-release` roots are validation artifacts and are not committed.
