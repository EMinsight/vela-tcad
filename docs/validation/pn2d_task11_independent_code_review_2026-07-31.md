# PN2D Task 11 independent code review

Date: 2026-07-31

Reviewed snapshot:
`ec003470cb100cea8285b9b05c3396503b0995c8`

Review mode: independent, read-only code review. The reviewer did not modify
the repository and did not use the independent scientific-review conclusion.

Verdict:

```text
APPROVE_WITH_CONDITIONS
```

No Critical or High severity implementation defect was found in the current
opt-in SG/Laux path. A separate production-default patch may be designed and
reviewed, but the current review does not authorize switching the default.

## Findings

### C1 - a PN2D BV template switch must be an atomic configuration change

Priority: P1; blocks direct default switching.

The current BV template contains a coherent triangle-source configuration:

```text
current_approximation = cell_reconstructed
source_mapping_mode = triangle_gss_gradqf_truncated
cell_reconstructed_midpoint_density = gss_logistic
```

The exact-lattice candidate runner changes the configuration atomically:

```text
current_approximation = element_edge_sg_gss_laux
source_mapping_mode = element_vertex_box_measure
cell_reconstructed_midpoint_density = bernoulli
```

Changing only `current_approximation` would produce a noncanonical
configuration and is rejected by
`validateImpactIonizationConfiguration()` in
`include/vela/equation/AssemblerUtils.h`.

A future template patch must change the complete group and add a fresh-render
end-to-end test.

### C2 - the phrase production default has two distinct scopes

Priority: P1; blocks an ambiguous default patch.

The PN2D BV template currently selects the triangle reconstruction explicitly.
The global C++ `ImpactIonizationModelConfig` defaults are a different coherent
legacy group:

```text
generation = carrier_density
drivingForce = electric_field
currentApproximation = mobility_density_gradient
sourceMappingMode = node_F_node_alpha_node_G
```

Changing only the C++ `currentApproximation` default would break omitted-field
and partial configurations by forming a noncanonical SG/Laux group.

The first production proposal should be limited to the PN2D BV template.
Changing the global C++ API defaults requires a separate compatibility and
migration contract with omitted-field tests.

### C3 - current tests prove the opt-in operator, not a new default

Priority: P2.

Existing tests cover:

- Tri3 geometry and orientation invariance;
- edge-node mapping;
- nonnegative obtuse support and exact-area conservation;
- source records and node mapping;
- automatic differentiation and finite-difference derivatives;
- full residual/Jacobian/source-only Jacobian consistency;
- diagnostic versus assembled source consistency.

A default-change patch must additionally cover:

1. freshly rendered `pn2d_bv` default fields;
2. explicit legacy `cell_reconstructed` configurations preserving their old
   behavior;
3. omitted-field C++ compatibility;
4. rollback to the previous template default;
5. default-render M0 and M2 end-to-end runs under the prospective contract.

## Confirmed implementation properties

- Residual assembly and diagnostic postprocessing use the same local current
  source path.
- Full assembly and source-only Jacobian both use
  `elementEdgeGssLauxAvalancheSourceIntegralsLocal`.
- Electron and hole continuity equations use the same combined source with
  consistent residual/Jacobian signs.
- The canonical configuration guard covers current approximation, source
  mapping, generation model, driving force, and QFP-gradient discretization.
- Newton and Gummel parsing both propagate `current_approximation`.
- BV and forward-IV mobility doping bases remain isolated:
  `net_doping` and `cell_reconstructed_total_impurity`, respectively.
- Commit `ec00347` contains five documentation files only. At review time the
  worktree contained only the pre-existing untracked `tmp/`; no generated
  simulation output or tracked modification was present.

## Required scope of the next code proposal

The next proposal may:

- change only the PN2D BV template's atomic SG/Laux configuration group;
- retain the C++ API defaults;
- preserve an explicit legacy template/configuration path;
- provide an explicit rollback;
- add the default-path tests and prospective M0/M2 runs listed above.

The next proposal may not:

- change the global C++ default implicitly;
- change only one SG/Laux configuration field;
- treat this review as approval of an unimplemented patch;
- weaken an existing numerical or scientific threshold.

## Evidence reviewed

- `configs/templates/pn2d_bv.template.json`
- `configs/templates/pn2d_iv.template.json`
- `include/vela/physics/ImpactIonizationModel.h`
- `include/vela/equation/AssemblerUtils.h`
- `include/vela/equation/ElementEdgeGssLauxAD.inl`
- `src/equation/CoupledDDAssembler.cpp`
- `src/solver/NewtonSolver.cpp`
- `src/solver/GummelSolver.cpp`
- `scripts/run_pn2d_bv_exact_lattice_process.py`
- `tests/test_element_edge_gss_laux_avalanche.cpp`
- `tests/test_impact_ionization.cpp`
- `docs/validation/pn2d_task11_regression_review_2026-07-31.md`
