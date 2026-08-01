# PN2D node-volume policy independent code review

Date: 2026-08-01

Reviewed state: HEAD `7623ee8` plus the current uncommitted worktree.

Review mode: independent, read-only code review. The reviewer did not modify
the repository and did not read or depend on the scientific-review verdict.

Verdict:

```text
APPROVE_WITH_CONDITIONS
```

This verdict approves implementing and re-reviewing a minimal atomic patch. It
does not approve the current worktree as a production-default change.

Current worktree authorized: **no**.

## Findings

### C1 - the actual default patch does not exist

Priority: P1; blocks current authorization.

- `configs/templates/pn2d_bv.template.json` does not render
  `mesh_geometry.node_volume_policy`.
- `avalanche_current_support_profile` still defaults to
  `legacy_cell_reconstructed`.
- The template version remains 2.
- Neither the template nor `scripts/generate_pn2d_config.py` appears in the
  current `git diff`.

There is no actual production-default patch to approve.

### C2 - SG/Laux and node-volume policy cannot currently be bound atomically

Priority: P1.

`BV_AVALANCHE_CURRENT_SUPPORT_PROFILES` contains only the three avalanche
current-support fields. The renderer updates only
`solver.impact_ionization`, and its validation checks only that group. It
cannot enforce:

- SG/Laux always paired with `mixed_voronoi`;
- the legacy profile always paired with explicit `barycentric`;
- rejection of SG/Laux+barycentric and legacy+mixed half-migrations.

### C3 - the acceptance aggregator is not bound to the actual default render

Priority: P1.

The contract candidate contains both mixed-Voronoi and SG/Laux, but
`review_pn2d_node_volume_policy_default_acceptance.py` accepts no template,
render manifest, or base-config argument. Its aggregate result depends on
curve/control reports, forward IV, and the CTest log. It can therefore report
`ready_for_independent_default_policy_reviews` while the checked-in template
is still legacy and has no `mesh_geometry` field.

The evidence demonstrates that explicitly prepared candidate runs are healthy;
it does not prove that the actual default renderer emits that candidate.

### C4 - existing tests cannot reject a missing or half-applied patch

Priority: P2.

- Template tests still expect the legacy default.
- Profile rollback/opt-in tests inspect only the three impact fields.
- There are no negative tests for SG/Laux+barycentric or legacy+mixed.
- The new acceptance unit test covers the typed slope outcome, not default
  render binding.
- The contract focused-release list lacks atomic default-render and rollback
  tests.

### C5 - the current worktree is not an atomic default change

Priority: P2.

The worktree contains unrelated Newton solver, runner, diagnostic, report, and
test changes, plus multiple untracked files. A future template patch must be
isolated and reviewed as its own scope; the current worktree must not be
treated as one atomic production-default commit.

## Required patch shape

One BV profile selector should atomically render:

| Profile | Node volume | Current approximation | Source mapping | Midpoint density |
|---|---|---|---|---|
| default candidate | mixed_voronoi | element_edge_sg_gss_laux | element_vertex_box_measure | bernoulli |
| legacy rollback | barycentric | cell_reconstructed | triangle_gss_gradqf_truncated | gss_logistic |

The patch must:

1. modify only the PN2D BV template, renderer/validator, and their tests;
2. change the BV default profile to SG/Laux;
3. have the same selector write `mesh_geometry.node_volume_policy`;
4. have the legacy profile explicitly write `barycentric`;
5. reject both half-migration combinations;
6. increment the template version from 2 to 3;
7. preserve the C++ `BoxGeometryBuilder::Options` barycentric default;
8. preserve omitted `mesh_geometry` as barycentric;
9. leave PN2D IV, generic solver defaults, and geometry algorithms unchanged.

## Required tests and binding

The patch review requires:

- default BV render asserts mixed-Voronoi and the complete SG/Laux group;
- one legacy override asserts explicit barycentric and the complete legacy
  group;
- both half-migration combinations fail closed;
- manifest checks cover profile, version, and full pairing without hidden
  opt-in;
- the IV default render remains unchanged;
- omitted, empty-object, and explicit-barycentric compatibility remains;
- the acceptance aggregator binds template hash, render manifest, base config,
  and both M0/M2 executed configs;
- mutation tests make a missing node policy, legacy default restoration,
  explicit hidden opt-in, or half-migration fail acceptance;
- atomic default and rollback tests join the focused Release gate.

Current verification remains healthy: 12/12 focused CTest and 9/9 Python tests
passed. In a split full run, tests 1-496 passed before the shell timeout and
tests 497-509 passed 13/13 on immediate rerun. These results validate the
current opt-in base, not a nonexistent default patch.

## Evidence reviewed

- `configs/templates/pn2d_bv.template.json`
- `scripts/generate_pn2d_config.py`
- `scripts/review_pn2d_node_volume_policy_default_acceptance.py`
- `src/simulation/ConfigParsing.cpp`
- `src/mesh/BoxGeometryBuilder.cpp`
- `include/vela/mesh/BoxGeometryBuilder.h`
- `tests/regression/test_pn2d_config_templates.py`
- `tests/regression/test_pn2d_node_volume_policy_default_acceptance.py`
- `tests/test_box_geometry.cpp`
- `tests/test_dc_sweep.cpp`
- `docs/validation/contracts/pn2d_node_volume_policy_default_acceptance_v1.json`
