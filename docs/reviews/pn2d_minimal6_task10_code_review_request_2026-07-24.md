# PN2D Minimal6 Task 10 code review request

Date: 2026-07-24

Review status: requested; focused and full regressions pass.

## Review scope

Review the scoped changes for:

- mesh coordinate-unit handling in native source integration;
- production-faithful triangle source selection in Phase F;
- electric-field mobility candidate template;
- high-field box replay and factorial attribution;
- Phase E residual/first-update diagnostics;
- impact factorization and independent verifier;
- source-unit isolated audit; and
- validation reports and plan ledger.

## Required code checks

1. `integrate_native_nodal_per_unit_depth()` must default legacy meshes to
   meters and convert declared micrometer coordinates exactly once.
2. `_triangle_source_per_cm_s()` must sum the raw production
   triangle-GSS local sources and convert `m^-1 s^-1` to `cm^-1 s^-1` by
   `1e-2`.
3. Geometric-zero local edges must throw if any source is nonzero.
4. The candidate JSON must differ from the baseline only at
   `solver.mobility.high_field_driving_force`.
5. Sentaurus projected currents must never be labeled native directed edges.
6. Verifiers must recompute arithmetic from CSV rows instead of trusting
   manifest summaries.
7. No generated build output may be staged.

## Validation result

| Validation | Result |
|---|---:|
| focused Python tests | 24 passed |
| high-risk C++ tests | 7 executables, all passed |
| Release build | passed |
| full CTest | 469/469 passed |
| `ascii_sources` | passed |
| `git diff --check` | passed |
| Task 3/4/5/6/source/Phase C verifiers | passed |
| final Phase F A/B verifiers | passed |
| final impact A/B verifier | passed |

## Known non-blocking boundary

The DCSweep `sg_avalanche_source_integral_total` diagnostic is a global-edge
SG proxy when the production assembler uses
`triangle_gss_gradqf_truncated`. Final Phase F no longer consumes that field
for impact parity; it reads the fixed-state operator-audit triangle source.
A future schema revision should rename or split the DCSweep proxy rather than
silently changing the existing CSV meaning.

No production formula change is included in the Task 8/9 delta.
