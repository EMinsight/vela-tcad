# Task 4 report: Recover potential gradients and electric-field vectors

## Status

Implemented the diagnostic-only Task 4 field inverse layer and committed exactly
the two requested Task 4 source/test files. No `include/` or `src/` file changed
relative to phase base `a5524cf`, no dependency was added, no remote run was
started, and the existing untracked `docs/validation/figures/` directory was
preserved.

## Commit

`1010e37 Add Minimal6 field inverse candidates`

Committed files:

- `scripts/pn2d_minimal6_diagnostics/inverse_fields.py`
- `tests/regression/test_pn2d_minimal6_inverse_fields.py`

The report itself was written after the commit and intentionally remains outside
that Task 4 code commit.

## TDD evidence

### Initial RED

The focused test was written before the implementation and run with:

```powershell
C:\Users\qzw\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe -m unittest tests.regression.test_pn2d_minimal6_inverse_fields -v
```

It failed for the intended reason:

```text
ModuleNotFoundError: No module named 'scripts.pn2d_minimal6_diagnostics.inverse_fields'
Ran 1 test in 0.000s
FAILED (errors=1)
```

### Initial GREEN

After the minimal implementation, the same command passed:

```text
Ran 5 tests in 0.003s
OK
```

### Review-fix RED and GREEN

Self-review found that the requirement says relative magnitude error is computed
only when the reference *exceeds* the declared floor. A boundary test at exactly
the floor was added first. RED produced one intended failure because the sample
was classified `VALID` rather than `BELOW_FLOOR`. The comparison was then changed
from `< floor` to `<= floor`; the focused suite returned to 5/5 GREEN.

## Fresh scoped verification

Only the requested Task 1 contract, Task 2 input, and Task 4 field modules were
run together:

```powershell
C:\Users\qzw\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe -m unittest tests.regression.test_pn2d_minimal6_inverse_contracts tests.regression.test_pn2d_minimal6_inverse_inputs tests.regression.test_pn2d_minimal6_inverse_fields -v
```

Result:

```text
Ran 24 tests in 16.126s
OK
```

Before commit, `git diff --cached --check` exited zero, the cached path list was
exactly the two requested Task 4 files, and
`git diff --name-only a5524cf -- include src` was empty.

## Delivered behavior

- `triangle_gradient()` implements the specified determinant formula for exact
  P1 gradients and rejects singular or non-finite triangles.
- `cell_to_node_vectors()` and `cell_to_edge_vectors()` area-weight adjacent cell
  vectors while retaining deterministic support identities and explicit normalized
  cell weights.
- `edge_scalar_difference()` uses the declared start-to-end orientation and
  rejects zero-length/non-finite edges.
- `mirror_vector()` applies the global Cartesian x-mirror convention
  `(Ex, Ey) -> (-Ex, Ey)`.
- `vector_error()` separates magnitude and direction status. It uses vector
  magnitude relative error only above the caller-declared reference floor, clamps
  the dot-product cosine before `acos`, and types geometric zero, below-floor,
  missing, non-finite, and undefined direction without converting any of them to
  ordinary zero error.
- `evaluate_field_candidates()` consumes canonical SI node observations plus
  explicit topology connectivity and emits the four fixed candidates:
  `triangle_minus_grad_psi`, `node_area_weighted_minus_grad_psi`,
  `edge_area_weighted_minus_grad_psi`, and
  `signed_edge_minus_delta_psi_over_h`.
- Candidate/reference comparisons use declared transforms on identical support:
  node electric field is explicitly averaged to cell or edge support, and the
  signed edge candidate is compared with the directed tangent projection.
- Classification uses only the immutable contract thresholds: median relative
  magnitude error `<= 2%` and median direction error `<= 1 degree`. No threshold
  or scale depends on bias, node, edge, cell, topology, discovery, or holdout.
- State, candidate, node, edge, cell, and adjacent-cell ordering is deterministic.

## Numerical and scientific self-review

- Triangle orientation is permitted in either winding because the gradient
  determinant retains its sign while reconstruction weights use absolute area.
- Determinants with absolute value `<= 1e-300`, zero-length edges, non-finite
  geometry, and non-finite derived gradients fail closed.
- `E = -grad(psi)` is explicit and covered on `f = 3*x - 4*y + 2`, yielding
  `grad(f) = (3, -4)` and `E = (-3, 4)`.
- Reversing a directed edge reverses its scalar difference sign; reconstruction
  of a physical tangent vector remains orientation-consistent.
- Reference geometric zeros and values at or below the floor are absent from
  relative-error statistics; zero candidate direction and zero reference
  direction are `DIRECTION_UNDEFINED`.
- Missing observations remain `None` with a typed status; incompatible topology

## Determinant overflow review fix (2026-07-21)

RED: the focused overflow-triangle regression failed because no `ValueError` was raised.
GREEN: after the finite determinant guard, the Task 4 suite passed 7/7 in 0.005 s.
Formatting cleanup commit `8019164 Normalize Task 4 determinant fix line endings`
normalizes the two fallback-patched lines; `git diff --check 3945a50..HEAD` is clean.
  or coordinate frames fail closed rather than being compared silently.

## Concerns

No blocking concern remains for Task 4. The evaluator intentionally requires
canonical SI observations and explicit directed edge/triangle topology; future
tasks must preserve those contracts rather than passing raw-unit or mixed-support
data. Task 5 was not started.

## Review-fix continuation (2026-07-21)

The review-fix work resumed at `1010e37` with only the requested Task 4 source
and test files modified, alongside unstaged Task 3/4 reports and the pre-existing
untracked `docs/validation/figures/` directory. The focused field suite was run
first and was already GREEN: 7 tests passed in 0.009 s. No new correction was
made, so no new RED result is claimed; the original RED/GREEN evidence above
remains the applicable test-first record for this behavior.

Debug/self-review confirmed that overflowing finite endpoint subtraction and
length are rejected before use; finite extreme vector directions are formed from
normalized components rather than an overflow-prone raw dot product; every
pre-existing non-`VALID` observation status is propagated unchanged; and only
the canonical field gates (`0.02`, `1.0`) are accepted. The adversarial tests
cover statuses, overflowing edges, huge and subnormal vectors, both exact gate
boundaries, just-outside gates, and a caller threshold override.

Fresh scoped GREEN verification used only Tasks 1, 2, and 4:

```powershell
C:\Users\qzw\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe -m unittest tests.regression.test_pn2d_minimal6_inverse_contracts tests.regression.test_pn2d_minimal6_inverse_inputs tests.regression.test_pn2d_minimal6_inverse_fields -v
```

Result: `Ran 26 tests in 24.674s` and `OK`.

Final review-fix commit: `3945a50 Harden Task 4 inverse field review fixes`.
It contains exactly `scripts/pn2d_minimal6_diagnostics/inverse_fields.py` and
`tests/regression/test_pn2d_minimal6_inverse_fields.py`; the report remains
