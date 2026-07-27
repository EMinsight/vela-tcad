# PN2D element-edge local forward AD source Jacobian

Date: 2026-07-27

Scope: opt-in `element_edge_sg_gss_laux` on Tri3 cells only.

## Outcome

Typed outcome: `jacobian_gate_failed`

The requested local forward-AD implementation is complete and its direct
cell-source derivative gate passes. Nonlinear PN2D authorization remains
stopped because the frozen assembled real-state source-block probe does not
meet the Task 12 `1e-8` true-relative gate on coarse7x3 `-20 V`. Dependency
contribution closure and residual/Jacobian configuration hashes are not yet
emitted by the probe.

No production selector, Van Overstraeten parameter, QFP-gradient default,
mobility default, geometry scale, or source scale changed.

## RED

The focused assembled test was changed from a nearly constant-alpha,
constant-mobility manufactured state to:

- native Van Overstraeten parameters;
- `quasi_fermi_gradient` avalanche drive;
- QFP-gradient high-field mobility; and
- a nonzero approximately 5--20 MV/m manufactured field.

With the previous cell-local `1e-7` central-FD Jacobian surrogate and an
independent outer `1e-6` central difference, the source-only true-relative
error was:

`1.0468418106e-7`

This exceeded the `1e-8` gate and established RED with nonzero source support.

## Minimal opt-in patch

The opt-in source now uses one scalar-templated local residual chain:

`statistics -> SG current -> high-field mobility ->`
`GSS/Laux cell current vector -> current norm ->`
`Van Overstraeten alpha -> vertex source`

For value assembly the scalar is `Real`. For Jacobian assembly it is a fixed
dual number with nine derivative slots:

`3 Tri3 vertices x (psi, electron QFP, hole QFP)`.

The nine local derivatives are inserted into both carrier-continuity rows
after the existing continuity-unit scaling. The prior internal finite
difference remains only on the separate triangle-GSS cell-local path.

The variable-intrinsic-density SG expression uses an algebraically equivalent
`expm1(delta-QFP / Vt)` factorization when:

- `abs(delta-QFP / Vt) < 50`; and
- neither endpoint Boltzmann exponent is clamped.

Otherwise it uses the previous two-term expression. Exact flat QFP retains the
production zero value, while the dual path evaluates the limiting derivative
instead of discarding it at the value short-circuit.

## Direct local-source verification

The direct test uses a Tri3 cell with:

- Van Overstraeten;
- Masetti high-field mobility driven by QFP gradient;
- an exact flat-QFP edge;
- a nonzero cell QFP-gradient drive through the third vertex; and
- identical Real and dual evaluations of the shared source function.

| Gate | Result |
|---|---:|
| Real versus dual source value, per vertex | `<= 1e-13`, pass |
| Outer FD `1e-6` true-relative | `2.6048566474e-7` |
| Outer FD `3e-7` true-relative | `2.344558440e-8` |
| Outer FD `1e-7` true-relative | `2.607230050e-9`, pass |

The error decreases across all three symmetric steps and the fine-step result
passes `1e-8`.

## Frozen assembled real-state evidence

The updated runner was applied without overwriting historical evidence. New
generated outputs are under:

`build-release/pn2d-task12-local-ad-20260727`

They are not committed.

| State | Outer step | Analytic norm | FD norm | Difference norm | True relative |
|---|---:|---:|---:|---:|---:|
| coarse7x3 `-10 V` | `1e-6` | `8.852681154e-15` | `1.141237211e-8` | `1.141237213e-8` | `1.000000002` |
| coarse7x3 `-10 V` | `1e-7` | `8.852681154e-15` | `1.141185591e-9` | `1.141185614e-9` | `1.000000020` |
| coarse7x3 `-10 V` | `1e-8` | `8.852681154e-15` | `1.135905002e-10` | `1.135905234e-10` | `1.000000204` |
| coarse7x3 `-20 V` | `1e-7` | `3.668330531e-12` | `3.805188302e-12` | `1.058459256e-12` | `2.781621229e-1` |
| coarse7x3 `-20 V` | `1e-8` | `3.668330531e-12` | `3.667519931e-12` | `3.889023939e-14` | `1.060161811e-2` |
| Minimal6 mirror `-20 V` | `1e-7` | `3.720579100e-16` | `1.324320536e-12` | `1.324324921e-12` | `1.000003311` |
| Minimal6 mirror `-20 V` | `1e-8` | `3.720579100e-16` | `1.324404892e-13` | `1.324444039e-13` | `1.000029558` |

The coarse `-10 V` and Minimal6 rows are analytic-near-zero. Their assembled
FD norms decrease approximately linearly with the perturbation and therefore
do not define a nonzero convergence plateau. Minimal6 reaches the `1e-12`
near-zero absolute gate at `1e-8`; coarse `-10 V` does not. Coarse `-20 V` is
nonzero and approaches the AD norm, but its best recorded true-relative error
is `1.06e-2`, not `1e-8`.

The direct shared-source test proves the local derivative implementation. It
does not justify relabeling the assembled real-state block as closed.

The assembler-use regression independently scatters the nine local derivatives
into the global source matrix. Its relative mismatch is `3.4e-16` without
scaling and `3.3e-16` with scaling enabled; the manufactured assembled FD
comparison is `7.084e-10`. These gates prove local-AD assembly and coordinate
scaling, but they do not replace dependency-level contribution emission.

A remaining diagnostic limitation is that record-level signed edge flux,
current vector, and alpha fields still come from the legacy observation path,
while electron, hole, and combined source integrals now all come from the
shared source core. The three source components close, but the intermediate
record is not yet a complete shared-core dependency trace. Configuration hashes
also remain unavailable.

## Focused verification

- element-edge GSS/Laux tests: 89 assertions in 7 cases passed;
- impact-ionization tests: 523 assertions in 42 cases passed;
- Newton/diagnostic tests: 728 assertions in 63 cases passed;
- `ascii_sources`: passed;
- `git diff --check`: passed before report creation.

## Stop condition

Do not authorize coarse self-consistent or fine PN2D nonlinear sweeps from
this patch. The next acceptable work is to make the real-state probe evaluate
the shared source residual directly, emit dependency contribution closure and
configuration hashes, and then repeat the frozen coarse gate without
subtracting numerically dominant full-continuity blocks.
