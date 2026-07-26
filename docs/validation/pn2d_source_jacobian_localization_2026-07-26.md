# PN2D source-only Jacobian localization

Date: 2026-07-26

Plan task: Task 12 of
`docs/superpowers/plans/2026-07-26-pn2d-high-bias-process-variable-jacobian-localization.md`

## Outcome

Typed outcome: `incomplete_analytic_derivative`

Nonlinear authorization: stopped. Tasks 15-19 remain unauthorized. Task 13
may proceed as the plan's independent native-observation evidence task.

The first unresolved dependency is not a missing Van Overstraeten scale or a
field, mobility, current, source, geometry, or voltage scale. The opt-in
element-edge source block in `CoupledDDAssembler::assembleJacobian` is a
cell-local central finite-difference surrogate for the complete source chain,
not an analytic derivative decomposed by dependency. It therefore cannot
satisfy the plan's independent analytic-versus-FD contribution and convergence
contract.

## Source chain

For each carrier the evaluated chain is:

`scaled state -> physical psi/QFP -> density -> low-field mobility ->`
`QFP-gradient high-field mobility -> SG signed edge flux ->`
`GSS/Laux cell current vector -> vector norm ->`
`Van Overstraeten alpha(norm(P1 QFP gradient)) ->`
`element-vertex box source -> node accumulation -> continuity scaling`.

The assembled source value and diagnostic replay use the same production
source evaluator. The Jacobian path perturbs `psi`, electron QFP, and hole QFP
at each cell vertex and central-differences the complete cell source. It does
not separately emit or sum analytic contributions for carrier statistics, SG
Bernoulli terms, mobility, current reconstruction, current norm, alpha,
mapping, accumulation, and scaling.

## Accepted configuration fix

The cell-local physical perturbation previously used

`1e-7 * max(1 V, abs(physical potential))`.

The assembler-level audit perturbs scaled coordinates, which corresponds to

`1e-7 * max(V0, abs(physical potential))`.

The opt-in source derivative now uses the latter rule. A focused contract test
covers both a near-zero physical potential and a 20 V potential. No production
default, model selector, Van Overstraeten parameter, QFP-gradient default, or
source scale changed.

## Frozen RED and post-fix evidence

The pre-fix frozen maximum source-only true-relative error was
`0.9640948767506723`.

Two post-fix raw roots were generated:

- `build-release/pn2d-task12-jacobian-stepfix-20260726-a`
- `build-release/pn2d-task12-jacobian-stepfix-20260726-b`

All 18 corresponding CSV files are byte-identical between the two roots.
Selected source-only rows at the frozen `1e-7` relative step are:

| Topology/bias/variant | Analytic norm | FD norm | Difference norm | True relative | psi relative | electron-QFP relative | hole-QFP relative |
|---|---:|---:|---:|---:|---:|---:|---:|
| coarse7x3 -10 V opt-in | 1.152481682e-9 | 1.152426816e-9 | 4.168605835e-13 | 3.617069060e-4 | 8.742555409e-10 | 1.272764689e-3 | 1.846765135e-4 |
| coarse7x3 -20 V opt-in | 3.805195125e-12 | 3.805188059e-12 | 2.613262474e-17 | 6.867617528e-6 | 1.489024279e-9 | 1.717137337e-6 | 2.004342545e-5 |
| Minimal6 mirror -20 V opt-in | 1.825509338e-12 | 1.825509337e-12 | 1.998767342e-21 | 1.094909404e-9 | 7.392520538e-9 | 1.103991174e-9 | 8.183515767e-12 |
| coarse7x3 -10 V production | 6.066834044e-2 | 6.066834043e-2 | 2.839060006e-11 | 4.679640131e-10 | 6.253658090e-11 | 6.715030822e-10 | 7.812244306e-10 |
| coarse7x3 -20 V production | 3.030023688e-6 | 3.030023689e-6 | 2.076511948e-15 | 6.853121165e-10 | 8.144603962e-10 | 4.968320372e-10 | 9.026384474e-10 |

The scaled-step fix closes the former Minimal6 -20 V failure at the selected
step. It does not close the nonzero coarse opt-in gate. The remaining coarse
difference is localized to the electron-QFP and hole-QFP columns; the `psi`
column passes.

Rows whose analytic and FD norms are at most `1e-12` satisfy the frozen
near-zero absolute gate. They are not relabeled as nonzero successes.

## Independent step convergence

The outer source-only FD was repeated at `1e-6`, `3e-7`, `1e-7`, `3e-8`, and
`1e-8`. The cell-local surrogate always retains its internal `1e-7`
perturbation.

| Case | Outer step | Analytic norm | FD norm | True relative |
|---|---:|---:|---:|---:|
| coarse7x3 -10 V opt-in | 1e-6 | 1.152481682e-9 | 1.141333828e-8 | 8.991219432e-1 |
| coarse7x3 -10 V opt-in | 3e-7 | 1.152481682e-9 | 3.423709014e-9 | 6.638360851e-1 |
| coarse7x3 -10 V opt-in | 1e-7 | 1.152481682e-9 | 1.152426816e-9 | 3.617069060e-4 |
| coarse7x3 -10 V opt-in | 3e-8 | 1.152481682e-9 | 4.067717052e-10 | 6.800738919e-1 |
| coarse7x3 -10 V opt-in | 1e-8 | 1.152481682e-9 | 5.364905210e-10 | 8.971865104e-1 |
| Minimal6 mirror -20 V opt-in | 1e-6 | 1.825509338e-12 | 1.331535237e-11 | 8.630926160e-1 |
| Minimal6 mirror -20 V opt-in | 3e-7 | 1.825509338e-12 | 4.217333003e-12 | 5.677335676e-1 |
| Minimal6 mirror -20 V opt-in | 1e-7 | 1.825509338e-12 | 1.825509337e-12 | 1.094909404e-9 |
| Minimal6 mirror -20 V opt-in | 3e-8 | 1.825509338e-12 | 2.073672496e-12 | 1.499825491e-1 |
| Minimal6 mirror -20 V opt-in | 1e-8 | 1.825509338e-12 | 4.743564632e-12 | 6.200868450e-1 |

There is no independent convergence plateau. Agreement at exactly `1e-7` is
expected because both paths then use the same perturbation policy; it is not
evidence for an analytic derivative.

The frozen Minimal6 base state is not at a zero-current or Van Overstraeten
switch branch. Its element current-vector magnitudes are approximately
`2.4e14` to `1.5e15 m^-2 s^-1`, and its QFP-gradient impact fields are
approximately `9.6e6` to `1.0e7 V/m`, below and separated from the
`4.0e7 V/m` Van Overstraeten parameter switch. The failed convergence is
therefore classified as an incomplete analytic derivative plus FD
resolution/truncation sensitivity, not as a proven physical branch
nonsmoothness.

## Exit gates

| Gate | Result |
|---|---|
| Nonzero analytic/FD true-relative <= 1e-8 | FAIL: 3.617069060e-4 on coarse7x3 -10 V opt-in |
| Near-zero absolute difference <= 1e-12 | PASS |
| Analytic dependency contribution sum <= 1e-12 | NOT AVAILABLE: no analytic dependency decomposition |
| Independent FD contribution sum <= 1e-12 | NOT AVAILABLE: no dependency-level FD decomposition |
| Diagnostic source equals residual source <= 1e-12 | same evaluator is used, but no independent hashed emission was produced |
| Residual/Jacobian configuration hashes identical | NOT AVAILABLE in the current C++ probe output |
| No unrelated/default formula changes | PASS |

The unavailable contribution and hash gates are not inferred from code-path
proximity. They remain unproven.

## Verification

- Scaled perturbation contract: 2 assertions passed.
- Existing element-edge source Jacobian focused test: 8 assertions passed.
- Newton solver and diagnostic tests: 728 assertions in 63 cases passed.
- A/B post-fix Jacobian CSVs: 18/18 byte-identical.

The next acceptable source-Jacobian implementation must provide a genuinely
analytic or otherwise independently convergent derivative for every dependency
stage, together with contribution closure and emitted configuration hashes.
No tolerance relaxation is authorized.
