# Sentaurus 2022 TransportModels DD/DG validation

Validation window: 2026-08-18 to 2026-08-20

Status: **complete with stated applicability limits**.  Both DD and DG now
contain exact 21-point Id-Vg and Id-Vd comparison lattices.  The baseline is
qualified for nonzero-bias Id-Vd and on-state Id-Vg; deep-off Id-Vg remains
below Vela's present nonlinear residual/current floor.

## Scope

This validation freezes the Sentaurus T-2022.03-SP2 `TransportModels`
MOSFET reference and establishes a controlled Vela comparison.  DD and DG use
the same imported TDR mesh, doping, materials, contacts, bias lattices, and
classical transport models.  The only DD-to-DG physics delta is the enabled
electron quantum-potential equation.

The reusable entry points are:

- `scripts/run_transportmodels_dd_dg_workflow.py` for materialization,
  execution, restart provenance, and exact-lattice trend comparison;
- `scripts/analyze_transportmodels_workflow.py` for bias-regime error
  analysis.

Generated evidence is under
`build-release/reference_tcad/transportmodels_sentaurus2022/vela_baseline`.

## Frozen Sentaurus evidence

The official T-2022.03-SP2 run contains four 21-point curves:

| Branch | Curve | Bias range | Final current (A/um) |
|---|---|---:|---:|
| DD | Id-Vg | -1.0 to 2.2 V | 1.81040907111e-3 |
| DG | Id-Vg | -1.0 to 2.2 V | 1.68607324102e-3 |
| DD | Id-Vd | 0.0 to 2.0 V | 8.08310849704e-4 |
| DG | Id-Vd | 0.0 to 2.0 V | 7.05525753105e-4 |

The neutral Vela TDR export contains 3315 nodes, 6456 triangles, seven
geometry regions, and four contacts.  The node count is larger than the
Sentaurus geometry inventory because interface-side nodes are split during
neutral export.

## Numerical path decision

The first DD attempt used the legacy edge projection of the quasi-Fermi
gradient.  It reached Vd=0.84 V but switched to a discontinuous current state
and failed at 0.8403125 V with `line_search_non_decrease`; the dominant
electron residual was in the drain-junction region.

The production candidate uses
`mobility.high_field_gradient_discretization = transport_cell_vector`.
This is the two-dimensional cell-vector interpretation of Sentaurus
`GradQuasiFermi`.  A separate electric-field-drive control also converged but
changed off-state current by many orders of magnitude, so it was not promoted
as the Sentaurus semantic match.

The Id-Vg dependency graph contains an explicit same-bias relaxation after
ramping Vd to 1.1 V.  Without this stage, the ramp endpoint carries a residual
quasi-Fermi gradient and reports about 5.95e-11 A/um at Vg=-1 V.  Re-solving
the identical bias reaches the numerical residual floor and reports about
1.58e-24 A/um.  The relaxed state, rather than the ramp endpoint, is therefore
the Id-Vg comparison seed.

## DD result

Workflow:
`build-release/reference_tcad/transportmodels_sentaurus2022/vela_baseline/workflow_dd_vector_run01`

All six stages pass: Id-Vg equilibrium, drain ramp, final-bias relaxation,
Id-Vg curve, independent Id-Vd equilibrium, and Id-Vd curve.  Both comparison
lattices contain all 21 exact Sentaurus biases and both trend checks pass.

| Curve / regime | Points | Median absolute log error (dex) | Maximum absolute log error (dex) | Maximum relative error |
|---|---:|---:|---:|---:|
| Id-Vg off, Vg <= -0.68 V | 3 | 8.82953 | 9.01373 | approximately 1.0 |
| Id-Vg transition, -0.52 <= Vg <= 0.12 V | 5 | 0.189391 | 0.207659 | 0.613091 |
| Id-Vg on, Vg >= 0.28 V | 13 | 0.0074512 | 0.0737019 | 0.184955 |
| Id-Vd, Vd > 0 V | 20 | 0.00310681 | 0.00374633 | 0.00866356 |

Endpoint agreement is 2.096% at Id-Vg Vg=2.2 V and 0.702% at Id-Vd
Vd=2.0 V.

The generic `comparison_status=pass` means exact lattice and monotonic trend
passed.  It is not an amplitude-accuracy waiver.  The qualified baseline is:

- qualified for nonzero-bias Id-Vd and on-state Id-Vg comparison;
- usable with an explicit 0.19-dex-scale model discrepancy in the transition
  region;
- not qualified for the three deep-off Id-Vg points, where the requested
  current is below the present nonlinear residual/current numerical floor.

## DG result

The DG workflow uses the identical DD configuration plus the electron
quantum-potential equation.  Its exact 21-point Id-Vg and Id-Vd lattices are
complete, and both monotonic-trend comparisons pass.

| Curve / regime | Points | Median absolute log error (dex) | Maximum absolute log error (dex) | Maximum relative error |
|---|---:|---:|---:|---:|
| Id-Vg off, Vg <= -0.68 V | 3 | 8.39509 | 8.51596 | 1.00000 |
| Id-Vg transition, -0.52 <= Vg <= 0.12 V | 5 | 0.131935 | 0.339151 | 0.542017 |
| Id-Vg on, Vg >= 0.28 V | 13 | 0.0366998 | 0.0820144 | 0.207854 |
| Id-Vd, Vd > 0 V | 20 | 0.0418032 | 0.0442337 | 0.107219 |

Endpoint agreement is 4.016% at Id-Vg Vg=2.2 V and 9.842% at Id-Vd
Vd=2.0 V.  At the Id-Vd endpoint, Sentaurus gives
7.05525753105e-4 A/um and Vela gives 7.74966247896e-4 A/um.

The DG qualified baseline is:

- qualified for the complete nonzero-bias Id-Vd sweep, whose maximum
  relative error is 10.72%;
- qualified for on-state Id-Vg shape and amplitude comparison, with at most
  0.0820 dex / 20.79% error over the defined on-state region;
- usable with an explicit 0.339-dex / 54.20% maximum discrepancy in the
  transition region;
- not qualified for the three deep-off Id-Vg points, where the Sentaurus
  current lies below Vela's present nonlinear residual/current floor.

### DG numerical convergence policy

The initial 20-outer-iteration workflow exhausted the electron
density-gradient outer fixed-point budget at Vd=0.3484375 V, while every
density-gradient inner solve converged.  Vector Aitken acceleration with
bounded relaxation completed most of the bias path.  The high-drain tail
required a provenance-preserving restart and a more conservative fixed
relaxation policy.

The final Vd=2.0 V point used fixed relaxation 1.5 and an outer budget of 200.
It converged at outer iteration 87 with a final raw quantum-potential change
of 3.728541065e-8 V.  This is a validation-tail policy, not a proposed global
Vela default; the equation, interface boundary condition, and convergence
tolerance were unchanged.

The final drain-tail workflow is
`workflow_dg_idvd_resume_1p9_direct2p0_fixed1p5_outer200_run01`.  Its manifest
and comparison status are both `pass`.

## Error definitions and method

All comparisons use the exact Sentaurus bias lattices; no interpolated point
is included in the regional statistics.  The absolute log error is
`abs(log10(abs(I_Vela) / abs(I_Sentaurus)))`.  Relative error is
`abs(I_Vela - I_Sentaurus) / abs(I_Sentaurus)`.  Id-Vd relative statistics
exclude Vd=0 V because both currents are zero.  Currents are total drain
current normalized by device width and reported in A/um.

The generic `comparison_status=pass` asserts lattice completeness and trend
agreement.  Amplitude qualification is governed by the regional statistics
and limits stated above.

## Final figures

![TransportModels device mesh](../progress_report_2026Q3/2026-08-19_transportmodels_dg_daily_report/figures/transportmodels_device_mesh.png)

![DG Id-Vg comparison](../progress_report_2026Q3/2026-08-19_transportmodels_dg_daily_report/figures/transportmodels_dg_idvg_comparison.png)

![DG Id-Vd comparison](../progress_report_2026Q3/2026-08-19_transportmodels_dg_daily_report/figures/transportmodels_dg_idvd_comparison.png)

The same directory also contains SVG versions and
`transportmodels_dg_figure_summary.json`, which records the chart-to-source
mapping and plotted summary metrics.

## Final evidence ledger

| Artifact | SHA-256 |
|---|---|
| DG 21-point Id-Vg candidate | `412DBE55FA53E774F82D44F214DD22AC8CD5B952AA36EEDE4248937E932C768D` |
| DG 21-point Id-Vd candidate | `FD1B5E7A4E9480ADCA8283797C4B228A4EABF34382AABC05E244A90DDEB77FD5` |
| DG Vd=2.0 V final state | `BDB95E6C7CC3FF538D3AE452F20C2172CEA4780DF06BF9065258ED922C4F1B53` |
| DG bias-regime analysis JSON | `748AA319E92A1AAA8C4D3D8E19A9ED067B15398499075352732098F2F7F913B1` |
| DG final workflow manifest | `14709AEDD4F6E9648454AC77EEEB585E98ABCD67E642F2EB5C9705C07542969B` |

## Validation conclusion

The TransportModels DD/DG comparison objective is complete.  Confidence is
**high** for evidence completeness, bias alignment, provenance, and the
reported calculations, and **medium** for cross-simulator amplitude agreement
because transition/deep-off Id-Vg retains an acknowledged model/numerical
discrepancy.  The result is ready to share with those caveats.  Within the
narrowed task scope, the next maintenance step is to freeze these candidates
as regression inputs and enforce separate on-state, transition, and Id-Vd
thresholds rather than a single whole-curve amplitude threshold.

## Final verification

An independent CSV audit confirmed 21 finite, unique, bias-aligned points in
each DG curve.  The largest reference/candidate bias-coordinate difference is
2.22e-16 V for Id-Vg and 1.11e-16 V for Id-Vd.  Recalculation independently
reproduced all regional values in the DG table, including the 10.7040% Id-Vd
maximum, 10.1040% Id-Vd median, and 9.8424% endpoint relative errors.

The focused Python regression suite reports `Ran 9 tests in 1.267s` and `OK`.
It covers split Id-Vg/Id-Vd candidate provenance, external drain-prefix
validation, checkpoint bridge biases, exact restart lattices, controlled DD/DG
physics deltas, regional analysis, and zero-bias exclusion.  Python byte-code
compilation also passes for the workflow, analysis, and plotting scripts.

Two focused C++ density-gradient/Newton tests pass: frozen quantum-potential
transport consistency and parsing of Sentaurus electron density-gradient
controls.  All PNGs were visually inspected after generation; labels, legends,
error panels, and mesh detail are visible without clipping.

## Historical execution ledger

The following checkpoints document earlier requested pauses.  They are
retained for provenance and are superseded by the final completion above.

### 2026-08-18 22:00 closeout

Closeout was recorded at `2026-08-18 22:00:03 +08:00` (Asia/Shanghai).
Execution is paused at this checkpoint; no later simulation stage was started.

No `vela_example_runner` process remains.  The most important generated
evidence hashes are:

| Artifact | SHA-256 |
|---|---|
| DD workflow manifest | `3CB75B60CEC7C683FEDAF05028A8B5C12765AA98CABDA089E9C26E49DE5F22DB` |
| DD bias-regime analysis | `47DC3E5B9B25629C6D515CEA5A73F8B171A9727829B6A0BACEA316A4921B9A64` |
| DG initial failed-workflow manifest | `205CC6798A0645028BDD72C79473FA19D59E76EB8EEDDC7C43F90F5568B714E9` |
| DG 0.30 V checkpoint | `339D0D4EDB4A415E4F1476F1814554E67466C4AB3864A8BF149F286B9242CD47` |
| DG 0.35 V checkpoint | `D2B19A5753F43B5E313E823082477136CEDD3A2576D903B3A247AEF38459C330` |
| DG 0.45 V checkpoint | `5DBE0C46F62DE4998D7345013C12B046F081D23965B01D8BA1E2B86369765B4E` |
| DG 0.50 V checkpoint | `60177597126F49FC041B9185F06BDF5DD2E897E86A6DAAB7E513DC34D4B2CFEB` |
| DG 1.10 V final checkpoint | `D174C4491C1C44596BC0BF538FDC960E8A0751C312E0B568639FF3B93510167F` |
| DG 0.55--1.10 V curve segment | `705DC6DE576C661D5672417B6D906BBBAFCC518180541C7A343CD6119EA11C55` |

The next resume point at that time was the DG 1.10 V final checkpoint.

### Historical verification at 2026-08-18 closeout

The focused regression command is:

```powershell
D:\msys64\ucrt64\bin\python.exe -m unittest `
  tests.regression.test_analyze_transportmodels_workflow `
  tests.regression.test_transportmodels_dd_dg_workflow -v
```

It covers the DD/DG controlled-delta contract, the six-stage restart graph,
the 2-D quasi-Fermi discretization, exact reference lattices, bias-regime
partitioning, and zero-bias exclusion for Id-Vd amplitude summaries.

Historical closeout result: `Ran 4 tests in 0.254s` and `OK`.

### 2026-08-18 22:10 requested pause

Execution briefly resumed from the frozen 1.10 V drain-ramp checkpoint using
the new provenance-preserving `--idvg-ramp-state dg=PATH` workflow entry.
The DG same-bias relaxation at Vg=-1.0 V and Vd=1.1 V completed successfully.
Its final-state SHA-256 is
`8050FD0DD0F7D0DEF9038ABBC3A437A4C7ECBECCB5BA98A129932FD7C28B515D`.

At the user's request, the first Id-Vg point (Vg=-0.84 V) was interrupted at
`2026-08-18 22:10:24 +08:00` while its density-gradient outer iteration was
still contracting.  It had not emitted a converged bias row or point state,
so it is not counted as completed evidence.  No `vela_example_runner` process
remains.  The next run must restart the Id-Vg curve from the completed
same-bias relaxation state in
`workflow_dg_aitken_resume_run01`; the workflow's `--resume` mode will reuse
the hashed relaxation stage and recompute the interrupted point.

After adding regression coverage for the external checkpoint provenance,
the focused suite reports `Ran 5 tests in 0.465s` and `OK`.

### 2026-08-19 requested pause

Execution was paused at `2026-08-19 15:35:34 +08:00` after the user
requested an immediate stop.  No `vela_example_runner` process remains.

The exact 21-point DG Id-Vg lattice is complete.  Its final-state SHA-256 is
`6992F8CC3B235A8BB2F92C264DE087E74C9BA192AA511B402E047B391A65CF5B`.
The DG Id-Vd equilibrium stage is also complete.  The Id-Vd curve has 14
completed checkpoints on the exact 0.1 V lattice, from 0.1 V through 1.4 V.
The 1.4 V checkpoint SHA-256 is
`3A968B07A1A6EE5ED31F058EFC98F51992D1756A929C1550C16B98A7735515BE`.

The interrupted 1.5 V solve had not emitted a point checkpoint and is not
counted as completed evidence.  Resume from the 1.4 V checkpoint and run only
the remaining exact biases 1.5--2.0 V.  The frozen outer-80 Id-Vd config has
SHA-256
`426C2EC4A8E0BFB361730D32B2AC8E8ED08A8371E9A7EBC7143E75F5220410E5`.
Those remaining points, merged candidate, analysis, figures, and report were
completed on 2026-08-20 as recorded above.  Work outside the narrowed
TransportModels DD/DG scope is intentionally not started here.
