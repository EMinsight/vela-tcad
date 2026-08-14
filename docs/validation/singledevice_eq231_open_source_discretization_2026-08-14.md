# SingleDevice Eq. 231 open-source discretization prototypes (2026-08-14)

## Scope

This audit implements and tests two opt-in alternatives to the qualified
`p1_direct` all-material density-gradient equation:

- `cvfem_full`: median-dual CVFEM integration of the complete expanded
  Laplacian, quadratic-gradient, and reaction terms, following the equation
  layout used by Charon-style control-volume finite elements.
- `oxide_boundary: devsim_wkb`: a DEVSIM/Garcia-Asenov-style oxide closure
  with an independent SiO2 density-gradient coefficient and a WKB penetration
  flux at transport/SiO2 interfaces.

Both remain experimental. The qualified default is unchanged.

## Unit oracle

The quasi-1D two-material MOSCAP strip test verifies:

- the analytic WKB penetration depth
  `hbar/sqrt(2*m_barrier*q*barrier)`; the DEVSIM example defaults
  `m_barrier=0.4*m0` and `barrier=3.15 V` give `0.1739 nm`;
- zero quantum potential without the oxide closure;
- a finite, symmetric interface response with the closure;
- the complete CVFEM operator preserves a zero manufactured state;
- direct API validation rejects unknown discretizations and non-positive WKB
  parameters;
- a small line-search step cannot report convergence while the nonlinear
  residual remains above tolerance.

Targeted results:

```text
test_density_gradient_quantum_potential: 62 assertions, 18 cases passed
test_newton_solver [density_gradient]:   38 assertions, 2 cases passed
```

After completing the Debug build, the broader CTest run passed 618 of 619
registered tests before the last concurrently scheduled contact-current
reporting test remained active for about 15 minutes. The run was stopped and
that exact remaining test was rerun alone; it passed all 16 assertions in
0.18 s. This points to a parallel CTest scheduling/process interaction rather
than a reproducible code failure.

## Fixed Sentaurus-state decomposition

The comparison uses the same imported linear and saturation endpoint states.
All residuals are evaluated before nonlinear iteration.

| endpoint | variant | max free residual | max node | stiffness L1 | quadratic L1 | reaction L1 | interface-source L1 |
|---|---|---:|---:|---:|---:|---:|---:|
| linear | `p1_direct` | 2080.199911 | 1848 | 116035.645 | 349648.854 | 242290.617 | 0 |
| linear | `cvfem_full` | 2080.199911 | 1848 | 116035.645 | 349648.854 | 242290.617 | 0 |
| linear | `cvfem_full + devsim_wkb` | 3851.482866 | 1848 | 116035.645 | 349648.854 | 102693.060 | 4730.663 |
| saturation | `p1_direct` | 2080.191965 | 1848 | 117383.049 | 351667.069 | 243923.146 | 0 |
| saturation | `cvfem_full` | 2080.191965 | 1848 | 117383.049 | 351667.069 | 243923.146 | 0 |
| saturation | `cvfem_full + devsim_wkb` | 3851.469037 | 1848 | 117383.049 | 351667.069 | 104546.292 | 4730.663 |

For affine Tri3 cells, the median-dual CVFEM flux and the Galerkin P1
stiffness give the same assembled fixed-state operator here. The quadratic
field is constant inside a Tri3 and both schemes assign one third of its cell
integral to each vertex. Consequently, replacing P1 by this complete CVFEM
form does not change the residual hotspot.

The WKB variant changes the SiO2 reaction coefficient and adds the interface
source, but increases the maximum fixed-state residual by about 85.15%. The
maximum remains node 1848 in the mirrored poly-reoxidation neighborhood. Its
direct interface-source contribution is zero, so the regression is caused by
the changed oxide volume equation rather than the face source acting directly
on that node.

## Endpoint oracle

The global quantum-potential inner solve must converge at both imported
endpoints before a complete self-consistent Id-Vg sweep is admitted.

| endpoint | variant | inner iterations | final inner residual | raw quantum change (V) | result |
|---|---|---:|---:|---:|---|
| linear | `p1_direct` | 500 | 1200.755171 | 1.389548 | fail |
| linear | `cvfem_full` | 500 | 1186.861301 | 1.412816 | fail |
| linear | `cvfem_full + devsim_wkb` | 114 | 139.473206 | 1.433381 | fail |
| saturation | `p1_direct` | 500 | 1167.735047 | 1.465372 | fail |
| saturation | `cvfem_full` | 500 | 1168.796687 | 1.429268 | fail |
| saturation | `cvfem_full + devsim_wkb` | 144 | 138.788982 | 1.433234 | fail |

The WKB runs previously stopped on a vanishing line-search update and were
incorrectly marked as converged despite residuals near 139. The global Eq. 231
solver now requires its residual tolerance as well as its update tolerance;
both endpoint runs correctly fail with
`electron_density_gradient_max_iterations`.

## Decision

Neither experimental path passes the two endpoint gates, so the full linear
and saturation Id-Vg curves are intentionally not run. The experiment rules
out ordinary Tri3 P1-versus-CVFEM volume integration as the missing mechanism,
and the DEVSIM WKB oxide closure is not a drop-in representation of the
Sentaurus all-material Eq. 231 interface contract for this structure.

The remaining implementation target is a material-side internal-interface
formulation: region-side equation traces or discontinuous auxiliary degrees of
freedom, with explicit Si/SiO2 and PolySilicon/SiO2 transmission conditions
and a consistent coupled Jacobian. That formulation should first reduce the
fixed-state node-1848 residual before either endpoint or a full Id-Vg curve is
rerun.

## Reproduction

Generate the three endpoint configurations with:

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass `
  -File scripts/make_singledevice_eq231_experimental_configs.ps1
```

Run a generated endpoint with `build-release/vela_example_runner.exe
--config <generated-json>`. Diagnostic CSV and summary files are written below
`reports/eq231_experimental_20260814/{p1,cvfem,cvfem_wkb}` in the generated
SingleDevice workspace.
