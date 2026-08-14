# SingleDevice Eq. 231 material/interface follow-up (2026-08-13)

## Outcome

The saturation failure is not a drift-diffusion or terminal quasi-Fermi-state
failure.  With the imported Sentaurus electron quantum potential frozen at
`Vg=2.2 V, Vd=1.1 V`, Vela converges in seven Newton iterations and reports
`1.70224338069e-3 A/um`, versus the Sentaurus reference
`1.7069265949e-3 A/um` (0.274% relative error).

The potential-form implementation now represents the material terms omitted by
the earlier `psi-phin` proxy:

- the O-2018.06 Silicon electron parameters `gamma=3.6`, `theta=0.5`,
  `xi=eta=1` and 300 K DOS mass `m/m0=1.0618016171622988`;
- the conduction-band drive from material affinity;
- the DOS-mass driving potential and a material-resolved `gamma/m` coefficient;
- material-row `gamma/m` coefficients while keeping the Eq. 231 normal
  interface flux unweighted by `gamma/m`;
- an opt-in global-insulator mode using the O-2018.06 SiO2/Si3N4 values
  `gamma=1`, `xi=eta=0`, and `m/m0=0.42`.

## Interface qualification

The global-insulator mode is deliberately not the default. Sentaurus Eq. 231
enforces continuity of the potential-like quantity
`Ec + Phi_m + Lambda` at internal interfaces. Vela now represents that
quantity as the continuous primary unknown and reconstructs `Lambda` from the
local material drive.

The qualified default still keeps the semiconductor-side homogeneous Neumann
interface condition until the all-material endpoint oracle passes.

At the saturation oracle, the default semiconductor-side model still fails its
inner solve after 500 iterations and moves by `24.803958 V`.  This isolates the
next implementation item to the inhomogeneous semiconductor-to-unsolved-region
step boundary in Sentaurus Eq. 233.  The current homogeneous Neumann treatment
is sufficient for the linear endpoint but is not accepted as saturation parity.

## Verification

- `test_density_gradient_quantum_potential.exe`: 15 assertions, 6 cases pass.
- `test_newton_solver.exe "[density_gradient]"`: 29 assertions, 2 cases pass.
- Frozen saturation endpoint: seven Newton iterations, 0.274% current error.
- Updated linear self-consistent endpoint: three quantum outer iterations,
  final raw change `2.87e-4 V`, and `0.229%` current error.

## Eq. 233 and global potential-like follow-up

The O-2018.06 manual equation was transcribed and implemented as the opt-in
`interface_boundary: sentaurus_step` mode.  The implementation includes the
stable step function `f(x)=(expm1(x)-x)/x^2`, its analytic derivative, the
unsolved-side DOS mass, solved-side gamma, BGN conduction-band share, and the
semiconductor-side normal driving-potential term.  Unit tests cover the
removable singularity, analytic derivative, boundary Jacobian, and JSON
controls.

This mode is not yet accepted for the original SingleDevice deck.  On the
frozen saturation DD state, its first inner solve still reaches the 500-step
limit.  The best tested discretization reduces the raw infinity-norm change
from the homogeneous-Neumann `24.803958 V` to `3.193316 V`, but does not meet
the `5e-4 V` gate.  Reversing the boundary normal worsens the change to
`6.358392 V`; omitting the reconstructed normal drive worsens it to
`6.150184 V`, so both alternatives are rejected.

The manual also establishes that the original deck enables
`eQuantumPotential` globally.  Therefore Eq. 233 is a useful solved/unsolved
boundary feature, but is not the exact interface contract for this case.  A
global potential-like prototype was added: it solves one continuous
`Lambda-(psi+affinity+DOS drive)` unknown while evaluating region-side Lambda.
The implementation now persists this continuous unknown as
`electron_quantum_potential_like_V` in restart CSV files while retaining
backward compatibility with Lambda-only restart files.

The apparent need for two interface Lambda restart values was rejected after
checking the raw per-region TDR export. Across all 216 shared region-pair node
occurrences, `eQuantumPotential`, `ConductionBandEnergy`, and their sum are
bitwise identical on the two exported sides. The Sentaurus node field uses a
single interface trace and resolves the material transition at adjacent
interior nodes. Vela follows the same conforming-node convention and evaluates
the material drive per cell away from the interface trace.

The direct weak form also corrected a dimensional error in the prototype:
with `C=gamma*hbar^2/(6*m*q)` and dimensionless Eq. 231 drive `w`, the reaction
coefficient is `2/C`, not `2*Vt/C`.

## Current endpoint result after continuous-state restart

- The one-material global Eq. 231 oracle converges, returns the continuous
  state, and reproduces the same solution after a restart from that state.
- The saturation endpoint (`Vg=2.2 V`, `Vd=1.1 V`) still reaches the 500-step
  quantum inner limit: residual `1.0538e6`, last accepted update `0.05 V`, raw
  Lambda change `3.5866 V`.
- The linear endpoint (`Vg=2.2 V`, `Vd=0.1 V`) also reaches the 500-step inner
  limit: residual `1.1854e6`, last accepted update `0.05 V`, raw Lambda change
  `2.0993 V`.

Therefore the requested restart representation is implemented and verified,
but the two endpoint oracles are not accepted. The remaining failure is in the
multi-material nonlinear/interface discretization, not missing region-side
restart data. The all-material mode remains experimental and is not promoted
to the default SingleDevice validation path.

## Re-derivation and manufactured/end-point rerun

The insulator specialization was re-derived directly from Eq. 231. Writing
`S=Phi/q`, `phin=-EFn/q`, and `Vt=1/beta/q`, the dimensionless drive is

`w = (-xi*phin - S + (eta-1)*psi)/Vt`.

Consequently the semiconductor uses `w=(-phin-S)/Vt`, while an insulator with
`xi=eta=0` uses `w=(-psi-S)/Vt`. Since
`S=-psi-chi-PhiDOS+Lambda`, the explicit electrostatic term cancels in the
insulator; Lambda remains a region-side reaction quantity and `S` is the
shared interface trace. The cell assembler now evaluates that reconstruction
with the cell-side material, rather than an ownership-selected nodal material.

Two-material manufactured tests verify the cancellation under a nonzero
linear electrostatic potential. A second manufactured solution verifies the
`theta=0.5` identity via `u=exp(w/2)`, and a regression prevents a vanishing
line-search scale from being reported as convergence. The exponential form is
retained as an opt-in manufactured-solution operator; `p1_direct` remains the
endpoint default because the imported oxide field contains element-scale
quantum-potential changes above 1 V.

The restart converter now removes the Sentaurus global band-energy origin
(`4.62577622745745 V` in both endpoint exports), includes the 300 K Silicon
DOS mass `1.0618016171622988`, and preserves the all-material quantum field.
After those corrections, the direct fixed-state maximum residual is
`7.84515e3` at the linear endpoint and `7.84514e3` at saturation; both endpoint
oracles still fail. The common hotspot is Silicon node 1739. Sentaurus exports
approximately zero Lambda there, but the Vela constant-affinity material drive
reconstructs `-22.7403 mV`, so the remaining residual is reaction dominated.
This narrows the unclosed feature to Sentaurus-equivalent doping/BGN-dependent
conduction-band drive, rather than the `xi=eta=0` transformation or interface
restart trace.

Current verification: 46 assertions in 13 density-gradient cases pass; 34
assertions in two Newton density-gradient cases pass; the restart converter
regression passes. Both `Vg=2.2 V` endpoint configurations execute and retain
the expected `electron_density_gradient_max_iterations` failure.

## Material-specific and BGN drive follow-up

The TDR fields independently establish that both Silicon and PolySilicon use
the same un-narrowed affinity in this deck:
`ElectronAffinity - 0.5*BandgapNarrowing = 4.072740384615385 eV` at every node
in both materials. Vela now stores that base value in the SingleDevice
material fixture and adds `0.5*DeltaEg` per cell vertex during Eq. 231
assembly. Material JSON also supports independent
`electron_quantum_gamma` and `electron_quantum_dos_mass_ratio`; Silicon and
PolySilicon are currently both `3.6` and `1.0618016171622988`, while SiO2 is
`1.0` and `0.42`.

The change removes the previous high-doping Silicon hotspot. The fixed-state
direct-P1 maximum residual is now `2.08020e3` for the linear endpoint and
`2.08019e3` for saturation. Both endpoint solves still fail the quantum inner
gate. The common maximum has moved to the SiO2 reoxidation region, where the
exported quantum correction is about `-2.39 V` and the direct expanded
`theta*|grad(w)|^2` term remains large. Therefore the requested BGN and
Si/Poly parameter development is implemented and verified; the remaining
endpoint mismatch is no longer led by the semiconductor material drive.
