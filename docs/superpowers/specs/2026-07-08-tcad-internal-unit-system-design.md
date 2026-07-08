# TCAD Internal Unit System Design

## Purpose

The current `unit_scaling` mode interprets public deck values in common TCAD
units, immediately converts them to SI, and then builds dimensionless solver
unknowns from those SI values. The requested change is to remove SI as the
internal physical intermediate for `unit_scaling`.

After this change, `unit_scaling` uses TCAD display units as the internal
physical unit system:

| Quantity | Legacy mode internal unit | `unit_scaling` internal unit |
| --- | --- | --- |
| length | m | um |
| area | m^2 | um^2 |
| volume | m^3 | um^3 |
| concentration | m^-3 | cm^-3 |
| sheet density | m^-2 | cm^-2 |
| mobility | m^2/(V s) | cm^2/(V s) |
| diffusivity | m^2/s | cm^2/s |
| electric field | V/m | V/cm |
| inverse length | m^-1 | cm^-1 |
| current density | A/m^2 | A/cm^2 |
| potential | V | V |
| temperature | K | K |
| energy | eV | eV |

The design keeps legacy SI behavior unchanged when `scaling` is omitted.

## Non-Goals

- Do not introduce a second public TCAD mode. The target is to change
  `scaling.mode = "unit_scaling"` semantics after tests and documentation are
  updated.
- Do not tune physical models or PN2D calibration while changing units.
- Do not preserve misleading SI-named output columns in paths that no longer
  compute SI internally without an explicit conversion at the output boundary.
- Do not rewrite unrelated solver algorithms.

## Architecture

Introduce a small unit-system layer, tentatively `PhysicalUnitSystem`. Store it
by value on `UnitScalingConfig` so `parseUnitScalingConfig` returns both the
mode and the active internal unit system as one small value object.

`UnitScalingConfig` should stop exposing `lengthToSI`, `concentrationToSI`,
`mobilityToSI`, and similar APIs as the main parser surface. Replace them with
internal-unit names:

- `lengthToInternal`
- `concentrationToInternal`
- `sheetDensityToInternal`
- `mobilityToInternal`
- `electricFieldToInternal`
- `inverseLengthToInternal`
- `surfaceFieldCoefficientToInternal`

For legacy mode these methods are identity because legacy internal units are
SI. For `unit_scaling` these methods are also identity for deck values that are
already written in TCAD units. This identity behavior is intentional: the calls
are a semantic boundary and a future extension point for any third unit mode.
Implementation code should still call them rather than bypassing the parser
surface, because direct reads tend to recreate hidden SI assumptions.

Add output-boundary helpers for display and compatibility:

- `internalLengthToMeters`
- `internalConcentrationToM3`
- `internalElectricFieldToVPerM`
- `internalCurrentDensityToAPerM2`
- `internalCurrentPerDeviceDepthToAPerUm`

These helpers are only for outputs, cross-mode comparisons, and compatibility
tests. They are not used as the default path into assemblers.

## Unit Constants

Each unit system needs factors from its internal units to SI only for deriving
dimensionally correct constants. Those factors live inside the unit-system
layer, not in parser code.

For the TCAD internal system:

- `length_m_per_internal = 1e-6`
- `area_m2_per_internal = 1e-12`
- `volume_m3_per_internal = 1e-18`
- `concentration_m3_per_internal = 1e6`
- `mobility_m2_per_V_s_per_internal = 1e-4`
- `field_V_per_m_per_internal = 1e2`
- `inverse_length_m_inv_per_internal = 1e2`
- `current_density_A_m2_per_internal = 1e4`

These factors should be centralized so equations do not scatter literals such
as `1e-6`, `1e4`, or `100`.

Preferred implementation path: keep the solver state and parsed physical
fields in the active internal unit system, but allow `UnitScalingSystem` to use
the centralized factors above to derive `lambda2`, `J0`, and `R0`. This keeps
SI conversion isolated inside the unit-system/reference-scale layer instead of
reintroducing SI as the physical state passed through parsers, assemblers, or
models.

Composite constants are required where TCAD base units mix `um` geometry with
`cm` density or mobility units:

| Term family | Internal operands | TCAD composite factor |
| --- | --- | --- |
| Volumetric charge, generation, recombination geometry | `cm^-3 * um^3` equivalent, including 2-D area times the unit-depth convention | `1e-12` |
| Interface, sheet, and Neumann charge geometry | `cm^-2 * um^2` equivalent, including 2-D edge length times the unit-depth convention | `1e-8` |
| Field from mesh coordinate differences | `V / um` to `V/cm` | `1e4` |
| Current or charge per internal device depth | per `um` display | identity for `unit_scaling`; explicit conversion only for legacy SI |

Do not let individual assemblers or diagnostics independently rederive these
numbers. Add named helpers on `PhysicalUnitSystem` for these composite factors
and use them in every charge, source, field, and terminal-current path.

## Scaling System

`UnitScalingSystem` should receive internal physical quantities and a
`PhysicalUnitSystem`.

References:

- `V0 = kT / q` in volts.
- `C0` is an internal concentration scale.
- `L0` is an internal length scale.
- `mu0` is an internal mobility scale.
- `D0 = mu0 * V0`, adjusted so its length unit matches the continuity geometry.
- `E0 = V0 / L0`, expressed in the internal electric-field unit.
- `J0` is expressed in the internal current-density unit.
- `R0` is expressed in the internal volumetric generation/recombination unit.
- `lambda2` is the dimensionless Poisson coefficient in the current internal
  unit system.

The implementation must not derive `lambda2`, `J0`, or `R0` by assuming SI
mesh length, SI concentration, or SI mobility. Compute them by converting only
the reference constants through the centralized unit-system factors inside
`UnitScalingSystem`. That keeps the change close to the existing scaling
architecture while preventing SI values from leaking into parser or assembler
state.

## Poisson Equation

The Poisson path is the highest-risk part because it couples charge density,
permittivity, geometry, and potential:

`-div(eps grad psi) = q * rho_number`

In TCAD internal units:

- geometry is in `um`,
- doping/fixed charge is in `cm^-3`,
- interface charge is in `cm^-2`,
- potential is in `V`.

The assembler must use unit-system factors to make volumetric and sheet charge
terms consistent with the geometric control volumes. The matrix and RHS
normalization should be documented in code next to the Poisson scaling spec.
Use the named composite factors from `PhysicalUnitSystem` for:

- Poisson main RHS volumetric doping and region fixed charge.
- Poisson fixed/interface charge RHS.
- Gummel/DD Poisson-with-carriers RHS.
- Coupled Newton Poisson residual and Jacobian rows.
- Neumann displacement boundaries.

Acceptance checks:

- A legacy SI deck and an equivalent `unit_scaling` deck produce the same
  physical potential after output conversion.
- A `0.01 um` square cell keeps geometry in `um^2` internally and does not
  construct `1e-16 m^2` node areas as the main geometry representation.

## Drift-Diffusion Equations

The drift-diffusion assemblers should operate on internal physical units and
dimensionless solver unknowns:

- `psi_hat = psi / V0`
- `phin_hat = phin / V0`
- `phip_hat = phip / V0`
- `n_hat = n / C0`
- `p_hat = p / C0`

Scharfetter-Gummel flux and recombination/source terms must be re-derived for
the internal unit system:

- Mobility is `cm^2/(V s)` in `unit_scaling`.
- Electric field and quasi-Fermi gradients are `V/cm`.
- Mesh edge lengths are `um`, so any flux expression using gradients must
  apply the centralized `cm_per_um = 1e-4` relation through the unit-system
  layer.
- Recombination rates and avalanche source integrals must have explicit
  internal units before residual normalization.
- The Poisson-charge composite factors used by the standalone Poisson path
  also apply to the Gummel and coupled Newton Poisson blocks. They should be
  shared through the same helper layer, not duplicated in each assembler.
- Any term that computes field-like quantities from coordinate differences
  must record the unit it has produced before applying output conversions. In
  `unit_scaling`, `abs(delta V) / edge_length` is `V/um`, not `V/m`.

The solver should still return `DDSolution` in internal physical units. For
legacy mode this means SI; for `unit_scaling` this means TCAD units. Consumers
must not assume `DDSolution` is SI without checking the associated
`UnitScalingConfig` or `PhysicalUnitSystem`.

## Materials And Models

Input parsing changes:

- Mesh coordinates in `unit_scaling` are stored as `um`, not converted to `m`.
- Doping and material concentrations are stored as `cm^-3`.
- Mobilities are stored as `cm^2/(V s)`.
- Impact ionization coefficients remain in `cm^-1` and `V/cm`.
- Surface mobility theta remains in `cm/V`.

Model implementations should use the same internal units as their inputs. When
a formula currently assumes SI, update either the formula constants or add a
unit-system coefficient at the model boundary.

Default model parameters are in scope for this change. Several defaults are
currently compiled as SI literals because old internal state was always SI.
Audit and migrate at least:

- Caughey-Thomas and Masetti mobility defaults in
  `include/vela/physics/MobilityModel.h`.
- Selberherr and Van Overstraeten impact-ionization defaults in
  `include/vela/physics/ImpactIonizationModel.h`.
- Hard-coded conversion helpers such as `inverseCmToInverseM` and
  `fieldVPerCmToVPerM` in `src/physics/ImpactIonizationModel.cpp`.

The preferred strategy is to keep canonical model defaults in documented TCAD
units where the source literature or Sentaurus-style parameters are already in
TCAD units, and convert to legacy SI only when constructing legacy-mode model
configs. This avoids treating an omitted JSON field differently from an
explicitly supplied equivalent field.

## Output Policy

Outputs should be explicit about display units.

Recommended policy:

- VTK potential and quasi-Fermi fields: `V`.
- VTK carrier and doping fields in `unit_scaling`: `cm^-3`.
- VTK electric field in `unit_scaling`: `V/cm`.
- VTK current-density vectors in `unit_scaling`: `A/cm^2`.
- DC sweep current display columns: keep `*_A_per_um`.
- If SI compatibility columns are retained, compute them only through explicit
  output conversion helpers and keep their names honest.
- Audit all diagnostics that compute fields or currents directly from mesh
  coordinate differences before applying fixed literals such as `100`, `1e4`,
  or `1e6`. This includes at least `DCSweep.cpp`,
  `ElectricFieldDiagnostics.*`, `ContactCurrent.*`, and the VTK current-density
  export path.

Restart-state CSV policy:

- Keep the existing restart-state header in the first implementation to avoid
  breaking restart readers during the unit-system migration.
- Document that restart-state values are in the active internal physical unit
  system, despite the historical `electrons_m3` and `holes_m3` column names.
- Require restart files to be read with the same scaling mode that wrote them.
  Cross-mode restart conversion is out of scope for this change.

## Compatibility And Migration

This is a breaking semantic change for decks using
`scaling.mode = "unit_scaling"`.

Migration work:

- Update `docs/config_schema.md`.
  - Replace the current statement that Poisson and drift-diffusion assemblers
    continue to receive SI values.
- Update `docs/development_poisson_unit_scaling.md` or replace it with a
  broader internal-unit-system note.
- Update example comments that mention SI normalization.
- Update reference TCAD tools only where they assumed Vela internal fields were
  SI. Candidate files include:
  - `src/tools/alpha_from_sentaurus_f.cpp`
  - `src/tools/pn2d_jacobian_block_audit.cpp`
  - `src/tools/sentaurus_import.cpp`
  - `src/tools/sentaurus_vanoverstraeten_fit.cpp`
  - `src/tools/vanoverstraeten_parameter_sweep.cpp`
- Audit internal diagnostic structures and Doxygen comments with SI suffixes or
  SI claims. Candidate fields include:
  - `ContactCurrentEdgeDiagnostic::edgeLength_m` and `edgeCouple_m`
  - `CoupledDDEdgeFluxDiagnostic::length_m`, `couple_m`,
    `netDopingAvg_m3`, `electricField_V_m`, and
    `electronMobility_m2_V_s`
  - `SweepTransportDiagnostics`
  - `ContinuityBalanceDiagnosticRow`
  - comments such as `DopingModel`'s "All concentrations are in SI units"

Legacy no-`scaling` decks must remain unchanged.

## Test Plan

Follow TDD. Add failing tests before production edits.

Initial red tests:

1. `UnitScalingConfig unit_scaling keeps TCAD values internal`
   - Mesh `x = 0.01` remains `0.01` internally.
   - Doping `1e17` remains `1e17` internally.
   - Mobility `1000` remains `1000` internally.

2. `UnitScalingSystem computes finite TCAD references`
   - With `L0 = 1 um`, `C0 = 1e17 cm^-3`, `mu0 = 1000 cm^2/(V s)`, all
     references are positive and dimensionless coefficients are finite.
   - Exercise `UnitScalingSystem::autoInputsFrom` directly in both legacy and
     `unit_scaling` contexts.

3. `Poisson unit_scaling matches equivalent legacy SI potential`
   - Build one tiny PN or charged slab case in both modes.
   - Compare final physical potential after explicit output conversion.

4. `Box geometry stays in internal units`
   - A `0.01 um` by `0.01 um` cell has area on the order of `1e-4 um^2`, not
     `1e-16 m^2`.

5. `Impact and mobility parsers do not preconvert unit_scaling parameters to SI`
   - `electron_A_m_inv = 1e6` remains the internal `cm^-1` value for TCAD mode.
   - `electron_mu_min_m2_V_s = 100` remains the internal `cm^2/(V s)` value for
     TCAD mode.

6. `Default model parameters migrate with the active unit system`
   - Build `unit_scaling` configs that omit optional mobility and impact
     coefficients.
   - Assert default mobility and impact coefficients equal the intended TCAD
     internal values, not the old SI literals.

7. `Poisson and DD charge factors are explicit`
   - Add small assembler-level tests that isolate volumetric charge and
     interface/sheet/Neumann charge contributions.
   - Assert the TCAD composite factors are applied (`1e-12` for volumetric
     `cm^-3 * um^3` equivalent and `1e-8` for sheet `cm^-2 * um^2`
     equivalent).

8. `Derived output columns are cross-mode consistent`
   - Compare equivalent legacy and `unit_scaling` sweeps for
     `*_A_per_um`, `max_electric_field_V_per_cm`, charge, capacitance, and VTK
     field/current-density diagnostics.
   - Include a case where field values are computed from mesh coordinate
     differences so `V/um` is not accidentally treated as `V/m`.

Regression commands:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
cmake --build build --parallel
ctest --test-dir build --output-on-failure -R "scaling|poisson|mesh|mobility|newton|gummel"
```

Run the full suite before claiming completion:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
ctest --test-dir build --output-on-failure
```

## Implementation Order

1. Add `PhysicalUnitSystem` and internal-unit parser methods.
2. Add tests proving `unit_scaling` no longer converts parser values to SI.
3. Audit and migrate compiled physical defaults and ad hoc conversion helpers
   in mobility and impact-ionization models.
4. Update mesh, doping, material, mobility, and impact parsing to use internal
   methods.
5. Update `UnitScalingSystem` to accept the unit system and compute references
   without SI assumptions leaking into solver state.
6. Consolidate automatic reference-scale input construction. Either make
   `PoissonSimulation.cpp` call `UnitScalingSystem::autoInputsFrom`, or list
   and update every duplicated derivation path explicitly.
7. Update Poisson assembly and tests, including explicit composite charge
   factors.
8. Update DD assembly and tests, including Gummel and coupled Newton Poisson
   blocks.
9. Audit diagnostics and post-processing paths that use coordinate differences
   or fixed conversion literals (`100`, `1e4`, `1e6`).
10. Update output naming/documentation, restart-state documentation, and
    internal diagnostic structure comments or field names.
11. Run focused and full test suites.

## Risks

- Poisson charge scaling can silently drift if permittivity, volume, and
  concentration factors are not derived together.
- SG flux can silently drift because mobility and geometry use different TCAD
  length bases (`cm` for mobility, `um` for mesh).
- Existing regression fixtures that compare against `current_total_A_per_um`
  may change if output conversion and internal current scaling are not kept
  equivalent.
- Restart-state files become ambiguous unless their scaling mode is preserved
  by convention or metadata.
- Compiled model defaults can silently remain in old SI units even after parser
  conversions are removed. Treat defaults as a separate migration target.
- Field and current diagnostics can double-convert when a value computed from
  TCAD-coordinate differences is passed to an SI-named helper such as
  `voltsPerMeterToVoltsPerCm`.
- Duplicate automatic reference-scale derivations can leave Poisson,
  Gummel/Newton, and post-processing on different internal scales if only one
  path is updated.

## Design Decisions

- Restart-state column names stay stable for the first implementation, but
  their documented meaning becomes active internal physical units.
- Existing user-facing JSON field names stay stable even when they contain
  historical SI suffixes. Under `unit_scaling`, their numeric interpretation is
  TCAD internal units.
- Keep common TCAD display columns such as `*_A_per_um`,
  `max_electric_field_V_per_cm`, and VTK `A/cm^2` current density.
- Preserve SI compatibility output only through explicit output conversion
  helpers. Do not feed converted SI values back into solver state or model
  evaluation.
