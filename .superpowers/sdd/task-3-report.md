# Task 3 Report: Independent Physics And Support-Conversion Library

## Status

Completed from clean base `1a88019a8b397d6c75545c26ca8159aba65bf5ee`.
The independent geometry, physics, and support-conversion library now has direct
coverage for every Task 3 edge contract, including the tracked Sentaurus
parameter provenance and an exact comparison against Vela production defaults.

## RED

The first required RED was run before changing implementation code:

```powershell
python -m unittest tests.regression.test_pn2d_minimal6_diagnostic_physics.DiagnosticPhysicsTest.test_edge_to_cell_rejects_native_avalanche_generation
```

Observed result:

```text
TypeError: edge_scalar_to_cells() got an unexpected keyword argument 'quantity'
Ran 1 test in 0.001s
FAILED (errors=1)
```

The expanded focused suite then failed during import because the independent
Vela production-default parser did not exist:

```text
ImportError: cannot import name 'parse_vela_van_overstraeten_defaults'
Ran 1 test in 0.000s
FAILED (errors=1)
```

The other newly direct tests characterize already-correct fail-closed behavior:
degenerate triangles, zero-length edges, an inverse alpha domain with no valid
candidate, inconsistent electron/hole `ni_eff` estimates without an
authoritative average, and partial-volume fractions outside `[0, 1]`.

## GREEN

Implemented the minimum missing behavior:

- `edge_scalar_to_cells(..., quantity=...)` rejects native Sentaurus
  `ImpactIonization` and Vela `AvalancheGeneration` before averaging.
- `parse_vela_van_overstraeten_defaults` reads the numeric production defaults
  from `include/vela/physics/ImpactIonizationModel.h`, converts the SI values to
  the cm-based units used by Sentaurus, and records the header path and SHA-256.
- No second unversioned coefficient table was added.

Focused result:

```text
....................
Ran 20 tests in 0.019s
OK
```

## Parameter provenance and comparison

Tracked source:

`reference_tcad/pn2d_sentaurus2018_minimal6/source/models.par`

SHA-256:

`b4b3ebfdefba530f756f3855d43d7d587720689771d8badc747b61439ed42742`

| Carrier | Coefficient | Sentaurus parsed | Vela production default after SI-to-cm conversion | Result |
|---|---:|---:|---:|---|
| electron | `a_low` | `7.03e5 1/cm` | `7.03e5 1/cm` | match |
| electron | `a_high` | `7.03e5 1/cm` | `7.03e5 1/cm` | match |
| electron | `b_low` | `1.231e6 V/cm` | `1.231e6 V/cm` | match |
| electron | `b_high` | `1.231e6 V/cm` | `1.231e6 V/cm` | match |
| hole | `a_low` | `1.582e6 1/cm` | `1.582e6 1/cm` | match |
| hole | `a_high` | `6.71e5 1/cm` | `6.71e5 1/cm` | match |
| hole | `b_low` | `2.036e6 V/cm` | `2.036e6 V/cm` | match |
| hole | `b_high` | `1.693e6 V/cm` | `1.693e6 V/cm` | match |
| shared | switch field | `4.0e5 V/cm` | `4.0e5 V/cm` | match |

The comparison status is `available`; there are no missing production
coefficients and no mismatches.

## Verification

Task 2 contracts:

```text
python -m unittest tests.regression.test_pn2d_minimal6_diagnostic_contracts
Ran 18 tests in 0.152s
OK
```

Task 2 compatibility/import controls:

```text
python -m compileall -q scripts\pn2d_minimal6_diagnostics
imports OK
```

C++ controls:

```text
build-release\test_impact_ionization.exe
All tests passed (511 assertions in 40 test cases)

build-release\test_cell_reconstructed_avalanche.exe
All tests passed (83 assertions in 14 test cases)
```

Formula gate control:

```text
python -m unittest tests.regression.test_pn2d_minimal6_fixed_state_audit.ReviewContractTests.test_complete_matrix_unique_keys_and_formula_gate -v
Ran 1 test in 0.078s
OK
```

`scripts/audit_pn2d_minimal6_fixed_state.py` remains unchanged with
`FORMULA_LIMIT = 5.0e-12`.

Repository hygiene:

```text
git diff --check
```

No output; exit code 0.

## Changed files

- `scripts/pn2d_minimal6_diagnostics/physics.py`
- `scripts/pn2d_minimal6_diagnostics/support.py`
- `tests/regression/test_pn2d_minimal6_diagnostic_physics.py`
- `.superpowers/sdd/task-3-report.md`

## Concerns

No parameter-coverage concern remains: Vela exposes the full comparable numeric
table. The production-default parser intentionally fails closed if those C++
declarations stop being explicit numeric initializers; such a source-format
change will require updating this independent audit control.

## Interrupted fix wave (2026-07-16)

This wave resumed at `2f178a8fc60eb5a6273279c92ea3af9e8ad0b7ee`
without discarding or reimplementing the existing uncommitted changes.

### RED

Before the interrupted implementation was added, the expanded focused suite
ran 21 tests and reported:

```text
Ran 21 tests
FAILED (failures=2, errors=1)
```

The two failures demonstrated that phonon mismatch and omitted conversion
provenance were not yet enforced. The error demonstrated that the Vela parser
did not yet expose `phonon_energy_eV` for the tracked comparison. These were
expected missing-behavior results, not fixture, import, or syntax failures.

### GREEN

The interrupted implementation added the minimum behavior needed by those
tests:

- parse Vela's numeric `phononEnergy` and retain the header path and SHA-256;
- compare `phonon_energy_eV` for both carriers;
- require recognized, non-native provenance for every node-to-cell and
  edge-to-cell average, rejecting omitted, unknown, `ImpactIonization`, and
  `AvalancheGeneration` provenance; and
- document the `infer_ni_eff` quasi-Fermi sign convention and inverse relations.

The recorded GREEN rerun was:

```text
Ran 21 tests in 0.020s
OK
```

### Phonon units and comparison scope

Sentaurus declares `hbarOmega = 0.063, 0.063` in eV, while Vela declares
`phononEnergy = 0.063` in eV. Unlike inverse-length and electric-field
coefficients, this value needs no SI-to-cm conversion. It is compared directly
for both carriers, with mismatch or unavailable returned on difference or
omission.

Vela's `aScale` and `bScale` are runtime diagnostic multipliers, not material
coefficients in the tracked Sentaurus block. Comparing them to `a(low/high)` or
`b(low/high)` would conflate controls with source parameters, so they remain
explicitly outside this comparison.

Vela's `referenceTemperature_K` and `temperature_K` are runtime inputs. The
tracked Sentaurus block documents the thermal factor in terms of `T0` and `T`
but supplies no numeric values for those inputs. Temperatures therefore remain
outside the direct default comparison. The common phonon energy is included
because both sources explicitly serialize it in eV.

### Fresh verification after resume

Focused physics/support:

```text
python -m unittest tests.regression.test_pn2d_minimal6_diagnostic_physics -v
Ran 21 tests in 0.020s
OK
```

Task 2 contracts and imports:

```text
python -m unittest tests.regression.test_pn2d_minimal6_diagnostic_contracts -v
Ran 18 tests in 0.205s
OK
python -m compileall -q scripts\pn2d_minimal6_diagnostics
imports OK
```

C++ controls:

```text
build-release\test_impact_ionization.exe
All tests passed (511 assertions in 40 test cases)
build-release\test_cell_reconstructed_avalanche.exe
All tests passed (83 assertions in 14 test cases)
```

Formula gate and repository hygiene:

```text
python -m unittest tests.regression.test_pn2d_minimal6_fixed_state_audit.ReviewContractTests.test_complete_matrix_unique_keys_and_formula_gate -v
Ran 1 test in 0.078s
OK
FORMULA_LIMIT = 5.0e-12
git diff --check
no output; exit code 0
```

### Provenance audit and self-review

A repository-wide callsite search found no consumers outside the focused test
module. Every successful `node_scalar_to_cells` or `edge_scalar_to_cells`
callsite passes an allowed provenance (`Potential` or `IonizationCoefficient`);
omitted and unknown values occur only in deliberate rejection tests. The
helpers fail closed when provenance is absent.

Production header provenance:
`include/vela/physics/ImpactIonizationModel.h`
SHA-256 `b1582b41074043207fc5dafcd4c73f95d3fc8b3af3eba4503cddb3113130c008`

Tracked Sentaurus provenance:
`reference_tcad/pn2d_sentaurus2018_minimal6/source/models.par`
SHA-256 `b4b3ebfdefba530f756f3855d43d7d587720689771d8badc747b61439ed42742`

Self-review found no correctness gap in the scoped diff. The allowlist is
intentionally narrow and requires an explicit code-and-test update when a new
intensive quantity is approved for averaging. The regex parser remains a
fail-closed audit over explicit numeric C++ initializers, not a general C++
parser.
