# Task 9 canonical inverse-input sealer report

## Outcome

Implemented a diagnostic-only, fail-closed adapter from the actual Task 8/9
schemas to three sealed `vela.pn2d_minimal6_inverse_input.v1` roots. The
adapter re-imports hash-bound/raw-sealed Sentaurus TDRs with the explicitly
supplied importer, validates canonical Minimal6 geometry and field contracts,
and atomically publishes only after `load_input_bundle()` accepts all three
staged roots.

Vela provides native nodal psi, electron/hole quasi-Fermi potential, and
carrier density only. Electric field, electron/hole current density,
electron/hole avalanche alpha, and impact-ionization generation remain blank
and therefore typed `MISSING_FIELD`. No edge/cell-to-node transform is
invented.

## TDD evidence

RED:

- Direct report/verifier execution from a non-repository cwd failed 1/1 with
  `ModuleNotFoundError: No module named 'scripts'`.
- The sealer contract failed because
  `scripts/pn2d_minimal6_diagnostics/inverse_input_sealer.py` did not exist;
  six dependent tests were skipped.
- Adding the sealer CLI to the direct-script test failed because
  `scripts/seal_pn2d_minimal6_inverse_inputs.py` did not exist.

GREEN:

- Sealer/CLI contract: 8/8 passed in 149.393 s.
- Existing inverse-input contract: 12/12 passed in 24.766 s.
- Report plus independent verifier: 5/5 passed in 42.193 s.
- State exporter/resume contract: 45/45 passed in 34.631 s.
- Scoped Python compilation passed.
- `git diff --exit-code a5524cf -- include src` passed with no output.

The sealer tests cover exact matrices, unit conversion, canonical coordinate
remapping, typed Vela gaps, byte-identical dual output, irrelevance of old
unhashed Task 8 field CSVs, source/member hash rejection, path escape,
incomplete/version/interpolation rejection, importer failure, imported
unit/component/mapping/geometry rejection, existing-output preservation, and
sealed-state tampering.

## Provenance and limitations

- Task 8 did not bind the Vela runner executable at execution time. The sealed
  Vela root labels its supplied runner hash as `post_hoc_observed`.
- Task 8 did not declare a remote Sentaurus release/binary. The base Sentaurus
  root labels this as `not_declared_by_source_manifest`; it does not infer the
  release from the supplemental run.
- The supplemental root requires exact `O-2018.06-SP2` provenance and all 40
  states passed before sealing.
- Current relevant production-source bytes are copied and hashed; the
  phase-base production no-change proof remains a separate explicit git gate.
- Missing Vela native nodal E/J/alpha/source evidence can only support
  `insufficient_data` or confounded conclusions, not an identified formula.

## Planned invocation

```powershell
python scripts\seal_pn2d_minimal6_inverse_inputs.py `
  --vela-sweep-root=<task8-vela-root> `
  --sentaurus-sweep-root=<task8-sentaurus-root> `
  --supplemental-root=<completed-task9-export-root> `
  --output-root=<sealed-input-root> `
  --importer=build-release\sentaurus_import.exe `
  --vela-executable=build-release\vela_example_runner.exe `
  --phase-base=a5524cf
```

Use `<sealed-input-root>\vela`, `sentaurus`, and `supplemental` as the three
report CLI roots. The report, verifier, and sealer CLIs now bootstrap the
repository path before importing `scripts.*`, so direct invocation works from
any cwd.
