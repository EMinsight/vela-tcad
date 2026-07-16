# Task 4 Report: Fixed-State Formula-Difference CLI And Counterfactual Engine

## Status
Implemented from clean base `658b6b2de297a6a61cbb5e936e227842748e23f7` without production C++ changes.

## Coverage audit
The historical partial had source integration, ledger rows, factor names, and closure/dominance helpers. Missing pieces were the strict residual name, executable substitutions, declared downstream invalidation, reverse-path interactions, ordered waterfall rows, exact audit-root validation, typed reconstructed-source labels, and two-run hash determinism.

## RED
- Baseline: 15 tests, one failure: `ValueError: residual record has the wrong name`. Fixed the producer; did not weaken Task 2 validation.
- Reverse interaction: expected `reverse_adjacent` was absent. Both path orders are now evaluated.
- DAG engine: import failed because `DependencyCounterfactualEngine` did not exist. The engine now restricts reads to declared parents and recomputes only changed plus downstream operators.
- CLI acceptance: -19 V `ni_eff` raised `OverflowError`. The Task 3 formula remains primary; a log-domain fallback emits finite log10 values and marks only unrepresentable linear values unavailable.

## GREEN
```text
python -B -m unittest tests.regression.test_pn2d_minimal6_formula_difference -v
Ran 17 tests in 8.119s
OK

python -B -m unittest tests.regression.test_pn2d_minimal6_diagnostic_contracts tests.regression.test_pn2d_minimal6_diagnostic_physics -v
Ran 39 tests in 0.220s
OK

git diff --check
exit 0, no output
```
Runtime and standalone Draft 2020-12 formula schema validation both pass.

## Acceptance evidence
The committed six-state fixture is augmented at runtime with deterministic `ImpactIonization`, velocities, electron/hole/mean ionization integrals, and temperature.

- Exact identities: `36/54/24`.
- Native SourceKind: `sentaurus`; both reconstructions: `derived`.
- Forward order is the declared eight-factor order; reverse is exact reverse.
- Maximum -12/-19 V closure error: `0.0 dex`.
- Maximum path-residual vs named Sentaurus semantics-residual difference: `2.4868995751603507e-14 dex`.
- Symmetric dominance spans both biases/topologies; the 25% residual gate removes the dominant factor.
- No synthetic interaction is emitted because no adjacent pair has a forward/reverse difference strictly above `0.3 dex`.

Two runs produced identical hashes:

| Artifact | SHA-256 |
|---|---|
| `quantity_ledger.csv` | `11dc0eec76dfac0e1f9240e1315028ba30c9b69951d925ab6ad7310f7628d388` |
| `factor_waterfall.csv` | `26fc1133c333276fb0b91e608d8ccec7b215446fb78c1eee54ca75d19c75e111` |
| `root_cause_summary.json` | `e34acac2500313d51e8527efa2a8a6cdb693194b2df83b4aaf85d01734cc70d5` |
| `root_cause_summary.md` | `1723ccd19357244e1590b62c9cfad077f7fcde34df68a1654c56e3affb9e0402` |

## Adversarial coverage
Missing fields, wrong units, inexact bias, clockwise topology, duplicate nodes, hash mutation, undeclared dependencies, and reconstructed-as-native labels all fail closed.

## Files
- `scripts/pn2d_minimal6_diagnostics/counterfactual.py`
- `scripts/diagnose_pn2d_minimal6_formula_difference.py`
- `tests/regression/test_pn2d_minimal6_formula_difference.py`
- `.superpowers/sdd/task-4-report.md`

## Concerns
The deterministic fixture closes every non-current control, so only current semantics is nonzero and no interaction is manufactured. Real inputs may activate other factors under the same DAG/gates. Linear `ni_eff` beyond IEEE-754 range is unavailable while finite log10 values remain recorded.
