# Task 4 Report: Fixed-State Formula-Difference CLI And Counterfactual Engine

> **Authoritative supersession (2026-07-16 terminal review):** The final section of this report supersedes every earlier causal-availability, factor-contribution, residual, dominance, interaction-count, and artifact-hash claim. Earlier RED/GREEN and provenance history remains historical evidence only; do not use its numerical conclusions.

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

## Review-fix verification (2026-07-16)

### Scoped review and TDD correction

The review-fix wave started from `f128fe1`. The inherited working diff covered the formula CLI, counterfactual engine, interaction plot, and focused regression. The plot change is intentional schema alignment from `factor_a/factor_b` to `first_factor/second_factor`; no unrelated production C++ or solver behavior changed.

The eight-stage formula engine has concrete density, gradient, mobility, current, field, alpha, partial-volume, and node-mapping operators. A real-fixture audit nevertheless found one false availability: at sketch `-19 V`, the raw impact-field replacement differs by `4462.635933946973 V/m`, but direct exported scalar alpha bypasses the field dependency and produced exactly `0.0 dex`. The audit also records that Vela DeMan coefficient provenance is unavailable. The CLI now keeps the field value in the ledger but marks `impact_driving_field` unavailable with `direct exported alpha lacks coefficient provenance for independent impact-field replay`; dominance therefore fails closed.

TDD evidence:

- RED: the real two-run CLI acceptance failed because `impact_driving_field` was `available` instead of `unavailable` (`1` test, `17.986s`).
- The first minimal placement exposed the preliminary-engine complete-map invariant before artifact creation. Root-cause tracing moved the availability guard after source-mapping normalization rather than weakening the engine.
- GREEN: targeted real CLI acceptance passed (`1/1`, `19.303s`).
- Fresh full formula suite passed (`23/23`, `59.004s`).
- Fresh Task 2 contracts plus Task 3 physics passed (`39/39`, `0.202s`).
- Fresh Task 1 state-export plus fixed-state audit/import boundary passed (`67/67`, `36.654s`).
- `compileall -q -f` exited `0`; direct imports of the CLI, counterfactual, plot, and regression modules printed `imports OK`.
- Configured Windows Git (`core.autocrlf=true`) reported `git diff --check` exit `0`.

### Provenance and adversarial evidence

Before creating the output directory, the CLI requires matching state/audit schemas, exact six-state association, equal Task 4 provenance, PASS replay/report command status, passing summary gates, the complete state-tree hash map, successful production replay verification, and byte-exact regeneration of `node_state.csv`, `edge_audit.csv`, and `triangle_audit.csv`.

Focused tests verify failure before artifacts for six CLI file mutations: missing velocity field, wrong `ImpactIonization` unit, reversed triangle topology, duplicate node ID, state-manifest mutation, and audit-CSV mutation. Separate gates reject missing/fabricated audit binding, undeclared dependencies, and a false native `SourceKind`. Reconstructed Sentaurus alpha-current and Vela alpha-flux sources remain `derived`; only raw Sentaurus `ImpactIonization` is native.

### Eight-factor and closure evidence

The synthetic operator test proves every parameterized stage can independently change the final source. The real fixture reports actual controls or a typed unavailable reason:

| Factor | Maximum absolute path contribution (dex) | Real-fixture status |
|---|---:|---|
| `ni_eff/BGN` | `0.013033088418356566` | available |
| `gradient_recovery` | `0.20912350759364517` | available |
| `mobility` | `0.0` | available; exact equal-input control |
| `current_semantics` | `0.17130567269615632` | available |
| `impact_driving_field` | `0.0` | unavailable; coefficient provenance is absent |
| `alpha_law` | `299.35969964281963` | available direct alpha-output replacement; not rankable while field is unavailable |
| `partial_volume` | `9.64327466553287e-17` | available C++/independent-Python control |
| `source_to_node_mapping` | `0.01176195976434569` | available |

- Exact identities: `36 node / 54 edge / 24 triangle`.
- Maximum eligible forward/reverse closure error: `5.684341886080802e-14 dex` (`<1e-10`).
- Maximum path-residual versus named `sentaurus_internal_semantics_residual` mismatch: `6.927791673660977e-14 dex`.
- Dominance: `insufficient_data`, reason `unavailable factors: impact_driving_field`; no `dominant_factor` is emitted.

Two fresh runs produced identical artifact hashes:

| Artifact | SHA-256 |
|---|---|
| `quantity_ledger.csv` | `b5551b5eadf6edf7bf7b41e8404745295ef2290ae80c08b655568f22983ba2bb` |
| `factor_waterfall.csv` | `51860e5e20327eb30ca8c8163674d8a0e419aecbd7957837417258c95372cc84` |
| `root_cause_summary.json` | `cef310b883ae89796fe0a36d7cc27154aff05534061626039aafa8bf174e138f` |
| `root_cause_summary.md` | `1723ccd19357244e1590b62c9cfad077f7fcde34df68a1654c56e3affb9e0402` |

### Interaction evidence and remaining concern

No adjacent pair in the real fixture crosses the strict `>0.3 dex` forward/reverse trigger, so no real interaction is fabricated. The focused interaction case uses `baseline=1`, `A=10`, `B=2`, and `A+B=50`, giving `log10(2.5) = 0.3979400086720376 dex`; it renders both `interaction.png` and `interaction.pdf` through the `first_factor/second_factor` schema.

The remaining limitation is explicit rather than hidden: direct exported alpha cannot independently identify the impact-driving-field effect without bound Vela coefficient provenance. The large direct alpha-output contribution must not be interpreted as a dominant causal factor; the fail-closed `insufficient_data` result is the reviewed outcome.

## Terminal causal-attribution review (authoritative, 2026-07-16)

This section supersedes the earlier Task 4 and review-fix causal results, including the `03391e1` hashes and the stale `299.35969964281963 dex` alpha attribution. Provenance and adversarial gates remain valid, but all earlier availability, contribution, residual, dominance, interaction-count, and hash values are historical.

### Corrected design

- Removed the target/preliminary calibration entirely. `source_to_node_mapping` weights are never scaled to force the Sentaurus alpha-current reconstruction.
- Preserved the actual conservative Vela and independent-Python mapping controls in the ledger. Their maximum partition/local conservation error is `0.0`.
- Marked `source_to_node_mapping` unavailable because Sentaurus mapping weights are not independently exported. Its unexplained total difference remains in the named residual.
- Marked `ni_eff/BGN` unavailable because raw carrier-density averaging is not an independently inferred, formula-bound `ni_eff/BGN` replay. Density interpolation remains visible only as ledger evidence.
- Marked both `impact_driving_field` and `alpha_law` unavailable because direct exported alpha confounds the field and law while bound Vela coefficient provenance is absent.
- Corrected `direction_rad` ledger units from `V/m` to `rad`.
- Kept the generic eight-operator engine and its fully supplied synthetic activation test unchanged.

### TDD and verification

- RED: the target-independence test changed mapping weights from `0.009763849402817335` to `0.16598543984789466` when only the target was multiplied by `17`; the real CLI also emitted `V/m` instead of `rad` (`2` tests, `29.052s`, two expected failures).
- GREEN: targeted conservative-mapping and real-CLI acceptance passed (`2/2`, `25.172s`).
- Fresh Task 4 formula suite passed (`24/24`, `63.305s`).
- Fresh Task 2 contracts plus Task 3 physics passed (`39/39`, `0.199s`).
- Fresh Task 1 state-export plus fixed-state audit/import boundary passed (`67/67`, `35.183s`).
- `compileall -q -f` exited `0`; direct imports of the CLI, counterfactual, plot, and regression modules printed `imports OK`.
- Configured Windows Git (`core.autocrlf=true`) reported `git diff --check` exit `0`.

### Authoritative factor availability

| Factor | Maximum absolute path contribution (dex) | Real-fixture status |
|---|---:|---|
| `ni_eff/BGN` | `0.0` | unavailable; raw density averaging is not inferred ni/BGN |
| `gradient_recovery` | `0.21386353869419972` | available |
| `mobility` | `0.0` | available; exact equal-input control |
| `current_semantics` | `0.22165744282050193` | available |
| `impact_driving_field` | `0.0` | unavailable; coefficient provenance absent |
| `alpha_law` | `0.0` | unavailable; exported alpha confounds field and law |
| `partial_volume` | `4.821637332766436e-17` | available C++/independent-Python control |
| `source_to_node_mapping` | `0.0` | unavailable; only conservative control weights exist |

All four unavailable factors contribute exactly zero on both paths. Dominance is `insufficient_data` with reason `unavailable factors: alpha_law, impact_driving_field, ni_eff/BGN, source_to_node_mapping`; no `dominant_factor` is emitted.

### Exact closure and named residual

- Exact identities: `36 node / 54 edge / 24 triangle`.
- Maximum eligible forward/reverse closure error: `0.0 dex`.
- Maximum path-residual versus named `sentaurus_internal_semantics_residual` mismatch: `0.0 dex`.
- Maximum mapping conservation error: `0.0`.
- Direction units: exactly `rad`.

| State | Named residual (dex) |
|---|---:|
| sketch `0 V` | `286.7681499114667` |
| sketch `-12 V` | `290.6705668672648` |
| sketch `-19 V` | `293.62234705748955` |
| mirror `0 V` | `286.89502314785227` |
| mirror `-12 V` | `290.7974401036503` |
| mirror `-19 V` | `293.74922029387506` |

The large residual is the reviewed result: missing independent causal inputs dominate the unexplained gap. It is not calibrated away or attributed to alpha, ni/BGN, field, or mapping.

### Determinism, provenance, and interactions

Two fresh runs produced identical authoritative hashes:

| Artifact | SHA-256 |
|---|---|
| `quantity_ledger.csv` | `50afbc13f947990469994fee26dc3d7bef6bf7d581bec737245dc13bbd85f132` |
| `factor_waterfall.csv` | `7874740f80eff42de9e7ffaedafa9a07477677f2d5fc7ac436c5a0a976ec33ad` |
| `root_cause_summary.json` | `45000b674741538d12983b4d43af3d8bae18ca0f5c69eb04067b4377d2966208` |
| `root_cause_summary.md` | `1723ccd19357244e1590b62c9cfad077f7fcde34df68a1654c56e3affb9e0402` |

Strict state/audit binding, production replay, byte-exact audit regeneration, six CLI file mutations, undeclared-dependency rejection, and false-native-kind rejection all remain covered. No real adjacent pair crosses the strict `>0.3 dex` trigger, so no real interaction is fabricated. The synthetic interaction contract remains `log10(2.5) = 0.3979400086720376 dex` and renders both PNG and PDF through `first_factor/second_factor`.

### Terminal concern

The package is integrity-complete but scientifically inconclusive. A causal ranking requires independently bound ni/BGN formulas, Vela alpha coefficients that permit field/law separation, and an exported Sentaurus source-mapping operator. Until those exist, the authoritative outcome is exact closure with a large named residual and `insufficient_data`.
