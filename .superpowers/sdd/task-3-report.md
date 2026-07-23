# Task 3 Report: Generalize the supplemental Sentaurus field matrix safely

## Status

Completed and committed as `0501c99` (`0501c99ef46cb07156e467f11638300147c44a85`).
Only the Task 3 exporter and its existing state-export regression module were
committed. The unrelated untracked `docs/validation/figures/` directory was
preserved.

## RED

Before production changes, the focused state-export module was run with the
new forty-state and version-provenance tests. It failed as intended:

```text
TypeError: validate_state_matrix() takes 1 positional argument but 2 were given
KeyError: 'sentaurus_version'
Ran 33 tests
FAILED (errors=2)
```

This demonstrated the missing declared-matrix API and absent version binding.

## GREEN

The exporter now:

- accepts a nonempty, finite, duplicate-free caller-declared matrix over valid
  minimal6 topologies while preserving the omitted-argument legacy six-state
  validator;
- records a sorted JSON `expected_matrix` at preparation, binds it to the
  prepared states, rejects mutation during execution, and threads it through
  recovered-archive validation;
- retains exact bias checking and rejects missing/changed states; and
- captures a nonempty `sdevice -version` result, records it in both manifest
  and states, and fails closed on missing or mixed versions.

Existing V2 member-hash, field-contract, recovery, and audit-provenance tests
remain in the focused state-export module and passed unchanged.

## Verification

```text
python tests\regression\test_pn2d_minimal6_state_export.py
Ran 33 tests in 18.431s
OK

python -m unittest discover -s tests\regression -p test_pn2d_minimal6_inverse_inputs.py
Ran 12 tests in 13.272s
OK

git diff --check
(no output)
```

Running the Task 2 input module directly is not its supported invocation: its
importlib fixture needs unittest discovery to establish the repository package
context. The prescribed discovery command above passed.

## Self-review

- The forty-state `(sketch, mirror) x (-1 .. -20 V)` fixture validates exact
  declared coverage and preserves the 1e-12 V final-bias tolerance.
- The default validator still requires exactly the legacy six states and keeps
  its existing `exact six-state matrix mismatch` error contract.
- `prepare_exports()` serializes canonical sorted matrix pairs; execution
  snapshots that record and restores it before writing a failed partial
  manifest if an executor mutates it.
- No remote Sentaurus run, V2 schema change, source/include change, generated
  output change, or Task 4 work was performed.
## Review-fix wave (2026-07-18)

### RED

After review, the expanded focused state-export suite reproduced three binding
failures before production edits:

```text
ValueError: recovered archive expected_matrix does not match caller declaration
AssertionError: ValueError not raised
AssertionError: ' O-2018.06-SP2 ' != 'O-2018.06-SP2'
Ran 38 tests
FAILED (2 failures, 1 error)
```

Those failures proved that caller order was treated as semantic, synchronized
in-memory mutation could replace the prepared contract, and padded versions
were retained verbatim.

### GREEN and commit

Commit `8ac29b5edd47494de306f064e4d285cda6a71c33` binds `run_exports()` to the
original on-disk prepared manifest, canonicalizes caller/recovered matrices
only after duplicate validation, and normalizes nonempty Sentaurus versions.
New-format recovery now requires matching manifest and per-state provenance;
legacy manifests without `expected_matrix` retain optional version behavior.

Focused verification:

```text
python tests\regression\test_pn2d_minimal6_state_export.py
Ran 38 tests in 19.999s
OK

git diff --check
(no output)
```

The added discriminating tests cover invalid, empty, nonfinite, duplicate, and
unknown-topology declarations; missing and extra declared states; generalized
recovery; caller ordering; missing, mixed, and whitespace versions; and a
synchronized in-memory matrix/state mutation before `run_exports()`.