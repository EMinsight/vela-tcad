# PN2D minimal6 Task 5 Design

## Goal

Make the fixed-state audit report fail closed on provenance and source-integrity
tampering, while preserving independent Python formulas and producing the
required synthetic report package.

## Design

`audit_pn2d_minimal6_fixed_state.py` remains the single report producer. Before
creating any output directory it validates the immutable Task 4 producer,
source commit, six replay identities, arguments, input hashes, committed output
hashes, and fresh replay output hashes. Vela aggregate columns are sums of raw
Task 4 CSV values; Python aggregates are computed independently and compared by
the existing hybrid gate. Geometric tiny-zero normalization is allowed only for
geometrically zero partial volumes.

The focused regression suite supplies adversarial mutations for each contract,
and the synthetic CLI is the end-to-end acceptance check. A Task 5 report records
RED, GREEN, replay, CLI, and independent review evidence. No real Sentaurus
state or Task 6 output is changed.

## Acceptance criteria

- Focused Task 5 tests pass with zero failures/errors.
- Synthetic report contains 36 node, 54 edge, and 24 triangle rows and 14 QA'd figures.
- Provenance replay is actually executed and reports PASS before artifacts are written.
- Maximum C++/Python formula hybrid error is below `5e-12`; state parity is below `1e-12`.
- Git diff is clean of whitespace errors.
