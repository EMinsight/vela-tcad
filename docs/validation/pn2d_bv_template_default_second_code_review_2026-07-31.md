# PN2D BV template-default second independent code review

Date: 2026-07-31

Final verdict:

- code and acceptance infrastructure: `APPROVE`;
- SG/Laux production-default change: `REJECT`.

## Scope

This review examined the actual template/profile patch, generator validation,
exact-lattice evidence writer, contract-domain analyzer, M0/M2 evaluator,
negative tests, and machine acceptance.

## Findings and disposition

1. The acceptance is failed and authorizes no default surface.  The template
   default must therefore be rolled back.  **Resolved.**
2. The atomic profile design is sound, rejects mixed/omitted combinations,
   preserves a single-command legacy rollback, and leaves global C++ defaults
   unchanged.  **Retained.**
3. Parity inputs initially needed stronger binding to Vela branch outputs and
   explicit Sentaurus aggregates.  **Resolved with path/hash checks and a
   self-consistent wrong-artifact negative test.**
4. Branch configs/IVs and render manifest/base config needed a complete hash
   chain.  **Resolved by recomputing branch hashes and regenerating the base
   config from the bound manifest.**
5. State determinism initially accepted partial matching dictionaries.
   **Resolved by requiring passed manifests, exact branch and full execution
   lattice coverage, and actual snapshot-file hash checks.**
6. Process probes initially required only existence and A/B hash identity.
   **Resolved by also requiring machine columns and every contract bias.**
7. CLI opt-in current-support metadata could report the legacy base fields.
   **Resolved by recording the effective CLI SG/Laux fields.**
8. Strict JSON could emit non-standard non-finite values.  **Resolved with
   finite metric checks and `allow_nan=false`.**

## Verification

- focused Python tests: 28/28 passed;
- Release CTest: 506/506 passed;
- final default render: legacy;
- SG/Laux explicit opt-in render: complete atomic profile;
- generated build and simulation output remains untracked.

The infrastructure is acceptable to retain, but the SG/Laux production
default is not approved.

The final read-only recheck found no unresolved P1/P2 item after the
evidence-chain fixes and rollback.
