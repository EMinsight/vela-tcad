# PN2D minimal6 Task 5 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Complete the PN2D minimal6 fixed-state audit gate and produce independently verifiable Task 5 evidence.

**Architecture:** Keep the Python audit as the report orchestrator, with strict pre-artifact C++ replay verification and separate Vela/Python source aggregates. Use the existing synthetic fixture and adversarial unittest suite as the executable contract.

**Tech Stack:** Python 3, unittest, NumPy, Pillow, Matplotlib, PyMuPDF, C++20/CMake/Ninja producer.

## Global Constraints

- Do not weaken provenance or skip replay when the producer is unavailable.
- Do not use Python-derived values to populate `vela_*` production columns.
- Do not start Task 6 or alter real Sentaurus state exports.
- Use the MSYS2 UCRT64 toolchain from `D:\msys64`.

### Task 1: Reproduce and isolate current RED

**Files:**
- Test: `tests/regression/test_pn2d_minimal6_fixed_state_audit.py`
- Code: `scripts/audit_pn2d_minimal6_fixed_state.py`

- [x] Run the focused unittest and record each failure category.
- [x] Build or configure the Task 4 producer without changing the contract.
- [x] Re-run only the affected tests to distinguish missing executable from logic defects.

### Task 2: Fix source aggregate contract

**Files:**
- Modify: `scripts/audit_pn2d_minimal6_fixed_state.py`
- Test: `tests/regression/test_pn2d_minimal6_fixed_state_audit.py`

- [x] Add a regression assertion using tolerance for independently summed raw C++ values.
- [x] Ensure `vela_*` sums raw CSV local-edge columns and Python sums Python columns.
- [x] Run the aggregate and formula-gate tests; keep the adversarial mutation failing closed.

### Task 3: Complete producer-backed provenance verification

**Files:**
- Modify: `scripts/audit_pn2d_minimal6_fixed_state.py`
- Test: `tests/regression/test_pn2d_minimal6_fixed_state_audit.py`

- [x] Verify executable/source/identity/argument/input/output hashes before `out.mkdir`.
- [x] Ensure error messages distinguish missing producer from forged hash where applicable.
- [x] Run all provenance mutation tests with the real producer available.

### Task 4: Run end-to-end Task 5 acceptance

**Files:**
- Create: `.superpowers/sdd/task-5-report.md`
- Generate ignored: `build-release/reference_tcad/pn2d_minimal6_synthetic_audit_task5_final_20260713/`

- [x] Run all 33 focused tests.
- [x] Run synthetic CLI and verify 36/54/24 rows plus 14 figures.
- [x] Run Tasks 1–3 regression groups and `git diff --check`.
- [x] Write the Task 5 report with RED/GREEN/replay/CLI evidence and PASS/PASS conclusion.
