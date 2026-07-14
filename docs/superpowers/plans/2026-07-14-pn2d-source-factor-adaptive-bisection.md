# PN2D Source-Volume-Factor Adaptive Bisection Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reproducibly bisect the PN2D diagnostic source-volume-factor stability boundary to width `0.001953125`, retain failure evidence, and publish the physical interpretation without changing defaults.

**Architecture:** Add one Python orchestration script built around pure deck-generation, curve-summary, and bracket-update functions. The script clones the existing `factor_0p921875.json` baseline, runs at most three candidates sequentially through the release runner, writes ignored artifacts, and uses the existing knee-shape formulas for comparison. Unit tests exercise invariants and bisection logic without running the expensive solver.

**Tech Stack:** Python 3 standard library, existing PN2D diagnostic modules, C++20 release `vela_example_runner.exe`, `unittest`, MSYS2 UCRT64.

## Global Constraints

- Change only `solver.impact_ionization.source_volume_factor` and candidate-specific output paths in generated decks.
- Keep `quasi_fermi_update_limit_V=0.05`, secant predictor, terminal-current gate, and electron p95 branch limit `2.0 dex` unchanged.
- Initial bracket is stable `0.921875`, unstable `0.9375`; stop at width `<=0.001953125` after at most three new runs.
- A candidate is stable only when an accepted converged row reaches exactly `-20 V` within the existing `0.01 V` reporting tolerance.
- Generated simulation artifacts remain under `build-release/` and are not committed.
- Do not change production defaults or claim a validated BV curve.

---

### Task 1: Pure Deck And Bracket Logic

**Files:**
- Create: `scripts/run_pn2d_bv_source_factor_bisection.py`
- Create: `tests/regression/test_pn2d_bv_source_factor_bisection.py`

**Interfaces:**
- Consumes: baseline JSON dictionary and candidate factor.
- Produces: `build_candidate_deck(base, factor, label) -> dict`, `candidate_label(factor) -> str`, and `update_bracket(lower, upper, factor, stable) -> tuple[float, float]`.

- [ ] **Step 1: Write failing pure-function tests**

Test that the first three midpoint labels are exact, stable/unstable bracket updates are monotone, and a generated deck differs from the baseline only at:

```text
solver.impact_ionization.source_volume_factor
output_csv
sweep.write_state_file
```

- [ ] **Step 2: Run the focused tests and observe RED**

Run:

```powershell
python -m unittest tests.regression.test_pn2d_bv_source_factor_bisection -v
```

Expected: import failure because the orchestration script does not exist.

- [ ] **Step 3: Implement minimal pure functions**

Use `Decimal(str(value))` for midpoint and label construction so the sequence is exactly `0.9296875`, then the midpoint of the updated bracket. Deep-copy with JSON round-trip, update only the factor and relative output filenames, and reject an invalid or nonshrinking bracket.

- [ ] **Step 4: Run focused tests and observe GREEN**

Expected: all pure-function tests pass.

- [ ] **Step 5: Commit Task 1**

```powershell
git add scripts/run_pn2d_bv_source_factor_bisection.py tests/regression/test_pn2d_bv_source_factor_bisection.py
git commit -m "Add PN2D source factor bisection driver"
```

### Task 2: Solver Execution And Evidence Summaries

**Files:**
- Modify: `scripts/run_pn2d_bv_source_factor_bisection.py`
- Modify: `tests/regression/test_pn2d_bv_source_factor_bisection.py`

**Interfaces:**
- Consumes: runner path, baseline deck, Sentaurus reference curve, output directory, and bracket limits.
- Produces: per-candidate `run_summary.json`, top-level `bisection_summary.json`, `bisection_summary.csv`, and `bisection_summary.md`.

- [ ] **Step 1: Write failing execution-summary tests**

Create synthetic CSV fixtures in a temporary directory and verify:

```python
assert summarize_candidate(csv_reaching_minus20)["stable"] is True
assert summarize_candidate(csv_stopping_early)["stable"] is False
assert summary["failure_reason"] == "line_search_non_decrease"
```

Also test `--dry-run` writes three candidate-independent invariants without invoking a solver.

- [ ] **Step 2: Run focused tests and observe RED**

Expected: missing summary and CLI functions.

- [ ] **Step 3: Implement solver and summary flow**

Run:

```text
vela_example_runner.exe --config <candidate.json>
```

with candidate directory as `cwd`, capture full stdout/stderr to files, parse converged/failure rows, calculate knee markers and maximum `-20..-10 V` log10 current error using `diagnose_pn2d_bv_knee_shape`, copy any emitted `newton_failure_diagnostics.json` path into the summary, and update the bracket only after evidence is complete.

- [ ] **Step 4: Run focused and neighboring analysis tests**

```powershell
python -m unittest `
  tests.regression.test_pn2d_bv_source_factor_bisection `
  tests.regression.test_pn2d_bv_compensated_sg_replay_orchestrator -v
```

Expected: all tests pass.

- [ ] **Step 5: Commit Task 2**

```powershell
git add scripts/run_pn2d_bv_source_factor_bisection.py tests/regression/test_pn2d_bv_source_factor_bisection.py
git commit -m "Record PN2D source factor bisection evidence"
```

### Task 3: Execute Three-Step Adaptive Bisection

**Files:**
- Generate only: `build-release/reference_tcad/pn2d_sentaurus2018/reports/qflim0p05_source_factor_bisection_20260714/`

**Interfaces:**
- Consumes: committed orchestration script, release runner, `factor_0p921875.json`, and Sentaurus reference CSV.
- Produces: final bracket and candidate solver/failure evidence.

- [ ] **Step 1: Verify the release runner and baseline inputs exist**

```powershell
Test-Path build-release/vela_example_runner.exe
Test-Path build-release/reference_tcad/pn2d_sentaurus2018/reports/qflim0p05_high_source_scan/factor_0p921875.json
Test-Path build-release/reference_tcad/pn2d_sentaurus2018/reference_curves/pn2d_sentaurus2018_bv_reference.csv
```

Expected: all `True`.

- [ ] **Step 2: Run the adaptive bisection**

```powershell
python scripts/run_pn2d_bv_source_factor_bisection.py `
  --out-dir build-release/reference_tcad/pn2d_sentaurus2018/reports/qflim0p05_source_factor_bisection_20260714
```

Expected: at most three candidates and final bracket width `<=0.001953125`. Solver nonzero exit is recorded as an unstable physical result and does not abort evidence generation.

- [ ] **Step 3: Run branch-divergence comparison**

Invoke `scripts/diagnose_pn2d_bv_factor_branch_divergence.py` with the stable baseline and every newly stable curve in increasing-factor order. If no new candidate reaches `-20 V`, compare the first failed transition against `factor_0p9375_newton_failure_diagnostics.json` in the final report instead.

- [ ] **Step 4: Audit evidence consistency**

Verify candidate factors equal bracket midpoints, every bracket update follows the recorded stable flag, and no generated deck changed a fixed baseline field.

### Task 4: Document And Verify The Physical Conclusion

**Files:**
- Modify: `docs/validation/pn2d_bv_validation.md`
- Modify: `docs/validation/pn2d_bv_current_progress_summary.md`
- Test: `tests/regression/test_pn2d_bv_source_factor_bisection.py`

**Interfaces:**
- Consumes: final bisection summary, knee summaries, and Newton failure diagnostics.
- Produces: durable current status and next-step decision.

- [ ] **Step 1: Update validation documents**

Record every candidate factor, stability outcome, deepest bias, failure reason, current error, knee markers, final bracket, and failure classification. State that this remains a diagnostic override and not a validated physical BV curve.

- [ ] **Step 2: Run final verification**

```powershell
python -m unittest `
  tests.regression.test_pn2d_bv_source_factor_bisection `
  tests.regression.test_pn2d_bv_compensated_sg_replay_orchestrator -v
git diff --check
```

Expected: all tests pass and `git diff --check` is silent.

- [ ] **Step 3: Commit documentation**

```powershell
git add docs/validation/pn2d_bv_validation.md docs/validation/pn2d_bv_current_progress_summary.md
git commit -m "Document PN2D source factor stability boundary"
```
