# PN2D Active-Region Branch/State Feedback Debug Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Find why Vela's active-region carrier branch/state is driven away from the Sentaurus avalanche support branch in high-bias PN2D BV.

**Architecture:** Treat the discrepancy as a feedback-loop problem, not a scalar avalanche-fit problem. The plan builds reproducible diagnostics that follow the chain `state -> SG flux -> alpha*|J| source -> continuity residual -> Newton update -> next state` on the same Sentaurus-active support edges.

**Tech Stack:** C++20 Vela release runner, Python diagnostic scripts under `scripts/`, existing coarse7x3 Sentaurus/Vela artifacts, CSV/VTK field exports, Catch2 plus Python unit tests where practical.

---

## Current Evidence To Preserve

- Release `psi_gradient_proxy` 1 V grid reached `0..-20 V` with 21 converged points.
- Local Sentaurus field exports currently cover `-18 V` and `-20 V`; `-19 V` is current-curve-only unless a new Sentaurus field export is produced.
- On Sentaurus-active top support, Vela `alpha*flux` is nearly absent:
  - `-18 V`: Vela/Sentaurus `alpha*flux sum = 3.226301e-72`
  - `-20 V`: Vela/Sentaurus `alpha*flux sum = 5.468599e-26`
- On the same support, branch/state deltas are large:
  - `-18 V`: `d(psi-phin)=-0.54104 V`, `log10(n_Vela/n_Sentaurus)=-9.08907`
  - `-20 V`: `d(psi-phin)=-0.438162 V`, `log10(n_Vela/n_Sentaurus)=-7.35997`

Do not overwrite these artifacts:

- `build-release/reference_tcad/pn2d_sentaurus2018_coarse7x3/reports/coarse_psi_gradient_proxy_1vgrid_20260705/`
- `build-release/reference_tcad/pn2d_sentaurus2018_coarse7x3/reports/coarse_vm_vector_compare/`

## File Structure

- Create: `scripts/diagnose_pn2d_active_region_branch_feedback.py`
  - One CLI that reads existing Vela/Sentaurus CSV/VTK diagnostics and emits normalized support, state, flux, source, residual-proxy, and rollback tables.
- Modify: `tests/test_reference_tcad_tools.py`
  - Add focused fixture tests for support selection, branch-delta math, and alpha-flux aggregation.
- Optional modify: `src/solver/NewtonSolver.cpp`
  - Only if existing probes cannot expose the needed first-step carrier-row deltas.
- Optional modify: `src/simulation/DCSweep.cpp`
  - Only if we need a default-off diagnostic dump for per-edge/per-node residual ownership during a real sweep.
- Create output directory per run:
  - `build-release/reference_tcad/pn2d_sentaurus2018_coarse7x3/reports/active_region_branch_feedback_YYYYMMDD/`

## Task 1: Freeze Inputs And Recompute The Baseline Tables

**Files:**
- Create output only under `build-release/reference_tcad/pn2d_sentaurus2018_coarse7x3/reports/active_region_branch_feedback_YYYYMMDD/`

- [ ] **Step 1: Verify exact input artifacts exist**

Run:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
Test-Path build-release\reference_tcad\pn2d_sentaurus2018_coarse7x3\reports\coarse_psi_gradient_proxy_1vgrid_20260705\coarse_psi_gradient_proxy_1vgrid_20260705.csv
Test-Path build-release\reference_tcad\pn2d_sentaurus2018_coarse7x3\reports\coarse_vm_vector_compare\sentaurus_multibias\sentaurus_-18v\field_manifest.json
Test-Path build-release\reference_tcad\pn2d_sentaurus2018_coarse7x3\reports\coarse_vm_vector_compare\sentaurus_multibias\sentaurus_-20v\field_manifest.json
```

Expected: all three commands print `True`.

- [ ] **Step 2: Record the artifact manifest**

Write `manifest.json` containing:

```json
{
  "vela_case": "coarse_psi_gradient_proxy_1vgrid_20260705",
  "sentaurus_field_biases": [-18.0, -20.0],
  "current_only_biases": [-19.0],
  "support_basis": "Sentaurus alpha*|J| on same imported mesh edges",
  "current_approximation": "psi_gradient_proxy"
}
```

- [ ] **Step 3: Copy only small summary CSVs**

Copy these existing CSVs into the new report directory with names prefixed `input_`:

```text
coarse_psi_gradient_proxy_1vgrid_20260705.csv
coarse_node_field_compare_psi_gradient_proxy_1vgrid.csv
coarse_current_support_compare_psi_gradient_proxy_1vgrid.csv
sentaurus_active_edge_branch_state_summary_psi_gradient_proxy_1vgrid.csv
```

Expected: copies are byte-identical to source files.

## Task 2: Add A Reusable Support/State Diagnostic Script

**Files:**
- Create: `scripts/diagnose_pn2d_active_region_branch_feedback.py`
- Test: `tests/test_reference_tcad_tools.py`

- [ ] **Step 1: Write failing tests for pure math helpers**

Add tests for:

```python
def branch_delta(row):
    return {
        "delta_psi_minus_phin": (row["vela_psi"] - row["vela_phin"]) - (row["sent_psi"] - row["sent_phin"]),
        "delta_phip_minus_psi": (row["vela_phip"] - row["vela_psi"]) - (row["sent_phip"] - row["sent_psi"]),
    }

def alpha_flux(e_flux, h_flux, e_alpha, h_alpha):
    return e_flux * e_alpha + h_flux * h_alpha
```

Test cases:

```python
def test_branch_delta_uses_absolute_branch_offsets():
    row = {
        "vela_psi": 1.0,
        "vela_phin": 0.7,
        "vela_phip": 1.4,
        "sent_psi": 1.0,
        "sent_phin": 0.9,
        "sent_phip": 1.2,
    }
    assert branch_delta(row)["delta_psi_minus_phin"] == pytest.approx(0.2)
    assert branch_delta(row)["delta_phip_minus_psi"] == pytest.approx(0.2)

def test_alpha_flux_combines_electron_and_hole_support():
    assert alpha_flux(2.0, 3.0, 5.0, 7.0) == pytest.approx(31.0)
```

- [ ] **Step 2: Run tests and confirm failure**

Run:

```powershell
python -m pytest tests\test_reference_tcad_tools.py -k "branch_delta or alpha_flux" -q
```

Expected: fail because helpers do not exist.

- [ ] **Step 3: Implement the helper functions**

Create helper functions in `scripts/diagnose_pn2d_active_region_branch_feedback.py`:

```python
def branch_delta(vela_psi: float, vela_phin: float, vela_phip: float,
                 sent_psi: float, sent_phin: float, sent_phip: float) -> dict[str, float]:
    return {
        "delta_psi_minus_phin_V": (vela_psi - vela_phin) - (sent_psi - sent_phin),
        "delta_phip_minus_psi_V": (vela_phip - vela_psi) - (sent_phip - sent_psi),
    }

def alpha_flux(electron_flux: float, hole_flux: float,
               electron_alpha: float, hole_alpha: float) -> float:
    return electron_flux * electron_alpha + hole_flux * hole_alpha
```

- [ ] **Step 4: Run tests and confirm pass**

Run:

```powershell
python -m pytest tests\test_reference_tcad_tools.py -k "branch_delta or alpha_flux" -q
```

Expected: pass.

## Task 3: Build Three Comparable Support Sets

**Files:**
- Modify: `scripts/diagnose_pn2d_active_region_branch_feedback.py`

- [ ] **Step 1: Add support selectors**

Implement CLI options:

```text
--support sentaurus_top
--support vela_top
--support overlap
--top 20
--biases -18,-20
```

Definitions:

- `sentaurus_top`: rank edges by Sentaurus `alpha_e*|Je|/q + alpha_h*|Jh|/q`.
- `vela_top`: rank edges by Vela `electron_alpha_m_inv*electron_flux_abs + hole_alpha_m_inv*hole_flux_abs`.
- `overlap`: intersection of top `N` Sentaurus and top `N` Vela edge ids at the same bias.

- [ ] **Step 2: Emit support membership CSV**

Output: `support_membership.csv`

Required columns:

```text
bias_V,support,rank,edge_id,node0,node1,edge_class,
sent_alpha_flux,vela_alpha_flux,vela_over_sent_alpha_flux,
electron_flux_over_sentaurus,hole_flux_over_sentaurus,
electron_alpha_over_sentaurus,hole_alpha_over_sentaurus
```

- [ ] **Step 3: Validate with existing evidence**

Run:

```powershell
python scripts\diagnose_pn2d_active_region_branch_feedback.py `
  --node-field-csv build-release\reference_tcad\pn2d_sentaurus2018_coarse7x3\reports\coarse_psi_gradient_proxy_1vgrid_20260705\coarse_node_field_compare_psi_gradient_proxy_1vgrid.csv `
  --support-csv build-release\reference_tcad\pn2d_sentaurus2018_coarse7x3\reports\coarse_psi_gradient_proxy_1vgrid_20260705\coarse_current_support_compare_psi_gradient_proxy_1vgrid.csv `
  --out-dir build-release\reference_tcad\pn2d_sentaurus2018_coarse7x3\reports\active_region_branch_feedback_20260705 `
  --biases -18,-20 `
  --top 20
```

Expected:

- `sentaurus_top` `-20 V` Vela/Sentaurus alpha-flux sum remains near `5.47e-26`.
- `sentaurus_top` `-18 V` Vela/Sentaurus alpha-flux sum remains near `3.23e-72`.

## Task 4: Quantify Branch/State On Each Support Set

**Files:**
- Modify: `scripts/diagnose_pn2d_active_region_branch_feedback.py`

- [ ] **Step 1: Add endpoint state aggregation**

For each selected edge, collect endpoint nodes and compute medians:

```text
delta_psi_minus_phin_V
delta_phip_minus_psi_V
log10_n_vela_over_sentaurus
log10_p_vela_over_sentaurus
```

- [ ] **Step 2: Emit `support_state_summary.csv`**

Required columns:

```text
bias_V,support,edge_count,node_count,
median_delta_psi_minus_phin_V,median_delta_phip_minus_psi_V,
median_log10_n_vela_over_sentaurus,median_log10_p_vela_over_sentaurus,
median_psi_diff_V,median_phin_diff_V,median_phip_diff_V
```

- [ ] **Step 3: Classify branch failure mode**

Add a textual `classification` field:

- `electron_branch_deficit`: `median_delta_psi_minus_phin_V < -0.05` and `median_log10_n < -0.5`
- `hole_branch_deficit`: `median_delta_phip_minus_psi_V > 0.05` and `median_log10_p < -0.5`
- `support_mismatch`: Sentaurus top and Vela top edge overlap fraction `< 0.25`
- `alpha_zeroed`: median Vela/Sentaurus alpha ratio `< 1e-6`
- `flux_zeroed`: median Vela/Sentaurus flux ratio `< 1e-6`

Expected for current evidence: Sentaurus-active support should classify at least `electron_branch_deficit`, `alpha_zeroed`, or `flux_zeroed` depending on bias.

## Task 5: Separate Alpha Failure From Flux/Carrier Failure

**Files:**
- Modify: `scripts/diagnose_pn2d_active_region_branch_feedback.py`

- [ ] **Step 1: Add four replay ratios**

For each Sentaurus-active support edge, compute:

```text
baseline_vela = alpha_vela * flux_vela
sent_alpha_vela_flux = alpha_sentaurus * flux_vela
vela_alpha_sent_flux = alpha_vela * flux_sentaurus
sent_full = alpha_sentaurus * flux_sentaurus
```

Use electron + hole combined values.

- [ ] **Step 2: Emit `support_replay_decomposition.csv`**

Required columns:

```text
bias_V,edge_id,
baseline_over_sent_full,
sent_alpha_vela_flux_over_sent_full,
vela_alpha_sent_flux_over_sent_full,
sent_full_reference,
limiting_factor
```

Classification rule:

- If `sent_alpha_vela_flux_over_sent_full << 1`, Vela flux/carrier state is insufficient even with Sentaurus alpha.
- If `vela_alpha_sent_flux_over_sent_full << 1`, Vela alpha/drive is insufficient even with Sentaurus current.
- If both are small, the branch has moved both state and driving force away from Sentaurus.

- [ ] **Step 3: Verify against `psi_gradient_proxy`**

Expected: both Vela alpha and Vela flux are very small on Sentaurus-active support, which means replacing current-density qF gradient with `psi` gradient is not a sufficient correction.

## Task 6: Trace The Feedback Through Continuity Residual Terms

**Files:**
- Prefer existing probe if available:
  - `scripts/diagnose_pn2d_bv_predictor_carrier_row_audit.py`
  - `scripts/diagnose_pn2d_bv_predictor_first_step_audit.py`
- Create wrapper if needed:
  - `scripts/diagnose_pn2d_active_region_feedback_rows.py`

- [ ] **Step 1: Create a Sentaurus-active mixed-state predictor**

For `-20 V`, create two initial states from Vela converged state:

1. `qf_only`: replace only `phin/phip` on Sentaurus-active endpoint nodes with Sentaurus values; reconstruct `n/p`.
2. `qf_plus_flux_support`: same state, plus diagnostic replay using Sentaurus active edge support for flux/source decomposition. Do not change production residual yet.

Output:

```text
predictor_qf_only_state.csv
predictor_qf_only_manifest.json
```

- [ ] **Step 2: Run first Newton step probe**

Use existing `newton_step_probe` or add a default-off runner mode only if the existing tool cannot consume the state CSV.

Required output:

```text
first_step_update_by_node.csv
```

Columns:

```text
node_id,support_class,
delta_phin_intended,delta_phin_newton,retained_phin_fraction,
delta_phip_intended,delta_phip_newton,retained_phip_fraction,
electron_residual_before,hole_residual_before,
electron_flux_term,hole_flux_term,
electron_impact_term,hole_impact_term
```

- [ ] **Step 3: Decide if rollback is residual-driven or globalization-driven**

Interpretation:

- If first full Newton step and carrier-only step both roll back the state, the carrier residual/Jacobian is the driver.
- If carrier-only keeps the state but full Newton rolls back, Poisson/coupling is the driver.
- If raw Newton keeps the state but line search rejects it, globalization/branch guard is the driver.

Expected from prior evidence: carrier-row flux/source balance likely treats Sentaurus-aligned branch as over-fluxed and rolls it back.

## Task 7: Inspect Carrier-Row Jacobian Coefficients On Active Endpoints

**Files:**
- Optional modify: `src/solver/NewtonSolver.cpp`
- Optional create: `scripts/diagnose_pn2d_active_region_jacobian_rows.py`

- [ ] **Step 1: Dump local row coefficients**

For active endpoint carrier rows, dump:

```text
row_id,node_id,carrier,
diag,abs_offdiag_sum,diag_dominance_ratio,
d_flux_d_self,d_flux_d_neighbor,
d_impact_d_self,d_impact_d_neighbor,
rhs,row_update
```

- [ ] **Step 2: Compare Vela baseline vs Sentaurus-like predictor**

Run at `-20 V`:

```powershell
python scripts\diagnose_pn2d_active_region_jacobian_rows.py `
  --baseline-state build-release\reference_tcad\pn2d_sentaurus2018_coarse7x3\reports\coarse_psi_gradient_proxy_1vgrid_20260705\coarse_psi_gradient_proxy_1vgrid_20260705_last_state.csv `
  --predictor-state build-release\reference_tcad\pn2d_sentaurus2018_coarse7x3\reports\active_region_branch_feedback_20260705\predictor_qf_only_state.csv `
  --support-membership build-release\reference_tcad\pn2d_sentaurus2018_coarse7x3\reports\active_region_branch_feedback_20260705\support_membership.csv `
  --out-dir build-release\reference_tcad\pn2d_sentaurus2018_coarse7x3\reports\active_region_branch_feedback_20260705\jacobian_rows
```

Expected: identify whether the active carrier rows are dominated by SG flux increase, missing impact feedback, or diagonal scaling.

## Task 8: Minimal Hypothesis Tests

**Files:**
- Diagnostic-only config files under the report directory.
- Do not modify production defaults.

- [ ] **Hypothesis A: Vela alpha/drive is the primary limiter**

Test: replay with Sentaurus alpha and Vela flux.

Pass condition for hypothesis: `sent_alpha_vela_flux_over_sent_full` moves close to `1`.

Fail condition: remains far below `1`.

- [ ] **Hypothesis B: Vela carrier/current state is the primary limiter**

Test: replay with Vela alpha and Sentaurus flux.

Pass condition for hypothesis: `vela_alpha_sent_flux_over_sent_full` moves close to `1`.

Fail condition: remains far below `1`.

- [ ] **Hypothesis C: Coupled carrier residual rolls back Sentaurus-like state**

Test: first Newton step on `qf_only` predictor.

Pass condition for hypothesis: retained QF shift fraction is much less than `1`, especially on Sentaurus-active endpoints.

- [ ] **Hypothesis D: Source feedback is too weak relative to SG flux increase**

Test: carrier-row term decomposition shows `impact_term / flux_delta` is too small to stabilize Sentaurus-like state.

Pass condition for hypothesis: impact scale needed to close active carrier rows is consistently above `1.5x`.

## Task 9: Acceptance Criteria For This Debug Pass

The debug pass is complete when the report answers all of these:

- [ ] On Sentaurus-active support, is the main deficit alpha, flux/current, carrier density, or support mismatch?
- [ ] Does the first Newton step move a Sentaurus-like active state toward or away from the low-current Vela branch?
- [ ] Which residual block causes the movement: carrier flux, impact source, Poisson coupling, or line search/globalization?
- [ ] Is `psi_gradient_proxy` ruled out as a corrective formula, or does it only expose path sensitivity?
- [ ] Is the next implementation target clear enough to write a small default-off experiment?

## Task 10: Final Report

**Files:**
- Create: `build-release/reference_tcad/pn2d_sentaurus2018_coarse7x3/reports/active_region_branch_feedback_20260705/active_region_branch_feedback_summary.md`

Required sections:

```markdown
# Active-Region Branch/State Feedback Debug Summary

## Inputs
## Current Curve Context
## Sentaurus-Active Support Result
## Vela-Active Support Result
## Alpha vs Flux Replay
## First Newton Step Rollback
## Carrier-Row Term Balance
## Classification
## Next Default-Off Experiment
```

Final classification must use one of:

- `alpha_drive_primary`
- `carrier_flux_state_primary`
- `carrier_row_feedback_primary`
- `globalization_branch_selection_primary`
- `mixed_alpha_and_state`

For the current evidence, expect `carrier_row_feedback_primary` or `mixed_alpha_and_state`, but do not assert either until Tasks 5-7 are complete.

## Execution Notes

- Use release runner for BV reproduction:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
.\build-release\vela_example_runner.exe --config <config.json>
```

- Keep all new diagnostics default-off.
- Do not tune Van Overstraeten coefficients in this debug pass.
- Do not replace production current-density formula in this debug pass.
- Do not use `psi_gradient_proxy` as an acceptance candidate; use it only as a diagnostic branch perturbation.

