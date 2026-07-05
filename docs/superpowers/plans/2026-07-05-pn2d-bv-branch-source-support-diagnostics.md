# PN2D BV Branch And Source-Support Diagnostics Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Determine whether the remaining PN2D BV Vela/Sentaurus high-bias mismatch is primarily caused by quasi-Fermi branch/state offsets, contact/neutral-region carrier-state drift, or carrier-row avalanche source-support feedback.

**Architecture:** Do not change solver physics first. Reuse existing Sentaurus exports, Vela VTK states, and Python diagnostics to build a bias-resolved evidence chain over `-15 V..-20 V`; only after the diagnostics identify a dominant mechanism should a follow-up implementation plan change C++ formulas or boundary/source mapping code.

**Tech Stack:** C++20, CMake/Ninja, MSYS2 UCRT64, Python standard library diagnostics, existing `build-release` PN2D reference artifacts, Sentaurus multibias TDR exports, and `vela_example_runner` output VTK files.

---

## Files And Responsibilities

- Use: `scripts/export_pn2d_bv_multibias_fields.py` to export missing high-bias Sentaurus field directories.
- Use: `scripts/diagnose_pn2d_bv_absolute_branch_offsets.py` to quantify `psi/phin/phip`, `psi-phin`, `phip-psi`, and density offsets by node class.
- Use: `scripts/diagnose_pn2d_bv_transition_edge_state.py` for selected active-path edge-level branch/drop comparisons.
- Use: `scripts/diagnose_pn2d_bv_sg_avalanche_edges.py` to generate per-bias Vela SG edge avalanche source dumps.
- Use: `scripts/diagnose_pn2d_bv_active_support_continuity_balance.py` to summarize transport, SRH, avalanche, and Sentaurus generation terms on active support nodes.
- Use: `scripts/diagnose_pn2d_bv_active_support_residual_proxy.py` to replay qF shifts and source substitutions against carrier-row residual proxies.
- Use: `scripts/diagnose_pn2d_bv_source_ownership_replay.py` only if an edge-local replay CSV is already available for the same support window.
- Modify: `docs/validation/pn2d_bv_validation.md` with a short result section after the diagnostics are run.

---

## Shared Paths

Use these PowerShell variables for every task:

```powershell
$env:Path = "D:\msys64\ucrt64\bin;D:\msys64\usr\bin;$env:Path"
$Root = "build-release\reference_tcad\pn2d_sentaurus2018"
$Vela = "$Root\vela"
$Sent = "$Root\sentaurus_multibias_highbias_20260705"
$Out = "$Root\reports\bv_branch_source_support_20260705"
$Mesh = "$Root\vela\mesh.json"
$Doping = "$Root\doping.csv"
$Support = "$Root\reports\sentaurus_default_bv_execution\thresholded_avalanche_support_m13p2\thresholded_avalanche_support_nodes.csv"
New-Item -ItemType Directory -Force $Out | Out-Null
```

If `$Mesh` does not exist, use `$Root\mesh.json`. Record the chosen path in the final report.

---

### Task 1: Prepare Comparable High-Bias Field Inputs

**Files:**
- Read: `build-release/reference_tcad/pn2d_sentaurus2018/sentaurus_vm_runs/pn2d_bv_codex_20260629/source`
- Create: `build-release/reference_tcad/pn2d_sentaurus2018/sentaurus_multibias_highbias_20260705`

- [ ] **Step 1: Confirm the Sentaurus TDR source contains the high-bias snapshots**

Run:

```powershell
Get-ChildItem "$Root\sentaurus_vm_runs\pn2d_bv_codex_20260629\source" -Filter "pn2d_bv_multibias_*_des.tdr" | Measure-Object
```

Expected: `Count` is `401`. If the count is lower than `401`, run the existing live Sentaurus VM refresh before continuing:

```powershell
python scripts\run_sentaurus_vm_reference.py pn2d --ssh-target sentaurus --source-dir reference_tcad\pn2d_sentaurus2018\source --local-output-dir "$Root\sentaurus_vm_runs" --remote-root ~/sentaurus_runs/vela_oracle --run-id pn2d_bv_branch_source_support_20260705 --stages bv
```

- [ ] **Step 2: Export Sentaurus high-bias field directories**

Run:

```powershell
python scripts\export_pn2d_bv_multibias_fields.py --source-dir "$Root\sentaurus_vm_runs\pn2d_bv_codex_20260629\source" --out-root $Sent --tdr-importer build-release\sentaurus_import.exe --start -15 --stop -20 --step -1 --final-bias -20 --intervals 400 --clean
```

Expected:

- `$Sent\sentaurus_-15v\fields` exists.
- `$Sent\sentaurus_-16v\fields` exists.
- `$Sent\sentaurus_-17v\fields` exists.
- `$Sent\sentaurus_-18v\fields` exists.
- `$Sent\sentaurus_-19v\fields` exists.
- `$Sent\sentaurus_-20v\fields` exists.
- `$Sent\export_summary.json` reports `count: 6`.

- [ ] **Step 3: Confirm exact-bias Vela VTK files exist**

Run:

```powershell
foreach ($b in @("-15","-16","-17","-18","-19","-20")) {
  Get-ChildItem $Vela -Filter "*.vtk" | Where-Object { $_.Name -match "_$b(?:\.0)?V\.vtk$" } | Select-Object -First 1 Name
}
```

Expected: one Vela VTK is printed for each bias. If any exact bias is missing, rerun the current `simulation_bv_minus20_avaljac.json` or the current promoted BV deck with VTK output enabled before continuing.

---

### Task 2: Quantify Branch And Carrier-State Offsets

**Files:**
- Read: `$Sent/sentaurus_<bias>v/fields`
- Read: `$Vela/*.vtk`
- Create: `$Out/absolute_branch_offsets`

- [ ] **Step 1: Run node-level branch offset summary**

Run:

```powershell
python scripts\diagnose_pn2d_bv_absolute_branch_offsets.py --mesh $Mesh --doping-csv $Doping --sentaurus-root $Sent --vela-vtk-root $Vela --out-dir "$Out\absolute_branch_offsets" --biases "-15,-16,-17,-18,-19,-20"
```

Expected:

- `$Out\absolute_branch_offsets\absolute_branch_offsets_nodes.csv`
- `$Out\absolute_branch_offsets\absolute_branch_offsets_summary.json`

- [ ] **Step 2: Classify the offset pattern**

Open `absolute_branch_offsets_summary.json` and record these facts in a scratch note:

- For `impact_active=true`, median and p90 of `delta_psi_minus_phin_V`.
- For `impact_active=true`, median and p90 of `delta_phip_minus_psi_V`.
- For `node_class=contact`, the same two metrics.
- For `doping_class=n` and `doping_class=p`, the same two metrics.
- For active nodes, `log10_vela_over_sentaurus_electron_density` and `log10_vela_over_sentaurus_hole_density`.

Decision rule:

- If active-node `delta_psi_minus_phin_V` or `delta_phip_minus_psi_V` grows monotonically from `-15 V` to `-20 V` and the density log ratio moves with it, classify the primary suspect as **branch/state feedback**.
- If contact nodes show a much larger offset than interior active nodes, classify the primary suspect as **contact carrier-state boundary anchoring**.
- If neutral n/p regions show large density offsets while active nodes do not, classify the primary suspect as **neutral-region carrier supply**.

- [ ] **Step 3: Add edge-level branch/drop evidence on the active path**

Use the active transition path CSV from the latest branch/localization report. Prefer:

```powershell
$PathEdges = "$Root\reports\sentaurus_default_bv_execution\focused_restart_m13p2_transition_edge_state\transition_path_edges.csv"
```

If that file does not exist, select the newest candidate and fail explicitly if none exists:

```powershell
if (-not (Test-Path $PathEdges)) {
  $PathEdges = Get-ChildItem "$Root\reports" -Recurse -Filter "*path*edge*.csv" -ErrorAction SilentlyContinue |
    Sort-Object LastWriteTime -Descending |
    Select-Object -First 1 -ExpandProperty FullName
}
if (-not $PathEdges) {
  throw "No active path edge CSV found under $Root\reports; run the transition-path localization diagnostic before this step."
}
```

Then run:

```powershell
python scripts\diagnose_pn2d_bv_transition_edge_state.py --path-edge-csv $PathEdges --mesh $Mesh --doping-csv $Doping --sentaurus-dir "$Sent\sentaurus_-20v" --vela-vtk-root $Vela --out-dir "$Out\transition_edge_state_m20" --bias -20 --carrier all
```

Expected:

- `$Out\transition_edge_state_m20\transition_edge_state.csv`
- Edge rows include `delta_phin_drop_V`, `delta_phip_drop_V`, `delta_electron_qf_field_V_m`, and density log ratios.

Decision rule:

- If endpoint absolute offsets are large but qF drops/fields are close, prioritize carrier density/branch level correction.
- If qF drops/fields are also far off on active edges, prioritize GradQF reconstruction/discretization.

---

### Task 3: Check Contact And Neutral-Region Carrier Supply

**Files:**
- Read: `$Out/absolute_branch_offsets/absolute_branch_offsets_nodes.csv`
- Create: `$Out/contact_neutral_state_summary.csv`

- [ ] **Step 1: Summarize offsets by node class and doping class**

Run:

```powershell
$Rows = Import-Csv "$Out\absolute_branch_offsets\absolute_branch_offsets_nodes.csv"
$Median = {
  param($values)
  $finite = @($values | Where-Object { $_ -ne $null -and $_ -ne "" } | ForEach-Object {[double]$_} | Sort-Object)
  if ($finite.Count -eq 0) { return "" }
  return $finite[[int](($finite.Count - 1) / 2)]
}
$Groups = $Rows | Group-Object bias_V,node_class,doping_class,impact_active
$Summary = foreach ($g in $Groups) {
  $items = $g.Group
  [pscustomobject]@{
    key = $g.Name
    count = $items.Count
    median_delta_psi_minus_phin_V = & $Median ($items | ForEach-Object {$_.delta_psi_minus_phin_V})
    median_delta_phip_minus_psi_V = & $Median ($items | ForEach-Object {$_.delta_phip_minus_psi_V})
    median_log10_n = & $Median ($items | ForEach-Object {$_.log10_vela_over_sentaurus_electron_density})
    median_log10_p = & $Median ($items | ForEach-Object {$_.log10_vela_over_sentaurus_hole_density})
  }
}
$Summary | Export-Csv "$Out\contact_neutral_state_summary.csv" -NoTypeInformation
```

Expected:

- `$Out\contact_neutral_state_summary.csv`

- [ ] **Step 2: Decide whether the carrier supply error is boundary-local or bulk-neutral**

Decision rule:

- If contact rows have a larger median carrier-density deficit than adjacent non-contact rows at the same bias, inspect ohmic contact phin/phip enforcement and minority relaxation knobs next.
- If n-neutral or p-neutral rows have the largest deficit while contact rows look aligned, inspect mobility/BGN/ni density reconstruction and neutral branch evolution.
- If only `impact_active=true` rows are deficient, focus on source feedback rather than contact supply.

---

### Task 4: Replay Source Injection And Carrier-Row Support

**Files:**
- Read: `$Support`
- Read: `$Sent/sentaurus_<bias>v`
- Read: `$Vela/*.vtk`
- Create: `$Out/sg_edges`
- Create: `$Out/active_support_balance_<bias>`
- Create: `$Out/active_support_residual_proxy_<bias>`

- [ ] **Step 1: Generate SG avalanche edge dumps at each high-bias point**

Run:

```powershell
New-Item -ItemType Directory -Force "$Out\sg_edges" | Out-Null
foreach ($b in @(-15,-16,-17,-18,-19,-20)) {
  $token = ("{0:g}" -f $b).Replace("-", "m").Replace(".", "p")
  $vtk = Get-ChildItem $Vela -Filter "*.vtk" | Where-Object { $_.Name -match "_$b(?:\.0)?V\.vtk$" } | Sort-Object LastWriteTime -Descending | Select-Object -First 1
  python scripts\diagnose_pn2d_bv_sg_avalanche_edges.py --vtk $vtk.FullName --mesh $Mesh --doping-csv $Doping --out-dir "$Out\sg_edges\bias_$token" --bias $b --top 2000 --material-ni-m3 1.4638914958767616e16 --bandgap-narrowing old_slotboom
}
```

Expected: each `$Out\sg_edges\bias_m*` directory contains `sg_avalanche_edges.csv` and `sg_avalanche_edge_summary.json`.

- [ ] **Step 2: Run active-support continuity balance at `-20 V`**

Run:

```powershell
python scripts\diagnose_pn2d_bv_active_support_continuity_balance.py --support-csv $Support --sg-edge-csv "$Out\sg_edges\bias_m20\sg_avalanche_edges.csv" --mesh $Mesh --doping-csv $Doping --sentaurus-dir "$Sent\sentaurus_-20v" --vela-vtk-root $Vela --out-dir "$Out\active_support_balance_m20" --bias -20 --material-ni-m3 1.4638914958767616e16 --bandgap-narrowing old_slotboom
```

Expected:

- `$Out\active_support_balance_m20\active_support_continuity_balance_nodes.csv`
- `$Out\active_support_balance_m20\active_support_continuity_balance_summary.json`

- [ ] **Step 3: Replay residual proxy with qF shifts matching the observed branch offsets**

Use the median active-support offsets from Task 2. If the offsets are close to the previously observed values, start with electron `+0.047 V` and hole `+0.056 V`:

```powershell
python scripts\diagnose_pn2d_bv_active_support_residual_proxy.py --support-csv $Support --sg-edge-csv "$Out\sg_edges\bias_m20\sg_avalanche_edges.csv" --mesh $Mesh --doping-csv $Doping --sentaurus-dir "$Sent\sentaurus_-20v" --vela-vtk-root $Vela --out-dir "$Out\active_support_residual_proxy_m20_qfshift" --bias -20 --material-ni-m3 1.4638914958767616e16 --bandgap-narrowing old_slotboom --transport-modes density_sg --electron-qf-shift-v 0.047 --hole-qf-shift-v 0.056 --qf-shift-scope support_nodes
```

Expected:

- `$Out\active_support_residual_proxy_m20_qfshift\active_support_residual_proxy_nodes.csv`
- `$Out\active_support_residual_proxy_m20_qfshift\active_support_residual_proxy_summary.json`

Decision rule:

- If qF shift replay substantially closes the carrier-row residual/source deficit, branch/state feedback is the dominant cause.
- If qF shift replay leaves the residual deficit mostly unchanged while source terms remain low, source support/current ownership is the dominant cause.

- [ ] **Step 4: Run source ownership replay only when a matching edge-local replay CSV exists**

Find a matching edge-local replay CSV:

```powershell
$EdgeReplayCsv = Get-ChildItem "$Root\reports" -Recurse -Filter "edge_local_source_replay*.csv" -ErrorAction SilentlyContinue |
  Where-Object { $_.FullName -match "m13p2|minus20|-20" } |
  Sort-Object LastWriteTime -Descending |
  Select-Object -First 1 -ExpandProperty FullName
$EdgeReplayCsv
```

If a CSV exists for the same `-20 V` state and `thresholded_avalanche_support_m13p2` support window, run:

```powershell
if ($EdgeReplayCsv) {
  python scripts\diagnose_pn2d_bv_source_ownership_replay.py --edge-replay-csv $EdgeReplayCsv --variant vela_psi_sentaurus_qf --bias -20 --support-class overlap --out-dir "$Out\source_ownership_replay_m20"
} else {
  Write-Host "No matching edge-local replay CSV found; skip source ownership replay and record the skip reason."
}
```

Expected:

- `$Out\source_ownership_replay_m20\source_ownership_replay_summary.csv`

Decision rule:

- If `vela_double_combined` or a carrier-split policy beats the production policy by a large margin, plan a follow-up source mapping experiment.
- If all ownership variants are poor, source ownership is not enough; return to branch/contact state.

---

### Task 5: Write The Evidence Summary

**Files:**
- Modify: `docs/validation/pn2d_bv_validation.md`

- [ ] **Step 1: Add a short section to the validation doc**

Append a section named `High-Bias Branch And Source-Support Diagnostic Pass` near the latest BV diagnostics. The section must cite these artifact roots:

- `build-release/reference_tcad/pn2d_sentaurus2018/sentaurus_multibias_highbias_20260705/`
- `build-release/reference_tcad/pn2d_sentaurus2018/reports/bv_branch_source_support_20260705/`

The section must contain one paragraph each for branch/state offsets, contact/neutral carrier-state evidence, source-support replay evidence, and final classification.

- [ ] **Step 2: Fill in the three finding bullets with numbers**

Required numbers:

- Active-node median `delta_psi_minus_phin_V` and `delta_phip_minus_psi_V` at `-20 V`.
- Active-node median electron/hole density log10 Vela/Sentaurus ratio at `-20 V`.
- Contact-node median offsets at `-20 V`.
- `active_support_residual_proxy_summary.json` before/after qF-shift implication.
- Any source ownership policy improvement if Task 4 Step 4 ran.

- [ ] **Step 3: Commit the documentation-only diagnostic record**

Run:

```powershell
git add docs/validation/pn2d_bv_validation.md
git commit -m "Document PN2D BV branch and source-support diagnostics"
```

---

## Stop Criteria

Stop and write a follow-up implementation plan when one of these is true:

- Branch/state offsets explain at least `0.5 decade` of the `-20 V` current gap through density/source replay.
- Contact rows show a distinct qF/carrier-density offset not present in interior active rows.
- Source ownership/support replay identifies a policy with a clearly better residual/source match than production.
- Required Sentaurus high-bias exports are missing and cannot be regenerated in the current environment.

Do not change Van Overstraeten coefficients, `alpha |J| / q` scaling, or SG current signs in this plan unless the diagnostics contradict the current formula audit.
