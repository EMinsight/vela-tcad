param(
    [string]$Root = "build-release/reference_tcad/singledevice_sentaurus2018/vela_import_fixedmaterials/vela"
)

$variants = @(
    @{ Name = "p1"; Discretization = "p1_direct"; OxideBoundary = "none"; TraceOffset = 0.0; SiWeight = 1.0; PolyWeight = 1.0; OxSiWeight = 1.0; OxPolyWeight = 1.0; InsulatorCornerWeight = 1.0 },
    @{ Name = "cvfem"; Discretization = "cvfem_full"; OxideBoundary = "none"; TraceOffset = 0.0; SiWeight = 1.0; PolyWeight = 1.0; OxSiWeight = 1.0; OxPolyWeight = 1.0; InsulatorCornerWeight = 1.0 },
    @{ Name = "p1_lambda"; Discretization = "p1_lambda_direct"; OxideBoundary = "none"; TraceOffset = 0.0; SiWeight = 1.0; PolyWeight = 1.0; OxSiWeight = 1.0; OxPolyWeight = 1.0; InsulatorCornerWeight = 1.0 },
    @{ Name = "gss_state"; Discretization = "gss_potentiallike_fitted"; OxideBoundary = "none"; TraceOffset = 0.0; SiWeight = 1.0; PolyWeight = 1.0; OxSiWeight = 1.0; OxPolyWeight = 1.0; InsulatorCornerWeight = 1.0 },
    @{ Name = "sentaurus_box"; Discretization = "sentaurus_box"; OxideBoundary = "none"; TraceOffset = 0.02012; SiWeight = 0.3613278292533479; PolyWeight = 1.0684933639683336; OxSiWeight = 2.6839079693374917; OxPolyWeight = 2.569027176700638; SiOffset = -0.00020247747279261268; PolyOffset = 0.01872581675079906; OxSiOffset = -0.0015039829729206406; OxPolyOffset = -0.0052046570150173915; InsulatorCornerWeight = 1.0 },
    @{ Name = "conservative_sqrt"; Discretization = "conservative_sqrt_fitted"; OxideBoundary = "none"; TraceOffset = 0.0; SiWeight = 1.0; PolyWeight = 1.0; OxSiWeight = 1.0; OxPolyWeight = 1.0; InsulatorCornerWeight = 1.0 },
    @{ Name = "gss_lambda"; Discretization = "gss_density_fitted"; OxideBoundary = "none"; TraceOffset = 0.0; SiWeight = 1.0; PolyWeight = 1.0; OxSiWeight = 1.0; OxPolyWeight = 1.0; InsulatorCornerWeight = 1.0 }
)

foreach ($endpoint in @("lin", "sat")) {
    $source = Join-Path $Root "${endpoint}_eq231_global_potentiallike_probe.json"
    $base = Get-Content -Raw $source | ConvertFrom-Json
    foreach ($variant in $variants) {
        $config = $base | ConvertTo-Json -Depth 100 | ConvertFrom-Json
        $quantum = $config.solver.electron_quantum_potential
        $quantum | Add-Member -NotePropertyName global_discretization `
            -NotePropertyValue $variant.Discretization -Force
        $quantum | Add-Member -NotePropertyName oxide_boundary `
            -NotePropertyValue $variant.OxideBoundary -Force
        $quantum | Add-Member `
            -NotePropertyName sentaurus_interface_insulator_half_jump_offset `
            -NotePropertyValue $variant.TraceOffset -Force
        $quantum | Add-Member -NotePropertyName `
            sentaurus_interface_silicon_half_jump_offset `
            -NotePropertyValue $(if ($variant.Name -eq "sentaurus_box") { -4.6008840569922854e-5 } else { 0.0 }) -Force
        $quantum | Add-Member -NotePropertyName `
            sentaurus_interface_polysilicon_half_jump_offset `
            -NotePropertyValue $(if ($variant.Name -eq "sentaurus_box") { 0.0026674992132365016 } else { 0.0 }) -Force
        $quantum | Add-Member -NotePropertyName `
            sentaurus_interface_silicon_reaction_weight `
            -NotePropertyValue $variant.SiWeight -Force
        $quantum | Add-Member -NotePropertyName `
            sentaurus_interface_polysilicon_reaction_weight `
            -NotePropertyValue $variant.PolyWeight -Force
        $quantum | Add-Member -NotePropertyName `
            sentaurus_interface_insulator_at_silicon_reaction_weight `
            -NotePropertyValue $variant.OxSiWeight -Force
        $quantum | Add-Member -NotePropertyName `
            sentaurus_interface_insulator_at_polysilicon_reaction_weight `
            -NotePropertyValue $variant.OxPolyWeight -Force
        $quantum | Add-Member -NotePropertyName `
            sentaurus_interface_silicon_reaction_offset_V `
            -NotePropertyValue $(if ($variant.ContainsKey("SiOffset")) { $variant.SiOffset } else { 0.0 }) -Force
        $quantum | Add-Member -NotePropertyName `
            sentaurus_interface_polysilicon_reaction_offset_V `
            -NotePropertyValue $(if ($variant.ContainsKey("PolyOffset")) { $variant.PolyOffset } else { 0.0 }) -Force
        $quantum | Add-Member -NotePropertyName `
            sentaurus_interface_insulator_at_silicon_reaction_offset_V `
            -NotePropertyValue $(if ($variant.ContainsKey("OxSiOffset")) { $variant.OxSiOffset } else { 0.0 }) -Force
        $quantum | Add-Member -NotePropertyName `
            sentaurus_interface_insulator_at_polysilicon_reaction_offset_V `
            -NotePropertyValue $(if ($variant.ContainsKey("OxPolyOffset")) { $variant.OxPolyOffset } else { 0.0 }) -Force
        $quantum | Add-Member -NotePropertyName `
            sentaurus_insulator_reentrant_corner_reaction_weight `
            -NotePropertyValue $variant.InsulatorCornerWeight -Force
        $quantum | Add-Member -NotePropertyName oxide_quantum_mass_ratio `
            -NotePropertyValue 0.14 -Force
        $quantum | Add-Member -NotePropertyName oxide_barrier_mass_ratio `
            -NotePropertyValue 0.4 -Force
        $quantum | Add-Member -NotePropertyName oxide_barrier_height_V `
            -NotePropertyValue 3.15 -Force
        $diagnosticRoot = Join-Path $Root `
            "reports/eq231_experimental_20260814/$($variant.Name)/$endpoint"
        $quantum.residual_diagnostic_prefix =
            [System.IO.Path]::GetFullPath($diagnosticRoot)
        $config.output_csv = "${endpoint}_eq231_$($variant.Name).csv"
        $config.log_file = "${endpoint}_eq231_$($variant.Name).log"
        # Both endpoint probes must consume the same checked-in material
        # contract.  Some older generated saturation probes referenced a
        # local "updated" copy and silently bypassed newly added parameters.
        $config.materials_file = "materials_sentaurus2018.json"
        $config.sweep.write_state_file =
            "${endpoint}_eq231_$($variant.Name)_state.csv"
        $target = Join-Path $Root "${endpoint}_eq231_$($variant.Name).json"
        $config | ConvertTo-Json -Depth 100 | Set-Content -Encoding UTF8 $target
        Write-Output $target

        if ($variant.Name -eq "sentaurus_box") {
            $curve = $config | ConvertTo-Json -Depth 100 | ConvertFrom-Json
            $curveQuantum = $curve.solver.electron_quantum_potential
            $curveQuantum.outer_max_iterations = 30
            $curveQuantum.max_iterations = 500
            $curveQuantum.max_update_V = 0.1
            $curveQuantum.damping = 0.5
            $curveQuantum.outer_acceleration = "none"
            $curveQuantum.outer_relaxation = 1.0
            $curveQuantum.outer_relaxation_min = 0.1
            $curveQuantum.outer_relaxation_max = 1.0
            $curveQuantum.relative_tolerance = 0.0
            $curveQuantum.absolute_tolerance_V = 1.0e-6
            $curveQuantum | Add-Member -NotePropertyName `
                outer_absolute_tolerance_V -NotePropertyValue 5.0e-4 -Force
            $curveQuantum.residual_diagnostic_prefix = ""
            $curveQuantum.residual_diagnostic_use_initial_state = $false
            $curve.sweep.start = -0.5
            $curve.sweep.stop = 2.2
            $curve.sweep.step = 0.135
            # Start from the durable low-gate-bias state.  The imported
            # Sentaurus endpoints are at Vg=2.2 V; solving the new nonlinear
            # interface closure from that high-inversion state excites an
            # unstable outer mode before continuation can begin.
            $curve.sweep.initial_state_file = $(
                if ($endpoint -eq "lin") {
                    "lin_self_consistent_final.csv"
                } else {
                    "sat_idvg_sentaurus_box_off_state.csv"
                })
            $curve.sweep.write_state_file =
                "${endpoint}_idvg_sentaurus_box_final_state.csv"
            $curve.output_csv = "${endpoint}_idvg_sentaurus_box_curve.csv"
            $curve.log_file = "${endpoint}_idvg_sentaurus_box_curve.log"
            if ($endpoint -eq "sat") {
                # The high-drain final frozen solve can reach a numerically
                # zero displayed residual and then reject another Newton step.
                # A 1e-6 relative tolerance accepts that already-converged DD
                # state without changing the 0.5 mV quantum outer criterion.
                $curve.solver.reltol = 1.0e-6
                $curve.solver.abstol = 1.0e-5
                # Start at half of the 21-point reference spacing through the
                # sensitive saturation subthreshold range.  The adaptive
                # sweep may grow back to the reference spacing once stable.
                $curve.sweep.step = 0.0675
                $curve.sweep.max_step = 0.135
            }
            $curveTarget = Join-Path $Root `
                "${endpoint}_idvg_sentaurus_box_curve.json"
            $curve | ConvertTo-Json -Depth 100 |
                Set-Content -Encoding UTF8 $curveTarget
            Write-Output $curveTarget

            if ($endpoint -eq "lin") {
                $offState = $curve | ConvertTo-Json -Depth 100 |
                    ConvertFrom-Json
                $offState.sweep.stop = -0.5
                $offState.sweep.write_state_file =
                    "lin_idvg_sentaurus_box_off_state.csv"
                $offState.output_csv =
                    "lin_idvg_sentaurus_box_off_state_point.csv"
                $offState.log_file =
                    "lin_idvg_sentaurus_box_off_state.log"
                $offTarget = Join-Path $Root `
                    "lin_idvg_sentaurus_box_off_state.json"
                $offState | ConvertTo-Json -Depth 100 |
                    Set-Content -Encoding UTF8 $offTarget
                Write-Output $offTarget
            } else {
                $drainRamp = $curve | ConvertTo-Json -Depth 100 |
                    ConvertFrom-Json
                ($drainRamp.contacts | Where-Object name -eq "gate").bias = -0.5
                $drainRamp.sweep.contact = "drain"
                $drainRamp.sweep.start = 0.1
                $drainRamp.sweep.stop = 1.1
                $drainRamp.sweep.step = 0.05
                $drainRamp.sweep.max_step = 0.05
                $drainRamp.sweep.initial_state_file =
                    "lin_idvg_sentaurus_box_off_state.csv"
                $drainRamp.sweep.write_state_file =
                    "sat_idvg_sentaurus_box_off_state.csv"
                $drainRamp.output_csv =
                    "sat_idvg_sentaurus_box_drain_ramp.csv"
                $drainRamp.log_file =
                    "sat_idvg_sentaurus_box_drain_ramp.log"
                $drainRampTarget = Join-Path $Root `
                    "sat_idvg_sentaurus_box_drain_ramp.json"
                $drainRamp | ConvertTo-Json -Depth 100 |
                    Set-Content -Encoding UTF8 $drainRampTarget
                Write-Output $drainRampTarget
            }
        }
    }
}
