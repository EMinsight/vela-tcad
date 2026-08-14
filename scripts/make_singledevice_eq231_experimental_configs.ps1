param(
    [string]$Root = "build-release/reference_tcad/singledevice_sentaurus2018/vela_import_fixedmaterials/vela"
)

$variants = @(
    @{ Name = "p1"; Discretization = "p1_direct"; OxideBoundary = "none" },
    @{ Name = "cvfem"; Discretization = "cvfem_full"; OxideBoundary = "none" },
    @{ Name = "cvfem_wkb"; Discretization = "cvfem_full"; OxideBoundary = "devsim_wkb" }
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
        $config.sweep.write_state_file =
            "${endpoint}_eq231_$($variant.Name)_state.csv"
        $target = Join-Path $Root "${endpoint}_eq231_$($variant.Name).json"
        $config | ConvertTo-Json -Depth 100 | Set-Content -Encoding UTF8 $target
        Write-Output $target
    }
}
