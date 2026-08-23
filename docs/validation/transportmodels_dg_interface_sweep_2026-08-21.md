# TransportModels DG Si/SiO2 interface comparison

Fixed-state audit at `Vg=1 V`, `Vd=2 V`. All three cases use the same mesh,
hybrid Vela-DD/Sentaurus-Q state, corrected material contract, and `sentaurus_box`
operator. The runner exit code `1` is expected because the diagnostic run is
deliberately limited to one inner and one outer iteration.

| Interface contract | Global L1 | Ratio to best | Max residual | Interface-cell L1 |
|---|---:|---:|---:|---:|
| Neutral continuous | 405275 | 1.0000 | 3391.67 | 52951.4 |
| Half-jump only | 405794 | 1.0013 | 3391.67 | 53243.9 |
| Full affine calibrated | 406365 | 1.0027 | 3391.67 | 53738.3 |

## Interpretation

- Best integrated residual: **Neutral continuous**.
- Best worst-node residual: **Neutral continuous**.
- The affine coefficients originated from the older SingleDevice mesh and are
  treated as a transferability experiment, not as universal physical constants.
- A model is eligible for the self-consistent stage only if it improves the global
  metric without creating a materially worse local hotspot.

Figure: `D:\code-repo\vela-tcad\build-release\reference_tcad\transportmodels_sentaurus2022\vela_baseline\dg_interface_fixed_state_sweep_2026-08-21\interface_comparison.png`
CSV: `D:\code-repo\vela-tcad\build-release\reference_tcad\transportmodels_sentaurus2022\vela_baseline\dg_interface_fixed_state_sweep_2026-08-21\interface_comparison.csv`
