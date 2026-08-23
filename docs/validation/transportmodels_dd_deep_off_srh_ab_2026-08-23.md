# TransportModels DD deep-off SRH model A/B

Bias: `Vg=-1 V`, `Vd=1.1 V`; exact Sentaurus psi, quasi-Fermi potentials, n, and p are held fixed.

| Variant | SRH integral (A/um) | Sentaurus ratio |
|---|---:|---:|
| fermi_old_slotboom | -1.054907008e-15 | 64.9453% |
| boltzmann_old_slotboom | -1.054908320e-15 | 64.9453% |
| fermi_no_bgn | -3.341043151e-16 | 20.5691% |
| boltzmann_no_bgn | -3.341044699e-16 | 20.5691% |
| fermi_old_slotboom_fermi_correction | -1.103251907e-15 | 67.9216% |
| fermi_corrected_bgn_sentaurus_ni | -1.656835045e-15 | 102.0029% |

## Isolated effects

- Generalized Fermi factors relative to Boltzmann: -0.0001%.
- OldSlotboom relative to no BGN: 215.7418%.
- Fermi BGN correction relative to OldSlotboom: 4.5829%.
- Sentaurus silicon ni relative to 1e10 cm^-3: 50.1774%.
- Sentaurus SRH-field integral versus terminal-current closure: 0.6351%.

Raw summary: `D:\code-repo\vela-tcad\build-release\reference_tcad\transportmodels_sentaurus2022\reports\idvg_dd_deep_off_fixed_state_20260823\srh_model_ab\summary.json`
