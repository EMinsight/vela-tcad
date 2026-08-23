# TransportModels DG band-drive audit

The Sentaurus base affinity is reconstructed as
`ElectronAffinity - 0.5 * BandgapNarrowing` at the Vg=1 V, Vd=2 V DG state.

| Region | Material | Sentaurus base affinity (eV) | Vela base affinity (eV) | Vela−Sentaurus (mV) | BGN range (eV) |
|---|---|---:|---:|---:|---:|
| R.Gateox | SiO2 | 0.9 | 0.95 | 50.000000 | 0–0 |
| R.PolyReox | SiO2 | 0.9 | 0.95 | 50.000000 | 0–0 |
| R.PolyReox_mirrored | SiO2 | 0.9 | 0.95 | 50.000000 | 0–0 |
| R.Substrate | Silicon | 4.07274038462 | 4.05 | -22.740385 | 0.0500034–0.362215 |
| R.Polygate | PolySilicon | 4.07274038462 | 4.05 | -22.740385 | 0.155223–0.155223 |
| R.Spacer | Nitride | 1.9 | 1.9 | -0.000000 | 0–0 |
| R.Spacer_mirrored | Nitride | 1.9 | 1.9 | -0.000000 | 0–0 |
