# Sentaurus 2018 BVmethods NMOS reference

This directory preserves the text inputs used by the Sentaurus Training
`BVmethods` NMOS example and the method mapping used by Vela validation. Binary
TDR/PLT/log outputs remain generated artifacts under `build*/reference_tcad/`
and are not checked in.

## Source inventory

| Checked-in file | Training workbench node | Method |
|---|---:|---|
| `source/bvmethods_nmos_sde.cmd` | 1 | symmetric NMOS structure and mesh |
| `source/aba_poisson_sdevice.cmd` | 3 | ABA, Poisson-only sweep |
| `source/aba_coupled_sdevice.cmd` | 4 | ABA postprocessing on coupled DD |
| `source/external_resistor_sdevice.cmd` | 5 | `1e7 ohm*um` series-resistor load line |
| `source/voltage_to_current_sdevice.cmd` | 6 | voltage ramp to 6 V, then current control |
| `source/continuation_sdevice.cmd` | 7 | Sentaurus continuation curve tracing |
| `source/transient_sdevice.cmd` | 8 | transient voltage ramp with resistor |
| `source/models.par` | 3--8 | shared SRH lifetime parameters |

The SDevice files retain the original contact names, physical models, target
values, step controls and current direction. File/output names were made
descriptive so the decks can coexist outside Sentaurus Workbench.

Run the structure deck first. It builds `bvmethods_nmos_half_msh.tdr`, mirrors
it about `x=0`, renames the mirrored drain contact to source, and writes
`bvmethods_nmos_msh.tdr`. Each SDevice deck expects that final mesh in its
working directory.

## Vela template mapping

- `configs/templates/bvmethods_nmos_external_resistor.template.json` maps the
  series-resistor load-line method.
- `configs/templates/bvmethods_nmos_voltage_to_current.template.json` maps the
  voltage-to-current boundary switch.
- IIC/ABA validation remains represented by the existing BV/IIC configuration
  workflow and the checked-in validation scripts.
- Sentaurus `Continuation` and `Transient` are archived as references only;
  Vela does not yet expose equivalent pseudo-arclength or transient BV modes.

The validated Sentaurus threshold voltages are `6.379791636 V` for the
external-resistor method and `6.383184201 V` for voltage-to-current control at
the example drain-current criterion. No mobility or avalanche coefficient
scaling is used by the Vela templates.
