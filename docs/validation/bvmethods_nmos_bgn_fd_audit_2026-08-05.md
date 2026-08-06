# BVmethods NMOS OldSlotboom/Fermi-Dirac audit (2026-08-05)

## Scope

This audit isolates the high-doping material/statistics chain used by the
Sentaurus 2018 `BVmethods` NMOS case: silicon intrinsic density, OldSlotboom
bandgap narrowing, Fermi-Dirac carrier density, and ideal Ohmic-contact charge
neutrality.  Mobility and avalanche parameters are held fixed.

## Root cause and implementation

The Sentaurus Device user guide states that simulations combining `Fermi` with
`EffectiveIntrinsicDensity(OldSlotboom)` apply an additional 300 K correction
to bandgap narrowing by default; `EffectiveIntrinsicDensity(NoFermi)` disables
it.  The correction is the sum, for donors and acceptors, of the difference
between the Fermi-Dirac and Maxwell-Boltzmann reduced Fermi energies:

`DeltaEg_Fermi = Vt300 * sum(F_1/2^-1(N/Ndos) - ln(N/Ndos))`.

Vela previously implemented the OldSlotboom term and Fermi density separately
but omitted this combined-model correction.  The new explicit configuration
`bandgap_narrowing.fermi_statistics_correction` adds it without changing
historical decks.  The BVmethods replay enables it and uses the extracted
Sentaurus silicon intrinsic density `1.4638914958767616e10 cm^-3`.

## Fixed-state closure

At the main source/drain contact (`Nd=5.1e20 cm^-3`, `Na=1e18 cm^-3`):

- raw OldSlotboom narrowing: `0.1539642324 eV`;
- Fermi-statistics correction: `0.1392000359 eV`;
- corrected total: `0.2931642683 eV`;
- value inferred independently from Sentaurus fields: `0.2949909650 eV`;
- corrected density prediction / Sentaurus density: `0.9937728522`.

Across the audited drain-contact nodes the corrected density ratio is
`0.99377` to `0.99516`; at the p-type substrate it is `1.01062`.

## Self-consistent 6.4 V replay

The corrected replay converged in five Newton iterations:

- drain current: `9.668456986e-6 A/um`;
- integrated avalanche current: `8.851535625e-6 A/um`;
- `Iava/Id = 0.9155065424`;
- n+ contact potential error: `+0.91335 mV` (previously about `+70.5 mV`);
- substrate potential error: `-0.27273 mV`.

The material/statistics correction therefore closes the ideal Ohmic-contact
built-in potential, but does not close current-type IIC.  On high-field edge
2226, the electron QF-drop, mobility, generalized-Einstein, and geometry ratios
remain near unity (`1.0045`, `1.0635`, `0.9995`, `1.0`), while the electron
density and current ratios remain `0.4592` and `0.5765`.  The next root-cause
target is the internal absolute quasi-Fermi branch and carrier population, not
OldSlotboom, mobility, or avalanche integration.

## Reproducibility

- Audit script: `scripts/audit_bvmethods_nmos_bgn_fd.py`
- Replay script: `scripts/run_bvmethods_nmos_iic_current_support.py`
- Replay output: `build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela_validation/btbt_e2_iic_fermi_bgn_fixed6p4_20260805`
- Contact audit: replay output under `bgn_fd_audit`
- Source ledger: replay output under `edge_ledger`
- SG transport audit: replay output under `sg_transport`
