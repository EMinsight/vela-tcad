# TransportModels DG fixed-state residual audit

Work point: Vg = 1.0 V, Vd = 2.0 V

Status: **pass**

## Fixed-state definition

The DD variables come from the converged Vela self-consistent 2 V state. The
electron quantum potential and continuous potential-like field come from the
Sentaurus 2022 DG state. The p1_direct Eq. 231 operator is assembled from this
initial state before applying a DG update.

## Main result

- Maximum free-node raw residual: `21010.7321813` at node
  `2` (-0.0667875 µm, 0.025 µm).
- Dominant component at that node: `reaction`.
- Global component L1 shares: stiffness 8.90%,
  gradient-squared 17.56%, reaction
  73.53%, explicit interface boundary
  0.00%.

## Region ranking

| Region | Material | Share of global free-node L1 | Interface share within region | Max cell |
|---|---|---:|---:|---:|
| R.Substrate | Si | 69.04% | 0.34% | 3245 |
| R.PolyReox | SiO2 | 9.90% | 70.43% | 670 |
| R.PolyReox_mirrored | SiO2 | 9.86% | 70.53% | 677 |
| R.Gateox | SiO2 | 6.43% | 21.93% | 114 |
| R.Spacer | Nitride | 2.03% | 25.29% | 6018 |

The raw residual is an integrated discrete equation value and is best used for
relative localization and component attribution. It is not a voltage error and
must not be compared directly with the Qn field error in mV.

## Provenance

- Config SHA-256: `BA1AB5586D48B0071031485BDE5FC6A655EF236FBB88A08B867D8A2180485E05`
- Fixed hybrid state SHA-256: `4665D7C01C4B87ABA3935B0E769909ACDA149ABBCA3E63D954DC7953971BB917`
- Node residual CSV SHA-256: `4B85F91EAAE1BABA1FD46C18AF0CA2626A4E6345B4AC52446340E6050D4373E5`
- Cell residual CSV SHA-256: `39DFE7EF11C6D5922781EFF8560E3A66AC9D7629B18A02C23E95BDB4CF1A3E04`
