# PN2D Minimal6 Sentaurus high-field driving-force control

Date: 2026-07-24

Status: `valid`

Primary outcome: `electric_field_best_supported_candidate`

Production formulas modified by this diagnostic: `false`

## Control design

Two isolated Sentaurus O-2018.06-SP2 control runs were generated, one for each
Minimal6 topology at -20 V. The control deck retained:

- `DopingDependence`;
- `SRH`;
- `Avalanche(VanOverstraeten)`;
- `EffectiveIntrinsicDensity(OldSlotboom)`; and
- the same mesh, material parameters, temperature, and contact bias.

The only mobility-model change was removal of `HighFieldSaturation`. Native
element electron and hole mobility were exported explicitly. Because the
remaining Masetti low-field branch depends on doping and temperature, its four
element values per topology were used as the Sentaurus-native low-field
coefficients for all 40 exact high-field states.

The documented Caughey-Thomas high-field law was then evaluated with three
candidate element driving fields:

1. exported native `eGradQuasiFermi/Element` or
   `hGradQuasiFermi/Element`;
2. the affine triangle gradient of exported node QFP; and
3. exported native `ElectricField/Element`.

## Native low-field coefficient comparison

| Carrier | N | Median Sentaurus/Vela low-field error (dex) | P95 (dex) | Maximum (dex) |
|---|---:|---:|---:|---:|
| electron | 160 | 0.139859 | 0.202939 | 0.202939 |
| hole | 160 | 0.0851601 | 0.123333 | 0.123333 |

This establishes that part of the earlier native element mobility residual is
already present before high-field saturation. It cannot be attributed solely
to the selected driving field.

## Direct final-mobility replay over all 40 states

| Carrier | Candidate drive | N | Median error (dex) | P95 (dex) | Maximum (dex) |
|---|---|---:|---:|---:|---:|
| electron | exported native QFP gradient | 160 | 0.497771 | 0.618552 | 0.641861 |
| electron | affine triangle QFP gradient | 160 | 0.0207907 | 0.107439 | 0.141608 |
| electron | native element electric field | 160 | 0.000833449 | 0.00130115 | 0.00134817 |
| hole | exported native QFP gradient | 160 | 0.0176517 | 0.0661291 | 0.0731585 |
| hole | affine triangle QFP gradient | 160 | 0.0176517 | 0.0661291 | 0.0731585 |
| hole | native element electric field | 160 | 0.000443255 | 0.00056937 | 0.000576072 |

## Inverted effective-field comparison

The final Sentaurus element mobility was inverted using the native low-field
coefficient and the same documented high-field exponent and saturation
velocity.

| Carrier | Field comparison | N | Median error (dex) | P95 (dex) | Maximum (dex) |
|---|---|---:|---:|---:|---:|
| electron | inverted versus native QFP gradient | 160 | 0.917208 | 1.10341 | 1.63502 |
| electron | inverted versus triangle QFP gradient | 160 | 0.0344636 | 0.220662 | 0.752279 |
| electron | inverted versus electric field | 160 | 0.00104502 | 0.00302827 | 0.00373147 |
| hole | inverted versus native QFP gradient | 160 | 0.0340003 | 0.215471 | 0.713104 |
| hole | inverted versus triangle QFP gradient | 160 | 0.0340003 | 0.215471 | 0.713104 |
| hole | inverted versus electric field | 160 | 0.000669557 | 0.00163904 | 0.00190564 |

## Interpretation for the current discrepancy

The strongest tested reconstruction of Sentaurus element mobility is:

`native Sentaurus low-field mobility + native element electric field`.

This result rejects treating exported electron
`eGradQuasiFermi/Element` as the internal high-field mobility drive. It also
shows that the affine triangle QFP gradient is a useful approximation but is
not the closest candidate once the native low-field coefficient is available.

The result refines, but does not overturn, the current-factor diagnosis:

- importing the Sentaurus potential/QFP state still removes about 0.72 dex of
  paired electron-current error and 0.85 dex of paired hole-current error;
- the self-consistent QFP state remains the dominant full-current discrepancy;
- native low-field mobility and element-to-edge support explain part of the
  remaining approximately 0.06 dex fixed-state current error; and
- all self-consistent sign mismatches remain localized to the very small
  central edge 1-5 current and follow the QFP differential-mode sign.

No SRH scaling change, continuity-source scaling change, production mobility
change, SG change, or QFP-equation change is justified by this control alone.

## Boundary of inference

The low-field control changes the solved carrier state, but the retained
Masetti coefficient is state-independent at fixed doping and 300 K. The
reconstruction is therefore a strong operator control, not a direct
observation of a proprietary internal Sentaurus option. A Vela production
change should first replay the complete box-edge current with native
low-field element mobility and electric-field saturation on the same imported
state.

## Evidence

- Evidence root A:
  `build-release/pn2d-minimal6-sentaurus-highfield-drive-20260724-a`
- Evidence root B:
  `build-release/pn2d-minimal6-sentaurus-highfield-drive-20260724-b`
- A/B directory diff: byte-identical.
- Carrier-element samples: 320.
- Task 3 native -20 V control runs: 2.
- Task 2 bias-invariance controls at -1, -10, and -20 V: 6.
- Bias-invariance roots: `build-release/pn2d-minimal6-lowfield-bias-invariance-20260724-a/b`.
- Sentaurus release: O-2018.06-SP2.
