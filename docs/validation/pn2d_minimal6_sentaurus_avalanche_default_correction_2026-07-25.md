# PN2D Minimal6 Sentaurus avalanche default correction

Date: 2026-07-25

Status: complete; independently verified; production formulas unchanged.

## Corrected outcome

Sentaurus Device O-2018.06-SP2 uses `GradQuasiFermi` as the default
driving-force selector for drift-diffusion avalanche generation. The
apparently electric-field-driven default observed in the earlier Minimal6
fixed-state replay is a mesh-specific contact fallback:

- by default, Sentaurus replaces the quasi-Fermi-potential gradient with the
  electric field in every element touching a contact;
- all four Minimal6 triangles touch either the Anode or the Cathode;
- consequently, implicit default, explicit `GradQuasiFermi`, and explicit
  `ElectricField` are exactly equivalent on this mesh unless
  `ComputeGradQuasiFermiAtContacts=UseQuasiFermi` is set.

The earlier electric-field fixed-state coefficient and source replay remains
numerically valid for the effective Minimal6 default. It must not be cited as
evidence that the global Sentaurus avalanche default is `ElectricField`.

## Manual contract

The local Sentaurus Device User Guide P-2019.03 establishes:

1. page 425: the default avalanche model is van Overstraeten-de Man and the
   default driving force is `GradQuasiFermi`;
2. page 380: inside elements touching a contact, the electric field replaces
   the quasi-Fermi gradient by default; setting
   `ComputeGradQuasiFermiAtContacts=UseQuasiFermi` requests the QFP gradient
   in every element;
3. page 426: the default avalanche current density is reconstructed from
   element-edge Scharfetter-Gummel currents; `AvalDensGradQF` selects the
   alternative `-q mu n grad(Phi)` current-density approximation;
4. page 439: `GradQuasiFermi` is the normal drift-diffusion avalanche
   default, while `ElectricField` is the plain-field option used for
   ionization-integral breakdown analysis.

Source manual:

`D:\工作\学习资料\TCAD软件手册\Sentaurus PDFManual 2019\data\sdevice_ug.pdf`

## Mesh contact-fallback proof

The four physical triangles have vertex sets:

| Element | Physical vertices | Contact touched |
| ---: | --- | --- |
| 0 | 0, 4, 5 | Anode through 0 and 4 |
| 1 | 0, 5, 1 | Anode through 0 |
| 2 | 1, 5, 2 | Cathode through 2 |
| 3 | 5, 3, 2 | Cathode through 2 and 3 |

There is no interior triangle on which the default `GradQuasiFermi` selector
can avoid the contact electric-field replacement.

## Experiment A: effective default

The first remote matrix contains two topologies, three biases
`-1/-10/-20 V`, and four branches:

- implicit `Avalanche(VanOverstraeten)`;
- explicit `Avalanche(VanOverstraeten GradQuasiFermi)`;
- explicit `Avalanche(VanOverstraeten ElectricField)`;
- explicit `GradQuasiFermi` plus `Math { AvalDensGradQF }`.

Across every parsed node, element, element-vertex, local edge, runtime
integral, CurrentPlot integral, and terminal current:

| Candidate relative to explicit GradQuasiFermi | Exact parsed match | Key result |
| --- | --- | --- |
| implicit default | yes | confirms default selector |
| explicit ElectricField | yes | confirms all-element contact fallback |
| GradQuasiFermi plus AvalDensGradQF | no | changes avalanche current-density support |

The first two equivalences are bitwise at the scientific `.plt` level for
both topologies. For the `AvalDensGradQF` control, the maximum changes over
the six states are:

| Quantity | Maximum difference |
| --- | ---: |
| potential | `6.53089e-6 V` |
| carrier density | `1.09714e-4 dex` |
| mobility | `9.64e-17 dex` |
| avalanche coefficient | `1.59e-15 dex` |
| node generation | `0.454116 dex` |
| device source integral | `0.360933 dex` |
| terminal current | `1.96947e-3` relative |

The zero coefficient and mobility differences, together with the source
change, confirm that `AvalDensGradQF` changes the current-density
approximation in the avalanche source. The small high-bias state and current
changes are nonlinear feedback from that source change.

## Experiment B: force true QFP gradients at contacts

The second matrix adds:

- `GradQuasiFermi` plus
  `ComputeGradQuasiFermiAtContacts=UseQuasiFermi`;
- the same branch plus `AvalDensGradQF`.

Sentaurus logs explicitly report both the Math parameter and
`driving force is Gradient Quasi Fermi` for electrons and holes.

Relative to the normal contact-fallback default, the forced-QFP branch gives:

| Bias | Potential max (V) | Density max (dex) | Mobility max (dex) | Alpha max (dex) | Generation max (dex) | Source-integral max (dex) | Terminal current relative |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| -10 | `3.06079e-3` | `0.0514191` | `0.0518441` | `2.65793` | `2.65693` | `1.18450` | `6.66779e-3` |
| -20 | `1.73772e-3` | `0.0291925` | `0.0293699` | `0.667478` | `0.667200` | `0.541665` | `5.38012e-3` |

This branch proves that the contact fallback is active, but it is not a pure
avalanche-drive experiment: the same Math setting also changes the default
high-field mobility drive in contact elements.

## Experiment C: mobility-isolated avalanche drive

The causal isolation disables high-field saturation mobility in both
branches, retains the same doping-dependent low-field mobility, sets
`ComputeGradQuasiFermiAtContacts=UseQuasiFermi` in both, and changes only:

- reference: `Avalanche(VanOverstraeten ElectricField)`;
- candidate: `Avalanche(VanOverstraeten GradQuasiFermi)`.

The generated decks are otherwise identical after output-name
normalization. Runtime logs confirm no high-field mobility model and the
declared avalanche driving force for each carrier.

The mirror and sketch results agree under their verified cell permutation:

| Bias | Mobility error (dex) | Alpha max (dex) | Generation max (dex) | Total source GradQF / E-field (dex) | Terminal current relative |
| ---: | ---: | ---: | ---: | ---: | ---: |
| -1 | `0` | `9.27086` | `48.5594` | `-4.97695` | `0` |
| -10 | `0` | `2.73497` | `2.73497` | `-0.564589` | `2.57967e-8` |
| -20 | `0` | `0.694370` | `0.694468` | `-0.227260` | `6.64978e-4` |

The `-1 V` ratios compare values near numerical zero and are retained only
as low-signal diagnostics. The `-10 V` and `-20 V` rows establish that the
avalanche driving force alone materially changes alpha and the avalanche
source while mobility is exactly unchanged.

At `-20 V`, the carrier-resolved device integrals are:

| Carrier | ElectricField source (A/um) | GradQuasiFermi source (A/um) | GradQF / E-field (dex) |
| --- | ---: | ---: | ---: |
| electron | `2.129501e-19` | `1.262123e-19` | `-0.227176` |
| hole | `1.283703e-22` | `5.163196e-23` | `-0.395546` |
| total | `2.130785e-19` | `1.262640e-19` | `-0.227260` |

## Scientific correction and production decision

The corrected dependency statement is:

1. Sentaurus global drift-diffusion avalanche default:
   `GradQuasiFermi`;
2. Sentaurus effective Minimal6 avalanche drive:
   electric field in all four cells because every cell touches a contact;
3. Sentaurus default avalanche current support:
   element-edge Scharfetter-Gummel reconstruction;
4. `AvalDensGradQF`:
   an alternative current-density approximation, not a request to change the
   avalanche coefficient driving force;
5. the earlier Vela electric-field fixed-state replay:
   valid for this Minimal6 target, but not general-mesh evidence.

No production formula or default is changed by this correction. The
`element_edge_sg_gss_laux` operator remains opt-in. A general-mesh production
decision requires at least one interior element not touching a contact and
must compare the Sentaurus default contact rule rather than globally forcing
either electric field or QFP gradient.

## Deterministic evidence

- Default/contact-fallback raw root:
  `build-release/pn2d-minimal6-sentaurus-avalanche-drive-controls-20260725`
- Forced-QFP contact raw root:
  `build-release/pn2d-minimal6-sentaurus-avalanche-contact-controls-20260725`
- Corrected combined comparison:
  `build-release/pn2d-minimal6-sentaurus-avalanche-corrected-comparison-20260725`
- Corrected combined independent verification:
  `build-release/pn2d-minimal6-sentaurus-avalanche-corrected-verification-20260725`
- Mobility-isolation raw root:
  `build-release/pn2d-minimal6-sentaurus-avalanche-mobility-isolation-20260725`
- Mobility-isolation comparison:
  `build-release/pn2d-minimal6-sentaurus-avalanche-mobility-isolation-comparison-20260725`
- Mobility-isolation independent verification:
  `build-release/pn2d-minimal6-sentaurus-avalanche-mobility-isolation-verification-20260725`

The corrected combined verifier checked 30 state comparisons and 12,180
individual physical-quantity comparisons. The mobility-isolation verifier
checked 6 state comparisons and 2,436 individual comparisons. Both report
`status: independently_verified`.
