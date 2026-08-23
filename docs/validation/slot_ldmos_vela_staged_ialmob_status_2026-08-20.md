# Slot-LDMOS Vela staged BVDS and Sentaurus IALMob status

Date: 2026-08-20

## Vela staged execution

| Stage | Result | Evidence |
|---|---|---|
| 00 equilibrium | Pass | Zero-bias state converged. |
| 01 1 Tohm unit control | Pass | At 1 V outer bias, inner voltage is 0.00837439824 V and drain current is 9.91625580e-13 A/um. |
| 02 avalanche off to 60 V | Pass | 15/15 targets converged; at 60 V outer bias, inner voltage is 0.897340756 V and drain current is 5.91026592e-11 A/um. |
| 03 IIC postprocess | Pass | At the 60 V load-line state, qG is 5.86051676e-11 A/um and maximum electric field is 3.46256e4 V/cm. |
| 04 self-consistent avalanche activation | Pass | The activated 1 V state converged; qG is 4.78650e-68 A/um. |
| 05 self-consistent avalanche to 60 V | Fail | First boundary probe at 0.0183743983 V fails at Newton iteration 2 with `line_search_non_decrease`; electron and hole residual blocks are 7.096/7.130. |
| 06 final external-resistor BVDS | Blocked | Requires a converged Stage 05 state. |

The Stage 05 highest Poisson-residual node is node 10237 at approximately
(0.018839, 0.366406) um. All incident triangles belong to the Silicon region;
it is not a PolySilicon or dielectric interface node.

## PolySilicon nontransport A/B control

PolySilicon is now electrostatic-only in the generated Vela bundle. Carrier
material parameters are omitted and exported node doping is retained only on
nodes incident to Silicon cells. The independent control bundle is
`build-release/reference_tcad/slot_ldmos_sentaurus2022/run01/vela_ready_poly_nontransport`.

| Metric | Original bundle | Poly nontransport | Relative change |
|---|---:|---:|---:|
| 60 V avalanche-off inner voltage (V) | 0.897341059 | 0.897340756 | -3.38e-7 |
| 60 V avalanche-off drain current (A/um) | 5.91026589e-11 | 5.91026592e-11 | +5.13e-9 |
| Stage 05 outcome | line-search failure | same line-search failure | no improvement |

Therefore PolySilicon transport semantics are corrected, but they do not cause
the leakage discrepancy or the Stage 05 nonlinear failure.

## Sentaurus IALMob on/off control

The no-IALMob controls were run with Sentaurus T-2022.03-SP2. Only
`Enormal(IALMob)` was removed. The process TDR, parameter file, high-field
mobility, SRH, Auger, avalanche, external resistor, bias schedule, and numerical
controls were unchanged. Both control decks returned exit code 0.

| Metric | IALMob on | IALMob off | Off - on |
|---|---:|---:|---:|
| 60 V avalanche-off drain current (A/um) | 6.587956014e-15 | 6.593721583e-15 | +0.087517% |
| BVDS at 1e-7 A/um (V) | 38.520901204 | 38.730622726 | +0.209721522 V (+0.544436%) |

IALMob has a measurable but sub-percent effect on BVDS. It cannot explain the
roughly 8.97e3 Vela/Sentaurus avalanche-off leakage ratio or the Vela Stage 05
line-search failure. A full IALMob implementation should therefore remain a
quantitative model-equivalence task, not the current nonlinear-solver blocker.

## Next development gate

The deterministic Stage 05 ablation at the first 0.0183743983 V inner-voltage
probe localized the failure to the finite-difference triangle-GSS avalanche
source Jacobian. Removing the avalanche residual converged, keeping the
self-consistent residual while freezing only its source Jacobian converged,
and bypassing Gummel still failed with the finite-difference source Jacobian.
This originally led the production Slot-LDMOS stages to use a self-consistent
residual with `source_jacobian: frozen`. A subsequent SDevice manual audit
corrected the interpretation: avalanche derivatives are enabled by default,
and the source deck does not contain `-AvalDerivatives`. The production
preparation now selects `source_jacobian: local_ad`; `frozen` is retained only
as an ablation. At the saved 15.7209505709 V state and hotspot node 10236, the
three-step JVP sweep reduced the maximum relative error from about 3.30e-2
(`frozen`) to 2.40e-5 (`local_ad`).
