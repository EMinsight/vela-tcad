# TransportModels DD deep-off self-consistent SRH fix

Bias: `Vg=-1 V`, `Vd=1.1 V`.

Applied changes: Sentaurus silicon intrinsic density `1.4638914958767616e10 cm^-3`, Fermi-corrected OldSlotboom BGN, and SRH `Nref=1e16 cm^-3` in Vela internal units.

| Metric | Result |
|---|---:|
| Vela Id (A/um) | 1.664369200e-15 |
| Sentaurus Id (A/um) | 1.634684064e-15 |
| Relative error | 1.8160% |
| Log-current error | 0.007816 dex |
| Newton iterations | 0 |
| Carrier-row violations | 0 |
| KCL residual (A/um) | 1.644936545e-18 |
| Numerical status | resolved |

A stricter `stall_residual_floor=1e-13` control retained global continuity closure and zero carrier-row violations, but its first Newton step failed with `line_search_non_decrease`. This bounds the remaining issue to nonlinear resolution below the accepted deep-off floor rather than the SRH material model.
