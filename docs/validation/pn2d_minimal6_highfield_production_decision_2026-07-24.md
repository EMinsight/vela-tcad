# PN2D Minimal6 high-field production decision

## Decision ledger

| Evidence | Result | Authorized action |
|---|---|---|
| Native element electric-field mobility replay | median below `0.001 dex`, all 320 element-carrier samples pass | retain evidence |
| Box-edge current replay | P95 below `0.0014 dex`, signs 100%, terminal error below `0.143%` | retain evidence |
| Fixed-state KCL after coefficient replacement | bounded failure, maximum `8.731080e-4` | require self-consistent solve |
| Factor attribution | `interaction_dominant`; element/global support is material | do not hide support difference in a fitted field scale |
| Imported-state residual | element-electric improves all 80 state/carrier pairs | causal residual evidence |
| Executable first step | global electric-field config improves all 48 paired updates | authorize config-only candidate sweep |
| Jacobian and boundary | `2.989497e-9`; boundary difference `0.0` | candidate is numerically admissible |
| SRH/source units | isolated factor audit passes; forward IV is insensitive | retain source-unit factor independently |

## Production action

Task 8 is authorized to run a separate, self-consistent 40-state candidate
with

`solver.mobility.high_field_driving_force = "electric_field"`.

This is a comparison-configuration change only.  It is not a change to the
Masetti/high-field formula, Scharfetter-Gummel current, or impact-ionization
formula.

The executable configuration still uses one global edge mobility, whereas
the best Sentaurus reconstruction uses coefficient-weighted native element
mobility.  Consequently, a successful candidate sweep can validate the drive
choice but cannot by itself prove element-support equivalence.  Any
element-supported production operator must be proposed as a separate
discretization design.

No fitted mobility, field scale, saturation velocity, or edge coefficient is
authorized.
