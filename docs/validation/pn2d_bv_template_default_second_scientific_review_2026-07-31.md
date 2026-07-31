# PN2D BV template-default second independent scientific review

Date: 2026-07-31

Verdict: `REJECT_DEFAULT_CHANGE`

## Scope

This read-only review examined the frozen M0/M2 prospective contract, both
independent runs at each level, same-grid Sentaurus-golden parity, closure,
state and IV determinism, and the unified machine decision.

## Findings

1. M0 strongly supports the SG/Laux operator on the same grid:
   `V_break` differs by 0.001 V, the effective curve maximum error is
   0.003749 dex, and closure ratios remain below 7.14e-6.
2. M2 remains directionally consistent: `V_break` differs by 0.011 V, the
   effective curve maximum error is 0.07883 dex, and closure passes.
3. M2 nevertheless fails the frozen knee contract.  The knee median error is
   0.05798 dex against 0.05, and neither curve has a required `V_slope`
   crossing in the frozen window.
4. IIC is incomplete in both M0 and M2.  Because IIC is postprocess-only, this
   is not direct evidence of avalanche source-feedback failure, but it blocks
   the predeclared causal-observation contract.
5. Cross-mesh convergence is observation-only and neither rescues nor adds to
   the same-grid decision.

## Required action

Restore the PN2D BV template default to `legacy_cell_reconstructed`; retain
SG/Laux as an explicit opt-in.  Do not change the global C++ default.  Preserve
all rejected-run evidence.
