# PN2D Minimal6 high-field residual and first-step follow-up

## Result

The mobility candidate is classified as `mobility_candidate_causal`.

The fixed-state element-electric reconstruction reduces the continuity
residual on all 80 topology/bias/carrier pairs.  Its median residual ratios
relative to Vela production are `0.869920` for electrons and `0.894376` for
holes.  The native-final element reference gives `0.868555` and `0.893641`;
the triangle-QFP replay gives `0.904986` and `0.925859`.  Constant mobility
increases the residual by median factors `7.22073` and `3.68784`.

The executable global-edge `electric_field` configuration was then compared
with the existing global-edge `quasi_fermi_gradient` configuration at
mirror/sketch x -1/-10/-20 V:

- carrier-only and fully coupled first steps were both evaluated;
- all 48 paired topology/bias/mode/carrier/node QFP updates decreased;
- electric/QFP absolute-update ratios range from `0.961340` to `0.990027`;
- all 90 analytic/finite-difference Jacobian blocks pass, with maximum
  relative difference `2.989497e-9`; and
- contact boundary rows remain identical, with maximum difference `0.0`.

The production edge-to-node residual replay closes at `3.308722e-24`.

## Scope

The element-electric, native-element, and triangle branches remain explicitly
typed as `box_operator_reconstruction`.  Their fixed-state residual evidence
must not be described as a native executable edge-mobility implementation.
The first-step and Jacobian evidence uses the existing executable global-edge
`electric_field` configuration.

All compared branches retain the same production source-unit snapshot.  No
mobility experiment was combined with an SRH or impact-source scaling change.

## Determinism

The `20260724-a` and `20260724-b` roots are byte-identical for every stable
CSV, report, and independent verification JSON.  Both independent
verifications pass with zero failures.

No production formula was changed.
