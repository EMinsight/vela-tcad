# PN2D Minimal6 current-factor attribution

## Result

The fixed imported-state experiment evaluated all 2 x 3 x 2 combinations of

- low-field mobility source: Vela or Sentaurus;
- drive: global-edge QFP, triangle QFP, or native element electric field; and
- support: one global edge or coefficient-weighted native elements.

For the binary Vela/global-QFP/global-edge to
Sentaurus/electric-field/native-element replacement, all six replacement
orders were evaluated on every one of the 400 valid carrier-edge samples.
The 4,800 factorial samples and 7,200 paired increments close to a maximum
`2.220446e-16 dex`, below the `1e-12 dex` gate.

The independently recomputed baseline exactly reproduces the recorded Vela
imported-state production-current branch, and the factorial target exactly
reproduces `sentaurus_lowfield_element_electric_field`.

## Ordinary active edges

| Carrier | Factor | Median abs Shapley (dex) | P95 (dex) | Current-weighted mean abs (dex) |
|---|---|---:|---:|---:|
| electron | low-field coefficient | 0.0336100 | 0.0603004 | 0.0403244 |
| electron | drive | 0.0193801 | 0.100148 | 0.0498389 |
| electron | support | 0.0445718 | 0.0854419 | 0.0683708 |
| hole | low-field coefficient | 0.0280090 | 0.0588656 | 0.0400342 |
| hole | drive | 0.0163687 | 0.0649597 | 0.0386445 |
| hole | support | 0.0439517 | 0.0803378 | 0.0673242 |

The unweighted ranking is `support > low_field > drive`; the
reference-current-weighted ranking is `support > drive > low_field`.  The
combined interaction remainder is larger than every individual weighted
factor (electron `0.155344 dex`, hole `0.150234 dex`), so the typed outcome is
`interaction_dominant`.  The rank disagreement is retained explicitly.

The central `1-5` near-zero-current tail is reported separately and is not
used to establish the ordinary-edge ranking.  No mobility, field scale, or
edge coefficient was fitted.

## Determinism

The `20260724-a` and `20260724-b` roots are byte-identical for every generated
CSV, report, and independent verification JSON.  Both independent
verifications pass with zero failures.

No production formula was changed.
