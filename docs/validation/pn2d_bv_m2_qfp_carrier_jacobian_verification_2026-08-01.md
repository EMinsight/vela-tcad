# PN2D M2 carrier-QFP residual and Jacobian verification

## Decision

The observation-only verification is complete.  The typed outcome is:

`phip_dominant__electron_flux_dominant__hole_flux_dominant__analytic_fd_inconsistent__both_qfp_updates_roll_back_from_sentaurus`

The formal all-block analytic/finite-difference gate remains failed, exactly as
declared, because the SRH/Auger relative comparison fails.  The seven-step
sensitivity result refines that outcome without changing the gate:

`formal_relative_gate_fails_only_at_srh_absolute_fd_floor`

The dynamically important Poisson, transport, SG-avalanche, and boundary/gauge
blocks pass the unchanged `5e-5` relative threshold.  Current evidence therefore
localizes the first meaningful mismatch to carrier transport response to QFP
inside the coupled solve; it does not identify a SG/Laux avalanche derivative
defect and does not authorize a production change.

## Frozen contract

- Mesh: balanced-junction M2, 115 Vela physical nodes.
- Biases: `-18`, `-19.5`, `-19.7`, and `-20 V`.
- States: Vela baseline, Sentaurus electron QFP only, Sentaurus hole QFP only,
  and both Sentaurus QFPs.
- Fixed-source operator: complete `element_edge_sg_gss_laux` with
  `element_vertex_box_measure`, `quasi_fermi_gradient`, and
  `postprocess_only` coupling.
- Residual terms: transport/SG flux, recombination, avalanche, gauge, and
  boundary, split by electron and hole equations.
- Jacobian states: Vela baseline and joint Sentaurus-QFP state.
- Main finite-difference step: `1e-7 V`, double symmetric.
- SRH sensitivity steps: `1e-5`, `3e-6`, `1e-6`, `3e-7`, `1e-7`, `3e-8`,
  and `1e-8 V` at `-19.5` and `-20 V`.
- Formal relative threshold: `5e-5`, unchanged.
- Two complete runs must be byte-identical.

No physical model, continuation schedule, production default, or acceptance
threshold was modified.

## Carrier-QFP source attribution

Hole QFP has the larger absolute frozen-source contribution at every declared
bias.  At `-20 V`, hole-QFP-only substitution recovers `0.991518` of the
joint-QFP log-source-error improvement, versus `0.859524` for electron QFP.

| Bias (V) | Electron-QFP source / Sentaurus | Electron recovery | Hole-QFP source / Sentaurus | Hole recovery |
|---:|---:|---:|---:|---:|
| -18.0 | 1.038639 | -0.552791 | 1.021354 | -0.883309 |
| -19.5 | 0.967082 | 1.609865 | 0.985487 | 1.938074 |
| -19.7 | 0.954803 | 1.192570 | 0.974784 | 1.444139 |
| -20.0 | 0.930731 | 0.859524 | 0.948076 | 0.991518 |

Recovery above one means that a one-carrier substitution overshoots the
improvement produced by replacing both QFPs; negative recovery means it moves
the source farther from Sentaurus.  This is attribution under a frozen state,
not proof that the hole continuity equation contains the implementation error.

## Carrier residual decomposition

At `-20 V`, joint-QFP substitution changes the interior residual mainly through
transport/SG flux:

| Carrier equation | Transport share | Avalanche share | Recombination share | Boundary/gauge share |
|---|---:|---:|---:|---:|
| Electron | 0.883898 | 0.116102 | 4.39e-8 | 0 |
| Hole | 0.892987 | 0.107013 | 4.05e-8 | 0 |

Transport wins at all four biases for both carriers.  The maximum absolute
term-sum closure error is `1.32349e-23`.

## Analytic versus finite-difference Jacobian

| Physical block | Worst relative difference | Bias (V) | State | Subblock | Formal gate |
|---|---:|---:|---|---|---|
| Poisson | 6.12734e-8 | -19.5 | baseline | phin column | pass |
| Transport | 5.15168e-8 | -18.0 | baseline | electron-phin | pass |
| SG avalanche | 9.05564e-10 | -18.0 | baseline | hole-phin | pass |
| Boundary/gauge | 3.55650e-10 | -19.5 | joint QFP | electron-phin | pass |
| SRH/Auger | 0.974853 | -19.5 | joint QFP | electron-phin | fail |

The SRH/Auger relative failure occurs at analytic and finite-difference matrix
norms of order `1e-15`.  Across seven perturbation steps, the largest absolute
difference is `2.70635e-15`, below the independently recorded `1e-13` absolute
floor classification threshold.  At `3e-6 V`, the best subblock comparison has
absolute difference `1.00713e-21` and relative difference `3.53218e-7`.
Smaller perturbations amplify subtraction cancellation.  The original formal
relative failure is retained rather than retrospectively reclassified as a
pass.

## First coupled update at -20 V

| Carrier | Initial residual | Trial / initial | Projection toward Sentaurus | Trial target-distance ratio |
|---|---:|---:|---:|---:|
| Electron | 8.36690e-8 | 4.18441 | -0.923986 | 1.14985 |
| Hole | 9.43464e-8 | 3.68000 | -0.889187 | 1.11458 |

Both carrier updates point away from the Sentaurus QFP target, increase their
carrier residual, and increase their distance from the target.  Together with
the residual decomposition, this establishes that the first meaningful
deviation is in carrier transport/QFP response inside the coupled solve rather
than in frozen avalanche-source evaluation.

## Determinism and limitations

All 148 node, term, update, Jacobian, and step-sensitivity artifacts are
byte-identical across two independent runs.  This experiment is a local
frozen-state plus first-update audit; it does not prove a unique source line,
establish causality for the final BV curve, or validate a correction.

The interactive HTML report passed canonical artifact validation and structural
verification.  Browser-level layout and source-dialog verification were not
run because no compatible local Chromium executable was available.

## Next read-only experiment

Keep SG/Laux unchanged.  At the carrier-residual hotspot support, decompose the
transport Jacobian into mobility, Bernoulli/GSS coefficient, QFP driving-force,
row-scaling, and contact-elimination contributions.  Finite-difference each
edge-level contribution for electron and hole equations on baseline and joint
QFP states.  Only after the first sign/scale mismatch is reproduced in that
smaller operator should an opt-in correction be designed.
