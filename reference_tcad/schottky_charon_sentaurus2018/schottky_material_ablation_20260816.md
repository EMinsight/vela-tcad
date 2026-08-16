# Schottky material-default ablation (2026-08-16)

## Outcome

The 1 V Vela-to-Sentaurus current error falls from 8.8886% to 0.6749% when the
Sentaurus 2018 candidate silicon material bundle is made explicit in Vela.  No
new physical model or solver feature was added.

| Case | Explicit cumulative changes | Vela current at 1 V (A/um) | Relative error at 1 V | Maximum log-current error |
|---|---|---:|---:|---:|
| Baseline | Vela built-in defaults | 1.082623822e-4 | 8.8886% | 0.478652 dex |
| `mu_n` | `mun=1417 cm2/(V s)` | 1.137182258e-4 | 4.2971% | 0.476749 dex |
| `mu_n_nc_eg` | previous plus `Nc=2.856665228e19 cm-3`, `Eg(300 K)=1.124159231 eV` | 1.142346185e-4 | 3.8625% | 0.419400 dex |
| `full_materials` | previous plus explicit `ni`, `mup`, `Nv`, affinity, permittivity and SRH lifetimes | 1.180222146e-4 | 0.6749% | 0.228313 dex |

All metrics compare the 24 nonzero Sentaurus points from 0.01 V through 1 V.
The comparison uses signed log-current interpolation, with linear interpolation
when adjacent candidate currents cross zero, matching `compare_reference_curves.py`.

## Interpretation

- Aligning electron mobility alone increases the 1 V current by 5.0395% and
  removes slightly more than half of the original error.
- Adding the 300 K conduction-band DOS and bandgap increases the 1 V current by
  a further 0.4541%.
- The remaining explicit material bundle increases it by another 3.3156%.  That
  bundle was intentionally frozen together in this pass; its `ni`, `Nv`, and
  `mup` contributions have not been individually assigned.
- The fully frozen curve is monotonic through the pseudo-arclength stage, with
  zero voltage backsteps, zero current decreases, and terminal-current
  consistency ratio 1 near 1 V.

The evidence therefore attributes the earlier 8.89% result primarily to
implicit material database differences, not to pseudo-arclength continuation.
The residual 0.675% at 1 V does not justify implementing a new Schottky boundary
or current-integration feature.  Such work should remain gated on a future
fixed-state or single-parameter ablation demonstrating a specific deficiency.

## Remaining boundary/current audit

The largest residual after full material alignment is at 0.181640361883 V:
Sentaurus gives `9.378776027e-10 A/um`, while interpolated Vela gives
`5.544124709e-10 A/um`.  This is a factor of 1.69166 (0.228313 dex).  The Vela
current there is electron dominated, and its stage-B terminal-current
consistency ratio is exactly 1 near 1 V.  Consequently:

- inconsistent terminal-current postprocessing is not supported as the source
  of the remaining error;
- the residual has the character of a low-current Schottky-emission prefactor
  or contact-state reference difference rather than missing high-field physics;
- the current evidence is insufficient to select a boundary-formula code
  change, so no such feature was implemented.

The next boundary-specific gate, if tighter low-current parity is required, is
to compare Sentaurus and Vela contact potential, electron density, and electron
normal flux at one low-bias state, or replay the same state through both contact
operators.

## Reproduction

Each case uses an independent stage-A voltage sweep and stage-B pseudo-arclength
restart under `vela/outputs/ablations/<case>/`.  Generated outputs remain
untracked.  The tracked inputs are:

- `vela/materials_mu_n_sentaurus2018.json`
- `vela/materials_mu_n_nc_eg_sentaurus2018.json`
- `vela/materials_full_sentaurus2018.json`
- the corresponding `simulation_iv_*.json` and
  `simulation_iv_*_arclength.json` decks

The stage-B decks set `stall_residual_floor=2e-9` because reloading the accepted
0.82 V state produced a residual of about 1.31e-9, just above the default 1e-9
stall threshold.  This is a restart acceptance adjustment, not a physical-model
change, and remains far below the configured nonlinear relative tolerance.
