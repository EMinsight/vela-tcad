# TransportModels Id-Vg deep-off discrepancy analysis

Date: 2026-08-21

Status: **diagnosis complete; no solver change made**.

## Technical summary

The Sentaurus current plateau near `1.6e-15 A/um` is not adequately described
as a numerical current floor. At `Vg=-1 V`, `99.9994%` of the Sentaurus drain
current is electron current, while the substrate carries an opposite hole
current with a magnitude ratio of `0.999634`.
Because the deck enables `SRH(DopingDep TempDependence)` and no Auger, BTBT, or
avalanche mechanism, this is a resolved SRH pair-generation current.

Vela does not reproduce that pair-generation plateau in either DD or DG. The
DG drain current at `-1 V` is `4.471109e-21 A/um`, while the
four-terminal KCL residual is `5.512245e-18 A/um`,
or `1232.9` times larger. The reported Vela drain
current is therefore below the numerical conservation resolution of this run.

## The discrepancy is confined to the deep-off source-dominated regime

| Vg (V) | Sentaurus DG Id (A/um) | Vela DG Id (A/um) | Relative error | Log error (dex) |
|---:|---:|---:|---:|---:|
| -1.00 | 1.636970e-15 | 4.471109e-21 | 99.9997% | 5.563625 |
| -0.84 | 1.637595e-15 | 4.471022e-21 | 99.9997% | 5.563800 |
| -0.68 | 2.162827e-15 | 4.999185e-16 | 76.8859% | 0.636123 |
| -0.52 | 4.282288e-14 | 3.868459e-14 | 9.6637% | 0.044138 |

At `Vg=-0.52 V`, ordinary channel transport has risen above the generation
plateau and the DG relative error falls to `9.66%`.

## Terminal-current decomposition identifies SRH generation

At `Vg=-1 V`:

- Sentaurus drain electron current: `1.636961e-15 A/um`.
- Sentaurus drain hole current: `9.458480e-21 A/um`.
- Sentaurus substrate hole current: `-1.636362e-15 A/um`.
- Sentaurus terminal KCL residual / drain current: `1.238e-15`.
- Vela drain electron current is exactly zero at the saved precision; its
  terminal current is set by a `-4.471109e-21 A/um` hole
  contribution.

The same approximately `1.6e-15 A/um` Sentaurus plateau is present in both its
DD and DG references, whereas both Vela DD and DG fall many orders lower. This
rules out the density-gradient equation as the primary deep-off cause.

## Verified controls

- The Sentaurus and Vela decks both enable SRH with Scharfetter doping
  dependence and the same stated `taumin`, `taumax`, `Nref`, `gamma`, and
  temperature exponents.
- Sentaurus does not enable Auger, band-to-band tunneling, avalanche, or a gate
  leakage model in this case.
- The Sentaurus and Vela comparison uses the same imported 3315-node topology.
- The discrepancy is present in DD as well as DG.

## Most likely Vela implementation gaps

1. **Cancellation in the Fermi-Dirac SRH numerator.** The Boltzmann path uses
   `ni^2 * expm1(deltaPhi/Vt)`, but the Fermi-Dirac path subtracts `n*p` and an
   independently reconstructed equilibrium product. In deep depletion these
   close quantities can lose the small net-generation signal.
2. **Incomplete generalized SRH formula.** Sentaurus generalizes SRH for Fermi
   statistics and quantization with carrier degeneracy factors in both the
   numerator and denominator. Vela reconstructs an equilibrium product for the
   numerator, but its SRH denominator remains `taup*(n+ni)+taun*(p+ni)`.
3. **Net-current resolution.** Vela's drift and diffusion diagnostics are each
   about `0.1 A/um` and cancel to the `1e-15 A/um` scale. At the two lowest gate
   biases the primary electron terminal current becomes exactly zero at saved
   precision, and the KCL residual exceeds the reported drain current.

## Secondary comparison controls

The final Vela regression sweeps the gate from `2.2 V` down to `-1 V`, while
Sentaurus initializes at `-1 V`, ramps the drain to `1.1 V`, and then sweeps the
gate upward. This is not the leading explanation because DD and DG show the
same missing plateau, but sweep direction, initial state, the `10 mV` DG outer
tolerance, and carrier residual tolerances must be matched before a final
off-state acceptance test.

## Recommended next experiment order

1. Add an SRH source audit that integrates generation over the silicon and
   reconciles it against electron and hole terminal currents at `Vg=-1 V`.
2. Implement a cancellation-free generalized Fermi SRH excess product and the
   Sentaurus degeneracy factors, with unit tests around equilibrium.
3. Require `abs(sum(I_contact))` to be at least one decade below the Id being
   compared; otherwise label the point unresolved rather than assigning a
   relative-error pass/fail.
4. Repeat forward and reverse DD sweeps with matched initialization and tighter
   tolerances; only then rerun DG.
5. Keep deep-off acceptance separate from transition/on-state acceptance and
   use both log-current error and a terminal-conservation criterion.

## Confidence and limitation

The classification of the Sentaurus plateau as SRH generation is high
confidence because it is supported by carrier-resolved terminal currents,
current conservation, and the enabled-physics deck. The exact share attributable
to Vela numerator cancellation versus the missing generalized denominator is
not yet measured; that split requires the controlled implementation A/B tests
listed above.
