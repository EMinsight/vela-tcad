# SingleDevice Eq. 231 Sentaurus-box Jacobian probe (2026-08-14)

> 2026-08-15 correction: `AverageBoxMethod` uses signed circumcentric measures
> on the obtuse SingleDevice cells.  The earlier positive mixed-Voronoi label
> is superseded by
> `singledevice_eq231_averagebox_polyoxide_central_closure_2026-08-15.md`.

## Outcome

The fixed-state mismatch is no longer attributable to an ordinary Tri3
quadratic-gradient integral.  Direct Sentaurus NewtonPlot perturbations show
that the existing `gss_potentiallike_fitted` edge operator reproduces the
Formula-0 oxide-interior Jacobian.  The remaining discrepancy is confined to
material-interface traces and the interface-adjacent control volume.

An opt-in `sentaurus_box` prototype was implemented.  It combines the verified
fitted edge flux, Sentaurus `MeasureCoefficients` mixed-Voronoi
element-vertex measures, and the maximum region-side Lambda trace in the box
reaction.  The linear fixed-state maximum residual fell from `952.158169` to
`18.350359`, but the imported endpoint still moved by `1.000394 V` and failed
the reference-state gate.  The saturation endpoint and both 21-point Id-Vg
curves were therefore not run.

## Direct residual/Jacobian oracle

The oracle loads the converged `Vg=2.2 V, Vd=0.1 V` Sentaurus state, perturbs
one `eQuantumPotential` node value inside the VM, and executes exactly one
Newton iteration of the quantum-potential equation.  The unperturbed
NewtonPlot RHS is at roundoff.  A separate `theta=0` parameter set isolates
the linear box and reaction terms.

At oxide node 1848, the current fitted operator and Sentaurus have the same
directed support: columns 1847, 1849, and 1852 are nonzero, while the strongly
downhill 1859/1861 edges are exponentially suppressed.  After correcting for
the actual TDX perturbation written to each field, corresponding derivatives
agree within about `0.1%-0.6%`.  The current fitted self derivative is
`6943.119`; the Sentaurus finite-difference value is approximately `6.95e3`.

At shared node 1847, perturbing only the oxide TDR occurrence has no effect;
perturbing the PolySilicon occurrence changes both exported copies and the
RHS.  This proves that the repeated TDR value is an owner/output trace, not an
independent oxide-side restart degree of freedom.  The off-diagonal linear and
fitted edge derivatives still agree, while the combined interface residual
and diagonal do not.  The missing term is therefore an interface trace/control
volume term rather than a bulk stiffness coefficient.

## GateOx interface-adjacent row

After the first prototype, the leading bulk-looking residuals moved to GateOx
nodes adjacent to the PolySilicon interface.  Node 2074 has fitted stiffness
`164.541792`, reaction `-177.975285`, and residual `-13.433493`.

Its two dominant edges terminate at shared node 2075.  The current oxide-side
half jump is `h=0.298765`.  Sentaurus's residual and finite-difference
derivative are simultaneously reproduced by `h approximately 0.3189`, which
corresponds to only about a `1.04 mV` region-side potential-like trace shift.
This explains why the same fitted formula is accurate at oxide-interior node
1848 but misses the interface-adjacent row: Sentaurus uses an internal
region-side interface trace that is not represented by the single exported
TDR value.

## Implementation and verification

- Added `global_discretization: sentaurus_box` as an experimental mode.
- Added mixed-Voronoi element-vertex reaction measures.
- Added a maximum material-side reaction-trace prototype with a consistent
  Jacobian.
- Added a two-material manufactured solution for the trace contract.
- `test_density_gradient_quantum_potential.exe`: 81 assertions in 21 cases
  pass.
- Linear imported endpoint: quantum inner solve converges in 25 iterations,
  but the one-step endpoint gate fails because the raw field change is
  `1.000394 V`.

## Decision

The prototype is retained only as a diagnostic.  The next core implementation
must represent a region-side interface trace (or an algebraically eliminated
equivalent) and its coupled Jacobian.  Further adjustment of bulk theta,
solver damping, ordinary P1/CVFEM volume integration, or a WKB boundary is not
supported by the new Jacobian evidence.

The original `SingleDevice` command enables plain `eQuantumPotential`; it does
not request the `Resolve` switch.  The O-2018.06 manual describes `Resolve`
only as a more accurate discretization for discontinuous band structures at
non-heterointerfaces and does not disclose either the default or resolved box
formula.  Consequently, substituting a guessed `Resolve` discretization would
not be an equivalent reproduction of this deck.

A direct TDR field cross-check also excludes a missing band-edge reconstruction
term as the approximately 1 mV trace offset.  At the PolySilicon/GateOx owner
vertices, the exported `ElectronAffinity` includes exactly half of the exported
OldSlotboom `BandgapNarrowing`; the Vela cell-side drive uses the same affinity,
BGN split, and the documented SiO2 electron DOS mass of 0.42.  The unresolved
offset is therefore a numerical interface-trace effect, not an omitted bulk
band/DOS material parameter.
