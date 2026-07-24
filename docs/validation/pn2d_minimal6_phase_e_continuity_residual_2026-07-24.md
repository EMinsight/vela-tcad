# PN2D Minimal6 Phase E continuity residual and branch localization

Date: 2026-07-24

Status: complete.

## Post-fix decision

The `1e-8` two-dimensional source-to-flux factor is now implemented once in
`PhysicalUnitSystem` and applied to SRH, impact generation, carrier-term
diagnostics, analytic/finite-difference source Jacobian paths, carrier
diagonal-floor scaling, and the Gummel continuity linearization. The legacy SI
path remains exactly `1.0`.

Two fresh `mirror/sketch x -1..-20 V` roots independently pass verification
and have identical hashes for all eight sealed outputs. The post-fix
classification is:

| Gate | Post-fix result |
|---|---:|
| maximum imported-state carrier residual | `2.107970162e-8` |
| maximum first carrier-only update at `-1 V` | `0.025149361 V` |
| carrier trial/baseline residual ratio at `-1 V` | `0.356490216` |
| maximum Jacobian block relative difference | `2.989497315e-9` |
| full-physics convergence at `-20 V` | `true` |
| boundary QFP maximum error | `1.776356839e-15 V` |
| boundary density maximum error | `4.426180435e-6 dex` |
| reconstructed legacy-to-production residual median ratio | `1.500562591e-4` |
| Sentaurus box-edge mobility median residual change | `0.124446891` |
| constant-mobility median residual change | `0.790334451` |

The original source-unit diagnosis is therefore confirmed and repaired.
Residual scaling is no longer the active cause. After the repair, the dominant
fixed-state discriminator moves to SG transport/mobility support. The
Sentaurus branch is still a box-operator reconstruction, not a native directed
edge-current observation, so this result does not yet justify a production
mobility or SG formula change.

The pre-fix localization evidence below is retained as the RED baseline.

## Pre-fix decision

Phase E identifies a reproducible Vela implementation-defect candidate in the
unit conversion between 2-D volumetric continuity sources and integrated SG
edge flux. It does not support a production mobility-formula change.

The fixed-state evidence gives the following classification:

| Candidate cause | Result |
|---|---|
| mobility coefficient | secondary; median residual change `1.149956e-5` for the Sentaurus box-edge branch |
| contact/boundary state | rejected; QFP maximum error `1.776357e-15 V` |
| boundary density conversion | rejected; maximum error `4.426180e-6 dex` |
| residual scaling | primary defect candidate; missing 2-D source-to-flux factor `1e-8` |
| analytic Jacobian implementation | rejected as the primary cause; maximum block relative difference `5.048694e-8` |
| first fixed-psi update | unstable because the unconverted source dominates; not an analytic/FD mismatch |
| high-bias nonlinear branch | separate secondary failure at `-20 V`, coupled to avalanche |

No production residual, Jacobian, mobility, SG, recombination, or impact
formula was changed in Phase E.

## Evidence contract

The audit covers the exact `mirror/sketch x -1..-20 V` lattice.

| Evidence | Count |
|---|---:|
| exact imported states | 40 |
| residual-waterfall rows | 480 |
| contact/boundary rows | 160 |
| Jacobian block rows | 200 |
| internal-node first-update rows | 80 |
| source-unit control rows | 160 |
| controlled branch experiments | 54 |

Each imported state uses Sentaurus `psi`, electron QFP, and hole QFP. Vela
recomputes carrier density with its Old Slotboom BGN implementation before the
continuity operator is evaluated.

The three mobility branches are:

1. native Vela production Masetti plus QFP-gradient high-field mobility;
2. coefficient-weighted Sentaurus box-edge mobility reconstruction; and
3. native Vela constant-mobility control.

The production and constant branches are recomputed by the C++ assembler. The
Sentaurus branch uses exact SG linearity with edge mobility. Its impact source
is reweighted on the alpha-current edge support; replay against the native
constant branch bounds that source reconstruction to `0.00359821` maximum
relative error.

## Residual sign and orientation closure

The directed edge convention is the production convention:

- add edge flux to `node0`;
- subtract edge flux from `node1`;
- add SRH to electron and hole rows;
- subtract impact generation from electron and hole rows; and
- replace contact carrier rows by QFP Dirichlet residuals.

After increasing only the diagnostic CSV precision to 17 digits, the
production edge-to-node replay closes to `2.194488e-16` maximum relative
error. The constant branch closes to `4.517416e-16`.

## First exact departure

Both topologies depart at the first exact state, `-1 V`.

| Topology | Bias | Maximum imported-state carrier residual | Maximum carrier-only first update | Trial/baseline carrier residual ratio |
|---|---:|---:|---:|---:|
| mirror | `-1 V` | `1.441532e-4` | `213.047104 V` | `1.100957e217` |
| sketch | `-1 V` | `1.441532e-4` | `213.047104 V` | `1.100957e217` |

The mirror and sketch values agree after the previously verified cell/support
mapping. There is no topology-specific first divergence.

At mirror `-1 V`, the representative electron row at node 5 is:

| Term | Normalized value |
|---|---:|
| SG divergence | `1.911105266e-8` |
| SRH as currently assembled | `-1.441680000e-4` |
| impact | approximately `-2.61e-62` |
| original residual | `-1.441488889e-4` |
| source-unit corrected residual | `1.910961098e-8` |

Thus the current SRH contribution is about four orders of magnitude larger
than the SG divergence at the first state. Mobility substitution changes this
residual only at the `1e-5` median relative level.

## Source-unit derivation

The TCAD internal unit system uses:

- concentration: `1 internal = 1e6 m^-3`;
- coordinate: `1 internal = 1e-6 m`; and
- mobility: `1 internal = 1e-4 m^2/(V s)`.

The 2-D volumetric source integral and integrated SG flux therefore differ by

`f_source = (1e6 * (1e-6)^2) / (1e6 * 1e-4) = 1e-8`.

In the current coupled residual, Poisson charge multiplies its explicit unit
conversion factor, but `R * node_area` and volumetric impact integrals are
added directly beside SG edge flux before all carrier rows are divided by
`C0*D0`. The analytic Jacobian mirrors the same residual convention, which is
why analytic/finite-difference agreement can pass while the physical residual
still has the wrong source/flux unit ratio.

Applying `1e-8` to SRH and impact only as an offline fixed-state control reduces
the absolute residual by a median factor of `1.5005626e-4`. This is a causal
RED result for a future production patch, not the patch itself.

## Boundary and Jacobian audit

| Gate | Result |
|---|---:|
| contact QFP maximum absolute error | `1.776357e-15 V` |
| contact BGN density maximum absolute error | `4.426180e-6 dex` |
| maximum Jacobian block relative difference | `5.048694e-8` |
| Jacobian classification | `jacobian_consistent` |
| boundary classification | `boundary_consistent` |

The Jacobian audit covers Poisson, transport, SRH/Auger, SG avalanche, and
Dirichlet/gauge blocks at all 40 states. It rejects a derivative coding error
as the primary low-bias cause. It does not validate the physical unit contract:
the analytic and finite-difference Jacobians differentiate the same
unconverted residual.

## Controlled branch experiments

Experiments were run at `-1`, `-10`, and `-20 V` on both topologies.

| Experiment | `-1 V` | `-10 V` | `-20 V` |
|---|---|---|---|
| fixed `psi`, carrier-only Newton iteration | fails after first huge step | fails after first huge step | fails after first huge step |
| coupled, avalanche disabled | converges | converges | converges |
| coupled, SRH disabled | converges | converges | fails before an accepted step |
| full physics | converges | converges | fails after one iteration |
| QFP homotopy, all five fractions | all converge | all converge | all fail |

At `-1 V`, disabling SRH reduces the coupled initial residual from
`2.279164665e-4` to `3.576378907e-8`. This independently confirms that the
low-bias imported-state residual is source dominated.

At `-20 V`, disabling avalanche changes the full-physics failure into a
converged solve. This is a separate high-bias avalanche/nonlinear-branch
problem. It does not explain the first `-1 V` QFP departure, where the impact
term is effectively zero.

## Phase E exit

Phase E passes its localization gate with:

- primary result: `volumetric_source_unit_conversion` implementation-defect
  candidate;
- secondary result: high-bias avalanche/nonlinear branch failure;
- rejected primary causes: mobility, contact state, boundary density
  conversion, and analytic Jacobian inconsistency; and
- remaining limitation: Sentaurus native directed-edge current is still
  unavailable, so the Sentaurus box-edge mobility branch remains an operator
  reconstruction.

The approximately `0.33 V` self-consistent QFP difference is therefore not
explained by the remaining mobility mismatch. The next allowed action is a
focused RED/GREEN production patch that applies the 2-D source/flux conversion
consistently in the carrier residual, term diagnostics, and all corresponding
Jacobian derivatives. Phase B/C gates and the full 40-state sweep must then be
rerun before accepting the patch.

## Deterministic evidence

Two independent roots were generated and verified:

Post-fix roots:

- `build-release/pn2d-minimal6-phase-e-continuity-sourcefix2-20260724-a`
- `build-release/pn2d-minimal6-phase-e-continuity-sourcefix2-20260724-b`

Both post-fix independent verifiers pass. All eight sealed output hashes match,
and the manifests are byte-identical with SHA-256
`dfb69a79e568d184e5686012fa3e6dba83f494de4c58e58a6bcb0f35a7429ab4`.

Pre-fix roots retained for the RED baseline:

- `build-release/pn2d-minimal6-phase-e-continuity-20260724-a`
- `build-release/pn2d-minimal6-phase-e-continuity-20260724-b`

Both independent verifiers pass. All eight sealed output hashes match, and the
two manifest files are byte-identical with SHA-256
`42515ed3da7fca2c351943adbd4c413eb6366d003e98a6712a8df168b9076d0e`.
