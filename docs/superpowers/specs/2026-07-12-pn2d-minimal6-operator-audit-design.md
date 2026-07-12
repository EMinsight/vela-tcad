# PN2D Minimal6 Fixed-State Operator Audit Design

## Purpose

Create a deliberately minimal PN2D mesh with six nodes and four triangles so
that Sentaurus and Vela local transport and avalanche operators can be audited
edge by edge and triangle by triangle. The fixture is diagnostic-only. Its
current-voltage behavior and apparent breakdown voltage are not physical BV
validation results and must not replace the existing PN2D reference meshes.

The audit uses Sentaurus-converged nodal states as immutable input to Vela.
Vela does not run Newton, Gummel, continuation, or an independent DC solve in
this workflow. This isolates discretization and current/source semantics from
state-branch and solver differences.

## Geometry And Canonical Nodes

Keep the current coarse PN2D dimensions and junction position:

```text
length = 2.0 um
height = 0.5 um
junction x = 1.0 um
```

The canonical node labels and coordinates are:

| node | x (um) | y (um) | doping |
|---:|---:|---:|---|
| 1 | 0.0 | 0.5 | Acceptor = 1e17 cm^-3 |
| 2 | 1.0 | 0.5 | Acceptor = Donor = 1e17 cm^-3 |
| 3 | 2.0 | 0.5 | Donor = 1e17 cm^-3 |
| 4 | 2.0 | 0.0 | Donor = 1e17 cm^-3 |
| 5 | 0.0 | 0.0 | Acceptor = 1e17 cm^-3 |
| 6 | 1.0 | 0.0 | Acceptor = Donor = 1e17 cm^-3 |

Nodes 2 and 6 intentionally retain the native compensated-junction semantics.
Their net signed doping is zero. The Anode is boundary edge 1-5 and the Cathode
is boundary edge 3-4. All four triangles belong to one silicon region.

Sentaurus internal node numbers are not authoritative. Every export is mapped
to the canonical labels by exact coordinate matching within `1e-12 um`.
Nearest-neighbor matching and interpolation are prohibited.

## Topology Variants

Both variants use counter-clockwise triangle orientation.

The sketch topology is:

```text
(1,5,2), (5,6,2), (2,6,4), (2,4,3)
```

The vertically mirrored topology is:

```text
(1,5,6), (1,6,2), (2,6,3), (6,4,3)
```

Each valid topology must contain exactly six nodes, four triangles, nine unique
edges, six boundary edges, three interior edges, and two contact edges. The
sketch and mirror variants are separate fixtures and artifact roots. Their
results are never merged before orientation-sensitivity ratios are computed.

## Authoritative Mesh Path

An explicit topology definition is the single source of truth. It defines
coordinates, triangle connectivity, silicon ownership, contact edges, nodal
doping, units, and a stable topology ID. The Sentaurus toolchain converts this
definition to a grid accepted by SDevice without changing the topology. Vela
then imports the exact Sentaurus TDR rather than constructing a parallel mesh.

The first implementation gate proves this conversion path in the Sentaurus VM.
The gate fails if Sentaurus rejects the explicit grid, adds a node, removes a
node, changes connectivity, flips a requested diagonal, changes a contact edge,
or changes a nodal doping value. No SDE/SNMesh automatic-remeshing fallback is
allowed. No separately generated Vela mesh is allowed.

## Bias States And Data Flow

For each topology, Sentaurus solves exactly these states:

```text
0 V, -12 V, -19 V
```

The requested bias must have an exact converged export. A nearest available
bias is not a valid substitute. Each run records commands, logs, convergence
status, source and output hashes, and the Sentaurus version.

The data flow is:

```text
explicit topology
  -> Sentaurus grid/TDR
  -> topology and doping contract gate
  -> Sentaurus fixed-bias solution
  -> complete state-field export
  -> Vela imports the same TDR and immutable state
  -> production C++ local operator evaluation
  -> independent Python mathematical reference
  -> joined audit tables and static figures
```

Vela must not alter `psi`, `phin`, `phip`, `n`, or `p` before evaluating local
operators. The imported-state parity check runs before any physics comparison.

## Required Sentaurus Fields

Each canonical node requires:

```text
ElectrostaticPotential
eQuasiFermiPotential
hQuasiFermiPotential
eDensity
hDensity
DonorConcentration
AcceptorConcentration
eMobility
hMobility
ElectricField vector
eCurrentDensity vector
hCurrentDensity vector
eAlphaAvalanche
hAlphaAvalanche
```

The field manifest must record region, component count, unit, mapping status,
and global-node mapping. Missing, scalar-only, incomplete, ambiguously mapped,
or incorrectly dimensioned fields fail the affected bias state.

## Audit Layers

### Node State

Write one wide row per topology, bias, and canonical node. The table preserves
the original Sentaurus values, converted SI values, values received by Vela,
absolute error, and relative error. Expected row count is:

```text
2 topologies * 3 biases * 6 nodes = 36 rows
```

### Edge Operators

Write one wide row per topology, bias, and canonical edge. Electron and hole
quantities remain in the same row so their endpoint and orientation conventions
cannot drift. Record:

- endpoint state, edge vector, length, canonical sign, and edge class;
- `Delta(phiFn)/h` and `Delta(phiFp)/h`;
- Bernoulli arguments, factors, and GSS midpoint densities;
- carrier mobilities;
- Vela production SG raw flux;
- PDF grad-qF flux;
- Vela reconstructed-vector projection and magnitude;
- Sentaurus current-vector projection and magnitude;
- carrier alpha, edge area proxy, and edge source contributions.

Expected row count is:

```text
2 topologies * 3 biases * 9 edges = 54 rows
```

### Triangle Operators

Write one wide row per topology, bias, and canonical triangle. Record:

- triangle orientation, signed area, and shape-function gradients;
- `grad(psi)`, `grad(phiFn)`, and `grad(phiFp)`;
- reconstructed electron and hole current vectors;
- electron and hole alpha;
- each local edge's Genius-truncated partial volume;
- `alpha * |J| / q` source terms;
- electron and hole nodal source partitions.

Expected row count is:

```text
2 topologies * 3 biases * 4 triangles = 24 rows
```

## Implementations Under Test

The Vela values must come from production C++ operator code. A dedicated
fixed-state diagnostic entry point may expose intermediate factors, but it must
call the same helpers used by residual assembly and diagnostics. Python is an
orchestrator and report generator, not a replacement implementation.

An independent Python reference computes triangle gradients, Bernoulli values,
SG fluxes, midpoint densities, and geometry directly from the documented
equations. It is an oracle for the C++ implementation and cannot be used to
populate columns labeled as Vela production results.

## Gates And Classification

Automatic gates are:

- coordinate error less than `1e-12 um`;
- exact topology counts and exact canonical connectivity;
- exact Anode and Cathode edge ownership;
- exact compensated doping semantics at nodes 2 and 6;
- imported Sentaurus-to-Vela state error less than `1e-12` under a documented
  absolute/relative hybrid norm so exact and near-zero values are well-defined;
- C++ versus Python geometry, gradient, Bernoulli, and SG formula relative
  error less than `5e-12`, using an absolute tolerance for values near zero;
- complete finite values for every required nonzero comparison;
- exact output row counts of 36, 54, and 24 with unique keys.

Zero values are explicitly classified and are not assigned artificial log
errors. Sentaurus-versus-Vela current and avalanche-source differences have no
automatic improvement threshold in the first audit. They are diagnostic
measurements because the products may use different current semantics.

The summary reports sketch/mirror ratios for each shared canonical quantity,
including integrated electron, hole, and total avalanche source. This measures
mesh-direction sensitivity separately from cross-product disagreement.

## Failure Behavior

- A topology conversion failure stops the workflow before state solving.
- A failed or inexact Sentaurus bias is reported and excluded; it is not
  replaced by another bias.
- A missing required field stops the affected topology/bias audit.
- Duplicate coordinates, duplicate triangle keys, unexpected edges, reversed
  orientation, unit mismatches, and non-finite values are fatal contract errors.
- Partial matrices are labeled incomplete and cannot be published as a complete
  six-state comparison.

The workflow never enables a QF clamp, source limiter, alpha scale, contact
fallback, or solver tuning. It does not modify Vela's default
`density_gradient` mode or the existing triangle-GSS experimental mode.

## Reports And Figures

The artifact root contains:

```text
manifest.json
node_state.csv
edge_audit.csv
triangle_audit.csv
summary.json
summary.md
figures/minimal6-topologies.{png,pdf}
figures/minimal6-edge-current-audit.{png,pdf}
figures/minimal6-triangle-source-audit.{png,pdf}
```

The manifest records topology definitions, bias states, model configuration,
tool versions, command status, file hashes, row counts, and gate status. Static
figures show canonical labels and avoid implying that the three bias samples
form a validated BV curve.

## Test Strategy

Unit tests cover canonical coordinate mapping, both topology edge sets,
counter-clockwise orientation, mirror equivalence, contact ownership, and
compensated-node doping. Formula tests cover triangle gradients, Bernoulli
limits, SG signs, midpoint density, partial-volume sums, and zero handling.

Contract tests reject missing fields, wrong vector component counts, duplicate
nodes, changed diagonals, inexact biases, wrong units, incomplete manifests,
and row-key duplication. An integration fixture exercises a synthetic
Sentaurus export through C++ replay and verifies the exact 36/54/24 output
matrix. Real Sentaurus VM execution is a separate validation gate and is not a
network-dependent CI requirement.

## Delivery Boundary

The first delivery is only the explicit-grid Sentaurus compatibility gate for
both topologies. Full C++ replay and reporting work begins only after both grids
round-trip through Sentaurus with unchanged topology. Final completion requires
all six exact topology/bias states and all automatic gates. No physical BV
conclusion is part of this design.
