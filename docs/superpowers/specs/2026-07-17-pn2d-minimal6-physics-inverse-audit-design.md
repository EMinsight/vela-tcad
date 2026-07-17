# PN2D Minimal6 physics inverse-audit design

Date: 2026-07-17

## Purpose

Produce a provenance-complete diagnostic report that compares and reverse-infers
the electrostatic potential, electric-field vector, electron and hole
quasi-Fermi-potential gradients, current densities, avalanche driving fields,
ionization coefficients, and impact-ionization generation used by Vela and
Sentaurus on the PN2D Minimal6 device.

This phase is diagnostic only. It may add analysis scripts, tests, immutable
report artifacts, and documentation, but it must not modify production C++
formulas or solver behavior. A later, separately approved task may use an
identified formula candidate as the basis for a test-driven production change.

## Baseline and constraints

Tasks 5-8 of the Minimal6 formula-difference and sweep plan are complete. The
fresh Task 8 roots provide exact sketch and mirror checkpoints across the common
`-1..-20 V` range. The Sentaurus exports include electrostatic potential,
electron and hole densities, electron and hole quasi-Fermi potentials,
electric-field vectors, electron and hole current-density vectors, avalanche
coefficients, and native impact-ionization generation.

Existing evidence establishes that electrostatic potential and peak electric
field are already close while the self-consistent quasi-Fermi branch, carrier
density, current, and source differ materially. That evidence motivates this
audit but does not predetermine the recovered formulas. The audit must preserve
typed `insufficient_data` or confounded outcomes when the exported data cannot
identify a unique operator.

The following constraints are mandatory:

- Inputs are immutable and hash-addressed.
- Sketch and mirror data are both required.
- No bias-, node-, edge-, cell-, carrier-, or topology-specific fitted scale is
  allowed in an identified formula.
- Raw fields are never overwritten by normalized or derived values.
- Relative errors exclude typed geometric zeros and values below declared
  numerical floors.
- Production C++ files must remain unchanged throughout this phase.

## Chosen approach

Use a layered physical inverse audit with one-factor and staged data
replacement. This approach was selected over two alternatives:

1. Direct black-box fitting could reduce aggregate error quickly, but it would
   confound units, support reconstruction, mobility, driving force, and source
   mapping.
2. Whole-state Sentaurus replay is useful as a localization control, but it
   cannot by itself identify an accurate formula.

The chosen approach first aligns support and units, then examines each physical
layer independently, and finally replays replacements through declared
dependencies. It favors causal and cross-bias evidence over aggregate curve
agreement.

## Architecture

The audit is split into five isolated components.

### 1. Provenance and input validation

The input loader validates the sealed Task 8 Vela and Sentaurus roots, strict
manifests, topology, exact bias, field component count, field unit, and file
hashes before reading values. It records the executable and tracked-source
hashes already bound into the source manifests.

If an independently required field is absent, the audit reports the missing
field before attempting inference. A targeted remote Sentaurus regeneration is
permitted for diagnostic fields such as electron and hole mobility. Any such
regeneration must preserve the device, mesh, physics, topology, bias matrix,
and Sentaurus version and must retain the generated deck, logs, raw TDR files,
imported fields, remote provenance, and hashes.

### 2. Canonical support and unit layer

All observations are represented by the key

```text
(solver, topology, bias_V, support_kind, support_id, quantity, component)
```

where `support_kind` is one of `node`, `edge`, `cell`, `contact`, or
`integrated`. Each row stores the raw value and unit, normalized SI value,
coordinate frame, orientation convention, source path and hash, and conversion
provenance.

The layer owns only transformations that can be independently verified:

- metres, centimetres, and micrometres;
- node, edge, and triangle identifiers;
- sketch-to-mirror coordinate and vector transformations;
- scalar, vector, projection, magnitude, and signed-flux semantics;
- node-to-cell, cell-to-node, node-to-edge, and cell-to-edge support recovery.

It must not embed a physical model or silently choose a candidate formula.

### 3. Candidate-formula evaluators

Each evaluator accepts canonical observations and returns derived observations,
validity masks, dimensional checks, and a dependency record. Evaluators remain
independent of report formatting.

### 4. Replacement and identifiability engine

The engine evaluates baseline, single replacement, staged replacement, reverse
replacement, and full replacement over an explicit dependency graph. It emits
the downstream metric change and closure residual for every replacement. It
does not infer causality from a final total that closes only through cancellation.

### 5. Report and verification package

The report builder consumes validated tables and emits machine-readable results,
static figures, a concise scientific narrative, and a manifest. A standalone
verifier replays key calculations from raw inputs and validates report hashes.

## Candidate formulas

### Electrostatic field

Compare the direct Sentaurus electric-field vector with:

- the constant triangular gradient `E = -grad(psi)`;
- area-weighted cell gradients reconstructed at nodes;
- area-weighted adjacent-cell gradients reconstructed at edges;
- the signed edge difference `-(psi_j - psi_i) / h`;
- projections of reconstructed vectors along the edge and current directions.

Magnitude and direction are assessed separately. Near-zero components do not
use component-wise relative error.

### Quasi-Fermi-potential gradients

Evaluate electron and hole quantities independently using:

- constant triangular `grad(phi_n)` and `grad(phi_p)`;
- signed edge differences;
- area-weighted node and edge reconstructions;
- vector magnitude and projection along the local current direction;
- effective gradients recovered from current,
  `|J_n|/(q*mu_n*n)` and `|J_p|/(q*mu_p*p)`.

The current-inverted form identifies a gradient only when mobility and carrier
density are independently available on a compatible support. Otherwise it is
reported as the confounded product `mu*|grad(phi_qf)|` or as an effective
mobility-gradient combination. The report must not promote it to a unique
gradient formula.

### Current density

Normalize Sentaurus nodal current-density vectors and Vela SG edge carrier
fluxes to compatible signed current units and compare:

- native SG carrier flux;
- quasi-Fermi-gradient current `q*mu*n*grad(phi_n)` and
  `q*mu*p*grad(phi_p)` under explicit carrier sign conventions;
- drift and diffusion terms separately;
- nodal current projected to an edge;
- edge currents reconstructed to a node or cell vector.

The audit distinguishes current magnitude, signed normal flux, and vector
current. Agreement in one semantic form cannot establish agreement in another.

### Avalanche driving force and alpha

Invert the configured piecewise Van Overstraeten coefficient only on valid,
non-saturated branches with known coefficient provenance. Compare the resulting
effective field with:

- `|E|`;
- `|grad(phi_n)|` and `|grad(phi_p)|`;
- the corresponding field projected along carrier current;
- the existing low-density interpolation candidates when their reference
  density is independently known.

Direct alpha agreement is a control, not proof of the driving force, when
multiple driving fields produce numerically indistinguishable coefficients.

### Impact-ionization generation and mapping

Reconstruct

```text
G_II = (alpha_n*|J_n| + alpha_p*|J_p|) / q
```

on compatible support, then independently apply triangle integration, edge
partial volumes, out-of-plane normalization, and node mapping. Compare native
Sentaurus generation, Sentaurus alpha-current reconstruction, Vela native
source, and Vela formula reconstruction without mixing volumetric and integrated
units.

## Replacement data flow

The declared upstream-to-downstream order is:

```text
potential/state
  -> gradient recovery
  -> mobility and carrier support
  -> current semantics
  -> avalanche driving force
  -> alpha law
  -> geometric integration
  -> source-to-node mapping
```

The potential/state layer is an immutable observed input in the primary audit.
Whole-state replay is retained only as a localization control. The causal matrix
then replaces, in order, gradient recovery, mobility, current semantics,
driving force, alpha, integration geometry, and source mapping. Each factor is
also replaced in isolation and in reverse order. Every result records which
observations are native, reconstructed, replaced, or unavailable.

## Discovery and validation split

The discovery set is the seven sketch checkpoints at `-1`, `-4`, `-8`, `-12`,
`-16`, `-19`, and `-20 V`, spanning low field, transition, and high field.
Formula ranking and any global constant selection use only the discovery set.

The holdout set contains the mirror topology and all remaining exact common
checkpoints. A formula fails identification if it needs holdout-specific
adjustment. Final tables report all exact common `-1..-20 V` checkpoints for
both topologies, even when individual samples are masked by a typed validity
rule.

## Identifiability and error handling

Each candidate receives exactly one primary classification:

- `identified`: one candidate satisfies all gates without local fitting;
- `consistent_nonunique`: multiple candidates cannot be distinguished on the
  available Minimal6 support;
- `confounded`: only a product or combined operator is identifiable;
- `insufficient_data`: required independent evidence or valid support is absent;
- `rejected`: units, direction, symmetry, cross-bias stability, or error gates
  fail.

Samples are separately typed as valid, geometric zero, below numerical floor,
missing field, incompatible support, invalid unit, direction undefined,
coefficient branch ambiguous, exponential underflow, or non-finite. Masking a
sample requires a recorded reason. Missing or invalid values never become zero.

## Acceptance gates

All thresholds apply to the holdout set as well as the combined final matrix.

### Electrostatic field

- Median magnitude error on valid high-field samples is at most 2%.
- Median vector-direction error is at most 1 degree.
- Sketch/mirror transformed results satisfy the same gates.

### Quasi-Fermi gradient and current

- Median absolute log10 magnitude error is at most 0.1 dex.
- The 95th percentile absolute log10 magnitude error is at most 0.3 dex.
- Median direction error on direction-defined samples is at most 5 degrees.
- No per-bias, per-support, per-carrier, or per-topology fitted scale is used.

### Alpha and impact generation

- Valid high-field integrated alpha/current generation differs by at most
  0.1 dex.
- Valid local node, edge, or cell support differs by at most 0.3 dex.
- A driving-force conclusion must remain stable across the configured
  low/high coefficient branch transition.

### Causal replacement

- A named replacement must reduce the error of its declared downstream
  quantity rather than only the final aggregate.
- Reverse replacement must restore the corresponding baseline discrepancy
  within the numerical replay tolerance.
- Full replacement closure and the sum of declared residual terms must close
  within `1e-10 dex` for finite nonzero states.

### Package integrity

- Unit and dimensional checks pass.
- Sketch/mirror invariance passes.
- Two identical runs produce byte-identical machine-readable outputs and
  rendered figures.
- A standalone verifier passes against the sealed inputs.
- The report manifest covers every input and output artifact.
- A final tracked diff against the phase-start baseline `a5524cf` confirms that
  production C++ files are unchanged.

Failure of an identification threshold does not fail the package when the
candidate is accurately classified as rejected, confounded, nonunique, or
insufficient. Schema, provenance, determinism, and classification consistency
remain mandatory.

## Report artifacts

The diagnostic root will contain:

- `input_manifest.json` with source roots, hashes, versions, and field inventory;
- normalized node, edge, cell, contact, and integrated CSV tables;
- `candidate_metrics.csv` and `candidate_classifications.json`;
- a one-factor and staged `replacement_matrix.csv`;
- per-bias and aggregate comparisons for potential, field, quasi-Fermi
  gradients, current, alpha, and generation;
- static PNG/PDF figures with a deterministic figure manifest;
- `physics_inverse_audit.json` as the authoritative machine-readable result;
- `physics_inverse_audit.md` as the scientific report;
- `verification.json` and a standalone verification entry point.

The report must state the recovered equation, sign convention, units, support,
validity domain, fitted constants if any, discovery/holdout results, and
identifiability classification for every promoted candidate. It ends with a
ranked list of candidates that could justify a later production-code task, but
does not implement them.

## Test strategy

Tests cover:

- exact unit conversions and dimensional rejection;
- analytic triangle gradients and edge projections;
- vector mirror transformations and orientation signs;
- node, edge, and cell support reconstructions;
- current and alpha inverse identities on synthetic non-degenerate fixtures;
- typed zero, floor, missing-data, ambiguity, and non-finite handling;
- single, staged, reverse, and full replacement closure;
- discovery/holdout separation and prohibition of local fitted scales;
- deterministic JSON, CSV, figure, and manifest generation;
- semantic replay by an independent verifier;
- an explicit guard that production C++ paths do not change in this phase.

## Completion criterion

This phase is complete when the report package passes all integrity and
reproducibility checks and assigns every candidate a defensible identifiability
classification. It is not necessary to force an `identified` result. An
evidence-backed `consistent_nonunique`, `confounded`, or `insufficient_data`
outcome is scientifically complete when the report states exactly which
additional independent field or mesh resolution would distinguish the
candidates.
