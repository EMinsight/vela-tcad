# PN2D Minimal6 Sentaurus native-edge Phase A

Date: 2026-07-24

Status: `schema-valid insufficient_data`

## Objective

Determine whether Sentaurus O-2018.06-SP2 can expose native directed electron
and hole currents on the nine Minimal6 finite-volume edges. This phase does
not modify Vela formulas and does not regenerate the 40-state lattice.

## Ordinary Plot result

An explicit SDevice `Plot` probe requested:

```text
eDensity/Edge
hDensity/Edge
eCurrentDensity/Edge/Vector
hCurrentDensity/Edge/Vector
```

SDevice rejected the command file before solving:

```text
offending input: Edge
syntax error !
```

Therefore `/Edge` is not an available SDevice `Plot` location in this
release.

## Runtime ReadFlux result

The documented CurrentPlot Tcl interface was loaded successfully. The probe
skipped the callback following the Poisson-only solve and called the runtime
API after:

```text
Coupled { Poisson Electron Hole }
```

The minimal call used an otherwise valid vertex field:

```tcl
set flux [$data ReadFlux $::des_data_edge "eDensity"]
```

SDevice returned:

```text
Unrecognized location for ReadFlux !
```

This agrees with the Device Data documentation: the actual `ReadFlux`
implementation works for vertex-based datasets only. It returns a
box-boundary surface integral of a chosen variable gradient divided by box
volume. It is not an accessor for individual directed carrier edge currents.

The carrier-field-name search is therefore terminated. Trying additional
names cannot overcome the rejected edge location.

## Runtime topology and box coefficients

The CurrentPlot Tcl mesh interface independently returned:

```text
dimension          2
runtime vertices  10
element vertices  12
edges               9
elements            4
regions             3
```

Runtime vertices 0-5 are the six physical semiconductor vertices. Vertices
6-9 duplicate contact coordinates for contact-region support. All nine bulk
edges directly use physical vertices 0-5, so no 13-to-9 edge collapse is
required.

The documented `ReadCoefficient` interface returned the following Sentaurus
box-method element-edge coefficients. The global coefficient is the sum of
the adjacent element contributions.

| Edge | Start -> end, um | Adjacent elements | Element contributions | Sum |
|---:|---|---|---|---:|
| 0 | (0,0.5) -> (0,0) | 0 | 1 | 1 |
| 1 | (0,0.5) -> (1,0) | 0,1 | 0, 0 | 0 |
| 2 | (0,0) -> (1,0) | 0 | 0.25 | 0.25 |
| 3 | (0,0.5) -> (1,0.5) | 1 | 0.25 | 0.25 |
| 4 | (1,0.5) -> (1,0) | 1,2 | 1, 1 | 2 |
| 5 | (1,0.5) -> (2,0.5) | 2 | 0.25 | 0.25 |
| 6 | (2,0.5) -> (1,0) | 2,3 | 0, 0 | 0 |
| 7 | (2,0) -> (1,0) | 3 | 0.25 | 0.25 |
| 8 | (2,0.5) -> (2,0) | 3 | 1 | 1 |

Edges 1 and 6 are the two right-triangle hypotenuses. Their contributions
are exactly zero in both adjacent elements. This directly confirms the
circumcentric finite-volume geometry: a constitutive edge current density
may be evaluated on those edges, but its integrated box flux has zero
geometric weight.

## Gate result

| Phase A condition | Result |
|---|---|
| Native edge fields through `Plot` | Failed: parser rejects `/Edge` |
| Runtime carrier flux on edge support | Failed: `ReadFlux` rejects edge location |
| Exact nine-edge topology | Passed |
| Sentaurus box coefficients | Passed |
| Native edge-current units | Not identifiable because no current is observable |
| Edge-current direction/sign | Not identifiable |
| Contact-current closure from native edges | Not executable |

Overall Phase A status is `schema-valid insufficient_data`. The native
directed-edge current required by the planned cross-solver comparison is not
observable through the documented SDevice output or CurrentPlot Tcl
interfaces in O-2018.06-SP2.

## Scientific consequence

The existing element-current direction result remains valid, and the
zero-weight diagonal geometry is now directly verified. However, node
current projection, element-vector projection, or a current recomputed from
vertex state and `ReadCoefficient` must be labeled reconstructed. None is a
native Sentaurus edge-current reference.

A production Vela current-formula change remains unauthorized. The next
valid route requires either:

1. a vendor-supported/internal interface exposing carrier-continuity
   element-edge residual contributions; or
2. a documented reconstruction of the Sentaurus box current from runtime
   coefficients and state, treated explicitly as an operator replay rather
   than a native observation.

## Evidence

Deterministic local evidence root:

`build-release/pn2d-minimal6-sentaurus-native-edge-phase-a-20260724-a`

The root contains the raw solver logs, minimal reproducer decks and Tcl
scripts, relevant official-manual excerpts, and a SHA-256 manifest.
