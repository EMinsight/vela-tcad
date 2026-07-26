# PN2D high-bias Sentaurus process-variable contract

Date: 2026-07-26

Plan task: Task 13 of
`docs/superpowers/plans/2026-07-26-pn2d-high-bias-process-variable-jacobian-localization.md`

## Outcome

Typed outcome: `native_process_variable_contract_available`

The coarse7x3 one-state contract is accepted at exact anode bias
`-19.95 V` on Sentaurus Device `O-2018.06-SP2`. This evidence may enter
Task 14. It does not override the Task 12 nonlinear stop and cannot authorize
a Vela self-consistent candidate.

## Raw roots and determinism

- `build-release/pn2d-task13-process-probe-20260726-a`
- `build-release/pn2d-task13-process-probe-20260726-b`

Both roots observed exactly one runtime callback at `-19.95 V`. Their deck,
mesh, parameter, and Tcl bundle SHA-256 values are identical. Their CurrentPlot
`.plt` files are byte-identical.

The raw TDR and solver logs contain run-path metadata and therefore have
different raw hashes. After local TDR parsing:

- 46 datasets were inventoried in each root;
- 53 exported files were compared;
- the only differing file was `metadata.json`, which contains the source path;
- all 52 physical, contact, field, and geometry exports were byte-identical;
- all 260 normalized `AVAL_PROBE_*` runtime records were byte-identical.

Runtime record counts were:

| Record | Count |
|---|---:|
| BEGIN / END | 1 / 1 |
| physical and contact-support vertices | 33 |
| semiconductor elements | 32 |
| element-local edges with `ReadCoefficient` | 96 |
| element-local vertices with `ReadMeasure` | 96 |
| carrier-split source integral | 1 |

The ordinary TDR contains 27 physical semiconductor vertices and 32 Tri3
elements. The six additional runtime vertices are contact-support duplicates;
they are not inserted into the physical-node field arrays.

## Accepted native fields

| Quantity | Centering | Components | Unit | Label |
|---|---|---:|---|---|
| potential, electron/hole QFP | vertex | 1 | V | `native_node` |
| electron/hole density | vertex | 1 | cm^-3 | `native_node` |
| electric field | vertex and element | 2 | V cm^-1 | `native_node`, `native_element` |
| electron/hole QFP gradient | vertex and element | 2 | V cm^-1 | `native_node`, `native_element` |
| electron/hole mobility | vertex and element | 1 | cm^2 V^-1 s^-1 | `native_node`, `native_element` |
| electron/hole current density | vertex and element | 2 | A cm^-2 | `native_node`, `native_element` |
| total current density | vertex | 2 | A cm^-2 | `native_node` |
| electron/hole velocity magnitude | vertex | 1 | cm s^-1 | `native_node` |
| electron/hole avalanche alpha | vertex | 1 | cm^-1 | `native_node` |
| electron/hole/total impact generation | vertex | 1 | cm^-3 s^-1 | `native_node` |
| electron/hole/mean ionization integral | vertex | 1 | 1 | `native_node` |
| doping, donor, acceptor, space charge | vertex | 1 | cm^-3 | `native_node` |
| SRH recombination | vertex | 1 | cm^-3 s^-1 | `native_node` |
| contact voltage/current/charge | contact | 1 | V, A, C | native contact observation |

The required electron and hole current-density vectors have two components on
both vertex and element support. Electron and hole alpha are present as scalar
vertex observations. No native element-centered alpha is claimed.

## CurrentPlot and documented runtime observations

The byte-identical CurrentPlot file contains:

- Cathode and Anode electron, hole, displacement, and total currents;
- electron and hole carrier-split avalanche-generation integrals;
- total avalanche-generation integral; and
- the runtime Tcl observation channel.

The documented runtime interface provides:

- element-local-edge `ReadCoefficient`;
- element-local-vertex `ReadMeasure`;
- element mobility and two-component E/QFP-gradient/current vectors; and
- independent carrier-split integration of native vertex generation with
  `ReadMeasure`.

These are native mesh coefficients, measures, fields, and integrals. A current
recomputed from vertex state plus `ReadCoefficient` remains
`operator_replay`; it is not a native directed-edge carrier current.

## Negative observation boundary

| Requested claim | Status |
|---|---|
| native directed-edge electron/hole current | `unsupported_native_edge` |
| `/Edge` ordinary Plot fields | parser rejection preserved from the frozen Phase A probe |
| `ReadFlux` as directed carrier edge current | unsupported location; not reinterpreted |
| native element-centered avalanche alpha | unsupported by this accepted contract |
| full nonlinear residual or Jacobian from ordinary Plot/CurrentPlot | unsupported |
| `CNormPrint` / `NewtonPlot` residual export | not probed because no locally documented minimal O-2018.06-SP2 syntax was available |

Box current reconstruction and SG replay are always labeled reconstructed.

## Provenance

The accepted A-root bundle hashes are:

| Input | SHA-256 |
|---|---|
| coarse7x3 mesh TDR | `11a358f95aba40d558588480b4a54ca751907322354f5907734dd4a795889ea7` |
| parameter file | `b4b3ebfdefba530f756f3855d43d7d587720689771d8badc747b61439ed42742` |
| exact process deck | `52b1b3bc52c7134ffb5ded8f0b273a8e5ecdec403358d105b2216a799c97564f` |
| runtime Tcl | `d4e9614d41706fb4a4dfa43fd5ed35513a33d7290be8c4203b19e0d4e6cf9730` |

The accepted byte-identical CurrentPlot SHA-256 is
`3e6de3d4925624b4b7ea8b264bd6f414bcccd1e3af6d62fdf33ebea5c3e563ae`.

Task 14 may now generate the declared exact high-bias lattice and paired
physics controls using this schema. Fine PN2D remains outside scope until the
coarse contract and all later gates authorize it.
