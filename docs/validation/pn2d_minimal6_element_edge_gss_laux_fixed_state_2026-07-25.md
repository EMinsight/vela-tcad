# PN2D Minimal6 element-edge GSS/Laux fixed-state audit

Date: 2026-07-25

## Scope

This audit evaluates an explicit opt-in avalanche source path on six exact
Minimal6 states:

- topologies: `mirror` and `sketch`;
- reverse biases: -1 V, -10 V, and -20 V;
- support per state: 6 nodes, 4 triangles, 9 global edges, and 12
  element-local directed edges;
- imported state: Sentaurus electrostatic potential, electron quasi-Fermi
  potential, and hole quasi-Fermi potential;
- recomputed state: Vela Old-Slotboom carrier density, Vela mobility,
  variable-intrinsic-density SG edge flux, GSS/Laux cell current vector,
  Van Overstraeten coefficient, and element-vertex source.

The Sentaurus directed-edge values remain a documented box-operator
reconstruction, not a native directed-edge observation. The new Vela path is
not the production default.

## Input correction

An initial run accidentally used the Vela self-consistent baseline QFP state.
At -20 V, its central electron QFP is approximately -10.009 V, whereas the
Sentaurus QFP is approximately -9.642 V. That mixed-state input generated a
false 4.5 to 6 dex current discrepancy.

The accepted run uses the original Sentaurus field export and replaces its
carrier-density columns with the independently verified Vela Old-Slotboom
density recomputation. The preceding 40-state density audit showed a maximum
recomputed-versus-exported Sentaurus density difference of 4.426180e-6 dex.

## Discrete operator

For each triangle and carrier, the diagnostic path:

1. evaluates all three directed variable-ni SG edge particle fluxes;
2. represents all three edge currents in the general operator; on this
   right-triangle mesh the zero-weight hypotenuse current is numerically
   inactive;
3. reconstructs one cell current vector with the GSS/Laux pair solve and box
   weighting;
4. evaluates the Van Overstraeten coefficient from the cell P1 electric field
   magnitude;
5. multiplies alpha by the reconstructed current magnitude;
6. distributes the generation with the three exact element-vertex box
   measures.

All 24 right-triangle hypotenuse partial volumes are exactly zero. They are
geometric integration weights. The implementation keeps a general three-edge
record, but the Minimal6 GSS/Laux weighted vector is exactly unchanged if a
zero-weight hypotenuse current is perturbed.

## Six-state comparison

The independently recomputed pooled errors are:

| Quantity | Electron median / max (dex) | Hole median / max (dex) |
|---|---:|---:|
| directed edge current | 0.069064 / 1.063346 | 0.047671 / 0.710160 |
| reconstructed cell current vector | 0.053081 / 0.220792 | 0.045517 / 0.084771 |
| Van Overstraeten alpha | 1.02e-12 / 1.55e-9 | 1.68e-12 / 2.56e-9 |
| accumulated node avalanche source | 0.058771 / 0.232730 | 0.047670 / 0.090251 |
| whole-device carrier source integral | 0.073549 / 0.197052 | 0.062718 / 0.076738 |

The large maximum directed-edge errors are the already identified small
central-edge mobility tails. They do not represent the pooled current or
source behavior.

Directly comparing one Vela element source with one Sentaurus
element-vertex value is not the continuity-vector contract. Sentaurus exports
a nodal generation value that is redistributed into every adjacent element
by `ReadMeasure`; Vela first forms a cell generation and then assembles its
three element-vertex contributions. The valid comparison is therefore the
sum over all adjacent elements at each physical node, followed by the
whole-device integral.

## Driving-force result

Using the cell quasi-Fermi gradient for alpha produced approximately 1.0 dex
electron and 1.7 dex hole pooled alpha errors, with much larger low-field
tails. Replacing only the alpha driving field by the P1 electric field
`|-grad(psi)|` reduces the maximum alpha errors to 1.55e-9 and 2.56e-9 dex.

This is a direct fixed-state result:

- the Vela and Sentaurus Van Overstraeten coefficient formulas agree;
- the Sentaurus default `Avalanche(VanOverstraeten)` coefficient in this
  experiment follows the electric-field support;
- the QFP variables still determine carrier density and SG current;
- QFP gradient is not the matching coefficient driver for this default
  Sentaurus branch.

## Verification

The standalone verifier checked:

- exactly 72 element-vertex rows and 6 states;
- exactly 24 zero hypotenuse partial volumes;
- alpha/current/source error gates;
- the identity
  `qg = alpha * |J| * element_vertex_measure`;
- accumulated physical-node sources;
- per-carrier whole-device source integrals.

The source identity maximum relative error is 4.04e-16 and the verifier
reported zero failures.

## Scientific conclusion

The fixed-state evidence supports the following decomposition:

1. imported Sentaurus potential/QFP plus Vela Old-Slotboom statistics
   reproduces the carrier state;
2. the remaining directed-edge and reconstructed-vector current errors are
   at the known Vela-versus-Sentaurus element-mobility scale;
3. GSS/Laux reduces exactly to its two positive-weight edges on Minimal6;
   a separate acute-scalene regression establishes that all three edges are
   active when all three partial volumes are positive;
4. electric-field-driven Van Overstraeten alpha agrees with Sentaurus to
   numerical precision;
5. accumulated node and integrated avalanche source errors reduce to the
   same mobility/current scale.

The evidence rejects changes to the Van Overstraeten coefficient formula. It
supports continuing with the opt-in element-edge current support and
electric-field coefficient driver, while retaining a general three-edge
operator for non-right-triangle meshes, before any production
default decision.

## Evidence

- Generator:
  `scripts/diagnose_pn2d_minimal6_element_edge_gss_laux_fixed_state.py`
- Independent verifier:
  `scripts/verify_pn2d_minimal6_element_edge_gss_laux_fixed_state.py`
- Evidence root:
  `build-release/pn2d-minimal6-element-edge-gss-laux-fixed-state-20260725`
- Detailed comparison: `fixed_state_comparison.csv`
- Accumulated node source: `node_source_comparison.csv`
- Whole-device source: `state_source_summary.csv`
- Independent result: `independent_verification.json`
