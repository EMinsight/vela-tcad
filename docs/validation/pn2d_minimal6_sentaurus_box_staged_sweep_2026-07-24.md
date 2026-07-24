# PN2D Minimal6 40-state Sentaurus box-current staged sweep

Date: 2026-07-24

## Scope

Extend the closed mirror/-1 V replacement experiment to the exact 40-state
lattice:

- topologies: `mirror` and `sketch`;
- reverse biases: -1 V through -20 V;
- support per state: 6 nodes, 4 triangles, 9 global edges, and 12
  element-local edges;
- replacement chain:
  `QFP -> carrier density -> element mobility -> box geometry`.

No production formula in `include/` or `src/` was modified.

## Input and support gates

The replay uses:

- Sentaurus and Vela node potential/QFP/density observations for all 40
  states;
- native Sentaurus element electron and hole mobility for all 160
  state-elements;
- Vela velocity-fix fixed-state edge current and inferred production edge
  mobility;
- the final Sentaurus terminal `.plt` current for every state; and
- the documented 2-D box coefficient and carrier-current equations validated
  by the preceding single-state audit.

Sentaurus element fields are stored in `region_cell_order`, which is not the
same as Vela triangle-id order for `sketch`. The mapping was independently
recovered by matching native element electric field to the P1 field
`-grad(psi)`:

| Topology | Vela triangle -> Sentaurus region-cell index | Maximum relative field residual |
|---|---|---:|
| mirror | 0->0, 1->1, 2->2, 3->3 | 2.959312e-16 |
| sketch | 0->0, 1->3, 2->2, 3->1 | 2.959312e-16 |

Using the uncorrected `sketch` order produced a false 4.04% carrier-terminal
error. Applying the field-verified permutation reduces the total-terminal
gate to the same approximately 1e-7 level as the original mirror/-1 V probe.

## Forty-state pooled result

There are five nonzero reference edges per carrier and state, giving 200
valid carrier-edge samples per carrier and stage. Exact-zero supports remain
typed and are excluded from dex statistics.

| Stage | Electron median / p95 / max (dex) | Electron sign | Hole median / p95 / max (dex) | Hole sign |
|---|---:|---:|---:|---:|
| Vela baseline | 5.427222 / 6.438195 / 6.718769 | 0.81 | 5.366981 / 6.439125 / 6.727331 | 0.81 |
| Sentaurus QFP | 5.398867 / 6.560457 / 7.237132 | 1.00 | 5.315772 / 6.205333 / 6.881592 | 1.00 |
| + Sentaurus density | 0.062732 / 0.961922 / 1.063451 | 1.00 | 0.057630 / 0.616202 / 0.710167 | 1.00 |
| + Sentaurus element mobility | 0 / 0 / 0 | 1.00 | 0 / 0 / 0 | 1.00 |
| + box geometry | 0 / 0 / 0 | 1.00 | 0 / 0 / 0 | 1.00 |
| recomputed-density control | 0.062727 / 0.961927 / 1.063455 | 1.00 | 0.057626 / 0.616206 / 0.710171 | 1.00 |

`mirror` and `sketch` give the same statistics after the cell-order
permutation, within floating-point roundoff.

## Paired incremental contribution

Values below are the per-edge paired median of
`previous abs-dex error - current abs-dex error`.

| Replacement step | Electron (dex) | Hole (dex) |
|---|---:|---:|
| Vela baseline -> Sentaurus QFP | +0.014939 | +0.017077 |
| QFP -> QFP + density | +5.344047 | +5.244281 |
| density -> element mobility | +0.062732 | +0.057630 |
| mobility -> geometry | 0.000000 | 0.000000 |

QFP replacement makes the direction correct on all 400 nonzero
carrier-edge samples but barely changes the magnitude. Carrier density is the
dominant closure step. Element mobility removes the remaining edge-dependent
residual, and geometry contributes zero.

## Bias trend

The table shows state-level median error for one topology; the other topology
is numerically identical after support mapping.

| Bias | Carrier | Vela baseline | Sentaurus QFP | + density | Full |
|---:|---|---:|---:|---:|---:|
| -1 V | electron | 5.156389 | 4.975977 | 0.093137 | 0 |
| -1 V | hole | 4.954712 | 4.781812 | 0.039681 | 0 |
| -5 V | electron | 5.437897 | 5.386237 | 0.041703 | 0 |
| -5 V | hole | 5.361590 | 5.310519 | 0.041990 | 0 |
| -10 V | electron | 5.485280 | 5.457856 | 0.065095 | 0 |
| -10 V | hole | 5.452378 | 5.425090 | 0.071364 | 0 |
| -15 V | electron | 5.508590 | 5.489893 | 0.083159 | 0 |
| -15 V | hole | 5.492349 | 5.473696 | 0.086520 | 0 |
| -20 V | electron | 6.074732 | 6.059103 | 0.094561 | 0 |
| -20 V | hole | 6.069584 | 6.053963 | 0.096749 | 0 |

The largest mobility-sensitive samples occur on the central 1-5 edge at
-20 V:

| Carrier | Reference current magnitude | Vela-mobility candidate magnitude | Error |
|---|---:|---:|---:|
| electron | 1.625904e-20 A/um | 1.881681e-19 A/um | 1.063451 dex |
| hole | 1.546089e-20 A/um | 7.932346e-20 A/um | 0.710167 dex |

These maxima are very small central-edge currents. They explain the p95/max
mobility tails but do not dominate terminal total current.

## Density self-consistency

For all 240 state-node samples, density recomputed from Sentaurus
electrostatic/QFP potentials and the effective intrinsic density recovered
from the Vela BGN state agrees with exported Sentaurus density to at most
`4.426180e-6 dex`.

The pooled current medians of the recomputed-density control differ from
direct Sentaurus-density replacement by approximately 4e-6 dex. Thus the
dominant density step is reproduced by the Vela BGN transformation and is
not an arbitrary fitted substitution.

## Closure and conservation gates

| Gate | Result |
|---|---:|
| exact states checked | 40 |
| stage edge samples checked | 4,320 |
| final nonzero carrier-edge values checked | 400 |
| final replay/reference maximum relative error | 0.0 |
| Vela baseline production/replay maximum relative difference | 5.295713e-14 |
| cell-mapping maximum relative field residual | 2.959312e-16 |
| carrier-specific terminal maximum relative error | 2.318564e-3 |
| total terminal maximum relative error | 1.137692e-7 |
| internal total-current KCL maximum relative error | 3.163255e-9 |
| independent verification failures | 0 |

The carrier-specific maximum occurs for the small Cathode hole current at
-20 V. Summing electron and hole currents restores the physically relevant
terminal conservation gate to `1.137692e-7`. The total-current KCL result
remains below `1e-8`.

## Scientific conclusion

The single-state conclusion is stable over all 40 exact states:

1. QFP replacement fixes current direction but not magnitude.
2. Density made self-consistent with the replaced potentials removes about
   5.24-5.34 dex and is the dominant closure step.
3. Sentaurus element mobility removes the remaining median 0.06 dex, with a
   high-bias central-edge tail up to 1.063 dex for a very small current.
4. The box/cotangent geometry contributes exactly zero error.
5. The reconstructed box current closes total terminal current and internal
   total-current KCL at approximately 1e-7 and 1e-9, respectively.
6. The apparent earlier `sketch` discrepancy was a cell-order mapping error,
   not a transport formula difference.

This establishes a deterministic Sentaurus-operator reconstruction across the
full lattice. It remains a reconstruction rather than a native directed-edge
observation, so the result supports diagnosis and counterfactual replacement
but does not alone authorize changing production current formulas.

## Evidence

- Evidence root:
  `build-release/pn2d-minimal6-sentaurus-box-staged-sweep-20260724-a`
- Full stage lattice: `stage_edge_samples.csv`
- Pooled summary: `stage_summary.csv`
- Paired contributions: `paired_contributions.csv`
- Per-state trend: `state_summary.csv`
- Density control: `density_recompute_control.csv`
- Mobility comparison: `mobility_comparison.csv`
- Geometry: `geometry_coefficients.csv`
- Sentaurus/Vela cell mapping: `cell_mapping.csv`
- Terminal closure: `terminal_closure.csv`
- Total-current KCL: `total_current_kcl.csv`
- Baseline production replay: `baseline_operator_crosscheck.csv`
- Independent verification: `independent_verification.json`

The independent verifier regenerated all pooled summaries and paired
contributions, independently recomputed contact and internal-node
conservation from the final edge table, checked the exact 40-state lattice,
and reported zero failures.
