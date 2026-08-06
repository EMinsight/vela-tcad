# BVmethods NMOS fixed-peak generation-support growth (2026-08-05)

## Scope and invariants

The accepted Fermi-Dirac SG states, mobility, Eparallel driving force,
Van Overstraeten coefficients, and avalanche source mapping were not changed.
This audit only compares the spatial support of electron avalanche generation
from 6.4 to 7.0 V.

Sentaurus node generation is integrated with the exact semiconductor P1
measure. Vela electron edge source is symmetrically projected to the same
nodes and divided by that same P1 measure. The projection conserves the Vela
electron-source integral at every bias. A second Vela-native edge-support
calculation is reported to expose threshold sensitivity to the projection.

## Threshold support

Absolute P1-node support-area ratios are:

| peak threshold | Vela/Sentaurus at 6.4 V | Vela/Sentaurus at 7.0 V |
|---:|---:|---:|
| 10% | 1.17638 | 1.00122 |
| 30% | 1.17512 | 1.09915 |
| 50% | 1.19552 | 1.25828 |
| 80% | 1.51554 | 1.25813 |

The common-node area-growth factors from 6.4 to 7.0 V are:

| peak threshold | Sentaurus factor | Vela P1-projected factor | Vela native-edge factor |
|---:|---:|---:|---:|
| 10% | 1.17522 | 1.00024 | 1.15839 |
| 30% | 1.08818 | 1.01783 | 1.06087 |
| 50% | 1.03786 | 1.09235 | 1.16016 |
| 80% | 1.15070 | 0.95525 | 1.14589 |

The exact threshold areas depend on edge-to-node projection, especially at
80%, but both representations give the same robust localization: Vela grows
too slowly in the 10--30% generation shoulder, while the 50% core is not too
narrow and the native 80% core growth is already aligned.

## Radial cumulative distribution

Using the Sentaurus peak node as a common center, the cumulative-source
fractions agree within 2.3 percentage points at every tested radius from
0.0025 to 0.2 um. At 7.0 V, Vela is 0.75 percentage points low inside
0.005 um, 2.27 points high inside 0.04 um, and 0.35 points low inside
0.08 um. Therefore the remaining support error is not a simple isotropic
radial-width error. It is an anisotropic shoulder redistribution along the
drain-body high-field corridor.

## Edge and cell localization

The ten largest negative-growth edges contain 57.8% of all qualified negative
edge-growth deficit. Their midpoints occupy the compact corridor

`x=0.07945--0.09711 um`, `y=0.00459--0.01868 um`,

downstream from the common peak near `(0.06621, 0.00137) um`.

The largest edge deficit is edge 2672, nodes 565--573, midpoint
`(0.083867, 0.009289) um`:

| quantity, 6.4 to 7.0 V | Sentaurus | Vela |
|---|---:|---:|
| endpoint generation growth | 2.06908 | 1.86367 |
| alpha growth | 1.16721 | 1.10531 |
| electron-current growth | 1.70142 | 1.68611 |
| Vela Eparallel growth | n/a | 1.09544 |

Its missing 7.0 V source relative to Sentaurus local growth is
`2.48976e17 m^-1 s^-1`. The current-growth difference is below 1%, whereas
the alpha-growth difference is about 5.3%; the deficit is therefore mainly the
shoulder Eparallel/alpha evolution, not the local electron-current slope.

The largest cell deficit is triangle 1628, nodes 564/573/565, centroid
`(0.082396, 0.008898) um`, containing edges 2672 and 2673. The next dominant
cells 1629, 1644, 1646, 1626, and 1642 remain in the same junction corridor.

## Decision

The fixed peak model is retained. The next implementation investigation should
freeze carrier currents and compare the carrier-specific Eparallel recovered
on edges 2663, 2666, 2672, 2675, 2681, and adjacent cells against Sentaurus
alpha-derived local driving fields. No global alpha coefficient or mobility
refit is justified by these data.

## Artifacts

- implementation: `scripts/audit_bvmethods_nmos_generation_support_growth.py`;
- output root:
  `build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela_validation/qf_vector_generation_support_growth_20260805`;
- threshold ledgers: `threshold_support_compare.csv` and
  `vela_native_edge_threshold_support.csv`;
- radial ledger: `radial_cumulative_source_compare.csv`;
- localization ledgers: `top_node_growth_deficits_6p4_7p0.csv`,
  `top_edge_growth_deficits_6p4_7p0.csv`, and
  `top_cell_growth_deficits_6p4_7p0.csv`.
