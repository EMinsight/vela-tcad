# PN2D avalanche-off SRH paired mesh convergence

Date: 2026-07-28

## Frozen inputs

- Starting commit contains `fa1c343`.
- BV mobility doping basis remains `net_doping`.
- Forward-IV mobility doping basis remains
  `cell_reconstructed_total_impurity`.
- Impact ionization is disabled; SRH and Old Slotboom BGN remain enabled.
- Sentaurus base source is
  `reference_tcad/pn2d_sentaurus2018_coarse7x3/source`.
- The same TDR-derived physical mesh and doping files are used by Vela at
  each level.

Ignored artifacts:

`build-release/pn2d-bv-off-srh-mesh-matrix-coarse-baseline-20260728/`

## Mesh and run evidence

| Level | Nodes | Triangles | Minimum edge (um) | Vela/Sentaurus points | log-current RMSE (dex) | max closure |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| M0 | 27 | 32 | 0.25 | 21/21 | `6.997e-5` | `7.482e-7` |
| M1 | 43 | 60 | 0.125 | 21/21 | `2.533e-2` | `1.684e-6` |
| M2 | 115 | 191 | 0.0625 | 21/21 | `1.118e-2` | `6.890e-7` |

All levels retain the `0..2 um` by `0..0.5 um` bounding box and `1 um2`
area. Anode and cathode remain on `x=0` and `x=2 um`.

The paired error does not improve monotonically from M0 to M2: M0 is already
the closest pair. Consequently the required evidence for
`mesh_resolution_difference` is absent.

## Failed mesh-independence gates

The M1-to-M2 positive-generation integral changes are:

| Bias (V) | Vela | Sentaurus |
| ---: | ---: | ---: |
| -1 | `49.98%` | `49.98%` |
| -5 | `72.19%` | `72.19%` |
| -10 | `11.02%` | `10.92%` |
| -15 | `12.98%` | `12.81%` |
| -20 | `23.05%` | `59.50%` |

Every anchor exceeds the `<2%` two-finest-level requirement.

The node-control-volume total-impurity dose changes by `-5.56%` on M1 and
`-8.33%` on M2 relative to M0, exceeding the `<0.1%` gate. The sealed SDE
places both constant P and N profiles on the `x=1 um` junction nodes. Their
continuum overlap has zero measure, but the imported nodal representation
assigns it a mesh-dependent control volume; junction refinement therefore
changes the discrete total-impurity dose and source support.

## Task 4 decision

Task 4 is not accepted. Both simulators converge and remain mutually close on
M0 and M2, but the study has not produced a mesh-independent source integral
or dose-preserving nested sequence. The RMSE also fails the monotonic
mesh-root-cause criterion.

Per the plan, stop before Task 5. No SRH production correction or default
change is authorized from this matrix. A follow-up must first define a
dose-preserving junction profile/support and obtain two finest levels with
`<2%` source change. That requires new Sentaurus inputs beyond the explicitly
authorized M0/M1/M2 upload set.
