# PN2D BV atomic SG/Laux + mixed-Voronoi default acceptance

Date: 2026-08-01

Status: passed; ready for the independent-review decision.

## Candidate and scope

The PN2D BV template version 3 renders one atomic default profile:

- `impact_ionization.current_approximation=element_edge_sg_gss_laux`;
- `impact_ionization.source_mapping_mode=element_vertex_box_measure`;
- `impact_ionization.cell_reconstructed_midpoint_density=bernoulli`;
- `mesh_geometry.node_volume_policy=mixed_voronoi`;
- `mesh_geometry.require_non_obtuse=true`.

The named `legacy_cell_reconstructed` rollback atomically restores the previous
three avalanche-support values together with `barycentric` and
`require_non_obtuse=false`. The global C++ defaults remain barycentric and do
not require a non-obtuse mesh. The PN2D IV template is unchanged and retains
impact ionization off.

## Prospective contract and bindings

The frozen contract is
`docs/validation/contracts/pn2d_node_volume_policy_atomic_default_acceptance_v2.json`
with SHA-256
`f7b2c356876892ab0293c244c0446ace3aab73ca2669692095600318b1dcc12c`.

Both M0 and M2 were rendered from template defaults without an avalanche
profile override. Every off, IIC, and on execution records
`current_support.origin=base_config`. The acceptance evaluator replayed both
render manifests and configs exactly and verified their hashes. It also bound
the unchanged PN2D IV template and the atomic rollback.

## Results

| Metric | M0 | M2 |
|---|---:|---:|
| Exact lattice, each branch/run | 29/29 | 29/29 |
| Independent run IV/state hashes | identical | identical |
| Off current RMSE vs Sentaurus | 0.0000888 dex | 0.0000387 dex |
| Avalanche-on 28-point RMSE | 0.001766 dex | 0.001923 dex |
| Avalanche-on maximum error | 0.003749 dex | 0.004647 dex |
| Knee RMSE | 0.002810 dex | 0.003051 dex |
| Absolute V_break error | 0.001 V | 0.001 V |
| Nonmonotonic intervals | none | none |
| Maximum on continuity closure ratio | 0.002200 | 0.00000192 |

IIC and avalanche-off have identical IV and state hashes at every requested
bias, confirming that postprocess-only impact ionization does not feed back
into the continuity equations.

The 201-point forward-IV guard passed. The Sentaurus anchor median relative
error is 0.2700%, the maximum is 0.4066%, and the maximum degradation relative
to barycentric is `2.42e-14`.

The final Release suite passed 512/512 tests. The release log SHA-256 is
`d9d4ea6328a4b21f41a7091e6014c6711cfa332a9bc69823ad44a117cb6754d3`.

## Evidence

Generated evidence is under
`build-release/pn2d-node-volume-atomic-default-acceptance-v2-20260801/`.
The final aggregate is `acceptance.json` with SHA-256
`b08b7dd7b5616cd7ad29697189fd2c7e4fd04d36f21a6d0016f21debea81c450`.

The acceptance applies only to the qualified non-obtuse PN2D BV M0/M2 Tri3
mesh family. It is not evidence for arbitrary obtuse meshes, other devices, or
a global solver default change.
