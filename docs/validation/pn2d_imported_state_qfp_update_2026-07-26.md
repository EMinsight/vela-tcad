# PN2D imported-state continuity residual and first-QFP-update audit

Date: 2026-07-26

Status: Task 7 evidence complete with failed Task 6 Jacobian and Task 7 cross-topology entry gates; Tasks 8-9 not entered.

## Decision

The typed outcome is `operator_improvement_without_qfp_causality`. Task 8 is not authorized for two independent reasons:

1. At the first material source-support departure (`-10 V`), Minimal6 worsens while coarse7x3 improves, so the same residual and first-update QFP error do not improve across topology classes.
2. The true source-isolated real-state `sg_avalanche` analytic/central-FD block has maximum nonzero relative difference `0.9640948767506723`, above the frozen `1e-8` gate.

The original Task 6 focused test compared full continuity matrices and could be dominated by SG transport. Review replaced it with an impact-only analytic difference and central FD of independently replayed source terms. A second review then caught that `rel_diff` used `max(1, fd_norm)`, turning small-source comparisons into absolute rather than relative errors. It now divides by the larger analytic/FD norm (and treats exact zero separately). The production block diagnostic also passes already isolated source/recombination pairs without subtracting the baseline twice. Fresh Van Overstraeten/Masetti/QFP evidence is therefore genuinely source-only but fails the required relative gate; no tolerance was changed.

Tasks 8 and 9 are recorded `not_entered_task6_jacobian_and_task7_causality_gates_failed`. The production default remains unchanged. `element_edge_sg_gss_laux` remains opt-in. No field, mobility, alpha, or geometry scale was fitted.

## Frozen inputs and topology provenance

The audit covers exact imported `psi`, electron QFP, and hole QFP at `-1`, `-10`, and `-20 V` for Minimal6 mirror/sketch and regenerated coarse7x3.

- All 27 coarse semiconductor node IDs and coordinates match the frozen probe exactly.
- All 32 coarse element IDs and local Tri3 vertex permutations match exactly.
- The six non-semiconductor probe vertices are excluded.
- Minimal6 and coarse contact-node sets and every imported field node ID are independently checked against each configured mesh.
- The runner, coarse log/mesh/doping/materials, every Minimal6 source config plus its directly referenced mesh/doping, and all imported Minimal6 fields are hash-sealed in each manifest.

The older 21-node/24-cell coarse export is not used. The matching frozen topology is 27 nodes/32 cells.

## Geometry scope

Task 6's constrained-obtuse RED proves that clipped circumcentric support was nonconservative. The opt-in helper now preserves nonnegative partials and exact total triangle area by uniform normalization. The patch does not have a native obtuse current/source oracle for the local vertex partition, so it is kept diagnostic and is not described as general operator parity.

Minimal6 and regenerated coarse7x3 are right-triangle controls without the pre-fix area-closure defect. The obtuse normalization therefore cannot itself explain their QFP updates.

## Frozen configurations

Production triangle:

- Van Overstraeten;
- QFP-gradient driving force;
- `cell_reconstructed` plus `gss_logistic`;
- `triangle_gss_gradqf_truncated`.

Opt-in candidate:

- Van Overstraeten;
- QFP-gradient driving force with `cell_gradient`;
- `element_edge_sg_gss_laux`;
- `element_vertex_box_measure`.

Both use unchanged Vela Masetti mobility with QFP-gradient high-field driving. No A/B or geometry scale is present. Box edge fluxes are a native-mobility reconstruction diagnostic, not a native Sentaurus directed edge-current observation and not a third first-update branch.

## Sealed v2 evidence

Two fresh roots are generated independently:

- `build-release/pn2d-imported-state-qfp-update-20260726-a`;
- `build-release/pn2d-imported-state-qfp-update-20260726-b`.

| Evidence | Rows |
|---|---:|
| topology/contact/state gate | 9 |
| residual decomposition | 468 |
| carrier-only and coupled first QFP updates | 600 |
| analytic/central-FD Jacobian blocks | 90 |
| paired avalanche-off and SRH-off controls | 936 |
| independently reproducible causality groups | 54 |

The v2 verifier does not import analyzer calculations. It independently checks the complete 117-JSON generated configuration lattice, frozen mobility/impact/contact/tolerance contracts, exact coarse `psi/eQFP/hQFP` values at all three biases, topology/contact mapping, the complete expected input-hash key set, paired controls, residual closure, causality directions, Jacobian gates, authorization, and typed outcome.

Its expected status is `pass: false` with the typed error `source-specific Jacobian gate failed`. A/B sealed outputs are byte-identical, complete provenance hashes validate, and both roots pass the full configuration/topology/path-binding checks. The false verifier result is the intended enforcement of the scientific gate, not a verifier malfunction.

| Gate | Result |
|---|---:|
| A/B sealed hashes | identical |
| input hashes | verified |
| residual closure maximum relative | `1.6654229802532582e-16` |
| boundary branch maximum absolute difference | `0` |
| avalanche-off branch maximum absolute difference | `0` |
| nonzero source Jacobian maximum relative | `0.9640948767506723` (fail) |
| near-zero source Jacobian maximum absolute | `4.6780474445207e-17` (pass) |
| cross-topology first-material causality | `false` |
| Task 8 authorized | `false` |

## Residual and first-update causality

Candidate/production L2 ratios below one improve.

| Topology | Bias | electron residual | hole residual | carrier-only electron | carrier-only hole | coupled electron | coupled hole |
|---|---:|---:|---:|---:|---:|---:|---:|
| Minimal6 mirror/sketch | -1 | `1.000000` | `1.000000` | `1.000000` | `1.000000` | `1.000000` | `1.000000` |
| Minimal6 mirror/sketch | -10 | `1.027156` | `1.027976` | `1.023601` | `1.024308` | `1.023601` | `1.024308` |
| Minimal6 mirror/sketch | -20 | `3.409014e-5` | `3.373837e-5` | `0.0297334` | `0.0294365` | `0.0297335` | `0.0294366` |
| coarse7x3 | -1 | `1.000000` | `1.000000` | `1.000000` | `1.000000` | `1.000000` | `1.000000` |
| coarse7x3 | -10 | `1.099286e-7` | `1.099887e-7` | `0.00158892` | `0.0140316` | `0.00106766` | `0.00924058` |
| coarse7x3 | -20 | `0.0180704` | `0.0181181` | `0.0120429` | `0.00857828` | `0.00891144` | `0.00639328` |

Exact equality at low-signal `-1 V` is classified `equal`, not as a failure. The independent rule identifies the first non-equal bias for every topology (`-10 V`) and then compares directions. Coarse improves in every residual/update group; Minimal6 worsens in every group. This is the decisive cross-topology contradiction.

The first maximum term departures are all avalanche terms:

- coarse7x3: `-10 V`, hole, node 10;
- Minimal6 mirror: `-10 V`, hole, node 1;
- Minimal6 sketch: `-10 V`, hole, node 5.

Paired avalanche-off branches are identical. Paired SRH-off branches retain the source-operator difference. Boundary rows are identical. This locates the first discriminator in avalanche source/current support, not SRH, contacts, or changed mobility configuration.

Because the initial state is exactly Sentaurus QFP, the Sentaurus-to-Vela reference vector is zero. Update direction relative to that zero vector is recorded `undefined_zero_reference`; `abs(delta_qfp)` is retained only as the post-update error magnitude.

## Residual units and scales

Normalized carrier terms come directly from C++ diagnostics. Physical columns multiply by the assembler continuity unit scale `C0*D0`. The analyzer-derived maximum of incident absolute SG flux, SRH, and avalanche is explicitly named `diagnostic_incident_term_scale`; it is not claimed to be an exported solver row scale. The runner currently does not export a separate solver row-scaling factor for this probe.

## Review disposition and Task 10 decision

Independent scientific and code reviews both found the masked Jacobian test blocking. The response is:

- source-only focused central FD added;
- real-state source/recombination block audit changed from subtracting full forward-FD matrices to direct isolated-term central FD;
- paired controls added for both variants;
- circular authorization removed;
- verifier independently derives the gate/outcome, validates both full configuration/provenance lattices, and intentionally fails the rejected source-Jacobian gate;
- input/topology/contact/state provenance sealed and independently checked;
- zero-reference QFP direction and first causal node/term recorded;
- analyzer-defined scale renamed honestly.

The applicable production decision remains: keep the candidate opt-in as a diagnostic; no default change and no self-consistent or fine-BV claim.

Generated `build-release` roots are not committed.
