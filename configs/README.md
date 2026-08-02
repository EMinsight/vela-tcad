# Vela simulation configuration assets

`templates/` contains versioned PN2D starting points. Forward IV and reverse BV
are separate because their mobility, impact-ionization, solver, and continuation
settings are intentionally different.

Use `scripts/generate_pn2d_config.py` instead of copying and editing a previous
run's JSON. The generator:

- accepts only parameters declared by the selected template;
- rejects invalid types, absolute paths by default, inconsistent sweep
  directions, and unqualified IV/BV physics combinations;
- writes stable, sorted JSON plus a separate `.manifest.json`;
- leaves the runnable Vela JSON free of template-only metadata.

Relative paths are interpreted from the generated simulation JSON's directory.
The schema in `schema/vela-simulation.schema.json` describes the common
machine-readable contract for rendered configurations.

## PN2D BV profiles

`pn2d_bv` template version 3 defaults to the qualified atomic
`element_edge_sg_gss_laux` profile. The generator resolves it to:

- `impact_ionization.current_approximation=element_edge_sg_gss_laux`;
- `impact_ionization.source_mapping_mode=element_vertex_box_measure`;
- `impact_ionization.cell_reconstructed_midpoint_density=bernoulli`;
- `mesh_geometry.node_volume_policy=mixed_voronoi`; and
- `mesh_geometry.require_non_obtuse=true`.

Use `avalanche_current_support_profile=legacy_cell_reconstructed` to restore
the complete cell-reconstructed, barycentric rollback configuration. Partial
mixes of the two profiles fail validation. This policy applies only to the
PN2D BV template; the PN2D IV template keeps impact ionization disabled and
the global parser default remains barycentric.
