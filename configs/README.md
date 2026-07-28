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
