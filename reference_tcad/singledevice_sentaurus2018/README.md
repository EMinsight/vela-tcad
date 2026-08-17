# Sentaurus 2018 SingleDevice nMOS reference

This fixture freezes the `GettingStarted/sdevice/SingleDevice` application-library
example from Sentaurus O-2018.06-SP2. The source decks are preserved verbatim.
Only the small, 21-point `Id-Vg` curves and their provenance are committed.

The completed baseline run used the original density-gradient, mobility, BGN,
and SRH equations. Large TDR, PLT, log, and restart artifacts are staged under
the ignored `build-release/reference_tcad/singledevice_sentaurus2018` tree.

The two output branches both start from the saved equilibrium state:

- `singledevice_idvg_lin_reference.csv`: `Vds=0.1 V`
- `singledevice_idvg_sat_reference.csv`: `Vds=1.1 V`

The original Save/Load dependency graph and both self-consistent Id-Vg branches
now pass with Vela's density-gradient outer coupling.  The current closeout,
derived metrics, low-current rule, and three-bias field/KCL audit are recorded
in `singledevice_validation_20260817.json`.  The older
`singledevice_validation_20260812.json` remains an immutable historical
snapshot of the previously partial state.

Reproduce the staged Vela workflow with `scripts/run_singledevice_workflow.py`,
then apply `scripts/close_singledevice_validation.py` and
`scripts/analyze_singledevice_closeout_fields.py`.  The gate-sweep numerical
settings are branch-specific: the linear branch uses a 0.1 V quasi-Fermi update
limit, while the already proven saturation branch retains 0.025 V.  This is a
solver-path setting only; no physics model was added for the closeout.

Sentaurus field exports can be converted to a Vela restart state with
`scripts/sentaurus_fields_to_restart.py`.  The converter clears Sentaurus's
insulator-only quasi-Fermi placeholders because Vela pins those algebraic rows
to zero.  The restart format preserves the electron quantum potential so that
frozen-potential and fully coupled checks can be separated.

Generate a credential-free run plan with:

```text
python scripts/run_singledevice_sentaurus_vm.py --run-id singledevice_check
```

A live run adds `--live`; authentication is supplied
by the caller's SSH configuration or agent and is never stored in the fixture.
