# Charon-derived Sentaurus 2018 Schottky reference

This fixture translates Charon's two-dimensional `n_diode` Schottky-contact
test into Sentaurus SDE/SDevice and then compares Vela on the exact imported
Sentaurus mesh. The source benchmark is a 1 um by 1 um, uniformly doped
`1e16 cm^-3` n-silicon slab with a 4.75 eV Schottky anode and an Ohmic
cathode at 300 K.

The Sentaurus O-2018.06-SP2 reference completed the full 0--1 V sweep. The
small committed curve contains the converged equilibrium point and all forward
points; the TDR, PLT, and logs remain under the ignored
`build-release/reference_tcad/schottky_charon_sentaurus2018/raw` tree and are
sealed by hashes in `schottky_charon_sentaurus2018_reference.json`.

Vela's pre-existing `dirichlet_barrier` Gummel model did not converge at 0 V
on the imported 697-node mesh after 150 iterations. The implemented minimum
parity feature is therefore the boundary equation used by the source deck:
electrostatic Schottky barrier pinning plus independent electron/hole
thermionic Robin fluxes in coupled Newton. No image-force lowering, tunnelling,
series resistance, AC, or high-field model is enabled.

The first acceptance range is intentionally bounded to 0.01--0.54 V. It
passes monotonic trend and a maximum 0.5-decade current-shape tolerance over
all 14 Sentaurus points in the range. A diagnostic extension reached 0.5625 V
and then stalled at 0.563125 V in the electron continuity block; 0.54--1 V is
kept as a separate numerical-continuation milestone instead of expanding this
feature slice.

Run the bounded Vela deck with:

```text
build/vela_example_runner.exe --config reference_tcad/schottky_charon_sentaurus2018/vela/simulation_iv.json
```

The deck uses `scaling.mode = unit_scaling`, so its thermionic velocity values
are in cm/s despite the stable historical `_m_per_s` field spelling.
