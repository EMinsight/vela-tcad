# SingleDevice Eq. 231 impurity-BGN and mass split (2026-08-14)

## Outcome

The Eq. 231-only OldSlotboom drive now depends on total impurity concentration
and no longer lets the inversion-layer carrier density dominate the bandgap
narrowing.  The quantum gradient coefficient mass is independently
configurable from the DOS mass, with legacy fallback to the DOS mass when the
new parameter is absent.  The generic drift-diffusion BGN call sites were not
changed.

For the Sentaurus 2018 SingleDevice material contract, the DOS mass remains
`1.0618016171622988` and the recovered coefficient mass is
`1.0906506732296395` for Silicon and PolySilicon.

## Fixed-state three-node oracle

The `sentaurus_box` linear and saturation probes were regenerated from the
same checked-in material file.  An old generated saturation source had pointed
to `materials_sentaurus2018_updated.json`; the generator now forces both
endpoints to use `materials_sentaurus2018.json`.

| endpoint | node 2101 residual | node 2100 residual | node 2044 residual | gate |
|---|---:|---:|---:|---|
| linear | `-5.18992e-7` | `-5.67389e-7` | `-5.96338e-7` | pass |
| saturation | `-4.95632e-7` | `-5.46068e-7` | `-5.75921e-7` | pass |

All six rows are below the required `1e-5` absolute residual.  This confirms
the cross-state impurity-BGN and coefficient-mass hypothesis without a
node-specific fit.

## Self-consistent endpoint gate

Both imported endpoint runs complete their classical Newton solve and their
quantum inner solve, but the one-step self-consistency gate still fails:

| endpoint | classical Newton iterations | quantum inner iterations | inner result | raw quantum change | result |
|---|---:|---:|---|---:|---|
| linear, Vg=2.2 V, Vd=0.1 V | 5 | 21 | converged | `1.000394 V` | fail |
| saturation, Vg=2.2 V, Vd=1.1 V | 7 | 21 | converged | `1.000394 V` | fail |

The remaining global fixed-state maximum is node 1793, with residual
`1.59751409` at the linear endpoint and `1.59751482` at saturation.  Node 1793
is an SiO2/Nitride shared corner (two `R.PolyReox` cells and four `R.Spacer`
cells), not one of the ordinary-Silicon Formula-0 rows closed by this task.
Consequently the 21-point Id-Vg curves remain gated.  The next localization
target is the insulator-insulator region-side trace/reaction closure at that
corner.

## Windows runtime packaging

The observed `libgcc_s_seh-1.dll` popup was a launch-environment issue rather
than a solver memory error.  MinGW builds now stage `libgcc_s_seh-1.dll`,
`libstdc++-6.dll`, `libwinpthread-1.dll`, and `libspdlog-1.17.dll` into the
build directory at CMake configure time.  Both `test_mobility.exe` and
`sentaurus_import.exe --help` run successfully from a PowerShell process that
does not prepend the UCRT64 directories to PATH.

## Verification

- density-gradient unit tests: 21 cases, 86 assertions passed;
- density-gradient Newton configuration tests: 1 case, 37 assertions passed;
- material/mesh tests: 21 cases, 70 assertions passed;
- mobility executable: 25 cases, 124 assertions passed without an injected
  MSYS2 PATH;
- SingleDevice Python fixture: 6 tests passed;
- `git diff --check`: passed.

