# SingleDevice Eq. 231 insulator-corner closure (2026-08-14)

## Outcome

The SiO2/Nitride shared-corner failure is closed for both imported
SingleDevice endpoint states. The correction is restricted to a pure
nontransport re-entrant vertex with exactly six incident triangles split 2:4
between two distinct materials. Ordinary nodes along the same material
interface are unchanged.

The same work also closes the associated state/output contract. The
`sentaurus_box` equation already used the maximum material-side band drive at
a shared vertex, but initialization, Dirichlet conversion, diagnostics, and
the returned/restart quantum potential used the CSV owner-side drive. At node
1793 those conventions differed by exactly 1 V. All conversions now use the
same maximum material-side drive selected by the equation.

## Direct oracle

Sentaurus Formula-0 finite differences at node 1793 show that all six Vela
off-diagonal flux/Jacobian entries have one common normalization relative to
Sentaurus. The inferred Sentaurus reaction diagonal is smaller than Vela's
generic control-volume reaction by the same factor obtained directly from the
fixed-state balance:

- Jacobian-derived reaction ratio: `0.971387583603058`;
- direct fixed-state closure ratio: `0.9713027870586473`.

The direct balance value is used by the SingleDevice experimental profile as
`sentaurus_insulator_reentrant_corner_reaction_weight`. The parameter defaults
to `1` and therefore does not alter generic simulations.

## Fixed-state cross-state result

| endpoint | node | old residual | closed residual |
|---|---:|---:|---:|
| linear, Vg=2.2 V, Vd=0.1 V | 1793 | `1.5975140900` | `0` |
| linear, mirrored corner | 3552 | `1.5975` | `1.61e-12` |
| saturation, Vg=2.2 V, Vd=1.1 V | 1793 | `1.5975148224` | `-1.91e-11` |
| saturation, mirrored corner | 3552 | `1.5975` | `2.69e-12` |

Control nodes 47, 1788, and 1794 on the ordinary SiO2/Nitride interface retain
their pre-existing near-zero residuals. The largest remaining fixed-state
residual is no longer an insulator/insulator row: it is about `0.447` at the
Si/SiO2 interface (nodes 2630/848 in the linear state).

## Endpoint effect

Both classical solves and both quantum inner solves converge. Correcting the
shared-node state/output convention reduces the one-step outer quantum change
from exactly `1.000000 V` to:

| endpoint | inner iterations | raw quantum change |
|---|---:|---:|
| linear | 19 | `0.000957 V` |
| saturation | 19 | `0.000926 V` |

These changes are still above the existing `0.5 mV` endpoint acceptance
threshold, so the 21-point Id-Vg curves remain gated. The next residual target
is the Si/SiO2 interface control-volume closure, not the SiO2/Nitride corner or
restart conversion.

## Verification

- density-gradient quantum tests: 22 cases, 93 assertions passed;
- Newton/configuration tests: 83 cases, 1161 assertions passed;
- manufactured 2:4 two-insulator corner detects the new weight and verifies a
  restart-consistent root;
- a mismatched output-owner manufactured test verifies maximum material-side
  drive conversion at a shared vertex;
- `git diff --check`: passed.
