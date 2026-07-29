# PN2D avalanche-off same-state SRH decomposition

Date: 2026-07-28

## Method

At `-1, -5, -10, -15, -20 V`, native Sentaurus `psi`, `n`, `p`,
electron/hole quasi-Fermi potentials, and SRH rate were read from the
27-node Sentaurus coarse7x3 mesh. The converted Vela mesh has the same
`mesh.json` SHA-256, so the recorded barycentric map is an identity-geometry
map rather than a fine-to-coarse projection.

Every target node records its source triangle and weights. There are no
duplicate target IDs, uncovered nodes, or implicit zero fills. Because the
TDR does not expose native effective intrinsic density, it is explicitly
inferred independently from the electron and hole density/quasi-Fermi
relations; both inferred values and their disagreement are retained in the
local table.

## Acceptance evidence

Artifacts:

`build-release/pn2d-bv-off-srh-same-state-coarse-baseline-20260728/`

| Gate | Result | Limit | Status |
| --- | ---: | ---: | --- |
| required anchors | 5/5 | 5/5 | pass |
| target-node coverage | 100% | 100% | pass |
| duplicate IDs | 0 | 0 | pass |
| implicit zero fills | 0 | 0 | pass |
| repeated same-state error | 0 | `1e-10` | pass |
| max nodal/element integration difference | `4.077e-16` | 1% | pass |
| minimum named-term assignment | `0.9999999999990323` | 90% | pass |

The element-integration conservation and node-permutation tests pass.

## Log-gap decomposition

All values below are signed dex contributions in the
Sentaurus-over-Vela direction.

| Bias (V) | state `n,p` | inferred `ni_eff`/BGN | local formula | support | same-state total mismatch |
| ---: | ---: | ---: | ---: | ---: | ---: |
| -1 | `-1.48e-8` | `-4.42e-6` | `+1.20e-6` | `0` | `1.20e-6` |
| -5 | `-4.31e-9` | `-4.41e-6` | `+1.18e-6` | `0` | `1.18e-6` |
| -10 | `-2.62e-9` | `-4.40e-6` | `+1.18e-6` | `0` | `1.18e-6` |
| -15 | `-2.20e-9` | `-4.40e-6` | `+1.17e-6` | `0` | `1.17e-6` |
| -20 | `-1.48e-4` | `-4.99e-6` | `+1.04e-6` | `0` | `1.04e-6` |

The source-dominant local Vela-operator versus native-Sentaurus rate mismatch
is only `1.20e-6 dex` or less. The same-state integrated mismatch is also
only `1.20e-6 dex` or less because coarse and native support are now the same.

## Decision gate

Classification at the Task 2 gate:

`state_difference`

Imported Sentaurus state plus the existing Vela SRH operator matches the
native Sentaurus integral within `0.05 dex` at every anchor. Task 3 is not
entered: there is no local-rate evidence for an SRH formula or parameter
mismatch. The superseded non-paired mesh result must not be carried forward.
