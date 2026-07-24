# PN2D Minimal6 Phase F self-consistent comparison

Date: 2026-07-24

Status: complete with typed outcome `model_difference`.

## Decision

The repaired production operator converges at every exact
`mirror/sketch x -1..-20 V` checkpoint. All 40 Vela states are classified as
`multiplication_like`; no transition fails and terminal-current magnitude
retains the expected reverse-bias order.

The provisional parity targets are not met. The first failed dependency metric
is the electron quasi-Fermi potential on internal nodes 1 and 5. Electrostatic
potential passes before that metric, so Phase F does not identify a Poisson
defect. The remaining density, mobility, directed-current, terminal-current,
and impact-source failures are downstream of the self-consistent QFP state.

No production mobility, Scharfetter-Gummel, impact-ionization, or Poisson
formula change is justified by Phase F.

## Exact sweep contract

Two independent candidate roots were generated from the repaired Release
binary:

- `build-release/pn2d-minimal6-phase-f-sweep-20260724-r1-a`
- `build-release/pn2d-minimal6-phase-f-sweep-20260724-r1-b`

Each root contains:

| Evidence | Count |
|---|---:|
| accepted Vela checkpoints | 40 |
| accepted Sentaurus checkpoints, including 0 V starts | 42 |
| rejected Vela transitions | 0 |
| accepted Vela segments | 40 |
| `multiplication_like` Vela branch classifications | 40 |

The comparison lattice excludes the two 0 V initial states and contains the
exact 40 common physical states. Interpolation remains forbidden.

An earlier background launch root,
`pn2d-minimal6-phase-f-sweep-20260724-a`, is retained but excluded. Both first
segments exited with Windows `0xC0000135` before solver execution because the
background process lacked the UCRT64 DLL path. The valid roots explicitly
inherit `D:/msys64/ucrt64/bin` and do not contain this infrastructure failure.

## Dependency-chain metrics

Thresholds are the pre-run provisional targets in the revised Phase B-G plan.

| Quantity | Count | Median | P95 | Maximum | Gate |
|---|---:|---:|---:|---:|---|
| electrostatic potential, V | 80 | `5.36593e-12` | `1.59494e-11` | `2.36565e-11` | pass |
| electron QFP, V | 80 | `0.0548737` | `0.0737589` | `0.0817169` | fail |
| hole QFP, V | 80 | `0.0620671` | `0.0821749` | `0.0891277` | fail |
| electron density, dex | 80 | `0.921833` | `1.23909` | `1.37278` | fail |
| hole density, dex | 80 | `1.04268` | `1.38047` | `1.49727` | fail |
| electron native-element mobility, dex | 160 | `0.465174` | `0.743976` | `0.775572` | fail |
| hole native-element mobility, dex | 160 | `0.292739` | `0.463503` | `0.488321` | fail |
| electron directed box-edge current, dex | 200 | `0.823474` | `2.22312` | `2.92481` | fail |
| hole directed box-edge current, dex | 200 | `0.939922` | `2.38723` | `2.88210` | fail |
| total terminal current, dex | 40 | `0.647139` | `0.821782` | `0.914590` | fail |
| integrated impact source, dex | 40 | `12.8959` | `13.0327` | `13.1197` | fail |

All 400 active directed-current samples retain typed signs and zero/reference
handling. The Vela carrier-edge convention is replayed from the shared raw
signed particle flux:

- electron box-current orientation: `+q * flux`;
- hole box-current orientation: `-q * flux`; and
- global-edge integration: dual-face length times `1e-6` unit depth.

An independent contact sum of the 400 recorded carrier-edge values reproduces
the Vela Anode current with maximum relative error `3.66527e-16`. Thus the
directed-current error is not an edge-to-terminal sign or integration artifact.

Sentaurus directed current remains explicitly labeled a
`box_operator_reconstruction`, not a native directed-edge observation.

## Scientific interpretation

The source-unit repair materially moved the self-consistent branch toward
Sentaurus: the earlier internal QFP discrepancy was approximately `0.33 V`,
whereas the Phase F medians are `0.0549 V` and `0.0621 V`. It did not reach the
frozen `0.01 V` median and `0.025 V` P95 targets.

The QFP errors predict order-one density errors through the exponential carrier
statistics. The observed approximately `0.92-1.04 dex` density medians are
therefore consistent with the self-consistent QFP displacement. This does not
contradict the fixed-imported-state BGN result: when Sentaurus electrostatic and
QFP potentials are imported, Vela recomputes all node densities within
`4.426181e-6 dex`.

The Phase F mobility rows compare different self-consistent QFP-gradient
states. They must not be interpreted as a new fixed-state mobility-formula
error. Phase D remains authoritative for same-state model/support separation.

The integrated impact source is approximately 13 dex below Sentaurus even
though total terminal current is lower by a median `0.647 dex`. This is
consistent with the exponential sensitivity of the avalanche coefficient to
the carrier driving field. The carrier-edge-to-terminal closure excludes a
missing current integration factor. Phase G must distinguish the remaining QFP
continuity/model difference before considering any impact-model change.

## Regression and determinism

The two Phase F analysis roots are:

- `build-release/pn2d-minimal6-phase-f-self-consistent-20260724-r3-a`
- `build-release/pn2d-minimal6-phase-f-self-consistent-20260724-r3-b`

Both contain 80 internal-node rows, 320 carrier-element mobility rows, 400
active directed-current rows, and 40 terminal/source rows. The independent
verifier passes all hashes and row counts, and the two manifests are identical
in content.

Phase B/C replay gates were independently rerun on both sealed Phase C roots:

| Gate | Result |
|---|---:|
| checked exact states | `40` |
| final valid carrier edges | `400` |
| Vela production replay maximum relative difference | `5.295713e-14` |
| imported-state density maximum error | `4.426180e-6 dex` |
| Sentaurus total terminal-current relative error | `1.137692e-7` |
| total-current KCL relative error | `3.163255e-9` |
| verifier failures | `0` |

## Phase F exit

Phase F exits as:

- primary typed result: `model_difference`;
- first failed metric: `electron_qfp`;
- preserved common prefix: all 40 exact states;
- branch result: 40/40 `multiplication_like`, no order failure;
- passed upstream metric: electrostatic potential;
- downstream failed metrics: QFP, density, mobility on differing state
  supports, directed current, terminal current, and impact source; and
- production decision: retain the Phase E source-unit patch, with no additional
  formula change from Phase F evidence.

Phase G may now perform the final production-decision ledger, full Release/CTest
regression, and independent scientific/code review.
