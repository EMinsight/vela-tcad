# PN2D Task 10 M0 stall root-cause audit

Date: 2026-07-31

## Outcome

The Vela M0 stall near `-17.1054 V` is resolved without changing the solver,
avalanche coefficients, SG/Laux operator, Newton tolerances, or production
defaults.

The failed mesh input preserved integrated impurity dose by assigning the
geometric junction nodes to the N region only. That creates a nonzero nodal
net doping at `x=1 um`. Because the nodal values are linearly interpolated
inside each triangle, the zero-net-doping contour no longer lies on the
geometric junction. This changes the high-field branch before avalanche
feedback is enabled.

The minimal corrected construction is:

```text
x < 1 um:       NA=1e17, ND=0
x = 1 um:       NA=5e16, ND=5e16
x > 1 um:       NA=0, ND=1e17       [cm^-3]
```

It preserves:

- zero net doping at the intended metallurgical junction;
- total impurity `ND+NA=1e17 cm^-3` at the junction;
- integrated total-impurity dose `1e17 cm^-3 um^2`;
- the existing M0 mesh (`27` nodes, `32` triangles).

## Reproduction and rejected-state capture

The original unconstrained continuation retained its last accepted parent at
`-17.105395800637609 V` and kept shrinking the attempted step to
`3.1616419704505461e-8 V`. It was terminated while continuing to retry, rather
than returning a natural terminal failure.

For a bounded reproduction, the same N-owned input was run with
`min_step=1e-5 V`. It stopped at `-17.10499053824512 V` after 16 rejected
attempts. The first rejection was captured as:

| Quantity | Value |
| --- | ---: |
| Parent bias | `-17.100000000000001 V` |
| Requested segment target | `-18 V` |
| Attempt target | `-17.150000000000002 V` |
| Attempted step | `-0.050000000000000711 V` |
| Initial residual norm | `1.4142135623730951` |
| Final residual norm | `1.7191292957326368e-6` |
| Newton iterations | `10` |
| Rejection | `line_search_non_decrease` |

The parent and initial-state hashes are both `a07a214d880ef9cf`. The optional
diagnostic therefore captures the true last accepted state, the exact initial
state presented to Newton, and the final rejected trial state.

## Controlled junction-node matrix

All five cases below use the same mesh, physics, solver, voltage lattice, and
SG/Laux self-consistent avalanche operator. Only the three nodes at `x=1 um`
change.

| Junction-node rule | Junction net (`cm^-3`) | Total dose (`cm^-3 um^2`) | `I(-17 V)` (`A/um`) | Result |
| --- | ---: | ---: | ---: | --- |
| old double species | `0` | `1.125e17` | `-5.87710e-17` | reaches `-20 V`, 0 rejects |
| N-owned | `+1e17` | `1.000e17` | `-7.12497e-15` | stalls near `-17.105 V` |
| balanced half | `0` | `1.000e17` | `-4.78734e-17` | reaches `-20 V`, 0 rejects |
| neutral zero | `0` | `0.875e17` | `-4.42410e-17` | reaches `-20 V`, 0 rejects |
| P-owned | `-1e17` | `1.000e17` | `-7.08058e-15` | stalls near `-17.105 V` |

The P-owned and N-owned cases fail at the same magnitude with opposite doping
sign. Both zero-net controls reach `-20 V`. This rules out a donor-specific
implementation defect and identifies the nonzero junction-node net doping as
the causal variable. Of the passing controls, only balanced-half also preserves
the intended total impurity.

## Physical-state comparison at -17 V

Representative values show that the N-owned input is already on a different
high-field/current branch:

| Quantity | old double species | N-owned | balanced half |
| --- | ---: | ---: | ---: |
| Junction potential (`V`) | `-8.5000` | `-2.31287` | `-8.5000` |
| Junction electron QFP (`V`) | `-8.10804` | `-2.05303` | `-8.10818` |
| Junction hole QFP (`V`) | `-8.88193` | `-2.58001` | `-8.88183` |
| Junction electron density (`m^-3`) | `5.106e9` | `7.143e11` | `4.333e9` |
| Junction hole density (`m^-3`) | `7.526e9` | `5.385e11` | `6.375e9` |
| Top impact field (`V/m`) | `3.24321e5` | `4.94418e5` | `3.24327e5` |
| Top alpha (`m^-1`) | `1.57957e4` | `5.82977e4` | `1.57968e4` |
| Integrated source | `1.59578e14` | `4.41558e16` | `1.35886e14` |

At `-17 V`, the N-owned current is about `121x` the old input and its integrated
source is about `277x`. Balanced-half restores the old physical branch while
removing the old input's 12.5% total-dose excess.

## First-failure branch decomposition

At the first rejected N-owned attempt (`-17.10 -> -17.15 V`):

- avalanche-off and IIC have identical residuals and Newton steps because the
  IIC source is observed but not fed back;
- their full-step trial combined residual is `1.32467e-3`;
- self-consistent avalanche-on changes the full coupled step and leaves a
  `5.52903e-3` trial combined residual;
- IIC and on both observe carrier impact-term L2 norm
  `2.03599e-5`, almost exactly balancing the SG flux L2 norm
  `2.03635e-5`;
- the largest impact term is at junction node `10`;
- off has a zero impact term.

The stall is therefore source-feedback-enabled, but the avalanche model is not
the upstream root cause: the invalid junction-node state first inflates the
field, carrier densities, SG current, and ionization source. Feeding that source
back then exposes the continuation frontier.

## Corrected Sentaurus M0 input

Sentaurus SDE O-2018.06-SP2 regenerated the balanced-half M0 TDR. Importing the
TDR confirms all three junction nodes have:

```text
ND = NA = 4.999999999999996e16 cm^-3
```

The imported mesh is byte-identical at the JSON level to the failed M0 mesh:

```text
mesh SHA-256:
c9aaf5f3130f2e1e78e399d155390ed8f19a306ff9ab5af4904230b5e328bc7e

Sentaurus TDR SHA-256:
999057b81108361fda6ffdf6c0d8cc40cdfc44e00b7f7b7f2b5b6d6d0a32bac0

imported doping SHA-256:
926dfc4de1d9322b43723ddcf087dd94d6af6c1b61e451e0434103f166ad1b11
```

## Exact-lattice confirmation

Two fresh independent Vela executions used the actual Sentaurus-imported
balanced M0 mesh and doping. In each execution:

- avalanche-off: `29/29` exact/converged points, 0 rejected attempts;
- IIC postprocess: `29/29` exact/converged points, 0 rejected attempts;
- SG/Laux avalanche-on: `29/29` exact/converged points, 0 rejected attempts;
- all branches reach `-20 V`;
- all 87 state hashes and all three IV hashes match across the two executions.

At `-20 V`, the currents are:

| Branch | Current (`A/um`) |
| --- | ---: |
| avalanche-off | `-6.24657e-17` |
| IIC postprocess | `-6.24657e-17` |
| SG/Laux avalanche-on | `-4.11185e-16` |

## Paired Sentaurus curve

The same balanced M0 TDR was then run in Sentaurus for off, IIC, on, and
on-with-avalanche-derivative observation. The normalized manifest contains four
complete 29-point branches, 176,668 field records, and 1,044 aggregate records:

```text
manifest SHA-256:
ed99ec9474c60931675947d8e4cd6ca8f1a04531e911b6732c0322bb439f2e5e
```

Sentaurus avalanche-on reaches `-20 V` with
`I=-3.94901e-16 A/um`, compared with Vela
`I=-4.11185e-16 A/um`. The paired curve metrics are:

| Metric | Observed |
| --- | ---: |
| Global median absolute log-current error | `0.00579 dex` |
| Global P95 absolute log-current error | `0.01755 dex` |
| Knee median absolute log-current error | `0.01755 dex` |
| Knee maximum absolute log-current error | `0.03010 dex` |
| Adjacent-slope RMSE | `0.08924 dex/V` |
| Sentaurus `V_break` | `-19.654 V` |
| Vela `V_break` | `-19.622 V` |
| Absolute `V_break` difference | `0.032 V` |
| Sentaurus/Vela `V_slope` | `-19.793/-19.828 V` |

The fail-closed parity wrapper reports `ill_conditioned_knee_metric`, not a
solver failure. Its strict composite result is held open by the already known
near-floor low-current nonmonotonicity and by closure columns that are not
present in these raw curve CSVs. The BV-effective knee-window current and
voltage agreement is nevertheless preserved by the corrected input. This
paired result closes the specific “Vela cannot reach -20 V on corrected M0”
blocker; it does not waive the separate closure or low-current contracts.

## Implementation

The change is deliberately split into input-contract and observation-only
parts:

1. `prepare_pn2d_bv_off_srh_mesh_matrix.py` now emits a narrow junction window
   containing equal half-concentration donor and acceptor species.
2. `DCSweep` accepts an optional
   `sweep.diagnostics.newton_history.rejected_state_directory`. When set, only
   rejected attempts receive parent/initial/final state snapshots and their
   paths are appended to `newton_attempts.csv`.
3. Controlled junction-ownership and first-failure audit scripts preserve all
   other inputs and generate machine-readable summaries.

The new diagnostic option defaults to empty. It does not affect residuals,
Jacobians, continuation, or accepted states when disabled.

## Decision

The M0 Vela blocker is closed. The prior M1-to-M2 mesh-independence result is a
separate, secondary gate and remains open because the earlier M1/M2 runs used
the superseded single-owner junction contract. No production default changes
are authorized by this result.
