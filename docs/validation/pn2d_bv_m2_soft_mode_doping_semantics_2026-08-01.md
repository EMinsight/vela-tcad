# PN2D BV M2 soft-mode projection and doping-semantics audit

Date: 2026-08-01

Status: P0-P2 complete; observation-only; production defaults unchanged.

Typed outcomes:

```text
carrier_block_linear_solve_decomposed
nodal_doping_topology_and_edge_average_exact__sentaurus_control_volume_not_exported
```

## Decision

The M2 knee discrepancy is not explained by a different nodal doping input,
triangle topology, compensated-node ownership, or edge endpoint-average
doping.  The two dominant carrier soft modes are created by modal cancellation
between transport stiffness and avalanche diagonal/cross-carrier stiffness.
The remaining geometry question is narrower: Vela uses barycentric node
volumes, while the local Sentaurus TDR export does not expose its node box
volumes for direct comparison.

This result does not authorize a production-default change, a doping change,
an avalanche-coefficient adjustment, or an acceptance-threshold change.

## P0: semantics correction and evidence freeze

The prior phrase "compensated-junction triangle doping smoothing" was
incorrect for the active BV configuration.  The frozen configuration uses:

- nodal donors and acceptors read from `node_doping_file`;
- nodal net doping `ND-NA` in Poisson fixed charge;
- `doping_concentration_basis=net_doping` for mobility;
- arithmetic mean of the two endpoint nodal net values for SG-edge mobility;
- default barycentric `area/3` node volumes; and
- no `cell_reconstructed_total_impurity` branch.

The pre-experiment inputs, configuration, prior result, and acceptance
thresholds are frozen in:

```text
docs/validation/contracts/pn2d_bv_m2_soft_mode_p0_freeze_v1.json
```

## P1: modal Jacobian and residual projection

The carrier-block probe was extended without changing the assembled equations.
For every singular pair `(u, v)`, it records:

```text
u^T J_transport v
u^T J_recombination v
u^T J_avalanche_diagonal v
u^T J_avalanche_cross v
u^T (-R_transport)
u^T (-R_recombination)
u^T (-R_avalanche)
```

It also projects the `no_cross_carrier`, `no_recombination`, `no_avalanche`,
and `transport_only` linear-solve steps onto the same full-system right
singular vectors.

Two complete runs over `-18`, `-19.5`, `-19.7`, and `-20 V`, for both Vela
baseline and joint Sentaurus-QFP states, are byte-identical.  Contract checks:

| Check | Maximum | Limit | Result |
|---|---:|---:|---|
| modal Jacobian relative closure | `5.96e-15` | `1e-10` | pass |
| modal RHS relative closure | `9.62e-14` | `1e-10` | pass |
| full linear closure | `7.19e-16` | `1e-8` | pass |
| row-scaled/full step difference | `2.10e-13` | `1e-8` | pass |
| SVD energy closure | `4.88e-15` | `1e-10` | pass |

At `-20 V`, the first two joint-QFP modes carry `58.70%` and `36.18%` of
the production step energy.  Their signed Jacobian projections, normalized by
the full singular value, are:

| Mode | transport | recombination | avalanche diagonal | avalanche cross | sum |
|---:|---:|---:|---:|---:|---:|
| 1 | `+9.1821` | `+2.91e-6` | `-3.2702` | `-4.9119` | `1.0000` |
| 2 | `+8.3493` | `+2.65e-6` | `-3.1866` | `-4.1627` | `1.0000` |

Thus transport by itself is not near-null.  Avalanche diagonal and
cross-carrier terms cancel `89.1%` and `88.0%` of the transport modal
stiffness, respectively, leaving the very small full singular values.  The
cross-carrier term is the larger avalanche contribution, but the diagonal and
cross terms jointly form the cancellation.  Recombination is negligible in
the modal Jacobian.

The modal RHS is also a cancellation.  At `-20 V`, transport contributes
`-18.64/-16.49` times the net RHS projection and avalanche contributes
`+18.61/+16.46` times; recombination supplies approximately `1.03` times the
small remainder.  This explains why small state changes can strongly change
the carrier step without a local sign error.

## P2: Sentaurus/Vela input and control-volume semantics

The original M2 Sentaurus mesh TDR was re-imported independently with
`compensated-doping-policy=reported`.  Results:

| Comparison | Result |
|---|---:|
| node IDs | exact, `115/115` |
| coordinates | maximum error `0 um` |
| unordered triangle topology | exact, `191/191` |
| donors | maximum relative error `0` |
| acceptors | maximum relative error `0` |
| derived nodal net doping | maximum relative error `0` |
| Sentaurus reported net versus `ND-NA` | maximum relative error `1.92e-15` |
| junction edge endpoint-average net doping | maximum relative error `0` |
| compensated nodes | `9` in both representations |

The Vela and Sentaurus-exported `doping.csv` files have the same SHA-256:

```text
c4fe644a164b90ea954b6f844326349224ac957f0b06da99d53851799e8b1165
```

Therefore the soft modes cannot be attributed to a different nodal doping
file, a hidden triangle total-impurity reconstruction, or different junction
connectivity.

Vela's active barycentric and diagnostic mixed-Voronoi volumes both conserve
the total `1e-12 m2` device area.  However, the six nodes selected by the two
dominant modes have mixed/barycentric volume ratios from `1.0` to `1.5`.
This includes the strongest p-side shoulder node at `(0.75, 0.5) um` and an
n-side boundary shoulder node at `(1.25, 0) um`, both with ratio `1.5`.

The TDR inventory contains coordinates, elements, donor/acceptor fields, and
reported net doping, but no node control-volume field.  Consequently this
audit proves input and topology parity but cannot claim exact Sentaurus/Vela
control-volume parity.

## Interpretation

The present evidence supports this narrower chain:

```text
joint-QFP state
  -> transport and avalanche modal stiffness nearly cancel
  -> two junction-shoulder near-null carrier modes
  -> large carrier-only Newton step
  -> self-consistent density/current/source deficit
```

Nodal doping determines where the depletion shoulders reside, but the values
and topology are identical between the two simulators.  It is therefore an
indirect spatial organizer, not a mismatched input.  A control-volume policy
difference remains testable because the soft-mode support includes nodes with
up to `1.5x` barycentric/mixed-Voronoi volume difference.

## Next bounded experiment

Keep SG/Laux and nodal doping unchanged.  Obtain the Sentaurus box volumes if
possible.  Otherwise, compare barycentric and mixed-Voronoi node volumes in a
frozen-state opt-in first-step experiment on the same four M2 states.  Measure
the two modal stiffnesses, modal RHS, carrier step, Poisson/carrier residuals,
and frozen integrated source.  A nodal-doping redistribution experiment is not
justified unless this control fails.

## Reproduction

```powershell
python scripts\run_pn2d_bv_m2_carrier_block_decomposition.py `
  --runner build-release\vela_example_runner.exe `
  --base-config build-release\pn2d-bv-m2-qfp-carrier-jacobian-verification-20260731\self_consistent_sg_laux_probe.json `
  --prior-root build-release\pn2d-bv-m2-qfp-carrier-jacobian-verification-20260731 `
  --output-root build-release\pn2d-bv-m2-soft-mode-component-projection-20260801 `
  --repeats 2

build-release\sentaurus_import.exe `
  --tdr build-release\pn2d-task10-balanced-m2-sentaurus-process-v2-20260731\avalanche_on\raw\pn2d_msh.tdr `
  --export-dir build-release\pn2d-bv-m2-doping-control-volume-audit-20260801\sentaurus_export `
  --compensated-doping-policy reported

python scripts\audit_pn2d_bv_m2_doping_control_volume_semantics.py `
  --vela-config build-release\pn2d-bv-m2-qfp-carrier-jacobian-verification-20260731\self_consistent_sg_laux_probe.json `
  --sentaurus-export build-release\pn2d-bv-m2-doping-control-volume-audit-20260801\sentaurus_export `
  --soft-modes build-release\pn2d-bv-m2-soft-mode-component-projection-20260801\dominant_singular_modes.csv `
  --output-root build-release\pn2d-bv-m2-doping-control-volume-audit-20260801
```

Generated simulation and audit outputs remain under `build-release/` and are
not committed.
