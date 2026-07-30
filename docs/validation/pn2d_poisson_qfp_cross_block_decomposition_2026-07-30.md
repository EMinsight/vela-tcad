# PN2D Poisson-QFP cross-block decomposition

Date: 2026-07-30

Outcome: `bidirectional_poisson_qfp_closed_loop_cause`

Production defaults changed: no

Task 8 authorized: no

## Question and controlled contract

The observation determines why a QFP-only carrier-block update points toward
the exact Sentaurus avalanche-on state while the complete Vela Newton update
points away from it near the BV knee.

For one production baseline residual and Jacobian at each bias, the unknowns
are partitioned as `x = [psi; qfp]`, with
`qfp = [phin; phip]`:

```text
J = [ A  B ] = [ J_psi_psi  J_psi_qfp ]
    [ C  D ]   [ J_qfp_psi  J_qfp_qfp ]
```

The diagnostic uses the frozen Sentaurus QFP substitution residual and solves:

1. independent blocks: `A dpsi = -rpsi`, `D dqfp = -rqfp`;
2. remove `B`: solve `A dpsi = -rpsi`, then
   `D dqfp = -rqfp - C dpsi`;
3. remove `C`: solve `D dqfp = -rqfp`, then
   `A dpsi = -rpsi - B dqfp`;
4. full Schur complement:
   `(D - C A^-1 B) dqfp = -rqfp - C A^-1(-rpsi)`;
5. the production full raw Newton solve and its configured capped/Poisson
   recorrected form.

All cases share the same Jacobian, scaling, contact constraints, state,
residual, and physical models. No candidate or production code path consumes
the counterfactual updates.

## Results

| Bias | Independent direction cosine | Full Schur direction cosine | Full capped direction cosine | Independent optimal projection scale | Schur/full max difference | Relative closure |
|---:|---:|---:|---:|---:|---:|---:|
| -19.7 V | 0.75664 | -0.44177 | -0.19941 | 2.7124e-5 | 3.9562e-13 V | 2.1966e-14 |
| -19.8 V | 0.75108 | -0.45824 | -0.21322 | 2.8129e-5 | 7.3193e-13 V | 4.5993e-14 |

Removing either cross block preserves the positive QFP direction of the
independent carrier block at both biases. Closing both directions through
`C A^-1 B` reverses it. The exact Schur reconstruction matches the production
raw full solve within `7.32e-13 V`, proving that the reversal is present before
configured step caps or Poisson recorrection. The capped step remains adverse,
so those post-solve operations attenuate but do not create the reversal.

The independent carrier block is not a usable correction: its unscaled unit
step has QFP RMSE near `1.19e4/1.15e4 V`, and its target-projection optimal
scale is only about `2.7e-5/2.8e-5`. This is direction evidence, not evidence
that decoupling or deleting a Jacobian block is numerically or physically
valid.

The strongest adverse cross-block projection is reproducibly located at node
12, `(x, y) = (1.0, 0.5)`, at both biases. Six contact nodes per bias preserve
zero target QFP change.

## Determinism and provenance

Two independent observation roots produced identical node and sparse-block
CSV hashes:

| Bias | Node CSV SHA-256 | Jacobian-block CSV SHA-256 |
|---:|---|---|
| -19.7 V | `2d4dd6d624685b5e2c4c1075d6c150502f6cba09ea6efcdc5445f267d980cccb` | `62c37c48ef67b44752a129533f6215664853064fe7d4bf2573fcf6cee53afb77` |
| -19.8 V | `0a596ad6bda801ba716fd65c3b34c1417553ded19f1ce652be28d099cca8f4fc` | `efea45c6776e9b71d00efe4c6f39ad7087edc12d4d506a7ca8b33f5eeb620836` |

Sealed inputs:

- base config SHA-256:
  `18df2306e6397c27f3b1ed14ddda4781b261c64aaf8a318529c24a2b7d5e2975`;
- Sentaurus process-chain input SHA-256:
  `d83f377d6377ac7fb881df48ba64c0a997ad1fcd3fb2665ce8ead108124b8c28`;
- Vela exact-state manifest SHA-256:
  `1c020a74510aa82566bdcbf970fa12cae47ffda826c93f0bf66b0ea508879f20`;
- diagnostic runner SHA-256:
  `e9d73111329ec79c74ea601cf5204fa231fa2a4e1d15ab4c3567f01d28f4c3be`.

Generated scorecard:
`build-release/pn2d-task7-cross-block-scorecard-20260730/acceptance.json`.

## Decision

The cross-bias causal classification is
`bidirectional_poisson_qfp_closed_loop_cause`. It localizes the first-update
direction reversal to the closed Poisson-QFP Jacobian loop; it does not yet
identify whether the responsible issue is model ownership, sign convention,
state scaling, conditioning, or a legitimate local linearization of different
physics.

The requested next observation has been completed in
`pn2d_schur_loop_source_decomposition_2026-07-30.md`. Transport and avalanche
loops independently sustain the adverse direction, while SRH/Auger does not;
directional finite differences pass. Task 8 remains prohibited pending
analytical model-ownership and independent-reference review.
