# PN2D WP7 paired process-chain analyzer

Date: 2026-07-30

Outcome: `density_qfp_feedback_cause`.

## Implemented contract

`scripts/analyze_pn2d_bv_process_chain.py` compares:

1. Sentaurus versus Vela avalanche-off;
2. Sentaurus versus Vela IIC/postprocess-only;
3. Sentaurus versus Vela avalanche-on;
4. Sentaurus IIC versus avalanche-on; and
5. Vela postprocess-only versus avalanche-on.

The fixed dependency order is state, density, drive, mobility, current,
alpha, generation, geometric source, residual/Jacobian, first Newton update,
and terminal current. The analyzer emits:

- `stage_summary.csv`;
- `support_summary.csv`;
- `hotspot_chain.csv`;
- `first_departure.json`;
- `source_terminal_closure.csv`;
- `newton_first_update.csv`;
- deterministic process-chain and hotspot SVGs; and
- `acceptance.json`.

It reports scalar signed/relative/log differences, vector component/magnitude
and angle differences, active-support overlap, hotspot, source-weighted
centroid, and cumulative 10/50/90% support.

## Focused acceptance

Synthetic fixtures inject an isolated error at each of all eleven stages.
Every injection is recovered as the first stage at both adjacent knee biases.
Row reordering and unrelated sub-floor tail values do not change the result.
A one-bias injection and a missing simulator both return
`insufficient_observation`.

A causal result also requires:

- every stage in all paired comparisons;
- matching non-implicit supports;
- source and terminal closure at `1e-12`; and
- the same earliest stage at two adjacent biases.

## Completed real-data gate

The current run binds the exact 29-point Sentaurus process manifest:

`build-release/pn2d-wp21-full-lattice-a-20260729/manifest.json`

with SHA-256
`190cb08f6c128ce64bdfd9bb8dfc6242bde95238234b7492ef5740b4fa2d3d15`.

The run-local Vela `max_iter=80` qualification supplies all 29 exact points
for avalanche-off, IIC/postprocess-only, and avalanche-on. Its process
manifest SHA-256 is
`b882ece81a9cd1e7633e5685adbdd1a9ffde8b4adf2d14dea4fbc2286d6ddf6d`.

Sentaurus `NewtonPlot(Error Residual Update)` and Vela
`newton_step_probe` were then evaluated on the same fixed-transition
contract: each knee target starts from the previous exact accepted state and
only the first Newton correction is observed. The six targets are
`-19.7, -19.8, -19.85, -19.9, -19.95, -20 V`.

Each simulator contributes 2,916 node records:

- per-equation L2-normalized residual signatures;
- potential and carrier-density first-update signatures; and
- three branches times six biases times 27 nodes times six quantities.

All 2,916 cross-simulator node coordinates agree exactly. Sentaurus and Vela
each have zero avalanche-off/IIC difference over their 972 corresponding
records. The full Sentaurus Jacobian matrix is not exported by the 2018.06
NewtonPlot interface; the contract records the native RHS and the observable
inverse action `delta_x = -J^-1 R` and explicitly sets matrix availability to
false.

The final analyzer output is under
`build-release/pn2d-wp7-process-chain-newton-complete-20260730`. All five
comparisons contain all eleven stages, all 203 source/terminal closure rows
pass, and `missing_stage_observations` is empty.

The accepted earliest departure is `state` for Sentaurus IIC versus
avalanche-on at adjacent biases `-19.7/-19.8 V`. The affected evidence is
minority-carrier quasi-Fermi/density feedback; for example, at physical node
14 and `-19.7 V`, the hole quasi-Fermi potential changes from
`-0.3529632468 V` to `-0.2646053763 V`, while the corresponding low-density
tail is deliberately excluded from density-stage causality by the frozen
active-support floor.

Typed outcome: `density_qfp_feedback_cause`. This satisfies the WP8/Task 7
entry condition, but it does not itself authorize a production-default
change.
